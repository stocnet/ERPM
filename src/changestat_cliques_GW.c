/**
 * @file
 * @brief  Change statistic for the ERPM term `cliques_GW` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `cliques_GW(lambda)`, which aggregates group-level k-cliques using a
 *  geometrically weighted series.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (groups = group mode, actors = actor mode)
 *  ------------------------------------------------------------
 *
 *  A bipartite network is assumed, with:
 *    - group mode  = structural groups (one membership per actor),
 *    - actor mode  = individuals whose membership edges are toggled.
 *
 *  If a group has size d (degree in the bipartite projection), its closed-form
 *  GW-clique contribution is:
 *
 *        S(d, λ) = λ * [1 - ((λ - 1) / λ)^d ]
 *
 *  Let r = (λ - 1) / λ.  
 *  A single membership toggle changes the group size from d_old to d_new, giving:
 *
 *        Δ = S(d_new, λ) − S(d_old, λ)
 *          = λ * ( r^{d_old} − r^{d_new} )
 *
 *  Only the group affected by the toggle contributes to the change.
 *
 *  ------------------------------------------------------------
 *  Implementation in {ergm} (one-toggle)
 *  ------------------------------------------------------------
 *
 *  - This function is called once per toggle.
 *  - `edgestate == 1` → the edge exists → toggle = deletion → d_new = d_old − 1  
 *    `edgestate == 0` → the edge does not exist → toggle = addition → d_new = d_old + 1
 *  - `ZERO_ALL_CHANGESTATS(0)` clears the output buffer for this toggle.
 *  - {ergm} accumulates the values from all toggles afterward.
 *
 *  Group vertex detection:
 *    The group node is always the endpoint located in the “group mode”
 *    (index > BIPARTITE). The other endpoint is always an actor.
 *
 *  Vectorisation over λ:
 *    INPUT_PARAM = [λ_0, r_0, λ_1, r_1, ..., λ_{J-1}, r_{J-1}]
 *    where r_j = (λ_j − 1) / λ_j.
 *
 *  Numerical notes:
 *    - λ ≥ 1 ⇒ r ∈ [0,1)
 *    - λ = 1 ⇒ r = 0, r^0 = 1, r^d = 0 for d > 0 (formula remains valid)
 *    - `dpow_double()` handles exponentiation via exponentiation-by-squaring (O(log d))
 *    - We use convention 0^0 = 1 for combinatorial consistency
 *
 *  Complexity:
 *    O(#lambdas) per toggle.  
 *    No adjacency traversal and no clique enumeration.
 *
 *  R interface:
 *    Provided by InitErgmTerm.cliques_GW:
 *      - validates λ and builds INPUT_PARAM = c(rbind(lambda, r))
 *      - emptynwstats = 0, coef.names vectorized
 *
 *  Output:
 *    For each λ_j:
 *
 *        CHANGE_STAT[j] += λ_j * ( r_j^{d_old} − r_j^{d_new} )
 *
 *  ------------------------------------------------------------
 *  @example Usage (R)
 *  ------------------------------------------------------------
 *  ```r
 *  library(ERPM)
 *
 *  # Example partition of actors into groups:
 *  part <- c(1,1,2,2,3,3)   # 6 actors, 3 groups
 *
 *  # Fit a model with one GW-clique term
 *  fit <- erpm(partition ~ cliques_GW(lambda = 2))
 *  summary(fit)
 *
 *  # Interpretation of one toggle:
 *  # Suppose the sampler proposes to add an actor to a group:
 *  #   d_old = 2
 *  #   d_new = 3
 *  #   λ = 2  →  r = (2−1)/2 = 1/2
 *  #
 *  #   Δ = 2 * ((1/2)^2 − (1/2)^3)
 *  #     = 2 * (1/4 − 1/8)
 *  #     = 1/4
 *  #
 *  # This is exactly what c_cliques_GW adds to CHANGE_STAT[0].
 *  ```
 *
 */



#include <math.h>
#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"

/**
 * @def DEBUG_CLIQUES_GW
 * @brief Enable verbose debugging output for ::c_cliques_GW.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 * - group vertex index and degrees before/after the toggle,
 * - λ and r values for each sub-term,
 * - intermediate powers r^d and the resulting Δ.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_CLIQUES_GW 0

/* -------------------------------------------------------------------------- */
/* Utility: fast exponentiation                                               */
/* -------------------------------------------------------------------------- */

/**
 * @brief Fast exponentiation for non-negative integer exponents.
 *
 * @details
 *  Computes base^exp using exponentiation-by-squaring (O(log exp)).
 *  Returns 1.0 when exp ≤ 0, implementing 0^0 = 1 for combinatorial
 *  consistency.
 *
 * @param base  The base in double precision.
 * @param exp   Non-negative integer exponent.
 * @return base^exp if exp > 0, else 1.0.
 */
static inline double dpow_double(double base, int exp){
  if(exp <= 0) return 1.0;  // convention base^0 = 1.0 for all base (including 0^0)

  double r = 1.0;
  double b = base;
  int e = exp;

  /* Exponentiation by squaring:
   * - invariant: r * b^e = base^exp at each step.
   * - when the least significant bit of e is 1, multiply r by current b.
   * - square b and shift e until e becomes 0.
   */
  while(e){
    if(e & 1) r *= b;
    b *= b;
    e >>= 1;
  }
  return r;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cliques_GW                                               */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cliques_GW(lambda)`.
 *
 * This is the {ergm} change-statistic function registered as
 * ::c_cliques_GW via ::C_CHANGESTAT_FN. It processes a single toggle
 * (tail, head) on the bipartite network:
 *
 * - It identifies the group-side vertex @c v2 (mode 2).
 * - It reads its current degree @c deg_old from the internal degree arrays.
 * - It infers @c deg_new depending on whether the edge is being added or removed.
 * - For each λ_j packed into INPUT_PARAM, it computes
 *   \f$\Delta_j = \lambda_j (r_j^{deg\_old} - r_j^{deg\_new})\f$
 *   and adds it to @c CHANGE_STAT[j].
 *
 * @param tail       Pointer to the tail vertex of the toggled edge.
 * @param head       Pointer to the head vertex of the toggled edge.
 * @param mtp        Pointer to the model term parameters (unused directly here,
 *                   but required by the macro signature).
 * @param nwp        Pointer to the current network-plus workspace (provides
 *                   OUT_DEG, IN_DEG, etc.).
 * @param edgestate  Indicator of the current edge state:
 *                   - 0 if the edge is absent (toggle = addition),
 *                   - 1 if the edge is present (toggle = deletion).
 *
 * @note
 * - The macro ::C_CHANGESTAT_FN provides the full prototype including
 *   the number of statistics (N_CHANGE_STATS) and the arrays CHANGE_STAT
 *   and INPUT_PARAM.
 * - The function assumes that the network is bipartite and that one endpoint
 *   of the dyad is always in mode 1 and the other in mode 2.
 */
C_CHANGESTAT_FN(c_cliques_GW){
  /* 1) Reset the output buffer for THIS toggle.
   *
   * {ergm} calls this function once per toggle and accumulates the
   * resulting CHANGE_STAT values outside this function. We therefore
   * explicitly zero the array at the beginning of each call.
   */
  ZERO_ALL_CHANGESTATS(0);  // equivalent to: memset(CHANGE_STAT, 0, N_CHANGE_STATS * sizeof(double))

  /* 2) Size of mode (actors).
   *
   * For bipartite networks, BIPARTITE > 0 equals the number of actors.
   * Mode 2 (groups) then occupies vertices in (BIPARTITE+1..N_NODES).
   * If BIPARTITE == 0 the network is not bipartite, but this term is
   * designed for the bipartite case only.
   */
  const int n1 = BIPARTITE;

  /* 3) Identify the group-side vertex v2 (mode 2).
   *
   * The toggle may be oriented as (actor -> group) or (group -> actor).
   * We select the endpoint with index strictly greater than n1 as the
   * group vertex.
   */
  Vertex v2 = (tail > (Vertex)n1) ? tail : head;

  /* 4) Degrees before and after the toggle.
   *
   * OUT_DEG and IN_DEG are supplied by the {ergm} engine. For bipartite
   * ERPM terms we sum them to get the actual size of the group in mode 2.
   * - deg_old: degree of the group before the toggle.
   * - delta:   +1 if we add an edge, -1 if we remove an edge.
   * - deg_new: degree after the toggle (must remain non-negative).
   */
  int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);
  int delta   = edgestate ? -1 : +1;   // edgestate=1 -> deletion, edgestate=0 -> addition
  int deg_new = deg_old + delta;

  /* Defensive note: in a well-formed chain, deg_new should never be negative.
   * If that ever happens, it indicates an inconsistency between the stored
   * edgestate and the actual adjacency structure.
   */

  /* 5) Optional debug traces. */
  #if DEBUG_CLIQUES_GW
    Rprintf("[c_cliques_GW] v2=%d | edgestate=%d | deg_old=%d -> deg_new=%d | nstats=%d\n",
            (int)v2, (int)edgestate, deg_old, deg_new, N_CHANGE_STATS);
  #endif

  /* 6) Loop over vectorised sub-terms (different λ_j values).
   *
   * The INPUT_PARAM array is packed as:
   *   [lambda_0, r_0, lambda_1, r_1, ..., lambda_{J-1}, r_{J-1}],
   * where:
   *   - lambda_j is the decay parameter λ_j,
   *   - r_j      = (lambda_j - 1) / lambda_j.
   *
   * For each j, we compute:
   *   term_old = r_j ^ deg_old,
   *   term_new = r_j ^ deg_new,
   *   Δ_j      = lambda_j * (term_old - term_new).
   */
  for(int j = 0; j < N_CHANGE_STATS; ++j){
    double lambda = INPUT_PARAM[2*j + 0];
    double r      = INPUT_PARAM[2*j + 1];  // r in [0,1) when lambda >= 1

    // Evaluate r^{deg_old} and r^{deg_new} using dpow_double().
    double term_old = dpow_double(r, deg_old);
    double term_new = dpow_double(r, deg_new);

    // Local change for this λ_j.
    double d = lambda * (term_old - term_new);

    // Accumulate in the j-th component of the change statistic.
    CHANGE_STAT[j] += d;

    #if DEBUG_CLIQUES_GW
      Rprintf("  j=%d | lambda=%.9g r=%.9g | r^deg_old=%.9g r^deg_new=%.9g | Δ=%.9g\n",
              j, lambda, r, term_old, term_new, d);
    #endif
  }
}



// /**
//  * @file    changestat_cliques_GW.c
//  * @brief   Change statistic {ergm} pour le terme ERPM `cliques_GW` (forme un-toggle)
//  *
//  * @details
//  *  Principe statistique :
//  *  - Pour un groupe de taille d, la contribution agrégée des k-cliques pondérées
//  *    géométriquement vaut S(d,λ) = λ[1 - ((λ-1)/λ)^d].
//  *  - Lors d’un toggle affectant un unique groupe v2, la variation locale est :
//  *        Δ = S(d_new,λ) - S(d_old,λ) = λ( r^{d_old} - r^{d_new} ),  r=(λ-1)/λ.
//  *
//  *  Implémentation {ergm} — forme "un-toggle" :
//  *  - Signature imposée : C_CHANGESTAT_FN(c_cliques_GW)(tail, head, mtp, nwp, edgestate)
//  *  - Chaque appel traite UN toggle (tail, head) dans l’état courant du réseau.
//  *  - edgestate = 1  → l’arête existe (toggle = RETRAIT) → d_new = d_old - 1.
//  *    edgestate = 0  → l’arête n’existe pas (toggle = AJOUT) → d_new = d_old + 1.
//  *  - ZERO_ALL_CHANGESTATS(0) : remet le tampon de sortie à zéro.
//  *    {ergm} additionne ensuite les contributions de tous les toggles.
//  *
//  *  Hypothèses et conventions internes :
//  *  - Réseau biparti : mode 1 = acteurs (1..BIPARTITE), mode 2 = groupes (BIPARTITE+1..N).
//  *  - Le toggle connecte toujours un acteur et un groupe → un seul groupe impacté par appel.
//  *  - Le "degré" d’un groupe est le degré total OUT_DEG[v2] + IN_DEG[v2] (représentation interne).
//  *  - Vectorisation sur j : le terme peut exposer plusieurs λ_j. Les entrées INPUT_PARAM sont
//  *    empaquetées par colonnes : [λ_0, r_0, λ_1, r_1, ..., λ_{J-1}, r_{J-1}], avec r_j = (λ_j-1)/λ_j.
//  *
//  *  Précisions numériques :
//  *  - λ ≥ 1 ⇒ r ∈ [0,1). Cas λ=1 : r=0, alors r^0=1 et r^d=0 si d>0. La variation Δ reste bien définie.
//  *  - dpow_double(base, exp) gère exp≥0 avec exponentiation rapide O(log exp).
//  *  - 0^0 renvoyé à 1 par convention (exp<=0 → 1.0) pour la cohérence avec la limite combinatoire.
//  *
//  *  Complexité :
//  *  - O(#lambdas) par toggle ; pas d’accès à la liste d’adjacence ni de parcours des voisins.
//  *
//  *  Entrées (côté R via InitErgmTerm.cliques_GW) :
//  *  - inputs = c(rbind(lambda, r)) de longueur 2*N_CHANGE_STATS.
//  *  - coef.names vectorisés, emptynwstats = 0.
//  *
//  *  Sortie :
//  *  - CHANGE_STAT[j] += Δ_j pour chaque λ_j.
//  *
//  *  Validation :
//  *  - Cohérence directe avec la forme fermée et tests unitaires : ajout/retrait donnent
//  *    Δ = λ( r^{d_old} - r^{d_new} ).
//  */

// #include <math.h>
// #include <R_ext/Print.h>
// #include "ergm_changestat.h"
// #include "ergm_storage.h"

// // Si 1, affiche des traces dans la console R pendant summary()/MCMC.
// #define DEBUG_CLIQUES_GW 0

// // -----------------------------------------------------------------------------
// // Exponentiation rapide (exponentiation by squaring) pour exp >= 0.
// // Retourne 1.0 si exp <= 0 (convention 0^0 = 1 pour la cohérence combinatoire).
// // -----------------------------------------------------------------------------
// static inline double dpow_double(double base, int exp){
//   if(exp <= 0) return 1.0;
//   double r = 1.0, b = base;
//   int e = exp;
//   while(e){
//     if(e & 1) r *= b;
//     b *= b;
//     e >>= 1;
//   }
//   return r;
// }

// C_CHANGESTAT_FN(c_cliques_GW){
//   // --- 1) Réinitialiser le tampon de sortie pour CE toggle. ---
//   ZERO_ALL_CHANGESTATS(0);  // memset(CHANGE_STAT, 0, nstats*sizeof(double))

//   // --- 2) Taille du mode 1 (acteurs). Si non-biparti, BIPARTITE vaut 0. ---
//   const int n1 = BIPARTITE;

//   // --- 3) Sélection du nœud côté groupes (mode 2) : indices > n1. ---
//   //     Le toggle peut être orienté dans un sens quelconque : on cible le sommet > n1.
//   Vertex v2 = (tail > (Vertex)n1) ? tail : head;

//   // --- 4) Degrés avant/après toggle (calcul local). ---
//   //     OUT_DEG/IN_DEG : tableaux internes des degrés sortants/entrants.
//   int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);
//   int delta   = edgestate ? -1 : +1;   // retrait:-1, ajout:+1
//   int deg_new = deg_old + delta;

//   // --- 5) Traces éventuelles. ---
//   #if DEBUG_CLIQUES_GW
//     Rprintf("[c_cliques_GW] v2=%d | edgestate=%d | deg_old=%d -> deg_new=%d | nstats=%d\n",
//             (int)v2, (int)edgestate, deg_old, deg_new, N_CHANGE_STATS);
//   #endif

//   // --- 6) Boucle sur les sous-termes vectorisés (plusieurs λ_j possibles). ---
//   // INPUT_PARAM : [λ_0, r_0, λ_1, r_1, ..., λ_{J-1}, r_{J-1}], r_j=(λ_j-1)/λ_j.
//   for(int j = 0; j < N_CHANGE_STATS; ++j){
//     double lambda = INPUT_PARAM[2*j + 0];
//     double r      = INPUT_PARAM[2*j + 1];  // ∈ [0,1) si λ>=1

//     // Δ_j = λ ( r^{deg_old} - r^{deg_new} ).
//     double term_old = dpow_double(r, deg_old);
//     double term_new = dpow_double(r, deg_new);
//     double d = lambda * (term_old - term_new);

//     CHANGE_STAT[j] += d;

//     #if DEBUG_CLIQUES_GW
//       Rprintf("  j=%d | lambda=%.9g r=%.9g | rold=%.9g rnew=%.9g | Δ=%.9g\n",
//               j, lambda, r, term_old, term_new, d);
//     #endif
//   }
// }
