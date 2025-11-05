/**
 * @file    changestat_cliques_GW.c
 * @brief   Change statistic {ergm} pour le terme ERPM `cliques_GW` (forme un-toggle)
 *
 * @details
 *  Principe statistique :
 *  - Pour un groupe de taille d, la contribution agrégée des k-cliques pondérées
 *    géométriquement vaut S(d,λ) = λ[1 - ((λ-1)/λ)^d].
 *  - Lors d’un toggle affectant un unique groupe v2, la variation locale est :
 *        Δ = S(d_new,λ) - S(d_old,λ) = λ( r^{d_old} - r^{d_new} ),  r=(λ-1)/λ.
 *
 *  Implémentation {ergm} — forme "un-toggle" :
 *  - Signature imposée : C_CHANGESTAT_FN(c_cliques_GW)(tail, head, mtp, nwp, edgestate)
 *  - Chaque appel traite UN toggle (tail, head) dans l’état courant du réseau.
 *  - edgestate = 1  → l’arête existe (toggle = RETRAIT) → d_new = d_old - 1.
 *    edgestate = 0  → l’arête n’existe pas (toggle = AJOUT) → d_new = d_old + 1.
 *  - ZERO_ALL_CHANGESTATS(0) : remet le tampon de sortie à zéro.
 *    {ergm} additionne ensuite les contributions de tous les toggles.
 *
 *  Hypothèses et conventions internes :
 *  - Réseau biparti : mode 1 = acteurs (1..BIPARTITE), mode 2 = groupes (BIPARTITE+1..N).
 *  - Le toggle connecte toujours un acteur et un groupe → un seul groupe impacté par appel.
 *  - Le "degré" d’un groupe est le degré total OUT_DEG[v2] + IN_DEG[v2] (représentation interne).
 *  - Vectorisation sur j : le terme peut exposer plusieurs λ_j. Les entrées INPUT_PARAM sont
 *    empaquetées par colonnes : [λ_0, r_0, λ_1, r_1, ..., λ_{J-1}, r_{J-1}], avec r_j = (λ_j-1)/λ_j.
 *
 *  Précisions numériques :
 *  - λ ≥ 1 ⇒ r ∈ [0,1). Cas λ=1 : r=0, alors r^0=1 et r^d=0 si d>0. La variation Δ reste bien définie.
 *  - dpow_double(base, exp) gère exp≥0 avec exponentiation rapide O(log exp).
 *  - 0^0 renvoyé à 1 par convention (exp<=0 → 1.0) pour la cohérence avec la limite combinatoire.
 *
 *  Complexité :
 *  - O(#lambdas) par toggle ; pas d’accès à la liste d’adjacence ni de parcours des voisins.
 *
 *  Entrées (côté R via InitErgmTerm.cliques_GW) :
 *  - inputs = c(rbind(lambda, r)) de longueur 2*N_CHANGE_STATS.
 *  - coef.names vectorisés, emptynwstats = 0.
 *
 *  Sortie :
 *  - CHANGE_STAT[j] += Δ_j pour chaque λ_j.
 *
 *  Validation :
 *  - Cohérence directe avec la forme fermée et tests unitaires : ajout/retrait donnent
 *    Δ = λ( r^{d_old} - r^{d_new} ).
 */

#include <math.h>
#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"

// Si 1, affiche des traces dans la console R pendant summary()/MCMC.
#define DEBUG_CLIQUES_GW 0

// -----------------------------------------------------------------------------
// Exponentiation rapide (exponentiation by squaring) pour exp >= 0.
// Retourne 1.0 si exp <= 0 (convention 0^0 = 1 pour la cohérence combinatoire).
// -----------------------------------------------------------------------------
static inline double dpow_double(double base, int exp){
  if(exp <= 0) return 1.0;
  double r = 1.0, b = base;
  int e = exp;
  while(e){
    if(e & 1) r *= b;
    b *= b;
    e >>= 1;
  }
  return r;
}

C_CHANGESTAT_FN(c_cliques_GW){
  // --- 1) Réinitialiser le tampon de sortie pour CE toggle. ---
  ZERO_ALL_CHANGESTATS(0);  // memset(CHANGE_STAT, 0, nstats*sizeof(double))

  // --- 2) Taille du mode 1 (acteurs). Si non-biparti, BIPARTITE vaut 0. ---
  const int n1 = BIPARTITE;

  // --- 3) Sélection du nœud côté groupes (mode 2) : indices > n1. ---
  //     Le toggle peut être orienté dans un sens quelconque : on cible le sommet > n1.
  Vertex v2 = (tail > n1) ? tail : head;

  // --- 4) Degrés avant/après toggle (calcul local). ---
  //     OUT_DEG/IN_DEG : tableaux internes des degrés sortants/entrants.
  int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);
  int delta   = edgestate ? -1 : +1;   // retrait:-1, ajout:+1
  int deg_new = deg_old + delta;

  // --- 5) Traces éventuelles. ---
  #if DEBUG_CLIQUES_GW
    Rprintf("[c_cliques_GW] v2=%d | edgestate=%d | deg_old=%d -> deg_new=%d | nstats=%d\n",
            (int)v2, (int)edgestate, deg_old, deg_new, N_CHANGE_STATS);
  #endif

  // --- 6) Boucle sur les sous-termes vectorisés (plusieurs λ_j possibles). ---
  // INPUT_PARAM : [λ_0, r_0, λ_1, r_1, ..., λ_{J-1}, r_{J-1}], r_j=(λ_j-1)/λ_j.
  for(int j = 0; j < N_CHANGE_STATS; ++j){
    double lambda = INPUT_PARAM[2*j + 0];
    double r      = INPUT_PARAM[2*j + 1];  // ∈ [0,1) si λ>=1

    // Δ_j = λ ( r^{deg_old} - r^{deg_new} ).
    double term_old = dpow_double(r, deg_old);
    double term_new = dpow_double(r, deg_new);
    double d = lambda * (term_old - term_new);

    CHANGE_STAT[j] += d;

    #if DEBUG_CLIQUES_GW
      Rprintf("  j=%d | lambda=%.9g r=%.9g | rold=%.9g rnew=%.9g | Δ=%.9g\n",
              j, lambda, r, term_old, term_new, d);
    #endif
  }
}