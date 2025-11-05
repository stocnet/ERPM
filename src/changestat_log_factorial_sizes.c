/**
 * @file    changestat_log_factorial_sizes.c
 * @brief   Change statistic {ergm} pour le terme ERPM `log_factorial_sizes` (forme un-toggle, non vectorisé)
 *
 * @details
 *  Définition :
 *    Stat = sum_{g ∈ mode-2} lgamma(deg(g)), avec convention f(0)=0.
 *
 *  Variations locales (toggle touchant un seul groupe v2) :
 *    - AJOUT  : deg_old = n → deg_new = n+1
 *               Δ = lgamma(n+1) - lgamma(n) = log(n)    si n ≥ 1 ; Δ = 0 si n = 0
 *    - RETRAIT: deg_old = n → deg_new = n-1
 *               Δ = lgamma(n-1) - lgamma(n) = -log(n-1) si n ≥ 2 ; Δ = 0 si n = 1
 *
 *  Implémentation {ergm} — forme "un-toggle" :
 *    - Signature imposée : C_CHANGESTAT_FN(c_log_factorial_sizes)(tail, head, mtp, nwp, edgestate)
 *    - Chaque appel traite UN toggle (tail, head) dans l’état courant du réseau.
 *    - edgestate = 1  → l’arête existe (toggle = RETRAIT).
 *      edgestate = 0  → l’arête n’existe pas (toggle = AJOUT).
 *    - ZERO_ALL_CHANGESTATS(0) réinitialise le tampon de sortie pour ce toggle.
 *
 *  Hypothèses et conventions :
 *    - Réseau biparti : mode 1 = acteurs (1..BIPARTITE), mode 2 = groupes (BIPARTITE+1..N).
 *    - Le toggle connecte un acteur et un groupe → un seul groupe impacté par appel.
 *    - Le degré du groupe est le degré total OUT_DEG[v2] + IN_DEG[v2] (représentation interne).
 *
 *  Complexité :
 *    - O(1) par toggle. Aucun parcours de voisins.
 *
 *  Entrées (depuis R via InitErgmTerm.log_factorial_sizes) :
 *    - Non vectorisé : N_CHANGE_STATS == 1, pas d’INPUT_PARAM.
 *
 *  Sortie :
 *    - CHANGE_STAT[0] += Δ.
 */

#include <math.h>
#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"

// Si 1, traces de debug pendant summary()/MCMC.
#define DEBUG_LOG_FACTORIAL 0

C_CHANGESTAT_FN(c_log_factorial_sizes){
  // --- 1) Réinitialiser le tampon de sortie pour CE toggle. ---
  ZERO_ALL_CHANGESTATS(0);

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

  #if DEBUG_LOG_FACTORIAL
    Rprintf("[c_log_factorial_sizes] v2=%d | edgestate=%d | deg_old=%d -> deg_new=%d\n",
            (int)v2, (int)edgestate, deg_old, deg_new);
  #endif

  // --- 5) Variation locale Δ selon les cas (sans appeler lgamma pour les bords). ---
  double d = 0.0;
  if(edgestate == 0){
    // AJOUT: Δ = log(deg_old) si deg_old >= 1 ; sinon 0.
    if(deg_old >= 1) d = log((double)deg_old);
  }else{
    // RETRAIT: Δ = -log(deg_old - 1) si deg_old >= 2 ; sinon 0.
    if(deg_old >= 2) d = -log((double)(deg_old - 1));
  }

  // --- 6) Accumulation (terme non vectorisé => N_CHANGE_STATS == 1). ---
  CHANGE_STAT[0] += d;

  #if DEBUG_LOG_FACTORIAL
    Rprintf("  Δ=%.9g\n", d);
  #endif
}