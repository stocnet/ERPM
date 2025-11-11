/**
 * @file    changestat_cov_ingroup.c
 * @brief   Change statistic {ergm} pour le terme ERPM `cov_ingroup`
 *
 * Définition :
 *   T = sum_{g in mode-2} [ n_g * (sum_{i in g} x_i) * 1[n_g in S] ].
 *
 * Variation locale pour un toggle (acteur i, groupe v2) :
 *   Soit X = sum_{j in g} x_j et n = deg(g) AVANT le toggle.
 *   AJOUT    : n' = n+1, X' = X + x_i  ->  Δ = (n'+0)*X' * 1[n' in S] - n*X * 1[n in S]
 *   RETRAIT  : n' = n-1, X' = X - x_i  ->  Δ = (n'+0)*X' * 1[n' in S] - n*X * 1[n in S]
 *
 * Hypothèses :
 *   - Réseau biparti : mode 1 = acteurs (1..BIPARTITE), mode 2 = groupes (>BIPARTITE).
 *   - Un toggle connecte exactement un acteur (v1<=n1) et un groupe (v2>n1).
 *   - On somme X par parcours local des voisins acteur du groupe v2.
 *
 * Inputs (depuis InitErgmTerm.cov_ingroup) dans INPUT_PARAM :
 *   idx 0 : n1   (double -> int)                nombre d'acteurs (mode 1)
 *   idx 1 : L    (double -> int)                taille de S (0 => toutes tailles)
 *   idx 2..2+L-1: sizes[k] (double -> int)      valeurs de S
 *   idx 2+L..    : x[1..n1] (double)            attributs des acteurs
 *
 * Sortie :
 *   CHANGE_STAT[0] += Δ.
 */

#include <math.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"

static inline int in_sizes_set(int n, const double *in, int L){
  if(L<=0) return 1; /* pas de filtre: toutes tailles acceptées */
  /* in[0..L-1] contient les tailles autorisées (stockées en double) */
  for(int k=0; k<L; ++k) if((int)in[k] == n) return 1;
  return 0;
}

C_CHANGESTAT_FN(c_cov_ingroup){
  ZERO_ALL_CHANGESTATS(0);

  /* ---------------- Lecture des inputs ---------------- */
  const int n1 = (int)INPUT_PARAM[0];
  const int L  = (int)INPUT_PARAM[1];

  const double *sizes = (L>0) ? (&INPUT_PARAM[2]) : NULL;
  const double *x     = (L>0) ? (&INPUT_PARAM[2+L]) : (&INPUT_PARAM[2]);

  /* ---------------- Identifier acteur v1 et groupe v2 ---------------- */
  const Vertex n1_lim = (Vertex)BIPARTITE; /* doit coïncider avec n1 */
  Vertex v2 = (tail > n1_lim) ? tail : head; /* groupe côté mode 2 */
  Vertex v1 = (tail > n1_lim) ? head : tail; /* acteur côté mode 1 */

  /* Sécurité : v1 doit être dans [1..n1] */
  if(v1 < (Vertex)1 || v1 > (Vertex)n1){
    CHANGE_STAT[0] += 0.0;
    return;
  }

  /* ---------------- Degré et somme d'attributs X pour le groupe v2 ---------------- */
  int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);

  /* Somme X = sum_{i in g} x_i en parcourant les voisins acteurs du groupe v2.
     Itérer sur les OUTEDGES de v2 suffit pour les réseaux non-orientés/bipartis. */
  double X = 0.0;
  Vertex h;
  Edge e;

  STEP_THROUGH_OUTEDGES(v2, e, h){
    if(h >= (Vertex)1 && h <= (Vertex)n1) X += x[(int)h - 1]; /* x indexé 0..n1-1, vertices 1..n */
  }
  /* Si des arêtes pouvaient être "entrantes" pour v2 (selon orientation),
     on pourrait dé-commenter le bloc suivant pour être exhaustif :
  STEP_THROUGH_INEDGES(v2, e, h){
    if(h >= (Vertex)1 && h <= (Vertex)n1) X += x[(int)h - 1];
  }
  */

  /* ---------------- Valeurs avant/après ---------------- */
  const int is_add   = (edgestate==0);             /* 1 si AJOUT, 0 si RETRAIT */
  const int n_new    = deg_old + (is_add ? +1 : -1);
  const double xi    = x[(int)v1 - 1];
  const double X_new = X + (is_add ? +xi : -xi);

  /* ---------------- Filtre de tailles S ---------------- */
  const int w_old = in_sizes_set(deg_old, sizes, L);
  const int w_new = in_sizes_set(n_new,  sizes, L);

  /* ---------------- Variation locale ---------------- */
  /* Δ = (n_new * X_new * w_new) - (deg_old * X * w_old) */
  double d = 0.0;
  d = (w_new ? ( (double)n_new * X_new ) : 0.0)
    - (w_old ? ( (double)deg_old * X   ) : 0.0);

  CHANGE_STAT[0] += d;
}