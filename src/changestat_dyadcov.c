// ============================================================================
// Fichier : changestat_dyadcov.c
// Terme   : dyadcov
// Stat    :
//
//  - Non normalisé (normalized = 0) :
//      T^{(k)}(p;Z) = sum_g sum_{C∈C_k(g)} prod_{i<j∈C} z_{ij}
//
//  - Normalisé (normalized = 1) :
//      T^{(k)}_norm(p;Z) = sum_g 1[n_g>=k] * (1/choose(n_g,k)) * sum_{C∈C_k(g)} prod_{i<j∈C} z_{ij}
//
//  où C_k(g) est l’ensemble des sous-ensembles C⊂A(g) de taille k.
//
// INPUT_PARAM : c(n1, k, normalized, Z[n1*n1])  (Z en ordre colonne-major R)
// ============================================================================

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>

#define DEBUG_DYADCOV 0
#define UNUSED_WARNING(x) (void)x

/* ----------------------------------------------------------------------------
 * choose_double(n,k) : combinaison nCk en double (n>=0, 0<=k<=n)
 * Implémentation stable pour k modéré.
 * -------------------------------------------------------------------------- */
static double choose_double(int n, int k){
  if(k < 0 || k > n) return 0.0;
  if(k == 0 || k == n) return 1.0;
  if(k > n - k) k = n - k; /* symétrie */
  double res = 1.0;
  for(int i = 1; i <= k; i++){
    res *= (double)(n - k + i);
    res /= (double)i;
  }
  return res;
}

/* ----------------------------------------------------------------------------
 * sum_cliques_k()
 * Somme des produits P(C;Z) sur toutes les cliques C de taille k dans un
 * vecteur d'acteurs donné.
 *
 * actors : indices 1..n1 des acteurs dans le groupe (taille ng)
 * ng     : nb d'acteurs dans le groupe
 * k      : taille de clique ciblée (k>=2)
 * Z      : matrice n1 x n1 en ordre colonne-major, Z[(j-1)*n1 + (i-1)] = z_ij
 * -------------------------------------------------------------------------- */
static double sum_cliques_k(const int *actors,
                            int ng, int k,
                            int n1,
                            const double *Z){

  if(k > ng) return 0.0;

  int *comb = (int*)Calloc(k, int);
  for(int i = 0; i < k; i++) comb[i] = i;

  double total = 0.0;

  while(1){
    /* Produit sur toutes les dyades i<j dans la clique courante */
    double prod = 1.0;
    for(int p = 0; p < k; p++){
      int idx_i = actors[ comb[p] ];   /* acteur (vertex) 1..n1 */
      int row   = idx_i - 1;           /* 0..n1-1 */

      for(int q = p + 1; q < k; q++){
        int idx_j = actors[ comb[q] ];
        int col   = idx_j - 1;
        int z_idx = col * n1 + row;    /* colonne-major */
        prod *= Z[z_idx];
      }
    }
    total += prod;

    /* Générer la combinaison suivante (lexicographique) */
    int pos = k - 1;
    while(pos >= 0 && comb[pos] == (ng - k + pos)) pos--;
    if(pos < 0) break; /* plus de combinaisons */

    comb[pos]++;
    for(int j = pos + 1; j < k; j++){
      comb[j] = comb[j - 1] + 1;
    }
  }

  Free(comb);
  return total;
}

/* ----------------------------------------------------------------------------
 * group_dyadcov_k()
 *
 * Calcule, pour un groupe g donné :
 *   - n_g : nombre d'acteurs reliés à g (via IN+OUT, dedoublonné)
 *   - S_g^{(k)}(Z) = sum_{C∈C_k(g)} prod_{i<j∈C} z_{ij}
 *
 * Renvoie S_g^{(k)}(Z) et, si n_g_out!=NULL, place n_g dans *n_g_out.
 * -------------------------------------------------------------------------- */
static double group_dyadcov_k(Vertex g,
                              int n1, int k,
                              const double *Z,
                              Network *nwp,
                              int *n_g_out){

  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char); /* init à 0 */
  Vertex h;
  Edge e;

  /* Marquage des acteurs voisins (OUT) */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      seen[idx] = 1;
    }
  }

  /* Marquage des acteurs voisins (IN) */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      seen[idx] = 1;
    }
  }

  /* Compter les acteurs dans le groupe g */
  int ng = 0;
  for(int i = 0; i < n1; i++){
    if(seen[i]) ng++;
  }
  if(n_g_out) *n_g_out = ng;

  /* Si ng < k, aucune clique possible -> contribution nulle */
  if(ng < k){
#if DEBUG_DYADCOV
    Rprintf("[dyadcov][group_dyadcov_k] g=%d ng=%d < k=%d -> 0\n",
            (int)g, ng, k);
#endif
    Free(seen);
    return 0.0;
  }

  /* Collecter les indices d’acteurs (1..n1) appartenant à g */
  int *actors = (int*)Calloc(ng, int);
  int idx = 0;
  for(int i = 0; i < n1; i++){
    if(seen[i]){
      actors[idx++] = i + 1;  /* vertex index 1-based */
    }
  }

  double sum = sum_cliques_k(actors, ng, k, n1, Z);

#if DEBUG_DYADCOV
  Rprintf("[dyadcov][group_dyadcov_k] g=%d ng=%d k=%d -> sum=%g\n",
          (int)g, ng, k, sum);
#endif

  Free(actors);
  Free(seen);

  return sum;
}

/* ----------------------------------------------------------------------------
 * Change statistic un-toggle pour dyadcov
 *
 * INPUT_PARAM layout:
 * [0]      = n1
 * [1]      = k
 * [2]      = normalized (0 ou 1)
 * [3..]    = Z[n1*n1] (colonne-major)
 * -------------------------------------------------------------------------- */
C_CHANGESTAT_FN(c_dyadcov){
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  const double *ip        = INPUT_PARAM;
  const int     n1        = (int)ip[0];
  int           k         = (int)ip[1];
  const int     normalized= (int)ip[2];
  const double *Z         = ip + 3;

#if DEBUG_DYADCOV
  Rprintf("[dyadcov] n1=%d k=%d normalized=%d\n", n1, k, normalized);
#endif

  if(k < 2){
    /* Sécurité : k<2 -> stat nulle */
    CHANGE_STAT[0] = 0.0;
    return;
  }

  /* Identifier les sommets acteur et groupe du toggle courant */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;
  Vertex group = (a <= (Vertex)n1) ? b : a;
  UNUSED_WARNING(actor);

  int ng_before = 0, ng_after = 0;

  /* Contribution du groupe avant / après basculement virtuel */
  double S_before = group_dyadcov_k(group, n1, k, Z, nwp, &ng_before);

  TOGGLE(a, b);  /* toggle virtuel */

  double S_after  = group_dyadcov_k(group, n1, k, Z, nwp, &ng_after);

  TOGGLE(a, b);  /* annuler le toggle */

  double delta;
  if(!normalized){
    /* Version brute : somme des poids de cliques */
    delta = S_after - S_before;
  } else {
    /* Version normalisée par choose(n_g,k) */
    double C_before = (ng_before >= k) ? choose_double(ng_before, k) : 0.0;
    double C_after  = (ng_after  >= k) ? choose_double(ng_after,  k) : 0.0;

    double T_before = (C_before > 0.0) ? (S_before / C_before) : 0.0;
    double T_after  = (C_after  > 0.0) ? (S_after  / C_after)  : 0.0;
    delta = T_after - T_before;
  }

  CHANGE_STAT[0] += delta;

#if DEBUG_DYADCOV
  Rprintf("[dyadcov] a=%d b=%d group=%d ng_before=%d ng_after=%d "
          "S_before=%g S_after=%g delta=%g\n",
          (int)a, (int)b, (int)group,
          ng_before, ng_after, S_before, S_after, delta);
#endif
}
