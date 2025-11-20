// ============================================================================
// Fichier : changestat_dyadcov_GW.c
// Terme   : dyadcov_GW
// Stat    :
//
//   T_GW(p;Z,lambda)
//     = sum_g sum_{k=2..n_g} a_k(lambda) * S_g^{(k)}(Z)
//     = sum_g sum_{k=2..n_g} (-1/lambda)^{k-2}
//                    * sum_{C∈C_k(g)} prod_{i<j∈C} z_{ij}
//
//   avec C_k(g) l’ensemble des sous-ensembles C⊂A(g) de taille k,
//   et S_g^{(k)}(Z) = sum_{C∈C_k(g)} prod_{i<j∈C} z_{ij}.
//
//   NB : a_k(lambda) = (-1/lambda)^{k-2} pour k>=2, ce qui donne pour lambda=2 :
//        a_2 = 1, a_3 = -1/2, a_4 = 1/4, ...
//
// INPUT_PARAM : c(n1, lambda, Z[n1*n1])  (Z en ordre colonne-major R)
// ============================================================================

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>

#define DEBUG_DYADCOV_GW 0
#define UNUSED_WARNING(x) (void)x

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
 * group_dyadcov_GW()
 *
 * Calcule, pour un groupe g donné :
 *   - n_g : nombre d'acteurs reliés à g (via IN+OUT, dédoublonné)
 *   - S_g^{GW}(Z,lambda)
 *       = sum_{k=2..n_g} a_k(lambda) * S_g^{(k)}(Z)
 *       avec a_k(lambda) = (-1/lambda)^{k-2}, k>=2.
 *
 * Renvoie S_g^{GW}(Z,lambda) et, si n_g_out!=NULL, place n_g dans *n_g_out.
 * -------------------------------------------------------------------------- */
static double group_dyadcov_GW(Vertex g,
                               int n1,
                               double lambda,
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

  /* Si ng < 2, aucune clique possible -> contribution nulle */
  if(ng < 2){
#if DEBUG_DYADCOV_GW
    Rprintf("[dyadcov_GW][group_dyadcov_GW] g=%d ng=%d < 2 -> 0\n",
            (int)g, ng);
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

  /* Somme géométriquement pondérée sur k = 2..ng */
  double sum_gw = 0.0;
  double factor = 1.0;          /* a_2(lambda) = 1 = (-1/lambda)^{0} */

  for(int k = 2; k <= ng; k++){
    double S_k = sum_cliques_k(actors, ng, k, n1, Z);
    sum_gw += factor * S_k;
    factor *= (-1.0 / lambda); /* a_{k+1} = a_k * (-1/lambda) */
#if DEBUG_DYADCOV_GW
    Rprintf("[dyadcov_GW][group_dyadcov_GW] g=%d ng=%d k=%d S_k=%g factor=%g\n",
            (int)g, ng, k, S_k, factor);
#endif
  }

  Free(actors);
  Free(seen);

#if DEBUG_DYADCOV_GW
  Rprintf("[dyadcov_GW][group_dyadcov_GW] g=%d ng=%d -> sum_gw=%g\n",
          (int)g, ng, sum_gw);
#endif

  return sum_gw;
}

/* ----------------------------------------------------------------------------
 * Change statistic un-toggle pour dyadcov_GW
 *
 * INPUT_PARAM layout:
 * [0]      = n1
 * [1]      = lambda
 * [2..]    = Z[n1*n1] (colonne-major)
 * -------------------------------------------------------------------------- */
C_CHANGESTAT_FN(c_dyadcov_GW){
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  const double *ip   = INPUT_PARAM;
  const int     n1   = (int)ip[0];
  const double  lambda = ip[1];
  const double *Z    = ip + 2;

#if DEBUG_DYADCOV_GW
  Rprintf("[dyadcov_GW] n1=%d lambda=%g\n", n1, lambda);
#endif

  if(lambda == 0.0){
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
  double S_before = group_dyadcov_GW(group, n1, lambda, Z, nwp, &ng_before);

  TOGGLE(a, b);  /* toggle virtuel */

  double S_after  = group_dyadcov_GW(group, n1, lambda, Z, nwp, &ng_after);

  TOGGLE(a, b);  /* annuler le toggle */

  double delta = S_after - S_before;
  CHANGE_STAT[0] += delta;

#if DEBUG_DYADCOV_GW
  Rprintf("[dyadcov_GW] a=%d b=%d group=%d ng_before=%d ng_after=%d "
          "S_before=%g S_after=%g delta=%g\n",
          (int)a, (int)b, (int)group,
          ng_before, ng_after, S_before, S_after, delta);
#endif
}
