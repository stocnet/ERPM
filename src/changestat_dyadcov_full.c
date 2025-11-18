// ============================================================================
// Fichier : changestat_dyadcov_full.c
// Terme   : dyadcov_full
// Stat    : T = sum_g 1[n_g in S] * sum_{i<j, i,j in g} z_{ij}
//           avec S = ensemble de tailles autorisées (optionnel)
// INPUT   : c(n1, L, sizes[L], Z[n1*n1])  (Z en ordre colonne-major R)
// ============================================================================

#include "ergm_changestat.h"
#include "ergm_storage.h"      // pour Calloc/Free
#include <R_ext/Print.h>

#define DEBUG_DYADCOV_FULL 0   /* mettre à 0 pour couper le debug */
#define UNUSED_WARNING(x) (void)x

/* ----------------------------------------------------------------------------
 * Vérifie si une taille n est incluse dans l'ensemble des tailles autorisées S.
 * Si L==0, toutes les tailles sont acceptées.
 * -------------------------------------------------------------------------- */
static inline int in_sizes(int n, int L, const double *sizes){
  if(L == 0) return 1;
  for(int i = 0; i < L; i++){
    if((int)sizes[i] == n) return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------------
 * Fonction locale group_dyadcov()
 * Évalue la contribution d’un groupe g :
 *   - identifie les acteurs reliés à g (via IN+OUT, dédoublonné) ;
 *   - calcule n_g = nombre d’acteurs dans g ;
 *   - si n_g >= 2 et n_g ∈ S, retourne:
 *       sum_{i<j, i,j in g} z_{ij}
 *     où Z est une matrice n1×n1 (ordre colonne-major comme en R).
 *
 * Z[ (j-1)*n1 + (i-1) ] = z_{ij}, i,j en {1,...,n1}
 * -------------------------------------------------------------------------- */
static double group_dyadcov(Vertex g,
                            int n1, int L, const double *sizes,
                            const double *Z,
                            Network *nwp){

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

  /* Groupe vide/singleton ou hors filtre de tailles -> contribution nulle */
  if(ng <= 1 || !in_sizes(ng, L, sizes)){
    #if DEBUG_DYADCOV_FULL
      Rprintf("[dyadcov_full][group_dyadcov] g=%d ng=%d -> 0 (vide/singleton/hors S)\n",
              (int)g, ng);
    #endif
    Free(seen);
    return 0.0;
  }

  /* Collecter les indices d’acteurs (1..n1) appartenant à g */
  int *actors = (int*)Calloc(ng, int);
  int k = 0;
  for(int i = 0; i < n1; i++){
    if(seen[i]){
      actors[k++] = i + 1; /* Vertex index 1-based */
    }
  }

  double sum = 0.0;

  /* Somme sur toutes les dyades i<j à l’intérieur du groupe g */
  for(int p = 0; p < ng; p++){
    int i = actors[p];             /* 1..n1 */
    int row = i - 1;               /* 0..n1-1 */

    for(int q = p + 1; q < ng; q++){
      int j = actors[q];           /* 1..n1 */
      int col = j - 1;             /* 0..n1-1 */

      /* Z en ordre colonne-major (style R) : Z[(j-1)*n1 + (i-1)] = z_ij */
      int idx = col * n1 + row;
      sum += Z[idx];
    }
  }

  #if DEBUG_DYADCOV_FULL
    Rprintf("[dyadcov_full][group_dyadcov] g=%d ng=%d -> sum=%g\n",
            (int)g, ng, sum);
  #endif

  Free(actors);
  Free(seen);

  return sum;
}

/* ----------------------------------------------------------------------------
 * Change statistic un-toggle pour dyadcov_full
 *
 * INPUT_PARAM layout:
 * [0]      = n1
 * [1]      = L
 * [2..1+L] = sizes[L]
 * [2+L..]  = Z[n1*n1] (colonne-major)
 * -------------------------------------------------------------------------- */
C_CHANGESTAT_FN(c_dyadcov_full){
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  const double *ip    = INPUT_PARAM;
  const int n1        = (int)ip[0];
  const int L         = (int)ip[1];
  const double *sizes = ip + 2;
  const double *Z     = ip + 2 + L;

  #if DEBUG_DYADCOV_FULL
    Rprintf("[dyadcov_full] n1=%d L=%d\n", n1, L);
  #endif

  /* Identifier les sommets acteur et groupe du toggle courant */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;
  Vertex group = (a <= (Vertex)n1) ? b : a;
  UNUSED_WARNING(actor);

  /* Contribution du groupe avant / après basculement virtuel */
  double F_before = group_dyadcov(group, n1, L, sizes, Z, nwp);

  TOGGLE(a, b);  /* toggle virtuel */

  double F_after  = group_dyadcov(group, n1, L, sizes, Z, nwp);

  TOGGLE(a, b);  /* annuler le toggle */

  CHANGE_STAT[0] += (F_after - F_before);

  #if DEBUG_DYADCOV_FULL
    Rprintf("[dyadcov_full] a=%d b=%d group=%d before=%g after=%g delta=%g\n",
            (int)a, (int)b, (int)group, F_before, F_after, (F_after - F_before));
  #endif
}