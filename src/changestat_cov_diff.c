// ============================================================================
// Fichier : changestat_cov_diff.c
// Terme   : cov_diff
// Stat    : T_k = sum_g sum_{S ⊂ g, |S|=k} (max_{i∈S} x_i - min_{i∈S} x_i)
// Option  : normalized = 1  => moyenne par groupe (division par C(n_g, k))
// INPUT   : c(n1, k, norm_flag, x[1..n1])
// ============================================================================

#include "ergm_changestat.h"
#include "ergm_storage.h"      // Calloc / Free
#include <R_ext/Print.h>
#include <math.h>

#define DEBUG_COV_DIFF 0
#define UNUSED_WARNING(x) (void)(x)

/* ----------------------------------------------------------------------------
 * sum_D_rec :
 *   Énumère toutes les combinaisons de taille k parmi ng indices
 *   (référencés dans idxs[0..ng-1]) et additionne pour chacune la
 *   différence max-min de x.
 *
 * pos   : profondeur (0..k)
 * start : index de départ dans idxs
 * k,ng  : paramètres de taille
 * idxs  : indices des acteurs du groupe (0-based pour x[])
 * x     : vecteur de covariée
 * comb  : tableau de taille k stockant les positions choisies dans idxs[]
 * acc   : accumulateur de la somme des D(S)
 * -------------------------------------------------------------------------- */
static void sum_D_rec(int pos, int start,
                      int k, int ng,
                      const int *idxs,
                      const double *x,
                      int *comb,
                      double *acc){

  if(pos == k){
    /* On a une combinaison complète dans comb[0..k-1] */
    double xmin = 0.0, xmax = 0.0;
    int first = 1;
    for(int t = 0; t < k; t++){
      int idx = idxs[ comb[t] ];   /* indice 0-based dans x[] */
      double val = x[idx];
      if(first){
        xmin = xmax = val;
        first = 0;
      }else{
        if(val < xmin) xmin = val;
        if(val > xmax) xmax = val;
      }
    }
    *acc += (xmax - xmin);
    return;
  }

  /* Choix du prochain élément dans idxs[start..ng-1] */
  for(int i = start; i <= ng - (k - pos); i++){
    comb[pos] = i;
    sum_D_rec(pos + 1, i + 1, k, ng, idxs, x, comb, acc);
  }
}

/* ----------------------------------------------------------------------------
 * group_covdiff :
 *   Calcule la contribution T_k(g) pour un groupe g :
 *     - extrait la liste des acteurs du groupe ;
 *     - si n_g < k => 0 ;
 *     - sinon somme sur toutes les k-cliques S de g la quantité D(S) ;
 *     - si norm_flag == 1, divise par C(n_g, k).
 *
 * n1        : taille du mode acteurs
 * k         : taille de clique
 * norm_flag : 0 = brut, 1 = by_group
 * x         : vecteur covariée (longueur >= n1)
 * nwp       : pointeur réseau
 * -------------------------------------------------------------------------- */
static double group_covdiff(Vertex g,
                            int n1,
                            int k,
                            int norm_flag,
                            const double *x,
                            Network *nwp){

  Edge e;
  Vertex h;
  int ng = 0;

  /* Marquage des acteurs déjà comptés (pour être robuste à des directions) */
  unsigned char *seen = (unsigned char *)Calloc(n1, unsigned char);  /* 0 par défaut */
  int *idxs          = (int *)Calloc(n1, int);                        /* indices 0-based dans x[] */

  /* OUT-neighbors g -> h */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        idxs[ng++] = idx;
      }
    }
  }

  /* IN-neighbors h -> g */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        idxs[ng++] = idx;
      }
    }
  }

  /* Réinitialisation du marquage et libération de seen */
  if(ng > 0){
    for(int i = 0; i < ng; i++){
      seen[ idxs[i] ] = 0;
    }
  }
  Free(seen);

  if(ng < k){
    #if DEBUG_COV_DIFF
      Rprintf("[cov_diff][group_covdiff] g=%d ng=%d < k=%d -> 0\n",
              (int)g, ng, k);
    #endif
    Free(idxs);
    return 0.0;
  }

  /* Somme des D(S) sur toutes les combinaisons de taille k */
  double sumD = 0.0;
  int *comb = (int *)Calloc(k, int);
  sum_D_rec(0, 0, k, ng, idxs, x, comb, &sumD);
  Free(comb);

  double res = sumD;

  /* Normalisation par groupe (moyenne sur les C(n_g, k) cliques) */
  if(norm_flag){
    double denom = CHOOSE(ng, k);
    if(denom > 0.0){
      res /= denom;
    }else{
      res = 0.0;
    }
  }

  #if DEBUG_COV_DIFF
    Rprintf("[cov_diff][group_covdiff] g=%d ng=%d k=%d norm=%d -> res=%g\n",
            (int)g, ng, k, norm_flag, res);
  #endif

  Free(idxs);
  return res;
}

/* ----------------------------------------------------------------------------
 * Change statistic un-toggle pour cov_diff
 *
 * INPUT_PARAM layout:
 * [0]      = n1
 * [1]      = k
 * [2]      = norm_flag (0 = brut, 1 = by_group)
 * [3..]    = x[1..n1]
 * -------------------------------------------------------------------------- */
C_CHANGESTAT_FN(c_cov_diff){
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  const double *ip       = INPUT_PARAM;
  const int n1           = (int)ip[0];
  const int k            = (int)ip[1];
  const int norm_flag    = (int)ip[2];
  const double *x        = ip + 3;

  #if DEBUG_COV_DIFF
    Rprintf("[cov_diff] n1=%d k=%d norm_flag=%d\n", n1, k, norm_flag);
  #endif

  /* Identifier l’acteur et le groupe du toggle courant */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;
  Vertex group = (a <= (Vertex)n1) ? b : a;
  UNUSED_WARNING(actor);

  /* Contribution avant/après toggle virtuel */
  double F_before = group_covdiff(group, n1, k, norm_flag, x, nwp);

  TOGGLE(a, b);  /* toggle virtuel */

  double F_after  = group_covdiff(group, n1, k, norm_flag, x, nwp);

  TOGGLE(a, b);  /* annuler le toggle */

  CHANGE_STAT[0] += (F_after - F_before);

  #if DEBUG_COV_DIFF
    Rprintf("[cov_diff] a=%d b=%d group=%d before=%g after=%g delta=%g\n",
            (int)a, (int)b, (int)group, F_before, F_after, (F_after-F_before));
  #endif
}
