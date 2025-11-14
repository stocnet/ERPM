// ============================================================================
// Fichier : changestat_cov_diff_GW.c
// Terme   : cov_diff_GW
// Stat    : T_GW = sum_{k>=2} (-1/lambda)^{k-1} c_k(cov,p)
//           avec
//             c_k = sum_g sum_{S ⊂ g, |S|=k} (max_{i∈S} x_i - min_{i∈S} x_i)
// INPUT   : [0] = n1
//           [1] = L (nombre de valeurs de lambda)
//           [2..(1+L)]   = lambda[1..L]
//           [2+L..]      = x[1..n1] (covariée acteurs, 0-based côté C)
// Remarque : un-toggle, variation locale sur le seul groupe affecté par (tail, head).
// ============================================================================

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>
#include <math.h>

#define DEBUG_COV_DIFF_GW 0
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
 * group_covdiff_allk :
 *   Pour un groupe g :
 *     - extrait la liste des acteurs du groupe (mode 1 : 1..n1),
 *     - si n_g < 2 -> tous les c_k = 0,
 *     - sinon, pour chaque k = 2..n_g, calcule
 *         c_k(g) = sum_{S ⊂ g, |S|=k} D(S),
 *       en appelant sum_D_rec().
 *
 * Résultat :
 *   - remplit *ng_out avec n_g (nombre d'acteurs dans le groupe),
 *   - remplit ck[2..n_g] avec c_k(g) (ck[0] et ck[1] laissés à 0).
 *
 * Les tableaux idxs[] et ck[] sont alloués en amont et passés par l'appelant :
 *   - idxs doit être de taille >= n1 (nous n'utilisons que les ng premières cases),
 *   - ck   doit être de taille >= (n1+1) pour indexer directement ck[k].
 * -------------------------------------------------------------------------- */
static void group_covdiff_allk(Vertex g,
                               int n1,
                               const double *x,
                               Network *nwp,
                               int *ng_out,
                               double *ck,
                               int n1_for_ck,
                               int *idxs){

  Edge e;
  Vertex h;
  int ng = 0;

  /* marquage local des acteurs du groupe (0..n1-1) */
  unsigned char *seen = (unsigned char *)Calloc(n1, unsigned char);

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

  /* nettoyage du marquage */
  if(ng > 0){
    for(int i = 0; i < ng; i++){
      seen[ idxs[i] ] = 0;
    }
  }
  Free(seen);

  *ng_out = ng;

  if(ng < 2){
    /* aucun c_k non nul */
    return;
  }

  if(n1_for_ck < n1+1){
    /* garde-fou : le buffer ck doit avoir au moins n1+1 cases */
    return;
  }

  /* initialiser ck[k] pour k=2..ng à 0 (ck[0], ck[1] restent 0) */
  for(int k = 2; k <= ng; k++){
    ck[k] = 0.0;
  }

  /* pour chaque taille k, sommer D(S) sur toutes les combinaisons de taille k */
  for(int k = 2; k <= ng; k++){
    double sumD = 0.0;
    int *comb = (int *)Calloc(k, int);
    sum_D_rec(0, 0, k, ng, idxs, x, comb, &sumD);
    Free(comb);
    ck[k] = sumD;
  }

#if DEBUG_COV_DIFF_GW
  Rprintf("[cov_diff_GW][group_covdiff_allk] g=%d ng=%d : ", (int)g, ng);
  for(int k = 2; k <= ng; k++){
    Rprintf(" c_%d=%g", k, ck[k]);
  }
  Rprintf("\n");
#endif
}

/* ----------------------------------------------------------------------------
 * Change statistic un-toggle pour cov_diff_GW
 *
 * INPUT_PARAM layout:
 * [0]      = n1
 * [1]      = L (nb de lambda)
 * [2..1+L] = lambda[1..L]
 * [2+L.. ] = x[1..n1]
 *
 * N_CHANGE_STATS = L (une stat par valeur de lambda).
 * -------------------------------------------------------------------------- */
C_CHANGESTAT_FN(c_cov_diff_GW){
  /* initialiser toutes les stats à 0 pour ce toggle */
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  const double *ip    = INPUT_PARAM;
  const int n1        = (int)ip[0];
  const int L         = (int)ip[1];
  const double *lambda = ip + 2;
  const double *x     = ip + 2 + L;

#if DEBUG_COV_DIFF_GW
  Rprintf("[cov_diff_GW] n1=%d L=%d | lambda[1]=%g\n", n1, L, (L > 0 ? lambda[0] : NA_REAL));
#endif

  /* Identifier l’acteur et le groupe du toggle courant */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;   /* pas utilisé mais gardé pour debug/robustesse */
  Vertex group = (a <= (Vertex)n1) ? b : a;
  UNUSED_WARNING(actor);

  /* buffers pour c_k avant/après :
   *   - taille n1+1 pour adresser ck[k] directement.
   *   - idxs pour les indices d'acteurs du groupe (0..n1-1).
   */
  double *ck_before = (double *)Calloc(n1 + 1, double);
  double *ck_after  = (double *)Calloc(n1 + 1, double);
  int    *idxs      = (int    *)Calloc(n1,     int);

  int ng_before = 0, ng_after = 0;

  /* c_k(g^-) avant toggle */
  group_covdiff_allk(group, n1, x, nwp, &ng_before, ck_before, n1 + 1, idxs);

  /* appliquer le toggle virtuel */
  TOGGLE(a, b);

  /* c_k(g^+) après toggle */
  group_covdiff_allk(group, n1, x, nwp, &ng_after, ck_after, n1 + 1, idxs);

  /* annuler le toggle */
  TOGGLE(a, b);

  /* Kmax local pour ce groupe : max(ng_before, ng_after) */
  int Kmax = ng_before > ng_after ? ng_before : ng_after;
  if(Kmax < 2){
    /* aucun effet (groupe de taille < 2) */
    Free(ck_before);
    Free(ck_after);
    Free(idxs);
    return;
  }

  /* calcul de Delta T_GW pour chaque lambda_l :
   *   Delta T_GW(l) = sum_{k=2..Kmax} (-1/lambda_l)^{k-1} * (c_k^+ - c_k^-)
   */
  for(int l = 0; l < L; l++){
    double lam = lambda[l];
    double base = -1.0 / lam;
    double pow  = base;   /* correspond à k=2 -> exposant (k-1)=1 */
    double deltaT = 0.0;

    for(int k = 2; k <= Kmax; k++){
      double ckB = ck_before[k];
      double ckA = ck_after[k];
      double dck = ckA - ckB;
      if(dck != 0.0){
        deltaT += pow * dck;
      }
      pow *= base;  /* pour le k+1 suivant */
    }

    CHANGE_STAT[l] += deltaT;

#if DEBUG_COV_DIFF_GW
    Rprintf("[cov_diff_GW] a=%d b=%d group=%d lambda[%d]=%g -> DeltaT=%g\n",
            (int)a, (int)b, (int)group, l+1, lam, deltaT);
#endif
  }

  Free(ck_before);
  Free(ck_after);
  Free(idxs);
}
