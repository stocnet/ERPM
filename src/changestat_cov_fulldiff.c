// ============================================================================
// Fichier : changestat_cov_fulldiff.c
// Terme   : cov_fulldiff
// Stat    : T = sum_g 1[n_g in S] * (x_g^max - x_g^min)
//           avec S = ensemble de tailles autorisées (optionnel)
// INPUT   : c(n1, L, sizes[L], x[1..n1])
// ============================================================================

#include "ergm_changestat.h"
#include "ergm_storage.h"      // pour Calloc/Free
#include <R_ext/Print.h>

#define DEBUG_COV_FULLDIFF 0   /* mettre à 0 pour couper le debug */
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
 * Fonction locale group_range()
 * Évalue la contribution d’un groupe g :
 *   - calcule n_g = nombre d’acteurs reliés à g (robuste: IN+OUT, dédoublonné) ;
 *   - calcule x_g^min et x_g^max sur ces acteurs ;
 *   - retourne (x_g^max - x_g^min) si n_g>1 et n_g ∈ S, sinon 0.
 *
 * x[0..n1-1] : valeur numérique de la covariée pour chaque acteur (mode 1).
 * -------------------------------------------------------------------------- */
static double group_range(Vertex g,
                          int n1, int L, const double *sizes,
                          const double *x,
                          Network *nwp){

  int ng = 0;
  double xmin = 0.0, xmax = 0.0;
  unsigned char first = 1;

  /* Marquage des acteurs déjà comptés pour éviter les doubles comptes
     (réseaux dirigés, arêtes réciproques, etc.). */
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char); /* Calloc => 0 */

  Vertex h;
  Edge e;

  /* OUT-neighbors g -> h */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        ng++;
        double val = x[idx];
        if(first){
          xmin = xmax = val;
          first = 0;
        }else{
          if(val < xmin) xmin = val;
          if(val > xmax) xmax = val;
        }
      }
    }
  }

  /* IN-neighbors h -> g */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        ng++;
        double val = x[idx];
        if(first){
          xmin = xmax = val;
          first = 0;
        }else{
          if(val < xmin) xmin = val;
          if(val > xmax) xmax = val;
        }
      }
    }
  }

  /* Nettoyage du tableau seen */
  STEP_THROUGH_OUTEDGES(g, e, h){ if(h <= (Vertex)n1) seen[(int)h-1] = 0; }
  STEP_THROUGH_INEDGES(g, e, h){ if(h <= (Vertex)n1) seen[(int)h-1] = 0; }
  Free(seen);

  /* Groupe vide ou singleton : pas de dispersion interne */
  if(ng <= 1){
    #if DEBUG_COV_FULLDIFF
      Rprintf("[cov_fulldiff][group_range] g=%d ng=%d -> 0 (vide/singleton)\n",
              (int)g, ng);
    #endif
    return 0.0;
  }

  /* Filtre de tailles S */
  if(!in_sizes(ng, L, sizes)){
    #if DEBUG_COV_FULLDIFF
      Rprintf("[cov_fulldiff][group_range] g=%d ng=%d -> hors S\n", (int)g, ng);
    #endif
    return 0.0;
  }

  double res = xmax - xmin;

  #if DEBUG_COV_FULLDIFF
    Rprintf("[cov_fulldiff][group_range] g=%d ng=%d xmin=%g xmax=%g -> res=%g\n",
            (int)g, ng, xmin, xmax, res);
  #endif

  return res;
}

/* ----------------------------------------------------------------------------
 * Change statistic un-toggle pour cov_fulldiff
 *
 * INPUT_PARAM layout:
 * [0]      = n1
 * [1]      = L
 * [2..1+L] = sizes[L]
 * [2+L..]  = x[1..n1]
 * -------------------------------------------------------------------------- */
C_CHANGESTAT_FN(c_cov_fulldiff){
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  const double *ip    = INPUT_PARAM;
  const int n1        = (int)ip[0];
  const int L         = (int)ip[1];
  const double *sizes = ip + 2;
  const double *x     = ip + 2 + L;

  #if DEBUG_COV_FULLDIFF
    Rprintf("[cov_fulldiff] n1=%d L=%d\n", n1, L);
  #endif

  /* Identifier les sommets acteur et groupe du toggle courant */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;
  Vertex group = (a <= (Vertex)n1) ? b : a;
  UNUSED_WARNING(actor);

  /* Contribution du groupe avant/après basculement virtuel */
  double F_before = group_range(group, n1, L, sizes, x, nwp);

  TOGGLE(a, b);  /* toggle virtuel */

  double F_after  = group_range(group, n1, L, sizes, x, nwp);

  TOGGLE(a, b);  /* annuler le toggle */

  CHANGE_STAT[0] += (F_after - F_before);

  #if DEBUG_COV_FULLDIFF
    Rprintf("[cov_fulldiff] a=%d b=%d group=%d before=%g after=%g delta=%g\n",
            (int)a, (int)b, (int)group, F_before, F_after, (F_after-F_before));
  #endif
}
