// =============================================================================
//  Fichier : changestat_cliques.c
//  Terme   : `cliques` — Somme des k-cliques d'acteurs via tailles des groupes
//
//  k >= 2 : Stat = sum_g C(n_g, k).
//  k == 1 : Stat = # { g : n_g == 1 }  (compte des groupes de taille 1).
//
//  Variation un-toggle :
//    k >= 2 :
//      +edge : Δ =  C(n_g,   k-1)
//      -edge : Δ = -C(n_g-1, k-1)
//    k == 1 :
//      +edge : Δ = (n_g==0) ? +1 : (n_g==1) ? -1 : 0
//      -edge : Δ = (n_g==2) ? +1 : (n_g==1) ? -1 : 0
//
//  Option d'échelle : scale_j = 1 ou 1/C(N1,k) (marche aussi pour k=1, denom = N1).
// =============================================================================

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>

#define DEBUG_CLIQUES 0

C_CHANGESTAT_FN(c_cliques){
  ZERO_ALL_CHANGESTATS(0);

  const int n1 = BIPARTITE; // taille du mode 1 (acteurs)

  // nœud "groupe" (mode 2 : index > n1)
  Vertex t = tail, h = head;
  Vertex v2 = (t > (Vertex)n1) ? t : h;

  // degré groupe avant toggle
  int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);

  const int is_add = edgestate ? 0 : 1;

  for(int j = 0; j < N_CHANGE_STATS; ++j){
    int    k     = (int)   INPUT_PARAM[2*j + 0];
    double scale = (double)INPUT_PARAM[2*j + 1];

    double delta = 0.0;

    if(k == 1){
      if(is_add){
        // 0->1 : +1 ; 1->2 : -1 ; sinon 0
        if(deg_old == 0) delta = +1.0;
        else if(deg_old == 1) delta = -1.0;
      }else{
        // 2->1 : +1 ; 1->0 : -1 ; sinon 0
        if(deg_old == 2) delta = +1.0;
        else if(deg_old == 1) delta = -1.0;
      }
    }else{
      if(is_add){
        if(deg_old >= k-1) delta = CHOOSE(deg_old, k-1);
      }else{
        if(deg_old-1 >= k-1 && deg_old >= 1) delta = -CHOOSE(deg_old - 1, k-1);
        else delta = 0.0;
      }
    }

    CHANGE_STAT[j] += scale * delta;

    #if DEBUG_CLIQUES
      Rprintf("[c_cliques] v2=%d, deg_old=%d, k=%d, add=%d, delta=%.6f, scale=%g, out=%.6f\n",
              (int)v2, deg_old, k, is_add, delta, scale, CHANGE_STAT[j]);
    #endif
  }
}
