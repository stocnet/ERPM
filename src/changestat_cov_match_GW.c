// ============================================================================
// Fichier : src/changestat_cov_match_GW.c
// Terme   : cov_match_GW
// Stat    : S = sum_g sum_r λ * (1 - r_λ^{ n_{g,r} }), r_λ=(λ-1)/λ
// Ciblé   : S^(κ) = sum_g λ * (1 - r_λ^{ n_{g,κ} })
// Normes  : by_group -> Σ_g [ Num(g)/Den(g) ], global -> S / [λ(1-r_λ^{N1})]
// Toggle  : ajout  -> Δ =  r_λ^{ m }
//           retrait-> Δ = -r_λ^{ m-1 }   où m = n_{g,r*} avant toggle
// Complexité : O(deg(g)) pour histogramme; Δ en O(1).
// ============================================================================

#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <math.h>
#include <string.h>

#define DEBUG_COV_MATCH_GW 0
#define UNUSED_WARNING(x) (void)x

static inline int code_of_actor(Vertex i, const double *z_codes){
  return (int)z_codes[(size_t)(i-1)]; // acteurs indexés 1..n1
}

/* Acteurs voisins uniques du groupe g */
static int neighbors_actors_of_group(Network *nwp, Vertex g, Vertex *actors, int n1){
  int cnt = 0;
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char); // 0-inited
  Vertex h; Edge e;

  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){ seen[idx]=1; actors[cnt++]=h; }
    }
  }
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){ seen[idx]=1; actors[cnt++]=h; }
    }
  }
  Free(seen);
  return cnt;
}

/* Histogramme des codes de modalités dans le groupe */
static int histogram_codes(const Vertex *actors, int na, const double *z_codes, int *codes, int *counts){
  int m = 0;
  for(int a=0; a<na; ++a){
    int code = code_of_actor(actors[a], z_codes);
    if(code<=0) continue; // ignore NA
    int found = 0;
    for(int j=0; j<m; ++j){
      if(codes[j]==code){ counts[j]++; found=1; break; }
    }
    if(!found){ codes[m]=code; counts[m]=1; m++; }
  }
  return m;
}

C_CHANGESTAT_FN(c_cov_match_GW){
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  const double *P = INPUT_PARAM;

  const int n1         = (int)P[0];
  const int K          = (int)P[1];
  const int norm_mode  = (int)P[2]; // 0 none, 1 by_group, 2 global
  const int has_kappa  = (int)P[3]; // 0/1
  const int kappa_code = (int)P[4];

  const double *lambdas = P + 5;
  const double *z_codes = P + 5 + K;

  // Repérage acteur/groupe
  Vertex t = tail, h = head;
  Vertex v_actor = (t <= n1) ? t : h;
  Vertex v_group = (t >  n1) ? t : h;
  if(v_actor<=0 || v_actor>n1) return;
  if(v_group<=n1) return;

  const int is_add = edgestate ? 0 : 1; // 1=ajout, 0=retrait

  // Infos locales du groupe
  int maxbuf = n1 < 8192 ? n1 : 8192;
  Vertex actors_buf[8192];
  int    codes_buf[8192];
  int    counts_buf[8192];

  int na = neighbors_actors_of_group(nwp, v_group, actors_buf, n1);
  if(na>maxbuf) na = maxbuf;
  int m = histogram_codes(actors_buf, na, z_codes, codes_buf, counts_buf);

  const int r_star = code_of_actor(v_actor, z_codes);
  int n_gr_old = 0;
  for(int j=0;j<m;++j) if(codes_buf[j]==r_star){ n_gr_old = counts_buf[j]; break; }
  const int n_g_old = na;

  // Pré-calcul n_{g,κ} si ciblé
  int n_gk_old = 0;
  if(has_kappa && kappa_code>0){
    for(int u=0; u<m; ++u) if(codes_buf[u]==kappa_code){ n_gk_old = counts_buf[u]; break; }
  }

  // Boucle sur λ
  for(int j=0; j<K; ++j){
    const double lambda = lambdas[j];
    const double rlam   = (lambda - 1.0) / lambda; // r_λ in [0,1)
    double delta_non_norm = 0.0;

    // Δ brute (non normalisée)
    // ajout:  r^m ; retrait: -r^{m-1}. Si r*=κ requis, sinon 0.
    if(has_kappa){
      if(r_star==kappa_code){
        if(is_add){
          delta_non_norm = pow(rlam, (double)n_gk_old);
        }else{
          // n_gk_old>=1 si l'arête existe
          delta_non_norm = -pow(rlam, (double)(n_gk_old-1));
        }
      }else{
        delta_non_norm = 0.0;
      }
    }else{
      if(is_add){
        delta_non_norm = pow(rlam, (double)n_gr_old);
      }else{
        delta_non_norm = -pow(rlam, (double)(n_gr_old-1));
      }
    }

    double delta = delta_non_norm;

    if(norm_mode==1){
      // by_group : Δ = (N_plus/D_plus) - (N_minus/D_minus)
      // Num(g) = Σ_r λ(1-r^{n_{g,r}})   (ou modalité κ seule)
      // Den(g) = λ(1-r^{n_g})
      double N_minus = 0.0;
      if(has_kappa){
        N_minus = lambda * (1.0 - pow(rlam, (double)n_gk_old));
      }else{
        for(int u=0; u<m; ++u){
          N_minus += lambda * (1.0 - pow(rlam, (double)counts_buf[u]));
        }
      }
      double D_minus = lambda * (1.0 - pow(rlam, (double)n_g_old));

      // Mises à jour locales
      const int n_g_new  = n_g_old + (is_add ? +1 : -1);
      double N_plus = N_minus;
      if(has_kappa){
        // Seule la cellule κ change : +λ r^{m} (ajout) ou -λ r^{m-1} (retrait)
        if(r_star==kappa_code){
          if(is_add)  N_plus += lambda * ( pow(rlam,(double)n_gk_old) );
          else        N_plus -= lambda * ( pow(rlam,(double)(n_gk_old-1)) );
        }
      }else{
        // Cellule r* change
        if(is_add)  N_plus += lambda * ( pow(rlam,(double)n_gr_old) );
        else        N_plus -= lambda * ( pow(rlam,(double)(n_gr_old-1)) );
      }
      double D_plus = lambda * (1.0 - pow(rlam, (double)n_g_new));

      double ratio_minus = (D_minus>0.0) ? (N_minus/D_minus) : 0.0;
      double ratio_plus  = (D_plus >0.0) ? (N_plus /D_plus ) : 0.0;

      delta = ratio_plus - ratio_minus;

    }else if(norm_mode==2){
      // global : diviser par λ(1-r^{N1}) (constante)
      const double Dglob = lambda * (1.0 - pow(rlam, (double)n1));
      delta = (Dglob>0.0) ? (delta_non_norm / Dglob) : 0.0;
    }

    CHANGE_STAT[j] += delta;
  }
}
