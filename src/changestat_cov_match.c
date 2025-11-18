// =============================================================================
/**
 * @file changestat_cov_match.c
 * @brief Change statistic (forme "un-toggle") pour le terme `cov_match`
 *
 * @details
 * Calcule, lors d’un toggle (ajout/retrait d’une arête acteur–groupe),
 * la variation locale de la statistique :
 *
 *   - Version brute (non normalisée) :
 *       S_k(B;c) = ∑_g ∑_r C(n_{g,r}, k)
 *
 *   - Version ciblée (category = κ) :
 *       S_k^{(κ)}(B;c) = ∑_g C(n_{g,κ}, k)
 *
 *   - Normalisation "by_group" :
 *       ∑_g [  (∑_r C(n_{g,r}, k)) / C(n_g, k)  ]      (ou C(n_{g,κ},k)/C(n_g,k) en ciblé)
 *
 *   - Normalisation "global" :
 *       (1 / C(N_actors, k)) * S_k(·)                   (ou S_k^{(κ)}(·) en ciblé)
 *
 * où n_{g,r} = # acteurs de modalité r dans le groupe g, n_g = taille du groupe g,
 * et C(a,b)=CHOOSE(a,b) avec la convention C(a,b)=0 si a<b.
 *
 * ## Contexte {ergm} (forme un-toggle)
 * Cf. l’exemple `squared_sizes` fourni. Même cycle de vie, même contrat.
 *
 * ## Entrée (depuis R → C, via InitErgmTerm.cov_match)
 * - INPUT_PARAM (double*) contient, dans l’ordre :
 *   [0]   = n1                        (taille du mode 1 = #acteurs)
 *   [1]   = K                         (nombre de valeurs de k vectorisées)
 *   [2]   = norm_mode                 (0 = none, 1 = by_group, 2 = global)
 *   [3]   = has_kappa                 (0 = non ciblé, 1 = ciblé)
 *   [4]   = kappa_code                (code entier de la modalité ciblée si has_kappa=1, sinon 0)
 *   [5..(5+K-1)]          = ks[j]     (valeurs k ≥ 1)
 *   [5+K .. 5+K+(n1-1)]   = z[1..n1]  (codes entiers des modalités par acteur, 0 si NA)
 *
 * - COEF NAMES : longueur K ⇒ N_CHANGE_STATS = K
 *
 * ## Sortie (par appel)
 * - CHANGE_STAT[j] reçoit Δ_j pour le toggle courant, selon k = ks[j] et la normalisation.
 *
 * ## Hypothèses
 * - Réseau biparti : acteurs 1..n1, groupes (n1+1)..N.
 * - Le toggle connecte exactement un acteur (mode 1) et un groupe (mode 2).
 *
 * ## Complexité
 * - O(deg(g)) pour construire les effectifs de modalités du groupe g au moment du toggle.
 *   (Optimisable en stockant des états par groupe ; suffisant dans un premier jet.)
 *
 * ## Remarques
 * - Identités utiles :
 *   * Ajout d’un acteur de modalité r* dans g :
 *       C(n_{g,r*}+1, k) - C(n_{g,r*}, k) = C(n_{g,r*}, k-1)
 *     donc Δ_non_norm = + C(n_{g,r*}^{old}, k-1)
 *   * Retrait :
 *       C(n_{g,r*}-1, k) - C(n_{g,r*}, k) = - C(n_{g,r*}-1, k-1)
 *     donc Δ_non_norm = - C(n_{g,r*}^{new}-1, k-1)  (avec new = old - 1)
 *
 */
// =============================================================================

#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <math.h>
#include <string.h>

#define DEBUG_COV_MATCH 0   /* mettre à 0 pour couper le debug */
#define UNUSED_WARNING(x) (void)x

// -----------------------------------------------------------------------------

static inline int code_of_actor(Vertex i, const double *z_codes){
  // z_codes est 1-indexé comme les sommets (acteurs 1..n1)
  return (int)z_codes[(size_t)(i-1)];
}

/* Retourne les acteurs voisins uniques du groupe g dans actors[]
 * Déduplique via seen[] pour éviter les doubles comptes OUT/IN. */
static int neighbors_actors_of_group(Network *nwp, Vertex g, Vertex *actors, int n1){
  int cnt = 0;
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char); /* Calloc -> 0 */

  Vertex h; Edge e;

  /* OUT-neighbors g -> h */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        actors[cnt++] = h;
      }
    }
  }

  /* IN-neighbors h -> g */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        actors[cnt++] = h;
      }
    }
  }

  Free(seen);
  return cnt;
}

// Construit un histogramme (codes[], counts[]) des fréquences par modalité dans le groupe.
static int histogram_codes(const Vertex *actors, int na, const double *z_codes, int *codes, int *counts){
  int m = 0; // nombre de modalités distinctes vues
  for(int a=0; a<na; ++a){
    int code = code_of_actor(actors[a], z_codes);
    if(code<=0) continue; // ignorer NA/indéfinis
    int found = 0;
    for(int j=0; j<m; ++j){
      if(codes[j]==code){ counts[j]++; found=1; break; }
    }
    if(!found){
      codes[m]=code; counts[m]=1; m++;
    }
  }
  return m;
}

C_CHANGESTAT_FN(c_cov_match){
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  // ----- Lecture des inputs packés par R -----
  const double *P = INPUT_PARAM;

  const int n1          = (int)P[0];
  const int K           = (int)P[1];
  const int norm_mode   = (int)P[2];  // 0 none, 1 by_group, 2 global
  const int has_kappa   = (int)P[3];  // 0/1
  const int kappa_code  = (int)P[4];

  const double *ks_d    = P + 5;
  const double *z_codes = P + 5 + K;  // longueur n1

  #if DEBUG_COV_MATCH
    Rprintf("[cov_match:init] n1=%d K=%d norm=%d has_kappa=%d kappa_code=%d ks{",
            n1, K, norm_mode, has_kappa, kappa_code);
    for(int jj=0;jj<K;++jj){ Rprintf("%s%d", (jj?",":""), (int)ks_d[jj]); }
    Rprintf("}\n");
  #endif

  // ----- Repérage acteur / groupe dans le toggle -----
  Vertex t = tail, h = head;
  Vertex v_actor = (t <= n1) ? t : h;
  Vertex v_group = (t >  n1) ? t : h;

  // Sanity minimal
  if(v_actor<=0 || v_actor>n1) return;
  if(v_group<=n1) return;

  // ----- Caractère du toggle -----
  // edgestate=1 ⇒ l’arête existe ⇒ toggle=retrait ;   edgestate=0 ⇒ toggle=ajout.
  const int is_add = edgestate ? 0 : 1;

  #if DEBUG_COV_MATCH
    Rprintf("[cov_match:toggle] actor=%d group=%d is_add=%d OUT_DEG[g]=%d IN_DEG[g]=%d\n",
            (int)v_actor, (int)v_group, is_add, (int)OUT_DEG[v_group], (int)IN_DEG[v_group]);
  #endif

  // ----- Préparer les infos locales sur le groupe affecté -----
  Vertex actors_buf[8192];
  int    codes_buf[ 8192];
  int    counts_buf[8192];

  int maxbuf = (n1 < 8192) ? n1 : 8192;
  int na = neighbors_actors_of_group(nwp, v_group, actors_buf, n1);
  if(na>maxbuf) na = maxbuf; // sécurité soft

  int m = histogram_codes(actors_buf, na, z_codes, codes_buf, counts_buf);

  // ----- Récupérer la modalité de l’acteur du toggle -----
  const int r_star = code_of_actor(v_actor, z_codes);
  int n_gr_old = 0;
  for(int j=0;j<m;++j) if(codes_buf[j]==r_star) { n_gr_old = counts_buf[j]; break; }

  const int n_g_old = na;

  #if DEBUG_COV_MATCH
    Rprintf("[cov_match:group] n_g_old=%d r*=%d n_gr_old=%d m=%d\n",
            n_g_old, r_star, n_gr_old, m);
    if(has_kappa){
      int n_gk_old_dbg = 0;
      if(kappa_code>0){
        for(int u=0; u<m; ++u) if(codes_buf[u]==kappa_code){ n_gk_old_dbg = counts_buf[u]; break; }
      }
      Rprintf("[cov_match:group] kappa=%d n_gk_old=%d\n", kappa_code, n_gk_old_dbg);
    }
  #endif

  // ----- Boucle sur les k vectorisés -----
  for(int j=0; j<K; ++j){
    const int k = (int)ks_d[j];

    /* Cas spécial: k==1, by_group, ciblé(κ) → somme_g 1_{n_{g,κ}≥1} */
    if(norm_mode==1 && has_kappa && k==1){
      int n_gk_old = 0;
      if(kappa_code>0){
        for(int u=0; u<m; ++u) if(codes_buf[u]==kappa_code){ n_gk_old = counts_buf[u]; break; }
      }
      int n_gk_new = n_gk_old;
      if(r_star==kappa_code){
        n_gk_new += (is_add ? +1 : -1);
        if(n_gk_new < 0) n_gk_new = 0;
      }
      double delta_ind = (n_gk_new>0) - (n_gk_old>0);
      CHANGE_STAT[j] += delta_ind;

      #if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=1][by_group][kappa] n_g_old=%d n_gk_old=%d n_gk_new=%d delta=%g\n",
                n_g_old, n_gk_old, n_gk_new, delta_ind);
      #endif
      continue;
    }

    double delta_non_norm = 0.0;

    if(has_kappa){
      int n_gk_old = 0;
      if(kappa_code>0){
        for(int u=0; u<m; ++u) if(codes_buf[u]==kappa_code){ n_gk_old = counts_buf[u]; break; }
      }
      if(is_add){
        delta_non_norm = (r_star==kappa_code) ? CHOOSE(n_gk_old, k-1) : 0.0;
      }else{
        delta_non_norm = (r_star==kappa_code) ? -CHOOSE((n_gk_old-1), k-1) : 0.0;
      }
      #if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=%d][kappa=%d] n_gk_old=%d delta_non_norm=%.6f\n",
                k, kappa_code, n_gk_old, delta_non_norm);
      #endif
    }else{
      if(is_add){
        delta_non_norm = CHOOSE(n_gr_old, k-1);
      }else{
        delta_non_norm = -CHOOSE((n_gr_old-1), k-1);
      }
      #if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=%d][non-ciblé] n_gr_old=%d delta_non_norm=%.6f\n",
                k, n_gr_old, delta_non_norm);
      #endif
    }

    // ----- Appliquer la normalisation -----
    double delta = delta_non_norm;

    if(norm_mode==1){
      // by_group : Δ = (N_plus/D_plus) - (N_minus/D_minus)
      double N_minus = 0.0;
      if(has_kappa){
        int n_gk_old = 0;
        if(kappa_code>0){
          for(int u=0; u<m; ++u) if(codes_buf[u]==kappa_code){ n_gk_old = counts_buf[u]; break; }
        }
        N_minus = CHOOSE(n_gk_old, k);
      }else{
        for(int u=0; u<m; ++u) N_minus += CHOOSE(counts_buf[u], k);
      }
      double D_minus = CHOOSE(n_g_old, k);

      double N_plus = N_minus + delta_non_norm;
      double D_plus = CHOOSE(n_g_old + (is_add?+1:-1), k);

      double ratio_minus = (D_minus>0.0) ? (N_minus/D_minus) : 0.0;
      double ratio_plus  = (D_plus >0.0) ? (N_plus /D_plus ) : 0.0;

      delta = ratio_plus - ratio_minus;

      #if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=%d][by_group] N-=%g D-=%g | N+=%g D+=%g | delta=%.6f\n",
                k, N_minus, D_minus, N_plus, D_plus, delta);
      #endif
    }
    else if(norm_mode==2){
      double Cglob = CHOOSE(n1, k);
      delta = (Cglob>0.0) ? (delta_non_norm / Cglob) : 0.0;
      #if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=%d][global] Cglob=%g delta_non_norm=%.6f -> delta=%.6f\n",
                k, Cglob, delta_non_norm, delta);
      #endif
    }

    CHANGE_STAT[j] += delta;
  }
}
