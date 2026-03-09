/**
 * @file changestat_cov_match.c
 * @brief Change statistic for the ERPM term `cov_match` (multi-toggle form).
 *
 * @author Jérémie Chichignoud
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
 *  `cov_match`, which measures within-group categorical agreement among actors
 *  in a bipartite partition network. The statistic can operate either in a
 *  non-targeted mode (all categories) or in a targeted mode focusing on a
 *  specific category κ, and supports several normalization schemes.
 *
 *  ------------------------------------------------------------
 *  Implementation requirements (multi-toggle / D_CHANGESTAT_FN)
 *  ------------------------------------------------------------
 *
 *  This change statistic is designed to support multi-toggle proposals
 *  (e.g. swap, split, merge) that are internally represented in \pkg{ergm}
 *  as a list of edge toggles.
 *
 *  Consequently:
 *
 *    - the statistic is implemented using the D_CHANGESTAT_FN API,
 *      which processes multiple toggles in sequence,
 *    - the R initializer `InitErgmTerm.cov_match` MUST set
 *      `d_func = TRUE`.
 *
 *  If `d_func = FALSE`, \pkg{ergm} will incorrectly call the
 *  one-toggle C entry point (`c_` interface), which may lead
 *  to undefined behaviour or segmentation faults.
 *
 *  ------------------------------------------------------------
 *  Compiled symbol naming convention
 *  ------------------------------------------------------------
 *
 *  The C entry point must be declared as:
 *
 *      D_CHANGESTAT_FN(d_cov_match)
 *
 *  which exposes the compiled symbol `d_cov_match`.
 *
 *  A symbol named `c_cov_match` must NOT be exposed with a
 *  D-signature, because \pkg{ergm} may resolve it as the
 *  one-toggle entry point.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  The bipartite network represents a partition structure:
 *
 *    - actor mode  : vertices corresponding to actors, each carrying
 *                    a categorical covariate c(i),
 *    - group mode  : vertices corresponding to structural groups.
 *
 *  Membership is encoded by edges between actors and groups.
 *
 *  For each group g:
 *
 *    - A(g)      = set of actors adjacent to g,
 *    - n_g       = |A(g)| = group size,
 *    - c(i)      = category code of actor i,
 *    - n_{g,r}   = number of actors in group g with category r.
 *
 *  For an integer k ≥ 1, define the combinatorial quantity:
 *
 *      S_k(B; c) = ∑_g ∑_r C(n_{g,r}, k),
 *
 *  where C(a,b) denotes the binomial coefficient (with
 *  C(a,b) = 0 whenever a < b).
 *
 *  Targeted version for a specific category κ:
 *
 *      S_k^{(κ)}(B; c) = ∑_g C(n_{g,κ}, k).
 *
 *  ------------------------------------------------------------
 *  Normalisation modes
 *  ------------------------------------------------------------
 *
 *  The statistic can be normalized in three different ways:
 *
 *  1) "none" (norm_mode = 0)
 *
 *        statistic = S_k(B; c)
 *
 *     or in targeted mode:
 *
 *        statistic = S_k^{(κ)}(B; c)
 *
 *
 *  2) "by_group" (norm_mode = 1)
 *
 *     Non-targeted:
 *
 *        ∑_g [ ( ∑_r C(n_{g,r}, k) ) / C(n_g, k) ]
 *
 *     Targeted:
 *
 *        ∑_g [ C(n_{g,κ}, k) / C(n_g, k) ]
 *
 *
 *  3) "global" (norm_mode = 2)
 *
 *     Non-targeted:
 *
 *        ∑_g [ ( ∑_r C(n_{g,r}, k) ) / n_g ]
 *
 *     Targeted:
 *
 *        ∑_g [ C(n_{g,κ}, k) / n_g ]
 *
 *     Groups with n_g = 0 contribute zero.
 *
 *
 *  Special case (k = 1, targeted, by_group):
 *
 *        statistic = ∑_g 1[n_{g,κ} ≥ 1]
 *
 *
 *  ------------------------------------------------------------
 *  Multi-toggle evaluation in \pkg{ergm}
 *  ------------------------------------------------------------
 *
 *  In multi-toggle proposals, several toggles may affect:
 *
 *    - different groups, or
 *    - the same group multiple times within a single proposal.
 *
 *  Each toggle must therefore be evaluated with respect to the
 *  CURRENT intermediate network state after applying all previous
 *  toggles in the proposal.
 *
 *  The standard \pkg{ergm} multi-toggle pattern is used:
 *
 *    for each toggle i:
 *
 *      1) read endpoints (TAIL(i), HEAD(i)),
 *      2) compute the edge state BEFORE the toggle,
 *      3) recompute the affected group contribution,
 *      4) accumulate the local change Δ,
 *      5) temporarily apply the toggle (TOGGLE_IF_MORE_TO_COME).
 *
 *  After processing all toggles, temporary modifications are reverted
 *  using UNDO_PREVIOUS_TOGGLES.
 *
 *
 *  ------------------------------------------------------------
 *  INPUT_PARAM layout (from InitErgmTerm.cov_match)
 *  ------------------------------------------------------------
 *
 *  The R initializer packs parameters as:
 *
 *      INPUT_PARAM = c(
 *        n1,
 *        K,
 *        norm_mode,
 *        has_kappa,
 *        kappa_code,
 *        ks[1:K],
 *        z_codes[1:n1]
 *      )
 *
 *  where:
 *
 *    - n1          = number of actors (size of actor mode),
 *    - K           = number of k values requested,
 *    - norm_mode   = normalization mode,
 *    - has_kappa   = indicator for targeted mode,
 *    - kappa_code  = targeted category code,
 *    - ks          = vector of k values,
 *    - z_codes     = categorical covariate codes for actors.
 *
 *
 *  At the C level (double* P = INPUT_PARAM):
 *
 *      P[0]          = n1
 *      P[1]          = K
 *      P[2]          = norm_mode
 *      P[3]          = has_kappa
 *      P[4]          = kappa_code
 *      P[5 .. 5+K-1] = ks[0 .. K-1]
 *      P[5+K .. ]    = z_codes[0 .. n1-1]
 *
 *  The term returns K statistics:
 *
 *      - N_CHANGE_STATS = K
 *      - one statistic component per value of k.
 *
 *
 *  ------------------------------------------------------------
 *  Computational complexity
 *  ------------------------------------------------------------
 *
 *  For each toggle:
 *
 *    - collecting the actors adjacent to a group is O(deg(g)),
 *    - constructing category counts is approximately
 *      O(deg(g) × m), where m is the number of distinct
 *      categories present in the group.
 *
 *  No caching is performed across toggles. This keeps the
 *  implementation simple and robust, which is acceptable
 *  for typical ERPM group sizes.
 */

#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <math.h>
#include <string.h>

/* -------------------------------------------------------------------------- */
/* Debug switches                                                             */
/* -------------------------------------------------------------------------- */

/**
 * @def DEBUG_COV_MATCH
 * @brief Enable verbose debugging output for ::d_cov_match.
 *
 * Set to 1 to print diagnostic info to the R console:
 *   - ntoggles (multi-toggle detection),
 *   - actor/group endpoints per toggle,
 *   - group size and key category counts,
 *   - per-k deltas.
 */
#define DEBUG_COV_MATCH 0
#define UNUSED_WARNING(x) (void)(x)

/* -------------------------------------------------------------------------- */
/* Helper: category code lookup                                               */
/* -------------------------------------------------------------------------- */

/**
 * @brief Retrieve the integer category code for a given actor.
 *
 * z_codes is stored as a double array of length n1:
 *   z_codes[i-1] = category code for actor vertex i (1..n1), 0 if NA/undefined.
 */
static inline int code_of_actor(Vertex i, const double *z_codes){
  return (int)z_codes[(size_t)(i-1)];
}

/* -------------------------------------------------------------------------- */
/* Helper: collect unique actor neighbours of a group                         */
/* -------------------------------------------------------------------------- */

/**
 * @brief Collect unique actor neighbours of a group vertex g (group mode).
 *
 * @param nwp     ergm network workspace.
 * @param g       group vertex index (> n1).
 * @param actors  output buffer of length n1 (Vertex).
 * @param n1      actor-mode size.
 *
 * @return number of unique actor neighbours written into actors[].
 *
 * Notes:
 * - We deduplicate using a temporary bitmap 'seen' of length n1.
 * - This function observes the CURRENT network state (including any temporary
 *   toggles already applied in a multi-toggle proposal).
 */
static int neighbors_actors_of_group(Network *nwp, Vertex g, Vertex *actors, int n1){
  int cnt = 0;
  unsigned char *seen = (unsigned char*)R_Calloc((size_t)n1, unsigned char);

  Vertex h; Edge e;

  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        actors[cnt++] = h;
      }
    }
  }

  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        actors[cnt++] = h;
      }
    }
  }

  R_Free(seen);
  return cnt;
}

/* -------------------------------------------------------------------------- */
/* Helper: histogram of category codes in a group                             */
/* -------------------------------------------------------------------------- */

/**
 * @brief Build a histogram (codes[], counts[]) from actor list actors[0..na-1].
 *
 * Simple O(na*m) implementation, sufficient for small/medium groups.
 */
static int histogram_codes(const Vertex *actors, int na,
                           const double *z_codes,
                           int *codes, int *counts){
  int m = 0;
  for(int a=0; a<na; ++a){
    int code = code_of_actor(actors[a], z_codes);
    if(code <= 0) continue; /* ignore NA/undefined */
    int found = 0;
    for(int j=0; j<m; ++j){
      if(codes[j] == code){
        counts[j]++; found = 1; break;
      }
    }
    if(!found){
      codes[m]  = code;
      counts[m] = 1;
      m++;
    }
  }
  return m;
}

/* -------------------------------------------------------------------------- */
/* Helper: compute group-level unnormalised quantity N(g) for a given k       */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute N(g) for the current group histogram:
 *   - targeted: N = C(n_{g,κ}, k)
 *   - non-targeted: N = sum_r C(n_{g,r}, k)
 */
static double group_N_for_k(int k,
                            int has_kappa,
                            int kappa_code,
                            const int *codes,
                            const int *counts,
                            int m){
  if(has_kappa){
    int n_gk = 0;
    if(kappa_code > 0){
      for(int u=0; u<m; ++u){
        if(codes[u] == kappa_code){ n_gk = counts[u]; break; }
      }
    }
    return CHOOSE(n_gk, k);
  }else{
    double s = 0.0;
    for(int u=0; u<m; ++u) s += CHOOSE(counts[u], k);
    return s;
  }
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_match (MULTI-TOGGLE)                                 */
/* -------------------------------------------------------------------------- */

/**
 * @brief Multi-toggle change statistic for `cov_match`.
 *
 * This is the D_ entrypoint called by ergm when `d_func=TRUE` is returned by
 * the R initializer (InitErgmTerm.cov_match).
 *
 * For each toggle i, we:
 *   - identify (actor, group),
 *   - recompute group membership + category histogram under current intermediate state,
 *   - compute per-k Δ according to normalization and targeted mode,
 *   - apply TOGGLE_IF_MORE_TO_COME(i) so subsequent toggles see updated state.
 */
D_CHANGESTAT_FN(d_cov_match){

#if DEBUG_COV_MATCH
  static int seen_multitoggle = 0;
  if(ntoggles > 1 && seen_multitoggle < 20){
    Rprintf("[cov_match] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen_multitoggle++;
  }
#endif

  /* Reset output buffer (this proposal). */
  ZERO_ALL_CHANGESTATS();

  /* Read packed inputs. */
  const double *P = INPUT_PARAM;

  const int n1         = (int)P[0];
  const int K          = (int)P[1];
  const int norm_mode  = (int)P[2]; /* 0 none, 1 by_group, 2 global */
  const int has_kappa  = (int)P[3]; /* 0/1 */
  const int kappa_code = (int)P[4];

  const double *ks_d    = P + 5;
  const double *z_codes = P + 5 + K;

  /* Scratch buffers sized by n1 (safe, deterministic).
   * - actors_buf length n1
   * - histogram buffers length n1 (worst-case distinct categories == group size)
   */
  Vertex *actors_buf = (Vertex*)R_Calloc((size_t)n1, Vertex);
  int    *codes_buf  = (int*)R_Calloc((size_t)n1, int);
  int    *counts_buf = (int*)R_Calloc((size_t)n1, int);

  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex t = TAIL(i);
    Vertex h = HEAD(i);

    /* Determine current edge state BEFORE toggling. */
    int edgestate = DIRECTED ? IS_OUTEDGE(t, h) : IS_UNDIRECTED_EDGE(t, h);
    int is_add    = edgestate ? 0 : 1;

    /* Identify actor/group endpoints (one in each mode). */
    Vertex v_actor = (t <= (Vertex)n1) ? t : h;
    Vertex v_group = (t >  (Vertex)n1) ? t : h;

    /* Guard: if malformed, skip safely. */
    if(v_actor <= 0 || v_actor > (Vertex)n1 || v_group <= (Vertex)n1){
      /* Still must apply toggle for consistency? No: this should never happen
       * in a correctly constrained bipartite proposal. We skip and do not toggle.
       */
#if DEBUG_COV_MATCH
      Rprintf("[cov_match][WARN] malformed toggle i=%d tail=%d head=%d (n1=%d)\n",
              i, (int)t, (int)h, n1);
#endif
      continue;
    }

    /* Recompute group membership in current intermediate state. */
    int na = neighbors_actors_of_group(nwp, v_group, actors_buf, n1);

    /* Build histogram of category codes for the group in current state. */
    int m = histogram_codes(actors_buf, na, z_codes, codes_buf, counts_buf);

    /* Category of the toggled actor. */
    const int r_star = code_of_actor(v_actor, z_codes);

    /* Count in actor's category r* within the group (old, before this toggle). */
    int n_gr_old = 0;
    if(r_star > 0){
      for(int u=0; u<m; ++u){
        if(codes_buf[u] == r_star){ n_gr_old = counts_buf[u]; break; }
      }
    }

    /* Current group size before this toggle. */
    const int n_g_old = na;

#if DEBUG_COV_MATCH
    Rprintf("[cov_match:toggle] i=%d tail=%d head=%d | actor=%d group=%d | edgestate=%d is_add=%d | n_g_old=%d r*=%d n_gr_old=%d\n",
            i, (int)t, (int)h, (int)v_actor, (int)v_group, edgestate, is_add, n_g_old, r_star, n_gr_old);
#endif

    /* Update each vectorized k. */
    for(int j=0; j<K; ++j){
      const int k = (int)ks_d[j];

      /* Special case:
       *   k=1, by_group, targeted:
       *     contrib(g) = 1[ n_{g,κ} >= 1 ]
       * so Δ = ind_new - ind_old.
       */
      if(norm_mode == 1 && has_kappa && k == 1){
        int n_gk_old = 0;
        if(kappa_code > 0){
          for(int u=0; u<m; ++u){
            if(codes_buf[u] == kappa_code){ n_gk_old = counts_buf[u]; break; }
          }
        }

        int n_gk_new = n_gk_old;
        if(r_star == kappa_code){
          n_gk_new += (is_add ? +1 : -1);
          if(n_gk_new < 0) n_gk_new = 0;
        }

        double delta = (n_gk_new > 0) - (n_gk_old > 0);
        CHANGE_STAT[j] += delta;

#if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=1][by_group][kappa] n_gk_old=%d n_gk_new=%d delta=%g\n",
                n_gk_old, n_gk_new, delta);
#endif
        continue;
      }

      /* Unnormalised local delta for THIS toggle (same identities as one-toggle). */
      double delta_non_norm = 0.0;

      if(has_kappa){
        int n_gk_old = 0;
        if(kappa_code > 0){
          for(int u=0; u<m; ++u){
            if(codes_buf[u] == kappa_code){ n_gk_old = counts_buf[u]; break; }
          }
        }

        if(is_add){
          delta_non_norm = (r_star == kappa_code) ? CHOOSE(n_gk_old, k-1) : 0.0;
        }else{
          delta_non_norm = (r_star == kappa_code) ? -CHOOSE((n_gk_old - 1), k-1) : 0.0;
        }
      }else{
        if(is_add){
          delta_non_norm = CHOOSE(n_gr_old, k-1);
        }else{
          delta_non_norm = -CHOOSE((n_gr_old - 1), k-1);
        }
      }

      /* Apply requested normalisation. */
      double delta = delta_non_norm;

      if(norm_mode == 1){
        /* by_group: ratio before vs after (group-level). */
        double N_minus = group_N_for_k(k, has_kappa, kappa_code, codes_buf, counts_buf, m);
        double D_minus = CHOOSE(n_g_old, k);

        /* After-toggle group size. */
        int n_g_new = n_g_old + (is_add ? +1 : -1);
        if(n_g_new < 0) n_g_new = 0;

        /* N_plus can be updated by delta_non_norm. */
        double N_plus = N_minus + delta_non_norm;
        double D_plus = CHOOSE(n_g_new, k);

        double ratio_minus = (D_minus > 0.0) ? (N_minus / D_minus) : 0.0;
        double ratio_plus  = (D_plus  > 0.0) ? (N_plus  / D_plus ) : 0.0;

        delta = ratio_plus - ratio_minus;
      }
      else if(norm_mode == 2){
        /* global: contrib(g) = N(g)/n_g (or 0 if n_g==0). */
        double N_minus = group_N_for_k(k, has_kappa, kappa_code, codes_buf, counts_buf, m);

        int n_g_new = n_g_old + (is_add ? +1 : -1);
        if(n_g_new < 0) n_g_new = 0;

        double N_plus = N_minus + delta_non_norm;

        double contrib_minus = (n_g_old > 0) ? (N_minus / (double)n_g_old) : 0.0;
        double contrib_plus  = (n_g_new > 0) ? (N_plus  / (double)n_g_new) : 0.0;

        delta = contrib_plus - contrib_minus;
      }

      CHANGE_STAT[j] += delta;

#if DEBUG_COV_MATCH
      Rprintf("[cov_match:k=%d] delta_non_norm=%.6f delta=%.6f cumul=%.6f\n",
              k, delta_non_norm, delta, CHANGE_STAT[j]);
#endif
    }

    /* Temporarily apply toggle so subsequent toggles see updated network state. */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* Undo all temporary toggles (restore original network state). */
  UNDO_PREVIOUS_TOGGLES(i);

  R_Free(counts_buf);
  R_Free(codes_buf);
  R_Free(actors_buf);
}
