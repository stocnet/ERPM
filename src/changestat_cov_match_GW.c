/**
 * @file changestat_cov_match_GW.c
 * @brief  Change statistic for the ERPM term `cov_match_GW` (multi-toggle form).
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
 *  `cov_match_GW`, which applies a geometrically weighted transform to
 *  group-level category counts for a categorical actor covariate.
 *
 *  ------------------------------------------------------------
 *  IMPORTANT (multi-toggle / D_CHANGESTAT_FN)
 *  ------------------------------------------------------------
 *  This changestat is implemented in the D_ API (multi-toggle):
 *
 *      D_CHANGESTAT_FN(d_cov_match_GW)
 *
 *  because ERPM/partition moves (swap/split/merge) can be represented as a
 *  list of membership toggles. Several toggles may touch the same group in a
 *  single proposal, so we must evaluate them sequentially on the intermediate
 *  state:
 *
 *    - For each toggle i, compute the local Δ using the CURRENT intermediate
 *      network state (i.e. after toggles 0..i-1 have been temporarily applied).
 *    - Then apply the toggle temporarily if more toggles remain:
 *        TOGGLE_IF_MORE_TO_COME(i)
 *    - At the end, undo all temporary toggles:
 *        UNDO_PREVIOUS_TOGGLES(i)
 *
 *  If you implement this as a one-toggle C_CHANGESTAT_FN and your proposal
 *  generates multi-toggle moves, the statistic will be wrong (or ergm may call
 *  the function with the wrong signature if d_func=TRUE is missing on the R side).
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  The network is bipartite:
 *    - actor mode  = vertices representing actors, each carrying a categorical
 *                    covariate c(i),
 *    - group mode  = vertices representing structural groups.
 *
 *  For each group vertex g in the group mode:
 *    - A(g)      = set of actors connected to g,
 *    - n_g       = |A(g)| = size of group g,
 *    - c(i)      = category code of actor i,
 *    - n_{g,r}   = number of actors in g with category r.
 *
 *  For a given decay parameter λ ≥ 1, define:
 *
 *      r_λ = (λ - 1) / λ  ∈ [0, 1).
 *
 *  In this implementation, the ERPM effect is encoded so that a membership
 *  toggle of an actor with category r* produces a local contribution for the
 *  targeted count n_{g,r*} of the form:
 *
 *    - addition (n_{g,r*} = m → m+1):
 *        Δ_non_norm =  r_λ^m,
 *    - deletion (n_{g,r*} = m → m-1):
 *        Δ_non_norm = -r_λ^{m-1}.
 *
 *  When a target category κ is specified, the transform is restricted to the
 *  single count n_{g,κ}. Otherwise, it aggregates over all categories present
 *  in the group.
 *
 *  Normalisation modes:
 *
 *    - "none" (norm_mode = 0):
 *        statistic ∝ ∑_g ∑_r GW(n_{g,r}; λ),
 *        local Δ given directly by Δ_non_norm.
 *
 *    - "by_group" (norm_mode = 1):
 *        For each group g, build:
 *
 *          Num(g) = ∑_r λ (1 - r_λ^{n_{g,r}})      (or λ (1 - r_λ^{n_{g,κ}})
 *                                                 in targeted version),
 *          Den(g) = λ (1 - r_λ^{n_g}),
 *
 *        then the contribution of g is:
 *
 *          Num(g) / Den(g),
 *
 *        so for a toggle in group g, the local change is:
 *
 *          Δ = [Num(g)_after / Den(g)_after] - [Num(g)_before / Den(g)_before].
 *
 *    - "global" (norm_mode = 2):
 *        A global normalisation rescales the non-normalised statistic by
 *        a group-size independent constant:
 *
 *          Δ = Δ_non_norm / [λ (1 - r_λ^{N_actors})],
 *
 *        where N_actors = n1 is the number of actors.
 *
 *  Only the group touched by the membership toggle is recomputed. All other
 *  groups are unaffected.
 *
 *  ------------------------------------------------------------
 *  Bipartite structure and actors/groups
 *  ------------------------------------------------------------
 *
 *  At the C level, the bipartite structure is encoded as:
 *    - the first n1 vertices are actors (actor mode),
 *    - the remaining vertices are groups (group mode).
 *
 *  Each membership toggle always connects:
 *    - exactly one actor vertex (index in 1..n1),
 *    - exactly one group vertex (index > n1).
 *
 *  ------------------------------------------------------------
 *  INPUT_PARAM layout
 *  ------------------------------------------------------------
 *
 *  The R initialiser packs the parameters as:
 *
 *      INPUT_PARAM = c(
 *        n1,
 *        K,
 *        norm_mode,
 *        has_kappa,
 *        kappa_code,
 *        lambdas[1:K],
 *        z_codes[1:n1]
 *      )
 *
 *  At C level:
 *
 *    P[0]          = n1
 *    P[1]          = K
 *    P[2]          = norm_mode
 *    P[3]          = has_kappa
 *    P[4]          = kappa_code
 *    P[5 .. 5+K-1] = lambdas[0 .. K-1]
 *    P[5+K .. ]    = z_codes[0 .. n1-1]
 *
 *  N_CHANGE_STATS is equal to K, one statistic per λ_j.
 */

#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <math.h>
#include <string.h>

/**
 * @def DEBUG_COV_MATCH_GW
 * @brief Enable verbose debugging output for ::d_cov_match_GW.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - actor and group vertices for each toggle,
 *  - group size and relevant category counts before the toggle,
 *  - intermediate non-normalised and normalised deltas per λ,
 *  - a one-line "MULTI-TOGGLE" banner when ntoggles > 1.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_COV_MATCH_GW 0

/* -------------------------------------------------------------------------- */
/* Helper: category code lookup                                               */
/* -------------------------------------------------------------------------- */

/**
 * @brief Retrieve the integer category code for a given actor.
 *
 * @details
 *  Actors are indexed 1..n1; z_codes is packed as a double array of length n1.
 *  Values <= 0 are interpreted as “no category” (e.g., NA) and ignored.
 */
static inline int code_of_actor(Vertex i, const double *z_codes){
  return (int)z_codes[(size_t)(i-1)];
}

/* -------------------------------------------------------------------------- */
/* Helper: neighbours of a group in the actor mode (unique)                   */
/* -------------------------------------------------------------------------- */

/**
 * @brief Collect unique actor neighbours of a group vertex g.
 *
 * @details
 *  This version is multi-toggle friendly and avoids O(n1) memset by using a
 *  stamp array:
 *    - stamp[idx] stores the last "mark" seen for actor idx.
 *    - To deduplicate within one call, we increment mark and compare.
 *
 * @param nwp    Network pointer.
 * @param g      Group vertex.
 * @param actors Output buffer (size at least n1).
 * @param n1     Actor-mode size.
 * @param stamp  int[n1] stamp array (persistent across calls).
 * @param mark   current mark value (incremented by caller per call).
 * @return number of unique actor neighbours written to actors[].
 */
static int neighbors_actors_of_group_stamped(Network *nwp, Vertex g,
                                             Vertex *actors, int n1,
                                             int *stamp, int mark){
  int cnt = 0;
  Vertex h;
  Edge e;

  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(stamp[idx] != mark){
        stamp[idx] = mark;
        actors[cnt++] = h;
      }
    }
  }

  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(stamp[idx] != mark){
        stamp[idx] = mark;
        actors[cnt++] = h;
      }
    }
  }

  return cnt;
}

/* -------------------------------------------------------------------------- */
/* Helper: histogram of category codes in a group                             */
/* -------------------------------------------------------------------------- */

/**
 * @brief Build a histogram of category codes for actors in a group.
 *
 * @details
 *  codes[] / counts[] are filled for distinct codes encountered among actors[].
 *  This is O(na * m) with m distinct categories in the group, which is fine for
 *  small/medium groups; can be optimised later if needed.
 *
 * @return number of distinct categories written (m).
 */
static int histogram_codes(const Vertex *actors, int na,
                           const double *z_codes,
                           int *codes, int *counts){
  int m = 0;
  for(int a=0; a<na; ++a){
    int code = code_of_actor(actors[a], z_codes);
    if(code <= 0) continue;
    int found = 0;
    for(int j=0; j<m; ++j){
      if(codes[j] == code){
        counts[j]++;
        found = 1;
        break;
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
/* Change statistic: cov_match_GW (multi-toggle)                              */
/* -------------------------------------------------------------------------- */

D_CHANGESTAT_FN(d_cov_match_GW){

#if DEBUG_COV_MATCH_GW
  static int seen_multi = 0;
  if(ntoggles > 1 && seen_multi < 20){
    Rprintf("[cov_match_GW] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen_multi++;
  }
#endif

  /* 1) Reset output buffer (this proposal may contain multiple toggles). */
  ZERO_ALL_CHANGESTATS();

  /* 2) Read packed parameters from INPUT_PARAM. */
  const double *P = INPUT_PARAM;

  const int n1         = (int)P[0];
  const int K          = (int)P[1];
  const int norm_mode  = (int)P[2]; /* 0 none, 1 by_group, 2 global */
  const int has_kappa  = (int)P[3]; /* 0/1 */
  const int kappa_code = (int)P[4];

  const double *lambdas = P + 5;
  const double *z_codes = P + 5 + K;

  /* 3) Allocate per-proposal working buffers (size n1). */
  Vertex *actors_buf = (Vertex*)R_Calloc((size_t)n1, Vertex);
  int    *codes_buf  = (int*)   R_Calloc((size_t)n1, int);
  int    *counts_buf = (int*)   R_Calloc((size_t)n1, int);
  int    *stamp      = (int*)   R_Calloc((size_t)n1, int);

  /* Stamp init: 0 means "never seen". */
  for(int i=0;i<n1;++i) stamp[i]=0;
  int mark = 1;

  /* 4) Process toggles sequentially on the intermediate state. */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex t = TAIL(i);
    Vertex h = HEAD(i);

    /* Determine current edge state BEFORE toggling (intermediate state). */
    int edgestate = DIRECTED ? IS_OUTEDGE(t, h) : IS_UNDIRECTED_EDGE(t, h);
    const int is_add = edgestate ? 0 : 1; /* 1=addition, 0=deletion */

    /* Identify actor and group endpoints. */
    Vertex v_actor = (t <= (Vertex)n1) ? t : h;
    Vertex v_group = (t >  (Vertex)n1) ? t : h;

    /* Minimal sanity checks: ensure we truly have an actor and a group. */
    if(v_actor <= 0 || v_actor > (Vertex)n1) goto toggle_apply;
    if(v_group <= (Vertex)n1) goto toggle_apply;

    /* Build group membership in the actor mode (unique actors). */
    int na = neighbors_actors_of_group_stamped(nwp, v_group, actors_buf, n1, stamp, mark++);
    const int n_g_old = na;

    /* Histogram of codes in the group (based on current intermediate state). */
    int m = histogram_codes(actors_buf, na, z_codes, codes_buf, counts_buf);

    /* Actor category code and its count in the group. */
    const int r_star = code_of_actor(v_actor, z_codes);

    int n_gr_old = 0;
    if(r_star > 0){
      for(int j=0; j<m; ++j){
        if(codes_buf[j] == r_star){
          n_gr_old = counts_buf[j];
          break;
        }
      }
    }

    /* Target category count (if applicable). */
    int n_gk_old = 0;
    if(has_kappa && kappa_code > 0){
      for(int j=0; j<m; ++j){
        if(codes_buf[j] == kappa_code){
          n_gk_old = counts_buf[j];
          break;
        }
      }
    }

#if DEBUG_COV_MATCH_GW
    Rprintf("[cov_match_GW:toggle] i=%d tail=%d head=%d | actor=%d group=%d | edgestate=%d is_add=%d | n_g_old=%d r*=%d n_gr_old=%d n_gk_old=%d (kappa=%d)\n",
            i, (int)t, (int)h, (int)v_actor, (int)v_group,
            edgestate, is_add, n_g_old, r_star, n_gr_old, n_gk_old, kappa_code);
#endif

    /* Main loop over lambdas (vectorised statistics). */
    for(int j=0; j<K; ++j){
      const double lambda = lambdas[j];
      const double rlam   = (lambda - 1.0) / lambda;

      /* 1) Non-normalised delta. */
      double delta_non_norm = 0.0;

      if(has_kappa){
        /* Targeted version: only the count n_{g,kappa} contributes. */
        if(r_star == kappa_code && r_star > 0){
          if(is_add){
            delta_non_norm = pow(rlam, (double)n_gk_old);
          }else{
            /* If edge exists, n_gk_old >= 1 for a valid deletion. */
            delta_non_norm = -pow(rlam, (double)(n_gk_old - 1));
          }
        }else{
          delta_non_norm = 0.0;
        }
      }else{
        /* Non-targeted: use count for actor's own category r*. */
        if(r_star > 0){
          if(is_add){
            delta_non_norm = pow(rlam, (double)n_gr_old);
          }else{
            delta_non_norm = -pow(rlam, (double)(n_gr_old - 1));
          }
        }else{
          /* Undefined category -> no contribution. */
          delta_non_norm = 0.0;
        }
      }

      /* 2) Apply normalisation. */
      double delta = delta_non_norm;

      if(norm_mode == 1){
        /* By-group normalisation:
         *   Δ = (Num_after/Den_after) - (Num_before/Den_before)
         */
        double N_minus = 0.0;

        if(has_kappa){
          N_minus = lambda * (1.0 - pow(rlam, (double)n_gk_old));
        }else{
          for(int u=0; u<m; ++u){
            N_minus += lambda * (1.0 - pow(rlam, (double)counts_buf[u]));
          }
        }

        const double D_minus = lambda * (1.0 - pow(rlam, (double)n_g_old));
        const int n_g_new = n_g_old + (is_add ? +1 : -1);

        /* Update Num for the single affected cell. */
        double N_plus = N_minus;

        if(has_kappa){
          if(r_star == kappa_code && r_star > 0){
            if(is_add){
              N_plus += lambda * pow(rlam, (double)n_gk_old);
            }else{
              N_plus -= lambda * pow(rlam, (double)(n_gk_old - 1));
            }
          }
        }else{
          if(r_star > 0){
            if(is_add){
              N_plus += lambda * pow(rlam, (double)n_gr_old);
            }else{
              N_plus -= lambda * pow(rlam, (double)(n_gr_old - 1));
            }
          }
        }

        const double D_plus = lambda * (1.0 - pow(rlam, (double)n_g_new));

        const double ratio_minus = (D_minus > 0.0) ? (N_minus / D_minus) : 0.0;
        const double ratio_plus  = (D_plus  > 0.0) ? (N_plus  / D_plus ) : 0.0;

        delta = ratio_plus - ratio_minus;

#if DEBUG_COV_MATCH_GW
        Rprintf("[cov_match_GW][by_group] i=%d j=%d lambda=%g N-=%g D-=%g N+=%g D+=%g delta=%g\n",
                i, j, lambda, N_minus, D_minus, N_plus, D_plus, delta);
#endif

      }else if(norm_mode == 2){
        /* Global normalisation:
         *   Δ = Δ_non_norm / [λ (1 - r_λ^{n1})]
         */
        const double Dglob = lambda * (1.0 - pow(rlam, (double)n1));
        delta = (Dglob > 0.0) ? (delta_non_norm / Dglob) : 0.0;

#if DEBUG_COV_MATCH_GW
        Rprintf("[cov_match_GW][global] i=%d j=%d lambda=%g Dglob=%g delta_non_norm=%g delta=%g\n",
                i, j, lambda, Dglob, delta_non_norm, delta);
#endif
      }

      CHANGE_STAT[j] += delta;
    }

toggle_apply:
    /* Temporarily apply this toggle so subsequent toggles see updated degrees. */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* 5) Undo temporary toggles to restore the original network state. */
  UNDO_PREVIOUS_TOGGLES(i);

  /* 6) Free buffers. */
  R_Free(stamp);
  R_Free(counts_buf);
  R_Free(codes_buf);
  R_Free(actors_buf);
}
