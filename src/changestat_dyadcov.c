/* =============================================================================
 * File    : changestat_dyadcov.c
 * Purpose : Change statistic for the ERPM term `dyadcov` (multi-toggle form).
 * Project : ERPM / ERGM extensions
 * =============================================================================
 *
 * This file is the MULTI-TOGGLE (D_CHANGESTAT_FN) counterpart of the historical
 * one-toggle implementation. It preserves the exact same statistic definition
 * but makes the change-statistic safe under proposals that consist of multiple
 * toggles (swap/split/merge decomposed into a list of edge toggles).
 *
 * ---------------------------------------------------------------------------
 * IMPORTANT (multi-toggle / D_CHANGESTAT_FN)
 * ---------------------------------------------------------------------------
 *
 * - ergm can propose moves containing multiple edge toggles.
 * - If the compiled changestat is one-toggle (C_CHANGESTAT_FN), ergm will call
 *   it repeatedly, once per toggle, on intermediate states that are managed by
 *   the engine.
 * - If you implement D_CHANGESTAT_FN, YOU must handle the full list of toggles
 *   consistently and return the aggregated change over the whole proposal.
 *
 * Design rule used here (same as squared_sizes):
 * - Process toggles sequentially.
 * - For each toggle i:
 *     1) compute the group contribution BEFORE the toggle under the current
 *        intermediate state (which already includes previous toggles applied),
 *     2) apply a VIRTUAL toggle to compute the AFTER contribution,
 *     3) undo the virtual toggle (back to the current intermediate state),
 *     4) accumulate Δ_i = after - before (with the chosen normalisation),
 *     5) temporarily apply the toggle for real IF there are more toggles coming,
 *        so that subsequent degrees/memberships are consistent.
 * - At the end: undo the temporarily-applied toggles and restore the original
 *   network state.
 *
 * This guarantees:
 * - Correct behaviour when multiple toggles hit the same group vertex.
 * - Correct behaviour when toggles interact (because “current state” is updated
 *   as we advance in the toggle list).
 *
 * ---------------------------------------------------------------------------
 * Statistic definition (unchanged)
 * ---------------------------------------------------------------------------
 *
 * See the long header in the user-provided reference version. In short:
 * - Bipartite network with:
 *     actor mode = 1..n1, group mode = n1+1..N
 * - Dyadic covariate matrix Z on actors (n1 x n1), column-major (R order).
 * - For a group g, for fixed clique size k>=2:
 *
 *     S_g^{(k)}(Z) = ∑_{C ∈ C_k(g)} ∏_{i<j∈C} (z_ij + z_ji).
 *
 * - Global statistics:
 *     norm_mode 0: ∑_g S_g
 *     norm_mode 1: ∑_g 1[n_g>=k] (1/n_g) S_g
 *     norm_mode 2: ∑_g 1[n_g>=k] (1/choose(n_g,k)) S_g
 *
 * ---------------------------------------------------------------------------
 * INPUT_PARAM layout (unchanged)
 * ---------------------------------------------------------------------------
 *   INPUT_PARAM = c(n1, k, norm_mode, as.vector(Z))
 *     ip[0]  = n1
 *     ip[1]  = k
 *     ip[2]  = norm_mode
 *     ip[3+] = Z (length n1*n1), column-major
 * ------------------------------------------------------------------------- */

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>

/**
 * @def DEBUG_DYADCOV
 * @brief Enable verbose debugging output for ::d_dyadcov.
 *
 * Set this macro to 1 to print detailed information to the R console during
 * summary() or MCMC runs:
 *  - current parameter values (n1, k, norm_mode),
 *  - per-toggle group sizes before/after,
 *  - S_before, S_after, delta for each toggle,
 *  - cumulative CHANGE_STAT.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_DYADCOV 0

/**
 * @def UNUSED_WARNING
 * @brief Macro to explicitly mark a parameter as intentionally unused.
 */
#define UNUSED_WARNING(x) (void)x

/* -------------------------------------------------------------------------- */
/* Clique enumeration: sum_cliques_k                                           */
/* -------------------------------------------------------------------------- */
/**
 * @brief Sum over all k-cliques of actors inside a given group.
 *
 * Given actor indices (1-based, in the actor mode) of length ng, enumerates all
 * k-subsets and sums:
 *    ∏_{i<j in clique} (z_ij + z_ji)
 *
 * Z is indexed as Z[(j-1)*n1 + (i-1)] (column-major, R convention).
 */
static double sum_cliques_k(const int *actors,
                            int ng, int k,
                            int n1,
                            const double *Z){

  if(k > ng) return 0.0;

  int *comb = (int*)R_Calloc(k, int);
  for(int i = 0; i < k; i++) comb[i] = i;

  double total = 0.0;

  while(1){
    double prod = 1.0;

    for(int p = 0; p < k; p++){
      int idx_i = actors[ comb[p] ]; /* 1..n1 */
      int row   = idx_i - 1;         /* 0..n1-1 */

      for(int q = p + 1; q < k; q++){
        int idx_j  = actors[ comb[q] ];
        int col    = idx_j - 1;
        int idx_ij = col * n1 + row; /* z_ij */
        int idx_ji = row * n1 + col; /* z_ji */
        prod *= (Z[idx_ij] + Z[idx_ji]);
      }
    }

    total += prod;

    int pos = k - 1;
    while(pos >= 0 && comb[pos] == (ng - k + pos)) pos--;
    if(pos < 0) break;

    comb[pos]++;
    for(int j = pos + 1; j < k; j++){
      comb[j] = comb[j - 1] + 1;
    }
  }

  R_Free(comb);
  return total;
}

/* -------------------------------------------------------------------------- */
/* Group-level functional: group_dyadcov_k                                     */
/* -------------------------------------------------------------------------- */
/**
 * @brief Compute S_g^{(k)}(Z) for a given group vertex g.
 *
 * Collects actor neighbours of g via OUT and IN edges (deduplicated), counts ng,
 * and returns 0 if ng<k; else returns sum_cliques_k over actor indices.
 */
static double group_dyadcov_k(Vertex g,
                              int n1, int k,
                              const double *Z,
                              Network *nwp,
                              int *n_g_out){

  unsigned char *seen = (unsigned char*)R_Calloc(n1, unsigned char);
  Vertex h;
  Edge e;

  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1) seen[(int)h - 1] = 1;
  }
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1) seen[(int)h - 1] = 1;
  }

  int ng = 0;
  for(int i = 0; i < n1; i++) if(seen[i]) ng++;
  if(n_g_out) *n_g_out = ng;

  if(ng < k){
#if DEBUG_DYADCOV
    Rprintf("[dyadcov][group] g=%d ng=%d < k=%d -> 0\n", (int)g, ng, k);
#endif
    R_Free(seen);
    return 0.0;
  }

  int *actors = (int*)R_Calloc(ng, int);
  int idx = 0;
  for(int i = 0; i < n1; i++) if(seen[i]) actors[idx++] = i + 1;

  double sum = sum_cliques_k(actors, ng, k, n1, Z);

#if DEBUG_DYADCOV
  Rprintf("[dyadcov][group] g=%d ng=%d k=%d -> sum=%g\n", (int)g, ng, k, sum);
#endif

  R_Free(actors);
  R_Free(seen);
  return sum;
}

/* -------------------------------------------------------------------------- */
/* Normalisation helper (group-level)                                          */
/* -------------------------------------------------------------------------- */
static inline double normalise_group_sum(double S, int ng, int k, int norm_mode){
  if(norm_mode == 0) return S;

  if(ng < k || ng <= 0) return 0.0;

  if(norm_mode == 1){
    /* "global" 1/ng */
    return S / (double)ng;
  }

  if(norm_mode == 2){
    /* "by_group" 1/choose(ng,k) */
    double denom = CHOOSE(ng, k);
    if(denom <= 0.0) return 0.0;
    return S / denom;
  }

  /* unknown mode: be safe */
  return S;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: dyadcov (multi-toggle)                                    */
/* -------------------------------------------------------------------------- */
/**
 * @brief Multi-toggle change statistic for dyadcov (D_CHANGESTAT_FN).
 *
 * Returns the aggregated Δ over the full list of toggles.
 *
 * NOTE:
 * - We do NOT use `edgestate` from the one-toggle signature.
 * - We rely on virtual toggling + sequential application to keep state consistent.
 */
D_CHANGESTAT_FN(d_dyadcov){

  ZERO_ALL_CHANGESTATS();

  /* Decode INPUT_PARAM: [n1, k, norm_mode, Z...] */
  const double *ip        = INPUT_PARAM;
  const int     n1        = (int)ip[0];
  const int     k         = (int)ip[1];
  const int     norm_mode = (int)ip[2];
  const double *Z         = ip + 3;

#if DEBUG_DYADCOV
  static int seen_multi = 0;
  if(ntoggles > 1 && seen_multi < 10){
    Rprintf("[dyadcov] MULTI-TOGGLE ntoggles=%d | n1=%d k=%d mode=%d\n",
            (int)ntoggles, n1, k, norm_mode);
    seen_multi++;
  }
#endif

  /* Degenerate k */
  if(k < 2){
    CHANGE_STAT[0] = 0.0;
    return;
  }

  /* Process toggles sequentially */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex t = TAIL(i);
    Vertex h = HEAD(i);

    /* Identify the affected group vertex (>n1). */
    Vertex group = (t > (Vertex)n1) ? t : h;

    /* Sanity: in a valid bipartite membership toggle, one endpoint is group. */
#if DEBUG_DYADCOV
    if(group <= (Vertex)n1){
      Rprintf("[dyadcov][WARN] toggle #%d has no group endpoint: tail=%d head=%d (n1=%d)\n",
              i, (int)t, (int)h, n1);
    }
#endif

    int ng_before = 0, ng_after = 0;

    /* Contribution BEFORE (current intermediate state) */
    double S_before = group_dyadcov_k(group, n1, k, Z, nwp, &ng_before);
    double T_before = normalise_group_sum(S_before, ng_before, k, norm_mode);

    /* Virtual toggle to compute AFTER */
    TOGGLE(t, h);
    double S_after = group_dyadcov_k(group, n1, k, Z, nwp, &ng_after);
    double T_after = normalise_group_sum(S_after, ng_after, k, norm_mode);
    TOGGLE(t, h);

    double delta = T_after - T_before;
    CHANGE_STAT[0] += delta;

#if DEBUG_DYADCOV
    Rprintf("[dyadcov][D] i=%d tail=%d head=%d group=%d | ng %d->%d | "
            "S %g->%g | T %g->%g | Δ=%g | cumul=%g\n",
            i, (int)t, (int)h, (int)group,
            ng_before, ng_after,
            S_before, S_after,
            T_before, T_after,
            delta, CHANGE_STAT[0]);
#endif

    /* Apply toggle for subsequent toggles (except the last one). */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* Restore original state */
  UNDO_PREVIOUS_TOGGLES(i);
}
