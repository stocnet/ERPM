// ============================================================================
// File    : src/changestat_cov_fulldiff.c
// Purpose : Change statistic for the ERPM term `cov_fulldiff` (MULTI-TOGGLE form)
// Project : ERPM / ERGM extensions
// ============================================================================

/**
 * @file changestat_cov_fulldiff.c
 * @brief  Change statistic for the ERPM term `cov_fulldiff` (multi-toggle form).
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
 *  `cov_fulldiff`, which measures, for each group in the group mode, the
 *  within-group dispersion of a numeric actor covariate via the range:
 *
 *      range_g(x) = x_g^max - x_g^min,
 *
 *  with an optional filter on group sizes.
 *
 *  ------------------------------------------------------------
 *  IMPORTANT (multi-toggle / D_CHANGESTAT_FN)
 *  ------------------------------------------------------------
 *
 *  This implementation is MULTI-TOGGLE safe:
 *   - It is declared with D_CHANGESTAT_FN(d_cov_fulldiff).
 *   - It supports proposals containing ntoggles >= 1 (swap/split/merge decomposed
 *     into a list of edge toggles).
 *
 *  Why it matters:
 *   - When ntoggles > 1, several toggles may affect the same group.
 *   - The statistic must be computed consistently with intermediate states.
 *   - Therefore we evaluate toggles sequentially, *temporarily applying* each
 *     toggle so later toggles see updated memberships.
 *
 *  Design choice (robustness over cleverness):
 *   - For each toggle i, we recompute the affected group contribution BEFORE and
 *     AFTER applying that toggle, using the current intermediate network state.
 *   - After processing all toggles, we undo ALL toggles in reverse order to restore
 *     the original network state.
 *
 *  This approach is simple and correct. It may be more expensive than a fully
 *  incremental min/max update, but avoids subtle bugs and stays local.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  The bipartite network represents a partition:
 *    - actor mode  = vertices representing actors (individuals),
 *    - group mode  = vertices representing structural groups.
 *
 *  Membership is encoded by edges between actors and groups.
 *  For a group vertex g in the group mode:
 *
 *    - A(g)       = set of actors adjacent to g (members of group g),
 *    - n_g        = |A(g)| = size of group g,
 *    - x_i        = numeric covariate value for actor i,
 *    - x_g^min    = min_{i in A(g)} x_i,
 *    - x_g^max    = max_{i in A(g)} x_i.
 *
 *  Let S be an optional set of allowed group sizes. The global statistic is:
 *
 *      T(p; x, S)
 *        = ∑_g 1[n_g ∈ S] * (x_g^max - x_g^min),
 *
 *  where:
 *    - p describes the bipartite membership pattern (partition),
 *    - the indicator 1[n_g ∈ S] is 1 if S is empty or n_g is in S,
 *      and 0 otherwise.
 *
 *  Groups of size 0 or 1 do not contribute because their internal dispersion
 *  is undefined or zero; this is handled explicitly in the implementation.
 *
 *  ------------------------------------------------------------
 *  Bipartite structure (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  At the C level, the bipartite structure is encoded via the boundary
 *  BIPARTITE:
 *
 *    - actor mode  : vertices 1 .. n1, where n1 = BIPARTITE,
 *    - group mode  : vertices > n1, representing groups.
 *
 *  This change statistic assumes:
 *    - the number of actors n1 (actor mode) is stored in INPUT_PARAM[0],
 *    - each membership toggle connects exactly one actor (vertex ≤ n1)
 *      and one group (vertex > n1).
 *
 *  The group vertex is detected as the endpoint in the group mode
 *  (vertex index > n1), and the actor vertex as the endpoint in
 *  the actor mode (vertex index ≤ n1).
 *
 *  ------------------------------------------------------------
 *  INPUT_PARAM layout (from InitErgmTerm.cov_fulldiff)
 *  ------------------------------------------------------------
 *
 *  The R initialiser packs the parameters into INPUT_PARAM as:
 *
 *    INPUT_PARAM = c(
 *      n1,
 *      L,
 *      sizes[1:L],
 *      x[1:n1]
 *    )
 *
 *  where:
 *    - n1         = number of actors (size of the actor mode),
 *    - L          = length of the size filter S,
 *    - sizes[ ]   = allowed group sizes (stored as doubles and cast to int),
 *    - x[ ]       = numeric covariate on actors, length n1.
 *
 *  In C, this becomes:
 *
 *    INPUT_PARAM[0]        = n1
 *    INPUT_PARAM[1]        = L
 *    INPUT_PARAM[2..1+L]   = sizes[0..L-1]
 *    INPUT_PARAM[2+L..]    = x[0..n1-1]
 *
 *  The term returns a single scalar statistic:
 *
 *    - N_CHANGE_STATS = 1,
 *    - CHANGE_STAT[0] is updated by the local Δ aggregated over toggles.
 *
 *  ------------------------------------------------------------
 *  Local change under a toggle
 *  ------------------------------------------------------------
 *
 *  A toggle flips membership of one actor in one group:
 *    - addition  : actor becomes member of the group,
 *    - deletion  : actor leaves the group.
 *
 *  Only the affected group g can change its contribution:
 *
 *      Δ = range_after(g) * 1[n_g_after ∈ S]
 *        - range_before(g) * 1[n_g_before ∈ S].
 *
 *  In multi-toggle mode, "before" and "after" are evaluated with respect to the
 *  current intermediate state right before and right after applying toggle i.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"      /* R_Calloc/R_Free */
#include <R_ext/Print.h>

/**
 * @def DEBUG_COV_FULLDIFF
 * @brief Enable verbose debugging output for ::d_cov_fulldiff.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - group sizes and ranges before and after each toggle,
 *  - inclusion or exclusion by size filter,
 *  - per-toggle Δ and cumulative Δ.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_COV_FULLDIFF 0   /* set to 0 to disable debug output */

/**
 * @def DEBUG_COV_FULLDIFF_MULTITOGGLE_PING
 * @brief Print a small ping when ntoggles>1 (even if DEBUG_COV_FULLDIFF=0).
 *
 * This is useful to confirm that the D_ entrypoint is used and that the proposal
 * path can generate multi-toggle moves.
 */
#define DEBUG_COV_FULLDIFF_MULTITOGGLE_PING 1

/**
 * @def UNUSED_WARNING
 * @brief Utility macro to explicitly mark unused parameters.
 */
#define UNUSED_WARNING(x) (void)x

/* -------------------------------------------------------------------------- */
/* Helper: group size filter                                                  */
/* -------------------------------------------------------------------------- */

static inline int in_sizes(int n, int L, const double *sizes){
  if(L == 0) return 1;
  for(int i = 0; i < L; i++){
    if((int)sizes[i] == n) return 1;
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Helper: compute range contribution for a single group                      */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute the covariate range contribution for a single group vertex.
 *
 * @details
 *  This function reconstructs actor membership for group g by traversing both
 *  out- and in-edges, deduplicating actors through a "seen" array.
 *
 *  It returns:
 *   - 0 if ng <= 1 (empty or singleton group),
 *   - 0 if size filter excludes ng,
 *   - otherwise (xmax - xmin) over actors in the group.
 */
static double group_range(Vertex g,
                          int n1, int L, const double *sizes,
                          const double *x,
                          Network *nwp){

  int ng = 0;
  double xmin = 0.0;
  double xmax = 0.0;
  unsigned char first = 1;

  unsigned char *seen = (unsigned char*)R_Calloc(n1, unsigned char);

  Vertex h;
  Edge e;

  /* OUT-neighbours: group → actor edges. */
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

  /* IN-neighbours: actor → group edges. */
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

  /* Free the buffer. (No need to reset marks: seen is local and freed.) */
  R_Free(seen);

  if(ng <= 1){
    return 0.0;
  }

  if(!in_sizes(ng, L, sizes)){
    return 0.0;
  }

  return (xmax - xmin);
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_fulldiff (multi-toggle)                              */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_fulldiff` (multi-toggle).
 *
 * @details
 *  - Declared as D_CHANGESTAT_FN(d_cov_fulldiff).
 *  - Processes toggles sequentially, applying each toggle to update the
 *    intermediate state.
 *  - Undoes all toggles in reverse order at the end.
 *
 *  This avoids relying on TOGGLE_IF_MORE_TO_COME/UNDO_PREVIOUS_TOGGLES macros,
 *  and remains correct even if:
 *    - several toggles hit the same group,
 *    - the same dyad is toggled multiple times in the proposal.
 */
D_CHANGESTAT_FN(d_cov_fulldiff){

#if DEBUG_COV_FULLDIFF_MULTITOGGLE_PING
  if(ntoggles > 1){
    static int seen = 0;
    if(seen < 25){
      Rprintf("[cov_fulldiff] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
      seen++;
    }
  }
#endif

  /* 1) Reset output buffer for the whole proposal. */
  ZERO_ALL_CHANGESTATS();

  /* 2) Read inputs from INPUT_PARAM. */
  const double *ip    = INPUT_PARAM;
  const int n1        = (int)ip[0];
  const int L         = (int)ip[1];
  const double *sizes = ip + 2;
  const double *x     = ip + 2 + L;

#if DEBUG_COV_FULLDIFF
  Rprintf("[cov_fulldiff][D] n1=%d L=%d ntoggles=%d\n", n1, L, (int)ntoggles);
#endif

  /* 3) Process each toggle i sequentially on the intermediate state. */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex t = TAIL(i);
    Vertex h = HEAD(i);

    /* Identify actor/group vertices using the bipartite split. */
    Vertex actor = (t <= (Vertex)n1) ? t : h;
    Vertex group = (t <= (Vertex)n1) ? h : t;

    /* Safety: if proposal is malformed (no group endpoint), just ignore.
     * (Should not happen for membership toggles in ERPM usage.) */
    if(group <= (Vertex)n1){
#if DEBUG_COV_FULLDIFF
      Rprintf("[cov_fulldiff][D][WARN] toggle #%d has no group endpoint: tail=%d head=%d (n1=%d)\n",
              i, (int)t, (int)h, n1);
#endif
      continue;
    }

    /* Contribution BEFORE applying this toggle (current intermediate state). */
    double F_before = group_range(group, n1, L, sizes, x, nwp);

    /* Apply this toggle on the intermediate network state. */
    TOGGLE(t, h);

    /* Contribution AFTER applying this toggle (updated intermediate state). */
    double F_after  = group_range(group, n1, L, sizes, x, nwp);

    /* Accumulate local delta. */
    CHANGE_STAT[0] += (F_after - F_before);

#if DEBUG_COV_FULLDIFF
    Rprintf("[cov_fulldiff][D] i=%d tail=%d head=%d | actor=%d group=%d | before=%g after=%g | d=%g | cumul=%g\n",
            i, (int)t, (int)h, (int)actor, (int)group, F_before, F_after, (F_after - F_before), CHANGE_STAT[0]);
#endif

    UNUSED_WARNING(actor);
  }

  /* 4) Undo ALL toggles in reverse order to restore original network state. */
  for(int j = (int)ntoggles - 1; j >= 0; --j){
    Vertex t = TAIL(j);
    Vertex h = HEAD(j);
    TOGGLE(t, h);
  }
}
