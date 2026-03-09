/**
 * @file changestat_cov_diff.c
 * @brief  Change statistic for the ERPM term `cov_diff` (multi-toggle form).
 * 
 * @author Jérémie Chichignoud
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
 *  `cov_diff`, which measures, for each group in the group mode, the
 *  dispersion of a numeric actor covariate over all k-subsets of actors
 *  inside that group.
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
 *    - A(g)  = set of actors adjacent to g (members of group g),
 *    - n_g   = |A(g)| = size of group g,
 *    - x_i   = numeric covariate value for actor i.
 *
 *  For a fixed integer k ≥ 2, define for each k-subset S ⊂ A(g),
 *
 *      D(S) = max_{i∈S} x_i - min_{i∈S} x_i.
 *
 *  The group-level contribution is:
 *
 *      T_k(g) = ∑_{S ⊂ A(g), |S| = k} D(S).
 *
 *  The global statistic is:
 *
 *      T_k(p; x)
 *        = ∑_g T_k(g),
 *
 *  where p encodes the bipartite membership structure.
 *
 *  Two normalisation modes are supported:
 *
 *    - by-group normalisation by the number of k-subsets:
 *
 *        T_k^{by\_group}(g) = T_k(g) / C(n_g, k),
 *
 *    - global normalisation by the group size:
 *
 *        T_k^{global}(g)   = T_k(g) / n_g.
 *
 *  The top-level statistic is then obtained by summing T_k(g),
 *  T_k^{by\_group}(g) or T_k^{global}(g) over all groups g.
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
 *  INPUT_PARAM layout (from InitErgmTerm.cov_diff)
 *  ------------------------------------------------------------
 *
 *  The R initialiser packs the parameters into INPUT_PARAM as:
 *
 *    INPUT_PARAM = c(
 *      n1,
 *      k,
 *      norm_mode,
 *      x[1:n1]
 *    )
 *
 *  where:
 *    - n1        = number of actors (size of the actor mode),
 *    - k         = size of subsets used in D(S),
 *    - norm_mode = 0 for raw T_k(g),
 *                  1 for per-group normalisation T_k(g)/C(n_g, k),
 *                  2 for global normalisation T_k(g)/n_g,
 *    - x[ ]      = numeric covariate on actors, length n1.
 *
 *  In C, this becomes:
 *
 *    INPUT_PARAM[0]     = n1
 *    INPUT_PARAM[1]     = k
 *    INPUT_PARAM[2]     = norm_mode
 *    INPUT_PARAM[3..]   = x[0..n1-1]
 *
 *  The term returns a single scalar statistic:
 *
 *    - N_CHANGE_STATS = 1,
 *    - CHANGE_STAT[0] is updated by the local Δ at each toggle.
 *
 *  ------------------------------------------------------------
 *  Local change under a toggle
 *  ------------------------------------------------------------
 *
 *  A membership toggle flips membership of one actor in one group:
 *    - addition  : actor becomes member of the group,
 *    - deletion  : actor leaves the group.
 *
 *  Only the affected group g can change its contribution:
 *
 *      Δ = T_after(g) - T_before(g),
 *
 *  where T(g) is either T_k(g), its per-group normalised version,
 *  or its global size-normalised version, depending on norm_mode.
 *  The helper function group_covdiff() recomputes the
 *  contribution for group g by:
 *
 *    1. reconstructing its membership in the actor mode (deduplicated),
 *    2. enumerating all k-subsets of its actors,
 *    3. computing D(S) = max(x_i) - min(x_i) for each subset S,
 *    4. summing over all subsets and optionally dividing by C(n_g, k)
 *       or by n_g.
 *
 *  ------------------------------------------------------------
 *  Multi-toggle / D_CHANGESTAT_FN (CRITICAL)
 *  ------------------------------------------------------------
 *
 *  This effect MUST support proposals decomposed into multiple toggles
 *  (swap/split/merge → a list of membership edge flips). When multiple
 *  toggles touch the same group, we must evaluate them in sequence under
 *  the correct intermediate state.
 *
 *  Therefore:
 *    - the change statistic is implemented using D_CHANGESTAT_FN,
 *    - toggles are processed sequentially,
 *    - after processing a toggle i, we temporarily apply it with
 *      TOGGLE_IF_MORE_TO_COME(i) so later toggles see updated degrees,
 *    - at the end, we restore the original network with UNDO_PREVIOUS_TOGGLES.
 *
 *  IMPORTANT:
 *  - We still compute "before/after" for each toggle via a *virtual* TOGGLE
 *    (apply, compute, undo) because group_covdiff() reconstructs neighbors
 *    from the current network state.
 *  - Then, TOGGLE_IF_MORE_TO_COME(i) applies the toggle *for real* only for
 *    i < ntoggles-1, so the intermediate state is correct.
 *
 *  ------------------------------------------------------------
 *  Complexity
 *  ------------------------------------------------------------
 *
 *  For a single group g of size n_g and fixed k:
 *    - neighbour reconstruction is O(n_g),
 *    - the number of k-subsets is C(n_g, k),
 *    - the total complexity is O(C(n_g, k) * k).
 *
 *  This implementation is intended for moderate group sizes or small k.
 *  There is no caching of per-group state across toggles.
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  The R initialiser (InitErgmTerm.cov_diff):
 *    - validates the numeric covariate on the actor mode,
 *    - sets k ≥ 2 and norm_mode ∈ {0,1,2},
 *    - packs n1, k, norm_mode and x into INPUT_PARAM,
 *    - sets emptynwstats and a single coef.name,
 *    - sets d_func = TRUE so ergm calls the D_ entrypoint.
 *
 *  ------------------------------------------------------------
 *  @example Usage (R)
 *  ------------------------------------------------------------
 *  @code{.r}
 *  library(ERPM)
 *
 *  # Example partition: 6 actors into 3 groups
 *  part <- c(1, 1, 2, 2, 3, 3)
 *
 *  # Numeric covariate on actors
 *  x <- c(1.0, 2.0, 0.5, 0.9, 1.5, 1.8)
 *
 *  # Raw cov_diff over all 2-subsets inside each group
 *  fit1 <- erpm(partition ~ cov_diff(attr = x, k = 2, normalized = FALSE))
 *  summary(fit1)
 *
 *  # Normalised version: average D(S) over all k-subsets per group
 *  fit2 <- erpm(partition ~ cov_diff(attr = x, k = 2, normalized = "by_group"))
 *  summary(fit2)
 *
 *  # Global version: average D(S) per actor in each group
 *  fit3 <- erpm(partition ~ cov_diff(attr = x, k = 2, normalized = "global"))
 *  summary(fit3)
 *
 *  # Internally, the MCMC proposal may produce multi-toggle moves. This effect
 *  # is implemented as a D_ changestat and remains consistent under such moves.
 *  @endcode
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"      /* R_Calloc / R_Free */
#include <R_ext/Print.h>
#include <math.h>

/**
 * @def DEBUG_COV_DIFF
 * @brief Enable verbose debugging output for ::d_cov_diff.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - group sizes for cov_diff,
 *  - values of k and norm_mode,
 *  - local contribution of the affected group,
 *  - multi-toggle traces (ntoggles, sequential deltas).
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_COV_DIFF 0

/**
 * @def UNUSED_VARIABLE
 * @brief Utility macro to explicitly mark unused parameters.
 *
 * @param x Parameter or variable that is intentionally unused in a
 *          particular compilation unit or function.
 */
#define UNUSED_VARIABLE(x) (void)(x)

/* -------------------------------------------------------------------------- */
/* Helper: recursive enumeration of k-subsets                                 */
/* -------------------------------------------------------------------------- */

/**
 * @brief Recursively sum D(S) = max(x_i) - min(x_i) over all k-subsets S.
 *
 * @details
 *  This function performs a depth-first enumeration of all combinations
 *  of size k among ng actor indices stored in @p idxs[0..ng-1]. The
 *  current partial combination is stored in @p comb, as positions in
 *  idxs[]. When a complete combination of size k is formed, it computes:
 *
 *      D(S) = max_{i∈S} x_i - min_{i∈S} x_i
 *
 *  and adds it to the accumulator @p acc.
 *
 *  The recursion uses:
 *    - @p pos   = current depth (number of elements chosen so far),
 *    - @p start = first index in idxs[] that may be chosen at this depth.
 *
 * @param pos    Current depth in the combination (0..k).
 * @param start  Start index in @p idxs for the next choice.
 * @param k      Target subset size.
 * @param ng     Number of available actors in the group.
 * @param idxs   Array of actor indices (0-based indices into x[]).
 * @param x      Numeric covariate vector for actors.
 * @param comb   Working array of length k storing chosen positions in idxs[].
 * @param acc    Pointer to the accumulator for the sum of D(S) over all S.
 */
static void sum_D_rec(int pos, int start,
                      int k, int ng,
                      const int *idxs,
                      const double *x,
                      int *comb,
                      double *acc){

  if(pos == k){
    /* We have a complete combination in comb[0..k-1]. */
    double xmin = 0.0, xmax = 0.0;
    int first = 1;
    for(int t = 0; t < k; t++){
      int idx = idxs[ comb[t] ];   /* 0-based index into x[] */
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

  /* Choose the next element among idxs[start..ng-1]. */
  for(int i = start; i <= ng - (k - pos); i++){
    comb[pos] = i;
    sum_D_rec(pos + 1, i + 1, k, ng, idxs, x, comb, acc);
  }
}

/* -------------------------------------------------------------------------- */
/* Helper: cov_diff contribution for a single group                           */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute the cov_diff contribution T_k(g) for a single group vertex.
 *
 * @details
 *  For a given group vertex g in the group mode, this function:
 *
 *    1. Reconstructs the set of member actors A(g) by:
 *         - traversing outgoing edges (group → actor),
 *         - traversing incoming edges (actor → group),
 *         - deduplicating actors using a "seen" array.
 *
 *    2. If n_g = |A(g)| < k, returns 0 (no valid k-subsets).
 *
 *    3. Otherwise, enumerates all k-subsets S ⊂ A(g) via sum_D_rec(),
 *       and accumulates:
 *
 *           T_k(g) = ∑_S D(S),
 *
 *       where D(S) = max_{i∈S} x_i - min_{i∈S} x_i.
 *
 *    4. If norm_mode == 1, divides T_k(g) by C(n_g, k) to obtain the
 *       group-level average over all k-subsets:
 *
 *           T_k^{by\_group}(g) = T_k(g) / C(n_g, k).
 *
 *       If norm_mode == 2, divides T_k(g) by n_g to obtain the group-level
 *       average per actor:
 *
 *           T_k^{global}(g) = T_k(g) / n_g.
 *
 *  The actor covariate x is indexed as:
 *    - x[idx] is the covariate of actor (idx + 1) in the actor mode,
 *      for idx in 0..n1-1.
 *
 * @param g          Group vertex in the group mode.
 * @param n1         Number of actors (size of the actor mode).
 * @param k          Subset size used in cov_diff.
 * @param norm_mode  0 for raw T_k(g),
 *                   1 for by-group normalisation T_k(g)/C(n_g, k),
 *                   2 for global normalisation T_k(g)/n_g.
 * @param x          Pointer to numeric covariate values (length ≥ n1).
 * @param nwp        Pointer to the network-plus workspace.
 *
 * @return The contribution of group g:
 *         - T_k(g)                  if norm_mode == 0,
 *         - T_k(g)/C(n_g, k)        if norm_mode == 1,
 *         - T_k(g)/n_g              if norm_mode == 2,
 *         - 0 if n_g < k.
 */
static double group_covdiff(Vertex g,
                            int n1,
                            int k,
                            int norm_mode,
                            const double *x,
                            Network *nwp){

  Edge e;
  Vertex h;
  int ng = 0;  /* number of actors in this group */

  /* "seen" marks actors already counted to avoid double counting. */
  unsigned char *seen = (unsigned char *)R_Calloc(n1, unsigned char);  /* 0-initialised */
  /* idxs holds 0-based indices into x[] for the actors in this group. */
  int *idxs          = (int *)R_Calloc(n1, int);

  /* OUT-neighbours: group → actor. */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        idxs[ng++] = idx;
      }
    }
  }

  /* IN-neighbours: actor → group. */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        idxs[ng++] = idx;
      }
    }
  }

  /* Reset seen marks for the actors actually present in the group. */
  if(ng > 0){
    for(int i = 0; i < ng; i++){
      seen[ idxs[i] ] = 0;
    }
  }
  R_Free(seen);

  /* Not enough actors to form a k-subset. */
  if(ng < k){
    #if DEBUG_COV_DIFF
      Rprintf("[cov_diff][group_covdiff] g=%d ng=%d < k=%d -> 0\n",
              (int)g, ng, k);
    #endif
    R_Free(idxs);
    return 0.0;
  }

  /* Sum D(S) over all k-subsets S of the group. */
  double sumD = 0.0;
  int *comb = (int *)R_Calloc(k, int);
  sum_D_rec(0, 0, k, ng, idxs, x, comb, &sumD);
  R_Free(comb);

  double res = sumD;

  /* Normalisation, if requested. */
  if(norm_mode == 1){
    /* By-group normalisation by C(n_g, k). */
    double denom = CHOOSE(ng, k);
    if(denom > 0.0){
      res /= denom;
    }else{
      res = 0.0;
    }
  }else if(norm_mode == 2){
    /* Global normalisation by group size n_g. */
    double denom = (double)ng;
    if(denom > 0.0){
      res /= denom;
    }else{
      res = 0.0;
    }
  }

  #if DEBUG_COV_DIFF
    Rprintf("[cov_diff][group_covdiff] g=%d ng=%d k=%d norm_mode=%d -> res=%g\n",
            (int)g, ng, k, norm_mode, res);
  #endif

  R_Free(idxs);
  return res;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_diff (multi-toggle)                                  */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_diff` (multi-toggle).
 *
 * @details
 *  This is the \pkg{ergm} change-statistic function registered as
 *  ::d_cov_diff via ::D_CHANGESTAT_FN. It computes the local change
 *  Δ in the cov_diff statistic for a multi-toggle proposal, i.e. a list
 *  of membership toggles between actors and groups.
 *
 *  The layout of INPUT_PARAM is:
 *
 *    INPUT_PARAM[0]     = n1         (number of actors)
 *    INPUT_PARAM[1]     = k          (subset size)
 *    INPUT_PARAM[2]     = norm_mode  (0 raw, 1 by-group, 2 global)
 *    INPUT_PARAM[3..]   = x[0..n1-1] (numeric covariate on actors)
 *
 *  For each toggle i (processed sequentially under the intermediate state):
 *    1. Identify the actor vertex and the group vertex using the
 *       bipartite boundary between actor mode and group mode.
 *    2. Evaluate the group contribution BEFORE the toggle via group_covdiff().
 *    3. Apply a *virtual* toggle (TOGGLE), evaluate AFTER, then undo it.
 *    4. Accumulate Δ_i = AFTER - BEFORE into CHANGE_STAT[0].
 *    5. Temporarily apply the toggle (TOGGLE_IF_MORE_TO_COME(i)) so that
 *       later toggles see updated membership degrees and neighborhoods.
 *
 *  At the end, undo the temporary toggles (UNDO_PREVIOUS_TOGGLES).
 *
 *  IMPORTANT:
 *  - The edgestate is not required because we explicitly do a virtual TOGGLE
 *    to compute the "after" state for the current toggle.
 */
D_CHANGESTAT_FN(d_cov_diff){

#if DEBUG_COV_DIFF
  static int seen = 0;
  if(ntoggles > 1 && seen < 10){
    Rprintf("[cov_diff] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen++;
  }
#endif

  /* 1) Reset the output buffer for THIS proposal. */
  ZERO_ALL_CHANGESTATS();

  /* 2) Read inputs from INPUT_PARAM. */
  const double *ip       = INPUT_PARAM;
  const int n1           = (int)ip[0];   /* number of actors (actor mode) */
  const int k            = (int)ip[1];   /* subset size for cov_diff */
  const int norm_mode    = (int)ip[2];   /* 0 raw, 1 by-group, 2 global */
  const double *x        = ip + 3;       /* actor covariate values */

#if DEBUG_COV_DIFF
  Rprintf("[cov_diff] n1=%d k=%d norm_mode=%d | BIPARTITE=%d\n",
          n1, k, norm_mode, (int)BIPARTITE);
#endif

  /* 3) Process toggles sequentially under the intermediate state. */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex a = TAIL(i);
    Vertex b = HEAD(i);

    /* Identify actor and group vertices (actor <= n1, group > n1). */
    Vertex actor = (a <= (Vertex)n1) ? a : b;
    Vertex group = (a <= (Vertex)n1) ? b : a;

    UNUSED_VARIABLE(actor);

#if DEBUG_COV_DIFF
    if(actor > (Vertex)n1 || group <= (Vertex)n1){
      Rprintf("[cov_diff][WARN] toggle #%d endpoints not (actor,group): tail=%d head=%d | actor=%d group=%d (n1=%d)\n",
              i, (int)a, (int)b, (int)actor, (int)group, n1);
    }
#endif

    /* BEFORE under current intermediate state. */
    double F_before = group_covdiff(group, n1, k, norm_mode, x, nwp);

    /* Virtual toggle to get AFTER. */
    TOGGLE(a, b);
    double F_after  = group_covdiff(group, n1, k, norm_mode, x, nwp);
    TOGGLE(a, b);

    double delta = (F_after - F_before);
    CHANGE_STAT[0] += delta;

#if DEBUG_COV_DIFF
    Rprintf("[D:d_cov_diff] i=%d tail=%d head=%d | group=%d | before=%g after=%g | Δ=%g | cumul=%g\n",
            i, (int)a, (int)b, (int)group, F_before, F_after, delta, CHANGE_STAT[0]);
#endif

    /* Apply this toggle for real if more toggles remain (intermediate state). */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* 4) Restore the original network state (undo temporary toggles). */
  UNDO_PREVIOUS_TOGGLES(i);
}
