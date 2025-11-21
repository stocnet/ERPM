/**
 * @file
 * @brief  Change statistic for the ERPM term `cov_diff` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
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
 *  where p encodes the bipartite membership structure. An optional
 *  per-group normalisation divides T_k(g) by C(n_g, k), the number
 *  of k-subsets of A(g):
 *
 *      T_k^norm(g) = T_k(g) / C(n_g, k),
 *
 *  so that:
 *
 *      T_k^norm(p; x) = ∑_g T_k^norm(g).
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
 *      norm_flag,
 *      x[1:n1]
 *    )
 *
 *  where:
 *    - n1        = number of actors (size of the actor mode),
 *    - k         = size of subsets used in D(S),
 *    - norm_flag = 0 for raw sum T_k(g),
 *                  1 for per-group normalisation T_k(g)/C(n_g, k),
 *    - x[ ]      = numeric covariate on actors, length n1.
 *
 *  In C, this becomes:
 *
 *    INPUT_PARAM[0]     = n1
 *    INPUT_PARAM[1]     = k
 *    INPUT_PARAM[2]     = norm_flag
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
 *  A single toggle flips membership of one actor in one group:
 *    - addition  : actor becomes member of the group,
 *    - deletion  : actor leaves the group.
 *
 *  Only the affected group g can change its contribution:
 *
 *      Δ = T_after(g) - T_before(g),
 *
 *  where T(g) is either T_k(g) or its normalised version, depending
 *  on norm_flag. The helper function group_covdiff() recomputes the
 *  contribution for group g by:
 *
 *    1. reconstructing its membership in the actor mode (deduplicated),
 *    2. enumerating all k-subsets of its actors,
 *    3. computing D(S) = max(x_i) - min(x_i) for each subset S,
 *    4. summing over all subsets and optionally dividing by C(n_g, k).
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
 *    - sets k ≥ 2 and norm_flag ∈ {0,1},
 *    - packs n1, k, norm_flag and x into INPUT_PARAM,
 *    - sets emptynwstats and a single coef.name.
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
 *  fit2 <- erpm(partition ~ cov_diff(attr = x, k = 2, normalized = TRUE))
 *  summary(fit2)
 *
 *  # Internally, each membership toggle between an actor and a group
 *  # calls c_cov_diff(), which:
 *  #   - recomputes the cov_diff contribution of the affected group
 *  #     before and after a virtual toggle,
 *  #   - applies the optional group-level normalisation,
 *  #   - updates CHANGE_STAT[0] by the difference.
 *  @endcode
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"      /* Calloc / Free */
#include <R_ext/Print.h>
#include <math.h>

/**
 * @def DEBUG_COV_DIFF
 * @brief Enable verbose debugging output for ::c_cov_diff.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - group sizes for cov_diff,
 *  - values of k and norm_flag,
 *  - local contribution of the affected group.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_COV_DIFF 0

/**
 * @def UNUSED_WARNING
 * @brief Utility macro to explicitly mark unused parameters.
 *
 * @param x Parameter or variable that is intentionally unused in a
 *          particular compilation unit or function.
 */
#define UNUSED_WARNING(x) (void)(x)

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
 *    4. If norm_flag == 1, divides T_k(g) by C(n_g, k) to obtain the
 *       group-level average over all k-subsets:
 *
 *           T_k^norm(g) = T_k(g) / C(n_g, k).
 *
 *  The actor covariate x is indexed as:
 *    - x[idx] is the covariate of actor (idx + 1) in the actor mode,
 *      for idx in 0..n1-1.
 *
 * @param g          Group vertex in the group mode.
 * @param n1         Number of actors (size of the actor mode).
 * @param k          Subset size used in cov_diff.
 * @param norm_flag  0 for raw T_k(g), 1 for per-group normalisation.
 * @param x          Pointer to numeric covariate values (length ≥ n1).
 * @param nwp        Pointer to the network-plus workspace.
 *
 * @return The contribution of group g:
 *         - T_k(g) if norm_flag == 0,
 *         - T_k(g)/C(n_g, k) if norm_flag == 1,
 *         - 0 if n_g < k.
 */
static double group_covdiff(Vertex g,
                            int n1,
                            int k,
                            int norm_flag,
                            const double *x,
                            Network *nwp){

  Edge e;
  Vertex h;
  int ng = 0;  /* number of actors in this group */

  /* "seen" marks actors already counted to avoid double counting. */
  unsigned char *seen = (unsigned char *)Calloc(n1, unsigned char);  /* 0-initialised */
  /* idxs holds 0-based indices into x[] for the actors in this group. */
  int *idxs          = (int *)Calloc(n1, int);

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
  Free(seen);

  /* Not enough actors to form a k-subset. */
  if(ng < k){
    #if DEBUG_COV_DIFF
      Rprintf("[cov_diff][group_covdiff] g=%d ng=%d < k=%d -> 0\n",
              (int)g, ng, k);
    #endif
    Free(idxs);
    return 0.0;
  }

  /* Sum D(S) over all k-subsets S of the group. */
  double sumD = 0.0;
  int *comb = (int *)Calloc(k, int);
  sum_D_rec(0, 0, k, ng, idxs, x, comb, &sumD);
  Free(comb);

  double res = sumD;

  /* Optional per-group normalisation by C(n_g, k). */
  if(norm_flag){
    double denom = CHOOSE(ng, k);
    if(denom > 0.0){
      res /= denom;
    }else{
      res = 0.0;
    }
  }

  #if DEBUG_COV_DIFF
    Rprintf("[cov_diff][group_covdiff] g=%d ng=%d k=%d norm_flag=%d -> res=%g\n",
            (int)g, ng, k, norm_flag, res);
  #endif

  Free(idxs);
  return res;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_diff (one-toggle)                                    */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_diff`.
 *
 * @details
 *  This is the {ergm} change-statistic function registered as
 *  ::c_cov_diff via ::C_CHANGESTAT_FN. It computes the local change
 *  Δ in the cov_diff statistic for a single membership toggle between
 *  an actor and a group.
 *
 *  The layout of INPUT_PARAM is:
 *
 *    INPUT_PARAM[0]     = n1         (number of actors)
 *    INPUT_PARAM[1]     = k          (subset size)
 *    INPUT_PARAM[2]     = norm_flag  (0 raw, 1 per-group normalisation)
 *    INPUT_PARAM[3..]   = x[0..n1-1] (numeric covariate on actors)
 *
 *  For each toggle:
 *    1. Identify the actor vertex and the group vertex using the
 *       bipartite boundary between actor mode and group mode.
 *    2. Evaluate the group contribution before the virtual toggle
 *       via group_covdiff().
 *    3. Apply a virtual toggle (TOGGLE) on the actor–group edge.
 *    4. Evaluate the group contribution after the virtual toggle.
 *    5. Undo the virtual toggle.
 *    6. Update:
 *
 *         CHANGE_STAT[0] += F_after - F_before.
 *
 *  The parameter @p edgestate is unused here, because the function
 *  explicitly performs a virtual TOGGLE to obtain the "after" state.
 *
 * @param tail       Tail vertex of the toggled edge.
 * @param head       Head vertex of the toggled edge.
 * @param mtp        Pointer to the model term parameters (unused here,
 *                   but required by the macro signature).
 * @param nwp        Pointer to the current network-plus workspace, providing
 *                   adjacency and degree information.
 * @param edgestate  Current state of the edge (unused in this implementation).
 */
C_CHANGESTAT_FN(c_cov_diff){
  /* 1) Reset the output buffer for THIS toggle. */
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  /* 2) Read inputs from INPUT_PARAM. */
  const double *ip       = INPUT_PARAM;
  const int n1           = (int)ip[0];   /* number of actors (actor mode) */
  const int k            = (int)ip[1];   /* subset size for cov_diff */
  const int norm_flag    = (int)ip[2];   /* 0 raw, 1 per-group normalisation */
  const double *x        = ip + 3;       /* actor covariate values */

  #if DEBUG_COV_DIFF
    Rprintf("[cov_diff] n1=%d k=%d norm_flag=%d\n", n1, k, norm_flag);
  #endif

  /* 3) Identify actor and group vertices for the current toggle.
   *
   * Actor vertices are 1..n1 (actor mode).
   * Group vertices are > n1 (group mode).
   */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;  /* actor vertex in the actor mode */
  Vertex group = (a <= (Vertex)n1) ? b : a;  /* group vertex in the group mode */
  UNUSED_WARNING(actor);

  /* 4) Evaluate group contribution before and after a virtual toggle. */
  double F_before = group_covdiff(group, n1, k, norm_flag, x, nwp);

  /* Apply a virtual toggle (single-edge API). */
  TOGGLE(a, b);

  double F_after  = group_covdiff(group, n1, k, norm_flag, x, nwp);

  /* Undo the virtual toggle to restore the original state. */
  TOGGLE(a, b);

  /* 5) Update change statistic: Δ = F_after - F_before. */
  CHANGE_STAT[0] += (F_after - F_before);

  #if DEBUG_COV_DIFF
    Rprintf("[cov_diff] a=%d b=%d group=%d before=%g after=%g delta=%g\n",
            (int)a, (int)b, (int)group, F_before, F_after, (F_after - F_before));
  #endif
}
