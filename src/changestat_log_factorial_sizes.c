/**
 * @file
 * @brief Change statistic for the ERPM term `log_factorial_sizes` (one-toggle, non-vectorised).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `log_factorial_sizes`, defined on a bipartite network with:
 *    - actor mode  = actor vertices,
 *    - group mode  = group vertices.
 *
 *  Each group-mode vertex represents a structural group. Its size is defined
 *  as the total degree of the vertex induced by membership edges from the
 *  actor mode.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (log-factorial of group sizes)
 *  ------------------------------------------------------------
 *
 *  Let:
 *    - G be the set of group-mode vertices,
 *    - deg(g) be the size (degree) of group g.
 *
 *  The statistic is:
 *
 *      Stat = sum_{g ∈ G} f(deg(g)),
 *      f(n) = log((n − 1)!) = lgamma(n),  with the convention f(0) = 0.
 *
 *  Using the identity lgamma(n+1) − lgamma(n) = log(n) for n ≥ 1, we obtain
 *  simple local updates when a single membership edge is toggled.
 *
 *  Local variations when a single toggle modifies one group g:
 *
 *    - Edge addition:
 *        deg_old = n  →  deg_new = n + 1
 *
 *        Δ = f(n + 1) − f(n)
 *          = lgamma(n + 1) − lgamma(n)
 *          = log(n)  if n ≥ 1
 *          = 0       if n = 0  (since f(0) = 0 and f(1) = 0)
 *
 *    - Edge deletion:
 *        deg_old = n  →  deg_new = n − 1
 *
 *        Δ = f(n − 1) − f(n)
 *          = lgamma(n − 1) − lgamma(n)
 *          = −log(n − 1)  if n ≥ 2
 *          = 0            if n = 1  (since f(1) = 0 and f(0) = 0)
 *
 *  The implementation uses these closed-form increments and never calls
 *  lgamma() inside the change statistic.
 *
 *  ------------------------------------------------------------
 *  Implementation in {ergm} (one-toggle)
 *  ------------------------------------------------------------
 *
 *  - A bipartite network is assumed, with:
 *      - actor mode  = the first BIPARTITE vertices,
 *      - group mode  = the remaining vertices.
 *
 *  - Each membership toggle connects exactly one actor-mode vertex and one
 *    group-mode vertex; only that group is affected.
 *
 *  - The macro ::C_CHANGESTAT_FN declares the function with the signature
 *    required by {ergm} and exposes:
 *      - N_CHANGE_STATS (here equal to 1),
 *      - CHANGE_STAT    (output buffer),
 *      - INPUT_PARAM    (unused here, non-vectorised term).
 *
 *  - For each call:
 *      1. The output buffer CHANGE_STAT is reset with ::ZERO_ALL_CHANGESTATS(0).
 *      2. The group-mode vertex v2 impacted by the toggle is identified.
 *      3. Its current size is read as:
 *           deg_old = OUT_DEG[v2] + IN_DEG[v2].
 *      4. The new size is obtained locally via:
 *           deg_new = deg_old + 1  for an addition,
 *           deg_new = deg_old − 1  for a deletion.
 *      5. The local increment Δ is computed using the closed-form formulas
 *         above, without evaluating lgamma().
 *      6. The scalar statistic is updated as:
 *           CHANGE_STAT[0] += Δ.
 *
 *  ------------------------------------------------------------
 *  Complexity
 *  ------------------------------------------------------------
 *
 *  - O(1) per toggle:
 *      - a constant number of degree reads,
 *      - a constant amount of arithmetic and at most one log() call.
 *  - No neighbour traversal, no scanning of other groups.
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  - The R-side initialiser (InitErgmTerm.log_factorial_sizes) sets:
 *      - N_CHANGE_STATS = 1 (non-vectorised term),
 *      - no INPUT_PARAM (the term has no hyper-parameters at the C level),
 *      - emptynwstats = 0.
 *
 *  ------------------------------------------------------------
 *  @example Usage (R)
 *  ------------------------------------------------------------
 *  @code{.r}
 *  library(ERPM)
 *
 *  # Example: partition of 5 actors into 2 groups
 *  part <- c(1, 1, 1, 2, 2)   # group sizes: 3 and 2
 *
 *  # Fit an ERPM with the log_factorial_sizes term
 *  fit <- erpm(partition ~ log_factorial_sizes)
 *  summary(fit)
 *
 *  # Interpretation of one toggle:
 *  # Suppose the sampler proposes to add an actor to a group of current size n = 2:
 *  #   deg_old = 2
 *  #   deg_new = 3
 *  #
 *  # The local increment is:
 *  #   Δ = lgamma(3) - lgamma(2)
 *  #     = log(2)
 *  #
 *  # This is exactly what c_log_factorial_sizes adds to CHANGE_STAT[1].
 *  @endcode
 *
 *  @test
 *  A self-test can:
 *    - build several small bipartite networks from known partitions,
 *    - compute the reference statistic as sum_g lgamma(deg(g)) with f(0) = 0,
 *    - call `summary()` on a model containing `log_factorial_sizes`,
 *    - verify equality between the reported statistic and the reference,
 *    - apply single-edge toggles and check that the observed Δ matches
 *      the closed-form formulas for additions and deletions.
 */

#include <math.h>
#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"

/**
 * @def DEBUG_LOG_FACTORIAL
 * @brief Enable verbose debugging output for ::c_log_factorial_sizes.
 *
 * Set this macro to 1 to print diagnostic traces to the R console during
 * `summary()` or MCMC runs:
 *  - index of the affected group vertex,
 *  - degree before and after the toggle,
 *  - edge state and resulting increment Δ.
 *
 * When set to 0, no debug output is produced.
 */
#define DEBUG_LOG_FACTORIAL 0

/* -------------------------------------------------------------------------- */
/* Change statistic: log_factorial_sizes                                      */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `log_factorial_sizes`.
 *
 * @details
 *  This {ergm} change-statistic function is registered as ::c_log_factorial_sizes
 *  via ::C_CHANGESTAT_FN. It implements the one-toggle update for the statistic
 *  defined by:
 *
 *      Stat = sum_{groups g} lgamma(deg(g)),  with f(0) = 0.
 *
 *  The network is assumed bipartite, with an actor mode and a group mode.
 *  Each toggle connects exactly one actor-mode vertex and one group-mode vertex.
 *
 *  For each call:
 *    - The function:
 *        1. Resets CHANGE_STAT to zero for this toggle.
 *        2. Identifies the group-mode vertex v2 affected by the membership toggle.
 *        3. Reads the current group size:
 *             deg_old = OUT_DEG[v2] + IN_DEG[v2].
 *        4. Determines whether the toggle is an addition or deletion from
 *           @p edgestate and sets:
 *             deg_new = deg_old + 1  if the edge is absent (addition),
 *             deg_new = deg_old − 1  if the edge is present (deletion).
 *        5. Computes the local change Δ using only log() and the closed-form:
 *
 *             addition (deg_old = n):
 *               Δ =  log(n)    if n ≥ 1
 *               Δ =  0         if n = 0
 *
 *             deletion (deg_old = n):
 *               Δ = −log(n − 1)  if n ≥ 2
 *               Δ =  0           if n = 1
 *
 *        6. Accumulates the change:
 *             CHANGE_STAT[0] += Δ.
 *
 *  The global statistic over all toggles is obtained by the {ergm} engine
 *  via accumulation of the per-toggle contributions.
 *
 * @param tail       Tail vertex of the toggled edge (actor or group).
 * @param head       Head vertex of the toggled edge (actor or group).
 * @param mtp        Pointer to the model term structure (unused directly here,
 *                   but required by the macro signature).
 * @param nwp        Pointer to the network-plus workspace (provides OUT_DEG,
 *                   IN_DEG, and other internals).
 * @param edgestate  Current state of the edge:
 *                   - 0 if the edge is absent (toggle = addition),
 *                   - 1 if the edge is present (toggle = deletion).
 *
 * @note
 *  - The actor mode occupies the first BIPARTITE vertices; the group mode
 *    consists of the remaining vertices.
 *  - The function assumes that each membership toggle involves exactly one
 *    actor-mode vertex and one group-mode vertex; only that group’s size
 *    changes.
 *  - The term is non-vectorised: N_CHANGE_STATS is always 1 and INPUT_PARAM
 *    is unused.
 */
C_CHANGESTAT_FN(c_log_factorial_sizes){
  /* 1) Reset the output buffer for THIS toggle.
   *
   * {ergm} accumulates the contributions from multiple calls; here we only
   * report the local increment Δ for the current toggle.
   */
  ZERO_ALL_CHANGESTATS(0);

  /* 2) Number of vertices in the actor mode.
   *
   * In a bipartite encoding, vertices 1..BIPARTITE belong to the actor mode,
   * while vertices with index > BIPARTITE belong to the group mode.
   */
  const int n1 = BIPARTITE;

  /* 3) Identify the affected group-mode vertex.
   *
   * The toggled edge always connects one actor-mode vertex and one group-mode
   * vertex. The group vertex is the endpoint whose index is strictly greater
   * than n1.
   */
  Vertex v2 = (tail > n1) ? tail : head;

  /* 4) Group sizes before and after the toggle (local computation).
   *
   * OUT_DEG and IN_DEG are the internal degree arrays. Their sum gives the
   * group size for vertex v2. The new size is obtained by adding or removing
   * one membership depending on the current edge state.
   */
  int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);
  int delta   = edgestate ? -1 : +1;   // edge present -> deletion (-1), edge absent -> addition (+1)
  int deg_new = deg_old + delta;

  #if DEBUG_LOG_FACTORIAL
    Rprintf("[c_log_factorial_sizes] v2=%d | edgestate=%d | deg_old=%d -> deg_new=%d\n",
            (int)v2, (int)edgestate, deg_old, deg_new);
  #endif

  /* 5) Local increment Δ using closed-form formulas.
   *
   * We avoid calling lgamma() and use:
   *   addition:  Δ = log(deg_old)     if deg_old >= 1, else 0
   *   deletion:  Δ = -log(deg_old-1)  if deg_old >= 2, else 0
   *
   * This is consistent with the definition Stat = sum_g lgamma(deg(g))
   * under the convention f(0) = 0.
   */
  double d = 0.0;
  if(edgestate == 0){
    /* Addition: edge is currently absent, we add it.
     *
     * deg_old = n, deg_new = n+1:
     *   Δ = lgamma(n+1) - lgamma(n)
     *     = log(n) for n >= 1,
     *     = 0      for n = 0.
     */
    if(deg_old >= 1) d = log((double)deg_old);
  }else{
    /* Deletion: edge is currently present, we remove it.
     *
     * deg_old = n, deg_new = n-1:
     *   Δ = lgamma(n-1) - lgamma(n)
     *     = -log(n-1) for n >= 2,
     *     = 0         for n = 1.
     */
    if(deg_old >= 2) d = -log((double)(deg_old - 1));
  }

  /* 6) Accumulate the contribution (single scalar statistic).
   *
   * The term is non-vectorised, so N_CHANGE_STATS == 1 and we always update
   * CHANGE_STAT[0].
   */
  CHANGE_STAT[0] += d;

  #if DEBUG_LOG_FACTORIAL
    Rprintf("  Δ=%.9g\n", d);
  #endif
}
