/**
 * @file
 * @brief  Change statistic for the ERPM term `cov_ingroup` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `cov_ingroup`, which couples group size and the sum of a numeric actor
 *  covariate inside each group, with an optional filter on group sizes.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  The network is bipartite:
 *    - actor mode  = vertices representing actors (individuals),
 *    - group mode  = vertices representing structural groups.
 *
 *  Membership is encoded as edges between actors and groups. For a group
 *  vertex g in the group mode:
 *
 *    - A(g)   = set of actors adjacent to g (members of group g),
 *    - n_g    = |A(g)| = group size,
 *    - x_i    = numeric covariate value for actor i.
 *
 *  Let S be a (possibly empty) set of allowed group sizes. Define
 *  the ingroup statistic:
 *
 *      T(p; x)
 *        = ∑_g  1[n_g ∈ S] * n_g * ∑_{i ∈ A(g)} x_i,
 *
 *  where:
 *    - p encodes the partition via the bipartite network,
 *    - 1[n_g ∈ S] = 1 if n_g ∈ S, else 0; if S is empty, all sizes
 *      are accepted (no filter).
 *
 *  Intuitively, for each group g, the contribution is the group size
 *  multiplied by the total sum of x inside the group, optionally masked
 *  by a size filter.
 *
 *  ------------------------------------------------------------
 *  Local change for a membership toggle
 *  ------------------------------------------------------------
 *
 *  A single toggle flips the membership of one actor i in one group g:
 *    - addition  : actor i joins group g,
 *    - deletion  : actor i leaves group g.
 *
 *  Let:
 *    - n      = n_g       (group size before the toggle),
 *    - X      = ∑_{j ∈ A(g)} x_j   (sum of x in g before the toggle),
 *    - x_i    = covariate value of the toggled actor i,
 *    - n'     = n ± 1     (group size after the toggle),
 *    - X'     = X ± x_i   (sum after the toggle),
 *    - w(n)   = 1[n ∈ S]  (or 1 for all n if S is empty).
 *
 *  Then the contribution of group g before and after the toggle is:
 *
 *      T_before(g) = n  * X  * w(n),
 *      T_after(g)  = n' * X' * w(n').
 *
 *  The local change contributed by group g to the global statistic is:
 *
 *      Δ = T_after(g) - T_before(g)
 *        = (n' * X' * w(n')) - (n * X * w(n)).
 *
 *  All other groups are unchanged by this toggle and contribute zero
 *  to the change statistic.
 *
 *  ------------------------------------------------------------
 *  Bipartite structure and actors/groups
 *  ------------------------------------------------------------
 *
 *  At the C level, the bipartite structure is encoded as:
 *    - the first BIPARTITE vertices are actors (actor mode),
 *    - the remaining vertices are groups (group mode).
 *
 *  This change statistic assumes:
 *    - the number of actors n1 stored in INPUT_PARAM[0] matches BIPARTITE,
 *    - each membership toggle connects exactly one actor (vertex ≤ n1)
 *      and one group (vertex > n1).
 *
 *  The group vertex is detected as the endpoint in the group mode
 *  (vertex index > BIPARTITE), and the actor vertex as the endpoint
 *  in the actor mode (vertex index ≤ BIPARTITE).
 *
 *  ------------------------------------------------------------
 *  INPUT_PARAM layout (from InitErgmTerm.cov_ingroup)
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
 *    - n1          = number of actors (size of the actor mode),
 *    - L           = length of the size filter S,
 *    - sizes[ ]    = allowed group sizes (as doubles, cast to int),
 *    - x[1..n1]    = numeric covariate values on the actor mode.
 *
 *  In C, this becomes:
 *
 *    INPUT_PARAM[0]          = n1
 *    INPUT_PARAM[1]          = L
 *    INPUT_PARAM[2..2+L-1]   = sizes[0..L-1]
 *    INPUT_PARAM[2+L..]      = x[0..n1-1]
 *
 *  Special case:
 *    - If L == 0, every group size is accepted (no filter).
 *
 *  The term returns a single scalar statistic:
 *
 *    - N_CHANGE_STATS = 1,
 *    - CHANGE_STAT[0] is updated by the local Δ at each toggle.
 *
 *  ------------------------------------------------------------
 *  Implementation notes
 *  ------------------------------------------------------------
 *
 *  For each toggle:
 *    1. Identify the actor vertex v1 and the group vertex v2.
 *    2. Read the current degree of v2 as:
 *         deg_old = OUT_DEG[v2] + IN_DEG[v2],
 *       which equals n_g before the toggle.
 *    3. Compute X = ∑_{i ∈ A(g)} x_i by traversing neighbours of v2
 *       in the actor mode. In the current implementation:
 *         - only outgoing edges from v2 are traversed, which is
 *           sufficient for undirected bipartite representations where
 *           group–actor memberships appear as out-edges from groups.
 *    4. Build the "after" quantities:
 *         n_new = deg_old ± 1,
 *         X_new = X ± x_i,
 *       depending on whether the toggle is an addition or deletion.
 *    5. Check the size filter S on n and n_new, via in_sizes_set().
 *    6. Compute:
 *         Δ = (n_new * X_new * w_new) - (deg_old * X * w_old),
 *       and add it to CHANGE_STAT[0].
 *
 *  The code does not toggle the edge inside the network; it uses the
 *  current adjacency combined with the sign derived from edgestate.
 *
 *  ------------------------------------------------------------
 *  Complexity
 *  ------------------------------------------------------------
 *
 *  For a single toggle involving group g:
 *    - retrieving deg_old is O(1),
 *    - building X is O(deg(g)) via neighbour traversal,
 *    - all scalar operations are O(1).
 *
 *  No persistent state is stored across toggles; recomputation is
 *  local to the affected group.
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  The R initialiser (InitErgmTerm.cov_ingroup):
 *    - validates the actor attribute x and the size filter S,
 *    - packs n1, L, sizes, and x into INPUT_PARAM,
 *    - sets emptynwstats = 0 and a single coef.name.
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
 *  # Numeric covariate on actors (e.g., a score or weight)
 *  x <- c(1.0, 2.0, 0.5, 1.5, 3.0, 2.5)
 *
 *  # Unfiltered cov_ingroup: all group sizes contribute
 *  fit1 <- erpm(partition ~ cov_ingroup(attr = x))
 *  summary(fit1)
 *
 *  # Filtered version: only groups with size 2 or 3 contribute
 *  fit2 <- erpm(partition ~ cov_ingroup(attr = x, size = c(2, 3)))
 *  summary(fit2)
 *
 *  # Internally, each membership toggle between an actor and a group
 *  # calls c_cov_ingroup(), which:
 *  #   - reconstructs the sum of x in the affected group,
 *  #   - computes the before/after quantities n * X,
 *  #   - applies the size filter,
 *  #   - updates CHANGE_STAT[0] by the local Δ.
 *  @endcode
 */

#include <math.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"

/* -------------------------------------------------------------------------- */
/* Helper: membership size filter                                             */
/* -------------------------------------------------------------------------- */

/**
 * @brief Check whether a group size belongs to the allowed size set S.
 *
 * @details
 *  The size filter S is encoded as an array of doubles @p in of length @p L.
 *  Each entry is cast to int and compared to @p n.
 *
 *  Special case:
 *    - If L <= 0, the filter is considered inactive and all sizes are
 *      accepted (the function returns 1).
 *
 * @param n   Group size to be tested.
 * @param in  Pointer to the array of allowed sizes (stored as doubles).
 * @param L   Number of allowed sizes in @p in.
 *
 * @return 1 if @p n is in S, or if L <= 0; 0 otherwise.
 */
static inline int in_sizes_set(int n, const double *in, int L){
  if(L<=0) return 1; /* no filter: all sizes accepted */
  /* in[0..L-1] stores allowed sizes as doubles, cast to int for comparison */
  for(int k=0; k<L; ++k){
    if((int)in[k] == n) return 1;
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_ingroup (one-toggle)                                 */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_ingroup`.
 *
 * @details
 *  This is the {ergm} change-statistic function registered as
 *  ::c_cov_ingroup via ::C_CHANGESTAT_FN. It computes the local change
 *  Δ for the ingroup covariate statistic when a single membership edge
 *  between an actor and a group is toggled.
 *
 *  The layout of INPUT_PARAM is:
 *
 *    INPUT_PARAM[0]        = n1               (number of actors, actor mode)
 *    INPUT_PARAM[1]        = L                (size of the size filter S)
 *    INPUT_PARAM[2..2+L-1] = sizes[0..L-1]    (allowed group sizes, as doubles)
 *    INPUT_PARAM[2+L..]    = x[0..n1-1]       (numeric covariate on actors)
 *
 *  The function:
 *    - identifies the actor vertex v1 and the group vertex v2 for the
 *      current toggle,
 *    - reconstructs the group size and sum of x in the group,
 *    - infers n_new and X_new after the hypothetical toggle,
 *    - applies the size filter before and after,
 *    - computes:
 *
 *          Δ = (n_new * X_new * w_new) - (deg_old * X * w_old),
 *
 *      and adds Δ to CHANGE_STAT[0].
 *
 *  The parameter @p edgestate indicates whether the edge currently exists:
 *    - edgestate = 0 → edge is absent (toggle = addition),
 *    - edgestate = 1 → edge is present (toggle = deletion).
 *
 *  The network itself is not modified; only the sign and magnitude of
 *  the change are computed based on the current state.
 *
 * @param tail       Tail vertex of the toggled edge.
 * @param head       Head vertex of the toggled edge.
 * @param mtp        Pointer to the model term parameters (unused here,
 *                   but required by the macro signature).
 * @param nwp        Pointer to the network-plus workspace, providing
 *                   adjacency and degree information.
 * @param edgestate  Current state of the edge:
 *                   - 0 if the edge is absent (toggle = addition),
 *                   - 1 if the edge is present (toggle = deletion).
 */
C_CHANGESTAT_FN(c_cov_ingroup){
  /* 1) Reset the output buffer for THIS toggle.
   *
   * {ergm} accumulates contributions across multiple toggles; this function
   * sets the local Δ for the current one.
   */
  ZERO_ALL_CHANGESTATS(0);

  /* 2) Read inputs from INPUT_PARAM. */
  const int n1 = (int)INPUT_PARAM[0];  /* number of actors in actor mode */
  const int L  = (int)INPUT_PARAM[1];  /* length of size filter S */

  /* sizes[] points to the allowed size set S, if any; x[] to actor covariates. */
  const double *sizes = (L>0) ? (&INPUT_PARAM[2])      : NULL;
  const double *x     = (L>0) ? (&INPUT_PARAM[2+L])    : (&INPUT_PARAM[2]);

  /* 3) Identify the actor vertex v1 and the group vertex v2 for this toggle.
   *
   * The boundary between actor mode and group mode is BIPARTITE, which
   * must coincide with n1.
   *
   *  - actor vertices: 1 .. n1,
   *  - group vertices: n1+1 .. N_NODES.
   */
  const Vertex n1_lim = (Vertex)BIPARTITE; /* should coincide with n1 */
  Vertex v2 = (tail > n1_lim) ? tail : head; /* group vertex in the group mode */
  Vertex v1 = (tail > n1_lim) ? head : tail; /* actor vertex in the actor mode */

  /* Safety check: v1 must be a valid actor index in [1..n1]. */
  if(v1 < (Vertex)1 || v1 > (Vertex)n1){
    CHANGE_STAT[0] += 0.0;
    return;
  }

  /* 4) Compute the group size and the sum of x in the group BEFORE the toggle.
   *
   *  - deg_old: degree of v2, equals n_g before the toggle.
   *  - X:       sum of x_i over all actors connected to v2 in the actor mode.
   */
  int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);

  double X = 0.0;
  Vertex h;
  Edge e;

  /* Traverse outgoing edges from the group vertex v2.
   *
   * For undirected bipartite representations where memberships appear as
   * edges from groups to actors, iterating over OUTEDGES is sufficient.
   */
  STEP_THROUGH_OUTEDGES(v2, e, h){
    if(h >= (Vertex)1 && h <= (Vertex)n1){
      /* x is stored as x[0..n1-1], vertex indices are 1-based. */
      X += x[(int)h - 1];
    }
  }
  /* If memberships could also appear as incoming edges to v2, one could
   * add a symmetric pass over INEDGES here (currently not needed in the
   * assumed representation).
   *
   * Example (commented out):
   *
   * STEP_THROUGH_INEDGES(v2, e, h){
   *   if(h >= (Vertex)1 && h <= (Vertex)n1) X += x[(int)h - 1];
   * }
   */

  /* 5) Compute "after-toggle" quantities.
   *
   * edgestate = 0 → addition:  n_new = n + 1, X_new = X + x_i
   * edgestate = 1 → deletion:  n_new = n - 1, X_new = X - x_i
   */
  const int is_add   = (edgestate==0);          /* 1 if addition, 0 if deletion */
  const int n_new    = deg_old + (is_add ? +1 : -1);
  const double xi    = x[(int)v1 - 1];
  const double X_new = X + (is_add ? +xi : -xi);

  /* 6) Evaluate the size filter S before and after the toggle. */
  const int w_old = in_sizes_set(deg_old, sizes, L);  /* 1 if deg_old ∈ S */
  const int w_new = in_sizes_set(n_new,  sizes, L);   /* 1 if n_new  ∈ S */

  /* 7) Compute the local change Δ and add it to CHANGE_STAT[0].
   *
   *    Δ = (n_new * X_new * w_new) - (deg_old * X * w_old).
   */
  double d = 0.0;
  d = (w_new ? ( (double)n_new * X_new ) : 0.0)
    - (w_old ? ( (double)deg_old * X   ) : 0.0);

  CHANGE_STAT[0] += d;
}
