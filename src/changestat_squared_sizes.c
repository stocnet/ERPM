/**
 * @file changestat_squared_sizes.c
 * @brief Change statistic for the ERPM term `squared_sizes` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `squared_sizes(sizes, pow)` on a bipartite network encoded as:
 *    - actor mode  = actor vertices,
 *    - group mode  = group vertices.
 *
 *  Each group-mode vertex represents a structural group. Its size is defined
 *  as the total degree of the vertex (sum of incoming and outgoing degrees)
 *  induced by membership edges from the actor mode.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (group sizes in a bipartite encoding)
 *  ------------------------------------------------------------
 *
 *  For each group vertex g, let:
 *    - deg_g be its total degree (its group size, i.e. number of members),
 *    - S be the set of admissible sizes (argument `sizes`),
 *    - pow be the scalar exponent (argument `pow`).
 *
 *  The contribution of group g to the statistic is:
 *
 *      f(deg_g) =
 *        { deg_g^{pow}  if deg_g ∈ S
 *        { 0            otherwise.
 *
 *  The full statistic is:
 *
 *      Stat = sum_over_groups f(deg_g).
 *
 *  A single membership toggle affects exactly one group g and changes its
 *  size from deg_old to deg_new = deg_old ± 1. The local change in the
 *  statistic is therefore:
 *
 *      Δ = f(deg_new) − f(deg_old).
 *
 *  The function implemented here computes this Δ and accumulates it
 *  in CHANGE_STAT[0].
 *
 *  ------------------------------------------------------------
 *  Implementation in {ergm} (one-toggle)
 *  ------------------------------------------------------------
 *
 *  - A bipartite network is assumed:
 *      - the actor mode occupies the first BIPARTITE vertices,
 *      - the group mode occupies the remaining vertices.
 *  - Each membership toggle connects exactly one actor-mode vertex and one
 *    group-mode vertex.
 *  - This change-statistic function:
 *      1. Identifies the group-mode vertex affected by the toggle.
 *      2. Reads its current size deg_old from OUT_DEG + IN_DEG.
 *      3. Derives deg_new = deg_old ± 1 from the @p edgestate flag.
 *      4. Reads the parameters:
 *           - pow     = INPUT_PARAM[0]        (scalar exponent),
 *           - K       = INPUT_PARAM[1]        (number of admissible sizes),
 *           - sizes_k = INPUT_PARAM[2 + k]    (k = 0..K-1).
 *      5. Computes:
 *
 *           Δ = f(deg_new) − f(deg_old),
 *
 *         where
 *
 *           f(d) = d^pow if d == sizes_k for some k, and 0 otherwise.
 *
 *         This Δ is then added to CHANGE_STAT[0].
 *
 *  - The macro ::C_CHANGESTAT_FN declares the function with the signature
 *    required by {ergm} and exposes:
 *      - N_CHANGE_STATS  (number of statistics, expected to be 1 here),
 *      - INPUT_PARAM     (packed parameters),
 *      - CHANGE_STAT     (output buffer).
 *
 *  Parameter packing:
 *  - The R-side initialiser (InitErgmTerm.squared_sizes) constructs:
 *
 *        INPUT_PARAM = c(pow, K, sizes_1, ..., sizes_K)
 *
 *    so that:
 *
 *        INPUT_PARAM[0] = pow        (integer >= 1)
 *        INPUT_PARAM[1] = K          (number of admissible sizes)
 *        INPUT_PARAM[2 + k] = sizes_{k+1} (k = 0..K-1).
 *
 *  ------------------------------------------------------------
 *  Complexity
 *  ------------------------------------------------------------
 *
 *  - O(K) par toggle:
 *      - un lookup de degré pour le groupe touché,
 *      - au plus deux scans du vecteur des tailles admissibles (deg_old / deg_new),
 *        avec au plus deux exponentiations (pour deg_old, deg_new).
 *  - Pas de parcours d’adjacence, pas de scan sur les autres groupes.
 *
 *  ------------------------------------------------------------
 *  R interface and usage
 *  ------------------------------------------------------------
 *
 *  On the R side, the corresponding ERPM term can be accessed via the wrapper.
 *
 *  @example Usage (R)
 *  @code{.r}
 *  library(ERPM)
 *
 *  # Example: partition of 6 actors into 3 groups
 *  part <- c(1, 1, 2, 2, 3, 3)
 *
 *  # Model with a squared_sizes term on group sizes equal to 2 or 3
 *  fit <- erpm(
 *    partition ~ squared_sizes(sizes = c(2L, 3L), pow = 2L)
 *  )
 *  summary(fit)
 *
 *  # Interpretation of one toggle:
 *  # Suppose a sampler proposes to add an actor to a group of current size 2:
 *  #   deg_old = 2
 *  #   deg_new = 3
 *  #
 *  # For sizes ∈ {2,3}, pow = 2:
 *  #   f(deg_old) = 2^2 = 4
 *  #   f(deg_new) = 3^2 = 9
 *  #
 *  #   Δ = f(deg_new) - f(deg_old) = 9 - 4 = 5
 *  #
 *  # This is exactly what c_squared_sizes adds to CHANGE_STAT[0].
 *  @endcode
 *
 *  @test
 *  A self-test can:
 *    - build small bipartite networks from known partitions,
 *    - compute the reference statistic by summing d^pow over group sizes
 *      belonging to the target set of sizes,
 *    - call `summary()` on an ERGM/ERPM model with `squared_sizes`,
 *    - compare the reported statistic to the reference value,
 *    - apply a series of single-edge toggles and verify that the incremental
 *      changes match Δ = f(deg_new) − f(deg_old).
 */

#include <math.h>
#include <R_ext/Print.h>        // Rprintf
#include "ergm_changestat.h"
#include "ergm_storage.h"

/**
 * @def DEBUG_SQUARED_SIZES
 * @brief Enable or disable verbose debugging for ::c_squared_sizes.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * for each toggle:
 *  - endpoints of the toggled edge,
 *  - index and degrees of the affected group vertex,
 *  - parameters (pow, K, sizes_k) and Δ value.
 *
 * When set to 0, no debug traces are emitted.
 */
#define DEBUG_SQUARED_SIZES 0

/* -------------------------------------------------------------------------- */
/* Utility: fast integer exponentiation                                       */
/* -------------------------------------------------------------------------- */

/**
 * @brief Fast exponentiation for integer base and non-negative exponent.
 *
 * @details
 *  Computes base^exp using exponentiation by squaring.
 *  This runs in O(log(exp)), which is typically faster and more
 *  numerically stable than a naive loop in O(exp) for larger exponents.
 *
 *  For convenience, if @p exp <= 0 the function returns 1.0, so in
 *  particular base^0 = 1.0 for any base.
 *
 * @param base Integer base.
 * @param exp  Non-negative integer exponent.
 * @return base^exp as a double if exp > 0, or 1.0 if exp <= 0.
 */
static inline double ipow_int(int base, int exp){
  if(exp <= 0) return 1.0;
  double r = 1.0, b = (double)base;
  while(exp){
    if(exp & 1) r *= b;  // when the current bit is set, multiply the accumulator
    b *= b;              // square the base for the next bit
    exp >>= 1;           // shift to the next bit
  }
  return r;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: squared_sizes                                            */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `squared_sizes(sizes, pow)`.
 *
 * @details
 *  This function is registered via ::C_CHANGESTAT_FN as ::c_squared_sizes.
 *  It processes a single membership toggle between an actor-mode vertex and
 *  a group-mode vertex in a bipartite network and computes the local change
 *  in the aggregated `squared_sizes` statistic:
 *
 *      T(y) = sum_{g in G} 1[deg(g) in sizes] * deg(g)^pow.
 *
 *  Actor / group modes:
 *  - The actor mode is represented by the first BIPARTITE vertices.
 *  - The group mode is represented by all remaining vertices.
 *
 *  Parameter packing (R side):
 *  - The R-side initializer constructs
 *
 *        INPUT_PARAM = c(pow, K, sizes_1, ..., sizes_K)
 *
 *    where:
 *      - pow       = common exponent (integer >= 1),
 *      - K         = number of admissible group sizes,
 *      - sizes_i   = i-th admissible group size (integer >= 1).
 *
 *  For each toggle:
 *  - The function:
 *      1. Zeros the CHANGE_STAT buffer (single component).
 *      2. Identifies the group-mode vertex affected by the toggle.
 *      3. Reads the group size before the toggle:
 *           deg_old = OUT_DEG[v2] + IN_DEG[v2].
 *      4. Determines if the toggle is an addition or a deletion via @p edgestate
 *         and sets:
 *           deg_new = deg_old + 1  (addition),
 *           deg_new = deg_old - 1  (deletion).
 *      5. Checks whether deg_old and deg_new belong to the set `sizes`.
 *      6. Computes:
 *
 *           Δ = 1[deg_new in sizes] * deg_new^pow
 *             − 1[deg_old in sizes] * deg_old^pow,
 *
 *         and adds Δ to CHANGE_STAT[0].
 */
C_CHANGESTAT_FN(c_squared_sizes){

  /* 1) Always reset the output buffer for THIS call.
   *
   * {ergm} sums the vectors returned by successive calls. Here we only
   * compute the local contribution of the current toggle.
   */
  ZERO_ALL_CHANGESTATS();

  /* 2) Number of vertices in the actor mode.
   *
   * In a bipartite network, BIPARTITE is the number of actor-mode vertices.
   * All vertices with index > BIPARTITE belong to the group mode.
   */
  const int n1 = BIPARTITE;

  /* 3) Read the endpoints of the current toggle. */
  Vertex t = tail;
  Vertex h = head;

  /* 4) Identify the group-mode vertex.
   *
   * The toggle involves exactly one actor-mode vertex and one group-mode vertex.
   * The group vertex is the one whose index is strictly greater than n1.
   */
  Vertex v2 = (t > n1) ? t : h;

  /* 5) Group size before and after the toggle (local computation).
   *
   * OUT_DEG and IN_DEG are internal degree arrays obtained from nwp.
   * The group size is the sum of outgoing and incoming degrees of v2.
   * The new size is obtained by incrementing or decrementing by 1, depending
   * on whether the edge is being added or removed.
   */
  int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);
  int delta   = edgestate ? -1 : +1;   // present -> deletion (-1), absent -> addition (+1)
  int deg_new = deg_old + delta;

  /* Defensive note: in a consistent chain, deg_new should stay non-negative.
   * Any negative value would signal an inconsistency between edgestate and
   * the internal network representation.
   */

  /* 6) Aggregated statistic over a set of admissible sizes.
   *
   * INPUT_PARAM layout (length = K + 2):
   *   INPUT_PARAM[0]         = pow        (exponent, integer >= 1)
   *   INPUT_PARAM[1]         = K          (number of admissible sizes)
   *   INPUT_PARAM[2..(K+1)]  = sizes_1..K (admissible group sizes)
   *
   * There is a single ERGM statistic component (N_CHANGE_STATS == 1),
   * which aggregates the contributions of all admissible sizes.
   */
  const int power = (int)INPUT_PARAM[0];
  const int K     = (int)INPUT_PARAM[1];

  int match_old = 0;
  int match_new = 0;

  for(int i = 0; i < K; ++i){
    int size_i = (int)INPUT_PARAM[2 + i];
    if(deg_old == size_i) match_old = 1;
    if(deg_new == size_i) match_new = 1;
  }

  double d = 0.0;

  if(match_new)
    d += ipow_int(deg_new, power);

  if(match_old)
    d -= ipow_int(deg_old, power);

  CHANGE_STAT[0] += d;

#if DEBUG_SQUARED_SIZES
  Rprintf("[C:c_squared_sizes] tail=%d head=%d | group=%d | edgestate=%d | "
          "deg_old=%d -> deg_new=%d | K=%d pow=%d | Δ=%.2f | cumul=%.2f\n",
          (int)t, (int)h, (int)v2, (int)edgestate,
          deg_old, deg_new, K, power, d, CHANGE_STAT[0]);
#endif
}