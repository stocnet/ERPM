/**
 * @file
 * @brief Change statistic for the ERPM term `squared_sizes` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `squared_sizes(from, to, pow)` on a bipartite network encoded as:
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
 *    - [from_j, to_j) and pow_j be the parameters of sub-term j.
 *
 *  The contribution of group g to sub-term j is:
 *
 *      f_j(deg_g) =
 *        { deg_g^{pow_j}  if deg_g ∈ [from_j, to_j)
 *        { 0              otherwise.
 *
 *  The full statistic for sub-term j is:
 *
 *      Stat_j = sum_over_groups f_j(deg_g).
 *
 *  A single membership toggle affects exactly one group g and changes its
 *  size from deg_old to deg_new = deg_old ± 1. The local change for sub-term j
 *  is therefore:
 *
 *      Δ_j = f_j(deg_new) − f_j(deg_old).
 *
 *  The function implemented here computes this Δ_j for all j and accumulates
 *  it in CHANGE_STAT[j].
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
 *      4. For each sub-term j with parameters (from_j, to_j, pow_j), computes:
 *
 *            Δ_j = [deg_new ∈ [from_j, to_j)] * deg_new^{pow_j}
 *                − [deg_old ∈ [from_j, to_j)] * deg_old^{pow_j},
 *
 *         and adds Δ_j to CHANGE_STAT[j].
 *
 *  - The macro ::C_CHANGESTAT_FN declares the function with the signature
 *    required by {ergm} and exposes:
 *      - N_CHANGE_STATS  (number of sub-terms),
 *      - INPUT_PARAM     (packed parameters),
 *      - CHANGE_STAT     (output buffer).
 *
 *  Parameter packing:
 *  - The R-side initialiser (InitErgmTerm.squared_sizes) constructs:
 *
 *        INPUT_PARAM = c(rbind(from, to, pow))
 *
 *    so that for sub-term j:
 *
 *        INPUT_PARAM[3*j + 0] = from_j  (inclusive lower bound)
 *        INPUT_PARAM[3*j + 1] = to_j    (exclusive upper bound)
 *        INPUT_PARAM[3*j + 2] = pow_j   (integer >= 1).
 *
 *  ------------------------------------------------------------
 *  Complexity
 *  ------------------------------------------------------------
 *
 *  - O(N_CHANGE_STATS) per toggle:
 *      - one degree lookup for the affected group,
 *      - two range tests and up to two exponentiations (for deg_old, deg_new)
 *        per sub-term.
 *  - No adjacency traversal, no scanning of other groups.
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
 *  # Model with a squared_sizes term on group sizes in [2, Inf)
 *  # (in practice, 'to' is finite and handled by the initialiser)
 *  fit <- erpm(
 *    partition ~ squared_sizes(from = 2L, to = Inf, pow = 2L)
 *  )
 *  summary(fit)
 *
 *  # Interpretation of one toggle:
 *  # Suppose a sampler proposes to add an actor to a group of current size 2:
 *  #   deg_old = 2
 *  #   deg_new = 3
 *  #
 *  # For from = 2, to = 4, pow = 2:
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
 *      in [from, to),
 *    - call `summary()` on an ERGM/ERPM model with `squared_sizes`,
 *    - compare the reported statistic to the reference value,
 *    - apply a series of single-edge toggles and verify that the incremental
 *      changes match Δ_j = f_j(deg_new) − f_j(deg_old) for each sub-term j.
 */

#include <math.h>
#include <R_ext/Print.h>        // Rprintf
#include "ergm_changestat.h"
#include "ergm_storage.h"

/**
 * @def IN_RANGE
 * @brief Test if a value lies in a half-open interval [a, b).
 *
 * @param x Value to test.
 * @param a Inclusive lower bound.
 * @param b Exclusive upper bound.
 */
#define IN_RANGE(x,a,b) ((x) >= (a) && (x) < (b))

/**
 * @def DEBUG_SQUARED_SIZES
 * @brief Enable or disable verbose debugging for ::c_squared_sizes.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * for each toggle:
 *  - endpoints of the toggled edge,
 *  - index and degrees of the affected group vertex,
 *  - parameters (from, to, pow) and per-sub-statistic Δ and cumulative values.
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
 *  Computes \f$\text{base}^{\text{exp}}\f$ using exponentiation by squaring.
 *  This runs in \f$O(\log(\text{exp}))\f$, which is typically faster and more
 *  numerically stable than a naive loop in \f$O(\text{exp})\f$ for larger
 *  exponents.
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
 * @brief Change statistic for the ERPM term `squared_sizes(from, to, pow)`.
 *
 * @details
 *  This function is registered via ::C_CHANGESTAT_FN as ::c_squared_sizes.
 *  It processes a single membership toggle between an actor-mode vertex and
 *  a group-mode vertex in a bipartite network and computes the local change
 *  in the `squared_sizes` statistic for each sub-term.
 *
 *  Actor / group modes:
 *  - The actor mode is represented by the first BIPARTITE vertices.
 *  - The group mode is represented by all remaining vertices.
 *
 *  For each call:
 *  - The function:
 *      1. Zeros the CHANGE_STAT buffer.
 *      2. Identifies the group-mode vertex @c v2 affected by the toggle.
 *      3. Reads the group size before the toggle:
 *           deg_old = OUT_DEG[v2] + IN_DEG[v2].
 *      4. Determines if the toggle is an addition or a deletion via @p edgestate
 *         and sets:
 *           deg_new = deg_old + 1  (addition),
 *           deg_new = deg_old - 1  (deletion).
 *      5. For each sub-term j, it reads:
 *           from_j = INPUT_PARAM[3*j + 0],
 *           to_j   = INPUT_PARAM[3*j + 1],
 *           pow_j  = INPUT_PARAM[3*j + 2],
 *         and computes:
 *
 *           Δ_j = [deg_new ∈ [from_j, to_j)] * deg_new^{pow_j}
 *               − [deg_old ∈ [from_j, to_j)] * deg_old^{pow_j}.
 *
 *         This Δ_j is then added to CHANGE_STAT[j].
 *
 *  The global statistic is obtained by {ergm} by accumulating the contributions
 *  from all calls to this function over all toggles.
 *
 * @param tail       Tail vertex of the toggled edge (actor or group).
 * @param head       Head vertex of the toggled edge (actor or group).
 * @param mtp        Pointer to the model term structure (contains INPUT_PARAM,
 *                   number of statistics, etc.).
 * @param nwp        Pointer to the network-plus workspace (provides degrees,
 *                   adjacency, and other internal structures).
 * @param edgestate  Current state of the edge:
 *                   - 0 if the edge is absent (toggle = addition),
 *                   - 1 if the edge is present (toggle = deletion).
 *
 * @note
 *  - This function assumes a bipartite representation where each membership
 *    connects one actor-mode vertex and one group-mode vertex.
 *  - The group size is inferred from the degree of the group vertex without
 *    explicitly toggling the edge inside this function (un-toggle form).
 */
C_CHANGESTAT_FN(c_squared_sizes){

  /* 1) Always reset the output buffer for THIS call.
   *
   * {ergm} sums the vectors returned by successive calls. Here we only
   * compute the local contribution of the current toggle.
   */
  ZERO_ALL_CHANGESTATS(0);  // memset(CHANGE_STAT, 0, N_CHANGE_STATS * sizeof(double))

  /* 2) Number of vertices in the actor mode.
   *
   * In a bipartite network, BIPARTITE is the number of actor-mode vertices.
   * All vertices with index > BIPARTITE belong to the group mode.
   */
  const int n1 = BIPARTITE;

  /* 3) Read the endpoints of the current toggle. */
  Vertex t = tail, h = head;

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

  /* 6) Optional debug logging. */
  #if DEBUG_SQUARED_SIZES
    Rprintf("[C:c_squared_sizes] tail=%d head=%d | group=%d | edgestate=%d | "
            "deg_old=%d -> deg_new=%d | N_STATS=%d\n",
            (int)t, (int)h, (int)v2, (int)edgestate, deg_old, deg_new, N_CHANGE_STATS);
  #endif

  /* 7) For each sub-term j, read (from, to, pow) and accumulate Δ_j. */
  for(int j = 0; j < N_CHANGE_STATS; ++j){
    // Parameters for sub-term j: [from_j, to_j, pow_j].
    int from  = (int)INPUT_PARAM[3*j + 0];
    int to    = (int)INPUT_PARAM[3*j + 1];
    int power = (int)INPUT_PARAM[3*j + 2];

    double d = 0.0;

    // Contribution of the new degree.
    if(IN_RANGE(deg_new, from, to))
      d += ipow_int(deg_new, power);

    // Remove the contribution of the old degree.
    if(IN_RANGE(deg_old, from, to))
      d -= ipow_int(deg_old, power);

    // Accumulate local change for this j.
    CHANGE_STAT[j] += d;

    #if DEBUG_SQUARED_SIZES
      Rprintf("  stat[%d]: from=%d to=%d pow=%d | Δ=%.2f | cumul=%.2f\n",
              j, from, to, power, d, CHANGE_STAT[j]);
    #endif
  }
}
