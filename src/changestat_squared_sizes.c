/**
 * @file changestat_squared_sizes.c
 * @brief Change statistic for the ERPM term `squared_sizes` (multi-toggle form).
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
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
 *  A membership toggle affects exactly one group g and changes its
 *  size from deg_old to deg_new = deg_old ± 1. The local change in the
 *  statistic is therefore:
 *
 *      Δ = f(deg_new) − f(deg_old).
 *
 *  In multi-toggle mode, multiple toggles can affect the same group.
 *  Therefore, this D_ change-statistic processes toggles sequentially,
 *  temporarily applying each toggle to keep subsequent degrees consistent,
 *  and then undoes the temporary toggles before returning.
 *
 *  ------------------------------------------------------------
 *  Implementation in \pkg{ergm} (multi-toggle)
 *  ------------------------------------------------------------
 *
 *  - The function is declared with ::D_CHANGESTAT_FN and receives:
 *      - ntoggles : number of toggles in the proposal,
 *      - tails/heads arrays: endpoints of each toggle.
 *
 *  - For each toggle i:
 *      1. Read endpoints (t, h).
 *      2. Determine current edge state edgestate BEFORE toggling.
 *      3. Identify the affected group-mode vertex v2.
 *      4. Read deg_old from current network state.
 *      5. Compute deg_new = deg_old ± 1 based on edgestate.
 *      6. Accumulate Δ = f(deg_new) − f(deg_old).
 *      7. Temporarily apply the toggle if more toggles remain, so that
 *         later toggles see updated degrees.
 *
 *  - At the end, undo all temporary toggles, restoring the original network.
 *
 *  ------------------------------------------------------------
 *  Parameter packing (R side)
 *  ------------------------------------------------------------
 *
 *  INPUT_PARAM = c(pow, K, sizes_1, ..., sizes_K)
 *
 *    INPUT_PARAM[0] = pow
 *    INPUT_PARAM[1] = K
 *    INPUT_PARAM[2 + k] = sizes_{k+1}, k = 0..K-1
 */

#include <math.h>
#include <R_ext/Print.h>        // Rprintf
#include "ergm_changestat.h"
#include "ergm_storage.h"

/**
 * @def DEBUG_SQUARED_SIZES
 * @brief Enable or disable verbose debugging for ::d_squared_sizes.
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
#define UNUSED_VARIABLE(x) (void)(x)

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
    if(exp & 1) r *= b;
    b *= b;
    exp >>= 1;
  }
  return r;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: squared_sizes (multi-toggle)                              */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `squared_sizes(sizes, pow)` (multi-toggle).
 *
 * @details
 *  See file header for the mathematical definition and the multi-toggle logic.
 *
 *  Key design point:
 *    - Multi-toggle proposals (swap/split/merge decomposed into multiple edge
 *      toggles) must be handled consistently when several toggles affect the
 *      same group node. We therefore evaluate toggles sequentially and
 *      temporarily apply them (TOGGLE_IF_MORE_TO_COME) so that degree arrays
 *      reflect the intermediate state.
 */
D_CHANGESTAT_FN(d_squared_sizes){

#if DEBUG_SQUARED_SIZES
  static int seen = 0;
  if(ntoggles > 1 && seen < 10){
    Rprintf("[squared_sizes] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen++;
  }
#endif

  /* 1) Reset output buffer. */
  ZERO_ALL_CHANGESTATS();

  /* 2) Actor-mode size in bipartite networks. */
  const int n1 = BIPARTITE;

  /* 3) Read parameters. */
  const int power = (int)INPUT_PARAM[0];
  const int K     = (int)INPUT_PARAM[1];

  /* 4) Process toggles sequentially. */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex t = TAIL(i);
    Vertex h = HEAD(i);

    /* Determine current edge state before toggling.
     *
     * - For directed networks, use IS_OUTEDGE(t,h).
     * - For undirected (including typical bipartite memberships), we can use
     *   IS_UNDIRECTED_EDGE(t,h), which is robust to endpoint order.
     */
    int edgestate = DIRECTED ? IS_OUTEDGE(t, h) : IS_UNDIRECTED_EDGE(t, h);

    /* Identify the group-mode vertex (index > n1). */
    Vertex v2 = (t > n1) ? t : h;

#if DEBUG_SQUARED_SIZES
    if(v2 <= n1){
      Rprintf("[d_squared_sizes][WARN] toggle #%d has no group endpoint: tail=%d head=%d (n1=%d)\n",
              i, (int)t, (int)h, n1);
    }
#endif

    /* Read current group size from degree arrays (current intermediate state). */
    int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);

    /* Present -> deletion (-1), absent -> addition (+1). */
    int delta   = edgestate ? -1 : +1;
    int deg_new = deg_old + delta;

    /* Membership test: does deg belong to admissible sizes? */
    int match_old = 0;
    int match_new = 0;

    for(int k = 0; k < K; ++k){
      int size_k = (int)INPUT_PARAM[2 + k];
      if(deg_old == size_k) match_old = 1;
      if(deg_new == size_k) match_new = 1;
    }

    /* Compute local delta for this toggle under current intermediate state. */
    double d = 0.0;
    if(match_new) d += ipow_int(deg_new, power);
    if(match_old) d -= ipow_int(deg_old, power);

    CHANGE_STAT[0] += d;

#if DEBUG_SQUARED_SIZES
    Rprintf("[D:d_squared_sizes] i=%d tail=%d head=%d | group=%d | edgestate=%d | "
            "deg_old=%d -> deg_new=%d | K=%d pow=%d | Δ=%.2f | cumul=%.2f\n",
            i, (int)t, (int)h, (int)v2, (int)edgestate,
            deg_old, deg_new, K, power, d, CHANGE_STAT[0]);
#endif

    /* Temporarily apply this toggle so subsequent toggles see updated degrees. */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* 5) Undo temporary toggles to restore the original network state. */
  UNDO_PREVIOUS_TOGGLES(i);
}
