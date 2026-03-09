/**
 * @file changestat_log_factorial_sizes.c
 * @brief Change statistic for the ERPM term `log_factorial_sizes` (MULTI-TOGGLE form).
 * 
 * @author Jérémie Chichignoud
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
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
 *  Implementation in \pkg{ergm} (MULTI-TOGGLE)
 *  ------------------------------------------------------------
 *
 *  Why multi-toggle?
 *  - Some proposal kernels (swap/split/merge) are decomposed into a list of
 *    edge toggles. In such a proposal, several toggles may touch the same
 *    group node. If we evaluate each toggle against the original degrees
 *    (without applying intermediate toggles), we get the wrong Δ.
 *
 *  Therefore:
 *  - This change-statistic is implemented as a D_CHANGESTAT_FN and processes
 *    the toggles sequentially.
 *  - After computing Δ for a toggle, we temporarily apply it so that degrees
 *    seen by subsequent toggles reflect the intermediate state.
 *  - At the end, we undo all temporary toggles, restoring the original network.
 *
 *  Mechanics:
 *    1) ZERO_ALL_CHANGESTATS()
 *    2) FOR_EACH_TOGGLE(i):
 *        - Read endpoints (TAIL/HEAD)
 *        - Determine edgestate BEFORE toggling
 *        - Identify affected group vertex v2
 *        - Read deg_old from current intermediate state
 *        - Compute Δ with the closed-form formulas
 *        - Accumulate CHANGE_STAT[0] += Δ
 *        - TOGGLE_IF_MORE_TO_COME(i)
 *    3) UNDO_PREVIOUS_TOGGLES(i)
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  - The R-side initialiser (InitErgmTerm.log_factorial_sizes) sets:
 *      - N_CHANGE_STATS = 1 (non-vectorised term),
 *      - no INPUT_PARAM (the term has no hyper-parameters at the C level),
 *      - emptynwstats = 0,
 *      - d_func = TRUE to select the D_ entrypoint.
 */

#include <math.h>
#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"

/**
 * @def DEBUG_LOG_FACTORIAL_SIZES
 * @brief Enable verbose debugging output for ::d_log_factorial_sizes.
 *
 * Set this macro to 1 to print diagnostic traces to the R console during
 * `summary()` or MCMC runs:
 *  - ntoggles (multi-toggle proposals),
 *  - endpoints of each toggle,
 *  - index of the affected group vertex,
 *  - degree before/after the toggle (intermediate state),
 *  - edge state and resulting increment Δ.
 *
 * When set to 0, no debug output is produced.
 */
#define DEBUG_LOG_FACTORIAL_SIZES 0
#define UNUSED_VARIABLE(x) (void)x

/* -------------------------------------------------------------------------- */
/* Change statistic: log_factorial_sizes (multi-toggle)                        */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `log_factorial_sizes` (multi-toggle).
 *
 * @details
 *  Implements the local change for:
 *      Stat = sum_{groups g} lgamma(deg(g)), with f(0)=0
 *  using closed-form O(1) increments and sequential toggle application.
 */
D_CHANGESTAT_FN(d_log_factorial_sizes){

#if DEBUG_LOG_FACTORIAL_SIZES
  static int seen = 0;
  if(ntoggles > 1 && seen < 10){
    Rprintf("[log_factorial_sizes] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen++;
  }
#endif

  /* 1) Reset output buffer for THIS proposal (list of toggles). */
  ZERO_ALL_CHANGESTATS();

  /* 2) Actor-mode size in bipartite networks. */
  const int n1 = BIPARTITE;

  /* 3) Process toggles sequentially (temporary apply to keep degrees consistent). */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex t = TAIL(i);
    Vertex h = HEAD(i);

    /* Determine current edge state before toggling. */
    int edgestate = DIRECTED ? IS_OUTEDGE(t, h) : IS_UNDIRECTED_EDGE(t, h);

    /* Identify the group-mode vertex (index > n1). */
    Vertex v2 = (t > n1) ? t : h;

#if DEBUG_LOG_FACTORIAL_SIZES
    if(v2 <= n1){
      Rprintf("[d_log_factorial_sizes][WARN] toggle #%d has no group endpoint: tail=%d head=%d (n1=%d)\n",
              i, (int)t, (int)h, n1);
    }
#endif

    /* Read current group size from degree arrays (intermediate state). */
    int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);

    /* Present -> deletion, absent -> addition. */
    int delta   = edgestate ? -1 : +1;
    int deg_new = deg_old + delta;
    UNUSED_VARIABLE(deg_new) ;

    /* Closed-form Δ:
     *  - addition:  Δ =  log(deg_old)      if deg_old >= 1, else 0
     *  - deletion:  Δ = -log(deg_old - 1)  if deg_old >= 2, else 0
     */
    double d = 0.0;
    if(edgestate == 0){
      if(deg_old >= 1) d = log((double)deg_old);
    }else{
      if(deg_old >= 2) d = -log((double)(deg_old - 1));
    }

    CHANGE_STAT[0] += d;

#if DEBUG_LOG_FACTORIAL_SIZES
    Rprintf("[D:d_log_factorial_sizes] i=%d tail=%d head=%d | group=%d | edgestate=%d | "
            "deg_old=%d -> deg_new=%d | Δ=%.12g | cumul=%.12g\n",
            i, (int)t, (int)h, (int)v2, (int)edgestate,
            deg_old, deg_new, d, CHANGE_STAT[0]);
#endif

    /* Temporarily apply this toggle so subsequent toggles see updated degrees. */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* 4) Undo temporary toggles to restore the original network state. */
  UNDO_PREVIOUS_TOGGLES(i);
}
