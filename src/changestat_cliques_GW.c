/**
 * @file changestat_cliques_GW.c
 * @brief Change statistic for the ERPM term `cliques_GW` (multi-toggle form).
 * 
 * @author Jérémie Chichignoud
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
 *  `cliques_GW(lambda)`, which aggregates group sizes through a geometrically
 *  weighted closed form.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (groups = group mode, actors = actor mode)
 *  ------------------------------------------------------------
 *
 *  A bipartite network is assumed, with:
 *    - actor mode  = vertices 1..BIPARTITE
 *    - group mode  = vertices (BIPARTITE+1)..N_NODES
 *
 *  If a group has size d (degree in the bipartite incidence), its closed-form
 *  GW contribution is:
 *
 *      S(d, λ) = λ * [1 - ((λ - 1) / λ)^d ]
 *
 *  Let r = (λ - 1) / λ.
 *  A membership toggle changes the group size from d_old to d_new, giving:
 *
 *      Δ = S(d_new, λ) − S(d_old, λ)
 *        = λ * ( r^{d_old} − r^{d_new} )
 *
 *  Only the group affected by the toggle contributes to the change.
 *
 *  ------------------------------------------------------------
 *  Implementation in \pkg{ergm} (multi-toggle / D_CHANGESTAT_FN)
 *  ------------------------------------------------------------
 *
 *  Why multi-toggle:
 *    - In ERPM workflows, MCMC proposals often represent a partition move
 *      (swap/split/merge) as a LIST of membership edge toggles.
 *    - If several toggles affect the same group within one proposal, degrees
 *      must be updated sequentially so each subsequent toggle sees the
 *      intermediate state.
 *
 *  Therefore this term is implemented as a D_ change-statistic:
 *    - It processes toggles sequentially.
 *    - It temporarily applies each toggle (TOGGLE_IF_MORE_TO_COME) so the degree
 *      arrays reflect the intermediate state.
 *    - At the end it undoes all temporary toggles (UNDO_PREVIOUS_TOGGLES),
 *      restoring the original network state before returning.
 *
 *  Vectorisation over λ:
 *    INPUT_PARAM = [λ_0, r_0, λ_1, r_1, ..., λ_{J-1}, r_{J-1}]
 *    where r_j = (λ_j − 1) / λ_j.
 *
 *  Output:
 *    For each λ_j:
 *      CHANGE_STAT[j] += λ_j * ( r_j^{d_old} − r_j^{d_new} )
 *
 *  R interface:
 *    Provided by InitErgmTerm.cliques_GW:
 *      - validates λ and builds INPUT_PARAM = c(rbind(lambda, r))
 *      - returns d_func=TRUE so ergm calls this D_ entrypoint
 */

#include <math.h>
#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"

/**
 * @def DEBUG_CLIQUES_GW
 * @brief Enable verbose debugging output for ::d_cliques_GW.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 * - toggle endpoints (tail/head),
 * - group vertex index and degrees before/after each toggle,
 * - λ and r values for each sub-term,
 * - intermediate powers r^d and the resulting Δ,
 * - ntoggles notice (multi-toggle detection).
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_CLIQUES_GW 0
#define UNUSED_VARIABLE(x) (void)x

/* -------------------------------------------------------------------------- */
/* Utility: fast exponentiation                                               */
/* -------------------------------------------------------------------------- */

/**
 * @brief Fast exponentiation for non-negative integer exponents.
 *
 * @details
 *  Computes base^exp using exponentiation-by-squaring (O(log exp)).
 *  Returns 1.0 when exp ≤ 0, implementing 0^0 = 1 for combinatorial
 *  consistency.
 *
 * @param base  The base in double precision.
 * @param exp   Non-negative integer exponent.
 * @return base^exp if exp > 0, else 1.0.
 */
static inline double dpow_double(double base, int exp){
  if(exp <= 0) return 1.0;  // convention base^0 = 1.0 for all base (including 0^0)

  double r = 1.0;
  double b = base;
  int e = exp;

  while(e){
    if(e & 1) r *= b;
    b *= b;
    e >>= 1;
  }
  return r;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cliques_GW (multi-toggle)                                 */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cliques_GW(lambda)` (multi-toggle).
 *
 * @details
 *  Processes a proposal composed of @c ntoggles membership toggles.
 *  For each toggle i:
 *   - determine current edge state BEFORE toggling,
 *   - identify the group-side endpoint v2,
 *   - read deg_old from the current intermediate state,
 *   - compute deg_new = deg_old ± 1,
 *   - for each λ_j accumulate Δ_j = λ_j( r_j^{deg_old} − r_j^{deg_new} ),
 *   - temporarily apply the toggle if more toggles remain.
 *
 *  At the end, undo all temporary toggles and return the accumulated changes.
 */
D_CHANGESTAT_FN(d_cliques_GW){

#if DEBUG_CLIQUES_GW
  static int seen = 0;
  if(ntoggles > 1 && seen < 10){
    Rprintf("[cliques_GW] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen++;
  }
#endif

  /* 1) Reset output buffer. */
  ZERO_ALL_CHANGESTATS();

  /* 2) Actor-mode size in bipartite networks. */
  const int n1 = BIPARTITE;

  /* 3) Process toggles sequentially. */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex t = TAIL(i);
    Vertex h = HEAD(i);

    /* Determine current edge state before toggling. */
    int edgestate = DIRECTED ? IS_OUTEDGE(t, h) : IS_UNDIRECTED_EDGE(t, h);

    /* Identify the group-mode endpoint (index > n1). */
    Vertex v2 = (t > n1) ? t : h;

#if DEBUG_CLIQUES_GW
    if(v2 <= n1){
      Rprintf("[d_cliques_GW][WARN] toggle #%d has no group endpoint: tail=%d head=%d (n1=%d)\n",
              i, (int)t, (int)h, n1);
    }
#endif

    /* Degrees before and after the toggle for this group. */
    int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);
    int delta   = edgestate ? -1 : +1;   /* present -> deletion (-1), absent -> addition (+1) */
    int deg_new = deg_old + delta;

#if DEBUG_CLIQUES_GW
    Rprintf("[d_cliques_GW] i=%d tail=%d head=%d | group=%d | edgestate=%d | deg_old=%d -> deg_new=%d\n",
            i, (int)t, (int)h, (int)v2, (int)edgestate, deg_old, deg_new);
#endif

    /* Loop over vectorized λ_j sub-terms.
     *
     * INPUT_PARAM is packed as:
     *   [lambda_0, r_0, lambda_1, r_1, ..., lambda_{J-1}, r_{J-1}]
     * where N_CHANGE_STATS == J.
     */
    for(int j = 0; j < N_CHANGE_STATS; ++j){
      double lambda = INPUT_PARAM[2*j + 0];
      double r      = INPUT_PARAM[2*j + 1];

      double term_old = dpow_double(r, deg_old);
      double term_new = dpow_double(r, deg_new);

      double d = lambda * (term_old - term_new);
      CHANGE_STAT[j] += d;

#if DEBUG_CLIQUES_GW
      Rprintf("  j=%d | lambda=%.9g r=%.9g | r^old=%.9g r^new=%.9g | Δ=%.9g | cumul=%.9g\n",
              j, lambda, r, term_old, term_new, d, CHANGE_STAT[j]);
#endif
    }

    /* Temporarily apply this toggle so subsequent toggles see updated degrees. */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* 4) Undo temporary toggles to restore the original network state. */
  UNDO_PREVIOUS_TOGGLES(i);
}
