/* ============================================================================
 * File    : src/changestat_cov_ingroup.c
 * Purpose : Change statistic for the ERPM term `cov_ingroup` (MULTI-TOGGLE form).
 * Project : ERPM / ERGM extensions
 * ============================================================================
 *
 * @file changestat_cov_ingroup.c
 * @brief  Change statistic for the ERPM term `cov_ingroup` (multi-toggle form).
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
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
 *  ------------------------------------------------------------
 *  IMPORTANT (multi-toggle / D_CHANGESTAT_FN)
 *  ------------------------------------------------------------
 *
 *  This term MUST support multi-toggle moves (swap/split/merge proposals)
 *  represented as a list of toggles in ergm’s MCMC.
 *
 *  Therefore:
 *    - the compiled change-statistic is implemented using D_CHANGESTAT_FN,
 *      with the symbol name: d_cov_ingroup
 *    - the R initializer MUST return d_func = TRUE so ergm calls the
 *      multi-toggle entrypoint with the correct signature.
 *
 *  Multi-toggle correctness constraint:
 *    - multiple toggles may affect the same group in the same proposal;
 *    - we must process toggles sequentially and TEMPORARILY apply them
 *      (TOGGLE_IF_MORE_TO_COME) so subsequent toggles see updated degrees
 *      and neighbourhoods;
 *    - we must UNDO_PREVIOUS_TOGGLES at the end to restore the original state.
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
 *  Special case:
 *    - If L == 0, every group size is accepted (no filter).
 *
 *  The term returns a single scalar statistic:
 *    - N_CHANGE_STATS = 1,
 *    - CHANGE_STAT[0] is updated by the local Δ at each toggle (and accumulated
 *      over all toggles in the proposal).
 *
 *  ------------------------------------------------------------
 *  Complexity
 *  ------------------------------------------------------------
 *
 *  For each toggle involving group g (under the current intermediate state):
 *    - retrieving deg_old is O(1),
 *    - building X is O(deg(g)) via neighbour traversal,
 *    - all scalar operations are O(1).
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  The R initialiser (InitErgmTerm.cov_ingroup):
 *    - validates the actor covariate vector and the size filter S,
 *    - packs n1, L, sizes, and x into INPUT_PARAM,
 *    - sets d_func = TRUE (REQUIRED for this D_ changestat),
 *    - sets emptynwstats = 0 and a single coef.name.
 *
 *  ------------------------------------------------------------
 *  Debugging
 *  ------------------------------------------------------------
 *
 *  Compile-time debug macro DEBUG_COV_INGROUP:
 *    - if set to 1, prints per-toggle diagnostics.
 *    - prints an explicit "MULTI-TOGGLE ntoggles=..." banner (limited to a few).
 *  This is intentionally compile-time (like squared_sizes) to avoid runtime cost.
 * ============================================================================ */

#include <math.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>

/**
 * @def DEBUG_COV_INGROUP
 * @brief Enable or disable verbose debugging output for ::d_cov_ingroup.
 *
 * When set to 1, the change-statistic function prints diagnostic information
 * to the R console for each toggle:
 *  - actor and group vertex indices,
 *  - degree and covariate sums before and after the toggle,
 *  - size-filter flags and local Δ,
 *  - the updated CHANGE_STAT[0].
 *
 * When set to 0, no debug output is produced.
 */
#define DEBUG_COV_INGROUP 0

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
  if(L <= 0) return 1; /* no filter: all sizes accepted */
  for(int k = 0; k < L; ++k){
    if((int)in[k] == n) return 1;
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Helper: sum covariate over neighbours of a group vertex                    */
/* -------------------------------------------------------------------------- */

/**
 * @brief Sum x over all actor neighbours of a group vertex v2.
 *
 * @details
 *  Orientation-agnostic: traverses both outgoing and incoming edges of v2 and
 *  accumulates x on actor-mode vertices 1..n1.
 *
 *  This must be correct under the CURRENT INTERMEDIATE STATE in multi-toggle
 *  mode (i.e., after some previous toggles have been temporarily applied).
 *
 * @param v2  Group vertex index (group mode).
 * @param x   Pointer to covariate values x[0..n1-1] on actors.
 * @param n1  Number of actors (size of actor mode).
 * @param nwp Network workspace pointer (required by ergm macros).
 *
 * @return Sum of x_i over all actor neighbours i of v2.
 */
static inline double SAFE_SUM_GROUP(Vertex v2, const double *x, int n1, Network *nwp){
  (void)nwp;

  double X = 0.0;
  Vertex h;
  Edge e;

  STEP_THROUGH_OUTEDGES(v2, e, h){
    if(h >= (Vertex)1 && h <= (Vertex)n1){
      X += x[(int)h - 1];
    }
  }

  STEP_THROUGH_INEDGES(v2, e, h){
    if(h >= (Vertex)1 && h <= (Vertex)n1){
      X += x[(int)h - 1];
    }
  }

  return X;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_ingroup (multi-toggle)                               */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_ingroup` (multi-toggle).
 *
 * @details
 *  Implemented as ::D_CHANGESTAT_FN(d_cov_ingroup).
 *
 *  For each toggle in the proposal (under the current intermediate state):
 *    1. Identify the actor vertex v1 and the group vertex v2.
 *    2. Read deg_old (current group size).
 *    3. Recompute X (current sum of x in the group).
 *    4. Infer n_new and X_new after applying this toggle.
 *    5. Apply the size filter on deg_old and n_new.
 *    6. Accumulate Δ = (n_new*X_new*w_new) - (deg_old*X*w_old).
 *    7. Temporarily apply the toggle so later toggles see updated state.
 *
 *  At the end: UNDO_PREVIOUS_TOGGLES to restore the original network state.
 */
D_CHANGESTAT_FN(d_cov_ingroup){

#if DEBUG_COV_INGROUP
  static int seen = 0;
  if(ntoggles > 1 && seen < 10){
    Rprintf("[cov_ingroup] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen++;
  }
#endif

  /* 1) Reset output buffer for the whole proposal. */
  ZERO_ALL_CHANGESTATS();

  /* 2) Read inputs from INPUT_PARAM. */
  const int n1 = (int)INPUT_PARAM[0];  /* number of actors (actor mode size) */
  const int L  = (int)INPUT_PARAM[1];  /* length of size filter S */

  const double *sizes = (L > 0) ? (&INPUT_PARAM[2])      : NULL;
  const double *x     = (L > 0) ? (&INPUT_PARAM[2 + L])  : (&INPUT_PARAM[2]);

  /* 3) Actor/group boundary. */
  const Vertex n1_lim = (Vertex)BIPARTITE; /* should coincide with n1 */

  /* 4) Process toggles sequentially (multi-toggle). */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex t = TAIL(i);
    Vertex h = HEAD(i);

    /* Determine current edge state BEFORE toggling (intermediate state). */
    int edgestate = DIRECTED ? IS_OUTEDGE(t, h) : IS_UNDIRECTED_EDGE(t, h);

    /* Identify endpoints: v2 = group, v1 = actor. */
    Vertex v2 = (t > n1_lim) ? t : h;
    Vertex v1 = (t > n1_lim) ? h : t;

    /* Safety: must toggle an actor-group edge. */
    if(v1 < (Vertex)1 || v1 > (Vertex)n1 || v2 <= n1_lim){
      /* Ignore invalid toggle; still apply temp toggle for consistency? No. */
#if DEBUG_COV_INGROUP
      Rprintf("[d_cov_ingroup][WARN] invalid toggle #%d: tail=%d head=%d (n1=%d bip=%d)\n",
              i, (int)t, (int)h, n1, (int)n1_lim);
#endif
      TOGGLE_IF_MORE_TO_COME(i);
      continue;
    }

    /* Current group size and ingroup sum under the intermediate state. */
    int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);
    double X    = SAFE_SUM_GROUP(v2, x, n1, nwp);

    /* Addition if absent; deletion if present. */
    const int is_add   = (edgestate == 0);
    const int n_new    = deg_old + (is_add ? +1 : -1);
    const double xi    = x[(int)v1 - 1];
    const double X_new = X + (is_add ? +xi : -xi);

    /* Apply size filter. */
    const int w_old = in_sizes_set(deg_old, sizes, L);
    const int w_new = in_sizes_set(n_new,  sizes, L);

    /* Local delta for this toggle under current intermediate state. */
    double d = 0.0;
    d = (w_new ? ((double)n_new * X_new) : 0.0)
      - (w_old ? ((double)deg_old * X)   : 0.0);

    CHANGE_STAT[0] += d;

#if DEBUG_COV_INGROUP
    Rprintf(
      "[D:d_cov_ingroup] i=%d tail=%d head=%d | v1=%d v2=%d | edgestate=%d is_add=%d | "
      "deg_old=%d -> n_new=%d | xi=%.6f | X=%.6f -> X_new=%.6f | "
      "w_old=%d w_new=%d | Δ=%.6f | cumul=%.6f\n",
      i, (int)t, (int)h, (int)v1, (int)v2, (int)edgestate, is_add,
      deg_old, n_new, xi, X, X_new, w_old, w_new, d, CHANGE_STAT[0]
    );
#endif

    /* Temporarily apply this toggle so subsequent toggles see updated state. */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* 5) Undo temporary toggles to restore the original network state. */
  UNDO_PREVIOUS_TOGGLES(i);
}
