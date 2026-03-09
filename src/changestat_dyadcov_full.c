/**
 * @file changestat_dyadcov_full.c
 * @brief  Change statistic for the ERPM term `dyadcov_full` (MULTI-TOGGLE form).
 * 
 * @author Jérémie Chichignoud
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
 *  `dyadcov_full`, which aggregates a dyadic covariate over all actor pairs
 *  inside each group, with an optional filter on group sizes.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  A bipartite network is assumed, with:
 *    - actor mode  = actor vertices,
 *    - group mode  = group vertices.
 *
 *  A numeric dyadic covariate matrix Z is defined on actor pairs:
 *    - dimension: n1 × n1, where n1 is the number of actors,
 *    - storage:   column-major order (R convention),
 *                 Z[(j-1)*n1 + (i-1)] = z_ij for actors i,j ∈ {1,…,n1}.
 *
 *  For a given group vertex g in the group mode, let:
 *    - A(g) be the set of actors connected to g (via incoming or outgoing edges),
 *    - n_g  = |A(g)| be the size of the group,
 *    - S be the (optional) set of allowed group sizes.
 *
 *  The group-level contribution is:
 *
 *      F_g(Z)
 *        = 1[ n_g ∈ S, n_g ≥ 2 ]
 *          * ∑_{i != j, i,j ∈ A(g)} z_ij.
 *
 *  The global statistic is:
 *
 *      T(p; Z)
 *        = ∑_g F_g(Z)
 *        = ∑_g 1[ n_g ∈ S, n_g ≥ 2 ]
 *            * ∑_{i != j, i,j ∈ A(g)} z_ij.
 *
 *  If Z is symmetric, this is equal to 2 × ∑_{i<j} z_ij: no division by 2 is
 *  applied in the implementation.
 *
 *  If S is empty (no size filter), all groups with n_g ≥ 2 contribute.
 *
 *  ------------------------------------------------------------
 *  Implementation outline (MULTI-TOGGLE / D_CHANGESTAT_FN)
 *  ------------------------------------------------------------
 *
 *  Why multi-toggle:
 *    - In ERPM / constrained bipartite proposals, one conceptual move
 *      (swap/split/merge) is often represented as a list of membership toggles.
 *    - Several toggles can affect the same group, so we must evaluate changes
 *      under the evolving intermediate state.
 *
 *  Strategy:
 *    - We compute, for each toggle, the group-level change:
 *
 *          Δ_i = F_after - F_before
 *
 *      where F_before is evaluated on the current intermediate network,
 *      and F_after is evaluated on the network with the toggle applied.
 *
 *    - For toggles 1..(ntoggles-1), we apply the toggle "for real" in the
 *      intermediate state using TOGGLE_IF_MORE_TO_COME(i), so subsequent toggles
 *      see updated memberships.
 *
 *    - For the last toggle, TOGGLE_IF_MORE_TO_COME(i) does not apply it (by design),
 *      so we apply a local virtual TOGGLE only to measure F_after, then undo it.
 *
 *    - At the end, UNDO_PREVIOUS_TOGGLES(i) restores the original network state
 *      by undoing the intermediate toggles applied via TOGGLE_IF_MORE_TO_COME.
 *
 *  ------------------------------------------------------------
 *  INPUT_PARAM layout
 *  ------------------------------------------------------------
 *
 *  The R-side initialiser passes the parameters as a flat numeric vector:
 *
 *      INPUT_PARAM = c(n1, L, sizes, as.vector(Z))
 *
 *  where:
 *    - n1    = number of actors (actor mode),
 *    - L     = length of the size filter vector `sizes`,
 *    - sizes = vector of allowed group sizes (may be empty),
 *    - Z     = n1 × n1 matrix in column-major order.
 *
 *  At the C level:
 *
 *      ip[0]          = n1
 *      ip[1]          = L
 *      ip[2 .. 1+L]   = sizes[0 .. L-1]
 *      ip[2+L .. end] = Z[0 .. (n1*n1-1)]
 *
 *  ------------------------------------------------------------
 *  Complexity and limitations
 *  ------------------------------------------------------------
 *
 *  For a group with n_g actors:
 *    - all actor pairs i,j distinct inside the group are considered,
 *    - the number of ordered pairs (i,j), i!=j, is n_g * (n_g - 1),
 *    - the cost per group is O(n_g^2).
 *
 *  The overall complexity per toggle is proportional to the number of
 *  actors in the affected group squared, plus the cost of scanning n1
 *  entries to reconstruct membership.
 *
 *  This term is intended for moderate group sizes and not for extremely
 *  large or dense bipartite components.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"      /* R_Calloc/R_Free */
#include <R_ext/Print.h>

/**
 * @def DEBUG_DYADCOV_FULL
 * @brief Enable verbose debugging output for ::d_dyadcov_full.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - multi-toggle detection (ntoggles>1),
 *  - toggle endpoints and touched group,
 *  - F_before / F_after / Δ per toggle.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_DYADCOV_FULL 0   /* set to 0 to disable debug output */

#define UNUSED_WARNING(x) (void)x

/* -------------------------------------------------------------------------- */
/* Size filter: in_sizes                                                      */
/* -------------------------------------------------------------------------- */

static inline int in_sizes(int n, int L, const double *sizes){
  if(L == 0) return 1;
  for(int i = 0; i < L; i++){
    if((int)sizes[i] == n) return 1;
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Group-level functional: group_dyadcov                                      */
/* -------------------------------------------------------------------------- */

static double group_dyadcov(Vertex g,
                            int n1, int L, const double *sizes,
                            const double *Z,
                            Network *nwp){

  /* Temporary bitmap of actor membership for this group (0/1 per actor). */
  unsigned char *seen = (unsigned char*)R_Calloc(n1, unsigned char); /* initialised to 0 */

  Vertex h;
  Edge e;

  /* Mark actor neighbours reachable via outgoing edges from g. */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      seen[idx] = 1;
    }
  }

  /* Mark actor neighbours reachable via incoming edges to g. */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      seen[idx] = 1;
    }
  }

  /* Count the number of actors in group g. */
  int ng = 0;
  for(int i = 0; i < n1; i++){
    if(seen[i]) ng++;
  }

  /* Empty/singleton group or filtered out by size -> zero contribution. */
  if(ng <= 1 || !in_sizes(ng, L, sizes)){
#if DEBUG_DYADCOV_FULL
    Rprintf("[dyadcov_full][group] g=%d ng=%d -> 0\n", (int)g, ng);
#endif
    R_Free(seen);
    return 0.0;
  }

  /* Collect the 1-based actor vertex indices belonging to g. */
  int *actors = (int*)R_Calloc(ng, int);
  int k = 0;
  for(int i = 0; i < n1; i++){
    if(seen[i]){
      actors[k++] = i + 1; /* actor vertex index in 1..n1 */
    }
  }

  double sum = 0.0;

  /* Sum over all ordered actor pairs i != j inside group g.
   *
   * We compute it as sum_{i<j} (z_ij + z_ji).
   */
  for(int p = 0; p < ng; p++){
    int i = actors[p];             /* 1..n1 */
    int row = i - 1;               /* 0..n1-1 */

    for(int q = p + 1; q < ng; q++){
      int j = actors[q];           /* 1..n1 */
      int col = j - 1;             /* 0..n1-1 */

      /* Z in column-major order (R style). */
      int idx_ij = col * n1 + row; /* z_ij */
      int idx_ji = row * n1 + col; /* z_ji */

      sum += Z[idx_ij] + Z[idx_ji];
    }
  }

#if DEBUG_DYADCOV_FULL
  Rprintf("[dyadcov_full][group] g=%d ng=%d -> sum=%g\n", (int)g, ng, sum);
#endif

  R_Free(actors);
  R_Free(seen);

  return sum;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: dyadcov_full (MULTI-TOGGLE)                               */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `dyadcov_full` (multi-toggle).
 *
 * @details
 *  This is the \pkg{ergm} change-statistic function registered as
 *  ::d_dyadcov_full via ::D_CHANGESTAT_FN.
 *
 *  The function processes the list of toggles sequentially. For each toggle,
 *  it computes the local change in the statistic for the unique touched group:
 *
 *      Δ = F_after − F_before
 *
 *  with:
 *    - F_before evaluated on the current intermediate network,
 *    - F_after evaluated on the network where the toggle is applied.
 *
 *  Intermediate state handling:
 *    - For toggles i < ntoggles-1:
 *        we apply the toggle via TOGGLE_IF_MORE_TO_COME(i), then evaluate F_after
 *        on the updated intermediate network.
 *    - For the last toggle:
 *        TOGGLE_IF_MORE_TO_COME(i) does not apply the toggle. We therefore apply
 *        a local virtual TOGGLE only to evaluate F_after, then immediately undo it.
 *
 *  Finally, UNDO_PREVIOUS_TOGGLES(i) restores the original network state by
 *  undoing the intermediate toggles applied by TOGGLE_IF_MORE_TO_COME.
 */
D_CHANGESTAT_FN(d_dyadcov_full){

#if DEBUG_DYADCOV_FULL
  static int seen_multi = 0;
  if(ntoggles > 1 && seen_multi < 20){
    Rprintf("[dyadcov_full] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen_multi++;
  }
#endif

  /* 1) Reset output buffer (vector length = N_CHANGE_STATS, here 1). */
  ZERO_ALL_CHANGESTATS();

  /* 2) Decode INPUT_PARAM layout. */
  const double *ip    = INPUT_PARAM;
  const int n1        = (int)ip[0];
  const int L         = (int)ip[1];
  const double *sizes = ip + 2;
  const double *Z     = ip + 2 + L;

  /* 3) Process toggles sequentially. */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex a = TAIL(i);
    Vertex b = HEAD(i);

    /* Identify actor/group endpoints for this toggle (robust to order). */
    Vertex actor = (a <= (Vertex)n1) ? a : b;
    Vertex group = (a <= (Vertex)n1) ? b : a;
    UNUSED_WARNING(actor);

    /* Evaluate before on the current intermediate state. */
    double F_before = group_dyadcov(group, n1, L, sizes, Z, nwp);

    /* Apply intermediate toggle if more toggles remain. */
    TOGGLE_IF_MORE_TO_COME(i);

    double F_after = 0.0;

    if(i < (int)ntoggles - 1){
      /* Toggle already applied to intermediate network. */
      F_after = group_dyadcov(group, n1, L, sizes, Z, nwp);
    } else {
      /* Last toggle: not applied by TOGGLE_IF_MORE_TO_COME. Do a local virtual toggle. */
      TOGGLE(a, b);
      F_after = group_dyadcov(group, n1, L, sizes, Z, nwp);
      TOGGLE(a, b);
    }

    double delta = (F_after - F_before);
    CHANGE_STAT[0] += delta;

#if DEBUG_DYADCOV_FULL
    Rprintf("[D:d_dyadcov_full] i=%d tail=%d head=%d | group=%d | before=%g after=%g Δ=%g | cumul=%g\n",
            i, (int)a, (int)b, (int)group, F_before, F_after, delta, CHANGE_STAT[0]);
#endif
  }

  /* 4) Restore original network state (undo intermediate toggles). */
  UNDO_PREVIOUS_TOGGLES(i);
}
