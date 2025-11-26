/**
 * @file
 * @brief  Change statistic for the ERPM term `dyadcov_full` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
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
 *  Implementation outline (one-toggle)
 *  ------------------------------------------------------------
 *
 *  The change statistic is computed in “one-toggle” form:
 *
 *    1. A single membership toggle connects one actor (actor mode) and
 *       one group (group mode).
 *
 *    2. For the affected group g:
 *         - reconstruct its actor members from the current network,
 *         - compute F_before = F_g(Z) under the current membership.
 *
 *    3. Apply a virtual toggle (::TOGGLE) for (actor, group):
 *         - reconstruct the actor members again,
 *         - compute F_after = F_g(Z) under the virtually updated membership.
 *
 *    4. Undo the virtual toggle (::TOGGLE) to restore the network.
 *
 *    5. The local change is:
 *
 *         Δ = F_after − F_before.
 *
 *    6. The scalar change statistic is updated as:
 *
 *         CHANGE_STAT[0] += Δ.
 *
 *  The {ergm} engine accumulates Δ over all toggles to obtain the
 *  total statistic during MCMC or summary evaluation.
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
#include "ergm_storage.h"      /* Calloc/Free */
#include <R_ext/Print.h>

/**
 * @def DEBUG_DYADCOV_FULL
 * @brief Enable verbose debugging output for ::c_dyadcov_full.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - group sizes and membership lists,
 *  - group-wise sums of dyadic covariates,
 *  - F_before, F_after and the local Δ.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_DYADCOV_FULL 0   /* set to 0 to disable debug output */

/**
 * @def UNUSED_WARNING
 * @brief Mark a parameter as intentionally unused.
 *
 * @param x Identifier of the unused variable.
 *
 * This macro is used to silence compiler warnings when a parameter is required
 * by the interface but not directly accessed in the implementation.
 */
#define UNUSED_WARNING(x) (void)x

/* -------------------------------------------------------------------------- */
/* Size filter: in_sizes                                                      */
/* -------------------------------------------------------------------------- */

/**
 * @brief Check whether a group size is contained in the allowed size set.
 *
 * @details
 *  The size filter S is represented by a vector @p sizes of length @p L.
 *  The effective rule is:
 *    - if L == 0, all sizes are accepted;
 *    - otherwise, a size @p n is accepted iff it matches one of the values
 *      in @p sizes (after casting to int).
 *
 * @param n      Group size to test.
 * @param L      Length of the filter vector @p sizes.
 * @param sizes  Pointer to a vector of length @p L containing allowed sizes.
 *
 * @return 1 if @p n is allowed, 0 otherwise.
 */
static inline int in_sizes(int n, int L, const double *sizes){
  if(L == 0) return 1;
  for(int i = 0; i < L; i++){
    if((int)sizes[i] == n) return 1;
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Group-level functional: group_dyadcov                                     */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute the dyadic covariate sum for one group vertex.
 *
 * @details
 *  For a group vertex @p g in the group mode, this function:
 *
 *    1. Identifies all actor neighbours connected to @p g via outgoing
 *       and incoming edges.
 *
 *    2. Counts them to obtain the group size n_g.
 *
 *    3. Applies the size filter:
 *         - if n_g <= 1, or
 *         - if n_g is not in the allowed size set S,
 *       then the contribution is 0.
 *
 *    4. Otherwise, collects the actor vertex indices in an array and
 *       computes:
 *
 *         ∑_{i != j, i,j ∈ A(g)} z_ij,
 *
 *       where Z is a n1 × n1 dyadic covariate matrix in column-major order.
 *
 *    For symmetric Z, this is equal to 2 × ∑_{i<j} z_ij. No division by 2 is
 *    applied.
 *
 * @param g       Group vertex whose actor members define the group.
 * @param n1      Number of actors (dimension of the actor mode).
 * @param L       Length of the size filter vector @p sizes.
 * @param sizes   Pointer to the vector of allowed group sizes (may have L=0).
 * @param Z       Pointer to the dyadic covariate matrix (n1*n1, column-major).
 * @param nwp     Pointer to the {ergm} Network structure (provides edges).
 *
 * @return The sum of z_ij over all ordered actor pairs i != j inside group
 *         @p g that satisfy the size filter, or 0.0 if n_g <= 1 or n_g is not
 *         allowed.
 */
static double group_dyadcov(Vertex g,
                            int n1, int L, const double *sizes,
                            const double *Z,
                            Network *nwp){

  /* Temporary bitmap of actor membership for this group (0/1 per actor). */
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char); /* initialised to 0 */

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
    Rprintf("[dyadcov_full][group_dyadcov] g=%d ng=%d -> 0 (empty/singleton/not in S)\n",
            (int)g, ng);
#endif
    Free(seen);
    return 0.0;
  }

  /* Collect the 1-based actor vertex indices belonging to g. */
  int *actors = (int*)Calloc(ng, int);
  int k = 0;
  for(int i = 0; i < n1; i++){
    if(seen[i]){
      actors[k++] = i + 1; /* actor vertex index in 1..n1 */
    }
  }

  double sum = 0.0;

  /* Sum over all ordered actor pairs i != j inside group g. */
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
  Rprintf("[dyadcov_full][group_dyadcov] g=%d ng=%d -> sum=%g\n",
          (int)g, ng, sum);
#endif

  Free(actors);
  Free(seen);

  return sum;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: dyadcov_full (one-toggle)                                */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `dyadcov_full`.
 *
 * @details
 *  This is the {ergm} change-statistic function registered as
 *  ::c_dyadcov_full via ::C_CHANGESTAT_FN. It implements the one-toggle
 *  update for the dyadic covariate statistic over actor pairs inside
 *  groups, with optional size filtering.
 *
 *  For a membership toggle (actor, group):
 *
 *    - It identifies the actor vertex (actor mode) and the group vertex
 *      (group mode) from the tail/head pair and the actor count n1.
 *
 *    - It computes F_before = F_g(Z) for the group in the current
 *      network using ::group_dyadcov().
 *
 *    - It applies a virtual toggle (::TOGGLE) for (actor, group) and
 *      computes F_after = F_g(Z) on the virtually updated network.
 *
 *    - It restores the network by toggling the edge back (::TOGGLE).
 *
 *    - The local change is:
 *
 *          Δ = F_after − F_before,
 *
 *      and this is added to the single scalar statistic:
 *
 *          CHANGE_STAT[0] += Δ.
 *
 * @param tail       Tail vertex of the toggled edge (actor or group).
 * @param head       Head vertex of the toggled edge (actor or group).
 * @param mtp        Pointer to the model term structure (unused directly
 *                   here but required by the macro signature).
 * @param nwp        Pointer to the network-plus workspace (used by
 *                   ::group_dyadcov to inspect edges).
 * @param edgestate  Current state of the edge:
 *                   - 0 if the edge is absent (toggle = addition),
 *                   - 1 if the edge is present (toggle = deletion).
 *
 * @note
 *  - The actor mode is determined by n1 from INPUT_PARAM; vertices with
 *    index ≤ n1 are actors, vertices with index > n1 are groups.
 *  - The term is non-vectorised: N_CHANGE_STATS == 1 and only
 *    CHANGE_STAT[0] is written.
 *  - The parameter @p edgestate is not used explicitly; the function
 *    relies on virtual toggling to obtain “before” and “after” values.
 */
C_CHANGESTAT_FN(c_dyadcov_full){
  /* 1) Reset the output buffer for THIS toggle. */
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  /* 2) Decode INPUT_PARAM layout:
   *    [0]      = n1
   *    [1]      = L
   *    [2..1+L] = sizes[L]
   *    [2+L..]  = Z[n1*n1] (column-major).
   */
  const double *ip    = INPUT_PARAM;
  const int n1        = (int)ip[0];
  const int L         = (int)ip[1];
  const double *sizes = ip + 2;
  const double *Z     = ip + 2 + L;

#if DEBUG_DYADCOV_FULL
  Rprintf("[dyadcov_full] n1=%d L=%d\n", n1, L);
#endif

  /* 3) Identify the actor and group vertices for this toggle. */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;
  Vertex group = (a <= (Vertex)n1) ? b : a;
  UNUSED_WARNING(actor);  /* actor is not used further explicitly. */

  /* 4) Group contribution BEFORE the virtual toggle. */
  double F_before = group_dyadcov(group, n1, L, sizes, Z, nwp);

  /* 5) Virtual toggle: temporarily change the membership edge. */
  TOGGLE(a, b);

  /* 6) Group contribution AFTER the virtual toggle. */
  double F_after = group_dyadcov(group, n1, L, sizes, Z, nwp);

  /* 7) Restore the original network by undoing the virtual toggle. */
  TOGGLE(a, b);

  /* 8) Compute and accumulate the local change Δ. */
  CHANGE_STAT[0] += (F_after - F_before);

#if DEBUG_DYADCOV_FULL
  Rprintf("[dyadcov_full] a=%d b=%d group=%d before=%g after=%g delta=%g\n",
          (int)a, (int)b, (int)group, F_before, F_after, (F_after - F_before));
#endif
}