/**
 * @file
 * @brief  Change statistic for the ERPM term `dyadcov_GW` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `dyadcov_GW(lambda)`, which aggregates a symmetrised dyadic covariate
 *  functional over all actor cliques inside each group, with a geometric
 *  weighting over the clique size.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  A bipartite network is assumed, with:
 *    - actor mode  = actor vertices,
 *    - group mode  = group vertices.
 *
 *  Each group vertex represents a structural group. Its members are the
 *  actors connected to it by membership edges (outgoing or incoming).
 *
 *  A numeric dyadic covariate matrix Z is defined on actor pairs:
 *    - dimension: n1 × n1, where n1 is the number of actors,
 *    - storage:   column-major order (R convention),
 *                 Z[(j-1)*n1 + (i-1)] = z_ij for actors i,j.
 *
 *  The matrix Z is not required to be symmetric: in general z_ij != z_ji.
 *  For each unordered actor pair {i,j} with i<j, the implementation uses
 *  the symmetrised value
 *
 *      w_ij = z_ij + z_ji,
 *
 *  in line with the definition of `dyadcov_full`.
 *
 *  For a given group g, let:
 *    - A(g) be the set of actors in group g,
 *    - n_g  = |A(g)| be the size of group g,
 *    - C_k(g) be the set of all k-subsets C ⊂ A(g).
 *
 *  For each clique C ∈ C_k(g), define the dyadic product based on w_ij:
 *
 *      P(C; Z) = ∏_{i<j ∈ C} w_ij
 *              = ∏_{i<j ∈ C} (z_ij + z_ji).
 *
 *  The group-level k-clique covariate functional is:
 *
 *      S_g^{(k)}(Z) = ∑_{C ∈ C_k(g)} P(C; Z).
 *
 *  The geometrically weighted dyadic covariate functional for a group is:
 *
 *      S_g^{GW}(Z, λ)
 *        = ∑_{k=2}^{n_g} a_k(λ) * S_g^{(k)}(Z),
 *
 *  with weights:
 *
 *      a_k(λ) = (-1 / λ)^{k-1},   for k ≥ 2.
 *
 *  The global statistic is:
 *
 *      T_GW(p; Z, λ)
 *        = ∑_g S_g^{GW}(Z, λ).
 *
 *  For λ = 2, the weights are:
 *      a_2 = 1, a_3 = -1/2, a_4 = 1/4, ...
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
 *         - compute n_g and S_g^{GW}(Z, λ) before the virtual toggle.
 *
 *    3. Apply a virtual toggle (::TOGGLE) for the edge (actor, group):
 *         - reconstruct the actor members again,
 *         - compute n_g and S_g^{GW}(Z, λ) after the virtual toggle.
 *
 *    4. Undo the virtual toggle (::TOGGLE again) to restore the network.
 *
 *    5. The local change is:
 *
 *         Δ = S_after − S_before.
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
 *      INPUT_PARAM = c(n1, lambda, as.vector(Z))
 *
 *  where:
 *    - n1     = number of actors (actor mode),
 *    - lambda = geometric decay parameter λ,
 *    - Z      = n1 × n1 matrix in column-major order.
 *
 *  At the C level:
 *
 *      ip[0]   = n1
 *      ip[1]   = lambda
 *      ip[2+]  = Z[0..(n1*n1-1)]
 *
 *  ------------------------------------------------------------
 *  Complexity and limitations
 *  ------------------------------------------------------------
 *
 *  For a group with n_g actors:
 *    - all clique sizes k = 2..n_g are considered,
 *    - the number of k-cliques is C(n_g, k),
 *    - each clique product involves k*(k-1)/2 unordered dyads, each mapped
 *      to a symmetrised value w_ij = z_ij + z_ji,
 *    - S_g^{GW}(Z, λ) is a sum of S_g^{(k)}(Z) weighted by a_k(λ).
 *
 *  This leads to combinatorial cost in n_g and k. The term is therefore
 *  intended for moderate group sizes and not for very large or dense
 *  actor sets.
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  The R initialiser (InitErgmTerm.dyadcov_GW):
 *    - validates n1, λ and the dimensions of Z,
 *    - flattens Z in column-major order,
 *    - sets N_CHANGE_STATS = 1 and emptynwstats = 0,
 *    - builds INPUT_PARAM as c(n1, lambda, as.vector(Z)).
 *
 *  ------------------------------------------------------------
 *  @example Usage (R)
 *  ------------------------------------------------------------
 *  @code{.r}
 *  library(ERPM)
 *
 *  # Example: 4 actors in 2 groups
 *  part <- c(1, 1, 2, 2)  # actors 1,2 in group 1; 3,4 in group 2
 *
 *  # Possibly non-symmetric dyadic covariate on the actor mode
 *  set.seed(1)
 *  Z <- matrix(runif(4 * 4), nrow = 4, ncol = 4)
 *  diag(Z) <- 0
 *
 *  # Geometrically weighted dyadic covariate effect
 *  fit <- erpm(partition ~ dyadcov_GW(lambda = 2, dyadcov = Z))
 *  summary(fit)
 *  @endcode
 *
 *  @test
 *  A self-test for this change statistic can:
 *    - construct small bipartite networks from known partitions;
 *    - compute T_GW(p; Z, λ) directly in R by:
 *         * enumerating all groups g and actor sets A(g),
 *         * enumerating all k-cliques C ⊂ A(g) for k ≥ 2,
 *         * computing w_ij = z_ij + z_ji and summing
 *           a_k(λ) * ∏_{i<j ∈ C} w_ij;
 *    - compare these reference values to:
 *         summary( erpm(partition ~ dyadcov_GW(...)) )$statistics;
 *    - apply explicit edge toggles and check that observed changes match
 *      the Δ produced by ::c_dyadcov_GW.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>

/**
 * @def DEBUG_DYADCOV_GW
 * @brief Enable verbose debugging output for ::c_dyadcov_GW.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - actor and group sizes per group,
 *  - intermediate S_g^{(k)}(Z) values and geometric weights a_k(λ),
 *  - S_before, S_after and the local Δ.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_DYADCOV_GW 0

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
/* Clique enumeration: sum_cliques_k                                         */
/* -------------------------------------------------------------------------- */

/**
 * @brief Sum of products over all k-cliques of actors inside a group.
 *
 * @details
 *  Given a vector of actor indices and a clique size k, this function:
 *    - enumerates all k-subsets C of the actors in the group,
 *    - for each C, computes:
 *
 *         P(C; Z) = ∏_{i<j ∈ C} w_ij,
 *
 *      where w_ij = z_ij + z_ji is the symmetrised dyadic covariate for
 *      the unordered pair {i,j}, using the matrix Z on the actor mode,
 *    - accumulates the total value:
 *
 *         ∑_{C} P(C; Z).
 *
 *  The matrix Z is indexed in column-major order (as in R):
 *
 *      Z[(j-1)*n1 + (i-1)] = z_ij.
 *
 * @param actors  Pointer to an array of length @p ng containing actor
 *                vertex indices in the actor mode (1-based).
 * @param ng      Number of actors in the group.
 * @param k       Clique size (k ≥ 2).
 * @param n1      Number of actors (dimension of the actor mode).
 * @param Z       Pointer to the dyadic covariate matrix (length n1*n1),
 *                stored in column-major order.
 *
 * @return The sum of products over all k-cliques of actors in this group.
 */
static double sum_cliques_k(const int *actors,
                            int ng, int k,
                            int n1,
                            const double *Z){

  if(k > ng) return 0.0;

  /* Working array for k-combinations of indices into actors[]. */
  int *comb = (int*)Calloc(k, int);
  for(int i = 0; i < k; i++) comb[i] = i;

  double total = 0.0;

  while(1){
    /* Product over all unordered actor dyads {i,j} in the current clique. */
    double prod = 1.0;
    for(int p = 0; p < k; p++){
      int idx_i = actors[ comb[p] ];   /* actor vertex index 1..n1 */
      int row   = idx_i - 1;           /* 0..n1-1 */

      for(int q = p + 1; q < k; q++){
        int idx_j = actors[ comb[q] ];
        int col   = idx_j - 1;

        int idx_ij = col * n1 + row;   /* z_ij   */
        int idx_ji = row * n1 + col;   /* z_ji   */
        double w_ij = Z[idx_ij] + Z[idx_ji];

        prod *= w_ij;
      }
    }
    total += prod;

    /* Generate the next k-combination in lexicographic order. */
    int pos = k - 1;
    while(pos >= 0 && comb[pos] == (ng - k + pos)) pos--;
    if(pos < 0) break; /* no more combinations */

    comb[pos]++;
    for(int j = pos + 1; j < k; j++){
      comb[j] = comb[j - 1] + 1;
    }
  }

  Free(comb);
  return total;
}

/* -------------------------------------------------------------------------- */
/* Group-level functional: group_dyadcov_GW                                  */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute S_g^{GW}(Z, λ) for a given group vertex.
 *
 * @details
 *  For a group vertex \p g this function:
 *
 *    1. Collects all distinct actor neighbours connected to g via outgoing
 *       and incoming edges.
 *
 *    2. Counts them (n_g) and optionally stores that value into @p *n_g_out
 *       if @p n_g_out is non-NULL.
 *
 *    3. If n_g < 2, no k-cliques exist and the contribution is 0.
 *
 *    4. Otherwise:
 *         - builds a list of actor vertex indices in the actor mode,
 *         - for each k from 2 to n_g:
 *             * computes S_g^{(k)}(Z) via ::sum_cliques_k(), using
 *               w_ij = z_ij + z_ji for each unordered pair {i,j},
 *             * multiplies it by the weight a_k(λ) = (-1/λ)^{k-1},
 *             * accumulates the result into S_g^{GW}(Z, λ).
 *
 * @param g        Group vertex whose actor members define the group.
 * @param n1       Number of actors (dimension of the actor mode).
 * @param lambda   Geometric decay parameter λ.
 * @param Z        Pointer to the dyadic covariate matrix (n1*n1, column-major).
 * @param nwp      Pointer to the {ergm} Network structure (provides edges).
 * @param n_g_out  If non-NULL, receives the number of actors in group g.
 *
 * @return S_g^{GW}(Z, λ) for group g.
 */
static double group_dyadcov_GW(Vertex g,
                               int n1,
                               double lambda,
                               const double *Z,
                               Network *nwp,
                               int *n_g_out){

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
  if(n_g_out) *n_g_out = ng;

  /* If the group has fewer than 2 actors, no cliques are possible. */
  if(ng < 2){
#if DEBUG_DYADCOV_GW
    Rprintf("[dyadcov_GW][group_dyadcov_GW] g=%d ng=%d < 2 -> 0\n",
            (int)g, ng);
#endif
    Free(seen);
    return 0.0;
  }

  /* Collect the 1-based actor vertex indices belonging to g. */
  int *actors = (int*)Calloc(ng, int);
  int idx = 0;
  for(int i = 0; i < n1; i++){
    if(seen[i]){
      actors[idx++] = i + 1;  /* actor vertex index in 1..n1 */
    }
  }

  /* Geometrically weighted sum over clique sizes k = 2..ng. */
  double sum_gw = 0.0;
  double factor = 1.0;          /* a_2(λ) = 1 = (-1/λ)^{1-1} */

  for(int k = 2; k <= ng; k++){
    double S_k = sum_cliques_k(actors, ng, k, n1, Z);
    sum_gw += factor * S_k;
#if DEBUG_DYADCOV_GW
    Rprintf("[dyadcov_GW][group_dyadcov_GW] g=%d ng=%d k=%d S_k=%g factor=%g\n",
            (int)g, ng, k, S_k, factor);
#endif
    factor *= (-1.0 / lambda); /* recurrence: a_{k+1} = a_k * (-1/λ) */
  }

  Free(actors);
  Free(seen);

#if DEBUG_DYADCOV_GW
  Rprintf("[dyadcov_GW][group_dyadcov_GW] g=%d ng=%d -> sum_gw=%g\n",
          (int)g, ng, sum_gw);
#endif

  return sum_gw;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: dyadcov_GW (one-toggle)                                  */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `dyadcov_GW(lambda)`.
 *
 * @details
 *  This is the {ergm} change-statistic function registered as
 *  ::c_dyadcov_GW via ::C_CHANGESTAT_FN. It implements the one-toggle
 *  update for the geometrically weighted dyadic covariate statistic on
 *  actor cliques inside each group, based on the symmetrised covariate
 *  w_ij = z_ij + z_ji for each unordered pair {i,j}.
 *
 *  For a membership toggle (actor, group):
 *
 *    - It identifies the actor vertex (actor mode) and the group vertex
 *      (group mode) from the tail/head pair and the actor count n1.
 *
 *    - It computes S_before = S_g^{GW}(Z, λ) for the group in the current
 *      network using ::group_dyadcov_GW().
 *
 *    - It applies a virtual toggle (::TOGGLE) for (actor, group) and
 *      computes S_after = S_g^{GW}(Z, λ) on the virtually updated network.
 *
 *    - It restores the network by toggling the edge back (::TOGGLE).
 *
 *    - The local change is:
 *
 *          Δ = S_after − S_before,
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
 *                   ::group_dyadcov_GW to inspect edges).
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
C_CHANGESTAT_FN(c_dyadcov_GW){
  /* 1) Reset the output buffer for THIS toggle.
   *
   * {ergm} accumulates contributions from multiple calls; here we only
   * report the local Δ for the current membership toggle.
   */
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  /* 2) Decode INPUT_PARAM layout: [n1, lambda, Z...]. */
  const double *ip     = INPUT_PARAM;
  const int     n1     = (int)ip[0];   /* number of actors (actor mode) */
  const double  lambda = ip[1];        /* geometric decay parameter λ   */
  const double *Z      = ip + 2;       /* dyadic covariate matrix       */

#if DEBUG_DYADCOV_GW
  Rprintf("[dyadcov_GW] n1=%d lambda=%g\n", n1, lambda);
#endif

  /* If λ = 0, the weights are not defined; the statistic is set to 0. */
  if(lambda == 0.0){
    CHANGE_STAT[0] = 0.0;
    return;
  }

  /* 3) Identify the actor and group vertices for this toggle.
   *
   * By construction:
   *   - exactly one endpoint has index ≤ n1 (actor),
   *   - the other endpoint has index > n1 (group).
   */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;
  Vertex group = (a <= (Vertex)n1) ? b : a;
  UNUSED_WARNING(actor);  /* currently unused beyond identification. */

  int ng_before = 0, ng_after = 0;

  /* 4) Group contribution BEFORE the virtual toggle. */
  double S_before = group_dyadcov_GW(group, n1, lambda, Z, nwp, &ng_before);

  /* 5) Virtual toggle: temporarily change the membership edge. */
  TOGGLE(a, b);

  /* 6) Group contribution AFTER the virtual toggle. */
  double S_after = group_dyadcov_GW(group, n1, lambda, Z, nwp, &ng_after);

  /* 7) Restore the original network by undoing the virtual toggle. */
  TOGGLE(a, b);

  /* 8) Compute and accumulate the local change Δ. */
  double delta = S_after - S_before;
  CHANGE_STAT[0] += delta;

#if DEBUG_DYADCOV_GW
  Rprintf("[dyadcov_GW] a=%d b=%d group=%d ng_before=%d ng_after=%d "
          "S_before=%g S_after=%g delta=%g\n",
          (int)a, (int)b, (int)group,
          ng_before, ng_after, S_before, S_after, delta);
#endif
}
