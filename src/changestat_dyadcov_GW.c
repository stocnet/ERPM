/**
 * @file changestat_dyadcov_GW.c
 * @brief  Change statistic for the ERPM term `dyadcov_GW` (MULTI-toggle form).
 * 
 * @author Jérémie Chichignoud
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
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
 *  Implementation outline (MULTI-toggle / D_CHANGESTAT_FN)
 *  ------------------------------------------------------------
 *
 *  IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
 *  - This term MUST support proposals that consist of multiple edge toggles
 *    (swap/split/merge decomposed into a list of toggles).
 *  - Therefore, the compiled change-statistic MUST be implemented using the
 *    D_CHANGESTAT_FN API (multi-toggle).
 *  - On the R side, we MUST advertise this to ergm by returning `d_func = TRUE`.
 *    Otherwise ergm will try to call the changestat as a one-toggle C_CHANGESTAT_FN,
 *    causing a signature mismatch and typically a segfault.
 *
 *  Multi-toggle logic used here:
 *    For each toggle i (tail/head):
 *      1) Compute S_before for the affected group under the CURRENT intermediate state.
 *      2) Virtually apply the toggle and compute S_after.
 *      3) Δ_i = S_after - S_before is added to CHANGE_STAT[0].
 *      4) If more toggles remain, we APPLY the toggle to the network state
 *         (so subsequent toggles see updated memberships), and later UNDO them
 *         all at the end of the function.
 *
 *  This is correct even if multiple toggles in the same proposal affect the
 *  same group vertex: each Δ_i is computed against the proper intermediate state.
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
 *    - builds INPUT_PARAM as c(n1, lambda, as.vector(Z)),
 *    - IMPORTANT: sets d_func=TRUE to use this D_ entrypoint.
 *
 *  ------------------------------------------------------------
 *  @test
 *  A self-test for this change statistic can:
 *    - validate summary() equivalences on explicit networks;
 *    - validate erpm() fits return finite coefficients;
 *    - trigger MCMC with a proposal that can emit multi-toggle steps, and
 *      observe C-level debug traces when DEBUG_DYADCOV_GW is enabled.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>

/**
 * @def DEBUG_DYADCOV_GW
 * @brief Enable verbose debugging output for ::d_dyadcov_GW.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - multi-toggle detection (ntoggles>1),
 *  - per-toggle endpoints and affected group,
 *  - S_before, S_after and the local Δ,
 *  - group sizes before/after (ng_before/ng_after).
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_DYADCOV_GW 0

/**
 * @def UNUSED_WARNING
 * @brief Mark a parameter as intentionally unused.
 *
 * @param x Identifier of the unused variable.
 */
#define UNUSED_WARNING(x) (void)x

/* -------------------------------------------------------------------------- */
/* Clique enumeration: sum_cliques_k                                           */
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
  int *comb = (int*)R_Calloc(k, int);
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

  R_Free(comb);
  return total;
}

/* -------------------------------------------------------------------------- */
/* Group-level functional: group_dyadcov_GW                                    */
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
 * @param nwp      Pointer to the \pkg{ergm} Network structure (provides edges).
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
  unsigned char *seen = (unsigned char*)R_Calloc(n1, unsigned char);
  Vertex h;
  Edge e;

  /* Mark actor neighbours reachable via outgoing edges from g. */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      seen[(int)h - 1] = 1;
    }
  }

  /* Mark actor neighbours reachable via incoming edges to g. */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      seen[(int)h - 1] = 1;
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
    Rprintf("[dyadcov_GW][group] g=%d ng=%d < 2 -> 0\n", (int)g, ng);
#endif
    R_Free(seen);
    return 0.0;
  }

  /* Collect the 1-based actor vertex indices belonging to g. */
  int *actors = (int*)R_Calloc(ng, int);
  int idx = 0;
  for(int i = 0; i < n1; i++){
    if(seen[i]) actors[idx++] = i + 1;
  }

  /* Geometrically weighted sum over clique sizes k = 2..ng. */
  double sum_gw = 0.0;
  double factor = 1.0; /* a_2(λ) = 1; recurrence factor *= (-1/λ) */

  for(int k = 2; k <= ng; k++){
    double S_k = sum_cliques_k(actors, ng, k, n1, Z);
    sum_gw += factor * S_k;

#if DEBUG_DYADCOV_GW
    Rprintf("[dyadcov_GW][group] g=%d ng=%d k=%d S_k=%g factor=%g\n",
            (int)g, ng, k, S_k, factor);
#endif

    factor *= (-1.0 / lambda);
  }

  R_Free(actors);
  R_Free(seen);

#if DEBUG_DYADCOV_GW
  Rprintf("[dyadcov_GW][group] g=%d ng=%d -> sum_gw=%g\n",
          (int)g, ng, sum_gw);
#endif

  return sum_gw;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: dyadcov_GW (MULTI-toggle)                                 */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `dyadcov_GW(lambda)` (multi-toggle).
 *
 * @details
 *  This is the \pkg{ergm} D_ change-statistic entrypoint registered as
 *  ::d_dyadcov_GW via ::D_CHANGESTAT_FN. It implements the multi-toggle
 *  update for the geometrically weighted dyadic covariate statistic on
 *  actor cliques inside each group, based on the symmetrised covariate
 *  w_ij = z_ij + z_ji for each unordered pair {i,j}.
 *
 *  The function processes toggles sequentially. For each toggle, the affected
 *  group is identified, and Δ_i is computed as:
 *
 *      Δ_i = S_after(intermediate + toggle_i) - S_before(intermediate)
 *
 *  where "intermediate" is the network state after applying previous toggles
 *  from the same proposal. This is achieved by:
 *    - computing S_before,
 *    - applying a virtual TOGGLE and computing S_after,
 *    - undoing that virtual TOGGLE,
 *    - then applying TOGGLE_IF_MORE_TO_COME(i) to carry the intermediate state
 *      forward when needed.
 *
 *  At the end, UNDO_PREVIOUS_TOGGLES(i) restores the original network.
 */
D_CHANGESTAT_FN(d_dyadcov_GW){

#if DEBUG_DYADCOV_GW
  static int seen_mt = 0;
  if(ntoggles > 1 && seen_mt < 20){
    Rprintf("[dyadcov_GW] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen_mt++;
  }
#endif

  /* 1) Reset output buffer for THIS proposal. */
  ZERO_ALL_CHANGESTATS();

  /* 2) Decode INPUT_PARAM layout: [n1, lambda, Z...]. */
  const double *ip     = INPUT_PARAM;
  const int     n1     = (int)ip[0];
  const double  lambda = ip[1];
  const double *Z      = ip + 2;

  /* Guard: λ must be non-zero (weights undefined at 0). */
  if(lambda == 0.0){
    CHANGE_STAT[0] = 0.0;
    return;
  }

  /* 3) Process toggles sequentially (multi-toggle). */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex tail = TAIL(i);
    Vertex head = HEAD(i);

    /* Identify the affected group vertex (index > n1). */
    Vertex actor = (tail <= (Vertex)n1) ? tail : head;
    Vertex group = (tail <= (Vertex)n1) ? head : tail;

#if DEBUG_DYADCOV_GW
    if(group <= (Vertex)n1){
      Rprintf("[dyadcov_GW][WARN] toggle #%d has no group endpoint: tail=%d head=%d (n1=%d)\n",
              (int)i, (int)tail, (int)head, n1);
    }
#endif
    UNUSED_WARNING(actor);

    int ng_before = 0, ng_after = 0;

    /* S_before under the current intermediate state. */
    double S_before = group_dyadcov_GW(group, n1, lambda, Z, nwp, &ng_before);

    /* Virtual toggle (compute after). */
    TOGGLE(tail, head);
    double S_after  = group_dyadcov_GW(group, n1, lambda, Z, nwp, &ng_after);
    TOGGLE(tail, head); /* restore */

    double delta = S_after - S_before;
    CHANGE_STAT[0] += delta;

#if DEBUG_DYADCOV_GW
    {
      int edgestate = DIRECTED ? IS_OUTEDGE(tail, head) : IS_UNDIRECTED_EDGE(tail, head);
      Rprintf("[dyadcov_GW][D] i=%d tail=%d head=%d | edgestate=%d | group=%d | "
              "ng_before=%d ng_after=%d | S_before=%g S_after=%g | Δ=%g | cumul=%g\n",
              (int)i, (int)tail, (int)head, edgestate, (int)group,
              ng_before, ng_after, S_before, S_after, delta, CHANGE_STAT[0]);
    }
#endif

    /* Carry the intermediate state forward if more toggles remain. */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* 4) Restore the original network state. */
  UNDO_PREVIOUS_TOGGLES(i);
}
