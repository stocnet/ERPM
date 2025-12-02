/**
 * @file changestat_dyadcov.c
 * @brief Change statistic for the ERPM term `dyadcov` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for an ERPM effect
 *  `dyadcov`, defined on a bipartite network with:
 *    - actor mode  = actor vertices,
 *    - group mode  = group vertices.
 *
 *  Each group-mode vertex represents a structural group, and its members are
 *  the actor-mode vertices connected to it by membership edges (in either
 *  direction, if the network is stored as directed).
 *
 *  A numeric dyadic covariate matrix Z is defined on the actor mode:
 *    - dimension: n1 × n1, where n1 is the number of actors,
 *    - indexing:  column-major order as in R,
 *      Z[(j-1)*n1 + (i-1)] = z_ij for actors i,j.
 *
 *  The matrix Z is allowed to be non-symmetric. For each unordered actor pair
 *  {i,j} with i<j, the effect uses the symmetric combination:
 *
 *      z_ij + z_ji
 *
 *  whenever both entries are available in the n1 × n1 block. If Z happens to
 *  be symmetric, this reduces to 2*z_ij for each pair {i,j}.
 *
 *  For a fixed clique size k ≥ 2, and for each group g, let:
 *    - A(g) be the set of actors in group g,
 *    - C_k(g) be the set of all k-subsets C ⊂ A(g),
 *    - for a given clique C, define:
 *
 *        P(C; Z) = ∏_{i<j ∈ C} (z_ij + z_ji)
 *
 *  Then the group-level dyadic covariate functional is:
 *
 *      S_g^{(k)}(Z) = ∑_{C ∈ C_k(g)} P(C; Z).
 *
 *  Two variants of the global statistic are supported.
 *
 *  ------------------------------------------------------------
 *  Non-normalised statistic
 *  ------------------------------------------------------------
 *
 *  The non-normalised statistic is:
 *
 *      T^{(k)}(p; Z) = ∑_g S_g^{(k)}(Z)
 *                    = ∑_g ∑_{C ∈ C_k(g)} ∏_{i<j ∈ C} (z_ij + z_ji).
 *
 *  All k-cliques of actors inside each group are enumerated, and their
 *  dyadic products (based on z_ij + z_ji) are summed.
 *
 *  ------------------------------------------------------------
 *  Normalised statistic (group-size normalisation)
 *  ------------------------------------------------------------
 *
 *  Let n_g = |A(g)| be the size of group g. The normalised version now uses
 *  a per-group factor 1 / n_g (rather than 1 / C(n_g, k)):
 *
 *      T_norm^{(k)}(p; Z) =
 *        ∑_g 1[n_g ≥ k] * (1 / n_g) * S_g^{(k)}(Z).
 *
 *  In other words, each group contribution is the clique-based sum S_g^{(k)}(Z)
 *  divided by the size n_g of the group (and groups with n_g < k contribute 0
 *  because S_g^{(k)}(Z) = 0 in that case).
 *
 *  ------------------------------------------------------------
 *  Implementation outline (one-toggle)
 *  ------------------------------------------------------------
 *
 *  A bipartite network is assumed with:
 *    - actor mode  = first n1 vertices,
 *    - group mode  = remaining vertices.
 *
 *  For each membership toggle between one actor and one group:
 *
 *    1. Identify the actor vertex (actor mode) and the group vertex
 *       (group mode) involved in the toggle.
 *
 *    2. For the affected group g:
 *         - reconstruct the set A(g) of actor members from the current
 *           network state,
 *         - compute n_g and S_before = S_g^{(k)}(Z) using the helper
 *           ::group_dyadcov_k().
 *
 *    3. Apply a virtual toggle (::TOGGLE) to temporarily update the network.
 *         - reconstruct A(g) again,
 *         - compute n_g and S_after = S_g^{(k)}(Z) on the virtually
 *           updated state.
 *
 *    4. Undo the virtual toggle (::TOGGLE again) to restore the original
 *       network state.
 *
 *    5. Depending on the normalisation flag:
 *         - If not normalised:
 *
 *              Δ = S_after − S_before.
 *
 *         - If normalised (group-size normalisation):
 *              - define n_g_before and n_g_after as the actor counts before
 *                and after the toggle,
 *              - define group-level normalised values:
 *
 *                  T_before = (1[n_g_before ≥ k] / n_g_before) * S_before
 *                  T_after  = (1[n_g_after  ≥ k] / n_g_after)  * S_after
 *
 *                with the convention that T_before/T_after are 0 if
 *                n_g_before or n_g_after is 0 or < k,
 *
 *              - and the local change is:
 *
 *                  Δ = T_after − T_before.
 *
 *    6. Update the scalar statistic:
 *
 *         CHANGE_STAT[0] += Δ.
 *
 *  The change statistic is “one-toggle”: each call reports only the local
 *  change Δ associated with a single toggle. The {ergm} engine accumulates
 *  these local changes over all toggles when evaluating the statistic.
 *
 *  ------------------------------------------------------------
 *  INPUT_PARAM layout
 *  ------------------------------------------------------------
 *
 *  The R-side initialiser (InitErgmTerm.dyadcov) passes parameters as:
 *
 *      INPUT_PARAM = c(n1, k, normalized, as.vector(Z))
 *
 *  where:
 *    - n1          = number of actors (size of actor mode),
 *    - k           = clique size (integer ≥ 2),
 *    - normalized  = 0 for raw sum, 1 for per-group normalisation 1/n_g,
 *    - Z           = numeric vector of length n1*n1 in column-major order.
 *
 *  At the C level the layout is:
 *
 *      ip[0]   = n1
 *      ip[1]   = k
 *      ip[2]   = normalized
 *      ip[3+]  = Z[0 .. n1*n1-1] (column-major)
 *
 *  ------------------------------------------------------------
 *  Complexity and limitations
 *  ------------------------------------------------------------
 *
 *  For a given group of size n_g:
 *    - The number of k-cliques is C(n_g, k),
 *    - For each clique, the product involves k*(k-1)/2 unordered pairs,
 *      and for each pair {i,j} the implementation reads z_ij and z_ji
 *      and uses their sum (z_ij + z_ji).
 *    - Therefore the complexity is roughly
 *      O(C(n_g, k) * k^2) with a constant factor reflecting two accesses
 *      per pair instead of one.
 *
 *  Since the implementation recomputes S_g^{(k)}(Z) from scratch before and
 *  after each toggle, this effect is not intended for very large groups
 *  or large k (combinatorial explosion of cliques).
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  The R-side term constructor:
 *    - validates n1, k, normalized, and the dimensions of Z,
 *    - flattens Z in column-major order,
 *    - sets N_CHANGE_STATS = 1 and emptynwstats = 0,
 *    - builds INPUT_PARAM as described above.
 *
 *  ------------------------------------------------------------
 *  @example Usage (R)
 *  ------------------------------------------------------------
 *  @code{.r}
 *  library(ERPM)
 *
 *  # Example: 4 actors, 2 groups (bipartite partition representation)
 *  part <- c(1, 1, 2, 2)  # actors 1,2 in group 1; actors 3,4 in group 2
 *
 *  # Dyadic covariate matrix on the actor mode (4 x 4)
 *  set.seed(1)
 *  Z <- matrix(runif(4 * 4), nrow = 4, ncol = 4)
 *  diag(Z) <- 0
 *
 *  # Non-normalised dyadic covariate statistic on 3-cliques,
 *  # using (z_ij + z_ji) for each unordered pair {i,j}
 *  fit_raw <- erpm(partition ~ dyadcov(k = 3, cov = Z, normalized = FALSE))
 *  summary(fit_raw)
 *
 *  # Normalised version (per-group scaling by 1 / n_g on 3-cliques)
 *  fit_norm <- erpm(partition ~ dyadcov(k = 3, cov = Z, normalized = TRUE))
 *  summary(fit_norm)
 *
 *  # Internally, each sampler toggle between an actor and a group calls
 *  # c_dyadcov(), which recomputes S_g^{(k)}(Z) for the affected group
 *  # before and after the virtual toggle, and updates CHANGE_STAT[0]
 *  # by Δ = (T_after - T_before).
 *  @endcode
 *
 *  @test
 *  A self-test for this change statistic can:
 *    - build small bipartite networks from known partitions,
 *    - compute T^{(k)}(p; Z) and T_norm^{(k)}(p; Z) directly in R by:
 *         - enumerating actor sets A(g) per group,
 *         - enumerating all k-subsets C ⊂ A(g),
 *         - forming products ∏_{i<j ∈ C} (z_ij + z_ji),
 *    - compare these reference values to:
 *         - summary( erpm(partition ~ dyadcov(...)) )$statistics,
 *    - apply single-edge toggles and check that observed changes match
 *      the local Δ reported by ::c_dyadcov.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>

/**
 * @def DEBUG_DYADCOV
 * @brief Enable verbose debugging output for ::c_dyadcov.
 *
 * Set this macro to 1 to print detailed information to the R console during
 * `summary()` or MCMC runs:
 *  - current parameter values (n1, k, normalized),
 *  - group sizes before and after the virtual toggle,
 *  - raw sums S_before, S_after and the resulting Δ.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_DYADCOV 0

/**
 * @def UNUSED_WARNING
 * @brief Macro to explicitly mark a parameter as intentionally unused.
 *
 * @param x Identifier of the variable that should be marked unused.
 *
 * This macro suppresses compiler warnings about unused variables by casting
 * them to void. It is used when a parameter is required by the interface
 * but not needed in the current implementation.
 */
#define UNUSED_WARNING(x) (void)x

/* -------------------------------------------------------------------------- */
/* Clique enumeration: sum_cliques_k                                         */
/* -------------------------------------------------------------------------- */
/**
 * @brief Sum over all k-cliques of actors inside a given group.
 *
 * @details
 *  Given a vector of actor indices \c actors (1-based vertex indices in the
 *  actor mode) of length \c ng, this function:
 *
 *    - enumerates all k-subsets C of these actors,
 *    - for each subset C, computes the product:
 *
 *          P(C; Z) = ∏_{i<j ∈ C} (z_ij + z_ji),
 *
 *      using the dyadic covariate matrix Z,
 *    - accumulates the total:
 *
 *          ∑_{C} P(C; Z).
 *
 *  The matrix Z is indexed as Z[(j-1)*n1 + (i-1)] for actors i,j, assuming
 *  column-major order (R convention).
 *
 * @param actors  Pointer to an array of length \c ng containing actor vertex
 *                indices in the actor mode (1-based).
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

  /* Working array for combinations of indices into actors[]. */
  int *comb = (int*)Calloc(k, int);
  for(int i = 0; i < k; i++) comb[i] = i;

  double total = 0.0;

  while(1){
    /* Product over all unordered pairs i<j within the current clique,
     * using (z_ij + z_ji) for each pair {i,j}. */
    double prod = 1.0;
    for(int p = 0; p < k; p++){
      int idx_i = actors[ comb[p] ];   /* actor vertex index in 1..n1 */
      int row   = idx_i - 1;           /* 0..n1-1 */

      for(int q = p + 1; q < k; q++){
        int idx_j = actors[ comb[q] ];
        int col   = idx_j - 1;
        int idx_ij = col * n1 + row;   /* z_ij, column-major */
        int idx_ji = row * n1 + col;   /* z_ji, column-major */
        double zij = Z[idx_ij];
        double zji = Z[idx_ji];
        prod *= (zij + zji);
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
/* Group-level functional: group_dyadcov_k                                   */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute S_g^{(k)}(Z) for a given group vertex.
 *
 * @details
 *  For a group-mode vertex \c g, this function:
 *
 *    1. Collects all neighbouring actor-mode vertices via both outgoing and
 *       incoming edges from g (deduplicated).
 *
 *    2. Let n_g be the number of such actors; optionally stores it into
 *       \c *n_g_out if \c n_g_out is non-NULL.
 *
 *    3. If n_g < k, returns 0 (no k-cliques are possible).
 *
 *    4. Otherwise:
 *         - builds the list of actor indices (1..n1) belonging to the group,
 *         - calls ::sum_cliques_k() to compute:
 *
 *              S_g^{(k)}(Z) = ∑_{C ∈ C_k(g)} ∏_{i<j ∈ C} (z_ij + z_ji).
 *
 * @param g        Group vertex whose actor members define the group.
 * @param n1       Number of actors (dimension of the actor mode).
 * @param k        Clique size (k ≥ 2).
 * @param Z        Pointer to the dyadic covariate matrix (n1*n1, column-major).
 * @param nwp      Pointer to the {ergm} Network structure (provides edges).
 * @param n_g_out  If non-NULL, receives the number of actors in group g.
 *
 * @return S_g^{(k)}(Z) for group g.
 */
static double group_dyadcov_k(Vertex g,
                              int n1, int k,
                              const double *Z,
                              Network *nwp,
                              int *n_g_out){

  /* Temporary bitmap of actor membership for this group. */
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char); /* initialised to 0 */
  Vertex h;
  Edge e;

  /* Mark actor neighbours reachable via outgoing edges of the group. */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      seen[idx] = 1;
    }
  }

  /* Mark actor neighbours reachable via incoming edges of the group. */
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

  /* If the group is smaller than k, no k-cliques exist. */
  if(ng < k){
#if DEBUG_DYADCOV
    Rprintf("[dyadcov][group_dyadcov_k] g=%d ng=%d < k=%d -> 0\n",
            (int)g, ng, k);
#endif
    Free(seen);
    return 0.0;
  }

  /* Collect the 1-based actor indices belonging to g. */
  int *actors = (int*)Calloc(ng, int);
  int idx = 0;
  for(int i = 0; i < n1; i++){
    if(seen[i]){
      actors[idx++] = i + 1;  /* vertex index in 1..n1 */
    }
  }

  double sum = sum_cliques_k(actors, ng, k, n1, Z);

#if DEBUG_DYADCOV
  Rprintf("[dyadcov][group_dyadcov_k] g=%d ng=%d k=%d -> sum=%g\n",
          (int)g, ng, k, sum);
#endif

  Free(actors);
  Free(seen);

  return sum;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: dyadcov (one-toggle)                                    */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `dyadcov`.
 *
 * @details
 *  This is the {ergm} change-statistic function registered as ::c_dyadcov via
 *  ::C_CHANGESTAT_FN. It implements the one-toggle update for the dyadic
 *  covariate statistic on k-cliques of actors inside each group.
 *
 *  The function:
 *    - assumes a bipartite network with an actor mode and a group mode,
 *    - receives:
 *        - n1        = number of actors (actor mode),
 *        - k         = clique size (k ≥ 2),
 *        - normalized= 0 or 1,
 *        - Z         = n1 × n1 dyadic covariate matrix on the actor mode,
 *    - recomputes the group-level sums S_g^{(k)}(Z) before and after a virtual
 *      toggle of the membership edge, for the unique affected group.
 *
 *  The local change Δ is then either:
 *    - S_after − S_before for the non-normalised statistic, or
 *    - (S_after / n_g_after) − (S_before / n_g_before) for the normalised
 *      statistic, with both terms set to 0 when n_g_before or n_g_after is 0
 *      or < k.
 *
 *  The result is accumulated into the single scalar statistic:
 *
 *      CHANGE_STAT[0] += Δ.
 *
 * @param tail       Tail vertex of the toggled edge (actor or group).
 * @param head       Head vertex of the toggled edge (actor or group).
 * @param mtp        Pointer to the model term structure (unused directly
 *                   here, but required by the macro signature).
 * @param nwp        Pointer to the network-plus workspace (used by the
 *                   group_dyadcov_k() helper to inspect edges).
 * @param edgestate  Current state of the edge:
 *                   - 0 if the edge is absent (toggle = addition),
 *                   - 1 if the edge is present (toggle = deletion).
 *
 * @note
 *  - The actor mode is identified via the parameter n1 coming from
 *    INPUT_PARAM; vertices with index ≤ n1 are actors, and vertices with
 *    index > n1 are groups.
 *  - The term is non-vectorised: N_CHANGE_STATS == 1 and only
 *    CHANGE_STAT[0] is written.
 *  - The function uses virtual toggling via ::TOGGLE to obtain “before”
 *    and “after” statistics for the same group in a consistent way.
 */
C_CHANGESTAT_FN(c_dyadcov){
  /* 1) Reset the output buffer for THIS toggle.
   *
   * {ergm} accumulates contributions from multiple calls; here we only
   * report the local change Δ for the current membership toggle.
   */
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);  /* edgestate is implicit in TOGGLE, but kept for signature consistency. */

  /* 2) Decode INPUT_PARAM layout: [n1, k, normalized, Z...]. */
  const double *ip         = INPUT_PARAM;
  const int     n1         = (int)ip[0];  /* number of actors */
  int           k          = (int)ip[1];  /* clique size      */
  const int     normalized = (int)ip[2];  /* 0 = raw, 1 = normalised (1/n_g) */
  const double *Z          = ip + 3;      /* dyadic covariate matrix */

#if DEBUG_DYADCOV
  Rprintf("[dyadcov] n1=%d k=%d normalized=%d\n", n1, k, normalized);
#endif

  /* Safety: k < 2 yields a degenerate statistic (no proper k-cliques). */
  if(k < 2){
    CHANGE_STAT[0] = 0.0;
    return;
  }

  /* 3) Identify the actor and group vertices involved in the toggle.
   *
   * By construction:
   *   - exactly one endpoint has index ≤ n1 (actor),
   *   - the other endpoint has index > n1 (group).
   */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;
  Vertex group = (a <= (Vertex)n1) ? b : a;
  UNUSED_WARNING(actor);  /* actor is not needed explicitly beyond this point. */

  int ng_before = 0, ng_after = 0;

  /* 4) Group contribution BEFORE the virtual toggle. */
  double S_before = group_dyadcov_k(group, n1, k, Z, nwp, &ng_before);

  /* 5) Virtual toggle: temporarily change the membership edge. */
  TOGGLE(a, b);

  /* 6) Group contribution AFTER the virtual toggle. */
  double S_after  = group_dyadcov_k(group, n1, k, Z, nwp, &ng_after);

  /* 7) Undo the virtual toggle to restore the original network state. */
  TOGGLE(a, b);

  /* 8) Compute the local change Δ depending on the normalisation flag. */
  double delta;
  if(!normalized){
    /* Raw version: sum of clique products. */
    delta = S_after - S_before;
  } else {
    /* Normalised version: group-size normalisation 1 / n_g. */
    double T_before = 0.0;
    double T_after  = 0.0;

    if(ng_before >= k && ng_before > 0){
      T_before = S_before / (double)ng_before;
    }
    if(ng_after >= k && ng_after > 0){
      T_after = S_after / (double)ng_after;
    }

    delta = T_after - T_before;
  }

  /* 9) Accumulate the scalar change statistic. */
  CHANGE_STAT[0] += delta;

#if DEBUG_DYADCOV
  Rprintf("[dyadcov] a=%d b=%d group=%d ng_before=%d ng_after=%d "
          "S_before=%g S_after=%g delta=%g\n",
          (int)a, (int)b, (int)group,
          ng_before, ng_after, S_before, S_after, delta);
#endif
}
