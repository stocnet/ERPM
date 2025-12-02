/**
 * @file changestat_cliques.c
 * @brief  Change statistic for the ERPM term `cliques(k, normalized)`.
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `cliques(k, normalized)`. The effect counts, for each group in the
 *  group mode, the number of actor k-cliques that are induced by group
 *  memberships in the actor mode.
 *
 *  Actor / group setup:
 *  - The network is bipartite:
 *    - the actor mode contains all actor vertices,
 *    - the group mode contains all group vertices.
 *  - Each actor belongs to exactly one group; edges encode membership
 *    between actor-mode and group-mode vertices.
 *
 *  Let n_g be the size of a group g (degree of the group vertex in
 *  the bipartite representation).
 *
 *  Statistic at the group level:
 *
 *  - For k >= 2:
 *
 *      cliques_k = sum_g C(n_g, k)
 *
 *    where C(n, k) is the binomial coefficient. Each group g contributes
 *    the number of k-size subsets of its actors, which correspond to
 *    k-cliques in the actor-mode projection.
 *
 *  - For k == 1:
 *
 *      cliques_1 = #{ g : n_g == 1 }
 *
 *    i.e. the number of groups of size exactly 1 (singleton groups).
 *
 *  Change statistic for a single membership toggle:
 *
 *  - For k >= 2 and a toggle that changes the size of a single group from
 *    n_g to n_g ± 1:
 *      - addition (is_add = 1):
 *          Δ =  C(n_g,   k-1)
 *      - deletion (is_add = 0):
 *          Δ = -C(n_g-1, k-1)
 *
 *  - For k == 1, using n_g = group size before the toggle:
 *      - addition:
 *          n_g == 0  →  Δ = +1
 *          n_g == 1  →  Δ = -1
 *          else      →  Δ =  0
 *      - deletion:
 *          n_g == 2  →  Δ = +1
 *          n_g == 1  →  Δ = -1
 *          else      →  Δ =  0
 *
 *  Scaling:
 *  - Each component j uses a scalar `scale_j` provided via INPUT_PARAM.
 *  - Typical choices:
 *      scale_j = 1
 *      scale_j = 1 / C(N_actors, k)   (global normalisation, valid also for k=1)
 *  - The final change for component j is:
 *
 *      CHANGE_STAT[j] += scale_j * Δ_j
 *
 *  Implementation outline:
 *  - The actor mode occupies vertex indices 1..BIPARTITE.
 *  - The group mode occupies vertex indices BIPARTITE+1..N_NODES.
 *  - For each toggle between an actor and a group:
 *      1) Identify the group-mode vertex.
 *      2) Read its current size n_g from OUT_DEG + IN_DEG.
 *      3) Determine whether the toggle is an addition or a deletion
 *         from `edgestate`.
 *      4) For each (k_j, scale_j) pair in INPUT_PARAM, compute Δ_j as
 *         above and accumulate it into CHANGE_STAT[j].
 *
 *  Complexity:
 *  - O(#stats) per toggle:
 *      - the group size is read once,
 *      - for k >= 2 we use CHOOSE(n_g, k-1) or CHOOSE(n_g-1, k-1),
 *      - no traversal of adjacency lists and no explicit clique enumeration.
 *
 *  R interface:
 *  - The corresponding R initialiser (InitErgmTerm.cliques) prepares:
 *      - a vector of (k_j, scale_j) pairs in INPUT_PARAM,
 *      - one statistic per requested k and normalisation,
 *      - emptynwstats = 0.
 *
 *  @note
 *  The design is strictly local: each toggle affects at most one group in the
 *  group mode, so only its size is needed to compute the change.
 *
 *  @example Usage (R)
 *  @code{.r}
 *  library(ERPM)
 *
 *  # Partition: 6 actors into 3 groups
 *  #   group 1: actors 1, 2
 *  #   group 2: actors 3, 4
 *  #   group 3: actors 5, 6
 *  part <- c(1, 1, 2, 2, 3, 3)
 *
 *  # Count 2-cliques (pairs of actors in the same group) without normalisation
 *  fit <- erpm(partition ~ cliques(k = 2, normalized = FALSE))
 *  summary(fit)
 *
 *  # Interpretation of a single toggle:
 *  # - Suppose a sampler proposes to add an actor to a group of size n_g = 2
 *  #   with k = 2.
 *  # - The change in the 2-clique count is:
 *  #       Δ = C(2, 1) = 2
 *  #   and this value is multiplied by the associated scale factor.
 *  @endcode
 *
 *  @test
 *  A self-test script should:
 *  - construct small bipartite networks from known partitions,
 *  - evaluate the statistic by:
 *      (i) calling `summary` on the ERGM with `cliques(k, normalized)` and
 *      (ii) computing the reference value from the group sizes `n_g`,
 *  - verify that single-edge toggles produce Δ consistent with the formulas:
 *      - k >= 2:  +C(n_g, k-1) for additions, -C(n_g-1, k-1) for deletions,
 *      - k == 1:  piecewise rules on n_g.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>

/**
 * @def DEBUG_CLIQUES
 * @brief Enable or disable verbose debugging output for ::c_cliques.
 *
 * When set to 1, the change-statistic function prints diagnostic information
 * to the R console for each toggle:
 *  - the group vertex index and its degree before the toggle,
 *  - the current k value and whether the toggle is an addition or deletion,
 *  - the raw Δ before scaling and the final accumulated CHANGE_STAT[j].
 *
 * When set to 0, no debug output is produced.
 */
#define DEBUG_CLIQUES 0

/* -------------------------------------------------------------------------- */
/* Change statistic: cliques(k, normalized)                                   */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cliques(k, normalized)`.
 *
 * @details
 *  This function is registered via ::C_CHANGESTAT_FN as ::c_cliques and is
 *  called once per proposed toggle between an actor and a group.
 *
 *  Actor / group modes:
 *  - Actor-mode vertices use indices 1..BIPARTITE.
 *  - Group-mode vertices use indices BIPARTITE+1..N_NODES.
 *
 *  For each toggle:
 *  - The function determines which endpoint is in the group mode.
 *  - It reads the group size before the toggle from OUT_DEG + IN_DEG.
 *  - It determines whether the toggle is an addition or a deletion from
 *    the `edgestate` flag.
 *  - For each statistic j, it reads:
 *      - k_j     = INPUT_PARAM[2*j + 0] (target clique size),
 *      - scale_j = INPUT_PARAM[2*j + 1] (scaling factor),
 *    and computes the local change Δ_j using the un-toggle formulas:
 *
 *      * k_j == 1:
 *          - addition:
 *              n_g == 0  →  Δ = +1
 *              n_g == 1  →  Δ = -1
 *              else      →  Δ =  0
 *          - deletion:
 *              n_g == 2  →  Δ = +1
 *              n_g == 1  →  Δ = -1
 *              else      →  Δ =  0
 *
 *      * k_j >= 2:
 *          - addition:
 *              if n_g >= k_j - 1:
 *                 Δ =  C(n_g,   k_j - 1)
 *          - deletion:
 *              if n_g - 1 >= k_j - 1 and n_g >= 1:
 *                 Δ = -C(n_g-1, k_j - 1)
 *              else:
 *                 Δ = 0
 *
 *    It then multiplies Δ_j by scale_j and accumulates the result in
 *    CHANGE_STAT[j].
 *
 * @note
 *  - The macro ::C_CHANGESTAT_FN injects the full signature, including
 *    access to CHANGE_STAT, INPUT_PARAM and N_CHANGE_STATS.
 *  - The term is designed for a bipartite representation of the partition
 *    where each actor belongs to exactly one group.
 */
C_CHANGESTAT_FN(c_cliques){
  /* 1) Reset the output buffer for this toggle.
   *
   * {ergm} accumulates the contributions of all toggles outside this
   * function, so we explicitly zero the vector here.
   */
  ZERO_ALL_CHANGESTATS(0);

  /* 2) Number of actor-mode vertices.
   *
   * In a bipartite ERGM, BIPARTITE holds the number of vertices in the
   * actor mode. All vertices with index > BIPARTITE are in the group mode.
   */
  const int n1 = BIPARTITE; // number of actors (actor mode)

  /* 3) Identify the group-mode vertex.
   *
   * The toggle always connects exactly one actor and one group.
   * Among (tail, head), the vertex with index > n1 is the group-side vertex.
   */
  Vertex t = tail, h = head;
  Vertex v2 = (t > (Vertex)n1) ? t : h;  // group-mode vertex

  /* 4) Group size before the toggle.
   *
   * OUT_DEG and IN_DEG are provided by {ergm}. For a bipartite membership
   * representation, the group size is the sum of outgoing and incoming
   * degrees of the group vertex.
   */
  int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);

  /* 5) Toggle type: addition vs deletion.
   *
   * edgestate == 0 → edge currently absent → we are adding a membership.
   * edgestate == 1 → edge currently present → we are removing a membership.
   */
  const int is_add = edgestate ? 0 : 1;

  /* 6) Loop over all statistics (possibly several k, scale pairs).
   *
   * INPUT_PARAM is organised as [k_0, scale_0, k_1, scale_1, ...].
   */
  for(int j = 0; j < N_CHANGE_STATS; ++j){
    int    k     = (int)   INPUT_PARAM[2*j + 0];  // target clique size
    double scale = (double)INPUT_PARAM[2*j + 1];  // scaling factor

    double delta = 0.0;  // raw change for this (k, scale) before scaling

    if(k == 1){
      /* Special case k == 1:
       *
       * We are counting groups of size exactly 1. Only transitions that
       * cross the threshold n_g ∈ {0,1,2} contribute.
       */
      if(is_add){
        // Addition: 0 -> 1 : +1 ; 1 -> 2 : -1 ; otherwise 0
        if(deg_old == 0)      delta = +1.0;
        else if(deg_old == 1) delta = -1.0;
      }else{
        // Deletion: 2 -> 1 : +1 ; 1 -> 0 : -1 ; otherwise 0
        if(deg_old == 2)      delta = +1.0;
        else if(deg_old == 1) delta = -1.0;
      }
    }else{
      /* Case k >= 2:
       *
       * We use the combinatorial identities for the number of k-cliques
       * induced by a group of size n_g:
       *
       *   cliques_k(g) = C(n_g, k).
       *
       * Adding or removing one actor in that group changes the count
       * according to:
       *
       *   addition: Δ = C(n_g,   k-1)
       *   deletion: Δ = -C(n_g-1, k-1)
       */
      if(is_add){
        if(deg_old >= k-1)
          delta = CHOOSE(deg_old, k-1);
      }else{
        if(deg_old-1 >= k-1 && deg_old >= 1)
          delta = -CHOOSE(deg_old - 1, k-1);
        else
          delta = 0.0;
      }
    }

    // Apply scaling and accumulate in the j-th statistic.
    CHANGE_STAT[j] += scale * delta;

    #if DEBUG_CLIQUES
      Rprintf("[c_cliques] group=%d, deg_old=%d, k=%d, add=%d, "
              "delta_raw=%.6f, scale=%g, stat_j=%.6f\n",
              (int)v2, deg_old, k, is_add, delta, scale, CHANGE_STAT[j]);
    #endif
  }
}
