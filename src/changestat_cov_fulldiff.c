/**
 * @file
 * @brief  Change statistic for the ERPM term `cov_fulldiff` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `cov_fulldiff`, which measures, for each group in the group mode, the
 *  within-group dispersion of a numeric actor covariate via the range:
 *
 *      range_g(x) = x_g^max - x_g^min,
 *
 *  with an optional filter on group sizes.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  The bipartite network represents a partition:
 *    - actor mode  = vertices representing actors (individuals),
 *    - group mode  = vertices representing structural groups.
 *
 *  Membership is encoded by edges between actors and groups.
 *  For a group vertex g in the group mode:
 *
 *    - A(g)       = set of actors adjacent to g (members of group g),
 *    - n_g        = |A(g)| = size of group g,
 *    - x_i        = numeric covariate value for actor i,
 *    - x_g^min    = min_{i in A(g)} x_i,
 *    - x_g^max    = max_{i in A(g)} x_i.
 *
 *  Let S be an optional set of allowed group sizes. The global statistic is:
 *
 *      T(p; x, S)
 *        = ∑_g 1[n_g ∈ S] * (x_g^max - x_g^min),
 *
 *  where:
 *    - p describes the bipartite membership pattern (partition),
 *    - the indicator 1[n_g ∈ S] is 1 if S is empty or n_g is in S,
 *      and 0 otherwise.
 *
 *  Groups of size 0 or 1 do not contribute because their internal dispersion
 *  is undefined or zero; this is handled explicitly in the implementation.
 *
 *  ------------------------------------------------------------
 *  Bipartite structure (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  At the C level, the bipartite structure is encoded via the boundary
 *  BIPARTITE:
 *
 *    - actor mode  : vertices 1 .. n1, where n1 = BIPARTITE,
 *    - group mode  : vertices > n1, representing groups.
 *
 *  This change statistic assumes:
 *    - the number of actors n1 (actor mode) is stored in INPUT_PARAM[0],
 *    - each membership toggle connects exactly one actor (vertex ≤ n1)
 *      and one group (vertex > n1).
 *
 *  The group vertex is detected as the endpoint in the group mode
 *  (vertex index > n1), and the actor vertex as the endpoint in
 *  the actor mode (vertex index ≤ n1).
 *
 *  ------------------------------------------------------------
 *  INPUT_PARAM layout (from InitErgmTerm.cov_fulldiff)
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
 *    - n1         = number of actors (size of the actor mode),
 *    - L          = length of the size filter S,
 *    - sizes[ ]   = allowed group sizes (stored as doubles and cast to int),
 *    - x[ ]       = numeric covariate on actors, length n1.
 *
 *  In C, this becomes:
 *
 *    INPUT_PARAM[0]        = n1
 *    INPUT_PARAM[1]        = L
 *    INPUT_PARAM[2..1+L]   = sizes[0..L-1]
 *    INPUT_PARAM[2+L..]    = x[0..n1-1]
 *
 *  The term returns a single scalar statistic:
 *
 *    - N_CHANGE_STATS = 1,
 *    - CHANGE_STAT[0] is updated by the local Δ at each toggle.
 *
 *  ------------------------------------------------------------
 *  Local change under a toggle
 *  ------------------------------------------------------------
 *
 *  A single toggle flips membership of one actor in one group:
 *    - addition  : actor becomes member of the group,
 *    - deletion  : actor leaves the group.
 *
 *  Only the affected group g can change its contribution:
 *
 *      Δ = range_after(g) * 1[n_g_after ∈ S]
 *        - range_before(g) * 1[n_g_before ∈ S],
 *
 *  where:
 *    - range(g)   = x_g^max - x_g^min,
 *    - n_g_before = group size before the toggle,
 *    - n_g_after  = group size after the toggle.
 *
 *  The helper function group_range() recomputes the contribution for
 *  group g by:
 *    1. reconstructing its membership in the actor mode (deduplicated),
 *    2. computing n_g, x_g^min and x_g^max,
 *    3. applying the size filter S and returning the contribution
 *       (range or 0).
 *
 *  ------------------------------------------------------------
 *  Complexity
 *  ------------------------------------------------------------
 *
 *  For a single toggle and affected group g:
 *    - neighbour traversal is O(deg(g)), including deduplication,
 *    - computing min and max is O(deg(g)),
 *    - the size filter check is O(L).
 *
 *  No per-group state is cached; the computation is local and
 *  recomputed for each toggle.
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  The R initialiser (InitErgmTerm.cov_fulldiff):
 *    - validates the numeric covariate on the actor mode,
 *    - sets up the size filter S,
 *    - packs n1, L, sizes and x into INPUT_PARAM,
 *    - sets emptynwstats = 0 and a single coef.name.
 *
 *  ------------------------------------------------------------
 *  @example Usage (R)
 *  ------------------------------------------------------------
 *  @code{.r}
 *  library(ERPM)
 *
 *  # Example partition: 6 actors into 3 groups
 *  part <- c(1, 1, 2, 2, 3, 3)
 *
 *  # Numeric covariate on actors (e.g., a score or continuous attribute)
 *  x <- c(1.2, 1.5, 0.3, 0.7, 2.0, 1.8)
 *
 *  # Range within each group, no size filter
 *  fit1 <- erpm(partition ~ cov_fulldiff(attr = x))
 *  summary(fit1)
 *
 *  # Only groups of size 3 contribute to the statistic
 *  fit2 <- erpm(partition ~ cov_fulldiff(attr = x, size = 3))
 *  summary(fit2)
 *
 *  # Internally, each membership toggle between an actor and a group
 *  # calls c_cov_fulldiff(), which:
 *  #   - recomputes the range-based contribution of the affected group
 *  #     before and after a virtual toggle,
 *  #   - applies the size filter,
 *  #   - updates CHANGE_STAT[0] by the difference.
 *  @endcode
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"      /* Calloc/Free */
#include <R_ext/Print.h>

/**
 * @def DEBUG_COV_FULLDIFF
 * @brief Enable verbose debugging output for ::c_cov_fulldiff.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - group sizes and ranges before and after the toggle,
 *  - inclusion or exclusion by size filter,
 *  - local contribution of the affected group.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_COV_FULLDIFF 0   /* set to 0 to disable debug output */

/**
 * @def UNUSED_WARNING
 * @brief Utility macro to explicitly mark unused parameters.
 *
 * @param x Parameter or variable that is intentionally unused in a
 *          particular compilation unit or function.
 */
#define UNUSED_WARNING(x) (void)x

/* -------------------------------------------------------------------------- */
/* Helper: group size filter                                                  */
/* -------------------------------------------------------------------------- */

/**
 * @brief Check whether a group size belongs to the allowed size set S.
 *
 * @details
 *  The size filter S is encoded as an array of doubles @p sizes of length @p L.
 *  Each entry is cast to int and compared to @p n.
 *
 *  Special case:
 *    - If L == 0, the filter is considered inactive and all sizes are
 *      accepted (the function returns 1).
 *
 * @param n      Group size to be checked.
 * @param L      Number of allowed sizes in @p sizes.
 * @param sizes  Pointer to the array of allowed sizes (stored as doubles).
 *
 * @return 1 if @p n is in S or if L == 0; 0 otherwise.
 */
static inline int in_sizes(int n, int L, const double *sizes){
  if(L == 0) return 1;
  for(int i = 0; i < L; i++){
    if((int)sizes[i] == n) return 1;
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Helper: compute range contribution for a single group                      */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute the covariate range contribution for a single group vertex.
 *
 * @details
 *  This function evaluates the contribution of one group in the group mode:
 *
 *    1. It reconstructs the set of actor neighbours A(g) by:
 *         - traversing outgoing edges (group → actor),
 *         - traversing incoming edges (actor → group),
 *         - deduplicating actors using a "seen" array.
 *
 *    2. It computes:
 *         - n_g      = number of actors in A(g),
 *         - x_g^min  = minimum covariate value over actors in A(g),
 *         - x_g^max  = maximum covariate value over actors in A(g).
 *
 *    3. It applies the size filter S:
 *         - if n_g <= 1, the contribution is 0 (no internal dispersion),
 *         - if n_g ∉ S, the contribution is 0,
 *         - otherwise the contribution is x_g^max - x_g^min.
 *
 *  The actor covariate is provided as an array:
 *    - x[i] stores the numeric covariate for actor i+1,
 *      where actors are vertices 1..n1 in the actor mode.
 *
 * @param g       Group vertex in the group mode.
 * @param n1      Number of actors (size of the actor mode).
 * @param L       Length of the size filter S.
 * @param sizes   Pointer to allowed group sizes S (as doubles).
 * @param x       Pointer to actor covariate values (length n1).
 * @param nwp     Pointer to the network-plus workspace.
 *
 * @return The contribution of group g:
 *         (x_g^max - x_g^min) if n_g>1 and n_g ∈ S, 0 otherwise.
 */
static double group_range(Vertex g,
                          int n1, int L, const double *sizes,
                          const double *x,
                          Network *nwp){

  int ng = 0;                /* number of actors in the group */
  double xmin = 0.0;         /* minimum covariate value in the group */
  double xmax = 0.0;         /* maximum covariate value in the group */
  unsigned char first = 1;   /* flag to initialise xmin/xmax on first actor */

  /* "seen" marks actors that have already been counted, to avoid
   * double counting in directed or reciprocated networks.
   * Allocated and zero-initialised by Calloc. */
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char);

  Vertex h;
  Edge e;

  /* OUT-neighbours: group → actor edges. */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;      /* actor index in 0..n1-1 */
      if(!seen[idx]){
        seen[idx] = 1;           /* mark actor as seen */
        ng++;                    /* increment group size */
        double val = x[idx];     /* covariate value for this actor */
        if(first){
          xmin = xmax = val;     /* initialisation on first actor */
          first = 0;
        }else{
          if(val < xmin) xmin = val;
          if(val > xmax) xmax = val;
        }
      }
    }
  }

  /* IN-neighbours: actor → group edges. */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;      /* actor index in 0..n1-1 */
      if(!seen[idx]){
        seen[idx] = 1;           /* mark actor as seen */
        ng++;                    /* increment group size */
        double val = x[idx];     /* covariate value for this actor */
        if(first){
          xmin = xmax = val;     /* initialisation on first actor */
          first = 0;
        }else{
          if(val < xmin) xmin = val;
          if(val > xmax) xmax = val;
        }
      }
    }
  }

  /* Optional debug output before freeing "seen". */
  #if DEBUG_COV_FULLDIFF
    Rprintf("[cov_fulldiff][group_range] g=%d ng=%d xmin=%g xmax=%g\n",
            (int)g, ng, xmin, xmax);
  #endif

  /* Reset seen marks and free the buffer. */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1) seen[(int)h - 1] = 0;
  }
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1) seen[(int)h - 1] = 0;
  }
  Free(seen);

  /* Empty group or singleton: no internal dispersion. */
  if(ng <= 1){
    #if DEBUG_COV_FULLDIFF
      Rprintf("[cov_fulldiff][group_range] g=%d ng=%d -> 0 (empty or singleton)\n",
              (int)g, ng);
    #endif
    return 0.0;
  }

  /* Apply size filter S. If group size is not allowed, contribution is 0. */
  if(!in_sizes(ng, L, sizes)){
    #if DEBUG_COV_FULLDIFF
      Rprintf("[cov_fulldiff][group_range] g=%d ng=%d -> excluded by size filter\n",
              (int)g, ng);
    #endif
    return 0.0;
  }

  /* Compute final contribution range = xmax - xmin. */
  double res = xmax - xmin;

  #if DEBUG_COV_FULLDIFF
    Rprintf("[cov_fulldiff][group_range] g=%d ng=%d xmin=%g xmax=%g -> res=%g\n",
            (int)g, ng, xmin, xmax, res);
  #endif

  return res;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_fulldiff (one-toggle)                                */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_fulldiff`.
 *
 * @details
 *  This is the {ergm} change-statistic function registered as
 *  ::c_cov_fulldiff via ::C_CHANGESTAT_FN. It computes the local change
 *  Δ in the range-based dispersion statistic for a single membership
 *  toggle between an actor and a group.
 *
 *  The layout of INPUT_PARAM is:
 *
 *    INPUT_PARAM[0]        = n1            (number of actors)
 *    INPUT_PARAM[1]        = L             (length of the size filter S)
 *    INPUT_PARAM[2..1+L]   = sizes[0..L-1] (allowed group sizes)
 *    INPUT_PARAM[2+L..]    = x[0..n1-1]    (numeric covariate on actors)
 *
 *  For each toggle:
 *    1. Identify the actor vertex and the group vertex using the
 *       bipartite boundary between actor mode and group mode.
 *    2. Evaluate the group contribution before the virtual toggle
 *       via group_range().
 *    3. Apply a virtual toggle (TOGGLE) on the actor–group edge.
 *    4. Evaluate the group contribution after the virtual toggle.
 *    5. Undo the virtual toggle.
 *    6. Update:
 *
 *         CHANGE_STAT[0] += F_after - F_before.
 *
 *  The parameter @p edgestate is unused here, because the function
 *  explicitly performs a virtual TOGGLE to obtain the "after" state.
 *
 * @param tail       Tail vertex of the toggled edge.
 * @param head       Head vertex of the toggled edge.
 * @param mtp        Pointer to the model term parameters (unused here,
 *                   but required by the macro signature).
 * @param nwp        Pointer to the current network-plus workspace, providing
 *                   adjacency and degree information.
 * @param edgestate  Current state of the edge (unused in this implementation).
 */
C_CHANGESTAT_FN(c_cov_fulldiff){
  /* 1) Reset the output buffer for THIS toggle. */
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  /* 2) Read inputs from INPUT_PARAM. */
  const double *ip    = INPUT_PARAM;
  const int n1        = (int)ip[0];       /* number of actors (actor mode) */
  const int L         = (int)ip[1];       /* length of size filter S */
  const double *sizes = ip + 2;           /* pointer to size filter array */
  const double *x     = ip + 2 + L;       /* pointer to actor covariates */

  #if DEBUG_COV_FULLDIFF
    Rprintf("[cov_fulldiff] n1=%d L=%d\n", n1, L);
  #endif

  /* 3) Identify actor and group vertices for the current toggle.
   *
   * Actor vertices are 1..n1 (actor mode).
   * Group vertices are > n1 (group mode).
   */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;  /* actor vertex in the actor mode */
  Vertex group = (a <= (Vertex)n1) ? b : a;  /* group vertex in the group mode */
  UNUSED_WARNING(actor);

  /* 4) Evaluate group contribution before and after a virtual toggle. */
  double F_before = group_range(group, n1, L, sizes, x, nwp);

  /* Apply a virtual toggle (single-edge API). */
  TOGGLE(a, b);

  double F_after  = group_range(group, n1, L, sizes, x, nwp);

  /* Undo the virtual toggle to restore the original state. */
  TOGGLE(a, b);

  /* 5) Update change statistic: Δ = F_after - F_before. */
  CHANGE_STAT[0] += (F_after - F_before);

  #if DEBUG_COV_FULLDIFF
    Rprintf("[cov_fulldiff] a=%d b=%d group=%d before=%g after=%g delta=%g\n",
            (int)a, (int)b, (int)group, F_before, F_after, (F_after - F_before));
  #endif
}
