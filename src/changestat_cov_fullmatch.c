/**
 * @file
 * @brief  Change statistic for the ERPM term `cov_fullmatch` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `cov_fullmatch`, which detects groups that are completely homogeneous
 *  with respect to a categorical actor covariate, with an optional filter
 *  on group sizes and an optional targeted category.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  The bipartite network encodes a partition:
 *    - actor mode  = vertices representing actors (individuals),
 *    - group mode  = vertices representing structural groups.
 *
 *  Membership is represented by edges between actors and groups.
 *  For a group vertex g in the group mode:
 *
 *    - A(g)      = set of actors adjacent to g (members of group g),
 *    - n_g       = |A(g)| = size of group g,
 *    - c_i       = integer code of the categorical covariate for actor i,
 *    - n_{g,r}   = number of actors in group g with category r.
 *
 *  Let S be an optional set of allowed group sizes.
 *
 *  Non-targeted version ("full match" on any category):
 *
 *      T(p; c)
 *        = ∑_g  1[n_g ∈ S] * 1[∃ r such that n_{g,r} = n_g],
 *
 *  i.e. each group contributes 1 if:
 *    - its size n_g is allowed by S (or S is empty), and
 *    - all actors in the group share the same category r,
 *      with no missing (NA) categories.
 *
 *  Targeted version (category κ):
 *
 *      T^{(κ)}(p; c)
 *        = ∑_g  1[n_g ∈ S] * 1[n_{g,κ} = n_g],
 *
 *  i.e. each group contributes 1 if:
 *    - its size n_g is allowed by S (or S is empty), and
 *    - all actors in the group have category κ and no NA.
 *
 *  In both variants, groups containing at least one NA category are
 *  never considered "full match".
 *
 *  ------------------------------------------------------------
 *  Bipartite structure (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  At the C level, the bipartite structure is encoded via the boundary
 *  BIPARTITE:
 *
 *    - actor mode   : vertices 1 .. n1, where n1 = BIPARTITE,
 *    - group mode   : vertices > n1, representing groups.
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
 *  INPUT_PARAM layout (from InitErgmTerm.cov_fullmatch)
 *  ------------------------------------------------------------
 *
 *  The R initialiser packs the parameters into INPUT_PARAM as:
 *
 *    INPUT_PARAM = c(
 *      n1,
 *      L,
 *      sizes[1:L],
 *      K,
 *      target,
 *      cats[1:n1]
 *    )
 *
 *  where:
 *    - n1          = number of actors (size of the actor mode),
 *    - L           = length of the size filter S,
 *    - sizes[ ]    = allowed group sizes (as doubles, cast to int),
 *    - K           = number of distinct categories encoded (max code),
 *    - target      = targeted category code (0 if non-targeted),
 *    - cats[ ]     = integer codes for categories on actors:
 *                     0 = NA / ignored, 1..K = valid categories.
 *
 *  In C, this becomes:
 *
 *    INPUT_PARAM[0]          = n1
 *    INPUT_PARAM[1]          = L
 *    INPUT_PARAM[2..1+L]     = sizes[0..L-1]
 *    INPUT_PARAM[2+L]        = K
 *    INPUT_PARAM[3+L]        = target
 *    INPUT_PARAM[4+L..]      = cats[0..n1-1]
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
 *  A single toggle flips the membership of one actor in one group:
 *    - addition  : actor becomes member of the group,
 *    - deletion  : actor leaves the group.
 *
 *  Only the affected group g can change its contribution:
 *
 *      Δ = full_flag_after(g) - full_flag_before(g),
 *
 *  where full_flag(g) ∈ {0,1} indicates whether group g satisfies:
 *    - the size filter S, and
 *    - the "full match" condition (on any category or on target κ).
 *
 *  The helper function group_flag() recomputes this flag for group g
 *  by:
 *    1. reconstructing its membership in the actor mode,
 *    2. building a frequency table of categories,
 *    3. checking the size filter and homogeneity conditions.
 *
 *  ------------------------------------------------------------
 *  Complexity
 *  ------------------------------------------------------------
 *
 *  For a single toggle and affected group g:
 *    - retrieving degrees is O(1),
 *    - building the "seen" array is O(n1) but effectively only touches
 *      neighbours of g,
 *    - neighbour traversal is O(deg(g)),
 *    - building counts over K categories is O(deg(g)) (bounded by K).
 *
 *  No per-group state is cached; the computation is local and
 *  recomputed for each toggle.
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  The R initialiser (InitErgmTerm.cov_fullmatch):
 *    - validates the actor covariate and the size filter S,
 *    - encodes categories into integer codes 0..K,
 *    - packs n1, L, sizes, K, target, and cats into INPUT_PARAM,
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
 *  # Categorical covariate on actors (e.g., a label or type)
 *  cat <- c("A", "A", "B", "B", "A", "C")
 *
 *  # Full homogeneous groups on any category, no size filter
 *  fit1 <- erpm(partition ~ cov_fullmatch(attr = cat))
 *  summary(fit1)
 *
 *  # Only groups of size 3 that are unanimously category "A"
 *  fit2 <- erpm(partition ~ cov_fullmatch(attr = cat,
 *                                         size = 3,
 *                                         category = "A"))
 *  summary(fit2)
 *
 *  # Internally, each membership toggle between an actor and a group
 *  # calls c_cov_fullmatch(), which:
 *  #   - recomputes the homogeneity flag for the affected group,
 *  #   - evaluates the size filter,
 *  #   - updates CHANGE_STAT[0] by the difference (after - before).
 *  @endcode
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"      /* Calloc/Free */
#include <R_ext/Print.h>
#include <string.h>            /* memset */

/**
 * @def DEBUG_COV_FULLMATCH
 * @brief Enable verbose debugging output for ::c_cov_fullmatch.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - group degrees and orientation,
 *  - actor category in the current toggle,
 *  - group-level flags before and after the virtual toggle.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_COV_FULLMATCH 0  /* set to 0 to disable debug output */

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
  if(L==0) return 1;
  for(int i=0; i<L; i++){
    if((int)sizes[i] == n) return 1;
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Helper: compute full-match flag for a single group                         */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute the "full match" indicator for a single group vertex.
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
 *         - cnt[r]   = count of actors with category r (for r=1..K),
 *         - has_na   = 1 if any actor has category code 0 (NA).
 *
 *    3. It applies the size filter S:
 *         - if n_g == 0 or n_g ∉ S, the result is 0.
 *
 *    4. It checks the full-match condition:
 *         - if target > 0:
 *             result = 1 if (no NA) and (cnt[target] == n_g),
 *             result = 0 otherwise,
 *         - if target == 0:
 *             result = 1 if (no NA) and (max_r cnt[r] == n_g),
 *             result = 0 otherwise.
 *
 *  The actor categories are provided as integer codes:
 *    - cats[i] = 0        → NA, excluded from matching,
 *    - cats[i] = 1..K     → valid category code.
 *
 *  The function uses an internal buffer for category counts. If K is
 *  small (≤ 1024) it uses a stack-allocated array, otherwise it allocates
 *  counts on the heap.
 *
 * @param g       Group vertex in the group mode.
 * @param n1      Number of actors (size of the actor mode).
 * @param L       Length of the size filter S.
 * @param sizes   Pointer to allowed group sizes S (as doubles).
 * @param K       Number of categories (maximum code).
 * @param target  Target category code (>0 for targeted variant, 0 otherwise).
 * @param cats    Pointer to actor category codes (length n1, 0..K).
 * @param nwp     Pointer to the network-plus workspace.
 *
 * @return 1.0 if group g is a "full match" under the given settings,
 *         0.0 otherwise.
 */
static double group_flag(Vertex g, int n1, int L, const double *sizes,
                         int K, int target, const double *cats,
                         Network *nwp){

  int ng = 0;            /* group size (number of actors in the group) */
  int has_na = 0;        /* flag indicating presence of NA categories */
  int cnt_max = 0;       /* maximum category count within this group */

  /* Decide whether to use stack or heap for category counts. */
  int use_stack_cnt = (K > 0 && K <= 1024);
  int cntK[1024];
  int *cnt = NULL;
  if(K > 0){
    cnt = use_stack_cnt ? cntK : (int*)Calloc(K, int);
    for(int i=0; i<K; i++) cnt[i] = 0;
  }

  /* ----------------------------------------------------------------------
   * Build the set of actor neighbours, deduplicated via `seen`.
   * -------------------------------------------------------------------- */
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char); /* Calloc => zeroed */
  Vertex h;
  Edge e;

  /* OUT-neighbours: group → actor edges. */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        ng++;                               /* increment group size */
        int cat = (int)cats[idx];          /* category for this actor */
        if(cat == 0){
          has_na = 1;                      /* mark presence of NA */
        }else if(cat >= 1 && cat <= K){
          int v = ++cnt[cat-1];            /* update category count */
          if(v > cnt_max) cnt_max = v;     /* track maximum count */
        }
      }
    }
  }

  /* IN-neighbours: actor → group edges. */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        ng++;                               /* increment group size */
        int cat = (int)cats[idx];          /* category for this actor */
        if(cat == 0){
          has_na = 1;                      /* mark presence of NA */
        }else if(cat >= 1 && cat <= K){
          int v = ++cnt[cat-1];            /* update category count */
          if(v > cnt_max) cnt_max = v;     /* track maximum count */
        }
      }
    }
  }

  /* Optional debug output on degrees and orientation. */
  #if DEBUG_COV_FULLMATCH
    Rprintf("[cov_fullmatch][deg] g=%d OUT_DEG=%d IN_DEG=%d DIRECTED=%d\n",
            (int)g, (int)OUT_DEG[g], (int)IN_DEG[g], (int)DIRECTED);
  #endif

  /* Empty group: never a full match. */
  if(ng == 0){
    #if DEBUG_COV_FULLMATCH
      Rprintf("[cov_fullmatch][group_flag] g=%d ng=0 -> res=0 (empty group)\n",
              (int)g);
    #endif
    if(cnt && !use_stack_cnt) Free(cnt);
    Free(seen);
    return 0.0;
  }

  /* Apply size filter S. If group size is not allowed, result is 0. */
  if(!in_sizes(ng, L, sizes)){
    #if DEBUG_COV_FULLMATCH
      Rprintf("[cov_fullmatch][group_flag] g=%d ng=%d -> excluded by size filter\n",
              (int)g, ng);
    #endif
    if(cnt && !use_stack_cnt) Free(cnt);
    Free(seen);
    return 0.0;
  }

  /* ----------------------------------------------------------------------
   * Compute final result:
   *  - targeted version: all actors must share the target category,
   *  - non-targeted version: some category must cover the entire group.
   *  Missing (NA) categories disqualify the group.
   * -------------------------------------------------------------------- */
  double res;
  if(target > 0){
    int c = has_na ? -1 : (K > 0 ? cnt[target-1] : 0);
    res = (!has_na && c == ng) ? 1.0 : 0.0;
    #if DEBUG_COV_FULLMATCH
      Rprintf("[cov_fullmatch][group_flag] g=%d ng=%d target=%d c=%d has_na=%d -> res=%.0f\n",
              (int)g, ng, target, c, has_na, res);
    #endif
  }else{
    res = (!has_na && cnt_max == ng) ? 1.0 : 0.0;
    #if DEBUG_COV_FULLMATCH
      Rprintf("[cov_fullmatch][group_flag] g=%d ng=%d cnt_max=%d has_na=%d -> res=%.0f\n",
              (int)g, ng, cnt_max, has_na, res);
    #endif
  }

  if(cnt && !use_stack_cnt) Free(cnt);

  /* Reset "seen" marks for this group and free the buffer. */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1) seen[(int)h - 1] = 0;
  }
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1) seen[(int)h - 1] = 0;
  }
  Free(seen);

  return res;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_fullmatch (one-toggle)                               */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_fullmatch`.
 *
 * @details
 *  This is the {ergm} change-statistic function registered as
 *  ::c_cov_fullmatch via ::C_CHANGESTAT_FN. It computes the local change
 *  Δ in the "full match" statistic for a single membership toggle between
 *  an actor and a group.
 *
 *  The layout of INPUT_PARAM is:
 *
 *    INPUT_PARAM[0]        = n1            (number of actors)
 *    INPUT_PARAM[1]        = L             (length of the size filter S)
 *    INPUT_PARAM[2..1+L]   = sizes[0..L-1] (allowed group sizes)
 *    INPUT_PARAM[2+L]      = K             (number of categories)
 *    INPUT_PARAM[3+L]      = target        (target category code, 0 if none)
 *    INPUT_PARAM[4+L..]    = cats[0..n1-1] (actor category codes 0..K)
 *
 *  For each toggle:
 *    1. Identify the actor vertex and the group vertex using the
 *       bipartite boundary between actor mode and group mode.
 *    2. Evaluate the group flag before the virtual toggle via group_flag().
 *    3. Apply a virtual toggle (TOGGLE) on the actor–group edge.
 *    4. Evaluate the group flag after the virtual toggle.
 *    5. Undo the virtual toggle.
 *    6. Update:
 *
 *         CHANGE_STAT[0] += full_after - full_before.
 *
 *  The parameter @p edgestate is unused here, because the function
 *  explicitly performs a virtual TOGGLE to obtain the "after" state.
 *
 * @param tail       Tail vertex of the toggled edge.
 * @param head       Head vertex of the toggled edge.
 * @param mtp        Pointer to the model term parameters (unused here,
 *                   but required by the macro signature).
 * @param nwp        Pointer to the network-plus workspace, providing
 *                   adjacency and degree information.
 * @param edgestate  Current state of the edge (unused in this implementation).
 */
C_CHANGESTAT_FN(c_cov_fullmatch){
  /* 1) Reset the output buffer for THIS toggle. */
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  /* 2) Read inputs from INPUT_PARAM. */
  const double *ip    = INPUT_PARAM;
  const int n1        = (int)ip[0];       /* number of actors (actor mode) */
  const int L         = (int)ip[1];       /* length of size filter S */
  const double *sizes = ip + 2;           /* pointer to size filter array */
  const int K         = (int)ip[2 + L];   /* number of categories */
  const int target    = (int)ip[3 + L];   /* target category code (0 if none) */
  const double *cats  = ip + 4 + L;       /* actor category codes (0..K) */

  #if DEBUG_COV_FULLMATCH
    Rprintf("[cov_fullmatch] n1=%d L=%d K=%d target=%d\n",
            n1, L, K, target);
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

  /* Optional diagnostic: presence of the edge before the virtual toggle. */
  #if DEBUG_COV_FULLMATCH
    int present_out = IS_OUTEDGE(actor, group);
    int present_in  = IS_INEDGE(group, actor);
    int actor_cat   = (int)cats[(int)actor - 1];
    Rprintf("[cov_fullmatch] actor=%d cat=%d group=%d present_out=%d present_in=%d\n",
            (int)actor, actor_cat, (int)group, present_out, present_in);
  #endif

  /* 4) Evaluate group flag before and after a virtual toggle. */
  double F_before = group_flag(group, n1, L, sizes, K, target, cats, nwp);

  /* Apply a virtual toggle (single-edge API). */
  TOGGLE(a, b);

  double F_after  = group_flag(group, n1, L, sizes, K, target, cats, nwp);

  /* Undo the virtual toggle to restore the original state. */
  TOGGLE(a, b);

  /* 5) Update change statistic: Δ = F_after - F_before. */
  CHANGE_STAT[0] += (F_after - F_before);

  #if DEBUG_COV_FULLMATCH
    Rprintf("[cov_fullmatch] a=%d b=%d group=%d before=%.0f after=%.0f delta=%.0f\n",
            (int)a, (int)b, (int)group, F_before, F_after, (F_after - F_before));
  #endif
}
