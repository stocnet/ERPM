/**
 * @file changestat_cov_fullmatch.c
 * @brief  Change statistic for the ERPM term `cov_fullmatch` (multi-toggle form).
 * 
 * @author Jérémie Chichignoud
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
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
 *  Multi-toggle (D_CHANGESTAT_FN) semantics
 *  ------------------------------------------------------------
 *
 *  In multi-toggle mode, multiple toggles can affect:
 *    - the same group (several actors moved in/out),
 *    - several groups.
 *
 *  The D_ change-statistic therefore processes toggles sequentially:
 *
 *    For each toggle i:
 *      1) identify the affected group g,
 *      2) compute F_before = group_flag(g) under the CURRENT intermediate state,
 *      3) apply the toggle to update the intermediate state,
 *      4) compute F_after  = group_flag(g) under the UPDATED intermediate state,
 *      5) accumulate Δ_i = F_after - F_before.
 *
 *  To keep degrees/adjacency consistent for subsequent toggles:
 *    - for i < ntoggles-1, we keep the toggle applied (TOGGLE_IF_MORE_TO_COME),
 *    - at the end, we undo all temporary toggles (UNDO_PREVIOUS_TOGGLES).
 *
 *  Special case: for the last toggle, TOGGLE_IF_MORE_TO_COME does nothing.
 *  We therefore apply a local virtual TOGGLE / TOGGLE to compute F_after
 *  for the last toggle only, without changing the final restored state.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"      /* R_Calloc/R_Free */
#include <R_ext/Print.h>
#include <string.h>            /* memset */

/**
 * @def DEBUG_COV_FULLMATCH
 * @brief Enable verbose debugging output for ::d_cov_fullmatch.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - group degrees and orientation,
 *  - group-level flags before and after each toggle,
 *  - MULTI-TOGGLE detection (ntoggles>1).
 *
 * When set to 0, the compiled code does not emit any debug traces.
 *
 * IMPORTANT:
 * - Do NOT leave this enabled for MLE fits unless you explicitly want spam.
 */
#define DEBUG_COV_FULLMATCH 0

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

static double group_flag(Vertex g, int n1, int L, const double *sizes,
                         int K, int target, const double *cats,
                         Network *nwp){

  int ng = 0;            /* group size (number of actors in the group) */
  int has_na = 0;        /* presence of NA categories */
  int cnt_max = 0;       /* maximum category count within this group */

  /* Stack vs heap for category counts. */
  int use_stack_cnt = (K > 0 && K <= 1024);
  int cntK[1024];
  int *cnt = NULL;
  if(K > 0){
    cnt = use_stack_cnt ? cntK : (int*)R_Calloc(K, int);
    for(int i=0; i<K; i++) cnt[i] = 0;
  }

  /* Deduplicate actor neighbours with a seen[] bitmap (actor mode only). */
  unsigned char *seen = (unsigned char*)R_Calloc(n1, unsigned char);

  Vertex h;
  Edge e;

  /* OUT-neighbours: group -> actor */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        ng++;
        int cat = (int)cats[idx];
        if(cat == 0){
          has_na = 1;
        }else if(cat >= 1 && cat <= K){
          int v = ++cnt[cat-1];
          if(v > cnt_max) cnt_max = v;
        }
      }
    }
  }

  /* IN-neighbours: actor -> group */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        ng++;
        int cat = (int)cats[idx];
        if(cat == 0){
          has_na = 1;
        }else if(cat >= 1 && cat <= K){
          int v = ++cnt[cat-1];
          if(v > cnt_max) cnt_max = v;
        }
      }
    }
  }

#if DEBUG_COV_FULLMATCH
  Rprintf("[cov_fullmatch][deg] g=%d OUT_DEG=%d IN_DEG=%d DIRECTED=%d\n",
          (int)g, (int)OUT_DEG[g], (int)IN_DEG[g], (int)DIRECTED);
#endif

  /* Empty group => never full match. */
  if(ng == 0){
#if DEBUG_COV_FULLMATCH
    Rprintf("[cov_fullmatch][group_flag] g=%d ng=0 -> res=0 (empty)\n", (int)g);
#endif
    if(cnt && !use_stack_cnt) R_Free(cnt);
    R_Free(seen);
    return 0.0;
  }

  /* Size filter S. */
  if(!in_sizes(ng, L, sizes)){
#if DEBUG_COV_FULLMATCH
    Rprintf("[cov_fullmatch][group_flag] g=%d ng=%d -> size filtered\n", (int)g, ng);
#endif
    if(cnt && !use_stack_cnt) R_Free(cnt);
    R_Free(seen);
    return 0.0;
  }

  /* Homogeneity / targeted check. */
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

  if(cnt && !use_stack_cnt) R_Free(cnt);

  /* Free buffers. */
  R_Free(seen);

  return res;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_fullmatch (multi-toggle)                              */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_fullmatch` (multi-toggle).
 *
 * @details
 *  This is the \pkg{ergm} multi-toggle change-statistic function registered as
 *  ::d_cov_fullmatch via ::D_CHANGESTAT_FN. It computes the local change
 *  Δ in the "full match" statistic for a *proposal consisting of multiple
 *  toggles* between actors and groups.
 *
 *  INPUT_PARAM layout:
 *    INPUT_PARAM[0]        = n1            (number of actors)
 *    INPUT_PARAM[1]        = L             (length of the size filter S)
 *    INPUT_PARAM[2..1+L]   = sizes[0..L-1] (allowed group sizes)
 *    INPUT_PARAM[2+L]      = K             (number of categories)
 *    INPUT_PARAM[3+L]      = target        (target category code, 0 if none)
 *    INPUT_PARAM[4+L..]    = cats[0..n1-1] (actor category codes 0..K)
 *
 *  Multi-toggle algorithm:
 *    - process toggles sequentially over the intermediate state;
 *    - for i < ntoggles-1, keep the toggle applied so subsequent toggles see
 *      the updated network (TOGGLE_IF_MORE_TO_COME);
 *    - for i == ntoggles-1, compute after-state via a local virtual TOGGLE,
 *      then undo immediately;
 *    - undo all temporary toggles at the end (UNDO_PREVIOUS_TOGGLES).
 */
D_CHANGESTAT_FN(d_cov_fullmatch){

#if DEBUG_COV_FULLMATCH
  static int seen_multi = 0;
  if(ntoggles > 1 && seen_multi < 10){
    Rprintf("[cov_fullmatch] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen_multi++;
  }
#endif

  /* Reset output buffer for the whole proposal. */
  ZERO_ALL_CHANGESTATS();

  /* Read inputs. */
  const double *ip    = INPUT_PARAM;
  const int n1        = (int)ip[0];
  const int L         = (int)ip[1];
  const double *sizes = ip + 2;
  const int K         = (int)ip[2 + L];
  const int target    = (int)ip[3 + L];
  const double *cats  = ip + 4 + L;

#if DEBUG_COV_FULLMATCH
  Rprintf("[cov_fullmatch] n1=%d L=%d K=%d target=%d\n", n1, L, K, target);
#endif

  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex a = TAIL(i);
    Vertex b = HEAD(i);

    /* Identify actor and group endpoints using n1 (actor mode boundary). */
    Vertex actor = (a <= (Vertex)n1) ? a : b;
    Vertex group = (a <= (Vertex)n1) ? b : a;

    /* Defensive: if toggle is not an actor-group edge, ignore safely. */
    if(actor > (Vertex)n1 || group <= (Vertex)n1){
#if DEBUG_COV_FULLMATCH
      Rprintf("[cov_fullmatch][WARN] non bipartite toggle: tail=%d head=%d (n1=%d)\n",
              (int)a, (int)b, n1);
#endif
      continue;
    }

#if DEBUG_COV_FULLMATCH
    int actor_cat = (int)cats[(int)actor - 1];
    int present_out = IS_OUTEDGE(actor, group);
    int present_in  = IS_INEDGE(group, actor);
    Rprintf("[cov_fullmatch] i=%d actor=%d cat=%d group=%d present_out=%d present_in=%d\n",
            i, (int)actor, actor_cat, (int)group, present_out, present_in);
#endif

    /* Before under current intermediate state. */
    double F_before = group_flag(group, n1, L, sizes, K, target, cats, nwp);

    double F_after = 0.0;

    if(i < (int)ntoggles - 1){
      /* Apply toggle and KEEP it for subsequent toggles. */
      TOGGLE_IF_MORE_TO_COME(i);
      F_after = group_flag(group, n1, L, sizes, K, target, cats, nwp);
    }else{
      /* Last toggle: TOGGLE_IF_MORE_TO_COME would do nothing, so do local virtual toggle. */
      TOGGLE(a, b);
      F_after = group_flag(group, n1, L, sizes, K, target, cats, nwp);
      TOGGLE(a, b);
    }

    CHANGE_STAT[0] += (F_after - F_before);

#if DEBUG_COV_FULLMATCH
    Rprintf("[cov_fullmatch] i=%d group=%d before=%.0f after=%.0f delta=%.0f cumul=%.0f\n",
            i, (int)group, F_before, F_after, (F_after - F_before), CHANGE_STAT[0]);
#endif
  }

  /* Undo all temporary toggles applied by TOGGLE_IF_MORE_TO_COME. */
  UNDO_PREVIOUS_TOGGLES(i);
}
