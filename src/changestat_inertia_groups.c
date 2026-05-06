/**
 * @file changestat_inertia_groups.c
 * @brief Change statistic for the ERPM inertial term 'inertia_groups'.
 * 
 * @author Jérémie Chichignoud
 *
 * @details
 *  This change statistic operates on a PLE (stacked block-diagonal) bipartite
 *  meta-network built by erpm_long().
 *
 *  What is counted:
 *    Among current groups (mode-2 vertices), count those whose exact actor
 *    membership set (restricted to the first n_eff actors of that block)
 *    matches a group observed in the past over a window of d = past_influence
 *    lags, optionally restricted by a size filter.
 *
 *  Depth semantics:
 *    A group contributes 1 iff its signature appears identically in AT LEAST
 *    ONE past partition across lags 1..d (union over lags, i.e. Sigma^(d)_t).
 *
 *  Variable block sizes:
 *    Blocks may have different numbers of actors. Actors are identified by
 *    position (actor i in block b corresponds to actor i in any other block).
 *    When a block's current size differs from a past partition size, only the
 *    first n_eff[b] = min(current size, all past sizes for that block) actors
 *    are considered. A warning is emitted by InitErgmTerm.inertia_groups.
 *
 *  Locality:
 *    This change statistic is local to a single bipartite edge toggle. Only
 *    the affected group can change its persistence flag for that toggle:
 *
 *      delta = flag_after(group) - flag_before(group)
 *
 * INPUT_PARAM layout (packed in R):
 *
 *    ip[0]                          n1_total
 *    ip[1]                          B        (number of blocks)
 *    ip[2]                          d        (past_influence)
 *    ip[3]                          L        (size filter count)
 *    ip[4 .. 3+L]                   sizes[0..L-1]
 *    ip[4+L .. 3+L+B]               actor_offsets[0..B-1]   (0-based)
 *    ip[4+L+B .. 3+L+2B]            n_eff[0..B-1]
 *    ip[4+L+2B .. 3+L+2B+n1_total]  group_to_block[0..n1_total-1]  (1..B)
 *    ip[4+L+2B+n1_total .. 3+L+2B+n1_total+B*d]  offsets[0..B*d-1]  (0-based)
 *    data blocks (block-major, lag-major):
 *      for each (b, lag): M, len_1, ids_1..., len_2, ids_2..., ...
 *      (ids are GLOBAL actor vertex ids 1..n1_total, truncated to n_eff[b])
 *
 *  Notes:
 *    - All values are stored as doubles and cast to int in C.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>
#include <string.h>

#define DEBUG_INERTIA_GROUPS 0
#define UNUSED_WARNING(x) (void)(x)

/**
 * @brief Check whether a size passes the optional filter.
 *
 * @param n Current group size.
 * @param L Length of the size filter vector.
 * @param sizes Array of allowed sizes.
 *
 * @return 1 if accepted, 0 otherwise.
 */
static inline int in_sizes(int n, int L, const int *sizes){
  if(L == 0) return 1;
  for(int i = 0; i < L; i++){
    if(sizes[i] == n) return 1;
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Helper: build sorted member list of a group within its block               */
/* -------------------------------------------------------------------------- */
/*
 * Collect neighbors of group vertex g that fall in the actor range
 * [actor_off+1 .. actor_off+n_eff] (global ids). The result is sorted
 * because we scan the seen[] array in ascending order.
 *
 * out_ids : int array of length >= n_eff (allocated by caller).
 * Returns : group size k (number of members found).
 */
static int collect_group_members_in_block(Vertex g, int actor_off, int n_eff,
                                           Network *nwp, int *out_ids){
  unsigned char *seen = (unsigned char*)R_Calloc(n_eff, unsigned char);
  Edge e; Vertex h;

  STEP_THROUGH_OUTEDGES(g, e, h){
    int idx = (int)h - actor_off - 1;   /* 0-based within block */
    if(idx >= 0 && idx < n_eff) seen[idx] = 1;
  }
  STEP_THROUGH_INEDGES(g, e, h){
    int idx = (int)h - actor_off - 1;
    if(idx >= 0 && idx < n_eff) seen[idx] = 1;
  }

  int k = 0;
  for(int i = 0; i < n_eff; i++){
    if(seen[i]) out_ids[k++] = actor_off + 1 + i;  /* global id (1-based) */
  }

  R_Free(seen);

#if DEBUG_INERTIA_GROUPS
  Rprintf("  [collect] g=%d actor_off=%d n_eff=%d -> k=%d ids=[",
          (int)g, actor_off, n_eff, k);
  for(int i = 0; i < k; i++) Rprintf("%d%s", out_ids[i], i < k-1 ? "," : "");
  Rprintf("]\n");
#endif

  return k;
}

/* -------------------------------------------------------------------------- */
/* Helper: check membership in one (b, lag) observed partition data block    */
/* -------------------------------------------------------------------------- */
/*
 * Returns 1 iff cur_ids (length cur_n, sorted global) matches at least one
 * observed group stored in the (b, lag) data block.
 *
 * base_offsets : 0-based index in ip of the first offset entry.
 */
static int matches_one_lag_block(const int *cur_ids, int cur_n,
                                  const double *ip,
                                  int b, int lag, int d,
                                  int base_offsets){
  int idx = (b - 1) * d + (lag - 1);           /* 0-based in offsets table */
  int pos = (int)ip[base_offsets + idx];        /* 0-based position in ip   */

#if DEBUG_INERTIA_GROUPS
  Rprintf("  [lag_match] b=%d lag=%d idx=%d pos=%d cur_n=%d cur=[",
          b, lag, idx, pos, cur_n);
  for(int i = 0; i < cur_n; i++) Rprintf("%d%s", cur_ids[i], i < cur_n-1 ? "," : "");
  Rprintf("]\n");
#endif

  if(pos < 0) return 0;

  int M = (int)ip[pos++];

#if DEBUG_INERTIA_GROUPS
  Rprintf("  [lag_match] M=%d past groups to scan\n", M);
#endif

  for(int j = 0; j < M; j++){
    int len = (int)ip[pos++];

#if DEBUG_INERTIA_GROUPS
    Rprintf("  [lag_match]   past_group[%d] len=%d ids=[", j, len);
    for(int u = 0; u < len; u++) Rprintf("%d%s", (int)ip[pos+u], u < len-1 ? "," : "");
    Rprintf("]\n");
#endif

    if(len == cur_n){
      int ok = 1;
      for(int u = 0; u < len; u++){
        if((int)ip[pos + u] != cur_ids[u]){ ok = 0; break; }
      }
      if(ok){
#if DEBUG_INERTIA_GROUPS
        Rprintf("  [lag_match]   -> MATCH at past_group[%d]\n", j);
#endif
        return 1;
      }
    }
    pos += len;
  }

#if DEBUG_INERTIA_GROUPS
  Rprintf("  [lag_match]   -> NO MATCH\n");
#endif
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Helper: compare current group to past groups for block b, ANY lag in 1..d */
/* -------------------------------------------------------------------------- */
/* Returns 1 iff cur_ids matches at least one group in at least one lag      */
/* (i.e. the signature belongs to the union Sigma^(d)_t = U_{l=1}^d Sigma_{t-l}) */
static int matches_any_past_groups_block(const int *cur_ids, int cur_n,
                                          const double *ip,
                                          int b, int d,
                                          int base_offsets){
  for(int lag = 1; lag <= d; lag++){
    if(matches_one_lag_block(cur_ids, cur_n, ip, b, lag, d, base_offsets)){
      return 1;
    }
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Helper: persistence flag for one group vertex                              */
/* -------------------------------------------------------------------------- */
static double group_flag(Vertex g, const double *ip,
                          const int *sizes_int, Network *nwp){

  const int n1_total     = (int)ip[0];
  const int B            = (int)ip[1];
  const int d            = (int)ip[2];
  const int L            = (int)ip[3];

  const int base_aoff    = 4 + L;
  const int base_neff    = 4 + L + B;
  const int base_gtb     = 4 + L + 2 * B;
  const int base_offsets = 4 + L + 2 * B + n1_total;

  UNUSED_WARNING(d);

  /* Look up which block this group belongs to (O(1) direct lookup) */
  int g_off = (int)g - (n1_total + 1);    /* 0-based among group vertices */
  int b     = (int)ip[base_gtb + g_off];  /* block index 1..B             */

  /* Per-block parameters */
  int actor_off = (int)ip[base_aoff + b - 1];  /* 0-based global offset */
  int n_eff     = (int)ip[base_neff + b - 1];  /* effective actor count  */

#if DEBUG_INERTIA_GROUPS
  Rprintf("[group_flag] g=%d g_off=%d -> b=%d actor_off=%d n_eff=%d\n",
          (int)g, g_off, b, actor_off, n_eff);
#endif

  /* Collect current group members (within block, sorted global ids) */
  int *ids = (int*)R_Calloc(n_eff, int);
  int ng   = collect_group_members_in_block(g, actor_off, n_eff, nwp, ids);

  if(ng == 0){
#if DEBUG_INERTIA_GROUPS
    Rprintf("[group_flag] g=%d -> empty group, flag=0\n", (int)g);
#endif
    R_Free(ids); return 0.0;
  }

  /* Size filter on current group */
  if(!in_sizes(ng, L, sizes_int)){
#if DEBUG_INERTIA_GROUPS
    Rprintf("[group_flag] g=%d ng=%d -> size filter REJECT, flag=0\n", (int)g, ng);
#endif
    R_Free(ids); return 0.0;
  }

  /* Match required in ANY lag (union Sigma^(d)_t), within the SAME block */
  int ok = matches_any_past_groups_block(ids, ng, ip, b, d, base_offsets);
  R_Free(ids);

#if DEBUG_INERTIA_GROUPS
  Rprintf("[group_flag] g=%d ng=%d -> all-lag match=%d flag=%.0f\n",
          (int)g, ng, ok, ok ? (double)ng : 0.0);
#endif

  /* Count actors in replicated groups (not groups): each persistent group
   * contributes its size ng to the statistic T_inertia = #{i : group_t(i) replicated}. */
  return ok ? (double)ng : 0.0;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: inertia_groups (one-toggle)                              */
/* -------------------------------------------------------------------------- */
C_CHANGESTAT_FN(c_inertia_groups){
  ZERO_ALL_CHANGESTATS(0);

  const double *ip = INPUT_PARAM;
  const int n1_total = (int)ip[0];
  const int L        = (int)ip[3];

  int *sizes_int = NULL;
  if(L > 0){
    sizes_int = (int*)R_Calloc(L, int);
    for(int i = 0; i < L; i++) sizes_int[i] = (int)ip[4 + i];
  }

  /* Identify actor and group vertices from the toggled edge */
  Vertex a  = tail, b_v = head;
  Vertex actor = (a <= (Vertex)n1_total) ? a  : b_v;
  Vertex group = (a <= (Vertex)n1_total) ? b_v : a;

#if DEBUG_INERTIA_GROUPS
  {
    const int B = (int)ip[1];
    const int d = (int)ip[2];
    Rprintf("[inertia_groups] toggle actor=%d group=%d | n1_total=%d B=%d d=%d L=%d\n",
            (int)actor, (int)group, n1_total, B, d, L);
  }
#endif

  UNUSED_WARNING(actor);

  /* Compute flag before toggle */
  double F_before = group_flag(group, ip, sizes_int, nwp);

  /* Virtual toggle */
  TOGGLE(a, b_v);

  /* Compute flag after toggle */
  double F_after = group_flag(group, ip, sizes_int, nwp);

  /* Undo toggle */
  TOGGLE(a, b_v);

  CHANGE_STAT[0] += (F_after - F_before);

#if DEBUG_INERTIA_GROUPS
  Rprintf("[inertia_groups] group=%d before=%.0f after=%.0f delta=%.0f\n",
          (int)group, F_before, F_after, (F_after - F_before));
#endif

  if(sizes_int) R_Free(sizes_int);
}
