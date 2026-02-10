// changestat_inertia_groups.c
/**
 * @file changestat_inertia_groups.c
 * @brief Change statistic for the ERPM inertial term `inertia_groups`.
 *
 * @details
 *  This change statistic is designed to be agnostic to the pseudo-longitudinal
 *  paradigm used by erpm_long():
 *    - PLS (sequential): a single bipartite network per time.
 *    - PLE (stacked): a block-diagonal ("stacked") bipartite meta-network.
 *
 *  IMPORTANT DESIGN RULE:
 *  The changestat does NOT implement any PLS/PLE branching logic.
 *  The distinction is handled entirely by InitErgmTerm.inertia_groups, which
 *  packs INPUT_PARAM so that the same C code applies in both cases.
 *
 *  What is counted:
 *    Among current groups (mode-2 vertices), count those whose exact actor
 *    membership set matches a group observed in the past over a window of
 *    d = past_influence lags, optionally restricted by a size filter.
 *
 *  IMPORTANT (depth semantics):
 *    A group contributes 1 iff its signature appears identically in
 *    ALL required observed past partitions across lags 1..d (intersection over lags).
 *
 *  Locality:
 *    This change statistic is local to a single bipartite edge toggle. Only the
 *    affected group can change its persistence flag for that toggle, so:
 *
 *      Δ = flag_after(group) - flag_before(group)
 *
 *    Note on b1part moves:
 *      Under the b1part constraint, a "move" of an actor from old group to new
 *      group is implemented as TWO toggles (remove old edge, add new edge). ERGM
 *      accounts for both impacted groups through the two toggles of the move.
 *
 *    where flag(group)=1 iff:
 *      - current group size passes optional size filter S, and
 *      - current actor set equals one of the past groups for the SAME BLOCK (PLE)
 *        or the only block (PLS), for EACH lag 1..d.
 *
 *  INPUT_PARAM packing (by InitErgmTerm.inertia_groups):
 *
 *    INPUT_PARAM = c(
 *      n1_total, n_block, G_block, B, d, L, sizes[1:L],
 *      offsets[1:(B*d)],  # 0-based offsets into INPUT_PARAM (double-coded ints)
 *      # data blocks appended in any deterministic order (here: block-major, lag-major):
 *      # for b=1..B:
 *      #   for lag=1..d:
 *      #     M,
 *      #       len_1, ids...
 *      #       len_2, ids...
 *      #       ...
 *    )
 *
 *  Notes:
 *    - All values are stored as doubles and cast to int in C.
 *    - Actor IDs stored in past groups are GLOBAL actor vertex ids in the current
 *      network's actor space (1..n1_total).
 *    - In PLS, we set B=1, n_block=n1_total, G_block=n1_total.
 *    - In PLE (stacked), we set n1_total = n_block*B, and groups are laid out
 *      blockwise in mode-2; the block index of a group vertex is derived from
 *      its vertex id using (n1_total, G_block, B).
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"      /* R_Calloc/R_Free */
#include <R_ext/Print.h>
#include <string.h>

#define DEBUG_INERTIA_GROUPS 0
#define UNUSED_WARNING(x) (void)(x)

/* -------------------------------------------------------------------------- */
/* Helper: size filter                                                        */
/* -------------------------------------------------------------------------- */
static inline int in_sizes(int n, int L, const int *sizes){
  if(L==0) return 1;
  for(int i=0; i<L; i++){
    if(sizes[i] == n) return 1;
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Helper: determine the block index of a group vertex                         */
/* -------------------------------------------------------------------------- */
/*
 * We assume the standard bipartite convention:
 * - actor vertices: 1..n1_total
 * - group vertices: (n1_total+1)..N
 *
 * For PLS:
 * - B = 1 => block = 1 for all groups.
 *
 * For PLE (stacked):
 * - groups are arranged blockwise, each block having G_block group vertices.
 * - group offset among groups is (g - (n1_total+1)) in [0 .. G_block*B - 1]
 * - block = floor(offset / G_block) + 1
 *
 * Returns block in 1..B. If anything looks inconsistent, returns 1 as a safe fallback.
 */
static inline int group_block_index(Vertex g, int n1_total, int G_block, int B){
  if(B <= 1) return 1;
  if(g <= (Vertex)n1_total) return 1;
  int off = (int)g - (n1_total + 1); /* 0-based among group vertices */
  if(off < 0) return 1;
  int b = (off / G_block) + 1;
  if(b < 1) b = 1;
  if(b > B) b = B;
  return b;
}

/* -------------------------------------------------------------------------- */
/* Helper: build sorted member list of a group (actors only)                  */
/* -------------------------------------------------------------------------- */
/*
 * Collect membership of group vertex g:
 * - deduplicate via seen[] using neighbor traversal
 * - then emit sorted actor ids by scanning seen[] from 1..n1_total
 *
 * Returns:
 * - out_ids : int array (allocated by caller) filled with sorted actor ids
 * - return value: group size (k)
 */
static int collect_group_members_sorted(Vertex g, int n1_total, Network *nwp, int *out_ids){
  unsigned char *seen = (unsigned char*)R_Calloc(n1_total, unsigned char);
  Edge e; Vertex h;

  /* Neighbors through OUT edges */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1_total){
      seen[(int)h - 1] = 1;
    }
  }
  /* Neighbors through IN edges */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1_total){
      seen[(int)h - 1] = 1;
    }
  }

  int k = 0;
  for(int i=1; i<=n1_total; i++){
    if(seen[i-1]){
      out_ids[k++] = i;
    }
  }

  R_Free(seen);
  return k;
}

/* -------------------------------------------------------------------------- */
/* Helper: check membership in one (b,lag) observed partition block           */
/* -------------------------------------------------------------------------- */
/*
 * Returns 1 iff cur_ids (length cur_n, sorted) matches at least one observed
 * group signature stored in the (b,lag) block.
 */
static int matches_one_lag_block(const int *cur_ids, int cur_n,
                                const double *ip,
                                int b, int lag,
                                int d, int L){
  int base_offsets = 6 + L;

  int idx = (b - 1) * d + (lag - 1);           /* 0-based in offsets table */
  int pos = (int)ip[base_offsets + idx];       /* 0-based position in ip */
  if(pos < 0) return 0;

  int M = (int)ip[pos++];
  for(int j=0; j<M; j++){
    int len = (int)ip[pos++];
    if(len == cur_n){
      int ok = 1;
      for(int u=0; u<len; u++){
        int pid = (int)ip[pos + u];
        if(pid != cur_ids[u]) { ok = 0; break; }
      }
      if(ok) return 1;
    }
    pos += len;
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Helper: compare current group to past groups (block b, ALL lags)           */
/* -------------------------------------------------------------------------- */
/*
 * Depth semantics:
 *   Require match in EVERY lag 1..d (intersection over lags).
 */
static int matches_all_past_groups_block(const int *cur_ids, int cur_n,
                                        const double *ip,
                                        int b, int d, int L){
  for(int lag=1; lag<=d; lag++){
    if(!matches_one_lag_block(cur_ids, cur_n, ip, b, lag, d, L)){
      return 0;
    }
  }
  return 1;
}

/* -------------------------------------------------------------------------- */
/* Helper: persistence flag for one group vertex                              */
/* -------------------------------------------------------------------------- */
static double group_flag(Vertex g,
                         const double *ip,
                         const int *sizes_int,
                         Network *nwp){

  const int n1_total = (int)ip[0];
  const int n_block  = (int)ip[1];
  const int G_block  = (int)ip[2];
  const int B        = (int)ip[3];
  const int d        = (int)ip[4];
  const int L        = (int)ip[5];

  UNUSED_WARNING(n_block);

  /* Determine which block this group belongs to (PLE) or 1 (PLS) */
  const int b = group_block_index(g, n1_total, G_block, B);

  /* Build current group member list */
  int *ids = (int*)R_Calloc(n1_total, int);
  int ng = collect_group_members_sorted(g, n1_total, nwp, ids);

  /* Empty group: never persistent */
  if(ng == 0){
    R_Free(ids);
    return 0.0;
  }

  /* Size filter on CURRENT group */
  if(!in_sizes(ng, L, sizes_int)){
    R_Free(ids);
    return 0.0;
  }

  /* Exact match required in ALL lags, within the SAME block */
  int ok = matches_all_past_groups_block(ids, ng, ip, b, d, L);
  R_Free(ids);

  return ok ? 1.0 : 0.0;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: inertia_groups (one-toggle)                              */
/* -------------------------------------------------------------------------- */
C_CHANGESTAT_FN(c_inertia_groups){
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(mtp);
  UNUSED_WARNING(edgestate);

  const double *ip = INPUT_PARAM;

  const int n1_total = (int)ip[0];
  const int d        = (int)ip[4];
  const int L        = (int)ip[5];

  UNUSED_WARNING(d);

  /* Copy size filter as ints for cheap comparisons */
  int *sizes_int = NULL;
  if(L > 0){
    sizes_int = (int*)R_Calloc(L, int);
    for(int i=0; i<L; i++) sizes_int[i] = (int)ip[6 + i];
  }

  #if DEBUG_INERTIA_GROUPS
    {
      const int n_block = (int)ip[1];
      const int G_block = (int)ip[2];
      const int B       = (int)ip[3];
      Rprintf("[inertia_groups] n1_total=%d n_block=%d G_block=%d B=%d d=%d L=%d\n",
              n1_total, n_block, G_block, B, d, L);
    }
  #endif

  /* Identify actor and group vertices */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1_total) ? a : b;
  Vertex group = (a <= (Vertex)n1_total) ? b : a;
  UNUSED_WARNING(actor);

  /* Compute before */
  double F_before = group_flag(group, ip, sizes_int, nwp);

  /* Virtual toggle */
  TOGGLE(a, b);

  /* Compute after */
  double F_after  = group_flag(group, ip, sizes_int, nwp);

  /* Undo toggle */
  TOGGLE(a, b);

  CHANGE_STAT[0] += (F_after - F_before);

  #if DEBUG_INERTIA_GROUPS
    Rprintf("[inertia_groups] group=%d before=%.0f after=%.0f delta=%.0f\n",
            (int)group, F_before, F_after, (F_after - F_before));
  #endif

  if(sizes_int) R_Free(sizes_int);
}
