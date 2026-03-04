/**
 * @file changestat_inertia_groups.c
 * @brief Change statistic for the ERPM inertial term 'inertia_groups'.
 *
 * @details
 * This file implements the C-side change statistic used by the ERPM term
 * 'inertia_groups'. The code is deliberately designed to remain independent 
 * of the longitudinal mode selected in erpm_long():
 *
 * - PLS (sequential): one bipartite network per time step.
 * - PLE (stacked): one block-diagonal bipartite meta-network.
 *
 * The change statistic does not branch on PLS vs PLE. The distinction is
 * handled entirely by InitErgmTerm.inertia_groups, which packs INPUT_PARAM
 * so that the same C code applies in both settings.
 *
 * What is counted:
 * For each current group (mode-2 vertex), we check whether its exact actor
 * membership matches a group observed in the past over a window of
 * d = past_influence lags. An optional size filter may restrict which
 * current groups are eligible.
 *
 * Depth semantics:
 * A group contributes 1 iff its signature appears identically in ALL
 * required past partitions across lags 1..d (intersection over lags).
 *
 * Locality:
 * The statistic is local to a single bipartite edge toggle. Only the
 * impacted group can change its persistence status for that toggle:
 *
 *   Δ = flag_after(group) - flag_before(group)
 *
 * Under the 'b1part' constraint, a reassignment of one actor corresponds
 * to two toggles (remove old edge, add new edge). ERGM evaluates both
 * toggles separately, so both affected groups are accounted for.
 *
 * INPUT_PARAM layout (packed in R):
 *
 *   c(
 *     n1_total, n_block, G_block, B, d, L,
 *     sizes[1:L],
 *     offsets[1:(B*d)],
 *     # followed by observed group signatures per (block, lag)
 *   )
 *
 * All values are stored as doubles and cast to int in C.
 * Actor IDs stored in past signatures are global actor vertex ids
 * (1..n1_total).
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
  if(L==0) return 1;
  for(int i=0; i<L; i++){
    if(sizes[i] == n) return 1;
  }
  return 0;
}

/**
 * @brief Determine the block index of a group vertex.
 *
 * In PLS (B=1), all groups belong to block 1.
 * In PLE, group vertices are arranged blockwise, each block
 * containing G_block groups.
 *
 * @param g Group vertex id.
 * @param n1_total Number of actor vertices.
 * @param G_block Number of groups per block.
 * @param B Number of blocks.
 *
 * @return Block index in 1..B.
 */
static inline int group_block_index(Vertex g, int n1_total, int G_block, int B){
  if(B <= 1) return 1;
  if(g <= (Vertex)n1_total) return 1;

  int off = (int)g - (n1_total + 1);
  if(off < 0) return 1;

  int b = (off / G_block) + 1;
  if(b < 1) b = 1;
  if(b > B) b = B;

  return b;
}

/**
 * @brief Collect sorted actor members of a group.
 *
 * The function traverses both IN and OUT edges of the group vertex,
 * records actor neighbors, and emits a sorted list of actor ids.
 *
 * @param g Group vertex.
 * @param n1_total Number of actor vertices.
 * @param nwp Network pointer.
 * @param out_ids Output array (allocated by caller).
 *
 * @return Group size.
 */
static int collect_group_members_sorted(Vertex g, int n1_total, Network *nwp, int *out_ids){
  unsigned char *seen = (unsigned char*)R_Calloc(n1_total, unsigned char);
  Edge e; Vertex h;

  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1_total){
      seen[(int)h - 1] = 1;
    }
  }
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

/**
 * @brief Check whether current membership matches one past lag.
 *
 * @param cur_ids Sorted current actor ids.
 * @param cur_n Current group size.
 * @param ip INPUT_PARAM.
 * @param b Block index.
 * @param lag Lag index.
 * @param d Number of lags.
 * @param L Size filter length.
 *
 * @return 1 if a matching past group is found, 0 otherwise.
 */
static int matches_one_lag_block(const int *cur_ids, int cur_n,
                                const double *ip,
                                int b, int lag,
                                int d, int L){
  int base_offsets = 6 + L;

  int idx = (b - 1) * d + (lag - 1);
  int pos = (int)ip[base_offsets + idx];
  if(pos < 0) return 0;

  int M = (int)ip[pos++];
  for(int j=0; j<M; j++){
    int len = (int)ip[pos++];
    if(len == cur_n){
      int ok = 1;
      for(int u=0; u<len; u++){
        if((int)ip[pos + u] != cur_ids[u]){
          ok = 0;
          break;
        }
      }
      if(ok) return 1;
    }
    pos += len;
  }
  return 0;
}

/**
 * @brief Require exact match across all lags (intersection semantics).
 *
 * @return 1 if the group matches in every lag, 0 otherwise.
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

/**
 * @brief Compute persistence flag for a single group vertex.
 *
 * A group is persistent (flag=1) iff:
 *  - it is non-empty,
 *  - it passes the size filter (if any),
 *  - its exact membership matches a past group in all required lags.
 *
 * @return 1.0 if persistent, 0.0 otherwise.
 */
static double group_flag(Vertex g,
                         const double *ip,
                         const int *sizes_int,
                         Network *nwp){

  const int n1_total = (int)ip[0];
  const int G_block  = (int)ip[2];
  const int B        = (int)ip[3];
  const int d        = (int)ip[4];
  const int L        = (int)ip[5];

  const int b = group_block_index(g, n1_total, G_block, B);

  int *ids = (int*)R_Calloc(n1_total, int);
  int ng = collect_group_members_sorted(g, n1_total, nwp, ids);

  if(ng == 0){
    R_Free(ids);
    return 0.0;
  }

  if(!in_sizes(ng, L, sizes_int)){
    R_Free(ids);
    return 0.0;
  }

  int ok = matches_all_past_groups_block(ids, ng, ip, b, d, L);
  R_Free(ids);

  return ok ? 1.0 : 0.0;
}

/**
 * @brief Change statistic for 'inertia_groups' (one-toggle version).
 *
 * The statistic evaluates the persistence flag of the impacted group
 * before and after a virtual toggle and returns the difference.
 */
C_CHANGESTAT_FN(c_inertia_groups){
  ZERO_ALL_CHANGESTATS(0);

  const double *ip = INPUT_PARAM;
  const int n1_total = (int)ip[0];
  const int L = (int)ip[5];

  int *sizes_int = NULL;
  if(L > 0){
    sizes_int = (int*)R_Calloc(L, int);
    for(int i=0; i<L; i++) sizes_int[i] = (int)ip[6 + i];
  }

  Vertex a = tail, b = head;
  Vertex group = (a <= (Vertex)n1_total) ? b : a;

  double F_before = group_flag(group, ip, sizes_int, nwp);

  TOGGLE(a, b);
  double F_after  = group_flag(group, ip, sizes_int, nwp);
  TOGGLE(a, b);

  CHANGE_STAT[0] += (F_after - F_before);

  if(sizes_int) R_Free(sizes_int);
}
