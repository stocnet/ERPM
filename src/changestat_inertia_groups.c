/**
 * @file changestat_inertia_groups.c
 * @brief Change statistic for the ERPM inertial term `inertia_groups`.
 *
 * @details
 *  This is a longitudinal (inertial) ERGM term intended to be used with
 *  erpm_long(). Past information is attached to the current network as
 *  network attributes before calling ergm().
 *
 *  The statistic counts, among current groups (group-mode vertices),
 *  those whose exact actor membership set matches at least one group observed
 *  in the past over a window of d = past_influence lags.
 *
 *  This change statistic is local to a single toggle (actor-group membership edge).
 *  Only the affected group can change its persistence status, so:
 *
 *    Δ = flag_after(group) - flag_before(group),
 *
 *  where flag(group) = 1 if:
 *    - group size passes optional size filter S, and
 *    - group actor set equals one of the past groups across any lag.
 *
 *  INPUT_PARAM is pre-packed by InitErgmTerm.inertia_groups as:
 *
 *    INPUT_PARAM = c(
 *      n1, d, L, sizes[1:L],
 *      # for each lag=1..d:
 *      M_lag,
 *        len_1, ids...
 *        len_2, ids...
 *        ...
 *    )
 *
 *  All values are stored as doubles and cast to int in C.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"      /* Calloc/Free */
#include <R_ext/Print.h>
#include <string.h>

#define DEBUG_INERTIA_GROUPS 0
#define UNUSED_WARNING(x) (void)x

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
/* Helper: build sorted member list of a group (actors only)                  */
/* -------------------------------------------------------------------------- */
/*
 * Collect membership of group vertex g:
 * - deduplicate via seen[] using neighbor traversal
 * - then emit sorted actor ids by scanning seen[] from 1..n1
 *
 * Returns:
 * - *out_ids : int array of length *out_n (allocated by caller)
 * - out_n    : group size
 */
static int collect_group_members_sorted(Vertex g, int n1, Network *nwp, int *out_ids){
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char);
  Edge e; Vertex h;

  /* Neighbors through OUT edges */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      seen[(int)h - 1] = 1;
    }
  }
  /* Neighbors through IN edges */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      seen[(int)h - 1] = 1;
    }
  }

  int k = 0;
  for(int i=1; i<=n1; i++){
    if(seen[i-1]){
      out_ids[k++] = i;
    }
  }

  Free(seen);
  return k;
}

/* -------------------------------------------------------------------------- */
/* Helper: compare current group to past groups in INPUT_PARAM                */
/* -------------------------------------------------------------------------- */
static int matches_any_past_group(const int *cur_ids, int cur_n,
                                 const double *ip, int n1, int d, int L){
  /* ip points to beginning of INPUT_PARAM. We must skip header and sizes. */
  int pos = 0;

  /* header: n1, d, L */
  pos += 3;

  /* sizes: L entries */
  pos += L;

  /* Now per lag blocks */
  for(int lag=1; lag<=d; lag++){
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
      /* skip ids */
      pos += len;
    }
  }

  UNUSED_WARNING(n1);
  return 0;
}

/* -------------------------------------------------------------------------- */
/* Helper: persistence flag for one group vertex                              */
/* -------------------------------------------------------------------------- */
static double group_flag(Vertex g,
                         int n1,
                         int d,
                         int L,
                         const int *sizes_int,
                         const double *ip,
                         Network *nwp){
  /* Build current group member list */
  int *ids = (int*)Calloc(n1, int);
  int ng = collect_group_members_sorted(g, n1, nwp, ids);

  /* Empty group: never persistent */
  if(ng == 0){
    Free(ids);
    return 0.0;
  }

  /* Size filter on CURRENT group */
  if(!in_sizes(ng, L, sizes_int)){
    Free(ids);
    return 0.0;
  }

  /* Exact match against any past group over any lag */
  int ok = matches_any_past_group(ids, ng, ip, n1, d, L);
  Free(ids);

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

  const int n1 = (int)ip[0];
  const int d  = (int)ip[1];
  const int L  = (int)ip[2];

  /* Copy size filter as ints for cheap comparisons */
  int *sizes_int = NULL;
  if(L > 0){
    sizes_int = (int*)Calloc(L, int);
    for(int i=0; i<L; i++) sizes_int[i] = (int)ip[3 + i];
  }

  #if DEBUG_INERTIA_GROUPS
    Rprintf("[inertia_groups] n1=%d d=%d L=%d\n", n1, d, L);
  #endif

  /* Identify actor and group vertices */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;
  Vertex group = (a <= (Vertex)n1) ? b : a;
  UNUSED_WARNING(actor);

  /* Compute before */
  double F_before = group_flag(group, n1, d, L, sizes_int, ip, nwp);

  /* Virtual toggle */
  TOGGLE(a, b);

  /* Compute after */
  double F_after  = group_flag(group, n1, d, L, sizes_int, ip, nwp);

  /* Undo toggle */
  TOGGLE(a, b);

  CHANGE_STAT[0] += (F_after - F_before);

  #if DEBUG_INERTIA_GROUPS
    Rprintf("[inertia_groups] group=%d before=%.0f after=%.0f delta=%.0f\n",
            (int)group, F_before, F_after, (F_after - F_before));
  #endif

  if(sizes_int) Free(sizes_int);
}
