/**
 * @file changestat_cliques.c
 * @brief  Change statistic for the ERPM term `cliques(k, normalized)` (multi-toggle).
 * 
 * @author Jérémie Chichignoud
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
 *  `cliques(k, normalized)`. The effect counts, for each group in the
 *  group mode, the number of actor k-cliques that are induced by group
 *  memberships in the actor mode.
 *
 *  IMPORTANT (multi-toggle):
 *  - Swap/split/merge proposals may be represented as a list of toggles.
 *  - Therefore we implement this effect with D_CHANGESTAT_FN and process toggles
 *    sequentially, temporarily applying each toggle so subsequent toggles see
 *    consistent intermediate degrees, then undo all toggles before returning.
 *
 *  Actor / group setup:
 *  - The network is bipartite:
 *    - actor mode = vertices 1..BIPARTITE,
 *    - group mode = vertices (BIPARTITE+1)..N_NODES.
 *  - Membership edges connect one actor to one group.
 *
 *  Let n_g be the size of a group g (degree of the group vertex in the bipartite
 *  representation). For each requested k:
 *
 *  Raw statistic:
 *    - For k >= 2:   cliques_k(y) = Σ_g C(n_g, k)
 *    - For k == 1:   cliques_1(y) = #{ g : n_g == 1 }
 *
 *  Group-size-normalized statistic:
 *    - For k >= 2:   cliques_k^{grp}(y) = Σ_g C(n_g, k) / n_g
 *    - For k == 1:   identical to raw (since n_g = 1 for contributing groups)
 *
 *  Change statistics (single toggle intuition):
 *    - Raw, k >= 2:
 *        addition: Δ_raw =  C(n_g,   k-1)
 *        deletion: Δ_raw = -C(n_g-1, k-1)
 *    - Raw, k == 1:
 *        addition: 0->1:+1 ; 1->2:-1 ; else:0
 *        deletion: 2->1:+1 ; 1->0:-1 ; else:0
 *    - Normalized, k >= 2:
 *        old_contrib = (n_g  >= k && n_g  > 0) ? C(n_g,  k) / n_g  : 0
 *        new_contrib = (n_g' >= k && n_g' > 0) ? C(n_g', k) / n_g' : 0
 *        Δ_grp = new_contrib - old_contrib
 *
 *  Scaling / mode flag via INPUT_PARAM:
 *    INPUT_PARAM = (k_1, scale_1, k_2, scale_2, ..., k_J, scale_J)
 *    - scale_j > 0 => raw mode
 *    - scale_j < 0 => group-size-normalized mode
 *    - |scale_j|   => multiplicative factor
 *
 *  For each toggle i, and each stat j:
 *    CHANGE_STAT[j] += |scale_j| * Δ_j
 *
 *  R interface:
 *    - InitErgmTerm.cliques MUST return d_func=TRUE to match this D_ signature.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>

/**
 * @def DEBUG_CLIQUES
 * @brief Enable or disable verbose debugging output for ::d_cliques.
 *
 * When set to 1, the change-statistic function prints diagnostic information:
 *  - multi-toggle detection (ntoggles),
 *  - per-toggle endpoints and group vertex,
 *  - degrees before/after, mode (raw/grp), Δ and cumulative stats.
 *
 * When set to 0, no debug output is produced.
 */
#define DEBUG_CLIQUES 0
#define UNUSED_VARIABLE(x) (void)x

/* -------------------------------------------------------------------------- */
/* Change statistic: cliques(k, normalized) (multi-toggle)                     */
/* -------------------------------------------------------------------------- */

D_CHANGESTAT_FN(d_cliques){

#if DEBUG_CLIQUES
  static int seen = 0;
  if(ntoggles > 1 && seen < 50){
    Rprintf("[cliques] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen++;
  }
#endif

  /* 1) Reset the output buffer for this proposal (all toggles). */
  ZERO_ALL_CHANGESTATS();

  /* 2) Actor-mode size. */
  const int n1 = BIPARTITE;

  /* 3) Process toggles sequentially (intermediate-state aware). */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex t = TAIL(i);
    Vertex h = HEAD(i);

    /* Determine current edge state BEFORE toggling, on the current intermediate
     * network state.
     *
     * - For directed networks, use IS_OUTEDGE(t,h).
     * - For undirected (typical memberships), use IS_UNDIRECTED_EDGE(t,h), which
     *   is robust to endpoint order.
     */
    int edgestate = DIRECTED ? IS_OUTEDGE(t, h) : IS_UNDIRECTED_EDGE(t, h);

    /* Identify the group-mode vertex (index > n1). */
    Vertex v2 = (t > (Vertex)n1) ? t : h;

#if DEBUG_CLIQUES
    if(v2 <= (Vertex)n1){
      Rprintf("[d_cliques][WARN] toggle #%d has no group endpoint: tail=%d head=%d (n1=%d)\n",
              i, (int)t, (int)h, n1);
    }
#endif

    /* Group size before this toggle, in the intermediate state. */
    int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);

    /* edgestate==0 => addition, edgestate==1 => deletion */
    const int is_add = edgestate ? 0 : 1;

    /* Loop over all statistics (vectorized k). */
    for(int j = 0; j < N_CHANGE_STATS; ++j){

      int    k         = (int)   INPUT_PARAM[2*j + 0];
      double scale_raw = (double)INPUT_PARAM[2*j + 1];

      int    use_grp   = (scale_raw < 0.0);
      double scale_abs = use_grp ? -scale_raw : scale_raw; /* >= 0 */

      double delta = 0.0;

      if(k == 1){
        /* k==1 counts singleton groups (size exactly 1).
         * Normalized variant is identical.
         */
        if(is_add){
          if(deg_old == 0)      delta = +1.0;
          else if(deg_old == 1) delta = -1.0;
        }else{
          if(deg_old == 2)      delta = +1.0;
          else if(deg_old == 1) delta = -1.0;
        }
      }else{
        /* k>=2 : raw vs group-size-normalized. */
        if(!use_grp){
          /* Raw (un-toggle formulas). */
          if(is_add){
            if(deg_old >= k-1) delta = CHOOSE(deg_old, k-1);
          }else{
            if(deg_old >= 1 && (deg_old - 1) >= (k - 1)){
              delta = -CHOOSE(deg_old - 1, k - 1);
            }
          }
        }else{
          /* Group-size-normalized: delta = new_contrib - old_contrib. */
          int n_old = deg_old;
          int n_new = is_add ? (n_old + 1) : (n_old - 1);

          double old_contrib = 0.0;
          double new_contrib = 0.0;

          if(n_old >= k && n_old > 0){
            old_contrib = CHOOSE(n_old, k) / (double)n_old;
          }
          if(n_new >= k && n_new > 0){
            new_contrib = CHOOSE(n_new, k) / (double)n_new;
          }

          delta = new_contrib - old_contrib;
        }
      }

      CHANGE_STAT[j] += scale_abs * delta;

#if DEBUG_CLIQUES
      Rprintf("[D:d_cliques] i=%d tail=%d head=%d | group=%d | edgestate=%d | add=%d | "
              "deg_old=%d | k=%d | mode=%s | delta=%.6f | scale_abs=%g | stat[%d]=%.6f\n",
              i, (int)t, (int)h, (int)v2, (int)edgestate, (int)is_add,
              deg_old, k, use_grp ? "grp" : "raw",
              delta, scale_abs, j, CHANGE_STAT[j]);
#endif
    }

    /* Temporarily apply this toggle so subsequent toggles see updated degrees. */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* 4) Undo temporary toggles to restore the original network state. */
  UNDO_PREVIOUS_TOGGLES(i);
}
