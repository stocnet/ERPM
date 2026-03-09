/**
 * @file changestat_cov_diff_GW.c
 * @brief  Change statistic for the ERPM term `cov_diff_GW` (multi-toggle form).
 * 
 * @author Jérémie Chichignoud
 *
 * @details
 *  This file implements the \pkg{ergm} change statistic for the ERPM effect
 *  `cov_diff_GW`, which applies a geometrically weighted transformation
 *  to the family of cov_diff statistics over all subset sizes k ≥ 2
 *  inside each group in the group mode.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  The bipartite network represents a partition:
 *    - actor mode  = vertices representing actors (individuals),
 *    - group mode  = vertices representing structural groups.
 *
 *  Membership is encoded by edges between actors and groups. Let:
 *    - A(g)        = set of actors adjacent to group vertex g,
 *    - n_g         = |A(g)|, size of group g,
 *    - x_i         = numeric covariate value for actor i in the actor mode.
 *
 *  For each integer k ≥ 2 and each group g, define:
 *
 *      c_k(g)
 *        = ∑_{S ⊂ A(g), |S| = k} D(S),
 *
 *  where:
 *
 *      D(S) = max_{i∈S} x_i - min_{i∈S} x_i.
 *
 *  The global cov_diff statistic of order k is:
 *
 *      c_k(p; x) = ∑_g c_k(g),
 *
 *  where p encodes the partition via the bipartite membership structure.
 *
 *  The geometrically weighted cov_diff statistic is defined as:
 *
 *      T_GW(p; x, λ)
 *        = ∑_{k≥2} (-1/λ)^{k-1} c_k(p; x),
 *
 *  for a decay parameter λ > 1. In practice the sum over k is finite,
 *  since groups have finite size. For multiple decay parameters
 *  λ_1, …, λ_L, the term is vectorised and returns L statistics, one
 *  per λ_ℓ.
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
 *      and one group (vertex > n1),
 *    - only the group in the group mode incident to the toggled edge
 *      can change its cov_diff_GW contribution.
 *
 *  ------------------------------------------------------------
 *  INPUT_PARAM layout (from InitErgmTerm.cov_diff_GW)
 *  ------------------------------------------------------------
 *
 *  The R initialiser packs the parameters into INPUT_PARAM as:
 *
 *    INPUT_PARAM = c(
 *      n1,
 *      L,
 *      lambda[1:L],
 *      x[1:n1]
 *    )
 *
 *  where:
 *    - n1         = number of actors in the actor mode,
 *    - L          = number of decay parameters λ_ℓ,
 *    - lambda[ ]  = vector of λ_ℓ values (λ_ℓ > 1),
 *    - x[ ]       = numeric covariate on actors, length n1.
 *
 *  In C, this becomes:
 *
 *    INPUT_PARAM[0]          = n1
 *    INPUT_PARAM[1]          = L
 *    INPUT_PARAM[2..(1+L)]   = lambda[0..L-1]     (0-based indexing in C)
 *    INPUT_PARAM[2+L..]      = x[0..n1-1]
 *
 *  The term returns L statistics:
 *
 *    - N_CHANGE_STATS = L,
 *    - CHANGE_STAT[l] corresponds to λ_l = lambda[l].
 *
 *  ------------------------------------------------------------
 *  Local change under a toggle
 *  ------------------------------------------------------------
 *
 *  A membership toggle flips membership of one actor in one group:
 *    - addition  : actor becomes a member of the group,
 *    - deletion  : actor leaves the group.
 *
 *  Only the affected group g can change its contribution. For this group,
 *  define:
 *
 *    c_k^-(g) = c_k(g) before the toggle,
 *    c_k^+(g) = c_k(g) after  the toggle.
 *
 *  For a given λ, the local change is:
 *
 *    Δ T_GW(λ)
 *      = ∑_{k=2}^{K_max} (-1/λ)^{k-1} ( c_k^+(g) - c_k^-(g) ),
 *
 *  where K_max is max(n_g^-, n_g^+), the maximum group size before or
 *  after the toggle. All c_k(g) are recomputed locally for the group
 *  before and after a virtual toggle:
 *
 *    1. Extract all actors in the group (actor mode) using adjacency.
 *    2. Enumerate all k-subsets for k = 2..n_g and sum D(S).
 *    3. Store the values c_k(g) into arrays ck_before[k] and ck_after[k].
 *    4. Combine them with the weights (-1/λ)^{k-1}.
 *
 *  ------------------------------------------------------------
 *  Implementation in \pkg{ergm} (multi-toggle / D_CHANGESTAT_FN)
 *  ------------------------------------------------------------
 *
 *  IMPORTANT:
 *  - This effect MUST support multi-toggle proposals (swap/split/merge), which
 *    ergm decomposes into a sequence of edge toggles (tails/heads arrays).
 *  - Therefore, the changestat is implemented with D_CHANGESTAT_FN (D_ entrypoint),
 *    and processes toggles sequentially.
 *
 *  Sequential evaluation rule (same pattern as squared_sizes):
 *    - For each toggle i:
 *        1) Compute the delta for THAT toggle under the current intermediate state.
 *        2) Temporarily apply the toggle (TOGGLE_IF_MORE_TO_COME) so subsequent
 *           toggles see updated degrees/adjs if they touch the same group.
 *    - At the end, UNDO_PREVIOUS_TOGGLES restores the original network state.
 *
 *  Key subtlety:
 *    - To compute "before" and "after" contributions for a toggle, we apply a
 *      virtual TOGGLE(t,h) twice (toggle-on then toggle-off). This is local and
 *      leaves the network unchanged.
 *    - Then, separately, TOGGLE_IF_MORE_TO_COME(i) applies the actual toggle
 *      persistently in the intermediate state.
 *
 *  ------------------------------------------------------------
 *  Complexity
 *  ------------------------------------------------------------
 *
 *  For a single group g of size n_g:
 *
 *    - neighbour reconstruction (actors in g) is O(n_g),
 *    - for each k, there are C(n_g, k) subsets,
 *    - D(S) for a given subset S is computed in O(k).
 *
 *  The total complexity per toggle is:
 *
 *    O( ∑_{k=2}^{n_g} C(n_g, k) * k ),
 *
 *  multiplied by the number of λ values L. This is intended for moderate
 *  group sizes or small k in practice. There is no caching of per-group
 *  state across toggles.
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  The R initialiser (InitErgmTerm.cov_diff_GW):
 *    - validates λ values (λ > 1) and the numeric covariate on actors,
 *    - encodes the actor covariate as a numeric vector,
 *    - packs n1, L, the λ vector and the covariate into INPUT_PARAM,
 *    - sets emptynwstats = numeric(L), L coefficient names,
 *    - MUST return d_func = TRUE so ergm calls the D_ changestat entrypoint.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <R_ext/Print.h>
#include <math.h>
#include <string.h> /* memset */

/**
 * @def DEBUG_COV_DIFF_GW
 * @brief Enable verbose debugging output for ::d_cov_diff_GW.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - actor and group sizes for the affected group,
 *  - local c_k(g) values before and after a toggle,
 *  - local Δ T_GW for each lambda,
 *  - multi-toggle markers (ntoggles).
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_COV_DIFF_GW 0
#define UNUSED_VARIABLE(x) (void)x

/**
 * @def UNUSED_WARNING
 * @brief Utility macro to explicitly mark unused parameters.
 *
 * @param x Parameter or variable intentionally unused in the function.
 */
#define UNUSED_WARNING(x) (void)(x)

/* -------------------------------------------------------------------------- */
/* Helper: recursive enumeration of k-subsets                                 */
/* -------------------------------------------------------------------------- */

/**
 * @brief Recursively sum D(S) = max(x_i) - min(x_i) over all k-subsets S.
 *
 * @details
 *  This function performs a depth-first enumeration of all combinations
 *  of size k among ng actor indices stored in @p idxs[0..ng-1]. The
 *  current partial combination is stored in @p comb, as positions in
 *  idxs[]. When a complete combination of size k is formed, it computes:
 *
 *      D(S) = max_{i∈S} x_i - min_{i∈S} x_i
 *
 *  and adds it to the accumulator @p acc.
 *
 *  The recursion uses:
 *    - @p pos   = current depth (number of elements chosen so far),
 *    - @p start = first index in idxs[] that may be chosen at this depth.
 *
 * @param pos    Current depth in the combination (0..k).
 * @param start  Start index in @p idxs for the next choice.
 * @param k      Target subset size.
 * @param ng     Number of available actors in the group.
 * @param idxs   Array of actor indices (0-based indices into x[]).
 * @param x      Numeric covariate vector for actors.
 * @param comb   Working array of length k storing chosen positions in idxs[].
 * @param acc    Pointer to the accumulator for the sum of D(S) over all S.
 */
static void sum_D_rec(int pos, int start,
                      int k, int ng,
                      const int *idxs,
                      const double *x,
                      int *comb,
                      double *acc){

  if(pos == k){
    /* We have a complete combination in comb[0..k-1]. */
    double xmin = 0.0, xmax = 0.0;
    int first = 1;
    for(int t = 0; t < k; t++){
      int idx = idxs[ comb[t] ];   /* 0-based index into x[] */
      double val = x[idx];
      if(first){
        xmin = xmax = val;
        first = 0;
      }else{
        if(val < xmin) xmin = val;
        if(val > xmax) xmax = val;
      }
    }
    *acc += (xmax - xmin);
    return;
  }

  /* Choose the next element among idxs[start..ng-1]. */
  for(int i = start; i <= ng - (k - pos); i++){
    comb[pos] = i;
    sum_D_rec(pos + 1, i + 1, k, ng, idxs, x, comb, acc);
  }
}

/* -------------------------------------------------------------------------- */
/* Helper: cov_diff contributions c_k(g) for a single group                   */
/* -------------------------------------------------------------------------- */

/**
 * @brief Compute cov_diff contributions c_k(g) for all k ≥ 2 in a group.
 *
 * @details
 *  For a given group vertex @p g in the group mode, this function:
 *
 *    1. Reconstructs the set of member actors A(g) by:
 *         - traversing outgoing edges (group → actor),
 *         - traversing incoming edges (actor → group),
 *         - deduplicating actors using a "seen" array.
 *
 *    2. Stores the 0-based indices of these actors (for x[]) into @p idxs.
 *
 *    3. If n_g = |A(g)| < 2, returns immediately with:
 *         - *ng_out = n_g,
 *         - all c_k(g) implicitly 0 via @p ck (caller responsibility).
 *
 *    4. Otherwise, for each k in {2, …, n_g}:
 *         - enumerates all k-subsets S of A(g) via sum_D_rec(),
 *         - computes c_k(g) = ∑_S D(S),
 *         - stores the result in ck[k].
 *
 *  The arrays @p ck and @p idxs are allocated by the caller:
 *    - idxs must have length ≥ n1 (we use only the first n_g entries),
 *    - ck   must have length ≥ n1_for_ck, with n1_for_ck ≥ n1+1 so that
 *      we can safely write ck[k] for 2 ≤ k ≤ n_g.
 *
 *  Indices are:
 *    - actors in the actor mode are 1..n1,
 *    - x[idx] is the covariate of actor (idx + 1) in the actor mode.
 *
 * @param g          Group vertex in the group mode.
 * @param n1         Number of actors (size of the actor mode).
 * @param x          Pointer to numeric covariate values for actors.
 * @param nwp        Pointer to the network-plus workspace.
 * @param ng_out     Output: number of actors n_g in group g.
 * @param ck         Output: array of length ≥ n1_for_ck storing c_k(g) in ck[k].
 * @param n1_for_ck  Declared length of the ck buffer (must be ≥ n1+1).
 * @param idxs       Working array for actor indices in group g.
 */
static void group_covdiff_allk(Vertex g,
                               int n1,
                               const double *x,
                               Network *nwp,
                               int *ng_out,
                               double *ck,
                               int n1_for_ck,
                               int *idxs){

  Edge e;
  Vertex h;
  int ng = 0;

  /* Local marking of actors in the group (0..n1-1) to avoid double counting. */
  unsigned char *seen = (unsigned char *)R_Calloc(n1, unsigned char);

  /* OUT-neighbours: group → actor. */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        idxs[ng++] = idx;
      }
    }
  }

  /* IN-neighbours: actor → group. */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        idxs[ng++] = idx;
      }
    }
  }

  /* Clear marks for the actors seen in this group and release the buffer. */
  if(ng > 0){
    for(int i = 0; i < ng; i++){
      seen[ idxs[i] ] = 0;
    }
  }
  R_Free(seen);

  *ng_out = ng;

  /* If the group has fewer than 2 actors, all c_k(g) are zero. */
  if(ng < 2){
    return;
  }

  /* Guard: ck buffer must be large enough to address ck[k] for k up to n1. */
  if(n1_for_ck < n1 + 1){
    return;
  }

  /* Initialise ck[k] for k = 2..ng to 0 (ck[0], ck[1] remain 0). */
  for(int k = 2; k <= ng; k++){
    ck[k] = 0.0;
  }

  /* For each k, sum D(S) over all k-subsets S of the group. */
  for(int k = 2; k <= ng; k++){
    double sumD = 0.0;
    int *comb = (int *)R_Calloc(k, int);
    sum_D_rec(0, 0, k, ng, idxs, x, comb, &sumD);
    R_Free(comb);
    ck[k] = sumD;
  }

#if DEBUG_COV_DIFF_GW
  Rprintf("[cov_diff_GW][group_covdiff_allk] g=%d ng=%d : ", (int)g, ng);
  for(int k = 2; k <= ng; k++){
    Rprintf(" c_%d=%g", k, ck[k]);
  }
  Rprintf("\n");
#endif
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_diff_GW (multi-toggle / D_CHANGESTAT_FN)             */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_diff_GW` (multi-toggle).
 *
 * @details
 *  This is the \pkg{ergm} multi-toggle change-statistic function registered as
 *  ::d_cov_diff_GW via ::D_CHANGESTAT_FN. It computes the total change
 *  Δ T_GW(λ_ℓ) for each λ_ℓ over a proposal that may contain multiple toggles.
 *
 *  The layout of INPUT_PARAM is:
 *
 *    INPUT_PARAM[0]          = n1         (number of actors)
 *    INPUT_PARAM[1]          = L          (number of lambda values)
 *    INPUT_PARAM[2..1+L]     = lambda[0..L-1]
 *    INPUT_PARAM[2+L..]      = x[0..n1-1] (numeric covariate on actors)
 *
 *  For each toggle i in the proposal:
 *    1. Identify the actor vertex and the group vertex using n1.
 *    2. Compute c_k^-(g) for all k ≥ 2 in the CURRENT intermediate state.
 *    3. Apply a local virtual toggle (TOGGLE) on the actor–group edge.
 *    4. Compute c_k^+(g) for all k ≥ 2.
 *    5. Undo the virtual toggle (TOGGLE again).
 *    6. Combine the differences with weights (-1/λ_ℓ)^{k-1} and accumulate
 *       into CHANGE_STAT[ℓ].
 *    7. Temporarily apply the actual toggle (TOGGLE_IF_MORE_TO_COME) so that
 *       later toggles see updated degrees/adjs.
 *
 *  At the end, UNDO_PREVIOUS_TOGGLES restores the original network.
 */
D_CHANGESTAT_FN(d_cov_diff_GW){

#if DEBUG_COV_DIFF_GW
  static int seen = 0;
  if(ntoggles > 1 && seen < 10){
    Rprintf("[cov_diff_GW] MULTI-TOGGLE ntoggles=%d\n", (int)ntoggles);
    seen++;
  }
#endif

  /* 1) Reset output buffer for the whole proposal. */
  ZERO_ALL_CHANGESTATS();

  /* 2) Read inputs from INPUT_PARAM. */
  const double *ip     = INPUT_PARAM;
  const int n1         = (int)ip[0];   /* number of actors in the actor mode */
  const int L          = (int)ip[1];   /* number of lambda values */
  const double *lambda = ip + 2;       /* lambda[0..L-1] */
  const double *x      = ip + 2 + L;   /* actor covariate values */

#if DEBUG_COV_DIFF_GW
  Rprintf("[cov_diff_GW] n1=%d L=%d | lambda[0]=%g\n",
          n1, L, (L > 0 ? lambda[0] : NA_REAL));
#endif

  /* 3) Allocate reusable buffers once (per proposal). */
  const int ck_len = n1 + 1;
  double *ck_before = (double *)R_Calloc(ck_len, double);
  double *ck_after  = (double *)R_Calloc(ck_len, double);
  int    *idxs      = (int    *)R_Calloc(n1,     int);

  /* 4) Process toggles sequentially. */
  int i = 0;
  FOR_EACH_TOGGLE(i){

    Vertex t = TAIL(i);
    Vertex h = HEAD(i);

    /* Identify actor and group vertices for the current toggle. */
    Vertex actor = (t <= (Vertex)n1) ? t : h;
    Vertex group = (t <= (Vertex)n1) ? h : t;

#if DEBUG_COV_DIFF_GW
    if(!(actor <= (Vertex)n1 && group > (Vertex)n1)){
      Rprintf("[cov_diff_GW][WARN] toggle #%d not (actor,group): tail=%d head=%d (n1=%d)\n",
              i, (int)t, (int)h, n1);
    }
#endif

    /* Clear ck buffers (defensive: avoids stale values for k > ng). */
    memset(ck_before, 0, (size_t)ck_len * sizeof(double));
    memset(ck_after,  0, (size_t)ck_len * sizeof(double));

    int ng_before = 0;
    int ng_after  = 0;

    /* c_k(g) before the toggle (current intermediate state). */
    group_covdiff_allk(group, n1, x, nwp, &ng_before, ck_before, ck_len, idxs);

    /* Virtual toggle on this dyad (compute "after"). */
    TOGGLE(t, h);
    group_covdiff_allk(group, n1, x, nwp, &ng_after, ck_after, ck_len, idxs);
    TOGGLE(t, h); /* undo virtual toggle */

    /* Local maximum group size for this group. */
    int Kmax = (ng_before > ng_after) ? ng_before : ng_after;
    if(Kmax >= 2){

      /* For each lambda_ℓ, compute Δ T_GW(λ_ℓ). */
      for(int l = 0; l < L; l++){
        double lam  = lambda[l];
        double base = -1.0 / lam;  /* weight base for k = 2 */
        double pw   = base;        /* (-1/lam)^(k-1) with k=2 => power 1 */
        double deltaT = 0.0;

        for(int k = 2; k <= Kmax; k++){
          double dck = ck_after[k] - ck_before[k];
          if(dck != 0.0){
            deltaT += pw * dck;
          }
          pw *= base;
        }

        CHANGE_STAT[l] += deltaT;

#if DEBUG_COV_DIFF_GW
        Rprintf("[D:cov_diff_GW] i=%d tail=%d head=%d | group=%d | ng=%d->%d | lambda[%d]=%g | Δ=%g | cumul=%g\n",
                i, (int)t, (int)h, (int)group, ng_before, ng_after, l+1, lam, deltaT, CHANGE_STAT[l]);
#endif
      }
    }

    /* Apply the actual toggle temporarily so later toggles see updated state. */
    TOGGLE_IF_MORE_TO_COME(i);
  }

  /* 5) Undo temporary toggles to restore the original network state. */
  UNDO_PREVIOUS_TOGGLES(i);

  /* 6) Release buffers. */
  R_Free(ck_before);
  R_Free(ck_after);
  R_Free(idxs);
}
