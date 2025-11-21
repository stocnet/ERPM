/**
 * @file
 * @brief  Change statistic for the ERPM term `cov_match` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `cov_match`, which counts monochromatic k-cliques of actors inside groups
 *  based on a categorical covariate, with several normalisation modes and an
 *  optional target category.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  A bipartite network is assumed, with:
 *    - actor mode  = actors carrying a categorical covariate c(i),
 *    - group mode  = structural groups, each connecting to multiple actors.
 *
 *  For each group g, let:
 *    - A(g)      = set of actors connected to g (actor membership),
 *    - n_g       = |A(g)| = group size,
 *    - c(i)      = category code of actor i,
 *    - n_{g,r}   = number of actors in group g with category r.
 *
 *  The basic combinatorial quantity for a given k ≥ 1 is:
 *
 *      S_k(B; c) = ∑_g ∑_r C(n_{g,r}, k),
 *
 *  where C(a,b)=CHOOSE(a,b) is the binomial coefficient with the convention
 *  C(a,b) = 0 if a < b.
 *
 *  Targeted version (category κ):
 *
 *      S_k^{(κ)}(B; c) = ∑_g C(n_{g,κ}, k).
 *
 *  Normalisation modes:
 *
 *   - "none" (norm_mode = 0):
 *        statistic = S_k(B; c)              (or S_k^{(κ)} in targeted mode).
 *
 *   - "by_group" (norm_mode = 1):
 *        non-target:
 *          ∑_g [ ( ∑_r C(n_{g,r}, k) ) / C(n_g, k) ],
 *        targeted:
 *          ∑_g [ C(n_{g,κ}, k) / C(n_g, k) ].
 *
 *   - "global" (norm_mode = 2):
 *        non-target:
 *          (1 / C(N_actors, k)) * S_k(B; c),
 *        targeted:
 *          (1 / C(N_actors, k)) * S_k^{(κ)}(B; c),
 *
 *  where N_actors = n1 is the number of actors.
 *
 *  Special case for k = 1, "by_group", targeted:
 *    The statistic reduces to:
 *
 *      ∑_g 1[ n_{g,κ} ≥ 1 ],
 *
 *    i.e. a count of groups that contain at least one actor in category κ.
 *
 *  ------------------------------------------------------------
 *  One-toggle implementation in {ergm}
 *  ------------------------------------------------------------
 *
 *  For each membership toggle between an actor and a group:
 *
 *    1. Identify:
 *         - the actor vertex in the actor mode,
 *         - the group vertex in the group mode.
 *
 *    2. Compute the current membership of the group:
 *         - list all actor neighbours of the group (incoming and outgoing),
 *         - deduplicate via a temporary bitmap,
 *         - count group size n_g and build an empirical histogram of
 *           category codes for that group.
 *
 *    3. From this local information, read:
 *         - n_{g,r*} for the actor involved in the toggle (r* = c(actor)),
 *         - optionally n_{g,κ} for the target category κ if requested.
 *
 *    4. Use binomial identities to compute the unnormalised local change:
 *         - addition of an actor with category r*:
 *             Δ_non_norm = + C(n_{g,r*}^{old}, k-1),
 *         - removal:
 *             Δ_non_norm = - C(n_{g,r*}^{new}-1, k-1),
 *           with new = old - 1.
 *
 *    5. Apply the requested normalisation:
 *         - none:        Δ = Δ_non_norm,
 *         - by_group:    Δ = (N_plus/D_plus) - (N_minus/D_minus),
 *         - global:      Δ = Δ_non_norm / C(N_actors, k).
 *
 *    6. Update:
 *
 *         CHANGE_STAT[j] += Δ
 *
 *       for each k = ks[j] in the vectorised set of k values.
 *
 *  Only the group touched by the membership toggle is recomputed. No other
 *  groups are inspected.
 *
 *  ------------------------------------------------------------
 *  INPUT_PARAM layout
 *  ------------------------------------------------------------
 *
 *  The R initialiser packs the parameters as:
 *
 *      INPUT_PARAM = c(
 *        n1,
 *        K,
 *        norm_mode,
 *        has_kappa,
 *        kappa_code,
 *        ks[1:K],
 *        z_codes[1:n1]
 *      )
 *
 *  where:
 *    - n1         = number of actors (actor mode),
 *    - K          = number of vectorised k values,
 *    - norm_mode  = 0 (none), 1 (by_group), 2 (global),
 *    - has_kappa  = 0 (no target category) or 1 (targeted),
 *    - kappa_code = integer code of the targeted category if has_kappa = 1,
 *    - ks         = integer vector of k ≥ 1,
 *    - z_codes    = integer codes of the categorical covariate on actors:
 *                     z_codes[i] = category code of actor i, 0 if NA.
 *
 *  At C level:
 *
 *    P[0]          = n1
 *    P[1]          = K
 *    P[2]          = norm_mode
 *    P[3]          = has_kappa
 *    P[4]          = kappa_code
 *    P[5 .. 5+K-1] = ks[0 .. K-1]
 *    P[5+K .. ]    = z_codes[0 .. n1-1]
 *
 *  N_CHANGE_STATS is equal to K, one component per k.
 *
 *  ------------------------------------------------------------
 *  Complexity and trade-offs
 *  ------------------------------------------------------------
 *
 *  For each toggle:
 *    - building the actor neighbourhood for the touched group is O(deg(group)),
 *    - building the histogram over distinct category codes is O(deg(group)),
 *    - the per-k updates reuse these counts and are O(K + m), where m is the
 *      number of distinct categories present in the group.
 *
 *  No state is cached across toggles. This is sufficient for small/medium-sized
 *  groups and suitable as a first implementation; further optimisation could
 *  store group-level summaries for reuse across toggles.
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  The R side (InitErgmTerm.cov_match):
 *    - validates the actor covariate and normalisation options,
 *    - encodes categories as integer codes z_codes,
 *    - constructs the vector ks of k values and sets N_CHANGE_STATS = length(ks),
 *    - builds INPUT_PARAM as above,
 *    - sets emptynwstats = 0 and appropriate coef.names.
 *
 *  ------------------------------------------------------------
 *  @example Usage (R)
 *  ------------------------------------------------------------
 *  @code{.r}
 *  library(ERPM)
 *
 *  # Example partition of 6 actors into 3 groups
 *  part <- c(1, 1, 2, 2, 3, 3)
 *
 *  # Categorical covariate on actors (e.g., colors)
 *  color <- c("red", "red", "blue", "blue", "red", "blue")
 *
 *  # Non-normalised cov_match with k = 2
 *  fit1 <- erpm(partition ~ cov_match(k = 2, attr = color))
 *  summary(fit1)
 *
 *  # Targeted and group-normalised version:
 *  #   counts, for each group, the proportion of pairs (i,j) with color == "red"
 *  fit2 <- erpm(partition ~ cov_match(k = 2,
 *                                     attr      = color,
 *                                     category  = "red",
 *                                     normalize = "by_group"))
 *  summary(fit2)
 *
 *  # Internally, each membership toggle between an actor and a group calls
 *  # c_cov_match() once, which recomputes the local category counts in the
 *  # affected group, applies the appropriate combinatorial formulas and
 *  # normalisation, then updates CHANGE_STAT[] for each requested k.
 *  @endcode
 */

#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <math.h>
#include <string.h>

/**
 * @def DEBUG_COV_MATCH
 * @brief Enable verbose debugging output for ::c_cov_match.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - actor and group participating in the toggle,
 *  - group size and category counts,
 *  - intermediate unnormalised and normalised deltas for each k.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_COV_MATCH 0   /* set to 0 to disable debug output */

/**
 * @def UNUSED_WARNING
 * @brief Mark a parameter as intentionally unused.
 *
 * @param x Identifier of the unused variable.
 *
 * This macro is used to silence compiler warnings when a parameter is
 * required by the interface but not accessed in the implementation.
 */
#define UNUSED_WARNING(x) (void)x

/* -------------------------------------------------------------------------- */
/* Helper: category code lookup                                               */
/* -------------------------------------------------------------------------- */

/**
 * @brief Retrieve the integer category code for a given actor.
 *
 * @details
 *  The array @p z_codes is indexed in the same way as actor vertices:
 *    - actor vertices are 1-based (1..n1),
 *    - z_codes is stored as a double array and cast to int on access.
 *  A value ≤ 0 is interpreted as “no category” (e.g., NA).
 *
 * @param i        Actor vertex index (1..n1).
 * @param z_codes  Pointer to the array of category codes (length n1).
 *
 * @return Integer category code for actor i, or 0 if undefined.
 */
static inline int code_of_actor(Vertex i, const double *z_codes){
  /* z_codes is 1-indexed in the same way as actor vertices (actors 1..n1). */
  return (int)z_codes[(size_t)(i-1)];
}

/* -------------------------------------------------------------------------- */
/* Helper: neighbours of a group in the actor mode                            */
/* -------------------------------------------------------------------------- */

/**
 * @brief Collect unique actor neighbours of a group vertex.
 *
 * @details
 *  For a group vertex @p g in the group mode, this function:
 *    - walks through outgoing edges from g and records any neighbour h that
 *      lies in the actor mode (h ≤ n1),
 *    - walks through incoming edges to g and records any neighbour h that
 *      lies in the actor mode,
 *    - deduplicates actor neighbours using a temporary bitmap @p seen,
 *    - stores the unique actor vertices in the provided buffer @p actors.
 *
 *  The buffer @p actors must be large enough to hold up to n1 vertices.
 *
 * @param nwp     Pointer to the {ergm} Network structure.
 * @param g       Group vertex whose actor neighbours are queried.
 * @param actors  Output buffer that will receive the actor vertex indices.
 * @param n1      Number of actors (and maximum number of neighbours).
 *
 * @return The number of unique actor neighbours written into @p actors.
 */
static int neighbors_actors_of_group(Network *nwp, Vertex g, Vertex *actors, int n1){
  int cnt = 0;
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char); /* Calloc -> 0 */

  Vertex h; Edge e;

  /* OUT-neighbours: group -> actor. */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        actors[cnt++] = h;
      }
    }
  }

  /* IN-neighbours: actor -> group. */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        actors[cnt++] = h;
      }
    }
  }

  Free(seen);
  return cnt;
}

/* -------------------------------------------------------------------------- */
/* Helper: histogram of category codes in a group                             */
/* -------------------------------------------------------------------------- */

/**
 * @brief Build a histogram of category codes for the actors of a group.
 *
 * @details
 *  Given a list of actor vertices in @p actors, this function:
 *    - looks up the integer category code for each actor via ::code_of_actor,
 *    - ignores codes ≤ 0 (e.g., NA),
 *    - builds arrays @p codes and @p counts of distinct categories and their
 *      frequencies, in no particular order.
 *
 *  The arrays @p codes and @p counts must be large enough to hold all distinct
 *  categories that can appear in the group (up to @p na).
 *
 * @param actors   Array of actor vertices belonging to the group.
 * @param na       Number of actors in @p actors.
 * @param z_codes  Pointer to the actor category code array.
 * @param codes    Output array receiving distinct category codes.
 * @param counts   Output array receiving counts per category code.
 *
 * @return The number of distinct categories written into @p codes and @p counts.
 */
static int histogram_codes(const Vertex *actors, int na, const double *z_codes, int *codes, int *counts){
  int m = 0; /* number of distinct categories seen */
  for(int a=0; a<na; ++a){
    int code = code_of_actor(actors[a], z_codes);
    if(code<=0) continue; /* ignore NA / undefined */
    int found = 0;
    for(int j=0; j<m; ++j){
      if(codes[j]==code){ counts[j]++; found=1; break; }
    }
    if(!found){
      codes[m]=code; counts[m]=1; m++;
    }
  }
  return m;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_match (one-toggle)                                   */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_match`.
 *
 * @details
 *  This is the {ergm} change-statistic function registered as
 *  ::c_cov_match via ::C_CHANGESTAT_FN. It implements the one-toggle
 *  update for the cov_match effect, which counts groups that are
 *  homogeneous (or partially homogeneous) with respect to a categorical
 *  actor covariate, possibly with:
 *
 *    - multiple values of k vectorised,
 *    - an optional target category,
 *    - three normalisation modes: none, by_group, global.
 *
 *  For a membership toggle between an actor and a group:
 *
 *    1. The function reads global parameters from INPUT_PARAM:
 *         - n1, K, norm_mode, has_kappa, kappa_code,
 *         - ks[0..K-1] (the vector of k values),
 *         - actor category codes z_codes[0..n1-1].
 *
 *    2. It identifies:
 *         - v_actor in the actor mode (index ≤ n1),
 *         - v_group in the group mode (index > n1).
 *
 *    3. It builds the list of current actor neighbours of v_group and
 *       constructs a histogram of category codes in that group.
 *
 *    4. It extracts:
 *         - n_g_old  = current group size,
 *         - n_{g,r*} = count of actors with the same category as v_actor,
 *         - n_{g,κ}  = count of actors in the target category κ (if any).
 *
 *    5. For each k = ks[j], it computes the unnormalised local change
 *       Δ_non_norm using binomial identities, then applies:
 *         - no normalisation (norm_mode = 0),
 *         - group-level normalisation (norm_mode = 1),
 *         - global normalisation (norm_mode = 2).
 *
 *    6. The result is accumulated in CHANGE_STAT[j].
 *
 *  The parameter @p edgestate indicates whether the edge currently exists:
 *    - edgestate = 0  → the edge is absent (toggle = addition),
 *    - edgestate = 1  → the edge is present (toggle = deletion).
 *
 *  In this implementation, @p edgestate is used only to infer the sign and
 *  direction of the local change, and not to explicitly toggle the network.
 *
 * @param tail       Tail vertex of the toggled edge (actor or group).
 * @param head       Head vertex of the toggled edge (actor or group).
 * @param mtp        Pointer to the model term structure (unused directly
 *                   here but required by the macro).
 * @param nwp        Pointer to the network-plus workspace (used to access
 *                   adjacency information and degrees).
 * @param edgestate  Current state of the edge:
 *                   - 0 if the edge is absent (toggle = addition),
 *                   - 1 if the edge is present (toggle = deletion).
 */
C_CHANGESTAT_FN(c_cov_match){
  /* 1) Reset the output buffer for THIS toggle.
   *
   * {ergm} accumulates contributions from multiple calls. Here we only
   * provide the local Δ for the current membership toggle.
   */
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  /* 2) Read packed inputs from R.
   *
   * Layout:
   *   P[0]   = n1
   *   P[1]   = K
   *   P[2]   = norm_mode (0 none, 1 by_group, 2 global)
   *   P[3]   = has_kappa (0/1)
   *   P[4]   = kappa_code
   *   P[5..] = ks[0..K-1]
   *   P[5+K..] = z_codes[0..n1-1]
   */
  const double *P = INPUT_PARAM;

  const int n1          = (int)P[0];
  const int K           = (int)P[1];
  const int norm_mode   = (int)P[2];  /* 0 none, 1 by_group, 2 global */
  const int has_kappa   = (int)P[3];  /* 0 or 1 */
  const int kappa_code  = (int)P[4];

  const double *ks_d    = P + 5;
  const double *z_codes = P + 5 + K;  /* length n1 */

  #if DEBUG_COV_MATCH
    Rprintf("[cov_match:init] n1=%d K=%d norm=%d has_kappa=%d kappa_code=%d ks{",
            n1, K, norm_mode, has_kappa, kappa_code);
    for(int jj=0;jj<K;++jj){ Rprintf("%s%d", (jj?",":""), (int)ks_d[jj]); }
    Rprintf("}\n");
  #endif

  /* 3) Identify actor and group vertices for this toggle.
   *
   * Exactly one endpoint lies in the actor mode (≤ n1), the other in the
   * group mode (> n1).
   */
  Vertex t = tail, h = head;
  Vertex v_actor = (t <= n1) ? t : h;
  Vertex v_group = (t >  n1) ? t : h;

  /* Minimal sanity checks: ensure we truly have one actor and one group. */
  if(v_actor<=0 || v_actor>n1) return;
  if(v_group<=n1) return;

  /* 4) Determine the nature of the toggle.
   *
   *   edgestate = 1 → edge exists → toggle = deletion,
   *   edgestate = 0 → edge absent → toggle = addition.
   */
  const int is_add = edgestate ? 0 : 1;

  #if DEBUG_COV_MATCH
    Rprintf("[cov_match:toggle] actor=%d group=%d is_add=%d OUT_DEG[g]=%d IN_DEG[g]=%d\n",
            (int)v_actor, (int)v_group, is_add, (int)OUT_DEG[v_group], (int)IN_DEG[v_group]);
  #endif

  /* 5) Build local information for the group affected by the toggle.
   *
   *  - actors_buf: actor members of the group,
   *  - codes_buf, counts_buf: category histogram for those actors.
   */
  Vertex actors_buf[8192];
  int    codes_buf[ 8192];
  int    counts_buf[8192];

  int maxbuf = (n1 < 8192) ? n1 : 8192;
  int na = neighbors_actors_of_group(nwp, v_group, actors_buf, n1);
  if(na>maxbuf) na = maxbuf; /* soft safety bound */

  int m = histogram_codes(actors_buf, na, z_codes, codes_buf, counts_buf);

  /* 6) Extract the category code and count for the actor in the toggle. */
  const int r_star = code_of_actor(v_actor, z_codes);
  int n_gr_old = 0;
  for(int j=0;j<m;++j) if(codes_buf[j]==r_star) { n_gr_old = counts_buf[j]; break; }

  const int n_g_old = na;

  #if DEBUG_COV_MATCH
    Rprintf("[cov_match:group] n_g_old=%d r*=%d n_gr_old=%d m=%d\n",
            n_g_old, r_star, n_gr_old, m);
    if(has_kappa){
      int n_gk_old_dbg = 0;
      if(kappa_code>0){
        for(int u=0; u<m; ++u) if(codes_buf[u]==kappa_code){ n_gk_old_dbg = counts_buf[u]; break; }
      }
      Rprintf("[cov_match:group] kappa=%d n_gk_old=%d\n", kappa_code, n_gk_old_dbg);
    }
  #endif

  /* 7) Main loop over vectorised k values. */
  for(int j=0; j<K; ++j){
    const int k = (int)ks_d[j];

    /* Special case:
     *   k == 1, norm_mode == by_group, targeted(κ):
     *   statistic = ∑_g 1[ n_{g,κ} ≥ 1 ].
     */
    if(norm_mode==1 && has_kappa && k==1){
      int n_gk_old = 0;
      if(kappa_code>0){
        for(int u=0; u<m; ++u) if(codes_buf[u]==kappa_code){ n_gk_old = counts_buf[u]; break; }
      }
      int n_gk_new = n_gk_old;
      if(r_star==kappa_code){
        n_gk_new += (is_add ? +1 : -1);
        if(n_gk_new < 0) n_gk_new = 0;
      }
      double delta_ind = (n_gk_new>0) - (n_gk_old>0);
      CHANGE_STAT[j] += delta_ind;

      #if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=1][by_group][kappa] n_g_old=%d n_gk_old=%d n_gk_new=%d delta=%g\n",
                n_g_old, n_gk_old, n_gk_new, delta_ind);
      #endif
      continue;
    }

    /* 8) Compute unnormalised delta for this k. */
    double delta_non_norm = 0.0;

    if(has_kappa){
      /* Targeted version: only counts in category κ matter. */
      int n_gk_old = 0;
      if(kappa_code>0){
        for(int u=0; u<m; ++u) if(codes_buf[u]==kappa_code){ n_gk_old = counts_buf[u]; break; }
      }
      if(is_add){
        delta_non_norm = (r_star==kappa_code) ? CHOOSE(n_gk_old, k-1) : 0.0;
      }else{
        delta_non_norm = (r_star==kappa_code) ? -CHOOSE((n_gk_old-1), k-1) : 0.0;
      }
      #if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=%d][kappa=%d] n_gk_old=%d delta_non_norm=%.6f\n",
                k, kappa_code, n_gk_old, delta_non_norm);
      #endif
    }else{
      /* Non-targeted version: use counts for the actor's own category r*. */
      if(is_add){
        delta_non_norm = CHOOSE(n_gr_old, k-1);
      }else{
        delta_non_norm = -CHOOSE((n_gr_old-1), k-1);
      }
      #if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=%d][non-targeted] n_gr_old=%d delta_non_norm=%.6f\n",
                k, n_gr_old, delta_non_norm);
      #endif
    }

    /* 9) Apply normalisation according to norm_mode. */
    double delta = delta_non_norm;

    if(norm_mode==1){
      /* by_group normalisation:
       *
       *   Δ = (N_plus / D_plus) - (N_minus / D_minus),
       *
       * where:
       *   - N_minus is either ∑_r C(n_{g,r}, k) or C(n_{g,κ}, k),
       *   - D_minus = C(n_g_old, k),
       *   - N_plus and D_plus are the same quantities after the toggle.
       */
      double N_minus = 0.0;
      if(has_kappa){
        int n_gk_old = 0;
        if(kappa_code>0){
          for(int u=0; u<m; ++u) if(codes_buf[u]==kappa_code){ n_gk_old = counts_buf[u]; break; }
        }
        N_minus = CHOOSE(n_gk_old, k);
      }else{
        for(int u=0; u<m; ++u) N_minus += CHOOSE(counts_buf[u], k);
      }
      double D_minus = CHOOSE(n_g_old, k);

      double N_plus = N_minus + delta_non_norm;
      double D_plus = CHOOSE(n_g_old + (is_add?+1:-1), k);

      double ratio_minus = (D_minus>0.0) ? (N_minus/D_minus) : 0.0;
      double ratio_plus  = (D_plus >0.0) ? (N_plus /D_plus ) : 0.0;

      delta = ratio_plus - ratio_minus;

      #if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=%d][by_group] N-=%g D-=%g | N+=%g D+=%g | delta=%.6f\n",
                k, N_minus, D_minus, N_plus, D_plus, delta);
      #endif
    }
    else if(norm_mode==2){
      /* global normalisation:
       *
       *   Δ = Δ_non_norm / C(N_actors, k).
       */
      double Cglob = CHOOSE(n1, k);
      delta = (Cglob>0.0) ? (delta_non_norm / Cglob) : 0.0;
      #if DEBUG_COV_MATCH
        Rprintf("[cov_match:k=%d][global] Cglob=%g delta_non_norm=%.6f -> delta=%.6f\n",
                k, Cglob, delta_non_norm, delta);
      #endif
    }

    /* 10) Accumulate the contribution for this k. */
    CHANGE_STAT[j] += delta;
  }
}
