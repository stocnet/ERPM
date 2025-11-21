/**
 * @file
 * @brief  Change statistic for the ERPM term `cov_match_GW` (one-toggle form).
 *
 * @details
 *  This file implements the {ergm} change statistic for the ERPM effect
 *  `cov_match_GW`, which applies a geometrically weighted transform to
 *  group-level category counts for a categorical actor covariate.
 *
 *  ------------------------------------------------------------
 *  Statistical principle (actor mode, group mode)
 *  ------------------------------------------------------------
 *
 *  The network is bipartite:
 *    - actor mode  = vertices representing actors, each carrying a categorical
 *                    covariate c(i),
 *    - group mode  = vertices representing structural groups.
 *
 *  For each group vertex g in the group mode:
 *    - A(g)      = set of actors connected to g,
 *    - n_g       = |A(g)| = size of group g,
 *    - c(i)      = category code of actor i,
 *    - n_{g,r}   = number of actors in g with category r.
 *
 *  For a given decay parameter λ ≥ 1, define:
 *
 *      r_λ = (λ - 1) / λ  ∈ [0, 1).
 *
 *  A generic geometrically weighted transform at group level can be written as:
 *
 *      S_g(λ; c) = ∑_r f_λ(n_{g,r}),
 *
 *  where f_λ(·) is a function of the category counts n_{g,r} that is
 *  chosen so that local changes can be expressed using simple powers of r_λ.
 *
 *  In this implementation, the ERPM effect is encoded so that a membership
 *  toggle of an actor with category r* produces a local contribution for the
 *  targeted count n_{g,r*} of the form:
 *
 *    - addition (n_{g,r*} = m → m+1):
 *        Δ_non_norm =  r_λ^m,
 *    - deletion (n_{g,r*} = m → m-1):
 *        Δ_non_norm = -r_λ^{m-1}.
 *
 *  When a target category κ is specified, the transform is restricted to the
 *  single count n_{g,κ}. Otherwise, it aggregates over all categories present
 *  in the group.
 *
 *  Normalisation modes:
 *
 *    - "none" (norm_mode = 0):
 *        statistic ∝ ∑_g ∑_r GW(n_{g,r}; λ),
 *        local Δ given directly by Δ_non_norm.
 *
 *    - "by_group" (norm_mode = 1):
 *        For each group g, build:
 *
 *          Num(g) = ∑_r λ (1 - r_λ^{n_{g,r}})      (or λ (1 - r_λ^{n_{g,κ}})
 *                                                 in targeted version),
 *          Den(g) = λ (1 - r_λ^{n_g}),
 *
 *        then the contribution of g is:
 *
 *          Num(g) / Den(g),
 *
 *        so for a toggle in group g, the local change is:
 *
 *          Δ = [Num(g)_after / Den(g)_after] - [Num(g)_before / Den(g)_before].
 *
 *    - "global" (norm_mode = 2):
 *        A global normalisation rescales the non-normalised statistic by
 *        a group-size independent constant:
 *
 *          Δ = Δ_non_norm / [λ (1 - r_λ^{N_actors})],
 *
 *        where N_actors = n1 is the number of actors.
 *
 *  Only the group touched by the membership toggle is recomputed. All other
 *  groups are unaffected.
 *
 *  ------------------------------------------------------------
 *  Bipartite structure and actors/groups
 *  ------------------------------------------------------------
 *
 *  At the C level, the bipartite structure is encoded as:
 *    - the first n1 vertices are actors (actor mode),
 *    - the remaining vertices are groups (group mode).
 *
 *  Each membership toggle always connects:
 *    - exactly one actor vertex (index in 1..n1),
 *    - exactly one group vertex (index > n1).
 *
 *  This convention is used to split the dyad (tail, head) into:
 *    - v_actor in the actor mode,
 *    - v_group in the group mode.
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
 *        lambdas[1:K],
 *        z_codes[1:n1]
 *      )
 *
 *  where:
 *    - n1          = number of actors,
 *    - K           = number of λ values (vectorised),
 *    - norm_mode   = 0 (none), 1 (by_group), 2 (global),
 *    - has_kappa   = 0 (no target category) or 1 (targeted),
 *    - kappa_code  = integer code of the targeted category (if has_kappa = 1),
 *    - lambdas     = vector of λ_j ≥ 1 (length K),
 *    - z_codes     = integer codes of the categorical covariate on actors:
 *                      z_codes[i] = category code of actor i, 0 if NA.
 *
 *  At C level:
 *
 *    P[0]          = n1
 *    P[1]          = K
 *    P[2]          = norm_mode
 *    P[3]          = has_kappa
 *    P[4]          = kappa_code
 *    P[5 .. 5+K-1] = lambdas[0 .. K-1]
 *    P[5+K .. ]    = z_codes[0 .. n1-1]
 *
 *  N_CHANGE_STATS is equal to K, one statistic per λ_j.
 *
 *  ------------------------------------------------------------
 *  Complexity
 *  ------------------------------------------------------------
 *
 *  For each membership toggle:
 *    - building the neighbour list for the touched group is O(deg(group)),
 *    - building the histogram of category codes is O(deg(group)),
 *    - the per-λ updates use only scalar operations (O(K)).
 *
 *  No state is cached across toggles. This is sufficient for small and
 *  medium-sized groups and can be optimised later if needed.
 *
 *  ------------------------------------------------------------
 *  R interface
 *  ------------------------------------------------------------
 *
 *  The R initialiser (InitErgmTerm.cov_match_GW):
 *    - validates λ > 1 and the actor-level covariate,
 *    - encodes the actor categories as integer codes z_codes,
 *    - packs the λ values into INPUT_PARAM,
 *    - sets emptynwstats = 0 and vectorised coef.names,
 *    - ensures N_CHANGE_STATS = length(lambdas).
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
 *  # Categorical covariate on actors
 *  color <- c("red", "red", "blue", "blue", "red", "blue")
 *
 *  # Geometrically weighted cov_match over all categories
 *  fit1 <- erpm(partition ~ cov_match_GW(lambda = 2, attr = color))
 *  summary(fit1)
 *
 *  # Targeted and group-normalised version:
 *  #   emphasises groups where the share of "red" actors is high, with
 *  #   a geometric down-weighting controlled by lambda.
 *  fit2 <- erpm(partition ~ cov_match_GW(lambda    = 2,
 *                                        attr      = color,
 *                                        category  = "red",
 *                                        normalize = "by_group"))
 *  summary(fit2)
 *
 *  # Internally, each membership toggle between an actor and a group
 *  # calls c_cov_match_GW() once, which:
 *  #   - reconstructs the category counts in the touched group,
 *  #   - computes the local Δ for each requested lambda,
 *  #   - applies the chosen normalisation mode,
 *  #   - updates the corresponding components of CHANGE_STAT[].
 *  @endcode
 */

#include <R_ext/Print.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include <math.h>
#include <string.h>

/**
 * @def DEBUG_COV_MATCH_GW
 * @brief Enable verbose debugging output for ::c_cov_match_GW.
 *
 * Set this macro to 1 to print diagnostic information to the R console
 * during `summary()` or MCMC runs:
 *  - actor and group vertices for the current toggle,
 *  - group size and category counts before the toggle,
 *  - intermediate non-normalised and normalised deltas per λ.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_COV_MATCH_GW 0

/**
 * @def UNUSED_WARNING
 * @brief Mark a variable as intentionally unused.
 *
 * @param x Identifier of the unused variable.
 *
 * This macro is used to silence compiler warnings when a parameter or
 * variable is required by the interface but not accessed in the code.
 */
#define UNUSED_WARNING(x) (void)x

/* -------------------------------------------------------------------------- */
/* Helper: category code lookup                                               */
/* -------------------------------------------------------------------------- */

/**
 * @brief Retrieve the integer category code for a given actor.
 *
 * @details
 *  The array @p z_codes is indexed consistently with actor vertices:
 *    - actor vertices are 1-based (1..n1),
 *    - z_codes is stored as a double array and cast to int on access.
 *  Values ≤ 0 are interpreted as “no category” (e.g., NA).
 *
 * @param i        Actor vertex index (1..n1).
 * @param z_codes  Pointer to the array of category codes (length n1).
 *
 * @return Integer category code for actor @p i, or 0 if undefined.
 */
static inline int code_of_actor(Vertex i, const double *z_codes){
  /* Actors are indexed 1..n1; convert to 0-based index in z_codes. */
  return (int)z_codes[(size_t)(i-1)]; // actors indexed 1..n1
}

/* -------------------------------------------------------------------------- */
/* Helper: neighbours of a group in the actor mode                            */
/* -------------------------------------------------------------------------- */

/**
 * @brief Collect unique actor neighbours of a given group vertex.
 *
 * @details
 *  For a group vertex @p g in the group mode, this function:
 *    - walks through outgoing edges from g and records actor neighbours,
 *    - walks through incoming edges to g and records actor neighbours,
 *    - deduplicates actor neighbours using a temporary bitmap @p seen,
 *    - writes the resulting list of actors into the @p actors buffer.
 *
 *  Only neighbours whose vertex index is ≤ n1 are considered actors.
 *
 * @param nwp     Pointer to the {ergm} Network structure.
 * @param g       Group vertex whose actor neighbours are queried.
 * @param actors  Output buffer that will receive actor vertex indices.
 * @param n1      Number of actors (size of the actor mode).
 *
 * @return The number of unique actor neighbours stored in @p actors.
 */
static int neighbors_actors_of_group(Network *nwp, Vertex g, Vertex *actors, int n1){
  int cnt = 0;
  /* Temporary bitmap indicating whether an actor index has been seen. */
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char); // 0-inited
  Vertex h;
  Edge e;

  /* Traverse outgoing edges from g (group -> actor). */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx]=1;
        actors[cnt++]=h;
      }
    }
  }

  /* Traverse incoming edges to g (actor -> group). */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx]=1;
        actors[cnt++]=h;
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
 * @brief Build a histogram of category codes for actors in a group.
 *
 * @details
 *  Given an array of actor vertices belonging to a group, this function:
 *    - looks up the integer category code for each actor via ::code_of_actor,
 *    - ignores codes ≤ 0 (e.g., NA),
 *    - accumulates counts per distinct category into @p codes and @p counts.
 *
 *  The arrays @p codes and @p counts must be large enough to hold all
 *  distinct categories that may appear (up to @p na).
 *
 * @param actors   Array of actor vertices in the group.
 * @param na       Number of actors in @p actors.
 * @param z_codes  Pointer to the actor category code array (length n1).
 * @param codes    Output array for distinct category codes.
 * @param counts   Output array for counts per category code.
 *
 * @return The number of distinct categories written into @p codes and @p counts.
 */
static int histogram_codes(const Vertex *actors, int na, const double *z_codes, int *codes, int *counts){
  int m = 0;
  for(int a=0; a<na; ++a){
    int code = code_of_actor(actors[a], z_codes);
    if(code<=0) continue; // ignore NA / undefined
    int found = 0;
    for(int j=0; j<m; ++j){
      if(codes[j]==code){
        counts[j]++;
        found=1;
        break;
      }
    }
    if(!found){
      codes[m]=code;
      counts[m]=1;
      m++;
    }
  }
  return m;
}

/* -------------------------------------------------------------------------- */
/* Change statistic: cov_match_GW (one-toggle)                                */
/* -------------------------------------------------------------------------- */

/**
 * @brief Change statistic for the ERPM term `cov_match_GW`.
 *
 * @details
 *  This is the {ergm} change-statistic function registered as
 *  ::c_cov_match_GW via ::C_CHANGESTAT_FN. It implements the one-toggle
 *  update for the geometrically weighted covariate-matching effect, with:
 *
 *    - potentially multiple λ values vectorised (one per statistic),
 *    - an optional targeted category κ,
 *    - three normalisation modes (none, by_group, global).
 *
 *  For each membership toggle between an actor and a group:
 *
 *    1. The function reads global parameters from INPUT_PARAM:
 *         - n1, K, norm_mode, has_kappa, kappa_code,
 *         - lambdas[0..K-1],
 *         - z_codes[0..n1-1] (actor categories).
 *
 *    2. It identifies:
 *         - v_actor as the endpoint in the actor mode (vertex index ≤ n1),
 *         - v_group as the endpoint in the group mode (vertex index > n1).
 *
 *    3. It constructs the membership of v_group in the actor mode by:
 *         - traversing incoming and outgoing edges touching v_group,
 *         - deduplicating actor neighbours,
 *         - building a local histogram of categories in the group.
 *
 *    4. It extracts:
 *         - n_g_old  = group size (number of actors in v_group),
 *         - n_{g,r*} = count of actors with the same category as v_actor,
 *         - n_{g,κ}  = count of actors with category κ (if has_kappa = 1).
 *
 *    5. For each λ_j in the vector:
 *         - computes r_λj = (λ_j - 1) / λ_j,
 *         - builds a non-normalised local delta Δ_non_norm as:
 *             * targeted:  based on n_{g,κ},
 *             * non-targeted: based on n_{g,r*},
 *           with addition/removal rules:
 *             - addition:  Δ_non_norm = r_λ^{m},
 *             - deletion:  Δ_non_norm = -r_λ^{m-1},
 *           where m is the relevant count before the toggle.
 *
 *    6. Applies normalisation:
 *         - none:      Δ = Δ_non_norm,
 *         - by_group:  Δ = (Num_after/Den_after) - (Num_before/Den_before),
 *         - global:    Δ = Δ_non_norm / [λ (1 - r_λ^{n1})].
 *
 *    7. Accumulates the result:
 *
 *          CHANGE_STAT[j] += Δ
 *
 *       for the j-th component corresponding to λ_j.
 *
 *  The parameter @p edgestate indicates whether the edge currently exists:
 *    - edgestate = 0 → the edge is absent (toggle = addition),
 *    - edgestate = 1 → the edge is present (toggle = deletion).
 *
 *  In this implementation, the network is not explicitly toggled; the
 *  computation is based on the current state and on the sign inferred
 *  from @p edgestate.
 *
 * @param tail       Tail vertex of the toggled edge (actor or group).
 * @param head       Head vertex of the toggled edge (actor or group).
 * @param mtp        Pointer to the model term parameters (unused directly
 *                   but required by the macro).
 * @param nwp        Pointer to the network-plus workspace, providing access
 *                   to adjacency and degree information.
 * @param edgestate  Current state of the edge:
 *                   - 0 if the edge is absent (toggle = addition),
 *                   - 1 if the edge is present (toggle = deletion).
 */
C_CHANGESTAT_FN(c_cov_match_GW){
  /* 1) Reset the output buffer for THIS toggle.
   *
   * {ergm} accumulates change-statistics across multiple toggles, but this
   * function computes the local Δ only for the current membership toggle.
   */
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  /* 2) Read packed parameters from INPUT_PARAM.
   *
   * Layout:
   *   P[0]   = n1
   *   P[1]   = K
   *   P[2]   = norm_mode (0 none, 1 by_group, 2 global)
   *   P[3]   = has_kappa (0/1)
   *   P[4]   = kappa_code
   *   P[5..] = lambdas[0..K-1]
   *   P[5+K..] = z_codes[0..n1-1]
   */
  const double *P = INPUT_PARAM;

  const int n1         = (int)P[0];
  const int K          = (int)P[1];
  const int norm_mode  = (int)P[2]; // 0 none, 1 by_group, 2 global
  const int has_kappa  = (int)P[3]; // 0/1
  const int kappa_code = (int)P[4];

  const double *lambdas = P + 5;
  const double *z_codes = P + 5 + K;

  /* 3) Identify the actor and group vertices involved in the toggle.
   *
   * Exactly one endpoint is in the actor mode (vertex index ≤ n1),
   * the other is in the group mode (vertex index > n1).
   */
  Vertex t = tail, h = head;
  Vertex v_actor = (t <= n1) ? t : h;
  Vertex v_group = (t >  n1) ? t : h;

  /* Minimal sanity checks: ensure we truly have an actor and a group. */
  if(v_actor<=0 || v_actor>n1) return;
  if(v_group<=n1) return;

  /* 4) Determine whether this toggle is an addition or a deletion.
   *
   *   edgestate = 1 → edge exists → toggle = deletion,
   *   edgestate = 0 → edge absent → toggle = addition.
   */
  const int is_add = edgestate ? 0 : 1; // 1=addition, 0=deletion

  #if DEBUG_COV_MATCH_GW
    Rprintf("[cov_match_GW:toggle] actor=%d group=%d is_add=%d OUT_DEG[g]=%d IN_DEG[g]=%d\n",
            (int)v_actor, (int)v_group, is_add,
            (int)OUT_DEG[v_group], (int)IN_DEG[v_group]);
  #endif

  /* 5) Build local information for the group affected by the toggle.
   *
   *  - actors_buf: list of actor neighbours of v_group,
   *  - codes_buf, counts_buf: histogram of category codes in that group.
   */
  int maxbuf = n1 < 8192 ? n1 : 8192;
  Vertex actors_buf[8192];
  int    codes_buf[8192];
  int    counts_buf[8192];

  int na = neighbors_actors_of_group(nwp, v_group, actors_buf, n1);
  if(na>maxbuf) na = maxbuf; /* soft cap to avoid overflow */

  int m = histogram_codes(actors_buf, na, z_codes, codes_buf, counts_buf);

  /* 6) Extract the category code and count for the actor in the toggle. */
  const int r_star = code_of_actor(v_actor, z_codes);
  int n_gr_old = 0;
  for(int j=0;j<m;++j){
    if(codes_buf[j]==r_star){
      n_gr_old = counts_buf[j];
      break;
    }
  }
  const int n_g_old = na;

  /* 7) Pre-compute n_{g,κ} if a target category κ is requested. */
  int n_gk_old = 0;
  if(has_kappa && kappa_code>0){
    for(int u=0; u<m; ++u){
      if(codes_buf[u]==kappa_code){
        n_gk_old = counts_buf[u];
        break;
      }
    }
  }

  #if DEBUG_COV_MATCH_GW
    Rprintf("[cov_match_GW:group] n_g_old=%d r*=%d n_gr_old=%d m=%d n_gk_old=%d (kappa=%d)\n",
            n_g_old, r_star, n_gr_old, m, n_gk_old, kappa_code);
  #endif

  /* 8) Main loop over λ values (vectorised statistics). */
  for(int j=0; j<K; ++j){
    const double lambda = lambdas[j];
    const double rlam   = (lambda - 1.0) / lambda; // r_λ in [0,1)
    double delta_non_norm = 0.0;

    /* 8.a Compute the non-normalised local change Δ_non_norm. */
    if(has_kappa){
      /* Targeted version: only the count n_{g,κ} contributes. */
      if(r_star==kappa_code){
        if(is_add){
          /* addition: m = n_gk_old → Δ = r_λ^{m} */
          delta_non_norm = pow(rlam, (double)n_gk_old);
        }else{
          /* deletion: m = n_gk_old → Δ = -r_λ^{m-1} (m≥1 if edge exists) */
          delta_non_norm = -pow(rlam, (double)(n_gk_old-1));
        }
      }else{
        /* Actor category is not κ → no contribution. */
        delta_non_norm = 0.0;
      }
    }else{
      /* Non-targeted: use the count for the actor's own category r*. */
      if(is_add){
        delta_non_norm = pow(rlam, (double)n_gr_old);
      }else{
        delta_non_norm = -pow(rlam, (double)(n_gr_old-1));
      }
    }

    double delta = delta_non_norm;

    /* 8.b Apply normalisation, if requested. */
    if(norm_mode==1){
      /* Group-level normalisation:
       *
       *   Δ = (Num_after / Den_after) - (Num_before / Den_before),
       *
       * where:
       *   Num(g) = ∑_r λ (1 - r_λ^{n_{g,r}})    (or λ (1 - r_λ^{n_{g,κ}})
       *                                             in targeted version),
       *   Den(g) = λ (1 - r_λ^{n_g}).
       */

      /* Compute Num_before and Den_before from local histogram. */
      double N_minus = 0.0;
      if(has_kappa){
        N_minus = lambda * (1.0 - pow(rlam, (double)n_gk_old));
      }else{
        for(int u=0; u<m; ++u){
          N_minus += lambda * (1.0 - pow(rlam, (double)counts_buf[u]));
        }
      }
      double D_minus = lambda * (1.0 - pow(rlam, (double)n_g_old));

      /* Group size after the toggle. */
      const int n_g_new  = n_g_old + (is_add ? +1 : -1);

      /* Update Num and Den locally to obtain Num_after, Den_after. */
      double N_plus = N_minus;
      if(has_kappa){
        /* Only the κ-cell changes. */
        if(r_star==kappa_code){
          if(is_add){
            N_plus += lambda * pow(rlam,(double)n_gk_old);
          }else{
            N_plus -= lambda * pow(rlam,(double)(n_gk_old-1));
          }
        }
      }else{
        /* Only the cell for r* changes. */
        if(is_add){
          N_plus += lambda * pow(rlam,(double)n_gr_old);
        }else{
          N_plus -= lambda * pow(rlam,(double)(n_gr_old-1));
        }
      }
      double D_plus = lambda * (1.0 - pow(rlam, (double)n_g_new));

      double ratio_minus = (D_minus>0.0) ? (N_minus/D_minus) : 0.0;
      double ratio_plus  = (D_plus >0.0) ? (N_plus /D_plus ) : 0.0;

      delta = ratio_plus - ratio_minus;

      #if DEBUG_COV_MATCH_GW
        Rprintf("[cov_match_GW][by_group] j=%d lambda=%g N-=%g D-=%g N+=%g D+=%g delta=%g\n",
                j, lambda, N_minus, D_minus, N_plus, D_plus, delta);
      #endif

    }else if(norm_mode==2){
      /* Global normalisation:
       *
       *   Δ = Δ_non_norm / [λ (1 - r_λ^{N_actors})],
       *
       * with N_actors = n1.
       */
      const double Dglob = lambda * (1.0 - pow(rlam, (double)n1));
      delta = (Dglob>0.0) ? (delta_non_norm / Dglob) : 0.0;

      #if DEBUG_COV_MATCH_GW
        Rprintf("[cov_match_GW][global] j=%d lambda=%g Dglob=%g delta_non_norm=%g delta=%g\n",
                j, lambda, Dglob, delta_non_norm, delta);
      #endif
    }

    /* 9) Accumulate the contribution for the j-th statistic. */
    CHANGE_STAT[j] += delta;
  }
}
