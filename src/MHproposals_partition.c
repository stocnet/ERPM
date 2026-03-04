/**
 * @file MHproposals_partition.c
 * @brief Metropolis-Hastings proposals for ERPM partition networks under '~b1part'.
 *
 * ERPM represents a partition as a bipartite 'membership' network:
 * - Actor-mode vertices are  1..BIPARTITE (n actors).
 * - Group-mode vertices are (BIPARTITE+1)..N_NODES (G groups, padded).
 * - Under '~b1part', each actor has exactly one incident membership edge.
 *
 * ERPM convention:
 * - The membership network must be undirected. If a directed network is detected,
 *   proposals stop with a hard error.
 *
 * Implemented proposals (legacy 'ergm' MH API):
 * - MH_ErpmToggleStep: move one actor from its current group to another group (2 toggles).
 * - MH_ErpmSwapStep  : swap memberships of two actors in different groups (4 toggles).
 * - MH_ErpmMix       : mixed proposal that selects ToggleStep/SwapStep by user weights.
 *
 * Legacy MH API notes (ergm_MHproposal.h):
 * - Proposals are implemented as MH_P_FN(MH_<Name>) and do not return a status.
 * - Initialization is detected when MHp->ntoggles == 0, and the proposal must set
 *   its fixed toggle count (2 or 4 here).
 * - A valid proposal must fill:
 *   MHp->ntoggles, Mtail[...], Mhead[...], and MHp->logratio.
 *
 * IMPORTANT (mixed proposal):
 * - MH_ErpmToggleStep and MH_ErpmSwapStep use an init branch:
 *     if(MHp->ntoggles==0){ MHp->ntoggles=fixed; return; }
 *   Therefore MH_ErpmMix MUST NOT call them directly, or the first dispatch would
 *   produce a no-move.
 * - MH_ErpmMix dispatches to per-iteration helpers:
 *     ErpmToggleStep_propose() / ErpmSwapStep_*propose()
 *   which do not have an init branch and always emit a proposal (or report "impossible").
 *
 * ---------------------------------------------------------------------------
 * Partition representation and invariants
 * ---------------------------------------------------------------------------
 * A partition of n actors is encoded as an undirected bipartite membership graph:
 *
 * - Actors: vertices 1..BIPARTITE.
 * - Groups: vertices (BIPARTITE+1)..N_NODES (padded group set).
 * - Membership: an edge (actor, group) indicates the actor's group assignment.
 *
 * Under '~b1part', the intended invariant is:
 * - each actor has degree exactly 1 (one membership edge).
 *
 * The sampler assumes '~b1part' is active during MCMC but still performs
 * defensive checks:
 * - directed membership networks are rejected (fatal error);
 * - corrupted states (actor without group) trigger MH_FAILED or error().
 *
 * Padded groups (degree 0) are allowed and can become active through moves.
 *
 * ---------------------------------------------------------------------------
 * Proposal semantics
 * ---------------------------------------------------------------------------
 *
 * 1) ToggleStep (2 toggles)
 * -------------------------
 * - Select an actor uniformly.
 * - Let group_old be its current group.
 * - Select group_new uniformly among all other group vertices
 *   (including padded groups).
 *
 * Proposal encoding:
 *   (actor, group_old) OFF
 *   (actor, group_new) ON
 *
 * Hastings ratio:
 * - symmetric proposal:
 *     actor uniform in both directions,
 *     target group uniform over constant set (G_total−1).
 * - therefore logratio = 0.
 *
 * 2) SwapStep (4 toggles)
 * -----------------------
 * - Requires at least two non-empty groups.
 * - Draw actor_i uniformly and read group_i.
 * - Draw actor_j uniformly until:
 *       actor_j != actor_i  AND  group_j != group_i.
 *
 * Proposal encoding:
 *   (actor_i, group_i) OFF → (actor_i, group_j) ON
 *   (actor_j, group_j) OFF → (actor_j, group_i) ON
 *
 * Hastings ratio:
 * - identical selection rule before and after swap,
 *   hence symmetric and logratio = 0.
 *
 * Structural property:
 * - SwapStep preserves group sizes exactly, so it only mixes states
 *   within the same size vector.
 *
 * If all actors belong to a single group, no swap is possible.
 *
 * 3) Mixed kernel (ErpmMix)
 * -------------------------
 * - Randomly selects ToggleStep or SwapStep according to user weights.
 * - If SwapStep is selected but impossible (e.g. only one non-empty group),
 *   the kernel falls back to ToggleStep instead of reporting MH_FAILED.
 *
 * Storage:
 * - persistent arrays stored in MHp->storage:
 *     move_codes[k] and cumprob[k].
 * - memory managed with R_Calloc / R_Free.
 *
 * ---------------------------------------------------------------------------
 * Interaction with ergm internals
 * ---------------------------------------------------------------------------
 * Uses ergm network macros and structures:
 * - BIPARTITE, N_NODES, IN_DEG[...] (reference only for SwapStep logic).
 * - IS_UNDIRECTED_EDGE(u,v) to test membership edges.
 *
 * Undirected requirement is checked through nwp->directed_flag.
 *
 * Failure behaviour:
 * - directed network → error().
 * - invalid '~b1part' state → MH_FAILED or error().
 * - no legal SwapStep:
 *     standalone SwapStep → MH_FAILED/error().
 *     inside Mix → fallback ToggleStep.
 *
 * ---------------------------------------------------------------------------
 * Debugging
 * ---------------------------------------------------------------------------
 * - DEBUG_ERPM_PROPOSALS enables verbose traces.
 * - DBG_max_print limits console output during long MCMC runs.
 * - Debugging is intended for development/selftests and should remain
 *   disabled in production builds.
 */

/**
 * @file MHproposals_partition.c
 * @brief Metropolis-Hastings proposals for ERPM partition networks under '~b1part'.
 *
 * ERPM represents a partition as a bipartite 'membership' network:
 * - Actor-mode vertices are  1..BIPARTITE (n actors).
 * - Group-mode vertices are (BIPARTITE+1)..N_NODES (G groups, padded).
 * - Under '~b1part', each actor has exactly one incident membership edge.
 *
 * ERPM convention:
 * - The membership network must be undirected. If a directed network is detected,
 *   proposals stop with a hard error.
 *
 * Implemented proposals (legacy 'ergm' MH API):
 * - MH_ErpmToggleStep: move one actor from its current group to another group (2 toggles).
 * - MH_ErpmSwapStep  : swap memberships of two actors in different groups (4 toggles).
 * - MH_ErpmMix       : mixed proposal that selects ToggleStep/SwapStep by user weights.
 *
 * Legacy MH API notes (ergm_MHproposal.h):
 * - Proposals are implemented as MH_P_FN(MH_<Name>) and do not return a status.
 * - Initialization is detected when MHp->ntoggles == 0, and the proposal must set
 *   its fixed toggle count (2 or 4 here).
 * - A valid proposal must fill:
 *   MHp->ntoggles, Mtail[...], Mhead[...], and MHp->logratio.
 *
 * IMPORTANT (mixed proposal):
 * - MH_ErpmToggleStep and MH_ErpmSwapStep use an init branch:
 *     if(MHp->ntoggles==0){ MHp->ntoggles=fixed; return; }
 *   Therefore MH_ErpmMix MUST NOT call them directly, or the first dispatch would
 *   produce a no-move.
 * - MH_ErpmMix dispatches to per-iteration helpers:
 *     ErpmToggleStep_propose() / ErpmSwapStep_*propose()
 *   which do not have an init branch and always emit a proposal (or report "impossible").
 *
 * ---------------------------------------------------------------------------
 * Partition representation and invariants
 * ---------------------------------------------------------------------------
 * ERPM encodes a partition of n actors into (padded) groups as an undirected,
 * bipartite membership graph.
 *
 * - Actor mode:
 *     vertices 1..BIPARTITE (n actors).
 * - Group mode:
 *     vertices (BIPARTITE+1)..N_NODES (G_total groups, padded).
 * - Membership:
 *     an undirected edge (actor, group) means "actor is assigned to group".
 *
 * Under the '~b1part' constraint, the intended invariant is:
 * - Each actor has degree exactly 1 (one and only one incident membership edge).
 *
 * This file assumes that '~b1part' (and the implied '~b1degrees' machinery inside
 * ergm) is active during MCMC. Nevertheless, defensive checks are used:
 * - A directed membership network is treated as fatal.
 * - A corrupted state (actor with no group) triggers MH_FAILED or a hard error,
 *   depending on whether MH_FAILED is available at compile time.
 *
 * Note on padding:
 * - Group vertices may include empty "padded" groups (degree 0).
 * - ToggleStep is allowed to target padded groups (this can create new non-empty
 *   groups and/or destroy singleton groups by moving their only actor away).
 *
 * ---------------------------------------------------------------------------
 * Proposal semantics (high-level)
 * ---------------------------------------------------------------------------
 * 1) ToggleStep (2 toggles)
 * ------------------------
 * - Choose an actor uniformly among all actors.
 * - Let group_old be its unique current group.
 * - Choose group_new uniformly among ALL group vertices except group_old,
 *   including padded empty groups.
 * - Propose reassignment actor: group_old -> group_new.
 *
 * Encoding:
 * - Toggle (actor, group_old) OFF
 * - Toggle (actor, group_new) ON
 *
 * Hastings correction:
 * - With this selection rule, the proposal is symmetric because:
 *     - actor is uniform in both directions;
 *     - group_new is uniform over a constant-size set (G_total - 1), independent
 *       of the number of non-empty groups;
 *   hence MHp->logratio = 0.
 *
 * 2) SwapStep (4 toggles) — WITHOUT IN_DEG
 * ---------------------------------------
 * - Choose two distinct actors uniformly, with the constraint that they currently
 *   belong to two different groups.
 * - Swap their memberships between these two groups.
 *
 * Operationally (UPDATED):
 * - Pre-check: if fewer than 2 non-empty groups exist, no legal swap exists.
 * - Draw actor_i uniformly in 1..BIPARTITE; read group_i.
 * - Draw actor_j uniformly in 1..BIPARTITE until:
 *     actor_j != actor_i AND group_j != group_i
 *   (NO bounded retry cap; see bias note below).
 *
 * Encoding:
 * - Toggle (actor_i, group_i) OFF, (actor_i, group_j) ON
 * - Toggle (actor_j, group_j) OFF, (actor_j, group_i) ON
 *
 * Hastings correction:
 * - With the same selection rule before and after the swap (uniform actors with the
 *   "different groups" constraint), the proposal is symmetric, hence logratio = 0.
 *
 * Bias note (why we avoid max_tries):
 * - A hard retry cap (e.g., 200) can incorrectly classify a rare-but-possible swap as
 *   "impossible" in highly unbalanced partitions, which subtly breaks exact symmetry.
 * - We therefore (i) pre-check P>=2, and (ii) sample until success without a small cap.
 *
 * Structural note:
 * - SwapStep preserves group sizes exactly, so it cannot connect states with different
 *   size vectors. It is intended as an intra-size mixing move.
 *
 * Impossibility:
 * - If all actors are in the same group, no legal swap exists.
 * - Standalone MH_ErpmSwapStep reports MH_FAILED (or error()) in that case.
 *
 * 3) Mix (ToggleStep / SwapStep) with SWAP fallback
 * ------------------------------------------------
 * - The mixed proposal samples which move to apply using user-provided weights.
 * - If the selected move is SWAP but no legal swap exists (e.g., only one non-empty
 *   group), Mix MUST FALL BACK to ToggleStep and MUST NOT propagate MH_FAILED for
 *   this specific reason (to keep the chain moving under mixed kernels).
 *
 * Storage:
 * - The mix uses persistent storage (MHp->storage) containing:
 *     - move_codes[k] (int)
 *     - cumprob[k] (double)
 * - Allocation uses R_Calloc/R_Free to match R's memory management conventions
 *   for long-lived proposal storage.
 *
 * ---------------------------------------------------------------------------
 * Interaction with ergm internals
 * ---------------------------------------------------------------------------
 * - This file uses ergm's Network representation and macros:
 *     - BIPARTITE, N_NODES, IN_DEG[...] from ergm_changestat.h (IN_DEG no longer
 *       used by SwapStep, but may still be used elsewhere / in reference code)
 *     - IS_UNDIRECTED_EDGE(u,v) to test membership edges
 * - The check "undirected network required" uses nwp->directed_flag because the
 *   DIRECTED macro from ergm_changestat_common.do_not_include_directly.h is not
 *   visible from this compilation unit.
 *
 * Failure behavior:
 * - Fatal configuration error:
 *     - directed membership network => error()
 * - Invalid/corrupted '~b1part' state:
 *     - missing membership edge for chosen actor => MH_FAILED (if available)
 *       or error().
 * - No legal SwapStep:
 *     - Standalone SwapStep => MH_FAILED (if available) or error().
 *     - Mix when SWAP selected => fallback ToggleStep (never MH_FAILED for this case).
 *
 * ---------------------------------------------------------------------------
 * Debugging
 * ---------------------------------------------------------------------------
 * - Compile-time flag DEBUG_ERPM_PROPOSALS controls verbose traces.
 * - A print limiter (DBG_max_print) prevents flooding the console during long
 *   MCMC runs.
 * - Debug output is intended for development and selftests; it should remain
 *   disabled by default in production builds. (UPDATED: default is now OFF)
 */

#include "ergm_MHproposals_degree.h"
#include "ergm_MHproposal.h" /* MH_FAILED */
#include "ergm_changestat.h"

#include <R_ext/Print.h>    /* Rprintf */
#include <R_ext/Error.h>    /* error  */
#include <R_ext/Memory.h>   /* R_Calloc, R_Free, R_alloc */
#include <math.h>           /* floor, log, isfinite */
#include <Rmath.h>          /* unif_rand */


/* -------------------------------------------------------------------------- */
/* Debugging                                                                  */
/* -------------------------------------------------------------------------- */
/* Allows print along the proposals execution */
#define DEBUG_ERPM_PROPOSALS 0
#define UNUSED_VARIABLE(x) (void)(x)

/* Print limiter (avoid flooding console during long MCMC). */
#if DEBUG_ERPM_PROPOSALS
  static int DBG_seen_toggle = 0;
  static int DBG_seen_swap   = 0;
  static int DBG_seen_mix    = 0;
  static const int DBG_max_print = 400;
#endif

#if DEBUG_ERPM_PROPOSALS
static inline void DBG_print_context_header(const char *who, Network *nwp){
  const Vertex n1 = BIPARTITE;
  const Vertex N  = N_NODES;
  const Vertex G  = N - n1;
  Rprintf("[ERPM][%s][CTX] n1=%d N=%d G_total=%d directed_flag=%d\n",
          who, (int)n1, (int)N, (int)G, (int)nwp->directed_flag);
}

static inline Vertex DBG_count_nonempty_groups(Network *nwp){
  Vertex P = 0;
  for(Vertex g = nwp->bipartite + 1; g <= nwp->nnodes; ++g){
    if(nwp->indegree[g] > 0) P++;
  }
  return P;
}
#endif

/* -------------------------------------------------------------------------- */
/* Guards                                                                     */
/* -------------------------------------------------------------------------- */
static inline void erpm_require_undirected(Network *nwp){
  if(nwp->directed_flag){
    error("ERPM partition proposals require an undirected membership network (DIRECTED=TRUE detected).");
  }
}

/* -------------------------------------------------------------------------- */
/* RNG helper                                                                 */
/* -------------------------------------------------------------------------- */
static inline Vertex UnifVertex1to(Vertex upper_bound){
  if(upper_bound < 1){
    error("ERPM proposals: UnifVertex1to() called with upper_bound < 1.");
  }
  return (Vertex)(1 + (Vertex)floor(unif_rand() * upper_bound));
}

/* -------------------------------------------------------------------------- */
/* Group counts (non-empty groups)                                            */
/* -------------------------------------------------------------------------- */
/* NOTE:
 * We use nwp->indegree[g] as a fast proxy for group size.
 * In an undirected bipartite membership network, a group vertex's indegree equals
 * its incident membership edge count (up to ergm's internal conventions).
 */
static inline Vertex erpm_count_nonempty_groups(Network *nwp){
  Vertex P = 0;
  for(Vertex g = BIPARTITE + 1; g <= N_NODES; ++g){
    if(nwp->indegree[g] > 0) ++P;
  }
  return P;
}

/* -------------------------------------------------------------------------- */
/* Membership decoding                                                        */
/* -------------------------------------------------------------------------- */
/* Slow-path decoder (O(G_total) scan). Used as a correctness baseline.
 * Performance-critical paths (SwapStep) can decode actor->group for all actors
 * once per propose to avoid repeated scans.
 */
static Vertex get_current_group_of_actor(Vertex actor_id, Network *nwp){
  erpm_require_undirected(nwp);

  const Vertex number_of_actors = BIPARTITE;
  if(actor_id < 1 || actor_id > number_of_actors) return 0;

  Vertex found_group = 0;

  for(Vertex group_id = number_of_actors + 1; group_id <= N_NODES; ++group_id){
    if(IS_UNDIRECTED_EDGE(actor_id, group_id)){
      if(found_group != 0){
#if DEBUG_ERPM_PROPOSALS
        if(DBG_seen_toggle < DBG_max_print){
          Rprintf("[ERPM][WARN] actor=%d has multiple group neighbors (at least %d and %d)\n",
                  (int)actor_id, (int)found_group, (int)group_id);
        }
#endif
        return found_group;
      }
      found_group = group_id;
    }
  }

#if DEBUG_ERPM_PROPOSALS
  if(found_group == 0 && DBG_seen_toggle < DBG_max_print){
    Rprintf("[ERPM][WARN] actor=%d has no group neighbor (invalid b1part state?)\n", (int)actor_id);
  }
#endif

  return found_group;
}

/* Decode actor->group for ALL actors in the current state.
 * Returns 1 on success, 0 if a corrupted b1part state is detected.
 *
 * Memory:
 * - Caller is expected to allocate ag[0..BIPARTITE] (Vertex array) with R_alloc
 *   or persistent storage. Index 0 is unused.
 */
static int erpm_decode_actor_groups(Vertex *ag, Network *nwp){
  for(Vertex a = 1; a <= BIPARTITE; ++a){
    Vertex g = get_current_group_of_actor(a, nwp);
    if(g == 0) return 0;
    ag[a] = g;
  }
  return 1;
}

/* -------------------------------------------------------------------------- */
/* MH_B1Part (reference from ergm)                                            */
/* -------------------------------------------------------------------------- */
MH_P_FN(MH_B1Part) {
  if (MHp->ntoggles == 0) { /* Initialize */
    MH_CondB1Degree(MHp, nwp);
    return;
  }

  MH_CondB1Degree(MHp, nwp);

  int dP = (IN_DEG[Mhead[1]] == 0) - (IN_DEG[Mhead[0]] == 1);

  if (dP) {
    Vertex P = 0;
    for (Vertex i = BIPARTITE + 1; i <= N_NODES; i++) {
      if (IN_DEG[i]) P++;
    }

    MHp->logratio += dP == -1 ? log(N_NODES - BIPARTITE - P + 1) : -log(N_NODES - BIPARTITE - P);
  }
}

/* -------------------------------------------------------------------------- */
/* Per-iteration helpers (no init branch)                                     */
/* -------------------------------------------------------------------------- */
static void ErpmToggleStep_propose(MHProposal *MHp, Network *nwp){
  erpm_require_undirected(nwp);

  const Vertex number_of_actors       = BIPARTITE;
  const Vertex number_of_groups_total = N_NODES - BIPARTITE;

  Vertex actor_id  = UnifVertex1to(number_of_actors);
  Vertex group_old = get_current_group_of_actor(actor_id, nwp);

  if(group_old == 0){
// We keep a MH_FAILED guard to be robust to ergm changes (to be remove later)
#ifdef MH_FAILED 
  #if DEBUG_ERPM_PROPOSALS
      if(DBG_seen_toggle < DBG_max_print){
        DBG_print_context_header("ToggleStep", nwp);
        Vertex P = DBG_count_nonempty_groups(nwp);
        Rprintf("[ERPM][ToggleStep][FAIL] MH_FAILED (actor has no group)\n");
        Rprintf("[ERPM][ToggleStep][FAIL] actor=%d | nonempty_groups(P)=%d\n",
                (int)actor_id, (int)P);
        DBG_seen_toggle++;
      }
  #endif
      MHp->ntoggles = MH_FAILED;
      MHp->logratio = 0.0;
  #else
      error("ERPM ToggleStep: invalid b1part state (actor has no group neighbor).");
#endif
    return;
  }

  // Draw group_new uniformly from ALL groups except group_old, without a retry loop. 
  if(number_of_groups_total < 2){
#ifdef MH_FAILED
    MHp->ntoggles = MH_FAILED;
    MHp->logratio = 0.0;
#else
    error("ERPM ToggleStep: need at least 2 group vertices to propose a move.");
#endif
    return;
  }

  const Vertex gold_rank = group_old - BIPARTITE;                /* 1..G_total */
  const Vertex r = UnifVertex1to(number_of_groups_total - 1);    /* 1..G_total-1 */
  const Vertex gnew_rank = (r >= gold_rank) ? (r + 1) : r;       /* skip old rank */
  const Vertex group_new = BIPARTITE + gnew_rank;

  MHp->ntoggles = 2;
  Mtail[0] = actor_id; Mhead[0] = group_old;
  Mtail[1] = actor_id; Mhead[1] = group_new;
  MHp->logratio = 0.0;
}

/**
 * Swap proposer (no init branch), returns whether a legal swap was emitted.
 *
 * Return value:
 * - 1 : success (ntoggles=4, toggles filled)
 * - 0 : no legal swap exists in the current state (e.g., all actors in same group)
 * - -1: corrupted b1part state encountered (caller may MH_FAILED or hard error)
 *
 * Rationale:
 * - Standalone MH_ErpmSwapStep should still signal failure when impossible.
 * - Mixed MH_ErpmMix wants a SWAP->TOGGLE fallback for the "impossible" case,
 *   without ever going through MH_FAILED for that reason.
 *
 * Bias/Perf:
 * - We first check that at least two non-empty groups exist (P>=2).
 * - We avoid bounded retry caps (max_tries) because they can misclassify a rare-but-possible
 *   swap as "impossible" in unbalanced partitions, breaking exact symmetry.
 * - We decode actor->group for all actors once per propose (O(n)) to avoid repeated O(G) scans.
 */
static int ErpmSwapStep_trypropose(MHProposal *MHp, Network *nwp){
  erpm_require_undirected(nwp);

  const Vertex number_of_actors = BIPARTITE;

  /* Fast impossibility pre-check: if fewer than 2 non-empty groups exist, no legal swap exists. */
  if(erpm_count_nonempty_groups(nwp) < 2){
    return 0;
  }

  /* Decode actor->group once to avoid repeated scans in the sampling loop. */
  Vertex *ag = (Vertex*) R_alloc((size_t)(number_of_actors + 1), sizeof(Vertex)); /* 1..n */
  if(!erpm_decode_actor_groups(ag, nwp)){
    return -1; /* corrupted state */
  }

  /* 1) Sample actor_i uniformly and read its group. */
  Vertex actor_i = UnifVertex1to(number_of_actors);
  Vertex group_i = ag[actor_i];
  if(group_i == 0) return -1;

  /* 2) Sample actor_j uniformly until it is in a different group (no small cap). */
  Vertex actor_j = 0;
  Vertex group_j = 0;

  for(;;){
    Vertex cand_actor = UnifVertex1to(number_of_actors);
    if(cand_actor == actor_i) continue;

    Vertex cand_group = ag[cand_actor];
    if(cand_group == 0) return -1;

    if(cand_group != group_i){
      actor_j = cand_actor;
      group_j = cand_group;
      break;
    }
  }

  /* 3) Encode swap (4 toggles). */
  MHp->ntoggles = 4;

  Mtail[0] = actor_i; Mhead[0] = group_i;
  Mtail[1] = actor_i; Mhead[1] = group_j;

  Mtail[2] = actor_j; Mhead[2] = group_j;
  Mtail[3] = actor_j; Mhead[3] = group_i;

  MHp->logratio = 0.0;
  return 1;
}

/* Backward-compatible wrapper: strict SwapStep behavior (standalone). */
static void ErpmSwapStep_propose_strict(MHProposal *MHp, Network *nwp){
  const int rc = ErpmSwapStep_trypropose(MHp, nwp);

  if(rc == 1) return;

  if(rc == 0){
#ifdef MH_FAILED
#if DEBUG_ERPM_PROPOSALS
    if(DBG_seen_swap < DBG_max_print){
      DBG_print_context_header("SwapStep", nwp);
      Rprintf("[ERPM][SwapStep][FAIL] MH_FAILED (no legal swap: <2 non-empty groups)\n");
      DBG_seen_swap++;
    }
#endif
    MHp->ntoggles = MH_FAILED;
    MHp->logratio = 0.0;
#else
    error("ERPM SwapStep: no legal swap exists (need at least two non-empty groups).");
#endif
    return;
  }

//rc == -1 : corrupted state 
// We keep a MH_FAILED guard to be robust to ergm changes (to be remove later)
#ifdef MH_FAILED
  #if DEBUG_ERPM_PROPOSALS
    if(DBG_seen_swap < DBG_max_print){
      DBG_print_context_header("SwapStep", nwp);
      Rprintf("[ERPM][SwapStep][FAIL] MH_FAILED (invalid b1part state encountered)\n");
      DBG_seen_swap++;
    }
  #endif
    MHp->ntoggles = MH_FAILED;
    MHp->logratio = 0.0;
#else
  error("ERPM SwapStep: invalid b1part state (actor has no group neighbor).");
#endif
}

/* -------------------------------------------------------------------------- */
/* ERPM Mix (ToggleStep / SwapStep)                                           */
/* -------------------------------------------------------------------------- */
#define ERPM_MOVE_TOGGLE 1
#define ERPM_MOVE_SWAP   2

typedef struct ErpmMixStorage {
  int K;
  int    *move_codes;
  double *cumprob;
} ErpmMixStorage;

static void erpm_mix_build_cumprob(ErpmMixStorage *st, MHProposal *MHp){
  if(!st)  error("ERPM ErpmMix: internal error (NULL storage).");
  if(!MHp) error("ERPM ErpmMix: internal error (NULL MHp).");

  if(MHp->iinputs == NULL)
    error("ERPM ErpmMix: missing iinputs (expected K + move codes).");

  const int K = (int)MHp->iinputs[0];
  if(K <= 0)
    error("ERPM ErpmMix: invalid K in iinputs[0] (must be > 0).");

  if(MHp->inputs == NULL)
    error("ERPM ErpmMix: missing inputs (expected weights).");

  st->K = K;
  st->move_codes = (int*)    R_Calloc((size_t)K, int);
  st->cumprob    = (double*) R_Calloc((size_t)K, double);

  double wsum = 0.0;

  for(int k = 0; k < K; ++k){
    const int code = (int)MHp->iinputs[1 + k];
    st->move_codes[k] = code;

    const double w = MHp->inputs[k];
    if(!(isfinite(w) && w > 0.0))
      error("ERPM ErpmMix: all weights must be finite and > 0.");

    if(code != ERPM_MOVE_TOGGLE && code != ERPM_MOVE_SWAP)
      error("ERPM ErpmMix: unknown move code in iinputs (expected TOGGLE/SWAP).");

    wsum += w;
  }

  double c = 0.0;
  for(int k = 0; k < K; ++k){
    c += MHp->inputs[k] / wsum;
    st->cumprob[k] = c;
  }
  st->cumprob[K - 1] = 1.0;
}

static int erpm_mix_sample_move(const ErpmMixStorage *st){
  const double u = unif_rand();
  for(int k = 0; k < st->K; ++k){
    if(u < st->cumprob[k]) return st->move_codes[k];
  }
  return st->move_codes[st->K - 1];
}

static void erpm_mix_free_storage(MHProposal *MHp){
  ErpmMixStorage *st = (ErpmMixStorage*) MHp->storage;
  if(!st) return;

  if(st->move_codes) R_Free(st->move_codes);
  if(st->cumprob)    R_Free(st->cumprob);
  R_Free(st);

  MHp->storage = NULL;
}

MH_I_FN(Mi_ErpmMix){
#if DEBUG_ERPM_PROPOSALS
  Rprintf("[ERPM][Mix][INIT] setting ntoggles=4 (max toggles)\n");
#endif

  ErpmMixStorage *st = (ErpmMixStorage*) R_Calloc(1, ErpmMixStorage);
  erpm_mix_build_cumprob(st, MHp);
  MHp->storage = st;

  /* Max toggles among component moves: SWAP is 4. */
  MHp->ntoggles = 4;
  MHp->logratio = 0.0;
}

MH_P_FN(MH_ErpmMix){
  erpm_require_undirected(nwp);

  /* Mi_ErpmMix must have run. Lazy init is treated as a wiring bug. */
  ErpmMixStorage *st = (ErpmMixStorage*) MHp->storage;
  if(st == NULL){
    error("ERPM ErpmMix: NULL storage (Mi_ErpmMix was not called).");
  }

  int move = erpm_mix_sample_move(st);

  if(move == ERPM_MOVE_TOGGLE){
    ErpmToggleStep_propose(MHp, nwp);
  } else if(move == ERPM_MOVE_SWAP){
    /* Try SWAP. If impossible, FALL BACK to TOGGLE and do NOT use MH_FAILED. */
    const int rc = ErpmSwapStep_trypropose(MHp, nwp);

    if(rc == 1){
      /* OK: swap emitted */
    } else if(rc == 0){
#if DEBUG_ERPM_PROPOSALS
      if(DBG_seen_mix < DBG_max_print){
        Rprintf("[ERPM][Mix][FALLBACK] SWAP impossible -> fallback TOGGLE\n");
        DBG_seen_mix++;
      }
#endif
      move = ERPM_MOVE_TOGGLE;  
      ErpmToggleStep_propose(MHp, nwp);
    } else {
      // rc == -1: corrupted state -> keep the usual failure behavior 
// We keep a MH_FAILED guard to be robust to ergm changes (to be remove later)
#ifdef MH_FAILED
  #if DEBUG_ERPM_PROPOSALS
        if(DBG_seen_mix < DBG_max_print){
          DBG_print_context_header("Mix", nwp);
          Rprintf("[ERPM][Mix][FAIL] SWAP encountered invalid b1part state -> MH_FAILED\n");
          DBG_seen_mix++;
        }
  #endif
        MHp->ntoggles = MH_FAILED;
        MHp->logratio = 0.0;
#else
      error("ERPM ErpmMix: invalid b1part state encountered while attempting SwapStep.");
#endif
    }
  } else {
    error("ERPM ErpmMix: internal error (unexpected move code).");
  }

#if DEBUG_ERPM_PROPOSALS
  if(DBG_seen_mix < DBG_max_print){
    Rprintf("[ERPM][Mix] move=%d | ntoggles=%d | logratio=%.6f\n",
            move, (int)MHp->ntoggles, MHp->logratio);
    DBG_seen_mix++;
  }
#endif
}

MH_F_FN(Mf_ErpmMix){
  erpm_mix_free_storage(MHp);
}

/* ========================================================================== */
/* MH_ErpmToggleStep                                                          */
/* ========================================================================== */
MH_P_FN(MH_ErpmToggleStep){
  erpm_require_undirected(nwp);

  if (MHp->ntoggles == 0) {
#if DEBUG_ERPM_PROPOSALS
    Rprintf("[ERPM][ToggleStep][INIT] setting ntoggles=2\n");
#endif
    MHp->ntoggles = 2;
    return;
  }

  ErpmToggleStep_propose(MHp, nwp);

#if DEBUG_ERPM_PROPOSALS
  if(MHp->ntoggles == 2 && DBG_seen_toggle < DBG_max_print){
    Rprintf("[ERPM][ToggleStep] actor=%d | group_old=%d -> group_new=%d | ntoggles=%d | "
            "toggles={(%d,%d),(%d,%d)} | logratio=%.6f\n",
            (int)Mtail[0], (int)Mhead[0], (int)Mhead[1], (int)MHp->ntoggles,
            (int)Mtail[0], (int)Mhead[0], (int)Mtail[1], (int)Mhead[1],
            MHp->logratio);
    DBG_seen_toggle++;
  }
#endif
}

/* ========================================================================== */
/* MH_ErpmSwapStep                                                            */
/* ========================================================================== */
MH_P_FN(MH_ErpmSwapStep){
  erpm_require_undirected(nwp);

  if (MHp->ntoggles == 0) {
#if DEBUG_ERPM_PROPOSALS
    Rprintf("[ERPM][SwapStep][INIT] setting ntoggles=4\n");
#endif
    MHp->ntoggles = 4;
    return;
  }

  /* Standalone behavior: strict => MH_FAILED (or error) if impossible. */
  ErpmSwapStep_propose_strict(MHp, nwp);

#if DEBUG_ERPM_PROPOSALS
  if(MHp->ntoggles == 4 && DBG_seen_swap < DBG_max_print){
    Rprintf("[ERPM][SwapStep] ntoggles=%d | toggles={(%d,%d),(%d,%d),(%d,%d),(%d,%d)} | logratio=%.6f\n",
            (int)MHp->ntoggles,
            (int)Mtail[0], (int)Mhead[0],
            (int)Mtail[1], (int)Mhead[1],
            (int)Mtail[2], (int)Mhead[2],
            (int)Mtail[3], (int)Mhead[3],
            MHp->logratio);
    DBG_seen_swap++;
  }
#endif
}
