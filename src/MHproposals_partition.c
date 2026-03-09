/**
 * @file MHproposals_partition.c
 * @brief Metropolis-Hastings proposals for ERPM partition networks under `~b1part`.
 * 
 * @author Jérémie Chichignoud
 *
 * @details
 * ERPM represents a partition as an undirected bipartite membership network:
 * - actor-mode vertices are `1..BIPARTITE`,
 * - group-mode vertices are `(BIPARTITE+1)..N_NODES`,
 * - under `~b1part`, each actor must have exactly one incident membership edge.
 *
 * This file implements the proposal side of that representation with the legacy
 * \pkg{ergm} MH API. The available moves are:
 * - `MH_ErpmToggleStep`: reassign one actor to another group (2 toggles),
 * - `MH_ErpmSwapStep`  : swap the groups of two actors (4 toggles),
 * - `MH_ErpmMergeStep` : move all actors of one non-empty group into another,
 * - `MH_ErpmSplitStep` : move a non-empty subset of one group into an empty group,
 * - `MH_ErpmMix`       : draw one move type from user weights, then apply it.
 *
 * A few conventions matter here:
 * - membership networks must be undirected; directed inputs are rejected,
 * - corrupted `~b1part` states are treated as proposal failures,
 * - standalone proposals are strict and fail when their move is infeasible,
 * - the mixed proposal keeps a state-independent first-stage draw and falls back
 *   to `TOGGLE` when the selected move is infeasible.
 *
 * The deliberate choice in `ErpmMix` is to interpret the user weights as
 * comparable attempt weights, not as exact observed move frequencies along the
 * chain. Feasibility is checked only after the move type has been drawn. This
 * keeps the first-stage mixture simple and avoids introducing a state-dependent
 * proposal distribution that would then need a matching Hastings correction in
 * `MHp->logratio`. The trade-off is straightforward: infeasible `SWAP`,
 * `MERGE`, or `SPLIT` draws are converted into `TOGGLE`s, so the observed number
 * of toggles is mechanically inflated whenever those moves are often infeasible.
 *
 * Legacy MH API reminders:
 * - proposals are written as `MH_P_FN(MH_<Name>)`,
 * - initialization is detected with `MHp->ntoggles == 0`,
 * - a valid proposal must fill `MHp->ntoggles`, `Mtail[...]`, `Mhead[...]`,
 *   and `MHp->logratio`.
 *
 * For the mixed kernel this has an important consequence: the standalone
 * proposal entry points all have an init branch, so `MH_ErpmMix` must never
 * dispatch to them directly. It uses helper routines with no init branch
 * instead, so every dispatch either emits a move or reports that the move is
 * infeasible in the current state.
 */

#include "ergm_MHproposals_degree.h"
#include "ergm_MHproposal.h" /* MH_FAILED */
#include "ergm_changestat.h"

#include <R_ext/Error.h>   /* error   */
#include <R_ext/Memory.h>  /* R_alloc, R_Calloc, R_Free */
#include <R_ext/Print.h>   /* Rprintf */
#include <Rmath.h>         /* unif_rand */
#include <math.h>          /* floor, log, isfinite */


/* ========================================================================= */
/* Debug helpers                                                             */
/* ========================================================================= */

/** @brief Enable verbose traces inside proposal code. */
#define DEBUG_ERPM_PROPOSALS 0

/** @brief Silence an intentionally unused variable. */
#define UNUSED_VARIABLE(x) (void)(x)

#if DEBUG_ERPM_PROPOSALS
/** @brief Maximum number of debug prints per move family. */
static const int DBG_max_print = 400;

/** @brief Debug print counters. */
static int DBG_seen_toggle = 0;
static int DBG_seen_swap   = 0;
static int DBG_seen_mix    = 0;

/**
 * @brief Print a compact description of the current partition network.
 *
 * @param who Short caller name for the log line.
 * @param nwp Current network.
 */
static inline void DBG_print_context_header(const char *who, Network *nwp){
  const Vertex n1 = BIPARTITE;
  const Vertex N  = N_NODES;
  const Vertex G  = N - n1;

  Rprintf("[ERPM][%s][CTX] n1=%d N=%d G_total=%d directed_flag=%d\n",
          who, (int)n1, (int)N, (int)G, (int)nwp->directed_flag);
}

/**
 * @brief Count non-empty groups for debug output.
 *
 * @param nwp Current network.
 * @return Number of non-empty group vertices.
 */
static inline Vertex DBG_count_nonempty_groups(Network *nwp){
  Vertex P = 0;
  for(Vertex g = nwp->bipartite + 1; g <= nwp->nnodes; ++g){
    if(nwp->indegree[g] > 0) P++;
  }
  return P;
}
#endif


/* ========================================================================= */
/* Small utilities                                                           */
/* ========================================================================= */
/**
 * @brief Abort if the membership network is directed.
 *
 * @param nwp Current network.
 */
static inline void erpm_require_undirected(Network *nwp){
  if(nwp->directed_flag){
    error("ERPM partition proposals require an undirected membership network (DIRECTED=TRUE detected).");
  }
}

/**
 * @brief Draw a uniform integer in `1..upper_bound`.
 *
 * @param upper_bound Inclusive upper bound.
 * @return A value in `1..upper_bound`.
 */
static inline Vertex UnifVertex1to(Vertex upper_bound){
  if(upper_bound < 1){
    error("ERPM proposals: UnifVertex1to() called with upper_bound < 1.");
  }
  return (Vertex)(1 + (Vertex)floor(unif_rand() * upper_bound));
}

/**
 * @brief Count the number of non-empty groups.
 *
 * @details
 * In this undirected bipartite membership representation, the indegree of a
 * group vertex is used as a fast proxy for its current size.
 *
 * @param nwp Current network.
 * @return Number of groups with size strictly greater than zero.
 */
static inline Vertex erpm_count_nonempty_groups(Network *nwp){
  Vertex P = 0;
  for(Vertex g = BIPARTITE + 1; g <= N_NODES; ++g){
    if(nwp->indegree[g] > 0) ++P;
  }
  return P;
}

/**
 * @brief Test whether a swap move is feasible.
 *
 * @param nwp Current network.
 * @return Non-zero iff at least two non-empty groups exist.
 */
static inline int erpm_can_swap(Network *nwp){
  return erpm_count_nonempty_groups(nwp) >= 2;
}

/**
 * @brief Test whether a merge move is feasible.
 *
 * @param nwp Current network.
 * @return Non-zero iff at least two non-empty groups exist.
 */
static inline int erpm_can_merge(Network *nwp){
  return erpm_count_nonempty_groups(nwp) >= 2;
}

/**
 * @brief Test whether a split move is feasible.
 *
 * @param nwp Current network.
 * @return Non-zero iff there is both a splittable group and an empty group.
 */
static inline int erpm_can_split(Network *nwp){
  int has_splittable_group = 0;
  int has_empty_group = 0;

  for(Vertex g = BIPARTITE + 1; g <= N_NODES; ++g){
    if(nwp->indegree[g] >= 2) has_splittable_group = 1;
    if(nwp->indegree[g] == 0) has_empty_group = 1;
    if(has_splittable_group && has_empty_group) return 1;
  }

  return 0;
}


/* ========================================================================= */
/* Membership decoding                                                       */
/* ========================================================================= */

/**
 * @brief Return the current group of one actor.
 *
 * @details
 * This is the slow-path decoder: it scans all group vertices and looks for the
 * unique incident membership edge of the actor. It is kept as the simple,
 * trustworthy baseline and is reused by the bulk decoder below.
 *
 * @param actor_id Actor vertex id in `1..BIPARTITE`.
 * @param nwp Current network.
 * @return Group vertex id, or `0` if no valid group could be found.
 */
static Vertex get_current_group_of_actor(Vertex actor_id, Network *nwp){
  erpm_require_undirected(nwp);

  if(actor_id < 1 || actor_id > BIPARTITE) return 0;

  Vertex found_group = 0;

  for(Vertex group_id = BIPARTITE + 1; group_id <= N_NODES; ++group_id){
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
    Rprintf("[ERPM][WARN] actor=%d has no group neighbor (invalid b1part state?)\n",
            (int)actor_id);
  }
#endif

  return found_group;
}

/**
 * @brief Decode the full actor-to-group map for the current state.
 *
 * @param ag Output array indexed by actor id. Index `0` is unused.
 * @param nwp Current network.
 * @return `1` on success, `0` if the state violates `~b1part`.
 */
static int erpm_decode_actor_groups(Vertex *ag, Network *nwp){
  for(Vertex a = 1; a <= BIPARTITE; ++a){
    Vertex g = get_current_group_of_actor(a, nwp);
    if(g == 0) return 0;
    ag[a] = g;
  }
  return 1;
}

/**
 * @brief Collect all actors currently assigned to a given group.
 *
 * @param group_id Group vertex id.
 * @param ag Actor-to-group map.
 * @param buf Output buffer of length at least `BIPARTITE`.
 * @param nwp Current network. Unused, kept for call-site symmetry.
 * @return Number of actors collected into `buf`.
 */
static Vertex erpm_collect_group_actors(Vertex group_id, Vertex *ag, Vertex *buf, Network *nwp){
  UNUSED_VARIABLE(nwp);

  Vertex k = 0;
  for(Vertex a = 1; a <= BIPARTITE; ++a){
    if(ag[a] == group_id){
      buf[k++] = a;
    }
  }
  return k;
}

/**
 * @brief Sample one empty group uniformly.
 *
 * @param exclude_g1 First excluded group id.
 * @param exclude_g2 Second excluded group id.
 * @param nwp Current network.
 * @return Empty group vertex id, or `0` if none exists.
 */
static Vertex erpm_sample_empty_group(Vertex exclude_g1, Vertex exclude_g2, Network *nwp){
  Vertex empty_count = 0;

  for(Vertex g = BIPARTITE + 1; g <= N_NODES; ++g){
    if(g == exclude_g1 || g == exclude_g2) continue;
    if(nwp->indegree[g] == 0) empty_count++;
  }
  if(empty_count == 0) return 0;

  Vertex r = UnifVertex1to(empty_count);
  for(Vertex g = BIPARTITE + 1; g <= N_NODES; ++g){
    if(g == exclude_g1 || g == exclude_g2) continue;
    if(nwp->indegree[g] == 0){
      r--;
      if(r == 0) return g;
    }
  }

  return 0;
}


/* ========================================================================= */
/* Reference proposal from ergm                                              */
/* ========================================================================= */
/**
 * @brief Reference copy of `MH_B1Part` from \pkg{ergm}.
 *
 * @details
 * This is kept here as a local reference for the state-dependent combinatorial
 * correction used by `B1Part`. It is not part of the ERPM mixed proposal logic,
 * but it remains useful context when comparing proposal designs under `~b1part`.
 */
MH_P_FN(MH_B1Part) {
  if (MHp->ntoggles == 0) {
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

    MHp->logratio += dP == -1 ? log(N_NODES - BIPARTITE - P + 1)
                              : -log(N_NODES - BIPARTITE - P);
  }
}


/* ========================================================================= */
/* Per-iteration move helpers                                                */
/* ========================================================================= */
/**
 * @brief Emit one toggle proposal.
 *
 * @details
 * One actor is drawn uniformly, its current group is decoded, and a distinct
 * destination group is drawn uniformly among all other group vertices, including
 * padded empty groups. This keeps the move simple and symmetric.
 *
 * @param MHp Proposal object to fill.
 * @param nwp Current network.
 */
static void ErpmToggleStep_propose(MHProposal *MHp, Network *nwp){
  erpm_require_undirected(nwp);

  const Vertex n_actors = BIPARTITE;
  const Vertex n_groups = N_NODES - BIPARTITE;

  Vertex actor_id  = UnifVertex1to(n_actors);
  Vertex group_old = get_current_group_of_actor(actor_id, nwp);

  if(group_old == 0){
#ifdef MH_FAILED
  #if DEBUG_ERPM_PROPOSALS
    if(DBG_seen_toggle < DBG_max_print){
      DBG_print_context_header("ToggleStep", nwp);
      Rprintf("[ERPM][ToggleStep][FAIL] MH_FAILED (actor has no group)\n");
      Rprintf("[ERPM][ToggleStep][FAIL] actor=%d | nonempty_groups(P)=%d\n",
              (int)actor_id, (int)DBG_count_nonempty_groups(nwp));
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

  if(n_groups < 2){
#ifdef MH_FAILED
    MHp->ntoggles = MH_FAILED;
    MHp->logratio = 0.0;
#else
    error("ERPM ToggleStep: need at least 2 group vertices to propose a move.");
#endif
    return;
  }

  const Vertex old_rank = group_old - BIPARTITE;
  const Vertex r = UnifVertex1to(n_groups - 1);
  const Vertex new_rank = (r >= old_rank) ? (r + 1) : r;
  const Vertex group_new = BIPARTITE + new_rank;

  MHp->ntoggles = 2;
  Mtail[0] = actor_id; Mhead[0] = group_old;
  Mtail[1] = actor_id; Mhead[1] = group_new;
  MHp->logratio = 0.0;
}

/**
 * @brief Try to emit one swap proposal.
 *
 * @details
 * This helper has no init branch. It is used both by the strict standalone
 * wrapper and by the mixed proposal. The move preserves group sizes exactly.
 *
 * @param MHp Proposal object to fill.
 * @param nwp Current network.
 * @return `1` on success, `0` if no legal swap exists, `-1` on corrupted state.
 */
static int ErpmSwapStep_trypropose(MHProposal *MHp, Network *nwp){
  erpm_require_undirected(nwp);

  if(!erpm_can_swap(nwp)){
    return 0;
  }

  Vertex *ag = (Vertex*) R_alloc((size_t)(BIPARTITE + 1), sizeof(Vertex));
  if(!erpm_decode_actor_groups(ag, nwp)){
    return -1;
  }

  Vertex actor_i = UnifVertex1to(BIPARTITE);
  Vertex group_i = ag[actor_i];
  if(group_i == 0) return -1;

  Vertex actor_j = 0;
  Vertex group_j = 0;

  for(;;){
    Vertex cand_actor = UnifVertex1to(BIPARTITE);
    if(cand_actor == actor_i) continue;

    Vertex cand_group = ag[cand_actor];
    if(cand_group == 0) return -1;

    if(cand_group != group_i){
      actor_j = cand_actor;
      group_j = cand_group;
      break;
    }
  }

  MHp->ntoggles = 4;
  Mtail[0] = actor_i; Mhead[0] = group_i;
  Mtail[1] = actor_i; Mhead[1] = group_j;
  Mtail[2] = actor_j; Mhead[2] = group_j;
  Mtail[3] = actor_j; Mhead[3] = group_i;
  MHp->logratio = 0.0;

  return 1;
}

/**
 * @brief Apply standalone swap semantics.
 *
 * @details
 * The standalone proposal stays strict: if a swap is infeasible, it reports
 * failure instead of silently changing move type.
 *
 * @param MHp Proposal object to fill.
 * @param nwp Current network.
 */
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

/**
 * @brief Try to emit one merge proposal.
 *
 * @details
 * Two distinct non-empty groups are drawn uniformly among non-empty groups.
 * Every actor currently in the source group is reassigned to the destination.
 *
 * @param MHp Proposal object to fill.
 * @param nwp Current network.
 * @return `1` on success, `0` if no legal merge exists, `-1` on corrupted state.
 */
static int ErpmMergeStep_trypropose(MHProposal *MHp, Network *nwp){
  erpm_require_undirected(nwp);

  if(!erpm_can_merge(nwp)){
    return 0;
  }

  Vertex *ag = (Vertex*) R_alloc((size_t)(BIPARTITE + 1), sizeof(Vertex));
  if(!erpm_decode_actor_groups(ag, nwp)){
    return -1;
  }

  Vertex P  = erpm_count_nonempty_groups(nwp);
  Vertex r1 = UnifVertex1to(P);
  Vertex r2 = UnifVertex1to(P - 1);

  Vertex g1 = 0;
  Vertex g2 = 0;

  for(Vertex g = BIPARTITE + 1; g <= N_NODES; ++g){
    if(nwp->indegree[g] > 0){
      if(--r1 == 0){ g1 = g; break; }
    }
  }
  for(Vertex g = BIPARTITE + 1; g <= N_NODES; ++g){
    if(nwp->indegree[g] > 0 && g != g1){
      if(--r2 == 0){ g2 = g; break; }
    }
  }

  if(g1 == 0 || g2 == 0){
    return 0;
  }

  Vertex *buf = (Vertex*) R_alloc((size_t)BIPARTITE, sizeof(Vertex));
  Vertex k = erpm_collect_group_actors(g2, ag, buf, nwp);
  if(k == 0){
    return 0;
  }

  MHp->ntoggles = (int)(2 * k);

  Vertex t = 0;
  for(Vertex j = 0; j < k; ++j){
    Vertex a = buf[j];
    Mtail[t] = a; Mhead[t] = g2; ++t; /* OFF */
    Mtail[t] = a; Mhead[t] = g1; ++t; /* ON  */
  }

  MHp->logratio = 0.0;
  return 1;
}

/**
 * @brief Apply standalone merge semantics.
 *
 * @param MHp Proposal object to fill.
 * @param nwp Current network.
 */
static void ErpmMergeStep_propose_strict(MHProposal *MHp, Network *nwp){
  const int rc = ErpmMergeStep_trypropose(MHp, nwp);

  if(rc == 1) return;

  if(rc == 0){
#ifdef MH_FAILED
    MHp->ntoggles = MH_FAILED;
    MHp->logratio = 0.0;
#else
    error("ERPM MergeStep: no legal merge exists (need at least two non-empty groups).");
#endif
    return;
  }

#ifdef MH_FAILED
  MHp->ntoggles = MH_FAILED;
  MHp->logratio = 0.0;
#else
  error("ERPM MergeStep: invalid b1part state (actor has no group neighbor).");
#endif
}

/**
 * @brief Try to emit one split proposal.
 *
 * @details
 * One splittable group and one empty destination group are drawn, then a
 * non-trivial subset of actors is moved to the empty group.
 *
 * @param MHp Proposal object to fill.
 * @param nwp Current network.
 * @return `1` on success, `0` if no legal split exists, `-1` on corrupted state.
 */
static int ErpmSplitStep_trypropose(MHProposal *MHp, Network *nwp){
  erpm_require_undirected(nwp);

  if(!erpm_can_split(nwp)){
    return 0;
  }

  Vertex *ag = (Vertex*) R_alloc((size_t)(BIPARTITE + 1), sizeof(Vertex));
  if(!erpm_decode_actor_groups(ag, nwp)){
    return -1;
  }

  Vertex *splittable = (Vertex*) R_alloc((size_t)(N_NODES + 1), sizeof(Vertex));
  Vertex S = 0;

  for(Vertex g = BIPARTITE + 1; g <= N_NODES; ++g){
    if(nwp->indegree[g] >= 2){
      splittable[S++] = g;
    }
  }
  if(S == 0){
    return 0;
  }

  Vertex g  = splittable[UnifVertex1to(S) - 1];
  Vertex gp = erpm_sample_empty_group(g, 0, nwp);
  if(gp == 0){
    return 0;
  }

  Vertex *buf = (Vertex*) R_alloc((size_t)BIPARTITE, sizeof(Vertex));
  Vertex k = erpm_collect_group_actors(g, ag, buf, nwp);
  if(k < 2){
    return 0;
  }

  Vertex m = UnifVertex1to(k - 1);

  for(Vertex j = 0; j < m; ++j){
    Vertex u   = j + (Vertex)floor(unif_rand() * (k - j));
    Vertex tmp = buf[j];
    buf[j] = buf[u];
    buf[u] = tmp;
  }

  MHp->ntoggles = (int)(2 * m);

  Vertex t = 0;
  for(Vertex j = 0; j < m; ++j){
    Vertex a = buf[j];
    Mtail[t] = a; Mhead[t] = g;  ++t; /* OFF */
    Mtail[t] = a; Mhead[t] = gp; ++t; /* ON  */
  }

  MHp->logratio = 0.0;
  return 1;
}

/**
 * @brief Apply standalone split semantics.
 *
 * @param MHp Proposal object to fill.
 * @param nwp Current network.
 */
static void ErpmSplitStep_propose_strict(MHProposal *MHp, Network *nwp){
  const int rc = ErpmSplitStep_trypropose(MHp, nwp);

  if(rc == 1) return;

  if(rc == 0){
#ifdef MH_FAILED
    MHp->ntoggles = MH_FAILED;
    MHp->logratio = 0.0;
#else
    error("ERPM SplitStep: no legal split exists (need a group of size >= 2 and one empty group).");
#endif
    return;
  }

#ifdef MH_FAILED
  MHp->ntoggles = MH_FAILED;
  MHp->logratio = 0.0;
#else
  error("ERPM SplitStep: invalid b1part state (actor has no group neighbor).");
#endif
}


/* ========================================================================= */
/* Mixed proposal storage                                                    */
/* ========================================================================= */

/** @brief Stable move code for `TOGGLE`. */
#define ERPM_MOVE_TOGGLE 1
/** @brief Stable move code for `SWAP`. */
#define ERPM_MOVE_SWAP   2
/** @brief Stable move code for `MERGE`. */
#define ERPM_MOVE_MERGE  3
/** @brief Stable move code for `SPLIT`. */
#define ERPM_MOVE_SPLIT  4

/**
 * @brief Persistent storage for `ErpmMix`.
 *
 * @details
 * The R initializer packs the user-requested move codes and weights. Here they
 * are normalized into cumulative probabilities for fast first-stage sampling.
 */
typedef struct ErpmMixStorage {
  int K;              /**< Number of move types in the mixture. */
  int    *move_codes; /**< Move code array of length `K`. */
  double *cumprob;    /**< Cumulative probabilities of length `K`. */
} ErpmMixStorage;

/**
 * @brief Install the canonical default mix.
 *
 * @details
 * The canonical fallback is the historical `toggle:2, swap:1` mixture. It is
 * used whenever the R-side packing is missing or invalid.
 *
 * @param st Storage object to populate.
 */
static void erpm_mix_set_default(ErpmMixStorage *st){
  if(!st) error("ERPM ErpmMix: internal error (NULL storage).");

  st->K = 2;
  st->move_codes = (int*)    R_Calloc((size_t)st->K, int);
  st->cumprob    = (double*) R_Calloc((size_t)st->K, double);

  st->move_codes[0] = ERPM_MOVE_TOGGLE;
  st->move_codes[1] = ERPM_MOVE_SWAP;

  st->cumprob[0] = 2.0 / 3.0;
  st->cumprob[1] = 1.0;
}

/**
 * @brief Build mixed-proposal storage from R-packed inputs.
 *
 * @param st Storage object to populate.
 * @param MHp Proposal object containing packed `iinputs` and `inputs`.
 * @return `1` on success, `0` if the packing is missing or invalid.
 */
static int erpm_mix_try_build_from_inputs(ErpmMixStorage *st, MHProposal *MHp){
  if(!st)  error("ERPM ErpmMix: internal error (NULL storage).");
  if(!MHp) error("ERPM ErpmMix: internal error (NULL MHp).");

  if(MHp->iinputs == NULL) return 0;
  if(MHp->inputs  == NULL) return 0;

  const int K = (int)MHp->iinputs[0];
  if(K <= 0 || K > 32) return 0;

  st->K = K;
  st->move_codes = (int*)    R_Calloc((size_t)K, int);
  st->cumprob    = (double*) R_Calloc((size_t)K, double);

  double wsum = 0.0;

  for(int k = 0; k < K; ++k){
    const int code  = (int)MHp->iinputs[1 + k];
    const double w  = MHp->inputs[k];

    if(code != ERPM_MOVE_TOGGLE &&
       code != ERPM_MOVE_SWAP   &&
       code != ERPM_MOVE_MERGE  &&
       code != ERPM_MOVE_SPLIT) return 0;

    if(!(isfinite(w) && w > 0.0)) return 0;

    st->move_codes[k] = code;
    wsum += w;
  }

  if(!(isfinite(wsum) && wsum > 0.0)) return 0;

  double c = 0.0;
  for(int k = 0; k < K; ++k){
    c += MHp->inputs[k] / wsum;
    st->cumprob[k] = c;
  }
  st->cumprob[K - 1] = 1.0;

  return 1;
}

/**
 * @brief Sample one move code from the mixed proposal.
 *
 * @param st Initialized mixed-proposal storage.
 * @return One of the `ERPM_MOVE_*` constants.
 */
static int erpm_mix_sample_move(const ErpmMixStorage *st){
  const double u = unif_rand();
  for(int k = 0; k < st->K; ++k){
    if(u < st->cumprob[k]) return st->move_codes[k];
  }
  return st->move_codes[st->K - 1];
}

/**
 * @brief Free `ErpmMix` storage.
 *
 * @param MHp Proposal object whose `storage` field owns an `ErpmMixStorage`.
 */
static void erpm_mix_free_storage(MHProposal *MHp){
  ErpmMixStorage *st = (ErpmMixStorage*) MHp->storage;
  if(!st) return;

  if(st->move_codes) R_Free(st->move_codes);
  if(st->cumprob)    R_Free(st->cumprob);
  R_Free(st);

  MHp->storage = NULL;
}


/* ========================================================================= */
/* Mixed proposal entry points                                               */
/* ========================================================================= */

/**
 * @brief Initialize storage for `MH_ErpmMix`.
 *
 * @details
 * The proposal arrays must be large enough for the biggest supported move, so
 * the init branch reserves up to `2 * BIPARTITE` toggles, which is enough for
 * a full merge.
 */
MH_I_FN(Mi_ErpmMix){
  if(MHp->storage != NULL){
    return;
  }

  ErpmMixStorage *st = (ErpmMixStorage*) R_Calloc(1, ErpmMixStorage);

  if(!erpm_mix_try_build_from_inputs(st, MHp)){
#if DEBUG_ERPM_PROPOSALS
    Rprintf("[ERPM][Mix][INIT] invalid packing -> fallback default (toggle:2 swap:1)\n");
#endif
    if(st->move_codes) R_Free(st->move_codes);
    if(st->cumprob)    R_Free(st->cumprob);
    st->move_codes = NULL;
    st->cumprob    = NULL;
    erpm_mix_set_default(st);
  }

  MHp->storage  = st;
  MHp->ntoggles = (int)(2 * BIPARTITE);
  MHp->logratio = 0.0;
}

/**
 * @brief Emit one proposal from the mixed kernel.
 *
 * @details
 * The move type is first drawn from the fixed user-level weights. Feasibility is
 * checked only afterwards. If the selected move is infeasible, the mixed kernel
 * falls back to `TOGGLE` instead of failing. This is the core policy choice of
 * the current implementation.
 */
MH_P_FN(MH_ErpmMix){
  erpm_require_undirected(nwp);

  ErpmMixStorage *st = (ErpmMixStorage*) MHp->storage;
  if(st == NULL){
    error("ERPM ErpmMix: NULL storage (Mi_ErpmMix was not called).");
  }

  int move = erpm_mix_sample_move(st);

  if(move == ERPM_MOVE_TOGGLE){
    ErpmToggleStep_propose(MHp, nwp);

  } else if(move == ERPM_MOVE_SWAP){
    const int rc = ErpmSwapStep_trypropose(MHp, nwp);

    if(rc == 1){
      /* move emitted */
    } else if(rc == 0){
#if DEBUG_ERPM_PROPOSALS
      if(DBG_seen_mix < DBG_max_print){
        Rprintf("[ERPM][Mix][FALLBACK] SWAP impossible -> TOGGLE\n");
        DBG_seen_mix++;
      }
#endif
      ErpmToggleStep_propose(MHp, nwp);
    } else {
#ifdef MH_FAILED
      MHp->ntoggles = MH_FAILED;
      MHp->logratio = 0.0;
#else
      error("ERPM ErpmMix: invalid b1part state encountered while attempting SwapStep.");
#endif
    }

  } else if(move == ERPM_MOVE_MERGE){
    const int rc = ErpmMergeStep_trypropose(MHp, nwp);

    if(rc == 1){
      /* move emitted */
    } else if(rc == 0){
#if DEBUG_ERPM_PROPOSALS
      if(DBG_seen_mix < DBG_max_print){
        Rprintf("[ERPM][Mix][FALLBACK] MERGE impossible -> TOGGLE\n");
        DBG_seen_mix++;
      }
#endif
      ErpmToggleStep_propose(MHp, nwp);
    } else {
#ifdef MH_FAILED
      MHp->ntoggles = MH_FAILED;
      MHp->logratio = 0.0;
#else
      error("ERPM ErpmMix: invalid b1part state encountered while attempting MergeStep.");
#endif
    }

  } else if(move == ERPM_MOVE_SPLIT){
    const int rc = ErpmSplitStep_trypropose(MHp, nwp);

    if(rc == 1){
      /* move emitted */
    } else if(rc == 0){
#if DEBUG_ERPM_PROPOSALS
      if(DBG_seen_mix < DBG_max_print){
        Rprintf("[ERPM][Mix][FALLBACK] SPLIT impossible -> TOGGLE\n");
        DBG_seen_mix++;
      }
#endif
      ErpmToggleStep_propose(MHp, nwp);
    } else {
#ifdef MH_FAILED
      MHp->ntoggles = MH_FAILED;
      MHp->logratio = 0.0;
#else
      error("ERPM ErpmMix: invalid b1part state encountered while attempting SplitStep.");
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

/**
 * @brief Free `MH_ErpmMix` storage.
 *
 * @param MHp Proposal object.
 */
MH_F_FN(Mf_ErpmMix){
  erpm_mix_free_storage(MHp);
}


/* ========================================================================= */
/* Standalone proposal entry points                                          */
/* ========================================================================= */

/**
 * @brief Standalone toggle proposal entry point.
 *
 * @details
 * Init reserves exactly 2 toggles. Subsequent calls emit one toggle move.
 */
MH_P_FN(MH_ErpmToggleStep){
  erpm_require_undirected(nwp);

  if(MHp->ntoggles == 0){
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

/**
 * @brief Standalone swap proposal entry point.
 *
 * @details
 * Init reserves 4 toggles. Subsequent calls apply the strict swap semantics.
 */
MH_P_FN(MH_ErpmSwapStep){
  erpm_require_undirected(nwp);

  if(MHp->ntoggles == 0){
#if DEBUG_ERPM_PROPOSALS
    Rprintf("[ERPM][SwapStep][INIT] setting ntoggles=4\n");
#endif
    MHp->ntoggles = 4;
    return;
  }

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

/**
 * @brief Standalone merge proposal entry point.
 *
 * @details
 * Init reserves enough room for the largest possible merge, namely moving all
 * actors of one group into another.
 */
MH_P_FN(MH_ErpmMergeStep){
  erpm_require_undirected(nwp);

  if(MHp->ntoggles == 0){
    MHp->ntoggles = (int)(2 * BIPARTITE);
    return;
  }

  ErpmMergeStep_propose_strict(MHp, nwp);
}

/**
 * @brief Standalone split proposal entry point.
 *
 * @details
 * Init reserves enough room for the largest possible split, namely moving
 * `BIPARTITE - 1` actors out of one group.
 */
MH_P_FN(MH_ErpmSplitStep){
  erpm_require_undirected(nwp);

  if(MHp->ntoggles == 0){
    MHp->ntoggles = (int)(2 * (BIPARTITE - 1));
    return;
  }

  ErpmSplitStep_propose_strict(MHp, nwp);
}
