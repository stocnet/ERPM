/**
 * @file MHproposals_partition.c
 * @brief Metropolis-Hastings proposals for ERPM partition networks under '~b1part'.
 *
 * ERPM represents a partition as an undirected bipartite "membership" network:
 * - Actor-mode vertices are  1..BIPARTITE (n actors).
 * - Group-mode vertices are (BIPARTITE+1)..N_NODES (G_total groups, padded).
 * - Under '~b1part', each actor has exactly one incident membership edge.
 *
 * Convention / safety:
 * - Membership networks MUST be undirected. Directed networks are rejected (fatal error).
 * - Corrupted '~b1part' states (actor without group) trigger MH_FAILED (if available)
 *   or a hard error.
 *
 * Implemented proposals (legacy ergm MH API):
 * - MH_ErpmToggleStep: move one actor to another group (2 toggles).
 * - MH_ErpmSwapStep  : swap memberships of two actors in different groups (4 toggles).
 * - MH_ErpmMix       : mixed proposal selecting ToggleStep/SwapStep by user weights.
 *
 * Legacy MH API notes (ergm_MHproposal.h):
 * - Proposals are MH_P_FN(MH_<Name>) and do not return a status.
 * - Initialization is detected when MHp->ntoggles == 0.
 * - A valid proposal must fill:
 *     MHp->ntoggles, Mtail[...], Mhead[...], MHp->logratio.
 *
 * IMPORTANT (mixed proposal wiring):
 * - MH_ErpmToggleStep and MH_ErpmSwapStep contain an init branch:
 *     if(MHp->ntoggles==0){ MHp->ntoggles=fixed; return; }
 *   Therefore MH_ErpmMix MUST NOT call them directly, or the first dispatch would
 *   produce a no-move.
 * - MH_ErpmMix dispatches to per-iteration helpers:
 *     ErpmToggleStep_propose() / ErpmSwapStep_trypropose()
 *   which have no init branch and always emit a proposal (or report impossibility).
 *
 * ---------------------------------------------------------------------------
 * Proposal semantics (summary)
 * ---------------------------------------------------------------------------
 * 1) ToggleStep (2 toggles)
 * - Choose an actor uniformly among 1..BIPARTITE.
 * - Read its current group_old (unique neighbor in group mode).
 * - Choose group_new uniformly among ALL group vertices except group_old,
 *   including padded empty groups.
 * - Encode:
 *     (actor, group_old) OFF
 *     (actor, group_new) ON
 * - Symmetric -> logratio = 0.
 *
 * 2) SwapStep (4 toggles)
 * - Requires at least two non-empty groups.
 * - Choose actor_i uniformly; read group_i.
 * - Choose actor_j uniformly until actor_j != actor_i AND group_j != group_i.
 * - Encode 4 toggles swapping memberships.
 * - Symmetric -> logratio = 0.
 * - Preserves group sizes exactly (intra-size mixing move).
 *
 * 3) ErpmMix (Toggle/Swap mixture)  [UPDATED]
 * - Randomly selects TOGGLE or SWAP according to user weights.
 * - If SWAP is selected but impossible (P<2 non-empty groups), ErpmMix falls back
 *   to TOGGLE (no MH_FAILED for this specific reason).
 * - ErpmMix fully relies on InitErgmProposal.ErpmMix packing:
 *     iinputs = c(K, move_codes...)
 *     inputs  = c(weights...)
 * - ErpmMix must be robust to missing/invalid iinputs/inputs and MUST fall back
 *   to the canonical mix (toggle:2, swap:1) instead of crashing.
 * - Only TOGGLE and SWAP are supported for now.
 *
 * ---------------------------------------------------------------------------
 * Interaction with ergm internals
 * ---------------------------------------------------------------------------
 * - Uses ergm Network macros/structures:
 *     - BIPARTITE, N_NODES, IN_DEG[...] (SwapStep uses nwp->indegree; IN_DEG kept for MH_B1Part ref)
 *     - IS_UNDIRECTED_EDGE(u,v) to test membership edges.
 * - Undirected requirement checked through nwp->directed_flag.
 *
 * Failure behavior:
 * - directed network => error()
 * - invalid '~b1part' state => MH_FAILED (if available) or error()
 * - no legal SwapStep:
 *     standalone SwapStep => MH_FAILED/error()
 *     inside Mix (when SWAP selected) => fallback ToggleStep
 *
 * ---------------------------------------------------------------------------
 * Debugging
 * ---------------------------------------------------------------------------
 * - DEBUG_ERPM_PROPOSALS enables verbose traces.
 * - DBG_max_print limits console output during long MCMC runs.
 * - Debugging is intended for development/selftests and should remain disabled
 *   by default in production builds.
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
#define DEBUG_ERPM_PROPOSALS 1
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
/* We keep a MH_FAILED guard to be robust to ergm changes (to be remove later) */
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

  /* Draw group_new uniformly from ALL groups except group_old, without a retry loop. */
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

  /* rc == -1 : corrupted state
   * We keep a MH_FAILED guard to be robust to ergm changes (to be remove later)
   */
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
  int    *move_codes;  /* length K */
  double *cumprob;     /* length K, last = 1.0 */
} ErpmMixStorage;

/* Canonical default mix used if user packing is missing/invalid. */
static void erpm_mix_set_default(ErpmMixStorage *st){
  if(!st) error("ERPM ErpmMix: internal error (NULL storage).");

  st->K = 2;
  st->move_codes = (int*)    R_Calloc((size_t)st->K, int);
  st->cumprob    = (double*) R_Calloc((size_t)st->K, double);

  st->move_codes[0] = ERPM_MOVE_TOGGLE;  /* weight 2 */
  st->move_codes[1] = ERPM_MOVE_SWAP;    /* weight 1 */

  st->cumprob[0] = 2.0/3.0;
  st->cumprob[1] = 1.0;
}

/* Validate user-provided packing and build cumulative probabilities.
 * Returns 1 on success, 0 if invalid (caller should use defaults).
 *
 * IMPORTANT (allocation discipline):
 * - This function allocates st->move_codes and st->cumprob on entry.
 * - If returning 0, caller MUST free these partial allocations (or reuse helper).
 */
static int erpm_mix_try_build_from_inputs(ErpmMixStorage *st, MHProposal *MHp){
  if(!st)  error("ERPM ErpmMix: internal error (NULL storage).");
  if(!MHp) error("ERPM ErpmMix: internal error (NULL MHp).");

  if(MHp->iinputs == NULL) return 0;
  if(MHp->inputs  == NULL) return 0;

  const int K = (int)MHp->iinputs[0];
  if(K <= 0) return 0;

  /* Defensive: K must not be absurd. */
  if(K > 32) return 0;

  /* Allocate. */
  st->K = K;
  st->move_codes = (int*)    R_Calloc((size_t)K, int);
  st->cumprob    = (double*) R_Calloc((size_t)K, double);

  /* Validate codes + weights and compute normalization. */
  double wsum = 0.0;

  for(int k = 0; k < K; ++k){
    const int code = (int)MHp->iinputs[1 + k];
    const double w = MHp->inputs[k];

    if(code != ERPM_MOVE_TOGGLE && code != ERPM_MOVE_SWAP) return 0;
    if(!(isfinite(w) && w > 0.0)) return 0;

    st->move_codes[k] = code;
    wsum += w;
  }
  if(!(isfinite(wsum) && wsum > 0.0)) return 0;

  /* Build cumulative probabilities. */
  double c = 0.0;
  for(int k = 0; k < K; ++k){
    c += MHp->inputs[k] / wsum;
    st->cumprob[k] = c;
  }
  st->cumprob[K - 1] = 1.0;

  return 1;
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
  /* If ergm re-calls Mi_* defensively, do not rebuild/allocate/print again. */
  if(MHp->storage != NULL){
    return;
  }

#if DEBUG_ERPM_PROPOSALS
  Rprintf("[ERPM][Mix][INIT] building storage from iinputs/inputs (or defaults)\n");
#endif

  ErpmMixStorage *st = (ErpmMixStorage*) R_Calloc(1, ErpmMixStorage);

  /* Try user config; if invalid, fallback to canonical default mix. */
  if(!erpm_mix_try_build_from_inputs(st, MHp)){
#if DEBUG_ERPM_PROPOSALS
    Rprintf("[ERPM][Mix][INIT] invalid packing -> fallback default (toggle:2 swap:1)\n");
#endif
    /* If try_build allocated partial memory then returned 0, we must free it. */
    if(st->move_codes) R_Free(st->move_codes);
    if(st->cumprob)    R_Free(st->cumprob);
    st->move_codes = NULL;
    st->cumprob    = NULL;

    erpm_mix_set_default(st);
  }

  MHp->storage = st;

  /* Max toggles among component moves: SWAP is 4. */
  MHp->ntoggles = 4;
  MHp->logratio = 0.0;
}

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
      /* OK */

    } else if(rc == 0){
      /* SWAP impossible -> fallback TOGGLE (never MH_FAILED for this case). */
#if DEBUG_ERPM_PROPOSALS
      if(DBG_seen_mix < DBG_max_print){
        Rprintf("[ERPM][Mix][FALLBACK] SWAP impossible -> TOGGLE\n");
        DBG_seen_mix++;
      }
#endif
      ErpmToggleStep_propose(MHp, nwp);

    } else {
      /* rc == -1: corrupted state */
#ifdef MH_FAILED
      MHp->ntoggles = MH_FAILED;
      MHp->logratio = 0.0;
#else
      error("ERPM ErpmMix: invalid b1part state encountered while attempting SwapStep.");
#endif
    }

  } else {
    /* Should not happen (storage validation), but keep hard error. */
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
