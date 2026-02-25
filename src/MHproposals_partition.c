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
 * - MH_ErpmSwapStep  : swap memberships of two actors in different non-empty groups (4 toggles).
 *
 * Legacy MH API notes (ergm_MHproposal.h):
 * - Proposals are implemented as MH_P_FN(MH_<Name>) and do not return a status.
 * - Initialization is detected when MHp->ntoggles == 0, and the proposal must set
 *   its fixed toggle count (2 or 4 here).
 * - A valid proposal must fill:
 *   MHp->ntoggles, Mtail[...], Mhead[...], and MHp->logratio.
 */

#include "ergm_MHproposals_degree.h"
#include "ergm_MHproposal.h"
#include "ergm_changestat.h"

#include <R_ext/Print.h>   /* Rprintf */
#include <R_ext/Error.h>   /* error()  */
#include <math.h>          /* floor, log */


/* -------------------------------------------------------------------------- */
/* Debugging                                                                  */
/* -------------------------------------------------------------------------- */
/**
 * @def DEBUG_ERPM_PROPOSALS
 * @brief Enable verbose debugging output.
 *
 * Set this macro to 1 to print detailed information to the R console during
 * summary() or MCMC runs.
 *
 * When set to 0, the compiled code does not emit any debug traces.
 */
#define DEBUG_ERPM_PROPOSALS 0
#define UNUSED_VARIABLE(x) (void)(x)

/* Print limiter (avoid flooding console during long MCMC). */
#if DEBUG_ERPM_PROPOSALS
  static int DBG_seen_toggle = 0;
  static int DBG_seen_swap   = 0;
  static const int DBG_max_print = 400;
#endif

/* -------------------------------------------------------------------------- */
/* Guards                                                                     */
/* -------------------------------------------------------------------------- */
/**
 * @brief Enforce the ERPM requirement that the membership network is undirected.
 * @param nwp Pointer to the current ergm Network.
 *
 * The ERPM partition representation is an undirected membership graph. Running
 * these proposals on a directed network indicates a programming/configuration
 * error and is treated as fatal.
 */
static inline void erpm_require_undirected(Network *nwp){
  /* Note: DIRECTED macro from ergm_changestat_common.do_not_include_directly.h
   * is not reachable from this compilation unit, so we use nwp->directed_flag.
   */
  if(nwp->directed_flag){
    error("ERPM partition proposals require an undirected membership network (DIRECTED=TRUE detected).");
  }
}

/* -------------------------------------------------------------------------- */
/* RNG helper                                                                 */
/* -------------------------------------------------------------------------- */
/**
 * @brief Draw a uniform integer in {1, ..., upper_bound}.
 * @param upper_bound Upper bound (must be >= 1).
 * @return Uniform draw in 1..upper_bound.
 */
static inline Vertex UnifVertex1to(Vertex upper_bound){
  return (Vertex)(1 + (Vertex)floor(unif_rand() * upper_bound));
}

/* -------------------------------------------------------------------------- */
/* Membership decoding                                                        */
/* -------------------------------------------------------------------------- */
/**
 * @brief Get the current group (mode-2 vertex id) of a given actor under '~b1part'.
 * @param actor_id Actor vertex id in 1..BIPARTITE.
 * @param nwp Pointer to the current ergm Network.
 * @return Group vertex id in (BIPARTITE+1)..N_NODES, or 0 if not found.
 *
 * Expected invariant under '~b1part': actor_id has exactly one group neighbor.
 *
 * Implementation notes:
 * - Robust approach with minimal assumptions: scan all group vertices and
 *   test whether the membership edge exists.
 * - This is O(number_of_groups_total) per call.
 * - Under ERPM padding, the total number of group vertices is typically close
 *   to the number of actors, so this remains acceptable for the two proposals
 *   implemented here.
 * - On an invalid state (no neighbor) we return 0; callers treat it as a no-move.
 * - If multiple neighbors are detected, we keep the first one and warn in debug
 *   mode; this indicates a corrupted '~b1part' state.
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

/* -------------------------------------------------------------------------- */
/* MH_B1Part                                                                  */
/* -------------------------------------------------------------------------- */
/**
 * @brief Reference implementation from ergm: B1Part proposal.
 *
 * Kept here as a comparison point: ergm's B1Part wraps MH_CondB1Degree and adds
 * a combinatorial Hastings correction when the set of reachable targets depends
 * on the number of non-empty groups P.
 *
 * ERPM's ErpmToggleStep uses a different selection rule (uniform over all
 * groups except the current one), so no such correction is needed there.
 *
 * In particular:
 * - If the target group is sampled from an (implicit) state-dependent set
 *   (e.g., "only non-empty groups"), then a Hastings correction like MH_B1Part
 *   may be required.
 * - ERPM ToggleStep samples uniformly from all padded groups except the current
 *   one; this set has constant size (number_of_groups_total - 1), so the proposal
 *   is symmetric and logratio = 0 is appropriate.
 */
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

/* ========================================================================== */
/* MH_ErpmToggleStep                                                          */
/* ========================================================================== */
/**
 * @brief ERPM ToggleStep: move one actor from its current group to another group.
 *
 * Move definition:
 * - Pick one actor uniformly over 1..BIPARTITE.
 * - Let group_old be the actor's current group (unique neighbor under '~b1part').
 * - Pick group_new uniformly over all group vertices except group_old, including
 *   potentially empty padded groups.
 * - Reassign the actor from group_old to group_new (2 toggles).
 *
 * Encoding (2 toggles):
 * - (actor_id, group_old) off
 * - (actor_id, group_new) on
 *
 * Hastings ratio:
 * - With this selection rule, the proposal is symmetric:
 *   * actor is uniform over actors, both forward and backward;
 *   * group_new is uniform over a constant-size set (all groups except group_old),
 *     independent of the current number of non-empty groups;
 *   therefore logratio = 0.
 *
 * Failure mode:
 * - If the state violates '~b1part' (actor has no group), the proposal returns a no-move.
 *
 * Notes on group creation/destruction:
 * - Moving an actor to a currently empty group creates a new non-empty group.
 * - Moving an actor out of a singleton group destroys that group (it becomes empty).
 * - Because the destination sampling set is constant-size, these events do not
 *   require a combinatorial Hastings adjustment (unlike ergm's MH_B1Part).
 */
MH_P_FN(MH_ErpmToggleStep){

  erpm_require_undirected(nwp);

  if (MHp->ntoggles == 0) {
#if DEBUG_ERPM_PROPOSALS
    Rprintf("[ERPM][ToggleStep][INIT] setting ntoggles=2\n");
#endif
    MHp->ntoggles = 2;
    return;
  }

  const Vertex number_of_actors      = BIPARTITE;            /* actor-mode vertices: 1..BIPARTITE */
  const Vertex number_of_groups_total = N_NODES - BIPARTITE; /* padded group-mode vertices */

  Vertex actor_id = UnifVertex1to(number_of_actors);

  Vertex group_old = get_current_group_of_actor(actor_id, nwp);

  if (group_old == 0) {
#if DEBUG_ERPM_PROPOSALS
    if(DBG_seen_toggle < DBG_max_print){
      Rprintf("[ERPM][ToggleStep][WARN] actor=%d has no current group => no move\n", (int)actor_id);
      DBG_seen_toggle++;
    }
#endif
    MHp->ntoggles = 0;
    MHp->logratio = 0.0;
    return;
  }

  Vertex group_new = group_old;
  while (group_new == group_old) {
    Vertex group_rank = UnifVertex1to(number_of_groups_total); /* rank in 1..number_of_groups_total */
    group_new = BIPARTITE + group_rank;                        /* map rank -> vertex id */
  }

  MHp->ntoggles = 2;
  Mtail[0] = actor_id; Mhead[0] = group_old;
  Mtail[1] = actor_id; Mhead[1] = group_new;

  MHp->logratio = 0.0;

#if DEBUG_ERPM_PROPOSALS
  if(DBG_seen_toggle < DBG_max_print){
    Rprintf("[ERPM][ToggleStep] actor=%d | group_old=%d -> group_new=%d | ntoggles=%d | "
            "toggles={(%d,%d),(%d,%d)} | logratio=%.6f\n",
            (int)actor_id, (int)group_old, (int)group_new, (int)MHp->ntoggles,
            (int)Mtail[0], (int)Mhead[0], (int)Mtail[1], (int)Mhead[1],
            MHp->logratio);
    DBG_seen_toggle++;
  }
#endif
}

/* ========================================================================== */
/* MH_ErpmSwapStep                                                            */
/* ========================================================================== */
/**
 * @brief ERPM SwapStep: swap memberships of two actors belonging to two distinct non-empty groups.
 *
 * Move definition:
 * - Build the set of non-empty groups (IN_DEG[group] > 0).
 * - Choose two distinct groups group_i, group_j uniformly from that set.
 * - Choose actor_i uniformly among actors currently in group_i.
 * - Choose actor_j uniformly among actors currently in group_j.
 * - Swap memberships of actor_i and actor_j between group_i and group_j (4 toggles).
 *
 * Encoding (4 toggles):
 * - (actor_i, group_i) off, (actor_i, group_j) on
 * - (actor_j, group_j) off, (actor_j, group_i) on
 *
 * Hastings ratio:
 * - With the same selection rule before and after the swap, the proposal is symmetric,
 *   hence logratio = 0.
 *
 * Structural note:
 * - SwapStep preserves group sizes exactly, so it cannot connect states with different
 *   size vectors. It is intended as an intra-size mixing move.
 * - Because groups are required to be non-empty by construction here, SwapStep cannot
 *   create or destroy groups (unlike ToggleStep).
 */
MH_P_FN(MH_ErpmSwapStep){

  erpm_require_undirected(nwp);

  if (MHp->ntoggles == 0) {
#if DEBUG_ERPM_PROPOSALS
    Rprintf("[ERPM][SwapStep][INIT] setting ntoggles=4\n");
#endif
    MHp->ntoggles = 4;
    return;
  }

  const Vertex number_of_actors       = BIPARTITE;
  const Vertex number_of_groups_total = N_NODES - BIPARTITE; /* unused currently; kept for symmetry with ToggleStep */

  /* ---------------------------------------------------------------------- */
  /* 1) Build list of non-empty groups (IN_DEG[group] > 0).                  */
  /* ---------------------------------------------------------------------- */
  Vertex number_of_nonempty_groups = 0;

  for(Vertex group_id = number_of_actors + 1; group_id <= N_NODES; ++group_id){
    if(IN_DEG[group_id] > 0) number_of_nonempty_groups++;
  }

  /* Need at least two non-empty groups to swap. */
  if(number_of_nonempty_groups < 2){
    MHp->ntoggles = 0;
    MHp->logratio = 0.0;
    return;
  }

  /* Allocate group list; automatically free by R. */
  Vertex *nonempty_group_ids = (Vertex*) R_alloc((size_t)number_of_nonempty_groups, sizeof(Vertex));

  Vertex write_index = 0;
  for(Vertex group_id = number_of_actors + 1; group_id <= N_NODES; ++group_id){
    if(IN_DEG[group_id] > 0) nonempty_group_ids[write_index++] = group_id;
  }

  /* ---------------------------------------------------------------------- */
  /* 2) Choose two distinct groups uniformly from the non-empty list.        */
  /* ---------------------------------------------------------------------- */
  Vertex index_i = UnifVertex1to(number_of_nonempty_groups) - 1; /* 0..count-1 */
  Vertex index_j = index_i;
  while(index_j == index_i){
    index_j = UnifVertex1to(number_of_nonempty_groups) - 1;
  }

  Vertex group_i = nonempty_group_ids[index_i];
  Vertex group_j = nonempty_group_ids[index_j];

  /* ---------------------------------------------------------------------- */
  /* 3) Choose one actor uniformly in each chosen group.                      */
  /*    We scan actors and pick the k-th member (rank sampling).              */
  /* ---------------------------------------------------------------------- */
  Vertex degree_group_i = (Vertex)IN_DEG[group_i];
  Vertex degree_group_j = (Vertex)IN_DEG[group_j];

  /* Defensive: should not happen given the non-empty filter, but keep it safe. */
  if(degree_group_i == 0 || degree_group_j == 0){
    MHp->ntoggles = 0;
    MHp->logratio = 0.0;
    return;
  }

  Vertex rank_in_group_i = UnifVertex1to(degree_group_i); /* rank in 1..degree_group_i */
  Vertex rank_in_group_j = UnifVertex1to(degree_group_j);

  Vertex actor_i = 0, actor_j = 0;

  Vertex seen_members = 0;
  for(Vertex actor_id = 1; actor_id <= number_of_actors; ++actor_id){
    if(IS_UNDIRECTED_EDGE(actor_id, group_i)){
      if(++seen_members == rank_in_group_i){ actor_i = actor_id; break; }
    }
  }

  seen_members = 0;
  for(Vertex actor_id = 1; actor_id <= number_of_actors; ++actor_id){
    if(IS_UNDIRECTED_EDGE(actor_id, group_j)){
      if(++seen_members == rank_in_group_j){ actor_j = actor_id; break; }
    }
  }

  if(actor_i == 0 || actor_j == 0 || actor_i == actor_j){
    /* Should not happen; fail safe. */
    MHp->ntoggles = 0;
    MHp->logratio = 0.0;
    return;
  }

  /* ---------------------------------------------------------------------- */
  /* 4) Encode swap (4 toggles).                                              */
  /* ---------------------------------------------------------------------- */
  MHp->ntoggles = 4;
  Mtail[0] = actor_i; Mhead[0] = group_i;
  Mtail[1] = actor_i; Mhead[1] = group_j;
  Mtail[2] = actor_j; Mhead[2] = group_j;
  Mtail[3] = actor_j; Mhead[3] = group_i;

  MHp->logratio = 0.0;

#if DEBUG_ERPM_PROPOSALS
  if(DBG_seen_swap < DBG_max_print){
    Rprintf("[ERPM][SwapStep] group_i=%d (deg=%d), group_j=%d (deg=%d) | "
            "actor_i=%d (rank=%d) <-> actor_j=%d (rank=%d) | "
            "ntoggles=%d | toggles={(%d,%d),(%d,%d),(%d,%d),(%d,%d)} | logratio=%.6f\n",
            (int)group_i, (int)degree_group_i, (int)group_j, (int)degree_group_j,
            (int)actor_i, (int)rank_in_group_i, (int)actor_j, (int)rank_in_group_j,
            (int)MHp->ntoggles,
            (int)Mtail[0], (int)Mhead[0],
            (int)Mtail[1], (int)Mhead[1],
            (int)Mtail[2], (int)Mhead[2],
            (int)Mtail[3], (int)Mhead[3],
            MHp->logratio);
    DBG_seen_swap++;
  }
#endif

  UNUSED_VARIABLE(number_of_groups_total);
}
