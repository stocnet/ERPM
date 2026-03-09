
# ==============================================================================
# File    : scripts/test/selftests/selftest_ergm_prop_swap-toggle_step.R
# Auteur : Jérémie Chichignoud - Cub'itech
# Purpose : Self-test (via erpm wrapper) for ERPM MH proposals:
#           - ErpmToggleStep (2 toggles)
#           - ErpmSwapStep   (4 toggles)
#
# Notes
#   - This script is meant to be run from the package root (DESCRIPTION present).
#   - It validates both 'wiring' (registration, simulate()) and basic properties
#     of the proposals under the 'b1part' constraint.
# ==============================================================================

options(ergm.loglik.warn_dyads = FALSE)
options(ERPM.zzz.verbose = TRUE)
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",    quietly = TRUE)) stop("Package 'ergm' requis.")
})

suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

# ------------------------------------------------------------------------------
# Optional: ERGM patch
# ------------------------------------------------------------------------------
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
}

# ------------------------------------------------------------------------------
# Load ERPM package - toremove
# ------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Run from the package root (DESCRIPTION) with devtools available.")
}

cat("=== SELFTEST ERPM: proposals via erpm() (ErpmToggleStep / ErpmSwapStep) ===\n")
cat("R:", paste(R.version$major, R.version$minor, sep="."), "\n")
cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")
cat("ERPM loaded.\n")

# ==============================================================================
# Run settings
# ==============================================================================
# RUN controls which phases execute. The default is to run everything.
# For quick iterations, disable expensive phases (simulate, disttest, connectivity).
RUN <- list(
  phase1_plumbing   = TRUE,
  phase2_simulate   = TRUE,
  phase3_disttest   = TRUE,
  phase4_fitcompare = TRUE,

  quiet_phase1     = FALSE,
  quiet_phase2     = FALSE,

  burnin           = 300,
  interval         = 1,
  samplesize       = 300,

  # Distribution probe (flat target; offset at 0 => constant loglik).
  # We measure the empirical target-choice law of ToggleStep from a fixed state.
  dist_nsteps       = 4000,
  dist_burnin       = 200,
  dist_interval     = 1,
  dist_tol_rel      = 0.20,
  dist_min_events   = 200
)

# ==============================================================================
# Data
# ==============================================================================
# Partitions are explicit integer vectors: one group label per actor.
# Groups are treated as labeled in this test script.
partitions <- list(
  P_small = c(1,1,2,2,2,3),                         # n=6
  P_mid   = c(1,2,2,3,3,3,4,5,5,5,6,6,6,6,6)        # n=15
)

# Node attributes for cov_match('bin_att') on the actor mode.
# build_bipartite_from_inputs must map these into vertex attributes.
nodes_small <- data.frame(
  label   = 1:6,
  bin_att = c("1","1","1","0","0","0"),
  stringsAsFactors = FALSE
)

nodes_mid <- data.frame(
  label   = 1:15,
  bin_att = c("1","1","1","0","0","0",  "1","0","1","0",  "1","0","1","0","1"),
  stringsAsFactors = FALSE
)

# ==============================================================================
# Helpers
# ==============================================================================

# ------------------------------------------------------------------------------
# .unwrap_network
# ------------------------------------------------------------------------------
# simulate() can return either:
#   - a 'network' directly, or
#   - a list of 'network' objects (length nsim), depending on call shape.
# This helper returns a single 'network' when possible, while never unwrapping
# an object that already inherits from 'network'.
.unwrap_network <- function(x) {
  if (inherits(x, "network")) return(x)
  if (is.list(x) && length(x) >= 1L && inherits(x[[1L]], "network")) return(x[[1L]])
  x
}

# ------------------------------------------------------------------------------
# .get_n1
# ------------------------------------------------------------------------------
# Extract actor-mode size (n1) from the 'bipartite' network attribute.
.get_n1 <- function(nw) {
  if (!inherits(nw, "network")) {
    stop(sprintf("Expected a 'network' object; got class=%s", paste(class(nw), collapse = ",")))
  }
  n1 <- network::get.network.attribute(nw, "bipartite")
  if (is.null(n1) || is.na(n1)) stop("Missing/invalid 'bipartite' network attribute.")
  as.integer(n1)
}

# ------------------------------------------------------------------------------
# .actor_degrees
# ------------------------------------------------------------------------------
# Actor degrees are computed from the edge list and must be exactly 1 under 'b1part'
# (one membership edge per actor).
.actor_degrees <- function(nw) {
  n1 <- .get_n1(nw)

  el <- network::as.edgelist(nw)
  if (is.null(el) || nrow(el) == 0L) return(rep.int(0L, n1))

  tails <- as.integer(el[, 1])
  heads <- as.integer(el[, 2])

  deg1 <- integer(n1)

  idx_t <- tails[tails >= 1L & tails <= n1]
  idx_h <- heads[heads >= 1L & heads <= n1]

  if (length(idx_t)) deg1[idx_t] <- deg1[idx_t] + 1L
  if (length(idx_h)) deg1[idx_h] <- deg1[idx_h] + 1L

  deg1
}

# ------------------------------------------------------------------------------
# .check_invariants_b1part
# ------------------------------------------------------------------------------
# Minimal structural sanity checks for a partition-as-bipartite encoding under
# the 'b1part' constraint:
#   - bipartite flag is TRUE and 'bipartite' attribute exists
#   - each actor has degree 1
#   - edgecount equals n1 (one membership edge per actor)
.check_invariants_b1part <- function(nw, tag = "") {
  if (!inherits(nw, "network")) {
    stop(sprintf("Invariant check%s: expected 'network', got class=%s",
                 if (nzchar(tag)) paste0(" [", tag, "]") else "",
                 paste(class(nw), collapse = ",")))
  }

  n1 <- .get_n1(nw)
  N  <- network::network.size(nw)

  if (N < n1 + 1L) {
    stop(sprintf("Invariant failed%s: network size N=%d < n1+1=%d.",
                 if (nzchar(tag)) paste0(" [", tag, "]") else "",
                 N, n1 + 1L))
  }

  if (!isTRUE(network::is.bipartite(nw))) {
    stop(sprintf("Invariant failed%s: network is not bipartite.",
                 if (nzchar(tag)) paste0(" [", tag, "]") else ""))
  }

  d1 <- .actor_degrees(nw)
  if (any(d1 != 1L)) {
    bad <- which(d1 != 1L)
    stop(sprintf(
      "Invariant failed%s: actor degrees not all 1. bad idx=%s | deg=%s",
      if (nzchar(tag)) paste0(" [", tag, "]") else "",
      paste(head(bad, 10), collapse = ","),
      paste(head(d1[bad], 10), collapse = ",")
    ))
  }

  m <- network::network.edgecount(nw)
  if (m != n1) {
    stop(sprintf(
      "Invariant failed%s: edgecount=%d but expected n1=%d (one membership edge per actor).",
      if (nzchar(tag)) paste0(" [", tag, "]") else "",
      m, n1
    ))
  }

  TRUE
}

# ------------------------------------------------------------------------------
# .make_translated_ergm_call
# ------------------------------------------------------------------------------
# Ask erpm() to translate 'partition ~ RHS' into an ergm() call, without fitting.
# This is used to validate that:
#   - the wrapper builds a bipartite network on the LHS
#   - constraints/terms are expanded as expected
.make_translated_ergm_call <- function(partition, rhs, nodes, verbose = FALSE) {
  if (!exists("erpm", mode = "function")) stop("erpm() not found (ERPM not loaded?).")

  user_formula <- as.formula(paste0("partition ~ ", rhs))
  environment(user_formula) <- list2env(list(partition = partition), parent = parent.frame())

  call <- erpm(
    user_formula,
    eval.call   = FALSE,
    verbose     = verbose,
    nodes       = nodes,
    dyads       = list(),
    constraints = NULL
  )

  if (!is.call(call) || !identical(call[[1L]], as.name("ergm"))) {
    stop("erpm(eval.call=FALSE) did not return an ergm() call.")
  }

  call
}

# ------------------------------------------------------------------------------
# .extract_formula_constraints_from_ergm_call
# ------------------------------------------------------------------------------
# Convenience: return translated formula + constraints argument (if present).
.extract_formula_constraints_from_ergm_call <- function(ergm_call) {
  fml <- ergm_call[[2L]]
  if (!inherits(fml, "formula")) stop("Unexpected: ergm_call[[2]] is not a formula.")

  al <- as.list(ergm_call)
  nm <- names(al)
  idx_constraints <- which(nm == "constraints")
  constraints <- if (length(idx_constraints) == 1L) al[[idx_constraints]] else NULL

  list(formula = fml, constraints = constraints)
}

# ------------------------------------------------------------------------------
# .sim_one_with_proposal
# ------------------------------------------------------------------------------
# Run a single simulate() call while forcing the proposal via MCMC.prop.
# We interpret 'burnin/interval/samplesize' as a total number of MH steps,
# and request only one final draw (nsim=1).
.sim_one_with_proposal <- function(ergm_formula, proposal_name,
                                  burnin, interval, samplesize,
                                  quiet = FALSE) {

  mcmc_prop <- as.formula(paste0('~ .select("', proposal_name, '") + b1part'))

  nsteps_total <- as.integer(burnin + samplesize * interval)
  if (nsteps_total < 1L) nsteps_total <- 1L

  ctrl <- control.simulate.formula(
    MCMC.prop         = mcmc_prop,
    MCMC.burnin       = nsteps_total,
    MCMC.interval     = 1,
    MCMC.packagenames = "ERPM"
  )

  if (!isTRUE(quiet)) {
    cat(sprintf("\n--- simulate() using proposal=%s ---\n", proposal_name))
    cat(sprintf("  requested: burnin=%d | interval=%d | samplesize=%d\n", burnin, interval, samplesize))
    cat(sprintf("  implemented: total_steps=%d (via MCMC.burnin) | nsim=1 | MCMC.interval=1\n", nsteps_total))
    cat("  MCMC.prop:", paste(deparse(mcmc_prop), collapse = " "), "\n")
  }

  sim <- try(
    simulate(
      ergm_formula,
      nsim    = 1,
      control = ctrl,
      verbose = !isTRUE(quiet)
    ),
    silent = TRUE
  )

  if (inherits(sim, "try-error")) {
    msg <- conditionMessage(attr(sim, "condition"))
    stop(sprintf("simulate() failed for proposal=%s: %s", proposal_name, msg))
  }

  sim_nw <- .unwrap_network(sim)
  if (!inherits(sim_nw, "network")) {
    stop(sprintf("Unexpected simulate() output type; expected a 'network', got class=%s",
                 paste(class(sim_nw), collapse=",")))
  }

  .check_invariants_b1part(sim_nw, tag = paste0("post-sim:", proposal_name))
  if (!isTRUE(quiet)) cat("  OK: invariants hold after simulation.\n")

  invisible(sim_nw)
}

# ------------------------------------------------------------------------------
# .get_actor_group_vec
# ------------------------------------------------------------------------------
# Decode a partition from the bipartite membership edges:
# returns an integer vector g[1..n1] where g[i] is the group vertex id (mode 2).
.get_actor_group_vec <- function(nw) {
  if (!inherits(nw, "network")) stop("Expected 'network' in .get_actor_group_vec().")

  n1 <- .get_n1(nw)
  el <- network::as.edgelist(nw)
  if (is.null(el) || nrow(el) == 0L) stop("Empty edgelist in a b1part network (invalid).")

  a <- integer(n1)
  for (k in seq_len(nrow(el))) {
    u <- el[k, 1]; v <- el[k, 2]
    if (u <= n1 && v > n1) a[u] <- v
    else if (v <= n1 && u > n1) a[v] <- u
  }
  if (any(a == 0L)) stop("Some actors have no group in edgelist decoding (invalid b1part).")
  a
}

# ------------------------------------------------------------------------------
# .decode_toggle_event
# ------------------------------------------------------------------------------
# If exactly one actor changed group between two networks, return:
#   list(i=<actor>, g_old=<group vertex>, g_new=<group vertex>)
# Otherwise return NULL.
.decode_toggle_event <- function(nw_prev, nw_next) {
  g_prev <- .get_actor_group_vec(nw_prev)
  g_next <- .get_actor_group_vec(nw_next)

  idx <- which(g_prev != g_next)
  if (length(idx) != 1L) return(NULL)

  i <- idx[1L]
  list(i = i, g_old = g_prev[i], g_new = g_next[i])
}

# ==============================================================================
# PHASE 1
# ==============================================================================
# 'Plumbing' checks:
#   - proposals are registered in ergm_proposal_table()
#   - erpm(eval.call=FALSE) produces a valid ergm() call
#   - the LHS network satisfies basic b1part invariants.
.run_phase1_plumbing <- function(quiet = FALSE) {
  cat("\n=== PHASE 1: ERPM PLUMBING ===\n")
  if (isTRUE(quiet)) cat("  [quiet]\n")

  tab <- ergm::ergm_proposal_table()
  tab_erpm <- tab[tab$Package == "ERPM", , drop = FALSE]

  if (!isTRUE(quiet)) {
    cat("  ergm proposal table (Package==ERPM):\n")
    print(tab_erpm)
  }

  has_toggle <- any(tab_erpm$Proposal == "ErpmToggleStep" & grepl("b1part", tab_erpm$Constraints))
  has_swap   <- any(tab_erpm$Proposal == "ErpmSwapStep"   & grepl("b1part", tab_erpm$Constraints))

  if (!has_toggle) stop("Missing proposal table row for ERPM::ErpmToggleStep under b1part.")
  if (!has_swap)   stop("Missing proposal table row for ERPM::ErpmSwapStep under b1part.")

  part <- partitions$P_small
  ergm_call <- .make_translated_ergm_call(part, rhs = 'cov_match("bin_att")', nodes = nodes_small, verbose = FALSE)

  x <- .extract_formula_constraints_from_ergm_call(ergm_call)
  fml <- x$formula
  constraints <- x$constraints

  if (!isTRUE(quiet)) {
    cat("\n  erpm() produced ergm() call:\n")
    cat("  ", paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "), "\n", sep = "")
    cat("  translated formula:\n")
    print(fml)
    cat("  constraints:\n")
    print(constraints)
  }

  nw0 <- eval(fml[[2L]], envir = environment(fml))
  if (!inherits(nw0, "network")) stop("Translated formula LHS is not a 'network'.")
  .check_invariants_b1part(nw0, tag = "initial")

  cat("  OK: ERPM proposals registered + erpm() translation + initial invariants.\n")
  invisible(list(ergm_formula = fml, constraints = constraints))
}

# ==============================================================================
# PHASE 2
# ==============================================================================
# Force each proposal via simulate() and check that invariants still hold.
# This is the minimal end-to-end test that the MH proposal entrypoints are
# reachable and do not break the partition encoding.
.run_phase2_simulate <- function(quiet = FALSE) {
  cat("\n=== PHASE 2: SIMULATE / INVARIANTS ===\n")
  if (isTRUE(quiet)) cat("  [quiet]\n")

  part <- partitions$P_mid
  ergm_call <- .make_translated_ergm_call(part, rhs = 'cov_match("bin_att")', nodes = nodes_mid, verbose = FALSE)
  x <- .extract_formula_constraints_from_ergm_call(ergm_call)
  fml <- x$formula

  nw0 <- eval(fml[[2L]], envir = environment(fml))
  .check_invariants_b1part(nw0, tag = "pre-sim")

  out_toggle <- .sim_one_with_proposal(
    fml, proposal_name = "ErpmToggleStep",
    burnin = RUN$burnin, interval = RUN$interval, samplesize = RUN$samplesize,
    quiet = quiet
  )

  out_swap <- .sim_one_with_proposal(
    fml, proposal_name = "ErpmSwapStep",
    burnin = RUN$burnin, interval = RUN$interval, samplesize = RUN$samplesize,
    quiet = quiet
  )

  if (!isTRUE(quiet)) {
    cat("\n--- Optional non-failing sanity: show edge list heads ---\n")
    cat("  pre  :", "\n"); print(utils::head(network::as.edgelist(nw0)))
    cat("  post ToggleStep:", "\n"); print(utils::head(network::as.edgelist(out_toggle)))
    cat("  post SwapStep  :", "\n"); print(utils::head(network::as.edgelist(out_swap)))
  }

  cat("\nPHASE 2 OK.\n")
  invisible(TRUE)
}

# ==============================================================================
# PHASE 3
# ==============================================================================
# Distribution sanity test for ToggleStep under a flat target:
#   - build a fixed state nw_fixed
#   - repeatedly draw one-step proposals from the same state
#   - check that the choice of new group is consistent with uniformity over
#     targets excluding the current group.
.run_phase3_disttest <- function() {
  cat("\n=== PHASE 3: PROPOSAL DISTRIBUTION TEST (ToggleStep) ===\n")
  cat("Target: flat via offset(cov_match('bin_att'), 0). We test proposal law at a fixed state.\n")

  part <- partitions$P_mid

  ergm_call <- .make_translated_ergm_call(part, rhs = 'cov_match("bin_att")', nodes = nodes_mid, verbose = FALSE)
  x <- .extract_formula_constraints_from_ergm_call(ergm_call)
  fml <- x$formula

  nw0 <- eval(fml[[2L]], envir = environment(fml))
  .check_invariants_b1part(nw0, tag = "dist-pre")

  n1 <- .get_n1(nw0)
  G  <- network::network.size(nw0) - n1
  if (G < 2L) stop("Need at least 2 groups for ToggleStep distribution test.")

  mcmc_prop <- ~ .select("ErpmToggleStep") + b1part

  one_step_ctrl <- control.simulate.formula(
    MCMC.prop         = mcmc_prop,
    MCMC.burnin       = 1,
    MCMC.interval     = 1,
    MCMC.packagenames = "ERPM"
  )

  nw_cur <- nw0
  if (RUN$dist_burnin > 0) {
    cat(sprintf("Burn-in chain (to pick a typical fixed state): %d steps\n", RUN$dist_burnin))
    for (b in seq_len(RUN$dist_burnin)) {
      out <- simulate(nw_cur ~ offset(cov_match("bin_att"), 0), nsim = 1, control = one_step_ctrl, verbose = FALSE)
      nw_cur <- .unwrap_network(out)
      .check_invariants_b1part(nw_cur, tag = "dist-burnin-chain")
    }
  }
  nw_fixed <- nw_cur
  .check_invariants_b1part(nw_fixed, tag = "dist-fixed")

  fml_flat <- nw_fixed ~ offset(cov_match("bin_att"), 0)

  cat(sprintf("Collecting %d independent 1-step moves from a fixed state...\n", RUN$dist_nsteps))

  counts <- matrix(0L, nrow = G, ncol = G)
  rownames(counts) <- paste0("g", seq_len(G))
  colnames(counts) <- paste0("g", seq_len(G))

  recorded <- 0L
  for (t in seq_len(RUN$dist_nsteps)) {
    out <- simulate(fml_flat, nsim = 1, control = one_step_ctrl, verbose = FALSE)
    nw_next <- .unwrap_network(out)
    .check_invariants_b1part(nw_next, tag = "dist-step")

    ev <- .decode_toggle_event(nw_fixed, nw_next)
    if (!is.null(ev)) {
      old_k <- ev$g_old - n1
      new_k <- ev$g_new - n1
      if (old_k >= 1L && old_k <= G && new_k >= 1L && new_k <= G) {
        counts[old_k, new_k] <- counts[old_k, new_k] + 1L
        recorded <- recorded + 1L
      }
    }
  }

  cat(sprintf("Recorded toggle events: %d (out of %d reps)\n", recorded, RUN$dist_nsteps))
  if (recorded < 200L) stop("Too few recorded events; increase dist_nsteps.")

  minN <- RUN$dist_min_events
  alpha <- 0.001

  tested <- 0L
  for (old_k in seq_len(G)) {
    v <- counts[old_k, ]
    v[old_k] <- 0L
    N <- sum(v)
    if (N < minN) next

    expected <- rep(N / (G - 1L), G - 1L)
    obs <- v[-old_k]

    chisq <- sum((obs - expected)^2 / expected)
    df <- (G - 1L) - 1L
    pval <- stats::pchisq(chisq, df = df, lower.tail = FALSE)

    if (!is.finite(pval)) {
      stop(sprintf("Chi-square test produced non-finite p-value for old_group=%d.", old_k))
    }
    if (pval < alpha) {
      stop(sprintf(
        "Distribution test failed for old_group=%d: chi^2=%.3f df=%d p=%.3g (alpha=%.3g)\nCounts=%s",
        old_k, chisq, df, pval, alpha, paste(v, collapse = ",")
      ))
    }

    tested <- tested + 1L
  }

  if (tested == 0L) {
    stop("No old-group had enough events for uniformity test (increase dist_nsteps or lower dist_min_events).")
  }

  cat(sprintf("OK: target choice consistent with uniformity at fixed state (tested groups=%d, alpha=%.3g).\n",
              tested, alpha))
  invisible(TRUE)
}

# ==============================================================================
# PHASE 4
# ==============================================================================
# Connectivity / ergodicity test for SwapStep only.
# This is a structural test, independent of the target distribution:
#   - enumerate all assignments {1..G}^n on a tiny case
#   - connect two states if one SwapStep transforms one into the other
#   - run BFS and count connected components
# Expected outcome:
#   SwapStep preserves the group-size vector, so the state space splits into
#   components indexed by that size vector (not ergodic globally).
.run_phase4_connectivity_swap <- function() {
  cat("\n=== PHASE 4: CONNECTIVITY TEST (SwapStep only) ===\n")

  n <- 6L
  G <- 4L
  cat(sprintf("Test case: n=%d actors | G=%d labeled groups\n", n, G))

  cat("Enumerating all partitions...\n")
  grids <- replicate(n, seq_len(G), simplify = FALSE)
  all_states <- as.matrix(do.call(expand.grid, grids))

  Nstates <- nrow(all_states)
  cat(sprintf("Total states: %d (= G^n)\n", Nstates))

  size_vec <- function(state) {
    tab <- tabulate(state, nbins = G)
    as.integer(tab)
  }
  sizes <- t(apply(all_states, 1, size_vec))

  can_swap <- function(s1, s2) {
    diff_idx <- which(s1 != s2)
    if (length(diff_idx) != 2L) return(FALSE)

    i <- diff_idx[1L]
    j <- diff_idx[2L]

    g1_i <- s1[i]; g1_j <- s1[j]
    g2_i <- s2[i]; g2_j <- s2[j]

    (g1_i != g1_j) &&
      (g2_i == g1_j) &&
      (g2_j == g1_i)
  }

  cat("Building adjacency graph under one SwapStep...\n")
  adj <- vector("list", Nstates)
  for (a in seq_len(Nstates)) {
    s1 <- all_states[a, ]
    neigh <- integer(0)
    for (b in seq_len(Nstates)) {
      if (a == b) next
      if (can_swap(s1, all_states[b, ])) neigh <- c(neigh, b)
    }
    adj[[a]] <- neigh
  }

  cat("Computing connected components (BFS)...\n")
  visited <- rep(FALSE, Nstates)
  comp_id <- integer(Nstates)
  comp <- 0L

  for (i in seq_len(Nstates)) {
    if (visited[i]) next

    comp <- comp + 1L
    queue <- i
    visited[i] <- TRUE
    comp_id[i] <- comp

    while (length(queue) > 0L) {
      v <- queue[1L]
      queue <- queue[-1L]
      for (u in adj[[v]]) {
        if (!visited[u]) {
          visited[u] <- TRUE
          comp_id[u] <- comp
          queue <- c(queue, u)
        }
      }
    }
  }

  cat(sprintf("Number of connected components under SwapStep: %d\n", comp))

  cat("Analyzing components vs size vectors...\n")
  size_str <- apply(sizes, 1, paste, collapse = "-")
  comp_sizes <- split(size_str, comp_id)
  comp_size_patterns <- lapply(comp_sizes, unique)
  all_size_patterns <- unique(size_str)

  cat(sprintf("Distinct size vectors in full space: %d\n", length(all_size_patterns)))

  for (k in seq_along(comp_size_patterns)) {
    cat(sprintf(
      "Component %d: %d states | distinct size vectors inside = %d\n",
      k,
      length(comp_sizes[[k]]),
      length(comp_size_patterns[[k]])
    ))
  }

  cat("\nConclusion:\n")
  cat("  SwapStep preserves group sizes, so it cannot connect states with different size vectors.\n")
  cat("  The space splits into components indexed by group-size configuration.\n")
  cat("  Therefore SwapStep alone is not ergodic on the full partition space.\n")

  invisible(list(
    n = n,
    G = G,
    Nstates = Nstates,
    n_components = comp,
    comp_id = comp_id,
    size_vectors = sizes
  ))
}

# ==============================================================================
# Main runner
# ==============================================================================
# Run phases according to RUN. Designed to be source-able and also runnable as
# a standalone script via Rscript.
run_all_tests_ergm_prop_swap_toggle <- function() {
  set.seed(42)

  if (isTRUE(RUN$phase1_plumbing)) {
    .run_phase1_plumbing(quiet = isTRUE(RUN$quiet_phase1))
  } else {
    cat("\n=== PHASE 1: ERPM PLUMBING ===\nSKIP\n")
  }

  if (isTRUE(RUN$phase2_simulate)) {
    .run_phase2_simulate(quiet = isTRUE(RUN$quiet_phase2))
  } else {
    cat("\n=== PHASE 2: SIMULATE / INVARIANTS ===\nSKIP\n")
  }

  if (isTRUE(RUN$phase3_disttest)) {
    .run_phase3_disttest()
  } else {
    cat("\n=== PHASE 3: PROPOSAL DISTRIBUTION TEST ===\nSKIP\n")
  }

  if (isTRUE(RUN$phase4_fitcompare)) {
    .run_phase4_connectivity_swap()
  } else {
    cat("\n=== PHASE 4: CONNECTIVITY TEST ===\nSKIP\n")
  }

  invisible(TRUE)
}

if (identical(environment(), globalenv())) {
  run_all_tests_ergm_prop_swap_toggle()
}

# ------------------------------------------------------------------------------
# Optional: disable ERGM patch 
# ------------------------------------------------------------------------------
if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()