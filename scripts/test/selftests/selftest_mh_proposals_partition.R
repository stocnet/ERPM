# ==============================================================================
# File    : scripts/test/selftests/selftest_mh_proposals_partition.R
# Auteur  : Jérémie Chichignoud - Cub'itech
# Purpose : Self-test for ERPM MH proposals (ErpmToggleStep / ErpmSwapStep)
#
# Goals
#   - PHASE 1 (INIT / PLUMBING)
#       * ensure InitErgmProposal.* are visible and return the expected proposal
#         descriptors (name/pkgname).
#       * ensure the proposals are actually accepted by ergm's control interface.
#
#   - PHASE 2 (SIMULATE / DYNAMICS)
#       * run short MCMC simulations using each proposal explicitly.
#       * verify invariants expected under ~ b1part partition representation:
#           - actor degrees stay == 1 (each actor belongs to exactly one group)
#           - network remains bipartite
#           - edge count stays == number of actors (n)
#
#   - PHASE 3 (MULTI-TOGGLE PROBE)
#       * observe that proposals generate multi-toggle moves:
#           - ToggleStep should produce ntoggles = 2
#           - SwapStep   should produce ntoggles = 4
#       * this phase is implemented as a "debug probe": it relies on optional C
#         logging enabled in MHproposal_partition.c (recommended).
#
# Notes
#   - This selftest does NOT validate acceptance rates or stationary behavior.
#     It only validates proposal wiring + structural invariants.
#   - If you do not see multi-toggle traces in PHASE 3, enable debug prints in
#     the C proposals (see instructions below) and recompile.
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
# Optional: ERGM patch hook (if you use it in your tree)
# ------------------------------------------------------------------------------
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
}

# ------------------------------------------------------------------------------
# Load ERPM package (dev)
# ------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Run from the package root (DESCRIPTION) with devtools available.")
}

# Wrapper ERPM (needed for build_bipartite_from_inputs in your tree)
if (!exists("build_bipartite_from_inputs", mode = "function")) {
  if (file.exists("R/erpm_build_bipartite.R")) {
    source("R/erpm_build_bipartite.R", local = FALSE)
  } else {
    stop("build_bipartite_from_inputs() missing and R/erpm_build_bipartite.R not found.")
  }
}

cat("=== SELFTEST ERPM: MH proposals (ErpmToggleStep / ErpmSwapStep) ===\n")
cat("R:", paste(R.version$major, R.version$minor, sep="."), "\n")
cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")

# ==============================================================================
# Run settings
# ==============================================================================
RUN <- list(
  phase1_init      = TRUE,
  phase2_simulate  = TRUE,
  phase3_multitog  = TRUE,

  quiet_phase1     = FALSE,
  quiet_phase2     = FALSE,

  # Simulation settings (keep small; this is a wiring/invariant test)
  burnin           = 200,
  interval         = 1,
  samplesize       = 200
)

# ==============================================================================
# Data (explicit partitions)
# ==============================================================================
partitions <- list(
  P1 = c(1,1,2,2,2,3),
  P2 = c(1,2,2,3,3,3,4,5,5,5,6,6,6,6,6)
)

.make_nodes <- function(part) {
  n <- length(part)
  set.seed(100 + n)
  data.frame(
    label  = paste0("A", seq_len(n)),
    age    = sample(20:60, n, TRUE),
    score  = round(runif(n, 0, 10), 2),
    stringsAsFactors = FALSE
  )
}

.make_nw <- function(part, nodes) {
  built <- build_bipartite_from_inputs(partition = part, nodes = nodes)
  if (is.list(built) && !is.null(built$network) && inherits(built$network, "network")) return(built$network)
  if (inherits(built, "network")) return(built)
  stop("build_bipartite_from_inputs() did not return a 'network'.")
}

# ==============================================================================
# Helpers: invariants under ~ b1part
# ==============================================================================
.get_n1 <- function(nw) {
  n1 <- network::get.network.attribute(nw, "bipartite")
  if (is.null(n1) || is.na(n1)) stop("Missing/invalid 'bipartite' network attribute.")
  as.integer(n1)
}

# actor degrees must be exactly 1 under b1part partition representation
.actor_degrees <- function(nw) {
  n1 <- .get_n1(nw)

  # Robust degree computation without relying on network.degree() exports.
  el <- network::as.edgelist(nw)
  if (is.null(el) || nrow(el) == 0L) {
    # No edges => all actor degrees 0 (will fail invariant, which is correct).
    return(rep.int(0L, n1))
  }

  # For undirected membership networks, each edge contributes 1 to each endpoint.
  # We only need actor-mode degrees: vertices 1..n1.
  tails <- as.integer(el[, 1])
  heads <- as.integer(el[, 2])

  deg1 <- integer(n1)

  # Count appearances of each actor vertex in either endpoint column.
  idx_t <- tails[tails >= 1L & tails <= n1]
  idx_h <- heads[heads >= 1L & heads <= n1]

  if (length(idx_t)) deg1[idx_t] <- deg1[idx_t] + 1L
  if (length(idx_h)) deg1[idx_h] <- deg1[idx_h] + 1L

  deg1
}

.check_invariants_b1part <- function(nw, tag = "") {
  n1 <- .get_n1(nw)
  N  <- network::network.size(nw)
  stopifnot(N >= n1 + 1L)

  # 1) bipartite attribute consistent
  if (!isTRUE(network::is.bipartite(nw))) {
    stop(sprintf("Invariant failed%s: network is not bipartite.", if (nzchar(tag)) paste0(" [", tag, "]") else ""))
  }

  # 2) actor degrees == 1
  d1 <- .actor_degrees(nw)
  if (any(d1 != 1L)) {
    bad <- which(d1 != 1L)
    msg <- sprintf(
      "Invariant failed%s: actor degrees not all 1. bad idx=%s | deg=%s",
      if (nzchar(tag)) paste0(" [", tag, "]") else "",
      paste(head(bad, 10), collapse = ","),
      paste(head(d1[bad], 10), collapse = ",")
    )
    stop(msg)
  }

  # 3) edge count == n1 (one membership edge per actor)
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

# ==============================================================================
# Helpers: building formulas in a controlled environment
# ==============================================================================
.make_f_nw <- function(nw, rhs = "edges") {
  f <- as.formula(paste0("nw ~ ", rhs))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  f
}

# ==============================================================================
# PHASE 1 — INIT / PLUMBING
# ==============================================================================
.run_phase1_init <- function(quiet = FALSE) {
  cat("\n=== PHASE 1: INIT / PLUMBING ===\n")
  if (isTRUE(quiet)) cat("  [quiet]\n")

  # Sanity: exported initializers exist
  req <- c("InitErgmProposal.ErpmToggleStep", "InitErgmProposal.ErpmSwapStep")
  for (fn in req) {
    if (!exists(fn, mode = "function")) stop("Missing function: ", fn)
  }

  part  <- partitions$P1
  nodes <- .make_nodes(part)
  nw    <- .make_nw(part, nodes)

  # Check returned descriptors
  p_toggle <- InitErgmProposal.ErpmToggleStep(nw)
  p_swap   <- InitErgmProposal.ErpmSwapStep(nw)

  if (!isTRUE(quiet)) {
    cat("  ToggleStep initializer returned:\n")
    print(p_toggle)
    cat("  SwapStep initializer returned:\n")
    print(p_swap)
  }

  stopifnot(is.list(p_toggle), identical(p_toggle$name, "ErpmToggleStep"), identical(p_toggle$pkgname, "ERPM"))
  stopifnot(is.list(p_swap),   identical(p_swap$name,   "ErpmSwapStep"),   identical(p_swap$pkgname,   "ERPM"))

  # Check invariants on the initial network
  .check_invariants_b1part(nw, tag = "initial")

  cat("  OK: proposal init descriptors + initial invariants.\n")
  invisible(TRUE)
}

# ==============================================================================
# PHASE 2 — SIMULATE / DYNAMICS (proposal wiring + invariants after MCMC)
# ==============================================================================
.sim_one <- function(nw, proposal_name, burnin, interval, samplesize, quiet = FALSE) {
  f <- .make_f_nw(nw, rhs = "edges")

  ctrl <- control.simulate.formula(
    MCMC.prop        = as.formula(paste0('~ .select("', proposal_name, '") + b1part')),
    MCMC.burnin      = burnin,
    MCMC.interval    = interval,
    MCMC.packagenames = "ERPM"
  )

  if (!isTRUE(quiet)) {
    cat(sprintf("\n--- simulate() with proposal=%s ---\n", proposal_name))
    cat(sprintf("  burnin=%d | interval=%d | (samplesize hint=%d)\n", burnin, interval, samplesize))
  }

  sim <- try(
    simulate(
      f,
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

  # simulate.formula returns a network (or list) as output.
  # Extract the simulated network robustly.
  sim_nw <- NULL
  if (inherits(sim, "network")) sim_nw <- sim
  else if (is.list(sim) && length(sim) >= 1L && inherits(sim[[1L]], "network")) sim_nw <- sim[[1L]]

  if (is.null(sim_nw)) stop("Unexpected simulate() output type; expected a 'network'.")

  .check_invariants_b1part(sim_nw, tag = paste0("post-sim:", proposal_name))

  if (!isTRUE(quiet)) {
    cat("  OK: invariants hold after MCMC.\n")
  }

  invisible(sim_nw)
}

.run_phase2_simulate <- function(quiet = FALSE) {
  cat("\n=== PHASE 2: SIMULATE / DYNAMICS ===\n")
  if (isTRUE(quiet)) cat("  [quiet]\n")

  part  <- partitions$P2
  nodes <- .make_nodes(part)
  nw    <- .make_nw(part, nodes)

  # Before MCMC
  .check_invariants_b1part(nw, tag = "pre-sim")

  # Run both proposals
  out1 <- .sim_one(
    nw, proposal_name = "ErpmToggleStep",
    burnin = RUN$burnin, interval = RUN$interval, samplesize = RUN$samplesize,
    quiet = quiet
  )

  out2 <- .sim_one(
    nw, proposal_name = "ErpmSwapStep",
    burnin = RUN$burnin, interval = RUN$interval, samplesize = RUN$samplesize,
    quiet = quiet
  )

  # Additional weak sanity: partitions changed at least sometimes (not guaranteed, but likely)
  # We only report, we do not fail if identical.
  if (!isTRUE(quiet)) {
    cat("\n--- Optional sanity: did something move? (non-failing) ---\n")
    cat("  pre edge list head:\n")
    print(utils::head(network::as.edgelist(nw)))
    cat("  post ToggleStep edge list head:\n")
    print(utils::head(network::as.edgelist(out1)))
    cat("  post SwapStep edge list head:\n")
    print(utils::head(network::as.edgelist(out2)))
  }

  cat("\nPHASE 2 OK.\n")
  invisible(TRUE)
}

# ==============================================================================
# PHASE 3 — MULTI-TOGGLE PROBE (requires C debug prints)
# ==============================================================================
.run_phase3_multitoggle_probe <- function() {
  cat("\n=== PHASE 3: MULTI-TOGGLE PROBE ===\n")
  cat("This phase is a debug probe.\n")
  cat("Enable C prints in MHproposal_partition.c (recommended):\n")
  cat("  - print MHp->ntoggles and (tail,head) pairs in MH_ErpmToggleStep / MH_ErpmSwapStep\n")
  cat("  - recompile\n\n")

  part  <- partitions$P1
  nodes <- .make_nodes(part)
  nw    <- .make_nw(part, nodes)

  # Very small runs, verbose=TRUE to surface C-side prints.
  local_ctrl <- list(burnin = 50, interval = 1, samplesize = 50)

  .sim_one(nw, "ErpmToggleStep", local_ctrl$burnin, local_ctrl$interval, local_ctrl$samplesize, quiet = FALSE)
  .sim_one(nw, "ErpmSwapStep",   local_ctrl$burnin, local_ctrl$interval, local_ctrl$samplesize, quiet = FALSE)

  cat("\nProbe done.\n")
  invisible(TRUE)
}

# ==============================================================================
# Main runner
# ==============================================================================
run_all_tests_mh_proposals <- function() {
  set.seed(42)

  if (isTRUE(RUN$phase1_init)) {
    .run_phase1_init(quiet = isTRUE(RUN$quiet_phase1))
  } else {
    cat("\n=== PHASE 1: INIT / PLUMBING ===\nSKIP\n")
  }

  if (isTRUE(RUN$phase2_simulate)) {
    .run_phase2_simulate(quiet = isTRUE(RUN$quiet_phase2))
  } else {
    cat("\n=== PHASE 2: SIMULATE / DYNAMICS ===\nSKIP\n")
  }

  if (isTRUE(RUN$phase3_multitog)) {
    .run_phase3_multitoggle_probe()
  } else {
    cat("\n=== PHASE 3: MULTI-TOGGLE PROBE ===\nSKIP\n")
  }

  invisible(TRUE)
}

if (identical(environment(), globalenv())) {
  run_all_tests_mh_proposals()
}

# ------------------------------------------------------------------------------
# Optional: disable ERGM patch hook
# ------------------------------------------------------------------------------
if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()