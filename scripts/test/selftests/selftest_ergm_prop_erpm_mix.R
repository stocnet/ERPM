# ==============================================================================
# File    : scripts/test/selftests/selftest_ergm_prop_erpm_mix.R
# Auteur : Jérémie Chichignoud - Cub'itech
# Purpose : Self-test (via erpm wrapper) for ERPM MH proposal:
#           - ErpmMix (mix of ErpmToggleStep / ErpmSwapStep under ~b1part)
#
# File purpose
#   - PHASE 1 (SIMULATE) : simulate() with ErpmMix + b1part invariants.
#   - PHASE 2 (PROBE)    : "1-step" probe: observe toggle + swap type moves
#                          (via simulate() with burnin=1 at each iteration).
#   - PHASE 3 (FIT)      : fits via erpm() (log_factorial_sizes, cov_match, dyadcov)
#                          with a dedicated control.ergm.
#   - PHASE 4 (DIST)     : "ergodicity + series" sanity under a flat target via simulate()
#                          (not erpm(save_networks), too fragile depending on versions).
#
# ==============================================================================

# --------------------------------------------------------------------------------------
# Init
# --------------------------------------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

options(ergm.loglik.warn_dyads = FALSE)
options(ERPM.zzz.verbose = TRUE)

suppressPackageStartupMessages({
  if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' required.")
  if (!requireNamespace("ergm",    quietly = TRUE)) stop("Package 'ergm' required.")
})
suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

# Optional ERGM patch
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
}

# Load ERPM
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Run from the package root (DESCRIPTION) with devtools available.")
}

cat("=== SELFTEST ERPM: ErpmMix proposal via erpm() / simulate() ===\n")
cat("R:", paste(R.version$major, R.version$minor, sep="."), "\n")
cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")
cat("ERPM loaded.\n")

# ======================================================================================
# Settings
# ======================================================================================
RUN <- list(
  phase1_simulate_invariants = TRUE,
  phase2_probe_steps         = TRUE,
  phase3_erpm_fits           = TRUE,
  phase4_dist_sanity         = TRUE,

  quiet = TRUE,

  # Phase 1 (simulate)
  burnin     = 400,
  interval   = 1,
  samplesize = 400,

  # Phase 2 (probe)
  probe_nsteps     = 200,
  probe_min_toggle = 25,
  probe_min_swap   = 25,

  # Phase 3 (fits)
  fit_seed          = 123,
  fit_burnin        = 50000,
  fit_interval      = 20,
  fit_samplesize    = 300000,
  fit_maxit         = 40,

  # Phase 4 (dist)
  dist_seed   = 777,
  dist_steps  = 20000,
  dist_burnin = 50000,
  dist_alpha  = 1e-3
)

# ======================================================================================
# Data for test
# ======================================================================================
partition_mid <- c(1,2,2,3,3,3,4,5,5,5,6,6,6,6,6)  # n=15

nodes_mid <- data.frame(
  label   = 1:15,
  bin_att = c("1","1","1","0","0","0",  "1","0","1","0",  "1","0","1","0","1"),
  stringsAsFactors = FALSE
)

Z1_mid <- matrix(
  c(
    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
    1,0,1,2,3,4,5,6,7,8, 9,10,11,12,13,
    2,1,0,1,2,3,4,5,6,7, 8, 9,10,11,12,
    3,2,1,0,1,2,3,4,5,6, 7, 8, 9,10,11,
    4,3,2,1,0,1,2,3,4,5, 6, 7, 8, 9,10,
    5,4,3,2,1,0,1,2,3,4, 5, 6, 7, 8, 9,
    6,5,4,3,2,1,0,1,2,3, 4, 5, 6, 7, 8,
    7,6,5,4,3,2,1,0,1,2, 3, 4, 5, 6, 7,
    8,7,6,5,4,3,2,1,0,1, 2, 3, 4, 5, 6,
    9,8,7,6,5,4,3,2,1,0, 1, 2, 3, 4, 5,
    10,9,8,7,6,5,4,3,2,1, 0, 1, 2, 3, 4,
    11,10,9,8,7,6,5,4,3,2, 1, 0, 1, 2, 3,
    12,11,10,9,8,7,6,5,4,3, 2, 1, 0, 1, 2,
    13,12,11,10,9,8,7,6,5,4, 3, 2, 1, 0, 1,
    14,13,12,11,10,9,8,7,6,5, 4, 3, 2, 1, 0
  ),
  nrow = 15, ncol = 15, byrow = TRUE
)
storage.mode(Z1_mid) <- "double"
dyads_mid <- list(Z1 = Z1_mid)

# ======================================================================================
# Helpers 
# ======================================================================================
.unwrap_network <- function(x) {
  if (inherits(x, "network")) return(x)
  if (is.list(x) && length(x) >= 1L && inherits(x[[1L]], "network")) return(x[[1L]])
  x
}

.get_n1 <- function(nw) {
  if (!inherits(nw, "network")) stop("Expected a 'network' object.")
  n1 <- network::get.network.attribute(nw, "bipartite")
  if (is.null(n1) || is.na(n1)) stop("Missing/invalid 'bipartite' network attribute.")
  as.integer(n1)
}

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
    stop(sprintf("Invariant failed%s: edgecount=%d but expected n1=%d.",
                 if (nzchar(tag)) paste0(" [", tag, "]") else "",
                 m, n1))
  }

  TRUE
}

.make_translated_ergm_call <- function(partition, rhs, nodes, dyads, verbose = FALSE) {
  if (!exists("erpm", mode = "function")) stop("erpm() not found (ERPM not loaded?).")

  user_formula <- as.formula(paste0("partition ~ ", rhs))
  environment(user_formula) <- list2env(list(partition = partition), parent = parent.frame())

  call <- erpm(
    user_formula,
    eval.call   = FALSE,
    verbose     = verbose,
    nodes       = nodes,
    dyads       = dyads,
    constraints = NULL
  )
  if (!is.call(call) || !identical(call[[1L]], as.name("ergm"))) {
    stop("erpm(eval.call=FALSE) did not return an ergm() call.")
  }
  call
}

.extract_constraints_from_ergm_call <- function(ergm_call) {
  al <- as.list(ergm_call)
  nm <- names(al)
  i <- which(nm == "constraints")
  if (length(i) == 1L) al[[i]] else NULL
}

.get_actor_group_vec <- function(nw) {
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

.classify_step <- function(nw_prev, nw_next) {
  g0 <- .get_actor_group_vec(nw_prev)
  g1 <- .get_actor_group_vec(nw_next)

  idx <- which(g0 != g1)
  if (length(idx) == 1L) return("toggle")
  if (length(idx) != 2L) return("other")

  i <- idx[1L]; j <- idx[2L]
  if (g1[i] == g0[j] && g1[j] == g0[i]) return("swap")
  "other"
}

# ======================================================================================
# Build translated base objects once
# ======================================================================================
ergm_call <- .make_translated_ergm_call(
  partition = partition_mid,
  rhs       = 'cov_match("bin_att")',
  nodes     = nodes_mid,
  dyads     = dyads_mid,
  verbose   = FALSE
)

fml <- ergm_call[[2L]]
constraints <- .extract_constraints_from_ergm_call(ergm_call)
if (is.null(constraints)) stop("Translated ergm() call has no constraints; expected ~b1part.")

nw0 <- eval(fml[[2L]], envir = environment(fml))
.check_invariants_b1part(nw0, tag = "initial")

# ======================================================================================
# PHASE 1: simulate() with ErpmMix + invariants (ctrl local)
# ======================================================================================
run_phase1_simulate_invariants <- function() {
  cat("\n=== PHASE 1: SIMULATE / INVARIANTS (ErpmMix) ===\n")

  ctrl <- control.simulate.formula(
    MCMC.prop         = ~ .select("ErpmMix"),
    MCMC.packagenames = "ERPM",
    MCMC.burnin       = RUN$burnin,
    MCMC.interval     = RUN$interval
  )

  # Flat target: stats=0 (here edges in offset=0).
  fml_flat <- nw0 ~ offset(edges, 0)

  sim1 <- simulate(
    fml_flat,
    nsim        = RUN$samplesize,
    constraints = constraints,
    control     = ctrl,
    verbose     = !RUN$quiet
  )

  # nsim>1 -> network.series
  nets <- try(as.list(sim1), silent = TRUE)
  if (inherits(nets, "try-error") || is.null(nets) || length(nets) < 1L) {
    stop("PHASE 1 failed: simulate() did not return a network.series.")
  }

  # Check invariants on a few samples (start/middle/end)
  idx <- unique(pmax(1L, pmin(length(nets), c(1L, floor(length(nets)/2), length(nets)))))
  for (k in idx) .check_invariants_b1part(nets[[k]], tag = sprintf("phase1-sample-%d", k))

  cat("PHASE 1 OK.\n")
  invisible(TRUE)
}

# ======================================================================================
# PHASE 2: one-step probe 
# ======================================================================================
run_phase2_probe_steps <- function() {
  cat("\n=== PHASE 2: ONE-STEP PROBE (ErpmMix emits toggle + swap) ===\n")
  cat(sprintf("Probe steps: %d | expected min: toggle>=%d, swap>=%d\n",
              RUN$probe_nsteps, RUN$probe_min_toggle, RUN$probe_min_swap))

  ctrl <- control.simulate.formula(
    MCMC.prop         = ~ .select("ErpmMix"),
    MCMC.packagenames = "ERPM",
    MCMC.burnin       = 1,
    MCMC.interval     = 1
  )

  nw_fixed <- nw0
  .check_invariants_b1part(nw_fixed, tag = "probe-fixed")

  count_toggle <- 0L
  count_swap   <- 0L
  count_other  <- 0L

  for (t in seq_len(RUN$probe_nsteps)) {
    out <- simulate(
      nw_fixed ~ offset(edges, 0),
      nsim        = 1,
      constraints = constraints,
      control     = ctrl,
      verbose     = FALSE
    )
    nw_next <- .unwrap_network(out)
    .check_invariants_b1part(nw_next, tag = "probe-step")

    typ <- .classify_step(nw_fixed, nw_next)
    if (typ == "toggle") count_toggle <- count_toggle + 1L
    else if (typ == "swap") count_swap <- count_swap + 1L
    else count_other <- count_other + 1L
  }

  cat(sprintf("Observed: toggle=%d | swap=%d | other=%d\n", count_toggle, count_swap, count_other))

  if (count_toggle < RUN$probe_min_toggle)
    stop(sprintf("Too few toggle-like steps observed (%d < %d).", count_toggle, RUN$probe_min_toggle))
  if (count_swap < RUN$probe_min_swap)
    stop(sprintf("Too few swap-like steps observed (%d < %d).", count_swap, RUN$probe_min_swap))

  cat("PHASE 2 OK.\n")
  invisible(list(toggle = count_toggle, swap = count_swap, other = count_other))
}

# ======================================================================================
# PHASE 3: fit tests via erpm() 
# ======================================================================================
.fit_one <- function(rhs) {
  user_formula <- as.formula(paste0("partition_mid ~ ", rhs))
  environment(user_formula) <- list2env(list(partition_mid = partition_mid), parent = parent.frame())

  ctrl <- control.ergm(
    MCMC.prop         = ~ .select("ErpmMix"),
    MCMC.packagenames = "ERPM",
    MCMC.burnin       = RUN$fit_burnin,
    MCMC.interval     = RUN$fit_interval,
    MCMC.samplesize   = RUN$fit_samplesize,
    MCMLE.maxit       = RUN$fit_maxit,
    parallel          = 0
  )

  erpm(
    user_formula,
    eval.call = TRUE,
    verbose   = !RUN$quiet,
    nodes     = nodes_mid,
    dyads     = dyads_mid,
    control   = ctrl
  )
}

run_phase3_erpm_fits <- function() {
  cat("\n=== PHASE 3: FIT TESTS via erpm() (log_factorial_sizes, cov_match, dyadcov) ===\n")
  set.seed(RUN$fit_seed)

  fit_log_factorial_sizes <- .fit_one("log_factorial_sizes")
  print(summary(fit_log_factorial_sizes))
  cat("  OK: fit log_factorial_sizes\n")

  # NOTE: adapt these arguments if your current API does not accept clique_size/normalized/normalize
  fit_covmatch <- .fit_one('cov_match("bin_att", clique_size=2, normalized="by_group")')
  print(summary(fit_covmatch))
  cat('  OK: fit cov_match("bin_att", k=2, by_group)\n')

  fit_dyadcov <- .fit_one('dyadcov("Z1", clique_size=2, normalize="global")')
  print(summary(fit_dyadcov))
  cat('  OK: fit dyadcov("Z1", k=2, global)\n')

  cat("PHASE 3 OK.\n")
  invisible(list(log_factorial_sizes = fit_log_factorial_sizes, covmatch = fit_covmatch, dyadcov = fit_dyadcov))
}

# ======================================================================================
# PHASE 4: distribution sanity (robust) via simulate()
#   Checks:
#     (A) ~b1part invariants held for the whole series
#     (B) the chain actually moves
#     (C) the mix does produce both "toggle" and "swap" moves
#     (D) (optional) "flat" uniformity in TOGGLE-ONLY (more reasonable)
# ======================================================================================
run_phase4_dist_sanity <- function() {
  cat("\n=== PHASE 4: DISTRIBUTION SANITY (ErpmMix) [simulate()] ===\n")
  set.seed(RUN$dist_seed)

  ctrl_mix <- control.simulate.formula(
    MCMC.prop         = ~ .select("ErpmMix"),
    MCMC.packagenames = "ERPM",
    MCMC.burnin       = RUN$dist_burnin,
    MCMC.interval     = 1
  )

  fml_flat <- nw0 ~ offset(edges, 0)

  sim <- simulate(
    fml_flat,
    nsim        = RUN$dist_steps,
    constraints = constraints,
    control     = ctrl_mix,
    verbose     = !RUN$quiet
  )

  nets <- try(as.list(sim), silent = TRUE)
  if (inherits(nets, "try-error") || is.null(nets) || length(nets) < 100L) {
    stop("Phase 4 failed: simulate() did not return a long enough network.series.")
  }

  n1 <- .get_n1(nets[[1L]])
  N  <- network::network.size(nets[[1L]])
  Gbig <- N - n1
  if (Gbig < 2L) stop("Need at least 2 groups for distribution sanity.")

  actors_to_track <- c(1L, 2L, 3L)

  P_active <- integer(length(nets))
  sizes_sd <- numeric(length(nets))
  gtraj <- matrix(NA_integer_, nrow = length(nets), ncol = length(actors_to_track),
                  dimnames = list(NULL, paste0("a", actors_to_track)))

  n_toggle <- 0L
  n_swap   <- 0L
  n_other  <- 0L

  nw_prev <- NULL
  for (t in seq_along(nets)) {
    nw_cur <- nets[[t]]
    .check_invariants_b1part(nw_cur, tag = sprintf("dist-sample-%d", t))

    gvec <- .get_actor_group_vec(nw_cur) - n1  # groups 1..Gbig
    tab <- tabulate(gvec, nbins = Gbig)
    P_active[t] <- sum(tab > 0L)
    sizes_sd[t] <- stats::sd(tab)

    for (k in seq_along(actors_to_track)) gtraj[t, k] <- gvec[actors_to_track[k]]

    if (!is.null(nw_prev)) {
      typ <- .classify_step(nw_prev, nw_cur)
      if (typ == "toggle") n_toggle <- n_toggle + 1L
      else if (typ == "swap") n_swap <- n_swap + 1L
      else n_other <- n_other + 1L
    }
    nw_prev <- nw_cur
  }

  cat(sprintf("  series length=%d | actors tracked=%s\n",
              length(nets), paste(actors_to_track, collapse = ",")))
  cat(sprintf("  move counts (between consecutive samples): toggle=%d | swap=%d | other=%d\n",
              n_toggle, n_swap, n_other))

  changes_per_actor <- integer(ncol(gtraj))
  for (k in seq_len(ncol(gtraj))) {
    changes_per_actor[k] <- sum(gtraj[-1L, k] != gtraj[-nrow(gtraj), k], na.rm = TRUE)
  }
  cat(sprintf("  tracked actor changes: %s\n",
              paste(paste0(colnames(gtraj), "=", changes_per_actor), collapse = " | ")))

  if (all(changes_per_actor == 0L)) stop("Phase 4 failed: tracked actors never moved; chain looks stuck.")
  if (n_toggle == 0L) stop("Phase 4 failed: no toggle-like move observed in the series.")
  if (n_swap   == 0L) stop("Phase 4 failed: no swap-like move observed in the series.")

  n_unique_P <- length(unique(P_active))
  cat(sprintf("  active groups P: min=%d max=%d unique=%d\n",
              min(P_active), max(P_active), n_unique_P))
  if (n_unique_P < 2L) {
    warning("Phase 4 warning: P_active never changed; toggle may not have hit empty/size-1 groups in this run.")
  }

  cat(sprintf("  group-size sd: mean=%.3f | min=%.3f | max=%.3f\n",
              mean(sizes_sd), min(sizes_sd), max(sizes_sd)))
  cat("PHASE 4 OK (ErpmMix): invariants + movement + toggle&swap present.\n")

  # ----------------------------
  # OPTIONAL: uniformity (flat target) in TOGGLE-ONLY
  # ----------------------------
  cat("\n  [OPTIONAL] TOGGLE-ONLY uniformity sanity (flat target)\n")
  ctrl_toggle <- control.simulate.formula(
    MCMC.prop         = ~ .select("ErpmToggleStep"),
    MCMC.packagenames = "ERPM",
    MCMC.burnin       = RUN$dist_burnin,
    MCMC.interval     = 1
  )

  sim_t <- simulate(
    fml_flat,
    nsim        = RUN$dist_steps,
    constraints = constraints,
    control     = ctrl_toggle,
    verbose     = FALSE
  )
  nets_t <- try(as.list(sim_t), silent = TRUE)

  counts <- NULL
  if (!inherits(nets_t, "try-error") && length(nets_t) >= 200L) {
    a <- actors_to_track[1L]
    counts <- integer(Gbig)
    for (nw_cur in nets_t) {
      gvec <- .get_actor_group_vec(nw_cur) - n1
      counts[gvec[a]] <- counts[gvec[a]] + 1L
    }
    obs <- counts
    Nobs <- sum(obs)
    exp <- rep(Nobs / Gbig, Gbig)
    chisq <- sum((obs - exp)^2 / exp)
    df <- Gbig - 1L
    pval <- stats::pchisq(chisq, df = df, lower.tail = FALSE)
    cat(sprintf("    toggle-only actor %d: chisq=%.2f df=%d p=%.3g\n", a, chisq, df, pval))
    if (!is.finite(pval) || pval < RUN$dist_alpha) {
      warning(sprintf("Toggle-only uniformity sanity: p=%.3g < alpha=%.3g.", pval, RUN$dist_alpha))
    }
  } else {
    warning("Toggle-only uniformity sanity skipped: too few samples returned by simulate().")
  }

  invisible(list(
    mix = list(toggle = n_toggle, swap = n_swap, other = n_other,
               P_active = P_active, sizes_sd = sizes_sd, gtraj = gtraj),
    toggle_only = list(counts_actor1 = counts)
  ))
}

# ======================================================================================
# Main run
# ======================================================================================
run_all_tests_erpm_mix <- function() {
  set.seed(1)
  cat("\n=== TEST ERPM: ErpmMix (toggle+swap) under ~b1part ===\n")

  if (isTRUE(RUN$phase1_simulate_invariants)) run_phase1_simulate_invariants()
  else cat("\n=== PHASE 1 ===\nSKIP\n")

  if (isTRUE(RUN$phase2_probe_steps)) run_phase2_probe_steps()
  else cat("\n=== PHASE 2 ===\nSKIP\n")

  if (isTRUE(RUN$phase3_erpm_fits)) run_phase3_erpm_fits()
  else cat("\n=== PHASE 3 ===\nSKIP\n")

  if (isTRUE(RUN$phase4_dist_sanity)) run_phase4_dist_sanity()
  else cat("\n=== PHASE 4 ===\nSKIP\n")

  cat("\nSELFTEST OK: ErpmMix wiring + invariants + probe + fits + dist sanity.\n")
  invisible(TRUE)
}

if (identical(environment(), globalenv())) {
  run_all_tests_erpm_mix()
}

# Optional ERGM patch: disable
if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()