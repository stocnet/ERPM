# ==============================================================================
# File    : scripts/test/selftests/selftest_ergm_prop_erpm_mix_dynamic.R
# Purpose : Self-test (via erpm wrapper) for ERPM MH proposal ErpmMix under ~b1part,
#           with dynamic moves/weights decoded by InitErgmProposal.ErpmMix.
# Run     : Rscript scripts/test/selftests/selftest_ergm_prop_erpm_mix_dynamic.R
#
#
# What this selftest enforces:
#   1) "Real chain": each PROBE is a single MCMC simulation (simulate() called
#      ONCE), then we classify transitions between successive states of the chain.
#
#   2) Requested proportions (on a flat target):
#      - we use any term offset with coefficient fixed at 0
#        => contribution to the log-likelihood always zero => acceptance ~1 (if proposal symmetric)
#        => accepted moves ≈ proposed moves.
#      - default (no args)  : ~ 2/3 toggle, 1/3 swap
#      - custom 4:1         : ~ 4/5 toggle, 1/5 swap
#      - invalid args       : canonical fallback (2/1)
#
# Important (metric):
#   - Here we measure *accepted* moves (difference between successive states).
#   - Since the objective is moves/weights decoding, we force a flat target via offset(...,0)
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

if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
}

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Run from the package root (DESCRIPTION) with devtools available.")
}

cat("=== SELFTEST ERPM: ErpmMix (dynamic args) under ~b1part (REAL CHAIN) ===\n")
cat("R:", paste(R.version$major, R.version$minor, sep="."), "\n")
cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")
cat("ERPM loaded.\n")

# ======================================================================================
# Settings
# ======================================================================================
RUN <- list(
  seed1 = 123,
  seed2 = 456,
  seed3 = 789,
  seed_fit = 1526,

  # Phase 1 probes: number of transitions (differences between successive saved states).
  probe_nsteps = 1200L,

  # Tolerances on observed toggle proportion (among toggle+swap, excluding 'stay'):
  tol_default = 0.06,   # expected 2/3
  tol_4_1     = 0.05,   # expected 4/5
  tol_invalid = 0.06    # expected 2/3
)

# ======================================================================================
# Data
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
.as_network_list <- function(x) {
  if (inherits(x, "network")) return(list(x))
  if (is.list(x) && length(x) >= 1L && all(vapply(x, inherits, logical(1), "network"))) return(x)
  if (is.list(x) && length(x) >= 1L && inherits(x[[1L]], "network")) return(x)
  stop("simulate() did not return a network or a list of networks.")
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
  n1 <- .get_n1(nw)
  N  <- network::network.size(nw)

  if (!isTRUE(network::is.bipartite(nw))) stop("Not bipartite ", tag)
  if (isTRUE(network::is.directed(nw))) stop("Directed network (expected undirected) ", tag)
  if (N < n1 + 1L) stop("N too small ", tag)

  el <- network::as.edgelist(nw)
  if (!is.null(el) && nrow(el) > 0L) {
    u <- as.integer(el[, 1]); v <- as.integer(el[, 2])
    ok <- (u <= n1 & v > n1) | (v <= n1 & u > n1)
    if (!all(ok)) stop("Found non-bipartite edges (actor-actor or group-group) ", tag)
  }

  d1 <- .actor_degrees(nw)
  if (any(d1 != 1L)) stop("Actor degrees != 1 ", tag)

  # Under b1part: each actor has degree 1 => edgecount == n1 
  m <- network::network.edgecount(nw)
  if (m != n1) stop("edgecount != n1 ", tag)

  TRUE
}

.get_actor_group_vec <- function(nw) {
  n1 <- .get_n1(nw)
  el <- network::as.edgelist(nw)
  if (is.null(el) || nrow(el) == 0L) stop("Empty edgelist (invalid).")

  a <- integer(n1)
  for (k in seq_len(nrow(el))) {
    u <- as.integer(el[k, 1]); v <- as.integer(el[k, 2])
    if (u <= n1 && v > n1) a[u] <- v
    else if (v <= n1 && u > n1) a[v] <- u
  }
  if (any(a == 0L)) stop("Some actors have no group (invalid).")
  a
}

.classify_step <- function(nw_prev, nw_next) {
  g0 <- .get_actor_group_vec(nw_prev)
  g1 <- .get_actor_group_vec(nw_next)

  idx <- which(g0 != g1)
  if (length(idx) == 0L) return("stay")    # rejected move (or no-op)
  if (length(idx) == 1L) return("toggle")  # one actor changed group

  if (length(idx) != 2L) return("other")   # should not happen if proposal correct

  i <- idx[1L]; j <- idx[2L]
  if (g1[i] == g0[j] && g1[j] == g0[i]) return("swap")
  "other"
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

# ======================================================================================
# Build base objects once (network + constraints)
# ======================================================================================
set.seed(RUN$seed1)

ergm_call0 <- .make_translated_ergm_call(
  partition = partition_mid,
  rhs       = 'cov_match("bin_att")',
  nodes     = nodes_mid,
  dyads     = dyads_mid,
  verbose   = FALSE
)

fml0 <- ergm_call0[[2L]]
constraints0 <- .extract_constraints_from_ergm_call(ergm_call0)
if (is.null(constraints0)) stop("Translated ergm() call has no constraints; expected ~b1part.")

nw0 <- eval(fml0[[2L]], envir = environment(fml0))
.check_invariants_b1part(nw0, tag = "[initial]")

# --------------------------------------------------------------------------------------
# Flat probe model :
#   - offset(log_factorial_sizes(), 0) fixes coefficient at 0 => likelihood contribution is 0
#   - hence Δ loglik = 0 always; with symmetric proposal => accept ~1
# --------------------------------------------------------------------------------------
fml_flat <- nw0 ~ offset(log_factorial_sizes(), 0)

# ======================================================================================
# Probe runner 
# ======================================================================================
.run_probe_chain <- function(ctrl, nsteps, label) {
  cat("\n--- Phase 1 / PROBE:", label, "(real chain) ---\n")

  out <- simulate(
    fml_flat,
    nsim        = as.integer(nsteps + 1L),
    constraints = constraints0,
    control     = ctrl,
    verbose     = FALSE
  )

  nwl <- .as_network_list(out)
  if (length(nwl) != nsteps + 1L) {
    stop(sprintf("simulate() returned %d networks; expected %d.",
                 length(nwl), nsteps + 1L))
  }

  ct <- 0L; cs <- 0L; cstay <- 0L; co <- 0L

  for (k in seq_len(nsteps)) {
    nw_prev <- nwl[[k]]
    nw_next <- nwl[[k + 1L]]

    .check_invariants_b1part(nw_prev, tag = paste0("[", label, " k=", k, " prev]"))
    .check_invariants_b1part(nw_next, tag = paste0("[", label, " k=", k, " next]"))

    typ <- .classify_step(nw_prev, nw_next)
    if (typ == "toggle") ct <- ct + 1L
    else if (typ == "swap") cs <- cs + 1L
    else if (typ == "stay") cstay <- cstay + 1L
    else co <- co + 1L
  }

  denom <- ct + cs
  p_toggle <- if (denom > 0L) ct / denom else NA_real_

  cat(sprintf("Observed: toggle=%d | swap=%d | stay=%d | other=%d\n", ct, cs, cstay, co))
  cat(sprintf("Observed (among move!=stay): p(toggle)=%.4f\n", p_toggle))

  invisible(list(toggle = ct, swap = cs, stay = cstay, other = co, p_toggle = p_toggle))
}

.assert_mix <- function(res, label, p_target, tol) {
  if (res$other != 0L) {
    stop(sprintf("%s: unexpected 'other' transitions observed (%d).", label, res$other))
  }

  # stays should be ~0 on the flat model; warn (don’t fail) if not.
  if (res$stay != 0L) {
    warning(sprintf("%s: %d 'stay' transitions (rejections/no-op) observed on flat model.",
                    label, res$stay))
  }

  denom <- res$toggle + res$swap
  if (denom <= 0L || !is.finite(res$p_toggle)) {
    stop(sprintf("%s: no classified moves (toggle+swap=0).", label))
  }

  if (abs(res$p_toggle - p_target) > tol) {
    stop(sprintf(
      "%s: p(toggle)=%.4f, expected %.4f ± %.4f (toggle=%d swap=%d stay=%d).",
      label, res$p_toggle, p_target, tol, res$toggle, res$swap, res$stay
    ))
  }

  TRUE
}

# ======================================================================================
# Phase 1: PROBES
# ======================================================================================
cat("\n================================================================================\n")
cat("Phase 1) PROBES: ErpmMix decoding + custom proportions (REAL CHAIN)\n")
cat("================================================================================\n")

# --- TEST 1: default ErpmMix (no args) -> expect ~2/3 toggle, ~1/3 swap
cat("\n=== Phase 1 / TEST 1: default ErpmMix (no args) ===\n")
ctrl_default <- control.simulate.formula(
  MCMC.prop         = ~ .select("ErpmMix"),
  MCMC.prop.args    = list(list()),         # empty args => default packing/fallback
  MCMC.packagenames = "ERPM",
  MCMC.burnin       = 100,
  MCMC.interval     = 1
)

set.seed(RUN$seed1)
res1 <- .run_probe_chain(ctrl_default, RUN$probe_nsteps, "default-mix")
.assert_mix(res1, "Default mix", p_target = 2/3, tol = RUN$tol_default)
cat("Phase 1 / TEST 1 OK.\n")

# --- TEST 2: custom proportions (toggle:4, swap:1) -> expect ~0.8 toggle
cat("\n=== Phase 1 / TEST 2: custom mix (toggle:4, swap:1) ===\n")
ctrl_4_1 <- control.simulate.formula(
  MCMC.prop         = ~ .select("ErpmMix"),
  MCMC.prop.args    = list(list(moves = c("toggle", "swap"), weights = c(4, 1))),
  MCMC.packagenames = "ERPM",
  MCMC.burnin       = 100,
  MCMC.interval     = 1
)

set.seed(RUN$seed2)
res2 <- .run_probe_chain(ctrl_4_1, RUN$probe_nsteps, "mix-4-1")
.assert_mix(res2, "Mix 4:1", p_target = 4/5, tol = RUN$tol_4_1)
cat("Phase 1 / TEST 2 OK.\n")

# --- TEST 3: invalid args -> must fall back to canonical default mix (2/1)
cat("\n=== Phase 1 / TEST 3: invalid args -> fallback default ===\n")
ctrl_invalid <- control.simulate.formula(
  MCMC.prop         = ~ .select("ErpmMix"),
  MCMC.prop.args    = list(list(moves = c("toggle", "nope"), weights = c(1, 1))),
  MCMC.packagenames = "ERPM",
  MCMC.burnin       = 100,
  MCMC.interval     = 1
)

set.seed(RUN$seed3)
res3 <- .run_probe_chain(ctrl_invalid, RUN$probe_nsteps, "invalid->default")
.assert_mix(res3, "Invalid->default", p_target = 2/3, tol = RUN$tol_invalid)
cat("Phase 1 / TEST 3 OK.\n")

cat("\nPhase 1 OK: default + custom 4:1 + invalid fallback, using a real chain.\n")

# ======================================================================================
# Quick end-to-end fits 
# ======================================================================================
cat("\n================================================================================\n")
cat("PHASE 2) FITS : ErpmMix end-to-end through erpm() \n")
cat("================================================================================\n")

.ctrl_fit_mix_4_1 <- function() {
  ergm::control.ergm(
    MCMC.prop         = ~ .select("ErpmMix"),
    MCMC.prop.args    = list(list(moves = c("toggle", "swap"), weights = c(4, 1))),
    MCMC.packagenames = "ERPM",
    MCMC.burnin       = 20000,
    MCMC.interval     = 1,
    MCMC.samplesize   = 40000,
    seed              = RUN$seed_fit
  )
}

.run_fit <- function(rhs, label) {
  cat("\n--- PHASE 2 / FIT:", label, "---\n")
  fml <- as.formula(paste0("partition ~ ", rhs))
  environment(fml) <- list2env(list(partition = partition_mid), parent = parent.frame())

  fit <- erpm(
    fml,
    eval.call    = TRUE,
    verbose      = TRUE,
    nodes        = nodes_mid,
    dyads        = dyads_mid,
    constraints  = NULL,               # erpm default => ~b1part
    control      = .ctrl_fit_mix_4_1()
  )

  print(summary(fit))

  if (!inherits(fit, "ergm")) stop("Expected an 'ergm' fit object.")
  cat("OK: fitted.\n")
  invisible(fit)
}

fit1 <- .run_fit('log_factorial_sizes()', "log_factorial_sizes")
fit2 <- .run_fit('cov_match("bin_att")', "cov_match(bin_att)")
fit3 <- .run_fit('dyadcov("Z1")', "dyadcov(Z1)")

# (Optional) cliques: can be close to constant under many swaps; still ok as an end-to-end smoke test.
# fit4 <- .run_fit('cliques(k=2)', "cliques(k=2)")

cat("\nPHASE 2 OK: fits completed under ~b1part using ErpmMix (toggle:4 swap:1).\n")
cat("\nSELFTEST OK.\n")

if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()