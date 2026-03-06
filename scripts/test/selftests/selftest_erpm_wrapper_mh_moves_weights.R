# ==============================================================================
# File    : scripts/test/selftests/selftest_erpm_wrapper_mh_moves_weights.R
# Purpose : Self-test for ERPM wrapper arguments `mh_moves` / `mh_weights`
#           under ~b1part.
#
# File purpose
#   - PHASE 1 (PROBES) :
#       validate wrapper-level injection of ErpmMix through erpm(...), then
#       observe the requested move families on a real chain generated from the
#       translated ergm call.
#
#   - PHASE 2 (VALIDATION) :
#       check that erpm(...) rejects invalid mh_moves / mh_weights inputs.
#
#   - PHASE 3 (FITS) :
#       run end-to-end fits through erpm(...) with wrapper-level MH injection,
#       for three representative ERPM terms:
#         * log_factorial_sizes
#         * cov_match("bin_att")
#         * dyadcov("Z1")
#
# Notes
#   - This selftest targets the erpm(...) interface, not the low-level proposal
#     initializer directly.
#   - Phase 1 still uses simulate(), but only after erpm(eval.call=FALSE) has
#     built the translated ergm() call. What is tested is the wrapper wiring.
#   - In 4-way settings, observed move frequencies are not expected to match
#     nominal weights exactly because ErpmMix keeps a state-independent first
#     draw and falls back to TOGGLE when a selected move is infeasible.
#
# Run
#   Rscript scripts/test/selftests/selftest_erpm_wrapper_mh_moves_weights.R
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
  if (!requireNamespace("ergm", quietly = TRUE))    stop("Package 'ergm' required.")
})
suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm, quietly = TRUE, warn.conflicts = FALSE)
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

cat("=== SELFTEST ERPM WRAPPER: mh_moves / mh_weights under ~b1part ===\n")
cat("R:", paste(R.version$major, R.version$minor, sep = "."), "\n")
cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")
cat("ERPM loaded.\n")

# ======================================================================================
# Settings
# ======================================================================================
RUN <- list(
  phase1_probes            = TRUE,
  phase2_input_validation  = TRUE,
  phase3_erpm_fits         = TRUE,

  quiet = TRUE,

  seed1 = 123,
  seed2 = 456,
  seed3 = 2468,
  seed4 = 8642,
  seed_fit = 1526,

  probe_nsteps_2way = 1600L,
  probe_nsteps_4way = 2200L,

  tol_default = 0.06,   # expected 2/3
  tol_4_1     = 0.05,   # expected 4/5

  fit_burnin     = 20000,
  fit_interval   = 1,
  fit_samplesize = 40000
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
`%||%` <- function(a, b) if (!is.null(a)) a else b

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
  if (length(idx) == 0L) return("stay")
  if (length(idx) == 1L) return("toggle")

  if (length(idx) == 2L) {
    i <- idx[1L]; j <- idx[2L]
    if (g1[i] == g0[j] && g1[j] == g0[i]) return("swap")
  }

  tab0 <- table(g0)
  tab1 <- table(g1)

  groups <- union(names(tab0), names(tab1))
  s0 <- setNames(integer(length(groups)), groups)
  s1 <- setNames(integer(length(groups)), groups)
  s0[names(tab0)] <- as.integer(tab0)
  s1[names(tab1)] <- as.integer(tab1)

  d <- s1 - s0
  nz <- which(d != 0L)

  if (length(nz) == 2L) {
    gA <- names(d)[nz[1L]]
    gB <- names(d)[nz[2L]]
    dA <- d[gA]
    dB <- d[gB]

    if (!(dA + dB == 0L)) return("other")
    if (abs(dA) <= 0L) return("other")

    loser <- if (dA < 0L) gA else gB
    win   <- if (dA > 0L) gA else gB

    if (s1[loser] == 0L) return("merge")
    if (s0[win] == 0L)   return("split")
  }

  "other"
}

.extract_constraints_from_ergm_call <- function(ergm_call) {
  al <- as.list(ergm_call)
  nm <- names(al)
  i <- which(nm == "constraints")
  if (length(i) == 1L) al[[i]] else NULL
}

.extract_control_from_ergm_call <- function(ergm_call) {
  al <- as.list(ergm_call)
  nm <- names(al)
  i <- which(nm == "control")
  if (length(i) == 1L) al[[i]] else NULL
}

.assert_erpmmix_injection <- function(ergm_call, expected_moves, expected_weights, label) {
  ctrl <- .extract_control_from_ergm_call(ergm_call)
  if (is.null(ctrl)) {
    stop(sprintf("%s: translated ergm() call has no control argument.", label))
  }

  prop_expr <- ctrl$MCMC.prop
  if (is.null(prop_expr)) {
    stop(sprintf("%s: control object has no MCMC.prop.", label))
  }

  prop_str <- paste(deparse(prop_expr, width.cutoff = 500L), collapse = " ")
  prop_str <- gsub("\\s+", "", prop_str)
  if (!identical(prop_str, "~.select(\"ErpmMix\")")) {
    stop(sprintf("%s: expected MCMC.prop = ~ .select(\"ErpmMix\"), got: %s", label, prop_str))
  }

  prop_args <- ctrl$MCMC.prop.args
  if (!is.list(prop_args) || length(prop_args) != 1L || !is.list(prop_args[[1L]])) {
    stop(sprintf("%s: invalid MCMC.prop.args structure.", label))
  }

  got_moves <- prop_args[[1L]]$moves
  got_weights <- prop_args[[1L]]$weights

  if (!identical(as.character(got_moves), as.character(expected_moves))) {
    stop(sprintf(
      "%s: injected moves mismatch.\n  expected: %s\n  got     : %s",
      label,
      paste(expected_moves, collapse = ", "),
      paste(got_moves, collapse = ", ")
    ))
  }

  if (!isTRUE(all.equal(as.numeric(got_weights), as.numeric(expected_weights)))) {
    stop(sprintf(
      "%s: injected weights mismatch.\n  expected: %s\n  got     : %s",
      label,
      paste(expected_weights, collapse = ", "),
      paste(got_weights, collapse = ", ")
    ))
  }

  TRUE
}

.make_erpm_call_with_wrapper_mix <- function(partition,
                                             rhs,
                                             nodes,
                                             dyads,
                                             mh_moves,
                                             mh_weights,
                                             verbose = FALSE) {
  if (!exists("erpm", mode = "function")) stop("erpm() not found (ERPM not loaded?).")

  user_formula <- as.formula(paste0("partition ~ ", rhs))
  environment(user_formula) <- list2env(list(partition = partition), parent = parent.frame())

  call <- erpm(
    user_formula,
    eval.call   = FALSE,
    verbose     = verbose,
    nodes       = nodes,
    dyads       = dyads,
    constraints = NULL,
    mh_moves    = mh_moves,
    mh_weights  = mh_weights
  )

  if (!is.call(call) || !identical(call[[1L]], as.name("ergm"))) {
    stop("erpm(eval.call=FALSE) did not return an ergm() call.")
  }

  call
}

.make_sim_control_from_ergm_call <- function(ergm_call,
                                             burnin = 100L,
                                             interval = 1L) {
  ctrl_ergm <- .extract_control_from_ergm_call(ergm_call)
  if (is.null(ctrl_ergm)) {
    stop("Translated ergm() call has no control argument.")
  }

  ergm::control.simulate.formula(
    MCMC.prop         = ctrl_ergm$MCMC.prop,
    MCMC.prop.args    = ctrl_ergm$MCMC.prop.args,
    MCMC.packagenames = ctrl_ergm$MCMC.packagenames %||% "ERPM",
    MCMC.burnin       = as.integer(burnin),
    MCMC.interval     = as.integer(interval)
  )
}

.run_probe_chain_from_erpm <- function(mh_moves,
                                       mh_weights,
                                       nsteps,
                                       label,
                                       burnin = 100L) {
  cat("\n--- Phase 1 / PROBE:", label, "(real chain via erpm wrapper) ---\n")

  ergm_call <- .make_erpm_call_with_wrapper_mix(
    partition  = partition_mid,
    rhs        = 'offset(log_factorial_sizes(), 0)',
    nodes      = nodes_mid,
    dyads      = dyads_mid,
    mh_moves   = mh_moves,
    mh_weights = mh_weights,
    verbose    = FALSE
  )

  .assert_erpmmix_injection(
    ergm_call        = ergm_call,
    expected_moves   = mh_moves,
    expected_weights = mh_weights,
    label            = label
  )

  fml <- ergm_call[[2L]]
  constraints0 <- .extract_constraints_from_ergm_call(ergm_call)
  if (is.null(constraints0)) {
    stop(sprintf("%s: translated ergm() call has no constraints.", label))
  }

  ctrl_sim <- .make_sim_control_from_ergm_call(
    ergm_call = ergm_call,
    burnin    = burnin,
    interval  = 1L
  )

  out <- suppressWarnings(simulate(
    fml,
    nsim        = as.integer(nsteps + 1L),
    constraints = constraints0,
    control     = ctrl_sim,
    verbose     = FALSE
  ))

  nwl <- .as_network_list(out)
  if (length(nwl) != nsteps + 1L) {
    stop(sprintf(
      "%s: simulate() returned %d networks; expected %d.",
      label, length(nwl), nsteps + 1L
    ))
  }

  ct <- 0L
  cs <- 0L
  cm <- 0L
  csp <- 0L
  cstay <- 0L
  co <- 0L

  for (k in seq_len(nsteps)) {
    nw_prev <- nwl[[k]]
    nw_next <- nwl[[k + 1L]]

    .check_invariants_b1part(nw_prev, tag = paste0("[", label, " k=", k, " prev]"))
    .check_invariants_b1part(nw_next, tag = paste0("[", label, " k=", k, " next]"))

    typ <- .classify_step(nw_prev, nw_next)
    if (typ == "toggle") ct <- ct + 1L
    else if (typ == "swap") cs <- cs + 1L
    else if (typ == "merge") cm <- cm + 1L
    else if (typ == "split") csp <- csp + 1L
    else if (typ == "stay") cstay <- cstay + 1L
    else co <- co + 1L
  }

  denom <- ct + cs + cm + csp
  p <- if (denom > 0L) {
    c(toggle = ct, swap = cs, merge = cm, split = csp) / denom
  } else {
    c(toggle = NA_real_, swap = NA_real_, merge = NA_real_, split = NA_real_)
  }

  cat(sprintf(
    "Observed: toggle=%d | swap=%d | merge=%d | split=%d | stay=%d | other=%d\n",
    ct, cs, cm, csp, cstay, co
  ))
  cat(sprintf(
    "Observed (among classified non-stay): p(toggle)=%.4f p(swap)=%.4f p(merge)=%.4f p(split)=%.4f\n",
    p[["toggle"]], p[["swap"]], p[["merge"]], p[["split"]]
  ))

  invisible(list(
    toggle = ct, swap = cs, merge = cm, split = csp,
    stay = cstay, other = co, p = p,
    ergm_call = ergm_call
  ))
}

.assert_mix_2way <- function(res, label, p_toggle_target, tol) {
  if (res$other != 0L) {
    stop(sprintf("%s: unexpected 'other' transitions observed (%d).", label, res$other))
  }
  if (res$merge != 0L || res$split != 0L) {
    stop(sprintf("%s: merge/split observed but they were not requested.", label))
  }

  denom <- res$toggle + res$swap
  if (denom <= 0L) {
    stop(sprintf("%s: no classified moves (toggle+swap=0).", label))
  }

  p_toggle <- res$toggle / denom
  if (abs(p_toggle - p_toggle_target) > tol) {
    stop(sprintf(
      "%s: p(toggle)=%.4f, expected %.4f ± %.4f (toggle=%d swap=%d stay=%d other=%d).",
      label, p_toggle, p_toggle_target, tol, res$toggle, res$swap, res$stay, res$other
    ))
  }

  TRUE
}

.assert_mix_4way_soft <- function(res_base, res_4way, label) {
  if (res_4way$other != 0L) {
    stop(sprintf("%s: unexpected 'other' transitions observed (%d).", label, res_4way$other))
  }

  denom4 <- res_4way$toggle + res_4way$swap + res_4way$merge + res_4way$split
  if (denom4 <= 0L) {
    stop(sprintf("%s: no classified moves in 4-way mix.", label))
  }

  if (res_4way$merge <= 0L) {
    stop(sprintf("%s: merge was never observed although requested with positive weight.", label))
  }
  if (res_4way$split <= 0L) {
    stop(sprintf("%s: split was never observed although requested with positive weight.", label))
  }

  if (res_4way$toggle <= res_4way$swap) {
    stop(sprintf("%s: toggle should remain more frequent than swap for the chosen weights.", label))
  }

  if ((res_4way$merge + res_4way$split) <= (res_base$merge + res_base$split)) {
    stop(sprintf("%s: merge/split activation did not increase versus the 2-way baseline.", label))
  }

  p4 <- res_4way$p[c("toggle", "swap", "merge", "split")]

  if (!(p4[["toggle"]] > p4[["swap"]])) {
    stop(sprintf("%s: expected p(toggle) > p(swap).", label))
  }

  if (!(p4[["merge"]] > 0 && p4[["split"]] > 0)) {
    stop(sprintf("%s: expected positive observed proportions for merge and split.", label))
  }

  if ((p4[["merge"]] + p4[["split"]]) >= 0.40) {
    stop(sprintf(
      "%s: merge+split proportion unexpectedly high on this dataset (%.4f).",
      label, p4[["merge"]] + p4[["split"]]
    ))
  }

  TRUE
}

.assert_erpm_error <- function(expr, pattern, label) {
  got <- try(eval.parent(substitute(expr)), silent = TRUE)
  if (!inherits(got, "try-error")) {
    stop(sprintf("%s: expected an error, but the call succeeded.", label))
  }
  msg <- conditionMessage(attr(got, "condition"))
  if (!grepl(pattern, msg, fixed = TRUE)) {
    stop(sprintf(
      "%s: error message mismatch.\nExpected pattern: %s\nGot: %s",
      label, pattern, msg
    ))
  }
  TRUE
}

.ctrl_fit_wrapper_mix_4way <- function() {
  ergm::control.ergm(
    MCMC.burnin     = RUN$fit_burnin,
    MCMC.interval   = RUN$fit_interval,
    MCMC.samplesize = RUN$fit_samplesize,
    seed            = RUN$seed_fit
  )
}

.run_fit_via_erpm <- function(rhs, label, mh_moves, mh_weights) {
  cat("\n--- Phase 3 / FIT:", label, "---\n")

  fml <- as.formula(paste0("partition_mid ~ ", rhs))
  environment(fml) <- list2env(list(partition_mid = partition_mid), parent = parent.frame())

  fit <- erpm(
    fml,
    eval.call   = TRUE,
    verbose     = !RUN$quiet,
    nodes       = nodes_mid,
    dyads       = dyads_mid,
    constraints = NULL,
    mh_moves    = mh_moves,
    mh_weights  = mh_weights,
    control     = .ctrl_fit_wrapper_mix_4way()
  )

  if (!inherits(fit, "ergm")) {
    stop(sprintf("%s: expected an 'ergm' fit object.", label))
  }

  print(summary(fit))
  cat("OK: fitted.\n")
  invisible(fit)
}

# ======================================================================================
# Sanity check on initial translated network
# ======================================================================================
set.seed(RUN$seed1)

ergm_call0 <- .make_erpm_call_with_wrapper_mix(
  partition  = partition_mid,
  rhs        = 'cov_match("bin_att")',
  nodes      = nodes_mid,
  dyads      = dyads_mid,
  mh_moves   = c("toggle", "swap"),
  mh_weights = c(2, 1),
  verbose    = FALSE
)

fml0 <- ergm_call0[[2L]]
nw0 <- eval(fml0[[2L]], envir = environment(fml0))
.check_invariants_b1part(nw0, tag = "[initial]")

# ======================================================================================
# Phase 1: PROBES through erpm(...)
# ======================================================================================
run_phase1_probes <- function() {
  cat("\n================================================================================\n")
  cat("Phase 1) PROBES: erpm(...) wrapper injection + practical activation\n")
  cat("================================================================================\n")

  cat("\n=== Phase 1 / TEST 1: erpm(..., mh_moves=c('toggle','swap'), mh_weights=c(2,1)) ===\n")
  set.seed(RUN$seed1)
  res1 <- .run_probe_chain_from_erpm(
    mh_moves   = c("toggle", "swap"),
    mh_weights = c(2, 1),
    nsteps     = RUN$probe_nsteps_2way,
    label      = "wrapper-mix-2-1"
  )
  .assert_mix_2way(res1, "Wrapper mix 2:1", p_toggle_target = 2/3, tol = RUN$tol_default)
  cat("Phase 1 / TEST 1 OK.\n")

  cat("\n=== Phase 1 / TEST 2: erpm(..., mh_moves=c('toggle','swap'), mh_weights=c(4,1)) ===\n")
  set.seed(RUN$seed2)
  res2 <- .run_probe_chain_from_erpm(
    mh_moves   = c("toggle", "swap"),
    mh_weights = c(4, 1),
    nsteps     = RUN$probe_nsteps_2way,
    label      = "wrapper-mix-4-1"
  )
  .assert_mix_2way(res2, "Wrapper mix 4:1", p_toggle_target = 4/5, tol = RUN$tol_4_1)
  cat("Phase 1 / TEST 2 OK.\n")

  cat("\n=== Phase 1 / TEST 3: erpm(..., mh_moves=toggle/swap/merge/split, mh_weights=6/2/1/1) ===\n")
  w4 <- c(toggle = 6, swap = 2, merge = 1, split = 1)

  set.seed(RUN$seed3)
  res3 <- .run_probe_chain_from_erpm(
    mh_moves   = names(w4),
    mh_weights = as.numeric(w4),
    nsteps     = RUN$probe_nsteps_4way,
    label      = "wrapper-mix-4way"
  )
  .assert_mix_4way_soft(res1, res3, "Wrapper mix 4-way")
  cat("Phase 1 / TEST 3 OK.\n")

  cat("\n=== Phase 1 / TEST 4: erpm(..., mh_moves=toggle/swap/merge/split, mh_weights=2/1/4/4) ===\n")
  w5 <- c(toggle = 2, swap = 1, merge = 4, split = 4)

  set.seed(RUN$seed4)
  res4 <- .run_probe_chain_from_erpm(
    mh_moves   = names(w5),
    mh_weights = as.numeric(w5),
    nsteps     = RUN$probe_nsteps_4way,
    label      = "wrapper-mix-merge-split-heavy"
  )

  if (res4$other != 0L) {
    stop(sprintf("Wrapper merge/split-heavy mix: unexpected 'other' transitions observed (%d).", res4$other))
  }
  if ((res4$merge + res4$split) <= (res3$merge + res3$split)) {
    stop("Wrapper merge/split-heavy mix: merge+split did not increase relative to the 4-way baseline.")
  }
  if (res4$merge <= 0L || res4$split <= 0L) {
    stop("Wrapper merge/split-heavy mix: merge and split should both be observed.")
  }
  cat("Phase 1 / TEST 4 OK.\n")

  cat("\nPhase 1 OK: wrapper injection + 2-way proportions + 4-way activation.\n")
  invisible(list(test1 = res1, test2 = res2, test3 = res3, test4 = res4))
}

# ======================================================================================
# Phase 2: wrapper-level validation errors
# ======================================================================================
run_phase2_input_validation <- function() {
  cat("\n================================================================================\n")
  cat("Phase 2) INPUT VALIDATION through erpm(...)\n")
  cat("================================================================================\n")

  cat("\n--- Phase 2 / TEST 1: mh_moves without mh_weights -> must fail ---\n")
  .assert_erpm_error(
    erpm(
      partition_mid ~ cov_match("bin_att"),
      eval.call   = FALSE,
      verbose     = FALSE,
      nodes       = nodes_mid,
      dyads       = dyads_mid,
      mh_moves    = c("toggle", "swap"),
      mh_weights  = NULL
    ),
    "`mh_moves` and `mh_weights` must be provided together or both left NULL.",
    "mh_moves only"
  )
  cat("Phase 2 / TEST 1 OK.\n")

  cat("\n--- Phase 2 / TEST 2: unsupported move name -> must fail ---\n")
  .assert_erpm_error(
    erpm(
      partition_mid ~ cov_match("bin_att"),
      eval.call   = FALSE,
      verbose     = FALSE,
      nodes       = nodes_mid,
      dyads       = dyads_mid,
      mh_moves    = c("toggle", "nope"),
      mh_weights  = c(1, 1)
    ),
    "`mh_moves` contains unsupported move(s): nope.",
    "unsupported move"
  )
  cat("Phase 2 / TEST 2 OK.\n")

  cat("\n--- Phase 2 / TEST 3: mismatched lengths -> must fail ---\n")
  .assert_erpm_error(
    erpm(
      partition_mid ~ cov_match("bin_att"),
      eval.call   = FALSE,
      verbose     = FALSE,
      nodes       = nodes_mid,
      dyads       = dyads_mid,
      mh_moves    = c("toggle", "swap"),
      mh_weights  = c(1, 1, 1)
    ),
    "`mh_moves` and `mh_weights` must have the same length.",
    "mismatched lengths"
  )
  cat("Phase 2 / TEST 3 OK.\n")

  cat("\nPhase 2 OK: wrapper input validation behaves as expected.\n")
  invisible(TRUE)
}

# ======================================================================================
# Phase 3: fits via erpm(...)
# ======================================================================================
run_phase3_erpm_fits <- function() {
  cat("\n================================================================================\n")
  cat("Phase 3) FITS through erpm(...)\n")
  cat("================================================================================\n")

  wfit <- c(toggle = 6, swap = 2, merge = 1, split = 1)

  fit1 <- .run_fit_via_erpm(
    rhs        = "log_factorial_sizes()",
    label      = "log_factorial_sizes",
    mh_moves   = names(wfit),
    mh_weights = as.numeric(wfit)
  )

  fit2 <- .run_fit_via_erpm(
    rhs        = 'cov_match("bin_att")',
    label      = 'cov_match("bin_att")',
    mh_moves   = names(wfit),
    mh_weights = as.numeric(wfit)
  )

  fit3 <- .run_fit_via_erpm(
    rhs        = 'dyadcov("Z1")',
    label      = 'dyadcov("Z1")',
    mh_moves   = names(wfit),
    mh_weights = as.numeric(wfit)
  )

  cat("\nPhase 3 OK: fits completed under erpm(...) with wrapper-injected ErpmMix.\n")
  invisible(list(
    log_factorial_sizes = fit1,
    cov_match = fit2,
    dyadcov = fit3
  ))
}

# ======================================================================================
# Main run
# ======================================================================================
run_all_selftests_erpm_wrapper_mh_moves_weights <- function() {
  if (isTRUE(RUN$phase1_probes)) run_phase1_probes()
  else cat("\n================================================================================\nPhase 1) SKIP\n================================================================================\n")

  if (isTRUE(RUN$phase2_input_validation)) run_phase2_input_validation()
  else cat("\n================================================================================\nPhase 2) SKIP\n================================================================================\n")

  if (isTRUE(RUN$phase3_erpm_fits)) run_phase3_erpm_fits()
  else cat("\n================================================================================\nPhase 3) SKIP\n================================================================================\n")

  cat("\nSELFTEST OK.\n")
  invisible(TRUE)
}

if (identical(environment(), globalenv())) {
  run_all_selftests_erpm_wrapper_mh_moves_weights()
}

if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()