# ==============================================================================
# File    : scripts/test/selftests/selftest_ergm_prop_erpm_mix_dynamic_toggle_swap_merge_split.R
# Purpose : Self-test (via erpm wrapper) for ERPM MH proposal ErpmMix under ~b1part,
#           validating dynamic decoding and practical activation of:
#             toggle / swap / merge / split
#
# What this selftest enforces:
#   1) "Real chain": each PROBE is a single MCMC simulation (simulate() called ONCE),
#      then transitions between successive states are classified.
#
#   2) Flat target:
#      - offset(log_factorial_sizes(), 0) => likelihood contribution always 0
#      - accepted moves are therefore driven mainly by the proposal kernel itself.
#
#   3) Classification:
#      - toggle : exactly 1 actor changes group
#      - swap   : exactly 2 actors change and memberships are swapped
#      - merge  : k>=1 actors moved from one group to another, and one group becomes empty
#      - split  : k>=1 actors moved from one group to a previously empty group
#      - stay   : identical successive states (rejection / impossible move / no-op)
#      - other  : anything else (MUST NOT occur)
#
#   4) Validation philosophy:
#      - for 2-way mixes (toggle/swap only), proportions are checked quantitatively;
#      - for 4-way mixes, the selftest is intentionally softer:
#          * all four moves must be decodable and observable,
#          * no "other" transition may occur,
#          * merge/split must appear when given positive weights,
#          * toggle must remain dominant over swap for the chosen weights,
#          * merge/split must remain minority moves on this dataset.
#
# Notes:
#   - merge/split feasibility depends on the current state of the chain.
#   - ErpmMix does not currently reweight or resample moves based on instantaneous feasibility.
#   - Therefore, in a 4-way mix, observed transition frequencies are NOT expected to match
#     nominal weights exactly.
#
# ==============================================================================

# --------------------------------------------------------------------------------------
# Init
# --------------------------------------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

options(ergm.loglik.warn_dyads = FALSE)
options(ERPM.zzz.verbose = TRUE)
options(Proposal.ErpmMix.debug = TRUE)

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

cat("=== SELFTEST ERPM: ErpmMix (dynamic args) under ~b1part (TOGGLE/SWAP/MERGE/SPLIT) ===\n")
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
  seed4 = 2468,
  seed5 = 8642,
  seed_fit = 1526,

  probe_nsteps_2way = 1600L,
  probe_nsteps_4way = 2200L,

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
    dA <- d[gA]; dB <- d[gB]

    if (!(dA + dB == 0L)) return("other")
    if (abs(dA) <= 0L) return("other")

    loser <- if (dA < 0L) gA else gB
    win   <- if (dA > 0L) gA else gB

    if (s1[loser] == 0L) return("merge")
    if (s0[win] == 0L)   return("split")
  }

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

.run_probe_chain <- function(ctrl, nsteps, label) {
  cat("\n--- Phase 1 / PROBE:", label, "(real chain) ---\n")

  out <- suppressWarnings(simulate(
    fml_flat,
    nsim        = as.integer(nsteps + 1L),
    constraints = constraints0,
    control     = ctrl,
    verbose     = FALSE
  ))

  nwl <- .as_network_list(out)
  if (length(nwl) != nsteps + 1L) {
    stop(sprintf("simulate() returned %d networks; expected %d.",
                 length(nwl), nsteps + 1L))
  }

  ct <- 0L; cs <- 0L; cm <- 0L; csp <- 0L; cstay <- 0L; co <- 0L

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
    stay = cstay, other = co, p = p
  ))
}

.assert_mix_2way <- function(res, label, p_toggle_target, tol) {
  if (res$other != 0L) {
    stop(sprintf("%s: unexpected 'other' transitions observed (%d).", label, res$other))
  }
  if (res$merge != 0L || res$split != 0L) {
    stop(sprintf("%s: merge/split observed but weights should be 0.", label))
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
    stop(sprintf("%s: merge was never observed although merge weight > 0.", label))
  }
  if (res_4way$split <= 0L) {
    stop(sprintf("%s: split was never observed although split weight > 0.", label))
  }

  if (res_4way$toggle <= res_4way$swap) {
    stop(sprintf("%s: toggle should remain more frequent than swap for the chosen weights.", label))
  }

  if ((res_4way$merge + res_4way$split) <= (res_base$merge + res_base$split)) {
    stop(sprintf("%s: merge/split activation did not increase versus the default mix.", label))
  }

  p4 <- res_4way$p[c("toggle", "swap", "merge", "split")]

  if (!(p4[["toggle"]] > p4[["swap"]])) {
    stop(sprintf("%s: expected p(toggle) > p(swap).", label))
  }

  if (!(p4[["merge"]] > 0 && p4[["split"]] > 0)) {
    stop(sprintf("%s: expected positive observed proportions for merge and split.", label))
  }

  if ((p4[["merge"]] + p4[["split"]]) >= 0.40) {
    stop(sprintf("%s: merge+split proportion unexpectedly high on this dataset (%.4f).",
                 label, p4[["merge"]] + p4[["split"]]))
  }

  TRUE
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

fml_flat <- nw0 ~ offset(log_factorial_sizes(), 0)

# ======================================================================================
# Phase 1: PROBES
# ======================================================================================
cat("\n================================================================================\n")
cat("Phase 1) PROBES: ErpmMix decoding + practical activation (REAL CHAIN)\n")
cat("================================================================================\n")

# --- TEST 1: default ErpmMix (no args) -> canonical toggle/swap only
cat("\n=== Phase 1 / TEST 1: default ErpmMix (no args) ===\n")
ctrl_default <- control.simulate.formula(
  MCMC.prop         = ~ .select("ErpmMix"),
  MCMC.prop.args    = list(list()),
  MCMC.packagenames = "ERPM",
  MCMC.burnin       = 100,
  MCMC.interval     = 1
)

set.seed(RUN$seed1)
res1 <- .run_probe_chain(ctrl_default, RUN$probe_nsteps_2way, "default-mix")
.assert_mix_2way(res1, "Default mix", p_toggle_target = 2/3, tol = RUN$tol_default)
cat("Phase 1 / TEST 1 OK.\n")

# --- TEST 2: custom 2-way (toggle:4, swap:1)
cat("\n=== Phase 1 / TEST 2: custom mix 2-way (toggle:4, swap:1) ===\n")
ctrl_4_1 <- control.simulate.formula(
  MCMC.prop         = ~ .select("ErpmMix"),
  MCMC.prop.args    = list(list(moves = c("toggle", "swap"), weights = c(4, 1))),
  MCMC.packagenames = "ERPM",
  MCMC.burnin       = 100,
  MCMC.interval     = 1
)

set.seed(RUN$seed2)
res2 <- .run_probe_chain(ctrl_4_1, RUN$probe_nsteps_2way, "mix-4-1")
.assert_mix_2way(res2, "Mix 4:1", p_toggle_target = 4/5, tol = RUN$tol_4_1)
cat("Phase 1 / TEST 2 OK.\n")

# --- TEST 3: invalid args -> fallback canonical default mix
cat("\n=== Phase 1 / TEST 3: invalid args -> fallback default ===\n")
ctrl_invalid <- control.simulate.formula(
  MCMC.prop         = ~ .select("ErpmMix"),
  MCMC.prop.args    = list(list(moves = c("toggle", "nope"), weights = c(1, 1))),
  MCMC.packagenames = "ERPM",
  MCMC.burnin       = 100,
  MCMC.interval     = 1
)

set.seed(RUN$seed3)
res3 <- .run_probe_chain(ctrl_invalid, RUN$probe_nsteps_2way, "invalid->default")
.assert_mix_2way(res3, "Invalid->default", p_toggle_target = 2/3, tol = RUN$tol_invalid)
cat("Phase 1 / TEST 3 OK.\n")

# --- TEST 4: custom 4-way mix -> soft validation only
cat("\n=== Phase 1 / TEST 4: custom 4-way (toggle/swap/merge/split) ===\n")
w4 <- c(toggle = 6, swap = 2, merge = 1, split = 1)

ctrl_4way <- control.simulate.formula(
  MCMC.prop         = ~ .select("ErpmMix"),
  MCMC.prop.args    = list(list(moves = names(w4), weights = as.numeric(w4))),
  MCMC.packagenames = "ERPM",
  MCMC.burnin       = 300,
  MCMC.interval     = 1
)

set.seed(RUN$seed4)
res4 <- .run_probe_chain(ctrl_4way, RUN$probe_nsteps_4way, "mix-4way")
.assert_mix_4way_soft(res1, res4, "Mix 4-way")
cat("Phase 1 / TEST 4 OK.\n")

# --- TEST 5: merge/split emphasized -> they must increase further
cat("\n=== Phase 1 / TEST 5: merge/split emphasized ===\n")
w5 <- c(toggle = 2, swap = 1, merge = 4, split = 4)

ctrl_ms <- control.simulate.formula(
  MCMC.prop         = ~ .select("ErpmMix"),
  MCMC.prop.args    = list(list(moves = names(w5), weights = as.numeric(w5))),
  MCMC.packagenames = "ERPM",
  MCMC.burnin       = 300,
  MCMC.interval     = 1
)

set.seed(RUN$seed5)
res5 <- .run_probe_chain(ctrl_ms, RUN$probe_nsteps_4way, "mix-merge-split-heavy")

if (res5$other != 0L) {
  stop(sprintf("Merge/split-heavy mix: unexpected 'other' transitions observed (%d).", res5$other))
}
if ((res5$merge + res5$split) <= (res4$merge + res4$split)) {
  stop("Merge/split-heavy mix: merge+split did not increase relative to the 4-way baseline.")
}
if (res5$merge <= 0L || res5$split <= 0L) {
  stop("Merge/split-heavy mix: merge and split should both be observed.")
}
cat("Phase 1 / TEST 5 OK.\n")

cat("\nPhase 1 OK: default + custom 2-way + invalid fallback + soft 4-way + merge/split-heavy.\n")

# ======================================================================================
# Phase 2: End-to-end fits (smoke tests)
# ======================================================================================
cat("\n================================================================================\n")
cat("PHASE 2) FITS : ErpmMix end-to-end through erpm() \n")
cat("================================================================================\n")

.ctrl_fit_mix_4way <- function() {
  ergm::control.ergm(
    MCMC.prop         = ~ .select("ErpmMix"),
    MCMC.prop.args    = list(list(moves = names(w4), weights = as.numeric(w4))),
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
    constraints  = NULL,
    control      = .ctrl_fit_mix_4way()
  )

  print(summary(fit))

  if (!inherits(fit, "ergm")) stop("Expected an 'ergm' fit object.")
  cat("OK: fitted.\n")
  invisible(fit)
}

fit1 <- .run_fit('log_factorial_sizes()', "log_factorial_sizes")
fit2 <- .run_fit('cov_match("bin_att")', "cov_match(bin_att)")
fit3 <- .run_fit('dyadcov("Z1")', "dyadcov(Z1)")

cat("\nPHASE 2 OK: fits completed under ~b1part using ErpmMix (4-way).\n")
cat("\nSELFTEST OK.\n")

if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()