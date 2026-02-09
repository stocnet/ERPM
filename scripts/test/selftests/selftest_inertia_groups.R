# ======================================================================================
# File    : scripts/test/selftests/selftest_inertia_groups_PLE.R
# Object  : Self-test (PLE only) for ERPM inertial term `inertia_groups`
# Run     : Rscript scripts/test/selftests/selftest_inertia_groups_PLE.R
# Notes   :
#   - Only PLE ("empile") is used in this selftest.
#   - No random generation: partitions, nodes, dyads are explicit numeric values.
#   - When calling erpm_long(): always provide nodes AND dyads, even if unused.
#   - When calling erpm_long(): verbose=TRUE and debug="deep".
#   - After each erpm_long(): print(out) (or print with a short label).
#   - Static effects exercised: cliques, cov_fullmatch, dyadcov
#   - Inertial effect: inertia_groups
#   - This version does NOT rely on blockdiag() (known unimplemented combo with b1part in your setup).
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule
# --------------------------------------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("network",   quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",      quietly = TRUE)) stop("Package 'ergm' requis.")
  if (!requireNamespace("devtools",  quietly = TRUE)) stop("Package 'devtools' requis.")
  if (!requireNamespace("rprojroot", quietly = TRUE)) stop("Package 'rprojroot' requis.")
})

suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

devtools::load_all(".")

# --------------------------------------------------------------------------------------
# Warnings capture helper (selftest)
# --------------------------------------------------------------------------------------
.warn_flush <- function() {
  invisible(warnings())
  invisible(NULL)
}

.capture_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = w)
}

# --------------------------------------------------------------------------------------
# Patch ERGM (optionnel, ne doit pas faire échouer le selftest)
# --------------------------------------------------------------------------------------
root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)

.patch_enabled <- FALSE
patch_path <- file.path(root, "scripts", "ergm_patch.R")
if (!file.exists(patch_path)) {
  message("[selftest_inertia_groups_PLE] scripts/ergm_patch.R introuvable: selftest continue sans patch.")
} else {
  source(patch_path, local = FALSE)
  if (exists("ergm_patch_enable", mode = "function")) {
    tryCatch(
      {
        ergm_patch_enable()
        .patch_enabled <- TRUE
      },
      error = function(e) {
        message("[selftest_inertia_groups_PLE] ergm_patch_enable() failed, continuing without patch.\n",
                "  reason: ", conditionMessage(e))
        .patch_enabled <- FALSE
      }
    )
  } else {
    message("[selftest_inertia_groups_PLE] ergm_patch_enable() introuvable: selftest continue sans patch.")
  }
}

cat("=== SELFTEST inertia_groups (PLE only) | NO blockdiag ===\n")

# ======================================================================================
# Global fixed scenario (explicit values)
# ======================================================================================

# Stacked (PLE) dimensions
B       <- 3L   # number of blocks (times)
n_block <- 4L   # actors per block
G_block <- 3L   # groups per block

# Total actors in stacked actor space (mode-1)
n1_total <- B * n_block

# Total vertices in stacked bipartite network
N_total <- n1_total + (B * G_block)

# Current partitions per block (time blocks inside the stacked meta-network)
P_cur_b1 <- c(1, 1, 2, 3)
P_cur_b2 <- c(1, 2, 2, 3)
P_cur_b3 <- c(1, 1, 3, 3)

# Past partitions per block by lag (lag=1 recent, lag=2 older)
P_b1_lag1 <- c(1, 1, 2, 3)
P_b1_lag2 <- c(1, 2, 2, 3)

P_b2_lag1 <- c(1, 2, 2, 3)
P_b2_lag2 <- c(1, 1, 2, 3)

P_b3_lag1 <- c(1, 1, 3, 3)
P_b3_lag2 <- c(1, 2, 2, 3)

# Container expected by InitErgmTerm.inertia_groups (PLE):
erpm_block_past_partitions <- list(
  list(P_b1_lag1, P_b1_lag2),
  list(P_b2_lag1, P_b2_lag2),
  list(P_b3_lag1, P_b3_lag2)
)

# ======================================================================================
# Shared helpers: stacking, oracle stats, and "no blockdiag" workaround for fits
# ======================================================================================

# --- Build stacked dyads (actor-actor) as block-diagonal n1_total x n1_total matrix
# Input: dyads_by_time = list(t=1..B) of list(fm=4x4, Z1=4x4)
.stack_actor_dyads_blockdiag <- function(dyads_by_time, name) {
  stopifnot(length(dyads_by_time) == B)
  M <- matrix(0, nrow = n1_total, ncol = n1_total)
  for (b in seq_len(B)) {
    Z <- dyads_by_time[[b]][[name]]
    if (!is.matrix(Z) || nrow(Z) != n_block || ncol(Z) != n_block) {
      stop(sprintf("dyads[%d]$%s must be a %dx%d matrix", b, name, n_block, n_block), call. = FALSE)
    }
    i0 <- (b - 1L) * n_block
    idx <- (i0 + 1L):(i0 + n_block)
    M[idx, idx] <- Z
  }
  M
}

.vec <- function(M) as.numeric(M)

# --- Oracle: current membership edges => per group member actor indices (global 1..n1_total)
# Groups are vertices (n1_total+1)..N_total, blockwise (G_block per block)
.get_group_members_global <- function(nw, gv) {
  nb <- network::get.neighborhood(nw, gv, type = "all")
  cur_ids <- sort(unique(as.integer(nb)))
  cur_ids <- cur_ids[cur_ids >= 1L & cur_ids <= n1_total]
  cur_ids
}

# --- Oracle inertia_groups (exogenous, PLE): strict intersection over lags 1..d
oracle_inertia_groups <- function(nw, d, size_filter_int) {
  cnt <- 0L
  for (b in seq_len(B)) {
    for (g in seq_len(G_block)) {
      gv <- n1_total + (b - 1L) * G_block + g
      cur_ids <- .get_group_members_global(nw, gv)
      if (length(cur_ids) == 0L) next
      if (length(size_filter_int) > 0L && !(length(cur_ids) %in% size_filter_int)) next

      ok_all_lags <- TRUE
      for (lag in seq_len(d)) {
        p_lag <- erpm_block_past_partitions[[b]][[lag]]
        grp_split <- split(seq_along(p_lag), as.integer(p_lag))
        grp_global <- lapply(grp_split, function(v) sort((b - 1L) * n_block + as.integer(v)))

        ok_lag <- FALSE
        for (u in seq_along(grp_global)) {
          if (length(grp_global[[u]]) == length(cur_ids) && all(grp_global[[u]] == cur_ids)) {
            ok_lag <- TRUE
            break
          }
        }
        if (!ok_lag) {
          ok_all_lags <- FALSE
          break
        }
      }
      if (ok_all_lags) cnt <- cnt + 1L
    }
  }
  cnt
}

# --- Oracle cliques(k=2): sum over groups choose(size,2)
oracle_cliques_k2 <- function(nw) {
  tot <- 0
  for (b in seq_len(B)) {
    for (g in seq_len(G_block)) {
      gv <- n1_total + (b - 1L) * G_block + g
      m <- length(.get_group_members_global(nw, gv))
      if (m >= 2L) tot <- tot + (m * (m - 1L)) / 2
    }
  }
  tot
}

# --- Oracle cov-based within-group sums: sum_{groups} sum_{i<j in group} M[i,j]
oracle_within_group_upper_sum <- function(nw, M_actor) {
  stopifnot(is.matrix(M_actor), nrow(M_actor) == n1_total, ncol(M_actor) == n1_total)
  tot <- 0
  for (b in seq_len(B)) {
    for (g in seq_len(G_block)) {
      gv <- n1_total + (b - 1L) * G_block + g
      ids <- .get_group_members_global(nw, gv)
      if (length(ids) < 2L) next
      ids <- as.integer(ids)
      for (ii in 1L:(length(ids) - 1L)) {
        for (jj in (ii + 1L):length(ids)) {
          tot <- tot + M_actor[ids[ii], ids[jj]]
        }
      }
    }
  }
  tot
}

# --- NO blockdiag workaround for fits (keep b1part only):
# Add a structural "almost-forbidden" edgecov on cross-block actor-group dyads.
# This DOES NOT use blockdiag(). It makes cross-block membership edges essentially impossible.
.build_forbidden_crossblock_edgecov <- function() {
  X <- matrix(0, nrow = N_total, ncol = N_total)

  # actor i in block ba, group vertex in block bg: forbid if ba != bg
  # Actors: 1..n1_total (global)
  # Groups: (n1_total+1)..N_total, blockwise
  for (ba in seq_len(B)) {
    a0 <- (ba - 1L) * n_block
    a_idx <- (a0 + 1L):(a0 + n_block)

    for (bg in seq_len(B)) {
      g0 <- n1_total + (bg - 1L) * G_block
      g_idx <- (g0 + 1L):(g0 + G_block)

      if (ba != bg) {
        # Put a huge negative value so offset(edgecov(X)) kills those edges.
        X[a_idx, g_idx] <- -1e6
        X[g_idx, a_idx] <- -1e6
      }
    }
  }
  X
}

# ======================================================================================
# SECTION 0) Build explicit inputs for Section 1 (manual meta-network + stacked dyads/monads)
# ======================================================================================

# Nodes (explicit, per time)
nodes <- list(
  data.frame(
    label  = c("A", "B", "C", "D"),
    gender = c(1, 1, 2, 1),
    age    = c(20, 22, 25, 30)
  ),
  data.frame(
    label  = c("FT", "AZ", "JI", "DO"),
    gender = c(2, 1, 2, 2),
    age    = c(10, 42, 25, 30)
  ),
  data.frame(
    label  = c("H", "Z", "S", "A"),
    gender = c(1, 1, 1, 1),
    age    = c(27, 26, 25, 28)
  )
)

# Dyads (explicit, per time)
dyads <- list(
  list(
    fm = matrix(c(
      0, 1, 1, 0,
      1, 0, 1, 0,
      1, 1, 0, 0,
      0, 0, 0, 0
    ), nrow = 4, byrow = TRUE),
    Z1 = matrix(c(
      0, 2, 3, 0,
      2, 0, 4, 0,
      3, 4, 0, 1,
      0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    fm = matrix(c(
      0, 1, 0, 1,
      1, 0, 1, 0,
      0, 1, 0, 1,
      1, 0, 1, 0
    ), nrow = 4, byrow = TRUE),
    Z1 = matrix(c(
      0, 5, 0, 2,
      5, 0, 1, 0,
      0, 1, 0, 3,
      2, 0, 3, 0
    ), nrow = 4, byrow = TRUE)
  ),
  list(
    fm = matrix(c(
      0, 0, 1, 1,
      0, 0, 1, 1,
      1, 1, 0, 0,
      1, 1, 0, 0
    ), nrow = 4, byrow = TRUE),
    Z1 = matrix(c(
      0, 1, 2, 3,
      1, 0, 4, 5,
      2, 4, 0, 6,
      3, 5, 6, 0
    ), nrow = 4, byrow = TRUE)
  )
)

# Partitions passed to erpm_long (explicit)
partitions <- list(
  c(1, 1, 2, 3),
  c(1, 2, 2, 3),
  c(1, 1, 3, 3)
)

# Stacked actor dyads (n1_total x n1_total) for manual meta-network
fm_stacked <- .stack_actor_dyads_blockdiag(dyads, "fm")
Z1_stacked <- .stack_actor_dyads_blockdiag(dyads, "Z1")


# ======================================================================================
# SECTION 1) summary() on a meta-network built by .erpm_long_empile_run() (ALL effects)
# Objective:
#   - build the stacked PLE meta-network using the exposed engine runner
#     .erpm_long_empile_run(...), i.e. exactly the internal path used by erpm_long();
#   - run summary() on: inertia_groups + cliques + cov_fullmatch + dyadcov
#   - verify each effect vs offline oracle (same nw)
# ======================================================================================

cat("\n================================================================================\n")
cat("SECTION 1) summary() on meta-network built by .erpm_long_empile_build_meta_nw() (ALL effects)\n")
cat("================================================================================\n")

# --------------------------------------------------------------------------------------
# 1) Build meta-network (PLE) via the exposed engine runner
# --------------------------------------------------------------------------------------

.f_empile_run <- get0(".erpm_long_empile_build_meta_nw", mode = "function", inherits = TRUE)
if (is.null(.f_empile_run) && "ERPM" %in% loadedNamespaces()) {
  .f_empile_run <- get0(".erpm_long_empile_build_meta_nw", envir = asNamespace("ERPM"),
                        mode = "function", inherits = FALSE)
}
if (is.null(.f_empile_run)) {
  stop("[SECTION 1] Cannot find exposed .erpm_long_empile_build_meta_nw() in current session.", call. = FALSE)
}

# RHS must be a pure ERGM RHS (already translated if needed). We only need it because the
# engine returns it; the build itself depends on inertial_present/past_influence.
rhs_all <- ~ inertia_groups(past_influence = 1) +
            cliques(k = 2) +
            cov_fullmatch("fm") +
            dyadcov("Z1")

# We want a "full" meta-network that supports checks for d=1 and d=2 on the SAME network,
# so we build with inertial_present=FALSE (=> idx_est = 1..T) and do not drop early blocks.
# This is purely for the summary/oracle validation; the erpm_long(estimation) behavior
# with inertial terms (dropping first d blocks) is tested in SECTION 2.
print(partitions)
print(rhs_all)
print(nodes)
print(dyads)
.emp <- tryCatch(
  .f_empile_run(
    partitions       = partitions,
    rhs              = rhs_all,
    inertial_present = TRUE,
    past_influence   = 1,
    nodes            = nodes,
    dyads            = dyads,
    group_labels     = NULL,
    directed         = FALSE,
    verbose          = TRUE
  ),
  error = function(e) {
    stop("[SECTION 1] .erpm_long_empile_build_meta_nw() failed: ", conditionMessage(e), call. = FALSE)
  }
)

cat("network.size types:\n")
tmp_nws <- lapply(partitions, function(p) .erpm_ple_partition_to_bipartite(p))
print(vapply(tmp_nws, function(nw) typeof(network::network.size(nw)), character(1)))

if (!is.list(.emp) || is.null(.emp$meta_nw) || !inherits(.emp$meta_nw, "network")) {
  stop("[SECTION 1] .erpm_long_empile_build_meta_nw() did not return a list with $meta_nw (class 'network').",
       call. = FALSE)
}

nw_meta <- .emp$meta_nw
print(.emp)

# --------------------------------------------------------------------------------------
# 1bis) Sanity-check: meta-nodes (monadic covariates) attached by the engine
# --------------------------------------------------------------------------------------
.gender <- network::get.network.attribute(nw_meta, "gender")
.age    <- network::get.network.attribute(nw_meta, "age")

if (is.null(.gender) || !is.matrix(.gender) || ncol(.gender) != 1L) {
  stop("[SECTION 1] meta-node 'gender' must be a (Nmeta x 1) matrix attached to meta_nw.", call. = FALSE)
}
if (is.null(.age) || !is.matrix(.age) || ncol(.age) != 1L) {
  stop("[SECTION 1] meta-node 'age' must be a (Nmeta x 1) matrix attached to meta_nw.", call. = FALSE)
}

cat(sprintf("[SECTION 1] meta-node gender: dim=%dx%d\n", nrow(.gender), ncol(.gender)))
cat(sprintf("[SECTION 1] meta-node age   : dim=%dx%d\n", nrow(.age), ncol(.age)))

cat("\n[SECTION 1] Meta-network built by .erpm_long_empile_build_meta_nw():\n")

# --------------------------------------------------------------------------------------
# 2) Normalize dyads contract for summary() terms (cov_fullmatch/dyadcov)
# --------------------------------------------------------------------------------------
# Your term initializers in this project sometimes read dyads from:
#   - nw %n% "dyads"  (list of vectors length nA^2), OR
#   - nw %n% "<name>" (matrix), depending on version.
#
# The PLE runner shown attaches dyads as matrices on network attributes named by nm.
# To make this selftest robust across both conventions, we ensure a canonical
# nw %n% "dyads" list exists, with vec(meta_matrix) for each dyad name.

.meta_nA <- B * n_block
if (is.null(.meta_nA) || is.na(.meta_nA)) stop("[SECTION 1] Internal error: meta_nA undefined.", call. = FALSE)

.fm_M <- network::get.network.attribute(nw_meta, "fm")
.Z1_M <- network::get.network.attribute(nw_meta, "Z1")

if (!is.matrix(.fm_M) || nrow(.fm_M) != .meta_nA || ncol(.fm_M) != .meta_nA) {
  stop("[SECTION 1] meta_nw %n% 'fm' must be a nA x nA matrix (here nA=12).", call. = FALSE)
}
if (!is.matrix(.Z1_M) || nrow(.Z1_M) != .meta_nA || ncol(.Z1_M) != .meta_nA) {
  stop("[SECTION 1] meta_nw %n% 'Z1' must be a nA x nA matrix (here nA=12).", call. = FALSE)
}

network::set.network.attribute(
  nw_meta, "dyads",
  list(fm = as.numeric(.fm_M), Z1 = as.numeric(.Z1_M))
)

dy_ret <- network::get.network.attribute(nw_meta, "dyads")
if (!is.list(dy_ret) || !all(c("fm", "Z1") %in% names(dy_ret))) {
  stop("[SECTION 1] Failed to enforce dyads list with names {fm,Z1} in nw %n% 'dyads'.", call. = FALSE)
}

# --------------------------------------------------------------------------------------
# 3) Ensure PLE past partitions attribute required by inertia_groups
# --------------------------------------------------------------------------------------
# The exposed runner may attach timeline_nws when inertial_present=TRUE, but your
# canonical inertia_groups expects erpm_block_past_partitions-style content.
# We attach it explicitly for this selftest’s oracle validation.

if (is.null(network::get.network.attribute(nw_meta, "erpm_block_past_partitions"))) {
  network::set.network.attribute(nw_meta, "erpm_block_past_partitions", erpm_block_past_partitions)
}

# --------------------------------------------------------------------------------------
# 4) Summary vs oracle
# --------------------------------------------------------------------------------------

case_specs <- list(
  list(label = "CASE S1: pi=1, size=NULL",    d = 1L, size = NULL),
  list(label = "CASE S2: pi=2, size=NULL",    d = 2L, size = NULL),
  list(label = "CASE S3: pi=1, size=2",       d = 1L, size = 2L),
  list(label = "CASE S4: pi=2, size=2",       d = 2L, size = 2L),
  list(label = "CASE S5: pi=1, size=c(1,3)",  d = 1L, size = c(1L, 3L))
)

oracle_all_effects <- function(nw, d, size_filter_int) {
  list(
    inertia_groups = oracle_inertia_groups(nw, d = d, size_filter_int = size_filter_int),
    cliques_k2     = oracle_cliques_k2(nw),
    cov_fullmatch  = oracle_within_group_upper_sum(nw, fm_stacked),
    dyadcov_Z1     = oracle_within_group_upper_sum(nw, Z1_stacked)
  )
}

for (k in seq_along(case_specs)) {
  spec <- case_specs[[k]]
  d_k <- spec$d
  size_k <- spec$size
  size_int <- if (is.null(size_k)) integer(0) else as.integer(size_k)

  cat("\n------------------------------------------------------------\n")
  cat(spec$label, "\n")

  .warn_flush()
  res <- .capture_warnings(
    summary(
      nw_meta ~
        inertia_groups(past_influence = d_k, size = size_k, debug = "deep") +
        cliques(k = 2) +
        cov_fullmatch("fm") +
        dyadcov("Z1")
    )
  )
  if (length(res$warnings)) cat(paste(res$warnings, collapse = "\n"), "\n")
  sum_obj <- res$value
  print(sum_obj)

  o <- oracle_all_effects(nw_meta, d = d_k, size_filter_int = size_int)
  svals <- as.numeric(sum_obj)

  cat(sprintf("[ORACLE] inertia_groups(pi=%d,size=%s) = %d\n",
              d_k, if (is.null(size_k)) "NULL" else paste(size_int, collapse = ","),
              as.integer(o$inertia_groups)))
  cat(sprintf("[ORACLE] cliques(k=2)                 = %g\n", o$cliques_k2))
  cat(sprintf("[ORACLE] cov_fullmatch('fm')          = %g\n", o$cov_fullmatch))
  cat(sprintf("[ORACLE] dyadcov('Z1')                = %g\n", o$dyadcov_Z1))

  cat(sprintf("[SUMMARY] inertia_groups              = %d\n", as.integer(svals[1L])))
  cat(sprintf("[SUMMARY] cliques(k=2)                = %g\n", svals[2L]))
  cat(sprintf("[SUMMARY] cov_fullmatch('fm')         = %g\n", svals[3L]))
  cat(sprintf("[SUMMARY] dyadcov('Z1')               = %g\n", svals[4L]))

  if (!isTRUE(all.equal(as.integer(svals[1L]), as.integer(o$inertia_groups))))
    stop(sprintf("Mismatch inertia_groups in %s", spec$label), call. = FALSE)
  if (!isTRUE(all.equal(svals[2L], o$cliques_k2)))
    stop(sprintf("Mismatch cliques(k=2) in %s", spec$label), call. = FALSE)
  if (!isTRUE(all.equal(svals[3L], o$cov_fullmatch)))
    stop(sprintf("Mismatch cov_fullmatch('fm') in %s", spec$label), call. = FALSE)
  if (!isTRUE(all.equal(svals[4L], o$dyadcov_Z1)))
    stop(sprintf("Mismatch dyadcov('Z1') in %s", spec$label), call. = FALSE)
}

cat("\nSECTION 1 OK: summary() matches oracle for inertia_groups + cliques + cov_fullmatch + dyadcov.\n")
# ======================================================================================
# SECTION 2) erpm_long() dry-run (PLE only) with ALL effects in one call
# Objective:
#   - validate call construction
#   - verify meta-network structure produced by erpm_long (dyads + monads) against expectations
# ======================================================================================

cat("\n================================================================================\n")
cat("SECTION 2) erpm_long() dry-run (PLE only) | ALL effects\n")
cat("================================================================================\n")

cat("\n--- DRY-RUN: inertia_groups + cliques + cov_fullmatch + dyadcov (pi=1) ---\n")

.warn_flush()
res <- .capture_warnings(
  erpm_long(
    partitions ~
      inertia_groups(past_influence = 1, debug = "deep") +
      cliques(k = 2) +
      cov_fullmatch("fm") +
      dyadcov("Z1"),
    nodes     = nodes,
    dyads     = dyads,
    mode      = "empile",
    verbose   = TRUE,
    debug     = "deep",
    eval.call = FALSE
  )
)
if (length(res$warnings)) cat(paste(res$warnings, collapse = "\n"), "\n")
out <- res$value
print(out)

# --- Validate meta-network built by erpm_long (expected kept times for pi=1: t=2,3 => B=2) ----
cat("\n[SECTION 2] Validating returned meta-network content (light checks)\n")

stopifnot(inherits(out, "erpm_long"))
stopifnot(is.list(out$network) || inherits(out$network, "network"))
nw_ret <- out$network
stopifnot(inherits(nw_ret, "network"))

# We expect (with pi=1) that erpm_long keeps times 2,3 => B_kept=2 => nA = 8, nG=6, N=14
na_ret <- network::get.network.attribute(nw_ret, "erpm_nA")
ng_ret <- network::get.network.attribute(nw_ret, "erpm_nG")
N_ret  <- network::network.size(nw_ret)

cat(sprintf("[SECTION 2] Returned: N=%d | erpm_nA=%d | erpm_nG=%d\n", N_ret, na_ret, ng_ret))

if (!identical(as.integer(na_ret), 8L) || !identical(as.integer(ng_ret), 6L) || !identical(as.integer(N_ret), 14L)) {
  stop("[SECTION 2] Unexpected stacked sizes from erpm_long() for pi=1 (expected N=14,nA=8,nG=6).", call. = FALSE)
}

# Check dyads exist and lengths match nA^2
dy_ret <- network::get.network.attribute(nw_ret, "dyads")
if (!is.list(dy_ret) || !all(c("fm", "Z1") %in% names(dy_ret))) {
  stop("[SECTION 2] Missing dyads list with names {fm,Z1} in nw %n% 'dyads'.", call. = FALSE)
}
if (length(dy_ret$fm) != 8L * 8L || length(dy_ret$Z1) != 8L * 8L) {
  stop("[SECTION 2] Dyads vectors have wrong length (expected nA^2=64).", call. = FALSE)
}

# Spot-check monadic attrs presence (x for actors, g for groups)
for (nm in c("gender", "age")) {
  M <- network::get.network.attribute(nw_ret, nm)
  if (is.null(M) || !is.matrix(M) || ncol(M) != 1L) {
    stop(sprintf("[SECTION 2] Missing meta-node '%s' as (Nmeta x 1) matrix.", nm), call. = FALSE)
  }
}

cat("[SECTION 2] OK: erpm_long dry-run returned a meta-network with expected dimensions + dyads + monads.\n")
cat("\nSECTION 2 OK (dry-run call executed).\n")

# ======================================================================================
# SECTION 3) FITS (NO blockdiag): direct ergm() on manual nw_meta with b1part only
# Objective:
#   - run minimal fits without blockdiag:
#       * inertia_groups only
#       * cliques only
#       * cov_fullmatch only
#       * dyadcov only
#       * all 4 effects together
#   - Keep cross-block edges essentially impossible WITHOUT using blockdiag()
# ======================================================================================

cat("\n================================================================================\n")
cat("SECTION 3) FITS (NO blockdiag): ergm() on manual meta-network + b1part only\n")
cat("================================================================================\n")

# Build the "forbid cross-block membership edges" edgecov matrix once
X_forbid <- .build_forbidden_crossblock_edgecov()

# Control (MPLE keeps this fast/stable; no MCMC proposal needed)
ctrl_mple <- ergm::control.ergm(
  MPLE.maxit = 50,
  MPLE.type = "glm"
)

.fit_one <- function(rhs, label) {
  cat("\n------------------------------------------------------------\n")
  cat(label, "\n")
  .warn_flush()

  # Add the "no cross-block" structural offset to keep the PLE logic sane without blockdiag()
  f <- as.formula(paste0(
    "nw_meta ~ offset(edgecov(X_forbid)) + ", rhs
  ))

  res <- .capture_warnings(
    tryCatch(
      ergm::ergm(
        f,
        constraints = ~b1part,
        # estimate = "MPLE",
        # control = ctrl_mple,
        verbose = TRUE
      ),
      error = function(e) e
    )
  )

  if (length(res$warnings)) cat(paste(res$warnings, collapse = "\n"), "\n")

  if (inherits(res$value, "error")) {
    cat("[FIT ERROR]\n")
    cat(conditionMessage(res$value), "\n")
    return(invisible(NULL))
  }

  fit <- res$value
  print(fit)
  invisible(fit)
}

# Individual fits
.fit_one(
  rhs = 'inertia_groups(past_influence = 1, debug = "deep")',
  label = "FIT 1) inertia_groups only (pi=1)"
)

.fit_one(
  rhs = "cliques(k = 2)",
  label = "FIT 2) cliques only (k=2)"
)

.fit_one(
  rhs = 'cov_fullmatch("fm")',
  label = "FIT 3) cov_fullmatch only ('fm')"
)

.fit_one(
  rhs = 'dyadcov("Z1")',
  label = "FIT 4) dyadcov only ('Z1')"
)

.fit_all <- .fit_one(
  rhs = paste(
    'inertia_groups(past_influence = 1, debug = "deep")',
    "cliques(k = 2)",
    'cov_fullmatch("fm")',
    'dyadcov("Z1")',
    sep = " + "
  ),
  label = "FIT 5) ALL effects: inertia_groups + cliques + cov_fullmatch + dyadcov"
)

cat("\nSELFTEST DONE.\n")

if (.patch_enabled && exists("ergm_patch_disable", mode = "function")) {
  try(ergm_patch_disable(), silent = TRUE)
}