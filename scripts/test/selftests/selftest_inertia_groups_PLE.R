# ==============================================================================
# File    : scripts/test/selftests/selftest_inertia_groups_PLE.R
# Auteur  : Jérémie Chichignoud - Cub'itech
# Purpose : Integration self-test for the ERPM inertial ERGM term `inertia_groups`
#           under the longitudinal stacked engine (PLE / "empile").
#
# Scope
#   - Tests the PLE longitudinal pipeline used by `erpm_long()`.
#   - Focuses on the exogenous inertial variant of `inertia_groups`.
#   - Exercises the meta-network construction and downstream ERGM evaluation.
#
# Constraints
#   - PLE ("empile") path only; the PLS engine is not tested here.
#   - Only exogenous inertial mode is considered.
#   - `blockdiag()` is intentionally excluded from constraints because the
#     combination `blockdiag + b1part` is currently known to be broken.
#   - All `summary()` / `ergm()` calls therefore use `constraints = ~ b1part`.
#
# Datasets
#   Dataset #1
#     - n = 4 actors, T = 3 partitions
#     - Dyads : {fm, Z1}
#     - Nodes : {label, gender, age}
#
#   Dataset #2
#     - n = 5 actors, T = 3 partitions
#     - Dyads : {Y, X1}
#     - Nodes : {id, sex, age_years}
#
# Test organization
#   1) DRY-RUN
#        Build the PLE meta-network and validate the `erpm()` call produced
#        by the wrapper (dataset #1 across several RHS scenarios).
#
#   2) SUMMARY
#        Build the meta-network via `erpm_long()`, run `summary()`, and compare
#        results with offline expected computations (datasets #1 and #2).
#
#   3) FIT
#        Run `erpm_long()` with estimation enabled over the same scenario grid
#        (datasets #1 and #2). Individual failures are logged but must not stop
#        execution of the whole test script.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Preamble (locale, packages, and reproducibility knobs)
# ------------------------------------------------------------------------------
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
options(keep.source = TRUE)
options(keep.source.pkgs = TRUE)
Sys.setenv(R_KEEP_PKG_SOURCE = "yes")
# devtools::load_all(".")


# ------------------------------------------------------------------------------
# Debug helpers and warning capture utilities
# ------------------------------------------------------------------------------
dbg <- TRUE
dbgcat <- function(...) if (dbg) cat("[selftest_inertia_groups_PLE][DEBUG] ", ..., "\n", sep = "")

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

.print_warnings <- function(w) {
  if (length(w)) cat(paste0("[WARN] ", w, collapse = "\n"), "\n")
}

# ------------------------------------------------------------------------------
# Optional ERGM patch: never fail the selftest if the patch is unavailable/broken
# ------------------------------------------------------------------------------
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

cat("=== SELFTEST inertia_groups (PLE only) | constraints: b1part only | NO blockdiag ===\n")

# ==============================================================================
# DATASET #1 (n=4, T=3): explicit nodes, dyads, and partitions
# ==============================================================================
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

# ==============================================================================
# DATASET #2 (n=5, T=3): explicit nodes, dyads, and partitions
# ==============================================================================
# Nodes (explicit, per time) — size = 5
nodes2 <- list(
  data.frame(
    id        = c("E", "F", "G", "H", "I"),
    sex       = c(1, 2, 1, 2, 1),
    age_years = c(21, 24, 29, 31, 26)
  ),
  data.frame(
    id        = c("KA", "LU", "MI", "NO", "PA"),
    sex       = c(2, 2, 1, 1, 2),
    age_years = c(35, 28, 22, 40, 33)
  ),
  data.frame(
    id        = c("Q", "R", "T", "U", "V"),
    sex       = c(1, 1, 2, 2, 1),
    age_years = c(27, 26, 34, 29, 31)
  )
)

# Dyads (explicit, per time) — 5x5
dyads2 <- list(
  list(
    Y = matrix(c(
      0, 1, 0, 1, 0,
      1, 0, 1, 0, 0,
      0, 1, 0, 1, 1,
      1, 0, 1, 0, 0,
      0, 0, 1, 0, 0
    ), nrow = 5, byrow = TRUE),
    X1 = matrix(c(
      0, 2, 0, 3, 0,
      2, 0, 1, 0, 0,
      0, 1, 0, 4, 2,
      3, 0, 4, 0, 1,
      0, 0, 2, 1, 0
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    Y = matrix(c(
      0, 1, 1, 0, 0,
      1, 0, 0, 1, 0,
      1, 0, 0, 1, 1,
      0, 1, 1, 0, 0,
      0, 0, 1, 0, 0
    ), nrow = 5, byrow = TRUE),
    X1 = matrix(c(
      0, 3, 2, 0, 0,
      3, 0, 0, 1, 0,
      2, 0, 0, 4, 5,
      0, 1, 4, 0, 0,
      0, 0, 5, 0, 0
    ), nrow = 5, byrow = TRUE)
  ),
  list(
    Y = matrix(c(
      0, 0, 1, 1, 0,
      0, 0, 1, 0, 1,
      1, 1, 0, 0, 1,
      1, 0, 0, 0, 0,
      0, 1, 1, 0, 0
    ), nrow = 5, byrow = TRUE),
    X1 = matrix(c(
      0, 1, 2, 3, 0,
      1, 0, 4, 0, 2,
      2, 4, 0, 5, 6,
      3, 0, 5, 0, 0,
      0, 2, 6, 0, 0
    ), nrow = 5, byrow = TRUE)
  )
)

# Partitions passed to erpm_long — length = 5
partitions2 <- list(
  c(1, 1, 2, 3, 3),
  c(1, 2, 2, 3, 1),
  c(2, 2, 3, 3, 1)
)

# ==============================================================================
# Offline helpers (expected values for summary cross-checks)
#
# Notation:
#   - cliques(k=2): sum_g choose(n_g, 2)
#   - cov_match(<attr>, k=2, normalized="none"): sum_g sum_r choose(n_{g,r}, 2)
#   - dyadcov(<name>, k=2, normalize=FALSE): sum_g sum_{i<j in g} (Zij + Zji)
#   - inertia_groups(pi=d): strict intersection across lags 1..d (PLE storage convention)
# ==============================================================================
.get_group_members_actor_ids <- function(nw, gv, n1) {
  nb <- network::get.neighborhood(nw, gv, type = "all")
  ids <- sort(unique(as.integer(nb)))
  ids[ids >= 1L & ids <= n1]
}

.expected_cliques_k2 <- function(nw) {
  n1 <- as.integer(nw %n% "bipartite")
  if (!is.finite(n1) || n1 <= 0L) stop("[expected] missing/invalid bipartite.")
  tot <- 0
  for (gv in (n1 + 1L):(2L * n1)) {
    m <- length(.get_group_members_actor_ids(nw, gv, n1))
    if (m >= 2L) tot <- tot + (m * (m - 1L)) / 2
  }
  tot
}

.expected_cov_match_k2_none <- function(nw, cov_attr) {
  n1 <- as.integer(nw %n% "bipartite")
  x <- network::get.vertex.attribute(nw, cov_attr)
  if (is.null(x)) stop("[expected] missing vertex attribute: ", cov_attr)
  x <- x[seq_len(n1)]
  tot <- 0
  for (gv in (n1 + 1L):(2L * n1)) {
    ids <- .get_group_members_actor_ids(nw, gv, n1)
    if (length(ids) < 2L) next
    tab <- table(x[ids], useNA = "no")
    if (!length(tab)) next
    tot <- tot + sum(tab * (tab - 1L) / 2)
  }
  tot
}

.expected_dyadcov_k2_raw <- function(nw, dyad_name) {
  n1 <- as.integer(nw %n% "bipartite")
  Z <- tryCatch(nw %n% dyad_name, error = function(e) NULL)

  if (is.null(Z)) {
    dy <- tryCatch(nw %n% "dyads", error = function(e) NULL)
    if (is.list(dy) && !is.null(dy[[dyad_name]])) Z <- dy[[dyad_name]]
  }
  if (is.null(Z) || !is.matrix(Z)) stop("[expected] dyad matrix not found: ", dyad_name)

  if (nrow(Z) < n1 || ncol(Z) < n1) stop("[expected] dyad dims too small for n1.")
  if (nrow(Z) != n1 || ncol(Z) != n1) Z <- Z[seq_len(n1), seq_len(n1), drop = FALSE]

  tot <- 0
  for (gv in (n1 + 1L):(2L * n1)) {
    ids <- .get_group_members_actor_ids(nw, gv, n1)
    if (length(ids) < 2L) next
    for (ii in 1L:(length(ids) - 1L)) {
      for (jj in (ii + 1L):length(ids)) {
        i <- ids[ii]
        j <- ids[jj]
        tot <- tot + (Z[i, j] + Z[j, i])
      }
    }
  }
  tot
}

.expected_inertia_groups <- function(nw, past_influence, size = NULL) {
  d <- as.integer(past_influence)
  if (!is.finite(d) || d < 1L) stop("[expected] past_influence must be >= 1.")

  n1 <- as.integer(nw %n% "bipartite")
  if (!is.finite(n1) || n1 <= 0L) stop("[expected] missing/invalid bipartite.")

  B <- tryCatch(as.integer(nw %n% "erpm_B"), error = function(e) NA_integer_)
  n_block <- tryCatch(as.integer(nw %n% "erpm_n"), error = function(e) NA_integer_)
  G_block <- tryCatch(as.integer(nw %n% "erpm_G"), error = function(e) NA_integer_)
  past <- tryCatch(nw %n% "erpm_block_past_partitions", error = function(e) NULL)

  if (!is.finite(B) || !is.finite(n_block) || !is.finite(G_block) || is.null(past)) {
    stop("[expected] missing inertia_groups PLE attributes: erpm_B/erpm_n/erpm_G/erpm_block_past_partitions.")
  }
  if (n1 != B * n_block) {
    stop("[expected] inconsistent sizes: bipartite n1=", n1, " but B*n=", B * n_block)
  }
  if (!is.list(past) || length(past) != B) stop("[expected] past container must have length B.")

  sizes_int <- if (is.null(size)) integer(0) else sort(unique(as.integer(size)))
  cnt <- 0L

  for (b in seq_len(B)) {
    if (!is.list(past[[b]]) || length(past[[b]]) < d) stop("[expected] past[[b]] must have at least d partitions.")
    for (g in seq_len(G_block)) {
      gv <- n1 + (b - 1L) * G_block + g
      cur_ids <- .get_group_members_actor_ids(nw, gv, n1)
      if (!length(cur_ids)) next
      if (length(sizes_int) && !(length(cur_ids) %in% sizes_int)) next

      ok_all_lags <- TRUE
      for (lag in seq_len(d)) {
        p_lag <- past[[b]][[lag]]
        if (!is.atomic(p_lag) || length(p_lag) != n_block) stop("[expected] past partition wrong shape.")
        groups <- split(seq_along(p_lag), as.integer(p_lag))
        groups_global <- lapply(groups, function(v) sort((b - 1L) * n_block + as.integer(v)))

        ok_lag <- FALSE
        for (u in seq_along(groups_global)) {
          if (length(groups_global[[u]]) == length(cur_ids) && all(groups_global[[u]] == cur_ids)) {
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

# ------------------------------------------------------------------------------
# Locate the PLE meta-network builder (internal engine entry point)
# ------------------------------------------------------------------------------
.f_empile_build <- get0(".erpm_long_empile_build_meta_nw", mode = "function", inherits = TRUE)
if (is.null(.f_empile_build) && "ERPM" %in% loadedNamespaces()) {
  .f_empile_build <- get0(".erpm_long_empile_build_meta_nw", envir = asNamespace("ERPM"),
                          mode = "function", inherits = FALSE)
}
if (is.null(.f_empile_build)) stop("[selftest] Cannot find .erpm_long_empile_build_meta_nw().", call. = FALSE)

# ==============================================================================
# SECTION 1) DRY-RUN (dataset #1)
#
# Goal:
#   - exercise the PLE builder directly and validate a few structural invariants
#   - check that erpm_long(eval.call=TRUE) returns a sane erpm() call for the same RHS
#
# Scenario grid (dataset #1), as requested:
#   1) partitions ok ; nodes=NULL ; dyads=NULL ; cliques
#   2) partitions ok ; nodes=NULL ; dyads=NULL ; inertia_groups
#   3) partitions ok ; nodes=NULL ; dyads=NULL ; cliques + inertia_groups
#   4) partitions ok ; nodes=filled ; dyads=NULL ; cov_match
#   5) partitions ok ; nodes=filled ; dyads=NULL ; inertia_groups
#   6) partitions ok ; nodes=filled ; dyads=NULL ; cov_match + inertia_groups
#   7) partitions ok ; nodes=filled ; dyads=filled ; dyadcov
#   8) partitions ok ; nodes=filled ; dyads=filled ; inertia_groups
#   9) partitions ok ; nodes=filled ; dyads=filled ; dyadcov + inertia_groups
# ==============================================================================
cat("\n================================================================================\n")
cat("SECTION 1) DRY-RUN: build meta-network + verify erpm() call (dataset #1)\n")
cat("================================================================================\n")

.run_dry_scenario <- function(label, rhs, partitions, nodes_arg, dyads_arg) {
  cat("\n------------------------------------------------------------\n")
  cat(label, "\n")
  dbgcat("rhs: ", deparse(rhs))
  dbgcat("nodes: ", if (is.null(nodes_arg)) "NULL" else "filled")
  dbgcat("dyads: ", if (is.null(dyads_arg)) "NULL" else "filled")

  inert <- grepl("\\binertia_groups\\b", paste(deparse(rhs), collapse = " "))
  d <- if (inert) 1L else 0L

  # Build the PLE meta-network via the engine builder (structural validation only)
  built <- tryCatch(
    .f_empile_build(
      partitions       = partitions,
      rhs              = rhs,
      inertial_present = inert,
      past_influence   = d,
      nodes            = nodes_arg,
      dyads            = dyads_arg,
      group_labels     = NULL,
      directed         = FALSE,
      verbose          = TRUE,
      debug            = TRUE
    ),
    error = function(e) e
  )

  if (inherits(built, "error")) {
    cat("[ENGINE ERROR]\n", conditionMessage(built), "\n")
    return(invisible(NULL))
  }

  nw <- built$meta_nw
  stopifnot(inherits(nw, "network"))
  n1 <- as.integer(nw %n% "bipartite")
  N  <- network::network.size(nw)

  dbgcat("meta network: N=", N, " | n1(bipartite)=", n1, " | selected_partition_indices={", paste(built$selected_partition_indices, collapse = ","), "}")

  # Sanity-check vertex timeblock attribute (expected for stacked networks)
  tb <- network::get.vertex.attribute(nw, "timeblock")
  if (is.null(tb) || length(tb) != N) stop("[DRY] missing/invalid vertex attr 'timeblock'.", call. = FALSE)

  if (!is.null(nodes_arg)) {
    # When nodes are provided, their covariates must be materialized as vertex attributes
    if (is.null(network::get.vertex.attribute(nw, "gender"))) stop("[DRY] missing vertex attr 'gender'.", call. = FALSE)
    if (is.null(network::get.vertex.attribute(nw, "age")))    stop("[DRY] missing vertex attr 'age'.", call. = FALSE)
  }

  if (!is.null(dyads_arg)) {
    # Dyads may be attached either directly (legacy) or under nw %n% "dyads" (preferred)
    z1 <- tryCatch(nw %n% "Z1", error = function(e) NULL)
    dy <- tryCatch(nw %n% "dyads", error = function(e) NULL)
    okZ1 <- is.matrix(z1) || (is.list(dy) && !is.null(dy[["Z1"]]))
    if (!okZ1) stop("[DRY] dyads filled but cannot find Z1 in nw %n% 'Z1' nor nw %n% 'dyads'[['Z1']].", call. = FALSE)
  }

  if (inert) {
    # inertia_groups in PLE requires a specific set of network attributes
    if (is.null(tryCatch(nw %n% "erpm_block_past_partitions", error = function(e) NULL)))
      stop("[DRY] missing nw %n% 'erpm_block_past_partitions' (required by inertia_groups).", call. = FALSE)
    if (!identical(as.character(tryCatch(nw %n% "erpm_mode", error = function(e) "")), "empile"))
      stop("[DRY] missing/invalid nw %n% 'erpm_mode' for PLE.", call. = FALSE)
  }

  # Verify the call object produced by erpm_long(eval.call=TRUE) for the same RHS
  call_erpm <- tryCatch(
    erpm_long(
      partitions ~ rhs,
      nodes     = nodes_arg,
      dyads     = dyads_arg,
      mode      = "empile",
      verbose   = TRUE,
      debug     = TRUE,
      eval.call = TRUE
    ),
    error = function(e) e
  )
  if (inherits(call_erpm, "error")) {
    cat("[ERPM_LONG DRY ERROR]\n", conditionMessage(call_erpm), "\n")
    return(invisible(NULL))
  }

  dbgcat("erpm_long returned call (eval.call=TRUE):")
  print(call_erpm)

  cat("[DRY] OK\n")
  invisible(list(meta_nw = nw, call = call_erpm))
}

# Build RHS expressions explicitly to keep term ordering deterministic
rhs1 <- quote(cliques(k = 2))
rhs2 <- quote(inertia_groups(past_influence = 1))
rhs3 <- quote(cliques(k = 2) + inertia_groups(past_influence = 1))
rhs4 <- quote(cov_match("gender", clique_size = 2, normalized = "none"))
rhs5 <- quote(inertia_groups(past_influence = 1))
rhs6 <- quote(cov_match("gender", clique_size = 2, normalized = "none") + inertia_groups(past_influence = 1))
rhs7 <- quote(dyadcov("Z1", clique_size = 2, normalize = FALSE))
rhs8 <- quote(inertia_groups(past_influence = 1))
rhs9 <- quote(dyadcov("Z1", clique_size = 2, normalize = FALSE) + inertia_groups(past_influence = 1))

.run_dry_scenario("S1.1) nodes=NULL; dyads=NULL; cliques",                      rhs1, partitions, NULL,  NULL)
.run_dry_scenario("S1.2) nodes=NULL; dyads=NULL; inertia_groups",              rhs2, partitions, NULL,  NULL)
.run_dry_scenario("S1.3) nodes=NULL; dyads=NULL; cliques + inertia_groups",    rhs3, partitions, NULL,  NULL)
.run_dry_scenario("S1.4) nodes=filled; dyads=NULL; cov_match",                 rhs4, partitions, nodes, NULL)
.run_dry_scenario("S1.5) nodes=filled; dyads=NULL; inertia_groups",            rhs5, partitions, nodes, NULL)
.run_dry_scenario("S1.6) nodes=filled; dyads=NULL; cov_match + inertia_groups",rhs6, partitions, nodes, NULL)
.run_dry_scenario("S1.7) nodes=filled; dyads=filled; dyadcov",                 rhs7, partitions, nodes, dyads)
.run_dry_scenario("S1.8) nodes=filled; dyads=filled; inertia_groups",          rhs8, partitions, nodes, dyads)
.run_dry_scenario("S1.9) nodes=filled; dyads=filled; dyadcov + inertia_groups",rhs9, partitions, nodes, dyads)

cat("\nSECTION 1 DONE.\n")

# ==============================================================================
# SECTION 3) SUMMARY vs OFFLINE COMPUTATION
#
# For each dataset:
#   - build the PLE meta-network via erpm_long(eval.call=FALSE)
#   - run summary(nw ~ terms, constraints = ~ b1part)
#   - compute expected values on the *exact same* network object and compare
#
# Scenario grid (for each dataset), as requested:
#   A) nodes=NULL; dyads=NULL; cliques
#   B) nodes=NULL; dyads=NULL; inertia_groups
#   C) nodes=NULL; dyads=NULL; cliques + inertia_groups
#   D) nodes=filled; dyads=NULL; cov_match
#   E) nodes=filled; dyads=NULL; inertia_groups
#   F) nodes=filled; dyads=NULL; cov_match + inertia_groups
#   G) nodes=filled; dyads=filled; dyadcov
#   H) nodes=filled; dyads=filled; inertia_groups
#   I) nodes=filled; dyads=filled; dyadcov + inertia_groups
#   J) nodes=filled; dyads=filled; cliques + cov_match + dyadcov + inertia_groups
# ==============================================================================
cat("\n================================================================================\n")
cat("SECTION 3) SUMMARY vs OFFLINE COMPUTATION (datas5====================================\n")

.run_summary_case <- function(label, partitions, nodes_arg, dyads_arg, rhs_expr,
                             cov_attr = NULL, dyad_name = NULL, past_influence = 1L) {
  cat("\n------------------------------------------------------------\n")
  cat(label, "\n")
  dbgcat("nodes: ", if (is.null(nodes_arg)) "NULL" else "filled",
         " | dyads: ", if (is.null(dyads_arg)) "NULL" else "filled",
         " | rhs: ", deparse(rhs_expr))

  # Build via erpm_long (we only need the constructed network, not a fit object)
  formula <- as.formula(call("~", quote(partitions), rhs_expr))
  
  .warn_flush()
  res <- .capture_warnings(
    tryCatch(
      erpm_long(
        formula   = formula,
        nodes     = nodes_arg,
        dyads     = dyads_arg,
        mode      = "empile",
        verbose   = TRUE,
        debug     = TRUE,
        eval.call = FALSE
      ),
      error = function(e) e
    )
  )
  .print_warnings(res$warnings)

  if (inherits(res$value, "error")) {
    cat("[ERPM_LONG ERROR]\n", conditionMessage(res$value), "\n")
    return(invisible(NULL))
  }

  out <- res$value
  stopifnot(is.list(out) || inherits(out, "erpm_long"))
  nw <- out$network
  stopifnot(inherits(nw, "network"))

  n1 <- as.integer(nw %n% "bipartite")
  N  <- network::network.size(nw)
  dbgcat("built nw: N=", N, " | n1=", n1)

  # Run summary() under b1part only (no blockdiag in this selftest)
  .warn_flush()
  sres <- .capture_warnings(
    tryCatch(
      summary(
        as.formula(call("~", nw, rhs_expr)),
        constraints = ~ b1part
      ),
      error = function(e) e
    )
  )
  .print_warnings(sres$warnings)
  if (inherits(sres$value, "error")) {
    cat("[SUMMARY ERROR]\n", conditionMessage(sres$value), "\n")
    return(invisible(NULL))
  }

  s <- sres$value
  print(s)

  # Compute offline expected stats, matching the RHS term order used above
  rhs_txt <- paste(deparse(rhs_expr), collapse = " ")

  want_cliques <- grepl("\\bcliques\\b", rhs_txt)
  want_cov     <- grepl("\\bcov_match\\b", rhs_txt)
  want_dyad    <- grepl("\\bdyadcov\\b", rhs_txt)
  want_inert   <- grepl("\\binertia_groups\\b", rhs_txt)

  o <- numeric(0)
  if (want_cliques) o <- c(o, .expected_cliques_k2(nw))
  if (want_cov)     o <- c(o, .expected_cov_match_k2_none(nw, cov_attr))
  if (want_dyad)    o <- c(o, .expected_dyadcov_k2_raw(nw, dyad_name))
  if (want_inert)   o <- c(o, .expected_inertia_groups(nw, past_influence = past_influence))

  sval <- as.numeric(s)

  cat("[OFFLINE]  ", paste(signif(o, 10), collapse = " | "), "\n")
  cat("[SUMMARY]  ", paste(signif(sval, 10), collapse = " | "), "\n")

  if (length(o) != length(sval)) {
    stop("[MISMATCH] expected length != summary length (check rhs ordering).", call. = FALSE)
  }
  if (!isTRUE(all.equal(sval, o))) {
    stop("[MISMATCH] summary != expected in: ", label, call. = FALSE)
  }

  cat("[SUMMARY vs OFFLINE] OK\n")
  invisible(list(out = out, nw = nw, summary = s))
}

# ---------------------------
# Dataset #1: case RHS blocks
# ---------------------------
rhs_1A <- quote(cliques(k = 2))
rhs_1B <- quote(inertia_groups(past_influence = 1))
rhs_1C <- quote(cliques(k = 2) + inertia_groups(past_influence = 1))
rhs_1D <- quote(cov_match("gender", clique_size = 2, normalized = "none"))
rhs_1E <- quote(inertia_groups(past_influence = 1))
rhs_1F <- quote(cov_match("gender", clique_size = 2, normalized = "none") + inertia_groups(past_influence = 1))
rhs_1G <- quote(dyadcov("Z1", clique_size = 2, normalize = FALSE))
rhs_1H <- quote(inertia_groups(past_influence = 1))
rhs_1I <- quote(dyadcov("Z1", clique_size = 2, normalize = FALSE) + inertia_groups(past_influence = 1))
rhs_1J <- quote(cliques(k = 2) + cov_match("gender", clique_size = 2, normalized = "none") +
                  dyadcov("Z1", clique_size = 2, normalize = FALSE) +
                  inertia_groups(past_influence = 1))

# Run dataset #1 grid
cat("\n--- DATASET #1 (n=4) ---\n")
.run_summary_case("S3.1A) nodes=NULL; dyads=NULL; cliques",
                  partitions, NULL, NULL, rhs_1A)

.run_summary_case("S3.1B) nodes=NULL; dyads=NULL; inertia_groups",
                  partitions, NULL, NULL, rhs_1B, past_influence = 1L)

.run_summary_case("S3.1C) nodes=NULL; dyads=NULL; cliques + inertia_groups",
                  partitions, NULL, NULL, rhs_1C, past_influence = 1L)

.run_summary_case("S3.1D) nodes=filled; dyads=NULL; cov_match(gender)",
                  partitions, nodes, NULL, rhs_1D, cov_attr = "gender")

.run_summary_case("S3.1E) nodes=filled; dyads=NULL; inertia_groups",
                  partitions, nodes, NULL, rhs_1E, past_influence = 1L)

.run_summary_case("S3.1F) nodes=filled; dyads=NULL; cov_match + inertia_groups",
                  partitions, nodes, NULL, rhs_1F, cov_attr = "gender", past_influence = 1L)

.run_summary_case("S3.1G) nodes=filled; dyads=filled; dyadcov(Z1)",
                  partitions, nodes, dyads, rhs_1G, dyad_name = "Z1")

.run_summary_case("S3.1H) nodes=filled; dyads=filled; inertia_groups",
                  partitions, nodes, dyads, rhs_1H, past_influence = 1L)

.run_summary_case("S3.1I) nodes=filled; dyads=filled; dyadcov + inertia_groups",
                  partitions, nodes, dyads, rhs_1I, dyad_name = "Z1", past_influence = 1L)

.run_summary_case("S3.1J) nodes=filled; dyads=filled; cliques + cov_match + dyadcov + inertia_groups",
                  partitions, nodes, dyads, rhs_1J, cov_attr = "gender", dyad_name = "Z1", past_influence = 1L)

# ---------------------------
# Dataset #2: case RHS blocks
#   - cov_match uses "sex"
#   - dyadcov uses "X1"
# ---------------------------
rhs_2A <- quote(cliques(k = 2))
rhs_2B <- quote(inertia_groups(past_influence = 1))
rhs_2C <- quote(cliques(k = 2) + inertia_groups(past_influence = 1))
rhs_2D <- quote(cov_match("sex", clique_size = 2, normalized = "none"))
rhs_2E <- quote(inertia_groups(past_influence = 1))
rhs_2F <- quote(cov_match("sex", clique_size = 2, normalized = "none") + inertia_groups(past_influence = 1))
rhs_2G <- quote(dyadcov("X1", clique_size = 2, normalize = FALSE))
rhs_2H <- quote(inertia_groups(past_influence = 1))
rhs_2I <- quote(dyadcov("X1", clique_size = 2, normalize = FALSE) + inertia_groups(past_influence = 1))
rhs_2J <- quote(cliques(k = 2) + cov_match("sex", clique_size = 2, normalized = "none") +
                  dyadcov("X1", clique_size = 2, normalize = FALSE) +
                  inertia_groups(past_influence = 1))

# Run dataset #2 grid
cat("\n--- DATASET #2 (n=5) ---\n")
.run_summary_case("S3.2A) nodes=NULL; dyads=NULL; cliques",
                  partitions2, NULL, NULL, rhs_2A)

.run_summary_case("S3.2B) nodes=NULL; dyads=NULL; inertia_groups",
                  partitions2, NULL, NULL, rhs_2B, past_influence = 1L)

.run_summary_case("S3.2C) nodes=NULL; dyads=NULL; cliques + inertia_groups",
                  partitions2, NULL, NULL, rhs_2C, past_influence = 1L)

.run_summary_case("S3.2D) nodes=filled; dyads=NULL; cov_match(sex)",
                  partitions2, nodes2, NULL, rhs_2D, cov_attr = "sex")

.run_summary_case("S3.2E) nodes=filled; dyads=NULL; inertia_groups",
                  partitions2, nodes2, NULL, rhs_2E, past_influence = 1L)

.run_summary_case("S3.2F) nodes=filled; dyads=NULL; cov_match + inertia_groups",
                  partitions2, nodes2, NULL, rhs_2F, cov_attr = "sex", past_influence = 1L)

.run_summary_case("S3.2G) nodes=filled; dyads=filled; dyadcov(X1)",
                  partitions2, nodes2, dyads2, rhs_2G, dyad_name = "X1")

.run_summary_case("S3.2H) nodes=filled; dyads=filled; inertia_groups",
                  partitions2, nodes2, dyads2, rhs_2H, past_influence = 1L)

.run_summary_case("S3.2I) nodes=filled; dyads=filled; dyadcov + inertia_groups",
                  partitions2, nodes2, dyads2, rhs_2I, dyad_name = "X1", past_influence = 1L)

.run_summary_case("S3.2J) nodes=filled; dyads=filled; cliques + cov_match + dyadcov + inertia_groups",
                  partitions2, nodes2, dyads2, rhs_2J, cov_attr = "sex", dyad_name = "X1", past_influence = 1L)

cat("\nSECTION 3 DONE.\n")

# ==============================================================================
# SECTION 4) FITS via erpm_long() (eval.call=FALSE)
#
# Notes:
#   - No explicit control/estimate tuning here: let ergm defaults apply.
#   - Only b1part is used internally (blockdiag is intentionally avoided here).
#   - This section is allowed to fail on some scenarios (separation, constant stats, etc.),
#     but the script must keep going: errors are logged and swallowed.
# ==============================================================================
cat("\n================================================================================\n")
cat("SECTION 4) FITS: erpm_long() eval.call=FALSE (datasets #1 and #2)\n")
cat("================================================================================\n")

.run_fit_case <- function(label, partitions, nodes_arg, dyads_arg, rhs_expr) {
  cat("\n------------------------------------------------------------\n")
  cat(label, "\n")
  dbgcat("nodes: ", if (is.null(nodes_arg)) "NULL" else "filled",
         " | dyads: ", if (is.null(dyads_arg)) "NULL" else "filled",
         " | rhs: ", deparse(rhs_expr))

  .warn_flush()
  res <- .capture_warnings(
    tryCatch(
      erpm_long(
        partitions ~ rhs_expr,
        nodes     = nodes_arg,
        dyads     = dyads_arg,
        mode      = "empile",
        verbose   = TRUE,
        debug     = NULL,
        eval.call = FALSE
      ),
      error = function(e) e
    )
  )
  .print_warnings(res$warnings)

  if (inherits(res$value, "error")) {
    cat("[FIT ERROR]\n", conditionMessage(res$value), "\n")
    return(invisible(NULL))
  }

  fit <- res$value
  print(fit)
  cat("[FIT] OK\n")
  invisible(fit)
}

# Dataset #1 fits
cat("\n--- DATASET #1 FITS (n=4) ---\n")
.run_fit_case("S4.1A) nodes=NULL; dyads=NULL; cliques",                      partitions, NULL,  NULL,  rhs_1A)
.run_fit_case("S4.1B) nodes=NULL; dyads=NULL; inertia_groups",              partitions, NULL,  NULL,  rhs_1B)
.run_fit_case("S4.1C) nodes=NULL; dyads=NULL; cliques + inertia_groups",    partitions, NULL,  NULL,  rhs_1C)
.run_fit_case("S4.1D) nodes=filled; dyads=NULL; cov_match",                 partitions, nodes, NULL,  rhs_1D)
.run_fit_case("S4.1E) nodes=filled; dyads=NULL; inertia_groups",            partitions, nodes, NULL,  rhs_1E)
.run_fit_case("S4.1F) nodes=filled; dyads=NULL; cov_match + inertia_groups",partitions, nodes, NULL,  rhs_1F)
.run_fit_case("S4.1G) nodes=filled; dyads=filled; dyadcov",                 partitions, nodes, dyads, rhs_1G)
.run_fit_case("S4.1H) nodes=filled; dyads=filled; inertia_groups",          partitions, nodes, dyads, rhs_1H)
.run_fit_case("S4.1I) nodes=filled; dyads=filled; dyadcov + inertia_groups",partitions, nodes, dyads, rhs_1I)
.run_fit_case("S4.1J) nodes=filled; dyads=filled; ALL",                     partitions, nodes, dyads, rhs_1J)

# Dataset #2 fits
cat("\n--- DATASET #2 FITS (n=5) ---\n")
.run_fit_case("S4.2A) nodes=NULL; dyads=NULL; cliques",                      partitions2, NULL,   NULL,   rhs_2A)
.run_fit_case("S4.2B) nodes=NULL; dyads=NULL; inertia_groups",              partitions2, NULL,   NULL,   rhs_2B)
.run_fit_case("S4.2C) nodes=NULL; dyads=NULL; cliques + inertia_groups",    partitions2, NULL,   NULL,   rhs_2C)
.run_fit_case("S4.2D) nodes=filled; dyads=NULL; cov_match",                 partitions2, nodes2, NULL,   rhs_2D)
.run_fit_case("S4.2E) nodes=filled; dyads=NULL; inertia_groups",            partitions2, nodes2, NULL,   rhs_2E)
.run_fit_case("S4.2F) nodes=filled; dyads=NULL; cov_match + inertia_groups",partitions2, nodes2, NULL,   rhs_2F)
.run_fit_case("S4.2G) nodes=filled; dyads=filled; dyadcov",                 partitions2, nodes2, dyads2, rhs_2G)
.run_fit_case("S4.2H) nodes=filled; dyads=filled; inertia_groups",          partitions2, nodes2, dyads2, rhs_2H)
.run_fit_case("S4.2I) nodes=filled; dyads=filled; dyadcov + inertia_groups",partitions2, nodes2, dyads2, rhs_2I)
.run_fit_case("S4.2J) nodes=filled; dyads=filled; ALL",                     partitions2, nodes2, dyads2, rhs_2J)

cat("\nSELFTEST DONE.\n")

# ------------------------------------------------------------------------------
# Cleanup: disable patch if it was enabled (best-effort only)
# ------------------------------------------------------------------------------
if (.patch_enabled && exists("ergm_patch_disable", mode = "function")) {
  try(ergm_patch_disable(), silent = TRUE)
}