
# ======================================================================================
# File    : scripts/test/selftests/selftest_inertia_groups.R
# Object  : Self-test (PLE only) for ERPM inertial term `inertia_groups`
# Run     : Rscript scripts/test/selftests/selftest_inertia_groups.R
#
# Notes
#   - PLE ("empile") only. No PLS. No inertial endogenous variant.
#   - Two explicit datasets:
#       * Dataset #1: n=4, T=3, dyads {fm, Z1}, nodes {label, gender, age}
#       * Dataset #2: n=5, T=3, dyads {Y, X1}, nodes {label, sex, age_years}
#   - The file is organized as:
#       SECTION 1) DRY-RUN: build meta-network + verify erpm() call (dataset #1)
#       SECTION 2) SUMMARY: summary() vs offline expected (datasets #1 and #2)
#                  + past_influence=2 + b1partblockdiag constraint regression
#       SECTION 3) FITS: run actual ergm fits (via erpm_long eval.call=TRUE + eval())
#       SECTION 4) b1partblockdiag: constraint structure + summary identity test
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

options(keep.source = TRUE)
options(keep.source.pkgs = TRUE)
Sys.setenv(R_KEEP_PKG_SOURCE = "yes")
# devtools::load_all(".")

# --------------------------------------------------------------------------------------
# Debug + warnings capture (selftest)
# --------------------------------------------------------------------------------------
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

# --------------------------------------------------------------------------------------
# Global pass/fail counters
# --------------------------------------------------------------------------------------
.st <- new.env(parent = emptyenv())
.st$ok   <- 0L
.st$fail <- 0L
.st_ok   <- function() { .st$ok   <- .st$ok   + 1L; invisible(NULL) }
.st_fail <- function() { .st$fail <- .st$fail + 1L; invisible(NULL) }

cat("=== SELFTEST inertia_groups (PLE only) ===\n")

# ======================================================================================
# DATASET #1 (n=4, T=3)
# ======================================================================================

# Nodes (explicit, per time)
nodes <- list(
  data.frame(
    label  = c("A", "B", "C", "D"),
    gender = c("H", "H", "F", "H"),
    age    = c(20, 22, 25, 30),
    stringsAsFactors = FALSE
  ),
  data.frame(
    label  = c("FT", "AZ", "JI", "DO"),
    gender = c("F", "H", "F", "F"),
    age    = c(10, 42, 25, 30),
    stringsAsFactors = FALSE
  ),
  data.frame(
    label  = c("H", "Z", "S", "A"),
    gender = c("H", "H", "H", "H"),
    age    = c(27, 26, 25, 28),
    stringsAsFactors = FALSE
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

# ======================================================================================
# DATASET #2 (n=5, T=3)
# ======================================================================================

# Nodes (explicit, per time) — size = 5
nodes2 <- list(
  data.frame(
    label     = c("E", "F", "G", "H", "I"),
    sex       = c("H", "F", "H", "F", "H"),
    age_years = c(21, 24, 29, 31, 26),
    stringsAsFactors = FALSE
  ),
  data.frame(
    label     = c("KA", "LU", "MI", "NO", "PA"),
    sex       = c("F", "F", "H", "H", "F"),
    age_years = c(35, 28, 22, 40, 33),
    stringsAsFactors = FALSE
  ),
  data.frame(
    label     = c("Q", "R", "T", "U", "V"),
    sex       = c("H", "H", "F", "F", "H"),
    age_years = c(27, 26, 34, 29, 31),
    stringsAsFactors = FALSE
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

# ======================================================================================
# Offline helpers (expected)
#   - cliques(k=2): sum_g choose(n_g, 2)
#   - cov_match("gender"/"sex", k=2, normalized="none"): sum_g sum_r choose(n_{g,r}, 2)
#   - dyadcov("Z1"/"X1", k=2, normalize=FALSE): sum_g sum_{i<j in g} (Zij + Zji)
#   - inertia_groups(pi=d): count actors whose current group membership matches
#     at least one past group signature in ANY lag in 1..d (union over lags)
# ======================================================================================

make_formula <- function(partitions_obj, rhs_expr) {
  env <- list2env(list(partitions = partitions_obj), parent = parent.frame())
  as.formula(call("~", quote(partitions), rhs_expr), env = env)
}

make_summary_formula <- function(nw, rhs_expr) {
  env <- list2env(list(nw = nw), parent = parent.frame())
  as.formula(call("~", quote(nw), rhs_expr), env = env)
}

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

  # Derive block structure from engine bookkeeping attributes (no erpm_B/n/G).
  sel_idx  <- tryCatch(as.integer(nw %n% "erpm_long.selected_partition_indices"),
                       error = function(e) NULL)
  nbr_by_t <- tryCatch(as.integer(nw %n% "erpm_long.nbr_actors_by_t"),
                       error = function(e) NULL)
  past     <- tryCatch(nw %n% "erpm_block_past_partitions", error = function(e) NULL)

  if (is.null(sel_idx) || is.null(nbr_by_t)) {
    stop("[expected] missing erpm_long.selected_partition_indices or erpm_long.nbr_actors_by_t.")
  }
  if (is.null(past)) {
    stop("[expected] missing erpm_block_past_partitions.")
  }

  B             <- length(sel_idx)
  n_b           <- nbr_by_t[sel_idx]                                  # actor count per block
  actor_offsets <- c(0L, cumsum(n_b))[seq_len(B)]                     # 0-based

  if (sum(n_b) != n1) {
    stop("[expected] inconsistent sizes: bipartite n1=", n1,
         " but sum(n_b)=", sum(n_b))
  }
  if (!is.list(past) || length(past) != B) {
    stop("[expected] past container must have length B=", B)
  }

  # Timeblock vertex attribute for group-to-block mapping
  N  <- network::network.size(nw)
  tb <- network::get.vertex.attribute(nw, "timeblock")
  if (is.null(tb) || length(tb) != N) stop("[expected] missing/invalid vertex attr 'timeblock'.")
  tb_groups <- as.integer(tb[(n1 + 1L):N])
  max_t     <- max(sel_idx)
  time_to_b <- integer(max_t)
  time_to_b[sel_idx] <- seq_len(B)
  group_to_block <- time_to_b[tb_groups]   # length n1, values 1..B

  # Per-block effective size (same truncation logic as InitErgmTerm)
  n_eff <- vapply(seq_len(B), function(b) {
    past_sizes <- vapply(past[[b]][seq_len(d)], length, integer(1))
    min(c(n_b[b], past_sizes))
  }, integer(1))

  sizes_int <- if (is.null(size)) integer(0) else sort(unique(as.integer(size)))
  cnt <- 0L

  for (b in seq_len(B)) {
    if (!is.list(past[[b]]) || length(past[[b]]) < d) {
      stop("[expected] past[[", b, "]] must have at least d=", d, " partitions.")
    }
    off_b  <- actor_offsets[b]   # 0-based actor offset for block b
    neff_b <- n_eff[b]
    G_b    <- n_b[b]              # group vertices per block (padded bipartite)

    for (g in seq_len(G_b)) {
      gv      <- n1 + (b - 1L) * G_b + g
      cur_ids <- .get_group_members_actor_ids(nw, gv, n1)
      # Restrict to actors in [off_b+1 .. off_b+neff_b]
      cur_ids <- sort(cur_ids[cur_ids >= off_b + 1L & cur_ids <= off_b + neff_b])
      if (!length(cur_ids)) next
      if (length(sizes_int) && !(length(cur_ids) %in% sizes_int)) next

      ok_any_lag <- FALSE
      for (lag in seq_len(d)) {
        p_lag  <- past[[b]][[lag]]
        p_eff  <- p_lag[seq_len(neff_b)]    # truncate to n_eff[b]
        groups <- split(seq_along(p_eff), as.integer(p_eff))
        # Global ids: off_b (0-based) + local (1-based)
        groups_global <- lapply(groups, function(v) sort(off_b + as.integer(v)))

        ok_lag <- any(vapply(groups_global, function(gg)
          length(gg) == length(cur_ids) && all(gg == cur_ids), logical(1)))

        if (ok_lag) { ok_any_lag <- TRUE; break }
      }

      if (ok_any_lag) cnt <- cnt + length(cur_ids)   # count actors, not groups
    }
  }

  cnt
}

# ======================================================================================
# Internal entry point lookup (engine)
# ======================================================================================
.f_empile_build <- get0(".erpm_long_empile_build_meta_nw", mode = "function", inherits = TRUE)
if (is.null(.f_empile_build) && "ERPM" %in% loadedNamespaces()) {
  .f_empile_build <- get0(".erpm_long_empile_build_meta_nw", envir = asNamespace("ERPM"),
                          mode = "function", inherits = FALSE)
}
if (is.null(.f_empile_build)) stop("[selftest] Cannot find .erpm_long_empile_build_meta_nw().", call. = FALSE)

# ======================================================================================
# SECTION 1) DRY-RUN (dataset #1)
# ======================================================================================

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
    .st_fail()
    return(invisible(NULL))
  }

  nw <- built$meta_nw
  stopifnot(inherits(nw, "network"))
  n1 <- as.integer(nw %n% "bipartite")
  N  <- network::network.size(nw)

  dbgcat("meta network: N=", N, " | n1(bipartite)=", n1,
         " | selected_partition_indices={", paste(built$selected_partition_indices, collapse = ","), "}")

  # light sanity checks
  tb <- network::get.vertex.attribute(nw, "timeblock")
  if (is.null(tb) || length(tb) != N) stop("[DRY] missing/invalid vertex attr 'timeblock'.", call. = FALSE)

  if (!is.null(nodes_arg)) {
    if (is.null(network::get.vertex.attribute(nw, "gender"))) stop("[DRY] missing vertex attr 'gender'.", call. = FALSE)
    if (is.null(network::get.vertex.attribute(nw, "age")))    stop("[DRY] missing vertex attr 'age'.", call. = FALSE)
  }

  if (!is.null(dyads_arg)) {
    z1 <- tryCatch(nw %n% "Z1", error = function(e) NULL)
    dy <- tryCatch(nw %n% "dyads", error = function(e) NULL)
    okZ1 <- is.matrix(z1) || (is.list(dy) && !is.null(dy[["Z1"]]))
    if (!okZ1) stop("[DRY] dyads filled but cannot find Z1 in nw %n% 'Z1' nor nw %n% 'dyads'[['Z1']].", call. = FALSE)
  }

  if (inert) {
    if (is.null(tryCatch(nw %n% "erpm_block_past_partitions", error = function(e) NULL)))
      stop("[DRY] missing nw %n% 'erpm_block_past_partitions' (required by inertia_groups).", call. = FALSE)
    if (!identical(as.character(tryCatch(nw %n% "erpm_mode", error = function(e) "")), "empile"))
      stop("[DRY] missing/invalid nw %n% 'erpm_mode' for PLE.", call. = FALSE)
  }

  # verify erpm() call returned by erpm_long(eval.call=TRUE)
  formula <- make_formula(partitions, rhs)

  call_erpm <- tryCatch(
    erpm_long(
      formula   = formula,
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
    .st_fail()
    return(invisible(NULL))
  }

  cat("[DRY] OK\n")
  .st_ok()
  invisible(list(meta_nw = nw, call = call_erpm))
}

# Scenario grid (dataset #1)
rhs1  <- quote(cliques(k = 2))
rhs2  <- quote(inertia_groups(past_influence = 1))
rhs3  <- quote(cliques(k = 2) + inertia_groups(past_influence = 1))
rhs4  <- quote(cov_match("gender", clique_size = 2, normalized = "none"))
rhs6  <- quote(cov_match("gender", clique_size = 2, normalized = "none") + inertia_groups(past_influence = 1))
rhs7  <- quote(dyadcov("Z1", clique_size = 2, normalize = FALSE))
rhs9  <- quote(dyadcov("Z1", clique_size = 2, normalize = FALSE) + inertia_groups(past_influence = 1))
rhs_d2 <- quote(inertia_groups(past_influence = 2))

.run_dry_scenario("S1.1)  nodes=NULL;   dyads=NULL;   cliques",                        rhs1,  partitions, NULL,  NULL)
.run_dry_scenario("S1.2)  nodes=NULL;   dyads=NULL;   inertia_groups(d=1)",            rhs2,  partitions, NULL,  NULL)
.run_dry_scenario("S1.3)  nodes=NULL;   dyads=NULL;   cliques + inertia_groups(d=1)",  rhs3,  partitions, NULL,  NULL)
.run_dry_scenario("S1.4)  nodes=filled; dyads=NULL;   cov_match",                      rhs4,  partitions, nodes, NULL)
.run_dry_scenario("S1.5)  nodes=filled; dyads=NULL;   inertia_groups(d=1)",            rhs2,  partitions, nodes, NULL)
.run_dry_scenario("S1.6)  nodes=filled; dyads=NULL;   cov_match + inertia_groups(d=1)",rhs6,  partitions, nodes, NULL)
.run_dry_scenario("S1.7)  nodes=filled; dyads=filled; dyadcov",                        rhs7,  partitions, nodes, dyads)
.run_dry_scenario("S1.8)  nodes=filled; dyads=filled; inertia_groups(d=1)",            rhs2,  partitions, nodes, dyads)
.run_dry_scenario("S1.9)  nodes=filled; dyads=filled; dyadcov + inertia_groups(d=1)",  rhs9,  partitions, nodes, dyads)
.run_dry_scenario("S1.10) nodes=NULL;   dyads=NULL;   inertia_groups(d=2)",            rhs_d2, partitions, NULL, NULL)

cat("\nSECTION 1 DONE.\n")

# ======================================================================================
# SECTION 2) SUMMARY vs OFFLINE COMPUTATION
# ======================================================================================

cat("\n================================================================================\n")
cat("SECTION 2) SUMMARY vs OFFLINE COMPUTATION (datasets #1 and #2)\n")
cat("================================================================================\n")

.run_summary_case <- function(label, partitions, nodes_arg, dyads_arg, rhs_expr,
                             cov_attr = NULL, dyad_name = NULL, past_influence = 1L) {
  cat("\n------------------------------------------------------------\n")
  cat(label, "\n")
  dbgcat("nodes: ", if (is.null(nodes_arg)) "NULL" else "filled",
         " | dyads: ", if (is.null(dyads_arg)) "NULL" else "filled",
         " | rhs: ", deparse(rhs_expr))

  formula <- make_formula(partitions, rhs_expr)

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
    .st_fail()
    return(invisible(NULL))
  }

  nw <- attr(res$value, "meta_nw")
  if (is.null(nw)) {
    cat("[ERPM_LONG ERROR] erpm_long() returned no 'meta_nw' attribute.\n")
    .st_fail()
    return(invisible(NULL))
  }
  stopifnot(inherits(nw, "network"))

  n1 <- as.integer(nw %n% "bipartite")
  N  <- network::network.size(nw)
  dbgcat("built nw: N=", N, " | n1=", n1)

  # summary() on that network, explicitly in a formula env that contains `nw`
  sum_formula <- make_summary_formula(nw, rhs_expr)

  .warn_flush()
  sres <- .capture_warnings(
    tryCatch(
      summary(sum_formula, constraints = ~ b1partblockdiag("timeblock")),
      error = function(e) e
    )
  )
  .print_warnings(sres$warnings)
  if (inherits(sres$value, "error")) {
    cat("[SUMMARY ERROR]\n", conditionMessage(sres$value), "\n")
    .st_fail()
    return(invisible(NULL))
  }

  s <- sres$value
  print(s)

  # offline expected vector in the same order as rhs_expr terms are written
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
    cat("[MISMATCH] expected length ", length(o), " != summary length ", length(sval), " (check rhs ordering).\n")
    .st_fail()
    return(invisible(NULL))
  }
  if (!isTRUE(all.equal(sval, o))) {
    cat("[MISMATCH] summary != expected in: ", label, "\n")
    .st_fail()
    return(invisible(NULL))
  }

  cat("[SUMMARY vs OFFLINE] OK\n")
  .st_ok()
  invisible(list(nw = nw, summary = s))
}

# ---------------------------
# Dataset #1: RHS blocks
# ---------------------------
rhs_1A  <- quote(cliques(k = 2))
rhs_1B  <- quote(inertia_groups(past_influence = 1))
rhs_1C  <- quote(cliques(k = 2) + inertia_groups(past_influence = 1))
rhs_1D  <- quote(cov_match("gender", clique_size = 2, normalized = "none"))
rhs_1F  <- quote(cov_match("gender", clique_size = 2, normalized = "none") + inertia_groups(past_influence = 1))
rhs_1G  <- quote(dyadcov("Z1", clique_size = 2, normalize = FALSE))
rhs_1I  <- quote(dyadcov("Z1", clique_size = 2, normalize = FALSE) + inertia_groups(past_influence = 1))
rhs_1J  <- quote(
  cliques(k = 2) +
    cov_match("gender", clique_size = 2, normalized = "none") +
    dyadcov("Z1", clique_size = 2, normalize = FALSE) +
    inertia_groups(past_influence = 1)
)
rhs_1K  <- quote(inertia_groups(past_influence = 2))

cat("\n--- DATASET #1 (n=4) ---\n")
.run_summary_case("S2.1A) nodes=NULL;   dyads=NULL;   cliques",
                  partitions, NULL, NULL, rhs_1A)

.run_summary_case("S2.1B) nodes=NULL;   dyads=NULL;   inertia_groups(d=1)",
                  partitions, NULL, NULL, rhs_1B, past_influence = 1L)

.run_summary_case("S2.1C) nodes=NULL;   dyads=NULL;   cliques + inertia_groups(d=1)",
                  partitions, NULL, NULL, rhs_1C, past_influence = 1L)

.run_summary_case("S2.1D) nodes=filled; dyads=NULL;   cov_match(gender)",
                  partitions, nodes, NULL, rhs_1D, cov_attr = "gender")

.run_summary_case("S2.1E) nodes=filled; dyads=NULL;   inertia_groups(d=1)",
                  partitions, nodes, NULL, rhs_1B, past_influence = 1L)

.run_summary_case("S2.1F) nodes=filled; dyads=NULL;   cov_match + inertia_groups(d=1)",
                  partitions, nodes, NULL, rhs_1F, cov_attr = "gender", past_influence = 1L)

.run_summary_case("S2.1G) nodes=filled; dyads=filled; dyadcov(Z1)",
                  partitions, nodes, dyads, rhs_1G, dyad_name = "Z1")

.run_summary_case("S2.1H) nodes=filled; dyads=filled; inertia_groups(d=1)",
                  partitions, nodes, dyads, rhs_1B, past_influence = 1L)

.run_summary_case("S2.1I) nodes=filled; dyads=filled; dyadcov + inertia_groups(d=1)",
                  partitions, nodes, dyads, rhs_1I, dyad_name = "Z1", past_influence = 1L)

.run_summary_case("S2.1J) nodes=filled; dyads=filled; ALL (d=1)",
                  partitions, nodes, dyads, rhs_1J, cov_attr = "gender", dyad_name = "Z1", past_influence = 1L)

.run_summary_case("S2.1K) nodes=NULL;   dyads=NULL;   inertia_groups(d=2)",
                  partitions, NULL, NULL, rhs_1K, past_influence = 2L)

# ---------------------------
# Dataset #2: RHS blocks
#   - cov_match uses "sex", dyadcov uses "X1"
# ---------------------------
rhs_2A  <- quote(cliques(k = 2))
rhs_2B  <- quote(inertia_groups(past_influence = 1))
rhs_2C  <- quote(cliques(k = 2) + inertia_groups(past_influence = 1))
rhs_2D  <- quote(cov_match("sex", clique_size = 2, normalized = "none"))
rhs_2F  <- quote(cov_match("sex", clique_size = 2, normalized = "none") + inertia_groups(past_influence = 1))
rhs_2G  <- quote(dyadcov("X1", clique_size = 2, normalize = FALSE))
rhs_2I  <- quote(dyadcov("X1", clique_size = 2, normalize = FALSE) + inertia_groups(past_influence = 1))
rhs_2J  <- quote(
  cliques(k = 2) +
    cov_match("sex", clique_size = 2, normalized = "none") +
    dyadcov("X1", clique_size = 2, normalize = FALSE) +
    inertia_groups(past_influence = 1)
)
rhs_2K  <- quote(inertia_groups(past_influence = 2))

cat("\n--- DATASET #2 (n=5) ---\n")
.run_summary_case("S2.2A) nodes=NULL;   dyads=NULL;   cliques",
                  partitions2, NULL, NULL, rhs_2A)

.run_summary_case("S2.2B) nodes=NULL;   dyads=NULL;   inertia_groups(d=1)",
                  partitions2, NULL, NULL, rhs_2B, past_influence = 1L)

.run_summary_case("S2.2C) nodes=NULL;   dyads=NULL;   cliques + inertia_groups(d=1)",
                  partitions2, NULL, NULL, rhs_2C, past_influence = 1L)

.run_summary_case("S2.2D) nodes=filled; dyads=NULL;   cov_match(sex)",
                  partitions2, nodes2, NULL, rhs_2D, cov_attr = "sex")

.run_summary_case("S2.2E) nodes=filled; dyads=NULL;   inertia_groups(d=1)",
                  partitions2, nodes2, NULL, rhs_2B, past_influence = 1L)

.run_summary_case("S2.2F) nodes=filled; dyads=NULL;   cov_match + inertia_groups(d=1)",
                  partitions2, nodes2, NULL, rhs_2F, cov_attr = "sex", past_influence = 1L)

.run_summary_case("S2.2G) nodes=filled; dyads=filled; dyadcov(X1)",
                  partitions2, nodes2, dyads2, rhs_2G, dyad_name = "X1")

.run_summary_case("S2.2H) nodes=filled; dyads=filled; inertia_groups(d=1)",
                  partitions2, nodes2, dyads2, rhs_2B, past_influence = 1L)

.run_summary_case("S2.2I) nodes=filled; dyads=filled; dyadcov + inertia_groups(d=1)",
                  partitions2, nodes2, dyads2, rhs_2I, dyad_name = "X1", past_influence = 1L)

.run_summary_case("S2.2J) nodes=filled; dyads=filled; ALL (d=1)",
                  partitions2, nodes2, dyads2, rhs_2J, cov_attr = "sex", dyad_name = "X1", past_influence = 1L)

.run_summary_case("S2.2K) nodes=NULL;   dyads=NULL;   inertia_groups(d=2)",
                  partitions2, NULL, NULL, rhs_2K, past_influence = 2L)

cat("\nSECTION 2 DONE.\n")

# ======================================================================================
# SECTION 3) FITS: erpm_long(eval.call=TRUE) + eval()
#   - Requires the ergm patch for ergm 4.10.x (replace bug).
#   - Some cases may fail due to near-constant statistics or separation on small data;
#     errors are logged and the selftest continues.
#   - On success, the coefficient estimate and convergence are reported.
# ======================================================================================

cat("\n================================================================================\n")
cat("SECTION 3) FITS: eval(erpm_long(..., eval.call=TRUE)) (datasets #1 and #2)\n")
cat("================================================================================\n")

if (!.patch_enabled) {
  cat("[SECTION 3 SKIPPED] ergm patch not available — fits require the patch.\n")
} else {

.run_fit_case <- function(label, partitions, nodes_arg, dyads_arg, rhs_expr) {
  cat("\n------------------------------------------------------------\n")
  cat(label, "\n")
  dbgcat("nodes: ", if (is.null(nodes_arg)) "NULL" else "filled",
         " | dyads: ", if (is.null(dyads_arg)) "NULL" else "filled",
         " | rhs: ", deparse(rhs_expr))

  formula <- make_formula(partitions, rhs_expr)

  # Step 1: build meta-network + get erpm() call
  .warn_flush()
  erpm_call <- .capture_warnings(
    tryCatch(
      erpm_long(
        formula      = formula,
        nodes        = nodes_arg,
        dyads        = dyads_arg,
        mode         = "empile",
        verbose      = FALSE,
        debug        = FALSE,
        eval.call    = TRUE,
        eval.loglik  = FALSE,
        constraints  = ~ b1partblockdiag("timeblock")
      ),
      error = function(e) e
    )
  )
  .print_warnings(erpm_call$warnings)

  if (inherits(erpm_call$value, "error")) {
    cat("[FIT ERROR (build)]\n", conditionMessage(erpm_call$value), "\n")
    .st_fail()
    return(invisible(NULL))
  }

  # Step 2: run the erpm() call, which runs ergm()
  fit_res <- .capture_warnings(
    tryCatch(eval(erpm_call$value), error = function(e) e)
  )
  .print_warnings(fit_res$warnings)

  if (inherits(fit_res$value, "error")) {
    msg <- conditionMessage(fit_res$value)
    if (grepl("essentially constant|degenerate|separation", msg, ignore.case = TRUE)) {
      cat("[FIT SKIP] degenerate data (constant statistic) —", label, "\n")
      return(invisible(NULL))
    }
    cat("[FIT ERROR (ergm)]\n", msg, "\n")
    .st_fail()
    return(invisible(NULL))
  }

  fit <- fit_res$value
  co  <- tryCatch(round(coef(fit), 4), error = function(e) NULL)
  cat("[FIT] OK | coef:", if (is.null(co)) "?" else paste(co, collapse = ", "), "\n")
  .st_ok()
  invisible(fit)
}

cat("\n--- DATASET #1 FITS (n=4) ---\n")
.run_fit_case("S3.1A) nodes=NULL;   dyads=NULL;   cliques",                       partitions,  NULL,  NULL,  rhs_1A)
.run_fit_case("S3.1B) nodes=NULL;   dyads=NULL;   inertia_groups(d=1)",           partitions,  NULL,  NULL,  rhs_1B)
.run_fit_case("S3.1C) nodes=NULL;   dyads=NULL;   cliques + inertia_groups(d=1)", partitions,  NULL,  NULL,  rhs_1C)
.run_fit_case("S3.1D) nodes=filled; dyads=NULL;   cov_match",                     partitions,  nodes, NULL,  rhs_1D)
.run_fit_case("S3.1E) nodes=filled; dyads=NULL;   inertia_groups(d=1)",           partitions,  nodes, NULL,  rhs_1B)
.run_fit_case("S3.1F) nodes=filled; dyads=NULL;   cov_match + inertia_groups",    partitions,  nodes, NULL,  rhs_1F)
.run_fit_case("S3.1G) nodes=filled; dyads=filled; dyadcov",                       partitions,  nodes, dyads, rhs_1G)
.run_fit_case("S3.1H) nodes=filled; dyads=filled; inertia_groups(d=1)",           partitions,  nodes, dyads, rhs_1B)
.run_fit_case("S3.1I) nodes=filled; dyads=filled; dyadcov + inertia_groups",      partitions,  nodes, dyads, rhs_1I)
.run_fit_case("S3.1J) nodes=filled; dyads=filled; ALL (d=1)",                     partitions,  nodes, dyads, rhs_1J)
.run_fit_case("S3.1K) nodes=NULL;   dyads=NULL;   inertia_groups(d=2)",           partitions,  NULL,  NULL,  rhs_1K)

cat("\n--- DATASET #2 FITS (n=5) ---\n")
.run_fit_case("S3.2A) nodes=NULL;   dyads=NULL;   cliques",                       partitions2, NULL,   NULL,   rhs_2A)
.run_fit_case("S3.2B) nodes=NULL;   dyads=NULL;   inertia_groups(d=1)",           partitions2, NULL,   NULL,   rhs_2B)
.run_fit_case("S3.2C) nodes=NULL;   dyads=NULL;   cliques + inertia_groups(d=1)", partitions2, NULL,   NULL,   rhs_2C)
.run_fit_case("S3.2D) nodes=filled; dyads=NULL;   cov_match",                     partitions2, nodes2, NULL,   rhs_2D)
.run_fit_case("S3.2E) nodes=filled; dyads=NULL;   inertia_groups(d=1)",           partitions2, nodes2, NULL,   rhs_2B)
.run_fit_case("S3.2F) nodes=filled; dyads=NULL;   cov_match + inertia_groups",    partitions2, nodes2, NULL,   rhs_2F)
.run_fit_case("S3.2G) nodes=filled; dyads=filled; dyadcov",                       partitions2, nodes2, dyads2, rhs_2G)
.run_fit_case("S3.2H) nodes=filled; dyads=filled; inertia_groups(d=1)",           partitions2, nodes2, dyads2, rhs_2B)
.run_fit_case("S3.2I) nodes=filled; dyads=filled; dyadcov + inertia_groups",      partitions2, nodes2, dyads2, rhs_2I)
.run_fit_case("S3.2J) nodes=filled; dyads=filled; ALL (d=1)",                     partitions2, nodes2, dyads2, rhs_2J)
.run_fit_case("S3.2K) nodes=NULL;   dyads=NULL;   inertia_groups(d=2)",           partitions2, NULL,   NULL,   rhs_2K)

} # end if (.patch_enabled)

cat("\nSECTION 3 DONE.\n")

# ======================================================================================
# SECTION 4) b1partblockdiag constraint regression
#   - Verifies constraint structure (dependence, implies, free_dyads class)
#   - Verifies summary(~b1partblockdiag) == summary(~b1part) on PLE meta-networks
#   - Regression test: formula env isolation (partitions2 not confused with global partitions)
# ======================================================================================

cat("\n================================================================================\n")
cat("SECTION 4) b1partblockdiag constraint regression\n")
cat("================================================================================\n")

.run_b1bd_case <- function(label, partitions, nodes_arg, dyads_arg, rhs_expr,
                           past_influence = 1L) {
  cat("\n------------------------------------------------------------\n")
  cat(label, "\n")

  formula <- make_formula(partitions, rhs_expr)

  # Build meta-network
  erpm_res <- tryCatch(
    erpm_long(formula = formula, nodes = nodes_arg, dyads = dyads_arg,
              mode = "empile", verbose = FALSE, debug = FALSE, eval.call = FALSE),
    error = function(e) e
  )
  if (inherits(erpm_res, "error")) {
    cat("[B1BD ERROR (build)]\n", conditionMessage(erpm_res), "\n")
    .st_fail(); return(invisible(NULL))
  }
  nw <- attr(erpm_res, "meta_nw")
  if (is.null(nw)) { cat("[B1BD ERROR] no meta_nw\n"); .st_fail(); return(invisible(NULL)) }

  # Test 1: constraint fields
  con <- tryCatch(
    InitErgmConstraint.b1partblockdiag(nw, list(attr = "timeblock")),
    error = function(e) e
  )
  if (inherits(con, "error")) {
    cat("[B1BD ERROR (constraint init)]\n", conditionMessage(con), "\n")
    .st_fail(); return(invisible(NULL))
  }
  if (!isTRUE(con$dependence)) {
    cat("[B1BD FAIL] dependence should be TRUE\n"); .st_fail(); return(invisible(NULL))
  }
  if (!identical(con$implies, c("b1degrees", "edges"))) {
    cat("[B1BD FAIL] implies mismatch\n"); .st_fail(); return(invisible(NULL))
  }
  if (!inherits(con$free_dyads, "rlebdm")) {
    cat("[B1BD FAIL] free_dyads is not rlebdm\n"); .st_fail(); return(invisible(NULL))
  }

  # Test 2: summary(~b1part) == summary(~b1partblockdiag("timeblock"))
  sum_formula <- make_summary_formula(nw, rhs_expr)
  s_b1 <- tryCatch(
    as.numeric(summary(sum_formula, constraints = ~ b1part)),
    error = function(e) e
  )
  s_bd <- tryCatch(
    as.numeric(summary(sum_formula, constraints = ~ b1partblockdiag("timeblock"))),
    error = function(e) e
  )
  if (inherits(s_b1, "error")) {
    cat("[B1BD ERROR (summary b1part)]\n", conditionMessage(s_b1), "\n")
    .st_fail(); return(invisible(NULL))
  }
  if (inherits(s_bd, "error")) {
    cat("[B1BD ERROR (summary b1partblockdiag)]\n", conditionMessage(s_bd), "\n")
    .st_fail(); return(invisible(NULL))
  }
  if (!isTRUE(all.equal(s_b1, s_bd))) {
    cat("[B1BD FAIL] summary(~b1part) != summary(~b1partblockdiag)\n")
    cat("  b1part :", s_b1, "\n")
    cat("  b1bd   :", s_bd, "\n")
    .st_fail(); return(invisible(NULL))
  }

  cat("[B1BD] OK | summary equal for both constraints\n")
  .st_ok()
  invisible(list(nw = nw, s_b1 = s_b1, s_bd = s_bd))
}

cat("\n--- Dataset #1 ---\n")
.run_b1bd_case("S4.1A) cliques",                   partitions,  NULL,  NULL,  rhs_1A)
.run_b1bd_case("S4.1B) inertia_groups(d=1)",        partitions,  NULL,  NULL,  rhs_1B, past_influence = 1L)
.run_b1bd_case("S4.1C) cliques + inertia_groups",   partitions,  NULL,  NULL,  rhs_1C, past_influence = 1L)
.run_b1bd_case("S4.1K) inertia_groups(d=2)",        partitions,  NULL,  NULL,  rhs_1K, past_influence = 2L)

cat("\n--- Dataset #2 ---\n")
.run_b1bd_case("S4.2A) cliques",                    partitions2, NULL,  NULL,  rhs_2A)
.run_b1bd_case("S4.2B) inertia_groups(d=1)",         partitions2, NULL,  NULL,  rhs_2B, past_influence = 1L)
.run_b1bd_case("S4.2K) inertia_groups(d=2)",         partitions2, NULL,  NULL,  rhs_2K, past_influence = 2L)

# Regression test: formula env isolation
# Verify that passing partitions2 (n=5) works correctly even when a global
# variable named 'partitions' (n=4) exists — regression for eval(lhs) scoping bug.
cat("\n--- Regression: formula env isolation ---\n")
.reg_formula <- make_formula(partitions2, rhs_2B)
.reg_res <- tryCatch(
  erpm_long(formula = .reg_formula, mode = "empile", verbose = FALSE,
            eval.call = FALSE, eval.loglik = FALSE),
  error = function(e) e
)
if (inherits(.reg_res, "error")) {
  cat("[REGRESSION FAIL] formula env isolation: ", conditionMessage(.reg_res), "\n")
  .st_fail()
} else {
  .reg_nw <- attr(.reg_res, "meta_nw")
  .reg_n1 <- as.integer(.reg_nw %n% "bipartite")
  if (.reg_n1 != 5L * 2L) {
    cat("[REGRESSION FAIL] expected n1=10 (5 actors * 2 selected blocks), got ", .reg_n1, "\n")
    .st_fail()
  } else {
    cat("[REGRESSION] formula env isolation OK (n1=", .reg_n1, ")\n", sep = "")
    .st_ok()
  }
}
rm(.reg_formula, .reg_res, .reg_nw, .reg_n1)

cat("\nSECTION 4 DONE.\n")

# ======================================================================================
# FINAL SUMMARY
# ======================================================================================

cat("\n================================================================================\n")
cat(sprintf("SELFTEST DONE — OK: %d | FAIL: %d | TOTAL: %d\n",
            .st$ok, .st$fail, .st$ok + .st$fail))
cat("================================================================================\n")

if (.st$fail > 0L) {
  message(sprintf("[SELFTEST] %d case(s) FAILED.", .st$fail))
}

if (.patch_enabled && exists("ergm_patch_disable", mode = "function")) {
  try(ergm_patch_disable(), silent = TRUE)
}