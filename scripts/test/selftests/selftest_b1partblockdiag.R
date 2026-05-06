# ==============================================================================
# File   : scripts/test/selftests/selftest_b1partblockdiag.R
# Object : Self-test for ERPM constraint `b1partblockdiag`
# Run    : Rscript -e "devtools::load_all('.', quiet=TRUE); source('scripts/test/selftests/selftest_b1partblockdiag.R')"
#
# SECTION 1) Constraint fields: dependence, implies, free_dyads class, attr
# SECTION 2) free_dyads structure × 4 configs (count + intra/inter-block)
# SECTION 3) Error conditions (non-contiguous blocks, non-bipartite)
# SECTION 4) Proposal B1PartBlockdiag (name + non-bipartite error)
# SECTION 5) PLE coherence via erpm_long(eval.call=FALSE)
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  if (!requireNamespace("network",  quietly = TRUE)) stop("Package 'network' required.")
  if (!requireNamespace("ergm",     quietly = TRUE)) stop("Package 'ergm' required.")
})
suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

# ---- counters ----------------------------------------------------------------
.st <- new.env(parent = emptyenv())
.st$ok   <- 0L
.st$fail <- 0L
.st_ok   <- function(label) { cat(sprintf("[OK]   %s\n", label)); .st$ok   <- .st$ok   + 1L; invisible(NULL) }
.st_fail <- function(label, reason = "") {
  cat(sprintf("[FAIL] %s%s\n", label, if (nchar(reason)) paste0(" — ", reason) else ""))
  .st$fail <- .st$fail + 1L; invisible(NULL)
}

cat("=== SELFTEST b1partblockdiag ===\n")

# ---- helpers -----------------------------------------------------------------

# Build a bipartite network with contiguous blocks.
# n_actors / n_groups: integer vectors of per-block sizes (same length).
.make_bip <- function(n_actors, n_groups, block_labels = seq_along(n_actors)) {
  n1 <- sum(n_actors)
  n2 <- sum(n_groups)
  nw <- network::network(n1 + n2, bipartite = n1, directed = FALSE)
  tb_actors <- rep(block_labels, times = n_actors)
  tb_groups <- rep(block_labels, times = n_groups)
  network::set.vertex.attribute(nw, "timeblock", c(tb_actors, tb_groups))
  nw
}

# Retrieve the constraint for a network.
.get_con <- function(nw) {
  InitErgmConstraint.b1partblockdiag(nw, list(attr = "timeblock"))
}

# Check free_dyads count and intra/inter-block structure.
.check_fd <- function(nw, n_actors, n_groups, label) {
  con <- tryCatch(.get_con(nw), error = function(e) e)
  if (inherits(con, "error")) {
    .st_fail(label, conditionMessage(con)); return(invisible(FALSE))
  }
  mat <- as.matrix(con$free_dyads)
  n1  <- sum(n_actors)
  n2  <- sum(n_groups)

  # Count check
  expected <- 2L * sum(as.integer(n_actors) * as.integer(n_groups))
  actual   <- as.integer(sum(mat))
  if (actual != expected) {
    .st_fail(label, sprintf("free_dyads count=%d expected=%d", actual, expected))
    return(invisible(FALSE))
  }

  # Intra-block TRUE, inter-block FALSE
  tb_a <- rep(seq_along(n_actors), times = n_actors)
  tb_g <- rep(seq_along(n_groups), times = n_groups)
  ok <- TRUE
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      exp_ij <- (tb_a[i] == tb_g[j])
      if (!identical(as.logical(mat[i, n1 + j]),   exp_ij) ||
          !identical(as.logical(mat[n1 + j, i]), exp_ij)) {
        .st_fail(label, sprintf("mat[%d,%d]=%s expected=%s",
                                i, n1 + j, mat[i, n1 + j], exp_ij))
        ok <- FALSE; break
      }
    }
    if (!ok) break
  }
  if (!ok) return(invisible(FALSE))

  .st_ok(label); invisible(TRUE)
}

# ==============================================================================
# SECTION 1 — Constraint fields
# ==============================================================================
cat("\n================================================================================\n")
cat("SECTION 1) Constraint fields\n")
cat("================================================================================\n")

{
  label <- "S1) fields: dependence, implies, free_dyads class, attr"
  nw <- .make_bip(c(3L, 2L), c(3L, 2L))
  con <- tryCatch(.get_con(nw), error = function(e) e)
  if (inherits(con, "error")) {
    .st_fail(label, conditionMessage(con))
  } else if (!isTRUE(con$dependence)) {
    .st_fail(label, "dependence is not TRUE")
  } else if (!identical(con$implies, c("b1degrees", "edges"))) {
    .st_fail(label, sprintf("implies = %s", paste(con$implies, collapse = ",")))
  } else if (!inherits(con$free_dyads, "rlebdm")) {
    .st_fail(label, sprintf("free_dyads class = %s", paste(class(con$free_dyads), collapse = "/")))
  } else if (!identical(con$attr, "timeblock")) {
    .st_fail(label, sprintf("attr = %s", con$attr))
  } else {
    .st_ok(label)
  }
}

# ==============================================================================
# SECTION 2 — free_dyads structure
# ==============================================================================
cat("\n================================================================================\n")
cat("SECTION 2) free_dyads structure\n")
cat("================================================================================\n")

# C1: B=1, n_actors=[4], n_groups=[4]  — free_dyads = 2*16 = 32
.check_fd(.make_bip(4L, 4L), 4L, 4L, "S2.C1) B=1 n=[4] | expected=32")

# C2: B=2 equal, n_actors=[3,3], n_groups=[3,3]  — free_dyads = 2*(9+9) = 36
.check_fd(.make_bip(c(3L, 3L), c(3L, 3L)), c(3L, 3L), c(3L, 3L),
          "S2.C2) B=2 equal n=[3,3] | expected=36")

# C3: B=2 unequal, n_actors=[2,4], n_groups=[2,4]  — free_dyads = 2*(4+16) = 40
.check_fd(.make_bip(c(2L, 4L), c(2L, 4L)), c(2L, 4L), c(2L, 4L),
          "S2.C3) B=2 unequal n=[2,4] | expected=40")

# C4: B=3, n_actors=[2,3,2], n_groups=[2,3,2]  — free_dyads = 2*(4+9+4) = 34
.check_fd(.make_bip(c(2L, 3L, 2L), c(2L, 3L, 2L)), c(2L, 3L, 2L), c(2L, 3L, 2L),
          "S2.C4) B=3 n=[2,3,2] | expected=34")

# ==============================================================================
# SECTION 3 — Error conditions
# ==============================================================================
cat("\n================================================================================\n")
cat("SECTION 3) Error conditions\n")
cat("================================================================================\n")

# Non-contiguous mode-1 (actors): timeblock = [1,2,1,2, 1,2,1,2]
{
  label <- "S3.1) non-contiguous mode-1 -> stop"
  nw <- network::network(8L, bipartite = 4L, directed = FALSE)
  network::set.vertex.attribute(nw, "timeblock", c(1L, 2L, 1L, 2L,  1L, 2L, 1L, 2L))
  res <- tryCatch(.get_con(nw), error = function(e) e)
  if (inherits(res, "error") && grepl("contig", conditionMessage(res), ignore.case = TRUE)) {
    .st_ok(label)
  } else if (inherits(res, "error")) {
    .st_ok(paste(label, "(error caught, different msg)"))
  } else {
    .st_fail(label, "no error raised for non-contiguous mode-1")
  }
}

# Non-contiguous mode-2 (groups): actors contiguous, groups non-contiguous
{
  label <- "S3.2) non-contiguous mode-2 -> stop"
  nw <- network::network(8L, bipartite = 4L, directed = FALSE)
  network::set.vertex.attribute(nw, "timeblock", c(1L, 1L, 2L, 2L,  1L, 2L, 1L, 2L))
  res <- tryCatch(.get_con(nw), error = function(e) e)
  if (inherits(res, "error") && grepl("contig", conditionMessage(res), ignore.case = TRUE)) {
    .st_ok(label)
  } else if (inherits(res, "error")) {
    .st_ok(paste(label, "(error caught, different msg)"))
  } else {
    .st_fail(label, "no error raised for non-contiguous mode-2")
  }
}

# Non-bipartite network
{
  label <- "S3.3) non-bipartite network -> stop"
  nw_dir <- network::network(6L, directed = FALSE)
  res <- tryCatch(
    InitErgmConstraint.b1partblockdiag(nw_dir, list(attr = "timeblock")),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    .st_ok(label)
  } else {
    .st_fail(label, "no error raised for non-bipartite network")
  }
}

# ==============================================================================
# SECTION 4 — Proposal B1PartBlockdiag
# ==============================================================================
cat("\n================================================================================\n")
cat("SECTION 4) Proposal B1PartBlockdiag\n")
cat("================================================================================\n")

# Normal case: returns name="B1Part"
{
  label <- "S4.1) InitErgmProposal.B1PartBlockdiag returns name='B1Part'"
  nw <- .make_bip(c(3L, 2L), c(3L, 2L))
  res <- tryCatch(InitErgmProposal.B1PartBlockdiag(NULL, nw), error = function(e) e)
  if (inherits(res, "error")) {
    .st_fail(label, conditionMessage(res))
  } else if (!identical(res$name, "B1Part")) {
    .st_fail(label, sprintf("name = '%s'", res$name))
  } else {
    .st_ok(label)
  }
}

# Non-bipartite: should error
{
  label <- "S4.2) InitErgmProposal.B1PartBlockdiag non-bipartite -> error"
  nw_dir <- network::network(6L, directed = FALSE)
  res <- tryCatch(InitErgmProposal.B1PartBlockdiag(NULL, nw_dir), error = function(e) e)
  if (inherits(res, "error")) {
    .st_ok(label)
  } else {
    .st_fail(label, "no error raised for non-bipartite")
  }
}

# ==============================================================================
# SECTION 5 — PLE coherence via erpm_long(eval.call=FALSE)
# ==============================================================================
cat("\n================================================================================\n")
cat("SECTION 5) PLE coherence\n")
cat("================================================================================\n")

# Dataset A: n=4, T=3, d=1  ->  B=2 blocks of 4+4, free_dyads = 2*(4*4+4*4) = 64
{
  label <- "S5.A) n=4 T=3 d=1 | B=2 blocks, free_dyads=64"
  parts_a <- list(c(1L, 1L, 2L, 3L), c(1L, 2L, 2L, 3L), c(1L, 1L, 3L, 3L))
  res <- tryCatch(
    erpm_long(parts_a ~ inertia_groups(past_influence = 1L),
              eval.call = FALSE, verbose = FALSE),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    .st_fail(label, conditionMessage(res))
  } else {
    nw_a <- attr(res, "meta_nw")
    con_a <- tryCatch(.get_con(nw_a), error = function(e) e)
    if (inherits(con_a, "error")) {
      .st_fail(label, conditionMessage(con_a))
    } else {
      actual <- as.integer(sum(as.matrix(con_a$free_dyads)))
      if (actual == 64L) {
        .st_ok(label)
      } else {
        .st_fail(label, sprintf("free_dyads count=%d expected=64", actual))
      }
    }
  }
}

# Dataset B: n=5, T=3, d=1  ->  B=2 blocks of 5+5, free_dyads = 2*(5*5+5*5) = 100
{
  label <- "S5.B) n=5 T=3 d=1 | B=2 blocks, free_dyads=100"
  parts_b <- list(c(1L, 1L, 2L, 3L, 3L), c(1L, 2L, 2L, 3L, 1L), c(2L, 2L, 3L, 3L, 1L))
  res <- tryCatch(
    erpm_long(parts_b ~ inertia_groups(past_influence = 1L),
              eval.call = FALSE, verbose = FALSE),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    .st_fail(label, conditionMessage(res))
  } else {
    nw_b <- attr(res, "meta_nw")
    con_b <- tryCatch(.get_con(nw_b), error = function(e) e)
    if (inherits(con_b, "error")) {
      .st_fail(label, conditionMessage(con_b))
    } else {
      actual <- as.integer(sum(as.matrix(con_b$free_dyads)))
      if (actual == 100L) {
        .st_ok(label)
      } else {
        .st_fail(label, sprintf("free_dyads count=%d expected=100", actual))
      }
    }
  }
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
cat("\n================================================================================\n")
cat(sprintf("SELFTEST DONE — OK: %d | FAIL: %d | TOTAL: %d\n",
            .st$ok, .st$fail, .st$ok + .st$fail))
cat("================================================================================\n")

if (.st$fail > 0L) {
  message(sprintf("[SELFTEST] %d case(s) FAILED.", .st$fail))
}
