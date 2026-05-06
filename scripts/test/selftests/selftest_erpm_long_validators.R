# ==============================================================================
# File    : scripts/test/selftests/selftest_erpm_long_validators.R
# Auteur  : Jérémie Chichignoud - Cub'itech
# Purpose : Self-test for erpm_long() input validators in PLE-only mode.
#
# Notes
#   - This script is meant to be run from the package root (DESCRIPTION present).
#   - It targets validator behavior only, using deliberately invalid inputs and
#     checking that erpm_long() fails with the expected error messages.
#   - It does not test the PLE engine construction, does not validate summary()
#     outputs, and does not run model fits (the wrapper is called with eval.call=FALSE).
#   - Two small datasets are used (n=4 and n=5, T=3) to cover basic shape variation.
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) stop("Package 'devtools' requis.")
  if (!requireNamespace("network",  quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",     quietly = TRUE)) stop("Package 'ergm' requis.")
})

# Load ERPM in dev mode (adjust if you prefer library(ERPM))
options(keep.source = TRUE)
options(keep.source.pkgs = TRUE)
Sys.setenv(R_KEEP_PKG_SOURCE = "yes")
# devtools::load_all(".")

cat("=== SELFTEST erpm_long validators (PLE only) ===\n")

# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------

.capture_error <- function(expr) {
  tryCatch(
    list(ok = TRUE,  err = NULL, val = eval.parent(substitute(expr))),
    error = function(e) list(ok = FALSE, err = conditionMessage(e), val = NULL)
  )
}

.expect_fail <- function(label, expr, pattern = NULL) {
  cat("\n------------------------------------------------------------\n")
  cat("[EXPECT FAIL] ", label, "\n", sep = "")
  r <- .capture_error(expr)
  if (r$ok) {
    cat("!! UNEXPECTED PASS\n")
    return(invisible(FALSE))
  }
  cat("got error: ", r$err, "\n", sep = "")
  if (!is.null(pattern) && !grepl(pattern, r$err)) {
    cat("!! ERROR PATTERN MISMATCH\n")
    cat("   expected pattern: ", pattern, "\n", sep = "")
    return(invisible(FALSE))
  }
  cat("OK\n")
  invisible(TRUE)
}

.expect_pass <- function(label, expr) {
  cat("\n------------------------------------------------------------\n")
  cat("[EXPECT PASS] ", label, "\n", sep = "")
  r <- .capture_error(expr)
  if (!r$ok) {
    cat("!! UNEXPECTED FAIL\n")
    cat("got error: ", r$err, "\n", sep = "")
    return(invisible(FALSE))
  }
  cat("OK\n")
  invisible(TRUE)
}

# We do NOT want to run engine or ergm. So we force eval.call=FALSE always.
.call_erpm_long <- function(formula, mode = "empile", nodes = NULL, dyads = NULL,
                            verbose = FALSE, debug = FALSE, seed = NULL, group_labels = NULL) {
  erpm_long(
    formula   = formula,
    mode      = mode,
    nodes     = nodes,
    dyads     = dyads,
    verbose   = verbose,
    debug     = debug,
    seed      = seed,
    group_labels = group_labels,
    eval.call = FALSE
  )
}

# --------------------------------------------------------------------------------------
# Baseline datasets (good)
# --------------------------------------------------------------------------------------
partitions1 <- list(
  c(1, 1, 2, 3),
  c(1, 2, 2, 3),
  c(1, 1, 3, 3)
)

nodes1 <- list(
  data.frame(label = c("A","B","C","D"), gender = c(1,1,2,1), age = c(20,22,25,30)),
  data.frame(label = c("FT","AZ","JI","DO"), gender = c(2,1,2,2), age = c(10,42,25,30)),
  data.frame(label = c("H","Z","S","A"), gender = c(1,1,1,1), age = c(27,26,25,28))
)

dyads1 <- list(
  list(
    fm = matrix(c(0,1,1,0, 1,0,1,0, 1,1,0,0, 0,0,0,0), nrow = 4, byrow = TRUE),
    Z1 = matrix(c(0,2,3,0, 2,0,4,0, 3,4,0,1, 0,0,1,0), nrow = 4, byrow = TRUE)
  ),
  list(
    fm = matrix(c(0,1,0,1, 1,0,1,0, 0,1,0,1, 1,0,1,0), nrow = 4, byrow = TRUE),
    Z1 = matrix(c(0,5,0,2, 5,0,1,0, 0,1,0,3, 2,0,3,0), nrow = 4, byrow = TRUE)
  ),
  list(
    fm = matrix(c(0,0,1,1, 0,0,1,1, 1,1,0,0, 1,1,0,0), nrow = 4, byrow = TRUE),
    Z1 = matrix(c(0,1,2,3, 1,0,4,5, 2,4,0,6, 3,5,6,0), nrow = 4, byrow = TRUE)
  )
)

partitions2 <- list(
  c(1, 1, 2, 3, 3),
  c(1, 2, 2, 3, 1),
  c(2, 2, 3, 3, 1)
)

nodes2 <- list(
  data.frame(id = c("E","F","G","H","I"), sex = c(1,2,1,2,1), age_years = c(21,24,29,31,26)),
  data.frame(id = c("KA","LU","MI","NO","PA"), sex = c(2,2,1,1,2), age_years = c(35,28,22,40,33)),
  data.frame(id = c("Q","R","T","U","V"), sex = c(1,1,2,2,1), age_years = c(27,26,34,29,31))
)

dyads2 <- list(
  list(
    Y  = matrix(c(0,1,0,1,0, 1,0,1,0,0, 0,1,0,1,1, 1,0,1,0,0, 0,0,1,0,0), nrow=5, byrow=TRUE),
    X1 = matrix(c(0,2,0,3,0, 2,0,1,0,0, 0,1,0,4,2, 3,0,4,0,1, 0,0,2,1,0), nrow=5, byrow=TRUE)
  ),
  list(
    Y  = matrix(c(0,1,1,0,0, 1,0,0,1,0, 1,0,0,1,1, 0,1,1,0,0, 0,0,1,0,0), nrow=5, byrow=TRUE),
    X1 = matrix(c(0,3,2,0,0, 3,0,0,1,0, 2,0,0,4,5, 0,1,4,0,0, 0,0,5,0,0), nrow=5, byrow=TRUE)
  ),
  list(
    Y  = matrix(c(0,0,1,1,0, 0,0,1,0,1, 1,1,0,0,1, 1,0,0,0,0, 0,1,1,0,0), nrow=5, byrow=TRUE),
    X1 = matrix(c(0,1,2,3,0, 1,0,4,0,2, 2,4,0,5,6, 3,0,5,0,0, 0,2,6,0,0), nrow=5, byrow=TRUE)
  )
)

# RHS expressions (good)
rhs_cliques <- quote(cliques(k = 2))
rhs_inert1  <- quote(inertia_groups(past_influence = 1))
rhs_combo   <- quote(cliques(k = 2) + inertia_groups(past_influence = 1))

# --------------------------------------------------------------------------------------
# Scenarios (~15) focused on validators
# --------------------------------------------------------------------------------------

# 0) PASS baseline (dataset #1)
.expect_pass(
  "P0) baseline ok (T=3, PLE, nodes/dyads ok, cliques)",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="empile", nodes=nodes1, dyads=dyads1)
)

# 1) FAIL: mode invalid
.expect_fail(
  "F1) mode invalid",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="banana", nodes=NULL, dyads=NULL),
  pattern = "Invalid mode"
)

# 2) FAIL: mode=PLS/sequential blocked
.expect_fail(
  "F2) mode=sequential (PLS) not implemented",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="sequential", nodes=NULL, dyads=NULL),
  pattern = "PLS|sequential"
)

# 3) FAIL: formula not formula
.expect_fail(
  "F3) formula not a formula",
  .call_erpm_long("not a formula", mode="empile"),
  pattern = "formula must be a formula"
)

# 4) FAIL: RHS empty
bad_rhs_empty <- as.formula("partitions1 ~ 1")
.expect_fail(
  "F4) RHS has no terms",
  .call_erpm_long(bad_rhs_empty, mode="empile"),
  pattern = "RHS must contain at least one term"
)

# 5) FAIL: LHS evaluates to non-list
bad_formula_lhs_scalar <- as.formula("c(1,2,3) ~ cliques(k=2)")
.expect_fail(
  "F5) LHS not a list of partitions",
  .call_erpm_long(bad_formula_lhs_scalar, mode="empile"),
  pattern = "LHS must evaluate to a non-empty list"
)

# 6) FAIL: T=1 rejected
partitions_T1 <- list(c(1,1,2,3))
.expect_fail(
  "F6) T=1 rejected (use erpm)",
  .call_erpm_long(partitions_T1 ~ rhs_cliques, mode="empile"),
  pattern = "T=1|Use erpm\\(\\)"
)

# 7) FAIL: partitions contains NA
partitions_NA <- partitions1
partitions_NA[[2]][2] <- NA
.expect_fail(
  "F7) partitions contains NA",
  .call_erpm_long(partitions_NA ~ rhs_cliques, mode="empile"),
  pattern = "contains NA"
)

# 8) FAIL: partitions element is list (not atomic)
partitions_bad_atomic <- partitions1
partitions_bad_atomic[[2]] <- list(1,2,3,4)
.expect_fail(
  "F8) partitions[[t]] not atomic",
  .call_erpm_long(partitions_bad_atomic ~ rhs_cliques, mode="empile"),
  pattern = "must be an atomic vector"
)

# 9) FAIL: verbose not scalar logical
.expect_fail(
  "F9) verbose not scalar logical",
  erpm_long(partitions1 ~ rhs_cliques, mode="empile", eval.call=FALSE, verbose="yes"),
  pattern = "verbose must be TRUE or FALSE"
)

# 10) FAIL: debug not scalar logical
.expect_fail(
  "F10) debug not scalar logical",
  erpm_long(partitions1 ~ rhs_cliques, mode="empile", eval.call=FALSE, debug="yes"),
  pattern = "debug must be TRUE or FALSE"
)

# 11) FAIL: nodes wrong length
nodes_bad_len <- nodes1[1:2]
.expect_fail(
  "F11) nodes length != T",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="empile", nodes=nodes_bad_len),
  pattern = "Invalid nodes|nodes length"
)

# 12) FAIL: nodes missing label
nodes_no_label <- nodes1
nodes_no_label[[1]] <- subset(nodes_no_label[[1]], select = -label)
.expect_fail(
  "F12) nodes[[1]] missing 'label'",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="empile", nodes=nodes_no_label),
  pattern = "must contain a 'label' column"
)

# 13) FAIL: nodes schema differs across time (rename a column)
nodes_schema_diff <- nodes1
names(nodes_schema_diff[[2]])[names(nodes_schema_diff[[2]]) == "gender"] <- "sex"
.expect_fail(
  "F13) nodes schema differs across time",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="empile", nodes=nodes_schema_diff),
  pattern = "column schema differs"
)

# 14) FAIL: nodes row count mismatch
nodes_row_mismatch <- nodes1
nodes_row_mismatch[[3]] <- nodes_row_mismatch[[3]][1:3, , drop = FALSE]
.expect_fail(
  "F14) nodes nrow != length(partition)",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="empile", nodes=nodes_row_mismatch),
  pattern = "has .* rows but partitions"
)

# 15) FAIL: dyads wrong length
dyads_bad_len <- dyads1[1:2]
.expect_fail(
  "F15) dyads length != T",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="empile", dyads=dyads_bad_len),
  pattern = "Invalid dyads|dyads must be a list of length T"
)

# 16) FAIL: dyads[[t]] not named
dyads_no_names <- dyads1
dyads_no_names[[1]] <- unname(dyads_no_names[[1]])
.expect_fail(
  "F16) dyads[[1]] not a NAMED list",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="empile", dyads=dyads_no_names),
  pattern = "must be a NAMED list"
)

# 17) FAIL: dyads matrix not numeric
dyads_non_numeric <- dyads1
dyads_non_numeric[[1]]$Z1 <- matrix(as.character(dyads_non_numeric[[1]]$Z1), nrow=4)
.expect_fail(
  "F17) dyads matrix non-numeric",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="empile", dyads=dyads_non_numeric),
  pattern = "must be numeric"
)

# 18) FAIL: dyads matrix wrong dims
dyads_wrong_dim <- dyads1
dyads_wrong_dim[[2]]$Z1 <- matrix(0, nrow=3, ncol=3)
.expect_fail(
  "F18) dyads matrix wrong dims",
  .call_erpm_long(partitions1 ~ rhs_cliques, mode="empile", dyads=dyads_wrong_dim),
  pattern = "expected 4x4|dim"
)

# 19) FAIL: inertia_groups past_influence too large for T
.expect_fail(
  "F19) inertia_groups past_influence >= T",
  .call_erpm_long(partitions1 ~ inertia_groups(past_influence = 3), mode="empile"),
  pattern = "need at least T|past_influence"
)

# 20) PASS: dataset #2 baseline (different column names; nodes validator expects 'label' so we pass nodes2=NULL)
.expect_pass(
  "P20) dataset #2 baseline ok with nodes=NULL dyads=NULL cliques",
  .call_erpm_long(partitions2 ~ rhs_cliques, mode="empile", nodes=NULL, dyads=NULL)
)

cat("\n=== VALIDATOR SELFTEST DONE ===\n")