# ==============================================================================
# MWE: simple usage of erpm(...) with mh_moves / mh_weights
# Purpose:
#   Show how to run ERPM fits:
#     1) without explicit moves
#     2) with toggle/swap
#     3) with toggle/swap/merge/split
#     4) with another combination
#
# Data:
#   - hard-coded partition
#   - hard-coded nodal attribute
#   - hard-coded dyadic attribute
#
# Tested effects:
#   - log_factorial_sizes()
#   - cov_match("bin_att")
#   - dyadcov("Z1")
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

options(ergm.loglik.warn_dyads = FALSE)
options(ERPM.zzz.verbose = TRUE)
options(Proposal.ErpmMix.debug = TRUE)

suppressPackageStartupMessages({
  library(network)
  library(ergm)
})

if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
}

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Run this script from the package root with devtools available.")
}

cat("=== ERPM MWE: mh_moves / mh_weights ===\n")
cat("R:", paste(R.version$major, R.version$minor, sep = "."), "\n")
cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")
cat("ERPM loaded.\n")

# ------------------------------------------------------------------------------
# Hard-coded data
# ------------------------------------------------------------------------------
partition <- c(1,2,2,3,3,3,4,5,5,5,6,6,6,6,6)

nodes <- data.frame(
  label   = 1:15,
  bin_att = c("1","1","1","0","0","0","1","0","1","0","1","0","1","0","1"),
  stringsAsFactors = FALSE
)

Z1 <- matrix(
  c(
    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
    1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,
    2,1,0,1,2,3,4,5,6,7,8,9,10,11,12,
    3,2,1,0,1,2,3,4,5,6,7,8,9,10,11,
    4,3,2,1,0,1,2,3,4,5,6,7,8,9,10,
    5,4,3,2,1,0,1,2,3,4,5,6,7,8,9,
    6,5,4,3,2,1,0,1,2,3,4,5,6,7,8,
    7,6,5,4,3,2,1,0,1,2,3,4,5,6,7,
    8,7,6,5,4,3,2,1,0,1,2,3,4,5,6,
    9,8,7,6,5,4,3,2,1,0,1,2,3,4,5,
    10,9,8,7,6,5,4,3,2,1,0,1,2,3,4,
    11,10,9,8,7,6,5,4,3,2,1,0,1,2,3,
    12,11,10,9,8,7,6,5,4,3,2,1,0,1,2,
    13,12,11,10,9,8,7,6,5,4,3,2,1,0,1,
    14,13,12,11,10,9,8,7,6,5,4,3,2,1,0
  ),
  nrow = 15, ncol = 15, byrow = TRUE
)
storage.mode(Z1) <- "double"

dyads <- list(Z1 = Z1)

# ------------------------------------------------------------------------------
# Baseline MCMC settings for fits
# ------------------------------------------------------------------------------
ctrl_base <- ergm::control.ergm(
  MCMC.burnin     = 20000,
  MCMC.interval   = 1,
  MCMC.samplesize = 40000,
  seed            = 123
)

# ==============================================================================
# 1) No explicit move: let ergm handle transitions
# ==============================================================================
cat("\n============================================================\n")
cat("CASE 1) No mh_moves / mh_weights\n")
cat("============================================================\n")

fit1_a <- erpm(
  partition ~ log_factorial_sizes(),
  nodes    = nodes,
  dyads    = dyads,
  control  = ctrl_base,
  verbose  = TRUE
)
print(summary(fit1_a))

fit1_b <- erpm(
  partition ~ cov_match("bin_att"),
  nodes    = nodes,
  dyads    = dyads,
  control  = ctrl_base,
  verbose  = TRUE
)
print(summary(fit1_b))

fit1_c <- erpm(
  partition ~ dyadcov("Z1"),
  nodes    = nodes,
  dyads    = dyads,
  control  = ctrl_base,
  verbose  = TRUE
)
print(summary(fit1_c))

# ==============================================================================
# 2) Toggle / swap moves
# ==============================================================================
cat("\n============================================================\n")
cat("CASE 2) mh_moves = c('toggle','swap')\n")
cat("============================================================\n")

fit2_a <- erpm(
  partition ~ log_factorial_sizes(),
  nodes      = nodes,
  dyads      = dyads,
  mh_moves   = c("toggle", "swap"),
  mh_weights = c(2, 1),
  control    = ctrl_base,
  verbose    = TRUE
)
print(summary(fit2_a))

fit2_b <- erpm(
  partition ~ cov_match("bin_att"),
  nodes      = nodes,
  dyads      = dyads,
  mh_moves   = c("toggle", "swap"),
  mh_weights = c(2, 1),
  control    = ctrl_base,
  verbose    = TRUE
)
print(summary(fit2_b))

fit2_c <- erpm(
  partition ~ dyadcov("Z1"),
  nodes      = nodes,
  dyads      = dyads,
  mh_moves   = c("toggle", "swap"),
  mh_weights = c(2, 1),
  control    = ctrl_base,
  verbose    = TRUE
)
print(summary(fit2_c))

# ==============================================================================
# 3) Toggle / swap / merge / split moves
# ==============================================================================
cat("\n============================================================\n")
cat("CASE 3) mh_moves = c('toggle','swap','merge','split')\n")
cat("============================================================\n")

fit3_a <- erpm(
  partition ~ log_factorial_sizes(),
  nodes      = nodes,
  dyads      = dyads,
  mh_moves   = c("toggle", "swap", "merge", "split"),
  mh_weights = c(6, 2, 1, 1),
  control    = ctrl_base,
  verbose    = TRUE
)
print(summary(fit3_a))

fit3_b <- erpm(
  partition ~ cov_match("bin_att"),
  nodes      = nodes,
  dyads      = dyads,
  mh_moves   = c("toggle", "swap", "merge", "split"),
  mh_weights = c(6, 2, 1, 1),
  control    = ctrl_base,
  verbose    = TRUE
)
print(summary(fit3_b))

fit3_c <- erpm(
  partition ~ dyadcov("Z1"),
  nodes      = nodes,
  dyads      = dyads,
  mh_moves   = c("toggle", "swap", "merge", "split"),
  mh_weights = c(6, 2, 1, 1),
  control    = ctrl_base,
  verbose    = TRUE
)
print(summary(fit3_c))

# ==============================================================================
# 4) Another combination
#    Example: favor merge / split
# ==============================================================================
cat("\n============================================================\n")
cat("CASE 4) another move combination\n")
cat("============================================================\n")

fit4_a <- erpm(
  partition ~ log_factorial_sizes(),
  nodes      = nodes,
  dyads      = dyads,
  mh_moves   = c("toggle", "swap", "merge", "split"),
  mh_weights = c(2, 1, 4, 4),
  control    = ctrl_base,
  verbose    = TRUE
)
print(summary(fit4_a))

fit4_b <- erpm(
  partition ~ cov_match("bin_att"),
  nodes      = nodes,
  dyads      = dyads,
  mh_moves   = c("toggle", "swap", "merge", "split"),
  mh_weights = c(2, 1, 4, 4),
  control    = ctrl_base,
  verbose    = TRUE
)
print(summary(fit4_b))

fit4_c <- erpm(
  partition ~ dyadcov("Z1"),
  nodes      = nodes,
  dyads      = dyads,
  mh_moves   = c("toggle", "swap", "merge", "split"),
  mh_weights = c(2, 1, 4, 4),
  control    = ctrl_base,
  verbose    = TRUE
)
print(summary(fit4_c))

cat("\nMWE finished.\n")

if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()