# ==============================================================================
# MWE minimal : cliques_GW(lambda) — summary + fit
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  library(devtools)  # charge le package local qui expose erpm()
  library(ergm)
})

# Charge le package local (expose erpm(), cliques_GW, build_bipartite_from_inputs, etc.)
devtools::load_all(".")

source("scripts/ergm_patch.R")
ergm_patch_enable()

# Partition de test (tailles = 1,4,2,3)
partition <- c(1,2,2,2,2,3,3,4,4,4)

# Fermé: S(n,λ)=λ*[1-((λ-1)/λ)^n]
S_gw <- function(n, lambda) lambda * (1 - ((lambda - 1) / lambda)^n)
expected_stat <- function(part, lambda) {
  sz <- as.integer(table(part))
  if (length(lambda) == 1L) sum(S_gw(sz, lambda)) else vapply(lambda, \(l) sum(S_gw(sz,l)), 0.0)
}

# ------------------------------------------------------------------------------
# 1) SUMMARY via dry-run erpm -> formule -> summary(formule, ~b1part)
# ------------------------------------------------------------------------------
lambda <- c(2, 4)

dry <- erpm(partition ~ cliques_GW(lambda = lambda), eval_call = FALSE, verbose = FALSE)
obs  <- as.numeric(summary(dry[[2]], constraints = ~ b1part))
exp  <- expected_stat(partition, lambda)

cat("[summary] observé :", paste(obs, collapse=", "),
    " | attendu :", paste(exp, collapse=", "), "\n")
stopifnot(isTRUE(all.equal(obs, exp, tol = 1e-10)))

# ------------------------------------------------------------------------------
# 2) FIT ERPM simple (MLE + logLik) pour λ = 2
# ------------------------------------------------------------------------------
fit <- erpm(partition ~ cliques_GW(lambda = 2),
            estimate    = "MLE",
            eval.loglik = TRUE,
            verbose     = TRUE)

print(summary(fit))