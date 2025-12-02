# ======================================================================================
# Fichier : scripts/test/user_tests/tests_single_effect_squared_sizes.R
#
# Objet   : Test single effect "squared_sizes" with different options
#
# ======================================================================================

Sys.setenv(LANG="fr_FR.UTF-8")
try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent=TRUE)
options(encoding="UTF-8")

# ----- packages -----------------------------------------------------------------------
# ----- Dépendances minimales ----------------------------------------------------------
suppressPackageStartupMessages({
  library(devtools)  # load_all(".")
  library(network)   # objet network + attribut 'bipartite'
  library(ergm)      # summary(), ergm(), control.ergm(), etc.
})
devtools::load_all(".")

if (!requireNamespace("Rglpk", quietly = TRUE)) {
  message("Avis: 'Rglpk' non installé. ergm utilisera 'lpSolveAPI'.")
}

options(ergm.loglik.warn_dyads = FALSE)

# ----- Active le patch {ergm} ---------------------------------------------------------
source("scripts/ergm_patch.R")
ergm_patch_enable()

# ----- Partitions de test -------------------------------------------------------------
partition_mix <- c(1, 2, 2, 3, 3, 3)
partition_balanced <- c(1, 1, 2, 2, 3, 3)
partition_full <- c(1, 1, 1, 1, 1, 1)
partition_singleton <- c(1, 2, 3, 4, 5, 6)

# ======================================================================================
# 1) STAT OBSERVÉE SANS FIT — simple cases + options
# ======================================================================================

# baseline test
dry <- erpm(partition_mix ~ squared_sizes, eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 14
dry <- erpm(partition_balanced ~ squared_sizes, eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 12
dry <- erpm(partition_full ~ squared_sizes, eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 36
dry <- erpm(partition_singleton ~ squared_sizes, eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 6

# with options on single size
dry <- erpm(partition_mix ~ squared_sizes(sizes=2), eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 4
dry <- erpm(partition_balanced ~ squared_sizes(sizes=2), eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 12
dry <- erpm(partition_full ~ squared_sizes(sizes=2), eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0
dry <- erpm(partition_singleton ~ squared_sizes(sizes=2), eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0

# with options on several sizes
dry <- erpm(partition_mix ~ squared_sizes(sizes=2:6), eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 13
dry <- erpm(partition_balanced ~ squared_sizes(sizes=2:6), eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 12
dry <- erpm(partition_full ~ squared_sizes(sizes=2:6), eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 36
dry <- erpm(partition_singleton ~ squared_sizes(sizes=2:6), eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0

# ======================================================================================
# 2) FIT MODEL
# ======================================================================================

ctrl_A <- control.ergm(
  CD.maxit        = 0,        # saute la phase contrastive
  MCMLE.maxit     = 10,
  MCMC.burnin     = 5000,
  MCMC.interval   = 1000,
  MCMC.samplesize = 1e4,
  force.main      = TRUE,
  parallel        = 0
)

make_nw_from_partition <- function(part,nodes) {
  built <- build_bipartite_from_inputs(
    partition = part,
    nodes     = nodes
    )
  built$network
}

# baseline case
set.seed(1)
nw <- make_nw_from_partition(partition_mix, NULL)
fit_ergm <- ergm(nw ~ squared_sizes,
                #  estimate="MLE", 
                #  control=ctrl_A,
                 constraints = ~b1part
                 )
print(summary(fit_ergm))

set.seed(1)
fit_erpm <- erpm(partition_mix ~ squared_sizes,
                #  estimate="MLE", 
                #  control=ctrl_A,
                 nodes = NULL
                 ) 
print(summary(fit_erpm))
fit_ergm$coefficients[1] - fit_erpm$coefficients[1]  # should be close to 0

ergm_patch_disable()