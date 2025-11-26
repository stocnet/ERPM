# ======================================================================================
# Fichier : scripts/test/user_tests/tests_single_effect_dyadcov_full.R
#
# Objet   : Test single effect dyadcov_full with different options
#
# ======================================================================================
Sys.setenv(LANG="fr_FR.UTF-8")
try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent=TRUE)
options(encoding="UTF-8")
options(ergm.loglik.warn_dyads=FALSE)
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

# ----- Active le patch {ergm} ---------------------------------------------------------
source("scripts/ergm_patch.R")
ergm_patch_enable()

# ----- Partitions de test -------------------------------------------------------------
partition_mix <- c(1, 2, 2, 3, 3, 3)
partition_balanced <- c(1, 1, 2, 2, 3, 3)
partition_full <- c(1, 1, 1, 1, 1, 1)
partition_singleton <- c(1, 2, 3, 4, 5, 6)

# ----- Binary attributes --------------------------------------------------------------
block_att <- matrix(c(1, 1, 1, 0, 0, 0,
                      1, 1, 1, 0, 0, 0,
                      1, 1, 1, 0, 0, 0,
                      0, 0, 0, 1, 1, 1,
                      0, 0, 0, 1, 1, 1,
                      0, 0, 0, 1, 1, 1), nrow=6, byrow=TRUE)
mix_att <- matrix(c(0, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0,
                    0, 1, 0, 0, 1, 0,
                    0, 0, 0, 0, 1, 1,
                    0, 0, 0, 1, 0, 1,
                    0, 1, 0, 0, 0, 0), nrow=6, byrow=TRUE)
nets_df <- list(block_att = block_att, mix_att = mix_att)

# ======================================================================================
# 1) STAT OBSERVÉE SANS FIT — simple cases + options
# ======================================================================================

# baseline test - block matrix -> error
dry <- erpm(partition_mix ~ dyadcov_full("block_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0+2+6=8 -> error
dry <- erpm(partition_balanced ~ dyadcov_full("block_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 2+0+2=4 -> error
dry <- erpm(partition_full ~ dyadcov_full("block_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 12 -> error
dry <- erpm(partition_singleton ~ dyadcov_full("block_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0 -> error

# baseline test - mixed matrix -> error (the matrix does not have to be symmetric)
dry <- erpm(partition_mix ~ dyadcov_full("mix_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0+2+4=6 -> error
dry <- erpm(partition_balanced ~ dyadcov_full("mix_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 2+0+1=3 -> error
dry <- erpm(partition_full ~ dyadcov_full("mix_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 12 -> error
dry <- erpm(partition_singleton ~ dyadcov_full("mix_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0 -> error

# option test - size = 2 TODO

# option test - size = 2:6 TODO

# ======================================================================================
# 2) FIT MODEL
# ======================================================================================

make_nw_from_partition <- function(part,dyads) {
  built <- build_bipartite_from_inputs(
    partition = part,
    dyads     = dyads)
  built$network
}

ctrl_A <- control.ergm(
  CD.maxit        = 0,        # saute la phase contrastive
  MCMLE.maxit     = 10,
  MCMC.burnin     = 5000,
  MCMC.interval   = 1000,
  MCMC.samplesize = 1e4,
  force.main      = TRUE,
  parallel        = 0
)

# baseline case TODO
set.seed(1)  
nw <- make_nw_from_partition(partition_balanced,nets_df)
fit_ergm <- ergm( nw ~ dyadcov_full("block_att"),
                  constraints = ~b1part, 
                  estimate="MLE", 
                  control=ctrl_A)
print(summary(fit_ergm))

set.seed(1)  
fit_erpm <- erpm(partition_balanced ~ dyadcov_full("block_att"),
                 dyads = nets_df,
                 estimate="MLE", 
                 control=ctrl_A) 
print(summary(fit_erpm))
print(fit_ergm$coefficients[1] - fit_erpm$coefficients[1])  # should be 0 with the call of the same seed for each case

