# ======================================================================================
# Fichier : scripts/test/user_tests/tests_single_effects.R
#
# Objet   : Test single effects one by one with different options
#
# ======================================================================================

# ----- packages -----------------------------------------------------------------------
# ----- Dépendances minimales ----------------------------------------------------------
suppressPackageStartupMessages({
  library(devtools)  # load_all(".")
  library(network)   # objet network + attribut 'bipartite'
  library(ergm)      # summary(), ergm(), control.ergm(), etc.
})
devtools::load_all(".")

# ----- Active le patch {ergm} ---------------------------------------------------------
#source("scripts/ergm_patch.R")
#ergm_patch_enable()

# ----- Partitions de test -------------------------------------------------------------
partition_mix <- c(1, 2, 2, 3, 3, 3)
partition_balanced <- c(1, 1, 2, 2, 3, 3)
partition_full <- c(1, 1, 1, 1, 1, 1)
partition_singleton <- c(1, 2, 3, 4, 5, 6)

# ======================================================================================
# 1) STAT OBSERVÉE SANS FIT — simple cases + options
# ======================================================================================

# baseline test
dry <- erpm(partition_mix ~ log_factorial_sizes, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be log(2)=0.6931472
dry <- erpm(partition_balanced ~ log_factorial_sizes, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 0
dry <- erpm(partition_full ~ log_factorial_sizes, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be log(2*3*4*5) = 4.787492
dry <- erpm(partition_singleton ~ log_factorial_sizes, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 0



# ======================================================================================
# 2) FIT MODEL
# ======================================================================================

set.seed(1)  # stabilise l’estimation si on utilise une estimation CD dans ergm

# create nw for ergm
n <- 6
nw <- network.initialize(n * 2, dir = FALSE, bip = n)

# baseline case #1
nw[cbind(seq_along(partition_mix), partition_mix + n)] <- 1
fit_ergm <- ergm(nw ~ log_factorial_sizes,
                 constraints = ~b1part)
summary(fit_ergm)

fit_erpm <- erpm(partition_mix ~ log_factorial_sizes) # should be around -0.1
summary(fit_erpm) 
fit_ergm$coefficients[1] - fit_erpm$coefficients[1]  # should be close to 0

# baseline case #2
nw[cbind(seq_along(partition_balanced), partition_balanced + n)] <- 1
fit_ergm <- ergm(nw ~ log_factorial_sizes,
                 constraints = ~b1part)
summary(fit_ergm)

fit_erpm <- erpm(partition_balanced ~ log_factorial_sizes) # should be around 0.4
summary(fit_erpm)
fit_ergm$coefficients[1] - fit_erpm$coefficients[1]  # should be close to 0
