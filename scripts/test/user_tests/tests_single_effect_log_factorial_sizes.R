# ======================================================================================
# Fichier : scripts/test/user_tests/tests_single_effect_log_factorial_sizes.R
#
# Objet   : Test single effect "log_factorial_sizes" with different options
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

options(ergm.loglik.warn_dyads=FALSE)

# ----- Active le patch {ergm} ---------------------------------------------------------
source("scripts/ergm_patch.R")
ergm_patch_enable()

# ----- Partitions de test -------------------------------------------------------------
partition_mix <- c(1, 2, 2, 3, 3, 3)
partition_balanced <- c(1, 1, 2, 2, 2, 3, 4, 4, 5) 
partition_full <- c(1, 1, 1, 1, 1, 1)
partition_singleton <- c(1, 2, 3, 4, 5, 6)

# ======================================================================================
# 1) STAT OBSERVÉE SANS FIT — simple cases + options
# ======================================================================================

# baseline test
dry <- erpm(partition_mix ~ log_factorial_sizes, eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be log(2)=0.6931472
dry <- erpm(partition_balanced ~ log_factorial_sizes, eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0
dry <- erpm(partition_full ~ log_factorial_sizes, eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be log(2*3*4*5) = 4.787492
dry <- erpm(partition_singleton ~ log_factorial_sizes, eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0



# ======================================================================================
# 2) FIT MODEL
# ======================================================================================

set.seed(1)  # stabilise l’estimation si on utilise une estimation CD dans ergm

make_nw_from_partition <- function(part) {
  n <- length(part)
  nw <- network.initialize(n * 2, dir = FALSE, bip = n)
  nw[cbind(seq_along(part), part + n)] <- 1
  nw
}


# create nw for ergm
# n <- 6
# nw <- network.initialize(n * 2, dir = FALSE, bip = n)

# baseline case #1
# ctrl_A <- control.ergm(
#   CD.maxit        = 0,        # saute la phase contrastive
#   MCMLE.maxit     = 10,
#   MCMC.burnin     = 5000,
#   MCMC.interval   = 1000,
#   MCMC.samplesize = 1e4,
#   force.main      = TRUE,
#   parallel        = 0
# )
# nw[cbind(seq_along(partition_mix), partition_mix + n)] <- 1
set.seed(1) 
nw <- make_nw_from_partition(partition_mix)
fit_ergm <- ergm( nw ~ log_factorial_sizes,
                  # estimate="MLE", 
                  # control=ctrl_A,
                  constraints = ~b1part
                  )
print(summary(fit_ergm))

set.seed(1) 
fit_erpm <- erpm(partition_mix ~ log_factorial_sizes) # should be around -0.1
print(summary(fit_erpm))
cat("[ERPM vs ERGM | 1]\n\t", sprintf("fit_ergm - fit_erpm = %f", fit_ergm$coefficients[1] - fit_erpm$coefficients[1]), "\n")  # should be close to 0

# baseline case #2
set.seed(1) 
nw <- make_nw_from_partition(partition_balanced)
fit_ergm <- ergm( nw ~ log_factorial_sizes,
                  # estimate="MLE",
                  # control=ctrl_A,
                  constraints = ~b1part
                  )
print(summary(fit_ergm))

set.seed(1) 
fit_erpm <- erpm( partition_balanced ~ log_factorial_sizes
                  # control=ctrl_A
                  ) # should be around 0.4
print(summary(fit_erpm))
# fit_ergm$coefficients[1] - fit_erpm$coefficients[1]  # should be close to 0
cat("[ERPM vs ERGM | 2]\n\t", sprintf("fit_ergm - fit_erpm = %f", fit_ergm$coefficients[1] - fit_erpm$coefficients[1]), "\n")  # should be close to 0
