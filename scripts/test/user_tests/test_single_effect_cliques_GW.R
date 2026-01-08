# ======================================================================================
# Fichier : scripts/test/user_tests/tests_single_effect_cliques_GW.R
#
# Objet   : Test single effect "cliques_GW" with different options
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

# ----- Active le patch \pkg{ergm} ---------------------------------------------------------
source("scripts/ergm_patch.R")
ergm_patch_enable()

# ----- Partitions de test -------------------------------------------------------------
partition_mix <- c(1, 2, 2, 3, 3, 3)
nodes_df <- data.frame(
  label  = sprintf("A%d", seq_along(partition_mix)),
  gender = sample(1:2, length(partition_mix), replace = TRUE),
  age    = sample(20:40, length(partition_mix), replace = TRUE),
  stringsAsFactors = FALSE
)

partition_balanced <- c(1, 1, 2, 2, 3, 3)
partition_full <- c(1, 1, 1, 1, 1, 1)
partition_singleton <- c(1, 2, 3, 4, 5, 6)


# ======================================================================================
# 1) STAT OBSERVÉE SANS FIT — simple cases + options
# ======================================================================================

# baseline test
dry <- erpm(partition_mix ~ cliques_GW, eval.call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 4.25
dry <- erpm(partition_balanced ~ cliques_GW, eval.call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 4.5
dry <- erpm(partition_full ~ cliques_GW, eval.call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 1.96875
dry <- erpm(partition_singleton ~ cliques_GW, eval.call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 6

# with options on lambda = 3
dry <- erpm(partition_mix ~ cliques_GW(lambda=3), eval.call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 4.778
dry <- erpm(partition_balanced ~ cliques_GW(lambda=3), eval.call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 5
dry <- erpm(partition_full ~ cliques_GW(lambda=3), eval.call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 2.736626
dry <- erpm(partition_singleton ~ cliques_GW(lambda=3), eval.call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 6

# with options on lambda < 1
#dry <- erpm(partition_mix ~ cliques_GW(lambda=0), eval.call = FALSE, verbose = TRUE)
# should be an error

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
    nodes     = nodes)
  built$network
}

# baseline case
nw <- make_nw_from_partition(partition_mix,nodes_df)
set.seed(1)  
fit_ergm <- ergm(nw ~ cliques_GW,
                 constraints = ~b1part,
                 estimate="MLE", 
                 control=ctrl_A)
print(summary(fit_ergm))

set.seed(1)  
fit_erpm <- erpm(partition_mix ~ cliques_GW,
                 nodes = nodes_df,
                 estimate="MLE", 
                 control=ctrl_A) 
print(summary(fit_erpm))
fit_ergm$coefficients[1] - fit_erpm$coefficients[1]  # should be close to 0

# option with lambda=3
nw <- make_nw_from_partition(partition_mix,nodes_df)
set.seed(1)  
fit_ergm <- ergm(nw ~ cliques_GW(lambda=3),
                 constraints = ~b1part,
                 estimate="MLE", 
                 control=ctrl_A)
print(summary(fit_ergm))

set.seed(1)  
fit_erpm <- erpm(partition_mix ~ cliques_GW(lambda=3),
                 nodes = nodes_df,
                 estimate="MLE", 
                 control=ctrl_A) 
print(summary(fit_erpm))
cat("[ERPM vs ERGM]\n\t", sprintf("fit_ergm - fit_erpm = %f", fit_ergm$coefficients[1] - fit_erpm$coefficients[1]), "\n")  # should be 0

