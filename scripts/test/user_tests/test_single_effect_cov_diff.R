# ======================================================================================
# Fichier : scripts/test/user_tests/tests_single_effect_cov_diff.R
#
# Objet   : Test single effect cov_diff with different options
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

# ----- Active le patch {ergm} ---------------------------------------------------------
source("scripts/ergm_patch.R")
ergm_patch_enable()

# ----- Partitions de test -------------------------------------------------------------
partition_mix <- c(1, 2, 2, 3, 3, 3)
partition_balanced <- c(1, 1, 2, 2, 3, 3)
partition_full <- c(1, 1, 1, 1, 1, 1)
partition_singleton <- c(1, 2, 3, 4, 5, 6)

# ----- Binary attributes --------------------------------------------------------------
bin_att <- c(1, 1, 1, 0, 0, 0)
cat_att <- c("A", "A", "B", "B", "C", "C")
nodes_df <- data.frame(label = 1:6, bin_att = bin_att, bin_cat = cat_att)

# ======================================================================================
# 1) STAT OBSERVÉE SANS FIT — simple cases + options
# ======================================================================================

# baseline test - binary attribute
dry <- erpm(partition_mix ~ cov_diff("bin_att"), 
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0
dry <- erpm(partition_balanced ~ cov_diff("bin_att"), 
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 1
dry <- erpm(partition_full ~ cov_diff("bin_att"), 
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 9
dry <- erpm(partition_singleton ~ cov_diff("bin_att"),
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0 

# baseline test - category attribute
#dry <- erpm(partition_mix ~ cov_diff("bin_cat"),
#            nodes = nodes_df,
#            eval_call = FALSE, verbose = TRUE) # should be an error


# option test - clique size = 3 -> error
dry <- erpm(partition_mix ~ cov_diff("bin_att", clique_size=3), 
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0
dry <- erpm(partition_balanced ~ cov_diff("bin_att", clique_size=3), 
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0
dry <- erpm(partition_full ~ cov_diff("bin_att", clique_size=3), 
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 18
dry <- erpm(partition_singleton ~ cov_diff("bin_att", clique_size=3),
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0 

# option test - clique size = 3 -> error
dry <- erpm(partition_mix ~ cov_diff("bin_att", normalized = T), 
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0
dry <- erpm(partition_balanced ~ cov_diff("bin_att", normalized = T), 
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0.5
dry <- erpm(partition_full ~ cov_diff("bin_att", normalized = T), 
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 1.5
dry <- erpm(partition_singleton ~ cov_diff("bin_att", normalized = T),
            nodes = nodes_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0 

# ======================================================================================
# 2) FIT MODEL
# ======================================================================================

set.seed(1)  # stabilise l’estimation si on utilise une estimation CD dans ergm

make_nw_from_partition <- function(part,nodes) {
  built <- build_bipartite_from_inputs(
    partition = part,
    nodes     = nodes)
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

# baseline case 
set.seed(1)  # stabilise l’estimation si on utilise une estimation CD dans ergm

nw <- make_nw_from_partition(partition_balanced,nodes_df)
fit_ergm <- ergm( nw ~ cov_diff("bin_att"),
                  constraints = ~b1part, 
                  estimate="MLE", 
                  control=ctrl_A)
print(summary(fit_ergm))

set.seed(1)  # stabilise l’estimation si on utilise une estimation CD dans ergm

fit_erpm <- erpm(partition_balanced ~ cov_diff("bin_att"),
                 nodes = nodes_df,
                 estimate="MLE", 
                 control=ctrl_A) 
print(summary(fit_erpm))
fit_ergm$coefficients[1] - fit_erpm$coefficients[1]  # should be close to 0

# option case TODO
