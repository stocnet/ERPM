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
# source("scripts/ergm_patch.R")
# ergm_patch_enable()

# ----- Partitions de test -------------------------------------------------------------
partition_mix <- c(1, 2, 2, 3, 3, 3)
partition_balanced <- c(1, 1, 2, 2, 3, 3)
partition_full <- c(1, 1, 1, 1, 1, 1)
partition_singleton <- c(1, 2, 3, 4, 5, 6)

# ----- Binary attributes --------------------------------------------------------------
bin_att <- c(1, 1, 1, 0, 0, 0)
cat_att <- c("A", "A", "B", "B", "C", "C")
nodes_df <- data.frame(label = 1:6, bin_att = bin_att, bin_cat = cat_att)

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
# 2) FIT MODELS
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

# groups+ squared_sizes
fit_erpm <- erpm(partition_mix ~ groups+ squared_sizes,
                 nodes = nodes_df,
                 dyads = nets_df,
                 estimate="MLE") 
print(summary(fit_erpm))

# groups+ cov_match + dyadcov
fit_erpm <- erpm(partition_mix ~ groups+ cov_match(clique_size=3, "bin_att") + dyadcov(clique_size=2, "block_att"),
                 nodes = nodes_df,
                 dyads = nets_df,
                 estimate="MLE", 
                 control=ctrl_A) 
print(summary(fit_erpm))

# groups+ cliques_GW + cov_diff
fit_erpm <- erpm(partition_balanced ~ groups+ cliques_GW + cov_diff("bin_att"),
                 nodes = nodes_df,
                 dyads = nets_df,
                 estimate="MLE", 
                 control=ctrl_A) 
print(summary(fit_erpm))

