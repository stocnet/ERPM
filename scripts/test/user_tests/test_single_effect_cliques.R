# ======================================================================================
# Fichier : scripts/test/user_tests/tests_single_effect_cliques.R
#
# Objet   : Test single effect "cliques" with different options
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
dry <- erpm(partition_mix ~ cliques, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 4
dry <- erpm(partition_balanced ~ cliques, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 3
dry <- erpm(partition_full ~ cliques, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 15
dry <- erpm(partition_singleton ~ cliques, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 0

# with options on clique size = 3
dry <- erpm(partition_mix ~ cliques(clique_size = 3), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 1
dry <- erpm(partition_balanced ~ cliques(clique_size = 3), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 0
dry <- erpm(partition_full ~ cliques(clique_size = 3), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 20
dry <- erpm(partition_singleton ~ cliques(clique_size = 3), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 0

# with options on other cliques sizes
dry <- erpm(partition_singleton ~ cliques(clique_size = 1), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 6
dry <- erpm(partition_mix ~ cliques(clique_size = 1), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 6 --> error?
dry <- erpm(partition_mix ~ cliques(clique_size = 7), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be an error/warning?

# with options normalized = TRUE
dry <- erpm(partition_mix ~ cliques(clique_size = 2, normalized = T), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 1/2 + 3/3 = 1.5 -> error ?
dry <- erpm(partition_mix ~ cliques(clique_size = 3, normalized = T), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # should be 1/3 -> error ?

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

set.seed(1)  # stabilise l’estimation si on utilise une estimation CD dans ergm

make_nw_from_partition <- function(part,nodes) {
  built <- build_bipartite_from_inputs(
    partition = part,
    nodes     = nodes)
  built$network
}

# baseline case
nw <- make_nw_from_partition(partition_mix,nodes_df)
fit_ergm <- ergm(nw ~ cliques,
                 constraints = ~b1part,
                 estimate="MLE", 
                 control=ctrl_A)
print(summary(fit_ergm))

fit_erpm <- erpm(partition_mix ~ cliques,
                 nodes = nodes_df,
                 estimate="MLE", 
                 control=ctrl_A) 
print(summary(fit_erpm))
fit_ergm$coefficients[1] - fit_erpm$coefficients[1]  # should be close to 0

# option with clique_size=3
nw <- make_nw_from_partition(partition_mix,nodes_df)
fit_ergm <- ergm(nw ~ cliques(clique_size=3),
                 constraints = ~b1part,
                 estimate="MLE", 
                 control=ctrl_A)
print(summary(fit_ergm))

fit_erpm <- erpm(partition_mix ~ cliques(clique_size=3),
                 nodes = nodes_df,
                 estimate="MLE", 
                 control=ctrl_A) 
print(summary(fit_erpm))
fit_ergm$coefficients[1] - fit_erpm$coefficients[1]  # should be close to 0

# option with normalized=T -> TODO