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
source("scripts/ergm_patch.R")
ergm_patch_enable()
options(ergm.loglik.warn_dyads=FALSE)

# ----- Partitions de test -------------------------------------------------------------
partition_mix       <- c(1, 2, 2, 3, 3, 3)  # tailles : 1, 2, 3
partition_balanced  <- c(1, 1, 2, 2, 3, 3)  # tailles : 2, 2, 2
partition_full      <- c(1, 1, 1, 1, 1, 1)  # taille : 6
partition_singleton <- c(1, 2, 3, 4, 5, 6)  # 6 singletons

# ----- Nodes pour les fits ------------------------------------------------------------
nodes_df <- data.frame(
  label = paste0("A", 1:6),
  stringsAsFactors = FALSE
)

# ======================================================================================
# 1) STAT OBSERVÉE SANS FIT — simple cases + options
# ======================================================================================
# Rappel : 
#   - pour k >= 2, cliques(k) = somme_g C(n_g, k)
#   - pour k = 1, cliques(1) = nombre de groupes de taille exactement 1
#   - normalized = TRUE : division par C(N1, k), où N1 = nb d'acteurs

# baseline test (k = 2)
dry <- erpm(partition_mix ~ cliques, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 4  (0 + 1 + 3)
dry <- erpm(partition_balanced ~ cliques, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 3  (1 + 1 + 1)
dry <- erpm(partition_full ~ cliques, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 15 (C(6,2))
dry <- erpm(partition_singleton ~ cliques, eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 0  (6 groupes de taille 1)

# with options on clique size = 3
dry <- erpm(partition_mix ~ cliques(clique_size = 3), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 1  (0 + 0 + 1)
dry <- erpm(partition_balanced ~ cliques(clique_size = 3), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 0
dry <- erpm(partition_full ~ cliques(clique_size = 3), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 20 (C(6,3))
dry <- erpm(partition_singleton ~ cliques(clique_size = 3), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 0

# with options on other cliques sizes
dry <- erpm(partition_singleton ~ cliques(clique_size = 1), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 6  (6 groupes de taille 1)
dry <- erpm(partition_mix ~ cliques(clique_size = 1), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 1  (un seul groupe de taille 1)
dry <- erpm(partition_mix ~ cliques(clique_size = 7), eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 0  (aucun groupe de taille >= 7)

# with options normalized = TRUE
# N1 = 6 pour partition_mix
dry <- erpm(partition_mix ~ cliques(clique_size = 2, normalized = TRUE),
            eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 4 / C(6,2) = 4/15 ≈ 0.2667

dry <- erpm(partition_mix ~ cliques(clique_size = 3, normalized = TRUE),
            eval_call = FALSE, verbose = TRUE)
summary(dry[[2]], constraints = ~ b1part) # attendu : 1 / C(6,3) = 1/20 = 0.05

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

make_nw_from_partition <- function(part, nodes) {
  built <- build_bipartite_from_inputs(
    partition = part,
    nodes     = nodes)
  built$network
}

# baseline case (k = 2)
set.seed(1)  # stabilise l’estimation si on utilise une estimation CD dans ergm
nw <- make_nw_from_partition(partition_mix, nodes_df)
fit_ergm <- ergm(nw ~ cliques,
                 constraints = ~ b1part,
                 estimate    = "MLE",
                 control     = ctrl_A)
print(summary(fit_ergm))

set.seed(1)
fit_erpm <- erpm(partition_mix ~ cliques,
                 nodes    = nodes_df,
                 estimate = "MLE",
                 control  = ctrl_A)
print(summary(fit_erpm))
fit_ergm$coefficients[1] - fit_erpm$coefficients[1]  # devrait être ≈ 0

# option with clique_size = 3
set.seed(1)
nw <- make_nw_from_partition(partition_mix, nodes_df)
fit_ergm <- ergm(nw ~ cliques(clique_size = 3),
                 constraints = ~ b1part,
                 estimate    = "MLE",
                 control     = ctrl_A)
print(summary(fit_ergm))

set.seed(1)
fit_erpm <- erpm(partition_mix ~ cliques(clique_size = 3),
                 nodes    = nodes_df,
                 estimate = "MLE",
                 control  = ctrl_A)
print(summary(fit_erpm))
fit_ergm$coefficients[1] - fit_erpm$coefficients[1]  # devrait être ≈ 0
