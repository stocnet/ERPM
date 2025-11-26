# ======================================================================================
# Fichier : scripts/test/user_tests/tests_single_effect_dyadcov.R
#
# Objet   : Test single effect dyadcov with different options
#
# ======================================================================================
Sys.setenv(LANG = "fr_FR.UTF-8")
try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE)
options(encoding = "UTF-8")
options(ergm.loglik.warn_dyads = FALSE)

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
partition_mix       <- c(1, 2, 2, 3, 3, 3)
partition_balanced  <- c(1, 1, 2, 2, 3, 3)
partition_full      <- c(1, 1, 1, 1, 1, 1)
partition_singleton <- c(1, 2, 3, 4, 5, 6)

# ----- Dyadic covariates (6 x 6) ------------------------------------------------------
block_att <- matrix(
  c(0, 1, 1, 0, 0, 0,
    1, 0, 1, 0, 0, 0,
    1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1,
    0, 0, 0, 1, 0, 1,
    0, 0, 0, 1, 1, 0),
  nrow = 6, byrow = TRUE
)

mix_att <- matrix(
  c(0, 1, 1, 0, 0, 1,
    1, 0, 1, 0, 0, 0,
    0, 1, 0, 0, 1, 0,
    0, 0, 0, 0, 1, 1,
    0, 0, 0, 1, 0, 1,
    0, 1, 0, 0, 0, 0),
  nrow = 6, byrow = TRUE
)

nets_df <- list(block_att = block_att, mix_att = mix_att)

# ======================================================================================
# 1) STAT OBSERVÉE SANS FIT — cas simples + options
# ======================================================================================

cat("========== 1) STAT OBSERVÉE SANS FIT — dyadcov (k=2) ==========\n")

# ---- baseline test - block matrix -------------------------------------------
dry <- erpm(
  partition_mix ~ dyadcov("block_att", clique_size = 2, normalized = FALSE),
  dyads     = nets_df,
  eval_call = FALSE,
  verbose   = TRUE
)
print(summary(dry[[2]], constraints = ~ b1part))
# attendu (comme dyadcov_full): 0 + 2 + 6 = 8

dry <- erpm(
  partition_balanced ~ dyadcov("block_att", clique_size = 2, normalized = FALSE),
  dyads     = nets_df,
  eval_call = FALSE,
  verbose   = TRUE
)
print(summary(dry[[2]], constraints = ~ b1part))
# attendu: 2 + 0 + 2 = 4

dry <- erpm(
  partition_full ~ dyadcov("block_att", clique_size = 2, normalized = FALSE),
  dyads     = nets_df,
  eval_call = FALSE,
  verbose   = TRUE
)
print(summary(dry[[2]], constraints = ~ b1part))
# attendu: 12

dry <- erpm(
  partition_singleton ~ dyadcov("block_att", clique_size = 2, normalized = FALSE),
  dyads     = nets_df,
  eval_call = FALSE,
  verbose   = TRUE
)
print(summary(dry[[2]], constraints = ~ b1part))
# attendu: 0

# ---- baseline test - mixed matrix (matrice non symétrique) ------------------
dry <- erpm(
  partition_mix ~ dyadcov("mix_att", clique_size = 2, normalized = FALSE),
  dyads     = nets_df,
  eval_call = FALSE,
  verbose   = TRUE
)
print(summary(dry[[2]], constraints = ~ b1part))
# attendu (comme dyadcov_full): 0 + 2 + 4 = 6

dry <- erpm(
  partition_balanced ~ dyadcov("mix_att", clique_size = 2, normalized = FALSE),
  dyads     = nets_df,
  eval_call = FALSE,
  verbose   = TRUE
)
print(summary(dry[[2]], constraints = ~ b1part))
# attendu: 2 + 0 + 1 = 3

dry <- erpm(
  partition_full ~ dyadcov("mix_att", clique_size = 2, normalized = FALSE),
  dyads     = nets_df,
  eval_call = FALSE,
  verbose   = TRUE
)
print(summary(dry[[2]], constraints = ~ b1part))
# attendu: 12

dry <- erpm(
  partition_singleton ~ dyadcov("mix_att", clique_size = 2, normalized = FALSE),
  dyads     = nets_df,
  eval_call = FALSE,
  verbose   = TRUE
)
print(summary(dry[[2]], constraints = ~ b1part))
# attendu: 0


# ======================================================================================
# 2) FIT MODEL — comparaison ergm vs erpm
# ======================================================================================

cat("========== 2) FIT MODEL — ergm vs erpm (dyadcov, k=2) ==========\n")

make_nw_from_partition <- function(part, dyads) {
  built <- build_bipartite_from_inputs(
    partition = part,
    dyads     = dyads
  )
  built$network
}

ctrl_A <- control.ergm(
  CD.maxit        = 0,      # saute la phase contrastive
  MCMLE.maxit     = 10,
  MCMC.burnin     = 5000,
  MCMC.interval   = 1000,
  MCMC.samplesize = 1e4,
  force.main      = TRUE,
  parallel        = 0
)

# ---- baseline case : partition_balanced + block_att -------------------------
set.seed(1)
nw <- make_nw_from_partition(partition_balanced, nets_df)

fit_ergm <- ergm(
  nw ~ dyadcov("block_att", clique_size = 2, normalized = FALSE),
#   estimate    = "MLE",
#   control     = ctrl_A,
  verbose = TRUE,
  constraints = ~ b1part
)
print(summary(fit_ergm))

set.seed(1)
fit_erpm <- erpm(
  partition_balanced ~ dyadcov("block_att", clique_size = 2, normalized = FALSE),
#   estimate = "MLE",
#   control  = ctrl_A,
  verbose = TRUE,
  dyads    = nets_df
)
print(summary(fit_erpm))

cat("Différence des coefficients (ergm - erpm) :\n")
print(fit_ergm$coefficients[1] - fit_erpm$coefficients[1])
# attendu: ≈ 0 si les chemins d'estimation et la définition de la stat coïncident
