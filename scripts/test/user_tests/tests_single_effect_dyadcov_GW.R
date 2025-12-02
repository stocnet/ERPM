# ======================================================================================
# Fichier : scripts/test/user_tests/tests_single_effect_dyadcov_GW.R
#
# Objet   : Test single effect dyadcov_GW avec différentes partitions et matrices
#
# ======================================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE)
options(encoding = "UTF-8")
options(ergm.loglik.warn_dyads = FALSE)

# ----- Packages -----------------------------------------------------------------------
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
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# ======================================================================================
# Partitions de test (n = 5)
# ======================================================================================

partition_mix       <- c(1, 1, 2, 2, 3)  # 2,2,1
partition_balanced  <- c(1, 2, 2, 3, 3)  # 1,2,2
partition_full      <- c(1, 1, 1, 1, 1)  # 5
partition_singleton <- 1:5               # 1,1,1,1,1


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
dry <- erpm(partition_mix ~ dyadcov_GW("block_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0+(2)+(6-8/2)=4
dry <- erpm(partition_balanced ~ dyadcov_GW("block_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be (2)+0+(2)=4
dry <- erpm(partition_full ~ dyadcov_GW("block_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 12-1/2x2x2^3+0...=4
dry <- erpm(partition_singleton ~ dyadcov_GW("block_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0

# baseline test - mixed matrix 
dry <- erpm(partition_mix ~ dyadcov_GW("mix_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0+2+(4-1/2x2)=5 
dry <- erpm(partition_balanced ~ dyadcov_GW("mix_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 2+0+1=3
dry <- erpm(partition_full ~ dyadcov_GW("mix_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 12-1/x8+0...=8
dry <- erpm(partition_singleton ~ dyadcov_GW("mix_att"), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0 

# baseline test - block matrix -> error
dry <- erpm(partition_mix ~ dyadcov_GW("block_att",lambda=3), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0+(2)+(6-8/3)=5.3333
dry <- erpm(partition_balanced ~ dyadcov_GW("block_att",lambda=3), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be (2)+0+(2)=4
dry <- erpm(partition_full ~ dyadcov_GW("block_att",lambda=3), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 12-1/3x2x2^3+0...=6.66667
dry <- erpm(partition_singleton ~ dyadcov_GW("block_att",lambda=3), 
            dyads = nets_df,
            eval_call = FALSE, verbose = TRUE)
print(summary(dry[[2]], constraints = ~ b1part)) # should be 0

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

# baseline case 
set.seed(1)  
nw <- make_nw_from_partition(partition_balanced,nets_df)
fit_ergm <- ergm( nw ~ dyadcov_GW("mix_att"),
                  constraints = ~b1part, 
                  estimate="MLE", 
                  control=ctrl_A)
print(summary(fit_ergm))

set.seed(1)  
fit_erpm <- erpm(partition_balanced ~ dyadcov_GW("mix_att"),
                 dyads = nets_df,
                 estimate="MLE", 
                 control=ctrl_A) 
print(summary(fit_erpm))
print(fit_ergm$coefficients[1] - fit_erpm$coefficients[1])  # should be 0 with the call of the same seed for each case


# option case 
set.seed(1)  
nw <- make_nw_from_partition(partition_balanced,nets_df)
fit_ergm <- ergm( nw ~ dyadcov_GW("mix_att", lambda=3),
                  constraints = ~b1part, 
                  estimate="MLE", 
                  control=ctrl_A)
print(summary(fit_ergm))

set.seed(1)  
fit_erpm <- erpm(partition_balanced ~ dyadcov_GW("mix_att", lambda=3),
                 dyads = nets_df,
                 estimate="MLE", 
                 control=ctrl_A) 
print(summary(fit_erpm))
print(fit_ergm$coefficients[1] - fit_erpm$coefficients[1])  # should be 0 with the call of the same seed for each case



# # ======================================================================================
# # Matrices dyadiques Z1, Z2 (5x5)
# # Z1 : symétrique, diag = 0
# # Z2 : non symétrique, diag = 0
# # ======================================================================================
# 
# Z1 <- matrix(
#   c(0.0, 1.0, 0.5, 0.3, 0.8,
#     1.0, 0.0, 1.2, 0.4, 0.2,
#     0.5, 1.2, 0.0, 0.9, 0.6,
#     0.3, 0.4, 0.9, 0.0, 1.1,
#     0.8, 0.2, 0.6, 1.1, 0.0),
#   nrow = 5, byrow = TRUE
# )
# 
# Z2 <- matrix(
#   c(0.0, 0.3, 0.7, 1.2, 0.5,
#     0.4, 0.0, 0.6, 0.9, 1.5,
#     1.1, 0.2, 0.0, 0.8, 0.3,
#     0.9, 1.4, 0.5, 0.0, 1.0,
#     0.2, 0.6, 1.3, 0.7, 0.0),
#   nrow = 5, byrow = TRUE
# )
# 
# nets_df <- list(Z1 = Z1, Z2 = Z2)
# 
# cat("Z1[1:5,1:5] =\n"); print(Z1)
# cat("\nZ2[1:5,1:5] =\n"); print(Z2)

# # ======================================================================================
# # 1) STAT OBSERVÉE SANS FIT — summary() via ERGM vs dry-run ERPM
# # ======================================================================================
# 
# make_nw_from_partition <- function(part, dyads) {
#   built <- build_bipartite_from_inputs(
#     partition = part,
#     dyads     = dyads
#   )
#   if (is.list(built) && !is.null(built$network)) return(built$network)
#   if (inherits(built, "network")) return(built)
#   stop("build_bipartite_from_inputs() n'a pas renvoyé un 'network'.")
# }
# 
# run_one_summary_check <- function(partition, dyads, dyad_name, lambda) {
#   cat("\n--- SUMMARY CHECK ---\n")
#   cat("partition : ", paste(partition, collapse = " "), "\n", sep = "")
#   cat("dyad      : ", dyad_name, ", lambda = ", lambda, "\n", sep = "")
# 
#   # Référence via network explicite
#   nw <- make_nw_from_partition(partition, dyads)
#   form_ref <- as.formula(paste0("nw ~ dyadcov_GW('", dyad_name, "', lambda = ", lambda, ")"))
#   summary_ref <- as.numeric(summary(form_ref, constraints = ~ b1part))
# 
#   # Via ERPM (dry-run)
#   form_erpm <- as.formula(paste0("partition ~ dyadcov_GW('", dyad_name, "', lambda = ", lambda, ")"))
#   dry <- erpm(
#     form_erpm,
#     eval_call = FALSE,
#     verbose   = TRUE,
#     dyads     = dyads
#   )
#   summary_erpm <- as.numeric(summary(dry[[2]], constraints = ~ b1part))
# 
#   cat(sprintf("[summary] dyadcov_GW('%s', lambda=%g) : ERPM-dry=%g | ERGM-ref=%g\n",
#               dyad_name, lambda, summary_erpm, summary_ref))
#   stopifnot(all.equal(summary_erpm, summary_ref, tol = 0))
# }
# 
# cat("\n=== PHASE 1 : Summary(nw) vs Summary(erpm-dry-run) [dyadcov_GW] ===\n")
# 
# # Partition mix
# run_one_summary_check(partition_mix,       nets_df, "Z1", 2)
# run_one_summary_check(partition_mix,       nets_df, "Z1", 3)
# run_one_summary_check(partition_mix,       nets_df, "Z2", 2)
# 
# # Partition décalée
# run_one_summary_check(partition_shifted,   nets_df, "Z1", 2)
# run_one_summary_check(partition_shifted,   nets_df, "Z1", 3)
# run_one_summary_check(partition_shifted,   nets_df, "Z2", 2)
# 
# # Partition full
# run_one_summary_check(partition_full,      nets_df, "Z1", 2)
# run_one_summary_check(partition_full,      nets_df, "Z1", 3)
# run_one_summary_check(partition_full,      nets_df, "Z2", 2)
# 
# # Partition singleton
# run_one_summary_check(partition_singleton, nets_df, "Z1", 2)
# run_one_summary_check(partition_singleton, nets_df, "Z1", 3)
# run_one_summary_check(partition_singleton, nets_df, "Z2", 2)
# 
# cat("\nFin PHASE 1 : toutes les égalités summary ERGM vs ERPM-dry-run ont été vérifiées.\n")
# 
# # ======================================================================================
# # 2) FIT MODEL — ergm() direct vs erpm()
# # ======================================================================================
# 
# cat("\n=== PHASE 2 : Fits courts (ergm vs erpm) [dyadcov_GW] ===\n")
# 
# # Baseline : partition_mix + Z1, lambda = 2
# set.seed(1)
# nw_ref_Z1 <- make_nw_from_partition(partition_mix, nets_df)
# fit_ergm_Z1_l2 <- ergm(
#   nw_ref_Z1 ~ dyadcov_GW("Z1", lambda = 2),
#   constraints = ~ b1part,
#   eval.loglik = TRUE
# )
# 
# set.seed(1)
# fit_erpm_Z1_l2 <- erpm(
#   partition_mix ~ dyadcov_GW("Z1", lambda = 2),
#   dyads       = nets_df,
#   eval.loglik = TRUE
# )
# 
# cat("\n--- summary(fit_ergm_Z1_l2) ---\n")
# print(summary(fit_ergm_Z1_l2))
# 
# cat("\n--- summary(fit_erpm_Z1_l2) ---\n")
# print(summary(fit_erpm_Z1_l2))
# 
# diff_Z1 <- fit_ergm_Z1_l2$coefficients[1] - fit_erpm_Z1_l2$coefficients[1]
# cat("\nDifférence coef (ergm - erpm) pour dyadcov_GW('Z1', lambda=2) : ", diff_Z1, "\n", sep = "")
# 
# # Baseline : partition_mix + Z2, lambda = 2
# set.seed(1)
# nw_ref_Z2 <- make_nw_from_partition(partition_mix, nets_df)
# fit_ergm_Z2_l2 <- ergm(
#   nw_ref_Z2 ~ dyadcov_GW("Z2", lambda = 2),
#   constraints = ~ b1part,
#   eval.loglik = TRUE
# )
# 
# set.seed(1)
# fit_erpm_Z2_l2 <- erpm(
#   partition_mix ~ dyadcov_GW("Z2", lambda = 2),
#   dyads       = nets_df,
#   eval.loglik = TRUE
# )
# 
# cat("\n--- summary(fit_ergm_Z2_l2) ---\n")
# print(summary(fit_ergm_Z2_l2))
# 
# cat("\n--- summary(fit_erpm_Z2_l2) ---\n")
# print(summary(fit_erpm_Z2_l2))
# 
# diff_Z2 <- fit_ergm_Z2_l2$coefficients[1] - fit_erpm_Z2_l2$coefficients[1]
# cat("\nDifférence coef (ergm - erpm) pour dyadcov_GW('Z2', lambda=2) : ", diff_Z2, "\n", sep = "")
# 
# cat("\nTest single effect dyadcov_GW terminé.\n")
# 
# # Désactive le patch si disponible
# if (exists("ergm_patch_disable", mode = "function")) {
#   ergm_patch_disable()
# }