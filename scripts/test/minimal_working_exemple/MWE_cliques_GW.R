# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cliques_GW.R
#
# Objet   : MWE pour l’effet ERPM `cliques_GW(lambda)` — calcule la stat observée
#            via dry-run `erpm()` + réalise un fit ERGM via `erpm()`.
#
# Chaîne : partition -> biparti -> traduction -> (dry) formule -> summary(formule)
#          partition -> biparti -> traduction -> formule -> ergm()
#
# Référence fermée : S(n,λ) = λ * [ 1 - ((λ-1)/λ)^n ];  stat = Σ_g S(n_g,λ)
# Domaine pour fit : λ > 1 (évite statistique quasi-constante)
# ==============================================================================

# ----- Préambule locale/UTF-8 -------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

# ----- Dépendances minimales ---------------------------------------------------
suppressPackageStartupMessages({
  library(devtools)  # load_all(".")
  library(network)
  library(ergm)
})

# ----- Charge le package local ERPM -------------------------------------------
devtools::load_all(".")

# ----- Active le patch {ergm} si présent --------------------------------------
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# ----- Partition de test -------------------------------------------------------
#   Groupes effectifs = {1,2,3,4} ; tailles = (1,4,2,3) ; N1 = 10 acteurs
partition <- c(1, 2, 2, 2, 2, 3, 3, 4, 4, 4)

cat("\nPartition :", paste(partition, collapse = ", "), "\n")
sizes <- as.integer(table(partition))
cat("Tailles des groupes (par ordre de label):", paste(sizes, collapse = ", "), "\n")
cat("Nombre d'acteurs (N1):", length(partition), " | Nombre de groupes:", length(sizes), "\n\n")

# ----- Helper forme fermée ----------------------------------------------------
S_gw <- function(n, lambda) lambda * (1 - ((lambda - 1) / lambda)^n)
stat_expected <- function(part, lambda) {
  sz <- as.integer(table(part))
  if (length(lambda) == 1L) {
    sum(S_gw(sz, lambda))
  } else {
    vapply(lambda, function(l) sum(S_gw(sz, l)), numeric(1))
  }
}

# ==============================================================================
# 1) STAT OBSERVÉE — dry-run erpm -> formule -> summary(formule, ~b1part)
# ==============================================================================

# -- Cas scalaire : lambda = 2
dry_l2 <- erpm(partition ~ cliques_GW(lambda = 2), eval_call = FALSE, verbose = TRUE)
fml_l2 <- dry_l2[[2]]
obs_l2 <- summary(fml_l2, constraints = ~ b1part)
exp_l2 <- stat_expected(partition, 2)

cat(sprintf("[MWE cliques_GW(λ=2)] summary observé = %s | attendu = %s\n",
            paste0(obs_l2, collapse=","), paste0(exp_l2, collapse=",")))
stopifnot(is.numeric(obs_l2), length(obs_l2) == 1L,
          isTRUE(all.equal(as.numeric(obs_l2), as.numeric(exp_l2), tol = 1e-10)))

# -- Cas vectorisé : lambda = c(1.25, 4)
dry_vec <- erpm(partition ~ cliques_GW(lambda = c(1.25, 4)), eval_call = FALSE, verbose = TRUE)
fml_vec <- dry_vec[[2]]
obs_vec <- summary(fml_vec, constraints = ~ b1part)
exp_vec <- stat_expected(partition, c(1.25, 4))

cat(sprintf("[MWE cliques_GW(λ=c(1.25,4))] summary observé = %s | attendu = %s\n",
            paste0(obs_vec, collapse=", "), paste0(exp_vec, collapse=", ")))
stopifnot(is.numeric(obs_vec), length(obs_vec) == 2L,
          isTRUE(all.equal(as.numeric(obs_vec), as.numeric(exp_vec), tol = 1e-10)))

# ==============================================================================
# 2) FIT ERGM COMPLET — lambda = 2 (λ>1 recommandé pour l'estimation)
# ==============================================================================

set.seed(1)
fit <- erpm(partition ~ cliques_GW(lambda = 2),
            eval_call   = TRUE,
            verbose     = TRUE,
            estimate    = "MLE",      # ou "CD" si MLE indispo/long
            eval.loglik = TRUE)

cat("\n--- summary(fit) (style ergm) ---\n")
print(summary(fit))

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)