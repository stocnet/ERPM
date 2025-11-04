# ======================================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_squared_sizes.R
#
# Objet   : MWE pour l’effet ERPM `squared_sizes()` — calcule la stat observée (summary)
#            via dry-run `erpm()` et calcule un fit ERGM via `erpm()`.
#
# Contexte: Chaîne ERPM → ERGM :
#            partition -> réseau biparti -> traduction -> (dry) formule -> summary(formule) / ergm()
#            partition -> réseau biparti -> traduction -> formule -> ergm()
#
# Résumé technique :
#   • `squared_sizes()` est implémenté côté ERGM (InitErgmTerm + change-stat C).
#   • La statistique vaut ∑_g |g|^2 sur le mode groupes du biparti.
#   • Pour le summary : `erpm(..., eval_call=FALSE)` renvoie un appel `ergm(formule, ...)` non évalué.
#     La formule capture `nw` dans son environnement ; `summary(formule, ...)` l’utilise.
# ======================================================================================

# ----- Préambule locale/UTF-8 ---------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

# ----- Dépendances minimales ----------------------------------------------------------
suppressPackageStartupMessages({
  library(devtools)  # load_all(".")
  library(network)   # objet network + attribut 'bipartite'
  library(ergm)      # summary(), ergm(), control.ergm(), etc.
})

# ----- Charge le package local ERPM ---------------------------------------------------
devtools::load_all(".")

# ----- Active le patch {ergm} ---------------------------------------------------------
source("scripts/ergm_patch.R")
ergm_patch_enable()

# ----- Partition de test --------------------------------------------------------------
#   Groupes = {1,2,3} ; tailles = (1,2,3) 
partition <- c(1, 2, 2, 3, 3, 3)

cat("\nPartition :", paste(partition, collapse = ", "), "\n")
sizes <- as.integer(table(partition))
cat("Tailles des groupes (par ordre de label):", paste(sizes, collapse = ", "), "\n")
cat("Nombre d'acteurs (N1):", length(partition), " | Nombre de groupes:", length(sizes), "\n\n")

# ======================================================================================
# 1) STAT OBSERVÉE SANS FIT
#    dry-run erpm -> extraction de la formule -> summary(formule, ~b1part)
# ======================================================================================

dry <- erpm(partition ~ squared_sizes(), eval_call = FALSE, verbose = TRUE)
fml <- dry[[2]]  # formule: nw ~ squared_sizes()

# Le réseau biparti 'nw' est dans l’environnement de la formule
obs <- summary(fml, constraints = ~ b1part)
cat("[MWE squared_sizes] summary observé =", obs, "\n\n")
stopifnot(is.numeric(obs), length(obs) == 1L, obs == 14)

# ======================================================================================
# 2) FIT ERGM COMPLET (ILLUSTRATIF)
#    Fit via erpm(), puis affichage standard. La stat observée du fit reste déterministe.
# ======================================================================================

set.seed(1) # stabilise l’estimation si on utilise une estimation CD dans ergm
fit <- erpm(partition ~ squared_sizes(),
            eval_call   = TRUE,
            verbose     = TRUE,
            estimate    = "MLE",     # MCMLE pour obtenir SE/logLik
            eval.loglik = TRUE,
            control     = list(MCMLE.maxit = 20))

# Vérification rapide de cohérence sur la formule du fit
obs_from_fit <- summary(fit$formula, constraints = ~ b1part)
stopifnot(obs_from_fit == 14)

cat("\n--- summary(fit) (style ergm) ---\n")
print(summary(fit))

cat("\n[OK] MWE squared_sizes: summary et fit vérifiés.\n")