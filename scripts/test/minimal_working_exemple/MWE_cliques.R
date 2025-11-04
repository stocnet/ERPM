# ======================================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cliques.R
#
# Objet   : MWE pour l’effet ERPM `cliques(k)` — calcule la stat observée (summary)
#            via dry-run `erpm()` + calcule un fit ERGM via `erpm()` non normalisé.
#
# Contexte: Ce script illustre la chaîne ERPM → ERGM :
#            partition -> réseau biparti -> traduction -> (dry) formule -> summary(formule) / ergm()
#            partition -> réseau biparti -> traduction -> formule -> ergm()
#
# Résumé technique :
#   • `cliques(k)` est implémenté côté ERGM (InitErgmTerm + change-stat C).
#   • Pour le summary : `erpm(..., eval_call=FALSE)` renvoie un appel `ergm(formule, ...)` non évalué.
#     La formule capture `nw` dans son environnement ; `summary(formule, ...)` l’utilise.
#   • k ≥ 1 pris en charge ; k=1 = nombre de groupes de taille 1.
#   • `normalized=TRUE` divise par C(N1,k), N1 = nb d’acteurs.
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
#   Groupes effectifs = {1,2,3,4} ; tailles = (1,4,2,3) ; N1 = 10 acteurs
partition <- c(1, 2, 2, 2, 2, 3, 3, 4, 4, 4)

cat("\nPartition :", paste(partition, collapse = ", "), "\n")
sizes <- as.integer(table(partition))
cat("Tailles des groupes (par ordre de label):", paste(sizes, collapse = ", "), "\n")
cat("Nombre d'acteurs (N1):", length(partition), " | Nombre de groupes:", length(sizes), "\n\n")

# ======================================================================================
# 1) STAT OBSERVÉE SANS FIT — 
#    dry-run erpm -> extraction de la formule -> summary(formule, ~b1part)
# ======================================================================================
k <- 2
dry <- erpm(partition ~ cliques(cliques_sizes = k), eval_call = FALSE, verbose = TRUE)
fml <- dry[[2]]  # formule: nw ~ cliques(clique_size=k, normalized=FALSE)

# Le réseau biparti 'nw' est dans l’environnement de la formule
obs <- summary(fml, constraints = ~ b1part)
cat(sprintf("[MWE cliques(k=%d)] summary observé = %s  | terme = %s\n\n",
            k, paste0(obs, collapse = ", "), paste0(names(obs), collapse = ", ")))

stopifnot(is.numeric(obs), length(obs) == 1L)

# ======================================================================================
# 2) STAT OBSERVÉE NORMALISÉE —
#    dry-run erpm avec normalized=TRUE, puis summary() sur la formule
# ======================================================================================

dry_norm <- erpm(partition ~ cliques(cliques_sizes = 2, normalized = TRUE),
                 eval_call = FALSE, verbose = TRUE)
fml_norm <- dry_norm[[2]]

obs_norm <- summary(fml_norm, constraints = ~ b1part)
cat("[MWE cliques(k=2, normalized=TRUE)] summary observé =", obs_norm, "\n\n")

stopifnot(is.numeric(obs_norm), length(obs_norm) == 1L)

# ======================================================================================
# 3) FIT ERGM COMPLET NON NORMALISÉ — k = 2
#    Un fit via erpm() est réalisé, puis un résumé standard du fit est affiché.
# ======================================================================================

set.seed(1)  # stabilise l’estimation si on utilise une estimation CD dans ergm
fit <- erpm(partition ~ cliques(cliques_sizes = 2),
            eval_call   = TRUE,
            verbose     = TRUE,
            estimate    = "MLE",
            eval.loglik = TRUE)

# Affiche un résumé du fit
cat("\n--- summary(fit) (style ergm) ---\n")
print(summary(fit))