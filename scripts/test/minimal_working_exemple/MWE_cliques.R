# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cliques.R
#
# Objet   : MWE pour l’effet ERPM `cliques(k)` — calcule la stat observée (summary)
#            via une formule directe sur le biparti + réalise un fit ERGM via `erpm()`.
#
# Chaîne ERPM → ERGM :
#   partition -> biparti (build_bipartite_from_inputs) -> summary(nw ~ ...) / erpm(partition ~ ...)
#
# Résumé technique :
#   • `cliques(k)` (k ≥ 1 ; k=1 = nombre de groupes de taille 1).
#   • `normalized=TRUE` divise par C(N1,k), N1 = nb d’acteurs.
#   • Summary : formule directe `nw ~ cliques(...)`, contrainte `~ b1part`.
#   • Fit : appel direct `erpm(partition ~ cliques(...), eval.loglik=TRUE)`.
# ==============================================================================

# ----- Préambule locale/UTF-8 -------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)   # 4) demandé

suppressPackageStartupMessages({
  library(devtools)   # expose erpm(), cliques(), build_bipartite_from_inputs, etc.
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

# ----- Construction explicite du biparti via le wrapper ERPM ------------------
# 1) demandé : construire le biparti depuis build_bipartite_from_inputs
bld <- build_bipartite_from_inputs(partition = partition)
nw  <- bld$network

# Petites fonctions de référence analytiques pour contrôle
expected_cliques_k1 <- function(part) {
  sz <- as.integer(table(part))
  sum(sz == 1L)
}
expected_cliques_k1_norm <- function(part) {
  num <- expected_cliques_k1(part)
  num / length(part)
}

# ==============================================================================
# 1) SUMMARY DIRECT — cas k = 1 et k = 2
# ==============================================================================

# ----- Cas k = 1 : doit renvoyer le nombre de groupes de taille 1 -------------
k1 <- 1L

stat_summary_cliques_k1 <- summary(
  nw ~ cliques(clique_size = k1),
  constraints = ~ b1part
)
stat_summary_cliques_k1_norm <- summary(
  nw ~ cliques(clique_size = k1, normalized = TRUE),
  constraints = ~ b1part
)

truth_k1      <- expected_cliques_k1(partition)
truth_k1_norm <- expected_cliques_k1_norm(partition)

cat(sprintf(
  "[summary] cliques(k=%d) (brut)       = %s | attendu = %s\n",
  k1,
  paste0(as.numeric(stat_summary_cliques_k1), collapse = ", "),
  truth_k1
))
cat(sprintf(
  "[summary] cliques(k=%d, normalized) = %s | attendu = %s\n\n",
  k1,
  paste0(as.numeric(stat_summary_cliques_k1_norm), collapse = ", "),
  format(truth_k1_norm)
))

stopifnot(
  length(stat_summary_cliques_k1)      == 1L,
  is.finite(stat_summary_cliques_k1),
  length(stat_summary_cliques_k1_norm) == 1L,
  is.finite(stat_summary_cliques_k1_norm)
)
stopifnot(
  isTRUE(all.equal(as.numeric(stat_summary_cliques_k1),
                   as.numeric(truth_k1),
                   tolerance = 1e-10))
)
stopifnot(
  isTRUE(all.equal(as.numeric(stat_summary_cliques_k1_norm),
                   as.numeric(truth_k1_norm),
                   tolerance = 1e-10))
)

# ----- Cas k = 2 : comme avant ------------------------------------------------
k <- 2L

# Summary non normalisé
stat_summary_cliques_k2 <- summary(
  nw ~ cliques(clique_size = k),
  constraints = ~ b1part
)
cat(sprintf("[summary] cliques(k=%d) = %s\n",
            k, paste0(as.numeric(stat_summary_cliques_k2), collapse = ", ")))

# Summary normalisé
stat_summary_cliques_k2_norm <- summary(
  nw ~ cliques(clique_size = k, normalized = TRUE),
  constraints = ~ b1part
)
cat(sprintf("[summary] cliques(k=%d, normalized=TRUE) = %s\n\n",
            k, paste0(as.numeric(stat_summary_cliques_k2_norm), collapse = ", ")))

# Garde-fous simples
stopifnot(length(stat_summary_cliques_k2)      == 1L, is.finite(stat_summary_cliques_k2))
stopifnot(length(stat_summary_cliques_k2_norm) == 1L, is.finite(stat_summary_cliques_k2_norm))

# ==============================================================================
# 2) FIT ERPM logLik — appel direct et minimal
# ==============================================================================
set.seed(1)

# 3) demandé :  eval.loglik=TRUE
fit_cliques_k2 <- erpm(
  partition ~ cliques(clique_size = 2),
  # estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE
)

cat("\n--- summary(fit_cliques_k2) (style ergm) ---\n")
print(summary(fit_cliques_k2))

# ----- Désactivation du patch {ergm} si activé --------------------------------
on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)
