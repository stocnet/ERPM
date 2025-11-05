# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_log_factorial_sizes.R
#
# Objet   : MWE pour l’effet ERPM `log_factorial_sizes` — calcule la stat observée
#            via dry-run `erpm()` + réalise un fit ERGM via `erpm()`.
#
# Chaîne : partition -> biparti -> traduction -> (dry) formule -> summary(formule)
#          partition -> biparti -> traduction -> formule -> ergm()
#
# Expression math  : stat = Σ_g log((n_g - 1)!) = Σ_g lgamma(n_g)
# Remarque : pas de paramètre pour cet effet. Contrainte ~b1part.
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
# Référence: Σ_g log((n_g - 1)!) = Σ_g lgamma(n_g)
stat_expected <- function(part) {
  sizes <- as.integer(table(part))
  sum(lgamma(sizes))
}

# ==============================================================================
# 1) STAT OBSERVÉE — dry-run erpm -> formule -> summary(formule, ~b1part)
# ==============================================================================

# Dry-run: partition LHS -> call ergm(...)
dry <- erpm(partition ~ log_factorial_sizes, eval_call = FALSE, verbose = TRUE)

# Extraire la formule {ergm} traduite et l’évaluer sous contrainte bipartite
fml <- dry[[2]]
obs <- summary(fml, constraints = ~ b1part)

# Comparaison à la valeur fermée
exp_val <- stat_expected(partition)

cat(sprintf("[MWE log_factorial_sizes] summary observé = %.12f | attendu = %.12f\n",
            as.numeric(obs), as.numeric(exp_val)))
stopifnot(is.numeric(obs), length(obs) == 1L,
          isTRUE(all.equal(as.numeric(obs), as.numeric(exp_val), tol = 1e-10)))

# ==============================================================================
# 1-bis) Explicitation (pas sûr que ça se dise) du réseau
# ==============================================================================
# Construction bipartie minimale si le helper projet n’est pas chargé
partition_to_bipartite_network_min <- function(labels, partition) {
  stopifnot(length(labels) == length(partition))
  nA <- length(partition)
  G  <- max(partition)
  inc <- matrix(0L, nrow = nA, ncol = G,
                dimnames = list(labels, paste0("G", seq_len(G))))
  inc[cbind(seq_len(nA), partition)] <- 1L
  nw <- network::network(inc, matrix.type = "bipartite", bipartite = nA, directed = FALSE)
  network::set.vertex.attribute(nw, "vertex.names", c(labels, colnames(inc)))
  nw
}

lbl <- paste0("A", seq_len(length(partition)))
nw_bip <- partition_to_bipartite_network_min(lbl, partition)

# Formule réseau explicite
obs_nw <- summary(nw_bip ~ log_factorial_sizes(), constraints = ~ b1part)
cat(sprintf("[Route réseau] summary(nw ~ log_factorial_sizes()) = %.12f\n", as.numeric(obs_nw)))
stopifnot(isTRUE(all.equal(as.numeric(obs_nw), as.numeric(exp_val), tol = 1e-10)))

# ==============================================================================
# 2) FIT ERGM COURT — estimation via erpm()
# ==============================================================================
set.seed(1)
oldopt <- options(ergm.loglik.warn_dyads = FALSE); on.exit(options(oldopt), add = TRUE)

fit <- erpm(partition ~ log_factorial_sizes,
            eval_call   = TRUE,
            verbose     = TRUE,
            estimate    = "CD",       # rapide et suffisant pour le MWE
            eval.loglik = FALSE,
            control     = list(MCMLE.maxit = 2, MCMC.samplesize = 500))

cat("\n--- summary(fit) (style ergm) ---\n")
print(summary(fit))

# ----- Nettoyage patch ---------------------------------------------------------
on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)