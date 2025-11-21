# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cov_fullmatch.R
# Objet   : MWE pour l’effet ERPM `cov_fullmatch`
# Chaîne  : partition -> biparti (wrapper) -> summary(nw ~ ...) / erpm(partition ~ ...)
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  library(devtools)  # expose erpm(), build_bipartite_from_inputs, etc.
  library(ergm)
})

# Charge le package local
devtools::load_all(".")

# Patch {ergm} si présent
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R"); ergm_patch_enable()
}

# ----------------------------- Données fixes -----------------------------------
# Partition: 9 acteurs, groupes {1,1,2,2,3,3,3,4,4} -> tailles = (2,2,3,2)
partition <- c(1,1, 2,2, 3,3,3, 4,4)
labels    <- paste0("A", seq_along(partition))

# Attribut catégoriel (dept)
dept  <- c("A","B", "A","A", "B","B","B", "A","A")
nodes <- data.frame(label = labels, dept = dept, stringsAsFactors = FALSE)

cat("\nPartition : ", paste(partition, collapse=" "), "\n", sep = "")
cat("Tailles groupes : ", paste(as.integer(table(partition))), "\n\n", sep = "")

# ---------------------- Biparti via wrapper (référence) ------------------------
# 1) demandé : construire le biparti avec build_bipartite_from_inputs()
bld <- build_bipartite_from_inputs(partition = partition, nodes = nodes)
nw  <- bld$network

# -------------------------- Summary de référence (direct) ----------------------
# Définition : nombre de groupes homogènes sur la covariée
summary_ref_allSizes <- as.numeric(summary(nw ~ cov_fullmatch("dept"), constraints = ~ b1part))
summary_ref_A_2to3   <- as.numeric(summary(nw ~ cov_fullmatch("dept", category = "A", size = 2:3),
                                           constraints = ~ b1part))

# ---------------------- Observés via dry-run erpm() ----------------------------
# 2) demandé : noms explicites + appels directs
dry_allSizes <- erpm(partition ~ cov_fullmatch("dept"),
                     eval_call = FALSE, verbose = FALSE, nodes = nodes)
summary_obs_allSizes <- as.numeric(summary(dry_allSizes[[2]], constraints = ~ b1part))

dry_A_2to3 <- erpm(partition ~ cov_fullmatch("dept", category = "A", size = 2:3),
                   eval_call = FALSE, verbose = FALSE, nodes = nodes)
summary_obs_A_2to3 <- as.numeric(summary(dry_A_2to3[[2]], constraints = ~ b1part))

cat(sprintf("[summary] cov_fullmatch(dept) : observé=%g | référence=%g\n",
            summary_obs_allSizes, summary_ref_allSizes))
stopifnot(all.equal(summary_obs_allSizes, summary_ref_allSizes, tol = 0))

cat(sprintf("[summary] cov_fullmatch(dept, A, 2:3) : observé=%g | référence=%g\n",
            summary_obs_A_2to3, summary_ref_A_2to3))
stopifnot(all.equal(summary_obs_A_2to3, summary_ref_A_2to3, tol = 0))

# ------------------------------- Fit court ---------------------------------
# 3) demandé : eval.loglik=TRUE
set.seed(1)
nw <- build_bipartite_from_inputs(partition = partition, nodes = nodes)
fit <- ergm(
  nw$network ~ cov_fullmatch("dept", category = "A", size = 1:3),
  constraints = ~ b1part,
  # estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE
)

set.seed(1)
fit_mle <- erpm(
  partition ~ cov_fullmatch("dept", category = "A", size = 1:3),
  # estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE,
  nodes       = nodes
)

cat("\n--- summary(fit_mle) ---\n")
print(summary(fit_mle))

# Nettoyage patch
on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)