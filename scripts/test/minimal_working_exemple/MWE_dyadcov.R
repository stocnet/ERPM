# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_dyadcov.R
# Objet   : MWE pour l’effet ERPM `dyadcov`
# Chaîne  : partition -> biparti ->
#           summary(nw ~ dyadcov("Z1", ...)) /
#           erpm(partition ~ dyadcov("Z1", ...))
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  library(devtools)
  library(ergm)
})

# Charge le package local ERPM
devtools::load_all(".")

# Patch {ergm} si présent
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# ----------------------------- Données fixes -----------------------------------
# Partition P1 : 15 acteurs, groupes {1,1,1,1,1,1,2,2,3,3,3,3,3,4,4}
partition <- c(1,1,1,1,1,1,
               2,2,
               3,3,3,3,3,
               4,4)
stopifnot(length(partition) == 15L)

labels <- paste0("N", seq_along(partition))
nodes  <- data.frame(label = labels, stringsAsFactors = FALSE)

# Covariée numérique -> matrice dyadique Z1 (15x15)
x  <- seq(0, by = 0.5, length.out = length(partition))
Z1 <- abs(outer(x, x, "-"))  # symétrique, diag = 0

cat("\nPartition : ", paste(partition, collapse = " "), "\n", sep = "")
cat("Tailles groupes : ", paste(as.integer(table(partition)), collapse = " "), "\n\n", sep = "")

cat("Z1[1:6,1:6] =\n")
print(round(Z1[1:6, 1:6], 3L))
cat("\n")

# ---------------------- Biparti via builder (référence) ------------------------
bld_ref <- build_bipartite_from_inputs(
  partition = partition,
  nodes     = nodes,
  dyads     = list(Z1 = Z1)
)
nw_ref <- bld_ref$network

# -------------------------- Summary de référence (direct) ----------------------
# Deux variantes illustratives :
#  - k = 2, non normalisé
#  - k = 2, normalisé
summary_ref_k2_nonnorm <- as.numeric(
  summary(nw_ref ~ dyadcov("Z1", clique_size = 2, normalized = FALSE),
          constraints = ~ b1part)
)

summary_ref_k2_norm <- as.numeric(
  summary(nw_ref ~ dyadcov("Z1", clique_size = 2, normalized = TRUE),
          constraints = ~ b1part)
)

# ---------------------- Observés via dry-run erpm() ----------------------------
# On laisse le wrapper reconstruire son biparti, mais on lui repasse Z1

dry_k2_nonnorm <- erpm(
  partition ~ dyadcov("Z1", clique_size = 2, normalized = FALSE),
  eval_call = FALSE,
  verbose   = FALSE,
  nodes     = nodes,
  dyads     = list(Z1 = Z1)
)
summary_obs_k2_nonnorm <- as.numeric(
  summary(dry_k2_nonnorm[[2]], constraints = ~ b1part)
)

dry_k2_norm <- erpm(
  partition ~ dyadcov("Z1", clique_size = 2, normalized = TRUE),
  eval_call = FALSE,
  verbose   = FALSE,
  nodes     = nodes,
  dyads     = list(Z1 = Z1)
)
summary_obs_k2_norm <- as.numeric(
  summary(dry_k2_norm[[2]], constraints = ~ b1part)
)

cat(sprintf("[summary] dyadcov('Z1', k=2, norm=FALSE) : observé=%g | référence=%g\n",
            summary_obs_k2_nonnorm, summary_ref_k2_nonnorm))
stopifnot(all.equal(summary_obs_k2_nonnorm, summary_ref_k2_nonnorm, tol = 0))

cat(sprintf("[summary] dyadcov('Z1', k=2, norm=TRUE)  : observé=%g | référence=%g\n",
            summary_obs_k2_norm, summary_ref_k2_norm))
stopifnot(all.equal(summary_obs_k2_norm, summary_ref_k2_norm, tol = 0))

# ------------------------------- Fit MLE court ---------------------------------
# On illustre un fit MLE avec la variante k=2 normalisée,
# qui est généralement mieux conditionnée sur ces exemples.

set.seed(1)
bld_fit <- build_bipartite_from_inputs(
  partition = partition,
  nodes     = nodes,
  dyads     = list(Z1 = Z1)
)
nw_fit <- bld_fit$network

cat("\n=== Fit de référence : ergm(nw ~ dyadcov('Z1', k=2, norm=TRUE)) ===\n")
fit_ref <- ergm(
  nw_fit ~ dyadcov("Z1", clique_size = 2, normalized = TRUE),
  constraints = ~ b1part,
#   estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE
)

set.seed(1)
cat("\n=== Fit via erpm(partition ~ dyadcov('Z1', k=2, norm=TRUE)) ===\n")
fit_erpm <- erpm(
  partition ~ dyadcov("Z1", clique_size = 2, normalized = TRUE),
#   estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE,
  nodes       = nodes,
  dyads       = list(Z1 = Z1)
)

cat("\n--- summary(fit_ref) ---\n")
print(summary(fit_ref))

cat("\n--- summary(fit_erpm) ---\n")
print(summary(fit_erpm))

# Nettoyage patch
on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)