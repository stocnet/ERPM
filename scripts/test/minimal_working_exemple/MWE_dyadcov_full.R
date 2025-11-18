# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_dyadcov_full.R
# Objet   : MWE pour l’effet ERPM `dyadcov_full`
# Chaîne  : partition -> biparti ->
#           summary(nw ~ dyadcov_full("Z1")) /
#           erpm(partition ~ dyadcov_full("Z1"))
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
# Partition P1 : 9 acteurs, groupes {1,1,1,1,2,3,3,3,4}
partition <- c(1,1,1,1, 2, 3,3,3, 4)
labels    <- paste0("N", seq_along(partition))
nodes     <- data.frame(label = labels, stringsAsFactors = FALSE)

# Covariée numérique -> matrice dyadique Z1 (9x9)
x  <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
Z1 <- abs(outer(x, x, "-"))  # symétrique, diag = 0

cat("\nPartition : ", paste(partition, collapse = " "), "\n", sep = "")
cat("Tailles groupes : ", paste(as.integer(table(partition)), collapse = " "), "\n\n", sep = "")

# ---------------------- Biparti via builder (référence) ------------------------
bld_ref <- build_bipartite_from_inputs(
  partition = partition,
  nodes     = nodes,
  dyads     = list(Z1 = Z1)
)
nw_ref <- bld_ref$network

# -------------------------- Summary de référence (direct) ----------------------
summary_ref_all <- as.numeric(
  summary(nw_ref ~ dyadcov_full("Z1"), constraints = ~ b1part)
)

summary_ref_2to3 <- as.numeric(
  summary(nw_ref ~ dyadcov_full("Z1", size = 2:3), constraints = ~ b1part)
)

# ---------------------- Observés via dry-run erpm() ----------------------------
# On laisse le wrapper reconstruire son biparti, mais on lui repasse Z1
dry_all <- erpm(
  partition ~ dyadcov_full("Z1"),
  eval_call = FALSE,
  verbose   = FALSE,
  nodes     = nodes,
  dyads     = list(Z1 = Z1)
)
summary_obs_all <- as.numeric(
  summary(dry_all[[2]], constraints = ~ b1part)
)

dry_2to3 <- erpm(
  partition ~ dyadcov_full("Z1", size = 2:3),
  eval_call = FALSE,
  verbose   = FALSE,
  nodes     = nodes,
  dyads     = list(Z1 = Z1)
)
summary_obs_2to3 <- as.numeric(
  summary(dry_2to3[[2]], constraints = ~ b1part)
)

cat(sprintf("[summary] dyadcov_full('Z1')           : observé=%g | référence=%g\n",
            summary_obs_all, summary_ref_all))
stopifnot(all.equal(summary_obs_all, summary_ref_all, tol = 0))

cat(sprintf("[summary] dyadcov_full('Z1', size=2:3) : observé=%g | référence=%g\n",
            summary_obs_2to3, summary_ref_2to3))
stopifnot(all.equal(summary_obs_2to3, summary_ref_2to3, tol = 0))

# ------------------------------- Fit MLE court ---------------------------------
set.seed(1)
bld_fit <- build_bipartite_from_inputs(
  partition = partition,
  nodes     = nodes,
  dyads     = list(Z1 = Z1)
)
nw_fit <- bld_fit$network

fit_ref <- ergm(
  nw_fit ~ dyadcov_full("Z1", size = 2:3),
  constraints = ~ b1part,
  estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE
)

set.seed(1)
fit_erpm <- erpm(
  partition ~ dyadcov_full("Z1", size = 2:3),
  estimate    = "MLE",
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