# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_dyadcov.R
# Objet   : MWE pour l’effet ERPM `dyadcov`
# Chaîne  : partition -> biparti ->
#           summary(nw ~ dyadcov("Z1", ...)) /
#           erpm(partition ~ dyadcov("Z1", ...))
#           + tests des nouvelles conventions (Z non symétrique, diag = 0)
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

# Covariée numérique -> matrice dyadique Z1 (15x15, symétrique, diag=0)
x  <- seq(0, by = 0.5, length.out = length(partition))
Z1 <- abs(outer(x, x, "-"))  # symétrique, diag = 0
diag(Z1) <- 0

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

# ==============================================================================
# 1) Summary de référence (Z symétrique) + dry-run erpm() (k = 2)
# ==============================================================================

summary_ref_k2_nonnorm <- as.numeric(
  summary(nw_ref ~ dyadcov("Z1", clique_size = 2, normalized = FALSE),
          constraints = ~ b1part)
)

summary_ref_k2_norm <- as.numeric(
  summary(nw_ref ~ dyadcov("Z1", clique_size = 2, normalized = TRUE),
          constraints = ~ b1part)
)

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

cat(sprintf("[Z1|summary] dyadcov('Z1', k=2, norm=FALSE) : observé=%g | référence=%g\n",
            summary_obs_k2_nonnorm, summary_ref_k2_nonnorm))
stopifnot(all.equal(summary_obs_k2_nonnorm, summary_ref_k2_nonnorm, tol = 0))

cat(sprintf("[Z1|summary] dyadcov('Z1', k=2, norm=TRUE)  : observé=%g | référence=%g\n",
            summary_obs_k2_norm, summary_ref_k2_norm))
stopifnot(all.equal(summary_obs_k2_norm, summary_ref_k2_norm, tol = 0))

# ==============================================================================
# 2) Nouvelles conventions : Z non symétrique (diag = 0) + symétrisation interne
# ==============================================================================

# Construction d'une matrice dyadique non symétrique Z2
# - base positive
# - asymétrie volontaire
# - diagonale forcée à 0 (contrainte du modèle)
base  <- abs(outer(x, rev(x), "-"))
Z2    <- base
Z2[lower.tri(Z2)] <- 2 * Z2[lower.tri(Z2)]   # asymétrie volontaire
diag(Z2) <- 0                                # convention: diag = 0

cat("Z2[1:6,1:6] (non symétrique, diag=0) =\n")
print(round(Z2[1:6, 1:6], 3L))
cat("\n")

# T^{(2)}(p;Z) = sum_g sum_{i<j in g} (z_ij+z_ji)
T2_manual <- function(part, Z) {
  stopifnot(length(part) == nrow(Z), nrow(Z) == ncol(Z))
  split_idx <- split(seq_along(part), part)
  total <- 0
  for (g in names(split_idx)) {
    idx <- split_idx[[g]]
    n_g <- length(idx)
    if (n_g < 2L) next
    for (ii in seq_len(n_g - 1L)) {
      i <- idx[ii]
      for (jj in (ii + 1L):n_g) {
        j <- idx[jj]
        total <- total + (Z[i, j] + Z[j, i])
      }
    }
  }
  total
}

# Variante normalisée : moyenne par groupe (division par choose(n_g, 2))
T2_manual_norm <- function(part, Z) {
  stopifnot(length(part) == nrow(Z), nrow(Z) == ncol(Z))
  split_idx <- split(seq_along(part), part)
  total <- 0
  for (g in names(split_idx)) {
    idx <- split_idx[[g]]
    n_g <- length(idx)
    if (n_g < 2L) next
    Sg <- 0
    for (ii in seq_len(n_g - 1L)) {
      i <- idx[ii]
      for (jj in (ii + 1L):n_g) {
        j <- idx[jj]
        Sg <- Sg + (Z[i, j] + Z[j, i])
      }
    }
    total <- total + Sg / choose(n_g, 2)
  }
  total
}

summary_ref_Z2_k2_nonnorm <- T2_manual(partition, Z2)
summary_ref_Z2_k2_norm    <- T2_manual_norm(partition, Z2)

bld_ref_Z2 <- build_bipartite_from_inputs(
  partition = partition,
  nodes     = nodes,
  dyads     = list(Z2 = Z2)
)
nw_ref_Z2 <- bld_ref_Z2$network

summary_ergm_Z2_k2_nonnorm <- as.numeric(
  summary(nw_ref_Z2 ~ dyadcov("Z2", clique_size = 2, normalized = FALSE),
          constraints = ~ b1part)
)

summary_ergm_Z2_k2_norm <- as.numeric(
  summary(nw_ref_Z2 ~ dyadcov("Z2", clique_size = 2, normalized = TRUE),
          constraints = ~ b1part)
)

dry_Z2_k2_nonnorm <- erpm(
  partition ~ dyadcov("Z2", clique_size = 2, normalized = FALSE),
  eval_call = FALSE,
  verbose   = FALSE,
  nodes     = nodes,
  dyads     = list(Z2 = Z2)
)
summary_erpm_Z2_k2_nonnorm <- as.numeric(
  summary(dry_Z2_k2_nonnorm[[2]], constraints = ~ b1part)
)

dry_Z2_k2_norm <- erpm(
  partition ~ dyadcov("Z2", clique_size = 2, normalized = TRUE),
  eval_call = FALSE,
  verbose   = FALSE,
  nodes     = nodes,
  dyads     = list(Z2 = Z2)
)
summary_erpm_Z2_k2_norm <- as.numeric(
  summary(dry_Z2_k2_norm[[2]], constraints = ~ b1part)
)

cat(sprintf("[Z2|manual vs ergm] dyadcov('Z2', k=2, norm=FALSE) : manual=%g | ergm=%g\n",
            summary_ref_Z2_k2_nonnorm, summary_ergm_Z2_k2_nonnorm))
stopifnot(all.equal(summary_ref_Z2_k2_nonnorm, summary_ergm_Z2_k2_nonnorm, tol = 1e-8))

cat(sprintf("[Z2|manual vs erpm] dyadcov('Z2', k=2, norm=FALSE) : manual=%g | erpm=%g\n",
            summary_ref_Z2_k2_nonnorm, summary_erpm_Z2_k2_nonnorm))
stopifnot(all.equal(summary_ref_Z2_k2_nonnorm, summary_erpm_Z2_k2_nonnorm, tol = 1e-8))

cat(sprintf("[Z2|manual vs ergm] dyadcov('Z2', k=2, norm=TRUE)  : manual=%g | ergm=%g\n",
            summary_ref_Z2_k2_norm, summary_ergm_Z2_k2_norm))
stopifnot(all.equal(summary_ref_Z2_k2_norm, summary_ergm_Z2_k2_norm, tol = 1e-8))

cat(sprintf("[Z2|manual vs erpm] dyadcov('Z2', k=2, norm=TRUE)  : manual=%g | erpm=%g\n",
            summary_ref_Z2_k2_norm, summary_erpm_Z2_k2_norm))
stopifnot(all.equal(summary_ref_Z2_k2_norm, summary_erpm_Z2_k2_norm, tol = 1e-8))

# ==============================================================================
# 3) Fit court (k = 2, normalisé) sur Z1
# ==============================================================================

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
  eval.loglik = TRUE,
  verbose     = TRUE
)

set.seed(1)
cat("\n=== Fit via erpm(partition ~ dyadcov('Z1', k=2, norm=TRUE)) ===\n")
fit_erpm <- erpm(
  partition ~ dyadcov("Z1", clique_size = 2, normalized = TRUE),
  eval.loglik = TRUE,
  verbose     = TRUE,
  nodes       = nodes,
  dyads       = list(Z1 = Z1)
)

cat("\n--- summary(fit_ref) ---\n")
print(summary(fit_ref))

cat("\n--- summary(fit_erpm) ---\n")
print(summary(fit_erpm))

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)