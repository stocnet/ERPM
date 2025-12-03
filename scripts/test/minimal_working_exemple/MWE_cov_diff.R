# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cov_diff.R
# Objet   : MWE minimal pour l’effet ERPM `cov_diff`
# Chaîne  : partition → biparti → summary(nw) / erpm(partition)
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  library(devtools)  # expose erpm(), build_bipartite_from_inputs
  library(ergm)
})

# Charger le package local
devtools::load_all(".")

# Patch {ergm} si disponible
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# ----------------------------- Données fixes -----------------------------------
# Partition simple (15 acteurs, 3 groupes de tailles 4,5,6)
partition <- c(
  rep(1, 4),
  rep(2, 5),
  rep(3, 6)
)
labels <- paste0("A", seq_along(partition))

# Attribut numérique : score
# g1 : variation faible, g2 : variation forte, g3 : rampe
score <- c(
  10, 11,  9, 10,                 # g1
   5,  7, 20,  9,  6,             # g2
   0,  2,  4,  6,  8, 10          # g3
)

nodes <- data.frame(
  label = labels,
  score = score,
  stringsAsFactors = FALSE
)

cat("\nPartition :  ", paste(partition, collapse = " "), "\n", sep = "")
cat("Tailles groupes : ", paste(as.integer(table(partition))), "\n\n", sep = "")

# ---------------------- Biparti via wrapper (référence) ------------------------
bld <- build_bipartite_from_inputs(partition = partition, nodes = nodes)
nw  <- bld$network

# -------------------------- Summary de référence -------------------------------
# k = 2 (non normalisé)
ref_k2_raw <- as.numeric(
  summary(nw ~ cov_diff("score", clique_size = 2),
          constraints = ~ b1part)
)

# k = 2 (normalisé par groupe : "by_group" / TRUE)
ref_k2_bygrp <- as.numeric(
  summary(nw ~ cov_diff("score", clique_size = 2, normalized = "by_group"),
          constraints = ~ b1part)
)

# k = 2 (normalisation globale par taille de groupe : "global")
ref_k2_glob <- as.numeric(
  summary(nw ~ cov_diff("score", clique_size = 2, normalized = "global"),
          constraints = ~ b1part)
)

# k = 3 (non normalisé)
ref_k3_raw <- as.numeric(
  summary(nw ~ cov_diff("score", clique_size = 3),
          constraints = ~ b1part)
)

# k = 3 (normalisé par groupe, pour vérifier k > 2)
ref_k3_bygrp <- as.numeric(
  summary(nw ~ cov_diff("score", clique_size = 3, normalized = "by_group"),
          constraints = ~ b1part)
)

# -------------------------- Summary via erpm() --------------------------------
# mêmes cas, en passant par la traduction ERPM → ERGM

# k = 2 (non normalisé)
dry_k2_raw <- erpm(
  partition ~ cov_diff("score", clique_size = 2),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_k2_raw <- as.numeric(
  summary(dry_k2_raw[[2]], constraints = ~ b1part)
)

# k = 2 (normalisé par groupe : TRUE doit se comporter comme "by_group")
dry_k2_bygrp <- erpm(
  partition ~ cov_diff("score", clique_size = 2, normalized = TRUE),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_k2_bygrp <- as.numeric(
  summary(dry_k2_bygrp[[2]], constraints = ~ b1part)
)

# k = 2 (normalisation globale explicite)
dry_k2_glob <- erpm(
  partition ~ cov_diff("score", clique_size = 2, normalized = "global"),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_k2_glob <- as.numeric(
  summary(dry_k2_glob[[2]], constraints = ~ b1part)
)

# k = 3 (non normalisé)
dry_k3_raw <- erpm(
  partition ~ cov_diff("score", clique_size = 3),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_k3_raw <- as.numeric(
  summary(dry_k3_raw[[2]], constraints = ~ b1part)
)

# k = 3 (normalisé par groupe)
dry_k3_bygrp <- erpm(
  partition ~ cov_diff("score", clique_size = 3, normalized = "by_group"),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_k3_bygrp <- as.numeric(
  summary(dry_k3_bygrp[[2]], constraints = ~ b1part)
)

# -------------------------- Vérifications --------------------------------------
cat(sprintf("[summary] cov_diff(score,k=2)                 : obs=%g | ref=%g\n",
            obs_k2_raw, ref_k2_raw))
stopifnot(all.equal(obs_k2_raw, ref_k2_raw, tol = 0))

cat(sprintf("[summary] cov_diff(score,k=2,by_group/TRUE)  : obs=%g | ref=%g\n",
            obs_k2_bygrp, ref_k2_bygrp))
stopifnot(all.equal(obs_k2_bygrp, ref_k2_bygrp, tol = 0))

cat(sprintf("[summary] cov_diff(score,k=2,global)         : obs=%g | ref=%g\n",
            obs_k2_glob, ref_k2_glob))
stopifnot(all.equal(obs_k2_glob, ref_k2_glob, tol = 0))

cat(sprintf("[summary] cov_diff(score,k=3)                 : obs=%g | ref=%g\n",
            obs_k3_raw, ref_k3_raw))
stopifnot(all.equal(obs_k3_raw, ref_k3_raw, tol = 0))

cat(sprintf("[summary] cov_diff(score,k=3,by_group)       : obs=%g | ref=%g\n",
            obs_k3_bygrp, ref_k3_bygrp))
stopifnot(all.equal(obs_k3_bygrp, ref_k3_bygrp, tol = 0))

# ------------------------------- Fit simple --------------------------------
# On ajoute squared_sizes() pour stabiliser la structure des groupes
set.seed(1)
nw2 <- build_bipartite_from_inputs(partition = partition, nodes = nodes)

ctrl_mle <- control.ergm(
  MCMLE.maxit     = 40,
  MCMC.samplesize = 30000
)

fit_ergm <- ergm(
  nw2$network ~ squared_sizes() + cov_diff("score", clique_size = 2),
  constraints = ~ b1part,
  # estimate    = "MLE",
  # control     = ctrl_mle,
  eval.loglik = TRUE
)

set.seed(1)
fit_erpm <- erpm(
  partition ~ squared_sizes() + cov_diff("score", clique_size = 2),
  # estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE,
  # control     = ctrl_mle,
  nodes       = nodes
)

cat("\n--- summary(fit_erpm) ---\n")
print(summary(fit_erpm))

cat("\n--- summary(fit_ergm) ---\n")
print(summary(fit_ergm))

if (exists("ergm_patch_disable")) ergm_patch_disable()
