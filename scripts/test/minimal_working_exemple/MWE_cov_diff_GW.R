# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cov_diff_GW.R
# Objet   : MWE minimal pour l’effet ERPM `cov_diff_GW`
# Chaîne  : partition → biparti → summary(nw) / erpm(partition) → fit MLE simple
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
# lambda simples et vectorisés, pour illustrer la pondération géométrique

# lambda = 2
ref_l2 <- as.numeric(
  summary(nw ~ cov_diff_GW("score", lambda = 2),
          constraints = ~ b1part)
)

# lambda = 3
ref_l3 <- as.numeric(
  summary(nw ~ cov_diff_GW("score", lambda = 3),
          constraints = ~ b1part)
)

# lambda = c(2,4) (vectorisation)
ref_l2_l4 <- as.numeric(
  summary(nw ~ cov_diff_GW("score", lambda = c(2, 4)),
          constraints = ~ b1part)
)

# -------------------------- Summary via erpm() --------------------------------
# mêmes trois cas, en passant par la traduction ERPM → ERGM

# lambda = 2
dry_l2 <- erpm(
  partition ~ cov_diff_GW("score", lambda = 2),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_l2 <- as.numeric(
  summary(dry_l2[[2]], constraints = ~ b1part)
)

# lambda = 3
dry_l3 <- erpm(
  partition ~ cov_diff_GW("score", lambda = 3),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_l3 <- as.numeric(
  summary(dry_l3[[2]], constraints = ~ b1part)
)

# lambda = c(2,4)
dry_l2_l4 <- erpm(
  partition ~ cov_diff_GW("score", lambda = c(2, 4)),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_l2_l4 <- as.numeric(
  summary(dry_l2_l4[[2]], constraints = ~ b1part)
)

# -------------------------- Vérifications summary ------------------------------
cat(sprintf("[summary] cov_diff_GW(score,lambda=2)      : obs=%g | ref=%g\n",
            obs_l2, ref_l2))
stopifnot(all.equal(obs_l2, ref_l2, tol = 0))

cat(sprintf("[summary] cov_diff_GW(score,lambda=3)      : obs=%g | ref=%g\n",
            obs_l3, ref_l3))
stopifnot(all.equal(obs_l3, ref_l3, tol = 0))

cat(sprintf("[summary] cov_diff_GW(score,lambda=c(2,4)) : obs=%s | ref=%s\n",
            paste(obs_l2_l4, collapse = ","), paste(ref_l2_l4, collapse = ",")))
stopifnot(all.equal(obs_l2_l4, ref_l2_l4, tol = 0))

# ------------------------------- Fit MLE simple --------------------------------
# Modèle jouet : terme structurel + cov_diff_GW(lambda=2)
set.seed(1)
bld2 <- build_bipartite_from_inputs(partition = partition, nodes = nodes)
nw2  <- bld2$network

ctrl_mle <- control.ergm(
  MCMLE.maxit     = 40,
  MCMC.samplesize = 30000
)

# Fit direct sur le biparti explicite
fit_ergm <- ergm(
  nw2 ~ cov_diff_GW("score", lambda = 2),
  constraints = ~ b1part,
  estimate    = "MLE",
  control     = ctrl_mle,
  eval.loglik = TRUE
)

# Fit via erpm() et traduction automatique
set.seed(1)
fit_erpm <- erpm(
  partition ~ cov_diff_GW("score", lambda = 2),
  estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE,
  control     = ctrl_mle,
  nodes       = nodes
)

cat("\n--- summary(fit_erpm) ---\n")
print(summary(fit_erpm))

cat("\n--- summary(fit_ergm) ---\n")
print(summary(fit_ergm))

if (file.exists("scripts/ergm_patch.R") && exists("ergm_patch_disable")) {
  ergm_patch_disable()
}