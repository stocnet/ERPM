# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cov_fulldiff.R
# Objet   : MWE minimal pour l’effet ERPM `cov_fulldiff`
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
# Partition simple (24 acteurs, 4 groupes de tailles 5,5,6,8)
partition <- c(
  rep(1, 5),
  rep(2, 5),
  rep(3, 6),
  rep(4, 8)
)
labels <- paste0("A", seq_along(partition))

# Attribut numérique : score
# g1 : variation faible, g2 : variation modérée,
# g3 : rampe, g4 : niveau élevé et étalé
score <- c(
  10, 11,  9, 10, 10,                 # g1
   5,  7, 12,  9,  6,                 # g2
   0,  1,  3,  2,  4,  6,             # g3
  50, 55, 60, 65, 70, 75, 80, 72      # g4
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
# cas 1 : toutes tailles >= 2
ref_all <- as.numeric(
  summary(nw ~ cov_fulldiff("score"),
          constraints = ~ b1part)
)

# cas 2 : tailles {5,6,8} (tous les groupes dans ce MWE, pour illustrer size comme filtre)
ref_mid <- as.numeric(
  summary(nw ~ cov_fulldiff("score", size = c(5,6,8)),
          constraints = ~ b1part)
)

# cas 3 : tailles >= 6 uniquement (groupes 3 et 4)
ref_big <- as.numeric(
  summary(nw ~ cov_fulldiff("score", size = 6:20),
          constraints = ~ b1part)
)

# -------------------------- Summary via erpm() --------------------------------
# mêmes trois cas, en passant par la traduction ERPM → ERGM

dry_all <- erpm(
  partition ~ cov_fulldiff("score"),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_all <- as.numeric(
  summary(dry_all[[2]], constraints = ~ b1part)
)

dry_mid <- erpm(
  partition ~ cov_fulldiff("score", size = c(5,6,8)),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_mid <- as.numeric(
  summary(dry_mid[[2]], constraints = ~ b1part)
)

dry_big <- erpm(
  partition ~ cov_fulldiff("score", size = 6:20),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_big <- as.numeric(
  summary(dry_big[[2]], constraints = ~ b1part)
)

# -------------------------- Vérifications --------------------------------------
cat(sprintf("[summary] cov_fulldiff(score)            : obs=%g | ref=%g\n",
            obs_all, ref_all))
stopifnot(all.equal(obs_all, ref_all, tol = 0))

cat(sprintf("[summary] cov_fulldiff(score,S={5,6,8})  : obs=%g | ref=%g\n",
            obs_mid, ref_mid))
stopifnot(all.equal(obs_mid, ref_mid, tol = 0))

cat(sprintf("[summary] cov_fulldiff(score,S=6:20)     : obs=%g | ref=%g\n",
            obs_big, ref_big))
stopifnot(all.equal(obs_big, ref_big, tol = 0))

# ------------------------------- Fit simple --------------------------------
set.seed(1)
nw2 <- build_bipartite_from_inputs(partition = partition, nodes = nodes)

fit_ergm <- ergm(
  nw2$network ~ cov_fulldiff("score"),
  constraints = ~ b1part,
  # estimate    = "MLE",
  # control     = control.ergm(
  #   MCMLE.maxit     = 40,
  #   MCMC.samplesize = 30000
  # ),
  eval.loglik = TRUE
)

set.seed(1)
fit_erpm <- erpm(
  partition ~ cov_fulldiff("score"),
  # estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE,
  # control     = control.ergm(
  #   MCMLE.maxit     = 40,
  #   MCMC.samplesize = 30000
  # ),
  nodes       = nodes
)

cat("\n--- summary(fit_erpm) ---\n")
print(summary(fit_erpm))

cat("\n--- summary(fit_ergm) ---\n")
print(summary(fit_ergm))

if (exists("ergm_patch_disable")) ergm_patch_disable()