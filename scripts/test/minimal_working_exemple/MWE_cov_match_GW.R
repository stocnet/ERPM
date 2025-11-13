# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cov_match_GW.R
# Objet   : MWE minimal pour l’effet ERPM `cov_match_GW`
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
# Partition simple (10 acteurs, 4 groupes)
partition <- c(1,1, 2,2, 3,3,3, 4,4,4)
labels    <- paste0("A", seq_along(partition))

# Attributs : sexe, dept
sexe <- c("F","H",
          "F","F",
          "H","H","H",
          "F","H","H")

dept <- c("RH","RH",
          "Info","Info",
          "RH","Info","Info",
          "RH","RH","Info")

nodes <- data.frame(
  label = labels,
  sexe  = sexe,
  dept  = dept,
  stringsAsFactors = FALSE
)

cat("\nPartition :  ", paste(partition, collapse = " "), "\n", sep = "")
cat("Tailles groupes : ", paste(as.integer(table(partition))), "\n\n", sep = "")

# ---------------------- Biparti via wrapper (référence) ------------------------
bld <- build_bipartite_from_inputs(partition = partition, nodes = nodes)
nw  <- bld$network

# -------------------------- Summary de référence -------------------------------
# cas 1 : lambda = 2, tout sexe
ref_l2 <- as.numeric(
  summary(nw ~ cov_match_GW("sexe", lambda = 2),
          constraints = ~ b1part)
)

# cas 2 : lambda = 2, normalisation par groupe
ref_l2_bg <- as.numeric(
  summary(nw ~ cov_match_GW("sexe", lambda = 2, normalized = "by_group"),
          constraints = ~ b1part)
)

# cas 3 : lambda = 2, catégorie cible F
ref_l2_F <- as.numeric(
  summary(nw ~ cov_match_GW("sexe", lambda = 2, category = "F"),
          constraints = ~ b1part)
)

# -------------------------- Summary via erpm() --------------------------------
# même trois cas, en passant par la traduction ERPM → ERGM

dry_l2 <- erpm(
  partition ~ cov_match_GW("sexe", lambda = 2),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_l2 <- as.numeric(
  summary(dry_l2[[2]], constraints = ~ b1part)
)

dry_l2_bg <- erpm(
  partition ~ cov_match_GW("sexe", lambda = 2, normalized = "by_group"),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_l2_bg <- as.numeric(
  summary(dry_l2_bg[[2]], constraints = ~ b1part)
)

dry_l2_F <- erpm(
  partition ~ cov_match_GW("sexe", lambda = 2, category = "F"),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs_l2_F <- as.numeric(
  summary(dry_l2_F[[2]], constraints = ~ b1part)
)

# -------------------------- Vérifications --------------------------------------
cat(sprintf("[summary] cov_match_GW(sexe,λ=2)           : obs=%g | ref=%g\n",
            obs_l2, ref_l2))
stopifnot(all.equal(obs_l2, ref_l2, tol = 0))

cat(sprintf("[summary] cov_match_GW(sexe,λ=2,by_group) : obs=%g | ref=%g\n",
            obs_l2_bg, ref_l2_bg))
stopifnot(all.equal(obs_l2_bg, ref_l2_bg, tol = 0))

cat(sprintf("[summary] cov_match_GW(sexe==F,λ=2)       : obs=%g | ref=%g\n",
            obs_l2_F, ref_l2_F))
stopifnot(all.equal(obs_l2_F, ref_l2_F, tol = 0))

# ------------------------------- Fit MLE simple --------------------------------
set.seed(1)
nw2 <- build_bipartite_from_inputs(partition = partition, nodes = nodes)

fit_ergm <- ergm(
  nw2$network ~ cov_match_GW("sexe", lambda = 2),
  constraints = ~ b1part,
  estimate    = "MLE",
  control     = control.ergm(
    MCMLE.maxit     = 40,
    MCMC.samplesize = 30000
  ),
  eval.loglik = TRUE
)

set.seed(1)
fit_erpm <- erpm(
  partition ~ cov_match_GW("sexe", lambda = 2),
  estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE,
  control     = control.ergm(
    MCMLE.maxit     = 40,
    MCMC.samplesize = 30000
  ),
  nodes       = nodes
)

cat("\n--- summary(fit_erpm) ---\n")
print(summary(fit_erpm))

cat("\n--- summary(fit_ergm) ---\n")
print(summary(fit_ergm))

ergm_patch_disable()