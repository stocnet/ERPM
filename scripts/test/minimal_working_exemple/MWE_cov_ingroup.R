# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cov_ingroup.R
# Objet   : MWE minimal pour l’effet ERPM `cov_ingroup`
# Chaîne  : partition -> biparti -> summary(nw~...) ; partition -> erpm(dry) ; partition -> erpm(fit)
# ==============================================================================

options(ergm.loglik.warn_dyads = FALSE)
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  library(devtools)
  library(ergm)
})

devtools::load_all(".")

source("scripts/ergm_patch.R")
ergm_patch_enable()

# ----------------------------- Données fixes -----------------------------------
partition <- c(1,1, 2,2, 3,3,3, 4,4)             # tailles = (2,2,3,2)
nodes <- data.frame(
  label = paste0("A", seq_along(partition)),
  age   = c(20,30, 40,50, 10,10,10, 5,5),
  dept  = c("A","B","A","C","A","A","B","C","A"),
  stringsAsFactors = FALSE
)

# -------------------- Biparti (wrapper) ----------------------------------------
nw <- build_bipartite_from_inputs(partition = partition, nodes = nodes)$network

# -------------------- SUMMARY de référence (direct) ----------------------------
exp_age_S23  <- as.numeric(summary(nw ~ cov_ingroup("age",  size = 2:3), constraints = ~ b1part))
exp_deptA_S3 <- as.numeric(summary(nw ~ cov_ingroup("dept", category = "A", size = 3), constraints = ~ b1part))
cat("[REF] age S{2,3} =", exp_age_S23, " | dept=='A' S{3} =", exp_deptA_S3, "\n")

# -------------------- SUMMARY via erpm(dry-run) --------------------------------
dry <- erpm(
  partition ~ cov_ingroup("age", size = 2:3) + cov_ingroup("dept", category = "A", size = 3),
  eval_call = FALSE, verbose = FALSE, nodes = nodes
)
obs <- as.numeric(summary(dry[[2]], constraints = ~ b1part))
cat("[DRY] stats =", paste(obs, collapse = ", "), "\n")

stopifnot(length(obs) == 2L, isTRUE(all.equal(obs, c(exp_age_S23, exp_deptA_S3), tol = 0)))

# -------------------- FIT (CD court et lisible) --------------------------------
set.seed(1)
fit <- erpm(
  partition ~ cov_ingroup("age", size = 2:3) + cov_ingroup("dept", category = "A", size = 3),
  estimate    = "CD",
  eval.loglik = TRUE,
  verbose     = TRUE,
  nodes       = nodes
)
print(summary(fit))

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)