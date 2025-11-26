# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_dyadcov_GW.R
# Objet   : MWE pour l’effet ERPM `dyadcov_GW`
# Chaîne  : partition -> biparti ->
#           summary(nw ~ dyadcov_GW("Zk", lambda=...)) /
#           erpm(partition ~ dyadcov_GW("Zk", lambda=...))
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  library(devtools)
  library(ergm)
})

# ----------------------------------------------------------------------
# Chargement du package local ERPM
# ----------------------------------------------------------------------
devtools::load_all(".")

# Patch {ergm} si présent
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# ==============================================================================
# 1. Données fixes : partition + matrices dyadiques Z1 / Z2
# ==============================================================================

# Partition P1 : 5 acteurs, groupes {1,1,2,2,3}
partition <- c(1L, 1L, 2L, 2L, 3L)
labels    <- paste0("N", seq_along(partition))
nodes     <- data.frame(label = labels, stringsAsFactors = FALSE)

cat("\nPartition : ", paste(partition, collapse = " "), "\n", sep = "")
cat("Tailles groupes : ", paste(as.integer(table(partition)), collapse = " "), "\n\n", sep = "")

# Matrice dyadique Z1 (5x5) - symétrique, diag = 0
Z1 <- matrix(
  c(
    0.0, 1.0, 0.5, 0.3, 0.8,
    1.0, 0.0, 1.2, 0.4, 0.2,
    0.5, 1.2, 0.0, 0.9, 0.6,
    0.3, 0.4, 0.9, 0.0, 1.1,
    0.8, 0.2, 0.6, 1.1, 0.0
  ),
  nrow = 5, ncol = 5, byrow = TRUE
)

# Matrice dyadique Z2 (5x5) - éventuellement non symétrique, diag = 0
Z2 <- matrix(
  c(
    0.0, 0.3, 0.7, 1.2, 0.5,
    0.4, 0.0, 0.6, 0.9, 1.5,
    1.1, 0.2, 0.0, 0.8, 0.3,
    0.9, 1.4, 0.5, 0.0, 1.0,
    0.2, 0.6, 1.3, 0.7, 0.0
  ),
  nrow = 5, ncol = 5, byrow = TRUE
)

stopifnot(all(diag(Z1) == 0), all(diag(Z2) == 0))

cat("Z1[1:5,1:5] =\n")
print(Z1)
cat("\nZ2[1:5,1:5] =\n")
print(Z2)
cat("\n")

# ==============================================================================
# 2. Biparti via builder + summary de référence
# ==============================================================================

bld_ref <- build_bipartite_from_inputs(
  partition = partition,
  nodes     = nodes,
  dyads     = list(Z1 = Z1, Z2 = Z2)
)
nw_ref <- bld_ref$network

# Helper : summary direct via ergm
summary_direct <- function(rhs_str) {
  f <- as.formula(paste0("nw_ref ~ ", rhs_str))
  environment(f) <- list2env(list(nw_ref = nw_ref), parent = parent.frame())
  as.numeric(summary(f, constraints = ~ b1part))
}

# Helper : summary via erpm(eval_call = FALSE), à la manière du selftest
summary_via_erpm <- function(rhs_str) {
  # Formule côté ERPM
  partition_local <- partition
  f_erpm <- as.formula(paste0("partition_local ~ ", rhs_str))
  environment(f_erpm) <- list2env(
    list(partition_local = partition, nodes = nodes),
    parent = parent.frame()
  )

  call_ergm <- erpm(
    f_erpm,
    eval_call = FALSE,
    verbose   = FALSE,
    nodes     = nodes,
    dyads     = list(Z1 = Z1, Z2 = Z2)
  )

  # Extraire la RHS de la formule traduite
  ergm_form <- call_ergm[[2L]]
  rhs_expr  <- ergm_form[[3L]]

  # Reconstruire un réseau via le builder
  bld2 <- build_bipartite_from_inputs(
    partition = partition,
    nodes     = nodes,
    dyads     = list(Z1 = Z1, Z2 = Z2)
  )
  nw2 <- bld2$network

  # Construire une formule nw2 ~ <rhs_expr> et appeler summary()
  f2 <- as.formula(bquote(nw2 ~ .(rhs_expr)))
  environment(f2) <- list2env(list(nw2 = nw2), parent = parent.frame())

  as.numeric(summary(f2, constraints = ~ b1part))
}

# Panel minimal de cas pour dyadcov_GW
cases <- c(
  "dyadcov_GW('Z1', lambda = 2)",
  "dyadcov_GW('Z1', lambda = 3)",
  "dyadcov_GW('Z2', lambda = 2)"
)

cat("=== PHASE 1 : Summary(nw) vs Summary(erpm-dry-run) [dyadcov_GW] ===\n")

for (rhs_str in cases) {
  s_ref  <- summary_direct(rhs_str)
  s_erpm <- summary_via_erpm(rhs_str)

  cat(sprintf("[summary] %s : observé=%g | référence=%g\n",
              rhs_str, s_erpm, s_ref))
  stopifnot(all.equal(s_erpm, s_ref, tol = 0))
}

# ==============================================================================
# 3. Fits courts : ergm direct vs erpm
# ==============================================================================

cat("\n=== PHASE 2 : Fits courts (ergm vs erpm) [dyadcov_GW] ===\n")

set.seed(1)
fit_ref_Z1_l2 <- ergm(
  nw_ref ~ dyadcov_GW("Z1", lambda = 2),
  constraints = ~ b1part,
  eval.loglik = TRUE,
  verbose     = FALSE
)

set.seed(1)
fit_erpm_Z1_l2 <- erpm(
  partition ~ dyadcov_GW("Z1", lambda = 2),
  eval.loglik = TRUE,
  verbose     = FALSE,
  nodes       = nodes,
  dyads       = list(Z1 = Z1, Z2 = Z2)
)

set.seed(1)
fit_ref_Z2_l2 <- ergm(
  nw_ref ~ dyadcov_GW("Z2", lambda = 2),
  constraints = ~ b1part,
  eval.loglik = TRUE,
  verbose     = FALSE
)

set.seed(1)
fit_erpm_Z2_l2 <- erpm(
  partition ~ dyadcov_GW("Z2", lambda = 2),
  eval.loglik = TRUE,
  verbose     = FALSE,
  nodes       = nodes,
  dyads       = list(Z1 = Z1, Z2 = Z2)
)

cat("\n--- summary(fit_ref_Z1_l2) ---\n")
print(summary(fit_ref_Z1_l2))

cat("\n--- summary(fit_erpm_Z1_l2) ---\n")
print(summary(fit_erpm_Z1_l2))

cat("\n--- summary(fit_ref_Z2_l2) ---\n")
print(summary(fit_ref_Z2_l2))

cat("\n--- summary(fit_erpm_Z2_l2) ---\n")
print(summary(fit_erpm_Z2_l2))

cat("\nMWE dyadcov_GW terminé sans erreur.\n")
on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)