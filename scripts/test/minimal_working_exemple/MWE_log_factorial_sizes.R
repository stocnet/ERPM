# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_log_factorial_sizes.R
# Objet   : MWE pour l’effet ERPM `log_factorial_sizes`
# Chaîne  : partition -> biparti (wrapper) -> summary de référence
#           partition -> erpm(dry-run) -> summary observé -> comparaison
#           partition -> erpm() -> summary(fit)
# Stat    : Σ_g log((n_g - 1)!) = Σ_g lgamma(n_g) ; contrainte ~ b1part
# ==============================================================================

options(ergm.loglik.warn_dyads = FALSE)
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  library(devtools)
  library(ergm)
})

devtools::load_all(".")

if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# ------------------------------------------------------------------------------ 
# Données fixes
# Groupes effectifs = {1,2,3,4} ; tailles = (1,4,2,3) ; N1 = 10
# ------------------------------------------------------------------------------ 
partition <- c(1, 2, 2, 2, 2, 3, 3, 4, 4, 4)
cat("Partition :", paste(partition, collapse = " "), "\n")
cat("Tailles groupes :", paste(as.integer(table(partition)), collapse = ", "), "\n\n")

# ------------------------------------------------------------------------------ 
# 1) Réseau biparti via wrapper (build_bipartite_from_inputs)
# ------------------------------------------------------------------------------ 
bld <- build_bipartite_from_inputs(partition = partition)
nw  <- bld$network

# ------------------------------------------------------------------------------ 
# 2) SUMMARY de référence direct sur le biparti
# ------------------------------------------------------------------------------ 
summary_reference_val     <- as.numeric(summary(nw ~ log_factorial_sizes(), constraints = ~ b1part))
summary_reference_closed  <- sum(lgamma(as.integer(table(partition))))
cat("[summary|référence] network ~ log_factorial_sizes() :", summary_reference_val, "\n")
cat("[summary|fermée   ] Σ_g lgamma(n_g)                :", summary_reference_closed, "\n\n")
stopifnot(isTRUE(all.equal(summary_reference_val, summary_reference_closed, tol = 0)))

# ------------------------------------------------------------------------------ 
# 3) SUMMARY observé via erpm(dry-run) -> formule -> summary(...)
# ------------------------------------------------------------------------------ 
dry <- erpm(partition ~ log_factorial_sizes, eval_call = FALSE, verbose = TRUE)
fml <- dry[[2]]
cons <- if (length(dry) >= 3L && inherits(dry[[3]], "formula")) dry[[3]] else ~ b1part

summary_observed_via_erpm_dry <- as.numeric(summary(fml, constraints = cons))
cat("[summary|erpm(dry)]                                 :", summary_observed_via_erpm_dry, "\n\n")
stopifnot(isTRUE(all.equal(summary_observed_via_erpm_dry, summary_reference_val, tol = 0)))

# ------------------------------------------------------------------------------ 
# 4) FIT via erpm() en  logLik (appel direct et concis)
# ------------------------------------------------------------------------------ 
set.seed(1)
fit_erpm_mle <- erpm(
  partition ~ log_factorial_sizes,
  # estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE
)

cat("\n--- summary(fit_erpm_mle) ---\n")
print(summary(fit_erpm_mle))

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)