# ======================================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_groups.R
#
# Objet   : MWE pour l’effet ERPM `groups(k)` — calcule la stat observée (summary)
#           et effectue un fit ERGM complet via `erpm()`, avec appels directs.
#
# Chaîne  : partition -> biparti (build_bipartite_from_inputs) -> summary(nw ~ ...)
#           partition -> erpm(partition ~ ...) -> fit MLE
#
# Résumé technique :
#   • `groups(k)` ≡ nombre de groupes de taille EXACTEMENT k.
#   • Le wrapper traduit `groups(k)` en {ergm} `b2degrange(from=k, to=k+1)`.
#   • On compare un summary direct sur `nw` et un summary via dry-run `erpm`.
# ======================================================================================

# ----- Préambule ----------------------------------------------------------------------
options(ergm.loglik.warn_dyads = FALSE)  # (4) demandé
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  library(devtools)   # expose erpm(), build_bipartite_from_inputs, etc.
  library(ergm)
})

devtools::load_all(".")

# Patch {ergm} optionnel si présent
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R"); ergm_patch_enable()
}

# ----- Données de test ---------------------------------------------------------------
# 4 groupes {1,2,3,4} ; tailles = (1,3,2,2)
partition <- c(1, 2, 2, 2, 3, 3, 4, 4)
k_group   <- 2  # groups(k) : taille EXACTE recherchée

cat("\nPartition :", paste(partition, collapse = " "), "\n")
cat("Tailles groupes :", paste(as.integer(table(partition)), collapse = ", "), "\n\n")

# ----- (1) Biparti via wrapper -------------------------------------------------------
bld <- build_bipartite_from_inputs(partition = partition)
nw  <- bld$network

# ======================================================================================
# 1) SUMMARY DIRECT sur le biparti (référence)
# ======================================================================================
stat_summary_ref_groups_k <- as.numeric(
  summary(nw ~ b2degrange(from = k_group, to = k_group + 1), constraints = ~ b1part)
)
cat(sprintf("[summary|ref] groups(k=%d) = %g\n", k_group, stat_summary_ref_groups_k))

# ======================================================================================
# 2) SUMMARY via erpm(dry-run) puis summary(formule) — vérification
# ======================================================================================
dry_erpm_groups_k <- erpm(partition ~ groups(k_group), eval_call = FALSE, verbose = FALSE)
stat_summary_via_dry_groups_k <- as.numeric(summary(dry_erpm_groups_k[[2]], constraints = ~ b1part))
cat(sprintf("[summary|erpm(dry)] groups(k=%d) = %g\n\n", k_group, stat_summary_via_dry_groups_k))

stopifnot(isTRUE(all.equal(stat_summary_via_dry_groups_k, stat_summary_ref_groups_k, tol = 0)))

# ======================================================================================
# 3) FIT ERPM COMPLET — (3) demandé : estimate="MLE", eval.loglik=TRUE
# ======================================================================================
set.seed(1)
fit_erpm_groups_k_mle <- erpm(
  partition ~ groups(k_group),
  estimate    = "MLE",
  eval.loglik = TRUE,
  verbose     = TRUE
)

cat("\n--- summary(fit_erpm_groups_k_mle) ---\n")
print(summary(fit_erpm_groups_k_mle))

# Patch off proprement si activé
on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)