# ==============================================================================
# Minimal Working Exemple : inertia_groups — summary + fit
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  library(devtools)
  library(ergm)
})

# ----------------------------------------------------------------------
# Helpers (demandés tels quels)
# ----------------------------------------------------------------------
.warn_flush <- function() {
  # Clear the global warnings buffer
  invisible(warnings())
  invisible(NULL)
}

.capture_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = w)
}

# ----------------------------------------------------------------------
# Charge le package local ERPM
# ----------------------------------------------------------------------
cat("\n=== [STEP] load_all('.') ===\n")
devtools::load_all(".")

# ----------------------------------------------------------------------
# Patch ergm si présent (avec nettoyage à la sortie)
# ----------------------------------------------------------------------
if (file.exists("scripts/ergm_patch.R")) {
  cat("\n=== [STEP] ergm_patch_enable() ===\n")
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# ----------------------------------------------------------------------
# Données fixes (8 acteurs, 5 partitions)
#   - Conçu pour que certaines signatures (taille>=2) persistent
#     sur des fenêtres (t-1,t-2) avec past_influence=2.
# ----------------------------------------------------------------------
cat("\n=== [STEP] data: partitions / nodes / dyads ===\n")

labels <- paste0("A", 1:8)

# 5 partitions strictes (8 acteurs -> 4 groupes)
p1 <- c(1, 1, 2, 2, 3, 3, 4, 4)  # signatures taille>=2 : {1,2}, {3,4}, {5,6}, {7,8}
p2 <- c(1, 2, 2, 2, 3, 3, 4, 4)  # {3,4}, {5,6}, {7,8} restent stables
p3 <- c(1, 2, 2, 2, 3, 3, 4, 4)  # identique p2 (renforce l'intersection)
p4 <- c(1, 2, 2, 2, 3, 4, 4, 4)  # changement : {5} singleton, {6,7,8} taille 3
p5 <- c(1, 2, 2, 2, 3, 4, 4, 4)  # identique p4

partitions <- list(p1, p2, p3, p4, p5)

cat("Partitions (T=5, N=8) :\n")
for (t in seq_along(partitions)) {
  cat(sprintf("  - p%d: %s\n", t, paste(partitions[[t]], collapse = " ")))
}
cat("\nTailles des groupes par temps :\n")
for (t in seq_along(partitions)) {
  tt <- table(partitions[[t]])
  cat(sprintf("  - p%d: %s\n", t, paste(as.integer(tt), collapse = " ")))
}

# Attributs noeudaux (plusieurs)
nodes <- data.frame(
  label  = labels,
  sector = c("S1","S1","S2","S2","S1","S3","S3","S2"),
  region = c("R1","R2","R1","R2","R1","R2","R1","R2"),
  score  = c(10, 12, 9, 15, 7, 11, 14, 8),
  stringsAsFactors = FALSE
)

cat("\nNodes (aperçu) :\n")
print(nodes)

# Dyadiques (2 matrices 8x8 déclarées "en brut")
# Convention : diag=0 ; Z1 symétrique ; Z2 non-symétrique (orientation utile)
Z1 <- matrix(c(
  0,1,2,1,3,2,4,3,
  1,0,1,2,2,3,3,4,
  2,1,0,1,3,2,4,3,
  1,2,1,0,2,3,3,4,
  3,2,3,2,0,1,2,1,
  2,3,2,3,1,0,1,2,
  4,3,4,3,2,1,0,1,
  3,4,3,4,1,2,1,0
), nrow = 8, byrow = TRUE)

Z2 <- matrix(c(
  0,2,0,1,0,3,0,1,
  0,0,2,0,1,0,3,0,
  1,0,0,2,0,1,0,3,
  3,1,0,0,2,0,1,0,
  0,3,1,0,0,2,0,1,
  1,0,3,1,0,0,2,0,
  0,1,0,3,1,0,0,2,
  2,0,1,0,3,1,0,0
), nrow = 8, byrow = TRUE)

dyads <- list(Z1 = Z1, Z2 = Z2)

cat("\nDyads :\n")
cat(sprintf("  - Z1: %dx%d (symétrique)\n", nrow(Z1), ncol(Z1)))
cat(sprintf("  - Z2: %dx%d (non-symétrique)\n", nrow(Z2), ncol(Z2)))

# ----------------------------------------------------------------------
# RHS (effets demandés)
#   - dyadcov : utilise Z1 (k=2 par défaut si votre API est dyadcov("Z1"))
#   - cov_match : utilise l'attribut 'sector'
#   - cliques : simple structure de groupes
#   - inertia_groups : past_influence=2 + size=2 + debug="deep"
# ----------------------------------------------------------------------
cat("\n=== [STEP] RHS ===\n")

rhs <- quote(
  dyadcov("Z1") +
  cov_match("sector", clique_size = 2, normalized = "by_group") +
  cliques(clique_size = 2, normalized = FALSE) +
  inertia_groups(past_influence = 2, size = 2, debug = "deep")
)

cat("RHS = "); print(rhs)

# ----------------------------------------------------------------------
# 1) SUMMARY (dry-run)
# ----------------------------------------------------------------------
cat("\n================================================================================\n")
cat("=== [SUMMARY] erpm_long(mode='empile', eval.call=FALSE) + summary(network) ===\n")
cat("================================================================================\n")

.warn_flush()
res_sum <- .capture_warnings(
  erpm_long(
    partitions ~ eval(rhs),
    nodes     = nodes,
    dyads     = dyads,
    mode      = "empile",
    verbose   = TRUE,
    debug     = "deep",
    eval.call = FALSE
  )
)

if (length(res_sum$warnings)) {
  cat("\n--- WARNINGS (summary) ---\n")
  cat(paste(res_sum$warnings, collapse = "\n"), "\n")
} else {
  cat("\n--- WARNINGS (summary) ---\n<none>\n")
}

out_sum <- res_sum$value
cat("\n--- PRINT out_sum ---\n")
print(out_sum)

# Hypothèse de wrapper : en dry-run, out[[2]] est le network construit
# (comme dans vos MWEs précédents).
nw_sum <- out_sum[[2]]

cat("\n--- summary(nw_sum ~ RHS, constraints=~b1part + blockdiag(timeblock)) ---\n")
# En PLE, blockdiag(timeblock) est attendu. On le met explicitement.
sum_vec <- summary(
  nw_sum ~ eval(rhs),
  constraints = ~ b1part + blockdiag("timeblock")
)
print(sum_vec)

# ----------------------------------------------------------------------
# 2) FIT
# ----------------------------------------------------------------------
cat("\n================================================================================\n")
cat("=== [FIT] erpm_long(mode='empile', eval.call=TRUE) ===\n")
cat("================================================================================\n")

# Contrôle ergm (optionnel)
# ctrl_fit <- control.ergm(
#   CD.maxit        = 0,
#   MCMLE.maxit     = 10,
#   MCMC.burnin     = 5000,
#   MCMC.interval   = 1000,
#   MCMC.samplesize = 1e4,
#   force.main      = TRUE,
#   parallel        = 0
# )

set.seed(1)
.warn_flush()
res_fit <- .capture_warnings(
  erpm_long(
    partitions ~ eval(rhs),
    nodes       = nodes,
    dyads       = dyads,
    mode        = "empile",
    verbose     = TRUE,
    debug       = "deep",
    eval.call   = TRUE,
    eval.loglik = TRUE
    # control   = ctrl_fit
  )
)

if (length(res_fit$warnings)) {
  cat("\n--- WARNINGS (fit) ---\n")
  cat(paste(res_fit$warnings, collapse = "\n"), "\n")
} else {
  cat("\n--- WARNINGS (fit) ---\n<none>\n")
}

fit <- res_fit$value
cat("\n--- PRINT fit ---\n")
print(fit)

cat("\n--- summary(fit) ---\n")
print(summary(fit))

cat("\n=== [DONE] MWE inertia_groups (PLE empilé) ===\n")
on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)