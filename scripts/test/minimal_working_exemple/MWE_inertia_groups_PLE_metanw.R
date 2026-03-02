# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_inertia_groups_PLE.R
# Objet   : MWE pour l’effet inertiel ERPM `inertia_groups` (PLE empilé)
# Chaîne  : partitions (5) -> erpm_long(mode="empile") ->
#           summary(nw ~ ...) / fit(erpm_long(...))
#
# Contraintes demandées :
#   - past_influence = 2
#   - 5 partitions, 8 acteurs
#   - RHS : dyadcov + cov_match + cliques + inertia_groups(size=2)
#   - plusieurs attributs noeudaux et dyadiques
#   - 1 summary (dry-run) + 1 fit
#   - capture + print des warnings (buffer global + withCallingHandlers)
#   - options erpm_long : mode="empile", verbose=TRUE, debug="deep"
#   - formule erpm_long au format : partitions ~ RHS
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) stop("Package 'devtools' requis.")
  if (!requireNamespace("network",  quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",     quietly = TRUE)) stop("Package 'ergm' requis.")
})

# Load ERPM in dev mode (adjust if you prefer library(ERPM))
options(keep.source = TRUE)
options(keep.source.pkgs = TRUE)
Sys.setenv(R_KEEP_PKG_SOURCE = "yes")
# devtools::load_all(".")

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
# Patch ergm si présent (avec nettoyage à la sortie)
# ----------------------------------------------------------------------
if (file.exists("scripts/ergm_patch.R")) {
  cat("\n=== [STEP] ergm_patch_enable() ===\n")
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}
# We do NOT want to run engine or ergm. So we force eval.call=FALSE always.
.call_erpm_long <- function(formula, mode = "empile", nodes = NULL, dyads = NULL,
                            verbose = FALSE, debug = FALSE, seed = NULL, group_labels = NULL) {
  erpm_long(
    formula   = formula,
    mode      = mode,
    nodes     = nodes,
    dyads     = dyads,
    verbose   = verbose,
    debug     = debug,
    seed      = seed,
    group_labels = group_labels,
    eval.call = FALSE
  )
}
# Nodes (explicit, per time)
nodes <- list(
  data.frame(
    label  = c("A", "B", "C"),
    gender = c(1, 1, 2),
    age    = c(20, 22, 25)
  ),
  data.frame(
    label  = c("A", "AZ", "JI"),
    gender = c(1, 2, 2),
    age    = c(20, 42, 30)
  ),
  data.frame(
    label  = c("A", "B", "Z"),
    gender = c(1, 1, 1),
    age    = c(20, 23, 2)
  )
)

# Dyads (explicit, per time)
dyads <- list(
  list(
    fm = matrix(c(
      0, 1, 1,
      1, 0, 1,
      1, 1, 0
    ), nrow = 3, byrow = TRUE),
    Z1 = matrix(c(
      0, 2, 3,
      2, 0, 4,
      3, 4, 0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    fm = matrix(c(
      0, 1, 0,
      1, 0, 1,
      0, 1, 0
    ), nrow = 3, byrow = TRUE),
    Z1 = matrix(c(
      0, 5, 0,
      5, 0, 1,
      0, 1, 0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    fm = matrix(c(
      0, 0, 1,
      0, 0, 1,
      1, 1, 0
    ), nrow = 3, byrow = TRUE),
    Z1 = matrix(c(
      0, 1, 2,
      1, 0, 4,
      2, 4, 0
    ), nrow = 3, byrow = TRUE)
  )
)

# Partitions passed to erpm_long (explicit)
partitions <- list(
  c(1, 1, 2),
  c(1, 2, 2),
  c(1, 1, 2)
)

# Nodes (explicit, per time)
nodes2 <- list(
  data.frame(
    label  = c("A", "B", "C"),
    gender = c(1, 1, 2),
    age    = c(20, 22, 25)
  ),
  data.frame(
    label  = c("A", "B"),
    gender = c(1, 1),
    age    = c(20, 22)
  ),
  data.frame(
    label  = c("A"),
    gender = c(1),
    age    = c(20)
  )
)

# Dyads (explicit, per time)
dyads2 <- list(
  list(
    fm = matrix(c(
      0, 1, 1,
      1, 0, 1,
      1, 1, 0
    ), nrow = 3, byrow = TRUE),
    Z1 = matrix(c(
      0, 2, 3,
      2, 0, 2,
      3, 2, 0
    ), nrow = 3, byrow = TRUE)
  ),
  list(
    fm = matrix(c(
      0, 1,
      1, 0
    ), nrow = 2, byrow = TRUE),
    Z1 = matrix(c(
      0, 5,
      5, 0
    ), nrow = 2, byrow = TRUE)
  ),
  list(
    fm = matrix(c(
      0
    ), nrow = 1, byrow = TRUE),
    Z1 = matrix(c(
      0
    ), nrow = 1, byrow = TRUE)
  )
)

# Partitions passed to erpm_long (explicit)
partitions2 <- list(
  c(1, 2, 2),
  c(1, 2),
  c(1)
)
# ----------------------------------------------------------------------
# RHS 
#   - dyadcov : utilise Z1 (k=2 par défaut si votre API est dyadcov("Z1"))
#   - cov_match : utilise l'attribut 'sector'
#   - cliques : simple structure de groupes
#   - inertia_groups : past_influence=2 + size=2 + debug="deep"
# ----------------------------------------------------------------------
cat("\n=== [STEP] RHS ===\n")

rhs1 <- quote(
  cliques(clique_size = 2, normalized = FALSE)
)

rhs2 <- quote(inertia_groups(past_influence = 2))

rhs3 <- quote(
  dyadcov("Z1") +
  cov_match("sector", clique_size = 2, normalized = "by_group") +
  cliques(clique_size = 2, normalized = FALSE) +
  inertia_groups(past_influence = 2, size = 2)
)

cat("RHS1 = "); print(rhs1)
cat("RHS2 = "); print(rhs2)
cat("RHS3 = "); print(rhs3)

# f_long <- stats::as.formula(paste0("partitions ~ ", rhs))
# environment(f_long) <- list2env(list(partitions = partitions), parent = parent.frame())

# ----------------------------------------------------------------------
# 1) SUMMARY (dry-run)
# ----------------------------------------------------------------------
cat("\n================================================================================\n")
cat("=== [SUMMARY] erpm_long(mode='empile', eval.call=FALSE) + summary(network) ===\n")
cat("================================================================================\n")

.call_erpm_long(partitions ~ cliques(clique_size = 2, normalized = FALSE),
                mode="empile",
                nodes=nodes,
                dyads=dyads)
.call_erpm_long(partitions ~ inertia_groups(past_influence = 1),
                mode="empile",
                nodes=nodes,
                dyads=dyads)

.call_erpm_long(partitions2 ~ cliques(clique_size = 2, normalized = FALSE),
                mode="empile",
                nodes=nodes2,
                dyads=dyads2)
.call_erpm_long(partitions2 ~ inertia_groups(past_influence = 1),
                mode="empile",
                nodes=nodes2,
                dyads=dyads2)

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