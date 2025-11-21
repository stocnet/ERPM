# ==============================================================================
# Fichier : scripts/test/selftests/selftest_cov_ingroup.R
# Objet   : Self-test rapide pour l’effet ERPM/ERGM `cov_ingroup`
# Exécution: Rscript scripts/test/selftests/selftest_cov_ingroup.R
# ==============================================================================

options(ergm.loglik.warn_dyads = FALSE)
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  library(network)
  library(ergm)
})

# Patch ERGM si disponible
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# Charger package
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Package non disponible.")
}

# Wrapper ERPM
if (!exists("erpm", mode = "function")) {
  if (file.exists("R/erpm_wrapper.R")) {
    source("R/erpm_wrapper.R", local = FALSE)
  } else {
    stop("R/erpm_wrapper.R introuvable.")
  }
}

if (!exists("build_bipartite_from_inputs", mode = "function")) {
  stop("build_bipartite_from_inputs() manquant.")
}

cat("=== TEST ERPM: cov_ingroup (version rapide) ===\n")

# ==============================================================================
# Données
# ==============================================================================

partitions <- list(
  P1 = c(1,2,3,3,3,4),
  P2 = c(1, 2,2, 3,3,3, 4, 5,5,5, 6,6,6,6,6),
  P3 = c(1,1,2,2,3)
)

.make_nodes <- function(part) {
  n <- length(part)
  set.seed(100 + n)
  data.frame(
    label  = paste0("N", seq_len(n)),
    age    = sample(20:60, n, TRUE),
    score  = round(runif(n, 0, 10), 1),
    gender = sample(c("F","M"), n, TRUE),
    dept   = sample(c("A","B","C"), n, TRUE, prob = c(0.5, 0.3, 0.2)),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Helpers
# ==============================================================================

make_nw <- function(part, nodes) {
  built <- build_bipartite_from_inputs(partition = part, nodes = nodes)
  if (is.list(built) && !is.null(built$network) && inherits(built$network, "network")) {
    return(built$network)
  }
  if (inherits(built, "network")) return(built)
  stop("build_bipartite_from_inputs() ne renvoie pas un 'network'.")
}

make_f <- function(nw, rhs) {
  f <- as.formula(paste0("nw ~ ", rhs))
  environment(f) <- list2env(list(nw = nw))
  f
}

summary_network <- function(part, nodes, rhs) {
  nw <- make_nw(part, nodes)
  f  <- make_f(nw, rhs)
  as.numeric(suppressMessages(summary(f, constraints = ~ b1part)))
}

summary_erpm <- function(part, nodes, rhs) {
  partition <- part
  f <- as.formula(paste0("partition ~ ", rhs))
  environment(f) <- list2env(list(partition = partition, nodes = nodes))

  call_ergm <- erpm(f, eval_call = FALSE, verbose = FALSE, nodes = nodes)
  ergm_form <- call_ergm[[2L]]
  rhs_e     <- ergm_form[[3L]]

  nw2 <- make_nw(part, nodes)
  f2  <- as.formula(bquote(nw2 ~ .(rhs_e)))
  environment(f2) <- list2env(list(nw2 = nw2))

  as.numeric(suppressMessages(summary(f2, constraints = ~ b1part)))
}

compare_case <- function(part, nodes, rhs) {
  s1 <- summary_network(part, nodes, rhs)
  s2 <- summary_erpm(part, nodes, rhs)
  cat(sprintf("[SUMMARY] n=%-3d RHS=%-35s net=%s | erpm=%s\n",
              length(part), rhs,
              paste(s1, collapse=","), paste(s2, collapse=",")))
  stopifnot(length(s1) == length(s2), all(s1 == s2))
}

# ==============================================================================
# Cas à tester
# ==============================================================================

cases <- c(
  "cov_ingroup('age')",
  "cov_ingroup('age', size=2:4)",
  "cov_ingroup('score')",
  "cov_ingroup('gender', category='F')",
  "cov_ingroup('dept', category='A')"
)

# ==============================================================================
# Phase unique : tests SUMMARY
# ==============================================================================

for (nm in names(partitions)) {
  part  <- partitions[[nm]]
  nodes <- .make_nodes(part)
  cat(sprintf("\n--- Partition %s --- n=%d ---\n", nm, length(part)))
  for (rhs in cases) compare_case(part, nodes, rhs)
}

cat("\nTous les tests cov_ingroup ont passé.\n")

if (exists("ergm_patch_disable")) ergm_patch_disable()