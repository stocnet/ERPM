# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cov_ingroup.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_ingroup`
# Exécution: Rscript scripts/test/selftests/selftest_cov_ingroup.R
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule
# --------------------------------------------------------------------------------------
options(ergm.loglik.warn_dyads = FALSE)
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",    quietly = TRUE)) stop("Package 'ergm' requis.")
})

suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

# Patch ERGM optionnel
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# Charger le package et le wrapper ERPM
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Le fichier DESCRIPTION n'existe pas ou devtools n'est pas installé.")
}
if (!exists("erpm", mode = "function")) {
  if (file.exists("R/erpm_wrapper.R")) {
    source("R/erpm_wrapper.R", local = FALSE)
  } else stop("erpm_wrapper.R introuvable.")
}
if (!exists("build_bipartite_from_inputs", mode = "function")) {
  stop("build_bipartite_from_inputs() indisponible. Il doit être exporté par R/erpm_wrapper.R.")
}

# --------------------------------------------------------------------------------------
# Logging local
# --------------------------------------------------------------------------------------
.get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  i <- grep("^--file=", a)
  if (length(i)) return(dirname(normalizePath(sub("^--file=", "", a[i[1]]))))
  fs <- sys.frames()
  ofiles <- vapply(fs, function(f) if (!is.null(f$ofile)) f$ofile else NA_character_, "")
  if (any(!is.na(ofiles))) {
    j <- which.max(nchar(ofiles))
    return(dirname(normalizePath(ofiles[j])))
  }
  normalizePath(getwd())
}

root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_cov_ingroup.log")
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
if (file.exists(log_path)) unlink(log_path, force = TRUE)
con_out <- file(log_path, open = "wt")
con_err <- file(log_path, open = "at")
sink(con_out, split = TRUE)
sink(con_err, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(close(con_err),        silent = TRUE)
  try(sink(),                silent = TRUE)
  try(close(con_out),        silent = TRUE)
  flush.console()
}, add = TRUE)
cat("==> Log:", log_path, "\n")

# ======================================================================================
# Données de test
# ======================================================================================
partitions <- list(
  P1 = c(1,2,2,3,3,3,4),
  P2 = c(1,1, 2,2,2, 3,3,3,3, 4,4, 5,5,5, 6,6,6,6,6,6),
  P3 = c(1,1,2,2,3)
)

.make_nodes_df_for_partition <- function(part) {
  n <- length(part)
  set.seed(123 + n)
  data.frame(
    label   = utils::head(LETTERS, n),
    age     = sample(20:60, n, replace = TRUE),
    score   = round(runif(n, min = 0, max = 10), 1),
    gender  = sample(c("F","M"), n, replace = TRUE),
    dept    = sample(c("A","B","C"), n, replace = TRUE, prob=c(0.5,0.3,0.2)),
    stringsAsFactors = FALSE
  )
}

# ======================================================================================
# Helpers réseau via builder du wrapper
# ======================================================================================

# Construit le biparti via build_bipartite_from_inputs(partition=..., nodes=...)
make_network_from_partition_via_builder <- function(partition_vec, nodes_df) {
  stopifnot(is.atomic(partition_vec), nrow(nodes_df) == length(partition_vec))
  built <- build_bipartite_from_inputs(partition = partition_vec, nodes = nodes_df)
  if (is.list(built) && !is.null(built$network)) return(built$network)
  if (inherits(built, "network")) return(built)
  stop("Le builder n’a pas renvoyé un objet 'network'.")
}

# Construit une formule `nw ~ <rhs>` avec capture de nw
make_formula_for_network_summary <- function(nw, rhs_txt) {
  f <- as.formula(paste0("nw ~ ", rhs_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  f
}

# ======================================================================================
# Fonctions explicitement nommées pour SUMMARY et ERPM
# ======================================================================================

# Summary direct sur le réseau construit via le builder
run_one_network_summary_case_for_cov_ingroup <- function(partition_vec, nodes_df, rhs_txt) {
  nw <- make_network_from_partition_via_builder(partition_vec, nodes_df)
  f  <- make_formula_for_network_summary(nw, rhs_txt)
  as.numeric(suppressMessages(summary(f, constraints = ~ b1part)))
}

# Summary côté ERPM: traduction via erpm(eval_call=FALSE), exécution du summary sur réseau équivalent
run_one_erpm_translated_summary_case_for_cov_ingroup <- function(partition_vec, nodes_df, rhs_txt) {
  partition <- partition_vec
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition, nodes = nodes_df), parent = parent.frame())

  call_ergm <- erpm(f, eval_call = FALSE, verbose = FALSE, nodes = nodes_df)

  # Extraire RHS depuis la formule traduite
  ergm_form <- call_ergm[[2L]]
  rhs_expr  <- ergm_form[[3L]]

  # Reconstruire le réseau via builder et exécuter summary avec la même RHS
  nw2 <- make_network_from_partition_via_builder(partition_vec, nodes_df)
  f2  <- as.formula(bquote(nw2 ~ .(rhs_expr)))
  environment(f2) <- list2env(list(nw2 = nw2), parent = parent.frame())

  # Contrainte: réutiliser celle de l'appel si fournie, sinon ~b1part
  call_args <- as.list(call_ergm)[-1L]
  cons <- call_args$constraints
  if (is.null(cons)) cons <- as.formula(~ b1part)

  as.numeric(suppressMessages(summary(f2, constraints = cons)))
}

# Vérifie l'équivalence exacte des deux summary
check_summary_equivalence_network_vs_erpm <- function(partition_vec, nodes_df, rhs_vec, tol = 0) {
  ok_all <- TRUE
  for (rhs in rhs_vec) {
    s_net  <- run_one_network_summary_case_for_cov_ingroup(partition_vec, nodes_df, rhs)
    s_erpm <- run_one_erpm_translated_summary_case_for_cov_ingroup(partition_vec, nodes_df, rhs)
    cat(sprintf("[SUMMARY-CHECK] n=%-3d RHS=%-50s net=%s  erpm=%s\n",
                length(partition_vec), rhs,
                paste(s_net,  collapse=","), paste(s_erpm, collapse=",")))
    if (any(!is.finite(s_net)) || any(!is.finite(s_erpm)) ||
        length(s_net) != length(s_erpm) ||
        !all(abs(s_net - s_erpm) <= tol)) {
      ok_all <- FALSE
      cat("  -> MISMATCH détecté.\n")
    }
  }
  ok_all
}

# ======================================================================================
# Panel de cas `cov_ingroup`
# ======================================================================================

cases_summary <- c(
  "cov_ingroup('age')",
  "cov_ingroup('age', size = 2:4)",
  "cov_ingroup('score')",
  "cov_ingroup('score', size = 3:6)",
  "cov_ingroup('gender', category='F')",
  "cov_ingroup('gender', category='M', size = 2:5)",
  "cov_ingroup('dept',   category='A')",
  "cov_ingroup('dept',   category='C', size = 2:4)"
)

# ======================================================================================
# Contrôles ERGM pour fitting
# ======================================================================================

ctrl <- control.ergm(
  MCMC.samplesize = 1000,
  MCMLE.maxit     = 20
)

# ======================================================================================
# Phase 1: SUMMARY comparatifs (réseau explicite vs ERPM traduit)
# ======================================================================================

run_phase1_summary_equivalence_checks <- function() {
  cat("=== PHASE 1 : Summary(nw via builder) vs Summary(ERPM-traduit) ===\n")
  total <- 0L; ok <- 0L
  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes_df_for_partition(part)
    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))
    res <- check_summary_equivalence_network_vs_erpm(part, nodes, cases_summary, tol = 0)
    total <- total + length(cases_summary)
    ok    <- ok + as.integer(res) * length(cases_summary)
  }
  cat(sprintf("\n=== Bilan Phase 1 : %d / %d checks OK ===\n", ok, total))
  if (ok < total) stop(sprintf("Summary mismatch sur %d cas.", total - ok))
  invisible(NULL)
}

# ======================================================================================
# Phase 2: Fits via ERPM
# ======================================================================================

run_one_erpm_fit_with_return <- function( partition_vec, 
                                          nodes_df, rhs_txt, 
                                          fit_name,
                                          estimate = NULL,
                                          control = NULL, 
                                          eval.loglik = TRUE
                                        ) {
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT %-20s] SKIP (erpm() indisponible)\n", fit_name))
    return(list(ok = NA, error = FALSE, coef = NA, fit = NULL))
  }
  set.seed(42)
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition_vec, nodes = nodes_df), parent = parent.frame())
  cat(sprintf("[ERPM-FIT %-20s] n=%-3d RHS=%s  | estimate=%s eval.loglik=%s\n",
              fit_name, length(partition_vec), rhs_txt, estimate, as.character(eval.loglik)))
  fit <- try(
    erpm(f, eval.loglik = eval.loglik,
            # estimate = estimate,  
            # control = control, 
            verbose = FALSE, 
            nodes = nodes_df),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    cat("  -> ERREUR fit:", as.character(fit), "\n")
    return(list(ok = FALSE, error = TRUE, coef = NA, fit = NULL))
  }
  ok_class <- inherits(fit, "ergm")
  cf <- try(stats::coef(fit), silent = TRUE)
  ok_coef <- !inherits(cf, "try-error") && all(is.finite(cf))
  cat(sprintf("  -> class(ergm)? %s | coef finies? %s | coef: %s\n",
              if (ok_class) "OK" else "KO",
              if (ok_coef) "OK" else "KO",
              if (ok_coef) paste(format(as.numeric(cf)), collapse=", ") else "NA"))
  list(ok = ok_class && ok_coef, error = FALSE, coef = if (ok_coef) cf else NA, fit = fit)
}

run_phase2_erpm_fits_and_print_summaries <- function() {
  cat("\n=== PHASE 2 : Fits erpm() ( + logLik) ===\n")

  rhs_list <- list(
    R1 = "cov_ingroup('age')",
    R2 = "cov_ingroup('gender', category='F', size=2:4)",
    R3 = "cov_ingroup('dept',   category='A')",
    R4 = "cov_ingroup('age') + cov_ingroup('gender', category='F', size=3:6)"
  )

  parts <- list(
    list(name="P1", part=partitions$P1),
    list(name="P2", part=partitions$P2)
  )

  fit_results <- list()
  for (px in parts) {
    nodes <- .make_nodes_df_for_partition(px$part)
    for (nm in names(rhs_list)) {
      key <- paste0(px$name, "_", nm)
      fit_results[[key]] <- run_one_erpm_fit_with_return(
        partition_vec = px$part,
        nodes_df      = nodes,
        rhs_txt       = rhs_list[[nm]],
        fit_name      = key,
        # estimate      = "CD",
        # control       = ctrl,
        eval.loglik   = TRUE
      )
    }
  }

  ok <- vapply(fit_results, function(x) isTRUE(x$ok), logical(1))
  n_ok <- sum(ok, na.rm = TRUE); n_tot <- sum(!is.na(ok))
  cat(sprintf("\n=== Bilan fits erpm() : %d / %d OK ===\n", n_ok, n_tot))

  cat("\n=== Résumés détaillés des fits ERPM réussis ===\n")
  for (nm in names(fit_results)) {
    fit_obj <- fit_results[[nm]]
    if (isTRUE(fit_obj$ok) && inherits(fit_obj$coef, "numeric") && !is.null(fit_obj$fit)) {
      cat(sprintf("\n--- Résumé fit %s ---\n", nm))
      print(summary(fit_obj$fit))
    }
  }

  if (n_ok < n_tot) stop(sprintf("Echec fits: %d KO", n_tot - n_ok))
  invisible(fit_results)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: cov_ingroup ===\n")
run_phase1_summary_equivalence_checks()
fit_results <- run_phase2_erpm_fits_and_print_summaries()

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)
cat("\nTous les tests cov_ingroup ont passé.\n")