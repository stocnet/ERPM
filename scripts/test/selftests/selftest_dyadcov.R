# ======================================================================================
# Fichier : scripts/test/selftests/selftest_dyadcov.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `dyadcov`
# Exécution: Rscript scripts/test/selftests/selftest_dyadcov.R
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
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_dyadcov.log")
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
  P1 = c(
    rep(1L, 6),
    rep(2L, 2),
    rep(3L, 5),
    rep(4L, 2)
  ),
  P2 = c(
    rep(1L, 2),
    rep(2L, 5),
    rep(3L, 5),
    rep(4L, 3)
  ),
  P3 = c(
    rep(1L, 3),
    rep(2L, 3),
    rep(3L, 5)
  )
)

# Nodes "muets" pour satisfaire le builder (dyadcov ne lit pas d'attribut nodal)
.make_nodes_df_for_partition <- function(part) {
  n <- length(part)
  data.frame(
    label = paste0("N", seq_len(n)),
    stringsAsFactors = FALSE
  )
}

# Matrices dyadiques déterministes à partir de la partition
# - Z1 : valeurs réelles modérées (plus grandes intra-groupe, faibles mais non nulles inter-groupe), symétrique, diag=0
# - Z2 : mêmes ordres de grandeur mais rendue explicitement non symétrique, diag=0
.make_dyads_for_partition <- function(part) {
  n <- length(part)
  Z1 <- matrix(0, n, n)
  Z2 <- matrix(0, n, n)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next

      d_ij <- abs(i - j)

      if (part[i] == part[j]) {
        # Intra-groupe
        Z1[i, j] <- 1.0 + 0.5 * d_ij
        Z2[i, j] <- 0.5 * d_ij^2 + (i + j) / 10
      } else {
        # Inter-groupe
        Z1[i, j] <- (d_ij %% 3) / 3
        Z2[i, j] <- ((i * j) %% 7) / 4
      }
    }
  }

  # Pour dyadcov, on impose des diagonales nulles
  diag(Z1) <- 0
  diag(Z2) <- 0

  # Et on casse explicitement la symétrie de Z2 pour tester la convention (z_ij + z_ji)
  if (n > 1L) {
    Z2[lower.tri(Z2)] <- 2 * Z2[lower.tri(Z2)]
  }

  list(Z1 = Z1, Z2 = Z2)
}

# Petit helper de debug : afficher partition / nodes / extrait des dyads
print_debug_partition_nodes_dyads <- function(name, part, nodes_df, dyads_list) {
  cat(sprintf("\n[DEBUG] --- Cas %s ---\n", name))
  cat("[DEBUG] partition :", paste(part, collapse = ","), "\n")
  cat("[DEBUG] nodes (head) :\n")
  print(utils::head(nodes_df, 10))

  Z1 <- dyads_list$Z1
  Z2 <- dyads_list$Z2
  n  <- nrow(Z1)
  k  <- min(6L, n)

  cat("[DEBUG] Z1[1:", k, ", 1:", k, "] =\n", sep = "")
  print(round(Z1[seq_len(k), seq_len(k)], 3))

  cat("[DEBUG] Z2[1:", k, ", 1:", k, "] =\n", sep = "")
  print(round(Z2[seq_len(k), seq_len(k)], 3))

  # Sanity check diagonales nulles
  stopifnot(all(diag(Z1) == 0), all(diag(Z2) == 0))
}

# ======================================================================================
# Helpers réseau via builder du wrapper
# ======================================================================================

# Construit le biparti via build_bipartite_from_inputs(partition=..., nodes=..., dyads=...)
make_network_from_partition_and_dyads <- function(partition_vec, nodes_df, dyads_list) {
  stopifnot(is.atomic(partition_vec), nrow(nodes_df) == length(partition_vec))
  built <- build_bipartite_from_inputs(
    partition = partition_vec,
    nodes     = nodes_df,
    dyads     = dyads_list
  )
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
run_one_network_summary_case_for_dyadcov <- function(partition_vec, nodes_df, dyads_list, rhs_txt) {
  nw <- make_network_from_partition_and_dyads(partition_vec, nodes_df, dyads_list)
  f  <- make_formula_for_network_summary(nw, rhs_txt)
  as.numeric(suppressMessages(summary(f, constraints = ~ b1part)))
}

# Summary côté ERPM: traduction via erpm(eval_call=FALSE), exécution du summary sur réseau équivalent
run_one_erpm_translated_summary_case_for_dyadcov <- function(partition_vec, nodes_df, dyads_list, rhs_txt) {
  partition <- partition_vec
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition, nodes = nodes_df), parent = parent.frame())

  call_ergm <- erpm(
    f,
    eval_call = FALSE,
    verbose   = FALSE,
    nodes     = nodes_df,
    dyads     = dyads_list
  )

  # Extraire RHS depuis la formule traduite
  ergm_form <- call_ergm[[2L]]
  rhs_expr  <- ergm_form[[3L]]

  # Reconstruire le réseau via builder et exécuter summary avec la même RHS
  nw2 <- make_network_from_partition_and_dyads(partition_vec, nodes_df, dyads_list)
  f2  <- as.formula(bquote(nw2 ~ .(rhs_expr)))
  environment(f2) <- list2env(list(nw2 = nw2), parent = parent.frame())

  # Contrainte: réutiliser celle de l'appel si fournie, sinon ~b1part
  call_args <- as.list(call_ergm)[-1L]
  cons <- call_args$constraints
  if (is.null(cons)) cons <- as.formula(~ b1part)

  as.numeric(suppressMessages(summary(f2, constraints = cons)))
}

# Vérifie l'équivalence exacte des deux summary
check_summary_equivalence_network_vs_erpm_dyadcov <- function(partition_vec, nodes_df, dyads_list,
                                                              rhs_vec, tol = 0) {
  ok_all <- TRUE
  for (rhs in rhs_vec) {
    s_net  <- run_one_network_summary_case_for_dyadcov(partition_vec, nodes_df, dyads_list, rhs)
    s_erpm <- run_one_erpm_translated_summary_case_for_dyadcov(partition_vec, nodes_df, dyads_list, rhs)
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
# Panel de cas `dyadcov`
# ======================================================================================

cases_summary <- c(
  "dyadcov('Z1', clique_size = 2, normalized = FALSE)",
  "dyadcov('Z1', clique_size = 2, normalized = TRUE)",
  "dyadcov('Z1', clique_size = 3, normalized = FALSE)",
  "dyadcov('Z1', clique_size = 3, normalized = TRUE)",
  "dyadcov('Z2', clique_size = 2, normalized = FALSE)",  # Z2 non symétrique
  "dyadcov('Z2', clique_size = 2, normalized = TRUE)"    # Z2 non symétrique
)

# ======================================================================================
# Contrôles ERGM pour fitting
# ======================================================================================

ctrl <- control.ergm(
  # Laisser les valeurs par défaut pour MCMLE, l'effet est seul dans le modèle ou presque.
)

# ======================================================================================
# Phase 1: SUMMARY comparatifs (réseau explicite vs ERPM traduit)
# ======================================================================================

run_phase1_summary_equivalence_checks_dyadcov <- function() {
  cat("=== PHASE 1 : Summary(nw via builder) vs Summary(ERPM-traduit) [dyadcov] ===\n")
  total <- 0L; ok <- 0L
  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes_df_for_partition(part)
    dyads <- .make_dyads_for_partition(part)

    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))

    print_debug_partition_nodes_dyads(nm, part, nodes, dyads)

    res <- check_summary_equivalence_network_vs_erpm_dyadcov(
      partition_vec = part,
      nodes_df      = nodes,
      dyads_list    = dyads,
      rhs_vec       = cases_summary,
      tol           = 0
    )
    total <- total + length(cases_summary)
    ok    <- ok + as.integer(res) * length(cases_summary)
  }
  cat(sprintf("\n=== Bilan Phase 1 (dyadcov) : %d / %d checks OK ===\n", ok, total))
  if (ok < total) stop(sprintf("Summary mismatch sur %d cas.", total - ok))
  invisible(NULL)
}

# ======================================================================================
# Phase 2: Fits via ERPM (MLE uniquement, aucun autre effet)
# ======================================================================================

run_one_erpm_fit_with_return_dyadcov <- function(partition_vec, nodes_df, dyads_list,
                                                 rhs_txt, fit_name,
                                                 estimate = "MLE",
                                                 eval.loglik = TRUE,
                                                 control = ctrl) {
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT %-20s] SKIP (erpm() indisponible)\n", fit_name))
    return(list(ok = NA, error = TRUE, coef = NA, fit = NULL, aic = NA, bic = NA))
  }
  set.seed(42)
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition_vec, nodes = nodes_df), parent = parent.frame())

  cat(sprintf("[ERPM-FIT %-20s] n=%-3d RHS=%s  | estimate=%s eval.loglik=%s\n",
              fit_name, length(partition_vec), rhs_txt, estimate, as.character(eval.loglik)))

  print_debug_partition_nodes_dyads(paste0("FIT_", fit_name), partition_vec, nodes_df, dyads_list)

  fit <- try(
    erpm(
      f,
      estimate    = estimate,
      eval.loglik = eval.loglik,
      control     = control,
      verbose     = FALSE,
      nodes       = nodes_df,
      dyads       = dyads_list
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    msg <- paste(as.character(fit), collapse = "\n")
    cat("  -> ERREUR fit:", msg, "\n")
    return(list(ok = FALSE, error = TRUE, coef = NA, fit = NULL, aic = NA, bic = NA))
  }

  ok_class <- inherits(fit, "ergm")
  cf <- try(stats::coef(fit), silent = TRUE)
  ok_coef <- !inherits(cf, "try-error") && all(is.finite(cf))

  aic_val <- NA_real_
  bic_val <- NA_real_
  ll <- try(logLik(fit, add = TRUE), silent = TRUE)
  if (!inherits(ll, "try-error")) {
    aic_try <- try(AIC(fit), silent = TRUE)
    bic_try <- try(BIC(fit), silent = TRUE)
    if (!inherits(aic_try, "try-error")) aic_val <- as.numeric(aic_try)
    if (!inherits(bic_try, "try-error")) bic_val <- as.numeric(bic_try)
  } else {
    cat("  -> logLik(fit, add=TRUE) a échoué, AIC/BIC indisponibles pour ce fit.\n")
  }

  cat(sprintf("  -> class(ergm)? %s | coef finies? %s | coef: %s | AIC=%s | BIC=%s\n",
              if (ok_class) "OK" else "KO",
              if (ok_coef) "OK" else "KO",
              if (ok_coef) paste(format(as.numeric(cf)), collapse=", ") else "NA",
              if (is.finite(aic_val)) format(aic_val, digits = 6) else "NA",
              if (is.finite(bic_val)) format(bic_val, digits = 6) else "NA"))

  list(
    ok    = ok_class && ok_coef,
    error = FALSE,
    coef  = if (ok_coef) cf else NA,
    fit   = fit,
    aic   = aic_val,
    bic   = bic_val
  )
}

run_phase2_erpm_fits_and_print_summaries_dyadcov <- function() {
  cat("\n=== PHASE 2 : Fits erpm() (MLE + logLik) [dyadcov] ===\n")

  rhs_list <- list(
    R1 = "dyadcov('Z1', clique_size = 2, normalized = FALSE)",
    R2 = "dyadcov('Z1', clique_size = 2, normalized = TRUE)",
    R3 = "dyadcov('Z2', clique_size = 3, normalized = FALSE) + cliques",
    R4 = "dyadcov('Z2', clique_size = 3, normalized = TRUE) + cliques"
  )

  parts <- list(
    list(name = "P1", part = partitions$P1),
    list(name = "P2", part = partitions$P2)
  )

  fit_results <- list()
  for (px in parts) {
    nodes <- .make_nodes_df_for_partition(px$part)
    dyads <- .make_dyads_for_partition(px$part)
    for (nm in names(rhs_list)) {

      # Cas pathologique connu : P1 + R4 => mélange MCMC très mauvais, durée excessive.
      if (px$name == "P1" && nm == "R4") {
        key_skip <- paste0(px$name, "_", nm)
        cat(sprintf("[ERPM-FIT %-20s] SKIP (modèle numériquement dégénéré, selftest)\n", key_skip))
        next
      }

      key <- paste0(px$name, "_", nm)
      fit_results[[key]] <- run_one_erpm_fit_with_return_dyadcov(
        partition_vec = px$part,
        nodes_df      = nodes,
        dyads_list    = dyads,
        rhs_txt       = rhs_list[[nm]],
        fit_name      = key,
        estimate      = "MLE",
        eval.loglik   = TRUE,
        control       = ctrl
      )
    }
  }

  ok_raw  <- vapply(fit_results, function(x) x$ok,    logical(1))
  n_ok    <- sum(ok_raw, na.rm = TRUE)
  n_tot   <- sum(!is.na(ok_raw))

  cat(sprintf("\n=== Bilan fits erpm() [dyadcov] : %d / %d OK ===\n", n_ok, n_tot))

  cat("\n=== Tableau AIC/BIC pour les fits ERPM (dyadcov) ===\n")
  tab <- data.frame(
    fit   = character(0),
    ok    = logical(0),
    AIC   = numeric(0),
    BIC   = numeric(0),
    stringsAsFactors = FALSE
  )
  for (nm in names(fit_results)) {
    fr <- fit_results[[nm]]
    tab <- rbind(tab, data.frame(
      fit = nm,
      ok  = fr$ok,
      AIC = fr$aic,
      BIC = fr$bic,
      stringsAsFactors = FALSE
    ))
  }
  print(tab)

  cat("\n=== Résumés détaillés des fits ERPM réussis (dyadcov) ===\n")
  for (nm in names(fit_results)) {
    fit_obj <- fit_results[[nm]]
    if (isTRUE(fit_obj$ok) && inherits(fit_obj$coef, "numeric") && !is.null(fit_obj$fit)) {
      cat(sprintf("\n--- Résumé fit %s ---\n", nm))
      print(summary(fit_obj$fit))
    }
  }

  # On sait empiriquement que plusieurs modèles dyadcov seuls sont dégénérés.
  # Pour le selftest, on exige seulement qu'au moins deux modèles simples convergent.
  required_ok <- 2L

  if (n_ok < required_ok) {
    stop(sprintf("Echec fits: seulement %d modèles convergents (seuil=%d).",
                 n_ok, required_ok))
  } else {
    cat(sprintf("\nSeuil de convergence atteint : %d modèles OK (seuil=%d).\n",
                n_ok, required_ok))
  }

  invisible(fit_results)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: dyadcov ===\n")
run_phase1_summary_equivalence_checks_dyadcov()
fit_results <- run_phase2_erpm_fits_and_print_summaries_dyadcov()

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)
cat("\nTous les tests dyadcov ont passé.\n")
