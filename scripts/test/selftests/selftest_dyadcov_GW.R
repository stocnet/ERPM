# ======================================================================================
# Fichier : scripts/test/selftests/selftest_dyadcov_GW.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `dyadcov_GW`
# Exécution: Rscript scripts/test/selftests/selftest_dyadcov_GW.R
#
# But du fichier
#   - PHASE 1 (SUMMARY) : valider summary() (réseau explicite vs ERPM traduit).
#   - PHASE 2 (ERPM FIT): valider qu'un fit erpm() passe et renvoie des coefs finis.
#   - PHASE 3 (MCMC)    : diagnostic "multi-toggle" (intérêt principal ici) :
#                         déclencher des étapes MCMC pouvant contenir ntoggles>1,
#                         et observer les traces C si DEBUG_DYADCOV_GW=1 côté C.
#
# Important (multi-toggle / D_CHANGESTAT_FN)
#   - dyadcov_GW est maintenant un terme multi-toggle (D_ changestat).
#   - Donc InitErgmTerm.dyadcov_GW doit renvoyer d_func=TRUE (sinon crash).
#   - La PHASE 3 ne “prouve” pas mathématiquement le multi-toggle ; elle sert à
#     déclencher le chemin MCMC où ergm peut proposer plusieurs toggles d'un coup
#     (selon MCMC.prop / contraintes) et à vérifier visuellement la présence
#     des logs de type "MULTI-TOGGLE ntoggles=...".
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
  if (!requireNamespace("rprojroot", quietly = TRUE)) stop("Package 'rprojroot' requis.")
})

suppressMessages(suppressPackageStartupMessages({
  library(network,  quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,     quietly = TRUE, warn.conflicts = FALSE)
  library(rprojroot,quietly = TRUE, warn.conflicts = FALSE)
}))

# Patch ERGM optionnel
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
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
if (!exists("InitErgmTerm.dyadcov_GW", mode = "function")) {
  stop("InitErgmTerm.dyadcov_GW introuvable après load_all().")
}

# ======================================================================================
# Réglages de run (point clé du fichier)
# ======================================================================================
RUN <- list(
  phase1_summary = FALSE,
  phase2_fit     = TRUE,
  phase3_mcmc    = FALSE,

  # "quiet" réduit la pollution console des phases 1/2 sans les supprimer.
  quiet_phase1   = FALSE,
  quiet_phase2   = FALSE
)

# ======================================================================================
# Logging local (fichier .log + console)
# ======================================================================================
root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_dyadcov_GW.log")
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
  P1 = c(1L, 1L, 2L, 2L, 3L),   # tailles groupes: 2,2,1
  P2 = c(1L, 2L, 2L, 3L, 3L),   # tailles groupes: 1,2,2
  P3 = c(1L, 1L, 2L, 3L)        # tailles groupes: 2,1,1
)

.make_nodes_df_for_partition <- function(part) {
  n <- length(part)
  data.frame(
    label = paste0("N", seq_len(n)),
    stringsAsFactors = FALSE
  )
}

.make_dyads_for_partition <- function(part) {
  n <- length(part)

  if (n == 5L) {
    Z1 <- matrix(
      c(
        0.0, 1.0, 0.5, 0.3, 0.8,
        1.0, 0.0, 1.2, 0.4, 0.2,
        0.5, 1.2, 0.0, 0.9, 0.6,
        0.3, 0.4, 0.9, 0.0, 1.1,
        0.8, 0.2, 0.6, 1.1, 0.0
      ),
      nrow = 5L, ncol = 5L, byrow = TRUE
    )

    Z2 <- matrix(
      c(
        0.0, 0.3, 0.7, 1.2, 0.5,
        0.4, 0.0, 0.6, 0.9, 1.5,
        1.1, 0.2, 0.0, 0.8, 0.3,
        0.9, 1.4, 0.5, 0.0, 1.0,
        0.2, 0.6, 1.3, 0.7, 0.0
      ),
      nrow = 5L, ncol = 5L, byrow = TRUE
    )

  } else if (n == 4L) {
    Z1 <- matrix(
      c(
        0.0, 0.9, 0.4, 0.7,
        0.9, 0.0, 1.0, 0.2,
        0.4, 1.0, 0.0, 0.6,
        0.7, 0.2, 0.6, 0.0
      ),
      nrow = 4L, ncol = 4L, byrow = TRUE
    )

    Z2 <- matrix(
      c(
        0.0, 0.5, 1.1, 0.3,
        0.8, 0.0, 0.4, 1.0,
        0.2, 1.3, 0.0, 0.9,
        1.2, 0.7, 0.1, 0.0
      ),
      nrow = 4L, ncol = 4L, byrow = TRUE
    )

  } else {
    stop("Taille de partition non supportée dans .make_dyads_for_partition(): n = ", n)
  }

  diag(Z1) <- 0
  diag(Z2) <- 0
  stopifnot(all(diag(Z1) == 0), all(diag(Z2) == 0))

  list(Z1 = Z1, Z2 = Z2)
}

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

  cat("[DEBUG] diag(Z1) =", paste(round(diag(Z1), 3), collapse = ","), "\n")
  cat("[DEBUG] diag(Z2) =", paste(round(diag(Z2), 3), collapse = ","), "\n")
}

# ======================================================================================
# Helpers réseau via builder du wrapper
# ======================================================================================

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

make_formula_for_network_summary <- function(nw, rhs_txt) {
  f <- as.formula(paste0("nw ~ ", rhs_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  f
}

# ======================================================================================
# SUMMARY helpers (réseau vs ERPM traduit)
# ======================================================================================

run_one_network_summary_case_for_dyadcov_GW <- function(partition_vec, nodes_df, dyads_list, rhs_txt) {
  nw <- make_network_from_partition_and_dyads(partition_vec, nodes_df, dyads_list)
  f  <- make_formula_for_network_summary(nw, rhs_txt)
  as.numeric(suppressMessages(summary(f, constraints = ~ b1part)))
}

run_one_erpm_translated_summary_case_for_dyadcov_GW <- function(partition_vec, nodes_df, dyads_list, rhs_txt) {
  partition <- partition_vec
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition, nodes = nodes_df), parent = parent.frame())

  call_ergm <- erpm(
    f,
    eval.call = FALSE,
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

check_summary_equivalence_network_vs_erpm_dyadcov_GW <- function(partition_vec, nodes_df, dyads_list,
                                                                 rhs_vec, tol = 0) {
  ok_all <- TRUE
  for (rhs in rhs_vec) {
    s_net  <- run_one_network_summary_case_for_dyadcov_GW(partition_vec, nodes_df, dyads_list, rhs)
    s_erpm <- run_one_erpm_translated_summary_case_for_dyadcov_GW(partition_vec, nodes_df, dyads_list, rhs)
    cat(sprintf("[SUMMARY-CHECK] n=%-3d RHS=%-40s net=%s  erpm=%s\n",
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
# Panel de cas `dyadcov_GW`
# ======================================================================================

cases_summary <- c(
  "dyadcov_GW('Z1', lambda = 2)",
  "dyadcov_GW('Z1', lambda = 3)",
  "dyadcov_GW('Z2', lambda = 2)",
  "dyadcov_GW('Z2', lambda = 1.5)"
)

# ======================================================================================
# Contrôles ERGM pour fitting
# ======================================================================================

ctrl <- control.ergm(
  # (laisser vide par défaut ici ; le wrapper peut injecter ses réglages)
)

# ======================================================================================
# Phase 1: SUMMARY comparatifs (réseau explicite vs ERPM)
# ======================================================================================

run_phase1_summary_equivalence_checks_dyadcov_GW <- function(quiet = FALSE) {
  cat("=== PHASE 1 : Summary(nw via builder) vs Summary(ERPM-traduit) [dyadcov_GW] ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  total <- 0L; ok <- 0L
  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes_df_for_partition(part)
    dyads <- .make_dyads_for_partition(part)

    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))

    if (!isTRUE(quiet)) print_debug_partition_nodes_dyads(nm, part, nodes, dyads)

    res <- check_summary_equivalence_network_vs_erpm_dyadcov_GW(
      partition_vec = part,
      nodes_df      = nodes,
      dyads_list    = dyads,
      rhs_vec       = cases_summary,
      tol           = 0
    )
    total <- total + length(cases_summary)
    ok    <- ok + as.integer(res) * length(cases_summary)
  }
  cat(sprintf("\n=== Bilan Phase 1 (dyadcov_GW) : %d / %d checks OK ===\n", ok, total))
  if (ok < total) stop(sprintf("Summary mismatch sur %d cas.", total - ok))
  invisible(NULL)
}

# ======================================================================================
# Phase 2: Fits via ERPM
# ======================================================================================

run_one_erpm_fit_with_return_dyadcov_GW <- function(partition_vec, nodes_df, dyads_list,
                                                    rhs_txt, fit_name,
                                                    eval.loglik = TRUE,
                                                    control = ctrl,
                                                    quiet = FALSE) {
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT %-20s] SKIP (erpm() indisponible)\n", fit_name))
    return(list(ok = NA, error = TRUE, coef = NA, fit = NULL, aic = NA, bic = NA))
  }
  set.seed(42)
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition_vec, nodes = nodes_df), parent = parent.frame())

  cat(sprintf("[ERPM-FIT %-20s] n=%-3d RHS=%s  | eval.loglik=%s\n",
              fit_name, length(partition_vec), rhs_txt, as.character(eval.loglik)))

  if (!isTRUE(quiet)) print_debug_partition_nodes_dyads(paste0("FIT_", fit_name), partition_vec, nodes_df, dyads_list)

  fit <- try(
    erpm(
      f,
      eval.loglik = eval.loglik,
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

run_phase2_erpm_fits_and_print_summaries_dyadcov_GW <- function(quiet = FALSE) {
  cat("\n=== PHASE 2 : Fits erpm() ( + logLik ) [dyadcov_GW] ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  fit_specs <- list(
    list(key = "P1_R1", part = partitions$P1, rhs = "dyadcov_GW('Z1', lambda = 2)"),
    list(key = "P1_R2", part = partitions$P1, rhs = "dyadcov_GW('Z2', lambda = 2)"),
    list(key = "P2_R2", part = partitions$P2, rhs = "dyadcov_GW('Z2', lambda = 2)"),
    list(key = "P3_R1", part = partitions$P3, rhs = "dyadcov_GW('Z1', lambda = 2)")
  )

  fit_results <- list()

  for (spec in fit_specs) {
    key   <- spec$key
    part  <- spec$part
    rhs   <- spec$rhs

    nodes <- .make_nodes_df_for_partition(part)
    dyads <- .make_dyads_for_partition(part)

    fit_results[[key]] <- run_one_erpm_fit_with_return_dyadcov_GW(
      partition_vec = part,
      nodes_df      = nodes,
      dyads_list    = dyads,
      rhs_txt       = rhs,
      fit_name      = key,
      eval.loglik   = TRUE,
      quiet         = quiet
    )
  }

  ok_raw <- vapply(fit_results, function(x) x$ok, logical(1))
  n_ok   <- sum(ok_raw, na.rm = TRUE)
  n_tot  <- sum(!is.na(ok_raw))

  cat(sprintf("\n=== Bilan fits erpm() [dyadcov_GW] : %d / %d OK ===\n", n_ok, n_tot))

  cat("\n=== Tableau AIC/BIC pour les fits ERPM (dyadcov_GW) ===\n")
  tab <- data.frame(
    fit = names(fit_results),
    ok  = vapply(fit_results, function(x) isTRUE(x$ok), logical(1)),
    AIC = vapply(fit_results, function(x) as.numeric(x$aic), numeric(1)),
    BIC = vapply(fit_results, function(x) as.numeric(x$bic), numeric(1)),
    stringsAsFactors = FALSE
  )
  print(tab)

  if (!isTRUE(quiet)) {
    cat("\n=== Résumés détaillés des fits ERPM réussis (dyadcov_GW) ===\n")
    for (nm in names(fit_results)) {
      fr <- fit_results[[nm]]
      if (isTRUE(fr$ok) && !is.null(fr$fit)) {
        cat(sprintf("\n--- Résumé fit %s ---\n", nm))
        print(summary(fr$fit))
      }
    }
  }

  if (n_ok < n_tot) stop(sprintf("Echec fits: %d KO", n_tot - n_ok))
  invisible(fit_results)
}

# ======================================================================================
# Phase 3: MCMC multi-toggle probe (diagnostic)
# ======================================================================================

.run_mcmc_multitoggle_probe_dyadcov_GW <- function(nw, dyad_key = "Z1", lambda = 2) {
  # Objectif:
  #   - déclencher la MCMC
  #   - potentiellement observer des pas multi-toggle (ntoggles>1)
  #   - si le C est compilé avec DEBUG_DYADCOV_GW=1, voir passer :
  #       "[dyadcov_GW] MULTI-TOGGLE ntoggles=..."
  #
  # Remarque:
  #   - selon la combinaison (contraintes/proposals), ergm peut rester en 1-toggle.
  #   - ce probe est un “smoke test” pour le chemin D_ changestat + logs.
  #
  # Pour augmenter les chances:
  #   - on utilise une proposal réputée susceptible d'agréger des toggles
  #     (ex: ~ sparse) et on laisse verbose=TRUE.
  ctrl <- control.simulate.formula(
    MCMC.burnin   = 1000,
    MCMC.interval = 1,
    MCMC.prop     = ~ sparse
  )

  rhs <- sprintf("dyadcov_GW('%s', lambda=%s)", dyad_key, format(lambda, digits = 6))
  f <- as.formula(paste0("nw ~ ", rhs))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())

  sim <- simulate(
    f,
    nsim    = 1,
    control = ctrl,
    verbose = TRUE,
    constraints = ~ b1part
  )

  print(sim)
  invisible(sim)
}

run_phase3_mcmc_probe <- function(part_probe, dyads_probe, quiet = FALSE) {
  cat("\n=== PHASE 3: MCMC MULTI-TOGGLE PROBE (dyadcov_GW) ===\n")
  cat("Objectif: déclencher le chemin multi-toggle (D_) et observer les traces C si activées.\n")
  cat("Pour activer les traces C: compiler avec DEBUG_DYADCOV_GW=1 dans changestat_dyadcov_GW.c.\n")
  cat("Note: selon MCMC.prop/contraintes, tu peux ne voir que des 1-toggle.\n\n")

  nodes <- .make_nodes_df_for_partition(part_probe)
  nw    <- make_network_from_partition_and_dyads(part_probe, nodes, dyads_probe)

  if (!isTRUE(quiet)) {
    cat("[PHASE 3] Réseau probe construit. n=", length(part_probe),
        " | groupes=", length(unique(part_probe)),
        " | tailles=", paste(sort(table(part_probe)), collapse=","), "\n", sep = "")
  }

  # Probe 1: Z1, lambda=2
  .run_mcmc_multitoggle_probe_dyadcov_GW(nw, dyad_key = "Z1", lambda = 2)

  cat("\nSi tu vois '[dyadcov_GW] MULTI-TOGGLE ntoggles=...' en console, probe OK.\n")
  invisible(TRUE)
}

# ======================================================================================
# Run principal
# ======================================================================================

run_all_tests_dyadcov_GW <- function() {
  set.seed(1)
  cat("=== TEST ERPM: dyadcov_GW ===\n")
  cat("R:", paste(R.version$major, R.version$minor, sep="."), "\n")
  cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")

  # PHASE 1
  if (isTRUE(RUN$phase1_summary)) {
    run_phase1_summary_equivalence_checks_dyadcov_GW(quiet = isTRUE(RUN$quiet_phase1))
  } else {
    cat("\n=== PHASE 1: SUMMARY ===\nSKIP (désactivée via RUN$phase1_summary = FALSE)\n")
  }

  # PHASE 2
  fit_results <- NULL
  if (isTRUE(RUN$phase2_fit)) {
    fit_results <- run_phase2_erpm_fits_and_print_summaries_dyadcov_GW(quiet = isTRUE(RUN$quiet_phase2))
  } else {
    cat("\n=== PHASE 2: ERPM FIT ===\nSKIP (désactivée via RUN$phase2_fit = FALSE)\n")
  }

  # PHASE 3
  if (isTRUE(RUN$phase3_mcmc)) {
    # On probe sur P1 (n=5) avec dyads hardcodées correspondantes.
    dyads_probe <- .make_dyads_for_partition(partitions$P1)
    run_phase3_mcmc_probe(part_probe = partitions$P1, dyads_probe = dyads_probe, quiet = FALSE)
  } else {
    cat("\n=== PHASE 3: MCMC MULTI-TOGGLE PROBE ===\nSKIP (désactivée via RUN$phase3_mcmc = FALSE)\n")
  }

  invisible(list(fit_results = fit_results))
}

# Exécution quand lancé en script
if (identical(environment(), globalenv())) {
  run_all_tests_dyadcov_GW()
}

# --------------------------------------------------------------------------------------
# Fin de script: on désactive le patch si présent
# --------------------------------------------------------------------------------------
if (exists("ergm_patch_disable", mode = "function")) {
  ergm_patch_disable()
}

cat("\nTous les tests dyadcov_GW ont passé.\n")