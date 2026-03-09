# ======================================================================================
# Fichier : scripts/test/selftests/selftest_dyadcov_full.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `dyadcov_full` (MULTI-TOGGLE)
# Auteur : Jérémie Chichignoud - Cub'itech
#
# Correctif principal (2026-02-18)
#   - Les checks "summary vs ref" affichaient parfois :
#       summary=40.4 ref=40.4 | ok=FALSE
#     alors que les valeurs étaient égales à l’arrondi.
#   - Cause : comparaison stricte tol=0 + affichage arrondi => faux négatifs (epsilon float).
#   - Fix : tolérance relative/absolue + affichage haute précision + diff explicite.
#
# Bonus
#   - Comptage Phase 1 corrigé : on compte par RHS, pas par bloc booléen global.
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule
# --------------------------------------------------------------------------------------
options(ergm.loglik.warn_dyads = FALSE)
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

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
  if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
}

# Charger le package et le wrapper ERPM
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Le fichier DESCRIPTION n'existe pas ou devtools n'est pas installé.")
}
if (!exists("erpm", mode = "function")) {
  if (file.exists("R/erpm.R")) {
    source("R/erpm.R", local = FALSE)
  } else {
    stop("erpm.R introuvable.")
  }
}
if (!exists("build_bipartite_from_inputs", mode = "function")) {
  stop("build_bipartite_from_inputs() indisponible. Il doit être exporté par R/erpm.R.")
}

# IMPORTANT: activer le debug R-side (initializer) si besoin
# options(ERPM.dyadcov_full.debug = TRUE)

# --------------------------------------------------------------------------------------
# Logging local
# --------------------------------------------------------------------------------------
root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)

log_path <- file.path(root, "scripts", "test", "selftests", "selftest_dyadcov_full.log")
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
# Réglages de run
# ======================================================================================
RUN <- list(
  phase1_summary = TRUE,
  phase2_fit     = TRUE,
  phase3_mcmc    = FALSE,

  quiet_phase1   = FALSE,
  quiet_phase2   = FALSE,

  # Tolérances pour les comparaisons float
  # - abs : pour éviter les faux négatifs autour de 1e-12
  # - rel : stable quand la stat devient grande
  tol_abs = 1e-10,
  tol_rel = 1e-10,

  # PHASE 3:
  mcmc_nsim        = 1,
  mcmc_burnin      = 2000,
  mcmc_interval    = 1,
  mcmc_samplesize  = 20000
)

.maybe_print <- function(x, quiet = FALSE) {
  if (!isTRUE(quiet)) print(x)
  invisible(NULL)
}

# ======================================================================================
# Données de test
# ======================================================================================

# Partitions de test
partitions <- list(
  P1 = c(
    rep(1L, 4),
    rep(2L, 1),
    rep(3L, 3),
    rep(4L, 1)
  ),
  P2 = c(
    rep(1L, 1),
    rep(2L, 3),
    rep(3L, 3),
    rep(4L, 2)
  ),
  P3 = c(
    rep(1L, 2),
    rep(2L, 1),
    rep(3L, 3)
  )
)

# Nodes "muets" pour satisfaire le builder
.make_nodes_df_for_partition <- function(part) {
  n <- length(part)
  data.frame(
    label = paste0("N", seq_len(n)),
    stringsAsFactors = FALSE
  )
}

# Matrices dyadiques déterministes à partir de la partition
# - Z1 : valeurs modérées, plus grandes intra-groupe
# - Z2 : valeurs plus dispersées
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
# Référence R (calcul exact) pour dyadcov_full
# ======================================================================================

# Calcule: sum_g 1[ng in S, ng>=2] * sum_{i!=j, i,j in g} z_ij
# IMPORTANT: la somme est sur des paires ordonnées (i,j), i!=j
.dyadcov_full_ref_from_partition <- function(part, Z, size = NULL) {
  stopifnot(is.matrix(Z), nrow(Z) == length(part), ncol(Z) == length(part))
  S <- NULL
  if (!is.null(size) && length(size)) {
    S <- sort(unique(as.integer(size)))
  }

  groups <- split(seq_along(part), part)
  tot <- 0

  for (ids in groups) {
    ng <- length(ids)
    if (ng < 2) next
    if (!is.null(S) && !(ng %in% S)) next

    s <- 0
    for (p in seq_len(ng)) {
      for (q in seq_len(ng)) {
        if (p == q) next
        s <- s + Z[ids[p], ids[q]]
      }
    }
    tot <- tot + s
  }

  as.numeric(tot)
}

# ======================================================================================
# Fonctions pour SUMMARY et ERPM
# ======================================================================================

run_one_network_summary_case_for_dyadcov_full <- function(partition_vec, nodes_df, dyads_list, rhs_txt) {
  nw <- make_network_from_partition_and_dyads(partition_vec, nodes_df, dyads_list)
  f  <- make_formula_for_network_summary(nw, rhs_txt)
  as.numeric(suppressMessages(summary(f, constraints = ~ b1part)))
}

run_one_erpm_translated_summary_case_for_dyadcov_full <- function(partition_vec, nodes_df, dyads_list, rhs_txt) {
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

  # RHS depuis la formule traduite
  ergm_form <- call_ergm[[2L]]
  rhs_expr  <- ergm_form[[3L]]

  nw2 <- make_network_from_partition_and_dyads(partition_vec, nodes_df, dyads_list)
  f2  <- as.formula(bquote(nw2 ~ .(rhs_expr)))
  environment(f2) <- list2env(list(nw2 = nw2), parent = parent.frame())

  call_args <- as.list(call_ergm)[-1L]
  cons <- call_args$constraints
  if (is.null(cons)) cons <- as.formula(~ b1part)

  as.numeric(suppressMessages(summary(f2, constraints = cons)))
}

# --------------------------------------------------------------------------------------
# Comparaison robuste (tolérance + affichage hi-precision)
# --------------------------------------------------------------------------------------
.is_close_num <- function(x, y, tol_abs, tol_rel) {
  if (!is.finite(x) || !is.finite(y)) return(FALSE)
  d <- abs(x - y)
  d <= (tol_abs + tol_rel * max(1, abs(x), abs(y)))
}

# Check supplémentaire: summary(.) doit matcher une référence R "exacte"
check_summary_vs_reference_R_dyadcov_full <- function(partition_vec, nodes_df, dyads_list, rhs_txt,
                                                     tol_abs = RUN$tol_abs, tol_rel = RUN$tol_rel) {
  # parse RHS minimal (cas supportés: ceux de cases_summary)
  Zname <- if (grepl("'Z2'", rhs_txt, fixed = TRUE)) "Z2" else "Z1"

  size <- NULL
  if (grepl("size", rhs_txt, fixed = TRUE)) {
    tmp <- gsub("^.*size\\s*=\\s*", "", rhs_txt)
    tmp <- gsub("\\).*$", "", tmp)
    size <- eval(parse(text = tmp), envir = baseenv())
  }

  nw <- make_network_from_partition_and_dyads(partition_vec, nodes_df, dyads_list)
  s_net <- as.numeric(suppressMessages(summary(make_formula_for_network_summary(nw, rhs_txt),
                                              constraints = ~ b1part)))

  ref <- .dyadcov_full_ref_from_partition(partition_vec, dyads_list[[Zname]], size = size)

  ok <- .is_close_num(s_net, ref, tol_abs = tol_abs, tol_rel = tol_rel)
  diff <- s_net - ref

  cat(sprintf(
    "[REF-CHECK] n=%-3d RHS=%-30s summary=%.17g ref=%.17g diff=%+.3e | ok=%s (abs=%g rel=%g)\n",
    length(partition_vec), rhs_txt, s_net, ref, diff, ok, tol_abs, tol_rel
  ))

  ok
}

check_summary_equivalence_network_vs_erpm_dyadcov_full <- function(partition_vec, nodes_df, dyads_list,
                                                                   rhs_vec,
                                                                   tol_abs = RUN$tol_abs, tol_rel = RUN$tol_rel) {
  ok_flags <- logical(length(rhs_vec))
  for (i in seq_along(rhs_vec)) {
    rhs <- rhs_vec[[i]]
    s_net  <- run_one_network_summary_case_for_dyadcov_full(partition_vec, nodes_df, dyads_list, rhs)
    s_erpm <- run_one_erpm_translated_summary_case_for_dyadcov_full(partition_vec, nodes_df, dyads_list, rhs)

    cat(sprintf("[SUMMARY-CHECK] n=%-3d RHS=%-30s net=%.17g  erpm=%.17g diff=%+.3e\n",
                length(partition_vec), rhs, s_net, s_erpm, (s_net - s_erpm)))

    ok_flags[[i]] <- .is_close_num(s_net, s_erpm, tol_abs = tol_abs, tol_rel = tol_rel)
    if (!ok_flags[[i]]) cat("  -> MISMATCH détecté.\n")
  }
  ok_flags
}

# ======================================================================================
# Panel de cas `dyadcov_full`
# ======================================================================================

cases_summary <- c(
  "dyadcov_full('Z1')",
  "dyadcov_full('Z1', size = 2)",
  "dyadcov_full('Z1', size = 2:3)",
  "dyadcov_full('Z2')",
  "dyadcov_full('Z2', size = 3)"
)

# ======================================================================================
# Contrôles ERGM pour fitting
# ======================================================================================

ctrl <- control.ergm(MCMLE.maxit = 10, MCMC.samplesize = 1e4)

# ======================================================================================
# Phase 1: SUMMARY comparatifs
# ======================================================================================

run_phase1_summary_equivalence_checks_dyadcov_full <- function() {
  cat("=== PHASE 1 : Summary(nw via builder) vs Summary(ERPM-traduit) + Référence R [dyadcov_full] ===\n")
  total <- 0L
  ok    <- 0L

  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes_df_for_partition(part)
    dyads <- .make_dyads_for_partition(part)

    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse = ",")))

    print_debug_partition_nodes_dyads(nm, part, nodes, dyads)

    # (A) net vs erpm traduit (compte par RHS)
    flagsA <- check_summary_equivalence_network_vs_erpm_dyadcov_full(
      partition_vec = part,
      nodes_df      = nodes,
      dyads_list    = dyads,
      rhs_vec       = cases_summary
    )
    total <- total + length(flagsA)
    ok    <- ok    + sum(flagsA)

    # (B) net vs référence R (compte par RHS)
    for (rhs in cases_summary) {
      total <- total + 1L
      ok    <- ok + as.integer(check_summary_vs_reference_R_dyadcov_full(part, nodes, dyads, rhs))
    }
  }

  cat(sprintf("\n=== Bilan Phase 1 (dyadcov_full) : %d / %d checks OK ===\n", ok, total))
  if (ok < total) stop(sprintf("Phase 1 KO: %d mismatch.", total - ok))
  invisible(NULL)
}

# ======================================================================================
# Phase 2: Fits via ERPM (inchangé)
# ======================================================================================

run_one_erpm_fit_with_return_dyadcov_full <- function(partition_vec, nodes_df, dyads_list,
                                                     rhs_txt, fit_name,
                                                     estimate    = NULL,
                                                     eval.loglik = TRUE,
                                                     control     = ctrl) {
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT %-20s] SKIP (erpm() indisponible)\n", fit_name))
    return(list(ok = NA, error = TRUE, coef = NA, fit = NULL, aic = NA, bic = NA))
  }

  set.seed(42)
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition_vec, nodes = nodes_df), parent = parent.frame())

  cat(sprintf("[ERPM-FIT %-20s] n=%-3d RHS=%s  | estimate=%s eval.loglik=%s\n",
              fit_name, length(partition_vec), rhs_txt,
              as.character(estimate), as.character(eval.loglik)))

  print_debug_partition_nodes_dyads(paste0("FIT_", fit_name), partition_vec, nodes_df, dyads_list)

  fit <- try(
    erpm(
      f,
      eval.loglik = eval.loglik,
      verbose     = TRUE,
      nodes       = nodes_df,
      dyads       = dyads_list,
      control     = control
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
  ok_coef <- !inherits(cf, "try-error") && length(cf) > 0L && all(is.finite(cf))

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
              if (ok_coef) paste(format(as.numeric(cf)), collapse = ", ") else "NA",
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

run_phase2_erpm_fits_and_print_summaries_dyadcov_full <- function() {
  cat("\n=== PHASE 2 : Fits erpm() ( + logLik ) [dyadcov_full] ===\n")

  rhs_list <- list(
    R1 = "dyadcov_full('Z1', size = 2:3)",
    R2 = "dyadcov_full('Z2', size = 3)"
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
      key <- paste0(px$name, "_", nm)
      fit_results[[key]] <- run_one_erpm_fit_with_return_dyadcov_full(
        partition_vec = px$part,
        nodes_df      = nodes,
        dyads_list    = dyads,
        rhs_txt       = rhs_list[[nm]],
        fit_name      = key,
        eval.loglik   = TRUE
      )
    }
  }

  ok_raw <- vapply(fit_results, function(x) x$ok, logical(1))
  n_ok   <- sum(ok_raw, na.rm = TRUE)
  n_tot  <- sum(!is.na(ok_raw))

  cat(sprintf("\n=== Bilan fits erpm() [dyadcov_full] : %d / %d OK ===\n", n_ok, n_tot))

  cat("\n=== Tableau AIC/BIC pour les fits ERPM (dyadcov_full) ===\n")
  tab <- data.frame(
    fit = names(fit_results),
    ok  = vapply(fit_results, function(x) isTRUE(x$ok), logical(1)),
    AIC = vapply(fit_results, function(x) x$aic, numeric(1)),
    BIC = vapply(fit_results, function(x) x$bic, numeric(1)),
    stringsAsFactors = FALSE
  )
  print(tab)

  cat("\n=== Résumés détaillés des fits ERPM réussis (dyadcov_full) ===\n")
  for (nm in names(fit_results)) {
    fit_obj <- fit_results[[nm]]
    if (isTRUE(fit_obj$ok) && !is.null(fit_obj$fit)) {
      cat(sprintf("\n--- Résumé fit %s ---\n", nm))
      print(summary(fit_obj$fit))
    }
  }

  if (n_ok < n_tot) stop(sprintf("Phase 2 KO: %d fits en échec.", n_tot - n_ok))
  invisible(fit_results)
}

# ======================================================================================
# Phase 3: MCMC multi-toggle probe (inchangé)
# ======================================================================================

.run_mcmc_multitoggle_probe <- function(nw, rhs_txt, cons, prop, label) {
  cat("\n--- PROBE:", label, "---\n")
  cat("RHS:", rhs_txt, "\n")
  cat("constraints:", deparse(cons), "\n")
  cat("MCMC.prop:", deparse(prop), "\n")
  cat("NOTE: si DEBUG_DYADCOV_FULL=1, tu dois voir côté C un message type:\n")
  cat("      [dyadcov_full] MULTI-TOGGLE ntoggles=...\n")

  ctrl <- control.simulate.formula(
    MCMC.burnin   = RUN$mcmc_burnin,
    MCMC.interval = RUN$mcmc_interval,
    MCMC.prop     = prop
  )

  f <- as.formula(paste0("nw ~ ", rhs_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())

  sim <- simulate(
    f,
    nsim        = RUN$mcmc_nsim,
    control     = ctrl,
    constraints = cons,
    verbose     = TRUE
  )

  print(sim)
  invisible(sim)
}

run_phase3_mcmc_multitoggle_dyadcov_full <- function() {
  cat("\n=== PHASE 3 : MCMC MULTI-TOGGLE PROBE [dyadcov_full] ===\n")
  cat("Objectif: déclencher ntoggles>1 dans le D_CHANGESTAT_FN.\n")
  cat("Si tu ne vois rien:\n")
  cat("  - DEBUG_DYADCOV_FULL est probablement à 0 (recompile nécessaire)\n")
  cat("  - ou bien le MCMC.prop choisi ne produit pas de multi-toggle sur ton setup\n")

  part  <- partitions$P1
  nodes <- .make_nodes_df_for_partition(part)
  dyads <- .make_dyads_for_partition(part)
  nw    <- make_network_from_partition_and_dyads(part, nodes, dyads)

  rhs  <- "dyadcov_full('Z1', size = 2:3)"
  cons <- as.formula(~ b1part)

  probes <- list(
    list(prop = ~ sparse,   label = "prop=~sparse"),
    list(prop = ~ default,  label = "prop=~default")
  )

  for (p in probes) {
    .run_mcmc_multitoggle_probe(nw, rhs, cons, p$prop, p$label)
  }

  cat("\nFIN PHASE 3.\n")
  cat("Attendu (si debug C activé): au moins un 'MULTI-TOGGLE ntoggles=...'.\n")
  invisible(TRUE)
}

# ======================================================================================
# Exécution
# ======================================================================================

run_all_tests_dyadcov_full <- function() {
  set.seed(1)
  cat("=== TEST ERPM: dyadcov_full (MULTI-TOGGLE) ===\n")
  cat("R:", paste(R.version$major, R.version$minor, sep = "."), "\n")
  cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")

  if (isTRUE(RUN$phase1_summary)) {
    run_phase1_summary_equivalence_checks_dyadcov_full()
  } else {
    cat("\n=== PHASE 1 : SUMMARY ===\nSKIP\n")
  }

  fit_results <- NULL
  if (isTRUE(RUN$phase2_fit)) {
    fit_results <- run_phase2_erpm_fits_and_print_summaries_dyadcov_full()
  } else {
    cat("\n=== PHASE 2 : FIT ===\nSKIP\n")
  }

  if (isTRUE(RUN$phase3_mcmc)) {
    run_phase3_mcmc_multitoggle_dyadcov_full()
  } else {
    cat("\n=== PHASE 3 : MCMC MULTI-TOGGLE PROBE ===\nSKIP\n")
  }

  invisible(list(fit_results = fit_results))
}

if (identical(environment(), globalenv())) {
  run_all_tests_dyadcov_full()
}

on.exit(try(if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable(), silent = TRUE), add = TRUE)
cat("\nTous les tests dyadcov_full ont passé.\n")