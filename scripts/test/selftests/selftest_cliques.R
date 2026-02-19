# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cliques.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cliques`
# Exécution: Rscript scripts/test/selftests/selftest_cliques.R
#
# But du fichier
#   - PHASE 1 (SUMMARY + TRANSLATION) : valider la statistique via summary(nw ~ cliques(...))
#     et vérifier la traduction erpm(...) -> ergm(...) (eval.call=FALSE).
#   - PHASE 2 (ERPM FIT) : valider qu'un fit erpm() passe et renvoie des coefs finis.
#   - PHASE 3 (MCMC)     : diagnostic "multi-toggle" (déclencher le chemin D_CHANGESTAT_FN).
#
# IMPORTANT (multi-toggle):
#   - Le terme `cliques` est maintenant un D_CHANGESTAT_FN (multi-toggle).
#   - Le test PHASE 3 vise à déclencher des propositions à plusieurs toggles.
#   - Pour VOIR des traces côté C, il faut compiler avec DEBUG_CLIQUES=1 dans
#     changestat_cliques.c, puis recompiler/reinstaller le package.
#
# Notes
#   - Les phases 1/2 peuvent spammer la console (prints/summaries).
#   - On peut réduire la sortie via RUN$quiet_phase1 / RUN$quiet_phase2.
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule: environnement et dépendances minimales
# --------------------------------------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",    quietly = TRUE)) stop("Package 'ergm' requis.")
})

suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

# Patch ERGM optionnel si besoin
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
}

# --------------------------------------------------------------------------------------
# Chargement ERPM et wrapper
# --------------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else {
  if (file.exists("R/erpm_wrapper.R")) source("R/erpm_wrapper.R", local = FALSE)
}

# Vérifs minimales des symboles requis
if (!exists("InitErgmTerm.cliques", mode = "function")) {
  stop("InitErgmTerm.cliques introuvable. Charge le package (devtools::load_all) ou assure-toi que le fichier R est présent.")
}
if (!exists("erpm", mode = "function")) {
  stop("erpm() indisponible. Charge le wrapper via devtools::load_all ou source('R/erpm_wrapper.R').")
}
if (!exists("build_bipartite_from_inputs", mode = "function")) {
  stop("build_bipartite_from_inputs() indisponible. Il doit être exporté par R/erpm_wrapper.R.")
}

# ======================================================================================
# Réglages de run (le point clé du fichier)
# ======================================================================================
RUN <- list(
  phase1_summary = TRUE,   # summary + référence + traduction erpm(eval.call=FALSE)
  phase2_fit     = TRUE,   # petits fits erpm()
  phase3_mcmc    = FALSE,   # probe MCMC pour déclencher multi-toggle

  quiet_phase1   = FALSE,
  quiet_phase2   = FALSE
)

# ======================================================================================
# Helpers généraux (chemin script + logging)
# ======================================================================================
.get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 1L) return(normalizePath(dirname(f), winslash = "/", mustWork = FALSE))
  if (!is.null(sys.frames()) && !is.null(sys.calls())) {
    for (i in rev(seq_along(sys.calls()))) {
      cf <- sys.frame(i)
      if (!is.null(cf$ofile)) return(normalizePath(dirname(cf$ofile), winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- .get_script_dir()
log_dir    <- file.path("scripts", "test", "selftests")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_path   <- file.path(log_dir, "selftest_cliques.log")

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

cat("==> Log: ", log_path, "\n")

# ======================================================================================
# Fonctions locales — utilitaires de référence
# ======================================================================================

# Construire un biparti depuis une partition via le builder du wrapper
make_network_from_partition_via_builder <- function(partition_vec) {
  stopifnot(is.atomic(partition_vec), length(partition_vec) >= 1L)
  built <- build_bipartite_from_inputs(partition = partition_vec)

  # Tolérance: selon versions, le builder peut renvoyer directement un network
  if (inherits(built, "network")) return(built)
  if (is.list(built) && !is.null(built$network) && inherits(built$network, "network")) return(built$network)

  stop("build_bipartite_from_inputs() n'a pas renvoyé un 'network' (ni built$network).")
}

# Taille des groupes d'une partition
group_sizes_from_partition <- function(part) as.integer(table(part))

# Valeur de référence analytique de la stat `cliques` pour une partition donnée
# - k == 1 : nombre de groupes de taille 1
# - k >= 2 : sum_g choose(n_g, k)
# - normalized :
#     * k == 1 : identique au cas brut,
#     * k >= 2 : somme_g choose(n_g, k) / n_g
expected_cliques_from_partition <- function(part, k = 2L, normalized = FALSE) {
  sz <- group_sizes_from_partition(part)

  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L) stop("'k' invalide pour expected_cliques_from_partition().")

  if (k == 1L) {
    num <- sum(sz == 1L)
    return(num)
  }

  num_raw <- sum(choose(sz, k))
  if (!isTRUE(normalized)) return(num_raw)

  contrib <- ifelse(sz > 0L, choose(sz, k) / sz, 0)
  sum(contrib)
}

# Normaliser la signature de cliques(...) pour pilotage summary/erpm
normalize_cliques_signature_args <- function(args = list()) {
  out <- list(k = 2L, normalized = FALSE)
  if (length(args)) {
    nm <- names(args)
    if (!is.null(nm) && length(nm)) {
      if ("k"           %in% nm) out$k          <- as.integer(args[["k"]])
      if ("clique_size" %in% nm) out$k          <- as.integer(args[["clique_size"]])
      if ("normalized"  %in% nm) out$normalized <- isTRUE(args[["normalized"]])
    }
  }
  out
}

# Vérifier que la traduction erpm(...) → ergm(...) contient bien `cliques(...)` avec les options pertinentes
check_translation_contains_cliques_with_args <- function(call_ergm, args = list()) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)

  if (!grepl("\\bcliques\\(", compact)) return(FALSE)

  def_k   <- 2L
  def_nrm <- FALSE

  checks <- logical(0)

  k_eff <- if (!is.null(args$clique_size)) as.integer(args$clique_size)
           else if (!is.null(args$k))      as.integer(args$k)
           else def_k
  if (!identical(k_eff, def_k)) {
    checks <- c(checks,
                grepl(paste0("k=", k_eff), compact, fixed = TRUE) |
                grepl(paste0("clique_size=", k_eff), compact, fixed = TRUE))
  }

  nrm_eff <- isTRUE(args$normalized)
  if (isTRUE(nrm_eff != def_nrm)) {
    checks <- c(checks, grepl("normalized=TRUE", compact, fixed = TRUE))
  }

  if (!length(checks)) return(TRUE)
  all(checks)
}

.maybe_print <- function(x, quiet = FALSE) {
  if (!isTRUE(quiet)) print(x)
  invisible(NULL)
}

# ======================================================================================
# PHASE 1 — SUMMARY + TRADUCTION
# ======================================================================================

run_one_summary_and_translation_case <- function(partition_vec, case_name, rhs_call_text, rhs_args, quiet = FALSE) {
  sig   <- normalize_cliques_signature_args(rhs_args)
  nw    <- make_network_from_partition_via_builder(partition_vec)
  truth <- expected_cliques_from_partition(partition_vec, k = sig$k, normalized = sig$normalized)

  f <- as.formula(paste0("nw ~ ", rhs_call_text))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  stat_val <- suppressMessages(as.numeric(summary(f)))

  ok_summary <- isTRUE(all.equal(unname(as.numeric(stat_val)),
                                 unname(as.numeric(truth)),
                                 tolerance = 1e-10))

  ok_translation <- NA
  if (exists("erpm", mode = "function")) {
    call_ergm <- erpm(f, eval.call = FALSE, verbose = FALSE)
    ok_translation <- check_translation_contains_cliques_with_args(call_ergm, args = rhs_args)
  }

  if (!isTRUE(quiet)) {
    cat(sprintf("\n[SUMMARY-CASE %-18s] part={%s}\n", case_name, paste(partition_vec, collapse=",")))
    cat(sprintf("  RHS         : %s\n", rhs_call_text))
    cat(sprintf("  summary(.)  : %s\n", paste(format(stat_val), collapse=", ")))
    cat(sprintf("  attendu     : %s\n", format(truth)))
    if (!is.na(ok_translation)) cat(sprintf("  translation : %s\n", if (ok_translation) "OK" else "KO"))
  }

  list(ok_summary = ok_summary, ok_translation = ok_translation, stat = stat_val, expected = truth)
}

run_summary_and_translation_panel_for_partition <- function(partition_vec, panel, quiet = FALSE) {
  res <- lapply(panel, function(cx) {
    out <- run_one_summary_and_translation_case(partition_vec, cx$name, cx$call_txt, cx$args, quiet = quiet)
    data.frame(
      case           = cx$name,
      ok_summary     = out$ok_summary,
      ok_translation = if (is.na(out$ok_translation)) NA else out$ok_translation,
      stat           = I(list(out$stat)),
      expected       = I(list(out$expected)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

run_all_summary_and_translation_tests <- function(partitions, cases, quiet = FALSE) {
  cat("\n=== PHASE 1: SUMMARY + TRANSLATION ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  all_results <- list()
  total_ok <- 0L
  total_n  <- 0L

  for (nm in names(partitions)) {
    cat(sprintf("\n--- Partition %s ---\n", nm))
    df <- run_summary_and_translation_panel_for_partition(partitions[[nm]], cases, quiet = quiet)
    all_results[[nm]] <- df

    total_ok <- total_ok + sum(df$ok_summary,     na.rm = TRUE)
    total_n  <- total_n  + sum(!is.na(df$ok_summary))

    total_ok <- total_ok + sum(df$ok_translation, na.rm = TRUE)
    total_n  <- total_n  + sum(!is.na(df$ok_translation))

    .maybe_print(df, quiet = quiet)
  }

  cat(sprintf("\n=== Bilan global (summary+translation) : %d / %d validations OK ===\n", total_ok, total_n))
  if (total_ok < total_n) stop(sprintf("Echec: %d validations KO", total_n - total_ok))
  invisible(all_results)
}

# ======================================================================================
# PHASE 2 — ERPM FIT : petits fits
# ======================================================================================

run_one_erpm_fit <- function(partition_vec, rhs, fit_name, eval.loglik = TRUE) {
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("\n[ERPM-FIT %-18s] SKIP (erpm() indisponible)\n", fit_name))
    return(list(ok = NA, error = FALSE, coef = NA, fit = NULL))
  }

  set.seed(42)
  partition <- partition_vec
  f <- as.formula(paste0("partition ~ ", rhs))
  environment(f) <- list2env(list(partition = partition), parent = parent.frame())

  cat(sprintf("\n[ERPM-FIT %-18s] part={%s}\n  RHS=%s\n", fit_name, paste(partition_vec, collapse=","), rhs))

  fit <- try(erpm(f, eval.loglik = eval.loglik, verbose = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) {
    cat(sprintf("  -> ERREUR (fit): %s\n", as.character(fit)[1]))
    return(list(ok = FALSE, error = TRUE, coef = NA, fit = NULL))
  }

  ok_class <- inherits(fit, "ergm")
  cf <- try(stats::coef(fit), silent = TRUE)
  ok_coef <- !inherits(cf, "try-error") && length(cf) > 0L && all(is.finite(cf))

  cat(sprintf("  -> class(ergm)? %s | coef finies? %s | coef: %s\n",
              if (ok_class) "OK" else "KO",
              if (ok_coef)  "OK" else "KO",
              if (ok_coef)  paste(format(as.numeric(cf)), collapse=", ") else "NA"))

  list(ok = ok_class && ok_coef, error = FALSE, coef = if (ok_coef) cf else NA, fit = fit)
}

run_all_erpm_fit_tests <- function(partitions, quiet = FALSE) {
  cat("\n=== PHASE 2: ERPM FIT (sanity check) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  fits <- list(
    F1_P1_k2  = run_one_erpm_fit(partitions$P1, "cliques()",                    "F1_P1_k2"),
    F2_P2_k3  = run_one_erpm_fit(partitions$P2, "cliques(k=3)",                 "F2_P2_k3"),
    F3_P3_k1n = run_one_erpm_fit(partitions$P3, "cliques(k=1, normalized=TRUE)", "F3_P3_k1n")
  )

  ok <- vapply(fits, function(x) if (is.na(x$ok)) NA else isTRUE(x$ok), logical(1))
  n_ok <- sum(ok, na.rm = TRUE)
  n_tot <- sum(!is.na(ok))

  cat(sprintf("\n=== Bilan fits erpm() : %d / %d OK ===\n", n_ok, n_tot))
  if (n_ok < n_tot) stop(sprintf("Echec fits: %d KO", n_tot - n_ok))

  if (!isTRUE(quiet)) {
    cat("\n=== Résumés détaillés des fits ERPM réussis ===\n")
    for (nm in names(fits)) {
      if (isTRUE(fits[[nm]]$ok) && !is.null(fits[[nm]]$fit)) {
        cat(sprintf("\n--- Résumé fit %s ---\n", nm))
        print(summary(fits[[nm]]$fit))
      }
    }
  }

  invisible(fits)
}

# ======================================================================================
# PHASE 3 — MCMC multi-toggle probe (diagnostic)
# ======================================================================================

.run_mcmc_multitoggle_probe <- function(nw) {
  # Objectif:
  #   - déclencher le code MCMC
  #   - observer (si DEBUG_CLIQUES=1 côté C) des traces type:
  #       "[cliques] MULTI-TOGGLE ntoggles=..."
  #
  # Important: on garde verbose=TRUE exprès.

  ctrl <- control.simulate.formula(
    MCMC.burnin   = 1000,
    MCMC.interval = 1,
    MCMC.prop     = ~ sparse
  )

  f <- nw ~ cliques()

  sim <- simulate(
    f,
    nsim    = 1,
    control = ctrl,
    verbose = TRUE
  )

  print(sim)
  invisible(sim)
}

run_phase3_mcmc <- function(part_probe) {
  cat("\n=== PHASE 3: MCMC MULTI-TOGGLE PROBE ===\n")
  cat("Objectif: déclencher des propositions multi-toggle (D_CHANGESTAT_FN).\n")
  cat("Pour voir des traces côté C: compiler changestat_cliques.c avec DEBUG_CLIQUES=1.\n\n")

  nw_probe <- make_network_from_partition_via_builder(part_probe)
  .run_mcmc_multitoggle_probe(nw_probe)

  cat("\nSi '[cliques] MULTI-TOGGLE ntoggles=...' apparaît (DEBUG_CLIQUES=1), test OK.\n")
  invisible(TRUE)
}

# ======================================================================================
# Jeu de tests
# ======================================================================================
partitions <- list(
  P1 = c(1, 2, 2, 3, 3, 3),           # tailles (1,2,3)
  P2 = c(1, 1, 2, 3, 3, 4, 4, 4),     # tailles (2,1,3)
  P3 = c(1, 1, 1, 2, 2, 3),           # tailles (3,2,1)
  P4 = c(1, 2, 3, 4, 5),              # 5 singletons
  P5 = rep(1, 6)                      # un seul groupe taille 6
)

cases <- list(
  list(name="k1_raw",     call_txt="cliques(k=1)",                  args=list(k=1, normalized=FALSE)),
  list(name="k1_norm",    call_txt="cliques(k=1, normalized=TRUE)", args=list(k=1, normalized=TRUE)),
  list(name="k2_default", call_txt="cliques()",                     args=list(k=2, normalized=FALSE)),
  list(name="k2_exp",     call_txt="cliques(k=2)",                  args=list(k=2, normalized=FALSE)),
  list(name="k2_norm",    call_txt="cliques(k=2, normalized=TRUE)", args=list(k=2, normalized=TRUE)),
  list(name="k3_raw",     call_txt="cliques(k=3)",                  args=list(k=3, normalized=FALSE)),
  list(name="k3_alias",   call_txt="cliques(clique_size=3)",        args=list(clique_size=3, normalized=FALSE))
)

# ======================================================================================
# Run principal
# ======================================================================================

run_all_tests_cliques <- function() {
  set.seed(42)

  cat("=== SELFTEST ERPM: cliques ===\n")
  cat("R:", paste(R.version$major, R.version$minor, sep="."), "\n")
  cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")

  # PHASE 1
  summary_results <- NULL
  if (isTRUE(RUN$phase1_summary)) {
    summary_results <- run_all_summary_and_translation_tests(
      partitions = partitions,
      cases      = cases,
      quiet      = isTRUE(RUN$quiet_phase1)
    )
  } else {
    cat("\n=== PHASE 1: SUMMARY + TRANSLATION ===\n")
    cat("SKIP (désactivée via RUN$phase1_summary = FALSE)\n")
  }

  # PHASE 2
  fit_results <- NULL
  if (isTRUE(RUN$phase2_fit)) {
    fit_results <- run_all_erpm_fit_tests(
      partitions = partitions,
      quiet      = isTRUE(RUN$quiet_phase2)
    )
  } else {
    cat("\n=== PHASE 2: ERPM FIT ===\n")
    cat("SKIP (désactivée via RUN$phase2_fit = FALSE)\n")
  }

  # PHASE 3
  if (isTRUE(RUN$phase3_mcmc)) {
    run_phase3_mcmc(part_probe = partitions$P1)
  } else {
    cat("\n=== PHASE 3: MCMC MULTI-TOGGLE PROBE ===\n")
    cat("SKIP (désactivée via RUN$phase3_mcmc = FALSE)\n")
  }

  invisible(list(summary_results = summary_results, fit_results = fit_results))
}

if (identical(environment(), globalenv())) {
  run_all_tests_cliques()
}

# --------------------------------------------------------------------------------------
# Fin de script: on désactive le patch si présent
# --------------------------------------------------------------------------------------
if (exists("ergm_patch_disable", mode = "function")) {
  ergm_patch_disable()
}