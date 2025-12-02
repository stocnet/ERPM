# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cliques.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cliques`
# Exécution: Rscript scripts/selftests/selftest_cliques.R
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule: environnement et dépendances minimales
# --------------------------------------------------------------------------------------
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

# Patch ERGM optionnel si besoin
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

options(ergm.loglik.warn_dyads = FALSE)

# --------------------------------------------------------------------------------------
# Helpers généraux (chemin script + logging)
# --------------------------------------------------------------------------------------
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
log_path   <- file.path(script_dir, "selftest_cliques.log")

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

# --------------------------------------------------------------------------------------
# Chargement ERPM et wrapper
# --------------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else {
  if (file.exists("R/erpm_wrapper.R")) {
    source("R/erpm_wrapper.R", local = FALSE)
  }
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
# Fonctions locales
# ======================================================================================

# Construire un biparti depuis une partition via le builder du wrapper
make_network_from_partition_via_builder <- function(partition_vec) {
  stopifnot(is.atomic(partition_vec), length(partition_vec) >= 1L)
  built <- build_bipartite_from_inputs(partition = partition_vec)
  built$network
}

# Taille des groupes d'une partition
group_sizes_from_partition <- function(part) as.integer(table(part))

# Valeur de référence analytique de la stat `cliques` pour une partition donnée
# - k == 1 : nombre de groupes de taille 1
# - k >= 2 : sum_g choose(n_g, k)
# - normalized :
#     * k == 1 : identique au cas brut (groupes de taille 1),
#     * k >= 2 : somme_g choose(n_g, k) / n_g (normalisation par taille de groupe).
expected_cliques_from_partition <- function(part, k = 2L, normalized = FALSE) {
  sz <- group_sizes_from_partition(part)
  if (k == 1L) {
    num <- sum(sz == 1L)
    if (!normalized) return(num)
    # Pour k=1, la normalisation par taille de groupe laisse la valeur inchangée.
    return(num)
  } else {
    num_raw <- sum(choose(sz, k))
    if (!normalized) return(num_raw)
    # Nouvelle normalisation par taille de groupe : ∑_g C(n_g, k) / n_g
    contrib <- ifelse(sz > 0L, choose(sz, k) / sz, 0)
    sum(contrib)
  }
}

# Normaliser la signature de cliques(...) pour pilotage summary/erpm
normalize_cliques_signature_args <- function(args = list()) {
  out <- list(k = 2L, normalized = FALSE)
  if (length(args)) {
    nm <- names(args)
    if (!is.null(nm) && length(nm)) {
      if ("k"          %in% nm) out$k          <- as.integer(args[["k"]])
      if ("clique_size"%in% nm) out$k          <- as.integer(args[["clique_size"]])
      if ("normalized" %in% nm) out$normalized <- isTRUE(args[["normalized"]])
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
    checks <- c(checks, grepl(paste0("k=", k_eff), compact, fixed = TRUE) |
                       grepl(paste0("clique_size=", k_eff), compact, fixed = TRUE))
  }

  nrm_eff <- isTRUE(args$normalized)
  if (isTRUE(nrm_eff != def_nrm)) {
    checks <- c(checks, grepl("normalized=TRUE", compact, fixed = TRUE))
  }

  if (!length(checks)) return(TRUE)
  all(checks)
}

# --------------------------------------------------------------------------------------
# SUMMARY + TRADUCTION : un cas
# --------------------------------------------------------------------------------------
run_one_summary_and_translation_case <- function(partition_vec, case_name, rhs_call_text, rhs_args) {
  sig  <- normalize_cliques_signature_args(rhs_args)
  nw   <- make_network_from_partition_via_builder(partition_vec)
  truth <- expected_cliques_from_partition(partition_vec, k = sig$k, normalized = sig$normalized)

  f <- as.formula(paste0("nw ~ ", rhs_call_text))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  stat_val <- suppressMessages(as.numeric(summary(f)))

  ok_summary <- isTRUE(all.equal(unname(as.numeric(stat_val)),
                                 unname(as.numeric(truth)),
                                 tolerance = 1e-10))

  ok_translation <- NA
  if (exists("erpm", mode = "function")) {
    call_ergm <- erpm(f, eval_call = FALSE, verbose = TRUE)
    ok_translation <- check_translation_contains_cliques_with_args(call_ergm, args = rhs_args)
  }

  cat(sprintf("\n[SUMMARY-CASE %-18s] part={%s}", case_name, paste(partition_vec, collapse=",")))
  cat(sprintf("\t  RHS         : %s", rhs_call_text))
  cat(sprintf("\t  summary(.)  : %s", paste(format(stat_val), collapse=", ")))
  cat(sprintf("\t  attendu     : %s", format(truth)))
  if (!is.na(ok_translation)) cat(sprintf("\t  translation : %s", if (ok_translation) "OK" else "KO"))
  cat("\n")

  list(ok_summary = ok_summary, ok_translation = ok_translation, stat = stat_val, expected = truth)
}

# Panel summary() + traduction sur une partition
run_summary_and_translation_panel_for_partition <- function(partition_vec, panel) {
  res <- lapply(panel, function(cx) {
    out <- run_one_summary_and_translation_case(partition_vec, cx$name, cx$call_txt, cx$args)
    data.frame(
      case          = cx$name,
      ok_summary    = out$ok_summary,
      ok_translation= if (is.na(out$ok_translation)) NA else out$ok_translation,
      stat          = I(list(out$stat)),
      expected      = I(list(out$expected)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

# --------------------------------------------------------------------------------------
# ERPM FIT : un cas
# --------------------------------------------------------------------------------------
run_one_erpm_fit <- function(partition_vec, rhs, fit_name,
                             eval.loglik = TRUE,
                             lhs_mode = c("partition","network")) {
  lhs_mode <- match.arg(lhs_mode)
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("\n[ERPM-FIT %-18s] SKIP (erpm() indisponible)\n", fit_name))
    return(list(ok = NA, error = FALSE, coef = NA, fit = NULL))
  }

  oldopt <- options(ergm.loglik.warn_dyads = FALSE); on.exit(options(oldopt), add = TRUE)
  set.seed(42)

  if (lhs_mode == "network") {
    nw <- make_network_from_partition_via_builder(partition_vec)
    f  <- as.formula(paste0("nw ~ ", rhs))
    environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  } else {
    partition <- partition_vec
    f <- as.formula(paste0("partition ~ ", rhs))
    environment(f) <- list2env(list(partition = partition), parent = parent.frame())
  }

  cat(sprintf("\n[ERPM-FIT %-18s] part={%s}\t RHS=%s\t LHS=%s\n",
              fit_name, paste(partition_vec, collapse=","), rhs, lhs_mode))

  fit <- try(
    erpm(f, eval.loglik = eval.loglik, verbose = FALSE),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    cat(sprintf("  -> ERREUR (fit): %s\n", as.character(fit)))
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

# ======================================================================================
# Jeu de tests
# ======================================================================================
partitions <- list(
  P1 = c(1, 2, 2, 3, 3, 3),           # tailles (1,2,3)
  P2 = c(1, 1, 2, 3, 3, 4, 4, 4),     # tailles (2,1,3)
  P3 = c(1, 1, 1, 2, 2, 3),           # tailles (3,2,1)
  P4 = c(1, 2, 3, 4, 5),              # toutes tailles 1 (P4: 5 singletons)
  P5 = rep(1, 6)                      # un seul groupe de taille 6
)

cases <- list(
  list(name="k1_raw",     call_txt="cliques(k=1)",                      args=list(k=1, normalized=FALSE)),
  list(name="k1_norm",    call_txt="cliques(k=1, normalized=TRUE)",     args=list(k=1, normalized=TRUE)),
  list(name="k2_default", call_txt="cliques()",                         args=list(k=2, normalized=FALSE)),
  list(name="k2_exp",     call_txt="cliques(k=2)",                      args=list(k=2, normalized=FALSE)),
  list(name="k2_norm",    call_txt="cliques(k=2, normalized=TRUE)",     args=list(k=2, normalized=TRUE)),
  list(name="k3_raw",     call_txt="cliques(k=3)",                      args=list(k=3, normalized=FALSE)),
  list(name="k3_alias",   call_txt="cliques(clique_size=3)",            args=list(clique_size=3, normalized=FALSE))
)

# ======================================================================================
# Phase 1 : tests summary() + traduction via erpm()
# ======================================================================================
run_all_summary_and_translation_tests <- function() {
  log_path <- file.path("scripts","test","selftests","selftest_cliques.log")
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(log_path)) unlink(log_path)

  sink(file = log_path, split = TRUE)
  con_msg <- file(log_path, open = "at")
  sink(con_msg, type = "message")
  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(close(con_msg),        silent = TRUE)
    try(sink(),                silent = TRUE)
    flush.console()
  }, add = TRUE)

  set.seed(42)
  cat("=== TEST ERPM: cliques — SUMMARY + TRANSLATION ===\n")

  all_results <- list(); total_ok <- 0L; total_n <- 0L
  for (nm in names(partitions)) {
    cat(sprintf("\n--- Partition %s ---\n", nm))
    df <- run_summary_and_translation_panel_for_partition(partitions[[nm]], cases)
    all_results[[nm]] <- df

    total_ok <- total_ok + sum(df$ok_summary,     na.rm = TRUE)
    total_n  <- total_n  + sum(!is.na(df$ok_summary))

    total_ok <- total_ok + sum(df$ok_translation, na.rm = TRUE)
    total_n  <- total_n  + sum(!is.na(df$ok_translation))

    print(df)
  }

  cat(sprintf("\n=== Bilan global (summary+translation) : %d / %d validations OK ===\n", total_ok, total_n))
  if (total_ok < total_n) stop(sprintf("Echec: %d validations KO", total_n - total_ok))
  invisible(all_results)
}

# ======================================================================================
# Phase 2 : petits fits via erpm()
# ======================================================================================
run_all_summary_translation_and_erpm_fit_tests <- function() {
  # Phase 1
  all_results <- run_all_summary_and_translation_tests()

  # Phase 2
  cat("\n=== PHASE 2 : Fits erpm() courts ===\n")
  ctrl <- list(MCMLE.maxit = 10, MCMC.samplesize = 1000)
  fit_results <- list()

  fit_results[["F1_P1_k2"]] <- run_one_erpm_fit(
    partition_vec = partitions$P1,
    rhs       = "cliques()",            # k=2 défaut
    fit_name  = "F1_P1_k2",
    eval.loglik = TRUE,
    lhs_mode  = "partition"
  )

  fit_results[["F2_P2_k3"]] <- run_one_erpm_fit(
    partition_vec = partitions$P2,
    rhs       = "cliques(k=3)",
    fit_name  = "F2_P2_k3",
    eval.loglik = TRUE,
    lhs_mode  = "partition"
  )

  fit_results[["F3_P3_k1n"]] <- run_one_erpm_fit(
    partition_vec = partitions$P3,
    rhs       = "cliques(k=1, normalized=TRUE)",
    fit_name  = "F3_P3_k1n",
    eval.loglik = TRUE,
    lhs_mode  = "partition"
  )

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
  invisible(list(summary_results = all_results, fit_results = fit_results))
}

# Point d'entrée
if (identical(environment(), globalenv())) {
  run_all_summary_translation_and_erpm_fit_tests()
}

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)
