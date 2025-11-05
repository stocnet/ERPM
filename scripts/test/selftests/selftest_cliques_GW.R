# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cliques_GW.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cliques_GW`
# Exécution: Rscript scripts/test/selftests/selftest_cliques_GW.R
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

# Patch ERGM optionnel si tu en as besoin pour le projet
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

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
log_path   <- file.path(script_dir, "selftest_cliques_GW.log")

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
# Chargements du package et des utilitaires ERPM
# --------------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else {
  if (!exists("InitErgmTerm.cliques_GW", mode = "function")) {
    stop("InitErgmTerm.cliques_GW introuvable. Exécute depuis le package (devtools::load_all) ou charge le package.")
  }
}

if (!exists("partition_to_bipartite_network", mode = "function")) {
  if (file.exists("R/functions_erpm_bip_network.R")) {
    source("R/functions_erpm_bip_network.R", local = FALSE)
  } else {
    stop("functions_erpm_bip_network.R introuvable. Place ce fichier sous R/.")
  }
}

if (!exists("erpm", mode = "function")) {
  if (file.exists("R/erpm_wrapper.R")) {
    source("R/erpm_wrapper.R", local = FALSE)
  } else {
    cat("[WARN] erpm() indisponible. Les tests de traduction et les fits via erpm() seront sautés.\n")
  }
}

# ======================================================================================
# Fonctions locales
# ======================================================================================

# Construire un biparti depuis une partition
.make_nw_from_partition <- function(part) {
  n <- length(part); stopifnot(n >= 1, is.atomic(part))
  lbl <- utils::head(LETTERS, n)
  partition_to_bipartite_network(labels = lbl, partition = part, attributes = list())
}

# Signature attendue pour la traduction
# Hypothèses par défaut usuelles des GW-terms: k=2, decay=0.5, fixed=TRUE
# Remplace .normalize_gw_signature(...)
.normalize_gw_signature <- function(args = list()) {
  out <- list(lambda = 2)
  if (length(args) && !is.null(args$lambda)) out$lambda <- as.numeric(args$lambda)
  out
}

# Remplace .check_translation_ok(...)
.check_translation_ok <- function(call_ergm, args = list()) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)
  if (!grepl("\\bcliques_GW\\(", compact)) return(FALSE)

  if (is.null(args$lambda)) return(TRUE)  # défaut

  lam <- as.numeric(args$lambda)
  # tolérance de formatage
  fmt1 <- function(x) sub("\\.?0+$","", format(x, trim=TRUE, scientific=FALSE))
  # si vecteur, on exige que chaque valeur apparaisse dans l’appel
  all(vapply(lam, function(v) {
    pat <- paste0("lambda=.*", fmt1(v))
    grepl(pat, compact)
  }, logical(1)))
}
# Exécuter un cas summary + traduction
.run_one_case <- function(part, name, call_txt, args) {
  nw <- .make_nw_from_partition(part)

  f <- as.formula(paste0("nw ~ ", call_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  stat_val <- suppressMessages(as.numeric(summary(f)))  # peut être long>=1

  ok_stat <- isTRUE(all(is.finite(stat_val)))

  ok_trad <- NA
  if (exists("erpm", mode = "function")) {
    call_ergm <- erpm(f, eval_call = FALSE, verbose = TRUE)
    ok_trad <- .check_translation_ok(call_ergm, args = args)
  }

  cat(sprintf("\n[CAS %-20s] part={%s}", name, paste(part, collapse=",")))
  cat(sprintf("\t  appel        : %s", call_txt))
  cat(sprintf("\t  summary(.)   : %s", paste(format(stat_val), collapse=", ")))
  if (!is.na(ok_trad)) cat(sprintf("\t  traduction   : %s", if (ok_trad) "OK" else "KO"))
  cat("\n")

  list(ok_stat = ok_stat, ok_trad = ok_trad, stat = stat_val)
}

# Panel sur une partition
.run_cases_for_partition <- function(part, panel) {
  res <- lapply(panel, function(cx) {
    out <- .run_one_case(part, cx$name, cx$call_txt, cx$args)
    data.frame(
      case    = cx$name,
      ok_stat = out$ok_stat,
      ok_trad = if (is.na(out$ok_trad)) NA else out$ok_trad,
      stat    = I(list(out$stat)),   # préserve le vecteur
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

# ======================================================================================
# Jeu de tests
# ======================================================================================
partitions <- list(
  P1 = c(1, 2, 2, 3, 3, 3),
  P2 = c(1, 1, 2, 3, 3, 4, 4, 4),
  P3 = c(1, 1, 1, 2, 2, 3),
  P4 = c(1, 2, 3, 4, 5),
  P5 = rep(1, 6)
)

# Remplace 'cases <- list(...)'
cases <- list(
  list(name="gw_default",  call_txt="cliques_GW()",                 args=list()),
  list(name="gw_lam2",     call_txt="cliques_GW(lambda=2)",         args=list(lambda=2)),
  list(name="gw_lam1_5",   call_txt="cliques_GW(lambda=1.5)",       args=list(lambda=1.5)),
  list(name="gw_lam_vec",  call_txt="cliques_GW(lambda=c(1.25,4))", args=list(lambda=c(1.25,4)))
)

# ======================================================================================
# Phase 1 : tests summary + traduction
# ======================================================================================
run_all_tests <- function() {
  log_path <- file.path("scripts","test","selftests","selftest_cliques_GW.log")
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
  cat("=== TEST ERPM: cliques_GW ===\n")

  all_results <- list(); total_ok <- 0L; total_n <- 0L
  for (nm in names(partitions)) {
    cat(sprintf("\n--- Partition %s ---\n", nm))
    df <- .run_cases_for_partition(partitions[[nm]], cases)
    all_results[[nm]] <- df

    total_ok <- total_ok + sum(df$ok_stat, na.rm = TRUE)
    total_n  <- total_n  + sum(!is.na(df$ok_stat))

    total_ok <- total_ok + sum(df$ok_trad, na.rm = TRUE)
    total_n  <- total_n  + sum(!is.na(df$ok_trad))

    print(df)
  }

  cat(sprintf("\n=== Bilan global (summary+trad) : %d / %d validations OK ===\n", total_ok, total_n))
  if (total_ok < total_n) stop(sprintf("Echec: %d validations KO", total_n - total_ok))
  invisible(all_results)
}

# ======================================================================================
# Phase 2 : petits fits via erpm() pour valider l'intégration
# ======================================================================================
.erpm_fit_one <- function(part, rhs, name,
                          estimate = "CD",
                          eval.loglik = FALSE,
                          control = list(MCMLE.maxit = 2, MCMC.samplesize = 500),
                          lhs_mode = c("partition","network")) {
  lhs_mode <- match.arg(lhs_mode)
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("\n[ERPM-FIT %-20s] SKIP (erpm() indisponible)\n", name))
    return(list(ok = NA, error = FALSE, coef = NA))
  }

  oldopt <- options(ergm.loglik.warn_dyads = FALSE); on.exit(options(oldopt), add = TRUE)
  set.seed(42)

  if (lhs_mode == "network") {
    nw <- .make_nw_from_partition(part)
    f  <- as.formula(paste0("nw ~ ", rhs))
    environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  } else {
    partition <- part
    f <- as.formula(paste0("partition ~ ", rhs))
    environment(f) <- list2env(list(partition = partition), parent = parent.frame())
  }

  cat(sprintf("\n[ERPM-FIT %-20s] part={%s}\t RHS=%s\t LHS=%s\n",
              name, paste(part, collapse=","), rhs, lhs_mode))

  fit <- try(
    erpm(f, estimate = estimate, eval.loglik = eval.loglik, control = control, verbose = FALSE),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    cat(sprintf("  -> ERREUR (fit): %s\n", as.character(fit)))
    return(list(ok = FALSE, error = TRUE, coef = NA))
  }

  ok_class <- inherits(fit, "ergm")
  cf <- try(stats::coef(fit), silent = TRUE)
  ok_coef <- !inherits(cf, "try-error") && all(is.finite(cf))
  cat(sprintf("  -> class(ergm)? %s | coef finies? %s | coef: %s\n",
              if (ok_class) "OK" else "KO",
              if (ok_coef) "OK" else "KO",
              if (ok_coef) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

  list(ok = ok_class && ok_coef, error = FALSE, coef = if (ok_coef) cf else NA)
}

run_all_tests_with_fits <- function() {
  # Phase 1
  all_results <- run_all_tests()

  # Phase 2
  cat("\n=== PHASE 2 : Fits erpm() courts ===\n")
  fit_results <- list()

  # Remplace les RHS des trois fits
  fit_results[["F1_P1_def"]] <- .erpm_fit_one(
    part   = partitions$P1,
    rhs    = "cliques_GW()",
    name   = "F1_P1_def",
    estimate    = "CD",
    eval.loglik = FALSE,
    control     = list(MCMLE.maxit = 2, MCMC.samplesize = 500),
    lhs_mode    = "partition"
  )

  fit_results[["F2_P2_lam3"]] <- .erpm_fit_one(
    part   = partitions$P2,
    rhs    = "cliques_GW(lambda=3)",
    name   = "F2_P2_lam3",
    estimate    = "CD",
    eval.loglik = FALSE,
    control     = list(MCMLE.maxit = 2, MCMC.samplesize = 500),
    lhs_mode    = "partition"
  )

  fit_results[["F3_P3_lam1"]] <- .erpm_fit_one(
    part   = partitions$P3,
    rhs    = "cliques_GW(lambda=1.1)",
    name   = "F3_P3_lam1.1",
    estimate    = "CD",
    eval.loglik = FALSE,
    control     = list(MCMLE.maxit = 2, MCMC.samplesize = 500),
    lhs_mode    = "partition"
  )

  ok <- vapply(fit_results, function(x) isTRUE(x$ok), logical(1))
  n_ok <- sum(ok, na.rm = TRUE); n_tot <- sum(!is.na(ok))
  cat(sprintf("\n=== Bilan fits erpm() : %d / %d OK ===\n", n_ok, n_tot))

  if (n_ok < n_tot) stop(sprintf("Echec fits: %d KO", n_tot - n_ok))
  invisible(list(summary_results = all_results, fit_results = fit_results))
}

# Point d'entrée
if (identical(environment(), globalenv())) {
  run_all_tests_with_fits()
}

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)