# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cliques_GW.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cliques_GW`
# Exécution: Rscript scripts/test/selftests/selftest_cliques_GW.R
#
# But du fichier
#   - PHASE 1 (SUMMARY) : valider la statistique via summary(nw ~ cliques_GW(...)).
#   - PHASE 2 (ERPM FIT): valider que erpm() construit un modèle et renvoie des coefs finis.
#   - PHASE 3 (MCMC)    : diagnostic "multi-toggle" (le point clé ici) pour vérifier
#                         que le changestat D_ est bien appelé en MCMC et qu'on
#                         ne segfault pas.
#
# Important
#   - Les phases 1/2 peuvent spammer la console (summaries, prints).
#   - Pour bosser sur la phase 3, on peut désactiver 1/2 via des flags.
#   - Pour voir les traces côté C:
#       * activer DEBUG_CLIQUES_GW=1 dans changestat_cliques_GW.c
#       * recompiler (INSTALL / devtools::load_all)
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

# Patch ERGM (optionnel)
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) {
    ergm_patch_enable()
  }
} else {
  message("[ergm_patch] scripts/ergm_patch.R introuvable, on continue sans patch")
}

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
if (!exists("InitErgmTerm.cliques_GW", mode = "function")) {
  stop("InitErgmTerm.cliques_GW introuvable. Charger le package (devtools::load_all) ou vérifier le fichier R.")
}
if (!exists("erpm", mode = "function")) {
  stop("erpm() indisponible. Charger le wrapper via devtools::load_all ou source('R/erpm_wrapper.R').")
}
if (!exists("build_bipartite_from_inputs", mode = "function")) {
  stop("build_bipartite_from_inputs() indisponible. Il doit être exposé par le wrapper.")
}

# ======================================================================================
# Réglages de run (le point clé du fichier)
# ======================================================================================
RUN <- list(
  phase1_summary = TRUE,   # TRUE = on valide summary() ; FALSE = on saute
  phase2_fit     = TRUE,   # TRUE = on fit un mini-modèle ; FALSE = on saute
  phase3_mcmc    = FALSE,   # TRUE = on lance le probe multi-toggle ; FALSE = on saute

  quiet_phase1   = FALSE,
  quiet_phase2   = FALSE
)

# ======================================================================================
# Helpers généraux (script dir + logging)
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

# ======================================================================================
# Fonctions locales — utilitaires de référence
# ======================================================================================

# Construire un biparti depuis une partition via le builder du wrapper
.make_network_from_partition_via_builder <- function(partition_vec) {
  stopifnot(is.atomic(partition_vec), length(partition_vec) >= 1L)
  built <- build_bipartite_from_inputs(partition = partition_vec)
  if (inherits(built, "network")) return(built)
  if (is.list(built) && !is.null(built$network) && inherits(built$network, "network")) return(built$network)

  # Fallback: essayer quelques clés usuelles si ton wrapper a varié
  if (is.list(built)) {
    for (nm in c("nw","net","graph","g","bip")) {
      if (!is.null(built[[nm]]) && inherits(built[[nm]], "network")) return(built[[nm]])
    }
  }
  stop("build_bipartite_from_inputs() n'a pas renvoyé de network exploitable.")
}

# Référence R: T_lambda(y) = sum_g lambda * (1 - r^deg(g)), r=(lambda-1)/lambda
.expected_T_lambda_from_partition <- function(part, lambda) {
  stopifnot(length(lambda) >= 1L)

  # group sizes from partition (counts)
  sizes <- as.integer(table(part))
  out <- numeric(length(lambda))

  for (j in seq_along(lambda)) {
    lam <- as.numeric(lambda[j])
    r   <- (lam - 1) / lam
    # contribution per group: lam*(1 - r^d)
    out[j] <- sum(lam * (1 - (r ^ sizes)))
  }
  out
}

# Vérification légère: la traduction erpm(...) contient cliques_GW et lambda=...
.check_translation_contains_cliques_GW_with_args <- function(call_ergm, args = list()) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)
  if (!grepl("\\bcliques_GW\\(", compact)) return(FALSE)

  if (is.null(args$lambda)) return(TRUE)

  lam <- as.numeric(args$lambda)
  fmt <- function(x) sub("\\.?0+$","", format(x, trim=TRUE, scientific=FALSE))
  all(vapply(lam, function(v) grepl(fmt(v), compact, fixed = TRUE), logical(1)))
}

.maybe_print <- function(x, quiet = FALSE) {
  if (!isTRUE(quiet)) print(x)
  invisible(NULL)
}

# ======================================================================================
# PHASE 1 — SUMMARY : summary(.) vs référence R (partition)
# ======================================================================================

.run_one_case_summary <- function(part, case, quiet = FALSE) {
  nw <- .make_network_from_partition_via_builder(part)

  f <- as.formula(paste0("nw ~ ", case$call_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  stat_val <- as.numeric(summary(f))  # vectorized possible

  # référence calculée sur la partition (en termes de tailles de groupes)
  lam <- if (is.null(case$args$lambda)) 2 else as.numeric(case$args$lambda)
  expected <- .expected_T_lambda_from_partition(part, lam)

  ok_stat <- isTRUE(all(is.finite(stat_val))) &&
            length(stat_val) == length(expected) &&
            isTRUE(all(abs(stat_val - expected) <= 1e-8))

  ok_trad <- NA
  if (exists("erpm", mode = "function")) {
    call_ergm <- erpm(f, eval.call = FALSE, verbose = FALSE)
    ok_trad <- .check_translation_contains_cliques_GW_with_args(call_ergm, args = case$args)
  }

  if (!isTRUE(quiet)) {
    cat(sprintf("  - %-14s | %-28s | stat=%s expected=%s | ok_stat=%s\n",
                case$name,
                case$call_txt,
                paste(format(stat_val, digits=8), collapse=", "),
                paste(format(expected, digits=8), collapse=", "),
                ok_stat))
  }

  list(case = case$name, signature = case$call_txt, ok_stat = ok_stat, ok_trad = ok_trad,
       stat = stat_val, expected = expected)
}

.run_phase1_summary <- function(partitions, cases, quiet = FALSE) {
  cat("\n=== PHASE 1: SUMMARY (statistique) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  all_res <- list()
  n_checks <- 0L
  n_ok <- 0L

  for (nm in names(partitions)) {
    cat("\n--- Partition ", nm, " ---\n", sep = "")
    one <- lapply(cases, function(cs) .run_one_case_summary(partitions[[nm]], cs, quiet = quiet))

    df <- data.frame(
      case      = vapply(one, `[[`, character(1), "case"),
      signature = vapply(one, `[[`, character(1), "signature"),
      ok_stat   = vapply(one, `[[`, logical(1), "ok_stat"),
      ok_trad   = vapply(one, function(x) if (is.na(x$ok_trad)) NA else isTRUE(x$ok_trad), logical(1)),
      stat      = I(lapply(one, `[[`, "stat")),
      expected  = I(lapply(one, `[[`, "expected")),
      stringsAsFactors = FALSE
    )

    .maybe_print(df, quiet = quiet)
    all_res[[nm]] <- df

    n_checks <- n_checks + sum(!is.na(df$ok_stat))
    n_ok     <- n_ok     + sum(df$ok_stat, na.rm = TRUE)

    n_checks <- n_checks + sum(!is.na(df$ok_trad))
    n_ok     <- n_ok     + sum(df$ok_trad, na.rm = TRUE)
  }

  cat(sprintf("\nBilan SUMMARY: %d / %d validations OK\n", n_ok, n_checks))
  if (n_ok < n_checks) stop(sprintf("Echec SUMMARY: %d validations KO", n_checks - n_ok))
  invisible(all_res)
}

# ======================================================================================
# PHASE 2 — ERPM FIT : vérifier qu'un fit passe et renvoie un coef fini
# ======================================================================================

.run_one_case_erpm_fit <- function(part, name, rhs_txt, quiet = FALSE) {
  if (!exists("erpm", mode = "function")) {
    if (!isTRUE(quiet)) cat(sprintf("  - %-16s | erpm() indisponible -> SKIP\n", name))
    return(list(name = name, ok = NA, coef = NA, fit = NULL))
  }

  partition <- part
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition), parent = parent.frame())

  res <- try(erpm(f, eval.loglik = TRUE, verbose = FALSE), silent = TRUE)
  if (inherits(res, "try-error")) {
    if (!isTRUE(quiet)) cat(sprintf("  - %-16s | FAIL (erpm erreur)\n", name))
    return(list(name = name, ok = FALSE, coef = NA, fit = NULL))
  }

  cf <- try(stats::coef(res), silent = TRUE)
  ok_coef  <- !(inherits(cf, "try-error")) && length(cf) > 0L && all(is.finite(cf))
  ok_class <- inherits(res, "ergm")
  ok <- isTRUE(ok_class && ok_coef)

  if (!isTRUE(quiet)) {
    cat(sprintf("  - %-16s | ok=%s | coef=%s\n",
                name, ok, if (ok_coef) paste(round(cf, 6), collapse = ", ") else "<NA>"))
  }

  list(name = name, ok = ok, coef = if (ok_coef) cf else NA, fit = res)
}

.run_phase2_fit <- function(part_ref, quiet = FALSE) {
  cat("\n=== PHASE 2: ERPM FIT (sanity check) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  fits <- list(
    FIT_default = .run_one_case_erpm_fit(part_ref, "gw_default", "cliques_GW()", quiet = quiet),
    FIT_lam2    = .run_one_case_erpm_fit(part_ref, "gw_lam2",    "cliques_GW(lambda=2)", quiet = quiet),
    FIT_lam3    = .run_one_case_erpm_fit(part_ref, "gw_lam3",    "cliques_GW(lambda=3)", quiet = quiet)
  )

  ok_fit <- vapply(fits, function(x) if (is.na(x$ok)) NA else isTRUE(x$ok), logical(1))
  n_ok  <- sum(ok_fit, na.rm = TRUE)
  n_tot <- sum(!is.na(ok_fit))

  cat(sprintf("\nBilan FIT: %d / %d OK\n", n_ok, n_tot))
  if (n_ok < n_tot) stop(sprintf("Echec FIT: %d KO", n_tot - n_ok))

  if (!isTRUE(quiet)) {
    cat("\n--- Résumés des fits OK ---\n")
    for (nm in names(fits)) {
      if (isTRUE(fits[[nm]]$ok)) {
        cat("\n#", nm, "\n")
        print(summary(fits[[nm]]$fit))
      }
    }
  }

  invisible(fits)
}

# ======================================================================================
# PHASE 3 — MCMC multi-toggle (diagnostic)
# ======================================================================================

.run_mcmc_multitoggle_probe <- function(nw) {
  # Ici le but n'est pas "faire une belle simu", c'est juste:
  #   - déclencher le code MCMC
  #   - observer les traces debug du changestat D_ (si DEBUG_CLIQUES_GW=1)
  #
  # On force un prop "sparse" (souvent multi-toggle en pratique côté ergm/ERGM),
  # et on garde verbose=TRUE pour ne pas cacher les logs.

  ctrl <- control.simulate.formula(
    MCMC.burnin   = 1000,
    MCMC.interval = 1,
    MCMC.prop     = ~ sparse
  )

  f <- nw ~ cliques_GW(lambda = c(1.5, 2, 4))

  sim <- simulate(
    f,
    nsim    = 1,
    control = ctrl,
    verbose = TRUE
  )

  print(sim)
  invisible(sim)
}

.run_phase3_mcmc <- function(part_probe) {
  cat("\n=== PHASE 3: MCMC MULTI-TOGGLE PROBE ===\n")
  cat("Objectif: vérifier que le terme cliques_GW passe en MCMC avec un changestat D_.\n")
  cat("Pour voir des traces C:\n")
  cat("  - mettre DEBUG_CLIQUES_GW=1 dans changestat_cliques_GW.c\n")
  cat("  - recompiler\n\n")

  nw <- .make_network_from_partition_via_builder(part_probe)
  .run_mcmc_multitoggle_probe(nw)

  cat("\nSi le run MCMC passe sans crash, la chaîne multi-toggle est OK.\n")
  invisible(TRUE)
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

cases <- list(
  list(name="gw_default", call_txt="cliques_GW()",                 args=list(lambda=NULL)),
  list(name="gw_lam2",    call_txt="cliques_GW(lambda=2)",         args=list(lambda=2)),
  list(name="gw_lam1_5",  call_txt="cliques_GW(lambda=1.5)",       args=list(lambda=1.5)),
  list(name="gw_lam_vec", call_txt="cliques_GW(lambda=c(1.25,4))", args=list(lambda=c(1.25,4)))
)

# ======================================================================================
# Run principal
# ======================================================================================

run_all_tests_cliques_GW <- function() {
  set.seed(42)

  cat("=== SELFTEST ERPM: cliques_GW ===\n")
  cat("R:", paste(R.version$major, R.version$minor, sep="."), "\n")
  cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")

  # PHASE 1
  summary_results <- NULL
  if (isTRUE(RUN$phase1_summary)) {
    summary_results <- .run_phase1_summary(
      partitions = partitions,
      cases      = cases,
      quiet      = isTRUE(RUN$quiet_phase1)
    )
  } else {
    cat("\n=== PHASE 1: SUMMARY ===\n")
    cat("SKIP (désactivée via RUN$phase1_summary = FALSE)\n")
  }

  # PHASE 2
  fit_results <- NULL
  if (isTRUE(RUN$phase2_fit)) {
    fit_results <- .run_phase2_fit(
      part_ref = partitions$P1,
      quiet    = isTRUE(RUN$quiet_phase2)
    )
  } else {
    cat("\n=== PHASE 2: ERPM FIT ===\n")
    cat("SKIP (désactivée via RUN$phase2_fit = FALSE)\n")
  }

  # PHASE 3
  if (isTRUE(RUN$phase3_mcmc)) {
    .run_phase3_mcmc(part_probe = partitions$P1)
  } else {
    cat("\n=== PHASE 3: MCMC MULTI-TOGGLE PROBE ===\n")
    cat("SKIP (désactivée via RUN$phase3_mcmc = FALSE)\n")
  }

  invisible(list(summary_results = summary_results, fit_results = fit_results))
}

# Exécution quand lancé en script
if (identical(environment(), globalenv())) {
  run_all_tests_cliques_GW()
}

# --------------------------------------------------------------------------------------
# Fin de script: on désactive le patch si présent (histoire de ne pas laisser traîner)
# --------------------------------------------------------------------------------------
if (exists("ergm_patch_disable", mode = "function")) {
  ergm_patch_disable()
}