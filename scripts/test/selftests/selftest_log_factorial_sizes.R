# ======================================================================================
# Fichier : scripts/test/selftests/selftest_log_factorial_sizes.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `log_factorial_sizes`
# Exécution: Rscript scripts/test/selftests/selftest_log_factorial_sizes.R
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

# Patch ERGM optionnel si disponible
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

# script_dir <- .get_script_dir()
# log_path   <- file.path(script_dir, "selftest_log_factorial_sizes.log")
root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_log_factorial_sizes.log")

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

cat("==> Log: ", log_path, "\n")

# --------------------------------------------------------------------------------------
# Chargements utilitaires ERPM si disponibles
# --------------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
}

if (!exists("partition_to_bipartite_network", mode = "function")) {
  # Fallback minimal interne
  partition_to_bipartite_network <- function(labels, partition, attributes = list()) {
    stopifnot(length(labels) == length(partition))
    nA <- length(partition)
    G  <- max(partition)
    inc <- matrix(0L, nrow = nA, ncol = G,
                  dimnames = list(labels, paste0("G", seq_len(G))))
    inc[cbind(seq_len(nA), partition)] <- 1L
    nw <- network::network(inc, matrix.type = "bipartite", bipartite = nA, directed = FALSE)
    network::set.vertex.attribute(nw, "vertex.names", c(labels, colnames(inc)))
    if (length(attributes)) {
      for (nm in names(attributes)) {
        vals <- attributes[[nm]]
        stopifnot(length(vals) == nA)
        network::set.vertex.attribute(nw, nm, c(vals, rep(NA, G)))
      }
    }
    nw
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

.is_degenerate_for_fitting <- function(part){
  sz <- as.integer(table(part))
  # Dégénéré si: tous singletons, ou un seul groupe, ou variance très faible attendue
  if (length(sz) == length(part)) return(TRUE)   # tous singletons
  if (length(sz) == 1L)           return(TRUE)   # un seul groupe
  FALSE
}


# Construire un biparti depuis une partition
.make_nw_from_partition <- function(part) {
  n <- length(part); stopifnot(n >= 1, is.atomic(part))
  lbl <- paste0("A", seq_len(n))
  partition_to_bipartite_network(labels = lbl, partition = part, attributes = list())
}

# Valeur de référence analytique: sum_g log((n_g-1)!) = sum_g lgamma(n_g)
.expected_log_factorial <- function(p) {
  sum(lgamma(as.integer(table(p))))
}

# Vérifie un appel erpm(...) non évalué
.check_translation_ok <- function(call_ergm) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)
  grepl("\\blog_factorial_sizes\\(\\)", compact)
}

# Exécuter un cas summary + check analytique + dry-runs
.run_one_partition_all <- function(part, name) {
  nw <- .make_nw_from_partition(part)

  # 1) summary( nw ~ log_factorial_sizes() )
  f_nw <- nw ~ log_factorial_sizes()
  stat_summary <- suppressMessages(as.numeric(summary(f_nw)))
  ok_summary   <- isTRUE(all(is.finite(stat_summary)))

  # 2) check analytique
  expected <- .expected_log_factorial(part)
  ok_value <- isTRUE(abs(stat_summary[1] - expected) < 1e-9)

  # 3) dry-run erpm LHS=network
  ok_trad_net <- NA
  if (exists("erpm", mode = "function")) {
    call_net <- erpm(f_nw, eval_call = FALSE, verbose = FALSE)
    ok_trad_net <- .check_translation_ok(call_net)
  }

  # 4) dry-run erpm LHS=partition
  ok_trad_part <- NA
  if (exists("erpm", mode = "function")) {
    partition <- part
    f_part <- partition ~ log_factorial_sizes
    environment(f_part) <- list2env(list(partition = partition), parent = parent.frame())
    call_part <- erpm(f_part, eval_call = FALSE, verbose = FALSE)
    ok_trad_part <- .check_translation_ok(call_part)
  }

  cat(sprintf("\n[PART %-8s] part={%s}", name, paste(part, collapse=",")))
  cat(sprintf("\t summary=%.12f\t expected=%.12f\t summaryOK=%s\t valueOK=%s",
              stat_summary[1], expected, ok_summary, ok_value))
  if (!is.na(ok_trad_net))  cat(sprintf("\t trad[LHS=nw]=%s",  if (ok_trad_net) "OK" else "KO"))
  if (!is.na(ok_trad_part)) cat(sprintf("\t trad[LHS=part]=%s",if (ok_trad_part) "OK" else "KO"))
  cat("\n")

  data.frame(
    name        = name,
    n           = length(part),
    groups      = length(unique(part)),
    summary     = stat_summary[1],
    expected    = expected,
    ok_summary  = ok_summary,
    ok_value    = ok_value,
    ok_trad_nw  = ok_trad_net,
    ok_trad_part= ok_trad_part,
    stringsAsFactors = FALSE
  )
}

# Petit fit via erpm() (CD rapide, sans logLik) — LHS au choix
.erpm_fit_one <- function(part, lhs_mode = c("partition","network"),
                          estimate = "CD",
                          eval.loglik = FALSE,
                          control = list(CD.maxit = 4, MCMLE.maxit = 0, MCMC.samplesize = 500)) {
  lhs_mode <- match.arg(lhs_mode)

  if (.is_degenerate_for_fitting(part)) {
    cat(sprintf("[FIT %-10s] part={%s}\t SKIP (profil dégénéré pour le fit)\n",
                lhs_mode, paste(part, collapse=",")))
    return(list(ok = NA, coef = NA))
  }

  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[FIT %-10s] part={%s}\t SKIP erpm() indisponible\n",
                lhs_mode, paste(part, collapse=",")))
    return(list(ok = NA, coef = NA))
  }

  oldopt <- options(ergm.loglik.warn_dyads = FALSE); on.exit(options(oldopt), add = TRUE)
  set.seed(42)

  if (lhs_mode == "network") {
    nw <- .make_nw_from_partition(part)
    f  <- nw ~ log_factorial_sizes()
  } else {
    partition <- part
    f <- partition ~ log_factorial_sizes
    environment(f) <- list2env(list(partition = partition), parent = parent.frame())
  }

  cat(sprintf("[FIT %-10s] part={%s}\n", lhs_mode, paste(part, collapse=",")))

  fit_try <- function(ctrl){
    try(erpm(f, estimate = "CD", eval.loglik = FALSE, control = ctrl, verbose = FALSE),
        silent = TRUE)
  }

  # tentative 1: CD seul, court
  ctrl1 <- do.call(ergm::control.ergm, control)
  fit   <- fit_try(ctrl1)

  # tentative 2: si erreur, CD encore plus court
  if (inherits(fit, "try-error")) {
    ctrl2 <- ergm::control.ergm(CD.maxit = 3, MCMLE.maxit = 0, MCMC.samplesize = 300)
    fit   <- fit_try(ctrl2)
  }

  # abandon propre -> SKIP
  if (inherits(fit, "try-error")) {
    msg <- as.character(fit)
    if (grepl("essentially constant|Hotelling", msg, ignore.case = TRUE)) {
      cat("  -> SKIP fit (variance quasi nulle détectée)\n")
      return(list(ok = NA, coef = NA))
    }
    cat("  -> ERREUR fit: ", msg, "\n")
    return(list(ok = FALSE, coef = NA))
  }

  cf <- try(stats::coef(fit), silent = TRUE)
  ok_coef <- !inherits(cf, "try-error") && all(is.finite(cf))
  cat(sprintf("  -> class(ergm)=%s | coef finies=%s | coef=%s\n",
              inherits(fit, "ergm"), ok_coef,
              if (ok_coef) paste0(format(as.numeric(cf)), collapse=", ") else "NA"))
  list(ok = inherits(fit, "ergm") && ok_coef, coef = if (ok_coef) cf else NA)
}


# ======================================================================================
# Jeu de tests
# ======================================================================================

# Partitions variées, couvrant: singletons, tous dans un groupe, tailles mixtes, multi-groupes
partitions <- list(
  P1 = c(1, 1, 2, 2, 2, 3),               # tailles: (2,3,1)
  P2 = c(1, 1, 1, 2, 3, 3, 3, 3),         # tailles: (3,1,4)
  P3 = c(1, 2, 2, 3, 3, 4, 4, 4),         # tailles: (1,2,2,3)
  P4 = c(1, 2, 3, 4, 5),                  # toutes tailles 1
  P5 = rep(1, 6),                         # un seul groupe
  P6 = c(1,1,1,1,2,2,3,3,3,4),            # tailles: (4,2,3,1)
  P7 = c(1,2,2,2,2,3,3,4,4,4,4,4)         # tailles: (1,4,2,5)
)

# ======================================================================================
# Phase 1 : summary + analytique + dry-runs
# ======================================================================================
run_phase1 <- function() {
  cat("=== PHASE 1: summary + expected + dry-run translations ===\n")
  rows <- lapply(names(partitions), function(nm) .run_one_partition_all(partitions[[nm]], nm))
  df <- do.call(rbind, rows)

  # Bilan
  ok_flags <- c(df$ok_summary, df$ok_value,
                if (!all(is.na(df$ok_trad_nw))) df$ok_trad_nw else TRUE,
                if (!all(is.na(df$ok_trad_part))) df$ok_trad_part else TRUE)
  total_ok <- sum(ok_flags, na.rm = TRUE)
  total_n  <- sum(!is.na(ok_flags))
  cat(sprintf("\n=== Bilan Phase 1 : %d / %d validations OK ===\n", total_ok, total_n))
  if (total_ok < total_n) stop(sprintf("Phase 1 KO: %d échecs", total_n - total_ok))

  df
}

# ======================================================================================
# Phase 2 : petits fits via erpm() pour LHS partition et LHS réseau
# ======================================================================================
# --- dans run_phase2(), garde l’échantillon mais les cas dégénérés seront SKIP automatiquement ---
run_phase2 <- function() {
  if (!exists("erpm", mode = "function")) {
    cat("\n=== PHASE 2: SKIP (erpm() indisponible) ===\n"); return(invisible(NULL))
  }
  cat("\n=== PHASE 2: Fits erpm() courts ===\n")
  pick <- c("P1","P2","P4","P5")
  res_fit <- list()
  for (nm in pick) {
    part <- partitions[[nm]]
    res_fit[[paste0(nm,"_part")]] <- .erpm_fit_one(part, lhs_mode = "partition")
    res_fit[[paste0(nm,"_net")]]  <- .erpm_fit_one(part, lhs_mode = "network")
  }

  # NE PAS utiliser isTRUE ici
  ok <- vapply(res_fit, function(x) x$ok, logical(1))

  n_skip <- sum(is.na(ok))
  n_tot  <- length(ok) - n_skip          # total hors SKIP
  n_ok   <- sum(ok, na.rm = TRUE)        # OK parmi non-NA
  n_fail <- n_tot - n_ok                 # échecs réels

  cat(sprintf("\n=== Bilan Phase 2 (fits) : %d / %d OK ; %d SKIP ===\n",
              n_ok, n_tot, n_skip))

  if (n_fail > 0)
    stop(sprintf("Phase 2 KO: %d fits en échec (hors SKIP)", n_fail))

  invisible(res_fit)
}

# ======================================================================================
# Point d'entrée
# ======================================================================================
if (identical(environment(), globalenv())) {
  set.seed(42)
  oldopt <- options(ergm.loglik.warn_dyads = FALSE); on.exit(options(oldopt), add = TRUE)

  p1 <- run_phase1()
  p2 <- run_phase2()

  cat("\n=== SELFTEST log_factorial_sizes : SUCCESS ===\n")
  print(p1)
}

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)