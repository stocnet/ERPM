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
options(ergm.loglik.warn_dyads = FALSE)

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
# Chargements utilitaires ERPM / wrapper
# --------------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
}

if (!exists("erpm", mode = "function") || !exists("build_bipartite_from_inputs", mode = "function")) {
  if (file.exists("R/erpm_wrapper.R")) {
    source("R/erpm_wrapper.R", local = FALSE)
  } else {
    cat("[WARN] erpm()/build_bipartite_from_inputs indisponibles. Fallbacks utilisés.\n")
  }
}

# Fallback minimal si wrapper absent
if (!exists("partition_to_bipartite_network", mode = "function")) {
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

# --------------------------------------------------------------------------------------
# Construction bipartie via wrapper ERPM
# --------------------------------------------------------------------------------------
.erpm_build_bipartite_network <- function(partition, nodes_df = NULL) {
  n <- length(partition)
  if (is.null(nodes_df)) {
    nodes_df <- data.frame(label = paste0("A", seq_len(n)), stringsAsFactors = FALSE)
  } else {
    stopifnot(nrow(nodes_df) == n)
    if (!"label" %in% names(nodes_df)) {
      nodes_df$label <- paste0("A", seq_len(n))
    }
  }
  attrs <- as.list(nodes_df[, setdiff(names(nodes_df), "label"), drop = FALSE])

  # 1) Tentatives via build_bipartite_from_inputs du wrapper
  if (exists("build_bipartite_from_inputs", mode = "function")) {
    builder <- get("build_bipartite_from_inputs")
    # Signature 1
    out <- try(builder(partition = partition, nodes = nodes_df), silent = TRUE)
    # Signature 2
    if (inherits(out, "try-error") || is.null(out)) {
      out <- try(builder(partition = partition, labels = nodes_df$label, attributes = attrs), silent = TRUE)
    }
    # Extraction d'un objet 'network'
    if (!inherits(out, "try-error") && !is.null(out)) {
      if (inherits(out, "network")) return(out)
      if (is.list(out)) {
        for (nm in c("network","nw","g","graph","bip","net")) {
          if (!is.null(out[[nm]]) && inherits(out[[nm]], "network")) return(out[[nm]])
        }
      }
    }
  }

  # 2) Fallback explicite
  partition_to_bipartite_network(labels = nodes_df$label, partition = partition, attributes = attrs)
}

# --------------------------------------------------------------------------------------
# Fonctions locales explicites (SUMMARY / ERPM)
# --------------------------------------------------------------------------------------

# Référence analytique: sum_g log((n_g-1)!) = sum_g lgamma(n_g)
summary_expected_log_factorial_value <- function(partition) {
  sum(lgamma(as.integer(table(partition))))
}

# Vérifie qu'un appel erpm(...) contient bien log_factorial_sizes()
erpm_check_translation_contains_lfs <- function(call_ergm) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)
  grepl("\\blog_factorial_sizes\\(\\)", compact)
}

# Exécute summary côté réseau biparti (construit via wrapper)
summary_on_bipartite_from_wrapper <- function(partition) {
  nw <- .erpm_build_bipartite_network(partition)
  as.numeric(suppressMessages(summary(nw ~ log_factorial_sizes())))
}

# Secoue la chaîne ERPM: dry-runs traduction pour LHS=network et LHS=partition
erpm_dryruns_check_translations <- function(partition) {
  if (!exists("erpm", mode = "function")) return(list(ok_net = NA, ok_part = NA))

  # LHS = network
  nw <- .erpm_build_bipartite_network(partition)
  call_net  <- erpm(nw ~ log_factorial_sizes(), eval_call = FALSE, verbose = FALSE)
  ok_net    <- erpm_check_translation_contains_lfs(call_net)

  # LHS = partition
  f_part <- partition ~ log_factorial_sizes
  environment(f_part) <- list2env(list(partition = partition), parent = parent.frame())
  call_part <- erpm(f_part, eval_call = FALSE, verbose = FALSE)
  ok_part   <- erpm_check_translation_contains_lfs(call_part)

  list(ok_net = ok_net, ok_part = ok_part)
}

# Un cas complet: summary + analytique + dry-runs
summary_run_one_partition_all <- function(partition, name) {
  stat_summary <- summary_on_bipartite_from_wrapper(partition)
  expected     <- summary_expected_log_factorial_value(partition)
  ok_summary   <- isTRUE(all(is.finite(stat_summary)))
  ok_value     <- isTRUE(abs(stat_summary[1] - expected) < 1e-9)

  tr <- erpm_dryruns_check_translations(partition)

  cat(sprintf("\n[SUMMARY %-8s] part={%s}\t summary=%.12f\t expected=%.12f\t summaryOK=%s\t valueOK=%s\t trad[LHS=nw]=%s\t trad[LHS=part]=%s\n",
              name, paste(partition, collapse=","), stat_summary[1], expected,
              ok_summary, ok_value,
              if (is.na(tr$ok_net)) "NA" else if (tr$ok_net) "OK" else "KO",
              if (is.na(tr$ok_part)) "NA" else if (tr$ok_part) "OK" else "KO"))

  data.frame(
    name          = name,
    n             = length(partition),
    groups        = length(unique(partition)),
    summary       = stat_summary[1],
    expected      = expected,
    ok_summary    = ok_summary,
    ok_value      = ok_value,
    ok_trad_nw    = tr$ok_net,
    ok_trad_part  = tr$ok_part,
    stringsAsFactors = FALSE
  )
}

# Fit via erpm() avec MLE + logLik, LHS = partition ou réseau
erpm_run_fit_one <- function(partition, lhs_mode = c("partition","network"),
                             estimate = "MLE", eval.loglik = TRUE,
                             control = control.ergm(MCMLE.maxit = 10, MCMC.samplesize = 3000)) {
  lhs_mode <- match.arg(lhs_mode)
  # profils dégénérés → SKIP
  sz <- as.integer(table(partition))
  if (length(sz) == length(partition) || length(sz) == 1L) {
    cat(sprintf("[ERPM-FIT %-10s] part={%s}\t SKIP (profil dégénéré)\n",
                lhs_mode, paste(partition, collapse=",")))
    return(list(ok = NA, error = FALSE, coef = NA, fit = NULL))
  }
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT %-10s] part={%s}\t SKIP erpm() indisponible\n",
                lhs_mode, paste(partition, collapse=",")))
    return(list(ok = NA, error = TRUE, coef = NA, fit = NULL))
  }

  set.seed(42)
  if (lhs_mode == "network") {
    nw <- .erpm_build_bipartite_network(partition)
    f  <- nw ~ log_factorial_sizes()
  } else {
    # LHS = partition, lier le bon symbole dans l'environnement
    f <- partition ~ log_factorial_sizes
    environment(f) <- list2env(list(partition = partition), parent = parent.frame())
  }

  cat(sprintf("[ERPM-FIT %-10s] part={%s}\n", lhs_mode, paste(partition, collapse=",")))

  res <- withCallingHandlers(
    try(erpm(f, estimate = estimate, eval.loglik = eval.loglik,
             control = control, verbose = FALSE), silent = TRUE),
    warning = function(w) { message <- conditionMessage(w); cat("  - WARNING: ", message, "\n", sep=""); invokeRestart("muffleWarning") }
  )
  if (inherits(res, "try-error")) {
    cat("  -> ERREUR fit: ", as.character(res), "\n", sep = "")
    return(list(ok = FALSE, error = TRUE, coef = NA, fit = NULL))
  }

  cf <- try(stats::coef(res), silent = TRUE)
  ok_coef  <- !(inherits(cf, "try-error")) && all(is.finite(cf))
  ok_class <- inherits(res, "ergm")
  cat(sprintf("  -> class(ergm)=%s | coef finies=%s | coef=%s\n",
              ok_class, if (ok_coef) "OK" else "KO",
              if (ok_coef) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

  list(ok = ok_class && ok_coef, error = FALSE, coef = if (ok_coef) cf else NA, fit = res)
}

# --------------------------------------------------------------------------------------
# Jeu de tests
# --------------------------------------------------------------------------------------
partitions <- list(
  P1 = c(1, 1, 2, 2, 2, 3),               # tailles: (2,3,1)
  P2 = c(1, 1, 1, 2, 3, 3, 3, 3),         # tailles: (3,1,4)
  P3 = c(1, 2, 2, 3, 3, 4, 4, 4),         # tailles: (1,2,2,3)
  P4 = c(1, 2, 3, 4, 5),                  # toutes tailles 1
  P5 = rep(1, 6),                         # un seul groupe
  P6 = c(1,1,1,1,2,2,3,3,3,4),            # tailles: (4,2,3,1)
  P7 = c(1,2,2,2,2,3,3,4,4,4,4,4)         # tailles: (1,4,2,5)
)

# --------------------------------------------------------------------------------------
# Phase 1 : summary + analytique + dry-runs
# --------------------------------------------------------------------------------------
summary_run_phase1_all <- function() {
  cat("=== PHASE 1: summary + expected + ERPM dry-runs ===\n")
  rows <- lapply(names(partitions), function(nm) summary_run_one_partition_all(partitions[[nm]], nm))
  df <- do.call(rbind, rows)

  # Bilan validations
  ok_flags <- c(df$ok_summary, df$ok_value,
                if (!all(is.na(df$ok_trad_nw))) df$ok_trad_nw else TRUE,
                if (!all(is.na(df$ok_trad_part))) df$ok_trad_part else TRUE)
  total_ok <- sum(ok_flags, na.rm = TRUE)
  total_n  <- sum(!is.na(ok_flags))
  cat(sprintf("\n=== Bilan Phase 1 : %d / %d validations OK ===\n", total_ok, total_n))
  if (total_ok < total_n) stop(sprintf("Phase 1 KO: %d échecs", total_n - total_ok))

  df
}

# --------------------------------------------------------------------------------------
# Phase 2 : fits MLE + logLik, pour LHS partition et LHS réseau
# --------------------------------------------------------------------------------------
erpm_run_phase2_fits <- function() {
  if (!exists("erpm", mode = "function")) {
    cat("\n=== PHASE 2: SKIP (erpm() indisponible) ===\n"); return(invisible(NULL))
  }
  cat("\n=== PHASE 2: Fits erpm() MLE + logLik ===\n")
  pick <- c("P1","P2","P4","P5")
  fit_results <- list()
  for (nm in pick) {
    part <- partitions[[nm]]
    fit_results[[paste0(nm,"_part")]] <- erpm_run_fit_one(part, lhs_mode = "partition")
    fit_results[[paste0(nm,"_net")]]  <- erpm_run_fit_one(part, lhs_mode = "network")
  }

    # Bilan des fits  — gère explicitement SKIP (= NA)
    ok_raw <- vapply(fit_results, function(x) x$ok, logical(1))
    n_skip <- sum(is.na(ok_raw))
    n_tot  <- length(ok_raw) - n_skip            # seulement les cas évalués
    n_ok   <- sum(ok_raw, na.rm = TRUE)          # OK parmi non-SKIP
    n_fail <- n_tot - n_ok                       # vrais échecs

    cat(sprintf("\n=== Bilan fits erpm() : %d / %d OK ; %d SKIP ===\n",
                n_ok, n_tot, n_skip))

    # Résumés détaillés des fits ERPM réussis
    cat("\n=== Résumés détaillés des fits ERPM réussis ===\n")
    for (nm in names(fit_results)) {
      fit_obj <- fit_results[[nm]]
      if (isTRUE(fit_obj$ok) && inherits(fit_obj$coef, "numeric")) {
        cat(sprintf("\n--- Résumé fit %s ---\n", nm))
        print(summary(fit_obj$fit))
      }
    }

    if (n_fail > 0) stop(sprintf("Echec fits: %d KO", n_fail))
    invisible(fit_results)
}

# --------------------------------------------------------------------------------------
# Point d'entrée
# --------------------------------------------------------------------------------------
if (identical(environment(), globalenv())) {
  set.seed(42)

  summary_results <- summary_run_phase1_all()
  fit_results     <- erpm_run_phase2_fits()

  cat("\n=== SELFTEST log_factorial_sizes : SUCCESS ===\n")
  print(summary_results)

  invisible(list(summary_results = summary_results, fit_results = fit_results))
}

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)
