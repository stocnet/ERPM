# ======================================================================================
# File    : scripts/test/selftests/selftest_log_factorial_sizes.R
# Auteur : Jérémie Chichignoud - Cub'itech
# Object  : Self-test autonome pour l'effet ERPM/ERGM `log_factorial_sizes`
#
# But du fichier
#   - PHASE 1 (SUMMARY) : valider la statistique via summary(nw ~ log_factorial_sizes(.)).
#   - PHASE 2 (ERPM FIT): valider que erpm() construit un modèle et renvoie des coefs finis.
#   - PHASE 3 (MCMC)    : diagnostic "multi-toggle" (déclencher D_CHANGESTAT_FN)
#                         et observer les traces debug côté C si activées.
#
# Important
#   - Les phases 1/2 peuvent spammer la console (print de dataframes + summaries).
#   - Pour bosser proprement sur la phase 3, on peut désactiver 1/2 via des flags.
#   - Le terme `log_factorial_sizes` n'a pas d'arguments.
#   - Le point critique: le changestat est en D_ (multi-toggle), donc InitErgmTerm
#     doit retourner d_func=TRUE, et ce test doit pouvoir déclencher un ntoggles>1.
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
  if (exists("ergm_patch_enable", mode = "function")) {
    ergm_patch_enable()
  }
}

# --------------------------------------------------------------------------------------
# Chargement package/terme + wrapper ERPM
# --------------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
}

if (!exists("InitErgmTerm.log_factorial_sizes", mode = "function")) {
  stop("InitErgmTerm.log_factorial_sizes introuvable après load_all().")
}

# Pour certains setups, le wrapper n'est pas attaché automatiquement.
if (!exists("erpm", mode = "function") || !exists("build_bipartite_from_inputs", mode = "function")) {
  if (file.exists("R/erpm.R")) {
    source("R/erpm.R", local = FALSE)
  } else {
    message("[WARN] R/erpm.R introuvable, certains checks ERPM seront SKIP.")
  }
}

# --------------------------------------------------------------------------------------
# Réglages de run (isoler facilement PHASE 3)
# --------------------------------------------------------------------------------------
RUN <- list(
  phase1_summary = TRUE,
  phase2_fit     = TRUE,
  phase3_mcmc    = TRUE,

  quiet_phase1   = FALSE,
  quiet_phase2   = FALSE
)

# --------------------------------------------------------------------------------------
# Fallback minimal si wrapper absent (conversion partition -> réseau biparti)
# --------------------------------------------------------------------------------------
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
# Construction bipartie via wrapper ERPM si présent
# --------------------------------------------------------------------------------------
.erpm_build_bipartite_network <- function(partition, nodes_df = NULL) {
  n <- length(partition)
  if (is.null(nodes_df)) {
    nodes_df <- data.frame(label = paste0("A", seq_len(n)), stringsAsFactors = FALSE)
  } else {
    stopifnot(nrow(nodes_df) == n)
    if (!"label" %in% names(nodes_df)) nodes_df$label <- paste0("A", seq_len(n))
  }
  attrs <- as.list(nodes_df[, setdiff(names(nodes_df), "label"), drop = FALSE])

  # 1) Tentatives via build_bipartite_from_inputs du wrapper
  if (exists("build_bipartite_from_inputs", mode = "function")) {
    builder <- get("build_bipartite_from_inputs")
    out <- try(builder(partition = partition, nodes = nodes_df), silent = TRUE)
    if (inherits(out, "try-error") || is.null(out)) {
      out <- try(builder(partition = partition, labels = nodes_df$label, attributes = attrs), silent = TRUE)
    }
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
# Référence analytique
# --------------------------------------------------------------------------------------
# Stat = sum_g lgamma(n_g) ; convention: lgamma(0)=0 (mais table(partition) exclut les vides)
summary_expected_log_factorial_value <- function(partition) {
  sum(lgamma(as.integer(table(partition))))
}

# --------------------------------------------------------------------------------------
# Vérification: traduction erpm(...) contient bien log_factorial_sizes (avec ou sans ())
# --------------------------------------------------------------------------------------
erpm_check_translation_contains_lfs <- function(call_ergm) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)
  # accepte log_factorial_sizes et log_factorial_sizes()
  grepl("\\blog_factorial_sizes(\\(\\))?", compact)
}

# --------------------------------------------------------------------------------------
# SUMMARY sur réseau biparti
# --------------------------------------------------------------------------------------
summary_on_bipartite_from_wrapper <- function(partition) {
  nw <- .erpm_build_bipartite_network(partition)
  # accepte la forme avec parenthèses (arglist vide), ou sans.
  as.numeric(suppressMessages(summary(nw ~ log_factorial_sizes())))
}

# --------------------------------------------------------------------------------------
# Dry-runs traduction ERPM: LHS=network et LHS=partition
# --------------------------------------------------------------------------------------
erpm_dryruns_check_translations <- function(partition) {
  if (!exists("erpm", mode = "function")) return(list(ok_net = NA, ok_part = NA))

  nw <- .erpm_build_bipartite_network(partition)
  call_net  <- erpm(nw ~ log_factorial_sizes(), eval.call = FALSE, verbose = FALSE)
  ok_net    <- erpm_check_translation_contains_lfs(call_net)

  f_part <- partition ~ log_factorial_sizes
  environment(f_part) <- list2env(list(partition = partition), parent = parent.frame())
  call_part <- erpm(f_part, eval.call = FALSE, verbose = FALSE)
  ok_part   <- erpm_check_translation_contains_lfs(call_part)

  list(ok_net = ok_net, ok_part = ok_part)
}

# --------------------------------------------------------------------------------------
# Un cas complet: summary + analytique + dry-runs
# --------------------------------------------------------------------------------------
summary_run_one_partition_all <- function(partition, name, quiet = FALSE) {
  stat_summary <- summary_on_bipartite_from_wrapper(partition)
  expected     <- summary_expected_log_factorial_value(partition)

  ok_summary <- isTRUE(all(is.finite(stat_summary)))
  ok_value   <- isTRUE(abs(stat_summary[1] - expected) < 1e-9)

  tr <- erpm_dryruns_check_translations(partition)

  if (!isTRUE(quiet)) {
    cat(sprintf(
      "\n[SUMMARY %-8s] part={%s}\t summary=%.12f\t expected=%.12f\t summaryOK=%s\t valueOK=%s\t trad[LHS=nw]=%s\t trad[LHS=part]=%s\n",
      name, paste(partition, collapse=","), stat_summary[1], expected,
      ok_summary, ok_value,
      if (is.na(tr$ok_net))  "NA" else if (tr$ok_net)  "OK" else "KO",
      if (is.na(tr$ok_part)) "NA" else if (tr$ok_part) "OK" else "KO"
    ))
  }

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

summary_run_phase1_all <- function(partitions, quiet = FALSE) {
  cat("\n=== PHASE 1: SUMMARY (statistique) + expected + ERPM dry-runs ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  rows <- lapply(names(partitions), function(nm) summary_run_one_partition_all(partitions[[nm]], nm, quiet = quiet))
  df <- do.call(rbind, rows)

  # Bilan validations (ok_trad_* optionnels)
  ok_flags <- c(df$ok_summary, df$ok_value,
                if (!all(is.na(df$ok_trad_nw)))   df$ok_trad_nw   else TRUE,
                if (!all(is.na(df$ok_trad_part))) df$ok_trad_part else TRUE)

  total_ok <- sum(ok_flags, na.rm = TRUE)
  total_n  <- sum(!is.na(ok_flags))
  cat(sprintf("\n=== Bilan Phase 1 : %d / %d validations OK ===\n", total_ok, total_n))
  if (total_ok < total_n) stop(sprintf("Phase 1 KO: %d échecs", total_n - total_ok))

  if (!isTRUE(quiet)) print(df)
  invisible(df)
}

# --------------------------------------------------------------------------------------
# PHASE 2 — Fit via erpm() (sanity)
# --------------------------------------------------------------------------------------
erpm_run_fit_one <- function(partition, lhs_mode = c("partition","network"),
                             eval.loglik = TRUE, quiet = FALSE) {
  lhs_mode <- match.arg(lhs_mode)

  # Profils dégénérés → SKIP (un seul groupe ou tous singletons)
  sz <- as.integer(table(partition))
  if (length(sz) == length(partition) || length(sz) == 1L) {
    if (!isTRUE(quiet)) cat(sprintf("[ERPM-FIT %-10s] part={%s}\t SKIP (profil dégénéré)\n",
                                    lhs_mode, paste(partition, collapse=",")))
    return(list(ok = NA, error = FALSE, coef = NA, fit = NULL))
  }

  if (!exists("erpm", mode = "function")) {
    if (!isTRUE(quiet)) cat(sprintf("[ERPM-FIT %-10s] part={%s}\t SKIP erpm() indisponible\n",
                                    lhs_mode, paste(partition, collapse=",")))
    return(list(ok = NA, error = TRUE, coef = NA, fit = NULL))
  }

  set.seed(42)
  if (lhs_mode == "network") {
    nw <- .erpm_build_bipartite_network(partition)
    f  <- nw ~ log_factorial_sizes()
  } else {
    f <- partition ~ log_factorial_sizes
    environment(f) <- list2env(list(partition = partition), parent = parent.frame())
  }

  if (!isTRUE(quiet)) cat(sprintf("[ERPM-FIT %-10s] part={%s}\n", lhs_mode, paste(partition, collapse=",")))

  res <- withCallingHandlers(
    try(
      erpm(f, eval.loglik = eval.loglik, verbose = FALSE),
      silent = TRUE
    ),
    warning = function(w) {
      if (!isTRUE(quiet)) cat("  - WARNING: ", conditionMessage(w), "\n", sep = "")
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(res, "try-error")) {
    if (!isTRUE(quiet)) cat("  -> ERREUR fit: ", as.character(res), "\n", sep = "")
    return(list(ok = FALSE, error = TRUE, coef = NA, fit = NULL))
  }

  cf <- try(stats::coef(res), silent = TRUE)
  ok_coef  <- !(inherits(cf, "try-error")) && length(cf) > 0L && all(is.finite(cf))
  ok_class <- inherits(res, "ergm")

  if (!isTRUE(quiet)) {
    cat(sprintf("  -> class(ergm)=%s | coef finies=%s | coef=%s\n",
                ok_class, if (ok_coef) "OK" else "KO",
                if (ok_coef) paste(format(as.numeric(cf)), collapse=", ") else "NA"))
  }

  list(ok = isTRUE(ok_class && ok_coef), error = FALSE, coef = if (ok_coef) cf else NA, fit = res)
}

erpm_run_phase2_fits <- function(partitions, quiet = FALSE) {
  if (!exists("erpm", mode = "function")) {
    cat("\n=== PHASE 2: SKIP (erpm() indisponible) ===\n")
    return(invisible(NULL))
  }
  cat("\n=== PHASE 2: Fits erpm() MLE + logLik ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  # Ne garder que des partitions non dégénérées (réduire les SKIP)
  pick <- intersect(names(partitions), c("P1","P2"))

  fit_results <- list()
  for (nm in pick) {
    part <- partitions[[nm]]
    fit_results[[paste0(nm,"_part")]] <- erpm_run_fit_one(part, lhs_mode = "partition", quiet = quiet)
    fit_results[[paste0(nm,"_net")]]  <- erpm_run_fit_one(part, lhs_mode = "network",   quiet = quiet)
  }

  ok_raw <- vapply(fit_results, function(x) x$ok, logical(1))
  n_skip <- sum(is.na(ok_raw))
  n_tot  <- length(ok_raw) - n_skip
  n_ok   <- sum(ok_raw, na.rm = TRUE)
  n_fail <- n_tot - n_ok

  cat(sprintf("\n=== Bilan fits erpm() : %d / %d OK ; %d SKIP ===\n", n_ok, n_tot, n_skip))

  if (!isTRUE(quiet)) {
    cat("\n=== Résumés détaillés des fits ERPM réussis ===\n")
    for (nm in names(fit_results)) {
      fit_obj <- fit_results[[nm]]
      if (isTRUE(fit_obj$ok) && inherits(fit_obj$fit, "ergm")) {
        cat(sprintf("\n--- Résumé fit %s ---\n", nm))
        print(summary(fit_obj$fit))
      }
    }
  }

  if (n_fail > 0) stop(sprintf("Echec fits: %d KO", n_fail))
  invisible(fit_results)
}

# --------------------------------------------------------------------------------------
# PHASE 3 — MCMC multi-toggle probe
# --------------------------------------------------------------------------------------
.run_mcmc_multitoggle_probe <- function(nw) {
  # Objectif: déclencher le chemin MCMC et provoquer des propositions multi-toggle.
  # Si DEBUG_LOG_FACTORIAL_SIZES=1 côté C, tu dois voir passer:
  #   "[log_factorial_sizes] MULTI-TOGGLE ntoggles=..."
  #
  # Important: on garde verbose=TRUE exprès.

  ctrl <- control.simulate.formula(
    MCMC.burnin   = 1000,
    MCMC.interval = 1,
    # Le kernel "sparse" est un bon candidat pour produire des listes de toggles.
    # (selon versions ergm, le contenu exact dépend des propositions activées).
    MCMC.prop     = ~ sparse
  )

  f <- nw ~ log_factorial_sizes

  sim <- simulate(
    f,
    nsim    = 1,
    control = ctrl,
    verbose = TRUE
  )

  print(sim)
  invisible(sim)
}

run_phase3_mcmc <- function(partition_probe) {
  cat("\n=== PHASE 3: MCMC MULTI-TOGGLE PROBE ===\n")
  cat("Objectif: déclencher le changestat D_ (multi-toggle) et observer ntoggles>1.\n")
  cat("Si tu ne vois rien, c'est souvent:\n")
  cat("  - DEBUG_LOG_FACTORIAL_SIZES=0 côté C (normal), ou\n")
  cat("  - le kernel n'a pas produit de multi-toggle sur ce run.\n\n")

  nw <- .erpm_build_bipartite_network(partition_probe)
  .run_mcmc_multitoggle_probe(nw)

  cat("\n[NOTE] Pour voir des traces, recompiler avec DEBUG_LOG_FACTORIAL_SIZES=1.\n")
  invisible(TRUE)
}

# --------------------------------------------------------------------------------------
# Jeu de données
# --------------------------------------------------------------------------------------
partitions <- list(
  P1 = c(1, 1, 2, 2, 2, 3),               # tailles: (2,3,1)
  P2 = c(1, 1, 1, 2, 3, 3, 3, 3),         # tailles: (3,1,4)
  P3 = c(1, 2, 2, 3, 3, 4, 4, 4),         # tailles: (1,2,2,3)
  P4 = c(1, 2, 3, 4, 5),                  # toutes tailles 1 (dégénéré)
  P5 = rep(1, 6),                         # un seul groupe (dégénéré)
  P6 = c(1,1,1,1,2,2,3,3,3,4),            # tailles: (4,2,3,1)
  P7 = c(1,2,2,2,2,3,3,4,4,4,4,4)         # tailles: (1,4,2,5)
)

# --------------------------------------------------------------------------------------
# Run principal
# --------------------------------------------------------------------------------------
run_all_tests_log_factorial_sizes <- function() {
  set.seed(42)

  cat("=== SELFTEST ERPM/ERGM: log_factorial_sizes ===\n")
  cat("R:", paste(R.version$major, R.version$minor, sep="."), "\n")
  cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")

  # PHASE 1
  summary_results <- NULL
  if (isTRUE(RUN$phase1_summary)) {
    summary_results <- summary_run_phase1_all(
      partitions = partitions,
      quiet      = isTRUE(RUN$quiet_phase1)
    )
  } else {
    cat("\n=== PHASE 1: SUMMARY ===\n")
    cat("SKIP (désactivée via RUN$phase1_summary = FALSE)\n")
  }

  # PHASE 2
  fit_results <- NULL
  if (isTRUE(RUN$phase2_fit)) {
    fit_results <- erpm_run_phase2_fits(
      partitions = partitions,
      quiet      = isTRUE(RUN$quiet_phase2)
    )
  } else {
    cat("\n=== PHASE 2: ERPM FIT ===\n")
    cat("SKIP (désactivée via RUN$phase2_fit = FALSE)\n")
  }

  # PHASE 3
  if (isTRUE(RUN$phase3_mcmc)) {
    run_phase3_mcmc(partition_probe = partitions$P1)
  } else {
    cat("\n=== PHASE 3: MCMC MULTI-TOGGLE PROBE ===\n")
    cat("SKIP (désactivée via RUN$phase3_mcmc = FALSE)\n")
  }

  invisible(list(summary_results = summary_results, fit_results = fit_results))
}

# Exécution quand lancé en script
if (identical(environment(), globalenv())) {
  run_all_tests_log_factorial_sizes()
}

# Fin: on désactive le patch si présent
if (exists("ergm_patch_disable", mode = "function")) {
  ergm_patch_disable()
}