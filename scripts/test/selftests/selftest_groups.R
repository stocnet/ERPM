# ======================================================================================
# Fichier : scripts/test/selftests/selftest_groups.R
# Objet   : Tests robustes pour l'effet ERPM `groups` → {ergm} `b2degrange`
# Exécution: Rscript scripts/test/selftests/selftest_groups.R
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule: environnement et dépendances minimales
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
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_groups.log")
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
# Helpers réseau via builder du wrapper
# ======================================================================================

# Construit le biparti via build_bipartite_from_inputs(partition=..., nodes=NULL)
make_network_from_partition_via_builder <- function(partition_vec) {
  stopifnot(is.atomic(partition_vec), length(partition_vec) >= 1L)
  built <- build_bipartite_from_inputs(partition = partition_vec, nodes = NULL)
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
# Fonctions utilitaires: normalisation et attentes
# ======================================================================================

# Normaliser la signature groups(...) → (from, to, texte attendu)
normalize_groups_signature <- function(args = list()) {
  if (length(args) == 0L) {
    return(list(from = 1L, to = Inf, text = "b2degrange(from=1,to=Inf)"))
  }
  nm <- names(args)
  if (length(args) == 1L && (is.null(nm) || isTRUE(nm[1L] == ""))) {
    k <- as.numeric(args[[1L]])
    stopifnot(length(k) == 1L, is.finite(k))
    return(list(from = as.integer(k), to = as.integer(k) + 1L,
                text = sprintf("b2degrange(from=%d,to=%d)", as.integer(k), as.integer(k)+1L)))
  }
  if (!is.null(nm) && all(c("from","to") %in% nm)) {
    f <- as.numeric(args[["from"]]); t <- args[["to"]]
    stopifnot(length(f) == 1L)
    f <- as.integer(f)
    txt <- sprintf("b2degrange(from=%s,to=%s)",
                   as.character(f),
                   if (is.infinite(t)) "Inf" else as.character(as.integer(t)))
    return(list(from = f, to = t, text = txt))
  }
  stop("Signature groups(...) inconnue pour ce test.")
}

# Nombre de groupes dont la taille est dans [from, to)
expected_groups_in_range <- function(part, from, to) {
  sz <- as.integer(table(part))
  sum(sz >= from & sz < to)
}

# Vérifier si une chaîne attendue est contenue dans un appel ergm traduit
check_translation_contains <- function(call_ergm, expected) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)
  grepl(expected, compact, fixed = TRUE)
}

# ======================================================================================
# Fonctions explicitement nommées pour SUMMARY et ERPM
# ======================================================================================

# Summary direct: nw construit via builder → summary(nw ~ b2degrange(...))
run_one_network_summary_case_for_groups <- function(partition_vec, groups_args) {
  nw <- make_network_from_partition_via_builder(partition_vec)
  ng <- normalize_groups_signature(groups_args)
  as.integer(suppressMessages(summary(as.formula(nw ~ b2degrange(from = ng$from, to = ng$to)),
                                      constraints = ~ b1part)))
}

# Summary côté ERPM: traduction via erpm(eval_call=FALSE) puis summary sur nw équivalent
run_one_erpm_translated_summary_case_for_groups <- function(partition_vec, call_txt) {
  f <- as.formula(paste0("partition ~ ", call_txt))
  environment(f) <- list2env(list(partition = partition_vec), parent = parent.frame())
  call_ergm <- erpm(f, eval_call = FALSE, verbose = FALSE)

  # Extraire RHS traduit
  ergm_form <- call_ergm[[2L]]
  rhs_expr  <- ergm_form[[3L]]

  # Reconstruire nw via builder et exécuter summary(nw ~ <RHS_traduit>)
  nw2 <- make_network_from_partition_via_builder(partition_vec)
  f2  <- as.formula(bquote(nw2 ~ .(rhs_expr)))
  environment(f2) <- list2env(list(nw2 = nw2), parent = parent.frame())

  # Contrainte de l'appel ou ~b1part
  call_args <- as.list(call_ergm)[-1L]
  cons <- call_args$constraints
  if (is.null(cons)) cons <- as.formula(~ b1part)

  as.integer(suppressMessages(summary(f2, constraints = cons)))
}

# Vérifie l’équivalence exacte des deux summary
check_summary_equivalence_network_vs_erpm_for_groups <- function(partition_vec, cases_df) {
  ok_all <- TRUE
  for (i in seq_len(nrow(cases_df))) {
    call_txt <- cases_df$call_txt[i]
    args     <- cases_df$args[[i]]
    s_net  <- run_one_network_summary_case_for_groups(partition_vec, args)
    s_erpm <- run_one_erpm_translated_summary_case_for_groups(partition_vec, call_txt)
    cat(sprintf("[SUMMARY-CHECK] n=%-3d RHS=%-30s net=%s  erpm=%s\n",
                length(partition_vec), call_txt, paste(s_net, collapse=","), paste(s_erpm, collapse=",")))
    if (!identical(s_net, s_erpm)) {
      ok_all <- FALSE
      cat("  -> MISMATCH détecté.\n")
    }
  }
  ok_all
}

# ======================================================================================
# Panel de partitions et cas groups(...)
# ======================================================================================
partitions <- list(
  P1 = c(1, 2, 2, 3, 3, 3),
  P2 = c(1, 1, 2, 3, 3, 4, 4, 4),
  P3 = c(1, 1, 1, 2, 2, 3),
  P4 = c(1, 2, 3, 4, 5),
  P5 = rep(1, 6)
)

cases <- list(
  list(name = "groups_all",         call_txt = "groups",                args = list()),
  list(name = "groups_exact_1",     call_txt = "groups(1)",             args = list(1)),
  list(name = "groups_exact_2",     call_txt = "groups(2)",             args = list(2)),
  list(name = "groups_exact_3",     call_txt = "groups(3)",             args = list(3)),
  list(name = "groups_range_2_4",   call_txt = "groups(from=2,to=4)",   args = list(from = 2, to = 4)),
  list(name = "groups_range_1_Inf", call_txt = "groups(from=1,to=Inf)", args = list(from = 1, to = Inf))
)
cases_df <- data.frame(
  name     = vapply(cases, `[[`, "", "name"),
  call_txt = vapply(cases, `[[`, "", "call_txt"),
  stringsAsFactors = FALSE
)
cases_df$args <- lapply(cases, `[[`, "args")

# ======================================================================================
# Phase 1: SUMMARY comparatifs (réseau explicite vs ERPM traduit)
# ======================================================================================
run_phase1_summary_equivalence_checks_for_groups <- function() {
  cat("=== PHASE 1 : Summary(nw via builder) vs Summary(ERPM-traduit) — groups ===\n")
  total <- 0L; ok <- 0L
  for (nm in names(partitions)) {
    part <- partitions[[nm]]
    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))
    res <- check_summary_equivalence_network_vs_erpm_for_groups(part, cases_df)
    total <- total + nrow(cases_df)
    ok    <- ok + as.integer(res) * nrow(cases_df)
  }
  cat(sprintf("\n=== Bilan Phase 1 : %d / %d checks OK ===\n", ok, total))
  if (ok < total) stop(sprintf("Summary mismatch sur %d cas.", total - ok))
  invisible(NULL)
}

# ======================================================================================
# Phase 2: Fits via ERPM (MLE + logLik)
# ======================================================================================

# Fit unique ERPM avec retour détaillé
run_one_erpm_fit_with_return_for_groups <- function(partition_vec, rhs_txt, fit_name,
                                                    estimate = "MLE", eval.loglik = TRUE,
                                                    control = control.ergm(MCMC.samplesize = 2000,
                                                                           MCMLE.maxit = 10)) {
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT %-20s] SKIP (erpm() indisponible)\n", fit_name))
    return(list(ok = NA, error = FALSE, coef = NA, fit = NULL))
    }
  set.seed(123)
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition_vec), parent = parent.frame())
  cat(sprintf("[ERPM-FIT %-20s] n=%-3d RHS=%s  | estimate=%s eval.loglik=%s\n",
              fit_name, length(partition_vec), rhs_txt, estimate, as.character(eval.loglik)))
  fit <- try(
    erpm(f, estimate = estimate, eval.loglik = eval.loglik, control = control, verbose = FALSE),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    cat("  -> ERREUR fit:", as.character(fit), "\n")
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

# Lancer plusieurs fits ERPM pour groups et afficher des résumés
run_phase2_erpm_fits_and_print_summaries_for_groups <- function() {
  cat("\n=== PHASE 2 : Fits erpm() (MLE + logLik) — groups ===\n")

  rhs_list <- list(
    R1 = "groups",
    R2 = "groups(1)",
    R3 = "groups(from=2,to=4)",
    R4 = "groups(3) + groups(from=1,to=Inf)"
  )

  parts <- list(
    list(name="P1", part=partitions$P1),
    list(name="P2", part=partitions$P2)
  )

  fit_results <- list()
  for (px in parts) {
    for (nm in names(rhs_list)) {
      key <- paste0(px$name, "_", nm)
      fit_results[[key]] <- run_one_erpm_fit_with_return_for_groups(
        partition_vec = px$part,
        rhs_txt       = rhs_list[[nm]],
        fit_name      = key,
        estimate      = "MLE",
        eval.loglik   = TRUE,
        control       = control.ergm(MCMC.samplesize = 3000, MCMLE.maxit = 10)
      )
    }
  }

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
  invisible(fit_results)
}

# ======================================================================================
# Exécution
# ======================================================================================
set.seed(42)
cat("=== TEST ERPM: groups → b2degrange ===\n")
run_phase1_summary_equivalence_checks_for_groups()
fit_results <- run_phase2_erpm_fits_and_print_summaries_for_groups()

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)
cat("\nTous les tests groups ont passé.\n")