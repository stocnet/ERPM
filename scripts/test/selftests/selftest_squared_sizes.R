# ======================================================================================
# Fichier : scripts/test/selftests/selftest_squared_sizes.R
# Objet   : Self-test autonome pour l'effet ERPM `squared_sizes`
# Exécution: Rscript scripts/test/selftests/selftest_squared_sizes.R
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

# Patch traçage/contournement ERGM si disponible
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) {
    ergm_patch_enable()
  } else {
    message("[ergm_patch] fonction ergm_patch_enable absente")
  }
} else {
  message("[ergm_patch] scripts/ergm_patch.R introuvable, on continue sans patch")
}

# --------------------------------------------------------------------------------------
# Local logging minimal
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
log_path   <- file.path(script_dir, "selftest_squared_sizes.log")
cat("[diag] script_dir =", script_dir, "\n")
cat("[diag] log_path   =", log_path,   "\n")
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

# --------------------------------------------------------------------------------------
# Chargement package/terme + wrapper ERPM
# --------------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  cat("[diag] devtools::load_all()...\n")
  devtools::load_all(quiet = TRUE)
  cat("[diag] load_all OK. Terme attendu présent ? ",
      exists("InitErgmTerm.squared_sizes", mode = "function"), "\n")
} else {
  stop("InitErgmTerm.squared_sizes introuvable. Exécuter depuis le package (load_all).")
}

if (!exists("build_bipartite_from_inputs", mode = "function") || !exists("erpm", mode = "function")) {
  if (file.exists("R/erpm_wrapper.R")) {
    cat("[diag] source R/erpm_wrapper.R\n")
    source("R/erpm_wrapper.R", local = FALSE)
  } else {
    cat("[diag] R/erpm_wrapper.R introuvable, on suppose erpm déjà chargé ailleurs\n")
  }
}
cat("[diag] fonctions wrapper dispo ? build_bipartite_from_inputs=",
    exists("build_bipartite_from_inputs", mode = "function"),
    " erpm=", exists("erpm", mode = "function"), "\n")
if (!exists("build_bipartite_from_inputs", mode = "function"))
  stop("build_bipartite_from_inputs indisponible. Charger R/erpm_wrapper.R.")

# ======================================================================================
# Fonctions locales — noms explicites: *summary* / *erpm*
# ======================================================================================

# Construction du biparti via le wrapper (exigé)
.build_bipartite_nw_via_wrapper <- function(part) {
  stopifnot(length(part) >= 1L)
  nodes <- data.frame(label = utils::head(LETTERS, length(part)),
                      stringsAsFactors = FALSE)

  if (!exists("build_bipartite_from_inputs", mode = "function")) {
    stop("build_bipartite_from_inputs indisponible (charger R/erpm_wrapper.R).")
  }
  builder <- get("build_bipartite_from_inputs", mode = "function")

  try_sigs <- list(
    list(expr = quote(builder(partition = part, nodes = nodes)),                                   tag = "partition+nodes"),
    list(expr = quote(builder(partition = part, labels = nodes$label, attributes = list())),        tag = "partition+labels+attributes"),
    list(expr = quote(builder(partition = part)),                                                   tag = "partition seul"),
    list(expr = quote(builder(partition = part, labels = nodes$label)),                             tag = "partition+labels")
  )

  last_err <- NULL
  for (sig in try_sigs) {
    cat("[diag] build_bipartite_from_inputs essai:", sig$tag, "\n")
    out <- try(eval(sig$expr), silent = TRUE)
    if (!inherits(out, "try-error") && !is.null(out)) {
      if (inherits(out, "network")) {
        cat("[diag] build OK (retour = network) via", sig$tag, "\n")
        return(out)
      }
      if (is.list(out)) {
        for (nm in c("network","nw","net","graph","g","bip")) {
          if (!is.null(out[[nm]]) && inherits(out[[nm]], "network")) {
            cat("[diag] build OK (retour = list$", nm, ") via", sig$tag, "\n")
            return(out[[nm]])
          }
        }
      }
    } else {
      last_err <- out
      cat("[diag] build FAIL via", sig$tag, ":", if (!is.null(out)) as.character(out)[1] else "<NA>", "\n")
    }
  }

  fm <- try(formals(builder), silent = TRUE)
  stop(paste0(
    "build_bipartite_from_inputs a échoué avec toutes les signatures testées.\n",
    "Formals: ", if (!inherits(fm, "try-error")) paste(names(fm), collapse = ", ") else "<inconnus>", "\n",
    "Dernière erreur: ", if (!is.null(last_err)) as.character(last_err)[1] else "<aucune>"
  ))
}

# Tailles de groupes depuis partition
.group_sizes_from_partition <- function(part) as.integer(table(part))

# Attendue: somme des tailles^pow filtrées par [from,to)
.expected_summary_squared_sizes_from_partition <- function(part, from = 1, to = Inf, pow = 2) {
  sz  <- .group_sizes_from_partition(part)
  idx <- (sz >= from) & (sz < to)
  sum((sz[idx])^pow)
}

# Identité pow2 sur réseau (contrôle)
.identity_pow2_all_summary_on_network <- function(nw) {
  stopifnot(inherits(nw, "network"))
  n1 <- network::get.network.attribute(nw, "bipartite")
  if (is.null(n1) || is.na(n1)) stop("Réseau non biparti.")
  n  <- network::network.size(nw)
  v2 <- seq.int(as.integer(n1) + 1L, n)
  deg2 <- vapply(v2, function(v) length(network::get.neighborhood(nw, v, type = "combined")), integer(1L))
  ecount <- network::network.edgecount(nw)
  ecount + 2L * sum(choose(deg2, 2))
}

# Normalisation de la signature
.normalize_squared_sizes_signature <- function(args = list()) {
  out <- list(from = 1, to = Inf, pow = 2, text = "squared_sizes(from=1,to=Inf,pow=2)")
  if (length(args)) {
    nm <- names(args)
    if (!is.null(nm) && length(nm)) {
      if ("from" %in% nm) out$from <- as.numeric(args[["from"]])
      if ("to"   %in% nm) out$to   <- as.numeric(args[["to"]])
      if ("pow"  %in% nm) out$pow  <- as.numeric(args[["pow"]])
    }
  }
  out$text <- sprintf("squared_sizes(from=%s,to=%s,pow=%s)",
                      if (is.infinite(out$from)) "Inf" else as.character(out$from),
                      if (is.infinite(out$to))   "Inf" else as.character(out$to),
                      as.character(out$pow))
  out
}

# Vérification que l'appel erpm traduit contient la bonne signature
.check_translation_ok_erpm <- function(call_ergm, term = "squared_sizes", args = list()) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)
  def <- list(from = 1, to = Inf, pow = 2)
  if (!grepl(paste0("\\b", term, "\\("), compact)) return(FALSE)
  checks <- logical(0)
  if (!is.null(args$from) && !identical(as.numeric(args$from), def$from)) {
    checks <- c(checks, grepl(paste0("from=", gsub("\\s+", "", as.character(args$from))), compact, fixed = TRUE))
  }
  if (!is.null(args$to) && !(is.infinite(args$to) && is.infinite(def$to)) &&
      !identical(as.numeric(args$to), def$to)) {
    checks <- c(checks, grepl(paste0("to=", gsub("\\s+", "", as.character(args$to))), compact, fixed = TRUE))
  }
  if (!is.null(args$pow) && !identical(as.numeric(args$pow), def$pow)) {
    checks <- c(checks, grepl(paste0("pow=", gsub("\\s+", "", as.character(args$pow))), compact, fixed = TRUE))
  }
  if (!length(checks)) return(TRUE)
  all(checks)
}

# ---------- Couche SUMMARY : un cas ----------
.run_one_case_summary_squared_sizes <- function(part, name, call_txt, args) {
  nw <- .build_bipartite_nw_via_wrapper(part)
  sg <- .normalize_squared_sizes_signature(args)

  expected_val <- .expected_summary_squared_sizes_from_partition(part, from = sg$from, to = sg$to, pow = sg$pow)

  f <- as.formula(paste0("nw ~ ", call_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  stat_val <- as.numeric(summary(f))

  ok_trad <- NA
  if (exists("erpm", mode = "function")) {
    call_ergm <- erpm(f, eval_call = FALSE, verbose = TRUE)
    ok_trad <- .check_translation_ok_erpm(call_ergm, term = "squared_sizes", args = args)
  }

  ok_identity <- NA
  if (isTRUE(sg$pow == 2 && sg$from == 1 && isTRUE(is.infinite(sg$to)))) {
    id_val <- .identity_pow2_all_summary_on_network(nw)
    ok_identity <- identical(unname(as.integer(stat_val)), unname(as.integer(id_val)))
  }

  cat(sprintf("\n[SUMMARY CAS %-18s] part={%s}", name, paste(part, collapse=",")))
  cat(sprintf("\tappel=%s", sg$text))
  cat(sprintf("\tsummary(.)=%s", format(stat_val)))
  cat(sprintf("\tattendu(part)=%s", format(expected_val)))
  if (!is.na(ok_identity)) cat(sprintf("\tidentité pow2=%s", if (ok_identity) "OK" else "KO"))
  if (!is.na(ok_trad))     cat(sprintf("\ttrad erpm=%s",    if (ok_trad) "OK" else "KO"))
  cat("\n")

  list(
    ok_stat   = identical(unname(as.integer(stat_val)), unname(as.integer(expected_val))),
    ok_trad   = ok_trad,
    ok_ident  = ok_identity,
    stat      = stat_val,
    expected  = expected_val
  )
}

# ---------- Couche SUMMARY : panneau ----------
.run_cases_for_partition_summary_squared_sizes <- function(part, panel) {
  res <- lapply(panel, function(cx) {
    out <- .run_one_case_summary_squared_sizes(part, cx$name, cx$call_txt, cx$args)
    data.frame(
      case     = cx$name,
      ok_stat  = out$ok_stat,
      ok_trad  = if (is.na(out$ok_trad)) NA else out$ok_trad,
      ok_ident = if (is.na(out$ok_ident)) NA else out$ok_ident,
      stat     = out$stat,
      expected = out$expected,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

# ---------- Couche ERPM FIT : un cas ----------
.run_one_case_erpm_fit_squared_sizes <- function(part, name, call_txt, ctrl) {
  cat(sprintf("\n[ERPM-FIT %s] début\n", name))

  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT %-18s] SKIP (erpm() indisponible)\n", name))
    return(list(ok = NA, error = TRUE, coef = NA, fit = NULL))
  }

  # Inspecte la signature d’erpm
  frm <- try(formals(erpm), silent = TRUE)
  if (inherits(frm, "try-error")) {
    cat("[diag] formals(erpm) -> ERREUR\n")
    frm_names <- character(0)
  } else {
    frm_names <- names(frm)
    cat("[diag] formals(erpm) noms =", paste(frm_names, collapse = ", "), "\n")
  }
  has_labels     <- "labels"     %in% frm_names
  has_attributes <- "attributes" %in% frm_names
  cat(sprintf("[diag] support labels=%s, attributes=%s\n", has_labels, has_attributes))

  # Formule ERPM avec partition en LHS
  f <- as.formula(paste0("part ~ ", call_txt))
  environment(f) <- list2env(list(part = part), parent = parent.frame())
  cat("[diag] formule ERPM =", paste(deparse(f, width.cutoff = 500L), collapse = " "), "\n")

  # Prépare arguments dynamiquement
  lab <- utils::head(LETTERS, length(part))
  arglist <- list(
    formula     = f,
    # estimate    = "MLE",
    # control     = ctrl,
    eval.loglik = TRUE,
    verbose     = TRUE
  )
  if (has_labels)     arglist$labels     <- lab
  if (has_attributes) arglist$attributes <- list()

  cat("[diag] arguments passés à erpm(): ", paste(names(arglist), collapse = ", "), "\n")

  # Dry-run pour voir l’appel ergm final
  arglist_dry <- arglist
  arglist_dry$eval_call <- FALSE
  dry <- try(do.call(erpm, arglist_dry), silent = TRUE)
  if (inherits(dry, "try-error")) {
    cat("[diag] dry-run erpm -> ERREUR:\n", as.character(dry)[1], "\n")
    return(list(ok = FALSE, error = TRUE, coef = NA, fit = NULL))
  }
  dry_s <- paste(deparse(dry, width.cutoff = 500L), collapse = " ")
  cat("[diag] dry-run appel ergm =", dry_s, "\n")
  dl <- as.list(dry); nm <- names(dl)
  has_constraints <- !is.null(nm) && any(nm == "constraints")
  cat("[diag] contraintes présentes dans l’appel dry-run: ", has_constraints, "\n")

  # Fit réel via erpm
  res <- try(do.call(erpm, arglist), silent = TRUE)
  if (inherits(res, "try-error")) {
    cat(sprintf("[ERPM-FIT %-18s] ERREUR évaluation: %s\n", name, as.character(res)[1]))
    tb <- try(utils::traceback(x = NULL, max.lines = 1), silent = TRUE)
    if (!inherits(tb, "try-error")) cat("[diag] traceback court affiché ci-dessus s’il existe.\n")
    return(list(ok = FALSE, error = TRUE, coef = NA, fit = NULL))
  }

  cf <- try(stats::coef(res), silent = TRUE)
  ok_coef  <- !(inherits(cf, "try-error")) && all(is.finite(cf))
  ok_class <- inherits(res, "ergm")
  cat(sprintf("[ERPM-FIT %-18s] coef finies: %s | coef=%s\n",
              name, if (ok_coef) "OK" else "KO",
              if (ok_coef) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

  list(ok = ok_class && ok_coef, error = FALSE, coef = if (ok_coef) cf else NA, fit = res)
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
  list(name="sq_all_pow2",    call_txt="squared_sizes",                   args=list(from=1, to=Inf, pow=2)),
  list(name="sq_2to5_pow2",   call_txt="squared_sizes(from=2,to=5)",      args=list(from=2, to=5,   pow=2)),
  list(name="sq_all_pow3",    call_txt="squared_sizes(pow=3)",            args=list(from=1, to=Inf, pow=3)),
  list(name="sq_1to2_pow2",   call_txt="squared_sizes(from=1,to=2)",      args=list(from=1, to=2,   pow=2)),
  list(name="sq_3toInf_pow2", call_txt="squared_sizes(from=3,to=Inf)",    args=list(from=3, to=Inf, pow=2))
)

# Contrôle MCMLE commun
ctrl_fit <- control.ergm(
  MCMLE.maxit     = 10,
  MCMC.samplesize = 1e4
)

# ======================================================================================
# Run principal — noms explicites
# ======================================================================================
run_all_tests_squared_sizes <- function() {
  # Logging centralisé côté repo
  log_path2 <- file.path("scripts","test","selftests","selftest_squared_sizes.log")
  dir.create(dirname(log_path2), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(log_path2)) unlink(log_path2)
  sink(file = log_path2, split = TRUE)
  con_msg <- file(log_path2, open = "at")
  sink(con_msg, type = "message")
  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(close(con_msg),        silent = TRUE)
    try(sink(),                silent = TRUE)
    flush.console()
  }, add = TRUE)

  set.seed(42)
  cat("=== TEST ERPM: squared_sizes ===\n")
  cat("[diag] R.version =", paste(R.version$major, R.version$minor, sep="."), "\n")
  cat("[diag] ergm version =", as.character(utils::packageVersion("ergm")), "\n")

  all_summary_results <- list()
  fit_results <- list()

  total_ok <- 0L; total_n <- 0L

  # ---------- Phase SUMMARY sur toutes partitions ----------
  for (nm in names(partitions)) {
    cat(sprintf("\n--- SUMMARY Partition %s ---\n", nm))
    df <- .run_cases_for_partition_summary_squared_sizes(partitions[[nm]], cases)
    all_summary_results[[nm]] <- df
    total_ok <- total_ok + sum(df$ok_stat,  na.rm = TRUE); total_n <- total_n + sum(!is.na(df$ok_stat))
    total_ok <- total_ok + sum(df$ok_ident, na.rm = TRUE); total_n <- total_n + sum(!is.na(df$ok_ident))
    total_ok <- total_ok + sum(df$ok_trad,  na.rm = TRUE); total_n <- total_n + sum(!is.na(df$ok_trad))
    print(df)
  }

  cat(sprintf("\n=== Bilan global SUMMARY : %d / %d validations OK ===\n", total_ok, total_n))
  if (total_ok < total_n) stop(sprintf("Echec SUMMARY: %d validations KO", total_n - total_ok))

  # ---------- Phase ERPM FIT sur quelques cas représentatifs ----------
  cat("\n=== PHASE ERPM FIT: estimation MLE + loglik ===\n")
  part_ref <- partitions$P1
  fit_results[["FIT_sq_all_pow2"]]  <- .run_one_case_erpm_fit_squared_sizes(part_ref, "sq_all_pow2",  "squared_sizes", ctrl_fit)
  fit_results[["FIT_sq_2to5_pow2"]] <- .run_one_case_erpm_fit_squared_sizes(part_ref, "sq_2to5_pow2", "squared_sizes(from=2,to=5)", ctrl_fit)
  fit_results[["FIT_sq_all_pow3"]]  <- .run_one_case_erpm_fit_squared_sizes(part_ref, "sq_all_pow3",  "squared_sizes(pow=3)", ctrl_fit)

  ok <- vapply(fit_results, function(x) isTRUE(x$ok), logical(1))
  n_ok <- sum(ok, na.rm = TRUE); n_tot <- sum(!is.na(ok))
  cat(sprintf("\n=== Bilan fits erpm() : %d / %d OK ===\n", n_ok, n_tot))

  cat("\n=== Résumés détaillés des fits ERPM réussis ===\n")
  for (nm in names(fit_results)) {
    fit_obj <- fit_results[[nm]]
    if (isTRUE(fit_obj$ok) && inherits(fit_obj$coef, "numeric")) {
      cat(sprintf("\n--- Résumé fit %s ---\n", nm))
      print(summary(fit_obj$fit))
    }
  }

  if (n_ok < n_tot) stop(sprintf("Echec fits: %d KO", n_tot - n_ok))

  invisible(list(summary_results = all_summary_results, fit_results = fit_results))
}

if (identical(environment(), globalenv())) run_all_tests_squared_sizes()
ergm_patch_disable()