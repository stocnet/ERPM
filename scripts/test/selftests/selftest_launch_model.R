# ==============================================================================
# Fichier : selftest_launch_model.R
# Objet   : Batterie de tests robuste pour launch_model() (summary | ergm | erpm)
# Auteur  : réécriture avec journalisation fiable et sortie sans couleurs
# ==============================================================================

# macOS/Linux : empêcher les séquences ANSI et forcer FR UTF-8
Sys.setenv(LANG = "fr_FR.UTF-8",
           CLI_NO_COLOR="1", NO_COLOR="1", R_CLI_NUM_COLORS="0", CRAYON_DISABLE="1", TERM="dumb")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))
options(crayon.enabled = FALSE, cli.num_colors = 1, warn = 1)

# --- Préambule ---------------------------------------------------------------
# Dossier de logs + duplication sortie console
log_dir <- "scripts/test"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file_out <- file.path(log_dir, "selftest_launch_model.log")

sink(log_file_out, split = TRUE)      # stdout -> fichier + console
sink(stdout(), type = "message")      # messages -> stdout (donc même fichier)

# Chargement strict des dépendances
.pkgs <- c("ergm", "network")
.miss <- .pkgs[!vapply(.pkgs, requireNamespace, FALSE, quietly = TRUE)]
if (length(.miss)) stop("Packages manquants: ", paste(.miss, collapse = ", "))
suppressPackageStartupMessages({
  lapply(.pkgs, require, character.only = TRUE)
})

# Sources locales nécessaires
source("scripts/local_source/settings.R",      local = FALSE)
source("scripts/local_source/init.R",          local = FALSE)
source("scripts/local_source/launcher.R",      local = FALSE)
source("R/erpm_wrapper.R",                     local = FALSE)
source("R/functions_erpm_bip_network.R",       local = FALSE)

# Initialisation du package ERPM
init_erpm(selftest = FALSE, verbose = TRUE)

# ------------------------------------------------------------------------------
#' @title .strip_ansi_file
#' @description Supprime les séquences ANSI (couleurs, gras, etc.) d’un fichier log.
#' @param path Chemin du fichier à nettoyer.
#' @return Invisible. Écrit le fichier nettoyé sur place.
.strip_ansi_file <- function(path) {
  if (!file.exists(path)) return(invisible())
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) return(character()))
  if (!length(txt)) return(invisible())
  if (requireNamespace("fansi", quietly = TRUE)) {
    txt <- fansi::strip_sgr(txt)
  } else if (requireNamespace("crayon", quietly = TRUE)) {
    txt <- crayon::strip_style(txt)
  } else {
    txt <- gsub("\\x1b\\[[0-9;]*[mK]", "", txt, perl = TRUE)
  }
  writeLines(txt, path)
  invisible()
}

# ------------------------------------------------------------------------------
#' @title .onexit_cleaners
#' @description Nettoyage complet à la fin d’exécution (désactive patchs ERGM,
#'   vide les environnements, ferme les sinks, supprime les séquences ANSI).
#' @return Invisible.
.onexit_cleaners <- function() {
  if (exists("ergm_patch_disable", mode = "function")) try(ergm_patch_disable(), silent = TRUE)
  if (exists("clean_global_env",  mode = "function")) try(clean_global_env(),  silent = TRUE)
  while (sink.number() > 0) sink(NULL)
  while (sink.number(type = "message") > 0) sink(NULL, type = "message")
  try(.strip_ansi_file(log_file_out), silent = TRUE)
}

CLEAN_ON_EXIT <- isTRUE(as.logical(Sys.getenv("ERPM_CLEAN_ON_EXIT", "false"))) || !interactive()
if (CLEAN_ON_EXIT) on.exit(.onexit_cleaners(), add = TRUE)

# Reproductibilité
set.seed(1234)

# ------------------------------------------------------------------------------
#' @title .with_null_device
#' @description Exécute du code graphique dans un device PDF temporaire (évite
#'   l’ouverture d’une fenêtre graphique pendant les tests).
#' @param code Expression R à exécuter.
#' @return Résultat de `code`.
.with_null_device <- function(code) {
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({ try(grDevices::dev.off(), silent = TRUE); unlink(tf) }, add = TRUE)
  force(code)
}

# ------------------------------------------------------------------------------
#' @title .summarize_ret
#' @description Affiche un résumé synthétique d’un résultat de `launch_model()`.
#' @param ret Résultat retourné par `launch_model()`.
#' @return Invisible. Écrit dans la console/log.
.summarize_ret <- function(ret) {
  if (!is.list(ret)) return(invisible(NULL))
  if (!is.null(ret$fit) && inherits(ret$fit, "ergm")) {
    k  <- tryCatch(length(stats::coef(ret$fit)), error = function(e) NA_integer_)
    try(suppressWarnings(logLik(ret$fit, add = TRUE)), silent = TRUE)
    ll <- tryCatch(as.numeric(stats::logLik(ret$fit)), error = function(e) NA_real_)
    cat(sprintf("   -> ergm fit | #coef=%s | logLik=%s\n",
                ifelse(is.na(k), "NA", k),
                ifelse(is.na(ll), "NA", sprintf("%.4f", ll))))
    tn <- tryCatch(ret$fit$termnames, error = function(e) NULL)
    if (length(tn)) cat("   -> termes: ", paste(tn, collapse = ", "), "\n", sep = "")
  } else if (!is.null(ret$result)) {
    cat(sprintf("   -> result class: %s\n", paste(class(ret$result), collapse = ",")))
  } else if (!is.null(ret$call_text)) {
    cat(sprintf("   -> call: %s\n", ret$call_text))
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
#' @title .is_degenerate
#' @description Détecte une dégénérescence d’un modèle ERGM.
#' @param ret Résultat d’un appel `launch_model("ergm", ...)`.
#' @return Booléen indiquant si le modèle est dégénéré.
.is_degenerate <- function(ret) {
  if (!is.list(ret) || is.null(ret$fit) || !inherits(ret$fit, "ergm")) return(FALSE)
  bad_coef <- tryCatch(any(!is.finite(stats::coef(ret$fit))), error = function(e) TRUE)
  msg <- tolower(paste(
    tryCatch(paste(ret$fit$optim$output, collapse = " "), error = function(e) ""),
    tryCatch(paste(ret$fit$optim$info,   collapse = " "), error = function(e) "")
  ))
  bad_mcmc <- grepl("non[ -]?converg|degener|infinite coefficient", msg)
  bad_coef || bad_mcmc
}

# ------------------------------------------------------------------------------
#' @title .run_test
#' @description Enveloppe d’exécution standardisée pour les tests unitaires.
#'   Gère les erreurs, warnings, temps d’exécution et validation optionnelle.
#' @param label Nom du test (affiché dans les logs).
#' @param expr Expression R à exécuter.
#' @param expect_error Si TRUE, considère qu’une erreur est attendue.
#' @param check_ok Fonction de validation additionnelle.
#' @param warnings_are_fail Si TRUE, les warnings font échouer le test.
#' @return Résultat de `expr`, invisible.
.run_test <- function(label, expr, expect_error = FALSE, check_ok = NULL, warnings_are_fail = TRUE) {
  cat("\n===== TEST: ", label, " =====\n", sep = "")
  t0 <- proc.time()[["elapsed"]]
  w <- character()
  out <- withCallingHandlers(
    tryCatch(force(expr),
             error = function(e) structure(list(`__err__` = TRUE, message = conditionMessage(e)), class = "test_error")),
    warning = function(m) {
      msg <- conditionMessage(m)
      harmless <- grepl("Rglpk.*falling back to.*lpSolveAPI", msg)
      if (!harmless && isTRUE(warnings_are_fail)) w <<- c(w, msg)
    }
  )
  dt <- proc.time()[["elapsed"]] - t0

  if (inherits(out, "test_error")) {
    if (isTRUE(expect_error)) {
      cat("✅ Erreur attendue capturée: ", out$message, "\n", sep = "")
    } else {
      cat("❌ ERREUR: ", out$message, "\n", sep = "")
    }
    cat(sprintf("   -> durée: %.2fs\n", dt))
    return(invisible(out))
  }

  ok <- TRUE
  if (length(w)) ok <- FALSE
  if (is.function(check_ok)) ok <- isTRUE(tryCatch(check_ok(out), error = function(e) FALSE))
  if (isTRUE(expect_error)) ok <- FALSE

  if (ok) {
    cat("✅ OK\n")
  } else {
    cat("❌ Échec\n")
    if (length(w)) cat("   -> warning(s): ", paste(w, collapse = " | "), "\n", sep = "")
  }
  .summarize_ret(out)
  cat(sprintf("   -> durée: %.2fs\n", dt))
  invisible(out)
}

# ------------------------------------------------------------------------------
#' @title create_demo
#' @description Crée les données de démonstration utilisées pour les tests.
#' @param seed Graine de reproductibilité.
#' @return Liste contenant un réseau biparti et une partition.
create_demo <- function(seed = 1234) {
  if (!exists("create_erpm_network", mode = "function")) {
    stop("create_erpm_network() introuvable. Vérifie 'R/functions_erpm_bip_network.R'.")
  }
  set.seed(seed)
  demo <- create_erpm_network()
  list(nw = demo$network, part = demo$partition)
}

.demo <- create_demo()
nw_demo   <- .demo$nw
part_demo <- .demo$part


# ==============================================================================
# TESTS SUMMARY (statistiques rapides)
# ==============================================================================

ret1  <- .run_test("SUMMARY | b2degrange simple", {
  launch_model("summary", effects = "b2degrange(from=2,to=3)")
})

ret2  <- .run_test("SUMMARY | multi-termes avec squared_sizes", {
  launch_model("summary", effects = "b2degrange(from=2,to=3) + squared_sizes(from=1,to=2)")
})

ret3  <- .run_test("SUMMARY | RHS par liste d'appels", {
  launch_model("summary",
               effects = list(
                 call("b2degrange", from = 2, to = 4),
                 call("b2star", 2)
               ))
})

ret4  <- .run_test("SUMMARY | dry-run pour inspecter la formule", {
  launch_model("summary", effects = "b2star(2) + b2degrange(from=1,to=3)", dry_run = TRUE)
})

# ==============================================================================
# TESTS ERGM (estimate doit être 'CD' ici + variantes)
# ==============================================================================

if (exists("ergm_patch_enable")) ergm_patch_enable(verbose = TRUE)

ret5  <- .run_test("ERGM | CD basique biparti (logLik)", {
  set.seed(1234)
  launch_model("ergm",
               effects     = "b2degrange(from=2,to=3)",
               estimate    = "CD",
               eval_loglik = TRUE)
})

ret6  <- .run_test("ERGM | CD deux termes + contrainte explicite ~b1part (logLik)", {
  set.seed(1234)
  launch_model("ergm",
               effects     = "b2degrange(from=2,to=3) + b2star(2)",
               constraints = ~ b1part,
               estimate    = "CD",
               eval_loglik = TRUE)
})

ret7  <- .run_test("ERGM | CD timeout court (pas de logLik)", {
  set.seed(1234)
  launch_model("ergm",
               effects     = "b2degrange(from=1,to=5)",
               estimate    = "CD",
               eval_loglik = FALSE,
               timeout     = 5)
})

ret8  <- .run_test("ERGM | CD + plot + verbose (logLik)", {
  set.seed(1234)
  .with_null_device({
    launch_model("ergm",
                 effects     = "b2degrange(from=2,to=5)",
                 estimate    = "CD",
                 eval_loglik = TRUE,
                 plot        = TRUE,
                 verbose     = TRUE)
  })
})

ret9  <- .run_test("ERGM | CD avec réseau fourni (nw_demo) (logLik)", {
  set.seed(1234)
  launch_model("ergm",
               effects     = "b2star(2)",
               nw          = nw_demo,
               estimate    = "CD",
               eval_loglik = TRUE)
})

ret10 <- .run_test("ERGM | CD effets en liste d'appels (logLik)", {
  set.seed(1234)
  launch_model("ergm",
               effects     = list(
                 call("b2degrange", from = 2, to = 3),
                 call("squared_sizes", from = 1, to = 2)
               ),
               estimate    = "CD",
               eval_loglik = TRUE)
})

ret11 <- .run_test("ERGM | dry-run pour inspecter l'appel ergm()", {
  launch_model("ergm",
               effects     = "b2degrange(from=2,to=3) + b2star(3)",
               estimate    = "CD",
               eval_loglik = FALSE,
               dry_run     = TRUE)
})

# Variante supplémentaire
ret11b <- .run_test("ERGM | MPLE + contrainte observed (pas de logLik)", {
  set.seed(1234)
  launch_model("ergm",
               effects     = "b2degrange(from=2,to=3)",
               constraints = ~ observed,
               estimate    = "MPLE",
               eval_loglik = FALSE)
})

# ==============================================================================
# TESTS ERPM (traduction + encapsulations)
# ==============================================================================

ret12 <- .run_test("ERPM | groups + squared_sizes en CD (logLik)", {
  set.seed(1234)
  launch_model("erpm",
               effects     = "groups(from=2,to=3) + squared_sizes(from=1,to=2)",
               estimate    = "CD",
               eval_loglik = TRUE)
})

ret13 <- .run_test("ERPM | dry-run traduction ERPM→ERGM (regex)", {
  out <- launch_model("erpm",
                      effects = "groups(from=1,to=3) + squared_sizes(from=1,to=2)",
                      dry_run = TRUE)
  call_txt <- tryCatch(out$call_text, error = function(e) NULL)
  if (length(call_txt)) {
    stopifnot(grepl("b2", paste(call_txt, collapse = " ")))
  } else {
    fml <- tryCatch(out$result, error = function(e) NULL)
    stopifnot(inherits(fml, "formula"))
    stopifnot(grepl("b2", paste(deparse(fml), collapse = " ")))
  }
  out
})

ret14 <- .run_test("ERPM | CD contrôle MCMC allégé (logLik)", {
  set.seed(1234)
  launch_model("erpm",
               effects     = "groups(from=2,to=4) + squared_sizes(from=1,to=3)",
               estimate    = "CD",
               eval_loglik = TRUE,
               control     = list(MCMLE.maxit = 3, MCMC.samplesize = 1000))
})

ret15 <- .run_test("ERPM | CD avec timeout (pas de logLik)", {
  set.seed(1234)
  launch_model("erpm",
               effects     = "groups(from=2,to=5) + squared_sizes(from=1,to=2)",
               estimate    = "CD",
               eval_loglik = FALSE,
               timeout     = 10)
})

ret16 <- .run_test("ERPM | CD avec réseau fourni (nw_demo) (logLik)", {
  set.seed(1234)
  launch_model("erpm",
               effects     = "groups(from=2,to=3)",
               nw          = nw_demo,
               estimate    = "CD",
               eval_loglik = TRUE)
})

ret17 <- .run_test("ERPM | CD contraintes auto (~b1part) implicites (logLik)", {
  set.seed(1234)
  launch_model("erpm",
               effects     = "squared_sizes(from=1,to=3)",
               estimate    = "CD",
               eval_loglik = TRUE)
})

# ==============================================================================
# CAS SPÉCIAUX / ROBUSTESSE
# ==============================================================================

ret18 <- .run_test(
  "ROBUSTESSE | ERGM terme inadapté au biparti (dégénérescence attendue)",
  {
    set.seed(1234)
    launch_model("ergm",
                 effects     = "triangles",
                 estimate    = "CD",
                 eval_loglik = FALSE)
  },
  expect_error = FALSE,
  check_ok     = .is_degenerate
)

ret19 <- .run_test("ROBUSTESSE | SUMMARY via build_rhs_from_effects()", {
  rhs_call <- build_rhs_from_effects(
    effects      = c("b2degrange", "squared_sizes"),
    effects_args = list(
      b2degrange    = list(from = 2, to = 4),
      squared_sizes = list(from = 1, to = 2)
    )
  )
  launch_model("summary", effects = rhs_call)
})

ret20 <- .run_test("ROBUSTESSE | ERGM CD contrôle custom (pas de logLik)", {
  set.seed(1234)
  launch_model("ergm",
               effects      = "b2degrange(from=2,to=3)",
               estimate     = "CD",
               eval_loglik  = FALSE,
               control      = list(MCMLE.maxit = 2, MCMC.samplesize = 800, MCMC.burnin = 5000))
})

# Résumé d'exécution
cat("\nTous les tests ont été exécutés.\n\n")

# ==============================================================================
# RÉCAPITULATIF GLOBAL
# ==============================================================================

all_ret <- list(
  ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10,
  ret11, ret11b, ret12, ret13, ret14, ret15, ret16, ret17, ret18, ret19, ret20
)

summary_table <- data.frame(
  test_id = sprintf("ret%-2d", seq_along(all_ret)),
  label   = c(
    "SUMMARY simple", "SUMMARY multi-termes", "SUMMARY RHS liste",
    "SUMMARY dry-run", "ERGM CD basique", "ERGM CD + contraintes",
    "ERGM CD timeout", "ERGM CD + plot", "ERGM CD + nw fourni",
    "ERGM CD effets liste", "ERGM dry-run", "ERGM MPLE + observed",
    "ERPM CD groups + squared_sizes", "ERPM dry-run traduction",
    "ERPM CD contrôle léger", "ERPM CD timeout", "ERPM CD + nw fourni",
    "ERPM CD contraintes auto", "ERGM terme inadapté", "SUMMARY via build_rhs",
    "ERGM CD contrôle custom"
  ),
  status  = vapply(all_ret, function(x) {
    if (inherits(x, "test_error")) {
      if (isTRUE(x$`__err__`)) "ERREUR" else "Inconnu"
    } else if (.is_degenerate(x)) {
      "Dégénéré"
    } else {
      "OK"
    }
  }, character(1)),
  stringsAsFactors = FALSE
)

cat("\n\n=========================\n")
cat("RÉSUMÉ GLOBAL DES TESTS\n")
cat("=========================\n\n")
print(summary_table, right = FALSE, row.names = FALSE)

n_ok  <- sum(summary_table$status == "OK")
n_err <- sum(summary_table$status == "ERREUR")
n_deg <- sum(summary_table$status == "Dégénéré")
cat("\nTotal OK :", n_ok, "/", nrow(summary_table), "\n")
cat("Total erreurs :", n_err, "\n")
cat("Total dégénérés :", n_deg, "\n")

# Code de sortie non nul si échecs (utile en CI)
if (!interactive()) {

  # Ferme proprement les sinks
  while (sink.number(type="message") > 0) sink(NULL, type="message")
  while (sink.number() > 0) sink(NULL)

  # Strip ANSI et, si voulu, remplace les symboles Unicode
  try({
    .strip_ansi_file(log_file_out)
    txt <- readLines(log_file_out, warn = FALSE)
    # Optionnel: ASCII only
    txt <- chartr("ℹ✔✖", "i+ x", txt)
    writeLines(txt, log_file_out)
  }, silent = TRUE)

  quit(save = "no", status = if ((n_err + n_deg) > 0) 1L else 0L)
}
