# ======================================================================================
# Fichier : scripts/test/test_groups.R
# Objet   : Tests robustes pour l'effet ERPM `groups` → {ergm} `b2degrange`
# Exécution: Rscript scripts/test/test_groups.R
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule: environnement et dépendances minimales
# --------------------------------------------------------------------------------------
# macOS/Linux: forcer UTF-8 pour affichages
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

# --- Logging console + fichier (indépendant du cwd) ---
get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 1L) return(normalizePath(dirname(f), winslash = "/", mustWork = FALSE))
  # fallback: si source() en interactif
  if (!is.null(sys.frames()) && !is.null(sys.calls())) {
    for (i in rev(seq_along(sys.calls()))) {
      cf <- sys.frame(i)
      if (!is.null(cf$ofile)) return(normalizePath(dirname(cf$ofile), winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
log_path   <- file.path(script_dir, "test_groups.log")

# réinit
if (file.exists(log_path)) unlink(log_path, force = TRUE)

# connexions dédiées
con_out <- file(log_path, open = "wt")  # sortie standard
con_err <- file(log_path, open = "at")  # messages/erreurs en append

sink(con_out, split = TRUE)             # stdout → fichier + console
sink(con_err, type = "message")         # messages → fichier

on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(close(con_err),        silent = TRUE)
  try(sink(),                silent = TRUE)
  try(close(con_out),        silent = TRUE)
  flush.console()
}, add = TRUE)

# Charger les fonctions biparti si non chargés
if (!exists("partition_to_bipartite_network", mode = "function")) {
  if (file.exists("R/functions_erpm_bip_network.R")) {
    source("R/functions_erpm_bip_network.R", local = FALSE)
  } else {
    stop("functions_erpm_bip_network.R introuvable. Place ce fichier sous R/.")
  }
}

# Charger le wrapper ERPM si non chargé
if (!exists("erpm", mode = "function")) {
  if (file.exists("R/erpm_wrapper.R")) {
    source("R/erpm_wrapper.R", local = FALSE)
  } else {
    stop("erpm_wrapper.R introuvable. Place ce fichier sous R/.")
  }
}

# ======================================================================================
# Fonctions locales
# ======================================================================================

#' Construire un réseau biparti depuis une partition
#' @param part vecteur de groupes par objet, ex. c(1,2,2,3,3,3).
#' @return objet `network` bipartite prêt pour {ergm}.
#' @keywords internal
make_nw_from_partition <- function(part) {
  n <- length(part)                           # nb objets
  stopifnot(n >= 1, is.atomic(part))          # garde-fous
  lbl <- utils::head(LETTERS, n)              # labels "A","B",...
  partition_to_bipartite_network(labels = lbl, partition = part, attributes = list())
}

#' Décompte attendu de groupes pour un intervalle [from, to)
#' @keywords internal
expected_groups_in_range <- function(part, from, to) {
  sz <- as.integer(table(part))               # tailles par groupe
  sum(sz >= from & sz < to)                   # filtre demi-ouvert
}

#' Normaliser la signature groups(...) en (from,to) numériques
#' @keywords internal
normalize_groups_signature <- function(args = list()) {
  if (length(args) == 0L) {
    return(list(from = 1, to = Inf, text = "b2degrange(from=1,to=Inf)"))
  }
  nm <- names(args)
  if (length(args) == 1L && (is.null(nm) || isTRUE(nm[1L] == ""))) {
    k <- as.numeric(args[[1L]])
    stopifnot(length(k) == 1L, is.finite(k))
    return(list(from = k, to = k + 1, text = sprintf("b2degrange(from=%g,to=%g)", k, k + 1)))
  }
  if (!is.null(nm) && all(c("from","to") %in% nm)) {
    f <- as.numeric(args[["from"]]); t <- as.numeric(args[["to"]])
    stopifnot(length(f) == 1L, length(t) == 1L, is.finite(f) || is.infinite(f), is.finite(t) || is.infinite(t))
    txt <- sprintf("b2degrange(from=%s,to=%s)",
                   if (is.infinite(f)) "Inf" else as.character(f),
                   if (is.infinite(t)) "Inf" else as.character(t))
    return(list(from = f, to = t, text = txt))
  }
  stop("Signature groups(...) inconnue pour ce test.")
}

#' Vérifier la traduction ERPM → ERGM dans l'appel erpm(..., eval_call=FALSE)
#' @keywords internal
check_translation_contains <- function(call_ergm, expected) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)
  grepl(expected, compact, fixed = TRUE)
}

# ======================================================================================
# Noyau de test
# ======================================================================================

#' Exécuter un sous-test pour une partition et une signature `groups(...)`
#' @keywords internal
run_one_case <- function(part, name, call_txt, args) {
  nw <- make_nw_from_partition(part)                         # réseau biparti
  ng <- normalize_groups_signature(args)                     # from/to attendus
  stat_val <- as.integer(summary(as.formula(nw ~ b2degrange(from = ng$from, to = ng$to))))
  truth    <- expected_groups_in_range(part, ng$from, ng$to)

  f <- as.formula(paste0("nw ~ ", call_txt))                 # formule ERPM
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  call_ergm <- erpm(f, eval_call = FALSE, verbose = TRUE)    # montre aussi les calls

  ok_trad <- check_translation_contains(call_ergm, gsub("\\s+", "", ng$text))

  cat(sprintf("\n[CAS %-18s] part={%s}", name, paste(part, collapse=",")))
  cat(sprintf("\t  attendu(groups)  : %s", ng$text))
  cat(sprintf("\t  summary(.)       : %d", stat_val))
  cat(sprintf("\t  attendu(part)    : %d", truth))
  cat(sprintf("\t  traduction erpm  : %s\n", if (ok_trad) "OK" else "KO"))

  list(
    ok_stat = identical(stat_val, truth),
    ok_trad = ok_trad,
    stat = stat_val,
    expected = truth,
    call_line = paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  )
}

#' Exécuter tous les cas pour une partition donnée
#' @keywords internal
run_cases_for_partition <- function(part, panel) {
  res <- lapply(panel, function(cx) {
    out <- run_one_case(part, cx$name, cx$call_txt, cx$args)
    data.frame(
      case     = cx$name,
      ok_stat  = out$ok_stat,
      ok_trad  = out$ok_trad,
      stat     = out$stat,
      expected = out$expected,
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

cases <- list(
  list(name = "groups_all",         call_txt = "groups",                args = list()),
  list(name = "groups_exact_1",     call_txt = "groups(1)",             args = list(1)),
  list(name = "groups_exact_2",     call_txt = "groups(2)",             args = list(2)),
  list(name = "groups_exact_3",     call_txt = "groups(3)",             args = list(3)),
  list(name = "groups_range_2_4",   call_txt = "groups(from=2,to=4)",   args = list(from = 2, to = 4)),
  list(name = "groups_range_1_Inf", call_txt = "groups(from=1,to=Inf)", args = list(from = 1, to = Inf))
)

# ======================================================================================
# Run principal
# ======================================================================================

#' Lancer l'ensemble des tests et produire un bilan
run_all_tests <- function() {
  # --- Logging console + fichier (dans une fonction => on.exit marche) ---
  log_path <- file.path("scripts","test","test_groups.log")
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(log_path)) unlink(log_path)

  sink(file = log_path, split = TRUE)          # stdout -> fichier + console
  con_msg <- file(log_path, open = "at")       # messages -> même fichier
  sink(con_msg, type = "message")

  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(close(con_msg),        silent = TRUE)
    try(sink(),                silent = TRUE)
    flush.console()
  }, add = TRUE)

  # --- Tests ---
  set.seed(42)
  cat("=== TEST ERPM: groups → b2degrange ===\n")

  all_results <- list(); total_ok <- 0L; total_n <- 0L
  for (nm in names(partitions)) {
    cat(sprintf("\n--- Partition %s ---\n", nm))
    df <- run_cases_for_partition(partitions[[nm]], cases)
    all_results[[nm]] <- df
    total_ok <- total_ok + sum(df$ok_stat) + sum(df$ok_trad)
    total_n  <- total_n  + 2L * nrow(df)
  }

  cat(sprintf("\n=== Bilan global : %d / %d validations OK ===\n", total_ok, total_n))
  if (total_ok < total_n) stop(sprintf("Echec: %d validations KO", total_n - total_ok))
  invisible(all_results)
}

if (identical(environment(), globalenv())) run_all_tests()