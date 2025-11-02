# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cliques.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cliques`
# Exécution: Rscript scripts/test/selftests/selftest_cliques.R
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

# --------------------------------------------------------------------------------------
# Chargements du package et des utilitaires ERPM
# --------------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else {
  # On peut fonctionner si le terme est déjà installé/chargé via library(ERPM)
  if (!exists("InitErgmTerm.cliques", mode = "function")) {
    stop("InitErgmTerm.cliques introuvable. Exécute depuis le package (devtools::load_all) ou charge le package.")
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
  }
  # sinon, on vivra sans vérification de traduction
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

# Tailles de groupes
.group_sizes_from_partition <- function(part) as.integer(table(part))

# N1 = #acteurs (taille du mode 1 dans le biparti)
.get_N1_from_nw <- function(nw) {
  n1 <- network::get.network.attribute(nw, "bipartite")
  if (is.null(n1) || is.na(n1)) stop("Réseau non biparti ou attribut 'bipartite' manquant.")
  as.integer(n1)
}

# Valeur attendue de la stat cliques (partition → valeur théorique)
# - k == 1 : nombre de groupes de taille 1
# - k >= 2 : sum_g choose(n_g, k)
# - normalized : division par choose(N1, k) (k=1 ⇒ denom = N1)
.expected_cliques <- function(part, k = 2L, normalized = FALSE) {
  sz <- .group_sizes_from_partition(part)
  if (k == 1L) {
    num <- sum(sz == 1L)
  } else {
    num <- sum(choose(sz, k))
  }
  if (!normalized) return(num)
  N1 <- length(part) # Dans notre construction, N1 = #acteurs = longueur de la partition
  denom <- choose(N1, k)
  if (denom == 0 || is.na(denom)) return(0)
  num / denom
}

# Normaliser la "signature" attendue côté traduction erpm
# Défauts ERPM/ERGM : k = 2, normalized = FALSE
.normalize_cliques_signature <- function(args = list()) {
  out <- list(k = 2L, normalized = FALSE, text = "cliques(k=2,normalized=FALSE)")
  if (length(args)) {
    nm <- names(args)
    if (!is.null(nm) && length(nm)) {
      if ("k"          %in% nm) out$k          <- as.integer(args[["k"]])
      if ("clique_size"%in% nm) out$k          <- as.integer(args[["clique_size"]])
      if ("normalized" %in% nm) out$normalized <- isTRUE(args[["normalized"]])
    }
  }
  out$text <- sprintf("cliques(k=%s,normalized=%s)", as.integer(out$k), if (out$normalized) "TRUE" else "FALSE")
  out
}

# Vérifier la présence des seuls arguments non-défaut dans l'appel ergm traduit
.check_translation_ok <- function(call_ergm, args = list()) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)

  # Par défaut on s'attend à "cliques(" (pas b2degrange), y compris k=1 (supporté côté C/R)
  if (!grepl("\\bcliques\\(", compact)) return(FALSE)

  # Défauts
  def_k   <- 2L
  def_nrm <- FALSE

  checks <- logical(0)

  # k présent seulement si ≠ défaut
  k_eff <- if (!is.null(args$clique_size)) as.integer(args$clique_size)
           else if (!is.null(args$k))      as.integer(args$k)
           else def_k
  if (!identical(k_eff, def_k)) {
    # On accepte soit cliques(k=...), soit cliques(clique_size=...)
    checks <- c(checks, grepl(paste0("k=", k_eff), compact, fixed = TRUE) |
                       grepl(paste0("clique_size=", k_eff), compact, fixed = TRUE))
  }

  # normalized présent seulement si TRUE (défaut = FALSE)
  nrm_eff <- isTRUE(args$normalized)
  if (isTRUE(nrm_eff != def_nrm)) {
    checks <- c(checks, grepl("normalized=TRUE", compact, fixed = TRUE))
  }

  # si aucun arg non-défaut → OK dès que le terme est présent
  if (!length(checks)) return(TRUE)
  all(checks)
}

# Exécuter un cas
.run_one_case <- function(part, name, call_txt, args) {
  nw <- .make_nw_from_partition(part)
  sig <- .normalize_cliques_signature(args)

  truth <- .expected_cliques(part, k = sig$k, normalized = sig$normalized)

  f <- as.formula(paste0("nw ~ ", call_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  stat_val <- as.numeric(summary(f))

  ok_trad <- NA
  if (exists("erpm", mode = "function")) {
    call_ergm <- erpm(f, eval_call = FALSE, verbose = TRUE)
    ok_trad <- .check_translation_ok(call_ergm, args = args)
  }

  cat(sprintf("\n[CAS %-18s] part={%s}", name, paste(part, collapse=",")))
  cat(sprintf("\t  appel          : %s", call_txt))
  cat(sprintf("\t  summary(.)     : %s", format(stat_val)))
  cat(sprintf("\t  attendu(part)  : %s", format(truth)))
  if (!is.na(ok_trad))     cat(sprintf("\t  traduction erpm : %s", if (ok_trad) "OK" else "KO"))
  cat("\n")

  list(
    # ok_stat   = isTRUE(all.equal(unname(as.numeric(stat_val)), unname(as.numeric(truth)))),
    ok_stat = isTRUE(all.equal(unname(as.numeric(stat_val)),
                           unname(as.numeric(truth)),
                           tolerance = 1e-10)),
    ok_trad   = ok_trad,
    stat      = stat_val,
    expected  = truth
  )
}

# Panel sur une partition
.run_cases_for_partition <- function(part, panel) {
  res <- lapply(panel, function(cx) {
    out <- .run_one_case(part, cx$name, cx$call_txt, cx$args)
    data.frame(
      case     = cx$name,
      ok_stat  = out$ok_stat,
      ok_trad  = if (is.na(out$ok_trad)) NA else out$ok_trad,
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
  P1 = c(1, 2, 2, 3, 3, 3),            # tailles = 1,2,3
  P2 = c(1, 1, 2, 3, 3, 4, 4, 4),      # tailles = 2,1,2,1,3  -> (après comptage: 2,1,2,1,3) -> tri: 1,1,2,2,3
  P3 = c(1, 1, 1, 2, 2, 3),            # tailles = 3,2,1
  P4 = c(1, 2, 3, 4, 5),               # tailles = 1,1,1,1,1
  P5 = rep(1, 6)                       # taille = 6 (un seul groupe)
)

# Quelques cas ciblés :
# - k=1 (non-normalisé / normalisé)
# - k=2 (défaut / explicite / normalisé)
# - k=3 (non-normalisé)
cases <- list(
  list(name="k1_raw",     call_txt="cliques(k=1)",                      args=list(k=1, normalized=FALSE)),
  list(name="k1_norm",    call_txt="cliques(k=1, normalized=TRUE)",     args=list(k=1, normalized=TRUE)),
  list(name="k2_default", call_txt="cliques()",                         args=list(k=2, normalized=FALSE)),
  list(name="k2_exp",     call_txt="cliques(k=2)",                      args=list(k=2, normalized=FALSE)),
  list(name="k2_norm",    call_txt="cliques(k=2, normalized=TRUE)",     args=list(k=2, normalized=TRUE)),
  list(name="k3_raw",     call_txt="cliques(k=3)",                      args=list(k=3, normalized=FALSE)),
  # variante avec alias clique_size=
  list(name="k3_alias",   call_txt="cliques(clique_size=3)",            args=list(clique_size=3, normalized=FALSE))
)

# ======================================================================================
# Run principal
# ======================================================================================
run_all_tests <- function() {
  # Logging aussi vers scripts/test/selftests/selftest_cliques.log
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
  cat("=== TEST ERPM: cliques ===\n")

  all_results <- list(); total_ok <- 0L; total_n <- 0L
  for (nm in names(partitions)) {
    cat(sprintf("\n--- Partition %s ---\n", nm))
    df <- .run_cases_for_partition(partitions[[nm]], cases)
    all_results[[nm]] <- df

    # ok_stat : 1 validation / cas
    total_ok <- total_ok + sum(df$ok_stat, na.rm = TRUE)
    total_n  <- total_n  + sum(!is.na(df$ok_stat))

    # traduction (si erpm dispo)
    total_ok <- total_ok + sum(df$ok_trad, na.rm = TRUE)
    total_n  <- total_n  + sum(!is.na(df$ok_trad))

    print(df)
  }

  cat(sprintf("\n=== Bilan global : %d / %d validations OK ===\n", total_ok, total_n))
  if (total_ok < total_n) stop(sprintf("Echec: %d validations KO", total_n - total_ok))
  invisible(all_results)
}

if (identical(environment(), globalenv())) run_all_tests()