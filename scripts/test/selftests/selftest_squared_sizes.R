# ======================================================================================
# Fichier : scripts/test/selftests/selftest_squared_sizes.R
# Objet   : Self-test autonome pour l'effet ERPM `squared_sizes`
# Auteur  : Jérémie Chichignoud - Cub'itech
#
# But du fichier
#   - PHASE 1 (SUMMARY) : valider la statistique via summary(nw ~ squared_sizes(...)).
#   - PHASE 2 (ERPM FIT): valider que erpm() construit un modèle et renvoie des coefs finis.
#   - PHASE 3 (MCMC)    : diagnostic "multi-toggle" (le seul truc qui nous intéresse quand
#                         on active le debug côté C et qu'on ne veut pas se faire noyer).
#
# Important
#   - Les phases 1/2 peuvent spammer la console (print de dataframes + summaries).
#   - Pour bosser proprement sur la phase 3, on peut désactiver 1/2 via des flags.
#   - Ici on structure volontairement le run pour que commenter / activer soit trivial.
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
  } else {
    message("[ergm_patch] fonction ergm_patch_enable absente")
  }
} else {
  message("[ergm_patch] scripts/ergm_patch.R introuvable, on continue sans patch")
}

# --------------------------------------------------------------------------------------
# Chargement package/terme + wrapper ERPM
# --------------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else {
  stop("Exécuter depuis la racine du package (DESCRIPTION) avec devtools disponible.")
}

if (!exists("InitErgmTerm.squared_sizes", mode = "function")) {
  stop("InitErgmTerm.squared_sizes introuvable après load_all().")
}

# Pour certains setups, le wrapper n'est pas attaché automatiquement.
# On tente un source direct en fallback (sans magie).
if (!exists("build_bipartite_from_inputs", mode = "function") || !exists("erpm", mode = "function")) {
  if (file.exists("R/erpm.R")) {
    source("R/erpm.R", local = FALSE)
  }
}

if (!exists("build_bipartite_from_inputs", mode = "function")) {
  stop("build_bipartite_from_inputs indisponible. Charger le wrapper ERPM.")
}

# ======================================================================================
# Réglages de run (le point clé du fichier)
# ======================================================================================
# Objectif: pouvoir isoler la PHASE 3 sans éditer 50 endroits.
RUN <- list(
  phase1_summary = TRUE,   # TRUE = on valide summary() ; FALSE = on saute
  phase2_fit     = TRUE,   # TRUE = on fit un mini-modèle ; FALSE = on saute
  phase3_mcmc    = TRUE,   # TRUE = on lance le probe multi-toggle ; FALSE = on saute

  # "quiet" réduit la pollution console des phases 1/2 sans les supprimer.
  # (PHASE 3, elle, doit rester verbeuse pour voir passer les traces debug.)
  quiet_phase1   = FALSE,
  quiet_phase2   = FALSE
)

# ======================================================================================
# Fonctions locales — utilitaires de référence
# ======================================================================================

# Construction du biparti via le wrapper (exigé)
.build_bipartite_nw_via_wrapper <- function(part) {
  stopifnot(length(part) >= 1L)

  nodes <- data.frame(
    label = utils::head(LETTERS, length(part)),
    stringsAsFactors = FALSE
  )

  builder <- get("build_bipartite_from_inputs", mode = "function")

  # Plusieurs signatures possibles selon les versions du wrapper.
  candidates <- list(
    quote(builder(partition = part, nodes = nodes)),
    quote(builder(partition = part, labels = nodes$label, attributes = list())),
    quote(builder(partition = part)),
    quote(builder(partition = part, labels = nodes$label))
  )

  last_err <- NULL
  for (expr in candidates) {
    out <- try(eval(expr), silent = TRUE)

    if (inherits(out, "try-error")) {
      last_err <- out
      next
    }

    if (inherits(out, "network")) return(out)

    if (is.list(out)) {
      for (nm in c("network","nw","net","graph","g","bip")) {
        if (!is.null(out[[nm]]) && inherits(out[[nm]], "network")) {
          return(out[[nm]])
        }
      }
    }
  }

  fm <- try(formals(builder), silent = TRUE)
  stop(paste0(
    "build_bipartite_from_inputs a échoué.\n",
    "Formals: ", if (!inherits(fm, "try-error")) paste(names(fm), collapse = ", ") else "<inconnus>", "\n",
    "Dernière erreur: ", if (!is.null(last_err)) as.character(last_err)[1] else "<aucune>"
  ))
}

# Tailles de groupes depuis partition
.group_sizes_from_partition <- function(part) as.integer(table(part))

# Attendue: sum_{g : size(g) ∈ sizes} size(g)^pow (scalaire agrégé)
.expected_summary_squared_sizes_from_partition <- function(part, sizes = NULL, pow = 2) {
  sz <- .group_sizes_from_partition(part)

  if (length(pow) != 1L) stop("'pow' doit être scalaire (longueur 1).")
  pow <- as.numeric(pow)

  n1 <- length(part)
  sizes_eff <- if (is.null(sizes) || length(sizes) == 0L) seq_len(as.integer(n1)) else as.integer(sizes)

  idx <- sz %in% sizes_eff
  if (!any(idx)) return(0)

  sum(sz[idx]^pow)
}

# Contrôle identité (uniquement pour pow=2, toutes tailles)
#   sum_g deg(g)^2 = m + 2 * sum_g C(deg(g), 2)
# où m = nombre d'arêtes (membership edges).
.identity_pow2_all_summary_on_network <- function(nw) {
  stopifnot(inherits(nw, "network"))

  n1 <- network::get.network.attribute(nw, "bipartite")
  if (is.null(n1) || is.na(n1)) stop("Réseau non biparti (attribut 'bipartite' manquant).")

  n  <- network::network.size(nw)
  v2 <- seq.int(as.integer(n1) + 1L, n)

  deg2 <- vapply(
    v2,
    function(v) length(network::get.neighborhood(nw, v, type = "combined")),
    integer(1L)
  )

  m <- network::network.edgecount(nw)
  m + 2L * sum(choose(deg2, 2))
}

# Normalisation simple (sizes/pow) + texte lisible pour affichage
.normalize_squared_sizes_args <- function(args = list()) {
  out <- list(sizes = NULL, pow = 2)

  if (length(args)) {
    nm <- names(args)
    if (!is.null(nm) && length(nm)) {
      if ("sizes" %in% nm) out$sizes <- args[["sizes"]]
      if ("pow"   %in% nm) out$pow   <- args[["pow"]]
    }
  }

  pieces <- character(0)
  if (!is.null(out$sizes)) pieces <- c(pieces, sprintf("sizes=c(%s)", paste(out$sizes, collapse = ",")))
  if (!identical(as.numeric(out$pow), 2)) pieces <- c(pieces, sprintf("pow=%s", out$pow))

  out$text <- if (length(pieces)) sprintf("squared_sizes(%s)", paste(pieces, collapse = ",")) else "squared_sizes"
  out
}

# Vérification légère: la traduction erpm(...) contient le terme et les args attendus
.check_translation_ok_erpm <- function(call_ergm, term = "squared_sizes", args = list()) {
  line    <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)

  if (!grepl(paste0("\\b", term, "\\("), compact)) return(FALSE)

  ok <- TRUE
  if (!is.null(args$sizes)) ok <- ok && grepl("sizes=", compact, fixed = TRUE)
  if (!is.null(args$pow) && !identical(as.numeric(args$pow), 2)) ok <- ok && grepl("pow=", compact, fixed = TRUE)

  ok
}

# Petit helper d'affichage: on évite les sorties inutiles si "quiet" est actif.
.maybe_print <- function(x, quiet = FALSE) {
  if (!isTRUE(quiet)) print(x)
  invisible(NULL)
}

# ======================================================================================
# PHASE 1 — SUMMARY : summary(.) vs référence partition
# ======================================================================================

.run_one_case_summary_squared_sizes <- function(part, case, quiet = FALSE) {
  nw <- .build_bipartite_nw_via_wrapper(part)

  sg <- .normalize_squared_sizes_args(case$args)

  expected_val <- .expected_summary_squared_sizes_from_partition(
    part,
    sizes = sg$sizes,
    pow   = sg$pow
  )

  # summary(nw ~ squared_sizes(...))
  f <- as.formula(paste0("nw ~ ", case$call_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  stat_val <- as.numeric(summary(f))

  # Traduction erpm() (optionnel)
  ok_trad <- NA
  if (exists("erpm", mode = "function")) {
    call_ergm <- erpm(f, eval.call = FALSE, verbose = FALSE)
    ok_trad <- .check_translation_ok_erpm(call_ergm, term = "squared_sizes", args = case$args)
  }

  # Identité (optionnel): uniquement pow=2 et tailles par défaut
  ok_identity <- NA
  if (isTRUE(as.numeric(sg$pow) == 2) && (is.null(sg$sizes) || length(sg$sizes) == 0L)) {
    id_val <- .identity_pow2_all_summary_on_network(nw)
    ok_identity <- identical(unname(as.integer(stat_val)), unname(as.integer(id_val)))
  }

  out <- list(
    case       = case$name,
    signature  = sg$text,
    stat       = stat_val,
    expected   = expected_val,
    ok_stat    = identical(unname(as.integer(stat_val)), unname(as.integer(expected_val))),
    ok_trad    = ok_trad,
    ok_ident   = ok_identity
  )

  if (!isTRUE(quiet)) {
    # Affichage compact: la ligne résume tout, inutile de spammer un gros dump.
    cat(sprintf("  - %-14s | %-28s | stat=%s expected=%s | ok_stat=%s\n",
                out$case, out$signature, out$stat, out$expected, out$ok_stat))
  }

  out
}

.run_cases_for_partition_summary_squared_sizes <- function(part, cases, quiet = FALSE) {
  out <- lapply(cases, function(cs) .run_one_case_summary_squared_sizes(part, cs, quiet = quiet))

  data.frame(
    case       = vapply(out, `[[`, character(1), "case"),
    signature  = vapply(out, `[[`, character(1), "signature"),
    ok_stat    = vapply(out, `[[`, logical(1),  "ok_stat"),
    ok_trad    = vapply(out, function(x) if (is.na(x$ok_trad))  NA else isTRUE(x$ok_trad),  logical(1)),
    ok_ident   = vapply(out, function(x) if (is.na(x$ok_ident)) NA else isTRUE(x$ok_ident), logical(1)),
    stat       = vapply(out, function(x) as.numeric(x$stat),     numeric(1)),
    expected   = vapply(out, function(x) as.numeric(x$expected), numeric(1)),
    stringsAsFactors = FALSE
  )
}

.run_phase1_summary <- function(partitions, cases, quiet = FALSE) {
  cat("\n=== PHASE 1: SUMMARY (statistique) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  all_summary <- list()
  n_checks <- 0L
  n_ok     <- 0L

  for (nm in names(partitions)) {
    cat("\n--- Partition ", nm, " ---\n", sep = "")
    df <- .run_cases_for_partition_summary_squared_sizes(partitions[[nm]], cases, quiet = quiet)

    # En mode normal, on affiche aussi le tableau complet (utile quand ça casse).
    .maybe_print(df, quiet = quiet)

    all_summary[[nm]] <- df

    # ok_stat toujours compté
    n_checks <- n_checks + sum(!is.na(df$ok_stat))
    n_ok     <- n_ok     + sum(df$ok_stat, na.rm = TRUE)

    # ok_ident et ok_trad sont optionnels (NA possible)
    n_checks <- n_checks + sum(!is.na(df$ok_ident))
    n_ok     <- n_ok     + sum(df$ok_ident, na.rm = TRUE)

    n_checks <- n_checks + sum(!is.na(df$ok_trad))
    n_ok     <- n_ok     + sum(df$ok_trad, na.rm = TRUE)
  }

  cat(sprintf("\nBilan SUMMARY: %d / %d validations OK\n", n_ok, n_checks))
  if (n_ok < n_checks) stop(sprintf("Echec SUMMARY: %d validations KO", n_checks - n_ok))

  invisible(all_summary)
}

# ======================================================================================
# PHASE 2 — ERPM FIT : vérifier qu'un fit passe et renvoie un coef fini
# ======================================================================================

.run_one_case_erpm_fit_squared_sizes <- function(part, name, call_txt, ctrl, quiet = FALSE) {
  if (!exists("erpm", mode = "function")) {
    if (!isTRUE(quiet)) cat(sprintf("  - %-16s | erpm() indisponible -> SKIP\n", name))
    return(list(name = name, ok = NA, coef = NA, fit = NULL))
  }

  f <- as.formula(paste0("part ~ ", call_txt))
  environment(f) <- list2env(list(part = part), parent = parent.frame())

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

.run_phase2_fit <- function(part_ref, ctrl_fit, quiet = FALSE) {
  cat("\n=== PHASE 2: ERPM FIT (sanity check) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  fits <- list(
    FIT_sq_all_pow2  = .run_one_case_erpm_fit_squared_sizes(part_ref, "sq_all_pow2",  "squared_sizes",            ctrl_fit, quiet = quiet),
    FIT_sq_2to4_pow2 = .run_one_case_erpm_fit_squared_sizes(part_ref, "sq_2to4_pow2", "squared_sizes(sizes=2:4)", ctrl_fit, quiet = quiet),
    FIT_sq_all_pow3  = .run_one_case_erpm_fit_squared_sizes(part_ref, "sq_all_pow3",  "squared_sizes(pow=3)",     ctrl_fit, quiet = quiet)
  )

  ok_fit <- vapply(fits, function(x) if (is.na(x$ok)) NA else isTRUE(x$ok), logical(1))
  n_fit_ok  <- sum(ok_fit, na.rm = TRUE)
  n_fit_tot <- sum(!is.na(ok_fit))

  cat(sprintf("\nBilan FIT: %d / %d OK\n", n_fit_ok, n_fit_tot))
  if (n_fit_ok < n_fit_tot) stop(sprintf("Echec FIT: %d KO", n_fit_tot - n_fit_ok))

  # En mode normal, on affiche les summaries (c'est ce qui pollue le plus).
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
  #   - observer dans la console les traces debug du changestat
  #     (ex: "[squared_sizes] MULTI-TOGGLE ntoggles=...")
  #
  # Important: on garde verbose=TRUE exprès, sinon on ne voit rien.

  ctrl <- control.simulate.formula(
    MCMC.burnin   = 1000,
    MCMC.interval = 1,
    MCMC.prop     = ~ sparse
  )

  f <- nw ~ squared_sizes

  sim <- simulate(
    f,
    nsim        = 1,
    control     = ctrl,
    verbose     = TRUE
  )

  print(sim)
  invisible(sim)
}

.run_phase3_mcmc <- function(part_probe) {
  cat("\n=== PHASE 3: MCMC MULTI-TOGGLE PROBE ===\n")
  cat("Objectif: voir passer les traces debug du changestat (multi-toggle).\n")
  cat("Si tu ne vois rien, c'est souvent que le define debug n'est pas actif côté C,\n")
  cat("ou que l'update ne déclenche pas le chemin multi-toggle.\n\n")

  nw_probe <- .build_bipartite_nw_via_wrapper(part_probe)
  .run_mcmc_multitoggle_probe(nw_probe)

  cat("\nSi '[squared_sizes] MULTI-TOGGLE ntoggles=...' apparaît, test OK.\n")
  invisible(TRUE)
}

# ======================================================================================
# Jeu de données + cas test
# ======================================================================================

partitions <- list(
  P1 = c(1, 2, 2, 3, 3, 3),
  P2 = c(1, 1, 2, 3, 3, 4, 4, 4),
  P3 = c(1, 1, 1, 2, 2, 3),
  P4 = c(1, 2, 3, 4, 5),
  P5 = rep(1, 6)
)

cases <- list(
  list(
    name     = "sq_all_pow2",
    call_txt = "squared_sizes",
    args     = list(sizes = NULL, pow = 2)
  ),
  list(
    name     = "sq_2to4_pow2",
    call_txt = "squared_sizes(sizes=c(2,3,4))",
    args     = list(sizes = c(2,3,4), pow = 2)
  ),
  list(
    name     = "sq_all_pow3",
    call_txt = "squared_sizes(pow=3)",
    args     = list(sizes = NULL, pow = 3)
  ),
  list(
    name     = "sq_1to2_pow2",
    call_txt = "squared_sizes(sizes=1:2)",
    args     = list(sizes = 1:2, pow = 2)
  ),
  list(
    name     = "sq_3_pow2",
    call_txt = "squared_sizes(sizes=3)",
    args     = list(sizes = 3, pow = 2)
  )
)

# Contrôle MCMLE commun (fit)
ctrl_fit <- control.ergm(
  MCMLE.maxit     = 10,
  MCMC.samplesize = 1e4
)

# ======================================================================================
# Run principal
# ======================================================================================

run_all_tests_squared_sizes <- function() {
  set.seed(42)

  cat("=== SELFTEST ERPM: squared_sizes ===\n")
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
      part_ref  = partitions$P1,
      ctrl_fit  = ctrl_fit,
      quiet     = isTRUE(RUN$quiet_phase2)
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
  run_all_tests_squared_sizes()
}

# --------------------------------------------------------------------------------------
# Fin de script: on désactive le patch si présent (histoire de ne pas laisser traîner)
# --------------------------------------------------------------------------------------
if (exists("ergm_patch_disable", mode = "function")) {
  ergm_patch_disable()
}