# ==============================================================================
# Fichier : scripts/test/selftests/selftest_cov_ingroup.R
# Objet   : Self-test autonome pour l’effet ERPM/ERGM `cov_ingroup` (multi-toggle ready)
# Exécution: Rscript scripts/test/selftests/selftest_cov_ingroup.R
#
# But du fichier
#   - PHASE 1 (SUMMARY) : valider la statistique via summary(nw ~ cov_ingroup(...))
#                         contre une référence pure R (partition + covariate).
#   - PHASE 2 (ERPM FIT): valider que erpm() construit un modèle et renvoie des coefs finis.
#   - PHASE 3 (MCMC)    : diagnostic "multi-toggle" (observer les traces debug C si activées).
#
# Notes importantes (multi-toggle)
#   - Le fait que le changestat soit D_ (multi-toggle) doit être visible même
#     si ergm propose parfois ntoggles=1. Le point ici: aucune hypothèse "1 toggle".
#   - Pour voir des traces côté C:
#       * activer DEBUG_COV_INGROUP=1 dans changestat_cov_ingroup.c
#       * recompiler le package
#     Puis lancer PHASE 3 (simulate) en verbose=TRUE.
#
#   - PHASE 3 n’est pas un test "statistique" : c’est un probe pour déclencher
#     la MCMC et voir passer les logs multitoggle.
# ==============================================================================

options(ergm.loglik.warn_dyads = FALSE)
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",    quietly = TRUE)) stop("Package 'ergm' requis.")
})

suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

# Patch ERGM si disponible
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
}

# Charger package
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Exécuter depuis la racine du package (DESCRIPTION) avec devtools disponible.")
}

# Wrapper ERPM
if (!exists("erpm", mode = "function")) {
  if (file.exists("R/erpm_wrapper.R")) {
    source("R/erpm_wrapper.R", local = FALSE)
  } else {
    stop("R/erpm_wrapper.R introuvable.")
  }
}

if (!exists("build_bipartite_from_inputs", mode = "function")) {
  stop("build_bipartite_from_inputs() manquant.")
}

cat("=== SELFTEST ERPM: cov_ingroup (multi-toggle) ===\n")
cat("R:", paste(R.version$major, R.version$minor, sep="."), "\n")
cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")

# ==============================================================================
# Réglages de run
# ==============================================================================
RUN <- list(
  phase1_summary = FALSE,
  phase2_fit     = TRUE,
  phase3_mcmc    = FALSE,

  quiet_phase1   = FALSE,
  quiet_phase2   = FALSE
)

# ==============================================================================
# Données
# ==============================================================================
partitions <- list(
  P1 = c(1,2,3,3,3,4),
  P2 = c(1, 2,2, 3,3,3, 4, 5,5,5, 6,6,6,6,6),
  P3 = c(1,1,2,2,3)
)

.make_nodes <- function(part) {
  n <- length(part)
  set.seed(100 + n)
  data.frame(
    label  = paste0("N", seq_len(n)),
    age    = sample(20:60, n, TRUE),
    score  = round(runif(n, 0, 10), 1),
    gender = sample(c("F","M"), n, TRUE),
    dept   = sample(c("A","B","C"), n, TRUE, prob = c(0.5, 0.3, 0.2)),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Helpers
# ==============================================================================

# Construction du biparti via le wrapper (exigé)
.make_nw <- function(part, nodes) {
  built <- build_bipartite_from_inputs(partition = part, nodes = nodes)
  if (is.list(built) && !is.null(built$network) && inherits(built$network, "network")) return(built$network)
  if (inherits(built, "network")) return(built)
  stop("build_bipartite_from_inputs() ne renvoie pas un 'network'.")
}

# Extraction x (numérique) pour la référence pure R
# - si cov est numérique: on prend nodes[[cov]]
# - si cov est catégoriel + category: on renvoie un indicateur 0/1
.get_x_ref <- function(nodes, cov, category = NULL) {
  if (!cov %in% names(nodes)) stop("Covariate introuvable dans nodes: ", cov)
  v <- nodes[[cov]]

  if (!is.null(category)) {
    xb <- as.integer(as.character(v) == category)
    xb[is.na(xb)] <- 0L
    return(as.double(xb))
  }

  x <- suppressWarnings(as.numeric(v))
  if (any(!is.finite(x))) stop("Covariate non numérique / contient NA/Inf: ", cov)
  as.double(x)
}

# Référence pure R:
#   T = sum_g 1[n_g ∈ S] * n_g * sum_{i in g} x_i
.expected_cov_ingroup <- function(part, nodes, cov, size = NULL, category = NULL) {
  x <- .get_x_ref(nodes, cov = cov, category = category)

  # groupes par partition
  groups <- split(seq_along(part), part)
  sizes_g <- vapply(groups, length, integer(1))

  S <- if (is.null(size) || length(size) == 0L) NULL else as.integer(size)

  tot <- 0.0
  for (g in seq_along(groups)) {
    ng <- sizes_g[g]
    if (!is.null(S) && !(ng %in% S)) next
    tot <- tot + ng * sum(x[groups[[g]]])
  }
  tot
}

# Build formula nw ~ rhs
.make_f_nw <- function(nw, rhs) {
  f <- as.formula(paste0("nw ~ ", rhs))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  f
}

# Build formula partition ~ rhs (for erpm wrapper)
.make_f_part <- function(part, nodes, rhs) {
  partition <- part
  f <- as.formula(paste0("partition ~ ", rhs))
  environment(f) <- list2env(list(partition = partition, nodes = nodes), parent = parent.frame())
  f
}

# SUMMARY sur network (ergm) + b1part
.summary_network <- function(part, nodes, rhs) {
  nw <- .make_nw(part, nodes)
  f  <- .make_f_nw(nw, rhs)
  as.numeric(suppressMessages(summary(f, constraints = ~ b1part)))
}

# SUMMARY via erpm() (traduction -> réseau)
.summary_erpm_translation <- function(part, nodes, rhs) {
  f <- .make_f_part(part, nodes, rhs)

  call_ergm <- erpm(f, eval.call = FALSE, verbose = FALSE, nodes = nodes)
  ergm_form <- call_ergm[[2L]]
  rhs_e     <- ergm_form[[3L]]

  nw2 <- .make_nw(part, nodes)
  f2  <- as.formula(bquote(nw2 ~ .(rhs_e)))
  environment(f2) <- list2env(list(nw2 = nw2), parent = parent.frame())

  as.numeric(suppressMessages(summary(f2, constraints = ~ b1part)))
}

.maybe_print <- function(x, quiet = FALSE) {
  if (!isTRUE(quiet)) print(x)
  invisible(NULL)
}

# Vérification simple de la "traduction" erpm(...) (présence du terme)
.check_translation_ok <- function(call_ergm, term = "cov_ingroup") {
  line    <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)
  grepl(paste0("\\b", term, "\\("), compact)
}

# ==============================================================================
# Cas à tester
# ==============================================================================
cases <- list(
  list(
    name = "age_all",
    rhs  = "cov_ingroup('age')",
    cov  = "age",
    size = NULL,
    category = NULL
  ),
  list(
    name = "age_2to4",
    rhs  = "cov_ingroup('age', size=2:4)",
    cov  = "age",
    size = 2:4,
    category = NULL
  ),
  list(
    name = "score_all",
    rhs  = "cov_ingroup('score')",
    cov  = "score",
    size = NULL,
    category = NULL
  ),
  list(
    name = "gender_F",
    rhs  = "cov_ingroup('gender', category='F')",
    cov  = "gender",
    size = NULL,
    category = "F"
  ),
  list(
    name = "dept_A_2to3",
    rhs  = "cov_ingroup('dept', category='A', size=c(2,3))",
    cov  = "dept",
    size = c(2,3),
    category = "A"
  )
)

# ==============================================================================
# PHASE 1 — SUMMARY : réseau vs référence pure R + traduction erpm
# ==============================================================================
.run_phase1_summary <- function(partitions, cases, quiet = FALSE) {
  cat("\n=== PHASE 1: SUMMARY (référence R vs réseau vs erpm-traduction) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  n_checks <- 0L
  n_ok     <- 0L

  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes(part)

    cat(sprintf("\n--- Partition %s --- n=%d ---\n", nm, length(part)))

    rows <- list()

    for (cs in cases) {
      ref <- .expected_cov_ingroup(part, nodes, cov = cs$cov, size = cs$size, category = cs$category)
      s1  <- .summary_network(part, nodes, cs$rhs)
      s2  <- .summary_erpm_translation(part, nodes, cs$rhs)

      ok1 <- isTRUE(all.equal(as.numeric(s1), as.numeric(ref), tolerance = 0))
      ok2 <- isTRUE(all.equal(as.numeric(s2), as.numeric(ref), tolerance = 0))

      # traduction (présence du terme dans le call)
      f_part   <- .make_f_part(part, nodes, cs$rhs)
      call_ergm <- erpm(f_part, eval.call = FALSE, verbose = FALSE, nodes = nodes)
      ok_trad  <- .check_translation_ok(call_ergm, term = "cov_ingroup")

      rows[[cs$name]] <- data.frame(
        case     = cs$name,
        rhs      = cs$rhs,
        ref      = ref,
        net      = s1,
        erpm     = s2,
        ok_net   = ok1,
        ok_erpm  = ok2,
        ok_trad  = ok_trad,
        stringsAsFactors = FALSE
      )

      if (!isTRUE(quiet)) {
        cat(sprintf("  - %-12s | ref=%s | net=%s (ok=%s) | erpm=%s (ok=%s) | trad_ok=%s\n",
                    cs$name, format(ref, digits=6),
                    format(s1, digits=6), ok1,
                    format(s2, digits=6), ok2,
                    ok_trad))
      }

      n_checks <- n_checks + 3L
      n_ok     <- n_ok + as.integer(ok1) + as.integer(ok2) + as.integer(ok_trad)
    }

    df <- do.call(rbind, rows)
    .maybe_print(df, quiet = quiet)
  }

  cat(sprintf("\nBilan SUMMARY: %d / %d validations OK\n", n_ok, n_checks))
  if (n_ok < n_checks) stop(sprintf("Echec SUMMARY: %d validations KO", n_checks - n_ok))
  invisible(TRUE)
}

# ==============================================================================
# PHASE 2 — ERPM FIT : sanity check (coefs finis)
# ==============================================================================
.run_one_fit <- function(part, nodes, rhs, quiet = FALSE) {
  if (!exists("erpm", mode = "function")) {
    if (!isTRUE(quiet)) cat("  - erpm() indisponible -> SKIP\n")
    return(list(ok = NA, fit = NULL, coef = NA))
  }

  f <- .make_f_part(part, nodes, rhs)

  res <- try(erpm(f, eval.loglik = TRUE, verbose = FALSE, nodes = nodes), silent = TRUE)
  if (inherits(res, "try-error")) {
    if (!isTRUE(quiet)) cat("  - FAIL (erpm erreur)\n")
    return(list(ok = FALSE, fit = NULL, coef = NA))
  }

  cf <- try(stats::coef(res), silent = TRUE)
  ok_coef  <- !(inherits(cf, "try-error")) && length(cf) > 0L && all(is.finite(cf))
  ok_class <- inherits(res, "ergm")

  ok <- isTRUE(ok_class && ok_coef)

  if (!isTRUE(quiet)) {
    cat(sprintf("  - RHS=%-45s | ok=%s | coef=%s\n",
                rhs, ok, if (ok_coef) paste(round(cf, 6), collapse = ", ") else "<NA>"))
  }

  list(ok = ok, fit = res, coef = if (ok_coef) cf else NA)
}

.run_phase2_fit <- function(part_ref, cases, quiet = FALSE) {
  cat("\n=== PHASE 2: ERPM FIT (sanity check) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  nodes <- .make_nodes(part_ref)

  # contrôle MCMLE minimal
  ctrl_fit <- control.ergm(
    MCMLE.maxit     = 10,
    MCMC.samplesize = 1e4
  )

  ok_all <- TRUE
  fits <- list()

  for (cs in cases) {
    # on ne fit pas tous les cas catégoriels si ça devient instable: ici on tente tout,
    # et on laisse le "ok" décider.
    out <- .run_one_fit(part_ref, nodes, cs$rhs, quiet = quiet)
    fits[[cs$name]] <- out
    if (!isTRUE(out$ok) && !is.na(out$ok)) ok_all <- FALSE
  }

  n_fit_ok  <- sum(vapply(fits, function(x) isTRUE(x$ok), logical(1)))
  n_fit_tot <- sum(!is.na(vapply(fits, function(x) x$ok, logical(1))))

  cat(sprintf("\nBilan FIT: %d / %d OK\n", n_fit_ok, n_fit_tot))
  if (n_fit_ok < n_fit_tot) stop(sprintf("Echec FIT: %d KO", n_fit_tot - n_fit_ok))

  # en mode normal, on affiche les summaries (ça pollue, donc optionnel)
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

# ==============================================================================
# PHASE 3 — MCMC multi-toggle (diagnostic)
# ==============================================================================
.run_mcmc_multitoggle_probe <- function(nw) {
  # Objectif:
  #   - déclencher le chemin MCMC
  #   - observer dans la console les traces debug du changestat D_
  #
  # NB: le fait d’avoir ntoggles>1 dépend des propositions utilisées par ergm.
  # Le point: le changestat ne doit PAS supposer ntoggles==1.
  ctrl <- control.simulate.formula(
    MCMC.burnin   = 1000,
    MCMC.interval = 1,
    # Proposition "sparse" est souvent suffisante pour déclencher des moves,
    # et peut générer des listes de toggles selon la config/proposal.
    MCMC.prop     = ~ sparse
  )

  f <- nw ~ cov_ingroup("age")

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
  cat("Objectif: déclencher la MCMC et observer les traces debug D_ côté C.\n")
  cat("Si tu ne vois rien, activer DEBUG_COV_INGROUP=1 et recompiler.\n\n")

  nodes <- .make_nodes(part_probe)
  nw    <- .make_nw(part_probe, nodes)
  .run_mcmc_multitoggle_probe(nw)

  cat("\nProbe terminé.\n")
  invisible(TRUE)
}

# ==============================================================================
# Run principal
# ==============================================================================
run_all_tests_cov_ingroup <- function() {
  set.seed(42)

  # PHASE 1
  if (isTRUE(RUN$phase1_summary)) {
    .run_phase1_summary(partitions = partitions, cases = cases, quiet = isTRUE(RUN$quiet_phase1))
  } else {
    cat("\n=== PHASE 1: SUMMARY ===\nSKIP\n")
  }

  # PHASE 2
  if (isTRUE(RUN$phase2_fit)) {
    .run_phase2_fit(part_ref = partitions$P1, cases = cases, quiet = isTRUE(RUN$quiet_phase2))
  } else {
    cat("\n=== PHASE 2: ERPM FIT ===\nSKIP\n")
  }

  # PHASE 3
  if (isTRUE(RUN$phase3_mcmc)) {
    .run_phase3_mcmc(part_probe = partitions$P1)
  } else {
    cat("\n=== PHASE 3: MCMC MULTI-TOGGLE PROBE ===\nSKIP\n")
  }

  invisible(TRUE)
}

if (identical(environment(), globalenv())) {
  run_all_tests_cov_ingroup()
}

# fin: désactiver patch si présent
if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()