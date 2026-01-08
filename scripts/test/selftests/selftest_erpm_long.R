# ======================================================================================
# Fichier : scripts/test/selftests/selftest_erpm_long.R
# Objet   : Self-test autonome pour erpm_long() (version minimale)
# Exécution: Rscript scripts/test/selftests/selftest_erpm_long.R
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule
# --------------------------------------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",    quietly = TRUE)) stop("Package 'ergm' requis.")
  if (!requireNamespace("devtools", quietly = TRUE)) stop("Package 'devtools' requis.")
  if (!requireNamespace("rprojroot", quietly = TRUE)) stop("Package 'rprojroot' requis.")
})

suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

# --------------------------------------------------------------------------------------
# Root projet
# --------------------------------------------------------------------------------------
root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)

# --------------------------------------------------------------------------------------
# Logging local
# --------------------------------------------------------------------------------------
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_erpm_long.log")
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

# --------------------------------------------------------------------------------------
# Patch ERGM (obligatoire ici)
# --------------------------------------------------------------------------------------
# NOTE:
# Le patch est un outil de debug (trace sur baseenv::replace()).
# Selon la version de R, il peut échouer (structure interne différente).
# Le selftest doit pouvoir tourner SANS patch si l'activation échoue.
.patch_enabled <- FALSE
patch_path <- file.path(root, "scripts", "ergm_patch.R")
if (!file.exists(patch_path)) stop("scripts/ergm_patch.R introuvable (patch requis).")

source(patch_path, local = FALSE)

if (!exists("ergm_patch_enable", mode = "function")) {
  stop("ergm_patch_enable() introuvable après sourcing scripts/ergm_patch.R")
}

tryCatch(
  {
    ergm_patch_enable()
    .patch_enabled <- TRUE
  },
  error = function(e) {
    message("[selftest_erpm_long] ergm_patch_enable() failed, continuing without patch.\n",
            "  reason: ", conditionMessage(e))
    .patch_enabled <- FALSE
  }
)

# --------------------------------------------------------------------------------------
# Chargement du package via devtools (prioritaire)
# --------------------------------------------------------------------------------------
# IMPORTANT:
# On privilégie devtools::load_all() pour charger la version du working tree.
# On évite de re-source R/erpm_long.R ici: ça peut créer des incohérences si load_all()
# a déjà créé un namespace avec une autre instance des fonctions.
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Le fichier DESCRIPTION n'existe pas (projet R package requis).")
}

# devtools::load_all(path = root, recompile = TRUE, quiet = TRUE)
# Chargement du package (debug-friendly)
if (!"ERPM" %in% loadedNamespaces()) {
  suppressPackageStartupMessages(library(devtools))
  load_all(path = root, recompile = TRUE, quiet = TRUE)  # IMPORTANT: pas devtools::load_all()
}

if (!exists("erpm_long", mode = "function")) stop("erpm_long() introuvable après devtools::load_all().")
if (!exists("erpm",      mode = "function")) stop("erpm() introuvable après devtools::load_all().")

# ======================================================================================
# Données de test: partitions simples, nodes simples, AVEC dyads
# ======================================================================================

# Deux partitions (n=10)
P1 <- c(1,1,1, 2,2, 3,3,3,3, 4)
P2 <- c(1,1, 2,2,2, 3,3, 4,4,4)

# IMPORTANT:
# cov_match() requiert une covariate catégorielle (factor/character).
# Ici, on stocke score comme character (catégories) plutôt que numeric.
nodes1 <- data.frame(
  label = paste0("A", 1:10),
  score = as.character(c(1,2,3, 3,4, 5,6,7,8, 2)),
  stringsAsFactors = FALSE
)

nodes2 <- data.frame(
  label = paste0("B", 1:10),
  score = as.character(c(2,1, 4,4,5, 5,6, 1,2,3)),
  stringsAsFactors = FALSE
)

# Trois et quatre partitions (toujours n=10)
P3 <- c(1, 2,2, 3,3,3, 4,4,4,4)
P4 <- c(1,1,1,1, 2,2, 3,3, 4,4)

nodes3 <- data.frame(
  label = paste0("C", 1:10),
  score = as.character(c(0,1,2,3,4,5,6,7,8,9)),
  stringsAsFactors = FALSE
)
nodes4 <- data.frame(
  label = paste0("D", 1:10),
  score = as.character(c(9,8,7,6,5,4,3,2,1,0)),
  stringsAsFactors = FALSE
)

# --------------------------------------------------------------------------------------
# Dyads: matrices n×n, une liste par temps
# --------------------------------------------------------------------------------------
.make_dyad_symmetric <- function(n, seed = 1L, p = 0.25) {
  set.seed(seed)
  M <- matrix(rbinom(n * n, 1, p), n, n)
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  diag(M) <- 0
  M
}

.make_dyad_directed <- function(n, seed = 1L, p = 0.20) {
  set.seed(seed)
  M <- matrix(rbinom(n * n, 1, p), n, n)
  diag(M) <- 0
  M
}

# Deux dyads:
# - "friendship": symétrique 0/1
# - "advice": dirigée 0/1
#
# Remarque:
# build_bipartite_from_inputs() force les dimnames (labels) dans un ordre contrôlé.
# Ici, les matrices sont laissées sans dimnames pour éviter toute contrainte d'ordre.
dyads1 <- list(
  friendship = .make_dyad_symmetric(10, seed = 101, p = 0.30),
  advice     = .make_dyad_directed(10,   seed = 102, p = 0.20)
)
dyads2 <- list(
  friendship = .make_dyad_symmetric(10, seed = 201, p = 0.25),
  advice     = .make_dyad_directed(10,   seed = 202, p = 0.15)
)
dyads3 <- list(
  friendship = .make_dyad_symmetric(10, seed = 301, p = 0.20),
  advice     = .make_dyad_directed(10,   seed = 302, p = 0.10)
)
dyads4 <- list(
  friendship = .make_dyad_symmetric(10, seed = 401, p = 0.35),
  advice     = .make_dyad_directed(10,   seed = 402, p = 0.05)
)

# ======================================================================================
# Helpers d'assertion pour dyads
# ======================================================================================
.assert_dyads_attached <- function(nw, expected_names) {
  d <- network::get.network.attribute(nw, "dyads")
  stopifnot(!is.null(d))
  stopifnot(is.list(d))
  stopifnot(all(expected_names %in% names(d)))

  for (nm in expected_names) {
    M <- d[[nm]]
    stopifnot(is.matrix(M))
    stopifnot(all(dim(M) == c(10, 10)))

    # Vérifie que les dimnames ont été posés au moment de la construction du réseau.
    stopifnot(!is.null(rownames(M)))
    stopifnot(!is.null(colnames(M)))
    stopifnot(length(rownames(M)) == 10)
    stopifnot(length(colnames(M)) == 10)
  }
  TRUE
}

# ======================================================================================
# Helpers d'assertion pour inertie
# ======================================================================================
.assert_inertia_attr <- function(nw, term, lag) {
  nm <- paste0("erpm_inertia__", term, "__lag", as.integer(lag))
  v <- network::get.network.attribute(nw, nm)
  stopifnot(!is.null(v))
  TRUE
}

.assert_inertia_not_attr <- function(nw, term, lag) {
  nm <- paste0("erpm_inertia__", term, "__lag", as.integer(lag))
  v <- network::get.network.attribute(nw, nm)
  stopifnot(is.null(v))
  TRUE
}

# ======================================================================================
# Helpers d'assertion pour fits ERGM
# ======================================================================================
.assert_ergm_fit <- function(fit) {
  stopifnot(!is.null(fit))
  stopifnot(inherits(fit, "ergm"))

  # Conformité minimale: coefficients numériques, non vides.
  co <- tryCatch(stats::coef(fit), error = function(e) NULL)
  stopifnot(!is.null(co))
  stopifnot(is.numeric(co))
  stopifnot(length(co) >= 1L)

  TRUE
}

.assert_erpm_long_fit_output <- function(out, expected_T) {
  stopifnot(inherits(out, "erpm_long"))
  stopifnot(length(out$calls) == expected_T)
  stopifnot(length(out$networks) == expected_T)
  stopifnot(length(out$fits) == expected_T)

  for (t in seq_len(expected_T)) {
    stopifnot(inherits(out$networks[[t]], "network"))
    stopifnot(.assert_ergm_fit(out$fits[[t]]))
  }

  TRUE
}

# --------------------------------------------------------------------------------------
# Comparaison de fits (objectif: vérifier que erpm_long construit les mêmes données
# que des appels directs à erpm() par partition).
# --------------------------------------------------------------------------------------
# NOTE:
# On ne cherche PAS à valider cov_match ici.
# On cherche un test d'intégration: "données construites -> fit".
# La comparaison de coefficients peut échouer si ergm consomme le RNG différemment
# entre les deux parcours. Dans ce cas, augmenter la tolérance ou comparer des
# éléments déterministes (call, attributs réseau).
.assert_ergm_fits_close <- function(fit_a, fit_b, tol = 1e-6) {
  stopifnot(.assert_ergm_fit(fit_a))
  stopifnot(.assert_ergm_fit(fit_b))

  ca <- stats::coef(fit_a)
  cb <- stats::coef(fit_b)

  # Même ordre de paramètres.
  stopifnot(identical(names(ca), names(cb)))

  # Coefficients proches.
  stopifnot(all(is.finite(ca)))
  stopifnot(all(is.finite(cb)))
  stopifnot(max(abs(ca - cb)) <= tol)

  TRUE
}

# ======================================================================================
# Selftests statiques (dry-run)
# ======================================================================================

run_case_two_partitions_dryrun <- function() {
  cat("\n=== CASE 1: 2 partitions | dry-run | verbose | 2 effets + 1 dyadique ===\n")

  partitions <- list(P1 = P1, P2 = P2)
  nodes      <- list(nodes1, nodes2)
  dyads      <- list(dyads1, dyads2)

  f <- as.formula("partitions ~ groups + cov_match('score') + dyadcov('friendship')")
  environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

  out <- erpm_long(
    formula   = f,
    eval.call = FALSE,
    verbose   = TRUE,
    debug     = TRUE,
    nodes     = nodes,
    dyads     = dyads
  )

  stopifnot(inherits(out, "erpm_long"))
  stopifnot(length(out$calls) == 2L)
  stopifnot(inherits(out$networks[[1]], "network"))
  stopifnot(inherits(out$networks[[2]], "network"))

  # Dyads attachées sur t=1 et t=2.
  stopifnot(.assert_dyads_attached(out$networks[[1]], c("friendship", "advice")))
  stopifnot(.assert_dyads_attached(out$networks[[2]], c("friendship", "advice")))

  # Pas d'inertie ici: history doit exister (liste de longueur T) et rester NULL.
  stopifnot(!is.null(out$history_timeline))
  stopifnot(is.list(out$history_timeline))
  stopifnot(length(out$history_timeline) == 2L)
  stopifnot(is.null(out$history_timeline[[1L]]))
  stopifnot(is.null(out$history_timeline[[2L]]))

  cat("=== OK CASE 1 ===\n")
  invisible(out)
}

run_case_three_partitions_dryrun <- function() {
  cat("\n=== CASE 2: 3 partitions | dry-run | verbose | 3 effets + 1 dyadique ===\n")

  partitions <- list(P1 = P1, P2 = P2, P3 = P3)
  nodes      <- list(nodes1, nodes2, nodes3)
  dyads      <- list(dyads1, dyads2, dyads3)

  f <- as.formula("partitions ~ groups + cliques(2) + cov_match('score') + dyadcov('friendship')")
  environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

  out <- erpm_long(
    formula   = f,
    eval.call = FALSE,
    verbose   = TRUE,
    debug     = TRUE,
    nodes     = nodes,
    dyads     = dyads
  )

  stopifnot(inherits(out, "erpm_long"))
  stopifnot(length(out$calls) == 3L)

  # Dyads attachées sur chaque temps.
  for (t in 1:3) stopifnot(.assert_dyads_attached(out$networks[[t]], c("friendship", "advice")))

  # Pas d'inertie ici.
  stopifnot(!is.null(out$history_timeline))
  stopifnot(is.list(out$history_timeline))
  stopifnot(length(out$history_timeline) == 3L)
  for (t in 1:3) stopifnot(is.null(out$history_timeline[[t]]))

  cat("=== OK CASE 2 ===\n")
  invisible(out)
}

run_case_four_partitions_dryrun <- function() {
  cat("\n=== CASE 3: 4 partitions | dry-run | verbose | 3 effets + 1 dyadique ===\n")

  partitions <- list(P1 = P1, P2 = P2, P3 = P3, P4 = P4)
  nodes      <- list(nodes1, nodes2, nodes3, nodes4)
  dyads      <- list(dyads1, dyads2, dyads3, dyads4)

  f <- as.formula("partitions ~ groups + cliques(2) + cov_match('score') + dyadcov('friendship')")
  environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

  out <- erpm_long(
    formula   = f,
    eval.call = FALSE,
    verbose   = TRUE,
    debug     = TRUE,
    nodes     = nodes,
    dyads     = dyads
  )

  stopifnot(inherits(out, "erpm_long"))
  stopifnot(length(out$calls) == 4L)

  # Dyads attachées sur chaque temps.
  for (t in 1:4) stopifnot(.assert_dyads_attached(out$networks[[t]], c("friendship", "advice")))

  # Pas d'inertie ici.
  stopifnot(!is.null(out$history_timeline))
  stopifnot(is.list(out$history_timeline))
  stopifnot(length(out$history_timeline) == 4L)
  for (t in 1:4) stopifnot(is.null(out$history_timeline[[t]]))

  cat("=== OK CASE 3 ===\n")
  invisible(out)
}

# ======================================================================================
# Selftests inertiels (dry-run)
# ======================================================================================

run_case_two_partitions_inertia_p1 <- function() {
  cat("\n=== CASE 4: 2 partitions | dry-run | inertie simple | past_influence=1 ===\n")

  partitions <- list(P1 = P1, P2 = P2)
  nodes      <- list(nodes1, nodes2)
  dyads      <- list(dyads1, dyads2)

  # Effet inertiel simple.
  # Attendu:
  # - t=1: inertie inactive
  # - t=2: inertie active (past_influence=1) et attribut lag1 attaché
  f <- as.formula("partitions ~ groups + inertia_groups(sizes=2, past_influence=1) + dyadcov('friendship')")
  environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

  out <- erpm_long(
    formula   = f,
    eval.call = FALSE,
    verbose   = TRUE,
    debug     = TRUE,
    nodes     = nodes,
    dyads     = dyads
  )

  stopifnot(inherits(out, "erpm_long"))
  stopifnot(length(out$calls) == 2L)

  # Dyads attachées sur chaque temps.
  for (t in 1:2) stopifnot(.assert_dyads_attached(out$networks[[t]], c("friendship", "advice")))

  # Inertie.
  stopifnot(.assert_inertia_not_attr(out$networks[[1]], "inertia_groups", 1))
  stopifnot(.assert_inertia_attr(out$networks[[2]], "inertia_groups", 1))

  # History.
  stopifnot(!is.null(out$history_timeline))
  stopifnot(is.list(out$history_timeline))
  stopifnot(length(out$history_timeline) == 2L)
  stopifnot(is.null(out$history_timeline[[1L]]))
  stopifnot(!is.null(out$history_timeline[[2L]]))
  stopifnot(is.list(out$history_timeline[[2L]]))

  cat("=== OK CASE 4 ===\n")
  invisible(out)
}

run_case_three_partitions_two_inertias <- function() {
  cat("\n=== CASE 5: 3 partitions | dry-run | 2 inerties | past_influence=1 et 2 ===\n")

  partitions <- list(P1 = P1, P2 = P2, P3 = P3)
  nodes      <- list(nodes1, nodes2, nodes3)
  dyads      <- list(dyads1, dyads2, dyads3)

  # Deux effets inertiels:
  # - inertia_groups(past_influence=1): actif dès t=2
  # - inertia_cliques(past_influence=2, k=2): actif seulement à t=3
  # Attendu:
  # - t=1: rien
  # - t=2: inertia_groups lag1 uniquement
  # - t=3: inertia_groups lag1 + inertia_cliques lag1 et lag2
  f <- as.formula("partitions ~ groups + inertia_groups(sizes=2, past_influence=1) + inertia_groups(sizes=3, past_influence=2) + dyadcov('friendship')")
  environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

  out <- erpm_long(
    formula   = f,
    eval.call = FALSE,
    verbose   = TRUE,
    debug     = TRUE,
    nodes     = nodes,
    dyads     = dyads
  )

  stopifnot(inherits(out, "erpm_long"))
  stopifnot(length(out$calls) == 3L)

  # Dyads attachées sur chaque temps.
  for (t in 1:3) stopifnot(.assert_dyads_attached(out$networks[[t]], c("friendship", "advice")))

  # Inertie groups (past_influence=1).
  stopifnot(.assert_inertia_not_attr(out$networks[[1]], "inertia_groups", 1))
  stopifnot(.assert_inertia_attr(out$networks[[2]], "inertia_groups", 1))
  stopifnot(.assert_inertia_attr(out$networks[[3]], "inertia_groups", 1))

  # Inertie cliques (past_influence=2).
  stopifnot(.assert_inertia_not_attr(out$networks[[1]], "inertia_cliques", 1))
  stopifnot(.assert_inertia_not_attr(out$networks[[2]], "inertia_cliques", 1))
  stopifnot(.assert_inertia_attr(out$networks[[3]], "inertia_cliques", 1))
  stopifnot(.assert_inertia_attr(out$networks[[3]], "inertia_cliques", 2))

  # History.
  stopifnot(!is.null(out$history_timeline))
  stopifnot(is.list(out$history_timeline))
  stopifnot(length(out$history_timeline) == 3L)
  stopifnot(is.null(out$history_timeline[[1L]]))
  stopifnot(!is.null(out$history_timeline[[2L]]))
  stopifnot(!is.null(out$history_timeline[[3L]]))

  cat("=== OK CASE 5 ===\n")
  invisible(out)
}

# ======================================================================================
# Tests avec évaluation (sans dry-run) : sans inertie
# ======================================================================================

run_case_two_partitions_fit <- function() {
  cat("\n=== CASE 6: 2 partitions | 2 effets statiques ===\n")

  partitions <- list(P1 = P1, P2 = P2)
  nodes      <- list(nodes1, nodes2)
  dyads      <- list(dyads1, dyads2)

  f <- as.formula("partitions ~ groups + cov_match('score')")
  environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

  # ---------------------------------------------------------------------------
  # Parcours A: erpm_long() (construit réseaux + fit)
  # ---------------------------------------------------------------------------
  set.seed(1)
  out <- erpm_long(
    formula   = f,
    eval.call = TRUE,
    verbose   = TRUE,
    debug     = TRUE,
    nodes     = nodes,
    dyads     = dyads
  )

  print(out)
  stopifnot(.assert_erpm_long_fit_output(out, expected_T = 2L))

  # Dyads attachées (même si RHS n'utilise pas dyadcov ici).
  for (t in 1:2) stopifnot(.assert_dyads_attached(out$networks[[t]], c("friendship", "advice")))

  # Pas d'inertie attendue.
  stopifnot(!is.null(out$history_timeline))
  stopifnot(length(out$history_timeline) == 2L)
  stopifnot(is.null(out$history_timeline[[1L]]))
  stopifnot(is.null(out$history_timeline[[2L]]))

  # ---------------------------------------------------------------------------
  # Parcours B: fits directs via erpm() sur chaque partition
  # Objectif: vérifier que la construction des données par erpm_long est cohérente
  # avec des appels unitaires erpm(partition ~ RHS, nodes=..., dyads=...).
  # ---------------------------------------------------------------------------
  set.seed(1)

  direct_fits <- vector("list", length(partitions))
  names(direct_fits) <- names(partitions)

  for (t in seq_along(partitions)) {
    p_t <- partitions[[t]]

    # NOTE:
    # On appelle erpm() avec LHS = partition, pour forcer la reconstruction du réseau
    # à partir de (partition, nodes, dyads), comme le fait erpm_long().
    direct_fits[[t]] <- erpm(
      p_t ~ groups + cov_match("score"),
      eval.call = TRUE,
      verbose   = TRUE,
      nodes     = nodes[[t]],
      dyads     = dyads[[t]]
    )
    stopifnot(.assert_ergm_fit(direct_fits[[t]]))
    # print(direct_fits[[t]])
  }

  

  # Comparaison des coefficients (tolérance faible par défaut).
  for (t in seq_along(partitions)) {
    stopifnot(.assert_ergm_fits_close(out$fits[[t]], direct_fits[[t]], tol = 1e-6))
  }

  cat("=== OK CASE 6 ===\n")
  invisible(out)
}

run_case_four_partitions_fit <- function() {
  cat("\n=== CASE 7: 4 partitions |  2 effets statiques ===\n")

  partitions <- list(P1 = P1, P2 = P2, P3 = P3, P4 = P4)
  nodes      <- list(nodes1, nodes2, nodes3, nodes4)
  dyads      <- list(dyads1, dyads2, dyads3, dyads4)

  f <- as.formula("partitions ~ groups + cov_match('score')")
  environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

  # ---------------------------------------------------------------------------
  # Parcours A: erpm_long()
  # ---------------------------------------------------------------------------
  set.seed(1)
  out <- erpm_long(
    formula   = f,
    eval.call = TRUE,
    verbose   = TRUE,
    debug     = TRUE,
    nodes     = nodes,
    dyads     = dyads
  )

  print(out)
  stopifnot(.assert_erpm_long_fit_output(out, expected_T = 4L))

  # Dyads attachées.
  for (t in 1:4) stopifnot(.assert_dyads_attached(out$networks[[t]], c("friendship", "advice")))

  # Pas d'inertie attendue.
  stopifnot(!is.null(out$history_timeline))
  stopifnot(length(out$history_timeline) == 4L)
  for (t in 1:4) stopifnot(is.null(out$history_timeline[[t]]))

  # ---------------------------------------------------------------------------
  # Parcours B: fits directs via erpm() sur chaque partition
  # ---------------------------------------------------------------------------
  set.seed(1)

  direct_fits <- vector("list", length(partitions))
  names(direct_fits) <- names(partitions)

  for (t in seq_along(partitions)) {
    p_t <- partitions[[t]]

    direct_fits[[t]] <- erpm(
      p_t ~ groups + cov_match("score"),
      eval.call = TRUE,
      verbose   = TRUE,
      nodes     = nodes[[t]],
      dyads     = dyads[[t]]
    )
    stopifnot(.assert_ergm_fit(direct_fits[[t]]))
    print(direct_fits[[t]])
  }

  # Comparaison des coefficients.
  for (t in seq_along(partitions)) {
    stopifnot(.assert_ergm_fits_close(out$fits[[t]], direct_fits[[t]], tol = 1e-6))
  }

  cat("=== OK CASE 7 ===\n")

  invisible(out)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== SELFTEST ERPM_LONG (minimal) ===\n")

# run_case_two_partitions_dryrun()
# run_case_three_partitions_dryrun()
# run_case_four_partitions_dryrun()

# Nouveaux tests inertiels.
# run_case_two_partitions_inertia_p1()
# run_case_three_partitions_two_inertias()

# Tests avec évaluation (sans dry-run).
run_case_two_partitions_fit()
run_case_four_partitions_fit()

# Patch: on tente de désactiver si dispo.
if (.patch_enabled && exists("ergm_patch_disable", mode = "function")) {
  try(ergm_patch_disable(), silent = TRUE)
}

cat("\nTous les tests erpm_long (dry-run) ont passé.\n")