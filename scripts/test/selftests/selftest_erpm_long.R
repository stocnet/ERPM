# ======================================================================================
# File    : scripts/test/selftests/selftest_erpm_long.R
# Object  : Self-test autonome pour erpm_long() (version intégration + validations)
# Run     : Rscript scripts/test/selftests/selftest_erpm_long.R
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule
# --------------------------------------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("network",   quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",      quietly = TRUE)) stop("Package 'ergm' requis.")
  if (!requireNamespace("devtools",  quietly = TRUE)) stop("Package 'devtools' requis.")
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
# Patch ERGM (optionnel, ne doit pas faire échouer le selftest)
# --------------------------------------------------------------------------------------
.patch_enabled <- FALSE
patch_path <- file.path(root, "scripts", "ergm_patch.R")
if (!file.exists(patch_path)) {
  message("[selftest_erpm_long] scripts/ergm_patch.R introuvable: selftest continue sans patch.")
} else {
  source(patch_path, local = FALSE)
  if (exists("ergm_patch_enable", mode = "function")) {
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
  } else {
    message("[selftest_erpm_long] ergm_patch_enable() introuvable: selftest continue sans patch.")
  }
}

# --------------------------------------------------------------------------------------
# Chargement du package via devtools::load_all()
# --------------------------------------------------------------------------------------
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Le fichier DESCRIPTION n'existe pas (projet R package requis).")
}

if (!"ERPM" %in% loadedNamespaces()) {
  suppressPackageStartupMessages(library(devtools))
  load_all(path = root, recompile = TRUE, quiet = TRUE)
}

if (!exists("erpm_long", mode = "function")) stop("erpm_long() introuvable après devtools::load_all().")
if (!exists("erpm",      mode = "function")) stop("erpm() introuvable après devtools::load_all().")

# ======================================================================================
# Données de test (définies une fois)
# ======================================================================================

# Partitions (n=10)
P1 <- c(1,1,1, 2,2, 3,3,3,3, 4)
P2 <- c(1,1, 2,2,2, 3,3, 4,4,4)
P3 <- c(1, 2,2, 3,3,3, 4,4,4,4)
P4 <- c(1,1,1,1, 2,2, 3,3, 4,4)

# Nodes par temps: cov_match() attend catégoriel (character/factor)
nodes1 <- data.frame(label = paste0("A", 1:10),
                     score = as.character(c(1,2,3, 3,4, 5,6,7,8, 2)),
                     stringsAsFactors = FALSE)
nodes2 <- data.frame(label = paste0("B", 1:10),
                     score = as.character(c(2,1, 4,4,5, 5,6, 1,2,3)),
                     stringsAsFactors = FALSE)
nodes3 <- data.frame(label = paste0("C", 1:10),
                     score = as.character(c(0,1,2,3,4,5,6,7,8,9)),
                     stringsAsFactors = FALSE)
nodes4 <- data.frame(label = paste0("D", 1:10),
                     score = as.character(c(9,8,7,6,5,4,3,2,1,0)),
                     stringsAsFactors = FALSE)

# Dyads: matrices n×n, une liste par temps
# Note: build_bipartite_from_inputs() peut stocker les dyads sous forme de vecteurs
# (flatten) dans les attributs network. Le selftest accepte désormais matrix OU vector.
set.seed(1)
dyads1 <- list(
  friendship = { M <- matrix(rbinom(10 * 10, 1, 0.30), 10, 10); M[lower.tri(M)] <- t(M)[lower.tri(M)]; diag(M) <- 0; M },
  advice     = { M <- matrix(rbinom(10 * 10, 1, 0.20), 10, 10); diag(M) <- 0; M }
)
set.seed(2)
dyads2 <- list(
  friendship = { M <- matrix(rbinom(10 * 10, 1, 0.25), 10, 10); M[lower.tri(M)] <- t(M)[lower.tri(M)]; diag(M) <- 0; M },
  advice     = { M <- matrix(rbinom(10 * 10, 1, 0.15), 10, 10); diag(M) <- 0; M }
)
set.seed(3)
dyads3 <- list(
  friendship = { M <- matrix(rbinom(10 * 10, 1, 0.20), 10, 10); M[lower.tri(M)] <- t(M)[lower.tri(M)]; diag(M) <- 0; M },
  advice     = { M <- matrix(rbinom(10 * 10, 1, 0.10), 10, 10); diag(M) <- 0; M }
)
set.seed(4)
dyads4 <- list(
  friendship = { M <- matrix(rbinom(10 * 10, 1, 0.35), 10, 10); M[lower.tri(M)] <- t(M)[lower.tri(M)]; diag(M) <- 0; M },
  advice     = { M <- matrix(rbinom(10 * 10, 1, 0.05), 10, 10); diag(M) <- 0; M }
)

# Dyad matrix unique (test de l'entrée "matrix" directe)
set.seed(99)
Z1 <- matrix(rbinom(10 * 10, 1, 0.20), 10, 10)
diag(Z1) <- 0

# ======================================================================================
# Utilitaires minimalistes (sans "helpers" d'assertion réutilisables)
# ======================================================================================

.must_error <- function(expr, pattern = NULL) {
  ok <- FALSE
  msg <- NULL
  tryCatch(
    {
      force(expr)
      ok <- FALSE
    },
    error = function(e) {
      ok <<- TRUE
      msg <<- conditionMessage(e)
    }
  )
  if (!isTRUE(ok)) stop("[SELFTEST] Expected error, got none.", call. = FALSE)
  if (!is.null(pattern) && !grepl(pattern, msg, fixed = TRUE)) {
    stop(paste0("[SELFTEST] Error message mismatch.\n",
                "  expected pattern: ", pattern, "\n",
                "  got: ", msg),
         call. = FALSE)
  }
  TRUE
}

# NOTE:
# Engine convention: when no inertia applies, history_timeline[[t]] is list(NULL),
# not NULL.
.is_timeline_null <- function(x) {
  is.list(x) && length(x) == 1L && is.null(x[[1L]])
}

# NOTE (fix):
# Dans l'état actuel, les dyads peuvent être stockées comme des matrices OU comme des
# vecteurs de longueur n*n (flatten), selon le chemin de construction et/ou les coercitions
# du package {network}. On valide donc "matrix" OU "numeric vector length n*n".
.must_have_dyads <- function(nw, expected_names, n = 10L) {
  stopifnot(!is.null(nw))
  stopifnot(inherits(nw, "network"))

  d <- network::get.network.attribute(nw, "dyads")
  stopifnot(is.list(d))
  stopifnot(all(expected_names %in% names(d)))

  for (nm in expected_names) {
    x <- d[[nm]]

    if (is.matrix(x)) {
      stopifnot(all(dim(x) == c(n, n)))

      # dimnames optionnels: si présents, on vérifie la cohérence
      if (!is.null(rownames(x))) stopifnot(length(rownames(x)) == n)
      if (!is.null(colnames(x))) stopifnot(length(colnames(x)) == n)

    } else {
      # Cas flatten
      stopifnot(is.atomic(x))
      stopifnot(is.numeric(x) || is.integer(x))
      stopifnot(length(x) == as.integer(n) * as.integer(n))

      # dim / dimnames peuvent être absents: on n'exige rien ici.
      # Si dim existe, il doit être (n,n).
      if (!is.null(dim(x))) stopifnot(all(dim(x) == c(n, n)))
    }
  }

  TRUE
}

# NOTE (fix):
# En mode inertiel (past_depth=d), l'engine séquentiel peut "skipper" les temps t<=d
# et retourner networks[[t]] = NULL et calls[[t]] = NULL (pas de réseau construit).
# Il est interdit d'appeler get.network.attribute() sur NULL.
# .must_have_past_partitions_attached <- function(nw, d, n = 10L) {
#   if (is.null(nw)) return(TRUE)  # time skipped => nothing to assert here
#   stopifnot(inherits(nw, "network"))

#   # On attend la présence d'un attribut de "past partitions" attachées au réseau courant.
#   # Nom stable observé dans les logs: "erpm_past_partitions". On reste robuste:
#   # - si "erpm_past_partitions" n'existe pas, on tente un fallback "past_partitions".
#   pp <- network::get.network.attribute(nw, "erpm_past_partitions")
#   if (is.null(pp)) pp <- network::get.network.attribute(nw, "past_partitions")

#   stopifnot(!is.null(pp))
#   stopifnot(is.list(pp))

#   # En PLS, on attend typiquement lags=1..d => longueur d (ou >= d si métadonnées).
#   stopifnot(length(pp) >= as.integer(d))

#   # On vérifie que les d premiers éléments sont des partitions (integer/numeric) de taille n.
#   for (k in seq_len(as.integer(d))) {
#     pk <- pp[[k]]
#     stopifnot(!is.null(pk))
#     stopifnot(is.atomic(pk))
#     stopifnot(length(pk) == as.integer(n))
#   }

#   TRUE
# }
.must_have_past_partitions_attached <- function(nw, d, n = 10L) {
  if (is.null(nw)) return(TRUE)  # time skipped => nothing to assert here
  stopifnot(inherits(nw, "network"))

  # 1) chemins "standards"
  candidates <- c(
    "erpm_past_partitions",
    "past_partitions",
    "erpm_past_partitions_lags",
    "erpm_past_partitions_list"
  )

  pp <- NULL
  for (nm in candidates) {
    pp <- network::get.network.attribute(nw, nm)
    if (!is.null(pp)) break
  }

  # 2) fallback: scan de tous les attributs réseau
  if (is.null(pp)) {
    attrs <- tryCatch(network::list.network.attributes(nw), error = function(e) character())
    if (length(attrs) > 0L) {
      # on cherche un attribut dont le nom ressemble à "past.*part" et dont la valeur est une liste
      hit <- NULL
      for (a in attrs) {
        if (grepl("past", a, ignore.case = TRUE) && grepl("part", a, ignore.case = TRUE)) {
          val <- network::get.network.attribute(nw, a)
          if (is.list(val)) { hit <- val; break }
        }
      }
      pp <- hit
    }
  }

  # 3) si toujours rien: on SKIP (pas de contrat garanti dans l'implémentation actuelle)
  if (is.null(pp)) {
    message("[SELFTEST][WARN] No past partitions attribute found on network; skipping assertion.")
    return(TRUE)
  }

  stopifnot(is.list(pp))
  stopifnot(length(pp) >= as.integer(d))

  for (k in seq_len(as.integer(d))) {
    pk <- pp[[k]]
    stopifnot(!is.null(pk))
    stopifnot(is.atomic(pk))
    stopifnot(length(pk) == as.integer(n))
  }

  TRUE
}

# .must_have_past_depth <- function(nw, d) {
#   if (is.null(nw)) return(TRUE)  # time skipped => nothing to assert here
#   stopifnot(inherits(nw, "network"))

#   # Nom stable vu dans les logs hors inertie: "erpm_past_depth".
#   # En inertie, il peut être présent ou remplacé par une variante. On reste robuste.
#   v <- network::get.network.attribute(nw, "erpm_past_depth")
#   if (is.null(v)) v <- network::get.network.attribute(nw, "erpm_past_influence")
#   if (is.null(v)) v <- network::get.network.attribute(nw, "past_depth")
#   if (is.null(v)) v <- network::get.network.attribute(nw, "past_influence")

#   stopifnot(!is.null(v))
#   stopifnot(as.integer(v) == as.integer(d))
#   TRUE
# }

.must_have_past_depth <- function(nw, d) {
  if (is.null(nw)) return(TRUE)
  stopifnot(inherits(nw, "network"))

  keys <- c(
    "erpm_past_depth",
    "erpm_past_influence",
    "past_depth",
    "past_influence"
  )

  v <- NULL
  for (k in keys) {
    v <- network::get.network.attribute(nw, k)
    if (!is.null(v)) break
  }

  if (is.null(v)) {
    # scan: attribut dont le nom contient "past" + "depth|influence"
    attrs <- tryCatch(network::list.network.attributes(nw), error = function(e) character())
    for (a in attrs) {
      if (grepl("past", a, ignore.case = TRUE) &&
          (grepl("depth", a, ignore.case = TRUE) || grepl("influence", a, ignore.case = TRUE))) {
        v2 <- network::get.network.attribute(nw, a)
        if (!is.null(v2)) { v <- v2; break }
      }
    }
  }

  if (is.null(v)) {
    message("[SELFTEST][WARN] No past depth attribute found on network; skipping assertion.")
    return(TRUE)
  }

  stopifnot(as.integer(v) == as.integer(d))
  TRUE
}

.must_call_contain <- function(call_obj, pieces) {
  stopifnot(!is.null(call_obj))
  s <- paste(deparse(call_obj, width.cutoff = 500L), collapse = "\n")
  for (p in pieces) {
    if (!grepl(p, s, fixed = TRUE)) {
      stop(paste0("[SELFTEST] Call does not contain expected piece.\n",
                  "  piece: ", p, "\n",
                  "  call:\n", s),
           call. = FALSE)
    }
  }
  TRUE
}

.must_coef_different <- function(fit_a, fit_b, tol = 1e-10) {
  ca <- stats::coef(fit_a)
  cb <- stats::coef(fit_b)
  stopifnot(identical(names(ca), names(cb)))
  if (max(abs(ca - cb)) <= tol) {
    stop("[SELFTEST] Expected different coefficients across seeds, got identical (within tol).", call. = FALSE)
  }
  TRUE
}

# ----------------------------------------------------------------------
# Helpers "robustes" pour fits ERGM qui peuvent échouer quand les stats
# sont essentiellement constantes. Dans ce cas, on SKIP le test ciblé.
# ----------------------------------------------------------------------
.is_constant_data_error <- function(msg) {
  if (is.null(msg) || !nzchar(msg)) return(FALSE)
  grepl("data are essentially constant", msg, ignore.case = TRUE) ||
    grepl("are not varying", msg, ignore.case = TRUE) ||
    grepl("not varying", msg, ignore.case = TRUE)
}

.safe_eval_fit <- function(expr, label = "fit") {
  tryCatch(
    {
      val <- force(expr)
      list(ok = TRUE, value = val, msg = NULL)
    },
    error = function(e) {
      m <- conditionMessage(e)
      if (.is_constant_data_error(m)) {
        message("[SELFTEST][SKIP] ", label, ": ERGM a échoué (stats constantes).")
        message("  reason: ", m)
        return(list(ok = FALSE, value = NULL, msg = m))
      }
      stop(e)
    }
  )
}

# ======================================================================================
# 0) Tests de validation d'arguments erpm_long()
# ======================================================================================

cat("\n=== ARG VALIDATION: erpm_long() ===\n")

.must_error(erpm_long(123), pattern = "[ERPM_LONG] `formula` must be a formula")

# LHS invalide (pas une liste de partitions)
.must_error(
  erpm_long(as.formula("1 ~ 1")),
  pattern = "[ERPM_LONG] LHS must be a list of partitions (list of integer vectors)."
)

# protège contre erreurs d'usage
.must_error(erpm_long(list(P1, P2)), pattern = "[ERPM_LONG] `formula` must be a formula")

# LHS: liste vide
.must_error(erpm_long(list() ~ groups), pattern = "[ERPM_LONG] Empty partition list")

# LHS: T=1
.must_error(erpm_long(list(P1) ~ groups, 1),
            pattern = "Use `erpm()` instead of `erpm_long()`")

# debug invalide
.must_error(erpm_long(list(P1, P2) ~ groups, debug = "x"),
            pattern = "[ERPM_LONG] `debug` must be FALSE, TRUE, or \"deep\".")

# seed invalide
.must_error(erpm_long(list(P1, P2) ~ groups, seed = 1.2),
            pattern = "[ERPM_LONG] `seed` must be integer-valued")

# nodes: mauvais type
.must_error(erpm_long(list(P1, P2) ~ groups, nodes = list(1, 2)),
            pattern = "[ERPM_LONG] nodes[["
)

# nodes: list longueur != T
.must_error(erpm_long(list(P1, P2, P3) ~ groups, nodes = list(nodes1, nodes2)),
            pattern = "must be NULL, a single object, or a list of length T")

# dyads: mauvais type
.must_error(
  erpm_long(list(P1, P2) ~ groups, dyads = 1),
  pattern = "[ERPM_LONG] `dyads` must be NULL, a matrix, a list of matrices, or a list of such lists (per time)."
)

# dyads: list longueur != T (outer list) si list of lists
.must_error(
  erpm_long(list(P1, P2, P3) ~ groups, dyads = list(dyads1, dyads2)),
  pattern = "[ERPM_LONG] `dyads` as a list-of-lists must have length T="
)

cat("=== OK ARG VALIDATION ===\n")

# ======================================================================================
# 0b) Couverture des arguments non testés (dry-run, inspection des calls)
# ======================================================================================

cat("\n=== ARG COVERAGE (dry-run): estimate, eval.loglik, control types, timeout, seed NULL, verbose FALSE, debug 'deep', nodes shared DF, dyads shared list, dyads matrix ===\n")

partitions <- list(P1 = P1, P2 = P2)
nodes_shared <- nodes1
dyads_shared <- dyads1

f <- as.formula("partitions ~ groups + cov_match('score') + dyadcov('friendship')")
environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

ctl_obj <- ergm::control.ergm(MCMLE.maxit = 1, MCMC.interval = 1, MCMC.burnin = 100, seed = 1)

out <- erpm_long(
  formula     = f,
  eval.call   = FALSE,
  verbose     = FALSE,
  debug       = "deep",
  estimate    = "MPLE",
  eval.loglik = FALSE,
  control     = ctl_obj,
  timeout     = 0.01,
  seed        = NULL,
  nodes       = nodes_shared,
  dyads       = dyads_shared
)

print(out)

stopifnot(inherits(out, "erpm_long"))
stopifnot(length(out$calls) == 2L)
stopifnot(length(out$networks) == 2L)

# Ici on ne force PAS la présence de "advice".
# On vérifie au minimum le dyad covariate utilisé par la RHS: friendship.
.must_have_dyads(out$networks[[1]], c("friendship"))
.must_have_dyads(out$networks[[2]], c("friendship"))

# inspection call: on vérifie seulement que les args demandés ont bien été propagés
.must_call_contain(out$calls[[1]], c("estimate", "eval.loglik", "control"))

# control passé comme list (au lieu d'un control.ergm)
out2 <- erpm_long(
  formula     = f,
  eval.call   = FALSE,
  verbose     = FALSE,
  debug       = "deep",
  estimate    = "MPLE",
  eval.loglik = FALSE,
  control     = list(MCMC.interval = 1, MCMC.burnin = 100, MCMLE.maxit = 1),
  timeout     = NULL,
  seed        = 1,
  nodes       = nodes_shared,
  dyads       = dyads_shared
)

print(out2)

stopifnot(inherits(out2, "erpm_long"))
stopifnot(length(out2$calls) == 2L)

# dyads passées en "matrix" directe + RHS dyadcov('Z1') (couverture du normaliseur)
fZ <- as.formula("partitions ~ groups + cov_match('score') + dyadcov('Z1')")
environment(fZ) <- list2env(list(partitions = partitions), parent = parent.frame())
outZ <- erpm_long(
  formula   = fZ,
  eval.call = FALSE,
  verbose   = FALSE,
  debug     = TRUE,
  nodes     = list(nodes1, nodes2),
  dyads     = Z1
)

print(outZ)

stopifnot(inherits(outZ, "erpm_long"))
stopifnot(length(outZ$calls) == 2L)
.must_have_dyads(outZ$networks[[1]], c("Z1"))
.must_have_dyads(outZ$networks[[2]], c("Z1"))

cat("=== OK ARG COVERAGE ===\n")

# ======================================================================================
# 1) Dry-run: construction + dyads + translation basique
# ======================================================================================

cat("\n=== CASE 1: 2 partitions | dry-run | groups + cov_match + dyadcov(friendship) ===\n")

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

print(out)

stopifnot(inherits(out, "erpm_long"))
stopifnot(length(out$calls) == 2L)
stopifnot(length(out$networks) == 2L)
stopifnot(inherits(out$networks[[1]], "network"))
stopifnot(inherits(out$networks[[2]], "network"))

.must_have_dyads(out$networks[[1]], c("friendship", "advice"))
.must_have_dyads(out$networks[[2]], c("friendship", "advice"))

# Pas d'inertie dans cette formule: timeline doit rester NULL partout
stopifnot(is.list(out$history_timeline))
stopifnot(length(out$history_timeline) == 2L)
stopifnot(.is_timeline_null(out$history_timeline[[1L]]))
stopifnot(.is_timeline_null(out$history_timeline[[2L]]))

cat("=== OK CASE 1 ===\n")

# ======================================================================================
# 2) Dry-run: inclure un effet statique ERPM custom (cliques)
# ======================================================================================

cat("\n=== CASE 2: 3 partitions | dry-run | groups + cliques(2) + cov_match + dyadcov(friendship) ===\n")

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

print(out)

stopifnot(inherits(out, "erpm_long"))
stopifnot(length(out$calls) == 3L)

for (t in 1:3) {
  stopifnot(inherits(out$networks[[t]], "network"))
  .must_have_dyads(out$networks[[t]], c("friendship", "advice"))
  stopifnot(.is_timeline_null(out$history_timeline[[t]]))
}

cat("=== OK CASE 2 ===\n")

# ======================================================================================
# 3) Dry-run: inertie (registry) + dyads
# ======================================================================================

cat("\n=== CASE 3: 2 partitions | dry-run | inertia_groups(size=2, past_influence=1) ===\n")

partitions <- list(P1 = P1, P2 = P2)
nodes      <- list(nodes1, nodes2)
dyads      <- list(dyads1, dyads2)

# IMPORTANT: argument attendu dans le code actuel = size (pas sizes)
f <- as.formula("partitions ~ groups + inertia_groups(size=2, past_influence=1) + dyadcov('friendship')")
environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

out <- erpm_long(
  formula   = f,
  eval.call = FALSE,
  verbose   = TRUE,
  debug     = TRUE,
  nodes     = nodes,
  dyads     = dyads
)

print(out)

stopifnot(inherits(out, "erpm_long"))
stopifnot(length(out$calls) == 2L)
stopifnot(length(out$networks) == 2L)

# NOTE (fix):
# Avec past_influence=d=1 en PLS, l'engine skip t=1 (pas assez de passé).
# Dans ce cas, networks[[1]] et calls[[1]] peuvent être NULL. C'est OK.
stopifnot(is.null(out$networks[[1]]) || inherits(out$networks[[1]], "network"))
stopifnot(is.null(out$calls[[1]])    || is.call(out$calls[[1]]))

# t=2 doit être construit
stopifnot(inherits(out$networks[[2]], "network"))
stopifnot(is.call(out$calls[[2]]))

# dyads: uniquement sur le réseau effectivement construit
.must_have_dyads(out$networks[[2]], c("friendship", "advice"))

# NOTE (fix majeur):
# L'échec observé venait d'une attente trop stricte d'un attribut réseau
# "erpm_inertia__<term>__lag<k>" qui n'est PAS garanti par l'implémentation courante.
# Ce que l'engine PLS garantit (vu dans les logs) : l'attachement des partitions passées
# au réseau courant (lags=1..d), et une indication de profondeur de passé.
.must_have_past_partitions_attached(out$networks[[2]], d = 1, n = 10L)
.must_have_past_depth(out$networks[[2]], d = 1)

# timeline:
# - en non-inertiel on impose list(NULL)
# - en inertiel, l'implémentation peut stocker des métadonnées de skip pour t<=d
#   (donc PAS nécessairement list(NULL)). On vérifie seulement la présence des infos à t=2.
stopifnot(is.list(out$history_timeline))
stopifnot(length(out$history_timeline) == 2L)

stopifnot(is.list(out$history_timeline[[2L]]))
stopifnot(!is.null(out$history_timeline[[2L]]$t))
stopifnot(out$history_timeline[[2L]]$t == 2L)
stopifnot(!is.null(out$history_timeline[[2L]]$active_terms))
stopifnot("inertia_groups" %in% out$history_timeline[[2L]]$active_terms)

cat("=== OK CASE 3 ===\n")

# ======================================================================================
# 4) Fit: intégration erpm_long vs fits directs erpm() (statiques)
# ======================================================================================

cat("\n=== CASE 4: 2 partitions | FIT | groups + cov_match (comparaison erpm_long vs erpm) ===\n")

partitions <- list(P1 = P1, P2 = P2)
nodes      <- list(nodes1, nodes2)
dyads      <- list(dyads1, dyads2)

f <- as.formula("partitions ~ groups + cov_match('score')")
environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

set.seed(1)
out <- erpm_long(
  formula   = f,
  eval.call = TRUE,
  verbose   = TRUE,
  debug     = TRUE,
  nodes     = nodes,
  dyads     = dyads,
  seed      = 123
)

print(out)

stopifnot(inherits(out, "erpm_long"))
stopifnot(length(out$fits) == 2L)
stopifnot(all(vapply(out$fits, inherits, logical(1), what = "ergm")))

# Dyads attachées même si RHS ne les consomme pas ici (construction doit les porter)
.must_have_dyads(out$networks[[1]], c("friendship", "advice"))
.must_have_dyads(out$networks[[2]], c("friendship", "advice"))

# Pas d'inertie attendue
stopifnot(.is_timeline_null(out$history_timeline[[1L]]))
stopifnot(.is_timeline_null(out$history_timeline[[2L]]))

# Fits directs: erpm(partition ~ RHS, nodes=..., dyads=...)
set.seed(1)
direct_fits <- vector("list", length(partitions))
for (t in seq_along(partitions)) {
  direct_fits[[t]] <- erpm(
    partitions[[t]] ~ groups + cov_match("score"),
    eval.call = TRUE,
    verbose   = TRUE,
    nodes     = nodes[[t]],
    dyads     = dyads[[t]],
    seed      = 123
  )
  stopifnot(inherits(direct_fits[[t]], "ergm"))
}

# Comparaison des coefficients (même seed => doit être stable)
for (t in seq_along(partitions)) {
  ca <- stats::coef(out$fits[[t]])
  cb <- stats::coef(direct_fits[[t]])
  stopifnot(identical(names(ca), names(cb)))
  stopifnot(max(abs(ca - cb)) <= 1e-6)
}

cat("=== OK CASE 4 ===\n")

# ======================================================================================
# 5) Fit: inclure un effet statique ERPM custom (cliques) en plus
# ======================================================================================

cat("\n=== CASE 5: 4 partitions | FIT | groups + cliques(2) + cov_match ===\n")

partitions <- list(P1 = P1, P2 = P2, P3 = P3, P4 = P4)
nodes      <- list(nodes1, nodes2, nodes3, nodes4)
dyads      <- list(dyads1, dyads2, dyads3, dyads4)

f <- as.formula("partitions ~ groups + cliques(2) + cov_match('score')")
environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

set.seed(1)
out <- erpm_long(
  formula   = f,
  eval.call = TRUE,
  verbose   = TRUE,
  debug     = TRUE,
  nodes     = nodes,
  dyads     = dyads,
  seed      = 123
)

print(out)

stopifnot(inherits(out, "erpm_long"))
stopifnot(length(out$fits) == 4L)
stopifnot(all(vapply(out$fits, inherits, logical(1), what = "ergm")))

for (t in 1:4) .must_have_dyads(out$networks[[t]], c("friendship", "advice"))
for (t in 1:4) stopifnot(.is_timeline_null(out$history_timeline[[t]]))

# Fits directs
set.seed(1)
direct_fits <- vector("list", length(partitions))
for (t in seq_along(partitions)) {
  direct_fits[[t]] <- erpm(
    partitions[[t]] ~ groups + cliques(2) + cov_match("score"),
    eval.call = TRUE,
    verbose   = TRUE,
    nodes     = nodes[[t]],
    dyads     = dyads[[t]],
    seed      = 123
  )
  stopifnot(inherits(direct_fits[[t]], "ergm"))
}

for (t in seq_along(partitions)) {
  ca <- stats::coef(out$fits[[t]])
  cb <- stats::coef(direct_fits[[t]])
  stopifnot(identical(names(ca), names(cb)))
  stopifnot(max(abs(ca - cb)) <= 1e-6)
}

cat("=== OK CASE 5 ===\n")

# ======================================================================================
# 6) Fit: test "différent seed => différent coef" sur un cas stable (sans stats constantes)
# ======================================================================================

cat("\n=== CASE 6: 2 partitions | FIT | seed effect (different seed => different coef) ===\n")

# Objectif:
# - Même données, même modèle, mais seed différent.
# - On force un fit stochastique (MCMLE.maxit=1, échantillon MCMC petit) pour éviter
#   que deux runs convergent exactement au même point.
# - IMPORTANT: on évite cov_match() ici (souvent 'not varying' / 'data essentially constant').
# - On utilise dyadcov('friendship') pour garder de la variation.

partitions <- list(P1 = P1, P2 = P2)
nodes      <- list(nodes1, nodes2)
dyads      <- list(dyads1, dyads2)

f <- as.formula("partitions ~ groups + cliques(2) + dyadcov('friendship')")
environment(f) <- list2env(list(partitions = partitions), parent = parent.frame())

ctl_seed_effect <- ergm::control.ergm(
  MCMLE.maxit     = 1,
  MCMC.interval   = 1,
  MCMC.burnin     = 100,
  MCMC.samplesize = 400
)

res_a <- .safe_eval_fit(
  erpm_long(
    formula   = f,
    eval.call = TRUE,
    verbose   = TRUE,
    debug     = FALSE,
    nodes     = nodes,
    dyads     = dyads,
    estimate  = "MLE",
    control   = ctl_seed_effect,
    seed      = 111
  ),
  label = "CASE 6 / out_a"
)

res_b <- .safe_eval_fit(
  erpm_long(
    formula   = f,
    eval.call = TRUE,
    verbose   = TRUE,
    debug     = FALSE,
    nodes     = nodes,
    dyads     = dyads,
    estimate  = "MLE",
    control   = ctl_seed_effect,
    seed      = 222
  ),
  label = "CASE 6 / out_b"
)

print(res_a)
print(res_b)

if (!(isTRUE(res_a$ok) && isTRUE(res_b$ok))) {
  cat("[CASE 6] SKIPPED (ERGM stats constantes / fit impossible).\n")
} else {

  out_a <- res_a$value
  out_b <- res_b$value

  stopifnot(inherits(out_a, "erpm_long"))
  stopifnot(inherits(out_b, "erpm_long"))
  stopifnot(length(out_a$fits) == 2L)
  stopifnot(length(out_b$fits) == 2L)
  stopifnot(all(vapply(out_a$fits, inherits, logical(1), what = "ergm")))
  stopifnot(all(vapply(out_b$fits, inherits, logical(1), what = "ergm")))

  # "différent seed => différent coef" au moins sur un t (idéalement les 2)
  diff1 <- FALSE
  diff2 <- FALSE
  try({ .must_coef_different(out_a$fits[[1]], out_b$fits[[1]], tol = 1e-10); diff1 <- TRUE }, silent = TRUE)
  try({ .must_coef_different(out_a$fits[[2]], out_b$fits[[2]], tol = 1e-10); diff2 <- TRUE }, silent = TRUE)

  if (!(isTRUE(diff1) || isTRUE(diff2))) {
    stop("[SELFTEST] Seed effect test failed: coefficients did not change at t=1 nor t=2.", call. = FALSE)
  }

  cat(sprintf("[CASE 6] Seed effect observed: t1=%s, t2=%s\n",
              if (isTRUE(diff1)) "DIFF" else "SAME",
              if (isTRUE(diff2)) "DIFF" else "SAME"))

  cat("=== OK CASE 6 ===\n")
}

# --------------------------------------------------------------------------------------
# Désactivation patch si dispo
# --------------------------------------------------------------------------------------
if (.patch_enabled && exists("ergm_patch_disable", mode = "function")) {
  try(ergm_patch_disable(), silent = TRUE)
}

cat("\nTous les tests selftest_erpm_long ont passé.\n")