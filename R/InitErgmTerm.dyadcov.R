# ==============================================================================
# Fichier : R/InitErgmTerm.dyadcov.R
# Terme   : dyadcov
# Stat    :
#   - normalized = FALSE :
#       T^{(k)}(p;Z) = sum_g sum_{C∈C_k(g)} prod_{i<j∈C} z_{ij}
#   - normalized = TRUE :
#       T^{(k)}_norm(p;Z) = sum_g 1[n_g>=k] * (1/choose(n_g,k)) *
#                           sum_{C∈C_k(g)} prod_{i<j∈C} z_{ij}
#
# INPUT_PARAM : c(n1, k, normalized, Z[n1*n1])  (Z en ordre colonne-major)
# Debug      : options(ERPM.dyadcov.debug = TRUE) pour activer les logs
# ==============================================================================

InitErgmTerm.dyadcov <- function(nw, arglist, ...) {
  termname <- "dyadcov"

  # -- debug helpers ------------------------------------------------------------
  dbg    <- isTRUE(getOption("ERPM.dyadcov.debug", TRUE))
  dbgcat <- function(...) if (dbg) cat("[dyadcov][DEBUG]", ..., "\n", sep = "")

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("dyadcov",                 "clique_size", "normalized"),
    vartypes      = c("matrix,character",        "numeric,integer", "logical"),
    defaultvalues = list(NULL,                   2,               FALSE),
    required      = c(TRUE,                      FALSE,           FALSE)
  )

  # ----- n1 : taille du mode acteurs -----
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L)
    stop(termname, ": réseau non biparti strict (attribut %n% 'bipartite' manquant ou invalide).")
  dbgcat("n1 = ", n1)

  # ----- récupérer la matrice dyadique (n1 x n1) -----
  dyad_raw   <- a$dyadcov
  dyad_label <- NULL

  if (is.character(dyad_raw) && length(dyad_raw) == 1L) {
    # cas : nom d'attribut %n% portant la matrice
    dyad_mat <- nw %n% dyad_raw
    dyad_label <- dyad_raw
    if (is.null(dyad_mat))
      stop(termname, ": attribut de niveau réseau inexistant: ", sQuote(dyad_raw), ".")
    dbgcat("dyadcov source = network attribute ", sQuote(dyad_label))
  } else {
    # cas : matrice passée littéralement
    dyad_mat   <- dyad_raw
    dyad_label <- "dyadcov"
    dbgcat("dyadcov source = matrix literal")
  }

  if (!is.matrix(dyad_mat))
    stop(termname, ": 'dyadcov' doit être une matrice ou le nom d'un attribut de niveau réseau.")

  nr <- nrow(dyad_mat)
  nc <- ncol(dyad_mat)

  if (nr < n1 || nc < n1)
    stop(termname, ": dimensions de la matrice dyadique (", nr, "x", nc,
         ") insuffisantes pour n1 = ", n1, ".")

  # On restreint, si besoin, au bloc supérieur gauche n1 x n1
  if (nr > n1 || nc > n1) {
    dyad_mat <- dyad_mat[seq_len(n1), seq_len(n1), drop = FALSE]
    dbgcat("dyadcov truncated to ", n1, "x", n1)
  }

  # ----- coercition numérique + FAIL-FAST NA -----
  if (!is.numeric(dyad_mat))
    stop(termname, ": la matrice dyadique doit être numérique.")

  if (anyNA(dyad_mat))
    stop(termname, ": NA non autorisé dans la matrice dyadique.")

  # ----- contrôle optionnel de symétrie -----
  tol <- 1e-8
  if (max(abs(dyad_mat - t(dyad_mat))) > tol) {
    stop(termname, ": la matrice dyadique doit être symétrique (différence > tolérance).")
  }

  dbgcat("dyadcov dim = ", paste(dim(dyad_mat), collapse = "x"),
         " | sample = ",
         paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L), collapse = ","))

  # ----- clique_size (k) -----
  k_raw <- a$clique_size
  if (is.null(k_raw) || length(k_raw) == 0L) {
    k <- 2L
  } else {
    if (!is.numeric(k_raw))
      stop(termname, ": 'clique_size' doit être numérique ou entier.")
    k <- as.integer(round(k_raw[1L]))
  }

  if (!is.finite(k) || k < 2L)
    stop(termname, ": 'clique_size' doit être un entier >= 2.")

  dbgcat("clique_size (k) = ", k)

  # ----- normalized (booléen) -----
  norm_raw <- a$normalized
  if (is.null(norm_raw) || length(norm_raw) == 0L) {
    normalized <- FALSE
  } else {
    normalized <- isTRUE(norm_raw[1L])
  }
  norm_flag <- if (normalized) 1 else 0
  dbgcat("normalized = ", normalized)

  # ----- nom de coefficient -----
  base_label <- sprintf("dyadcov[%s]_k%d", dyad_label, k)
  coef.name  <- if (normalized) paste0(base_label, "_norm") else base_label
  dbgcat("coef.name = ", coef.name)

  # ----- construire INPUT_PARAM -----
  # as.double(matrice) => ordre colonne-major compatible avec l'indexation C
  inputs <- c(
    as.double(n1),
    as.double(k),
    as.double(norm_flag),
    as.double(dyad_mat)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " k=", k, " normalized=", norm_flag,
         " | Z[1:6]=", paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L),
                             collapse = ","))

  list(
    name         = "dyadcov",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, k, normalized, Z[n1*n1]
    dependence   = TRUE,
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}
