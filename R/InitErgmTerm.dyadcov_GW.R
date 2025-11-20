# ==============================================================================
# Fichier : R/InitErgmTerm.dyadcov_GW.R
# Terme   : dyadcov_GW
# Stat    :
#   T_GW(p;Z, lambda)
#     = sum_g sum_{k=2..n_g} a_k(lambda) * S_g^{(k)}(Z)
#     = sum_g sum_{k=2..n_g} (-1/lambda)^{k-2}
#                      * sum_{C∈C_k(g)} prod_{i<j∈C} z_{ij}
#
# INPUT_PARAM : c(n1, lambda, Z[n1*n1])  (Z en ordre colonne-major)
# Debug      : options(ERPM.dyadcov_GW.debug = TRUE) pour activer les logs
# ==============================================================================

InitErgmTerm.dyadcov_GW <- function(nw, arglist, ...) {
  termname <- "dyadcov_GW"

  # -- debug helpers ------------------------------------------------------------
  dbg    <- isTRUE(getOption("ERPM.dyadcov_GW.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[dyadcov_GW][DEBUG]", ..., "\n", sep = "")

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("dyadcov",              "lambda"),
    vartypes      = c("matrix,character",     "numeric"),
    defaultvalues = list(NULL,                2),
    required      = c(TRUE,                   FALSE)
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

  # ----- coercition numérique + FAIL-FAST NA -----------------------------------
  if (!is.numeric(dyad_mat))
    stop(termname, ": la matrice dyadique doit être numérique.")

  if (anyNA(dyad_mat))
    stop(termname, ": NA non autorisé dans la matrice dyadique.")

  # ----- contrôle optionnel de symétrie ----------------------------------------
  tol <- 1e-8
  if (max(abs(dyad_mat - t(dyad_mat))) > tol) {
    stop(termname, ": la matrice dyadique doit être symétrique (différence > tolérance).")
  }

  dbgcat("dyadcov dim = ", paste(dim(dyad_mat), collapse = "x"),
         " | sample = ",
         paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L), collapse = ","))

  # ----- lambda ---------------------------------------------------------------
  lambda_raw <- a$lambda
  if (is.null(lambda_raw) || length(lambda_raw) == 0L) {
    lambda <- 2
  } else {
    if (!is.numeric(lambda_raw))
      stop(termname, ": 'lambda' doit être numérique.")
    lambda <- as.numeric(lambda_raw[1L])
  }

  if (!is.finite(lambda) || lambda <= 0)
    stop(termname, ": 'lambda' doit être un réel strictement positif (et typiquement > 1).")

  dbgcat("lambda = ", format(lambda, digits = 6L))

  # ----- nom de coefficient ---------------------------------------------------
  # On encode lambda de façon compacte dans le nom
  lambda_tag <- gsub("[^0-9\\.eE\\-]+", "_", format(lambda, digits = 4L))
  base_label <- sprintf("dyadcov_GW[%s]_lambda%s", dyad_label, lambda_tag)
  coef.name  <- base_label
  dbgcat("coef.name = ", coef.name)

  # ----- construire INPUT_PARAM -----------------------------------------------
  # as.double(matrice) => ordre colonne-major compatible avec l'indexation C
  inputs <- c(
    as.double(n1),
    as.double(lambda),
    as.double(dyad_mat)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " lambda=", format(lambda, digits = 6L),
         " | Z[1:6]=", paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L),
                             collapse = ","))

  list(
    name         = "dyadcov_GW",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, lambda, Z[n1*n1]
    dependence   = TRUE,
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}
