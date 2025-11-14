# ==============================================================================
# Fichier : R/InitErgmTerm.cov_diff_GW.R
# Terme   : cov_diff_GW
# Stat    : T_GW = sum_{k>=2} (-1/lambda)^{k-1} c_k(cov,p)
#           avec c_k = sum_g sum_{S ⊂ g, |S|=k} (max_{i∈S} x_i - min_{i∈S} x_i)
# Règle   : NA interdits dans la covariée (fail-fast)
# INPUT_PARAM : c(n1, L, lambda[1..L], x[1..n1])
# Debug   : options(ERPM.cov_diff_GW.debug = TRUE) pour activer les logs
# ==============================================================================

InitErgmTerm.cov_diff_GW <- function(nw, arglist, ...) {
  termname <- "cov_diff_GW"

  # -- debug helpers ------------------------------------------------------------
  dbg    <- isTRUE(getOption("ERPM.cov_diff_GW.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[cov_diff_GW][DEBUG]", ..., "\n", sep = "")

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                              "lambda"),
    vartypes      = c("character,numeric,logical,vector", "numeric,vector"),
    defaultvalues = list(NULL,                             2),
    required      = c(TRUE,                                FALSE)
  )

  # ----- n1 : taille du mode acteurs -----
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) stop(termname, ": réseau non biparti strict.")
  dbgcat("n1 = ", n1)

  # ----- récupérer le vecteur covarié acteurs (longueur >= n1) -----
  cov_raw <- a$cov
  if (is.character(cov_raw) && length(cov_raw) == 1L) {
    cov_vec   <- network::get.vertex.attribute(nw, cov_raw)
    cov_label <- cov_raw
    if (is.null(cov_vec))
      stop(termname, ": attribut inexistant: ", sQuote(cov_raw), ".")
    dbgcat("cov source = vertex attribute ", sQuote(cov_label))
  } else {
    cov_vec   <- cov_raw
    cov_label <- "cov"
    dbgcat("cov source = vector literal")
  }

  if (length(cov_vec) < n1)
    stop(termname, ": longueur de la covariée < n1.")

  cov_vec <- cov_vec[seq_len(n1)]
  dbgcat("cov length = ", length(cov_vec),
         " | head = ", paste(utils::head(as.character(cov_vec), 6L), collapse = ","))

  # ----- contrainte numérique + FAIL-FAST NA -----
  if (is.logical(cov_vec) || is.integer(cov_vec)) {
    cov_vec <- as.numeric(cov_vec)
  }
  if (!is.numeric(cov_vec))
    stop(termname, ": la covariée doit être convertissable en numérique.")

  if (anyNA(cov_vec))
    stop(termname, ": NA non autorisé dans la covariée du mode acteurs.")

  # ----- lambda : numérique, lambda > 1, éventuellement vectorisé ---------------
  lambda_raw <- a$lambda
  if (!is.numeric(lambda_raw) || length(lambda_raw) < 1L)
    stop(termname, ": 'lambda' doit être numérique (scalaire ou vecteur).")

  lambda_vec <- as.double(lambda_raw)
  if (any(!is.finite(lambda_vec)))
    stop(termname, ": 'lambda' doit être fini.")
  if (any(lambda_vec <= 1))
    stop(termname, ": toutes les valeurs de 'lambda' doivent être > 1.")

  L <- length(lambda_vec)
  dbgcat("lambda_vec = {", paste(signif(lambda_vec, 5L), collapse = ", "), "} (L = ", L, ")")

  # ----- noms de coefficients ---------------------------------------------------
  coef.names <- sprintf(
    "cov_diff_GW[%s]_lambda%.5g",
    cov_label,
    lambda_vec
  )
  dbgcat("coef.names = ", paste(coef.names, collapse = " | "))

  # ----- construire INPUT_PARAM ------------------------------------------------
  # Layout : [0] = n1, [1] = L, [2..(1+L)] = lambda[1..L], [2+L..] = x[1..n1]
  inputs <- c(
    as.double(n1),
    as.double(L),
    as.double(lambda_vec),
    as.double(cov_vec)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L,
         " | lambda[1:6]=", paste(utils::head(signif(lambda_vec, 5L), 6L), collapse = ","),
         " | cov[1:6]=",    paste(utils::head(signif(cov_vec,    5L), 6L), collapse = ","))

  list(
    name         = "cov_diff_GW",
    coef.names   = coef.names,
    inputs       = inputs,              # n1, L, lambda[L], x[n1]
    dependence   = TRUE,
    minval       = -Inf,                # combinaison alternée -> peut être négative
    maxval       = Inf,
    emptynwstats = numeric(L)
  )
}
