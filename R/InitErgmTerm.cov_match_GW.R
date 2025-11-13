# ==============================================================================
# Fichier : R/InitErgmTerm.cov_match_GW.R
# Terme   : cov_match_GW(cov, lambda = 2, category = NULL,
#                        normalized = c("none","by_group","global"))
# Stat    : non normalisée :
#             S_GW = sum_g sum_r  λ * (1 - r_λ^{ n_{g,r} })
#           ciblée :
#             S_GW^{(κ)} = sum_g  λ * (1 - r_λ^{ n_{g,κ} })
#           by_group :
#             sum_g [ Num(g) / Den(g) ] avec
#               Num(g)=sum_r λ(1-r_λ^{n_{g,r}})  (ou λ(1-r_λ^{n_{g,κ}}) ciblé)
#               Den(g)=λ(1-r_λ^{n_g})
#           global :
#             [ sum_g Num(g) ] / [ λ(1-r_λ^{N1}) ]
# ==============================================================================

#' ERGM term: cov_match_GW
#' @param cov character Nom d'un attribut acteur (factor/character).
#' @param lambda numeric>1  Scalaire ou vecteur. Poids géométrique (r_λ=(λ-1)/λ).
#' @param category character|NULL Modalité ciblée (facultatif).
#' @param normalized character|logical "none"|"by_group"|"global".
#' @export
InitErgmTerm.cov_match_GW <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_match_GW"

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov","lambda","category","normalized"),
    vartypes      = c("character","numeric","character","logical,character"),
    defaultvalues = list(NULL,        2,        NULL,       "none"),
    required      = c(TRUE,           FALSE,    FALSE,      FALSE)
  )

  # -- biparti strict
  n1 <- tryCatch(nw %n% "bipartite", error = function(e) NA_integer_)
  if (!is.numeric(n1) || !is.finite(n1) || n1 <= 0)
    ergm_Init_stop(sQuote(termname), ": réseau non biparti ou attribut %n% 'bipartite' manquant/invalide.")

  # -- normalized
  normalized <- a$normalized
  if (is.logical(normalized)) normalized <- if (isTRUE(normalized)) "by_group" else "none"
  normalized <- match.arg(tolower(as.character(normalized)), c("none","by_group","global"))
  norm_mode  <- switch(normalized, none=0L, by_group=1L, global=2L)

  # -- lambda
  lambdas <- as.double(a$lambda)
  if (!length(lambdas)) lambdas <- 2
  if (any(!is.finite(lambdas)) || any(lambdas <= 1))
    ergm_Init_stop(sQuote(termname), ": 'lambda' doit être > 1 (numérique, fini).")
  lambdas <- as.double(unique(lambdas))
  K <- length(lambdas)

  # -- attribut catégoriel et catégorie ciblée
  covname   <- a$cov
  if (!(is.character(covname) && length(covname)==1L))
    ergm_Init_stop(sQuote(termname), ": 'cov' doit être le nom d'un attribut acteur (factor/character).")

  ia <- seq_len(n1)
  vals <- network::get.vertex.attribute(nw, covname)
  if (is.null(vals))
    ergm_Init_stop(sQuote(termname), ": attribut inexistant: ", sQuote(covname), ".")

  f <- as.factor(vals[ia])
  category <- a$category
  if (!is.null(category) && !(category %in% levels(f)))
    levels(f) <- c(levels(f), category)

  z <- as.integer(f); z[!is.finite(z)] <- 0L  # 0 = NA/absent
  kappa_code <- if (is.null(category)) 0L else as.integer(match(category, levels(f)))
  has_kappa  <- as.double(as.integer(kappa_code > 0))
  cov_label  <- if (is.null(category)) covname else paste0(covname,"==",category)

  # -- inputs packés pour le C
  #   [0] = n1
  #   [1] = K
  #   [2] = norm_mode
  #   [3] = has_kappa
  #   [4] = kappa_code
  #   [5..(5+K-1)]       = lambdas
  #   [5+K .. 5+K+n1-1]  = z[1..n1]
  inputs <- c(
    as.double(n1),
    as.double(K),
    as.double(norm_mode),
    has_kappa,
    as.double(kappa_code),
    lambdas,
    as.double(z)
  )

  # -- noms des coefficients
  suffix_norm <- switch(normalized, none="", by_group="_bygrp", global="_glob")
  fmt_lambda  <- function(x) sub("\\.?0+$", "", format(x, trim=TRUE))
  coef.names  <- paste0("cov_match_GW[", cov_label, "]_l", vapply(lambdas, fmt_lambda, ""), suffix_norm)

  list(
    name         = "cov_match_GW",   # doit matcher C_CHANGESTAT_FN(c_cov_match_GW)
    coef.names   = coef.names,       # longueur = K
    inputs       = inputs,
    dependence   = TRUE,
    emptynwstats = 0
  )
}
