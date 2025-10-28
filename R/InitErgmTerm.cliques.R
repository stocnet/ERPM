#' ERGM term: cliques
#' Terme composite simple: Δ = w * ( ΔTriangles + ΔTwoStars ) autour du toggle.
#' Arguments:
#'   w (num) poids global, défaut 1.
InitErgmTerm.cliques <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cliques"

  a <- check.ErgmTerm(
    nw, arglist,
    directed = NULL, bipartite = NULL,
    varnames = c("w"),
    vartypes = c("numeric"),
    defaultvalues = list(1),
    required = c(FALSE)
  )

  w <- a$w
  if (length(w) != 1L) ergm_Init_stop(sQuote(termname), ": 'w' doit être scalaire.")

  # Un seul coefficient
  coef.names <- "cliques"

  # Inputs: passer w
  inputs <- c(w)

  # Réseau vide: pas de 2-stars ni triangles => 0
  en0 <- 0

  list(
    name = "cliques",
    coef.names = coef.names,
    inputs = inputs,
    dependence = TRUE,
    emptynwstats = en0
  )
}
