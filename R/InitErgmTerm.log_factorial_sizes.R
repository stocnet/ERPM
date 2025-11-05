# ==============================================================================
# Fichier : R/InitErgmTerm.log_factorial_sizes.R
# Terme   : log_factorial_sizes (sans arguments)
# Stat    : sum_{g in mode-2} lgamma(deg_g), avec f(0)=0
# Usage   : summary(nw ~ log_factorial_sizes) | ergm(nw ~ log_factorial_sizes)
#           erpm(nw ~ log_factorial_sizes)   | erpm(partition ~ log_factorial_sizes)
# ==============================================================================

#' ERGM term: log_factorial_sizes
#'
#' Somme, sur les nœuds du mode 2 (groupes), de lgamma(deg(g)).
#' Aucun argument utilisateur. Réseau strictement biparti, non orienté.
#'
#' @export
InitErgmTerm.log_factorial_sizes <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "log_factorial_sizes"

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,             # exige un biparti
    varnames      = character(0),     # aucun argument
    vartypes      = character(0),
    defaultvalues = list(),
    required      = logical(0)
  )

  # Garde-fous cohérents avec tes autres inits
  if (isTRUE(nw %n% "directed"))
    ergm_Init_stop(sQuote(termname), ": utiliser un biparti non orienté (edges entre modes).")

  coef.names   <- termname
  inputs       <- NULL
  emptynwstats <- 0                   # réseau vide -> 0

  list(
    name         = termname,          # doit matcher c_log_factorial_sizes
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,
    emptynwstats = emptynwstats
  )
}