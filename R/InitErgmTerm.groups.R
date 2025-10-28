#' ERGM term: groups
#' Compte les sommets dont le degré total (in+out) ∈ [from,to).
#' Delta-stat = somme des variations sur tail et head lors d'un toggle.
#' Arguments:
#'   from (num), to (num, Inf autorisé).
#' Retour:
#'   list(name="groups", inputs= c(rbind(from,to)), dependence=TRUE, emptynwstats=…)
InitErgmTerm.groups <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "groups"

  a <- check.ErgmTerm(
    nw, arglist,
    directed = NULL, bipartite = NULL,
    varnames = c("from", "to"),
    vartypes = c("numeric", "numeric"),
    defaultvalues = list(NULL, Inf),
    required = c(TRUE, FALSE)
  )

  from <- a$from
  to   <- a$to

  # Harmonisation des longueurs
  if (length(to) == 1L && length(from) > 1L)      to   <- rep(to,   length(from))
  else if (length(from) == 1L && length(to) > 1L) from <- rep(from, length(to))
  else if (length(from) != length(to)) {
    ergm_Init_stop("Les arguments de ", sQuote(termname),
                   " doivent avoir la même longueur ou l’un être de longueur 1.")
  }

  # to = Inf -> borne ouverte au-delà de la taille du graphe
  to <- ifelse(is.infinite(to), pmax(from, network.size(nw)) + 1, to)
  if (any(from >= to)) ergm_Init_stop(sQuote(termname), ": exigence ", sQuote("from < to"), " violée.")

  if (length(from) == 0L) return(NULL)

  # Noms des coefficients
  coef.names <- ifelse(
    to >= network.size(nw) + 1,
    paste0("groups_", from, "+"),
    paste0("groups_", from, "to", to)
  )

  # Inputs empilés par paires
  inputs <- c(rbind(from, to))

  # Stat du réseau vide: deg=0 partout
  # => pour chaque intervalle, compte = n si 0 ∈ [from,to) sinon 0.
  en0 <- integer(length(from))
  en0[ from <= 0 & 0 < to ] <- network.size(nw)

  list(
    name = "groups",
    coef.names = coef.names,
    inputs = inputs,
    dependence = TRUE,
    emptynwstats = en0
  )
}
