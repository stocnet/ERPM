#' ERGM term: squared_sizes
#' Somme, sur les sommets dont le degré total ∈ [from,to), de (deg^pow).
#' Delta-stat = somme des variations sur tail et head lors d'un toggle.
#' Arguments:
#'   from (num), to (num, Inf autorisé), pow (num entier >=1 ; vectorisable).
InitErgmTerm.squared_sizes <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "squared_sizes"

  a <- check.ErgmTerm(
    nw, arglist,
    directed = NULL, bipartite = NULL,   # OK pour biparti aussi, terme lit deg total
    varnames = c("from", "to", "pow"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, Inf, 2),
    required = c(TRUE, FALSE, FALSE)
  )

  from <- a$from
  to   <- a$to
  pow  <- a$pow

  # Harmonisation des longueurs: recycle pow si besoin
  if (length(to) == 1L && length(from) > 1L)      to   <- rep(to,   length(from))
  else if (length(from) == 1L && length(to) > 1L) from <- rep(from, length(to))
  else if (length(from) != length(to)) {
    ergm_Init_stop("Les arguments de ", sQuote(termname),
                   " doivent avoir la même longueur ou l’un être de longueur 1.")
  }
  if (length(pow) == 1L && length(from) > 1L) pow <- rep(pow, length(from))
  if (length(pow) != length(from)) {
    ergm_Init_stop(sQuote(termname), ": la longueur de 'pow' doit égaler celle de 'from'/'to' ou être 1.")
  }

  to <- ifelse(is.infinite(to), pmax(from, network.size(nw)) + 1, to)
  if (any(from >= to)) ergm_Init_stop(sQuote(termname), ": exigence ", sQuote("from < to"), " violée.")
  if (any(pow < 1 | pow != as.integer(pow))) ergm_Init_stop(sQuote(termname), ": 'pow' doit être entier >= 1.")

  if (length(from) == 0L) return(NULL)

  coef.names <- ifelse(
    to >= network.size(nw) + 1,
    paste0("squared_sizes_", from, "+", ifelse(pow != 1, paste0("_pow", pow), "")),
    paste0("squared_sizes_", from, "to", to, ifelse(pow != 1, paste0("_pow", pow), ""))
  )

  inputs <- c(rbind(from, to, pow))

  # Réseau vide: deg=0, donc contribution = 0^pow = 0 si 0 ∈ [from,to), sinon 0.
  # => toujours 0.
  en0 <- numeric(length(from))

  list(
    name = "squared_sizes",
    coef.names = coef.names,
    inputs = inputs,
    dependence = TRUE,
    emptynwstats = en0
  )
}
