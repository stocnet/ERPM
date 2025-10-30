# ==============================================================================
# Fichier : InitErgmTerm.squared_sizes.R
# Fonction : InitErgmTerm.squared_sizes()
# Utilité : Déclarer le terme ERGM 'squared_sizes' (mode 2 uniquement)
# ==============================================================================

#' ERGM term: squared_sizes
#'
#' @title Somme des puissances des tailles de groupes
#' @description Calcule la somme, sur les nœuds du **mode 2** (groupes) uniquement,
#' de \code{deg^pow} pour les nœuds dont le degré \code{deg} appartient à l'intervalle
#' \code{[from, to)}. Par défaut, tous les groupes sont inclus avec \code{pow = 2}.
#'
#' @details
#' Le réseau doit être biparti \{mode1=acteurs\}–\{mode2=groupes\}. La taille d'un groupe
#' est le degré du nœud correspondant dans le mode 2. Le terme est vectorisable :
#' on peut fournir plusieurs intervalles \code{from/to} et puissances \code{pow}.
#'
#' @param from \code{numeric} entier(s) ≥ 0, borne(s) inférieure(s) incluse(s). Défaut : \code{1}.
#' @param to \code{numeric} entier(s), borne(s) supérieure(s) exclusive(s). \code{Inf} autorisé. Défaut : \code{Inf}.
#' @param pow \code{numeric} entier(s) ≥ 1, puissance(s) appliquée(s) au degré. Défaut : \code{2}.
#'
#' @return Liste d'initialisation pour {ergm}.
#'
#' @examples
#' \dontrun{
#'   # Tous les groupes, somme des carrés des tailles
#'   summary(nw ~ squared_sizes + b1part)
#'
#'   # Groupes de taille 2 à 4, puissance 3
#'   summary(nw ~ squared_sizes(from=2, to=5, pow=3) + b1part)
#' }
#' @export
InitErgmTerm.squared_sizes <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "squared_sizes"

  # -- Validation base et exigences de structure --------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed   = NULL,
    bipartite  = TRUE,                      # exige un réseau biparti
    varnames   =    c(    "from",     "to",       "pow"),
    vartypes   =    c(    "numeric",  "numeric",  "numeric"),
    defaultvalues = list( 1,          Inf,        2),        # par défaut : tous les groupes
    required   =    c(    FALSE,      FALSE,      FALSE)
  )

  from <- a$from; to <- a$to; pow <- a$pow

  # -- Si la taille des vecteurs ne sont pas les mêmes -----------------------------------------
  if (length(to)   == 1L && length(from) >  1L) to   <- rep(to,   length(from))  # Crée un vecteur to à partir du scalaire : exemple d'appel squared_sizes(from=c(1,2,4), to=5, pow=2) -> squared_sizes(from=c(1,2,4), to=c(5,5,5), pow=2)
  if (length(from) == 1L && length(to)   >  1L) from <- rep(from, length(to))    # Crée un vecteur from à partir du scalaire : exemple d'appel squared_sizes(from=c(1), to=c(5,2,4), pow=2) -> squared_sizes(from=c(1,1,1), to=c(5,2,4), pow=2)
  if (length(pow)  == 1L && length(from) >  1L) pow  <- rep(pow,  length(from))  # Crée un vecteur pow à partir du scalaire : exemple d'appel squared_sizes(from=c(1,2,4), to=c(5,2,4), pow=2) -> squared_sizes(from=c(1,2,4), to=c(5,2,4), pow=c(2,2,2)

  if (length(from) != length(to) || length(pow) != length(from))
    ergm_Init_stop(sQuote(termname), ": longueurs incompatibles pour from/to/pow.")

  # -- Contraintes sur pow -------------------------------------------------------
  if (any(pow < 1 | pow != as.integer(pow)))
    ergm_Init_stop(sQuote(termname), ": 'pow' doit être un entier >= 1.")

  # -- Borne supérieure : remplacer Inf par (n+1) pour respecter [from, to) -----
  to <- ifelse(is.infinite(to), pmax(from, network.size(nw)) + 1, to)
  # n1 <- network.bipartite(nw)  # nb de sommets du mode 1 (acteurs)
  # to <- ifelse(is.infinite(to), pmax(from, n1) + 1, to)

  # -- Cohérence des intervalles -------------------------------------------------
  if (any(from >= to))
    ergm_Init_stop(sQuote(termname), ": from < to requis.")

  # -- Cas trivial : aucun intervalle -------------------------------------------
  if (length(from) == 0L) return(NULL)

  # -- Noms de coefficients lisibles --------------------------------------------
  coef.names <- ifelse(
    to >= network.size(nw) + 1,
    paste0("squared_sizes_", from, "+", ifelse(pow != 1, paste0("_pow", pow), "")),
    paste0("squared_sizes_", from, "to", to, ifelse(pow != 1, paste0("_pow", pow), ""))
  )

  # -- Entrées compactées pour le C : (from, to, pow) en colonnes puis aplati ----
  inputs <- c(rbind(as.integer(from), as.integer(to), as.integer(pow)))

  # -- Retour standard d'initialisation -----------------------------------------
  list(
    name         = "squared_sizes",   # doit matcher CHANGESTAT_FN(c_squared_sizes) mais sans le c_ qui est parse par ergm
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,              # dépend des toggles
    emptynwstats = numeric(length(from))  # réseau vide -> 0
  )
}