# ==============================================================================
# Fichier : InitErgmTerm.cliques.R
# Fonction : InitErgmTerm.cliques()
# Utilité : Déclarer le terme ERGM 'cliques' (générique k >= 2) sur biparti
# ==============================================================================

#' ERGM term: cliques
#'
#' @title Somme des k-cliques d'acteurs via les tailles de groupes (mode 2)
#' @description
#' Pour un réseau biparti {acteurs (mode 1)}–{groupes (mode 2)}, la projection
#' acteur–acteur forme, pour chaque groupe g de taille \eqn{n_g}, une clique
#' complète \eqn{K_{n_g}}. Le nombre total de k-cliques d'acteurs vaut donc
#' \eqn{\sum_g \binom{n_g}{k}}. Ce terme calcule cette statistique **sans**
#' matérialiser la projection, en ne touchant qu'au nœud de groupe impacté par
#' chaque toggle.
#'
#' @details
#' - Réseau **biparti requis**. Les tailles de groupes sont les degrés des
#'   nœuds du mode 2.
#' - Pour \code{k = 1}, la statistique est le **nombre de groupes de taille 1**
#'   (équivalent au nombre de nœuds de degré 1 dans la projection acteur–acteur).
#' - Vectorisable en \code{k}. Option \code{normalized=TRUE} pour renvoyer une
#'   proportion par rapport au maximum possible \eqn{\binom{N_1}{k}}, où
#'   \eqn{N_1} = nombre d'acteurs.
#'
#' @param k \code{numeric} entier(s) \eqn{\ge 1}. Taille(s) de clique(s) d'acteurs.
#' @param normalized \code{logical} scalaire. Si \code{TRUE}, renvoie la proportion
#'   \eqn{\frac{\sum_g \binom{n_g}{k}}{\binom{N_1}{k}}}. Défaut : \code{FALSE}.
#'
#' @return Liste d'initialisation pour {ergm}.
#'
#' @examples
#' \dontrun{
#'   # Paires de co-membres (k=2)
#'   summary(nw ~ cliques(k=2) + b1part)
#'
#'   # Triangles de co-membres (k=3), normalisé
#'   summary(nw ~ cliques(k=3, normalized=TRUE) + b1part)
#'
#'   # Plusieurs k en même temps
#'   summary(nw ~ cliques(k=c(2,3,4)) + b1part)
#' }
#' @export
InitErgmTerm.cliques <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cliques"
  # (alias clique_size -> k inchangé)

  # --- 0) Normalisation des arguments utilisateurs ---------------------------
  # Supporter :
  #   - cliques(1)             -> k = 1
  #   - cliques(clique_size=1) -> k = 1
  #   - cliques(k=1)           -> k = 1
  if (length(arglist) == 1L) {
    nm <- names(arglist)
    if (is.null(nm) || isTRUE(nm[1L] == "")) {
      # argument positionnel unique : cliques(1)
      arglist <- list(k = arglist[[1L]])
    }
  }
  if (!is.null(names(arglist)) && "clique_size" %in% names(arglist)) {
    arglist[["k"]] <- arglist[["clique_size"]]
    arglist[["clique_size"]] <- NULL
  }

  # --- 1) Récupération robuste de N1 (mode 1) ---
  n1 <- network::get.network.attribute(nw, "bipartite")
  if (is.null(n1) || is.na(n1)) {
    ergm_Init_stop(sQuote(termname), ": réseau non biparti ou attribut 'bipartite' manquant.")
  }
  n1 <- as.integer(n1)
  if (n1 <= 1L) ergm_Init_stop(sQuote(termname), ": biparti invalide (N1 <= 1).")

  # --- 2) Spécification des args attendus ---
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("k", "normalized"),
    vartypes      = c("numeric", "logical"),
    defaultvalues = list(2, FALSE),
    required      = c(FALSE, FALSE)
  )

  k  <- a$k
  nz <- a$normalized

  # --- 3) Validation ---   # <- k >= 1 désormais
  if (length(k) < 1L) ergm_Init_stop(sQuote(termname), ": spécifiez au moins une valeur de k.")
  if (any(k != as.integer(k) | k < 1)) ergm_Init_stop(sQuote(termname), ": 'k' doit contenir des entiers >= 1.")
  if (length(nz) != 1L || is.na(nz)) ergm_Init_stop(sQuote(termname), ": 'normalized' doit être un booléen scalaire.")

  # --- 4) Échelle (normalisation optionnelle) ---
  scale <- rep(1, length(k))
  if (isTRUE(nz)) {
    denom <- choose(n1, k)              # marche aussi pour k=1 (denom = N1)
    denom[denom == 0 | is.na(denom)] <- Inf
    scale <- 1 / denom
  }

  # --- 5) Noms + inputs C (k, scale) ---
  coef.names <- if (isTRUE(nz)) paste0("cliques_k", k, "_norm") else paste0("cliques_k", k)
  inputs <- c(rbind(as.integer(k), as.double(scale)))

  list(
    name         = "cliques",
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,
    emptynwstats = numeric(length(k))
  )
}