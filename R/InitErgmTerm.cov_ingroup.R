# ==============================================================================
# Fichier : InitErgmTerm.cov_ingroup.R
# Terme   : cov_ingroup(cov, size = NULL, category = NULL)
# Stat    : T = sum_g [ n_g * (sum_{i in g} x_i ) * 1[n_g in S] ]
#           avec x_i numérique, ou x_i = 1[c_i == category] si 'category' est fourni.
#           S = ensemble de tailles autorisées (size), NULL => toutes tailles.
# ==============================================================================

#' ERGM term: cov_ingroup
#'
#' \deqn{T(B) = \sum_g n_g \Big(\sum_{i\in g} x_i\Big)\, \mathbf{1}[\,n_g\in S\,],}
#' où x_i est un attribut acteur, n_g la taille du groupe, S l’ensemble des tailles retenues.
#'
#' @param cov character|numeric  Nom d'un attribut acteur (numérique ou qualitatif),
#'                               ou un vecteur numérique de longueur |A| (mode 1).
#' @param size integer|numeric|NULL  Ensemble de tailles S (ex. c(2,3,4)). NULL => toutes.
#' @param category character|NULL     Modalité ciblée si attribut qualitatif.
#' @export
InitErgmTerm.cov_ingroup <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_ingroup"

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                 "size",  "category"),
    vartypes      = c("numeric,character",   "numeric","character"),
    defaultvalues = list(NULL,               NULL,     NULL),
    required      = c(TRUE,                  FALSE,    FALSE)
  )

  # --- 0) Taille du mode acteurs (n1) -----------------------------------------
  n1 <- tryCatch(nw %n% "bipartite", error = function(e) NA_integer_)
  if (!is.numeric(n1) || !is.finite(n1) || n1 <= 0)
    ergm_Init_stop(sQuote(termname), ": réseau non biparti ou attribut %n% 'bipartite' manquant/invalide.")

  # --- 1) Construire le vecteur acteur x (longueur n1) ------------------------
  cov      <- a$cov
  category <- a$category

  get_actor_cov <- function(nw, cov, category = NULL) {
    n1 <- as.integer(nw %n% "bipartite")

    # Par convention du builder: acteurs en tête (1..n1). Fallback si nécessaire.
    ia <- seq_len(n1)
    vn <- network::network.vertex.names(nw)
    if (length(vn) >= n1) {
      ia_guess <- which(!grepl("^G\\d+$", vn))
      if (length(ia_guess) == n1) ia <- ia_guess
    }

    if (is.character(cov) && length(cov) == 1L) {
      vals <- network::get.vertex.attribute(nw, cov)
      if (is.null(vals))
        ergm_Init_stop(sQuote(termname), ": attribut inexistant: ", sQuote(cov), ".")
      x <- vals[ia]

      if (!is.null(category)) {
        xb <- as.integer(as.character(x) == category); xb[is.na(xb)] <- 0L
        return(list(x = as.double(xb), cov_label = paste0(cov, "==", category)))
      } else {
        x_num <- suppressWarnings(as.numeric(x))
        if (any(!is.finite(x_num))) {
          bad <- which(!is.finite(x_num))[1]
          ergm_Init_stop(sQuote(termname),
                         ": attribut numérique contient NA/NaN/Inf côté acteurs. ",
                         "Exemple: vertex=", vn[ia[bad]], ", valeur=", as.character(x[bad]))
        }
        return(list(x = as.double(x_num), cov_label = cov))
      }

    } else {
      x_num <- suppressWarnings(as.numeric(cov))
      if (any(!is.finite(x_num)))
        ergm_Init_stop(sQuote(termname), ": vecteur 'cov' contient NA/NaN/Inf.")
      if (length(x_num) < n1)
        ergm_Init_stop(sQuote(termname), ": longueur(cov) < |A| = ", n1, ".")
      if (!is.null(category))
        ergm_Init_stop(sQuote(termname), ": 'category' ne s'applique pas quand 'cov' est numérique direct.")
      return(list(x = as.double(x_num[seq_len(n1)]), cov_label = "cov"))
    }
  }

  ax <- get_actor_cov(nw, cov = cov, category = category)
  x <- ax$x
  cov_label <- ax$cov_label

  # --- 2) Normaliser le critère de tailles S (argument 'size') ----------------
  S <- a$size
  if (is.null(S) || length(S) == 0L) {
    sizes <- integer(0)                      # S = toutes tailles
  } else {
    S <- unique(as.integer(S))
    if (any(!is.finite(S)) || any(S < 1L))
      ergm_Init_stop(sQuote(termname), ": 'size' doit contenir des entiers >= 1.")
    sizes <- sort(S)
  }

  # --- 3) Inputs pour le C ----------------------------------------------------
  # INPUT_PARAM layout:
  #   [1]     = n1
  #   [2]     = L = length(sizes)
  #   [3..]   = sizes (L entrées, éventuellement L=0)
  #   [ .. ]  = x_1..x_n1
  L <- length(sizes)
  inputs <- c(
    as.double(n1),
    as.double(L),
    as.double(sizes),
    as.double(x)
  )

  # --- 4) Nom du coefficient ---------------------------------------------------
  pretty_sizes <- if (L == 0L) "all" else paste0("S{", paste(sizes, collapse = ","), "}")
  coef.names <- paste0("cov_ingroup[", cov_label, "]_", pretty_sizes)

  # --- 5) Spécification du terme ----------------------------------------------
  list(
    name         = "cov_ingroup",    # doit matcher C_CHANGESTAT_FN(c_cov_ingroup)
    coef.names   = coef.names,
    inputs       = inputs,           # n1, L, sizes[L], x[n1]
    dependence   = TRUE,
    emptynwstats = 0
  )
}