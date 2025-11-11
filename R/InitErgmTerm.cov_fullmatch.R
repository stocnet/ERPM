# ==============================================================================
# Fichier : R/InitErgmTerm.cov_fullmatch.R
# Terme   : cov_fullmatch
# Règle   : NA interdits dans la covariée (fail-fast)
# INPUT_PARAM : c(n1, L, sizes[L], K, target, cats[n1])
# ==============================================================================

InitErgmTerm.cov_fullmatch <- function(nw, arglist, ...) {
  termname <- "cov_fullmatch"

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                              "size",             "category"),
    vartypes      = c("character,numeric,logical,vector", "numeric,integer",  "character,numeric,logical"),
    defaultvalues = list(NULL,                        NULL,              NULL),
    required      = c(TRUE,                           FALSE,             FALSE)
  )

  # if (!is.null(arglist) && "size" %in% names(arglist)) {
  #   sz_expr <- arglist[["size"]]
  #   if (is.language(sz_expr) && is.call(sz_expr)) {
  #     if (identical(sz_expr[[1L]], as.name("c")) && length(sz_expr) == 1L) {
  #       stop("cov_fullmatch: 'size' vide (c()) interdit. Utilisez NULL pour toutes tailles.")
  #     }
  #   }
  # }

  # ----- n1 : taille du mode acteurs -----
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) stop(termname, ": réseau non biparti strict.")

  # ----- récupérer le vecteur covarié acteurs (longueur >= n1) -----
  cov_raw <- a$cov
  if (is.character(cov_raw) && length(cov_raw) == 1L) {
    cov_vec   <- network::get.vertex.attribute(nw, cov_raw)
    cov_label <- cov_raw
    if (is.null(cov_vec)) stop(termname, ": attribut inexistant: ", sQuote(cov_raw), ".")
  } else {
    cov_vec   <- cov_raw
    cov_label <- "cov"
  }
  if (length(cov_vec) < n1) stop(termname, ": longueur de la covariée < n1.")
  cov_vec <- cov_vec[seq_len(n1)]

  # ----- FAIL-FAST NA -----
  if (anyNA(cov_vec)) stop(termname, ": NA non autorisé dans la covariée du mode acteurs.")

  # ----- normalisation en catégories entières 1..K -----
  if (is.logical(cov_vec)) cov_vec <- as.integer(cov_vec)
  f        <- factor(cov_vec)           # sans NA à ce stade
  levels_f <- levels(f)
  K        <- length(levels_f)
  if (K == 0L) stop(termname, ": aucune modalité valide trouvée.")
  cats     <- as.integer(f)             # 1..K

  # ----- category ciblée -> target (0 si non fourni) -----
  target <- 0L
  if (!is.null(a$category)) {
    cat_val <- a$category
    # harmoniser le type pour matcher levels
    if (is.logical(cat_val)) cat_val <- as.integer(cat_val)
    ix <- match(as.character(cat_val), levels_f, nomatch = 0L)
    if (ix == 0L) {
      warning(termname, ": 'category' non trouvée dans les modalités; cible ignorée.")
    } else {
      target <- as.integer(ix)
      cov_label <- paste0(cov_label, "==", levels_f[target])
    }
  }

  # ----- filtre de tailles S -----
  sizes <- a$size
  print(sizes)
  if (is.null(sizes)) {
    # NULL => toutes tailles (comportement historique)
    L <- 0L
    sizes_vec <- numeric(0)
  } else {
    if (!is.numeric(sizes))
      stop("cov_fullmatch: 'size' doit être numérique.")
    if (length(sizes) == 0L)
      stop("cov_fullmatch: 'size' vide (integer(0)) interdit. Utilisez NULL pour toutes tailles.")
    sizes <- as.integer(round(sizes))
    if (any(!is.finite(sizes)) || any(sizes <= 0L))
      stop("cov_fullmatch: 'size' doit contenir des entiers positifs.")
    L <- length(sizes)
    sizes_vec <- as.double(sizes)
  }
  # sizes <- a$size
  # if (is.null(sizes) || length(sizes) == 0L) {
  #   L <- 0L
  #   sizes_vec <- numeric(0)
  # } else {
  #   if (!is.numeric(sizes)) stop(termname, ": 'size' doit être numérique.")
  #   sizes <- as.integer(round(sizes))
  #   if (any(sizes <= 0L)) stop(termname, ": 'size' doit contenir des entiers positifs.")
  #   L <- length(sizes)
  #   sizes_vec <- as.double(sizes)
  # }

  # ----- nom de coefficient (rétrocompatible) -----
  size_label <- if (L > 0L) sprintf("_S{%s}", paste(sizes, collapse = ",")) else "_all"
  coef.name  <- sprintf("cov_fullmatch[%s]%s", cov_label, size_label)

  # ----- construire INPUT_PARAM -----
  inputs <- c(
    as.double(n1),
    as.double(L),
    sizes_vec,
    as.double(K),
    as.double(target),
    as.double(cats)
  )

  list(
    name         = "cov_fullmatch",
    coef.names   = coef.name,
    inputs       = inputs,
    dependence   = TRUE,
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}

# # ==============================================================================
# # Fichier : InitErgmTerm.cov_fullmatch.R
# # Terme   : cov_fullmatch(cov, size = NULL, category = NULL)
# # Stat    : T = sum_g 1[ n_g in S ] * 1[ exists r in L : n_{g,r} = n_g ]
# #           Variante ciblée (category = kappa) :
# #           T_kappa = sum_g 1[ n_g in S ] * 1[ n_{g,kappa} = n_g ]
# #           S = ensemble de tailles autorisées (size), NULL => toutes tailles.
# # ==============================================================================

# #' ERGM term: cov_fullmatch
# #'
# #' \deqn{T(B;c,S) = \sum_g \mathbf{1}[\,n_g\in S\,] \,\mathbf{1}\big[\exists r\in L : n_{g,r}=n_g\big].}
# #' \deqn{T_\kappa(B;c,S) = \sum_g \mathbf{1}[\,n_g\in S\,] \,\mathbf{1}\big[n_{g,\kappa}=n_g\big].}
# #'
# #' @param cov character|vector  Nom d'un attribut acteur (qualitatif),
# #'                              ou un vecteur (numérique/texte) de longueur |A| à traiter comme catégoriel.
# #' @param size integer|numeric|NULL  Ensemble de tailles S (ex. c(2,3,4)). NULL => toutes tailles.
# #' @param category character|NULL     Modalité ciblée si on veut la version ciblée.
# #' @export
# InitErgmTerm.cov_fullmatch <- function(nw, arglist, ..., version = packageVersion("ergm")) {
#   termname <- "cov_fullmatch"

#   a <- check.ErgmTerm(
#     nw, arglist,
#     directed      = NULL,
#     bipartite     = TRUE,
#     varnames      = c("cov",               "size",  "category"),
#     vartypes      = c("character,numeric", "numeric","character"),
#     defaultvalues = list(NULL,             NULL,     NULL),
#     required      = c(TRUE,                FALSE,    FALSE)
#   )

#   # --- 0) Taille du mode acteurs (n1) -----------------------------------------
#   n1 <- tryCatch(nw %n% "bipartite", error = function(e) NA_integer_)
#   if (!is.numeric(n1) || !is.finite(n1) || n1 <= 0)
#     ergm_Init_stop(sQuote(termname), ": réseau non biparti ou attribut %n% 'bipartite' manquant/invalide.")

#   # Acteurs supposés en tête 1..n1, avec fallback sur les noms si besoin
#   ia <- seq_len(n1)
#   vn <- network::network.vertex.names(nw)
#   if (length(vn) >= n1) {
#     ia_guess <- which(!grepl("^G\\d+$", vn))
#     if (length(ia_guess) == n1) ia <- ia_guess
#   }

#   cov      <- a$cov
#   category <- a$category

#   # --- 1) Construire le vecteur des catégories acteurs -> entiers 1..K --------
#   get_actor_cats <- function(nw, cov, category = NULL, n1, ia) {
#     # Récupération de l'attribut ou d'un vecteur fourni
#     if (is.character(cov) && length(cov) == 1L) {
#       vals <- network::get.vertex.attribute(nw, cov)
#       if (is.null(vals))
#         ergm_Init_stop(sQuote(termname), ": attribut inexistant: ", sQuote(cov), ".")
#       v <- vals[ia]
#       cov_label <- cov
#     } else {
#       v <- cov
#       if (length(v) < n1)
#         ergm_Init_stop(sQuote(termname), ": longueur(cov) < |A| = ", n1, ".")
#       v <- v[seq_len(n1)]
#       cov_label <- "cov"
#     }

#     # Transformer en facteur catégoriel
#     f <- as.factor(as.character(v))
#     # NA explicites -> modalité 0 côté C
#     cats <- as.integer(f)
#     cats[is.na(cats)] <- 0L

#     levels_f <- levels(f)
#     K <- length(levels_f)

#     if (!is.null(category)) {
#       kappa <- match(category, levels_f, nomatch = NA_integer_)
#       if (is.na(kappa))
#         ergm_Init_stop(sQuote(termname), ": 'category' inconnue pour ", sQuote(cov_label), ".")
#       cov_label <- paste0(cov_label, "==", category)
#       target <- as.integer(kappa)
#     } else {
#       target <- 0L
#     }

#     list(cats = as.double(cats), K = as.integer(K), target = as.integer(target),
#          cov_label = cov_label)
#   }

#   ac <- get_actor_cats(nw, cov, category, n1 = n1, ia = ia)
#   cats      <- ac$cats     # doubles côté C, valeurs entières
#   K         <- ac$K
#   target    <- ac$target   # 0 => non ciblé, sinon index 1..K
#   cov_label <- ac$cov_label

#   # --- 2) Normaliser le critère de tailles S (argument 'size') ----------------
#   S <- a$size
#   if (is.null(S) || length(S) == 0L) {
#     sizes <- integer(0)                      # S = toutes tailles
#   } else {
#     S <- unique(as.integer(S))
#     if (any(!is.finite(S)) || any(S < 1L))
#       ergm_Init_stop(sQuote(termname), ": 'size' doit contenir des entiers >= 1.")
#     sizes <- sort(S)
#   }

#   # --- 3) Inputs pour le C ----------------------------------------------------
#   # INPUT_PARAM layout:
#   #   [1]     = n1
#   #   [2]     = L = length(sizes)
#   #   [3..]   = sizes (L entrées, éventuellement L=0)
#   #   [next]  = K = nb de modalités (hors 0)
#   #   [next]  = target (0 => non ciblé, >0 => index de la modalité cible)
#   #   [ .. ]  = cats_1..cats_n1  (indices 0..K)
#   L <- length(sizes)
#   inputs <- c(
#     as.double(n1),
#     as.double(L),
#     as.double(sizes),
#     as.double(K),
#     as.double(target),
#     as.double(cats)
#   )

#   # --- 4) Nom du coefficient ---------------------------------------------------
#   pretty_sizes <- if (L == 0L) "all" else paste0("S{", paste(sizes, collapse = ","), "}")
#   coef.names <- paste0("cov_fullmatch[", cov_label, "]_", pretty_sizes)

#   # --- 5) Spécification du terme ----------------------------------------------
#   list(
#     name         = "cov_fullmatch",  # doit matcher C_CHANGESTAT_FN(c_cov_fullmatch)
#     coef.names   = coef.names,
#     inputs       = inputs,           # n1, L, sizes[L], K, target, cats[n1]
#     dependence   = TRUE,
#     emptynwstats = 0
#   )
# }