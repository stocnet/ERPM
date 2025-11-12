# ==============================================================================
# Fichier : R/InitErgmTerm.cov_fullmatch.R
# Terme   : cov_fullmatch
# Règle   : NA interdits dans la covariée (fail-fast)
# INPUT_PARAM : c(n1, L, sizes[L], K, target, cats[n1])
# Debug    : options(ERPM.cov_fullmatch.debug = TRUE) pour activer les logs
# ==============================================================================

InitErgmTerm.cov_fullmatch <- function(nw, arglist, ...) {
  termname <- "cov_fullmatch"

  # -- debug helpers ------------------------------------------------------------
  dbg <- isTRUE(getOption("ERPM.cov_fullmatch.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[cov_fullmatch][DEBUG]", ..., "\n", sep = "")

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                              "size",             "category"),
    vartypes      = c("character,numeric,logical,vector", "numeric,integer",  "character,numeric,logical"),
    defaultvalues = list(NULL,                             NULL,               NULL),
    required      = c(TRUE,                                FALSE,              FALSE)
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
    if (is.null(cov_vec)) stop(termname, ": attribut inexistant: ", sQuote(cov_raw), ".")
    dbgcat("cov source = vertex attribute ", sQuote(cov_label))
  } else {
    cov_vec   <- cov_raw
    cov_label <- "cov"
    dbgcat("cov source = vector literal")
  }
  if (length(cov_vec) < n1) stop(termname, ": longueur de la covariée < n1.")
  cov_vec <- cov_vec[seq_len(n1)]
  dbgcat("cov length = ", length(cov_vec), " | head = ", paste(utils::head(as.character(cov_vec), 6L), collapse = ","))

  # ----- FAIL-FAST NA -----
  if (anyNA(cov_vec)) stop(termname, ": NA non autorisé dans la covariée du mode acteurs.")

  # ----- normalisation en catégories entières 1..K -----
  if (is.logical(cov_vec)) cov_vec <- as.integer(cov_vec)
  f        <- factor(cov_vec)           # sans NA à ce stade
  levels_f <- levels(f)
  K        <- length(levels_f)
  if (K == 0L) stop(termname, ": aucune modalité valide trouvée.")
  cats     <- as.integer(f)             # 1..K
  dbgcat("K = ", K, " | levels = {", paste(levels_f, collapse = ","), "}")

  # ----- category ciblée -> target (0 si non fournie) -----
  has_category <- {
    x <- a$category
    !is.null(x) && length(x) == 1L && !is.na(x) && nzchar(as.character(x))
  }
  dbgcat("has_category = ", has_category)

  target <- 0L
  if (has_category) {
    cat_val <- a$category
    if (is.logical(cat_val)) cat_val <- as.integer(cat_val)
    ix <- match(as.character(cat_val), levels_f, nomatch = 0L)
    if (ix == 0L) {
      warning(termname, ": 'category' non trouvée dans les modalités; cible ignorée.", call. = FALSE)
      dbgcat("category ", sQuote(as.character(a$category)), " not found -> target=0 (ignored)")
    } else {
      target <- as.integer(ix)
      cov_label <- paste0(cov_label, "==", levels_f[target])
      dbgcat("category target = ", target, " -> label suffix = ", levels_f[target])
    }
  }

  # ----- filtre de tailles S -----
  sizes <- a$size
  if (is.null(sizes)) {
    L <- 0L
    sizes_vec <- numeric(0)
    dbgcat("size filter = <ALL> (L=0)")
  } else {
    if (!is.numeric(sizes))
      stop("cov_fullmatch: 'size' doit être numérique.")
    if (length(sizes) == 0L)
      stop("cov_fullmatch: 'size' vide (integer(0)) interdit. Utilisez NULL pour toutes tailles.")
    sizes <- as.integer(round(sizes))
    if (any(!is.finite(sizes)) || any(sizes <= 0L))
      stop("cov_fullmatch: 'size' doit contenir des entiers positifs.")
    sizes <- sort(unique(sizes))
    L <- length(sizes)
    sizes_vec <- as.double(sizes)
    dbgcat("size filter = {", paste(sizes, collapse = ","), "} (L=", L, ")")
  }

  # ----- nom de coefficient (rétrocompatible) -----
  size_label <- if (L > 0L) sprintf("_S{%s}", paste(sizes, collapse = ",")) else "_all"
  coef.name  <- sprintf("cov_fullmatch[%s]%s", cov_label, size_label)
  dbgcat("coef.name = ", coef.name)

  # ----- construire INPUT_PARAM -----
  inputs <- c(
    as.double(n1),
    as.double(L),
    sizes_vec,
    as.double(K),
    as.double(target),
    as.double(cats)
  )
  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L, " K=", K, " target=", target, " | cats[1:6]=",
         paste(utils::head(as.integer(cats), 6L), collapse = ","))

  list(
    name         = "cov_fullmatch",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, L, sizes[L], K, target, cats[n1]
    dependence   = TRUE,
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}