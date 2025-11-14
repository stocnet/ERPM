# ==============================================================================
# Fichier : R/InitErgmTerm.cov_fulldiff.R
# Terme   : cov_fulldiff
# Stat    : T = sum_g 1[n_g in S] * (x_g^max - x_g^min)
# Règle   : NA interdits dans la covariée (fail-fast)
# INPUT_PARAM : c(n1, L, sizes[L], x[1..n1])
# Debug    : options(ERPM.cov_fulldiff.debug = TRUE) pour activer les logs
# ==============================================================================

InitErgmTerm.cov_fulldiff <- function(nw, arglist, ...) {
  termname <- "cov_fulldiff"

  # -- debug helpers ------------------------------------------------------------
  dbg    <- isTRUE(getOption("ERPM.cov_fulldiff.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[cov_fulldiff][DEBUG]", ..., "\n", sep = "")

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                              "size"),
    vartypes      = c("character,numeric,logical,vector", "numeric,integer"),
    defaultvalues = list(NULL,                             NULL),
    required      = c(TRUE,                                FALSE)
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
    if (is.null(cov_vec))
      stop(termname, ": attribut inexistant: ", sQuote(cov_raw), ".")
    dbgcat("cov source = vertex attribute ", sQuote(cov_label))
  } else {
    cov_vec   <- cov_raw
    cov_label <- "cov"
    dbgcat("cov source = vector literal")
  }

  if (length(cov_vec) < n1)
    stop(termname, ": longueur de la covariée < n1.")

  cov_vec <- cov_vec[seq_len(n1)]
  dbgcat("cov length = ", length(cov_vec),
         " | head = ", paste(utils::head(as.character(cov_vec), 6L), collapse = ","))

  # ----- coercition numérique + FAIL-FAST NA -----
  if (is.logical(cov_vec) || is.integer(cov_vec)) {
    cov_vec <- as.numeric(cov_vec)
  }
  if (!is.numeric(cov_vec))
    stop(termname, ": la covariée doit être coercible en numérique.")

  if (anyNA(cov_vec))
    stop(termname, ": NA non autorisé dans la covariée du mode acteurs.")

  # ----- filtre de tailles S -----
  sizes <- a$size
  if (is.null(sizes)) {
    L         <- 0L
    sizes_vec <- numeric(0)
    dbgcat("size filter = <ALL> (L=0)")
    size_label <- "_all"
  } else {
    if (!is.numeric(sizes))
      stop(termname, ": 'size' doit être numérique.")
    if (length(sizes) == 0L)
      stop(termname, ": 'size' vide (integer(0)) interdit. Utilisez NULL pour toutes tailles.")

    sizes <- as.integer(round(sizes))
    if (any(!is.finite(sizes)) || any(sizes <= 0L))
      stop(termname, ": 'size' doit contenir des entiers positifs.")

    sizes     <- sort(unique(sizes))
    L         <- length(sizes)
    sizes_vec <- as.double(sizes)
    dbgcat("size filter = {", paste(sizes, collapse = ","), "} (L=", L, ")")
    size_label <- sprintf("_S{%s}", paste(sizes, collapse = ","))
  }

  # ----- nom de coefficient -----
  coef.name <- sprintf("cov_fulldiff[%s]%s", cov_label, size_label)
  dbgcat("coef.name = ", coef.name)

  # ----- construire INPUT_PARAM -----
  inputs <- c(
    as.double(n1),
    as.double(L),
    sizes_vec,
    as.double(cov_vec)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L,
         " | cov[1:6]=", paste(utils::head(signif(cov_vec, 5L), 6L), collapse = ","))

  list(
    name         = "cov_fulldiff",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, L, sizes[L], x[n1]
    dependence   = TRUE,
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}
