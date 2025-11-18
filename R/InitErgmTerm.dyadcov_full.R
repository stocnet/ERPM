# ==============================================================================
# Fichier : R/InitErgmTerm.dyadcov_full.R
# Terme   : dyadcov_full
# Stat    : T = sum_g 1[n_g in S] * sum_{i<j, i,j in g} z_{ij}
# Règle   : NA interdits dans la matrice dyadique (fail-fast)
# INPUT_PARAM : c(n1, L, sizes[L], Z[n1*n1])  (Z en ordre colonne-major)
# Debug    : options(ERPM.dyadcov_full.debug = TRUE) pour activer les logs
# ==============================================================================

InitErgmTerm.dyadcov_full <- function(nw, arglist, ...) {
  termname <- "dyadcov_full"

  # -- debug helpers ------------------------------------------------------------
  dbg    <- isTRUE(getOption("ERPM.dyadcov_full.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[dyadcov_full][DEBUG]", ..., "\n", sep = "")

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("dyadcov",                         "size"),
    vartypes      = c("matrix,character",               "numeric,integer"),
    defaultvalues = list(NULL,                           NULL),
    required      = c(TRUE,                              FALSE)
  )

  # ----- n1 : taille du mode acteurs -----
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L)
    stop(termname, ": réseau non biparti strict (attribut %n% 'bipartite' manquant ou invalide).")
  dbgcat("n1 = ", n1)

  # ----- récupérer la matrice dyadique (n1 x n1) -----
  dyad_raw <- a$dyadcov
  dyad_label <- NULL

  if (is.character(dyad_raw) && length(dyad_raw) == 1L) {
    # cas : nom d'attribut %n% portant la matrice
    dyad_mat <- nw %n% dyad_raw
    dyad_label <- dyad_raw
    if (is.null(dyad_mat))
      stop(termname, ": attribut de niveau réseau inexistant: ", sQuote(dyad_raw), ".")
    dbgcat("dyadcov source = network attribute ", sQuote(dyad_label))
  } else {
    # cas : matrice passée littéralement
    dyad_mat <- dyad_raw
    dyad_label <- "dyadcov"
    dbgcat("dyadcov source = matrix literal")
  }

  if (!is.matrix(dyad_mat))
    stop(termname, ": 'dyadcov' doit être une matrice ou le nom d'un attribut de niveau réseau.")

  nr <- nrow(dyad_mat)
  nc <- ncol(dyad_mat)

  if (nr < n1 || nc < n1)
    stop(termname, ": dimensions de la matrice dyadique (", nr, "x", nc,
         ") insuffisantes pour n1 = ", n1, ".")

  # On restreint, si besoin, au bloc supérieur gauche n1 x n1
  if (nr > n1 || nc > n1) {
    dyad_mat <- dyad_mat[seq_len(n1), seq_len(n1), drop = FALSE]
    dbgcat("dyadcov truncated to ", n1, "x", n1)
  }

  # ----- coercition numérique + FAIL-FAST NA -----
  if (!is.numeric(dyad_mat))
    stop(termname, ": la matrice dyadique doit être numérique.")

  if (anyNA(dyad_mat))
    stop(termname, ": NA non autorisé dans la matrice dyadique.")

  # ----- contrôle optionnel de symétrie -----
  # petite tolérance numérique
  tol <- 1e-8
  if (max(abs(dyad_mat - t(dyad_mat))) > tol) {
    stop(termname, ": la matrice dyadique doit être symétrique (différence > tolérance).")
  }

  dbgcat("dyadcov dim = ", paste(dim(dyad_mat), collapse = "x"),
         " | sample = ",
         paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L), collapse = ","))

  # ----- filtre de tailles S -----
  sizes <- a$size
  if (is.null(sizes)) {
    L         <- 0L
    sizes_vec <- numeric(0)
    dbgcat("size filter = <ALL> (L=0)")
    size_label <- "_all"
  } else {
    if (!is.numeric(sizes))
      stop(termname, ": 'size' doit être numérique ou entier.")
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
  coef.name <- sprintf("dyadcov_full[%s]%s", dyad_label, size_label)
  dbgcat("coef.name = ", coef.name)

  # ----- construire INPUT_PARAM -----
  # as.double(matrice) => ordre colonne-major compatible avec l'indexation C
  inputs <- c(
    as.double(n1),
    as.double(L),
    sizes_vec,
    as.double(dyad_mat)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L,
         " | Z[1:6]=", paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L),
                             collapse = ","))

  list(
    name         = "dyadcov_full",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, L, sizes[L], Z[n1*n1]
    dependence   = TRUE,
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}