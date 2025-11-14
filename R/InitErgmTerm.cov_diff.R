# ==============================================================================
# Fichier : R/InitErgmTerm.cov_diff.R
# Terme   : cov_diff
# Stat    : T_k = sum_g sum_{S ⊂ g, |S|=k} (max_{i∈S} x_i - min_{i∈S} x_i)
# Option  : normalized = TRUE  => moyenne par groupe (division par C(n_g, k))
# Règle   : NA interdits dans la covariée (fail-fast)
# INPUT_PARAM : c(n1, k, norm_flag, x[1..n1])
# Debug   : options(ERPM.cov_diff.debug = TRUE) pour activer les logs
# ==============================================================================

InitErgmTerm.cov_diff <- function(nw, arglist, ...) {
  termname <- "cov_diff"

  # -- debug helpers ------------------------------------------------------------
  dbg    <- isTRUE(getOption("ERPM.cov_diff.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[cov_diff][DEBUG]", ..., "\n", sep = "")

  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                              "clique_size",            "normalized"),
    vartypes      = c("character,numeric,logical,vector", "numeric,integer",        "logical,character,numeric"),
    defaultvalues = list(NULL,                             2,                         FALSE),
    required      = c(TRUE,                                FALSE,                     FALSE)
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

  # ----- contrainte numérique + FAIL-FAST NA -----
  if (is.logical(cov_vec) || is.integer(cov_vec)) {
    cov_vec <- as.numeric(cov_vec)
  }
  if (!is.numeric(cov_vec))
    stop(termname, ": la covariée doit être convertissable en numérique.")

  if (anyNA(cov_vec))
    stop(termname, ": NA non autorisé dans la covariée du mode acteurs.")

  # ----- clique_size : entier k >= 2 -----
  k_raw <- a$clique_size
  if (length(k_raw) != 1L || !is.numeric(k_raw))
    stop(termname, ": 'clique_size' doit être un scalaire numérique.")
  k <- as.integer(round(k_raw))
  if (!is.finite(k) || k < 2L)
    stop(termname, ": 'clique_size' doit être un entier >= 2.")
  dbgcat("clique_size k = ", k)

  # ----- normalized : FALSE / TRUE (none / by_group) ---------------------------
  norm_raw <- a$normalized
  norm_flag <- 0L
  norm_label <- "raw"

  if (is.null(norm_raw)) {
    norm_flag  <- 0L
    norm_label <- "raw"
  } else if (is.logical(norm_raw)) {
    norm_flag  <- if (isTRUE(norm_raw)) 1L else 0L
    norm_label <- if (norm_flag == 1L) "by_group" else "raw"
  } else if (is.character(norm_raw) && length(norm_raw) == 1L) {
    nr <- match.arg(norm_raw, c("none", "by_group"))
    norm_flag  <- if (nr == "by_group") 1L else 0L
    norm_label <- nr
  } else if (is.numeric(norm_raw) && length(norm_raw) == 1L) {
    norm_flag  <- if (norm_raw != 0) 1L else 0L
    norm_label <- if (norm_flag == 1L) "by_group" else "raw"
  } else {
    stop(termname, ": 'normalized' doit être logique, numérique ou 'none'/'by_group'.")
  }

  dbgcat("normalized = ", norm_label, " (flag=", norm_flag, ")")

  # ----- nom de coefficient ----------------------------------------------------
  coef.name <- sprintf(
    "cov_diff[%s]_k%d%s",
    cov_label, k,
    if (norm_flag == 1L) "_norm" else ""
  )
  dbgcat("coef.name = ", coef.name)

  # ----- construire INPUT_PARAM ------------------------------------------------
  inputs <- c(
    as.double(n1),
    as.double(k),
    as.double(norm_flag),
    as.double(cov_vec)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " k=", k, " norm_flag=", norm_flag,
         " | cov[1:6]=", paste(utils::head(signif(cov_vec, 5L), 6L), collapse = ","))

  list(
    name         = "cov_diff",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, k, norm_flag, x[n1]
    dependence   = TRUE,
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}
