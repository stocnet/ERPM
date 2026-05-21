# ==============================================================================
# File    : R/InitErgmTerm.cov_fullmatch.R
# Purpose : Declare the ERGM term 'cov_fullmatch' for bipartite networks
#           (group-level unanimity on an actor covariate).
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM term: cov_fullmatch (group-level unanimity)
#' @name InitErgmTerm.cov_fullmatch
#' @aliases cov_fullmatch
#' @note InitErgmTerm.cov_fullmatch.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' \code{cov_fullmatch} is an ERGM term for bipartite networks that counts
#' groups whose actors are unanimously homogeneous on a categorical covariate.
#' The bipartite network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side that represents groups).
#' }
#'
#' Each actor in the actor mode carries a covariate value, which is mapped
#' internally to integer categories \eqn{1,\dots,K}. For each group node, we
#' consider the set of adjacent actors (its group of members) and check whether
#' all actors in that group share the same category. An optional size filter
#' restricts the set of group sizes that contribute to the statistic, and an
#' optional \code{category} argument targets a single covariate level.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic.
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - This term MUST support multi-toggle moves (swap/split/merge) represented
#'   as a list of edge toggles in ERGM's MCMC.
#' - Therefore, the compiled change-statistic must be implemented using the
#'   D_CHANGESTAT_FN API (multi-toggle).
#' - On the R side, we MUST advertise this to ergm by returning `d_func = TRUE`.
#'   Otherwise ergm will call the one-toggle entrypoint (C_CHANGESTAT_FN) with
#'   the wrong signature (=> crash / segfault).
#'
#' Compiled symbol naming convention:
#' - Recommended: implement the C function as `d_cov_fullmatch` via
#'   D_CHANGESTAT_FN(d_cov_fullmatch).
#' - Keeping a legacy one-toggle `c_cov_fullmatch` is fine for backward
#'   compatibility, but the ERGM term should be configured to call the D_ entrypoint.
#'
#' The R initializer below:
#' \itemize{
#'   \item enforces that the network is bipartite via \code{nw \%n\% "bipartite"};
#'   \item extracts and encodes an actor-mode covariate as integer categories;
#'   \item normalizes the optional \code{size} filter into a sorted set of
#'         positive integers;
#'   \item encodes an optional targeted category into an integer code
#'         \eqn{\texttt{target} \in \{0,\dots,KEst-ce\}} (with \code{0} = "no target");
#'   \item packs a compact \code{inputs} vector for the C layer.
#' }
#'
#' @template erpm-initergmterm-args
#'
#' @export
InitErgmTerm.cov_fullmatch <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_fullmatch"

  # Debug helpers (controlled by an R option):
  #   options(ERPM.cov_fullmatch.debug = TRUE/FALSE)
  dbg <- isTRUE(getOption("ERPM.cov_fullmatch.debug", TRUE))
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

  # ----- 1) Actor-mode size n1 -----------------------------------------------
  # n1 is the number of actors, as stored in the bipartite network attribute.
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) ergm_Init_stop(sQuote(termname), ": strictly bipartite network required.")
  dbgcat("n1 = ", n1)

  # ----- 2) Extract the actor covariate vector (length >= n1) -----------------
  # The covariate can be given as:
  # - the name of a vertex attribute (preferred);
  # - a literal vector. In all cases we keep only the first n1 entries for actors.
  cov_raw <- a$cov
  if (is.character(cov_raw) && length(cov_raw) == 1L) {
    cov_vec   <- network::get.vertex.attribute(nw, cov_raw)
    cov_label <- cov_raw
    if (is.null(cov_vec)) ergm_Init_stop(sQuote(termname), ": attribut inexistant: ", sQuote(cov_raw), ".")
    dbgcat("cov source = vertex attribute ", sQuote(cov_label))
  } else {
    cov_vec   <- cov_raw
    cov_label <- "cov"
    dbgcat("cov source = vector literal")
  }
  if (length(cov_vec) < n1) ergm_Init_stop(sQuote(termname), ": covariate length < n1.")
  cov_vec <- cov_vec[seq_len(n1)]
  dbgcat("cov length = ", length(cov_vec), " | head = ", paste(utils::head(as.character(cov_vec), 6L), collapse = ","))

  # ----- 3) Fail-fast on missing values ---------------------------------------
  # Any NA on the actor mode is rejected, since the change-statistic expects
  # well-defined categories for all actors.
  if (anyNA(cov_vec)) ergm_Init_stop(sQuote(termname), ": NA values are not allowed in the actor-mode covariate.")

  # ----- 4) Normalize to integer category codes 1..K --------------------------
  if (is.logical(cov_vec)) cov_vec <- as.integer(cov_vec)
  f        <- factor(cov_vec)     # no NA at this stage
  levels_f <- levels(f)
  K        <- length(levels_f)
  if (K == 0L) ergm_Init_stop(sQuote(termname), ": no valid modality found.")
  cats     <- as.integer(f)       # 1..K
  dbgcat("K = ", K, " | levels = {", paste(levels_f, collapse = ","), "}")

  # ----- 5) Targeted category -> 'target' code (0 if not provided) ------------
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
      warning(termname, ": 'category' not found among the modalities; target ignored.", call. = FALSE)
      dbgcat("category ", sQuote(as.character(a$category)), " not found -> target=0 (ignored)")
    } else {
      target <- as.integer(ix)
      cov_label <- paste0(cov_label, "==", levels_f[target])
      dbgcat("category target = ", target, " -> label suffix = ", levels_f[target])
    }
  }

  # ----- 6) Size filter S (argument 'size') -----------------------------------
  sizes <- a$size
  if (is.null(sizes)) {
    L <- 0L
    sizes_vec <- numeric(0)
    dbgcat("size filter = <ALL> (L=0)")
  } else {
    if (!is.numeric(sizes))
      ergm_Init_stop(sQuote(termname), ": 'size' must be numeric.")
    if (length(sizes) == 0L)
      ergm_Init_stop(sQuote(termname), ": empty 'size' (integer(0)) is not allowed. Use NULL for all sizes.")
    sizes <- as.integer(round(sizes))
    if (any(!is.finite(sizes)) || any(sizes <= 0L))
      ergm_Init_stop(sQuote(termname), ": 'size' must contain positive integers.")
    sizes <- sort(unique(sizes))
    L <- length(sizes)
    sizes_vec <- as.double(sizes)
    dbgcat("size filter = {", paste(sizes, collapse = ","), "} (L=", L, ")")
  }

  # ----- 7) Coefficient name (kept backward compatible) -----------------------
  size_label <- if (L > 0L) sprintf("_S{%s}", paste(sizes, collapse = ",")) else "_all"
  coef.name  <- sprintf("cov_fullmatch[%s]%s", cov_label, size_label)
  dbgcat("coef.name = ", coef.name)

  # ----- 8) Build INPUT_PARAM for the C layer ---------------------------------
  # Layout (double vector):
  #   n1, L, sizes[L], K, target, cats[n1]
  inputs <- c(
    as.double(n1),
    as.double(L),
    sizes_vec,
    as.double(K),
    as.double(target),
    as.double(cats)
  )
  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L, " K=", K, " target=", target,
         " | cats[1:6]=", paste(utils::head(as.integer(cats), 6L), collapse = ","))

  # ----- 9) Return the ERGM term specification --------------------------------
  # IMPORTANT:
  # - `d_func = TRUE` is REQUIRED because the compiled changestat is D_ (multi-toggle).
  list(
    name         = "cov_fullmatch",
    coef.names   = coef.name,
    inputs       = inputs,
    dependence   = TRUE,
    d_func       = TRUE,   # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}