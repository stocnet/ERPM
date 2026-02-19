# ==============================================================================
# File    : R/InitErgmTerm.cov_fulldiff.R
# Purpose : Declare the ERGM term 'cov_fulldiff' for bipartite networks
#           (within-group covariate range), MULTI-TOGGLE SAFE.
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM term: cov_fulldiff (within-group range of a numeric covariate)
#' @name InitErgmTerm.cov_fulldiff
#' @aliases cov_fulldiff
#' @note InitErgmTerm.cov_fulldiff.R
#'
#' @description
#' \code{cov_fulldiff} is an ERGM term for bipartite networks that measures,
#' for each group, the dispersion of a numeric actor covariate through the
#' max-min range, optionally restricted to a subset of group sizes.
#'
#' The bipartite network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side, representing groups).
#' }
#'
#' Each actor in the actor mode carries a numeric covariate value \eqn{x_i}. For
#' each group node, we look at the actors adjacent to that group and compute
#' the within-group range \eqn{x_g^{\max} - x_g^{\min}}. An optional \code{size}
#' filter restricts which group sizes contribute to the statistic.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic.
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - This term is intended to support multi-toggle proposals (swap/split/merge
#'   represented as a list of toggles).
#' - Therefore, the compiled change-statistic MUST be implemented using the
#'   D_CHANGESTAT_FN API (multi-toggle), e.g. `d_cov_fulldiff`.
#' - On the R side, we MUST advertise this to ergm by returning `d_func = TRUE`.
#'   Otherwise ergm will call the one-toggle entrypoint (C_CHANGESTAT_FN) with the
#'   wrong signature, which typically crashes (segfault).
#'
#' Compiled symbol naming convention:
#' - Recommended: implement the C function as `d_cov_fulldiff` via
#'   D_CHANGESTAT_FN(d_cov_fulldiff).
#' - Avoid exposing a symbol named `c_cov_fulldiff` with a D-signature.
#'
#' The R initializer below:
#' \itemize{
#'   \item enforces that the network is bipartite via \code{nw \%n\% "bipartite"};
#'   \item extracts an actor-level covariate from a vertex attribute or from a
#'         literal vector, coercing it to numeric and failing fast on \code{NA};
#'   \item normalizes the optional \code{size} argument into a sorted set of
#'         strictly positive integers;
#'   \item packs the actor covariate and the size filter into a compact
#'         \code{INPUT_PARAM} vector for the C layer.
#' }
#'
#' @export
InitErgmTerm.cov_fulldiff <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_fulldiff"

  # Debug helpers:
  # - dbg: logical flag, controlled by an R option;
  # - dbgcat(): emits prefixed debug messages when dbg is TRUE.
  #
  # NOTE: This is R-side debug only. C-side debug remains controlled by compile-time
  # macros in changestat_cov_fulldiff.c.
  dbg    <- isTRUE(getOption("ERPM.cov_fulldiff.debug", TRUE))
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

  # ----- 1) Actor-mode size n1 -----------------------------------------------
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) stop(termname, ": strictly bipartite network required.")
  dbgcat("n1 = ", n1)

  # ----- 2) Extract actor covariate (length >= n1) ---------------------------
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
    stop(termname, ": covariate length < n1.")

  cov_vec <- cov_vec[seq_len(n1)]
  dbgcat("cov length = ", length(cov_vec),
         " | head = ", paste(utils::head(as.character(cov_vec), 6L), collapse = ","))

  # ----- 3) Numeric coercion + fail-fast on NA -------------------------------
  if (is.logical(cov_vec) || is.integer(cov_vec)) {
    cov_vec <- as.numeric(cov_vec)
  }
  if (!is.numeric(cov_vec))
    stop(termname, ": the covariate must be coercible to numeric.")

  if (anyNA(cov_vec))
    stop(termname, ": NA values are not allowed in the actor-mode covariate.")

  # ----- 4) Size filter S (argument 'size') ----------------------------------
  sizes <- a$size
  if (is.null(sizes)) {
    L         <- 0L
    sizes_vec <- numeric(0)
    dbgcat("size filter = <ALL> (L=0)")
    size_label <- "_all"
  } else {
    if (!is.numeric(sizes))
      stop(termname, ": 'size' must be numeric.")
    if (length(sizes) == 0L)
      stop(termname, ": empty 'size' (integer(0)) is not allowed. Use NULL for all sizes.")

    sizes <- as.integer(round(sizes))
    if (any(!is.finite(sizes)) || any(sizes <= 0L))
      stop(termname, ": 'size' must contain positive integers.")

    sizes     <- sort(unique(sizes))
    L         <- length(sizes)
    sizes_vec <- as.double(sizes)
    dbgcat("size filter = {", paste(sizes, collapse = ","), "} (L=", L, ")")
    size_label <- sprintf("_S{%s}", paste(sizes, collapse = ","))
  }

  # ----- 5) Coefficient name --------------------------------------------------
  coef.name <- sprintf("cov_fulldiff[%s]%s", cov_label, size_label)
  dbgcat("coef.name = ", coef.name)

  # ----- 6) Build INPUT_PARAM for the C layer --------------------------------
  # Layout:
  #   c(n1, L, sizes[1:L], x[1:n1])
  inputs <- c(
    as.double(n1),
    as.double(L),
    sizes_vec,
    as.double(cov_vec)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L,
         " | cov[1:6]=", paste(utils::head(signif(cov_vec, 5L), 6L), collapse = ","))

  # ----- 7) Return ERGM term specification -----------------------------------
  # IMPORTANT:
  # - d_func = TRUE tells ergm to call the D_ change-statistic entrypoint
  #   (multi-toggle signature).
  list(
    name         = "cov_fulldiff",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, L, sizes[L], x[n1]
    dependence   = TRUE,
    d_func       = TRUE,        # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    minval       = 0,
    maxval       = Inf,
    emptynwstats = 0
  )
}