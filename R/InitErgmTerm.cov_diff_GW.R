# ==============================================================================
# File    : R/InitErgmTerm.cov_diff_GW.R
# Purpose : Declare the ERGM term 'cov_diff_GW' for bipartite networks
#           (multi-toggle / D_CHANGESTAT_FN).
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM term: cov_diff_GW (geometrically weighted range over k-actor subsets)
#'
#' @name InitErgmTerm.cov_diff_GW
#' @aliases cov_diff_GW
#' @note InitErgmTerm.cov_diff_GW.R
#'
#' @description
#' \code{cov_diff_GW} is an ERGM term for bipartite networks that builds a
#' geometrically weighted combination of \code{cov_diff}-type statistics over
#' all subset sizes \eqn{k \ge 2} inside each group. The bipartite network is
#' interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side, representing groups).
#' }
#'
#' Each actor in the actor mode carries a numeric covariate value \eqn{x_i}. For
#' each group and each subset size \eqn{k \ge 2}, \code{cov_diff_GW} considers
#' all \eqn{k}-actor subsets inside that group and computes on each subset the
#' max-min range of the covariate. These \code{cov_diff}-type contributions are
#' then combined over all \eqn{k \ge 2} using a geometric weight controlled by
#' \eqn{\lambda > 1}.
#'
#' For a given value of \eqn{\lambda > 1}, let
#' \eqn{c_k} denote the \code{cov_diff}-style statistic of order \eqn{k}, i.e.
#' the sum of ranges over all \eqn{k}-actor subsets in all groups. The term
#' \code{cov_diff_GW} is defined as
#' \deqn{
#'   T_{\mathrm{GW}}(\lambda)
#'   =
#'   \sum_{k \ge 2} \left(-\frac{1}{\lambda}\right)^{k-1} c_k.
#' }
#' The initializer supports a vector of \eqn{\lambda} values and returns one
#' scalar statistic per value.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic.
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - This term MUST support multi-toggle proposals (swap/split/merge decomposed
#'   into multiple toggles).
#' - Therefore, the compiled change-statistic MUST be implemented using the
#'   D_CHANGESTAT_FN API (multi-toggle).
#' - On the R side, we MUST advertise this to ergm by returning `d_func = TRUE`.
#'   Otherwise ergm will call the changestat as a one-toggle C_CHANGESTAT_FN,
#'   causing a signature mismatch (often a segfault).
#'
#' Compiled symbol naming convention:
#' - Implement the C function as `d_cov_diff_GW` via D_CHANGESTAT_FN(d_cov_diff_GW).
#' - Avoid exposing a symbol named `c_cov_diff_GW` with a D-signature.
#'
#' Debugging output for the initializer can be enabled via:
#' \preformatted{
#'   options(ERPM.cov_diff_GW.debug = TRUE)
#' }
#'
#' @param nw A \pkg{network} object.
#' @param arglist A named list of term arguments constructed by \pkg{ergm}.
#'   Expected components include \code{cov} (vertex attribute name or numeric
#'   vector) and \code{lambda} (numeric scalar or vector with values > 1).
#' @param ... Passed through by \pkg{ergm}; not used.
#' @param version ERGM API version; not used.
#'
#' @export
InitErgmTerm.cov_diff_GW <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cov_diff_GW"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.cov_diff_GW.debug = TRUE/FALSE)
  dbg    <- isTRUE(getOption("ERPM.cov_diff_GW.debug", TRUE))
  dbgcat <- function(...) if (dbg) cat("[cov_diff_GW][DEBUG]", ..., "\n", sep = "")

  # ---------------------------------------------------------------------------
  # Parse and validate term arguments via ergm helpers
  # ---------------------------------------------------------------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("cov",                              "lambda"),
    vartypes      = c("character,numeric,logical,vector", "numeric,vector"),
    defaultvalues = list(NULL,                             2),
    required      = c(TRUE,                                FALSE)
  )

  # ---------------------------------------------------------------------------
  # 1) Actor-mode size n1
  # ---------------------------------------------------------------------------
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) {
    ergm_Init_stop(sQuote(termname), ": strictly bipartite network required (valid 'bipartite' attribute).")
  }
  dbgcat("n1 = ", n1)

  # ---------------------------------------------------------------------------
  # 2) Extract actor covariate (cov)
  # ---------------------------------------------------------------------------
  # The covariate may be:
  # - a vertex attribute name (character scalar),
  # - a literal vector (numeric/integer/logical).
  cov_raw <- a$cov
  if (is.character(cov_raw) && length(cov_raw) == 1L) {
    cov_vec   <- network::get.vertex.attribute(nw, cov_raw)
    cov_label <- cov_raw
    if (is.null(cov_vec)) {
      ergm_Init_stop(sQuote(termname), ": attribut inexistant: ", sQuote(cov_raw), ".")
    }
    dbgcat("cov source = vertex attribute ", sQuote(cov_label))
  } else {
    cov_vec   <- cov_raw
    cov_label <- "cov"
    dbgcat("cov source = vector literal")
  }

  if (length(cov_vec) < n1) {
    ergm_Init_stop(sQuote(termname), ": covariate length < n1.")
  }

  # Keep actor mode only (first n1 entries).
  cov_vec <- cov_vec[seq_len(n1)]
  dbgcat("cov length = ", length(cov_vec),
         " | head = ", paste(utils::head(as.character(cov_vec), 6L), collapse = ","))

  # Numeric coercion + fail-fast on NA.
  if (is.logical(cov_vec) || is.integer(cov_vec)) cov_vec <- as.numeric(cov_vec)
  if (!is.numeric(cov_vec)) {
    ergm_Init_stop(sQuote(termname), ": the covariate must be coercible to numeric.")
  }
  cov_vec <- as.double(cov_vec)
  if (anyNA(cov_vec)) {
    ergm_Init_stop(sQuote(termname), ": NA values are not allowed in the actor-mode covariate.")
  }

  # ---------------------------------------------------------------------------
  # 3) Lambda (vectorized), must be finite and strictly > 1
  # ---------------------------------------------------------------------------
  lambda_raw <- a$lambda
  if (!is.numeric(lambda_raw) || length(lambda_raw) < 1L) {
    ergm_Init_stop(sQuote(termname), ": 'lambda' must be numeric (scalar or vector).")
  }
  lambda_vec <- as.double(lambda_raw)
  if (any(!is.finite(lambda_vec))) {
    ergm_Init_stop(sQuote(termname), ": 'lambda' must be finite.")
  }
  if (any(lambda_vec <= 1)) {
    ergm_Init_stop(sQuote(termname), ": all values of 'lambda' must be > 1.")
  }
  L <- length(lambda_vec)
  dbgcat("lambda_vec = {", paste(signif(lambda_vec, 6L), collapse = ", "), "} (L = ", L, ")")

  # ---------------------------------------------------------------------------
  # 4) Coefficient names
  # ---------------------------------------------------------------------------
  coef.names <- sprintf("cov_diff_GW[%s]_lambda%.6g", cov_label, lambda_vec)
  dbgcat("coef.names = ", paste(coef.names, collapse = " | "))

  # ---------------------------------------------------------------------------
  # 5) INPUT_PARAM for the C layer
  # ---------------------------------------------------------------------------
  # Layout:
  #   INPUT_PARAM = c(
  #     n1,
  #     L,
  #     lambda[1:L],
  #     x[1:n1]
  #   )
  inputs <- c(
    as.double(n1),
    as.double(L),
    as.double(lambda_vec),
    as.double(cov_vec)
  )

  dbgcat("inputs len=", length(inputs),
         " | n1=", n1, " L=", L,
         " | lambda head=", paste(utils::head(signif(lambda_vec, 6L), 6L), collapse = ","),
         " | cov head=",    paste(utils::head(signif(cov_vec,    6L), 6L), collapse = ","))

  # ---------------------------------------------------------------------------
  # 6) Return ERGM term specification
  # ---------------------------------------------------------------------------
  # IMPORTANT:
  # - `d_func = TRUE` is REQUIRED (multi-toggle / D_CHANGESTAT_FN).
  list(
    name         = "cov_diff_GW",
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,
    d_func       = TRUE,        # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = numeric(L)
  )
}