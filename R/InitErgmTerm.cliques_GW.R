# ==============================================================================
# File    : InitErgmTerm.cliques_GW.R
# Purpose : Declare the ERGM term 'cliques_GW' for bipartite networks
#           (geometrically weighted group sizes).
# Project : ERPM / ERGM extensions
# ============================================================================

#' ERGM term: cliques_GW (geometrically weighted group sizes)
#' @name InitErgmTerm.cliques_GW
#' @aliases cliques_GW
#' @note InitErgmTerm.cliques_GW.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' \code{cliques_GW} is an ERGM term for bipartite actor–group networks that
#' aggregates group sizes through a geometrically weighted transform. The
#' bipartite network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode}, whose size is given by \code{nw \%n\% "bipartite"};
#'   \item a \emph{group mode}, consisting of the remaining nodes that represent
#'         groups.
#' }
#'
#' For each group node \eqn{g} in the group mode, let \eqn{n_g} be its degree
#' (the number of adjacent actors). For a given \eqn{\lambda \ge 1}, define
#' \deqn{
#'   S(n_g, \lambda)
#'   =
#'   \lambda \Big[1 - r_\lambda^{\,n_g}\Big],
#'   \qquad
#'   r_\lambda = \frac{\lambda - 1}{\lambda}.
#' }
#' The \code{cliques_GW} term computes the statistic
#' \deqn{
#'   T_\lambda(y)
#'   =
#'   \sum_{g \in G} S(n_g, \lambda)
#'   =
#'   \sum_{g \in G} \lambda \Big[ 1 - r_\lambda^{\,n_g} \Big],
#' }
#' where \eqn{G} is the set of group-mode nodes. Intuitively, each group
#' contributes a geometrically weighted function of its size, with
#' \eqn{r_\lambda \in [0, 1)} whenever \eqn{\lambda > 1}.
#'
#' The initializer is vectorized in \code{lambda}: each value \eqn{\lambda_j}
#' produces one scalar statistic \eqn{T_{\lambda_j}(y)} and one corresponding
#' coefficient.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic.
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - This term MUST support multi-toggle proposals (swap/split/merge proposals
#'   decomposed into a list of membership toggles).
#' - Therefore the compiled change-statistic MUST be implemented using the
#'   D_CHANGESTAT_FN API (multi-toggle).
#' - On the R side, we MUST advertise this to ergm by returning `d_func = TRUE`.
#'   Otherwise ergm will call the one-toggle C_CHANGESTAT_FN entrypoint, causing
#'   a signature mismatch and typically a segfault.
#'
#' Compiled symbol naming convention:
#' - Recommended: implement the C function as `d_cliques_GW` via
#'   D_CHANGESTAT_FN(d_cliques_GW).
#' - Avoid exposing a symbol named `c_cliques_GW` with a D-signature.
#'
#' R initializer responsibilities:
#' \itemize{
#'   \item enforce bipartite network via \code{check.ErgmTerm()};
#'   \item accept one or several values of \code{lambda};
#'   \item validate domain: finite and \eqn{\lambda >= 1};
#'   \item precompute \eqn{r_\lambda = (\lambda - 1) / \lambda};
#'   \item pack \code{lambda} and \code{r_lambda} into \code{inputs} as interleaved pairs.
#' }
#'
#' INPUT_PARAM layout (C side):
#' \deqn{
#'   \text{INPUT\_PARAM}
#'   =
#'   (\lambda_1, r_{\lambda_1},
#'    \lambda_2, r_{\lambda_2},
#'    \dots,
#'    \lambda_J, r_{\lambda_J}).
#' }
#'
#' @note
#' For \eqn{\lambda = 1}, \eqn{r_\lambda = 0} implies:
#' - \eqn{S(0,1)=0}
#' - \eqn{S(n,1)=1} for any \eqn{n \ge 1}
#' so the statistic reduces to the number of non-empty groups.
#'
#' @keywords ERGM term bipartite groups cliques geometric
#' @md
#'
#' @template erpm-initergmterm-args
#'
#' @export
InitErgmTerm.cliques_GW <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cliques_GW"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.cliques_GW.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.cliques_GW.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[cliques_GW][DEBUG]", ..., "\n", sep = "")

  dbgcat("InitErgmTerm.cliques_GW called with args: ",
         paste(names(arglist), collapse = ", "))

  # ---------------------------------------------------------------------------
  # Base validation and structural requirements
  # ---------------------------------------------------------------------------
  # - Enforce bipartite network (actor mode / group mode).
  # - Declare allowed arguments: lambda.
  # - Provide default: lambda = 2.
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("lambda"),
    vartypes      = c("numeric"),
    defaultvalues = list(2),
    required      = c(FALSE)
  )

  lambda <- a$lambda
  if (length(lambda) == 0L) {
    dbgcat("lambda is empty -> returning NULL term")
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Domain validation for lambda
  # ---------------------------------------------------------------------------
  # - finite values only
  # - lambda >= 1
  if (any(!is.finite(lambda))) {
    ergm_Init_stop(sQuote(termname), ": 'lambda' must be finite (no NA/NaN/Inf).")
  }
  if (any(lambda < 1)) {
    ergm_Init_stop(sQuote(termname), ": 'lambda' must be >= 1.")
  }

  # ---------------------------------------------------------------------------
  # Precompute r = (lambda - 1) / lambda
  # ---------------------------------------------------------------------------
  r <- (lambda - 1) / lambda
  dbgcat("lambda = ", paste(lambda, collapse = ", "),
         " | r = ", paste(signif(r, 8), collapse = ", "))

  # ---------------------------------------------------------------------------
  # Coefficient names (vectorized)
  # ---------------------------------------------------------------------------
  pretty_num <- function(x) {
    s <- formatC(x, digits = 6, format = "fg", flag = "#")
    sub("\\.$", "", s)
  }
  coef.names <- paste0("cliques_GW_lambda", vapply(lambda, pretty_num, ""))
  dbgcat("coef.names = ", paste(coef.names, collapse = ", "))

  # ---------------------------------------------------------------------------
  # Pack INPUT_PARAM as interleaved (lambda_j, r_j) pairs
  # ---------------------------------------------------------------------------
  # Layout:
  #   [2*j + 0] = lambda_j
  #   [2*j + 1] = r_j
  inputs <- c(rbind(as.double(lambda), as.double(r)))
  dbgcat("inputs length = ", length(inputs),
         " (", length(lambda), " lambdas)")

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  # IMPORTANT:
  # - `d_func = TRUE` tells ergm to call the multi-toggle (D_) changestat entrypoint.
  # - Without it, ergm assumes a one-toggle C_ changestat and will call the function
  #   with the wrong signature if you compiled only a D_ function (=> segfault).
  list(
    name         = "cliques_GW",
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,
    d_func       = TRUE,                    # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    emptynwstats = numeric(length(lambda))
  )
}