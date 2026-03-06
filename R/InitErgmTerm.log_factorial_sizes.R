# ==============================================================================
# File    : R/InitErgmTerm.log_factorial_sizes.R
# Purpose : Declare the ERGM term 'log_factorial_sizes' for bipartite networks
#           (aggregation over the group mode) — MULTI-TOGGLE READY (D_ changestat).
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM term: log_factorial_sizes (group-mode degrees)
#'
#' @name InitErgmTerm.log_factorial_sizes
#' @aliases log_factorial_sizes
#' @note InitErgmTerm.log_factorial_sizes.R
#'
#' @description
#' \code{log_factorial_sizes} is an ERGM term for bipartite networks that sums
#' the log-factorial of group sizes. The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side used as \code{\%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each node in the group mode, we compute the degree \eqn{\deg(g)} (number
#' of adjacent actors) and accumulate
#' \deqn{
#'   \sum_{g \in \text{group mode}} \log\Gamma(\deg(g)),
#' }
#' with the convention \code{lgamma(0) = 0}, so that empty groups do not
#' contribute.
#'
#' This term has no user-visible arguments and is intended to be used directly
#' in \pkg{ergm} formulas or through the ERPM wrapper.
#'
#' @param nw A bipartite \pkg{network} object.
#' @param arglist A list of arguments passed by \pkg{ergm} to the initializer.
#'   This term expects no user-visible arguments, so \code{arglist} should be empty.
#' @param ... Further arguments passed by \pkg{ergm}; not used.
#' @param version ERGM API version. Defaults to \code{packageVersion("ergm")}.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic.
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - This term MUST support multi-toggle proposals (swap/split/merge decomposed
#'   into multiple edge toggles).
#' - Therefore the compiled change-statistic is implemented using the
#'   D_CHANGESTAT_FN API (multi-toggle).
#' - On the R side, we MUST advertise this to ergm by returning `d_func = TRUE`.
#'   Otherwise ergm will assume a one-toggle C_CHANGESTAT_FN and call the
#'   compiled symbol with the wrong signature (=> crash/segfault).
#'
#' Compiled symbol naming convention:
#' - Recommended: implement the C function as `d_log_factorial_sizes` via
#'   D_CHANGESTAT_FN(d_log_factorial_sizes).
#' - Avoid exposing a symbol named `c_log_factorial_sizes` with a D-signature:
#'   ergm may resolve it as the one-toggle entrypoint and crash.
#'
#' The R initializer below:
#' \itemize{
#'   \item enforces that the network is bipartite and undirected;
#'   \item declares that the term is dependent (\code{dependence = TRUE});
#'   \item defines a single scalar statistic whose empty-network value is 0;
#'   \item declares `d_func = TRUE` (multi-toggle changestat).
#' }
#'
#' The ERGM infrastructure will call the C change-statistic whenever a toggle
#' affects an edge between an actor and a group. The C code is responsible for
#' updating the log-factorial sum for the group node at the group mode.
#'
#' @section Mathematical definition:
#' Let:
#' \itemize{
#'   \item \eqn{A} denote the set of actor-mode nodes;
#'   \item \eqn{G} denote the set of group-mode nodes;
#'   \item \eqn{y} the adjacency matrix between actors and groups;
#'   \item \eqn{\deg(g)} the degree of group node \eqn{g \in G}, i.e. the number
#'         of actors in that group.
#' }
#' The statistic is:
#' \deqn{
#'   T(y) = \sum_{g \in G} \log\Gamma(\deg(g)), \quad \log\Gamma(0) := 0.
#' }
#'
#' @section Usage:
#' Typical usage with \pkg{ergm}:
#' \preformatted{
#'   summary(nw ~ log_factorial_sizes)
#'   summary(nw ~ log_factorial_sizes())
#'   ergm(nw ~ log_factorial_sizes)
#' }
#'
#' When using the ERPM wrapper:
#' \preformatted{
#'   erpm(nw ~ log_factorial_sizes)
#'   erpm(partition ~ log_factorial_sizes)
#' }
#'
#' @note
#' The network must be strictly bipartite and undirected:
#' \itemize{
#'   \item the actor mode is identified by \code{nw \%n\% "bipartite"};
#'   \item the group mode is the complement of the actor mode;
#'   \item \code{nw \%n\% "directed"} must not be \code{TRUE}.
#' }
#' If these conditions are not met, the initializer will fail fast via
#' \code{ergm_Init_stop()}.
#'
#' @section Tests:
#' Self-tests for \code{log_factorial_sizes} construct small bipartite networks
#' with known group sizes and compare:
#' \itemize{
#'   \item the ERGM summary \code{summary(nw ~ log_factorial_sizes)};
#'   \item a direct evaluation of \code{sum(lgamma(group_sizes))} with
#'         \code{lgamma(0) <- 0}.
#' }
#' These tests also check that multi-toggle proposals exercise the D_ entrypoint.
#'
#' @keywords ERGM term bipartite groups factorial
#' @md
#'
#' @export
InitErgmTerm.log_factorial_sizes <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "log_factorial_sizes"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.log_factorial_sizes.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.log_factorial_sizes.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[log_factorial_sizes][DEBUG]", ..., "\n", sep = "")

  dbgcat("InitErgmTerm.log_factorial_sizes called with args: ",
         paste(names(arglist), collapse = ", "))

  # Run standard ERGM term checks:
  # - enforce bipartite network;
  # - no user arguments;
  # - let ergm handle all other generic validations.
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,         # no explicit restriction here, see extra guard below
    bipartite     = TRUE,         # require actor-mode/group-mode encoding
    varnames      = character(0), # no user-facing arguments
    vartypes      = character(0),
    defaultvalues = list(),
    required      = logical(0)
  )
  dbgcat("check.ErgmTerm ok")

  # Extra guard for direction: require an undirected bipartite network.
  if (isTRUE(nw %n% "directed")) {
    ergm_Init_stop(sQuote(termname),
                   ": use an undirected bipartite network (actor-group edges).")
  }

  # Single scalar statistic with a fixed, argument-free name.
  coef.names <- termname

  # No numeric INPUT_PARAM passed to the C code for this term.
  inputs <- NULL

  # Statistic on the empty network:
  # - group-mode degrees are all zero;
  # - each contributes lgamma(0) := 0;
  # - sum is therefore 0.
  emptynwstats <- 0

  # ---------------------------------------------------------------------------
  # IMPORTANT: multi-toggle (D_) entrypoint declaration
  # ---------------------------------------------------------------------------
  # - We implement the compiled changestat as D_CHANGESTAT_FN(d_log_factorial_sizes)
  # - Therefore we must return d_func = TRUE so ergm calls the D_ signature.
  dbgcat("returning d_func = TRUE (multi-toggle)")

  list(
    name         = termname,  # must match the registered term name ("log_factorial_sizes")
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,
    d_func       = TRUE,      # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    emptynwstats = emptynwstats
  )
}