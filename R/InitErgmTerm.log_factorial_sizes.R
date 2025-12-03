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
#' in {ergm} formulas or through the ERPM wrapper.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic, declared in
#' the compiled code under the name \code{c_log_factorial_sizes}. The R
#' initializer below:
#' \itemize{
#'   \item enforces that the network is bipartite and undirected;
#'   \item declares that the term is dependent (\code{dependence = TRUE});
#'   \item defines a single scalar statistic whose empty-network value is 0.
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
#' Typical usage with {ergm}:
#' \preformatted{
#'   summary(nw ~ log_factorial_sizes)
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
#' @examples
#' \dontrun{
#'   library(network)
#'   library(ergm)
#'
#'   # Build a small bipartite network: 4 actors, 3 groups
#'   n_actors <- 4
#'   n_groups <- 3
#'   n_total  <- n_actors + n_groups
#'
#'   adj <- matrix(0, n_total, n_total)
#'   # Actors = 1..4, Groups = 5..7
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'   adj[3, 6] <- adj[6, 3] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'   # Group 7 remains empty
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw \%n\% "bipartite" <- n_actors  # actor mode size
#'
#'   # Inspect the statistic
#'   summary(nw ~ log_factorial_sizes)
#'
#'   # Fit a simple ERGM with the term
#'   fit <- ergm(nw ~ log_factorial_sizes)
#'   summary(fit)
#' }
#'
#' @section Tests:
#' Self-tests for \code{log_factorial_sizes} construct small bipartite networks
#' with known group sizes and compare:
#' \itemize{
#'   \item the ERGM summary \code{summary(nw ~ log_factorial_sizes)};
#'   \item a direct evaluation of \code{sum(lgamma(group_sizes))} with
#'         \code{lgamma(0) <- 0}.
#' }
#' These tests also check that toggling an actor–group tie updates the statistic
#' by the correct local increment as defined in the C change-statistic.
#'
#' @keywords ERGM term bipartite groups factorial
#' @md
#'
#' @export
InitErgmTerm.log_factorial_sizes <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "log_factorial_sizes"

  # Run standard ERGM term checks:
  # - enforce bipartite network;
  # - no user arguments;
  # - let {ergm} handle all other generic validations.
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,         # no explicit restriction here, see extra guard below
    bipartite     = TRUE,         # require a bipartite representation (actor mode / group mode)
    varnames      = character(0), # no user-facing arguments
    vartypes      = character(0),
    defaultvalues = list(),
    required      = logical(0)
  )

  # Extra guard for direction: require an undirected bipartite network.
  # This keeps the semantics consistent with "actors in groups" on an undirected bipartite graph.
  if (isTRUE(nw %n% "directed"))
    ergm_Init_stop(sQuote(termname), ": use an undirected bipartite network (actor–group edges).")

  # Single scalar statistic with a fixed, argument-free name.
  coef.names   <- termname

  # No numeric INPUT_PARAM passed to the C code for this term.
  inputs       <- NULL

  # Statistic on the empty network:
  # - group-mode degrees are all zero;
  # - each contributes lgamma(0) := 0;
  # - sum is therefore 0.
  emptynwstats <- 0

  # Return the ERGM term specification expected by {ergm}.
  # The field `name` must match the C symbol `c_log_factorial_sizes`.
  list(
    name         = termname,   # must match c_log_factorial_sizes in the C changestats
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,       # statistic depends on the network configuration
    emptynwstats = emptynwstats
  )
}
