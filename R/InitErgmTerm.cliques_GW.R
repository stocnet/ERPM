#' ERGM term: cliques_GW (geometrically weighted group sizes)
#'
#' @note InitErgmTerm.cliques_GW.R
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
#' The term is implemented as a native ERGM C change-statistic, declared in the
#' compiled code under a symbol compatible with \code{name = "cliques_GW"}. The
#' R initializer:
#' \itemize{
#'   \item enforces that the network is bipartite via \code{check.ErgmTerm()};
#'   \item accepts one or several values of \code{lambda};
#'   \item validates the numerical domain of \code{lambda};
#'   \item precomputes the ratio \eqn{r_\lambda = (\lambda-1)/\lambda};
#'   \item packs \code{lambda} and \code{r_\lambda} into a compact
#'         \code{INPUT_PARAM} layout for the C layer.
#' }
#'
#' The \code{INPUT_PARAM} vector passed to the C layer has the layout
#' \deqn{
#'   \text{INPUT\_PARAM}
#'   =
#'   (\lambda_1, r_{\lambda_1},
#'    \lambda_2, r_{\lambda_2},
#'    \dots,
#'    \lambda_J, r_{\lambda_J}),
#' }
#' where \eqn{\lambda_j} are the user-specified parameters and
#' \eqn{r_{\lambda_j} = (\lambda_j - 1)/\lambda_j}.
#'
#' For \eqn{\lambda = 1}, the definition \eqn{r_\lambda = 0} implies that
#' \eqn{S(n_g, 1) = 1} for any \eqn{n_g \ge 1}, and \eqn{S(0, 1) = 0}, so the
#' statistic reduces to the number of non-empty groups in the group mode.
#'
#' @section Mathematical definition:
#' Let:
#' \itemize{
#'   \item \eqn{A} be the set of actor-mode nodes;
#'   \item \eqn{G} be the set of group-mode nodes;
#'   \item \eqn{B} the bipartite adjacency matrix between actors and groups,
#'         with \eqn{B_{ig} = 1} if actor \eqn{i} belongs to group \eqn{g};
#'   \item \eqn{A_g = \{ i \in A : B_{ig} = 1 \}} the set of actors in group
#'         \eqn{g};
#'   \item \eqn{n_g = |A_g|} the size of group \eqn{g}.
#' }
#' For a given \eqn{\lambda \ge 1}, define
#' \deqn{
#'   S(n_g, \lambda)
#'   =
#'   \lambda \Big[ 1 - \Big(\frac{\lambda - 1}{\lambda}\Big)^{n_g} \Big],
#' }
#' and the statistic
#' \deqn{
#'   T_\lambda(y)
#'   =
#'   \sum_{g \in G} S(n_g, \lambda).
#' }
#' When multiple values \eqn{\lambda_1,\dots,\lambda_J} are supplied, the ERGM
#' term returns the vector
#' \eqn{(T_{\lambda_1}(y), \dots, T_{\lambda_J}(y))}.
#'
#' @section Usage:
#' Typical usage with {ergm} on a bipartite actor–group network \code{nw}:
#' \preformatted{
#'   # Single lambda
#'   summary(nw ~ cliques_GW(lambda = 2))
#'
#'   # Multiple lambda values (vectorized)
#'   summary(nw ~ cliques_GW(lambda = c(1.5, 2, 4)))
#'
#'   # Fit an ERGM with a cliques_GW term
#'   fit <- ergm(nw ~ cliques_GW(lambda = 2))
#'   summary(fit)
#' }
#'
#' When using the ERPM wrapper on a partition-based workflow, the term can be
#' invoked indirectly as:
#' \preformatted{
#'   erpm(partition ~ cliques_GW(lambda = 2))
#'   erpm(partition ~ cliques_GW(lambda = c(2, 3)))
#' }
#' provided that the wrapper builds a consistent actor–group bipartite network
#' from the partition.
#'
#' @note
#' The network is interpreted as a bipartite actor–group graph:
#' \itemize{
#'   \item the actor mode size is given by \code{nw \%n\% "bipartite"} and must
#'         be a strictly positive integer;
#'   \item the group mode consists of the remaining nodes and represents groups;
#'   \item the \code{cliques_GW} statistic depends only on the degrees of
#'         group-mode nodes and the actor–group incidence structure.
#' }
#'
#' The initializer is vectorized in \code{lambda}. Domain checks ensure that
#' \code{lambda} is finite and at least 1. Values \eqn{\lambda > 1} yield
#' \eqn{r_\lambda \in (0, 1)} and therefore a strictly decaying geometric
#' profile as group size increases.
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
#'   # Adjacency matrix: actors 1..4, groups 5..7
#'   adj <- matrix(0, n_total, n_total)
#'
#'   # Group 5: actors {1, 2}
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'
#'   # Group 6: actors {2, 3, 4}
#'   adj[2, 6] <- adj[6, 2] <- 1
#'   adj[3, 6] <- adj[6, 3] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'
#'   # Group 7: singleton {1}
#'   adj[1, 7] <- adj[7, 1] <- 1
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw %n% "bipartite" <- n_actors  # actor mode size
#'
#'   # Inspect the geometric cliques statistic for different lambdas
#'   summary(nw ~ cliques_GW(lambda = c(1.5, 2, 4)))
#'
#'   # Fit a simple ERGM with the term
#'   fit <- ergm(nw ~ cliques_GW(lambda = 2))
#'   summary(fit)
#' }
#'
#' @section Tests:
#' Self-tests for \code{cliques_GW} construct small bipartite actor–group
#' networks with known group sizes, then:
#' \itemize{
#'   \item compute group-mode degrees \eqn{n_g} and evaluate
#'         \eqn{T_\lambda(y) = \sum_g \lambda \big[1 - r_\lambda^{n_g}\big]}
#'         directly in R for several values of \eqn{\lambda};
#'   \item compare these reference values with
#'         \code{summary(nw ~ cliques_GW(lambda = lambda_vec))};
#'   \item check that toggling an actor–group tie modifies the statistic by the
#'         increment predicted by the local change in the degree of the affected
#'         group(s).
#' }
#'
#' @keywords ERGM term bipartite groups cliques geometric
#' @md
#'
#' @export
InitErgmTerm.cliques_GW <- function(nw, arglist, ..., version = packageVersion("ergm")) {
    termname <- "cliques_GW"

    # Run standard ERGM term checks:
    # - enforce bipartite network;
    # - expect a numeric 'lambda' argument (scalar or vector);
    # - let {ergm} handle generic validations.
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
    if (length(lambda) == 0L) return(NULL)

    # Domain validation for lambda:
    # - finite values only;
    # - lambda >= 1 (lambda > 1 gives r in (0, 1), lambda = 1 gives r = 0).
    if (any(!is.finite(lambda)))
        ergm_Init_stop(sQuote(termname), ": 'lambda' must be finite (no NA/NaN/Inf).")
    if (any(lambda < 1))
        ergm_Init_stop(sQuote(termname), ": 'lambda' must be >= 1.")

    # Precompute r = (lambda - 1) / lambda on the R side.
    # This is passed to the C code alongside the original lambda values.
    r <- (lambda - 1) / lambda

    # Helper to generate readable coefficient names, e.g.:
    #   cliques_GW_lambda2, cliques_GW_lambda1.5, cliques_GW_lambda0.125
    pretty_num <- function(x) {
        s <- formatC(x, digits = 4, format = "fg", flag = "#")
        sub("\\.$", "", s)  # strip trailing "." if present ("2." -> "2")
    }
    coef.names <- paste0("cliques_GW_lambda", vapply(lambda, pretty_num, ""))

    # Pack INPUT_PARAM as interleaved (lambda_j, r_j) pairs, column-wise.
    # Layout:
    #   [2*j + 0] = lambda_j
    #   [2*j + 1] = r_j
    inputs <- c(rbind(as.double(lambda), as.double(r)))

    # Return the ERGM term specification expected by {ergm}.
    # The field 'name' must match the C change-statistic symbol 'cliques_GW'.
    list(
        name         = "cliques_GW",
        coef.names   = coef.names,
        inputs       = inputs,
        dependence   = TRUE,
        emptynwstats = numeric(length(lambda))
    )
}
