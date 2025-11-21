# ==============================================================================
# File    : InitErgmTerm.squared_sizes.R
# Purpose : Declare the ERGM term 'squared_sizes' for bipartite networks
#           (aggregation over the group mode).
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM term: squared_sizes (group-mode degrees raised to a power)
#'
#' @file InitErgmTerm.squared_sizes.R
#'
#' @description
#' \code{squared_sizes} is an ERGM term for bipartite networks that aggregates
#' powers of group sizes. The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each node in the group mode with degree \eqn{\deg(g)} (number of adjacent
#' actors), the term can accumulate contributions of the form
#' \deqn{
#'   \deg(g)^{\text{pow}},
#' }
#' restricted to degree intervals \eqn{[\text{from}, \text{to})}. The term is
#' vectorized: providing vectors \code{from}, \code{to} and \code{pow} yields
#' one statistic component per triplet \eqn{(\text{from}_k, \text{to}_k, \text{pow}_k)}.
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic, declared in the
#' compiled code under the name \code{c_squared_sizes}. The R initializer:
#' \itemize{
#'   \item enforces that the network is bipartite (actor mode vs group mode);
#'   \item parses and recycles the user arguments \code{from}, \code{to} and \code{pow};
#'   \item checks that \code{pow} is integer and at least 1;
#'   \item replaces infinite upper bounds \code{to = Inf} by a finite value
#'         while preserving the \eqn{[\text{from}, \text{to})} semantics;
#'   \item builds a compact \code{inputs} vector for the C code:
#'         \code{c(from_1, to_1, pow_1, ..., from_K, to_K, pow_K)}.
#' }
#'
#' The ERGM infrastructure will call the C change-statistic whenever a toggle
#' affects an actor–group edge. On the C side, the code is responsible for
#' updating each component \eqn{T_k(y)} based on how the degrees of group-mode
#' nodes change.
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
#' For each component \eqn{k} (corresponding to \code{from[k]}, \code{to[k]},
#' \code{pow[k]}), the statistic is
#' \deqn{
#'   T_k(y)
#'   =
#'   \sum_{g \in G}
#'     \mathbf{1}\big[\deg(g) \in [\text{from}_k, \text{to}_k)\big]\,
#'     \deg(g)^{\text{pow}_k}.
#' }
#' The overall statistic is the vector \eqn{(T_1(y), \dots, T_K(y))}.
#'
#' @section Usage:
#' Typical usage with {ergm} (constraints explicit):
#' \preformatted{
#'   # Sum of squared group sizes over all groups
#'   summary(nw ~ squared_sizes, constraints = ~ b1part)
#'
#'   # Groups with size in [2, 5), cubic power
#'   summary(
#'     nw ~ squared_sizes(from = 2, to = 5, pow = 3),
#'     constraints = ~ b1part
#'   )
#' }
#'
#' When using the ERPM wrapper:
#' \preformatted{
#'   # On an explicit bipartite network
#'   erpm(nw ~ squared_sizes)
#'
#'   # From a partition vector (actor -> group id)
#'   erpm(partition ~ squared_sizes)
#' }
#'
#' @param from numeric, integer(s) \eqn{\ge 0}. Lower bounds of the degree
#'   intervals, inclusive. Defaults to \code{1}. Can be a scalar or a vector.
#' @param to numeric, integer(s). Upper bounds of the degree intervals,
#'   exclusive. \code{Inf} is allowed to represent "no upper bound".
#'   Defaults to \code{Inf}. Can be a scalar or a vector.
#' @param pow numeric, integer(s) \eqn{\ge 1}. Powers applied to the group
#'   degrees. Defaults to \code{2}. Can be a scalar or a vector.
#'
#' @return
#' A standard {ergm} term initialization list with components:
#' \itemize{
#'   \item \code{name}         = \code{"squared_sizes"};
#'   \item \code{coef.names}   = character vector, one name per (from,to,pow) triplet;
#'   \item \code{inputs}       = integer vector \code{c(from_1, to_1, pow_1, ..., from_K, to_K, pow_K)};
#'   \item \code{dependence}   = \code{TRUE};
#'   \item \code{emptynwstats} = numeric vector of zeros (one per component).
#' }
#'
#' @note
#' \itemize{
#'   \item The network must be bipartite and interpreted as actors versus groups.
#'   \item The actor mode is identified by \code{nw \%n\% "bipartite"}.
#'   \item The group mode is the complement of the actor mode.
#'   \item The initializer enforces bipartiteness via \code{check.ErgmTerm(..., bipartite = TRUE)}.
#'   \item The term is declared as dependent (\code{dependence = TRUE}), since it
#'         depends on the current pattern of actor–group memberships.
#' }
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
#'   nw %n% "bipartite" <- n_actors  # actor-mode size
#'
#'   # Sum of squared group sizes over all groups
#'   summary(nw ~ squared_sizes, constraints = ~ b1part)
#'
#'   # Vectorized example: two components
#'   #  - squared sizes for groups of size in [1, +Inf)
#'   #  - cubic sizes for groups of size in [3, 6)
#'   summary(
#'     nw ~ squared_sizes(from = c(1, 3),
#'                        to   = c(Inf, 6),
#'                        pow  = c(2, 3)),
#'     constraints = ~ b1part
#'   )
#'
#'   # ERGM fit using squared_sizes
#'   fit <- ergm(nw ~ squared_sizes, constraints = ~ b1part)
#'   summary(fit)
#' }
#'
#' @test
#' Self-tests for \code{squared_sizes} construct small bipartite networks with
#' known group degrees and compare:
#' \itemize{
#'   \item the ERGM summary \code{summary(nw ~ squared_sizes(...), constraints = ~ b1part)};
#'   \item a direct computation of
#'         \code{colSums( (deg_g %in% [from_k, to_k)) * deg_g^pow_k )} for each component;
#' }
#' They also verify that toggling a single actor–group edge changes each
#' component \eqn{T_k(y)} by the expected local increment, in agreement with
#' the C change-statistic \code{c_squared_sizes}.
#'
#' @keywords ERGM term bipartite groups degree power
#' @md
#'
#' @export
InitErgmTerm.squared_sizes <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "squared_sizes"

  # ---------------------------------------------------------------------------
  # Base validation and structural requirements
  # ---------------------------------------------------------------------------
  # - Enforce a bipartite network (actor mode / group mode).
  # - Declare the allowed arguments: from, to, pow.
  # - Provide default values: from = 1, to = Inf, pow = 2.
  a <- check.ErgmTerm(
    nw, arglist,
    directed   = NULL,
    bipartite  = TRUE,                     # require a bipartite actor–group representation
    varnames   =    c("from",   "to",      "pow"),
    vartypes   =    c("numeric","numeric", "numeric"),
    defaultvalues = list(1,      Inf,      2),   # default: all groups, squared sizes
    required   =    c(FALSE,    FALSE,     FALSE)
  )

  from <- a$from
  to   <- a$to
  pow  <- a$pow

  # ---------------------------------------------------------------------------
  # Vector length harmonization
  # ---------------------------------------------------------------------------
  # Recycle scalar arguments so that from, to and pow have the same length.
  # Examples:
  #   squared_sizes(from = c(1, 2, 4), to = 5,      pow = 2)
  #   squared_sizes(from = 1,          to = c(5, 2), pow = 2)
  #   squared_sizes(from = c(1, 2, 4), to = c(5, 2, 4), pow = 2)
  if (length(to)   == 1L && length(from) > 1L) to   <- rep(to,   length(from))
  if (length(from) == 1L && length(to)   > 1L) from <- rep(from, length(to))
  if (length(pow)  == 1L && length(from) > 1L) pow  <- rep(pow,  length(from))

  if (length(from) != length(to) || length(pow) != length(from)) {
    ergm_Init_stop(sQuote(termname), ": incompatible lengths for from/to/pow.")
  }

  # ---------------------------------------------------------------------------
  # Constraints on pow
  # ---------------------------------------------------------------------------
  # Require pow >= 1 and integer-valued.
  if (any(pow < 1 | pow != as.integer(pow))) {
    ergm_Init_stop(sQuote(termname), ": 'pow' must be an integer >= 1.")
  }

  # ---------------------------------------------------------------------------
  # Replace Inf in 'to' with a finite upper bound
  # ---------------------------------------------------------------------------
  # To keep the [from, to) semantics while providing a finite upper bound to the
  # C code, replace Inf with:
  #   max(from, network.size(nw)) + 1
  # so that degrees smaller than this value behave as "no upper limit" in
  # practice.
  to <- ifelse(
    is.infinite(to),
    pmax(from, network.size(nw)) + 1L,
    to
  )
  # Historical alternative (commented out, not used anymore):
  # n1 <- network.bipartite(nw)  # number of actor-mode vertices
  # to <- ifelse(is.infinite(to), pmax(from, n1) + 1, to)

  # ---------------------------------------------------------------------------
  # Interval sanity checks
  # ---------------------------------------------------------------------------
  # Each component must define a valid interval [from_k, to_k) with from_k < to_k.
  if (any(from >= to)) {
    ergm_Init_stop(sQuote(termname), ": from < to is required for all components.")
  }

  # ---------------------------------------------------------------------------
  # Trivial case: no intervals => no statistic
  # ---------------------------------------------------------------------------
  if (length(from) == 0L) return(NULL)

  # ---------------------------------------------------------------------------
  # Coefficient names
  # ---------------------------------------------------------------------------
  # Build readable coefficient names, encoding:
  # - the degree interval [from, to) or [from, +Inf);
  # - the power pow, when pow != 1.
  coef.names <- ifelse(
    to >= network.size(nw) + 1L,
    paste0(
      "squared_sizes_", from, "+",
      ifelse(pow != 1, paste0("_pow", pow), "")
    ),
    paste0(
      "squared_sizes_", from, "to", to,
      ifelse(pow != 1, paste0("_pow", pow), "")
    )
  )

  # ---------------------------------------------------------------------------
  # Compact INPUT_PARAM vector for the C change-statistic
  # ---------------------------------------------------------------------------
  # Pack (from_k, to_k, pow_k) by rows, then flatten:
  #   inputs = c(from_1, to_1, pow_1, from_2, to_2, pow_2, ...)
  inputs <- c(rbind(as.integer(from), as.integer(to), as.integer(pow)))

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  # - name         : must match the C symbol (without the 'c_' prefix parsed by {ergm}).
  # - coef.names   : one label per (from,to,pow) triplet.
  # - inputs       : numeric vector passed to the C change-statistic as INPUT_PARAM.
  # - dependence   : TRUE since the term depends on the network configuration.
  # - emptynwstats : empty network -> all components equal to 0.
  list(
    name         = "squared_sizes",          # must match CHANGESTAT_FN(c_squared_sizes) (without 'c_')
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,
    emptynwstats = numeric(length(from))     # empty network => 0 for each component
  )
}