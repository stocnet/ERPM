#' ERGM term: cliques (k-actor cliques via group sizes)
#'
#' @file InitErgmTerm.cliques.R
#'
#' @description
#' \code{cliques} is an ERGM term for bipartite actor–group networks that counts
#' k-actor cliques induced by group memberships. The bipartite network is
#' interpreted as:
#' \itemize{
#'   \item an \emph{actor mode}, whose size is given by \code{nw %n% "bipartite"};
#'   \item a \emph{group mode}, consisting of the remaining nodes that represent
#'         groups.
#' }
#'
#' Each group node in the group mode has a degree \eqn{n_g} (number of adjacent
#' actors). Interpreting each group as forming a complete clique among its
#' actors, the total number of k-actor cliques is
#' \deqn{
#'   T_k(y) = \sum_{g \in \text{group mode}} \binom{n_g}{k}.
#' }
#' The term \code{cliques} computes this statistic directly from group sizes,
#' without explicitly materializing the actor–actor projection.
#'
#' The initializer supports:
#' \itemize{
#'   \item \eqn{k \ge 1}, where:
#'     \itemize{
#'       \item for \eqn{k \ge 2}, \eqn{T_k(y)} is the number of k-actor cliques;
#'       \item for \eqn{k = 1}, \eqn{T_1(y)} is the number of groups of size 1
#'             (equivalently, the number of actors that belong to exactly one
#'             singleton group in the projection semantics).
#'     }
#'   \item an optional normalization factor that scales \eqn{T_k(y)} by the
#'         maximum possible \eqn{\binom{N_A}{k}}, where \eqn{N_A} is the
#'         actor-mode size.
#' }
#'
#' @details
#' The term is implemented as a native ERGM C change-statistic, declared in the
#' compiled code under a symbol compatible with \code{name = "cliques"}. The R
#' initializer:
#' \itemize{
#'   \item enforces that the network is bipartite via \code{nw %n% "bipartite"};
#'   \item accepts one or several values of \code{k} (vectorized interface);
#'   \item optionally applies a global normalization by \eqn{\binom{N_A}{k}};
#'   \item packs \code{k} and the associated scaling factors into a compact
#'         \code{INPUT_PARAM} layout for the C layer.
#' }
#'
#' Internally, the actor mode is identified by \code{nw %n% "bipartite"} and the
#' group mode is defined as the complementary set of nodes. The C
#' change-statistic only updates the group node(s) whose membership changes
#' after a toggle, using an un-toggle computation on the degree of the affected
#' group(s).
#'
#' The initializer is vectorized in \code{k}: each entry \eqn{k_j} produces one
#' scalar statistic \eqn{T_{k_j}(y)} (possibly normalized) and one coefficient
#' name.
#'
#' The \code{INPUT_PARAM} vector passed to the C layer has the layout
#' \deqn{
#'   \text{INPUT\_PARAM}
#'   =
#'   (k_1, \text{scale}_1, k_2, \text{scale}_2, \dots, k_J, \text{scale}_J),
#' }
#' where:
#' \itemize{
#'   \item \eqn{k_j} are the requested clique sizes;
#'   \item \eqn{\text{scale}_j} is either 1 (non-normalized case) or
#'         \eqn{1 / \binom{N_A}{k_j}} in the normalized case.
#' }
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
#' For a given integer \eqn{k \ge 1}, define
#' \deqn{
#'   T_k(y)
#'   =
#'   \sum_{g \in G} \binom{n_g}{k}.
#' }
#' When \code{normalized = FALSE}, the ERGM term returns the vector
#' \eqn{(T_{k_1}(y),\dots,T_{k_J}(y))} for all requested values \eqn{k_1,\dots,k_J}.
#'
#' When \code{normalized = TRUE}, let \eqn{N_A = |A|} be the actor-mode size and
#' define the normalized statistic
#' \deqn{
#'   T_k^{\mathrm{norm}}(y)
#'   =
#'   \frac{T_k(y)}{\binom{N_A}{k}}
#' }
#' with the convention that \eqn{1 / \binom{N_A}{k} = 0} whenever
#' \eqn{\binom{N_A}{k} = 0}, so that the corresponding contribution is zero.
#'
#' @section Usage:
#' Typical usage with {ergm} on a bipartite network \code{nw}:
#' \preformatted{
#'   # Count 2-actor cliques induced by groups
#'   summary(nw ~ cliques(k = 2))
#'
#'   # Count 3-actor cliques, normalized by the maximum possible number
#'   summary(nw ~ cliques(k = 3, normalized = TRUE))
#'
#'   # Vectorized k: 2-, 3- and 4-actor cliques in a single term
#'   summary(nw ~ cliques(k = c(2, 3, 4)))
#' }
#'
#' When using the ERPM wrapper on a partition-based workflow, the term can be
#' invoked indirectly as:
#' \preformatted{
#'   erpm(partition ~ cliques(k = 2))
#'   erpm(partition ~ cliques(k = c(2, 3), normalized = TRUE))
#' }
#' provided that the wrapper builds a consistent actor–group bipartite network
#' from the partition.
#'
#' @note
#' The network must be strictly bipartite:
#' \itemize{
#'   \item the actor mode size is given by \code{nw %n% "bipartite"} and must be
#'         a strictly positive integer;
#'   \item the group mode consists of the remaining nodes and represents groups;
#'   \item the \code{cliques} term only depends on degrees of group-mode nodes
#'         and the actor–group incidence.
#' }
#'
#' The initializer is tolerant with respect to the argument name for
#' \code{k}: it accepts positional usage \code{cliques(2)}, the legacy
#' \code{clique_size} name, as well as the explicit \code{k} argument. All
#' values of \code{k} must be integers greater than or equal to 1.
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
#'   # Inspect the number of 2-actor and 3-actor cliques
#'   summary(nw ~ cliques(k = c(2, 3)))
#'
#'   # Normalized version for k = 2
#'   summary(nw ~ cliques(k = 2, normalized = TRUE))
#'
#'   # Fit a simple ERGM with the term
#'   fit <- ergm(nw ~ cliques(k = 2))
#'   summary(fit)
#' }
#'
#' @test
#' Self-tests for \code{cliques} construct small bipartite networks with known
#' group sizes and compare:
#' \itemize{
#'   \item the ERGM summary
#'         \code{summary(nw ~ cliques(k = k_vec, normalized = FALSE))};
#'   \item a direct evaluation of
#'         \eqn{\sum_{g} \binom{n_g}{k}} from the group-mode degrees.
#' }
#' Additional checks verify that:
#' \itemize{
#'   \item the normalized version matches
#'         \eqn{\sum_{g} \binom{n_g}{k} / \binom{N_A}{k}} for each \eqn{k};
#'   \item toggling an actor–group tie changes the statistic by exactly the
#'         increment implied by the local change in the size of the affected
#'         group(s);
#'   \item the vectorized interface over multiple values of \code{k} returns
#'         consistent coefficients and statistics.
#' }
#'
#' @keywords ERGM term bipartite groups cliques
#' @md
#'
#' @export
InitErgmTerm.cliques <- function(nw, arglist, ..., version = packageVersion("ergm")) {
  termname <- "cliques"

  # Normalize user arguments so that the initializer consistently sees 'k':
  # - cliques(1)             -> k = 1
  # - cliques(clique_size=1) -> k = 1
  # - cliques(k=1)           -> k = 1
  if (length(arglist) == 1L) {
    nm <- names(arglist)
    if (is.null(nm) || isTRUE(nm[1L] == "")) {
      # Single positional argument: cliques(1)
      arglist <- list(k = arglist[[1L]])
    }
  }
  if (!is.null(names(arglist)) && "clique_size" %in% names(arglist)) {
    # Backward compatibility: accept 'clique_size' and rename it to 'k'
    arglist[["k"]] <- arglist[["clique_size"]]
    arglist[["clique_size"]] <- NULL
  }

  # Retrieve the actor-mode size N_A from the bipartite attribute.
  # This is the number of actors; the remaining nodes form the group mode.
  n1 <- network::get.network.attribute(nw, "bipartite")
  if (is.null(n1) || is.na(n1)) {
    ergm_Init_stop(sQuote(termname), ": réseau non biparti ou attribut 'bipartite' manquant.")
  }
  n1 <- as.integer(n1)
  if (n1 <= 1L) ergm_Init_stop(sQuote(termname), ": biparti invalide (N1 <= 1).")

  # Run standard ERGM term checks:
  # - enforce bipartite network;
  # - expect arguments 'k' and 'normalized';
  # - let {ergm} handle generic validations.
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("k", "normalized"),
    vartypes      = c("numeric", "logical"),
    defaultvalues = list(2, FALSE),
    required      = c(FALSE, FALSE)
  )

  k  <- a$k
  nz <- a$normalized

  # Validate k and normalized:
  # - require at least one value of k;
  # - enforce integer k >= 1;
  # - enforce scalar logical normalized.
  if (length(k) < 1L)
    ergm_Init_stop(sQuote(termname), ": specify at least one value of k.")
  if (any(k != as.integer(k) | k < 1))
    ergm_Init_stop(sQuote(termname), ": 'k' must contain integers >= 1.")
  if (length(nz) != 1L || is.na(nz))
    ergm_Init_stop(sQuote(termname), ": 'normalized' must be a scalar boolean.")

  # Compute normalization factors:
  # - if normalized = FALSE: scale = 1;
  # - if normalized = TRUE : scale = 1 / choose(N_A, k),
  #   with protection against zero denominators.
  scale <- rep(1, length(k))
  if (isTRUE(nz)) {
    denom <- choose(n1, k)              # also valid for k = 1 (denom = N_A)
    denom[denom == 0 | is.na(denom)] <- Inf
    scale <- 1 / denom
  }

  # Build coefficient names and INPUT_PARAM for the C layer:
  # - one coefficient per k;
  # - INPUT_PARAM = (k_1, scale_1, k_2, scale_2, ...).
  coef.names <- if (isTRUE(nz)) paste0("cliques_k", k, "_norm") else paste0("cliques_k", k)
  inputs <- c(rbind(as.integer(k), as.double(scale)))

  # Return the ERGM term specification expected by {ergm}.
  # The field 'name' must match the C change-statistic symbol 'cliques'.
  list(
    name         = "cliques",
    coef.names   = coef.names,
    inputs       = inputs,
    dependence   = TRUE,
    emptynwstats = numeric(length(k))
  )
}