# ==============================================================================
# File    : R/InitErgmTerm.dyadcov_GW.R
# Term    : dyadcov_GW
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM term: dyadcov_GW (geometrically weighted dyadic covariate)
#'
#' @name InitErgmTerm.dyadcov_GW
#' @aliases dyadcov_GW
#' @note InitErgmTerm.dyadcov_GW.R
#'
#' @description
#' \code{dyadcov_GW} is an ERGM term for bipartite networks that builds a
#' geometrically weighted aggregation of a dyadic covariate
#' \eqn{Z = (z_{ij})} over actor cliques within each group.
#'
#' The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each group in the group mode, we consider its incident actors, all
#' \eqn{k}-subsets of these actors, and the products of \emph{symmetrised}
#' dyadic covariates \eqn{z_{ij} + z_{ji}} over each subset. These clique
#' contributions are then combined across clique sizes \eqn{k \ge 2} with a
#' geometric weight controlled by \eqn{\lambda}.
#'
#' @details
#' The core statistic for a given partition \eqn{p}, dyadic covariate \eqn{Z},
#' and scalar \eqn{\lambda > 0} is:
#'
#' @section Definition:
#'   T_GW(p; Z, lambda)
#'     = sum_g sum_{k = 2..n_g} a_k(lambda) * S_g^{(k)}(Z)
#'     = sum_g sum_{k = 2..n_g} (-1/lambda)^{k-1}
#'                      * sum_{C in C_k(g)} prod_{i<j in C} (z_{ij} + z_{ji})
#'
#' where:
#' \itemize{
#'   \item \eqn{g} ranges over group-mode nodes (groups);
#'   \item \eqn{n_g} is the number of actors attached to group \eqn{g};
#'   \item \eqn{C_k(g)} is the set of all actor subsets \eqn{C} of size \eqn{k}
#'         within group \eqn{g};
#'   \item \eqn{S_g^{(k)}(Z)} is the raw \eqn{k}-clique sum
#'         \eqn{\sum_{C \in C_k(g)} \prod_{i<j \in C} (z_{ij}+z_{ji})};
#'   \item \eqn{a_k(\lambda) = (-1/\lambda)^{k-1}} is the geometric weight
#'         applied to clique size \eqn{k}.
#' }
#'
#' The dyadic covariate \eqn{Z} is a real actor-by-actor matrix of
#' size \eqn{n_A \times n_A}, where \eqn{n_A} is the size of the actor mode
#' \code{n1 = nw \%n\% "bipartite"}. Symmetry is \strong{not} required:
#' both \eqn{z_{ij}} and \eqn{z_{ji}} are used for each unordered pair
#' \eqn{\{i,j\}} when available. The initializer:
#' \itemize{
#'   \item enforces that the network is bipartite and retrieves the actor-mode size;
#'   \item accepts \code{dyadcov} as either a literal matrix or the name of a
#'         network-level attribute (\code{nw \%n\% "..."});
#'   \item truncates \code{dyadcov} to its top-left \code{n1 x n1} block if
#'         larger than required;
#'   \item checks that \code{dyadcov} is numeric and free of \code{NA};
#'   \item optionally reports (in debug mode) how far the matrix is from being
#'         symmetric, but does not require symmetry.
#' }
#'
#' @section Implementation and change-statistic:
#' The term is implemented as a native ERGM C change-statistic, exposed under the
#' symbol \code{c_dyadcov_GW}. The R initializer below:
#' \itemize{
#'   \item packages \code{n1}, \code{lambda} and the flattened \code{dyadcov}
#'         matrix into \code{INPUT_PARAM};
#'   \item declares the term as dependent (\code{dependence = TRUE}) with
#'         unbounded support (\code{minval = -Inf}, \code{maxval = Inf});
#'   \item sets the empty-network statistic to \code{0}.
#' }
#'
#' On each toggle of an actor–group edge, the C change-statistic recomputes the
#' local contribution for the unique group touched by the toggle, updating the
#' geometrically weighted sum over clique contributions
#' \eqn{\{S_g^{(k)}(Z)\}_{k \ge 2}} based on the symmetrised dyadic covariate
#' \eqn{z_{ij}+z_{ji}}.
#'
#' @section Arguments:
#' The initializer is not called directly by users; it is invoked automatically
#' by {ergm} when the term \code{dyadcov_GW(...)} appears on the right-hand
#' side of a model formula.
#'
#' @param dyadcov matrix or character. Either:
#'   \itemize{
#'     \item a numeric matrix of size at least \code{n1 x n1}, where
#'           \code{n1 = nw \%n\% "bipartite"} is the actor-mode size; or
#'     \item the name of a network-level attribute containing such a matrix
#'           (retrieved as \code{nw \%n\% dyadcov}).
#'   }
#'   In both cases the matrix is truncated, if necessary, to its top-left
#'   \code{n1 x n1} block.
#' @param lambda numeric scalar. Geometric weight parameter \eqn{\lambda > 0}.
#'   Typical usage sets \eqn{\lambda > 1} to ensure decaying weights in
#'   \eqn{a_k(\lambda) = (-1/\lambda)^{k-1}}, but the initializer only requires
#'   \eqn{\lambda > 0} and finiteness.
#'
#' @return
#' A standard {ergm} term initialization list with components:
#' \itemize{
#'   \item \code{name}         = \code{"dyadcov_GW"};
#'   \item \code{coef.names}   = a single coefficient name encoding the
#'         dyadic covariate label and the value of \code{lambda};
#'   \item \code{inputs}       = numeric vector
#'         \code{c(n1, lambda, as.double(Z))};
#'   \item \code{dependence}   = \code{TRUE};
#'   \item \code{minval}       = \code{-Inf};
#'   \item \code{maxval}       = \code{Inf};
#'   \item \code{emptynwstats} = \code{0}.
#' }
#'
#' @note
#' \itemize{
#'   \item The network must be bipartite and interpreted as actors versus groups.
#'   \item The actor mode is identified by \code{nw \%n\% "bipartite"} and must
#'         be a strictly positive integer.
#'   \item The dyadic matrix must be numeric and free of \code{NA} values.
#'         Any \code{NA} triggers a fail-fast error. Symmetry is \emph{not}
#'         required: the C code always uses the symmetrised combination
#'         \code{z_ij + z_ji} for each unordered pair \code{\{i,j\}}.
#'   \item Debug logging is controlled by the global option
#'         \code{options(ERPM.dyadcov_GW.debug = TRUE/FALSE)}. When \code{TRUE},
#'         the initializer prints diagnostic messages prefixed by
#'         \code{"[dyadcov_GW][DEBUG]"} and may report how far the dyadic
#'         matrix is from symmetry.
#' }
#'
#' @examples
#' \dontrun{
#'   library(network)
#'   library(ergm)
#'
#'   # -------------------------------------------------------------------------
#'   # Build a toy bipartite network: 4 actors, 2 groups
#'   # -------------------------------------------------------------------------
#'   n_actors <- 4
#'   n_groups <- 2
#'   n_total  <- n_actors + n_groups
#'
#'   adj <- matrix(0, n_total, n_total)
#'
#'   # Actors = 1..4, Groups = 5..6
#'   # Group 5 has actors 1, 2, 3
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'   adj[3, 5] <- adj[5, 3] <- 1
#'   # Group 6 has actors 3, 4
#'   adj[3, 6] <- adj[6, 3] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw %n% "bipartite" <- n_actors  # actor-mode size
#'
#'   # -------------------------------------------------------------------------
#'   # Build a (possibly non-symmetric) dyadic covariate on actors
#'   # -------------------------------------------------------------------------
#'   Z <- matrix(
#'     c(
#'       0, 1, 2, 3,
#'       10, 0, 4, 5,
#'       20, 40, 0, 6,
#'       30, 50, 60, 0
#'     ),
#'     n_actors, n_actors, byrow = TRUE
#'   )
#'   nw %n% "Z_example" <- Z
#'
#'   # -------------------------------------------------------------------------
#'   # Inspect the dyadcov_GW statistic
#'   # -------------------------------------------------------------------------
#'   # Geometrically weighted aggregation with lambda = 2
#'   summary(
#'     nw ~ dyadcov_GW("Z_example", lambda = 2),
#'     constraints = ~ b1part
#'   )
#' }
#'
#' @section Tests:
#' Self-tests for \code{dyadcov_GW} (not shown here) typically:
#' \itemize{
#'   \item construct small bipartite networks with a known partition of actors
#'         into groups;
#'   \item define dyadic matrices \code{Z} (symmetric or not) where clique
#'         products based on \code{z_ij + z_ji} and their geometric weights
#'         \eqn{(-1/\lambda)^{k-1}} can be computed explicitly;
#'   \item compare \code{summary(nw ~ dyadcov_GW(...), constraints = ~ b1part)}
#'         with a direct R implementation of
#'         \code{sum_g sum_{k >= 2} (-1/lambda)^{k-1} * S_g^{(k)}(Z)};
#'   \item apply explicit edge toggles and check that observed changes match
#'         the Δ produced by \code{c_dyadcov_GW}.
#' }
#'
#' @keywords ERGM term bipartite dyadic covariate geometrically-weighted cliques
#' @md
InitErgmTerm.dyadcov_GW <- function(nw, arglist, ...) {
  termname <- "dyadcov_GW"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.dyadcov_GW.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.dyadcov_GW.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[dyadcov_GW][DEBUG]", ..., "\n", sep = "")

  # ---------------------------------------------------------------------------
  # Base ERGM term validation and argument parsing
  # ---------------------------------------------------------------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("dyadcov",              "lambda"),
    vartypes      = c("matrix,character",     "numeric"),
    defaultvalues = list(NULL,                2),
    required      = c(TRUE,                   FALSE)
  )

  # ---------------------------------------------------------------------------
  # Actor-mode size (n1) from the bipartite attribute
  # ---------------------------------------------------------------------------
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L)
    stop(termname, ": réseau non biparti strict (attribut %n% 'bipartite' manquant ou invalide).")
  dbgcat("n1 = ", n1)

  # ---------------------------------------------------------------------------
  # Retrieve the dyadic matrix (n1 x n1)
  # ---------------------------------------------------------------------------
  dyad_raw   <- a$dyadcov
  dyad_label <- NULL

  if (is.character(dyad_raw) && length(dyad_raw) == 1L) {
    # Case: name of a network-level attribute storing the matrix
    dyad_mat   <- nw %n% dyad_raw
    dyad_label <- dyad_raw
    if (is.null(dyad_mat))
      stop(termname, ": attribut de niveau réseau inexistant: ", sQuote(dyad_raw), ".")
    dbgcat("dyadcov source = network attribute ", sQuote(dyad_label))
  } else {
    # Case: matrix passed literally as an argument
    dyad_mat   <- dyad_raw
    dyad_label <- "dyadcov"
    dbgcat("dyadcov source = matrix literal")
  }

  if (!is.matrix(dyad_mat))
    stop(termname, ": 'dyadcov' doit être une matrice ou le nom d'un attribut de niveau réseau.")

  nr <- nrow(dyad_mat)
  nc <- ncol(dyad_mat)

  if (nr < n1 || nc < n1)
    stop(termname, ": dimensions de la matrice dyadique (", nr, "x", nc,
         ") insuffisantes pour n1 = ", n1, ".")

  # If larger than needed, restrict to the top-left n1 x n1 block
  if (nr > n1 || nc > n1) {
    dyad_mat <- dyad_mat[seq_len(n1), seq_len(n1), drop = FALSE]
    dbgcat("dyadcov truncated to ", n1, "x", n1)
  }

  # ---------------------------------------------------------------------------
  # Numeric coercion and fail-fast on NA
  # ---------------------------------------------------------------------------
  if (!is.numeric(dyad_mat))
    stop(termname, ": la matrice dyadique doit être numérique.")

  if (anyNA(dyad_mat))
    stop(termname, ": NA non autorisé dans la matrice dyadique.")

  # ---------------------------------------------------------------------------
  # Optional symmetry diagnostics (non-blocking)
  # ---------------------------------------------------------------------------
  tol <- 1e-8
  max_asym <- max(abs(dyad_mat - t(dyad_mat)))
  if (dbg && max_asym > tol) {
    dbgcat("warning: dyadcov matrix is not symmetric, max |Z - t(Z)| = ",
           signif(max_asym, 5L))
  }

  dbgcat("dyadcov dim = ", paste(dim(dyad_mat), collapse = "x"),
         " | sample = ",
         paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L), collapse = ","))

  # ---------------------------------------------------------------------------
  # lambda (geometric weight parameter)
  # ---------------------------------------------------------------------------
  lambda_raw <- a$lambda
  if (is.null(lambda_raw) || length(lambda_raw) == 0L) {
    lambda <- 2
  } else {
    if (!is.numeric(lambda_raw))
      stop(termname, ": 'lambda' doit être numérique.")
    lambda <- as.numeric(lambda_raw[1L])
  }

  if (!is.finite(lambda) || lambda <= 0)
    stop(termname, ": 'lambda' doit être un réel strictement positif (et typiquement > 1).")

  dbgcat("lambda = ", format(lambda, digits = 6L))

  # ---------------------------------------------------------------------------
  # Coefficient name
  # ---------------------------------------------------------------------------
  # Encode lambda compactly into the coefficient label.
  lambda_tag <- gsub("[^0-9\\.eE\\-]+", "_", format(lambda, digits = 4L))
  base_label <- sprintf("dyadcov_GW[%s]_lambda%s", dyad_label, lambda_tag)
  coef.name  <- base_label
  dbgcat("coef.name = ", coef.name)

  # ---------------------------------------------------------------------------
  # Build INPUT_PARAM for the C side
  # ---------------------------------------------------------------------------
  # Layout:
  #   inputs = c(
  #     as.double(n1),
  #     as.double(lambda),
  #     as.double(Z[1]), ..., as.double(Z[n1*n1])
  #   )
  # where as.double(matrix) uses column-major order, consistent with C indexing.
  inputs <- c(
    as.double(n1),
    as.double(lambda),
    as.double(dyad_mat)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " lambda=", format(lambda, digits = 6L),
         " | Z[1:6]=", paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L),
                             collapse = ","))

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  list(
    name         = "dyadcov_GW",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, lambda, then Z[n1*n1]
    dependence   = TRUE,
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}
