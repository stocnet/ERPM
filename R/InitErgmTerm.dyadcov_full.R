# ==============================================================================
# File    : R/InitErgmTerm.dyadcov_full.R
# Term    : dyadcov_full
# Project : ERPM / ERGM extensions
# ==============================================================================
# Statistic:
#   T = sum_g 1[n_g in S] * sum_{i<j, i,j in g} z_{ij}
#
# INPUT_PARAM layout (C side):
#   c(n1, L, sizes[L], Z[n1*n1])  with Z in column-major order
#
# Debugging:
#   options(ERPM.dyadcov_full.debug = TRUE) to enable debug logs
# ==============================================================================

#' ERGM term: dyadcov_full (within-group dyadic covariate sums)
#'
#' @file InitErgmTerm.dyadcov_full.R
#'
#' @description
#' \code{dyadcov_full} is an ERGM term for bipartite networks that aggregates a
#' symmetric dyadic covariate \eqn{Z = (z_{ij})} over pairs of actors within
#' each group. The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each group in the group mode, we consider all unordered pairs of actors
#' attached to that group and sum the dyadic covariate values \eqn{z_{ij}}. An
#' optional size filter restricts the sum to groups whose size belongs to a
#' specified set \eqn{S}.
#'
#' @details
#' Let:
#' \itemize{
#'   \item \eqn{A} be the set of actor-mode nodes;
#'   \item \eqn{G} be the set of group-mode nodes;
#'   \item \eqn{y} be the bipartite adjacency between actors and groups;
#'   \item \eqn{n_g} be the number of actors attached to group \eqn{g \in G};
#'   \item \eqn{Z = (z_{ij})_{i,j \in A}} be a symmetric numeric matrix.
#' }
#'
#' The statistic implemented by \code{dyadcov_full} is:
#'
#' @code{
#'   T(p; Z, S) =
#'     sum_g 1[n_g in S] * sum_{i<j, i,j in g} z_{ij},
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{S} is a set of admissible group sizes (if \code{size = NULL},
#'         \eqn{S} is implicitly "all positive sizes");
#'   \item the inner sum runs over all unordered pairs of distinct actors
#'         \eqn{i, j} that share group \eqn{g}.
#' }
#'
#' The dyadic covariate \eqn{Z} is provided as an actor-by-actor matrix of size
#' at least \code{n1 x n1}, where \code{n1 = nw \%n\% "bipartite"} is the size
#' of the actor mode. The initializer:
#' \itemize{
#'   \item enforces that the network is bipartite and retrieves \code{n1};
#'   \item accepts \code{dyadcov} either as a literal matrix or as the name of a
#'         network-level attribute containing the matrix;
#'   \item truncates the matrix to its top-left \code{n1 x n1} block if it is
#'         larger than required;
#'   \item checks that the matrix is numeric, free of \code{NA}, and symmetric
#'         up to a small tolerance;
#'   \item parses \code{size} into a sorted vector of distinct positive integers,
#'         or uses the convention "all sizes" when \code{size = NULL}.
#' }
#'
#' @section Implementation and INPUT_PARAM:
#' On the C side, the term is implemented by the change-statistic
#' \code{c_dyadcov_full}. The R initializer below packages the following into
#' the numeric \code{INPUT_PARAM} vector:
#'
#' @code{
#'   INPUT_PARAM = c(
#'     n1,             # actor-mode size
#'     L,              # number of targeted group sizes (length of S)
#'     sizes[1:L],     # sorted distinct positive group sizes
#'     as.double(Z)    # flattened n1 x n1 matrix in column-major order
#'   )
#' }
#'
#' For each toggle of an actor–group edge, the C change-statistic recomputes
#' the within-group dyadic covariate sum for the unique group touched by the
#' toggle, respecting the size filter \eqn{S} when present.
#'
#' @section Arguments:
#' The initializer is not called directly by users. It is invoked automatically
#' by {ergm} when the term \code{dyadcov_full(...)} appears on the right-hand
#' side of a model formula. The user-facing arguments are:
#'
#' @param dyadcov matrix or character. Either:
#'   \itemize{
#'     \item a numeric symmetric matrix of size at least \code{n1 x n1}, where
#'           \code{n1 = nw \%n\% "bipartite"} is the actor-mode size; or
#'     \item the name of a network-level attribute containing such a matrix
#'           (retrieved via \code{nw \%n\% dyadcov}).
#'   }
#'   In both cases the matrix is truncated, if needed, to its top-left
#'   \code{n1 x n1} block.
#' @param size numeric or integer vector of positive group sizes to retain.
#'   If \code{NULL} (the default), all positive group sizes are included
#'   (\eqn{S} = all sizes). If non-\code{NULL}, the vector is rounded to
#'   integers, filtered to positive finite values, and deduplicated. An empty
#'   non-\code{NULL} vector is rejected.
#'
#' @return
#' A standard {ergm} term initialization list with components:
#' \itemize{
#'   \item \code{name}         = \code{"dyadcov_full"};
#'   \item \code{coef.names}   = a single coefficient name encoding the dyadic
#'         covariate label and the size filter;
#'   \item \code{inputs}       = numeric vector
#'         \code{c(n1, L, sizes[L], as.double(Z))};
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
#'   \item The dyadic matrix must be numeric, symmetric (up to a small
#'         tolerance) and free of \code{NA} values. Any \code{NA} triggers a
#'         fail-fast error.
#'   \item When \code{size = NULL}, no size filter is applied and all groups
#'         contribute. When \code{size} is non-\code{NULL}, only groups with
#'         size in the set \eqn{S} contribute.
#'   \item Debug logging is controlled by the global option
#'         \code{options(ERPM.dyadcov_full.debug = TRUE/FALSE)}. When
#'         \code{TRUE}, the initializer prints diagnostic messages prefixed by
#'         \code{"[dyadcov_full][DEBUG]"}.
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
#'   # Group 5 has actors 1, 2, 3  (size 3)
#'   adj[1, 5] <- adj[5, 1] <- 1
#'   adj[2, 5] <- adj[5, 2] <- 1
#'   adj[3, 5] <- adj[5, 3] <- 1
#'   # Group 6 has actors 3, 4    (size 2)
#'   adj[3, 6] <- adj[6, 3] <- 1
#'   adj[4, 6] <- adj[6, 4] <- 1
#'
#'   nw <- network(adj, directed = FALSE, matrix.type = "adjacency")
#'   nw %n% "bipartite" <- n_actors  # actor-mode size
#'
#'   # -------------------------------------------------------------------------
#'   # Define a symmetric dyadic covariate on actors
#'   # -------------------------------------------------------------------------
#'   Z <- matrix(
#'     c(
#'       0, 1, 2, 3,
#'       1, 0, 4, 5,
#'       2, 4, 0, 6,
#'       3, 5, 6, 0
#'     ),
#'     n_actors, n_actors, byrow = TRUE
#'   )
#'   nw %n% "Z_example" <- Z
#'
#'   # -------------------------------------------------------------------------
#'   # Example 1: all group sizes (size = NULL)
#'   # -------------------------------------------------------------------------
#'   summary(
#'     nw ~ dyadcov_full("Z_example"),
#'     constraints = ~ b1part
#'   )
#'
#'   # -------------------------------------------------------------------------
#'   # Example 2: restrict to groups of size 3 only
#'   # -------------------------------------------------------------------------
#'   summary(
#'     nw ~ dyadcov_full("Z_example", size = 3),
#'     constraints = ~ b1part
#'   )
#'
#'   # Fit a simple ERGM with dyadcov_full and a size filter
#'   fit <- ergm(
#'     nw ~ dyadcov_full("Z_example", size = c(2, 3)),
#'     constraints = ~ b1part
#'   )
#'   summary(fit)
#'
#'   # Example call through the ERPM wrapper (network already bipartite)
#'   # erpm(nw ~ dyadcov_full("Z_example", size = c(2, 3)))
#' }
#'
#' @test
#' Self-tests for \code{dyadcov_full} (not shown here) typically:
#' \itemize{
#'   \item construct small bipartite networks with a known partition of actors
#'         into groups and specified size sets \eqn{S};
#'   \item define simple symmetric dyadic matrices \code{Z} for which
#'         \code{sum_g 1[n_g in S] * sum_{i<j, i,j in g} z_{ij}} can be computed
#'         analytically in R;
#'   \item compare \code{summary(nw ~ dyadcov_full(...), constraints = ~ b1part)}
#'         to a direct implementation of the formula above;
#'   \item verify that toggling a single actor–group edge changes the statistic
#'         exactly by the difference between the "before" and "after" within-group
#'         dyadic sums for the affected group, in agreement with the C
#'         change-statistic \code{c_dyadcov_full}.
#' }
#'
#' @keywords ERGM term bipartite dyadic covariate group-sum
#' @md
InitErgmTerm.dyadcov_full <- function(nw, arglist, ...) {
  termname <- "dyadcov_full"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.dyadcov_full.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.dyadcov_full.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[dyadcov_full][DEBUG]", ..., "\n", sep = "")

  # ---------------------------------------------------------------------------
  # Base ERGM term validation and argument parsing
  # ---------------------------------------------------------------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("dyadcov",                         "size"),
    vartypes      = c("matrix,character",               "numeric,integer"),
    defaultvalues = list(NULL,                           NULL),
    required      = c(TRUE,                              FALSE)
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
  # Optional symmetry check on the dyadic matrix
  # ---------------------------------------------------------------------------
  # Small numerical tolerance for symmetric check
  tol <- 1e-8
  if (max(abs(dyad_mat - t(dyad_mat))) > tol) {
    stop(termname, ": la matrice dyadique doit être symétrique (différence > tolérance).")
  }

  dbgcat("dyadcov dim = ", paste(dim(dyad_mat), collapse = "x"),
         " | sample = ",
         paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L), collapse = ","))

  # ---------------------------------------------------------------------------
  # Size filter S: optional set of admissible group sizes
  # ---------------------------------------------------------------------------
  sizes <- a$size
  if (is.null(sizes)) {
    # No filter: all group sizes are allowed (L = 0 encodes "no explicit filter").
    L         <- 0L
    sizes_vec <- numeric(0)
    dbgcat("size filter = <ALL> (L=0)")
    size_label <- "_all"
  } else {
    if (!is.numeric(sizes))
      stop(termname, ": 'size' doit être numérique ou entier.")
    if (length(sizes) == 0L)
      stop(termname, ": 'size' vide (integer(0)) interdit. Utilisez NULL pour toutes tailles.")

    sizes <- as.integer(round(sizes))
    if (any(!is.finite(sizes)) || any(sizes <= 0L))
      stop(termname, ": 'size' doit contenir des entiers positifs.")

    sizes     <- sort(unique(sizes))
    L         <- length(sizes)
    sizes_vec <- as.double(sizes)
    dbgcat("size filter = {", paste(sizes, collapse = ","), "} (L=", L, ")")
    size_label <- sprintf("_S{%s}", paste(sizes, collapse = ","))
  }

  # ---------------------------------------------------------------------------
  # Coefficient name
  # ---------------------------------------------------------------------------
  coef.name <- sprintf("dyadcov_full[%s]%s", dyad_label, size_label)
  dbgcat("coef.name = ", coef.name)

  # ---------------------------------------------------------------------------
  # Build INPUT_PARAM for the C side
  # ---------------------------------------------------------------------------
  # Layout:
  #   inputs = c(
  #     as.double(n1),
  #     as.double(L),
  #     sizes_vec[1:L],
  #     as.double(Z[1]), ..., as.double(Z[n1*n1])
  #   )
  # where as.double(matrix) uses column-major order, consistent with C indexing.
  inputs <- c(
    as.double(n1),
    as.double(L),
    sizes_vec,
    as.double(dyad_mat)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L,
         " | Z[1:6]=", paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L),
                             collapse = ","))

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  list(
    name         = "dyadcov_full",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, L, sizes[L], Z[n1*n1]
    dependence   = TRUE,
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}