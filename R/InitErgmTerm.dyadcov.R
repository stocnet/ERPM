# ==============================================================================
# File    : R/InitErgmTerm.dyadcov.R
# Term    : dyadcov
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM term: dyadcov (clique-based dyadic covariate)
#'
#' @file InitErgmTerm.dyadcov.R
#'
#' @description
#' \code{dyadcov} is an ERGM term for bipartite networks that aggregates a
#' symmetric dyadic covariate \eqn{Z = (z_{ij})} over cliques of actors within
#' each group. The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each group node in the group mode, we consider the set of adjacent actors
#' and all \eqn{k}-cliques of actors within that group. For a given clique
#' \eqn{C} of size \eqn{k}, the contribution is the product of dyadic covariates
#' over all actor pairs in the clique:
#' \deqn{
#'   \prod_{i<j,\, i,j \in C} z_{ij}.
#' }
#'
#' Two variants are supported:
#' \itemize{
#'   \item \code{normalized = FALSE} (raw sum over cliques);
#'   \item \code{normalized = TRUE} (average over all actor subsets of size
#'         \eqn{k} within each group).
#' }
#'
#' @details
#' The dyadic covariate \eqn{Z} is a real, symmetric, actor-by-actor matrix
#' of size \eqn{n_A \times n_A}, where \eqn{n_A} is the size of the actor mode
#' (\code{n1 = nw \%n\% "bipartite"}). The initializer:
#' \itemize{
#'   \item enforces bipartiteness and retrieves the actor-mode size \eqn{n_A};
#'   \item accepts \code{dyadcov} either as a literal matrix or as the name of a
#'         network-level attribute (\code{nw \%n\% "..."});
#'   \item truncates \code{dyadcov} to its top-left \code{n1 x n1} block if larger;
#'   \item checks that \code{dyadcov} is numeric, free of \code{NA}, and symmetric
#'         up to a small tolerance;
#'   \item parses \code{clique_size} into an integer \eqn{k \ge 2};
#'   \item parses \code{normalized} into a logical flag.
#' }
#'
#' The core statistic (for a fixed \eqn{k} and scalar \code{normalized}) can be
#' written in terms of the partition \eqn{p} and the dyadic covariate \eqn{Z}:
#'
#' @code
#'   - normalized = FALSE :
#'       T^{(k)}(p; Z)
#'         = sum_g sum_{C in C_k(g)} prod_{i<j in C} z_{ij}
#'
#'   - normalized = TRUE :
#'       T^{(k)}_norm(p; Z)
#'         = sum_g 1[n_g >= k] * (1 / choose(n_g, k)) *
#'             sum_{C in C_k(g)} prod_{i<j in C} z_{ij}
#'
#' where:
#' \itemize{
#'   \item \eqn{g} runs over group-mode nodes (groups);
#'   \item \eqn{n_g} is the number of actors adjacent to group \eqn{g};
#'   \item \eqn{C_k(g)} is the set of all subsets \eqn{C \subseteq A(g)} of size
#'         \eqn{k}, where \eqn{A(g)} is the actor set attached to group \eqn{g}.
#' }
#'
#' @section Implementation and change-statistic:
#' The term is implemented as a native ERGM C change-statistic, exposed under the
#' symbol \code{c_dyadcov}. The R initializer below:
#' \itemize{
#'   \item packages \code{n1}, \code{k}, \code{normalized} and the flattened
#'         \code{dyadcov} matrix into \code{INPUT_PARAM};
#'   \item declares the term as dependent (\code{dependence = TRUE}) with no
#'         finite \code{minval}/\code{maxval};
#'   \item sets the empty-network statistic to \code{0}.
#' }
#'
#' On each toggle of an actor–group edge, the C change-statistic recomputes the
#' local contribution for the unique group affected by the toggle, by updating
#' the clique-based sums according to the new set of actors in that group.
#'
#' @section Arguments:
#' The initializer is not called directly by users; it is invoked automatically
#' by {ergm} when the term \code{dyadcov(...)} appears on the right-hand side of
#' a model formula.
#'
#' @param dyadcov matrix or character. Either:
#'   \itemize{
#'     \item a numeric symmetric matrix of size at least \code{n1 x n1}, where
#'           \code{n1 = nw \%n\% "bipartite"} is the actor-mode size; or
#'     \item the name of a network-level attribute containing such a matrix
#'           (retrieved as \code{nw \%n\% dyadcov}).
#'   }
#'   In both cases the matrix is truncated, if necessary, to its top-left
#'   \code{n1 x n1} block.
#' @param clique_size numeric or integer scalar. Target clique size \eqn{k}.
#'   It is rounded and coerced to an integer and must satisfy \eqn{k \ge 2}.
#'   Defaults to \code{2}.
#' @param normalized logical scalar. When \code{FALSE}, the statistic is a raw
#'   sum over all \eqn{k}-cliques per group. When \code{TRUE}, each group's
#'   contribution is normalized by \code{choose(n_g, k)} (and groups with
#'   \eqn{n_g < k} contribute \code{0}). Defaults to \code{FALSE}.
#'
#' @return
#' A standard {ergm} term initialization list with components:
#' \itemize{
#'   \item \code{name}         = \code{"dyadcov"};
#'   \item \code{coef.names}   = a single coefficient name encoding
#'         \code{dyadcov} label, \code{k} and the normalization flag;
#'   \item \code{inputs}       = numeric vector
#'         \code{c(n1, k, normalized_flag, as.double(Z))};
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
#'         be a positive integer.
#'   \item The dyadic matrix must be numeric, symmetric (up to a small tolerance)
#'         and free of \code{NA} values.
#'   \item Debug logging is controlled by the global option
#'         \code{options(ERPM.dyadcov.debug = TRUE/FALSE)}. When \code{TRUE},
#'         the initializer prints diagnostic messages prefixed by
#'         \code{"[dyadcov][DEBUG]"}.
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
#'   # Build a symmetric dyadic covariate on actors
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
#'   # Inspect the dyadcov statistic
#'   # -------------------------------------------------------------------------
#'   # k = 2, raw sum over all cliques of size 2 (i.e. pairs) within groups
#'   summary(nw ~ dyadcov("Z_example", clique_size = 2, normalized = FALSE),
#'           constraints = ~ b1part)
#'
#'   # k = 3, normalized by the number of 3-actor subsets in each group
#'   summary(nw ~ dyadcov("Z_example", clique_size = 3, normalized = TRUE),
#'           constraints = ~ b1part)
#'
#'   # -------------------------------------------------------------------------
#'   # Example with the ERPM wrapper (network already bipartite)
#'   # -------------------------------------------------------------------------
#'   # erpm(nw ~ dyadcov("Z_example", clique_size = 2))
#' }
#'
#' @test
#' Self-tests for \code{dyadcov} (not shown here) typically:
#' \itemize{
#'   \item build small bipartite networks with a known partition of actors into
#'         groups;
#'   \item define simple symmetric dyadic matrices \code{Z} so that clique
#'         products can be computed analytically;
#'   \item compare \code{summary(nw ~ dyadcov(...), constraints = ~ b1part)} with
#'         a direct implementation of
#'         \code{sum_g sum_{C in C_k(g)} prod_{i<j in C} Z[i,j]} (raw or normalized);
#'   \item verify that toggling a single actor–group edge changes the statistic
#'         by the local difference between the "before" and "after" clique sums,
#'         in agreement with the C change-statistic \code{c_dyadcov}.
#' }
#'
#' @keywords ERGM term bipartite dyadic covariate cliques
#' @md
InitErgmTerm.dyadcov <- function(nw, arglist, ...) {
  termname <- "dyadcov"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.dyadcov.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.dyadcov.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[dyadcov][DEBUG]", ..., "\n", sep = "")

  # ---------------------------------------------------------------------------
  # Base ERGM term validation and argument parsing
  # ---------------------------------------------------------------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("dyadcov",                 "clique_size",      "normalized"),
    vartypes      = c("matrix,character",        "numeric,integer",  "logical"),
    defaultvalues = list(NULL,                   2,                  FALSE),
    required      = c(TRUE,                      FALSE,              FALSE)
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
  # Optional symmetry check
  # ---------------------------------------------------------------------------
  tol <- 1e-8
  if (max(abs(dyad_mat - t(dyad_mat))) > tol) {
    stop(termname, ": la matrice dyadique doit être symétrique (différence > tolérance).")
  }

  dbgcat("dyadcov dim = ", paste(dim(dyad_mat), collapse = "x"),
         " | sample = ",
         paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L), collapse = ","))

  # ---------------------------------------------------------------------------
  # clique_size (k)
  # ---------------------------------------------------------------------------
  k_raw <- a$clique_size
  if (is.null(k_raw) || length(k_raw) == 0L) {
    k <- 2L
  } else {
    if (!is.numeric(k_raw))
      stop(termname, ": 'clique_size' doit être numérique ou entier.")
    k <- as.integer(round(k_raw[1L]))
  }

  if (!is.finite(k) || k < 2L)
    stop(termname, ": 'clique_size' doit être un entier >= 2.")

  dbgcat("clique_size (k) = ", k)

  # ---------------------------------------------------------------------------
  # normalized (boolean flag)
  # ---------------------------------------------------------------------------
  norm_raw <- a$normalized
  if (is.null(norm_raw) || length(norm_raw) == 0L) {
    normalized <- FALSE
  } else {
    normalized <- isTRUE(norm_raw[1L])
  }
  norm_flag <- if (normalized) 1 else 0
  dbgcat("normalized = ", normalized)

  # ---------------------------------------------------------------------------
  # Coefficient name
  # ---------------------------------------------------------------------------
  # Encode the dyadic matrix label, the clique size, and the normalization flag.
  base_label <- sprintf("dyadcov[%s]_k%d", dyad_label, k)
  coef.name  <- if (normalized) paste0(base_label, "_norm") else base_label
  dbgcat("coef.name = ", coef.name)

  # ---------------------------------------------------------------------------
  # Build INPUT_PARAM for the C side
  # ---------------------------------------------------------------------------
  # Layout:
  #   inputs = c(
  #     as.double(n1),
  #     as.double(k),
  #     as.double(normalized_flag),
  #     as.double(Z[1]), ..., as.double(Z[n1*n1])
  #   )
  # where as.double(matrix) uses column-major order, consistent with C indexing.
  inputs <- c(
    as.double(n1),
    as.double(k),
    as.double(norm_flag),
    as.double(dyad_mat)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " k=", k, " normalized=", norm_flag,
         " | Z[1:6]=", paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L),
                             collapse = ","))

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  list(
    name         = "dyadcov",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, k, normalized flag, then Z[n1*n1]
    dependence   = TRUE,
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}