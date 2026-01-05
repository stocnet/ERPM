# ==============================================================================
# File    : R/InitErgmTerm.dyadcov_full.R
# Term    : dyadcov_full
# Project : ERPM / ERGM extensions
# ==============================================================================
# Statistic:
#   T = sum_g 1[n_g in S] * sum_{i != j, i,j in g} z_{ij}
#
# INPUT_PARAM layout (C side):
#   c(n1, L, sizes[L], Z[n1*n1])  with Z in column-major order
#
# Debugging:
#   options(ERPM.dyadcov_full.debug = TRUE) to enable debug logs
# ==============================================================================

#' ERGM term: dyadcov_full (within-group dyadic covariate sums)
#'
#' @name InitErgmTerm.dyadcov_full
#' @aliases dyadcov_full
#' @note InitErgmTerm.dyadcov_full.R
#'
#' @description
#' \code{dyadcov_full} is an ERGM term for bipartite networks that aggregates a
#' numeric dyadic covariate \eqn{Z = (z_{ij})} over *ordered* pairs of actors that
#' are adjacent to the same group.
#'
#' The network is interpreted as:
#' \itemize{
#'   \item an \emph{actor mode} (the side identified by \code{nw \%n\% "bipartite"});
#'   \item a \emph{group mode} (the complementary side of the bipartite graph).
#' }
#'
#' For each group node, we consider the set of adjacent actors and sum \eqn{z_{ij}}
#' over all ordered pairs \eqn{(i,j)} with \eqn{i \ne j} within that group.
#'
#' An optional size filter \code{size} restricts which groups contribute by
#' requiring the group size \eqn{n_g} to be in the provided set.
#'
#' @details
#' The dyadic covariate \eqn{Z} is a real actor-by-actor matrix of size
#' \eqn{n_A \times n_A}, where \eqn{n_A} is the size of the actor mode
#' (\code{n1 = nw \%n\% "bipartite"}).
#'
#' The initializer:
#' \itemize{
#'   \item enforces bipartiteness and retrieves the actor-mode size \eqn{n_A};
#'   \item accepts \code{dyadcov} either as a literal matrix or as the name of a
#'         network-level attribute (\code{nw \%n\% "..."});
#'   \item truncates \code{dyadcov} to its top-left \code{n1 x n1} block if larger;
#'   \item checks that \code{dyadcov} is numeric and free of \code{NA};
#'   \item parses \code{size} into a sorted, distinct set of positive integers.
#' }
#'
#' @section Implementation and change-statistic:
#' The term is implemented as a native ERGM C change-statistic, exposed under the
#' symbol \code{c_dyadcov_full}. The R initializer below:
#' \itemize{
#'   \item packages \code{n1}, \code{L}, \code{sizes[1:L]} and the flattened
#'         \code{dyadcov} matrix into \code{INPUT_PARAM};
#'   \item declares the term as dependent (\code{dependence = TRUE}) with no
#'         finite \code{minval}/\code{maxval};
#'   \item sets the empty-network statistic to \code{0}.
#' }
#'
#' For each toggle of an actor–group edge, the C change-statistic recomputes
#' the within-group dyadic covariate sum for the unique group touched by the
#' toggle, respecting the size filter \eqn{S} when present.
#'
#' @section Arguments:
#' The initializer is not called directly by users; it is invoked automatically
#' by {ergm} when the term \code{dyadcov_full(...)} appears on the right-hand side
#' of a model formula.
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
#' @param size optional numeric/integer vector. If provided, only groups whose
#'   actor-degree is in \code{size} contribute to the statistic.
#'
#' @return
#' A standard {ergm} term initialization list with components:
#' \itemize{
#'   \item \code{name}         = \code{"dyadcov_full"};
#'   \item \code{coef.names}   = a single coefficient name encoding
#'         \code{dyadcov} label and the size filter;
#'   \item \code{inputs}       = numeric vector
#'         \code{c(n1, L, sizes[1:L], as.double(Z))};
#'   \item \code{dependence}   = \code{TRUE};
#'   \item \code{minval}       = \code{-Inf};
#'   \item \code{maxval}       = \code{Inf};
#'   \item \code{emptynwstats} = \code{0}.
#' }
#'
#' @keywords ERGM term bipartite dyadic covariate
#' @md
InitErgmTerm.dyadcov_full <- function(nw, arglist, ...) {
  termname <- "dyadcov_full"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  dbg    <- isTRUE(getOption("ERPM.dyadcov_full.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[dyadcov_full][DEBUG]", ..., "\n", sep = "")

  # ---------------------------------------------------------------------------
  # Base ERGM term validation and argument parsing
  # ---------------------------------------------------------------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("dyadcov",                 "size"),
    vartypes      = c("matrix,character",        "numeric,integer"),
    defaultvalues = list(NULL,                   NULL),
    required      = c(TRUE,                      FALSE)
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
    # Case: name of a dyadic matrix.
    #
    # Since 2026-01, the ERPM bipartite builder stores all dyadic matrices under
    # a single network attribute `%n% "dyads"` (a named list). For backward
    # compatibility, we still support the historical behavior where the matrix
    # is stored directly as a top-level `%n%` attribute.
    dyads_list <- nw %n% "dyads"

    if (is.list(dyads_list) && !is.null(dyads_list[[dyad_raw]])) {
      dyad_mat   <- dyads_list[[dyad_raw]]
      dyad_label <- paste0("dyads$", dyad_raw)
      dbgcat("dyadcov source = nw %n% 'dyads'$", sQuote(dyad_raw))
    } else {
      dyad_mat   <- nw %n% dyad_raw
      dyad_label <- dyad_raw
      if (is.null(dyad_mat))
        stop(termname, ": matrice dyadique introuvable: ",
             sQuote(dyad_raw),
             " (cherché dans nw %n% 'dyads' puis dans nw %n% ", sQuote(dyad_raw), ").")
      dbgcat("dyadcov source = network attribute ", sQuote(dyad_label))
    }
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

  dbgcat("dyadcov dim = ", paste(dim(dyad_mat), collapse = "x"),
         " | sample = ",
         paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L), collapse = ","))

  # ---------------------------------------------------------------------------
  # Parse size filter (S)
  # ---------------------------------------------------------------------------
  sizes_raw <- a$size
  if (is.null(sizes_raw) || length(sizes_raw) == 0L) {
    sizes_vec <- numeric(0L)
    L <- 0L
    size_label <- ""
  } else {
    if (!is.numeric(sizes_raw))
      stop(termname, ": 'size' doit être numérique ou entier.")
    iv <- as.integer(round(sizes_raw))
    if (any(!is.finite(sizes_raw)) || any(iv <= 0L) || !isTRUE(all.equal(sizes_raw, iv))) {
      stop(termname, ": 'size' doit contenir des entiers positifs.")
    }
    iv <- sort(unique(iv))
    sizes_vec <- as.double(iv)
    L <- length(iv)
    size_label <- paste0("_size", paste(iv, collapse = "_"))
  }

  dbgcat("size filter: L=", L,
         if (L) paste0(" | sizes=", paste(sizes_vec, collapse = ",")) else " | none")

  # ---------------------------------------------------------------------------
  # Coefficient name
  # ---------------------------------------------------------------------------
  coef.name <- sprintf("dyadcov_full[%s]%s", dyad_label, size_label)
  dbgcat("coef.name = ", coef.name)

  # ---------------------------------------------------------------------------
  # Build INPUT_PARAM for the C side
  # ---------------------------------------------------------------------------
  inputs <- c(
    as.double(n1),
    as.double(L),
    sizes_vec,
    as.double(dyad_mat)
  )

  dbgcat("inputs summary: len=", length(inputs),
         " | n1=", n1, " L=", L,
         " | Z[1:6]=",
         paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L), collapse = ","))

  # ---------------------------------------------------------------------------
  # Standard ERGM term initialization return value
  # ---------------------------------------------------------------------------
  list(
    name         = "dyadcov_full",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, L, sizes[1:L], then Z[n1*n1]
    dependence   = TRUE,
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}
