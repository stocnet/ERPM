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
#   options(ERPM.dyadcov_full.debug = TRUE) to enable debug logs (R side)
#
# IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#   - dyadcov_full is intended to support multi-toggle proposals (swap/split/merge)
#     decomposed into a list of toggles.
#   - Therefore, the compiled change-statistic MUST be implemented as D_CHANGESTAT_FN,
#     and the R initializer MUST advertise this by returning `d_func = TRUE`.
#   - Symbol naming convention:
#       C side: D_CHANGESTAT_FN(d_dyadcov_full)
#       R side: d_func=TRUE and name="dyadcov_full"
# ==============================================================================

#' ERGM term: dyadcov_full (within-group dyadic covariate sums)
#'
#' @name InitErgmTerm.dyadcov_full
#' @aliases dyadcov_full
#' @note InitErgmTerm.dyadcov_full.R
#' @author Jérémie Chichignoud - Cub'itech
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
#'   \item supports the ERPM convention where dyadic matrices are stored in
#'         \code{nw \%n\% "dyads"} as a named list (e.g. \code{"Z1"}, \code{"Z2"});
#'   \item truncates \code{dyadcov} to its top-left \code{n1 x n1} block if larger;
#'   \item checks that \code{dyadcov} is numeric and free of \code{NA};
#'   \item parses \code{size} into a sorted, distinct set of positive integers.
#' }
#'
#' @section Implementation and change-statistic:
#' The term is implemented as a native ERGM C change-statistic in MULTI-TOGGLE form,
#' exposed under the symbol \code{d_dyadcov_full} via \code{D_CHANGESTAT_FN}.
#'
#' IMPORTANT:
#' \itemize{
#'   \item The initializer MUST return \code{d_func = TRUE}.
#'   \item Otherwise, \pkg{ergm} will assume a one-toggle \code{C_CHANGESTAT_FN}
#'         entrypoint and call the function with the wrong signature.
#' }
#'
#' The R initializer below:
#' \itemize{
#'   \item packages \code{n1}, \code{L}, \code{sizes[1:L]} and the flattened
#'         \code{dyadcov} matrix into \code{INPUT_PARAM};
#'   \item declares the term as dependent (\code{dependence = TRUE});
#'   \item sets the empty-network statistic to \code{0}.
#' }
#'
#' @section Arguments:
#' The initializer is not called directly by users; it is invoked automatically
#' by \pkg{ergm} when the term \code{dyadcov_full(...)} appears on the right-hand side
#' of a model formula.
#'
#' @param nw A \pkg{network} object.
#' @param arglist A named list of term arguments. Expected components:
#'   \itemize{
#'     \item \code{dyadcov}: matrix or character. Either a numeric matrix of size at least
#'           \code{n1 x n1} (with \code{n1 = nw \%n\% "bipartite"}), or the name of a
#'           network-level attribute containing such a matrix (retrieved as
#'           \code{nw \%n\% dyadcov}). Preferably, use ERPM convention:
#'           \code{nw \%n\% "dyads"} as a named list and pass \code{"Z1"} / \code{"Z2"}.
#'     \item \code{size}: optional numeric/integer vector. If provided, only groups whose
#'           actor-degree is in \code{size} contribute to the statistic.
#'   }
#' @param ... Passed through by \pkg{ergm}; not used.
#'
#' @return A standard \pkg{ergm} term initialization list.
#'
#' @keywords ERGM term bipartite dyadic covariate
#' @md
#' @export
InitErgmTerm.dyadcov_full <- function(nw, arglist, ...) {
  termname <- "dyadcov_full"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  dbg    <- isTRUE(getOption("ERPM.dyadcov_full.debug", TRUE))
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
  if (is.na(n1) || n1 <= 0L) {
    stop(termname, ": strictly bipartite network required (attribut %n% 'bipartite' manquant ou invalide).")
  }
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
      if (is.null(dyad_mat)) {
        stop(
          termname, ": dyadic matrix not found: ", sQuote(dyad_raw),
          " (looked up in nw %n% 'dyads', then in nw %n% ", sQuote(dyad_raw), ")."
        )
      }
      dbgcat("dyadcov source = network attribute ", sQuote(dyad_label))
    }
  } else {
    # Case: matrix passed literally as an argument
    dyad_mat   <- dyad_raw
    dyad_label <- "dyadcov"
    dbgcat("dyadcov source = literal matrix")
  }

  if (!is.matrix(dyad_mat)) {
    stop(termname, ": 'dyadcov' must be a matrix or the name of a network-level attribute.")
  }

  nr <- nrow(dyad_mat)
  nc <- ncol(dyad_mat)

  if (nr < n1 || nc < n1) {
    stop(termname, ": dyadic matrix dimension (", nr, "x", nc, ") insuffisantes pour n1 = ", n1, ".")
  }

  # If larger than needed, restrict to the top-left n1 x n1 block
  if (nr > n1 || nc > n1) {
    dyad_mat <- dyad_mat[seq_len(n1), seq_len(n1), drop = FALSE]
    dbgcat("dyadcov truncated to ", n1, "x", n1)
  }

  # ---------------------------------------------------------------------------
  # Numeric coercion and fail-fast on NA
  # ---------------------------------------------------------------------------
  if (!is.numeric(dyad_mat)) stop(termname, ": dyadic matrix must be numeric.")
  if (anyNA(dyad_mat))       stop(termname, ": NA values are not allowed in the dyadic matrix.")

  dbgcat("dyadcov dim = ", paste(dim(dyad_mat), collapse = "x"),
         " | sample = ",
         paste(utils::head(signif(as.numeric(dyad_mat), 5L), 6L), collapse = ","))

  # ---------------------------------------------------------------------------
  # Parse size filter (S)
  # ---------------------------------------------------------------------------
  sizes_raw <- a$size
  if (is.null(sizes_raw) || length(sizes_raw) == 0L) {
    sizes_vec  <- numeric(0L)
    L          <- 0L
    size_label <- ""
  } else {
    if (!is.numeric(sizes_raw)) stop(termname, ": 'size' must be numeric ou entier.")
    iv <- as.integer(round(sizes_raw))
    if (any(!is.finite(sizes_raw)) || any(iv <= 0L) || !isTRUE(all.equal(sizes_raw, iv))) {
      stop(termname, ": 'size' must contain positive integers.")
    }
    iv <- sort(unique(iv))
    sizes_vec  <- as.double(iv)
    L          <- length(iv)
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
  # IMPORTANT:
  # - `d_func = TRUE` tells ergm to call the multi-toggle (D_) changestat entrypoint.
  # - Without it, ergm assumes a one-toggle C_ changestat and will call the function
  #   with the wrong signature if you compiled only a D_ function (=> crash).
  list(
    name         = "dyadcov_full",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, L, sizes[1:L], then Z[n1*n1]
    dependence   = TRUE,
    d_func       = TRUE,        # <-- REQUIRED (multi-toggle D_CHANGESTAT_FN)
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}