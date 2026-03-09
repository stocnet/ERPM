# ==============================================================================
# File    : R/InitErgmTerm.dyadcov_GW.R
# Term    : dyadcov_GW
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM term: dyadcov_GW (geometrically weighted dyadic covariate over cliques)
#'
#' @name InitErgmTerm.dyadcov_GW
#' @aliases dyadcov_GW
#' @note InitErgmTerm.dyadcov_GW.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' \code{dyadcov_GW} is an ERGM term for bipartite networks that aggregates a
#' dyadic covariate matrix \eqn{Z=(z_{ij})} over actor cliques within each group,
#' using geometric weights controlled by \code{lambda}.
#'
#' The dyadic covariate \eqn{Z} is a real actor-by-actor matrix of
#' size \eqn{n_A \times n_A}, where \eqn{n_A} is the size of the actor mode
#' \code{n1 = nw \%n\% "bipartite"}. Symmetry is not required: both orientations
#' \eqn{z_{ij}} and \eqn{z_{ji}} are used through the symmetrised form
#' \eqn{z_{ij}+z_{ji}}.
#'
#' @section Implementation and change-statistic:
#' The term is implemented as a native ERGM C change-statistic.
#'
#' IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#' - This term MUST support multi-toggle proposals (swap/split/merge decomposed
#'   into multiple edge toggles).
#' - Therefore the compiled change-statistic is implemented as a D_ entrypoint:
#'     D_CHANGESTAT_FN(d_dyadcov_GW)
#' - On the R side, we MUST set `d_func = TRUE`, otherwise ergm will try to call
#'   a one-toggle C_ entrypoint and you risk a signature mismatch / segfault.
#'
#' Parameter packing:
#'   INPUT_PARAM = c(n1, lambda, as.vector(Z))
#' where as.vector(Z) uses R's column-major order (consistent with C indexing).
#'
#' @param nw A \pkg{network} object.
#' @param arglist A named list of term arguments. Expected components include
#'   \code{dyadcov} (matrix or character) and \code{lambda} (numeric scalar).
#' @param ... Passed through by \pkg{ergm}; not used.
#'
#' @details
#' Term arguments are passed via \code{arglist} by \pkg{ergm}:
#' \itemize{
#' \item \code{dyadcov}: matrix or character. Either a numeric matrix, or the name of a
#' network-level attribute containing such a matrix (retrieved as \code{nw \%n\% dyadcov}).
#' In ERPM usage, the matrix can also come from \code{nw \%n\% "dyads"} (a named list),
#' i.e. \code{(nw \%n\% "dyads")[[dyadcov]]}.
#' \item \code{lambda}: numeric scalar. Geometric weight parameter \eqn{\lambda > 0}.
#' }
#'
#' @keywords ERGM term bipartite dyadic covariate geometrically-weighted cliques
#' @md
#' @export
InitErgmTerm.dyadcov_GW <- function(nw, arglist, ...) {
  termname <- "dyadcov_GW"

  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  # Global option:
  #   options(ERPM.dyadcov_GW.debug = TRUE/FALSE)
  # When TRUE, the initializer prints diagnostic messages to the console.
  dbg    <- isTRUE(getOption("ERPM.dyadcov_GW.debug", TRUE))
  dbgcat <- function(...) if (dbg) cat("[dyadcov_GW][DEBUG]", ..., "\n", sep = "")

  # ---------------------------------------------------------------------------
  # Base ERGM term validation and argument parsing
  # ---------------------------------------------------------------------------
  a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("dyadcov",         "lambda"),
    vartypes      = c("matrix,character","numeric"),
    defaultvalues = list(NULL,            2),
    required      = c(TRUE,              FALSE)
  )

  # ---------------------------------------------------------------------------
  # Actor-mode size (n1) from the bipartite attribute
  # ---------------------------------------------------------------------------
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L)
    stop(termname, ": strictly bipartite network required (attribut %n% 'bipartite' manquant ou invalide).")
  dbgcat("n1 = ", n1)

  # ---------------------------------------------------------------------------
  # Retrieve the dyadic matrix (n1 x n1)
  # ---------------------------------------------------------------------------
  dyad_raw   <- a$dyadcov
  dyad_label <- NULL
  dyad_mat   <- NULL

  if (is.character(dyad_raw) && length(dyad_raw) == 1L) {
    # Case A: name of a network-level attribute storing the matrix (ERGM standard)
    dyad_mat <- nw %n% dyad_raw
    if (!is.null(dyad_mat)) {
      dyad_label <- dyad_raw
      dbgcat("dyadcov source = network attribute ", sQuote(dyad_label))
    } else {
      # Case B (ERPM): look for a named list of dyadic matrices in nw %n% "dyads"
      dyads_list <- nw %n% "dyads"
      if (is.list(dyads_list) && !is.null(dyads_list[[dyad_raw]])) {
        dyad_mat   <- dyads_list[[dyad_raw]]
        dyad_label <- paste0("dyads$", dyad_raw)
        dbgcat("dyadcov source = dyads list, key=", sQuote(dyad_raw))
      } else {
        stop(termname, ": nonexistent network-level attribute and dyads[[",
             sQuote(dyad_raw), "]] not found")
      }
    }
  } else {
    # Case C: matrix passed literally as an argument
    dyad_mat   <- dyad_raw
    dyad_label <- "dyadcov"
    dbgcat("dyadcov source = matrix literal")
  }

  if (!is.matrix(dyad_mat))
    stop(termname, ": 'dyadcov' must be a matrix or the name of a network-level attribute.")

  nr <- nrow(dyad_mat)
  nc <- ncol(dyad_mat)
  if (nr < n1 || nc < n1)
    stop(termname, ": dimensions de la matrice dyadique (", nr, "x", nc,
         ") insuffisantes pour n1 = ", n1, ".")

  # Truncate to actor-mode block
  if (nr != n1 || nc != n1) {
    dyad_mat <- dyad_mat[seq_len(n1), seq_len(n1), drop = FALSE]
    dbgcat("dyadcov truncated to ", n1, "x", n1)
  }

  # Numeric coercion + NA guard
  if (is.logical(dyad_mat) || is.integer(dyad_mat)) dyad_mat <- as.numeric(dyad_mat)
  if (!is.numeric(dyad_mat))
    stop(termname, ": the dyadic matrix must be numeric (or coercible to numeric).")
  if (anyNA(dyad_mat))
    stop(termname, ": NA values are not allowed in the dyadic matrix.")

  # Optional symmetry diagnostics (no requirement)
  if (dbg) {
    asym <- max(abs(dyad_mat - t(dyad_mat)))
    dbgcat("max |Z - t(Z)| = ", format(asym, digits = 6L))
  }

  # ---------------------------------------------------------------------------
  # Lambda handling
  # ---------------------------------------------------------------------------
  lambda <- as.double(a$lambda)
  if (!length(lambda) || is.na(lambda)) lambda <- 2
  if (length(lambda) != 1L)
    stop(termname, ": 'lambda' must be a scalar.")
  if (!is.finite(lambda) || lambda <= 0)
    stop(termname, ": 'lambda' must be a strictly positive real number (and typically > 1).")
  dbgcat("lambda = ", format(lambda, digits = 6L))

  # ---------------------------------------------------------------------------
  # Coefficient name
  # ---------------------------------------------------------------------------
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
  # IMPORTANT:
  # - `d_func = TRUE` tells ergm to call the multi-toggle (D_) changestat entrypoint.
  # - Without it, ergm assumes a one-toggle C_ changestat and will call the function
  #   with the wrong signature if you compiled only a D_ function (=> segfault).
  list(
    name         = "dyadcov_GW",
    coef.names   = coef.name,
    inputs       = inputs,      # n1, lambda, then Z[n1*n1]
    dependence   = TRUE,
    d_func       = TRUE,        # <-- REQUIRED for D_CHANGESTAT_FN (multi-toggle)
    minval       = -Inf,
    maxval       = Inf,
    emptynwstats = 0
  )
}