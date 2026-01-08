################################################################################
# FILE: R/erpm_long_helpers.R
################################################################################
#' ERPM longitudinal helpers: small utilities for erpm_long() (internal)
#'
#' @name erpm_long_helpers
#' @note erpm_long_helpers.R
#'
#' @description
#' This file defines a collection of small, internal helper functions used by
#' \code{erpm_long()} and its execution engine.
#'
#' These helpers are intentionally:
#' \itemize{
#'   \item stateless and side-effect free;
#'   \item focused on one narrow task (type checks, formatting, indexing);
#'   \item reusable across the longitudinal pipeline.
#' }
#'
#' They are separated from the main logic of \code{erpm_long()} to:
#' \itemize{
#'   \item improve readability of the main algorithm;
#'   \item reduce cognitive load when debugging or extending the code;
#'   \item make implicit assumptions explicit (e.g. expected input shapes).
#' }
#'
#' All functions are internal (non-exported) and prefixed with \code{.erpm_long_}.
#'
#' @keywords ERPM ERGM longitudinal internal helpers
NULL

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Check whether an object is a list of data.frames
#'
#' This is used to validate inputs such as `nodes` that may be provided either
#' as a single data.frame (recycled over time) or as a list of per-time data.frames.
#'
#' @param x Any R object.
#' @return Logical scalar. TRUE if `x` is a non-empty list and all elements
#'         are data.frames.
#' @noRd
.erpm_long_is_list_of_df <- function(x) {
  is.list(x) &&
    length(x) &&
    all(vapply(x, is.data.frame, logical(1)))
}

#' Check whether an object is a list of dyadic inputs
#'
#' In erpm_long(), dyads are expected to be supplied as:
#' \itemize{
#'   \item a single list of matrices (recycled over time), or
#'   \item a list of such lists, one per time point.
#' }
#'
#' This helper only checks the *outer* structure: list of lists.
#'
#' @param x Any R object.
#' @return Logical scalar. TRUE if `x` is a non-empty list whose elements
#'         are themselves lists.
#' @noRd
.erpm_long_is_list_of_dyads <- function(x) {
  is.list(x) &&
    length(x) &&
    all(vapply(x, is.list, logical(1)))
}

#' Extract the time-\eqn{t} element of a possibly time-indexed input
#'
#' Many inputs in erpm_long() (nodes, dyads, group_labels) can be provided as:
#' \itemize{
#'   \item NULL,
#'   \item a single object (used at all time points),
#'   \item a list of length T, providing one object per time point.
#' }
#'
#' This helper centralizes the logic for safely extracting the correct
#' element at time `t`, and for emitting a clear error if the length is invalid.
#'
#' @param x Input object (NULL, single object, or list).
#' @param t Integer time index (1-based).
#' @param T Total number of time points.
#' @param what Character string used in error messages to identify the input.
#' @return The object to be used at time `t`, or NULL.
#' @noRd
.erpm_long_get_t <- function(x, t, T, what = "input") {
  # NULL stays NULL: no input provided.
  if (is.null(x)) return(NULL)

  # Non-list or data.frame inputs are treated as time-invariant
  # and reused at every t.
  if (!is.list(x) || is.data.frame(x)) return(x)

  # If a list is provided, it must match the number of time points.
  if (length(x) != T) {
    stop(sprintf(
      "[ERPM_LONG] `%s` must be NULL, a single object, or a list of length T=%d.",
      what, T
    ), call. = FALSE)
  }

  # Extract the element for time t.
  x[[t]]
}

#' Convert a RHS expression into a single-line character string
#'
#' Used purely for verbose/debug output, to make logs compact and readable.
#' Long expressions are deparsed with a large width cutoff and collapsed
#' onto a single line.
#'
#' @param expr A language object (typically a RHS expression).
#' @return Character string representation on one line.
#' @noRd
.erpm_long_rhs_oneline <- function(expr) {
  paste(deparse(expr, width.cutoff = 500L), collapse = " ")
}

#' Combine a list of ERGM terms into a single RHS expression
#'
#' Given a list of calls/symbols representing ERGM terms, this helper
#' reconstructs a sum-of-terms expression using `+`.
#'
#' Special cases:
#' \itemize{
#'   \item empty list   -> NULL (no RHS);
#'   \item length 1     -> the term itself (no wrapping);
#'   \item length > 1   -> Reduce with `+`.
#' }
#'
#' @param terms List of language objects (calls or symbols).
#' @return A language object suitable for use on the RHS of a formula, or NULL.
#' @noRd
.erpm_long_combine_terms <- function(terms) {
  if (!length(terms)) return(NULL)
  if (length(terms) == 1L) return(terms[[1L]])

  Reduce(function(x, y) call("+", x, y), terms)
}

#' Compact formatter for named integer tables
#'
#' This helper is used only for verbose output to summarize group size tables
#' in a compact, stable textual form (e.g. "1:3 2:5 3:2").
#'
#' @param x A named integer vector or table.
#' @return Character string.
#' @noRd
.erpm_long_tab_str <- function(x) {
  paste(sprintf("%s:%s", names(x), as.integer(x)), collapse = " ")
}

#' Build a human-readable summary of a partition
#'
#' This helper summarizes a partition vector by reporting:
#' \itemize{
#'   \item total number of actors;
#'   \item number of distinct groups;
#'   \item sorted group sizes.
#' }
#'
#' It is intended for logging/debugging only and has no effect on model logic.
#'
#' @param p Integer vector encoding a partition (group id per actor).
#' @return Character string describing the partition.
#' @noRd
.erpm_long_partition_info <- function(p) {
  gid <- as.integer(p)
  sizes <- sort(table(gid))

  sprintf(
    "n=%d | groups=%d | sizes={%s}",
    length(gid),
    length(unique(gid)),
    .erpm_long_tab_str(sizes)
  )
}