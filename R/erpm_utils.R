################################################################################
# FILE: R/erpm_utils.R
################################################################################
#' ERPM utilities: small helpers used across the wrapper
#' @name erpm_utils
#' @note erpm_utils.R
#'
#' @description
#' This file contains small generic helpers used by ERPM :
#' \itemize{
#'   \item a null-coalescing operator \code{\%||\%};
#'   \item compact deparsing utilities for logging and error messages.
#' }
#'
#' @keywords ERPM ERGM wrapper utilities

# ============================================================================
# Small generic helpers
# ============================================================================

#' Null-coalescing operator
#'
#' Return `a` if not `NULL`, else `b`.
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' One-line deparser (compact logging)
#' @noRd
.oneline <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

#' Compact whitespace for readable logs
#'
#' Replace runs of whitespace by a single space and trim ends.
#' This is preferred for console logs.
#' @noRd
.compact_ws <- function(s) trimws(gsub("\\s+", " ", s))

#' Remove all whitespace from a string
#'
#' Keep this helper for strict comparisons in tests if needed.
#' Avoid using it for user-facing logs.
#' @noRd
.tight <- function(s) gsub("\\s+", "", s)