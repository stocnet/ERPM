################################################################################
# FILE: R/erpm_utils.R
################################################################################
#' ERPM utilities: small helpers used across the wrapper
#'
#' @name erpm_utils
#' @note erpm_utils.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' This file contains small generic helpers used by ERPM:
#' \itemize{
#'   \item a null-coalescing operator \code{\%||\%};
#'   \item compact deparsing utilities for logging and error messages.
#' }
#'
#' The helpers here are intentionally lightweight and domain-agnostic.
#'
#' @keywords ERPM ERGM wrapper utilities
NULL
################################################################################

# ==============================================================================
# Small generic helpers
# ==============================================================================

#' Null-coalescing operator
#'
#' Return \code{a} if not \code{NULL}, else \code{b}.
#'
#' @param a First value.
#' @param b Fallback value.
#' @return \code{a} when non-NULL, otherwise \code{b}.
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' One-line deparser (compact logging)
#'
#' @param x Object to deparse.
#' @return Single-line character representation.
#' @noRd
.oneline <- function(x) {
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}

#' Compact whitespace for readable logs
#'
#' Replace runs of whitespace by a single space and trim ends.
#' This is preferred for console logs.
#'
#' @param s Character string.
#' @return Compacted character string.
#' @noRd
.compact_ws <- function(s) {
  trimws(gsub("\\s+", " ", s))
}

#' Remove all whitespace from a string
#'
#' Keep this helper for strict comparisons in tests if needed.
#' Avoid using it for user-facing logs.
#'
#' @param s Character string.
#' @return Character string without whitespace.
#' @noRd
.tight <- function(s) {
  gsub("\\s+", "", s)
}