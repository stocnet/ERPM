# ################################################################################
# # FILE: R/erpm_long_common_rhs.R
# ################################################################################
# #' ERPM longitudinal common RHS helpers (internal)
# #'
# #' @name erpm_long_common_rhs
# #' @note erpm_long_common_rhs.R
# #'
# #' @description
# #' Shared RHS parsing/formatting utilities for \code{erpm_long()}.
# #'
# #' @keywords ERPM ERGM longitudinal internal helpers
# NULL

# #' Split a sum-of-terms RHS recursively
# #' @noRd
# .erpm_long_split_sum_terms <- function(expr) {
#   if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
#     return(c(.erpm_long_split_sum_terms(expr[[2]]),
#              .erpm_long_split_sum_terms(expr[[3]])))
#   }
#   list(expr)
# }

# #' Combine a list of ERGM terms into a single RHS expression
# #' @noRd
# .erpm_long_combine_terms <- function(terms) {
#   if (!length(terms)) return(NULL)
#   if (length(terms) == 1L) return(terms[[1L]])
#   Reduce(function(x, y) call("+", x, y), terms)
# }

# #' Convert a RHS expression into a single-line character string
# #' @noRd
# .erpm_long_rhs_oneline <- function(expr) {
#   paste(deparse(expr, width.cutoff = 500L), collapse = " ")
# }

# #' Extract a term/function name from a RHS item
# #' @noRd
# .erpm_long_term_name <- function(tt) {
#   if (is.symbol(tt)) return(as.character(tt))
#   if (is.call(tt) && is.symbol(tt[[1L]])) return(as.character(tt[[1L]]))
#   NA_character_
# }

# #' Flatten an expression tree to collect all calls matching a predicate
# #' @noRd
# .erpm_long_collect_calls <- function(expr, pred) {
#   out <- list()

#   walk <- function(x) {
#     if (is.call(x)) {
#       if (isTRUE(pred(x))) out[[length(out) + 1L]] <<- x
#       for (i in seq_along(x)) walk(x[[i]])
#     } else if (is.pairlist(x) || is.list(x)) {
#       for (i in seq_along(x)) walk(x[[i]])
#     }
#     invisible(NULL)
#   }

#   walk(expr)
#   out
# }

# #' Collect dyadic covariate names referenced on the RHS (string literal only)
# #'
# #' Conservative parser: only literal strings, no evaluation
# #'
# #' Recognized patterns:
# #' - dyadcov_*("name") and dyadcov("name") (prefix "^dyadcov")
# #' - cov_fullmatch("name")
# #'
# #' @noRd
# .erpm_long_rhs_dyad_names <- function(rhs_expr) {
#   if (is.null(rhs_expr)) return(character(0))

#   is_target <- function(cc) {
#     if (!is.call(cc)) return(FALSE)
#     if (!is.symbol(cc[[1L]])) return(FALSE)
#     fn <- as.character(cc[[1L]])
#     startsWith(fn, "dyadcov") || identical(fn, "cov_fullmatch")
#   }

#   calls <- .erpm_long_collect_calls(rhs_expr, is_target)
#   if (!length(calls)) return(character(0))

#   nm <- character(0)
#   for (cc in calls) {
#     if (length(cc) >= 2L) {
#       a1 <- cc[[2L]]
#       if (is.character(a1) && length(a1) == 1L && nzchar(a1)) nm <- c(nm, a1)
#     }
#   }

#   unique(nm)
# }

# #' Infer a single dyad name from the RHS if possible
# #' @noRd
# .erpm_long_infer_single_dyad_name <- function(rhs_expr) {
#   nm <- .erpm_long_rhs_dyad_names(rhs_expr)
#   if (length(nm) == 1L) nm else NULL
# }