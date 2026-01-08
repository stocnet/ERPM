################################################################################
# FILE: R/erpm_long_inertia_static_summary.R
################################################################################
#' ERPM longitudinal inertia: static summary computation on past networks
#'
#' @name erpm_long_inertia_static_summary
#' @note erpm_long_inertia_static_summary.R
#'
#' @description
#' This file defines a low-level helper used by \code{erpm_long()} to compute
#' static ERGM summary statistics on *past* networks.
#'
#' The function implemented here is only used when an inertial term does NOT
#' provide its own \code{build_attr()} method. In that case, inertia is handled
#' by:
#' \enumerate{
#'   \item taking one past network \code{nw_prev};
#'   \item building a temporary ERGM formula \code{nw_prev ~ <static terms>};
#'   \item calling \code{summary()} under the ERPM bipartite constraint
#'         \code{~ b1part};
#'   \item returning the resulting statistics as a named numeric vector.
#' }
#'
#' This mechanism allows inertial effects to reuse existing ERGM terms
#' transparently, without duplicating changestat logic in the longitudinal
#' layer.
#'
#' @keywords ERPM ERGM longitudinal inertia internal
NULL

# =============================================================================
# Helpers
# =============================================================================

#' Compute summary statistics for a list of ERGM term calls on one network
#'
#' @description
#' Internal helper used by \code{erpm_long()} to evaluate ERGM statistics on a
#' previously constructed network. This is typically invoked for inertial
#' effects that rely on *static summaries* of past partitions rather than
#' custom-built attributes.
#'
#' The function:
#' \enumerate{
#'   \item combines a list of ERGM RHS terms into a single RHS expression;
#'   \item builds a temporary formula of the form \code{nw ~ RHS};
#'   \item evaluates \code{summary()} under constraint \code{~ b1part};
#'   \item returns a named numeric vector of statistics.
#' }
#'
#' If no static terms are provided, an empty named numeric vector is returned.
#'
#' @param nw_prev A \pkg{network} object corresponding to a past partition.
#' @param static_terms A list of ERGM RHS calls (each a call or symbol).
#' @param env0 The original evaluation environment of the ERPM formula.
#'
#' @return A named numeric vector of ERGM summary statistics.
#'
#' @noRd
.erpm_long_compute_static_summary_one <- function(nw_prev, static_terms, env0) {

  # ---------------------------------------------------------------------------
  # Trivial case: no static terms requested
  # ---------------------------------------------------------------------------
  # This happens when an inertial effect declares no summaries to compute.
  # We return an empty named numeric vector for consistency.
  if (!length(static_terms)) {
    return(setNames(numeric(0), character(0)))
  }

  # ---------------------------------------------------------------------------
  # Combine multiple RHS terms into a single ERGM RHS expression
  # ---------------------------------------------------------------------------
  # If only one term is present, keep it as-is.
  # Otherwise, reduce the list using '+' calls to mimic a standard ERGM RHS.
  rhs_static <- if (length(static_terms) == 1L) {
    static_terms[[1L]]
  } else {
    Reduce(function(x, y) call("+", x, y), static_terms)
  }

  # ---------------------------------------------------------------------------
  # Build an evaluation environment where 'nw' refers to the past network
  # ---------------------------------------------------------------------------
  # We intentionally avoid touching the global environment.
  # The parent environment remains env0 so that symbols and user-defined
  # objects used in the terms can still be resolved.
  eval_env_prev <- list2env(list(nw = nw_prev), parent = env0)

  # Construct the temporary formula: nw ~ <static RHS>
  f_prev <- as.formula(bquote(nw ~ .(rhs_static)))
  environment(f_prev) <- eval_env_prev

  # ---------------------------------------------------------------------------
  # Compute ERGM summary under the ERPM bipartite constraint
  # ---------------------------------------------------------------------------
  # We suppress messages because ergm::summary() can be verbose and because
  # non-dyad-independent constraints (b1part) often trigger warnings that are
  # expected in this context.
  s_prev <- try(
    suppressMessages(
      summary(f_prev, constraints = ~ b1part)
    ),
    silent = TRUE
  )

  # Propagate a clean error if summary() failed.
  if (inherits(s_prev, "try-error")) {
    msg <- conditionMessage(attr(s_prev, "condition"))
    stop(sprintf("[ERPM_LONG] summary() failed: %s", msg), call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Normalize statistic names and coerce to numeric
  # ---------------------------------------------------------------------------
  # ERGM summary names can sometimes be missing or empty (edge cases).
  # We defensively generate generic names to guarantee a well-formed result.
  sn <- names(s_prev)
  if (is.null(sn) || any(!nzchar(sn))) {
    sn <- paste0("stat_", seq_along(s_prev))
  }

  setNames(as.numeric(s_prev), sn)
}