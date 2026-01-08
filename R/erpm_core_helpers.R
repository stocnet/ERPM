################################################################################
# FILE: R/erpm_core_helpers.R
################################################################################
#' ERPM core helpers: shared internals for wrappers
#' @name erpm_core_helpers
#' @note erpm_core_helpers.R
#'
#' @description
#' This file contains small internal helpers shared by multiple ERPM entrypoints,
#' notably \code{erpm()} and \code{erpm_long()}.
#'
#' The intent is to avoid duplication while keeping public APIs stable.
#'
#' @keywords ERPM ERGM internal helpers
NULL

# ============================================================================
# Shared core helpers (internal)
# ============================================================================

#' Parse and normalize the user formula (internal helper)
#' @noRd
.erpm_parse_formula <- function(formula) {
  if (!inherits(formula, "formula"))
    stop("Expected a `lhs ~ ...` formula with lhs = partition OR bipartite network.")

  env0 <- environment(formula)
  if (is.null(env0)) env0 <- parent.frame()

  list(
    env0             = env0,
    lhs_expr         = formula[[2]],
    rhs_expr         = formula[[3]],
    # For logs: keep a compact but readable, single-line representation.
    user_formula_str = .compact_ws(.oneline(formula))
  )
}

#' Evaluate LHS safely (internal helper)
#' @noRd
.erpm_eval_lhs <- function(lhs_expr, eval_env) {
  tryCatch(eval(lhs_expr, envir = eval_env), error = function(e) e)
}

#' Keep original formula while updating its environment (internal helper)
#'
#' Used when the LHS is already a network and no rebuild is needed.
#' @noRd
.erpm_keep_formula <- function(formula, eval_env, lhs_val) {
  new_formula <- formula
  if (!(inherits(lhs_val, "error"))) environment(new_formula) <- eval_env
  new_formula
}

#' Build the final \pkg{ergm} call (internal helper)
#' @noRd
.erpm_build_ergm_call <- function(formula,
                                  constraints,
                                  estimate,
                                  eval.loglik,
                                  control,
                                  verbose_arg_missing,
                                  verbose) {
  call_args <- list(
    as.name("ergm"),
    formula,
    constraints = constraints
  )

  if (!is.null(estimate))    call_args$estimate    <- estimate
  if (!is.null(eval.loglik)) call_args$eval.loglik <- eval.loglik

  # Pass control as an object to avoid hidden state in the evaluation environment.
  # ergm() accepts a control.ergm object directly.
  if (!is.null(control)) call_args$control <- control

  # Propagation rétroactive : si verbose a été explicitement fourni à erpm(),
  # on le transmet aussi à ergm(verbose = <valeur>).
  if (!isTRUE(verbose_arg_missing)) call_args$verbose <- verbose

  as.call(call_args)
}

#' Compact logging for translation (internal helper)
#' @noRd
.erpm_log_translation <- function(user_formula_str,
                                 final_str,
                                 estimate,
                                 eval.loglik,
                                 control,
                                 constraints_str = "~ b1part") {
  cat(sprintf("[ERPM] call initial : erpm(%s) -> call final : %s\n",
              user_formula_str, final_str))
  cat("\t Constraints:  ", constraints_str, "\n", sep = "")
  cat("\t Options: estimate=",
      if (is.null(estimate)) "NULL" else estimate,
      ", eval.loglik=",
      if (is.null(eval.loglik)) "NULL" else eval.loglik,
      ", control=",
      if (is.null(control)) "NULL" else class(control)[1L],
      "\n", sep = "")
}

#' Translate a RHS expression only (internal helper)
#'
#' Returns the translated RHS expression (call/symbol), without touching any formula.
#' This is useful for erpm_long(): history summaries must be computed on translated terms
#' (e.g., groups() -> b2degrange()).
#' @noRd
.erpm_translate_rhs_expr <- function(rhs_expr,
                                    env_eval,
                                    effect_rename_map = c(),
                                    wrap_with_proj1   = c(),
                                    wrap_with_B       = c()) {
  .erpm_translate_rhs_pipeline(
    rhs_expr   = rhs_expr,
    rename_map = effect_rename_map,
    wrap_proj1 = wrap_with_proj1,
    wrap_B     = wrap_with_B,
    env_eval   = env_eval
  )
}