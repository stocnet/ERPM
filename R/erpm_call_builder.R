################################################################################
# FILE: R/erpm_call_builder.R
################################################################################
#' ERPM call builder: resolve LHS, build control, compose ergm() call, and evaluate
#'
#' @name erpm_call_builder
#' @note erpm_call_builder.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' This file contains the orchestration helpers used by \code{erpm()} after
#' validation and RHS translation.
#'
#' The module is deliberately structured around four responsibilities:
#' \itemize{
#'   \item resolve the LHS and create the working formula/environment;
#'   \item construct the effective \code{control.ergm} object;
#'   \item build the final \code{ergm()} call and emit compact verbose logs;
#'   \item either return the call or evaluate it with standardized ERPM errors.
#' }
#'
#' This keeps the public wrapper short while preserving historical behavior.
#'
#' @keywords ERPM ERGM call builder control evaluation
NULL
################################################################################

# ==============================================================================
# LHS resolution and working formula
# ==============================================================================

#' Evaluate the LHS safely (internal helper)
#'
#' @param lhs_expr LHS expression from the user formula.
#' @param eval_env Environment where it should be evaluated.
#' @return Either the evaluated object or the caught error object.
#' @noRd
.erpm_eval_lhs <- function(lhs_expr, eval_env) {
  tryCatch(eval(lhs_expr, envir = eval_env), error = function(e) e)
}

#' Keep the original formula while updating its environment (internal helper)
#'
#' @param formula Original formula.
#' @param eval_env Target environment.
#' @param lhs_val Evaluated LHS value.
#' @return Formula.
#' @noRd
.erpm_keep_formula <- function(formula, eval_env, lhs_val) {
  new_formula <- formula
  if (!(inherits(lhs_val, "error"))) environment(new_formula) <- eval_env
  new_formula
}

#' Resolve the LHS and prepare a working formula for ergm()
#'
#' @param lhs_expr LHS expression.
#' @param rhs_expr RHS expression.
#' @param env0 Original formula environment.
#' @param nodes Optional node table.
#' @param dyads Optional dyadic inputs.
#' @param group_labels Optional group labels.
#' @param original_formula Original user formula.
#'
#' @return List with \code{lhs_kind}, \code{eval_env}, and \code{formula}.
#' @noRd
.erpm_prepare_formula_for_ergm <- function(lhs_expr,
                                           rhs_expr,
                                           env0,
                                           nodes,
                                           dyads,
                                           group_labels,
                                           original_formula) {
  lhs_val <- .erpm_eval_lhs(lhs_expr, env0)

  if (!(inherits(lhs_val, "error")) &&
      is.atomic(lhs_val) &&
      !inherits(lhs_val, "network")) {

    built <- build_bipartite_from_inputs(
      partition    = lhs_val,
      nodes        = nodes,
      dyads        = dyads,
      group_labels = group_labels
    )

    nw2 <- built$network
    eval_env <- list2env(list(nw = nw2), parent = env0)

    new_formula <- as.formula(bquote(nw ~ .(rhs_expr)))
    environment(new_formula) <- eval_env

    return(list(
      lhs_kind = "partition",
      eval_env = eval_env,
      formula  = new_formula
    ))
  }

  if (!(inherits(lhs_val, "error")) && inherits(lhs_val, "network")) {
    bip <- tryCatch(
      network::get.network.attribute(lhs_val, "bipartite"),
      error = function(e) NULL
    )

    if (is.null(bip) || is.na(bip)) {
      stop("LHS network is not bipartite or missing `%n% 'bipartite'` attribute.")
    }

    eval_env <- list2env(list(nw = lhs_val), parent = env0)

    new_formula <- as.formula(bquote(nw ~ .(rhs_expr)))
    environment(new_formula) <- eval_env

    return(list(
      lhs_kind = "network",
      eval_env = eval_env,
      formula  = new_formula
    ))
  }

  eval_env <- env0
  new_formula <- .erpm_keep_formula(original_formula, eval_env, lhs_val)

  list(
    lhs_kind = if (inherits(lhs_val, "network")) "network" else "unknown",
    eval_env = eval_env,
    formula  = new_formula
  )
}

# ==============================================================================
# Control / call construction
# ==============================================================================

#' Build the final \pkg{ergm} call (internal helper)
#'
#' @param formula Working formula.
#' @param constraints Effective constraints formula.
#' @param estimate Optional estimation mode.
#' @param eval.loglik Optional evaluation flag.
#' @param control Optional control.ergm object.
#' @param verbose Verbose value.
#'
#' @return Call object.
#' @noRd
.erpm_build_ergm_call <- function(formula,
                                  constraints,
                                  estimate,
                                  eval.loglik,
                                  control,
                                  verbose) {
  call_args <- list(
    as.name("ergm"),
    formula,
    constraints = constraints,
    verbose = verbose
  )

  if (!is.null(estimate))    call_args$estimate    <- estimate
  if (!is.null(eval.loglik)) call_args$eval.loglik <- eval.loglik
  if (!is.null(control))     call_args$control     <- control

  as.call(call_args)
}

#' Build the effective control.ergm object (internal helper)
#'
#' @param control User control input.
#' @param new_formula Translated formula.
#' @param constraints Effective constraints.
#' @param mh_moves Optional ERPM MH move names.
#' @param mh_weights Optional ERPM MH move weights.
#'
#' @return \code{control.ergm} object or NULL.
#' @noRd
.erpm_build_control <- function(control,
                                new_formula,
                                constraints,
                                mh_moves = NULL,
                                mh_weights = NULL) {
  ctrl <- if (is.null(control)) {
    NULL
  } else if (inherits(control, "control.ergm")) {
    control
  } else {
    do.call(ergm::control.ergm, as.list(control))
  }

  if (!is.null(mh_moves) && !is.null(mh_weights)) {
    ctrl_args <- if (is.null(ctrl)) list() else as.list(ctrl)

    ctrl_args$MCMC.prop <- as.formula(~ .select("ErpmMix"))
    ctrl_args$MCMC.prop.args <- list(list(
      moves   = mh_moves,
      weights = mh_weights
    ))

    ctrl <- do.call(ergm::control.ergm, ctrl_args)
  }

  if (is.null(ctrl)) return(NULL)

  k <- length(summary(new_formula, constraints = constraints))
  if (!is.null(ctrl$init) && length(ctrl$init) != k) {
    ctrl$init <- NULL
  }

  ctrl
}

#' Compact logging for translation / call construction (internal helper)
#'
#' @param user_formula_str Original user formula as one line.
#' @param final_str Compact translated \code{ergm(...)} string.
#' @param estimate Estimate mode.
#' @param eval.loglik eval.loglik value.
#' @param control Effective control object.
#' @param constraints_str Constraints deparsed as one line.
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

#' Build the full call bundle used by erpm() (internal helper)
#'
#' @param formula Translated working formula.
#' @param eval_env Evaluation environment.
#' @param constraints Effective constraints formula.
#' @param estimate Optional estimate mode.
#' @param eval.loglik Optional eval.loglik flag.
#' @param control User control input.
#' @param timeout Optional timeout, kept for symmetry with evaluation.
#' @param seed Optional seed, kept for symmetry with evaluation.
#' @param mh_moves Optional ERPM MH move names.
#' @param mh_weights Optional ERPM MH move weights.
#' @param verbose Verbosity flag.
#' @param user_formula_str Original user call as a one-line string.
#'
#' @return List with \code{control} and \code{ergm_call}.
#' @noRd
.erpm_build_call_bundle <- function(formula,
                                    eval_env,
                                    constraints,
                                    estimate,
                                    eval.loglik,
                                    control,
                                    timeout,
                                    seed,
                                    mh_moves,
                                    mh_weights,
                                    verbose,
                                    user_formula_str) {
  ctrl <- .erpm_build_control(
    control     = control,
    new_formula = formula,
    constraints = constraints,
    mh_moves    = mh_moves,
    mh_weights  = mh_weights
  )

  ergm_call <- .erpm_build_ergm_call(
    formula     = formula,
    constraints = constraints,
    estimate    = estimate,
    eval.loglik = eval.loglik,
    control     = ctrl,
    verbose     = verbose
  )

  if (isTRUE(verbose)) {
    final_fun <- as.call(list(as.name("ergm"), formula))
    final_str <- .compact_ws(.oneline(final_fun))

    .erpm_log_translation(
      user_formula_str = user_formula_str,
      final_str        = final_str,
      estimate         = estimate,
      eval.loglik      = eval.loglik,
      control          = ctrl,
      constraints_str  = paste(deparse(constraints, width.cutoff = 500L), collapse = " ")
    )
  }

  list(
    control   = ctrl,
    ergm_call = ergm_call,
    eval_env  = eval_env,
    timeout   = timeout,
    seed      = seed
  )
}

# ==============================================================================
# Evaluation
# ==============================================================================

#' Evaluate or return call with ERPM error formatting (internal helper)
#'
#' @param ergm_call The \code{ergm()} call.
#' @param eval.call Logical: evaluate or return call.
#' @param timeout Optional timeout in seconds.
#' @param seed Optional RNG seed.
#' @param eval_env Evaluation environment.
#' @param user_formula_str User call string for error formatting.
#'
#' @return Either the call or the evaluated result.
#' @noRd
.erpm_eval_or_return <- function(ergm_call,
                                 eval.call,
                                 timeout,
                                 seed,
                                 eval_env,
                                 user_formula_str) {
  if (!isTRUE(eval.call)) return(ergm_call)

  if (!is.null(seed)) {
    if (!(is.numeric(seed) && length(seed) == 1L && is.finite(seed))) {
      stop("[ERPM] `seed` must be a single finite numeric value (integer-like) or NULL.", call. = FALSE)
    }
    si <- as.integer(round(seed))
    if (!isTRUE(all.equal(seed, si))) {
      stop("[ERPM] `seed` must be integer-valued (e.g., 1, 2, 42).", call. = FALSE)
    }
    seed <- si
  }

  res <- try({
    if (is.null(timeout)) {
      if (!is.null(seed)) set.seed(seed)
      eval(ergm_call, envir = eval_env)
    } else {
      R.utils::withTimeout(
        {
          if (!is.null(seed)) set.seed(seed)
          eval(ergm_call, envir = eval_env)
        },
        timeout   = as.numeric(timeout),
        onTimeout = "silent"
      )
    }
  }, silent = TRUE)

  if (inherits(res, "try-error")) {
    msg <- paste0(
      "[ERPM ERROR]\n",
      "  user call : erpm(", user_formula_str, ")\n",
      "  ergm call : ", paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "), "\n",
      "  message   : ", conditionMessage(attr(res, "condition"))
    )
    stop(msg, call. = FALSE)
  }

  res
}