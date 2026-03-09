################################################################################
# FILE: R/erpm.R
################################################################################
#' ERPM main wrapper: validate inputs, build/translate formula, and optionally fit
#'
#' @name erpm_file
#' @note erpm.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' This file provides the public \code{erpm()} wrapper around \pkg{ergm}.
#'
#' The wrapper is intentionally kept short and orchestration-only. It:
#' \enumerate{
#'   \item validates and normalizes user inputs;
#'   \item resolves the LHS (partition vector vs pre-built bipartite network),
#'         building a bipartite network when needed;
#'   \item translates the RHS ERPM syntax into an \pkg{ergm} RHS;
#'   \item builds a standard \code{ergm()} call with explicit constraints/control;
#'   \item either returns that call (dry-run) or evaluates it.
#' }
#'
#' All user-facing console messages are in English for consistency with \pkg{ergm}.
#'
#' @keywords ERPM ERGM wrapper bipartite translation
NULL
################################################################################

# ==============================================================================
# Public API
# ==============================================================================

#' ERPM main wrapper: translate and optionally fit with \pkg{ergm}
#'
#' @param formula ERGM/ERPM formula. LHS must evaluate either to a partition
#'   vector or to a bipartite \pkg{network} object.
#' @param eval.call Logical. If \code{FALSE}, return the constructed
#'   \code{ergm()} call instead of evaluating it.
#' @param verbose Logical verbosity flag.
#' @param estimate Optional \pkg{ergm} estimation mode.
#' @param eval.loglik Passed to \code{ergm()} when non-NULL.
#' @param control Optional \code{control.ergm} object or list forwarded to
#'   \code{ergm::control.ergm()}.
#' @param timeout Optional timeout in seconds for evaluation.
#' @param seed Optional RNG seed.
#' @param nodes Optional node table or named list coercible to a node table when
#'   the LHS is a partition vector.
#' @param dyads Optional named list of dyadic matrices; for backward compatibility,
#'   a single matrix is also accepted when the RHS contains exactly one dyadic term.
#' @param group_labels Optional group labels passed to
#'   \code{build_bipartite_from_inputs()}.
#' @param constraints Optional constraints formula. Defaults to \code{~ b1part}.
#' @param mh_moves Optional ERPM MH move names forwarded to \code{ErpmMix}.
#' @param mh_weights Optional ERPM MH move weights forwarded to \code{ErpmMix}.
#'
#' @return Either the evaluated \code{ergm()} result or the unevaluated call
#'   when \code{eval.call = FALSE}.
#'
#' @examples
#' \dontrun{
#'   partition <- c(1, 1, 2, 2, 3)
#'   fit <- erpm(partition ~ groups(2) + cliques(k = 2), eval.call = FALSE)
#' }
#'
#' @keywords ERPM ERGM bipartite wrapper
#' @export
erpm <- function(formula,
                 eval.call    = TRUE,
                 verbose      = FALSE,
                 estimate     = NULL,
                 eval.loglik  = NULL,
                 control      = NULL,
                 timeout      = NULL,
                 seed         = NULL,
                 nodes        = NULL,
                 dyads        = list(),
                 group_labels = NULL,
                 constraints  = NULL,
                 mh_moves     = NULL,
                 mh_weights   = NULL) {

  # ---------------------------------------------------------------------------
  # 1) Validate inputs + normalize wrapper settings
  # ---------------------------------------------------------------------------
  validated <- .erpm_validate_inputs(
    formula       = formula,
    eval.call     = eval.call,
    verbose       = verbose,
    estimate      = estimate,
    eval.loglik   = eval.loglik,
    control       = control,
    timeout       = timeout,
    seed          = seed,
    nodes         = nodes,
    dyads         = dyads,
    group_labels  = group_labels,
    constraints   = constraints,
    mh_moves      = mh_moves,
    mh_weights    = mh_weights
  )

  if (isTRUE(validated$verbose) &&
      !is.null(validated$mh_moves) &&
      !is.null(validated$mh_weights)) {
    warning(
      paste(
        "[ERPM] `mh_moves` / `mh_weights` activate the mixed proposal `ErpmMix`.",
        "The supplied weights are interpreted as move attempt weights, not as exact observed frequencies.",
        "Feasibility is checked only after a move type has been drawn, and infeasible",
        "SWAP / MERGE / SPLIT draws are redirected to TOGGLE.",
        "As a result, requested proportions for MERGE and SPLIT may not be respected exactly,",
        "and the observed number of TOGGLE moves can be mechanically inflated."
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 2) Resolve/build LHS network and working formula
  # ---------------------------------------------------------------------------
  resolved <- .erpm_prepare_formula_for_ergm(
    lhs_expr         = validated$lhs_expr,
    rhs_expr         = validated$rhs_expr,
    env0             = validated$env0,
    nodes            = validated$nodes,
    dyads            = validated$dyads,
    group_labels     = validated$group_labels,
    original_formula = formula
  )

  # ---------------------------------------------------------------------------
  # 3) Translate RHS ERPM syntax into an ergm RHS
  # ---------------------------------------------------------------------------
  translated_formula <- .erpm_translate_formula_rhs(
    formula           = resolved$formula,
    env_eval          = resolved$eval_env,
    effect_rename_map = c(),
    wrap_with_proj1   = c(),
    wrap_with_B       = c(),
    verbose           = validated$verbose,
    user_formula_str  = validated$user_formula_str
  )

  # ---------------------------------------------------------------------------
  # 4) Build the final ergm() call
  # ---------------------------------------------------------------------------
  built_call <- .erpm_build_call_bundle(
    formula          = translated_formula,
    eval_env         = resolved$eval_env,
    constraints      = validated$constraints,
    estimate         = validated$estimate,
    eval.loglik      = validated$eval.loglik,
    control          = validated$control,
    timeout          = validated$timeout,
    seed             = validated$seed,
    mh_moves         = validated$mh_moves,
    mh_weights       = validated$mh_weights,
    verbose          = validated$verbose,
    user_formula_str = validated$user_formula_str
  )

  if (isTRUE(validated$verbose) && !isTRUE(validated$eval.call)) {
    cat(
      "\t dry-run ergm call : ",
      paste(deparse(built_call$ergm_call, width.cutoff = 500L), collapse = " "),
      "\n",
      sep = ""
    )
  }

  # ---------------------------------------------------------------------------
  # 5) Evaluate or return
  # ---------------------------------------------------------------------------
  .erpm_eval_or_return(
    ergm_call        = built_call$ergm_call,
    eval.call        = validated$eval.call,
    timeout          = validated$timeout,
    seed             = validated$seed,
    eval_env         = resolved$eval_env,
    user_formula_str = validated$user_formula_str
  )
}