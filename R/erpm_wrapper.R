################################################################################
# FILE: R/erpm_wrapper.R
################################################################################
#' ERPM main wrapper: translate ERPM formulas to \pkg{ergm} calls and optionally fit
#' @name erpm_wrapper
#' @note erpm_wrapper.R
#'
#' @description
#' This module provides the main ERPM wrapper around \pkg{ergm} to:
#' \enumerate{
#'   \item Build a bipartite network from a partition, with optional node and dyadic inputs.
#'   \item Translate ERPM RHS terms (e.g., \code{groups}, \code{cov_match}, \code{cliques})
#'         into \pkg{ergm} terms, with optional encapsulations (\code{Proj1}, \code{B}).
#'   \item Compose a standard call \code{ergm(nw ~ <translated RHS>, constraints = ~ b1part, ...)}.
#'   \item Either return the call (dry-run) or evaluate it and return the fitted model.
#' }
#'
#' All user-facing console messages are in English for consistency with \pkg{ergm}.
#'
#' @note
#' This wrapper assumes that the user has attached the \pkg{ergm} and \pkg{network} packages
#' (or that they are available in the search path) and that the \code{b1part}
#' constraint is meaningful for the constructed network.
#'
#' @note
#' The main wrapper is exercised in self-tests and MWEs under \code{scripts/test}
#' by comparing ERPM-based fits to direct \pkg{ergm} calls on constructed bipartite networks.
#'
#' @keywords ERPM ERGM wrapper bipartite translation

# ============================================================================
# Main ERPM → \pkg{ergm} wrapper
# ============================================================================

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

# ============================================================================
# Small internal pipeline helpers (testable units)
# ============================================================================

#' Parse and normalize the user formula (internal helper)
#' @noRd
.erpm_parse_formula <- function(formula) {
  if (!inherits(formula, "formula"))
    stop("Expected a `lhs ~ ...` formula with lhs = partition OR bipartite network.")

  env0 <- environment(formula) %||% parent.frame()

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

#' Resolve LHS: build network if needed, and build a unique eval_env (internal helper)
#' @noRd
.erpm_resolve_lhs <- function(lhs_val, rhs_expr, env0, nodes, dyads, group_labels = NULL) {
  # Case A: LHS is a partition vector (atomic and not a network)
  if (!(inherits(lhs_val, "error")) &&
      is.atomic(lhs_val) && !inherits(lhs_val, "network")) {

    # Forward group_labels to the builder (optional; no impact if NULL).
    built <- build_bipartite_from_inputs(
      partition    = lhs_val,
      nodes        = nodes,
      dyads        = dyads,
      group_labels = group_labels
    )
    nw2   <- built$network

    # Unique environment used for RHS validation and final evaluation.
    eval_env <- list2env(list(nw = nw2), parent = env0)

    new_formula <- as.formula(bquote(nw ~ .(rhs_expr)))
    environment(new_formula) <- eval_env

    return(list(
      lhs_kind    = "partition",
      eval_env    = eval_env,
      new_formula = new_formula
    ))
  }

  # Case B: LHS is already a network
  if (!(inherits(lhs_val, "error")) && inherits(lhs_val, "network")) {
    bip <- tryCatch(network::get.network.attribute(lhs_val, "bipartite"),
                    error = function(e) NULL)
    if (is.null(bip) || is.na(bip))
      stop("LHS network is not bipartite or missing `%n% 'bipartite'` attribute.")
  }

  # Default: do not fabricate a placeholder formula.
  # The public wrapper will keep the original formula via .erpm_keep_formula().
  list(
    lhs_kind = if (inherits(lhs_val, "network")) "network" else "unknown",
    eval_env = env0
  )
}

# The default constructor above is awkward if left as-is.
# Replace with an explicit helper returning the original formula, while keeping the
# same semantics (no rebuild).
#' @noRd
.erpm_keep_formula <- function(formula, eval_env, lhs_val) {
  new_formula <- formula
  if (!(inherits(lhs_val, "error"))) environment(new_formula) <- eval_env
  new_formula
}

#' Validate one translated term (internal helper)
#' @noRd
.erpm_validate_translated_term <- function(term_call, env_eval) {
  if (is.symbol(term_call)) term_call <- as.call(list(term_call))
  if (!is.call(term_call)) return(term_call)

  f <- term_call[[1L]]

  # Validation for b2degrange(from,to) produced by groups(...)
  if (is.symbol(f) && identical(f, as.name("b2degrange"))) {
    al <- as.pairlist(as.list(term_call)[-1L])
    from <- eval(al$from, envir = env_eval)
    to   <- if (is.null(al$to)) Inf else eval(al$to, envir = env_eval)

    if (is.finite(from) && from < 0L)
      stop(sprintf("groups(from|k): 'from'/'k' must be >= 0. Got: %d", from))

    if (!(is.infinite(to) || (is.finite(to) && to > from))) {
      to_str <- if (is.infinite(to)) "Inf" else as.character(to)
      stop(sprintf("groups(from,to): requires 'from' < 'to'. Got: from=%d, to=%s",
                   from, to_str))
    }
  }

  term_call
}

#' Translate the RHS through a clear pipeline (internal helper)
#' @noRd
.erpm_translate_rhs_pipeline <- function(rhs_expr,
                                        rename_map,
                                        wrap_proj1,
                                        wrap_B,
                                        env_eval) {
  # 1) split
  rhs_terms <- .erpm_split_sum_terms(rhs_expr)

  # 2) translate
  translated <- lapply(
    rhs_terms,
    .erpm_translate_one_term,
    rename_map = rename_map,
    wrap_proj1 = wrap_proj1,
    wrap_B     = wrap_B
  )

  # 3) validate
  translated <- lapply(translated, .erpm_validate_translated_term, env_eval = env_eval)

  # 4) recombine
  if (length(translated) == 1L) translated[[1L]]
  else Reduce(function(x, y) call("+", x, y), translated)
}

#' Translate RHS and inject into formula (internal helper)
#' @noRd
.erpm_translate_rhs <- function(new_formula,
                               eval_env,
                               effect_rename_map,
                               wrap_with_proj1,
                               wrap_with_B,
                               verbose,
                               user_formula_str) {
  rhs_expr <- new_formula[[3]]

  new_rhs <- tryCatch(
    .erpm_translate_rhs_pipeline(
      rhs_expr,
      rename_map = effect_rename_map,
      wrap_proj1 = wrap_with_proj1,
      wrap_B     = wrap_with_B,
      env_eval   = eval_env
    ),
    error = function(e) {
      if (isTRUE(verbose)) {
        cat(sprintf("[ERPM] call initial : erpm(%s)\n", user_formula_str))
        cat("\t error during translation: ", conditionMessage(e), "\n", sep = "")
      }
      stop(e)
    }
  )

  new_formula[[3]] <- new_rhs
  environment(new_formula) <- eval_env
  new_formula
}

#' Build control object (internal helper)
#' @noRd
.erpm_build_control <- function(control, new_formula) {
  if (is.null(control)) return(NULL)

  ctrl <- if (inherits(control, "control.ergm")) control
  else do.call(ergm::control.ergm, as.list(control))

  # Guard: drop `init` when its length does not match the number of stats.
  k <- length(summary(new_formula, constraints = ~ b1part))
  if (!is.null(ctrl$init) && length(ctrl$init) != k) ctrl$init <- NULL

  ctrl
}

#' Evaluate or return call with ERPM error formatting (internal helper)
#' @noRd
.erpm_eval_or_return <- function(ergm_call, eval.call, timeout, eval_env, user_formula_str) {
  if (!isTRUE(eval.call)) return(ergm_call)

  res <- try({
    if (is.null(timeout)) {
      eval(ergm_call, envir = eval_env)
    } else {
      R.utils::withTimeout(
        eval(ergm_call, envir = eval_env),
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

# ============================================================================
# Public wrapper
# ============================================================================

#' ERPM main wrapper: translate and optionally fit with \pkg{ergm}
#'
#' @export
erpm <- function(formula,
                 eval.call    = TRUE,
                 verbose      = TRUE,
                 estimate     = NULL,
                 eval.loglik  = NULL,
                 control      = NULL,
                 timeout      = NULL,
                 nodes        = NULL,
                 dyads        = list(),
                 group_labels = NULL) {

  verbose_arg_missing <- missing(verbose)

  # --- 0) Normalize options ---------------------------------------------------
  if (!is.null(estimate)) {
    estimate <- match.arg(estimate, c("MLE", "CD", "MPLE", "MCMLE"))
    if (identical(estimate, "MCMLE")) estimate <- "MLE"
  }

  # --- 1) Parse formula and build invariants ---------------------------------
  input <- .erpm_parse_formula(formula)
  env0  <- input$env0

  # Single evaluation env is introduced by resolving the LHS.
  lhs_val <- .erpm_eval_lhs(input$lhs_expr, env0)

  # Resolve LHS and create the unique eval_env if partition.
  resolved <- .erpm_resolve_lhs(
    lhs_val, input$rhs_expr, env0, nodes, dyads,
    group_labels = group_labels
  )

  if (identical(resolved$lhs_kind, "partition")) {
    eval_env    <- resolved$eval_env
    new_formula <- resolved$new_formula
  } else {
    eval_env    <- env0
    new_formula <- .erpm_keep_formula(formula, eval_env, lhs_val)
  }

  # --- 2) Translate RHS (pipeline) -------------------------------------------
  effect_rename_map <- c()
  wrap_with_proj1   <- c()
  wrap_with_B       <- c()

  new_formula <- .erpm_translate_rhs(
    new_formula       = new_formula,
    eval_env          = eval_env,
    effect_rename_map = effect_rename_map,
    wrap_with_proj1   = wrap_with_proj1,
    wrap_with_B       = wrap_with_B,
    verbose           = verbose,
    user_formula_str  = input$user_formula_str
  )

  # --- 3) Constraints and control --------------------------------------------
  constraint_expression <- as.formula(~ b1part)
  ctrl <- .erpm_build_control(control, new_formula)

  # --- 4) Build ergm call -----------------------------------------------------
  ergm_call <- .erpm_build_ergm_call(
    formula             = new_formula,
    constraints         = constraint_expression,
    estimate            = estimate,
    eval.loglik         = eval.loglik,
    control             = ctrl,
    verbose_arg_missing = verbose_arg_missing,
    verbose             = verbose
  )

  # --- 5) Logging ------------------------------------------------------------
  if (isTRUE(verbose)) {
    final_fun <- as.call(list(as.name("ergm"), new_formula))
    final_str <- .compact_ws(.oneline(final_fun))
    .erpm_log_translation(
      user_formula_str = input$user_formula_str,
      final_str        = final_str,
      estimate         = estimate,
      eval.loglik      = eval.loglik,
      control          = ctrl,
      constraints_str  = "~ b1part"
    )
  }

  # --- 6) Evaluate or return -------------------------------------------------
  if (!isTRUE(eval.call) && isTRUE(verbose)) {
    cat("\t dry-run ergm call : ",
        paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "),
        "\n", sep = "")
  }

  .erpm_eval_or_return(
    ergm_call         = ergm_call,
    eval.call         = eval.call,
    timeout           = timeout,
    eval_env          = eval_env,
    user_formula_str  = input$user_formula_str
  )
}