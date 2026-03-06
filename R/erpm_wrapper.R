################################################################################
# FILE: R/erpm_wrapper.R
################################################################################
#' ERPM main wrapper: translate ERPM formulas to \pkg{ergm} calls and optionally fit
#' @name erpm_wrapper
#' @note erpm_wrapper.R
#'
#' @description
#' This file provides the main ERPM wrapper around \pkg{ergm}. It is responsible for:
#' \enumerate{
#'   \item resolving the formula LHS (partition vector vs pre-built bipartite network);
#'   \item building a bipartite network from a partition when needed, with optional
#'         node and dyadic inputs;
#'   \item translating ERPM RHS terms (e.g., \code{groups}, \code{cov_match}, \code{cliques})
#'         into \pkg{ergm} terms, including optional wrappers (\code{Proj1}, \code{B});
#'   \item assembling a standard \code{ergm()} call with explicit constraints and control;
#'   \item optionally injecting an ERPM mixed MH proposal through \code{mh_moves}/\code{mh_weights};
#'   \item either returning the call (dry-run) or evaluating it and returning the fitted model.
#' }
#'
#' All user-facing console messages are in English for consistency with \pkg{ergm}.
#'
#' @note
#' The argument \code{constraints} is supported:
#' \itemize{
#'   \item if NULL, historical behavior is preserved (\code{constraints = ~ b1part});
#'   \item otherwise, it must be a constraints formula forwarded to \code{ergm()}.
#' }
#' This extension is required for PLE meta-networks where block-diagonal constraints may
#' need to be enforced.
#'
#' The arguments \code{mh_moves} and \code{mh_weights} are optional:
#' \itemize{
#'   \item if both are left \code{NULL}, \pkg{ergm} keeps its usual transition logic;
#'   \item if both are provided and valid, \code{erpm()} injects
#'         \code{MCMC.prop = ~ .select("ErpmMix")} together with aligned
#'         \code{MCMC.prop.args}.
#' }
#'
#' @keywords ERPM ERGM wrapper bipartite translation
################################################################################

# ============================================================================
# Bootstrap (dev script support)
# ============================================================================
if (!exists(".erpm_parse_formula", mode = "function") &&
    file.exists("R/erpm_core_helpers.R")) {
  source("R/erpm_core_helpers.R", local = FALSE)
}

# ============================================================================
# Small internal pipeline helpers (testable units)
# ============================================================================

#' Resolve LHS and build a dedicated evaluation environment (internal helper)
#' @param lhs_val Evaluated LHS value (partition vector or network).
#' @param rhs_expr RHS expression.
#' @param env0 Original calling environment.
#' @param nodes Optional node table.
#' @param dyads Optional dyadic inputs.
#' @param group_labels Optional group labels.
#' @return List describing the resolved LHS kind and the updated formula/env.
#' @noRd
.erpm_resolve_lhs <- function(lhs_val, rhs_expr, env0, nodes, dyads, group_labels = NULL) {
  if (!(inherits(lhs_val, "error")) &&
      is.atomic(lhs_val) && !inherits(lhs_val, "network")) {

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
      lhs_kind    = "partition",
      eval_env    = eval_env,
      new_formula = new_formula
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
      lhs_kind    = "network",
      eval_env    = eval_env,
      new_formula = new_formula
    ))
  }

  list(
    lhs_kind = if (inherits(lhs_val, "network")) "network" else "unknown",
    eval_env = env0
  )
}

#' Validate one translated term (internal helper)
#' @param term_call A translated term call.
#' @param env_eval Evaluation environment for term arguments.
#' @return The validated term call (unchanged) or stops on invalid args.
#' @noRd
.erpm_validate_translated_term <- function(term_call, env_eval) {
  if (is.symbol(term_call)) term_call <- as.call(list(term_call))
  if (!is.call(term_call)) return(term_call)

  f <- term_call[[1L]]

  if (is.symbol(f) && identical(f, as.name("b2degrange"))) {
    al <- as.pairlist(as.list(term_call)[-1L])
    from <- eval(al$from, envir = env_eval)
    to   <- if (is.null(al$to)) Inf else eval(al$to, envir = env_eval)

    if (is.finite(from) && from < 0L) {
      stop(sprintf("groups(from|k): 'from'/'k' must be >= 0. Got: %d", from))
    }

    if (!(is.infinite(to) || (is.finite(to) && to > from))) {
      to_str <- if (is.infinite(to)) "Inf" else as.character(to)
      stop(sprintf(
        "groups(from,to): requires 'from' < 'to'. Got: from=%d, to=%s",
        from, to_str
      ))
    }
  }

  term_call
}

#' Translate the RHS through a strict pipeline (internal helper)
#' @param rhs_expr RHS expression.
#' @param rename_map Term rename map.
#' @param wrap_proj1 Wrap terms with Proj1().
#' @param wrap_B Wrap terms with B().
#' @param env_eval Evaluation environment.
#' @return A single RHS expression (call or symbol).
#' @noRd
.erpm_translate_rhs_pipeline <- function(rhs_expr,
                                         rename_map,
                                         wrap_proj1,
                                         wrap_B,
                                         env_eval) {
  rhs_terms <- .erpm_split_sum_terms(rhs_expr)

  translated <- lapply(
    rhs_terms,
    .erpm_translate_one_term,
    rename_map = rename_map,
    wrap_proj1 = wrap_proj1,
    wrap_B     = wrap_B,
    env_eval   = env_eval
  )

  translated <- lapply(
    translated,
    .erpm_validate_translated_term,
    env_eval = env_eval
  )

  if (length(translated) == 1L) translated[[1L]]
  else Reduce(function(x, y) call("+", x, y), translated)
}

#' Translate RHS and inject into formula (internal helper)
#' @param new_formula Working formula.
#' @param eval_env Evaluation environment.
#' @param effect_rename_map Rename map.
#' @param wrap_with_proj1 Wrapper flags.
#' @param wrap_with_B Wrapper flags.
#' @param verbose Verbose flag.
#' @param user_formula_str User formula rendered as a single string.
#' @return Updated formula with translated RHS.
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
#' @param control User control input.
#' @param new_formula Translated formula.
#' @param constraints Effective constraints formula.
#' @param mh_moves Optional ERPM MH move names.
#' @param mh_weights Optional ERPM MH move weights.
#' @return control.ergm object (or NULL).
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

  # Keep init length aligned with the translated formula under the effective constraints.
  k <- length(summary(new_formula, constraints = constraints))
  if (!is.null(ctrl$init) && length(ctrl$init) != k) {
    ctrl$init <- NULL
  }

  ctrl
}

#' Evaluate or return call with ERPM error formatting (internal helper)
#' @param ergm_call The ergm() call.
#' @param eval.call Logical: evaluate or return call.
#' @param timeout Optional timeout in seconds.
#' @param seed Optional RNG seed.
#' @param eval_env Evaluation environment.
#' @param user_formula_str User call string (for error messages).
#' @return Either the call or the evaluated result.
#' @noRd
.erpm_eval_or_return <- function(ergm_call, eval.call, timeout, seed, eval_env, user_formula_str) {
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

# ============================================================================
# Public wrapper
# ============================================================================

#' ERPM main wrapper: translate and optionally fit with \pkg{ergm}
#'
#' (doc omitted here for brevity — keep your current Rd block in your tree)
#'
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

  # Resolve constraints (default remains ~ b1part)
  if (is.null(constraints)) {
    constraints <- as.formula(~ b1part)
  } else {
    if (!(inherits(constraints, "formula") && length(constraints) >= 2L)) {
      stop(
        "[ERPM] `constraints` must be a formula like `~ b1part` or `~ b1part + blockdiag(timeblock)`.",
        call. = FALSE
      )
    }
  }

  verbose_arg_missing <- missing(verbose)

  # --- 0) Normalize options ---------------------------------------------------
  if (!is.null(estimate)) {
    estimate <- match.arg(estimate, c("MLE", "CD", "MPLE", "MCMLE"))
    if (identical(estimate, "MCMLE")) estimate <- "MLE"
  }

  if (!is.null(seed)) {
    if (!(is.numeric(seed) && length(seed) == 1L && is.finite(seed))) {
      stop("[ERPM] `seed` must be a single finite numeric value (integer-like) or NULL.", call. = FALSE)
    }
    seed_i <- as.integer(round(seed))
    if (!isTRUE(all.equal(seed, seed_i))) {
      stop("[ERPM] `seed` must be integer-valued (e.g., 1, 2, 42).", call. = FALSE)
    }
    seed <- seed_i
  }

  # Validate optional MH mix arguments early. When both are NULL, historical
  # behavior is preserved and ergm keeps its own proposal selection logic.
  mh_spec <- .erpm_validate_mh_mix_inputs(
    mh_moves   = mh_moves,
    mh_weights = mh_weights
  )
  mh_moves   <- mh_spec$mh_moves
  mh_weights <- mh_spec$mh_weights

  # When the mixed ERPM proposal is requested, the user-provided weights should
  # be read as attempt weights. They are used for the first-stage move draw,
  # but the observed move frequencies along the chain can differ because
  # infeasible SWAP / MERGE / SPLIT draws fall back to TOGGLE.
  if (isTRUE(verbose) && !is.null(mh_moves) && !is.null(mh_weights)) {
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

  # --- 1) Parse formula and build invariants ---------------------------------
  input <- .erpm_parse_formula(formula)
  env0  <- input$env0

  # --- Dyads normalization ----------------------------------------------------
  if (is.matrix(dyads)) {

    .rhs_dyad_names <- function(rhs_expr) {
      if (is.null(rhs_expr)) return(character(0))

      out <- character(0)

      walk <- function(x) {
        if (is.call(x)) {
          if (is.symbol(x[[1L]])) {
            fn <- as.character(x[[1L]])
            if (startsWith(fn, "dyadcov") || identical(fn, "cov_fullmatch")) {
              if (length(x) >= 2L) {
                a1 <- x[[2L]]
                if (is.character(a1) && length(a1) == 1L && nzchar(a1)) {
                  out <<- c(out, a1)
                }
              }
            }
          }
          for (i in seq_along(x)) walk(x[[i]])
        } else if (is.pairlist(x) || is.list(x)) {
          for (i in seq_along(x)) walk(x[[i]])
        }
        invisible(NULL)
      }

      walk(rhs_expr)
      unique(out)
    }

    nm <- .rhs_dyad_names(input$rhs_expr)
    if (length(nm) != 1L) {
      stop(
        "[ERPM] `dyads` was provided as a matrix, but the RHS does not contain exactly one dyadic name.\n",
        "  Expected something like: dyadcov_full(\"X\") with a unique X.\n",
        "  Fix: pass `dyads = list(X = M)` or ensure the RHS contains one unique dyad name.",
        call. = FALSE
      )
    }

    dyads <- setNames(list(dyads), nm)
  }

  lhs_val <- .erpm_eval_lhs(input$lhs_expr, env0)

  resolved <- .erpm_resolve_lhs(
    lhs_val      = lhs_val,
    rhs_expr     = input$rhs_expr,
    env0         = env0,
    nodes        = nodes,
    dyads        = dyads,
    group_labels = group_labels
  )

  if (identical(resolved$lhs_kind, "partition") ||
      identical(resolved$lhs_kind, "network")) {
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
  ctrl <- .erpm_build_control(
    control     = control,
    new_formula = new_formula,
    constraints = constraints,
    mh_moves    = mh_moves,
    mh_weights  = mh_weights
  )

  # --- 4) Build ergm call -----------------------------------------------------
  ergm_call <- .erpm_build_ergm_call(
    formula             = new_formula,
    constraints         = constraints,
    estimate            = estimate,
    eval.loglik         = eval.loglik,
    control             = ctrl,
    verbose_arg_missing = verbose_arg_missing,
    verbose             = verbose
  )

  # --- 5) Logging -------------------------------------------------------------
  if (isTRUE(verbose)) {
    final_fun <- as.call(list(as.name("ergm"), new_formula))
    final_str <- .compact_ws(.oneline(final_fun))
    .erpm_log_translation(
      user_formula_str = input$user_formula_str,
      final_str        = final_str,
      estimate         = estimate,
      eval.loglik      = eval.loglik,
      control          = ctrl,
      constraints_str  = paste(deparse(constraints, width.cutoff = 500L), collapse = " ")
    )
  }

  # --- 6) Evaluate or return --------------------------------------------------
  if (!isTRUE(eval.call) && isTRUE(verbose)) {
    cat(
      "\t dry-run ergm call : ",
      paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "),
      "\n",
      sep = ""
    )
  }

  .erpm_eval_or_return(
    ergm_call        = ergm_call,
    eval.call        = eval.call,
    timeout          = timeout,
    seed             = seed,
    eval_env         = eval_env,
    user_formula_str = input$user_formula_str
  )
}