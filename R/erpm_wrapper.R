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
#' A new argument `constraints` is supported.
#' - If provided, it is forwarded to ergm().
#' - Otherwise the historical default is used: constraints = ~ b1part.
#' This change is required for PLE/stacked meta-networks where blockdiag() must be enforced.
#'
#' @keywords ERPM ERGM wrapper bipartite translation

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

#' Resolve LHS: build network if needed, and build a unique eval_env (internal helper)
#' @noRd
.erpm_legacy_resolve_lhs <- function(lhs_val, rhs_expr, env0, nodes, dyads, group_labels = NULL) {
  if (!(inherits(lhs_val, "error")) &&
      is.atomic(lhs_val) && !inherits(lhs_val, "network")) {

    built <- build_bipartite_from_inputs(
      partition    = lhs_val,
      nodes        = nodes,
      dyads        = dyads,
      group_labels = group_labels
    )
    nw2   <- built$network
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
    bip <- tryCatch(network::get.network.attribute(lhs_val, "bipartite"),
                    error = function(e) NULL)
    if (is.null(bip) || is.na(bip))
      stop("LHS network is not bipartite or missing `%n% 'bipartite'` attribute.")

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
#' @noRd
.erpm_legacy_validate_translated_term <- function(term_call, env_eval) {
  if (is.symbol(term_call)) term_call <- as.call(list(term_call))
  if (!is.call(term_call)) return(term_call)

  f <- term_call[[1L]]

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
.erpm_legacy_translate_rhs_pipeline <- function(rhs_expr,
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

  translated <- lapply(translated, .erpm_legacy_validate_translated_term, env_eval = env_eval)

  if (length(translated) == 1L) translated[[1L]]
  else Reduce(function(x, y) call("+", x, y), translated)
}

#' Translate RHS and inject into formula (internal helper)
#' @noRd
.erpm_legacy_translate_rhs <- function(new_formula,
                                      eval_env,
                                      effect_rename_map,
                                      wrap_with_proj1,
                                      wrap_with_B,
                                      verbose,
                                      user_formula_str) {
  rhs_expr <- new_formula[[3]]

  new_rhs <- tryCatch(
    .erpm_legacy_translate_rhs_pipeline(
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
.erpm_legacy_build_control <- function(control, new_formula, constraints) {
  if (is.null(control)) return(NULL)

  ctrl <- if (inherits(control, "control.ergm")) control
  else do.call(ergm::control.ergm, as.list(control))

  # -------------------------------------------------------------------------
  # CHANGE (justified):
  # Previously we hard-coded constraints=~b1part when computing k=number of stats.
  # With the new `constraints` argument, we must compute k under the *effective*
  # constraints, otherwise we may incorrectly drop/init or keep an incompatible init.
  # This is backward compatible because constraints defaults to ~b1part.
  # -------------------------------------------------------------------------
  k <- length(summary(new_formula, constraints = constraints))
  if (!is.null(ctrl$init) && length(ctrl$init) != k) ctrl$init <- NULL

  ctrl
}

#' Evaluate or return call with ERPM error formatting (internal helper)
#' @noRd
.erpm_legacy_eval_or_return <- function(ergm_call, eval.call, timeout, seed, eval_env, user_formula_str) {
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
#' @noRd
.erpm_legacy_wrapper <- function(formula,
                                 eval.call    = TRUE,
                                 verbose      = TRUE,
                                 debug        = FALSE,
                                 estimate     = NULL,
                                 eval.loglik  = NULL,
                                 control      = NULL,
                                 timeout      = NULL,
                                 seed         = NULL,
                                 nodes        = NULL,
                                 dyads        = list(),
                                 group_labels = NULL,
                                 constraints  = NULL) {

  # -------------------------------------------------------------------------
  # CHANGE (justified):
  # New argument `constraints`:
  # - if NULL: keep historical behavior constraints = ~ b1part
  # - else: must be a formula like ~ b1part + blockdiag(timeblock)
  # This is the minimal extension required to let erpm_long (PLE) enforce blockdiag.
  # -------------------------------------------------------------------------
  if (is.null(constraints)) {
    constraints <- as.formula(~ b1part)
  } else {
    if (!(inherits(constraints, "formula") && length(constraints) >= 2L)) {
      stop("[ERPM] `constraints` must be a formula like `~ b1part` or `~ b1part + blockdiag(timeblock)`.", call. = FALSE)
    }
  }

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

  # --- 1) Parse formula and build invariants ---------------------------------
  input <- .erpm_parse_formula(formula)
  env0  <- input$env0

  # ---  Dyads normalization -----------------------------------
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
                if (is.character(a1) && length(a1) == 1L && nzchar(a1)) out <<- c(out, a1)
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

  resolved <- .erpm_legacy_resolve_lhs(
    lhs_val, input$rhs_expr, env0, nodes, dyads,
    group_labels = group_labels
  )

  if (identical(resolved$lhs_kind, "partition") || identical(resolved$lhs_kind, "network")) {
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

  new_formula <- .erpm_legacy_translate_rhs(
    new_formula       = new_formula,
    eval_env          = eval_env,
    effect_rename_map = effect_rename_map,
    wrap_with_proj1   = wrap_with_proj1,
    wrap_with_B       = wrap_with_B,
    verbose           = verbose,
    user_formula_str  = input$user_formula_str
  )

  # --- 3) Constraints and control --------------------------------------------
  ctrl <- .erpm_legacy_build_control(control, new_formula, constraints = constraints)

  # --- 4) Build ergm call -----------------------------------------------------
  ergm_call <- .erpm_build_ergm_call(
    formula             = new_formula,
    constraints         = constraints,
    estimate            = estimate,
    eval.loglik         = eval.loglik,
    control             = ctrl,
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
      # CHANGE (justified): constraints string must reflect the effective constraints.
      constraints_str  = paste(deparse(constraints, width.cutoff = 500L), collapse = " ")
    )
  }

  # --- 6) Evaluate or return -------------------------------------------------
  # Propagate debug flag to constraint-level option so InitErgmConstraint.b1partblockdiag
  # prints its diagnostic messages during the ergm() call below.
  .old_b1bd_dbg <- getOption("ERPM.b1partblockdiag.debug", FALSE)
  if (isTRUE(debug)) options(ERPM.b1partblockdiag.debug = TRUE)
  on.exit(options(ERPM.b1partblockdiag.debug = .old_b1bd_dbg), add = TRUE)

  if (!isTRUE(eval.call) && isTRUE(verbose)) {
    cat("\t dry-run ergm call : ",
        paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "),
        "\n", sep = "")
  }

  .erpm_legacy_eval_or_return(
    ergm_call         = ergm_call,
    eval.call         = eval.call,
    timeout           = timeout,
    seed              = seed,
    eval_env          = eval_env,
    user_formula_str  = input$user_formula_str
  )
}