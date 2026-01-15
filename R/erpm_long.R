# =====================================================================
# FILE: R/erpm_long.R
# =====================================================================
#' ERPM longitudinal wrapper: sequential ERGM fits with history covariates
#'
#' @name erpm_long
#' @note erpm_long.R
#'
#' @description
#' This file implements the longitudinal wrapper \code{erpm_long()}.
#' It runs a sequence of ERPM/ERGM fits over time by:
#' \enumerate{
#'   \item evaluating the LHS into a list of partitions (length \eqn{T});
#'   \item detecting which RHS terms are inertial (history-dependent) using an internal registry;
#'   \item building a bipartite membership network for each time \eqn{t};
#'   \item attaching history covariates as network attributes for each required lag
#'         (one attribute per lag, with no aggregation across lags);
#'   \item calling \code{erpm()} at each \eqn{t} with a time-specific RHS composed of
#'         static terms plus the inertial terms that are active at that time.
#' }
#'
#' @keywords ERPM ERGM longitudinal
NULL

# ============================================================================
# Bootstrap (dev script support)
# ============================================================================
# Some selftests may source this file directly.
# Ensure shared helpers are available.
if (!exists(".erpm_parse_formula", mode = "function") &&
    file.exists("R/erpm_core_helpers.R")) {
  # When this file is sourced standalone (outside the package load process),
  # pull in the shared formula parsing helpers from the repository.
  source("R/erpm_core_helpers.R", local = FALSE)
}

if (!exists(".erpm_long_dbg", mode = "function") &&
    file.exists("R/erpm_long_inertia_utils.R")) {
  # Low-level debug/printing helpers used by the inertia subsystem.
  source("R/erpm_long_inertia_utils.R", local = FALSE)
}

if (!exists(".erpm_long_detect_inertial_calls", mode = "function") &&
    file.exists("R/erpm_long_inertia_detect.R")) {
  # RHS parser that splits static terms vs inertial terms using the registry.
  source("R/erpm_long_inertia_detect.R", local = FALSE)
}

if (!exists(".erpm_long_compute_static_summary_one", mode = "function") &&
    file.exists("R/erpm_long_inertia_static_summary.R")) {
  # Fallback helper to compute static summaries on a past network when needed.
  source("R/erpm_long_inertia_static_summary.R", local = FALSE)
}

if (!exists(".erpm_long_inertia_registry", mode = "list") &&
    file.exists("R/erpm_long_inertia_registry.R")) {
  # Declarative registry of inertial terms (how many lags, how to build attributes, etc.).
  source("R/erpm_long_inertia_registry.R", local = FALSE)
}

# erpm_long helpers/engine (dev script support)
if (!exists(".erpm_long_run", mode = "function") &&
    file.exists("R/erpm_long_engine.R")) {
  # Main engine: loop over t, build networks, attach inertia, call erpm().
  source("R/erpm_long_engine.R", local = FALSE)
}
if (!exists(".erpm_long_get_t", mode = "function") &&
    file.exists("R/erpm_long_helpers.R")) {
  # Shared small helpers for indexing per-time inputs (nodes/dyads/labels).
  source("R/erpm_long_helpers.R", local = FALSE)
}

#' ERPM longitudinal wrapper: sequential fits with inertial network attributes
#'
#' Fits a sequence of ERPM models over \eqn{t = 1,\dots,T} from a list of partitions on the LHS.
#' The RHS is shared across time, but \emph{inertial terms} are activated only when enough past
#' partitions exist, based on each term's \code{past_influence}.
#'
#' The workflow is:
#' \enumerate{
#'   \item Parse \code{formula} and evaluate its LHS in the formula environment to obtain a list
#'         of partitions \code{parts} (length \eqn{T}).
#'   \item Split the shared RHS into \emph{static terms} (always active) and \emph{inertial terms}
#'         (activation depends on time) using an internal registry (see
#'         \code{R/erpm_long_inertia_specs.R}).
#'   \item For each time point \code{t}:
#'         \enumerate{
#'           \item build a padded bipartite membership \pkg{network} from \code{parts[[t]]} via
#'                 \code{build_bipartite_from_inputs()}, optionally using time-varying
#'                 \code{nodes}, \code{dyads}, and \code{group_labels};
#'           \item determine which inertial terms are active (a term with \code{past_influence=d}
#'                 is active only when \code{t > d});
#'           \item when active, compute the required lag-specific summaries on past partitions and
#'                 attach \eqn{d} distinct network attributes to the current network (one per lag,
#'                 no aggregation across lags);
#'           \item construct \code{RHS(t)} as static terms plus inertial terms active at \code{t};
#'           \item call \code{erpm()} on \code{nw ~ RHS(t)} either in dry-run mode
#'                 (\code{eval.call=FALSE}) or fit mode (\code{eval.call=TRUE}).
#'         }
#' }
#'
#' The returned object preserves the per-time networks, calls, fitted models, and a timeline
#' describing inertial activation and attribute attachment.
#'
#' @param formula A formula \code{parts ~ <ERPM terms>}, where \code{parts} evaluates to a list of
#'   partitions (one per time point). Each partition must be a non-empty atomic vector of group ids.
#' @param eval.call Logical. If TRUE, evaluate each resulting \code{ergm()} call. If FALSE, return the
#'   unevaluated calls (dry-run) for all time points.
#' @param verbose Logical. If TRUE, print a compact per-time log plus a header summarizing the shared RHS,
#'   options, and inertial detection. When explicitly provided, its value is also forwarded to
#'   \code{ergm(verbose = ...)} via \code{erpm()}.
#' @param debug Logical or character. If FALSE, no extra diagnostics. If TRUE, emit basic diagnostics.
#'   If \code{"deep"}, emit basic diagnostics plus deep data checks (e.g. network summaries and
#'   nodecov/nodefactor input readiness on the translated RHS).
#' @param estimate Character or NULL. Forwarded to \code{erpm()} (and then to \code{ergm(estimate = ...)}),
#'   after the same light normalization used by \code{erpm()}.
#' @param eval.loglik Logical or NULL. Forwarded to \code{erpm()} (and then to \code{ergm(eval.loglik = ...)}).
#' @param control A list, a \code{control.ergm} object, or NULL. Forwarded to \code{erpm()} at each time point.
#' @param timeout Numeric seconds or NULL. If set, each evaluation is run under
#'   \code{R.utils::withTimeout(onTimeout="silent")} in \code{erpm()}.
#' @param seed Integer or NULL. If provided, validated as integer-like and forwarded to \code{erpm()} for each
#'   time point for reproducibility.
#' @param nodes Optional node data. May be:
#'   \itemize{
#'     \item NULL (auto labels),
#'     \item a single \code{data.frame} shared across time,
#'     \item a list of \code{data.frame}s (one per time point).
#'   }
#' @param dyads Optional dyadic inputs. Accepts matrices or lists of matrices and is normalized internally to
#'   a time-indexed list structure consistent with the longitudinal engine. These inputs are attached as
#'   dedicated network attributes during network construction.
#' @param group_labels Optional group vertex labels. May be time-invariant or time-indexed; forwarded to
#'   \code{build_bipartite_from_inputs()} for readable group labels in the padded bipartite networks.
#'
#' @return An object of class \code{"erpm_long"} with components:
#'   \itemize{
#'     \item \code{calls}: list of per-time \code{ergm()} calls (or extracted calls from fits);
#'     \item \code{fits}: list of per-time fitted models (or NULLs in dry-run mode);
#'     \item \code{networks}: list of per-time networks used for fitting;
#'     \item \code{history_timeline}: list describing inertial activation and attached attributes per time.
#'   }
#'
#' @examples
#' \dontrun{
#' parts <- list(
#'   c(1, 1, 2, 2, 3),
#'   c(1, 2, 2, 3, 3),
#'   c(1, 1, 1, 2, 3)
#' )
#'
#' # Dry-run: inspect per-time calls and inertial activation
#' out_calls <- erpm_long(parts ~ groups + inertia_term(past_influence = 2), eval.call = FALSE)
#'
#' # Fit sequentially with shared RHS, inertial term activates when eligible
#' out_fits <- erpm_long(parts ~ groups + inertia_term(past_influence = 2), estimate = "MLE")
#' }
#'
#' @note
#' Inertial terms are detected and handled via \code{R/erpm_long_inertia_specs.R}. A term requesting
#' \code{past_influence = d} is active only when \code{t > d}. When active, the wrapper materializes
#' \eqn{d} distinct lag-specific network attributes on the current network (one per lag), without
#' combining information across lags.
#'
#' @note
#' This wrapper delegates model translation, constraint handling, and fitting to \code{erpm()} at each
#' time point, so static and longitudinal workflows share the same call structure and safeguards.
#'
#' @export
erpm_long <- function(formula,
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
                      group_labels = NULL) {

  # ---------------------------------------------------------------------------
  # 0) Basic checks + standard parsing
  # ---------------------------------------------------------------------------
  if (!inherits(formula, "formula")) {
    # We require an R formula so we can evaluate the LHS in its environment and
    # keep the RHS as an expression to be analyzed (static vs inertial terms).
    stop("[ERPM_LONG] `formula` must be a formula: partitions ~ rhs", call. = FALSE)
  }

  # Normalize debug flag:
  # - FALSE: off
  # - TRUE : basic debug traces
  # - "deep": basic + deep data dumps (nodecov/nodefactor inputs, etc.)
  if (!(isFALSE(debug) || isTRUE(debug) || (is.character(debug) && identical(debug, "deep")))) {
    stop("[ERPM_LONG] `debug` must be FALSE, TRUE, or \"deep\".", call. = FALSE)
  }
  debug_deep <- is.character(debug) && identical(debug, "deep")
  debug_any  <- isTRUE(debug) || debug_deep

  .dbg <- function(...) {
    # Local debug printer.
    # Kept minimal so we can pass it down to the engine without pulling global state.
    if (isTRUE(debug_any)) cat("[ERPM_LONG][debug] ", ..., "\n", sep = "")
  }

  # If `seed` is not NULL, ensure it is compatible with set.seed().
  # Keep the normalized integer value (so downstream code never has to repeat checks).
  if (!is.null(seed)) {
    if (!(is.numeric(seed) && length(seed) == 1L && is.finite(seed))) {
      stop("[ERPM_LONG] `seed` must be a single finite numeric value (integer-like) or NULL.", call. = FALSE)
    }
    seed_i <- as.integer(round(seed))
    if (!isTRUE(all.equal(seed, seed_i))) {
      stop("[ERPM_LONG] `seed` must be integer-valued (e.g., 1, 2, 42).", call. = FALSE)
    }
    seed <- seed_i
  }

  input <- .erpm_parse_formula(formula)
  env0  <- input$env0
  rhs_expr <- input$rhs_expr

  # Evaluate the LHS in the original formula environment.
  # We expect a list of partitions, one per time point.
  parts <- try(eval(input$lhs_expr, envir = env0), silent = TRUE)
  if (inherits(parts, "try-error")) {
    stop("[ERPM_LONG] LHS must evaluate to a list of partitions.", call. = FALSE)
  }
  if (!is.list(parts)) {
    stop("[ERPM_LONG] LHS must be a list of partitions (list of integer vectors).", call. = FALSE)
  }

  # Number of time points.
  T <- length(parts)
  if (T < 1L) stop("[ERPM_LONG] Empty partition list.", call. = FALSE)
  if (T == 1L) {
    # If there is only one partition, longitudinal logic is unnecessary and
    # inertia cannot be defined in a meaningful way.
    stop("[ERPM_LONG] A single partition was provided. Use `erpm()` instead of `erpm_long()`.", call. = FALSE)
  }

  for (t in seq_len(T)) {
    p <- parts[[t]]

    # Each partition must be a non-empty atomic vector (typically integer-like group ids).
    if (is.null(p) || !is.atomic(p) || length(p) < 1L) {
      stop(sprintf("[ERPM_LONG] partitions[[%d]] must be a non-empty atomic vector.", t), call. = FALSE)
    }
  }

  # erpm_long delegates the actual ERGM call construction and fitting to erpm().
  if (!exists("erpm", mode = "function")) {
    stop("[ERPM_LONG] erpm() is not available in this session.", call. = FALSE)
  }
  # The longitudinal wrapper requires the bipartite builder because each time step
  # starts from a partition and is converted into a membership network.
  if (!exists("build_bipartite_from_inputs", mode = "function")) {
    stop("[ERPM_LONG] build_bipartite_from_inputs() is required for erpm_long().", call. = FALSE)
  }

  # Inertia helpers must be available for this version of erpm_long().
  if (!exists(".erpm_long_detect_inertial_calls", mode = "function")) {
    stop("[ERPM_LONG] missing inertia helpers. Source R/erpm_long_inertia_specs.R.", call. = FALSE)
  }
  if (!exists(".erpm_long_run", mode = "function")) {
    stop("[ERPM_LONG] missing engine helpers. Source R/erpm_long_engine.R.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 1) Helpers for per-time inputs
  # ---------------------------------------------------------------------------
  if (!is.null(nodes) && is.list(nodes) && !is.data.frame(nodes) && !.erpm_long_is_list_of_df(nodes)) {
    # `nodes` can be:
    # - NULL (no node table),
    # - one data.frame (shared across time),
    # - a list of data.frames (one per time point).
    stop("[ERPM_LONG] `nodes` must be a data.frame or a list of data.frames.", call. = FALSE)
  }

  # Normalize dyads early.
  # Accept: matrix, named list of matrices, list(T) of matrices, list(T) of named lists.
  dyads <- .erpm_long_normalize_dyads_input(dyads, rhs_expr = rhs_expr, T = T)

  # ---------------------------------------------------------------------------
  # 2) Detect inertial calls once (RHS shared, activation depends on t)
  # ---------------------------------------------------------------------------
  # The RHS expression is the same for all t, but inertial terms can be inactive
  # at early times depending on how many past partitions they require.
  det <- .erpm_long_detect_inertial_calls(rhs_expr, env0, debug = debug_any)
  static_terms     <- det$static_terms
  inertial_calls <- det$inertial_calls

  # ---------------------------------------------------------------------------
  # 2b) User-facing notification: inertial terms that can never activate
  # ---------------------------------------------------------------------------
  # If an inertial term requests past_influence=d, it is active only when t > d.
  # Therefore, for a timeline of length T, the term is never active if T <= d.
  if (isTRUE(det$enabled) && length(inertial_calls)) {
    never_active <- character(0)
    info_lines   <- character(0)

    for (cc in inertial_calls) {
      nm <- .erpm_long_term_name(cc)
      spec <- .erpm_long_inertia_registry[[nm]]
      if (is.null(spec) || !is.function(spec$past_influence)) next

      d <- spec$past_influence(cc, env0)
      d <- as.integer(d)

      start_t <- d + 1L
      if (start_t > T) {
        never_active <- c(never_active, sprintf("%s(past_influence=%d, T=%d)", nm, d, T))
      }

      info_lines <- c(info_lines, sprintf("  - %s: past_influence=%d -> active for t >= %d", nm, d, start_t))
    }

    if (length(never_active)) {
      warning(
        paste0(
          "[ERPM_LONG] Some inertial terms will never be applied with the provided partitions.\n",
          "  Rule: a term with past_influence=d is active only when t > d.\n",
          sprintf("  Here: T=%d, so terms with d >= %d are never active.\n", T, T),
          "  Affected term(s): ", paste(unique(never_active), collapse = ", "), "\n",
          "  Fix: provide more partitions or reduce past_influence."
        ),
        call. = FALSE
      )
    }

    if (isTRUE(verbose) && length(info_lines)) {
      cat("[ERPM_LONG] inertial activation schedule (per term):\n")
      cat(paste(info_lines, collapse = "\n"), "\n", sep = "")
    }
  }

  # ---------------------------------------------------------------------------
  # 3) Verbose header
  # ---------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    # Human-friendly header so the user can see what will happen before the loop runs.
    cat("============================================================\n")
    cat("[ERPM_LONG] Start\n")
    cat(sprintf("[ERPM_LONG] T = %d partitions\n", T))
    cat(sprintf("[ERPM_LONG] eval.call = %s | estimate = %s | eval.loglik = %s | seed = %s\n",
                if (isTRUE(eval.call)) "TRUE" else "FALSE",
                if (is.null(estimate)) "NULL" else as.character(estimate),
                if (is.null(eval.loglik)) "NULL" else as.character(eval.loglik),
                if (is.null(seed)) "NULL" else as.character(seed)))

    # Print the shared RHS once (useful for debugging formula parsing).
    cat("[ERPM_LONG] RHS (shared): ", .erpm_long_rhs_oneline(rhs_expr), "\n", sep = "")
    if (isTRUE(det$enabled)) {
      cat("[ERPM_LONG] inertial terms detected: ", paste(det$inertial_names, collapse = ", "), "\n", sep = "")
      # Important behavioral point: we do not pool information across lags, each lag is explicit.
      cat("[ERPM_LONG] note: activation is per-term via past_influence; lags are not reduced.\n")
    } else {
      cat("[ERPM_LONG] inertial terms detected: none\n")
    }

    if (debug_deep) {
      cat("[ERPM_LONG] debug level: deep\n")
    } else if (debug_any) {
      cat("[ERPM_LONG] debug level: basic\n")
    }

    cat("------------------------------------------------------------\n")
  }

  # ---------------------------------------------------------------------------
  # 4) Run
  # ---------------------------------------------------------------------------
  .erpm_long_run(
    parts          = parts,
    T              = T,
    env0           = env0,
    rhs_expr       = rhs_expr,
    static_terms   = static_terms,
    inertial_calls = inertial_calls,
    det            = det,
    eval.call      = eval.call,
    verbose        = verbose,
    debug          = debug,     # propagate raw debug value (FALSE/TRUE/"deep")
    .dbg           = .dbg,
    estimate       = estimate,
    eval.loglik    = eval.loglik,
    control        = control,
    timeout        = timeout,
    seed           = seed,
    nodes          = nodes,
    dyads          = dyads,
    group_labels   = group_labels
  )
}