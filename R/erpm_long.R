################################################################################
# FILE: R/erpm_long.R
################################################################################
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
#' This function:
#' \enumerate{
#'   \item Interprets the LHS of a formula as a list of partitions \code{parts} of length \eqn{T},
#'         where each element is a partition vector at time \eqn{t};
#'   \item Splits the shared RHS into \emph{static terms} (always active) and \emph{inertial terms}
#'         (activated only when enough past partitions exist), using an internal declarative
#'         registry (see \code{R/erpm_long_inertia_specs.R});
#'   \item Iterates over time points \eqn{t = 1,\dots,T} and, at each step:
#'         \enumerate{
#'           \item Builds a padded bipartite membership \pkg{network} from the current partition
#'                 (via \code{build_bipartite_from_inputs()}, with optional \code{nodes}, \code{dyads},
#'                 and \code{group_labels});
#'           \item For each inertial term present in the RHS, reads its \code{past_influence = d}
#'                 argument and activates the term only when \eqn{t > d};
#'           \item When an inertial term is active, computes the term-specific required summaries on
#'                 past partitions independently for each lag \eqn{\ell \in \{1,\dots,d\}}, and attaches
#'                 \eqn{d} distinct network attributes to the current network (one attribute per lag,
#'                 with no aggregation across lags);
#'           \item Constructs a time-specific RHS \code{RHS(t)} consisting of static terms plus the inertial
#'                 terms that are active at time \eqn{t};
#'           \item Calls \code{erpm()} on the current network and \code{RHS(t)}, either as a dry-run
#'                 (\code{eval.call = FALSE}) or as an evaluated \pkg{ergm} fit (\code{eval.call = TRUE}).
#'         }
#'   \item Returns the list of \eqn{T} time-specific results (calls or fitted models), reflecting the fact
#'         that inertial terms are intentionally absent at early time points until their
#'         \code{past_influence} requirement is met.
#' }
#'
#' @param formula A formula \code{parts ~ <ERPM terms>}, where \code{parts} evaluates to a list of
#'   partitions (one per time point). Each partition is an atomic vector of group ids.
#' @param eval.call Logical. If TRUE, evaluate each resulting \code{ergm()} call. If FALSE, return the
#'   unevaluated calls (dry-run) for all time points.
#' @param verbose Logical. If TRUE, print a compact per-time translation log and the effective options.
#'   When explicitly provided, its value is also forwarded to \code{ergm(verbose = ...)} via \code{erpm()}.
#' @param debug Logical. If TRUE, emit additional diagnostic messages about inertial-term activation and
#'   attribute construction.
#' @param estimate Character or NULL. Forwarded to \code{erpm()} (and then to \code{ergm(estimate = ...)}),
#'   after the same light normalization used by \code{erpm()} (e.g., mapping \code{"MCMLE"} to \code{"MLE"}).
#' @param eval.loglik Logical or NULL. Forwarded to \code{erpm()} (and then to \code{ergm(eval.loglik = ...)}).
#' @param control A list, a \code{control.ergm} object, or NULL. Forwarded to \code{erpm()} for each time
#'   point, with the same safeguards regarding incompatible \code{init} lengths.
#' @param timeout Numeric seconds or NULL. If set, each evaluation is run under \code{R.utils::withTimeout()}
#'   with \code{onTimeout="silent"}.
#' @param nodes Optional \code{data.frame} for actor attributes and labels, used when building each bipartite
#'   network from a partition.
#' @param dyads Optional named list of \eqn{n\times n} matrices to attach to each network as a dedicated network
#'   attribute (currently \code{\%n\% "dyads"}), used when building each bipartite network from a partition.
#' @param group_labels Optional character vector (or NULL) forwarded to \code{build_bipartite_from_inputs()} to
#'   set readable group vertex labels in each padded bipartite network.
#'
#' @return A list of length \eqn{T}. If \code{eval.call=TRUE}, each element is a fitted \pkg{ergm} model for
#'   time \eqn{t}. If \code{eval.call=FALSE}, each element is the unevaluated \code{ergm()} call produced for
#'   time \eqn{t}.
#'
#' @examples
#' \dontrun{
#'   parts <- list(
#'     c(1, 1, 2, 2, 3),
#'     c(1, 2, 2, 3, 3),
#'     c(1, 1, 1, 2, 3)
#'   )
#'
#'   # Dry-run: inspect per-time translated calls and inertial activation
#'   calls <- erpm_long(parts ~ groups + inertia_term(past_influence = 2), eval.call = FALSE)
#'
#'   # Fit sequentially with shared RHS and inertial activation when eligible
#'   fits <- erpm_long(parts ~ groups + inertia_term(past_influence = 2), estimate = "MLE")
#' }
#'
#' @note
#' Inertial terms are detected and handled via \code{R/erpm_long_inertia_specs.R}. Each inertial term can
#' request \code{past_influence = d} lags, which are materialized as \eqn{d} separate network attributes on
#' the current network (one per lag), without combining information across lags.
#'
#' @note
#' This wrapper delegates model construction, translation, and constraint handling to \code{erpm()} at each
#' time point, ensuring a consistent call structure and shared safeguards across static and longitudinal use.
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

  .dbg <- function(...) {
    # Local debug printer.
    # Kept minimal so we can pass it down to the engine without pulling global state.
    if (isTRUE(debug)) cat("[ERPM_LONG][debug] ", ..., "\n", sep = "")
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
  if (!is.null(dyads) && is.list(dyads) && .erpm_long_is_list_of_dyads(dyads)) {
    # ok: list of lists (per time), e.g. dyads[[t]][["X"]] is a matrix.
  } else if (!is.null(dyads) && is.list(dyads)) {
    # ok: a single dyads list shared across time.
  } else if (!is.null(dyads)) {
    stop("[ERPM_LONG] `dyads` must be a list (single) or a list of lists (per time).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 2) Detect inertial calls once (RHS shared, activation depends on t)
  # ---------------------------------------------------------------------------
  # The RHS expression is the same for all t, but inertial terms can be inactive
  # at early times depending on how many past partitions they require.
  det <- .erpm_long_detect_inertial_calls(rhs_expr, env0, debug = debug)
  static_terms     <- det$static_terms
  inertial_calls <- det$inertial_calls

  # ---------------------------------------------------------------------------
  # 3) Verbose header
  # ---------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    # Human-friendly header so the user can see what will happen before the loop runs.
    cat("============================================================\n")
    cat("[ERPM_LONG] Start\n")
    cat(sprintf("[ERPM_LONG] T = %d partitions\n", T))
    cat(sprintf("[ERPM_LONG] eval.call = %s | estimate = %s | eval.loglik = %s\n",
                if (isTRUE(eval.call)) "TRUE" else "FALSE",
                if (is.null(estimate)) "NULL" else as.character(estimate),
                if (is.null(eval.loglik)) "NULL" else as.character(eval.loglik)))

    # Print the shared RHS once (useful for debugging formula parsing).
    cat("[ERPM_LONG] RHS (shared): ", .erpm_long_rhs_oneline(rhs_expr), "\n", sep = "")
    if (isTRUE(det$enabled)) {
      cat("[ERPM_LONG] inertial terms detected: ", paste(det$inertial_names, collapse = ", "), "\n", sep = "")
      # Important behavioral point: we do not pool information across lags, each lag is explicit.
      cat("[ERPM_LONG] note: activation is per-term via past_influence; lags are not reduced.\n")
    } else {
      cat("[ERPM_LONG] inertial terms detected: none\n")
    }
    cat("------------------------------------------------------------\n")
  }

  # ---------------------------------------------------------------------------
  # 4) Run
  # ---------------------------------------------------------------------------
  # Delegate the full time loop to the engine so `erpm_long()` stays a small
  # orchestration wrapper (parse, validate, then run).
  .erpm_long_run(
    parts         = parts,
    T             = T,
    env0          = env0,
    rhs_expr      = rhs_expr,
    static_terms    = static_terms,
    inertial_calls= inertial_calls,
    det           = det,
    eval.call     = eval.call,
    verbose       = verbose,
    debug         = debug,
    .dbg          = .dbg,
    estimate      = estimate,
    eval.loglik   = eval.loglik,
    control       = control,
    timeout       = timeout,
    nodes         = nodes,
    dyads         = dyads,
    group_labels  = group_labels
  )
}