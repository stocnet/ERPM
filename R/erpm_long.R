################################################################################
# FILE: R/erpm_long.R
################################################################################

#' ERPM pseudo-longitudinal wrapper (PLE-only)
#'
#' @name erpm_long
#' @note erpm_long.R
#'
#' @description
#' This module implements a PLE-only (stacked/empile) pseudo-longitudinal wrapper.
#' The user provides a list of partitions on the LHS (length \eqn{T}) and a shared
#' RHS specification. The wrapper:
#' \enumerate{
#'   \item validates inputs and resolves derived settings (e.g., inertial presence,
#'         maximum \code{past_influence});
#'   \item builds a stacked bipartite meta-network (engine);
#'   \item delegates fitting (or dry-run call generation) to \code{erpm()}.
#' }
#'
#' The meta-network is constructed so that all covariates needed by RHS terms are
#' attached prior to calling \code{erpm()}.
#'
#' @keywords ERPM ERGM longitudinal wrapper PLE

# ==============================================================================
# Internal helpers (PLE)
# ==============================================================================

#' Normalize dyadic inputs for PLE timeline mode (internal helper)
#'
#' @param dyads NULL or a list of length T, each element a named list (name -> matrix).
#' @param T Number of partitions (timeline length).
#'
#' @return NULL or a list with fields:
#' \itemize{
#'   \item \code{mode}: always \code{"timeline"}
#'   \item \code{dyads}: named list; for each name, a list of length T of matrices
#' }
#'
#' @noRd
.erpm_ple_normalize_dyads <- function(dyads, T) {
  if (is.null(dyads)) return(NULL)

  # Expected: dyads = list(t=1..T) of list(name -> matrix)
  if (!is.list(dyads) || length(dyads) != T) {
    stop("[ERPM_PLE] dyads must be a list of length T (timeline mode).")
  }

  # Union of dyad names across time
  nms <- unique(unlist(lapply(dyads, names)))
  if (!length(nms)) stop("[ERPM_PLE] dyads timeline has no named attributes.")

  # Build dyads_data[[name]][[t]] = matrix
  dyads_data <- setNames(vector("list", length(nms)), nms)
  for (nm in nms) {
    mats <- vector("list", T)
    for (t in seq_len(T)) {
      Mt <- dyads[[t]][[nm]]
      if (is.null(Mt)) stop(sprintf("[ERPM_PLE] dyads[%d] missing '%s'.", t, nm))
      if (!is.matrix(Mt)) stop(sprintf("[ERPM_PLE] dyads[%d]$%s must be a matrix.", t, nm))
      mats[[t]] <- Mt
    }
    dyads_data[[nm]] <- mats
  }

  list(mode = "timeline", dyads = dyads_data)
}

#' Detect inertial terms on RHS and extract max past influence (internal helper)
#'
#' @description
#' Conservative detection: only \code{inertia_groups(...)} is recognized here.
#' The maximum \code{past_influence} across occurrences is returned.
#'
#' @param rhs RHS expression from the user formula.
#'
#' @return A list with fields:
#' \itemize{
#'   \item \code{inertial_present} (logical)
#'   \item \code{d} (integer): max past influence (0 if absent)
#' }
#'
#' @noRd
.erpm_long_detect_inertial <- function(rhs) {
  tt <- terms(rhs)
  labels <- attr(tt, "term.labels")

  if (is.null(labels) || !length(labels)) {
    return(list(inertial_present = FALSE, d = 0L))
  }

  idx <- grep("^inertia_groups\\b", labels)
  if (!length(idx)) {
    return(list(inertial_present = FALSE, d = 0L))
  }

  dmax <- 1L

  for (lab in labels[idx]) {
    expr <- try(parse(text = lab)[[1L]], silent = TRUE)
    if (inherits(expr, "try-error") || !is.call(expr)) next
    if (!identical(as.character(expr[[1L]]), "inertia_groups")) next

    args <- as.list(expr)[-1L]

    if (!length(args)) {
      d_i <- 1L
    } else if ("past_influence" %in% names(args)) {
      d_i <- suppressWarnings(as.integer(round(eval(args[["past_influence"]], parent.frame()))))
      if (is.na(d_i)) d_i <- 1L
    } else if ("pi" %in% names(args)) {
      d_i <- suppressWarnings(as.integer(round(eval(args[["pi"]], parent.frame()))))
      if (is.na(d_i)) d_i <- 1L
    } else if ("d" %in% names(args)) {
      d_i <- suppressWarnings(as.integer(round(eval(args[["d"]], parent.frame()))))
      if (is.na(d_i)) d_i <- 1L
    } else {
      d_i <- suppressWarnings(as.integer(round(eval(args[[1L]], parent.frame()))))
      if (is.na(d_i)) d_i <- 1L
    }

    if (d_i < 1L) stop("[ERPM_LONG] inertia_groups: past_influence must be >= 1.")
    if (d_i > dmax) dmax <- d_i
  }

  list(inertial_present = TRUE, d = dmax)
}

#' Validate erpm_long inputs and resolve derived settings (internal helper)
#'
#' @description
#' Centralizes all validation and resolution needed before building the meta-network:
#' \itemize{
#'   \item formula shape and mode (PLE-only);
#'   \item LHS evaluation (list of partitions);
#'   \item inertial detection and \code{past_influence} bounds.
#' }
#'
#' @param formula User formula: \code{partitions_list ~ rhs}.
#' @param mode Only \code{"empile"} / \code{"PLE"} supported.
#' @param eval.call Kept for signature compatibility (not used here).
#' @param verbose Kept for signature compatibility (not used here).
#' @param debug Kept for signature compatibility (not used here).
#' @param estimate Kept for signature compatibility (not used here).
#' @param eval.loglik Kept for signature compatibility (not used here).
#' @param control Kept for signature compatibility (not used here).
#' @param timeout Kept for signature compatibility (not used here).
#' @param seed Kept for signature compatibility (not used here).
#' @param nodes Kept for signature compatibility (not used here).
#' @param dyads Kept for signature compatibility (not used here).
#' @param group_labels Kept for signature compatibility (not used here).
#'
#' @return A list of resolved objects (partitions, rhs, inertial flags, etc.).
#'
#' @noRd
.erpm_long_validate_inputs <- function(formula, mode, eval.call, verbose, debug,
                                      estimate, eval.loglik, control, timeout, seed,
                                      nodes, dyads, group_labels) {
  if (!inherits(formula, "formula")) stop("[ERPM_LONG] formula must be a formula.")
  mode <- match.arg(mode)

  if (!identical(mode, "empile") && !identical(mode, "PLE")) {
    stop("[ERPM_LONG] Only PLE/empile mode is supported (PLS removed; V2 only).")
  }

  # LHS/RHS extraction
  lhs <- formula[[2L]]
  rhs <- formula[[3L]]

  if (is.null(lhs) || is.null(rhs)) {
    stop("[ERPM_LONG] formula must be of the form: partitions_list ~ terms")
  }

  # Evaluate partitions from the calling environment
  partitions <- eval(lhs, envir = parent.frame())
  if (!is.list(partitions) || !length(partitions)) {
    stop("[ERPM_LONG] LHS must evaluate to a non-empty list of partitions.")
  }
  T <- length(partitions)

  # Inertial detection (currently: inertia_groups)
  inert <- .erpm_long_detect_inertial(rhs)
  inertial_present <- isTRUE(inert$inertial_present)
  d <- as.integer(inert$d)

  # past_influence must leave at least one estimable block
  if (inertial_present && d >= T) {
    stop(sprintf("[ERPM_LONG] past_influence=%d but T=%d: need d <= T-1.", d, T))
  }

  list(
    lhs             = lhs,
    rhs             = rhs,
    partitions      = partitions,
    T               = T,
    inertial_present = inertial_present,
    past_influence  = d
  )
}

# ==============================================================================
# Public API
# ==============================================================================

#' Pseudo-longitudinal ERPM wrapper (PLE-only)
#'
#' @param formula ERGM/ERPM formula. LHS must evaluate to a list of partitions
#'   (length \eqn{T}). RHS is a shared model specification.
#' @param mode Kept for compatibility; only \code{"empile"} / \code{"PLE"} supported.
#' @param eval.call If TRUE, returns the \code{erpm()} call instead of evaluating.
#' @param verbose Verbosity.
#' @param debug Reserved for compatibility (passed through to \code{erpm()}).
#' @param estimate Passed through to \code{erpm()}.
#' @param eval.loglik Passed through to \code{erpm()}.
#' @param control Passed through to \code{erpm()}.
#' @param timeout Passed through to \code{erpm()}.
#' @param seed Passed through to \code{erpm()}.
#' @param nodes NULL or list length T of data.frames of monadic covariates.
#' @param dyads NULL or dyadic covariates (engine-dependent format).
#' @param group_labels Optional group labels (passed to the engine when relevant).
#'
#' @return Result of \code{erpm()} (or the call if \code{eval.call=TRUE}).
#'   The built meta-network is attached as \code{attr(out, "meta_nw")} when evaluated.
#'
#' @export
erpm_long <- function(formula,
                      mode = c("empile", "PLE"),
                      eval.call = TRUE,
                      verbose = FALSE,
                      debug = NULL,
                      estimate = NULL,
                      eval.loglik = NULL,
                      control = NULL,
                      timeout = NULL,
                      seed = NULL,
                      nodes = NULL,
                      dyads = NULL,
                      group_labels = NULL) {

  # ---------------------------------------------------------------------------
  # 1) Validate inputs + resolve derived settings
  # ---------------------------------------------------------------------------
  validated <- .erpm_long_validate_inputs(
    formula       = formula,
    mode          = mode,
    eval.call     = eval.call,
    verbose       = verbose,
    debug         = debug,
    estimate      = estimate,
    eval.loglik   = eval.loglik,
    control       = control,
    timeout       = timeout,
    seed          = seed,
    nodes         = nodes,
    dyads         = dyads,
    group_labels  = group_labels
  )

  rhs              <- validated$rhs
  partitions       <- validated$partitions
  inertial_present <- validated$inertial_present
  d                <- validated$past_influence

  # ---------------------------------------------------------------------------
  # 2) Build meta-network (engine)
  # ---------------------------------------------------------------------------
  built <- .erpm_long_empile_build_meta_nw(
    partitions       = partitions,
    rhs              = rhs,
    inertial_present = inertial_present,
    past_influence   = d,
    nodes            = nodes,
    dyads            = dyads,
    group_labels     = group_labels,
    directed         = FALSE,
    verbose          = verbose
  )
  meta_nw <- built$meta_nw

  # ---------------------------------------------------------------------------
  # 3) Compose erpm() call on the meta-network
  # ---------------------------------------------------------------------------
  call_erpm <- as.call(c(
    list(as.name("erpm"), as.formula(call("~", meta_nw, rhs))),
    if (!is.null(eval.call))    list(eval.call = eval.call)     else list(),
    if (!is.null(verbose))      list(verbose = verbose)         else list(),
    if (!is.null(debug))        list(debug = debug)             else list(),
    if (!is.null(estimate))     list(estimate = estimate)       else list(),
    if (!is.null(eval.loglik))  list(eval.loglik = eval.loglik) else list(),
    if (!is.null(control))      list(control = control)         else list(),
    if (!is.null(timeout))      list(timeout = timeout)         else list(),
    if (!is.null(seed))         list(seed = seed)               else list()
  ))

  if (isTRUE(eval.call)) return(call_erpm)

  out <- eval(call_erpm, envir = parent.frame())

  # Attach meta-network for debugging / introspection
  attr(out, "meta_nw") <- meta_nw
  attr(out, "erpm_long.idx_est") <- built$idx_est
  attr(out, "erpm_long.d") <- if (inertial_present) d else 0

  out
}