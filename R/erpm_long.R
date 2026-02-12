################################################################################
# FILE: R/erpm_long.R
################################################################################
#' ERPM pseudo-longitudinal wrapper (PLE-only)
#'
#' @name erpm_long
#' @note erpm_long.R
#'
#' @description
#' This file provides \code{erpm_long()}, a pragmatic pseudo-longitudinal wrapper
#' that operates in PLE ("empile") mode only.
#'
#' The user supplies:
#' \itemize{
#'   \item a timeline of partitions on the LHS (a list of length \eqn{T \ge 2});
#'   \item a shared RHS ERPM/ERGM specification (applied to all estimated blocks).
#' }
#'
#' The wrapper then:
#' \enumerate{
#'   \item validates inputs and derives settings (notably inertial presence and the
#'         effective maximum \code{past_influence});
#'   \item builds a stacked bipartite meta-network (PLE engine);
#'   \item delegates call construction and optional evaluation to \code{erpm()}.
#' }
#'
#' In PLE, the first \code{d} time points are not estimated when inertial terms with
#' depth \code{d} are present, since those blocks do not have enough past context.
#'
#' @keywords ERPM ERGM longitudinal wrapper PLE empile
################################################################################

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
#' @examples
#' \dontrun{
#'   partitions <- list(c(1,1,2), c(1,2,2), c(2,2,1))
#'   fit <- erpm_long(partitions ~ cliques(k=2), eval.call=FALSE)
#' }
#'
#' @keywords ERPM ERGM longitudinal PLE
#' @export
erpm_long <- function(formula,
                      mode = "empile",
                      eval.call = TRUE,
                      verbose = FALSE,
                      debug = FALSE,
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

  if (isTRUE(verbose)) {
    tlabs <- attr(terms(rhs), "term.labels")
    .erpm_long_vcat(verbose, sprintf("[ERPM_LONG] mode=PLE | T=%d | eval.call=%s", length(partitions), as.character(isTRUE(eval.call))))
    .erpm_long_vcat(verbose, sprintf("[ERPM_LONG] RHS=%s", if (length(tlabs)) paste(tlabs, collapse = " + ") else "<empty>"))
    .erpm_long_vcat(verbose, sprintf("[ERPM_LONG] inertial_present=%s | past_influence(d)=%d",
                                     as.character(inertial_present), if (inertial_present) d else 0L))
  }

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
    verbose          = verbose,
    debug            = debug
  )
  meta_nw <- built$meta_nw

  # ---------------------------------------------------------------------------
  # 3) Compose erpm() call on the meta-network
  # ---------------------------------------------------------------------------
  call_erpm <- as.call(c(
    list(as.name("erpm"), as.formula(call("~", meta_nw, rhs))),
    if (!is.null(eval.call))    list(eval.call = eval.call)     else list(),
    if (!is.null(verbose))      list(verbose = verbose)         else list(),
    if (isTRUE(debug))          list(debug = TRUE)              else list(),
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
  attr(out, "erpm_long.selected_partition_indices") <- built$selected_partition_indices
  attr(out, "erpm_long.d") <- if (inertial_present) d else 0

  out
}