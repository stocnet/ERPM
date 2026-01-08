################################################################################
# FILE: R/erpm_long_inertia_detect.R
################################################################################
#' ERPM longitudinal inertia detection helper (internal)
#'
#' @name erpm_long_inertia_detect
#' @note erpm_long_inertia_detect.R
#'
#' @description
#' This file implements the detection logic for inertial (longitudinal) terms
#' appearing in the RHS of an \code{erpm_long()} formula.
#'
#' Its role is purely *syntactic and declarative*:
#' \enumerate{
#'   \item split the shared RHS expression into individual ERGM terms;
#'   \item identify which of these terms correspond to registered inertial effects
#'         (as declared in \code{.erpm_long_inertia_registry});
#'   \item separate inertial terms from standard (static) ERGM terms.
#' }
#'
#' IMPORTANT:
#' \itemize{
#'   \item This file does \emph{not} decide when inertial terms are active in time.
#'         Temporal activation is handled later, based on \code{past_influence}.
#'   \item No model fitting or network manipulation happens here.
#'   \item The output is a lightweight "plan" describing how the RHS should be
#'         interpreted by the longitudinal engine.
#' }
#'
#' @keywords ERPM ERGM longitudinal inertia internal
NULL

# =============================================================================
# Helpers
# =============================================================================

#' Detect inertial terms in a shared RHS expression
#'
#' This internal helper inspects the RHS of an \code{erpm_long()} formula and
#' determines which terms are inertial (longitudinal) and which are standard
#' static ERGM terms.
#'
#' Conceptually, it performs a *single pass* over the RHS and builds a small
#' descriptor used later by the longitudinal engine:
#' \itemize{
#'   \item inertial terms are identified by name, using
#'         \code{.erpm_long_inertia_registry};
#'   \item non-inertial terms are kept unchanged and will be included at every
#'         time step;
#'   \item inertial terms are preserved as calls so their arguments
#'         (e.g. \code{past_influence}) can be evaluated later.
#' }
#'
#' @param rhs_expr An unevaluated RHS expression (as obtained from formula parsing).
#' @param env0     The environment in which the formula was defined (used only
#'                 for debugging messages at this stage).
#' @param debug    Logical; if TRUE, emit diagnostic messages.
#'
#' @return A list with four components:
#' \itemize{
#'   \item \code{enabled}: TRUE if at least one inertial term is present;
#'   \item \code{inertial_calls}: list of inertial term calls;
#'   \item \code{inertial_names}: character vector of inertial term names;
#'   \item \code{static_terms}: list of non-inertial RHS terms.
#' }
#'
#' @noRd
.erpm_long_detect_inertial_calls <- function(rhs_expr, env0, debug = FALSE) {

  # ---------------------------------------------------------------------------
  # 1) Split RHS into individual terms
  # ---------------------------------------------------------------------------
  # The RHS may be a nested sum (a + b + c). We recursively flatten it into
  # a simple list of terms so that each effect can be inspected independently.
  terms <- .erpm_long_split_sum_terms(rhs_expr)

  # Extract the function/term name of each RHS item.
  # For example:
  #   groups              -> "groups"
  #   inertia_groups(...) -> "inertia_groups"
  fnames <- vapply(terms, .erpm_long_term_name, character(1))

  # Drop malformed or unidentifiable terms (defensive programming).
  fnames <- fnames[!is.na(fnames)]

  # ---------------------------------------------------------------------------
  # 2) Identify inertial terms by registry lookup
  # ---------------------------------------------------------------------------
  # An inertial term is defined purely by its name being present in the
  # inertia registry. No argument evaluation is done here.
  inertial_names <- intersect(unique(fnames),
                              names(.erpm_long_inertia_registry))

  # Flag indicating whether *any* inertial effect is present at all.
  enabled <- length(inertial_names) > 0L

  # Optional debug output for developers.
  .erpm_long_dbg(
    debug,
    "[ERPM_LONG|DEBUG] detect inertial:",
    paste(inertial_names, collapse = ", ")
  )

  # ---------------------------------------------------------------------------
  # 3) Fast exit if no inertial terms are present
  # ---------------------------------------------------------------------------
  # In the common case where the RHS is purely static, we simply return all
  # terms as base terms and avoid any additional processing downstream.
  if (!enabled) {
    return(list(
      enabled        = FALSE,
      inertial_calls = list(),
      inertial_names = character(),
      static_terms     = terms
    ))
  }

  # ---------------------------------------------------------------------------
  # 4) Separate inertial calls from base (static) terms
  # ---------------------------------------------------------------------------
  inertial_calls <- list()
  static_terms     <- list()

  for (tt in terms) {

    # Determine the name of the current term.
    nm <- .erpm_long_term_name(tt)

    if (nm %in% inertial_names) {
      # This term is inertial.
      # If it was written as a bare symbol (e.g. inertia_groups),
      # normalize it into a call object for uniform downstream handling.
      if (is.symbol(tt)) tt <- as.call(list(tt))

      inertial_calls[[length(inertial_calls) + 1L]] <- tt
    } else {
      # This term is static and will be included at every time step.
      static_terms[[length(static_terms) + 1L]] <- tt
    }
  }

  # ---------------------------------------------------------------------------
  # 5) Return inertia detection summary
  # ---------------------------------------------------------------------------
  list(
    enabled        = TRUE,
    inertial_calls = inertial_calls,
    inertial_names = inertial_names,
    static_terms     = static_terms
  )
}