################################################################################
# FILE: R/erpm_validate_inputs.R
################################################################################
#' ERPM input validation: nodes, dyads, and optional MH mix arguments
#' @name erpm_validate_inputs
#' @note erpm_validate_inputs.R
#'
#' @description
#' This file provides internal validators used by the ERPM wrapper:
#' \itemize{
#'   \item label column detection for node data frames;
#'   \item validation of node attribute tables;
#'   \item validation of dyadic n×n matrices;
#'   \item validation of optional \code{erpm()} MH mix arguments
#'         (\code{mh_moves}, \code{mh_weights}).
#' }
#'
#' These helpers are used by the bipartite builder and by the main wrapper.
#'
#' @keywords ERPM ERGM wrapper validation

# ============================================================================
# Node / dyad validation helpers
# ============================================================================

#' Guess a label column in a node data frame
#' @noRd
.erpm_get_label_col <- function(nodes, prefer = "label") {
  stopifnot(is.data.frame(nodes))
  nms <- trimws(names(nodes))
  aliases <- c(prefer, "label", "nom", "name", "id")
  # Avoid duplicate "label" entries when prefer = "label"
  hit <- intersect(aliases, nms)
  if (length(hit)) hit[1L] else NULL
}

#' Validate a node data frame
#' @noRd
.erpm_check_nodes_df <- function(nodes) {
  stopifnot(is.data.frame(nodes))
  lab <- .erpm_get_label_col(nodes)  # may be NULL if no obvious label column
  if (!is.null(lab) && anyDuplicated(nodes[[lab]])) {
    stop(sprintf("nodes$%s contains duplicates.", lab))
  }
  invisible(TRUE)
}

#' Validate dyadic n×n matrices
#' @noRd
.erpm_check_dyads <- function(dyads, n, labels) {
  if (length(dyads) == 0L) return(invisible(TRUE))
  stopifnot(is.list(dyads))

  if (is.null(names(dyads)) || any(!nzchar(names(dyads)))) {
    stop("dyads must be a *named* list of n×n matrices (e.g., list(X = M)).")
  }

  for (nm in names(dyads)) {
    M <- dyads[[nm]]
    stopifnot(is.matrix(M), nrow(M) == n, ncol(M) == n)
    if (!is.null(rownames(M)) && !is.null(colnames(M))) {
      if (!identical(rownames(M), labels) || !identical(colnames(M), labels)) {
        stop(sprintf("dyads['%s']: row/colnames must match the actor order.", nm))
      }
    }
  }
  invisible(TRUE)
}

# ============================================================================
# Optional MH mix validation for erpm(...)
# ============================================================================

#' Validate optional mh_moves / mh_weights passed to erpm()
#'
#' @details
#' This validator keeps the wrapper behavior explicit:
#' - if both arguments are NULL, nothing is injected and ergm keeps its default
#'   proposal logic;
#' - if one is provided without the other, we stop early because the intent is
#'   incomplete;
#' - if both are provided, they must already form a clean, aligned specification
#'   that can be forwarded to ErpmMix.
#'
#' The R-side initializer of ErpmMix remains permissive by design, but at the
#' erpm() entry point we keep validation stricter so user mistakes are caught
#' immediately.
#'
#' @param mh_moves Optional character vector of move names.
#' @param mh_weights Optional numeric vector of move weights.
#' @return A normalized list with elements \code{mh_moves} and \code{mh_weights}.
#' @noRd
.erpm_validate_mh_mix_inputs <- function(mh_moves = NULL, mh_weights = NULL) {
  both_null <- is.null(mh_moves) && is.null(mh_weights)
  if (both_null) {
    return(list(mh_moves = NULL, mh_weights = NULL))
  }

  if (is.null(mh_moves) || is.null(mh_weights)) {
    stop(
      "[ERPM] `mh_moves` and `mh_weights` must be provided together or both left NULL.",
      call. = FALSE
    )
  }

  allowed_moves <- c("toggle", "swap", "merge", "split")

  # --- Validate mh_moves on its own ------------------------------------------
  if (!is.character(mh_moves)) {
    stop(
      "[ERPM] `mh_moves` must be a character vector like c(\"toggle\", \"swap\").",
      call. = FALSE
    )
  }

  mh_moves <- tolower(trimws(mh_moves))
  if (length(mh_moves) < 1L || any(!nzchar(mh_moves))) {
    stop(
      "[ERPM] `mh_moves` must contain at least one non-empty move name.",
      call. = FALSE
    )
  }

  bad_moves <- setdiff(unique(mh_moves), allowed_moves)
  if (length(bad_moves) > 0L) {
    stop(
      sprintf(
        "[ERPM] `mh_moves` contains unsupported move(s): %s. Allowed values are: %s.",
        paste(bad_moves, collapse = ", "),
        paste(allowed_moves, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Keep the wrapper strict and predictable: duplicated names are usually a user
  # mistake, not a meaningful specification.
  if (anyDuplicated(mh_moves)) {
    stop(
      "[ERPM] `mh_moves` must not contain duplicates.",
      call. = FALSE
    )
  }

  # --- Validate mh_weights on its own ----------------------------------------
  if (!is.numeric(mh_weights)) {
    stop(
      "[ERPM] `mh_weights` must be a numeric vector like c(2, 1).",
      call. = FALSE
    )
  }

  if (length(mh_weights) < 1L) {
    stop(
      "[ERPM] `mh_weights` must contain at least one value.",
      call. = FALSE
    )
  }

  if (any(!is.finite(mh_weights)) || any(mh_weights <= 0)) {
    stop(
      "[ERPM] `mh_weights` must contain only finite positive values.",
      call. = FALSE
    )
  }

  # --- Joint validation -------------------------------------------------------
  if (length(mh_moves) != length(mh_weights)) {
    stop(
      sprintf(
        "[ERPM] `mh_moves` and `mh_weights` must have the same length. Got %d and %d.",
        length(mh_moves), length(mh_weights)
      ),
      call. = FALSE
    )
  }

  list(
    mh_moves   = mh_moves,
    mh_weights = as.numeric(mh_weights)
  )
}

# ################################################################################
# # FILE: R/erpm_validate_inputs.R
# ################################################################################
# #' ERPM input validation: nodes and dyads
# #' @name erpm_validate_inputs
# #' @note erpm_validate_inputs.R
# #'
# #' @description
# #' This file provides internal validators used by the ERPM wrapper:
# #' \itemize{
# #'   \item label column detection for node data frames;
# #'   \item validation of node attribute tables;
# #'   \item validation of dyadic n×n matrices.
# #' }
# #'
# #' These helpers are used by the bipartite builder.
# #'
# #' @keywords ERPM ERGM wrapper validation

# # ============================================================================
# # Node / dyad validation helpers
# # ============================================================================

# #' Guess a label column in a node data frame
# #' @noRd
# .erpm_get_label_col <- function(nodes, prefer = "label") {
#   stopifnot(is.data.frame(nodes))
#   nms <- trimws(names(nodes))
#   aliases <- c(prefer, "label", "nom", "name", "id")
#   # Avoid duplicate "label" entries when prefer = "label"
#   hit <- intersect(aliases, nms)
#   if (length(hit)) hit[1L] else NULL
# }

# #' Validate a node data frame
# #' @noRd
# .erpm_check_nodes_df <- function(nodes) {
#   stopifnot(is.data.frame(nodes))
#   lab <- .erpm_get_label_col(nodes)  # may be NULL if no obvious label column
#   if (!is.null(lab) && anyDuplicated(nodes[[lab]]))
#     stop(sprintf("nodes$%s contains duplicates.", lab))
#   invisible(TRUE)
# }

# #' Validate dyadic n×n matrices
# #' @noRd
# .erpm_check_dyads <- function(dyads, n, labels) {
#   if (length(dyads) == 0L) return(invisible(TRUE))
#   stopifnot(is.list(dyads))

#   if (is.null(names(dyads)) || any(!nzchar(names(dyads)))) {
#     stop("dyads must be a *named* list of n×n matrices (e.g., list(X = M)).")
#   }

#   for (nm in names(dyads)) {
#     M <- dyads[[nm]]
#     stopifnot(is.matrix(M), nrow(M) == n, ncol(M) == n)
#     if (!is.null(rownames(M)) && !is.null(colnames(M))) {
#       if (!identical(rownames(M), labels) || !identical(colnames(M), labels))
#         stop(sprintf("dyads['%s']: row/colnames must match the actor order.", nm))
#     }
#   }
#   invisible(TRUE)
# }