################################################################################
# FILE: R/erpm_validate_inputs.R
################################################################################
#' ERPM input validation: nodes and dyads
#' @name erpm_validate_inputs
#' @note erpm_validate_inputs.R
#'
#' @description
#' This file provides internal validators used by the ERPM wrapper:
#' \itemize{
#'   \item label column detection for node data frames;
#'   \item validation of node attribute tables;
#'   \item validation of dyadic n×n matrices.
#' }
#'
#' These helpers are used by the bipartite builder.
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
  if (!is.null(lab) && anyDuplicated(nodes[[lab]]))
    stop(sprintf("nodes$%s contains duplicates.", lab))
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
      if (!identical(rownames(M), labels) || !identical(colnames(M), labels))
        stop(sprintf("dyads['%s']: row/colnames must match the actor order.", nm))
    }
  }
  invisible(TRUE)
}