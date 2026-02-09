################################################################################
# FILE: R/erpm_long_validate.R
# OBJECT: Validation utilities for PLE-only erpm_long()
# NOTES :
#   - Centralises all user-facing validation.
#   - No PLS logic.
#   - Fails early and explicitly when meta-network construction is impossible.
################################################################################

# ------------------------------------------------------------------------------
# Validate partitions input
# ------------------------------------------------------------------------------
.erpm_long_validate_partitions <- function(partitions) {
  if (!is.list(partitions) || !length(partitions)) {
    stop("[ERPM_LONG] partitions must be a non-empty list.")
  }

  nA_ref <- length(partitions[[1L]])
  if (nA_ref == 0L) {
    stop("[ERPM_LONG] partitions[[1]] is empty.")
  }

  for (t in seq_along(partitions)) {
    p <- partitions[[t]]
    if (!is.vector(p) || is.list(p)) {
      stop(sprintf("[ERPM_LONG] partitions[[%d]] must be an atomic vector.", t))
    }
    if (length(p) != nA_ref) {
      stop(sprintf(
        "[ERPM_LONG] partitions[[%d]] length (%d) != partitions[[1]] length (%d).",
        t, length(p), nA_ref
      ))
    }
    if (any(is.na(p))) {
      stop(sprintf("[ERPM_LONG] partitions[[%d]] contains NA values.", t))
    }
  }

  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Validate nodes (monadic covariates)
# ------------------------------------------------------------------------------
.erpm_long_validate_nodes <- function(nodes, partitions) {
  if (is.null(nodes)) return(invisible(TRUE))

  T <- length(partitions)
  if (!is.list(nodes) || length(nodes) != T) {
    stop("[ERPM_LONG] nodes must be NULL or a list of length T (one data.frame per partition).")
  }

  nA <- length(partitions[[1L]])

  ref_cols <- NULL
  ref_label <- NULL

  for (t in seq_len(T)) {
    df <- nodes[[t]]
    if (!is.data.frame(df)) stop(sprintf("[ERPM_LONG] nodes[[%d]] must be a data.frame.", t))
    if (nrow(df) != nA) {
      stop(sprintf("[ERPM_LONG] nodes[[%d]] has %d rows but partitions have %d actors.", t, nrow(df), nA))
    }
    if (!("label" %in% colnames(df))) {
      stop(sprintf("[ERPM_LONG] nodes[[%d]] must contain a 'label' column.", t))
    }
    cov_names <- setdiff(colnames(df), "label")
    if (!length(cov_names)) {
      stop(sprintf("[ERPM_LONG] nodes[[%d]] must have at least one covariate column besides 'label'.", t))
    }

    if (is.null(ref_cols)) {
      ref_cols <- colnames(df)
      ref_label <- as.character(df[["label"]])
      if (anyNA(ref_label) || any(ref_label == "")) stop("[ERPM_LONG] nodes[[1]]$label contains NA/empty values.")
    } else {
      if (!identical(colnames(df), ref_cols)) {
        stop(sprintf("[ERPM_LONG] nodes[[%d]] column names differ from nodes[[1]].", t))
      }
      lab <- as.character(df[["label"]])
      if (!identical(lab, ref_label)) {
        stop(sprintf("[ERPM_LONG] nodes[[%d]]$label differs from nodes[[1]]$label. Actor ordering must be stable across time.", t))
      }
    }
  }

  invisible(TRUE)
}
# .erpm_long_validate_nodes <- function(nodes, partitions) {
#   if (is.null(nodes)) return(invisible(TRUE))

#   T <- length(partitions)

#   if (!is.list(nodes) || length(nodes) != T) {
#     stop("[ERPM_LONG] nodes must be NULL or a list of length T (one data.frame per partition).")
#   }

#   nA <- length(partitions[[1L]])

#   ref_names <- NULL
#   for (t in seq_len(T)) {
#     df <- nodes[[t]]
#     if (!is.data.frame(df)) {
#       stop(sprintf("[ERPM_LONG] nodes[[%d]] must be a data.frame.", t))
#     }
#     if (nrow(df) != nA) {
#       stop(sprintf(
#         "[ERPM_LONG] nodes[[%d]] has %d rows but partitions have %d actors.",
#         t, nrow(df), nA
#       ))
#     }

#     if (is.null(ref_names)) {
#       ref_names <- colnames(df)
#       if (!length(ref_names)) {
#         stop("[ERPM_LONG] nodes data.frames must have at least one column.")
#       }
#     } else {
#       if (!identical(colnames(df), ref_names)) {
#         stop(sprintf(
#           "[ERPM_LONG] nodes[[%d]] column names differ from nodes[[1]].",
#           t
#         ))
#       }
#     }
#   }

#   invisible(TRUE)
# }

# ------------------------------------------------------------------------------
# Validate dyads (dyadic covariates)
# ------------------------------------------------------------------------------
.erpm_long_validate_dyads <- function(dyads, partitions) {
  if (is.null(dyads)) return(invisible(TRUE))

  T <- length(partitions)
  nA <- length(partitions[[1L]])

  check_mat <- function(M, label, t) {
    if (!is.matrix(M)) {
      stop(sprintf("[ERPM_LONG] dyads '%s' at t=%d is not a matrix.", label, t))
    }
    if (nrow(M) != nA || ncol(M) != nA) {
      stop(sprintf(
        "[ERPM_LONG] dyads '%s' at t=%d has dim %dx%d; expected %dx%d.",
        label, t, nrow(M), ncol(M), nA, nA
      ))
    }
  }

  # Single matrix (ambiguous except T=1)
  if (is.matrix(dyads)) {
    if (T != 1L) {
      stop("[ERPM_LONG] dyads given as a matrix but T>1. Use a list of length T.")
    }
    check_mat(dyads, "Z", 1L)
    return(invisible(TRUE))
  }

  if (!is.list(dyads)) {
    stop("[ERPM_LONG] dyads must be NULL, a matrix, or a list.")
  }

  # Unnamed list of matrices: interpreted as single attribute over time
  if (is.null(names(dyads))) {
    if (length(dyads) != T) {
      stop("[ERPM_LONG] dyads list must have length T.")
    }
    for (t in seq_len(T)) {
      check_mat(dyads[[t]], "Z", t)
    }
    return(invisible(TRUE))
  }

  # Named list: one attribute per name
  for (nm in names(dyads)) {
    x <- dyads[[nm]]
    if (is.matrix(x)) {
      if (T != 1L) {
        stop(sprintf(
          "[ERPM_LONG] dyads[['%s']] is a matrix but T>1. Use list of length T.",
          nm
        ))
      }
      check_mat(x, nm, 1L)
    } else if (is.list(x)) {
      if (length(x) != T) {
        stop(sprintf(
          "[ERPM_LONG] dyads[['%s']] must be a list of length T.", nm
        ))
      }
      for (t in seq_len(T)) {
        check_mat(x[[t]], nm, t)
      }
    } else {
      stop(sprintf(
        "[ERPM_LONG] dyads[['%s']] must be a matrix or a list of matrices.", nm
      ))
    }
  }

  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Validate past_influence feasibility
# ------------------------------------------------------------------------------
.erpm_long_validate_past_influence <- function(T, inertial_present, d) {
  if (!inertial_present) return(invisible(TRUE))
  if (d < 0L) {
    stop("[ERPM_LONG] past_influence must be >= 0.")
  }
  if (d >= T) {
    stop(sprintf(
      "[ERPM_LONG] past_influence=%d but T=%d: need d <= T-1.",
      d, T
    ))
  }
  invisible(TRUE)
}

# ################################################################################
# # FILE: R/erpm_long_validators.R
# ################################################################################
# #' ERPM longitudinal validators (internal)
# #'
# #' @name erpm_long_validators
# #' @note erpm_long_validators.R
# #'
# #' @description
# #' Small validators used by \code{erpm_long()} for:
# #' \itemize{
# #'   \item reproducibility seed normalization;
# #'   \item time-indexed nodes input integrity.
# #' }
# #'
# #' @keywords ERPM ERGM longitudinal internal helpers
# NULL

# # ============================================================================
# # Seed validator (internal)
# # ============================================================================

# #' Validate and normalize a seed for reproducible evaluation
# #'
# #' Ensures compatibility with \code{set.seed()} and returns a normalized integer.
# #'
# #' @noRd
# .erpm_long_validate_seed <- function(seed) {
#   if (is.null(seed)) return(NULL)

#   if (!(is.numeric(seed) && length(seed) == 1L && is.finite(seed))) {
#     stop("[ERPM_LONG] `seed` must be a single finite numeric value (integer-like) or NULL.", call. = FALSE)
#   }

#   si <- as.integer(round(seed))
#   if (!isTRUE(all.equal(seed, si))) {
#     stop("[ERPM_LONG] `seed` must be integer-valued (e.g., 1, 2, 42).", call. = FALSE)
#   }

#   # Base R expects a non-negative integer in practice.
#   # (R stores RNG seed as integer vector; negative is not a valid user seed here.)
#   if (si < 0L) {
#     stop("[ERPM_LONG] `seed` must be >= 0.", call. = FALSE)
#   }

#   # Be explicit about the usual integer range used by R RNG.
#   # This guards against accidental double seeds > 2^31-1.
#   if (si > .Machine$integer.max) {
#     stop(sprintf("[ERPM_LONG] `seed` must be <= %d.", .Machine$integer.max), call. = FALSE)
#   }

#   si
# }

# # ============================================================================
# # Nodes validator (internal)
# # ============================================================================

# #' Validate nodes input for erpm_long()
# #'
# #' Accepted:
# #' - NULL
# #' - one data.frame (shared across time)
# #' - list of length T of data.frame (time-indexed)
# #'
# #' This validator checks structure only; size consistency vs partitions is checked
# #' at engine level (PLS/PLE) when the partition length is known for a given t.
# #'
# #' @noRd
# # .erpm_long_validate_nodes_input <- function(nodes, T = NULL) {

# #   if (is.null(nodes)) return(invisible(TRUE))

# #   # Shared nodes: single data.frame
# #   if (is.data.frame(nodes)) {
# #     ok <- try(.erpm_check_nodes_df(nodes), silent = TRUE)
# #     if (inherits(ok, "try-error")) {
# #       stop(
# #         paste0("[ERPM_LONG] Invalid `nodes` data.frame: ",
# #                conditionMessage(attr(ok, "condition"))),
# #         call. = FALSE
# #       )
# #     }
# #     return(invisible(TRUE))
# #   }

# #   if (is.list(nodes) && !is.data.frame(nodes)) {
# #     if (length(nodes) && (is.null(names(nodes)) || any(!nzchar(names(nodes))))) {
# #         stop("[ERPM_LONG] `nodes` as list must be named (e.g., list(colors=..., shapes=...)).", call. = FALSE)
# #       }
# #       if (length(nodes) && !all(vapply(nodes, is.atomic, logical(1)))) {
# #         stop("[ERPM_LONG] `nodes` list must contain only atomic vectors.", call. = FALSE)
# #       }
# #       return(invisible(TRUE))
# #   }

# #   # Time-indexed nodes: list(T) of data.frame
# #   if (!is.list(nodes)) {
# #     stop("[ERPM_LONG] `nodes` must be NULL, a data.frame, or a list of data.frames.", call. = FALSE)
# #   }

# #   if (!is.null(T)) {
# #     if (!(length(nodes) == T)) {
# #       stop(sprintf("[ERPM_LONG] `nodes` as a list must have length T=%d.", T), call. = FALSE)
# #     }
# #   }

# #   for (t in seq_along(nodes)) {
# #     nt <- nodes[[t]]
# #     if (is.null(nt)) next  # allow explicit NULL per time point
# #     if (is.list(nt) && !is.data.frame(nt)) {
# #       if (length(nt) && (is.null(names(nt)) || any(!nzchar(names(nt))))) {
# #         stop(sprintf("[ERPM_LONG] nodes[[%d]] as list must be named (e.g., list(colors=..., shapes=...)).", t), call. = FALSE)
# #       }
# #       if (length(nt) && !all(vapply(nt, is.atomic, logical(1)))) {
# #         stop(sprintf("[ERPM_LONG] nodes[[%d]] list must contain only atomic vectors.", t), call. = FALSE)
# #       }
# #       next
# #     }
# #     ok <- try(.erpm_check_nodes_df(nt), silent = TRUE)
# #     if (inherits(ok, "try-error")) {
# #       stop(
# #         paste0("[ERPM_LONG] Invalid nodes[[", t, "]] data.frame: ",
# #                conditionMessage(attr(ok, "condition"))),
# #         call. = FALSE
# #       )
# #     }
# #   }

# #   invisible(TRUE)
# # }
# .erpm_long_validate_nodes_input <- function(nodes, T = NULL) {

#   if (is.null(nodes)) return(invisible(TRUE))

#   # Shared nodes: single data.frame
#   if (is.data.frame(nodes)) {
#     ok <- try(.erpm_check_nodes_df(nodes), silent = TRUE)
#     if (inherits(ok, "try-error")) {
#       stop(
#         paste0("[ERPM_LONG] Invalid `nodes` data.frame: ",
#                conditionMessage(attr(ok, "condition"))),
#         call. = FALSE
#       )
#     }
#     return(invisible(TRUE))
#   }

#   # From here: must be a list (either shared list-of-vectors OR list(T) time-indexed)
#   if (!is.list(nodes)) {
#     stop("[ERPM_LONG] `nodes` must be NULL, a data.frame, a named list of vectors, or a list of such objects (per time).", call. = FALSE)
#   }

#   # -------------------------------------------------------------------------
#   # Shared nodes: named list of vectors (e.g., list(colors=..., shapes=...))
#   #
#   # IMPORTANT:
#   # Do NOT misclassify a time-indexed list(T) as a shared attribute list.
#   # If T is known and length(nodes)==T, we treat it as time-indexed.
#   # -------------------------------------------------------------------------
#   is_time_indexed <- !is.null(T) && length(nodes) == T

#   if (!is_time_indexed) {
#     # Accept empty list() as a degenerate shared nodes spec (will become label-only later)
#     if (length(nodes) == 0L) return(invisible(TRUE))

#     if (is.null(names(nodes)) || any(!nzchar(names(nodes)))) {
#       stop("[ERPM_LONG] `nodes` as list must be named (e.g., list(colors=..., shapes=...)).", call. = FALSE)
#     }
#     if (!all(vapply(nodes, is.atomic, logical(1)))) {
#       stop("[ERPM_LONG] `nodes` list must contain only atomic vectors.", call. = FALSE)
#     }
#     return(invisible(TRUE))
#   }

#   # -------------------------------------------------------------------------
#   # Time-indexed nodes: list(T) of data.frame OR named list-of-vectors OR NULL
#   # -------------------------------------------------------------------------
#     for (t in seq_along(nodes)) {
#     nt <- nodes[[t]]
#     if (is.null(nt)) next  # allow explicit NULL per time point

#     # Allow per-time list inputs:
#     #   - named list of atomic vectors (colors/shapes/...)
#     #   - OR list(actors=<df>, groups=<df>) where each df is validated
#     if (is.list(nt) && !is.data.frame(nt)) {

#       # Case A: list(actors=df, groups=df) (or a subset)
#       if (!is.null(names(nt)) && any(names(nt) %in% c("actors", "groups"))) {
#         bad_names <- setdiff(names(nt), c("actors", "groups"))
#         if (length(bad_names) > 0L) {
#           stop(sprintf(
#             "[ERPM_LONG] nodes[[%d]] list may only contain 'actors'/'groups' when using data.frame mode (found: %s).",
#             t, paste(bad_names, collapse = ", ")
#           ), call. = FALSE)
#         }

#         if (!is.null(nt$actors)) {
#           if (!is.data.frame(nt$actors)) {
#             stop(sprintf("[ERPM_LONG] nodes[[%d]]$actors must be a data.frame or NULL.", t), call. = FALSE)
#           }
#           ok <- try(.erpm_check_nodes_df(nt$actors), silent = TRUE)
#           if (inherits(ok, "try-error")) {
#             stop(
#               paste0("[ERPM_LONG] Invalid nodes[[", t, "]]$actors data.frame: ",
#                      conditionMessage(attr(ok, "condition"))),
#               call. = FALSE
#             )
#           }
#         }

#         if (!is.null(nt$groups)) {
#           if (!is.data.frame(nt$groups)) {
#             stop(sprintf("[ERPM_LONG] nodes[[%d]]$groups must be a data.frame or NULL.", t), call. = FALSE)
#           }
#           ok <- try(.erpm_check_nodes_df(nt$groups), silent = TRUE)
#           if (inherits(ok, "try-error")) {
#             stop(
#               paste0("[ERPM_LONG] Invalid nodes[[", t, "]]$groups data.frame: ",
#                      conditionMessage(attr(ok, "condition"))),
#               call. = FALSE
#             )
#           }
#         }

#         next
#       }

#       # Case B: named list of atomic vectors (colors/shapes/...)
#       if (length(nt) && (is.null(names(nt)) || any(!nzchar(names(nt))))) {
#         stop(sprintf("[ERPM_LONG] nodes[[%d]] as list must be named (e.g., list(colors=..., shapes=...)).", t), call. = FALSE)
#       }
#       if (length(nt) && !all(vapply(nt, is.atomic, logical(1)))) {
#         stop(sprintf("[ERPM_LONG] nodes[[%d]] list must contain only atomic vectors.", t), call. = FALSE)
#       }
#       next
#     }

#     # data.frame case
#     ok <- try(.erpm_check_nodes_df(nt), silent = TRUE)
#     if (inherits(ok, "try-error")) {
#       stop(
#         paste0("[ERPM_LONG] Invalid nodes[[", t, "]] data.frame: ",
#                conditionMessage(attr(ok, "condition"))),
#         call. = FALSE
#       )
#     }
#   }
#   # for (t in seq_along(nodes)) {
#   #   nt <- nodes[[t]]
#   #   if (is.null(nt)) next  # allow explicit NULL per time point

#   #   # Allow per-time named list-of-vectors
#   #   if (is.list(nt) && !is.data.frame(nt)) {
#   #     if (length(nt) && (is.null(names(nt)) || any(!nzchar(names(nt))))) {
#   #       stop(sprintf("[ERPM_LONG] nodes[[%d]] as list must be named (e.g., list(colors=..., shapes=...)).", t), call. = FALSE)
#   #     }
#   #     if (length(nt) && !all(vapply(nt, is.atomic, logical(1)))) {
#   #       stop(sprintf("[ERPM_LONG] nodes[[%d]] list must contain only atomic vectors.", t), call. = FALSE)
#   #     }
#   #     next
#   #   }

#   #   # data.frame case
#   #   if (!is.data.frame(nt)) {
#   #     stop(sprintf("[ERPM_LONG] nodes[[%d]] must be a data.frame, a named list of vectors, or NULL.", t), call. = FALSE)
#   #   }

#   #   ok <- try(.erpm_check_nodes_df(nt), silent = TRUE)
#   #   if (inherits(ok, "try-error")) {
#   #     stop(
#   #       paste0("[ERPM_LONG] Invalid nodes[[", t, "]] data.frame: ",
#   #              conditionMessage(attr(ok, "condition"))),
#   #       call. = FALSE
#   #     )
#   #   }
#   # }

#   invisible(TRUE)
# }