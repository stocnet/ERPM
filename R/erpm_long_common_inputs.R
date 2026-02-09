# ################################################################################
# # FILE: R/erpm_long_common_inputs.R
# ################################################################################
# #' ERPM longitudinal common input helpers (internal)
# #'
# #' @name erpm_long_common_inputs
# #' @note erpm_long_common_inputs.R
# #'
# #' @description
# #' Shared helpers for time-indexed inputs in \code{erpm_long()}:
# #' \itemize{
# #'   \item dyads normalization;
# #'   \item per-time extraction for (nodes, dyads, group_labels);
# #'   \item compact formatting helpers.
# #' }
# #'
# #' @keywords ERPM ERGM longitudinal internal helpers
# NULL

# # ==============================================================================
# # Predicates / structure checks
# # ==============================================================================

# #' Check whether an object is a list of matrices (non-nested)
# #' @noRd
# .erpm_long_is_list_of_matrices <- function(x) {
#   is.list(x) &&
#     length(x) &&
#     all(vapply(x, is.matrix, logical(1)))
# }

# #' Check whether an object is a list of dyadic inputs (list of lists)
# #' @noRd
# .erpm_long_is_list_of_dyads <- function(x) {
#   is.list(x) &&
#     length(x) &&
#     all(vapply(x, is.list, logical(1)))
# }

# #' Check whether an object is a *named* list of matrices (non-nested)
# #' @noRd
# .erpm_long_is_named_list_of_matrices <- function(x) {
#   is.list(x) &&
#     length(x) &&
#     !is.null(names(x)) &&
#     all(nzchar(names(x))) &&
#     all(vapply(x, is.matrix, logical(1)))
# }

# # ==============================================================================
# # RHS parsing helpers (dyads)
# # ==============================================================================

# #' Extract all dyadic covariate names referenced in RHS (dyadcov*, cov_fullmatch*)
# #'
# #' Notes:
# #' - We treat dyad terms as those whose first argument is a character scalar naming
# #'   a dyadic matrix attached to the network.
# #' - This is used to decide whether a bare matrix dyads input can be auto-named.
# #'
# #' @noRd
# .erpm_long_extract_dyad_names <- function(rhs_expr) {
#   out <- character(0)

#   # Terms that reference a dyadic matrix name as their first argument (string).
#   dyad_terms <- c(
#     "dyadcov",
#     "dyadcov_full",
#     "dyadcov_GW",
#     "cov_fullmatch"
#   )

#   walk <- function(e) {
#     if (is.null(e)) return(invisible(NULL))

#     if (is.call(e)) {
#       fn <- e[[1L]]
#       if (is.symbol(fn) && as.character(fn) %in% dyad_terms) {
#         if (length(e) >= 2L &&
#             is.character(e[[2L]]) &&
#             length(e[[2L]]) == 1L &&
#             nzchar(e[[2L]])) {
#           out <<- c(out, e[[2L]])
#         }
#       }
#       for (i in seq_along(e)) walk(e[[i]])
#       return(invisible(NULL))
#     }

#     if (is.pairlist(e) || is.list(e)) {
#       for (i in seq_along(e)) walk(e[[i]])
#     }

#     invisible(NULL)
#   }

#   walk(rhs_expr)
#   unique(out)
# }

# # ==============================================================================
# # Dyads normalization
# # ==============================================================================

# #' Normalize dyads input for erpm_long()
# #'
# #' Accepted user inputs:
# #' - NULL or list(): no dyads
# #' - matrix: if RHS has exactly one dyad name => list(<name>=M);
# #'           if RHS has no dyad term => list("__erpm_unused_dyads__"=M);
# #'           if RHS has multiple dyad names => error (ambiguous).
# #' - named list of matrices: used as-is (shared across time)
# #' - list of length 1 containing a single matrix (unnamed): treated like matrix (rule above)
# #' - list of length T of matrices: per-time dyads, each wrapped as list(<name>=M) (rule above)
# #' - list of length T of (matrices OR lists of matrices): per-time dyads; each element
# #'   is normalized to a named list of matrices (possibly empty).
# #'
# #' Output invariant (critical):
# #' - Either a named list of matrices (shared),
# #' - Or a list(T) of named lists of matrices (time-indexed).
# #'
# #' @noRd
# .erpm_long_normalize_dyads_input <- function(dyads, rhs_expr, T, mode) {
#   stopifnot(mode %in% c("sequentiel", "empile"))
#   if (is.null(dyads) || (is.list(dyads) && length(dyads) == 0L)) return(NULL)
#   # if (is.list(dyads) && length(dyads) == 0L) return(list())

#   needed   <- .erpm_long_extract_dyad_names(rhs_expr)
#   n_needed <- length(needed)

#   # Decide how to name a bare dyads matrix when possible.
#   infer_single_name_or_error <- function() {
#     if (n_needed == 1L) return(needed[[1L]])
#     if (n_needed == 0L) return(NULL)
#     stop(
#       "[ERPM_LONG] `dyads` was provided as a matrix (or unnamed single matrix), but the RHS references multiple dyadic term names.\n",
#       "  Fix: pass `dyads` as a named list (shared) or a list(T) of named lists.\n",
#       call. = FALSE
#     )
#   }

#   # Placeholder naming when RHS does not consume dyads:
#   # users may pass dyads "just in case" without breaking.
#   placeholder_name <- function(i = NULL) {
#     if (is.null(i)) "__erpm_unused_dyads__" else paste0("__erpm_unused_dyads__", i, "__")
#   }

#   # ---------------------------------------------------------------------------
#   # Case 1: dyads is a bare matrix
#   # ---------------------------------------------------------------------------
#   if (is.matrix(dyads)) {
#     if (n_needed == 0L) return(setNames(list(dyads), placeholder_name()))
#     nm <- infer_single_name_or_error()
#     return(setNames(list(dyads), nm))
#   }

#   if (!is.list(dyads)) {
#     stop("[ERPM_LONG] `dyads` must be NULL, a matrix, a list of matrices, or a list of such lists (per time).", call. = FALSE)
#   }

#   # ---------------------------------------------------------------------------
#   # Case 2: dyads is a list of matrices (shared OR time-indexed)
#   # ---------------------------------------------------------------------------
#   if (.erpm_long_is_list_of_matrices(dyads)) {

#     unnamed <- is.null(names(dyads)) || any(!nzchar(names(dyads)))

#     # Unnamed singleton list(list(M)) -> treat like matrix.
#     if (unnamed && length(dyads) == 1L) {
#       if (n_needed == 0L) return(setNames(list(dyads[[1L]]), placeholder_name()))
#       nm <- infer_single_name_or_error()
#       return(setNames(list(dyads[[1L]]), nm))
#     }

#     # Per-time list(T) of matrices (unnamed): infer name if exactly one dyad name; else placeholders if none.
#     if (unnamed && length(dyads) == T) {
#       if (n_needed == 0L) {
#         return(lapply(dyads, function(M) setNames(list(M), placeholder_name())))
#       }
#       nm <- infer_single_name_or_error()
#       return(lapply(dyads, function(M) setNames(list(M), nm)))
#     }

#     # Shared named list of matrices: keep as shared.
#     if (.erpm_long_is_named_list_of_matrices(dyads)) {
#       return(dyads)
#     }

#     # Unnamed multi-matrix list:
#     # - If RHS needs none: accept and auto-name placeholders.
#     # - Else: require explicit names.
#     if (n_needed == 0L) {
#       out <- dyads
#       names(out) <- vapply(seq_along(out), function(i) placeholder_name(i), character(1))
#       return(out)
#     }

#     stop(
#       "[ERPM_LONG] `dyads` is a list of matrices but is not a valid named list.\n",
#       "  Fix: pass a named list like `list(block_att = M, mix_att = M2)`,\n",
#       "  or pass a single matrix / singleton list when the RHS contains exactly one dyad name.",
#       call. = FALSE
#     )
#   }

#   # ---------------------------------------------------------------------------
#   # Case 3: dyads is a list-of-lists (per-time)
#   # ---------------------------------------------------------------------------
#   if (.erpm_long_is_list_of_dyads(dyads)) {
#     if (length(dyads) != T) {
#       stop(sprintf("[ERPM_LONG] `dyads` as a list-of-lists must have length T=%d.", T), call. = FALSE)
#     }

#     out <- vector("list", T)

#     for (t in seq_len(T)) {
#       dt <- dyads[[t]]

#       if (is.null(dt) || (is.list(dt) && length(dt) == 0L)) {
#         out[[t]] <- list()
#         next
#       }

#       # dyads[[t]] provided as a bare matrix
#       if (is.matrix(dt)) {
#         if (n_needed == 0L) {
#           out[[t]] <- setNames(list(dt), placeholder_name())
#         } else {
#           nm <- infer_single_name_or_error()
#           out[[t]] <- setNames(list(dt), nm)
#         }
#         next
#       }

#       if (!is.list(dt)) {
#         stop(sprintf("[ERPM_LONG] dyads[[%d]] must be a matrix, a list of matrices, or empty.", t), call. = FALSE)
#       }

#       if (!all(vapply(dt, is.matrix, logical(1)))) {
#         stop(sprintf("[ERPM_LONG] dyads[[%d]] must contain only matrices.", t), call. = FALSE)
#       }

#       unnamed_t <- is.null(names(dt)) || any(!nzchar(names(dt)))

#       if (unnamed_t) {
#         # Unnamed singleton -> infer (or placeholder if RHS needs none)
#         if (length(dt) == 1L) {
#           if (n_needed == 0L) {
#             out[[t]] <- setNames(list(dt[[1L]]), placeholder_name())
#           } else {
#             nm <- infer_single_name_or_error()
#             out[[t]] <- setNames(list(dt[[1L]]), nm)
#           }
#         } else {
#           # Unnamed multi-list:
#           if (n_needed == 0L) {
#             tmp <- dt
#             names(tmp) <- vapply(seq_along(tmp), function(i) placeholder_name(i), character(1))
#             out[[t]] <- tmp
#           } else {
#             stop(sprintf("[ERPM_LONG] dyads[[%d]] must be a named list of matrices.", t), call. = FALSE)
#           }
#         }
#       } else {
#         # out[[t]] <- dt
#       }
#     }

#     # if (identical(mode, "empile")) {

#     #   # Shared dyads: must be matrices with identical dims/order across all times
#     #   if (is.list(dyads) && length(dyads)==T && all(vapply(dyads,is.list,TRUE))) {

#     #       stop(
#     #         "[ERPM_LONG|PLE] Shared dyads are not allowed unless explicitly time-indexed.\n",
#     #         "Fix: pass dyads as a list(T) of named lists of matrices.",
#     #         call. = FALSE
#     #       )
#     #   }
#     # }

#     return(out)
#   }

#   stop(
#     "[ERPM_LONG] `dyads` structure not recognized.\n",
#     "  Accepted: matrix, named list of matrices, list(T) of matrices, or list(T) of (matrices / named lists of matrices).",
#     call. = FALSE
#   )

# }


# # ==============================================================================
# # Nodes normalization (list-of-vectors -> data.frame)
# # ==============================================================================

# #' Normalize nodes input for erpm_long(): data.frame or named list of vectors
# #' @noRd
# # .erpm_long_normalize_nodes_input <- function(nodes, parts, T) {
# #   if (is.null(nodes)) return(NULL)

# #   # Helper: named list of vectors -> data.frame
# #   as_nodes_df <- function(x, n, fallback_labels) {
# #     if (is.data.frame(x)) return(x)

# #     if (!is.list(x)) {
# #       stop("[ERPM_LONG] `nodes` must be a data.frame, a named list of vectors, or a list(T) of those.", call. = FALSE)
# #     }

# #     if (length(x) == 0L) {
# #       # allow empty -> will become label-only df
# #       return(data.frame(label = fallback_labels, stringsAsFactors = FALSE))
# #     }

# #     if (is.null(names(x)) || any(!nzchar(names(x)))) {
# #       stop("[ERPM_LONG] `nodes` as a list must be a *named* list (e.g., list(colors=..., shapes=...)).", call. = FALSE)
# #     }

# #     if (!all(vapply(x, is.atomic, logical(1)))) {
# #       stop("[ERPM_LONG] `nodes` list must contain only atomic vectors.", call. = FALSE)
# #     }

# #     lens <- vapply(x, length, integer(1))
# #     if (any(lens != n)) {
# #       stop(sprintf("[ERPM_LONG] `nodes` list elements must all have length n=%d.", n), call. = FALSE)
# #     }

# #     df <- as.data.frame(x, stringsAsFactors = FALSE)

# #     # Ensure there is a label column (reuse existing if present)
# #     labcol <- .erpm_get_label_col(df)
# #     if (is.null(labcol)) {
# #       df$label <- fallback_labels
# #     } else {
# #       # standardize `label` column name
# #       df$label <- as.character(df[[labcol]])
# #     }

# #     df
# #   }

# #   # Shared nodes
# #   if (is.data.frame(nodes)) return(nodes)

# #   # Time-indexed nodes: list(T)
# #   if (is.list(nodes) && length(nodes) == T && !is.data.frame(nodes)) {
# #     n0 <- length(parts[[1L]])
# #     out <- vector("list", T)
# #     for (t in seq_len(T)) {
# #       n_t <- length(parts[[t]])
# #       if (n_t != n0) {
# #         stop("[ERPM_LONG] partitions must have constant length across time to normalize nodes.", call. = FALSE)
# #       }
# #       fallback_labels <- sprintf("A%d", seq_len(n_t))
# #       if (is.null(nodes[[t]])) {
# #         out[[t]] <- NULL
# #       } else {
# #         out[[t]] <- as_nodes_df(nodes[[t]], n = n_t, fallback_labels = fallback_labels)
# #       }
# #     }
# #     return(out)
# #   }

# #   # Shared list-of-vectors
# #   n <- length(parts[[1L]])
# #   fallback_labels <- sprintf("A%d", seq_len(n))
# #   as_nodes_df(nodes, n = n, fallback_labels = fallback_labels)
# # }
# .erpm_long_normalize_nodes_input <- function(nodes, parts, T, mode) {

#   if (is.null(nodes)) return(NULL)

#   # ... garde ton code existant au-dessus ...

#   as_nodes_df <- function(x, n, fallback_labels) {
#     if (is.data.frame(x)) return(x)

#     if (!is.list(x)) {
#       stop("[ERPM_LONG] `nodes` must be a data.frame, a named list of vectors, or a list(T) of those.",
#            call. = FALSE)
#     }

#     if (length(x) == 0L) {
#       return(data.frame(label = fallback_labels, stringsAsFactors = FALSE))
#     }

#     if (is.null(names(x)) || any(!nzchar(names(x)))) {
#       stop("[ERPM_LONG] `nodes` as a list must be a *named* list (e.g., list(colors=..., shapes=...)).",
#            call. = FALSE)
#     }

#     if (!all(vapply(x, is.atomic, logical(1)))) {
#       stop("[ERPM_LONG] `nodes` list must contain only atomic vectors.",
#            call. = FALSE)
#     }

#     lens <- vapply(x, length, integer(1))
#     if (any(lens != n)) {
#       stop(sprintf("[ERPM_LONG] `nodes` list elements must all have length n=%d.", n),
#            call. = FALSE)
#     }

#     df <- as.data.frame(x, stringsAsFactors = FALSE)
#     labcol <- .erpm_get_label_col(df)
#     if (is.null(labcol)) {
#       df$label <- fallback_labels
#     } else {
#       df$label <- as.character(df[[labcol]])
#     }
#     df
#   }

#   add_label_if_missing <- function(df, fallback_labels) {
#     if (is.null(df)) return(NULL)
#     if (!is.data.frame(df)) {
#       stop("[ERPM_LONG] `nodes` actors/groups must be data.frames or NULL.", call. = FALSE)
#     }
#     labcol <- .erpm_get_label_col(df)
#     if (is.null(labcol)) {
#       df$label <- fallback_labels
#     } else {
#       df$label <- as.character(df[[labcol]])
#     }
#     df
#   }

#   # if (is.data.frame(nodes)) return(nodes)

#   if (identical(mode, "empile") && is.data.frame(nodes)) {
#     if (!("id" %in% names(nodes))) stop( "[ERPM_LONG|PLE] Shared `nodes` must contain an actor identifier column (e.g. `id`).", call. = FALSE)
#     return(nodes)
#   }

#   if (is.list(nodes) && length(nodes) == T && !is.data.frame(nodes)) {
#     n0 <- length(parts[[1L]])
#     out <- vector("list", T)

#     for (t in seq_len(T)) {
#       n_t <- length(parts[[t]])
#       if (n_t != n0) {
#         stop("[ERPM_LONG] partitions must have constant length across time to normalize nodes.",
#              call. = FALSE)
#       }

#       fallback_actor_labels <- sprintf("A%d", seq_len(n_t))

#       if (is.null(nodes[[t]])) {
#         out[[t]] <- NULL
#         next
#       }

#       nt <- nodes[[t]]

#       if (is.list(nt) && !is.data.frame(nt) &&
#           !is.null(names(nt)) && any(names(nt) %in% c("actors", "groups"))) {

#         bad_names <- setdiff(names(nt), c("actors", "groups"))
#         if (length(bad_names) > 0L) {
#           stop(sprintf(
#             "[ERPM_LONG] nodes[[%d]] list may only contain 'actors'/'groups' (found: %s).",
#             t, paste(bad_names, collapse = ", ")
#           ), call. = FALSE)
#         }

#         actors_df <- add_label_if_missing(nt$actors, fallback_actor_labels)

#         # fallback group labels: based on provided df rows if any, otherwise empty
#         groups_df <- NULL
#         if (!is.null(nt$groups)) {
#           fallback_group_labels <- sprintf("G%d", seq_len(nrow(nt$groups)))
#           groups_df <- add_label_if_missing(nt$groups, fallback_group_labels)
#         }

#         out[[t]] <- list(actors = actors_df, groups = groups_df)
#         next
#       }

#       # Existing behavior: data.frame OR named list-of-vectors
#       out[[t]] <- as_nodes_df(nt, n = n_t, fallback_labels = fallback_actor_labels)
#     }

#     return(out)
#   }

#   # ... garde le reste de ton code existant (cas shared list-of-vectors, etc.) ...
#   as_nodes_df(nodes, n = length(parts[[1L]]), fallback_labels = sprintf("A%d", seq_len(length(parts[[1L]]))))
# }


# #' Extract the time-t element of a possibly time-indexed input
# #'
# #' IMPORTANT (dyads invariant):
# #' After normalization, `dyads` is either:
# #'   - a named list of matrices (shared across time), OR
# #'   - a list of length T of named lists of matrices (time-indexed).
# #' Therefore, `get_t(..., what="dyads")` must only choose between these two cases.
# #'
# #' @noRd
# .erpm_long_get_t <- function(x, t, T, what = "input") {
#   if (is.null(x)) return(NULL)
#   if (is.list(x) && length(x) == 0L) return(NULL)

#   # Non-list inputs (or data.frames) are treated as shared objects.
#   if (!is.list(x) || is.data.frame(x)) return(x)

#   # Dyads: rely on the post-normalization invariant.
#   if (identical(what, "dyads")) {
#     # Time-indexed dyads: list(T) of dyad-lists
#     if (length(x) == T && all(vapply(x, is.list, logical(1)))) {
#       return(x[[t]])
#     }
#     # Shared dyads: named list of matrices (or empty handled above)
#     return(x)
#   }

#   # Generic time-indexed case (nodes, group_labels, etc.)
#   if (length(x) != T) {
#     stop(sprintf(
#       "[ERPM_LONG] `%s` must be NULL, a single object, or a list of length T=%d.",
#       what, T
#     ), call. = FALSE)
#   }

#   x[[t]]
# }

# # ==============================================================================
# # Compact formatting helpers
# # ==============================================================================

# #' Compact formatter for named integer tables
# #' @noRd
# .erpm_long_tab_str <- function(x) {
#   paste(sprintf("%s:%s", names(x), as.integer(x)), collapse = " ")
# }

# #' Build a human-readable summary of a partition
# #' @noRd
# .erpm_long_partition_info <- function(p) {
#   gid <- as.integer(p)
#   sizes <- sort(table(gid))

#   sprintf(
#     "n=%d | groups=%d | sizes={%s}",
#     length(gid),
#     length(unique(gid)),
#     .erpm_long_tab_str(sizes)
#   )
# }