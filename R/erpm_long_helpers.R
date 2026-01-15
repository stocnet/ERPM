# =====================================================================
# FILE: R/erpm_long_helpers.R
# =====================================================================
#' ERPM longitudinal helpers: small utilities for erpm_long() (internal)
#'
#' @name erpm_long_helpers
#' @note erpm_long_helpers.R
#'
#' @description
#' Internal helpers for erpm_long(): input normalization, indexing, and small formatters used by the longitudinal pipeline.
#'
#' These helpers are intentionally:
#' \itemize{
#'   \item stateless and side-effect free;
#'   \item focused on one narrow task (type checks, formatting, indexing);
#'   \item reusable across the longitudinal pipeline.
#' }
#'
#' They are separated from the main logic of \code{erpm_long()} to:
#' \itemize{
#'   \item improve readability of the main algorithm;
#'   \item reduce cognitive load when debugging or extending the code;
#'   \item make implicit assumptions explicit (e.g. expected input shapes).
#' }
#'
#' All functions are internal (non-exported) and prefixed with \code{.erpm_long_}.
#'
#' @keywords ERPM ERGM longitudinal internal helpers
NULL

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Check whether an object is a list of data.frames
#'
#' This is used to validate inputs such as `nodes` that may be provided either
#' as a single data.frame (recycled over time) or as a list of per-time data.frames.
#'
#' @param x Any R object.
#' @return Logical scalar. TRUE if `x` is a non-empty list and all elements
#'         are data.frames.
#' @noRd
.erpm_long_is_list_of_df <- function(x) {
  is.list(x) &&
    length(x) &&
    all(vapply(x, is.data.frame, logical(1)))
}

#' Check whether an object is a list of dyadic inputs
#'
#' In erpm_long(), dyads are expected to be supplied as:
#' \itemize{
#'   \item a single list of matrices (recycled over time), or
#'   \item a list of such lists, one per time point.
#' }
#'
#' This helper only checks the *outer* structure: list of lists.
#'
#' @param x Any R object.
#' @return Logical scalar. TRUE if `x` is a non-empty list whose elements
#'         are themselves lists.
#' @noRd
.erpm_long_is_list_of_dyads <- function(x) {
  is.list(x) &&
    length(x) &&
    all(vapply(x, is.list, logical(1)))
}

#' Check whether an object is a list of matrices (non-nested)
#' @noRd
.erpm_long_is_list_of_matrices <- function(x) {
  is.list(x) &&
    length(x) &&
    all(vapply(x, is.matrix, logical(1)))
}

#' Collect dyadic covariate names referenced on the RHS (string literal only)
#'
#' Conservative parser: only literal strings, no evaluation
#'
#' Recognized patterns:
#' - dyadcov_*("name") and dyadcov("name") (prefix "^dyadcov")
#' - cov_fullmatch("name") (known dyadic covariate term in ERPM/ERGM extensions)
#'
#' @param rhs_expr RHS expression (language object).
#' @return Character vector of referenced dyad names (unique).
#' @noRd
.erpm_long_rhs_dyad_names <- function(rhs_expr) {
  if (is.null(rhs_expr)) return(character(0))

  is_target <- function(cc) {
    if (!is.call(cc)) return(FALSE)
    if (!is.symbol(cc[[1L]])) return(FALSE)
    fn <- as.character(cc[[1L]])
    startsWith(fn, "dyadcov") || identical(fn, "cov_fullmatch")
  }

  # Local collector (duplicated from deep debug helpers to avoid tight coupling).
  out <- list()
  walk <- function(x) {
    if (is.call(x)) {
      if (isTRUE(is_target(x))) out[[length(out) + 1L]] <<- x
      for (i in seq_along(x)) walk(x[[i]])
    } else if (is.pairlist(x) || is.list(x)) {
      for (i in seq_along(x)) walk(x[[i]])
    }
    invisible(NULL)
  }
  walk(rhs_expr)

  if (!length(out)) return(character(0))

  nm <- character(0)
  for (cc in out) {
    if (length(cc) >= 2L) {
      a1 <- cc[[2L]]
      if (is.character(a1) && length(a1) == 1L && nzchar(a1)) nm <- c(nm, a1)
    }
  }

  unique(nm)
}

#' Infer a single dyad name from the RHS if possible
#' @noRd
.erpm_long_infer_single_dyad_name <- function(rhs_expr) {
  nm <- .erpm_long_rhs_dyad_names(rhs_expr)
  if (length(nm) == 1L) nm else NULL
}

#' Normalize dyads input for erpm_long()
#'
#' Accepted user inputs:
#' - NULL or list(): no dyads
#' - matrix: converted to list(<inferred_name> = M) if RHS has exactly one dyad name
#' - named list of matrices: used as-is (shared across time)
#' - list of length T of matrices: treated as per-time dyads, each wrapped as list(<name>=M)
#' - list of length T of lists of matrices: treated as per-time dyads, validated
#'
#' @param dyads User dyads input.
#' @param rhs_expr RHS expression used to infer a single dyad name when needed.
#' @param T Number of time points.
#' @return Either:
#'   - a single named list of matrices, OR
#'   - a list of length T, each element a named list of matrices
#' @noRd
.erpm_long_normalize_dyads_input <- function(dyads, rhs_expr, T) {

  # NULL or empty list -> empty dyads shared.
  if (is.null(dyads)) return(list())
  if (is.list(dyads) && length(dyads) == 0L) return(list())

  infer_name <- function() {
    nm <- .erpm_long_infer_single_dyad_name(rhs_expr)
    if (is.null(nm)) {
      stop(
        "[ERPM_LONG] `dyads` was provided as a matrix (or unnamed single matrix), but the RHS does not contain exactly one dyadic term name.\n",
        "  Expected something like: dyadcov_full(\"X\") with a unique X.\n",
        "  Fix: pass `dyads = list(X = M)` or ensure the RHS contains one unique dyad name.",
        call. = FALSE
      )
    }
    nm
  }

  # Case A: dyads is a matrix -> wrap with inferred name
  if (is.matrix(dyads)) {
    nm <- infer_name()
    return(setNames(list(dyads), nm))
  }

  # From here dyads must be a list-like object.
  if (!is.list(dyads)) {
    stop("[ERPM_LONG] `dyads` must be NULL, a matrix, a list of matrices, or a list of such lists (per time).", call. = FALSE)
  }

  # Case B: list of matrices (non-nested)
  if (.erpm_long_is_list_of_matrices(dyads)) {

    # If length == T, interpret as per-time single dyad matrix.
    if (length(dyads) == T) {
      nm <- infer_name()
      return(lapply(dyads, function(M) setNames(list(M), nm)))
    }

    # Otherwise interpret as shared dyads list.
    if (is.null(names(dyads)) || any(!nzchar(names(dyads)))) {
      if (length(dyads) == 1L) {
        nm <- infer_name()
        return(setNames(list(dyads[[1L]]), nm))
      }
      stop(
        "[ERPM_LONG] `dyads` is a list of matrices but is not named.\n",
        "  Fix: pass a named list like `list(block_att = M, mix_att = M2)`.",
        call. = FALSE
      )
    }

    return(dyads)
  }

  # Case C: list of lists (per-time dyads)
  if (.erpm_long_is_list_of_dyads(dyads)) {
    if (length(dyads) != T) {
      stop(sprintf("[ERPM_LONG] `dyads` as a list-of-lists must have length T=%d.", T), call. = FALSE)
    }

    out <- vector("list", T)
    for (t in seq_len(T)) {
      dt <- dyads[[t]]
      if (is.null(dt) || (is.list(dt) && length(dt) == 0L)) {
        out[[t]] <- list()
        next
      }
      if (is.matrix(dt)) {
        # Defensive: allow a matrix inside list-of-lists.
        nm <- infer_name()
        out[[t]] <- setNames(list(dt), nm)
        next
      }
      if (!is.list(dt)) {
        stop(sprintf("[ERPM_LONG] dyads[[%d]] must be a list of matrices (or empty).", t), call. = FALSE)
      }
      if (!all(vapply(dt, is.matrix, logical(1)))) {
        stop(sprintf("[ERPM_LONG] dyads[[%d]] must contain only matrices.", t), call. = FALSE)
      }
      if (is.null(names(dt)) || any(!nzchar(names(dt)))) {
        if (length(dt) == 1L) {
          nm <- infer_name()
          out[[t]] <- setNames(list(dt[[1L]]), nm)
        } else {
          stop(sprintf("[ERPM_LONG] dyads[[%d]] must be a named list of matrices.", t), call. = FALSE)
        }
      } else {
        out[[t]] <- dt
      }
    }
    return(out)
  }

  stop(
    "[ERPM_LONG] `dyads` structure not recognized.\n",
    "  Accepted: matrix, named list of matrices, list(T) of matrices, or list(T) of named lists of matrices.",
    call. = FALSE
  )
}

#' Extract the time-\eqn{t} element of a possibly time-indexed input
#'
#' Many inputs in erpm_long() (nodes, dyads, group_labels) can be provided as:
#' \itemize{
#'   \item NULL,
#'   \item a single object (used at all time points),
#'   \item a list of length T, providing one object per time point.
#' }
#'
#' IMPORTANT (dyads special-case):
#' `dyads` can be supplied as a *single named list of matrices* (shared across time).
#' That object is a list, but it is NOT time-indexed. We treat it as a single object.
#'
#' This helper centralizes the logic for safely extracting the correct
#' element at time `t`, and for emitting a clear error if the length is invalid.
#'
#' @param x Input object (NULL, single object, or list).
#' @param t Integer time index (1-based).
#' @param T Total number of time points.
#' @param what Character string used in error messages to identify the input.
#' @return The object to be used at time `t`, or NULL.
#' @noRd
.erpm_long_get_t <- function(x, t, T, what = "input") {
  # NULL stays NULL: no input provided.
  if (is.null(x)) return(NULL)

  if (is.list(x) && length(x) == 0L) return(NULL)

  # Non-list or data.frame inputs are treated as time-invariant
  # and reused at every t.
  if (!is.list(x) || is.data.frame(x)) return(x)

  # Special-case: `dyads` may be a single named list of matrices shared across time.
  # That object is a list, but it is not a time-indexed list(T).
  if (identical(what, "dyads") && .erpm_long_is_list_of_matrices(x)) {
    return(x)
  }

  # If a list is provided, it must match the number of time points.
  if (length(x) != T) {
    stop(sprintf(
      "[ERPM_LONG] `%s` must be NULL, a single object, or a list of length T=%d.",
      what, T
    ), call. = FALSE)
  }

  # Extract the element for time t.
  x[[t]]
}

#' Convert a RHS expression into a single-line character string
#'
#' Used purely for verbose/debug output, to make logs compact and readable.
#' Long expressions are deparsed with a large width cutoff and collapsed
#' onto a single line.
#'
#' @param expr A language object (typically a RHS expression).
#' @return Character string representation on one line.
#' @noRd
.erpm_long_rhs_oneline <- function(expr) {
  paste(deparse(expr, width.cutoff = 500L), collapse = " ")
}

#' Combine a list of ERGM terms into a single RHS expression
#'
#' Given a list of calls/symbols representing ERGM terms, this helper
#' reconstructs a sum-of-terms expression using `+`.
#'
#' Special cases:
#' \itemize{
#'   \item empty list   -> NULL (no RHS);
#'   \item length 1     -> the term itself (no wrapping);
#'   \item length > 1   -> Reduce with `+`.
#' }
#'
#' @param terms List of language objects (calls or symbols).
#' @return A language object suitable for use on the RHS of a formula, or NULL.
#' @noRd
.erpm_long_combine_terms <- function(terms) {
  if (!length(terms)) return(NULL)
  if (length(terms) == 1L) return(terms[[1L]])

  Reduce(function(x, y) call("+", x, y), terms)
}

#' Compact formatter for named integer tables
#'
#' This helper is used only for verbose output to summarize group size tables
#' in a compact, stable textual form (e.g. "1:3 2:5 3:2").
#'
#' @param x A named integer vector or table.
#' @return Character string.
#' @noRd
.erpm_long_tab_str <- function(x) {
  paste(sprintf("%s:%s", names(x), as.integer(x)), collapse = " ")
}

#' Build a human-readable summary of a partition
#'
#' This helper summarizes a partition vector by reporting:
#' \itemize{
#'   \item total number of actors;
#'   \item number of distinct groups;
#'   \item sorted group sizes.
#' }
#'
#' It is intended for logging/debugging only and has no effect on model logic.
#'
#' @param p Integer vector encoding a partition (group id per actor).
#' @return Character string describing the partition.
#' @noRd
.erpm_long_partition_info <- function(p) {
  gid <- as.integer(p)
  sizes <- sort(table(gid))

  sprintf(
    "n=%d | groups=%d | sizes={%s}",
    length(gid),
    length(unique(gid)),
    .erpm_long_tab_str(sizes)
  )
}

# ---------------------------------------------------------------------------
# Deep debug helpers (nodecov/nodefactor readiness diagnostics)
# ---------------------------------------------------------------------------

#' Flatten an expression tree to collect all calls matching a predicate
#' @noRd
.erpm_long_collect_calls <- function(expr, pred) {
  out <- list()

  walk <- function(x) {
    if (is.call(x)) {
      if (isTRUE(pred(x))) out[[length(out) + 1L]] <<- x
      for (i in seq_along(x)) walk(x[[i]])
    } else if (is.pairlist(x) || is.list(x)) {
      for (i in seq_along(x)) walk(x[[i]])
    }
    invisible(NULL)
  }

  walk(expr)
  out
}

#' Extract node attribute names referenced by nodecov()/nodefactor() in a RHS expression
#'
#' We keep this intentionally conservative:
#' - Only accept literal strings: nodecov("age") / nodefactor('gender').
#' - Ignore dynamic expressions so deep debug does not accidentally evaluate user code.
#' @noRd
.erpm_long_node_attr_names_from_rhs <- function(rhs_expr) {
  if (is.null(rhs_expr)) return(character(0))

  is_target <- function(cc) {
    is.call(cc) &&
      is.symbol(cc[[1L]]) &&
      as.character(cc[[1L]]) %in% c("nodecov", "nodefactor")
  }

  calls <- .erpm_long_collect_calls(rhs_expr, is_target)
  if (!length(calls)) return(character(0))

  attrs <- character(0)
  for (cc in calls) {
    # nodecov("age") / nodefactor("gender") => first argument after function name.
    if (length(cc) >= 2L) {
      a1 <- cc[[2L]]
      if (is.character(a1) && length(a1) == 1L) attrs <- c(attrs, a1)
    }
  }

  unique(attrs)
}

#' Print missingness diagnostics for a node attribute by mode (actors vs groups)
#' @noRd
.erpm_long_deep_dump_node_attr <- function(nw, attr_name, max_show = 6L) {
  n1 <- as.integer(network::get.network.attribute(nw, "bipartite"))
  n  <- network::network.size(nw)
  i_actor <- seq_len(n1)
  i_group <- if (n1 < n) (n1 + 1L):n else integer(0)

  v <- try(network::get.vertex.attribute(nw, attr_name), silent = TRUE)
  if (inherits(v, "try-error") || is.null(v)) {
    cat(sprintf("[ERPM_LONG][deep] node attr '%s': MISSING (no such vertex attribute)\n", attr_name))
    return(invisible(NULL))
  }

  # Guarantee length for safety in weird cases.
  if (length(v) < n) v <- c(v, rep(NA, n - length(v)))

  vA <- v[i_actor]
  vG <- if (length(i_group)) v[i_group] else NULL

  naA <- sum(is.na(vA))
  naG <- if (is.null(vG)) 0L else sum(is.na(vG))

  cat(sprintf("[ERPM_LONG][deep] node attr '%s':\n", attr_name))
  cat(sprintf("  - actors (1..%d): NA=%d/%d\n", n1, naA, length(vA)))
  if (length(i_group)) {
    cat(sprintf("  - groups (%d..%d): NA=%d/%d\n", n1 + 1L, n, naG, length(vG)))
  } else {
    cat("  - groups: none\n")
  }

  # Small value preview (non-NA), to catch type/pathologies quickly.
  showA <- head(vA[!is.na(vA)], max_show)
  showG <- if (is.null(vG)) NULL else head(vG[!is.na(vG)], max_show)

  if (length(showA)) {
    cat("  - actors sample (non-NA): ", paste(showA, collapse = ", "), "\n", sep = "")
  } else {
    cat("  - actors sample (non-NA): (none)\n")
  }

  if (length(i_group)) {
    if (length(showG)) {
      cat("  - groups sample (non-NA): ", paste(showG, collapse = ", "), "\n", sep = "")
    } else {
      cat("  - groups sample (non-NA): (none)\n")
    }
  }

  invisible(NULL)
}

#' Deep debug entrypoint: inspect which node attributes are required by nodecov/nodefactor
#'
#' This uses the translated RHS (ergm-side names) so it matches what ergm() will see.
#' @noRd
.erpm_long_deep_debug_nodecov_inputs <- function(nw, rhs_translated) {
  attrs <- .erpm_long_node_attr_names_from_rhs(rhs_translated)
  if (!length(attrs)) {
    cat("[ERPM_LONG][deep] no nodecov()/nodefactor() terms detected in translated RHS.\n")
    return(invisible(NULL))
  }

  cat("[ERPM_LONG][deep] nodecov/nodefactor attributes referenced: ",
      paste(attrs, collapse = ", "), "\n", sep = "")

  for (a in attrs) .erpm_long_deep_dump_node_attr(nw, a)

  invisible(NULL)
}