################################################################################
# FILE: R/erpm_validators.R
################################################################################
#' ERPM validators: input checks and normalized settings for erpm()
#'
#' @name erpm_validators
#' @note erpm_validators.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' This file implements the validation layer used by \code{erpm()} before any
#' bipartite network construction or \code{ergm()} call assembly happens.
#'
#' The validator is organized as:
#' \itemize{
#'   \item \strong{small helpers}: formula parsing and dyadic-name extraction;
#'   \item \strong{atomic validators}: one helper per argument family;
#'   \item \strong{normalizers}: strict-but-retrocompatible option normalization;
#'   \item \strong{orchestrator}: a single entry point returning normalized settings.
#' }
#'
#' The intent is to keep \code{erpm()} short while preserving historical behavior
#' and early-failing on malformed inputs with explicit messages.
#'
#' @keywords ERPM ERGM validators wrapper inputs
NULL
################################################################################

# ==============================================================================
# Small helpers
# ==============================================================================

#' Parse and normalize the user formula (internal helper)
#'
#' @param formula User formula.
#' @return List containing the original formula environment, LHS, RHS, and a
#'   compact string representation for logs.
#' @noRd
.erpm_parse_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("Expected a `lhs ~ ...` formula with lhs = partition OR bipartite network.")
  }

  env0 <- environment(formula)
  if (is.null(env0)) env0 <- parent.frame()

  list(
    env0             = env0,
    lhs_expr         = formula[[2L]],
    rhs_expr         = formula[[3L]],
    user_formula_str = .compact_ws(.oneline(formula))
  )
}

#' Extract dyadic names referenced in the RHS (internal helper)
#'
#' @param rhs_expr RHS expression.
#' @return Character vector of unique dyadic names referenced through
#'   \code{dyadcov*} or \code{cov_fullmatch(...)} terms.
#' @noRd
.erpm_rhs_dyad_names <- function(rhs_expr) {
  if (is.null(rhs_expr)) return(character(0))

  out <- character(0)

  walk <- function(x) {
    if (is.call(x)) {
      if (is.symbol(x[[1L]])) {
        fn <- as.character(x[[1L]])
        if (startsWith(fn, "dyadcov") || identical(fn, "cov_fullmatch")) {
          if (length(x) >= 2L) {
            a1 <- x[[2L]]
            if (is.character(a1) && length(a1) == 1L && nzchar(a1)) {
              out <<- c(out, a1)
            }
          }
        }
      }
      for (i in seq_along(x)) walk(x[[i]])
    } else if (is.pairlist(x) || is.list(x)) {
      for (i in seq_along(x)) walk(x[[i]])
    }
    invisible(NULL)
  }

  walk(rhs_expr)
  unique(out)
}

# ==============================================================================
# Atomic validators / normalizers
# ==============================================================================

#' Normalize constraints for erpm() (internal helper)
#'
#' @param constraints User constraints value.
#' @return Constraints formula.
#' @noRd
.erpm_validate_constraints <- function(constraints) {
  if (is.null(constraints)) {
    return(as.formula(~ b1part))
  }

  if (!(inherits(constraints, "formula") && length(constraints) >= 2L)) {
    stop(
      "[ERPM] `constraints` must be a formula like `~ b1part` or `~ b1part + blockdiag(timeblock)`.",
      call. = FALSE
    )
  }

  constraints
}

#' Normalize estimate for erpm() (internal helper)
#'
#' @param estimate User estimate value.
#' @return Normalized estimate or NULL.
#' @noRd
.erpm_validate_estimate <- function(estimate) {
  if (is.null(estimate)) return(NULL)

  estimate <- match.arg(estimate, c("MLE", "CD", "MPLE", "MCMLE"))
  if (identical(estimate, "MCMLE")) estimate <- "MLE"
  estimate
}

#' Normalize seed for erpm() (internal helper)
#'
#' @param seed User seed value.
#' @return Integer seed or NULL.
#' @noRd
.erpm_validate_seed <- function(seed) {
  if (is.null(seed)) return(NULL)

  if (!(is.numeric(seed) && length(seed) == 1L && is.finite(seed))) {
    stop("[ERPM] `seed` must be a single finite numeric value (integer-like) or NULL.", call. = FALSE)
  }

  seed_i <- as.integer(round(seed))
  if (!isTRUE(all.equal(seed, seed_i))) {
    stop("[ERPM] `seed` must be integer-valued (e.g., 1, 2, 42).", call. = FALSE)
  }

  seed_i
}

#' Guess a label column in a node data.frame
#'
#' @param nodes Node data.frame.
#' @param prefer Preferred label-column name.
#' @return Character scalar or NULL.
#' @noRd
.erpm_get_label_col <- function(nodes, prefer = "label") {
  stopifnot(is.data.frame(nodes))
  nms <- trimws(names(nodes))
  aliases <- c(prefer, "label", "nom", "name", "id")
  hit <- intersect(aliases, nms)
  if (length(hit)) hit[1L] else NULL
}

#' Validate a node data.frame
#'
#' @param nodes Node data.frame.
#' @return TRUE invisibly.
#' @noRd
.erpm_check_nodes_df <- function(nodes) {
  stopifnot(is.data.frame(nodes))
  lab <- .erpm_get_label_col(nodes)
  if (!is.null(lab) && anyDuplicated(nodes[[lab]])) {
    stop(sprintf("nodes$%s contains duplicates.", lab))
  }
  invisible(TRUE)
}

#' Validate dyadic n x n matrices
#'
#' @param dyads Named list of matrices.
#' @param n Actor count.
#' @param labels Actor labels.
#' @return TRUE invisibly.
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

#' Validate optional mh_moves / mh_weights passed to erpm()
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

  if (anyDuplicated(mh_moves)) {
    stop(
      "[ERPM] `mh_moves` must not contain duplicates.",
      call. = FALSE
    )
  }

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

#' Normalize dyads input for erpm() (internal helper)
#'
#' @param dyads User dyads input.
#' @param rhs_expr RHS expression.
#' @return Dyads unchanged if already a list; otherwise a named list when a
#'   single matrix can be normalized retrocompatibly.
#' @noRd
.erpm_normalize_dyads_input <- function(dyads, rhs_expr) {
  if (!is.matrix(dyads)) return(dyads)

  nm <- .erpm_rhs_dyad_names(rhs_expr)
  if (length(nm) != 1L) {
    stop(
      "[ERPM] `dyads` was provided as a matrix, but the RHS does not contain exactly one dyadic name.\n",
      "  Expected something like: dyadcov_full(\"X\") with a unique X.\n",
      "  Fix: pass `dyads = list(X = M)` or ensure the RHS contains one unique dyad name.",
      call. = FALSE
    )
  }

  setNames(list(dyads), nm)
}

# ==============================================================================
# Orchestrator
# ==============================================================================

#' Validate and normalize erpm() inputs
#'
#' @param formula See \code{erpm()}.
#' @param eval.call See \code{erpm()}.
#' @param verbose See \code{erpm()}.
#' @param estimate See \code{erpm()}.
#' @param eval.loglik See \code{erpm()}.
#' @param control See \code{erpm()}.
#' @param timeout See \code{erpm()}.
#' @param seed See \code{erpm()}.
#' @param nodes See \code{erpm()}.
#' @param dyads See \code{erpm()}.
#' @param group_labels See \code{erpm()}.
#' @param constraints See \code{erpm()}.
#' @param mh_moves See \code{erpm()}.
#' @param mh_weights See \code{erpm()}.
#'
#' @return List of normalized settings used downstream by \code{erpm()}.
#' @noRd
.erpm_validate_inputs <- function(formula,
                                  eval.call,
                                  verbose,
                                  estimate,
                                  eval.loglik,
                                  control,
                                  timeout,
                                  seed,
                                  nodes,
                                  dyads,
                                  group_labels,
                                  constraints,
                                  mh_moves,
                                  mh_weights) {
  parsed <- .erpm_parse_formula(formula)

  constraints <- .erpm_validate_constraints(constraints)
  estimate    <- .erpm_validate_estimate(estimate)
  seed        <- .erpm_validate_seed(seed)
  dyads       <- .erpm_normalize_dyads_input(dyads, parsed$rhs_expr)

  mh_spec <- .erpm_validate_mh_mix_inputs(
    mh_moves   = mh_moves,
    mh_weights = mh_weights
  )

  list(
    env0             = parsed$env0,
    lhs_expr         = parsed$lhs_expr,
    rhs_expr         = parsed$rhs_expr,
    user_formula_str = parsed$user_formula_str,
    eval.call        = eval.call,
    verbose          = verbose,
    estimate         = estimate,
    eval.loglik      = eval.loglik,
    control          = control,
    timeout          = timeout,
    seed             = seed,
    nodes            = nodes,
    dyads            = dyads,
    group_labels     = group_labels,
    constraints      = constraints,
    mh_moves         = mh_spec$mh_moves,
    mh_weights       = mh_spec$mh_weights
  )
}