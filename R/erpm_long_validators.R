################################################################################
# FILE: R/erpm_long_validators.R
################################################################################
#' ERPM long validators: input checks and derived settings for erpm_long()
#'
#' @name erpm_long_validators
#' @note erpm_long_validators.R
#'
#' @description
#' This file implements the validation layer used by \code{erpm_long()} before any
#' meta-network construction happens.
#'
#' The validator is organized as:
#' \itemize{
#'   \item \strong{usage helpers}: reusable, user-facing format expectations;
#'   \item \strong{atomic validators}: one helper per argument (easy to test);
#'   \item \strong{coherence checks}: cross-argument consistency rules;
#'   \item \strong{orchestrator}: a single entry point returning normalized settings.
#' }
#'
#' The intent is to fail early with explicit messages whenever the inputs would
#' produce ambiguous or broken meta-networks (wrong T, schema mismatch, bad dyads
#' dimensions, invalid past influence).
#'
#' @keywords ERPM ERGM longitudinal validators PLE
################################################################################

# ==============================================================================
# Usage blocks (error messages)
# ==============================================================================

#' Usage string for the nodes argument (internal helper)
#' @return Character scalar, multi-line usage hint.
#' @noRd
.erpm_long_usage_nodes <- function() {
  paste(
    "Expected: nodes = NULL or list(data.frame) of length T, one data.frame per time.",
    "Each data.frame must have:",
    "  - a 'label' column (character or atomic),",
    "  - at least one covariate column besides 'label',",
    "  - nrow(nodes[[t]]) == length(partitions[[t]]).",
    "",
    "Example:",
    "nodes <- list(",
    "  data.frame(label=c('A','B','C','D'), gender=c(1,1,2,1), age=c(20,22,25,30)),",
    "  data.frame(label=c('FT','AZ','JI','DO'), gender=c(2,1,2,2), age=c(10,42,25,30)),",
    "  data.frame(label=c('H','Z','S','A'),  gender=c(1,1,1,1), age=c(27,26,25,28))",
    ")",
    sep = "\n"
  )
}

#' Usage string for the dyads argument (internal helper)
#' @return Character scalar, multi-line usage hint.
#' @noRd
.erpm_long_usage_dyads <- function() {
  paste(
    "Expected: dyads = NULL or list(list(matrix)) of length T.",
    "Format: dyads[[t]] is a NAMED list of dyadic matrices (e.g., fm=..., Z1=...).",
    "Constraints:",
    "  - length(dyads) == T",
    "  - names(dyads[[t]]) must be non-empty",
    "  - each dyads[[t]][[name]] must be a numeric matrix",
    "  - each matrix must be square with dim nA_t x nA_t where nA_t = length(partitions[[t]]).",
    "",
    "Example:",
    "dyads <- list(",
    "  list(fm=matrix(..., nrow=4, byrow=TRUE), Z1=matrix(..., nrow=4, byrow=TRUE)),",
    "  list(fm=matrix(..., nrow=4, byrow=TRUE), Z1=matrix(..., nrow=4, byrow=TRUE)),",
    "  list(fm=matrix(..., nrow=4, byrow=TRUE), Z1=matrix(..., nrow=4, byrow=TRUE))",
    ")",
    sep = "\n"
  )
}

# ==============================================================================
# Small helpers
# ==============================================================================

.erpm_long_stop <- function(...) stop(sprintf(...), call. = FALSE)

.erpm_long_is_scalar_logical <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

.erpm_long_is_scalar_character_or_null <- function(x) {
  is.null(x) || (is.character(x) && length(x) == 1L && !is.na(x))
}

# ==============================================================================
# Argument validators (one per argument)
# ==============================================================================

# ---- mode --------------------------------------------------------------------

.erpm_long_validate_mode <- function(mode) {
  # Accepts: PLE/PLS/empile/sequential (case-insensitive), with synonyms:
  # PLE == empile ; PLS == sequential.
  if (is.null(mode)) .erpm_long_stop("[ERPM_LONG] mode must be provided.")

  if (length(mode) != 1L || is.na(mode)) {
    .erpm_long_stop("[ERPM_LONG] mode must be a single string among: PLE, empile, PLS, sequential.")
  }

  m <- tolower(as.character(mode))

  if (m %in% c("ple", "empile")) {
    return(list(mode_user = mode, mode_norm = "PLE"))
  }

  if (m %in% c("pls", "sequential")) {
    .erpm_long_stop("[ERPM_LONG] PLS/sequential mode is not implemented yet. Use PLE/empile.")
  }

  .erpm_long_stop("[ERPM_LONG] Invalid mode='%s'. Allowed: PLE, empile, PLS, sequential.", as.character(mode))
}

# ---- formula + partitions + rhs ----------------------------------------------

.erpm_long_validate_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    .erpm_long_stop("[ERPM_LONG] formula must be a formula.")
  }
  if (length(formula) < 3L) {
    .erpm_long_stop("[ERPM_LONG] formula must be of the form: partitions_list ~ terms")
  }
  invisible(TRUE)
}

.erpm_long_eval_partitions_from_lhs <- function(formula) {
  lhs <- formula[[2L]]
  if (is.null(lhs)) .erpm_long_stop("[ERPM_LONG] Missing LHS in formula.")

  partitions <- eval(lhs, envir = parent.frame())

  if (!is.list(partitions) || !length(partitions)) {
    .erpm_long_stop("[ERPM_LONG] LHS must evaluate to a non-empty list of partitions.")
  }

  if (length(partitions) == 1L) {
    .erpm_long_stop("[ERPM_LONG] Only one partition provided (T=1). Use erpm() instead of erpm_long().")
  }

  # Per your rule: list of partitions; atomic; no NA.
  for (t in seq_along(partitions)) {
    p <- partitions[[t]]
    if (is.null(p) || is.list(p) || !is.atomic(p)) {
      .erpm_long_stop("[ERPM_LONG] partitions[[%d]] must be an atomic vector (not a list).", t)
    }
    if (!length(p)) {
      .erpm_long_stop("[ERPM_LONG] partitions[[%d]] is empty.", t)
    }
    if (anyNA(p)) {
      .erpm_long_stop("[ERPM_LONG] partitions[[%d]] contains NA values.", t)
    }
  }

  partitions
}

.erpm_long_validate_rhs <- function(formula) {
  rhs <- formula[[3L]]
  if (is.null(rhs)) .erpm_long_stop("[ERPM_LONG] Missing RHS in formula.")

  tt <- terms(as.formula(call("~", rhs)))
  term_labels <- attr(tt, "term.labels")

  if (is.null(term_labels) || !length(term_labels)) {
    .erpm_long_stop("[ERPM_LONG] RHS must contain at least one term (use '+').")
  }

  rhs
}

# ---- inertial detection + past_influence -------------------------------------

.erpm_long_detect_inertial <- function(rhs) {
  # Conservative: recognizes inertia_groups(...) only (your current behavior).
  tt <- terms(as.formula(call("~", rhs)))
  term_labels <- attr(tt, "term.labels")

  idx <- grep("^inertia_groups\\b", term_labels)
  if (!length(idx)) return(list(inertial_present = FALSE, d = 0L))

  dmax <- 1L

  for (lab in term_labels[idx]) {
    expr <- try(parse(text = lab)[[1L]], silent = TRUE)
    if (inherits(expr, "try-error") || !is.call(expr)) next
    if (!identical(as.character(expr[[1L]]), "inertia_groups")) next

    args <- as.list(expr)[-1L]

    # Default past_influence=1 if missing
    d_i <- 1L

    # Allowed aliases (based on your current wrapper behavior)
    if (length(args)) {
      if ("past_influence" %in% names(args)) {
        d_i <- suppressWarnings(as.integer(round(eval(args[["past_influence"]], parent.frame()))))
      } else if ("pi" %in% names(args)) {
        d_i <- suppressWarnings(as.integer(round(eval(args[["pi"]], parent.frame()))))
      } else if ("d" %in% names(args)) {
        d_i <- suppressWarnings(as.integer(round(eval(args[["d"]], parent.frame()))))
      } else {
        d_i <- suppressWarnings(as.integer(round(eval(args[[1L]], parent.frame()))))
      }
      if (is.na(d_i)) d_i <- 1L
    }

    if (d_i < 1L) {
      .erpm_long_stop("[ERPM_LONG] inertia_groups: past_influence must be >= 1.")
    }

    dmax <- max(dmax, d_i)
  }

  list(inertial_present = TRUE, d = as.integer(dmax))
}

.erpm_long_validate_past_influence_vs_T <- function(T, inertial_present, d) {
  if (!isTRUE(inertial_present)) return(invisible(TRUE))
  if (is.na(d) || d < 1L) {
    .erpm_long_stop("[ERPM_LONG] past_influence must be >= 1 when inertial terms are present.")
  }
  if (T < (d + 1L)) {
    .erpm_long_stop("[ERPM_LONG] T=%d but past_influence=%d: need at least T >= past_influence + 1.", T, d)
  }
  invisible(TRUE)
}

# ---- simple flags ------------------------------------------------------------

.erpm_long_validate_eval_call <- function(eval.call) {
  if (!.erpm_long_is_scalar_logical(eval.call)) {
    .erpm_long_stop("[ERPM_LONG] eval.call must be TRUE or FALSE.")
  }
  invisible(TRUE)
}

.erpm_long_validate_verbose <- function(verbose) {
  if (!.erpm_long_is_scalar_logical(verbose)) {
    .erpm_long_stop("[ERPM_LONG] verbose must be TRUE or FALSE.")
  }
  invisible(TRUE)
}

.erpm_long_validate_debug <- function(debug) {
  if (!.erpm_long_is_scalar_logical(debug)) {
    .erpm_long_stop("[ERPM_LONG] debug must be TRUE or FALSE.")
  }
  invisible(TRUE)
}

# ---- seed --------------------------------------------------------------------

.erpm_long_validate_seed <- function(seed) {
  if (is.null(seed)) return(invisible(TRUE))

  if (length(seed) != 1L || is.na(seed)) {
    .erpm_long_stop("[ERPM_LONG] seed must be NULL or a scalar value compatible with set.seed().")
  }

  old <- NULL
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv)
  }, add = TRUE)

  ok <- !inherits(try(set.seed(seed), silent = TRUE), "try-error")
  if (!ok) {
    .erpm_long_stop("[ERPM_LONG] seed is not compatible with set.seed(): got class=%s.", paste(class(seed), collapse = "/"))
  }

  invisible(TRUE)
}

# ---- group_labels ------------------------------------------------------------

.erpm_long_validate_group_labels <- function(group_labels) {
  if (!.erpm_long_is_scalar_character_or_null(group_labels)) {
    .erpm_long_stop("[ERPM_LONG] group_labels must be NULL or a single character string.")
  }
  invisible(TRUE)
}

# ---- nodes -------------------------------------------------------------------

.erpm_long_validate_nodes <- function(nodes, partitions) {
  if (is.null(nodes)) return(invisible(TRUE))

  T <- length(partitions)

  if (!is.list(nodes) || length(nodes) != T) {
    .erpm_long_stop("[ERPM_LONG] Invalid nodes.\n\n%s", .erpm_long_usage_nodes())
  }

  cols_ref <- NULL

  for (t in seq_len(T)) {
    df <- nodes[[t]]
    if (!is.data.frame(df)) {
      .erpm_long_stop("[ERPM_LONG] Invalid nodes: nodes[[%d]] must be a data.frame.\n\n%s", t, .erpm_long_usage_nodes())
    }

    nA_t <- length(partitions[[t]])
    if (nrow(df) != nA_t) {
      .erpm_long_stop(
        "[ERPM_LONG] Invalid nodes: nodes[[%d]] has %d rows but partitions[[%d]] has %d actors.\n\n%s",
        t, nrow(df), t, nA_t, .erpm_long_usage_nodes()
      )
    }

    if (!("label" %in% colnames(df))) {
      .erpm_long_stop("[ERPM_LONG] Invalid nodes: nodes[[%d]] must contain a 'label' column.\n\n%s", t, .erpm_long_usage_nodes())
    }

    cov_names <- setdiff(colnames(df), "label")
    if (!length(cov_names)) {
      .erpm_long_stop("[ERPM_LONG] Invalid nodes: nodes[[%d]] must have at least one covariate besides 'label'.\n\n%s", t, .erpm_long_usage_nodes())
    }

    # Enforce schema consistency (same columns, same order) across time
    if (is.null(cols_ref)) cols_ref <- colnames(df)
    if (!identical(colnames(df), cols_ref)) {
      .erpm_long_stop("[ERPM_LONG] Invalid nodes: column schema differs at t=%d (must match t=1).", t)
    }

    lab <- df[["label"]]
    if (anyNA(lab) || any(!nzchar(as.character(lab)))) {
      .erpm_long_stop("[ERPM_LONG] Invalid nodes: nodes[[%d]]$label contains NA/empty values.", t)
    }
  }

  invisible(TRUE)
}

# ---- dyads -------------------------------------------------------------------

.erpm_long_validate_dyads <- function(dyads, partitions) {
  if (is.null(dyads)) return(invisible(TRUE))

  T <- length(partitions)

  if (!is.list(dyads) || length(dyads) != T) {
    .erpm_long_stop("[ERPM_LONG] Invalid dyads.\n\n%s", .erpm_long_usage_dyads())
  }

  for (t in seq_len(T)) {
    dt <- dyads[[t]]

    if (!is.list(dt) || !length(dt) || is.null(names(dt)) || any(!nzchar(names(dt)))) {
      .erpm_long_stop("[ERPM_LONG] Invalid dyads: dyads[[%d]] must be a NAMED list of matrices.\n\n%s", t, .erpm_long_usage_dyads())
    }

    nA_t <- length(partitions[[t]])

    for (nm in names(dt)) {
      M <- dt[[nm]]
      if (!is.matrix(M)) {
        .erpm_long_stop("[ERPM_LONG] Invalid dyads: dyads[[%d]][['%s']] must be a matrix.\n\n%s", t, nm, .erpm_long_usage_dyads())
      }
      if (!is.numeric(M)) {
        .erpm_long_stop("[ERPM_LONG] Invalid dyads: dyads[[%d]][['%s']] must be numeric.", t, nm)
      }
      if (nrow(M) != nA_t || ncol(M) != nA_t) {
        .erpm_long_stop(
          "[ERPM_LONG] Invalid dyads: dyads[[%d]][['%s']] has dim %dx%d; expected %dx%d (nA_t=%d).",
          t, nm, nrow(M), ncol(M), nA_t, nA_t, nA_t
        )
      }
      if (any(!is.finite(M))) {
        .erpm_long_stop("[ERPM_LONG] Invalid dyads: dyads[[%d]][['%s']] contains non-finite values.", t, nm)
      }
    }
  }

  invisible(TRUE)
}

# ==============================================================================
# Inter-input coherence checks
# ==============================================================================

.erpm_long_validate_coherence <- function(partitions, nodes, dyads) {
  # At this point:
  # - partitions: list length T>=2
  # - nodes: NULL or list length T with consistent schema and row counts
  # - dyads: NULL or list length T with per-time square matrices matching partitions[[t]]
  # This function keeps any cross-argument rules in one place.

  T <- length(partitions)

  if (!is.null(nodes) && length(nodes) != T) {
    .erpm_long_stop("[ERPM_LONG] Coherence error: nodes length != T.")
  }

  if (!is.null(dyads) && length(dyads) != T) {
    .erpm_long_stop("[ERPM_LONG] Coherence error: dyads length != T.")
  }

  invisible(TRUE)
}

# ==============================================================================
# Orchestrator: main validator called by erpm_long()
# ==============================================================================

.erpm_long_validate_inputs <- function(formula,
                                      mode,
                                      eval.call,
                                      verbose,
                                      debug,
                                      estimate,
                                      eval.loglik,
                                      control,
                                      timeout,
                                      seed,
                                      nodes,
                                      dyads,
                                      group_labels) {

  # Formula shape
  .erpm_long_validate_formula(formula)

  # Mode normalization + (PLS stop)
  m <- .erpm_long_validate_mode(mode)

  # LHS/RHS
  partitions <- .erpm_long_eval_partitions_from_lhs(formula)
  rhs        <- .erpm_long_validate_rhs(formula)

  # Flags
  .erpm_long_validate_eval_call(eval.call)
  .erpm_long_validate_verbose(verbose)
  .erpm_long_validate_debug(debug)

  # seed + group_labels
  .erpm_long_validate_seed(seed)
  .erpm_long_validate_group_labels(group_labels)

  # nodes/dyads formats
  .erpm_long_validate_nodes(nodes, partitions)
  .erpm_long_validate_dyads(dyads, partitions)

  # inertial detection + past_influence constraints
  inert <- .erpm_long_detect_inertial(rhs)
  inertial_present <- isTRUE(inert$inertial_present)
  d <- as.integer(inert$d)

  .erpm_long_validate_past_influence_vs_T(T = length(partitions), inertial_present = inertial_present, d = d)

  # coherence across inputs
  .erpm_long_validate_coherence(partitions, nodes, dyads)

  # estimate/eval.loglik/control/timeout: pass-through by design
  list(
    lhs              = formula[[2L]],
    rhs              = rhs,
    partitions       = partitions,
    T                = length(partitions),
    inertial_present = inertial_present,
    past_influence   = if (inertial_present) d else 0L,
    mode_user        = mode,
    mode_norm        = m$mode_norm
  )
}