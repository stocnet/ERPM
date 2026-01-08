################################################################################
# FILE: R/erpm_long_inertia_utils.R
################################################################################
#' ERPM longitudinal inertia : helpers 
#' @name erpm_long_inertia_utils
#' @note erpm_long_inertia_utils.R
#' 
#' @description
#' This file groups small internal helpers used by longitudinal ERPM 
#' code to:
#'   - parse a RHS expression into individual terms (robustly, even when nested);
#'   - extract stable term names for user-facing error messages;
#'   - safely evaluate user-provided arguments (integers, named options);
#'   - derive and validate history summaries from an ERPM bipartite network:
#'       partition extraction, group signatures, co-membership matrices.
#'
#' Conventions:
#'   - All functions are internal and start with ".erpm_long_".
#'   - Errors are explicit and user-facing (prefixed with \code{ERPM_LONG}).
#'   - The "bipartite" network is assumed to follow build_bipartite_from_inputs()
#'     conventions: first mode are actors, second mode are groups.
#'
#' @keywords ERPM ERGM longitudinal inertia internal helpers
NULL

# =============================================================================
# Helpers
# =============================================================================

#' Split a sum-of-terms RHS recursively (local copy for robustness)
#'
#' Internal utility that turns an expression like `a + b + f(x)` into a flat list
#' of RHS items, while preserving non-`+` expressions as singletons.
#' This is duplicated locally to avoid tight coupling with other wrapper files.
#'
#' @param expr An R language object (typically a parsed RHS expression).
#' @return A list of RHS items (language objects).
#' @noRd
.erpm_long_split_sum_terms <- function(expr) {
  # If the expression is a call to `+`, split left and right recursively.
  if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
    return(c(.erpm_long_split_sum_terms(expr[[2]]),
             .erpm_long_split_sum_terms(expr[[3]])))
  }
  # Base case: a single term, returned as a length-1 list.
  list(expr)
}

#' Extract a term/function name from a RHS item
#'
#' Used for consistent, readable error messages. Accepts both symbols (`edges`)
#' and calls (`inertia(past_influence=2)`), returning the function/symbol name.
#'
#' @param tt A RHS item as a language object.
#' @return A character scalar (term name) or NA_character_ if unknown.
#' @noRd
.erpm_long_term_name <- function(tt) {
  # Bare symbol: "edges", "nodematch", etc.
  if (is.symbol(tt)) return(as.character(tt))
  # Call: first element is the function name symbol.
  if (is.call(tt) && is.symbol(tt[[1L]])) return(as.character(tt[[1L]]))
  # Fallback when structure is unexpected.
  NA_character_
}

#' Evaluate an integer argument (positional) from an inertial call
#'
#' Validates the presence of a positional argument, evaluates it in env0, and
#' enforces "finite numeric scalar" + "integer-valued" + lower bound.
#' This prevents silent coercions that would make longitudinal logic ambiguous.
#'
#' @param call A term call (language object), e.g. inertia(2).
#' @param env0 Evaluation environment for user expressions.
#' @param pos 1-based position among call arguments (excluding function name).
#' @param name User-facing name for error messages.
#' @param min Minimum allowed integer value (inclusive).
#' @return An integer scalar.
#' @noRd
.erpm_long_eval_int_arg <- function(call, env0, pos, name, min = 0L) {
  stopifnot(is.call(call)) # Defensive: this helper expects a call.

  # Extract only arguments (drop function name at position 1).
  args <- as.list(call)[-1L]

  # Enforce required positional argument.
  if (length(args) < pos) {
    stop(sprintf("[ERPM_LONG] %s(): missing argument '%s' (pos=%d).",
                 .erpm_long_term_name(call), name, pos),
         call. = FALSE)
  }

  # Evaluate the user expression in the provided environment (may reference objects).
  x_expr <- args[[pos]]
  x_val  <- try(eval(x_expr, envir = env0), silent = TRUE)
  if (inherits(x_val, "try-error")) {
    stop(sprintf("[ERPM_LONG] %s(): cannot evaluate '%s'.",
                 .erpm_long_term_name(call), name),
         call. = FALSE)
  }

  # Require a single finite numeric value (reject vectors, NA, Inf).
  if (!(is.numeric(x_val) && length(x_val) == 1L && is.finite(x_val))) {
    stop(sprintf("[ERPM_LONG] %s(): '%s' must be a finite numeric scalar.",
                 .erpm_long_term_name(call), name),
         call. = FALSE)
  }

  # Convert to integer only if it is truly integer-valued (e.g. 2.0 OK, 2.5 not OK).
  xi <- as.integer(round(x_val))
  if (!isTRUE(all.equal(x_val, xi))) {
    stop(sprintf("[ERPM_LONG] %s(): '%s' must be integer-valued.",
                 .erpm_long_term_name(call), name),
         call. = FALSE)
  }

  # Enforce lower bound to keep downstream history logic well-defined.
  if (xi < as.integer(min)) {
    stop(sprintf("[ERPM_LONG] %s(): '%s' must be >= %d.",
                 .erpm_long_term_name(call), name, as.integer(min)),
         call. = FALSE)
  }

  xi
}

#' Evaluate named argument `past_influence` (default=1)
#'
#' past_influence = d means:
#' - the effect becomes active only when t > d
#' - we compute and attach d summaries, one per lag 1..d
#'
#' @param call A term call (language object) that may include past_influence=...
#' @param env0 Evaluation environment for user expressions.
#' @param default Default lag depth if not specified by the user.
#' @return An integer scalar d >= 1.
#' @noRd
.erpm_long_get_past_influence <- function(call, env0, default = 1L) {
  stopifnot(is.call(call)) # This helper expects a call with named args.

  # Extract arguments and their names (may be NULL for purely positional calls).
  args <- as.list(call)[-1L]
  nm <- names(args)

  # If user did not specify past_influence, use default.
  if (is.null(nm) || !("past_influence" %in% nm)) {
    d <- as.integer(default)
  } else {
    # Evaluate in the user's environment to allow symbols/expressions.
    d_expr <- args[["past_influence"]]
    d_val  <- try(eval(d_expr, envir = env0), silent = TRUE)
    if (inherits(d_val, "try-error")) {
      stop(sprintf("[ERPM_LONG] %s(): cannot evaluate 'past_influence'.",
                   .erpm_long_term_name(call)),
           call. = FALSE)
    }

    # Require a single finite numeric scalar.
    if (!(is.numeric(d_val) && length(d_val) == 1L && is.finite(d_val))) {
      stop(sprintf("[ERPM_LONG] %s(): 'past_influence' must be a finite numeric scalar.",
                   .erpm_long_term_name(call)),
           call. = FALSE)
    }

    # Enforce integer-valued semantics.
    d <- as.integer(round(d_val))
    if (!isTRUE(all.equal(d_val, d))) {
      stop(sprintf("[ERPM_LONG] %s(): 'past_influence' must be integer-valued.",
                   .erpm_long_term_name(call)),
           call. = FALSE)
    }
  }

  # Longitudinal history depth must be at least 1 (lag-0 is not "past").
  if (d < 1L) {
    stop(sprintf("[ERPM_LONG] %s(): 'past_influence' must be >= 1.",
                 .erpm_long_term_name(call)),
         call. = FALSE)
  }

  as.integer(d)
}

#' Debug print helper (internal)
#'
#' Centralizes conditional debug prints so callers do not sprinkle `if (debug)`
#' everywhere. Always returns invisibly to avoid interfering with pipelines.
#'
#' @param debug Logical flag controlling output.
#' @param ... Items passed to cat().
#' @return Invisible NULL.
#' @noRd
.erpm_long_dbg <- function(debug, ...) {
  if (isTRUE(debug)) cat(..., "\n", sep = "")
  invisible(NULL)
}

#' Safe check: is this a bipartite membership network of the ERPM construction
#'
#' Minimal structural check used before assuming ERPM-specific vertex layout.
#' We require a {network} object and a numeric "bipartite" attribute.
#'
#' @param nw A candidate network object.
#' @return TRUE/FALSE.
#' @noRd
.erpm_long_is_erpm_bipartite <- function(nw) {
  inherits(nw, "network") &&
    is.numeric(network::get.network.attribute(nw, "bipartite"))
}

#' Extract an actor->group partition vector from a bipartite membership network. (could be used in very specific cases)
#'
#' Convention assumed (as in build_bipartite_from_inputs()):
#' - vertices 1..n1 are actors
#' - vertices (n1+1)..(n1+n2) are groups
#' - each actor has exactly one incident membership edge to its group
#'
#' Returned partition uses group ids in 1..n2.
#'
#' @param nw_prev A bipartite membership {network} built by ERPM.
#' @return Integer vector p of length n1, where \code{p[i]} is the group index (1..n2).
#' @noRd
.erpm_long_extract_partition_from_network <- function(nw_prev) {
  # Reject early if the object cannot be an ERPM bipartite network.
  if (!.erpm_long_is_erpm_bipartite(nw_prev)) {
    stop("[ERPM_LONG] cannot extract partition: not a bipartite network.", call. = FALSE)
  }

  # In {network}, attribute "bipartite" stores the size of mode-1 (actors).
  n1 <- as.integer(network::get.network.attribute(nw_prev, "bipartite"))
  if (!is.finite(n1) || n1 < 1L) {
    stop("[ERPM_LONG] cannot extract partition: invalid 'bipartite' attribute.", call. = FALSE)
  }

  # Convert the network to an edgelist matrix with two columns (tail, head).
  el <- network::as.edgelist(nw_prev)
  if (!is.matrix(el) || ncol(el) != 2L) {
    stop("[ERPM_LONG] cannot extract partition: invalid edgelist.", call. = FALSE)
  }

  # Normalize undirected orientation so first column is actor when possible.
  # This makes downstream filters simpler and stable across representations.
  a <- el[, 1L]
  b <- el[, 2L]
  swap <- (a > n1) & (b <= n1) # group->actor edges that should be flipped to actor->group
  if (any(swap)) {
    tmp <- a[swap]; a[swap] <- b[swap]; b[swap] <- tmp
  }

  # Keep only edges actor->group (actor in 1..n1, group > n1).
  # Any extra edges (if present) are ignored rather than silently misused.
  keep <- (a >= 1L & a <= n1) & (b > n1)
  a <- a[keep]
  b <- b[keep]

  if (!length(a)) {
    stop("[ERPM_LONG] cannot extract partition: no actor->group edges found.", call. = FALSE)
  }

  # Each actor should appear exactly once.
  # We enforce this strictly because inertial summaries assume single membership.
  tab <- table(a)
  if (any(tab != 1L)) {
    stop("[ERPM_LONG] cannot extract partition: actors do not have exactly one membership edge.", call. = FALSE)
  }

  # Map group vertex id -> group index 1..n2 by sorting unique group vertices.
  # Sorting ensures stable group indices across runs given the same network.
  gverts <- sort(unique(b))
  gid_map <- setNames(seq_along(gverts), as.character(gverts))

  # Build partition vector p where position is actor id and value is group index.
  p <- integer(n1)
  p[a] <- unname(gid_map[as.character(b)])

  # If an actor index is missing, its p value stays 0 and we fail explicitly.
  if (any(p < 1L)) {
    stop("[ERPM_LONG] cannot extract partition: missing memberships for some actors.", call. = FALSE)
  }

  p
}

#' Partition -> list of groups as integer actor index vectors
#'
#' Converts a partition vector p into a list where each element is the set of
#' actors belonging to the same group. Group labels are taken from p values.
#'
#' @param p Integer-like partition vector of length n (actors).
#' @return A list of integer vectors (actor indices), one per group.
#' @noRd
.erpm_long_groups_from_partition <- function(p) {
  p <- as.integer(p)         # Normalize storage type early.
  split(seq_along(p), p)     # Group actors by their partition label.
}

#' Groups -> stable signatures (character), e.g. "1,3,5"
#'
#' Used to compare groups across time by their member set. Sorting ensures that
#' equivalent groups yield identical signatures regardless of original order.
#'
#' @param groups List of integer vectors (actor indices).
#' @return Character vector of signatures, one per group.
#' @noRd
.erpm_long_group_signatures <- function(groups) {
  vapply(groups, function(v) paste(sort(as.integer(v)), collapse = ","), character(1))
}

#' Build a co-membership matrix (pairs) from a partition.
#'
#' Returns an n x n integer matrix with 1 if same group, 0 otherwise, diag=0.
#'
#' @param p Integer-like partition vector of length n.
#' @return An n x n integer matrix (0/1), with a zero diagonal.
#' @noRd
.erpm_long_comembership_matrix <- function(p) {
  p <- as.integer(p)             # Ensure comparisons are stable.
  n <- length(p)                # n actors.
  M <- outer(p, p, "==")         # Logical matrix: TRUE when in same group.
  diag(M) <- FALSE               # Exclude self-pairs by convention.
  storage.mode(M) <- "integer"   # Store as 0/1 for compactness and C interop.
  M
}

#' Evaluate named argument `size` if present; else return NULL.
#' The semantics are left to the changestat. Here we only forward user intent.
#'
#' Accepts numeric/integer input and rounds to integer. Returns NULL when the
#' argument is missing or explicitly set to NULL by the user.
#'
#' @param call A term call that may include size=...
#' @param env0 Evaluation environment for user expressions.
#' @return Integer vector/scalar (depending on user input) or NULL.
#' @noRd
.erpm_long_get_size_filter <- function(call, env0) {
  stopifnot(is.call(call)) # This helper expects a call.

  # Extract arguments and check if "size" is present among named arguments.
  args <- as.list(call)[-1L]
  nm <- names(args)
  if (is.null(nm) || !("size" %in% nm)) return(NULL)

  # Evaluate user expression in the provided environment.
  s_expr <- args[["size"]]
  s_val  <- try(eval(s_expr, envir = env0), silent = TRUE)
  if (inherits(s_val, "try-error")) {
    stop(sprintf("[ERPM_LONG] %s(): cannot evaluate 'size'.",
                 .erpm_long_term_name(call)),
         call. = FALSE)
  }

  # NULL means "no filter" and is forwarded as-is.
  if (is.null(s_val)) return(NULL)

  # Numeric is accepted and coerced to integer by rounding (caller defines semantics).
  if (is.numeric(s_val)) {
    s_val <- as.integer(round(s_val))
    return(s_val)
  }

  # Reject other types to avoid ambiguous downstream behavior.
  stop(sprintf("[ERPM_LONG] %s(): 'size' must be numeric/integer (or NULL).",
               .erpm_long_term_name(call)),
       call. = FALSE)
}