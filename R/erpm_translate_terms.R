################################################################################
# FILE: R/erpm_translate_terms.R
################################################################################
#' ERPM  translation helpers: normalize and translate RHS terms
#' @name erpm_translate_terms
#' @note erpm_translate_terms.R
#'
#' @description
#' This file implements the ERPM→{ergm} translation helpers:
#' \itemize{
#'   \item RHS syntactic splitting by `+`;
#'   \item argument normalization for ERPM-specific terms (e.g. \code{groups(...)});
#'   \item single-term translation with optional wrapping.
#' }
#'
#' @keywords ERPM ERGM translation

# ============================================================================
# ERPM term normalizers and splitters
# ============================================================================

#' Normalize `groups(...)` arguments
#' @noRd
.erpm_normalize_groups_args <- function(args_list) {
  nm <- names(args_list)

  # -- helpers ---------------------------------------------------------------
  .deparse1 <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  .eval_num_scalar <- function(x, env = parent.frame()) {
    if (is.numeric(x) && length(x) == 1L) return(x)
    vx <- try(eval(x, envir = env), silent = TRUE)
    if (inherits(vx, "try-error")) return(NULL)
    if (is.numeric(vx) && length(vx) == 1L) return(vx)
    NULL
  }
  .as_from_int <- function(x, env = parent.frame()) {
    if (is.symbol(x) && identical(x, as.name("Inf")))
      stop("groups(from,to): 'from' cannot be Inf.")
    v <- .eval_num_scalar(x, env)
    if (is.null(v) || !is.finite(v))
      stop(sprintf("groups(from): must be a finite integer >= 0. Got: %s", .deparse1(x)))
    iv <- as.integer(round(v))
    if (!isTRUE(all.equal(v, iv))) stop(sprintf("groups(from): integer required. Got: %s", format(v)))
    if (iv < 0L) stop("groups(from): must be >= 0.")
    iv
  }
  .as_to_val <- function(x, env = parent.frame()) {
    # Accept both symbol Inf and numeric Inf.
    if ((is.symbol(x) && identical(x, as.name("Inf"))) ||
        (is.numeric(x) && length(x) == 1L && is.infinite(x))) {
      return(quote(Inf))
    }
    v <- .eval_num_scalar(x, env)
    if (is.null(v) || !is.finite(v))
      stop(sprintf("groups(to): must be a finite integer or Inf. Got: %s", .deparse1(x)))
    iv <- as.integer(round(v))
    if (!isTRUE(all.equal(v, iv))) stop(sprintf("groups(to): integer or Inf required. Got: %s", format(v)))
    iv
  }
  # -------------------------------------------------------------------------

  # No args: groups ≡ [1, Inf)
  if (length(args_list) == 0L) return(list(from = 1L, to = quote(Inf)))

  # Alias support: groups(size = k)
  if (!is.null(nm) && "size" %in% nm && !("from" %in% nm) && !("to" %in% nm)) {
    args_list <- list(args_list[["size"]]); nm <- NULL
  }

  # Single positional: groups(k) ≡ [k, k+1)
  if (length(args_list) == 1L && (is.null(nm) || isTRUE(nm[1L] == ""))) {
    k <- .as_from_int(args_list[[1L]], parent.frame())
    return(list(from = k, to = k + 1L))
  }

  # Named pair: groups(from=..., to=...)
  if (!is.null(nm) && all(c("from","to") %in% nm)) {
    from <- .as_from_int(args_list[["from"]], parent.frame())
    to   <- .as_to_val  (args_list[["to"]],   parent.frame())

    # Resolve numeric to for check; keep quote(Inf) if Inf.
    to_num <- if (is.language(to)) Inf else to
    if (!(is.infinite(to_num) || (is.finite(to_num) && to_num > from))) {
      to_str <- if (is.language(to)) "Inf" else as.character(to_num)
      stop(sprintf("groups(from,to): requires 'from' < 'to'. Got: from=%d, to=%s", from, to_str))
    }
    return(list(from = from, to = to))
  }

  stop("groups(): use `groups`, `groups(k)`/`groups(size=k)`, or `groups(from=..,to=..)`. ")
}

#' Split a sum-of-terms RHS recursively
#' @noRd
.erpm_split_sum_terms <- function(expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
    return(c(.erpm_split_sum_terms(expr[[2]]),
             .erpm_split_sum_terms(expr[[3]])))
  }
  list(expr)
}

# ============================================================================
# Translation result object (call + meta)
# ============================================================================

#' Standard translator return object: {call, meta}
#'
#' @description
#' Translator functions may return either:
#'   - a raw call (legacy behavior), or
#'   - an object of class "erpm_tr" with fields {call, meta}.
#'
#' The wrapper will always use the call as the effective ergm term.
#' Metadata is attached as an attribute on the call and is ignored by ergm().
#'
#' @noRd
.erpm_tr_result <- function(call, meta = NULL) {
  stopifnot(is.call(call) || is.symbol(call))
  if (is.null(meta)) meta <- list()
  structure(list(call = call, meta = meta), class = "erpm_tr")
}

# ============================================================================
# Single-term translator ERPM → \pkg{ergm}
# ============================================================================

# Special-case translators.
#
# Rationale:
# Most terms are not constructed by the wrapper. They are implemented as ergm terms
# via InitErgmTerm.* and can be passed through (optionally renamed/wrapped).
#
# Only a small subset needs syntactic sugar at the wrapper level (e.g. groups()).

#' Translate `groups(...)` to ergm's `b2degrange(from,to)`
#' @noRd
.erpm_tr_groups <- function(fun_sym, args_list, rename_map, wrap_proj1, wrap_B, env_eval) {
  gt <- .erpm_normalize_groups_args(args_list)
  out_call <- as.call(list(as.name("b2degrange"), from = gt$from, to = gt$to))

  .erpm_tr_result(
    call = out_call,
    meta = list(
      term_name_original = "groups",
      term_name_final    = "b2degrange",
      wrappers_applied   = character(),
      notes              = sprintf(
        "normalized groups -> [from=%s, to=%s)",
        as.character(gt$from),
        if (is.language(gt$to)) "Inf" else as.character(gt$to)
      )
    )
  )
}

#' Translate `cliques(...)` by normalizing aliases only (keeps ERPM term name)
#' @noRd
.erpm_tr_cliques <- function(fun_sym, args_list, rename_map, wrap_proj1, wrap_B, env_eval) {
  # Keep a zero-arg call unchanged.
  if (length(args_list) == 0L) {
    return(.erpm_tr_result(
      call = as.call(list(as.name("cliques"))),
      meta = list(
        term_name_original = "cliques",
        term_name_final    = "cliques",
        wrappers_applied   = character(),
        notes              = "zero-arg call"
      )
    ))
  }

  al <- as.pairlist(args_list)

  # Alias: clique_size -> k (only if k is not already provided)
  if (!is.null(names(al))) {
    if (!is.null(al$clique_size) && is.null(al$k)) {
      al$k <- al$clique_size
      al$clique_size <- NULL
    }
  }

  out_call <- as.call(c(as.name("cliques"), as.list(al)))

  .erpm_tr_result(
    call = out_call,
    meta = list(
      term_name_original = "cliques",
      term_name_final    = "cliques",
      wrappers_applied   = character(),
      notes              = "normalized aliases: clique_size -> k"
    )
  )
}

# ============================================================================
# Dispatch table (term specs): translate + optional validate + optional deps
# ============================================================================

#' Term spec table: translate + optional validate + optional deps
#'
#' @description
#' Each entry can define:
#'   - translate: function(fun_sym, args_list, ..., env_eval) -> call OR erpm_tr
#'   - validate : optional function(args_list, env_eval) for early argument checks
#'   - deps     : optional dependency checks (bipartite, dyads keys, etc.)
#'
#' @noRd
.erpm_term_specs <- list(
  groups = list(
    translate = .erpm_tr_groups,
    validate  = NULL,
    deps      = list(requires_bipartite = TRUE)
  ),
  cliques = list(
    translate = .erpm_tr_cliques,
    validate  = NULL,
    deps      = list(requires_bipartite = TRUE)
  )
)

# Backward-compat alias (kept in case other internal code expects this name).
#' @noRd
.erpm_term_translators <- lapply(.erpm_term_specs, `[[`, "translate")

#' Apply spec-level dependency checks and optional early validation
#' @noRd
.erpm_apply_term_spec_checks <- function(spec, fname, args_list, env_eval) {
  deps <- spec$deps

  # Dependency: requires a bipartite `nw` (when available in env).
  if (!is.null(deps) && isTRUE(deps$requires_bipartite)) {
    nw <- try(get("nw", envir = env_eval), silent = TRUE)
    if (!inherits(nw, "try-error")) {
      bip <- try(network::get.network.attribute(nw, "bipartite"), silent = TRUE)
      if (inherits(bip, "try-error") || is.null(bip) || is.na(bip)) {
        stop(sprintf("%s(): requires a bipartite `nw` in evaluation environment.", fname))
      }
    }
  }

  # Dependency: requires named dyads attached to `nw` as `%n% "dyads"` list.
  if (!is.null(deps) && !is.null(deps$requires_dyads)) {
    nw <- try(get("nw", envir = env_eval), silent = TRUE)
    if (!inherits(nw, "try-error")) {
      dy <- try(network::get.network.attribute(nw, "dyads"), silent = TRUE)
      if (inherits(dy, "try-error") || is.null(dy)) dy <- list()
      need <- deps$requires_dyads
      miss <- setdiff(need, names(dy))
      if (length(miss)) {
        stop(sprintf("%s(): missing required dyads: %s", fname, paste(miss, collapse = ", ")))
      }
    }
  }

  # Optional early validation hook.
  if (is.function(spec$validate)) {
    spec$validate(args_list, env_eval = env_eval)
  }

  invisible(TRUE)
}

# ============================================================================
# Generic translator (rename + optional wrappers)
# ============================================================================

#' Translate a single ERPM term into an \pkg{ergm} term
#'
#' @description
#' Returns a call. If metadata is available (from a spec translator), it is
#' attached as attribute "erpm_meta" on the returned call.
#'
#' @noRd
.erpm_translate_one_term <- function(term_call,
                                     rename_map,
                                     wrap_proj1 = character(),
                                     wrap_B     = character(),
                                     env_eval   = parent.frame()) {
  # Turn bare symbols into zero-arg calls to unify processing.
  if (is.symbol(term_call)) term_call <- as.call(list(term_call))
  if (!is.call(term_call)) return(term_call)

  fun_sym   <- term_call[[1L]]
  args_list <- as.list(term_call)[-1L]

  # Determine term/function name.
  fname <- if (is.symbol(fun_sym)) as.character(fun_sym) else deparse(fun_sym)[1L]

  # --- 1) Special cases dispatch via specs -----------------------------------
  spec <- .erpm_term_specs[[fname]]
  if (!is.null(spec)) {
    .erpm_apply_term_spec_checks(spec, fname, args_list, env_eval)

    tr_res <- spec$translate(
      fun_sym, args_list,
      rename_map = rename_map,
      wrap_proj1 = wrap_proj1,
      wrap_B     = wrap_B,
      env_eval   = env_eval
    )

    # New style: {call, meta}
    if (is.list(tr_res) && inherits(tr_res, "erpm_tr")) {
      out_call <- tr_res$call
      attr(out_call, "erpm_meta") <- tr_res$meta
      return(out_call)
    }

    # Legacy style: raw call
    return(tr_res)
  }

  # --- 2) Generic path: rename + optional wrappers ---------------------------
  fname2 <- if (fname %in% names(rename_map)) rename_map[[fname]] else fname
  out <- as.call(c(as.name(fname2), args_list))

  if (fname2 %in% wrap_B)     out <- call("B",     call("~", out), form = "nonzero")
  if (fname2 %in% wrap_proj1) out <- call("Proj1", call("~", out))

  out
}