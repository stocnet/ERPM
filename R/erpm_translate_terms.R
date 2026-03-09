################################################################################
# FILE: R/erpm_translate_terms.R
################################################################################
#' ERPM term translation helpers: normalize and translate RHS terms
#'
#' @name erpm_translate_terms
#' @note erpm_translate_terms.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' This file implements the RHS translation layer used by \code{erpm()}.
#'
#' In the current ERPM API, translation is deliberately minimal:
#' \itemize{
#'   \item \code{groups(...)} is syntactic sugar and is translated to
#'         \code{b2degrange(from, to)};
#'   \item all other ERPM effects are expected to be handled by their own
#'         \code{InitErgmTerm.*} definitions and are therefore passed through
#'         unchanged, apart from optional generic renaming/wrapping hooks kept
#'         for backward compatibility.
#' }
#'
#' The module is organized as:
#' \itemize{
#'   \item term splitters and argument normalizers;
#'   \item a small translator-result class carrying optional metadata;
#'   \item single-term translation;
#'   \item RHS pipeline reconstruction back into a formula.
#' }
#'
#' @keywords ERPM ERGM term translation groups
NULL
################################################################################

# ==============================================================================
# Term normalizers and splitters
# ==============================================================================

#' Normalize `groups(...)` arguments
#'
#' IMPORTANT (env handling):
#' This normalizer may need to evaluate user-provided expressions, e.g.
#' \preformatted{
#'   k <- 3
#'   erpm(partition ~ groups(k))
#' }
#'
#' Therefore, evaluation must occur in the wrapper evaluation environment
#' (\code{env_eval}), whose parent is the user's formula environment.
#' Using \code{parent.frame()} here is fragile because it depends on the
#' internal call stack.
#'
#' @param args_list Raw argument list from the call.
#' @param env_eval Evaluation environment.
#' @return List with entries \code{from} and \code{to}.
#' @noRd
.erpm_normalize_groups_args <- function(args_list, env_eval) {
  nm <- names(args_list)

  .deparse1 <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

  .eval_num_scalar <- function(x, env = env_eval) {
    if (is.numeric(x) && length(x) == 1L) return(x)
    vx <- try(eval(x, envir = env), silent = TRUE)
    if (inherits(vx, "try-error")) return(NULL)
    if (is.numeric(vx) && length(vx) == 1L) return(vx)
    NULL
  }

  .as_from_int <- function(x, env = env_eval) {
    if (is.symbol(x) && identical(x, as.name("Inf"))) {
      stop("groups(from,to): 'from' cannot be Inf.")
    }
    v <- .eval_num_scalar(x, env)
    if (is.null(v) || !is.finite(v)) {
      stop(sprintf("groups(from): must be a finite integer >= 0. Got: %s", .deparse1(x)))
    }
    iv <- as.integer(round(v))
    if (!isTRUE(all.equal(v, iv))) {
      stop(sprintf("groups(from): integer required. Got: %s", format(v)))
    }
    if (iv < 0L) stop("groups(from): must be >= 0.")
    iv
  }

  .as_to_val <- function(x, env = env_eval) {
    if ((is.symbol(x) && identical(x, as.name("Inf"))) ||
        (is.numeric(x) && length(x) == 1L && is.infinite(x))) {
      return(quote(Inf))
    }

    v <- .eval_num_scalar(x, env)
    if (is.null(v) || !is.finite(v)) {
      stop(sprintf("groups(to): must be a finite integer or Inf. Got: %s", .deparse1(x)))
    }

    iv <- as.integer(round(v))
    if (!isTRUE(all.equal(v, iv))) {
      stop(sprintf("groups(to): integer or Inf required. Got: %s", format(v)))
    }

    iv
  }

  if (length(args_list) == 0L) {
    return(list(from = 1L, to = quote(Inf)))
  }

  if (!is.null(nm) && "size" %in% nm && !("from" %in% nm) && !("to" %in% nm)) {
    args_list <- list(args_list[["size"]])
    nm <- NULL
  }

  if (length(args_list) == 1L && (is.null(nm) || isTRUE(nm[1L] == ""))) {
    k <- .as_from_int(args_list[[1L]], env_eval)
    return(list(from = k, to = k + 1L))
  }

  if (!is.null(nm) && all(c("from", "to") %in% nm)) {
    from <- .as_from_int(args_list[["from"]], env_eval)
    to   <- .as_to_val(args_list[["to"]], env_eval)

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
#'
#' @param expr RHS expression.
#' @return List of individual term expressions.
#' @noRd
.erpm_split_sum_terms <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], as.name("+"))) {
    return(c(
      .erpm_split_sum_terms(expr[[2L]]),
      .erpm_split_sum_terms(expr[[3L]])
    ))
  }
  list(expr)
}

# ==============================================================================
# Translation result object
# ==============================================================================

#' Standard translator return object: \code{{call, meta}}
#'
#' @param call Call or symbol to be used as effective ergm term.
#' @param meta Optional list of metadata.
#' @return Object of class \code{"erpm_tr"}.
#' @noRd
.erpm_tr_result <- function(call, meta = NULL) {
  stopifnot(is.call(call) || is.symbol(call))
  if (is.null(meta)) meta <- list()
  structure(list(call = call, meta = meta), class = "erpm_tr")
}

# ==============================================================================
# Single-term translation
# ==============================================================================

#' Translate `groups(...)` to ergm's `b2degrange(from, to)`
#'
#' @param fun_sym Function symbol.
#' @param args_list Raw argument list.
#' @param rename_map Optional rename map.
#' @param wrap_proj1 Optional Proj1 wrappers.
#' @param wrap_B Optional B wrappers.
#' @param env_eval Evaluation environment.
#' @return \code{erpm_tr} object.
#' @noRd
.erpm_tr_groups <- function(fun_sym, args_list, rename_map, wrap_proj1, wrap_B, env_eval) {
  gt <- .erpm_normalize_groups_args(args_list, env_eval = env_eval)
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

#' Wrapper-level term specification table
#'
#' @description
#' Only terms that need wrapper-level syntactic sugar should appear here.
#' In the current design this is limited to \code{groups(...)}.
#'
#' @noRd
.erpm_term_specs <- list(
  groups = list(
    translate = .erpm_tr_groups,
    validate  = NULL,
    deps      = list(requires_bipartite = TRUE)
  )
)

#' Backward-compat alias for translator access
#' @noRd
.erpm_term_translators <- lapply(.erpm_term_specs, `[[`, "translate")

#' Apply spec-level dependency checks and optional validation
#'
#' @param spec Term spec entry.
#' @param fname Term name.
#' @param args_list Raw argument list.
#' @param env_eval Evaluation environment.
#' @return TRUE invisibly.
#' @noRd
.erpm_apply_term_spec_checks <- function(spec, fname, args_list, env_eval) {
  deps <- spec$deps

  if (!is.null(deps) && isTRUE(deps$requires_bipartite)) {
    nw <- try(get("nw", envir = env_eval), silent = TRUE)
    if (!inherits(nw, "try-error")) {
      bip <- try(network::get.network.attribute(nw, "bipartite"), silent = TRUE)
      if (inherits(bip, "try-error") || is.null(bip) || is.na(bip)) {
        stop(sprintf("%s(): requires a bipartite `nw` in evaluation environment.", fname))
      }
    }
  }

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

  if (is.function(spec$validate)) {
    spec$validate(args_list, env_eval = env_eval)
  }

  invisible(TRUE)
}

#' Translate a single ERPM term into an \pkg{ergm} term
#'
#' @param term_call Call or symbol representing one RHS term.
#' @param rename_map Optional rename map.
#' @param wrap_proj1 Terms to wrap with \code{Proj1(~ ...)}.
#' @param wrap_B Terms to wrap with \code{B(~ ..., form = "nonzero")}.
#' @param env_eval Evaluation environment.
#' @return Call or symbol.
#' @noRd
.erpm_translate_one_term <- function(term_call,
                                     rename_map,
                                     wrap_proj1 = character(),
                                     wrap_B     = character(),
                                     env_eval   = parent.frame()) {
  if (is.symbol(term_call)) term_call <- as.call(list(term_call))
  if (!is.call(term_call)) return(term_call)

  fun_sym   <- term_call[[1L]]
  args_list <- as.list(term_call)[-1L]
  fname <- if (is.symbol(fun_sym)) as.character(fun_sym) else deparse(fun_sym)[1L]

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

    if (is.list(tr_res) && inherits(tr_res, "erpm_tr")) {
      out_call <- tr_res$call
      attr(out_call, "erpm_meta") <- tr_res$meta
      return(out_call)
    }

    return(tr_res)
  }

  fname2 <- if (fname %in% names(rename_map)) rename_map[[fname]] else fname
  out <- as.call(c(as.name(fname2), args_list))

  if (fname2 %in% wrap_B)     out <- call("B", call("~", out), form = "nonzero")
  if (fname2 %in% wrap_proj1) out <- call("Proj1", call("~", out))

  out
}

# ==============================================================================
# RHS validation / pipeline
# ==============================================================================

#' Validate one translated term
#'
#' @param term_call Translated term call.
#' @param env_eval Evaluation environment.
#' @return Unchanged term call.
#' @noRd
.erpm_validate_translated_term <- function(term_call, env_eval) {
  if (is.symbol(term_call)) term_call <- as.call(list(term_call))
  if (!is.call(term_call)) return(term_call)

  f <- term_call[[1L]]

  if (is.symbol(f) && identical(f, as.name("b2degrange"))) {
    al <- as.pairlist(as.list(term_call)[-1L])
    from <- eval(al$from, envir = env_eval)
    to   <- if (is.null(al$to)) Inf else eval(al$to, envir = env_eval)

    if (is.finite(from) && from < 0L) {
      stop(sprintf("groups(from|k): 'from'/'k' must be >= 0. Got: %d", from))
    }

    if (!(is.infinite(to) || (is.finite(to) && to > from))) {
      to_str <- if (is.infinite(to)) "Inf" else as.character(to)
      stop(sprintf(
        "groups(from,to): requires 'from' < 'to'. Got: from=%d, to=%s",
        from, to_str
      ))
    }
  }

  term_call
}

#' Translate an RHS expression through the strict pipeline
#'
#' @param rhs_expr RHS expression.
#' @param rename_map Term rename map.
#' @param wrap_proj1 Terms to wrap with Proj1.
#' @param wrap_B Terms to wrap with B.
#' @param env_eval Evaluation environment.
#' @return Reconstructed RHS expression.
#' @noRd
.erpm_translate_rhs_pipeline <- function(rhs_expr,
                                         rename_map,
                                         wrap_proj1,
                                         wrap_B,
                                         env_eval) {
  rhs_terms <- .erpm_split_sum_terms(rhs_expr)

  translated <- lapply(
    rhs_terms,
    .erpm_translate_one_term,
    rename_map = rename_map,
    wrap_proj1 = wrap_proj1,
    wrap_B     = wrap_B,
    env_eval   = env_eval
  )

  translated <- lapply(
    translated,
    .erpm_validate_translated_term,
    env_eval = env_eval
  )

  if (length(translated) == 1L) {
    translated[[1L]]
  } else {
    Reduce(function(x, y) call("+", x, y), translated)
  }
}

#' Translate a RHS expression only
#'
#' @param rhs_expr RHS expression.
#' @param env_eval Evaluation environment.
#' @param effect_rename_map Optional rename map.
#' @param wrap_with_proj1 Optional Proj1 wrappers.
#' @param wrap_with_B Optional B wrappers.
#' @return Translated RHS expression.
#' @noRd
.erpm_translate_rhs_expr <- function(rhs_expr,
                                     env_eval,
                                     effect_rename_map = c(),
                                     wrap_with_proj1   = c(),
                                     wrap_with_B       = c()) {
  .erpm_translate_rhs_pipeline(
    rhs_expr   = rhs_expr,
    rename_map = effect_rename_map,
    wrap_proj1 = wrap_with_proj1,
    wrap_B     = wrap_with_B,
    env_eval   = env_eval
  )
}

#' Translate the RHS of a working formula
#'
#' @param formula Working formula.
#' @param env_eval Evaluation environment.
#' @param effect_rename_map Optional rename map.
#' @param wrap_with_proj1 Optional Proj1 wrappers.
#' @param wrap_with_B Optional B wrappers.
#' @param verbose Verbosity flag.
#' @param user_formula_str Original user formula as one line.
#' @return Updated formula.
#' @noRd
.erpm_translate_formula_rhs <- function(formula,
                                        env_eval,
                                        effect_rename_map,
                                        wrap_with_proj1,
                                        wrap_with_B,
                                        verbose,
                                        user_formula_str) {
  rhs_expr <- formula[[3L]]

  new_rhs <- tryCatch(
    .erpm_translate_rhs_pipeline(
      rhs_expr   = rhs_expr,
      rename_map = effect_rename_map,
      wrap_proj1 = wrap_with_proj1,
      wrap_B     = wrap_with_B,
      env_eval   = env_eval
    ),
    error = function(e) {
      if (isTRUE(verbose)) {
        cat(sprintf("[ERPM] call initial : erpm(%s)\n", user_formula_str))
        cat("\t error during translation: ", conditionMessage(e), "\n", sep = "")
      }
      stop(e)
    }
  )

  formula[[3L]] <- new_rhs
  environment(formula) <- env_eval
  formula
}