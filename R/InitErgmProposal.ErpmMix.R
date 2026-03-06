# ==============================================================================
# File    : R/InitErgmProposal.ErpmMix.R
# Purpose : Register the ERPM MCMC proposal 'ErpmMix' for ~b1part partitions.
# ==============================================================================

#' ERGM proposal: ErpmMix
#'
#' @name InitErgmProposal.ErpmMix
#' @aliases ErpmMix
#'
#' @description
#' ErpmMix is a Metropolis-Hastings proposal that draws one partition move type
#' at each MCMC iteration and then applies the corresponding ERPM step.
#'
#' Supported moves:
#' - "toggle" : ErpmToggleStep (2 toggles)
#' - "swap"   : ErpmSwapStep   (4 toggles)
#' - "merge"  : ErpmMergeStep  (merge one non-empty group into another)
#' - "split"  : ErpmSplitStep  (split a non-empty group toward an empty one)
#'
#' @details
#' User-facing arguments:
#' - moves   : character vector in {"toggle","swap","merge","split"}
#' - weights : positive numeric vector of same length as moves
#'
#' Important ergm note:
#' - When this function is called through ergm's `.select("ErpmMix")`, the object
#'   passed as `arguments` usually contains meta-entries such as `constraints`
#'   and `reference`, while the user payload from `MCMC.prop.args[[i]]` is stored
#'   in the first unnamed element: `arguments[[1]]`.
#'
#' Decoding policy:
#' - If user args are missing/empty, use the default move set
#'   moves=c("toggle","swap","merge","split") with weights=c(2,1,0,0).
#' - If `moves` is missing/empty, use defaults.
#' - If `moves` contains unknown entries, do not error: fall back to the canonical
#'   default mix returned to C, namely toggle:2 and swap:1.
#' - If `moves` contains duplicates, treat it as invalid input and fall back to the
#'   same canonical default mix.
#' - If `weights` is missing, use default weights when `moves` is exactly the default
#'   move set; otherwise use a vector of 1s aligned with the selected moves.
#' - If `weights` is invalid (non-numeric, wrong length, non-finite, <=0), fall back
#'   to the canonical default mix.
#'
#' Packing convention to C:
#' - iinputs = c(K, move_codes...)
#' - inputs  = c(weights...)
#' where move codes are stable integers:
#'   toggle = 1, swap = 2, merge = 3, split = 4
#'
#' Notes:
#' - This R initializer only packs the user-requested move mixture.
#' - Feasibility is handled later by the C proposal code. In particular, the move
#'   draw itself stays state-independent; if a selected move is impossible in the
#'   current network state, the C backend falls back to TOGGLE.
#' - As a consequence, observed move frequencies along the chain do not have to
#'   match the nominal user weights exactly, especially when merge/split are often
#'   infeasible on the visited states.
#'
#' @param arguments A list of proposal arguments (may be empty).
#' @param nw A \pkg{network} object (unused here; required by ergm API).
#'
#' @return A list describing the compiled proposal for \pkg{ergm}.
#' @export
InitErgmProposal.ErpmMix <- function(arguments, nw) {
  # ---------------------------------------------------------------------------
  # Debug helpers
  # ---------------------------------------------------------------------------
  dbg    <- isTRUE(getOption("Proposal.ErpmMix.debug", FALSE))
  dbgcat <- function(...) if (dbg) cat("[InitProposal][ErpmMix][DEBUG] ", ..., "\n", sep = "")

  .arg_names_chr <- function(x) {
    nm <- names(x)
    if (is.null(nm) || length(nm) == 0L) return("<none>")
    nm[nm == ""] <- "<unnamed>"
    paste(nm, collapse = ", ")
  }

  .fmt_dbg_value <- function(x) {
    if (is.null(x)) return("<NULL>")
    if (length(x) == 0L) return("<empty>")
    paste(utils::capture.output(dput(x)), collapse = "")
  }

  # Keep 'nw' for API compatibility (unused).
  if (!is.null(nw)) { }  # no-op

  # -----------------------------
  # Defaults
  # -----------------------------
  default_moves   <- c("toggle", "swap", "merge", "split")
  default_weights <- c(2, 1, 0, 0)

  .fallback_default <- function(reason = "unspecified") {
    dbgcat("fallback                = canonical default mix (toggle:2, swap:1)")
    dbgcat("fallback reason         = ", reason)

    codes   <- c(1L, 2L)      # toggle, swap
    weights <- c(2, 1)        # canonical default mix 2:1

    list(
      name    = "ErpmMix",
      pkgname = "ERPM",
      inputs  = as.numeric(weights),
      iinputs = as.integer(c(length(codes), codes))
    )
  }

  # -----------------------------
  # Normalize arguments container
  # -----------------------------
  if (is.null(arguments)) arguments <- list()
  if (!is.list(arguments)) return(.fallback_default("`arguments` is not a list"))

  # With ergm's .select() wrapper, user payload is typically stored in arguments[[1]].
  user_args <- list()
  if (length(arguments) >= 1L && is.list(arguments[[1L]])) {
    user_args <- arguments[[1L]]
  }

  dbgcat("raw argument names      = {", .arg_names_chr(arguments), "}")

  if (length(arguments) >= 1L && is.list(arguments[[1L]])) {
    dbgcat("user payload source     = arguments[[1]]")
    dbgcat("user payload names      = {", .arg_names_chr(arguments[[1L]]), "}")
    dbgcat("user payload$moves      = ", .fmt_dbg_value(arguments[[1L]]$moves))
    dbgcat("user payload$weights    = ", .fmt_dbg_value(arguments[[1L]]$weights))
  } else {
    dbgcat("user payload source     = <missing or not a list>")
  }

  # Extract user fields
  moves   <- user_args$moves
  weights <- user_args$weights

  # -----------------------------
  # Decode moves
  # -----------------------------
  if (is.null(moves) || length(moves) < 1L) {
    moves <- default_moves
  } else {
    if (!is.character(moves)) return(.fallback_default("`moves` is not a character vector"))
    moves <- tolower(trimws(moves))
    moves <- moves[nzchar(moves)]
    if (length(moves) < 1L) return(.fallback_default("`moves` is empty after trimming"))
  }

  allowed <- c("toggle", "swap", "merge", "split")

  # Unknown move names are treated as invalid user input: keep behavior robust
  # and return the canonical default mix instead of partially decoding.
  if (any(!moves %in% allowed)) return(.fallback_default("unsupported move name"))

  # Duplicated move names are treated the same way. The intent is ambiguous and
  # we prefer a predictable fallback over ad hoc aggregation.
  if (any(duplicated(moves))) return(.fallback_default("duplicated move name"))

  # -----------------------------
  # Decode weights
  # -----------------------------
  if (is.null(weights)) {
    # For the full default move set, preserve the intended baseline:
    # toggle=2, swap=1, merge=0, split=0.
    # For any custom move subset, default to equal user-level attempt weights.
    if (identical(moves, default_moves)) weights <- default_weights
    else weights <- rep(1, length(moves))
  }

  if (!is.numeric(weights) || length(weights) != length(moves)) {
    return(.fallback_default("invalid weight vector length/type"))
  }
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    return(.fallback_default("weights must be finite and > 0"))
  }

  move_code <- function(m) switch(m,
    toggle = 1L,
    swap   = 2L,
    merge  = 3L,
    split  = 4L
  )

  codes <- vapply(moves, move_code, integer(1))
  K <- length(codes)

  dbgcat("decoded moves           = ", paste(moves, collapse = ", "))
  dbgcat("decoded weights         = ", paste(format(weights), collapse = ", "))
  dbgcat("packed iinputs          = ", paste(c(K, codes), collapse = ", "))

  list(
    name    = "ErpmMix",   # Must match MH_ErpmMix in C
    pkgname = "ERPM",
    inputs  = as.numeric(weights),
    iinputs = as.integer(c(K, codes))
  )
}


# # ==============================================================================
# # File    : R/InitErgmProposal.ErpmMix.R
# # Purpose : Register the ERPM MCMC proposal 'ErpmMix' for ~b1part partitions.
# # ==============================================================================

# #' ERGM proposal: ErpmMix
# #'
# #' @name InitErgmProposal.ErpmMix
# #' @aliases ErpmMix
# #'
# #' @description
# #' ErpmMix is a Metropolis-Hastings proposal that draws one partition move type
# #' at each MCMC iteration and then applies the corresponding ERPM step.
# #'
# #' Supported moves:
# #' - "toggle" : ErpmToggleStep (2 toggles)
# #' - "swap"   : ErpmSwapStep   (4 toggles)
# #' - "merge"  : ErpmMergeStep  (merge one non-empty group into another)
# #' - "split"  : ErpmSplitStep  (split a non-empty group toward an empty one)
# #'
# #' @details
# #' User-facing arguments:
# #' - moves   : character vector in {"toggle","swap","merge","split"}
# #' - weights : positive numeric vector of same length as moves
# #'
# #' Important ergm note:
# #' - When this function is called through ergm's `.select("ErpmMix")`, the object
# #'   passed as `arguments` usually contains meta-entries such as `constraints`
# #'   and `reference`, while the user payload from `MCMC.prop.args[[i]]` is stored
# #'   in the first unnamed element: `arguments[[1]]`.
# #'
# #' Decoding policy:
# #' - If user args are missing/empty, use the default move set
# #'   moves=c("toggle","swap","merge","split") with weights=c(2,1,0,0).
# #' - If `moves` is missing/empty, use defaults.
# #' - If `moves` contains unknown entries, do not error: fall back to the canonical
# #'   default mix returned to C, namely toggle:2 and swap:1.
# #' - If `moves` contains duplicates, treat it as invalid input and fall back to the
# #'   same canonical default mix.
# #' - If `weights` is missing, use default weights when `moves` is exactly the default
# #'   move set; otherwise use a vector of 1s aligned with the selected moves.
# #' - If `weights` is invalid (non-numeric, wrong length, non-finite, <=0), fall back
# #'   to the canonical default mix.
# #'
# #' Packing convention to C:
# #' - iinputs = c(K, move_codes...)
# #' - inputs  = c(weights...)
# #' where move codes are stable integers:
# #'   toggle = 1, swap = 2, merge = 3, split = 4
# #'
# #' Notes:
# #' - This R initializer only packs the user-requested move mixture.
# #' - Feasibility is handled later by the C proposal code. In particular, the move
# #'   draw itself stays state-independent; if a selected move is impossible in the
# #'   current network state, the C backend falls back to TOGGLE.
# #' - As a consequence, observed move frequencies along the chain do not have to
# #'   match the nominal user weights exactly, especially when merge/split are often
# #'   infeasible on the visited states.
# #'
# #' @param arguments A list of proposal arguments (may be empty).
# #' @param nw A \pkg{network} object (unused here; required by ergm API).
# #'
# #' @return A list describing the compiled proposal for \pkg{ergm}.
# #' @export
# InitErgmProposal.ErpmMix <- function(arguments, nw) {
#   # ---------------------------------------------------------------------------
#   # Debug helpers
#   # ---------------------------------------------------------------------------
#   dbg    <- isTRUE(getOption("Proposal.ErpmMix.debug", FALSE))
#   dbgcat <- function(...) if (dbg) cat("[InitProposal][ErpmMix][DEBUG] ", ..., "\n", sep = "")

#   # Keep 'nw' for API compatibility (unused).
#   if (!is.null(nw)) { }  # no-op

#   # -----------------------------
#   # Defaults
#   # -----------------------------
#   default_moves   <- c("toggle", "swap", "merge", "split")
#   default_weights <- c(2, 1, 0, 0)

#   .fallback_default <- function() {
#     codes   <- c(1L, 2L)      # toggle, swap
#     weights <- c(2, 1)        # canonical default mix 2:1

#     list(
#       name    = "ErpmMix",
#       pkgname = "ERPM",
#       inputs  = as.numeric(weights),
#       iinputs = as.integer(c(length(codes), codes))
#     )
#   }

#   # -----------------------------
#   # Normalize arguments container
#   # -----------------------------
#   if (is.null(arguments)) arguments <- list()
#   if (!is.list(arguments)) return(.fallback_default())

#   # With ergm's .select() wrapper, user payload is typically stored in arguments[[1]].
#   user_args <- list()
#   if (length(arguments) >= 1L && is.list(arguments[[1L]])) {
#     user_args <- arguments[[1L]]
#   }

#   dbgcat("names(arguments) = {", paste(names(arguments), collapse = ", "), "}")
#   dbgcat("user_args (arguments[[1]]) = ", paste(capture.output(str(user_args)), collapse = " "))

#   # Extract user fields
#   moves   <- user_args$moves
#   weights <- user_args$weights

#   # -----------------------------
#   # Decode moves
#   # -----------------------------
#   if (is.null(moves) || length(moves) < 1L) {
#     moves <- default_moves
#   } else {
#     if (!is.character(moves)) return(.fallback_default())
#     moves <- tolower(trimws(moves))
#     moves <- moves[nzchar(moves)]
#     if (length(moves) < 1L) return(.fallback_default())
#   }

#   allowed <- c("toggle", "swap", "merge", "split")

#   # Unknown move names are treated as invalid user input: keep behavior robust
#   # and return the canonical default mix instead of partially decoding.
#   if (any(!moves %in% allowed)) return(.fallback_default())

#   # Duplicated move names are treated the same way. The intent is ambiguous and
#   # we prefer a predictable fallback over ad hoc aggregation.
#   if (any(duplicated(moves))) return(.fallback_default())

#   # -----------------------------
#   # Decode weights
#   # -----------------------------
#   if (is.null(weights)) {
#     # For the full default move set, preserve the intended baseline:
#     # toggle=2, swap=1, merge=0, split=0.
#     # For any custom move subset, default to equal user-level attempt weights.
#     if (identical(moves, default_moves)) weights <- default_weights
#     else weights <- rep(1, length(moves))
#   }

#   if (!is.numeric(weights) || length(weights) != length(moves)) return(.fallback_default())
#   if (any(!is.finite(weights)) || any(weights <= 0)) return(.fallback_default())

#   move_code <- function(m) switch(m,
#     toggle = 1L,
#     swap   = 2L,
#     merge  = 3L,
#     split  = 4L
#   )

#   codes <- vapply(moves, move_code, integer(1))
#   K <- length(codes)

#   dbgcat("decoded moves   = ", paste(moves, collapse = ", "))
#   dbgcat("decoded weights = ", paste(format(weights), collapse = ", "))
#   dbgcat("packed iinputs  = ", paste(c(K, codes), collapse = ", "))

#   list(
#     name    = "ErpmMix",   # Must match MH_ErpmMix in C
#     pkgname = "ERPM",
#     inputs  = as.numeric(weights),
#     iinputs = as.integer(c(K, codes))
#   )
# }