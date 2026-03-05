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
#' ErpmMix is a Metropolis-Hastings proposal that randomly selects a move type
#' at each MCMC iteration and then applies the corresponding ERPM partition step.
#'
#' Supported moves (for now):
#' - "toggle" : ErpmToggleStep (2 toggles)
#' - "swap"   : ErpmSwapStep   (4 toggles)
#'
#' @details
#' User-facing arguments:
#' - moves   : character vector in {"toggle","swap"}
#' - weights : positive numeric vector of same length as moves
#'
#' Important ergm note :
#' - When this function is called through ergm's `.select("ErpmMix")`, the object
#'   passed as `arguments` typically contains meta-entries such as `constraints`
#'   and `reference`, and the user payload from `MCMC.prop.args[[i]]` is stored
#'   in the first unnamed element: `arguments[[1]]`.
#'
#' Decoding policy :
#' - If user args are missing/empty, use the default mix: moves=c("toggle","swap"),
#'   weights=c(2,1).
#' - If `moves` is missing/empty, use defaults.
#' - If `moves` contains unknown entries, DO NOT error: fall back to defaults.
#' - If `weights` is missing, use defaults (aligned with the selected moves).
#' - If `weights` is invalid (non-numeric, wrong length, non-finite, <=0), fall back
#'   to defaults.
#'
#' Packing convention to C:
#' - iinputs = c(K, move_codes...)
#' - inputs  = c(weights...)
#' where move_codes are stable integers:
#'   toggle = 1, swap = 2
#'
#' Note:
#' - Only toggle/swap are supported for now; additional moves can be added later
#'   by extending `allowed` + `move_code`.
#'
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

  # Keep 'nw' for API compatibility (unused).
  if (!is.null(nw)) { }  # no-op

  # -----------------------------
  # Defaults 
  # -----------------------------
  default_moves   <- c("toggle", "swap")
  default_weights <- c(2, 1)

  .fallback_default <- function() {
    codes <- c(1L, 2L)
    list(
      name    = "ErpmMix",
      pkgname = "ERPM",
      inputs  = as.numeric(default_weights),
      iinputs = as.integer(c(length(codes), codes))
    )
  }

  # -----------------------------
  # Normalize arguments container
  # -----------------------------
  if (is.null(arguments)) arguments <- list()
  if (!is.list(arguments)) return(.fallback_default())

  # ergm meta wrapper: user payload is typically in arguments[[1]].
  user_args <- list()
  if (length(arguments) >= 1L && is.list(arguments[[1L]])) {
    user_args <- arguments[[1L]]
  }

  dbgcat("names(arguments) = {", paste(names(arguments), collapse = ", "), "}")
  dbgcat("user_args (arguments[[1]]) = ", paste(capture.output(str(user_args)), collapse = " "))

  # Extract user fields
  moves   <- user_args$moves
  weights <- user_args$weights

  # -----------------------------
  # Decode moves
  # -----------------------------
  if (is.null(moves) || length(moves) < 1L) {
    moves <- default_moves
  } else {
    if (!is.character(moves)) return(.fallback_default())
    moves <- tolower(trimws(moves))
    moves <- moves[nzchar(moves)]
    if (length(moves) < 1L) return(.fallback_default())
  }

  allowed <- c("toggle", "swap")

  # If any unknown move is requested: do NOT partially accept; fall back to canonical.
  if (any(!moves %in% allowed)) return(.fallback_default())

  # If user repeats a move, treat as mistake and default.
  if (any(duplicated(moves))) return(.fallback_default())

  # -----------------------------
  # Decode weights
  # -----------------------------
  if (is.null(weights)) {
    # If moves are exactly defaults, keep default weights; else use 1s.
    if (identical(moves, default_moves)) weights <- default_weights
    else weights <- rep(1, length(moves))
  }

  if (!is.numeric(weights) || length(weights) != length(moves)) return(.fallback_default())
  if (any(!is.finite(weights)) || any(weights <= 0)) return(.fallback_default())

  move_code <- function(m) switch(m,
    toggle = 1L,
    swap   = 2L
  )

  codes <- vapply(moves, move_code, integer(1))
  K <- length(codes)

  dbgcat("decoded moves   = ", paste(moves, collapse = ", "))
  dbgcat("decoded weights = ", paste(format(weights), collapse = ", "))
  dbgcat("packed iinputs  = ", paste(c(K, codes), collapse = ", "))

  list(
    name    = "ErpmMix",   # Must match MH_ErpmMix in C
    pkgname = "ERPM",
    inputs  = as.numeric(weights),
    iinputs = as.integer(c(K, codes))
  )
}