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
#' Packing convention to C:
#' - iinputs = c(K, move_codes...)
#' - inputs  = c(weights...)
#' where move_codes are stable integers:
#'   toggle = 1, swap = 2
#'
#' @param arguments A list of proposal arguments (may be empty). Supported keys:
#'   - moves: character vector in {"toggle","swap"}
#'   - weights: positive numeric vector of same length as moves
#' @param nw A \pkg{network} object (unused here; required by ergm API).
#'
#' @return A list describing the compiled proposal for \pkg{ergm}.
#' @export
InitErgmProposal.ErpmMix <- function(arguments, nw) {
  # Keep 'nw' for API compatibility (unused).
  if (!is.null(nw)) { }  # no-op, avoids linters

  if (is.null(arguments)) arguments <- list()
  if (!is.list(arguments))
    stop("[ERPM] ErpmMix: `arguments` must be a list (ergm internal).", call. = FALSE)

  # Defaults used when `.select("ErpmMix")` provides no arguments.
  moves   <- arguments$moves
  weights <- arguments$weights

  if (is.null(moves))   moves   <- c("toggle", "swap")
  if (is.null(weights)) weights <- c(2, 1)

  if (!is.character(moves) || length(moves) < 1L)
    stop("[ERPM] ErpmMix: `moves` must be a non-empty character vector.", call. = FALSE)

  moves <- tolower(trimws(moves))

  allowed <- c("toggle", "swap")
  if (any(!moves %in% allowed)) {
    bad <- unique(moves[!moves %in% allowed])
    stop("[ERPM] ErpmMix: unsupported moves: ", paste(bad, collapse = ", "),
         ". Allowed: ", paste(allowed, collapse = ", "), call. = FALSE)
  }

  if (!is.numeric(weights) || length(weights) != length(moves))
    stop("[ERPM] ErpmMix: `weights` must be numeric with same length as `moves`.", call. = FALSE)

  if (any(!is.finite(weights)) || any(weights <= 0))
    stop("[ERPM] ErpmMix: all `weights` must be finite and > 0.", call. = FALSE)

  move_code <- function(m) switch(m,
    toggle = 1L,
    swap   = 2L
  )

  codes <- vapply(moves, move_code, integer(1))
  K <- length(codes)

  list(
    name    = "ErpmMix",  # Must match MH_ErpmMix in C
    pkgname = "ERPM",
    inputs  = as.numeric(weights),
    iinputs = as.integer(c(K, codes))
  )
}