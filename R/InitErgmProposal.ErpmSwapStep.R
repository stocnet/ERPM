# ==============================================================================
# File    : R/InitErgmProposal.ErpmSwapStep.R
# Purpose : Register the ERPM MCMC proposal 'ErpmSwapStep' for partition moves
#           under the ~ b1part constraint (bipartite actor-group network).
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM proposal: ErpmSwapStep
#'
#' @name InitErgmProposal.ErpmSwapStep
#' @aliases ErpmSwapStep
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' ErpmSwapStep is a Metropolis-Hastings proposal used in ERPM when a partition
#' is represented as a bipartite actor-group membership network under the
#' constraint ~ b1part.
#'
#' The bipartite network is interpreted as:
#' \itemize{
#'   \item mode 1: actors (vertices 1..n, where n = nw %n% 'bipartite');
#'   \item mode 2: groups (vertices n+1..N), including possible empty groups.
#' }
#'
#' The move selects two distinct actors i and j and swaps their group
#' memberships. In the bipartite encoding, this corresponds to a 4-toggle move:
#' removing edges (i, g_i) and (j, g_j), and adding (i, g_j) and (j, g_i).
#'
#' Since actor pairs are sampled uniformly and the swap is constructed
#' deterministically, the proposal is symmetric. The Hastings correction
#' is therefore zero (logratio = 0 in the C implementation).
#'
#' @details
#' This function only registers the proposal with ergm. The actual move logic
#' is implemented in C in MH_ErpmSwapStep.
#'
#' ErpmSwapStep preserves the size of every group. As a consequence, used
#' alone, it is not ergodic on the full partition space: it only explores
#' states sharing the same group-size vector.
#'
#' @return
#' A list describing the proposal for ergm, including:
#' \itemize{
#'   \item name: the compiled proposal name (must match the C symbol suffix);
#'   \item pkgname: the package providing the C code;
#'   \item inputs, iinputs: optional parameters (unused here).
#' }
#'
#' @keywords ERGM MCMC proposal ERPM b1part swap
#' @export
InitErgmProposal.ErpmSwapStep <- function(nw, ...){
  list(
    name    = "ErpmSwapStep",  # Must match MH_ErpmSwapStep in C.
    pkgname = "ERPM",
    inputs  = NULL,
    iinputs = NULL
  )
}