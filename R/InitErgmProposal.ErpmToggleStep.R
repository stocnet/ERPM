# ==============================================================================
# File    : R/InitErgmProposal.ErpmToggleStep.R
# Purpose : Register the MCMC proposal 'ErpmToggleStep' for ERPM partitions
#           represented as bipartite actor → group membership networks.
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM proposal: ErpmToggleStep
#'
#' @name InitErgmProposal.ErpmToggleStep
#' @aliases ErpmToggleStep
#' @note InitErgmProposal.ErpmToggleStep.R
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' ErpmToggleStep defines a Metropolis-Hastings proposal used when sampling
#' partitions encoded as bipartite actor-group networks under the constraint
#' ~ b1part.
#'
#' Conceptually, a toggle step reassigns one actor to a different group.
#' In the bipartite representation, this corresponds to a 2-toggle move:
#' turning off the current membership edge (i, g_old) and turning on
#' (i, g_new).
#'
#' The proposal:
#' - draws an actor uniformly among all actors,
#' - draws a new group uniformly among all groups except the current one,
#'   including empty padded groups.
#'
#' Because both directions are sampled with the same probability,
#' the proposal is symmetric and the Hastings log-ratio is zero.
#'
#' @details
#' This R initializer only declares the compiled proposal to ergm.
#' The actual move logic is implemented in C in MH_ErpmToggleStep.
#'
#' It is meant to be used together with constraints = ~ b1part.
#'
#' @return
#' A list understood by ergm, specifying:
#' - name: the compiled proposal identifier,
#' - pkgname: the package providing the C implementation,
#' - inputs / iinputs: optional parameters (unused here).
#'
#' @keywords ERGM proposal MCMC ERPM b1part partition
#' @md
#' @export
InitErgmProposal.ErpmToggleStep <- function(nw, ...){
  list(
    name    = "ErpmToggleStep",  # Must match MH_ErpmToggleStep in C
    pkgname = "ERPM",
    inputs  = NULL,
    iinputs = NULL
  )
}