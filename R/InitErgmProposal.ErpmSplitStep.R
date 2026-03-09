# ==============================================================================
# File    : R/InitErgmProposal.ErpmSplitStep.R
# Purpose : Register the ERPM MCMC proposal 'ErpmSplitStep' for partition moves
#           under the ~ b1part constraint (bipartite actor-group network).
# Project : ERPM / ERGM extensions
# ==============================================================================

#' ERGM proposal: ErpmSplitStep
#'
#' @name InitErgmProposal.ErpmSplitStep
#' @aliases ErpmSplitStep
#' @author Jérémie Chichignoud - Cub'itech
#'
#' @description
#' ErpmSplitStep is a Metropolis-Hastings proposal used in ERPM when a partition
#' is represented as a bipartite actor-group membership network under the
#' constraint ~ b1part.
#'
#' The move selects a non-empty group g, picks a subset of its actors, selects an
#' empty group vertex g', and moves the selected subset from g to g'. In the
#' bipartite encoding, this corresponds to a multi-toggle move of size 2*|S|:
#' for each actor i in S, remove (i, g) and add (i, g').
#'
#' @details
#' This function only registers the proposal with ergm. The actual move logic
#' is implemented in C in MH_ErpmSplitStep.
#'
#' @return A list describing the proposal for ergm.
#' @keywords ERGM MCMC proposal ERPM b1part split
#' @export
InitErgmProposal.ErpmSplitStep <- function(nw, ...){
  list(
    name    = "ErpmSplitStep",  # Must match MH_ErpmSplitStep in C.
    pkgname = "ERPM",
    inputs  = NULL,
    iinputs = NULL
  )
}