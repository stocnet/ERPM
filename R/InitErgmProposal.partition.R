#  File R/InitErgmProposal.partition.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#' @templateVar name B1Part
#' @aliases InitErgmProposal.B1Part
#' @title MHp for b1part constraints
#'
#' @description MHp for `constraints= ~b1part`. At this time, this is
#'   just the [`CondB1Degree`][CondB1Degree-ergmProposal] proposal
#'   with a combinatoric adjustment.
#'
#' @template ergmProposal-general
#' @importFrom ergm ergm_Init_stop
#' @importFrom network is.bipartite
NULL
InitErgmProposal.B1Part <- function(arguments, nw) {
  proposal <- list(name = "B1Part", inputs=NULL)
  if (!is.bipartite(nw)) # Really, this should never trigger, since the InitErgmConstraint function should check.
    ergm_Init_stop("The B1Part proposal function does not work with a non-bipartite network.")

  proposal
}
