#  File R/InitErgmConstraint.partition.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#' @templateVar name b1part
#' @title Sample over partitioning of Mode 1 nodes identified by edges to Mode 2 nodes
#'
#' @description For bipartite networks, preserve the degree for the
#'   first mode of each vertex at 1, while the second mode is used to
#'   define a partitioning, with an appropriate combinatoric adjustment, if necessary.
#'
#' @usage
#' # b1part
#'
#' @template ergmConstraint-general
#'
#' @concept bipartite
InitErgmConstraint.b1part <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, bipartite = TRUE)
  list(dependence = TRUE, implies = c("b1degrees", "edges"))
}
