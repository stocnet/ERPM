.onLoad <- function(libname, pkgname) {
  .RegisterProposals()
}

#' @importFrom ergm ergm_proposal_table
.RegisterProposals <- function() {
  ergm_proposal_table("c", "Bernoulli", "&b1part",  0, "random", "B1Part")
}

.onUnload <- function(libpath) {
  library.dynam.unload("ERPM", libpath)
}
