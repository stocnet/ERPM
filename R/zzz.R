# .onLoad <- function(libname, pkgname) {
#   # Charge le namespace + la DLL d'ergm avant tout symbole MH/_*
#   requireNamespace("ergm", quietly = TRUE)
#   # Enregistre le proposal si l'API est dispo (évite les erreurs au chargement)
#   if (is.function(getNamespaceExports("ergm")[["ergm_proposal_table"]]) ||
#       exists("ergm_proposal_table", where = asNamespace("ergm"), inherits = FALSE)) {
#     try(.RegisterProposals(), silent = TRUE)
#   }
# }

# #' @importFrom ergm ergm_proposal_table
# .RegisterProposals <- function() {
#   # Classe "c", référence Bernoulli, contrainte &b1part, prior 0, nom "random", tag "B1Part"
#   ergm_proposal_table("c", "Bernoulli", "&b1part", 0, "random", "B1Part")
# }

# .onUnload <- function(libpath) {
#   library.dynam.unload("ERPM", libpath)
# }

.onLoad <- function(libname, pkgname) {
  #requireNamespace("ergm")  # charge le namespace + sa DLL
  .RegisterProposals()
}

#' @importFrom ergm ergm_proposal_table
.RegisterProposals <- function() {
  ergm_proposal_table("c", "Bernoulli", "&b1part",  0, "random", "B1Part")
}

.onUnload <- function(libpath) {
  library.dynam.unload("ERPM", libpath)
}
