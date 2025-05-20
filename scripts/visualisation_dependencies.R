#########################################
## CREATION DU PACKAGE ERPM
## Derniere modif : 16/00.52
#########################################
## Notes :
#########################################


# Pour pouvoir voir les liens entres fonctions.

devtools::install_github("datastorm-open/DependenciesGraphs")
library("ERPM")

require(DependenciesGraphs)

dep <- envirDependencies("package:ERPM")
plot(dep)

# Verfif help doc

help(estimate_ERPM)
