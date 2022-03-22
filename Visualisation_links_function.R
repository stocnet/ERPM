#########################################
## CREATION DU PACKAGE ERPM
## Derniere modif : 16/01/22
#########################################
## Notes :
#########################################


# Pour pouvoir voir les liens entres fonctions.

devtools::install_github("datastorm-open/DependenciesGraphs")
library('erpm')

require(DependenciesGraphs)

dep <- envirDependencies('package:erpm')
plot(dep)

# Verfif help doc

help(estimate_ERPM)
