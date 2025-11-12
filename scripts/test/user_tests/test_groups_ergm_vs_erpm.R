
Sys.setenv(LANG="fr_FR.UTF-8")
try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent=TRUE)
options(encoding="UTF-8")

suppressPackageStartupMessages({
  library(devtools)  # load_all(".")
  library(network)   # objet network + attribut 'bipartite'
  library(ergm)      # summary(), ergm(), control.ergm(), etc.
})
devtools::load_all(".")

source("scripts/ergm_patch.R")
ergm_patch_enable()

options(ergm.loglik.warn_dyads=FALSE)

# 1) Même réseau dans les deux cas
partition_mix <- c(1,2,2,3,3,3)

built <- build_bipartite_from_inputs(partition = partition_mix)
nw2   <- built$network  # réseau biparti du wrapper


set.seed(1)

# 2) Appel ergm  sur le même réseau
message("[ERGM] Appel ergm sur le même réseau")
fit_ergm <- ergm(
  nw2 ~ b2degrange(1, Inf),
  constraints = ~ b1part,
  estimate    = "MLE",
  eval.loglik = TRUE
)


message("[ERPM] Appel erpm strictement équivalent")
# 3) Appel erpm  équivalent
fit_erpm <- erpm(
  partition_mix ~ groups,
  estimate    = "MLE",
  eval.loglik = TRUE
)


# 5) Vérif rapide
print(summary(fit_ergm))
print(summary(fit_erpm))
fit_ergm$coef - fit_erpm$coef