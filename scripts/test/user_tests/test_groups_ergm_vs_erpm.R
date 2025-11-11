
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

# # create nw for ergm
# n <- 6
# nw <- network.initialize(n * 2, dir = FALSE, bip = n)

# # baseline case #1
# nw[cbind(seq_along(partition_mix), partition_mix + n)] <- 1


# 2) Contrôles identiques
ctrl <- control.ergm(
  # init              = "MPLE",
  # CD.nsteps         = 20,        # ou 0 pour sauter CD
  # MCMC.burnin       = 2e5,
  # MCMC.interval     = 5e3,
  # MCMC.samplesize   = 5e3,
  # MCMLE.termination = "Hummel",
  # MCMLE.maxit       = 60
)

set.seed(1)

# 3) Appel {ergm} de référence sur le même réseau
message("[ERGM] Appel {ergm} de référence sur le même réseau")
fit_ergm <- ergm(
  nw2 ~ b2degrange(1, Inf),
  constraints = ~ b1part,
  estimate    = "MLE",
  eval.loglik = TRUE,
  control     = ctrl
)


message("[ERPM] Appel wrapper strictement équivalent")
# 4) Appel wrapper strictement équivalent
fit_erpm <- erpm(
  partition_mix ~ groups,
  estimate    = "MLE",
  eval.loglik = TRUE,
  control     = ctrl
)


# 5) Vérif rapide
print(summary(fit_ergm))
print(summary(fit_erpm))
fit_ergm$coef - fit_erpm$coef