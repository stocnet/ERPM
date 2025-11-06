# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cov_ingroup.R
# Objet   : MWE pour l’effet ERPM `cov_ingroup`
# Chaîne  : partition -> biparti -> traduction -> (dry) formule -> summary
#           partition -> biparti -> traduction -> formule -> ergm()
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  library(devtools)  # load_all(".")
  library(network)
  library(ergm)
})

# Charge le package local
devtools::load_all(".")

# Patch {ergm} si présent
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R"); ergm_patch_enable()
}

# ----------------------------- Données fixes -----------------------------------
# Partition: 9 acteurs, groupes {1,1,2,2,3,3,3,4,4} -> tailles = (2,2,3,2)
partition <- c(1,1, 2,2, 3,3,3, 4,4)
labels    <- paste0("A", seq_along(partition))

# Attributs: un numérique (age), un catégoriel (dept)
age  <- c(20,30, 40,50, 10,10,10, 5,5)                 # choix volontairement simples
dept <- c("A","B", "A","C", "A","A","B", "C","A")

nodes <- data.frame(label = labels, age = age, dept = dept, stringsAsFactors = FALSE)

cat("\nPartition : ", paste(partition, collapse=" "), "\n", sep="")
cat("Tailles groupes : ", paste(as.integer(table(partition)), collapse=", "), "\n\n", sep="")

# ---------------------- Réseau biparti explicite (référence) -------------------
nA <- length(partition); G <- max(partition)
inc <- matrix(0L, nrow = nA, ncol = G, dimnames = list(labels, paste0("G",1:G)))
inc[cbind(seq_len(nA), partition)] <- 1L

nw <- network::network(inc, matrix.type = "bipartite", bipartite = nA, directed = FALSE)
network::set.vertex.attribute(nw, "vertex.names", c(labels, colnames(inc)))
network::set.vertex.attribute(nw, "age",  c(nodes$age,  rep(NA_real_, G)))
network::set.vertex.attribute(nw, "dept", c(nodes$dept, rep(NA_character_, G)))

# -------------------------- Cas 1 et Cas 2 : attendus --------------------------
# Définition du terme (implémenté): T = sum_g n_g * sum_{i in g} x_i * 1[n_g ∈ S]
# Cas 1: x_i = age_i, S = {2,3}  -> toutes les tailles de cette partition
# Cas 2: x_i = 1[dept_i == "A"], S = {3}

exp1 <- as.numeric(summary(nw ~ cov_ingroup("age",  size = 2:3), constraints = ~ b1part))
exp2 <- as.numeric(summary(nw ~ cov_ingroup("dept", category = "A", size = 3), constraints = ~ b1part))

# ----------------------- Observés via dry-run erpm() ---------------------------
dry1 <- erpm(partition ~ cov_ingroup("age",  size = 2:3), eval_call = FALSE, verbose = TRUE, nodes = nodes)
obs1 <- as.numeric(summary(dry1[[2]], constraints = ~ b1part))

dry2 <- erpm(partition ~ cov_ingroup("dept", category = "A", size = 3), eval_call = FALSE, verbose = TRUE, nodes = nodes)
obs2 <- as.numeric(summary(dry2[[2]], constraints = ~ b1part))

cat(sprintf("[MWE cov_ingroup] age, S={2,3} : summary=%g | ref=%g\n", obs1, exp1))
stopifnot(all.equal(obs1, exp1, tol = 0))

cat(sprintf("[MWE cov_ingroup] dept=='A', S={3} : summary=%g | ref=%g\n", obs2, exp2))
stopifnot(all.equal(obs2, exp2, tol = 0))

# ------------------------------- Fit CD court ----------------------------------
ctrl_cd <- control.ergm(
  init.method      = "zeros",
  CD.maxit         = 30,
  CD.samplesize    = 20000,
  CD.nsteps        = 3,
  CD.NR.maxit      = 0,
  CD.conv.min.pval = 1e-12,
  MCMLE.maxit      = 0,
  parallel         = 0
)

set.seed(1)
fit <- erpm(
  partition ~ cov_ingroup("age", size = 2:3) + cov_ingroup("dept", category = "A", size = 3),
  estimate    = "CD",
  eval.loglik = TRUE,
  control     = ctrl_cd,
  verbose     = TRUE,
  nodes       = nodes
)
mcmc.diagnostics(fit)
cat("\n--- summary(fit) ---\n")
print(summary(fit))

# Nettoyage patch
on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)