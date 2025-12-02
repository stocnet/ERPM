# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_cov_match.R
# Objet   : MWE minimal pour l’effet ERPM `cov_match`
# Chaîne  : partition → biparti → summary(nw) / erpm(partition)
#           avec vérification des normalisations none / by_group / global
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  library(devtools)  # expose erpm(), build_bipartite_from_inputs
  library(ergm)
})

# Charger le package local
devtools::load_all(".")

# Patch {ergm} si disponible
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# ----------------------------- Données fixes -----------------------------------
# Partition simple (10 acteurs)
partition <- c(1, 1,
               2, 2,
               3, 3, 3,
               4, 4, 4)
labels    <- paste0("A", seq_along(partition))

# Attributs : sexe, dept
sexe <- c("F", "H",
          "F", "F",
          "H", "H", "H",
          "F", "H", "H")

dept <- c("RH",  "RH",
          "Info","Info",
          "RH",  "Info","Info",
          "RH",  "RH",  "Info")

nodes <- data.frame(
  label = labels,
  sexe  = sexe,
  dept  = dept,
  stringsAsFactors = FALSE
)

cat("\nPartition :       ", paste(partition, collapse = " "), "\n", sep = "")
cat("Tailles groupes : ", paste(as.integer(table(partition))), "\n\n", sep = "")

# ---------------------- Biparti via wrapper (référence) ------------------------
bld <- build_bipartite_from_inputs(partition = partition, nodes = nodes)
nw  <- bld$network

# -------------------------- Summary de référence -------------------------------
# k = 2, normalisation "none"
ref_k2 <- as.numeric(
  summary(nw ~ cov_match("sexe", clique_size = 2),
          constraints = ~ b1part)
)

# k = 2, normalisation "by_group"
ref_k2_bg <- as.numeric(
  summary(nw ~ cov_match("sexe", clique_size = 2, normalized = "by_group"),
          constraints = ~ b1part)
)

# k = 2, normalisation "global"
ref_k2_glob <- as.numeric(
  summary(nw ~ cov_match("sexe", clique_size = 2, normalized = "global"),
          constraints = ~ b1part)
)

# k = 1, normalisation "by_group" + catégorie ciblée
ref_k1_F_bg <- as.numeric(
  summary(nw ~ cov_match("sexe", clique_size = 1, category = "F",
                         normalized = "by_group"),
          constraints = ~ b1part)
)

# -------------------------- Summary via erpm() ---------------------------------
# k = 2, normalisation "none"
dry_k2 <- erpm(
  partition ~ cov_match("sexe", clique_size = 2),
  eval_call = FALSE,
  verbose   = FALSE,
  nodes     = nodes
)
obs_k2 <- as.numeric(
  summary(dry_k2[[2]], constraints = ~ b1part)
)

# k = 2, normalisation "by_group"
dry_k2_bg <- erpm(
  partition ~ cov_match("sexe", clique_size = 2, normalized = "by_group"),
  eval_call = FALSE,
  verbose   = FALSE,
  nodes     = nodes
)
obs_k2_bg <- as.numeric(
  summary(dry_k2_bg[[2]], constraints = ~ b1part)
)

# k = 2, normalisation "global"
dry_k2_glob <- erpm(
  partition ~ cov_match("sexe", clique_size = 2, normalized = "global"),
  eval_call = FALSE,
  verbose   = FALSE,
  nodes     = nodes
)
obs_k2_glob <- as.numeric(
  summary(dry_k2_glob[[2]], constraints = ~ b1part)
)

# k = 1, normalisation "by_group" + catégorie ciblée
dry_k1_F_bg <- erpm(
  partition ~ cov_match("sexe", clique_size = 1, category = "F",
                        normalized = "by_group"),
  eval_call = FALSE,
  verbose   = FALSE,
  nodes     = nodes
)
obs_k1_F_bg <- as.numeric(
  summary(dry_k1_F_bg[[2]], constraints = ~ b1part)
)

# -------------------------- Vérifications --------------------------------------
cat(sprintf("[summary] cov_match(sexe,k=2)                 : obs=%g | ref=%g\n",
            obs_k2, ref_k2))
stopifnot(all.equal(obs_k2, ref_k2, tol = 0))

cat(sprintf("[summary] cov_match(sexe,k=2,by_group)       : obs=%g | ref=%g\n",
            obs_k2_bg, ref_k2_bg))
stopifnot(all.equal(obs_k2_bg, ref_k2_bg, tol = 0))

cat(sprintf("[summary] cov_match(sexe,k=2,global)         : obs=%g | ref=%g\n",
            obs_k2_glob, ref_k2_glob))
stopifnot(all.equal(obs_k2_glob, ref_k2_glob, tol = 0))

cat(sprintf("[summary] cov_match(sexe==F,k=1,by_group)    : obs=%g | ref=%g\n",
            obs_k1_F_bg, ref_k1_F_bg))
stopifnot(all.equal(obs_k1_F_bg, ref_k1_F_bg, tol = 0))

# ------------------------------- Fit simple ------------------------------------
set.seed(1)
bld2 <- build_bipartite_from_inputs(partition = partition, nodes = nodes)
nw2  <- bld2$network

# Fit de référence (appel direct à ergm)
fit_ref <- ergm(
  nw2 ~ cov_match("sexe", clique_size = 2),
  constraints = ~ b1part,
  eval.loglik = TRUE
)

# Fit via erpm() sur la même statistique non normalisée
set.seed(1)
fit_erpm <- erpm(
  partition ~ cov_match("sexe", clique_size = 2),
  eval.loglik = TRUE,
  verbose     = TRUE,
  nodes       = nodes
)

cat("\n--- summary(fit_ref) ---\n")
print(summary(fit_ref))

cat("\n--- summary(fit_erpm) ---\n")
print(summary(fit_erpm))

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)