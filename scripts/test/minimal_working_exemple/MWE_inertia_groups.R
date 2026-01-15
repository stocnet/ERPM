# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_erpm_long_inertia_groups.R
# Objet   : MWE pour comparer erpm_long() et ergm() sur un cas longitudinal
# Chaîne  : 3 partitions -> réseaux bipartis ->
#           (A) erpm_long(...) -> 3 fits (t=1..3)
#           (B) ergm(...) manuels : t=1..2 (cliques), t=3 (cliques + inertia_groups)
# Notes   :
#   - past_influence = 2 => effet inertiel activé uniquement à t = 3 (car t <= L => off).
#   - La consommation des seeds peut diverger entre pipelines. On compare "à peu près"
#     (coefficients proches, logLik proche), pas forcément identiques bit-à-bit.
# ==============================================================================

Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  library(devtools)
  library(network)
  library(ergm)
})

# ----------------------------------------------------------------------
# Chargement du package local ERPM
# ----------------------------------------------------------------------
devtools::load_all(".")

# Patch ergm si présent
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

cat("=== MWE erpm_long vs ergm : cliques + inertia_groups(past_influence=2) ===\n")

# ==============================================================================
# 1) Données : 3 partitions (un peu plus complexes)
# ==============================================================================

# N = 12, partitions plus riches (splits/merges/mouvements)
P1 <- c(
  1,1,1,     # {1,2,3}
  2,2,       # {4,5}
  3,3,       # {6,7}
  4,4,       # {8,9}
  5,         # {10}
  6,6        # {11,12}
)

# P2: split du groupe 1, merge partiel, déplacement d'un singleton
P2 <- c(
  1,1,       # {1,2} (ex-groupe1)
  2,         # {3}   (ex-groupe1 devient singleton)
  2,2,       # {4,5} reste
  3,         # {6}   (split du groupe3)
  4,         # {7}   (split du groupe3)
  4,4,       # {8,9} reste
  5,5,       # {10,11} (merge {10} avec {11} qui quitte {11,12})
  6          # {12} (devient singleton)
)

# P3: re-merge/swap plus fort pour activer inertia_groups sur des signatures variées
P3 <- c(
  1,         # {1} singleton
  2,2,2,     # {2,3,4} merge (2 vient de 1, 3 vient de 2, 4 vient de 2)
  3,3,       # {5,6}
  4,4,4,     # {7,8,9} merge
  5,         # {10} singleton (ex-merge avec 11 se défait)
  6,6        # {11,12} re-merge
)

parts <- list(P1 = P1, P2 = P2, P3 = P3)

cat("\nPartitions:\n")
for (nm in names(parts)) {
  p <- parts[[nm]]
  tab <- sort(table(p))
  cat(sprintf("  %s : %s\n", nm, paste(p, collapse = " ")))
  cat(sprintf("       groupes=%d | tailles=%s\n",
              length(tab),
              paste(as.integer(tab), collapse = ",")))
}

# Nodes (fixes dans le temps)
set.seed(1)
nodes <- data.frame(
  label = paste0("N", seq_along(P1)),
  age   = sample(20:60, length(P1), TRUE),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2) Helpers : builder + attache inertie (lag1, lag2) pour le pipeline manuel
# ==============================================================================

if (!exists("build_bipartite_from_inputs", mode = "function")) {
  stop("build_bipartite_from_inputs() manquant. Charge ERPM correctement.")
}

.groups_from_partition <- function(p) split(seq_along(p), as.integer(p))

.group_signatures <- function(groups) {
  vapply(groups, function(v) paste(sort(as.integer(v)), collapse = ","), character(1))
}

.attach_inertia_groups_lag <- function(nw_t, lag, p_prev) {
  termname <- "inertia_groups"
  n1 <- as.integer(nw_t %n% "bipartite")
  if (!is.finite(n1) || n1 < 1L) stop("nw_t: bipartite invalide.")
  p_prev <- as.integer(round(p_prev))
  if (length(p_prev) != n1) stop("p_prev: longueur != n1.")

  groups <- .groups_from_partition(p_prev)
  sig    <- .group_signatures(groups)
  sz     <- vapply(groups, length, integer(1))

  obj <- list(
    type       = "group_signature_set",
    lag        = as.integer(lag),
    signatures = as.character(sig),
    sizes      = as.integer(sz)
  )

  attr_name <- paste0("erpm_inertia__", termname, "__lag", lag)
  network::set.network.attribute(nw_t, attr_name, obj)

  nw_t
}

make_nw <- function(part) {
  bld <- build_bipartite_from_inputs(partition = part, nodes = nodes, dyads = list())
  if (is.list(bld) && !is.null(bld$network) && inherits(bld$network, "network")) return(bld$network)
  if (inherits(bld, "network")) return(bld)
  stop("build_bipartite_from_inputs() ne renvoie pas un 'network'.")
}

# Réseaux du pipeline manuel
nw1 <- make_nw(P1)
nw2 <- make_nw(P2)
nw3 <- make_nw(P3)

# past_influence = 2 => pour t=3, on attache lag1 = P2 et lag2 = P1
nw3 <- .attach_inertia_groups_lag(nw3, lag = 1L, p_prev = P2)
nw3 <- .attach_inertia_groups_lag(nw3, lag = 2L, p_prev = P1)

# ==============================================================================
# 3) Panel : effets (pour affichage)
# ==============================================================================

rhs_static <- "cliques(2)"
rhs_inert  <- "inertia_groups(past_influence=2)"

cat("\nRHS:\n")
cat("  t=1 :", rhs_static, "\n")
cat("  t=2 :", rhs_static, "\n")
cat("  t=3 :", paste(rhs_static, "+", rhs_inert), "\n\n")

# ==============================================================================
# 4) (A) erpm_long : 3 fits (sans extraction de calls, sans re-fit)
# ==============================================================================

if (!exists("erpm_long", mode = "function")) {
  stop("erpm_long() manquant.")
}

set.seed(1)
res_long <- erpm_long(
  parts ~ cliques(2) + inertia_groups(past_influence = 2),
  eval.call     = TRUE,
  eval.loglik   = TRUE,
  seed          = 1,  
  verbose       = TRUE,
  debug         = TRUE,
  nodes         = nodes,
  dyads         = list()
)

# Extraction robuste des fits produits par erpm_long()
extract_fits <- function(x) {
  if (is.list(x) && !is.null(x$models)) return(x$models)
  if (is.list(x) && !is.null(x$fits))   return(x$fits)
  if (is.list(x) && !is.null(x$models_fitted)) return(x$models_fitted)
  if (is.list(x) && all(vapply(x, function(z) inherits(z, "ergm"), logical(1)))) return(x)
  stop("Impossible d’extraire la liste des fits depuis erpm_long(). Regarde str(res_long).")
}

fits_long <- extract_fits(res_long)

cat("\n--- erpm_long : structure ---\n")
print(str(res_long, max.level = 2))
cat("\n--- erpm_long : fits détectés ---\n")
cat(sprintf("nfits=%d | classes=%s\n",
            length(fits_long),
            paste(unique(vapply(fits_long, function(z) class(z)[1], character(1))), collapse = ",")))

if (length(fits_long) != 3L) {
  stop(sprintf("Attendu 3 fits via erpm_long(), obtenu %d.", length(fits_long)))
}

# ==============================================================================
# 5) (B) ergm manuels : 3 fits (écrits ici, pas issus de erpm_long)
# ==============================================================================

# ctrl <- control.ergm(
#   MCMLE.maxit     = 3,
#   MCMC.interval   = 1024,
#   MCMC.burnin     = 1024,
#   MCMC.samplesize = 1024,
#   seed            = 1
# )

set.seed(1)
fit_manual_1 <- ergm(
  nw1 ~ cliques(2),
  constraints = ~ b1part,
  eval.loglik = TRUE,
#   control     = NULL,
  verbose     = TRUE
)

set.seed(1)
fit_manual_2 <- ergm(
  nw2 ~ cliques(2),
  constraints = ~ b1part,
  eval.loglik = TRUE,
#   control     = NULL,
  verbose     = TRUE
)

set.seed(1)
fit_manual_3 <- ergm(
  nw3 ~ cliques(2) + inertia_groups(past_influence = 2),
  constraints = ~ b1part,
  eval.loglik = TRUE,
#   control     = NULL,
  verbose     = TRUE
)

fits_manual <- list(fit_manual_1, fit_manual_2, fit_manual_3)

# ==============================================================================
# 6) Comparaisons : coefficients + logLik (tolérances)
# ==============================================================================

coef_comp <- function(a, b) {
  ca <- coef(a)
  cb <- coef(b)
  alln <- union(names(ca), names(cb))
  ca2 <- ca[alln]; cb2 <- cb[alln]
  ca2[is.na(ca2)] <- NA_real_
  cb2[is.na(cb2)] <- NA_real_

  data.frame(
    term = alln,
    coef_long   = as.numeric(ca2),
    coef_manual = as.numeric(cb2),
    diff        = as.numeric(ca2 - cb2),
    row.names   = NULL
  )
}

loglik_safe <- function(fit) {
  out <- try(as.numeric(logLik(fit)), silent = TRUE)
  if (inherits(out, "try-error")) return(NA_real_)
  out
}

cat("\n=== COMPARAISON erpm_long (fits) vs ergm (fits manuels) ===\n")
for (t in 1:3) {
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("[t=%d]\n", t))

  cat("--- summary(erpm_long fit) ---\n")
  print(summary(fits_long[[t]]))

  cat("\n--- summary(ergm manuel) ---\n")
  print(summary(fits_manual[[t]]))

  tab <- coef_comp(fits_long[[t]], fits_manual[[t]])
  cat("\n--- Coefficients (aligned) ---\n")
  print(tab)

  ll_long   <- loglik_safe(fits_long[[t]])
  ll_manual <- loglik_safe(fits_manual[[t]])
  cat(sprintf("\nlogLik: long=%s | manual=%s | diff=%s\n",
              format(ll_long), format(ll_manual), format(ll_long - ll_manual)))

  tol_coef <- 1e-2
  tol_ll   <- 1e-1

  ok_coef <- all(is.na(tab$diff) | abs(tab$diff) <= tol_coef)
  ok_ll   <- is.na(ll_long) || is.na(ll_manual) || abs(ll_long - ll_manual) <= tol_ll

  cat(sprintf("Check soft: coef<=%g ? %s | logLik<=%g ? %s\n",
              tol_coef, if (ok_coef) "OK" else "NO",
              tol_ll, if (ok_ll) "OK" else "NO"))
}

cat("\nMWE erpm_long vs ergm terminé.\n")
on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)