# ======================================================================================
# Fichier : scripts/test/selftests/selftest_dyadcov_retrocompat_erpm_matrix.R
# Objet   : Self-test rétrocompatibilité erpm() : dyads en list (legacy) vs matrice
#           (nouvelle compatibilité) + effets dyadcov / dyadcov_GW
# Exécution: Rscript scripts/test/selftests/selftest_dyadcov_retrocompat_erpm_matrix.R
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule
# --------------------------------------------------------------------------------------
options(ergm.loglik.warn_dyads = FALSE)
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",    quietly = TRUE)) stop("Package 'ergm' requis.")
})

suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

# Patch ERGM optionnel
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# Charger le package et le wrapper ERPM
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Le fichier DESCRIPTION n'existe pas ou devtools n'est pas installé.")
}
if (!exists("erpm", mode = "function")) {
  if (file.exists("R/erpm_wrapper.R")) {
    source("R/erpm_wrapper.R", local = FALSE)
  } else stop("erpm_wrapper.R introuvable.")
}
if (!exists("build_bipartite_from_inputs", mode = "function")) {
  stop("build_bipartite_from_inputs() indisponible. Il doit être exporté par R/erpm_wrapper.R.")
}

# --------------------------------------------------------------------------------------
# Logging local
# --------------------------------------------------------------------------------------
.get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  i <- grep("^--file=", a)
  if (length(i)) return(dirname(normalizePath(sub("^--file=", "", a[i[1]]))))
  fs <- sys.frames()
  ofiles <- vapply(fs, function(f) if (!is.null(f$ofile)) f$ofile else NA_character_, "")
  if (any(!is.na(ofiles))) {
    j <- which.max(nchar(ofiles))
    return(dirname(normalizePath(ofiles[j])))
  }
  normalizePath(getwd())
}

root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_dyadcov_retrocompat_erpm_matrix.log")
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
if (file.exists(log_path)) unlink(log_path, force = TRUE)
con_out <- file(log_path, open = "wt")
con_err <- file(log_path, open = "at")
sink(con_out, split = TRUE)
sink(con_err, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(close(con_err),        silent = TRUE)
  try(sink(),                silent = TRUE)
  try(close(con_out),        silent = TRUE)
  flush.console()
}, add = TRUE)
cat("==> Log:", log_path, "\n")

# ======================================================================================
# Données de test
# ======================================================================================

partitions <- list(
  P1 = c(
    rep(1L, 6),
    rep(2L, 2),
    rep(3L, 5),
    rep(4L, 2)
  ),
  P2 = c(
    rep(1L, 2),
    rep(2L, 5),
    rep(3L, 5),
    rep(4L, 3)
  ),
  P3 = c(
    rep(1L, 3),
    rep(2L, 3),
    rep(3L, 5)
  )
)

# Nodes muets pour satisfaire le builder
.make_nodes_df_for_partition <- function(part) {
  n <- length(part)
  data.frame(
    label = paste0("N", seq_len(n)),
    stringsAsFactors = FALSE
  )
}

# ======================================================================================
# Matrices dyadiques écrites en dur (diagonale = 0)
# ======================================================================================

# ----- Partition P1 (n = 15) -----

Z1_P1 <- matrix(c(
  0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5,
  1.25, 0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25,
  1.5, 1.25, 0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4,
  1.75, 1.5, 1.25, 0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75,
  2, 1.75, 1.5, 1.25, 0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5,
  2.25, 2, 1.75, 1.5, 1.25, 0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25,
  2.5, 2.25, 2, 1.75, 1.5, 1.25, 0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3,
  2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75,
  3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 0, 1.25, 1.5, 1.75, 2, 2.25, 2.5,
  3.25, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 0, 1.25, 1.5, 1.75, 2, 2.25,
  3.5, 3.25, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 0, 1.25, 1.5, 1.75, 2,
  3.75, 3.5, 3.25, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 0, 1.25, 1.5, 1.75,
  4, 3.75, 3.5, 3.25, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 0, 1.25, 1.5,
  4.25, 4, 3.75, 3.5, 3.25, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 0, 1.25,
  4.5, 4.25, 4, 3.75, 3.5, 3.25, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 0
), nrow = 15, ncol = 15, byrow = TRUE)

Z2_P1 <- matrix(c(
  0, 0.5, 0.8, 1.1, 1.4, 1.7, 2, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1, 4.4,
  0.7, 0, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3, 3.3, 3.6, 3.9, 4.2, 4.5,
  0.9, 1.2, 0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7, 4, 4.3, 4.6,
  1.1, 1.4, 1.7, 0, 1.7, 2, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1, 4.4, 4.7,
  1.3, 1.6, 1.9, 2.2, 0, 2.1, 2.4, 2.7, 3, 3.3, 3.6, 3.9, 4.2, 4.5, 4.8,
  1.5, 1.8, 2.1, 2.4, 2.7, 0, 2.5, 2.8, 3.1, 3.4, 3.7, 4, 4.3, 4.6, 4.9,
  1.7, 2, 2.3, 2.6, 2.9, 3.2, 0, 2.9, 3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 5,
  1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7, 0, 3.3, 3.6, 3.9, 4.2, 4.5, 4.8, 5.1,
  2.1, 2.4, 2.7, 3, 3.3, 3.6, 3.9, 4.2, 0, 3.7, 4, 4.3, 4.6, 4.9, 5.2,
  2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 0, 4.1, 4.4, 4.7, 5, 5.3,
  2.5, 2.8, 3.1, 3.4, 3.7, 4, 4.3, 4.6, 4.9, 5.2, 0, 4.5, 4.8, 5.1, 5.4,
  2.7, 3, 3.3, 3.6, 3.9, 4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 0, 5, 5.3, 5.6,
  2.9, 3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 5, 5.3, 5.6, 5.9, 6.2, 0, 5.5, 5.8,
  3.1, 3.4, 3.7, 4, 4.3, 4.6, 4.9, 5.2, 5.5, 5.8, 6.1, 6.4, 6.7, 0, 6,
  3.3, 3.6, 3.9, 4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6, 6.3, 6.6, 6.9, 7.2, 0
), nrow = 15, ncol = 15, byrow = TRUE)

dyads_P1 <- list(Z1 = Z1_P1, Z2 = Z2_P1)

# ----- Partition P2 (n = 15) -----

Z1_P2 <- matrix(c(
  0, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1,
  0.8, 0, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1,
  0.9, 1, 0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2,
  1, 1.1, 1.2, 0, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3,
  1.1, 1.2, 1.3, 1.4, 0, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4,
  1.2, 1.3, 1.4, 1.5, 1.6, 0, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5,
  1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 0, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
  1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
  1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 0, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8,
  1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 0, 2.5, 2.6, 2.7, 2.8, 2.9,
  1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 0, 2.7, 2.8, 2.9, 3,
  1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 0, 2.9, 3, 3.1,
  1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 0, 3.1, 3.2,
  2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 0, 3.3,
  2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 0
), nrow = 15, ncol = 15, byrow = TRUE)

Z2_P2 <- matrix(c(
  0, 0.65, 0.85, 1.1, 1.4, 1.75, 2.15, 2.6, 3.1, 3.65, 4.25, 4.9, 5.6, 6.35, 7.15,
  0.75, 0, 0.9, 1.15, 1.45, 1.8, 2.2, 2.65, 3.15, 3.7, 4.3, 4.95, 5.65, 6.4, 7.2,
  0.95, 1.15, 0, 1.2, 1.5, 1.85, 2.25, 2.7, 3.2, 3.75, 4.35, 5, 5.7, 6.45, 7.25,
  1.15, 1.3, 1.5, 0, 1.55, 1.9, 2.3, 2.75, 3.25, 3.8, 4.4, 5.05, 5.75, 6.5, 7.3,
  1.35, 1.45, 1.65, 1.9, 0, 1.95, 2.35, 2.8, 3.3, 3.85, 4.45, 5.1, 5.8, 6.55, 7.35,
  1.55, 1.6, 1.8, 2.05, 2.35, 0, 2.4, 2.85, 3.35, 3.9, 4.5, 5.15, 5.85, 6.6, 7.4,
  1.75, 1.75, 1.95, 2.2, 2.5, 2.85, 0, 2.9, 3.4, 3.95, 4.55, 5.2, 5.9, 6.65, 7.45,
  1.95, 1.9, 2.1, 2.35, 2.65, 3, 3.4, 0, 3.45, 4, 4.6, 5.25, 5.95, 6.7, 7.5,
  2.15, 2.05, 2.25, 2.5, 2.8, 3.15, 3.55, 4, 0, 4.05, 4.65, 5.3, 6, 6.75, 7.55,
  2.35, 2.2, 2.4, 2.65, 2.95, 3.3, 3.7, 4.15, 4.65, 0, 4.7, 5.35, 6.05, 6.8, 7.6,
  2.55, 2.35, 2.55, 2.8, 3.1, 3.45, 3.85, 4.3, 4.8, 5.35, 0, 5.4, 6.1, 6.85, 7.65,
  2.75, 2.5, 2.7, 2.95, 3.25, 3.6, 4, 4.45, 4.95, 5.5, 6.1, 0, 6.15, 6.9, 7.7,
  2.95, 2.65, 2.85, 3.1, 3.4, 3.75, 4.15, 4.6, 5.1, 5.65, 6.25, 6.9, 0, 6.95, 7.75,
  3.15, 2.8, 3, 3.25, 3.55, 3.9, 4.3, 4.75, 5.25, 5.8, 6.4, 7.05, 7.75, 0, 7.8,
  3.35, 2.95, 3.15, 3.4, 3.7, 4.05, 4.45, 4.9, 5.4, 5.95, 6.55, 7.2, 7.9, 8.65, 0
), nrow = 15, ncol = 15, byrow = TRUE)

dyads_P2 <- list(Z1 = Z1_P2, Z2 = Z2_P2)

# ----- Partition P3 (n = 11) -----
# Z1_P3 : symétrique
# Z2_P3 : asymétrique (Z2[i,j] != Z2[j,i] en général)

Z1_P3 <- matrix(c(
  0, 1.1, 1.4, 1.7, 2, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8,
  1.1, 0, 1.1, 1.4, 1.7, 2, 2.3, 2.6, 2.9, 3.2, 3.5,
  1.4, 1.1, 0, 1.1, 1.4, 1.7, 2, 2.3, 2.6, 2.9, 3.2,
  1.7, 1.4, 1.1, 0, 1.1, 1.4, 1.7, 2, 2.3, 2.6, 2.9,
  2, 1.7, 1.4, 1.1, 0, 1.1, 1.4, 1.7, 2, 2.3, 2.6,
  2.3, 2, 1.7, 1.4, 1.1, 0, 1.1, 1.4, 1.7, 2, 2.3,
  2.6, 2.3, 2, 1.7, 1.4, 1.1, 0, 1.1, 1.4, 1.7, 2,
  2.9, 2.6, 2.3, 2, 1.7, 1.4, 1.1, 0, 1.1, 1.4, 1.7,
  3.2, 2.9, 2.6, 2.3, 2, 1.7, 1.4, 1.1, 0, 1.1, 1.4,
  3.5, 3.2, 2.9, 2.6, 2.3, 2, 1.7, 1.4, 1.1, 0, 1.1,
  3.8, 3.5, 3.2, 2.9, 2.6, 2.3, 2, 1.7, 1.4, 1.1, 0
), nrow = 11, ncol = 11, byrow = TRUE)

Z2_P3 <- matrix(c(
  0, 0.77, 1.33, 2.03, 2.87, 3.85, 4.97, 6.23, 7.63, 9.17, 10.85,
  0.56, 0, 1.47, 2.17, 3.01, 3.99, 5.11, 6.37, 7.77, 9.31, 10.99,
  0.84, 0.98, 0, 2.31, 3.15, 4.13, 5.25, 6.51, 7.91, 9.45, 11.13,
  1.12, 1.4, 1.96, 0, 3.29, 4.27, 5.39, 6.65, 8.05, 9.59, 11.27,
  1.4, 1.82, 2.38, 3.08, 0, 4.41, 5.53, 6.79, 8.19, 9.73, 11.41,
  1.68, 2.24, 2.8, 3.5, 4.34, 0, 5.67, 6.93, 8.33, 9.87, 11.55,
  1.96, 2.66, 3.22, 3.92, 4.76, 5.74, 0, 7.07, 8.47, 10.01, 11.69,
  2.24, 3.08, 3.64, 4.34, 5.18, 6.16, 7.28, 0, 8.61, 10.15, 11.83,
  2.52, 3.5, 4.06, 4.76, 5.6, 6.58, 7.7, 8.96, 0, 10.29, 11.97,
  2.8, 3.92, 4.48, 5.18, 6.02, 7, 8.12, 9.38, 10.78, 0, 12.11,
  3.08, 4.34, 4.9, 5.6, 6.44, 7.42, 8.54, 9.8, 11.2, 12.74, 0
), nrow = 11, ncol = 11, byrow = TRUE)

dyads_P3 <- list(Z1 = Z1_P3, Z2 = Z2_P3)

# ======================================================================================
# Réduction P3 (pour fits rapides)
# ======================================================================================

# Sous-échantillonne partition + dyads en gardant la cohérence des indices.
.subsample_partition_and_dyads <- function(part, dyads, idx) {
  idx <- as.integer(idx)
  stopifnot(all(idx >= 1L), all(idx <= length(part)))
  part2 <- as.integer(part[idx])

  dy2 <- lapply(dyads, function(Z) {
    if (!is.matrix(Z)) stop("dyads doit contenir des matrices.")
    Z[idx, idx, drop = FALSE]
  })

  list(part = part2, dyads = dy2)
}

# Crée une version "petite" de P3 (n=8) : indices 1..8
P3_small <- .subsample_partition_and_dyads(partitions$P3, dyads_P3, idx = 1:8)
partitions$P3s <- P3_small$part
dyads_P3s <- P3_small$dyads

# Sélecteur de matrices dyadiques prédéfinies
.make_dyads_for_partition <- function(part) {
    if (length(part) == length(partitions$P1) && identical(as.integer(part), partitions$P1)) {
        return(dyads_P1)
    }
    if (length(part) == length(partitions$P2) && identical(as.integer(part), partitions$P2)) {
        return(dyads_P2)
    }
    if (length(part) == length(partitions$P3s) && identical(as.integer(part), partitions$P3s)) {
        return(dyads_P3s)
    }
  stop("Aucune matrice dyadique prédéfinie pour cette partition.")
}

# Sanity check: diagonales nulles
stopifnot(all(diag(dyads_P1$Z1) == 0),
          all(diag(dyads_P1$Z2) == 0),
          all(diag(dyads_P2$Z1) == 0),
          all(diag(dyads_P2$Z2) == 0),
          all(diag(dyads_P3$Z1) == 0),
          all(diag(dyads_P3$Z2) == 0))

# ======================================================================================
# Helpers réseau via builder du wrapper
# ======================================================================================

make_network_from_partition_and_dyads <- function(partition_vec, nodes_df, dyads_list) {
  stopifnot(is.atomic(partition_vec), nrow(nodes_df) == length(partition_vec))
  built <- build_bipartite_from_inputs(
    partition = partition_vec,
    nodes     = nodes_df,
    dyads     = dyads_list
  )
  if (is.list(built) && !is.null(built$network)) return(built$network)
  if (inherits(built, "network")) return(built)
  stop("Le builder n’a pas renvoyé un objet 'network'.")
}

# ======================================================================================
# Debug partition / nodes / dyads
# ======================================================================================

print_debug_partition_nodes_dyads <- function(name, part, nodes_df, dyads_list) {
  cat(sprintf("\n[DEBUG] --- Cas %s ---\n", name))
  cat("[DEBUG] partition :", paste(part, collapse = ","), "\n")
  cat("[DEBUG] nodes (head) :\n")
  print(utils::head(nodes_df, 10))

  Z1 <- dyads_list$Z1
  Z2 <- dyads_list$Z2
  n  <- nrow(Z1)
  k  <- min(6L, n)

  cat("[DEBUG] Z1[1:", k, ", 1:", k, "] =\n", sep = "")
  print(round(Z1[seq_len(k), seq_len(k)], 3))

  cat("[DEBUG] Z2[1:", k, ", 1:", k, "] =\n", sep = "")
  print(round(Z2[seq_len(k), seq_len(k)], 3))

  stopifnot(all(diag(Z1) == 0), all(diag(Z2) == 0))
}

# ======================================================================================
# Helpers "erpm -> rhs_expr -> summary" (pour comparer legacy vs matrice)
# ======================================================================================

# (A) Erpm: retourne le rhs_expr traduit + contraintes, sans fit
get_erpm_translated_rhs_expr <- function(partition_vec, nodes_df, dyads_arg, rhs_txt, verbose = FALSE) {
  partition <- partition_vec
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition, nodes = nodes_df), parent = parent.frame())

  call_ergm <- erpm(
    f,
    eval.call = FALSE,
    verbose   = verbose,
    nodes     = nodes_df,
    dyads     = dyads_arg
  )

  call_args <- as.list(call_ergm)[-1L]
  cons <- call_args$constraints
  if (is.null(cons)) cons <- as.formula(~ b1part)

  ergm_form <- call_ergm[[2L]]
  rhs_expr  <- ergm_form[[3L]]

  list(rhs_expr = rhs_expr, constraints = cons, call_ergm = call_ergm)
}

# (B) Summary sur un réseau construit via builder, mais RHS venant de erpm()
summary_from_translated_rhs <- function(nw, rhs_expr, constraints) {
  f2 <- as.formula(bquote(nw ~ .(rhs_expr)))
  environment(f2) <- list2env(list(nw = nw), parent = parent.frame())
  as.numeric(suppressMessages(summary(f2, constraints = constraints)))
}

# ======================================================================================
# Cas test RHS : dyadcov & dyadcov_GW (corrigés selon les InitErgmTerm)
# ======================================================================================

cases_retro <- list(
  # dyadcov : un seul nom dyadique => matrice doit marcher
  list(tag = "dyadcov_Z1_k2",
       rhs = "dyadcov('Z1', clique_size = 2, normalize = FALSE)",
       expects_matrix_ok = TRUE,
       expects_error = FALSE),
  list(tag = "dyadcov_Z2_k3_bygrp",
       rhs = "dyadcov('Z2', clique_size = 3, normalize = 'by_group')",
       expects_matrix_ok = TRUE,
       expects_error = FALSE),

  # dyadcov_GW : l'InitTerm attend dyadcov= et lambda= (pas decay/normalized)
  list(tag = "dyadcov_GW_Z1_lambda05",
       rhs = "dyadcov_GW(dyadcov = 'Z1', lambda = 0.5)",
       expects_matrix_ok = TRUE,
       expects_error = FALSE),

  # Cas volontairement ambigu (2 matrices) : matrice DOIT échouer (nouveau stop explicite),
  # legacy list doit marcher (rétrocompatibilité)
  list(tag = "ambiguous_two_dyads",
       rhs = "dyadcov('Z1', clique_size = 2, normalize = FALSE) + dyadcov('Z2', clique_size = 2, normalize = FALSE)",
       expects_matrix_ok = FALSE,
       expects_error = TRUE)
)

# ======================================================================================
# Phase 1 : Comparer traduction RHS + summary (legacy list vs matrice)
# ======================================================================================

run_phase1_retrocompat_summary <- function() {
  cat("=== PHASE 1 : Rétrocompatibilité erpm() (dyads list vs matrice) via summary ===\n")

  n_ok <- 0L
  n_tot <- 0L

  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes_df_for_partition(part)
    dyads <- .make_dyads_for_partition(part)

    # réseau de référence (builder)
    nw <- make_network_from_partition_and_dyads(part, nodes, dyads)

    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))
    print_debug_partition_nodes_dyads(nm, part, nodes, dyads)

    for (cs in cases_retro) {
      rhs_txt <- cs$rhs
      tag     <- cs$tag

      cat(sprintf("\n[RETRO-SUMMARY] Case=%s | RHS=%s\n", tag, rhs_txt))

      # --- Legacy: dyads list (doit rester OK dans tous les cas utiles) ---
      tr_list <- try(get_erpm_translated_rhs_expr(
        partition_vec = part, nodes_df = nodes, dyads_arg = dyads,
        rhs_txt = rhs_txt, verbose = FALSE
      ), silent = TRUE)

      if (inherits(tr_list, "try-error")) {
        msg <- paste(as.character(tr_list), collapse = "\n")
        cat("  -> LEGACY(list) ERREUR:\n", msg, "\n")
        stop(sprintf("Rétrocompatibilité cassée: erpm() échoue avec dyads=list() sur case=%s partition=%s", tag, nm))
      }

      s_list <- summary_from_translated_rhs(nw, tr_list$rhs_expr, tr_list$constraints)
      if (any(!is.finite(s_list))) {
        stop(sprintf("Legacy(list) summary non-fini sur case=%s partition=%s", tag, nm))
      }
      cat(sprintf("  legacy(list) summary = %s\n", paste(s_list, collapse = ",")))

      # --- New: dyads matrix ---
      # On passe la matrice correspondant au nom attendu quand c'est non-ambigu.
      # Sinon (cas ambiguous_two_dyads), on passe Z1 arbitrairement, et on attend une erreur.
      dyads_mat <- dyads$Z1

      tr_mat <- try(get_erpm_translated_rhs_expr(
        partition_vec = part, nodes_df = nodes, dyads_arg = dyads_mat,
        rhs_txt = rhs_txt, verbose = FALSE
      ), silent = TRUE)

      n_tot <- n_tot + 1L

      if (isTRUE(cs$expects_error)) {
        if (!inherits(tr_mat, "try-error")) {
          cat("  -> MATRIX: attendu une erreur, mais erpm() a accepté.\n")
          stop(sprintf("Cas matrix ambigu aurait dû échouer: case=%s partition=%s", tag, nm))
        } else {
          cat("  -> MATRIX: erreur attendue (OK).\n")
          n_ok <- n_ok + 1L
          next
        }
      }

      if (inherits(tr_mat, "try-error")) {
        msg <- paste(as.character(tr_mat), collapse = "\n")
        cat("  -> MATRIX ERREUR inattendue:\n", msg, "\n")
        stop(sprintf("erpm() n'accepte pas dyads=matrice alors que le RHS est non-ambigu: case=%s partition=%s", tag, nm))
      }

      s_mat <- summary_from_translated_rhs(nw, tr_mat$rhs_expr, tr_mat$constraints)
      if (any(!is.finite(s_mat))) {
        stop(sprintf("Matrix summary non-fini sur case=%s partition=%s", tag, nm))
      }
      cat(sprintf("  matrix summary = %s\n", paste(s_mat, collapse = ",")))

      # --- Critère rétrocompatibilité: summary identique ---
      # On exige égalité numérique stricte ici car summary() est déterministe.
      if (length(s_list) != length(s_mat) || any(abs(s_list - s_mat) > 0)) {
        stop(sprintf("Mismatch legacy(list) vs matrix sur case=%s partition=%s", tag, nm))
      }

      # Bonus: vérifier que le RHS traduit est identique
      d_list <- paste(deparse(tr_list$rhs_expr), collapse = " ")
      d_mat  <- paste(deparse(tr_mat$rhs_expr),  collapse = " ")
      if (!identical(d_list, d_mat)) {
        cat("  [WARN] RHS traduit diffère entre list et matrix (peut être OK si noms/coefs changent mais stats identiques).\n")
        cat("         list :", d_list, "\n")
        cat("         mat  :", d_mat,  "\n")
      }

      cat("  -> OK (legacy == matrix)\n")
      n_ok <- n_ok + 1L
    }
  }

  cat(sprintf("\n=== Bilan Phase 1 : %d / %d checks OK ===\n", n_ok, n_tot))
  invisible(NULL)
}

# ======================================================================================
# Phase 2 : Fits erpm() (legacy list vs matrice) sur un sous-ensemble stable
# ======================================================================================

run_one_erpm_fit <- function(partition_vec, nodes_df, dyads_arg, rhs_txt, fit_name) {
  set.seed(42)
  partition <- partition_vec
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition, nodes = nodes_df), parent = parent.frame())

  cat(sprintf("[ERPM-FIT %-26s] n=%-3d RHS=%s\n", fit_name, length(partition_vec), rhs_txt))

  fit <- try(
    erpm(
      f,
      verbose = FALSE,
      nodes   = nodes_df,
      dyads   = dyads_arg
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    msg <- paste(as.character(fit), collapse = "\n")
    cat("  -> ERREUR fit:", msg, "\n")
    return(list(ok = FALSE, fit = NULL, coef = NA))
  }

  ok_class <- inherits(fit, "ergm")
  cf <- try(stats::coef(fit), silent = TRUE)
  ok_coef <- !inherits(cf, "try-error") && all(is.finite(cf))

  cat(sprintf("  -> class(ergm)? %s | coef finies? %s | coef: %s\n",
              if (ok_class) "OK" else "KO",
              if (ok_coef) "OK" else "KO",
              if (ok_coef) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

  list(ok = ok_class && ok_coef, fit = fit, coef = cf)
}

run_phase2_retrocompat_fits <- function() {
  cat("\n=== PHASE 2 : Fits erpm() (legacy list vs matrice) ===\n")

  # Fits non-ambigus (un seul nom dyadique dans le RHS)
  rhs_fits <- list(
    F1 = "dyadcov('Z1', clique_size = 2, normalize = FALSE) + cliques",
    F2 = "dyadcov('Z2', clique_size = 3, normalize = 'by_group') + cliques",
    # dyadcov_GW: dyadcov= + lambda= (pas decay/normalized)
    F3 = "dyadcov_GW(dyadcov = 'Z1', lambda = 0.5) + cliques"
  )

  combos <- list(
    list(p = "P1", f = "F1"),
    list(p = "P2", f = "F1"),
    list(p = "P3s", f = "F2"),
    list(p = "P3s", f = "F3")
  )

  for (cb in combos) {
    part  <- partitions[[cb$p]]
    nodes <- .make_nodes_df_for_partition(part)
    dyads <- .make_dyads_for_partition(part)

    rhs_txt <- rhs_fits[[cb$f]]
    tag <- paste0(cb$p, "_", cb$f)

    cat(sprintf("\n--- FIT COMBO %s ---\n", tag))
    print_debug_partition_nodes_dyads(paste0("FIT_", tag), part, nodes, dyads)

    # Legacy list
    fit_list <- run_one_erpm_fit(part, nodes, dyads, rhs_txt, paste0(tag, "_legacy_list"))
    if (!isTRUE(fit_list$ok)) stop(sprintf("Fit legacy(list) a échoué: %s", tag))

    # New matrix: choisir matrice correspondant au RHS (heuristique simple)
    # - dyadcov('Z2', ...) => Z2
    # - dyadcov_GW(... 'Z1' ...) => Z1
    use_mat <- if (grepl("'Z2'", rhs_txt, fixed = TRUE)) dyads$Z2 else dyads$Z1

    fit_mat <- run_one_erpm_fit(part, nodes, use_mat, rhs_txt, paste0(tag, "_new_matrix"))
    if (!isTRUE(fit_mat$ok)) stop(sprintf("Fit matrix a échoué: %s", tag))

    # Comparaison faible: même longueur de coef + valeurs finies.
    # (Ne force pas égalité numérique stricte: MCMC, seeds, etc.)
    cf1 <- as.numeric(fit_list$coef)
    cf2 <- as.numeric(fit_mat$coef)

    if (length(cf1) != length(cf2) || any(!is.finite(cf1)) || any(!is.finite(cf2))) {
      stop(sprintf("Coefficients invalides ou dimensions différentes: %s", tag))
    }

    cat(sprintf("[RETRO-FIT] %s OK (legacy(list) + matrix) | dim=%d\n", tag, length(cf1)))
  }

  cat("\n=== Phase 2 OK: fits legacy(list) et matrix passent sur un sous-ensemble stable ===\n")
  invisible(NULL)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: rétrocompat dyads=list vs dyads=matrix | dyadcov + dyadcov_GW ===\n")
run_phase1_retrocompat_summary()
run_phase2_retrocompat_fits()

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)
cat("\nTous les tests de rétrocompatibilité dyadcov/dyadcov_GW ont passé.\n")