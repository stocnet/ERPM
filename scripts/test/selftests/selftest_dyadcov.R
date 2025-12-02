# ======================================================================================
# Fichier : scripts/test/selftests/selftest_dyadcov.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `dyadcov`
# Exécution: Rscript scripts/test/selftests/selftest_dyadcov.R
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
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_dyadcov.log")
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

# Sélecteur de matrices dyadiques prédéfinies
.make_dyads_for_partition <- function(part) {
  if (length(part) == length(partitions$P1) && identical(as.integer(part), partitions$P1)) {
    return(dyads_P1)
  }
  if (length(part) == length(partitions$P2) && identical(as.integer(part), partitions$P2)) {
    return(dyads_P2)
  }
  if (length(part) == length(partitions$P3) && identical(as.integer(part), partitions$P3)) {
    return(dyads_P3)
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

make_formula_for_network_summary <- function(nw, rhs_txt) {
  f <- as.formula(paste0("nw ~ ", rhs_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  f
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
# Référence R directe de la stat dyadcov
# ======================================================================================

dyadcov_reference_from_partition <- function(partition_vec, Z, k, normalized = FALSE) {
  n <- length(partition_vec)
  groups <- sort(unique(partition_vec))
  total <- 0
  for (g in groups) {
    idx_g <- which(partition_vec == g)
    n_g   <- length(idx_g)
    if (n_g < k) next

    S_g <- 0
    comb_mat <- utils::combn(idx_g, k)
    if (!is.matrix(comb_mat)) comb_mat <- matrix(comb_mat, nrow = k)

    for (col in seq_len(ncol(comb_mat))) {
      clique <- comb_mat[, col]
      prod_C <- 1
      if (k >= 2L) {
        for (p in 1:(k - 1L)) {
          i <- clique[p]
          for (q in (p + 1L):k) {
            j <- clique[q]
            prod_C <- prod_C * (Z[i, j] + Z[j, i])
          }
        }
      }
      S_g <- S_g + prod_C
    }

    if (!normalized) {
      total <- total + S_g
    } else {
      total <- total + S_g / n_g
    }
  }
  as.numeric(total)
}

run_phase0_analytic_checks_dyadcov <- function() {
  cat("=== PHASE 0 : Vérifications analytiques directes [dyadcov] ===\n")

  part  <- partitions$P1
  nodes <- .make_nodes_df_for_partition(part)
  dyads <- .make_dyads_for_partition(part)
  Z1    <- dyads$Z1

  nw <- make_network_from_partition_and_dyads(part, nodes, dyads)

  tests <- list(
    list(Zname = "Z1", Z = Z1, k = 2L, normalized = FALSE),
    list(Zname = "Z1", Z = Z1, k = 2L, normalized = TRUE),
    list(Zname = "Z1", Z = Z1, k = 3L, normalized = FALSE),
    list(Zname = "Z1", Z = Z1, k = 3L, normalized = TRUE)
  )

  for (ts in tests) {
    rhs <- sprintf("dyadcov('%s', clique_size = %d, normalized = %s)",
                   ts$Zname, ts$k, if (ts$normalized) "TRUE" else "FALSE")
    f   <- make_formula_for_network_summary(nw, rhs)
    val_summary <- as.numeric(suppressMessages(summary(f, constraints = ~ b1part)))
    val_ref     <- dyadcov_reference_from_partition(part, ts$Z, ts$k, ts$normalized)

    cat(sprintf("[ANALYTIC-CHECK] k=%d normalized=%s  summary=%g  ref=%g  diff=%g\n",
                ts$k, ts$normalized, val_summary, val_ref, val_summary - val_ref))

    if (!is.finite(val_summary) || !is.finite(val_ref) ||
        abs(val_summary - val_ref) > 1e-8) {
      stop(sprintf("Mismatch analytique dyadcov: k=%d normalized=%s diff=%g",
                   ts$k, ts$normalized, val_summary - val_ref))
    }
  }

  cat("=== Phase 0 OK: définition dyadcov (normalisation 1/n_g) validée sur cas simples ===\n\n")
  invisible(NULL)
}

# ======================================================================================
# Fonctions SUMMARY / ERPM
# ======================================================================================

run_one_network_summary_case_for_dyadcov <- function(partition_vec, nodes_df, dyads_list, rhs_txt) {
  nw <- make_network_from_partition_and_dyads(partition_vec, nodes_df, dyads_list)
  f  <- make_formula_for_network_summary(nw, rhs_txt)
  as.numeric(suppressMessages(summary(f, constraints = ~ b1part)))
}

run_one_erpm_translated_summary_case_for_dyadcov <- function(partition_vec, nodes_df, dyads_list, rhs_txt) {
  partition <- partition_vec
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition, nodes = nodes_df), parent = parent.frame())

  call_ergm <- erpm(
    f,
    eval_call = FALSE,
    verbose   = TRUE,
    nodes     = nodes_df,
    dyads     = dyads_list
  )

  ergm_form <- call_ergm[[2L]]
  rhs_expr  <- ergm_form[[3L]]

  nw2 <- make_network_from_partition_and_dyads(partition_vec, nodes_df, dyads_list)
  f2  <- as.formula(bquote(nw2 ~ .(rhs_expr)))
  environment(f2) <- list2env(list(nw2 = nw2), parent = parent.frame())

  call_args <- as.list(call_ergm)[-1L]
  cons <- call_args$constraints
  if (is.null(cons)) cons <- as.formula(~ b1part)

  as.numeric(suppressMessages(summary(f2, constraints = cons)))
}

check_summary_equivalence_network_vs_erpm_dyadcov <- function(partition_vec, nodes_df, dyads_list,
                                                              rhs_vec, tol = 0) {
  ok_all <- TRUE
  for (rhs in rhs_vec) {
    s_net  <- run_one_network_summary_case_for_dyadcov(partition_vec, nodes_df, dyads_list, rhs)
    s_erpm <- run_one_erpm_translated_summary_case_for_dyadcov(partition_vec, nodes_df, dyads_list, rhs)
    cat(sprintf("[SUMMARY-CHECK] n=%-3d RHS=%-50s net=%s  erpm=%s\n",
                length(partition_vec), rhs,
                paste(s_net,  collapse=","), paste(s_erpm, collapse=",")))
    if (any(!is.finite(s_net)) || any(!is.finite(s_erpm)) ||
        length(s_net) != length(s_erpm) ||
        !all(abs(s_net - s_erpm) <= tol)) {
      ok_all <- FALSE
      cat("  -> MISMATCH détecté.\n")
    }
  }
  ok_all
}

cases_summary <- c(
  "dyadcov('Z1', clique_size = 2, normalized = FALSE)",
  "dyadcov('Z1', clique_size = 2, normalized = TRUE)",
  "dyadcov('Z1', clique_size = 3, normalized = FALSE)",
  "dyadcov('Z1', clique_size = 3, normalized = TRUE)",
  "dyadcov('Z2', clique_size = 2, normalized = FALSE)",
  "dyadcov('Z2', clique_size = 2, normalized = TRUE)"
)

run_phase1_summary_equivalence_checks_dyadcov <- function() {
  cat("=== PHASE 1 : Summary(nw via builder) vs Summary(ERPM-traduit) [dyadcov] ===\n")
  total <- 0L; ok <- 0L
  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes_df_for_partition(part)
    dyads <- .make_dyads_for_partition(part)

    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))

    print_debug_partition_nodes_dyads(nm, part, nodes, dyads)

    res <- check_summary_equivalence_network_vs_erpm_dyadcov(
      partition_vec = part,
      nodes_df      = nodes,
      dyads_list    = dyads,
      rhs_vec       = cases_summary,
      tol           = 0
    )
    total <- total + length(cases_summary)
    ok    <- ok + as.integer(res) * length(cases_summary)
  }
  cat(sprintf("\n=== Bilan Phase 1 (dyadcov) : %d / %d checks OK ===\n", ok, total))
  if (ok < total) stop(sprintf("Summary mismatch sur %d cas.", total - ok))
  invisible(NULL)
}

# ======================================================================================
# Phase 2: Fits via ERPM
# ======================================================================================

run_one_erpm_fit_with_return_dyadcov <- function(partition_vec, nodes_df, dyads_list,
                                                 rhs_txt, fit_name) {
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT %-20s] SKIP (erpm() indisponible)\n", fit_name))
    return(list(ok = NA, error = TRUE, coef = NA, fit = NULL, aic = NA, bic = NA))
  }
  set.seed(42)
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition_vec, nodes = nodes_df), parent = parent.frame())

  cat(sprintf("[ERPM-FIT %-20s] n=%-3d RHS=%s\n",
              fit_name, length(partition_vec), rhs_txt))

  print_debug_partition_nodes_dyads(paste0("FIT_", fit_name), partition_vec, nodes_df, dyads_list)

  fit <- try(
    erpm(
      f,
      verbose = FALSE,
      nodes   = nodes_df,
      dyads   = dyads_list
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    msg <- paste(as.character(fit), collapse = "\n")
    cat("  -> ERREUR fit:", msg, "\n")
    return(list(ok = FALSE, error = TRUE, coef = NA, fit = NULL, aic = NA, bic = NA))
  }

  ok_class <- inherits(fit, "ergm")
  cf <- try(stats::coef(fit), silent = TRUE)
  ok_coef <- !inherits(cf, "try-error") && all(is.finite(cf))

  aic_val <- NA_real_
  bic_val <- NA_real_
  ll <- try(logLik(fit, add = TRUE), silent = TRUE)
  if (!inherits(ll, "try-error")) {
    aic_try <- try(AIC(fit), silent = TRUE)
    bic_try <- try(BIC(fit), silent = TRUE)
    if (!inherits(aic_try, "try-error")) aic_val <- as.numeric(aic_try)
    if (!inherits(bic_try, "try-error")) bic_val <- as.numeric(bic_try)
  } else {
    cat("  -> logLik(fit, add=TRUE) a échoué, AIC/BIC indisponibles pour ce fit.\n")
  }

  cat(sprintf("  -> class(ergm)? %s | coef finies? %s | coef: %s | AIC=%s | BIC=%s\n",
              if (ok_class) "OK" else "KO",
              if (ok_coef) "OK" else "KO",
              if (ok_coef) paste(format(as.numeric(cf)), collapse=", ") else "NA",
              if (is.finite(aic_val)) format(aic_val, digits = 6) else "NA",
              if (is.finite(bic_val)) format(bic_val, digits = 6) else "NA"))

  list(
    ok    = ok_class && ok_coef,
    error = FALSE,
    coef  = if (ok_coef) cf else NA,
    fit   = fit,
    aic   = aic_val,
    bic   = bic_val
  )
}

run_phase2_erpm_fits_and_print_summaries_dyadcov <- function() {
  cat("\n=== PHASE 2 : Fits erpm() [dyadcov] ===\n")

  rhs_list <- list(
    R1 = "dyadcov('Z2', clique_size = 3, normalized = FALSE) + cliques",
    R2 = "dyadcov('Z1', clique_size = 2, normalized = FALSE) + cliques",
    R3 = "dyadcov('Z1', clique_size = 2, normalized = TRUE) + cliques"
  )

  parts <- list(
    list(name = "P1", part = partitions$P1),
    list(name = "P2", part = partitions$P2),
    list(name = "P3", part = partitions$P3)
  )

  fit_results <- list()
  # On garde les deux fits stables sur P1/P2 et deux fits sur P3
  combos <- list(
    list(p = "P1", r = "R1"),
    list(p = "P2", r = "R1"),
    list(p = "P3", r = "R2"),
    list(p = "P3", r = "R3")
  )

  for (cb in combos) {
    px <- parts[[which(vapply(parts, function(x) x$name == cb$p, logical(1)))]]
    nodes <- .make_nodes_df_for_partition(px$part)
    dyads <- .make_dyads_for_partition(px$part)
    nmfit <- paste0(px$name, "_", cb$r)
    fit_results[[nmfit]] <- run_one_erpm_fit_with_return_dyadcov(
      partition_vec = px$part,
      nodes_df      = nodes,
      dyads_list    = dyads,
      rhs_txt       = rhs_list[[cb$r]],
      fit_name      = nmfit
    )
  }

  ok_raw  <- vapply(fit_results, function(x) x$ok,    logical(1))
  n_ok    <- sum(ok_raw, na.rm = TRUE)
  n_tot   <- sum(!is.na(ok_raw))

  cat(sprintf("\n=== Bilan fits erpm() [dyadcov] : %d / %d OK ===\n", n_ok, n_tot))

  cat("\n=== Tableau AIC/BIC pour les fits ERPM (dyadcov) ===\n")
  tab <- data.frame(
    fit   = character(0),
    ok    = logical(0),
    AIC   = numeric(0),
    BIC   = numeric(0),
    stringsAsFactors = FALSE
  )
  for (nm in names(fit_results)) {
    fr <- fit_results[[nm]]
    tab <- rbind(tab, data.frame(
      fit = nm,
      ok  = fr$ok,
      AIC = fr$aic,
      BIC = fr$bic,
      stringsAsFactors = FALSE
    ))
  }
  print(tab)

  cat("\n=== Résumés détaillés des fits ERPM réussis (dyadcov) ===\n")
  for (nm in names(fit_results)) {
    fit_obj <- fit_results[[nm]]
    if (isTRUE(fit_obj$ok) && inherits(fit_obj$coef, "numeric") && !is.null(fit_obj$fit)) {
      cat(sprintf("\n--- Résumé fit %s ---\n", nm))
      print(summary(fit_obj$fit))
    }
  }

  # Ici tu peux réajuster le seuil selon ce que tu veux imposer
  required_ok <- n_tot

  if (n_ok < required_ok) {
    stop(sprintf("Echec fits: seulement %d modèles convergents (seuil=%d).",
                 n_ok, required_ok))
  } else {
    cat(sprintf("\nSeuil de convergence atteint : %d modèles OK (seuil=%d).\n",
                n_ok, required_ok))
  }

  invisible(fit_results)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: dyadcov ===\n")
run_phase0_analytic_checks_dyadcov()
run_phase1_summary_equivalence_checks_dyadcov()
fit_results <- run_phase2_erpm_fits_and_print_summaries_dyadcov()

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)
cat("\nTous les tests dyadcov ont passé.\n")