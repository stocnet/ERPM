# ==============================================================================
# File    : scripts/test/selftests/selftest_cov_diff.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_diff` (multi-toggle)
# Exécution: Rscript scripts/test/selftests/selftest_cov_diff.R
#
# But du fichier
#   - PHASE 1 (SUMMARY) : valider la statistique via summary(nw ~ cov_diff(...))
#                         vs calcul analytique explicite (combn).
#   - PHASE 2 (EQUIV)   : valider summary() réseau explicite vs réseau ERPM-traduit.
#   - PHASE 3 (FIT)     : sanity check erpm() (coefs finis) avec un terme structurel.
#   - PHASE 4 (MCMC)    : probe multi-toggle (le point clé après migration D_CHANGESTAT_FN).
#
# Important
#   - Le probe multi-toggle (PHASE 4) n'a pas pour but de “bien simuler” :
#     on veut juste déclencher le chemin MCMC et voir passer (si activé)
#     les traces debug C (ex: "[cov_diff] MULTI-TOGGLE ntoggles=...").
#   - Pour voir ces traces:
#       * côté C: mettre DEBUG_COV_DIFF = 1 dans changestat_cov_diff.c
#       * recompiler le package (load_all(recompile=TRUE) ou install)
# ==============================================================================

# --------------------------------------------------------------------------------------
# Préambule
# --------------------------------------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

options(ergm.loglik.warn_dyads = FALSE)

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
  if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
}

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Le fichier DESCRIPTION n'existe pas ou devtools n'est pas installé.")
}

# --------------------------------------------------------------------------------------
# Logging local
# --------------------------------------------------------------------------------------
root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_cov_diff.log")
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

# --------------------------------------------------------------------------------------
# Chargement utilitaires ERPM / wrapper
# --------------------------------------------------------------------------------------
if (!exists("partition_to_bipartite_network", mode = "function")) {
  if (file.exists("R/functions_erpm_bip_network.R")) {
    source("R/functions_erpm_bip_network.R", local = FALSE)
  }
}
if (!exists("erpm", mode = "function") || !exists("build_bipartite_from_inputs", mode = "function")) {
  if (file.exists("R/erpm_wrapper.R")) {
    source("R/erpm_wrapper.R", local = FALSE)
  } else {
    cat("[WARN] erpm()/build_bipartite_from_inputs indisponibles. Certaines étapes seront sautées.\n")
  }
}

# ======================================================================================
# Réglages de run (activation/désactivation des phases)
# ======================================================================================
RUN <- list(
  phase1_summary_expected = FALSE,
  phase2_equiv_summary    = FALSE,
  phase3_erpm_fit         = TRUE,
  phase4_mcmc_probe       = FALSE,

  quiet_phase1 = FALSE,
  quiet_phase2 = FALSE,
  quiet_phase3 = FALSE
)

# ======================================================================================
# Données de test pour cov_diff (alignées sur cub_test cov_diff)
# ======================================================================================

# PartA : 3 groupes (4,5,6) = 15 sommets
partA <- c(
  rep(1, 4),
  rep(2, 5),
  rep(3, 6)
)

nodesA <- data.frame(
  label = paste0("A", seq_along(partA)),
  score = c(
    10, 11,  9, 10,             # g1 (variation faible)
     5,  7, 20,  9,  6,         # g2 (variation forte)
     0,  2,  4,  6,  8, 10      # g3 (croissant régulier)
  ),
  stringsAsFactors = FALSE
)

# PartB : 4 groupes (4,4,5,5) = 18 sommets
partB <- c(
  rep(1, 4),
  rep(2, 4),
  rep(3, 5),
  rep(4, 5)
)

nodesB <- data.frame(
  label = paste0("B", seq_along(partB)),
  score = c(
    10, 10, 11,  9,        # g1 (homogène)
     0, 40, 80,  5,        # g2 (très étalé)
     3,  4,  5,  6,  7,    # g3 (rampe régulière)
     0, 20, 40, 60, 80     # g4 (étalement fort)
  ),
  stringsAsFactors = FALSE
)

# PartC : 3 groupes (4,5,4) = 13 sommets
partC <- c(
  rep(1, 4),
  rep(2, 5),
  rep(3, 4)
)

nodesC <- data.frame(
  label = paste0("C", seq_along(partC)),
  score = c(
    5, 5, 5, 5,        # g1 -> dispersion nulle
    1, 2, 4, 3, 5,     # g2 -> dispersion modérée
    0, 0,100, 0        # g3 -> gros outlier
  ),
  stringsAsFactors = FALSE
)

# Panel générique de partitions pour les tests d'équivalence summary
partitions <- list(
  A = partA,
  B = partB,
  C = partC
)

.make_nodes_numeric <- function(part) {
  n <- length(part)
  set.seed(200 + n)
  data.frame(
    label = paste0("N", seq_len(n)),
    score = sample(0:20, n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# ======================================================================================
# Helpers réseau biparti et summaries
# ======================================================================================

.erpm_build_bipartite_nw <- function(part, nodes) {
  stopifnot(length(part) == nrow(nodes))
  attrs <- as.list(nodes[ , setdiff(names(nodes),"label"), drop=FALSE])

  if (exists("build_bipartite_from_inputs", mode = "function")) {
    builder <- get("build_bipartite_from_inputs")
    out <- try(builder(partition = part, nodes = nodes), silent = TRUE)
    if (inherits(out, "try-error") || is.null(out)) {
      out <- try(builder(partition = part, labels = nodes$label, attributes = attrs), silent = TRUE)
    }
    if (!inherits(out, "try-error") && !is.null(out)) {
      if (inherits(out, "network")) return(out)
      if (is.list(out)) {
        for (nm in c("network","nw","g","graph","bip","net")) {
          if (!is.null(out[[nm]]) && inherits(out[[nm]], "network")) return(out[[nm]])
        }
      }
    }
  }

  if (exists("partition_to_bipartite_network", mode = "function")) {
    return(partition_to_bipartite_network(labels = nodes$label, partition = part, attributes = attrs))
  }

  stop("Aucun constructeur biparti valide n'a produit un objet 'network'.")
}

.formula_nw <- function(nw, rhs_txt) {
  f <- as.formula(paste0("nw ~ ", rhs_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  f
}

summary_on_bipartite_network <- function(part, nodes, rhs_txt) {
  nw <- .erpm_build_bipartite_nw(part, nodes)
  f  <- .formula_nw(nw, rhs_txt)
  as.numeric(suppressMessages(summary(f)))
}

summary_on_erpm_translation <- function(part, nodes, rhs_txt) {
  if (!exists("erpm", mode = "function")) return(NA_real_)

  partition <- part
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition, nodes = nodes), parent = parent.frame())

  call_ergm <- erpm(f, eval.call = FALSE, verbose = FALSE, nodes = nodes)

  ergm_form <- call_ergm[[2L]]
  rhs_expr  <- ergm_form[[3L]]

  nw <- .erpm_build_bipartite_nw(part, nodes)
  f2 <- as.formula(bquote(nw ~ .(rhs_expr)))
  environment(f2) <- list2env(list(nw = nw), parent = parent.frame())

  call_args <- as.list(call_ergm)[-1L]
  cons <- call_args$constraints
  if (is.null(cons)) cons <- as.formula(~ b1part)

  as.numeric(suppressMessages(summary(f2, constraints = cons)))
}

# ======================================================================================
# Attentes analytiques pour cov_diff
# ======================================================================================

# Stat non normalisée :
#   T_k = sum_g sum_{S in C_k(g)} (max(x_S) - min(x_S))
# Stat normalisée par groupe :
#   T_k^{by_group} = sum_g [ (1/choose(n_g,k)) * sum_{S in C_k(g)} (...) ]
# Stat normalisée globale :
#   T_k^{global} = sum_g [ (1/n_g) * sum_{S in C_k(g)} (...) ]  (si n_g>=k)
expected_cov_diff <- function(part, x, clique_size = 2L, normalized = FALSE) {
  stopifnot(length(part) == length(x))
  k <- as.integer(clique_size)
  if (k < 2L) stop("expected_cov_diff: clique_size doit être >= 2.")

  norm_mode <- 0L
  if (is.logical(normalized)) {
    norm_mode <- if (isTRUE(normalized)) 1L else 0L
  } else if (is.character(normalized) && length(normalized) == 1L) {
    normalized <- tolower(normalized)
    norm_mode <- if (normalized == "by_group") 1L else if (normalized == "global") 2L else 0L
  } else if (is.numeric(normalized) && length(normalized) == 1L) {
    norm_mode <- as.integer(normalized)
  }

  gid <- as.integer(part)
  split_idx <- split(seq_along(gid), gid)
  tot <- 0

  for (ix in split_idx) {
    ng <- length(ix)
    if (ng < k) next

    v <- x[ix]
    idx_mat <- utils::combn(ng, k)
    diffs <- apply(idx_mat, 2L, function(idr) {
      vals <- v[idr]
      max(vals) - min(vals)
    })
    Sg <- sum(diffs)

    if (norm_mode == 0L) {
      tot <- tot + Sg
    } else if (norm_mode == 1L) {
      tot <- tot + Sg / choose(ng, k)
    } else if (norm_mode == 2L) {
      tot <- tot + Sg / ng
    } else {
      stop("expected_cov_diff: norm_mode invalide")
    }
  }

  tot
}

# ======================================================================================
# Phase 1: Summary — attentes explicites (analytique vs summary)
# ======================================================================================

run_phase1_summary_expected_cov_diff <- function(quiet = FALSE) {
  cat("=== PHASE 1 : Summary cov_diff avec attentes explicites ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  # ---------- Cas A ----------
  exp_A_k2_raw  <- expected_cov_diff(partA, nodesA$score, clique_size = 2L, normalized = FALSE)
  s_A_k2_raw    <- summary_on_bipartite_network(partA, nodesA, "cov_diff('score', clique_size = 2)")
  if (!isTRUE(quiet)) cat(sprintf("[A] cov_diff('score',k=2)              obtenu=%g  attendu=%g\n", s_A_k2_raw, exp_A_k2_raw))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_k2_raw), as.numeric(exp_A_k2_raw))))

  exp_A_k2_norm <- expected_cov_diff(partA, nodesA$score, clique_size = 2L, normalized = TRUE)
  s_A_k2_norm   <- summary_on_bipartite_network(partA, nodesA, "cov_diff('score', clique_size = 2, normalized = TRUE)")
  if (!isTRUE(quiet)) cat(sprintf("[A] cov_diff('score',k=2,norm=TRUE)    obtenu=%g  attendu=%g\n", s_A_k2_norm, exp_A_k2_norm))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_k2_norm), as.numeric(exp_A_k2_norm))))

  exp_A_k2_glob <- expected_cov_diff(partA, nodesA$score, clique_size = 2L, normalized = "global")
  s_A_k2_glob   <- summary_on_bipartite_network(partA, nodesA, "cov_diff('score', clique_size = 2, normalize = 'global')")
  if (!isTRUE(quiet)) cat(sprintf("[A] cov_diff('score',k=2,glob)         obtenu=%g  attendu=%g\n", s_A_k2_glob, exp_A_k2_glob))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_k2_glob), as.numeric(exp_A_k2_glob))))

  exp_A_k3_raw  <- expected_cov_diff(partA, nodesA$score, clique_size = 3L, normalized = FALSE)
  s_A_k3_raw    <- summary_on_bipartite_network(partA, nodesA, "cov_diff('score', clique_size = 3)")
  if (!isTRUE(quiet)) cat(sprintf("[A] cov_diff('score',k=3)              obtenu=%g  attendu=%g\n", s_A_k3_raw, exp_A_k3_raw))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_k3_raw), as.numeric(exp_A_k3_raw))))

  # ---------- Cas B ----------
  exp_B_k2_raw  <- expected_cov_diff(partB, nodesB$score, clique_size = 2L, normalized = FALSE)
  s_B_k2_raw    <- summary_on_bipartite_network(partB, nodesB, "cov_diff('score', clique_size = 2)")
  if (!isTRUE(quiet)) cat(sprintf("[B] cov_diff('score',k=2)              obtenu=%g  attendu=%g\n", s_B_k2_raw, exp_B_k2_raw))
  stopifnot(isTRUE(all.equal(as.numeric(s_B_k2_raw), as.numeric(exp_B_k2_raw))))

  exp_B_k2_norm <- expected_cov_diff(partB, nodesB$score, clique_size = 2L, normalized = TRUE)
  s_B_k2_norm   <- summary_on_bipartite_network(partB, nodesB, "cov_diff('score', clique_size = 2, normalized = TRUE)")
  if (!isTRUE(quiet)) cat(sprintf("[B] cov_diff('score',k=2,norm=TRUE)    obtenu=%g  attendu=%g\n", s_B_k2_norm, exp_B_k2_norm))
  stopifnot(isTRUE(all.equal(as.numeric(s_B_k2_norm), as.numeric(exp_B_k2_norm))))

  exp_B_k2_glob <- expected_cov_diff(partB, nodesB$score, clique_size = 2L, normalized = "global")
  s_B_k2_glob   <- summary_on_bipartite_network(partB, nodesB, "cov_diff('score', clique_size = 2, normalize = 'global')")
  if (!isTRUE(quiet)) cat(sprintf("[B] cov_diff('score',k=2,glob)         obtenu=%g  attendu=%g\n", s_B_k2_glob, exp_B_k2_glob))
  stopifnot(isTRUE(all.equal(as.numeric(s_B_k2_glob), as.numeric(exp_B_k2_glob))))

  exp_B_k3_raw  <- expected_cov_diff(partB, nodesB$score, clique_size = 3L, normalized = FALSE)
  s_B_k3_raw    <- summary_on_bipartite_network(partB, nodesB, "cov_diff('score', clique_size = 3)")
  if (!isTRUE(quiet)) cat(sprintf("[B] cov_diff('score',k=3)              obtenu=%g  attendu=%g\n", s_B_k3_raw, exp_B_k3_raw))
  stopifnot(isTRUE(all.equal(as.numeric(s_B_k3_raw), as.numeric(exp_B_k3_raw))))

  # ---------- Cas C ----------
  exp_C_k2_raw  <- expected_cov_diff(partC, nodesC$score, clique_size = 2L, normalized = FALSE)
  s_C_k2_raw    <- summary_on_bipartite_network(partC, nodesC, "cov_diff('score', clique_size = 2)")
  if (!isTRUE(quiet)) cat(sprintf("[C] cov_diff('score',k=2)              obtenu=%g  attendu=%g\n", s_C_k2_raw, exp_C_k2_raw))
  stopifnot(isTRUE(all.equal(as.numeric(s_C_k2_raw), as.numeric(exp_C_k2_raw))))

  exp_C_k2_norm <- expected_cov_diff(partC, nodesC$score, clique_size = 2L, normalized = TRUE)
  s_C_k2_norm   <- summary_on_bipartite_network(partC, nodesC, "cov_diff('score', clique_size = 2, normalized = TRUE)")
  if (!isTRUE(quiet)) cat(sprintf("[C] cov_diff('score',k=2,norm=TRUE)    obtenu=%g  attendu=%g\n", s_C_k2_norm, exp_C_k2_norm))
  stopifnot(isTRUE(all.equal(as.numeric(s_C_k2_norm), as.numeric(exp_C_k2_norm))))

  exp_C_k2_glob <- expected_cov_diff(partC, nodesC$score, clique_size = 2L, normalized = "global")
  s_C_k2_glob   <- summary_on_bipartite_network(partC, nodesC, "cov_diff('score', clique_size = 2, normalize = 'global')")
  if (!isTRUE(quiet)) cat(sprintf("[C] cov_diff('score',k=2,glob)         obtenu=%g  attendu=%g\n", s_C_k2_glob, exp_C_k2_glob))
  stopifnot(isTRUE(all.equal(as.numeric(s_C_k2_glob), as.numeric(exp_C_k2_glob))))

  exp_C_k3_raw  <- expected_cov_diff(partC, nodesC$score, clique_size = 3L, normalized = FALSE)
  s_C_k3_raw    <- summary_on_bipartite_network(partC, nodesC, "cov_diff('score', clique_size = 3)")
  if (!isTRUE(quiet)) cat(sprintf("[C] cov_diff('score',k=3)              obtenu=%g  attendu=%g\n", s_C_k3_raw, exp_C_k3_raw))
  stopifnot(isTRUE(all.equal(as.numeric(s_C_k3_raw), as.numeric(exp_C_k3_raw))))

  cat("\n=== Phase 1 OK ===\n")
  invisible(TRUE)
}

# ======================================================================================
# Phase 2: Summary — équivalence réseau explicite vs ERPM-traduit
# ======================================================================================

cases_equiv_diff <- c(
  "cov_diff('score', clique_size = 2)",
  "cov_diff('score', clique_size = 2, normalized = TRUE)",
  "cov_diff('score', clique_size = 2, normalize = 'global')",
  "cov_diff('score', clique_size = 3)"
)

check_summary_equivalence_cov_diff <- function(part, nodes, rhs_vec) {
  for (rhs in rhs_vec) {
    s_net  <- summary_on_bipartite_network(part, nodes, rhs)
    s_erpm <- summary_on_erpm_translation(part, nodes, rhs)
    cat(sprintf("[EQUIV-DIFF] n=%-3d RHS=%-55s net=%s | erpm=%s\n",
                length(part), rhs, paste(s_net, collapse=","), paste(s_erpm, collapse=",")))
    if (length(s_net) != length(s_erpm)) stop("Longueur de statistique différente.")
    if (!all(is.finite(s_net)) || !all(is.finite(s_erpm))) stop("Stat non finie.")
    if (!isTRUE(all.equal(as.numeric(s_net), as.numeric(s_erpm)))) {
      stop(sprintf("Mismatch summary net vs ERPM pour RHS=%s", rhs))
    }
  }
  TRUE
}

run_phase2_summary_equiv_cov_diff <- function(quiet = FALSE) {
  cat("\n=== PHASE 2 : Summary(nw) vs Summary(ERPM-traduit) pour cov_diff ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  # Cas aléatoires (scores simulés)
  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes_numeric(part)
    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))
    ok <- check_summary_equivalence_cov_diff(part, nodes, cases_equiv_diff)
    if (!ok) stop("Equivalence summary échouée.")
  }

  # Cas déterministes A,B,C sur 'score'
  cat("\n--- Partition A (déterministe) ---\n")
  stopifnot(check_summary_equivalence_cov_diff(partA, nodesA, cases_equiv_diff))

  cat("\n--- Partition B (déterministe) ---\n")
  stopifnot(check_summary_equivalence_cov_diff(partB, nodesB, cases_equiv_diff))

  cat("\n--- Partition C (déterministe) ---\n")
  stopifnot(check_summary_equivalence_cov_diff(partC, nodesC, cases_equiv_diff))

  cat("=== Phase 2 OK ===\n")
  invisible(TRUE)
}

# ======================================================================================
# Phase 3: Fits courts via erpm() — sanity check
# ======================================================================================

.with_warning_capture <- function(expr) {
  warnings <- character()
  val <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = warnings)
}

run_fit_cov_diff <- function(part, nodes, rhs, tag, quiet = FALSE) {
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT-COVDIFF %-12s] SKIP (erpm() indisponible)\n", tag))
    return(list(ok = NA, fit = NULL, coef = NA))
  }

  f <- as.formula(paste0("partition ~ ", rhs))
  environment(f) <- list2env(list(partition = part, nodes = nodes), parent = parent.frame())

  s_obs <- try(summary_on_erpm_translation(part, nodes, rhs), silent = TRUE)
  if (!inherits(s_obs, "try-error") && !isTRUE(quiet)) {
    cat(sprintf("[ERPM-FIT-COVDIFF %-12s] stat_observee=%s\n", tag, paste(format(s_obs), collapse=", ")))
  }

  res <- .with_warning_capture(
    try(erpm(f, eval.loglik = TRUE, verbose = FALSE, nodes = nodes), silent = TRUE)
  )
  fit   <- res$value
  warns <- res$warnings

  if (length(warns) && !isTRUE(quiet)) {
    cat(sprintf("[ERPM-FIT-COVDIFF %-12s] WARNINGS (%d):\n", tag, length(warns)))
    for (w in unique(warns)) cat("  - ", w, "\n", sep = "")
  }

  if (inherits(fit, "try-error")) {
    msg <- as.character(fit)
    cat(sprintf("[ERPM-FIT-COVDIFF %-12s] ERREUR: %s\n", tag, msg))
    return(list(ok = FALSE, fit = NULL, coef = NA))
  }

  cf <- try(stats::coef(fit), silent = TRUE)
  ok <- !(inherits(cf, "try-error")) && all(is.finite(cf))
  cat(sprintf("[ERPM-FIT-COVDIFF %-12s] coef finies: %s | coef=%s\n",
              tag, if (ok) "OK" else "KO",
              if (ok) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

  list(ok = ok, fit = fit, coef = cf)
}

run_phase3_erpm_fits_cov_diff <- function(quiet = FALSE) {
  cat("\n=== PHASE 3 : Fits erpm() (sanity check) pour cov_diff (avec squared_sizes) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  fits <- list(
    A_k2 = run_fit_cov_diff(
      partA, nodesA,
      "squared_sizes() + cov_diff('score', clique_size = 2)",
      "A_k2", quiet = quiet
    ),
    A_k2_norm = run_fit_cov_diff(
      partA, nodesA,
      "squared_sizes() + cov_diff('score', clique_size = 2, normalized = TRUE)",
      "A_k2_norm", quiet = quiet
    ),
    A_k2_glob = run_fit_cov_diff(
      partA, nodesA,
      "squared_sizes() + cov_diff('score', clique_size = 2, normalize = 'global')",
      "A_k2_glob", quiet = quiet
    ),
    B_k2 = run_fit_cov_diff(
      partB, nodesB,
      "squared_sizes() + cov_diff('score', clique_size = 2)",
      "B_k2", quiet = quiet
    )
  )

  ok    <- vapply(fits, function(x) isTRUE(x$ok), logical(1))
  n_ok  <- sum(ok, na.rm = TRUE)
  n_tot <- sum(!is.na(ok))
  cat(sprintf("\n=== Bilan fits erpm() cov_diff : %d / %d OK ===\n", n_ok, n_tot))

  if (!isTRUE(quiet)) {
    cat("\n=== Résumés des fits ERPM cov_diff réussis ===\n")
    for (nm in names(fits)) {
      fx <- fits[[nm]]
      if (isTRUE(fx$ok) && inherits(fx$fit, "ergm")) {
        cat(sprintf("\n--- Résumé fit %s ---\n", nm))
        print(summary(fx$fit))
      }
    }
  }

  if (n_ok < n_tot) stop(sprintf("Echec fits cov_diff: %d KO", n_tot - n_ok))
  invisible(fits)
}

# ======================================================================================
# Phase 4: MCMC multi-toggle probe (diagnostic)
# ======================================================================================

.run_mcmc_multitoggle_probe <- function(nw) {
  # Objectif: déclencher la MCMC, idéalement avec des propositions multi-toggles.
  # On force un prop "sparse" (souvent multi-toggle dans les backends ERGM),
  # mais ce n'est pas garanti sur tous les setups.
  #
  # Pour voir les traces côté C:
  #   - mettre DEBUG_COV_DIFF = 1 dans changestat_cov_diff.c
  #   - recompiler le package
  ctrl <- control.simulate.formula(
    MCMC.burnin   = 2000,
    MCMC.interval = 1,
    MCMC.prop     = ~ sparse
  )

  f <- nw ~ cov_diff("score", clique_size = 2)

  sim <- simulate(
    f,
    nsim    = 1,
    control = ctrl,
    verbose = TRUE
  )

  print(sim)
  invisible(sim)
}

run_phase4_mcmc_probe <- function() {
  cat("\n=== PHASE 4 : MCMC MULTI-TOGGLE PROBE (cov_diff) ===\n")
  cat("Objectif: voir passer des traces debug du changestat D_ (si activées côté C).\n")
  cat("Si tu ne vois rien:\n")
  cat("  - DEBUG_COV_DIFF=0 (normal) => aucune trace\n")
  cat("  - ou les propositions ne sont pas multi-toggle sur ton setup\n\n")

  # Un réseau “pas trop trivial”, sinon la MCMC peut être pauvre.
  nw_probe <- .erpm_build_bipartite_nw(partB, nodesB)
  .run_mcmc_multitoggle_probe(nw_probe)

  cat("\nSi '[cov_diff] MULTI-TOGGLE ntoggles=...' apparaît, migration multi-toggle OK.\n")
  invisible(TRUE)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: cov_diff (multi-toggle) ===\n")

if (isTRUE(RUN$phase1_summary_expected)) {
  run_phase1_summary_expected_cov_diff(quiet = isTRUE(RUN$quiet_phase1))
} else {
  cat("\n=== PHASE 1 : SUMMARY EXPECTED ===\nSKIP\n")
}

if (isTRUE(RUN$phase2_equiv_summary)) {
  run_phase2_summary_equiv_cov_diff(quiet = isTRUE(RUN$quiet_phase2))
} else {
  cat("\n=== PHASE 2 : EQUIV SUMMARY ===\nSKIP\n")
}

res_fits_cov_diff <- NULL
if (isTRUE(RUN$phase3_erpm_fit)) {
  res_fits_cov_diff <- run_phase3_erpm_fits_cov_diff(quiet = isTRUE(RUN$quiet_phase3))
} else {
  cat("\n=== PHASE 3 : ERPM FIT ===\nSKIP\n")
}

if (isTRUE(RUN$phase4_mcmc_probe)) {
  run_phase4_mcmc_probe()
} else {
  cat("\n=== PHASE 4 : MCMC PROBE ===\nSKIP\n")
}

if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()

cat("\nTous les tests cov_diff ont passé.\n")