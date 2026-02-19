# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cov_match.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_match` (multi-toggle)
# Exécution: Rscript scripts/test/selftests/selftest_cov_match.R
#
# But du fichier
#   - PHASE 1 (SUMMARY) : valider la statistique via summary(nw ~ cov_match(...)).
#   - PHASE 2 (EQUIV)   : valider summary(nw) == summary(ERPM-traduit).
#   - PHASE 3 (FIT)     : valider que erpm() construit un modèle et renvoie des coefs finis.
#   - PHASE 4 (MCMC)    : diagnostic "multi-toggle" (le seul truc qui nous intéresse quand
#                         on active le debug côté C et qu'on ne veut pas se faire noyer).
#
# Important (multi-toggle / D_CHANGESTAT_FN)
#   - Ce selftest suppose que le changestat C de cov_match est compilé en D_ (multi-toggle),
#     i.e. D_CHANGESTAT_FN(d_cov_match), et que InitErgmTerm.cov_match annonce d_func=TRUE.
#   - Si tu vois des segfaults ou "C function not found", c’est typiquement un mismatch
#     entre C_ vs D_ (ou un symbole mal nommé).
#
# Important (Phase 4)
#   - La Phase 4 ne "prouve" pas formellement que ergm propose toujours des multi-toggles :
#     elle sert de PROBE. On cherche à voir passer, dans la console, des traces de debug
#     du changestat côté C pour des propositions ntoggles>1.
#   - Pour voir ces traces, il faut compiler avec DEBUG_COV_MATCH=1 dans changestat_cov_match.c,
#     puis recompiler/recharger le package.
# ======================================================================================

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
} else{
  stop("Exécuter depuis la racine du package (DESCRIPTION) avec devtools disponible.")
}

# --------------------------------------------------------------------------------------
# Logging local
# --------------------------------------------------------------------------------------
root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_cov_match.log")
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
# Réglages de run (point clé)
# ======================================================================================
# Objectif: pouvoir isoler la phase MCMC multitoggle sans éditer 50 endroits.
RUN <- list(
  phase1_summary_expected = FALSE,
  phase2_summary_equiv    = FALSE,
  phase3_erpm_fits        = TRUE,
  phase4_mcmc_probe       = FALSE,

  quiet_phase1            = FALSE,
  quiet_phase2            = FALSE,
  quiet_phase3            = FALSE
)

# ======================================================================================
# Données de test
# ======================================================================================

# Panel de partitions variées
partitions <- list(
  A = c(1,1, 2,2,2, 3,3,3,3, 4),               # tailles: 2,3,4,1
  B = c(1,1,1, 2,2, 3,3,3, 4,4, 5),             # tailles: 3,2,3,2,1
  C = c(1,1,1, 2,2,2, 3, 4)                     # tailles: 3,3,1,1
)

# Génère un jeu d'attributs contrôlé par partition
.make_nodes <- function(part) {
  n <- length(part)
  set.seed(100 + n)
  data.frame(
    label  = paste0("N", seq_len(n)),
    sexe   = sample(c("F","H"), n, replace = TRUE),
    dept   = sample(c("IT","RH","V"), n, replace = TRUE, prob = c(0.4,0.3,0.3)),
    grade  = sample(c("G1","G2","G3"), n, replace = TRUE, prob = c(0.4,0.4,0.2)),
    score  = sample(1:5, n, replace = TRUE),  # numérique (doit être rejeté par cov_match)
    stringsAsFactors = FALSE
  )
}

# Cas spécifiques déterministes (reprennent tes exemples)
partA  <- c(1,1, 2,2,2, 3,3,3,3, 4)
nodesA <- data.frame(
  label = LETTERS[1:length(partA)],
  sexe  = c("F","F",  "H","F","H",  "H","H","H","F",  "H"),
  dept  = c("RH","RH","V","V","V",  "IT","IT","IT","IT","V"),
  stringsAsFactors = FALSE
)

partC  <- c(1,1,1, 2,2,2, 3, 4)
nodesC <- data.frame(
  label = paste0("C", seq_along(partC)),
  sexe  = c("F","H","H",  "H","H","H",  "F","H"),
  grade = c("G1","G2","G2","G1","G1","G3","G1","G2"),
  stringsAsFactors = FALSE
)

# Cas déterministe dédié au fit B_k2_RHbg (évite stat == 0 et constance)
partB_RHbg <- c(1,1,1, 2,2, 3,3)  # tailles: 3,2,2
nodesB_RHbg <- data.frame(
  label = paste0("BR", seq_along(partB_RHbg)),
  sexe  = c("F","H","F",  "H","H",  "F","H"),
  dept  = c("RH","RH","V",  "RH","V",  "V","V"),  # un seul groupe avec paires RH
  stringsAsFactors = FALSE
)

# ======================================================================================
# Helpers réseau biparti et summaries
# ======================================================================================

# Constructeur biparti robuste
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

# Fabrique formule nw ~ <rhs>
.formula_nw <- function(nw, rhs_txt) {
  f <- as.formula(paste0("nw ~ ", rhs_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  f
}

# Summary côté réseau biparti explicite
summary_on_bipartite_network <- function(part, nodes, rhs_txt, constraints = ~ b1part) {
  nw <- .erpm_build_bipartite_nw(part, nodes)
  f  <- .formula_nw(nw, rhs_txt)
  as.numeric(suppressMessages(summary(f, constraints = constraints)))
}

# Summary côté ERPM (traduction)
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
  cons <- call_args$constraints; if (is.null(cons)) cons <- as.formula(~ b1part)
  as.numeric(suppressMessages(summary(f2, constraints = cons)))
}

# ======================================================================================
# Attentes analytiques pour cov_match (réplique la sémantique implémentée)
# ======================================================================================

.choose <- function(n, k) if (n >= k) choose(n, k) else 0

# normalized = "global"
#   T_k(B; c) = sum_g S_k(g) / n_g  (ou C(n_{g,κ},k) / n_g en ciblé)
#   avec les groupes de taille n_g = 0 contribuant 0.
expected_cov_match <- function(part, vals, k = 2L, category = NULL,
                               normalized = c("none","by_group","global")) {
  normalized <- match.arg(normalized)
  gid <- as.integer(part)
  split_idx <- split(seq_along(gid), gid)

  stat_group <- function(ix) {
    v   <- vals[ix]
    n_g <- length(ix)

    if (!is.null(category)) {
      n_k <- sum(v == category, na.rm = TRUE)
      if (k == 1L && normalized == "by_group") {
        # Cas spécial: somme_g 1_{n_{g,κ} >= 1}
        return(as.numeric(n_k >= 1L))
      }
      clq <- .choose(n_k, k)
    } else {
      tbl <- table(v, useNA = "no")
      clq <- sum(vapply(as.integer(tbl), .choose, numeric(1), k = k))
    }

    if (normalized == "none") return(clq)

    if (normalized == "by_group") {
      den <- .choose(n_g, k)
      if (den == 0) return(0)
      return(clq / den)
    }

    den <- n_g
    if (den == 0) return(0)
    return(clq / den)
  }

  sum(vapply(split_idx, stat_group, numeric(1)))
}

# ======================================================================================
# Phase 1: Summary — attentes numériques explicites + erreurs attendues
# ======================================================================================

run_phase1_summary_expected <- function(quiet = FALSE) {
  cat("=== PHASE 1 : Summary avec attentes explicites ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  # 1) Cas A (déterministe)
  exp_A_k2 <- expected_cov_match(partA, nodesA$sexe, k = 2, normalized = "none")
  s_A_k2   <- summary_on_bipartite_network(partA, nodesA, "cov_match('sexe', clique_size = 2)")
  cat(sprintf("[A] cov_match('sexe',k=2)            obtenu=%g  attendu=%g\n", s_A_k2, exp_A_k2))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_k2), as.numeric(exp_A_k2))))

  exp_A_k3 <- expected_cov_match(partA, nodesA$sexe, k = 3, normalized = "none")
  s_A_k3   <- summary_on_bipartite_network(partA, nodesA, "cov_match('sexe', clique_size = 3)")
  cat(sprintf("[A] cov_match('sexe',k=3)            obtenu=%g  attendu=%g\n", s_A_k3, exp_A_k3))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_k3), as.numeric(exp_A_k3))))

  exp_A_F2 <- expected_cov_match(partA, nodesA$sexe, k = 2, category = "F", normalized = "none")
  s_A_F2   <- summary_on_bipartite_network(partA, nodesA, "cov_match('sexe', clique_size = 2, category = 'F')")
  cat(sprintf("[A] cov_match('sexe==F',k=2)         obtenu=%g  attendu=%g\n", s_A_F2, exp_A_F2))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_F2), as.numeric(exp_A_F2))))

  exp_A_bg2 <- expected_cov_match(partA, nodesA$sexe, k = 2, normalized = "by_group")
  s_A_bg2   <- summary_on_bipartite_network(partA, nodesA, "cov_match('sexe', clique_size = 2, normalized = 'by_group')")
  cat(sprintf("[A] cov_match('sexe',k=2,by_group)   obtenu=%.6f  attendu=%.6f\n", s_A_bg2, exp_A_bg2))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_bg2), as.numeric(exp_A_bg2))))

  exp_A_gl2 <- expected_cov_match(partA, nodesA$sexe, k = 2, normalized = "global")
  s_A_gl2   <- summary_on_bipartite_network(partA, nodesA, "cov_match('sexe', clique_size = 2, normalized = 'global')")
  cat(sprintf("[A] cov_match('sexe',k=2,global)     obtenu=%.7f  attendu=%.7f\n", s_A_gl2, exp_A_gl2))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_gl2), as.numeric(exp_A_gl2))))

  # 2) Cas C (k=1 autorisés en by_group)
  exp_C_k1_bg <- expected_cov_match(partC, nodesC$sexe, k = 1, normalized = "by_group")
  s_C_k1_bg   <- summary_on_bipartite_network(partC, nodesC, "cov_match('sexe', clique_size = 1, normalized = 'by_group')")
  cat(sprintf("[C] cov_match('sexe',k=1,by_group)   obtenu=%g  attendu=%g\n", s_C_k1_bg, exp_C_k1_bg))
  stopifnot(isTRUE(all.equal(as.numeric(s_C_k1_bg), as.numeric(exp_C_k1_bg))))

  exp_C_g1_bg <- expected_cov_match(partC, nodesC$grade, k = 1, category = "G1", normalized = "by_group")
  s_C_g1_bg   <- summary_on_bipartite_network(partC, nodesC, "cov_match('grade', clique_size = 1, category = 'G1', normalized = 'by_group')")
  cat(sprintf("[C] cov_match('grade==G1',k=1,bg)    obtenu=%g  attendu=%g\n", s_C_g1_bg, exp_C_g1_bg))
  stopifnot(isTRUE(all.equal(as.numeric(s_C_g1_bg), as.numeric(exp_C_g1_bg))))

  exp_C_abs <- expected_cov_match(partC, nodesC$grade, k = 2, category = "G999", normalized = "none")
  s_C_abs   <- summary_on_bipartite_network(partC, nodesC, "cov_match('grade', clique_size = 2, category = 'G999')")
  cat(sprintf("[C] cov_match('grade==G999',k=2)     obtenu=%g  attendu=%g\n", s_C_abs, exp_C_abs))
  stopifnot(isTRUE(all.equal(as.numeric(s_C_abs), as.numeric(exp_C_abs))))

  # --------------------------------------------------------------------
  # BLOC D'ERREURS ATTENDUES POUR k=1, normalized != 'by_group'
  # --------------------------------------------------------------------
  cat("\n--- Vérification des erreurs attendues pour k=1 ---\n")

  err <- NULL
  tryCatch({
    summary_on_bipartite_network(partC, nodesC, "cov_match('sexe', clique_size = 1)")
  }, error = function(e) err <<- conditionMessage(e))
  if (is.null(err)) {
    stop("Erreur attendue non levée pour cov_match('sexe', clique_size = 1) (normalized='none').")
  } else {
    cat("[OK] Erreur attendue pour k=1, normalized='none':\n     ", err, "\n")
  }

  err <- NULL
  tryCatch({
    summary_on_bipartite_network(partC, nodesC, "cov_match('sexe', clique_size = 1, normalized = 'global')")
  }, error = function(e) err <<- conditionMessage(e))
  if (is.null(err)) {
    stop("Erreur attendue non levée pour cov_match('sexe', clique_size = 1, normalized='global').")
  } else {
    cat("[OK] Erreur attendue pour k=1, normalized='global':\n     ", err, "\n")
  }

  cat("\n=== Phase 1 OK ===\n")
  invisible(TRUE)
}

# ======================================================================================
# Phase 2: Summary — équivalence réseau explicite vs ERPM-traduit
# ======================================================================================

cases_equiv <- c(
  "cov_match('sexe', clique_size = 2)",
  "cov_match('sexe', clique_size = c(2,3))",
  "cov_match('dept', clique_size = 2, category='RH')",
  "cov_match('dept', clique_size = 2, normalized = 'by_group')",
  "cov_match('grade', clique_size = 1, category='G1', normalized='by_group')"
)

check_summary_equivalence <- function(part, nodes, rhs_vec) {
  for (rhs in rhs_vec) {
    s_net  <- summary_on_bipartite_network(part, nodes, rhs)
    s_erpm <- summary_on_erpm_translation(part, nodes, rhs)
    cat(sprintf("[EQUIV] n=%-3d RHS=%-60s net=%s | erpm=%s\n",
                length(part), rhs, paste(s_net, collapse=","), paste(s_erpm, collapse=",")))
    if (length(s_net) != length(s_erpm)) stop("Longueur de statistique différente.")
    if (!all(is.finite(s_net)) || !all(is.finite(s_erpm))) stop("Stat non finie.")
    if (!isTRUE(all.equal(as.numeric(s_net), as.numeric(s_erpm)))) {
      stop(sprintf("Mismatch summary net vs ERPM pour RHS=%s", rhs))
    }
  }
  TRUE
}

run_phase2_summary_equiv <- function(quiet = FALSE) {
  cat("\n=== PHASE 2 : Summary(nw) vs Summary(ERPM-traduit) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes(part)
    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))
    ok <- check_summary_equivalence(part, nodes, cases_equiv)
    if (!ok) stop("Equivalence summary échouée.")
  }

  cat("\n--- Partition A (déterministe) ---\n")
  stopifnot(check_summary_equivalence(partA, nodesA, cases_equiv[1:3]))

  cat("\n--- Partition C (déterministe) ---\n")
  stopifnot(check_summary_equivalence(partC, nodesC, cases_equiv[c(1,3,5)]))

  cat("=== Phase 2 OK ===\n")
  invisible(TRUE)
}

# ======================================================================================
# Phase 3: Fits courts via erpm() — MLE + loglik
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

run_fit <- function(part, nodes, rhs, tag, quiet = FALSE) {
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT %-12s] SKIP (erpm() indisponible)\n", tag))
    return(list(ok = NA, fit = NULL, coef = NA))
  }

  f <- as.formula(paste0("partition ~ ", rhs))
  environment(f) <- list2env(list(partition = part, nodes = nodes), parent = parent.frame())

  s_obs <- try(summary_on_erpm_translation(part, nodes, rhs), silent = TRUE)
  if (!inherits(s_obs, "try-error") && !isTRUE(quiet)) {
    cat(sprintf("[ERPM-FIT %-12s] stat_observee=%s\n", tag, paste(format(s_obs), collapse=", ")))
  }

  res <- .with_warning_capture(
    try(erpm(f, eval.loglik = TRUE, verbose = FALSE, nodes = nodes), silent = TRUE)
  )
  fit   <- res$value
  warns <- res$warnings
  if (length(warns) && !isTRUE(quiet)) {
    cat(sprintf("[ERPM-FIT %-12s] WARNINGS (%d):\n", tag, length(warns)))
    for (w in unique(warns)) cat("  - ", w, "\n", sep = "")
  }

  if (inherits(fit, "try-error")) {
    msg <- as.character(fit)
    cat(sprintf("[ERPM-FIT %-12s] ERREUR: %s\n", tag, msg))
    return(list(ok = FALSE, fit = NULL, coef = NA))
  }

  cf <- try(stats::coef(fit), silent = TRUE)
  ok <- !(inherits(cf, "try-error")) && all(is.finite(cf))
  cat(sprintf("[ERPM-FIT %-12s] coef finies: %s | coef=%s\n",
              tag, if (ok) "OK" else "KO",
              if (ok) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

  list(ok = ok, fit = fit, coef = cf)
}

run_phase3_erpm_fits <- function(quiet = FALSE) {
  cat("\n=== PHASE 3 : Fits erpm() (MLE + loglik) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  fits <- list(
    A_k2       = run_fit(partA, nodesA,
                         "cov_match('sexe', clique_size = 2)", "A_k2", quiet = quiet),
    A_k2F      = run_fit(partA, nodesA,
                         "cov_match('sexe', clique_size = 2, category = 'F')", "A_k2F", quiet = quiet),
    A_k2_bg    = run_fit(partA, nodesA,
                         "cov_match('sexe', clique_size = 2, normalized = 'by_group')", "A_k2_bg", quiet = quiet),
    B_k23_H    = run_fit(partitions$B, .make_nodes(partitions$B),
                         "cov_match('sexe', clique_size = c(2,3), category = 'H')", "B_k23_H", quiet = quiet),
    B_k2_RHbg  = run_fit(partB_RHbg, nodesB_RHbg,
                         "cov_match('dept', clique_size = 2, category = 'RH', normalized = 'by_group')", "B_k2_RHbg", quiet = quiet)
  )

  ok <- vapply(fits, function(x) isTRUE(x$ok), logical(1))
  n_ok  <- sum(ok, na.rm = TRUE)
  n_tot <- sum(!is.na(ok))
  cat(sprintf("\n=== Bilan fits erpm() : %d / %d OK ===\n", n_ok, n_tot))

  if (n_ok < n_tot) stop(sprintf("Echec fits: %d KO", n_tot - n_ok))

  if (!isTRUE(quiet)) {
    cat("\n=== Résumés des fits ERPM réussis ===\n")
    for (nm in names(fits)) {
      fx <- fits[[nm]]
      if (isTRUE(fx$ok) && inherits(fx$fit, "ergm")) {
        cat(sprintf("\n--- Résumé fit %s ---\n", nm))
        print(summary(fx$fit))
      }
    }
  }

  invisible(fits)
}

# ======================================================================================
# Phase 4: MCMC multi-toggle probe (diagnostic)
# ======================================================================================

.run_mcmc_multitoggle_probe <- function(nw, rhs_txt) {
  # Objectif: déclencher le MCMC, et observer dans la console les traces debug du C
  # (si DEBUG_COV_MATCH=1 au build), notamment des indices que ntoggles>1 se produit.
  #
  # Notes:
  # - verbose=TRUE est voulu.
  # - MCMC.prop=~ sparse est un bon candidat pour produire des moves avec plusieurs toggles.
  ctrl <- control.simulate.formula(
    MCMC.burnin   = 1000,
    MCMC.interval = 1,
    MCMC.prop     = ~ sparse
  )

  f <- as.formula(paste0("nw ~ ", rhs_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())

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
  cat("\n=== PHASE 4 : MCMC MULTI-TOGGLE PROBE ===\n")
  cat("Objectif: voir passer des traces debug du changestat D_ (multi-toggle).\n")
  cat("Pré-requis: compiler changestat_cov_match.c avec DEBUG_COV_MATCH=1, puis recharger le package.\n\n")

  # On prend un cas simple qui bouge (éviter stat const/0).
  nw <- .erpm_build_bipartite_nw(partA, nodesA)

  # RHS volontairement simple (k=2) + by_group (stable numériquement)
  rhs <- "cov_match('sexe', clique_size = 2, normalized = 'by_group')"

  .run_mcmc_multitoggle_probe(nw, rhs_txt = rhs)

  cat("\nSi tu vois des traces style '[cov_match] MULTI-TOGGLE ntoggles=...' (côté C), c'est OK.\n")
  invisible(TRUE)
}

# ======================================================================================
# Run principal
# ======================================================================================

run_all_tests_cov_match <- function() {
  set.seed(1)

  cat("=== TEST ERPM: cov_match (multi-toggle) ===\n")
  cat("R:", paste(R.version$major, R.version$minor, sep="."), "\n")
  cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")

  # PHASE 1
  if (isTRUE(RUN$phase1_summary_expected)) {
    run_phase1_summary_expected(quiet = isTRUE(RUN$quiet_phase1))
  } else {
    cat("\n=== PHASE 1 : Summary ===\nSKIP\n")
  }

  # PHASE 2
  if (isTRUE(RUN$phase2_summary_equiv)) {
    run_phase2_summary_equiv(quiet = isTRUE(RUN$quiet_phase2))
  } else {
    cat("\n=== PHASE 2 : Equivalence ===\nSKIP\n")
  }

  # PHASE 3
  res_fits <- NULL
  if (isTRUE(RUN$phase3_erpm_fits)) {
    res_fits <- run_phase3_erpm_fits(quiet = isTRUE(RUN$quiet_phase3))
  } else {
    cat("\n=== PHASE 3 : Fits erpm() ===\nSKIP\n")
  }

  # PHASE 4
  if (isTRUE(RUN$phase4_mcmc_probe)) {
    run_phase4_mcmc_probe()
  } else {
    cat("\n=== PHASE 4 : MCMC probe ===\nSKIP\n")
  }

  invisible(list(fits = res_fits))
}

if (identical(environment(), globalenv())) {
  run_all_tests_cov_match()
}

# --------------------------------------------------------------------------------------
# Fin: on désactive le patch si présent
# --------------------------------------------------------------------------------------
if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()
cat("\nTous les tests cov_match ont passé.\n")