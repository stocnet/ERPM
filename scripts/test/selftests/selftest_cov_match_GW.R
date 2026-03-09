# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cov_match_GW.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_match_GW`
# Auteur : Jérémie Chichignoud - Cub'itech
#
# But du fichier
#   - PHASE 1 (SUMMARY) : valider la statistique via summary(nw ~ cov_match_GW(...))
#                         contre des attentes analytiques explicites (normalized='none').
#   - PHASE 2 (EQUIV)   : valider summary(nw) == summary(ERPM-traduit) sur un panel de RHS.
#   - PHASE 3 (FIT)     : vérifier que erpm() construit un modèle et renvoie des coefs finis.
#   - PHASE 4 (MCMC)    : diagnostic "multi-toggle" (objectif: déclencher le D_ changestat).
#
# IMPORTANT (multi-toggle / D_CHANGESTAT_FN):
#   - Le terme est désormais compilé en D_CHANGESTAT_FN(d_cov_match_GW).
#   - InitErgmTerm.cov_match_GW doit renvoyer d_func=TRUE (sinon crash possible).
#   - La PHASE 4 tente de déclencher ntoggles>1; selon le proposal/contraintes,
#     ce n'est pas garanti. Si tu ne vois pas de traces "MULTI-TOGGLE" alors:
#       * active DEBUG_COV_MATCH_GW=1 côté C, recompile,
#       * augmente burnin / interval,
#       * et/ou utilise un proposal qui fait des moves multi-arêtes (selon ton setup).
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
  # Recompile=TRUE utile ici car on vient de basculer en D_ changestat.
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else{
  stop("Le fichier DESCRIPTION n'existe pas ou devtools n'est pas installé.")
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
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_cov_match_GW.log")
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
  if (file.exists("R/erpm.R")) {
    source("R/erpm.R", local = FALSE)
  } else {
    cat("[WARN] erpm()/build_bipartite_from_inputs indisponibles. Certaines étapes seront sautées.\n")
  }
}

# ======================================================================================
# Réglages de run (point clé)
# ======================================================================================
RUN <- list(
  phase1_summary = TRUE,
  phase2_equiv   = TRUE,
  phase3_fit     = TRUE,
  phase4_mcmc    = TRUE,

  quiet_phase1   = FALSE,
  quiet_phase2   = FALSE,
  quiet_phase3   = FALSE
)

# ======================================================================================
# Données de test
# ======================================================================================

partitions <- list(
  # A : 3 groupes (2,3,2) = 7 sommets
  A = c(1,1, 2,2,2, 3,3),

  # B : 6 groupes (5,5,4,4,3,3) = 24 sommets
  B = c(rep(1, 5), rep(2, 5), rep(3, 4), rep(4, 4), rep(5, 3), rep(6, 3)),

  # C : 3 groupes (2,2,3) = 7 sommets
  C = c(1,1, 2,2, 3,3,3)
)

.make_nodes <- function(part) {
  n <- length(part)
  set.seed(100 + n)
  data.frame(
    label  = paste0("N", seq_len(n)),
    sexe   = sample(c("F","H"), n, replace = TRUE),
    dept   = sample(c("IT","RH","V"), n, replace = TRUE, prob = c(0.4,0.3,0.3)),
    grade  = sample(c("G1","G2","G3"), n, replace = TRUE, prob = c(0.4,0.4,0.2)),
    score  = sample(1:5, n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# Cas déterministes (reprend tes patterns, inchangés)
partA  <- partitions$A
nodesA <- data.frame(
  label = LETTERS[1:length(partA)],
  sexe  = c("F","H", "H","F","F", "H","H"),
  dept  = c("RH","RH", "V","IT","IT", "RH","V"),
  stringsAsFactors = FALSE
)

partC  <- partitions$C
nodesC <- data.frame(
  label = paste0("C", seq_along(partC)),
  sexe  = c("F","H", "H","H", "F","H","F"),
  grade = c("G1","G2", "G1","G3", "G1","G1","G2"),
  stringsAsFactors = FALSE
)

# ======================================================================================
# Helpers réseau biparti et summaries
# ======================================================================================

.erpm_build_bipartite_nw <- function(part, nodes) {
  stopifnot(length(part) == nrow(nodes))
  attrs <- as.list(nodes[, setdiff(names(nodes), "label"), drop=FALSE])

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
  as.numeric(suppressMessages(summary(f, constraints = ~ b1part)))
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
# Attentes analytiques pour cov_match_GW (sans normalisation)
# ======================================================================================

.S_GW <- function(n, lambda) {
  r <- (lambda - 1) / lambda
  lambda * (1 - r^n)
}

expected_cov_match_GW <- function(part, vals, lambda, category = NULL) {
  gid <- as.integer(part)
  split_idx <- split(seq_along(gid), gid)

  lambda <- as.numeric(lambda)
  out <- numeric(length(lambda))

  for (h in seq_along(lambda)) {
    lam <- lambda[h]
    tot <- 0
    for (ix in split_idx) {
      v <- vals[ix]
      if (is.null(category)) {
        tbl <- table(v, useNA = "no")
        for (n in as.integer(tbl)) {
          if (n > 0) tot <- tot + .S_GW(n, lam)
        }
      } else {
        n_k <- sum(v == category, na.rm = TRUE)
        if (n_k > 0) tot <- tot + .S_GW(n_k, lam)
      }
    }
    out[h] <- tot
  }
  out
}

# ======================================================================================
# Phase 1: Summary — attentes explicites sans normalisation
# ======================================================================================

run_phase1_summary_expected <- function() {
  cat("=== PHASE 1 : Summary cov_match_GW avec attentes explicites (normalized='none') ===\n")

  exp_A_l2 <- expected_cov_match_GW(partA, nodesA$sexe, lambda = 2)
  s_A_l2   <- summary_on_bipartite_network(partA, nodesA, "cov_match_GW('sexe', lambda = 2)")
  cat(sprintf("[A] cov_match_GW('sexe',λ=2)          obtenu=%g  attendu=%g\n", s_A_l2, exp_A_l2))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_l2), as.numeric(exp_A_l2))))

  exp_A_l3 <- expected_cov_match_GW(partA, nodesA$sexe, lambda = 3)
  s_A_l3   <- summary_on_bipartite_network(partA, nodesA, "cov_match_GW('sexe', lambda = 3)")
  cat(sprintf("[A] cov_match_GW('sexe',λ=3)          obtenu=%g  attendu=%g\n", s_A_l3, exp_A_l3))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_l3), as.numeric(exp_A_l3))))

  exp_A_l_vec <- expected_cov_match_GW(partA, nodesA$sexe, lambda = c(1.5, 2, 4))
  s_A_l_vec   <- summary_on_bipartite_network(partA, nodesA, "cov_match_GW('sexe', lambda = c(1.5, 2, 4))")
  cat(sprintf("[A] cov_match_GW('sexe',λ=c(1.5,2,4)) obtenu=%s  attendu=%s\n",
              paste(format(s_A_l_vec), collapse=", "),
              paste(format(exp_A_l_vec), collapse=", ")))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_l_vec), as.numeric(exp_A_l_vec))))

  exp_A_F_l2 <- expected_cov_match_GW(partA, nodesA$sexe, lambda = 2, category = "F")
  s_A_F_l2   <- summary_on_bipartite_network(partA, nodesA, "cov_match_GW('sexe', lambda = 2, category = 'F')")
  cat(sprintf("[A] cov_match_GW('sexe==F',λ=2)       obtenu=%g  attendu=%g\n", s_A_F_l2, exp_A_F_l2))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_F_l2), as.numeric(exp_A_F_l2))))

  exp_C_abs <- expected_cov_match_GW(partC, nodesC$grade, lambda = 2, category = "G999")
  s_C_abs   <- summary_on_bipartite_network(partC, nodesC, "cov_match_GW('grade', lambda = 2, category = 'G999')")
  cat(sprintf("[C] cov_match_GW('grade==G999',λ=2)   obtenu=%g  attendu=%g\n", s_C_abs, exp_C_abs))
  stopifnot(isTRUE(all.equal(as.numeric(s_C_abs), as.numeric(exp_C_abs))))

  cat("\n=== Phase 1 OK ===\n")
  invisible(NULL)
}

# ======================================================================================
# Phase 2: Summary — équivalence réseau explicite vs ERPM-traduit
# ======================================================================================

cases_equiv_GW <- c(
  "cov_match_GW('sexe', lambda = 2)",
  "cov_match_GW('sexe', lambda = c(1.5, 2, 4))",
  "cov_match_GW('sexe', lambda = 2, category = 'F')",
  "cov_match_GW('sexe', lambda = 2, normalized = 'by_group')",
  "cov_match_GW('sexe', lambda = 2, normalized = 'global')",
  "cov_match_GW('grade', lambda = 2, category = 'G1', normalized = 'by_group')"
)

check_summary_equivalence_GW <- function(part, nodes, rhs_vec) {
  for (rhs in rhs_vec) {
    s_net  <- summary_on_bipartite_network(part, nodes, rhs)
    s_erpm <- summary_on_erpm_translation(part, nodes, rhs)
    cat(sprintf("[EQUIV-GW] n=%-3d RHS=%-60s net=%s | erpm=%s\n",
                length(part), rhs, paste(s_net, collapse=","), paste(s_erpm, collapse=",")))
    if (length(s_net) != length(s_erpm)) stop("Longueur de statistique différente.")
    if (!all(is.finite(s_net)) || !all(is.finite(s_erpm))) stop("Stat non finie.")
    if (!isTRUE(all.equal(as.numeric(s_net), as.numeric(s_erpm)))) {
      stop(sprintf("Mismatch summary net vs ERPM pour RHS=%s", rhs))
    }
  }
  TRUE
}

run_phase2_summary_equiv_GW <- function() {
  cat("\n=== PHASE 2 : Summary(nw) vs Summary(ERPM-traduit) pour cov_match_GW ===\n")
  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes(part)
    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))
    ok <- check_summary_equivalence_GW(part, nodes, cases_equiv_GW)
    if (!ok) stop("Equivalence summary échouée.")
  }

  cat("\n--- Partition A (déterministe) ---\n")
  stopifnot(check_summary_equivalence_GW(partA, nodesA, cases_equiv_GW[c(1,2,3,4,5)]))

  cat("\n--- Partition C (déterministe) ---\n")
  stopifnot(check_summary_equivalence_GW(partC, nodesC, cases_equiv_GW[c(1,3,6)]))

  cat("=== Phase 2 OK ===\n")
  invisible(NULL)
}

# ======================================================================================
# Phase 3: Fits courts via erpm() — MLE + loglik
# ======================================================================================

ctrl_mle <- control.ergm(
  MCMLE.maxit     = 20,
  MCMC.burnin     = 10000,
  MCMC.interval   = 2000,
  MCMC.samplesize = 10000,
  parallel        = 0
)

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

run_fit_GW <- function(part, nodes, rhs, tag) {
  if (!exists("erpm", mode = "function")) {
    cat(sprintf("[ERPM-FIT-GW %-12s] SKIP (erpm() indisponible)\n", tag))
    return(list(ok = NA, fit = NULL, coef = NA))
  }

  f <- as.formula(paste0("partition ~ ", rhs))
  environment(f) <- list2env(list(partition = part, nodes = nodes), parent = parent.frame())

  s_obs <- try(summary_on_erpm_translation(part, nodes, rhs), silent = TRUE)
  if (!inherits(s_obs, "try-error")) {
    cat(sprintf("[ERPM-FIT-GW %-12s] stat_observee=%s\n", tag, paste(format(s_obs), collapse=", ")))
  } else {
    cat(sprintf("[ERPM-FIT-GW %-12s] stat_observee=NA (%s)\n", tag,
                conditionMessage(attr(s_obs, "condition"))))
  }

  res <- .with_warning_capture(
    try(erpm(f,
             eval.loglik = TRUE,
             # estimate = "MLE",
             # control  = ctrl_mle,
             verbose = FALSE,
             nodes = nodes),
        silent = TRUE)
  )

  fit   <- res$value
  warns <- res$warnings
  if (length(warns)) {
    cat(sprintf("[ERPM-FIT-GW %-12s] WARNINGS (%d):\n", tag, length(warns)))
    for (w in unique(warns)) cat("  - ", w, "\n", sep = "")
  }

  if (inherits(fit, "try-error")) {
    msg <- as.character(fit)
    cat(sprintf("[ERPM-FIT-GW %-12s] ERREUR: %s\n", tag, msg))
    return(list(ok = FALSE, fit = NULL, coef = NA))
  }

  cf <- try(stats::coef(fit), silent = TRUE)
  ok <- !(inherits(cf, "try-error")) && all(is.finite(cf))
  cat(sprintf("[ERPM-FIT-GW %-12s] coef finies: %s | coef=%s\n",
              tag, if (ok) "OK" else "KO",
              if (ok) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

  list(ok = ok, fit = fit, coef = cf)
}

run_phase3_erpm_fits_GW <- function() {
  cat("\n=== PHASE 3 : Fits erpm() (MLE + loglik) pour cov_match_GW ===\n")

  fits <- list(
    A_l2    = run_fit_GW(partA, nodesA, "cov_match_GW('sexe', lambda = 2)",                 "A_l2"),
    A_l24   = run_fit_GW(partA, nodesA, "cov_match_GW('sexe', lambda = c(2,4))",            "A_l24"),
    A_l2_F  = run_fit_GW(partA, nodesA, "cov_match_GW('sexe', lambda = 2, category = 'F')", "A_l2_F"),
    B_l23_H = run_fit_GW(partitions$B, .make_nodes(partitions$B),
                         "cov_match_GW('sexe', lambda = c(2,3), category = 'H')",          "B_l23_H")
  )

  ok    <- vapply(fits, function(x) if (is.na(x$ok)) NA else isTRUE(x$ok), logical(1))
  n_ok  <- sum(ok, na.rm = TRUE)
  n_tot <- sum(!is.na(ok))
  cat(sprintf("\n=== Bilan fits erpm() cov_match_GW : %d / %d OK ===\n", n_ok, n_tot))

  if (!isTRUE(RUN$quiet_phase3)) {
    cat("\n=== Résumés des fits ERPM cov_match_GW réussis ===\n")
    for (nm in names(fits)) {
      fx <- fits[[nm]]
      if (isTRUE(fx$ok) && inherits(fx$fit, "ergm")) {
        cat(sprintf("\n--- Résumé fit %s ---\n", nm))
        print(summary(fx$fit))
      }
    }
  }

  if (n_ok < n_tot) stop(sprintf("Echec fits cov_match_GW: %d KO", n_tot - n_ok))
  invisible(fits)
}

# ======================================================================================
# Phase 4: MCMC multi-toggle probe (diagnostic)
# ======================================================================================

.run_mcmc_multitoggle_probe_GW <- function(nw, rhs_txt) {
  # Objectif:
  #   - déclencher le chemin simulate()/MCMC
  #   - observer les traces C si DEBUG_COV_MATCH_GW=1
  #   - surtout: voir passer la bannière "MULTI-TOGGLE ntoggles=..."
  #
  # Notes:
  #   - Le fait d'avoir ntoggles>1 dépend du proposal/contraintes.
  #   - Le setup ci-dessous est volontairement minimal et peut nécessiter
  #     ajustement selon ton environnement ergm/ERPM.

  ctrl <- control.simulate.formula(
    MCMC.burnin   = 2000,
    MCMC.interval = 1,
    MCMC.prop     = ~ sparse
  )

  f <- .formula_nw(nw, rhs_txt)

  sim <- simulate(
    f,
    nsim        = 1,
    constraints = ~ b1part,
    control     = ctrl,
    verbose     = TRUE
  )

  print(sim)
  invisible(sim)
}

run_phase4_mcmc_probe_GW <- function() {
  cat("\n=== PHASE 4 : MCMC MULTI-TOGGLE PROBE (cov_match_GW) ===\n")
  cat("Objectif: voir passer les traces debug du changestat D_ (multi-toggle).\n")
  cat("Si tu ne vois pas 'MULTI-TOGGLE', ce n'est pas forcément KO: ça dépend\n")
  cat("du proposal (ntoggles peut rester à 1). Active DEBUG_COV_MATCH_GW=1 côté C\n")
  cat("et/ou change MCMC.prop si besoin.\n\n")

  nw_probe <- .erpm_build_bipartite_nw(partA, nodesA)

  # Un RHS simple, non-normalisé, pour éviter les effets numériques instables.
  .run_mcmc_multitoggle_probe_GW(nw_probe, "cov_match_GW('sexe', lambda = 2, normalized = 'none')")

  cat("\nSi '[cov_match_GW] MULTI-TOGGLE ntoggles=...' apparaît, probe multi-toggle OK.\n")
  invisible(TRUE)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: cov_match_GW ===\n")

if (isTRUE(RUN$phase1_summary)) {
  run_phase1_summary_expected()
} else {
  cat("\n=== PHASE 1 : SUMMARY ===\nSKIP\n")
}

if (isTRUE(RUN$phase2_equiv)) {
  run_phase2_summary_equiv_GW()
} else {
  cat("\n=== PHASE 2 : EQUIV ===\nSKIP\n")
}

res_fits_GW <- NULL
if (isTRUE(RUN$phase3_fit)) {
  res_fits_GW <- run_phase3_erpm_fits_GW()
} else {
  cat("\n=== PHASE 3 : FIT ===\nSKIP\n")
}

if (isTRUE(RUN$phase4_mcmc)) {
  run_phase4_mcmc_probe_GW()
} else {
  cat("\n=== PHASE 4 : MCMC MULTI-TOGGLE PROBE ===\nSKIP\n")
}

if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()

cat("\nTous les tests cov_match_GW ont passé.\n")