# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cov_fulldiff.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_fulldiff` (MULTI-TOGGLE)
# Exécution: Rscript scripts/test/selftests/selftest_cov_fulldiff.R
#
# But du fichier
#   - PHASE 1 (SUMMARY) : valider la statistique via summary(nw ~ cov_fulldiff(...))
#                         contre une référence analytique (range max-min par groupe).
#   - PHASE 2 (EQUIV)   : valider l'équivalence "réseau explicite" vs "ERPM-traduit".
#   - PHASE 3 (MCMC)    : déclencher le chemin MULTI-TOGGLE côté C (D_CHANGESTAT_FN)
#                         et observer les traces debug quand ntoggles > 1.
#
# Important (multi-toggle / D_CHANGESTAT_FN)
#   - cov_fulldiff est désormais implémenté côté C en D_CHANGESTAT_FN (multi-toggle).
#   - Le test PHASE 3 vise à vérifier qu'une proposition contrainte (b1part) produit
#     bien des moves à plusieurs toggles (typiquement 2: delete+add).
#   - Pour voir des traces, activer le debug côté C :
#       #define DEBUG_COV_FULLDIFF 1
#     puis recompiler (load_all(recompile=TRUE) suffit en général).
#
# Notes pratiques
#   - Les phases 1/2 peuvent spammer la console (summaries + tableaux).
#   - Pour bosser "proprement" sur la phase 3, désactiver 1/2 via RUN$phaseX.
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

# Charge le package (et recompile le code natif si nécessaire)
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else{
  stop("Le fichier DESCRIPTION n'existe pas ou devtools n'est pas installé.")
}

# --------------------------------------------------------------------------------------
# Logging local (identique esprit à tes autres selftests)
# --------------------------------------------------------------------------------------
root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_cov_fulldiff.log")
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
# Chargement utilitaires ERPM / wrapper (fallbacks)
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
# Réglages de run
# ======================================================================================
RUN <- list(
  phase1_summary = FALSE,
  phase2_equiv   = TRUE,
  phase3_mcmc    = FALSE,

  quiet_phase1   = FALSE,
  quiet_phase2   = FALSE
)

# ======================================================================================
# Données de test pour cov_fulldiff
# ======================================================================================

# PartA : 4 groupes (5,8,6,5) = 24 acteurs
partA  <- c(rep(1, 5), rep(2, 8), rep(3, 6), rep(4, 5))
nodesA <- data.frame(
  label = paste0("A", seq_along(partA)),
  score = c(
    10, 11,  9, 10, 12,          # g1
     5,  7, 20, 15,  9,  6, 30, 12,# g2
     0,  2,  4,  6,  8, 10,       # g3
    15, 16, 14, 17, 15           # g4
  ),
  stringsAsFactors = FALSE
)

# PartB : 6 groupes (6,6,5,4,3,2) = 26 acteurs
partB <- c(rep(1, 6), rep(2, 6), rep(3, 5), rep(4, 4), rep(5, 3), rep(6, 2))
nodesB <- data.frame(
  label = paste0("B", seq_along(partB)),
  score = c(
    10, 10, 11,  9, 10, 10,      # g1
     0, 40, 80,  5, 60,100,      # g2
     3,  4,  5,  6,  7,          # g3
    20, 21, 19, 20,              # g4
     5,  7,  6,                  # g5
     0, 10                       # g6
  ),
  stringsAsFactors = FALSE
)

# PartC : 4 groupes (5,7,8,4) = 24 acteurs
partC <- c(rep(1, 5), rep(2, 7), rep(3, 8), rep(4, 4))
nodesC <- data.frame(
  label = paste0("C", seq_along(partC)),
  score = c(
     5, 5, 5, 5, 5,              # g1 constant
     1, 2, 4, 3, 5, 6, 4,         # g2
     0, 0, 0, 0, 0, 0, 0,100,     # g3 outlier
     7, 8, 6, 7                   # g4
  ),
  stringsAsFactors = FALSE
)

partitions <- list(A = partA, B = partB, C = partC)

.make_nodes_numeric <- function(part) {
  n <- length(part)
  set.seed(100 + n)
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
  attrs <- as.list(nodes[, setdiff(names(nodes), "label"), drop = FALSE])

  if (exists("build_bipartite_from_inputs", mode = "function")) {
    builder <- get("build_bipartite_from_inputs")

    # Plusieurs signatures possibles selon tes versions du wrapper.
    candidates <- list(
      quote(builder(partition = part, nodes = nodes)),
      quote(builder(partition = part, labels = nodes$label, attributes = attrs)),
      quote(builder(partition = part)),
      quote(builder(partition = part, labels = nodes$label))
    )

    last_err <- NULL
    for (expr in candidates) {
      out <- try(eval(expr), silent = TRUE)
      if (inherits(out, "try-error")) { last_err <- out; next }

      if (inherits(out, "network")) return(out)

      if (is.list(out)) {
        for (nm in c("network","nw","g","graph","bip","net")) {
          if (!is.null(out[[nm]]) && inherits(out[[nm]], "network")) return(out[[nm]])
        }
      }
    }

    fm <- try(formals(builder), silent = TRUE)
    stop(paste0(
      "build_bipartite_from_inputs a échoué.\n",
      "Formals: ", if (!inherits(fm, "try-error")) paste(names(fm), collapse = ", ") else "<inconnus>", "\n",
      "Dernière erreur: ", if (!is.null(last_err)) as.character(last_err)[1] else "<aucune>"
    ))
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

summary_on_bipartite_network <- function(part, nodes, rhs_txt, constraints = NULL) {
  nw <- .erpm_build_bipartite_nw(part, nodes)
  f  <- .formula_nw(nw, rhs_txt)
  as.numeric(suppressMessages(summary(f, constraints = constraints)))
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
# Attentes analytiques pour cov_fulldiff
# ======================================================================================

# Pour chaque groupe g de taille n_g >= 2 :
#   D_g = max(x_i : i dans g) - min(x_i : i dans g)
# Stat T = sum_g D_g, éventuellement filtrée par un ensemble de tailles S.
expected_cov_fulldiff <- function(part, x, size = NULL) {
  gid <- as.integer(part)
  split_idx <- split(seq_along(gid), gid)
  tot <- 0
  for (ix in split_idx) {
    ng <- length(ix)
    if (ng < 2L) next
    if (!is.null(size) && !(ng %in% size)) next
    v <- x[ix]
    rng <- range(v)
    tot <- tot + (rng[2L] - rng[1L])
  }
  tot
}

# ======================================================================================
# Phase 1: Summary — attentes explicites (analytique vs summary)
# ======================================================================================

run_phase1_summary_expected_fulldiff <- function() {
  cat("=== PHASE 1 : Summary cov_fulldiff avec attentes explicites ===\n")

  exp_A_all <- expected_cov_fulldiff(partA, nodesA$score, size = NULL)
  s_A_all   <- summary_on_bipartite_network(partA, nodesA, "cov_fulldiff('score')")
  cat(sprintf("[A] cov_fulldiff('score')             obtenu=%g  attendu=%g\n", s_A_all, exp_A_all))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_all), as.numeric(exp_A_all))))

  exp_A_mid <- expected_cov_fulldiff(partA, nodesA$score, size = c(5,6,8))
  s_A_mid   <- summary_on_bipartite_network(partA, nodesA, "cov_fulldiff('score', size = c(5,6,8))")
  cat(sprintf("[A] cov_fulldiff('score',S={5,6,8})   obtenu=%g  attendu=%g\n", s_A_mid, exp_A_mid))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_mid), as.numeric(exp_A_mid))))

  exp_A_8   <- expected_cov_fulldiff(partA, nodesA$score, size = 8)
  s_A_8     <- summary_on_bipartite_network(partA, nodesA, "cov_fulldiff('score', size = 8)")
  cat(sprintf("[A] cov_fulldiff('score',S={8})       obtenu=%g  attendu=%g\n", s_A_8, exp_A_8))
  stopifnot(isTRUE(all.equal(as.numeric(s_A_8), as.numeric(exp_A_8))))

  exp_B_all <- expected_cov_fulldiff(partB, nodesB$score, size = NULL)
  s_B_all   <- summary_on_bipartite_network(partB, nodesB, "cov_fulldiff('score')")
  cat(sprintf("[B] cov_fulldiff('score')             obtenu=%g  attendu=%g\n", s_B_all, exp_B_all))
  stopifnot(isTRUE(all.equal(as.numeric(s_B_all), as.numeric(exp_B_all))))

  exp_B_6   <- expected_cov_fulldiff(partB, nodesB$score, size = 6)
  s_B_6     <- summary_on_bipartite_network(partB, nodesB, "cov_fulldiff('score', size = 6)")
  cat(sprintf("[B] cov_fulldiff('score',S={6})       obtenu=%g  attendu=%g\n", s_B_6, exp_B_6))
  stopifnot(isTRUE(all.equal(as.numeric(s_B_6), as.numeric(exp_B_6))))

  exp_B_56  <- expected_cov_fulldiff(partB, nodesB$score, size = c(5,6))
  s_B_56    <- summary_on_bipartite_network(partB, nodesB, "cov_fulldiff('score', size = c(5,6))")
  cat(sprintf("[B] cov_fulldiff('score',S={5,6})     obtenu=%g  attendu=%g\n", s_B_56, exp_B_56))
  stopifnot(isTRUE(all.equal(as.numeric(s_B_56), as.numeric(exp_B_56))))

  exp_C_all <- expected_cov_fulldiff(partC, nodesC$score, size = NULL)
  s_C_all   <- summary_on_bipartite_network(partC, nodesC, "cov_fulldiff('score')")
  cat(sprintf("[C] cov_fulldiff('score')             obtenu=%g  attendu=%g\n", s_C_all, exp_C_all))
  stopifnot(isTRUE(all.equal(as.numeric(s_C_all), as.numeric(exp_C_all))))

  exp_C_ge6 <- expected_cov_fulldiff(partC, nodesC$score, size = 6:20)
  s_C_ge6   <- summary_on_bipartite_network(partC, nodesC, "cov_fulldiff('score', size = 6:20)")
  cat(sprintf("[C] cov_fulldiff('score',S=6:20)      obtenu=%g  attendu=%g\n", s_C_ge6, exp_C_ge6))
  stopifnot(isTRUE(all.equal(as.numeric(s_C_ge6), as.numeric(exp_C_ge6))))

  cat("\n=== Phase 1 OK ===\n")
  invisible(NULL)
}

# ======================================================================================
# Phase 2: Summary — équivalence réseau explicite vs ERPM-traduit
# ======================================================================================

cases_equiv_fulldiff <- c(
  "cov_fulldiff('score')",
  "cov_fulldiff('score', size = 2)",
  "cov_fulldiff('score', size = c(4,5,6,8))",
  "cov_fulldiff('score', size = 6:20)"
)

check_summary_equivalence_fulldiff <- function(part, nodes, rhs_vec) {
  for (rhs in rhs_vec) {
    s_net  <- summary_on_bipartite_network(part, nodes, rhs)
    s_erpm <- summary_on_erpm_translation(part, nodes, rhs)

    cat(sprintf("[EQUIV-FULLDIFF] n=%-3d RHS=%-40s net=%s | erpm=%s\n",
                length(part), rhs, paste(s_net, collapse=","), paste(s_erpm, collapse=",")))

    if (length(s_net) != length(s_erpm)) stop("Longueur de statistique différente.")
    if (!all(is.finite(s_net)) || !all(is.finite(s_erpm))) stop("Stat non finie.")
    if (!isTRUE(all.equal(as.numeric(s_net), as.numeric(s_erpm)))) {
      stop(sprintf("Mismatch summary net vs ERPM pour RHS=%s", rhs))
    }
  }
  TRUE
}

run_phase2_summary_equiv_fulldiff <- function() {
  cat("\n=== PHASE 2 : Summary(nw) vs Summary(ERPM-traduit) pour cov_fulldiff ===\n")

  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes_numeric(part)

    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))

    ok <- check_summary_equivalence_fulldiff(part, nodes, cases_equiv_fulldiff)
    if (!ok) stop("Equivalence summary échouée.")
  }

  cat("\n--- Partition A (déterministe) ---\n")
  stopifnot(check_summary_equivalence_fulldiff(partA, nodesA, cases_equiv_fulldiff[c(1,3,4)]))

  cat("\n--- Partition B (déterministe) ---\n")
  stopifnot(check_summary_equivalence_fulldiff(partB, nodesB, cases_equiv_fulldiff[c(1,2,3)]))

  cat("\n--- Partition C (déterministe) ---\n")
  stopifnot(check_summary_equivalence_fulldiff(partC, nodesC, cases_equiv_fulldiff[c(1,3,4)]))

  cat("=== Phase 2 OK ===\n")
  invisible(NULL)
}

# ======================================================================================
# Phase 3: MCMC multi-toggle probe (objectif: ntoggles > 1 côté C)
# ======================================================================================

# Le but ici n'est PAS de faire un fit MLE propre, juste:
#   - déclencher la MCMC (simulate)
#   - sous contrainte b1part, provoquer des moves "réaffectation" (2 toggles)
#   - observer les traces debug dans la console si DEBUG_COV_FULLDIFF=1 côté C.
.run_mcmc_multitoggle_probe <- function(nw) {
  # On préfère des settings courts et bruyants (verbose=TRUE).
  ctrl <- control.simulate.formula(
    MCMC.burnin   = 2000,
    MCMC.interval = 1,
    MCMC.prop     = ~ b1part
  )

  f <- nw ~ cov_fulldiff("score")

  # NB: la contrainte b1part est cruciale pour obtenir des moves multi-toggle
  # (swap: supprimer ancienne appartenance + ajouter nouvelle).
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

run_phase3_mcmc_multitoggle <- function() {
  cat("\n=== PHASE 3 : MCMC MULTI-TOGGLE PROBE (cov_fulldiff) ===\n")
  cat("Objectif: voir passer les traces debug du changestat D_ quand ntoggles>1.\n")
  cat("Pré-requis: recompiler avec DEBUG_COV_FULLDIFF=1 dans changestat_cov_fulldiff.c\n")
  cat("et utiliser une contrainte b1part.\n\n")

  # Réseau simple mais non-trivial: on reprend A (déterministe)
  nw <- .erpm_build_bipartite_nw(partA, nodesA)

  # Important: la covariate doit exister sur le mode acteur (les n1 premiers sommets).
  # Selon ton builder, l'attribut peut déjà être présent; sinon on le force ici.
  n1 <- as.integer(nw %n% "bipartite")
  if (is.na(n1) || n1 <= 0L) stop("Réseau non biparti (bipartite manquant).")

  # Si l'attribut n'existe pas, on le pose proprement:
  v <- network::get.vertex.attribute(nw, "score")
  if (is.null(v)) {
    # convention: attribut défini sur tous les sommets; NA sur le mode groupe
    network::set.vertex.attribute(nw, "score", c(nodesA$score, rep(NA_real_, network::network.size(nw) - n1)))
  }

  .run_mcmc_multitoggle_probe(nw)

  cat("\nSi tu vois dans la console des lignes du style:\n")
  cat("  [cov_fulldiff] MULTI-TOGGLE ntoggles=2\n")
  cat("alors le test multi-toggle est OK.\n")

  invisible(TRUE)
}

# ======================================================================================
# Exécution
# ======================================================================================
set.seed(1)
cat("=== TEST ERPM: cov_fulldiff (MULTI-TOGGLE) ===\n")

if (isTRUE(RUN$phase1_summary)) {
  if (isTRUE(RUN$quiet_phase1)) cat("[PHASE 1] mode quiet\n")
  run_phase1_summary_expected_fulldiff()
} else {
  cat("\n=== PHASE 1 : SUMMARY ===\nSKIP\n")
}

if (isTRUE(RUN$phase2_equiv)) {
  if (isTRUE(RUN$quiet_phase2)) cat("[PHASE 2] mode quiet\n")
  run_phase2_summary_equiv_fulldiff()
} else {
  cat("\n=== PHASE 2 : EQUIV ===\nSKIP\n")
}

if (isTRUE(RUN$phase3_mcmc)) {
  run_phase3_mcmc_multitoggle()
} else {
  cat("\n=== PHASE 3 : MCMC MULTI-TOGGLE ===\nSKIP\n")
}

if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()

cat("\nTous les tests cov_fulldiff ont passé.\n")