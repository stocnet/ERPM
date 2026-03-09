# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cov_fullmatch.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_fullmatch`
# Auteur : Jérémie Chichignoud - Cub'itech
#
# Objectifs (structure "squared_sizes"-style, mais plus compact)
#   - PHASE 1 (EXPECTED): valider la statistique via summary(nw ~ cov_fullmatch(...))
#                         avec des attentes numériques explicites (cas analytiques).
#   - PHASE 2 (EQUIV)   : vérifier l’équivalence summary(nw) vs summary(ERPM-traduit)
#                         sur un panel de partitions + covariées variées.
#   - PHASE 3 (FIT)     : vérifier qu’un fit via erpm() passe et renvoie un coef fini.
#   - PHASE 4 (MCMC)    : "probe" multi-toggle: déclencher le code MCMC et observer
#                         (si activé côté C) les traces indiquant ntoggles>1.
#
# Notes importantes
#   - cov_fullmatch doit compter le NOMBRE de groupes homogènes (et non la somme des tailles).
#   - Les phases 1/2 peuvent être bruyantes; on peut les désactiver via RUN.
#   - La PHASE 4 n’est pas un test statistique "fort": c’est un test de chemin d’exécution
#     (multitoggle) destiné à valider une implémentation D_CHANGESTAT_FN si vous l’avez.
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule
# --------------------------------------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",    quietly = TRUE)) stop("Package 'ergm' requis.")
  if (!requireNamespace("utils",   quietly = TRUE)) stop("Package 'utils' requis.")
})

suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

# Patch ERGM optionnel (si utilisé dans le projet)
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable", mode = "function")) {
    ergm_patch_enable()
  } else {
    message("[ergm_patch] fonction ergm_patch_enable absente")
  }
}

# Chargement package depuis la racine (DESCRIPTION)
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Exécuter depuis la racine du package (DESCRIPTION) avec devtools disponible.")
}

# --------------------------------------------------------------------------------------
# Logging local (fichier .log à côté du selftest)
# --------------------------------------------------------------------------------------
root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_cov_fullmatch.log")
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
# Chargement utilitaires ERPM (fallback)
# --------------------------------------------------------------------------------------
# Fallback si la fonction n'est pas exposée par le wrapper
if (!exists("partition_to_bipartite_network", mode = "function")) {
  if (file.exists("R/functions_erpm_bip_network.R")) {
    source("R/functions_erpm_bip_network.R", local = FALSE)
  }
}

# Wrapper ERPM (erpm + builder biparti)
if (!exists("erpm", mode = "function") || !exists("build_bipartite_from_inputs", mode = "function")) {
  if (file.exists("R/erpm.R")) {
    source("R/erpm.R", local = FALSE)
  } else {
    cat("[WARN] erpm()/build_bipartite_from_inputs indisponibles. Certaines phases seront sautées.\n")
  }
}

# ======================================================================================
# Réglages de run (clé du fichier)
# ======================================================================================
RUN <- list(
  phase1_expected = FALSE,
  phase2_equiv    = FALSE,
  phase3_fit      = TRUE,
  phase4_mcmc     = FALSE,

  quiet_phase1    = FALSE,
  quiet_phase2    = FALSE,
  quiet_phase3    = FALSE
)

# ======================================================================================
# Données de test
# ======================================================================================

# Panel de partitions variées (pour équivalence)
partitions <- list(
  P1 = c(1,2,2,3,3,3,4),                              # petite
  P2 = c(1,1, 2,2,2, 3,3,3,3, 4,4, 5,5,5, 6,6,6,6,6), # moyenne
  P3 = c(1,1,2,2,3)                                   # très simple
)

# Génère un jeu d'attributs contrôlé par partition
.make_nodes <- function(part) {
  n <- length(part)
  set.seed(123 + n)
  data.frame(
    label   = utils::head(LETTERS, n),
    age     = sample(20:60, n, replace = TRUE),                           # numérique
    gender  = sample(c("F", "H"), n, replace = TRUE),                     # binaire
    dept    = sample(c("Info", "RH", "Rech", "Vent"), n, replace = TRUE,  # catégoriel
                     prob = c(0.35, 0.25, 0.25, 0.15)),
    stringsAsFactors = FALSE
  )
}

# ======================================================================================
# Helpers pour construction du réseau biparti et summaries
# ======================================================================================

# Construire un biparti depuis build_bipartite_from_inputs (wrapper ERPM),
# avec fallback explicite vers partition_to_bipartite_network si besoin.
.erpm_build_bipartite_nw <- function(part, nodes) {
  stopifnot(length(part) == nrow(nodes))

  attrs <- as.list(nodes[, setdiff(names(nodes), "label"), drop = FALSE])

  # 1) via build_bipartite_from_inputs si dispo
  if (exists("build_bipartite_from_inputs", mode = "function")) {
    builder <- get("build_bipartite_from_inputs", mode = "function")

    candidates <- list(
      quote(builder(partition = part, nodes = nodes)),
      quote(builder(partition = part, labels = nodes$label, attributes = attrs)),
      quote(builder(partition = part, labels = nodes$label)),
      quote(builder(partition = part))
    )

    last_err <- NULL
    for (expr in candidates) {
      out <- try(eval(expr), silent = TRUE)
      if (inherits(out, "try-error") || is.null(out)) {
        last_err <- out
        next
      }
      if (inherits(out, "network")) return(out)
      if (is.list(out)) {
        for (nm in c("network", "nw", "net", "graph", "g", "bip")) {
          if (!is.null(out[[nm]]) && inherits(out[[nm]], "network")) return(out[[nm]])
        }
      }
    }

    # Si on arrive ici, le builder a échoué
    fm <- try(formals(builder), silent = TRUE)
    cat("[WARN] build_bipartite_from_inputs a échoué; fallback si possible.\n")
    if (!inherits(fm, "try-error")) cat("[WARN] formals(builder): ", paste(names(fm), collapse = ", "), "\n", sep = "")
    if (!is.null(last_err)) cat("[WARN] dernière erreur: ", as.character(last_err)[1], "\n", sep = "")
  }

  # 2) fallback partition_to_bipartite_network si dispo
  if (exists("partition_to_bipartite_network", mode = "function")) {
    return(partition_to_bipartite_network(
      labels     = nodes$label,
      partition  = part,
      attributes = attrs
    ))
  }

  stop("Aucun constructeur biparti valide n'a produit un objet 'network'.")
}

# Construire formule `nw ~ <rhs>` avec `nw` capturé dans l'env
.formula_nw <- function(nw, rhs_txt) {
  f <- as.formula(paste0("nw ~ ", rhs_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  f
}

# summary() sur réseau biparti explicite
summary_on_bipartite_network <- function(part, nodes, rhs_txt) {
  nw <- .erpm_build_bipartite_nw(part, nodes)
  f  <- .formula_nw(nw, rhs_txt)
  as.numeric(suppressMessages(summary(f)))
}

# summary() en passant par la traduction erpm() (dry-run -> formule ergm -> summary)
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

# Vérif stricte: summary(nw) == summary(ERPM-traduit)
check_summary_equivalence <- function(part, nodes, rhs_vec, tol = 0) {
  ok_all <- TRUE
  for (rhs in rhs_vec) {
    s_net  <- summary_on_bipartite_network(part, nodes, rhs)
    s_erpm <- summary_on_erpm_translation(part, nodes, rhs)

    cat(sprintf("[CHECK] n=%d RHS=%-45s net=%s erpm=%s\n",
                length(part), rhs, paste(s_net, collapse = ","), paste(s_erpm, collapse = ",")))

    if (any(!is.finite(s_net)))  stop("summary réseau non fini.")
    if (any(!is.finite(s_erpm))) stop("summary ERPM non fini.")
    if (length(s_net) != length(s_erpm)) stop("Longueur de vecteur de stats différente.")

    if (!all(abs(s_net - s_erpm) <= tol)) {
      ok_all <- FALSE
      cat(sprintf("  -> MISMATCH au-delà de tol=%g\n", tol))
    }
  }
  ok_all
}

# ======================================================================================
# Cas analytiques avec attentes explicites pour cov_fullmatch
# Définition attendue: compte le NOMBRE de groupes homogènes (pas la somme des tailles)
# ======================================================================================

# 1) Tie 2v2 vs groupe homogène 3
case_tie <- list(
  part = c(1,1,1,1,  2,2,2),                         # tailles: 4 et 3
  val  = c("A","A","B","B",  "Z","Z","Z"),
  checks = list(
    list(rhs = "cov_fullmatch('val')",               expect = 1), # seul groupe taille 3 est homogène
    list(rhs = "cov_fullmatch('val', size = 4)",     expect = 0), # filtre isole le 4 → non homogène
    list(rhs = "cov_fullmatch('val', category='Z')", expect = 1), # groupe 3 x 'Z'
    list(rhs = "cov_fullmatch('val', category='A')", expect = 0)  # aucun groupe tout-'A'
  )
)

# 2) Tailles 1,2,3 toutes homogènes
case_sizes <- list(
  part = c(1, 2,2, 3,3,3),
  val  = c("X","Y","Y","Z","Z","Z"),
  checks = list(
    list(rhs = "cov_fullmatch('val')",                  expect = 3), # 3 groupes homogènes
    list(rhs = "cov_fullmatch('val', size = 1)",        expect = 1), # seul singleton
    list(rhs = "cov_fullmatch('val', size = c(1,3))",   expect = 2), # groupes {1} et {3}
    list(rhs = "cov_fullmatch('val', category='Y')",    expect = 1), # le groupe {Y,Y}
    list(rhs = "cov_fullmatch('val', category='Z')",    expect = 1)  # le groupe {Z,Z,Z}
  )
)

# 3) Numérique dense: seuls les singletons comptent
set.seed(42)
case_dense <- list(
  part = c(1,1,1, 2,2,2, 3,3,3, 4, 5),               # 3 groupes de 3, 2 singletons
  val  = c(1,2,3, 4,5,6, 7,8,9, 10, 11),             # tous distincts
  checks = list(
    list(rhs = "cov_fullmatch('val')",           expect = 2), # les 2 singletons
    list(rhs = "cov_fullmatch('val', size = 3)", expect = 0)  # filtre supprime les singletons
  )
)

# --------------------------------------------------------------------------------------
# Cas volontairement en erreur — la logique métier exige une erreur.
# Laisse commenté par défaut. Décommente pour tester le fail-fast côté InitErgmTerm.
# --------------------------------------------------------------------------------------
# 4) Catégorie absente → DOIT PRODUIRE UNE ERREUR
# case_cat_absent <- list(
#   part = c(1,1, 2,2,2, 3,3,3,3),
#   val  = c("A","A", "B","B","B", "C","C","C","C"),
#   checks = list(
#     list(rhs="cov_fullmatch('val', category='ZZZ')",  expect=0)
#   )
# )
# 5) Valeurs NA dans la covariée → DOIT PRODUIRE UNE ERREUR
# case_na <- list(
#   part = c(1,1,1, 2,2,2, 3,3,3,3, 4),
#   val  = c(NA,NA,NA,  "A",NA,"A",  "B","B","B","B",  "C")
# )
# 6) size = integer(0) → DOIT PRODUIRE UNE ERREUR
# summary_on_bipartite_network(..., "cov_fullmatch('val', size = c())")

# ======================================================================================
# Contrôles ERGM pour fitting via erpm()
# ======================================================================================
ctrl_mle <- control.ergm(
  MCMLE.maxit      = 10,
  MCMC.samplesize  = 3000
)

# ======================================================================================
# PHASE 1: Summary — attentes numériques explicites
# ======================================================================================
run_phase1_summary_expected <- function(quiet = FALSE) {
  cat("\n=== PHASE 1 : SUMMARY (attentes explicites) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  run_one <- function(part, val, rhs, expect) {
    nodes <- data.frame(label = paste0("A", seq_along(part)), val = val, stringsAsFactors = FALSE)

    s_net <- summary_on_bipartite_network(part, nodes, rhs)
    if (!isTRUE(quiet)) {
      cat(sprintf("[EXPECT][nw]   RHS=%-40s obtenu=%s attendu=%s\n",
                  rhs, paste0(s_net, collapse = ","), expect))
    }

    if (length(s_net) != 1L || !is.finite(s_net)) stop("summary(nw) non scalaire ou non fini.")
    if (!isTRUE(all.equal(as.numeric(s_net), as.numeric(expect)))) {
      stop(sprintf("Mismatch summary(nw) RHS=%s : obtenu=%s attendu=%s",
                   rhs, as.numeric(s_net), as.numeric(expect)))
    }

    s_erpm <- summary_on_erpm_translation(part, nodes, rhs)
    if (!isTRUE(quiet)) {
      cat(sprintf("[EXPECT][erpm] RHS=%-40s obtenu=%s attendu=%s\n",
                  rhs, paste0(s_erpm, collapse = ","), expect))
    }

    if (is.na(s_erpm)) stop("summary(ERPM-traduit) a retourné NA (erpm indisponible ?).")
    if (!isTRUE(all.equal(as.numeric(s_erpm), as.numeric(expect)))) {
      stop(sprintf("Mismatch summary(ERPM) RHS=%s : obtenu=%s attendu=%s",
                   rhs, as.numeric(s_erpm), as.numeric(expect)))
    }

    invisible(TRUE)
  }

  for (ck in case_tie$checks)   run_one(case_tie$part,   case_tie$val,   ck$rhs, ck$expect)
  for (ck in case_sizes$checks) run_one(case_sizes$part, case_sizes$val, ck$rhs, ck$expect)
  for (ck in case_dense$checks) run_one(case_dense$part, case_dense$val, ck$rhs, ck$expect)

  cat("=== PHASE 1 OK ===\n")
  invisible(TRUE)
}

# ======================================================================================
# PHASE 2: Summary — équivalence réseau explicite vs ERPM-traduit
# ======================================================================================
cases_equiv <- c(
  "cov_fullmatch('gender')",
  "cov_fullmatch('gender', size = 2:4)",
  "cov_fullmatch('dept', category='Info')",
  "cov_fullmatch('dept', category='Rech', size = 2:5)"
)

run_phase2_summary_equiv <- function(quiet = FALSE) {
  cat("\n=== PHASE 2 : SUMMARY EQUIV (nw vs ERPM-traduit) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes(part)

    if (!isTRUE(quiet)) {
      cat(sprintf("\n--- Partition %s --- n=%d | groupes=%d | tailles: %s\n",
                  nm, length(part), length(unique(part)), paste(sort(table(part)), collapse = ",")))
    }

    ok <- check_summary_equivalence(part, nodes, cases_equiv, tol = 0)
    if (!ok) stop("Equivalence summary échouée.")
  }

  cat("=== PHASE 2 OK ===\n")
  invisible(TRUE)
}

# ======================================================================================
# PHASE 3: Fits courts via erpm() — coefficients finis
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

.group_diag_table <- function(part, nodes, col = "val") {
  stopifnot(col %in% names(nodes))
  gid  <- as.integer(part)
  vals <- nodes[[col]]
  split_idx <- split(seq_along(gid), gid)

  df <- do.call(rbind, lapply(names(split_idx), function(g) {
    idx <- split_idx[[g]]
    v   <- vals[idx]
    data.frame(
      group       = as.integer(g),
      size        = length(idx),
      n_unique    = length(unique(v)),
      homogeneous = as.integer(length(unique(v)) == 1L),
      values      = paste(v, collapse = ","),
      stringsAsFactors = FALSE
    )
  }))
  df[order(df$group), , drop = FALSE]
}

.print_fit_diagnostic <- function(tag, part, nodes, rhs_txt, col = "val") {
  cat("\n--- DIAGNOSTIC -------------------------------------------------\n")
  cat(sprintf("[TAG] %s\n", tag))
  cat(sprintf("[RHS] %s\n", rhs_txt))
  cat(sprintf("[Partition] n=%d | groupes=%d | tailles: %s\n",
              length(part), length(unique(part)), paste(sort(table(part)), collapse = ",")))

  s_obs <- tryCatch(summary_on_erpm_translation(part, nodes, rhs_txt),
                    error = function(e) { cat("[OBS] erreur summary(ERPM): ", conditionMessage(e), "\n", sep = ""); NA_real_ })
  cat("[Observed stat via summary(ERPM-traduit)] ", if (is.finite(s_obs)) s_obs else "NA", "\n", sep = "")

  diag_df <- .group_diag_table(part, nodes, col = col)
  cat("[Group diag] head:\n")
  print(utils::head(diag_df, 12L), row.names = FALSE)
  cat("---------------------------------------------------------------\n\n")
}

run_phase3_erpm_fits <- function(quiet = FALSE) {
  cat("\n=== PHASE 3 : FITS ERPM (coef finis) ===\n")
  if (isTRUE(quiet)) cat("  [mode quiet] sortie console réduite\n")

  if (!exists("erpm", mode = "function")) {
    cat("[PHASE 3] SKIP: erpm() indisponible\n")
    return(invisible(list(ok = NA)))
  }

  run_fit <- function(part, nodes, rhs, tag, col = "val") {
    f <- as.formula(paste0("partition ~ ", rhs))
    environment(f) <- list2env(list(partition = part, nodes = nodes), parent = parent.frame())

    res <- .with_warning_capture(
      try(erpm(f, eval.loglik = TRUE, verbose = FALSE, nodes = nodes, control = ctrl_mle), silent = TRUE)
    )

    fit <- res$value
    warns <- unique(res$warnings)

    if (length(warns) && !isTRUE(quiet)) {
      cat(sprintf("[ERPM-FIT %-10s] WARNINGS (%d)\n", tag, length(warns)))
      for (w in warns) cat("  - ", w, "\n", sep = "")
    }

    if (inherits(fit, "try-error")) {
      cat(sprintf("[ERPM-FIT %-10s] FAIL: %s\n", tag, as.character(fit)))
      .print_fit_diagnostic(tag, part, nodes, rhs, col = col)
      return(list(ok = FALSE, fit = NULL, coef = NA))
    }

    cf <- try(stats::coef(fit), silent = TRUE)
    ok_coef  <- !(inherits(cf, "try-error")) && length(cf) > 0L && all(is.finite(cf))
    ok_class <- inherits(fit, "ergm")
    ok <- isTRUE(ok_coef && ok_class)

    cat(sprintf("[ERPM-FIT %-10s] ok=%s | coef=%s\n",
                tag, ok,
                if (ok_coef) paste(format(as.numeric(cf)), collapse = ", ") else "NA"))

    if (!ok) .print_fit_diagnostic(tag, part, nodes, rhs, col = col)
    list(ok = ok, fit = fit, coef = if (ok_coef) cf else NA)
  }

  nodes_sizes <- data.frame(label = paste0("S", seq_along(case_sizes$part)), val = case_sizes$val, stringsAsFactors = FALSE)
  nodes_tie   <- data.frame(label = paste0("T", seq_along(case_tie$part)),   val = case_tie$val,   stringsAsFactors = FALSE)
  nodes_dense <- data.frame(label = paste0("D", seq_along(case_dense$part)), val = case_dense$val, stringsAsFactors = FALSE)

  fits <- list(
    ALL      = run_fit(case_sizes$part, nodes_sizes, "cov_fullmatch('val')",                "ALL"),
    S1       = run_fit(case_sizes$part, nodes_sizes, "cov_fullmatch('val', size = 1)",      "S1"),
    S1_3     = run_fit(case_sizes$part, nodes_sizes, "cov_fullmatch('val', size = c(1,3))", "S1_3"),
    TIE_ALL  = run_fit(case_tie$part,   nodes_tie,   "cov_fullmatch('val')",                "TIE_ALL"),
    TIE_S4   = run_fit(case_tie$part,   nodes_tie,   "cov_fullmatch('val', size = 4)",      "TIE_S4"),
    DENSE_ALL= run_fit(case_dense$part, nodes_dense, "cov_fullmatch('val')",                "DENSE_ALL"),
    DENSE_S3 = run_fit(case_dense$part, nodes_dense, "cov_fullmatch('val', size = 3)",      "DENSE_S3")
  )

  ok <- vapply(fits, function(x) isTRUE(x$ok), logical(1))
  n_ok <- sum(ok); n_tot <- length(ok)

  cat(sprintf("\n=== Bilan FITS: %d / %d OK ===\n", n_ok, n_tot))
  if (n_ok < n_tot) stop(sprintf("Echec FITS: %d KO", n_tot - n_ok))

  # En mode normal: résumés des fits
  if (!isTRUE(quiet)) {
    cat("\n=== Résumés fits OK ===\n")
    for (nm in names(fits)) {
      if (isTRUE(fits[[nm]]$ok)) {
        cat("\n---", nm, "---\n")
        print(summary(fits[[nm]]$fit))
      }
    }
  }

  invisible(fits)
}

# ======================================================================================
# PHASE 4: MCMC multi-toggle probe (diagnostic)
# ======================================================================================
# But:
#   - Déclencher le code MCMC sur un réseau biparti.
#   - Observer (si DEBUG côté C est activé) des traces indiquant un appel multi-toggle
#     (ntoggles > 1), comme vous l’avez fait pour squared_sizes.
#
# Important:
#   - Ce probe ne "force" pas mathématiquement ntoggles>1: il dépend du proposal.
#   - On choisit un proposal réputé pouvoir générer des multi-toggles selon les setups.
#   - Si vous n’observez jamais ntoggles>1, changez MCMC.prop vers un proposal multi-toggle
#     de votre codebase, ou augmentez l’intensité/longueur de la simulation.
# ======================================================================================

.run_mcmc_multitoggle_probe <- function(nw, rhs_txt) {
  ctrl <- control.simulate.formula(
    MCMC.burnin   = 1000,
    MCMC.interval = 1,
    MCMC.prop     = ~ sparse
  )

  f <- .formula_nw(nw, rhs_txt)

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
  cat("Objectif: déclencher le code MCMC et observer des traces debug multi-toggle côté C.\n")
  cat("Si vous ne voyez rien:\n")
  cat("  - soit le debug n'est pas activé dans le changestat C,\n")
  cat("  - soit le proposal ne génère pas de moves multi-toggle sur ce setup.\n\n")

  # Cas probe: on prend un cas simple avec covariée 'val'
  part  <- case_sizes$part
  nodes <- data.frame(label = paste0("P", seq_along(part)), val = case_sizes$val, stringsAsFactors = FALSE)
  nw <- .erpm_build_bipartite_nw(part, nodes)

  # RHS minimal: cov_fullmatch('val')
  rhs <- "cov_fullmatch('val')"

  cat("[PROBE] network.size=", network::network.size(nw), " | edges=", network::network.edgecount(nw), "\n", sep = "")
  cat("[PROBE] RHS=", rhs, "\n", sep = "")
  .run_mcmc_multitoggle_probe(nw, rhs)

  cat("\n[PROBE] Si des traces type 'MULTI-TOGGLE ntoggles=...' apparaissent: chemin multitoggle OK.\n")
  invisible(TRUE)
}

# ======================================================================================
# Run principal
# ======================================================================================
run_all_tests_cov_fullmatch <- function() {
  set.seed(1)

  cat("=== SELFTEST ERPM: cov_fullmatch ===\n")
  cat("R:", paste(R.version$major, R.version$minor, sep = "."), "\n")
  cat("ergm:", as.character(utils::packageVersion("ergm")), "\n")

  if (isTRUE(RUN$phase1_expected)) {
    run_phase1_summary_expected(quiet = isTRUE(RUN$quiet_phase1))
  } else {
    cat("\n=== PHASE 1 === SKIP (RUN$phase1_expected=FALSE)\n")
  }

  if (isTRUE(RUN$phase2_equiv)) {
    run_phase2_summary_equiv(quiet = isTRUE(RUN$quiet_phase2))
  } else {
    cat("\n=== PHASE 2 === SKIP (RUN$phase2_equiv=FALSE)\n")
  }

  if (isTRUE(RUN$phase3_fit)) {
    run_phase3_erpm_fits(quiet = isTRUE(RUN$quiet_phase3))
  } else {
    cat("\n=== PHASE 3 === SKIP (RUN$phase3_fit=FALSE)\n")
  }

  if (isTRUE(RUN$phase4_mcmc)) {
    run_phase4_mcmc_probe()
  } else {
    cat("\n=== PHASE 4 === SKIP (RUN$phase4_mcmc=FALSE)\n")
  }

  invisible(TRUE)
}

# Exécution quand lancé en script
if (identical(environment(), globalenv())) {
  run_all_tests_cov_fullmatch()
}

# --------------------------------------------------------------------------------------
# Fin: on désactive le patch si présent
# --------------------------------------------------------------------------------------
if (exists("ergm_patch_disable", mode = "function")) {
  ergm_patch_disable()
}

cat("\nTous les tests cov_fullmatch ont passé.\n")