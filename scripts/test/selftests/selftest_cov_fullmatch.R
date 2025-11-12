# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cov_fullmatch.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_fullmatch`
# Exécution: Rscript scripts/test/selftests/selftest_cov_fullmatch.R
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

# Patch ERGM optionnel (si utilisé dans le projet)
if (file.exists("scripts/ergm_patch.R")) {
    source("scripts/ergm_patch.R")
    ergm_patch_enable()
}

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
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
# Chargement utilitaires ERPM
# --------------------------------------------------------------------------------------
# Fallback si la fonction n'est pas exposée par le wrapper
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
# Données de test
# ======================================================================================

# Panel de partitions variées
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
    gender  = sample(c("F","H"), n, replace = TRUE),                      # binaire
    dept    = sample(c("Info","RH","Rech","Vent"), n, replace = TRUE,     # catégoriel
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
  attrs <- as.list(nodes[ , setdiff(names(nodes),"label"), drop=FALSE])
  builder <- NULL
  if (exists("build_bipartite_from_inputs", mode = "function")) builder <- get("build_bipartite_from_inputs")
  if (!is.null(builder)) {
    # Essaye plusieurs signatures courantes pour robustesse.
    out <- try({ builder(partition = part, nodes = nodes) }, silent = TRUE)
    if (inherits(out, "try-error") || is.null(out)) {
      out <- try({ builder(partition = part, labels = nodes$label, attributes = attrs) }, silent = TRUE)
    }
    if (!inherits(out, "try-error") && !is.null(out)) return(out)
  }
  if (exists("partition_to_bipartite_network", mode = "function")) {
    return(partition_to_bipartite_network(labels = nodes$label, partition = part, attributes = attrs))
  }
  stop("Aucun constructeur biparti disponible.")
}

# Redéfinition robuste du constructeur biparti pour garantir un objet 'network'
.erpm_build_bipartite_nw <- function(part, nodes) {
  stopifnot(length(part) == nrow(nodes))
  attrs <- as.list(nodes[ , setdiff(names(nodes),"label"), drop=FALSE])
  # 1) Tentative via build_bipartite_from_inputs
  if (exists("build_bipartite_from_inputs", mode = "function")) {
    builder <- get("build_bipartite_from_inputs")
    # Signature 1: (partition, nodes)
    out <- try(builder(partition = part, nodes = nodes), silent = TRUE)
    # Signature 2: (partition, labels, attributes)
    if (inherits(out, "try-error") || is.null(out)) {
      out <- try(builder(partition = part, labels = nodes$label, attributes = attrs), silent = TRUE)
    }
    # Extraction robuste du network si l'API retourne une liste
    if (!inherits(out, "try-error") && !is.null(out)) {
      if (inherits(out, "network")) return(out)
      if (is.list(out)) {
        cand_names <- c("network","nw","g","graph","bip","net")
        for (nm in cand_names) {
          if (!is.null(out[[nm]]) && inherits(out[[nm]], "network")) return(out[[nm]])
        }
      }
    }
  }
  # 2) Fallback explicite
  if (exists("partition_to_bipartite_network", mode = "function")) {
    return(partition_to_bipartite_network(labels = nodes$label, partition = part, attributes = attrs))
  }
  stop("Aucun constructeur biparti valide n'a produit un objet 'network'.")
}

# Construire formule `nw ~ <rhs>` avec `nw` capturé dans l'env
.formula_nw <- function(nw, rhs_txt) {
  f <- as.formula(paste0("nw ~ ", rhs_txt))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  f
}

# Exécuter un summary côté réseau biparti explicite
summary_on_bipartite_network <- function(part, nodes, rhs_txt) {
  nw <- .erpm_build_bipartite_nw(part, nodes)
  f  <- .formula_nw(nw, rhs_txt)
  as.numeric(suppressMessages(summary(f)))
}

# Exécuter un summary côté ERPM (LHS = partition) via wrapper, en réutilisant le RHS traduit
summary_on_erpm_translation <- function(part, nodes, rhs_txt) {
  if (!exists("erpm", mode = "function")) return(NA_real_)
  partition <- part
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition, nodes = nodes), parent = parent.frame())

  call_ergm <- erpm(f, eval_call = FALSE, verbose = FALSE, nodes = nodes)
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

# Vérif stricte: summary(nw) == summary(ERPM-traduit) sur plusieurs RHS
check_summary_equivalence <- function(part, nodes, rhs_vec, tol=0) {
  ok_all <- TRUE
  for (rhs in rhs_vec) {
    s_net  <- summary_on_bipartite_network(part, nodes, rhs)
    s_erpm <- summary_on_erpm_translation(part, nodes, rhs)
    cat(sprintf("[CHECK] part(n=%d) RHS=%-50s  net=%s  erpm=%s\n",
                length(part), rhs, paste(s_net, collapse=","), paste(s_erpm, collapse=",")))
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
# Définition: compte le NOMBRE de groupes homogènes (et non la somme des tailles)
# ======================================================================================

# 1) Tie 2v2 vs groupe homogène 3
case_tie <- list(
  part = c(1,1,1,1,  2,2,2),                         # tailles: 4 et 3
  val  = c("A","A","B","B",  "Z","Z","Z"),
  checks = list(
    list(rhs="cov_fullmatch('val')",                  expect=1),         # seul groupe taille 3 est homogène
    list(rhs="cov_fullmatch('val', size = 4)",        expect=0),         # filtre isole le 4 → non homogène
    list(rhs="cov_fullmatch('val', category='Z')",    expect=1),         # groupe 3 x 'Z'
    list(rhs="cov_fullmatch('val', category='A')",    expect=0)          # aucun groupe tout-'A'
  )
)

# 2) Tailles 1,2,3 toutes homogènes
case_sizes <- list(
  part = c(1, 2,2, 3,3,3),
  val  = c("X","Y","Y","Z","Z","Z"),
  checks = list(
    list(rhs="cov_fullmatch('val')",                  expect=3),         # 3 groupes homogènes
    list(rhs="cov_fullmatch('val', size = 1)",        expect=1),         # seul singleton
    list(rhs="cov_fullmatch('val', size = c(1,3))",   expect=2),         # groupes {1} et {3}
    list(rhs="cov_fullmatch('val', category='Y')",    expect=1),         # le groupe {Y,Y}
    list(rhs="cov_fullmatch('val', category='Z')",    expect=1)          # le groupe {Z,Z,Z}
  )
)

# 3) Numérique dense: seuls les singletons comptent
set.seed(42)
case_dense <- list(
  part = c(1,1,1, 2,2,2, 3,3,3, 4, 5),               # 3 groupes de 3, 2 singletons
  val  = c(1,2,3, 4,5,6, 7,8,9,  10, 11),            # tous distincts
  checks = list(
    list(rhs="cov_fullmatch('val')",                  expect=2),         # les 2 singletons
    list(rhs="cov_fullmatch('val', size = 3)",        expect=0)          # filtre supprime les singletons
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
#     list(rhs="cov_fullmatch('val', category='ZZZ')",  expect=0) # intentionnellement invalide
#   )
# )

# 5) Valeurs NA dans la covariée → DOIT PRODUIRE UNE ERREUR
# case_na <- list(
#   part = c(1,1,1, 2,2,2, 3,3,3,3, 4),
#   val  = c(NA,NA,NA,  "A",NA,"A",  "B","B","B","B",  "C")
# )

# 6) size = integer(0) → DOIT PRODUIRE UNE ERREUR
# Exemple d’appel à décommenter dans run_phase1_expected():
# summary_on_bipartite_network(case_sizes$part, data.frame(label=..., val=...), "cov_fullmatch('val', size = c())")

# ======================================================================================
# Contrôles ERGM pour fitting via erpm()
# ======================================================================================

ctrl_mle <- control.ergm(
  MCMLE.maxit      = 10,
  MCMC.samplesize = 3000
  # MCMC.burnin      = 5000,
  # MCMC.interval    = 1000,
  # force.main       = TRUE,
  # parallel         = 0,
  # seed             = 1
)

# ======================================================================================
# Phase 1: Summary — attentes numériques explicites
# ======================================================================================

run_phase1_summary_expected <- function() {
  cat("=== PHASE 1 : Summary avec attentes explicites ===\n")

  run_one <- function(part, val, rhs, expect) {
    nodes <- data.frame(label = paste0("A", seq_along(part)), val = val, stringsAsFactors = FALSE)
    s <- summary_on_bipartite_network(part, nodes, rhs)
    cat(sprintf("[EXPECT] RHS=%-45s  obtenu=%s  attendu=%s\n", rhs, paste0(s, collapse=","), expect))
    if (length(s)!=1L || !is.finite(s)) stop("summary non scalaire ou non fini.")
    if (!isTRUE(all.equal(as.numeric(s), as.numeric(expect)))) {
      stop(sprintf("Mismatch sur RHS=%s : obtenu=%s attendu=%s", rhs, as.numeric(s), as.numeric(expect)))
    }
    s2 <- summary_on_erpm_translation(part, nodes, rhs)
    if (!isTRUE(all.equal(as.numeric(s2), as.numeric(expect)))) {
      stop(sprintf("Mismatch ERPM sur RHS=%s : obtenu=%s attendu=%s", rhs, as.numeric(s2), as.numeric(expect)))
    }
  }

  for (ck in case_tie$checks)   run_one(case_tie$part,   case_tie$val,   ck$rhs, ck$expect)
  for (ck in case_sizes$checks) run_one(case_sizes$part, case_sizes$val, ck$rhs, ck$expect)
  for (ck in case_dense$checks) run_one(case_dense$part, case_dense$val, ck$rhs, ck$expect)

  # -------------------------------
  # BLOCS D'ERREURS À TESTER À LA DEMANDE
  # Décommentez un bloc pour valider le fail-fast correspondant.
  # -------------------------------

  # 1) Category absente -> erreur attendue
  # err <- NULL
  # tryCatch({
  #   summary_on_bipartite_network(case_sizes$part,
  #                       data.frame(label=paste0("S",seq_along(case_sizes$part)),
  #                                  val=case_sizes$val, stringsAsFactors=FALSE),
  #                       "cov_fullmatch('val', category='ZZZ')")
  # }, error = function(e) err <<- e$message)
  # if (is.null(err) || !grepl("cov_fullmatch: 'category' absente des modalités", err))
  #   stop("Erreur attendue sur category absente non levée côté summary(nw).")

  # 2) NA dans la covariée -> erreur attendue
  # part_na <- c(1,1,1, 2,2,2, 3,3,3,3, 4)
  # val_na  <- c(NA,NA,NA, "A",NA,"A",  "B","B","B","B",  "C")
  # err <- NULL
  # tryCatch({
  #   summary_on_bipartite_network(part_na,
  #                       data.frame(label=paste0("N",seq_along(part_na)),
  #                                  val=val_na, stringsAsFactors=FALSE),
  #                       "cov_fullmatch('val')")
  # }, error = function(e) err <<- e$message)
  # if (is.null(err) || !grepl("cov_fullmatch: NA non autorisé", err))
  #   stop("Erreur attendue sur NA non levée côté summary(nw).")

  # 3) size = integer(0) -> erreur attendue
  # err <- NULL
  # tryCatch({
  #   summary_on_bipartite_network(case_sizes$part,
  #                       data.frame(label=paste0("S",seq_along(case_sizes$part)),
  #                                  val=case_sizes$val, stringsAsFactors=FALSE),
  #                       "cov_fullmatch('val', size = c())")
  # }, error = function(e) err <<- e$message)
  # if (is.null(err) || !grepl("cov_fullmatch: 'size' vide (integer\\(0\\)) interdit", err))
  #   stop("Erreur attendue sur size=integer(0) non levée côté summary(nw).")

  cat("=== Phase 1 OK ===\n")
  invisible(NULL)
}

# ======================================================================================
# Phase 2: Summary — équivalence réseau explicite vs ERPM-traduit
# ======================================================================================

cases_equiv <- c(
  "cov_fullmatch('gender')",
  "cov_fullmatch('gender', size = 2:4)",
  "cov_fullmatch('dept', category='Info')",
  "cov_fullmatch('dept', category='Rech', size = 2:5)"
)

run_phase2_summary_equiv <- function() {
  cat("\n=== PHASE 2 : Summary(nw) vs Summary(ERPM-traduit) ===\n")
  for (nm in names(partitions)) {
    part  <- partitions[[nm]]
    nodes <- .make_nodes(part)
    cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))
    ok <- check_summary_equivalence(part, nodes, cases_equiv, tol = 0)
    if (!ok) stop("Equivalence summary échouée.")
  }
  cat("=== Phase 2 OK ===\n")
  invisible(NULL)
}

# ======================================================================================
# Phase 3: Fits courts via erpm() — exécution et coefficients finis
# ======================================================================================

# ======================================================================================
# Helpers diagnostics pour Phase 3
# ======================================================================================

# Capture et mémorise les warnings pendant une évaluation
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

# Tableau diagnostique par groupe: taille et homogénéité de la covariée 'val'
.group_diag_table <- function(part, nodes) {
  stopifnot("val" %in% names(nodes))
  gid <- as.integer(part)
  vals <- nodes$val
  split_idx <- split(seq_along(gid), gid)
  df <- do.call(rbind, lapply(names(split_idx), function(g) {
    idx <- split_idx[[g]]
    v   <- vals[idx]
    data.frame(
      group      = as.integer(g),
      size       = length(idx),
      n_unique   = length(unique(v)),
      homogeneous= as.integer(length(unique(v)) == 1L),
      values     = paste(v, collapse = ","),
      stringsAsFactors = FALSE
    )
  }))
  df[order(df$group), , drop=FALSE]
}

# Impression compacte d’un diagnostic complet
.print_fit_diagnostic <- function(tag, part, nodes, rhs_txt) {
  cat("\n--- DIAGNOSTIC -------------------------------------------------\n")
  cat(sprintf("[TAG] %s\n", tag))
  cat(sprintf("[RHS] %s\n", rhs_txt))
  cat(sprintf("[Partition] n=%d | groupes=%d | tailles: %s\n",
              length(part), length(unique(part)),
              paste(sort(table(part)), collapse=",")))

  # Stat observée côté summary(ERPM-traduit)
  s_obs <- NA_real_
  err_s <- NULL
  tryCatch({
    s_obs <- summary_on_erpm_translation(part, nodes, rhs_txt)
  }, error = function(e) err_s <<- conditionMessage(e))
  cat("[Observed stat via summary(ERPM-traduit)] ",
      if (is.finite(s_obs)) as.character(s_obs) else paste0("NA", if (!is.null(err_s)) paste0(" (", err_s, ")") else ""),
      "\n", sep="")

  # RHS/constraints issus du dry-run erpm()
  partition <- part
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition, nodes = nodes), parent = parent.frame())
  call_ergm <- try(erpm(f, eval_call = FALSE, verbose = FALSE, nodes = nodes), silent = TRUE)
  if (!inherits(call_ergm, "try-error")) {
    ergm_form <- call_ergm[[2L]]
    rhs_expr  <- try(ergm_form[[3L]], silent = TRUE)
    call_args <- as.list(call_ergm)[-1L]
    cons      <- call_args$constraints
    cat("[Dry-run] constraints = ",
        if (!is.null(cons)) deparse(cons) else "~ b1part",
        "\n", sep="")
    cat("[Dry-run] RHS expr    = ",
        if (!inherits(rhs_expr, "try-error")) paste(deparse(rhs_expr), collapse=" ") else "?", "\n", sep="")
  } else {
    cat("[Dry-run] Impossible d’obtenir l’appel erpm(): ",
        as.character(call_ergm), "\n", sep="")
  }

  # Aperçu homogénéité par groupe
  diag_df <- .group_diag_table(part, nodes)
  cat("[Group diag] head:\n")
  print(utils::head(diag_df, 10L), row.names = FALSE)
  cat("---------------------------------------------------------------\n\n")
}

# ======================================================================================
# Phase 3: Fits via erpm() — version MLE + loglik
# ======================================================================================

run_phase3_erpm_fits <- function() {
  cat("\n=== PHASE 3 : Fits erpm() (MLE + loglik) ===\n")

  run_fit <- function(part, nodes, rhs, tag) {
    if (!exists("erpm", mode = "function")) {
      cat(sprintf("[ERPM-FIT %-14s] SKIP (erpm() indisponible)\n", tag));
      return(list(ok = NA, error = TRUE, coef = NA, fit = NULL))
    }
    f <- as.formula(paste0("partition ~ ", rhs))
    environment(f) <- list2env(list(partition = part, nodes = nodes), parent = parent.frame())

    # Stat observée avant fit
    s_obs <- NA_real_
    err_obs <- NULL
    s_obs <- tryCatch(summary_on_erpm_translation(part, nodes, rhs),
                      error = function(e){ err_obs <<- conditionMessage(e); NA_real_ })
    if (!is.na(s_obs)) {
      cat(sprintf("[ERPM-FIT %-14s] stat_observee=%s\n", tag, format(s_obs)))
    } else if (!is.null(err_obs)) {
      cat(sprintf("[ERPM-FIT %-14s] stat_observee=NA (%s)\n", tag, err_obs))
    }

    # Exécution avec capture des warnings
    res <- .with_warning_capture(
      try(erpm(f, estimate = "MLE", eval.loglik = TRUE, control = ctrl_mle,
               verbose = FALSE, nodes = nodes), silent = TRUE)
    )
    fit <- res$value
    warns <- res$warnings
    if (length(warns)) {
      cat(sprintf("[ERPM-FIT %-14s] WARNINGS (%d):\n", tag, length(warns)))
      for (w in unique(warns)) cat("  - ", w, "\n", sep="")
    }

    # Cas erreur d’exécution
    if (inherits(fit,"try-error")) {
      cat(sprintf("[ERPM-FIT %-14s] ERREUR: %s\n", tag, as.character(fit)))
      .print_fit_diagnostic(tag, part, nodes, rhs)
      return(list(ok = FALSE, error = TRUE, coef = NA, fit = NULL))
    }

    # Vérification coefficients
    cf <- try(stats::coef(fit), silent = TRUE)
    ok_coef <- !(inherits(cf, "try-error")) && all(is.finite(cf))
    ok_class <- inherits(fit, "ergm")
    cat(sprintf("[ERPM-FIT %-14s] coef finies: %s | coef=%s\n",
                tag, if (ok_coef) "OK" else "KO",
                if (ok_coef) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

    if (!(ok_coef && ok_class)) {
      sm <- try(suppressMessages(summary(fit)), silent = TRUE)
      if (!inherits(sm, "try-error")) {
        cat(sprintf("[ERPM-FIT %-14s] summary(fit) extract:\n", tag))
        if (!is.null(sm$coefs)) {
          print(sm$coefs)
        } else if (!is.null(sm$coefficients)) {
          print(sm$coefficients)
        } else {
          str(sm, max.level = 1)
        }
      }
      .print_fit_diagnostic(tag, part, nodes, rhs)

      # Tolérance cas extrême attendu
      if (is.finite(s_obs) && s_obs == 0) {
        cat(sprintf("[ERPM-FIT %-14s] CAS EXTRÊME ATTENDU: stat_observee=0 -> PASS technique\n", tag))
        return(list(ok = TRUE, error = FALSE, coef = if (ok_coef) cf else NA, fit = fit))
      }
      return(list(ok = FALSE, error = FALSE, coef = if (ok_coef) cf else NA, fit = fit))
    }

    list(ok = TRUE, error = FALSE, coef = cf, fit = fit)
  }

  fit_results <- list()

  # sizes: 1,2,3
  fit_results[["ALL"]]   <- run_fit(case_sizes$part,
                            data.frame(label=paste0("S",seq_along(case_sizes$part)),
                                       val=case_sizes$val, stringsAsFactors=FALSE),
                            "cov_fullmatch('val')", "ALL")
  fit_results[["S1"]]    <- run_fit(case_sizes$part,
                            data.frame(label=paste0("S",seq_along(case_sizes$part)),
                                       val=case_sizes$val, stringsAsFactors=FALSE),
                            "cov_fullmatch('val', size = c(1))", "S1")
  fit_results[["S1_3"]]  <- run_fit(case_sizes$part,
                            data.frame(label=paste0("S",seq_along(case_sizes$part)),
                                       val=case_sizes$val, stringsAsFactors=FALSE),
                            "cov_fullmatch('val', size = c(1,3))", "S1_3")

  # tie 2v2 et filtre size=4
  fit_results[["TIE_ALL"]] <- run_fit(case_tie$part,
                            data.frame(label=paste0("T",seq_along(case_tie$part)),
                                       val=case_tie$val, stringsAsFactors=FALSE),
                            "cov_fullmatch('val')", "TIE_ALL")
  fit_results[["TIE_S4"]]  <- run_fit(case_tie$part,
                            data.frame(label=paste0("T",seq_along(case_tie$part)),
                                       val=case_tie$val, stringsAsFactors=FALSE),
                            "cov_fullmatch('val', size = c(4))", "TIE_S4")

  # numérique dense
  fit_results[["DENSE_ALL"]] <- run_fit(case_dense$part,
                            data.frame(label=paste0("D",seq_along(case_dense$part)),
                                       val=case_dense$val, stringsAsFactors=FALSE),
                            "cov_fullmatch('val')", "DENSE_ALL")
  fit_results[["DENSE_S3"]]  <- run_fit(case_dense$part,
                            data.frame(label=paste0("D",seq_along(case_dense$part)),
                                       val=case_dense$val, stringsAsFactors=FALSE),
                            "cov_fullmatch('val', size = c(3))", "DENSE_S3")

  # Bilan global + résumés détaillés
  ok <- vapply(fit_results, function(x) isTRUE(x$ok), logical(1))
  n_ok <- sum(ok, na.rm = TRUE); n_tot <- sum(!is.na(ok))
  cat(sprintf("\n=== Bilan fits erpm() : %d / %d OK ===\n", n_ok, n_tot))

  cat("\n=== Résumés détaillés des fits ERPM réussis ===\n")
  for (nm in names(fit_results)) {
    fit_obj <- fit_results[[nm]]
    if (isTRUE(fit_obj$ok) && inherits(fit_obj$coef, "numeric")) {
      cat(sprintf("\n--- Résumé fit %s ---\n", nm))
      print(summary(fit_obj$fit))
    }
  }

  if (n_ok < n_tot) stop(sprintf("Echec fits: %d KO", n_tot - n_ok))

  invisible(list(fit_results = fit_results))
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: cov_fullmatch ===\n")
run_phase1_summary_expected()
run_phase2_summary_equiv()

res_fits <- run_phase3_erpm_fits()

ergm_patch_disable()
cat("\nTous les tests cov_fullmatch ont passé.\n")
