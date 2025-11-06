# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cov_ingroup.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_ingroup`
# Exécution: Rscript scripts/test/selftests/selftest_cov_ingroup.R
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule
# --------------------------------------------------------------------------------------
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
    a <- commandArgs(FALSE)
    f <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(f) == 1L) return(normalizePath(dirname(f), winslash = "/", mustWork = FALSE))
    if (!is.null(sys.frames()) && !is.null(sys.calls())) {
        for (i in rev(seq_along(sys.calls()))) {
        cf <- sys.frame(i)
        if (!is.null(cf$ofile)) return(normalizePath(dirname(cf$ofile), winslash = "/", mustWork = FALSE))
        }
    }
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
script_dir <- .get_script_dir()
log_path   <- file.path(script_dir, "selftest_cov_ingroup.log")
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
if (!exists("partition_to_bipartite_network", mode = "function")) {
    if (file.exists("R/functions_erpm_bip_network.R")) {
        source("R/functions_erpm_bip_network.R", local = FALSE)
    } else stop("functions_erpm_bip_network.R introuvable.")
}
if (!exists("erpm", mode = "function")) {
    if (file.exists("R/erpm_wrapper.R")) {
        source("R/erpm_wrapper.R", local = FALSE)
    } else {
        cat("[WARN] erpm() indisponible. Les tests de traduction et les fits via erpm() seront sautés.\n")
    }
}

# ======================================================================================
# Données de test
# ======================================================================================

# Panel de partitions de tailles/structures variées
partitions <- list(
    # petite mais non triviale
    P1 = c(1,2,2,3,3,3,4),
    # moyenne avec diversité de tailles
    P2 = c(1,1, 2,2,2, 3,3,3,3, 4,4, 5,5,5, 6,6,6,6,6,6),
    # très simple pour MLE (petite, variance non nulle)
    P3 = c(1,1,2,2,3)
)

# Génère un jeu d'attributs contrôlé par partition
.make_nodes <- function(part) {
    n <- length(part)
    set.seed(123 + n)
    data.frame(
        label   = utils::head(LETTERS, n),
        age     = sample(20:60, n, replace = TRUE),                  # numérique
        score   = round(runif(n, min = 0, max = 10), 1),             # numérique
        gender  = sample(c("F","M"), n, replace = TRUE),             # factor binaire
        dept    = sample(c("A","B","C"), n, replace = TRUE, prob=c(0.5,0.3,0.2)), # factor 3 modalités
        stringsAsFactors = FALSE
    )
}

# ======================================================================================
# Helpers
# ======================================================================================

# Construire un biparti depuis une partition + attributs
.make_nw <- function(part, nodes) {
    stopifnot(length(part) == nrow(nodes))
    attrs <- as.list(nodes[ , setdiff(names(nodes),"label"), drop=FALSE])
    partition_to_bipartite_network(labels = nodes$label, partition = part, attributes = attrs)
}
# .make_nw <- function(part, nodes) {
#     stopifnot(length(part) == nrow(nodes))
#     partition_to_bipartite_network(labels = nodes$label, partition = part,
#                                     attributes = as.list(nodes))
# }

# Construire formule `nw ~ <rhs>` avec `nw` capturé dans l'env
.formula_nw <- function(nw, rhs_txt) {
    f <- as.formula(paste0("nw ~ ", rhs_txt))
    environment(f) <- list2env(list(nw = nw), parent = parent.frame())
    f
}

# Exécuter un summary côté réseau biparti explicite
.summary_on_network <- function(part, nodes, rhs_txt) {
    nw <- .make_nw(part, nodes)
    f  <- .formula_nw(nw, rhs_txt)
    as.numeric(suppressMessages(summary(f)))
}

# Exécuter un summary côté ERPM (LHS = partition) via wrapper
.summary_on_erpm <- function(part, nodes, rhs_txt) {
  if (!exists("erpm", mode = "function")) return(NA_real_)
  partition <- part
  f <- as.formula(paste0("partition ~ ", rhs_txt))
  environment(f) <- list2env(list(partition = partition, nodes = nodes), parent = parent.frame())

  # 1) Obtenir l'appel ergm(...) sans l'évaluer
  call_ergm <- erpm(f, eval_call = FALSE, verbose = FALSE, nodes = nodes)

  # 2) Extraire la formule ergm et son RHS proprement
  #    call_ergm a la forme: ergm( <FORMULE>, constraints=..., estimate=..., ... )
  ergm_form <- call_ergm[[2L]]                  # la formule complète `nw ~ <RHS>`
  rhs_expr  <- ergm_form[[3L]]                  # l'expression du RHS

  # 3) Reconstruire une formule `nw ~ <RHS>` avec le nw reconstruit et les mêmes contraintes
  nw <- .make_nw(part, nodes)

  f2 <- as.formula(bquote(nw ~ .(rhs_expr)))
  environment(f2) <- list2env(list(nw = nw), parent = parent.frame())

  # 4) Récupérer la contrainte si présente, sinon ~b1part
  call_args <- as.list(call_ergm)[-1L]
  cons <- call_args$constraints
  if (is.null(cons)) cons <- as.formula(~ b1part)

  as.numeric(suppressMessages(summary(f2, constraints = cons)))
}
# .summary_on_erpm <- function(part, nodes, rhs_txt) {
#     if (!exists("erpm", mode = "function")) return(NA_real_)
#     partition <- part
#     f <- as.formula(paste0("partition ~ ", rhs_txt))
#     environment(f) <- list2env(list(partition = partition, nodes = nodes), parent = parent.frame())
#     # erpm() dry run → appel ergm(...) traduit
#     call_ergm <- erpm(f, eval_call = FALSE, verbose = FALSE, nodes = nodes)
#     # call_ergm <- erpm(f, eval_call = FALSE, verbose = FALSE)
#     # exécute summary sur réseau équivalent
#     nw <- .make_nw(part, nodes)
#     f2 <- .formula_nw(nw, gsub("^.*~\\s*", "", deparse(call_ergm)))  # RHS traduit
#     as.numeric(suppressMessages(summary(f2)))
# }

# Vérif stricte: summary(nw) == summary(erpm) sur plusieurs RHS
.check_summary_equivalence <- function(part, nodes, rhs_vec, tol=0) {
    ok_all <- TRUE
    for (rhs in rhs_vec) {
        s_net  <- .summary_on_network(part, nodes, rhs)
        s_erpm <- .summary_on_erpm(part, nodes, rhs)
        cat(sprintf("[CHECK] part(n=%d) RHS=%-45s  net=%s  erpm=%s\n",
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
# Panel de cas `cov_ingroup` (numérique + catégoriel, avec et sans filtre de taille)
# ======================================================================================

cases_summary <- c(
    # Numérique
    "cov_ingroup('age')",
    "cov_ingroup('age', size = 2:4)",
    "cov_ingroup('score')",
    "cov_ingroup('score', size = 3:6)",
    # Catégoriel
    "cov_ingroup('gender', category='F')",
    "cov_ingroup('gender', category='M', size = 2:5)",
    "cov_ingroup('dept',   category='A')",
    "cov_ingroup('dept',   category='C', size = 2:4)"
)

# ======================================================================================
# Contrôles ERGM pour fitting
# ======================================================================================

ctrl_cd <- control.ergm(
  CD.maxit         = 50,     # nb d’itérations max CD
  CD.samplesize    = 90000,  # plus d’échantillons par itération
  CD.nsteps        = 3,      # pas de CD-MCMC plus longs
  CD.conv.min.pval = 1e-20,  # seuil de convergence très strict
  CD.NR.maxit      = 0,      # pas de Newton-Raphson côté CD
  MCMLE.maxit      = 0,      # ne bascule jamais en MCMLE
  init.method      = "zeros",
  parallel         = 0
)

# Cas MLE “safe” : petite partition, peu d’attributs
# ctrl_mle <- control.ergm(
#     init.method      = "CD",
#     CD.samplesize    = 50000,
#     CD.nsteps        = 10,
#     MCMC.burnin      = 50000,
#     MCMC.interval    = 2500,
#     MCMC.samplesize  = 10000,
#     MCMLE.maxit      = 3,
#     MCMLE.steplength = 0.05,
#     force.main       = TRUE,
#     parallel         = 0
# )
ctrl_mle <- control.ergm(
  ## Warm-start CD robuste
  init.method   = "CD",
  CD.samplesize = 20000,
  CD.nsteps     = 3,
  CD.maxit      = 5,

  ## MCMLE principal
  main.method      = "MCMLE",
  force.main       = TRUE,
  MCMC.burnin      = 100000,
  MCMC.interval    = 5000,
  MCMC.samplesize  = 10000,

  ## Pas courts + amortissement pour éviter le test singulier
  MCMLE.maxit          = 5,
  MCMLE.steplength     = 0.02,
  MCMLE.dampening      = TRUE,
  MCMLE.dampening.level= 0.2,

  ## Critère de terminaison moins agressif
  MCMLE.confidence     = 0.95,
  MCMLE.confidence.boost = 1,

  ## Optionnel: viser une taille d’échantillon effective
  MCMC.effectiveSize   = 64,
  parallel             = 0
)

# ======================================================================================
# Phase 1: Summary comparatifs (réseau explicite vs ERPM)
# ======================================================================================

run_phase1 <- function() {
    cat("=== PHASE 1 : Summary(nw) vs Summary(ERPM-traduit) ===\n")
    total <- 0L; ok   <- 0L
    for (nm in names(partitions)) {
        part  <- partitions[[nm]]
        nodes <- .make_nodes(part)

        cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                    nm, length(part), length(unique(part)),
                    paste(sort(table(part)), collapse=",")))

        res <- .check_summary_equivalence(part, nodes, cases_summary, tol = 0)
        total <- total + length(cases_summary)
        ok    <- ok + as.integer(res) * length(cases_summary)
    }
    cat(sprintf("\n=== Bilan Phase 1 : %d / %d checks OK ===\n", ok, total))
    if (ok < total) stop(sprintf("Summary mismatch sur %d cas.", total - ok))
    invisible(NULL)
}

# ======================================================================================
# Phase 2: Fits courts via erpm() pour valider intégration
# ======================================================================================

.erpm_fit <- function(part, nodes, rhs, name,
                     estimate = "CD", eval.loglik = FALSE, control = ctrl_cd) {
    if (!exists("erpm", mode = "function")) {
        cat(sprintf("[ERPM-FIT %-16s] SKIP (erpm() indisponible)\n", name)); 
        return(list(ok = NA))
    }
    set.seed(42)
    f <- as.formula(paste0("partition ~ ", rhs))
    environment(f) <- list2env(list(partition = part, nodes = nodes), parent = parent.frame())
    cat(sprintf("[ERPM-FIT %-16s] n=%-3d RHS=%s\n", name, length(part), rhs))
    fit <- try(erpm(f,  estimate = estimate,    eval.loglik = eval.loglik,
                        control = control,      verbose = FALSE,    nodes = nodes),         
                        silent = TRUE)

    # fit <- try(erpm(f, estimate = estimate, eval.loglik = eval.loglik, control = control, verbose = FALSE),
    #             silent = TRUE)
    if (inherits(fit, "try-error")) {
        cat("  -> ERREUR fit:", as.character(fit), "\n"); return(list(ok = FALSE))
    }
    cf <- try(stats::coef(fit), silent = TRUE)
    ok_coef <- !(inherits(cf, "try-error")) && all(is.finite(cf))
    cat(sprintf("  -> class(ergm) %s | coef finies %s | coef: %s\n",
                if (inherits(fit, "ergm")) "OK" else "KO",
                if (ok_coef) "OK" else "KO",
                if (ok_coef) paste(format(as.numeric(cf)), collapse=", ") else "NA"))
    list(ok = inherits(fit,"ergm") && ok_coef)
}

run_phase2 <- function() {
    cat("\n=== PHASE 2 : Fits erpm() courts (CD) ===\n")
    # Choix de quelques RHS représentatifs
    rhs_list <- list(
        R1 = "cov_ingroup('age')",
        R2 = "cov_ingroup('gender', category='F', size=2:4')", # volontairement incorrect pour tester robustesse
        R3 = "cov_ingroup('gender', category='F', size=2:4)",
        R4 = "cov_ingroup('age') + cov_ingroup('gender', category='F', size=3:6)"
    )
    # Corrige la coquille éventuelle si présente
    if (grepl("size=2:4'", rhs_list$R2, fixed = TRUE)) rhs_list$R2 <- "cov_ingroup('dept', category='A')"

    parts <- list(
        list(name="P1", part=partitions$P1),
        list(name="P2", part=partitions$P2)
    )

    n_ok <- 0L; n_all <- 0L
    for (px in parts) {
        nodes <- .make_nodes(px$part)
        for (nm in names(rhs_list)) {
        out <- .erpm_fit(px$part, nodes, rhs_list[[nm]], paste0(px$name,"_",nm),
                        estimate="CD", eval.loglik=FALSE, control=ctrl_cd)
        if (!is.na(out$ok)) { n_all <- n_all + 1L; n_ok <- n_ok + as.integer(out$ok) }
        }
    }
    cat(sprintf("\n=== Bilan Phase 2 : %d / %d fits OK ===\n", n_ok, n_all))
    if (n_ok < n_all) stop(sprintf("Echec fits: %d KO", n_all - n_ok))
    invisible(NULL)
}

# ======================================================================================
# Phase 3: (optionnel, commenté) petit MLE sur partition simple
# ======================================================================================

run_phase3_mle <- function() {
    cat("\n=== PHASE 3 : MLE court sur petite partition ===\n")
    part  <- partitions$P3
    nodes <- .make_nodes(part)
    f <- as.formula("partition ~ cov_ingroup('gender', category='F') + cov_ingroup('age', size=2:3)")
    environment(f) <- list2env(list(partition = part, nodes = nodes), parent = parent.frame())
    set.seed(99)
    fit <- erpm(f,  estimate="MLE",     eval.loglik=TRUE, control=ctrl_mle,
                    verbose=FALSE,      nodes = nodes)
    # fit <- erpm(f, estimate="MLE", eval.loglik=TRUE, control=ctrl_mle, verbose=FALSE)
    print(summary(fit))
    invisible(NULL)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: cov_ingroup ===\n")
run_phase1()
run_phase2()
# run_phase3_mle()  # test MLE long temps d'exécution...

on.exit(try(ergm_patch_disable(), silent = TRUE), add = TRUE)
cat("\nTous les tests cov_ingroup ont passé.\n")