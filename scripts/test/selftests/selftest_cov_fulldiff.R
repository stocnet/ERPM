# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cov_fulldiff.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_fulldiff`
# Exécution: Rscript scripts/test/selftests/selftest_cov_fulldiff.R
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
# Données de test pour cov_fulldiff (partitions un peu plus grandes)
# ======================================================================================

# PartA : 4 groupes (5,8,6,5) = 24 sommets
partA  <- c(
    rep(1, 5),
    rep(2, 8),
    rep(3, 6),
    rep(4, 5)
)
nodesA <- data.frame(
    label = paste0("A", seq_along(partA)),
    score = c(
        # g1 (variation faible)
        10, 11,  9, 10, 12,
        # g2 (variation forte)
         5,  7, 20, 15,  9,  6, 30, 12,
        # g3 (rampe modérée)
         0,  2,  4,  6,  8, 10,
        # g4 (niveau moyen)
        15, 16, 14, 17, 15
    ),
    stringsAsFactors = FALSE
)

# PartB : 6 groupes (6,6,5,4,3,2) = 26 sommets
partB <- c(
    rep(1, 6),
    rep(2, 6),
    rep(3, 5),
    rep(4, 4),
    rep(5, 3),
    rep(6, 2)
)
nodesB <- data.frame(
    label = paste0("B", seq_along(partB)),
    score = c(
        # g1 (quasi constant)
        10, 10, 11,  9, 10, 10,
        # g2 (très étalé)
         0, 40, 80,  5, 60,100,
        # g3 (rampe)
         3,  4,  5,  6,  7,
        # g4 (deux niveaux)
        20, 21, 19, 20,
        # g5 (petit groupe)
         5,  7,  6,
        # g6 (groupe de taille 2)
         0, 10
    ),
    stringsAsFactors = FALSE
)

# PartC : 4 groupes (5,7,8,4) = 24 sommets
partC <- c(
    rep(1, 5),
    rep(2, 7),
    rep(3, 8),
    rep(4, 4)
)
nodesC <- data.frame(
    label = paste0("C", seq_along(partC)),
    score = c(
        # g1 constant -> D=0
         5, 5, 5, 5, 5,
        # g2 variation modérée
         1, 2, 4, 3, 5, 6, 4,
        # g3 énorme outlier
         0, 0, 0, 0, 0, 0, 0,100,
        # g4 petite variation
         7, 8, 6, 7
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
    call_ergm <- erpm(f, eval_call = FALSE, verbose = FALSE, nodes = nodes)
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

    # ---------- Cas A (n=24, tailles 5,8,6,5) ----------
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

    # ---------- Cas B (n=26, tailles 6,6,5,4,3,2) ----------
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

    # ---------- Cas C (n=24, tailles 5,7,8,4) ----------
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

    # Cas déterministes A,B,C sur 'score'
    cat("\n--- Partition A (déterministe) ---\n")
    stopifnot(check_summary_equivalence_fulldiff(partA, nodesA,
                                                 cases_equiv_fulldiff[c(1,3,4)]))
    cat("\n--- Partition B (déterministe) ---\n")
    stopifnot(check_summary_equivalence_fulldiff(partB, nodesB,
                                                 cases_equiv_fulldiff[c(1,2,3)]))
    cat("\n--- Partition C (déterministe) ---\n")
    stopifnot(check_summary_equivalence_fulldiff(partC, nodesC,
                                                 cases_equiv_fulldiff[c(1,3,4)]))

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

run_fit_fulldiff <- function(part, nodes, rhs, tag) {
    if (!exists("erpm", mode = "function")) {
        cat(sprintf("[ERPM-FIT-FULLDIFF %-12s] SKIP (erpm() indisponible)\n", tag))
        return(list(ok = NA, fit = NULL, coef = NA, expected_error = FALSE))
    }

    f <- as.formula(paste0("partition ~ ", rhs))
    environment(f) <- list2env(list(partition = part, nodes = nodes), parent = parent.frame())

    s_obs <- try(summary_on_erpm_translation(part, nodes, rhs), silent = TRUE)
    if (!inherits(s_obs, "try-error")) {
        cat(sprintf("[ERPM-FIT-FULLDIFF %-12s] stat_observee=%s\n", tag, paste(format(s_obs), collapse=", ")))
    } else {
        cat(sprintf("[ERPM-FIT-FULLDIFF %-12s] stat_observee=NA (%s)\n", tag, conditionMessage(attr(s_obs, "condition"))))
    }

    res <- .with_warning_capture(
        try(erpm(f, estimate = "MLE", eval.loglik = TRUE, control = ctrl_mle,
                 verbose = FALSE, nodes = nodes), silent = TRUE)
    )
    fit   <- res$value
    warns <- res$warnings
    if (length(warns)) {
        cat(sprintf("[ERPM-FIT-FULLDIFF %-12s] WARNINGS (%d):\n", tag, length(warns)))
        for (w in unique(warns)) cat("  - ", w, "\n", sep = "")
    }

    if (inherits(fit, "try-error")) {
        msg <- as.character(fit)
        cat(sprintf("[ERPM-FIT-FULLDIFF %-12s] ERREUR: %s\n", tag, msg))
        return(list(ok = FALSE, fit = NULL, coef = NA, expected_error = FALSE))
    }

    cf <- try(stats::coef(fit), silent = TRUE)
    ok <- !(inherits(cf, "try-error")) && all(is.finite(cf))
    cat(sprintf("[ERPM-FIT-FULLDIFF %-12s] coef finies: %s | coef=%s\n",
                tag, if (ok) "OK" else "KO",
                if (ok) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

    list(ok = ok, fit = fit, coef = cf, expected_error = FALSE)
}

run_phase3_erpm_fits_fulldiff <- function() {
    cat("\n=== PHASE 3 : Fits erpm() (MLE + loglik) pour cov_fulldiff ===\n")

    fits <- list(
        A_all  = run_fit_fulldiff(
            partA, nodesA,
            "cov_fulldiff('score')",
            "A_all"
        ),
        A_mid  = run_fit_fulldiff(
            partA, nodesA,
            "cov_fulldiff('score', size = c(5,6,8))",
            "A_mid"
        ),
        B_all  = run_fit_fulldiff(
            partB, nodesB,
            "cov_fulldiff('score')",
            "B_all"
        )
    )

    ok    <- vapply(fits, function(x) isTRUE(x$ok), logical(1))
    n_ok  <- sum(ok, na.rm = TRUE)
    n_tot <- sum(!is.na(ok))
    cat(sprintf("\n=== Bilan fits erpm() cov_fulldiff : %d / %d OK ===\n", n_ok, n_tot))

    cat("\n=== Résumés des fits ERPM cov_fulldiff réussis ===\n")
    for (nm in names(fits)) {
        fx <- fits[[nm]]
        if (isTRUE(fx$ok) && inherits(fx$fit, "ergm")) {
            cat(sprintf("\n--- Résumé fit %s ---\n", nm))
            print(summary(fx$fit))
        }
    }

    if (n_ok < n_tot) stop(sprintf("Echec fits cov_fulldiff: %d KO", n_tot - n_ok))
    invisible(fits)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: cov_fulldiff ===\n")
run_phase1_summary_expected_fulldiff()
run_phase2_summary_equiv_fulldiff()
res_fits_fulldiff <- run_phase3_erpm_fits_fulldiff()

if (exists("ergm_patch_disable")) ergm_patch_disable()

cat("\nTous les tests cov_fulldiff ont passé.\n")