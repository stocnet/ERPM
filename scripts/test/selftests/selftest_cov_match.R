# ======================================================================================
# Fichier : scripts/test/selftests/selftest_cov_match.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_match`
# Exécution: Rscript scripts/test/selftests/selftest_cov_match.R
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
        score  = sample(1:5, n, replace = TRUE),  # numérique
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
summary_on_bipartite_network <- function(part, nodes, rhs_txt) {
    nw <- .erpm_build_bipartite_nw(part, nodes)
    f  <- .formula_nw(nw, rhs_txt)
    as.numeric(suppressMessages(summary(f)))
}

# Summary côté ERPM (traduction)
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
# Attentes analytiques pour cov_match (réplique la sémantique implémentée)
# ======================================================================================

.choose <- function(n, k) if (n >= k) choose(n, k) else 0

expected_cov_match <- function(part, vals, k = 2L, category = NULL, normalized = c("none","by_group","global")) {
    normalized <- match.arg(normalized)
    gid <- as.integer(part)
    split_idx <- split(seq_along(gid), gid)
    N1 <- length(vals)

    stat_group <- function(ix) {
        v <- vals[ix]; n_g <- length(ix)
        if (!is.null(category)) {
        n_k <- sum(v == category, na.rm = TRUE)
        if (k == 1L && normalized == "by_group") {
            # Cas spécial observé: somme_g 1_{n_{g,κ} >= 1}
            return(as.numeric(n_k >= 1L))
        }
        clq <- .choose(n_k, k)
        } else {
        # somme des cliques homogènes toutes catégories
        tbl <- table(v, useNA = "no")
        clq <- sum(vapply(as.integer(tbl), .choose, numeric(1), k = k))
        }
        if (normalized == "none") return(clq)
        if (normalized == "by_group") {
        den <- .choose(n_g, k)
        if (den == 0) return(0)
        return(clq / den)
        }
        # normalized == "global"
        den <- .choose(N1, k)
        if (den == 0) return(0)
        return(clq / den)
    }

    sum(vapply(split_idx, stat_group, numeric(1)))
}

# ======================================================================================
# Phase 1: Summary — attentes numériques explicites + erreurs attendues
# ======================================================================================

run_phase1_summary_expected <- function() {
    cat("=== PHASE 1 : Summary avec attentes explicites ===\n")

    # 1) Cas A (déterministe comme tes logs)
    # - k=2 non normalisé
    exp_A_k2 <- expected_cov_match(partA, nodesA$sexe, k = 2, normalized = "none")
    s_A_k2   <- summary_on_bipartite_network(partA, nodesA, "cov_match('sexe', clique_size = 2)")
    cat(sprintf("[A] cov_match('sexe',k=2)            obtenu=%g  attendu=%g\n", s_A_k2, exp_A_k2))
    stopifnot(isTRUE(all.equal(as.numeric(s_A_k2), as.numeric(exp_A_k2))))

    # - k=3 non normalisé
    exp_A_k3 <- expected_cov_match(partA, nodesA$sexe, k = 3, normalized = "none")
    s_A_k3   <- summary_on_bipartite_network(partA, nodesA, "cov_match('sexe', clique_size = 3)")
    cat(sprintf("[A] cov_match('sexe',k=3)            obtenu=%g  attendu=%g\n", s_A_k3, exp_A_k3))
    stopifnot(isTRUE(all.equal(as.numeric(s_A_k3), as.numeric(exp_A_k3))))

    # - ciblage κ='F', k=2
    exp_A_F2 <- expected_cov_match(partA, nodesA$sexe, k = 2, category = "F", normalized = "none")
    s_A_F2   <- summary_on_bipartite_network(partA, nodesA, "cov_match('sexe', clique_size = 2, category = 'F')")
    cat(sprintf("[A] cov_match('sexe==F',k=2)         obtenu=%g  attendu=%g\n", s_A_F2, exp_A_F2))
    stopifnot(isTRUE(all.equal(as.numeric(s_A_F2), as.numeric(exp_A_F2))))

    # - normalisation by_group, k=2
    exp_A_bg2 <- expected_cov_match(partA, nodesA$sexe, k = 2, normalized = "by_group")
    s_A_bg2   <- summary_on_bipartite_network(partA, nodesA, "cov_match('sexe', clique_size = 2, normalized = 'by_group')")
    cat(sprintf("[A] cov_match('sexe',k=2,by_group)   obtenu=%.6f  attendu=%.6f\n", s_A_bg2, exp_A_bg2))
    stopifnot(isTRUE(all.equal(as.numeric(s_A_bg2), as.numeric(exp_A_bg2))))

    # - normalisation globale, k=2
    exp_A_gl2 <- expected_cov_match(partA, nodesA$sexe, k = 2, normalized = "global")
    s_A_gl2   <- summary_on_bipartite_network(partA, nodesA, "cov_match('sexe', clique_size = 2, normalized = 'global')")
    cat(sprintf("[A] cov_match('sexe',k=2,global)     obtenu=%.7f  attendu=%.7f\n", s_A_gl2, exp_A_gl2))
    stopifnot(isTRUE(all.equal(as.numeric(s_A_gl2), as.numeric(exp_A_gl2))))

    # 2) Cas C (k=1 autorisés + erreurs attendues)
    # - k=1, normalized='by_group' : chaque groupe non vide contribue 1
    exp_C_k1_bg <- expected_cov_match(partC, nodesC$sexe, k = 1, normalized = "by_group")
    s_C_k1_bg   <- summary_on_bipartite_network(partC, nodesC, "cov_match('sexe', clique_size = 1, normalized = 'by_group')")
    cat(sprintf("[C] cov_match('sexe',k=1,by_group)   obtenu=%g  attendu=%g\n", s_C_k1_bg, exp_C_k1_bg))
    stopifnot(isTRUE(all.equal(as.numeric(s_C_k1_bg), as.numeric(exp_C_k1_bg))))

    # - k=1 ciblé κ='G1' + by_group : somme_g 1_{n_{g,κ} >= 1}
    exp_C_g1_bg <- expected_cov_match(partC, nodesC$grade, k = 1, category = "G1", normalized = "by_group")
    s_C_g1_bg   <- summary_on_bipartite_network(partC, nodesC, "cov_match('grade', clique_size = 1, category = 'G1', normalized = 'by_group')")
    cat(sprintf("[C] cov_match('grade==G1',k=1,bg)    obtenu=%g  attendu=%g\n", s_C_g1_bg, exp_C_g1_bg))
    stopifnot(isTRUE(all.equal(as.numeric(s_C_g1_bg), as.numeric(exp_C_g1_bg))))

    # - catégorie absente => 0
    exp_C_abs <- expected_cov_match(partC, nodesC$grade, k = 2, category = "G999", normalized = "none")
    s_C_abs   <- summary_on_bipartite_network(partC, nodesC, "cov_match('grade', clique_size = 2, category = 'G999')")
    cat(sprintf("[C] cov_match('grade==G999',k=2)     obtenu=%g  attendu=%g\n", s_C_abs, exp_C_abs))
    stopifnot(isTRUE(all.equal(as.numeric(s_C_abs), as.numeric(exp_C_abs))))

    # --------------------------------------------------------------------
    # BLOC D'ERREURS ATTENDUES POUR k=1, normalized != 'by_group'
    # --------------------------------------------------------------------
    cat("\n--- Vérification des erreurs attendues pour k=1 ---\n")

    # 1) k=1, normalized = 'none' (par défaut) -> doit ERREUR
    err <- NULL
    tryCatch({
        summary_on_bipartite_network(partC, nodesC, "cov_match('sexe', clique_size = 1)")
    }, error = function(e) err <<- conditionMessage(e))
    if (is.null(err)) {
        stop("Erreur attendue non levée pour cov_match('sexe', clique_size = 1) (normalized='none').")
    } else {
        cat("[OK] Erreur attendue pour k=1, normalized='none':\n     ", err, "\n")
    }

    # 2) k=1, normalized = 'global' -> doit ERREUR
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
    invisible(NULL)
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
        cat(sprintf("[EQUIV] n=%-3d RHS=%-55s net=%s | erpm=%s\n",
                    length(part), rhs, paste(s_net, collapse=","), paste(s_erpm, collapse=",")))
        if (length(s_net) != length(s_erpm)) stop("Longueur de statistique différente.")
        if (!all(is.finite(s_net)) || !all(is.finite(s_erpm))) stop("Stat non finie.")
        if (!isTRUE(all.equal(as.numeric(s_net), as.numeric(s_erpm)))) {
        stop(sprintf("Mismatch summary net vs ERPM pour RHS=%s", rhs))
        }
    }
    TRUE
}

run_phase2_summary_equiv <- function() {
    cat("\n=== PHASE 2 : Summary(nw) vs Summary(ERPM-traduit) ===\n")
    for (nm in names(partitions)) {
        part  <- partitions[[nm]]
        nodes <- .make_nodes(part)
        cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                    nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))
        ok <- check_summary_equivalence(part, nodes, cases_equiv)
        if (!ok) stop("Equivalence summary échouée.")
    }
    # Cas déterministes A et C
    cat("\n--- Partition A (déterministe) ---\n")
    stopifnot(check_summary_equivalence(partA, nodesA, cases_equiv[1:3]))
    cat("\n--- Partition C (déterministe) ---\n")
    stopifnot(check_summary_equivalence(partC, nodesC, cases_equiv[c(1,3,5)]))
    cat("=== Phase 2 OK ===\n")
    invisible(NULL)
}

# ======================================================================================
# Phase 3: Fits courts via erpm() — MLE + loglik
# ======================================================================================

ctrl_mle <- control.ergm(
    MCMLE.maxit   = 20,
    MCMC.burnin   = 10000,
    MCMC.interval = 2000,
    parallel      = 0
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

run_fit <- function(part, nodes, rhs, tag) {
    if (!exists("erpm", mode = "function")) {
        cat(sprintf("[ERPM-FIT %-12s] SKIP (erpm() indisponible)\n", tag))
        return(list(ok = NA, fit = NULL, coef = NA, expected_error = FALSE))
    }

    f <- as.formula(paste0("partition ~ ", rhs))
    environment(f) <- list2env(list(partition = part, nodes = nodes), parent = parent.frame())

    # Stat observée
    s_obs <- try(summary_on_erpm_translation(part, nodes, rhs), silent = TRUE)
    if (!inherits(s_obs, "try-error")) {
        cat(sprintf("[ERPM-FIT %-12s] stat_observee=%s\n", tag, paste(format(s_obs), collapse=", ")))
    } else {
        cat(sprintf("[ERPM-FIT %-12s] stat_observee=NA (%s)\n", tag, conditionMessage(attr(s_obs, "condition"))))
    }
    obs_zero <- (!inherits(s_obs, "try-error") &&
                all(is.finite(s_obs)) &&
                length(s_obs) == 1L &&
                as.numeric(s_obs) == 0)

    # Exécution avec capture des warnings
    res <- .with_warning_capture(
        try(erpm(f, estimate = "MLE", eval.loglik = TRUE, control = ctrl_mle,
                verbose = FALSE, nodes = nodes), silent = TRUE)
    )
    fit   <- res$value
    warns <- res$warnings
    if (length(warns)) {
        cat(sprintf("[ERPM-FIT %-12s] WARNINGS (%d):\n", tag, length(warns)))
        for (w in unique(warns)) cat("  - ", w, "\n", sep = "")
    }

    # Cas erreur d’exécution
    if (inherits(fit, "try-error")) {
        msg <- as.character(fit)

        # Cas prévu: stat_observee = 0 + matrice avec diag négative
        if (obs_zero && grepl("Matrix .* negative elements on the diagonal", msg)) {
        cat(sprintf("[ERPM-FIT %-12s] ERREUR ATTENDUE (stat_observee=0): %s\n", tag, msg))
        return(list(ok = TRUE, fit = NULL, coef = NA, expected_error = TRUE))
        }

        cat(sprintf("[ERPM-FIT %-12s] ERREUR: %s\n", tag, msg))
        return(list(ok = FALSE, fit = NULL, coef = NA, expected_error = FALSE))
    }

    # Vérification coefficients
    cf <- try(stats::coef(fit), silent = TRUE)
    ok <- !(inherits(cf, "try-error")) && all(is.finite(cf))
    cat(sprintf("[ERPM-FIT %-12s] coef finies: %s | coef=%s\n",
                tag, if (ok) "OK" else "KO",
                if (ok) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

    list(ok = ok, fit = fit, coef = cf, expected_error = FALSE)
}

run_phase3_erpm_fits <- function() {
    cat("\n=== PHASE 3 : Fits erpm() (MLE + loglik) ===\n")

    fits <- list(
        A_k2       = run_fit(partA, nodesA, "cov_match('sexe', clique_size = 2)", "A_k2"),
        A_k2F      = run_fit(partA, nodesA, "cov_match('sexe', clique_size = 2, category = 'F')", "A_k2F"),
        A_k2_bg    = run_fit(partA, nodesA, "cov_match('sexe', clique_size = 2, normalized = 'by_group')", "A_k2_bg"),
        B_k23_H    = run_fit(partitions$B, .make_nodes(partitions$B), "cov_match('sexe', clique_size = c(2,3), category = 'H')", "B_k23_H"),
        B_k2_RHbg  = run_fit(partitions$B, .make_nodes(partitions$B), "cov_match('dept', clique_size = 2, category = 'RH', normalized = 'by_group')", "B_k2_RHbg")
    )

    ok <- vapply(fits, function(x) isTRUE(x$ok), logical(1))
    n_ok  <- sum(ok, na.rm = TRUE)
    n_tot <- sum(!is.na(ok))
    cat(sprintf("\n=== Bilan fits erpm() : %d / %d OK ===\n", n_ok, n_tot))

    cat("\n=== Résumés des fits ERPM réussis ===\n")
    for (nm in names(fits)) {
        fx <- fits[[nm]]
        if (isTRUE(fx$ok) && inherits(fx$fit, "ergm")) {
        cat(sprintf("\n--- Résumé fit %s ---\n", nm))
        print(summary(fx$fit))
        }
    }

    if (n_ok < n_tot) stop(sprintf("Echec fits: %d KO", n_tot - n_ok))
    invisible(fits)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: cov_match ===\n")
run_phase1_summary_expected()
run_phase2_summary_equiv()
res_fits <- run_phase3_erpm_fits()

ergm_patch_disable()
cat("\nTous les tests cov_match ont passé.\n")