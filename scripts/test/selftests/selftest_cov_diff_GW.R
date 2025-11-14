# ==============================================================================
# Fichier : scripts/test/selftests/selftest_cov_diff_GW.R
# Objet   : Self-test autonome pour l'effet ERPM/ERGM `cov_diff_GW`
# Exécution: Rscript scripts/test/selftests/selftest_cov_diff_GW.R
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
    ergm_patch_enable()
}

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
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
log_path <- file.path(root, "scripts", "test", "selftests", "selftest_cov_diff_GW.log")
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
# Données de test pour cov_diff_GW (reprennent les cas A/B/C de cov_diff)
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
    set.seed(300 + n)
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
# Attentes analytiques pour cov_diff et cov_diff_GW
# ======================================================================================

# Stat de base pour cov_diff (k fixé)
expected_cov_diff <- function(part, x, clique_size = 2L, normalized = FALSE) {
    stopifnot(length(part) == length(x))
    k <- as.integer(clique_size)
    if (k < 2L) stop("expected_cov_diff: clique_size doit être >= 2.")

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

        if (!normalized) {
            tot <- tot + Sg
        } else {
            denom <- choose(ng, k)
            tot <- tot + Sg / denom
        }
    }

    tot
}

# Stat globale T_GW(p;x;lambda) = sum_{k=2..Kmax} (-1/lambda)^{k-1} c_k
expected_cov_diff_GW <- function(part, x, lambda) {
    stopifnot(length(part) == length(x))
    lam <- as.numeric(lambda)
    if (any(lam <= 1)) stop("expected_cov_diff_GW: lambda doit être > 1.")
    Kmax <- max(table(part))
    if (Kmax < 2L) return(rep(0, length(lam)))

    # Pré-calcul des c_k
    ck <- vapply(2:Kmax, function(k) {
        expected_cov_diff(part, x, clique_size = k, normalized = FALSE)
    }, numeric(1))
    k_vec <- 2:Kmax

    sapply(lam, function(l) {
        weights <- (-1 / l)^(k_vec - 1L)
        sum(weights * ck)
    })
}

# ======================================================================================
# Phase 1: Summary — attentes explicites (analytique vs summary)
# ======================================================================================

run_phase1_summary_expected_cov_diff_GW <- function() {
    cat("=== PHASE 1 : Summary cov_diff_GW avec attentes explicites ===\n")

    lambdas <- c(1.5, 2, 4)

    # ---------- Cas A ----------
    for (l in lambdas) {
        expA <- expected_cov_diff_GW(partA, nodesA$score, lambda = l)
        sA   <- summary_on_bipartite_network(partA, nodesA,
                    sprintf("cov_diff_GW('score', lambda = %g)", l))
        cat(sprintf("[A] cov_diff_GW('score', lambda=%g)  obtenu=%g  attendu=%g\n",
                    l, sA, expA))
        stopifnot(isTRUE(all.equal(as.numeric(sA), as.numeric(expA))))
    }

    # ---------- Cas B ----------
    for (l in lambdas) {
        expB <- expected_cov_diff_GW(partB, nodesB$score, lambda = l)
        sB   <- summary_on_bipartite_network(partB, nodesB,
                    sprintf("cov_diff_GW('score', lambda = %g)", l))
        cat(sprintf("[B] cov_diff_GW('score', lambda=%g)  obtenu=%g  attendu=%g\n",
                    l, sB, expB))
        stopifnot(isTRUE(all.equal(as.numeric(sB), as.numeric(expB))))
    }

    # ---------- Cas C ----------
    for (l in lambdas) {
        expC <- expected_cov_diff_GW(partC, nodesC$score, lambda = l)
        sC   <- summary_on_bipartite_network(partC, nodesC,
                    sprintf("cov_diff_GW('score', lambda = %g)", l))
        cat(sprintf("[C] cov_diff_GW('score', lambda=%g)  obtenu=%g  attendu=%g\n",
                    l, sC, expC))
        stopifnot(isTRUE(all.equal(as.numeric(sC), as.numeric(expC))))
    }

    # ---------- Test vectorisé lambda = c(2,4) sur A ----------
    expA_v <- expected_cov_diff_GW(partA, nodesA$score, lambda = c(2, 4))
    sA_v   <- summary_on_bipartite_network(partA, nodesA,
                    "cov_diff_GW('score', lambda = c(2,4))")
    cat(sprintf("[A] cov_diff_GW('score', lambda=c(2,4))  obtenu=%s  attendu=%s\n",
                paste(sA_v, collapse = ","),
                paste(expA_v, collapse = ",")))
    stopifnot(isTRUE(all.equal(as.numeric(sA_v), as.numeric(expA_v))))

    cat("\n=== Phase 1 OK ===\n")
    invisible(NULL)
}

# ======================================================================================
# Phase 2: Summary — équivalence réseau explicite vs ERPM-traduit
# ======================================================================================

cases_equiv_diff_GW <- c(
    "cov_diff_GW('score', lambda = 2)",
    "cov_diff_GW('score', lambda = 3)",
    "cov_diff_GW('score', lambda = c(2,4))"
)

check_summary_equivalence_cov_diff_GW <- function(part, nodes, rhs_vec) {
    for (rhs in rhs_vec) {
        s_net  <- summary_on_bipartite_network(part, nodes, rhs)
        s_erpm <- summary_on_erpm_translation(part, nodes, rhs)
        cat(sprintf("[EQUIV-GW] n=%-3d RHS=%-45s net=%s | erpm=%s\n",
                    length(part), rhs,
                    paste(s_net, collapse=","), paste(s_erpm, collapse=",")))
        if (length(s_net) != length(s_erpm)) stop("Longueur de statistique différente.")
        if (!all(is.finite(s_net)) || !all(is.finite(s_erpm))) stop("Stat non finie.")
        if (!isTRUE(all.equal(as.numeric(s_net), as.numeric(s_erpm)))) {
            stop(sprintf("Mismatch summary net vs ERPM pour RHS=%s", rhs))
        }
    }
    TRUE
}

run_phase2_summary_equiv_cov_diff_GW <- function() {
    cat("\n=== PHASE 2 : Summary(nw) vs Summary(ERPM-traduit) pour cov_diff_GW ===\n")

    # Cas aléatoires (scores simulés)
    for (nm in names(partitions)) {
        part  <- partitions[[nm]]
        nodes <- .make_nodes_numeric(part)
        cat(sprintf("\n--- Partition %s ---  n=%d | groupes=%d | tailles: %s\n",
                    nm, length(part), length(unique(part)), paste(sort(table(part)), collapse=",")))
        ok <- check_summary_equivalence_cov_diff_GW(part, nodes, cases_equiv_diff_GW)
        if (!ok) stop("Equivalence summary échouée pour cov_diff_GW.")
    }

    # Cas déterministes A,B,C sur 'score'
    cat("\n--- Partition A (déterministe) ---\n")
    stopifnot(check_summary_equivalence_cov_diff_GW(partA, nodesA, cases_equiv_diff_GW))

    cat("\n--- Partition B (déterministe) ---\n")
    stopifnot(check_summary_equivalence_cov_diff_GW(partB, nodesB, cases_equiv_diff_GW))

    cat("\n--- Partition C (déterministe) ---\n")
    stopifnot(check_summary_equivalence_cov_diff_GW(partC, nodesC, cases_equiv_diff_GW))

    cat("=== Phase 2 OK ===\n")
    invisible(NULL)
}

# ======================================================================================
# Phase 3: Fits courts via erpm() — MLE + loglik (avec terme structurel)
# ======================================================================================

ctrl_mle <- control.ergm(
    MCMLE.maxit     = 20,
    MCMC.burnin     = 10000,
    MCMC.interval   = 1000,
    MCMC.samplesize = 50000,
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

run_fit_cov_diff_GW <- function(part, nodes, rhs, tag) {
    if (!exists("erpm", mode = "function")) {
        cat(sprintf("[ERPM-FIT-COVGW %-12s] SKIP (erpm() indisponible)\n", tag))
        return(list(ok = NA, fit = NULL, coef = NA, expected_error = FALSE))
    }

    f <- as.formula(paste0("partition ~ ", rhs))
    environment(f) <- list2env(list(partition = part, nodes = nodes), parent = parent.frame())

    s_obs <- try(summary_on_erpm_translation(part, nodes, rhs), silent = TRUE)
    if (!inherits(s_obs, "try-error")) {
        cat(sprintf("[ERPM-FIT-COVGW %-12s] stat_observee=%s\n", tag, paste(format(s_obs), collapse=", ")))
    } else {
        cat(sprintf("[ERPM-FIT-COVGW %-12s] stat_observee=NA (%s)\n",
                    tag, conditionMessage(attr(s_obs, "condition"))))
    }

    res <- .with_warning_capture(
        try(erpm(f, estimate = "MLE", eval.loglik = TRUE, control = ctrl_mle,
                 verbose = FALSE, nodes = nodes), silent = TRUE)
    )
    fit   <- res$value
    warns <- res$warnings
    if (length(warns)) {
        cat(sprintf("[ERPM-FIT-COVGW %-12s] WARNINGS (%d):\n", tag, length(warns)))
        for (w in unique(warns)) cat("  - ", w, "\n", sep = "")
    }

    if (inherits(fit, "try-error")) {
        msg <- as.character(fit)
        cat(sprintf("[ERPM-FIT-COVGW %-12s] ERREUR: %s\n", tag, msg))
        return(list(ok = FALSE, fit = NULL, coef = NA, expected_error = FALSE))
    }

    cf <- try(stats::coef(fit), silent = TRUE)
    ok <- !(inherits(cf, "try-error")) && all(is.finite(cf))
    cat(sprintf("[ERPM-FIT-COVGW %-12s] coef finies: %s | coef=%s\n",
                tag, if (ok) "OK" else "KO",
                if (ok) paste(format(as.numeric(cf)), collapse=", ") else "NA"))

    list(ok = ok, fit = fit, coef = cf, expected_error = FALSE)
}

run_phase3_erpm_fits_cov_diff_GW <- function() {
    cat("\n=== PHASE 3 : Fits erpm() (MLE + loglik) pour cov_diff_GW (avec squared_sizes) ===\n")

    fits <- list(
        A_l2 = run_fit_cov_diff_GW(
            partA, nodesA,
            "cov_diff_GW('score', lambda = 2)",
            "A_l2"
        ),
        A_l2_l4 = run_fit_cov_diff_GW(
            partA, nodesA,
            "cov_diff_GW('score', lambda = c(2,4))",
            "A_l2_l4"
        ),
        B_l2 = run_fit_cov_diff_GW(
            partB, nodesB,
            "cov_diff_GW('score', lambda = 2)",
            "B_l2"
        )
    )

    ok    <- vapply(fits, function(x) isTRUE(x$ok), logical(1))
    n_ok  <- sum(ok, na.rm = TRUE)
    n_tot <- sum(!is.na(ok))
    cat(sprintf("\n=== Bilan fits erpm() cov_diff_GW : %d / %d OK ===\n", n_ok, n_tot))

    cat("\n=== Résumés des fits ERPM cov_diff_GW réussis ===\n")
    for (nm in names(fits)) {
        fx <- fits[[nm]]
        if (isTRUE(fx$ok) && inherits(fx$fit, "ergm")) {
            cat(sprintf("\n--- Résumé fit %s ---\n", nm))
            print(summary(fx$fit))
        }
    }

    if (n_ok < n_tot) stop(sprintf("Echec fits cov_diff_GW: %d KO", n_tot - n_ok))
    invisible(fits)
}

# ======================================================================================
# Exécution
# ======================================================================================

set.seed(1)
cat("=== TEST ERPM: cov_diff_GW ===\n")
run_phase1_summary_expected_cov_diff_GW()
run_phase2_summary_equiv_cov_diff_GW()
res_fits_cov_diff_GW <- run_phase3_erpm_fits_cov_diff_GW()

if (exists("ergm_patch_disable")) ergm_patch_disable()

cat("\nTous les tests cov_diff_GW ont passé.\n")