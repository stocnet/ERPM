# ==============================================================================
# File    : scripts/test/selftests/selftest_erpm_long_vs_erpm_static_PLS.R
# Object  : Vérifier (T=2, PLS, statique) que erpm_long() == 2x erpm()
#           pour les effets : cov_ingroup, squared_sizes, dyadcov
#
# Usage   : Rscript scripts/test/selftests/selftest_erpm_long_vs_erpm_static_PLS.R
# Run from: racine du package
# ==============================================================================

# --------------------------------------------------------------------------------------
# Préambule
# --------------------------------------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))
options(ergm.loglik.warn_dyads = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("network",   quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",      quietly = TRUE)) stop("Package 'ergm' requis.")
  if (!requireNamespace("devtools",  quietly = TRUE)) stop("Package 'devtools' requis.")
  if (!requireNamespace("rprojroot", quietly = TRUE)) stop("Package 'rprojroot' requis.")
})

suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

devtools::load_all(".")

if (!exists("erpm_long", mode = "function")) stop("erpm_long() introuvable.")
if (!exists("erpm",      mode = "function")) stop("erpm() introuvable.")

# --------------------------------------------------------------------------------------
# Patch ERGM (optionnel, ne doit pas faire échouer le selftest)
# --------------------------------------------------------------------------------------
root <- tryCatch(
  rprojroot::find_root(rprojroot::is_r_package),
  error = function(e) getwd()
)

.patch_enabled <- FALSE
patch_path <- file.path(root, "scripts", "ergm_patch.R")
if (!file.exists(patch_path)) {
  message("[selftest_erpm_long] scripts/ergm_patch.R introuvable: selftest continue sans patch.")
} else {
  source(patch_path, local = FALSE)
  if (exists("ergm_patch_enable", mode = "function")) {
    tryCatch(
      {
        ergm_patch_enable()
        .patch_enabled <- TRUE
      },
      error = function(e) {
        message("[selftest_erpm_long] ergm_patch_enable() failed, continuing without patch.\n",
                "  reason: ", conditionMessage(e))
        .patch_enabled <- FALSE
      }
    )
  } else {
    message("[selftest_erpm_long] ergm_patch_enable() introuvable: selftest continue sans patch.")
  }
}

# --------------------------------------------------------------------------------------

.stopf <- function(...) stop(sprintf(...), call. = FALSE)

.fit_sig <- function(fit) {
  if (!inherits(fit, "ergm")) .stopf("Objet inattendu (pas un 'ergm'): %s", paste(class(fit), collapse = "/"))
  list(coef = stats::coef(fit), ll = as.numeric(stats::logLik(fit)))
}

.must_same <- function(a, b, label) {
  ok1 <- isTRUE(all.equal(a$coef, b$coef, tolerance = 1e-10))
  ok2 <- isTRUE(all.equal(a$ll,   b$ll,   tolerance = 1e-10))
  if (!ok1 || !ok2) {
    cat("\n[FAIL] ", label, "\n", sep = "")
    cat("  coef(A): ", paste(signif(a$coef, 12), collapse = ", "), "\n", sep = "")
    cat("  coef(B): ", paste(signif(b$coef, 12), collapse = ", "), "\n", sep = "")
    cat("  logLik(A): ", signif(a$ll, 12), "\n", sep = "")
    cat("  logLik(B): ", signif(b$ll, 12), "\n", sep = "")
    .stopf("Mismatch: %s", label)
  }
  cat("[OK] ", label, "\n", sep = "")
  TRUE
}

.get_long_fits <- function(out, T) {
  if (inherits(out, "erpm_long") && !is.null(out$fits) && length(out$fits) == T) {
    if (all(vapply(out$fits, inherits, logical(1), what = "ergm"))) return(out$fits)
  }
  fits <- list()
  rec <- function(x) {
    if (inherits(x, "ergm")) fits[[length(fits) + 1L]] <<- x
    if (is.list(x)) for (i in seq_along(x)) rec(x[[i]])
  }
  rec(out)
  if (length(fits) < T) .stopf("erpm_long(): attendu %d fits 'ergm', obtenu %d.", T, length(fits))
  fits[seq_len(T)]
}

# ------------------------------- Données test ---------------------------------
n1 <- 10L
P1 <- c(1,1,1, 2,2, 3,3,3,3, 4)
P2 <- c(1,1, 2,2,2, 3,3, 4,4,4)

x <- c(0.2, 1.5, -0.3, 2.0, 0.7, 1.1, -0.5, 0.0, 0.9, 1.8)
nodes_df <- data.frame(id = seq_len(n1), x = x)

set.seed(123)
Z <- matrix(rnorm(n1 * n1), n1, n1)
diag(Z) <- 0
Z <- (Z + t(Z)) / 2

# Seed unique (erpm_long n'accepte qu'un seed scalaire).
seed0 <- 424242L

rhs_list <- list(
  cov_ingroup   = 'cov_ingroup("x")',
  squared_sizes = 'squared_sizes',
  dyadcov       = 'dyadcov("Z1")'
)

# ------------------------------- Runner ---------------------------------------
run_one <- function(effect_name, rhs) {
  cat("\n============================================================\n")
  cat("=== TEST: ", effect_name, " | RHS = ", rhs, " ===\n", sep = "")
  cat("============================================================\n")

  partitions <- list(P1 = P1, P2 = P2)

  f_long <- stats::as.formula(paste0("partitions ~ ", rhs))
  environment(f_long) <- list2env(list(partitions = partitions), parent = parent.frame())

  # erpm_long(): on passe une matrice nue même si RHS ne consomme pas de dyads.
  set.seed(seed0)
  out_long <- erpm_long(
    formula   = f_long,
    mode      = "sequentiel",
    eval.call = TRUE,
    verbose   = TRUE,
    debug     = FALSE,
    nodes     = nodes_df,
    dyads     = Z,
    seed      = seed0
  )
  long_fits <- .get_long_fits(out_long, T = 2L)

  # erpm(): on reste strict -> dyads doit être une liste nommée si matrice fournie.
  dyads_erpm <- if (identical(effect_name, "dyadcov")) list(Z1 = Z) else list(dummy = Z)

  f1 <- stats::as.formula(paste0("P1 ~ ", rhs))
  environment(f1) <- environment()
  f2 <- stats::as.formula(paste0("P2 ~ ", rhs))
  environment(f2) <- environment()

  set.seed(seed0)
  fit1 <- erpm(
    f1,
    eval.call = TRUE,
    verbose   = TRUE,
    nodes     = nodes_df,
    dyads     = dyads_erpm,
    seed      = seed0
  )

  set.seed(seed0)
  fit2 <- erpm(
    f2,
    eval.call = TRUE,
    verbose   = TRUE,
    nodes     = nodes_df,
    dyads     = dyads_erpm,
    seed      = seed0
  )

  .must_same(.fit_sig(long_fits[[1]]), .fit_sig(fit1), paste0(effect_name, " | t=1 (erpm_long vs erpm)"))
  .must_same(.fit_sig(long_fits[[2]]), .fit_sig(fit2), paste0(effect_name, " | t=2 (erpm_long vs erpm)"))

  invisible(TRUE)
}

for (nm in names(rhs_list)) run_one(nm, rhs_list[[nm]])

cat("\n=== ALL TESTS PASSED ===\n")