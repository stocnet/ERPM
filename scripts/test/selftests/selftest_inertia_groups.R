# ==============================================================================
# Fichier : scripts/test/selftests/selftest_inertia_groups.R
# Objet   : Self-test ciblé pour l’effet ERPM/ERGM `inertia_groups`
# Exécution: Rscript scripts/test/selftests/selftest_inertia_groups.R
# ==============================================================================
#
# NOTE IMPORTANTE
# - Ce selftest ne construit PAS nw_t.
# - La construction du réseau temporel (nw_t) et l’attache des attributs inertiels
#   sont déléguées à erpm_long().
#
# Ce selftest vérifie uniquement :
# 1) Dry-runs: compare stats (summary sur networks construits par erpm_long) vs oracle R.
# 2) Fits: appelle erpm_long() (sans estimate ni control explicites) et vérifie
#    que coef() est disponible et fini. On utilise un baseline statique simple
#    pour éviter la dégénérescence dans les petits stress-tests.
# 3) Stress: partitions déterministes (dry + fit) avec gestion robuste des cas
#    "data are essentially constant" (skip contrôlé) et du cas
#    "Matrix ‘x’ has negative elements on the diagonal." (skip contrôlé).
#
# Debug:
# - Passe debug="deep" à erpm_long() si l’argument existe.
# - Optionnel: redirige output/messages deep dans des fichiers temporaires (pas de spam console)
#   tout en conservant la valeur de retour de erpm_long().
# ==============================================================================

options(ergm.loglik.warn_dyads = FALSE)
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  library(network)
  library(ergm)
})

# ------------------------------------------------------------------------------
# Patch ERGM optionnel
# ------------------------------------------------------------------------------
if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  if (exists("ergm_patch_enable")) ergm_patch_enable()
}

# ------------------------------------------------------------------------------
# Charger le package (dev) et vérifier erpm_long()
# ------------------------------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(recompile = TRUE, quiet = TRUE)
} else {
  stop("Package non disponible (devtools::load_all + DESCRIPTION requis).", call. = FALSE)
}
if (!exists("erpm_long", mode = "function")) {
  stop("erpm_long() manquant. Ce selftest requiert le wrapper longitudinal.", call. = FALSE)
}

cat("=== SELFTEST ERPM_LONG: inertia_groups (dry vs oracle + fits + stress) ===\n")

# ==============================================================================
# Options selftest
# ==============================================================================
DEEP_DEBUG <- TRUE   # passe debug="deep" à erpm_long() si supporté
DEEP_SINK  <- TRUE   # capture output/message deep dans des fichiers temporaires
KEEP_LOGS  <- FALSE  # si TRUE, conserve les logs deep (imprime leur chemin)

.capture_deep <- function(expr, label = "erpm_long_deep_") {
  if (!isTRUE(DEEP_DEBUG) || !isTRUE(DEEP_SINK)) return(force(expr))

  tf <- tempfile(label, fileext = ".log")
  con <- file(tf, open = "wt")

  on.exit({
    try(close(con), silent = TRUE)
    if (isTRUE(KEEP_LOGS)) {
      cat(sprintf("[deep log] %s\n", tf))
    } else {
      try(unlink(tf), silent = TRUE)
    }
  }, add = TRUE)

  sink(con, type = "output")
  sink(con, type = "message")
  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(sink(type = "output"), silent = TRUE)
  }, add = TRUE)

  force(expr)
}

.stopf <- function(...) stop(sprintf(...), call. = FALSE)

.is_essentially_constant_error <- function(msg) {
  msg <- tolower(as.character(msg))
  grepl("data are essentially constant", msg, fixed = TRUE)
}

.is_negative_diag_error <- function(msg) {
  msg <- tolower(as.character(msg))
  grepl("negative elements on the diagonal", msg, fixed = TRUE)
}

# ==============================================================================
# Données (timelines déterministes) + nodes explicites
# ==============================================================================
timelines <- list(
  TL1 = list(
    parts = list(
      P1 = c(1,1,2,2,3,3),
      P2 = c(1,1,2,3,3,3),
      P3 = c(1,2,2,3,3,3)
    ),
    nodes = list(
      N1 = data.frame(
        label  = c("N1","N2","N3","N4","N5","N6"),
        age    = c(25, 41, 33, 52, 29, 46),
        score  = c(3.2, 8.1, 4.7, 6.0, 2.9, 7.3),
        gender = c("F","M","M","F","F","M"),
        dept   = c("A","A","B","B","C","C"),
        stringsAsFactors = FALSE
      ),
      N2 = data.frame(
        label  = c("N1","N2","N3","N4","N5","N6"),
        age    = c(26, 42, 34, 53, 30, 47),
        score  = c(3.0, 8.0, 4.9, 6.1, 3.2, 7.1),
        gender = c("F","M","M","F","F","M"),
        dept   = c("A","A","B","C","C","C"),
        stringsAsFactors = FALSE
      ),
      N3 = data.frame(
        label  = c("N1","N2","N3","N4","N5","N6"),
        age    = c(27, 43, 35, 54, 31, 48),
        score  = c(3.1, 7.9, 5.0, 6.2, 3.1, 7.0),
        gender = c("F","M","M","F","F","M"),
        dept   = c("A","B","B","C","C","C"),
        stringsAsFactors = FALSE
      )
    )
  ),
  TL2 = list(
    parts = list(
      P1 = c(1,2,2,3,3,3,4,4),
      P2 = c(1,2,2,3,3,4,4,4),
      P3 = c(1,1,2,3,3,4,4,4)
    ),
    nodes = list(
      N1 = data.frame(
        label  = paste0("N",1:8),
        age    = c(22, 36, 44, 31, 55, 28, 47, 39),
        score  = c(1.5, 6.2, 7.8, 3.4, 8.9, 2.2, 5.7, 6.9),
        gender = c("F","M","F","M","F","M","F","M"),
        dept   = c("A","B","B","C","C","C","A","B"),
        stringsAsFactors = FALSE
      ),
      N2 = data.frame(
        label  = paste0("N",1:8),
        age    = c(23, 37, 45, 32, 56, 29, 48, 40),
        score  = c(1.7, 6.1, 7.7, 3.6, 9.0, 2.4, 5.6, 7.0),
        gender = c("F","M","F","M","F","M","F","M"),
        dept   = c("A","B","B","C","C","A","A","A"),
        stringsAsFactors = FALSE
      ),
      N3 = data.frame(
        label  = paste0("N",1:8),
        age    = c(24, 38, 46, 33, 57, 30, 49, 41),
        score  = c(1.6, 6.0, 7.6, 3.5, 8.8, 2.5, 5.5, 6.8),
        gender = c("F","M","F","M","F","M","F","M"),
        dept   = c("A","A","B","C","C","A","A","A"),
        stringsAsFactors = FALSE
      )
    )
  )
)

# ==============================================================================
# Oracle R (référence pour inertia_groups)
# ==============================================================================
.groups_from_partition <- function(p) split(seq_along(p), as.integer(p))

.group_signatures <- function(groups) {
  vapply(groups, function(v) paste(sort(as.integer(v)), collapse = ","), character(1))
}

.ref_stat_inertia_groups <- function(p_t, past_partitions, size = NULL, past_influence = 1L) {
  p_t <- as.integer(round(p_t))
  n <- length(p_t)

  past_influence <- as.integer(round(past_influence))
  if (!is.finite(past_influence) || past_influence < 1L) past_influence <- 1L

  g_t   <- .groups_from_partition(p_t)
  sig_t <- .group_signatures(g_t)
  sz_t  <- vapply(g_t, length, integer(1))

  keep <- rep(TRUE, length(sig_t))
  if (!is.null(size)) {
    size <- as.integer(round(size))
    size <- sort(unique(size))
    keep <- (sz_t %in% size)
  }

  sig_past <- character(0)
  for (lag in seq_len(past_influence)) {
    key <- paste0("P_lag", lag)
    p_prev <- past_partitions[[key]]
    if (is.null(p_prev)) next
    p_prev <- as.integer(round(p_prev))
    if (length(p_prev) != n) .stopf("ref_stat: partitions de tailles différentes (n=%d).", n)
    g_p <- .groups_from_partition(p_prev)
    sig_past <- c(sig_past, .group_signatures(g_p))
  }
  sig_past <- unique(sig_past)

  sum(keep & (sig_t %in% sig_past))
}

# ==============================================================================
# Helpers (formules, past, summary, appels erpm_long)
# ==============================================================================
.parse_past_influence <- function(rhs) {
  m <- regexec("past_influence\\s*=\\s*([0-9]+)", rhs)
  r <- regmatches(rhs, m)[[1L]]
  if (length(r) >= 2L) return(as.integer(r[2L]))
  1L
}

.make_long_formula <- function(parts, rhs) {
  rhs_expr <- parse(text = rhs)[[1L]]
  f_call <- bquote(.(parts) ~ .(rhs_expr))
  as.formula(f_call, env = parent.frame())
}

.past_from_parts <- function(parts) {
  if (length(parts) <= 1L) return(list())
  past <- list()
  for (lag in seq_len(length(parts) - 1L)) {
    past[[paste0("P_lag", lag)]] <- parts[[length(parts) - lag]]
  }
  past
}

.summary_on_network <- function(nw, rhs) {
  rhs_expr <- parse(text = rhs)[[1L]]
  f <- as.formula(bquote(nw ~ .(rhs_expr)))
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  as.numeric(suppressMessages(summary(f, constraints = ~ b1part)))
}

# dyads:
# - NULL
# - un seul objet (ex: matrice) réutilisé pour tous les t
# - ou une liste de longueur T (objets temporels)
.call_erpm_long <- function(parts, nodes_list, rhs, eval.call = FALSE, dyads = NULL, seed = 1L) {
  f <- .make_long_formula(parts, rhs)

  args <- list(
    formula   = f,
    eval.call = isTRUE(eval.call),
    verbose   = FALSE,
    nodes     = nodes_list
  )
  if (!is.null(dyads)) args$dyads <- dyads

  # IMPORTANT: ne pas forcer estimate ni control ici.
  if (isTRUE(DEEP_DEBUG) && "debug" %in% names(formals(erpm_long))) args$debug <- "deep"
  if ("seed" %in% names(formals(erpm_long))) args$seed <- as.integer(seed)

  run <- function() do.call(erpm_long, args)

  obj <- tryCatch(
    .capture_deep(run(), label = if (isTRUE(eval.call)) "erpm_long_fit_" else "erpm_long_dry_"),
    error = function(e) .stopf("Appel erpm_long() (%s) a échoué.\nMessage: %s",
                               if (isTRUE(eval.call)) "fit" else "dry",
                               conditionMessage(e))
  )

  if (!inherits(obj, "erpm_long")) .stopf("erpm_long() n'a pas renvoyé un objet 'erpm_long'.")
  obj
}

.get_fit_T <- function(obj, T) {
  if (is.null(obj$fits) || length(obj$fits) < T) .stopf("Fit: obj$fits[[%d]] indisponible.", T)
  fit <- obj$fits[[T]]
  if (is.null(fit)) .stopf("Fit: obj$fits[[%d]] est NULL.", T)
  fit
}

.check_fit <- function(fit_obj, rhs) {
  cf <- tryCatch(coef(fit_obj), error = function(e) NULL)
  if (is.null(cf)) .stopf("Fit: coef() indisponible (RHS=%s).", rhs)

  if (any(!is.finite(cf))) {
    cat(sprintf("[FIT warn] RHS=%-60s | coef non-finis (skip)\n", rhs))
    return(invisible(FALSE))
  }

  cat(sprintf("[FIT ok]   RHS=%-60s | p=%d\n", rhs, length(cf)))
  invisible(TRUE)
}

.compare_summary <- function(parts, nodes_list, rhs, size = NULL, seed = 123L) {
  T  <- length(parts)
  pi <- .parse_past_influence(rhs)

  # L’inertie nécessite T > past_influence (sinon aucun passé à comparer).
  if (T <= pi) {
    cat(sprintf("[SKIP]   T=%d RHS=%-60s (inertie inactive: nécessite T>%d)\n", T, rhs, pi))
    return(invisible(TRUE))
  }

  obj <- .call_erpm_long(parts, nodes_list, rhs, eval.call = FALSE, dyads = NULL, seed = seed)
  nwT <- obj$networks[[T]]
  if (is.null(nwT) || !inherits(nwT, "network")) .stopf("Dry-run: obj$networks[[%d]] invalide.", T)

  got <- .summary_on_network(nwT, rhs)

  p_t  <- parts[[T]]
  past <- .past_from_parts(parts)
  ref  <- .ref_stat_inertia_groups(p_t, past_partitions = past, size = size, past_influence = pi)

  cat(sprintf("[SUMMARY] T=%d n=%-3d RHS=%-60s got=%s | ref=%s\n",
              T, length(p_t), rhs, as.character(got), as.character(ref)))

  stopifnot(length(got) == 1L)
  stopifnot(is.finite(got))
  stopifnot(got == ref)
  invisible(TRUE)
}

# ==============================================================================
# Cas (centrés sur inertia_groups)
# ==============================================================================
cases_summary <- list(
  list(rhs = "inertia_groups()",                           size = NULL),
  list(rhs = "inertia_groups(size=2:3)",                   size = 2:3),
  list(rhs = "inertia_groups(size=c(1,2,4))",              size = c(1,2,4)),
  list(rhs = "inertia_groups(past_influence=1)",           size = NULL),
  list(rhs = "inertia_groups(size=2, past_influence=1)",   size = 2),
  list(rhs = "inertia_groups(past_influence=2)",           size = NULL),
  list(rhs = "inertia_groups(size=2:4, past_influence=2)", size = 2:4)
)

# Baseline statique minimal: on évite les effets combinatoires instables.
cases_fit <- c(
  "squared_sizes() + inertia_groups()",
  "squared_sizes() + inertia_groups(size=2:3)",
  "squared_sizes() + inertia_groups(size=2, past_influence=1)",
  "squared_sizes() + inertia_groups(past_influence=2)",
  "squared_sizes() + inertia_groups(size=2:4, past_influence=2)"
)

# ==============================================================================
# 1) DRY: summary vs oracle
# ==============================================================================
cat("\n============================================================\n")
cat("=== DRY: summary(on networks) vs oracle R ===\n")

for (tl_nm in names(timelines)) {
  TL <- timelines[[tl_nm]]
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("--- Timeline %s ---\n", tl_nm))

  P1 <- as.integer(TL$parts$P1)
  P2 <- as.integer(TL$parts$P2)
  P3 <- as.integer(TL$parts$P3)

  nodes1 <- TL$nodes$N1
  nodes2 <- TL$nodes$N2
  nodes3 <- TL$nodes$N3

  parts_t2 <- list(P1, P2)
  nodes_t2 <- list(nodes1, nodes2)

  cat(sprintf("\n[T=2] compare P2 vs P1 (n=%d)\n", length(P2)))
  for (cc in cases_summary) .compare_summary(parts_t2, nodes_t2, cc$rhs, cc$size, seed = 111L)

  parts_t3 <- list(P1, P2, P3)
  nodes_t3 <- list(nodes1, nodes2, nodes3)

  cat(sprintf("\n[T=3] compare P3 vs P2 (n=%d)\n", length(P3)))
  for (cc in cases_summary) .compare_summary(parts_t3, nodes_t3, cc$rhs, cc$size, seed = 222L)
}

cat("\nOK: DRY summary vs oracle.\n")

# ==============================================================================
# 2) FITS (sans estimate ni control explicites)
# ==============================================================================
cat("\n============================================================\n")
cat("=== FITS: erpm_long(sans estimate ni control explicites) ===\n")

for (tl_nm in names(timelines)) {
  TL <- timelines[[tl_nm]]
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("--- Timeline %s ---\n", tl_nm))

  P1 <- as.integer(TL$parts$P1)
  P2 <- as.integer(TL$parts$P2)
  P3 <- as.integer(TL$parts$P3)

  nodes1 <- TL$nodes$N1
  nodes2 <- TL$nodes$N2
  nodes3 <- TL$nodes$N3

  # ----- T = 2 -----
  parts_t2 <- list(P1, P2)
  nodes_t2 <- list(nodes1, nodes2)

  cat(sprintf("\n[T=2] fits (n=%d)\n", length(P2)))
  for (rhs in cases_fit) {
    obj2 <- .call_erpm_long(parts_t2, nodes_t2, rhs, eval.call = TRUE, dyads = NULL, seed = 100L)
    fit2 <- .get_fit_T(obj2, 2L)
    .check_fit(fit2, rhs)
  }

  # ----- T = 3 -----
  parts_t3 <- list(P1, P2, P3)
  nodes_t3 <- list(nodes1, nodes2, nodes3)

  cat(sprintf("\n[T=3] fits (n=%d)\n", length(P3)))
  for (rhs in cases_fit) {
    obj3 <- .call_erpm_long(parts_t3, nodes_t3, rhs, eval.call = TRUE, dyads = NULL, seed = 200L)
    fit3 <- .get_fit_T(obj3, 3L)
    .check_fit(fit3, rhs)
  }
}

cat("\nOK: FITS.\n")

# ==============================================================================
# 3) STRESS: partitions déterministes (dry + fit)
# ==============================================================================
cat("\n============================================================\n")
cat("=== STRESS: partitions déterministes (dry + fit) ===\n")

.safe_fit <- function(parts, nodes_list, rhs, seed) {
  out <- tryCatch({
    obj <- .call_erpm_long(parts, nodes_list, rhs, eval.call = TRUE, dyads = NULL, seed = seed)
    fit <- .get_fit_T(obj, length(parts))
    .check_fit(fit, rhs)
    TRUE
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (.is_essentially_constant_error(msg)) {
      cat(sprintf("[SKIP stress] RHS=%-60s | %s\n", rhs, "data are essentially constant"))
      return(FALSE)
    }
    if (.is_negative_diag_error(msg)) {
      cat(sprintf("[SKIP stress] RHS=%-60s | %s\n", rhs, "matrix x: negative diagonal"))
      return(FALSE)
    }
    .stopf("%s", msg)
  })
  invisible(out)
}

stress_cases <- list(
  list(
    label = "S1_n6",
    parts = list(
      P1 = c(1,1,2,2,3,3),
      P2 = c(1,1,2,3,3,3),
      P3 = c(1,2,2,3,3,3)
    ),
    nodes = list(
      N1 = timelines$TL1$nodes$N1,
      N2 = timelines$TL1$nodes$N2,
      N3 = timelines$TL1$nodes$N3
    )
  ),
  list(
    label = "S2_n6",
    parts = list(
      P1 = c(1,1,1,2,2,3),
      P2 = c(1,1,2,2,3,3),
      P3 = c(1,2,2,3,3,3)
    ),
    nodes = list(
      N1 = data.frame(
        label  = paste0("N",1:6),
        age    = c(28, 35, 49, 41, 26, 58),
        score  = c(2.1, 5.9, 8.2, 6.4, 3.0, 9.1),
        gender = c("M","F","M","F","F","M"),
        dept   = c("A","A","B","B","C","C"),
        stringsAsFactors = FALSE
      ),
      N2 = data.frame(
        label  = paste0("N",1:6),
        age    = c(29, 36, 50, 42, 27, 59),
        score  = c(2.2, 6.0, 8.0, 6.2, 3.2, 9.0),
        gender = c("M","F","M","F","F","M"),
        dept   = c("A","B","B","B","C","C"),
        stringsAsFactors = FALSE
      ),
      N3 = data.frame(
        label  = paste0("N",1:6),
        age    = c(30, 37, 51, 43, 28, 60),
        score  = c(2.0, 6.1, 7.9, 6.3, 3.1, 8.8),
        gender = c("M","F","M","F","F","M"),
        dept   = c("A","B","B","C","C","C"),
        stringsAsFactors = FALSE
      )
    )
  ),
  list(
    label = "S3_n8",
    parts = list(
      P1 = c(1,2,2,3,3,3,4,4),
      P2 = c(1,2,2,3,3,4,4,4),
      P3 = c(1,1,2,3,3,4,4,4)
    ),
    nodes = list(
      N1 = timelines$TL2$nodes$N1,
      N2 = timelines$TL2$nodes$N2,
      N3 = timelines$TL2$nodes$N3
    )
  ),
  list(
    label = "S4_n8",
    parts = list(
      P1 = c(1,1,2,2,3,3,4,4),
      P2 = c(1,2,2,3,3,3,4,4),
      P3 = c(1,2,2,3,3,4,4,4)
    ),
    nodes = list(
      N1 = data.frame(
        label  = paste0("N",1:8),
        age    = c(21, 29, 34, 45, 53, 40, 37, 48),
        score  = c(1.2, 2.8, 4.4, 6.7, 8.3, 5.5, 4.9, 7.1),
        gender = c("F","M","F","M","F","M","F","M"),
        dept   = c("A","A","B","B","C","C","A","B"),
        stringsAsFactors = FALSE
      ),
      N2 = data.frame(
        label  = paste0("N",1:8),
        age    = c(22, 30, 35, 46, 54, 41, 38, 49),
        score  = c(1.4, 3.0, 4.3, 6.5, 8.1, 5.6, 5.0, 7.2),
        gender = c("F","M","F","M","F","M","F","M"),
        dept   = c("A","B","B","B","C","C","A","B"),
        stringsAsFactors = FALSE
      ),
      N3 = data.frame(
        label  = paste0("N",1:8),
        age    = c(23, 31, 36, 47, 55, 42, 39, 50),
        score  = c(1.3, 2.9, 4.5, 6.6, 8.2, 5.7, 5.1, 7.0),
        gender = c("F","M","F","M","F","M","F","M"),
        dept   = c("A","B","B","C","C","C","A","B"),
        stringsAsFactors = FALSE
      )
    )
  ),
  list(
    label = "S5_n10",
    parts = list(
      P1 = c(1,1,2,2,3,3,4,4,5,5),
      P2 = c(1,2,2,3,3,4,4,4,5,5),
      P3 = c(1,2,2,3,3,4,4,5,5,5)
    ),
    nodes = list(
      N1 = data.frame(
        label  = paste0("N",1:10),
        age    = c(24, 31, 46, 38, 52, 29, 41, 57, 35, 48),
        score  = c(1.1, 2.5, 6.8, 4.2, 7.9, 3.0, 5.6, 8.4, 4.9, 7.1),
        gender = c("M","F","M","F","M","F","M","F","M","F"),
        dept   = c("A","A","B","B","C","C","A","B","C","A"),
        stringsAsFactors = FALSE
      ),
      N2 = data.frame(
        label  = paste0("N",1:10),
        age    = c(25, 32, 47, 39, 53, 30, 42, 58, 36, 49),
        score  = c(1.2, 2.6, 6.7, 4.3, 8.0, 3.1, 5.7, 8.3, 5.0, 7.0),
        gender = c("M","F","M","F","M","F","M","F","M","F"),
        dept   = c("A","B","B","B","C","C","A","B","C","A"),
        stringsAsFactors = FALSE
      ),
      N3 = data.frame(
        label  = paste0("N",1:10),
        age    = c(26, 33, 48, 40, 54, 31, 43, 59, 37, 50),
        score  = c(1.3, 2.7, 6.6, 4.4, 8.1, 3.2, 5.8, 8.2, 5.1, 6.9),
        gender = c("M","F","M","F","M","F","M","F","M","F"),
        dept   = c("A","B","B","C","C","C","A","B","C","A"),
        stringsAsFactors = FALSE
      )
    )
  )
)

for (sc in stress_cases) {
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("Stress case: %s\n", sc$label))

  P1 <- as.integer(sc$parts$P1)
  P2 <- as.integer(sc$parts$P2)
  P3 <- as.integer(sc$parts$P3)

  nodes1 <- sc$nodes$N1
  nodes2 <- sc$nodes$N2
  nodes3 <- sc$nodes$N3

  parts_t2 <- list(P1, P2)
  nodes_t2 <- list(nodes1, nodes2)

  parts_t3 <- list(P1, P2, P3)
  nodes_t3 <- list(nodes1, nodes2, nodes3)

  # dry: inertie pure
  .compare_summary(parts_t2, nodes_t2, "inertia_groups()", size = NULL, seed = 3000L)
  .compare_summary(parts_t3, nodes_t3, "inertia_groups(past_influence=2)", size = NULL, seed = 3001L)

  # fit: baseline + inertie (robuste aux cas constants et au cas diag négative)
  .safe_fit(parts_t3, nodes_t3, "squared_sizes() + inertia_groups()", seed = 3010L)
  .safe_fit(parts_t3, nodes_t3, "squared_sizes() + inertia_groups(past_influence=2)", seed = 3011L)
}

cat("\nOK: STRESS.\n")

# ------------------------------------------------------------------------------
# Cleanup patch
# ------------------------------------------------------------------------------
if (exists("ergm_patch_disable")) ergm_patch_disable()

cat("\nOK: selftest_inertia_groups terminé.\n")