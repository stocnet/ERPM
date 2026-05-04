# ==============================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_inertia_groups.R
# Objet   : MWE pour l'effet inertiel ERPM `inertia_groups` (PLE empilé)
#
# L'effet inertia_groups compte les acteurs dont le groupe au temps t est
# EXACTEMENT reproduit (même ensemble de membres) dans au moins un des d
# pas de temps précédents — sémantique d'UNION sur les lags 1..d.
#
#   T_inertia(y, d) = #{i, t>d : groupe_t(i) = groupe_{t-l}(i) pour au moins un l in 1..d}
#
# Chaîne ERPM PLE :
#   partitions (liste T) -> erpm_long(mode="empile") ->
#   réseau bipartite méta -> summary(méta ~ ...) / fit(erpm_long(...))
#
# Structure :
#   EXEMPLE 1 — données codées en dur (T=3, n=4)
#     - vérification analytique pour d=1 et d=2
#     - summary via erpm_long dry-run + stopifnot
#     - fit erpm_long + print(summary)
#   EXEMPLE 2 — données aléatoires (T=3, n=6, K=3, set.seed=42)
#     - mêmes vérifications pour d=1
#     - fit erpm_long + print(summary)
# ==============================================================================

options(ergm.loglik.warn_dyads = FALSE)
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  library(devtools)
  library(ergm)
})

devtools::load_all(".")

if (file.exists("scripts/ergm_patch.R")) {
  source("scripts/ergm_patch.R")
  ergm_patch_enable()
}

# ==============================================================================
# Fonction de référence analytique
# ==============================================================================

# Compte les acteurs dont le groupe au temps t (t > d) est identique à leur
# groupe en au moins un des d pas précédents (UNION des lags 1..d).
# "Identique" signifie : même ensemble de membres exactement.
ref_inertia_groups <- function(partitions, d = 1L) {
  T <- length(partitions)
  d <- as.integer(d)
  cnt <- 0L
  for (t in seq(d + 1L, T)) {
    cur <- partitions[[t]]
    n_t <- length(cur)
    for (i in seq_len(n_t)) {
      g_cur <- sort(which(cur == cur[i]))
      for (lag in seq_len(d)) {
        prev <- partitions[[t - lag]]
        if (length(prev) != n_t) next
        if (identical(g_cur, sort(which(prev == prev[i])))) {
          cnt <- cnt + 1L
          break
        }
      }
    }
  }
  cnt
}

# ==============================================================================
# EXEMPLE 1 — Données codées en dur
# ==============================================================================
cat("\n")
cat("================================================================================\n")
cat("EXEMPLE 1 — Données codées en dur (T=3, n=4)\n")
cat("================================================================================\n")

# Partitions :
#   t=1 : {1,2} | {3,4}    (2 groupes de taille 2)
#   t=2 : {1,2} | {3} | {4} (groupe {1,2} conservé, {3,4} scindé)
#   t=3 : {1,2} | {3,4}    (retour à t=1)
#
# Valeurs attendues :
#   d=1 : t=2 → acteurs 1,2 inertiels (+2) ; t=3 → acteurs 1,2 (+2)      = 4
#   d=2 : t=3 → acteurs 1,2 (lag 1) + acteurs 3,4 (lag 2, groupe {3,4}   = 4
#               identique à t=1) [union sur lags : chaque acteur compté 1 fois]
partitions1 <- list(
  c(1, 1, 2, 2),
  c(1, 1, 2, 3),
  c(1, 1, 2, 2)
)

cat("\nPartitions :\n")
for (t in seq_along(partitions1)) {
  sz <- as.integer(table(partitions1[[t]]))
  cat(sprintf("  t=%d : %s  (tailles : %s)\n",
              t, paste(partitions1[[t]], collapse=" "), paste(sz, collapse=", ")))
}

ref1_d1 <- ref_inertia_groups(partitions1, d = 1L)
ref1_d2 <- ref_inertia_groups(partitions1, d = 2L)
cat(sprintf("\nRéférences analytiques :\n"))
cat(sprintf("  inertia_groups(d=1) = %g\n", ref1_d1))
cat(sprintf("  inertia_groups(d=2) = %g\n", ref1_d2))

# ---- Summary via erpm_long (dry-run) ----
cat("\n--- Summary via erpm_long(eval.call=FALSE) ---\n")

out1_d1 <- erpm_long(
  partitions1 ~ inertia_groups(past_influence = 1L),
  eval.call   = FALSE,
  verbose     = FALSE,
  constraints = ~ b1partblockdiag("timeblock")
)
meta_nw1 <- attr(out1_d1, "meta_nw")

obs1_d1 <- as.numeric(summary(
  meta_nw1 ~ inertia_groups(past_influence = 1L),
  constraints = ~ b1partblockdiag("timeblock")
))
cat(sprintf("  inertia_groups(d=1) : obs=%-6g  ref=%-6g  %s\n",
            obs1_d1, ref1_d1,
            if (isTRUE(all.equal(obs1_d1, ref1_d1, tol = 0))) "OK" else "*** MISMATCH ***"))
stopifnot(isTRUE(all.equal(obs1_d1, ref1_d1, tol = 0)))

out1_d2 <- erpm_long(
  partitions1 ~ inertia_groups(past_influence = 2L),
  eval.call   = FALSE,
  verbose     = FALSE,
  constraints = ~ b1partblockdiag("timeblock")
)
obs1_d2 <- as.numeric(summary(
  attr(out1_d2, "meta_nw") ~ inertia_groups(past_influence = 2L),
  constraints = ~ b1partblockdiag("timeblock")
))
cat(sprintf("  inertia_groups(d=2) : obs=%-6g  ref=%-6g  %s\n",
            obs1_d2, ref1_d2,
            if (isTRUE(all.equal(obs1_d2, ref1_d2, tol = 0))) "OK" else "*** MISMATCH ***"))
stopifnot(isTRUE(all.equal(obs1_d2, ref1_d2, tol = 0)))

# ---- Fit ERPM ----
# Note : eval.loglik=FALSE contourne un bug dans ergm 4.10.x (replace() appelé
# avec un builtin) qui affecte uniquement le calcul de la log-vraisemblance.
# Les coefficients estimés ne sont pas affectés.
cat("\n--- Fit : erpm_long(partitions1 ~ inertia_groups(d=1), eval.loglik=FALSE) ---\n")
erpm_call1 <- erpm_long(
  partitions1 ~ inertia_groups(past_influence = 1L),
  eval.call   = TRUE,
  eval.loglik = FALSE,
  verbose     = FALSE,
  constraints = ~ b1partblockdiag("timeblock")
)
set.seed(1)
fit1 <- eval(erpm_call1)
cat(sprintf("  Coefficient inertia_groups(d=1) : %.4f\n", coef(fit1)))
cat("\n")
print(summary(fit1))

# ==============================================================================
# EXEMPLE 2 — Données aléatoires
# ==============================================================================
cat("\n")
cat("================================================================================\n")
cat("EXEMPLE 2 — Données aléatoires (T=3, n=6, K=3, seed=42)\n")
cat("================================================================================\n")

set.seed(42)
n2 <- 6L
T2 <- 3L
K2 <- 3L
partitions2 <- lapply(seq_len(T2), function(t) sample(seq_len(K2), n2, replace = TRUE))

# Garantir que chaque groupe est représenté à chaque pas
for (t in seq_len(T2)) {
  missing <- setdiff(seq_len(K2), partitions2[[t]])
  for (k in missing) partitions2[[t]][sample(n2, 1L)] <- k
}
partitions2 <- lapply(partitions2, function(p) as.integer(factor(p)))

cat("\nPartitions :\n")
for (t in seq_along(partitions2)) {
  sz <- as.integer(table(partitions2[[t]]))
  cat(sprintf("  t=%d : %s  (tailles : %s)\n",
              t, paste(partitions2[[t]], collapse=" "), paste(sz, collapse=", ")))
}

ref2_d1 <- ref_inertia_groups(partitions2, d = 1L)
cat(sprintf("\nRéférence analytique :\n"))
cat(sprintf("  inertia_groups(d=1) = %g\n", ref2_d1))

# ---- Summary via erpm_long (dry-run) ----
cat("\n--- Summary via erpm_long(eval.call=FALSE) ---\n")

out2_d1 <- erpm_long(
  partitions2 ~ inertia_groups(past_influence = 1L),
  eval.call   = FALSE,
  verbose     = FALSE,
  constraints = ~ b1partblockdiag("timeblock")
)
obs2_d1 <- as.numeric(summary(
  attr(out2_d1, "meta_nw") ~ inertia_groups(past_influence = 1L),
  constraints = ~ b1partblockdiag("timeblock")
))
cat(sprintf("  inertia_groups(d=1) : obs=%-6g  ref=%-6g  %s\n",
            obs2_d1, ref2_d1,
            if (isTRUE(all.equal(obs2_d1, ref2_d1, tol = 0))) "OK" else "*** MISMATCH ***"))
stopifnot(isTRUE(all.equal(obs2_d1, ref2_d1, tol = 0)))

# ---- Fit ERPM ----
cat("\n--- Fit : erpm_long(partitions2 ~ inertia_groups(d=1), eval.loglik=FALSE) ---\n")
erpm_call2 <- erpm_long(
  partitions2 ~ inertia_groups(past_influence = 1L),
  eval.call   = TRUE,
  eval.loglik = FALSE,
  verbose     = FALSE,
  constraints = ~ b1partblockdiag("timeblock")
)
set.seed(42)
fit2 <- eval(erpm_call2)
cat(sprintf("  Coefficient inertia_groups(d=1) : %.4f\n", coef(fit2)))
cat("\n")
print(summary(fit2))

cat("\n=== [DONE] MWE inertia_groups ===\n")

if (exists("ergm_patch_disable", mode = "function")) {
  try(ergm_patch_disable(), silent = TRUE)
}
