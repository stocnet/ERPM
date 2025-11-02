# ====================================================================================== 
# Fichier : cub_test.R
# Fonction :
# Utilité : test pour les nouvelles fonctions ERPM
# ====================================================================================== 

# macOS/Linux
Sys.setenv(LANG="fr_FR.UTF-8")
invisible(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"))

# ====================================================================================== 
# ======================================== INIT ======================================== 
suppressMessages(library(ergm))

if (!exists(".__debug_loaded", envir = .GlobalEnv)) {
  source("scripts/local_source/debug.R", local = FALSE)
}

source("scripts/local_source/init.R", local = FALSE) # init l'environnement (source de settings et du launcher)
source("R/erpm_wrapper.R", local = FALSE) # n'est pas dans scripts/local_source/ donc on le charge indépendemment
# debug_source_file("R/erpm_wrapper.R") 


init_erpm(selftest=FALSE, verbose=FALSE)

if (exists("log_msg")) {
  log_msg("INFO", "Démarrage du script pour ERPM")
} else {
  message("Démarrage du script pour ERPM")
}

if (exists("ergm_patch_enable")) ergm_patch_enable( verbose = VERBOSE )


# ====================================================================================== 
# ======================================== RUN ========================================= 

# =========================
# Partitions de test
# =========================
partitions <- list(
  list(name = "P1", part = c(1,1,2)),
  list(name = "P2", part = c(1,1,1)),
  list(name = "P3", part = c(1,1,1,1,2,2,3)),
  list(name = "P4", part = c(1,1,1,1,2,2,2)),
  list(name = "P5", part = c(1,1,1,1,2,2,2,2,3)),
  list(name = "P6", part = c(1,1,1,1,2,3,4,5,5))
)

# =========================
# Effets à tester (ton format exact)
# =========================
cases <- list(
  # Acceptés
  list(name = "clq_default",         effects = "cliques"),
  list(name = "clq_k2_F",            effects = "cliques(clique_size = 2, normalized = FALSE)"),
  list(name = "clq_k2_T",            effects = "cliques(clique_size = 2, normalized = TRUE)"),
  list(name = "clq_normF_defaultk",  effects = "cliques(normalized = FALSE)"),
  list(name = "clq_normT_defaultk",  effects = "cliques(normalized = TRUE)"),
  list(name = "clq_k2_only",         effects = "cliques(clique_size = 2)")
  # NB: les cas k=1 sont déplacés dans err_cases plus bas
)

results <- list()

launch_model(engine="erpm", effects="cliques",                         dry_run=TRUE)
launch_model(engine="erpm", effects="cliques()",                       dry_run=TRUE)
launch_model(engine="erpm", effects="cliques(normalized=FALSE)",       dry_run=TRUE)
launch_model(engine="erpm", effects="cliques(normalized=TRUE)",        dry_run=TRUE)
launch_model(engine="erpm", effects="cliques(clique_size=3)",          dry_run=TRUE)
launch_model(engine="erpm", effects="cliques(clique_size=1)",          dry_run=FALSE)  

# =========================
# Boucle principale : pour chaque partition × chaque effet
# =========================
for (pp in partitions) {
  if (exists("log_msg")) log_msg("INFO", sprintf("==== Partition %s : %s ====", pp$name, paste(pp$part, collapse=","))) else
    message(sprintf("==== Partition %s : %s ====", pp$name, paste(pp$part, collapse=",")))

  for (cx in cases) {
    # 1) SUMMARY (stat ERGM déterministe)
    if (exists("log_msg")) log_msg("INFO", sprintf("[summary] %s | %s", pp$name, cx$name))
    sres <- try(
      launch_model(
        engine      = "summary",
        effects     = cx$effects,
        partition   = pp$part,
        dry_run     = FALSE,
        verbose     = FALSE,     # évite le bruit en console, logs gardent la trace
        constraints = ~ b1part
      ),
      silent = TRUE
    )
    results[[sprintf("%s__%s__summary", pp$name, cx$name)]] <- sres

    # 2) ERGM (fit via wrapper ERPM → ERGM), comme dans ton format
    if (exists("log_msg")) log_msg("INFO", sprintf("[ergm/CD] %s | %s", pp$name, cx$name))
    fres <- try(
      launch_model(
        engine      = "erpm",
        effects     = cx$effects,
        partition   = pp$part,
        dry_run     = FALSE,
        estimate    = "CD",
        eval_loglik = TRUE,
        control     = list(MCMLE.maxit = 3, MCMC.samplesize = 1000),
        timeout     = NULL,
        constraints = ~ b1part
      ),
      silent = TRUE
    )
    results[[sprintf("%s__%s__fit", pp$name, cx$name)]] <- fres
  }
}

# =========================
# Cas erreurs attendus (k < 2)
# =========================
err_cases <- list(
  list(name = "clq_k1_default", effects = "cliques(clique_size = 1)"),
  list(name = "clq_k1_F",       effects = "cliques(clique_size = 1, normalized = FALSE)")
)

for (pp in partitions) {
  for (ex in err_cases) {
    if (exists("log_msg")) log_msg("INFO", sprintf("[ERROR EXPECTED] %s | %s", pp$name, ex$name))
    res <- try(
      launch_model(
        engine      = "summary",
        effects     = ex$effects,
        partition   = pp$part,
        dry_run     = FALSE,
        verbose     = FALSE,
        constraints = ~ b1part
      ),
      silent = TRUE
    )
    results[[sprintf("%s__%s__summary", pp$name, ex$name)]] <- res

    res2 <- try(
      launch_model(
        engine      = "erpm",
        effects     = ex$effects,
        partition   = pp$part,
        dry_run     = FALSE,
        estimate    = "CD",
        eval_loglik = TRUE,
        control     = list(MCMLE.maxit = 3, MCMC.samplesize = 1000),
        timeout     = NULL,
        constraints = ~ b1part
      ),
      silent = TRUE
    )
    results[[sprintf("%s__%s__fit", pp$name, ex$name)]] <- res2
  }
}

# cases <- list(
#   # Acceptés
#   list(name = "clq_default",         effects = "cliques"),                                       # -> edges()
#   list(name = "clq_k2_F",            effects = "cliques(clique_size = 2, normalized = FALSE)"),  # -> edges()
#   list(name = "clq_k2_T",            effects = "cliques(clique_size = 2, normalized = TRUE)"),   # -> edges() + normalisation post
#   list(name = "clq_normF_defaultk",  effects = "cliques(normalized = FALSE)"),                   # -> edges()
#   list(name = "clq_normT_defaultk",  effects = "cliques(normalized = TRUE)"),                    # -> edges() + normalisation post
#   list(name = "clq_k2_only",         effects = "cliques(clique_size = 2)"),                      # -> edges()
#   list(name = "clq_k1_default",      effects = "cliques(clique_size = 1)"),                      # -> b2degrange(0,0)
#   list(name = "clq_k1_F",            effects = "cliques(clique_size = 1, normalized = FALSE)")   # -> b2degrange(0,0)
# )

# results <- list()

# for (cx in cases) {
#   log_msg("INFO", sprintf("ERPM run: %s -> %s", cx$name, cx$effects))
#   results[[cx$name]] <- launch_model(
#     engine      = "erpm",
#     effects     = cx$effects,
#     dry_run     = FALSE,
#     estimate    = "CD",
#     eval_loglik = TRUE,
#     control     = list(MCMLE.maxit = 3, MCMC.samplesize = 1000),
#     timeout     = NULL
#   )
# }

# ======================================================================================= 
# ======================================== CLEAN ======================================== 

# Désactive le patch ERGM si actif
if (exists("ergm_patch_disable")) ergm_patch_disable( verbose = VERBOSE )

# Vide le buffer stdout
flush.console()  

# Log
log_msg("INFO", "Fin du programme -- Nettoyage de l'environnement global")

# Nettoie l'environnement
clean_global_env(verbose = VERBOSE)