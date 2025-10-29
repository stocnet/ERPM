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
if (!exists(".__debug_loaded", envir = .GlobalEnv)) {
  source("scripts/local_source/debug.R", local = FALSE)
}

source("scripts/local_source/init.R", local = FALSE) # init l'environnement (source de settings et du launcher)
source("R/erpm_wrapper.R", local = FALSE) # n'est pas dans scripts/local_source/ donc on le charge indépendemment
# debug_source_file("R/erpm_wrapper.R") # n'est pas dans scripts/local_source/ donc on le charge indépendemment


init_erpm(selftest=FALSE, verbose=FALSE)

if (exists("log_msg")) {
  log_msg("INFO", "Démarrage du script pour ERPM")
} else {
  message("Démarrage du script pour ERPM")
}

if (exists("ergm_patch_enable")) ergm_patch_enable( verbose = VERBOSE )


# ====================================================================================== 
# ======================================== RUN ========================================= 


cases <- list(
  list(name="groups_all",        effects="groups"),                  # ≈ b2degrange(1, Inf)
  list(name="groups_exact_3",    effects="groups(3)"),               # ≈ b2degrange(3, 4)
  list(name="groups_interval",   effects="groups(from=2,to=4)")      # ≈ b2degrange(2, 4)
)

results <- list()

for (cx in cases) {
  log_msg("INFO", sprintf("ERPM run: %s -> %s", cx$name, cx$effects))
  results[[cx$name]] <- launch_model(
    engine      = "erpm",
    effects     = cx$effects,
    dry_run     = FALSE,
    estimate    = "CD",
    eval_loglik = TRUE,
    control     = list(MCMLE.maxit = 3, MCMC.samplesize = 1000),
    timeout     = NULL
  )
}


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