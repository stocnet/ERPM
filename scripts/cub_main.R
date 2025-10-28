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

# ret <- launch_model(
#   engine      = "ergm",
#   effects     = "b2degrange(from=2,to=3) + squared_sizes(from=1,to=2)", #b2degrange #group
#   dry_run     = FALSE,
#   estimate    = "CD",#"MLE",#
#   eval_loglik = FALSE,         # pas de bridge sampling
#   control     = list(MCMLE.maxit = 3, MCMC.samplesize = 1000), # override léger
#   timeout     = 60             # coupe au bout de 60s
# )

ret <- launch_model(
  engine      = "summary",#"ergm",#"erpm",#
  effects     = "groups(from=2,to=3)", #b2degrange #group
  dry_run     = FALSE,#TRUE,#
  estimate    = "CD",#"MLE",#
  eval_loglik = FALSE,         # pas de bridge sampling
  control     = list(MCMLE.maxit = 3, MCMC.samplesize = 1000), # override léger
  timeout     = 60             # coupe au bout de 60s
)

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