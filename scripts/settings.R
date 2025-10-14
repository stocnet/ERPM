# ==============================================================================
# Fichier : settings.R
# Fonction : 
# Utilité : Variables et fonctions globales du projet ERPM
# ==============================================================================
if (!exists(".__settings_loaded", envir = .GlobalEnv)) {

    Sys.setlocale("LC_CTYPE", "UTF-8")

    # Variables globales
    # VERBOSE <- TRUE             # Affiche des infos facultatives
    # DEBUG_MODE <- FALSE         # Affiche des infos de debug
    # PROJECT_ROOT <- getwd()     # Récupère le chemin du répertoire courant
    if (!exists("VERBOSE"))         assign("VERBOSE",       TRUE,                           envir = .GlobalEnv)
    if (!exists("DEBUG_MODE"))      assign("DEBUG_MODE",    FALSE,                          envir = .GlobalEnv)
    if (!exists("PROJECT_ROOT"))    assign("PROJECT_ROOT",  getwd(),                        envir = .GlobalEnv)
    if (!exists("CRAN_MIRROR"))     assign("CRAN_MIRROR",  "https://cloud.r-project.org",   envir = .GlobalEnv)

    # Fonctions/Variables Globales
    source("scripts/colors.R")
    source("scripts/logging.R")
    source("scripts/import_dependencies.R")

    # --- Marqueur interne pour indiquer que le module est chargé ---
    assign(".__settings_loaded", TRUE, envir = .GlobalEnv)
}