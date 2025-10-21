# ==============================================================================
# Fichier : settings.R
# Fonction : 
# Utilité : Variables et fonctions globales du projet ERPM
# ==============================================================================
if (!exists(".__settings_loaded", envir = .GlobalEnv)) {

    
    # --- Détection automatique du répertoire racine du projet ---
    # (fonctionne même si on fait source("scripts/cub_main.R"))
    this_file <- tryCatch({
        # cas 1 : exécuté via source()
        if (!is.null(sys.calls()) && length(sys.calls()) > 0) {
        src <- sys.calls()
        call <- as.character(src[[length(src)]])
        # Récupère le chemin complet si possible
        if (grepl("source", call)) {
            normalizePath(sub('.*source\\(\"?([^\"]+)\"?\\).*', "\\1", call))
        } else {
            NA
        }
        } else NA
    }, error = function(e) NA)

    if (!is.na(this_file)) {
        PROJECT_ROOT <- dirname(dirname(this_file))
    } else {
        PROJECT_ROOT <- getwd()
    }

    Sys.setlocale("LC_CTYPE", "UTF-8")

    # Variables globales
    if (!exists("VERBOSE"))         assign("VERBOSE",       TRUE,                           envir = .GlobalEnv)
    if (!exists("DEV_MODE"))        assign("DEV_MODE",      FALSE,                          envir = .GlobalEnv)
    if (!exists("PROJECT_ROOT"))    assign("PROJECT_ROOT",  getwd(),                        envir = .GlobalEnv)
    if (!exists("CRAN_MIRROR"))     assign("CRAN_MIRROR",   "https://cloud.r-project.org",  envir = .GlobalEnv)
    if (!exists("INC_DIR"))         assign("INC_DIR",       "scripts/local_source/",        envir = .GlobalEnv)

    # Fonctions/Variables Globales
    source(sprintf("%scolors.R", INC_DIR))
    source(sprintf("%slogging.R", INC_DIR))
    source(sprintf("%simport_dependencies.R", INC_DIR))

    # --- Marqueur interne pour indiquer que le module est chargé ---
    assign(".__settings_loaded", TRUE, envir = .GlobalEnv)
}