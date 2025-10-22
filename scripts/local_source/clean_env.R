# ==============================================================================
# Fichier : clean_env.R
# Fonction :
# Utilité : 
# ==============================================================================

if(!exists(".__clean_env_loaded", envir = .GlobalEnv)){

    # ==============================================================================
    #' Nettoyage de l'environnement global
    #'
    #' Cette fonction supprime toutes les variables globales sauf celles de la
    #' "liste blanche", désactive le patch ERGM et libère la mémoire.
    #'
    #' @param keep_vars Optionnel : vecteur de noms de variables à conserver.
    #' @param verbose Logique : afficher un message de confirmation.
    #' @export
    clean_global_env <- function(keep_vars = NULL, verbose = TRUE) {
        if (is.null(keep_vars)) {
            if (exists("DEV_MODE", envir = .GlobalEnv) && DEV_MODE) {
            keep_vars <- c()
            } else {
            keep_vars <- c(
                "VERBOSE", "DEV_MODE", "PROJECT_ROOT", "CRAN_MIRROR", "INC_DIR",
                "log_msg", "red", "green", "yellow", "blue", "bold"
            )
            }
        }
        
        # Supprime le reste
        vars_to_remove <- setdiff(ls(envir = .GlobalEnv, all.names = TRUE), keep_vars)
        if (length(vars_to_remove) > 0) rm(list = vars_to_remove, envir = .GlobalEnv)
        
        # Nettoyage mémoire
        gc()
        if (isTRUE(VERBOSE)) cat("\n✔ Environnement nettoyé avec succès.\n")
    }

    assign("clean_global_env",      clean_global_env,   envir = .GlobalEnv)
    assign(".__clean_env_loaded",   TRUE,               envir = .GlobalEnv)
}
