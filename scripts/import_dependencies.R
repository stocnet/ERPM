# ==============================================================================
# Fichier : import_dependencies.R
# Fonction : install_if_missing()
# Utilité : Installe et charge les bibliothèques R nécessaires, avec log et couleur
# ==============================================================================
if(!exists(".__install_and_import_if_missing", envir = .GlobalEnv)){
     #assign("install_if_missing", function(pkgs, stop_rule=TRUE, log_dir="log") {
     install_and_import_if_missing <- function(pkgs, stop_rule=TRUE, log_dir = "log") {
        
        # Charge settings si pas déjà présent
        if (!exists(".__settings_loaded", envir = .GlobalEnv)) {
            tryCatch(source("scripts/settings.R"), error=function(e) VERBOSE <<- TRUE)
        }

        # Charge couleurs si absentes
        if (!exists(".__colors_loaded", envir = .GlobalEnv)) {
            tryCatch(source("scripts/colors.R"), error=function(e) {
                green  <<- function(x) paste0("\033[32m", x, "\033[0m")
                red    <<- function(x) paste0("\033[31m", x, "\033[0m")
                yellow <<- function(x) paste0("\033[33m", x, "\033[0m")
            })
        }

        p_loaded <- logical(length(pkgs))  # Initialise un vecteur logique

        # On charge toutes les bibliothèques demandées depuis la source standard
        for (i in seq_along(pkgs)) {
            p <- pkgs[i]

            # Vérifie si le package est installé, sinon on l'installe
            if (!requireNamespace(p, quietly = TRUE)) {
                message(sprintf("Installation du module manquant : %s", p))

                tryCatch({
                    install.packages(p, repos=CRAN_MIRROR)
                    if (exists("log_msg")) log_msg("SUCCESS", paste("Installation réussie pour", p))
                    else msg("SUCCESS", paste("Installation réussie pour", p))

                }, error=function(e) {
                    if (exists("log_msg")) log_msg("ERROR", paste("Échec installation de", p, ":", e$message))
                    else cat("ERROR", paste("Échec installation de", p, ":", e$message))
                
                })
            }

            # Charge le package après installation
            if (!suppressPackageStartupMessages(require(p, character.only = TRUE))) {
                p_loaded[i] <- FALSE
                if (stop_rule) stop(paste("Échec du chargement de", p))

            } else {
                p_loaded[i] <- TRUE
                if (exists("log_msg")) log_msg("SUCCESS", paste("Module", p, "chargé"))
                else cat("SUCCESS", paste("Module", p, "chargé"))

            }
        }

        # Résumé visuel si VERBOSE = TRUE
        if (exists("VERBOSE") && VERBOSE) {

            for (i in seq_along(pkgs)) {
                col <- if (p_loaded[i]) green else red
                cat(col(sprintf(" - %s : %s\n", pkgs[i], if (p_loaded[i]) "OK" else "Échec")))
            }
        }
        
        invisible(p_loaded)
    }#, envir = .GlobalEnv))

    assign(".__install_and_import_if_missing", TRUE, envir = .GlobalEnv)

}