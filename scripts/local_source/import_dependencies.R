# ==============================================================================
# Fichier : import_dependencies.R
# Fonction : install_if_missing()
# Utilité : Installe et charge les bibliothèques R nécessaires, avec log et couleur
# ==============================================================================
if(!exists(".__import_dependencies_loaded", envir = .GlobalEnv)){

    # Fonction
    install_and_import_if_missing <- function(pkgs, stop_rule=TRUE, log_dir = "log") {
    
        # Charge settings si pas déjà présent
        if (!exists(".__settings_loaded", envir = .GlobalEnv)) {   
            tryCatch(source(sprintf("%ssettings.R", INC_DIR)), error=function(e) VERBOSE <<- TRUE)
        }

        # Charge couleurs si absentes
        if (!exists(".__colors_loaded", envir = .GlobalEnv)) {
            tryCatch(source(sprintf("%scolors.R", INC_DIR)), error=function(e) {
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
                    else cat("SUCCESS", paste("Installation réussie pour", p))

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
                #cat(col(sprintf("%s : %s\n", pkgs[i], if (p_loaded[i]) "OK" else "Échec")))
            }
        }
        
        return(list(packages=pkgs, packages_loaded=p_loaded))
    } # install_and_import_if_missing

    # Fonction
    export_dependencies_from_description_file <- function(file="DESCRIPTION", optional=FALSE){
        desc <- read.dcf(file)

        # Depends et Imports
        depends <- unlist(strsplit(desc[1, "Depends"], ","))
        imports <- unlist(strsplit(desc[1, "Imports"], ","))

        pkgs <- trimws(c(depends, imports))
        
        # Extraire la version R si mentionnée
        r_dep <- pkgs[grepl("^R\\s*\\(>=\\s*[0-9.]+\\)$", pkgs)]
        r_version <- if(length(r_dep)) gsub("^R\\s*\\(>=\\s*([0-9.]+)\\)$", "\\1", r_dep) else NULL
        pkgs <- pkgs[!grepl("^R", pkgs)]  # enlever R du vecteur de packages

        # Extraire les versions des packages (ex: "ergm (>= 4.8.0)")
        pkg_versions <- sapply(pkgs, function(x){
            m <- regmatches(x, regexec("^([a-zA-Z0-9.]+)\\s*\\(>=\\s*([0-9.]+)\\)$", x))[[1]]
            if(length(m) == 3) return(list(name=m[2], version=m[3]))
            else return(list(name=x, version=NULL))
        }, simplify = FALSE)

        # Ne garder que les noms pour installer
        base_pkgs <- sapply(pkg_versions, `[[`, "name")

        return(list(packages=base_pkgs, pkg_versions=pkg_versions, r_version=r_version))
    }

    assign("install_and_import_if_missing",                 install_and_import_if_missing,              envir = .GlobalEnv)
    assign("export_dependencies_from_description_file",     export_dependencies_from_description_file,  envir = .GlobalEnv)
    assign(".__import_dependencies_loaded",                 TRUE,                                       envir = .GlobalEnv)
}
