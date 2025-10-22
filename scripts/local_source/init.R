# ==============================================================================
# Fichier : init.R
# Fonction : init_erpm
# Utilité : Initialiser l'environnement ERPM, charger les settings, packages et le patch ERGM
# ==============================================================================

if(!exists(".__init_loaded", envir = .GlobalEnv)){

  #' Initialisation de l'environnement ERPM
  #'
  #' Cette fonction prépare l'environnement de travail pour le projet ERPM :
  #' - Charge les settings globaux.
  #' - Vérifie et installe les packages nécessaires listés dans DESCRIPTION.
  #' - Charge automatiquement les modules via \pkg{devtools} si disponible.
  #' - Ajoute le patch pour corriger le bug \code{type 'builtin' d'indice incorrect} dans \pkg{ergm}.
  #' - Exécute un self-test du patch si \code{selftest = TRUE}.
  #'
  #' @param settings_file Chemin vers le fichier de settings (par défaut \code{"scripts/local_source/settings.R"}).
  #' @param ergm_patch_file Chemin vers le fichier contenant le patch ERGM.
  #' @param description_file Fichier DESCRIPTION contenant les dépendances du projet.
  #' @param verbose Logique : afficher les messages d'information ou non (par défaut TRUE).
  #' @param selftest Logique : exécuter le self-test du patch ERGM après initialisation.
  #'
  #' @return \code{TRUE} invisiblement si l'initialisation a réussi.
  #' @examples
  #' \dontrun{
  #' init_erpm(verbose = TRUE, selftest = FALSE)
  #' }
  #' @export
  init_erpm <- function(
    ergm_patch_file = "scripts/ergm_patch.R",
    description_file = "DESCRIPTION",
    verbose = TRUE, 
    selftest = FALSE
  ) {
    # Garde-fou : normalement settings est défini avant le call de init_erpm
    if (!exists(".__settings_loaded", envir = .GlobalEnv)) {
      source("scripts/local_source/settings.R", local = FALSE)
    }

    # ==== Vérifie que toutes les variables sont bien là ====
    required_globals <- c("VERBOSE", "DEV_MODE", "PROJECT_ROOT")
    missing_globals <- required_globals[!sapply(required_globals, exists, envir = .GlobalEnv)]
    if (length(missing_globals) > 0) {
      stop(paste("Les variables globales suivantes ne sont pas définies dans", settings_file, ":", 
                 paste(missing_globals, collapse = ", ")))
    }

    # ==== Récupère les packages depuis DESCRIPTION ====
    deps_info <- export_dependencies_from_description_file(description_file)

    # ==== Version de R ====
    if (!is.null(deps_info$r_version) && getRversion() < deps_info$r_version) {
      stop(sprintf("R >= %s requis par %s", deps_info$r_version, description_file))
    }

    # ==== Modules à charger ====
    base_pkgs <- unique(c(deps_info$packages, "parallel", "igraph", "RColorBrewer", "devtools"))
    installation_status <- install_and_import_if_missing(base_pkgs, TRUE)

    # ==== Vérifie si 'devtools' a été chargé avec succès ====
    devtools_index <- which(installation_status$packages == "devtools")
    if (length(devtools_index) == 0 || !installation_status$packages_loaded[devtools_index]) {
      msg <- "'devtools' n'est pas disponible — chargement manuel des fichiers R."
      if (exists("log_msg")) log_msg("WARN", msg) else warning(msg)
      r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
      invisible(sapply(r_files, source, local = .GlobalEnv))
    } else {
      if (exists("log_msg")) log_msg("INFO", "'devtools' détecté — utilisation de devtools::load_all().")
      else message("'devtools' détecté — utilisation de devtools::load_all().")
      devtools::load_all()
    }

    # ==== Ajout du patch ergm ====
    source(ergm_patch_file, local = TRUE)

    # ==== Self-test du patch si demandé ====
    if(selftest){
      res_tests <- ergm_patch_selftest(run_diagnostics = FALSE, verbose = verbose)
    }

    invisible(TRUE)
  }

  # --- Attribution à l'environnement global ---
  assign("init_erpm", init_erpm, envir = .GlobalEnv)
  assign(".__init_loaded", TRUE, envir = .GlobalEnv)
}
