# ==============================================================================
# Fichier : cub_test.R
# Fonction :
# Utilité : test pour les nouvelles fonctions ERPM
# ==============================================================================

# Initialise l'environnement du projet
# if (!exists(".ERPMenv", envir = .GlobalEnv))
# .ERPMenv <- new.env(parent = baseenv())

with(.GlobalEnv, {

  # ======================================== INIT ======================================== 
  # ==== Charge les settings et force les variables dans .GlobalEnv ====
  source("scripts/settings.R", local = FALSE)

  # ==== Vérifie que toutes les variables sont bien là ====
  required_globals <- c("VERBOSE", "DEV_MODE", "PROJECT_ROOT")
  missing_globals <- required_globals[!sapply(required_globals, exists, envir = .GlobalEnv)]
  if (length(missing_globals) > 0) {
      stop(red(paste( "Les variables globales suivantes ne sont pas définies dans scripts/settings.R:", paste(missing_globals, collapse = ", "))))
  }
  
  # ==== Récupère les packages depuis DESCRIPTION ====
  deps_info <- export_dependencies_from_description_file("DESCRIPTION")

  # ==== Version de R ====
  if (!is.null(deps_info$r_version)) {
    if (getRversion() < deps_info$r_version) {
      stop(sprintf("R >= %s requis par DESCRIPTION", deps_info$r_version))
    }
  }

  # ==== Modules à charger ( pour le moment : pas de gestion des versions des modules) ====
  base_pkgs <- deps_info$packages
  base_pkgs <- unique(c(base_pkgs, "parallel", "igraph", "RColorBrewer", "devtools")) # Ajouter des packages supplémentaires nécessaires
  installation_status <- install_and_import_if_missing(base_pkgs, TRUE)

  # ==== Vérifie si 'devtools' a été chargé avec succès ====
  devtools_index <- which(installation_status$packages == "devtools")

  if (length(devtools_index) == 0 || !installation_status$packages_loaded[devtools_index]) {
    msg <- "'devtools' n'est pas disponible — chargement manuel des fichiers R."
    if (exists("log_msg")) log_msg("WARN", msg) else warning(msg)

    # Fallback : charge les fichiers du dossier R manuellement
    r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
    invisible(sapply(r_files, source, local = .GlobalEnv))

  } else {
    if (exists("log_msg")) log_msg("INFO", "'devtools' détecté — utilisation de devtools::load_all().")
    else message("'devtools' détecté — utilisation de devtools::load_all().")
    devtools::load_all()
  }

  # ==== Ajout du patch ergm ====
  source("scripts/ergm_patch.R", local = TRUE)
  res_tests <- ergm_patch_selftest(run_diagnostics = FALSE, verbose = TRUE) # selftest pour tester le patch sur ergm

  ergm_patch_enable() # Appliquer le patch pour la suite

  log_msg("INFO", "Démarrage du script pour ERPM")

  # ======================================== RUN ======================================== 
  # Création du réseau et récupération de la partition
  res <- create_erpm_network()
  log_erpm_network(res)

  # # Extraction des données de réseau
  nw <- res$network
  partition <- res$partition
  na_nodes <- names(which(is.na(partition))) # Objets sans groupe
  
  # # Plots the bipartite network with the partition clusters.
  plot_partition_clusters(nw)

  # Estimation du modèle ergm
  options(error = recover) # Pour le debug
  ergm_model <- ergm(nw ~ b2degrange(from=2, to=3))

  # ======================================== CLEAN ======================================== 
  log_msg("INFO", "Fin du programme -- Nettoyage de l'environnement global")

  # Désactive le patch 
  ergm_patch_enable()

  # Liste blanche : variables à conserver
  if (exists("DEV_MODE")){
    keep_vars <- c("")
  } else {
    keep_vars <- c(
      "VERBOSE", "DEV_MODE", "PROJECT_ROOT", "CRAN_MIRROR", "INC_DIR",
      "log_msg", "red", "green", "yellow", "blue", "bold",
      ".__settings_loaded"
    )
  }


  # Supprime toutes les autres variables globales
  all_vars <- ls(envir = .GlobalEnv, all.names = TRUE)
  vars_to_remove <- setdiff(all_vars, keep_vars)
  if (length(vars_to_remove) > 0) {
    rm(list = vars_to_remove, envir = .GlobalEnv)
  }

  gc()  # Nettoyage mémoire
  cat("\n✔ Environnement nettoyé avec succès.\n")

})