# ====================================================================================== 
# Fichier : cub_test.R
# Fonction :
# Utilité : test pour les nouvelles fonctions ERPM
# ====================================================================================== 

#with(.GlobalEnv, {


  # ====================================================================================== 
  # ======================================== INIT ======================================== 
  source("scripts/local_source/settings.R", local = FALSE)
  source("scripts/local_source/init.R", local = FALSE)
  

  init_erpm(selftest=FALSE, verbose=FALSE)

  if (exists("log_msg")) {
    log_msg("INFO", "Démarrage du script pour ERPM")
  } else {
    message("Démarrage du script pour ERPM")
  }


  if (exists("ergm_patch_enable")) ergm_patch_enable(verbose = VERBOSE)



  # ====================================================================================== 
  # ======================================== RUN ========================================= 
  # Création du réseau et récupération de la partition
  res <- create_erpm_network()
  log_erpm_network(res)

  # Extraction des données de réseau
  nw <- res$network
  partition <- res$partition
  na_nodes <- names(which(is.na(partition))) # Objets sans groupe

  # Plot du réseau bipartite avec les clusters
  plot_partition_clusters(nw)

  # Estimation du modèle ERGM
  options(error = recover) # Pour le debug

  # Exécute le modèle ERGM et capture la sortie
  ergm_model <- log_capture({
    ergm(nw ~ b2degrange(from = 2, to = 3), verbose = 0)
  }, type = "INFO")

  # Résumé automatique du modèle ERGM
  #summary_ergm_model(ergm_model)
  summary_ergm_model(ergm_model, nw = nw, partition = partition)

  # Désactive le patch ERGM si actif
  if (exists("ergm_patch_disable")) ergm_patch_disable(verbose = VERBOSE)

  



  # ======================================================================================= 
  # ======================================== CLEAN ======================================== 
  # Vide le buffer stdout
  flush.console()  

  # Log
  log_msg("INFO", "Fin du programme -- Nettoyage de l'environnement global")

  # Nettoie l'environnement
  clean_global_env(verbose = VERBOSE)

#})