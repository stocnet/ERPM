# ==============================================================================
# Fichier : logging.R
# Fonction : log_msg
# Utilité : Outils de logging pour l'execution du code
# ==============================================================================
if (!exists(".__logging_loaded", envir = .GlobalEnv)) {

  # Dossier de logs global
  LOG_DIR <- "logs"
  if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
  
  # Nom du fichier de log avec horodatage
  timestamp_file <- format(Sys.time(), "%Y%m%d_%H%M%S")
  LOG_FILE <- file.path(LOG_DIR, paste0("cubi_", timestamp_file, ".log"))
  
  # Fonction globale de log
  log_msg <- function(type, msg, file = LOG_FILE) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    line <- sprintf("[%s] [%s] %s\n", timestamp, toupper(type), msg)
    
    # Écrit dans le fichier de log
    cat(line, file = file, append = TRUE)
    
    # Affiche aussi dans la console si VERBOSE activé
    if (exists("VERBOSE") && VERBOSE) {
      if (type == "ERROR") cat(red(line))
      else if (type == "WARN") cat(yellow(line))
      else cat(green(line))
    }
  }

  log_msg("INFO", sprintf("Fichier de log créé : %s", LOG_FILE))

  source("scripts/colors.R")
  cat("\n", yellow(sprintf("\n→ Log complet enregistré dans : %s\n", LOG_FILE)))

  # --- Marqueur interne pour indiquer que le module est chargé ---
  assign("log_file",          LOG_FILE,  envir = .GlobalEnv)
  assign("log_msg",           log_msg,   envir = .GlobalEnv)
  assign(".__logging_loaded", TRUE,      envir = .GlobalEnv)

}