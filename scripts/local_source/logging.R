# ==============================================================================
# Fichier : logging.R
# Fonction : log_msg
# Utilité : Outils de logging pour l'execution du code
# ==============================================================================
if (!exists(".__logging_loaded", envir = .GlobalEnv)) {

  # ==============================================================================
  #' Fonction de log général
  #'
  #' Cette fonction enregistre un message dans le fichier de log global et l'affiche
  #' dans la console si le mode VERBOSE est activé. Les messages peuvent être de type
  #' \code{INFO}, \code{WARN}, \code{ERROR}, ou \code{SUCCESS}.
  #'
  #' @param type Type de message (\code{character}) : "INFO", "WARN", "ERROR" ou "SUCCESS".
  #' @param msg Message (\code{character}) à enregistrer.
  #' @param file Chemin vers le fichier de log (\code{character}). Par défaut, utilise \code{LOG_FILE}.
  #'
  #' @return Renvoie le chemin complet du fichier de log utilisé (\code{character}).
  #' @examples
  #' log_msg("INFO", "Début du traitement")
  #' log_msg("WARN", "Attention, données manquantes")
  #' log_msg("ERROR", "Impossible de charger le package")
  #' log_msg("SUCCESS", "Traitement terminé avec succès")
  log_msg <- function(type, msg, file = LOG_FILE) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    line <- sprintf("[%s] [%s] %s\n", timestamp, toupper(type), msg)
    
    # Écrit dans le fichier de log
    cat(line, file = file, append = TRUE)
    
    # Affiche aussi dans la console si VERBOSE activé
    if (exists("VERBOSE") && VERBOSE) {
      if (type == "ERROR") cat(bold(red(line)))
      else if (type == "WARN") cat(bold(yellow(line)))
      else if (type == "INFO") cat(bold(cyan(line)))
      else if (type == "SUCCESS") cat(bold(green(line)))
    }
    
    LOG_FILE # Renvoie le nom du fichier log
  }

  # ==============================================================================
  #' Capture et log de toute la sortie d'une expression R
  #'
  #' Cette fonction exécute une expression R tout en capturant sa sortie standard
  #' et ses messages d'avertissement/erreur, puis les envoie au logger via
  #' \code{log_msg()}. Utile pour capturer la sortie de fonctions générant beaucoup
  #' de messages console comme \code{ergm()}.
  #'
  #' @param expr Expression R à exécuter (\code{expression} ou \code{{} block}).
  #' @param type Type de message (\code{character}) pour le log : "INFO" par défaut.
  #'
  #' @return Retourne la valeur renvoyée par l'expression (\code{any}).
  #' @examples
  #' log_capture({
  #'   print("Ceci sera redirigé vers le log")
  #'   message("Message interne également loggé")
  #' })
  log_capture <- function(expr, type = "INFO", verbose = TRUE) {
    # Conserve le call si verbose = FALSE
    call_text <- deparse(substitute(expr))

    # Redirige stdout et messages
    con <- textConnection("tmp_output", "w", local = TRUE)
    sink(con, type = "output")
    sink(con, type = "message")
    
    on.exit({
      sink(type = "message")
      sink(type = "output")
      close(con)
      
      # Si verbose TRUE, log complet, sinon juste le call
      output_to_log <- if(isTRUE(verbose)) paste(tmp_output, collapse = "\n") else paste("Call:", call_text)
      
      if (exists("log_msg")) log_msg(type, output_to_log)
      else cat(output_to_log, "\n")
    }, add = TRUE)
    
    force(expr)
  }

  # Dossier de logs global
  LOG_DIR <- "logs"
  if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
  
  # Nom du fichier de log avec horodatage
  timestamp_file <- format(Sys.time(), "%Y%m%d_%H%M%S")
  LOG_FILE <- file.path(LOG_DIR, paste0("cubi_", timestamp_file, ".log"))
  
  log_msg("INFO", sprintf("Fichier de log créé : %s", LOG_FILE))

  # --- Marqueur interne pour indiquer que le module est chargé ---
  assign("log_file",          LOG_FILE,     envir = .GlobalEnv)
  assign("log_msg",           log_msg,      envir = .GlobalEnv)
  assign("log_capture",       log_capture,  envir = .GlobalEnv)
  assign(".__logging_loaded", TRUE,      envir = .GlobalEnv)

}