# ==============================================================================
# Fichier : logging.R
# Fonction : log_msg, log_msgf, safe_sprintf, log_capture
# Utilité  : Outils de logging pour l'exécution du code
# ==============================================================================

if (!exists(".__logging_loaded", envir = .GlobalEnv)) {

  # ----------------------------------------------------------------------------
  #' sprintf sécurisé
  #'
  #' Évalue `sprintf(fmt, ...)` avec repli automatique si une conversion flottante
  #' (`%f`, `%.xf`) échoue. En cas d'erreur, remplace tous les `%f` par `%s` et
  #' convertit les arguments non numériques en chaînes.
  #'
  #' @param fmt Chaîne de format (`character(1)`).
  #' @param ... Arguments à formater.
  #'
  #' @return Une `character(1)` formattée. En ultime recours, un message d'erreur
  #'         explicite est renvoyé.
  #' @examples
  #' safe_sprintf("x=%.2f", 3.14159)
  #' safe_sprintf("x=%.2f y=%.2f", "NA", 2.3)  # -> remplace par %s pour le premier
  safe_sprintf <- function(fmt, ...) {
    # Capture des args pour double tentative
    args <- list(...)

    # Tentative 1 : sprintf direct
    tryCatch(
      do.call(sprintf, c(list(fmt), args)),
      error = function(e) {
        # Tentative 2 : repli — remplacer tous les %f (et %.xf) par %s
        fmt2  <- gsub("%(\\.[0-9]+)?f", "%s", fmt, perl = TRUE)

        # Convertir en texte tout ce qui n'est pas un numérique scalaire fini
        args2 <- lapply(
          args,
          function(x) {
            if (is.numeric(x) && length(x) == 1L && is.finite(x)) x else as.character(x)
          }
        )

        # Dernière tentative, sinon message d'erreur textuel
        tryCatch(
          do.call(sprintf, c(list(fmt2), args2)),
          error = function(e2) paste("<<format error:", conditionMessage(e), "fmt=", fmt, ">>")
        )
      }
    )
  }

  # ----------------------------------------------------------------------------
  #' Logger principal
  #'
  #' Enregistre un message dans un fichier de log et, si `VERBOSE==TRUE`, l'affiche
  #' en console avec couleurs (si disponibles). Accepte tout objet en `msg`.
  #'
  #' @param type Type de message (`character(1)`): "INFO", "WARN", "ERROR", "SUCCESS".
  #' @param msg  Message à enregistrer. Tout type d'objet accepté; sera converti
  #'             en texte de manière robuste.
  #' @param file Chemin du fichier de log (`character(1)`). Par défaut `LOG_FILE`.
  #'
  #' @return Le chemin du fichier de log utilisé.
  #' @examples
  #' log_msg("INFO", "Début du traitement")
  #' log_msg("ERROR", try(stop("Oups"), silent=TRUE))
  log_msg <- function(type, msg, file = LOG_FILE) {
    # -- Normalisation du type de message
    type2 <- tryCatch(toupper(as.character(type)[1L]), error = function(e) "INFO")

    # -- Conversion robuste en texte ------------------------------------------
    .to_string <- function(x) {
      # NULL explicite
      if (is.null(x)) return("<NULL>")

      # try-error : préfixer pour le distinguer visuellement
      if (inherits(x, "try-error")) return(paste0("ERROR: ", as.character(x)))

      # Atomique scalaire
      if (is.atomic(x) && length(x) == 1L) return(as.character(x))

      # Atomique vectoriel court (couper à 50 items)
      if (is.atomic(x) && length(x) > 1L)
        return(paste(utils::head(as.character(x), 50L), collapse = ", "))

      # Fallback structuré : str() puis print() au besoin
      s <- paste(capture.output(utils::str(x, give.head = FALSE, vec.len = 20L)), collapse = " ")
      if (!nzchar(s)) s <- paste(capture.output(print(x)), collapse = " ")
      s
    }

    # Protéger la stringification à tout prix
    msg2 <- tryCatch(
      if (is.character(msg) && length(msg) == 1L) msg else .to_string(msg),
      error = function(e) paste0("<<log_msg stringification error: ", conditionMessage(e), ">>")
    )

    # -- Construction de la ligne de log --------------------------------------
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    line <- sprintf("[%s] [%s] %s\n", timestamp, type2, msg2)

    # -- Écriture fichier (fallback console si 'file' invalide) ----------------
    file_ok <- is.character(file) && length(file) == 1L && nzchar(file)
    if (file_ok) {
      # Créer le répertoire si besoin
      dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
      # Écriture silencieuse (protéger contre E/S)
      try(cat(line, file = file, append = TRUE), silent = TRUE)
    } else {
      cat(line)
    }

    # -- Echo console optionnel avec couleurs ---------------------------------
    if (isTRUE(get0("VERBOSE", ifnotfound = FALSE))) {
      # Sans dépendance dure à {crayon} :
      safe_color <- function(f, x) { ff <- get0(f); if (is.function(ff)) ff(x) else x }
      colored <- switch(
        type2,
        "ERROR"   = safe_color("bold", safe_color("red",    line)),
        "WARN"    = safe_color("bold", safe_color("yellow", line)),
        "INFO"    = safe_color("bold", safe_color("cyan",   line)),
        "SUCCESS" = safe_color("bold", safe_color("green",  line)),
        line
      )
      cat(colored)
    }

    # Retourner le chemin du fichier (ou chaîne vide si absent)
    file %||% ""
  }

  # ----------------------------------------------------------------------------
  #' Logger formaté (style printf)
  #'
  #' Comme \code{log_msg()}, mais accepte un format et des arguments à la façon
  #' de \code{sprintf}. Sécurisé via \code{safe_sprintf()} pour éviter les erreurs
  #' de formatage lorsqu'un argument n'est pas numérique.
  #'
  #' @param type Type de message (`character(1)`).
  #' @param fmt  Format \code{sprintf} (`character(1)`).
  #' @param ...  Arguments à interpoler dans \code{fmt}.
  #' @param file Chemin du fichier de log (`character(1)`).
  #'
  #' @return Le chemin du fichier de log utilisé.
  #' @examples
  #' log_msgf("INFO", "x=%.3f, y=%s", 3.14159, "ok")
  #' log_msgf("INFO", "obs=%.6f exp=%.6f", "try-error: ..", 0.123456) # -> %f -> %s pour obs
  log_msgf <- function(type, fmt, ..., file = LOG_FILE) {
    # Construire la ligne via le printf sécurisé
    formatted <- safe_sprintf(fmt, ...)
    # Déléguer l'écriture au logger principal
    log_msg(type, formatted, file = file)
  }

  # ----------------------------------------------------------------------------
  #' Capturer et logger la sortie d'une expression
  #'
  #' Exécute une expression en capturant sa sortie standard et ses messages,
  #' puis injecte le résultat dans le logger. Utile pour encapsuler des appels
  #' bruyants comme \code{ergm()}.
  #'
  #' @param expr    Expression R à évaluer (passée non évaluée).
  #' @param type    Type de message pour le log (`character(1)`) ; "INFO" par défaut.
  #' @param verbose Si TRUE, log de la sortie complète ; sinon log du seul appel.
  #'
  #' @return La valeur renvoyée par \code{expr}.
  #' @examples
  #' log_capture({
  #'   print("Hello")
  #'   message("World")
  #' })
  log_capture <- function(expr, type = "INFO", verbose = TRUE) {
    # Conserver une trace compacte de l'appel si verbose=FALSE
    call_text <- deparse(substitute(expr))

    # Ouvrir une connexion de texte pour capturer stdout+messages
    con <- textConnection("tmp_output", "w", local = TRUE)
    sink(con, type = "output")
    sink(con, type = "message")

    # Toujours restaurer l'état des sinks, même si expr échoue
    on.exit({
      sink(type = "message")
      sink(type = "output")
      close(con)

      # Choix du contenu loggé
      output_to_log <- if (isTRUE(verbose)) {
        paste(tmp_output, collapse = "\n")
      } else {
        paste("Call:", call_text)
      }

      # Écrire via log_msg (ou cat si non disponible)
      if (exists("log_msg")) log_msg(type, output_to_log) else cat(output_to_log, "\n")
    }, add = TRUE)

    # Évaluer l'expression
    force(expr)
  }

  # ----------------------------------------------------------------------------
  # Initialisation du répertoire et du fichier de log
  # ----------------------------------------------------------------------------
  LOG_DIR <- "logs"                                      # Dossier par défaut
  if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)

  timestamp_file <- format(Sys.time(), "%Y%m%d_%H%M%S")  # Horodatage fichier
  LOG_FILE <- file.path(LOG_DIR, paste0("cubi_", timestamp_file, ".log"))

  # Premier message pour tracer l'initialisation
  log_msg("INFO", sprintf("Fichier de log créé : %s", LOG_FILE))

  # Exposer les fonctions/utilitaires dans l'environnement global
  assign("log_file",          LOG_FILE,    envir = .GlobalEnv)
  assign("log_msg",           log_msg,     envir = .GlobalEnv)
  assign("log_msgf",          log_msgf,    envir = .GlobalEnv)
  assign("safe_sprintf",      safe_sprintf,envir = .GlobalEnv)
  assign("log_capture",       log_capture, envir = .GlobalEnv)
  assign(".__logging_loaded", TRUE,        envir = .GlobalEnv)
}