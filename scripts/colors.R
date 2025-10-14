# ==============================================================================
# Fichier : colors.R
# Fonction : 
# Utilité : Définit des couleurs pour le projet
# ==============================================================================

if (!exists(".__colors_loaded", envir = .GlobalEnv)) {

  assign("green",  function(x) paste0("\033[32m", x, "\033[0m"), envir = .GlobalEnv)
  assign("red",    function(x) paste0("\033[31m", x, "\033[0m"), envir = .GlobalEnv)
  assign("yellow", function(x) paste0("\033[33m", x, "\033[0m"), envir = .GlobalEnv)
  assign("blue",   function(x) paste0("\033[34m", x, "\033[0m"), envir = .GlobalEnv)
  assign("bold",   function(x) paste0("\033[1m",  x, "\033[0m"), envir = .GlobalEnv)

  assign(".__colors_loaded", TRUE, envir = .GlobalEnv)
}