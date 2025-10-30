# ==============================================================================
# Fichier : init.R
# Fonction : init_erpm
# Utilité : Initialiser l'environnement ERPM, charger les settings, packages,
#           recompiler src/ au besoin, et activer le patch ERGM
# ==============================================================================

if(!exists(".__init_loaded", envir = .GlobalEnv)){

  #' Charger le code R du package et recompiler src/ si nécessaire
  #'
  #' Cette fonction centralise le chargement dynamique du package ERPM pendant le
  #' développement. Elle :
  #' - détecte tous les fichiers `.c`, `.cc`, `.cpp` dans `src/` ;
  #' - compare leurs dates de modification avec la DLL compilée (`ERPM.so`) ;
  #' - recompile automatiquement si un fichier natif est plus récent ou absent ;
  #' - recharge le code R via `devtools::load_all()` si disponible, sinon `pkgbuild` ou `R CMD SHLIB`.
  #'
  #' @param pkg_name Nom du package (utilisé pour la DLL, défaut `"ERPM"`).
  #' @param src_dir Répertoire contenant les fichiers C/C++ (défaut `"src"`).
  #' @param quiet Si TRUE, réduit le bruit d’affichage pendant la compilation.
  #' @return TRUE si le code a été chargé avec succès ; stoppe en cas d’échec de compilation.
  erpm_load_with_recompile <- function(pkg_name = "ERPM", src_dir = "src", quiet = TRUE) {
    requireNamespace("ergm", quietly = TRUE)

    # Extension dynamique selon le système (.so ou .dll)
    dynext <- .Platform$dynlib.ext

    # Chemin attendu de la librairie compilée
    dll_path <- file.path(src_dir, paste0(pkg_name, dynext))

    # Liste des fichiers source C/C++
    cfiles <- list.files(src_dir, pattern = "\\.(c|cc|cpp)$", full.names = TRUE)

    # --- Fonction interne : déterminer si une recompilation est nécessaire ---
    need_rebuild <- function() {
      if (!length(cfiles)) return(FALSE)                      # Aucun fichier source
      if (!file.exists(dll_path)) return(TRUE)                # Pas encore de DLL compilée
      max(file.info(cfiles)$mtime, na.rm = TRUE) > file.info(dll_path)$mtime  # Au moins un .c plus récent
    }

    # --- 1) Si devtools est disponible ---
    if (requireNamespace("devtools", quietly = TRUE)) {

      # Recompile seulement si besoin, sinon simple rechargement
      devtools::load_all(recompile = need_rebuild(), quiet = quiet)
      return(TRUE)
    }

    # --- 2) Si devtools absent : recharger manuellement tous les R/*.R ---
    r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
    invisible(lapply(r_files, function(f) try(source(f, local = .GlobalEnv), silent = TRUE)))

    # --- 3) Si pkgbuild est disponible : l’utiliser pour compiler les DLL ---
    if (requireNamespace("pkgbuild", quietly = TRUE)) {
      requireNamespace("ergm", quietly = TRUE)
      if (need_rebuild()) pkgbuild::compile_dll(path = ".", quiet = quiet)

      # Charger la DLL si non encore chargée
      if (file.exists(dll_path) && !isTRUE(dll_path %in% names(getLoadedDLLs())))
        try(dyn.load(dll_path), silent = TRUE)

      return(TRUE)
    }

    # --- 4) Fallback : appel direct à R CMD SHLIB ---
    if (need_rebuild()) {

      # Récupère l’inclusion des headers ergm (utile pour les changestats)
      inc_ergm <- try(system.file("include", package = "ergm"), silent = TRUE)
      inc_flag <- if (!inherits(inc_ergm, "try-error") && nzchar(inc_ergm))
        sprintf('-I"%s"', inc_ergm) else ""

      # Construction de la commande de compilation
      cmd <- sprintf(
        'R CMD SHLIB %s -o %s %s',
        paste(shQuote(cfiles), collapse = " "),
        shQuote(dll_path),
        inc_flag
      )
      
      # Exécution et contrôle d’erreur
      status <- system(cmd)
      if (status != 0) stop("Échec compilation native: ", cmd)
    }

    # --- 5) Chargement dynamique final si la DLL existe ---
    if (file.exists(dll_path) && !isTRUE(dll_path %in% names(getLoadedDLLs())))
      try(dyn.load(dll_path), silent = TRUE)

    TRUE
  }


  #' Initialisation de l'environnement ERPM
  #'
  #' - Charge les settings globaux.
  #' - Vérifie/installe les packages (DESCRIPTION).
  #' - Charge le code R et recompile src/ si nécessaire.
  #' - Ajoute le patch ergm et exécute un self-test optionnel.
  #'
  #' @param ergm_patch_file Chemin du patch ERGM.
  #' @param description_file Fichier DESCRIPTION.
  #' @param verbose Affichages.
  #' @param selftest Exécuter le self-test du patch.
  #' @return TRUE (invisible).
  #' @export
  init_erpm <- function(
    ergm_patch_file = "scripts/ergm_patch.R",
    description_file = "DESCRIPTION",
    verbose = TRUE,
    selftest = FALSE
  ){
    # Settings
    if (!exists(".__settings_loaded", envir = .GlobalEnv)) {
      source("scripts/local_source/settings.R", local = FALSE)
    }

    # Garde-fous
    required_globals <- c("VERBOSE", "DEV_MODE", "PROJECT_ROOT", "INC_DIR", "CRAN_MIRROR")
    missing_globals <- required_globals[!sapply(required_globals, exists, envir = .GlobalEnv)]
    if (length(missing_globals) > 0) {
      stop(
        "Variables globales manquantes dans scripts/local_source/settings.R : ",
        paste(missing_globals, collapse = ", ")
      )
    }

    # Lanceur
    if (!exists(".__launcher_loaded", envir = .GlobalEnv)) {
      source("scripts/local_source/launcher.R", local = FALSE)
    }

    # Dépendances
    deps_info <- export_dependencies_from_description_file(description_file)
    if (!is.null(deps_info$r_version) && getRversion() < deps_info$r_version) {
      stop(sprintf("R >= %s requis par %s", deps_info$r_version, description_file))
    }
    base_pkgs <- unique(c(deps_info$packages, "parallel", "igraph", "RColorBrewer", "devtools", "R.utils", "pkgbuild"))
    install_and_import_if_missing(base_pkgs, TRUE)

    # Chargement + recompilation conditionnelle
    erpm_load_with_recompile(pkg_name = "ERPM", src_dir = "src", quiet = !isTRUE(verbose))

    # Patch ergm
    source(ergm_patch_file, local = TRUE)
    if (isTRUE(selftest)) ergm_patch_selftest(run_diagnostics = FALSE, verbose = verbose)

    invisible(TRUE)
  }

  assign("erpm_load_with_recompile", erpm_load_with_recompile, envir = .GlobalEnv)
  assign("init_erpm", init_erpm, envir = .GlobalEnv)
  assign(".__init_loaded", TRUE, envir = .GlobalEnv)
}