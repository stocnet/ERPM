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
  erpm_load_with_recompile <- function(pkg_name = "ERPM",
                                      src_dir = "src",
                                      quiet   = TRUE,
                                      install_verbose = FALSE) {
    requireNamespace("ergm", quietly = TRUE)

    dynext  <- .Platform$dynlib.ext
    dll_path<- file.path(src_dir, paste0(pkg_name, dynext))
    cfiles  <- list.files(src_dir, pattern = "\\.(c|cc|cpp)$", full.names = TRUE)

    need_rebuild <- function() {
      if (!length(cfiles)) return(FALSE)
      if (!file.exists(dll_path)) return(TRUE)
      max(file.info(cfiles)$mtime, na.rm = TRUE) > file.info(dll_path)$mtime
    }

    if (requireNamespace("devtools", quietly = TRUE)) {

      if (isTRUE(install_verbose)) {
        # Build verbeux + logs détaillés
        # old <- Sys.getenv(c("MAKEFLAGS","R_KEEP_PKG_SOURCE","PKG_CFLAGS","PKG_CPPFLAGS"), unset = NA)
        # on.exit({
        #   # restore
        #   mapply(Sys.setenv,
        #         names(old),
        #         ifelse(is.na(old), "", old))
        # }, add = TRUE)

        old <- Sys.getenv(c("MAKEFLAGS","R_KEEP_PKG_SOURCE","PKG_CFLAGS","PKG_CPPFLAGS"),
                  unset = NA_character_)

        on.exit({
          # rétablir celles qui existaient
          to_set <- as.list(old[!is.na(old)])
          if (length(to_set)) do.call(Sys.setenv, to_set)
          # retirer celles qui n'existaient pas
          to_unset <- names(old)[is.na(old)]
          if (length(to_unset)) Sys.unsetenv(to_unset)
        }, add = TRUE)

        Sys.setenv(
          MAKEFLAGS            = "-j1",
          R_KEEP_PKG_SOURCE    = "yes",
          PKG_CFLAGS           = "-std=c11 -O0 -g3 -Wall -Wextra -pedantic -ferror-limit=0 -fmacro-backtrace-limit=0",
          PKG_CPPFLAGS         = "-DDEBUG_COV_FULLMATCH=1 -fdiagnostics-color=never"
        )

        devtools::install(
          ".", upgrade = "never", quiet = FALSE,
          args = c("--no-clean-on-error","--preclean","-v")
        )
        # recharge le code R après install pour garder l’environnement dev
        devtools::load_all(recompile = FALSE, quiet = quiet)
        return(TRUE)
      } else {
        # Flux rapide: recompile si nécessaire, sinon charge
        Sys.setenv(
          MAKEFLAGS    = "-j1",
          PKG_CFLAGS   = "-std=c11 -O2 -Wall",
          PKG_CPPFLAGS = ""
        )
        devtools::load_all(recompile = need_rebuild(), quiet = quiet)
        return(TRUE)
      }
    }

    # Fallbacks identiques à ta version…
    r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
    invisible(lapply(r_files, function(f) try(source(f, local = .GlobalEnv), silent = TRUE)))

    if (requireNamespace("pkgbuild", quietly = TRUE)) {
      if (need_rebuild()) pkgbuild::compile_dll(path = ".", quiet = quiet)
      if (file.exists(dll_path) && !isTRUE(dll_path %in% names(getLoadedDLLs())))
        try(dyn.load(dll_path), silent = TRUE)
      return(TRUE)
    }

    # SHLIB direct (inchangé)
    if (need_rebuild()) {
      rincl   <- sprintf('-I"%s"', R.home("include"))
      ergminc <- system.file("include", package = "ergm")
      eincl   <- if (nzchar(ergminc)) sprintf('-I"%s"', ergminc) else ""
      cmd <- sprintf('R CMD SHLIB %s -o %s %s %s',
                    paste(shQuote(cfiles), collapse = " "),
                    shQuote(dll_path), rincl, eincl)
      status <- system(cmd)
      if (status != 0) stop("Échec compilation native: ", cmd)
    }
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
    selftest = FALSE, 
    install_verbose = FALSE
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
    erpm_load_with_recompile( pkg_name        = "ERPM",
                              src_dir         = "src",
                              quiet           = !isTRUE(verbose),
                              install_verbose = install_verbose)

    # Patch ergm
    source(ergm_patch_file, local = TRUE)
    if (isTRUE(selftest)) ergm_patch_selftest(run_diagnostics = FALSE, verbose = verbose)

    invisible(TRUE)
  }

  assign("erpm_load_with_recompile", erpm_load_with_recompile, envir = .GlobalEnv)
  assign("init_erpm", init_erpm, envir = .GlobalEnv)
  assign(".__init_loaded", TRUE, envir = .GlobalEnv)
}