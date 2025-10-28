# ==============================================================================
# Fichier : debug.R
# Fonction : debug_source_file()
# Utilité : Diagnostiquer une erreur au chargement d’un fichier R (ex. lors d’un source())
# ==============================================================================

if (!exists(".__debug_loaded", envir = .GlobalEnv)) {

    #' Diagnostiquer une erreur lors du sourcing d’un fichier R
    #'
    #' Lit un fichier R, le parse avec conservation des références source, puis
    #' évalue les expressions **une par une** dans un environnement isolé pour
    #' identifier la première qui échoue.
    #'
    #' Affiche avant chaque évaluation un **aperçu lisible du code** (premières lignes
    #' et/ou type d’expression). Si le fichier est un gros bloc `if (...) { ... }`,
    #' la fonction le déplie pour exécuter chaque sous-instruction.
    #'
    #' À la fin, un résumé global indique le statut, le nombre d’expressions exécutées
    #' et les objets créés.
    #'
    #' @param file Chemin du fichier R à diagnostiquer.
    #' @param encoding Encodage du fichier (défaut \code{"UTF-8"}).
    #' @param preview_lines Nombre maximal de lignes de code à afficher pour chaque expression.
    #' @return Invisiblement, une liste \code{list(success, evaluated, failed_index, env)}.
    #' @examples
    #' \dontrun{
    #'   res <- debug_source_file("R/erpm_wrapper.R")
    #'   ls(res$env)
    #' }
    #' @export
    debug_source_file <- function(file = "R/erpm_wrapper.R",
                                encoding = "UTF-8",
                                preview_lines = 4L) {
    # Helper d'encodage sûr pour la console
    .u <- function(x) enc2native(enc2utf8(x))

    exprs <- parse(file, keep.source = TRUE, encoding = encoding)
    diag_env <- new.env(parent = baseenv())

    writeLines(.u(sprintf("\n=== Diagnostic du fichier : %s ===", file)))

    # Déplie les blocs if / { ... }
    decompose_block <- function(expr, env) {
        if (is.null(expr) || !is.call(expr)) return(list(expr))
        head <- expr[[1L]]
        if (identical(head, as.name("{"))) {
        as.list(expr)[-1L]
        } else if (identical(head, as.name("if"))) {
        cond <- try(eval(expr[[2L]], envir = env), silent = TRUE)
        branch <- if (!inherits(cond, "try-error") && isTRUE(cond)) {
            expr[[3L]]
        } else if (length(expr) >= 4L) expr[[4L]] else NULL
        if (is.null(branch)) return(list())
        if (is.call(branch) && identical(branch[[1L]], as.name("{"))) {
            as.list(branch)[-1L]
        } else list(branch)
        } else {
        list(expr)
        }
    }

    top <- if (length(exprs) >= 1L) exprs[[1L]] else NULL
    to_eval <- if (length(exprs) == 1L && is.call(top) &&
                    (identical(top[[1L]], as.name("{")) || identical(top[[1L]], as.name("if")))) {
        decompose_block(top, diag_env)
    } else {
        as.list(exprs)
    }

    n <- length(to_eval)
    evaluated <- 0L
    failed_index <- NA_integer_

    for (i in seq_len(n)) {
        expr <- to_eval[[i]]

        # Crée un aperçu du code
        code_preview <- paste0(capture.output(print(expr)), collapse = "\n")
        code_lines <- strsplit(code_preview, "\n")[[1]]
        if (length(code_lines) > preview_lines)
        code_lines <- c(code_lines[1:preview_lines],
                        sprintf("... (+%d lignes)", length(code_lines) - preview_lines))

        writeLines(.u(sprintf(
        "\n-- Évaluation de l’expression %d/%d --\n%s",
        i, n, paste0(">> ", code_lines, collapse = "\n")
        )))

        ok <- try(eval(expr, envir = diag_env), silent = TRUE)

        if (inherits(ok, "try-error")) {
        writeLines(.u(sprintf("✗ ERREUR sur l’expression #%d", i)))
        sr <- attr(expr, "srcref")
        if (!is.null(sr)) {
            writeLines(.u(sprintf("Lignes %d:%d → %d:%d", sr[[1]], sr[[2]], sr[[3]], sr[[4]])))
        }
        failed_index <- i
        break
        } else {
        evaluated <- evaluated + 1L
        }
    }

    success <- is.na(failed_index)
    objs <- ls(envir = diag_env, all.names = TRUE)

    writeLines(.u("\n=== Résumé diagnostic ==="))
    writeLines(.u(sprintf("Statut        : %s", if (success) "SUCCÈS" else "ÉCHEC")))
    writeLines(.u(sprintf("Évaluées      : %d / %d", evaluated, n)))
    if (!success)
        writeLines(.u(sprintf("Première erreur à l’expression #%d", failed_index)))
    writeLines(.u(sprintf("Objets créés  : %d", length(objs))))
    if (length(objs))
        writeLines(.u(sprintf(" - %s", paste(objs, collapse = ", "))))

    invisible(list(success = success,
                    evaluated = evaluated,
                    failed_index = failed_index,
                    env = diag_env))
    }


    # --- Attributions à l’environnement global ---
    assign("debug_source_file",         debug_source_file,      envir = .GlobalEnv)
    assign(".__debug_loaded",           TRUE,                   envir = .GlobalEnv)

}