# ==============================================================================
# Fichier : launcher.R
# Fonction : launch_model()
# Utilité : Lancer rapidement un test de stat/fit (summary | ergm | erpm) avec 1+ effets
# ==============================================================================

if (!exists(".__launcher_loaded", envir = .GlobalEnv)) {

  # ============================================================
  # Dépendances (packages + sources locales) — helpers
  # ============================================================

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  .cands <- function(rel) unique(c(
    rel,
    file.path("scripts", "local_source", rel),
    file.path("R", rel),
    file.path(.root, rel),
    file.path(.root, "scripts", "local_source", rel),
    file.path(.root, "R", rel)
  ))

  .ensure_pkgs <- function(pkgs) {
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss)) {
      message(sprintf("Installation des packages manquants : %s", paste(miss, collapse = ", ")))
      install.packages(miss, repos = getOption("repos", "https://cloud.r-project.org"))
    }
    invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(require(p, character.only = TRUE))))
  }

  .modify_list <- function(a, b) { a <- as.list(a); b <- as.list(b); for (nm in names(b)) a[[nm]] <- b[[nm]]; a }

  .source_if_found <- function(paths, local = FALSE, required = FALSE) {
    for (p in paths) if (file.exists(p)) {
      tryCatch({ source(p, local = local); return(invisible(TRUE)) },
               error = function(e) message(sprintf("Erreur source(\"%s\"): %s", p, conditionMessage(e))))
    }
    if (required) stop(sprintf("Introuvable : %s", paste(paths, collapse = " | ")))
    invisible(FALSE)
  }

  .head_lines <- function(x, n = 8L) if (length(x) <= n) x else c(x[seq_len(n)], sprintf("... (+%d lignes)", length(x) - n))
  .fmt_call <- function(x) {
    if (is.call(x) && identical(x[[1]], as.name("~")) && length(x) >= 3L) {
      lhs <- x[[2]]; if (inherits(lhs, "network")) x[[2]] <- as.name("nw")
    }
    out <- tryCatch(deparse(x, width.cutoff = 500L), error = function(e) sprintf("<deparse: %s>", conditionMessage(e)))
    paste(out, collapse = "  ")
  }
  .fmt_rhs <- function(rhs) paste(deparse(rhs, width.cutoff = 1000L), collapse = " ")
  .log_info <- function(txt) { if (exists("log_msg")) log_msg("INFO", txt) else message(txt) }

  .brief_result <- function(obj, engine = NULL) {
    if (identical(engine, "summary")) {
      co <- utils::capture.output(print(obj))
      if (!length(co)) co <- utils::capture.output(print.default(obj))
      return(paste(.head_lines(co), collapse = "\n"))
    }
    if (inherits(obj, "ergm")) {
      ll <- tryCatch(as.numeric(stats::logLik(obj)), error = function(e) NA_real_)
      k  <- tryCatch(length(stats::coef(obj)),       error = function(e) NA_integer_)
      return(sprintf("ergm fit: logLik=%.4f ; #coef=%s", ll, ifelse(is.na(k), "NA", k)))
    }
    if (is.call(obj)) return(paste("call:", .fmt_call(obj)))
    sprintf("Objet de classe: %s", paste(class(obj), collapse = ","))
  }

  .summarize_fit <- function(fit, nw, partition, verbose = TRUE) {
    if (exists("summary_ergm_model")) {
      send_to_log <- exists("log_msg")
      msg <- summary_ergm_model(fit, log = send_to_log, nw = nw, partition = partition)
      if (isTRUE(verbose) && !send_to_log) cat(msg, "\n")
      invisible(msg)
    } else {
      fit_brief <- .brief_result(fit, engine = "ergm")
      .log_info(paste("Résultat (ergm):", fit_brief))
      if (isTRUE(verbose)) cat(fit_brief, "\n")
      invisible(fit_brief)
    }
  }

  .normalize_wrappers <- function(expr) {
    if (!is.language(expr)) return(expr)
    fix_one <- function(sym, node) {
      if (is.call(node) && identical(node[[1]], as.name(sym))) {
        if (length(node) >= 2L) {
          arg1 <- node[[2]]
          if (!(is.call(arg1) && identical(arg1[[1]], as.name("~")))) node[[2]] <- as.call(list(as.name("~"), arg1))
        }
      }
      node
    }
    if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
      expr[[2]] <- .normalize_wrappers(expr[[2]])
      expr[[3]] <- .normalize_wrappers(expr[[3]])
      return(expr)
    }
    expr <- fix_one("Proj1", expr); expr <- fix_one("B", expr); expr
  }

  # ============================================================
  # Construction réseau / formule / effets
  # ============================================================

  # IMPORTANT: utilise systématiquement build_bipartite_from_inputs() du wrapper
  coerce_to_bipartite <- function(nw = NULL, partition = NULL, nodes = NULL, dyads = list(),
                                  verbose = TRUE, plot = FALSE) {
    if (!is.null(nw)) {
      bip <- tryCatch(nw %n% "bipartite", error = function(e) NULL)
      if (is.null(bip) || is.na(bip)) stop("nw fourni n'est pas biparti ou manque l'attribut %n% 'bipartite'.")
      # Partition best-effort depuis l’adjacence acteur→groupe
      n <- as.integer(bip)
      vn <- network::network.vertex.names(nw); if (is.null(vn)) vn <- seq_len(network::network.size(nw))
      actors <- vn[seq_len(n)]
      groups <- vn[seq.int(n + 1L, n + max(1L, length(vn) - n))]
      A <- as.matrix(nw, matrix.type = "adjacency")                 # recommandé

    

      part <- vapply(seq_len(n), function(i) {
        gi <- which(A[actors[i], groups] > 0L)
        if (length(gi)) as.integer(sub("^G", "", groups[gi[1L]])) else NA_integer_
      }, integer(1))
      return(list(network = nw, partition = part))
    }

    # Cas principal: construire via le wrapper
    if (is.null(partition)) stop("Fournir `partition` ou `nw`.")
    if (!exists("build_bipartite_from_inputs", mode = "function"))
      stop("build_bipartite_from_inputs() introuvable. Sourcez R/erpm_wrapper.R avant.")
    built <- build_bipartite_from_inputs(partition = partition, nodes = nodes, dyads = dyads)
    if (exists("log_erpm_network")) log_erpm_network(built, verbose = verbose)
    if (isTRUE(plot)) try(plot_partition_clusters(built$network), silent = TRUE)
    list(network = built$network, partition = built$partition)
  }

  build_formula_from_rhs <- function(rhs, nw = NULL) {
    rhs_expr <- if (is.character(rhs)) parse(text = rhs)[[1]] else rhs
    f <- as.formula(bquote(nw ~ .(rhs_expr)))
    if (!is.null(nw)) environment(f) <- list2env(list(nw = nw), parent = parent.frame())
    f
  }

  build_effect_call <- function(effect, effect_args = list()) {
    if (is.call(effect)) return(effect)
    if (is.character(effect) && length(effect) == 1L) return(as.call(c(as.name(effect), effect_args)))
    stop("`effect` doit être un nom d'effet (caractère) ou un appel R (call).")
  }

  build_rhs_from_effects <- function(effects, effects_args = list()) {
    if (is.character(effects) && length(effects) == 1L) return(parse(text = effects)[[1]])
    if (is.language(effects)) return(effects)
    if (is.list(effects) && length(effects) >= 1L && all(vapply(effects, is.call, logical(1)))) {
      if (length(effects) == 1L) return(effects[[1L]])
      return(Reduce(function(x, y) call("+", x, y), effects))
    }
    if (is.character(effects) && length(effects) >= 1L) {
      calls <- vector("list", length(effects))
      named <- length(names(effects_args)) > 0
      for (i in seq_along(effects)) {
        eff <- effects[i]
        args_i <- if (named) effects_args[[eff]] %||% list() else effects_args[[i]] %||% list()
        calls[[i]] <- build_effect_call(eff, args_i)
      }
      if (length(calls) == 1L) return(calls[[1L]])
      return(Reduce(function(x, y) call("+", x, y), calls))
    }
    stop("Format `effects` non reconnu.")
  }

  # ============================================================
  # Lancer un test/stat/fit (summary | ergm | erpm) avec 1+ effets
  # ============================================================

  launch_model <- function(engine = c("summary", "ergm", "erpm"),
                           effects,
                           effects_args = list(),
                           partition = NULL,
                           nodes = NULL,
                           dyads = list(),
                           nw = NULL,
                           dry_run = FALSE,
                           verbose = TRUE,
                           plot = FALSE,
                           constraints = NULL,
                           estimate = c("MLE", "MPLE", "CD"),
                           eval_loglik = FALSE,
                           control = list(),
                           timeout = NULL) {

    estimate <- match.arg(estimate)
    engine   <- match.arg(engine)

    # Réseau: si absent, construire via le wrapper
    if (is.null(nw)) {
      built <- coerce_to_bipartite(nw = NULL, partition = partition, nodes = nodes, dyads = dyads,
                                   verbose = verbose, plot = plot)
      nw <- built$network; partition <- built$partition
    } else {
      built <- coerce_to_bipartite(nw = nw, partition = NULL, verbose = verbose, plot = FALSE)
      partition <- built$partition
    }

    # RHS et wrappers
    rhs <- build_rhs_from_effects(effects, effects_args)
    rhs <- .normalize_wrappers(rhs)
    rhs_text <- .fmt_call(rhs)

    # Contrainte par défaut bipartite
    if (is.null(constraints)) {
      bip <- tryCatch(nw %n% "bipartite", error = function(e) NULL)
      if (!is.null(bip) && !is.na(bip)) constraints <- as.formula(~ b1part)
    }
    constraints_text <- if (is.null(constraints)) "NULL" else .fmt_call(constraints)

    .log_info(sprintf("launch_model(engine=%s, dry_run=%s) | RHS=%s | constraints=%s",
                      engine, dry_run, rhs_text, constraints_text))

    # ---------- summary ----------
    if (engine == "summary") {
      f <- build_formula_from_rhs(rhs, nw = nw)
      call_text <- .fmt_call(f)
      if (dry_run) {
        .log_info(paste("Résultat (dry-run/summary): summary(", call_text, ")"))
        return(list(engine=engine, mode="dryrun", rhs=rhs, rhs_text=rhs_text,
                    constraints=constraints, constraints_text=constraints_text,
                    call=f, call_text=call_text, fit=NULL, result=f,
                    network=nw, partition=partition))
      }
      out <- if (is.null(constraints)) summary(f) else summary(f, constraints = constraints)
      out_brief <- .brief_result(out, engine = "summary")
      .log_info(paste("Résultat (summary):\n", out_brief))
      return(list(engine=engine, mode="summary", rhs=rhs, rhs_text=rhs_text,
                  constraints=constraints, constraints_text=constraints_text,
                  call=f, call_text=call_text, fit=NULL, result=out,
                  network=nw, partition=partition))
    }

    # Prépare contrôles par défaut
    ctrl_defaults <- list(MCMLE.maxit = 5L, MCMC.samplesize = 2000L, MCMC.burnin = 10000L, MCMC.interval = 1000L)
    ctrl_list <- .modify_list(ctrl_defaults, control)
    ctrl <- do.call(ergm::control.ergm, ctrl_list)

    # ---------- ergm ----------
    if (engine == "ergm") {
      f <- build_formula_from_rhs(rhs, nw = nw)
      call_eval <- as.call(list(as.name("ergm"), f, constraints = constraints))
      call_text <- .fmt_call(call_eval)
      if (dry_run) {
        .log_info(paste("Résultat (dry-run/ergm):", call_text))
        return(list(engine=engine, mode="dryrun", rhs=rhs, rhs_text=rhs_text,
                    constraints=constraints, constraints_text=constraints_text,
                    call=call_eval, call_text=call_text, fit=NULL, result=NULL,
                    network=nw, partition=partition))
      }
      .run <- function() {
        if (is.null(constraints)) ergm::ergm(f, estimate = estimate, eval.loglik = eval_loglik, control = ctrl)
        else                      ergm::ergm(f, constraints = constraints, estimate = estimate, eval.loglik = eval_loglik, control = ctrl)
      }
      fit <- if (is.null(timeout)) .run() else R.utils::withTimeout(.run(), timeout = as.numeric(timeout), onTimeout = "silent")
      if (!inherits(fit, "ergm")) {
        .log_info("Fit non disponible ou interrompu."); 
        return(list(engine=engine, mode="fit_failed", rhs=rhs, rhs_text=rhs_text,
                    constraints=constraints, constraints_text=constraints_text,
                    call=NULL, call_text=call_text, fit=fit, result=fit,
                    network=nw, partition=partition))
      }
      .summarize_fit(fit, nw, partition, verbose)
      return(list(engine=engine, mode="fit", rhs=rhs, rhs_text=rhs_text,
                  constraints=constraints, constraints_text=constraints_text,
                  call=NULL, call_text=call_text, fit=fit, result=fit,
                  network=nw, partition=partition))
    }

    # ---------- erpm ----------
    if (!exists("erpm")) stop("La fonction `erpm()` est introuvable. Sourcez R/erpm_wrapper.R.")
    f <- build_formula_from_rhs(rhs, nw = nw)

    if (dry_run) {
      call_erpm <- erpm(f, eval_call = FALSE, verbose = verbose, estimate = estimate)
      call_text <- .fmt_call(call_erpm)
      .log_info(paste("Résultat (dry-run/erpm):", call_text))
      return(list(engine=engine, mode="dryrun", rhs=rhs, rhs_text=rhs_text,
                  constraints=constraints, constraints_text=constraints_text,
                  call=call_erpm, call_text=call_text, fit=NULL, result=NULL,
                  network=nw, partition=partition))
    }

    # garde-fou init incohérent
    if (!is.null(control$init)) {
      k <- length(summary(build_formula_from_rhs(rhs, nw = nw)))
      if (length(control$init) %in% c(1L, k-1L, k+1L)) control$init <- NULL
      ctrl <- do.call(ergm::control.ergm, .modify_list(ctrl_defaults, control))
    }

    fit <- erpm(f, eval_call = TRUE, verbose = verbose, estimate = estimate,
                eval.loglik = eval_loglik, control = ctrl, timeout = timeout)

    if (!inherits(fit, "ergm")) {
      .log_info("Fit non disponible ou interrompu.")
      return(list(engine=engine, mode="fit_failed", rhs=rhs, rhs_text=rhs_text,
                  constraints=constraints, constraints_text=constraints_text,
                  call=NULL, call_text=NULL, fit=fit, result=fit,
                  network=nw, partition=partition))
    }

    .summarize_fit(fit, nw, partition, verbose)
    return(list(engine=engine, mode="fit", rhs=rhs, rhs_text=rhs_text,
                constraints=constraints, constraints_text=constraints_text,
                call=NULL, call_text=NULL, fit=fit, result=fit,
                network=nw, partition=partition))
  }

  # ============================================================
  # Chargement des dépendances à l'exécution
  # ============================================================

  .ensure_pkgs(c("network", "ergm", "R.utils"))
  invisible(try(.ensure_pkgs(c("igraph", "intergraph", "RColorBrewer")), silent = TRUE))
  .root <- getwd()

  # Charger le wrapper avant usage pour build_bipartite_from_inputs()
  .source_if_found(.cands("R/erpm_wrapper.R"), required = TRUE)
  .source_if_found(.cands("scripts/local_source/settings.R"))
  .source_if_found(.cands("scripts/local_source/colors.R"))
  .source_if_found(.cands("scripts/local_source/logging.R"))
  .source_if_found(.cands("scripts/local_source/ergm_utils.R"))
  .source_if_found(.cands("R/functions_erpm_bip_network.R"))

  # ============================================================
  # Exports
  # ============================================================

  assign("coerce_to_bipartite",        coerce_to_bipartite,        envir = .GlobalEnv)
  assign("build_formula_from_rhs",     build_formula_from_rhs,     envir = .GlobalEnv)
  assign("build_effect_call",          build_effect_call,          envir = .GlobalEnv)
  assign("build_rhs_from_effects",     build_rhs_from_effects,     envir = .GlobalEnv)
  assign("launch_model",               launch_model,               envir = .GlobalEnv)
  assign(".__launcher_loaded",         TRUE,                       envir = .GlobalEnv)
}

# # ==============================================================================
# # Fichier : launcher.R
# # Fonction : launch_model()
# # Utilité : Lancer rapidement un test de stat/fit (summary | ergm | erpm) avec 1+ effets
# # ==============================================================================

# if (!exists(".__launcher_loaded", envir = .GlobalEnv)) {

#     # ============================================================
#     # Dépendances (packages + sources locales) — helpers
#     # ============================================================

#     #' Opérateur coalescent simple
#     #'
#     #' Retourne \code{a} si non-NULL, sinon \code{b}.
#     #'
#     #' @param a,b Objets R quelconques.
#     #' @return \code{a} si non-NULL, sinon \code{b}.
#     #' @keywords internal
#     `%||%` <- function(a, b) if (!is.null(a)) a else b

#     #' Générer des chemins candidats pour un fichier relatif
#     #'
#     #' Construit une liste de chemins potentiels (relatifs/absolus) où chercher
#     #' un fichier, en couvrant les répertoires usuels du projet.
#     #'
#     #' @param rel Chemin relatif (ex. \code{"R/erpm_wrapper.R"}).
#     #' @return Un vecteur de chemins candidats.
#     #' @keywords internal
#     .cands <- function(rel) unique(c(
#         rel,
#         file.path("scripts", "local_source", rel),
#         file.path("R", rel),
#         file.path(.root, rel),
#         file.path(.root, "scripts", "local_source", rel),
#         file.path(.root, "R", rel)
#     ))

#     #' Charger/installer des packages si besoin
#     #'
#     #' Vérifie la présence des packages, les installe au besoin, puis les charge
#     #' silencieusement. Retourne de façon invisible la liste des \code{require()}.
#     #'
#     #' @param pkgs Vecteur de caractères des noms de packages.
#     #' @return (invisible) Liste de résultats de \code{require()}.
#     #' @keywords internal
#     .ensure_pkgs <- function(pkgs) {
#         miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
#         if (length(miss)) {
#             message(sprintf("Installation des packages manquants : %s", paste(miss, collapse = ", ")))
#             install.packages(miss, repos = getOption("repos", "https://cloud.r-project.org"))
#         }
#         invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(require(p, character.only = TRUE))))
#     }

#     #' Fusionner deux listes (les clés de `b` priment)
#     #'
#     #' Convertit `a` et `b` en listes, puis remplace/ajoute dans `a` toutes les
#     #' entrées nommées présentes dans `b`. Les noms de `b` ont toujours priorité.
#     #'
#     #' @param a,b Objets convertibles en liste (list, environment, pairlist, etc.).
#     #' @return Une liste résultant de `a` modifiée par `b`.
#     #' @examples
#     #' .modify_list(list(x=1, y=2), list(y=9, z=3))
#     #' # $x 1 ; $y 9 ; $z 3
#     #' @keywords internal
#     .modify_list <- function(a, b) {
#     a <- as.list(a); b <- as.list(b)
#     for (nm in names(b)) a[[nm]] <- b[[nm]]
#     a
#     }

#     #' Source un fichier si présent (essaie plusieurs chemins)
#     #'
#     #' Tente \code{source()} sur la première occurrence existante parmi une liste de chemins.
#     #' Peut lever une erreur si \code{required = TRUE} et qu'aucun chemin n'existe.
#     #'
#     #' @param paths    Vecteur de chemins candidats.
#     #' @param local    Passé à \code{source(local = ...)}.
#     #' @param required Si TRUE, stoppe si aucun fichier n'est trouvé.
#     #' @return TRUE si un fichier a été sourcé, FALSE sinon.
#     #' @keywords internal
#     .source_if_found <- function(paths, local = FALSE, required = FALSE) {
#         for (p in paths) {
#             if (file.exists(p)) {
#                 tryCatch({
#                     source(p, local = local)
#                     return(invisible(TRUE))
#                 }, error = function(e) {
#                     message(sprintf("Erreur lors de source(\"%s\") : %s", p, conditionMessage(e)))
#                 })
#             }
#         }
#         if (required) stop(sprintf("Introuvable : %s", paste(paths, collapse = " | ")))
#         invisible(FALSE)
#     }

#     # ============================================================
#     # Helpers logging/format
#     # ============================================================

#     #' Résumé bref et typé d'un objet
#     #'
#     #' Produit un résumé compact selon le cas :
#     #' 1) \code{engine = "summary"} : capture \code{print(obj)} puis tronque.
#     #' 2) Objet \code{ergm} : \code{"ergm fit: logLik=<valeur> ; #coef=<k>"}.
#     #' 3) Objet de type \code{call} : \code{"call: <appel>"} via \code{.fmt_call()}.
#     #' 4) Autres : \code{"Objet de classe: <classes>"}.
#     #'
#     #' @param obj Objet à résumer.
#     #' @param engine Optionnel, \code{"summary"} pour forcer la capture/tri.
#     #' @return \code{character(1)} résumé textuel.
#     #' @keywords internal
#     .brief_result <- function(obj, engine = NULL) {
#         if (identical(engine, "summary")) {
#             co <- utils::capture.output(print(obj))
#             if (length(co) == 0L) {
#                 # Fallbacks si aucune sortie n'est capturée
#                 if (is.atomic(obj) && !is.null(names(obj))) {
#                     co <- sprintf("%s = %s", names(obj), as.character(obj))
#                 } else if (is.atomic(obj)) {
#                     co <- toString(obj)
#                 } else {
#                     co <- utils::capture.output(print.default(obj))
#                     if (length(co) == 0L) co <- "<aucune sortie capturée>"
#                 }
#             }
#             return(paste(.head_lines(co), collapse = "\n"))
#         }
#         if (inherits(obj, "ergm")) {
#             ll <- tryCatch(as.numeric(stats::logLik(obj)), error = function(e) NA_real_)
#             k  <- tryCatch(length(stats::coef(obj)),       error = function(e) NA_integer_)
#             return(sprintf("ergm fit: logLik=%.4f ; #coef=%s", ll, ifelse(is.na(k), "NA", k)))
#         }
#         if (is.call(obj)) return(paste("call:", .fmt_call(obj)))
#         sprintf("Objet de classe: %s", paste(class(obj), collapse = ","))
#     }

#     #' Formater un appel R sur une seule ligne
#     #'
#     #' Convertit une expression/forme/appel (\code{call}) en texte sur une seule ligne.
#     #' Si c’est une formule et que le LHS est un \code{network}, remplace par \code{nw}
#     #' pour éviter d’imprimer l’objet complet.
#     #'
#     #' @param x Expression/forme/appel R.
#     #' @return \code{character(1)} ligne formatée.
#     #' @keywords internal
#     .fmt_call <- function(x) {
#         if (is.call(x) && identical(x[[1]], as.name("~")) && length(x) >= 3L) {
#             lhs <- x[[2]]
#             if (inherits(lhs, "network")) x[[2]] <- as.name("nw")
#         }
#         out <- tryCatch(deparse(x, width.cutoff = 500L),
#                         error = function(e) sprintf("<erreur de deparse: %s>", conditionMessage(e)))
#         paste(out, collapse = "  ")
#     }

#     #' Formater un RHS (termes d'effets) sur une seule ligne
#     #'
#     #' @param rhs Appel R représentant le RHS combiné (avec des '+').
#     #' @return \code{character(1)} RHS compact.
#     #' @keywords internal
#     .fmt_rhs <- function(rhs) paste(deparse(rhs, width.cutoff = 1000L), collapse = " ")

#     #' Extraire la tête d'un vecteur de lignes avec indication de troncature
#     #'
#     #' @param x Vecteur de lignes.
#     #' @param n Nombre maximal de lignes à retourner.
#     #' @return Les \code{n} premières lignes (avec \code{"... (+K lignes)"} si troncature).
#     #' @keywords internal
#     .head_lines <- function(x, n = 8L) {
#         if (length(x) <= n) x else c(x[seq_len(n)], sprintf("... (+%d lignes)", length(x) - n))
#     }

#     #' Log d'information interne
#     #'
#     #' Écrit un message INFO via \code{log_msg()} si présent, sinon \code{message()}.
#     #'
#     #' @param txt Message à consigner.
#     #' @return \code{invisible(NULL)}.
#     #' @keywords internal
#     .log_info <- function(txt) {
#         if (exists("log_msg")) log_msg("INFO", txt) else message(txt)
#     }

#     #' Afficher un fit via summary_ergm_model() si dispo (sinon fallback)
#     #'
#     #' @param fit Objet \code{ergm}.
#     #' @param nw Réseau (optionnel pour le résumé enrichi).
#     #' @param partition Partition (optionnelle pour le résumé enrichi).
#     #' @param verbose Afficher sur console si TRUE (si pas de log).
#     #' @return \code{invisible(character)} message résumé.
#     #' @keywords internal
#     .summarize_fit <- function(fit, nw, partition, verbose = TRUE) {
#         if (exists("summary_ergm_model")) {
#             # Si log_msg() existe, summary_ergm_model écrira au log ; sinon on cat()
#             send_to_log <- exists("log_msg")
#             msg <- summary_ergm_model(fit, log = send_to_log, nw = nw, partition = partition)
#             if (isTRUE(verbose) && !send_to_log) cat(msg, "\n")
#             invisible(msg)
#         } else {
#             # Fallback : résumé bref
#             fit_brief <- .brief_result(fit, engine = "ergm")
#             .log_info(paste("Résultat (ergm):", fit_brief))
#             if (isTRUE(verbose)) cat(fit_brief, "\n")
#             invisible(fit_brief)
#         }
#     }

#     # ============================================================
#     # Normalisation des wrappers (Proj1 / B)
#     # ============================================================

#     #' Normaliser les wrappers \code{Proj1()} / \code{B()} sans \verb{~}
#     #'
#     #' Transforme récursivement un RHS d'effets pour que les appels \code{Proj1(x)}
#     #' ou \code{B(x)} deviennent \code{Proj1(~ x)} ou \code{B(~ x)} si l'argument
#     #' n'est pas déjà une formule. Ne modifie rien si la forme est correcte.
#     #'
#     #' @param expr Expression/call représentant un RHS (avec \code{+} possibles).
#     #' @return Expression normalisée.
#     #' @keywords internal
#     .normalize_wrappers <- function(expr) {
#         if (!is.language(expr)) return(expr)

#         fix_one <- function(call_sym, node) {
#             if (is.call(node) && identical(node[[1]], as.name(call_sym))) {
#                 if (length(node) >= 2L) {
#                     arg1 <- node[[2]]
#                     if (!(is.call(arg1) && identical(arg1[[1]], as.name("~")))) {
#                         node[[2]] <- as.call(list(as.name("~"), arg1))
#                     }
#                 }
#             }
#             node
#         }

#         if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
#             expr[[2]] <- .normalize_wrappers(expr[[2]])
#             expr[[3]] <- .normalize_wrappers(expr[[3]])
#             return(expr)
#         }

#         expr <- fix_one("Proj1", expr)
#         expr <- fix_one("B",     expr)
#         expr
#     }

#     # ============================================================
#     # Construction réseau / formule / effets
#     # ============================================================

#                 #' Construire un réseau bipartite depuis une partition
#                 #'
#                 #' Construit un \code{network} bipartite à partir d'une partition fournie,
#                 #' ou utilise \code{create_erpm_network()} si \code{partition} est \code{NULL}.
#                 #' Peut tracer le réseau et consigner des infos via \code{log_erpm_network()}.
#                 #'
#                 #' @param partition  Vecteur des groupes (entiers/factor) ou \code{NULL}.
#                 #' @param labels     Noms des objets (défaut : \code{LETTERS[...]}).
#                 #' @param attributes Liste d'attributs par objet.
#                 #' @param verbose    Affichages/logs.
#                 #' @param plot       Tracer le réseau/partition.
#                 #' @return \code{list(network = nw, partition = part)}.
#                 # build_bipartite_from_partition <- function(partition = NULL,
#                 #                                            labels    = NULL,
#                 #                                            attributes = list(),
#                 #                                            verbose   = TRUE,
#                 #                                            plot      = FALSE) {
#                 #     if (is.null(partition)) {
#                 #         if (exists("log_msg")) log_msg("INFO", "Aucune partition en entrée — utilisation de la partition de démonstration.")
#                 #         else if (isTRUE(verbose)) message("INFO: aucune partition en entrée — utilisation de la partition de démonstration.")
#                 #         res  <- create_erpm_network()
#                 #         nw   <- res$network
#                 #         part <- res$partition
#                 #     } else {
#                 #         stopifnot(is.atomic(partition), length(partition) >= 2)
#                 #         if (is.null(labels)) labels <- LETTERS[seq_along(partition)]
#                 #         stopifnot(length(labels) == length(partition))
#                 #         nw   <- partition_to_bipartite_network(labels, partition, attributes)
#                 #         part <- partition
#                 #     }

#                 #     if (exists("log_erpm_network")) {
#                 #         log_erpm_network(list(network = nw, partition = part), verbose = verbose)
#                 #     } else if (isTRUE(verbose)) {
#                 #         cat(sprintf("Réseau : %d sommets, %d arêtes\n",
#                 #                     network::network.size(nw), network::network.edgecount(nw)))
#                 #     }

#                 #     if (isTRUE(plot)) {
#                 #         try(plot_partition_clusters(nw), silent = TRUE)
#                 #     }

#                 #     list(network = nw, partition = part)
#                 # }

#     # Unifier l'entrée: partition | network biparti
#     coerce_to_bipartite <- function(nw = NULL,
#                                     partition = NULL,
#                                     nodes = NULL,
#                                     dyads = list(),
#                                     verbose = TRUE,
#                                     plot = FALSE) {
#         if (!is.null(nw)) {
#             if (!is_bipartite_network(nw))
#             stop("nw fourni n'est pas biparti (attribut %n% 'bipartite' manquant ou invalide).")
#             # Partition best-effort à partir des arêtes acteur→groupe
#             n <- nw %n% "bipartite"
#             vn <- network::network.vertex.names(nw)
#             actors <- vn[seq_len(n)]
#             groups <- vn[seq.int(n + 1L, n + max(1L, length(vn) - n))]
#             A <- network::as.sociomatrix(nw)#A <- sna::as.matrix.network.adjacency(nw)
#             part <- vapply(actors, function(a) {
#             g <- groups[which(A[a, groups] > 0L)]
#             if (length(g) == 0L) NA_integer_ else as.integer(sub("^G", "", g[1L]))
#             }, integer(1))
#             return(list(network = nw, partition = part))
#         }

#         # sinon construire depuis partition (+ attributs)
#         if (is.null(partition)) {
#             if (exists("log_msg")) log_msg("INFO","Aucune partition — réseau démo.")
#             res <- create_erpm_network()
#             nw  <- res$network; part <- res$partition
#             return(list(network = nw, partition = part))
#         }

#         built <- build_bipartite_from_inputs(partition = partition, nodes = nodes, dyads = dyads)
#         if (exists("log_erpm_network")) log_erpm_network(built, verbose = verbose)
#         if (isTRUE(plot)) try(plot_partition_clusters(built$network), silent = TRUE)
#         list(network = built$network, partition = built$partition)
#     }

#     #' Construire une formule \code{nw ~ <rhs>} en liant \code{nw} dans l'environnement
#     #'
#     #' @param rhs RHS (expression \code{call} ou chaîne à parser).
#     #' @param nw  Objet \code{network} à lier dans l'environnement de la formule.
#     #' @return Une \code{formula} \code{nw ~ <rhs>} dont l'environnement contient \code{nw}.
#     #' @keywords internal
#     build_formula_from_rhs <- function(rhs, nw = NULL) {
#         rhs_expr <- if (is.character(rhs)) parse(text = rhs)[[1]] else rhs
#         f <- as.formula(bquote(nw ~ .(rhs_expr)))
#         if (!is.null(nw)) {
#             environment(f) <- list2env(list(nw = nw), parent = parent.frame())
#         }
#         f
#     }

#     #' Construire un appel d'effet unique (brique élémentaire)
#     #'
#     #' @param effect      Nom de l'effet (ex. \code{"b2degrange"}) ou appel \code{call()}.
#     #' @param effect_args Liste nommée d'arguments si \code{effect} est un caractère.
#     #' @return Un objet \code{call} (ex. \code{b2degrange(from=2,to=3)}).
#     #' @keywords internal
#     build_effect_call <- function(effect, effect_args = list()) {
#         if (is.call(effect)) return(effect)
#         if (is.character(effect) && length(effect) == 1L) {
#             return(as.call(c(as.name(effect), effect_args)))
#         }
#         stop("`effect` doit être un symbole d'effet (caractère) ou un appel R (call).")
#     }

#     #' Construire un RHS à partir de 1+ effets
#     #'
#     #' Accepte : chaîne avec \code{'+'}, expression/call, liste de \code{call}s,
#     #' ou vecteur de noms d'effets avec dictionnaire d'arguments.
#     #'
#     #' @param effects      Description des effets (string/call/list/char vector).
#     #' @param effects_args Liste d’arguments (nommée par effet ou parallèle).
#     #' @return Un \code{call} représentant le RHS combiné par \code{'+'}.
#     #' @keywords internal
#     build_rhs_from_effects <- function(effects, effects_args = list()) {
#         if (is.character(effects) && length(effects) == 1L) {
#             return(parse(text = effects)[[1]])
#         }
#         if (is.language(effects)) {
#             return(effects)
#         }
#         if (is.list(effects) && length(effects) >= 1L && all(vapply(effects, is.call, logical(1)))) {
#             if (length(effects) == 1L) return(effects[[1L]])
#             return(Reduce(function(x, y) call("+", x, y), effects))
#         }
#         if (is.character(effects) && length(effects) >= 1L) {
#             calls <- vector("list", length(effects))
#             named <- length(names(effects_args)) > 0
#             for (i in seq_along(effects)) {
#                 eff <- effects[i]
#                 args_i <- if (named) effects_args[[eff]] %||% list() else effects_args[[i]] %||% list()
#                 calls[[i]] <- build_effect_call(eff, args_i)
#             }
#             if (length(calls) == 1L) return(calls[[1L]])
#             return(Reduce(function(x, y) call("+", x, y), calls))
#         }
#         stop("Format `effects` non reconnu. Utilisez chaîne/expression/liste d'appels ou vecteur de noms.")
#     }

#     # ============================================================
#     # Lancer un test/stat/fit (summary | ergm | erpm) avec 1+ effets
#     # ============================================================

#     #' Lancer un calcul (summary | ergm | erpm) avec 1+ effets
#     #'
#     #' Point d'entrée unique pour :
#     #' - calculer une statistique via \code{summary(nw ~ <RHS>)},
#     #' - fitter un modèle \pkg{ergm} (\code{engine="ergm"}),
#     #' - traduire puis fitter via le wrapper \code{erpm} (\code{engine="erpm"}).
#     #'
#     #' Journalise l'appel (moteur, RHS, contraintes, dry-run) et un résumé du résultat.
#     #'
#     #' @param engine \code{"summary"} | \code{"ergm"} | \code{"erpm"}.
#     #' @param effects Description des effets (string/call/list) ; voir \code{build_rhs_from_effects()}.
#     #' @param effects_args Liste d'arguments pour les effets (nommée ou parallèle).
#     #' @param partition Partition; si \code{NULL}, utilise \code{create_erpm_network()}.
#     #' @param nw Objet \code{network}; si \code{NULL}, construit depuis \code{partition}.
#     #' @param dry_run \code{TRUE} = retourne l'appel/la formule sans exécuter.
#     #' @param verbose Affichages utilisateur.
#     #' @param plot Trace le réseau/partition.
#     #' @param constraints Contrainte ERGM ; par défaut \code{~ b1part} si réseau biparti.
#     #' @param estimate Méthode d'estimation \code{c("CD","MCMLE")} (par défaut \code{"CD"} pour rapide).
#     #' @param eval_loglik \code{logical}. Si \code{TRUE}, évalue la log-vraisemblance (bridge sampling).
#     #'   Utile surtout avec \code{estimate="MCMLE"} ; coûteux.
#     #' @param control Liste passée à \code{control.ergm()} (ex. \code{list(MCMLE.maxit=3, MCMC.samplesize=1000)}).
#     #'   Les valeurs fournies écrasent les valeurs par défaut "rapides".
#     #' @param timeout Nombre de secondes pour interrompre proprement le fit (via \pkg{R.utils}).
#     #'
#     #' @return Une liste : \code{engine}, \code{mode}, \code{rhs} (+texte), \code{constraints} (+texte),
#     #'   \code{call}/\code{fit}/\code{result}, \code{network}, \code{partition}.
#     #' @examples
#     #' # 1) Statistiques rapides (summary)
#     #' launch_model("summary", effects="b2degree(d=2)")
#     #'
#     #' # 2) ERGM rapide (CD) avec contrainte bipartite auto
#     #' launch_model("ergm", effects="b2degree(d=2) + Proj1(nodematch('gender'))")
#     #'
#     #' # 3) ERPM → ERGM (renommage + encapsulations) en MCMLE light
#     #' launch_model("erpm",
#     #'   effects     = "groups(from=2,to=3) + squared_sizes(from=1,to=2)",
#     #'   estimate    = "MCMLE",
#     #'   control     = list(MCMLE.maxit=3, MCMC.samplesize=1000),
#     #'   eval_loglik = FALSE,
#     #'   timeout     = 60)
#     #' @export
#     launch_model <- function(   engine = c("summary", "ergm", "erpm"),
#                                 effects,
#                                 effects_args = list(),
#                                 partition = NULL,
#                                 nodes = NULL, 
#                                 dyads = list(),
#                                 nw = NULL,
#                                 dry_run = FALSE,
#                                 verbose = TRUE,
#                                 plot = FALSE,
#                                 constraints = NULL,
#                                 estimate = c("MLE", "MPLE", "CD"),
#                                 eval_loglik = FALSE,
#                                 control = list(),
#                                 timeout = NULL) {

#         estimate <- match.arg(estimate)
#         engine <- match.arg(engine)                                   # Valide/choisit le moteur

#                 # -- Réseau : construit à partir de la partition si besoin ----
#                 # if (is.null(nw)) {
#                 #     built <- build_bipartite_from_partition(                  # Crée un réseau biparti démo ou depuis `partition`
#                 #         partition = partition,
#                 #         labels    = NULL,
#                 #         attributes = list(),
#                 #         verbose   = verbose,
#                 #         plot      = plot
#                 #     )
#                 #     nw        <- built$network                                # Récupère le network créé
#                 #     partition <- built$partition                              # Récupère la partition associée
#                 # }

#         # -- Réseau : déduire/constuire depuis partition|nw + attributs ----
#         if (is.null(nw)) {
#         built <- coerce_to_bipartite(nw = NULL,
#                                     partition = partition,
#                                     nodes = nodes,
#                                     dyads = dyads,
#                                     verbose = verbose,
#                                     plot = plot)
#         nw        <- built$network
#         partition <- built$partition
#         } else {
#         built <- coerce_to_bipartite(nw = nw, partition = NULL, verbose = verbose, plot = FALSE)
#         # nw conservé tel quel, partition best-effort
#         partition <- built$partition
#         }

#         # -- RHS : combine 1+ effets, puis normalise Proj1/B ---------
#         rhs <- build_rhs_from_effects(effects, effects_args)          # Construit l'appel RHS (somme d’effets)
#         rhs <- .normalize_wrappers(rhs)                               # Corrige Proj1(x) -> Proj1(~ x), etc.
#         rhs_text <- .fmt_call(rhs)                                    # Version texte pour logs

#         # -- Contraintes par défaut si réseau biparti -----------------
#         if (is.null(constraints)) {
#             bip <- tryCatch(nw %n% "bipartite", error = function(e) NULL) # Lit l’attribut bipartite du network
#             if (!is.null(bip) && !is.na(bip))                           # Si biparti, applique contrainte standard
#                 constraints <- as.formula(~ b1part)
#         }
#         constraints_text <- if (is.null(constraints)) "NULL" else .fmt_call(constraints)

#         # -- Log d’entrée ---------------------------------------------
#         .log_info(sprintf("launch_model(engine=%s, dry_run=%s) | RHS=%s | constraints=%s",
#                           engine, dry_run, rhs_text, constraints_text))

#         # ========================= ENGINE: summary ====================
#         if (engine == "summary") {
#             f <- build_formula_from_rhs(rhs, nw = nw)                 # Formule `nw ~ <rhs>` (nw lié dans env)
#             call_text <- .fmt_call(f)                                 # Pour affichage/log

#             if (dry_run) {                                            # Mode dry-run : ne calcule pas
#                 .log_info(paste("Résultat (dry-run/summary): ", "summary(", call_text, ")"))
#                 return(list(
#                     engine = engine, mode = "dryrun",
#                     rhs = rhs, rhs_text = rhs_text,
#                     constraints = constraints, constraints_text = constraints_text,
#                     call = f, call_text = call_text,
#                     fit = NULL, result = f,
#                     network = nw, partition = partition
#                 ))
#             }

#             out <- summary(f)                                         # Calcule les stats de summary()
#             out_brief <- .brief_result(out, engine = "summary")       # Résumé compact pour logs
#             .log_info(paste("Résultat (summary):\n", out_brief))
#             # NB : on n'affiche pas systématiquement sur console pour éviter l'entrelacement

#             return(list(                                              # Retour standardisé
#                 engine = engine, mode = "summary",
#                 rhs = rhs, rhs_text = rhs_text,
#                 constraints = constraints, constraints_text = constraints_text,
#                 call = f, call_text = call_text,
#                 fit = NULL, result = out,
#                 network = nw, partition = partition
#             ))
#         }

#         # ========================== ENGINE: ergm ======================
#         if (engine == "ergm") {
#             f <- build_formula_from_rhs(rhs, nw = nw)                 # Formule `nw ~ <rhs>` pour ergm
#             call_eval <- as.call(list(as.name("ergm"), f, constraints = constraints))  # Appel ergm(...)
#             call_text <- .fmt_call(call_eval)

#             if (dry_run) {                                            # Mode dry-run : renvoie l’appel
#                 .log_info(paste("Résultat (dry-run/ergm):", call_text))
#                 return(list(
#                     engine = engine, mode = "dryrun",
#                     rhs = rhs, rhs_text = rhs_text,
#                     constraints = constraints, constraints_text = constraints_text,
#                     call = call_eval, call_text = call_text,
#                     fit = NULL, result = NULL,
#                     network = nw, partition = partition
#                 ))
#             }

#             # Contrôles par défaut (rapides) fusionnés avec ceux de l’utilisateur
#             ctrl_defaults <- list(
#                 MCMLE.maxit     = 5L,
#                 MCMC.samplesize = 2000L,
#                 MCMC.burnin     = 10000L,
#                 MCMC.interval   = 1000L
#             )
#             ctrl_list <- .modify_list(ctrl_defaults, control)
#             ctrl <- do.call(ergm::control.ergm, ctrl_list)

#             # --- Fit réel via ergm::ergm() avec garde-fous ---
#             .run_ergm <- function() {
#                 if (is.null(constraints)) {
#                     ergm::ergm(f, estimate = estimate, eval.loglik = eval_loglik, control = ctrl)
#                 } else {
#                     ergm::ergm(f, constraints = constraints, estimate = estimate, eval.loglik = eval_loglik, control = ctrl)
#                 }
#             }

#             fit <- if (is.null(timeout)) {
#                 .run_ergm()
#             } else {
#                 R.utils::withTimeout(.run_ergm(), timeout = as.numeric(timeout), onTimeout = "silent")
#             }

#             # fit <- if (is.null(constraints)) {                         # Fit réel via ergm::ergm()
#             #     ergm::ergm(f) 
#             # } else {
#             #     ergm::ergm(f, constraints = constraints)
#             # }

#             if (!inherits(fit, "ergm")) {
#                 .log_info("Fit non disponible ou interrompu (timeout/erreur). Aucun résumé.")
#                 return(list(
#                     engine=engine, mode="fit_failed",
#                     rhs=rhs, rhs_text=rhs_text,
#                     constraints=constraints, constraints_text=constraints_text,
#                     call=NULL, call_text=call_text,
#                     fit=fit, result=fit,
#                     network=nw, partition=partition
#                 ))
#             }
#             .summarize_fit(fit, nw = nw, partition = partition, verbose = verbose)  # Résumé enrichi si dispo

#             return(list(                                              # Retour standardisé
#                 engine = engine, mode = "fit",
#                 rhs = rhs, rhs_text = rhs_text,
#                 constraints = constraints, constraints_text = constraints_text,
#                 call = NULL, call_text = call_text,
#                 fit = fit, result = fit,
#                 network = nw, partition = partition
#             ))
#         }

#         # =========================== ENGINE: erpm =====================
#         if (!exists("erpm"))                                          # Sécurité : wrapper requis
#             stop("La fonction `erpm()` est introuvable. Sourcez R/erpm_wrapper.R.")

#         f <- build_formula_from_rhs(rhs, nw = nw)                     # Formule pour wrapper erpm

#         if (dry_run) {                                                # Dry-run : ne fit pas
#             call_erpm <- erpm(  f,                      # Demande l’appel traduit sans évaluer
#                                 eval_call   = FALSE,
#                                 verbose     = verbose,
#                                 estimate    = estimate
#                             )  
#             call_text <- .fmt_call(call_erpm)
#             .log_info(paste("Résultat (dry-run/erpm):", call_text))
#             return(list(
#                 engine = engine, mode = "dryrun",
#                 rhs = rhs, rhs_text = rhs_text,
#                 constraints = constraints, constraints_text = constraints_text,
#                 call = call_erpm, call_text = call_text,
#                 fit = NULL, result = NULL,
#                 network = nw, partition = partition
#             ))
#         }

#         # fit <- erpm(f, eval_call = TRUE, verbose = verbose)           # Fit réel via erpm (traduit → ergm)
        
#         ctrl_defaults <- list(
#                 MCMLE.maxit     = 5L,
#                 MCMC.samplesize = 2000L,
#                 MCMC.burnin     = 10000L,
#                 MCMC.interval   = 1000L
#             )
#         ctrl_list <- .modify_list(ctrl_defaults, control)

#         # juste avant do.call(ergm::control.ergm, ctrl_list)
#         if (!is.null(control$init)) {
#             k <- length(summary(build_formula_from_rhs(rhs, nw = nw)))
#             if (length(control$init) %in% c(1L, k-1L, k+1L)) control$init <- NULL
#         }
#         ctrl <- do.call(ergm::control.ergm, ctrl_list)
#         fit <- erpm(
#             f,
#             eval_call   = TRUE,
#             verbose     = verbose,
#             estimate    = estimate,
#             eval.loglik = eval_loglik,
#             control     = ctrl,
#             timeout     = timeout
#         )
        
#         if (!inherits(fit, "ergm")) {
#             .log_info("Fit non disponible ou interrompu (timeout/erreur). Aucun résumé.")
#             return(list(
#                 engine=engine, mode="fit_failed",
#                 rhs=rhs, rhs_text=rhs_text,
#                 constraints=constraints, constraints_text=constraints_text,
#                 call=NULL, call_text=call_text,
#                 fit=fit, result=fit,
#                 network=nw, partition=partition
#             ))
#         }

#         .summarize_fit(fit, nw = nw, partition = partition, verbose = verbose)  # Résumé enrichi si dispo

#         return(list(                                                  # Retour standardisé
#             engine = engine, mode = "fit",
#             rhs = rhs, rhs_text = rhs_text,
#             constraints = constraints, constraints_text = constraints_text,
#             call = NULL, call_text = NULL,
#             fit = fit, result = fit,
#             network = nw, partition = partition
#         ))
#     }

#     # ============================================================
#     # Chargement des dépendances à l'exécution
#     # ============================================================

#     .ensure_pkgs(c("network", "ergm", "R.utils"))
#     invisible(try(.ensure_pkgs(c("igraph", "intergraph", "RColorBrewer")), silent = TRUE))
#     .root <- getwd()

#     .source_if_found(.cands("scripts/local_source/settings.R"))
#     .source_if_found(.cands("scripts/local_source/colors.R"))
#     .source_if_found(.cands("scripts/local_source/logging.R"))
#     .source_if_found(.cands("scripts/local_source/ergm_utils.R"))
#     .source_if_found(.cands("R/functions_erpm_bip_network.R"))
#     .source_if_found(.cands("R/erpm_wrapper.R"))

#     # ============================================================
#     # Exports
#     # ============================================================

#     assign("coerce_to_bipartite", coerce_to_bipartite, envir = .GlobalEnv)
#     assign("build_formula_from_rhs",         build_formula_from_rhs,         envir = .GlobalEnv)
#     assign("build_effect_call",              build_effect_call,              envir = .GlobalEnv)
#     assign("build_rhs_from_effects",         build_rhs_from_effects,         envir = .GlobalEnv)
#     assign("launch_model",                   launch_model,                   envir = .GlobalEnv)
#     assign(".__launcher_loaded",             TRUE,                           envir = .GlobalEnv)
# }
