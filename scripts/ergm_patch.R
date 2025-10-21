# ==============================================================================
# Fichier : ergm_patch.R
# Fonction : activer ou désactiver le patch sur base::replace()
# Utilité : contourner le bug "type 'builtin' d'indice incorrect" dans ergm.fit()
# ==============================================================================

if(!exists(".__ergm_patch_loaded", envir = .GlobalEnv)){

  #' Active un patch pour contourner un bug dans ergm.fit()
  #'
  #' Cette fonction trace et redéfinit temporairement la fonction \code{base::replace()}
  #' afin de contourner le bug "type 'builtin' d'indice incorrect" présent dans certaines
  #' versions de \pkg{ergm}. Le patch modifie le comportement de \code{replace()} pour
  #' traiter correctement les fonctions passées en argument.
  #'
  #' @details
  #' - Supprime tout trace existant sur \code{ergm.fit} ou \code{replace()}.
  #' - Active le tracing global.
  #' - Ajoute un trace sur \code{base::replace()} avec une logique de correction.
  #' - Vérifie le bon fonctionnement via un test minimal.
  #'
  #' @return Aucun. Affiche un message de confirmation.
  #' @examples
  #' ergm_patch_enable()
  #' @export
  ergm_patch_enable <- function() {
    # Nettoyage des traces précédentes
    if ("ergm" %in% loadedNamespaces()) {
      try(untrace("ergm.fit", where = asNamespace("ergm")), silent = TRUE)
    }
    try(untrace("replace", where = baseenv()), silent = TRUE)

    # Activation du tracing global
    tracingState(on = TRUE)

    # Patch via trace sur base::replace : correction des appels erronés
    trace(
      what  = "replace",
      where = baseenv(),
      tracer = quote({
        if (is.function(list)) list <- list(x)  # si un test de type is.na est passé, le corriger
      }),
      print = FALSE
    )

    # Vérification de bon fonctionnement du patch
    stopifnot(identical(replace(c(1, NA, 3), is.na, 0), c(1, 0, 3)))

    message("[ergm_patch_enable] patch appliqué avec succès.")
  }

  #' Désactive le patch sur base::replace()
  #'
  #' Supprime le traçage ajouté par \code{ergm_patch_enable()} et désactive le tracing global.
  #' @return Aucun. Affiche un message de confirmation.
  #' @examples
  #' ergm_patch_disable()
  #' @export
  ergm_patch_disable <- function() {
    try(untrace("replace", where = baseenv()), silent = TRUE)  # retire le trace
    tracingState(on = FALSE)                                   # désactive le tracing global
    message("[ergm_patch_disable] patch retiré et tracing désactivé.")
  }

  #' Auto-test du patch ergm
  #'
  #' Cette fonction exécute une série de tests automatiques sur différents types
  #' de réseaux (unipartite, bipartite, avec attributs, etc.) afin de vérifier le bon
  #' fonctionnement du patch sur \code{base::replace()}.
  #'
  #' @param run_unipartite Logique : exécuter le test uniparti simple.
  #' @param run_bipartite Logique : exécuter les tests bipartis.
  #' @param run_attributes Logique : inclure les attributs de nœuds dans un test uniparti.
  #' @param run_bip2star Logique : exécuter le test sur \code{b2star(2)} biparti.
  #' @param run_diagnostics Logique : générer des diagnostics MCMC si TRUE.
  #' @param seed Entier : graine aléatoire pour la reproductibilité.
  #' @param verbose Logique : afficher les messages d'exécution.
  #' @param patch Logique : activer automatiquement le patch avant les tests.
  #'
  #' @return Une liste contenant les résultats et éventuelles erreurs des tests.
  #' @examples
  #' ergm_patch_selftest(run_unipartite = TRUE, run_bipartite = TRUE)
  #' @export
  ergm_patch_selftest <- function(
                                  run_unipartite   = TRUE,
                                  run_bipartite    = TRUE,
                                  run_attributes   = TRUE,
                                  run_bip2star     = TRUE,
                                  run_diagnostics  = FALSE,
                                  seed             = 123,
                                  verbose          = TRUE,
                                  patch            = TRUE
                                  ){
    # Vérifie que les packages nécessaires sont présents
    stopifnot(requireNamespace("network", quietly = TRUE),
              requireNamespace("ergm", quietly = TRUE))

    # Active le patch si demandé
    if(patch) ergm_patch_enable()

    say <- function(...) if (isTRUE(verbose)) cat(...)

    results <- list()
    set.seed(seed)  # pour la reproductibilité des tirages aléatoires

    # --- utilitaire d'exécution et de capture propre ---------------------------
    # Cette fonction exécute un bloc expr, capture erreurs et temps d'exécution
    .run_case <- function(name, expr){
      t0 <- proc.time()
      out <- tryCatch(
        {
          val <- force(expr)
          list(ok = TRUE, value = val, error = NULL,
               time_sec = as.numeric((proc.time()-t0)["elapsed"]))
        },
        error = function(e){
          list(ok = FALSE, value = NULL, error = conditionMessage(e),
               time_sec = as.numeric((proc.time()-t0)["elapsed"]))
        }
      )
      results[[name]] <<- out  # stocke le résultat dans l'environnement parent
      
      # Affichage du résultat du test
      if (verbose) {
        if (out$ok)
          say(sprintf("✅ %-18s (%.2fs)\n", paste0(name, ":"), out$time_sec))
        else
          say(sprintf("❌ %-18s (%.2fs) -> %s\n", paste0(name, ":"), out$time_sec, out$error))
      }
      invisible(out)
    }

    # ------------------------- CAS 1 : UNIPARTI --------------------------------
    if (run_unipartite){
      .run_case("unipartite_edges+triangles", {
        # Crée un petit réseau non orienté à 8 nœuds
        nw <- network::network.initialize(8, directed = FALSE)
        # Ajoute une chaîne d’arêtes
        network::add.edges(nw, c(1,2,3,4,5), c(2,3,4,5,6))
        # Ajuste un modèle ERGM simple
        fit <- ergm::ergm(nw ~ edges + triangles, verbose = if (verbose) 2 else 0)
        # Retourne les résultats principaux
        list(
          coef  = stats::coef(fit),
          ll    = as.numeric(stats::logLik(fit)),
          summ  = utils::capture.output(print(summary(fit)))
        )
      })
    }

    # ------------------------- CAS 2 : ATTRIBUTS -------------------------------
    if (run_attributes){
      .run_case("unipartite_with_nodematch", {
        # Crée un graphe uniparti de 10 nœuds
        nw <- network::network.initialize(10, directed = FALSE)
        pairs <- utils::combn(1:10, 2)
        keep  <- which(stats::runif(ncol(pairs)) < 0.25)
        if (length(keep)) network::add.edges(nw, pairs[1,keep], pairs[2,keep])

        # Ajoute un attribut de groupe
        network::set.vertex.attribute(nw, "grp", sample(c("A","B"), network::network.size(nw), TRUE))
        # Ajuste un modèle avec un terme nodematch
        fit <- ergm::ergm(nw ~ edges + nodematch("grp"), verbose = if (verbose) 2 else 0)
        list(coef = stats::coef(fit), ll = as.numeric(stats::logLik(fit)))
      })
    }

    # ------------------------- CAS 3 : BIPARTI sans patch ----------------------
    if (run_bipartite){
      .run_case("bipartite_b2degrange_NO_PATCH", {
        # Réseau biparti 3x3 (objets 1..3, groupes 4..6)
        nw <- network::network.initialize(6, directed = FALSE, bipartite = 3)
        network::add.edges(nw, c(1,1,2,3), c(4,5,5,6))
        # Pas de patch ici -> on s’attend à une erreur sur certaines versions d’ergm
        fit <- ergm::ergm(nw ~ b2degrange(from = 2, to = 3), verbose = if (verbose) 2 else 0)
        list(coef = stats::coef(fit), ll = as.numeric(stats::logLik(fit)))
      })
    }

    # ------------------------- CAS 4 : BIPARTI avec patch ----------------------
    if (run_bipartite){
      .run_case("bipartite_b2degrange_WITH_PATCH", {
        # Active le patch temporairement
        if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
        on.exit({
          if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()
        }, add = TRUE)

        # Même réseau biparti, cette fois avec patch actif
        nw <- network::network.initialize(6, directed = FALSE, bipartite = 3)
        network::add.edges(nw, c(1,1,2,3), c(4,5,5,6))
        fit <- ergm::ergm(nw ~ b2degrange(from = 2, to = 3), verbose = if (verbose) 2 else 0)

        # Diagnostics facultatifs
        diag <- NULL
        if (isTRUE(run_diagnostics)) {
          diag <- tryCatch(utils::capture.output(ergm::mcmc.diagnostics(fit)),
                           error = function(e) paste("diag_error:", conditionMessage(e)))
        }

        list(coef = stats::coef(fit),
             ll = as.numeric(stats::logLik(fit)),
             diagnostics = diag)
      })
    }

    # ------------------------- CAS 5 : BIPARTI b2star(2) -----------------------
    if (run_bip2star){
      .run_case("bipartite_b2star2_WITH_PATCH", {
        if (exists("ergm_patch_enable", mode = "function")) ergm_patch_enable()
        on.exit({
          if (exists("ergm_patch_disable", mode = "function")) ergm_patch_disable()
        }, add = TRUE)

        # Crée un réseau biparti 4x4 et ajoute quelques liens
        nw <- network::network.initialize(8, directed = FALSE, bipartite = 4)
        tails <- c(1,1,2,3,4)
        heads <- c(5,6,6,7,8)
        network::add.edges(nw, tails, heads)

        # Ajuste un modèle ERGM b2star(2)
        fit <- ergm::ergm(nw ~ b2star(2), verbose = if (verbose) 2 else 0)
        list(coef = stats::coef(fit), ll = as.numeric(stats::logLik(fit)))
      })
    }

    # ------------------------- Résumé final ------------------------------------
    if (verbose) {
      say("\n===== RÉSUMÉ DES TESTS =====\n")
      for (nm in names(results)) {
        r <- results[[nm]]
        if (!isTRUE(r$ok)) {
          say(sprintf("• %-28s : ❌ %s\n", nm, r$error))
        } else {
          co <- r$value$coef
          ll <- r$value$ll
          say(sprintf("• %-28s : ✅ ll=%.4f; coef=%s\n",
                      nm, ifelse(is.null(ll), NaN, ll),
                      if (is.null(co)) "—" else paste(round(co, 4), collapse = ", ")))
        }
      }
    }

    invisible(results)
  }

  # Attribution des fonctions à l’environnement global
  assign("ergm_patch_enable",    ergm_patch_enable,    envir = .GlobalEnv)
  assign("ergm_patch_disable",   ergm_patch_disable,   envir = .GlobalEnv)
  assign("ergm_patch_selftest",  ergm_patch_selftest,  envir = .GlobalEnv)
  assign(".__ergm_patch_loaded", TRUE,                 envir = .GlobalEnv)
}
