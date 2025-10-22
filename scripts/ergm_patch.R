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
  ergm_patch_enable <- function(verbose = TRUE) {
    if(!ERGM_PATCHED){
      # Nettoyage des traces précédentes
      if ("ergm" %in% loadedNamespaces()) {
        suppressWarnings(suppressMessages(try(untrace("ergm.fit", where = asNamespace("ergm")), silent = TRUE)))
      }
      suppressWarnings(suppressMessages(try(untrace("replace", where = baseenv()), silent = TRUE)))

      tracingState(on = TRUE)

      # Patch via trace sur base::replace : correction des appels erronés
      suppressWarnings(suppressMessages(
        trace(
          what  = "replace",
          where = baseenv(),
          tracer = quote({
            if (is.function(list)) list <- list(x)
          }),
          print = FALSE
        )
      ))

      stopifnot(identical(replace(c(1, NA, 3), is.na, 0), c(1, 0, 3)))

      assign("ERGM_PATCHED", TRUE, envir = .GlobalEnv)
      if(isTRUE(verbose)) message("[ergm_patch_enable] patch appliqué avec succès. Tracing sur baseenv::replace() activé.")
    }
  }

  #' Désactive le patch sur base::replace()
  #'
  #' Supprime le traçage ajouté par \code{ergm_patch_enable()} et désactive le tracing global.
  #' @return Aucun. Affiche un message de confirmation.
  #' @examples
  #' ergm_patch_disable()
  #' @export
  ergm_patch_disable <- function(verbose = TRUE) {
  if(ERGM_PATCHED){
    suppressWarnings(suppressMessages(try(untrace("replace", where = baseenv()), silent = TRUE)))
    tracingState(on = FALSE)
    
    assign("ERGM_PATCHED", FALSE, envir = .GlobalEnv)
    if(isTRUE(verbose)) message("[ergm_patch_disable] patch retiré avec succès. Tracing sur baseenv::replace() desactivé..")
  }
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
    stopifnot(requireNamespace("network", quietly = TRUE),
              requireNamespace("ergm", quietly = TRUE))

    say <- function(...) if (isTRUE(verbose)) cat(...)

    say("\n\n===================== SELF-TEST PATCH ERGM =================================\n")

    # Active le patch une seule fois au début si demandé
    if(patch) ergm_patch_enable(verbose = verbose)
    on.exit({
      # Désactive le patch à la fin du self-test
      if(patch) ergm_patch_disable(verbose = verbose)
      say("\n===================== FIN SELF-TEST PATCH ERGM =================================\n\n")
    }, add = TRUE)

    results <- list()
    set.seed(seed)

    .run_case <- function(name, expr){
      t0 <- proc.time()
      out <- tryCatch(
        {
          val <- NULL
          if(verbose){
            # laisse ergm afficher normalement
            val <- force(expr)
          } else {
            # capture toute la sortie si verbose = FALSE
            val <- suppressWarnings(suppressMessages(
              capture.output(val <- force(expr), type = "output")
            ))
          }
          list(ok = TRUE, value = val, error = NULL,
              time_sec = as.numeric((proc.time()-t0)["elapsed"]))
        },
        error = function(e){
          list(ok = FALSE, value = NULL, error = conditionMessage(e),
              time_sec = as.numeric((proc.time()-t0)["elapsed"]))
        }
      )
      results[[name]] <<- out
      if(verbose){
        if(out$ok) say(sprintf("✅ %-28s (%.2fs)\n", paste0(name, ":"), out$time_sec))
        else say(sprintf("❌ %-28s (%.2fs) -> %s\n", paste0(name, ":"), out$time_sec, out$error))
      }
      invisible(out)
    }

    # ------------------------- CAS 1 : UNIPARTI -------------------------------
    if(run_unipartite){
      .run_case("unipartite_edges+triangles", {
        nw <- network::network.initialize(8, directed = FALSE)
        network::add.edges(nw, c(1,2,3,4,5), c(2,3,4,5,6))
        fit <- if(verbose){ 
                  ergm::ergm(nw ~ edges + triangles, verbose = 2)
                } else {
                  suppressWarnings(
                    suppressMessages(
                      ergm::ergm(nw ~ edges + triangles, verbose = 0)
                    )
                  )
                }
        list(coef = stats::coef(fit), ll = as.numeric(stats::logLik(fit)),
            summ = utils::capture.output(print(summary(fit))))
      })
    }

    # ------------------------- CAS 2 : ATTRIBUTS -------------------------------
    if(run_attributes){
      .run_case("unipartite_with_nodematch", {
        nw <- network::network.initialize(10, directed = FALSE)
        pairs <- utils::combn(1:10, 2)
        keep <- which(stats::runif(ncol(pairs)) < 0.25)
        if(length(keep)) network::add.edges(nw, pairs[1,keep], pairs[2,keep])
        network::set.vertex.attribute(nw, "grp", sample(c("A","B"), network::network.size(nw), TRUE))
        fit <- if(verbose){ 
                  ergm::ergm(nw ~ edges + nodematch("grp"), verbose = 2)
                } else {
                  suppressWarnings(
                    suppressMessages(
                      ergm::ergm(nw ~ edges + nodematch("grp"), verbose = 0)
                    )
                  )
                }
        list(coef = stats::coef(fit), ll = as.numeric(stats::logLik(fit)))
      })
    }


    # ------------------------- CAS 4 : BIPARTI ----------------------
    if(run_bipartite){
      .run_case("bipartite_b2degrange_WITH_PATCH", {
        nw <- network::network.initialize(6, directed = FALSE, bipartite = 3)
        network::add.edges(nw, c(1,1,2,3), c(4,5,5,6))
        fit <- if(verbose){ 
                  ergm::ergm(nw ~ b2degrange(from = 2, to = 3), verbose = 2)
                } else {
                  suppressWarnings(
                    suppressMessages(
                      ergm::ergm(nw ~ b2degrange(from = 2, to = 3), verbose = 0)
                    )
                  )
                }
        diag <- NULL
        if(isTRUE(run_diagnostics)) diag <- tryCatch(utils::capture.output(ergm::mcmc.diagnostics(fit)), error=function(e) paste("diag_error:", conditionMessage(e)))
        list(coef = stats::coef(fit), ll = as.numeric(stats::logLik(fit)), diagnostics = diag)
      })
    }

    # ------------------------- CAS 5 : BIPARTI b2star(2) -----------------------
    if(run_bip2star){
      .run_case("bipartite_b2star2_WITH_PATCH", {
        nw <- network::network.initialize(8, directed = FALSE, bipartite = 4)
        tails <- c(1,1,2,3,4)
        heads <- c(5,6,6,7,8)
        network::add.edges(nw, tails, heads)
        fit <- if(verbose){ 
                  ergm::ergm(nw ~ b2star(2), verbose = 2)
                } else {
                  suppressWarnings(
                    suppressMessages(
                      ergm::ergm(nw ~ b2star(2), verbose = 0)
                    )
                  )
                }
        list(coef = stats::coef(fit), ll = as.numeric(stats::logLik(fit)))
      })
    }

    # ------------------------- Résumé final ------------------------------------
    say("\n===== RÉSUMÉ DES TESTS =====\n")
    for(nm in names(results)){
      r <- results[[nm]]
      if(!isTRUE(r$ok)){
        say(sprintf("• %-28s : ❌ %s\n", nm, r$error))
      } else if(!is.list(r$value)){
        # si r$value n'est pas une liste, on l'affiche brute
        say(sprintf("• %-28s : ⚠ valeur inattendue : %s\n", nm, paste(r$value, collapse=", ")))
      } else {
        co <- r$value$coef
        ll <- r$value$ll
        say(sprintf("• %-28s : ✅ ll=%.4f; coef=%s\n", nm, ifelse(is.null(ll), NaN, ll),
                    if(is.null(co)) "—" else paste(round(co, 4), collapse=", ")))
      }
    }
    

    invisible(results)
  }


  # Attribution des fonctions à l’environnement global
  assign("ERGM_PATCHED",          FALSE,                envir = .GlobalEnv)
  assign("ergm_patch_enable",     ergm_patch_enable,    envir = .GlobalEnv)
  assign("ergm_patch_disable",    ergm_patch_disable,   envir = .GlobalEnv)
  assign("ergm_patch_selftest",   ergm_patch_selftest,  envir = .GlobalEnv)
  assign(".__ergm_patch_loaded",  TRUE,                 envir = .GlobalEnv)
}
