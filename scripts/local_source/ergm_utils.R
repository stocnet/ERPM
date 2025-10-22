# ==============================================================================
# Fichier : ergm_utils.R
# Fonction : summary_ergm_model()
# Utilité : Fournit un résumé automatique et loggué pour un objet de classe ERGM
# ==============================================================================

if(!exists(".__ergm_utils_loaded", envir = .GlobalEnv)) {

  #' Résumé automatique d'un modèle ERGM
  #'
  #' Cette fonction affiche un résumé clair d'un objet de classe \code{ergm},
  #' incluant la log-vraisemblance, les coefficients, la structure du réseau
  #' (si fournie) et la répartition de la partition (si fournie).
  #' Les erreurs de calcul ou coefficients non-finis sont gérées proprement.
  #'
  #' @param ergm_model Objet de classe \code{ergm}.
  #' @param log Logical : si TRUE, envoie le résumé vers \code{log_msg} sinon vers la console.
  #' @param nw Optionnel : objet réseau (de type \code{network}) à résumer.
  #' @param partition Optionnel : vecteur de partition à afficher.
  #' @return Invisible, affiche le résumé.
  #' @examples
  #' \dontrun{
  #' nw <- network::network.initialize(6, directed = FALSE)
  #' network::add.edges(nw, c(1,2,3), c(2,3,4))
  #' fit <- ergm::ergm(nw ~ edges)
  #' summary_ergm_model(fit, nw = nw, partition = c(1,1,2,2,3,3))
  #' }
  summary_ergm_model <- function(ergm_model, log = TRUE, nw = NULL, partition = NULL) {

    # --- Vérifications de base ---
    if (!inherits(ergm_model, "ergm")) {
      warning("L'objet fourni n'est pas de classe 'ergm'.")
      return(invisible(NULL))
    }

    # --- Récupération des infos du modèle ---
    coefs <- tryCatch(stats::coef(ergm_model), error = function(e) NULL)
    ll    <- tryCatch(as.numeric(stats::logLik(ergm_model)), error = function(e) NA)
    summ  <- tryCatch(utils::capture.output(print(summary(ergm_model))),
                      error = function(e) "Résumé non disponible")

    # --- Informations réseau (si fournies) ---
    net_info <- NULL
    if (!is.null(nw)) {
      n <- tryCatch(network::network.size(nw), error = function(e) NA_integer_)
      m <- tryCatch(network::network.edgecount(nw), error = function(e) NA_integer_)
      bip <- tryCatch(nw %n% "bipartite", error = function(e) NULL)
      if (is.null(bip)) bip <- NA_integer_

      if (!is.na(bip) && !is.na(n)) {
        net_info <- sprintf("RÉSEAU : %d sommets, %d arêtes, biparti = TRUE (|U|=%d, |V|=%d)",
                            n, m, bip, n - bip)
      } else {
        net_info <- sprintf("RÉSEAU : %d sommets, %d arêtes, biparti = FALSE/NA", n, m)
      }
    }

    # --- Informations sur la partition (si fournie) ---
    part_info <- NULL
    if (!is.null(partition)) {
      ng  <- length(unique(stats::na.omit(partition)))
      tab <- tryCatch(table(partition, useNA = "ifany"), error = function(e) NULL)
      dist <- if (!is.null(tab)) paste(names(tab), as.integer(tab), sep = ":", collapse = ", ") else "—"
      part_info <- sprintf("PARTITION : %d groupes | répartition [%s]", ng, dist)
    }

    # --- Construction du message final ---
    msg <- paste0(
      "\n\n===== RÉSUMÉ DU MODÈLE ERGM =====\n",
      if (!is.null(net_info)) paste0(net_info, "\n") else "",
      if (!is.null(part_info)) paste0(part_info, "\n") else "",
      "Log-likelihood = ", ifelse(is.na(ll), "NA", round(ll, 4)), "\n",
      "Coefficients  = ", if(!is.null(coefs)) paste(round(coefs, 4), collapse = ", ") else "—", "\n",
      paste(summ, collapse = "\n")
    )

    # --- Sortie ou log ---
    if (isTRUE(log) && exists("log_msg")) {
      log_msg("INFO", msg)
    } else {
      cat(msg, "\n")
    }

    invisible(msg)
  }

  # --- Attribution à l'environnement global ---
  assign("summary_ergm_model", summary_ergm_model, envir = .GlobalEnv)
  assign(".__ergm_utils_loaded", TRUE, envir = .GlobalEnv)
}
