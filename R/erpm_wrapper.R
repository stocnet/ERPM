# ==============================================================================
# Fichier : erpm_wrapper.R
# Fonctions : .normalize_groups_args(), .split_sum_terms(), .translate_one_term(), erpm()
# Utilité : Traduire des formules ERPM en appels {ergm} (renommage + encapsulations)
# ==============================================================================

if (!exists(".__erpm_wrapper_loaded", envir = .GlobalEnv)) {

  # ============================================================================
  # Helpers internes
  # ============================================================================

  .oneline <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  .tight   <- function(s) gsub("\\s+", "", s)

  is_bipartite_network <- function(x) {
    inherits(x, "network") && !is.null(x %n% "bipartite") && is.finite(x %n% "bipartite")
  }

  .check_nodes_df <- function(nodes) {
    stopifnot(is.data.frame(nodes), "label" %in% names(nodes))
    if (anyDuplicated(nodes$label)) stop("nodes$label contient des doublons.")
    invisible(TRUE)
  }

  .check_dyads <- function(dyads, n, labels) {
    if (length(dyads) == 0L) return(invisible(TRUE))
    stopifnot(is.list(dyads))
    for (nm in names(dyads)) {
      M <- dyads[[nm]]
      stopifnot(is.matrix(M), nrow(M) == n, ncol(M) == n)
      # Optionnel mais utile : aligner par noms si fournis
      if (!is.null(rownames(M)) && !is.null(colnames(M))) {
        if (!identical(rownames(M), labels) || !identical(colnames(M), labels))
          stop(sprintf("dyads['%s']: row/colnames doivent matcher nodes$label.", nm))
      }
    }
    invisible(TRUE)
  }


  #' Normaliser les arguments de groups(...) vers b2degrange(from, to)
  #'
  #' Convertit les variantes d'appel ERPM `groups`, `groups(taille)`,
  #' et `groups(from=.., to=..)` en une paire (from, to) prête pour
  #' `b2degrange`. Les conventions sont :
  #' - `groups`            → [1, Inf)
  #' - `groups(taille)`    → [taille, taille+1)
  #' - `groups(from,to)`   → [from, to)
  #'
  #' @param args_list Liste d’arguments capturés depuis l’appel `groups(...)`.
  #' @return list(from=?, to=?), où `to` peut être `Inf` ou un appel `s+1`.
  #' @examples
  #' .normalize_groups_args(list())                           # 1, Inf
  #' .normalize_groups_args(list(2))                          # 2, 4 (to = 3 en sortie finale)
  #' .normalize_groups_args(list(from=2, to=5))               # 2, 5
  #' @keywords internal
  .normalize_groups_args <- function(args_list) {
    # Cas 0 : groups -> [1, Inf)
    if (length(args_list) == 0L)
      return(list(from = 1, to = quote(Inf)))

    nm <- names(args_list)

    # Cas 1 : groups(taille) -> [taille, taille+1)
    if (length(args_list) == 1L && (is.null(nm) || isTRUE(nm[1L] == ""))) {
      s <- args_list[[1L]]
      # Si taille est un numérique scalaire fini, calculer s+1 maintenant pour l’affichage propre
      to_val <- if (is.numeric(s) && length(s) == 1L && is.finite(s)) s + 1 else call("+", s, 1)
      return(list(from = s, to = to_val))
    }

    # Cas 2 : groups(from=..., to=...)
    if (!is.null(nm) && all(c("from","to") %in% nm)) {
      return(list(from = args_list[["from"]], to = args_list[["to"]]))
    }

    # Sinon : signature invalide
    stop("groups(): utilisez soit `groups`, soit `groups(taille)`, soit `groups(from=.., to=..)`.")  
  }

  #' Normaliser les arguments de cliques(...) vers (k, normalized)
  #'
  #' Conventions/par défauts :
  #' - `cliques` ou `cliques()`                 → k = 2, normalized = FALSE
  #' - `cliques(normalized=TRUE|FALSE)`         → k = 2, normalized = <bool>
  #' - `cliques(clique_size=k)`                 → k = <k>, normalized = FALSE
  #' - `cliques(k)` (un seul argument nu)       → k = <k>, normalized = FALSE
  #'
  #' @param args_list Liste d’arguments capturés depuis l’appel `cliques(...)`.
  #' @return list(k = <int>, normalized = <logical>)
  #' @keywords internal
  .normalize_cliques_args <- function(args_list) {
    k     <- 2L
    norm  <- FALSE

    if (length(args_list) == 0L) {
      return(list(k = k, normalized = norm))
    }

    nm <- names(args_list)

    # Forme abrégée : cliques(3)
    if (length(args_list) == 1L && (is.null(nm) || isTRUE(nm[1L] == ""))) {
      k <- as.integer(args_list[[1L]])
      return(list(k = k, normalized = norm))
    }

    # Nominal : clique_size=..., normalized=...
    if (!is.null(nm)) {
      if ("clique_size" %in% nm) k    <- as.integer(args_list[["clique_size"]])
      if ("normalized"  %in% nm) norm <- isTRUE(args_list[["normalized"]])
      return(list(k = k, normalized = norm))
    }

    # Sinon, tenter d'être robustes (ex: ordre des params non nommés)
    # 1er arg = k éventuel, 2e = normalized éventuel
    if (length(args_list) >= 1L) k    <- as.integer(args_list[[1L]])
    if (length(args_list) >= 2L) norm <- isTRUE(args_list[[2L]])
    list(k = k, normalized = norm)
  }

  #' Découper récursivement les termes d'un RHS
  #'
  #' Décompose une expression du côté droit d'une formule (RHS) contenant des
  #' additions de termes (ex. `edges + triangles + nodematch("grp")`) en une liste
  #' d'appels élémentaires, afin de pouvoir transformer chaque terme indépendamment.
  #'
  #' @param expr Une expression R correspondant au RHS d'une formule.
  #' @return Une liste d'appels (objets de type \code{call}) représentant les termes unitaires.
  #' @examples
  #' .split_sum_terms(quote(edges + triangles + nodematch("grp")))
  #' @keywords internal
  .split_sum_terms <- function(expr) {
    # Si c’est une somme, descendre récursivement à gauche et à droite
    if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
      return(c(
        .split_sum_terms(expr[[2]]),
        .split_sum_terms(expr[[3]])
      ))
    }
    # Sinon, paquetiser en liste à un élément
    list(expr)
  }

  build_bipartite_from_inputs <- function(partition = NULL,
                                        nodes     = NULL,
                                        dyads     = list()) {
    # Cas réseau déjà construit géré ailleurs
    stopifnot(!is.null(partition), is.atomic(partition), length(partition) >= 1L)

    # labels d’acteurs
    if (is.null(nodes)) {
      n      <- length(partition)
      labels <- sprintf("A%d", seq_len(n))
      nodes  <- data.frame(label = labels, stringsAsFactors = FALSE)
    } else {
      .check_nodes_df(nodes)
      if (nrow(nodes) != length(partition))
        stop("nrow(nodes) doit égaler length(partition).")
      labels <- as.character(nodes$label)
    }
    n <- length(labels)

    # valider dyads
    .check_dyads(dyads, n, labels)

    # groupes effectifs 1..G
    G       <- max(partition, na.rm = TRUE)
    g_names <- sprintf("G%d", seq_len(G))
    all_v   <- c(labels, g_names)

    # matrice adjacence bipartie
    adj <- matrix(0L, n + G, n + G,
                  dimnames = list(all_v, all_v))
    idx_actor  <- match(labels, all_v)
    idx_group  <- n + partition
    adj[cbind(idx_actor, idx_group)] <- 1L
    adj[cbind(idx_group, idx_actor)] <- 1L

    # réseau
    nw <- network::network(adj, directed = FALSE, matrix.type = "adjacency")
    network::set.network.attribute(nw, "bipartite", n)
    network::set.vertex.attribute(nw, "vertex.names", all_v)

    # attributs nodaux: copier colonnes de nodes (sauf label) sur les n premiers sommets
    if (ncol(nodes) > 1L) {
      for (a in setdiff(names(nodes), "label")) {
        vals <- nodes[[a]]
        network::set.vertex.attribute(nw, a, c(vals, rep(NA, G)))
      }
    }

    # covariables dyadiques: poser en attributs réseau (%n%) pour edgecov()
    if (length(dyads)) {
      for (nm in names(dyads)) {
        # edgecov opère sur le sous-graphe acteur-acteur via Proj1/B → on stocke n×n
        M <- dyads[[nm]]
        dimnames(M) <- list(labels, labels) # garantir ordre
        nw %n% nm <- M
      }
    }

    list(network = nw, partition = partition,
        actor_labels = labels, group_labels = g_names)
  }

  #' Traduire un terme ERPM vers un terme {ergm}
  #'
  #' Renomme un effet "ERPM-friendly" en son équivalent \pkg{ergm} via un
  #' dictionnaire, puis applique au besoin des encapsulations de forme
  #' \code{Proj1(~ ...)} et/ou \code{B(~ ..., form = "nonzero")}.
  #'
  #' @param term_call Un appel de fonction (objet \code{call}), p.ex. \code{cov_match("grp")}.
  #' @param rename_map Nom -> nouveau_nom, p.ex. \code{c(groups = "b2degrange", cov_match = "nodematch")}.
  #' @param wrap_proj1 Noms d’effets à encapsuler dans \code{Proj1(~ ...)}.
  #' @param wrap_B Noms d’effets à encapsuler dans \code{B(~ ..., form = "nonzero")}.
  #' @return Un objet \code{call} transformé, prêt à être inséré dans une formule \pkg{ergm}.
  #' @examples
  #' .translate_one_term(quote(cov_match("grp")),
  #'                    rename_map = c(cov_match = "nodematch"),
  #'                    wrap_proj1 = "nodematch")
  #' @keywords internal
  .translate_one_term <- function(term_call,
                                  rename_map,
                                  wrap_proj1 = character(),
                                  wrap_B     = character()) {

    # Accepter un symbole nu: groups  -> call(groups)
    if (is.symbol(term_call)) term_call <- as.call(list(term_call))
    # Passer au travers si non-call
    if (!is.call(term_call)) return(term_call)

    # Extraire symbole de fonction et arguments
    fun_sym   <- term_call[[1L]]
    args_list <- as.list(term_call)[-1L]

    # Cas : groups(...) -> b2degrange(from,to)
    if (is.symbol(fun_sym) && identical(fun_sym, as.name("groups"))) {
      gt <- .normalize_groups_args(args_list)  # Calculer from/to selon les cas
      # Construire l’appel b2degrange(from=?, to=?)
      return(as.call(list(as.name("b2degrange"), from = gt$from, to = gt$to)))
    }

    # Cas : cliques(...) -> cliques(clique_size=?, normalized=?)
    if (is.symbol(fun_sym) && identical(fun_sym, as.name("cliques"))) {
      # Valeurs par défaut ERPM/ERGM
      k    <- 2L
      norm <- FALSE

      # Lire args, en acceptant k, clique_size et normalized, + cliques(3) positionnel
      if (length(args_list) == 1L && (is.null(names(args_list)) || names(args_list)[1L] == "")) {
        k <- as.integer(args_list[[1L]])
      } else {
        al <- as.pairlist(args_list)
        if (!is.null(al$k))           k    <- as.integer(al$k)
        if (!is.null(al$clique_size)) k    <- as.integer(al$clique_size)
        if (!is.null(al$normalized))  norm <- isTRUE(al$normalized)
      }

      # Toujours renvoyer le terme cliques {ergm}, même pour k=1
      return(as.call(list(as.name("cliques"),
                          clique_size = k,
                          normalized  = norm)))
    }

    # Renommage générique via dictionnaire
    fname <- if (is.symbol(fun_sym)) as.character(fun_sym) else deparse(fun_sym)[1L]
    fname <- if (fname %in% names(rename_map)) rename_map[[fname]] else fname
    out <- as.call(c(as.name(fname), args_list))  # Recomposer l’appel

    # Encapsulations éventuelles
    if (fname %in% wrap_B)     out <- call("B",     call("~", out), form = "nonzero")
    if (fname %in% wrap_proj1) out <- call("Proj1", call("~", out))
    out
  }

  #' Wrapper ERPM → ERGM (traduction + fit optionnel)
  #'
  #' Traduit une formule “ERPM-friendly” (\code{groups}, \code{cov_match}, …) en
  #' termes \pkg{ergm} (\code{b2degrange}, \code{nodematch}, …), applique les
  #' encapsulations nécessaires (\code{Proj1(~ ...)}, \code{B(~ ..., form="nonzero")}),
  #' construit \code{ergm(~ nw + <RHS>, constraints=~b1part)}, puis évalue l'appel
  #' si demandé.
  #'
  #' @param formula Formule \code{nw ~ <termes ERPM>}.
  #' @param eval_call \code{logical}. \code{TRUE} (défaut) pour évaluer \code{ergm()}, \code{FALSE} pour renvoyer l'appel.
  #' @param verbose \code{logical}. Affiche la formule traduite et les contraintes.
  #' @param estimate \code{c("CD","MCMLE")}. Méthode d'estimation transmise à \code{ergm()}.
  #' @param eval.loglik \code{logical}. Évaluer la log-vraisemblance (bridge sampling).
  #' @param control Liste ou objet \code{control.ergm}. Si liste, convertie via \code{control.ergm()}.
  #' @param timeout \code{numeric}. Délai (s) pour interrompre l'évaluation, via \pkg{R.utils}.
  #'
  #' @return Un objet \code{ergm} si \code{eval_call=TRUE}, sinon un \code{call}.
  #' @examples
  #' \dontrun{
  #'   # Juste voir la traduction
  #'   erpm(nw ~ groups(from=2,to=3) + cov_match("grp"), eval_call=FALSE)
  #' }
  #' @export
  erpm <- function( formula, 
                  eval_call = TRUE, 
                  verbose = TRUE,
                  estimate = c("CD","MPLE","MLE","MCMLE"),
                  eval.loglik = TRUE,
                  control = NULL,
                  timeout = NULL,
                  # nouveaux canaux d'entrée optionnels (sans casser l'API)
                  nodes = NULL,
                  dyads = list()) {

    estimate <- match.arg(estimate)
    if (identical(estimate, "MCMLE")) estimate <- "MLE"
    if (!inherits(formula, "formula"))
      stop("The input should be a formula of the form `lhs ~ ...` with lhs = partition OU network biparti.")

    # --- Détection du LHS : partition vs network biparti ---
    env0    <- environment(formula) %||% parent.frame()
    lhs     <- formula[[2]]
    rhs_expr<- formula[[3]]

    # Évalue prudemment le LHS dans l’env de la formule, sinon parent.frame()
    lhs_val <- tryCatch(eval(lhs, envir = env0), error = function(e) e)

    # Cas A) LHS est une partition → construire le biparti + attributs
    if (!(inherits(lhs_val, "error")) &&
        is.atomic(lhs_val) && !inherits(lhs_val, "network")) {

      built <- build_bipartite_from_inputs(partition = lhs_val, nodes = nodes, dyads = dyads)
      nw2   <- built$network

      # Recréer une formule `nw ~ <rhs>` en liant `nw` dans l'environnement
      new_formula <- as.formula(bquote(nw ~ .(rhs_expr)))
      environment(new_formula) <- list2env(list(nw = nw2), parent = env0)

    } else {
      # Cas B) LHS déjà un network ou un symbole pointant dessus
      # Vérifie bipartite si possible
      if (!(inherits(lhs_val, "error")) && inherits(lhs_val, "network")) {
        bip <- tryCatch(lhs_val %n% "bipartite", error = function(e) NULL)
        if (is.null(bip) || is.na(bip)) stop("Le réseau LHS n'est pas biparti ou l'attribut %n% 'bipartite' est manquant.")
      }
      # Conserver la formule telle quelle, mais s’assurer que l’env contient l’objet
      new_formula <- formula
      if (!(inherits(lhs_val, "error"))) {
        environment(new_formula) <- env0
      }
    }

    # À partir d’ici, continuer exactement comme avant, en remplaçant `formula` par `new_formula`.
    orig_formula <- new_formula
    rhs_expr     <- new_formula[[3]]

    effect_rename_map <- c(
      cov_match     = "nodematch",
      cov_diff      = "absdiff",
      dyadcov       = "edgecov"
    )
    wrap_with_proj1 <- c("nodematch", "absdiff", "edgecov")
    wrap_with_B     <- c("edgecov")

    rhs_terms <- .split_sum_terms(rhs_expr)
    translated_terms <- lapply(
      rhs_terms,
      .translate_one_term,
      rename_map = effect_rename_map,
      wrap_proj1 = wrap_with_proj1,
      wrap_B     = wrap_with_B
    )

    new_rhs <- if (length(translated_terms) == 1L) translated_terms[[1L]] else
      Reduce(function(x, y) call("+", x, y), translated_terms)

    env <- environment(new_formula)
    new_formula[[3]] <- new_rhs
    if (!is.null(env)) environment(new_formula) <- env

    # 5) Contrainte par défaut bipartite
    constraint_expression <- as.formula(~ b1part)

    # 6) Contrôle ergm : accepter control.ergm existant, sinon construire
    ctrl <- if (inherits(control, "control.ergm")) control
            else if (is.null(control)) ergm::control.ergm()
            else do.call(ergm::control.ergm, as.list(control))

    # Placer le contrôle dans l’environnement d’éval et ne passer qu’un symbole dans l’appel
    ctrl_sym  <- as.name(sprintf(".ctrl_%08x", as.integer(runif(1, 0, .Machine$integer.max))))
    eval_env  <- if (!is.null(env)) env else parent.frame()
    assign(as.character(ctrl_sym), ctrl, envir = eval_env)

    # 7) Construire l’appel ergm(...)
    ergm_call <- as.call(list(
      as.name("ergm"),
      new_formula,
      constraints = constraint_expression,
      estimate    = estimate,
      eval.loglik = eval.loglik,
      control     = ctrl_sym
    ))

    # Journal compact du mapping complet
    if (isTRUE(verbose)) {
      init_str  <- .tight(.oneline(orig_formula))                # ex. nw~groups(3)
      final_fun <- as.call(list(as.name("ergm"), new_formula))   # ex. ergm(nw~b2degrange(...))
      final_str <- .tight(.oneline(final_fun))
      cat(sprintf("[ERPM] call initial : erpm(%s) -> call final : %s\n",
                  init_str, final_str))
      cat("\t Contraintes:  ~ b1part\n")
      cat("\t Options: estimate=", estimate,
          ", eval.loglik=", eval.loglik,
          ", control=", as.character(ctrl_sym), "\n", sep = "")
    }

    # 8) Évaluation ou renvoi du call
    if (!isTRUE(eval_call)) {
      if (isTRUE(verbose)) {
        cat("\t dry-run ergm call : ",
            paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "),
            "\n", sep = "")
      }
      return(ergm_call)
    }

    # Évaluation réelle
    if (is.null(timeout)) {
      return(eval(ergm_call, envir = eval_env))
    } else {
      return(R.utils::withTimeout(
        eval(ergm_call, envir = eval_env),
        timeout   = as.numeric(timeout),
        onTimeout = "silent"
      ))
    }
  }

  # ============================================================================
  # Export des fonctions
  # ============================================================================
  assign(".normalize_cliques_args", .normalize_cliques_args, envir = .GlobalEnv)
  assign(".normalize_groups_args",  .normalize_groups_args,  envir = .GlobalEnv)
  assign(".split_sum_terms",        .split_sum_terms,        envir = .GlobalEnv)
  assign(".translate_one_term",     .translate_one_term,     envir = .GlobalEnv)
  assign("erpm",                    erpm,                    envir = .GlobalEnv)
  assign(".__erpm_wrapper_loaded",  TRUE,                    envir = .GlobalEnv)
}