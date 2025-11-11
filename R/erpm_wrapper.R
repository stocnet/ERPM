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

  #' Test minimal bipartite
  #' @keywords internal
  is_bipartite_network <- function(x) {
    inherits(x, "network") && !is.null(x %n% "bipartite") && is.finite(x %n% "bipartite")
  }

  .get_label_col <- function(nodes, prefer = "label") {
    stopifnot(is.data.frame(nodes))
    nms <- trimws(names(nodes))
    aliases <- c(prefer, "label", "nom", "name", "id")
    hit <- intersect(aliases, nms)
    if (length(hit)) hit[1L] else NULL
  }

  .check_nodes_df <- function(nodes) {
    stopifnot(is.data.frame(nodes))
    lab <- .get_label_col(nodes)           # peut être NULL
    if (!is.null(lab) && anyDuplicated(nodes[[lab]]))
      stop(sprintf("nodes$%s contient des doublons.", lab))
    invisible(TRUE)
  }
  # .check_nodes_df <- function(nodes) {
  #   stopifnot(is.data.frame(nodes), "label" %in% names(nodes))
  #   if (anyDuplicated(nodes$label)) stop("nodes$label contient des doublons.")
  #   invisible(TRUE)
  # }

  #' Validation des matrices dyadiques n×n
  #' @keywords internal
  .check_dyads <- function(dyads, n, labels) {
    if (length(dyads) == 0L) return(invisible(TRUE))
    stopifnot(is.list(dyads))
    for (nm in names(dyads)) {
      M <- dyads[[nm]]
      stopifnot(is.matrix(M), nrow(M) == n, ncol(M) == n)
      if (!is.null(rownames(M)) && !is.null(colnames(M))) {
        if (!identical(rownames(M), labels) || !identical(colnames(M), labels))
          stop(sprintf("dyads['%s']: row/colnames doivent matcher l’ordre des acteurs.", nm))
      }
    }
    invisible(TRUE)
  }

  #' Normaliser les arguments de groups(...) vers b2degrange(from, to)
  #'
  #' Convertit les variantes d'appel ERPM `groups`, `groups(taille)`,
  #' et `groups(from=.., to=..)` en une paire (from, to) prête pour `b2degrange`.
  #'
  #' Conventions :
  #' - `groups`            → [1, Inf)
  #' - `groups(taille)`    → [taille, taille+1)
  #' - `groups(from,to)`   → [from, to)
  #'
  #' @param args_list Liste d’arguments capturés depuis l’appel `groups(...)`.
  #' @return list(from=?, to=?), où `to` peut être `Inf` ou un appel `s+1`.
  #' @examples
  #' .normalize_groups_args(list())                 # -> from=1,   to=Inf
  #' .normalize_groups_args(list(2))                # -> from=2,   to=3
  #' .normalize_groups_args(list(from=2, to=5))     # -> from=2,   to=5
  #' @keywords internal
  .normalize_groups_args <- function(args_list) {
    # Cas 0 : groups -> [1, Inf)
    if (length(args_list) == 0L)
      return(list(from = 1, to = quote(Inf)))

    nm <- names(args_list)

    # Cas 1 : groups(taille) -> [taille, taille+1)
    if (length(args_list) == 1L && (is.null(nm) || isTRUE(nm[1L] == ""))) {
      s <- args_list[[1L]]
      # Si taille est un numérique scalaire fini, calculer s+1 maintenant
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
  #' Conventions :
  #' - `cliques()`                         → k=2, normalized=FALSE
  #' - `cliques(normalized=TRUE|FALSE)`    → k=2, normalized=<bool>
  #' - `cliques(clique_size=k)`            → k=<k>, normalized=FALSE
  #' - `cliques(k)`                        → k=<k>, normalized=FALSE
  #'
  #' @param args_list Liste d’arguments capturés depuis l’appel `cliques(...)`.
  #' @return list(k = <int>, normalized = <logical>)
  #' @keywords internal
  .normalize_cliques_args <- function(args_list) {
    k     <- 2L
    norm  <- FALSE

    if (length(args_list) == 0L) return(list(k = k, normalized = norm))

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

    # Robuste aux non nommés
    if (length(args_list) >= 1L) k    <- as.integer(args_list[[1L]])
    if (length(args_list) >= 2L) norm <- isTRUE(args_list[[2L]])
    list(k = k, normalized = norm)
  }

  #' Découper récursivement les termes d'un RHS
  #'
  #' Décompose `edges + triangles + nodematch("grp")` en liste de calls unitaires.
  #' @param expr Une expression R correspondant au RHS d'une formule.
  #' @return Liste de calls.
  #' @examples
  #' .split_sum_terms(quote(edges + triangles + nodematch("grp")))
  #' @keywords internal
  .split_sum_terms <- function(expr) {
    if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
      return(c(.split_sum_terms(expr[[2]]),
               .split_sum_terms(expr[[3]])))
    }
    list(expr)
  }

  #' Construire un biparti à partir d’une partition et d’inputs optionnels
  #' @keywords internal
  build_bipartite_from_inputs <- function(partition = NULL,
                                          nodes     = NULL,
                                          dyads     = list()) {
    stopifnot(!is.null(partition), is.atomic(partition), length(partition) >= 1L)

    # # 1) Labels acteurs
    # if (is.null(nodes)) {
    #   n      <- length(partition)
    #   labels <- sprintf("A%d", seq_len(n))
    #   nodes  <- data.frame(label = labels, stringsAsFactors = FALSE)
    # } else {
    #   .check_nodes_df(nodes)
    #   if (nrow(nodes) != length(partition))
    #     stop("nrow(nodes) doit égaler length(partition).")
    #   labels <- as.character(nodes$label)
    # }
    # n <- length(labels)
    # 1) Labels acteurs
    if (is.null(nodes)) {
      n      <- length(partition)
      labels <- sprintf("A%d", seq_len(n))
      nodes  <- data.frame(label = labels, stringsAsFactors = FALSE)
    } else {
      .check_nodes_df(nodes)
      if (nrow(nodes) != length(partition))
        stop("nrow(nodes) doit égaler length(partition).")
      labcol <- .get_label_col(nodes)  # peut être NULL

      if (is.null(labcol)) {
        labels <- sprintf("A%d", seq_len(nrow(nodes)))
        nodes$label <- labels
      } else {
        labels <- as.character(nodes[[labcol]])
        if (anyNA(labels) || any(!nzchar(labels))) stop("labels vides/NA non autorisés.")
        if (anyDuplicated(labels)) stop("labels en double non autorisés.")
        if (!("label" %in% names(nodes))) nodes$label <- labels
      }
    }
    n <- length(labels)


    # 2) Valider dyads
    .check_dyads(dyads, n, labels)

    # 3) Groupes 1..G
    G       <- max(partition, na.rm = TRUE)
    g_names <- sprintf("G%d", seq_len(G))
    all_v   <- c(labels, g_names)

    # 4) Matrice adjacency bipartie
    adj <- matrix(0L, n + G, n + G, dimnames = list(all_v, all_v))
    idx_actor  <- match(labels, all_v)
    idx_group  <- n + partition
    adj[cbind(idx_actor, idx_group)] <- 1L
    adj[cbind(idx_group, idx_actor)] <- 1L

    # 5) Réseau
    nw <- network::network(adj, directed = FALSE, matrix.type = "adjacency")
    network::set.network.attribute(nw, "bipartite", n)
    network::set.vertex.attribute(nw, "vertex.names", all_v)

    # 6) Attributs nodaux (évite doublons)
    if (ncol(nodes) > 1L) {
      for (a in setdiff(names(nodes), c("label"))) {   # exclut 'label' interne
        vals <- nodes[[a]]
        network::set.vertex.attribute(nw, a, c(vals, rep(NA, G)))
      }
    }

    # 7) dyads: inchangé, mais garantit l’ordre selon `labels`
    if (length(dyads)) {
      for (nm in names(dyads)) {
        M <- dyads[[nm]]
        dimnames(M) <- list(labels, labels)
        nw %n% nm <- M
      }
    }

    list(network = nw, partition = partition,
        actor_labels = labels, group_labels = g_names)
  }

  #   # 6) Attributs nodaux
  #   if (ncol(nodes) > 1L) {
  #     for (a in setdiff(names(nodes), "label")) {
  #       vals <- nodes[[a]]
  #       network::set.vertex.attribute(nw, a, c(vals, rep(NA, G)))
  #     }
  #   }

  #   # 7) Covariables dyadiques n×n (pour edgecov() après Proj1/B)
  #   if (length(dyads)) {
  #     for (nm in names(dyads)) {
  #       M <- dyads[[nm]]
  #       dimnames(M) <- list(labels, labels) # garantir l’ordre
  #       nw %n% nm <- M
  #     }
  #   }

  #   list(network = nw, partition = partition,
  #        actor_labels = labels, group_labels = g_names)
  # }

  #' Traduire un terme ERPM vers un terme {ergm}
  #'
  #' Renomme via dictionnaire, puis encapsule avec Proj1/B si requis.
  #' @param term_call call, p.ex. cov_match("grp")
  #' @param rename_map dict Nom->nouveau_nom
  #' @param wrap_proj1 noms à encapsuler par Proj1(~ ...)
  #' @param wrap_B noms à encapsuler par B(~ ..., form="nonzero")
  #' @return call prêt pour une formule {ergm}
  #' @keywords internal
  .translate_one_term <- function(term_call,
                                  rename_map,
                                  wrap_proj1 = character(),
                                  wrap_B     = character()) {

    # Symbole nu : groups -> call(groups)
    if (is.symbol(term_call)) term_call <- as.call(list(term_call))
    # Non-call : retourner tel quel
    if (!is.call(term_call)) return(term_call)

    # Extraire symbole et args
    fun_sym   <- term_call[[1L]]
    args_list <- as.list(term_call)[-1L]

    # Cas : groups(...) -> b2degrange(from,to)
    if (is.symbol(fun_sym) && identical(fun_sym, as.name("groups"))) {
      gt <- .normalize_groups_args(args_list)
      return(as.call(list(as.name("b2degrange"), from = gt$from, to = gt$to)))
    }

    # Cas : cliques(...) -> cliques(clique_size=?, normalized=?)
    if (is.symbol(fun_sym) && identical(fun_sym, as.name("cliques"))) {
      k <- 2L; norm <- FALSE
      if (length(args_list) == 1L && (is.null(names(args_list)) || names(args_list)[1L] == "")) {
        k <- as.integer(args_list[[1L]])
      } else {
        al <- as.pairlist(args_list)
        if (!is.null(al$k))           k    <- as.integer(al$k)
        if (!is.null(al$clique_size)) k    <- as.integer(al$clique_size)
        if (!is.null(al$normalized))  norm <- isTRUE(al$normalized)
      }
      return(as.call(list(as.name("cliques"),
                          clique_size = k,
                          normalized  = norm)))
    }

    # Renommage générique via dictionnaire
    fname <- if (is.symbol(fun_sym)) as.character(fun_sym) else deparse(fun_sym)[1L]
    fname <- if (fname %in% names(rename_map)) rename_map[[fname]] else fname
    out <- as.call(c(as.name(fname), args_list))  # Recomposer

    # Encapsulations éventuelles
    if (fname %in% wrap_B)     out <- call("B",     call("~", out), form = "nonzero")
    if (fname %in% wrap_proj1) out <- call("Proj1", call("~", out))
    out
  }

  #' Wrapper ERPM → ERGM (traduction + fit optionnel)
  #'
  #' Traduit une formule “ERPM-friendly” (groups, cov_match, …) en termes {ergm}
  #' (b2degrange, nodematch, …), ajoute Proj1/B si nécessaire, construit
  #' ergm(nw ~ <RHS>, constraints=~b1part), puis exécute l'appel si demandé.
  #'
  #' @param formula Formule `lhs ~ <termes ERPM>`. `lhs` = partition OU network biparti.
  #' @param eval_call logical. TRUE pour évaluer {ergm}, FALSE pour retourner le call.
  #' @param verbose logical. Affiche la traduction et les contraintes.
  #' @param estimate c("CD","MPLE","MLE","MCMLE"). "MCMLE" est normalisé en "MLE".
  #' @param eval.loglik logical. Log-vraisemblance via bridge sampling.
  #' @param control list ou control.ergm.
  #' @param timeout numeric. Délai max en secondes (via R.utils::withTimeout).
  #' @param nodes data.frame optionnel (colonnes nodales, inclut `label`).
  #' @param dyads list optionnelle de matrices n×n (edgecov sur Proj1).
  #' @return Objet {ergm} si eval_call=TRUE, sinon un call.
  #' @examples
  #' @export
  erpm <- function(formula, 
                   eval_call   = TRUE, 
                   verbose     = TRUE,
                   estimate    = c("CD","MPLE","MLE","MCMLE"),
                   eval.loglik = TRUE,
                   control     = NULL,
                   timeout     = NULL,
                   nodes       = NULL,
                   dyads       = list()) {

    estimate <- match.arg(estimate)
    if (identical(estimate, "MCMLE")) estimate <- "MLE"
    if (!inherits(formula, "formula"))
      stop("Attendu une formule `lhs ~ ...` avec lhs = partition OU network biparti.")

    # --- 0) Détecter la nature du LHS ----------------------------------------
    env0     <- environment(formula) %||% parent.frame()
    lhs      <- formula[[2]]
    rhs_expr <- formula[[3]]

    # Évaluer prudemment le LHS dans l'environnement de la formule
    lhs_val <- tryCatch(eval(lhs, envir = env0), error = function(e) e)

    # --- 1) Cas A : LHS est une partition -> construire biparti --------------
    if (!(inherits(lhs_val, "error")) &&
        is.atomic(lhs_val) && !inherits(lhs_val, "network")) {

      built <- build_bipartite_from_inputs(partition = lhs_val, nodes = nodes, dyads = dyads)
      nw2   <- built$network

      # Recréer une formule `nw ~ <rhs>` en liant `nw` dans l'environnement
      new_formula <- as.formula(bquote(nw ~ .(rhs_expr)))
      environment(new_formula) <- list2env(list(nw = nw2), parent = env0)

    } else {
      # --- 2) Cas B : LHS est déjà un network biparti -> NE PAS reconstruire ---
      if (!(inherits(lhs_val, "error")) && inherits(lhs_val, "network")) {
        # Vérification de l'attribut bipartite
        bip <- tryCatch(lhs_val %n% "bipartite", error = function(e) NULL)
        if (is.null(bip) || is.na(bip))
          stop("Le réseau LHS n'est pas biparti ou l'attribut %n% 'bipartite' est manquant.")
      }
      # Conserver la formule telle quelle et s’assurer que l’environnement contient l’objet
      new_formula <- formula
      if (!(inherits(lhs_val, "error"))) environment(new_formula) <- env0
    }

    # --- 3) Traduction des termes du RHS -------------------------------------
    orig_formula <- new_formula
    rhs_expr     <- new_formula[[3]]

    effect_rename_map <- c(
      cov_match = "nodematch",
      cov_diff  = "absdiff",
      dyadcov   = "edgecov"
      # Remarque : log_factorial_sizes et cliques restent inchangés
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

    new_rhs <- if (length(translated_terms) == 1L) translated_terms[[1L]]
               else Reduce(function(x, y) call("+", x, y), translated_terms)

    env <- environment(new_formula)
    new_formula[[3]] <- new_rhs
    if (!is.null(env)) environment(new_formula) <- env

    # --- 4) Contraintes -------------------------------------------------------
    constraint_expression <- as.formula(~ b1part)

    # --- 5) Contrôle ergm -----------------------------------------------------
    ctrl <- if (inherits(control, "control.ergm")) control
            else if (is.null(control)) ergm::control.ergm()
            else do.call(ergm::control.ergm, as.list(control))

    # Garde fou pour init pas de même longueur que le call
    k <- length(summary(new_formula, constraints = ~ b1part))
    if (!is.null(ctrl$init) && length(ctrl$init) != k) ctrl$init <- NULL

    # Placer le contrôle dans l’environnement d’éval et ne passer qu’un symbole
    ctrl_sym  <- as.name(sprintf(".ctrl_%08x", as.integer(runif(1, 0, .Machine$integer.max))))
    eval_env  <- if (!is.null(env)) env else parent.frame()
    assign(as.character(ctrl_sym), ctrl, envir = eval_env)

    

    # --- 6) Construire l’appel ergm(...) -------------------------------------
    ergm_call <- as.call(list(
      as.name("ergm"),
      new_formula,
      constraints = constraint_expression,
      estimate    = estimate,
      eval.loglik = eval.loglik,
      control     = ctrl_sym
    ))

    # --- 7) Journalisation compacte ------------------------------------------
    if (isTRUE(verbose)) {
      init_str  <- .tight(.oneline(orig_formula))
      final_fun <- as.call(list(as.name("ergm"), new_formula))
      final_str <- .tight(.oneline(final_fun))
      cat(sprintf("[ERPM] call initial : erpm(%s) -> call final : %s\n", init_str, final_str))
      cat("\t Contraintes:  ~ b1part\n")
      cat("\t Options: estimate=", estimate,
          ", eval.loglik=", eval.loglik,
          ", control=", as.character(ctrl_sym), "\n", sep = "")
    }

    # --- 8) Évaluer ou renvoyer le call --------------------------------------
    if (!isTRUE(eval_call)) {
      if (isTRUE(verbose)) {
        cat("\t dry-run ergm call : ",
            paste(deparse(ergm_call, width.cutoff = 500L), collapse = " "),
            "\n", sep = "")
      }
      return(ergm_call)
    }

    if (is.null(timeout)) {
      eval(ergm_call, envir = eval_env)
    } else {
      R.utils::withTimeout(
        eval(ergm_call, envir = eval_env),
        timeout   = as.numeric(timeout),
        onTimeout = "silent"
      )
    }
  }

  # ============================================================================
  # Export des fonctions
  # ============================================================================
  assign("build_bipartite_from_inputs", build_bipartite_from_inputs,  envir = .GlobalEnv)
  assign(".normalize_cliques_args",     .normalize_cliques_args,      envir = .GlobalEnv)
  assign(".normalize_groups_args",      .normalize_groups_args,       envir = .GlobalEnv)
  assign(".split_sum_terms",            .split_sum_terms,             envir = .GlobalEnv)
  assign(".translate_one_term",         .translate_one_term,          envir = .GlobalEnv)
  assign("erpm",                        erpm,                         envir = .GlobalEnv)
  assign(".__erpm_wrapper_loaded",      TRUE,                         envir = .GlobalEnv)
}