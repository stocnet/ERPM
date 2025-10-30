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

    # Cas spécial : groups(...) -> b2degrange(from,to)
    if (is.symbol(fun_sym) && identical(fun_sym, as.name("groups"))) {
      gt <- .normalize_groups_args(args_list)  # Calculer from/to selon les cas
      # Construire l’appel b2degrange(from=?, to=?)
      return(as.call(list(as.name("b2degrange"), from = gt$from, to = gt$to)))
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
                    estimate = c("CD","MPLE", "MLE", "MCMLE"),
                    eval.loglik = TRUE,
                    control = NULL,
                    timeout = NULL) {

    # Valider arguments
    estimate <- match.arg(estimate)
    if (identical(estimate, "MCMLE")) estimate <- "MLE"

    if (!inherits(formula, "formula"))
      stop("The input should be a formula of the form `nw ~ ...`.")

    # Conserver la formule d’origine pour affichage
    orig_formula <- formula
    rhs_expr     <- formula[[3]]

    # Dictionnaires de traduction/encapsulation
    effect_rename_map <- c(
      cov_match     = "nodematch",
      cov_diff      = "absdiff",
      dyadcov       = "edgecov"
    )
    wrap_with_proj1 <- c("nodematch", "absdiff", "edgecov")
    wrap_with_B     <- c("edgecov")

    # 1) Découper le RHS en termes unitaires
    rhs_terms <- .split_sum_terms(rhs_expr)

    # 2) Traduire chaque terme
    translated_terms <- lapply(
      rhs_terms,
      .translate_one_term,
      rename_map = effect_rename_map,
      wrap_proj1 = wrap_with_proj1,
      wrap_B     = wrap_with_B
    )

    # 3) Recomposer le RHS avec des '+'
    new_rhs <- if (length(translated_terms) == 1L) translated_terms[[1L]] else
      Reduce(function(x, y) call("+", x, y), translated_terms)

    # 4) Recréer la formule avec le RHS traduit en préservant l’environnement
    env <- environment(formula)
    new_formula <- formula
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
  assign(".normalize_groups_args",  .normalize_groups_args,  envir = .GlobalEnv)
  assign(".split_sum_terms",        .split_sum_terms,        envir = .GlobalEnv)
  assign(".translate_one_term",     .translate_one_term,     envir = .GlobalEnv)
  assign("erpm",                    erpm,                    envir = .GlobalEnv)
  assign(".__erpm_wrapper_loaded",  TRUE,                    envir = .GlobalEnv)
}