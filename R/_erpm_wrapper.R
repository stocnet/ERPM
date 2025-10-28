# # ==============================================================================
# # Fichier : erpm_wrapper.R
# # Fonctions : split_sum_terms(), translate_one_term(), erpm()
# # Utilité : Traduire des formules ERPM en appels {ergm} (renommage + encapsulations)
# # ==============================================================================

# if (!exists(".__erpm_wrapper_loaded", envir = .GlobalEnv)) {

#   #' Wrapper ERPM → ERGM
#   #'
#   #' Transforme une formule utilisant une grammaire “ERPM-friendly”
#   #' (p.ex. \code{groups}, \code{cov_match}, …) en une formule compatible \pkg{ergm}
#   #' (p.ex. \code{b2degrange}, \code{nodematch}, …), applique les encapsulations
#   #' nécessaires (p.ex. \code{Proj1(~ ...)}, \code{B(~ ..., form="nonzero")}), puis
#   #' construit l’appel \code{ergm()} avec la contrainte \code{~ b1part}.
#   #'
#   #' @param formula Formule de la forme \code{nw ~ <termes ERPM>}.
#   #' @param eval_call Logique. TRUE (défaut) pour évaluer \code{ergm()}, FALSE pour renvoyer l’appel.
#   #' @param verbose Logique. TRUE pour afficher la formule traduite et les contraintes.
#   #'
#   #' @return Un objet \code{ergm} si \code{eval_call = TRUE}, sinon un objet \code{call}.
#   #' @examples
#   #' \dontrun{
#   #'   erpm(nw ~ groups(from=2,to=3) + cov_match("grp"), eval_call = FALSE)
#   #'   fit <- erpm(nw ~ groups(from=2,to=3) + cov_match("grp"))
#   #'   summary(fit)
#   #' }
#   #' @export
#   erpm <- function(formula, eval_call = TRUE, verbose = TRUE) {
#     if (!inherits(formula, "formula")) stop("The input should be a formula of the form `nw ~ ...`.")

#     network_symbol <- formula[[2]]
#     rhs_expr       <- formula[[3]]

#     effect_rename_map <- c(
#       groups        = "b2degrange",
#       cov_match     = "nodematch",
#       cov_diff      = "absdiff",
#       dyadcov       = "edgecov",
#       squared_sizes = "squared_sizes"
#     )

#     wrap_with_proj1 <- c("nodematch", "absdiff", "edgecov")
#     wrap_with_B     <- c("edgecov")

#     rhs_terms <- split_sum_terms(rhs_expr)

#     translated_terms <- lapply(
#       rhs_terms,
#       translate_one_term,
#       rename_map = effect_rename_map,
#       wrap_proj1 = wrap_with_proj1,
#       wrap_B     = wrap_with_B
#     )

#     new_rhs <- if (length(translated_terms) == 1L) {
#       translated_terms[[1L]]
#     } else {
#       Reduce(function(x, y) call("+", x, y), translated_terms)
#     }

#     transformed_formula   <- call("~", network_symbol, new_rhs)
#     constraint_expression <- call("~", as.name("b1part"))

#     ergm_call <- as.call(list(as.name("ergm"), transformed_formula, constraints = constraint_expression))

#     if (isTRUE(verbose)) {
#       cat("\n[ERPM] Formule transformée:\n  "); print(transformed_formula)
#       cat("[ERPM] Contraintes:\n  ");          print(constraint_expression)
#     }

#     if (isTRUE(eval_call)) {
#       return(eval(ergm_call, envir = parent.frame()))
#     } else {
#       return(ergm_call)
#     }
#   }

#    #' Découper récursivement les termes d'un RHS
#   #'
#   #' Décompose une expression du côté droit d'une formule (RHS) contenant des
#   #' additions de termes (ex. `edges + triangles + nodematch("grp")`) en une liste
#   #' d'appels élémentaires, afin de pouvoir transformer chaque terme indépendamment.
#   #'
#   #' @param expr Une expression R correspondant au RHS d'une formule.
#   #' @return Une liste d'appels (objets de type \code{call}) représentant les termes unitaires.
#   #' @examples
#   #' split_sum_terms(quote(edges + triangles + nodematch("grp")))
#   split_sum_terms <- function(expr) {
#     if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
#       c(
#         split_sum_terms(expr[[2]]),
#         split_sum_terms(expr[[3]])
#       )
#     } else {
#       list(expr)
#     }
#   }

#   #' Traduire un terme ERPM vers un terme {ergm}
#   #'
#   #' Renomme un effet "ERPM-friendly" en son équivalent \pkg{ergm} via un
#   #' dictionnaire, puis applique au besoin des encapsulations de forme
#   #' \code{Proj1(~ ...)} et/ou \code{B(~ ..., form = "nonzero")}.
#   #'
#   #' @param term_call Un appel de fonction (objet \code{call}), p.ex. \code{cov_match("grp")}.
#   #' @param rename_map Nom -> nouveau_nom, p.ex. \code{c(groups = "b2degrange", cov_match = "nodematch")}.
#   #' @param wrap_proj1 Noms d’effets à encapsuler dans \code{Proj1(~ ...)}.
#   #' @param wrap_B Noms d’effets à encapsuler dans \code{B(~ ..., form = "nonzero")}.
#   #' @return Un objet \code{call} transformé, prêt à être inséré dans une formule \pkg{ergm}.
#   #' @examples
#   #' translate_one_term(quote(cov_match("grp")),
#   #'                    rename_map = c(cov_match = "nodematch"),
#   #'                    wrap_proj1 = "nodematch")
#   translate_one_term <- function(term_call,
#                                  rename_map,
#                                  wrap_proj1 = character(),
#                                  wrap_B     = character()) {
#     if (!is.call(term_call)){
#       return(term_call)
#     } 

#     original_fun <- as.character(term_call[[1]])
#     args_list    <- as.list(term_call)[-1]

#     new_fun_name <- if (original_fun %in% names(rename_map)) {
#       rename_map[[original_fun]]
#     } else {
#       original_fun
#     }

#     renamed_call <- as.call(c(as.name(new_fun_name), args_list))

#     if (new_fun_name %in% wrap_B) {
#       renamed_call <- call("B", call("~", renamed_call), form = "nonzero")
#     }
#     if (new_fun_name %in% wrap_proj1) {
#       renamed_call <- call("Proj1", call("~", renamed_call))
#     }

#     renamed_call
#   }

#   # --- Fonctions statiques ---
#   assign("split_sum_terms",        split_sum_terms,        envir = .GlobalEnv)
#   assign("translate_one_term",     translate_one_term,     envir = .GlobalEnv)
  
#   # --- Attributions à l’environnement global ---
#   assign("erpm",                   erpm,                   envir = .GlobalEnv)
#   assign(".__erpm_wrapper_loaded", TRUE,                   envir = .GlobalEnv)
# }


# # erpm <- function(formula) {
# #   if (!inherits(formula, "formula")) stop("The input should be a formula.")

# #   lhs <- formula[[2]] # nw
# #   rhs <- formula[[3]] # formula

# #   # --- 1. Renaming table ---
# #   rename_map <- c(
# #     groups        = "b2degrange",
# #     cov_match     = "nodematch",
# #     cov_diff      = "absdiff",
# #     dyadcov       = "edgecov",
# #     squared_sizes = "squared_sizes"
# #   )

# #   # --- 2. Proj1(~ ...) encapsulation list ---
# #   needs_proj1 <- c("nodematch", "absdiff", "edgecov")

# #   # --- 3. B(~ ..., form = "nonzero") encapsulation list ---
# #   needs_B <- c("edgecov")  # Exemple

# #   # --- Find individual terms ---
# #   extract_terms <- function(expr) {
# #     if (is.call(expr) && expr[[1]] == as.name("+")) {
# #       c(extract_terms(expr[[2]]), extract_terms(expr[[3]]))
# #     } else {
# #       list(expr)
# #     }
# #   }

# #   terms <- extract_terms(rhs)
# #   transformed_terms <- list()

# #   for (term in terms) {
# #     if (is.call(term)) {
# #       original_fun <- as.character(term[[1]])
# #       args <- as.list(term)[-1]

# #       # 1. Renaming
# #       renamed_fun <- if (original_fun %in% names(rename_map)) rename_map[[original_fun]] else original_fun
# #       new_term <- as.call(c(as.name(renamed_fun), args))

# #       # 2. B encapsulation
# #       if (renamed_fun %in% needs_B) {
# #         new_term <- call("B", call("~", new_term), form = "nonzero")
# #       }

# #       # 3. Proj1 encapsulation
# #       if (renamed_fun %in% needs_proj1) {
# #         new_term <- call("Proj1", call("~", new_term))
# #       }

# #       transformed_terms[[length(transformed_terms) + 1]] <- new_term
# #     } else {
# #       message(sprintf(" Unknown ERPM effect : %s", deparse(term)))
# #       transformed_terms[[length(transformed_terms) + 1]] <- term
# #     }
# #   }

# #   # Reconstruc the formula
# #   new_rhs <- if (length(transformed_terms) == 1) {
# #     transformed_terms[[1]]
# #   } else {
# #     Reduce(function(x, y) call("+", x, y), transformed_terms)
# #   }
# #   form <- list(
# #     formula = call("~", new_rhs),
# #     constraints = call("~", as.name("b1part")))
  
# #   #full_formula <- as.call(list(as.name("~"), as.name("nw"), form$formula[[2]]))
# #   full_formula <- as.call(list(as.name("~"), substitute(lhs), form$formula[[2]]))

# #   # build the final call
# #   return(call("ergm", full_formula, constraints = form$constraints))

# #   # Make the call if wanted
# #   # eval(call("ergm", full_formula, constraints = form$constraints)) 
# # }