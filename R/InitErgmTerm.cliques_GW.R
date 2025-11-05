# ==============================================================================
# Fichier : InitErgmTerm.cliques_GW.R
# Terme   : cliques_GW(lambda = 2)
# Stat    : Somme_g λ[1 - ((λ-1)/λ)^{deg_g}], deg_g = taille du groupe
# Vectorisé sur 'lambda'
# ==============================================================================

#' ERGM term: cliques_GW
#'
#' @param lambda numeric, >= 1, scalaire ou vecteur. Défaut 2.
#' @export
InitErgmTerm.cliques_GW <- function(nw, arglist, ..., version = packageVersion("ergm")) {
    termname <- "cliques_GW"

    a <- check.ErgmTerm(
    nw, arglist,
    directed      = NULL,
    bipartite     = TRUE,
    varnames      = c("lambda"),
    vartypes      = c("numeric"),
    defaultvalues = list(2),
    required      = c(FALSE)
    )

    lambda <- a$lambda
    if (length(lambda) == 0L) return(NULL)

    # -- Validation domaine
    if (any(!is.finite(lambda)))
        ergm_Init_stop(sQuote(termname), ": 'lambda' non fini (NA/NaN/Inf) interdit.")
    if (any(lambda < 1))
        ergm_Init_stop(sQuote(termname), ": 'lambda' doit être >= 1.")

    # -- Pré-calcule r = (lambda-1)/lambda côté R
    r <- (lambda - 1) / lambda

    # -- Utilitaire pour produire des noms de coefficients comme : cliques_GW_lambda2, cliques_GW_lambda1.5, cliques_GW_lambda0.125
    pretty_num <- function(x){
        s <- formatC(x, digits = 4, format = "fg", flag = "#")  # Format général, 6 digits, sans . en tant que dernier caractère
        sub("\\.$", "", s)                                      # Supprime le point final si présent ("2." → "2").
    }
    coef.names <- paste0("cliques_GW_lambda", vapply(lambda, pretty_num, ""))

    # -- Entrées aplaties par colonnes: (lambda_j, r_j)
    inputs <- c(rbind(as.double(lambda), as.double(r)))

    list(
        name         = "cliques_GW",                 # doit matcher c_cliques_GW (sans le c_)
        coef.names   = coef.names,
        inputs       = inputs,                        # [2*j+0]=lambda_j, [2*j+1]=r_j
        dependence   = TRUE,                          # dépend des toggles
        emptynwstats = numeric(length(lambda))        # réseau vide -> 0
    )
}