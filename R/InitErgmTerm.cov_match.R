# ==============================================================================
# Fichier : InitErgmTerm.cov_match.R
# Terme   : cov_match(cov, clique_size = 2, category = NULL, normalized = c("none","by_group","global"))
# Stat    : non normalisée :
#             S_k(B;c)          = sum_g sum_r C(n_{g,r}, k)
#           ciblée (category=κ) :
#             S_k^{(κ)}(B;c)    = sum_g C(n_{g,κ}, k)
#           by_group :
#             sum_g [ S_k(g)/C(n_g,k) ]       (ou C(n_{g,κ},k)/C(n_g,k) en ciblé)
#           global :
#             (1/C(N1,k)) * S_k(·)
# ==============================================================================

#' ERGM term: cov_match
#'
#' @param cov character|factor|numeric  Nom d'un attribut acteur (factor/character) ou vecteur (longueur |A|).
#' @param clique_size integer|numeric   Taille(s) k des cliques monochromatiques (k >= 1). Vectorisable.
#' @param category character|NULL       Modalité ciblée (si attribut catégoriel).
#' @param normalized character|logical  "none"|"by_group"|"global" ; TRUE≡"by_group", FALSE≡"none".
#' @export
InitErgmTerm.cov_match <- function(nw, arglist, ..., version = packageVersion("ergm")) {
    termname <- "cov_match"   
    options(erpm.debug.cov_match_init = FALSE)

    a <- check.ErgmTerm(
        nw, arglist,
        directed      = NULL,
        bipartite     = TRUE,
        varnames      = c("cov","clique_size","category","normalized"),
        vartypes      = c("numeric,character","numeric","character","logical,character"),
        defaultvalues = list(NULL,            2,            NULL,       "none"),
        required      = c(TRUE,               FALSE,        FALSE,      FALSE)
    )

    # --- 0) Biparti ------------------------------------------------------
    n1 <- tryCatch(nw %n% "bipartite", error = function(e) NA_integer_)
    if (!is.numeric(n1) || !is.finite(n1) || n1 <= 0)
        ergm_Init_stop(sQuote(termname), ": réseau non biparti ou attribut %n% 'bipartite' manquant/invalide.")

    # --- 1) Normaliser  ------------------------------------------------
    normalized <- a$normalized
    if (is.logical(normalized)) {
        normalized <- if (isTRUE(normalized)) "by_group" else "none"
    }
    normalized <- match.arg(tolower(as.character(normalized)), c("none","by_group","global"))
    norm_mode  <- switch(normalized, none=0L, by_group=1L, global=2L)

    # --- 2) Normaliser les tailles k (clique_size) -------------------------------
    ks <- as.integer(a$clique_size)
    if (any(!is.finite(ks)) || any(ks < 1L)) ergm_Init_abort(sQuote(termname), ": 'clique_size' doit contenir des entiers >= 1.")

    .allow_k1_nn <- isTRUE(getOption("ERPM.allow.k1.nonnormalized", FALSE))

    # Interdit k=1 si non normalisé ou global, sauf si override explicite
    if (any(ks == 1L) && (normalized %in% c("none","global")) && !.allow_k1_nn) {
    ergm_Init_abort(
        sQuote(termname),
        ": cov_match(..., clique_size=1) avec normalized='", normalized,
        "' est constant pour les fits ERGM. ",
        "Utilisez k>=2, normalized='by_group', ou offset(...). ",
        "Pour forcer: options(ERPM.allow.k1.nonnormalized=TRUE)."
    )
    }

    ks <- unique(ks)
    K  <- length(ks)

    # --- 3) Construire modalités par acteur -------------------------
    cov      <- a$cov
    category <- a$category

    get_actor_codes <- function(nw, cov, category = NULL) {
        n1 <- as.integer(nw %n% "bipartite")
        ia <- seq_len(n1)  # indices des acteurs
        vn <- network::network.vertex.names(nw)
        if (length(vn) >= n1) {
        ia_guess <- which(!grepl("^G\\d+$", vn))
        if (length(ia_guess) == n1) ia <- ia_guess
        }

        if (is.character(cov) && length(cov) == 1L) {
        vals <- network::get.vertex.attribute(nw, cov)
        if (is.null(vals))
            ergm_Init_stop(sQuote(termname), ": attribut inexistant: ", sQuote(cov), ".")
        x <- vals[ia]

        # Forcer en facteur pour coder proprement les modalités
        f <- as.factor(x)
        z <- as.integer(f)          # codes 1..R, NA -> NA
        z[!is.finite(z)] <- 0L      # 0 = "absent/indéfini"

        kappa_code <- 0L
        cov_label  <- cov
        if (!is.null(category)) {
            # Si catégorie absente, l'ajouter aux levels pour obtenir une fréquence nulle sans erreur
            if (!(category %in% levels(f))) {
                levels(f) <- c(levels(f), category)
                # Recalculer les codes avec les nouveaux levels
                z <- as.integer(f)
                z[!is.finite(z)] <- 0L
            }
            kappa_code <- as.integer(match(category, levels(f)))
            cov_label  <- paste0(cov, "==", category)
        }

        # Retourne les codes, la catégorie éventuelle, et les niveaux pour debug
        return(list(
            z          = as.double(z),
            kappa_code = as.double(kappa_code),
            cov_label  = cov_label,
            levels     = levels(f)
        ))

        } else {
        # Vecteur numérique direct
        # Ici, 'category' ne s'applique pas
        x_num <- suppressWarnings(as.numeric(cov))
        if (any(!is.finite(x_num))) ergm_Init_stop(sQuote(termname), ": vecteur 'cov' contient NA/NaN/Inf.")
        if (length(x_num) < n1)     ergm_Init_stop(sQuote(termname), ": longueur(cov) < |A| = ", n1, ".")
        # Pour un attribut numérique, 'cov_match' n'a pas de sens (besoin de classes).
        ergm_Init_stop(sQuote(termname), ": 'cov_match' requiert un attribut catégoriel (factor/character).")
        }
    }

    ax <- get_actor_codes(nw, cov = cov, category = category)
    z_codes    <- ax$z          # double*, codes 0 (NA) ou 1..R
    kappa_code <- ax$kappa_code  # 0 si non ciblé
    cov_label  <- ax$cov_label
    levs       <- ax$levels %||% character(0)

    has_kappa <- as.double(as.integer(kappa_code > 0))

    # ===================== DEBUG LOCAL (AJOUT) =====================
    # Activer avec : options(erpm.debug.cov_match_init = TRUE)
    DEBUG_COV_MATCH_INIT <- isTRUE(getOption("erpm.debug.cov_match_init", FALSE)) || FALSE
    if (DEBUG_COV_MATCH_INIT) {
        cat(sprintf("[Init:%s] n1=%d | normalized=%s (mode=%d) | K=%d | ks={%s}\n",
                    termname, n1, normalized, norm_mode, K, paste(ks, collapse=",")))
        cat(sprintf("[Init:%s] cov=%s | has_kappa=%s | kappa_code=%s\n",
                    termname, deparse(substitute(cov)), as.logical(has_kappa), as.integer(kappa_code)))
        if (length(levs)) {
        map <- paste(sprintf("%d:%s", seq_along(levs), levs), collapse=", ")
        cat(sprintf("[Init:%s] levels map: %s\n", termname, map))
        zi <- as.integer(z_codes)
        zi[!is.finite(zi)] <- 0L
        tab <- as.integer(table(factor(zi, levels=0:length(levs))))
        cat(sprintf("[Init:%s] freq codes (0=NA): {%s}\n",
                    termname, paste(tab, collapse=", ")))
        }
    }
    # ===============================================================

    # --- 4) Inputs pour le C -----------------------------------------------------
    # Layout INPUT_PARAM :
    #   [0] = n1
    #   [1] = K
    #   [2] = norm_mode  (0 none, 1 by_group, 2 global)
    #   [3] = has_kappa  (0/1)
    #   [4] = kappa_code (0 si non ciblé)
    #   [5..(5+K-1)] = ks (k>=1)
    #   [5+K .. 5+K+(n1-1)] = z[1..n1] (codes modalités par acteur)
    inputs <- c(
        as.double(n1),
        as.double(K),
        as.double(norm_mode),
        has_kappa,
        as.double(kappa_code),
        as.double(ks),
        z_codes
    )

    # --- 5) Noms des coefficients -----------------------------------------------
    # Ex.: cov_match[sex]_k2, cov_match[sex==F]_k3_bygrp, etc.
    suffix_norm <- switch(normalized, none="", by_group="_bygrp", global="_glob")
    coef.names  <- paste0("cov_match[", cov_label, "]_k", ks, suffix_norm)

    # --- 6) Spécification du terme ----------------------------------------------
    list(
        name         = "cov_match",   # doit matcher C_CHANGESTAT_FN(c_cov_match)
        coef.names   = coef.names,    # longueur = K
        inputs       = inputs,        # tel que ci-dessus
        dependence   = TRUE,
        emptynwstats = 0
    )
}
