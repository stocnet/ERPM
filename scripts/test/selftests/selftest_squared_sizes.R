# ======================================================================================
# Fichier : scripts/test/selftests/selftest_squared_sizes.R
# Objet   : Self-test autonome pour l'effet ERPM `squared_sizes`
# Exécution: Rscript scripts/test/selftests/selftest_squared_sizes.R
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule: environnement et dépendances minimales
# --------------------------------------------------------------------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
    if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' requis.")
    if (!requireNamespace("ergm",    quietly = TRUE)) stop("Package 'ergm' requis.")
})

suppressMessages(suppressPackageStartupMessages({
    library(network, quietly = TRUE, warn.conflicts = FALSE)
    library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))

#' Obtenir le répertoire du script en cours d'exécution
#'
#' Cette fonction renvoie le chemin absolu du répertoire contenant
#' le script R actuellement exécuté, que celui-ci soit lancé avec
#' `Rscript`, `source()`, ou en session interactive.
#'
#' @details
#' - Si le script est lancé avec `Rscript`, le chemin est extrait de `--file=`.
#' - Si le script est exécuté via `source()`, la fonction parcourt la pile
#'   d'appels pour retrouver l'origine (`ofile`).
#' - En dernier recours, elle renvoie le répertoire de travail courant (`getwd()`).
#'
#' @return Un chemin absolu (chaîne de caractères) normalisé avec des barres obliques `/`.
#' @examples
#' .get_script_dir()
#'
#' @keywords internal
.get_script_dir <- function() {
    a <- commandArgs(FALSE)
    f <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(f) == 1L) return(normalizePath(dirname(f), winslash = "/", mustWork = FALSE))
    if (!is.null(sys.frames()) && !is.null(sys.calls())) {
        for (i in rev(seq_along(sys.calls()))) {
        cf <- sys.frame(i)
        if (!is.null(cf$ofile)) return(normalizePath(dirname(cf$ofile), winslash = "/", mustWork = FALSE))
        }
    }
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- .get_script_dir()
log_path   <- file.path(script_dir, "selftest_squared_sizes.log")

# réinit
if (file.exists(log_path)) unlink(log_path, force = TRUE)

# connexions dédiées
con_out <- file(log_path, open = "wt")  # sortie standard
con_err <- file(log_path, open = "at")  # messages/erreurs en append

sink(con_out, split = TRUE)             # stdout → fichier + console
sink(con_err, type = "message")         # messages → fichier

on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(close(con_err),        silent = TRUE)
    try(sink(),                silent = TRUE)
    try(close(con_out),        silent = TRUE)
    flush.console()
}, add = TRUE)

# Charger le terme ERGM squared_sizes si non chargé
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    # si tu es dans le package, charge tout (inclut C/registration si déjà compilé)
    devtools::load_all(quiet = TRUE)
} else {
    stop("InitErgmTerm.squared_sizes introuvable. Place le fichier sous R/ ou exécute depuis le package (load_all).")
}

# Charger les fonctions biparti si non chargées
if (!exists("partition_to_bipartite_network", mode = "function")) {
    if (file.exists("R/functions_erpm_bip_network.R")) {
        source("R/functions_erpm_bip_network.R", local = FALSE)
    } else {
        stop("functions_erpm_bip_network.R introuvable. Place ce fichier sous R/.")
    }
}

# Charger le wrapper ERPM si tu veux aussi inspecter la traduction
if (!exists("erpm", mode = "function")) {
  if (file.exists("R/erpm_wrapper.R")) {
    source("R/erpm_wrapper.R", local = FALSE)
  } else {
    stop("erpm_wrapper.R introuvable. Place ce fichier sous R/.")
  }
}

# ======================================================================================
# Fonctions locales
# ======================================================================================

#' Construire un réseau biparti depuis une partition
#' @param part vecteur de groupes par objet, ex. c(1,2,2,3,3,3).
#' @return objet `network` bipartite prêt pour {ergm}.
#' @keywords internal
.make_nw_from_partition <- function(part) {
    n <- length(part)
    stopifnot(n >= 1, is.atomic(part))
    lbl <- utils::head(LETTERS, n)
    partition_to_bipartite_network(labels = lbl, partition = part, attributes = list())
}

#' Calculer les tailles de groupes à partir d'une partition
#'
#' Cette fonction calcule le nombre d'éléments dans chaque groupe
#' à partir d'un vecteur de partition.
#'
#' @param part Un vecteur indiquant, pour chaque élément, le groupe auquel il appartient.
#'
#' @return Un vecteur d'entiers représentant la taille de chaque groupe,
#'         ordonné selon les identifiants de groupes présents dans `part`.
#'
#' @examples
#' .group_sizes_from_partition(c(1, 2, 2, 3, 3, 3))
#' # Renvoie : 1 2 3
#'
#' @keywords internal
.group_sizes_from_partition <- function(part) {
    as.integer(table(part))
}

#' Calculer la somme des tailles de groupes élevées à une puissance donnée
#'
#' Cette fonction calcule la somme des tailles de groupes (issues d'une partition)
#' élevées à la puissance `pow`, en ne considérant que les groupes dont la taille
#' est comprise entre `from` (inclus) et `to` (exclus).
#'
#' @param part Un vecteur de partition indiquant l'appartenance de chaque élément à un groupe.
#' @param from Taille minimale des groupes à inclure (valeur par défaut : 1).
#' @param to Taille maximale des groupes à inclure (valeur par défaut : \code{Inf}).
#' @param pow Puissance à laquelle élever les tailles de groupes (valeur par défaut : 2).
#'
#' @return Un scalaire numérique correspondant à la somme des tailles de groupes sélectionnés
#'         élevées à la puissance spécifiée.
#'
#' @examples
#' .expected_squared_sizes(c(1, 2, 2, 3, 3, 3))
#' # Calcule 1^2 + 2^2 + 3^2 = 14
#'
#' .expected_squared_sizes(c(1, 1, 2, 2, 2), from = 2, pow = 3)
#' # Calcule uniquement les groupes de taille >= 2 : 3^3 + 2^3 = 35
#'
#' @keywords internal
.expected_squared_sizes <- function(part, from = 1, to = Inf, pow = 2) {
    sz <- .group_sizes_from_partition(part)
    idx <- (sz >= from) & (sz < to)
    sum( (sz[idx])^pow )
}

#' Calculer la statistique d'identité au carré sur un réseau biparti
#'
#' Cette fonction évalue la statistique utilisée pour vérifier l'effet
#' \texttt{squared\_sizes} sur un réseau biparti.  
#' Elle combine le nombre total d’arêtes et la somme des combinaisons
#' possibles entre les degrés du mode 2.
#'
#' @param nw Un objet \code{network} biparti (issu du package \pkg{network}).
#'
#' @details
#' La fonction :
#' \enumerate{
#'   \item identifie la frontière entre les deux modes du graphe ;
#'   \item calcule les degrés des nœuds du second mode ;
#'   \item renvoie la somme du nombre d’arêtes et de deux fois la somme des combinaisons \code{choose(deg2, 2)}.
#' }
#'
#' @return Un scalaire numérique correspondant à la valeur totale de la statistique.
#'
#' @examples
#' nw <- network.initialize(6, bipartite = 3)
#' add.edges(nw, tail = c(1,2,3), head = c(4,5,6))
#' .identity_pow2_all(nw)
#'
#' @keywords internal
.identity_pow2_all <- function(nw) {
    # Garde-fous explicites
    stopifnot(inherits(nw, "network"))
    n_attr <- network::get.network.attribute(nw, "bipartite")
    if (is.null(n_attr) || is.na(n_attr)) stop("Réseau non biparti ou attribut 'bipartite' manquant.")
    n1 <- as.integer(n_attr)
    n  <- network::network.size(nw)
    if (!(n1 >= 0L && n1 < n)) stop("Attribut 'bipartite' incohérent avec la taille du réseau.")

    # Mode 2
    v2 <- seq.int(n1 + 1L, n)

    # Degrés du mode 2 sans conversion en matrice
    deg2 <- vapply(
    v2,
    function(v) length(network::get.neighborhood(nw, v, type = "combined")),
    integer(1L)
    )

    # Nombre d'arêtes sans summary()
    ecount <- network::network.edgecount(nw)

    # Identité: sum_g d_g^2 = edges + 2 * sum_g C(d_g, 2)
    ecount + 2L * sum(choose(deg2, 2))
}

#' Normaliser les arguments de l'effet \code{squared_sizes}
#'
#' Cette fonction prépare et complète la liste d’arguments associée à l’effet
#' \code{squared_sizes}, en appliquant des valeurs par défaut et en produisant
#' une représentation textuelle standardisée.
#'
#' @param args Une liste d’arguments optionnels pouvant contenir
#'        \code{from}, \code{to} et \code{pow}.
#'
#' @details
#' - Les valeurs par défaut sont : \code{from = 1}, \code{to = Inf}, \code{pow = 2}.  
#' - Les arguments fournis remplacent sélectivement ces valeurs.  
#' - Un champ \code{text} est ajouté pour fournir une version textuelle
#'   complète de la configuration (ex. : \code{"squared_sizes(from=1,to=Inf,pow=2)"}).
#'
#' @return Une liste contenant les éléments \code{from}, \code{to}, \code{pow}
#'         et \code{text}.
#'
#' @examples
#' .normalize_squared_sizes_signature(list(from = 2, pow = 3))
#' # Renvoie une liste équivalente à :
#' # $from = 2, $to = Inf, $pow = 3,
#' # $text = "squared_sizes(from=2,to=Inf,pow=3)"
#'
#' @keywords internal
.normalize_squared_sizes_signature <- function(args = list()) {
    # défaut : from=1, to=Inf, pow=2
    out <- list(from = 1, to = Inf, pow = 2, text = "squared_sizes(from=1,to=Inf,pow=2)")
    if (length(args)) {
    nm <- names(args)
    if (!is.null(nm) && length(nm)) {
        if ("from" %in% nm) out$from <- as.numeric(args[["from"]])
        if ("to"   %in% nm) out$to   <- as.numeric(args[["to"]])
        if ("pow"  %in% nm) out$pow  <- as.numeric(args[["pow"]])
    }
    }
    txt <- sprintf("squared_sizes(from=%s,to=%s,pow=%s)",
                    if (is.infinite(out$from)) "Inf" else as.character(out$from),
                    if (is.infinite(out$to))   "Inf" else as.character(out$to),
                    as.character(out$pow))
    out$text <- txt
    out
}

#' Vérifier la présence d'un motif dans un appel \code{ergm} traduit
#'
#' Cette fonction teste si une chaîne attendue (\code{expected})
#' est bien présente dans la version compacte (sans espaces)
#' d’un appel \code{ergm} déparsé.
#'
#' @param call_ergm Un appel \R{} (objet de type \code{call}) correspondant
#'        à la traduction générée pour \code{ergm()}.
#' @param expected Une chaîne de caractères à rechercher dans l’appel.
#'
#' @details
#' La fonction :
#' \enumerate{
#'   \item convertit l’appel en texte avec \code{deparse()};
#'   \item supprime tous les espaces et sauts de ligne ;
#'   \item vérifie si la chaîne \code{expected} est contenue dans le texte.
#' }
#'
#' @return Un booléen : \code{TRUE} si la sous-chaîne est trouvée, sinon \code{FALSE}.
#'
#' @examples
#' call <- quote(ergm(nw ~ b2degrange(from=1,to=5)))
#' .check_translation_contains(call, "b2degrange(from=1,to=5)")
#' # Renvoie TRUE
#'
#' @keywords internal
.check_translation_ok <- function(call_ergm, term = "squared_sizes", args = list()) {
    line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
    compact <- gsub("\\s+", "", line)

    # défauts reconnus par InitErgmTerm.squared_sizes
    def <- list(from = 1, to = Inf, pow = 2)

    # 1) le terme doit être présent
    if (!grepl(paste0("\\b", term, "\\("), compact)) return(FALSE)

    # 2) on n’exige la présence QUE des arguments ≠ défauts
    checks <- logical(0)
    if (!is.null(args$from) && !identical(as.numeric(args$from), def$from)) {
    checks <- c(checks, grepl(paste0("from=", gsub("\\s+", "", as.character(args$from))), compact, fixed = TRUE))
    }
    if (!is.null(args$to) && !(is.infinite(args$to) && is.infinite(def$to)) &&
        !identical(as.numeric(args$to), def$to)) {
    checks <- c(checks, grepl(paste0("to=", gsub("\\s+", "", as.character(args$to))), compact, fixed = TRUE))
    }
    if (!is.null(args$pow) && !identical(as.numeric(args$pow), def$pow)) {
    checks <- c(checks, grepl(paste0("pow=", gsub("\\s+", "", as.character(args$pow))), compact, fixed = TRUE))
    }

    # si aucun arg non-défaut → OK dès que le terme est présent
    if (!length(checks)) return(TRUE)
    all(checks)
}
# .check_translation_contains <- function(call_ergm, expected) {
#   line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
#   compact <- gsub("\\s+", "", line)
#   grepl(expected, compact, fixed = TRUE)
# }

# ======================================================================================
# Noyau de test
# ======================================================================================

#' Exécuter un cas de test pour l'effet \code{squared_sizes}
#'
#' Construit un réseau biparti à partir d'une partition, normalise la signature
#' \code{squared_sizes}, calcule la valeur thérorique, évalue la statistique
#' via \code{summary()} d'ERGM, vérifie éventuellement la traduction \code{erpm}
#' et teste l'identité de contrôle pour \code{pow = 2} et \code{[from,to) = [1,Inf)}.
#' Un récapitulatif est imprimé sur la console.
#'
#' @param part Vecteur de partition définissant l'appartenance des nœuds.
#' @param name Nom court du cas de test.
#' @param call_txt Terme ERPM à évaluer dans la formule, ex. \code{"squared_sizes(from=1,to=Inf,pow=2)"}.
#' @param args Liste d'arguments pour \code{squared_sizes} (peut contenir \code{from}, \code{to}, \code{pow}).
#'
#' @details
#' Étapes principales :
#' \enumerate{
#'   \item \code{nw <- .make_nw_from_partition(part)}
#'   \item \code{sg <- .normalize_squared_sizes_signature(args)}
#'   \item valeur thérorique : \code{.expected_squared_sizes(part, sg$from, sg$to, sg$pow)}
#'   \item Stat ERGM : \code{summary(as.formula(paste0("nw ~ ", call_txt)))}
#'   \item Traduction (optionnelle) : si \code{erpm()} existe, \code{.check_translation_contains()}
#'   \item Identité de contrôle : si \code{pow=2} et intervalle complet, comparaison avec \code{.identity_pow2_all(nw)}
#' }
#'
#' @return Une liste avec :
#' \describe{
#'   \item{\code{ok_stat}}{Booléen. \code{summary()} égale à la valeur thérorique.}
#'   \item{\code{ok_trad}}{Booléen ou \code{NA}. Traduction \code{erpm} conforme à la signature attendue.}
#'   \item{\code{ok_ident}}{Booléen ou \code{NA}. Identité de contrôle valide pour \code{pow=2}.}
#'   \item{\code{stat}}{Valeur numérique de \code{summary()}.}
#'   \item{\code{expected}}{Valeur numérique attendue.}
#' }
#'
#' @examples
#' .run_one_case(c(1,2,2,3,3,3), "pow2_full", "squared_sizes(from=1,to=Inf,pow=2)",
#'              list(from=1, to=Inf, pow=2))
#'
#' @seealso \code{\link{.normalize_squared_sizes_signature}},
#'          \code{\link{.expected_squared_sizes}},
#'          \code{\link{.identity_pow2_all}},
#'          \code{\link{.check_translation_contains}}
#' @keywords internal
.run_one_case <- function(part, name, call_txt, args) {
    nw <- .make_nw_from_partition(part)
    sg <- .normalize_squared_sizes_signature(args)

    # Valeur attendue depuis la partition
    truth <- .expected_squared_sizes(part, from = sg$from, to = sg$to, pow = sg$pow)

    # Valeur via {ergm}
    f <- as.formula(paste0("nw ~ ", call_txt))
    environment(f) <- list2env(list(nw = nw), parent = parent.frame())
    stat_val <- as.numeric(summary(f))

    # (Optionnel) inspection traduction via erpm si dispo
    ok_trad <- NA
    if (exists("erpm", mode = "function")) {
        call_ergm <- erpm(f, eval_call = FALSE, verbose = TRUE)
        ok_trad <- .check_translation_ok(call_ergm, term = "squared_sizes", args = args)
    }

    # Identité de contrôle pour pow=2 et [1,Inf)
    ok_identity <- NA
    if (isTRUE(sg$pow == 2 && isTRUE(sg$from == 1) && isTRUE(is.infinite(sg$to)))) {
        id_val <- .identity_pow2_all(nw)
        ok_identity <- identical(unname(as.integer(stat_val)), unname(as.integer(id_val)))
    }

    cat(sprintf("\n[CAS %-18s] part={%s}", name, paste(part, collapse=",")))
    cat(sprintf("\t  appel          : %s", sg$text))
    cat(sprintf("\t  summary(.)     : %s", format(stat_val)))
    cat(sprintf("\t  attendu(part)  : %s", format(truth)))
    if (!is.na(ok_identity)) cat(sprintf("\t  identité pow2   : %s", if (ok_identity) "OK" else "KO"))
    if (!is.na(ok_trad))     cat(sprintf("\t  traduction erpm : %s", if (ok_trad) "OK" else "KO"))
    cat("\n")

    list(
        ok_stat    = identical(unname(as.integer(stat_val)), unname(as.integer(truth))),
        ok_trad    = ok_trad,
        ok_ident   = ok_identity,
        stat       = stat_val,
        expected   = truth
    )
}

#' Exécuter une série de cas de test pour une partition donnée (\code{squared_sizes})
#'
#' Cette fonction applique \code{.run_one_case()} à chaque cas de test défini
#' dans un panneau (\code{panel}) pour une partition donnée, puis regroupe
#' les résultats dans un tableau récapitulatif.
#'
#' @param part Un vecteur de partition (entiers ou facteurs) indiquant
#'        l’appartenance des nœuds à des groupes.
#' @param panel Une liste de cas de test, où chaque élément doit contenir :
#'        \itemize{
#'          \item \code{name} : nom du cas de test (chaîne) ;
#'          \item \code{call_txt} : texte de la formule ERPM à tester ;
#'          \item \code{args} : liste d’arguments de l’effet (\code{from}, \code{to}, \code{pow}, etc.).
#'        }
#'
#' @details
#' Pour chaque cas du panneau :
#' \enumerate{
#'   \item Exécute \code{.run_one_case(part, name, call_txt, args)} ;
#'   \item Extrait les principaux indicateurs de validation :
#'         \code{ok_stat}, \code{ok_trad}, \code{ok_ident}, \code{stat}, \code{expected} ;
#'   \item Convertit le tout en \code{data.frame} pour un affichage et un traitement faciles.
#' }
#'
#' @return Un \code{data.frame} combiné contenant une ligne par cas de test et les colonnes :
#' \code{case}, \code{ok_stat}, \code{ok_trad}, \code{ok_ident}, \code{stat}, \code{expected}.
#'
#' @examples
#' panel <- list(
#'   list(name = "full_pow2", call_txt = "squared_sizes(from=1,to=Inf,pow=2)", args = list(from=1, to=Inf, pow=2)),
#'   list(name = "range_pow3", call_txt = "squared_sizes(from=2,to=5,pow=3)", args = list(from=2, to=5, pow=3))
#' )
#' .run_cases_for_partition(c(1,2,2,3,3,3), panel)
#'
#' @seealso \code{\link{.run_one_case}}
#' @keywords internal
.run_cases_for_partition <- function(part, panel) {
    res <- lapply(panel, function(cx) {
            out <- .run_one_case(part, cx$name, cx$call_txt, cx$args)
            data.frame(
                case      = cx$name,
                ok_stat   = out$ok_stat,
                ok_trad   = if (is.na(out$ok_trad)) NA else out$ok_trad,
                ok_ident  = if (is.na(out$ok_ident)) NA else out$ok_ident,
                stat      = out$stat,
                expected  = out$expected,
                stringsAsFactors = FALSE
            )
        })
  do.call(rbind, res)
}

# ======================================================================================
# Jeu de tests
# ======================================================================================
partitions <- list(
    P1 = c(1, 2, 2, 3, 3, 3),
    P2 = c(1, 1, 2, 3, 3, 4, 4, 4),
    P3 = c(1, 1, 1, 2, 2, 3),
    P4 = c(1, 2, 3, 4, 5),
    P5 = rep(1, 6)
)

cases <- list(
    list(name="sq_all_pow2",    call_txt="squared_sizes",                   args=list(from=1, to=Inf, pow=2)),
    list(name="sq_2to5_pow2",   call_txt="squared_sizes(from=2,to=5)",      args=list(from=2, to=5,   pow=2)),
    list(name="sq_all_pow3",    call_txt="squared_sizes(pow=3)",            args=list(from=1, to=Inf, pow=3)),
    list(name="sq_1to2_pow2",   call_txt="squared_sizes(from=1,to=2)",      args=list(from=1, to=2,   pow=2)),
    list(name="sq_3toInf_pow2", call_txt="squared_sizes(from=3,to=Inf)",    args=list(from=3, to=Inf, pow=2))
)

# ======================================================================================
# Run principal
# ======================================================================================

#' Exécuter l'ensemble de tests \code{ERPM: squared_sizes}
#'
#' Lance tous les tests pour l'effet \code{squared_sizes} :
#' (i) comparaison statistique \code{summary()} vs valeur thérorique,
#' (ii) identité de contrôle (si applicable),
#' (iii) vérification de la traduction \code{erpm} -> \code{ergm} (si disponible).
#' Gère le logging, agrège les résultats, et échoue si des validations sont en erreur.
#'
#' @details
#' - Journalisation dans \code{scripts/test/selftests/selftest_squared_sizes.log} avec redirection de
#'   \code{stdout} et \code{message}.  
#' - \code{set.seed(42)} pour la reproductibilité.  
#' - Pour chaque partition définie dans \code{partitions} et chaque cas de \code{cases} :
#'   exécute \code{.run_cases_for_partition()}, imprime le \code{data.frame} de résultats,
#'   et cumule les validations :
#'   \code{ok_stat} (toujours), \code{ok_ident} (si testable), \code{ok_trad} (si \code{erpm} présent).  
#' - Interrompt avec \code{stop()} si au moins une validation échoue.
#' - Suppose \code{partitions} et \code{cases} disponibles dans l'environnement.
#'
#' @return
#' Une liste nommée (invisible) de \code{data.frame}, un par partition, contenant
#' les colonnes \code{case}, \code{ok_stat}, \code{ok_trad}, \code{ok_ident}, \code{stat}, \code{expected}.
#'
#' @examples
#' res <- run_all_tests()
#' # file.show("scripts/test/selftests/selftest_squared_sizes.log")
#'
#' @seealso \code{\link{.run_cases_for_partition}}, \code{\link{.run_one_case}},
#'          \code{\link{.normalize_squared_sizes_signature}},
#'          \code{\link{.expected_squared_sizes}},
#'          \code{\link{.check_translation_contains}}
#' @keywords internal
run_all_tests <- function() {
    # Logging dans scripts/test/selftests/selftest_squared_sizes.log
    log_path <- file.path("scripts","test", "selftests", "selftest_squared_sizes.log")
    dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(log_path)) unlink(log_path)

    sink(file = log_path, split = TRUE)
    con_msg <- file(log_path, open = "at")
    sink(con_msg, type = "message")

    on.exit({
        try(sink(type = "message"), silent = TRUE)
        try(close(con_msg),        silent = TRUE)
        try(sink(),                silent = TRUE)
        flush.console()
    }, add = TRUE)

    set.seed(42)
    cat("=== TEST ERPM: squared_sizes ===\n")

    all_results <- list(); total_ok <- 0L; total_n <- 0L
    for (nm in names(partitions)) {
        cat(sprintf("\n--- Partition %s ---\n", nm))
        df <- .run_cases_for_partition(partitions[[nm]], cases)
        all_results[[nm]] <- df
        # ok_stat : 1 validation / cas
        total_ok <- total_ok + sum(df$ok_stat, na.rm = TRUE)
        total_n  <- total_n  + sum(!is.na(df$ok_stat))
        # identité (si applicable)
        total_ok <- total_ok + sum(df$ok_ident, na.rm = TRUE)
        total_n  <- total_n  + sum(!is.na(df$ok_ident))
        # traduction (si erpm dispo)
        total_ok <- total_ok + sum(df$ok_trad, na.rm = TRUE)
        total_n  <- total_n  + sum(!is.na(df$ok_trad))
        print(df)
    }

    cat(sprintf("\n=== Bilan global : %d / %d validations OK ===\n", total_ok, total_n))
    if (total_ok < total_n) stop(sprintf("Echec: %d validations KO", total_n - total_ok))
    invisible(all_results)
}

if (identical(environment(), globalenv())) run_all_tests()