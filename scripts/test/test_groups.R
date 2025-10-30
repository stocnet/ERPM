# ======================================================================================
# Fichier : scripts/test/test_groups.R
# Objet   : Tests robustes pour l'effet ERPM `groups` → {ergm} `b2degrange`
# Exécution: Rscript scripts/test/test_groups.R
# ======================================================================================

# --------------------------------------------------------------------------------------
# Préambule: environnement et dépendances minimales
# --------------------------------------------------------------------------------------
# macOS/Linux: forcer UTF-8 pour affichages
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
#' get_script_dir()
#'
#' @keywords internal
get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 1L) return(normalizePath(dirname(f), winslash = "/", mustWork = FALSE))
  # fallback: si source() en interactif
  if (!is.null(sys.frames()) && !is.null(sys.calls())) {
    for (i in rev(seq_along(sys.calls()))) {
      cf <- sys.frame(i)
      if (!is.null(cf$ofile)) return(normalizePath(dirname(cf$ofile), winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
log_path   <- file.path(script_dir, "test_groups.log")

# Réinit
if (file.exists(log_path)) unlink(log_path, force = TRUE)

# pipe pour rediriger en parallèle stdout dans un fichier log
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

# Charger les fonctions biparti si non chargés
if (!exists("partition_to_bipartite_network", mode = "function")) {
  if (file.exists("R/functions_erpm_bip_network.R")) {
    source("R/functions_erpm_bip_network.R", local = FALSE)
  } else {
    stop("functions_erpm_bip_network.R introuvable. Place ce fichier sous R/.")
  }
}

# Charger le wrapper ERPM si non chargé
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
make_nw_from_partition <- function(part) {
  n <- length(part)                           # nb objets
  stopifnot(n >= 1, is.atomic(part))          # garde-fous
  lbl <- utils::head(LETTERS, n)              # labels "A","B",...
  partition_to_bipartite_network(labels = lbl, partition = part, attributes = list())
}

#' Compter le nombre de groupes dont la taille est comprise dans un intervalle donné
#'
#' Cette fonction calcule, à partir d'une partition, le nombre de groupes
#' dont la taille appartient à l'intervalle demi-ouvert \code{[from, to)}.
#'
#' @param part Un vecteur indiquant, pour chaque élément, le groupe auquel il appartient.
#' @param from Taille minimale de groupe à inclure (valeur par défaut : 1).
#' @param to Taille maximale de groupe à inclure (exclue, valeur par défaut : \code{Inf}).
#'
#' @return Un entier représentant le nombre de groupes dont la taille vérifie
#'         \code{from <= taille < to}.
#'
#' @examples
#' expected_groups_in_range(c(1, 2, 2, 3, 3, 3), from = 2, to = 4)
#' # Trois groupes au total : tailles 1, 2, 3 → seuls 2 et 3 sont dans [2,4), résultat = 2
#'
#' @keywords internal
expected_groups_in_range <- function(part, from, to) {
  sz <- as.integer(table(part))               # tailles par groupe
  sum(sz >= from & sz < to)                   # filtre demi-ouvert
}

#' Normaliser les arguments de l'effet \code{groups}
#'
#' Cette fonction standardise la signature d’appel d’un effet \code{groups()}
#' afin de produire une version cohérente et complète des arguments
#' \code{from} et \code{to}, ainsi qu’une représentation textuelle normalisée.
#'
#' @param args Une liste d’arguments (nommés ou non) passée à \code{groups()}.
#'        Peut contenir :
#'        \itemize{
#'          \item aucun argument — valeurs par défaut : \code{from = 1, to = Inf} ;
#'          \item un argument numérique seul — interprété comme une borne inférieure \code{from = k}, \code{to = k + 1} ;
#'          \item deux arguments nommés \code{from} et \code{to}.
#'        }
#'
#' @details
#' La fonction :
#' \enumerate{
#'   \item définit des valeurs par défaut si aucun argument n’est fourni ;
#'   \item interprète un seul argument non nommé comme une borne inférieure ;
#'   \item vérifie la validité numérique et la cohérence des bornes ;
#'   \item génère un champ \code{text} sous la forme \code{"b2degrange(from=...,to=...)"}.
#' }
#'
#' @return Une liste contenant :
#' \describe{
#'   \item{\code{from}}{borne inférieure du degré (numérique).}
#'   \item{\code{to}}{borne supérieure du degré (numérique ou \code{Inf}).}
#'   \item{\code{text}}{chaîne de caractères décrivant la signature normalisée.}
#' }
#'
#' @examples
#' normalize_groups_signature()
#' # $from = 1, $to = Inf, $text = "b2degrange(from=1,to=Inf)"
#'
#' normalize_groups_signature(list(2))
#' # $from = 2, $to = 3, $text = "b2degrange(from=2,to=3)"
#'
#' normalize_groups_signature(list(from = 1, to = 5))
#' # $from = 1, $to = 5, $text = "b2degrange(from=1,to=5)"
#'
#' @keywords internal
normalize_groups_signature <- function(args = list()) {
  if (length(args) == 0L) {
    return(list(from = 1, to = Inf, text = "b2degrange(from=1,to=Inf)"))
  }
  nm <- names(args)
  if (length(args) == 1L && (is.null(nm) || isTRUE(nm[1L] == ""))) {
    k <- as.numeric(args[[1L]])
    stopifnot(length(k) == 1L, is.finite(k))
    return(list(from = k, to = k + 1, text = sprintf("b2degrange(from=%g,to=%g)", k, k + 1)))
  }
  if (!is.null(nm) && all(c("from","to") %in% nm)) {
    f <- as.numeric(args[["from"]]); t <- as.numeric(args[["to"]])
    stopifnot(length(f) == 1L, length(t) == 1L, is.finite(f) || is.infinite(f), is.finite(t) || is.infinite(t))
    txt <- sprintf("b2degrange(from=%s,to=%s)",
                   if (is.infinite(f)) "Inf" else as.character(f),
                   if (is.infinite(t)) "Inf" else as.character(t))
    return(list(from = f, to = t, text = txt))
  }
  stop("Signature groups(...) inconnue pour ce test.")
}

#' Vérifier si une chaîne attendue est contenue dans un appel \code{ergm} traduit
#'
#' Cette fonction teste la présence d’une sous-chaîne donnée dans la version
#' compacte (sans espaces ni sauts de ligne) d’un appel \code{ergm()}.
#'
#' @param call_ergm Un objet de type \code{call} représentant l’appel \code{ergm} généré.
#' @param expected Une chaîne de caractères à rechercher dans la version compacte de l’appel.
#'
#' @details
#' L’appel est d’abord converti en texte avec \code{deparse()}, puis
#' tous les espaces et retours à la ligne sont supprimés avant recherche.
#'
#' @return Un booléen : \code{TRUE} si la sous-chaîne \code{expected} est trouvée, sinon \code{FALSE}.
#'
#' @examples
#' call <- quote(ergm(nw ~ b2degrange(from=1,to=5)))
#' check_translation_contains(call, "b2degrange(from=1,to=5)")
#' # Retourne TRUE
#'
#' @keywords internal
check_translation_contains <- function(call_ergm, expected) {
  line <- paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  compact <- gsub("\\s+", "", line)
  grepl(expected, compact, fixed = TRUE)
}

# ======================================================================================
# Noyau de test
# ======================================================================================

#' Exécuter un cas de test ERPM pour l'effet "groups"
#'
#' Construit un réseau biparti à partir d'une partition, calcule la statistique
#' de référence avec \code{summary( ~ b2degrange )}, évalue la valeur théorique via
#' \code{expected_groups_in_range()}, génère l'appel traduit avec \code{erpm()}
#' sans exécution, vérifie la traduction (\code{check_translation_contains()}),
#' affiche un récapitulatif, puis renvoie les résultats.
#'
#' @param part Vecteur de partition (entiers ou facteurs) définissant l'appartenance des nœuds.
#' @param name Nom court du cas de test (chaîne).
#' @param call_txt Terme ERPM à tester, sous forme de texte pour la formule (ex. \code{"groups(from=1,to=5)"}).
#' @param args Liste d'arguments pour la signature \code{groups} (au moins \code{from}, \code{to}).
#'
#' @details
#' Étapes :
#' \enumerate{
#'   \item \code{make_nw_from_partition(part)} pour obtenir le réseau biparti.
#'   \item \code{normalize_groups_signature(args)} pour normaliser \code{from}/\code{to}.
#'   \item Statistique ERGM: \code{summary(nw ~ b2degrange(from=...,to=...))}.
#'   \item valeur théorique: \code{expected_groups_in_range(part, from, to)}.
#'   \item Traduction ERPM: \code{erpm(f, eval_call = FALSE, verbose = TRUE)} avec \code{f <- as.formula(paste0("nw ~ ", call_txt))}.
#'   \item Vérification de traduction: \code{check_translation_contains(call_ergm, gsub("\\s+", "", ng$text))}.
#'   \item Impression d'un récapitulatif sur la console.
#' }
#'
#' @return Une liste avec :
#' \describe{
#'   \item{\code{ok_stat}}{Booléen, TRUE si \code{summary} = valeur théorique.}
#'   \item{\code{ok_trad}}{Booléen, TRUE si la traduction ERPM→ERGM contient la signature attendue.}
#'   \item{\code{stat}}{Valeur entière renvoyée par \code{summary}.}
#'   \item{\code{expected}}{Valeur entière attendue selon la partition.}
#'   \item{\code{call_line}}{Chaîne compactée de l'appel traduit.}
#' }
#'
#' @examples
#' run_one_case(c(1,2,2,3,3,3), "demo", "groups(from=1,to=3)", list(from=1, to=3))
#'
#' @keywords internal
run_one_case <- function(part, name, call_txt, args) {
  nw <- make_nw_from_partition(part)                         # réseau biparti
  ng <- normalize_groups_signature(args)                     # from/to attendus
  stat_val <- as.integer(summary(as.formula(nw ~ b2degrange(from = ng$from, to = ng$to))))
  truth    <- expected_groups_in_range(part, ng$from, ng$to)

  f <- as.formula(paste0("nw ~ ", call_txt))                 # formule ERPM
  environment(f) <- list2env(list(nw = nw), parent = parent.frame())
  call_ergm <- erpm(f, eval_call = FALSE, verbose = TRUE)    # montre aussi les calls

  ok_trad <- check_translation_contains(call_ergm, gsub("\\s+", "", ng$text))

  cat(sprintf("\n[CAS %-18s] part={%s}", name, paste(part, collapse=",")))
  cat(sprintf("\t  attendu(groups)  : %s", ng$text))
  cat(sprintf("\t  summary(.)       : %d", stat_val))
  cat(sprintf("\t  attendu(part)    : %d", truth))
  cat(sprintf("\t  traduction erpm  : %s\n", if (ok_trad) "OK" else "KO"))

  list(
    ok_stat = identical(stat_val, truth),
    ok_trad = ok_trad,
    stat = stat_val,
    expected = truth,
    call_line = paste(deparse(call_ergm, width.cutoff = 500L), collapse = " ")
  )
}

#' Exécuter une série de cas de test ERPM pour une partition donnée
#'
#' Cette fonction applique \code{run_one_case()} à l’ensemble des cas décrits
#' dans un ensemble de cas de test  (\code{panel}), pour une partition donnée, et
#' agrège les résultats sous forme de tableau.
#'
#' @param part Un vecteur de partition (entiers ou facteurs) définissant
#'        l’appartenance des nœuds à des groupes.
#' @param panel Une liste de cas de test, où chaque élément doit contenir :
#'        \itemize{
#'          \item \code{name} : nom du cas de test ;
#'          \item \code{call_txt} : texte de l’appel ERPM à tester ;
#'          \item \code{args} : liste d’arguments pour la signature.
#'        }
#'
#' @details
#' Pour chaque cas du ensemble de cas de test, la fonction :
#' \enumerate{
#'   \item exécute \code{run_one_case(part, name, call_txt, args)} ;
#'   \item extrait les résultats principaux (\code{ok_stat}, \code{ok_trad}, \code{stat}, \code{expected}) ;
#'   \item compile l’ensemble sous forme de \code{data.frame}.
#' }
#'
#' @return Un \code{data.frame} contenant une ligne par cas testé, avec les colonnes :
#' \code{case}, \code{ok_stat}, \code{ok_trad}, \code{stat}, \code{expected}.
#'
#' @examples
#' panel <- list(
#'   list(name = "cas1", call_txt = "groups(from=1,to=3)", args = list(from=1, to=3)),
#'   list(name = "cas2", call_txt = "groups(from=2,to=4)", args = list(from=2, to=4))
#' )
#' run_cases_for_partition(c(1,2,2,3,3,3), panel)
#'
#' @keywords internal
run_cases_for_partition <- function(part, panel) {
  res <- lapply(panel, function(cx) {
    out <- run_one_case(part, cx$name, cx$call_txt, cx$args)
    data.frame(
      case     = cx$name,
      ok_stat  = out$ok_stat,
      ok_trad  = out$ok_trad,
      stat     = out$stat,
      expected = out$expected,
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
  list(name = "groups_all",         call_txt = "groups",                args = list()),
  list(name = "groups_exact_1",     call_txt = "groups(1)",             args = list(1)),
  list(name = "groups_exact_2",     call_txt = "groups(2)",             args = list(2)),
  list(name = "groups_exact_3",     call_txt = "groups(3)",             args = list(3)),
  list(name = "groups_range_2_4",   call_txt = "groups(from=2,to=4)",   args = list(from = 2, to = 4)),
  list(name = "groups_range_1_Inf", call_txt = "groups(from=1,to=Inf)", args = list(from = 1, to = Inf))
)

# ======================================================================================
# Run principal
# ======================================================================================

#' Exécuter l'ensemble de tests \code{ERPM: groups -> b2degrange}
#'
#' Lance l'ensemble des tests pour l'effet \code{groups} en vérifiant :
#' (i) l'égalité des statistiques calculées par \code{summary(~ b2degrange)} et la valeur théorique,
#' (ii) la traduction \code{erpm} -> \code{ergm}.  
#' Gère le logging vers fichier et console, agrège les résultats, et échoue si des validations sont en erreur.
#'
#' @details
#' - Journalisation dans \code{scripts/test/test_groups.log} avec redirection de \code{stdout} et \code{message}.  
#' - \code{set.seed(42)} pour la reproductibilité.  
#' - Pour chaque partition de \code{partitions} et chaque cas de \code{cases} :
#'   exécute \code{run_cases_for_partition()}, cumule les validations, et imprime un bilan global.  
#' - Interrompt avec \code{stop()} si au moins une validation échoue.  
#' - Suppose \code{partitions} et \code{cases} disponibles dans l'environnement appelant.
#'
#' @return
#' Une liste nommée (invisible) de \code{data.frame}, un par partition,
#' contenant les colonnes \code{case}, \code{ok_stat}, \code{ok_trad}, \code{stat}, \code{expected}.
#'
#' @examples
#' # Exécuter tous les tests et consulter le log :
#' res <- run_all_tests()
#' # file.show("scripts/test/test_groups.log")
#'
#' @seealso \code{\link{run_cases_for_partition}}, \code{\link{run_one_case}},
#'          \code{\link{normalize_groups_signature}}, \code{\link{check_translation_contains}}
#' @keywords internal
run_all_tests <- function() {
  # --- Logging console + fichier (dans une fonction => on.exit marche) ---
  log_path <- file.path("scripts","test","test_groups.log")
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(log_path)) unlink(log_path)

  sink(file = log_path, split = TRUE)          # stdout -> fichier + console
  con_msg <- file(log_path, open = "at")       # messages -> même fichier
  sink(con_msg, type = "message")

  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(close(con_msg),        silent = TRUE)
    try(sink(),                silent = TRUE)
    flush.console()
  }, add = TRUE)

  # --- Tests ---
  set.seed(42)
  cat("=== TEST ERPM: groups → b2degrange ===\n")

  all_results <- list(); total_ok <- 0L; total_n <- 0L
  for (nm in names(partitions)) {
    cat(sprintf("\n--- Partition %s ---\n", nm))
    df <- run_cases_for_partition(partitions[[nm]], cases)
    all_results[[nm]] <- df
    total_ok <- total_ok + sum(df$ok_stat) + sum(df$ok_trad)
    total_n  <- total_n  + 2L * nrow(df)
  }

  cat(sprintf("\n=== Bilan global : %d / %d validations OK ===\n", total_ok, total_n))
  if (total_ok < total_n) stop(sprintf("Echec: %d validations KO", total_n - total_ok))
  invisible(all_results)
}

if (identical(environment(), globalenv())) run_all_tests()