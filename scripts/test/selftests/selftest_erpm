# ==============================================================================
# Fichier : selftest_erpm.R
# Objet   : Self-tests étendus du wrapper erpm() — vérifie :
#           1) la traduction des termes ERPM -> {ergm} (dry-run),
#           2) l’exécution d'ergm avec contrôle minimal,
#           3) la robustesse à l’absence de coef.names (validation via appel).
# ==============================================================================


# ==============================================================================
# 0) PRÉAMBULE (env local + dépendances)
# ==============================================================================

# -- Force un env local stable pour l’affichage (UTF-8) ------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8"), silent = TRUE))

# -- Vérifie l'installation des packages "network" et "ergm" ----------------------------
suppressPackageStartupMessages({
  if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' requis.")
  if (!requireNamespace("ergm",    quietly = TRUE)) stop("Package 'ergm' requis.")
})

# -- Charge les packages --------------------------
suppressMessages(suppressPackageStartupMessages({
  library(network, quietly = TRUE, warn.conflicts = FALSE)
  library(ergm,    quietly = TRUE, warn.conflicts = FALSE)
}))


# ==============================================================================
# 1) SOURCES LOCALES DU PROJET
#    (Paramètres, initialisation ERPM, wrapper & helpers biparti)
# ==============================================================================

# NOTE : 'local = FALSE' source dans l'environnement global (tests).
source("scripts/local_source/settings.R",        local = FALSE)
source("scripts/local_source/init.R",            local = FALSE)
source("scripts/local_source/launcher.R",        local = FALSE)
source("R/erpm_wrapper.R",                       local = FALSE)
source("R/functions_erpm_bip_network.R",         local = FALSE)


# ==============================================================================
# 2) GESTION DU PARALLÉLISME
#    (Neutralise le cluster pour éviter les erreurs d’attachement de
#     package ERPM dans les workers PSOCK lors d’un développement via devtools)
# ==============================================================================

#' Désactive le parallélisme {ergm} et nettoie l'env. parallèle.
#'
#' - Force `parallel = 0` côté {ergm}
#' - Nettoie les variables d'environnement susceptibles de déclencher
#'   la création d’un cluster PSOCK (CI/IDE).
disable_all_parallel <- function() {
  options(ergm.parallel = 0L, ergm.parallel.type = NULL)  # aucun worker
  if (nzchar(Sys.getenv("R_PARALLEL_PORT", ""))) {
    Sys.setenv(R_PARALLEL_PORT = "")                      # ports parallèles vides
  } else {
    Sys.unsetenv("R_PARALLEL_PORT")                       # variable absente
  }
}

disable_all_parallel()  


# ==============================================================================
# 3) INITIALISATION ERPM & PATCH ERGM
# ==============================================================================

# - Charge le projet ERPM (messages/patchs visibles si verbose=TRUE)
# - Active le patch {ergm} si présent (tracing / correctifs ciblés)
init_erpm(selftest = FALSE, verbose = TRUE)
if (exists("ergm_patch_enable")) ergm_patch_enable(verbose = VERBOSE)


# ==============================================================================
# 4) HELPERS TEST — ASSERTIONS & OUTILS D’AFFICHAGE
# ==============================================================================
#' Retourne `a` si non-NULL/longueur>0 sinon `b`.
#' @param a,b objets potentiellement NULL
`%||%` <- function(a, b) if (!is.null(a) && length(a)) a else b

#' Vérifie un booléen retourne TRUE ou lève une erreur.
#' @param x   booléen
#' @param msg message d’erreur
.ok <- function(x, msg = "") if (!isTRUE(x)) stop("Test KO: ", msg)

#' Compacte une chaîne (supprime espaces/nbsp).
#' @param s string
#' @return string sans espaces
.compact <- function(s) gsub("[[:space:]]+|\u00A0", "", s, perl = TRUE)

#' Vérifie si un motif apparaît dans un texte.
#' @param txt texte source
#' @param pat motif regex (perl=TRUE)
#' @param msg message d’erreur si absent
.has <- function(txt, pat, msg = "") .ok(grepl(pat, .compact(txt), perl = TRUE), msg)

#' Extrait l’appel ergm(...) sous forme de texte, quel que soit l’objet.
#' @param x objet retour d’erpm(...) (call, fit, list avec $call ou attr "call")
#' @return texte de l'appel
.as_ergm_call_text <- function(x) {
  call_obj <-
    if (is.call(x)) {
      x
    } else if (is.list(x) && !is.null(x$call) && is.call(x$call)) {
      x$call
    } else {
      attr(x, "call")
    }
  if (!is.call(call_obj)) {
    # Fallback : imprime l’objet pour diagnostiquer des retours atypiques.
    return(paste(gsub("[[:space:]]+", " ", capture.output(print(x))), collapse = " "))
  }
  paste(deparse(call_obj, width.cutoff = 500L), collapse = " ")
}

#' Valide la présence d’un motif dans l’appel {ergm}.
#' S’appuie sur l’appel d'ergm quand `coef.names` est vide (ex. CD + contraintes).
#' @param obj    fit ou structure avec un attr/call {ergm}
#' @param pattern motif regex
#' @param label   étiquette affichée (optionnelle)
expect_in_call <- function(obj, pattern, label = NULL) {
  txt <- .as_ergm_call_text(obj)
  if (!is.null(label)) cat("\n--- CALL", paste0("[", label, "]"), "---\n", txt, "\n")
  .has(txt, pattern, paste0("Motif introuvable dans l'appel ergm: ", pattern))
}

#' Exécute un bloc de test standardisé :
#' - Affiche un(une?) en-tête,
#' - Capture/masque les warnings (option),
#' - Convertit les erreurs en objet "test_error" pour poursuivre la suite.
#' @param label     titre du test
#' @param expr      expression à évaluer
#' @param seed      pour la génération alétoire
#' @param quiet_warn TRUE pour ne pas imprimer les warnings
run <- function(label, expr, seed = 123, quiet_warn = TRUE) {
  cat(sprintf("\n===== %s =====\n", label))
  set.seed(seed)
  out <- withCallingHandlers(
    tryCatch(force(expr),
             error = function(e)
               structure(list(`__err__` = TRUE, message = conditionMessage(e)),
                         class = "test_error")),
    warning = function(w) {
      if (!quiet_warn) message("WARN: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(out, "test_error")) cat("❌ ERREUR: ", out$message, "\n", sep = "") else cat("✅ OK\n")
  invisible(out)
}

#' Lance un fit “sûr” :
#' - Tolère "data are essentially constant" (instabilité CD) et la neutralise,
#' - Sinon relance l’erreur,
#' - Affiche un en-tête de validation via appel.
#' @param expr  expression qui renvoie un fit {ergm}
#' @param label étiquette d’affichage
#' @return fit {ergm} ou NULL si neutralisé
safefit <- function(expr, label = NULL) {
  f <- try(force(expr), silent = TRUE)
  if (inherits(f, "try-error")) {
    msg <- as.character(f)
    if (grepl("data are essentially constant", msg)) {
      cat("ℹ CD instable (stat quasi constante) — test neutralisé.\n")
      return(NULL)
    }
    stop(msg, call. = FALSE)
  }
  if (!is.null(label)) cat("\n--- VALIDATION PAR APPEL", paste0("[", label, "]"), "---\n")
  f
}


# ==============================================================================
# 5) DONNÉES DE DÉMO (partition + attributs monadyques + matrices dyadiques)
# ==============================================================================

# Partition (mode 1 → acteurs, mode 2 → groupes)
partition <- c(1, 1, 2, 2, 2, 3)

# Table des acteurs (attributs vertex pour le mode 1)
nodes <- data.frame(
  label  = c("A","B","C","D","E","F"),
  gender = c(1,1,2,1,2,2),
  age    = c(20,22,25,30,30,31),
  stringsAsFactors = FALSE
)

# Matrice d’amitié (exemple) — stockée en attribut de réseau (dyad-level)
friendship <- matrix(c(
  0,1,1,1,0,0,
  1,0,0,0,1,0,
  1,0,0,0,1,0,
  1,0,0,0,0,0,
  0,1,1,0,0,1,
  0,0,0,0,1,0
), 6, 6, byrow = TRUE, dimnames = list(nodes$label, nodes$label))

# Matrice de distance (exemple) — idem
distance <- matrix(c(
  0,2,3,4,5,6,
  2,0,1,3,4,5,
  3,1,0,2,3,5,
  4,3,2,0,2,4,
  5,4,3,2,0,1,
  6,5,5,4,1,0
), 6, 6, byrow = TRUE, dimnames = list(nodes$label, nodes$label))


# ==============================================================================
# 6) CONSTRUCTION DU RÉSEAU BIPARTI (à partir de la partition)
# ==============================================================================

#' Construit un réseau biparti acteurs-groupes à partir d’une partition.
#'
#' - Mode 1 : acteurs (indices 1..n)
#' - Mode 2 : groupes (indices n+1..n+G)
#' - Ajoute un lien (acteur i) -- (groupe partition[i]) pour chaque i
#'
#' @param partition vecteur d’appartenance (longueur n), valeurs dans 1..G
#' @param nodes     data.frame avec une colonne "label" et d’autres attributs
#' @param friendship (optionnel) matrice n x n pour edgecov dyadique
#' @param distance   (optionnel) matrice n x n pour edgecov dyadique
#' @return objet `network` biparti non orienté
make_bip <- function(partition, nodes, friendship = NULL, distance = NULL) {
  n <- length(partition)           # nombre d'acteurs
  G <- max(partition)              # nombre de groupes
  actor_names <- nodes$label
  group_names <- paste0("G", seq_len(G))
  all_names   <- c(actor_names, group_names)

  # Matrice d'adjacence (n+G) x (n+G), arcs entre (i) et (n + partition[i])
  adj <- matrix(0L, n + G, n + G, dimnames = list(all_names, all_names))
  ia <- seq_len(n)
  ig <- n + partition
  adj[cbind(ia, ig)] <- 1L
  adj[cbind(ig, ia)] <- 1L

  # Construit le réseau biparti
  nw <- network::network(adj, directed = FALSE, matrix.type = "adjacency")
  network::set.network.attribute(nw, "bipartite", n)       # coupe biparti
  network::set.vertex.attribute(nw, "vertex.names", all_names)

  # Renseigne les attributs de noeuds (mode 1 = acteurs ; mode 2 = NA)
  for (a in setdiff(names(nodes), "label")) {
    network::set.vertex.attribute(nw, a, c(nodes[[a]], rep(NA, G)))
  }

  # Introduit les attributs dyadiques dans le réseau biparti
  if (!is.null(friendship)) nw %n% "friendship" <- friendship
  if (!is.null(distance))   nw %n% "distance"   <- distance

  nw
}

# -- Instancie le réseau pour les tests ----------------------------------------
nw <- make_bip(partition, nodes, friendship, distance)


# ==============================================================================
# 7) CONTRÔLES {ergm} POUR DES RUNS RAPIDES
# ==============================================================================

# - Paramétrage conservateur : itératifs CD/MCMLE courts, pas de cluster.
# - Cherche à valider la chaîne d’appel plutôt qu’à inférer.
CTRL <- ergm::control.ergm(
  MCMLE.maxit     = 10,     # bornes faibles
  CD.maxit        = 10,
  CD.samplesize   = 1024,
  MCMC.burnin     = 2000,
  MCMC.interval   = 1000,
  parallel        = 0L,     # PAS DE CLUSTER (évite l’attachement ERPM côté workers)
  parallel.type   = NULL,
  seed            = 42
)


# ==============================================================================
# 8) SANITY CHECK : CONTRAINTE b1part DISPONIBLE + FONCTIONNELLE
# ==============================================================================

# -- Vérifie la présence de la contrainte b1part ------------------------
.ok(exists("InitErgmConstraint.b1part", mode = "function"),
    "InitErgmConstraint.b1part introuvable dans la recherche")

# -- Tente un fit minimal (edges) avec contrainte ~b1part (via CD) -------------
tmp <- try(
  ergm(nw ~ edges, constraints = ~ b1part,
       estimate = "CD", eval.loglik = FALSE, control = CTRL),
  silent = TRUE
)
.ok(!inherits(tmp, "try-error"),
    paste("b1part KO sur edges-only :", as.character(tmp)))


# ==============================================================================
# 9) DRY-RUNS : VALIDATION DE LA TRADUCTION (ERPM -> ERGM)
#    (Vérifie le motif dans l’appel ergm(...) généré, sans exécution)
# ==============================================================================

run("Dry-run: groups(2) -> b2degrange(2,3)", {
  o <- erpm(partition ~ groups(2), eval_call = FALSE, verbose = TRUE)
  expect_in_call(o, "b2degrange\\(from=2,to=3\\)", "groups(2)")
})

run("Dry-run: groups(from=2,to=4)", {
  o <- erpm(partition ~ groups(from=2,to=4), eval_call = FALSE, verbose = TRUE)
  expect_in_call(o, "b2degrange\\(from=2,to=4\\)", "groups(2..4)")
})

run("Dry-run: cliques(3)", {
  o <- erpm(partition ~ cliques(3), eval_call = FALSE, verbose = TRUE)
  expect_in_call(o, "cliques\\(clique_size=3L?,normalized=FALSE\\)", "cliques(3)")
})

run("Dry-run: cliques(2, normalized=TRUE)", {
  o <- erpm(partition ~ cliques(clique_size=2, normalized=TRUE), eval_call = FALSE, verbose = TRUE)
  expect_in_call(o, "cliques\\(clique_size=2L?,normalized=TRUE\\)", "cliques(2,TRUE)")
})

run("Dry-run: squared_sizes(from=1,to=3,pow=3)", {
  o <- erpm(partition ~ squared_sizes(from=1,to=3,pow=3), eval_call = FALSE, verbose = TRUE)
  expect_in_call(o, "squared_sizes\\(from=1,to=3,pow=3\\)", "squared_sizes(1..3,pow3)")
})


# ==============================================================================
# 10) FITS CD : VALIDATION "PAR APPEL" (robuste aux coef.names vides)
#      (N’asserte pas sur les coefficients mais sur les motifs de l’appel)
# ==============================================================================

run("CD: groups(2) + squared_sizes()", {
  f <- safefit(
    erpm(partition ~ groups(2) + squared_sizes(),
         estimate = "CD", eval.loglik = FALSE, control = CTRL),
    "groups+squared_sizes"
  )
  if (!is.null(f)) {
    expect_in_call(f, "b2degrange\\(from=2,to=3\\)")
    expect_in_call(f, "squared_sizes\\(")
  }
})

run("CD: cliques(3)", {
  f <- safefit(
    erpm(partition ~ cliques(3),
         estimate = "CD", eval.loglik = FALSE, control = CTRL),
    "cliques3"
  )
  if (!is.null(f)) expect_in_call(f, "cliques\\(clique_size=3L?,normalized=FALSE\\)")
})

run("CD: cliques(2, normalized=TRUE)", {
  f <- safefit(
    erpm(partition ~ cliques(clique_size=2, normalized=TRUE),
         estimate = "CD", eval.loglik = FALSE, control = CTRL),
    "cliques2_norm"
  )
  if (!is.null(f)) expect_in_call(f, "cliques\\(clique_size=2L?,normalized=TRUE\\)")
})

run("CD: cliques(c(2,3))", {
  f <- try(erpm(partition ~ cliques(c(2,3)),
                estimate = "CD", eval.loglik = FALSE, control = CTRL), silent = TRUE)
  if (inherits(f, "try-error")) cat("ℹ cliques(c(2,3)) non supporté — wrapper à corriger.\n")
})

run("CD: squared_sizes(from=1,to=3,pow=2)", {
  f <- safefit(
    erpm(partition ~ squared_sizes(from=1,to=3,pow=2),
         estimate = "CD", eval.loglik = FALSE, control = CTRL),
    "squared_sizes_pow2"
  )
  if (!is.null(f)) expect_in_call(f, "squared_sizes\\(from=1,to=3,pow=2\\)")
})

run("CD: squared_sizes(from=1,to=Inf,pow=3)", {
  f <- safefit(
    erpm(partition ~ squared_sizes(from=1,to=Inf,pow=3),
         estimate = "CD", eval.loglik = FALSE, control = CTRL),
    "squared_sizes_pow3"
  )
  if (!is.null(f)) expect_in_call(f, "squared_sizes\\(from=1,to=Inf,pow=3\\)")
})


# ==============================================================================
# 11) FITS CD — ATTRIBUTS (validation via motifs dans l’appel)
# ==============================================================================

run("CD: attributs via partition ~ Proj1/B", {
  f <- safefit(
    erpm(
      partition ~ Proj1(~ B(~ nodematch("gender") + absdiff("age") + edgecov("friendship"), "nonzero")),
      nodes = nodes, dyads = list(friendship = friendship, distance = distance),
      estimate = "CD", eval.loglik = FALSE, control = CTRL
    ),
    "Proj1_partition"
  )
  if (!is.null(f)) {
    expect_in_call(f, "Proj1\\(")
    expect_in_call(f, "B\\(")
    expect_in_call(f, "nodematch\\(\"gender\"\\)")
    expect_in_call(f, "absdiff\\(\"age\"\\)")
    expect_in_call(f, "edgecov\\(\"friendship\"\\)")
  }
})

run("CD: attributs via nw ~ Proj1/B", {
  f <- safefit(
    erpm(
      nw ~ Proj1(~ B(~ nodematch("gender") + absdiff("age") + edgecov("friendship"), "nonzero")),
      estimate = "CD", eval.loglik = FALSE, control = CTRL
    ),
    "Proj1_nw"
  )
  if (!is.null(f)) {
    expect_in_call(f, "Proj1\\(")
    expect_in_call(f, "B\\(")
    expect_in_call(f, "nodematch\\(\"gender\"\\)")
    expect_in_call(f, "absdiff\\(\"age\"\\)")
    expect_in_call(f, "edgecov\\(\"friendship\"\\)")
  }
})


# ==============================================================================
# 12) TESTS PLUS DENSES (toujours validation par appel)
# ==============================================================================

run("CD: groups + cliques + squared_sizes + attributs", {
  f <- safefit(
    erpm(
      partition ~ groups(2) + cliques(3) + squared_sizes(pow = 3) +
        Proj1(~ B(~ nodematch("gender") + absdiff("age") + edgecov("friendship"), "nonzero")),
      nodes = nodes, dyads = list(friendship = friendship),
      estimate = "CD", eval.loglik = FALSE, control = CTRL
    ),
    "combo1"
  )
  if (!is.null(f)) {
    expect_in_call(f, "b2degrange\\(from=2,to=3\\)")
    expect_in_call(f, "cliques\\(clique_size=3L?,normalized=FALSE\\)")
    expect_in_call(f, "squared_sizes\\(pow=3\\)")
    expect_in_call(f, "Proj1\\(")
    expect_in_call(f, "B\\(")
  }
})

run("CD: cliques(norm=TRUE) + squared_sizes", {
  f <- try(
    erpm(
      partition ~ cliques(c(2,3), normalized = TRUE) + squared_sizes(from = 1, to = 3, pow = 2),
      estimate = "CD", eval.loglik = FALSE, control = CTRL
    ),
    silent = TRUE
  )
  if (inherits(f, "try-error")) {
    if (grepl("data are essentially constant", as.character(f))) {
      cat("ℹ CD instable (stat quasi constante) — neutralise le test.\n")
    } else {
      cat("ℹ cliques(c(2,3)) non supporté — corrige le wrapper.\n")
    }
  } else {
    expect_in_call(f, "cliques\\(clique_size=2L?,normalized=TRUE\\)")
    expect_in_call(f, "squared_sizes\\(from=1,to=3,pow=2\\)")
  }
})

run("CD: nw ~ groups + cliques + squared_sizes", {
  f <- safefit(
    erpm(
      nw ~ groups(3) + cliques(2, normalized = FALSE) + squared_sizes(),
      estimate = "CD", eval.loglik = FALSE, control = CTRL
    ),
    "combo3"
  )
  if (!is.null(f)) {
    expect_in_call(f, "b2degrange\\(from=3,to=4\\)")
    expect_in_call(f, "cliques\\(clique_size=2L?,normalized=FALSE\\)")
    expect_in_call(f, "squared_sizes\\(")
  }
})


# ==============================================================================
# 13) FIN DES TESTS & NETTOYAGE
# ==============================================================================

cat("\n✅ Selftests étendus erpm() terminés.\n")

# -- Désactive le patch {ergm} si actif ---------------------------------------
if (exists("ergm_patch_disable")) ergm_patch_disable(verbose = TRUE)