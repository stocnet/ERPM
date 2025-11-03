# ====================================================================================== 
# Fichier : cub_test.R
# Fonction :
# Utilité : test pour les nouvelles fonctions ERPM
# ====================================================================================== 

# macOS/Linux
Sys.setenv(LANG="fr_FR.UTF-8")
invisible(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"))

# ====================================================================================== 
# ======================================== INIT ======================================== 
suppressMessages(library(ergm))

if (!exists(".__debug_loaded", envir = .GlobalEnv)) {
  source("scripts/local_source/debug.R", local = FALSE)
}

source("scripts/local_source/init.R", local = FALSE) # init l'environnement (source de settings et du launcher)
source("R/erpm_wrapper.R", local = FALSE) # n'est pas dans scripts/local_source/ donc on le charge indépendemment
# debug_source_file("R/erpm_wrapper.R") 


init_erpm(selftest=FALSE, verbose=FALSE)

if (exists("log_msg")) {
  log_msg("INFO", "Démarrage du script pour ERPM")
} else {
  message("Démarrage du script pour ERPM")
}

if (exists("ergm_patch_enable")) ergm_patch_enable( verbose = VERBOSE )


# ====================================================================================== 
# ======================================== RUN ========================================= 


# cases <- list(
#   # Acceptés
#   list(name = "clq_default",         effects = "cliques"),                                       # -> edges()
#   list(name = "clq_k2_F",            effects = "cliques(clique_size = 2, normalized = FALSE)"),  # -> edges()
#   list(name = "clq_k2_T",            effects = "cliques(clique_size = 2, normalized = TRUE)"),   # -> edges() + normalisation post
#   list(name = "clq_normF_defaultk",  effects = "cliques(normalized = FALSE)"),                   # -> edges()
#   list(name = "clq_normT_defaultk",  effects = "cliques(normalized = TRUE)"),                    # -> edges() + normalisation post
#   list(name = "clq_k2_only",         effects = "cliques(clique_size = 2)"),                      # -> edges()
#   list(name = "clq_k1_default",      effects = "cliques(clique_size = 1)"),                      # -> b2degrange(0,0)
#   list(name = "clq_k1_F",            effects = "cliques(clique_size = 1, normalized = FALSE)")   # -> b2degrange(0,0)
# )

# results <- list()

# for (cx in cases) {
#   log_msg("INFO", sprintf("ERPM run: %s -> %s", cx$name, cx$effects))
#   results[[cx$name]] <- launch_model(
#     engine      = "erpm",
#     effects     = cx$effects,
#     dry_run     = FALSE,
#     estimate    = "CD",
#     eval_loglik = TRUE,
#     control     = list(MCMLE.maxit = 3, MCMC.samplesize = 1000),
#     timeout     = NULL
#   )
# }

# Données
partition <- c(1, 1, 2, 2, 2, 3)

nodes <- data.frame(
  label  = c("A","B","C","D","E","F"),
  gender = c(1, 1, 2, 1, 2, 2),
  age    = c(20,22,25,30,30,31),
  stringsAsFactors = FALSE
)

friendship <- matrix(c(
  0,1,1,1,0,0,
  1,0,0,0,1,0,
  1,0,0,0,1,0,
  1,0,0,0,0,0,
  0,1,1,0,0,1,
  0,0,0,0,1,0
), 6, 6, byrow = TRUE, dimnames = list(nodes$label, nodes$label))

distance <- matrix(c(
  0,2,3,4,5,6,
  2,0,1,3,4,5,
  3,1,0,2,3,5,
  4,3,2,0,2,4,
  5,4,3,2,0,1,
  6,5,5,4,1,0
), 6, 6, byrow = TRUE, dimnames = list(nodes$label, nodes$label))

# --- Construire nw (avec attributs) ---
{
  # Vérifie que le package 'network' est bien disponible
  stopifnot(requireNamespace("network", quietly = TRUE))

  # n = nombre d’acteurs ; G = nombre de groupes distincts
  n <- length(partition)
  G <- max(partition)

  # Noms des acteurs et des groupes
  actor_names <- nodes$label
  group_names <- paste0("G", seq_len(G))
  all_names   <- c(actor_names, group_names)  # noms de tous les sommets

  # Initialise la matrice d’adjacence (0 partout)
  adj <- matrix(0L, n + G, n + G, dimnames = list(all_names, all_names))

  # Indices des acteurs et de leurs groupes
  ia <- seq_len(n)        # indices 1..n pour les acteurs
  ig <- n + partition     # indices des groupes associés dans la matrice

  # Relie chaque acteur à son groupe dans la matrice
  adj[cbind(ia, ig)] <- 1L
  adj[cbind(ig, ia)] <- 1L

  # Crée un objet 'network' biparti non orienté à partir de la matrice
  nw <- network::network(adj, directed = FALSE, matrix.type = "adjacency")

  # Spécifie la taille du premier mode (nombre d’acteurs)
  network::set.network.attribute(nw, "bipartite", n)

  # Enregistre les noms des sommets
  network::set.vertex.attribute(nw, "vertex.names", all_names)

  # Ajoute les attributs monadiques (gender, age, …)
  for (a in setdiff(names(nodes), "label")) {
    # On complète avec NA pour les sommets représentant les groupes
    network::set.vertex.attribute(nw, a, c(nodes[[a]], rep(NA, G)))
  }

  # Ajoute les matrices dyadiques comme attributs réseau (edgecov utilisera ça)
  nw %n% "friendship" <- friendship
  nw %n% "distance"   <- distance
}

# --- Construire nw2 (sans attributs) ---
{
  n <- length(partition)              # nb d'acteurs (mode 1)
  G <- max(partition)                 # nb de groupes distincts (mode 2)

  actor_names <- nodes$label          # noms des acteurs (A,B,C, …)
  group_names <- paste0("G", seq_len(G)) # noms des groupes (G1,G2,…,GG)
  all_names   <- c(actor_names, group_names) # noms de tous les sommets (acteurs puis groupes)

  # Matrice d’adjacence vide (0) de taille (n+G)×(n+G) avec noms de lignes/colonnes
  adj <- matrix(0L, n + G, n + G, dimnames = list(all_names, all_names))

  # Indices des acteurs dans 'all_names' (ici 1..n, mais on reste robuste)
  ia <- match(actor_names, all_names)

  # Indices, dans 'all_names', des groupes auxquels appartient chaque acteur
  # Exemple: partition = c(1,1,2,2,2,3)  →  "G1","G1","G2","G2","G2","G3"
  ig <- match(paste0("G", partition), all_names)

  # Pose les arêtes bipartites ACTEUR->GROUPE (demi-matrice supérieure)
  adj[cbind(ia, ig)] <- 1L

  # Symétrise (graphe non orienté) : GROUPE->ACTEUR
  adj[cbind(ig, ia)] <- 1L

  # Construit l’objet `network` non orienté à partir de la matrice d’adjacence
  nw2 <- network::network(adj, directed = FALSE, matrix.type = "adjacency")

  # Déclare le graphe biparti et indique la taille du premier mode (acteurs)
  network::set.network.attribute(nw2, "bipartite", n)

  # Enregistre les étiquettes des sommets (utile pour inspection/affichage)
  network::set.vertex.attribute(nw2, "vertex.names", all_names)
}


# A1. Groupes non vides
erpm(partition ~ groups + squared_sizes, estimate = "MLE", eval.loglik = FALSE)

# A2. Groupes de taille exactement 2
erpm(partition ~ groups(2))

# A3. Groupes de taille dans [2,4)
erpm(partition ~ groups(from = 2, to = 4))

# A4. Cliques k=2 (défaut)
erpm(partition ~ cliques())

# A5. Cliques k=3
erpm(partition ~ cliques(3))

# A6. Cliques normalisées k=2
erpm(partition ~ cliques(clique_size = 2, normalized = TRUE))

# A7. Somme des tailles^2
erpm(partition ~ squared_sizes())

# B1. Effets acteurs simples via Proj1/B
erpm(
  partition ~ Proj1(~ B(~ nodematch("gender") + absdiff("age") + edgecov("friendship"), "nonzero")),
  nodes = nodes, 
  dyads = list(friendship = friendship, distance = distance)
)

# B2. Modèle combiné: groups + cliques + squared_sizes + covariés acteurs
erpm(
  partition ~ groups() + cliques(2) + squared_sizes(pow=3) +
               Proj1(~ B(~ nodematch("gender") + absdiff("age") + edgecov("friendship"), "nonzero")),
  nodes = nodes, 
  dyads = list(friendship = friendship)
)

# C1. Avec nw (attributs disponibles)
erpm(
  nw ~ groups() + cliques(3) + squared_sizes() +
       Proj1(~ B(~ nodematch("gender") + absdiff("age") + edgecov("friendship"), "nonzero")),
  estimate = "CD", 
  eval.loglik = FALSE
)

# C2. Avec nw2 (aucun attribut → pas de covariés acteurs)
erpm(
  nw2 ~ groups() + cliques(2) + squared_sizes(),
  estimate = "CD", 
  eval.loglik = FALSE
)

# C3. Exemple minimal structurel sur nw
erpm(nw ~ groups() + cliques(3))

# C4. Exemple minimal structurel sur nw2
erpm(nw2 ~ squared_sizes() + groups(from = 2, to = 4))

# ======================================================================================= 
# ======================================== CLEAN ======================================== 

# Désactive le patch ERGM si actif
if (exists("ergm_patch_disable")) ergm_patch_disable( verbose = VERBOSE )

# Vide le buffer stdout
flush.console()  

# Log
if (exists("log_msg")) log_msg("INFO", "Fin du programme -- Nettoyage de l'environnement global")

# Nettoie l'environnement
clean_global_env(verbose = VERBOSE)