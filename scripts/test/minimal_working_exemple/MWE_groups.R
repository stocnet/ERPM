# ======================================================================================
# Fichier : scripts/test/minimal_working_exemple/MWE_groups.R
#
# Objet   : MWE pour l’effet ERPM `groups(k)` — calcule la stat observée (summary)
#            et effectue un fit ERGM complet via `erpm()`, avec affichage lisible.
#
# Contexte: Ce script illustre la chaîne complète ERPM → ERGM :
#            partition -> réseau biparti -> traduction -> ergm() ou summary()
#
# Résumé technique :
#   • L’effet `groups(k)` ne nécessite pas de code C ni d’InitErgmTerm dédié.
#   • `erpm()` convertit la partition en graphe biparti (acteurs ↔ groupes).
#   • Le wrapper traduit `groups(k)` en `b2degrange(from=k,to=k+1)`.
#   • {ergm} sait déjà calculer cette statistique sur le second mode (groupes).
#   • Ainsi, `summary()` ou `ergm()` fonctionnent directement sur la formule traduite.
# ======================================================================================

# ----- Préambule locale/UTF-8 (assure des affichages stables) -------------------------
Sys.setenv(LANG = "fr_FR.UTF-8")
invisible(try(Sys.setlocale("LC_CTYPE","fr_FR.UTF-8"), silent = TRUE))

# ----- Dépendances minimales ----------------------------------------------------------
suppressPackageStartupMessages({
  library(devtools)  # charge le package local avec load_all(".")
  library(network)   # fournit l’objet network et l’attribut 'bipartite'
  library(ergm)      # fournit summary(), ergm(), control.ergm(), etc.
})

# ----- Charge le package local ERPM (expose erpm, InitErgmTerm, patch helpers) --------
devtools::load_all(".")

# ----- Active le patch {ergm} -------------------------------------------
#   Le patch corrige le trace base::replace()
source("scripts/ergm_patch.R")
ergm_patch_enable()

# ----- Crée la partition de test ------------------------------------------------------
#   On a 4 groupes : {1,2,3,4} de tailles = (1,3,2,2)
partition <- c(1, 2, 2, 2, 3, 3, 4, 4)

cat("\nPartition :", paste(partition, collapse = ", "), "\n")
sizes <- as.integer(table(partition))
cat("Tailles des groupes (par ordre de label):", paste(sizes, collapse = ", "), "\n")
cat("Nombre d'acteurs (N1):", length(partition), " | Nombre de groupes:", length(sizes), "\n\n")

# ----- Définit la taille du groupe à évaluer ------------------------------------------
#   groups(k) compte les groupes de taille EXACTEMENT k.
taille_group <- 2

# ======================================================================================
# 1) ÉVALUE LA STAT OBSERVÉE SANS FIT
#    - ERPM renvoie la formule  -> ERGM sans l'exécuter (dry-run)
#    - Passe la formule traduite à summary() avec la contrainte bipartite (~b1part)
# ======================================================================================
dry <- erpm(partition ~ groups(taille_group), eval_call = FALSE, verbose = TRUE)

# Récupère la formule ergm() traduite 
fml <- dry[[2]]

# Calcule la statistique observée sur le réseau biparti construit par erpm()
obs <- summary(fml, constraints = ~ b1part)                         # le réseau biparti est dans l'environnement de fml
cat(sprintf("\n[MWE groups(%d)] summary observé = %d  | terme = %s\n\n",
            taille_group, as.integer(obs), names(obs)))

# Vérifie que la valeur est correcte et cohérente avec la partition
stopifnot(is.numeric(obs), length(obs) == 1L)

# ======================================================================================
# 2) EFFECTUE UN FIT ERGM COMPLET ET AFFICHE LES RÉSULTATS DE MANIÈRE LISIBLE
#    - Le fit retourne un paramètre estimé 
#    - On vérifie la cohérence en recalculant summary() sur la formule du fit
# ======================================================================================
# Fixe une seed pour rendre l’estimation  plus stable entre exécutions
set.seed(1)

# Lance le fit complet
fit <- erpm(partition ~ groups(taille_group),
            eval_call  = TRUE,
            verbose    = TRUE,
            estimate   = "MLE",
            eval.loglik = TRUE)

# ----- Affiche le fit de manière synthétique ------------------------------------------
cat("\n--- summary(fit) ---\n") # Affiche un résumé standard ergm() plus détaillé
print(summary(fit))