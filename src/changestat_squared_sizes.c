// =============================================================================
/**
 * @file changestat_squared_sizes.c
 * @brief Change statistic C pour le terme `squared_sizes` (forme "un-toggle")
 *
 * @details
 * Calcule la variation de la somme, **sur les nœuds du mode 2 (groupes)**,
 * de `(deg^pow)` pour les nœuds dont le degré total `deg ∈ [from, to)`,
 * lors **du traitement d’un toggle** (ajout ou retrait d’arête).
 *
 * ## Contexte {ergm} (forme un-toggle)
 * - La macro `C_CHANGESTAT_FN(f)` déclare une fonction C de signature :
 *   ```
 *   void f(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate)
 *   ```
 *   Chaque appel traite **un seul toggle** `(tail, head)` de l’état courant.
 *   - `edgestate = 1`  → l’arête existe (toggle = retrait)
 *   - `edgestate = 0`  → l’arête n’existe pas (toggle = ajout)
 *
 * - À **chaque appel**, `ZERO_ALL_CHANGESTATS(...)` remet `CHANGE_STAT[0..N_CHANGE_STATS-1]` à zéro.
 *   Le moteur {ergm} **accumule ensuite** les contributions de tous les appels pour produire
 *   la statistique finale (summary ou MCMC).
 *
 * ## Entrée (depuis R → C, via InitErgmTerm.squared_sizes)
 * - `N_CHANGE_STATS` :
 *    nombre de statistiques élémentaires renvoyées par ce terme.
 *    - Défini automatiquement par la longueur de `coef.names` (créée dans `InitErgmTerm`).
 *    - Exemple :
 *      - un seul triplet `(from,to,pow)` ⇒ `N_CHANGE_STATS = 1`
 *      - plusieurs intervalles `(from,to)` ⇒ `N_CHANGE_STATS = length(from)`
 *
 * - `INPUT_PARAM` :
 *    vecteur mémoire contigu de taille `3*N_CHANGE_STATS` contenant, pour chaque sous-terme j :
 *    ```
 *    [3*j + 0] = from  (entier inclusif)
 *    [3*j + 1] = to    (entier exclusif)
 *    [3*j + 2] = pow   (entier >= 1)
 *    ```
 *    (Construit côté R par :
 *    `inputs <- c(rbind(as.integer(from), as.integer(to), as.integer(pow)))`)
 *
 * ## Sortie (par appel)
 * - Écrit dans `CHANGE_STAT[j]` la variation Δⱼ due **au toggle courant** :
 *   ```
 *   Δⱼ = f(deg_new) - f(deg_old),
 *   f(d) = d^pow si d ∈ [from,to), sinon 0
 *   ```
 *   - On ne met à jour que le sommet **côté mode 2 (groupes)** impacté.
 *   - {ergm} additionne ensuite tous les vecteurs produits pour calculer la statistique finale.
 *
 * ## Hypothèses
 * - Réseau biparti : le mode 1 (acteurs) = `1..BIPARTITE`, mode 2 (groupes) = `(BIPARTITE+1)..N`.
 * - Le toggle connecte un acteur et un groupe (un seul groupe est modifié à chaque fois).
 * - Le degré total utilisé est `OUT_DEG[v] + IN_DEG[v]` (représentation interne non dirigée).
 *
 * ## Complexité
 * - O(N_CHANGE_STATS) par appel (un-toggle)
 *
 * ## Remarques
 * - `ZERO_ALL_CHANGESTATS(...)` doit être appelé **au début** : chaque appel ne renvoie
 *   que sa contribution locale, pas la somme globale.
 * - `edgestate ? -1 : +1` :
 *   - si l’arête existe ⇒ retrait ⇒ degré -1
 *   - sinon ⇒ ajout ⇒ degré +1
 * - L’intervalle `[from,to)` est inclusif/exclusif strict.
 * - Validation (bornes, bipartition, Inf→n+1) faite côté R.
 *
 * ## Journalisation
 * - `Rprintf(...)` est utilisé ici pour le debug à chaud, visible pendant `summary()` / MCMC.
 */
// =============================================================================

#include <math.h>
#include <R_ext/Print.h>        // Rprintf
#include "ergm_changestat.h"
#include "ergm_storage.h"

// Test d’appartenance à l’intervalle [a,b)
#define IN_RANGE(x,a,b) ((x) >= (a) && (x) < (b))

// Si vaut 1, permet de print dans la console
#define DEBUG_SQUARED_SIZES 0

// -----------------------------------------------------------------------------
// Exponentiation rapide 
//   Idée : écriture binaire de exp, boucle en O(log exp) (bien plus rapide
//   et stable qu’une boucle naïve en O(exp)).
// -----------------------------------------------------------------------------
static inline double ipow_int(int base, int exp){
  if(exp <= 0) return 1.0;
  double r = 1.0, b = (double)base;
  while(exp){
    if(exp & 1) r *= b;  // si le bit courant est 1, on multiplie le résultat
    b *= b;              // on met à jour la puissance en carrant
    exp >>= 1;           // on passe au bit suivant
  }
  return r;
}

/*
 * Forme C "un-toggle" : traite UN toggle (tail, head, edgestate).
 *
 * Signature imposée par {ergm} :
 *   #define C_CHANGESTAT_FN(a) \
 *     void a (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate)
 *
 * Paramètres :
 *   - tail, head  : les sommets du toggle courant (indices 1..N).
 *   - mtp         : descripteur du terme (inclut nstats, inputparams, etc.).
 *   - nwp         : pointeur vers le réseau (degrés, voisinages, ...).
 *   - edgestate   : 1 si l’arête (tail,head) existe déjà (toggle = retrait), 0 sinon (toggle = ajout).
 *
 * Effets :
 *   - Remplit CHANGE_STAT[0..N_CHANGE_STATS-1] avec la contribution de CE toggle uniquement.
 *   - Ne modifie pas explicitement les degrés : on calcule deg_new = deg_old ± 1 localement.
 *     (Dans la forme multi-toggle, les macros TOGGLE/UNDO mettent à jour nwp ; ici non.)
 */
C_CHANGESTAT_FN(c_squared_sizes){

  // --- 1) Toujours réinitialiser le tampon de sortie pour CET appel. ---
  //     {ergm} fera l’accumulation globale en dehors de cette fonction.
  ZERO_ALL_CHANGESTATS(0);  // memset(CHANGE_STAT, 0, nstats*sizeof(double))

  // --- 2) Taille du mode 1 (acteurs). Si non-biparti, BIPARTITE vaut 0. ---
  const int n1 = BIPARTITE;

  // --- 3) Lecture des extrémités du toggle passé par {ergm}. ---
  Vertex t = tail, h = head;

  // --- 4) Sélection du nœud côté groupes (mode 2) : indices > n1. ---
  //     Le toggle peut être dans n’importe quel sens : on choisit le sommet > n1.
  Vertex v2 = (t > n1) ? t : h;

  // --- 5) Degré total du groupe avant/après toggle (calcul local pour la delta-stat). ---
  //     OUT_DEG/IN_DEG : pointeurs sur nwp->outdegree / nwp->indegree.
  int deg_old = (int)(OUT_DEG[v2] + IN_DEG[v2]);
  int delta   = edgestate ? -1 : +1;   // si l’arête existe, on la retire → -1 ; sinon on l’ajoute → +1
  int deg_new = deg_old + delta;

  // --- 6) Logs debug 
  #if DEBUG_SQUARED_SIZES
    Rprintf("[C:c_squared_sizes] tail=%d head=%d | v2=%d (mode2) | edgestate=%d | deg_old=%d -> deg_new=%d | N_STATS=%d\n",
            (int)t, (int)h, (int)v2, (int)edgestate, deg_old, deg_new, N_CHANGE_STATS);
  #endif

  // --- 7) Pour chaque sous-terme vectorisé j : lire (from,to,pow) et accumuler Δ_j. ---
  for(int j = 0; j < N_CHANGE_STATS; ++j){
    // INPUT_PARAM = mtp->inputparams (double*) construit côté R (aplati par colonnes)
    int from  = (int)INPUT_PARAM[3*j + 0];
    int to    = (int)INPUT_PARAM[3*j + 1];
    int power = (int)INPUT_PARAM[3*j + 2];

    double d = 0.0;
    if(IN_RANGE(deg_new, from, to)) d += ipow_int(deg_new, power);
    if(IN_RANGE(deg_old, from, to)) d -= ipow_int(deg_old, power);

    // On additionne dans le buffer DE CE SEUL appel.
    // Le moteur {ergm} additionnera ensuite ce buffer avec ceux des autres appels.
    CHANGE_STAT[j] += d;

    // --- 8) Logs debug 
    #if DEBUG_SQUARED_SIZES
      Rprintf("  stat[%d]: from=%d to=%d pow=%d | Δ=%.2f | cumul=%.2f\n",
              j, from, to, power, d, CHANGE_STAT[j]);
    #endif
  }
}
