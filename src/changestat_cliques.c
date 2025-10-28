// changestats_cliques.c
// Term: cliques — exemple composite minimaliste pour illustrer la
//        variation combinée "triangles + 2-stars" pondérée par un unique poids.
// ATTENTION: ce terme n’est pas un clone fidèle d’un terme {ergm}. Il sert
//            uniquement à corriger les erreurs de compilation et de signe
//            observées dans le code fourni.
// Entrée (INPUT_PARAM):
//   [0] = w  (poids double appliqué à toute la variation)
// Définition opérationnelle:
//   Δ = w * ( ΔTriangles(head,tail) + ΔTwoStars(head,tail) ),
//   où ΔTriangles = nb de voisins communs de head et tail avec signe du toggle,
//       ΔTwoStars = (deg(head)-edgestate) + (deg(tail)-edgestate) avec signe.
// Complexité: O(deg(head) + deg(tail)).
#include "ergm_changestat.h"
#include "ergm_storage.h"

C_CHANGESTAT_FN(c_cliques){
  const double w = (N_INPUT_PARAMS >= 1) ? INPUT_PARAM[0] : 1.0;
  const int add = edgestate ? -1 : +1; // +1 si on ajoute, -1 si on retire

  // --- Triangles (voisins communs) ---
  // On parcourt les voisins de head et on teste la présence d'un lien avec tail.
  // Pour des graphes non orientés, IS_UNDIRECTED_EDGE est adapté.
  // Pour orientés, ce n'est qu'une approximation.
  Vertex n3;
  Edge e;
  int commons = 0;

  STEP_THROUGH_OUTEDGES(head, e, n3){ commons += IS_UNDIRECTED_EDGE(n3, tail); }
  STEP_THROUGH_INEDGES(head,  e, n3){ commons += IS_UNDIRECTED_EDGE(n3, tail); }

  double d = 0.0;
  d += w * add * (double)commons;

  // --- 2-stars centrées sur head et tail ---
  // Variation pour un toggle: (deg(head) + deg(tail)) avec signe add.
  const int taildeg = (int)(OUT_DEG[tail] + IN_DEG[tail] - edgestate);
  const int headdeg = (int)(OUT_DEG[head] + IN_DEG[head] - edgestate);
  d += w * add * (double)(taildeg + headdeg);

  CHANGE_STAT[0] += d;
}
