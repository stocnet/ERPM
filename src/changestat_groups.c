// c_groups.c
// Term: groups — variation du nombre de sommets dont le degré total
//                appartient à un intervalle [from, to) (compte 1 par sommet).
// Remplace l'ancienne version boguée (sizeof(id)*pow, etc.).
// Entrée (INPUT_PARAM par statistique j):
//   [2*j + 0] = from  (entier, inclusif)
//   [2*j + 1] = to    (entier, exclusif)
// Sortie:
//   CHANGE_STAT[j] += Δ, où Δ = Σ_{v ∈ {tail, head}} ( 1_{new∈[from,to)} - 1_{old∈[from,to)} ).
// Complexité: O(N_CHANGE_STATS).
#include "ergm_changestat.h"
#include "ergm_storage.h"

#define IN_RANGE(x,a,b) ((x) >= (a) && (x) < (b))

C_CHANGESTAT_FN(c_groups){
  const int delta = edgestate ? -1 : +1;

  const int taildeg_old = (int)(OUT_DEG[tail] + IN_DEG[tail]);
  const int headdeg_old = (int)(OUT_DEG[head] + IN_DEG[head]);

  const int taildeg_new = taildeg_old + delta;
  const int headdeg_new = headdeg_old + delta;

  for(int j = 0; j < N_CHANGE_STATS; ++j){
    const int from = (int)INPUT_PARAM[2*j + 0];
    const int to   = (int)INPUT_PARAM[2*j + 1];

    double d = 0.0;

    if(IN_RANGE(taildeg_new, from, to)) d += 1.0;
    if(IN_RANGE(taildeg_old, from, to)) d -= 1.0;

    if(IN_RANGE(headdeg_new, from, to)) d += 1.0;
    if(IN_RANGE(headdeg_old, from, to)) d -= 1.0;

    CHANGE_STAT[j] += d;
  }
}
