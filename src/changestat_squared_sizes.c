// changestat_squared_sizes.c
// Term: squared_sizes — variation de la somme, sur des sommets dont le degré total
//        est dans [from, to), de deg^power, lors du toggle d’une arête.
// Entrée (INPUT_PARAM par statistique j):
//   [3*j + 0] = from  (entier, inclusif)
//   [3*j + 1] = to    (entier, exclusif)
//   [3*j + 2] = power (entier >= 1)
// Sortie:
//   CHANGE_STAT[j] += Δ, où Δ = Σ_sommets affectés ( f(new_deg) - f(old_deg) )
//   avec f(d) = d^power si d ∈ [from, to), sinon 0.
// Hypothèses:
//   - Graphes éventuellement dirigés. Le "degré total" = IN_DEG + OUT_DEG.
//   - Le toggle affecte exactement les deux sommets tail et head.
// Complexité: O(N_CHANGE_STATS).
#include <math.h>
#include "ergm_changestat.h"
#include "ergm_storage.h"

#define IN_RANGE(x,a,b) ((x) >= (a) && (x) < (b))

static inline double ipow_int(int base, int exp){
  // puissance entière sans erreurs d’arrondi, exp >= 0
  if(exp <= 0) return 1.0;
  double r = 1.0, b = (double)base;
  int e = exp;
  while(e){
    if(e & 1) r *= b;
    b *= b;
    e >>= 1;
  }
  return r;
}

C_CHANGESTAT_FN(c_squared_sizes){
  // delta = +1 si on ajoute l’arête, -1 si on la retire
  const int delta = edgestate ? -1 : +1;

  // degrés totaux actuels
  const int taildeg_old = (int)(OUT_DEG[tail] + IN_DEG[tail]);
  const int headdeg_old = (int)(OUT_DEG[head] + IN_DEG[head]);

  // degrés après toggle
  const int taildeg_new = taildeg_old + delta;
  const int headdeg_new = headdeg_old + delta;

  for(int j = 0; j < N_CHANGE_STATS; ++j){
    const int from  = (int)INPUT_PARAM[3*j + 0];
    const int to    = (int)INPUT_PARAM[3*j + 1];
    const int power = (int)INPUT_PARAM[3*j + 2];

    double d = 0.0;

    // contribution du tail
    if(IN_RANGE(taildeg_new, from, to)) d += ipow_int(taildeg_new, power);
    if(IN_RANGE(taildeg_old, from, to)) d -= ipow_int(taildeg_old, power);

    // contribution du head
    if(IN_RANGE(headdeg_new, from, to)) d += ipow_int(headdeg_new, power);
    if(IN_RANGE(headdeg_old, from, to)) d -= ipow_int(headdeg_old, power);

    CHANGE_STAT[j] += d;
  }
}
