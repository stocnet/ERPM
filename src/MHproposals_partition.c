/*  File src/MHproposal_partition.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#include "ergm_MHproposals_degree.h"

MH_P_FN(MH_B1Part) {

  if (MHp->ntoggles == 0) { /* Initialize */
    MH_CondB1Degree(MHp, nwp);
    return;
  }

  MH_CondB1Degree(MHp, nwp);

  int dP = (IN_DEG[Mhead[1]] == 0) - (IN_DEG[Mhead[0]] == 1);

  if (dP) {
    Vertex P = 0;
    for (Vertex i = BIPARTITE + 1; i <= N_NODES; i++) {
      if (IN_DEG[i]) P++;
    }

    MHp->logratio += dP == -1 ? log(N_NODES - BIPARTITE - P + 1) : -log(N_NODES - BIPARTITE - P);
  }
}
