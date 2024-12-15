#define STRICT_Wt_HEADERS
#include "ergm_wtchangestat_operator.h"
#include "ergm_changestat_operator.h"
#include "ergm_storage.h"

typedef struct {
  WtModel *m;
  Vertex *t, *h;
  double *w;
} StoreWtModelAndWtChanges;

// Proj1: Mode 1 projection
I_CHANGESTAT_FN(i_on_proj1_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  ALLOC_STORAGE(1, StoreWtModelAndWtChanges, storage);
  // No need to allocate it: we are only storing a pointer to a model.
  storage->m = WtModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state, pnwp, FALSE);

  /* WtSELECT_C_OR_D_BASED_ON_SUBMODEL(storage->m); */
  WtDELETE_IF_UNUSED_IN_SUBMODEL(x_func, storage->m);
  WtDELETE_IF_UNUSED_IN_SUBMODEL(z_func, storage->m);

  storage->t = R_Calloc(BIPARTITE-1, Vertex);
  storage->h = R_Calloc(BIPARTITE-1, Vertex);
  storage->w = R_Calloc(BIPARTITE-1, double);
}

/* D_CHANGESTAT_FN(d_on_proj1_net){ */
/*   GET_AUX_STORAGE(WtNetwork, pnwp); */
/*   GET_STORAGE(WtModel, m); */

/*   WtChangeStats(ntoggles, tails, heads, weights, pnwp, m); */

/*   memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double)); */
/* } */

C_CHANGESTAT_FN(c_on_proj1_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  GET_STORAGE(StoreWtModelAndWtChanges, storage);

  int echange = edgestate ? -1 : 1;

  unsigned int nt = 0;
  EXEC_THROUGH_FINEDGES(head, e, tail2, {
      if(tail!=tail2){
        storage->t[nt] = MIN(tail, tail2);
        storage->h[nt] = MAX(tail, tail2);
        storage->w[nt] = WtGETWT(tail, tail2, pnwp) + echange;
        nt++;
      }
    });

  WtChangeStats(nt, storage->t, storage->h, storage->w, pnwp, storage->m);

  memcpy(CHANGE_STAT, storage->m->workspace, N_CHANGE_STATS*sizeof(double));
}

Z_CHANGESTAT_FN(z_on_proj1_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  GET_STORAGE(StoreWtModelAndWtChanges, storage)

  WtZStats(pnwp, storage->m, FALSE);

  memcpy(CHANGE_STAT, storage->m->workspace, N_CHANGE_STATS*sizeof(double));
}

X_CHANGESTAT_FN(x_on_proj1_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  GET_STORAGE(StoreWtModelAndWtChanges, storage);
  ModelTerm *_mymtp = mtp;
  WtSEND_X_SIGNAL_INTO(pnwp, storage->m, NULL, _mymtp->dstats, type, data);
}

F_CHANGESTAT_FN(f_on_proj1_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  GET_STORAGE(StoreWtModelAndWtChanges, storage);

  R_Free(storage->t);
  R_Free(storage->h);
  R_Free(storage->w);
  WtModelDestroy(pnwp, storage->m);
}



I_CHANGESTAT_FN(i__proj1_net){
  WtNetwork *pnwp = AUX_STORAGE = WtNetworkInitialize(NULL, NULL, NULL, 0, BIPARTITE, DIRECTED, FALSE, FALSE, 0, NULL);

  EXEC_THROUGH_NET_EDGES_PRE(tail1, head, e1, {
      EXEC_THROUGH_FINEDGES(head, e2, tail2, {
          if(tail1 < tail2)
            WtSETWT(tail1, tail2, WtGETWT(tail1, tail2, pnwp) + 1, pnwp);
        });
    });
}

U_CHANGESTAT_FN(u__proj1_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  int echange = edgestate ? -1 : 1;

  EXEC_THROUGH_FINEDGES(head, e, tail2, {
      if(tail!=tail2)
	WtSETWT(tail, tail2, WtGETWT(tail, tail2, pnwp) + echange, pnwp);
    });
}

F_CHANGESTAT_FN(f__proj1_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  WtNetworkDestroy(pnwp);
  AUX_STORAGE = NULL;
}
