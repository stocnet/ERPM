#include "ergm_changestat.h"
#include "ergm_storage.h"

C_CHANGESTAT_FN(c_nodecovrange) {
  double oldmin = R_PosInf, oldmax = R_NegInf, newmin = R_PosInf, newmax = R_NegInf;

  EXEC_THROUGH_EDGES(tail, e, u, {
      oldmin = MIN(oldmin, INPUT_PARAM[u-1]);
      oldmax = MAX(oldmax, INPUT_PARAM[u-1]);

      if(!edgestate || u!=head){
        newmin = MIN(newmin, INPUT_PARAM[u-1]);
        newmax = MAX(newmax, INPUT_PARAM[u-1]);
      }
    });

  if(!edgestate){
    newmin = MIN(newmin, INPUT_PARAM[head-1]);
    newmax = MAX(newmax, INPUT_PARAM[head-1]);
  }

  CHANGE_STAT[0] += (isinf(newmax) ? 0 : newmax-newmin) - (isinf(oldmax) ? 0 : oldmax-oldmin);

  oldmin = R_PosInf; oldmax = R_NegInf; newmin = R_PosInf; newmax = R_NegInf;

  EXEC_THROUGH_EDGES(head, e, u, {
      oldmin = MIN(oldmin, INPUT_PARAM[u-1]);
      oldmax = MAX(oldmax, INPUT_PARAM[u-1]);

      if(!edgestate || u!=tail){
        newmin = MIN(newmin, INPUT_PARAM[u-1]);
        newmax = MAX(newmax, INPUT_PARAM[u-1]);
      }
    });

  if(!edgestate){
    newmin = MIN(newmin, INPUT_PARAM[tail-1]);
    newmax = MAX(newmax, INPUT_PARAM[tail-1]);
  }

  CHANGE_STAT[0] += (isinf(newmax) ? 0 : newmax-newmin) - (isinf(oldmax) ? 0 : oldmax-oldmin);
}



C_CHANGESTAT_FN(c_nodeocovrange) {
  double oldmin = R_PosInf, oldmax = R_NegInf, newmin = R_PosInf, newmax = R_NegInf;

  EXEC_THROUGH_OUTEDGES(tail, e, u, {
      oldmin = MIN(oldmin, INPUT_PARAM[u-1]);
      oldmax = MAX(oldmax, INPUT_PARAM[u-1]);

      if(!edgestate || u!=head){
        newmin = MIN(newmin, INPUT_PARAM[u-1]);
        newmax = MAX(newmax, INPUT_PARAM[u-1]);
      }
    });

  if(!edgestate){
    newmin = MIN(newmin, INPUT_PARAM[head-1]);
    newmax = MAX(newmax, INPUT_PARAM[head-1]);
  }

  CHANGE_STAT[0] += (isinf(newmax) ? 0 : newmax-newmin) - (isinf(oldmax) ? 0 : oldmax-oldmin);
}


C_CHANGESTAT_FN(c_nodeicovrange) {
  double oldmin = R_PosInf, oldmax = R_NegInf, newmin = R_PosInf, newmax = R_NegInf;

  EXEC_THROUGH_INEDGES(head, e, u, {
      oldmin = MIN(oldmin, INPUT_PARAM[u-1]);
      oldmax = MAX(oldmax, INPUT_PARAM[u-1]);

      if(!edgestate || u!=tail){
        newmin = MIN(newmin, INPUT_PARAM[u-1]);
        newmax = MAX(newmax, INPUT_PARAM[u-1]);
      }
    });

  if(!edgestate){
    newmin = MIN(newmin, INPUT_PARAM[tail-1]);
    newmax = MAX(newmax, INPUT_PARAM[tail-1]);
  }

  CHANGE_STAT[0] += (isinf(newmax) ? 0 : newmax-newmin) - (isinf(oldmax) ? 0 : oldmax-oldmin);
}


C_CHANGESTAT_FN(c_b1covrange) {
  double oldmin = R_PosInf, oldmax = R_NegInf, newmin = R_PosInf, newmax = R_NegInf;

  EXEC_THROUGH_OUTEDGES(tail, e, u, {
      oldmin = MIN(oldmin, INPUT_PARAM[u-1 - BIPARTITE]);
      oldmax = MAX(oldmax, INPUT_PARAM[u-1 - BIPARTITE]);

      if(!edgestate || u!=head){
        newmin = MIN(newmin, INPUT_PARAM[u-1 - BIPARTITE]);
        newmax = MAX(newmax, INPUT_PARAM[u-1 - BIPARTITE]);
      }
    });

  if(!edgestate){
    newmin = MIN(newmin, INPUT_PARAM[head-1 - BIPARTITE]);
    newmax = MAX(newmax, INPUT_PARAM[head-1 - BIPARTITE]);
  }

  CHANGE_STAT[0] += (isinf(newmax) ? 0 : newmax-newmin) - (isinf(oldmax) ? 0 : oldmax-oldmin);
}


I_CHANGESTAT_FN(i_nodefactordistinct) {
  unsigned int ncats = IINPUT_PARAM[0];
  
  ALLOC_STORAGE(ncats*N_NODES, unsigned int, freqs);

  EXEC_THROUGH_NET_EDGES(tail, head, e, {
      if(IINPUT_PARAM[head]) freqs[(tail-1)*ncats + IINPUT_PARAM[head]]++;
      if(IINPUT_PARAM[tail]) freqs[(head-1)*ncats + IINPUT_PARAM[tail]]++;
    });  
}

C_CHANGESTAT_FN(c_nodefactordistinct) {
  unsigned int ncats = IINPUT_PARAM[0];
  GET_STORAGE(unsigned int, freqs);
  int change = edgestate ? -1 : +1;
  if(IINPUT_PARAM[head]) {
    unsigned int oldfreq = freqs[(tail-1)*ncats + IINPUT_PARAM[head]];
    unsigned int newfreq = oldfreq + change;
    CHANGE_STAT[0] += (newfreq!=0) - (oldfreq!=0);
  }

  if(IINPUT_PARAM[tail]) {
    unsigned int oldfreq = freqs[(head-1)*ncats + IINPUT_PARAM[tail]];
    unsigned int newfreq = oldfreq + change;
    CHANGE_STAT[0] += (newfreq!=0) - (oldfreq!=0);
  }
}


U_CHANGESTAT_FN(u_nodefactordistinct) {
  unsigned int ncats = IINPUT_PARAM[0];
  GET_STORAGE(unsigned int, freqs);
  int change = edgestate ? -1 : +1;
  if(IINPUT_PARAM[head])
    freqs[(tail-1)*ncats + IINPUT_PARAM[head]] += change;
  if(IINPUT_PARAM[tail])
    freqs[(head-1)*ncats + IINPUT_PARAM[tail]] += change;
}




// ---------- TEST --------------

I_CHANGESTAT_FN(i_b1factordistinct) {
  unsigned int ncats = IINPUT_PARAM[0];
  
  ALLOC_STORAGE(ncats*(N_NODES - BIPARTITE), unsigned int, freqs);
  
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
  //if(IINPUT_PARAM[head]) freqs[(tail-1)*ncats + IINPUT_PARAM[head]]++;
    if(IINPUT_PARAM[tail]) freqs[(head-1-BIPARTITE)*ncats + IINPUT_PARAM[tail]]++;
  });  
}

C_CHANGESTAT_FN(c_b1factordistinct) {
  unsigned int ncats = IINPUT_PARAM[0];
  GET_STORAGE(unsigned int, freqs);
  int change = edgestate ? -1 : +1;
  
  if(IINPUT_PARAM[tail]) {
    unsigned int oldfreq = freqs[(head-1-BIPARTITE)*ncats + IINPUT_PARAM[tail]];
    unsigned int newfreq = oldfreq + change;
    CHANGE_STAT[0] += (newfreq!=0) - (oldfreq!=0);
  }
}


U_CHANGESTAT_FN(u_b1factordistinct) {
  unsigned int ncats = IINPUT_PARAM[0];
  GET_STORAGE(unsigned int, freqs);
  int change = edgestate ? -1 : +1;
  if(IINPUT_PARAM[tail])
    freqs[(head-1-BIPARTITE)*ncats + IINPUT_PARAM[tail]] += change;
}

