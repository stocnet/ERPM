// #include "ergm_changestat.h"
// //#include "changestats.h"
// #include "ergm_storage.h"
// #include "ergm_dyad_hashmap.h"
// #include "ergm_edgelist.h"


// // A macro indicating whether x is in [from,to)
// #define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))
// #define SQUARE(a,p,res); res=a; for(int x=1;x<p; x++){res=res*a;}

// /*****************
//  changestat: d_idegrange
// *****************/
// C_CHANGESTAT_FN(c_squared_sizes) {
//     int j, echange, res, p;
//     Vertex *id, *od;

//     id=IN_DEG;
//     od=OUT_DEG;

//     /* *** don't forget tail -> head */
//     echange=edgestate ? -1 : +1;
//     Vertex headdeg = od[head] + id[head]; 
//     // Vertex taildeg = od[tail] + id[tail];
//     // SQUARE(2, 2, res)
//     // Rprintf("val: %d, pw: %d, res: %d\n", 2, 2, res); 
    
//     for(j = 0; j < N_CHANGE_STATS; j++) {
//         Vertex from = INPUT_PARAM[3*j], to = INPUT_PARAM[3*j+1];
//         int power = INPUT_PARAM[3*j+2];

//         if (FROM_TO((headdeg + 1), from, to)){
//             // Rprintf("change stat : + %d - %d\n", pow(pow(headdeg + 1, power),power), pow(pow(headdeg, power),power));
//             SQUARE((headdeg + 1), power, res)
//             CHANGE_STAT[j] += res;
//             if (headdeg>=from){
//                 SQUARE((headdeg), power, res)
//                 CHANGE_STAT[j] -= res;
//             }
//         }
//         else if (headdeg + 1 == to)
//         {
//             SQUARE(headdeg, power, res)
//             CHANGE_STAT[j] -= res;
//         }
//     }
// }

// // S_CHANGESTAT_FN(s_groups) {
// //   /* int  echange;
// //      Vertex taild, headd=0, *id, *od; */
// //   Vertex *id, *od;

// //   id=IN_DEG;
// //   od=OUT_DEG;

// //   /* *** don't forget tail -> head */
// //   CHANGE_STAT[0] = 0.0;
// //   for(Vertex tail=1; tail <= N_NODES; tail++){
// //     if(od[tail] + id[tail] == 2)
// //       CHANGE_STAT[0] ++;
// //   }
// // }

// C_CHANGESTAT_FN(c_cliques) {

//   /* *** don't forget tail -> head */
//     Vertex node3;
//     Edge change = 0, e;
//     Rprintf("In cliques");
//     /* edgestate is 1 if edge exists and will disappear
//        edgestate is 0 if edge DNE and will appear */

//     // -3 * triangles
//     Vertex num_clq = INPUT_PARAM[0];

//     STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
//       change += IS_UNDIRECTED_EDGE(node3,tail);
//     }
//     STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
//       change += IS_UNDIRECTED_EDGE(node3,tail);
//     }
//     CHANGE_STAT[0] += change * (edgestate ? num_clq : num_clq);


//     // +1 * 2-stars

//     Vertex taild = OUT_DEG[tail] + IN_DEG[tail] - edgestate;
//     Vertex headd = OUT_DEG[head] + IN_DEG[head] - edgestate;
//     change = taild + headd;
//     CHANGE_STAT[0] += (edgestate ?  -change : change);
//     Rprintf("end cliques");

// }

// // Plus utilisé -> remplacer par ergm directement
// C_CHANGESTAT_FN(c_groups) {
//   int j, echange;
//   Vertex *id, *od;

//   id=IN_DEG;
//   od=OUT_DEG;

//   /* *** don't forget tail -> head */
//     echange=edgestate ? -1:+1;
//     Vertex taildeg = od[tail] + id[tail], headdeg = od[head] + id[head];
//     for(j = 0; j < N_CHANGE_STATS; j++) {
//         Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1], pow = INPUT_PARAM[2*j+2];

//         if(FROM_TO(taildeg+headdeg, from, to)){
//             Rprintf("change stat : +%d\n", sizeof(id)*pow);
//             CHANGE_STAT[j] += sizeof(id)*pow;
//         }
//     }

//     // STEP_THROUGH_OUTEDGES(head, e, node3){
//     //     Rprintf("head: %d, node3: %d, e: %d\n", head, node3, e);
//     // }
//     Rprintf("tail: %d, head: %d, taildeg: %d, headdeg: %d, size id: %d\n", tail, head, taildeg, headdeg, sizeof(id));
// }


// #undef FROM_TO
