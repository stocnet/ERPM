// ============================================================================
// Fichier : changestat_cov_fullmatch.c
// Terme   : cov_fullmatch
// Stat    : T = sum_g 1[n_g in S] * 1[exists r : n_{g,r} = n_g]
//           Variante ciblée: 1[n_{g,kappa} = n_g]
// ============================================================================

#include "ergm_changestat.h"
#include "ergm_storage.h"      // pour Calloc/Free
#include <R_ext/Print.h>
#include <string.h>            // memset

#define DEBUG_COV_FULLMATCH 0  /* mettre à 0 pour couper le debug */
#define UNUSED_WARNING(x) (void)x

/* ----------------------------------------------------------------------------
 * Vérifie si une taille n est incluse dans l'ensemble des tailles autorisées S.
 * Si L==0, toutes les tailles sont acceptées.
 * -------------------------------------------------------------------------- */
static inline int in_sizes(int n, int L, const double *sizes){
  if(L==0) return 1;
  for(int i=0;i<L;i++) if((int)sizes[i]==n) return 1;
  return 0;
}

/* ----------------------------------------------------------------------------
 * Fonction locale group_flag()
 * Évalue la contribution d’un groupe g :
 *   - calcule n_g = nombre d’acteurs reliés à g (robuste: IN+OUT, dédoublonné) ;
 *   - calcule les effectifs par catégorie ;
 *   - retourne 1 si le groupe est « full match » (entièrement homogène).
 *
 * cats[1..n1] : index de modalité 0..K (0 = NA) pour chaque acteur.
 * Si target>0, teste la modalité cible (n_{g,target} == n_g).
 * Sinon, teste s’il existe une catégorie r telle que n_{g,r} == n_g.
 * -------------------------------------------------------------------------- */
static double group_flag(Vertex g, int n1, int L, const double *sizes,
                         int K, int target, const double *cats,
                         Network *nwp){

  int ng = 0, has_na = 0;
  int cnt_max = 0;

  /* Allocation sécurisée : pile si petit K */
  int use_stack_cnt = (K > 0 && K <= 1024);
  int cntK[1024];
  int *cnt = NULL;
  if(K > 0){
    cnt = use_stack_cnt ? cntK : (int*)Calloc(K, int);
    for(int i=0;i<K;i++) cnt[i]=0;
  }

  /* ----------------------------------------------------------------------
   * Marquage des acteurs déjà comptés pour éviter les doubles comptes
   * (utile si réseau non dirigé ou si arêtes réciproques existent).
   * On alloue n1 octets et on remettra à 0 uniquement pour les voisins vus.
   * -------------------------------------------------------------------- */
  unsigned char *seen = (unsigned char*)Calloc(n1, unsigned char); /* Calloc => 0 */
  Vertex h; Edge e;

  /* OUT-neighbors g -> h */
  STEP_THROUGH_OUTEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        ng++;
        int cat = (int)cats[idx];
        if(cat==0){
          has_na = 1;
        }else if(cat>=1 && cat<=K){
          int v = ++cnt[cat-1];
          if(v>cnt_max) cnt_max = v;
        }
      }
    }
  }

  /* IN-neighbors h -> g */
  STEP_THROUGH_INEDGES(g, e, h){
    if(h <= (Vertex)n1){
      int idx = (int)h - 1;
      if(!seen[idx]){
        seen[idx] = 1;
        ng++;
        int cat = (int)cats[idx];
        if(cat==0){
          has_na = 1;
        }else if(cat>=1 && cat<=K){
          int v = ++cnt[cat-1];
          if(v>cnt_max) cnt_max = v;
        }
      }
    }
  }

  /* (Debug) degrés et orientation */
  #if DEBUG_COV_FULLMATCH
    Rprintf("[cov_fullmatch][deg] g=%d OUT_DEG=%d IN_DEG=%d DIRECTED=%d\n",
            (int)g, (int)OUT_DEG[g], (int)IN_DEG[g], (int)DIRECTED);
  #endif

  /* Groupe vide : jamais « full match » */
  if(ng == 0){
    #if DEBUG_COV_FULLMATCH
      Rprintf("[cov_fullmatch][group_flag] g=%d ng=0 -> res=0 (groupe vide)\n", (int)g);
    #endif
    if(cnt && !use_stack_cnt) Free(cnt);
    Free(seen);
    return 0.0;
  }

  /* Vérification du filtre de tailles S */
  if(!in_sizes(ng, L, sizes)){
    #if DEBUG_COV_FULLMATCH
      Rprintf("[cov_fullmatch][group_flag] g=%d ng=%d -> hors S\n", (int)g, ng);
    #endif
    if(cnt && !use_stack_cnt) Free(cnt);
    Free(seen);
    return 0.0;
  }

  /* ----------------------------------------------------------------------
   * Calcul du résultat :
   *  - Variante ciblée : tous les acteurs du groupe partagent la modalité cible.
   *  - Variante non ciblée : il existe une modalité couvrant tout le groupe.
   * -------------------------------------------------------------------- */
  double res;
  if(target>0){
    int c = has_na ? -1 : (K>0 ? cnt[target-1] : 0);
    res = (!has_na && c==ng) ? 1.0 : 0.0;
    #if DEBUG_COV_FULLMATCH
      Rprintf("[cov_fullmatch][group_flag] g=%d ng=%d target=%d c=%d has_na=%d -> res=%.0f\n",
              (int)g, ng, target, c, has_na, res);
    #endif
  }else{
    res = (!has_na && cnt_max==ng) ? 1.0 : 0.0;
    #if DEBUG_COV_FULLMATCH
      Rprintf("[cov_fullmatch][group_flag] g=%d ng=%d cnt_max=%d has_na=%d -> res=%.0f\n",
              (int)g, ng, cnt_max, has_na, res);
    #endif
  }

  if(cnt && !use_stack_cnt) Free(cnt);

  /* Remise à zéro des marques pour ce groupe, via les mêmes parcours. */
  STEP_THROUGH_OUTEDGES(g, e, h){ if(h <= (Vertex)n1) seen[(int)h-1] = 0; }
  STEP_THROUGH_INEDGES(g, e, h){ if(h <= (Vertex)n1) seen[(int)h-1] = 0; }
  Free(seen);

  return res;
}

/* ----------------------------------------------------------------------------
 * Change statistic un-toggle pour cov_fullmatch
 * -------------------------------------------------------------------------- */
C_CHANGESTAT_FN(c_cov_fullmatch){
  ZERO_ALL_CHANGESTATS(0);
  UNUSED_WARNING(edgestate);

  // INPUT_PARAM layout:
  // [0]=n1,
  // [1]=L,
  // [2..(1+L)]=sizes[L],
  // [2+L]=K,
  // [3+L]=target,
  // [4+L..]=cats[n1]  (cats est 1..n1, valeurs entières 0..K)

  
  const double *ip    = INPUT_PARAM;
  const int n1        = (int)ip[0];
  const int L         = (int)ip[1];
  const double *sizes = ip + 2;
  const int K         = (int)ip[2+L];
  const int target    = (int)ip[3+L];
  const double *cats  = ip + 4 + L;

  #if DEBUG_COV_FULLMATCH
    Rprintf("[cov_fullmatch] n1=%d L=%d K=%d target=%d\n", n1, L, K, target);
  #endif

  /* ----------------------------------------------------------------------
   * Identifier les sommets acteur et groupe du toggle courant
   * -------------------------------------------------------------------- */
  Vertex a = tail, b = head;
  Vertex actor = (a <= (Vertex)n1) ? a : b;
  Vertex group = (a <= (Vertex)n1) ? b : a;

  UNUSED_WARNING(actor);

  /* Diagnostic : état de l’arête avant le toggle (robuste: OUT || IN) */
  #if DEBUG_COV_FULLMATCH
    int present_out = IS_OUTEDGE(actor, group);
    int present_in  = IS_INEDGE(group, actor);
    int actor_cat   = (int)cats[(int)actor - 1];
    Rprintf("[cov_fullmatch] actor=%d cat=%d group=%d present_out=%d present_in=%d\n",
            (int)actor, actor_cat, (int)group, present_out, present_in);
  #endif

  /* ----------------------------------------------------------------------
   * Évaluer la contribution du groupe avant/après un toggle virtuel
   * -------------------------------------------------------------------- */
  double F_before = group_flag(group, n1, L, sizes, K, target, cats, nwp);

  // Appliquer le toggle virtuellement (mono-toggle API)
  TOGGLE(a, b);

  double F_after  = group_flag(group, n1, L, sizes, K, target, cats, nwp);

  // Annuler le toggle
  TOGGLE(a, b);

  /* ----------------------------------------------------------------------
   * Mise à jour du change statistic
   * -------------------------------------------------------------------- */
  CHANGE_STAT[0] += (F_after - F_before);

  #if DEBUG_COV_FULLMATCH
    Rprintf("[cov_fullmatch] a=%d b=%d group=%d before=%.0f after=%.0f delta=%.0f\n",
            (int)a, (int)b, (int)group, F_before, F_after, (F_after-F_before));
  #endif
}
