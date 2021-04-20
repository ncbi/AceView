/*  file: intrinsictree.c
 *  Author: Ulrich Sauvage (ulrich@kaa.crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 20 14:05 1998 (fw)
 * Created: Thu Jun 23 12:24:12 1994 (ulrich)
 *-------------------------------------------------------------------
 */

/* $Id: intrinsictree.c,v 1.4 2020/05/30 16:50:34 mieg Exp $ */

/***************************************************************/

#include "acedb.h"
#include "keyset.h"
#include "topology.h"
#include "interval.h"

#ifndef NON_GRAPHIC
#include "graph.h"
#include "chrono.h"
#endif /* !NON_GRAPHIC */

/***************************************************************/

static INTTREE look ;

typedef struct { KEYSET next, zero ; int n, g ;} TR;
typedef struct { unsigned char tdis ; int i, j ;} TREEDIST ;

/***************************************************************/
/********************  Display  functions  *********************/
/***************************************************************/
#ifdef NON_GRAPHIC
void intrinsicTreeDestroy(void)
{
}
#else  /* !NON_GRAPHIC */
static Graph
  intTreeGraph = 0 ;
static Array
  intTreeCurrGraph = 0 ;

void intrinsicTreeDestroy(void)
{ int i ;
  Graph *gp ;

  if (intTreeCurrGraph)
    { i = arrayMax(intTreeCurrGraph) ;
      gp = arrp(intTreeCurrGraph, 0, Graph) - 1 ;
      while(gp++, i--)
	if (graphExists(*gp))
	  { graphActivate(*gp) ;
	    graphDestroy() ;
	  }
      arrayDestroy(intTreeCurrGraph) ;
    }
  if (intTreeGraph)
    if (graphExists(intTreeGraph))
      { graphActivate(intTreeGraph) ;
	graphDestroy() ;
      }
  intTreeGraph = 0 ;
} /* intrinsicTreeDestroy */

/***************************************************************/

static void localPreserve(void)
{ 
  if (!intTreeCurrGraph)
    intTreeCurrGraph = arrayCreate(16, Graph) ;
  array(intTreeCurrGraph, arrayMax(intTreeCurrGraph), Graph) = intTreeGraph ;
  intTreeGraph = 0 ;
} /* localPreserve */

/***************************************************************/

static MENUOPT intMapGraphMenu[] = {
  { graphDestroy,	"Quit" },
  { help,		"Help" },
  { graphPrint,		"Print" },
  { localPreserve,	"Preserve" },
  { 0, 0 }
};
#endif /* !NON_GRAPHIC */

/***************************************************************/
/************************* Sort Order **************************/
/***************************************************************/

static int treeOrder (const void *a, const void *b)
{ int x ;

  if (a == b)
    messcrash("Sort is stupid") ;
  x =  ((const TREEDIST *) b)->tdis > ((const TREEDIST *) a)->tdis ? 1 : 
     (((const TREEDIST *) b)->tdis == ((const TREEDIST *) a)->tdis ? 0 : -1) ;
  if (x == 0)
    x = ((const TREEDIST *) b)->i - ((const TREEDIST *) a)->i ;
  if (x == 0)
    x = ((const TREEDIST *) b)->j - ((const TREEDIST *) a)->j ;
  return
    x < 0 ? -1 : 1 ;
} /* treeOrder */

/***************************************************************/

static int intTreeOrder (const void *a, const void *b)
{ return 
    arr(look->tree,*(const KEY*)a,TR).n 
    - arr(look->tree,*(const KEY*)b,TR).n  ;
} /* intTreeOrder */

/***************************************************************/
/************************  Tree  Plot  *************************/
/***************************************************************/

static int intTreePlotDef(int i, int x, int fath, int prev)
{ 
  int j, ksmax, delta = 100, epsilon = 200 ;
  char *v , *vp, *v0 = 0 ;

  if (fath != -1)
    { v = v0 + (i > fath ? i + fath * look->nd : fath + i * look->nd) ;
      if (look->distances && assFind(look->distances, v, &vp))
	delta = vp - v0 ;
    }
  if (prev != -1)
    { v = v0 + (i > prev ? i + prev * look->nd : prev + i * look->nd ) ;
      if (look->distances && assFind(look->distances, v, &vp))
	epsilon = vp - v0 ;
    }
  if (arr(look->tree, i, TR).zero)
    { int z = x + 1 ;
      int k = i ;
      j = keySetMax(arr(look->tree, i, TR).zero) ;
      while(j--)
	k = intTreePlotDef(keySet(arr(look->tree, i, TR).zero, j), z++, i, k) ;
    }

#ifndef NON_GRAPHIC
  if (look->display)
    graphText
      ( messprintf
       ("%s     %d     %d",     
	name(keySet(look->def, i)), delta, epsilon
	), 2 + 8 * x, 2 + look->yTree) ;
#else
  if (0)  messprintf
       ("%s     %d     %d",     
	name(keySet(look->def, i)), delta, epsilon
	) ;
#endif /* !NON_GRAPHIC */

  array(look->cOrder, look->yTree++, int) = i ;
  ksmax = keySetMax((look->trd).ks) ;
  keySet((look->trd).ks, ksmax) = keySet(look->def, i) ;
  keySet((look->trd).x, ksmax)  = x ;
  return i ;
} /* intTreePlotDef */

/***************************************************************/

static void intTreePlotNext(int i, int x, int fath)
{ int n, j, k ;
  KEYSET s = arr(look->tree, i, TR).next ;

  arr(look->tree, i, TR).n = 1 ; /* to prevent looping */
  if (s)
    { n = 0 ;
      for(j = 0 ; j < keySetMax(s) ; j++)
	{ k = keySet(s, j) ;
	  if(!arr(look->tree, k, TR).n)
	    intTreePlotNext(k, x + n++, i) ;
	}
    }
  look->prev = intTreePlotDef(i, x, fath, look->prev) ;
} /* intTreePlotNext */

/***************************************************************/

static int intTreeLength(int k)
{ KEYSET s = arrp(look->tree, k, TR)->next ;
  int i, l = 0, n ;
  i = s ? keySetMax(s) : 0 ;
  arr(look->tree, k, TR).n = 1 ; /* To prevent looping */
  while(i--)
    { if(!arr(look->tree, keySet(s, i), TR).n) /* not the father */
	if ((n = intTreeLength(keySet(s, i))) > l)
	  l = n ;
    }
  arr(look->tree, k, TR).n = ++l ; /* ++ for self */
  return l ;
} /* intTreeLength */

/***************************************************************/

static int intTreeEnd(int k)
{ KEYSET s = arrp(look->tree, k, TR)->next ;
  int i, n ;
  int j = -1 , l = -1 , myself = arr(look->tree, k, TR).n ;

  i = s ? keySetMax(s) : 0 ;
  if (myself == 1)
    return k ;
  while(i--)
    { n = arr(look->tree, keySet(s,i), TR).n ;
      if (n < myself && n > l)
	{ l = n ;
	  j = keySet(s, i) ;
	}
    }
  if (j == -1)
    messcrash("Bad tree in intTreeEnd, sorry") ;
  return intTreeEnd(j) ;
}

/***************************************************************/

static void intTreeClearN(void)
{ int i ;
  TR *trp ;

  i = look->nd ;
  trp = arrp(look->tree, 0, TR) - 1 ;
  while(trp++, i--)
    trp->n = 0 ;
} /* intTreeClearN */

/***************************************************************/

static void intTreeSort(int k)
{ KEYSET s = arrp(look->tree,k,TR)->next ;
  int i = s ? keySetMax(s) : 0 ;

  arr(look->tree, k, TR).n *= -1 ; /* To prevent looping */
  while(i--)
    if(arr(look->tree, keySet(s, i), TR).n > 0) /* not the father */
      intTreeSort(keySet(s, i)) ;
	
  if (s) arraySort(s, intTreeOrder) ;
} /* intTreeSort */

/***************************************************************/

static int intTreeFindHead(int x)
{ int i ;

  intTreeLength(x) ;
  i = intTreeEnd(x) ;
  intTreeClearN() ;
  intTreeLength(i) ; /* Move backwards */
  intTreeSort(i) ;
  intTreeClearN() ;
  return i ;
} /* intTreeFindHead */

/***************************************************************/

static void intTreePlot(int nn)
{ int first = 0, j, k, segmax = 0 ;
  TR *destrp ;

#ifndef NON_GRAPHIC
  if (look->display)
    {
      if(graphActivate(intTreeGraph))
	{ graphPop() ;
	  graphClear() ;
	  if (look->map)
	    graphText (messprintf("Map : %s", name(look->map)),  1., .2) ;
	}
      else
	{ intTreeGraph = graphCreate(TEXT_FULL_SCROLL, "Intrinsic tree", 0.180,0.25,0.4,0.65) ;
	  graphMenu(intMapGraphMenu) ;
	  graphHelp("Deficiency_map") ;
	  if (look->map)
	    graphText (messprintf("Map : %s", name(look->map)),  1., .2) ;
	  graphTextBounds(80, look->nd + 5) ;
	}
    }
#endif /* !NON_GRAPHIC */
  look->segment = arrayCreate(nn, TREE_DEF) ;
  for (j = 1 ; j <= nn ; j++)
    { 
      for (k = 0; k < look->nd && arr(look->tree, k, TR).g != j; k++)
	;
      if (k == look->nd)
	continue ;
      first = intTreeFindHead(k) ;
      (look->trd).ks = keySetCreate() ;
      (look->trd).x  = keySetCreate() ;
      look->prev = - 1 ;
      intTreePlotNext(first, 0, -1) ;
      array(look->segment, segmax++, TREE_DEF) = look->trd ;
    }
  j = arrayMax(look->tree) ;
  destrp = arrp(look->tree, j - 1, TR) + 1 ;
  while(destrp--, j--)
    { keySetDestroy(destrp->zero) ;
      keySetDestroy(destrp->next) ;
    }
  arrayDestroy(look->tree) ;
  arrayDestroy(look->tabledis) ;
  assDestroy(look->distances) ;
#ifndef NON_GRAPHIC
  if (look->display)
    { graphTextBounds(80, look->yTree + 5) ;
      graphRedraw() ;
    }
#endif /* !NON_GRAPHIC */

  return;
} /* intTreePlot */

/***************************************************************/
/*********************** Make Distances ************************/
/***************************************************************/

static void makeShortDistances(int limit)
{ int i, j, jmax, ii, jj ;
  int u, v , dummy, n, ND ;
  char *assp, *assp0 = 0, *valcp ;
  KEY *ksp, *kspi, *kspj ;
  KEYSET ks = 0, *ksi, *ksj ;
  Associator dista ;
  Array table ;

  if (assExists(look->distances))
    messcrash ("Double allocation de look->distances") ;
  if (look->display)
    look->distances = assBigCreate(20 * look->nd) ;
  else
    look->distances = 0 ;
  look->tabledis = arrayCreate(10 * look->nd, TREEDIST) ;/* Tableau contenant les distances !infinie */
  dista = look->distances ;
  table = look->tabledis ;
  n = 0 ;/* indice dans le tableau tabledis */
  i = ND = look->nd ;
  ksi = arrp(look->mind, i - 1, KEYSET) + 1 ;
  while(ksi--, i--)
    { if (!(j = keySetMax (*ksi)) || keySet(*ksi, 0) & F_FLAG)
	continue ; /* Pour obtenir le keyset des voisins de i */
      ks = keySetReCreate(ks) ;
      jj = 0 ;
      kspi = arrp(*ksi, j - 1, KEY) + 1 ;
      while(kspi--, j--)
	{ ksj = arrp(look->dinm, *kspi, KEYSET) ;
	  if ((ii = keySetMax(*ksj)))
	    { kspj = arrp(*ksj, ii - 1, KEY) ;
	      while(ii--)
		{ if (*kspj < i) /* Pour ne pas avoir a faire d(i, i) et d(j, i) */
		    keySet(ks, jj++) = *kspj ;
		  kspj-- ;
		}
	    }
	  else messerror ("tableau vide dans makeShortDistances") ;
	}
      keySetSort(ks) ;
      keySetCompress(ks) ;
      jmax = keySetMax(ks) ;
      if (!jmax)
	continue ; /* i e pas de voisin en dessous => passe a la def suivante */
      ksp = arrp(ks, jmax - 1, KEY) + 1 ;
      while(ksp--, jmax-- && !(keySet(*ksi,0) & F_FLAG))
	{ ksj = arrp(look->mind, *ksp, KEYSET) ;
	  if (!keySetMax (*ksj) || keySet(*ksj, 0) & F_FLAG) /* test sur max normalement inutile */
	    continue ;
	  j = *ksp ;
	  valcp = assp0 + (i + ND * j) ;
	  ii = keySetMax(*ksi) ;
	  jj = keySetMax(*ksj) ;
	  kspi = arrp(*ksi, ii - 1, KEY) ;
	  kspj = arrp(*ksj, jj - 1, KEY) ;
	  u = ii + jj ;
	  v = 0 ;
	  while(ii && jj)
	    {
	      if (*kspi == *kspj)
		{ v++ ;
		  kspi-- ; ii-- ;
		  kspj-- ; jj-- ;
		}
	      else if (*kspi > *kspj)
		{ ii-- ;
		  kspi-- ;
		}
	      else
		{ jj-- ;
		  kspj-- ;
		}
	    }
	  u -= (2 * v) ;
	  if (!u && !v) u = 1 ;
	  dummy = 120 * u ;
	  dummy /= (u + v) ;
	  if (v == keySetMax(*ksj)) /* alors ksj inclus ou egale ksi */
	    dummy = 0 ;
	  else if (v == keySetMax(*ksi)) /* alors ksi inclus ksj */
	    { dummy = 0 ;
	      keySet(*ksi, 0) |= F_ZERO ;
	      if (dista)
		assInsert(dista, valcp, assp0) ; /* association avec dummy = 0 donc assp0 directement */
	      array(table, n, TREEDIST).tdis = (unsigned char)dummy ;
	      arr(table, n, TREEDIST).i = j ; /* car c'est i qu'on supprime au */
	      arr(table, n, TREEDIST).j = i ; /* retour dans le calcul de l'arbre */
	      dummy = 100 ; /* Pour etre au-dessus de la limite */
	      n++ ;
	    }
	  assp = assp0 + dummy ;
	  if (!dummy) /* ksj identique (ou inclus) a ksi */
	    keySet(*ksj, 0) |= F_ZERO ;
	  if (dummy <= limit)
	    { 
	      if (dista)
		assInsert(dista, valcp, assp) ;
	      array(table, n, TREEDIST).tdis = (unsigned char) dummy ;
	      arr(table, n, TREEDIST).i = i ;
	      arr(table, n, TREEDIST).j = j ;
	      n++ ;
	    }
	}
    }
  keySetDestroy(ks) ;
} /* makeShortDistances */

/***************************************************************/

static void intProporEdw(Array *intpro)
{ int i = look->nm, j, k = 0, l = 0, *ip ;
  KEYSET *ksp ;
  KEY *keyp ;

  *intpro = arrayCreate(i, int) ;
  ip = arrayp(*intpro, i - 1, int) ;
  ksp = arrp(look->dinm, i - 1, KEYSET) + 1 ;
  while (ksp--, i--)
    { j = keySetMax(*ksp) ;
      k = 0 ;
      keyp = arrp(*ksp, 0, KEY) - 1 ;
      while(keyp++, j--)
	(*keyp & F_PERE) ? k++ : l++ ;
      if (!k && !l)
	*ip-- = 50 ;
      else
	*ip-- = (100 * k) / (k + l) ;
    }
} /* intProporEdw */

/***************************************************************/

static void makeDistancesEdwards(int limit)
{ int i, j, jmax, ii, jj ;
  int u, v, w, dummy, n, ND ;
  char *assp, *assp0 = 0, *valcp ;
  KEY *ksp, *kspi, *kspj ;
  KEYSET ks = 0, *ksi, *ksj ;
  Associator dista ;
  Array table, propor = 0 ;

  look->distances = assBigCreate(20 * look->nd) ;
  look->tabledis = arrayCreate(10 * look->nd, TREEDIST) ;/* Tableau contenant les distances !infinie */
  dista = look->distances ;
  table = look->tabledis ;
  intProporEdw(&propor) ;
  n = 0 ;/* indice dans le tableau tabledis */
  i = ND = look->nd ;
  ksi = arrp(look->mind, i - 1, KEYSET) + 1 ;
  while(ksi--, i--)
    { if (keySet(*ksi, 0) & F_ZV_FLAG)
	continue ;
      j = keySetMax(*ksi) ; /* Pour obtenir le keyset des voisins de i */
      ks = keySetReCreate(ks) ;
      jj = 0 ;
      kspi = arrp(*ksi, j - 1, KEY) + 1 ;
      while(kspi--, j--)
	{ ksj = arrp(look->dinm, *kspi & F_WHO, KEYSET) ;
	  ii = keySetMax(*ksj) ;
	  kspj = arrp(*ksj, ii - 1, KEY) ;
	  while(ii--)
	    { if ((*kspi & *kspj & F_PEME) && ((*kspj & F_WHO) < i)) /* Pour ne pas avoir a faire d(i, i) et d(j, i) */
		keySet(ks, jj++) = *kspj & F_WHO ;
	      kspj-- ;
	    }
	}
      keySetSort(ks) ;
      keySetCompress(ks) ;
      jmax = keySetMax(ks) ;
      if (!jmax)
	continue ; /* i e pas de voisin en dessous => passe a la def suivante */
      ksp = arrp(ks, jmax - 1, KEY) + 1 ;
      while(ksp--, jmax-- && !(keySet(*ksi,0) & F_ZV_FLAG))
	{ ksj = arrp(look->mind, *ksp, KEYSET) ;
	  if (keySet(*ksj, 0) & F_ZV_FLAG)
	    continue ;
	  j = *ksp ;
	  valcp = assp0 + (i + ND * j) ;
	  ii = keySetMax(*ksi) ;
	  jj = keySetMax(*ksj) ;
	  kspi = arrp(*ksi, ii - 1, KEY) ;
	  kspj = arrp(*ksj, jj - 1, KEY) ;
	  u = v = 0 ;
	  while(ii && jj)
	    {
	      if ((*kspi & F_WHO) == (*kspj & F_WHO))
		{ 
		  if (*kspi & *kspj & F_PEME)
		    v += 200 ;
		  else
		    u += 200 ;
		  kspi-- ; ii-- ;
		  kspj-- ; jj-- ;
		}
	      else if ((*kspi & F_WHO) > (*kspj & F_WHO))
		{ ii-- ;
		  w = arr(propor, *kspi & F_WHO, int) ;
		  if (*kspi & F_PERE) /* inconue => 1 chance sur deux d'etre 1 ou autre */
		    { v += w ;
		      u += 100 - w ;
		    }
		  else
		    { v += 100 - w ;
		      u += w ;
		    }
		  kspi-- ;
		}
	      else
		{ jj-- ;
		  w = arr(propor, *kspj & F_WHO, int) ;
		  if (*kspj & F_PERE) /* inconue => 1 chance sur deux d'etre 1 ou autre */
		    { v += w ;
		      u += 100 - w ;
		    }
		  else
		    { v += 100 - w ;
		      u += w ;
		    }
		  kspj-- ;
		}
	    }
	  while(ii--)
	    { w = arr(propor, *kspi & F_WHO, int) ;
	      if (*kspi & F_PERE) /* inconue => 1 chance sur deux d'etre 1 ou autre */
		{ v += w ;
		  u += 100 - w ;
		}
	      else
		{ v += 100 - w ;
		  u += w ;
		}
	      kspi-- ;
	    }
	  while(jj--)
	    { w = arr(propor, *kspj & F_WHO, int) ;
	      if (*kspj & F_PERE) /* inconue => 1 chance sur deux d'etre 1 ou autre */
		{ v += w ;
		  u += 100 - w ;
		}
	      else
		{ v += 100 - w ;
		  u += w ;
		}
	      kspj-- ;
	    }
	  if (!u && !v) u = 1 ;
	  dummy = 120 * u ;
	  dummy /= (u + v) ;
	  if (!dummy) /* ksj identique (ou inclus) a ksi */
	    keySet(*ksj, 0) |= F_ZERO ;
	  if (dummy <= limit)
	    { assp = assp0 + dummy ;
	      assInsert(dista, valcp, assp) ;
	      array(table, n, TREEDIST).tdis = (unsigned char) dummy ;
	      arr(table, n, TREEDIST).i = i ;
	      arr(table, n, TREEDIST).j = j ;
	      n++ ;
	    }
	}
    }
  keySetDestroy(ks) ;
  arrayDestroy(propor) ;
}

/***************************************************************/

static void intCalPropor(Array *intpro)
{ int i = look->nd, j, k = 0, *ip ;
  KEYSET *ksp ;
  KEY *keyp ;

  *intpro = arrayCreate(i, int) ;
  ip = arrayp(*intpro, i - 1, int) ;
  ksp = arrp(look->mind, i - 1, KEYSET) + 1 ;
  while (ksp--, i--)
    { j = keySetMax(*ksp) ;
      k = 0 ;
      keyp = arrp(*ksp, 0, KEY) - 1 ;
      while(keyp++, j--)
	if (!(*keyp & F_NO)) 
	  k++ ;
      *ip-- = (100 * k) / look->nd ;
    }
} /* intCalPropor */

/***************************************************************/

static void makeDistances(int limit, int choix)
{ int i, j, k, jmax, ii, jj, dati, datj ;
  int u, v , dummy, n, ND ;
  char *assp, *assp0 = 0, *valcp ;
  KEY *ksp, *kspi, *kspj ;
  KEYSET ks = 0, *ksi, *ksj ;
  Associator dista ;
  Array table, propor = 0 ;

  look->distances = assBigCreate(20 * look->nd) ;
  look->tabledis = arrayCreate(10 * look->nd, TREEDIST) ;/* Tableau contenant les distances !infinie */
  dista = look->distances ;
  table = look->tabledis ;
  n = 0 ;/* indice dans le tableau tabledis */
  i = ND = look->nd ;
  if (choix == 4)
    intCalPropor(&propor) ;
  ksi = arrp(look->mind, i - 1, KEYSET) + 1 ;
  while(ksi--, i--)
    { if (keySet(*ksi, 0) & F_ZV_FLAG)
	continue ;
      j = keySetMax(*ksi) ; /* Pour obtenir le keyset des voisins de i */
      ks = keySetReCreate(ks) ;
      jj = 0 ;
      kspi = arrp(*ksi, j - 1, KEY) + 1 ;
      while(kspi--, j--)
	{ if (*kspi & F_NO) /* on ne prend que si il existe un vrai YES en commun */
	    continue ;
	  ksj = arrp(look->dinm, *kspi & F_WHO, KEYSET) ;
	  ii = keySetMax(*ksj) ;
	  kspj = arrp(*ksj, ii - 1, KEY) ;
	  while(ii--)
	    { if (!(*kspj & F_NO) && ((*kspj & F_WHO) < i)) /* Pour ne pas avoir a faire */
		keySet(ks, jj++) = *kspj & F_WHO ; /* d(i, i) et d(j, i) et ceux qui n'ont pas */
	      kspj-- ;                                      /* au moins un YES en commun */
	    }
	}
      keySetSort(ks) ;
      keySetCompress(ks) ;
      jmax = keySetMax(ks) ;
      if (!jmax)
	continue ; /* i e pas de voisin en dessous => passe a la def suivante */
      ksp = arrp(ks, jmax - 1, KEY) + 1 ;
      while(ksp--, jmax-- && !(keySet(*ksi,0) & F_ZERO))
	{ ksj = arrp(look->mind, *ksp, KEYSET) ;
	  if (keySet(*ksj, 0) & F_ZV_FLAG)
	    continue ;
	  j = *ksp ;
	  valcp = assp0 + (i + ND * j) ;
	  ii = keySetMax(*ksi) ;
	  jj = keySetMax(*ksj) ;
	  kspi = arrp(*ksi, ii - 1, KEY) ;
	  dati = (*kspi & F_NO) ? 0 : 1 ;
	  kspj = arrp(*ksj, jj - 1, KEY) ;
	  datj = (*kspj & F_NO) ? 0 : 1 ;
	  u = v = 0 ;
	  while(ii && jj)
	    {
	      if ((*kspi & F_NON_DIAG) || (*kspj & F_NON_DIAG))
		{
		  if (*kspi & F_NON_DIAG)
		    { kspi-- ; i-- ;
		      dati = (*kspi & F_NO) ? 0 : 1 ;
		    }
		  if (*kspj & F_NON_DIAG)
		    { kspj-- ; j-- ;
		      datj = (*kspj & F_NO) ? 0 : 1 ;
		    }
		  continue ;
		}
	      if ((*kspi & F_WHO) == (*kspj & F_WHO))
		{ switch(dati + datj)
		    { case 0: /*  no - no  */
			break ;
		      case 1: /* yes - no  */
			u += 200 ;
			break ;
		      case 2: /* yes - yes */
			v += 200 ;
			break ;
		      }
		  kspi-- ; ii-- ;
		  dati = (*kspi & F_NO) ? 0 : 1 ;
		  kspj-- ; jj-- ;
		  datj = (*kspj & F_NO) ? 0 : 1 ;
		}
	      else if ((*kspi & F_WHO) > (*kspj & F_WHO))
		{ switch(dati)
		    { case 0: /*  no - inc */
			if (choix == 3)
			  u += 100 ;
			if (choix == 4)
			  u += arr(propor, i, int) ;
			break ;
		      case 1: /* yes - inc */
			switch(choix)
			  { case 1:
			      u += 200 ;
			      break ;
			    case 2:
			      u += 100 ;
			      break ;
			    case 3:
			      u += 50 ;
			      v += 50 ;
			      break ;
			    case 4:
			      v += arr(propor, i, int) ;
			      u += 100 - arr(propor, i, int) ;
			      break ;
			    default:
			      break ;
			    }
			break ;
		      }
		  ii-- ;
		  kspi-- ;
		  dati = (*kspi & F_NO) ? 0 : 1 ;
		}
	      else
		{ switch(datj)
		    { case 0: /*  no - inc */
			if (choix == 3)
			  u += 100 ;
			if (choix == 4)
			  u += arr(propor, j, int) ;
			break ;
		      case 1: /* yes - inc */
			switch(choix)
			  { case 1:
			      u += 200 ;
			      break ;
			    case 2:
			      u += 100 ;
			      break ;
			    case 3:
			      u += 50 ;
			      v += 50 ;
			      break ;
			    case 4:
			      v += arr(propor, j, int) ;
			      u += 100 - arr(propor, j, int) ;
			      break ;
			    default:
			      break ;
			    }
			break ;
		      }
		  jj-- ;
		  kspj-- ;
		  datj = (*kspj & F_NO) ? 0 : 1 ;
		}
	    }
	  k = i ;
	  if (jj)
	    { if (ii)
		messcrash("Probleme dans makeDistances") ;
	      ii = jj ;
	      dati = datj ;
	      kspi = kspj ;
	      k = j ;
	    }
	  while(ii--)
	    { switch(dati)
		{ case 0: /*  no - inc */
		    if (choix == 3)
		      u += 100 ;
		    if (choix == 4)
		      u += arr(propor, k, int) ;
		    break ;
		  case 1: /* yes - inc */
		    switch(choix)
		      { case 1:
			  u += 200 ;
			  break ;
			case 2:
			  u += 100 ;
			  break ;
			case 3:
			  u += 50 ;
			  v += 50 ;
			  break ;
			case 4:
			  v += arr(propor, k, int) ;
			  u += 100 - arr(propor, k, int) ;
			  break ;
			default:
			  break ;
			}
		    break ;
		  }
	      kspi-- ;
	      dati = (*kspi & F_NO) ? 0 : 1 ;
	    }
	  if (!u && !v) u = 1 ;
	  dummy = 120 * u ;
	  dummy /= (u + v) ;
	  if (v == 200 * keySetMax(*ksj)) /* alors ksj inclus ou egale ksi */
	    dummy = 0 ;
	  else if (v == 200 * keySetMax(*ksi)) /* alors ksi inclus ksj */
	    { dummy = 0 ;
	      keySet(*ksi, 0) |= F_ZERO ;
	      assInsert(dista, valcp, assp0) ; /* association a 0 */
	      array(table, n, TREEDIST).tdis = (unsigned char)dummy ;
	      arr(table, n, TREEDIST).i = j ; /* car c'est i qu'on supprime au */
	      arr(table, n, TREEDIST).j = i ; /* retour dans le calcul de l'arbre */
	      dummy = 100 ; /* Pour etre au-dessus de la limite */
	      n++ ;
	    }
	  assp = assp0 + dummy ;
	  if (!dummy) /* ksj identique (ou inclus) a ksi */
	    keySet(*ksj, 0) |= F_ZERO ;
	  if (dummy <= limit)
	    { assInsert(dista, valcp, assp) ;
	      array(table, n, TREEDIST).tdis = (unsigned char) dummy ;
	      arr(table, n, TREEDIST).i = i ;
	      arr(table, n, TREEDIST).j = j ;
	      n++ ;
	    }
	}
    }
  keySetDestroy(ks) ;
  if (choix == 4)
    arrayDestroy(propor) ;
} /* makeDistances */

/***************************************************************/
/************************  Make  Tree  *************************/
/***************************************************************/

static int makeIntTree(void)
{ Array tree ;
  int nmax, narb, inarb, newi, newj, valg ;
  int valni, valnj , izer ;
  KEYSET newKs ;

  look->tree = arrayCreate(look->nd, TR) ;
  tree = look->tree ;
  array(tree, look->nd - 1, TR).g = 0 ;
  nmax = arrayMax(look->tabledis) ;
  arraySort(look->tabledis, treeOrder) ;
  narb = 0 ;  /* number of small tree */
  inarb = 0 ;   /* indice de l'arbre */
  while(nmax--)
    { newi = arr(look->tabledis, nmax, TREEDIST).i ;
      newj = arr(look->tabledis, nmax, TREEDIST).j ;
      if (!(arr(look->tabledis, nmax, TREEDIST).tdis)) /* Suppression des identiques */
	{ arr(tree, newj, TR).g = - 1 ;
	  newKs = arr(tree, newi, TR).zero ;
	  if (!arrayExists(newKs))
	    { newKs = keySetCreate() ;
	      arr(tree, newi, TR).zero = newKs ;
	    }
	  izer = keySetMax(newKs) ;
	  if (!(arr(tree, newi, TR).g))
	    { arr(tree, newi, TR).g = ++inarb ;
	      narb++ ;
	    }
	  keySet(newKs, izer) = newj ;
	  continue ;
	}
      if (keySet(arr(look->mind, newi, KEYSET), 0) & F_ZERO)
	continue ;                                 /* identiques dans un autre arbre */
      if (keySet(arr(look->mind, newj, KEYSET), 0) & F_ZERO)
	continue ;
      if (!(arr(tree, newi, TR).g))
	{ newKs = keySetCreate() ;
	  arr(tree, newi, TR).next = newKs ;
	  keySet(arr(tree, newi, TR).next, 0) = newj ;
	  if (!(arr(tree, newj, TR).g))
	    { newKs = keySetCreate() ;
	      arr(tree, newj, TR).next = newKs ;
	      keySet(arr(tree, newj, TR).next, 0) = newi ;
	      arr(tree, newi, TR).g = ++inarb ;
	      arr(tree, newj, TR).g = inarb ;
	      narb++ ;
	    }
	  else
	    { 
	      if (!(arr(tree, newj, TR).next)) /* il est possible d'etre dans un arbre */
		{ valnj = 0 ;                  /* sans avoir de next, si on a un zero */
		  newKs = keySetCreate() ;
		  arr(tree, newj, TR).next = newKs ;
		}
	      else
		valnj = keySetMax(arr(tree, newj, TR).next) ;
	      keySet(arr(tree, newj, TR).next, valnj) = newi ;
	      arr(tree, newi, TR).g = arr(tree, newj, TR).g ;
	    } 
	}
      else
	{
	  if (arr(tree, newj, TR).g)
	    {
	      if (arr(tree, newi, TR).g != arr(tree, newj, TR).g)
		{ narb-- ;
		  if (!(arr(tree, newi, TR).next))
		    { valni = 0 ;
		      newKs = keySetCreate() ;
		      arr(tree, newi, TR).next = newKs ;
		    }
		  else
		    valni = keySetMax(arr(tree, newi, TR).next) ;
		  if (!(arr(tree, newj, TR).next))
		    { valnj = 0 ;
		      newKs = keySetCreate() ;
		      arr(tree, newj, TR).next = newKs ;
		    }
		  else
		    valnj = keySetMax(arr(tree,newj,TR).next) ;
		  keySet(arr(tree, newi, TR).next, valni) = newj ;
		  keySet(arr(tree, newj, TR).next, valnj) = newi ;
		  if (arr( tree, newi, TR).g < arr(tree, newj, TR).g)
		    { valnj = arr(tree, newj, TR).g ;
		      valni = look->nd ;
		      valg = arr(tree, newi, TR).g ;
		      while (valni--)
			{ if ((arr(tree, valni, TR).g) == valnj)
			    arr(tree, valni, TR).g = valg ;
			}
		    }
		  else
		    { valni = arr(tree, newi, TR).g ;
		      valnj = look->nd ;
		      valg = arr(tree, newj, TR).g ;
		      while (valnj--)
			{ if ((arr(tree, valnj, TR).g) == valni)
			    arr(tree, valnj, TR).g = valg ;
			}
		    }
		}
	    }
	  else
	    {
	      if (!(arr(tree, newi, TR).next))
		{ valni = 0 ;
		  newKs = keySetCreate() ;
		  arr(tree, newi, TR).next = newKs ;
		}
	      else
		valni = keySetMax(arr(tree, newi, TR).next) ;
	      keySet(arr(tree, newi, TR).next, valni) = newj ;
	      newKs = keySetCreate() ;
	      arr(tree, newj, TR).next = newKs ;
	      keySet(arr(tree, newj, TR).next, 0) = newi ;
	      arr(tree, newj, TR).g = arr(tree, newi, TR).g ;
	    }
	}
    }
  newi = look->nd ;
  while (newi--)
    { if (!arr(tree, newi, TR).g)
	{ narb++ ;
	  arr(tree, newi, TR).g = ++inarb ;
	}
    }
  return inarb ;
} /* makeIntTree */

/***************************************************************/

static Array intRetourneSeg(void)
{ int i = arrayMax(look->segment), j, k ;
  TREE_DEF *ksmp, *ksmq ;
  KEY *keys, *keyx ;
  Array maille = arrayCreate(i, TREE_DEF) ;

  if (!i) return maille ;
  ksmp = arrp(look->segment, i - 1, TREE_DEF) + 1 ;
  ksmq = arrayp(maille, i - 1, TREE_DEF) + 1 ;
  while(ksmp--, ksmq--, i--)
    { j = keySetMax(ksmp->ks) ;
      keys = arrp(ksmp->ks, j - 1, KEY) ;
      keyx = arrp(ksmp->x, j - 1, KEY) ;
      ksmq->ks = keySetCreate() ;
      ksmq->x = keySetCreate() ;
      k = 0 ;
      while(j--)
	{ keySet(ksmq->ks, k) = *keys-- ;
	  keySet(ksmq->x, k++) = *keyx-- ;
	}
      keySetDestroy(ksmp->ks) ;
      keySetDestroy(ksmp->x) ;
    }
  return maille ;
} /* intRetourneSeg */

/***************************************************************/

static void intDestroySeg()
{ int i ;
  TREE_DEF *ksm ;

  i = arrayMax(look->segment) ;
  ksm = arrp(look->segment, 0, TREE_DEF) - 1 ;
  while(ksm++, i--)
    { keySetDestroy(ksm->ks) ;
      keySetDestroy(ksm->x) ;
    }
  arrayDestroy(look->segment) ;
} /* intDestroySeg */

/***************************************************************/

void intCptSupFlag(Array donnee, KEY flag)
{ int i, j ;
  KEYSET *ksp ;
  KEY *mkp ;

  i = arrayMax(donnee) ;
  ksp = arrp(donnee, i - 1, KEYSET) + 1 ;
  while(ksp--, i--)
    { if ((j = keySetMax(*ksp)))
	{ mkp = arrp(*ksp, j - 1, KEY) ;
	  while(j--)
	    *mkp-- &= ~flag ;
	}
    }
} /* intCptSupFlag */

/***************************************************************/

int intCptTree(DEFCPT defCptLook)
{ int nbTree, dlimit, mchoix , n ;

  dlimit = defCptLook->dlimit ;
  if (dlimit <= 0 || dlimit > 100)
    messcrash("Bad value for distance limit in intCptTree") ;
  dlimit *= 12 ; dlimit /= 10 ;
  if (dlimit == 120) dlimit = 119 ;
  if ((mchoix = defCptLook->method) > 4) /* attention en fait M_EDWARD */
    messcrash("Bad value for method in intCptTree") ;
  look = (INTTREE)messalloc(sizeof(struct INTTREESTUFF)) ;
  look->magic = INTTREEMAG ;
  look->yTree = 0 ;
  look->display = defCptLook->display ;
  look->map = defCptLook->selectedMap ;
  look->def = defCptLook->def ;
  look->genes = defCptLook->mar ;
  look->cOrder = defCptLook->colOrder ;
  look->mind = defCptLook->marInDef ;
  look->dinm = defCptLook->defInMar ;
  if ((look->nd = arrayMax(look->mind)) != keySetMax(look->def) || 
      (look->nm = arrayMax(look->dinm)) != keySetMax(look->genes))
    messcrash("uncompatible dimension of arrays data") ;
  if (arrayMax(look->cOrder) != look->nd)
    messcrash("uncompatible dimension of array Order") ;
  if (!look->nd || !look->nm)
    messcrash("ND or NM is null") ;
  if (!(defCptLook->whatDis & F_G_ZERO))
    { intCptSupFlag(look->mind, F_ZERO) ;
      intCptSupFlag(look->dinm, F_ZERO) ;
    }
  if (!(defCptLook->whatDis & F_E_YN) && mchoix != M_EDWARD)   /* attention choix des distances a revoir */
    makeShortDistances(dlimit) ;
  else if (mchoix == M_EDWARD)
    makeDistancesEdwards(dlimit) ;
  else
    makeDistances(dlimit, defCptLook->whatDis & 7) ;
  if (!arrayMax(look->tabledis))
    { arrayDestroy(look->tabledis) ;
      look->magic = 0 ;
      messfree(look) ;
      return 0 ;
    }
  nbTree = makeIntTree() ;
  intTreePlot(nbTree) ;  /* destroys look->tree and look->tabledis */
  n = arrayMax(look->segment) ;
  if (mchoix == M_EDWARD)
    { intDestroySeg() ;
      intCptSupFlag(look->mind, F_ZERO) ;
      intCptSupFlag(look->dinm, F_ZERO) ;
    }
  else if(n)
    defCptLook->maillon = intRetourneSeg() ;

  arrayDestroy(look->segment) ; 
  look->magic = 0 ;
  messfree(look) ;
  return n ;
} /* intCptTree */

/***************************************************************/
/************************* eof *********************************/
/***************************************************************/
