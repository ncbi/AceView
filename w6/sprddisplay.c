/*  File: sprddisplay.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
    To get an automatic display of a multimap
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 17:47 1998 (fw)
 * Created: Mon Jun 21 00:00:09 1993 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: sprddisplay.c,v 1.2 2007/03/16 18:27:17 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "spreaddisp_.h"

#include "a.h"
#include "bs.h"
#include "tags.h"
#include "sysclass.h"
#include "classes.h"
#include "systags.h"
#include "query.h"
#include "pick.h"
#include "lex.h"


static BOOL doMultiMapDisplay2 (char *title, KEYSET ks, int cl, 
				KEY anchor_group, KEY anchor, int minMap)
{ OBJ obj ;
  Stack s = 0 ;
  KEYSET loci = 0 , loci1 = 0 ; 
  Array aa = 0, flipA ; int i, j, n = 0 ;

  aa = arrayCreate(100, BSunit) ;
  loci = keySetCreate() ; n = 0 ;
  flipA = arrayCreate(12, int) ;

  i = keySetMax(ks) ;
  while(i--)
    if ((obj = bsCreate(keySet(ks, i))))
      { arrayMax(aa) = 0 ;
	if (bsFindTag(obj, _Contains) &&
	    bsFlatten(obj, 2, aa))
	  for (j = 1 ; j < arrayMax(aa); j += 2)
	    keySet(loci, n++) = arr(aa, j, BSunit).k ;
	if (bsFindTag(obj, _Flipped))
	  array(flipA,i,int) = 1 ;
	bsDestroy(obj) ;
      }
  keySetSort(loci) ;
  keySetCompress(loci) ;

  s = stackCreate(100) ;
  if (!anchor)
    { loci1 = loci ;
      loci = query(loci1, "COUNT Map >= 1") ;  
      
      if (ks)
	{
	  pushText(s,"Title ") ;
	  catText(s, "Multi-map ") ;
	  if (title)
	      catText (s, title) ;
	  catText (s, "\n\n") ;
	  catText(s,"Colonne 1\n  Width 12 \n Visible \n Class Locus\n") ;
	  catText(s,messprintf("Condition COUNT Map >= %d\n\n", minMap)) ;
	  
	  for (i = 0; i <  keySetMax(ks) ; i++)
	    { catText(s, "Colonne ") ;
	      catText(s, messprintf("%d\n", i + 2)) ;
	      catText(s, messprintf("Subtitle %s\n", name(keySet(ks,i)))) ;
	      catText(s, "Float\n From 1 \n Tag Map = \"") ;
	      catText(s, name(keySet(ks,i))) ;
	      catText(s, "\"") ;
	      catText(s," # Position") ;
	      catText(s," \nVisible \nOptional") ;
	      if (array(flipA,i,int))
		catText(s,"\nFlipMap") ;
	      catText(s, "\n\n") ;
	    }
	  
	  spreadDispCreate(FALSE) ;
	  { 
	    SPREAD spread = currentSpread("multiMapDisplay") ;
	    
	    spreadDoReadDefinitions (spread, 0, 0, s, "", FALSE) ; /* will close f */
	    spread->quitWithoutConfirmation = TRUE ;
	    
	    /* try to restore data from disk */
	    spreadRecomputeKeySet (spread, loci, 0, 0, 0) ;
	    spreadDisplayData(spread) ;
	    
	    /* try to remember data on disk */
	    spreadMapCreate() ;
	  }
	}
    }
  else
    { loci1 = loci ;
      loci = query(loci1, messprintf(">%s", name(anchor_group))) ;
      
      if (ks)
	{
	  pushText(s,"Title ") ;
	  catText(s, "\"Homology-map ") ;
	  if (title)
	    catText (s, title) ;
	  catText (s, "\"\n\n") ;
	  catText(s,"Colonne 1\n  Width 12 \n Visible \n Class ") ;
	  catText (s, pickClass2Word (cl)) ;
	  catText (s, "\n") ;
	  
	  for (i = 0; i <  keySetMax(ks) ; i++)
	    { catText(s, "Colonne ") ;
	      catText(s, messprintf("%d\n", 2*i + 2)) ;
	      catText(s, messprintf("Subtitle Locus on %s\n", name(keySet(ks,i)))) ;
	      catText(s,"Class ") ;
	      catText (s, pickClass2Word (cl)) ;
	      catText (s, "\n From 1 \n Tag ") ;
	      catText (s, name (anchor)) ;
	      catText (s, "\n ") ;
	      catText(s,"Condition Map = \"") ;
	      catText(s,name(keySet(ks,i))) ;
	      catText(s,"\"\n") ;
	      catText(s," \nVisible \nOptional") ;
	      catText(s, "\n\n") ;

	      catText(s, "Colonne ") ;
	      catText(s, messprintf("%d\n", 2*i + 3)) ;
	      catText(s, messprintf("Subtitle %s\n", name(keySet(ks,i)))) ;
	      catText(s,messprintf("Float\n From %d \n Tag Map = \"", 2*i + 2)) ;
	      catText(s,name(keySet(ks,i))) ;
	      catText(s,"\" # Position") ;
	      catText(s," \nVisible \nOptional") ;
	      if (array(flipA,i,int))
		catText(s,"\nFlipMap") ;
	      catText(s, "\n\n") ;
	    }
	  
	  spreadDispCreate(FALSE) ;
	  { 
	    SPREAD spread = currentSpread("multiMapDisplay") ;
	    
	    spreadDoReadDefinitions (spread, 0, 0, s, "", FALSE) ; /* will close f */
	    spread->quitWithoutConfirmation = TRUE ;
	    
	    /* try to restore data from disk */
	    spreadRecomputeKeySet (spread, loci, 0, 0, 0) ;
	    spreadDisplayData(spread) ;
	    
	    /* try to remember data on disk */
	    spreadMapCreate() ;
	  }
	}
    }

/*   do NOT destroy ks, it belongs to keysetdisplay */
  keySetDestroy(loci) ;
  if (keySetExists(loci1))
    keySetDestroy(loci1) ;
  arrayDestroy(aa) ;
  arrayDestroy(flipA) ;
  stackDestroy(s) ;
  return TRUE ;
}

/********************************************************************************/

BOOL multiMapDisplay (KEY key, KEY from, BOOL reuse) 
{ OBJ obj ;
  KEY map, anchor = 0 , anchor_group;
  int n = 0 , cl, minMap = 2 ;
  KEYSET ks = keySetCreate() ;
  char *cp ;

  if (!key)
    return FALSE ;

  if (class(key) != _VMultiMap ||
      !(obj = bsCreate(key)))
    goto abort ;

  bsGetData(obj, _Min, _Int, &minMap) ;
  if (minMap <1) minMap = 1 ;
  if (bsGetKey(obj, _Map, &map))
    do { 	  
	 keySet(ks, n++) = map ;
       } while (bsGetKey(obj, _bsDown, &map)) ;
  if (minMap >n)
    minMap = n ;
  if (bsGetData (obj, _Anchor, _Text, &cp) &&
      (cl = pickWord2Class (cp)) &&
      bsGetData (obj, _bsRight, _Text, &cp) &&
      lexword2key (cp, &anchor_group, _VSystem)  &&
      bsGetData (obj, _bsRight, _Text, &cp) &&
      lexword2key (cp, &anchor, _VSystem)) ;
  else
    { cl = 0 ; anchor = 0 ; anchor_group = 0 ; }
  bsDestroy(obj) ;
  return 
    doMultiMapDisplay2 (name(key), ks, cl, anchor_group, anchor, minMap) ;
 abort:  
  keySetDestroy(ks) ;
  display (key, 0, TREE) ;
  return FALSE ;
}

/********************************************************************************/

void spreadTableDisplay (const char *tableName, const char *parms)
{ KEY key = 0 ;
  FILE *f = 0 ;

  if (!lexword2key (tableName, &key, _VTable) &&
      ! (f = filopen(tableName,"","r")))
    return ;
  
  spreadDispCreate(FALSE) ;
  { 
    SPREAD spread = currentSpread("spreadTableDisplay") ;
    
    spreadDoReadDefinitions (spread, key, f, 0, parms, FALSE) ; /* will close f */
    spreadDoRecompute (spread) ;
    spreadDisplayData (spread) ;
  }
}

BOOL doMultiMapDisplay (KEYSET ks)
{ return doMultiMapDisplay2 (0, ks, 0, 0, 0, 2) ;
}

/***********************************************************************/
/************************ eof ******************************************/
