/*  File: keysetdump.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Needed in acequery and xace
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 19 10:30 1998 (fw)
 * * Aug 19 17:36 1998 (rd): changed condition on dumping keys from
 *	(!pickList[class(*kp)].protected && !lexIsKeyInvisible(*kp))
 *      to (lexIsKeyVisible (*kp)) and prevented lexaddkey() on parsing
 *	protected keys, silently.
 * * Nov 12 14:23 1991 (mieg): added heap system
 * * Oct 23 18:12 1991 (mieg): added keySetParse
 * Created: Mon Oct 21 19:55:34 1991 (mieg)
 *-------------------------------------------------------------------
 */


/* $Id: keysetdump.c,v 1.8 2009/02/25 15:18:31 mieg Exp $ */

#include "acedb.h"
#include "keyset.h"
#include "bs.h"	  /* need bsKeySet */
#include "lex.h"  /* needed in keysetdump and keysetorder */
#include "pick.h" /* needed in keysetdump */
#include "a.h"    /* needed in keySetParse */
#include "freeout.h"
#include "sysclass.h"
#include "java.h"
#include "dump.h"

/************************************************************/
  /* Heap sort
   * source is an array to be partially sorted,
   * nn is the desired number of sorted objects,
   * order is the order function
   * min and max are extrema of the sorting functions
   */

/**************************************************************/
/************ Dedicated heap function *************************/

KEYSET  keySetHeap(KEYSET source, unsigned int nn, int (*order)(KEY *, KEY *))
{ KEYSET h, result ;
  int i , k = 1 , max , pos ;
  KEY  new ;
    
  if(!source || !arrayMax(source))
    return keySetCreate() ;
  
  if(!arrayExists(source))
    messcrash("keySetHeap received a non-magic source keySet") ;

  messStatus ("Heaping") ;
    /* Note that KEY 0 is used as maximal element */
  if(nn > arrayMax(source))
    nn = arrayMax(source) ;
  while( ( 1 << ++k) <= nn) ;
  max = 1 << k ;
  h = arrayCreate(max, KEY) ;
  arrayMax(h) = max ; 
  
  i = arrayMax(source) ;
  while (i --)
    { new = keySet(source, i) ;

        /*  HEAP_PUSH(new) */
      { pos = 1 ;
	while( ! arr(h, pos, KEY) ||
	      (order(&new,arrp(h, pos, KEY)) == -1 ) )
	  { if (pos>>1)
	      arr(h, pos>>1, KEY) = arr(h, pos, KEY) ;
	    arr(h, pos, KEY) = new ;
	    if ( (pos <<= 1) >= max )
	      break ;
	    if( ! arr(h, pos + 1, KEY) ||
	       ( order( arrp(h, pos, KEY), arrp(h, pos + 1 ,KEY)) == -1) )
	      pos++ ;
	  }
      }
    }

/* Test print out 
  for(i=0 ; i < arrayMax(h) ; i++)
    printf("%d %s \n", i, name(arr(h,i,KEY))) ;
*/  
              /* To dump,
	       * Create a minimal key
	       * and push it repeatedly
	       * collecting the top key 
	       */
  result = arrayCreate(max, KEY) ;
  arrayMax(result) = max - 1 ;
  
  i = max - 1 ;  /* pos 0 is not filled with any relevant data */
  new = KEYMAKE(0, lexMax(0)) ;
  while(i--)
    { arr(result, i, KEY) = arr(h, 1, KEY) ;
        /*  HEAP_PUSH(minimal key) */
      { pos = 1 ;
	while( TRUE ) /* I pretend to push a minimal key */
	  { if (pos>>1)
	      arr(h, pos>>1, KEY) = arr(h, pos, KEY) ;
	    arr(h, pos, KEY) = new ;
	    if ( (pos <<= 1) >= max )
	      break ;
	    if( arr (h, pos, KEY) == new ||
	       (order( arrp(h, pos, KEY), arrp(h, pos + 1 ,KEY)) == -1))
	      pos++ ;
	  }
      }
    }

  if(arrayMax(source) < arrayMax(result))
    (arrayMax(result) = arrayMax(source)) ;
  
  keySetDestroy(h) ;

  return result ;
}

/*************************************/

KEYSET keySetNeighbours(KEYSET oldSet)
{
  KEYSET new = keySetCreate() , nks ; KEY key, *kp ;
  int  i = keySetMax(oldSet) , j1 = 0, j2;
 
  i = keySetMax(oldSet) ;
  while(i--)
    { nks = (KEYSET) bsKeySet(keySet(oldSet, i)) ;
      if (nks)
	{ j2 = keySetMax(nks) ; kp = arrp(nks,0,KEY) - 1 ;
	  while(kp++, j2--)
	    if (!pickXref(class(*kp)) && 
		iskey(*kp) == 2 &&
		class(*kp) != _VUserSession)
	      keySet(new, j1++) = *kp ;
	  keySetDestroy(nks) ;
	}
      if (messIsInterruptCalled())
	break ;
    }  
  i = keySetMax(oldSet) ;
  while(i--)
    {
      key = keySet(oldSet, i) ;
      if (class(key) == _VKeySet &&
	  (nks = arrayGet (key,KEY,"k")))
	{ 
	  j2 = keySetMax(nks) ; kp = arrp(nks,0,KEY) - 1 ;
	  while(kp++, j2--)
	    if (!pickXref(class(*kp)) && 
		iskey(*kp) == 2 &&
		class(*kp) != _VUserSession)
	      keySet(new, j1++) = *kp ;
	  keySetDestroy(nks) ;
	}
    }
  keySetSort(new) ;
  keySetCompress(new) ;
  keySetDestroy (oldSet) ;
  return new ;
}

/**************************************************************/
/**************************************************************/

       /* first order classes */
     /* then compares KEYs acoording to their Word meaning*/

int keySetAlphaOrder (const void *a, const void *b)
{
  register int i = class(*(KEY*)a), j = class(*(KEY*)b) ;

  if (i != j)   
    return lexstrcmp (className(*(KEY *)a),
                   className(*(KEY *)b)  ) ;
  else
    return lexstrcmp(name(*(KEY *)a),
                   name(*(KEY *)b)  );
}

/**************************************************************/
BOOL keySetDump(FILE *f, Stack buffer, KEYSET s)
{
  int i = keySetMax(s) ;
  KEY *kp  = arrp(s,0,KEY)  - 1 ; 

  if(s->size != sizeof(KEY))
    return FALSE ;
  if (f)
    { while (kp++, i--)
	if (lexIsKeyVisible(*kp))
	  fprintf(f, "%s %s\n",
		  className(*kp),
		  freeprotect(name(*kp)) ) ;
    }
  else if (stackExists(buffer))
    { while (kp++, i--)
	if (lexIsKeyVisible(*kp))
	  catText(buffer,
		  messprintf("%s %s\n",
			     className(*kp),
			     freeprotect(name(*kp)))) ;
    }
  else
    { while (kp++, i--)
	if (lexIsKeyVisible(*kp))
	  freeOutf ("%s %s\n",
		    className(*kp),
		    freeprotect(name(*kp))) ;
    }
  return TRUE ;
}

/************************************/

BOOL javaDumpKeySet(KEY key) 
{ KEYSET s = arrayGet(key,KEY,"k") ;
  int i ;
  KEY *kp;
  BOOL first = TRUE ;

  if (!s) 
    return FALSE;

  freeOutf("?Keyset?%s?\t",freejavaprotect(name(key)));
 
  i = keySetMax(s) ;
  kp = arrp(s,0,KEY) - 1 ;
  while (kp++,i--)
    if (lexIsKeyVisible(*kp)) 
      { 
	if (!first) 
	  freeOut("\t") ;
	first = FALSE ;
	freeOutf("?%s?%s?\n",className(*kp),
		 freejavaprotect(name(*kp)));
      }
  
  freeOut("\n");

  arrayDestroy(s);
  return TRUE ;
}	

/************************************/
/* keyset dumpfunc is new interface for dump functions */
BOOL keySetDumpFunc (FILE *f, Stack buffer, KEY key)
{ 
  KEYSET s = arrayGet(key, KEY, "k");
  BOOL result;

  if (!s)
    return FALSE;

  result = keySetDump(f, buffer, s);
  arrayDestroy(s);
  return result;
}

/************************************/

BOOL keySetParse (int level, KEY key)
{
  char *cp = 0 , c;
  KEYSET kS = keySetCreate() ;
  int  classe ; KEY k ;

  while (freecard(level) && (cp = freewordcut(":-/\" ",&c)))
    {
      classe = pickWord2Class(cp) ;
      if (!classe)
	goto abort ;
      if (pickList[classe].protected) /* can't create protected keys */
	continue ;
      freestep(':') ;
      if((cp = freeword()))
	{ lexaddkey(cp,&k,classe) ;
	  keySetInsert(kS,k) ;
	}
      else
	goto abort ;
    }
  if (keySetMax(kS))
    arrayStore(key,kS,"k") ;
  keySetDestroy(kS) ;
  return TRUE ;
  
 abort:
  messerror ("keySetParse error at  Line %7d  KeySet %s, bad word %s\n", 
	    freestreamline(level), name(key),  cp ? cp : "void line") ;
  keySetDestroy(kS) ;
  return FALSE ;
}

/**************************************************************/

BOOL keySetAceDump (FILE *f, KEYSET s)
{ 
  int i ;
  
  for (i = 0 ; i < keySetMax(s) ; ++i)
    dumpKey (keySet(s,i), f, 0) ; /* dumps in f OR s OR freeout */
  return TRUE ;
}

/*********************************************/
/*********************************************/


/* read a .ace formatted keyset listing as generated by the List command
   and return it as a keyset data structure
   original by Detlef Wolf, DKFZ
   new version 951017 from Detlef adpated in minor ways by rd

   accept  Locus a;b;c;Sequence e;f  // all on one or multiple lines
   first line: Keyset toto is optional
*/

KEYSET keySetRead (int level) 
{ KEYSET ks = keySetCreate() ;
  int n, table = 0 , tt, nbad = 0 ;
  KEY key;
  char *word ;
    
  while (freecard (level))   /* closes file itself */
    { 
      while ((word = freeword()))
        { 
	  if ((tt = pickWord2Class (word)))
	    { table = tt ;
	    freestep (':') ;
	    word = freeword() ;
	    }
	  if (tt) table = tt ;
	  if (!table)
            { ++nbad ; continue ; }
          if (!word)
            { freeOutf ("// No object name, line %d", 
			freestreamline (level)) ;
              continue ;
            }
          if (table == _VKeySet && !keySetMax(ks))
            continue ;
          if (lexword2key (word, &key, table))
            keySet(ks, keySetMax(ks)) = key ;
          else
            ++nbad ;
        }
    }
  freeclose (level) ; /* a guarantee */
  if (nbad)
    freeOutf ("// %d objects not recognized", nbad) ;
				/* finally sort new keyset */
  n = keySetMax(ks) ;
  keySetSort(ks) ;
  keySetCompress(ks) ;
  if (n - keySetMax(ks)) 
    freeOutf ("// %d duplicate objects removed", n - keySetMax(ks)) ;
  return ks ;
}


/*************************************************************/
/*************************************************************/

void keySetKill (KEYSET ks) 
{
  int i ;
  KEY key ;

  i = keySetExists (ks) ? keySetMax (ks) : 0 ;
  while (i--)
    { 
      key = keySet (ks, i) ;
      if (class(key))
	switch ( pickType(key))
	  {
	  case 'A':
	    arrayKill (key) ;
	    break ;
	  case 'B':
	    { OBJ obj ;
	    if ((obj = bsUpdate (key)))
	      bsKill (obj) ;	      
	    }
	    break ;
	  }
    }
}

/**************************************************************/
/**************************************************************/
