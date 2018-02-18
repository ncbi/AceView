/*  File: lexalpha.c
 *  Author: Jean Thierry-Mieg (mieg@kaa.crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jul  9 11:28 2003 (rnc)
 * * Sep  1 23:24 1998 (rd): restructure keySetAlphaHeap() a bit
 *	remove keySetCopy if keySetmax == 1 - not correct if key not visible
 *	follow aliases if keySetMax() < LIMIT ; this always happened if > LIMIT
 *	  because lexSetStatus() follows aliases.  Also it is desirable.
 *	for > LIMIT, transfer responsibility for lexIsVisibleKey() check to
 *	  alphaMake(), and change key names to avoid errors on old caches.
 *        This is cleaner, and more efficient.
 *	alphaMake() always makes a keyset - safer, free, and perhaps needed now
 *	added lexMaxVisible() for command.c
 * Created: Tue Nov 22 11:34:29 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: lexalpha.c,v 1.6 2009/02/25 15:18:31 mieg Exp $ */

#define LIMIT 300
#define MAXCLASS 256
#include "acedb.h"
#include "lex.h"
#include "a.h"
#include "whooks/sysclass.h"
#include "pick.h"
#include "utils.h"

static KEYSET alphaSet[MAXCLASS] ;
static void alphaGet(int c) ;

static int classAlphaOrder (const void *a, const void *b)
{
  int c1 = *(const int*)a , c2 = *(const int*) b ;
  return 
    lexstrcmp (className(KEYMAKE(c1, 0)), className(KEYMAKE(c2,0))) ;
}

  /* return the alphabetic first nn keys of ks */
KEYSET lexAlphaClassList (BOOL redo)
{
  static KEYSET ks = 0 ;
  int i ;

  if (redo)
    keySetDestroy (ks) ;
  if (!ks)
    {
      ks = keySetCreate () ;
      i = MAXCLASS ;
      while (i--) keySet(ks, i) = i ;
      arraySort (ks, classAlphaOrder) ;
    }
  return ks ;
}

KEYSET keySetAlphaHeap (KEYSET ks, unsigned int nn)
{
  char ff[MAXCLASS] ; /* flags the class present in ks */
  unsigned int i, c , c1, j ;
  KEYSET ksa ;
  KEY *kp, *kpa, k ;
  KEYSET alphaClassList ;

  if (!keySetExists(ks))
    return 0 ;   /* may be messcrash ? */

  if (keySetMax(ks) == 0 || nn == 0) /* nothing to sort */
    return keySetCreate();    /* return empty keyset */

  memset(ff, 0, MAXCLASS) ;

  if (keySetMax(ks) < LIMIT)	/* do directly */
    { 
      ksa = arrayCreate (keySetMax(ks), KEY) ;
      keySet (ksa, keySetMax(ks)) = 0 ; /* make room, is that needed after arrayCreate ? */ 
   /* follow aliases and remove invisible keys */
      kp = arrp (ks, 0, KEY) - 1 ;
      i = keySetMax (ks) ; j = 0 ;
      while (kp++, i--)
	{ k = lexAliasOf(*kp) ;	/* follow aliases */
	  if (lexIsKeyVisible(k))
	    arr (ksa, j++,KEY) = k ; /* use a macro, space is guaranteed */
	}
      keySetMax(ksa) = j ;      
				/* then sort and return */
      arraySort(ksa,keySetAlphaOrder) ;
      arrayCompress (ksa) ;  /* avoid duplicates */
      return ksa ;
    }

	/* keySetMax > LIMIT
	   set flags and then loop through presorted keysets = alphaSet[c]
	*/


  kp = arrp (ks, 0, KEY) ;
  for (i = keySetMax(ks) ; i-- ; ++kp)
    { c = class(*kp) ;
      if (!ff[c])
	{ ff[c] = 1 ;
	  lexClearClassStatus(c, ALPHASTATUS) ;
	}
      if (*kp)			/* michel, dec 98 avoid *kp = 0 */
	lexSetStatus (*kp, ALPHASTATUS) ; /* does follow aliases */
      else
	/* messcrash ("bad key 0") */ ;  /* RD 990107: should messcrash really */
    }

  if (nn > keySetMax(ks))
    nn = keySetMax(ks) ;

  ksa = arrayCreate (nn, KEY) ;
  keySet(ksa, nn - 1) = 0 ; /* make room */
  kpa = arrp(ksa, 0, KEY) ;
  j = 0 ;

  alphaClassList = lexAlphaClassList (FALSE) ;
  for (c1 = 0 ; c1 < MAXCLASS ; c1++)
    {
      c = keySet(alphaClassList, c1) ;
      if (ff[c])
	{ alphaGet(c) ; /* creates alphaSet */
	kp = arrp(alphaSet[c], 0, KEY) ;
	for (i = keySetMax (alphaSet[c]) ; i-- ; ++kp)
	  if ((lexGetStatus(*kp) & ALPHASTATUS))
	    { *kpa++ = *kp ;
	    if (++j == nn)
	      goto ok ;  /* break all loops */
	    }
	}
    }
ok:
  keySetMax(ksa) = j ;
  return ksa ;
}

/******************************************/
/******************************************/

static BOOL firstPass = TRUE ;
static char isMarked[MAXCLASS] ; /* class changed - rebuild needed */
static char isTouched[MAXCLASS] ; /* lexalpha changed - needs saving */

static void alphaMake (int c) 
{
  int i, max = lexMax(c) ;
  KEYSET ks ;
  KEY kk, *kp ;

  if (max)
    {		/* make a keySet of visible keys in the class */
      ks = arrayCreate (max, KEY) ;
      keySet(ks, max - 1) = 0 ;	/* make space */
      for (i = max, kp = arrp(ks, 0, KEY), kk = KEYMAKE(c,0) ; i-- ; kk++)
	if (lexIsKeyVisible (kk))
	  *kp++ = kk ;
      keySetMax(ks) = kp - arrp(ks,0,KEY) ;

				/* sort it */
      arraySort (ks, keySetAlphaOrder) ;
    }
  else
    ks = keySetCreate() ;	/* make an empty keyset */

  keySetDestroy (alphaSet[c]) ;
  alphaSet[c] = ks ;
  isTouched[c] = 1 ;
  isMarked[c] = 0;
}
  
/***********************/

static void alphaRead (int c) 
{ KEY key ;

  lexaddkey(messprintf("_lexalpha%d",c), &key, _VVoc) ;
  alphaSet [c] = arrayGet (key, KEY, "k") ;
  isTouched [c] = 0 ;
}

/***********************/

static void alphaGet (int c) 
{ int i ;

  if (firstPass)
    { firstPass = FALSE ;
      i = MAXCLASS ;
      while (i--)
	{ isTouched[i] = 0 ;
	  isMarked [i] = 0 ;
	  alphaSet[i] = 0 ;
	}
    }
  
  if (isMarked[c])  /* obsolete on disk */
    alphaMake (c) ;
  else
    { if (!alphaSet[c])
	alphaRead (c) ;
      if (!alphaSet[c])
	alphaMake (c) ;
    }
}
  
/***********************/
/***********************/

void lexAlphaMakeAll (void)
{ int c = MAXCLASS ;
  while (c--)
    if (lexMax(c))
      alphaMake (c) ;
}

/***********************/

void lexAlphaSaveAll (void)
{ int c = MAXCLASS ;
  KEY key ;
  KEYSET ks ;

  while (c--)
    { lexaddkey(messprintf("_lexalpha%d",c), &key, _VVoc) ;
	/* rd - must call lexaddkey() before isMarked[], else on
	   c == _VVoc it can change the value !!!!
	*/
      if (isMarked[c])
	alphaMake (c) ;
      ks = alphaSet [c] ;
      if (ks && isTouched [c]) 
	{ arrayStore (key, ks,  "k") ;
	  isTouched[c] = 0 ;
	}
    }
}

/***********************/

void lexAlphaShutDown (void)
{
  int i = 256 ;
  
  while (i--)
    arrayDestroy (alphaSet[i]) ;
}

/***********************/

void lexAlphaMark (int c)
{ 
  if (firstPass)
    { int i = MAXCLASS ;
      firstPass = FALSE ;
      while (i--)
	{ isTouched[i] = 0 ;
	  isMarked [i] = 0 ;
	  alphaSet[i] = 0 ;
	}
    }

  if (c >= 0 && c < MAXCLASS)
    isMarked [c] = 1 ;		/* RD 980901 - no need to destroy here 
				   and potentially bad for lexAlphaIterator*() */
}

/***********************/

void lexAlphaClear (void)
{ int i = MAXCLASS ;
  while (i--)
    { isTouched[i] = 0 ;
      isMarked [i] = 0 ;
      if (keySetExists(alphaSet[i]))
	keySetDestroy (alphaSet[i]) ;
    }
}

/***********************/

int lexMaxVisible (int c)
{
  alphaGet(c) ;
  return keySetMax(alphaSet[c]) ;
}

/***********************/
/***********************/

#ifdef LEX_ALPHA_ITERATOR	/* needs more work */

struct LexAlphaIteratorStruct
{
  int c ;
  int i ;
  unsigned char mask ;
} ;

LexAlphaIterator *lexAlphaIteratorCreate (KEY classe, AC_HANDLE h)
{
  LexAlphaIterator lat ;
  unsigned char mask ;

  if (!pickIsA (&classe, &mask))
    return 0 ;

  lat = (LexAlphaIterator) halloc (sizeof(LexAlphaIteratorStruct), h) ;
  lat->c = classe ;
  lat->i = -1 ;
  lat->mask = mask ;

  alphaGet(lat->c) ;

  return lat ;
}

KEY lexAlphaIteratorNext (LexAlphaIterator lat)
{
  KEY k ;

  if (!alphaSet[lat->c])
    messcrash ("lost alphaSet[%s] in lexAlphaIteratorNext()", 
	       pickClass2Word(lat->c)) ;

  while (TRUE)
    { ++lat->i ;
      if (lat->i >= arrayMax(alphaSet[lat->c]))
	return 0 ;
      k = arr(alphaSet[lat->c],i,KEY) ;
      if (!mask || mask == (mask & q->isMask)) /* q is private to lexsubs4! */
	return k ;
    }
  return 0 ;			/* for compiler happiness */
}

#endif

/***********************/

unsigned int keySetCountVisible (KEYSET kset)
{
  unsigned int max = 0, i ;

  if (kset && keySetMax(kset))
    { KEY k, *kp = arrp (kset, 0, KEY) - 1 ;
      
      i = keySetMax(kset) ;
      while (kp++, i--)
	{ 
	  k = lexAliasOf(*kp) ;	/* follow aliases */
	  if (lexIsKeyVisible(k)) max++ ;
	}
    }
  return max;
}

/***********************/
/***********************/

int keyNamesAverage (KEYSET kSet, int *goodEntries)
{
  int length = 0;
  const char *cp ;
  int		i;
  
  *goodEntries = 0;
  /* count the number of usable items in this list 
     and determine the average length for names in the list */
  if (kSet)
    {
      KEY *kp = arrp (kSet, 0, KEY) ;
      
      i = keySetMax(kSet) ;
      while (kp && i--)
	{
	  if (lexIsKeyVisible(*kp))
	    {
	      *goodEntries += 1 ;
	      cp = 0 ;
	      if (nextName (*kp, &cp) && *cp != '\177')
		length += strlen (cp);
	    }
	  kp++;
	}
      if (*goodEntries)
	length = length / *goodEntries;
    }
  return length;
}

/*********************** end of file *********************/
