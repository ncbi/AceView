/*  File: class.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
    This file deals with the class hierarchy
    The MainClasses are listed and initialised in w5/tags.c
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 19 11:11 1998 (fw)
 * Created: Sun Oct 11 22:35:44 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: class.c,v 1.6 2014/11/30 03:21:21 mieg Exp $ */


#include <wh/acedb.h>

#include <wh/menu.h>
#include <wh/a.h>
#include <wh/bs.h>
#include <wh/bindex.h>
#include <whooks/systags.h>
#include <whooks/sysclass.h>
#include <wh/pick.h>
#include <wh/lex_bl_.h>		/* since need lexmark() */
#include <wh/lex.h>
#include <wh/query.h>

typedef struct CTstructure *CTP ;
typedef struct CTstructure 
 { KEY key ;
   CTP sta ;
   BOOL isTopClass ;
   int color ; 
   int generation , x, y, len , box ;
 } CT ;


/**************************************************************/
/**************** Constraints verifications ******************/
/*************************************************************/

static BOOL checkLoops (KEYSET ks, int top)
{
  OBJ obj ;
  static int nWarnings = 0 ;
  static KEYSET children = 0 ;
  KEY child, parent ;
  int i, index ;
  BOOL thisOne ;

  if (!lexword2key (pickClass2Word(top), &parent, _VClass))
    { if (nWarnings++ < 3)
	messout ("Can't find the Class entry %s of class %d, this is a programming bug, sorry ", 
		 pickClass2Word(top), top) ;
      if (nWarnings++ == 3)
	messout("Further similar warnings suppressed") ;
      return FALSE ;
    }
  keySet(ks, 0) = KEYMAKE(_VClass, 0) ;	/* for keySetOrder */
  children = keySetReCreate (children) ;

/* strategy is to move to the parent set only those children 
   all of whose parents are already there.  Keep both sets as
   regular sorted keysets.
*/

  while (parent)
    { keySetInsert (ks, parent) ;

      if (!(obj = bsCreate (parent)))
	{ messout ("Can't find the Class object for %s:%s", 
		   className(parent), name(parent)) ;
	  return FALSE ;
	}
      if (bsGetKey (obj, _Is_a_superclass_of, &child)) do
	if (!keySetFind (children, child, &index))
	  keySetInsert (children, child) ;
      while (bsGetKey (obj, _bsDown, &child)) ;
      bsDestroy (obj) ;

      parent = 0 ;
      for (i = 0 ; i < keySetMax(children) ; ++i)
	{ child = keySet(children, i) ;
	  if (!(obj = bsCreate(child)))
	    { messout ("Can't find the Class object for %s:%s", 
		       className(child), name(child)) ;
	      return FALSE ;
	    }
	  thisOne = TRUE ;
	  if (bsGetKey (obj, _Is_a_subclass_of, &parent)) do
	    if (!keySetFind (ks, parent, &index))
	      thisOne = FALSE ;
	  while (bsGetKey (obj, _bsDown, &parent)) ;
	  bsDestroy (obj) ;
	  if (thisOne)
	    { keySetRemove (children, child) ;
	      parent = child ;
	      break ;
	    }
	  else
	    parent = 0 ;
	}
    }

  keySet (ks, 0) = top ;

  if (keySetMax (children))
    return FALSE ;
  else
    return TRUE ;
}

static int findInStack(Stack s, char *text)
{ char *cp ;
  int i = 0 ;

  stackCursor(s, 0) ;
  while (i++, cp = stackNextText(s))
    if (!strcmp(cp, text))
      return i ;
  pushText(s, text) ;
  return i ;
}

static int checkFilters(KEYSET ks, int nn, BOOL makeIt)
{
  int i, n, max = keySetMax(ks) ;
  KEY classe, isa, ft ; 
  Stack s, new ; 
  BOOL modif = FALSE ;
  char *cp, *cq , buffer[32] ;
  OBJ obj, obj1 ;

  if (max <= 2 )
    return TRUE ;

  new = stackCreate(20) ;
  stackTextOnly(new); /* since we want to save it */
  pushText(new, "toto") ;
  
      /* asociate a mask to each filter */
  for (i = 1 ; i < max ; i++)
    if ((obj = bsCreate(classe = keySet(ks, i))))
      {
	if (bsGetData(obj, _Filter, _Text, &cp))
	  { 
	    /* There is a bug in the 
	       query code which references off the end of the string
	       passed to condCheckSyntax, and triggers electric fence.
	       It doesn't otherwise cause problems, and the code in
	       V. strange, so I out this work-around in.... SRK */
	    char *cp1 = messalloc(strlen(cp)+200);
	    n = findInStack( new, cp) ;  /* will pushText if not already there */
	    strcpy(cp1, cp);
	    if (!condCheckSyntax(cp1))
	      {	
		bsDestroy(obj) ;
		messfree(cp1);
		return FALSE ;
	      }
	    messfree(cp1);
	    
	    if (n > 9999) /* since there is a dummy entry on the stack */
	      { messout("More than 8 eight filters in class %s",
			pickClass2Word(nn)) ;
		bsDestroy(obj) ;
		return FALSE ;
	      }
	  }
	bsDestroy(obj) ;
      }
  
    /* check if one should rehash */
  
  sprintf(buffer,"__filters%d", nn) ;
  lexaddkey(buffer, &ft, _VCalcul) ;
  
  s = pickList[nn].filters = stackGet(ft) ;
  if (!s)
    { modif = TRUE ;
      if (makeIt)
	s = pickList[nn].filters = new ;
    }
  else
    { stackCursor(s, 0) ;
      stackCursor(new, 0) ;
      while ((cp = stackNextText(s)))
	{ if (!(cq = stackNextText(new)) ||
	      strcmp(cq, cp))
	    { modif = TRUE ;
	      break ;
	    }
	}
      if (stackNextText(new))
	modif = TRUE ;

      if (makeIt)
	{ stackDestroy(s) ;
	  pickList[nn].filters = new ;
	}
    }
      
  for (i = 1 ; i < max ; i++)
    if ((obj = bsUpdate(classe = keySet(ks, i))))
      { n = 0 ;
	if (bsGetData(obj, _Filter, _Text, &cp))
	  n = 1 << (findInStack( new, cp) - 2) ;  /* n == 1 == 1<< 00 on first filter */

	if (bsGetKey(obj, _Is_a_subclass_of, &isa))
	  do
	    { int ii ;
	      obj1 = bsCreate(isa) ;
	      if (bsGetData(obj1, _Mask, _Int, &ii))
		n |= ii ;
	      else if (bsFindTag(obj1, _Is_a_subclass_of))
		{ messout("Class %s cannot find previous mask of %s", name(classe),name(isa)) ;
		  bsDestroy(obj1) ;
		  bsDestroy(obj) ;
		  return FALSE ;
		}
	      bsDestroy(obj1) ;
	    } while (bsGetKey(obj, _bsDown, &isa)) ;  /* do while construct */

	bsAddData (obj, _Mask, _Int, &n) ;
	n = keySet(ks, 0) ;
	bsAddData(obj, _Belongs_to_class, _Int, &n) ;
	bsSave(obj) ;
      }
  
  if (!arrayExists (pickList[nn].conditions))
    pickList[nn].conditions = arrayCreate(8, void *) ;
  queryIsAInit (pickList[nn].filters, 
		pickList[nn].conditions) ;

  return modif ? 2 : 1 ;
}

static BOOL classSortClass(int i)
{ KEY k = 0 ;
  OBJ obj ;
  Array a = pickList[i].conditions ;
  Stack s = pickList[i].filters ;
  char buffer[32] ;

  if (!stackExists(s) || !arrayMax(a))
    return FALSE ;

  messStatus (messprintf ("Sorting %s\n", pickClass2Word(i))) ;

  while (lexNext(i, &k))
    if ((obj = bsCreate(k)))
      { bIndexObject (0, obj) ; bsDestroy (obj) ; }
  lexmark(i) ;
  sprintf(buffer,"__filters%d", i) ;
  lexaddkey(buffer, &k, _VCalcul) ;
  stackStore(k, s) ;
  return TRUE ;
}

BOOL classSortAll(void)
{ KEYSET ks = 0 ;
  int i ;

/*
  if (!checkClasses(FALSE))
    return FALSE ;
*/
    
  i = 256 ;
  while (i--)
    if (pickList[i].type == 'B')
      { ks = keySetReCreate(ks) ;
	if (checkLoops (ks, i) && 
	    checkFilters(ks, i, TRUE) == 2)
	  classSortClass(i) ;
      }
  keySetDestroy(ks) ;
  return TRUE ;
}
  

/*************************************************************/
/*************************************************************/

