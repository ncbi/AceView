/*  File: flag.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: flag.c,v 1.4 2007/03/30 23:57:15 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr  3 14:36 1997 (rd)
 * Created: Mon Jan 20 23:53:24 1997 (rd)
 *-------------------------------------------------------------------
 */

#include "acedb.h"
#include "array.h"
#include "flag.h"
#include "lex.h"

static DICT *masterDict = 0 ;
static FLAGSET flagValue[32] ;
static Array dictList ;
static BOOL charOK[256] ;

/************************************************************/

DICT* flagSetDict (const char *set)
{ 
  int i = -1 ;

				/* general initialisation if necessary */
  if (!masterDict)
    { AC_HANDLE h = handleCreate() ;

      masterDict = dictCreate (16) ;
      dictList = arrayCreate (16, DICT*) ;
      for (i = 0 ; i < MAXFLAG ; ++i)
	flagValue[i] = 1 << i ;
				/* build charOK[] character set */
      for (i = 0 ; i < 26 ; ++i)
	{ charOK['a'+i] = TRUE ;
	  charOK['A'+i] = TRUE ;
	}
      for (i = 0 ; i < 10 ; ++i)
	charOK['0'+i] = TRUE ;
      charOK['_'] = TRUE ;
      
      for (i = 0 ; i < MAXFLAG-1 ; ++i)
	flag ("Default", strnew (messprintf ("bit%d", i+1), h)) ;

      handleDestroy (h) ;
    }

  if (dictAdd (masterDict, set, &i))
    { DICT *d = dictCreate (32) ;
      int len = strlen (set) ;
      KEY key ;
      
      array(dictList, i, DICT*) = d ;
      for (key = 0 ; key < lexMax(0) ; ++key)
	if (!strncmp (name(key), "Flag", 4) &&
	    !strncmp (name(key)+4, set, len) &&
	    name(key)[4+len] == '#')
	  dictAdd (d, name(key)+5+len, 0) ;
    }

  if (i == -1)
    return 0 ;

  return arr(dictList, i, DICT*) ;
}

/************************************************************/

/* s has form "flagstring1|flagstring2|...|flagstringN" */

FLAGSET flag (char *set, char *s)
{ 
  FLAGSET res = 0 ;
  int i ;
  char *t ;
  static char *lastSet = 0 ;
  static DICT *dict ;

  if (set != lastSet)		/* efficiency */
    { dict = flagSetDict (set) ;
      lastSet = set ;
    }

  if (!dict || !s)
    messcrash ("bad arguments to flag()") ;

  if (dictFind (dict, s, &i))	/* shortcut attempt for simple cases */
    return flagValue[i] ;

  while (*s)
    { for (t = s ; *t && *t != '|' ; ++t) /* look for separator */
	if (!charOK[(int)*t])
	  { messerror ("Bad character %c (0x%x) in flag", *t, *t) ;
	    return 0 ;
	  }
      if (*t)
	*t = 0 ;
      else
	t = 0 ;

      if (dictAdd (dict, s, &i)) /* process this string */
	{ KEY dummy ;
	  lexaddkey (messprintf ("Flag%s#%s", set, s), &dummy, 0) ;
	  if (i >= MAXFLAG - 1)
	    messcrash ("flag set %s has overflowed, with %s", set, s) ;
	}
      res |= flagValue[i] ;

      if (t)			/* restore and cycle if necessary */
	{ *t = '|' ;
	  s = t+1 ;
	}
      else 
	break ;
    }

  return res ;
}

/************************************************************/

char *flagNameFromDict (DICT *dict, FLAGSET flag)
{ 
  int i ;
  static Stack s = 0 ;

  if (!dict)
    return 0 ;

  s = stackReCreate (s, 128) ;

  for (i = 1 ; i <= dictMax(dict) ; ++i)
    if (flag & flagValue[i])
      { if (!stackEmpty(s))
	  catText (s, "|") ;
	catText (s, dictName (dict, i)) ;
      }

  return stackText (s, 0) ;
}

/************************************************************/

char *flagName (char *set, FLAGSET flag)
{ 
  DICT *dict = flagSetDict (set) ;

  return flagNameFromDict (dict, flag) ;
}

/************************************************************/
 
