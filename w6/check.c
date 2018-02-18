/*  File: check.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 31 13:43 1994 (srk)
 * Created: Wed Sep 16 13:33:07 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: check.c,v 1.3 2014/11/10 18:05:09 mieg Exp $ */

#include <ctype.h>

#include "acedb.h"
#include "systags.h"
#include "check.h"


BOOL checkRegExpSyntax (char *cp)
{
  if (!cp || ! *cp)
    return FALSE ;

  return TRUE ;
}

BOOL checkText (char *text, KEY cst)
{ char c, *cp, *cq = text ;

  cp = name (cst) ;

  if (!strncmp(cp, "FORMAT ",7)) 
    { cp += 7 ;
      
      while ((c = *cp++))
	switch (c)
	  {
	  case 'a':
	    if (!isalpha((int)*cq))
	      return FALSE ;
	    cq++ ;
	    break ;
	  case 'n':
	    if (!isdigit((int)*cq))
	      return FALSE ;
	    cq++ ;
	    break ;
	  default:
	    if (c != *cq)
	      return FALSE ;
	    cq++ ;
	    break;  
	  }
    }
  return TRUE ;
}

/************************************************/

BOOL checkKey (KEY new, KEY cst)
{ 

  
  
  return TRUE ;
}

/************************************************/

BOOL checkData (void *xp, KEY type, KEY cst)
{ int i ; 
  float x = 0 ;
  mytime_t tm ;

  switch (type)
    {
    case _Int: 
      i = *(int*) xp ;
      x = i ;
      break ;
    case _Float: 
      i = *(float*) xp ;
      x = i ;
      break ;
    case _DateType: 
      tm = *(mytime_t*) xp ;
      x = tm ;
      break ;
    default:
      messcrash ("badtype in checkData x=%f", x) ;
    }

  
  
  
  return TRUE ;
}

/************************************************/
/************************************************/
