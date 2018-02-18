/*  File: pepdisp.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: generic methods for sequence displays
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 27 09:45 1998 (edgrif)
 * * Jun 24 13:37 1998 (edgrif): Fixed initialisation bug in methodInitialise,
        isDone was not being set. Made methodInitialise an internal routine.
	Inserted missing call to methodInitialise in method routine.
 * Created: Wed Jun  1 15:05:00 1994 (rd)
 *-------------------------------------------------------------------
 */

/*  $Id: method.c,v 1.2 2005/08/24 23:02:04 mieg Exp $  */

#include "acedb.h"
#include "method.h"
#include "lex.h"
#include "classes.h"
#include "bs.h"
#include "tags.h"
#include "systags.h"

/****** global information array *******/

Array methodInfo = 0 ;

/****** keys used here - do they need to be global? *****/

static KEY _Width ;
static KEY _Frame_sensitive ;
static KEY _Strand_sensitive ;
static KEY _Show_up_strand ;  /* mieg */
static KEY _Blastn ;
static KEY _Blixem_X ;
static KEY _Blixem_N ;
static KEY _Blixem_P ;
static KEY _Belvu ;
static KEY _Bumpable ;
static KEY _Score_bounds ;
static KEY _Score_by_offset ;
static KEY _Score_by_width ;
static KEY _Score_by_histogram ;
static KEY _BlastN ;
static KEY _Percent ;
static KEY _Right_priority ;
static KEY _Web_Right_priority ;
static KEY _EMBL_dump ;


/* Internal routines.                                                      */
static void methodInitialise (void) ;




/*                                                                           */
/* This is an internal routine, code using the method package does not need  */
/* to worry about initialisation, all interface routines call this routine   */
/* to make sure the package is initialised.                                  */
/*                                                                           */
static void methodInitialise (void)
{ 
  static int isDone = FALSE ;

  if (isDone == FALSE)
    {
    lexaddkey ("Width", &_Width, 0) ;
    lexaddkey ("Frame_sensitive", &_Frame_sensitive, 0) ;
    lexaddkey ("Strand_sensitive", &_Strand_sensitive, 0) ;
    lexaddkey ("Show_up_strand", &_Show_up_strand, 0) ; /* mieg */
    lexaddkey ("Blastn", &_Blastn, 0) ;
    lexaddkey ("Blixem_X", &_Blixem_X, 0) ;
    lexaddkey ("Blixem_N", &_Blixem_N, 0) ;
    lexaddkey ("Blixem_P", &_Blixem_P, 0) ;
    lexaddkey ("Belvu", &_Belvu, 0) ;
    lexaddkey ("Bumpable", &_Bumpable, 0) ;
    lexaddkey ("Score_bounds", &_Score_bounds, 0) ;
    lexaddkey ("Score_by_offset", &_Score_by_offset, 0) ;
    lexaddkey ("Score_by_width", &_Score_by_width, 0) ;
    lexaddkey ("Score_by_histogram", &_Score_by_histogram, 0) ;
    lexaddkey ("BlastN", &_BlastN, 0) ;
    lexaddkey ("Percent", &_Percent, 0) ;
    lexaddkey ("Right_priority", &_Right_priority, 0) ;
    lexaddkey ("Web_right_priority", &_Web_Right_priority, 0) ;

    isDone = TRUE  ;
    }

  return ;
}



/* methodAdd ()
 * General routine to handle cached information on methods in methodInfo
 * Will always ensure METHOD *method(method) exists
 *
 * Return values:
 * 0 if previously DONE
 * 1 if new and not in database (set to defaults)
 * 2 if new and read from database
 */

int methodAdd (KEY view, KEY method)
{
  OBJ obj ;
  char *text ;
  METHOD *meth ;


  methodInitialise() ;				  /* Make sure we are initialised. */


  method = lexAliasOf(method); /* use canonical key */

  if (class(method) != _VMethod)
    messcrash ("Bad class of argument to addMethod()") ;

  if (!methodInfo)
    methodInfo = arrayCreate (lexMax(_VMethod), METHOD*) ;

  if (!(meth = array (methodInfo, KEYKEY(method), METHOD*)))
    meth = array (methodInfo, KEYKEY(method), METHOD*) = (METHOD*) messalloc (sizeof(METHOD)) ;

  if (meth->flags & METHOD_DONE)
    return 0 ;

  if (!(obj = bsCreate (method))) 
    return 1 ;

  meth->flags &= METHOD_CALCULATED ;	/* reset to defaults */
  meth->flags |= METHOD_DONE ;
  meth->col = 0 ;
  meth->width = 0 ;
  meth->symbol = 0 ;
  meth->priority = 0 ;
  meth->min = meth->max = 0 ;

  if (bsGetKeyTags (obj, _Colour, &meth->col))
    meth->col -= _WHITE ;
  bsGetData (obj, _Width, _Float, &meth->width) ;
  bsGetData (obj, _Right_priority, _Float, &meth->priority) ;
  if (view)
    bsGetData (obj, _Web_Right_priority, _Float, &meth->priority) ;
  if (bsGetData (obj, _Symbol, _Text, &text))
    meth->symbol = *text ;
  if (bsFindTag (obj, _Frame_sensitive))
    meth->flags |= METHOD_FRAME_SENSITIVE ;
  if (bsFindTag (obj, _Strand_sensitive))
    meth->flags |= METHOD_STRAND_SENSITIVE ;
  if (bsFindTag (obj, _Show_up_strand))
    meth->flags |= METHOD_SHOW_UP_STRAND ;
  if (bsFindTag (obj, _BlastN))
    meth->flags |= METHOD_BLASTN ;
  if (bsFindTag (obj, _Blixem_X))
    meth->flags |= METHOD_BLIXEM_X ;
  if (bsFindTag (obj, _Blixem_N))
    meth->flags |= METHOD_BLIXEM_N ;
  if (bsFindTag (obj, _Blixem_P))
    meth->flags |= METHOD_BLIXEM_P ;
  if (bsFindTag (obj, _Belvu))
    meth->flags |= METHOD_BELVU ;
  if (bsFindTag (obj, _Bumpable))
    meth->flags |= METHOD_BUMPABLE ;
  if (bsFindTag (obj, _EMBL_dump))
    meth->flags |= METHOD_EMBL_DUMP ;
  if (bsFindTag (obj, _Percent))
    { meth->flags |= METHOD_PERCENT ;
      meth->min = 25 ;		/* so actual display is linear */
      meth->max = 100 ;		/* can override explicitly */
    }
  if (bsGetData (obj, _Score_bounds, _Float, &meth->min) &&
      bsGetData (obj, _bsRight, _Float, &meth->max))
    {
      if (bsFindTag (obj, _Score_by_offset))
	meth->flags |= METHOD_SCORE_BY_OFFSET ;
      else if (bsFindTag (obj, _Score_by_width))
	meth->flags |= METHOD_SCORE_BY_WIDTH ;
      else if (bsFindTag (obj, _Score_by_histogram))
	{ meth->flags |= METHOD_SCORE_BY_HIST ;
        meth->histBase = meth->min ;
	bsGetData (obj, _bsRight, _Float, &meth->histBase) ;
	}
    }


  bsDestroy (obj) ;
  return 2 ;
}

METHOD *method (KEY view, KEY key)
{
  METHOD *meth ;

  methodInitialise() ;				  /* Make sure we are initialised. */

  key = lexAliasOf(key) ;
  if (!key) return 0 ;
  if (!methodInfo ||
      KEYKEY(key) >= arrayMax(methodInfo) ||
      !(meth = arr(methodInfo, KEYKEY(key), METHOD *)))
    { methodAdd (view, key) ;
      meth = arr(methodInfo, KEYKEY(key), METHOD *) ;
    }
  return meth ;
}


void methodSet (char *name, int col, unsigned int flags, float priority,
		float width, char symbol, float min, float max)
{
  KEY key ;
  METHOD *meth ;

  methodInitialise() ;				  /* Make sure we are initialised. */

  lexaddkey (name, &key, _VMethod) ;
  if (methodAdd (0, key) != 1)			  /* only set here if not in database */
    return ;

  meth = arr (methodInfo, KEYKEY(key), METHOD*) ; /* made in methodAdd */

  meth->col = col ;

  /* I have no idea why this is set, it can probably be removed, I can't see */
  /* how it is used.                                                         */
  meth->flags = flags | METHOD_INTERNAL ;

  meth->priority = priority ;
  meth->width = width ;
  meth->symbol = symbol ;
  meth->min = min ;
  meth->max = max ;
}

/**************************** end of file *******************************/
 
 
 
 
 
 
