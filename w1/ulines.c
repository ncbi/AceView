/*  File: arraysub.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *              routines to break a long text in lines
 *              separated out from arraysubs.c becuase not-rentrant
 *               (part of regular.h - the header for libfree.a)
 * Exported functions:
 *              See Header file: array.h (includes lots of macros)
 * HISTORY:
 * Last edited: Dec  4 11:12 1998 (fw)
 * * Nov  1 16:11 1996 (srk)
 *		-	MEM_DEBUG code clean-up 
 *                      (some loose ends crept in from WIN32)
 *		-	int (*order)(void*,void*) prototypes used 
 *                                            uniformly throughout
 * * Jun  5 00:48 1996 (rd)
 * * May  2 22:33 1995 (mieg): killed the arrayReport at 20000
      otherwise it swamps 50 Mbytes of RAM
 * * Jan 21 16:25 1992 (mieg): messcrash in uArrayCreate(size<=0)
 * * Dec 17 11:40 1991 (mieg): stackTokeniseTextOn() tokeniser.
 * * Dec 12 15:45 1991 (mieg): Stack magic and stackNextText
 * Created: Thu Dec 12 15:43:25 1989 (mieg)
 *-------------------------------------------------------------------
 */

#include "regular.h"
#include "array.h"

/****************** routines to set text into lines ***********************/

static char *linesText ;
static Array lines ;
static Array textcopy ;
static int kLine, popLine ;

/**********/

int uLinesText (char *text, int width)
{
  char *cp,*bp ;
  int i ;
  int nlines = 0 ;
  int length = strlen (text) ;
  int safe = length + 2*(length/(width > 0 ? width : 1) + 1) ; /* mieg: avoid zero divide */
  static int isFirst = TRUE ;

  if (isFirst)
    { isFirst = FALSE ;
      lines = arrayCreate(16,char*) ;
      textcopy = arrayCreate(128,char) ;
    }

  linesText = text ;
  array(textcopy,safe,char) = 0 ;   /* ensures textcopy is long enough */

  if (!*text)
    { nlines = popLine = kLine = 0 ;
      array(lines,0,char*) = 0 ;
      return 0 ;
    }

  cp = textcopy->base ;
  nlines = 0 ;
  for (bp = text ; *bp ; ++bp)
    { array(lines,nlines++,char*) = cp ;
      for (i = 0 ; (*cp = *bp) && *cp != '\n' ; ++i, ++cp, ++bp)
        if (i == width)		/* back up to last space */
          { while (i--)
	      { --bp ; --cp ;
		if (*cp == ' ' || *cp == ',' || *cp == ';')
		  goto eol ;
	      }
	    cp += width ;	/* no coma or spaces in whole line ! */
	    bp += width ;
	    break ;
	  }
eol:  if (!*cp)
        break ;
      if (*cp != '\n')
	++cp ;
      *cp++ = 0 ;
    }
  kLine = 0 ;			/* reset for uNextLine() */
  popLine = nlines ;
  array(lines,nlines,char*) = 0 ; /* 0 terminate */

  return nlines ;
 }

char *uNextLine (char *text)
 {
   if (text != linesText)
     messout ("Warning : uNextLine being called with bad context") ;
   return array(lines,kLine++,char*) ;
 }

char *uPopLine (char *text)
 {
   if (text != linesText)
     messout ("Warning : uPopLine being called with bad context") ;
   if (popLine)
     return array(lines,--popLine,char*) ;
   else return 0;
 }

char **uBrokenLines (char *text, int width)
{
  uLinesText (text,width) ;
  return (char**)lines->base ;
}

char *uBrokenText (char *text, int width)
{
  char *cp ;

  uLinesText (text,width) ;
  uNextLine (text) ;
  while ((cp = uNextLine (text)))
    *(cp-1) = '\n' ;
  return arrp(textcopy,0,char) ;
}

