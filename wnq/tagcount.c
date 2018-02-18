/*  File: tagcount.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: tagcount.c,v 1.3 2010/06/15 17:55:03 mieg Exp $
 * Description: little program to count usages of tags in each class
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  2 15:54 1998 (edgrif)
 * * Dec  2 15:53 1998 (edgrif): Make declaration of main correct, add
 *              code to record build time of this module.
 * Created: Mon Sep 29 13:26:50 1997 (rd)
 *-------------------------------------------------------------------
 */

#define DEFINE_OBJ
typedef struct sobj *OBJ ;
#include "acedb.h"
#include "bs_.h"
#include "version.h"

static Array a ;

static void addInfo (BS bs)
{
  if (!class(bs->key) && bs->key > _LastN)
    ++array(a,bs->key,int) ;
  if (bs->right)
    addInfo (bs->right) ;
  if (bs->down)
    addInfo (bs->down) ;
}

static void printUsage (void)
{
  fprintf (stdout, 
	   "Usage: tagcount [options] <acedbdir>\n"
	   "   Options :\n"
	   "       -u   only user classes\n"
	   "       -s   only system classes\n"
	   "       -a   all classes\n"
	   "   if no option is given, user classes are shown.\n") ;
}


/* Defines a routine to return the compile date of this file.                */
UT_MAKE_GETCOMPILEDATEROUTINE()


/************************************************************/

int main (int argc, char **argv)
{
  int c, i, n ;
  KEY k ;
  OBJ o ;
  BOOL IsSysClassShown = FALSE;
  BOOL IsUserClassShown = TRUE;

  --argc ; ++argv ;
  if (!argc)
    { 
      printUsage ();
      return (EXIT_FAILURE);
    }

  if (argv[0][0] == '-')
    {
      /* an option is specified */
      if (argv[0][1] == 'u')
	{
	  IsSysClassShown = FALSE;
	  IsUserClassShown = TRUE;
	}
      else if (argv[0][1] == 's')
	{
	  IsSysClassShown = TRUE;
	  IsUserClassShown = FALSE;
	}
      else if (argv[0][1] == 'a')
	{
	  IsSysClassShown = TRUE;
	  IsUserClassShown = TRUE;
	}
      else
	{
	  fprintf (stderr, "tagcount: illegal option %s\n", *argv);
	  printUsage ();
	  return (EXIT_FAILURE);
	}
      --argc ; ++argv ;
    }

  if (!argc)
    { 
      printUsage ();
      return (EXIT_FAILURE);
    }

  aceInit (*argv) ;

  a = arrayCreate (5000, int) ;

  for (c = 0 ; c < 256 ; ++c)
    { 
      if (pickList[c].type != 'B')
	continue;
      if ((IsSysClassShown && pickList[c].protected) ||
	  (IsUserClassShown && !pickList[c].protected))
	{
	  a = arrayReCreate (a, 5000, int) ;
	  k = 0 ;
	  n = 0 ;
	  while (lexNext(c,&k))
	    if ((o = bsCreate (k)))
	      { 
		addInfo (o->root) ;
		bsDestroy (o) ;
		++n ;
	      }
	  printf ("%s\tOBJECT_COUNT\t%d\n", pickClass2Word(c), n) ;
	  for (i = _LastN ; i < arrayMax (a) ; ++i)
	    if (arr(a,i,int))
	      printf ("%s\t%s\t%d\n", 
		      pickClass2Word(c), name(i), arr(a,i,int)) ;
	}
    }
  return (EXIT_SUCCESS) ;
}

/******************* end of file **********************/
 
