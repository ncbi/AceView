/*  File: alignmentdisp.c
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: graphical alignment display functions
 * Exported functions:
 *              alignDisplay
 * HISTORY:
 * Last edited: Nov 19 15:18 1998 (fw)
 * Created: Thu Nov 19 11:41:03 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: alignmentdisp.c,v 1.1.1.1 2002/07/19 20:22:58 sienkiew Exp $ */

#include "acedb.h"
#include "graph.h"
#include "display.h"

#include "alignment_.h"

/************************************************************/

#ifdef JUNK

static int ALIGN_MAGIC ;

static ALIGNMENT *alignFromGraph (void)
{ 
  ALIGNMENT *al ;
  if (!(graphAssFind (&ALIGN_MAGIC, &al)))
    messcrash ("Failed to find alignment from graph") ;
  return al ;
}



static void alignTreeShow (void)
{ 
  ALIGNMENT *al = alignFromGraph () ;

  display (al->key, 0, TREE) ;
}

static void alignDisplayDump (void)
{ 
  ALIGNMENT *al = alignFromGraph () ;
  FILE *fil = filqueryopen (0, 0, "aln", "w", "Alignment file name") ;

  if (!fil)
    return ;
  alignDumpKey (al->key, fil) ;
  filclose (fil) ;
}

static void alignDisplayKeySetDump (void)
{ 
  FILE *fil ;
  KEYSET ks ;

  if (!keySetActive (&ks, 0))
    { messout ("Can't find active keyset") ;
      return ;
    }
  if (!(fil = filqueryopen (0, 0, "aln", "w", "Alignment file name")))
    return ;
  alignDumpKeySet (ks, fil) ;
  filclose (fil) ;
}
 


static MENUOPT menu[] = {
  {graphDestroy, "Quit"},
  {graphPrint, "Print"},
  {alignTreeShow, "Show as tree"},
  {alignDisplayDump, "Dump"},
  {alignDisplayKeySetDump, "Dump KeySet"},
   {0, 0}
} ;


BOOL alignDisplay_old (KEY key, KEY from, BOOL isOldGraph)
{ 
  ALIGNMENT *al ;
  ALIGN_COMP *c ;
  int i, xmax = 0, ymax = 0, x, y, maxlen = 0 ;
  int line = 0 ;

  if (!(al = alignGet (key)))
    { display (key, key, TREE) ;
      return FALSE ;
    }

  if (isOldGraph)
    graphClear () ;
  else
    { displayCreate (DtAlignment) ;
      graphMenu (menu) ;
    }

  graphAssociate (&ALIGN_MAGIC, al) ;

				/* now draw things */

  graphText (name(key), 0, line) ;
  line += 2 ;

  c = arrp(al->comp, 0, ALIGN_COMP) ;
  for (i = arrayMax(al->comp) ; i-- ; ++c)
    { x = strlen (messprintf ("%d - %d", c->start, c->end)) ;
      if (x > xmax)
	xmax = x ;
      y = strlen (name (c->key)) ;
      if (y > ymax)
	ymax = y ;
    }
  ymax += 2 ;
  xmax += ymax + 2 ;

  c = arrp(al->comp, 0, ALIGN_COMP) ;
  for (i = arrayMax(al->comp) ; i-- ; ++c)
    { graphText (name(c->key), 0, line) ;
      graphText (messprintf ("%d - %d", c->start, c->end), ymax, line) ;
      if (c->seq)
	{ graphText (c->seq, xmax, line) ;
	  if (strlen (c->seq) > maxlen)
	    maxlen = strlen (c->seq) ;
	}
      ++line ;
    }
  
  graphTextBounds (xmax+maxlen, line) ;
  graphRedraw () ;
  return TRUE ;
}
#endif /* JUNK */

BOOL alignDisplay (KEY key, KEY from, BOOL isOldGraph)
{ 
#if !defined(NO_POPEN)

  ALIGNMENT *al ;
  ALIGN_COMP *c ;
  int i, x, xmax = 0 ;
  FILE *pipe ;
  char *domain, *cp, *script ;
  char belvu_script[] = "wscripts/belvu_script" ;

  extern char *findCommand (char *command, char **retp);
  extern Stack BelvuMatchStack;

  if (!(al = alignGet (key)))
    { display (key, key, TREE) ;
      return FALSE ;
    }

  /* Open pipe to belvu */
  /* Can't call callScriptPipe since we want write access and buildCommand 
     is not an exported function in the graph library */
  if ((cp = filName (belvu_script, 0, "x")))
    script = cp ;
  else
    script = belvu_script ;
  printf("Calling \"%s %s\"", script, name(key));
  fflush(stdout);
  pipe = (FILE *)popen(messprintf("%s %s", script, name(key)), "w");

  /* Find longest name string */
  c = arrp(al->comp, 0, ALIGN_COMP) ;
  for (i = arrayMax(al->comp) ; i-- ; ++c)
    { x = strlen (messprintf ("%s/%d-%d", name(c->key), c->start, c->end)) ;
      if (x > xmax)
	xmax = x ;
    }

  domain = messalloc(xmax+1);

  c = arrp(al->comp, 0, ALIGN_COMP) ;
  for (i = arrayMax(al->comp) ; i-- ; ++c)
    { sprintf (domain, "%s/%d-%d", name(c->key), c->start, c->end) ;
      fprintf (pipe, "%-*s ", xmax, domain) ;
      if (c->seq)
	fprintf (pipe, "%s", c->seq) ;
      fprintf (pipe, "\n") ;
    }
  
  messfree(domain);

  /* MatchFooter */
  if (stackExists(BelvuMatchStack)) {
    char *cp;

    fprintf(pipe, "# matchFooter\n");

    stackCursor(BelvuMatchStack, 0) ;
    cp = stackNextText(BelvuMatchStack);
    printf(" with match %s\n", cp);

    fprintf(pipe, "%s\n", cp);
    fprintf(pipe, "%s\n", stackNextText(BelvuMatchStack));
    while ((cp = stackNextText (BelvuMatchStack)))
      fprintf(pipe, "%s ", cp);
    stackDestroy(BelvuMatchStack);
  }
  else printf("\n");
  fflush(stdout);

  fprintf (pipe, "%c\n", EOF) ; /* To close the pipe, sigh */
  fflush(pipe);
  /* pclose is no good - waits till belvu is finished. */

  return TRUE ;
#else
  return FALSE ;
#endif  /* not NO_POPEN */
} /* alignDisplay */

/*************************** eof ****************************/
