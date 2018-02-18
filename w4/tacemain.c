/*  File: tacemain.c
 *  Author: Jean Thierry-Mieg (mieg@crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: main fucntion for a text-only acedb interface
 *              
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  2 16:22 1998 (edgrif)
 * * Dec  2 16:22 1998 (edgrif): Corrected decl. of main, added code
 *              to record build time of this module.
 * * Jun  7 16:49 1996 (srk)
 * * May 15 12:37 1992 (mieg): This file used to be called querymain
 * *                  I import an adapt parseInput from DKFZ
 * Created: Fri May 15 12:27:28 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: tacemain.c,v 1.11 2012/04/05 03:06:46 mieg Exp $ */

#include "acedb.h"
#include "session.h"
#include "command.h"
#include "banner.h"
#include "version.h"

/**************************************************/
/******* ACEDB non graphic user interface  ********/
/**************************************************/

/* Defines a routine to return the compile date of this file.                */
UT_MAKE_GETCOMPILEDATEROUTINE()


int main (int argc, const char **argv)
     /* main function used for tace and tacembly */
{
  int level ;
  BOOL myIsInteractive = TRUE ;
  int ii,  outDelay = 0 ;
  const char *db = NULL;
  const char *prompt = 0, *buff = 0 ;

  messErrorInit(argv[0]) ;			  /* Record program name for crash messages. */

  if (getCmdLineOption (&argc, argv, "-noprompt", 0))
    myIsInteractive = FALSE ;
  if (getCmdLineOption (&argc, argv, "-no_prompt", 0))
    myIsInteractive = FALSE ;
  if (getCmdLineOption (&argc, argv, "-prompt", &prompt))
    set_user_prompt (prompt) ;
  if (getCmdLineOption (&argc, argv, "-out_delay", &buff) &&
      (ii = atoi (buff)) &&
      ii > 0
      )
    outDelay = ii ;
  if (!getCmdLineOption (&argc, argv, "-db", &db) && argc>=2)
    db = argv[1] ;

  setbuf (stdout, NULL) ;
  setbuf (stderr, NULL) ;
  aceInit( db );

#ifdef ACEMBLY
  bannerWrite (bannerMainStrings ("tacembly", FALSE, TRUE)) ;
#else
  bannerWrite (bannerMainStrings ("tace", FALSE, FALSE)) ;
#endif
  bannerWrite (bannerDataInfoStrings ()); /* must come after aceInit() */
    
  level = freesetfile(stdin,"") ;  
#ifdef ACEMBLY
  commandExecute (level, FALSE, myIsInteractive, stdout, 0, 7, 0) ;
#else
  commandExecute (level, FALSE, myIsInteractive, stdout, 0, 3, 0) ;
#endif

  /* finish the sesion and clean up */
  if (isWriteAccess() && 
      (!myIsInteractive || messQuery ("You did not save your work, should I ?")))
    aceQuit (TRUE);
  else
    aceQuit (FALSE);

  {
    int mx ;
    messAllocMaxStatus (&mx) ; 
    fprintf (stderr, "// %s done: max memory %d Mb\n", timeShowNow(), mx) ;
  }

  printf ("\n// A bientot \n\n") ;
  if (outDelay > 0 && outDelay < 120)
    sleep (outDelay) ; /* hack around nfs bugs, i hope this way that the next acedb sees a synch disk */
  return(EXIT_SUCCESS) ;
}

/**************************************************/
/****************** eof ***************************/
