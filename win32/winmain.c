/*  File: winmain.c
 *  Author: Richard Bruskiewich	(rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: WIN32 specific variants of ACEDB initialisation,
 *              main window order taking and termination routines
 * Exported functions:
 * HISTORY:
 * Last edited: Jun 9 14:54 1996 (rbrusk): see winace.cpp for WIN32 wormClose()
 *  -   wormClose() here renamed to DoIQuit() but kept here for convenience 
 * * Jun 30 14:23 1995 (rbrusk):
 *      -   Extracted worminit() and wormclose() from xacemain, modified for WIN32
 * Created: Tue Oct 22 12:51:04 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: winmain.c,v 1.3 2011/10/09 03:13:14 mieg Exp $ */

#include "win32.h"

#include "acedb.h"
#include "graph.h"
#include "session.h"
#include "../wh/menu.h"
#include "systags.h"
#include "display.h"
#include "graph_.h"

extern VoidRoutine messcrashroutine ; /* in messubs.c */
extern MENUOPT quovadis[] ;

/******************* INITIALISATION ***************************/
/**************************************************************/

MENU mainMenu ;		/* global so that session, wcs can change */
Graph mainGraph ;

extern void winCrash (void);	/* see below */

// extern void*openMonitor(const char *) ;
// extern BOOL closeMonitor(void *) ;
// static void *wormProgress ;

void wormInit (int argc, char **argv )
{
  /* The order of the initialisations is crucial */

  // Progress Monitor to feedback ACEDB session
  // completion activities to user
  // wormProgress = openMonitor("Initializing WinAce") ;

  freeinit () ;			/* must come before graphInit */
  freespecial ("\n\t\"/%\\@") ;	/* No subshells ($ removed) */
  graphInit (&argc, argv) ;	/* before calls to messcrash() */
  sessionInit (&argc, argv) ;

 /* must come after successful session initialisation */
  messcrashroutine = winCrash ;  /* see WIN32\mainfrm.cpp */ 

  // closeMonitor(wormProgress) ; // ending monitoring
}

/*************************************************************/
/****************** MAIN LOOP **********************/

extern Graph pickCreate (BOOL chooseSession) ;

void askOrders(void)
{
	mainGraph = pickCreate (FALSE) ; 

	mainMenu = menuInitialise ("Main menu", (MENUSPEC*)quovadis) ;
#ifdef WCS
	menuSuppress (mainMenu, "WCS annotate") ;
#endif
	if (sessionCheckUserName (FALSE))
		menuSetFlags (menuItem (mainMenu, "Save"), MENUFLAG_DISABLED) ;
	else    
	{
		menuSuppress (mainMenu, "Write Access") ;
		menuSuppress (mainMenu, "Add-Alias-Rename") ;
	}
	graphNewMenu (mainMenu) ;

	// No graphLoop(FALSE) here for WIN32, because WIN32 has
	// an implicit event loop built into the top level of the program
}

BOOL DoIQuit (void)
{
	/* Returns a value to WinAceApp for consideration */
	if (  thisSession.session != 1 &&
		 !messQuery ("Do you really want to quit WinAce?") )
		return FALSE;
	else
		return TRUE;
}

/*************************************************/
/****************** TERMINATION ******************/

void wormClose2 (void)	/* distinct from wormClose since called by messcrash */
{
  graphCleanUp () ;  /* kills all open windows, forces displayed objects to cache */

  graphWaitCursor(TRUE) ;

  if (isWriteAccess() && !messQuery ("You did not save your work, should I ?"))
    sessionClose (FALSE) ;
  else
    sessionClose (TRUE) ;

#ifdef XCLIENT
  xClientClose() ;
#endif

  /* better kill the wait cursor here
     since main graph is dying now */
  graphWaitCursor(FALSE) ;

  graphFinish () ;
  filtmpcleanup () ;		/* deletes any temporary files made */

#if defined(MEM_DEBUG)
  handleCleanUp() ;			/* Do this last, since it frees ALL memory */
#endif

  /* just return to caller, for further cleanup before exit() */
}


/******************  ACEDB  ***********************/
/**************************************************/

/*
    no main() routine... entry point handled within 
	MFC constructed WinApp in WinAce.cpp
*/

/**************** End of File ******************/
 
 
