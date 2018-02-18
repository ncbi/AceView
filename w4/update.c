/*  File: update.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Parses the official acedb updates.
 * Exported functions:
 *              updateData
 *              updateDoAction
 * HISTORY:
 * Last edited: Jan  8 11:30 1999 (fw)
 * * Jan 20 10:16 1992 (mieg): interactive access to alias, 
   should be done differently really.
 * Created: Fri Oct 25 13:35:14 1991 (mieg)
 *-------------------------------------------------------------------
 */

/*$Id: update.c,v 1.3 2006/12/16 05:04:34 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "main.h"		/* pickDraw */
#include "session.h"
#include "parse.h"
#include "model.h"
#include "lex.h"		/* for lexAlphaMakeAll */
#include "bs.h"			/* for global BOOL flags */

#include "vmap.h"		/* for vMapMakeAll() */

/*************************************************************/

static int lineBox ;

/*************************************************************/
/*************************************************************/

static Graph updateGraph = 0;
static int line = 3 ;
static void upMessage(char *cp)
{
  Graph old = graphActive() ;
  if(graphActivate(updateGraph))
    { graphText (cp, 1, line++) ;
      graphTextBounds (80,line+2) ;
      graphBoxDraw (0,-1,-1) ;
    }
  graphActivate(old) ;
}
    
/******************************************/

static void updateAlign (void)
{
  extern void gMapMakeAll(void) ;
  extern void pMapMakeAll(void) ;
  extern void cMapMakeAll(void) ;

  upMessage ("Making physical maps - a few minutes") ;
  pMapMakeAll () ;

  upMessage ("Making genetic maps - a minute or two") ;
  vMapMakeAll () ;
  gMapMakeAll () ;

  upMessage ("Making physical chromo maps - a minute or two") ;
  cMapMakeAll () ;

  upMessage ("Sorting the lexiques - a minute or two") ;
  lexAlphaMakeAll () ;

  upMessage ("I save to disk") ;
  sessionClose (TRUE);

  pickDraw () ;			/* updates title bar */

  upMessage ("Update complete") ;
  if (sessionDbName())
    messout ("We are now using data release %s %d-%d",
	     sessionDbName(),
	     thisSession.mainDataRelease, thisSession.subDataRelease) ;
  else
    messout ("We are now using data release %d-%d",
	     thisSession.mainDataRelease, thisSession.subDataRelease) ;
}

/******************************************/

static BOOL updateCheckFile (FILE *fil)
{ 
  int mainRelNum = -1, subRelNum = -1 ;

  if (fscanf (fil, "// acedb update %d %d\n",
	      &mainRelNum, &subRelNum) != 2)
    { messout ("The first line of the update does not have the "
	       "correct format.  Abandoning update.") ;
      return FALSE ;
    }

  if (mainRelNum != thisSession.mainDataRelease ||
      subRelNum != thisSession.subDataRelease+1)
    { 
      messout  ("Update %d.%d does not match the current state",
		mainRelNum, subRelNum) ;
      return FALSE ;
    }
  return TRUE ;
}

/******************************************/

static BOOL updateParseFile (FILE *fil)
{
  BOOL b ;

  messStatus (messprintf ("Adding update %d", 
			    thisSession.subDataRelease + 1)) ;
  upMessage ("Update started. This may take a few minutes please wait.") ;

  XREF_DISABLED_FOR_UPDATE = TRUE ;
  b = parseFile (fil, lineBox, 0) ;  /* Does close fil */
  XREF_DISABLED_FOR_UPDATE = FALSE ;

   /* mieg: DO NOT:    filremove (fileName, "") ;  remove file after successful update
    * this would be an ugly side effect.
    */
  return b ;
}

/******************************************/

void updateDoAction (BOOL doAll)
{
  VoidRoutine previousCrashRoutine ;
  extern VoidRoutine messcrashroutine ;
  extern void simpleCrashRoutine (void);
  FILE *fil = 0 ;
  char fileName[256] , * fn;
  BOOL doAlign = FALSE ;
  char timeBuf[25] ;

#if !defined(MACINTOSH)
  if (getenv ("ACEDB_DATA"))
    strcpy (fileName, getenv("ACEDB_DATA")) ;
  else
#endif
  if (!filName ("rawdata", 0, "r"))
    { messout ("Sorry, can't find rawdata directory") ;
      return ;
    }
  else
    strcpy (fileName, filName ("rawdata", 0, "r")) ;

  if (isWriteAccess())
    { messout("You still have write access\n"
	      "First, save your earlier work") ;
      return ;
    }

  if (!sessionGainWriteAccess())
    { 
      messout ("Sorry, you cannot gain write access\n"
	       "Update impossible.") ;
      return ;
    }

  fn = fileName + strlen(fileName) ;

  previousCrashRoutine = messcrashroutine ;
  messcrashroutine = simpleCrashRoutine ; /* don't allow crash recovery */

  upMessage ("Reading in model file") ;
  if (!getNewModels ())		/* needs write access */
    messExit ("There were errors while reading the models during the update") ;

  while(TRUE)
    { if (!isWriteAccess())
	sessionGainWriteAccess() ;
      if (sessionDbName())
	strcpy (fn, messprintf ("/update.%s.%d-%d",
				sessionDbName(),
				thisSession.mainDataRelease,
				thisSession.subDataRelease + 1)) ;
      else
	strcpy (fn, messprintf ("/update.%d-%d",
			      thisSession.mainDataRelease,
			      thisSession.subDataRelease + 1)) ;
      upMessage (messprintf("Looking for file %s",fileName)) ;
      upMessage (timeShow(timeNow(), timeBuf, 25)) ;
      if (!filName (fileName, "", "r"))
	{ if (!doAlign)
	    messout ("Could not find file %s",fileName) ;
	  else
	    upMessage (messprintf ("Could not find file %s",
				   fileName)) ;
	  break ;
	}
      if (!(fil = filopen (fileName, "", "r"))) /* should not fail */
	{ messcrash ("code bug:filName/filopen inconsistency - please report") ;
	  break ;
	}
      
      messStatus("Updating, please wait.") ;
      if (!updateCheckFile (fil))
	{ filclose(fil) ;
	  break ;
	}

      if (!updateParseFile (fil))  /* does close fil */
	messcrash ("There were errors while parsing the "
		   "official data: I exit.") ;

      thisSession.subDataRelease++ ;
      doAlign = TRUE ;
      sessionClose(TRUE) ;
      if (!doAll && !messQuery
	  (messprintf ("We are now using data release %d-%d\n%s",
		       thisSession.mainDataRelease,
		       thisSession.subDataRelease,
		       "Do you want to read the next update file ? ")))
	break ;
    }
  upMessage (timeShow(timeNow(), timeBuf, 25)) ;
  if (doAlign)
    { if (!isWriteAccess())
	sessionGainWriteAccess() ;
      updateAlign () ;
    }
  sessionClose(TRUE) ;  
  upMessage (timeShow(timeNow(), timeBuf, 25)) ;
  messcrashroutine = previousCrashRoutine; /* reestablish crash recovery */
} /* updateDoAction */

/******************************************/

static void updateAction(void)
{ updateDoAction(FALSE) ;
}

static void allUpdateAction(void)
{ updateDoAction(TRUE) ;
}

/******************************************/

static void updateDestroy (void)
{
  updateGraph = 0 ;
}

/*********************************************/
/************* public routine ****************/

static MENUOPT updateMenu[]=
{
  { graphDestroy,	"Quit" },
  { help,		"Help" },
  { updateAction,	"Next Update" },
  { allUpdateAction,	"All Updates" },
  { 0,0 }
};

/********************************************/

void updateData(void)
{
  if (!graphActivate(updateGraph))
    { updateGraph = displayCreate(DtUpdate) ;
      graphMenu (updateMenu);
      graphTextBounds (80,100) ;
      graphRegister(DESTROY,updateDestroy) ;
      graphColor (BLACK) ;
      graphButtons (updateMenu, 1, 1, 50) ;
      graphText ("Line:", 2, 2.5) ;
      lineBox = graphBoxStart () ;
      *parseLineText = 0 ;
      graphTextPtr (parseLineText, 8, 2.5, 60) ;
      graphBoxEnd () ;
      line = 4 ;
      if (sessionDbName())
	upMessage (messprintf ("We are now using data release %s %d-%d",
			     sessionDbName(),
			     thisSession.mainDataRelease,
			     thisSession.subDataRelease)) ;
      else
	upMessage (messprintf ("We are now using data release %d-%d", /*  */
			     thisSession.mainDataRelease,
			     thisSession.subDataRelease)) ;
      graphRedraw () ;
    }
  else
    graphPop() ;
}
 
/*********************************************/
/*********************************************/

 
