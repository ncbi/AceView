/*  File: sessiondisp.c
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: interactive session control
 * Exported functions:
 *              void sessionControl(void)
 * HISTORY:
 * Last edited: Mar 16 16:56 1999 (fw)
 * * Dec  3 14:45 1998 (edgrif): Change calls to new interface to aceversion.
 * * Nov 23 16:57 1998 (fw): functions moved here from session.c
 * Created: Mon Nov 23 16:57:22 1998 (fw)
 *-------------------------------------------------------------------
 */

/************************************************************/

#include "acedb.h"

#include "display.h"
#include "main.h"
#include "lex.h"
#include "lex_sess_.h"
#include "systags.h"		/* for _This_session tag */
#include "bs.h"
#include "pick.h"
#include "disk.h"
#include "disk_.h"   /* defines BLOCKHEADER */
#include "session.h"
#include "session_.h"		/* uses BLOCKHEADER */

/************************************************************/
#define SESS_COL_READLOCKED ORANGE
#define SESS_COL_PERMANENT  GREEN
#define SESS_COL_SELECTED   RED
#define SESS_COL_KEPT_ALIVE YELLOW
#define SESS_COL_DESTROYED  LIGHTGRAY
#define SESS_COL_CURRENT    PALEBLUE

static void sessDispDraw(void);
static void sessDispDestroy(void);
static void sessDispSetStartSession(void);
static void sessDispFixSession(void);
static void sessDispUnFixSession(void);
static void sessDispDestroyLastSession(void);
static void sessDispToggleDeadShown (void);
static void sessionPick(int box);
/************************************************************/

static Graph sessionSelectionGraph = 0 ;
static KEY sessionChosen = 0 ;
static int selectedSessionBox = 0 ;
static Array box2sess = 0;	/* arrayp of ST */
static BOOL IsDeadSessionShown = FALSE;

static MENUOPT sessionMenu[] = {
  { graphDestroy,		"Quit"},
  { help,			"Help"},
  { menuSpacer,			""},
  { sessDispSetStartSession,	"Start Session"},
  { sessDispFixSession,		"Fix Session"},
  { sessDispUnFixSession,	"Unfix Session"},
  { sessDispDestroyLastSession,	"Destroy Last Session"}, 
  { sessDispToggleDeadShown,		"Show/Hide dead sessions"},
  { 0,0 }
} ;
/************************************************************/



/*************************************************************/
/**************** Interactive session control ****************/
/*************************************************************/

void sessionControl(void)
{ 
  sessionChosen = 0;

  if (graphActivate(sessionSelectionGraph))
    graphPop() ;
  else
    sessionSelectionGraph =  displayCreate(DtSession) ;
  if (! sessionSelectionGraph) /* happens if display.wrm missed the config */
    sessionSelectionGraph = graphCreate (TEXT_FULL_SCROLL, "Session Control", .41, .21, .64, .64) ;
  graphMenu(sessionMenu) ;

  graphRegister(DESTROY, sessDispDestroy) ;
  graphRegister(RESIZE, sessDispDraw) ;
  graphRegister(PICK,  sessionPick) ;

  sessionChangeRegister (sessDispDraw);

  sessDispDraw() ;
} /* sessionControl */

/*************************************************************/

static void sessDispDestroy(void)
{ 
  sessionChangeRegister (0);

  sessionSelectionGraph =  0 ;
  arrayDestroy(box2sess); 
  box2sess = 0;

  return;
} /* sessDispDestroy */

/*************************************************************/

static void sessDispSetStartSession(void)
{ 
  if(!sessionChosen)
   { messout("First pick a live session") ;
     return ;
   }

  if (isWriteAccess())
    { messout("You still have write access\n"
	      "First save your previous work") ;
      return ;
    }
  
  if (writeAccessPossible() && 
      !messQuery ("By starting again in a previous session "
		  "you will lose the ability to re-gain write access, "
		  "even by going back to the current session.\n"
		  "Continue ?"))
    return;
		  


  pickPopMain() ;
  graphCleanUp () ;		/* kill all but main graph */

  sessionDoClose () ;		/* To empty all caches */
  sessionForbidWriteAccess();	/* can never gain write access again */
  
  lexClear() ;

  sessionStart(sessionChosen);	/* link up the current session to
				   the selected session */
  pickInit () ;                 /* needed, used to be part of sessionInit */
  sessionInitialize ();
  bIndexInit(BINDEX_AUTO) ;           /* will do only from disk */
  sessionControl();		/* re-display session control
				   it was killed by graphCleanUp */

  return;
} /* sessDispSetStartSession */

/************************************************************/  

static int stBoxColour (ST *st)
{
  int box_colour;
  
  /* colour the box according to the status of the session */
  if (st->isReadlocked)
    box_colour = SESS_COL_READLOCKED;
  else if(st->isPermanent) 
    box_colour = SESS_COL_PERMANENT;
  else if (st->isDestroyed)
    box_colour = SESS_COL_DESTROYED;
  else
    box_colour = SESS_COL_KEPT_ALIVE;
  
  if (st->key == _This_session)
    box_colour = SESS_COL_CURRENT;

  return box_colour;
} /* stBoxColour */
/************************************************************/  

static void sessionPick(int box)
{ 
  ST *st = NULL ;

  if (!box2sess)
    messcrash ("No boxes drawn yet for sessionPick");

  /* redraw the old box in its original colours */
  if (selectedSessionBox)
    {
      st = arrp(box2sess, selectedSessionBox, ST);
      graphBoxDraw(selectedSessionBox, BLACK, stBoxColour(st));
    }

  if (box == 0)
    {
      /* click on background unselects the currently selected box */
      sessionChosen = 0;
      selectedSessionBox = 0;
      return;
    }

  /* box is not the background */

  st = arrayp(box2sess, box, ST); /* use of arrayp protects against
				     boxnumbers going beyond the 
				     bos2sess index */

  if (st->number)		/* clicked on a session object box 
				   (not just a blank ST struct) */
    {
      if (!st->isDestroyed && st->key != _This_session)
	sessionChosen = st->key ;
      else
	sessionChosen = 0 ;
      
      /* draw the selection red */
      graphBoxDraw(box, BLACK, SESS_COL_SELECTED) ;
      
      /* double click displays the session object */
      if (box == selectedSessionBox && st->key != 1)
	display(st->key, 0, 0) ;
      
      selectedSessionBox = box ;
    }
  else
    {
      /* we clicked on a box that wasn't a session object's box */
      sessionChosen = 0 ;
      selectedSessionBox = 0;
    }

  return;
} /* sessionPick */

/************************************************************/  
/*
 this would allow complete versioning, by removing te
 isOtherSesssion blocking of write access, but it
is both too complex to be useful, and bugged in
a complex way when destroying former sessions on various barnches

it may be usefull however to recover from a fatal crash
*/
static void sessDispDestroyLastSession (void)
{ KEY father, grandFather ;
  ST *st, *st1 ;
  int i ;
  extern VoidRoutine messcrashroutine ;   
  OBJ Session ;
  int  a1 = 0, a3 = 0, b3 = 0, c3 = 0, as = 0 ;
  BLOCK theSuperBlock ;
  Array sessionTree;
  char timeBuf[25] ;

  if (isWriteAccess())
    { messout("You still have write access\n"
	      "First save your previous work") ;
      return ;
    }
 
  if (!I_own_the_process () || 
      !sessionFilName ("wspec/passwd", "wrm", "a"))
    { messout
	("Sorry, to destroy a session you must own the "
	 "acedb executable and the file wspec/passwd.wrm") ;
      return ;
    }

  father = thisSession.from ; /* only allowed murder */
  if (!father)
    { messout ("First select a living session.") ;
      return ;
    }

  /***** find the predecessor of the last live session 
	 (grandFather of thisSession) *******************/

  sessionTree = sessionTreeCreate(FALSE);

  for (i = 0, st = arrp(sessionTree,0,ST);
       i < arrayMax(sessionTree) ; st++, i++)
    if (st->key == father)
      break ;

  if (i >= arrayMax(sessionTree))
    { 
      messout("Failed sorry") ;
      sessionTreeDestroy (sessionTree);
      return ;
    }

  grandFather = st->ancester ; /* mieg 2001. was father ; */
  for (i = 0, st1 = arrp(sessionTree,0,ST);
       i < arrayMax(sessionTree) ; st1++, i++)
    if (st1->key == st->ancester)  /* mieg 2001. was father ; */
      break ;

  if (i >= arrayMax(sessionTree))
    { 
      messout("Failed sorry") ;
      sessionTreeDestroy (sessionTree);
      return ;
    }

  if (!grandFather || !st1 || iskey(st1->key) != 2 || st1->isDestroyed)
    { messout ("You cannot destroy the last session, "
	       "because the previous one is dead") ;
      sessionTreeDestroy (sessionTree);
      return ;
    }

  sessionTreeDestroy (sessionTree);

  /***** found grandFather ******/

  if (!messQuery(messprintf("%s%s%s%s",
 "Do you really want to destroy the last session,\n",
 "this is an irreversible operation which is designed to help ",
 "recover from some rare fatal disk errors\n",
 "It is probably not useful in other cases."))) 
    return ;

  if (!messQuery(messprintf("%s%s%s",
 "Are you really sure ?\n\n",
 "If you proceed, acedb will destroy the last session, exit, ",
 "and you will have to restart the code")))
      return ;

  if (!checkWriteAccess())
    return ;

  /* On we go, no possible recovery ! */
  messcrashroutine = 0 ;  
  sessionForbidWriteAccess(); /* so we can never re-gain write access,
				 clears the readlock and
				 redraws if nec to change buttons etc. */

  Session = bsCreate(grandFather) ;
  if(!Session)
    messcrash("sessionDestroy() - Sorry, I cannot find grandfather '%s'",
	      name(grandFather));

  if(!bsGetData(Session,_GlobalLex,_Int,&a1) ||
     !bsGetData(Session,_Session,_Int,&as)||
     !bsGetData(Session,_VocLex,_Int,&a3) ||
     !bsGetData(Session,_bsRight,_Int,&b3)  ||
     !bsGetData(Session,_bsRight,_Int,&c3)  ||
         /* Hasher only introduced in release 1.4  !bsGetData(Session,_bsRight,_Int,&h0)  || */
     !bsGetData(Session,_CodeRelease,
	    _Int,&thisSession.mainCodeRelease) ||
     !bsGetData(Session,_bsRight,
		_Int,&thisSession.subCodeRelease) ||
     !bsGetData(Session,_DataRelease,
	    _Int,&thisSession.mainDataRelease) ||
     !bsGetData(Session,_bsRight,
	    _Int,&thisSession.subDataRelease) ) 
    messcrash("Anomaly while recovering info from previous session, sorry") ;

  bsDestroy(Session) ;

  /* proceed to change the superblock */
  /* decrement session number */
      
  diskblockread (&theSuperBlock,1) ;

  theSuperBlock.gAddress = thisSession.gAddress = a1 ;
  theSuperBlock.mainRelease = aceGetVersion() ;
  theSuperBlock.subCodeRelease = aceGetRelease() ;

#ifdef ACEDB5   
  theSuperBlock.session =  as ;
#else
  theSuperBlock.h.session =  as ;
#endif
  diskWriteSuperBlock(&theSuperBlock) ;	/* writes and zeros the BATS */
  messdump ("DESTROYED SESSION %d %s\n", 
	    thisSession.session-1, timeShow(timeNow(), timeBuf, 25)) ;
  messout("Last session saved, I will now exit, please restart acedb") ;

  filtmpcleanup () ;		/* deletes any temporary files made */
				/* NOTE, readlocks already cleared */
  exit (0) ;
} /* sessDispDestroyLastSession */

/****************************************************************/

static void sessDispFixSession(void)
{ 
  OBJ Session ;
  char *cp = 0 ;
  ST *st;

  st = arrp (box2sess, selectedSessionBox, ST);

  if (!st->number)
    {
      messout ("No session selected.");
      return;
    }

  if (st->isDestroyed)
    {
      messout ("Cannot fix a dead session.") ;
      return ;
    }

  if (st->isPermanent)
    { messout ("Selected session is already permanent.") ;
      return ;
    }

  if (iskey(sessionChosen) != 2)	/* check that it is still in a 
				   cache or on disk */
    { 
      messout ("Cache inconsistency:\n"
	       "This session no longer exists!") ;
      sessDispDraw() ;		/* re-draw */
      return ;
    }

  if (!checkWriteAccess())	/* will output a message */
    return ;

  /* I save once here to let the system a chance to crash if need be */
  saveAll() ; 
  Session = bsUpdate(sessionChosen) ;

  if (bsFindTag(Session, _Destroyed_by_session))
    { 
      messerror("sessiontree inconsistency: "
		"selected session should not be dead!") ;
      bsDestroy(Session) ;
      sessDispDraw() ;
      return ;
    }

  if (st->number != 1 && !bsFindTag(Session, _Up_linked_to))
    /* for the first session it is OK not to have an uplink */
    { 
      messout("I cannot rename a session created before code 1-6") ;
      bsDestroy(Session) ;
      return ;
    }

  cp = 0 ;
  bsGetData(Session, _Session_Title,_Text,&cp) ;
  
  if (messPrompt("Rename this session and make it permanent",cp ? cp : "","t") &&
      ( cp = freeword()))
    { if(strlen(cp) > 80)
      *(cp + 80) = 0 ;
    bsAddTag(Session, _Permanent_session) ;
    bsAddData(Session, _Session_Title, _Text, cp) ;
    }
  bsSave(Session) ;		/* also destroys the OBJ */
  /* This save will rewrite Session to disk without altering the BAT */
  saveAll() ;


  sessDispDraw() ;
} /* sessDispFixSession */

/****************************************************************/

static void sessDispUnFixSession(void)
{ 
  OBJ Session ;
  ST *st;

  st = arrp (box2sess, selectedSessionBox, ST);

  if (!st->number)
    {
      messout ("No session selected.");
      return;
    }

  if (st->isDestroyed)
    { 
      messout ("Cannot unfix a dead session.");
      return;
    }

  if (!st->isPermanent)
    { 
      messout ("Can only unfix a permanent session.");
      return;
    }

  if (st->isReadlocked)
    { 
      messout ("Somebody is still working on this session\nCannot unfix this session now.");
      return;
    }

  if (iskey(sessionChosen) != 2)	/* check that it is still in a 
				   cache or on disk */
    { 
      messout ("Cache inconsistency:\n"
	       "This session no longer exists!") ;
      sessDispDraw() ;		/* re-draw */
      return ;
    }

  if (!checkWriteAccess())	/* will output a message */
    return ;

  /* I save once here to let the system a chance to crash if need be */
  saveAll() ; 
  Session = bsUpdate(sessionChosen) ;

  /* this check should not be necessary, but we better make sure,
     in case the database was updated under our feet or so */
  
  if (!bsFindTag(Session, _Destroyed_by_session))
    {
      if (bsFindTag(Session,_Permanent_session))
	{
	  /* remove the tag that was accessed last - the permanent-tag */
	  bsRemove(Session);
	  bsSave(Session);	/* also destroys the OBJ */
	  /* This save will rewrite Session to disk 
	     without altering the BAT */
	  saveAll() ; 
	}
      else
	{			/* oops, it is no longer permanent */
	  messerror("sessiontree inconsistency: "
		    "selected session should be permanent!") ;
	  bsDestroy(Session) ;
	}
    }
  else
    {				/* oops it has since been destroyed */
      messerror("sessiontree inconsistency: "
		"selected session should not be dead!") ;
      bsDestroy(Session) ;
    }

  sessDispDraw() ;		/* re-draw, if we had detected an
				   inconsistency, these out-of-date
				   objects will now have been re-read */

  return;
} /* sessDispUnFixSession */

/****************************************************************/

static void sessDispToggleDeadShown (void)
{
  IsDeadSessionShown = IsDeadSessionShown ? FALSE : TRUE ;
  sessDispDraw() ;
} /* sessDispToggleDeadShown */

/****************************************************************/

static void sessDispDrawBranch(Array sessionTree, ST* st, 
			       Array xStart, Array xEnd, int y0)
{ 
  int i, gen = st->generation;
  int treeMax = arrayMax(sessionTree) ;
  ST *st1, *stnew ;
  int box;

  /* Draw self */
  if (gen)
    {
      box = graphBoxStart() ;
      stnew = arrayp(box2sess, box, ST);
      memcpy (stnew, st, sizeof(ST));

      st->len = strlen(name(st->key)) ;

      if(!array(xStart, gen, int))
	/* nothing drawn on level `g' yet */
	{
	  if (gen < 2)
	    /* start at the left edge */
	    array(xStart, gen, int) = 4 ;
	  else
	    /* start below the rightmost box of the previous generation */
	    array(xStart, gen, int) = array(xStart, gen-1, int) ;
	}
      else
	/* something drawn at that level yet */
	{
	  /* start this box just right of the end of the
	     rightmost box at that level */
	  array(xStart, gen, int) = array(xEnd, gen, int) + 2 ;
	}

      array(xEnd, gen, int) = array(xStart, gen, int) + st->len;

      st->x = array(xStart, gen, int) ;
      st->y = y0 + 3*gen ;
      graphText(name(st->key), st->x, st->y ) ;
      graphBoxEnd() ;

      /* maintain index of box-number to session object mapping */
      
      graphBoxDraw(box, BLACK, stBoxColour(st));

      /* draw the connecting lines */
      for (i = 0, st1 = arrp(sessionTree,0,ST); i < treeMax; st1++, i++)
	if (st->ancester == st1->key)
	  graphLine(st1->x + st1->len/2, st1->y + 1.0, 
		    st->x + st->len/2 ,  st->y - .0) ;
    }

  /* Recursion */
  for (i = 0, st1 = arrp(sessionTree, treeMax - 1,ST); 
       i < treeMax; st1--, i++)
     if (st1->sta == st)
       sessDispDrawBranch(sessionTree, st1, xStart, xEnd, y0) ;

  return;
} /* sessDispDrawBranch */

/****************************************************************/

static void sessDispDraw() 
{
  int i, max = 0 , maxGen = 0 , n, y0 ;
  int box;
  ST *st ;
  Array levelSessBoxStarts = arrayCreate(12, int) ;
  Array levelSessBoxEnds = arrayCreate(12, int) ;
  Array sessionTree = sessionTreeCreate(IsDeadSessionShown);

  if (!graphActivate(sessionSelectionGraph))
    messcrash ("sessDispDraw() - Can't activate graph.");
      
  graphClear() ;

  graphButtons(sessionMenu, 3.0, 2.0, 40) ;
  
  y0 = 9 ;
  graphText("Double click on any previous session for text info.", 3, y0++) ;
  graphText("Restart from any live session (in read only mode)", 3, y0++) ;
  graphText("Fix any live session, to prevent its future destruction", 3, y0++) ;
  y0++ ;
  graphText("Colours: ", 3, y0) ;

  box = graphBoxStart ();
  graphText("live (permanent)", 12, y0);
  graphBoxEnd();
  graphBoxDraw(box, BLACK, SESS_COL_PERMANENT);

  box = graphBoxStart ();
  graphText("live (still in use)", 31, y0);
  graphBoxEnd();
  graphBoxDraw(box, BLACK, SESS_COL_READLOCKED);

  y0 += 2;

  box = graphBoxStart ();
  graphText("live (kept alive)", 12, y0);
  graphBoxEnd();
  graphBoxDraw(box, BLACK, SESS_COL_KEPT_ALIVE);

  box = graphBoxStart ();
  graphText("dead", 31, y0);
  graphBoxEnd();
  graphBoxDraw(box, BLACK, SESS_COL_DESTROYED);

  box = graphBoxStart ();
  graphText("selected", 37, y0++);
  graphBoxEnd();
  graphBoxDraw(box, BLACK, SESS_COL_SELECTED);

  graphPop() ;

  /* maps box numbers to session tree objects */
  box2sess = arrayReCreate(box2sess, 50, ST);

  n = arrayMax (sessionTree) ;
  for (i = 0, st = arrp(sessionTree, n - 1, ST); i < n; st--, i++)
    { 
      if (st->generation == 1)
	sessDispDrawBranch(sessionTree, st,
			  levelSessBoxStarts, levelSessBoxEnds, y0) ;
      if (st->generation > maxGen)
	maxGen = st->generation ;
    }

  max = 66;			/* width of writing */
  for (i = 0; i < arrayMax(levelSessBoxEnds); i++)
    if (max < (arr(levelSessBoxEnds, i, int) + 3))
      max = (arr(levelSessBoxEnds, i, int) + 3);
 
  graphTextBounds (max, y0 + 3 * maxGen + 2) ;
  graphRedraw() ;

  arrayDestroy (levelSessBoxStarts);
  arrayDestroy (levelSessBoxEnds);
  sessionTreeDestroy (sessionTree);

  return;
} /* sessDispDraw */

/************************* eof ******************************/
