/*  File: action.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **    Calls an action application.
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  3 14:44 1998 (edgrif)
 * * Dec  3 14:44 1998 (edgrif): Change calls to new interface to aceversion.
 * * Dec  2 15:56 1998 (edgrif): Fix minor bug in acedb version message.
 * * (I assume these changes made by rbrusk)
 *		-	'externalAsynchroneCommand' must return a value
 *		-	'SIGUSR1' : undefined in WIN32; user can user windows mgr "End Task"
 * * May  7 17:40 1996 (rd)
 * * May 12 17:39 1993 (mieg): pick_me_to_call now expects: command parameters
 * Created: Fri Jan 10 23:40:16 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: action.c,v 1.4 2015/08/18 23:24:07 mieg Exp $ */

#include "acedb.h"
#include "aceversion.h"
#include "bs.h"
#include "lex.h"
#include "display.h"   /* self declaration */
#include "systags.h"
#include "sysclass.h"
#include "tags.h"
#include "pick.h"
#include "mytime.h"
#include "call.h"
#include "session.h"

void externalDisplay (KEY key)  
{
  OBJ obj = bsCreate(key) ;
  char *dir = 0, *prog = 0, *args = 0 ;

  if (!obj)
    messcrash ("Pick_me_to_call called on unkown obj %s", name(key)) ;
  if (!bsGetData (obj, _Pick_me_to_call, _Text, &prog))
    { messout ("Sorry, a program and its parameters are expected here") ;
      goto abort ;
    }
 
  bsGetData (obj, _bsRight, _Text, &args) ;  /* get the parameters */
  callCdScript (dir, prog, args) ; 
  printf ("acedb call: %s %s %s\n", dir ? dir : "", prog , args ?  args: "") ;

 abort:
  bsDestroy (obj) ;
}

/**************************************************************/
/**************************************************************/

char* exportStackToMail (Stack s) /* called from droso code */
{ FILE *f = 0 ;
  char *cp, *nom = 0 ;
 
  f = filtmpopen (&nom, "w") ;
  
  if (!f)
    { messprintf ("Sorry, i cannot create a tmp file") ;
      return 0 ;
    }

  if (stackExists(s))
    { stackCursor (s, 0) ;
      while ((cp = stackNextText(s)))
	fprintf (f, "%s\n", cp) ;
    }
  filclose (f) ;
  fflush (f) ;
  return nom ;
}

/****************/

void acedbMailComments(void)
{ Stack s ;
  char * nm ;
  if (!messQuery ("Do you want to send a comment to the authors of this program"))
    return ;
  s = stackCreate (100) ;
  pushText(s, messprintf("Code release:   %s,  %s\n",
			 aceGetVersionString(), aceGetLinkDateString())) ;

  nm = exportStackToMail (s) ;
  stackDestroy (s) ;
  if (nm)
    callScript ("acedb_mailer", 
		messprintf("\"mieg@crbm.cnrs-mop.fr rd@sanger.ac.uk\" %s", nm)) ;
  
}

void acedbMailer (KEY key, KEYSET ks, Stack sText)
{ char *cp ;
  Stack  sm = stackCreate(50) ;
  int i, nmails = 0 ;
  OBJ obj ;
  int manip = 0 ;

  if (!keySetExists(ks))
    { ks = keySetCreate() ;  /* I will destroy at end */
      manip = 10 ;
    }
  if (key)
    { keySet(ks, keySetMax(ks)) = key ;
      manip++ ;
    }
  
  for (i = 0 ; i < keySetMax(ks) ; i++)
    { if (pickType(keySet(ks,i)) != 'B')
	continue ;
      obj = bsCreate(keySet(ks,i)) ;
  
      if (!obj)
	continue ;
      if (bsGetData(obj, _E_mail, _Text,&cp))
	do
	  { 
	    if (nmails++)
	      catText(sm, ", ") ;
	    else
	      pushText(sm,"\'") ;
	    catText(sm, cp) ;
	  } while (bsGetData(obj, _bsDown, _Text,&cp)) ;
      if (stackMark(sm) > 800)
	{ messout ("%d recipients is enough, I stop at %s",
		   nmails, name(keySet(ks,i))) ;
	  break ;
	}

      bsDestroy(obj) ;
    }
  catText (sm,"\'  ") ;
  catText (sm,exportStackToMail(sText) ) ;
  catText (sm," & ") ;
  
  if (!nmails)
    messout("No E_mail address in this keySet, sorry") ;
  else
    callScript ("acedb_mailer", stackText(sm,0)) ;

  stackDestroy(sm) ;

  if (manip >= 10)
    keySetDestroy(ks) ;
  else if (manip == 1)
    keySetMax(ks)-- ;
}

/**************************************************************/

/* call an external shell command and print output in a text_scroll window
   from Erik Sonnhammer, 92-01-12
*/

static void doExternalFileDisplay (const char *title, FILE *f, Stack s, BOOL isPipe)
{
  int line = 2, level = 0, max = 60 ;
  char *cp ;
  Graph old = graphActive() ; 
  Stack s1 = 0 ;

  graphCreate (TEXT_SCROLL, title ? title : "Action Display" , 0, 0, 0.7, 0.5);

  if (title)
    graphText(title, 1,1) ; 
  if (stackExists(s))
    {
      s1 = stackCreate (1000) ;
      stackCursor(s, 0) ;
      while ((cp = stackNextText (s)))
	{ catText(s1, cp) ; catText(s1, "\n") ; }
      level = freesettext (stackText(s1, 0),"") ;
    }
  else if (f)
    level = isPipe ? freesetpipe(f,"")  : freesetfile(f,"") ;

  if (level)
    {
      freespecial ("\n\t\"") ;
      while (line++, (cp = freecard(level))) /* closes the file */
	{ if (max < strlen(cp))
	  max = strlen(cp) ;
	graphText(cp, 0, line) ;
	}
      freeclose (level) ;
    }
  stackDestroy (s1) ;
  graphTextBounds (max, line);
  graphRedraw() ;
  graphActivate (old) ;
}

void externalFileDisplay (const char *title, FILE *f, Stack s)
{ doExternalFileDisplay (title, f, s, FALSE) ; }

void externalPipeDisplay (const char *title, FILE *f, Stack s)
{ doExternalFileDisplay (title, f, s, TRUE) ; }

void externalCommand (const char *command)
{ 
  char *cp, *cp0 = command && *command ? strnew(command, 0) : 0 ;
  
  if (cp0)
    {
      cp = cp0 ;
      while (*cp && *cp != ' ') cp++ ;
      *cp++ = 0 ; /* first space or terminal 0 */
    }
#if !defined(MACINTOSH)
  if (cp0 && *cp0)
    doExternalFileDisplay (command, callScriptPipe(command, cp), 0, TRUE) ;
#endif
  
  messfree (cp0) ;
}

#if defined(MACINTOSH)
BOOL externalAsynchroneCommand (const char *command, char *parms,
                                void *look, void(*g)(FILE *f, void *lk))
{ return ; }
#else
#include <signal.h>

/* This system is very simple minded, 
   We create a communication file xx->fName
   then call the command in the background
   and test for the existence of xx->fName.out 
   when a SIGUSER1 signal is received.
    If the outfile exists, the registered function is
   called on that output file, and cleanup is done
   */
   
    
static Array xxA = 0 ;

typedef struct xxLookStruct 
{ 
  int i ; /* index in xxA, 0: reusable */
  void (*g)(FILE *f, void *lk) ;
  char *command ;
  char *fName ;
  void *look ;
} XX ;

static void xxDestroy (XX *xx)
{ 	  /* I think that xx->look should be destroyed in the calling code */
  xx->i = 0 ;  /* slot can be reused */
  messfree(xx->command) ;
  filremove (xx->fName, "out") ;
  filtmpremove (xx->fName) ;
}

static void externalCompletion (int result)
{ FILE *f = 0 ;
  XX *xx ;
  int i ;

  if (!arrayExists (xxA))
    return ;

#if !defined(WIN32)
  signal (SIGUSR1, &externalCompletion) ;  /* could be called only once */
#endif

  for (i = 0 ; i < arrayMax(xxA) ; i++)
    {
      xx = arrp (xxA, i, XX) ;
      if (xx->i &&                               /* not destroyed */
	  filName (xx->fName, "out", "r") &&     /* no messout if absent  */
	  (f = filopen (xx->fName, "out", "r"))) /* would messout if absent */
	{
	  if (xx->g)  /* must destroy xx->look */
	    (*(xx->g))(f, xx->look) ;
	  xxDestroy (xx) ;
	}
    }
    
}

BOOL externalAsynchroneCommand (const char *command, const char *parms,
                                void *look, void(*g)(FILE *f, void *lk))
{ int i ;
  char *fName ;
  XX  *xx ; 
  FILE *f ;
  Stack s = 0 ;
  int pid = getpid() ;

  if (!xxA)
    xxA = arrayCreate (20, XX) ;
  f = filtmpopen (&fName, "w") ;
  if (!f)
    return FALSE ;
  fprintf(f, "// Call to: %s\n", command) ;
  filclose (f) ;

  for (i = 1 ; i < arrayMax(xxA) ; i++) /* start at i = 1 ! */
    if (!array (xxA, i, XX).i)
      break ;
  xx = arrayp (xxA, i, XX) ;  
  xx->i = i ;
  xx->command = strnew(command, 0) ;
  xx->look = look ;
  xx->g = g ;
  xx->fName = fName ;

#if !defined(WIN32)
  signal (SIGUSR1, &externalCompletion) ;  /* could be called only once */
#endif

  s  = stackCreate (250) ;
  pushText (s, messprintf(" %d ", pid)) ;
  catText (s, xx->fName) ;
  catText (s, "  ") ;
  if (*parms) catText (s, parms) ;
  catText (s, " &") ;
  
  if (
      (filName(messprintf("wscripts/%s", command),0,"x") || /* check name here since */
      filName(messprintf("%s", command),0,"x")) &&
      !callScript (command, stackText (s, 0))  /* because of &background, callScript always returns TRUE */
      )
    { 
      stackDestroy (s) ;
      return TRUE ;
    }
  else
    {
      messout ("Cannot find script %s", command) ;  
      if (xx->i && xx->g)
	(*(xx->g))(0, xx->look) ;
      xxDestroy (xx) ;
      return FALSE ;
    }
  return TRUE ;
}

#endif /* MACINTOSH */

/***************************************************************/
/***************************************************************/

typedef struct LOOKSTUFF
  { int   magic;        /* == MAGIC */
    KEY   key ;
    Graph graph ;
  } *LOOK ;

static MENUOPT actionMenu[] = {
  {graphDestroy, "Quit"},
  {help,"Help"},
  {graphPrint,"Print"},
  {displayPreserve,"Preserve"},
   {0, 0}
} ;

static int MAGIC = 723544 ;	/* use also for graphAssociate pointer */

#define LOOKGET(name)     LOOK look ; \
                          if (!graphAssFind (&MAGIC, &look)) \
		            messcrash ("graph not found in %s", name) ; \
			  if (!look) \
                            messcrash ("%s received a null pointer", name) ; \
                          if (look->magic != MAGIC) \
                            messcrash ("%s received a wrong pointer", name)

static BOOL actionExecute (LOOK look)
{
  graphText(messprintf ("Action on %s not yet code sorry", 
			name(look->key)), 3, 3) ;
  graphRedraw() ;
  return TRUE ;
}

static void actionDestroy (void)
{
  LOOKGET("actionDestroy") ;

  look->magic = 0 ;
  messfree (look) ;
  
  graphAssRemove (&MAGIC) ;
}

BOOL actionDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  LOOK look ;
  OBJ Action = 0 ;
  int _VAction = 0 ;
  KEY akey ;

  lexword2key("Action", &akey, _VMainClasses) ;
  _VAction = KEYKEY(akey) ;

  if (!key || !_VAction || class(key) != _VAction ||
      !(Action = bsCreate(key)))
    return FALSE ;

  bsDestroy(Action) ;

  look=(LOOK)messalloc(sizeof(struct LOOKSTUFF)) ;
  look->key = key ;
  look->magic = MAGIC;

  if (isOldGraph)
    { actionDestroy () ;
      graphClear () ;
      graphGoto (0,0) ;
      graphRetitle (name(key)) ;
      graphAssRemove (&MAGIC) ;
    }
  else 
    { if (!displayCreate (DtAction)) 
	{ messfree (look) ;
	  return FALSE ;
	}
      
      graphRetitle (name(key)) ;
      graphRegister (DESTROY, actionDestroy) ;
      graphMenu (actionMenu) ;
    }
  graphAssociate (&MAGIC, look) ;
  look->graph = graphActive () ;
  actionExecute (look) ;
  
  return TRUE ;
}

/**************************************************************/
