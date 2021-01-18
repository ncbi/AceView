/*  File: xclient.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@kaa.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 30 18:34 1998 (fw)
 * Created: Tue Jun 28 14:13:30 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: xclient.c,v 1.5 2020/05/30 16:50:36 mieg Exp $ */

#include "acedb.h"
#include "aceclient.h"
#include "parse.h"
#include "pick.h"
#include "sysclass.h"
#include "session.h"
#include "call.h"
#include "dump.h"
#include "lex.h"
#include "bs.h"
/* #include "classes.h" */
#include "client.h"
#include "table.h"
#include "main.h"

static char server[255] ; 
static char *journalFileName = 0 ; 
static int clientCount = 0 , serverCount = 0 ;

static int timeOut = 300 ;

static ace_handle *handle = 0 ;
static u_long port = DEFAULT_PORT;
static BOOL serverIsAllowed  = FALSE ; 
static void xClientCommit (void) ;
static KEYSET xClientGet (KEY key, char *quer, char *grp, BOOL getNeighbours) ;
static KEYSET xClientQuery (KEY key, char *quer, char *grp, BOOL getNeighbours) ;
static void  xClientSave (KEY key, int action) ;
static void  xClientAlias (KEY key, char *oldname, BOOL keepOldName) ;
static TABLE* xClientTableMaker (char *quer) ;
void xClientInit (void) ;

/***************************************************************************/
/***************************************************************************/

char *xClientServerName(void)
{ return server ; }

static BOOL getServerCount(void)
{
  /* asks server for its update-data number */
  BOOL ok = FALSE ;
  Stack sr = 0 ;
  int  encore = 0,erc, dummy, level, n ;
  char *cp ;
  unsigned char *reponse = 0 ;

  if (!serverIsAllowed)
    return FALSE ;
  
  if (!handle)
    xClientInit() ;
  if (!handle)
    return 0 ;

  erc = askServerBinary (handle, "Data_Version", &reponse, &dummy, &encore, 0) ;

  if (!reponse || erc)
    return FALSE ;
 
  sr = stackReCreate (sr, 1000) ;
  if (reponse)
    {
      catText (sr, (char*)reponse) ;
      free (reponse) ; /* NOT messfree */
    }
  stackCursor(sr, 0) ;
  level = freesettext (stackText(sr,0), 0) ;
  freespecial("\n\"\\") ;
  while (freecard(level))
    {
      while ((cp = freeword()))
	{
	  if (!strcmp(cp,"version"))
	    { 
	      if (!freeint(&n))
		goto fin ;
	      serverCount = n ;
	      ok = TRUE ;
	    }
	}
    }
fin:
  freeclose (level) ;
  stackDestroy (sr) ;
  
  clientCount = thisSession.subDataRelease ;
  return ok ;
}

static BOOL xClientState (int n)
{
  switch (n)
    {
    case 0: serverIsAllowed = FALSE ; externalServer = 0 ; break ;
    case 1: 
      if (!serverIsAllowed) xClientInit () ;
      break ;
    default: break ;
    }
  return serverIsAllowed ;
}

void xClientClose (void)
{ 
  closeServer (handle) ;
  handle = 0 ; 
}

void xClientLost (void)
{ 
  xClientClose () ;
  mainPickReport (2) ; /* report in main window */
}


static BOOL xClientRefresh (void)
{
  register int i ;
  getServerCount() ;
  
  if (!serverIsAllowed) return FALSE ; /* cannot get serverCount */
  if (!serverCount ||   /* continuous refresh */
      serverCount > clientCount)
    for (i = 0 ; i < 256 ; i++)
      lexClearClassStatus(i, SERVERSTATUS) ;   
  thisSession.subDataRelease = serverCount ;
  return TRUE ;
}


void xClientInit(void)
{ static int firstPass = 0 ;

  char *cp ;
  int i, level ; 
  FILE *fil=NULL;

  if (firstPass==0)
    { firstPass = 1 ;

      port = 0 ; *server = 0 ;
      cp = sessionFilName ("wspec/server", "wrm", "r") ;
      if (!cp ||
	  !(fil = filopen (cp, 0, "r")))
	messExit ("initXClient cannot open server configuration file: %s", cp) ;
      level = freesetfile(fil, 0) ;
      while (freecard(level))
	{ cp = freeword() ;
	  if (cp)
	    { 
	      if (!strcmp (cp,"EXTERNAL_SERVER") &&
		  (cp = freeword()) )
		strncpy (server, cp, 254) ;
	      else if (!strcmp (cp,"EXTERNAL_PORT") &&
		       freeint (&i))
		port = i ;
	    }
	}
      freeclose (level) ;
    }

  if (*server && port && !handle)
    handle = openServer (server, port, timeOut) ;
  if (!handle)
    { if (firstPass==1) messout (
			 "Connection to server %s, port %lu  failed, sorry, please check wspec/server.wrm", 
			 server, port) ;
    }
  else
    {
      externalServer =  xClientQuery ;  /* success */
      externalServerState = xClientState ;
      serverIsAllowed = TRUE ;
      externalSaver = xClientSave ;
      externalCommit = xClientCommit ;
      externalTableMaker = xClientTableMaker ;
      externalServerName = xClientServerName ;
      externalAlias = xClientAlias ;
      if (firstPass == 1)
	{
	  xClientRefresh() ;
	  fprintf(stderr, "Connection to external Server %s, port %lu established\n",
		  server, port) ;
	}
    }
  firstPass = 2 ;
}

/***************************************************************************/
/* If key == -1, export the keyset oldKeySetForServer to the server
 * this is useful before a show (to obtain all data at once
 * or before a table, to run on the active table
 */

static KEYSET xClientGet (KEY key, char *quer, char *grp, BOOL getNeighbours)
{ int   i, erc ;
  char *cp ;
  int  encore = 0, dummy, 
       count = -1 ;		/* init to avoid warning (?) */
  KEYSET kk = keySetCreate (), oldSet = 0 ;
  KEY  *kp ;
  Stack sq = 0, sr = 0 ;
  BOOL pasFini = TRUE ;
  static int nerr = 0 ;
  BOOL show = FALSE, pleaseSave = FALSE ;
  unsigned char *reponse = 0 ;
  extern BOOL   X_CLIENT_PARSING ;  /* in lexsubs4.c */
  extern KEYSET  oldSetForServer ;
  static KEYSET requiredSet = 0 ;

  requiredSet = keySetReCreate (requiredSet) ;
  if(!serverIsAllowed)
    return kk ;

  /* do not ask for whole database with a name like *?* */
  if (key == -1)
    { key = 0 ; oldSet =  oldSetForServer ;
      if (oldSet && !keySetMax(oldSet)) return kk ;
      show = getNeighbours ; getNeighbours = FALSE ;
    }
  if (key == -4) { key = 0 ; pleaseSave = TRUE ; }
  if (quer || grp)
    { count = key ; key = 0 ; }
  if (key)
    {
      if ((lexGetStatus(key) & SERVERSTATUS) || !dumpClassPossible(class(key), 'a'))
	return kk ;
      show = TRUE ;
      cp = name(key) ;
      if (!cp || !*cp) return kk ;
      cp-- ; while (*++cp) if (*cp != '*' && *cp != '?') break ;
      if (!*cp) return kk ;
      keySet (requiredSet,0) = key ;
    }
  else
    if (!pleaseSave && !oldSet && ((!quer && !grp) || (quer && !*quer) || (grp && !*grp)))
      return kk ;
  
  if (nerr > 20)
    { nerr = 0 ; xClientLost() ; }

  messStatus ("Reading from server") ;
  sq = stackReCreate (sq, 1000) ;
  sr = stackReCreate (sr, 1000) ;

  /*
    if (oldSet && !quer && !grp && !show) ;
    then i am just passing my active set from tableMaker
    */

  if (oldSet)
    { int tt1 = 0, tt2, j ;
      pushText(sq, "KeySet-Read = ") ;
      i = keySetMax(oldSet) ; kp = arrp(oldSet, 0, KEY) - 1 ;
      j = 0 ;
      while (kp++, i--)
	{ tt2 = class(*kp) ;
	  if (pickList[tt2].protected)
	    continue ;
	  if (!KEYKEY(*kp))
	    continue ;
	  if (!quer && !grp && show)  /* anticipate on these objects */
	    {
	      if (lexGetStatus(*kp) & SERVERSTATUS) 
		continue ;
	      else
		keySet (requiredSet,j++) = *kp ;
	    }
	  if (tt2 != tt1)
	    catText(sq,messprintf("%s ",className(*kp))) ;
	  tt1 = tt2 ;
	  catText(sq,messprintf("%s;",freeprotect(name(*kp)))) ;
	}
      if (!quer && !grp && show && !j) return kk ;
    }
  else if (quer)
    { pushText(sq,"Query ") ;
    catText (sq, quer) ;
    }
  else if (pleaseSave)
    { pushText(sq,"Save ") ;
    pasFini = FALSE ;
     }
  else if (grp)
    { pushText(sq,"Grep ") ;
    catText (sq, grp) ;
    }
  else if (getNeighbours)
    pushText(sq,messprintf("Query Find %s \"%s\" ; NEIGHBOURS \n", 
			  className(key), name(key))) ;
  else
    pushText(sq,messprintf("Query Find %s \"%s\" \n", 
			  className(key), name(key))) ;
  
lao:
  if (!handle)
    xClientInit() ;
  if (!handle)
    return kk ;
suiteReponse:
  erc = askServerBinary (handle, stackText(sq,0), &reponse, &dummy, &encore, 0) ;

  if (erc || !reponse)
    { handle = 0 ;
      fprintf(stderr, "call_ace_server, xCientGet,  return with : %i\n",erc);
      nerr++ ;
      return kk ;
    }

 
  if (reponse)
    {
      catText (sr, (char*)reponse) ;
      free (reponse) ; /* NOT messfree */
    }
  if (encore) 
    { stackClear(sq) ; pushText (sq,"encore") ; goto suiteReponse ;}
  
  if (pasFini)
    { 
      pasFini = FALSE ;
      stackClear(sq) ;
      if (oldSet)
	{ oldSet = 0 ; /* was not allocated here, is just a pointer really */
	  if (quer)
	    {  pasFini = TRUE ;
	    pushText(sq,"Query ") ;
	    catText (sq, quer) ;
	    }
	  else if (grp)
	    {  pasFini = TRUE ;
	    pushText(sq,"Grep ") ;
	    catText (sq, grp) ;
	    }
	  else if (show)
	    pushText(sq,"Show -a") ; 
	  else
	    return kk ;  /* just made oldks the active list on the server */
	}
      else if (show)
	pushText(sq,"Show -a") ;
      else
	{ pushText(sq,"List -a") ;
	  if (count) catText (sq, messprintf(" -c %d ", count)) ;
	}
      goto lao ;
    }
 
  stackCursor(sr, 0) ;

  X_CLIENT_PARSING = TRUE ;
  if (show)  /* server did a show */
    while((cp = stackNextText(sr)))
      { kk = keySetReCreate (kk) ;
        sessionGainWriteAccess() ;
        parseBuffer (cp, kk) ;
	i = keySetMax(kk) ; kp = arrp(kk, 0, KEY) - 1 ;
	while (kp++, i--)
	  if (*kp) lexSetStatus (*kp, SERVERSTATUS) ; /* prevent recursions */
      }
   else
     { 
       int cla, j = 0, level = freesettext (stackText(sr,0), "") ;
       kk = keySetReCreate (kk) ; i = 0 ;
       while (freecard(level))
	 {
	   cp = freeword() ;
	   if (!cp || !(cla = pickWord2Class(cp)))
	     continue ;
	   freenext () ; freestep(':') ; freenext() ;
	   cp = freeword() ;
	   if (!cp) continue ;
	   lexaddkey (cp, &key, cla) ;
	   if (i++) keySet(kk, j++) = key ; /* jump keyset answer */
	 }
     }
  if (key && getNeighbours && show &&  !(lexGetStatus(key) & TOUCHSTATUS))
    {
      KEYSET ks1 = bsKeySet(key) ; /* server should have exported all these */
      if (ks1)
	{
	  keySetDestroy (requiredSet) ;
	  requiredSet = ks1 ;
	}
    }
  if (show)
    {
      i = keySetMax(requiredSet) ; kp = arrp(requiredSet, 0, KEY) - 1 ;
      while (kp++, i--)
	if (*kp) lexSetStatus (*kp, SERVERSTATUS) ; /* since i asked explicitelly */
    }
  X_CLIENT_PARSING = FALSE ;
  keySetSort (kk) ;
  keySetCompress (kk) ;
  return kk ;
}

/***************************************************************************/

static BOOL xClientWrite (char *buffer)
{
  int  erc, encore = 3,   /* ACE_IN */
       dummy;
  static int nerr = 0 ;
  unsigned char *reponse = 0 ;
  static int nn = 0 ;

  if(!serverIsAllowed)
    return FALSE ; 

  messStatus ("Writing to server") ;
  if (nerr == 20)
     { messout("More than 20 errors when writing to server, i close the connection") ;
       nerr++ ;
     }
  if (nerr > 20)
     { xClientClose() ; return FALSE ; }

  if (!handle && nn < 10)
    { nn++ ;
      handle = openServer (server, port, timeOut) ;
    }
  if (!handle)
    return FALSE ; 
  if (nn > 0) nn-- ;
  erc = askServerBinary (handle, buffer, &reponse, &dummy, &encore, 0) ;

  if (erc != -3 || !reponse)
    { handle = 0 ;
      fprintf(stderr, "call_ace_server return with : %i\n",erc);
      nerr++ ;
      return FALSE ;
    }

  /*
  s = stackReCreate(s,1000) ;
  pushText (s, reponse) ;
  stackCursor(s, 0) ;
  */
  if (reponse) free (reponse) ; /* NOT messfree */
  
  return TRUE ;
}

/***************************************************************************/
/* if quer or grp, key is the number of lines needed 
   key == -1  -> query from keySetForServer
       == -2  -> disable server
       == -3  -> reenable
       == -4  -> save //  see above
       */
 
static KEYSET xClientQuery (KEY key, char *quer, char *grp, BOOL getNeighbours)
{ 
  KEYSET kk = 0 ;
  extern BOOL   X_CLIENT_PARSING ;  /* in lexsubs4.c */
  int i ;

  if (key == -2) { xClientState(0) ; return 0 ; }
  if (key == -3) { xClientState(1) ; return 0 ; }
  if (key == -4) { xClientGet (key, 0, 0, 0) ; return 0 ; }
  if (quer|| grp)invokeDebugger() ;
  if (!serverIsAllowed)
    return keySetCreate () ;
  if (quer)
    return xClientGet (key, quer, 0, getNeighbours) ;
  else if (grp)
    return xClientGet (key, 0, grp, getNeighbours) ; 
  else if (key == -1)
    return xClientGet (key, 0, 0, getNeighbours) ; 
  else if (key)
    {
      if(KEYKEY (key) && 
	 !X_CLIENT_PARSING &&
	 class(key) > 4 &&
	 dumpClassPossible(class(key), 'a') &&
	 !(lexGetStatus(key) & SERVERSTATUS) &&
	 !(lexGetStatus(key) & TOUCHSTATUS) &&
	 class(key) != _VMainClasses)
	{ 
	  invokeDebugger() ;
	  kk = xClientGet (key, 0, 0, getNeighbours) ;
	  keySetDestroy (kk) ;
	}
      else
	{
	  i = KEYKEY (key) ;
	  i = !X_CLIENT_PARSING ;
          i = class(key) ;
	  i = dumpClassPossible(class(key), 'a') ;
	  i = !(lexGetStatus(key) & SERVERSTATUS) ;
          i = !(lexGetStatus(key) & TOUCHSTATUS) ;
	  i =  class(key) ;
	  i = _VMainClasses ;
	}
    }
  return 0 ;
}

static TABLE* xClientTableMaker (char *quer)
{ 
  Stack sr = 0, sq = 0 ;
  TABLE *tt = 0 ;
  int  nerr = 0, encore = 0,erc, dummy, level ;

  unsigned char *reponse = 0 ;

  if (!serverIsAllowed || !quer || !*quer)
    return 0 ;
  
  if (nerr > 20)
    { nerr = 0 ; xClientLost() ; }

  if (!handle)
    xClientInit() ;
  if (!handle)
    return 0 ;

  mainActivity ("Reading from server") ;
  sq = stackReCreate (sq, 1000) ;
  sr = stackReCreate (sr, 1000) ;
  pushText (sq, quer) ;

suiteReponse:
  erc = askServerBinary (handle, stackText(sq,0), &reponse, &dummy, &encore, 0) ;

  if (!reponse || (erc && (erc + encore)))
    { handle = 0 ;
      fprintf(stderr, "call_ace_server, xCientGet,  return with : %i\n",erc);
      nerr++ ;
      stackDestroy (sq) ; stackDestroy (sr) ;
      return 0 ;
    }

 
  if (reponse)
    {
      catText (sr, (char*)reponse) ;
      free (reponse) ; /* NOT messfree */
    }
  if (encore) 
    { stackClear(sq) ; pushText (sq,"encore") ; goto suiteReponse ;}
  
  stackCursor(sr, 0) ;
  fprintf(stderr, stackText(sr,0)) ;
  level = freesettext (stackText(sr,0), 0) ;
  freespecial("\n\"/\\") ;
  tt = tableDoParse (level, 0) ;
  freeclose (level) ;
  stackDestroy (sr) ;

  return tt ;
}

/***************************************************************************/
/***************************************************************************/

static Stack getAceFile (level)
{ static Stack s = 0 ;
  char *cp, cutter ;
  int nn = 0 ;

  s = stackCreate (30000) ;

  while (freecard (level))
    { cp = freewordcut("\n" ,&cutter) ;
      if (cp && *cp)
	{ 
	  nn++ ;
	  catText (s, cp) ;
	}
      catText (s, "\n") ;
    }
  if (nn)
    return s ;
  stackDestroy(s) ;
  return 0 ;
}

static Stack getJournal(void)
{
  int level ;
  FILE *f = 0 ;

  if (!journalFileName)
    return 0 ;
  f = filopen (journalFileName,"","r") ;
  if (!f) return 0 ;
  level = freesetfile(f, "") ;
  freespecial ("\n\"/\\\t@") ;
  return getAceFile (level) ;
}

static void  xClientCommit (void)
{
  static BOOL failed = FALSE ;
  Stack s = 0 ;
  char *cp ;
  BOOL ok = TRUE ;

  if (failed) return ; /* all commit shoud succeed */

 
    
  s = getJournal () ;
  if (s)
    {
      if (messQuery("Do you want to send your editions back to the server"))
	{
	  stackCursor (s, 0) ;
	  ok = TRUE ;
	  while ((cp = stackNextText(s)))
	    ok &= xClientWrite (cp) ;
	}
      else
	ok = FALSE ;
    }

  stackDestroy (s) ;
  if (ok) /* success, throw away journal */
    {
      /*
	SAVE
	uncomment this line if you wish to force a save on the server
	this command issues a Save request
	but this may slow the code quite a bit

	xClientQuery ((KEY)(-4), 0, 0, 0) ;
	*/
      if (journalFileName)
	{
	  system(messprintf("touch %s.done ; cat %s >> %s.done",
			    journalFileName,journalFileName,journalFileName)) ;
	  filremove (journalFileName,"") ;
	}
      journalFileName = 0 ;
      return ;
    }
  failed = TRUE ;
}

/***************************************************************************/
/***************************************************************************/
static char* oldaliasname = 0 ;
static void  xClientSave (KEY key, int action) 
{
  static char* tmpName1 = 0 ;
  char *cp, *tmpName2 = 0 ;
  FILE *f ; 
  static BOOL isKilling = FALSE ;
  extern BOOL X_CLIENT_PARSING, isXrefing ;

  if (X_CLIENT_PARSING || !isXrefing  || (isKilling && action != 5) ||
      !dumpClassPossible (class(key),'a') )
    return ;

  if (!journalFileName) 
    { cp = filName("clientUpdate","","a") ;
    if (cp) journalFileName = strnew(cp, 0) ;
    }

  if (!journalFileName) 
    return ;

  switch (action)
    {
    case 1:  /* save a B object, phase 1: old version */
      if (tmpName1) { filtmpremove(tmpName1) ; messfree (tmpName1) ;tmpName1 = 0 ;}
      f = filtmpopen (&tmpName1,"w") ;   
      if (!f && tmpName1)  { filtmpremove(tmpName1) ; tmpName1 = 0 ; return ; }
      cp = tmpName1 ; tmpName1 = strnew(cp, 0) ;
      dumpKey(key, f, 0) ;
      filclose(f) ;
      break ;
    case 2:  /* save a B object, phase 2: new version */
      if (!tmpName1) return ;
      f = filtmpopen (&tmpName2,"w") ;   
      if (!f && tmpName2)  { filtmpremove(tmpName2) ; tmpName2 = 0 ; return ; }
      dumpKey(key, f, 0) ; 
      filclose(f) ;

      callScript("/bin/env acediff", messprintf(" %s %s >> %s", tmpName1, tmpName2, journalFileName)) ;
      filtmpremove(tmpName1) ; messfree(tmpName1) ; tmpName1 = 0 ;
      filtmpremove(tmpName2) ; tmpName2 = 0 ;
      break ;
    case 3:  /* dump */
      f = filopen (journalFileName,"","a") ;
      if (f)
	{ fprintf(f,"\n") ;
	  dumpKey(key, f, 0) ; 
	  fprintf(f,"\n") ; 
	}
      filclose(f) ;
      break ;
    case 4:  /* kill */
      isKilling = TRUE ;
    case 41:  /* kill a obj, single call */
      /* do NOT kill a server object
      f = filopen (journalFileName,"","a") ;
      if (f)
	fprintf(f,"\n-D %s %s\n\n",className(key), name(key)) ;
      filclose(f) ;
      */
      break ;
    case 5:  /* kill done */
      isKilling = FALSE ;
      break ;
    case 6:  /* aliasing */
      f = filopen (journalFileName,"","a") ;
      if (f)
	fprintf(f,"\n-A %s %s %s\n\n",className(key), oldaliasname, name(key)) ;
      filclose(f) ;
    case 7:  /* renaming */ /* needed in the temp-gene fixing idea, may confilct with otehr clients */
      f = filopen (journalFileName,"","a") ;
      if (f)
	fprintf(f,"\n-R %s %s %s\n\n",className(key), oldaliasname, name(key)) ;
      filclose(f) ;
    default:
      break ;
    }
}

static void  xClientAlias (KEY key, char *oldname, BOOL keepOldName) 
{
  oldaliasname = oldname ;
  xClientSave (key, keepOldName ? 6 : 7 ) ;
}

/***************************************************************************/
/***************************************************************************/
