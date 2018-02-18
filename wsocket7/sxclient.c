/*  File: sxclient.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 2000
 *-------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Implements the client/server version of xace. The only
 *              external function is xClientInit() which is called from
 *              xclientmain.c. Other functions are set up as callbacks
 *              by xClientInit() so that xaceclient will make calls to
 *              the server rather than to its own local copy of the
 *              database when it needs new data.
 *              
 * Exported functions: xClientInit()
 * HISTORY:
 * Last edited: Sep 21 14:06 2000 (edgrif)
 * Created: Wed Apr 26 13:25:31 2000 (edgrif)
 * CVS info:   $Id: sxclient.c,v 1.2 2006/06/22 20:21:40 mieg Exp $
 *-------------------------------------------------------------------
 */

#include <wh/acedb.h>
#include <wh/aceio.h>

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
#include <wh/aceclient.h>
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

#include <wh/parse.h>
#include <wh/pick.h>
#include <whooks/sysclass.h>
#include <wh/session.h>
#include <wh/call.h>
#include <wh/dump.h>
#include <wh/lex.h>
#include <wh/bs.h>
#include <wh/client.h>
#include <wh/table.h>
#include <wh/main.h>
#include <wh/dbpath.h>
#include <wh/xclient.h>
#include <wsocket/sclientlib.h>

void *handle = 0 ;
static u_long port = DEFAULT_PORT;
static BOOL serverIsAllowed  = FALSE ; 
static void xClientCommit (void) ;
static KEYSET xClientGet (KEY key, char *quer, char *grp, BOOL getNeighbours) ;
static KEYSET xClientQuery (KEY key, char *quer, char *grp, BOOL getNeighbours) ;
static void  xClientSave (KEY key, int action) ;
static void  xClientAlias (KEY key, char *oldname, BOOL keepOldName) ;
static TABLE* xClientTableMaker (char *quer) ;



static char* oldaliasname = 0 ;
static char server[255] ; 
static char *journalFileName = 0 ; 
static int clientCount = 0 , serverCount = 0 ;
static int timeOut = 300 ;

/* uuuggghhh, I hate this, but nothing else to do but go along with globals  */
/* for now...                                                                */
Client client = NULL ;
int time_out_unused = 0 ;


/***************************************************************************/
/***************************************************************************/

char *xClientServerName(void)
{
  return server ;
}

static BOOL getServerCount(void)
{
  /* asks server for its update-data number */
  BOOL ok = FALSE ;
  Stack sr = 0 ;
  int  encore = 0,erc, dummy, level, n ;
  char *cp ;
  char *reponse = 0 ;


  if (!serverIsAllowed)
    return FALSE ;
  
  if (!handle)
    xClientInit() ;
  if (!handle)
    return 0 ;

  erc = askServerBinary (handle, "Data_Version", &reponse, &dummy, &encore, 0, time_out_unused) ;

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
  freespecial("\n\\") ;
  while (freecard(level))
    {
      while ((cp = freeword()))				    /* does a freeclose */
	{
	  if (!strcmp(cp,"version"))
	    { 
	      if (!freeint(&n))
		{
		  freeclose (level) ;
		  goto fin ;
		}
	      serverCount = n ;
	      ok = TRUE ;
	    }
	}
    }
fin:
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
{
  static int firstPass = 0 ;
  char *cp ;
  int i, level ; 
  FILE *fil=NULL;
  int sock ;

  if (firstPass==0)
    { 
      firstPass = 1 ;
      port = 0 ; *server = 0 ;


      /* Set up our client state data.                                           */
      /* (if userid != NULL then we prompt for passwd)                           */
      client = (Client)messalloc(sizeof(ClientStruct)) ;
      client->userid = client->nonce = client->passwd_hash = client->nonce_hash = NULL ;


      /* Hack for now....                                                    */
      /* Somehow this will have to be passed in or obtained via dialogs...   */
      /*                                                                     */
      client->userid = "edgrif" ;
      client->passwd_hash = "38a6557860893d95c7e738c0e5cd233d" ;



      cp = dbPathStrictFilName ("wspec", "server", "wrm", "r", 0) ;
      if (!cp || !(fil = filopen (cp, 0, "r")))
	messExit ("initXClient cannot open server configuration file: %s", cp) ;
      messfree(cp);

      level = freesetfile(fil, 0) ;
      while (freecard(level))
	{
	  cp = freeword() ;
	  if (cp)
	    { 
	      if (!strcmp(cp,"EXTERNAL_SERVER") && (cp = freeword()))
		strncpy (server, cp, 254) ;
	      else if (!strcmp(cp,"EXTERNAL_PORT") && freeint (&i))
		port = i ;
	    }
	}
    }

  if (*server && port && !handle)
    handle = openServer (server, port, timeOut, &sock, client) ;
  if (!handle)
    {
      if (firstPass==1)
	messout("Connection to server %s, port %lu  failed, sorry, please check wspec/server.wrm", 
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
  char *reponse = 0 ;
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
  erc = askServerBinary (handle, stackText(sq,0), &reponse, &dummy, &encore, 0, time_out_unused) ;

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
        parseBufferInternal(cp, kk, TRUE) ;
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
  char *reponse = 0 ;
  static int nn = 0 ;
  int sock ;

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
    {
      nn++ ;
      
      /* cripes, the code needs to create another load of stuff here for client*/
      /* etc. better make client global for now until I understand this...     */
      handle = openServer (server, port, timeOut, &sock, client) ;
    }
  if (!handle)
    return FALSE ; 
  if (nn > 0) nn-- ;
  erc = askServerBinary (handle, buffer, &reponse, &dummy, &encore, 0, time_out_unused) ;

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
  int  nerr = 0, encore = 0,erc, dummy;
  ACEIN fi = 0 ;
  char *reponse = 0 ;
  char *unused ;

  if (!serverIsAllowed || !quer || !*quer)
    return 0 ;
  
  if (nerr > 20)
    { nerr = 0 ; xClientLost() ; }

  if (!handle)
    xClientInit() ;
  if (!handle)
    return 0 ;

  messStatus ("Reading from server") ;
  sq = stackReCreate (sq, 1000) ;
  sr = stackReCreate (sr, 1000) ;
  pushText (sq, quer) ;

suiteReponse:
  erc = askServerBinary (handle, stackText(sq,0), &reponse, &dummy, &encore, 0, time_out_unused) ;

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

  fi = aceInCreateFromText (stackText(sr,0), "", 0);
  aceInSpecial(fi,"\n/\\") ;
  tt = tableCreateFromAceIn (fi,
			     "(table)",			    /* Do we know the table name ?? */
			     &unused);
  aceInDestroy (fi) ;

  stackDestroy (sr) ;

  return tt ;
} /* xClientTableMaker */

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
  freespecial ("\n/\\\t@") ;
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

static void  xClientSave (KEY key, int action) 
{
  static char* tmpName1 = NULL ;
  static BOOL isKilling = FALSE ;
  char *tmpName2 = NULL ;
  FILE *f ; 
  ACEOUT dump_out;


  if (X_CLIENT_PARSING 
      || (isKilling && action != 5)
      || !dumpClassPossible (class(key),'a') )
    return ;

  if (!journalFileName) 
    {
      journalFileName = filGetName("clientUpdate","","a", 0) ;
      if (!journalFileName) 
	return ;
    }

  switch (action)
    {
    case 1:  /* save a B object, phase 1: old version */

      if (tmpName1) { filtmpremove(tmpName1) ; tmpName1 = NULL ;}
      f = filtmpopen (&tmpName1,"w") ;   
      if (!f && tmpName1)  { filtmpremove(tmpName1) ; tmpName1 = NULL ; return ; }
      filclose(f) ;		/* just interested in tmpname */

      dump_out = aceOutCreateToFile (tmpName1, "w", 0);
      if (!dump_out)  { filtmpremove(tmpName1) ; tmpName1 = NULL ; return ; }
      dumpKey(dump_out, key);
      aceOutDestroy (dump_out);
      break ;
    case 2:  /* save a B object, phase 2: new version */
      if (!tmpName1) return ;
      f = filtmpopen (&tmpName2,"w") ;   
      if (!f && tmpName2)  { filtmpremove(tmpName2) ; tmpName2 = NULL ; return ; }
      filclose(f) ;		/* just interested in tmpname */

      dump_out = aceOutCreateToFile (tmpName2, "w", 0);
      if (!dump_out)  { filtmpremove(tmpName2) ; tmpName2 = NULL ; return ; }
      dumpKey(dump_out, key);
      aceOutDestroy (dump_out);


#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
      /* I don't know what this /bin/env stuff was supposed to do, it can be */
      /* used to modify the commands environment but this call does not do   */
      /* that... it's just wierd and seems to produce path problems for me.. */
      callScript("/bin/env acediff", messprintf(" %s %s >> %s", tmpName1, tmpName2,
						journalFileName)) ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

      callScript("acediff", messprintf(" %s %s >> %s", tmpName1, tmpName2,
						journalFileName)) ;

      filtmpremove(tmpName1), tmpName1 = NULL ;
      filtmpremove(tmpName2), tmpName2 = NULL ;
      break ;
    case 3:  /* dump */
      dump_out = aceOutCreateToFile (journalFileName, "a", 0);
      if (dump_out)
	{ 
	  aceOutf (dump_out, "\n");
	  dumpKey(dump_out, key);
	  aceOutf (dump_out, "\n");
	  aceOutDestroy (dump_out);
	}
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
