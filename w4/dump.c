/*  File: dump.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov  6 16:52 1998 (rd)
 * * Aug 18 10:00 1998 (rd):
 *	     - split dumpPossible() into lexIsKeyVisible() and dumpClassPossible
 *	     - remove dumpKey1() as unnecessary
 * * Aug  3 11:29 1998 (rd): 
 *           - user can dump to any directory. Also changed the error format for
 *           - errors form dump  to class \t name \t error mess.
 * * Apr  3 16:18 1997 (rd)
 * * Jun 21 17:05 1992 (mieg): arranged  so the dump file can be read back
 * Created: Sun Feb 16 01:31:41 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: dump.c,v 1.14 2007/10/29 20:35:29 mieg Exp $ */

#include "acedb.h"

#include "lex.h"
#include "a.h"
#include "bs.h"
#include "dump.h"
#include "pick.h"
#include "session.h"
#include "sysclass.h"
#include "mytime.h"
#include "query.h"
#include "freeout.h"
#include "classes.h"
#include "java.h"
#include <setjmp.h>

#if defined(applec)
#include "MPWIncludes.h"
#endif

static Array doNotDumpKeySet = 0 ;

/******************************************************************/
/******************************************************************/

static FILE* f = 0 ;
static Stack s = 0 ;
DumpFunc dumpFunc[256] ;
static char dumpDir[DIR_BUFFER_SIZE] = "";         /* use this to choose where to dump too */

/*********************************************************/

static BOOL dumpSetDirectory (void)
{
  char fname[FIL_BUFFER_SIZE]="";
  char dname[DIR_BUFFER_SIZE]="";
  FILE *f = 0;

  f = filqueryopen (dname, fname, "", "rd", 
		    "Choose directory for dump files");

  if (!f)
    return FALSE;
  fclose (f) ;

  sprintf (dumpDir,"%s/",dname);

  return TRUE;
}


/*********************************************************/

/* RD 980817: split dumpPossible() into dumpKeyPossible() and
   dumpClassPossible().  I want dumpKeyPossible to be used for all key
   listing, e.g. in keysetdisplay, list and show commands in tace.
   dumpClassPossible will be used by the complete dump code.

   NB I allow protected classes (including XREF and Model classes) to
   be seen by dumpKeyPossible(), so we can see them in keysets and
   display them, and also list and show them from tace.  But they
   still fail dumpClassPossible() so will not be in full dumps.  This
   is a change.  Is that OK?

   RD 980819: dumpKeyPossible() is the same as lexIsKeyVisible(), not
   surprisingly, so I supressed it.
*/

BOOL dumpClassPossible (int classe, char style)
{
  classe &= 0xff ;

  if (pickList[classe].type == 'A')
    return dumpFunc[classe] ? TRUE : FALSE ;

  if (pickList[classe].private)  /* private to kernel, nobody's business  */
    return FALSE ;

  /* protected classes can not be parsed, so can't be dumped in ace style 
   *  XREF classed are created on the fly by the kernel
   */
  if ((style == 'a' || style == 'A' || style == 0) && 
      ( pickXref(classe) || pickList[classe].protected))
    return FALSE ;

  if (!lexMaxVisible(classe))
    return FALSE ;

  return TRUE ;
}

/*********************************************************/

static BOOL dumpKey2 (KEY key, char style, COND cond)
{ 
  OBJ obj ;
  int classe = 0, oldc = 0 ;
  BOOL notEmpty = FALSE ;

  if (doNotDumpKeySet && keySetFind (doNotDumpKeySet, key, 0))
    return FALSE ;

  if (!lexIsKeyVisible (key))
    return FALSE ;

		/* test for tag now, to avoid writing an object name line.  
		   cond is a compiled query built in the calling routine
		   queryFind3 uses bIndexFind() and only makes obj if necessary
		*/
  obj = 0 ;
  if (cond && !queryFindLocalise3 (cond, &obj, key))
    { if (obj) bsDestroy (obj) ;
      return FALSE ;
    }

  if (pickType(key) == 'A' && !(style == 'j' || style == 'J' || style == 'C' ))
    style = 0 ;	/* can only dump arrays in .ace, java, C */

				/* write object name line for .ace files */
  if (style == 0 || style == 'a' || style == 'A' )
    {
      notEmpty = TRUE ;
      if (f)
	fprintf (f,"%s %s\n",	   
		 className(key), freeprotect (name(key))) ;
      else if (s)
	pushText(s, messprintf ("%s %s\n",	   
				className(key), freeprotect (name(key)))) ;
      else
	freeOutf ("%s %s\n", className(key), freeprotect (name(key))) ;
    }

  switch (pickType (key))
    {
    case 'A': 
      switch (style)
	{
	case 'C':    /* C dumping */
          {
	  if (f || s) break; /* C dump knows only freeOut */
	  classe = class(key) ;
	  notEmpty = TRUE ; 
	  if (iskey(key) == 2)
	    { 
	      if (classe == _VDNA)
		{  dnaDumpKeyCstyle (key); }
	      else if (classe == _VPeptide)
		 { peptideDumpKeyCstyle (key); }
	      /*
		else if (classe == _VKeySet)
		{ notEmpty = TRUE ; javaDumpKeySet(key); }
		else if (classe == _VLongText)
		{ notEmpty = TRUE ; javaDumpLongText(key); }
	      */
	    }
	  else
	    { 
	      char buf[2] ;
	      char *cp = name(key) ;
	      
	      buf[0] = 'n' ; buf[0] = class(key) ;	/* bug? */
	      freeOutBinary (buf, 2) ;
	      freeOutBinary (cp, strlen(cp) + 1) ; 
	      buf[0] = '#' ; 
	      freeOutBinary (buf,1) ; 
	    }
	  break ;
          }
	case 'j': case 'J':   /* java dumping */
          {
	  if (f || s) break; /* Java dump knows only freeOut */
	  classe = class(key) ;
	  if (iskey(key) == 2)
	    {
	      if (classe == _VDNA)
		{ javaDumpDNA(key); }
	      else if (classe == _VPeptide)
		{ javaDumpPeptide(key); }
	      else if (classe == _VKeySet)
		{ notEmpty = TRUE ; javaDumpKeySet(key); }
	      else if (classe == _VLongText)
		{ notEmpty = TRUE ; javaDumpLongText(key); }
	    }
	  else
	    switch (style)
	      {
	      case 'j':
		notEmpty = TRUE ; 
		freeOutf ("?%s?%s?\n",className(key),
			  freejavaprotect(name(key))) ;
		break ;
	      case 'J':
		notEmpty = TRUE ; 
		if (oldc != classe)
		  { 
		    oldc = classe ;
		    freeOutf ("?%s", className (key)) ;
		  }
		else
		  freeOutf ("#", className (key)) ;	/* bug? need fmt? */
		if (iskey(key) == 2)
		  freeOutf ("?%s?\n",
			    freejavaprotect(name(key))) ;
		else
		  freeOutf("!%s?\n",
				  freejavaprotect(name(key))) ;
	      default:
		break ;
	      }
	  break ;
          }
	default:
	  if (dumpFunc[class(key)])
	    {
	      if (iskey (key) == 2)
		notEmpty = TRUE ;
	      dumpFunc[class(key)] (f, s, key) ;
	    }
	  break ;
	}
      break ;
    case 'B':
      if (obj || (obj = bsCreate (key)))
	{
	  notEmpty = TRUE ;
	  switch (style)
	    { 

	    case 0: case 'a': case 'A': /* ace format */
	      if (class(key) == _VModel)
		{ if (f) fprintf (f, "// a model can not be dumped in .ace format\n") ;
		  else if (s) pushText (s, "// a model can not be dumped in .ace format\n") ;
		  else freeOut ("// a model can not be dumped in .ace format\n") ;
		}
	      else
		bsAceDump (obj, f, s, (cond != 0) ? TRUE : FALSE) ;
	      break ;

	    default:
	      niceDump (obj, style) ;
	      break ;
	    }
	}
      else 
        {
	switch (style)
	  { 
	  case 'j':
	    notEmpty = TRUE ; 
	    freeOutf ("?%s?%s?\n",className(key),
		      freejavaprotect(name(key))) ;
	    break ;
	  case 'J':
	    notEmpty = TRUE ; 
	    if (oldc != classe)
	      { 
		oldc = classe ;
		freeOutf ("?%s", className (key)) ;
	      }
	    else
	      freeOutf ("#", className (key)) ;
	    if (iskey(key) == 2)
	      freeOutf ("?%s?\n",
			      freejavaprotect(name(key))) ;
	    else
	      freeOutf ("!%s?\n",
			      freejavaprotect(name(key))) ;
	    break ;
	  case 'C':
	    notEmpty = TRUE ; 
	    { 
	      char buf[2] ;
	      char *cp = name(key) ;

	      buf[0] = 'n' ; buf[1] = class(key) ;
	      freeOutBinary (buf, 2) ;
	      freeOutBinary (cp, strlen(cp) + 1) ;
	      buf[0] = '\n' ; 
	      buf[1] = '#' ; 
	      freeOutBinary (buf,2) ; 
	    }
	    break ;
          case 'h':
            {
            notEmpty = TRUE ;
	    freeOut (messprintf ("%s %s\n", className (key), name(key)));
            break ;
            }
	  default:
	    
	    break ;
	  }
        }
      break ;
    }
  if (obj) bsDestroy (obj) ;

		/* always leave a blank line after dumped object */
  if (notEmpty)
    {
      if (f)
	fputc ('\n',f) ;
      else if (s)
	pushText (s,"") ;
      else
	freeOut ("\n") ;
    }

  return TRUE ;
}

/*********************************************/
/************* public routine ****************/
/*********************************************/

/* command.c interface, dumps into freeOut */

BOOL dumpKeyBeautifully (KEY key, char style, COND cond) 
{
  f = 0 ; s = 0 ;
  return dumpKey2 (key, style, cond) ;
}

/* explicit interface, dumps where stated, 0,0 means freeOut */
BOOL dumpKey (KEY key, FILE* ff, Stack ss) 
{
  f = ff ; s = ss ;
  return dumpKey2 (key, 0, 0) ;
}


/*********************************************/
/*********************************************/

#ifndef NON_GRAPHIC

#include "graph.h"
#include "display.h"

static MENUOPT dumpMenu[]=
{
  { graphDestroy,"Quit"},
  { help,"Help"},
  { 0,0}
};

/******************************************/

static Graph dumpGraph = 0;

/******************************************/

#endif /* NON_GRAPHIC */

static Array trapKeys = 0, trapMessages ;

static BOOL dumpClass (int t, int *nk, int *no, int *ne, KEY *keyp) 
{
  FILE *f1 = f ;
  jmp_buf jmpBuf ;
  jmp_buf *jmpBufp = &jmpBuf ; /* mieg , pb on alpha compiler */
  static jmp_buf *oldJmpBuf ;	/* static to avoid volatile problems */
	/* NB following must be volatile to ensure there when jump back */
  volatile KEY kVolatile = 0 ;
  int nn = 0 ; /* number of bytes dumped */
  *nk = *no = *ne = 0 ;

  s = stackReCreate (s, 100000) ;
  f1 = f ; f = 0 ; /* forces dump on the stack */
  while (lexNext(t,keyp))
    { 
      (*nk)++ ;
      
      stackClear (s) ;
      kVolatile = *keyp ;
      if (!setjmp (jmpBuf))
	{ oldJmpBuf = messCatchCrash (jmpBufp) ;
	  if (dumpKey2 (*keyp,0,0))
	    (*no)++ ;
	}
      else			/* an error caught in messCrash() */
	{ if (!trapKeys)
	    { trapKeys = arrayCreate(32,KEY);
	      trapMessages = arrayCreate(32,char*);
	    }
	  (*ne)++ ;
	  array(trapKeys,arrayMax(trapKeys),KEY) = kVolatile ;
	  array(trapMessages,arrayMax(trapMessages),char*) = 
	    strnew (messCaughtMessage(), 0) ;
	}

      messCatchCrash (oldJmpBuf) ; /* restore value for enclosing function */
      nn += stackMark (s) ;
      fprintf (f1, "%s\n", stackText (s, 0)) ;
      if (nn > 1<<30)
	break ;
    }
  stackDestroy (s) ;
  f = f1 ;
  return (nn > 1<<30) ? TRUE : FALSE ; /* class is finished */
}

/******************************************/
#ifndef NON_GRAPHIC

static void dumpDestroy(void)
{
  dumpGraph = 0 ;
}

#endif /* NON_GRAPHIC */
/******************************************/

static FILE* dumpOpenOutput(BOOL newDump, char *classname, int nn)
{ 
  char name[DIR_BUFFER_SIZE+FIL_BUFFER_SIZE] ;
  static char letter ;
  static char date[12] ;
  char timeBuf[25] ;

  if (newDump)
    { letter = 'A' ;
      strcpy (date, timeShow (timeParse ("today"), timeBuf, 25)) ;

      if (classname)
	sprintf (name, "%sdump_%s_%c_%s.%d", dumpDir, date, letter, classname, nn) ;
      else
	sprintf (name, "%sdump_%s_%c.%d", dumpDir, date, letter, nn) ;
      while (filName (name, "ace", "r"))
	{ ++letter ;
	  if (classname)
	    sprintf (name, "%sdump_%s_%c_%s.%d", dumpDir, date, letter, classname, nn) ;
	  else
	    sprintf (name, "%sdump_%s_%c.%d", dumpDir, date, letter, nn) ;
	}
    }
  else
    if (classname)
      sprintf (name, "%sdump_%s_%c_%s.%d", dumpDir, date, letter, classname, nn) ;
    else
      sprintf (name, "%sdump_%s_%c.%d", dumpDir,date, letter, nn);
	
#ifdef NON_GRAPHIC
  if (!(f = filopen(name,"ace","w")))
    freeOut("Dump could not open the output file, Sorry");
#else
  if (!(f = filopen(name,"ace","w")))
    graphText ("Dump could not open the output file, Sorry", 5, 5) ;
  else
    if (newDump && !classname)
      graphText (messprintf ("Ready to dump into : %s", name), 2, 3) ;
  graphRedraw () ;
#endif /* NON_GRAPHIC */
  return f ;
}

/********************************************/

static void dump1Check (void)
{
  FILE *f=0;
  int i;

  if (trapKeys)
    {
      messout ("%d errors occurred during dump", arrayMax(trapKeys)) ;
      messdump ("DUMP ERROR: %d errors on dumping\n", arrayMax(trapKeys)) ;
#ifndef NON_GRAPHIC
      f = filqueryopen (0, 0, "err", "w", 
			"Choose file to dump the errors to") ;
      if (!f)
	messout ("Unable to open error file. Therefore writing to log.wrm") ;
#endif      
      if (f)
	fprintf (f, "%d errors on dumping\n", arrayMax(trapKeys)) ;

      for (i = 0 ; i < arrayMax(trapKeys) ; i++)
	{ 
	  if (f)
	    fprintf (f,"%s\t%s\t%s\n",
		     className(arr(trapKeys,i,KEY)),
		     name(arr(trapKeys,i,KEY)),
		     arr(trapMessages, i, char*)) ;
	  else
	    messdump ("%s\t%s\t%s\n",
		      className(arr(trapKeys,i,KEY)),
		      name(arr(trapKeys,i,KEY)),
		      arr(trapMessages, i, char*)) ;

	  messfree (arr(trapMessages, i, char*)) ;
	}

#ifdef NON_GRAPHIC
      messout ("Errors written to log.wrm");
#endif
      arrayDestroy (trapKeys) ;
      arrayDestroy (trapMessages) ;
      fclose (f) ;
    }
}

void dumpAll(void)
{ 
  dumpTimeStamps = messQuery ("Do you want to dump timestamps?") ;
  dumpComments = messQuery ("Do you want to dump comments?") ;
  if (!dumpSetDirectory ())
    return;
  dumpAllNonInteractive (0, FALSE) ;
}

void dumpAllNonInteractive (char *dir, BOOL split)
{
  int nk, ne, no , t , total = 0, nfile = 0 ;
  KEY cKey, last ;
#ifndef NON_GRAPHIC
  int box, line = 0 ;

  if(!graphActivate(dumpGraph))
    { dumpGraph =  displayCreate(DtDump) ;
      graphMenu (dumpMenu);
      graphTextBounds (80,100) ;
      graphRegister(DESTROY,dumpDestroy) ;
      graphColor (BLACK) ;
    }
  else
    { graphPop() ;
      graphClear();
    }
#endif /* NON_GRAPHIC */

  if (dir)
    strcpy (dumpDir, dir) ;

  if (!split)
    { f = dumpOpenOutput(TRUE, 0, ++nfile);
      if (!f) 
      goto end;
    }
  else 
    f = 0;

#ifndef NON_GRAPHIC
  
  graphText(messprintf("Session %d", thisSession.session),
	    4, 2) ;
  graphText ("Any object in the keyset called DoNotDump will be jumped", 
	     4, 6) ;
  graphRedraw () ;
  if (!messQuery("Should I proceed ?"))
    return ;

  line = 8 ;
  
  graphText ("Dumping classes ", 4, line++)  ;
  graphRedraw () ;
#endif /* NON_GRAPHIC */ 
  
  if (lexword2key ("DoNotDump",&cKey,_VKeySet))
    doNotDumpKeySet  = arrayGet(cKey, KEY, "k") ;
  for (t = 5; t<256 ; t++)
    if (lexMax(t) &&
	dumpClassPossible(t, 'a'))
      {
	BOOL more = TRUE ;
	last = 0 ;

	while (more)
	  {
	    more = FALSE ; 
	    if (split || total > 20000)
	      { 
		BOOL first;
		if (split) 
		  nfile = 0;
		if (f) 
		  { 
		    fprintf(f,"\n\n // End of this dump file \n\n") ;
		    filclose(f);
		    total = 0;
		    first = FALSE;
		  }
		else
		  first = TRUE;
#ifndef NON_GRAPHIC
		graphText("Starting a new dump file", 8, line += 3) ;
#endif
		f = dumpOpenOutput(first, split ? pickClass2Word(t) : 0, ++nfile) ;
		if (!f) 
		  break;
	      }
	    fprintf (f, "\n\n // Class %s \n\n", pickClass2Word(t))  ;
	    /* some classes may be larger than 2GB, we should measure the exported
	       size then dumpclass would say more = TRUE and we would loop
               this code is not yet written, mieg feb 2003 */ 
	    nk = no = ne = 0 ;
	    while (dumpClass (t, &nk, &no, &ne, &last))
	      {
		if (f) 
		  { 
		    fprintf(f,"\n\n // End of this dump file \n\n") ;
		    filclose(f);
		    total = 0;
		  }
		f = dumpOpenOutput(FALSE, split ? pickClass2Word(t) : 0, ++nfile) ;
		if (!f) 
		  break;
	      }
	    total += no ;
#ifndef NON_GRAPHIC
	    box = graphBoxStart() ;
	    graphText (messprintf ("Class %s : %d keys, %d obj, %d errors",
				   pickClass2Word(t),
				   nk,no, ne) , 4, line++)  ;
	    graphBoxEnd () ;
	    graphBoxDraw(box,BLACK,WHITE) ;
#endif /* NON_GRAPHIC */
	  }
      }
 end:
  if (f)
    { fprintf(f,"\n\n // End of this dump file \n\n") ;
      filclose(f) ;
    }
  else
#ifdef NON_GRAPHIC
    freeOut("Failed to open dump file, aborting dump.\n");
#else
    graphText("Failed to open dump file, cannot complete dump.", 8, line);
#endif 
  keySetDestroy (doNotDumpKeySet) ;
#ifndef NON_GRAPHIC
  graphRedraw() ;
#endif

  dump1Check();
  dumpTimeStamps = FALSE ;
  dumpDir[0] = '\0';
}

/******************************************/
/******************************************/
/**********************************************************/
/* dump/parse binary buffer as text, finished by empty line */

/* encode aceBinary array into ASCII by mapping to '_' + half-byte */
void aceBinaryDoDump (unsigned const char *cp, int n)
{
  register int i, j ;
  unsigned char buffer [62] ;

  buffer[60] = '\n' ;
  buffer[61] = 0 ;

  for (i = j = 0 ; i < n ; cp++, i++)
    { 
      buffer[j++] = (((*cp & 0xf0) >> 4) & 0xf) + '_' ;
      buffer[j++] = ((*cp & 0xf)  & 0xf) + '_' ;
      if (j == 60)
	{
	  freeOut ((char*)buffer) ;
	  j = 0 ;
	}
    }
  /* dump remnant */
  buffer[j++] = '\n' ;
  buffer[j++] = 0 ;
  freeOut ((char*)buffer) ;
  /* dump empty line */
  freeOut ("\n") ;
 
  return;
} /* aceBinaryDoDump */

/*****************************/

BOOL aceBinaryDump (FILE* f, Stack buf, KEY k) 
{
  Array a = 0 ;
  int level = 0 ;

  if (! (a = arrayGet(k, unsigned char, "c"))) 
    return FALSE ;
 
  if (f)
    level = freeOutSetFile (f) ;
  else if (buf)
    level = freeOutSetStack (buf) ;

  aceBinaryDoDump(arrp (a, 0, unsigned char), arrayMax(a)) ;

  if (level)
    freeOutClose(level) ;
  arrayDestroy(a);

  return TRUE ;
} /* aceBinaryDump */

/************************************/

BOOL aceBinaryParse (int level, KEY key)
{
  char *cp ;
  unsigned char *cq ;
  Array a = 0 ;
  int N = 1000000, n = 0 ; 
  BOOL ok = FALSE ; /* safe bet */

  a = arrayCreate (N, unsigned char) ;
  while (freecard(level) && (cp = freeword()))
    {
      while (*cp && *(cp+1))
	{
	  cq = arrayp (a, n++, unsigned char) ;
	  *cq = (((*cp - '_') & 0xf) << 4) | ((*(cp+1) - '_') & 0xf) ;
	  cp += 2 ;
	}
    }
  if (n)
    {
      arrayStore (key, a, "c") ;
      ok = TRUE ;
    }
  arrayDestroy(a) ;
  return ok ;
} /* aceBinaryParse */

/**********************************************************/
/**********************************************************/
/**********************************************************/
