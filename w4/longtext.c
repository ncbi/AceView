/*  File: longtext.c
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
 * Last edited: Mar 11 15:11 1999 (fw)
 * * Oct 15 15:18 1992 (mieg): Editing 
 * Created: Fri Mar 13 15:06:53 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: longtext.c,v 1.11 2016/02/04 16:56:16 mieg Exp $ */

#include "acedb.h"
#include "a.h"

#include "lex.h"
#include <ctype.h>
#include "sysclass.h"
#include "classes.h"
#include "freeout.h"
#include "pick.h"  /* needed in keysetdump */
#include "java.h"
#include "query.h"
#include "longtext.h"

/*****************************************************************/
/******************     non graphic routines *********************/
/*****************************************************************/

/* we could compress on the fly in dump, 
   but this is not longtetx specific 

  FILE *pipe ;
  extern FILE* popen() ;
  a = arrayCreate(stackMark(s), char) ;
  cp = cp0 = arrayp(a,0,char) ;
  stackCursor(s,0) ;
  pipe = popen ("compress", "rw") ;
  if (pipe)
    { register int i ; 
      while (cp = stackNextText(s))
	fprintf(pipe, "%s\n", cp) ;
      fprintf(pipe, "\n***LongTextEnd***\n") ;
      while((i = getc(pipe)) != EOF)
	*cp++ = i ;
      pclose(pipe) ;
      arrayMax(a) = cp - cp0 ;
      arrayStore(key, a,"c") ;
      stackDestroy(s) ;
      arrayDestroy(a) ;
      return TRUE ;
    }
 */    

BOOL longTextDump(FILE *f, Stack buffer, KEY k)
{ char  *cp ;
  Stack s = stackGet(k);

  if (s)
    {
      stackCursor(s,0) ;
      if (f)
	{ while ((cp = stackNextText(s)))
	  fprintf(f, "%s\n", cp) ;
	}
      else if(stackExists(buffer))
	{ while ((cp = stackNextText(s)))
	  { catText(buffer, cp) ;
	  catText(buffer,"\n") ;
	  }
	}
      else
	{ while ((cp = stackNextText(s)))
	  { freeOut (cp) ;
	  freeOut ("\n") ;
	  }
	}
    }
      

  if (f)
    fprintf(f, "***LongTextEnd***\n") ;
  else if(stackExists(buffer))
    catText(buffer, "***LongTextEnd***\n") ;
  else
    freeOut ("***LongTextEnd***\n") ;

  stackDestroy(s) ;
  return TRUE ;
}  

/************************************/

BOOL javaDumpLongText(KEY key) 
{
  char *cp;
  Stack s = stackGet(key);

  if (!s)
    return FALSE ;
  stackCursor(s,0) ;

  freeOut(messprintf("?LongText?%s?\t?txt?",freejavaprotect(name(key))));
  while ((cp = stackNextText(s)) )
    { freeOut(freejavaprotect(cp));
      freeOut("\\n");
    }
  freeOut("?\n\n");
  stackDestroy (s) ;
  return TRUE ;
}

/************************************/

BOOL bsMinilibDumpLongText(FILE *f, Stack buffer, KEY k)
{ char  *cp ;
  Stack s = stackGet(k);

  if (s)
    {
      stackCursor(s,0) ;
      if (f)
	{
	  while ((cp = stackNextText(s)))
	  fprintf(f, " %s\n", cp) ;
	}
      else if(stackExists(buffer))
	{ while ((cp = stackNextText(s)))
	  { 
	  catText(buffer," ");
	  catText(buffer, cp) ;
	  catText(buffer,"\n") ;
	  }
	}
      else
	{ while ((cp = stackNextText(s)))
	  { 
	  freeOut (" ");
	  freeOut (cp) ;
	  freeOut ("\n") ;
	  }
	}
    }

  stackDestroy(s) ;
  return TRUE ;
}  


/************************************/

BOOL longTextParse (int level, KEY key)
{
  char *cp = 0 ;
  Stack s = stackCreate(500) ;
  BOOL foundEnd = FALSE ;

#ifdef ACEDB4
  stackTextOnly(s); /* since we save it. */
#endif

  freespecial ("\n\\") ;
  while ((cp = freecard(level)))
    if (!strcmp(cp,"***LongTextEnd***"))
      { foundEnd = TRUE ;
	break ; 
      }
    else
      pushText(s, cp) ;

  if (!foundEnd)
    { messerror ("no ***LongTextEnd*** found when parsing LongText %s",
		 name(key)) ;
      return FALSE ;
    }
 
  if (stackMark(s))
    stackStore(key,s) ;

  stackDestroy(s) ;
  return TRUE ;
}

/************************************/
/* filters on template or AND's all the members of the stack */
BOOL longTextGrep (KEY key, char *template, Stack sWhat)
{
  char *cp = 0 ;
  int nn ;
  Stack s ;
  static Stack s1 = 0 ;
  BOOL found = FALSE ;

  if (class(key) != _VLongText &&
      class(key) != _VFicheView)
    return FALSE ;
      
  s = stackGet(key) ;
  if (!stackExists(s))
    return FALSE ;
  
  s1 = stackReCreate(s1, stackMark(s)) ;

  stackCursor(s, 0) ;
  while ((cp = stackNextText(s)))
    { catText(s1, cp) ;
      catText(s1, " ") ;
    }
  
  if (template)
    found = pickMatch(stackText(s1,0), template) ;
  else if (sWhat)
    {
      found = TRUE ; 
      stackCursor(sWhat,0) ;
      while (found && 
	     (template = stackNextText (sWhat)))
	{
	  nn = pickMatch(stackText(s1,0), template) ;
	  if (!nn)
	    found = FALSE ;
	}
    }
  stackDestroy(s) ;
  return found ;
}

/************************************/

static KEYSET longGrepDo (char *template, BOOL reMap)
{ 
  KEY key = 0 , key1 ;
  KEYSET ks = keySetCreate() ;
  int n = 0 ;
  
  while(lexNext(_VLongText, &key))
    { 
      if (longTextGrep(key, template, 0))
	{
	  BOOL ok = FALSE ;
	  if (reMap)
	    {  
	      int c ; /* c = _VPaper ; */
	      for (c = 0 ; c < 256 ; c++)
		{
		  key1 = key ;
		  if (lexReClass(key, &key1, c))
		    {
		      ok = TRUE ;
		      keySet(ks, n++) = key1 ;
		    }
		}
	    }
	  if (! ok)
	    keySet(ks, n++) = key ;
	}
    }
  keySetSort (ks) ;  /* needed because of ReClass calls */
  keySetCompress(ks) ;
  return ks ;
}


int longTextFilter (KEYSET ks, Stack s, BOOL reMap)
{ 
  KEY key, key1 ;
  int i, j ;
  
  for (i = j = 0 ; i < keySetMax (ks) ; i++)
    {
      key = keySet (ks, i) ;
      if (longTextGrep(key, 0, s))
	{
	  if (!reMap || !lexReClass(key, &key1, _VPaper))
	    key1 = key ;
	  keySet (ks, j++) = key1 ;
	}
    }
  keySetMax (ks) = j ;
  return j ;
}


KEYSET longGrep (char *template)
{
  return longGrepDo (template, TRUE) ;
}

KEYSET longGrepPlain (char *template)
{
  KEY paper = 0, abstract ;
  KEYSET ks = keySetCreate() ;
  int n = 0 ;
  KEY _Gene = str2tag ("Gene") ;
  KEY _Abstract = str2tag ("Abstract") ;
  
  
  while(lexNext(_VPaper, &paper))
    { 
      if (!keyFindTag (paper, _Gene) || 
	  ! (abstract = keyGetKey (paper, _Abstract)) ||
	  strcmp (name(paper), name(abstract)))
	continue ;
      if (longTextGrep(abstract, template, 0))
	keySet(ks, n++) = abstract ;
    }
  keySetSort (ks) ;  /* needed because of ReClass calls */
  keySetCompress(ks) ;
  return ks ;
}

/*****************************************************************/
/****************** end non graphic routines *********************/
/*****************************************************************/

#ifndef NON_GRAPHIC

#include "display.h"

typedef struct LOOKSTUFF
  { int   magic;        /* == MAGIC */
    KEY   key ;
    Graph graph ;
    int   activeBox,  editing ;
    Stack textStack ;
    int myTextStack ;
    KEYSET box2key ;
    Array segs,      	/* array of SEG's from the obj */
          boxIndex ;    /* if >256 a SEG, else something special */
  } *LOOK ;

typedef struct
  { KEY key ;
    float x, dx ;
    unsigned int flag ; /* I use shift right operators */
  } SEG ;
#define segFormat "k2fi"

static void longTextDestroy (void) ;
static void longTextRecover (void) ;
#ifdef ACEDB_MAILER
static void longTextMail (void) ;
#endif
static void longTextRead (void) ;
static void longTextDecorate (void) ;
BOOL longTextParse (int level, KEY key) ;
BOOL longTextDisplay (KEY key, KEY from, BOOL isOldGraph) ;
static void generalGrep(void) ;

static MENUOPT longTextMenu[] = 
            { {graphDestroy, "Quit"},
		{help,"Help"},
		{graphPrint,"Print"},
		{displayPreserve,"Preserve"},
		{longTextDecorate, "Search"},
#if defined(JUNK)
#if !defined(WIN32)
		{longTextEditor,"Edit"},  /* used to work, but no longer */
#endif
#endif
		{longTextRecover,"Save"},
#ifdef ACEDB_MAILER
		{longTextMail,"Mail"},
#endif
		{longTextRead,"Read"},
		{generalGrep,"General Grep"},
		{0, 0}
            } ;

static int MAGIC = 464565 ;	/* also use address as graphAss handle */
#define LOOKGET(name)     LOOK look ; \
                          if (!graphAssFind (&MAGIC, &look)) \
		            messcrash ("graph not found in %s",name) ; \
                          if (look->magic != MAGIC) \
                            messcrash ("%s received a wrong pointer",name)


/**************************************************************/

static void longTextPick(int box)
{
  LOOKGET("longTextPick") ;

  if (!box || ! keySetExists(look->box2key) || box >= keySetMax(look->box2key))
    return ;
  
  if (box == look->activeBox)         /* a second hit - follow it */
    display (keySet(look->box2key, box), 0, 0) ;
  else
    { if (look->activeBox)
	graphBoxDraw(box, BLACK, YELLOW) ;
      graphBoxDraw(box, BLACK, RED) ;
    }

  look->activeBox = box ;
}

/************************************/

static void longTextRead (void)
{
  char *cp = 0 ; int level ;
  Stack s ; 
  FILE *f ;

  LOOKGET("longTextRead") ;

  if (!(f = filqueryopen (0, 0, "ace", "w", "File to export the results")))
    return ;

  if (look->editing &&
      messQuery ("You did not save the edited form of the former text, should I")
      )
    longTextRecover () ;    
  
  look->editing = FALSE ;
  s = stackCreate(500) ;
  level = freesetfile(f,"") ;
  freespecial ("\n\t") ;
  while ((cp = freecard(level)))
    { if (!strcmp(cp,"***LongTextEnd***"))
	break ; 
      else
	pushText(s, cp) ;
    }
  stackDestroy(look->textStack) ;
  look->textStack = s ;
}

/**************************************************************/

static void longTextDecorate (void)
{ int x = 1 , y = 4, maxx = 0, box ;
  char *cp , c ;
  Stack s ; 
  KEY key ;
  int myClassList[5] ;
  LOOKGET("longTextDecorate") ;

  myClassList[0] = _VLocus ; /* SGI need constants for at declaration */
  myClassList[1] = _VGene ; /* SGI need constants for at declaration */
  myClassList[2] = _VSequence ; /* SGI need constants for at declaration */
  myClassList[3] = _VAllele ; /* SGI need constants for at declaration */
  myClassList[4] = -1 ; /* SGI need constants for at declaration */

  
  messStatus("Searching the text") ;
  graphClear() ;
  graphButtons (longTextMenu, 3, 1.5, 60) ;

  s = look->textStack ;
  if (stackExists(s) != look->myTextStack)
    { messout("longTextDraw lost its stack, sorry") ;
      return ;
    }
  
  look->box2key = keySetReCreate(look->box2key) ;
  stackCursor(s, 0) ;
  while ((cp = stackNextText(s)))
    { int level = freesettext(cp,"") ;
      freespecial("") ;
      freecard(level) ;
      while (*cp && *cp == ' ')	/* maintain indenting after decorate */
	{ ++x ; ++cp ; }
      while((cp = freewordcut(" .,/;:[]{}()!@#$%^&*+=\\|\"\'~",&c)) || c)
	{
	  if (cp)
	    { 
	      int *classp ;
	      
	      for (classp = myClassList ; *classp != -1 ; classp++)
		if (lexword2key(cp, &key, *classp++))
		  { keySet(look->box2key, box = graphBoxStart()) = key ;
		    graphText (cp, x, y) ;
		    graphBoxEnd() ;
		    graphBoxDraw(box, BLACK, YELLOW) ;
		    goto done ;
		  }

	      graphText (cp, x, y) ;
	    done:
	      x += strlen(cp) ;
	      
	      if (x > maxx)
		maxx = x ;
	    }
	  if (c)
	    graphText(messprintf("%c", c), x++, y) ;
	}
      x = 1 ;
      y++;
    }
  graphTextBounds( maxx + 2 , y+3) ;
  graphRedraw() ;
}

/**************************************************************/

static void longTextDraw (LOOK look, BOOL reuse)
{ int x = 1 , y = 4, maxx = 0  ;
  char *cp ;
  Stack s = look->textStack ;
  
  graphClear() ;
  graphButtons (longTextMenu, 3, 1.5, 60) ;

  if (stackExists(s) != look->myTextStack)
    { messout("longTextDraw lost its stack, sorry") ;
      return ;
    }
  stackCursor(s, 0) ;
  while ((cp = stackNextText(s)))
    {
      graphText (cp, x, y) ;
      x += strlen(cp) + 1 ;
      if (x > maxx)
	maxx = x ;
      x = 1 ;
      y++;
    }
  graphTextBounds( maxx + 2 , y+3) ;
  graphRedraw() ;
}

/**********************************/

#if defined(JUNK)
#if !defined (WIN32)
static void longTextEditor (void) 
{ FILE *f = 0 ;
  static int n = 1, ok = 12 ;
  Stack  s ;
  char *cp ;
  LOOKGET("longTextEditor") ;
  
  s = look->textStack ;
  if (!stackExists(s))
    return ;
  while((filName(messprintf("/tmp/acedb.editor.%d",n),"","r")) ||
	(filName(messprintf("/tmp/acedb.editor.%d.done",n),"","r")))
    n++ ;
  f = filopen(messprintf("/tmp/acedb.editor.%d",n),"","w") ;
  
  if (!f)
    { messprintf("Sorry, i cannot create the file %s",
		  messprintf("/tmp/acedb.editor.%d",n),"w") ;
      return ;
    }

  stackCursor(s, 0) ;
  while ((cp = stackNextText(s)))
    fprintf(f, "%s\n", cp) ;
  filclose(f) ;
 
  if ((ok = callScript (messprintf("acedb_editor /tmp/acedb.editor.%d &",n))) != -1)
    look->editing = n++ ;
  displayPreserve() ;
}
#endif
#endif

/**********************************/

#ifdef ACEDB_MAILER
static void longTextMail (void) 
{ void *v ; KEYSET ks ;
  LOOKGET("longTextMail") ;
 
  if (!stackExists(look->textStack))
    return ;
  
  if (!keySetActive(&ks, &v))
    { messout ("First select a keySet containing the recipients") ;
      return ;
    }
  acedbMailer(0, ks, look->textStack) ;
}
#endif

/**********************************/

static void longTextRecover (void) 
{ FILE *f ;
  int level, nn = 0 ;
  Stack s ;
  char *cp ;
  LOOKGET("longTextRecover") ;
  
  if (!look->editing)
    return ;
 lao:
  f = filopen(messprintf("/tmp/acedb.editor.%d.done",look->editing),"","r") ;
  
  if (!f)
    { if (!messQuery (
"To store the edited text, first exit the text editor.\n Are you ready ?"))
	return ;
      else
	if (nn++ > 1 &&
	    messQuery("do you want to abandon the edition"))
	  goto done ;
	goto lao ;
    }

  level = freesetfile(f,"") ;
  freespecial ("\n\t") ;
  if (longTextParse(level, look->key) &&
      (s = stackGet(look->key)))
      { stackClear(look->textStack) ;
	stackCursor(s, 0) ;
	while ((cp = stackNextText(s)))
	  pushText(look->textStack,cp) ;
	stackDestroy(s) ;
      }
  else
    messout("I failed to recover the file /tmp/acedb.editor.%d.done",
	    look->editing) ;
 
done:
  look->editing = 0 ;
  longTextDraw (look, 0) ;
}

/**********************************/

BOOL longTextDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  LOOK look=(LOOK)messalloc(sizeof(struct LOOKSTUFF)) ;
  KEY  curr = 0 ;
  
  look->magic = MAGIC;

  if (key && class(key) != _VLongText)
    goto abort ;

  if (! (look->textStack = stackGet(key)))
    goto abort ;

  look->key = key ;

  look->boxIndex = arrayCreate (64,SEG*) ;
  look->activeBox = 0 ;
  keySetDestroy(look->box2key) ;

  if (isOldGraph)
    { 
      longTextDestroy () ;
      graphClear () ;
      graphGoto (0,0) ;
      graphRetitle (name (key)) ;
      graphAssRemove (&MAGIC) ;
    }
  else 
    { if (!displayCreate(DtLongText)) 
	goto abort ;
      
      graphRetitle (name(key)) ;
      graphRegister (DESTROY, longTextDestroy) ;
      graphRegister (PICK, longTextPick) ;
/*
      graphRegister (KEYBOARD, longTextKbd) ;
      graphRegister (MIDDLE_DOWN, longTextMiddleDown) ;
 */
      
      graphMenu (longTextMenu) ;
    }

  look->graph = graphActive() ;
  look->myTextStack = stackExists(look->textStack) ;
  graphAssociate (&MAGIC, look) ;
  longTextDraw (look, curr) ;

  return TRUE ;

abort :
  messfree (look) ;
  return FALSE ;
}

/************************************************************/
/***************** Registered routines *********************/

static void longTextDestroy (void)
{
  LOOKGET("longTextDestroy") ;

  if (look->editing &&
      messQuery ("You did not save the edited form of this text, should I")
      )
    longTextRecover () ;    
  
  arrayDestroy(look->segs) ;
  arrayDestroy(look->boxIndex) ;
  stackDestroy(look->textStack) ;
  keySetDestroy(look->box2key) ;

  look->magic = 0 ;
  messfree (look) ;
  
  graphAssRemove (&MAGIC) ;
}

/**************************************************************/
/**************************************************************/


   /* general LongText Analyser */

static void generalGrep(void)
{
  KEY  _Abstract, *ppp, key, *kp, ltk;
  char *cp, cutter;
  Stack s ;
  OBJ obj = 0 ;
  int i, level ,  nX = 0 , nLT = 0 , xn ;
  KEYSET xs = 0 , ks=0, ks1=0 ;
  static char dName[DIR_BUFFER_SIZE], filName[FIL_BUFFER_SIZE] ;
  FILE *g, *f ; 
  int myClassList[5] ;
  void *dummy ;

  myClassList[0] = _VLocus ; /* SGI need constants for at declaration */
  myClassList[1] = _VGene ; /* SGI need constants for at declaration */
  myClassList[2] = _VSequence ; /* SGI need constants for at declaration */
  myClassList[3] = _VAllele ; /* SGI need constants for at declaration */
  myClassList[4] = -1 ; /* SGI need constants for at declaration */

  if (!lexword2key("Abstract",&_Abstract,0))
    { messout ("Tag Abstract is absent from this database, nothing I can report, sorry") ;
    return ;
    }
  if (!keySetActive(&ks1, &dummy))
    { messout ("First select a keyset containing Papers with abstracts") ;
      return ;
    }
   
  messout ("This button searches the text of the abstracts "
	   "of the papers in the active keyset for "
	   "names known in class: Gene Locus Sequence Allele ") ;

  ks = query(ks1,"CLASS Paper ; Abstract") ;
  if (!keySetMax(ks))
    { messout ("First select a keyset containing papers with abstracts") ;
      keySetDestroy (ks) ;
      return ;
    }
  f = filqueryopen(dName, filName, "ace", "w","I need a place to export the Grep results") ;
  g = filopen(filName, "bad","w") ;
  if (!f || !g)
    { messout("Got a cancel file, i quit") ;
      if (f) filclose(f) ;
      if (g) filclose(g) ;
      keySetDestroy (ks) ;
      return ;
    }
  
  i = keySetMax(ks) ;
  ppp = arrp(ks, 0, KEY) - 1 ;
  while (ppp++, i--)
    { 
      obj = bsCreate(*ppp) ;
      if (!obj) continue ;
      ltk = 0 ;
      bsGetKey (obj,_Abstract,&ltk) ;
      bsDestroy (obj) ;
      if (!ltk) continue ;

      s = stackGet(ltk) ;
      if (!stackExists(s))
	continue ;
      nLT++ ;
      stackCursor(s,0) ;
      xs = keySetReCreate(xs) ; xn = 0 ;

      while (s && (cp = stackNextText(s)))
	{ level = freesettext(cp,"") ;
	  	  
	  freespecial("") ;
	  freecard(level) ;
	  	  
	  while((cp = freewordcut(" .,/;:[]{}()!@#$%^&*+=\\|\"\'~\n",&cutter)) || cutter)
	    if (cp)
	      { 
		int *classp ;
		
		for (classp = myClassList ; *classp != -1 ; classp++)
		  if (*classp && lexword2key(cp, &key, *classp))
		    { keySet(xs, xn++) = key ;
		      break ;
		    }
		if (islower((int)*cp) &&
		    islower((int)*(cp+1)) &&
		    islower((int)*(cp+2))&&
		    (*(cp+3) == '-') &&
		    isdigit((int)*(cp+4)) &&
		    !lexword2key(cp, &key, _VLocus))
		  fprintf(g, "Paper \"%s\"\nGene \"%s\"\n\n", name(*ppp), cp) ;
	      }
	}

      if (xn)
	{ keySetSort(xs) ;
	  keySetCompress(xs) ;
	  
	  fprintf(f, "\nPaper %s\n", name(*ppp)) ;
	  xn = keySetMax(xs) ;
	  nX += xn ;
	  kp = arrp(xs, 0, KEY) - 1 ;
	  while (kp++, xn--)
	    fprintf(f, "%s %s \n", className(*kp), name(*kp)) ;
	  
	}
      stackDestroy(s) ;
    }

  messout ("Found %d XREF in %d LongTexts \n", nX, nLT) ;

  keySetDestroy(ks) ;

  filclose(f) ;
  filclose(g) ;
  keySetDestroy(xs) ;
}

/**************************************************************/

#endif  /* NON_GRAPHIC */

/**************************************************************/
/**************************************************************/
