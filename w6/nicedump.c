/*  File: nicedump.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 *  Lincoln Stein added java dump
 *  Doug Bigwood and John Barnett changes to handle html format
 * Last edited: Nov 23 16:24 1998 (fw)
 * * Jul 16 00:11 1998 (rd): lots of cleaning up
 	new timestamps for Java (never used in human or perl)
	down looping not recursion
	perl layout simplified (no freeOutxy) and IL \n "fix" undone
		(it's bug was because of bad interaction with freeOutxy)
	removed dumpQueries static because it was always set TRUE -
		Java always dumps ATTACHes to recursion level one,
		human and perl never do
	isModel set correctly for NEW_MODELS
 * * Dec  8 13:53 1994 (mieg): 

 * Created: Thu Jan 28 13:50:53 1993 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: nicedump.c,v 1.10 2014/10/08 14:33:12 mieg Exp $ */

#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct bblock *BBLOCK ;
typedef BBLOCK BP ;
#define DEFINE_OBJ
typedef struct sobj *OBJ ;

#include "acedb.h"
#include "bs_.h" 
#include "dump.h"
#include "java.h"
#include "model.h"		/* for isModel() */

#define LINE_LENGTH 120

static BOOL isModele ;
extern BOOL dumpComments ;
extern BOOL dumpTimeStamps ;

/********************************/
    
static char* niceTextUsual (BS bs, char *timeBuf)
{
  if (isModele)
    return name (bs->key) ;

  if (bs->key <= _LastC)
    return bsText(bs) ;
  else if (bs->key == _Int)
    return messprintf ("%d", bs->n.i) ;
  else if (bs->key == _Float)
    return messprintf ("%1.7g", bs->n.f) ;
  else if (bs->key == _DateType)
    return timeShow (bs->n.time, timeBuf, 25) ;
  else
    return name (bs->key) ;
}

/*************/

static char *niceExpandText (KEY key, char *timeBuf)
{
  KEY tag = pickExpandTag(key) ;
  OBJ obj ;
  char *text = 0 ;

  if (tag && (obj = bsCreate(key)))
    { if (bsGetKeyTags (obj,tag,0))
	text = niceTextUsual (obj->curr, timeBuf) ;
      bsDestroy (obj) ;
    }
  return text ;
}

/************************************/

static char* niceTextJava (BS bs, int style)
{ 
  static KEY old =  0, oldTime = 0 ;
  static Stack s = 0 ;
  char timeBuf[25] ;

  if (style < 2) old = oldTime = 0 ;
  if (!bs) return 0 ; /* allow to simply reset old to zero */

  if (!s) s = stackCreate(128) ;
  stackClear (s) ;

  if (style == 0)  /* old style, allways start with a ? */
    pushText (s, "?") ;

  if (isModele)
    {
      catText (s, class(bs->key) == _VModel ? "Model" : "tag") ;
      catText (s, messprintf ("?%s?", freejavaprotect(name(bs->key)))) ;
    }
  else if (bs->key <= _LastC)
    {
      if (!old || bs->key != old)
	catText (s, "txt") ;
      catText (s, messprintf ("?%s?",freejavaprotect(bsText(bs)))) ;
    }
  else if (bs->key == _Int)
    {
      if (!old || bs->key != old)
	catText (s, "int") ;
      catText (s, messprintf ("?%d?", bs->n.i)) ;
    }
  else if (bs->key == _Float)
    {
      if (!old || bs->key != old)
	catText (s, "float") ;
      catText (s, messprintf ("?%1.7g?", bs->n.f)) ;
    }
  else if (bs->key == _DateType)
    {
      if (!old || bs->key != old)
	catText (s, "date") ;
      catText (s, messprintf ("?%s?", timeShowJava (bs->n.time, timeBuf, 25))) ;
    }
  else if (!class(bs->key)) /* a tag */
    {
      if (!old || class(old))
	catText (s, "tag") ;
      catText (s, messprintf ("?%s?", freejavaprotect(name(bs->key)))) ;
    }
  else
    {
      if (!old || class(bs->key) != class(old))
	catText (s, messprintf ("%s",freejavaprotect(className(bs->key)))) ;
      catText (s, messprintf ("?%s?", freejavaprotect(name(bs->key)))) ;
    }
  if (dumpTimeStamps && bs->timeStamp &&
      (!oldTime || bs->timeStamp != oldTime ))
    catText (s, name(bs->timeStamp)) ;

  old = bs->key ; oldTime = bs->timeStamp ;
  return stackText (s,0) ;
}

/*************/
static char *qm = 0 ; /* character or sequence of chars to use as a quote mark */
static char* niceTextPerl(BS bs)
{
  char *title ;
  char timeBuf[25] ;

  if (isModele)
    return name (bs->key);

  if (bs->key <= _LastC)
    return messprintf ("ty=>tx,va=>%s%s%s", qm, bsText(bs), qm);
  else if (bs->key == _Int)
    return messprintf ("ty=>in,va=>%d", bs->n.i);
  else if (bs->key == _Float)
    return messprintf ("ty=>fl,va=>%1.7g", bs->n.f);
  else if (bs->key == _DateType)
    return messprintf ("ty=>dt,va=>%s%s%s", qm, timeShow (bs->n.time, timeBuf, 25), qm);
  else if (iskey(bs->key) == 2)
    {
      if ((title = niceExpandText(bs->key, 0)))
	return 
	  messprintf ("ty=>ob,cl=>%s%s%s,va=>%s%s%s,ti=>%s%s%s", 
		      qm, className(bs->key), qm, qm, name(bs->key), qm, qm, title, qm);
      else 
	return 
	  messprintf ("ty=>ob,cl=>%s%s%s,va=>%s%s%s", 
		      qm, className(bs->key), qm, qm, name(bs->key), qm);
    }
  else if (!class(bs->key)) /* a tag */
    return messprintf ("ty=>tg,va=>%s%s%s", qm, name(bs->key), qm);
  else
    return messprintf ("ty=>ob,cl=>%s%s%s,va=>%s%s%s,mt=>1", 
		       qm, className(bs->key), qm, qm, name(bs->key), qm);
}

/**************/
				/* recursive routine */
				/* RD 980716 simpler layout */
static void niceBSPerl (BS bs, BS bsm, int position, int level)
{
  static int attachRecursionLevel = 0 ;
  int i ;

  while (bs)
    {
      bsModelMatch (bs, &bsm) ;
      for (i=position; i<level; i++)	/* indent */
	freeOut("\t");

      freeOut ("{") ;
      freeOut (niceTextPerl(bs)) ;

      if (bs->right) 
	{ freeOut (",Pn=>[\n") ;
	  niceBSPerl(bs->right, bsModelRight(bsm), 0 /*level */ , level+1) ;
	  freeOut ("]") ;
	}


      else if (!attachRecursionLevel &&
	       bsm && bsm->right && class(bsm->right->key) == _VQuery)
	{ 
	  char *cq, *question = name(bsm->right->key) ;
	  OBJ obj2 = 0 ; BS bs1 = bs ;
	  KEYSET ks1 ; KEY key1 ;
	  int ii ;
	  BOOL wholeObj = !strcmp("WHOLE",question) ;
	  
	  cq = question + strlen(question) - 1 ;
	  while (cq > question && *cq != ';')
	    cq-- ;
	  if (*cq == ';')
	    cq++ ;
	  wholeObj = !strcmp("WHOLE",cq) ;
	  
	  /* position k1 upward on first true non-tag key */
	  while (bs1->up && 
		 ( pickType(bs1->key) != 'B' ||
		  !class(bs1->key) || bsIsComment(bs1)))
	    bs1 = bs1->up ;
	  ks1 = queryKey (bs1->key, !wholeObj ? question : "IS *") ;
	  if (wholeObj || queryFindLocalise (0, 0, cq)) 
	    for (ii = 0 ; ii < keySetMax (ks1) ; ii++)
	      { key1 = keySet (ks1, ii) ;
		obj2 = bsCreate(key1) ;
		if (!obj2)
		  continue ;
		attachRecursionLevel++ ;
		if (wholeObj || queryFindLocalise (obj2, key1, 0)) 
		  {
		    freeOut (",Pn=>[\n") ;
		    niceBSPerl (obj2->curr->right, 
				bsModelRight(obj2->modCurr), 
				0 /*level */, level + 1) ;  
		    freeOut ("]") ;
		  }
		bsDestroy (obj2) ;
		attachRecursionLevel-- ;
	      }
	  keySetDestroy (ks1) ;
	}
      freeOut ("}") ;

      bs = bs->down ;
      if (!bs)
	break ;			/* exit from loop */

      freeOut (",\n") ;
    }
}

/**************/
            /* recursive routine */
static int niceBSJava (BS bs, BS bsm, int position, int level, int style) 
{
  static int attachRecursionLevel = 0 ;
  int i;

  while (bs)
    {
      bsModelMatch (bs, &bsm) ;

				/* tab over to the correct level */
      for (i=position;i<level;i++)
	freeOut("\t");
      freeOut(niceTextJava(bs, style));

      if (bs->right)
	level=niceBSJava (bs->right, bsModelRight(bsm), level,level+1, style);

      else if (!attachRecursionLevel &&
	       bsm && bsm->right && class(bsm->right->key) == _VQuery)
	{ 
	  char *cq, *question = name(bsm->right->key) ;
	  OBJ obj2 = 0 ; BS bs1 = bs ;
	  KEYSET ks1 ; KEY key1 ;
	  int ii ;
	  BOOL wholeObj = !strcmp("WHOLE",question) ;
	  
	  cq = question + strlen(question) - 1 ;
	  while (cq > question && *cq != ';')
	    cq-- ;
	  if (*cq == ';')
	    cq++ ;
	  wholeObj = !strcmp("WHOLE",cq) ;
	  
	  /* position k1 upward on first true non-tag key */
	  while (bs1->up && 
		 ( pickType(bs1->key) != 'B' ||
		  !class(bs1->key) || bsIsComment(bs1)))
	    bs1 = bs1->up ;
	  ks1 = queryKey (bs1->key, !wholeObj ? question : "IS *") ;
	  if (wholeObj || queryFindLocalise (0, 0, cq)) 
	    for (ii = 0 ; ii < keySetMax (ks1) ; ii++)
	      { key1 = keySet (ks1, ii) ;
		obj2 = bsCreate(key1) ;
		if (!obj2)
		  continue ;
		attachRecursionLevel++ ;
		if (wholeObj || queryFindLocalise (obj2, key1, 0)) 
		  level=niceBSJava (obj2->curr->right, 
				    bsModelRight(obj2->modCurr), 
				    level, level + 1, style) ;
		bsDestroy (obj2) ;
		attachRecursionLevel-- ;
	      }
	  keySetDestroy (ks1) ;
	}
    
      bs = bs->down ;
      if (!bs)
	break ;			/* exit from loop */

      freeOut("\n");
      position = 0;
  }

  level--;
  return level;
}

/**************/
            /* recursive routine */
static int niceBSC (BS bs, BS bsm, int position)
{
  int i, level ;
  char direction ;
  char *cp, buf[8] ;

  level = position + 1;
  while (bs)
    {
      bsModelMatch (bs, &bsm) ;
	
      i = level - position ;
      switch (i)
	{
	case 0: direction = '.' ; break ;
	case 1: direction = '>' ; break ;
	case -1: direction = '<' ; break ;
	default:
	  buf[0] = 'l' ;
	  while (i < -255)
	    { buf[1] = 255 ; freeOutBinary (buf,2) ; i += 255 ; }
	  if (i < 0)
	    { buf[1] = -i ; freeOutBinary (buf,2) ; } 
	  direction = '.' ;
	  break ;
	}
      buf[0] = direction ; 
      position = level ;

      /* dump cell itself */
      if (bs->key <= _LastC)
	{
	  buf[1] = 't' ; 
	  freeOutBinary (buf,2) ; 
	  cp = bsText(bs) ;
	  freeOutBinary (cp, strlen(cp) + 1) ;
	}
      else if (bs->key == _Int)
	{
	  int ii = bs->n.i ;
	  
	  buf[1] = 'i' ; 
	  freeOutBinary (buf,2) ; 
	  freeOutBinary ((char*)&ii, 4) ;
	}
      else if (bs->key == _Float)
	{ 
	  float f = bs->n.f ;
	  
	  buf[1] = 'f' ; 
	  freeOutBinary (buf,2) ; 
	  freeOutBinary ((char*) &f, 4) ;
	}
      else if (bs->key == _DateType)
	{ 
	  int ii = bs->n.time ;
	  
	  buf[1] = 'd' ; 
	  freeOutBinary (buf,2) ; 
	  freeOutBinary ((char*) &ii, 4) ;
	}
      else if (iskey(bs->key) == 2)
	{
	  char *title = 0 ;

	  if ((title = niceExpandText(bs->key, 0)))
	    {
	      buf[1] = 't' ;
	      freeOutBinary (buf,2) ; 
	      cp = title ;
	      freeOutBinary (cp, strlen(cp) + 1) ;
	    }
	  else 
	    {
	      buf[1] = 'K' ; buf[2] = class(bs->key) ;
	      freeOutBinary (buf,3) ; 
	      cp = name(bs->key) ;
	      freeOutBinary (cp, strlen(cp) + 1) ;
	    }
	}
      else if (class(bs->key)) /* an emtpty object */
	{ 
	  buf[1] = 'k' ; buf[2] = class(bs->key) ;
	  freeOutBinary (buf,3) ; 
	  cp = name(bs->key) ;
	  freeOutBinary (cp, strlen(cp) + 1) ;
	}
      else if (!class(bs->key)) /* a tag */
	{ 
	  int ii = bs->key ; 
	  
	  buf[1] = 'g' ; 
	  freeOutBinary (buf,2) ; 
	  freeOutBinary ((char*) &ii, 4) ;
	}
      
      /* dump right */
      if (bsm && bs->right)
	{
	  /* this !..@ construction is non compatible with vahan's wacec, now removed */
	if (bsIsType(bsm->right))
          freeOutBinary( "!", 1);      /* entering subtype */

	position = niceBSC (bs->right, bsModelRight(bsm), position) ;

        if (bsIsType(bsm->right))
          freeOutBinary( "@", 1);      /* leaving subtype */

	}

      /* move down, dump in the cuurent while loop */
      bs = bs->down ;
      if (!bs)
	break ;			/* exit from loop */
      
      freeOutBinary("\n", 1);
  }
  /*
    buf[0] = '<' ; 
    freeOutBinary (buf,1) ; 
  */
  return position ;
}

/**************/

static int niceBShuman (BS bs, int x, int y)	/* recursive routine */
{
  int length, xPlus, yMax ;
  char *text, *cp ;
  char buf100[100] ;

  for ( ; bs ; bs = bs->down)
    {
      if (!dumpComments && class(bs->key) == _VComment)
	continue ;

      if (!(text = niceExpandText (bs->key, buf100)))
	text = niceTextUsual (bs, buf100) ;
   
      if (x < LINE_LENGTH - 18)
	length = LINE_LENGTH - x ;
      else
	length = 18 ;

      xPlus = 0 ;
      yMax = y ;
  
      uLinesText (text, length) ;
      if ((cp = uNextLine(text)))
	{ 
	  if (strlen(cp) > xPlus)
	    xPlus = strlen(cp) ;
	  if ( x + xPlus > LINE_LENGTH )
	    { yMax++ ; x = LINE_LENGTH - length ;}
	  freeOutxy (cp,x,yMax++) ;	/* write out this node */
	}
      while ((cp = uNextLine(text))) /* indent following lines */
	{ freeOutxy (cp,x+2,yMax++) ;
	  if (strlen(cp)+2 > xPlus)
	    xPlus = strlen(cp)+2 ;
	}
  
      xPlus += x ;
      xPlus = xPlus + 6 - xPlus%4 ;
  
      if (bs->right)		/* to the right at same y */
	y = niceBShuman (bs->right, xPlus, y) ;
  
      if (yMax > y)
	y = yMax ;
    }
  
  return y ;
}

/***************************************************************************/

void niceDump (OBJ obj, char style)
{ 
  static BOOL isFirst = TRUE ;
  BS oldDown = 0 ;
  char buf100[25] ;

  if (isFirst)
    { qm = getenv("ACEQM"); /* quote mark for perl output */
      if (!qm) 
	{ qm = messalloc (2);
	  qm[0] = 1; qm[1] = 0; /* default character is ^A */
	}
      isFirst = FALSE ;
    }

  isModele = (class(obj->root->key) == _VModel) || isModel(obj->root->key) ;

  if (obj->curr)
    switch (style)
      {
      case 'h': case 'H':
	freeOut (messprintf ("%s %s", className (obj->key),
			     niceTextUsual (obj->root, buf100))) ;
	if (obj->curr != obj->root)
	  {
	    oldDown = obj->curr->down ; obj->curr->down = 0 ;
	    niceBShuman (obj->curr,2,1) ;
	    obj->curr->down = oldDown ;
	  }
	else
	  niceBShuman (obj->curr->right,2,1) ;
	freeOutxy ("",0,1) ; /* space line */
	break ;

      case 'j': case 'J':
	if (obj->curr != obj->root)
	  {
	    freeOut (niceTextJava(obj->root, style == 'J' ? 1 : 0)) ; 
	    oldDown = obj->curr->down ; obj->curr->down = 0 ;
	    niceBSJava (obj->curr, obj->modCurr, 0,1, style == 'J' ? 2 : 0) ;
	    obj->curr->down = oldDown ;
	  }
	else
	  {
	    oldDown = obj->curr->down ; obj->curr->down = 0 ;
	    niceBSJava (obj->curr, obj->modCurr, 0,0, style == 'J' ? 1 : 0) ;
	    obj->curr->down = oldDown ;
	  }
	freeOut("\n") ;
	break;

      case 'p':
	if (!obj->curr->right) break ;
	freeOut ("{");
	freeOut (niceTextPerl (obj->root)) ;
	freeOut (",Pn=>[\n");
	niceBSPerl (obj->curr->right, bsModelRight(obj->modCurr), 0, 1) ;
	freeOut ("]}\n");
	break ;

      case 'C':
	{ 
	  char buf[2] ; 
	  char *cp = name(obj->key) ;

	  buf[0] = 'N' ; buf[1] = class(obj->key) ;
	  freeOutBinary (buf, 2) ;
	  freeOutBinary (cp, strlen (cp) + 1) ;
	  if (obj->curr->right) 
	    niceBSC (obj->curr->right, bsModelRight(obj->modCurr), 0) ;
	  buf[0] = '#' ; 
	  freeOutBinary (buf,1) ; 
	}
	break ;

#ifdef JUNK   /* not yet done */
>       case 'b': /* Boulder */
>       -name=AC4
>         -class=Sequence
>         FEATURES=....
>         HOMOL=...
> 
> 
>       niceText = niceTextBoulder;
>       freeOut(niceText(obj->root));
>       niceBSBoulder (obj->root->key, obj->curr->right, obj->modCurr->right, 0,1);
>       freeOut("\n");
>       break;
> 
#endif

      default:
	messcrash ("Bad style %c in niceDump", style) ;
      }
}

/************************************************************************/
/************************************************************************/
