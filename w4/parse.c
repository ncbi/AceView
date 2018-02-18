/*  File: parse.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 * -------------------------------------------------------------------
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
 * along with this program; if n t, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: ace file parser
 * Exported functions: parseControl, parseFile, parseBuffer.
 * HISTORY:
 * Last edited: Jan 18 09:53 1999 (fw)
 * * Jul 13 03:15 1998 (rd)
 * * May 12 15:19 1993 (mieg): added protected mode
 * * Apr 26 15:31 1993 (cgc): removed -T Top option
 * * Apr 25 23:25 1993 (cgc): added -R replace option: relpos, _bsHere
 * * Jun  5 16:18 1992 (mieg): added NON_GRAPHIC to parse in non-graphic mode
 * * Mar  4 02:48 1992 (rd): added -R option for renaming objects
 * Created: Sun Jan 12 17:52:35 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: parse.c,v 1.25 2014/11/19 22:20:16 mieg Exp $ */

#include "acedb.h"
#include "freeout.h"		/* for freeOutf */

#include "a.h"			/* for arraykill */
#include "bs.h"
#include "lex.h"
#include "systags.h"
#include "sysclass.h"
#include "pick.h"
#include "mytime.h"		/* for timeParse() */
#include "query.h"		/* queryCheckConstraints() */

#ifndef NON_GRAPHIC
#include "call.h"
#include "session.h"
#endif /* !NON_GRAPHIC */

#define FREESPECIAL  freespecial ("\n\"/\\\t@") ;
static Associator pipeAss = 0 ; /* handle pclose */

#ifndef NON_GRAPHIC
#include "display.h"
#include "main.h"

static void readAll (void) ;
static void
  parseDraw(void),
  openFile(void), closeFile(void),
  readItem(void), readLine(void), skipItem(void) ;

static MENUOPT parseOpts[] = {
  { openFile, "Open  ace File"},
#if defined(MACINTOSH)
  { readAll,  "Read all"},			     /* MAC */
#elif defined(WIN32)
  { readAll,  "Read all (Use \"Esc\" to break in)"}, /* WINDOWS */
#else
  { readAll,  "Read all (hit key F4 to break in)"},  /* UNIX-X11 */
#endif /* UNIX */
  { readItem, "Read one Item"},
  { readLine, "Read one Line"},
  { skipItem, "Skip to next Item"},
#ifndef NON_GRAPHIC
  { graphDestroy, "Quit"},
#endif
  {0,0}
} ;

static Graph parseGraph = 0 ;
static int buttonBox = 0 ;
#endif /* !NON_GRAPHIC */

BOOL parseKeepGoing = FALSE ; /* set to TRUE in tacemain */
ParseFunc parseFunc[256] ;  /* == MAXTABLE of lexsubs.c */
extern BOOL  isKernelParsing  ; /* mieg 2006_02_17 kernel communication with objcache.c */

/*********** stuff local to the parser routine and initialiser **********/
/************************************************************************/

/*
* internal state of parseLine() analyzer
*/

enum ParseState 
	{
	DONE, 
		/* XQX
		* Indicates that we saw EOF on the file we are parsing.
		* I'm not sure why it gets a state, since you won't
		* call parseLine any more after this.
		*/
	OUTSIDE, 
		/* XQX
		* OUTSIDE means parseLine is ready to see the start of
		* an object.  It will consume blank lines, comment lines,
		* and skip over '!' objects looking for the start of an
		* object.  It also recognizes and processes -D, -R, and 
		* -A lines.
		*/
	REF_BLOCKED,
		/* XQX
		* After parseLine handles the start of an object, it moves
		* to this state.  This state really just does some simple
		* initialization that could have been done at the end of
		* OUTSIDE.  When the initialization is done, it falls
		* directly into INSIDE to parse data.
		*/
	INSIDE 
		/* XQX
		* This is the state when we are reading an object.  Each
		* line is another tag, possibly with data.  parseLine
		* remains in this state until it sees the end of an
		* object.
		*/
	} ;

char parseLineText[64] ;	/* global in parse.h */

#ifndef NON_GRAPHIC
static char itemText[64], nerrorText[64] ;
static char directory[DIR_BUFFER_SIZE], filename[FIL_BUFFER_SIZE] ;
static int  lineBox, itemBox, nerrorBox ;
#endif /* !NON_GRAPHIC */

static int parseShow = 0 ;

static enum ParseState state = DONE ; /* start off not parsing */
static BOOL isError ;
static int nob = 0;		/* number of objects parsed so far */
static int nerror = 0 ;
static char *obname = 0 ;
static FILE *parseFil = 0 ;
static long filLength = 0 ;
static int level = 0 ;  /* freecard level */
int parseErrorsLastParse;	/* hack to sneak a number out of deep nesting */

/* all these statics are because of callback implementation */

#ifdef PROGBOX			/* Frank's progress bar */
static int progBox ;
#define PROG_X 5
#define PROG_Y 2
#define PROG_LEN 30

/* code to draw - for parseDraw() */
   graphRectangle (PROG_X, PROG_Y, PROG_X + PROG_LEN + 0.5, PROG_Y+1) ;
   progBox = graphBoxStart () ;
   graphFillRectangle (PROG_X, PROG_Y, PROG_X+0.5, PROG_Y+1) ;
   graphBoxEnd () ;

/* code to update - for showState() */
   if (parseFil && progBox && filLength)
     { float position = ftell(parseFil) / (float) filLength ;
       graphBoxShift (progBox, PROG_X + position*PROG_LEN, PROG_Y) ;
     }
#endif /* PROGBOX */

/********************************************/

static void parseArray (KEY key, int classe)
{
  if (classe > 0 && classe < 256 && parseFunc[classe])
    { if (parseFunc[classe] (level,key))
	{ FREESPECIAL ;/* may change freespecial */
	  return ;
	}
    }
  else
#ifndef NON_GRAPHIC
    messout ("No parser available for array class %d", classe) ;
#else
  freeOutf ("// No parser available for array class %d", classe) ;
#endif /* !NON_GRAPHIC */
  if (!parseKeepGoing)
    isError = TRUE ;    /* skip to end of this reference */
  ++nerror ;
  while (freecard(level) && freeword()) ;
  FREESPECIAL ; /* restore */
}

/*******************************************/

#ifndef NON_GRAPHIC
static void showState (void)
{
  if (itemBox)
    { strcpy (itemText, messprintf ("%d  %.40s", nob, obname ? obname : "")) ;
      graphBoxDraw (itemBox, BLACK, WHITE) ;
    }
  if (lineBox)
    { if (state != DONE)
	{ strcpy (parseLineText,
		  messprintf ("%d", freestreamline (level))) ;
	  if (parseFil && filLength)
	    { long pos = ftell (parseFil) ;
	      strcat (parseLineText,
		      messprintf (" (%2.0f%%)", (100.0 * pos) / filLength)) ;
	    }
	}
      else			/* can't know where I am! */
	*parseLineText = 0 ;
      graphBoxDraw (lineBox, BLACK, WHITE) ;
    }
  if (nerrorBox)
    { strcpy (nerrorText, messprintf ("%d", nerror)) ;
      graphBoxDraw (nerrorBox, BLACK, WHITE) ;
    }
}
#endif /* !NON_GRAPHIC */

/********************************************/

BOOL overRideParseProtection = FALSE ; /* used in kernel */

static void parseLine (int isSkip, KEYSET ks, KEYSET ksDirty, BOOL silent) /* deals with one entry */
{
  static int recursion= 0 ; /* around calls that would reenter this code */
  int   table, ix ;
  float fx ;
  mytime_t maDate ;
  char  *word, dummyChar, *card, *obname2 = 0 ;
  KEY   key, tag, type, relpos ;
  static BOOL isDelete, isRename, isAlias, isBang ;
  static OBJ obj = 0 ;
  static KEY  ref ;
  static char *errtext = 0 ;
  BOOL constraintError = FALSE ;
#if defined(NEW_MODELS)
  extern BOOL READING_MODELS ;
#endif

  isError = FALSE ;
  if (isSkip == 2)		/* as in parseDestroy */
    goto done ;
  if (isSkip)
    { if (state == DONE)
        goto done ;
      if (state == OUTSIDE)
	return ;
      goto abort ;
    }

  switch (state)
    {
    case DONE :
      {
      goto done ;
      }

    case OUTSIDE :
      {
      isBang = FALSE ;
      while (TRUE)  /* look for  next object */
	{ if (!(card = freecard(level)))
	    goto done ;
	  if (!isBang)
	    { word = freeword() ;
	      if (word && *word != '!')
		break ;
	      if (word && *word == '!') /* jump paragraph */
		isBang = TRUE ;
	    }
	  else
	    { word = freeword() ;
	      if (!word)
		isBang = FALSE ;  /* end of jumped paragraph */
	    }
	}
#ifndef NON_GRAPHIC
      if (parseGraph && parseShow == 3)
	{ strncpy (parseLineText, card, 60) ;
	  graphBoxDraw (lineBox, BLACK, WHITE) ;
	}
#endif /* !NON_GRAPHIC */

      isDelete = isRename = isAlias = FALSE ;
      if (!strcmp (word, "-D"))	/* symbol for deletion */
	{ isDelete = TRUE ;
	  if (!(word = freeword()))
	    { errtext = "No class after -D" ;
	      goto abort ;
	    }
	}
      else if (!strcmp (word, "-R"))	/* symbol for renaming */
	{ isRename = TRUE ;
	  if (!(word = freeword()))
	    { errtext = "No class after -R" ;
	      goto abort ;
	    }
	}
      else if (!strcmp (word, "-A"))	/* symbol for aliasing */
	{ isAlias = isRename = TRUE ;
	if (1 || !(word = freeword()))    /* mieg jan 2002, aliasing is a nightmare */
	    { errtext = "No class after -A" ;
		  errtext = "-A Aliasing is now forbidden, use -R rename in place" ;
	      goto abort ;
	    }
	}

      /* recognize > as fasta format for DNA */
      if (*word == '>')
	{ 
	  table = pickWord2Class("DNA") ;
	  /*
	  * obname is global - it remembers the object name for when this function
	  * is called for the next line
	  */
	  if (*(word + 1))
	    obname = word + 1 ;
	  else
	    obname = freeword () ;
	  goto fasta ;
	}

      /* XQX
      * word is a class name - get the matching class number
      */
      table = pickWord2Class(word) ;
          /* Could search classes to parse subclasses */

      if (pickList[table].type != 'A'
	  && pickList[table].type != 'B')
	{ errtext = messprintf ("Bad class name : %s", word ? word : "NULL") ;
	  goto abort ;
	}

      if (!table)
	{ errtext = messprintf("Unrecognised Class %s ", word) ;
	  goto abort ;
	}

      if (pickList[table].protected && !overRideParseProtection)
	{ errtext = "Trying to read/alter a protected class" ;
	  goto abort ;
	}

#if !defined(NEW_MODELS)
		/* old models */
      if (pickList[table].type == 'B') /* a kludge to check existence of the model */
	{
	  if ((obj = bsCreate (KEYMAKE (table,0))))
	    bsDestroy (obj) ;
	  else
	    {
	      errtext = messprintf ("Missing model of class %s", word) ;
	      goto abort ;
	    }
	}
#else
	      /* new models */
      if (!READING_MODELS && pickList[table].type == 'B' && !pickList[table].model)
	{
	  errtext = messprintf ("Missing model of class %s", word) ;
	  goto abort ;
	}
#endif

      freestep(':') ;               /* optional colon after class */
      obname = freeword() ; /* fasta may get it glued on > char */
	  /*
	  * obname is global - it remembers the object name for when this function
	  * is called for the next line
	  */
    fasta:
      if (!obname || !strlen(obname2=lexcleanup(obname,0)))
	{ errtext = "while parsing fasta format :  No object name" ; 
	messfree (obname2) ;
	  goto abort ;
	}
      messfree (obname2) ;
      if (isRename)
	{ ref = 0 ;
	  lexword2key (obname, &ref, table) ;
	  ++nob ;
#ifndef NON_GRAPHIC
	  switch (parseShow)
	    {
	    case 1:
	      if (!(nob%100))
		showState () ;
	      if (!(nob%5000))
		graphRedraw () ;
	      break ;
	    case 2: case 3:
	      showState () ;
	      break ;
	    default:
	      break ;
	    }
#endif /* !NON_GRAPHIC */
				/* always check the syntax of the ace file */
	  if (!(obname = freeword()))
	    { errtext = "No new name for rename/alias (-R or -A)" ;
	      goto abort ;
	    }
	  if (!ref)
	    return ; /* if obname unknown, do nothing */
	  { int oldLevel  = level ; /* mieg: feb 93, To avoid a bug,
				       because level gets altered by lexAlias call
				       which itself call the present code.
	                               Notice that if bsFuse were reentrant
				       this would cause the bug again
                                       */
	    recursion++ ;
	    if (recursion >= 2)
	      messcrash ("Multiple recursive call to lexalias inside parser") ;
	    if (!lexAlias (&ref, word, FALSE, isAlias))
	      { errtext = "Rename failed" ;
		recursion-- ;
		goto abort ;
	      }

	    state = OUTSIDE ;
	    level = oldLevel ;
	    recursion-- ;
	  }
	  return ;
	}

      if (isDelete)
	{
	  if (!lexword2key(obname, &ref, table))
	    return ;
	  if (lexDestroyAlias (ref)) /* true if this is an alias */
	    return ;
	  switch (pickType (ref))
	    {
	    case 'A':
	      arrayKill (ref) ;
	      break ;
	    case 'B':
	      if ((obj = bsUpdate(ref)))
		bsKill(obj) ;
	      break ;
	    default:
	      errtext = "Can't handle because unknown class type" ;
	      goto abort ;
	    }
	  state = OUTSIDE ;
	  return ;
	}

      /* XQX
      * we know the class number and the object name.
      * find or create the object
      */
      lexaddkey (obname, &ref, table) ;
      ++nob ;

#ifndef NON_GRAPHIC
	  switch (parseShow)
	    {
	    case 1:
	      if (!(nob%100))
		showState () ;
	      if (!(nob%5000))
		graphRedraw () ;
	      break ;
	    case 2: case 3:
	      showState () ;
	      break ;
	    default:
	      break ;
	    }
#endif /* !NON_GRAPHIC */

      switch (pickType (ref))
	{
	case 'A':
	  /* XQX
          * if we have a keyset, add the newly created object to
	  * the keyset
	  */
	  if (ks)
	    keySet (ks, keySetMax (ks)) = ref ;

	  /* XQX
	  * this is an A object, so it has no tags.  It does have
	  * a special parser, so go call it.
	  */
	  parseArray (ref,table) ;
	  return ;

	case 'B':
	  break ;		/* handled by REF_BLOCKED */

	default:
	  errtext = "Can't handle because unknown class type" ;
	  goto abort ;
	}
      state = REF_BLOCKED ;
      } /* end case OUTSIDE */

    case REF_BLOCKED :			/* note fall through */
      {
      if (!(obj = bsUpdate (ref)))
	{ errtext = messprintf ("// Object %s is locked - free before continuing",
		   obname) ;
	  goto abort ;
	}
      state = INSIDE ;
      return ;
      } /* end case REF_BLOCKED */

    case INSIDE :
      {
      if (!(card = freecard(level)))
	goto done ;
#ifndef NON_GRAPHIC
      if (parseGraph && parseShow == 3)
	{ strncpy (parseLineText, card, 60) ;
	  graphBoxDraw (lineBox, BLACK, WHITE) ;
	}
#endif /* !NON_GRAPHIC */
      if (!(word = freeword ()))  /* white line separates objects */
	{
	  /* XQX
          * this null word means we hit the blank line at the end of the
	  * object.  close out the object and save it.  we are done
	  */
	  if (pickList[class(ref)].constraints &&
	      !queryCheckConstraints (obj, ref,pickList[class(ref)].constraints))
	    { errtext = messprintf ("Class:Object %s:%s does not satisfy constraints", name (ref)) ;
	    constraintError = TRUE ; /* prevent jumping next object */
	      goto abort ;
	    }

	  if (ks) keySet (ks, keySetMax (ks)) = ref ;
	  if (ksDirty && bsIsDirty (obj))
	    keySet (ksDirty, keySetMax (ksDirty)) = ref ;
	  bsSave (obj) ;
	  state = OUTSIDE ;
	  /* XQX
          * end of object parsing
          */
	  return ;
	}

      if (*word == '!')
        {
        /* XQX
        * this ! was in front of a tag, so we just ignore it.
        */
	return ;
        }

      if (!strcmp (word, "-D"))	/* symbol for deletion */
	{ 
	  /*
          * -D means we are going to delete a tag or a value after a tag.
          */
          isDelete = TRUE ;
	  if (!(word = freeword()))
	    { errtext = "No tag after -D" ;
	      goto abort ;
	    }
	}
      else
	isDelete = FALSE ;

      if (!lexword2key (word,&tag,0)) /* find tag */
	{ errtext = messprintf ("Unknown tag \"%s\"",word) ;
	  goto abort ;
	}
      freenext() ; freestep (':') ; /* optional colon after tag */

      /* XQX
      * set the object cursor back to the root of the object
      */
      bsGoto(obj, 0) ;		/* pop out of subobjs */

      /* XQX
      * Move the object cursor to the tag;  returns false if you can't
      * have that tag in that object.  This implicitly creates the
      * tag and all the tags that come before it in the hierarchy.
      */
      if (bsAddTag (obj,tag))	/* subobjects come on the same line */
	while ((word = freeword ()))
	  {
          /* XQX
	  * We have read the tag name and moved the object cursor to
	  * that location in the object.  Now, each iteration through 
	  * this loop picks up one more column of data and adds it to 
	  * the object.
	  *
	  * As far as I can tell, bsAddKey() and bsAddData() must advance
	  * the object cursor after depositing the data.
          */
	  scanSubObj:
	    relpos = _bsRight ;
	    if (*word == '-' && strlen(word) == 2)
	      switch (word[1])
		{
		case 'C':	/* standard comment */
		  if (!(word = freeword()))
		    { errtext = "No Comment after -C" ;
		      goto abort ;
		    }
		  bsAddComment (obj, word, 'C') ;
		  if (!(word = freeword()))
		    continue ;
		  if (strcmp(word,"-O")) /* only a timestamp possible */
		    { errtext = "Extra item on -C line" ;
		      goto abort ;
		    }
				/* note fall-through if -O after -C */
		case 'O':	/* timestamp (UserSession class) */
		  if (!(word = freeword()))
		    { errtext = "No UserSession key after -O" ;
		      goto abort ;
		    }
		  lexaddkey (word, &key, _VUserSession) ;
		  bsSetTimeStamp (obj, key) ;
		  continue ;

		case 'R':	/* Replace */
		  if (!(word = freeword()))
		    continue ;
		  relpos = _bsHere ;
		}

	    /*
	    * fetch the type of the tree node that we are about to modify.
 	    */
	    if (!(type = bsType (obj, relpos)))
	      { errtext = messprintf ("class:obj %s:\"%s\" tag %s : %s - off end of model\n",
				      className(bsKey(obj)), bsName(obj), name (tag), word) ;
		goto abort ;
	      }

	    if (class(type) == _VQuery)
	      break ;

	    if
#ifdef NEW_MODELS
		/* new models */
               (class(type) == _VModel) /* a type */
#else
		/* old models */
               (class(type) && type == KEYMAKE(class(type), 1)) /* a type */
#endif /* !NEW_MODELS */
	      {
		if (!bsPushObj(obj))
		  {
                    /*
                    * we have more columns of data than the model allows 
		    */
		    errtext = messprintf ("tag %s : %s - more data than model allows, "
					  "addObj failed\n", name(tag), word) ;
		    goto abort ;
		  }
		goto scanSubObj ;
              }

	      if
#ifdef NEW_MODELS
                (class(type) && (table = type, pickIsA(&table, 0), table))
#else
		((table = class(type)))
#endif /* !NEW_MODELS */
		{
		  char *word2 = 0 ;
 		  /* XQX
		  * The data type is a class number.  That means the data 
		  * element expected is an object name.
		  */
		  if (!strlen(word2=lexcleanup(word,0))) /* refuse blank strings: "    " */
		    { errtext = messprintf ("Blank object names names are illicit, sorry") ;
		    messfree (word2) ;
		    goto abort ;
		    }

		  /* XQX
                  * ???
		  */
		  if (pickKnown (table) &&
		      !lexword2key (word2,&key,table))
		    { errtext = messprintf ("tag %s : class:%s %s - unkown.  You can not create an "
					    "entry in a \"Known\" class in an indirect way (see "
					    "options.wrm)\n",
					    name (tag),pickClass2Word(table), word2) ;
		    messfree (word2) ;
		    goto abort ;
		    }
                  
		  /* XQX
		  * create a key for the referenced object and store it 
		  * in the current object at the cursor
		  */
		  lexaddkey (word,&key,table) ;
		  messfree (word2) ;
		  bsAddKey (obj,relpos,key) ;
		}
	      else if (type == _Int)
		{
 		  /* XQX
		  * The model expects an integer here.
		  */
		  if (sscanf (word,"%d%c",&ix,&dummyChar) == 1)
		    bsAddData (obj,relpos,_Int,&ix) ;
		  else
		  {
		    errtext = messprintf ("%s not an integer after %s\n", word, name(tag)) ;
		    goto abort ;
		  }
		}
	      else if (type == _Float)
		{
 		  /* XQX
		  * The model expects a floating point number.
		  * Notice that we scan as a double, but the database only
		  * stores float.  That must be what this range check is for.
		  */
		  double bigfloat = 0 ;
		  if (sscanf (word,"%lf%c",&bigfloat,&dummyChar) == 1 &&
		      bigfloat < 1e40 && bigfloat > -1e40 )
		    {
		      char buf[128] ;

		      if  (bigfloat <= ACE_FLT_MIN  && bigfloat >= -ACE_FLT_MIN)
			bigfloat = 0 ;
		      sprintf (buf, "%1.7g", bigfloat) ; sscanf (buf, "%f", &fx) ;
		      bsAddData (obj,relpos,_Float,&fx) ;
		    }
		  else
		    {
		      errtext = messprintf ("%s not a float after %s\n",
					    word,name(tag)) ;
		      goto abort ;
		    }
		}
	      else if (type == _DateType)
		{
 		  /* XQX
		  * the model expects a date/time here.  see wdoc/time-format
		  * for syntax of the input
		  */
		  if ((maDate = timeParse (word)))
		    bsAddData (obj,relpos,_DateType,&maDate) ;
		  else
		  {
		    errtext = messprintf ("%s not a DateType after %s\n",
					  word,name(tag)) ;
		    goto abort ;
		  }
		}
	      else if (type < _LastC)
		{
 		  /* XQX
		  * ???
		  */
		  bsAddData (obj,relpos,type,word) ;
		}
	      else          /* a tag */
		{
		  /* XQX
		  * The model expects a tag here.  
		  * because the  .ace file can say "tag1 tag2" 
		  * it is sufficent to say tag2 directly
		  * and tag1 is created implicitely to conform to the schema
		  */
		  if (!lexword2key (word,&tag,0)) /* find tag */
		    {
		      errtext = messprintf ("Unknown tag \"%s\" in class ",word, pickClass2Word(table)) ;
		      goto abort ;
		    }
		  if (!bsAddTag(obj, tag))
		    {
		      errtext = messprintf ("Cannot add Tag \"%s\" in class:obj %s:%s at this position"
					    ,word, className(bsKey(obj)), bsName(obj)) ;
		      goto abort ;
		    }
		}
	  }
      else
	{
	  errtext = messprintf ("Tag %s does not fit model of class:obj %s:%s \n"
				, name(tag), className(bsKey(obj)), bsName(obj)) ;
	  goto abort ;
	}

      if (isDelete)
	bsPrune (obj) ;

      /* XQX
      * move cursor back to root of object so it isn't pointing some
      * random place if we accidentally use it somewhere
      */
      bsGoto(obj, 0) ;

      /*
      * done parsing a line of tag data
      */
      return ;
      } /* end case INSIDE */
    }

abort:
  if (errtext)
    {
      if (!parseKeepGoing)
	messout ("Parse error near line %d in %s:%.25s : %s\n",
		 freestreamline(level) - 1, obj ?  className(ref) : "" ,  obj ?  name(ref) : "No obj" , errtext) ;
      else
	freeOutf ("// Parse error near line %d in %s:%.25s : %s\n",
		  freestreamline(level) - 1, obj ?  className(ref) : "" , obj ?  name(ref) : "No obj" , errtext) ;
      
      messdump ("Parse error near line %d in %s:%.25s : %s\n",
		freestreamline(level) - 1, obj ?  className(ref) : "" ,  obj ?  name(ref) : "No obj" , errtext) ;
    }
  if (!constraintError)
    while (freecard(level) && freeword()) ; /* skip to end of this reference */
  bsDestroy (obj) ;
  state = OUTSIDE ;
  if (!parseKeepGoing)
    isError = TRUE ;
  ++nerror ;
#ifndef NON_GRAPHIC
  showState () ;
  return ;
#else   /* In non graphic mode, fall through to DONE on errors */
  if (!isError)
    return ;
#endif


done :
  if (ks) keySet (ks, keySetMax (ks)) = ref ;
  if (ksDirty && bsIsDirty (obj))
    keySet (ksDirty, keySetMax (ksDirty)) = ref ;
  bsSave (obj) ;
  freeclose (level) ; /* skip to end of this stream */
  parseFil = 0 ;
  state = DONE ;		/* must come before showState */
#ifndef NON_GRAPHIC		/* GRAPHIC */
  if (!recursion && /* the recursion is private to the code, not for the end user */
      isSkip < 2)  /* in case parseDestroy, the graph is gone, thus show state fails */
    { showState () ;
      if (parseGraph)
	{ parseOpts->text = "Open  File" ;
	  parseOpts->f = openFile ;
	  if (buttonBox < 0)
	    { buttonBox = - buttonBox ;
	      graphBoxClear (buttonBox + 2) ;
	    }
	  parseDraw () ;
	}
    }
#else  /* NON_GRAPHIC */
  if (!recursion && !silent)
    freeOutf("// %d objects read with %d errors\n", nob, nerror) ;
#endif /* NON_GRAPHIC */
  parseErrorsLastParse = nerror;
  return ;
}

/****************************************************************/

#ifndef NON_GRAPHIC

static void readAll (void)
{
  /* NII to NWG assess the efficiency of the hash system */
  /* extern int NII, NIC, NWC, NWF, NWG  ; */

  parseShow = 1 ;
  do { parseLine (FALSE, 0, 0, FALSE) ;
     } while (
	      !messIsInterruptCalled() &&
	      !isError && state != REF_BLOCKED && state != DONE) ;
  /*  printf ("\nNII = %d, NIC =%d, NWF = %d NWG = %d NWC = %d\n\n",
      NII, NIC, NWF, NWG, NWC) ;  */
}

static void readItem (void)
{
  parseShow = 2 ;
  do { parseLine (FALSE, 0, 0, FALSE) ;
     } while (state == INSIDE) ;
  parseLine (FALSE, 0, 0, FALSE) ;		/* get next item's name */
}

static void readLine (void)
{
  parseShow = 3 ;
  parseLine (FALSE, 0, 0, FALSE) ;
}

static void skipItem (void)
{
  parseShow = 3 ;
  parseLine (TRUE, 0, 0, FALSE) ;
  parseLine (FALSE, 0, 0, FALSE) ;
}

static void openFile (void)
{
  char *end, cprs = 0 ;
  int n ;
  BOOL isPipe = FALSE ;

  /* write access may be lost doing a save */

  if (!sessionGainWriteAccess())
    {
      /* may occur is somebody else grabbed it */
      messout ("Sorry, you no longer have write access");
      graphDestroy () ;
      return ;
    }

  if (!*directory)
    {
      if (getenv ("ACEDB_DATA"))
	strcpy (directory, getenv("ACEDB_DATA")) ;
      else if (filName ("rawdata", 0, "r"))
	strcpy (directory, filName ("rawdata", 0, "r")) ;
    }

  end = "ace*" ;
  parseFil = filqueryopen (directory, filename, end, "r",
		      "Which ace file do you wish to parse ?") ;
  if (!parseFil)
    return ;
  if (!graphExists (parseGraph) || /* maybe the window is gone */
      !*filename)
    {
      if (pipeAss && assFind (pipeAss, parseFil, 0))
	{ pclose (parseFil) ; assRemove (pipeAss, parseFil) ; }
      else
	filclose(parseFil) ;
      parseFil = 0 ;
      return ;
    }

  /* autorecognition */
  n = strlen(filename) ;
  if (n > 3 && !strcmp(filename+n-3, ".gz")) cprs = 'g' ;
  else if (n > 2 && !strcmp(filename+n-2, ".Z")) cprs = 'Z' ;
  else if (n > 4 && !strcmp(filename+n-4, ".zip")) cprs = 'P' ;
  switch (cprs)
    {
#if defined(WIN32)
    case 'P' :
      if (pipeAss && assFind (pipeAss, parseFil, 0))
	{ pclose (parseFil) ; assRemove (pipeAss, parseFil) ; }
      else
	filclose(parseFil) ;
      parseFil = 0 ;
      parseFil = callScriptPipe ("pkunzip",messprintf("-c %s/%s",
						 directory, filename)) ;
      if (!pipeAss) pipeAss = assCreate () ; assInsert (pipeAss, parseFil, 0) ;
      isPipe = TRUE ;
      break ;
#else  /* UNIX */
    case 'Z' :
      if (pipeAss && assFind (pipeAss, parseFil, 0))
	{ pclose (parseFil) ; assRemove (pipeAss, parseFil) ; }
      else
	filclose(parseFil) ;
      parseFil = 0 ;
      parseFil = callScriptPipe ("zcat",messprintf(" %s/%s",
					      directory, filename)) ;
      if (!pipeAss) pipeAss = assCreate () ; assInsert (pipeAss, parseFil, 0) ;
      isPipe = TRUE ;
      break ;
    case 'g' :
      if (pipeAss && assFind (pipeAss, parseFil, 0))
	{ pclose (parseFil) ; assRemove (pipeAss, parseFil) ; }
      else
	filclose(parseFil) ;
      parseFil = 0 ;
      parseFil = callScriptPipe ("gunzip",messprintf(" -c %s/%s",
					      directory, filename)) ;
      if (!parseFil)
	parseFil = callScriptPipe ("gzip", messprintf(" -dc %s/%s",
						   directory, filename)) ;
      if (!pipeAss) pipeAss = assCreate () ; assInsert (pipeAss, parseFil, 0) ;
      isPipe = TRUE ;
      break ;
#endif /* UNIX */
    default: break ;
    }
  if (!parseFil)
    return ;

  /*#ifndef NON_GRAPHIC*/
  /* fseek does not work on stdin or pipes */
  if (!isPipe && !fseek (parseFil, 0L, 2))
    { filLength = ftell (parseFil) ;
      if (fseek (parseFil, 0L, 0))
	messcrash ("Can't seek back to original position in openFile") ;
    }
  else
    /*#endif */
    filLength = 0 ;

  if (isPipe)
    level = freesetpipe (parseFil, "") ;
  else
    level = freesetnamedfile (filename, parseFil, "") ;
  FREESPECIAL ;

  parseOpts->f = closeFile ;
  parseOpts->text = "Close File" ;
  if (buttonBox > 0)
    { buttonBox = - buttonBox ;
      graphBoxClear (buttonBox + 2) ;
    }

  messdump ("// Parsing file %s" SUBDIR_DELIMITER_STR "%s\n",
	    directory, filename) ;
  freeOutf ("// Parsing file %s" SUBDIR_DELIMITER_STR "%s\n",
	    directory, filename) ;
  state = OUTSIDE ;
  nob = 0 ;
  nerror = 0 ;
  parseDraw () ;
}

static void closeFile (void)
{
  state = DONE ;
  freeclose (level) ;
  parseFil = 0 ;
/*  *filename = 0 ;	want to leave name of last file read */
#ifndef NON_GRAPHIC
  if (parseGraph)
    { parseOpts->text = "Open  File" ;
      parseOpts->f = openFile ;
      if (buttonBox < 0)
	{ buttonBox = - buttonBox ;
	  graphBoxClear (buttonBox + 2) ;
	}
      parseDraw () ;
    }
#endif
}

/********************************************/

static void parseDraw (void)
{ int line = 9 ;
  if (!graphActivate (parseGraph))
    return ;
  graphPop () ;
  graphClear () ;

  buttonBox = graphButtons (parseOpts, 1, 3, 55) ;

#if !defined (THINK_C)
  if (*filename)
    { graphText ("Filename:",1,1) ;
      graphText (filename,11,1) ;
    }
#endif

  graphText ("Item:",1,line) ;
  itemBox = graphBoxStart () ;
  graphTextPtr (itemText,7,line,64) ;
  graphBoxEnd () ;

  line += 1.5 ;
  graphText ("Line:",1,line) ;
  lineBox = graphBoxStart () ;
  graphTextPtr (parseLineText,7,line,64) ;
  graphBoxEnd () ;

  line += 1.5 ;
  graphText ("Errors:",1,line) ;
  nerrorBox = graphBoxStart () ;
  graphTextPtr (nerrorText,9,line,64) ;
  graphBoxEnd () ;

  showState () ;

  graphRedraw () ;
}

/***************************************/

static void parseDestroy (void)
{
  parseGraph = 0 ;
  if (parseFil)	/* close it correctly */
    { state = DONE ;
      parseLine (2, 0, 0, FALSE) ;
      parseOpts->text = "Open  File" ;
      parseOpts->f = openFile ;
    }
}

/********************************************/

void parseControl (void)
{
  if (!checkWriteAccess())
    return ;

  if (graphActivate (parseGraph))
    { graphPop () ;
      return ;
    }

  if (!*directory)
    {
      if (getenv ("ACEDB_DATA"))
	strcpy (directory, getenv("ACEDB_DATA")) ;
      else if (filName ("rawdata", 0, "r"))
	strcpy (directory, filName ("rawdata", 0, "r")) ;
    }

  parseGraph = displayCreate (DtAce_Parser) ;
  displayPreserve () ; /* mieg 2002, avoid usage by other tools */
  graphTextBounds (60, 10) ;	/* needed to for text box sizing */
  graphRegister (DESTROY, parseDestroy) ;


  *itemText = 0 ;  /* harmless, but i am quite sure i needed that in a peculiar crash */
  *parseLineText = 0 ;
  *nerrorText = 0 ;
  nob = 0 ; nerror = 0 ; state = DONE ;
  parseDraw () ;
}

#else
void parseControl (void) { }
#endif  /* ndef NON_GRAPHIC */

/* RD - make parseFile() and parseBuffer() reentrant */

static BOOL doParseFile (char *fileName, FILE *fileHandle, int newLineBox, KEYSET ks, KEYSET ksDirty, BOOL isPipe, BOOL silent)
{
  int line, oldNerror = nerror, oldNob = nob, oldParseShow = parseShow ;
  long oldLen = filLength ;
  BOOL oldIsKernelParsing ;
  FILE *oldFil = parseFil ;
#ifndef NON_GRAPHIC
  long curr = 0 ;
  int oldItemBox = itemBox, oldLineBox = lineBox, oldNerrorBox = nerrorBox ;


  itemBox = nerrorBox = 0 ;
  lineBox = newLineBox ;
#endif
  nerror = nob = 0 ;

  parseFil = fileHandle ;
  cacheSaveAll () ;    /* any object touched before is flushed */
  oldIsKernelParsing = isKernelParsing ; /* reentrant */
  isKernelParsing = TRUE ; /* if an error occurs, the file should be re-parsed */
#if defined(MACINTOSH)
  setvbuf ( parseFil, nil, _IOFBF, 65536 ) ;
#endif

#ifndef NON_GRAPHIC
  /* fseek does not work on stdin */
  if (!isPipe)
    { curr = ftell (parseFil) ;
      if (!fseek (parseFil, 0L, 2))
	{ filLength = ftell (parseFil) ;
	  if (fseek (parseFil, curr, 0))
	    messcrash ("Can't seek back to original position in parseFile") ;
	}
      else
	filLength = 0 ;
    }
  else
#endif
    filLength = 0 ;

  if (isPipe)
    { if (!pipeAss) pipeAss = assCreate () ;
      assInsert (pipeAss, parseFil, 0) ;
      level = freesetpipe (parseFil, "") ;
    }
  else
    { if (pipeAss && assFind (pipeAss, parseFil, 0))
	assRemove (pipeAss, parseFil) ;
      level = freesetnamedfile (fileName, parseFil, "") ;
    }
  FREESPECIAL ;

  state = OUTSIDE ;
  parseShow = 1 ;
  line = 0 ;
  do { 
       parseLine (FALSE, ks, ksDirty, silent) ;
       if (++line > 100000)
	 {
	   /* mieg 2006_02_17
	    * when parsing very large files it is better if there
	    * are not too many dirty cache-objects so their memory
	    * can be shared 
	    */
	   line = 0 ;
	   cacheSaveAll () ;    /* any object touched before is flushed */
	 }
     } while (
	      !messIsInterruptCalled() &&
	      !isError && state != REF_BLOCKED && state != DONE) ;

  isKernelParsing = oldIsKernelParsing ; /* if an error file should be re-parsed */
  nerror = oldNerror ; nob = oldNob ; parseShow = oldParseShow ;
#ifndef NON_GRAPHIC
  itemBox = oldItemBox ; lineBox = oldLineBox ; nerrorBox = oldNerrorBox ;
#endif
  parseFil = oldFil ; filLength = oldLen ;
  freeclose (level) ;		/* to make sure it is closed */
  if (ks)
    { keySetSort(ks) ;
      keySetCompress (ks) ;
#ifndef NON_GRAPHIC
      keySetShow (ks, 0) ; /* michel, mars 2000, mieg 2003 */
#endif
    }
  return (state == DONE && nerror == 0) ? TRUE : FALSE ;
}

BOOL parseFile (FILE *fileHandle, int newLineBox, KEYSET ks)
{
  KEYSET ksDirty = 0 ;
  return doParseFile (0, fileHandle,newLineBox, ks, ksDirty, FALSE, FALSE) ;
}

BOOL parseNamedFile (char *fileName, FILE *fileHandle, int newLineBox, KEYSET ks)
{
  KEYSET ksDirty = 0 ;
  return doParseFile (fileName, fileHandle,newLineBox, ks, ksDirty, FALSE, FALSE) ;
}

BOOL parsePipe (FILE *fileHandle, int newLineBox, KEYSET ks)
{
  KEYSET ksDirty = 0 ;
  return doParseFile (0, fileHandle,newLineBox, ks, ksDirty, TRUE, FALSE) ;
}

/* mieg, jul 98, the extra param silent is to avoid messout's in parseBuffer */
static BOOL doParseLevel (int mylevel, KEYSET ks,  KEYSET ksDirty, BOOL silent)
{
  int oldNerror = nerror, oldNob = nob, oldParseShow = parseShow ;
  FILE *oldFil = parseFil ;
  int oldlevel = level ;
#ifndef NON_GRAPHIC
  int oldItemBox = itemBox, oldLineBox = lineBox, oldNerrorBox = nerrorBox ;

  lineBox = itemBox = nerrorBox = 0 ;
#endif
  level = mylevel ;
  nerror = nob = 0 ; parseFil = 0 ; parseShow = 0 ;

  state = OUTSIDE ;
  parseShow = 0 ;		/* no status drawing */

  do { parseLine (FALSE, ks, ksDirty, TRUE) ;  /* silent */
     } while (!isError && state != REF_BLOCKED && state != DONE) ;
  if (isError)
    { messout ("// Some sort of parse error") ;
      goto abort ;
    }
  if (state == REF_BLOCKED)
    { messout ("Parsing blocked by locked object") ;
      goto abort ;
    }
  nerror = oldNerror ; nob = oldNob ; parseShow = oldParseShow ; parseFil = oldFil ;
#ifndef NON_GRAPHIC
  itemBox = oldItemBox ; lineBox = oldLineBox ; nerrorBox = oldNerrorBox ;
#endif
  freeclose (level) ;		/* to make sure it is closed */
  if (ks)
    { keySetSort(ks) ;
      keySetCompress (ks) ;
    }
  level = oldlevel ;
  return TRUE ;

 abort:
  state = DONE ;
  parseLine (FALSE, ks, ksDirty, TRUE) ;
  nerror = oldNerror ; nob = oldNob ; parseShow = oldParseShow ; parseFil = oldFil ;
#ifndef NON_GRAPHIC
  itemBox = oldItemBox ; lineBox = oldLineBox ; nerrorBox = oldNerrorBox ;
#endif
  freeclose (level) ;		/* to make sure it is closed */
  if (ks)
    { keySetSort(ks) ;
      keySetCompress (ks) ;
    }
  level = oldlevel ;
  return FALSE ;
}

BOOL parseLevel (int mylevel, KEYSET ks)
{
  KEYSET ksDirty = 0 ;
  return doParseLevel (mylevel, ks, ksDirty, FALSE) ;
}

BOOL parseBuffer (char *text, KEYSET ks)
{
  BOOL t;
  KEYSET ksDirty = 0 ;
  if (!text || !*text)
    return FALSE ;
  level = freesettext (text,"") ;
  FREESPECIAL ;
  t = doParseLevel (level, ks, ksDirty, TRUE) ;
  return t;
}

/* convenience call private to AceC */
BOOL parseAceCBuffer (const char *text, Stack s, void *ace_out, KEYSET ks)
{
  BOOL t;
  int level_out;
  KEYSET ksDirty = 0 ;

  if (!text || !*text)
    return FALSE ;

  level = freesettext (text,"") ; 
  FREESPECIAL ;
  level_out = freeOutSetStack(s) ;
  t = doParseLevel (level, ks, ksDirty, TRUE) ; 
  freeOutClose(level_out);

  return t;
}

void parseOneFile (FILE *fil, KEYSET ks)
{
#ifndef NON_GRAPHIC
  BOOL ParseGraphWasOpen = FALSE ;

  if (!checkWriteAccess())
    return ;

  if (!graphActivate (parseGraph))
    {
      parseGraph = displayCreate (DtAce_Parser) ;
      graphTextBounds (60, 10);	/* needed to for text box sizing */
      graphRegister (DESTROY, parseDestroy) ;


      *itemText = 0 ;  /* harmless, but i am quite sure i needed that in a peculiar crash */
      *parseLineText = 0 ;
      *nerrorText = 0 ;
      nob = 0 ; nerror = 0 ; state = DONE ;

      parseDraw () ;
    }
  else
    {
      ParseGraphWasOpen = TRUE ;
      graphPop () ;
    }

  parseFile(fil,lineBox,ks) ;

  /* For cosmetic reasons: assume the user
     does not need/want to keep it open? */
  if (graphActivate (parseGraph) &&
      !ParseGraphWasOpen )
    graphDestroy() ;
#endif /* !NON_GRAPHIC */
}



/***** end of file ****/
