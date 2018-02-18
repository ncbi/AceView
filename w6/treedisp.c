/*  File: treedisp.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: display/editor for generic tree objects
 * Exported functions:
 	treeDisplay
 * HISTORY:
 * Last edited: May 17 14:12 1999 (edgrif)
 * * May 17 14:12 1999 (edgrif): Set block mode for choose tag window.
 * * Sep 28 14:24 1998 (fw): make timestamp boxes inactive, this 
                    prevents seg-faults and invalid boxNum errors
 * -	extra 'bs->up' tests in lookKbd() to prevent NULL pointer references?
 * * Aug 12 16:56 1996 (srk)
 * * Jun  3 12:42 1996 (rd)
 * * May  8 20:22 1993 (mieg): ON/OFF tag expansion
 * * May  8 20:21 1993 (mieg): updtateMode, to prevent expandText or ON/OFF
 * * Sep 20 01:05 1992 (rd): fixPath() handles new UNIQUE correctly
 * * Aug 16 10:59 1992 (mieg): tagChooser
 * * May  3 22:29 1992 (rd): redid MENU_PRESERVE (nl1 version lost, sorry)
 * * Feb 11 13:09 1992 (mieg): Introduced graphCompletionEntry
 * Created: Wed Oct 16 17:41:46 1991 (rd)
 *-------------------------------------------------------------------
 *
 * $Id: treedisp.c,v 1.33 2016/03/29 22:37:14 mieg Exp $
 */


#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct bblock *BBLOCK ;
typedef BBLOCK BP ;
#define DEFINE_OBJ
typedef struct sobj *OBJ ;
#include "acedb.h"

#include "bs_.h"
#include "cache.h"
#include "dump.h"
#include "biblio.h"
#include "pref.h"
#include "bitset.h"		/* to remember open tags */
#include "fingerp.h"		/* for fpDisplay() fw-980928 */
#include "parse.h"		/* for parseBuffer() fw-980928 */
#include "help.h"		/* for helpOn() fw-981103 */
#include "main.h"		/* for ksetClassComplete */
#include "spreaddisp.h"		/* for spreadTableDisplay() */
#include "tree.h"

/************************************************************/

typedef struct LOOKSTUFF
  { int   magic ;        /* == MAGIC */
    KEY   key ;
    OBJ   objet ;
    Graph graph ;
    int   activebox ;
    Array content ;       /* content of each box: BS */
    Array attachedObjs ;
    Associator bs2box ;
    BS    bsWait ;	/* for use during update */
    AC_HANDLE handleWait ;
    Stack textStack ;	/* for typing in things */
    int   classe ;         /* for completion mechanism */
    BOOL  updateMode, tagChooserMode, justKnownTags ;
    Array longTextArray ;
  } *LOOK ;

static Array bbSetArray = 0 ;
static BOOL autoPreserve = TRUE ;
static BOOL treeJustKnownTags = TRUE ;  /* global toggle, controls tag chooser */
static int MAGIC = 27182 ;

static FREEOPT treeMenu[] = {
  { 8,"Tree display" },
  { 99,"Quit" },
  { 98,"Help"} ,
  { 13,"Print"} ,
  {  0,""},
  /* 
     { 20,"Preserve"} ,           useless, just touch a node
     { 21,"Set AutoPreserve"} ,   useless
     { 22,"Unset AutoPreserve"} , useless
     { 23,"Biblio"} ,         see keySet menu
     { 24,"Neighbours"} ,     see keySet query
     { 14,"Export data"} ,    see keySet menu
     { 25,"Contract"} ,       useless, just touch a tag
     { 26,"Expand"},          useless, just touch a tag
     */
  {  9,"Update"} ,
  { 17,"Toggle Timestamps"} ,
  { 18,"Toggle TagCount"} ,
  {  5,"Show as text"} ,
} ;

static FREEOPT updateMenu[] = {
  {  9,"Tree Edit"} ,
  { 99,"Quit"} ,
  { 98,"Help"} ,
  {  3,"Save"} ,
  {  7,"Add Comment"} ,
  { 10,"Cancel"} ,
  {  8,"Delete"} ,
  { 11,"Edit"} ,
  { 25,"Contract"} ,
  { 26,"Expand"} 
} ;


extern void  lookBiblio	       (OBJ objet, KEY k) ;

static void  lookDestroy       (void) ;
static void  lookRedraw        (void) ;
static void  lookMenu 	       (KEY k) ;
static void  lookPick          (int k) ;
static void  lookKbd           (int k) ;
static void  lookDraw          (LOOK look);
static void  updatePick        (int k) ;
static void  updateKbd	       (int k) ;
static void  defuse            (BS bs) ;
static void  addComment        (char *string) ;
static void  addData           (char *string) ;
static void  addKeyText        (char *string) ;
static void  addKey            (KEY key) ;
static void  editEntry	       (char *text) ;
static void  treeDump          (LOOK look) ;
static int   countSiblings(BS bs) ;
static void  detach1 (LOOK look) ;
static void  treeDispMailer (KEY key) ;	/* at end of file */
static char* drawBStext (BS bs) ;

static BOOL  showTimeStamps = FALSE ;
static BOOL  showTagCount = FALSE ;

#define LINE_LENGTH 120

#define LOOKGET(name)     LOOK look ; \
     			  if (!graphAssFind (&MAGIC, &look)) \
		            messcrash ("%s can't find graph",name) ; \
                          if (look->magic != MAGIC) \
                            messcrash ("%s received a wrong pointer",name)

/**************************************************************/

static int tagColour = BLACK ;

static void colourBox (int box, BS bs, BOOL isActive)
{
  if (!bs) return ;

  if (isActive)
    {
      if (bs->size & (MODEL_FLAG | UNIQUE_FLAG))
	{
	  if (bs->size & SUBTYPE_FLAG)
	    graphBoxDraw (box, WHITE, GREEN) ;
	  else
	    graphBoxDraw (box, WHITE, BLUE) ;
	}
      else if (bsIsTag(bs))
	graphBoxDraw (box, WHITE, tagColour) ;
      else
	graphBoxDraw (box, WHITE, BLACK) ;
    }
  else
    {
      if (bs->size & MODEL_FLAG)
	{
	  if (bs->size & SUBTYPE_FLAG)
	    graphBoxDraw (box, BLACK, LIGHTGREEN) ;
	  else if (bsIsTag(bs))
	    graphBoxDraw (box, tagColour, CYAN) ;
	  else
	    graphBoxDraw (box, BLACK, CYAN) ;
	}
      else
	{
	  if (bsIsComment(bs))
	    graphBoxDraw (box, BLACK, LIGHTGRAY) ;
	  else if (bsIsTag(bs))
	    graphBoxDraw (box, tagColour, WHITE) ;
	  else
	    graphBoxDraw (box, BLACK, WHITE) ;
	}
    }
}

static void colourActive (LOOK look, BOOL isOn)
{
  if (look->activebox > 0)
    colourBox (look->activebox, 
	       arr(look->content,look->activebox,BS), isOn) ;
}

/************************************/

static void clearFlags (BS bs)
{
  bs->size = 0 ;
  if (bs->right)
    clearFlags (bs->right) ;
  if (bs->down)
    clearFlags (bs->down) ;
}

static void setOnFlag (BS bs)
{
  if (bs->right)
    { if (countSiblings(bs->right) > 6  ||
	  (pickExpandTag (bs->right->key) && countSiblings(bs->right) > 3))
	bs->size |= ON_FLAG ;
      setOnFlag (bs->right) ;
    }
  if (bs->down)
    setOnFlag (bs->down) ;
}

static void setOnFlagModelOnly (BS bs)
{
  if (bs->right)
    { if ((bs->size & MODEL_FLAG) &&
	  (countSiblings(bs->right) > 6  ||
	   (pickExpandTag (bs->right->key) && countSiblings(bs->right) > 3)))
	bs->size |= ON_FLAG ;
      setOnFlagModelOnly (bs->right) ;
    }
  if (bs->down)
    setOnFlagModelOnly (bs->down) ;
}

static void getOnFlagSet (BS bs, BitSet bb)
{
  if (bs->right)
    { 
      if (!class(bs->key) && bs->key >= _Date)
	{
	  if (bs->size & ON_FLAG)
	    bitSet (bb, bs->key) ;
	  else
	    bitUnSet (bb, bs->key) ;
	}
      getOnFlagSet (bs->right, bb) ;
    }
  if (bs->down)
    getOnFlagSet (bs->down, bb) ;
}

static void setOnFlagSet (BS bs, BitSet bb)
{
  if (bs->right)
    { 
      if (!class(bs->key) && bs->key >= _Date)
	{
	  if (bit(bb,bs->key))
	    bs->size |= ON_FLAG ;
	  else
	    bs->size &= ~ON_FLAG ;
	}
      setOnFlagSet (bs->right, bb) ;
    }
  if (bs->down)
    setOnFlagSet (bs->down, bb) ;
}

/************************************/

BOOL treeDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  LOOK look=(LOOK)messalloc(sizeof(struct LOOKSTUFF));

  look->magic = MAGIC;
  look->key = key;
  look->justKnownTags = treeJustKnownTags ;

  if (pickType (key) != 'B' || !(look->objet = bsCreateCopy(key)))
    goto abort ;
  clearFlags (look->objet->root) ;
  setOnFlag (look->objet->root) ;
  
  tagColour = prefInt ("TAG_COLOUR_IN_TREE_DISPLAY") ; 
  if (!tagColour) tagColour = BLACK ;

  if (isOldGraph)
    {
      lookDestroy () ;
      graphClear () ;
      graphGoto (0,0) ;
      graphRetitle (messprintf("%s: %s", className(key),name (key))) ;
    }
  else
    { 
      if (!displayCreate (TREE)) 
	goto abort ;
      graphRetitle (messprintf("%s: %s", className(key),name (key))) ;
      graphRegister (DESTROY, lookDestroy) ;
      graphRegister (PICK, lookPick) ;
      graphRegister (KEYBOARD, lookKbd) ;
      graphRegister (RESIZE, lookRedraw) ;
      graphFreeMenu (lookMenu, treeMenu) ;
    }

  look->graph = graphActive() ;
  graphAssociate (&MAGIC, look) ;

  look->content = arrayCreate(32,BS);
  look->bs2box = assCreate () ;

  look->activebox = 0 ;
  look->updateMode = FALSE ;
  look->tagChooserMode = FALSE ;
  lookDraw (look) ;

  return TRUE ;
  
abort:
  bsDestroy (look->objet) ;
  messfree (look) ;
  return FALSE ;
}

/************************************************************/

static void lookDestroy (void)
{
  LOOKGET("lookDestroy") ;

  if (look->bsWait && (look->bsWait->size & ADD_KEY_FLAG))
    displayUnBlock() ;	
  detach1 (look) ;
  if (look->objet && isCacheModified(look->objet->x) &&
      ( !strcmp("Tag-chooser", name(look->key)) ||
      messQuery (messprintf("%s not saved, Save ?",name(look->key)))))
    { defuse (look->objet->root) ;
      bsSave (look->objet);
    }
  else
    bsDestroy (look->objet);

  arrayDestroy (look->content) ;
  assDestroy (look->bs2box) ;
  if (look->textStack)
    stackDestroy (look->textStack) ;
  if (look->longTextArray)
    arrayDestroy (look->longTextArray);

  look->magic = 0 ;
  messfree (look) ;

  graphAssRemove (&MAGIC) ;
}

/************************************************************/

static void lookRedraw (void)
{
  LOOKGET("lookRedraw") ;

  lookDraw(look) ;
}

/*********************************************************/

static void lookPick (int box)
{ BS bs ;
  LOOKGET("lookPick") ;

  if (autoPreserve)
    displayPreserve () ;

  if (!box)
    { look->activebox = 0 ;
      return ;
    }
  else if ( box >= arrayMax (look->content))
    {
      messerror("LookPick received a wrong box number");
      return ;
    }

  bs = arr(look->content,box,BS) ;

  if (!bs)
    /* timestamp boxes have NULL bs attached, so they are inactive */
    return;

  /* the text of the currently picked box */
  graphPostBuffer (drawBStext (bs)) ; /* mhmp 11/04/97 */

  bs = arr(look->content,look->activebox,BS);

  if (box == look->activebox) /* a second hit - follow it */
    {  if (bs->size & ON_FLAG) /* if hidden, show branch */
	{ bs->size ^= ON_FLAG ;
	  lookDraw(look) ;
	}
      else if (bs->key == _Pick_me_to_call &&
	       !look->updateMode && !look->tagChooserMode)
	externalDisplay(look->key) ;
      else if (bs->key == _Gel && !look->updateMode &&
	       !look->tagChooserMode)
	fpDisplay (look->key) ;
      else if (bs->key == _E_mail &&  !look->updateMode &&
	       !look->tagChooserMode)
	treeDispMailer (look->key) ;
      else if (class(bs->key))
	display (bs->key,look->key, 0) ;
      else if (bs->right && !look->updateMode &&
	       !look->tagChooserMode)
	{ bs->size ^= ON_FLAG ;
	  lookDraw(look) ;
	}
    }
  else                            /* change now */
    { 
      colourActive (look, FALSE) ;
      look->activebox = box ;
      colourActive (look, TRUE) ;
    }
}

/*****************/

static void lookKbd (int k)
{ char *vp ;
  BS bs ;
  int box ;
  int nn = 3 ;
  LOOKGET("lookKbd") ;

  if (isDisplayBlocked())
    return ;

  if (autoPreserve)
    displayPreserve () ;

  if (!look->activebox) return ;
  colourActive (look, FALSE) ;
  if (!(bs = arr(look->content,look->activebox,BS)))
    return ;

retry: /* may be a hidden bs->box is hit, try nn times */
  switch (k)
    {
    case LEFT_KEY :
      while (bs->up && bs->up->down == bs)
	bs = bs->up ;
      if (bs->up && bs->up->up)
	bs = bs->up ;
      break ;
    case RIGHT_KEY :
      if (bs->right)
	bs = bs->right ;
      break ;
    case DOWN_KEY :
      if (bs->down)
	bs = bs->down ;
      else while (bs->up && bs->up->down == bs)
	bs = bs->up ;
      break ;
    case UP_KEY :
      if (bs->up && bs->up->down == bs)
	bs = bs->up ;
      else while (bs->down)
	bs = bs->down ;
      break ;
    case SPACE_KEY :
      if (bs->right)
	bs = bs->right ;
      else if (bs->down)
	bs = bs->down ;
      else
	while (bs->up && bs->up->up)
	  { while (bs->up && bs->up->down == bs)
	      bs = bs->up ;
	    if (bs->up && !bs->up->up)
	      break ;
	    bs = bs->up ;
	    if (bs->down)
	      { bs = bs->down ;
		break ;
	      }
	  }
      break ;
    }

  if (assFind (look->bs2box, bs, &vp))
    {
      box = assInt(vp) ;
      colourBox (box, bs, TRUE) ;
      if (box != look->activebox &&
	  class(bs->key) && iskey(bs->key) == 2)
	display (bs->key,look->key,0) ;
      
      look->activebox = box ;
    }
  else if (nn-- >0) goto retry ;
}

/*************************************************************/

void treeUpdate (void)
{
  BitSet bb = 0 ;
  LOOKGET("lookMenu") ;

  if(!checkWriteAccess())
      return ;

  if (!KEYKEY(look->key))
    { messout ("Sorry, you can not update models this way") ;
      return ;
    }
  if (pickList[class(look->key) & 255 ].protected)
    { messout ("Sorry, this is a protected class, you cannot update it interactively") ;
      return ;
    }
  detach1 (look) ;
  bsDestroy (look->objet) ;
  if (!(look->objet = bsUpdate (look->key)))
    { look->objet = bsCreateCopy (look->key) ;
      lookDraw (look) ;
      return ;
    }
  cacheMark (look->objet->x) ;
  look->updateMode = TRUE ;

  clearFlags (look->objet->root) ;
  bsFuseModel (look->objet, 0, look->justKnownTags) ;
  setOnFlag (look->objet->root) ; /* RD put back after FuseModel */
  if (bbSetArray)   /* mieg, inertia in open tags */
    {
      bb = array(bbSetArray, class(look->key), BitSet) ;
      if (bb)
	setOnFlagSet (look->objet->root, bb) ;
    }

  look->activebox = 0 ;
  graphPop() ;
  lookDraw (look) ;
  graphRegister (PICK, updatePick) ;
  graphRegister (KEYBOARD, updateKbd) ;
  graphFreeMenu (lookMenu, updateMenu) ;
  graphHelp ("Update") ;
}

static void lookMenu (KEY k)
{
  BS bs ;
  BitSet bb = 0 ;
  LOOKGET("lookMenu") ;

  /* if at end then check gActive is correct first */
  if (autoPreserve)
    displayPreserve () ;

  switch ((int) k)
    {

    case 3 :                        /* save */
      if (look->bsWait)
	{ messout ("Finish what you have started first, or cancel") ;
	  return ;
	}
      defuse (look->objet->root) ;
      if (!bbSetArray)   /* mieg, inertia in open tags */
	bbSetArray = arrayCreate (256, BitSet) ;
      bb = array(bbSetArray, class(look->key), BitSet) ;
      if (!bb)
	bb =  array(bbSetArray, class(look->key), BitSet) = bitSetCreate(lexMax(0),0) ;
      getOnFlagSet (look->objet->root, bb) ;

      bsSave (look->objet) ;
      look->objet = bsCreateCopy (look->key) ;
      /* If the object exists only in the cache, and is empty, it will
	 evaporate, if so put it back there (see displayKey in newkey.c) */
      if (!look->objet)
	{ OBJ obj = bsUpdate(look->key) ;   /* to give it a cache entry */
	  bsDestroy (obj) ;
	  look->objet = bsCreateCopy (look->key) ;
	}  
      clearFlags (look->objet->root) ;
      setOnFlagSet (look->objet->root, bb) ;

      look->activebox = 0 ;
      look->updateMode = FALSE ;
      lookDraw (look) ;
      graphRegister (PICK, lookPick) ;
      graphRegister (KEYBOARD, lookKbd) ;
      graphFreeMenu (lookMenu, treeMenu) ;
      graphHelp ("Tree") ;
      break ;
      
    case 5 :
      if (!look->activebox)
	{ messout ("Please first select a node") ;
	  return ;
	}
      bs = arr(look->content,look->activebox,BS) ;
      if (iskey(bs->key) == 2 )
	{ 
	  if(pickType(bs->key) == 'B')
	    display(bs->key, look->key, TREE) ;
	  else
	    messout("Sorry, I cannot display %s as a tree", name(bs->key)) ;
	}
      else
	messout("No data associated to %s", name(bs->key)) ;
	
      break ;

    case 7 :		/* add comment */
      if (look->bsWait)
	{ messout ("Please finish first what you have started, or cancel") ;
	  return ;
	}
      if (!look->activebox)
	{ messout ("Sorry, you must first select a node to add to (single-click)") ;
	  return ;
	}
      bs = arr(look->content,look->activebox,BS) ;
      if (bs->size & MODEL_FLAG)
	{ messout (
 "You must add comments to the original object (black and white stuff)") ;
	  return ;
	}
      if (bsIsComment(bs))
	{ messout ("You can't add comments to a comment") ;
	  return ;
	}
      bs->size |= ADD_COMMENT_FLAG ;
      look->bsWait = bs ;
      look->handleWait = cacheStoreHandle (look->objet->x) ;
      lookDraw (look) ;
      break ;

    case 8 :		/* delete */
      if (look->bsWait)
	{ messout ("Finish what you have started first, or cancel") ;
	  return ;
	}
      if (!look->activebox)
	{ messout ("You must first select a node to prune (single-click)") ;
	  return ;
	}
      bs = arr(look->content,look->activebox,BS) ;
      if (bs->size & MODEL_FLAG)
	{ messout ("You can only delete parts of the original object") ;
	  return ;
	}
      look->objet->curr = bs ;
      look->objet->modCurr = (bs->bt) ? bs->bt->bsm : 0 ;

      defuse (look->objet->root) ;
      bsRemove (look->objet) ;
      bsFuseModel (look->objet, 0, look->justKnownTags) ;
      setOnFlagModelOnly (look->objet->root) ;

      lookDraw (look) ;
      break ;

    case 9 :		/* update */
      displayPreserve() ;
      treeUpdate () ;
      break ;

    case 10 :		/* cancel operation */
      if (!look->bsWait)
	return ;
      if (1 ||  /* mieg 2004/11/25 */
	  (look->bsWait->size & ADD_KEY_FLAG))
	displayUnBlock() ;		/* clears display block */
      look->bsWait->size &= ~WAITING_FLAGS ;   /* clears edit and fake node flags */
      look->bsWait = 0 ;
      look->handleWait = 0 ;
      lookDraw (look) ;
      break ;

    case 11:		/* edit */
      if (look->bsWait)
	{ messout ("Finish what you have started first, or cancel") ;
	  return ;
	}
      if (!look->activebox)
	{ messout ("You must first select a node to edit (single-click)") ;
	  return ;
	}
      bs = arr(look->content,look->activebox,BS) ;
      if (bs->size & MODEL_FLAG)
	{ messout ("You can only edit parts of the original object") ;
	  return ;
	}
      if (bsIsTag(bs))
	{ messout ("You can not edit tags - only object names or data") ;
	  return ;
	}
      look->handleWait = cacheStoreHandle (look->objet->x) ;
      look->bsWait = bs ;
      bs->size |= EDIT_NODE_FLAG ;
      lookDraw (look) ;
      break ;

    case 13 :
      graphPrint() ; break ;

    case 17:
      showTimeStamps = showTimeStamps ? FALSE : TRUE ;
      look->activebox = 0;	/* don't carry box-selection 
				   over to new display style */
      lookDraw (look) ;
      break ;

    case 18:
      showTagCount = showTagCount ? FALSE : TRUE ;
      look->activebox = 0;	/* don't carry box-selection 
				   over to new display style */
      lookDraw (look) ;
      break ;


    case 25 :			/* contract */
      if (!look->activebox || 
	  !(bs = arr(look->content,look->activebox,BS)))
	{ messout ("You must first select a node to edit (single-click)") ;
	  return ;
	}
      if (!(bs->size & ON_FLAG))
	{ bs->size |= ON_FLAG ;
	  lookDraw (look) ;
	}
      else
	messout ("That node was already contracted") ;
      break ;

    case 26 :			/* expand */
      if (!look->activebox || 
	  !(bs = arr(look->content,look->activebox,BS)))
	{ messout ("You must first select a node to edit (single-click)") ;
	  return ;
	}
      if (bs->size & ON_FLAG)
	{ bs->size ^= ON_FLAG ;
	  lookDraw (look) ;
	}
      else
	messout ("That node was already expanded") ;
      break ;

#ifdef JUNK_ACEDB_MAILER
    case 14:
      treeDump (look) ; break ;

    case 20 :
      displayPreserve() ; break ;

    case 21 :
      autoPreserve = TRUE ; break ;

    case 22 :
      autoPreserve = FALSE ; break ;

    case 23 :
      biblioKey (look->key) ; break ;

    case 24 :
      { KEYSET nks =  bsKeySet(look->key) , new = keySetCreate() ;
	KEY *kp ;
	Graph  g = graphActive () ;
	int j1 = 0 , j2 ;
	if (nks)
	  { j2 = keySetMax(nks) ; kp = arrp(nks,0,KEY) - 1 ;
	    while(kp++, j2--)
	      if (!pickXref(class(*kp)) && iskey(*kp) == 2)
		keySet(new, j1++) = *kp ;
	    keySetDestroy(nks) ;
	  }
	displayCreate(DtKeySet) ;
	graphRetitle(messprintf("Neighbours of %s", name(look->key))) ;
	keySetShow (new,0) ;
	keySetSelect () ;
	graphActivate (g) ;
      }
      break ;
      /* now that we can paste cut, just paste the adress to your favorite mailer */
    case 31 :  /* Mailer */
      { 
	OBJ  obj = bsCreate(look->key) ;
	char *mailAddress ;
	graphHelp ("Mailer") ;

	if (!obj)
	  break ;
	if (!bsGetData(obj, _E_mail, _Text,&mailAddress))
	 {
	   if ( messPrompt  ("Please specify an email address","","w"))
	     { mailAddress = freeword() ;
	       	       
	       { OBJ obj1 = bsUpdate(look->key) ;
		 if (!obj1)
		   goto fin ;
		   		 
		 if (!bsAddData(obj1, _E_mail, _Text, mailAddress))
		   { bsDestroy(obj1) ;
		     messout ("Sorry, I can't save your address, no E_mail Tag in this class") ;
		     goto fin ;
		   }
		 bsSave(obj1) ;
	       }
	     }
	   else
	     goto fin ;
	 }
	else
	acedbMailer(look->key, 0,0) ;
      fin:
	bsDestroy(obj) ;
	break ;
      }
#endif
	
    case 98 : 
      graphHelp("Tree");
      { int nn = 0 ;
      Stack h = stackCreate (50) ;
      KEYSET ks1 = keySetCreate () ;
      KEY k1 ;
      
      if (look->activebox) 
	{
	  bs = arr(look->content,look->activebox,BS) ;
	  while (bs)   /* look upward for a tag */
	    { keySet (ks1, nn++) = bs->key ;
	    if (!class (bs->key) && bs->key >= 50)
	      break ;
	    bs = bs->up ;
	    } 
	}
      
      pushText (h, messprintf("%s_%s",graphHelp(0),className(look->key))) ;
      while (nn--)
	{ catText (h, "_") ;
	k1 = keySet(ks1,nn) ;
	catText (h, class(k1) ? className(k1) : name(k1)) ;
	}
      
      helpOn (stackText (h, 0)) ;
      stackDestroy (h) ;
      keySetDestroy (ks1) ;
      }
      break ;

    case 99 :
      graphDestroy () ; break ;
    }
}

/****************************************************/

static void treeDump (LOOK look)
{
  FILE *fil ;
  static char dirName[DIR_BUFFER_SIZE] = "" ;
  static char fileName[FIL_BUFFER_SIZE] = "object" ;

  if (! look) return ;  /* to please te compiler */
  if (!lexIsKeyVisible (look->key))
    { messout("Sorry, I don't know how to dump %s\n", name(look->key));
      return;
    }
  if (!(fil = filqueryopen (dirName, fileName, "ace", "w",
			    "Where do you wish to dump this object ?")))
    { 
      return ;
    }
  fprintf (fil,"// data dumped from tree display\n\n") ;
  dumpKey (look->key, fil, 0) ;
  filclose (fil) ;
  messout ("Object %s written to %s/%s.ace\n", 
	   name(look->key), dirName, fileName) ;
}

/*************************************************************/
/*************** drawing package *****************************/
/*******************************************************/

static void attach (void) ;
static void detach (void) ;
static void attachRecursively (void) ;
static void attachHelp (void) ;

static MENUOPT attachMenu [] = {
  {attach, "Attach"},
  {detach, "Detach"},
  {attachRecursively, "Attach Recursively"},
  {attachHelp, "Help on attach"},
   {0, 0}
} ;

static LOOK drawLook ;
static int xmax, ymax ;
static BOOL isModel ;

/******************* attach subpackage ***********************/

static void doDetach(BS bs)
{
  if (bs->right)
    doDetach (bs->right) ;
  if (bs->down)
    doDetach (bs->down) ;
  if (bs->right && (bs->right->size & ATTACH_FLAG))
    { bsTreePrune (bs->down) ; bs->down = 0 ; }
}

static void detach1 (LOOK look)
{ int i ;
  Array a = look->attachedObjs ;
  
  if (arrayExists(a))
    { doDetach (look->objet->root) ;
      i = arrayMax (a) ;
      while (i--)
	bsDestroy (array (a, i, OBJ)) ;
    }
  arrayDestroy (look->attachedObjs) ;
}
	
static void detach (void)
{ int box ;
  BS bs ;
  LOOKGET ("detach") ;
  
  box = look->activebox ;
  if (!box)
    { messout ("Please first select a node") ;
      return ;
    }
  else if ( box >= arrayMax (look->content))
    { messerror("Sorry, detach received a wrong box number");
      return ;
    }

  bs = arr(look->content,box,BS) ;

  if (bs && bs->right)
    doDetach (bs->right) ;
 /* do not detach down */
  if (bs && bs->right && (bs->right->size & ATTACH_FLAG))
    bs->right = 0 ;

  lookDraw (look) ;
}

static void doAttach (OBJ mainObj, BS bs, BOOL whole, Array a, KEYSET ks, 
		      int recursive, char *question)
{ OBJ obj ;
  int i ;
  int dummy ;
  KEY key ;
  KEYSET ks1 = 0 ;
  BS bs1 = bs , bs2 = bs ;

  while (bs1->up && 
	 ( pickType(bs1->key) != 'B' ||
	  !class(bs1->key) || class(bs1->key) == _VComment ||
	  class(bs1->key) == _VUserSession))
    bs1 = bs1->up ;
	
  if (bs->right) goto more ;

  ks1 = queryKey (bs1->key, question ? question : "IS *") ;
  for (i = 0 ; i < keySetMax (ks1) ; i++)
    { key = keySet (ks1, i) ;
      if ((recursive < 2 || !keySetFind (ks, key, &dummy)) &&
	  (obj = bsCreateCopy (key)))
      { 
	if (whole || queryFindLocalise (obj, key, 0))
	  { 
	    bs1 = bsTransferBranch (mainObj, obj, whole) ;
	    if (bs2 == bs)
	      bs2->right = bs1 ;
	    else  
	      bs2->down = bs1 ;

	    bs1->up = bs2 ;
	    bsDestroy (obj) ; /* destroy the now stripped tmp obj */
	  }
	else
	  { bsDestroy (obj) ;
	    goto more ;
	  }
	clearFlags (bs1) ;
	setOnFlag (bs1) ;
	bs1->size |= ATTACH_FLAG ;
	array (a, arrayMax (a), OBJ) = obj ;	
	if (ks && recursive > 1)
	  keySetInsert (ks, key) ;
      }
    }
  keySetDestroy (ks1) ;
 more:
  if (bs->right)
    if (recursive || !(bs->right->size & ATTACH_FLAG))
      doAttach (mainObj, bs->right, whole, a, ks, recursive ? recursive + 1 : 0, question) ;
  if (bs->down)
    doAttach (mainObj, bs->down, whole, a, ks, recursive, question) ;
}

static void attach1 (int box, BOOL recursive) 
{ 
  BS bs, bs1, down ;
  int type, tcl ; 
  KEY tag = 0 ;
  char *cp ;
  Stack sta = 0 ;
  KEYSET ks = 0 ;
  LOOKGET("attach") ;

  if (!box)
    { messout ("Please first select a node") ;
      return ;
    }
  else if ( box >= arrayMax (look->content))
    { messerror("Sorry, attach received a wrong box number");
      return ;
    }

  bs = arr(look->content,box,BS) ;
  bs1 = bs ;
  while (bs1->right && !class(bs1->key)) bs1 = bs1->right ;
  if (pickType(bs1->key) != 'B')
    { messout ("First object to the right is not a Tree. I can't attach, sorry") ;
      return ;
    }
  sta = stackCreate (50) ;
  if (treeChooseTagFromModel (&type, &tcl, class(bs1->key), &tag, sta, 0))
    { if (!arrayExists(look->attachedObjs))
	look->attachedObjs = arrayCreate (12, OBJ) ;
      if (recursive)
	ks = keySetCreate () ; /* against recursive attach runaway */
      down = bs->down ;
      bs->down = 0 ;
      cp = stackText (sta, 0) ;
      if (!strcmp ("WHOLE", cp))
	doAttach (look->objet, bs, 1, look->attachedObjs, ks, recursive, 0) ;
      else
	{ queryFindLocalise (0, 0, cp) ; /* Initialise queryFind */
	  doAttach (look->objet, bs, 0, look->attachedObjs, ks, recursive, cp) ;
	}
      bs->size &= ~ON_FLAG ;
      stackDestroy (sta) ;
      keySetDestroy (ks) ;
      bs->down = down ;
    }
  lookDraw (look) ;
}

static void attach (void) 
{ LOOKGET ("attach") ;

  attach1(look->activebox, FALSE) ;
}

static void attachRecursively (void) 
{ LOOKGET ("attachRecursivelly") ;
 
  attach1(look->activebox, TRUE) ;
}

static void attachQuery (OBJ obj, BS bs, KEY qr)
{ char *cp0, *cp = name (qr) ;
  LOOKGET ("attachQuery") ;

  if (!arrayExists(look->attachedObjs))
    look->attachedObjs = arrayCreate (12, OBJ) ;

  cp0 = cp ;
  cp += strlen (cp0) ;
  while (*(cp - 1 ) != ';' && cp > cp0) cp-- ;
  if (!strcmp ("WHOLE", cp))
    doAttach (obj, bs, 1, look->attachedObjs, 0, 0, 0) ;
  else
    { queryFindLocalise (0, 0, cp) ; /* Initialise queryFind */
      doAttach (obj, bs, 0, look->attachedObjs, 
		0, FALSE, cp0) ;
    }
}

static void attachHelp (void) 
{ helpOn ("attach") ;
}

/*******************************************************************/
static FREEOPT directTreeActionMenu[] = {
  {  4, "directTreeAction menu"},
  {'t', "Show as text"},
  {'g', "Show graphics"},
  {'z', ""},
  {'d', "Eliminate this node and everything on its right"}
} ;

static void directTreeAction (KEY key, int box)
{
  BS bs, bs1 ;
  Stack s, s1, t ;
  OBJ obj ;
  char timeBuf[25] ;
  LOOKGET ("directTreeAction") ;
  
  if (arrayExists(look->content) && box > 0 &&
      box < arrayMax(look->content))
    bs = array (look->content, box, BS) ;
  else
    return ;
  if (!checkWriteAccess ())
    return ;
  if (pickList[class(look->key) & 255 ].protected)
    { messout ("Sorry, this is a protected class, you cannot update it interactively") ;
    return ;
    }
  if (!class(look->key))
    return ;
  obj = bsUpdate(look->key) ;
  if (!obj) return ;
  bsSave(obj) ;
  switch (key)
    {
    case 't':
      displayPreserve () ;
      display (bs->key, look->key, TREE) ;
      break ;
    case 'g':
      display (bs->key, look->key, 0) ;
      break ;
    case 'd': /* eliminate */
      s = stackCreate (64) ; s1 = stackCreate (1024) ; 
      t = stackCreate (64) ;
      push (s, bs, BS) ; /* but not on t */
      bs1 = bs ; bs = bs->up ;
      while (bs)
	{ 
	  if (bs->right == bs1) 
	    { push (s, bs, BS) ; push (t, bs, BS) ; }
	  bs1 = bs ;
	  bs = bs->up ; 
	}
      pushText (s1, messprintf("%s %s\n-D  ", className(look->key),
		freeprotect(name(look->key)))) ;	  
      bs = pop (s, BS) ; /* jump root node */
      while (!stackAtEnd(s))
	{
	  bs = pop (s, BS) ;
	  if (class(bs->key) == _VUserSession)
	    continue ;
	  if (class(bs->key) == _VComment)
	     {
	       catText (s1, " -C ") ;
	       catText (s1,  freeprotect (name(bs->key))) ;
	     }
	  else if (class(bs->key))
	    catText (s1,  freeprotect (name(bs->key))) ;
	  else if (bs->key <= _LastC)
	    catText (s1,  freeprotect (bsText(bs))) ;
	  else if (bs->key == _Int)
	    catText (s1, messprintf ("%d", bs->n.i)) ;
	  else if (bs->key == _Float)
	    catText (s1, messprintf ("%g", bs->n.f)) ;
	  else if (bs->key == _DateType)
	    catText (s1, timeShow (bs->n.time, timeBuf, 25)) ;
	  else
	    catText (s1, name (bs->key)) ;
	  catText (s1, " ") ;
	} 
      catText (s1, "\n") ;
      /* now reestablish upper node to counteract rd's bizare system  */
      bs = pop (t, BS) ; /* jump root node */
      while (!stackAtEnd(t))
	{
	  bs = pop (t, BS) ;
	  if (class(bs->key) == _VUserSession)
	    continue ;
	  if (class(bs->key) == _VComment)
	     {
	       catText (s1, " -C ") ;
	       catText (s1,  freeprotect (name(bs->key))) ;
	     }
	  else if (class(bs->key))
	    catText (s1,  freeprotect (name(bs->key))) ;
	  else if (bs->key <= _LastC)
	    catText (s1,  freeprotect (bsText(bs))) ;
	  else if (bs->key == _Int)
	    catText (s1, messprintf ("%d", bs->n.i)) ;
	  else if (bs->key == _Float)
	    catText (s1, messprintf ("%g", bs->n.f)) ;
	  else if (bs->key == _DateType)
	    catText (s1, timeShow (bs->n.time, timeBuf, 25)) ;
	  else
	    catText (s1, name (bs->key)) ;
	  catText (s1, " ") ;
	}
      parseBuffer (stackText (s1,0), 0) ;
      stackDestroy (s) ;
      stackDestroy (s1);
      bsDestroy (look->objet) ;
      look->objet = bsCreateCopy(look->key) ; 
      clearFlags (look->objet->root) ;
      setOnFlag (look->objet->root) ;

      lookDraw (look) ;
      break ;
    default: break ;
    }
}

/*******************************************************************/

static int treeClassComplete (char *cp, int len)
{
  if (drawLook->classe)
    return ksetClassComplete (cp, len, drawLook->classe) ;
  else
    return 0 ;
}

static int drawTextEntry (int classe, int len, 
			  int x, int y, void (*fn)(char*))
{
  int box = 
    graphCompScrollEntry (treeClassComplete, 
			  stackText(drawLook->textStack, 0), 
			  len, 32, x, y, fn) ;
  array (drawLook->content, box, BS) = 0 ;
  if (x + 32 > xmax)	/* RD 9402 - 32 chars appear on screen */
    xmax = x + 32 ;
  drawLook->classe = classe ;
  return y+1 ;
}

static BOOL dontExpand = FALSE ;

static char *expandText (KEY key)
{
  KEY tag = pickExpandTag(key) ;
  OBJ obj ;
  char *text = 0 ;

  if (tag && (obj = bsCreate(key)))
    { if (bsGetKeyTags (obj,tag,0))
	text = drawBStext (obj->curr) ;
      bsDestroy (obj) ;
    }

  return text ;
}

static char* drawBStext (BS bs)
{
  char *text ;
  static char timeBuf[25] ;

  if (bs->size & SUBTYPE_FLAG)
    return name (bs->key - 1) ;
  else if (bs->size & MODEL_FLAG)
    return KEYKEY (bs->key) ? name (bs->key) : messprintf ("?%s", className(bs->key)) ;  /* use data tag names in model sections */
  else if (!dontExpand && (text = expandText(bs->key)))
    return text ;
  else if (isModel)
    return name (bs->key) ;
  else if (bs->key <= _LastC)
    return bsText(bs) ;
  else if (bs->key == _Int)
    return messprintf ("%d", bs->n.i) ;
  else if (bs->key == _Float)
    return messprintf ("%g", bs->n.f) ;
  else if (bs->key == _DateType)
    return timeShow (bs->n.time, timeBuf, 25) ;
  else
    return name (bs->key) ;
}

static int countSiblings (BS bs)
{ int n = 0 ;
  
  while (bs)
    { if (class(bs->key) != _VUserSession)
	n++ ; 
      bs = bs->down ; 
    }
  return n ;
}

static void drawTriangle(BS bs, int x, int y)
{
  int n = countSiblings(bs->right) ;
  graphColor(BLUE) ;
  graphText(messprintf("  -----> %d ",n), x, y) ;
  graphColor(BLACK) ;
}

static void treeDispTagCount (int cl, BS bs, BOOL mask)	/* recursive routine */
{
  
  while (bs)
    {
      if (!class (bs->key) && bs->key >= _Date)
	{
	  bs->tagCount = (0xffffffff & bIndexTagCount (cl, bs->key)) ;
	  /* -1, dont know in bindex
	   *  0, not evaluated
	   */
	}
      else if (bs->key != _UNIQUE)
	mask = TRUE ;
      if (mask) 
	bs->tagCount = -1 ;
      if (bs->right)
	treeDispTagCount (cl, bs->right, mask) ;
      bs = bs->down ;
    }
} /* treeDispTagCount */


static int drawBS (OBJ obj, BS bs, BS bsm, int x, int y)	/* recursive routine */
{
  int yMe = y ;
  int xPlus = 0 ;
  char *text, *cp, *vp ;
  int box ;
  int oldTextFormat = PLAIN_FORMAT ;

  bsModelMatch (bs, &bsm) ; /* will fail on Time stamp,
			       ignore returned value */

  if (bs->size & ATTACH_FLAG && bs->up && bs->up->key == bs->key)
    goto suite ; /* because, i don t draw the attached node itself */

  if (bs->size & ADD_DATA_FLAG)
    { y = drawTextEntry (0, bs->key < _LastC ? 2500 : 32, x, y, addData) ;
    }
  else if (bs->size & ADD_KEY_FLAG)
    { y = drawTextEntry (class(bs->key), 2500, x, y, addKeyText) ;
      displayBlock (addKey, "Or you can type in the name.  "
		    "If you press the TAB key you will autocomplete "
		    "or get a list of possible entries.") ;
    }

  text = drawBStext (bs) ;

  yMe = y ;  /* must do this here again: y can change */

  if (bs->size & EDIT_NODE_FLAG)		/* edit entry */
    { stackTextForceFeed (drawLook->textStack, strlen(text) + 2560) ;
      stackClear (drawLook->textStack) ;
      pushText (drawLook->textStack, text) ;
      y = drawTextEntry (class(bs->key), strlen(text) + 2500,
			 x, y, editEntry) ;
    }
  else if (!drawLook->updateMode && 
	   !drawLook->tagChooserMode &&
	   class(bs->key) == _VLongText &&
	   KEYKEY(bs->key) && iskey(bs->key) == 2)
    {
      text = messprintf ("<see below %d>", arrayMax(drawLook->longTextArray)+1) ;
      graphText (text, x, yMe++) ;
      if (strlen(text) > xPlus)
	xPlus = strlen(text) ;
      array(drawLook->longTextArray,
	    arrayMax(drawLook->longTextArray), BS) = bs ;
    }
  else
    {
      box = graphBoxStart() ;
      array (drawLook->content, box, BS) = bs ;
      vp = (char *)0 + box ;
      assInsert (drawLook->bs2box, bs, vp) ;
      if ((iskey (bs->key) == 2 && class(bs->key) != _VText) ||
	  ((bs->size & ON_FLAG) && bs->right) )
	oldTextFormat = graphTextFormat(BOLD) ;
      if (bs->key == _Greek)
	oldTextFormat = graphTextFormat(GREEK) ;
      uLinesText (text, LINE_LENGTH) ;
      if ((cp = uNextLine(text)))
	{ graphText (cp,x,yMe++) ;	/* write out this node */
	  if (strlen(cp) > xPlus)
	    xPlus = strlen(cp) ;
	  while ((cp = uNextLine(text)))	/* indent following lines */
	    { graphText (cp,x+2,yMe++) ;
	      if (strlen(cp)+2 > xPlus)
		xPlus = strlen(cp)+2 ;
	    }
	}

      if ((bs->size & ON_FLAG) && bs->right )
	drawTriangle(bs, x+xPlus, yMe-1) ;
      if ((iskey (bs->key) == 2 && class(bs->key) != _VText) ||
	  ( (bs->size & ON_FLAG) && bs->right  ) ||
	  bs->key == _Greek )
	graphTextFormat(oldTextFormat) ;
      graphBoxEnd() ;
      graphBoxInfo (box, bs->key, 0) ;
      if (box == drawLook->activebox)
	colourBox (box, bs, TRUE) ;
      else
	colourBox (box, bs, FALSE) ;
      if (!drawLook->updateMode && class(bs->key) != _VUserSession)
	{
	  BS bs1 = drawLook->attachedObjs ? bs : 0 ;
	  while (bs1)
	    { if (bs1->size & ATTACH_FLAG) break ; bs1 = bs1->up ; }
	  /* for 4.7h don't provided the [eliminate this node ...] action */
	  if (!bs1)
	    graphBoxFreeMenu (box, directTreeAction, directTreeActionMenu) ;
	  /* yes i want it   */
	}
      if (bs->right && class(bs->key) &&   /* try to align things to the right of columns of keys */
	  (bs->down || (bs->up && bs->up->down == bs)))
	xPlus = ((xPlus + 3)/4) * 4 ;
    }

 suite:
  xPlus += x ;
  if (xPlus > xmax)
    xmax = xPlus ;
  xPlus += 2 ; /* spacing */

  if (bs->size & ADD_COMMENT_FLAG)	/* add comment box */
    y = drawTextEntry (_VComment, 2500, xPlus, y, addComment) ;
  else if (!dontExpand && showTimeStamps)
    { if (bs->timeStamp)
	{
	  if ( (bs->size & ON_FLAG) && bs->right ) /* if collapsed */
	    xPlus += 9+2;	/* length of drawTriangle()-string+2 */

	  box = graphBoxStart() ;
	  graphText (name(bs->timeStamp), xPlus, y) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box, BLACK, LIGHTGRAY) ;
	  array (drawLook->content, box, BS) = 0 ;
				/* do not attach box to a bs */
	}
      y = yMe ; xPlus = x+2 ;	/* force vertical layout */
    }
  else if (!dontExpand && showTagCount && bsIsTag(bs) && bs->tagCount >= 0)
    { 
      if ( (bs->size & ON_FLAG) && bs->right ) /* if collapsed */
	xPlus += 9+2;	/* length of drawTriangle()-string+2 */
      
      box = graphBoxStart() ;
      graphText (messprintf ("%d", bs->tagCount), xPlus, y) ;
      graphBoxEnd () ;
      graphBoxDraw (box, BLACK, bs->tagCount > 0 ? GREEN : RED) ;
      array (drawLook->content, box, BS) = 0 ;
				/* do not attach box to a bs */

      y = yMe ; xPlus = x+2 ;	/* force vertical layout */
    }
  if (bs->right && !(bs->size & ON_FLAG))
    { /* mieg, may 30 96, try to lay out in a more vertical way
         without breaking inside data values 
	 */
      static int mieglayout = 0;
      if (!mieglayout)
	mieglayout = prefValue("HORIZONTAL_TREE") ? 1:-1;
      if ((mieglayout == -1) && 
	  bsIsTag(bs) && bs->key != _XREF && (bsIsTag(bs->right) || 
					      ( class (bs->right->key) && countSiblings (bs->right) > 1)))
	{
	  if (yMe > y)
	    y = yMe ;	
	  y = drawBS (obj, bs->right, bsm ? bsModelRight(bsm) : 0, x + 2, y) ;  /* below, at x + 2 */
	}
      else
	y = drawBS (obj, bs->right, bsm ? bsModelRight(bsm) : 0, xPlus, y) ;  /* to the right at same y */
    }
  
  if (!bs->right && !(bs->size & ON_FLAG) &&
      bsm && bsm->right &&
      class (bsm->right->key) == _VQuery &&
      KEYKEY (bsm->right->key) != 1)  /* a hack, because treedisp overloads 1 as #type */
    { attachQuery (obj, bs, bsm->right->key) ;
      if (bs->right)
	y = drawBS (obj, bs->right, bsm ? bsModelRight (bsm) : 0 , xPlus, y) ;  /* to the right at same y */
    }

  if (yMe > y)
    y = yMe ;

  if (bs->down)
    y = drawBS (obj, bs->down, bsm, x, y) ;		/* below at new y location */

  return y ;
}

static Array tableList = 0 , class2table = 0 ;

static FREEOPT* initTableList(int nn)
{ FILE *f ;
  int cc , level ;
  char *cp, *cq ; 
  BOOL inside = FALSE ;
  int n1 = 0, n2 = 1 , nt = 0 , k ;

  if (arrayExists(tableList))
    goto done ;

  tableList = arrayCreate(10, FREEOPT) ;
  class2table = arrayCreate(64, int) ;

  if (!filName ("wspec/table.menu", "wrm", "r"))
    return 0;
  if (!(f = filopen ("wspec/table.menu", "wrm", "r")))
    return 0 ;

  level = freesetfile(f,"") ;
  
  while(freecard(level))
    {
      cq = freepos() ; cp = freeword() ;
      if (!cp)
	{ if (inside)
	    { inside = FALSE ;
	      array(tableList,n1, FREEOPT).key = nt ;
	      array(tableList,n2++, FREEOPT).key = 0 ;
	    }
	}
      else
	{
	  if (!inside)
	    { if ((cc = pickWord2Class(cp)))
		{ inside = TRUE ;
		  nt = 0 ;
		  array(class2table, cc, int) = n1 = n2 ;
		  array(tableList,n2++, FREEOPT).text = pickClass2Word(cc) ;
		}
	    }
	  else
	    { nt++ ;
	      array(tableList,n2, FREEOPT).key = nt ;
	      array(tableList,n2++, FREEOPT).text = cp = messalloc(strlen(cq) + 1 ) ;
	      strcpy(cp, cq) ;
	    }
	}
    }

 done:
    { k = array(class2table, nn, int) ;
      if (k)
	return arrayp(tableList, k, FREEOPT) ;
      else
	return  NULL ;
    }
}

static void tables(void)
{ FREEOPT* ff = 0 ;
  int i ; KEY k ;
  LOOKGET("tables") ;

  i = array(class2table,class(look->key), int) ;
  if (i)
    ff = arrayp(tableList,i,FREEOPT) ;
  if (ff && graphSelect(&k, ff))
    { char * cp = freekey2text(k, ff) , *cq ;
      
      cq = cp ;  /* last word is file name */
      while(*cq) cq++ ;
      while(cq > cp && *(cq - 1) != ' ') cq-- ;
      
      spreadTableDisplay (cq, name(look->key)) ;
    }
}

         /* same but with a different menu interface */
static void tablePick(KEY k, int box)
{ FREEOPT* ff = 0 ;
  int i ; 
  LOOKGET("tables") ;

  i = array(class2table,class(look->key), int) ;
  if (i)
    ff = arrayp(tableList,i,FREEOPT) ;
  if (ff)
    { char * cp = freekey2text(k, ff) , *cq ;
      
      cq = cp ;  /* last word is file name */
      while(*cq) cq++ ;
      if (cp < cq) cq-- ;
      while(cq > cp && *cq == ' ') *cq-- = 0 ;
      while(cq > cp && *(cq - 1) != ' ') cq-- ;
      
      spreadTableDisplay(cq, name(look->key)) ;
    }
}
/*  myLongTextShow  is not used, you have my blessing to delete it, mieg.
static void myLongTextShow (void)
{
  int box = graphBoxAt (graphEventX, graphEventY, 0, 0) ;
  BS bs ;
  LOOKGET ("myLongTextShow") ;
 
  bs = arr(look->content, box, BS) ;
  if (bs && bs->key)
    display (bs->key, look->key, 0) ;
}
*/
static float drawLongText (float ymax)
{
  int box ;
  BS bs ;
  int i, x ;
  char *cp ;
  Stack text, path ;

  if (!arrayMax (drawLook->longTextArray))
    return ymax ;

  path = stackCreate(64) ;
  for (i = 0; i < arrayMax (drawLook->longTextArray); i++)
    { bs = arr(drawLook->longTextArray, i, BS) ;
      if (!(text = stackGet (bs->key)))
	continue ;
				/* draw a spacer line */
      ymax += 1.0 ;
      graphLine (0, ymax, 1000, ymax) ;
      ymax += 1.0 ;
				/* draw the path to this node */
      stackClear (path) ;
      do
	{ push (path, bs, BS) ;
	  while (bs->up->down == bs) bs = bs->up ;
	  bs = bs->up ;
	} while (bs->up) ;	/* until reach root */
      x = 0 ; 
      cp = messprintf ("%d.  ", i+1) ;
      graphText (cp, x, ymax) ;
      x += strlen(cp) ;
      while (!stackEmpty (path))
	{ cp = drawBStext (bs = pop (path, BS)) ;
	  if (stackEmpty (path))
	    { box = graphBoxStart() ;
	      array (drawLook->content, box, BS) = bs ;
	      graphTextFormat (BOLD) ;
	    }
	  graphText (cp, x, ymax) ;
	  x += strlen(cp) + 1 ;
	}
      graphBoxEnd () ;
      graphTextFormat (PLAIN_FORMAT) ;

      ymax += 2 ;
      if (xmax < 12) xmax = 12 ;
      stackCursor (text, 0) ;
      while ((cp = stackNextText(text)))
	{ 
	  char *cq ;
	  uLinesText (cp, LINE_LENGTH) ;  
	 while ((cq = uNextLine(cp)))
	   {
	     graphText (cq,1,ymax++) ;
	     if (1 + strlen (cq) > xmax)
	       xmax = 1 + strlen(cq) ;
	   }
	}
      stackDestroy (text) ;
    }

  stackDestroy (path) ;
  return ymax ;
}

static void knownTagsButton (void)
{
  LOOKGET ("knownTagButton") ;
  treeJustKnownTags = look->justKnownTags = look->justKnownTags ? FALSE : TRUE ;
  defuse (look->objet->root) ;
  bsFuseModel (look->objet, 0, treeJustKnownTags) ;
  lookDraw (look) ;
}

static void saveButton (void)
{
  lookMenu (3) ;
}

static void updateButton (void)
{
  lookMenu (9) ;
}

static void biblioButton (void)
{
  LOOKGET ("biblioButton") ;
  biblioKey (look->key) ;
}

static void lookDraw (LOOK look)
{ FREEOPT* ff = 0 ;
  int box, dx, offSet = 59 ;
  float x1, y1, x2, y2, ddx = 1 ;

  if (!graphActivate (look->graph))
    return ;
  if (look->activebox)
    graphBoxDim (look->activebox, &x1, &y1, &x2, &y2) ;
  else 
    y1 = -1 ;

  graphClear () ;
  arrayMax(look->content) = 0 ;
  assClear (look->bs2box) ;

  graphBoxInfo (0, look->key, 0) ;
  box = graphBoxStart() ;
  array (look->content, box, BS) = look->objet->root ;
  graphText (name(look->key),0,0) ;
  graphBoxEnd () ;

  isModel = class(look->key) == _VModel  || !KEYKEY(look->key)  ?
    TRUE : FALSE ; 

  if (showTagCount)
    {
      int cl = isModel ? pickWord2Class (name(look->key)+1) : class (look->key) ;
      if (cl)
	treeDispTagCount (cl, look->objet->root->right, FALSE) ;
    }

  if (look->objet->root && look->objet->root->right)
    { drawLook = look ;
      drawLook->textStack = stackReCreate (drawLook->textStack, 2560) ;
      drawLook->longTextArray
	= arrayReCreate (drawLook->longTextArray, 8, BS);
      xmax = 0 ;
      dontExpand = look->updateMode  || look->tagChooserMode ;
      ymax = drawBS (look->objet, look->objet->root->right, 
		     bsModelRoot(look->objet)->right, 2, 4) ;
      ymax = drawLongText (ymax) ; /* show any long text that occurred */
    }
  else
    { xmax = 1 ; 
      ymax = 2 ;		/* for buttons */
    }

if (xmax <  offSet) xmax =  offSet + 1 ;
#if !defined(macintosh)
  dx = 4 + ddx ;
  graphButton("Quit", graphDestroy, offSet - dx, .4) ;
  offSet -= dx + ddx ;
#endif

  if (!look->updateMode && !look->tagChooserMode)
    {
      if (writeAccessPossible())
	{
	  dx = 6 + ddx ;
	  box = graphButton ("Update", updateButton, offSet - dx, .4) ;
	  offSet -= dx + ddx ; ;
	}
      dx = 9 + ddx ;
      box = graphButton ("Attach...", attach, offSet - dx, .4) ;
      if (xmax < 52) xmax = 52 ;
      graphBoxMenu (box, attachMenu) ; 
      offSet -= dx + ddx ;
    }
  else
    {
      dx = 4 + ddx ;
      box = graphButton ("Save", saveButton, offSet - dx, .4) ;
      offSet -= dx + ddx ; ;

      if (treeJustKnownTags)
	{
	  dx = 13 + ddx ;
	  box = graphButton ("Show all tags", knownTagsButton, offSet - dx, .4) ;
	}
      else
	{
	  dx = 18 + ddx ;
	  box = graphButton ("Limit to known tags", knownTagsButton, offSet - dx, .4) ;
	}
      offSet -= dx + ddx ;
    }

  if ((ff = initTableList(class(look->key))))
    { 
      dx = 9 + ddx ;
      box = graphButton("Tables...", tables, offSet - dx, 3.4) ;
      offSet -= dx + ddx ; ;
      graphBoxFreeMenu(box, (FreeMenuFunction) tablePick, ff) ;
    }

  if (biblioKeyPossible(look->key))
    {
      dx = 6 + ddx ;
      graphButton("Biblio", biblioButton,  offSet - dx, .4) ; 
      offSet -= dx + ddx ; 
    }

  graphTextBounds (xmax,ymax) ;
  if (y1 > 0)
    graphGoto (x1, y1) ;
  graphRedraw() ;
  pickRememberDisplaySize ("TREE") ;

  if (look->activebox >= arrayMax(look->content))
    look->activebox = 0 ;

  isModel = FALSE ;
}

/************************************************************/
/************************************************************/

/* Updating:
   The plan is to build a new tree containing the union of the model
   and the current object.  We will rebuild every time an action changes the
   structure.  
   Parts of the model not yet represented in the object will be coloured
   blue.  The XREF, UNIQUE etc bits of the model will not be shown. Class
   slots (?Class) will appear in cyan at the bottom of each column of objects,
   so you can add new members (or replace if unique).
   slots (#Class) will not be shown.
   If you double-click on a cyan tag it will be added to your object.
   If you double click on a (cyan) ?Class then you will be prompted for a name,
   but you can also pick something from that class elsewhere in the system.  If
   you give a name that is a new lex entry it will ask for confirmation before
   adding.
   When displayBlocked you can not display any new objects, nor get the biblio
   etc.

   We will use bs->size to hold the information about whether an entry is
   in the model or not.  Adding and subtracting the model is done by two
   routines bsFuseModel() in bssubs.c and defuse() below.
   Codes for ->size are:
   	1	part of model (painted cyan)
	2	unique part of original tree
	4	fake node to add data (TextEntry)
	8	fake node to add comment (TextEntry)
	16	fake node to add key (TextEntry and DisplayBlock)
	32      edit current node - rather than enter a new one
	64	special node to add Type (TextEntry and DisplayBlock)
	   all from 4 on are in fact mutually exclusive - there
	   should be at most one at any one time (set to bsWait).
       128	ON_FLAG  to contract/uncontract
       256	ATTACH_FLAG
       512	new node this session
   We disable the keyboard routine (used for data entry), and enable a new
   pick routine.  The menu should change to contain Save (returns to standard
   form), Delete, and Add comment.  We should replace Save in the standard
   menu by Update.  Perhaps the standard menu could also contain Add comment?
*/

static void fixPath (LOOK look, BS bs) ;
extern void bsSubFuseType(OBJ obj, BS bs) ;

static void updatePick (int box)
{
  BS bs ;

  LOOKGET("updatePick") ;

  treeDump (0) ; /* to please te compiler */

  if (autoPreserve)
    displayPreserve () ;
  if (!box)
    { colourActive (look, FALSE) ;
      look->activebox = 0 ;
      return ;
    }
  else if (box >= arrayMax (look->content))
    { messerror("updatePick received a wrong box number %d",box);
      return ;
    }

  if (box == look->activebox)         /* a second hit - follow it */
    { bs = arr(look->content,look->activebox,BS) ;

      if (!bs)	      /* not a bs box, e.g. inside a textEntry */
	return ;

      if (bs->size & ON_FLAG) /* if hidden, show branch */
	{ bs->size ^= ON_FLAG ;
	  lookDraw(look) ;
	}
      else if (look->bsWait)
	messout ("Finish what you are doing first (or cancel)") ;
      else if (isDisplayBlocked())	/* from different update window */
	messout ("First deal with the other object you are waiting on") ;
      else if (bsIsTag(bs))
	{
	  if (bs->size & MODEL_FLAG)
	    fixPath (look, bs) ;
	  else
	    { 
	      bs->size |= ON_FLAG ;
	      lookDraw(look) ;
	    }
	}
      else if (!(bs->size & MODEL_FLAG)) /* in original object */
	{
	  look->bsWait = bs ;
	  look->handleWait = cacheStoreHandle (look->objet->x) ;
	  bs->size |= EDIT_NODE_FLAG ;	/* edit flag */
	  lookDraw (look) ;
	}
      else if (bsIsTag(bs->up) || !(bs->up->size & MODEL_FLAG))
	{
	  if (bs->size & SUBTYPE_FLAG)	/* a type */
	    {
	      bsSubFuseType (look->objet, bs) ;
	      lookDraw (look) ;
	    }
	  else
	    { 
	      look->bsWait = bs ;  
	      look->handleWait = cacheStoreHandle (look->objet->x) ;
	      if (class(bs->key))
		bs->size |= ADD_KEY_FLAG ;	/* displayBlock() + textEntry */
	      else
		bs->size |= ADD_DATA_FLAG ;	/* textEntry only */
	      lookDraw (look) ;
	    }
	}
      else
        messout ("Sorry - you must be next to a tag "
		 "or an existing part of the object") ;
    }
  else                              /* change highlighted entry */
    {
      colourActive (look, FALSE) ;
      look->activebox = box ;
      bs = arr(look->content,box,BS) ;
      if (!bs)		/* not a bs box */
	look->activebox = 0 ;
      colourActive (look, TRUE) ;
      if (bs)
	graphPostBuffer (drawBStext (bs)) ;
    }
}

static void updateKbd (int k)
{ char *vp ;
  BS bs ;
  int box ;
  LOOKGET("updateKbd") ;

  if (!look->activebox)
    return ;
  if (!(bs = arr(look->content,look->activebox,BS)))
    return ;

  switch (k)
    {
    case RETURN_KEY:
      updatePick (look->activebox) ;
      break ;
    case LEFT_KEY :
      while (bs->up->down == bs)
	bs = bs->up ;
      if (bs->up->up)
	bs = bs->up ;
      break ;
    case RIGHT_KEY :
      if (bs->right)
	bs = bs->right ;
      break ;
    case DOWN_KEY :
      if (bs->down)
	bs = bs->down ;
      else while (bs->up->down == bs)
	bs = bs->up ;
      break ;
    case UP_KEY :
      if (bs->up->down == bs)
	bs = bs->up ;
      else while (bs->down)
	bs = bs->down ;
      break ;
    case '\t': case SPACE_KEY :
      if (bs->right)
	bs = bs->right ;
      else if (bs->down)
	bs = bs->down ;
      else
	while (bs->up->up)
	  { while (bs->up->down == bs)
	      bs = bs->up ;
	    if (!bs->up->up)
	      break ;
	    bs = bs->up ;
	    if (bs->down)
	      { bs = bs->down ;
		break ;
	      }
	  }
      break ;
    default: return ;
    }

  if (assFind (look->bs2box, bs, &vp) &&
      (box = vp - (char *) 0) &&
      box != look->activebox)
    updatePick (box) ;
}

/********************************/

static void defuse (BS bs)	/* recursive */
{
  BS temp ;

  bs->size = (bs->size & (NEW_FLAG | ON_FLAG)) ;

  while (bs->right && (bs->right->size & MODEL_FLAG))
    { temp = bs->right ;
      bs->right = bs->right->down ;
      if (bs->right)
	bs->right->up = bs ;
      temp->up = 0 ;
      temp->down = 0 ;
      bsTreePrune (temp) ;
    }
  if (bs->right)
    defuse (bs->right) ;

  while (bs->down && (bs->down->size & MODEL_FLAG))
    { temp = bs->down ;
      bs->down = bs->down->down ;
      if (bs->down)
	bs->down->up = bs ;
      temp->up = 0 ;
      temp->down = 0 ;
      bsTreePrune (temp) ;
    }
  if (bs->down)
    defuse (bs->down) ;
}

/***********************************/

static void fixPath (LOOK look, BS bs)
{
  BS start = bs ;

  while (bs->size & MODEL_FLAG)	/* code for belonging to model */
    { bs->size &= ~MODEL_FLAG ;
      bs->size |= NEW_FLAG ;
      bs->timeStamp = sessionUserKey() ;
      while (bs->up && bs->up->down == bs)	/* go to top of column */
	bs = bs->up ;
      if (bs->up)		/* go to entry in previous column */
	bs = bs->up ;
    }

  defuse (look->objet->root) ;	/* leaves new stuff */

  bs = start ;
  for (bs = start ; bs->up ; bs = bs->up)
    { if (bs->bt && bs->bt->bsm && 
	  bs->bt->bsm->n.key & UNIQUE_BIT) /* UNIQUE */
	{ while (bs->down)
	    { look->objet->curr = bs->down ;
	      look->objet->modCurr = bs->bt->bsm ;
	      bsRemove (look->objet) ;
	    }
	  while (bs->up->down == bs && !bsIsComment(bs->up))
	    { look->objet->curr = bs->up ;
	      look->objet->modCurr = bs->bt->bsm ;
	      bsRemove (look->objet) ;
	    }
	  }
      while (bs->up->down == bs)	/* go to top of column */
	bs = bs->up ;
    }

  bsFuseModel (look->objet, 0, treeJustKnownTags) ;
  setOnFlagModelOnly (look->objet->root) ;
  look->objet->flag |= OBJFLAG_TOUCHED ;

  lookDraw (look) ;
}

/*********************************/

static void addComment (char* text)
{ 
  BS bs ;
  LOOKGET("addComment") ;

  if (!look->bsWait)
    messcrash ("screwup in addComment - no bs to add to") ;

  look->objet->curr = look->bsWait ;	/* point to add to */
  bsAddComment (look->objet, text, 'C') ;
  
  look->bsWait->size &= ~ADD_COMMENT_FLAG ;
  bs = look->bsWait ;
  look->bsWait = 0 ;
  fixPath (look, bs) ;
}

static void addKey (KEY key)
{ KEY old ;
  BS bs, bsm ;
  OBJ obj ;
  LOOKGET("addKey") ;

  if (!key)		/* cancellation */
    return ;

  if (!(bs = look->bsWait))
    messcrash ("screwup in addKey - bs is missing") ;

/* mieg, modified, sept 12.93 ***/

  if (bsIsComment (bs))		/* no matching bsm */
    {
      if (class(key) != class(bs->key))
	{
	  messout ("Sorry, that is not a comment you selected") ;
	  return ;
	}
      else
	bs->key = key ;
    }
  else
    { 
      old = bs->key ;
      bs->key = key ;		/* try to enter the key */
      if (!bs->bt || !(bsm = bs->bt->bsm))
	messcrash ("screwup in addKey - bsm is missing") ;
      if (!bsModelMatch (bs, &bsm))
	{ 
	  messout ("Sorry - that key is not in the right class") ;
	  bs->key = old ;
	  return ;
	}
    }
                          /* X-ref */
  if (pickXref(class(key)) || KEYKEY(bsm->n.key))	
    { 
      obj = look->objet ;
      if (!obj->xref) obj->xref = stackCreate (64) ;
      push (obj->xref, key, KEY) ;

      if (pickXref(class(key)))
        push (obj->xref, (_Quoted_in), KEY) ;
      else
	push (obj->xref, KEYKEY(bsm->n.key), KEY) ;
    }

  bs->size |= NEW_FLAG ;
  bs->timeStamp = sessionUserKey() ;
  bs->size &= ~(ADD_KEY_FLAG | EDIT_NODE_FLAG) ; /* must do this since display block is freed on return */
  ++look->activebox ;
  look->bsWait = 0 ;
  look->handleWait = 0 ;

  fixPath (look, bs) ;		/* clear "model" flags and redisplay */
}

static void addKeyText (char *text)
{ 
  BS bs, bs1 ;
  KEY key ;
  int table ;
  LOOKGET("addKeyText") ;
  
  if (!(bs = look->bsWait))
    messcrash ("Screwup in addKeyText - bs missing") ;

  table = class(bs->key) ;
  if (!*text || *text == '*' || !lexIsGoodName (text))
    {
      while(*text) 
	*text++ = 0 ;  /* necessary to set totally to 0 for graphTextEntry */
      graphRedraw() ;
      return ;
    }

  if (!lexword2key(text,&key,table))
    {
      if (table != _VText && 
	  !messQuery 
	  (messprintf("Unknown name - do you want to create a new %s called:%s", 
		      className(KEYMAKE(table,0)), text)))
	{ 
	  while(*text) 
	    *text++ = 0 ;  /* necessary to set totally to 0 for graphTextEntry */ 
	  graphRedraw() ;
	  return ;
	}
      else
	lexaddkey (text,&key,table) ;
    }
  bs1 = bs ;
  while (bs1->up && bs1->up->down == bs1)
    {
      bs1 = bs1->up ;
      if (bs1->key == key)
	{
	  strcpy (text, "sorry, this name is present above in the same column") ;
	  graphRedraw() ;
	  return ;
	}
    }
  display (key,0,0) ;		/* calls addKey via displayBlock() */
}

static void addData (char* text)
{
  BS bs ;
  int i ;
  float f ;
  mytime_t tm ;
  LOOKGET("addData") ;

  if (!(bs = look->bsWait))
    messcrash ("screwup in addData - no bs to add to") ;

  if (bs->key < _LastC)
    { bs->bt->cp = strnew (text, look->handleWait) ;
    }
  else if (bs->key == _Int)
    {
      if (sscanf (text,"%d",&i))
	bs->n.i = i ;
      else
	{ messout ("Sorry - not a good integer") ; return ; }
    }
  else if (bs->key == _Float)
    {
      if (sscanf (text,"%f",&f))
	bs->n.f = f ;
      else
	{ messout ("Sorry - not a good floating point number") ; return ; }
    }
  else if (bs->key == _DateType)
    {
      if ((tm = timeParse (text)))
	bs->n.time = tm ;
      else
	{ messout ("Sorry - not a good date") ; return ; }
    }
  else
    messcrash ("Bad data type %d in treedisp addData",bs->key) ;

  bs->size |= NEW_FLAG ;
  bs->timeStamp = sessionUserKey() ;
  bs->size &= ~(ADD_DATA_FLAG | EDIT_NODE_FLAG) ;	/* only get here if we added OK */
  look->bsWait = 0 ;
  look->handleWait =  0 ;
  ++look->activebox ;

  fixPath (look, bs) ;
}

/*******************/

static void editEntry (char *text)
{
  BS bs, bsm=0 ;
  KEY key ;
  OBJ obj ;
  int table ;
  LOOKGET("editEntry") ;

  if (!*text || !lexIsGoodName(text))
    { while(*text) 
	  *text++ = 0 ;  /* necessary to set totally to 0 for graphTextEntry */
      graphRedraw() ;
      return ;
    }

  if (!(bs = look->bsWait))
    messcrash ("screwup in editEntry - no bs to edit") ;
  if (bs->key < _LastN)
    addData (text) ;
  else
    { table = class(bs->key) ;
      if (!lexword2key(text,&key,table))
	{
	  if (!strcmp(text,"*") || !strcmp(text,"?"))
	    { messout("Please do not use ? or * as the name of an object, it would confuse subsequent searches") ;
	    return ;
	    }
	  else if (!pickXref(table) && !messQuery ("Unknown name - do you want to create a new object? "))
	    return ;
	  else
	    lexaddkey (text,&key,table) ;
	}

      if (pickXref(table))	/* especially comments! */
	{ obj = look->objet ;
	  if (!obj->xref) obj->xref = stackCreate (64) ;
	  push (obj->xref, bs->key, KEY) ;
	  push (obj->xref, (_Quoted_in | DELETE_BIT), KEY) ;
	}
      else 
	{ if (!bs->bt || !(bsm = bs->bt->bsm))
	    messcrash ("screwup in editEntry - bsm is missing") ;
	  if (KEYKEY(bsm->n.key))		/* must delete XREF */
	    { obj = look->objet ;
	      if (!obj->xref) obj->xref = stackCreate (64) ;
	      push (obj->xref, bs->key, KEY) ;
	      push (obj->xref, (KEYKEY(bsm->n.key) | DELETE_BIT), 
		    KEY) ;
	    }
	}

      addKey (key) ;
    }
}

/************************************************************/
/************************************************************/

static void chooseTagCancel(void)
{ LOOKGET("choosetagCancel") ;

  if (look->bsWait)
     return ;
  
  graphUnBlock(FALSE) ;
}

static void choosePick(int box)
{ BS bs ;
  LOOKGET("choosePick") ;
  
  if (look->bsWait)
     return ;

  if (!box)
    return ;
  else if (box >= arrayMax (look->content))
    { messerror("updatePick received a wrong box number %d",box);
      return ;
    }

  if (box == look->activebox)         /* a second hit - follow it */
    { bs = arr(look->content,look->activebox,BS) ;
      if (!bs)	      /* not a bs box, e.g. inside a textEntry */
	return ;
      if (bs == look->objet->root)
	return ;
      if (isDisplayBlocked())	/* from different update window */
	{ messout ("First deal with the other object you are waiting on") ;
	  return ;
	}
      if (bs->size & SUBTYPE_FLAG)        /* a type */
	{ bsSubFuseType (look->objet, bs) ;
	  lookDraw (look) ;
	}
      else 
	{ while (bs->size & MODEL_FLAG)	/* as in fixPath, but no UNIQUE problems */
	    { bs->size &= ~MODEL_FLAG ;
	      while (bs->up && bs->up->down == bs)
		bs = bs->up ;
	      if (bs->up)
		bs = bs->up ;
	    }
	  graphUnBlock(TRUE) ;
	}
     }
  else                              /* change highlighted entry */
    {
      colourActive (look, FALSE) ;
      look->activebox = box ;
      bs = arr(look->content,box,BS) ;
      if (!bs)		/* not a bs box */
	look->activebox = 0 ;
      colourActive (look, TRUE) ;
    }
}

/************/

static MENUOPT chooseMenu[]=
  {
   {chooseTagCancel , "Cancel"},
   {help, "Help"},
   {graphPrint, "Print"},
    {0, 0}
} ;

/************/

static LOOK chooseCreate (KEY key, OBJ objet, Graph oldGraph)
{ 
  LOOK look=(LOOK)messalloc(sizeof(struct LOOKSTUFF)) ;

   if (graphExists(oldGraph))
    { graphActivate(oldGraph) ;
      graphPop() ;
      lookDestroy () ;
      graphClear () ;
      graphGoto (0,0) ;
    }
  else 
    if (!displayCreate (TREE)) 
      return 0 ;


  graphSetBlockMode(GRAPH_BLOCKING) ;			    /* User must answer tree choose dialog. */
  graphRetitle (messprintf("Tag chooser: Class %s", className(key))) ;
  graphRegister (DESTROY, lookDestroy) ;
  graphRegister (PICK, choosePick) ;
  graphMenu (chooseMenu) ;
  graphHelp("Tag-chooser") ;

  look->magic = MAGIC;
  look->key = key;
  look->objet = objet ;
  look->tagChooserMode = TRUE ;

  look->graph = graphActive() ;
  graphAssociate (&MAGIC, look) ;

  look->content = arrayReCreate(look->content, 32,BS);
  look->bs2box = assReCreate (look->bs2box) ;
  if (look->textStack)
    stackDestroy (look->textStack) ;

  lookDraw (look) ;
  graphButtons (chooseMenu, 10, 1, 55) ;
  if (xmax < 65)
    graphTextBounds (65, ymax) ;

  return look ;
}

/************/

BOOL treeChooseTagFromModel(int *type, int *targetClass, int classe, 
			    KEY *tagp, Stack s, int continuation)
{ Graph oldGraph = graphActive() ;
  Graph myGraph = 0 ;
  Stack s1 ; Array toto ;
  BOOL lastTag = FALSE, contNeeded = FALSE ;
  LOOK look ;
  KEY key ; int i ;
  OBJ obj ;
  BS bs ; char *cp ;
  
  lexaddkey("Tag-chooser", &key, classe) ;
  obj = bsUpdate (key) ;

  while (bsGetKeyTags (obj, _bsRight, 0))
    bsRemove (obj) ;
  if (continuation && *tagp)
    bsAddTag (obj, *tagp) ;
  cacheMark (obj->x) ;
  bsFuseModel (obj, *tagp, treeJustKnownTags ? continuation : 0) ;
  
  look = chooseCreate(key, obj, myGraph) ;
  /* setOnFlag (look->objet->root) ;  */
  if (!look)
    { bsKill(obj) ;
      return FALSE ;
    }

  myGraph = look->graph ;

  messStatus("Please choose a Tag") ;

  if (!graphBlock())
    { look->objet = 0 ; /* no silly questions  in graphDestroy */
      look->tagChooserMode = FALSE ;
      if (graphActivate(myGraph))
	graphDestroy() ;
      graphActivate(oldGraph) ;
      graphPop() ;
      if (obj->magic)
	bsKill(obj) ;
      return FALSE ;
    }
  
  look->objet = 0 ; /* no silly questions in graphDestroy */
  look->tagChooserMode = FALSE ;
  if (graphActivate(myGraph))
    graphDestroy() ;
  graphActivate(oldGraph) ;
  graphPop() ;

  defuse(obj->root) ;
  bs= obj->root ;
  lastTag = FALSE ;
  
  s1 = stackCreate(32) ;
  while (bs->right)
    bs = bs->right ;
  if (!bs->up)
    messcrash("No bs->up in treeChooseTagFromModel") ;
  
  if (bsIsTag(bs))
    { *type = 'b' ;
      lastTag = TRUE ;
      *tagp = bs->key ;
      pushText (s1, name(bs->key)) ;
    }
  else if (bs->key == _Text)
    *type = 't' ;
  else if (bs->key == _Int)
    *type = 'i' ;
  else if (bs->key == _Float)
    *type = 'f' ;
  else if (bs->key == _DateType)
    *type = 'd' ;
  else if (class(bs->key))
    { *type = 'k' ;  
      *targetClass = class(bs->key) ;
    }

  while (bs = bs->up , bs->up)
    { if (bsIsTag(bs))
	{ if (!lastTag)
	    pushText (s1, name(bs->key)) ;
	  *tagp = bs->key ;
	  lastTag = TRUE ;
	}
      else 
	{ if (lastTag)
	    pushText (s1, "HERE  #") ;
	  else
	    pushText (s1, "HERE") ;	    
	  contNeeded = TRUE ; 
	  break ;
	}
    }

  stackCursor(s1, 0) ;
  toto = arrayCreate(32, char*) ;
  i = 0 ;
  while ((cp = stackNextText(s1)))
    array(toto, i++, char*) = cp ;
  
  if (i--)
    pushText(s, arr(toto,i,char*)) ;
  while(i--)
    catText(s, arr(toto,i,char*)) ;
  stackDestroy(s1) ;
  arrayDestroy(toto) ;

  bsKill(obj) ;

  if (contNeeded && continuation != 2)
    { messout (
	"Except from a derived column of the table maker, You must choose a rooted tag") ;
      return FALSE ;
    }
  return TRUE ;
}

/*********************************************************/

static void treeDispMailer (KEY key)
{
  OBJ  obj = bsCreate(key), obj1 ;
	
  if (!obj)
    return ;
  if (!bsGetData (obj, _E_mail, _Text, 0))
    { if (!messPrompt  ("Please specify an email address", "", "w"))
	goto fin ;
      if (!(obj1 = bsUpdate (key)))
	goto fin ;
	
      if (!bsAddData (obj1, _E_mail, _Text, freeword()))
	{ bsDestroy (obj1) ;
	  messout ("Sorry, I can't save your address") ;
	  goto fin ;
	}
      bsSave (obj1) ;
    }

  acedbMailer (key, 0, 0) ;

 fin:
  bsDestroy (obj) ;
}

/************************ end of file ********************/
