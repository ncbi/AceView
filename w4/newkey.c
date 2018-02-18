/*  File: newkey.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: user interface to edit operations 
 *              at level of whole objects:
 *			new object
 *			alias/fuse (separate these ?)
 *			rename
 *			delete
 *              included in graphical ace compilations
 * Exported functions:
 *              void addKeys (void)
 * HISTORY:
 * Last edited: Nov 24 14:45 1998 (fw)
 * * Dec  1 23:55 1996 (rd)
 * * Jun  6 18:36 1996 (rd)
 * * Feb  8 00:04 1992 (rd): added rename, delete and changed alias
 * * Jan 13 02:31 1992 (mieg): alias.
 * * Dec  2 19:02 1991 (mieg): graphClassEntry
 * Created: Mon Dec  2 19:01:40 1991 (mieg)
 *-------------------------------------------------------------------
 */
/************************************************************/

#include "acedb.h"
#include "graph.h"
#include "lex.h"
#include "display.h"
#include "pick.h"
#include "session.h"
#include "bs.h"
#include "systags.h"
#include "a.h"
#include "main.h"
#include "tree.h"

/************************************************************/

static Graph addGraph = 0 ;
static int currClassBox = 0 ;
static int currClass = 0 ;
static char nameText[128] ;
static Array box2class ;
static Array classList = 0 ;
static int classListBox = -1 ;
static char classNameBuf[32] = {'\0'} ;

/************************************************************/

static void findAndDisplayKey (void) ;
static void editNewKey (void) ;
static void (*preferredAction)(void) = editNewKey ;

/************************************************************/

static void doAction (void)
{ (*preferredAction)() ; }

static BOOL checkName (void)
{ char *cp ;

  if (!checkWriteAccess ())	/* may be lost after graph was created */
    { messout ("Sorry, you do not have Write Access") ; return FALSE ; }
  if (!currClass)
    { messout ("Please, first pick a class") ; return FALSE ; }
  if (!*nameText)
    { if (!messPrompt 
      ("Please enter the name of the object to be edited", "","w"))
	return FALSE ; 
      freenext () ;
      cp = freepos() ;
      strncpy (nameText, cp, 127) ;
    }

  if (!*nameText)
    { messout ("The name can not be empty") ; return FALSE ; }

  cp = nameText ; cp-- ;
  while (*cp++)
    if (*cp == '*' || *cp == '?')
      { messout("Please do not use * or ? in a new name") ;
	return FALSE ;
      }
  return TRUE ;
}

/******************* action functions ******************/

static KEY makeKey (void)
{ KEY key ;

  if (!checkName()) return 0 ;
  if (!lexaddkey (nameText, &key, currClass))
    messout ("%s %s already exists.", 		     
	     name(currClass), name(key)) ;
  return key ;
}

static void displayKey (KEY key)
{
  if (pickType(key) == 'B')
    { OBJ obj = bsUpdate(key) ;	/* to give it a cache entry */
      bsDestroy (obj) ;
      display (key, 0, TREE) ;
    }
  else
    display (key, 0, 0) ; /* A type cannot be shown as TREE */  
}

/****************/

static void findAndDisplayKey (void)
{
  KEY key ;

  preferredAction = findAndDisplayKey ;
  if (!checkName()) return ;

  if (!lexword2key (nameText, &key, currClass) &&
      !(messQuery
	(messprintf ("%s : %s does not yet exist\n"
		     "Do you want to Create it ?",
		    name (currClass), nameText)) &&
	(key = makeKey())))
    return ;

  displayKey (key) ;
}

static void makeAndDisplayKey (void)
{
  KEY key = makeKey () ;

  preferredAction = makeAndDisplayKey ;
  if (key)
    displayKey (key) ;
}

static void editNewKey (void)
{
  static Graph treeGraph = 0 ;
  KEY key ; 
  OBJ obj ;
  KEY class = currClass;
  unsigned char dummy;

  preferredAction = editNewKey ;
  if (!checkName()) return ;
  if (!lexword2key (nameText, &key, currClass))
    { if (!messQuery
	  (messprintf
	   ("%s : %s does not exist, do you want to create it ?", 
	    name (currClass), nameText)))
	return ;
      else
	lexaddkey (nameText, &key, currClass) ;
    }
  
  pickIsA(&class, &dummy);
  if (pickType(KEYMAKE(class, 0)) != 'B')
    { makeKey () ;
      messout ("Only Tree objects can be edited here, sorry") ; 
      return ;
    }
  
  obj = bsUpdate(key) ;
  if (graphActivate (treeGraph))
    treeDisplay (key, 0, TRUE) ;
  else
    { display (key, 0, TREE) ;
      treeGraph = graphActive() ;
    }
  bsDestroy(obj) ;  /* All this is to see empty objects */
  treeUpdate () ;
}

#ifdef JUNK
static void makeAlias (void)
{
  KEY key ;
  char *aliasName ;
  int classe = currClass ; 
 
  preferredAction = makeAlias ;
  if (!checkName())
    return ;
  if (!lexword2key (nameText, &key, classe)) 
    { messout ("Old name %s is unknown", nameText) ; return ; }
  if (!messPrompt ("Alias will result in both the old and new names being "
		    "recognized, but the name you type now will be the "
		    "one the database writes out (canonical name).  If it is "
		    "the name of an existing object then the two will be fused.  "
		    "Cancel and start again if you wish the name you "
		    "originally typed to be the canonical name."
		    "\n\nGive alias name:",
		    nameText, "t"))
    return ;

  aliasName = freeword () ;
  if (!*aliasName)
    { messout ("The name can not be empty") ; return ; }
  
  if (lexAlias (&key, aliasName, TRUE, TRUE))
    display (key, 0, TREE) ;
  else
    messout ("Alias failed, sorry") ;
}
#endif

static void changeName (void)
{
  KEY key ;
  char *newName ;
  int classe = currClass ; 
 
  preferredAction = changeName ;
  if (!checkName())
    return ;
  if (!lexword2key (nameText, &key, classe)) 
    { messout ("Old name %s is unknown", nameText) ; return ; }
  if (!messPrompt ("Rename is a hard rename function.  Unless you are just "
		    "changing case, the old name will no "
		    "longer apply to this or any other current object (you can "
		    "reuse it for a new object).  If the new name belongs to "
		    "an existing object the two objects will be fused.  "
		    "For both names to apply, cancel "
		    "now and use the Alias menu entry.  "
		    "\n\nGive new name:",
		    nameText, "t"))
    return ;

  newName = freeword () ;
  if (!*newName)
    { messout ("The name cannot be empty") ; return ; }
  if (!lexAlias (&key, newName, TRUE, FALSE))
    messout ("Rename failed, sorry") ;
}

/************************* delete functions **********************/

static void deleteKey (KEY key)
{ OBJ obj ;

  switch (pickType (key))
    {
    case 'A':
      arrayKill (key) ;		/* neat, heh! */
      break ;
    case 'B':
      if ((obj = bsUpdate (key)))
	bsKill (obj) ;
      else
	messout ("Object %s is locked", name (key)) ;
      break ;
    }
}

static void deleteObject (void)
{
  KEY key ;
  int classe = currClass ; 
 
  preferredAction = deleteObject ;
  if (!checkName()) return ;
   if (!lexword2key (nameText, &key, classe))
    messout ("Object %s is unknown", nameText) ;
  else if (messQuery (messprintf ("Delete all data associated %s ?",
				   name(key))))
    deleteKey (key) ;
}

/* Now part of the key set menu
static void deleteKeyset (void)
{
  KEYSET kset ;
  int i ;

  if (!keySetActive (&kset, 0))
    messout ("There is no active keyset") ;
  else if (messQuery (messprintf ("Delete all data from the %d objects "
				  "in the active keyset?", keySetMax (kset))))
    for (i = 0 ; i < keySetMax(kset) ; ++i)
      deleteKey (keySet(kset,i)) ;
}
*/

/****************************************/

static BOOL classListMenuCall (KEY k, int box)
{

	if (currClassBox)
		graphBoxDraw (currClassBox, BLACK, WHITE) ;

	currClassBox = box ;
	currClass = k ;
	strncpy(classNameBuf,name(currClass),31) ;
	classNameBuf[31] = '\0' ;
	graphBoxDraw(classListBox,BLACK,CYAN) ;

	return TRUE ;
}
static void classPick (void)
{ 
	KEY key ;
	if (graphSelect (&key, arrp (classList, 0, FREEOPT)))
		classListMenuCall (key, 0) ;
	graphBoxDraw(classListBox,BLACK,CYAN) ;
}

static void pick (int box, double x)
{ 
  if (!box)
    return ;

  if (box < arrayMax(box2class))
  {
      if (currClassBox)
        graphBoxDraw (currClassBox, BLACK, WHITE) ;

      currClassBox = box ;
      graphBoxDraw (currClassBox, WHITE, BLACK) ;
      currClass = arr(box2class,currClassBox,int) ;
      strncpy(classNameBuf,name(currClass),31) ;
      classNameBuf[31] = '\0' ;
      graphBoxDraw(classListBox,BLACK,CYAN) ;
  }
}

/******************* pick Object button *************/

static void blockPickObject (KEY key)
{
  int i ;
  KEY classKey ;
  unsigned char mask ;

  if (currClassBox)
      graphBoxDraw (currClassBox, BLACK, WHITE) ;
  currClassBox = 0 ;
  currClass = class(key) ;

  for (i = 0 ; i < arrayMax(box2class) ; ++i)
    { classKey = arr(box2class, i, int) ;
      pickIsA (&classKey, &mask) ;
      if (!mask && classKey == currClass)
      { 
          currClassBox = i ;
          graphBoxDraw (i, WHITE, BLACK) ;
      }
    }
  strncpy(classNameBuf,className(key),31) ;
  classNameBuf[31] = '\0' ;
  graphBoxDraw(classListBox,BLACK,CYAN) ;
  strncpy (nameText, name(key), 127) ;
  nameText[127] = 0 ;
  graphCompletionEntry (0, nameText, 0, 0, 0, 0) ;
}

static void pickObject (void)
{ displayBlock (blockPickObject, 0) ;
}

/****************************************************/

static void destroy (void)
{
  addGraph = 0 ;
  currClassBox = 0 ;
}

static MENUOPT newMenu[] = {
  {graphDestroy, "Quit"},
  {help, "Help"},
  {makeAndDisplayKey, "Create"},
  {findAndDisplayKey, "Display"}, 
  {editNewKey, "Edit"},
  /*   {makeAlias, "Alias-Fuse"},  suppressed, mieg feb 2002, acedb alias system is globally buged and utterly useless */
  {changeName, "Rename-Fuse"},
  {deleteObject, "Delete"},
/*   {deleteKeyset, "Delete KeySet"}, now in keyset menu */
   {0, 0}
}  ;

/*******************************************************/

static int newKeyClassComplete(char *cp, int len)
{
  if (currClass)
    return ksetClassComplete (cp, len , currClass) ;
  else
    return 0 ;
}

/******************* main function ****************/

void addKeys (void)
{
  int i, ncol = 3 ;
  int x, sep = 0 ;
  float y ;
  int nclasses = pickVocList[0].key ;
  int classPickBox ;
  int graphWidth, graphHeight ;

  if(!checkWriteAccess())
    { messout("Sorry, you do not have Write Access");
      return ;
    }

  if (graphActivate (addGraph))
    { graphPop () ;
      return ;
    }
  
  addGraph = displayCreate(DtAlias) ;
  graphRegister (DESTROY, destroy) ;
  graphRegister (PICK, pick) ;
  graphMenu (newMenu) ;
  
  box2class = arrayReCreate (box2class, 32, int) ;

  for (i = 1 ; i <= nclasses ; ++i)
    if (strlen(pickVocList[i].text) > sep)
      sep = strlen (pickVocList[i].text) ;
  sep += 2 ;
  
  classListCreate (&classList) ;
  classPickBox = graphButton ("Class:", classPick, 16, 2.2) ;
  graphBoxFreeMenu(classPickBox, (FreeMenuFunction) classListMenuCall, 
		   arrp (classList, 0, FREEOPT) ) ;
  classListBox = graphBoxStart() ;
  graphTextPtr(classNameBuf, 24, 2.35, 31) ;
  graphBoxEnd() ;
  graphBoxDraw(classListBox,BLACK,CYAN) ;

  graphButton ("Pick object", pickObject, 2.0, 2.2) ;

  graphFitBounds (&graphWidth, &graphHeight) ;
  graphButtons (newMenu, 2, 4.5, graphWidth) ;
  x = 0 ; y = 6.5 + (graphWidth < 60 ? 2 : 0);
  for (i = 1 ; i <= nclasses ; ++i) 
    { KEY class1 = pickVocList[i].key ;
      unsigned char dummy ;

      pickIsA(&class1, &dummy) ;
      
      if (pickList[class1 & 255].type != 'B')
	continue ;
      if (pickList[class1].protected)
	continue ;
      array (box2class,graphBoxStart(),int) = pickVocList[i].key ;
      graphText (name(pickVocList[i].key), 2.0 + (float)(x*sep), y) ;
      graphBoxEnd () ;
      if (!(++x%ncol))
	{ x = 0 ; y += 1 ; }
    }

  graphText ("Name: ", 2.0, 0.5) ;
  *nameText = 0 ;
  graphCompletionEntry (newKeyClassComplete, 
			nameText, 127, 9.0, 0.5, 
			TEXT_ENTRY_FUNC_CAST doAction) ;

  currClass = currClassBox = 0 ; 
  graphFitBounds (&graphWidth, &graphHeight) ;
  graphButtons (newMenu, 2, y + 2, graphWidth) ;
  graphRedraw () ;
} /* addKeys */

/************************** eof *****************************/
 
