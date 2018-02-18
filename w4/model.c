/*  File: model.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **  Reads the models of the ACeDB program.                   
 * Exported functions:
     readModels
 * HISTORY:
 * Last edited: Mar  4 16:24 1999 (fw)
 * * Oct 21 15:01 1998 (edgrif): Replace messout + exit() with messExit.
 * Created: Wed Nov  6 13:27:39 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: model.c,v 1.9 2015/08/18 14:33:42 mieg Exp $ */

#define DEFINE_OBJ
typedef struct sobj *OBJ ;
#include "acedb.h"

#include "bs_.h"
#include "cache.h"		/* used just in duplicate model  */
#include "parse.h"
#include "model.h"

/**************************************************/

static KEYSET checkModelSyntax (void) ;
#define MTEST FALSE
static VoidRoutine modelChangeRoutine = 0;
 
/**************************************************/

VoidRoutine modelChangeRegister (VoidRoutine func)
{ 
  VoidRoutine old = modelChangeRoutine ; 
  modelChangeRoutine = func ; 
  return old ;
} /* modelChangeRegister */

/**************************************************/

BOOL isModel (KEY key)		/* old models */
{ 
#if !defined(NEW_MODELS)
 return (KEYKEY(key) == 0) ;
#else
 return (class(key) == _VModel) ;
#endif
} /* isModel */

/**********************************************************************/

static int line ;
static Stack s1 = 0 ;

static void cleanVisibility (void)
{ KEY dummy ;
  OBJ obj ;
  int i ;
  KEYSET ks ;

  if (!lexword2key("Class", &dummy, _VClass))
    return ;
  ks = query (0, "Find Class Visibility") ;
 
  for (i = 0 ; i < keySetMax (ks) ; i ++)
    { if ((obj = bsUpdate (keySet(ks,i))))
      { if (bsFindTag (obj, _Visibility))
	  bsRemove (obj) ;
        bsSave (obj) ;
	}
    }
  i = 256 ;
  while (i--) 
    pickList[i].visible = 0 ;

  return;
} /* cleanVisibility */


static void cleanOldModels (void)
{
  KEY key = 0 ;
  OBJ obj ;
  extern BOOL READING_MODELS ;   /* from bssubs.c */

  READING_MODELS = TRUE ;

  while (lexNext (_VModel, &key))
    { if (iskey(key) == 2 && (obj = bsUpdate (key)))
	bsKill (obj) ;
      lexkill(key) ;
    }  
  READING_MODELS = FALSE ;

  return;
} /* cleanOldModels */

static int getLine (Array aNode, int *plevel)	/* returns FALSE at eof */
{
  static FILE *fil ;
  static char *sysModelFileName ;
  static Array levelCount ;
  static Stack sword ;		/* holds arbitrarily long string */
  char c = 0, *cp, *word, buf[2], inQuotes = 0, old = 0 ;
  KEY dummy ;
  static int pass = 1 ;
  int iNode, counter, inword, eow, eol, level=0 ;
#ifndef NEW_MODELS
  int table ;
#endif 
  Stack s ;

  buf[1] = 0 ;
pass2:
      iNode = counter = 0 ;
      inword =  eow = FALSE ;

  if (!fil)
    { if (pass == 1)
	{ if (!(fil =  filtmpopen (&sysModelFileName,"w")))
	    messExit ("Can't open sysmodel file to write : %s (%s)",
		      sysModelFileName, messSysErrorText()) ;
	  s = sysModels () ;
	  stackCursor (s, 0) ;
	  while ((cp = stackNextText(s)))
	    fprintf(fil, "%s\n", cp) ;
	  stackDestroy (s) ;
	  filclose (fil) ;
	  if (!(fil =  filopen (sysModelFileName,0, "r")))
	    messExit ("Can't reopen /tmp/sysmodels.wrm file (%s)",
		      messSysErrorText()) ;
	}
      else
	{
	  cp = sessionFilName ("wspec/models","wrm","r") ;
	  if (!cp || !(fil =  filopen (cp, 0,"r")))
	    messExit ("Can't open model file : wspec/models.wrm (%s)",
		      messSysErrorText()) ;
	}
      sword = stackCreate (64) ;
      levelCount = arrayCreate (8,int) ;
      array(levelCount,0,int) = 0 ;
      line = 1 ;
    }

  arrayMax(aNode) = 0 ; 
  for (eol = FALSE ; !eol && !feof(fil) ; ++counter)
    { old = c ;
      switch (c = (getc (fil) | inQuotes))
	{
	case '\t': 
	  while (++counter%8) ;
	  --counter ; eow = TRUE ; break ;
        case '/':
          while((c = getc(fil)) && c != '\n' && !feof(fil) ) ;
	  /* empty line and fall through */
	case '\n': case (char)EOF: /* case ('\n' | 128): */
	  eol = TRUE ; eow = TRUE ; inQuotes = FALSE ; break ;
	case ' ':
	  eow = TRUE ; break ;
        default: /* this includes the inQuotes modified characters */
	  c &= 127 ;
	  if (old != '\\' && c == '\"') 
	    inQuotes = inQuotes ? 0 : 128 ; /* toggle */
	  if (!inword)
	    { if (!iNode)	/* find level */
		{ for (level = 0 ; level < arrayMax(levelCount) ; level++)
		    if (arr(levelCount,level,int) >= counter)
		      break ;
		  if (level == arrayMax(levelCount) ||
		      arr(levelCount,level,int) > counter)
		    messcrash ("Model file: tree screwup line %d level %d",
			       line, level) ;
		  *plevel = level ;
		}
	      array(levelCount,level,int) = counter ;
	      ++level ;
	      stackClear (sword) ;
	      inword = TRUE ; eow = FALSE ;
	    }
	  old = buf[0] = c ;  /* a trick because push char fails on alpha */
	  catText (sword,buf) ;
	}
      if (inword && eow)			/* end of word */
	{ 
	  word = stackText(sword,0) ;
	  
	  switch (*word)
	    {
	    case '?':
	    case '#':
	    case '@':
	      if (!strlen(word+1))
		messcrash ("Model file: no space allowed behind a %c line %d level %d",
			   *word, line, level) ;
	      break ;
	    default:
	      break ;
	    }
	  if ( *word == '@' ) /* a query */
	    { word++ ;
	      if (*word == '\"')
		{ 
		  char *cp ;
		  word++ ;
		  if (strlen(word))
		    { 
		      cp = word + strlen(word) - 1 ;
		      if (*cp == '\"')
			*cp = 0 ;
		    }
		}	
	      if (lexMax (_VQuery) < 2)  /* a hack: KEYKEY = 1 breaks treedisp */
		lexaddkey("(NULL Query)", &dummy, _VQuery) ;
	      lexaddkey(word, arrayp(aNode,iNode,KEY), _VQuery) ;
	    }
#ifndef NEW_MODELS
	  else if ((*word == '?' || *word == '#') &&  /* Class pointer or type */
		   (table = pickMakeClass(word+1)))
	    /* Note that ?System is forbidden 
	     * This is because implicit creation of a tag would create havock 
	     */

	    { array(aNode,iNode,KEY) = KEYMAKE(table,*word == '#' ? 1 : 0) ;
	      if(!lexMax(table))  /* New class being defined */
		lexaddkey(0,&dummy,table) ;
	    }
#else
	  else if (*word == '?')   /* Class pointer */
	    {
              pickMakeClass(word+1) ;  /* declare the class */
	      if (iNode)
		lexaddkey (word+1, &dummy, _VClass) ;
	      else  /* Declaration of a new class */
		{
		  *word = '#' ;
		  lexaddkey (word, &dummy, _VModel) ; /* interpret as a Model */
		}
	      array(aNode,iNode,KEY) = dummy ;
	    }
	  else if (*word == '#')   /* Constructed type */
	    {
	      lexaddkey (word, &dummy, _VModel) ;
              array(aNode,iNode,KEY) = dummy ;
	    }
#endif
	  else
	    lexaddkey (word,arrayp(aNode,iNode,KEY),0) ;

	  ++iNode ;
	  inword = FALSE ;
	}
    }

#ifdef JUNK	
  if (MTEST)      /* keep stderr here since it is for debugging */
  { fprintf (stderr,"%d : ",*plevel) ;
    for (iNode = 0 ; iNode < arrayMax(aNode) ; ++iNode)
      fprintf (stderr,"%s ",name(arr(aNode,iNode,KEY))) ;
    fprintf (stderr,"\n") ;
  }
#endif

  if (feof (fil))
    { filclose (fil) ;
      fil = 0 ;
      arrayDestroy (levelCount) ;
      stackDestroy (sword) ;
      
      if (pass++ == 1)
      {
	filtmpremove (sysModelFileName) ;
	goto pass2 ;
      }
      pass = 1 ;
      return FALSE ;
    }
  filremove ("/tmp/sysmodels","wrm") ;
  ++line ;

  return TRUE ;
} /* getLine */

/*******************************************************/

static void getSubClasses(void)
{
  char *cp = sessionFilName ("wspec/subclasses","wrm","r") ;
  FILE *f = cp ? filopen(cp, 0, "r") : 0 ;
  KEY dummy ;

  if (lexword2key ("Class", &dummy, _VClass))
    {			/* fails during initial database construction */
      KEYSET old = 
	query (0, "FIND Class Is_a_subclass_of OR Filter OR Mask") ;
      OBJ obj ;
      int i  = keySetMax (old) ;

      while (i--)
	if ((obj = bsUpdate(keySet(old,i))))
	  {
	    if (bsFindTag (obj, _Is_a_subclass_of))
	      bsRemove(obj) ;
	    if (bsFindTag (obj, _Filter))
	      bsRemove(obj) ;
	    if (bsFindTag (obj, _Mask))
	      bsRemove(obj) ;
	    bsSave (obj) ;
	  }
    }
  if (f)
    {
      overRideParseProtection = TRUE ;
      parseFile (f, 0, 0) ;  /* will close f */
      overRideParseProtection = FALSE ;
    }

  return;
} /* getSubClasses */

/*******************************************************/

static void removeConstraints (void)
{
  OBJ obj ;
  KEY key = 0 ;
  
  while (lexNext (_VClass, &key))
    if ((obj = bsUpdate (key)))
      { if (bsFindTag (obj, _Constraints))
	  bsRemove(obj) ;
	bsSave (obj) ;
      }
  
  queryClearConstraints () ;

  return;
} /* removeConstraints */

/*********/

static void getConstraints (void)
{
  char *cp = sessionFilName ("wspec/constraints","wrm","r") ;
  FILE *f = cp ? filopen(cp, 0, "r") : 0 ;
  
  removeConstraints () ;

  if (f)
    {
      overRideParseProtection = TRUE ;
      parseFile (f, 0, 0) ;  /* will close f */
      overRideParseProtection = FALSE ;
      pickRegisterConstraints () ;
    }
} /* getConstraints */

/******************************************************************/
#ifndef NEW_MODELS

static void cleanUpModels (BS bs)
{
  if (bs->right)
    cleanUpModels (bs->right) ;
  if (bs->down)
    cleanUpModels (bs->down) ;
  if (class(bs->key) > 1) /* not on tags */
    { char buf[256] ;
      switch (KEYKEY (bs->key))
	{
	case 0:
	  sprintf(buf, "?%s",className(bs->key)) ;
	  lexaddkey (buf, &bs->key, _VModel) ;
	  break ;
	case 1:   
	  sprintf(buf, "#%s",className(bs->key)) ;
	  lexaddkey (buf, &bs->key, _VModel) ;
	  break ;
	}
    }

  return;
} /* cleanUpModels */

/******************************************************************/

static void aliasModels (void)
{
  KEY k = 0, k1 ;
  char *cp, buf [256] ;
  
  while (lexNext(_VModel, &k))
    { cp = name (k) ;
      if (*cp == '#')
	{ sprintf(buf, "?%s",cp + 1) ;
	  if (lexword2key (buf, &k1, _VModel))
	      lexAlias (&k, buf, FALSE, TRUE) ;	}
    }

  return;
} /* aliasModels */

/******************************************************************/

static void duplicateModel(OBJ Model)
{
  KEY model ;
  OBJ obj ;

  if(class(Model->key) == _VModel)  /* prevent loops */
    return ; 
  lexaddkey(name(Model->key), &model, _VModel) ;
  obj = bsUpdate(model) ;
  bsTreePrune(obj->root) ;
  obj->root = bsTreeCopy(Model->root, cacheStoreHandle (obj->x)) ;
  obj->root->key = model ;
  cacheMark(obj->x) ;
  cleanUpModels (obj->root->right) ;
  obj->flag |= OBJFLAG_TOUCHED ;
  bsSave(obj) ;

  return;
} /* duplicateModel */


/********************************************************************/
/* Verify that the obj0 model is the same as the Model class model */

#ifdef NOT_USED_ANYWHERE_BUT_LEAVE_FOR_NOW__RD_980821

static BOOL treeCompare(BS one, BS two)
{
  int  fault = 0;
  static char buf[256]; /* recursive routine, save memory */

  if (class(one->key) > 1 && (KEYKEY(one->key) == 0 || KEYKEY(one->key) == 1))
    {  sprintf(buf, "?%s",className(one->key)) ;
       if (strcmp(buf, name(two->key)) != 0)
	fault++;
    }
  else
    { if (one->key != two->key)
	fault++;
    }
  
  if (one->down)
    { if (two->down)
	{ if (!treeCompare(one->down, two->down))
	    fault++;
	}
    else
      fault++;
    }
  else
    { if (two->down)
	fault++;
    }
  
  if (one->right)
    { if (two->right)
	{ if (!treeCompare(one->right, two->right))
	    fault++;
	}
      else
	fault++;
    }
  else
    { if (two->right)
	fault++;
    }

  return (fault == 0); 
} /* treeCompare */



static BOOL verifyModel(int class)
{ KEY modelKey, zeroKey;
  OBJ zero, model;
  BOOL ok = TRUE;

  if ((pickList[class].type == 'B') && pickList[class].model)
    { zeroKey = KEYMAKE(class, 0);
      modelKey = pickList[class].model;
      
      zero = bsCreate(zeroKey);
      model = bsCreate(modelKey);
      
      if (zero && zero->root->right &&
	  model && model->root->right && 
	  !treeCompare(zero->root->right, model->root->right))
	ok = FALSE;
      
      bsDestroy(zero);
      bsDestroy(model);
    }

  return ok;
} /* verifyModel */

void verifyModels(void)
{
  int i;

  for (i = 1 ; i < 256 ; ++i)
    if (!verifyModel(i))
      messout("Failed to verify class %s.\n", className(KEYMAKE(i,0)));

  return;
} /* verifyModels */
#endif /* NOT_USED */


#endif   /* !NEW_MODELS */

/******************************************************************/
/* updates the models 
 * called by update.c and readModels() 
 */

BOOL getNewModels (void)   
{
  Array aNode,aDown ;
  KEY   key = 0 ;
#ifdef NEW_MODELS
  KEY  oldKey = 0 ;
#endif
  OBJ   obj = 0 ;
  BS    bs=0 ;
  int   i,level=0 ;
  int   table = 0 ;
  VoidRoutine previousCrashRoutine ;
  extern VoidRoutine messcrashroutine ;
  extern void simpleCrashRoutine (void);
  extern BOOL READING_MODELS;	/* from bssubs.c */

  if (sizeof(void *) != STACK_ALIGNMENT)
    messcrash ("This is a %d bits machine, but the program was compiled for a %d bits machine, please select the correct ACEDB_MACHINE and recompile", 8*sizeof(void *), 8*STACK_ALIGNMENT) ;

  if (!isWriteAccess())
    { messout ("Sorry, you do not have Write Access");
      return FALSE ;
    }
  pickInit () ;

  previousCrashRoutine = messcrashroutine ;
  messcrashroutine = simpleCrashRoutine; /* don't allow crash recovery */

  cleanOldModels () ;
  cleanVisibility () ;
  sysClassOptions () ;
 
  aNode = arrayCreate (4,KEY) ;
  aDown = arrayCreate (8,BS) ;	/* available down pointers */
  s1 = stackCreate(500) ;	/* to store queries and so on */
  pushText(s1, "dummy") ;	/* to move stackMark(s1) away from 0 and 1  */
  READING_MODELS = TRUE ;

  while (getLine (aNode,&level))
    if (arrayMax(aNode))
      { i = 0 ;
	key = arr(aNode,i,KEY) ;
	if
#ifndef NEW_MODELS
	  (!KEYKEY(key) || KEYKEY(key) == 1) /* new class */
#else
	  /* (class(key) == _VClass || */ ( class(key) == _VModel)
#endif
	  { i++ ;
	    if (level)
	      messcrash ("%s starts line but not new class at line %d",
			 name(key), line) ;
	    if (obj)	/* close it */
	      { 
#ifdef NEW_MODELS
		oldKey = bsKey (obj) ;
#endif
		if (MTEST)
		  bsTreeDump (obj->root) ;
		cacheMark (obj->x) ;
#ifndef NEW_MODELS
		duplicateModel(obj) ;
#endif
		obj->flag |= OBJFLAG_TOUCHED ;
		bsSave (obj) ;
#ifndef NEW_MODELS
		bsMakePaths (table) ;
#else
		bsMakePaths (oldKey) ;
#endif
	      }
#ifndef NEW_MODELS
	    table = class(key) ; key = KEYMAKE(table,0) ;
#else
	    table = KEYKEY(key)  ;
#endif
	    obj = bsUpdate (key) ;
	    if (!obj)
	      messExit ("Can not form model for %s at line %d",
			name(key), line) ;
	    bs = obj->root ;
	    bs->right = BSalloc() ;  /* very first addition is right */
	    bs->right->up = bs ;
	    bs = bs->right ;
	    bs->key = arr (aNode,i++,KEY) ;
	    ++level ; /* hocus pocus */
	    array(aDown,level++,BS) = bs ;
	  }
	else if (level <= 0 || level > arrayMax(aDown))
	  messExit ("Model screwup: impossible level line %d",
		    line) ;
	else
	  { bs = arr(aDown,level,BS) ;
	    bs->down = BSalloc() ;     /* first addition is down */
	    bs->down->up = bs ;
	    bs = bs->down ;
	    bs->key = arr (aNode,i++,KEY) ;
	    array(aDown,level++,BS) = bs ;
	  }
	    
        while (i < arrayMax(aNode))
	  { bs->right = BSalloc () ;	/* rest are right */
	    bs->right->up = bs ;
	    bs = bs->right ;
	    bs->key = arr (aNode,i++,KEY) ;
	    array(aDown,level++,BS) = bs ;
	  }
	arrayMax(aDown) = level ;
      }
  
  if (obj)
    { if (MTEST)
	bsTreeDump (obj->root) ;
      cacheMark (obj->x) ;
#ifndef NEW_MODELS
      duplicateModel (obj) ;
#endif
      obj->flag |= OBJFLAG_TOUCHED ;
#ifdef NEW_MODELS
      oldKey = bsKey (obj) ;
#endif
      bsSave (obj) ;
#ifndef NEW_MODELS
		bsMakePaths (table) ;
#else
		bsMakePaths (oldKey) ;
#endif
    }
	
  getSubClasses() ;
  getConstraints () ;
  READING_MODELS = FALSE ;
  arrayDestroy(aDown) ;
  arrayDestroy(aNode) ;
  stackDestroy(s1) ;
#ifndef NEW_MODELS
  aliasModels () ;
#endif
  pickCheckConsistency() ;


  lexAlphaClassList (TRUE) ; /* re alpha order class names */
  /* used to remake class lists etc. when models change */
  if (modelChangeRoutine)
    (*modelChangeRoutine)();

  {
    KEYSET ks = checkModelSyntax () ;
    messdump ("Read models got %d models\n", keySetMax(ks)) ;
    keySetDestroy (ks) ;
  }

  bIndexNewModels() ;


  messcrashroutine = previousCrashRoutine;

  return TRUE ;
} /* getNewModels */

/*******************************************************/
/***** Authorisation checks ****************************/

BOOL readModels (void)
{
  /* Only the acedb administrator is entitled to proceed
   * I check this by testing write permission on database/ACEDB.wrm
   */

  /* First you must own the running process */

  if (getuid() != geteuid() || 
      !sessionFilName ("database/ACEDB","wrm","w"))
    { messout ("Sorry, only the administrator owning "
	       "the $(ACEDB) directory can modify the models");
      return FALSE ;
    }

  if (thisSession.session != 1)
    /* if we're re-reading the models not during the bootstrapping
       of an empty database we may need to ask a few more questions */
    {
      if (isWriteAccess())
	{
	  if (messQuery ("You did not save the current session: "
			 "should I ?"))
	    sessionClose (TRUE) ;
	  else
	    return FALSE ;		/* don't want to save, can't proceed */
	}

      if (!messQuery ("Watch out ! An erroneous modification of "
		      "the models may screw up the system.\n"
		      "Do you want to proceed?"))
	return FALSE ;

      /* we either got here without write access or 
	 sessionClose dropped it ... */

      if (!sessionGainWriteAccess ())
	/* ... we are therefore guaranteed that the write access 
	   will change. That will re-draw and change buttons/menus
	   if necessary */
	{ 
	  messout("Sorry, you cannot gain write access.  ",
		  "I cannot update the models.") ;
	  return FALSE ;
	}
    }

  getNewModels() ;
  
#ifdef JUNK
  if (thisSession.session != 1) /* close(1) would crash inside session 1 */
    sessionClose (TRUE) ;   /* save */
#endif

  return TRUE ;			/* else it will have crashed */
} /* readModels */

/**********************************************************************/
/**************** syntax check - everything now done here *************/

typedef struct
{ KEY fromModel ;
  KEY fromTag ;
  KEY toModel ;
  KEY toTag ;
} XrefInfo ;

static Array xinfo ;		/* all XREF info */

/* #define DEBUG 1 */
static int depth = 1 ;

static BOOL atRoot ;

static void checkRecurse (KEY model, BS bs, BOOL inTags,
			  Associator tag2right)
{
  char *zero = 0 ;

#ifdef DEBUG
  { int i ;
    for (i = 2*depth ; i-- ;)
      putchar (' ') ;
    printf ("%s\n", name(bs->key)) ;
  }
#endif

  if (atRoot)
    atRoot = FALSE ;
  else if (bs->key == _UNIQUE)
    { if (bs->down)
	messerror ("Model error: DOWN from UNIQUE in model %s", name(model)) ;
      if (!bs->right)
	messerror ("Model error: UNIQUE in model %s with nothing to the right",
		   name(model)) ;
    }
  else if (bs->key == _REPEAT)
    { if (bs->right)
	messerror ("Model error:  something right of REPEAT in model %s",
		   name(model)) ;
      if (bs->down)
	messerror ("Model error: something down from REPEAT in model %s",
		   name(model)) ;
      return ;
    }
  else if (bs->key == _COORD)
    { if (bs->down)
	messerror ("Model error: something down from COORD in model %s",
		   name(model)) ;
      if (bs->up->key != _Int && bs->up->key != _Float)
	messerror ("Model error: COORD not preceeded by Int or Float in model %s",
		   name(model)) ;
      return ;
    }
  else if (class(bs->key) == _VQuery && !atRoot)
    { if (bs->right)
	messerror ("Model error: something right of attach query %s in model %s",
		   name(bs->key), name(model)) ;
      if (bs->down)
	messerror ("Model error: something down from attach query %s in model %s",
		   name(bs->key), name(model)) ;
      return ;
    }
  else if (bs->key == _ANY)
    { if (bs->down)
	messerror ("Model error: ANY with something DOWN in model %s") ;
    }
  else if (class(bs->key) == _VModel) /* check subtypes */
    { KEY key ;

      lexword2key (name(bs->key), &key, _VModel) ;
      if (key != bs->key)	/* then bs->key is #..., a type */
	{ OBJ obj ;

	  if (bs->right)
	    messerror ("Model error: something right of subtype %s in model %s",
		     name(bs->key), name(model)) ;
	  if (bs->down)
	    messerror ("Model error: something down from subtype %s in model %s",
		     name(bs->key), name(model)) ;
	  if (!(obj = bsCreate (bs->key)))
	    messerror ("Model error: no model for subtype %s in model %s",
		     name(bs->key), name(model)) ;
	  else
	    { if (!obj->root->right)
		messerror ("Model error: empty model for subtype %s in model %s",
			 name(bs->key), name(model)) ;
	      bsDestroy (obj) ;
	    }
	}
      inTags = FALSE ;
    }
  else if (bs->key < _LastN || class (bs->key))
    inTags = FALSE ; 
  else				/* a tag */
    { if (!inTags)
	messerror ("Model error: tag %s outside rooted subtree in model %s",
		   name(bs->key), name(model)) ;
      if (assFind (tag2right, zero + bs->key, 0))
	messerror ("Model error: duplicate tag %s in model %s",
		   name(bs->key), name(model)) ;
      else if (!bs->right || (bs->right->key == _UNIQUE && !bs->right->right))
	assInsert (tag2right, zero + bs->key, zero) ;
      else if (bs->right->key != _UNIQUE)
	assInsert (tag2right, zero + bs->key, zero + bs->right->key) ;
      else
	assInsert (tag2right, zero + bs->key, zero + bs->right->right->key) ;
    }

  if (!bs->right)		/* done */
    return ;

  if (bs->right->key == _XREF)
    { 
      if (class(bs->key) != _VModel && class(bs->key) != _VClass)
	messerror ("Model error: XREF from %s not a class in model %s",
		   name(bs->key), name(model)) ;
      if (bs->right->down)
	messerror ("Model error: DOWN from XREF %s in model %s",
		   name(bs->key), name(model)) ;
      if (!bs->right->right || class(bs->right->right->key) ||
	       bs->right->right->key < _LastN)
	{ messerror ("Model error: no tag after XREF from %s in model %s",
		     name(bs->key), name(model)) ;
	  return ;
	}
      if (bs->right->right->down)
	messerror ("Model error: down from tag after XREF from %s in model %s",
		   name(bs->key), name(model)) ;
      { XrefInfo *inf = arrayp(xinfo, arrayMax(xinfo), XrefInfo) ;
	inf->fromModel = model ;
	if (bs->up->key == _UNIQUE)
	  inf->fromTag = bs->up->up->key ;
	else
	  inf->fromTag = bs->up->key ;
#ifndef NEW_MODELS
	inf->toModel = bs->key ;
#else
	if (class(bs->key) == _VModel)
	  inf->toModel = bs->key ;
	else
	  lexword2key(messprintf("#%s",name(bs->key)), &(inf->toModel), _VModel) ;
	if (!inf->toModel) invokeDebugger() ;
#endif
	inf->toTag = bs->right->right->key ;
      }
      bs = bs->right->right ;		/* eat up XREF and tag */
      if (!bs->right)
	return ;
    }

  bs = bs->right ;
  ++depth ;
  if (bs->down) do
    { if (bs->key < _LastN || class (bs->key))
	messerror ("Model error: non-tag %s in column in model %s",
		   name(bs->key), name(model)) ;
      checkRecurse (model, bs, inTags, tag2right) ;
    } while ((bs = bs->down)) ;
  else
    checkRecurse (model, bs, inTags, tag2right) ;
  --depth ;

  return;
} /* checkRecurse */

static KEYSET checkModelSyntax (void)
{ 
  int i ;
  KEY model ;
  OBJ obj ;
  Associator tag2right ;
  Array tagsInModel ;	/* of Arrays */
  KEYSET ks = keySetCreate() ;
  char *zero = 0 , *cp ;
  int errorCount ;

  tagsInModel = arrayCreate (256, Array) ;
  xinfo = arrayCreate (256, XrefInfo) ;

  errorCount = messErrorCount() ;

  for (i = 1 ; i < lexMax(_VModel) ; ++i)
    { model = KEYMAKE(_VModel, i) ;
      if (iskey(model) != 2)
	continue ;
      keySet(ks, keySetMax(ks)) = model ;
      obj = bsCreate (model) ;
      if (obj->root->right)
	{ tag2right = array (tagsInModel, i, Associator) = assCreate () ;
#ifdef DEBUG
	  printf ("%s\n", name(model)) ;
#endif
	  atRoot = TRUE ;
	  checkRecurse (model, obj->root, TRUE, tag2right) ;
	}
      bsDestroy (obj) ; /* mieg */
    }

  for (i = 0 ; i < arrayMax (xinfo) ; ++i)
    { XrefInfo *inf = arrp(xinfo, i, XrefInfo) ;
      KEY toModel ;
#ifdef DEBUG
      printf ("%d: in model %s, tag %s XREF's to tag %s in model %s\n", 
	      i, name(inf->fromModel), name(inf->fromTag), name(inf->toTag), name(inf->toModel)) ;
#endif
      toModel = inf->toModel ;
      if ((tag2right = array (tagsInModel, KEYKEY(toModel), Associator)))
	{ if (!assFind (tag2right, zero + inf->toTag, &cp))
	    messerror ("Model error: In model %s, tag %s XREF's to tag %s in class %s\n"
		       "             can't find tag %s in model %s",
		       name(inf->fromModel), name(inf->fromTag), name(inf->toTag), name(inf->toModel),
		       name(inf->toTag), name(toModel)) ;
	  else 
	    { 
/*           there is a problem when using this paragraph because of XREF embedded in contructed types
              KEY backClass = cp - zero ;
	      
	      if (strcasecmp (name(backClass) + 1, name(inf->fromModel) + 1))
		messerror ("Model error: In model %s, tag %s XREF's to tag %s in class %s\n"
			 "             reverse tag %s in model %s points to %s, not %s",
			 name(inf->fromModel), name(inf->fromTag), name(inf->toTag), name(inf->toModel),
			 name(inf->toTag), name(toModel), name(backClass), name(inf->fromModel)) ;
*/
	    }
	}
      else
	messerror ("Model error: In model %s, tag %s XREF's to tag %s in class %s\n"
		   "             can't find model %s",
		   name(inf->fromModel), name(inf->fromTag), name(inf->toTag), name(inf->toModel),
		   name(toModel)) ;
    }


  /* Any errors in the model file and we exit..no save attempted.            */
  if (messErrorCount() != errorCount) messExit ("Errors in model file") ;


  /* clean up memory */
  for (i = 1 ; i < lexMax(_VModel) ; ++i)
    assDestroy (array (tagsInModel, i, Associator)) ;
  arrayDestroy (xinfo) ;
  arrayDestroy (tagsInModel) ;

  return ks ;
} /* checkModelSyntax */


/********************************************************/

static void showModelNodeJaq (BS bs, int x, int *y,
			      KEY lastTag, int dLastTag)
{ 
  AC_HANDLE h ;

  if (!bs) return ;
  if (bs->key == _XREF)
    { if (!(bs = bs->right)) return ;
      if (!(bs = bs->right)) return ;
    }
  if (bs->key == _UNIQUE)
    if (!(bs = bs->right)) return ;
  if (bs->key == _REPEAT)
    return ;

  h = handleCreate() ;

  while (bs)
    { char *qname ;
      char *tname ;

      if (x == 1 && *y == 1)
	qname = "NAME" ;
      else if (!class(bs->key) && bs->key > _LastN)
	{ qname = strnew (messprintf ("EXISTS %s", name(bs->key)), h) ;
	  lastTag = bs->key ;
	  dLastTag = 0 ;
	}
      else
	{ if (!lastTag)
	    qname = strnew (messprintf ("[%d]", dLastTag+1), h) ;
	  else if (dLastTag > 0)
	    qname = strnew (messprintf ("%s[%d]", name(lastTag), dLastTag+1), h) ;
	  else
	    qname = name(lastTag) ; 
	  ++dLastTag ;
	}

      if (!class(bs->key))
	tname = name(bs->key) ;
      else if (pickType(bs->key) == 'A')
	tname = strnew (messprintf ("@%s", className(bs->key)), h) ;
      else if (KEYKEY(bs->key) == 1)
	tname = strnew (messprintf ("#%s", className(bs->key)), h) ;
      else
	tname = name(bs->key) ;

#ifndef NEW_MODELS_X
      if (class(bs->key) == _VModel)
	freeOutf ("%d %d |#%s|%s|\n", *y, x, name(bs->key)+1, qname) ;
      else
#endif
	freeOutf ("%d %d |%s|%s|\n", *y, x, tname, qname) ;
      showModelNodeJaq (bs->right, x + strlen(tname) + 1, y, lastTag, dLastTag) ;
      if ((bs = bs->down))
	++*y ;
    }

  handleDestroy (h) ;
} /* showModelNodeJaq */


void showModelJaq (KEY key)	/* key is in class Model */
{ 
  OBJ obj ;
  int x, y ;

#ifndef NEW_MODELS
  key = KEYMAKE(pickWord2Class(name(key)+1), 0) ;
#endif /* !NEW_MODELS */

  obj = bsCreate (key) ;

  if (!obj) return ;
  
  x = 1 ; y = 1 ;
  showModelNodeJaq (obj->root, x, &y, 0, 0) ;
  
  bsDestroy (obj) ;

  return;
} /* showModelJaq */

/********************************************************/
/********************************************************/


 
 
 
