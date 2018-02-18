/*  File: picksubs.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
     pickInit
       must be called once at start of program.
     pickWord2Class()
       returns a class number
     pickMatch matches a name to a template 
       accepts the wild char * (anything) and ? (single character).

 * HISTORY:
 * Last edited: Dec 14 14:50 1998 (fw)
 * * Sep 28 14:53 1992 (mieg): a fix because one tag has a - rather than _
 * * Apr  2 13:41 1992 (mieg): Simplified customization by defining
   the new wspec files (sys)options and displays.
 * Created: Thu Apr  2 13:41:10 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */

#include <ctype.h>
#include "acedb.h"
#include "pick.h"
#include "lex.h"
#include "session.h"
#include "bs.h"
#include "query.h"
#include "whooks/systags.h"
#include "whooks/sysclass.h"
#include "whooks/tags.h"
#include "utils.h"


#ifndef NON_GRAPHIC
#include "display.h"
#endif /* !NON_GRAPHIC */

extern void classInit (void) ;
extern void sysClassInit (void) ;
extern void sysClassOptions(void) ; 
extern void classSortAll(void) ;
extern BOOL READING_MODELS ;

/******************************************************************/

FREEOPT pickVocList[256] = { {0,  "Choose from"}} ;

/*
* pickList is a description of every class known to the system. It
* is indexed by class number.  
* class number 0 is the _VSystem class, it contains the acedb Tags
*/
PICKLIST pickList[256] ;

/******************************************************************/

static Stack pickStack = 0 ;
static Array clA = 0 ;
static Array classMasks = 0 ;
static Array classComposite = 0;
typedef struct { KEY mainClasse ; unsigned char mask ; }  CLASS_MASK ;

/********************************************/

char *className(KEY kk)  
{     
  return stackText(pickStack, pickList[class(kk)].name)  ;
}

/********************************************/

char *pickClass2Word (int classe)
{     
  return className(KEYMAKE(classe,0)) ;
}

/********************************************/

int superClass (KEY cl)
{ 
  pickIsA (&cl, 0) ;
  return cl ;
}

/********************************************/

/* Given the key that equates to "?<model>" (e.g. "?Sequence")
 * as found in the model for a class, return the class, returns 0 on failure.
 * Seems surprising this doesn't already exist, perhaps it does and I just
 * couldn't find it. */
int classModel2Class(KEY class_model_key)
{
  int class = 0 ;
  char *class_name ;

  class_name = name(class_model_key) ;
  if (*class_name == '?')
    {
      class_name++ ;
      class = pickWord2Class(class_name) ;
    }

  return class ;
}


/********************************************/
Array pickIsComposite(KEY classe)
{
  int c;

  if (classe > 256)
    {
      if (class(classe) != _VClass)
	messcrash("Non-class passed to pickIsComposite");
      c = KEYKEY(classe);
    }
  else
    c = classe;

  return array(classComposite, c, Array);
}

BOOL pickBelongsToComposite(KEY classe, KEY key)
{
  Array a = pickIsComposite(classe);
  int i;

  if (classe > 256)
    classe = KEYKEY(classe);
  
  if (!a)
    return classe == class(key);
  
  for (i=0; i<arrayMax(a); i++)
    if (KEYKEY(array(a, i, KEY)) == class(key))
      return TRUE;
  
  return FALSE;
}
     
/********************************************/

BOOL pickIsA (KEY *classp, unsigned char *maskp)
{
  KEY key ;

  if (maskp) *maskp = 0x00 ;  /* Normal situation */
  if (*classp > 256)
    { 
#ifdef OLDJUNK
      if ((obj = bsCreate(*classp)))
	{
	  if (bsGetData (obj,_Belongs_to_class, _Int, &topClass))
	    { n = 0 ;
	      bsGetData (obj,_Mask, _Int, &n) ;
	      *maskp = (unsigned char) n ;
	      *classp = topClass ;
	    }
	  bsDestroy(obj) ;
	}
#endif
      if (class(*classp) == _VClass && classMasks &&
	  (key = KEYKEY (*classp)) && key < arrayMax(classMasks))
	{ 
	  if (maskp) *maskp = arr(classMasks, key, CLASS_MASK).mask ;
	  *classp = arr(classMasks, key, CLASS_MASK).mainClasse ;
	  return TRUE ;
	}
      *classp = 0 ;
      return FALSE ;
    }

  return TRUE ;
}

/******************************************************************/

static void pickGetClassNames(void)
{ 
  int  k , v ;
  KEY key = 0 ;

  pickStack = stackReCreate(pickStack, 200) ;
  pushText(pickStack,"System") ; /* prevents name = 0 */


  clA = arrayCreate(20, int) ;
  v = 0 ;

  while(lexNext (_VMainClasses, &key))
    { k = KEYKEY(key) ;
      if (*name(key) != '_')
	pickList[k].name = stackMark(pickStack) ;
      pushText(pickStack, name(key) );
      array(clA, v++, int) = k ;
    }
} 

/******************************************************************/

static KEYSET classAppearance = 0 ;
static FREEOPT classOptions[] =
{
  {14, "class options"},
  {'h', "Hidden"},
  {'v', "Visible"},
  {'a', "Array"},
  {'b', "Btree"},
  {'x', "XREF"},
  {'d', "Display"},
  {'t', "Title"},
  {'s', "Symbol"},
  {'r', "Rename"},
  {'c', "CaseSensitive"},
  {'p', "Protected"},
  {'S', "Sybase"},
  {'k', "Known"},
  {'m', "Map"}
  } ;

int pickMakeClass (const char* cp)
{
  KEY key, treeKey, model, class, classeKey ;
  int k, k1 ;

  lexaddkey (cp, &key, _VMainClasses) ;   k = KEYKEY (key) ;
  lexaddkey (cp, &classeKey, _VClass) ;   k1 = KEYKEY (classeKey) ;
  if (!classMasks) classMasks = arrayCreate (128, CLASS_MASK) ;
  array (classMasks, k1, CLASS_MASK).mainClasse = k ;

  if (!classComposite) classComposite = arrayCreate (256, Array) ;
  array (classComposite, k1, Array) = NULL ;

  if (_VModel)
    {
#ifdef NEW_MODELS
      lexaddkey (messprintf("#%s", name(key)), &model, _VModel) ;
#else
      lexaddkey (messprintf("?%s", name(key)), &model, _VModel) ;
#endif
      pickList[k].model = model ;
    }

  if (!pickList[k].name || !pickList[k].classe ||
      strcmp(stackText(pickStack, pickList[k].name),cp))
    { pickList[k].name = stackMark(pickStack) ;
      pushText(pickStack, cp) ;
      lexaddkey(cp, &class, _VClass) ;
      lexAlias(&class, cp, FALSE, FALSE) ; /* rename it to fix case */

      pickList[k].type = 'B' ;	/* default is B */
      pickList[k].visible = 'v' ; /* for  B default is visible */
      /* attention, model.c redeclares them as V */
      lexaddkey("TREE", &treeKey, _VDisplay) ;
      pickList[k].displayKey = treeKey ; /* for B default is TREE */
      pickList[k].classe = classeKey ;
      pickList[k].isCaseSensitive = FALSE;
      pickList[k].updateNames = FALSE;
    }

      /* all B classes should be listed in _VModel with the constructed types */

  return k ;
}
	
/******************************************************************/

static void pickGetClassProperties (BOOL pass)
{
  int  k, nc, i ;
  int line = 0 ;
  char *cp = sessionFilName ("wspec/options","wrm", "r") ;
  FILE * fil = cp ? filopen(cp, 0,"r") : 0 ;
  KEY tag , option, dummy ;

  if (!fil)
    messcrash("pickInit cannot find Class definition file wspec/options.wrm") ;
  if (!classAppearance)
    { nc = 0 ; classAppearance = keySetCreate() ; }
  else
    nc = keySetMax(classAppearance) ;
 
  freespecial ("\n\t\"\\") ;
  while(line ++ ,freeread(fil))
                                                /* read #define */
    if ((cp = freeword()) && *cp++ == '_' && *cp++ == 'V')
      {
	if(!*cp)
	  messcrash("Isolated _V while parsing line %d of wspec/options.wrm", line ) ;

	if (!strcmp (cp, "System"))
	  k = 0 ;
	else if (!lexword2key(cp, &dummy, _VMainClasses))
	  continue ; /* disregard option of classes without models */
	else
	  k = pickMakeClass (cp) ;

	for (i = 0 ; i < nc ; ++i)
	  if (keySet(classAppearance, i) == k)
	    break ;
	if (i == nc)		/* not found */
	  keySet(classAppearance, nc++) = k ;

	while (TRUE)
	  {
	    freenext() ;
	    if (!freestep('-')) 
	      {
		if ((cp = freeword()))
		  messcrash ("In %options.wrm line %d, no - at start of option %s", 
			     line, freepos()) ;
		else
		  break ;
	      }
	    if (!freekey(&option, classOptions))
	      messcrash ("In options.wrm line %d, unknown option %s", line, freepos()) ;
	    
	    switch (option) 
	      {
	      case 'h':
		pickList[k].visible = 'H' ;  /* default is visible */
		break ;
	      case 'v':
		pickList[k].visible = 'V' ;  /* default is visible */
		break ;
	      case 'a':
		pickList[k].type = 'A' ;
		break ;
	      case 'b':
		pickList[k].type = 'B' ;
		break ;
	      case 'x':
		pickList[k].Xref = TRUE ;
		pickList[k].visible = 'H' ; 
		pickList[k].type = 'B' ;
		break ;
	      case 'd':                       /* Display type (enum) */
		freenext() ;
		cp = freeword() ;
#ifndef NON_GRAPHIC		/* in non-graph mode we're not interested in displays */
		if (pass && !lexword2key (cp, &(pickList[k].displayKey), _VDisplay))
		  messcrash ("Bad display type %s in line %d of options.wrm", 
			     cp, line) ;
#endif
		break ;
	      case 't':
                               /* Tag for Name expansion  */
		pickList[k].tag = 0 ;
		if ((cp = freeword()))
		  { 
		    lexaddkey(cp, &tag, _VSystem) ; 
		    pickList[k].tag = tag ;
		  }
		break ;
	      case 's':
                               /* Tag for Symbol in maps  */
		pickList[k].symbol = 0 ;
		if ((cp = freeword()))
		  { 
		    lexaddkey(cp, &tag, _VSystem) ; 
		    pickList[k].symbol = tag ;
		  }
		break ;
	      case 'r':   /* other class name */
		if ((cp = freeword()))
		 { if (!pickList[k].alias)  /* Always keep the compilation value */
		     pickList[k].alias = pickList[k].name ;
		   pickList[k].name = stackMark(pickStack) ;
		   pushText(pickStack, cp) ;
		 }
		break ;
	      case 'c':
		/* notice that we cannot backtrack on a casesensitive class */
		pickList[k].isCaseSensitive = TRUE ;
		break ;
	      case 'p':
		pickList[k].protected = TRUE ;
		break ;
	      case 'S':
		pickList[k].sybase = TRUE ;
		break ;
	      case 'k':
		pickList[k].known = TRUE ;
		break ;
	      case 'm':
		freenext() ;
		cp = freeword() ;
		if (!lexword2key (cp, &(pickList[k].mapKey), _VDisplay))
		  messcrash ("Bad map type %s in line %d of options.wrm\n", cp, line) ;
		break ;
	      }
	  }
      }
  filclose ( fil ) ;
} /* pickGetClassProperties */

/******************************************************************/
   /* called from whooks/sysclass.c */ 
void pickSetClassOption (const char *nom, char type,  char v, 
			    const char *disp, BOOL protected, BOOL private,
			    BOOL case_sensitive, BOOL update_names )
{
  int i, k, nc ;

  if (!classAppearance)
    { nc = 0 ; classAppearance = keySetCreate() ; }
  else
    nc = keySetMax(classAppearance) ;

  if (!strcmp (nom, "System"))
    k = 0 ;
  else
    k = pickMakeClass (nom) ;

  for (i = 0 ; i < nc ; ++i)
    if (keySet(classAppearance, i) == k)
      break ;
  if (i == nc)		/* not found */
    keySet(classAppearance, nc++) = k ;
  
  if (!pickList[k].visible || pickList[k].visible == ace_lower(pickList[k].visible))
    pickList[k].visible = v ;
  pickList[k].type = type ;
  if (type == 'X')
    { pickList[k].Xref = TRUE ;
      pickList[k].type = 'B' ;
    }
  pickList[k].isCaseSensitive = case_sensitive;
  pickList[k].updateNames = update_names;

  if (disp) lexaddkey (disp, &(pickList[k].displayKey), _VDisplay) ;
  pickList[k].protected = protected ;  /* do not parse */
  pickList[k].private = private ;      /* do not dump */
#ifdef NEW_MODELS
  if (pickList[k].type == 'B') 
    lexaddkey (messprintf("#%s", nom), &pickList[k].model, _VModel) ;
#else
  lexaddkey (messprintf("?%s", nom), &pickList[k].model, _VModel) ;
#endif
  lexaddkey (nom, &pickList[k].classe, _VClass) ;
} /* pickSetClassOption */

/******************************************************************/

  /* Try name, then alias = compilation name */
int  pickWord2Class(const char *word)     /* returns class of word */
{
  register int i = 256 ;

  while(--i && lexstrcmp(word,stackText(pickStack, pickList[i].name))) ;
  if (!i)
    { i = 256 ;
      while(--i && 
	    lexstrcmp(word, stackText(pickStack, pickList[i].alias)));
    }
 
  return i ;
} /* pickWord2Class */

/**************************************************************/

static void pickDefaultProperties(void)
{
  KEY treeKey ;
  int j, k ;
  
  lexaddkey("TREE", &treeKey, _VDisplay) ;
      /* clA is used to have the Visible Class in main menu in
       * same order as they appear in classes.wrm.
       */
  for (j = 0 ; j < arrayMax(clA) ; j++)
    { k = array(clA, j, int) ;
      if (pickList[k].name)
	{
	  if (!pickList[k].type)            /* default is B */
	    pickList[k].type = 'B' ;
	  if (!pickList[k].displayKey)
	    {
	      if (pickList[k].type == 'B')
		pickList[k].displayKey = treeKey ;   /* for B default is TREE */
	      else
		pickList[k].displayKey = 0 ;   /* for A default is hidden */
	    }
	}
    }

  arrayDestroy(clA);

  return;
} /* pickDefaultProperties */

/******************************************************************/

/*
Actually this function should be here, but i moved it to whooks/class.c
because that file is NON_GRAPHIC_VARIANT and i need graphWindowSize
sorry, this is not very clean, meig may 97


void pickRememberDisplaySize (char *display)
{
  float sx, sy, sw, sh ;
  if (graphWindowSize (&sx, &sy, &sw, &sh))
    pickSetDisplaySize (display, sx, sy, sw, sh) ;
}

*/

/*******************************************************/

void pickRegisterConstraints (void)
{ OBJ obj ;
  Stack s = 0 ;
  char *cp ;
  KEY key = 0 ;
  KEY mainClasse = 0 ;
  
  while (lexNext (_VClass, &key))
    { 
      mainClasse = 0 ;
      if ((obj = bsCreate (key)))
	{ if (bsGetData (obj, _Constraints, _Text, &cp))
	    { if (!s)
		s = stackCreate (50) ;
	      pushText (s, messprintf ("(%s)", cp)) ;
	      while (bsGetData (obj, _bsDown, _Text, &cp))
		catText (s, messprintf ("AND ( %s )", cp)) ;
	      mainClasse = key ;
	      pickIsA (&mainClasse, 0) ;
	    }
	  bsDestroy (obj) ;
	}
      if (s && stackMark (s) && mainClasse &&
	  !queryConstraintsInit (stackText(s,0), mainClasse))
	messcrash ("Syntax error in file wspec/constraints.wrm, class %s\n: %s\n", 
		   pickClass2Word (mainClasse), stackText(s,0)) ;
      if (s) stackClear (s) ;
    }
  stackDestroy (s) ;
} /* pickRegisterConstraints */

/********************************************************/
/******************* Public Routine *********************/
/********************************************************/

/* pickPreInit is a duplicate of pickInit
it could certainly be simplified but this works like that
so i don t care (mieg, oct 18,94)
*/


void pickPreInit(void)
{ 
  BOOL old = READING_MODELS ;
  register int k ;
  
  READING_MODELS = TRUE ;

  k = 256 ;
  while(k--)
    pickList[k].type = 0;
  pickList[_VGlobal].type = 'A' ;
  pickList[_VVoc].type = 'A' ;

  sysClassInit () ;
  pickList[_VClass].type = 'B' ; /* needed for pickMakeClass() to work */
  pickList[_VDisplay].type = 'B' ;
  classInit () ;

  pickGetClassNames() ;

  pickMakeClass ("Session") ;

  pickGetClassProperties(FALSE) ;	/* options.wrm */
  sysClassOptions() ; 

  pickDefaultProperties() ;

  READING_MODELS = old ;

  return;
} /* pickPreInit */



void pickInit(void)
{
  BOOL old = READING_MODELS ;
  register int k ;
  
  READING_MODELS = TRUE ;

  if (classAppearance)
    keySetReCreate(classAppearance);

  k = 256 ;
  while(k--)
    { 
      memset(&pickList[k],0,sizeof(pickList[k]));
      pickList[k].type = 0;
      pickList[k].name = 0;
      pickList[k].alias = 0;
      pickList[k].visible = 0;
      pickList[k].mask = 0;
      pickList[k].model = 0;
      pickList[k].displayKey = 0;
      pickList[k].protected = 0;
      pickList[k].Xref = FALSE;
    }
  pickList[_VGlobal].type = 'A' ;
  pickList[_VVoc].type = 'A' ;

  sysClassInit () ;
  pickList[_VClass].type = 'B' ; /* needed for pickMakeClass() to work */
  pickList[_VDisplay].type = 'B' ;
  classInit () ;

  pickGetClassNames() ;

#ifndef NON_GRAPHIC
  pickDefaultDisplays() ;
  pickGetDisplayTypes() ;
  pickSetGraphTypes () ;
#endif

  pickGetClassProperties(TRUE) ;
  sysClassOptions() ;

  pickDefaultProperties() ;
  pickRegisterConstraints () ;

  READING_MODELS = old ;

  return;
} /* pickInit */

/******************************************/


/********************************************************/

static void pickCheckClassDefinitions(void)
{ 
  KEY metaClass = 0, _Buried ;
  int n = 256 ; 
  OBJ obj; 

#ifdef NEW_MODELS
  if (!lexword2key("#Class", &metaClass, _VModel))
    return ;   /* i.e. the Class model does not yet exists   */
#else
  if (iskey(KEYMAKE(_VClass,0)) != 2)
    return ;   /* i.e. the Class model does not yet exists   */
#endif
  lexaddkey ("Buried", &_Buried, 0) ;
  while (n--)
    { if (!pickList[n].name)
	continue ;

      lexaddkey(pickClass2Word(n), &metaClass, _VClass) ;
      obj = bsUpdate(metaClass) ;

	  /* too early, .model is not yet esatablished */
	  if (!pickList[n].visible)
	    {
	      if (pickList[n].type == 'B' &&  /* for  B default is visible */
		  iskey(pickList[n].model)) /* iff the model exists */
		pickList[n].visible = 'v' ;   
	      else
		pickList[n].visible = 'h' ;   /* for  A default is hidden */
	    }

      if (ace_lower(pickList[n].visible) == 'v')
	bsAddTag(obj, _Visible) ;   
      else if (pickList[n].visible == 'H')
	bsAddTag(obj, _Hidden) ; /* show in the triangle */
      else
	bsAddTag(obj, _Buried) ;    /* show nowhere */

#ifdef ACEDB4
      /* KLUDGE, mieg sept 96
	 there is clearly a bug here, 
	 consisder a subclass sub of class c
	 sub should not have an entry in pickList and mainclasses
         furthermore, without the kludge, Belongs to class 
	 is reset in each newe session to sub then to c
	 implying always some locked obj, which imply
	 taht cache 1 has to get write access silently
         to prevent complete locking

	 this kludge avoids the problem, but the code
	 is not clean as we know in this area
	 */
      if (!bsFindTag(obj, _Is_a_subclass_of)) /* KLUDGE */
	/* compiled_as is used in bindex.c before sep 1998
	 * Belongs_to_class is used in pickIsA before sep 1998
         * it is therefore essential to position them correctly
	 * untill ACEDB5
	 */
	{ 
	  bsAddData(obj, _Compiled_as,  _Text, 
		  stackText(pickStack, 
			    pickList[n].alias ? pickList[n].alias : pickList[n].name
			    )) ;
          bsAddData(obj, _bsRight, _Int, &n) ;
	  bsAddData(obj, _Belongs_to_class, _Int, &n) ;
	}
#endif
      bsSave(obj) ;
    }
} /* pickCheckClassDefinitions */

/**********************************************************************/

static void  pickRegisterClassMasks (void)
{
  KEY key = 0 ;
  int i, n, k ;
  OBJ obj = 0 ;

  while (lexNext (_VClass, &key))
    {
      obj = bsCreate (key) ; k = KEYKEY (key) ;
      if (obj)
	{
	  if (bsGetData(obj, _Belongs_to_class, _Int, &n))
	    array(classMasks, k, CLASS_MASK).mainClasse  = n ;
	  if (bsGetData (obj, _Mask, _Int, &n))
	    array(classMasks, k, CLASS_MASK).mask  = n ;

	  if (bsFindTag(obj, str2tag("Composite_of")))
	    {
	      Array a = arrayCreate(10, BSunit);
	      KEYSET aa = keySetCreate () ;
	      if (bsFlatten(obj, 1, a))
		{
		  for (i = 0 ; i < arrayMax (a) ; i++)
		    keySet (aa, i) = arr (a, i, BSunit).k ;
		    
		  array(classComposite, k, Array) = aa ;
		  /* Do not make any real objects in this class */
		  pickList[k].protected = TRUE;
		}
	      else
		arrayDestroy(aa);
	      arrayDestroy (a) ;
	    }
	  bsDestroy (obj) ;
	}
    }      
} /* pickRegisterClassMasks */

/**********************************************************************/

void pickVisibility(void)
{
  KEY  *classp, metaClass ;
  OBJ obj = 0 ;
  KEYSET subclasses = 0 ;
  int v = 0, i ;

  i= keySetMax(classAppearance) ;
  if(!i)
    return ;

  subclasses = query(0, ">?Class Is_a_subclass_of") ;

  classp = arrp(classAppearance,0,KEY) - 1 ;
  while(classp++, i--)
    if (lexword2key(pickClass2Word(*classp), &metaClass, _VClass) &&
	( obj = bsCreate(metaClass)) )
      { if (bsFindTag(obj, _Visible))
	  {
	    v++ ;
	    pickVocList[v].text = name(metaClass) ;
	    pickVocList[v].key = metaClass ;	  
	  }
	
	bsDestroy(obj) ;
	if (v >= 254)
	  break ;
      }

  i= keySetMax(subclasses) ;
  
  classp = arrp(subclasses,0,KEY) - 1 ;
  while(classp++, i--)
    if ((obj = bsCreate(*classp)))
	{ if (bsFindTag(obj, _Visible))
	    {
	      v++ ;
	      pickVocList[v].text = name(*classp) ;
	      pickVocList[v].key = *classp ;	  
	    }

	  bsDestroy(obj) ;
	  if (v >= 254)
	    break ;
	}

  keySetDestroy(subclasses) ;
  pickVocList[0].key = v ;
  pickVocList[0].text =  "Choose a class" ;
  if (v < 1)
    messout("This database has no visible class.\n"
            "Most probably, you should study and edit\n"
	    "the files wspec/options.wrm and wspec/models.wrm,\n",
	    "then use the option Read-Models from main menu.") ;
}

/********************************************************/

void pickCheckConsistency(void)
{
  pickCheckClassDefinitions() ;
  pickRegisterClassMasks () ;
  pickVisibility() ;
  classSortAll() ;
} /* pickCheckConsistency */

/*************************************************************/
/*                   messStatus 
   A routine that is called by variuos functions which might
   take longer than normal or require the user interact with
   a specific window.
   The message is intended to be displayed in the main window
   where the user can check what's going on in case (s)he has
   lost track of what other windows are doing.

   To be moved to an ace-general  module, like acesubs.c or so, 
   It to be in mainpick.c (the graphical bit) and in 
   session.c (undef'd for non-graph compilations) 

   Currently only graphical code that goes thorugh 
   acedbGraphInit() will call messStatusRegister() in order
   to register messStatusDisplayInWindow() (in mainpick.c)
*/
/*************************************************************/

OutRoutine messStatusRoutine = 0;

void messStatus (char* text)
{
  if (messStatusRoutine)
    (*messStatusRoutine)(text);

  return ;
} /* messStatus */

OutRoutine messStatusRegister (OutRoutine func)
{
  OutRoutine old = messStatusRoutine ;
  messStatusRoutine = func ; 
  return old;
} /* messStatusRegister */


/********************* eof ***********************************/
