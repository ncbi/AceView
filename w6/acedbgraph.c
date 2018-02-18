/*  File: acedbgraph.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Initialises the graph acedb interface which allows acedb
 *              to alter the graph package behaviour for certain required
 *              acedb functions (window position/printer preference etc).
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  4 12:02 1999 (fw)
 * Created: Thu Sep 24 13:55:44 1998 (edgrif)
 *-------------------------------------------------------------------
 */

#include "acedbgraph.h"					    /* ace -> graph header. */
#include "graphAcedbInterface.h"			    /* graph -> ace header. */
#include "acedb.h"					    /* For messStatus & register */
#include "main.h"					    /* For messStatusDisplayInWindow() */
#include "display.h"					    /* For displayCreate */
#include "disptype.h"					    /* For DtChrono etc. */
#include "session.h"					    /* for getLogin */
#include "freeout.h"					    /* for freeOutNN */
#include "colcontrol_.h"				    /* Sadly we need this for the
							       VIEWWINDOW structure. */
#include "lex.h"					    /* for lexaddkey */
#include "sysclass.h"					    /* for _VMainClasses */
#include "systags.h"					    /* for _Float */
#include "bs.h"						    /* For lots of bs calls. */


/* The following is lifted from graphxt.c, it would seem that for non-acedb  */
/* code the wcs stuff is never included, we always get the below local void  */
/* functions. For acedb code the real wcs functions are only included (from  */
/* here and in wcscomm.c if WCS is defined, currently this is not defined    */
/* in our makefiles....in other words only the void routines below are       */
/* ever really used by us.                                                   */
/*                                                                           */

#if defined(WCS) && !defined(NO_WCS)

#include "wcomm.h"

#else

static void wcommInit (Display *d, Window w) { }
static void wcommFinish (void) { }
static BOOL WPropMessage (XEvent *ev) { return FALSE ;}

#endif




/* This MUST be static so it persists between calls.                         */
/* This array has been moved from graphprint.c with the following comment    */
/* explaining why the array had to be changed from the default one in        */
/* graphprint.c                                                              */
/*                                                                           */
/* in ACEDB people got confused with OK  being the top 
     of the menu in that window, as  every other ACEDB 
     window has Quit there that does nothing and people 
     ended up with unwanted printouts */
static MENUOPT pdMenu[] =
{
  {graphDestroy, "Quit"},
  {pdExecute, " OK "},
  {0, 0}
} ;


static char colorDefault (void) ;
static void getFonts(FILE **fil, Stack fontNameStack, int *nFonts) ;
static void classPrint(FILE *fil, KEY key) ;
static void GetXlibFonts(int *level) ;

static void setAnon(VIEWWINDOW view, COLCONTROL control) ;
static BOOL setSaveView(MENU menu) ;
static void setSaveViewMsg(MAPCONTROL map, VIEWWINDOW view) ;
static void setViewMenu(MENUSPEC **menu) ;


/* Document routines moved from viewedit                                     */


/* Interface routines for colcontrol.c                                       */
static void resetkey(KEY *key) ;
static void resetsave(VoidCOLINSTANCEOBJ *savefunc) ;
static void menufree(int box, FREEOPT *viewmenu) ;
static void maplocate(OBJ init, VoidCOLINSTANCEOBJ *savefunc, MAPCONTROL map) ;
static void scale(OBJ init, VoidCOLINSTANCEOBJ *savefunc, MAPCONTROL map) ;
static void space(OBJ init, VoidCOLINSTANCEOBJ *spacefunc, SPACERPRIV private) ;
static void addmap(char *name, KEY *key) ;

/* These routines were in colcontrol.c and have been shifted here.           */
static KEY controlSaveView(MAPCONTROL map);
static BOOL controlSetView(MAPCONTROL map, KEY view);
static Array addToMenu(Array menu, AC_HANDLE handle, KEY view) ;
static void viewMenuFunc(KEY key, int box) ;
static void mapLocatorSave(COLINSTANCE instance, OBJ obj) ;
static void mapScaleSave(COLINSTANCE instance, OBJ obj) ;
static void spacerSave(COLINSTANCE instance, OBJ obj) ;




/* This routine initialises the graph package and the graph/acedb interface, */
/* it encapsulates the two so that they occur in the right order and with    */
/* the correct parameters. (In fact the order is not important at the moment */
/* but this may change later).                                               */
void acedbAppGraphInit(int *argcptr, char **argv)
  {

  /* Initialise the graphics package and register any overloaded functions.  */
  graphInit (argcptr, argv) ;	/* before calls to messcrash() */
							    /* Is this messcrash comment right ? */

  /* Initialise the acedb/graph package interface so that graph can do acedb */
  /* specific displays, this must follow immediately after graphInit().      */
  acedbGraphInit() ;

  return ;
  }


/* This routine is called by acedb applications which use the graph package  */
/* to set up the graph package so that it will do acedb specific things.     */
/*                                                                           */
void acedbGraphInit(void)
  {
  char *login ;

  setGraphAcedbMainActivity(messStatus) ;

  /* Register stuff for creating various types of graph but in acedb style.  */
  setGraphAcedbDisplayCreate(displayCreate, DtChrono, DtColControl, DtFile_Chooser) ;


  /* Register stuff for graphprint                                           */
  login = strnew(getLogin(TRUE), 0) ;
  setGraphAcedbPDMenu(&pdMenu[0], colorDefault, login) ;

  /* Register stuff for graphps                                              */
  setGraphAcedbPS(thisSession.name, getFonts) ;


  /* Set style of graphSelect window.                                        */
  graphSetSelectDisplayMode(GRAPH_SELECT_PLAIN) ;


  /* Register boxInfo classname print routine.                               */
  setGraphAcedbClassPrint(classPrint) ;

  setGraphAcedbXFont(GetXlibFonts) ;

  setGraphAcedbXWcomm(wcommInit, wcommFinish, WPropMessage) ;


  /* Register viewedit routines.                                             */
  setGraphAcedbView(setAnon, setSaveView, setSaveViewMsg, setViewMenu) ;

  /* Register colcontrol routines.                                           */
  setGraphAcedbColControl(controlSetView, resetkey, resetsave,
			  menufree, maplocate, scale, space, addmap) ;

  messStatusRegister (messStatusDisplayInWindow);

  return ;
  }



/* Get PS fonts specific to an acedb database (used in printing files).*/
/*                                                                     */
static char colorDefault ()
{ FILE *fil ;
  int level ;
  char *cp ;
  BOOL isColor = FALSE ;

  if ((fil = filopen ("wspec/psfonts","wrm","r")))
    { level = freesetfile(fil,"") ;
      freespecial("\n\"/\\") ;
      
      while (freecard(level))
	while ((cp = freeword ()))
	  if (!strcmp(cp,"COLOUR") || !strcmp (cp,"COLOR"))
	    isColor = TRUE ;
      /* filclose (fil) ; already closed by freecard at eof */
    }
  return isColor ? 'c' : 'p' ;
}


/* Get PS fonts specific to an acedb database (used in saving files). */
/*                                                                    */
static void getFonts(FILE **fil, Stack fontNameStack, int *nFonts)
{

  if ((*fil = filopen("wspec/psfonts","wrm","r")))
    {
    int level = freesetfile(*fil,"") ;
    char *cp ;

    freespecial("\n\"/\\") ;
      
    while (freecard(level)) /* will close fil */
      while ((cp = freeword ()))
	if (!strcmp(cp,"COLOUR") || !strcmp (cp,"COLOR")) ; /* do nothing */
	else
	  {
	  pushText(fontNameStack,cp) ;
	  (*nFonts)++ ;
	  }
      /* filclose (fil) ; already closed by freecard at eof */
    }

  return ;
  }


/* This code has been extracted from graphsub.c as it is acedb specific.     */
/* Lincoln Stein made some adjustments to the code so I've kept his date     */
/* stamp below.                                                              */
/*                                                                           */
static void classPrint(FILE *fil, KEY key)
  {
  /* 2 Mar 1998 LS, so that boxes can go to gifaceserver */
  /* 2 Mar 1998 LS */
  if (fil != NULL) fprintf (fil, "  %s:%s", className(key), freeprotect(name(key))) ;
  else freeOutf ("  %s:%s", className(key), freeprotect(name(key))) ;

  return ;
  }



/* This code was extracted from graphxlib.c, it attempts to find the acedb   */
/* fonts file and sets this as input for freesubs calls in the graph         */
/* package.                                                                  */
/*                                                                           */
static void GetXlibFonts(int *level)
  {
  FILE *f = filopen ("wspec/xfonts","wrm","r") ;

  if (!f)
    {
    messerror("cannot open the Acedb Xfonts file wspec/xfonts.wrm."
	      " The MIT default fonts will be used instead") ;
    }
  else
    {
    *level = freesetfile(f,"") ;
    freespecial ("\n/\\") ;
    }

  return ;
  }



/*****************************************************************************/
/* Viewedit functions.                                                       */

/* Small function for viewedit, sets altered views to be anonymous.          */
/*                                                                           */
static void setAnon(VIEWWINDOW view, COLCONTROL control)
  {
  int i;

  /* if we've edited and not saved, make views anonymous */
  if (view->dirty)
    for (i=0; i<arrayMax(control->maps); i++)
      arr(control->maps, i, MAPCONTROL)->viewKey = 0;

  return ;
  }


static BOOL setSaveView(MENU menu)
  {
  BOOL result = FALSE ;

  if (writeAccessPossible())
    {
    menuUnsetFlags(menuItem(menu,"Save view"), MENUFLAG_DISABLED);
    menuUnsetFlags(menuItem(menu,"Save view as default"), MENUFLAG_DISABLED);
    result = TRUE ;
    }

  return result ;
  }

static void setSaveViewMsg(MAPCONTROL map, VIEWWINDOW view)
  {
  if (map->viewKey && !view->dirty)
    graphText(messprintf("View stored in: %s", name(map->viewKey)), 1, 1);
  else
    graphText("View not stored.", 1, 1);

  return ;
  }


/* All of these routines seem to be specific to acedb code...                */
/* They are all set in motion by a menu...so the menu is the key...          */

/* only called by other acedb-only routines...see below...                   */
/* Need to make currentView,viewWindowDraw  non-static,                      */
/* otherwise looks ok...                                                     */
/*                                                                           */
static void viewSaveView(KEY tag)
{
  VIEWWINDOW view = currentView("viewSaveView");
  MAPCONTROL map = view->currentMap;
  KEY viewKey;
  OBJ obj;
  KEY _View;
  
  lexaddkey("View", &_View, 0);
  
  viewKey = controlSaveView(map);
  if (!viewKey)
    return;

  view->dirty = FALSE;
  viewWindowDraw(view);
  
  obj = bsUpdate(map->key);
  if (bsAddTag(obj, _View))
    bsAddKey(obj, _View, viewKey);
  else
    messerror("No tag View in object %s - cannot complete save", 
	      name(map->key));

  if (tag)
    { if (bsAddTag(obj, tag))
	bsAddKey(obj, tag, viewKey);
      else  
	messerror("No tag %s in object %s - cannot complete save", 
		  name(tag),
		  name(map->key));
    }
  bsSave(obj);
}


/* only referenced in the acedb menu....only uses above function.            */
static void viewSaveMap(MENUITEM m)
{ viewSaveView(0);
}

/* only referenced in the acedb menu....only uses above funciton.            */
static void viewSaveAsDefault(MENUITEM m)
{ KEY _Default_view; 
  lexaddkey("Default_view", &_Default_view, 0);
  viewSaveView(_Default_view);
}

/* only referenced in the acedb menu....                                     */
static void loadNamedView(MENUITEM m)
{
  KEY key, _VView;
  VIEWWINDOW view = currentView("loadNamedView");
  MAPCONTROL map = view->currentMap;
  char buffer[41];

  lexaddkey("View", &key, _VMainClasses);
  _VView = KEYKEY(key);
  
  if (!messPrompt("Load view named :", "", "t"))
    return;
  
  strncpy(buffer, freeword(), 40);
  if (!lexword2key(buffer, &key, _VView))
    { messout(messprintf("Can't find %s", buffer));
      return;
    }

  (void)controlSetView(map, key);

  view->dirty = FALSE;
  viewWindowDraw(view);
}


/*                                                                           */
/*                                                                           */
/* This is a very dangerous interface in that there is no checking of whether*/
/* we run off the end of the MENUSPEC array, this should only be altered in  */
/* conjunction with the function makeViewMenu in viewedit.c.                 */
/*                                                                           */
static void setViewMenu(MENUSPEC **curr)
  {
  MENUSPEC *i = *curr ;					    /* easier to follow. */

  i++ ;							    /* Go past the current entry. */
  i->f = (MENUFUNCTION)help ;
  i->text = "Help" ;
 
  i++ ;
  i->f = (MENUFUNCTION)menuSpacer ;
  i->text = "" ;

  i++ ;
  i->f = loadNamedView ;
  i->text = "Load named view" ;

  i++ ;
  i->f = viewSaveMap ;
  i->text = "Save view" ;

  i++ ;
  i->f = viewSaveAsDefault ;
  i->text = "Save view as default" ;

  *curr = i ;						    /* Update pointer to last entry. */

  return ;
  }


/*****************************************************************************/
/* ColControl routines.                                                      */


/* Trivial, but helps retain semantics of original colcontrol code where     */
/* this was #def'd for ACEDB only.                                           */
/*                                                                           */
static void resetkey(KEY *key)
  {
  *key = 0 ;

  return ;
  }

static void resetsave(VoidCOLINSTANCEOBJ *savefunc)
  {
  savefunc = NULL ;

  return ;
  }


static void menufree(int box, FREEOPT *viewmenu)
  {

  graphBoxFreeMenu(box, viewMenuFunc, viewmenu);

  return ;
  }


static void maplocate(OBJ init, VoidCOLINSTANCEOBJ *savefunc, MAPCONTROL map)
  {
  KEY _Magnification, _Projection_lines_on;

  *savefunc = mapLocatorSave ;
 
  if (init)
    {
    lexaddkey("Magnification", &_Magnification, 0);
    lexaddkey("Projection_lines_on", &_Projection_lines_on, 0);
      
    if (bsFindTag(init, _Projection_lines_on))
      map->hasProjectionLines = TRUE;
    else
      map->hasProjectionLines = FALSE;

    if (bsGetData(init, _Magnification, _Float, &map->magConf))
      map->mag = map->magConf;
    }

  }

static void scale(OBJ init, VoidCOLINSTANCEOBJ *savefunc, MAPCONTROL map)
  {
  KEY _Cursor_unit, _Cursor_on, _Scale_unit;
  float f;

  *savefunc = mapScaleSave ;

  if (init)
    {
    lexaddkey("Cursor_unit", &_Cursor_unit, 0);
    lexaddkey("Cursor_on", &_Cursor_on, 0);
    lexaddkey("Scale_unit", &_Scale_unit, 0);
    if (bsGetData(init, _Cursor_unit, _Float, &f))
      map->cursor.unit = f;
    if (bsGetData(init, _Scale_unit, _Float, &f))
      map->scaleUnit = f;
    if (bsFindTag(init, _Cursor_on))
      map->hasCursor = TRUE;
    }

  return ;
  }


static void space(OBJ init, VoidCOLINSTANCEOBJ* spacefunc, SPACERPRIV private)
  {
  KEY _Spacer_width, _Spacer_colour;
  float w;

  *spacefunc = spacerSave;

  if (init)
    {
    lexaddkey("Spacer_width", &_Spacer_width, 0);
    lexaddkey("Spacer_colour", &_Spacer_colour, 0);
    if (bsGetData(init, _Spacer_width, _Float, &w))
      private->width = w;
    if (bsFindTag(init, _Spacer_colour))
      private->colour = controlGetColour(init);
    }

  return ;
  }


static void addmap(char *name, KEY *key)
  {

  lexaddkey(name, key, 0) ;

  return ;
  }


/* These routines from from colcontrol.c, they are acedb specific.           */

/* Read a view and make the columns and put them in map. */
static BOOL controlSetView(MAPCONTROL map, KEY view)
  {
  OBJ View  = bsCreate(view);
  KEY _Columns, _Submenus, _Cambridge, _Nobuttons, _HideHeader, colKey;
  COLPROTO proto=0;
  int i, disp;
  BSMARK mark = 0;
  COLCONTROL control = map->control;
  COLINSTANCE instance;
  char *s;

  lexaddkey("Columns", &_Columns, 0);
  lexaddkey("Submenus", &_Submenus, 0);
  lexaddkey("Cambridge", &_Cambridge, 0);
  lexaddkey("No_buttons", &_Nobuttons, 0);
  lexaddkey("Hide_header", &_HideHeader, 0);

  if (!bsFindTag(View, map->viewTag) || !bsGetData(View, _Columns, _Text, &s))
    { messout("Cannot use View object %s, Tag %s not set, "
	      "or no columns found.", 
	      name(view), name(map->viewTag));
      bsDestroy(View);
      return FALSE;
    }

/* Delete any existing column instances for this map */   
  for (i=0; i<arrayMax(control->instances); i++)
    { instance = arr(control->instances, i, COLINSTANCE);
      if (instance->map == map)
	{ controlDeleteInstance(control, instance);
	  i--;
	}
    }

  map->viewKey = view;

  map->submenus = bsFindTag(View, _Submenus);
  map->cambridgeOptions = bsFindTag(View, _Cambridge);
  map->noButtons = bsFindTag(View, _Nobuttons);
  control->hideHeader = bsFindTag(View, _HideHeader);

  if (bsGetData(View, _Columns, _Text, &s))
    do
      { 
	mark = bsMark(View, mark);
	if  (bsGetData(View, _bsRight, _Int, &disp) && bsPushObj(View))
	  {
	    if (bsGetKeyTags(View, _bsRight, &colKey))
	      for (i = 0; i<arrayMax(map->protoArray); i++)
		{ proto = arrp(map->protoArray, i, struct ProtoStruct); 
		  if (proto->key == colKey)
		    { if (proto->unique && instanceExists(map, proto))
			{ messout("Cannot create two columns of type %s",
				  name(colKey));
			  break;
			}
		      instance = controlInstanceCreate(proto,
						       map,
						       disp,
						       View,
						       s); 
		      if (instance)
			array(map->control->instances, 
			      arrayMax(map->control->instances),
			      COLINSTANCE) = instance ;
		      break;
		    }
		}		
	    if (! proto) 
	      messout("Column prototype missing for %s", name(colKey));
	  }
	bsGoto(View, mark); 
      }
    while (bsGetData(View, _bsDown, _Text, &s));
  
  bsDestroy(View);
  bsMarkFree(mark);
  if (arrayMax(control->maps) > 1)
    control->hideHeader = FALSE; /* don't hide header on multiple maps - too
				    confusing */
  return TRUE;
  }


static KEY controlSaveView(MAPCONTROL map)
{ KEY _VView, _Columns, _Submenus, _Cambridge, _Nobuttons, _HideHeader;
  KEY protoTag;
  COLCONTROL control = map->control;
  COLINSTANCE instance; 
  char buffer[41];
  KEY key;
  int i;
  BSMARK mark = 0 ;
  OBJ obj;
  
  lexaddkey("View", &key, _VMainClasses);
  _VView = KEYKEY(key);
  lexaddkey("Columns" , &_Columns, 0);
  lexaddkey("Submenus", &_Submenus, 0);
  lexaddkey("Cambridge", &_Cambridge, 0);
  lexaddkey("No_buttons", &_Nobuttons, 0);
  lexaddkey("Hide_header", &_HideHeader, 0);

  if (!checkWriteAccess())
      return 0;

  if (!messPrompt(messprintf("Save config of %s as:", map->name),
		  map->viewKey ? name(map->viewKey) : "", "t"))
    return 0;

  strncpy(buffer, freeword(), 40);
  if (strlen(buffer) == 0)
    return 0;

  if (lexword2key(buffer, &key, _VView) && 
      !messQuery(messprintf("Overwrite existing %s?", buffer)))
    return 0;

  
  lexaddkey(buffer, &key, _VView);
  
  if (!key || !(obj = bsUpdate(key)))
    return 0;

  map->viewKey = key;

  bsKill(obj) ; obj = bsUpdate (key) ; /* kill off old version */

  bsAddTag(obj, map->viewTag);

  if (map->submenus)
    bsAddTag(obj, _Submenus);
  if (map->cambridgeOptions)
    bsAddTag(obj, _Cambridge);
  if (map->noButtons)
    bsAddTag(obj, _Nobuttons);
  if (control->hideHeader)
    bsAddTag(obj, _HideHeader);

  for (i=0; i<arrayMax(control->instances); i++)
    { instance = arr(control->instances, i, COLINSTANCE);
      if (instance->map != map)
	continue ;
      mark = bsMark(obj, mark) ;
      if (!bsAddData (obj, _Columns, _Text, instance->name) ||
	  !bsAddData (obj, _bsRight, _Int, &instance->displayed) ||
	  !bsPushObj (obj))
	{ messerror ("Invalid View model - need Columns Text Int #Column") ;
	  break ;
	}
      lexaddkey (instance->proto->name, &protoTag, 0) ;
      if (!bsAddTag (obj, protoTag))
	{ messerror ("Column type %s not a tag in ?Column -- can't store",
		     instance->proto->name) ;
	  bsGoto (obj, 0) ;
	  bsAddData (obj, _Columns, _Text, instance->name) ;
	  bsRemove (obj) ;
	  continue ;
	}
      if (instance->save)
	(*instance->save)(instance, obj) ;
      bsGoto(obj, mark) ;
    }
  bsSave (obj) ;

  bsMarkFree (mark) ;

  return key;
}


int controlGetColour(OBJ obj)
{ BSMARK mark = bsMark(obj, 0);
  KEY col;
  int result = WHITE;

  if (bsPushObj(obj) && bsGetKeyTags(obj, _bsRight, &col))
    result = col - _WHITE;
  
  bsGoto(obj, mark);
  bsMarkFree(mark);
  return result;
}

void controlSetColour(OBJ obj, int colour)      
{ BSMARK mark = bsMark(obj, 0);
  if (bsPushObj(obj))
    bsAddTag(obj, _WHITE + colour);
  bsGoto(obj, mark);

  bsMarkFree(mark);
}



static Array addToMenu(Array menu, AC_HANDLE handle, KEY view)
{ char *viewnam;
  KEY _Name;
  int i, j;
  OBJ v;

  lexaddkey("Name", &_Name, 0);
  
  if (menu)
    { for (i = 0; i<arrayMax(menu); i++) /* ignore duplicates */
	if ( array(menu, i, FREEOPT).key == view)
	  return menu;
    }
  else
    { menu = arrayHandleCreate(10, FREEOPT, handle);
      array(menu, 0, FREEOPT).text = "Views";
    }
  
  if (!view) /* view = 0 means call control window */
    viewnam = "View control";
  else
    { char *s;
      if (!(v = bsCreate(view)) || !bsGetData(v, _Name, _Text, &s))
	s = name(view); /* default */
      viewnam = handleAlloc(0, handle, 1+strlen(s));
      strcpy(viewnam, s);
      bsDestroy(v); 
    }

  j = arrayMax(menu);
  array(menu, j, FREEOPT).text = viewnam;
  array(menu, j, FREEOPT).key = view;
  array(menu, 0, FREEOPT).key = j;
  
  return menu;
}





void controlMakeColumns(MAPCONTROL map, COLDEFAULT colInit, KEY override, BOOL needMinimal)
{ 
  /* Add a map to the current window.
     If a view is avialable, then use it to display columns, else 
     colInitString says which columns should be instantiated,
     a view supplied as arg takes priority over any found in the map object.
     Columns are added at the RHS of any existing in the window */

  OBJ o ;
  KEY view;
  KEY def = 0;
  KEY min = 0;
  KEY _View, _Default_view, _Minimal_view, _From_map;
  Array menu = 0;
  KEY key = map->key;


  
  lexaddkey("View", &_View, 0);
  lexaddkey("Default_view", &_Default_view, 0);
  lexaddkey("Minimal_view", &_Minimal_view, 0);
  lexaddkey("From_map", &_From_map, 0);


  menu = addToMenu(menu, map->handle, 0); /* call vieweditor, for mac */
  
  while (key && (o = bsCreate(key)))
    { 
      if (bsGetKey(o, _Default_view, &view))/* look for Default_view */
	{ if (!def) /* first default we find, we use */
	    def = view;
	  menu = addToMenu(menu, map->handle, view);
	}
      
      if (bsGetKey(o, _Minimal_view, &view))/* look for Default_view */
	{ if (!min) /* first default we find, we use */
	    min = view;
	  menu = addToMenu(menu, map->handle, view);
	}
      
      if (bsGetKey(o, _View, &view))
	do { 
	  menu = addToMenu(menu, map->handle, view);
	}  while (bsGetKey(o, _bsDown, &view));
    
      if (!bsGetKey(o, _From_map, &key)) /* do any we inherit from */
	key = 0;
      bsDestroy(o);
    }

  if (menu)
    map->viewMenu = arrayp(menu, 0, FREEOPT);

  if (!override)
    { if (needMinimal && min)
	override = min;
      else
	override = def;
    }
  if (!override  || !controlSetView(map, override)) 
    /* if there's a default view use it */
    controlReadConfig(map, colInit);

}

static void viewMenuFunc(KEY key, int box)
{
  COLCONTROL control = currentColControl("viewMenuFunc");

  if (key == 0)
    viewWindowCreate(0);
  else
    {
      controlSetView(control->currentMap, key);
      if (control->viewWindow) 
	{
	  /* destroy all memory associated with the viewWindow */
	  handleDestroy(control->viewWindow->handle);
	  /* kill the viewWindow struct and cause viewWindowFinalise() */
	  messfree (control->viewWindow);/* destroys the window */
	}
      controlDraw();
    }

  return;
} /* viewMenuFunc */



static void mapLocatorSave(COLINSTANCE instance, OBJ obj)
{ KEY _Magnification, _Projection_lines_on;

  lexaddkey("Magnification", &_Magnification, 0);
  lexaddkey("Projection_lines_on", &_Projection_lines_on, 0);
  
  if (instance->map->magConf != 0.0)
    bsAddData(obj, _Magnification, _Float, &instance->map->mag);
  if (instance->map->hasProjectionLines)
    bsAddTag(obj, _Projection_lines_on);
}




static void mapScaleSave(COLINSTANCE instance, OBJ obj)
{  KEY _Cursor_unit, _Cursor_on, _Scale_unit;
   MAPCONTROL map = instance->map;

   lexaddkey("Cursor_unit", &_Cursor_unit, 0);
   lexaddkey("Cursor_on", &_Cursor_on, 0);
   lexaddkey("Scale_unit", &_Scale_unit, 0);
   
   bsAddData(obj, _Cursor_unit, _Float, &map->cursor.unit);
   bsAddData(obj, _Scale_unit, _Float, &map->scaleUnit);
   if (map->hasCursor)
     bsAddTag(obj, _Cursor_on);
}


static void spacerSave(COLINSTANCE instance, OBJ obj)
{ KEY _Spacer_width, _Spacer_colour;
  SPACERPRIV private = instance->private;

  lexaddkey("Spacer_Width", &_Spacer_width, 0);
  lexaddkey("Spacer_colour", &_Spacer_colour, 0);

  bsAddData(obj, _Spacer_width, _Float, &private->width);
  bsAddTag(obj, _Spacer_colour);
  controlSetColour(obj, private->colour);
}

