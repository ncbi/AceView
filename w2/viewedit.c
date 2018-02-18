/*  File:viewedit.c
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 * Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: view object editor.
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  4 11:59 1999 (fw)
 * * Oct 13 14:40 1998 (edgrif): Replace ACEDB defs with function calls.
 * Created: Wed 1 Mar 1995 (srk)
 *-------------------------------------------------------------------
 */

/*  $Id: viewedit.c,v 1.2 2003/03/17 05:03:42 mieg Exp $  */

#include "graph_.h"
#include "../wh/menu.h"
#include "colcontrol_.h"

/************************************************************/

static void viewWindowFinalise(void *p);
static void viewWindowPick(int box);
static MENUSPEC *makeViewMenu() ;

/************************************************************/

static magic_t GRAPH2VIEWWINDOW_ASSOC = "VIEWWINDOW";
static magic_t VIEWWINDOW_MAGIC = "VIEWWINDOW";

/************************************************************/


VIEWWINDOW currentView (char *caller)
{
  VIEWWINDOW view ;

  if (!(graphAssFind(&GRAPH2VIEWWINDOW_ASSOC, &view)))
    messcrash("%s() could not find VIEWWINDOW on graph", caller);
  if (!view)
    messcrash("%s() received NULL VIEWWINDOW pointer", caller);
  if (view->magic != &VIEWWINDOW_MAGIC)
    messcrash("%s() received non-magic VIEWWINDOW pointer", caller);
  
  return view;
} /* currentView */


static void viewWindowFinalise (void *p)
{
  VIEWWINDOW view = (VIEWWINDOW)p;
  COLCONTROL control = view->control;
  Graph old = graphActive();
  
  graphActivate (view->graph) ;
  control->viewWindow = 0 ;

  /* ACEDB-GRAPH INTERFACE: acedb needs to mark altered but not saved views. */
  if (getGraphAcedbViewAnon() != NULL) (getGraphAcedbViewAnon())(view, control) ;

  graphDestroy () ;
  graphActivate (old) ;
  
  return ;
} /* viewWindowFinalise */


static void viewWindowQuit(void)
{ 
  VIEWWINDOW view = currentView("viewWindowQuit");

  controlDrawControl(view->control);
  handleDestroy (view->handle);	/* clear storage inside view struct */
  messfree (view);		/* will cause viewWindowFinalise */

  return;
} /* viewWindowQuit */


static void viewWindowMenuQuit(MENUITEM m)
{
  viewWindowQuit();

  return;
} /* viewWindowMenuQuit */


static void viewDeleteInstance(MENUITEM m)
{
  VIEWWINDOW view = currentView("viewDeleteInstance");

  if (view->currentInstance)
    { controlDeleteInstance(view->control, view->currentInstance);
      view->currentInstance = 0;
      view->editInstance = FALSE;
      view->dirty = TRUE;
      viewWindowDraw(view);
    }

  return;
} /* viewDeleteInstance */


static void viewConfigInstance(MENUITEM m)
{ 
  VIEWWINDOW view = currentView("viewConfigInstance");

  if (view->currentInstance && view->currentInstance->configure)
    if ((*view->currentInstance->configure)(view->currentInstance))
      { view->dirty = TRUE;
	viewWindowDraw(view);
      }

}
  
/* Instance names must be unique in map, to save properly */

static char *nameDisambiguate(MAPCONTROL map, char *namein)
{ COLCONTROL control = map->control;
  COLINSTANCE instance;
  static char buf[200];
  int i;
  char *cp;
  
  if (namein != buf)
    strcpy(buf, namein);
  
  for (cp = buf; *cp; cp++)
    if (*cp == ' ')
      *cp = '_';
  
  for (i = 0; i<arrayMax(control->instances); i++)
    { instance = arr(control->instances, i, COLINSTANCE);
      if (instance->map == map &&
	  strcmp(instance->name, buf) == 0)
	{
	  strcat(buf, "+");
	  return nameDisambiguate(map, buf);
	}
    }
	
  return buf;  
}




static void helpOnColumn(MENUITEM m)
{
  VIEWWINDOW view = currentView("helpOnColumn");
  char *help;

  if (view->currentInstance)
    help = view->currentInstance->proto->helpText;
  else if (view->currentProto)
    help = view->currentProto->helpText;
  else
    return;
  
  if (!help)
    help = "Sorry, no help available here.\n";
 
  graphUnMessage();
  graphMessage(help);
}

void editDone(char *buff)
{ VIEWWINDOW view = currentView("editDone");
  COLINSTANCE instance = view->currentInstance;
  
  view->editInstance = FALSE;
  if (strcmp(instance->name, buff) != 0)
    { char *name = nameDisambiguate(instance->map, buff);
      view->dirty = TRUE;
      messfree(instance->name);
      instance->name = strnew (name, instance->handle);
    }

  viewWindowDraw(view);
}


static void viewEditMenuFunc(KEY key, int box)
{
  VIEWWINDOW view = currentView("viewEditMenuFunc");


  /* ACEDB-GRAPH INTERFACE: acedb has its own colcontrol routine.            */
  if (getGraphAcedbGetSetView() != NULL) (getGraphAcedbGetSetView())(view->currentMap, key);


  controlDrawControl(view->control);
  view->dirty = FALSE;
  view->currentInstance = 0;
  view->currentProto = 0;
  view->editInstance = FALSE;
  viewWindowDraw(view);
}


static void nonFunc(void)
{ return;
}


static void viewSelectMap(KEY key, int box)
{
  VIEWWINDOW view = currentView("viewSelectMap");
  COLCONTROL control = view->control;

  view->currentMap = arr(control->maps, (int)key, MAPCONTROL);
  view->currentProto = 0;
  view->currentInstance = 0;
  view->editInstance = FALSE;
  viewWindowDraw(view);
}

void viewWindowDraw(VIEWWINDOW view)
{ int box, i;
  COLCONTROL control = view->control;
  COLINSTANCE instance;
  MAPCONTROL map = view->currentMap;
  COLPROTO proto;
  BOOL multi = arrayMax(control->maps)>1;
  MENU menu = view->menu;
  float iOff = 35;
  float y = 5;
  BOOL result ;

  if (graphActivate(view->graph))
    graphPop();
  else
    return;

  graphClear();

  if (view->currentInstance && view->currentInstance->configure)
    menuUnsetFlags(menuItem(menu, "Configure column"), MENUFLAG_DISABLED);
  else
    menuSetFlags(menuItem(menu, "Configure column"), MENUFLAG_DISABLED);

  if (view->currentInstance)
    menuUnsetFlags(menuItem(menu, "Delete column"), MENUFLAG_DISABLED);
  else
    menuSetFlags(menuItem(menu, "Delete column"), MENUFLAG_DISABLED);

  if (view->currentInstance || view->currentProto)
    menuUnsetFlags(menuItem(menu, "Help on column"), MENUFLAG_DISABLED);
  else
    menuSetFlags(menuItem(menu, "Help on column"), MENUFLAG_DISABLED);


  /* ACEDB-GRAPH INTERFACE: acedb has its own save menu items, if they are   */
  /* registered use them, otherwise use graph defaults.                      */
  result = FALSE ;
  if (getGraphAcedbSaveView() != NULL) result = (getGraphAcedbSaveView())(menu) ;
  if (result == FALSE)
    {
    menuSetFlags(menuItem(menu, "Save view"), MENUFLAG_DISABLED);
    menuSetFlags(menuItem(menu, "Save view as default"), MENUFLAG_DISABLED);
    }



  graphNewMenu(view->menu);

  if (multi)
    { box = graphButton("Maps...", nonFunc, iOff+5, 1);
      graphBoxFreeMenu(box, viewSelectMap, view->mapMenu); 
    }
  else
    { box = graphButton("Views..", nonFunc, iOff+5, 1);
      if (map->viewMenu)
	graphBoxFreeMenu(box, viewEditMenuFunc, map->viewMenu);
    }

  graphColouredButton("Re-draw", controlDrawControl, control,
		      BLACK, WHITE, iOff+18, 1);

  graphButton("Quit", viewWindowQuit, iOff+30, 1);
  
  assClear(view->boxToArrow);
  assClear(view->boxToInstance);
  assClear(view->boxToProto);
  assClear(view->boxToButton);

  graphLine(iOff-2, 0, iOff-2, 200);
  graphLine(iOff-2, 4, 100, 4);
  graphLine(0, 7, iOff-2, 7);


  if (multi)
    { box = graphBoxStart();
      graphText(messprintf("Current map: %s", map->name), 1, 1);
      graphBoxEnd();
      graphBoxDraw(box, BLACK, map->colour);
    }
  else
    {
    /* ACEDB-GRAPH INTERFACE: acedb displays its own save view messages.     */
    if (getGraphAcedbSaveViewMsg()) (getGraphAcedbSaveViewMsg())(map, view) ;
    }

  
  graphText(messprintf("Column type: %s", view->currentInstance ?  
		       view->currentInstance->proto->name : ""),
	    iOff, y++);
  if (multi)
    { box = graphBoxStart();
      graphText(messprintf("Map: %s", view->currentInstance ?
			 view->currentInstance->map->name : ""),
	      iOff, y++);
      graphBoxEnd();
      if (view->currentInstance)
	graphBoxDraw(box, BLACK, view->currentInstance->map->colour);
    }

  view->submenusBox = graphBoxStart();
  graphArc(2, 3.5, 0.7, 0, 360);
  if (view->currentMap->submenus)
    graphFillArc(2, 3.5, 0.4, 0, 360);
  graphBoxEnd(); 
  graphBoxDraw(view->submenusBox, BLACK, TRANSPARENT);
  graphText("Use submenus", 3.2, 3);

  view->cambridgeBox = graphBoxStart();
  graphArc(2, 4.5, 0.7, 0, 360);
  if (view->currentMap->cambridgeOptions)
    graphFillArc(2, 4.5, 0.4, 0, 360);
  graphBoxEnd(); 
  graphBoxDraw(view->cambridgeBox, BLACK, TRANSPARENT);
  graphText("Cambridge options", 3.2, 4);

  view->hideBox = graphBoxStart();
  graphArc(2, 5.5, 0.7, 0, 360);
  if (view->control->hideHeader)
    graphFillArc(2, 5.5, 0.4, 0, 360);
  graphBoxEnd(); 
  graphBoxDraw(view->hideBox, BLACK, TRANSPARENT);
  graphText("Hide headers", 3.2, 5);
 
  view->suppressBox = 0;
  view->hideBox = 0;
  if (multi)
    { view->suppressBox = graphBoxStart();
      graphArc(2, 5.5, 0.7, 0, 360);
      if (view->currentMap->suppressed)
	graphFillArc(2, 5.5, 0.4, 0, 360);
      graphBoxEnd(); 
      graphBoxDraw(view->suppressBox, BLACK, TRANSPARENT);
      graphText("Suppress map", 3.2, 5);
    }
  else
    { view->hideBox = graphBoxStart();
      graphArc(2, 5.5, 0.7, 0, 360);
      if (view->control->hideHeader)
	graphFillArc(2, 5.5, 0.4, 0, 360);
      graphBoxEnd(); 
      graphBoxDraw(view->hideBox, BLACK, TRANSPARENT);
      graphText("Hide headers", 3.2, 5);
    }

  box = graphBoxStart();
  graphLine(iOff+1, y+1.0, iOff+2, y+1.0);
  graphLine(iOff+1.5, y+0.75, iOff+2, y+1.0);
  graphLine(iOff+1.5, y+1.25, iOff+2, y+1.0);
  graphBoxEnd();
  y++;

  assInsert(view->boxToArrow, assVoid(box), assVoid(0));

  for(i=0; i<arrayMax(control->instances); i++)
    { 
      instance = array(control->instances, i,  COLINSTANCE);
      if (!instance->map->suppressed)
	{ 
	  box = graphBoxStart();
	  graphArc(iOff+3.3, y+0.5, 0.7, 0, 360);
	  if (instance->displayed)
	    graphFillArc(iOff+3.3, y+0.5, 0.4, 0, 360);
	  graphBoxEnd(); 
	  graphBoxDraw(box, BLACK, TRANSPARENT);
	  assInsert(view->boxToButton, assVoid(box), assVoid(i));

	  box = graphBoxStart();
	  graphLine(iOff+1, y+1.0, iOff+2, y+1.0);
	  graphLine(iOff+1.5, y+0.75, iOff+2, y+1.0);
	  graphLine(iOff+1.5, y+1.25, iOff+2, y+1.0);
	  graphBoxEnd();
	  assInsert(view->boxToArrow, assVoid(box), assVoid(i+1));
	  
	  
	  if (instance == view->currentInstance && view->editInstance)
	    {
	      strcpy(view->buffer, instance->name);
	      box = graphTextEntry(view->buffer, EDITLEN, 
				   iOff+4.5, y, editDone);
	    }
	  else
	    {
	      box = graphBoxStart();
	      graphText(instance->name, iOff+4.5, y);
	      graphBoxEnd();
	      if (instance == view->currentInstance)
		graphBoxDraw(box, WHITE, BLACK);
	      else if (multi)
		graphBoxDraw(box, BLACK, instance->map->colour);
	    }
	  assInsert(view->boxToInstance, assVoid(box), assVoid(i));
	  y++;
	}
    }


  if (!map->suppressed)
    {
      graphText("Available Column Types", 4, 8);
      
      for(y=10, i=0; i<arrayMax(map->protoArray); i++, y++)
	{
	  proto = arrp(map->protoArray, i, struct ProtoStruct);
	  box = graphBoxStart();
	  if (proto->unique && instanceExists(map, proto))
	    graphTextFormat(ITALIC); /* Can't have another one */
	  graphText(proto->name, 3, y);
	  graphTextFormat(PLAIN_FORMAT);
	  graphBoxEnd();
	  assInsert(view->boxToProto, assVoid(box), proto);
	  if (proto == view->currentProto)
	    graphBoxDraw(box, WHITE, BLACK); 
	}
    }

  graphRedraw();
}

static void viewWindowPick(int box)
{
  VIEWWINDOW view = currentView("viewWindowPick") ;
  COLCONTROL control = view->control ;
  COLINSTANCE instance;
  COLPROTO proto;
  int i, insertBefore;
  void *b;

  if (box == 0) return; 

  if (box == view->submenusBox)
    { view->currentMap->submenus = !view->currentMap->submenus;
      view->dirty = TRUE;
    }
  else if (box == view->cambridgeBox)
    { view->currentMap->cambridgeOptions = !view->currentMap->cambridgeOptions;
      view->dirty = TRUE;
    }
  else if (box == view->hideBox)
    { view->control->hideHeader = !view->control->hideHeader;
      view->dirty = TRUE;
    }
  else if (box == view->suppressBox)
    { view->currentMap->suppressed = !view->currentMap->suppressed;
      view->dirty = TRUE;
    }
  else if (assFind(view->boxToButton, assVoid(box), &b))
    { instance = array(control->instances, assInt(b), COLINSTANCE);
      instance->displayed = !instance->displayed;
      view->dirty = TRUE;
    }
  else if (assFind(view->boxToProto, assVoid(box), &proto))
    { /* if can't be picked, ignore */
      if (!proto->unique || !instanceExists(view->currentMap, proto))
	{
	  view->currentInstance = 0;
	  view->editInstance = FALSE;
	  view->currentProto = proto;
	}
    }
  else if (assFind(view->boxToInstance, assVoid(box), &b))
    { 
      if (view->currentInstance == array(control->instances, assInt(b),
					 COLINSTANCE))
	{ view->editInstance = TRUE;
	}
      else
	{ view->currentInstance = array(control->instances, assInt(b),
					COLINSTANCE);
	  view->currentInstanceIndex = assInt(b);
	  view->currentProto = 0;
	  view->editInstance = FALSE;
	}
    }
  else if (assFind(view->boxToArrow, assVoid(box), &insertBefore) 
      && (view->currentInstance || view->currentProto))
    { /* Pressed an Arrow, if the current box is an instance, move it here,
	 else make an instance from the prototype an put that here */

      if (view->currentInstance)
	{ /* current is an instance, pull it out and close up the space */
	  int i;
	  instance = view->currentInstance;
	  if (insertBefore > view->currentInstanceIndex)
	    insertBefore--; /* co-ords change */
	  for (i=view->currentInstanceIndex; 
	       i<arrayMax(control->instances)-1;
	       i++)
	    arr(control->instances, i, COLINSTANCE) = 
	      arr(control->instances, i+1, COLINSTANCE);
	  arrayMax(control->instances) = arrayMax(control->instances) - 1;
	}
      else 
	{ /* current is a prototype, make an instance */
	  instance = 
	    controlInstanceCreate(view->currentProto,
				  view->currentMap,
				  TRUE,
				  0,
				  nameDisambiguate(view->currentMap,
						   view->currentProto->name)
				  );
	  view->currentProto = 0;
	}
      /* bump up the others to make room for the new one. */
      /* then insert it */
      if (instance)
	{ for (i=arrayMax(control->instances)-1; i >= insertBefore; i--)
	    { array(control->instances, i+1, COLINSTANCE) =
		array(control->instances, i, COLINSTANCE);
	    }
	  view->currentInstanceIndex = insertBefore;
	  array(control->instances, insertBefore, COLINSTANCE) = instance ;
	  view->currentInstance = instance;
	  view->editInstance = FALSE;
	  view->dirty = TRUE;
	}
    }
  
  viewWindowDraw(view);
  return;
  
}



/* This function constructs the viewedit menu of functions.                  */
/*                                                                           */
static MENUSPEC *makeViewMenu()
{
  static MENUSPEC menu_items[15] ;			    /* Must be static. */
  MENUSPEC *result, *curr ;

  result = curr = &menu_items[0] ;

  curr->f = viewWindowMenuQuit ;			    /* First item always 'quit' */
  curr->text = "Quit" ;

  /* ACEDB-GRAPH INTERFACE: acedb needs to add some functions to the view    */
  /* menu, call the function to do this if registered.                       */
  /* This call will use the current pointer and leave it pointing to its     */
  /* last entry.                                                             */
  if (getGraphAcedbViewMenu() != NULL) (getGraphAcedbViewMenu())(&curr) ;

  curr++ ;

  /* OK, now add the rest of the standard view functions.                    */
  curr->f = (MENUFUNCTION)menuSpacer ;
  curr->text = "" ;
  curr++ ;
  curr->f = viewConfigInstance ;
  curr->text = "Configure column" ;
  curr++ ;
  curr->f = viewDeleteInstance ;
  curr->text = "Delete column" ;
  curr++ ;
  curr->f = helpOnColumn ;
  curr->text = "Help on column" ;
  curr++ ;
  curr->f = NULL ;
  curr->text = NULL ;

  return(result) ;
} /* makeViewMenu */





void viewWindowCreate(void *dummy)
{
  COLCONTROL control = currentColControl("viewWindowCreate");
  VIEWWINDOW view;
  BOOL multi = arrayMax(control->maps)>1;
  int i;
  Graph old = graphActive();
  Array mapMenu;

  /* if map window already has a control, just redraw it to pop it */
  if (control->viewWindow) 
    { viewWindowDraw(control->viewWindow);
      graphActivate(old);
      return;
    }

  view = (VIEWWINDOW)halloc(sizeof(struct ViewWindowStruct),
			    control->handle);
  blockSetFinalise (view, viewWindowFinalise);
  /* NB, this causes proper cleanup and destruction of window purely
     by freeing the parents handle. To kill the view-window and all
     its associated memory, a messfree is sufficient */

  view->magic = &VIEWWINDOW_MAGIC;
  view->handle = handleHandleCreate(control->handle);


  mapMenu = arrayHandleCreate(10, FREEOPT, view->handle);
  array(mapMenu, 0, FREEOPT).text = "ViewEdit maps";
  array(mapMenu, 0, FREEOPT).key = arrayMax(control->maps);
  view->mapMenu = arrayp(mapMenu, 0, FREEOPT);
  for (i=0; i<arrayMax(control->maps); i++)
    { array(mapMenu, i+1, FREEOPT).text =
	arr(control->maps, i, MAPCONTROL)->name;
      array(mapMenu, i+1, FREEOPT).key = i; /* abuse key as int */
    }

  view->menu = menuInitialise("ViewEdit menu", makeViewMenu());
  menuSetFlags(menuItem(view->menu, "Configure column"), MENUFLAG_DISABLED);
  menuSetFlags(menuItem(view->menu, "Delete column"), MENUFLAG_DISABLED);
  if (multi)
    { menuSuppress(view->menu, "Save view");
      menuSuppress(view->menu, "Save view as default");
    }
  
 
  view->currentMap = array(control->maps, 0, MAPCONTROL);
  view->control = control;


  /* ACEDB-GRAPH INTERFACE: if acedb has registered a different display      */
  /* routine then call it, otherwise use the default.                        */
  if (getGraphAcedbDisplayCreate() != NULL)
    view->graph = (getGraphAcedbDisplayCreate())(getGraphAcedbViewName()) ;
  else
    view->graph = graphCreate(TEXT_FIT, "Display Chooser", 0.2, 0.2, 0.8, 0.5);


  view->boxToArrow = assHandleCreate(view->handle);
  view->boxToInstance = assHandleCreate(view->handle);
  view->boxToProto = assHandleCreate(view->handle);
  view->boxToButton = assHandleCreate(view->handle);
  view->currentInstance = 0;
  view->currentProto = 0;
  view->editInstance = FALSE;
  view->dirty = FALSE;
  graphAssociate (&GRAPH2VIEWWINDOW_ASSOC, view) ;
  graphRegister (PICK, viewWindowPick) ;

  control->viewWindow = view;
  
  viewWindowDraw(view);

  graphActivate(old); /* otherwise the callback code gets confused */
} /* viewWindowCreate */


/************************** eof *****************************/
