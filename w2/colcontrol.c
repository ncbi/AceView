/*  File:colcontrol.c
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 * Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: package for generic map drawing containing:
 		columns cotrol
		zoom, middle button scroll, locator
		findScaleUnits
		multiple Map support
		column data abstraction
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  4 11:36 1999 (fw)
 * * Oct 13 14:30 1998 (edgrif): Replaced ACEDB defines with function
 *              calls. Removed some acedb specific code to acedbgraph.c
 * Created: Wed Mar 20 13:40:35 1992 (srk)
 *-------------------------------------------------------------------
 */

/*  $Id: colcontrol.c,v 1.4 2007/12/10 19:26:34 mieg Exp $ */

#include "graph_.h"

#include "colcontrol_.h"

/************************************************************/

magic_t GRAPH2COLCONTROL_ASSOC = "COLCONTROL";
magic_t GRAPH2COLCONTROL_LOCALS_ASSOC = "COLCONTROL_LOCALS";

static magic_t COLCONTROL_MAGIC = "COLCONTROL";
static magic_t MAPCONTROL_MAGIC = "MAPCONTROL";

/************************************************************/

static void controlDestroyed(void);
static void controlLeftDrag(double x, double y);
static void controlLeftUp(double x, double y);
static void controlPick(int box, double x, double y);
static void controlSelect(int box, double x, double y);
static void controlKybd(int k);


static void controlMiddleDrag(double x, double y);
static void controlMiddleUp(double x, double y);
static void controlMiddleDown(double x, double y);

/* Used to assign colours to maps */
static int colourOrder[] = { GREEN, YELLOW, CYAN, MAGENTA, LIGHTGRAY };

static COLOUROPT ColButton = { viewWindowCreate, 0, BLACK, WHITE, "Views...", 0 };

/************************************************************/


/* redo the configure routine to be non blocking and use RD's editor routines */
static void controlConfigOK(void)
{
  COLINSTANCE instance;
  void *locals;

  if(graphCheckEditors(graphActive(),TRUE))
    { 
      graphAssFind(&GRAPH2COLCONTROL_ASSOC, &instance);
      graphAssFind(&GRAPH2COLCONTROL_LOCALS_ASSOC, &locals); 
      (*(instance->configFinal))(instance, locals, TRUE); /* FALSE in KO case */
      
      instance->configGraph = 0;
      graphDestroy();
      
      controlDrawControl(instance->map->control);
    } 
} /* controlConfigOK */

static void controlConfigKO(void)
{
  COLINSTANCE instance;
  void *locals;

  if(!graphAssFind(&GRAPH2COLCONTROL_ASSOC, &instance))
    printf("Could not find instance associated with GRAPH2COLCONTROL_ASSOC\n");
  if(!graphAssFind(&GRAPH2COLCONTROL_LOCALS_ASSOC, &locals))
    printf("Could not find locals associated with GRAPH2COLCONTROL_LOCALS_ASSOC\n"); 
  (*(instance->configFinal))(instance, locals, FALSE); /* FALSE in KO case */

  instance->configGraph = 0;
  graphDestroy();
} /* controlConfigKO */


static void controlConfigApply(void)
{
  COLINSTANCE instance;
  void *locals;

  if(!graphAssFind(&GRAPH2COLCONTROL_ASSOC, &instance))
    printf("Could not find instance associated with GRAPH2COLCONTROL_ASSOC\n");
  if(!graphAssFind(&GRAPH2COLCONTROL_LOCALS_ASSOC, &locals))
    printf("Could not find locals associated with GRAPH2COLCONTROL_LOCALS_ASSOC\n"); 
  (*(instance->configFinal))(instance, locals, TRUE); /* FALSE in KO case */

  controlDrawControl(instance->map->control);

/*  instance->configGraph = 0;
  graphDestroy();*/
} /* controlConfigApply */


Graph controlCreateConfig(COLINSTANCE instance, void *locals,char *text, float width, float height) 
{
  Graph g;
  int width2,height2;

  if (instance->configGraph)
    {
      graphActivate(instance->configGraph);
      graphPop();
      return 0;
    }
  
  g = graphCreate(TEXT_SCROLL,text,0.5,0.4, width, height); 

  instance->configGraph = g;

  graphAssociate(&GRAPH2COLCONTROL_ASSOC, instance);
  graphAssociate(&GRAPH2COLCONTROL_LOCALS_ASSOC, locals);

  graphRegister(DESTROY,controlConfigKO);

  graphFitBounds(&width2,&height2);
  graphButton("Apply",controlConfigApply,(float)(width2/2.0)-7.0,height2 - 2.0);
  graphButton("OK",controlConfigOK,(float)(width2/2.0),height2 - 2.0);
  graphButton("Cancel",graphDestroy,(float)(width2/2.0)+4.0,height2 - 2.0);

  return g;
}



/* There used to be on colControlCreate which had different function         */
/* signatures and code according to whether it was compiled with ACEDB or    */
/* not. I have now split this into two functions with some duplicated code   */
/* which will be gradually be unified....                                    */
/*                                                                           */
/* This is the form for ACEDB type code that is going to create its own      */
/* display.                                                                  */
/*                                                                           */
/* NOTE that this routine is not going to work well unless the application   */
/* has registered a displayCreate routine to set up the graph, the default   */
/* action is to attempt to display something.                                */
/*                                                                           */
COLCONTROL colControlCreate(BOOL isoldgraph, char *name, char *dispName)
{
  int gw, gh;
  AC_HANDLE handle = handleCreate();
  COLCONTROL control = (COLCONTROL)halloc(sizeof (struct ColControlStruct), handle);

  if (isoldgraph) 
    {
    controlDestroyed();
    graphClear();
    control->graph = graphActive();
    }
  else 
    {
    /* ACEDB-GRAPH INTERFACE: if display function registered, call it,       */
    /* otherwise do our own, note that if no function is registered then an  */
    /* educated guess is made about window size.                             */
    if (getGraphAcedbDisplayCreate() != NULL)
      control->graph = (getGraphAcedbDisplayCreate())(dispName) ;
    else
      control->graph = graphCreate(TEXT_FIT, name, 0.5, 0.5, 0, 0);

    if (!control->graph)
      {
      handleDestroy(handle);
      return 0;
      }
    graphRegister(RESIZE, controlDraw);
    graphRegister(DESTROY, controlDestroyed);
    graphRegister(KEYBOARD, controlKybd);
    graphRegister(PICK, controlPick);
    graphRegister(MIDDLE_DOWN, controlMiddleDown);
    }
  
  graphRetitle(name);					    /* Different code */

  graphAssociate(&GRAPH2COLCONTROL_ASSOC, control);
  /* link this map to the window */
      
  control->magic = &COLCONTROL_MAGIC ;
  control->handle = handle;
  control->instances = arrayHandleCreate (50, COLINSTANCE, handle);
  control->boxIndex = arrayHandleCreate(100, COLINSTANCE, handle);
  control->boxIndex2 = arrayHandleCreate(100, void *, handle);
  control->maps = arrayHandleCreate (4, MAPCONTROL, handle) ;
  control->mapTransitions = arrayHandleCreate(100, float, handle);
  control->bottomBoxes = assHandleCreate(handle);
  control->currentMap = 0;
  control->activeBox = 0;

  control->activeKey = 0;				    /* Different code */

  control->activeInstance = 0;
  control->viewWindow = 0;
  control->hideHeader = FALSE;
  graphFitBounds(&gw, &gh); /* returns ints */
  control->graphWidth = gw; /* make floats */
  control->graphHeight = gh;

  return control ;
  }


/* This routine is essentially the same as the one above with a few          */
/* differences - labelled in the routine above.                              */
/* Only coltest.c uses this form of the call within the entire CVS_ACEDB*/
/* tree....if someone outside uses it they can simply change their call */
/* to colControlBasicCreate and all should work.                        */
/*                                                                           */
COLCONTROL colControlBasicCreate(BOOL isoldgraph, int type, char *name, 
				 float x, float y, float w, float h)
{
  int gw, gh;
  AC_HANDLE handle = handleCreate();
  COLCONTROL control = (COLCONTROL)halloc(sizeof (struct ColControlStruct), handle);

  if (isoldgraph) 
    {
      controlDestroyed();
      graphClear();
      control->graph = graphActive();
    }
  else 
    { 
      control->graph = graphCreate(type, name, x, y, w, h);
      
      if (!control->graph)
	{
	  handleDestroy(handle);
	  return 0;
	}
      graphRegister(RESIZE, controlDraw);
      graphRegister(DESTROY, controlDestroyed);
      graphRegister(KEYBOARD, controlKybd);
      graphRegister(PICK, controlPick);
      graphRegister(MIDDLE_DOWN, controlMiddleDown);
    }
  
  graphAssociate(&GRAPH2COLCONTROL_ASSOC, control);
  /* link this map to the window */
      
  control->magic = &COLCONTROL_MAGIC ;
  control->handle = handle;
  control->instances = arrayHandleCreate (50, COLINSTANCE, handle);
  control->boxIndex = arrayHandleCreate(100, COLINSTANCE, handle);
  control->boxIndex2 = arrayHandleCreate(100, void *, handle);
  control->maps = arrayHandleCreate (4, MAPCONTROL, handle) ;
  control->mapTransitions = arrayHandleCreate(100, float, handle);
  control->bottomBoxes = assHandleCreate(handle);
  control->currentMap = 0;
  control->activeBox = 0;

  control->activeInstance = 0;
  control->viewWindow = 0;
  control->hideHeader = FALSE;
  graphFitBounds(&gw, &gh); /* returns ints */
  control->graphWidth = gw; /* make floats */
  control->graphHeight = gh;

  return control;
  }




 

/* Do the common parts of map creation. The caller must fill in lots more    */
/* fields in the map structure before calling controlAddMap()                */
/*                                                                           */
MAPCONTROL mapControlCreate(COLCONTROL control, void (*finaliseFunc)(void *p))
{ 
  AC_HANDLE handle = handleHandleCreate(control->handle);
  MAPCONTROL map;

  map = (MAPCONTROL)halloc(sizeof(struct MapControlStruct), handle);
  blockSetFinalise (map, finaliseFunc);
  
  map->handle = handle;
  map->magic = &MAPCONTROL_MAGIC;
  map->name = ""; /* default */
  map->menu = 0; /* default */
  

  /* ACEDB-GRAPH INTERFACE: set acedb specific map field. */
  if (getGraphAcedbResetKey() != NULL) (getGraphAcedbResetKey())(&(map->viewKey)) ;


  map->beforeDraw = 0;
  map->drawHeader = 0;
  map->headerPick = 0;
  map->convertTrans = 0;
  map->keyboard = 0;
  map->topMargin = 0;
  map->permConversionRoutines = arrayHandleCreate(50, 
						  struct ConversionRecord,
						  handle);
  map->transConversionRoutines = arrayHandleCreate(50,
						   struct ConversionRecord,
						   handle);
  map->buttons = 0;
  map->buttonsAddHere = &map->buttons;
  map->menusOnButtons = arrayHandleCreate(10, MENUOPT *, handle);

  map->suppressed = FALSE;
  map->submenus = FALSE;
  map->cambridgeOptions = TRUE;
  map->noButtons = FALSE;
  map->cursor.unit = 1; /* allows calls to set cursor without a scale column */

  return map;
}



COLCONTROL currentColControl(char *caller)
{
  COLCONTROL control;

  if (!(graphAssFind(&GRAPH2COLCONTROL_ASSOC, &control)))
    messcrash("%s() could not find COLCONTROL on graph", caller);
  if (!control)
    messcrash("%s() received NULL COLCONTROL pointer", caller);
  if (control->magic != &COLCONTROL_MAGIC)
    messcrash("%s() received non-magic COLCONTROL pointer", caller);
  
  return control;
} /* currentColControl */

MAPCONTROL currentMapControl(void)
{
  MAPCONTROL map = currentColControl("currentMapControl")->currentMap;
  char *caller = "currentMapControl";

  if (!map)
    messcrash("%s() received NULL MAPCONTROL pointer", caller);
  if (map->magic != &MAPCONTROL_MAGIC)
    messcrash("%s() received non-magic MAPCONTROL pointer", caller);

  return map;
} /* currentMapControl */


void controlAddMap(COLCONTROL control, MAPCONTROL map)
{
  COLPROTO proto, *protop;
  int c,m,i;

  /* find an unused Colour */
  for (c=0; colourOrder[c] != LIGHTGRAY; c++)
    { for (m=0; m<arrayMax(control->maps); m++)
	if (arr(control->maps, m, MAPCONTROL)->colour == colourOrder[c])
	    break;
      if (m == arrayMax(control->maps))
	break;
    }
  map->colour = colourOrder[c];
 
  array(control->maps, arrayMax(control->maps), MAPCONTROL) = map;
  map->control = control;
  controlAddButton(map, &ColButton, 0);
  
  /* the prototypes have map-private data fields, so we first make a private 
     copy of the prototypes, we put this in an array for easy access. */

  /* we also put the key correponding to the name in both the original,
   and copies, to match them again in conteolSetView and controlReadConfig */
  map->protoArray = arrayHandleCreate(30, struct ProtoStruct, map->handle);
  for (protop = map->prototypes; *protop; protop++)
    {

    /* ACEDB-GRAPH INTERFACE: If acedb is registered record the name and key */
    /* of this map.                                                          */
    if (getGraphAcedbAddMap() != NULL) (getGraphAcedbAddMap())((*protop)->name, &(*protop)->key) ;


    array(map->protoArray, arrayMax(map->protoArray), struct ProtoStruct) = **protop;
    }

  for (i = 0; i<arrayMax(map->protoArray); i++)
    { proto = arrp(map->protoArray, i, struct ProtoStruct);
      /* call init routines, if they have them. */
      if (proto->init) 
	(*proto->init)(proto, map);

      /* now register any permanent conversion routines */
 
      if (proto->conRoutine)
	{ struct ConversionRecord *c;
	  int i;
	  
	  for(i=0; i<arrayMax(map->permConversionRoutines); i++)
	    { c = arrp(map->permConversionRoutines, i, 
		       struct ConversionRecord );
	      if (c->convertRoutine == proto->conRoutine && 
		  c->params == proto->convertParams) 
		goto done;
	    }
	  /* new one here */	  
	  c = arrayp(map->permConversionRoutines,
		     arrayMax(map->permConversionRoutines),
		     struct ConversionRecord);
	  c->convertRoutine = proto->conRoutine;
	  c->params = proto->convertParams;

	done:
	  proto->convertResults = &c->results;
	}
     } 
  map->needConvert = TRUE;
  if (!control->currentMap) /* first one */
    control->currentMap = map ;
  
  if (map->topMargin > control->realTopMargin) 
    control->realTopMargin = map->topMargin;
}


void controlReadConfig(MAPCONTROL map, COLDEFAULT colInit)
{
  int i= 0;
  
  while (colInit[i].proto)
    { int j;
      COLPROTO proto;
      COLINSTANCE instance = 0;
      /* have a reference to the original proto array, find
	 correponding entry in the per-map copy */
      for (j=0; j<arrayMax(map->protoArray); j++)
	{ proto = arrp(map->protoArray, j, struct ProtoStruct);
	  if (proto->key == (colInit[i].proto)->key)
	    instance = controlInstanceCreate(proto,
					     map,
					     colInit[i].displayed,
					     0,
					     colInit[i].name);
	}
      if (instance)
	array(map->control->instances, 
	      arrayMax(map->control->instances),
	      COLINSTANCE) = instance ;

      i++;
    }
  
}	      


/* delete the instance from the map */
void controlDeleteInstance(COLCONTROL control, COLINSTANCE instance)
{
  int i,j;
  Graph old = graphActive();
  
  if (control->activeInstance == instance)
    {
    control->activeInstance = 0;

    /* ACEDB-GRAPH INTERFACE:  reset acedb specific field.                   */
    if (getGraphAcedbResetKey() != NULL) (getGraphAcedbResetKey())(&(control->activeKey)) ;

    }

  for (i=0,j=0; i<arrayMax(control->instances); i++)
    { COLINSTANCE instance1 = arr(control->instances, i, COLINSTANCE);
      if (instance1 != instance)
	arr(control->instances, j++, COLINSTANCE) = instance1;
    }

  arrayMax(control->instances) = j;

  if(instance->configGraph &&
     old != instance->configGraph)
    { /* remove configure menu if it exists */
      graphActivate(instance->configGraph);
      graphDestroy();
      graphActivate(old);
    }

  handleDestroy(instance->handle);
}

/* delete a map from the display, if it's the last, destroy the window. */
void  controlDeleteMap(COLCONTROL control, MAPCONTROL map)
{ int i, j;

  if (arrayMax(control->maps) == 1)
    { controlDestroy(control);
      return;
    }

  for (i=0; i<arrayMax(control->instances); i++)
    { COLINSTANCE instance = arr(control->instances, i, COLINSTANCE);
      if (instance->map == map)
	{ controlDeleteInstance(control, instance);
	  i--; /* above deletes one */
	}
    }
  
  for (i=0, j=0; i<arrayMax(control->maps); i++)
    { MAPCONTROL map1 = arr(control->maps, i, MAPCONTROL);
      if (map1 != map)
	arr(control->maps, j++, MAPCONTROL) = map1;
    }
  arrayMax(control->maps) = j;

  if (control->currentMap == map) /* pick another current map if needed */
    control->currentMap = arr(control->maps, j-1, MAPCONTROL);

  handleDestroy(map->handle);
  
  controlDrawControl(control);
}
  

BOOL instanceExists(MAPCONTROL map, COLPROTO proto)
{ /* returns true if an instance of proto already exists on map */
  int i;
  COLCONTROL control = map->control;
  COLINSTANCE instance;

  for (i=0; i<arrayMax(control->instances); i++)
    { instance =arr(control->instances, i, COLINSTANCE);
      if (instance->proto == proto && instance->map == map)
      return TRUE;
    }

  return FALSE;
}	




COLINSTANCE controlInstanceCreate(COLPROTO proto,
				  MAPCONTROL map,
				  BOOL displayed,
				  OBJ init,
				  char *name)
{
  COLINSTANCE instance;
  AC_HANDLE handle = handleHandleCreate(map->handle);
  
  instance = (COLINSTANCE) handleAlloc(proto->destroy, 
				       handle,
				       sizeof(struct InstanceStruct));
  instance->map = map;
  instance->handle = handle;
  instance->proto = proto;
  instance->name = handleAlloc(0, handle, strlen(name)+1);
  strcpy(instance->name, name);
  instance->displayed = displayed;
  instance->drawInit = 0; /* Next few function are optional */
  instance->configure = 0; /* Init so the create routine needn't */


  /* ACEDB-GRAPH INTERFACE: set acedb specific map field.                    */
  if (getGraphAcedbResetSave() != NULL) (getGraphAcedbResetSave())(&(instance->save)) ;


  instance->conversionRegister = 0;
  instance->setSelectBox = 0;
  instance->unSelectBox = 0;
  instance->doColour = 0;
  instance->pick = 0;
	      
  if (!(*proto->create)(instance, init))
    { handleDestroy(instance->handle) ;	
      return 0;
    }

  return instance ;
}
      
/* call this to destroy graph and everything - public */
/* can also just destroy graph */
void controlDestroy(COLCONTROL control)
{ 
  if (graphActivate(control->graph))
    { 
      graphDestroy();
      control->graph = 0 ;
    }      
}

/* called by graph Destroy, and by create code when it puts an new control */
/* on an existing graph */
static void controlDestroyed(void)
{
  COLCONTROL control = currentColControl("controlDestroyed");
  Graph old = graphActive() ;
  int i;
  
  for (i=0; i<arrayMax(control->instances); i++)
    { COLINSTANCE instance = arr(control->instances, i, COLINSTANCE);
      if(instance->configGraph &&
	 instance->configGraph != old)
	{ /* remove configure menu if it exists */
	  graphActivate(instance->configGraph);
	  graphDestroy();
	  graphActivate(old);
	}      
    }
  graphAssRemove(&GRAPH2COLCONTROL_ASSOC);
  handleDestroy(control->handle);
}


BOOL controlSetColByName(MAPCONTROL map, char *name, int mode)
{ COLCONTROL control = map->control;
  int i;
  COLINSTANCE instance;

  for(i=0; i<arrayMax(control->instances); i++)
    { instance = arr(control->instances, i, COLINSTANCE);
      if (instance->map == map && !strcmp(instance->name, name))
	switch (mode)
	  { 
	  case 0: 
	    instance->displayed = FALSE;
	    return TRUE;
	    
	  case 1:
	    instance->displayed = TRUE;
	    return TRUE;
	    
	  case 2:
	    instance->displayed = !instance->displayed;
	    return TRUE;
	    
	  default: 
	    messcrash("Illegal mode in controlSetColByName");
	  }
    }

  return FALSE;
}

void controlDrawControl(void *c)
{
  COLCONTROL control = (COLCONTROL) c;
  MAPCONTROL map;
  COLINSTANCE instance;
  BOOL active;
  float offset = 0;
  float oldOffset;
  int i, gw, gh;
  int firstFromBox = 0;

  if (graphActivate(control->graph))
    graphPop();
  else
    return;

/* The pathalogical case if no columns crashes later on unless we trap here */
  if (arrayMax(control->instances) == 0)
    return;

  control->boxIndex = arrayReCreate(control->boxIndex, 100, COLINSTANCE);
  control->boxIndex2 = arrayReCreate(control->boxIndex2, 100, void *);
  assClear(control->bottomBoxes);
  control->activeBox = 0;
  control->fromBox = 0;

  controlConvertRegister(control); /* transient conversion */

  graphFitBounds(&gw, &gh);
  control->graphWidth = gw;
  control->graphHeight = gh;

  graphClear () ;

  /* draw the header */
  if (control->hideHeader)
    control->topMargin = 1;/* au lieu de 0  mhmp 23.10.97 */
  else
    { control->topMargin = (*control->currentMap->drawHeader)(control->currentMap) ;
      control->graphHeight -= 0.6; /* space for map colour bar */
    }

  control->halfGraphHeight = 0.5 * ( control->graphHeight -
					 control->topMargin );
  control->lastHeaderBox = graphLastBox();

  for (i=0; i<arrayMax(control->maps); i++)
    { map = array(control->maps, i, MAPCONTROL);
      if (map->beforeDraw) (*map->beforeDraw)(map);
      if (map->convertTrans) (*map->convertTrans)(control, map);
    }
  
  for(i=0; i<arrayMax(control->instances); i++)
    { instance = array(control->instances, i, COLINSTANCE);
      if (instance->drawInit)
	(*instance->drawInit)(instance);
    }

  active = 
    (array(control->instances, 0, COLINSTANCE)->map == control->currentMap);

  for(i=0; i<arrayMax(control->instances); i++)
    { instance = array(control->instances, i, COLINSTANCE);
      if (instance->displayed && !instance->map->suppressed)
	{
	  oldOffset = offset;
          if (!control->hideHeader && 
	      active != (instance->map != control->currentMap)) 
	    /* only when status changes */
	    {
              if ((active = (instance->map != control->currentMap)))
		graphColor(PALEGRAY);
	      else
                graphColor(WHITE);
	      graphFillRectangle(offset, control->topMargin,
                                 offset+200, control->graphHeight);

	      graphColor(BLACK);
	    }
	  
	  (*instance->draw)(instance, &offset);
	  if (!firstFromBox && control->fromBox)
	    firstFromBox = control->fromBox;
	  
	  if (!control->hideHeader)
	    { char *cp;
	      char tmp[2];
	      int bottomBox = graphBoxStart();
	      int oh = graphTextHeight(0.5);
	      assInsert(control->bottomBoxes,
			assVoid(bottomBox),
			instance);
	      graphLine(oldOffset, control->graphHeight, 
			oldOffset, control->graphHeight+0.6);
	      graphLine(offset, control->graphHeight,
			offset, control->graphHeight+0.6);
	      for (cp = instance->name ; 
		   *cp && oldOffset<offset-1; 
		   cp++, oldOffset += 0.8)
		{ tmp[0] = *cp;
		  tmp[1] = 0;
		  graphText(tmp, oldOffset+0.5, control->graphHeight);
		}
	      graphBoxEnd();
	      graphTextHeight(oh);
	      graphBoxDraw(bottomBox, BLACK, instance->map->colour);
	    }
	  
	}
      array(control->mapTransitions, i, float) = offset;
    }  

  graphColor (WHITE) ;
  graphFillRectangle(offset, control->topMargin,
		     offset+200, control->graphHeight);
  graphColor (BLACK) ;

  if (control->currentMap->menu) 
    graphMenu(control->currentMap->menu);

  controlOverlay() ;
  graphRedraw () ;


  /* ACEDB-GRAPH INTERFACE:                                                  */
  if (getGraphAcedbResetKey() != NULL) (getGraphAcedbResetKey())(&(control->from)) ;

  if (firstFromBox)
    controlSelectBox(firstFromBox);

  return;
} /* controlDrawControl */


void controlDraw(void)
{
  COLCONTROL control = currentColControl("controlDraw");

  controlDrawControl((void *) control);
} /* controlDraw */


void controlConvertRegister(COLCONTROL control)
{ COLINSTANCE instance;
  MAPCONTROL map;
  int i;
  
  for (i=0; i<arrayMax(control->maps); i++)
    { map = array(control->maps, i, MAPCONTROL);
      arrayMax(map->transConversionRoutines) = 0;
    }
  
  for(i=0; i<arrayMax(control->instances); i++)
    { instance = array(control->instances, i, COLINSTANCE);
      if (instance->conversionRegister)
	(*instance->conversionRegister)(control, instance);
    }
}

void **controlConvert(Array conversionRoutines,
		      MAPCONTROL map, 
		      void *(*conRoutine)(), 
		      void *params)
{ struct ConversionRecord *c;
  int i;

  for(i=0; i<arrayMax(map->transConversionRoutines); i++)
    { c = arrp(map->transConversionRoutines, i, struct ConversionRecord );
      if (c->convertRoutine == conRoutine && c->params == params) 
	return &c->results ;
    }
  
  c = arrayp(map->transConversionRoutines,
	     arrayMax(map->transConversionRoutines),
	     struct ConversionRecord);
  c->convertRoutine = conRoutine;
  c->params = params;
  return &c->results;
}

void controlConvertTrans(MAPCONTROL map)
{ int i;
  struct ConversionRecord *c;
  
  for (i=0; i<arrayMax(map->transConversionRoutines); i++)
    { c = arrp(map->transConversionRoutines, i, struct ConversionRecord);
      c->results = (*c->convertRoutine)(map, c->params);
    }
}

void controlConvertPerm(MAPCONTROL map)
{ int i;
  struct ConversionRecord *c;
  
  for (i=0; i<arrayMax(map->permConversionRoutines); i++)
    { c = arrp(map->permConversionRoutines, i, struct ConversionRecord);
      c->results = (*c->convertRoutine)(map, c->params);
    }
}
      

void controlAddButton(MAPCONTROL map, COLOUROPT *button, MENUOPT *menu)
{
  COLOUROPT **p = map->buttonsAddHere;
 
  *p = (COLOUROPT *)halloc(sizeof(COLOUROPT), map->handle);
  **p = *button; /* copies whole struct */
  (*p)->next = *p;
  map->buttonsAddHere = &((*p)->next);
  
  array(map->menusOnButtons, arrayMax(map->menusOnButtons), MENUOPT *) = menu;

  return;
} /* controlAddButton */


float controlDrawButtons(MAPCONTROL map, float x, float y, float max)
{ 
  int i;
  int box = graphColouredButtons(map->buttons, x, &y, max);
  MENUOPT *menu;
  float y2 ;
  
  for (i=0; i<arrayMax(map->menusOnButtons); i++, box++)
    { menu = arr(map->menusOnButtons, i, MENUOPT *);
      if (menu) 
	graphBoxMenu(box, menu);
      
      
      /* ACEDB-GRAPH INTERFACE: acedb may need to free a popped up menu.     */
      if (getGraphAcedbFreemenu() != NULL)
	{
	  if (i == 0 && map->viewMenu) (getGraphAcedbFreemenu())(box, map->viewMenu) ;
	}
    }
  
  graphBoxDim (--box, 0, 0, 0, &y2) ;
  return y2 ;
}


static double dragOldy;

static void controlLeftDrag(double x, double y)
{
  COLCONTROL control = currentColControl("controlLeftDrag");

  graphXorLine(0, dragOldy, control->graphWidth, dragOldy);
  dragOldy = y;
  graphXorLine(0, y, control->graphWidth, y);

  return;
} /* controlLeftDrag */

static void controlLeftUp(double x, double y)
{
  COLCONTROL control = currentColControl("controlLeftUp");

  graphXorLine(0, dragOldy, control->graphWidth, dragOldy);
  graphRegister(LEFT_DRAG, 0);
  graphRegister(LEFT_UP, 0);
  
  return;
} /* controlLeftUp */

static void controlPick(int box, double x, double y)
{
  COLCONTROL control = currentColControl("controlPick");
  COLINSTANCE instance;
  float x1,x2,y1,y2;

  if (box && box <= control->lastHeaderBox) /* picked a header box, despatch to map */
    { if (control->currentMap->headerPick)
	(*control->currentMap->headerPick)(control->currentMap, box, x, y);
      return;
    }

  if (box && assFind(control->bottomBoxes, assVoid(box), &instance))
    /* a bottom coloured box, click here to configure */
    {
    if (instance->configure)
      {
      if ((*instance->configure)(instance))
	{ 
	/* ACEDB-GRAPH INTERFACE: set acedb specific map field.              */
	if (getGraphAcedbResetKey() != NULL) (getGraphAcedbResetKey())(&(instance->map->viewKey)) ;

	controlDrawControl(control);
	}
      }
    return;
    }

  /* Second click: follow it */
  if (box && box == control->activeBox)
    { instance = array(control->boxIndex, box, COLINSTANCE);
      if (instance && instance->followBox)
	(*instance->followBox)(instance, box, x, y);
      return;
    }
  
/*  if (!box)
    controlMiddleUp(x,y);  il rm 12/6/97 srk added for some reason 4/3/97 ??? */

  instance = array(control->boxIndex, box, COLINSTANCE);
  if (instance && instance->pick)
    { (*instance->pick)(instance, box, x, y);
      return;
    }

  controlSelect(box, x, y);

  if (box == 0 && control->needRuler)
    { graphBoxDim(0, &x1, &y1, &x2, &y2);
      y += y1;
      dragOldy = y;
      graphXorLine(0, y, control->graphWidth, y);
      graphRegister(LEFT_DRAG, (GraphFunc)controlLeftDrag);
      graphRegister(LEFT_UP,   (GraphFunc)controlLeftUp);
    }
  
  return;
} /* controlPick */

static void controlSelect(int box, double x, double y)
{
  COLCONTROL control = currentColControl("controlSelect");
  COLINSTANCE instance;
  MAPCONTROL newMap = control->currentMap;
  BOOL redrawOld = FALSE, redrawNew = FALSE, redrawMap = FALSE;
  int i;

  /* first deselect any previous box */

  if (control->activeBox)
    { instance = array(control->boxIndex, control->activeBox, COLINSTANCE);
      if (instance && instance->unSelectBox) 
	redrawOld = (*instance->unSelectBox)(instance, control->activeBox);
    }

  /* Now set the new one as selected */
  control->activeInstance = 0;

  /* ACEDB-GRAPH INTERFACE: set acedb specific field.                        */
  if (getGraphAcedbResetKey() != NULL) (getGraphAcedbResetKey())(&(control->activeKey)) ;


  if (box == 0 && x > 0) /* Pick in no box, can still change map*/
    { /* switch current map if necessary */
      for(i=0; i<arrayMax(control->instances); i++)
	{ if (array(control->mapTransitions, i, float) > x)
	    { newMap = array(control->instances, i, COLINSTANCE)->map;
	      break;
	    }
	}
    }
  else
     { 
       instance = array(control->boxIndex, box, COLINSTANCE);
       if (instance && instance->setSelectBox)
	 { control->activeInstance = instance;
	   redrawNew = (*instance->setSelectBox)(instance, box, x, y);
	   newMap = instance->map;
	 }
     }

  redrawMap = newMap != control->currentMap;
  control->currentMap = newMap;

  if (redrawOld || redrawNew || redrawMap) 
    controlDraw();
  else
    controlOverlay();

  return;
} /* controlSelect */


void controlOverlay(void)
{
  COLCONTROL control = currentColControl("controlOverlay");
  COLINSTANCE instance;
  int i;

  control->activeBox = 0;
  for (i=0; i<arrayMax(control->boxIndex); i++)
    { instance = array(control->boxIndex, i, COLINSTANCE);
      if (instance && instance->doColour)
	(*instance->doColour)(instance, i);
      /* this call sets control->activeBox */
    }
  
  return;
} /* controlOverlay */


void controlUnselectAll(void)
/* Public function to unhighlight the currently highlit thing */
{
  controlSelect(0, -1, 0);
}

void controlSelectBox (int box)
/* Public function to highlight the thing in box */
{ 
  controlSelect (box, -1, 0);
}

static void controlKybd(int k)
{
  MAPCONTROL map = currentMapControl();

  if (map->keyboard)
    (*(map->keyboard))(k);
}

void controlPrint(void)
{
  COLCONTROL control = currentColControl("controlPrint");
  MAPCONTROL map = control->currentMap;
  float oldMin = map->min;
  float oldMax = map->max;
  float oldCentre = map->centre;
  BOOL oldHH = control->hideHeader;
  float min, max, mag;
  int nx, ny;

  graphFitBounds (&nx, &ny) ;
  
  if (!messPrompt("Please state the zone you wish to print",
	     messprintf("%g   %g", map->min, map->max), "ffz"))
    return ;
  freefloat(&min) ; freefloat(&max) ;
  if (min >= max)
    return ;
  
  map->min = min;
  map->max = max;
  map->centre = (max+min)/2.0;
  control->hideHeader = TRUE;

  mag = map->mag > 0 ? map->mag : -map->mag;
  
  graphBoundsPrint (nx + 0.2, 1.05 * (max - min) * mag + 5, controlDraw) ;
  
  map->min = oldMin;
  map->max = oldMax;
  map->centre = oldCentre;
  control->hideHeader = oldHH;

  controlDrawControl(control);
}

/************************************/
/* a few utility functions (public) */
/************************************/

/* set/read box to private data mapping */
void controlRegBox(COLINSTANCE instance, int box, void *private)
{ COLCONTROL control = instance->map->control;

  array(control->boxIndex, box, COLINSTANCE) = instance;
  array(control->boxIndex2, box, void *) = private;
}

void *controlBoxRegd(COLINSTANCE instance, int box)
{ return arr(instance->map->control->boxIndex2, box, void *);
}

/* utility func which helps truncate column drawing:
   coord is y value in map coords, slop is how far above/below that we will
   draw, ret is given translated graph co-ord, truncated if needed.
   returns FALSE is truncation done, TRUE otherwise . */
BOOL controlCheckY(MAPCONTROL map, float coord, float slop, float *ret)
{ COLCONTROL control = map->control;
  float gco = MAP2GRAPH(map, coord);
  
  if (gco < (control->topMargin + slop + 0.1))
    { if (ret)
	*ret = control->topMargin + 0.1;
      return FALSE;
    }

  if (gco > (control->graphHeight - (0.1 + slop)))
    { if (ret)
	*ret = control->graphHeight - 0.1; 
      return FALSE;
    }

  if (ret)
    *ret = gco;
  return TRUE;
}

void controlTruncatedLine(COLCONTROL control, float x1, float y1, 
			  float x2, float y2)
     /* x2, y2 is assumed to be between topMargin and graphHeight */
     /* x1, y1 can be outside, draw part if line bewteen tm and gh */
{ float t = control->topMargin;
  float b = control->graphHeight;

  if (y1 > t && y1 < b)
    graphLine(x1, y1, x2, y2);
  else if (y1 <= t)
    { float crossx = x1 + (x2-x1)*(control->topMargin-y1)/(y2-y1);
      graphLine(crossx, control->topMargin, x2, y2);
    }
  else
    { float crossx = x1 + (x2-x1)*(y1-control->graphHeight)/(y1-y2);
      graphLine(crossx, control->graphHeight, x2, y2);
    }
}



/*************************************************************************/
/* A few columns for use anywhere, which depend only on map information. */
/*************************************************************************/

static void mapLocatorDrawInit(COLINSTANCE instance)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control ;
  float y0, yn;

  y0 = 1 + control->topMargin ;         
  yn = control->graphHeight - 1;        

  if (map->mag < 0)
    { map->thumb.offset = map->max ;
      map->thumb.fac = (y0 - yn) / (map->max - map->min) ;
    }
  else
    { map->thumb.offset = map->min ;
      map->thumb.fac = (yn - y0) / (map->max - map->min) ;
    }

  map->thumb.halfwidth = 0.5 * map->thumb.fac *
		      (control->graphHeight-control->topMargin-2.0) / map->mag;

  map->thumb.x = 0 ;
  map->thumb.drawn = FALSE;
  map->cursor.drawn = FALSE; /* In case scale instance removed */

}


static void mapLocatorDraw (COLINSTANCE instance, float *offset)
{ float y;
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control ;
  float t, b, box1, box2;

  *offset += 1;

  map->thumb.x = *offset;
  map->thumb.drawn = TRUE;

  box1 = graphBoxStart () ;
  graphFillRectangle (*offset - 0.25, MAP2WHOLE(map, map->min),
		      *offset + 0.25, MAP2WHOLE(map, map->max));
  y = MAP2WHOLE(map, map->centre) ;
  graphLine (*offset - 0.5, y, *offset + 0.5, y) ; /* to give correct width */

  y = MAP2WHOLE(map, map->centre);
  box2 = graphBoxStart();
  
  t = y - map->thumb.halfwidth;
  b = y + map->thumb.halfwidth;
  if (t < control->topMargin)
    t = control->topMargin;
  if (b > control->graphHeight)
    b = control->graphHeight;

  if (b > control->topMargin && 
      t < control->graphHeight)
    graphRectangle (*offset - 0.5, t, *offset + 0.5, b) ;

  graphBoxEnd();  /* box2 */
  graphBoxDraw (box2, DARKGREEN, GREEN);  
  
  graphBoxEnd();  /* box1 */

  if (isGifDisplay) /* fall down on big bar */
    {
      map->thumb.box = box1 ;
      graphBoxSetPick (box2, FALSE) ; 
    }
  else  /* do not react on left button, use mid button to recenter */
    {
       map->thumb.box = box2 ;
       graphBoxSetPick (box1, FALSE) ; 
    }
  array(control->boxIndex, map->thumb.box, COLINSTANCE) = instance;

  if (map->cursor.drawn && map->hasProjectionLines)
    { controlTruncatedLine(control, *offset,
			   y - map->thumb.halfwidth,
			   map->cursor.x,
			   control->topMargin+1);
      controlTruncatedLine(control, *offset,
			   y + map->thumb.halfwidth,
			   map->cursor.x,
			   control->graphHeight-1);
    }

  *offset += 1;

}

static void mapThumbDrag (float *x, float *y, BOOL isDone)
{
  COLCONTROL control = currentColControl("mapThumbDrag");
  MAPCONTROL map = control->thumbMap;
  float top, bottom;

  /* fix x */  
  *x = map->thumb.x - 0.5;

  /* stop the user pulling the thumb bar out of the extent */
  top = MAP2WHOLE(map, map->min);
  bottom = MAP2WHOLE(map, map->max);
  if (map->thumb.fac < 0 )
    { float tmp = bottom; 
      bottom = top;
      top = tmp;
    }
  if (*y < top)
    *y = top;
  if ((*y + 2 * map->thumb.halfwidth) > bottom)
    *y = bottom - 2 * map->thumb.halfwidth;
  
  if (isDone)
    { map->centre = WHOLE2MAP(map, *y + map->thumb.halfwidth);
      controlDraw();
    }
  
}


static void mapLocatorPick(COLINSTANCE instance,
			     int box, 
			     double x, 
			     double y)
{
  instance->map->control->thumbMap = instance->map;
  if (!isGifDisplay)
    graphBoxDrag(box, mapThumbDrag);
  else   
    {  
      double top = MAP2WHOLE(instance->map, instance->map->min);
      instance->map->centre =
	WHOLE2MAP (instance->map, y +  top) ;
      controlDraw() ;  
    }
  
} /* mapLocatorPick */
 
struct configMapLocal{
  float magConf;
  BOOL lines;
};


static BOOL mapLocatorConfigure (COLINSTANCE instance)
{
  MAPCONTROL map = instance->map;
  struct configMapLocal *cf = (struct configMapLocal *) messalloc(sizeof(struct configMapLocal));
  float line = 2.0;

  if(controlCreateConfig(instance,cf,"Configure Locator column",0.5,0.15)){
    /* initialise data */
    cf->magConf = map->magConf;
    cf->lines =  map->hasProjectionLines;
 
    graphFloatEditor("Magnification :",&cf->magConf,4.0,line++,0);
    graphToggleEditor("Show projection Lines",&cf->lines,4.0,line++);
    graphRedraw();
  }

  return FALSE;
} /* mapLocatorConfigure */


static void mapLocatorConfigFinal (COLINSTANCE instance, void *locals, BOOL ok)
{
  MAPCONTROL map = instance->map;
  struct configMapLocal *cf = locals;
  
  if(ok)
    {
      map->magConf = cf->magConf;
      if (map->magConf != 0.0)
	map->mag = map->magConf;
      map->hasProjectionLines = cf->lines;
    }
  else
    messfree(cf);

  return;
} /* mapLocatorConfigFinal */



static BOOL mapLocatorCreate (COLINSTANCE instance, OBJ init)
{ 
  MAPCONTROL map = instance->map;

  instance->draw = mapLocatorDraw;
  instance->drawInit = mapLocatorDrawInit;
  instance->pick = mapLocatorPick;
  instance->configure = mapLocatorConfigure;
  instance->configFinal = mapLocatorConfigFinal;

  map->hasProjectionLines = TRUE; /* default */
  map->magConf = map->mag ;	 /*  map->magConf = 0;  mhmp 24.10.97*/


  /* ACEDB-GRAPH INTERFACE: Call acedb specific code if registered to get    */
  /* magnification and projection.                                           */
  if (getGraphAcedbMapLocate() != NULL)
    (getGraphAcedbMapLocate())(init, &(instance->save), map) ;


  return TRUE ;
} /* mapLocatorCreate */




static void mapWhole (void *m)
{ 
  BOOL isNeg ;
  MAPCONTROL map = (MAPCONTROL) m; 


  map->centre = (map->max + map->min) / 2.0 ;  
  isNeg = (map->mag < 0) ;
  map->mag = (map->control->graphHeight - map->control->topMargin - 2.0) /
    			(map->max - map->min) ;
  
  if (isNeg)
    map->mag = -map->mag ;

  controlDraw() ;

  return;
} /* mapWhole */


static void mapZoomIn (void *m)
{ 
  MAPCONTROL map = (MAPCONTROL) m;

  map->mag *= 2 ;
  controlDraw() ;

  return;
} /* mapZoomIn */

static void mapZoomOut (void *m)
{
  MAPCONTROL map = (MAPCONTROL) m; 

  map->mag /= 2 ; 
  controlDraw() ;

  return;
} /* mapZoomOut */


static void mapLocatorInit(COLPROTO proto, MAPCONTROL map)
{
  COLOUROPT *button;
  
  button = (COLOUROPT *) messalloc(sizeof(COLOUROPT));

  button->text = "Whole";
  button->f = mapWhole;
  button->arg = (void *) map;
  button->fg = BLACK;
  button->bg = WHITE;
  button->next = 0;
  controlAddButton(map, button, 0);

  button->text = "Zoom in";
  button->f = mapZoomIn;
  button->arg = (void *) map;
  button->fg = BLACK;
  button->bg = WHITE;
  button->next = 0;
  controlAddButton(map, button, 0);

  button->text = "Zoom out";
  button->f = mapZoomOut;
  button->arg = (void *) map;
  button->fg = BLACK;
  button->bg = WHITE;
  button->next = button; /* Terminator */
  controlAddButton(map, button, 0);
  
  messfree(button);

  return;
} /* mapLocatorInit */

  

struct ProtoStruct mapLocatorColumn = {
  mapLocatorInit,
  mapLocatorCreate,
  0,
  "Locator",
  0,
  TRUE,
  0,
  0,
  "The locator column displays a green slider bar which can be used to "
    "control the position of the map display.\n\n"
      "The magnification configuration option can be used to set a specific "
	"factor between map co-ordinates and screen distances in character "
	  "size units. If a configuration is saved with a magnification set, "
	    "that will be used for initial display in preference to defaults "
	      "or information derived from the displayed object.\n\n"
		"The \"Show projection lines\" toggle controls the lines "
		  "drawn between the green bar on the locator and the ends "
		    "of the scale, if it exists.\n"
};

/************************************************************/
/*   scroller, to be used from giface/netscape  */


static void mapScrollerDrawInit(COLINSTANCE instance)
{
}


static void mapScrollerDraw(COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control ;
  int i, box ;
  float t, b, y1, y2 ;
  Array aa = 0 ;
  *offset += 1.5;
  
  box = graphBoxStart () ;
  

  t = GRAPH2MAP(map, control->topMargin + 1) ;
  t = MAP2WHOLE (map, t) ;
   if (t < control->topMargin + .2)
     t =  control->topMargin + .2 ;
  b = GRAPH2MAP(map, control->graphHeight-1) ;
  b = MAP2WHOLE (map, b) ;
  if (b > control->graphHeight)
    b = control->graphHeight ;

  graphFillRectangle (*offset - 0.5, t, *offset + 0.5, b) ;
  
  aa = arrayCreate (8, float) ; i = 0 ;
  y1 = control->topMargin + 1. ;
  y2 = control->topMargin + 2. ;
  
  if (y1 < t )
    {
      array (aa, i++, float) = *offset ;
      array (aa, i++, float) = y1 ;
      array (aa, i++, float) = *offset + .5 ;
      array (aa, i++, float) = y2 ;
      array (aa, i++, float) = *offset - .5 ;
      array (aa, i++, float) = y2 ;
      array (aa, i++, float) = *offset ;
      array (aa, i++, float) = y1 ;
      graphPolygon (aa) ;
      graphLine (*offset, t, *offset, y2) ;
    }

  arrayDestroy (aa) ;
  aa = arrayCreate (8, float) ;
  y1 = control->graphHeight - 1. ;
  y2 = control->graphHeight - 2. ;
  
  i = 0 ;
  if (y1 > b )
    {
      array (aa, i++, float) = *offset ;
      array (aa, i++, float) = y1 ;
      array (aa, i++, float) = *offset + .5 ;
      array (aa, i++, float) = y2 ;
      array (aa, i++, float) = *offset - .5 ;
      array (aa, i++, float) = y2 ;
      array (aa, i++, float) = *offset ;
      array (aa, i++, float) = y1 ;
      graphPolygon (aa) ;
      graphLine (*offset, b, *offset, y2) ;
    }
  arrayDestroy (aa) ;
  graphBoxEnd();

  graphBoxDraw (box, DARKGREEN, WHITE);
  /* graphBoxInfo (box, 0, "BUTTON:\"Click here to recentre\"") ; */
  array(control->boxIndex, box, COLINSTANCE) = instance;

  *offset += 1;
}

static void mapScrollerInit(COLPROTO proto, MAPCONTROL map)
{ 
}

  
static void mapScrollerPick(COLINSTANCE instance,
			     int box, 
			     double x, 
			     double y)
{
  instance->map->centre = WHOLE2MAP (instance->map, 
				     y + instance->map->control->topMargin + 1) ;
  controlDraw() ;
}
 
static BOOL mapScrollerCreate(COLINSTANCE instance, OBJ init)
{ 
  instance->draw = mapScrollerDraw;
  instance->drawInit = mapScrollerDrawInit;
  instance->pick = mapScrollerPick;

  return TRUE ;
}


struct ProtoStruct mapScrollerColumn = {
  mapScrollerInit,
  mapScrollerCreate,
  0,
  "Scroller",
  0,
  0,
  0,
  0,
  "The scroller column displays a scroll bar which is "
   " sensitive to the left mouse button and hance may by "
    " used from Netscpape. \n\n"
};

/************************************************************/


void mapControlCursorSet (MAPCONTROL map, float x)
{
  if (x > 0)
    map->cursor.val = 0.5 + x / map->cursor.unit ;
  else
    map->cursor.val = -0.5 + x / map->cursor.unit ;

  if (map->hasCursor && !map->cursor.box) /* off screen before, draw it now. */
    controlDrawControl(map->control);
  else
    mapControlCursorShift (map) ;
}


void mapControlCursorShift (MAPCONTROL map)
{
  float x1, x2, y1, y2, z = mapControlCursorPos(map) ;

  if (!map->cursor.box)
    return ;
  graphBoxDim (map->cursor.box, &x1, &y1, &x2, &y2) ;
  strcpy (map->cursor.text, messprintf ("%.*f", map->cursor.resolution, z)) ;
  y1 = MAP2GRAPH(map, mapControlCursorPos(map)) ;
  graphBoxShift (map->cursor.box, x1, y1-0.5) ;
}

static void mapControlFindScaleUnit (MAPCONTROL map, float *u, float *sub)
{
  float cutoff = 5 / map->mag ;
  float unit = *u ;
  float subunit = *u ;

  if (cutoff < 0)
    cutoff = -cutoff ;

  while (unit < cutoff)
    { unit *= 2 ;
      subunit *= 5 ;
      if (unit >= cutoff)
	break ;
      unit *= 2.5000001 ;	/* safe rounding */
      if (unit >= cutoff)
	break ;
      unit *= 2 ;
      subunit *= 2 ;
    }
  subunit /= 10 ;
  if (subunit > *sub)
    *sub = subunit ;
  *u = unit ;
}



static void mapScaleDrawInit(COLINSTANCE instance)
{
  MAPCONTROL map = instance->map;

  map->cursor.x = 0;
  map->cursor.drawn = FALSE;
  map->thumb.drawn = FALSE; /* In case locator instance removed */
}


static void mapScaleDraw (COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  float unit = map->scaleUnit;
  float subunit = unit/10.0;
  float x, y ;
  int resolution, max = 0 ;
  char *cp = 0 ;
  COLCONTROL control = map->control ;

  map->cursor.x = *offset+1.5;
  map->cursor.drawn = TRUE;
  map->cursor.box = 0;

  graphColor (BLACK) ;
  mapControlFindScaleUnit (map, &unit, &subunit) ;
  if (unit >= 1)
    resolution = 0 ;
  else if (unit >= .1)
    resolution = 1 ;
  else if (unit >= .01)
    resolution = 2 ;
  else if (unit >= .001)
    resolution = 3 ;
  else 
    resolution = 0 ;

  x = GRAPH2MAP(map, control->topMargin+1) ;
  x = unit * ((((x>=0)^(map->mag<0))?1:0) + (int)(x/unit)) ;
  while ((y = MAP2GRAPH(map, x)) < control->graphHeight - 0.95)
    { graphLine (*offset+0.5,y,*offset+1.5,y) ;
      cp = messprintf ("%-4.*f", resolution, x) ;
      graphText (cp, *offset+2, y-0.5) ;
      if (strlen(cp)+3 > max)
	max = 3+strlen(cp) ;
      if (map->mag >0)
	x += unit ;
      else
	x -= unit;
    }

  x = GRAPH2MAP(map, control->topMargin+1) ;
  x = subunit * ((((x>=0)^(map->mag<0))?1:0) + (int)(x/subunit)) ;
  while ((y = MAP2GRAPH(map, x)) < control->graphHeight - 1.5)
    { graphLine (*offset+1.0,y,*offset+1.5,y) ;
      if (map->mag >0)
	x += subunit ;
      else
	x -= subunit;
    }
  
  graphLine (*offset+1.5, control->topMargin+1, 
	     *offset+1.5, control->graphHeight-1 ) ;


  if (map->thumb.drawn && map->hasProjectionLines)
    { controlTruncatedLine(control, map->thumb.x, 
			   MAP2WHOLE(map, map->centre) - map->thumb.halfwidth,
			   *offset+1.5, control->topMargin+1) ;
      controlTruncatedLine(control, map->thumb.x, 
			   MAP2WHOLE(map, map->centre) + map->thumb.halfwidth,
			   *offset+1.5, control->graphHeight-1.0) ;
    }

  if (map->hasCursor && !control->hideHeader)
    {
      float z = mapControlCursorPos(map) ;
      float y = MAP2GRAPH(map, z) ;	/* Jean - your x here was a bug */
      float thumbx = map->thumb.drawn ? map->thumb.x : 0;
      
      if (y > control->topMargin+1 && y < control->graphHeight -1)
	{
	  if (map->cursor.unit >= .99)
	    map->cursor.resolution = 0 ;
	  else if (map->cursor.unit >= .099)
	    map->cursor.resolution = 1 ;
	  else if (map->cursor.unit >= .0099)
	    map->cursor.resolution = 2 ;
	  else if (map->cursor.unit >= .00099)
	    map->cursor.resolution = 3 ;
	  else
	    map->cursor.resolution = 0 ;
	  
	  strcpy (map->cursor.text, 
		  messprintf ("%.*f", map->cursor.resolution, z)) ;
	  map->cursor.box = graphBoxStart() ;
	  graphLine (thumbx, y, control->graphWidth+1, y);
	  map->cursor.pickBox = graphBoxStart() ;
	  
	  if (max < 3.0+strlen(map->cursor.text))
	    max = 3.0+strlen(map->cursor.text);
	  graphColor (LIGHTGREEN) ;
	  graphFillRectangle (*offset+0.5, y-0.5, 
			      *offset+max-0.5, y+0.5) ;
	  graphColor (BLACK) ;
	  graphTextPtr (map->cursor.text, *offset+0.5, 
			y-0.5, strlen(map->cursor.text)) ;
	  graphBoxEnd () ;
	  graphBoxEnd () ;
	  graphBoxSetPick (map->cursor.box, FALSE) ; /* only pick on .pickBox */
	  array(control->boxIndex, map->cursor.pickBox, COLINSTANCE) = instance;
	  graphBoxDraw (map->cursor.box, BLACK, TRANSPARENT) ;
	}
    }

  *offset += max ;
}

static void mapControlCursorDrag (float *x, float *y, BOOL isDone)
{
  COLCONTROL control = currentColControl("mapControlCursorDrag");
  MAPCONTROL map = control->thumbMap;

  *x = map->cursor.x - 1.0;
  if (*y < control->topMargin + 0.5)
    *y = control->topMargin + 0.5;
  if (*y > control->graphHeight - 1.5)
    *y = control->graphHeight - 1.5;
  
  if (isDone)
    mapControlCursorSet (map, GRAPH2MAP(map, *y+0.5)) ;
}


static void mapScalePick(COLINSTANCE instance,
			   int box, 
			   double x, 
			   double y)
{ 
  instance->map->control->thumbMap = instance->map;
  graphBoxDrag(box, mapControlCursorDrag);
}

struct configLocalsName{
  float scale,cursor;
  BOOL showCursor;
};

static BOOL mapScaleConfigure(COLINSTANCE instance)
{ MAPCONTROL map = instance->map;
  struct configLocalsName *cf = (struct configLocalsName *) messalloc(sizeof(struct configLocalsName));
  float line = 2.0;
  /* internal representation depends on cursor unit */
  
  if(controlCreateConfig(instance,cf,"Configure Scale column",0.5,0.15)){
    /* initialise data */
    cf->scale = map->scaleUnit;
    cf->cursor = map->cursor.unit;
    cf->showCursor = map->hasCursor;

    graphFloatEditor("Scale unit:",&cf->scale,4.0,line++,0);
    graphFloatEditor("Cursor unit:",&cf->cursor,4.0,line++,0);
    graphToggleEditor("Show cursor",&cf->showCursor,4.0,line++);
    graphRedraw();
  }
  return FALSE;
}
static void mapScaleConfigFinal(COLINSTANCE instance, void *locals, BOOL ok){
  struct configLocalsName *cf = locals;
  MAPCONTROL map = instance->map;
  float old = mapControlCursorPos(map);

  if(ok){
    map->hasCursor = cf->showCursor;
    map->scaleUnit = cf->scale;
    if (map->scaleUnit <0.0001)
      map->scaleUnit = 0.0001;
    map->cursor.unit = cf->cursor;
    if (map->cursor.unit < 0.001)
      map->cursor.unit = 0.001;
    
    if (old > 0)
      map->cursor.val = 0.5 + old / map->cursor.unit ;
    else
      map->magConf = 0.0;
  }
  else
    messfree(cf);
}


static BOOL mapScaleCreate(COLINSTANCE instance, OBJ init)
{
  MAPCONTROL map = instance->map;

  instance->draw = mapScaleDraw;
  instance->drawInit = mapScaleDrawInit;
  instance->pick = mapScalePick;
  instance->configure = mapScaleConfigure;
  instance->configFinal = mapScaleConfigFinal;

  map->cursor.val = 0; 
  map->cursor.unit = 1.0;
  map->scaleUnit = 0.01;
  map->hasCursor = FALSE;


  /* ACEDB-GRAPH INTERFACE: Call acedb specific code if registered to get    */
  /* scale information.                                                      */
  if (getGraphAcedbScale() != NULL)
    (getGraphAcedbScale())(init, &(instance->save), map) ;


  return TRUE ;
}

struct ProtoStruct mapScaleColumn = {
  0,
  mapScaleCreate,
  0,
  "Scale",
  0,
  TRUE, 
  0,
  0,
  "The scale column shows a ruler scale with map co-ordinates, and "
    "optionally a cursor and indicator line.\n\n"
      "The \"cursor unit\" configuration option controls the distance by "
	"which the up and down arrows move the cursor.\n\n"
	  "The \"scale unit\" configuration option controls the smallest "
	    "difference between two scale ticks which can be drawn."
};
  
/************************************************************************/


static void spacerDraw(COLINSTANCE instance, float *offset)
{ SPACERPRIV private = instance->private;
  COLCONTROL control = instance->map->control;
  
  if (private->colour != WHITE) /* white == transparent */
    { graphColor(private->colour);
      graphFillRectangle(*offset, control->topMargin,
			 *offset + private->width, control->graphHeight);
      graphColor(BLACK);
    }
  *offset += private->width;
}

struct spacerLocal{
  float width;
  int colour;
};

static BOOL spacerConfigure(COLINSTANCE instance)
{ SPACERPRIV private = instance->private;
  struct spacerLocal *cf = (struct spacerLocal *) messalloc(sizeof(struct spacerLocal));
  float line = 2.0;

  if(controlCreateConfig(instance,cf,"Configure spacer column",0.5,0.15)){
    /* initialise data */
    cf->width = private->width;
    cf->colour = private->colour;

    graphFloatEditor("Column width:",&cf->width,4.0,line++,0);
    graphColourEditor(" "," ",&cf->colour,4.0,line);
    graphRedraw();
  }
  return FALSE;
}
static void spacerConfigFinal(COLINSTANCE instance, void *locals, BOOL ok){
  SPACERPRIV private = instance->private;
  struct spacerLocal *cf = locals;

  if(ok){
    private->colour = cf->colour;
    private->width = cf->width;
  }
  else
    messfree(cf);
}
    
static BOOL spacerCreate(COLINSTANCE instance, OBJ init)
  {
  SPACERPRIV private = (SPACERPRIV)handleAlloc(0, instance->handle, sizeof(struct SpacerPriv));

  instance->private = private;
  instance->draw = spacerDraw;
  instance->configure = spacerConfigure;
  instance->configFinal = spacerConfigFinal;

  private->width = 1.0;
  private->colour = WHITE;

  
  /* ACEDB-GRAPH INTERFACE: Call acedb specific code if registered to get    */
  /* space information.                                                      */
  if (getGraphAcedbSpace() != NULL)
    (getGraphAcedbSpace())(init, &(instance->save), private) ;


  return TRUE;
  }



struct ProtoStruct spacerColumn = {
  0,
  spacerCreate,
  0,
  "Spacer",
  0,
  FALSE,
  0,
  0,
  "The spacer column provides space in the map display, of user specified "
    "width and colour. The width is specified in units equal to the width  of "
      "a single character on the display"
}; 



/*************************************************/
/************* middle button for thumb **********/

static double oldx, oldy, oldDy ;
static BOOL dragFast ;

static void controlMiddleDrag (double x, double y) 
{
  MAPCONTROL map = currentMapControl();

  if(dragFast)
    { graphXorLine (0, oldy - oldDy, map->thumb.x, oldy - oldDy) ;
      graphXorLine (0, oldy + oldDy, map->thumb.x, oldy + oldDy) ;
    }
  else
    graphXorLine (map->thumb.x, oldy, map->control->graphWidth, oldy) ;

  oldy = y ;

  if(dragFast)
    { oldDy *= exp ((x - oldx) / 25.) ;
      oldx = x ;
      graphXorLine (0, y - oldDy, map->thumb.x, y - oldDy) ;
      graphXorLine (0, y + oldDy, map->thumb.x, y + oldDy) ;
    }
  else
    graphXorLine (map->thumb.x, y, map->control->graphWidth, y) ;
}

static void controlMiddleUp (double x, double y) 
{ float x1,x2,y1,y2 ;
  MAPCONTROL map = currentMapControl();

  if (dragFast)
    { graphBoxDim (map->thumb.box, &x1, &y1, &x2, &y2) ;
      map->mag *= (y2 - y1) / (2. * oldDy) ;
      map->centre = WHOLE2MAP(map, y) ;
    }
  else
    map->centre = GRAPH2MAP(map,y) ;

  controlDraw();
}

static void controlMiddleDown (double x, double y) 
{ float x1 = 0, x2 = 0, y1 = 0 , y2 = 1 ;
  MAPCONTROL map = currentMapControl();

  if (map->thumb.box)
    graphBoxDim (map->thumb.box, &x1, &y1, &x2, &y2) ;
  oldDy = (y2 - y1) / 2. ;
  
  dragFast = (x < map->thumb.x) ;

  if (dragFast)
    { graphXorLine (0, y - oldDy, map->thumb.x, y - oldDy) ;
      graphXorLine (0, y + oldDy, map->thumb.x, y + oldDy) ;
    }
  else
    graphXorLine (map->thumb.x, y, map->control->graphWidth, y) ;
   
  oldx = x ;
  oldy = y ;
  graphRegister (MIDDLE_DRAG, controlMiddleDrag) ;	/* must redo */
  graphRegister (MIDDLE_UP, controlMiddleUp) ;
}

  
  
  
 
 
 
 
 
 
 
 
