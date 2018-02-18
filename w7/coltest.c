/*  File: coltest.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Test column control and graph.
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 17 16:15 1998 (edgrif)
 * * Sep 17 16:14 1998 (edgrif): Had to change the colControlCreate to
 *              a new routine colControlBasicCreate because of graph
 *              restructuring changes.
 * Created: Thu Sep 17 16:13:56 1998 (edgrif)
 *-------------------------------------------------------------------
 */
/* $Id: coltest.c,v 1.1.1.1 2002/07/19 20:22:59 sienkiew Exp $ */

#include "graph.h"
#include "colcontrol.h"


/* all per-map-private data goes here */
typedef struct ExampleLook {
  int from;
  MAPCONTROL map;
} *EXAMPLELOOK;

extern struct ProtoStruct exampleColumn; /* forward decl */

static COLPROTO columnsTab[] = {
  &mapLocatorColumn,
  &exampleColumn,
  &mapScaleColumn,
  0
};

static COLDEFAULT columnSpec = {
 { &mapLocatorColumn, TRUE, "Locator" },
 { &exampleColumn, TRUE, "Example" },
 { &mapScaleColumn, FALSE, "Scale" },
 { 0, 0, 0}
};
 /* scale is initially hidden */

static MENUOPT exampleMenu[] = {
  graphDestroy, "Quit",
  0,0
};

static void exampleDrawHeader(MAPCONTROL map)
{ controlDrawButtons(map, 5, 0.5, map->control->graphWidth);
}
 
main(int argc, char *argv[])
{ 


  COLCONTROL control;
  MAPCONTROL map;
  EXAMPLELOOK look;


  graphInit(&argc, argv);

  control = colControlBasicCreate(FALSE, TEXT_FIT, "example", 0.5, 0.5, 0, 0);
  
  map = mapControlCreate(control, 0); 
    /* 2nd param is finalise routine, if reqd */
  

/* this block is required */  
  map->prototypes = columnsTab; /* the column types we can display */

/* this block is optional */
  map->menu = exampleMenu; /* window menu whenever our map is current */
  map->name = "Test";      /* defaults to "" */
  map->drawHeader = exampleDrawHeader; /* draw the header */
  map->topMargin = 3; /* how much space the header needs */

/* other optional fields:

   map->beforeDraw   called on each map before drawing columns
   map->headerpick   called on click on header area
   map->convertTrans display list stuff
   map->keyboard     called on current map when key pressed
   
*/

  /* This memory will evaporate automatically when the map does */
  look  = (EXAMPLELOOK)handleAlloc(0, map->handle, sizeof(struct ExampleLook));
  
  map->look = look;
  look->map = map;
  
  map->max = 10;
  map->min = -10;
  map->centre = 0;

  controlAddMap(control, map);
  
  controlConvertPerm(map); /* make display lists */
    
  controlReadConfig(map, columnSpec); /* column spec says what to make */

  look->from = 4; /* highlight box displaying this first */
  controlDraw(); /* finally do it */
  look->from = 0; /* don't highlight on future draws */

  graphLoop (FALSE) ;
}

/*********************************************************/
/** convert function: may be shared between many columns */
/*********************************************************/
typedef struct DrawCommand {
  int number;
  float coord;
} *DRAWCOMMAND;

static void *exampleConvert(MAPCONTROL map, void *params)
{ Array dl = arrayHandleCreate(100, struct DrawCommand, map->handle);
  int i;
  int disp = 7;

  /* noddy */
  for (i=0; i<10; i++)
    { DRAWCOMMAND draw = arrayp(dl, i, struct DrawCommand);
      draw->coord = i;
      draw->number = disp++;
    }

  return dl;
}

/******************************/
/*** definition of a column ***/
/******************************/
typedef struct ExamplePriv {
  Array dl; /* drawing list */
  int selno; /* number of selected thing */
} * EXAMPLEPRIV;

static void exampleDraw(COLINSTANCE instance, float *offset)
{ EXAMPLEPRIV private = instance->private;
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  EXAMPLELOOK look = map->look;
  Array dl  = private->dl;
  int i, ibox;
  DRAWCOMMAND draw;
  float y;
  int max = 0, len;

  for(i=0; i<arrayMax(dl); i++)
    { draw = arrayp(dl, i, struct DrawCommand);
      y = MAP2GRAPH(map, draw->coord);
      if (y > control->topMargin && y < control->graphHeight-1)
	{ ibox = graphBoxStart();
	  array(control->boxIndex, ibox, COLINSTANCE) = instance;
	  /* NB the next line depends on not extending dl after drawing */
	  array(control->boxIndex2, ibox, void *) = draw;
	  graphText(messprintf("%d", draw->number), *offset, y);
	  len = strlen(messprintf("%d", draw->number));
	  if (len > max) 
	    max = len;
	  graphBoxEnd();

/* this box gets selected */
	  if (draw->number == look->from)
	    control->fromBox = ibox;
	}
    }
    
  *offset += max+1;
}
 


static BOOL exampleSetSelect(COLINSTANCE instance, int box, double x, double y)
/* Called when a box is left clicked. Does no drawing, simply sets up the
   columns state so that it knows what is selected. Note that the column
   code keeps track of which instance has the selected box in 
   control->activeInstance */
{ 
  EXAMPLEPRIV private = instance->private;
  COLCONTROL control = instance->map->control;
  DRAWCOMMAND draw =  (DRAWCOMMAND) array(control->boxIndex2, box, void *);

  private->selno = draw->number;
  return FALSE; /* return TRUE if re-draw is needed */
}

static BOOL exampleUnselect(COLINSTANCE instance, int box)
/* Called if a box is unselected, and no re-draw is done, must either
   fix of the window, or return TRUE for a re-draw. Doesn't need to 
   change select state, the column control does this by zeroing 
   control->activeInstance */
{ graphBoxDraw(box, BLACK, WHITE);
  return FALSE;
}

static void exampleDoColour(COLINSTANCE instance, int box)
/* Called for every box */
{ COLCONTROL control = instance->map->control;
  EXAMPLEPRIV private = instance->private;
  DRAWCOMMAND draw =  (DRAWCOMMAND) array(control->boxIndex2, box, void *);

  if (control->activeInstance == instance && /* always check this */
      private->selno == draw->number)
    { graphBoxDraw(box, WHITE, RED);
      control->activeBox = box; /* compulsary */
    }
}

static void exampleFollowBox(COLINSTANCE instance, int box, double x, double y)
/* called when selected box is clicked again, intended to open a new display
   for the double-clicked object */
     
{   COLCONTROL control = instance->map->control;
    DRAWCOMMAND draw =  (DRAWCOMMAND) array(control->boxIndex2, box, void *);
    messout (" double clicked on box %d", draw->number);
}

static BOOL exampleCreate(COLINSTANCE instance, OBJ init)
/* in non-ACEDB case init is always zero */
{ EXAMPLEPRIV private = (EXAMPLEPRIV) handleAlloc(0, instance->handle,
						  sizeof(struct ExamplePriv));

  /* this block needed */
  instance->draw = exampleDraw;
/* this block optional */
  instance->private  = private;
  instance->setSelectBox = exampleSetSelect;
  instance->unSelectBox = exampleUnselect;
  instance->doColour = exampleDoColour;
  instance->followBox = exampleFollowBox;
  
  /* also allowed here: 
     instance->beforeDraw   called before any draw routines
     instance->configure    change config of flexible columns
     */
  
  /* Powerfull magic to pick up results from our registered convert function */
  private->dl = *(instance->proto->convertResults);
		  
		  
  return TRUE; /* can return FALSE if not possible to create */
}

static void exampleButtonPressed(void *arg)
{ messout ("Ding dong");
}


static COLOUROPT exampleButton[] = {
  exampleButtonPressed, 0, /* 0 is arg to function */
  BLACK, WHITE,
  "A Button", 0 };

static MENUOPT menuOnButton[] = { 
  graphDestroy, "Quit-from-button",
  0, 0};

static void exampleInit(MAPCONTROL map, COLCONTROL control)
/* called once per prototype from controlAddMap, 
   chief job is to call control addButton*/
{
  controlAddButton(map, exampleButton,  menuOnButton);
}


struct ProtoStruct exampleColumn = { /* it all links in via here */
  exampleInit, /* called from controlAddMap for each prototype */
  exampleCreate, /* called to make an instance */
  0, /* becomes finalisation routine for all instance structures */
  "example", /*  name */
  0, /* void * specialisation, can have two column prototypes  
	identical, except for name and specialisation behave subtly differently
        access via instance->proto->specialisation in create routine */
  FALSE, /* unique, is set only one instance of column allowed */
  exampleConvert, /* call this routine for our display list */
  0 /* with this parameter */
};

