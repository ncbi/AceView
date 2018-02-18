/*  File: mapcontrol.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: package for generic map drawing containing:
 		columns control
		zoom, middle button scroll, locator
		findScaleUnits
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 10:41 1998 (fw)
 * * (mieg) zoom_set_g at least should be part of the map structure
            anyway i disable 
 * * Jul 23 14:51 1998 (edgrif): Add fmap.h public header for function defs
 *      removed from map.h
 * * Dec 24 1997 (rd): made MapColRec local.  Must use mapInsertCol().
 * * Mar 12 18:29 1995 (rd): added MapColRec2 for dymanic specs
 	based on priorities, and changed to copy cols, not use 
	colSwitch.
 * Created: Wed Nov 18 00:40:35 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: mapcontrol.c,v 1.6 2016/11/23 19:15:12 mieg Exp $ */

#include "acedb.h"
#include "key.h"                /* mapKbdScroom keys */
#include "map.h"		/* declare opaque LookStruct */
				/* to be completed by users 
				   of the mapPackage */

#include "fmap.h"		/* public fmap function headers */

/************************************************************/
                     /* public mapPackage globals */

MAP mapGifMap = 0 ;	       /* bypasses graphrequirement for giface */
int mapGraphWidth, mapGraphHeight; /* globals, only valid while drawing */
float topMargin, halfGraphHeight ;

magic_t MAP2LOOK_ASSOC = "MAP2LOOK" ;	/* find LOOK pointer for the map */

/************************************************************/

typedef struct {
  BOOL   isOn ;
  char   *name ;
  MapColDrawFunc func;
  float  priority ;
} MapColRec ;

static magic_t GRAPH2MAP_ASSOC = "GRAPH2MAP"; /* find the MAP on gActive */
static magic_t MAP_MAGIC = "MAP_MAGIC";       /* verify MAP-pointer */

static BOOL isMagPossible = TRUE ;
static MAP selectedMap = 0 ;

#define NAME_SPACE 20

/************************************************************/

static void mapThumbCalc (MAP map) ;

/************************************************************/

MAP currentMap(char *caller)
{
  MAP map;

  if (!(map = mapGifMap) &&
      !(graphAssFind(&GRAPH2MAP_ASSOC, &map)))
    messcrash("%s() could not find MAP on graph", caller);
  if (!map)
    messcrash("%s() received NULL MAP pointer", caller);
  if (map->magic != &MAP_MAGIC)
    messcrash("%s() received non-magic MAP pointer", caller);

  return map;
} /* currentMap */

/******************************/

MAP mapCreate (VoidRoutine draw)
{
  MAP map ;

  map = (MAP) messalloc (sizeof (struct MapStruct)) ;

  map->handle = handleCreate ();
  map->magic = &MAP_MAGIC ;
  map->draw = draw ;

  map->cols = arrayHandleCreate (32, MapColRec, map->handle) ;

  if (!mapGifMap)
    { graphAssociate (&GRAPH2MAP_ASSOC, map) ;
      graphRegister (MIDDLE_DOWN, mapMiddleDown) ;
    }
  selectedMap = map ;

  return map ;
} /* mapCreate */


void mapAttachToGraph (MAP map)
{
  graphAssociate (&GRAPH2MAP_ASSOC, map) ;
  graphRegister (MIDDLE_DOWN, mapMiddleDown) ;

  return;
} /* mapAttachToGraph */


void uMapDestroy (MAP map)
{
  if (map->magic != &MAP_MAGIC)
    messcrash ("uMapDestroy received mon-magic MAP pointer");

  map->magic = 0 ;
  graphAssRemove (&GRAPH2MAP_ASSOC) ;
  handleDestroy (map->handle) ;
  messfree (map);
} /* uMapDestroy */

/*****************************************/

static int mapColOrder (const void *aa, const void *bb)
{
  const MapColRec *a = (const MapColRec*)aa, *b = (const MapColRec*)bb ;
  if (a->priority > b->priority) return 1 ;
  else if (a->priority < b->priority) return -1 ;
  else return strcmp (a->name, b->name) ;
} /* mapColOrder */


BOOL mapInsertCol (MAP map, float priority, BOOL isOn, char *name, 
		   MapColDrawFunc func)
{
  int i ;
  MapColRec r ;

  if (!map)
    return FALSE ;

  for (i = 0 ; i < arrayMax(map->cols) ; ++i)
    if (strcmp (name, arrp(map->cols, i, MapColRec)->name) == 0)
      /* there already */
      return FALSE ;

  r.priority = priority ;
  r.isOn = isOn ;
  r.name = strnew (name, map->handle) ;
  r.func = func ;

  arrayInsert (map->cols, &r, mapColOrder) ;
  return TRUE ;
} /* mapInsertCol */

/**********************************************************/
/************* code to choose display columns *************/

BOOL mapColSetByName (char *text, int isOn)
{
  int i ;
  BOOL debug = FALSE ;
  MapColRec *col ;
  MAP map = currentMap("mapColSetByName") ;

  col = arrp(map->cols,0,MapColRec) ;
  for (i = 0 ; i < arrayMax(map->cols) ; ++i, ++col)
    {
      if (debug)
	printf("%d:: %s\n", i, col->name) ;
      if (strcasecmp (col->name, text) == 0)
	{ 
	  if (isOn == 2)   /* toggle */
	    col->isOn ^= TRUE ;
	  else if (isOn == -1)   /* query */
	    return col->isOn ;
	  else if (isOn)
	    col->isOn = TRUE ;
	  else
	    col->isOn = FALSE ;
	  return TRUE ;
	}
    }

  return FALSE ;		/* not found */
} /* mapColSetByName */


/************/

void mapColToggle (KEY key, int box)
{
  MAP map = currentMap("mapColToggle") ;

  arr(map->cols,key,MapColRec).isOn ^= TRUE ;
  (map->draw)() ;

  return;
} /* mapColToggle */


void mapColMenu (int box)
/****************************************************************
 * Function to build a menu of map-columns of the
 * map associated with the active graph. This menu
 * will be attached to the box with the given number
 * (Taken that we've just drawn a button with that box-number)
 ****************************************************************/
{
  static Associator cols2options = 0 ;
  int n;
  FREEOPT *options ;
  MAP map = currentMap("mapColMenu") ;

  if (!cols2options)
    cols2options = assCreate () ;
  if (assFind (cols2options, map->cols, &options) &&
      arrayMax(map->cols) == options->key)
    { graphBoxFreeMenu (box, mapColToggle, options) ;
      return ;
    }

  options = (FREEOPT*) messalloc ((arrayMax(map->cols)+1) * sizeof(FREEOPT)) ; 
  options->key = arrayMax(map->cols) ;
  options->text = "Columns" ;

  for (n = 0 ; n < arrayMax(map->cols) ; ++n)
    { options[n+1].text = arrp(map->cols,n,MapColRec)->name ;
      options[n+1].key = n ;
    }

  graphBoxFreeMenu (box, mapColToggle, options) ;
  assInsert (cols2options, map->cols, options) ;

  return;
} /* mapColMenu */

/************************************************************/

static void colPick (int i)
{
  int old ; 
  MAP map = currentMap("colPick") ;

  selectedMap = map ;
  if (i < 2)
    return ;

  old = arrp(map->cols,i-2,MapColRec)->isOn ; /* background and button */

  if (old)
    graphBoxDraw (i, BLACK, WHITE) ;
  else
    graphBoxDraw (i, WHITE, BLACK) ;
  arrp(map->cols,i-2,MapColRec)->isOn = old ? FALSE : TRUE;

  return;
} /* colPick */

/***********************************/

static MENUOPT simpleMenu[] = {
  { graphDestroy,	"Quit" },
  { help,		"Help" },
  { 0, 0 }
} ;

void mapColControl (void)
{
  int i, c, nx, ny, x, y, box ;
  MapColRec *col ;
  char *cp, *cq, old ;
  MAP map = currentMap("mapColControl") ;

  col = arrp(map->cols,0,MapColRec) ;

  graphClear () ;
  graphButton ("Return to display", map->draw, 5, 1) ;
  graphText ("Display options:", 1, 3) ;

  graphFitBounds (&nx, &ny) ;
  nx /= NAME_SPACE ; if (!nx) nx = 1 ;
  nx = 1 ; /* mieg hard coded, to get the form of the display fixed */
  x = 2 ; y = 4 ;
  for (c = 0, i = 0 ; c < arrayMax(map->cols) ; ++c, ++col)
    { box = graphBoxStart() ;
      cp = cq = col->name ;
      while (*cq && *cq != ':') cq++ ;
      old = *cq ; if (old) *cq = 0 ;/* mieg, segmentation fault on solaris */
      graphText (cp, x, y) ;
      graphBoxEnd() ;
      if (old) *cq = old ; /* mieg, segmentation fault on solaris */
      if (old)
	graphText (cq, x+ strlen(cp) + .4 , y) ;
      
      if (col->isOn)
	graphBoxDraw (box, BLACK, PALEBLUE) ;
      else
	graphBoxDraw (box, BLACK, WHITE) ;
      if ((++i) % nx)		/* increment i here */
	x += NAME_SPACE ;
      else
	{ x = 2 ;
	  ++y ;
	}
    }
  graphRegister (PICK, colPick) ;
  graphMenu (simpleMenu) ;
  graphRedraw () ;
  
  return;
} /* mapColControl */

/************************************/

void mapDrawColumns (MAP map)	/* topMargin must be set */
{
  float offset = 0 ;
  float oldOffset = 0 ;
  float maxOffset = 0 ;
  float oldPriority = -1000000 ;
  MapColRec *col;
  int i, frameCol, frame = 0 ;
  LOOK look ;

  graphFitBounds (&mapGraphWidth, &mapGraphHeight) ;

  halfGraphHeight = 0.5 * (mapGraphHeight - topMargin) ;

  if (!(look = fMapGifLook) &&
      !graphAssFind (&MAP2LOOK_ASSOC, &look))
    messcrash ("mapDrawColumns() - Can't find look for this map.\n"
	       "The programmer may have forgotten to associate "
	       "MAP2LOOK_ASSOC with the look on this graph.") ;

  selectedMap = map ;
  mapThumbCalc (map) ;
  map->cursor.box = 0 ;
  map->thumb.box = 0 ;
  map->thumb.x = 0 ;

  arraySort (map->cols, mapColOrder) ;

  col = arrp (map->cols, 0, MapColRec) ;

  frameCol = 0 ;
  map->isFrame = FALSE ;
  for (i = 0 ; i < arrayMax(map->cols) ; ++i) 
    {
      /* check start of frame-sensitive region */
      if (col[i].isOn &&
	  map->fixFrame && !map->isFrame &&
	  (col[i].func == fMapShow3FrameTranslation ||
	   col[i].func == fMapShowORF))
	{ frameCol = i ;
	  map->isFrame = TRUE ;
	  oldPriority = -1000000 ; /* so no overlap */
	}
				/* check end of frame-sensitive region */
      if (map->isFrame && ((col[i].func == fMapShowGeneTranslation) ||
			   (col[i].func == fMapShowDownGeneTranslation) ||
			   (col[i].func == fMapShowUpGeneTranslation))) 
	{ if (frameCol && frame < 2)
	    { ++(*map->fixFrame) ;
	      ++frame ;
	      i = frameCol ;	/* loop */
	    }
	  else
	    { *map->fixFrame -= 2 ; /* restore original value */
	      map->isFrame = FALSE ;
	    }
	  oldPriority = -1000000 ; /* so no overlap */
	}
				/* decide whether or not to overlap */
      if (col[i].priority < oldPriority + 0.01001)
	offset = oldOffset ;
      else
	oldOffset = offset ;
				/* draw the column */
      if (col[i].isOn)
	{ map->activeColName = col[i].name ;
	  if (col[i].func)
	    (*col[i].func)(look, &offset) ;
	  if (offset < maxOffset)
	    offset = maxOffset ;
	  else
	    maxOffset = offset ;
	}

      oldPriority = col[i].priority ;
    }

  return;
} /* mapDrawColumns */

/******************************************************/
/************ zoom and display stuff ******************/ 

void mapWhole (void)
{ 
  BOOL isNeg ;
  MAP map = currentMap("mapWhole") ; 

  selectedMap = map ;
  graphFitBounds (&mapGraphWidth, &mapGraphHeight) ;

  halfGraphHeight = 0.5 * (mapGraphHeight - topMargin) ;

  map->centre = (map->max + map->min) / 2 ;  

  isNeg = (map->mag < 0) ;

  /* Impenetrable fudge factors....(EG), the "- 5" and the "1.05" just seem  */
  /* to be hacks to make things come out right....                           */
#ifdef ED_G_NEVER_INCLUDE_THIS_COD
  map->mag = (mapGraphHeight - topMargin - 5) /  (1.05 * (map->max - map->min)) ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
  map->mag = (mapGraphHeight - 3 - topMargin - 3) /  ((map->max - map->min)) ;

  if (isNeg) map->mag = -map->mag ;

  (map->draw) () ;

  return;
} /* mapWhole */


void mapZoomIn (void)
{ 
  MAP map = currentMap("mapZoomIn") ;

  map->zoom_set_G++ ;

  selectedMap = map ;
  map->mag *= 2 ;
  (map->draw) () ;

  return;
} /* mapZoomIn */


void mapZoomOut (void)
{
  MAP map = currentMap("mapZoomOut") ; 

  if (TRUE) /* map->zoom_set_G > 0, very annoying i cannot move from a gene to a cosmid */
    {
      map->zoom_set_G-- ;
      
      selectedMap = map ;
      map->mag /= 2 ; 
      (map->draw)() ;
    }
  return;
} /* mapZoomOut */

static void mapPageUp (void)
{
  float height  ;
  MAP map = currentMap("mapZoomOut") ; 

  selectedMap = map ; 

  graphFitBounds (&mapGraphWidth, &mapGraphHeight) ;
  height = 
    GRAPH2MAP (map, mapGraphHeight)
    -  GRAPH2MAP (map, topMargin) ;  
  map->centre -= height/2.2 ;

  (map->draw)() ;
  return ;
}

static void mapPageDown (void)
{
  float height  ;
  MAP map = currentMap("mapZoomOut") ; 

  selectedMap = map ; 

  graphFitBounds (&mapGraphWidth, &mapGraphHeight) ;
  height = 
    GRAPH2MAP (map, mapGraphHeight)
    -  GRAPH2MAP (map, topMargin) ;  
  map->centre += height/2.2 ;

  (map->draw)() ;
  return ;
}

void mapKbdScroom (int k)
{
  switch (k)
    {
    case LEFT_KEY :
      mapZoomIn () ;
      break ;
    case RIGHT_KEY :
      mapZoomOut () ;
      break ;
    case PAGE_UP_KEY:
      mapPageUp () ;
      mapPageUp () ;
      break ;
    case UP_KEY : 
      mapPageUp () ;
      break ;
    case PAGE_DOWN_KEY: 
      mapPageDown () ;
      mapPageDown () ;
      break ;
    case DOWN_KEY :
      mapPageDown () ;
      break ;
    }
  return ; 
} /* mapKbdScroom */

/*************************************************/
/************* middle button for thumb **********/

static double oldx, oldy, oldDy ;
static BOOL dragFast ;

static void mapMiddleDrag (double x, double y) 
{
  MAP map = currentMap ("mapMiddleDrag") ;

  if(dragFast)
    { graphXorLine (0, oldy - oldDy, map->thumb.x, oldy - oldDy) ;
      graphXorLine (0, oldy + oldDy, map->thumb.x, oldy + oldDy) ;
    }
  else
    graphXorLine (map->thumb.x, oldy, mapGraphWidth, oldy) ;

  oldy = y ;

  if(dragFast)
    { oldDy *= exp ((x - oldx) / 25.) ;
      oldx = x ;
      graphXorLine (0, y - oldDy, map->thumb.x, y - oldDy) ;
      graphXorLine (0, y + oldDy, map->thumb.x, y + oldDy) ;
    }
  else
    graphXorLine (map->thumb.x, y, mapGraphWidth, y) ;

  return;
} /* mapMiddleDrag */

/**************/

void mapNoMag (void)
{ isMagPossible = FALSE ; }

/**************/
static void mapMiddleUp (double x, double y) 
{
  float x1,x2,y1,y2 ;
  MAP map = currentMap("mapMiddleUp") ;

  if (dragFast)
    { graphBoxDim (map->thumb.box, &x1, &y1, &x2, &y2) ;
      if (isMagPossible)
        map->mag *= (y2 - y1) / (2. * oldDy) ;
      map->centre = WHOLE2MAP(map, y) ;
    }
  else
    map->centre = GRAPH2MAP(map,y) ;

  (map->draw) () ;

  return;
} /* mapMiddleUp */

void mapMiddleDown (double x, double y)
{
  float x1,x2,y1,y2 ;
  MAP map = currentMap("mapMiddleDown") ;

  selectedMap = map ;
  if (map->thumb.box)
    { graphBoxDim (map->thumb.box, &x1, &y1, &x2, &y2) ;
      oldDy = (y2 - y1) / 2. ;
    }
  else
    oldDy = 0;
  
  dragFast = (oldDy && x < map->thumb.x) ; /* thumb->box == 0 if locator not shown */

  if (dragFast)
    { graphXorLine (0, y - oldDy, map->thumb.x, y - oldDy) ;
      graphXorLine (0, y + oldDy, map->thumb.x, y + oldDy) ;
    }
  else
    graphXorLine (map->thumb.x, y, mapGraphWidth, y) ;
   
  oldx = x ;
  oldy = y ;
  graphRegister (MIDDLE_DRAG, mapMiddleDrag) ;	/* must redo */
  graphRegister (MIDDLE_UP, mapMiddleUp) ;

  return;
} /* mapMiddleDown */

/****************************************************/
/********************** scale ***********************/

void mapFindScaleUnit (MAP map, float *u, float *sub)
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

  return;
} /* mapFindScaleUnit */

/***********************************************/
/************** whole map locator **************/

static void mapThumbCalc (MAP map) /* sets up whole transformation */
{
  float y0, yn;
  float  delta = map->max - map->min ;

  y0 = 3 + topMargin ;
  yn = mapGraphHeight - 3 ;

  if (delta == 0)
    delta = 1 ;

  if (map->mag < 0)
    { map->thumb.offset = map->max ;
      map->thumb.fac = (y0 - yn) / delta ;
    }
  else
    { map->thumb.offset = map->min ;
      map->thumb.fac = (yn - y0) / delta ;
    }

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  map->thumb.halfwidth = 0.5 * map->thumb.fac *
				(mapGraphHeight-topMargin) / map->mag ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
  /* This is still not right, the half width ends up being too small by a    */
  /* fraction, I think some of the above can be rationalised.                */
  map->thumb.halfwidth = 0.5 * map->thumb.fac *
				(yn - y0) / map->mag ;

  map->thumb.x = 0 ;

  return;
} /* mapThumbCalc */

/*******************************/

void mapThumbDrag (float *x, float *y, BOOL isDone)
{
  static float oldX ;
  static BOOL isDragging = FALSE ; 
  float top, bottom ;
  MAP map = currentMap("scaleBoxDrag") ;

  if (isDragging)
    *x = oldX ;   /* fix x */  
  else
    { oldX = *x ; 
      isDragging = TRUE ;
    }

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
    { 
      isDragging = FALSE ;
      map->centre = WHOLE2MAP(map, *y + map->thumb.halfwidth);
      (*map->draw)() ;
    }

  return;
} /* mapThumbDrag */


  /* look not used but kept because this is a column display */
void mapShowLocator (LOOK look, float *offset)
{
  MAP map = currentMap("mapShowLocator") ;
  float t, b, y, y1, y2 ;
  Array aa = 0 ;
  int box1, box2, i = 0 ;

  *offset += 1.5 ;

  y = MAP2WHOLE(map, map->centre) ;

  box1 = graphBoxStart () ;  /* global exterior box */

  /* insures correct size of pickable box */

  y1 = topMargin + 1. ;
  y2 = topMargin + 2. ;

  t = y - map->thumb.halfwidth;
  b = y + map->thumb.halfwidth;
  if (t < topMargin)
    t = topMargin;
  if (b > mapGraphHeight)
    b = mapGraphHeight;

  if (t > y2)
    {
      aa = arrayReCreate (aa, 8, float) ; i = 0 ;
      /* Top green triangle */
      graphColor (GREEN) ;
      array (aa, i++, float) = *offset ;
      array (aa, i++, float) = y1 ;
      array (aa, i++, float) = *offset + .7 ;
      array (aa, i++, float) = y2 ;
      array (aa, i++, float) = *offset - .7 ;
      array (aa, i++, float) = y2 ;
      array (aa, i++, float) = *offset ;
      array (aa, i++, float) = y1 ;
      graphPolygon (aa) ;
      
      /* Thick black line above thumb box */
      graphColor (BLACK) ;
      graphFillRectangle (*offset-0.25, y2,
			  *offset+0.25, t) ;
    }

  y1 = mapGraphHeight - 1. ;
  y2 = mapGraphHeight - 2. ;
  
  i = 0 ;
  if (b < y2)
    {
      aa = arrayReCreate (aa, 8, float) ;
       /* Bottom green triangle */
      graphColor (GREEN) ;
      array (aa, i++, float) = *offset ;
      array (aa, i++, float) = y1 ;
      array (aa, i++, float) = *offset + .7 ;
      array (aa, i++, float) = y2 ;
      array (aa, i++, float) = *offset - .7 ;
      array (aa, i++, float) = y2 ;
      array (aa, i++, float) = *offset ;
      array (aa, i++, float) = y1 ;
      graphPolygon (aa) ;  

      /* Thick black line under thumb box */
      graphColor (BLACK) ;
      graphFillRectangle (*offset-0.25, b,
			  *offset+0.25, y2) ;
    }

  /* now the box-dragable thumb box */
  map->thumb.x = *offset;

  box2 = graphBoxStart();  
  if (b > topMargin && 
      t < mapGraphHeight)
    graphRectangle (*offset - 0.5, t, *offset + 0.5, b) ;

  graphBoxEnd();   /* box2 */
  graphBoxDraw (box2, DARKGREEN, GREEN);

  graphBoxEnd();   /* box1 */
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

  *offset += 1 ;
  
  arrayDestroy (aa) ;

  return;
} /* mapShowLocator */

/*********************************************/
/***************** cursor ********************/

void mapCursorCreate (MAP map, float unit, float x0)
{
  map->cursor.unit = unit ;
  map->cursor.val = x0 / unit ;

  return;
} /* mapCursorCreate */


void mapCursorSet (MAP map, float x)
{
  if (x > 0)
    map->cursor.val = 0.5 + x / map->cursor.unit ;
  else
    map->cursor.val = -0.5 + x / map->cursor.unit ;
  mapCursorShift (map) ;

  return;
} /* mapCursorSet */

static float yCursor ;

void mapCursorDraw (MAP map, float x)
{
  float z = mapCursorPos(map) ;
  float y = MAP2GRAPH(map, z) ;	/* Jean - your x here was a bug */

  if (map->flip)
    z = - z ; 
  if (map->cursor.unit < .011)
    strcpy (map->cursor.text, messprintf ("%.2f",z)) ;
  else
    strcpy (map->cursor.text, messprintf ("%.0f",z)) ;
  map->cursor.box = graphBoxStart() ;
  graphLine (map->thumb.x, y, mapGraphWidth+1, y) ;
  map->cursor.pickBox = graphBoxStart() ;
  graphColor (LIGHTGREEN) ;
  graphFillRectangle (x+0.5, y-0.5, x+2.5+strlen(map->cursor.text), y+0.5) ;
  graphColor (BLACK) ;
  graphTextPtr (map->cursor.text, x+0.5, y-0.5, strlen(map->cursor.text)) ;
  graphBoxEnd () ;
  graphBoxEnd () ;
  graphBoxSetPick (map->cursor.box, FALSE) ; /* only pick on .pickBox */
  graphBoxDraw (map->cursor.box, BLACK, TRANSPARENT) ;
  yCursor = y ;

  return;
} /* mapCursorDraw */


void mapCursorShift (MAP map)
{
  float x1, x2, y1, y2, z = mapCursorPos(map) ;

  if (!map->cursor.box)
    return ;
  graphBoxDim (map->cursor.box, &x1, &y1, &x2, &y2) ;
  if (map->flip)
    z = - z ; 
  if (map->cursor.unit < .011)
    strcpy (map->cursor.text, messprintf ("%.2f",z)) ;
  else
    strcpy (map->cursor.text, messprintf ("%.0f",z)) ;
  y1 = MAP2GRAPH(map, mapCursorPos(map)) ;
  graphBoxShift (map->cursor.box, x1, y1-0.5) ;
  yCursor = y1 ;

  return;
} /* mapCursorShift */


void mapCursorDrag (float *x, float *y, BOOL isDone)
{
  static float oldX ;
  static BOOL isDragging = FALSE ;

  if (isDragging)
    *x = oldX ;
  else
    { oldX = *x ;
      isDragging = TRUE ;
    }

  if (isDone)
    { float newpos ;
      MAP map = currentMap("mapCursorDrag") ;

      newpos = GRAPH2MAP(map, *y+0.5) ;
      mapCursorSet (map, newpos) ;
      isDragging = FALSE ;
    }

  return;
} /* mapCursorDrag */

/***************************************************/

	/* look not used but kept because this is a column display */

void mapShowScale (LOOK look, float *offset)
{
  float unit = 0.01 ;
  float subunit = 0.001 ;
  float x, xx, y ;
  int resolution = 0, max = 0 ;
  char *cp = 0 ;
  MAP map = currentMap("mapShowScale") ;

  mapFindScaleUnit (map, &unit, &subunit) ;
  if (unit >= 1)
    resolution = 0 ;
  else if (unit >= .1)
    resolution = 1 ;
  else if (unit >= .01)
    resolution = 2 ;
  else if (unit >= .001)
    resolution = 3 ;

  x = GRAPH2MAP(map, topMargin+1) ;
  x = subunit * (((x>0)?1:0) + (int)(x/subunit)) ;
  while ((y = MAP2GRAPH(map, x)) < mapGraphHeight - 1)
    { graphLine (*offset+1.0,y,*offset+1.5,y) ;
      x += subunit ;
    }

  x = GRAPH2MAP(map, topMargin+1) ;
  x = unit * (((x>0)?1:0) + (int)(x/unit)) ;
  while ((y = MAP2GRAPH(map, x)) < mapGraphHeight - 1)
    { graphLine (*offset+0.5,y,*offset+1.5,y) ;
      xx = map->flip ? -x : x ;
      cp = messprintf ("%-4.*f", resolution, xx) ;
      graphText (cp, *offset+2, y-0.5) ;
      if (strlen(cp) > max)
	max = strlen(cp) ;
      x += unit ;
    }

  
  graphLine (*offset+1.5, topMargin+1, 
	     *offset+1.5, mapGraphHeight-0.5 ) ;
  if (map->thumb.x)
    { graphLine (map->thumb.x, 
		 MAP2WHOLE(map, map->centre) - map->thumb.halfwidth,
		 *offset+1.5, topMargin+1) ;
      graphLine (map->thumb.x, 
		 MAP2WHOLE(map, map->centre) + map->thumb.halfwidth,
		 *offset+1.5, mapGraphHeight-0.5) ;
    }

  mapCursorDraw (map, *offset) ;

  *offset += 2 + max ;

  return;
} /* mapShowScale */

/*************** mapPrint stuff back from the graphPackage **************/

void mapPrint (void)
{ 
  float min, max, tmp, oldCentre, oldMag ;
  int   origin = 0, nx = 1, ny, imin, imax ;
  void  *fMapLook = 0 ;
  MAP map = currentMap("mapPrint") ; 

  selectedMap = map ;
  graphFitBounds (&nx, &ny) ;

  oldCentre = map->centre ;
  oldMag = map->mag ;

  if (!(fMapLook = fMapGifLook) &&
      !graphAssFind (&MAP2LOOK_ASSOC, &fMapLook) &&
      fMapFindZone(fMapLook, &imin, &imax, &origin))
    { min = imin ; max = imax ; }
  else
    { min = map->min ; 
      max = map->max ; 
      origin = 0 ; 
    }
				/* convert to user visible coords */
  if (map->mag < 0)
    { tmp = min ;
      min = map->max - max + 1 - origin ; 
      max = map->max - min - origin ;
    }
  else
    { min -= origin - 1 ;
      max -= origin ;
    }
				/* ask for bounds */
#ifdef WORM
  if (!messPrompt("Give nx, mag, min, max",
		  messprintf("%d %.2f %.2f %.2f", 
			     nx, mag, min, max), 
		  "ifffz"))
    return ;
  freeint (&nx) ; freefloat (&map->mag) ;
#else
  if (!messPrompt("Please state the zone you wish to print",
	     messprintf("%g   %g", min, max), "ffz"))
    return ;
#endif
  freefloat(&min) ; freefloat(&max) ;
  if (min >= max)
    { map->mag = oldMag ;
      return ;
    }
  				/* convert back */
  if (map->mag < 0)
    { tmp = min ;
      min = map->max - origin - max ;
      max = map->max - origin - tmp + 1 ;
    }
  else
    { min += origin - 1 ;
      max += origin ;
    }

  if (min < map->min)
    min = map->min ;
  if (max > map->max)
    max = map->max ;

  map->centre = (max + min) / 2 ;
 
  graphBoundsPrint (nx + 0.2, 
		    1.05 * (max - min) * (map->mag > 0 ? map->mag : - map->mag) + topMargin + 5,
		    map->draw) ;
  
  map->centre = oldCentre ;
  map->mag = oldMag ;
  (map->draw)() ;

  return;
} /* mapPrint */

/***************************************************************/
/************ code to draw buttons at bottom of page ***********/

static void mapFlipButton (void *arg)
{
  MapColRec *col = (MapColRec*)arg ;
  MAP map = currentMap ("mapFlipButton") ;

  col->isOn = !col->isOn ;
  (map->draw)() ;

  return;
} /* mapFlipButton */


void mapDrawColumnButtons (MAP map)
{
  int i ;
  int box ;
  float h ;
  MapColRec *col ;
  COLOUROPT *buttons ;

  /* use ColouredButtons to get background and to pass arg = MapColRec pointer */
  buttons = (COLOUROPT*) messalloc (sizeof(COLOUROPT)*(arrayMax(map->cols)+1)) ;
  for (i = 0 ; i < arrayMax(map->cols) ; ++i)
    { col = arrp(map->cols, i, MapColRec) ;
      buttons[i].text = col->name ;
      buttons[i].arg = col ;
      buttons[i].f = mapFlipButton ;
      buttons[i].fg = BLACK ;
      buttons[i].bg = col->isOn ? LIGHTGRAY : WHITE ;
    }

  /* draw dummy version off screen to find height */
  box = graphBoxStart () ;
  h = 0 ;
  graphColouredButtons (buttons, -1000, &h, -1000+mapGraphWidth) ;
  graphBoxEnd () ;
  graphBoxDim (box,0,0,0,&h) ;

  if (h >= mapGraphHeight)
    { messout ("Sorry, window not tall enough to show column buttons") ;
      return ;
    }

  /* white out background */
  graphColor (WHITE) ;
  graphFillRectangle (0, mapGraphHeight-h, mapGraphWidth+1, mapGraphHeight+1) ;

  graphColor (BLACK) ;
  h = mapGraphHeight-h  ;
  graphColouredButtons (buttons, 0, &h, mapGraphWidth) ;

  return;
} /* mapDrawColumnButtons */

/************** end of file ***********/
 

#ifdef CODE_NEVER_USED
BOOL mapActive(KEYSET *setp, void** mapp) 
{ 
  if (selectedMap && selectedMap->magic == &MAP_MAGIC
      && selectedMap->activeSet)
    { if (setp) *setp = selectedMap->activeSet (selectedMap->look) ;
      if (mapp) *mapp = selectedMap ;
      return TRUE ;
    }
  selectedMap = 0 ;
  if (setp) *setp = 0 ; 
  if (mapp) *mapp = 0 ;
  return FALSE ;
}
#endif /* CODE_NEVER_USED */

/***************************************/
