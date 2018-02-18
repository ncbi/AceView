/*  File: map.h
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: map.h - public interface for the mapcontrol.c package
 *              This is not to be confused with Simon Kelleys
 *              ColControl, which is seperate and mutually exlcusive
 *              with this one.
 *              This mapPackage is used for the fmap display
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 10:41 1998 (fw)
 * * Jul 23 14:53 1998 (edgrif): Remove redeclarations of fmap
 *     - modules must include fmap.h public header instead.
 * Created: Sat Nov 21 17:13:13 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: map.h,v 1.2 2007/03/21 21:03:14 mieg Exp $ */

#ifndef _MAP_H
#define _MAP_H

#include "acedb.h"

/************************************************************/

typedef struct LookStruct *LOOK; /* forward definition
				    This is an incomplete type so far 
				    and needs to be completed. */

typedef void   (*MapColDrawFunc)(LOOK look, float *pos) ;

/************************************************************/


typedef struct MapStruct {
  magic_t *magic ;		/* == &MAP_MAGIC */
  Array cols ;
  float min, max ;
  float mag, centre ;
  float fastScroll ;		/* x coord of locator */
  int *fixFrame ;
  char *activeColName ;
  BOOL flip ;    /* To get the coordinates going upwards */
  struct 
    { int box ;
      float fac, offset, halfwidth, x ;
    } thumb ;
  struct
    { int box, pickBox ;
      int val ;
      float unit ;
      char text[32] ;
    } cursor ;
  VoidRoutine draw ;
  LOOK look ;   /* for call back */
  KEYSET (*activeSet)(LOOK look) ;  /* should return keySet(look->seg->key) */
  AC_HANDLE handle ;
  BOOL isFrame ;		/* for when showing 3 frames */ 
  int zoom_set_G ;  /* disabled ! Holds current zoom setting */		
} *MAP ;

/************************************************************/
                     /* public mapPackage globals */

extern MAP mapGifMap ;	/* used to bypass need for graph in giface ops */
extern int mapGraphWidth, mapGraphHeight ;
extern float halfGraphHeight, topMargin ;

extern magic_t MAP2LOOK_ASSOC ;	/* find LOOK pointer for the map */

/************************************************************/

#define MAP2GRAPH(map,x) \
  (map->mag * ((x) - map->centre) + halfGraphHeight + topMargin)
#define GRAPH2MAP(map,x) \
  (((x) - halfGraphHeight - topMargin) / map->mag + map->centre)
#define MAP2WHOLE(map,x) \
  (3 + topMargin + ((x) - map->thumb.offset)*map->thumb.fac)
#define WHOLE2MAP(map,x) \
  (map->thumb.offset + ((x) - 3 - topMargin)/map->thumb.fac)

/************************************************************/

MAP currentMap(char *caller);	/* find MAP on graph */
MAP mapCreate (VoidRoutine draw) ; /* use mapInsertCol() to add cols now */
void uMapDestroy (MAP z) ;
#define mapDestroy(z) ((z) ? uMapDestroy(z), (z)=0, TRUE : FALSE)
void mapAttachToGraph (MAP map) ; /* needed for fmap create order control */

BOOL mapInsertCol (MAP map, float priority, BOOL init, char *name, 
		   MapColDrawFunc func) ;
	/* if name present already, return FALSE - don't change isOn
	   else insert at priority level
	*/
BOOL mapColSetByName (char *text, int isOn) ;
	/* isOn = -1: return present state, 
           	   0: set off, 
		   1: set on,
		   2: toggle 
	*/
void mapColControl (void) ;
void mapDrawColumns (MAP map) ;
void mapDrawColumnButtons (MAP map) ;


void mapWhole (void) ;
void mapZoomIn (void) ;
void mapZoomOut (void) ;
void mapKbdScroom (int k) ; /* the 4 arrows, pageUp/pageDown */

void mapColMenu (int box) ;
void mapPrint (void) ;

void mapFindScaleUnit (MAP map, float *u, float *sub) ;
void mapShowLocator (LOOK look, float *offset) ;
void mapShowScale (LOOK look, float *offset) ;
void mapThumbDrag (float *x, float *y, BOOL isDone) ;

void mapCursorCreate (MAP map, float unit, float x0) ;
void mapCursorSet (MAP map, float x) ;
void mapCursorDraw (MAP map, float x) ;
void mapCursorShift (MAP map) ;
void mapCursorDrag (float *x, float *y, BOOL isDone) ;
void mapCursorDrawLine (MAP map) ;

void mapNoMag (void) ;          /* prevents rescaling by horizontal mouse motion */
void mapMiddleDown (double x, double y) ;  /* callback for mouse middle down */
#define mapCursorPos(map) ((map)->cursor.val * (map)->cursor.unit)

/*****************************/

#ifdef CODE_NEVER_USED
 /* returns the objects in the active map */
BOOL mapActive(KEYSET *setp, void** mapp)  ;
#endif /* CODE_NEVER_USED */

#endif  /* _MAP_H */

/*************************** eof ****************************/
