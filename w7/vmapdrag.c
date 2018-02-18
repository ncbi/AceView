/*  File: vmapdrag[192z.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  8 11:44 1999 (fw)
 * Created: Fri Apr 23 11:24:45 1993 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: vmapdrag.c,v 1.3 2015/08/17 16:12:02 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "vmap_.h"

/****** geneDragging evaluation *************/
    /* JTM 6 12 91 */
static BOOL dragTopBreakPoint = TRUE, dragLocus ;
static KEY draggedGene = 0 ;
static VerticalMap dragLook = 0 ; /* mhmp 04.12.98 pour un seul Drag ...*/
static SEG dragSeg ;
static GraphFunc oldMiddleDown = 0 , oldMiddleUp = 0, oldMiddleDrag = 0 ;
static float oldDy, oldy, oldx , xBase, yBase ; 
static  int nx, ny , draggedBox ;
static float xCursor ;
#define DRAGFASTLIMIT 6
static Array dragUndoArray = 0 ;

typedef struct { Array segs ;
		 SEG oldseg, newseg ;
	       } DU ;

/************************************************************/

void vMapDragUndo (void)
     /* private to vMapPackage */
{
  int i ;
  DU *dup = 0 ;
  VerticalMap look = currentVerticalMap("dragUndo") ;

  if (!dragUndoArray)
    return ;
  i = arrayMax(dragUndoArray) ;
  while (i--)
    { dup = arrayp(dragUndoArray, i, DU) ;
      if (dup->segs == look->segs)
	break ;
    }
  if (i == -1)
    { messout("No further undo available") ;
      return ;
    }
  arraySort (look->segs, vMapOrder) ;  /* because inclusion of map data may unsort */
  arrayRemove (look->segs, &(dup->newseg), vMapOrder) ;
  arrayInsert (look->segs, &(dup->oldseg), vMapOrder) ;
  vMapDraw (look, dup->newseg.key) ;
  
  for (; i + 1 < arrayMax(dragUndoArray) ; dup++, i++)
    { *dup = *(dup+1) ;
    }
  arrayMax(dragUndoArray)-- ;
}
  

static void vMapMiddleUpGene (double x, double y) 
{
  int dummy ;
  DU* dup ;

   dragSeg.flag |= FLAG_MOVED ;

   if (!dragUndoArray)
     dragUndoArray = arrayCreate(20, DU) ;
   dup = arrayp(dragUndoArray, arrayMax(dragUndoArray), DU) ;

   dup->segs = dragLook->segs ;
   dup->oldseg = dragSeg ;
   
   if (arrayFind (dragLook->segs, &dragSeg, &dummy, vMapOrder))
     messcrash ("vmapmidleup");
   if(!dragLocus)
     { float
	 x =  GRAPH2MAP(dragLook->map, y) ,
	 x1 = dragSeg.x ,
	 x2 = dragSeg.x + dragSeg.dx ;
       if (dragTopBreakPoint)
	 { dragSeg.x = x ;
	   dragSeg.dx = x2 - x ;
	 }
       else
	 { dragSeg.x = x1 ;
	   dragSeg.dx = x - x1 ;
	 }
       if(dragSeg.dx < 0)
	 { dragSeg.dx = - dragSeg.dx ;
	   dragSeg.x -= dragSeg.dx ;
	 }
     }
   else
     { dragSeg.x = GRAPH2MAP(dragLook->map, y) ;
       /* dragLook->map->centre + (y - ny2 - topMargin)  / dragLook->map->mag ; */
       dragSeg.dx = 2 * oldDy  / dragLook->map->mag ;
     }
   
   graphRegister (MIDDLE_DOWN, oldMiddleDown) ;
   graphRegister (MIDDLE_DRAG, oldMiddleDrag) ;
   graphRegister (MIDDLE_UP, oldMiddleUp) ;
   

   arrayInsert (dragLook->segs, &dragSeg, vMapOrder) ;
   dup->newseg = dragSeg ;

   vMapDraw (dragLook, draggedGene) ;

   dragLook = 0 ;
   dragSeg.key = 0 ;
}

/************/

static void vMapMiddleDragGene (double x, double y) 
{  /* erase */
  if(dragLocus)
    { graphXorLine (DRAGFASTLIMIT, oldy , nx, oldy ) ;
      if(oldDy > .1)
	{ graphXorLine (xCursor, oldy + oldDy , xBase, oldy) ;
	  graphXorLine (xCursor, oldy - oldDy , xBase, oldy) ;
	  graphXorLine (DRAGFASTLIMIT, oldy - oldDy , nx, oldy - oldDy) ;
	  graphXorLine (DRAGFASTLIMIT, oldy + oldDy , nx, oldy + oldDy) ;
	}
    }
  else  /* rearrangenmnet */
    { graphXorLine (DRAGFASTLIMIT, oldy - oldDy , nx, oldy - oldDy) ;
      if(oldDy) 
	graphXorLine (DRAGFASTLIMIT, oldy + oldDy , nx, oldy + oldDy) ;
    }
  /* modify coordinates */ 
  oldy = y ;
  if(oldDy && oldx != -99999)
    oldDy *= exp ((x - oldx) / 25.) ; 
  if(oldx != -99999)
    xBase *=  exp ((x - oldx) / 25.) ;
  yBase = y ;
  oldx = x ;
  
  /* redraw */
  if(dragLocus)
    { graphXorLine (DRAGFASTLIMIT, y , nx, y ) ;
      graphBoxShift(draggedBox, xBase, yBase) ;
      
      if(oldDy > .1)
	{ graphXorLine (xCursor, y + oldDy , xBase, y) ;
	  graphXorLine (xCursor, y - oldDy , xBase, y) ;
	  graphXorLine (DRAGFASTLIMIT, y - oldDy , nx, y - oldDy) ;
	  graphXorLine (DRAGFASTLIMIT, y + oldDy , nx, y + oldDy) ;
	}
    }
  else  /* rearrangenmnet */
    { graphXorLine (DRAGFASTLIMIT, y - oldDy , nx, y - oldDy) ;
      if(oldDy) 
	graphXorLine (DRAGFASTLIMIT, y + oldDy , nx, y + oldDy) ;
    }
}
  
/************/

static void vMapMiddleDownRear (double x, double y) 
{
  if (y < oldy)  /* erase and start dragging top line */
    { graphXorLine (DRAGFASTLIMIT, oldy - oldDy , nx, oldy - oldDy) ;
      dragTopBreakPoint = TRUE ;
    }
  else  /* erase and start dragging bottom line */
    { graphXorLine (DRAGFASTLIMIT, oldy + oldDy , nx, oldy + oldDy) ;
      dragTopBreakPoint = FALSE ;
    }
  oldDy = 0 ;
  oldy = y ;
  oldx = x ;
  graphXorLine (DRAGFASTLIMIT, y, nx, y) ;
}

/************/

void vMapDragButton (void)
     /* private to vMapPackage */
{ 
   KEY key ;
   SEG *seg1 ;
   float y , x1, x2 ;
   VerticalMap look = currentVerticalMap("vMapDragButton") ;
 
   if (dragLook) /* mhmp 04.12.98 */
     return ;

   graphFitBounds (&nx, &ny) ;
   seg1 = arrp(look->segs, 
	      arr(look->boxIndex,look->activeBox,int),
	      SEG) ;
   if(!seg1)
     { messout("First pick an object on the map") ;
       return ;
     }
   dragSeg = arr(look->segs, 
		 arr(look->boxIndex,look->activeBox,int),
		 SEG) ;

   
   draggedBox = look->activeBox ;
   graphBoxDim(look->map->cursor.box, &xCursor, &y, & x1, &x2) ;
   graphBoxDim(draggedBox, &xBase, &yBase, & x1, &x2) ;
   key = dragSeg.key ;
   if(dragSeg.flag & FLAG_ANY_LOCUS)
     dragLocus = TRUE ;
   else if (dragSeg.flag & FLAG_ANY_INTERVAL)
     dragLocus = FALSE ;
   else
     { messout("First pick a gene or a locus or an interval") ;
       return ;
     }
   
   dragLook = look ;
   y = MAP2GRAPH(look->map,dragSeg.x) ;
   oldDy = dragSeg.dx * look->map->mag / 2 ;
   oldx = -99999 ;
   
   draggedGene = key ;
   dragLook = look ;
   
   arraySort (look->segs, vMapOrder) ;  /* because inclusion of map data may unsort */
   if (!arrayRemove (look->segs, &dragSeg, vMapOrder))
     messcrash ("Confusion in vMapDragButton");
     
   if(dragLocus)
     { graphXorLine (DRAGFASTLIMIT, y , nx, y ) ;
       if(oldDy > .1)
	 { graphXorLine (xCursor, y + oldDy , xBase, y) ;
	   graphXorLine (xCursor, y - oldDy , xBase, y) ;
	   graphXorLine (DRAGFASTLIMIT, y - oldDy , nx, y - oldDy) ;
	   graphXorLine (DRAGFASTLIMIT, y + oldDy , nx, y + oldDy) ;
	 }
       oldMiddleDown = graphRegister (MIDDLE_DOWN,  vMapMiddleDragGene) ;
     }
   else  /* rearrangenmnet */
     { y += oldDy ;
       graphXorLine (DRAGFASTLIMIT, y - oldDy , nx, y - oldDy) ;
       if(oldDy) 
	 graphXorLine (DRAGFASTLIMIT, y + oldDy , nx, y + oldDy) ;
       oldMiddleDown = graphRegister (MIDDLE_DOWN, vMapMiddleDownRear) ;
     }
  oldy = y ;
   
  oldMiddleDrag = graphRegister (MIDDLE_DRAG, vMapMiddleDragGene) ;
  oldMiddleUp = graphRegister (MIDDLE_UP, vMapMiddleUpGene) ;
}

/***************************************************************/
/*************************** end of file ***********************/
