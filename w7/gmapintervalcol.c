/*  File: gmapintervalcol.c
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Interval display column for the genetic map
 * Exported functions:
 *     gMapRDIntervalColumn
 *     gMapJTMIntervalColumn 
 *     gMapChromBandsColumn
 * HISTORY:
 * Last edited: Dec 21 14:28 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: gmapintervalcol.c,v 1.2 2015/08/18 23:24:11 mieg Exp $ */

#include "gmap.h"
#include "bump.h"
#include "query.h"

enum intervalSpecialisations { intervalJTM, intervalRD, chromBand };

typedef struct {
  COND query;
  char *text;
  int colour;
} TAGCOLOUR;

typedef struct IntPrivStruct {
  char *query;
  KEYSET keyset;
  BOOL displayNames;
  BOOL displayIntervals;
  BOOL mapShowMultiple;
  BOOL noNeighbours;
  int width;
  int pne, pe, nne, ne;
  Array colours;
  Associator colourLookup;
  KEY symbolTag;
  float namOffset;
  int isNam; /* This column displays objects twice, use these to
		    keep track of which is active */
  float selx, seldx; /* same key can have more than one seg */
  Array segs;
} *INTPRIV;



static void intervalCursorToTop (int box)
{
  COLCONTROL control = currentColControl("intervalCursorToTop");
  MAPCONTROL map = arr(control->boxIndex, box, COLINSTANCE)->map;
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *);
  
  if (map->mag>0)
    mapControlCursorSet (map, seg->x - seg->dx) ;
  else
    mapControlCursorSet (map, seg->x + seg->dx) ;
    
  return;
} /* intervalCursorToTop */


static void intervalCursorToBottom (int box)
{
  COLCONTROL control = currentColControl("intervalCursorToBottom");
  MAPCONTROL map = arr(control->boxIndex, box, COLINSTANCE)->map;
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *);
  
  if (map->mag>0)
    mapControlCursorSet (map, seg->x + seg->dx) ;
  else
    mapControlCursorSet (map, seg->x - seg->dx) ;

  return;
} /* intervalCursorToBottom */

#define fabsf(x) ((x)<0.0 ? -(x) : (x))

static void intervalTopToCursor (int box)
{
  COLCONTROL control = currentColControl("intervalTopToCursor");
  COLINSTANCE instance = arr(control->boxIndex, box, COLINSTANCE);
  MAPCONTROL map = instance->map; 
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *);
  INTPRIV private = instance->private;
  float newTop, length, top, bottom; 
  GMAPSEG save;
  int i;

  save = *seg;
  
  newTop = mapControlCursorPos(map) ;
  if (map->mag>0)
    { top = save.x - save.dx;
      bottom = save.x + save.dx;
    }
  else
    { top = save.x + save.dx;
      bottom = save.x - save.dx;
    } 
  mapControlCursorSet (map, top);
  
  length = (bottom - newTop)/2.0; 
  save.dx = fabsf(length);
  save.x = newTop + length;
  save.flag |= FLAG_MOVED ;
  
  arrayRemove(private->segs, seg, gMapOrder); /* to maintain order */
  arrayInsert(private->segs, &save, gMapOrder); /* to maintain order */
  arrayFind(private->segs, &save, &i, gMapOrder);
  /* must update this, in order for controlSelectBox to set correct values */
  arr(control->boxIndex2, box, void *) = arrp(private->segs, i, GMAPSEG);

  controlSelectBox(box);
  controlDraw();

  return;
} /* intervalTopToCursor */


static void intervalBottomToCursor (int box)
{
  COLCONTROL control = currentColControl("intervalBottomToCursor");
  COLINSTANCE instance = arr(control->boxIndex, box, COLINSTANCE);
  MAPCONTROL map = instance->map; 
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *);
  INTPRIV private = instance->private;
  float newBottom, length, top, bottom; 
  GMAPSEG save;
  int i;

  save = *seg;

  newBottom = mapControlCursorPos(map) ;
  if (map->mag>0)
    { top = save.x - save.dx;
      bottom = save.x + save.dx;
    }
  else
    { top = save.x + save.dx;
      bottom = save.x - save.dx;
    }
  mapControlCursorSet (map, bottom);
  
  length = (newBottom - top)/2.0; 
  save.dx = fabsf(length);
  save.x = top + length; 
  save.flag |= FLAG_MOVED ;

  arrayRemove(private->segs, seg, gMapOrder); /* to maintain order */
  arrayInsert(private->segs, &save, gMapOrder); /* to maintain order */
  arrayFind(private->segs, &save, &i, gMapOrder);
  /* must update this, in order for controlSelectBox to set correct values */
  arr(control->boxIndex2, box, void *) = arrp(private->segs, i, GMAPSEG);

  controlSelectBox(box);
  controlDraw() ;

  return;
} /* intervalBottomToCursor */

static void intervalHighlightSeg(int box)
{
  COLCONTROL control = currentColControl("intervalHighlightSeg");
  MAPCONTROL map = arr(control->boxIndex, box, COLINSTANCE)->map;
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *);

  gMapHighlightKey(map, seg->key);

  return;
} /* intervalHighlightSeg */

MENUOPT gMapIntervalMenu[] = {
  { (VoidRoutine)intervalHighlightSeg, "Highlight" },
  { (VoidRoutine)intervalCursorToTop, "Cursor to top" },
  { (VoidRoutine)intervalCursorToBottom, "Cursor to bottom" },
  { (VoidRoutine)intervalTopToCursor, "Move top to cursor" },
  { (VoidRoutine)intervalBottomToCursor, "Move bottom to cursor" },
  { 0, 0 }
} ;

static void japanese (char *text, float x, float y)
{
  int n = strlen(text) ;

  graphTextUp (text, x, y + 0.65*n) ;
#ifdef OLDCODE
  a[1] = 0 ;
  for (i = 0 ; i < n ; ++i)
    { a[0] = *text++ ;
      graphText (a, x, y) ;
      y += 0.65 ;
    }
#endif
}


static void rearrLine (int flag, float x, float y1, float y2)
{
  if (flag & FLAG_DEFICIENCY)
    graphLine (x, y1, x, y2) ;
  else if (flag & FLAG_DUPLICATION)
    { graphLine (x-0.2, y1, x-0.2, y2) ;
      graphLine (x+0.2, y1, x+0.2, y2) ;
      graphColor (GRAY) ;
      graphFillRectangle (x-0.2, y1, x+0.2, y2) ;
      graphColor (BLACK) ;
    }
  else
    graphFillRectangle (x-0.15, y1, x+0.15, y2) ;
}

static BOOL intervalIsNam(COLINSTANCE instance, int box)
{ INTPRIV private = instance->private;
  float x1, x2, y1, y2;

  graphBoxDim(box, &x1, &y1, &x2, &y2);
  
  return private->namOffset < x1;
}

static BOOL intervalSetSelect(COLINSTANCE instance, 
			      int box, double x, double y)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *); 
  GeneticMap look = (GeneticMap)map->look;
  INTPRIV private = instance->private;

  private->isNam = intervalIsNam(instance, box);
  private->selx = seg->x;
  private->seldx = seg->dx;
  
  gMapUnselect(map);
  look->selectKey = seg->key;
  look->friendsInfoValid = TRUE;
  look->neighboursInfoValid = TRUE;
  look->x1 = seg->x - seg->dx;
  look->x2 = seg->x + seg->dx;
  (void)gMapNeighbours(map, &look->neighbours, seg->key);
  (void)gMapPositive(map, &look->positive, seg->key);
  (void)gMapNegative(map, &look->negative, seg->key);
  control->activeKey = seg->key; /* inter-map communication */  
    
  return FALSE; /* never need to redraw */
}

static BOOL intervalUnselect(COLINSTANCE instance, int box)
{ 
  GeneticMap look = (GeneticMap)(instance->map->look);

  *look->messageText = 0;

  if (look->messageBox)
    graphBoxDraw(look->messageBox, -1, -1);
  
  return FALSE;
  
}

static void intervalFollowBox(COLINSTANCE instance, 
			      int box, double x, double y)
{
  INTPRIV private = instance->private;
  GMAPSEG *seg = (GMAPSEG *) arr(instance->map->control->boxIndex2, 
				 box, void *); 
  BOOL isMap = (class(seg->key) == _VMap);

  if (isMap)
    { if (!private->mapShowMultiple ||
	  !gMapFollowMap (instance->map->control, instance->map, seg->key))
	{ displayPreserve();
	  display(seg->key, instance->map->key, 0);
	}
    }
  else
    display(seg->key, instance->map->key, TREE);
}

static void intervalDoColour(COLINSTANCE instance, int box)
{ 
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *); 
  GeneticMap look = (GeneticMap)map->look;
  INTPRIV private = instance->private;
  KEY key;
  int flags;

  if (!seg) 
    { graphBoxDraw(box, BLACK, TRANSPARENT); /* no seg, just refresh circle */
      return;
    }
  
  key = seg->key;
  
  if (key == look->selectKey && 
      instance == control->activeInstance &&
      private->selx == seg->x && private->seldx == seg->dx &&
      intervalIsNam(instance, box) == private->isNam)
    /* do highlight */
    { graphBoxDraw(box, BLACK, CYAN);
      strncpy(look->messageText, name(seg->key), 100);
      strcat(look->messageText, messprintf(" %.2f", seg->x - seg->dx));
      strcat(look->messageText, messprintf(" %.2f", seg->x + seg->dx));
	
      if (look->messageBox) 
	graphBoxDraw(look->messageBox, -1, -1);
      control->activeBox = box;
    }
  else
    { flags = gMapOverlap(map, key, seg->x - seg->dx, seg->x + seg->dx);
      if (keySetFind(look->highlight, key, 0))
	graphBoxDraw(box, BLACK, MAGENTA);
      else if (flags & FLAG_STRESSED)
	graphBoxDraw(box, graphContrast(private->pe), private->pe);
      else if (flags & FLAG_ANTI_STRESSED)
	graphBoxDraw(box, graphContrast(private->ne), private->ne);
      else if (flags & FLAG_ANTI_RELATED)
	graphBoxDraw(box, graphContrast(private->nne), private->nne);
      else if (flags & FLAG_RELATED)
	 graphBoxDraw(box, graphContrast(private->pne), private->pne);
      else if (!private->noNeighbours && gMapIsNeighbour(map, key))
	graphBoxDraw(box, BLACK, PALECYAN);   
      else if (key == look->selectKey) /* still show self */
	graphBoxDraw(box, BLACK, PALECYAN);   
      else
	{ void *col;
	  if (private->colours && 
	      !intervalIsNam(instance, box) &&
	      assFind(private->colourLookup, assVoid(box), &col))
	    graphBoxDraw(box, BLACK, assInt(col));
	  else
	    graphBoxDraw (box, BLACK, WHITE) ;
	}
    }
}

static void RDIntervalDraw(COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)map->look;
  BUMP  bump = bumpCreate (control->graphWidth, 0) ;
  float y1, y2, z1, z2, x, n ;
  int   i, ibox, ix, cmax = 0 ;
  GMAPSEG   *seg ;
  INTPRIV private = instance->private;
  int firstSeg, lastSeg, countDir;
  KEY symb;

  if (map->mag > 0)
    { firstSeg = 0;
      lastSeg = arrayMax(private->segs);
      countDir = 1;
    }
  else
    { firstSeg = arrayMax(private->segs)-1 ;
      lastSeg = -1;
      countDir = -1;
    }
  
  *offset += 1.25 ;

  graphTextHeight (0.8) ;
  private->namOffset = 0.0; /* not needed, set to incocuous value */
 
  for (i = firstSeg ; i != lastSeg ; i += countDir)
    { seg = arrp(private->segs,i,GMAPSEG) ;
      y1 = MAP2GRAPH(map,seg->x - seg->dx) ;
      y2 = MAP2GRAPH(map,seg->x + seg->dx) ;
      if (y1>y2)
	{ float tmp = y1; y1 = y2; y2 = tmp; }
      
      if (!(seg->flag & FLAG_ANY_INTERVAL))
	continue;
      if (y2 < control->topMargin || y1 > control->graphHeight-1)
	continue;
      if (private->keyset && !keySetFind(private->keyset, seg->key, 0))
	continue ;
      if (keySetFind(look->hidden, seg->key, 0))
	continue ;

      symb = gMapKeyOnTag(seg->key, private->symbolTag);
      n = private->displayNames ? 0.65*strlen(name(symb)) : 0;
      ix = 0 ; 
      if (y1 < control->topMargin) y1 = control->topMargin ;
      if (y2 > control->graphHeight) y2 = control->graphHeight ;
      if (y2 > y1+4+n)
	{ z1 = y1 ; z2 = y2 ; }
      else if (y1 > control->topMargin+n)
	{ z1 = y1-n ; z2 = y2 ; }
      else
	{ z1 = y1 ; z2 = y2+n ; }
      bumpItem (bump,1,(z2-z1+1)+0.2,&ix,&z1) ;
      if (private->width && ix >= private->width)
	ix = private->width;
      
      if (cmax < ix)
	cmax = ix;
       
      remarkRegister(map, seg->key, seg->x);
      controlRegBox(instance, ibox = graphBoxStart(), seg);
      graphBoxInfo (ibox, seg->key, 0); /* WWW */
      
      x = *offset + ix ;
      if (y1 > control->topMargin) 
	graphLine (x-0.25, y1, x+0.25, y1) ; 
      if (y2 < control->graphHeight)
	graphLine (x-0.25, y2, x+0.25, y2) ; 
      if (private->displayNames)
	{
	  if (y2 > y1+4+n)
	    { japanese (name(symb), x-0.5, y1+2) ;
	      rearrLine (seg->flag, x, y1, y1+2) ;
	      rearrLine (seg->flag, x, y1+2+n, y2) ;
	    }
	  else if (y1 > control->topMargin+n)
	    { japanese (name(symb), x-0.5, y1-n) ;
	      rearrLine (seg->flag, x, y1, y2) ;
	    }
	  else 
	    { rearrLine (seg->flag, x, y1, y2) ;
	      japanese (name(symb), x-0.5, y2) ;
	    }
	}
      graphBoxEnd () ;
      if (map->submenus)
	graphBoxMenu (ibox, gMapIntervalMenu) ; 
      
      if(seg->key == control->from)
	control->fromBox = ibox;
    }
  
  graphTextFormat (PLAIN_FORMAT) ;
  graphTextHeight (0.0) ;
  
  *offset += cmax+1 ;
  bumpDestroy (bump) ;
}

static void JTMIntervalDraw (COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)map->look;
  BUMP  bump = bumpCreate (control->graphWidth, 0) ;
  float y1, y2, x, n;
  int   i, ibox, ix, ix1, nmax = 0, cmax = 0;
  GMAPSEG   *seg ;
  INTPRIV private = instance->private;
  void *v ;
  Associator anyIntAss = assCreate() ;
  float anyIntervalOffset = *offset += 1.25 ;
  int firstSeg, lastSeg, countDir;
  KEY symb;

  if (map->mag > 0)
    { firstSeg = 0;
      lastSeg = arrayMax(private->segs);
      countDir = 1;
    }
  else
    { firstSeg = arrayMax(private->segs) -1;
      lastSeg = -1;
      countDir = -1;
    }

  graphTextHeight (0.8) ;

  for (i = firstSeg ; i !=  lastSeg ; i += countDir)
    { seg = arrp(private->segs,i,GMAPSEG) ;
      y1 = MAP2GRAPH(map,seg->x - seg->dx) ;
      y2 = MAP2GRAPH(map,seg->x + seg->dx) ;
      if (y1>y2)
	{ float tmp = y1; y1 = y2; y2 = tmp; }
      
      if (!(seg->flag & FLAG_ANY_INTERVAL))
	continue;
      if (y2 < control->topMargin || y1 > control->graphHeight-1)
	continue;
      if (private->keyset &&!keySetFind(private->keyset, seg->key, 0))
	continue ;
      if (keySetFind(look->hidden, seg->key, 0))
	continue ;

      ix = 0 ; 
      if (y1 < control->topMargin) y1 = control->topMargin ;
      if (y2 > control->graphHeight) y2 = control->graphHeight ;
      bumpItem (bump,1,(y2-y1+1)+0.2,&ix,&y1) ;
      if (private->width && ix >= private->width)
	ix = private->width-1;

      if (cmax < ix)
	cmax = ix;
       
      remarkRegister(map, seg->key, seg->x);
      controlRegBox(instance, ibox = graphBoxStart(), seg);
      graphBoxInfo (ibox, seg->key, 0) ; /* WWW */
      
      x = *offset + ix ;
      if (y1 > control->topMargin) 
	graphLine (x - 0.25, y1, x + 0.25, y1) ; 
      else
	{ graphLine (x - 0.25, control->topMargin + .5,
		     x, control->topMargin) ; 
	  graphLine (x + 0.25, control->topMargin + .5, 
		     x, control->topMargin) ; 
	}
      if (y2 < control->graphHeight)
	graphLine (x - 0.25, y2, x+0.25, y2) ; 
      else
	{ graphLine (x - 0.25, control->graphHeight - .5, 
		     x, control->graphHeight) ; 
	  graphLine (x + 0.25, control->graphHeight - .5, 
		     x, control-> graphHeight) ; 
	}
      rearrLine (seg->flag, x, y1, y2) ;
      graphBoxEnd () ;
      graphBoxDraw (ibox, BLACK, TRANSPARENT) ;
      
      v = (char*)(&anyIntAss) + ix ;
      assInsert(anyIntAss, seg, v) ;
      
      if (seg->key == control->from)
	control->fromBox = ibox;
      /* if we're doing names, they overwrite frombox and get priority */
      
      if (map->submenus)
	graphBoxMenu (ibox, gMapIntervalMenu) ;
    }
  
  *offset += cmax+1 ;
  bumpDestroy (bump) ;
  bump = bumpCreate (1, 0) ;

  private->namOffset = *offset+1.0;
  *offset += 1.25 ;

  graphTextHeight (0.8) ;

  if (private->displayNames)
    for (i = firstSeg ; i != lastSeg ; i += countDir)
      { seg = arrp(private->segs,i,GMAPSEG) ;
	y1 = MAP2GRAPH(map,seg->x - seg->dx) ;
	y2 = MAP2GRAPH(map,seg->x + seg->dx) ;
	if (y1>y2)
	  { float tmp = y1; y1 = y2; y2 = tmp; }

	if (!(seg->flag & FLAG_ANY_INTERVAL))
	  continue;
	if (y2 < control->topMargin || y1 > control->graphHeight-1)
	  continue;
	if (private->keyset && !keySetFind(private->keyset, seg->key, 0))
	  continue;
	if (keySetFind(look->hidden, seg->key, 0))
	  continue ;

	ibox = graphBoxStart();
	graphBoxInfo (ibox, seg->key, 0) ; /* WWW */

	array(control->boxIndex, ibox, COLINSTANCE) = instance;
	array(control->boxIndex2, ibox, void *) = (void *)seg;
	symb = gMapKeyOnTag(seg->key, private->symbolTag);
	n = 0.8*strlen(name(symb)) ; 
	if (n > nmax)
	  nmax = n ;
	ix = 0 ; 
	if (y1 < control->topMargin) y1 = control->topMargin ;
	if (y2 > control->graphHeight) y2 = control->graphHeight ;
	bumpTest (bump,1,1,&ix,&y1) ;
	x = *offset + ix ;
	if(y1 < y2 + 3 &&
	   y1 < control->graphHeight-1)
	  { bumpRegister (bump,1,1,&ix,&y1) ;
	    graphText(name(symb), x, y1) ; 
	  }
	graphBoxEnd () ;
	
	if (map->submenus)
	  graphBoxMenu (ibox, gMapIntervalMenu) ;
	  
	if (seg->key == control->from)
	  control->fromBox = ibox;
	
	if (anyIntAss && assFind (anyIntAss, seg, &v))
	  { ix1 = (char*)v - (char*)(&anyIntAss) ;
	    y1 += .5;
	    if (y1 < control->graphHeight-1 &&
		y1 < y2)
	      { ibox = graphBoxStart();
		array(control->boxIndex, ibox, COLINSTANCE) = instance;
		array(control->boxIndex2, ibox, void *) = 0;
		graphCircle(anyIntervalOffset + ix1, y1, .4) ;
		graphBoxEnd();
		graphBoxSetPick(ibox, FALSE);
	      }
	  }
      }
  

  *offset += nmax + 1 ;
  bumpDestroy (bump) ;
  assDestroy(anyIntAss);

  graphTextFormat (PLAIN_FORMAT) ;
  graphTextHeight (0.0) ;

  
}

static void chromBandsDraw (COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)map->look;
  BUMP  bump = bumpCreate (control->graphWidth, 0) ;
  float y1, y2, x, n;
  int   i, j, ibox, ix, nmax = 0, cmax = 0;
  GMAPSEG   *seg ;
  INTPRIV private = instance->private;
  int firstSeg, lastSeg, countDir;
  KEY symb;
  OBJ obj;
  
  private->colourLookup = assReCreate(private->colourLookup);

  if (map->mag > 0)
    { firstSeg = 0;
      lastSeg = arrayMax(private->segs);
      countDir = 1;
    }
  else
    { firstSeg = arrayMax(private->segs) -1;
      lastSeg = -1;
      countDir = -1;
    }

  graphTextHeight (0.8) ;
  if(private->displayIntervals)
    { for (i = firstSeg ; i !=  lastSeg ; i += countDir)
      { seg = arrp(private->segs,i,GMAPSEG) ;
        y1 = MAP2GRAPH(map,seg->x - seg->dx) ;
	y2 = MAP2GRAPH(map,seg->x + seg->dx) ;
	if (y1>y2)
	  { float tmp = y1; y1 = y2; y2 = tmp; }
	
	if (!(seg->flag & FLAG_ANY_INTERVAL))
	  continue;
	if (y2 < control->topMargin || y1 > control->graphHeight-1)
	  continue;
	if (private->keyset &&!keySetFind(private->keyset, seg->key, 0))
	  continue ;
	if (keySetFind(look->hidden, seg->key, 0))
	continue ;
	
	ix = 0 ; 
	if (y1 < control->topMargin) y1 = control->topMargin ;
	if (y2 > control->graphHeight) y2 = control->graphHeight ;
	bumpItem (bump,1,(y2-y1),&ix,&y1) ;
	if (private->width && ix >= private->width)
	  ix = private->width-1;
	
	if (cmax < ix)
	  cmax = ix;
	
	remarkRegister(map, seg->key, seg->x);
	controlRegBox(instance, ibox = graphBoxStart(), seg);
	graphBoxInfo (ibox, seg->key, 0) ; /* WWW */
						
	x = *offset + (ix * 1.5) ;
	if (y1 > control->topMargin) 
	  graphLine (x, y1, x+1, y1) ; 
	if (y2 < control->graphHeight)
	  graphLine (x, y2, x+1, y2) ; 
	if (seg->flag & FLAG_CENTROMERE)
	  { graphLine (x, y1, x+1, y2) ;
	  graphLine (x+1, y1, x, y2) ;
	}
	else
	  { graphLine (x, y1, x, y2) ;
	  graphLine (x+1, y1, x+1, y2) ;
	  }
	graphBoxEnd () ;
	
	obj = 0;
	for (j=0; j<arrayMax(private->colours); j++)
	  { COND q = arr(private->colours, j, TAGCOLOUR).query;
	  if (q && queryFind3(q, &obj, seg->key))
	    { assInsert(private->colourLookup, assVoid(ibox),
			assVoid(arr(private->colours, j, TAGCOLOUR).colour));
	    break;
	    }
	  }
	bsDestroy(obj);
	
	if (seg->key == control->from)
	  control->fromBox = ibox;
	/* if we're doing names, they overwrite frombox and get priority */
	
	if (map->submenus)
	  graphBoxMenu (ibox, gMapIntervalMenu) ;
      }
  
    *offset += (cmax * 1.5) + 1;
    bumpDestroy (bump) ;
    } /* end if display Intervals */
  bump = bumpCreate (1, 0) ;

  private->namOffset = *offset+0.8;
  *offset += 1.0;

  graphTextHeight (0.8) ;

  if (private->displayNames)
    for (i = firstSeg ; i != lastSeg ; i += countDir)
      { seg = arrp(private->segs,i,GMAPSEG) ;
	y1 = MAP2GRAPH(map,seg->x - seg->dx) ;
	y2 = MAP2GRAPH(map,seg->x + seg->dx) ;
	if (y1>y2)
	  { float tmp = y1; y1 = y2; y2 = tmp; }

	if (!(seg->flag & FLAG_ANY_INTERVAL))
	  continue;
	if (y2 < control->topMargin || y1 > control->graphHeight-1)
	  continue;
	if (private->keyset && !keySetFind(private->keyset, seg->key, 0))
	  continue;
	if (keySetFind(look->hidden, seg->key, 0))
	  continue ;

	symb = gMapKeyOnTag(seg->key, private->symbolTag);
	ix = 0 ; 
	if (y1 < control->topMargin) y1 = control->topMargin ;
	if (y2 > control->graphHeight) y2 = control->graphHeight ;
	bumpTest (bump,1,1,&ix,&y1) ;
	x = *offset + ix ;
	if(y1 > y2 + 3 ||
	   y1 > control->graphHeight-1)
	  continue ;
	n = 0.8*strlen(name(symb)) ; 
	if (n > nmax)
	  nmax = n ;

	ibox = graphBoxStart();
	graphBoxInfo (ibox, seg->key, 0) ; /* WWW */

	array(control->boxIndex, ibox, COLINSTANCE) = instance;
	array(control->boxIndex2, ibox, void *) = (void *)seg;
	
	bumpRegister (bump,1,1,&ix,&y1) ;
	graphText(name(symb), x, y1) ; 

	graphBoxEnd () ;
	
	if (map->submenus)
	  graphBoxMenu (ibox, gMapIntervalMenu) ;
	  
	if (seg->key == control->from)
	  control->fromBox = ibox;
      }
  

  *offset += nmax + 1 ;
  bumpDestroy (bump) ;

  graphTextFormat (PLAIN_FORMAT) ;
  graphTextHeight (0.0) ;

}

static void intPrivDestroy(void *p)
{ INTPRIV private = (INTPRIV)p;
  int i;
  
  if (private->keyset) 
    keySetDestroy(private->keyset);

  if (private->colours)
    { for (i=0; i<arrayMax(private->colours); i++)
	condDestroy(arr(private->colours, i, TAGCOLOUR).query);
      arrayDestroy(private->colours);
    }
}

static void buildQueries(INTPRIV private)
{ int i;
  char *t;
  COND c;
  
  if (!private->colours)
    return;
  
  for(i=0; i<arrayMax(private->colours); i++)
    { condDestroy(array(private->colours, i, TAGCOLOUR).query); /* OK zero */
      if (!(t = array(private->colours, i, TAGCOLOUR).text) || 
	  !condConstruct(t, &c))
	c = 0;
      array(private->colours, i, TAGCOLOUR).query = c;
    }
}

#define NLINES 15
struct configLocals {
  char query[280],tag[30];
  int colWidth;
  BOOL displayNames,showIntervals,showMultMap,displayIntervals;
  int colour[4];
  char *temp[NLINES];
  int colour2[NLINES];
};

static BOOL intervalConfigure(COLINSTANCE instance)
{ INTPRIV private = instance->private;
  int i;
  float line = 15.0,height;

  struct  configLocals *cf = (struct configLocals *) messalloc(sizeof(struct configLocals));

  if(private->colours)
    height = 0.8;
  else
    height = 0.3;
  if(controlCreateConfig(instance,cf,"Configure interval column",0.8,height)){
    
  /*initialise data*/
    if(private->query)
      strcpy(cf->query,private->query);
    cf->colWidth = private->width;
    if(private->symbolTag)
      sprintf(cf->tag,"%s",name(private->symbolTag));
    if(private->colours){
      for(i=0; i!=NLINES; i++){
	cf->temp[i] = (char *) messalloc(sizeof(char)*80);
	if(i < arrayMax(private->colours)){
	  strcpy(cf->temp[i],arrp(private->colours, i, TAGCOLOUR)->text);
	  cf->colour2[i] = arrp(private->colours, i, TAGCOLOUR)->colour;
	}
	else
	  cf->colour2[i] = WHITE; 
      }
    }

    cf->displayIntervals = private->displayIntervals;
    cf->displayNames = private->displayNames;
    cf->showIntervals = !private->noNeighbours;
    cf->showMultMap =private->mapShowMultiple; 
    cf->colour[0] = private->pne;
    cf->colour[1] = private->pe;
    cf->colour[2] = private->nne;
    cf->colour[3] = private->ne;

  /* draw the configuration options */

    graphTextEditor("Query for display:",cf->query,280,4,2,0);
    graphIntEditor("Column width Restriction:",&cf->colWidth,4.0,4.0,0);
    graphTextEditor("Symbol Tag:",cf->tag,30,40.0,4.0,0);
    
    if(instance->draw == chromBandsDraw)         /* only display for chrom band */
      graphToggleEditor("Display interval boxes",&cf->displayIntervals,4.0,9.0);
    graphToggleEditor("Display names",&cf->displayNames,4.0,10.0);
    graphToggleEditor("Show intervals as neighbours",&cf->showIntervals,4.0,11.0);
    graphToggleEditor("Show multiple maps",&cf->showMultMap,4.0,12.0);
    
    graphColourEditor("Positive data, no error."," ",&cf->colour[0],36.0,10.0);
    graphColourEditor("Positive data, error in co-ordinates."," ",&cf->colour[1],36.0,11.1);
    graphColourEditor("Negative data, no error."," ",&cf->colour[2],36.0,12.2);
    graphColourEditor("Negative data, error in coordinates."," ",&cf->colour[3],36.0,13.3);
    
    if(private->colours){
      graphText("Queries:", 1, line);
      graphText("Colours:", 76, line++);
      for(i = 0; i!=NLINES; i++)
	{
	  graphTextEditor(" ",cf->temp[i],70,2.0,line,0);
	  graphColourEditor(" "," ",&cf->colour2[i],77.0,line);
	  line++;line++;
	}
    }
  }
  graphRedraw();
  return FALSE;
}

static void gmapintervalFinal(COLINSTANCE instance, void *locals, BOOL ok)
{ struct configLocals *cf = locals;
  INTPRIV private = instance->private;
  int i,n;

  if(ok){

    private->displayNames =cf->displayNames;
    private->displayIntervals =cf->displayIntervals;
    private->mapShowMultiple =cf->showMultMap;
    private->noNeighbours = !cf->showIntervals;

    private->pne =cf->colour[0];
    private->pe= cf->colour[1];
    private->nne = cf->colour[2];
    private->ne = cf->colour[3];
    private->width = cf->colWidth;

    if (private->keyset)
      keySetDestroy(private->keyset);
    if (private->query)
      messfree(private->query);

    if (strlen(cf->query) != 0)
      { private->query = handleAlloc(0, 
				     instance->handle, 
				     1+strlen(cf->query));
	strcpy(private->query, cf->query);
	private->keyset = queryKey(instance->map->key, private->query);
      }
    else
      { private->keyset = 0;
	private->query = 0;
      }
    
    if (strlen(cf->tag) != 0)
      lexaddkey(cf->tag, &private->symbolTag, 0);
    else
      private->symbolTag = 0;

    if (private->colours)
      { /* first delete old, cond will be destroyed by buildqueries */
	char *s;
	for (i=0; i<arrayMax(private->colours); i++)
	  if ((s = array(private->colours, i, TAGCOLOUR).text))
	    { messfree(s);
	      array(private->colours, i, TAGCOLOUR).text = 0;
	      array(private->colours, i, TAGCOLOUR).colour = WHITE;
	    }
	for (n = 0, i=0; i!=NLINES; i++)
	  { if (strlen(cf->temp[i]) != 0)
	      { s = handleAlloc(0, instance->handle, 
				1+strlen (cf->temp[i]));
		strcpy(s, cf->temp[i]);
		array(private->colours, n, TAGCOLOUR).text = s;
		array(private->colours, n, TAGCOLOUR).colour = cf->colour2[i];
		n++;
	      }
	  }
	buildQueries(private);
      }
  }
  else{
    if(private->colours){
      for (i=0; i!=NLINES; i++)
	messfree(cf->temp[i]);
    }
    messfree(cf);
  } 
}

static void intervalSave(COLINSTANCE instance, OBJ obj)
{
  INTPRIV private = instance->private;
  int i;

  bsPushObj(obj); /* into Interval_col_conf */

  if (private->query)
    bsAddData(obj, str2tag("Query"), _Text, private->query);
  
  if (private->displayNames)
    bsAddTag(obj, str2tag("Names_on"));

  if (!private->displayIntervals && instance->draw == chromBandsDraw)
    bsAddTag(obj, str2tag("No_interval_boxes"));

  if (private->noNeighbours)
    bsAddTag(obj, str2tag("No_neighbours"));

  if (private->mapShowMultiple)
    bsAddTag(obj, str2tag("Show_multiple"));

  
  if (private->width)
    bsAddData(obj, str2tag("Width"), _Int, &private->width);
  
  if (private->colours)
    for (i=0; i<arrayMax(private->colours); i++)
      { char *s = array(private->colours, i, TAGCOLOUR).text;
	if (s)
	  bsAddData(obj, str2tag("Colours"), _Text, s);
	controlSetColour(obj, arr(private->colours, i, TAGCOLOUR).colour);
      }

  bsAddTag(obj, str2tag("Pne")); controlSetColour(obj, private->pne);
  bsAddTag(obj, str2tag("Pe")); controlSetColour(obj, private->pe);
  bsAddTag(obj, str2tag("Nne")); controlSetColour(obj, private->nne);
  bsAddTag(obj, str2tag("Ne")); controlSetColour(obj, private->ne);
  
  if (private->symbolTag)
    bsAddData(obj, str2tag("Symbol"), _Text, name(private->symbolTag));  
}

      
	  
static BOOL intervalCreate(COLINSTANCE instance, OBJ init)
{ 
  char *s1;
  int w;
  INTPRIV private;

  private = (INTPRIV)halloc(sizeof(struct IntPrivStruct), instance->handle);
  blockSetFinalise (private, intPrivDestroy);

  switch (assInt(instance->proto->specialisation))
    {
    case intervalRD:
      instance->draw = RDIntervalDraw;
      private->colours = 0;
      break;
    case intervalJTM:
      instance->draw = JTMIntervalDraw;
      private->colours = 0;
      break;
    case chromBand: 
      instance->draw = chromBandsDraw;
      /* private->colours is needed to free the conds in intPrivDestroy */
      /* to avoid premature freeing, we free it there */
      private->colours = arrayCreate(10, TAGCOLOUR);
      private->colourLookup = assHandleCreate(instance->handle);
      break;
    } 
  instance->setSelectBox = intervalSetSelect; 
  instance->unSelectBox = intervalUnselect; 
  instance->doColour = intervalDoColour;
  instance->followBox = intervalFollowBox;
  instance->configure = intervalConfigure;
  instance->save = intervalSave;
  instance->private = private;
  instance->configFinal = gmapintervalFinal;
  private->keyset = 0; 
  private->displayNames = FALSE;
  private->displayIntervals = TRUE;
  private->mapShowMultiple = FALSE;
  private->noNeighbours = FALSE;
  private->width = 0;
  private->query = 0;
  private->symbolTag = 0;
  private->pne = GREEN;
  private->pe = RED;
  private->nne = LIGHTBLUE;
  private->ne = RED;
  private->segs = *(instance->proto->convertResults);

  if (init)
    { bsPushObj(init); /* into intervalColCOnf */
      if (bsGetData(init, str2tag("Query"), _Text, &s1))
	{ private->query = handleAlloc(0, instance->handle, 1+strlen(s1));
	  strcpy(private->query, s1);
	  private->keyset = queryKey(instance->map->key, private->query);
	}
      
      private->displayNames = bsFindTag(init, str2tag("Names_on"));
      private->noNeighbours = bsFindTag(init, str2tag("No_neighbours"));
      private->mapShowMultiple = bsFindTag(init, str2tag("Show_multiple"));

      if (bsGetData(init, str2tag("Width"), _Int, &w))
	private->width = w;
      
      if (private->colours && bsGetData(init, str2tag("Colours"), _Text, &s1))
	{ do 
	    { int n = arrayMax(private->colours);
	      char *new = (char *)handleAlloc(0, instance->handle,
					      1+strlen(s1));
	      strcpy(new, s1);
	      array(private->colours, n, TAGCOLOUR).text = new;
	      array(private->colours, n, TAGCOLOUR).query = 0;
	      array(private->colours, n, TAGCOLOUR).colour = 
		controlGetColour(init);
	    } while (bsGetData(init, _bsDown, _Text, &s1));
	  buildQueries(private);
	}

      if (bsFindTag(init, str2tag("Pne")))
	private->pne = controlGetColour(init);
      if (bsFindTag(init, str2tag("Pe")))
	private->pe = controlGetColour(init);
      if (bsFindTag(init, str2tag("Nne")))
	private->nne = controlGetColour(init);
      if (bsFindTag(init, str2tag("Ne")))
	private->ne = controlGetColour(init);

      if (bsGetData (init, str2tag("Symbol"), _Text, &s1))
	lexaddkey(s1, &private->symbolTag, 0);

      if(instance->draw == chromBandsDraw)
	private->displayIntervals = !(bsFindTag(init, str2tag("No_interval_boxes")));
    }
  else
    { /* No init, use default name and leave private->keyset zero */
      private->displayNames = TRUE;
    }
  
  return TRUE;
}

struct ProtoStruct gMapRDIntervalColumn = {
  0,
  intervalCreate,
  0,
  "Interval_RD",
  (void *) intervalRD,
  FALSE,
  gMapConvert,
  0,
  "The interval_RD column can display any interval object (ie one with "
  "Left and Right or Multi_ends tags in the #Map_postion). The set of such "
  "objects in a map is found by scanning all objects two to the right of the "
  "Contains tag. This "
  "set is then filtered using the query in the configuration to give the set "
  "of objects to draw. \n\n"
  "Intervals are drawn by default as bold vertical "
  "lines. If the tag Deletetion if found a narrow line is used instead; if "
  "the tag Duplication is found a wide gray line is used. If \"Display "
  "names\" is selected, the names of objects are written, vertically, on "
  "objects.\n\n"
  "The relationship of interval objects to other objects on the map, "
  "described "
  "by the Positive and Negative tags, is shown using four colours that can "
  "be picked in the configuration panel. \"Positive data, no error\" is used "
  "when an interval is recorded as containing a point object and their "
  "repspective co-ordinates agree with this. If the co-ordinates do not "
  "agree, the colour for \"Positive data, error in co-ordinates\" is used. "
  "Similarly, if an interval is recorded as not containing a point, one of " 
  "the "
  "Negative data colours is used, depending on the co-ordinates. "

};

struct ProtoStruct gMapJTMIntervalColumn = {
  0,
  intervalCreate,
  0,
  "Interval_JTM",
  (void *) intervalJTM,
  FALSE,
  gMapConvert,
  0,
  "The interval_JTM column can display any interval object (ie one with "
  "Left and Right or Multi_ends tags in the #Map_postion). The set of such "
  "objects in a map is found by scanning all objects two to the right of the "
  "Contains tag. This "
  "set is then filtered using the query in the configuration to give the set "
  "of objects to draw. \n\n"
  "Intervals are drawn as vertical "
  "lines. If \"Display "
  "names\" is selected, the names of objects are written, horizontally, next "
  "to the lines. If the name of an object can be written within the limits "
  "of the interval it represents, a circle is drawn on the vertical line at "
  "the position of the name, to indicate this. Picking a vertical line will "
  "colour the corresponding name in light green, and vice versa.\n\n"
  "The relationship of interval objects to other objects on the map, "
  "described "
  "by the Positive and Negative tags, is shown using four colours that can "
  "be picked in the configuration panel. \"Positive data, no error\" is used "
  "when an interval is recorded as containing a point object and their "
  "repspective co-ordinates agree with this. If the co-ordinates do not "
  "agree, the colour for \"Positive data, error in co-ordinates\" is used. "
  "Similarly, if an interval is recorded as not containing a point, one of " 
  "the "
  "Negative data colours is used, depending on the co-ordinates. "
};

struct ProtoStruct gMapChromBandsColumn = {
  0,
  intervalCreate,
  0,
  "Interval_SRK",
  (void *) chromBand,
  FALSE,
  gMapConvert,
  0,
 "The interval_SRK column can display any interval object (ie one with "
  "Left and Right or Multi_ends tags in the #Map_postion). The set of such "
  "objects in a map is found by scanning all objects two to the right of the "
  "Contains tag. This "
  "set is then filtered using the query in the configuration to give the set "
  "of objects to draw. \n\n"
  "Intervals are drawn as coloured vertical "
  "bars. If \"Display "
  "names\" is selected, the names of objects are written, horizontally, next "
  "to the lines. Picking a vertical bar will "
  "colour the corresponding name in light green, and vice versa.\n\n"
  "The relationship of interval objects to other objects on the map, "
  "described "
  "by the Positive and Negative tags, is shown using four colours that can "
  "be picked in the configuration panel. \"Positive data, no error\" is used "
  "when an interval is recorded as containing a point object and their "
  "repspective co-ordinates agree with this. If the co-ordinates do not "
  "agree, the colour for \"Positive data, error in co-ordinates\" is used. "
  "Similarly, if an interval is recorded as not containing a point, one of " 
  "the "
  "Negative data colours is used, depending on the co-ordinates.\n\n"
  "The colours of the bars are controlled by looking for the presence "
  "of certain tags in the objects they represent. What these tags are and the "
  "colours they produce is controlled by the column's configuration, which "
  "has a table of tags and colours. To change a tag click on it to turn it "
  "yellow, then type; deleteing the tag name will remove that line from the "
  "table. Choose colours either by clicking with the left mouse button, to "
  "cycle through colours, or using a menu on the colour chooser box. "
  "If an object has more than one of the tags in the table, the one highest "
  "up the table takes priority when determining its colour."
};



















 
 
 
