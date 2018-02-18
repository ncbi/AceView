/*  File: cmapdisp.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: physical unit display of a chromosome
 	based on gmap
 * Exported functions: cmapDisplay()
 * HISTORY:
 * Last edited: Nov 19 13:13 1998 (fw)
 * Created: Thu Jan  9 22:54:23 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: cmapdisp.c,v 1.4 2015/08/18 23:24:11 mieg Exp $ */

#include "acedb.h"
#include "keyset.h"
#include "graph.h"
#include "key.h"
#include "lex.h"
#include "bs.h"
#include "a.h"
#include "systags.h"
#include "classes.h"
#include "sysclass.h"
#include "tags.h"
#include "session.h" 

#include "disptype.h"
#include "display.h"
#include "query.h"
#include "grid.h"
#include "bump.h"

/* the basic units are Mb - left end is 0 - calculate length from contigs
*/

/*******************************************/

typedef struct {
  float bin ;
  float wfac ;			/* width factor */
  int colour ;
  int offset ;
  Array coords, ctgs, seqlen, count ;
  Graph graph ;
  char name[32] ;
  int isLine ;
} COLUMN ;

typedef struct LOOKSTUFF
  { int   magic ;	/* == MAGIC */
    KEY   key ;		/* a chromosome */
    int   activeBox ;
    unsigned int flag ;
    float centre, mag ;	
    float max ;
    Array segs,      	/* info on things to be drawn */
          boxes,	/* a BOX structure of each box */
          neighbours ;	/* 0-terminated lists of neighbours */
    Associator ctgAss ;
    Array columns ;		/* of COLUMN* */
    AC_HANDLE handle ;
  } *LOOK ;

static int MAGIC = 361718 ;

typedef struct
  { KEY key ;
    float x1, x2 ;
    float p1, p2 ;
    int flag ;
  } SEG ;

typedef struct
  { SEG* seg ;		/* pointer OK: segs fixed when drawing */
    int  neigh ;	/* index in look->neighbours */
    int  col ;
  } BOX ;

/* p1, p2 are used for contigs: cmap <-> pmap
                   for in_situ data for in_situ clones
*/

				/* seg->flag flags */
#define FLAG_HIGHLIGHT		0x0001
#define FLAG_PROCESSED		0x0002
#define FLAG_INSITU		0x0004
#define FLAG_HYB		0x0008
#define FLAG_DELETING		0x0010

#define LOOKGET(name)     LOOK look ; \
                          if (!graphAssFind (&MAGIC, &look)) \
		            messcrash ("graph not found in %s",name) ; \
                          if (look->magic != MAGIC) \
                            messcrash ("%s received a wrong pointer",name)

static void cMapDestroy (void) ;
static void cMapRecalculate(void) ;
static void cMapPick (int box, double x , double y) ;
static void cMapDragCursor (float *x, float *y, BOOL isDone) ;
static void cMapMiddleDown (double x, double y) ;
static void cMapDrag (double x, double y) ;
static void cMapUp (double x, double y) ;
static void cMapDraw (LOOK look, KEY key) ;
static BOOL cMapConvert (LOOK look, BOOL isForce) ;
static void cMapResize (void) ;
static void cMapSelect (LOOK look, int box) ;
static void cMapFollow (LOOK look, double x, double y) ;
static void cMapHighlight (void) ;
static void cMapAdd (void) ;
static void cMapSubtract (void) ;
static void cMapWhole (void) ;
static void cMapFullGenome (void) ;
static void cMapPartGenome (void) ;
static void cMapChangeSymbolSize (void) ;
static void cMapSetSloppyBump (void) ;
static void cMapSetMag (void) ;
static void cMapColFromKeySet (void) ;
static void cMapColFromFile (void) ;
static void cMapColEdit (void) ;
static void cMapColDraw (LOOK look, COLUMN *col) ;
static void cMapHardColumns (void) ;

static BOOL isCoverage = FALSE ;
static void changeCoverage (void) { isCoverage = 1 - isCoverage ; }

static MENUOPT cMapMenu[] =
              { {graphDestroy, "Quit"},
		{help,"Help"},
		{graphPrint,"Print"},
		{displayPreserve,"Preserve"},
		{cMapRecalculate,"Recalculate"},
		{cMapHighlight,"Highlight Selected Objects"},
		{cMapAdd,"Add Selected Objects"},
		{cMapSubtract,"Subtract Selected Objects"},
		{cMapFullGenome,"Full Genome Distribution"},
		{cMapPartGenome, "II, III, X Distribution"},
		{cMapChangeSymbolSize,"Change Symbol Size"},
		{cMapSetSloppyBump,"Change bump sloppiness"},
		{changeCoverage,"Coverage <-> bump"},
		{cMapSetMag,"Set magnification"},
		{cMapColFromKeySet, "Create histogram from keyset"},
		{cMapColFromFile, "Create histogram from file"},
		{cMapColEdit, "Edit histogram columns"},
	        {cMapHardColumns, "Cheat - direct column set"},
		 {0, 0}
	      } ;

static int	cursorBox ;

/**********************************/

     /* Recalculate all cMaps */

void cMapMakeAll(void)
{
  KEY chromosome = 0, cacheKey ;	/* 0 primes lexNext() */
  struct LOOKSTUFF lookStuff ;
  KEY key = 0, _VcMap;

  lexword2key("cMap", &key, _VMainClasses);
  _VcMap = KEYKEY(key);
  if (!key) return ;
  

  while (lexNext (_VMap, &chromosome))
    { lookStuff.key = chromosome ;
      lookStuff.segs = 0 ;
      if (lexReClass (chromosome, &cacheKey, _VcMap)) /* force recalculation */
	arrayKill (cacheKey) ;
      if (cMapConvert (&lookStuff, TRUE))
	arrayDestroy (lookStuff.segs) ;
      if (messIsInterruptCalled ())
	break ;
    }
}

/**********************/

BOOL cMapDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  AC_HANDLE handle = handleCreate () ;
  LOOK look=(LOOK)halloc(sizeof(struct LOOKSTUFF), handle) ;

 /* I hope either key or from is a chromo */
  if (key && class(key) != _VMap)
    { 
      KEY tmp = from ;
      from = key ;
      if (class(tmp) == _VMap)
	key = tmp ;
      else
	key = 0 ;
    }

/* if not try to get it using from */

  if (!key && from)
    { OBJ obj = bsCreate(from) ;
      if (obj)
	{ 
	  bsGetKey(obj, _Map, &key) ;
	  bsDestroy(obj) ;
	}
    }

  if (!key)
    { if (from)
	display (from, 0, TREE) ;
      goto abort ;
    }

  look->handle = handle ;
  look->activeBox = 0 ;
  look->magic = MAGIC;
  look->flag = 0 ;
  look->max = 15.0 ;		/* arbitrary min */
  look->key = key ;
  look->segs = arrayHandleCreate (128, SEG, handle) ;
  look->boxes = arrayHandleCreate (64, BOX, handle) ;
  look->neighbours = arrayHandleCreate (32, int, handle) ;
  look->columns = arrayHandleCreate (8, COLUMN*, handle) ;
  look->ctgAss = assHandleCreate (handle) ;

  if (!cMapConvert (look, FALSE) ||
      !look->max)
    goto abort ;
  
  if (isOldGraph)
    { graphRetitle (name (key)) ;
      cMapDestroy () ;
      graphAssRemove (&MAGIC) ;
    }
  else 
    { if (!displayCreate(CMAP))
	goto abort ;
      graphRetitle (name(key)) ;
      
      graphRegister (RESIZE,(GraphFunc)cMapResize) ;
      graphRegister (DESTROY, cMapDestroy) ;
      graphRegister (PICK, (GraphFunc)cMapPick) ;
      graphRegister (MIDDLE_DOWN, (GraphFunc)cMapMiddleDown) ;
      graphRegister (MIDDLE_DRAG, (GraphFunc) cMapDrag) ;
      graphRegister (MIDDLE_UP, (GraphFunc) cMapUp) ;
      graphMenu (cMapMenu) ;
    }

  graphAssociate (&MAGIC, look) ;

  cMapWhole () ;	/* sets ->centre, ->mag and calls Draw() */
  return TRUE ;

abort :
  messfree (look) ;
  return FALSE ;
}

/*********************/

static void cMapDestroy (void)
{
  LOOKGET("cMapDestroy") ;

  graphAssRemove (&MAGIC) ;
  look->magic = 0 ;
  messfree (look->handle) ;
}

/**********************************/

static void cMapResize (void)
{
  LOOKGET("cMapResize") ;

  cMapDraw (look, 0) ;
}

/***********************************/

static void cMapPick (int box, double x, double y) 
{
  LOOKGET("cMapPick") ;

  if (!box)
    return ;

  if (box == cursorBox)
    graphBoxDrag (cursorBox, cMapDragCursor) ;
  else if (box == look->activeBox)
    cMapFollow (look, x, y) ;
  else
    cMapSelect (look, box) ;
}

/*********************/

static void cMapClear (void)
{
  KEY curr ;
  LOOKGET("cMapClear") ;
  
  if (look->activeBox && arrp(look->boxes,look->activeBox,BOX)->seg)
    curr = arrp(look->boxes,look->activeBox,BOX)->seg->key ;
  else 
    curr = 0 ;
  cMapConvert (look, FALSE) ;
  cMapDraw (look, curr) ;
}

static void cMapRecalculate (void)
{
  KEY cacheKey ;
  KEY key = 0, _VcMap;
  LOOKGET("cMapRecalculate") ;

  lexaddkey("cMap", &key, _VMainClasses);
  _VcMap = KEYKEY(key);
  if (!key) return ;

  if (lexReClass (look->key, &cacheKey, _VcMap))
    arrayKill (cacheKey) ;

  cMapClear () ;
}

static void cMapSelect (LOOK look, int box) 
{
  int i, j ;

  if (look->activeBox)
    { graphBoxDraw (look->activeBox, DARKBLUE, LIGHTBLUE) ;
      if ((i = arrp(look->boxes, box, BOX)->neigh))
	while ((j = arr(look->neighbours, i--, int)))
	  graphBoxDraw (j, DARKBLUE, LIGHTBLUE) ;
    }
  if (arrp(look->boxes, box, BOX)->seg)
    { look->activeBox = box ;
      graphBoxDraw (box, DARKBLUE, RED) ;
      if ((i = arrp(look->boxes, box, BOX)->neigh))
	while ((j = arr(look->neighbours, i--, int)))
	  graphBoxDraw (j, DARKBLUE, LIGHTRED) ;
    }
  else
    look->activeBox = 0 ;
}

static void cMapFollow (LOOK look, double x, double y)
{
  display (arrp(look->boxes, look->activeBox, BOX)->seg->key,
		 look->key, 0) ;
}

/**************************************************/
/**************** drawing info ********************/

static int	nx, ny ;	/* window dimensions */
static float	yCentre ;	/* ny/2 */
static float	yLength ;	/* length of picture */
static float	topMargin = 2 ;	/* space at top for buttons etc */
static float	bottomMargin = 1 ; /* space at bottom */
static float	xCursor = 5 ;	/* midline of mini-chromosome */
static float	xScale = 10 ;	/* scale bar text LHS */
static float	xContig = 14 ;	/* midline of contig boxes */
static int	xItem = 16 ;	/* left edge of item field */
static int	xInSitu = 30 ;	/* for insitu clones */

static float	symbolSize = 1.0 ; /* size of little squares */

static BOOL isFullGenome = FALSE ;

#define MAP2GRAPH(look,x) \
  (yCentre  +  look->mag * ((x) - look->centre))
#define GRAPH2MAP(look,x) \
  (((x) - yCentre) / look->mag + look->centre)

/********* start off with some utility routines ********/

static void getNxNy(void)
{
  graphFitBounds (&nx, &ny) ;
  yLength = (ny - topMargin - bottomMargin) ;
  yCentre = topMargin + 0.5*yLength ;
}

/***************************************/

static void cMapWhole (void)
{ 
  LOOKGET("cMapWhole") ; 
  getNxNy() ;
  look->centre = 0.5 * look->max ;
  look->mag = yLength / look->max ;
  cMapDraw (look, 0) ;
}

static void cMapZoomIn (void)
{ 
  LOOKGET("cMapZoomIn") ; 
  look->mag *= 2 ; 
  cMapDraw (look, 0) ;
}

static void cMapZoomOut (void)
{ 
  LOOKGET("cMapZoomOut") ; 
  look->mag /= 2 ; 
  cMapDraw (look, 0) ;
}

static void cMapChangeSymbolSize (void)
{
  if (messPrompt ("Change size of little square symbols to",
		   messprintf ("%f", symbolSize), "fz"))
    freefloat (&symbolSize) ;
}

static void cMapSetMag (void)
{
  LOOKGET("cMapZoomIn") ; 
  if (graphPrompt ("Give new magnification ", messprintf ("%.2f", look->mag), "fz"))
    freefloat (&look->mag) ;
  cMapDraw (look, 0) ;
}

/**************************************************************/
/***************** dragging code - middle button **************/

static double	oldy, oldDy, oldx;
static BOOL	dragFast ;
#define DRAGFASTLIMIT xScale

static void cMapDrag (double x, double y) 
{
  if (dragFast)
    { graphXorLine (0, oldy - oldDy, DRAGFASTLIMIT, oldy - oldDy) ;
      graphXorLine (0, oldy + oldDy, DRAGFASTLIMIT, oldy + oldDy) ;
    }
  else
    graphXorLine (DRAGFASTLIMIT, oldy, nx, oldy) ;

  oldy = y ;

  if (dragFast)
    { oldDy *= exp ((x - oldx) / 25.) ;
      oldx = x ;
      graphXorLine (0, y - oldDy, DRAGFASTLIMIT, y - oldDy) ;
      graphXorLine (0, y + oldDy, DRAGFASTLIMIT, y + oldDy) ;
    }
  else
    graphXorLine (DRAGFASTLIMIT, y, nx, y) ;
}

static void cMapUp (double x, double y) 
{ 
  float x1,x2,y1,y2 ;
  LOOKGET("cMapUp") ;

  if (dragFast)
    { graphBoxDim (cursorBox, &x1, &y1, &x2, &y2) ;
      look->mag *= (y2 - y1) / (2. * oldDy) ;
      look->centre = look->max * (y - topMargin) / yLength ;
    }
  else
    look->centre +=  (y - 0.5 - yCentre) / look->mag ;
  cMapDraw (look, 0) ;
}

static void cMapMiddleDown (double x, double y) 
{  
  float x1,x2,y1,y2 ;
  LOOKGET("cMapMiddleDown") ;

  getNxNy () ; 

  graphBoxDim (cursorBox, &x1, &y1, &x2, &y2) ;
  oldDy = (y2 - y1) / 2. ;

  dragFast = (x < DRAGFASTLIMIT) ? TRUE : FALSE ;

  if(dragFast)
    { graphXorLine (0, y - oldDy, DRAGFASTLIMIT, y - oldDy) ;
      graphXorLine (0, y + oldDy, DRAGFASTLIMIT, y + oldDy) ;
    }
  else
    graphXorLine (DRAGFASTLIMIT, y, nx, y) ;
   
  oldx = x ;
  oldy = y ;
  graphRegister (MIDDLE_DRAG, (GraphFunc) cMapDrag) ;
  graphRegister (MIDDLE_UP, (GraphFunc) cMapUp) ;
}

static void cMapDragCursor (float *x, float *y, BOOL isDone)
{
  if (isDone)
    { float x1,y1,x2,y2 ;
      LOOKGET("cMapDragCursor") ;

      getNxNy() ;
      graphBoxDim (cursorBox, &x1, &y1, &x2, &y2) ;
      look->centre = look->max * 
	((*y + 0.5*(y2-y1)) - topMargin) / yLength ;
      cMapDraw (look, 0) ;
    }
  else
    *x = xCursor - 0.5 ;
}

static float cmapSloppyness = 0 ;
static void cMapSetSloppyBump (void)
{
  if (messPrompt ("Set sloppiness of bump system - 0 for nosloppy",
		   messprintf ("%f", cmapSloppyness), "fz"))
    freefloat (&cmapSloppyness) ;
}

/**************************************************/
/**************************************************/

static void drawScale (LOOK look)
{
  float cutoff = 5 / look->mag ;
  float unit = 0.01 ;
  float subunit = 0.001 ;
  float x, y, start, end, oldxScale = xScale ;

  while (unit < cutoff)
    { unit *= 2 ;
      subunit *= 5 ;
      if (unit >= cutoff)
	break ;
      unit *= 2.5 ;
      if (unit >= cutoff)
	break ;
      unit *= 2 ;
      subunit *= 2 ;
    }

  if (isFullGenome)
    { graphTextHeight (0.8) ;
      xScale = xContig-0.5 ;
    }

  start = GRAPH2MAP(look, topMargin) ;
  if (start < 0)
    { start = 0 ;
/*      graphLine (xScale-1.5, MAP2GRAPH(look,start), 
		 xScale-0.5, MAP2GRAPH(look,start)) ; */
    }
  end = GRAPH2MAP(look, ny-bottomMargin) ;
  if (end > look->max)
    { end = look->max ;
/*      graphLine (xScale-1.5, MAP2GRAPH(look,end),
		 xScale-0.5, MAP2GRAPH(look,end)) ; */
    }
      
  x = unit * (((start > 0) ? 1 : 0) + (int)(start/unit)) ;
  while (x <= end)
    { y = MAP2GRAPH(look, x) ;
      if (isFullGenome)
	{
	  graphLine (xScale-1.5,y,xScale-0.5,y) ;
	  if (unit >= 1)
	    graphText (messprintf ("%4.0f",x),xScale-4.2,y-.3) ;
	  else if (unit >= .1)
	    graphText (messprintf ("%4.1f",x),xScale-4.2,y-.3) ;
	  else if (unit >= .01)
	    graphText (messprintf ("%4.2f",x),xScale-4.2,y-.3) ;
	  else if (unit >= .001)
	    graphText (messprintf ("%4.3f",x),xScale-4.2,y-.3) ;
	}
      else
	{
	  graphLine (xScale-1.5,y,xScale-0.5,y) ;
	  if(unit >= 1)
	    graphText (messprintf ("%-4.0f",x),xScale,y-.5) ;
	  else if(unit >= .1)
	    graphText (messprintf ("%-4.1f",x),xScale,y-.5) ;
	  else if(unit >= .01)
	    graphText (messprintf ("%-4.2f",x),xScale,y-.5) ;
	  else if(unit >= .001)
	    graphText (messprintf ("%-4.3f",x),xScale,y-.5) ;
	}
      x += unit ;
    }

  if (!isFullGenome)
    { x = subunit * (((start>=0)?1:0) + (int)(start/subunit)) ;
      while (x <= end)
	{ y = MAP2GRAPH(look,x) ;
	  graphLine (xScale-1.0,y,xScale-0.5,y) ;
	  x += subunit ;
	}
    }

  graphLine (xScale-0.5, MAP2GRAPH(look,start), 
	     xScale-0.5, MAP2GRAPH(look,end)) ;

  if (isFullGenome)
    { xScale = oldxScale ;
      graphTextHeight (0.0) ;
    }
}

static void drawChromosome (LOOK look)
{
  int i ;

#define MAP2CHROM(x) \
  (topMargin + yLength * (x) / look->max)

  graphFillRectangle (xCursor - 0.25, topMargin, 
		      xCursor + 0.25, ny - bottomMargin) ;

  cursorBox = graphBoxStart() ;
  arrayp(look->boxes,cursorBox,BOX)->seg = 0 ;
  graphRectangle (xCursor - 0.5, 
		  MAP2CHROM(GRAPH2MAP(look, topMargin)), 
		  xCursor + 0.5, 
		  MAP2CHROM(GRAPH2MAP(look, ny - bottomMargin))) ;
  graphBoxEnd () ;
  graphBoxDraw (cursorBox,DARKGREEN,GREEN) ;

  graphColor (DARKGRAY) ;
  graphLine (xCursor, MAP2CHROM(GRAPH2MAP(look,topMargin)), 
	     xScale-0.5, topMargin ) ;
  graphLine (xCursor, MAP2CHROM(GRAPH2MAP(look,ny - bottomMargin)),
	     xScale-0.5, ny - bottomMargin) ;
  graphColor (BLACK) ;

  graphTextHeight (0.75) ;
  for (i = 0 ; i <= 10 ; ++i)
    graphText (messprintf ("%3d%%",10*i),
	       1, topMargin + yLength * i / 10.0 - 0.25) ;
  graphTextHeight (0) ;
}

/***********************************************/

static void cMapDrawBoxFinish (int ibox, BOX *box)
{
  graphBoxEnd () ;
  if (box->seg && box->seg->flag & FLAG_HIGHLIGHT)
    graphBoxDraw (ibox, BLACK, MAGENTA) ;
  else
    graphBoxDraw (ibox, BLACK, box->col) ;
}

static void cMapAddNeigh (LOOK look, BOX *box, int ibox)
{
  if (!box->neigh)
    { array(look->neighbours, arrayMax(look->neighbours), int) = 0 ;
      box->neigh = arrayMax(look->neighbours) ;
    }
  array(look->neighbours, arrayMax(look->neighbours), int) = ibox ;
}

/***********************************************/

static MENUOPT buttonOpts[] = {
  {cMapWhole, "Whole"}, 
  {cMapZoomIn, "Zoom In"},
  {cMapZoomOut, "Zoom Out"},
  {cMapClear, "Clear"},
  {0, 0}} ;

static void cMapDraw (LOOK look, KEY curr)
{
  int iseg, ibox ;
  SEG *seg ;
  float screenMin, screenMax ;
  int x, xs ;
  float y, y1s, y2s ;
  BUMP keyBump, inSituBump ;
  BOX *box ;

  look->neighbours = arrayReCreate (look->neighbours, 32, int) ;
  look->boxes = arrayReCreate (look->boxes, 64, BOX) ;
  look->activeBox = 0 ;

  if (!isFullGenome)
    { graphClear () ;
      graphColor (BLACK) ;
      getNxNy () ;
      if (yLength < 6)
	{ messout ("Sorry, this window is too small for a chromo map") ;
	  return ;
	}
      drawScale (look) ;
      drawChromosome (look) ;
      graphText (name(look->key), 1, 0.7) ;
    }
  

  screenMin = GRAPH2MAP(look, topMargin) ;
  screenMax = GRAPH2MAP(look, ny - bottomMargin) ;

  keyBump = bumpCreate ((nx-xItem)/symbolSize, 20) ;
  bumpSetSloppy(keyBump, cmapSloppyness) ;
  inSituBump = bumpCreate (nx-xInSitu, 20) ;
  bumpSetSloppy(inSituBump, cmapSloppyness) ;

  for (iseg = 0 ; iseg < arrayMax(look->segs) ; ++iseg)
    { seg = arrp(look->segs, iseg, SEG) ;
      if (seg->x2 < screenMin || seg->x1 > screenMax)
	continue ;
      ibox = graphBoxStart () ;
      box = arrayp(look->boxes, ibox, BOX) ;
      box->seg = seg ;
      box->neigh = 0 ;
      if (class(seg->key) == _VContig)
	{
	  graphRectangle (xContig - 0.5, MAP2GRAPH(look, seg->x1), 
			  xContig + 0.5, MAP2GRAPH(look, seg->x2)) ;
	  box->col = YELLOW ;
	  cMapDrawBoxFinish (ibox, box) ;
	}
      
      else if (class(seg->key) == _VClone || class(seg->key) == _VGene ||
	     class(seg->key) == _VSequence)
	{
	  x = 0 ;
	  y = MAP2GRAPH(look,0.5*(seg->x1+seg->x2)) ;
	  if (isFullGenome)
	    { bumpItem (keyBump, 1, 0.04, &x, &y) ;
	      graphLine (xItem + x, y, xItem + x + 1, y) ;
	    }
	  else if (isCoverage)
	    graphFillRectangle (xItem, MAP2GRAPH(look,seg->x1),
				xItem+symbolSize, MAP2GRAPH(look,seg->x2)) ;
	  else
	    { bumpItem (keyBump, 1, symbolSize * 0.75, &x, &y) ;
	      graphRectangle (xItem + x*symbolSize, y, 
			      xItem + (x+0.8)*symbolSize, y + 0.6*symbolSize) ;
	    }
	  box->col = LIGHTBLUE ;
	  cMapDrawBoxFinish (ibox, box) ;
	  if (!isFullGenome && (seg->flag & FLAG_INSITU))
	    { xs = 0 ;
	      y1s = MAP2GRAPH(look, look->max * seg->p1) ;
	      y2s = MAP2GRAPH(look, look->max * seg->p2) ;
	      bumpItem (inSituBump, 1, y2s - y1s + 1.0, &xs, &y1s) ;
	      xs += xInSitu ;
	      graphLine (xs, y1s, x + xItem + 0.8, y) ;
	      graphLine (xs, y2s, x + xItem + 0.8, y + 0.6) ;
	      ibox = graphBoxStart() ;
	      cMapAddNeigh (look, box, ibox) ; /* to old box */
	      box = arrayp(look->boxes, ibox, BOX) ;
	      box->seg = seg ;
	      box->neigh = 0 ;
	      cMapAddNeigh (look, box, ibox-1) ;
	      graphRectangle (xs, y1s, xs + 0.5, y2s) ;
	      box->col = LIGHTBLUE ;
	      cMapDrawBoxFinish (ibox, box) ;
	    }
	}
      else
	graphBoxEnd() ;
    }
  
  bumpDestroy (keyBump) ;
  bumpDestroy (inSituBump) ;

				/* draw histogram columns */
  { int i ;
    for (i = 0 ; i < arrayMax (look->columns) ; ++i)
      cMapColDraw (look, arr(look->columns, i, COLUMN*)) ;
  }

  if (!isFullGenome)
    graphButtons (buttonOpts, 5, 0.5, nx) ; /* at end of mapDraw() */
  else
    drawScale (look) ;

  graphRedraw () ;
}

/*********************************************************/
/********************* conversion code *******************/

#define PMAP_FAC  (1830.0 / 1000000.0)

static int cMapOrder (const void *a, const void *b) /* for arraySort() call */
{
  const SEG *seg1 = (const SEG*)a, *seg2 = (const SEG*)b ;
  float diff = seg1->x1 - seg2->x1 ;

  if (diff > 0)
    return 1 ;
  else if (diff < 0)
    return -1 ;
  else
    return 0 ;
}

/********************/

typedef struct
  { int pmap ;
    float vmap ;
    SEG *seg ;
    KEY key ;
  } BLIT ;

static int blitOrder (const void *a1, const void *a2) /* for arraySort() call */
{
  const BLIT *b1 = (const BLIT*)a1, *b2 = (const BLIT*)a2 ;
  float diff = b1->vmap - b2->vmap ;

  if (diff > 0) return 1 ;
  else if (diff < 0) return -1 ;
  else return 0 ;
}

static BOOL cMapConvert (LOOK look, BOOL isForce)
{
  KEY contig, clone, key, _VcMap;
  OBJ Chrom, Contig, Clone ;
  int i, min, max, x1, x2 ;
  float vaxmap, total = 0.0, end ;
  SEG *seg ;
  int iseg = 0 ;
  BLIT *blit, *blitOld ;
  int iblit = 0 ;
  static Array blits = 0 ;

  if (!(Chrom = bsCreate (look->key)))
    return FALSE ;

  lexword2key("cMap", &key, _VMainClasses);
  _VcMap = KEYKEY(key);

  if (!isForce && lexReClass (look->key, &key, _VcMap))
    { arrayDestroy (look->segs) ;
      if ((look->segs = arrayGet (key, SEG, "kkkkfi")))
	{ end = 0.0 ;
	  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
	    { seg = arrp(look->segs, i, SEG) ;
	      if (seg->x2 > end) end = seg->x2 ;
	    }
	  look->max = end ;
	  return TRUE ;
	}
    }
  look->segs = arrayReCreate (look->segs, 64, SEG) ;
  blits = arrayReCreate (blits,64,BLIT) ;

  if (bsGetKey (Chrom, _Contig, &contig)) do
    { if (!(Contig = bsCreate (contig)))
	continue ;
      seg = arrayp (look->segs, iseg, SEG) ;
      seg->key = contig ;
      seg->flag = 0 ;
      min = 1000000 ;
      max = -1000000 ;
      if (bsGetKey (Contig, _Clone, &clone)) do
	{ if (!(Clone = bsCreate (clone)))
	    continue ;
	  if (bsGetKey (Clone, _pMap, &key) &&
	      bsGetData (Clone, _bsRight, _Int, &x1) &&
	      bsGetData (Clone, _bsRight, _Int, &x2))
	    { if (x1 < min) min = x1 ;
	      if (x2 > max) max = x2 ;
	      if (bsGetData (Clone, _Vaxmap, _Float, &vaxmap))
		{ blit = arrayp (blits, iblit, BLIT) ;
		  blit->pmap = 0.5 * (x1 + x2) ;
		  blit->vmap = vaxmap ;
		  blit->seg = seg ;
		  blit->key = clone ;
		  ++iblit ;
		}
	    }
	  bsDestroy (Clone) ;
	} while (bsGetKey (Contig, _bsDown, &clone)) ;
      seg->p1 = min ;
      seg->p2 = max ;
      seg->x1 = 0.0 ;
      if (max > min)
	seg->x2 = (max-min) * PMAP_FAC ;
      else
	seg->x2 = 0.0 ;
      total += seg->x2 ;
      bsDestroy (Contig) ;
      ++iseg ;
    } while (bsGetKey (Chrom, _bsDown, &contig)) ;

  bsDestroy (Chrom) ;

  arraySort (blits, blitOrder) ;
  blitOld = 0 ;
  end = 0.0 ;
  for (iblit = 0 ; iblit < arrayMax(blits) ; ++iblit)
    { blit = arrp(blits, iblit, BLIT) ;
      seg = blit->seg ;
      if ((!blitOld || blitOld->seg != seg) &&
	  !(seg->flag & FLAG_PROCESSED))
        { end += 0.1 ;
	  seg->x1 = end ;
	  end += seg->x2 ;
	  seg->x2 = end ;
	  end += 0.1 ;
	  seg->flag |= FLAG_PROCESSED ;
	}
      blitOld = blit ;
    }
  look->max = end ;
  
  arraySort (look->segs, cMapOrder) ;
  lexaddkey (name(look->key), &key, _VcMap) ;
  arrayStore (key, look->segs, "kkkkfi") ;
  return TRUE ;
}

static LOOK aLook ;

static SEG* cMapAddSeg (KEY key, int ictg, float x1, float x2)
{
  SEG *seg = arrayp (aLook->segs, arrayMax (aLook->segs), SEG) ;
  SEG *ctgSeg = arrp (aLook->segs, ictg, SEG) ;
  
  seg->key = key ;
  seg->x1 = ctgSeg->x1 + PMAP_FAC*(x1 - ctgSeg->p1) ;
  seg->x2 = ctgSeg->x1 + PMAP_FAC*(x2 - ctgSeg->p1) ;
  seg->flag = 0 ;

  return seg ;
}

static void cMapAdd (void)	/* add keys from active keyset */
{ 
  int i, k, x1, x2 ;
  int iseg;
  void *ictg;
  SEG *seg ;
  KEY key, contig, clone ;
  OBJ obj, Clone ;
  KEYSET kSet = 0 ;
  static Array map = 0 ;
  static Array units = 0 ;
  GRIDMAP *grid ;
 
  LOOKGET("cMapAdd") ;

  if (!keySetActive(&kSet, 0))
    { messout("First select a keySet window, thank you.") ;
      return ;
    }

  map = arrayReCreate (map, 8, GRIDMAP) ;
  look->ctgAss = assReCreate (look->ctgAss) ;

  for (iseg = 0 ; iseg < arrayMax(look->segs) ; ++iseg)
    { seg = arrp(look->segs, iseg, SEG) ;
      if (class(seg->key) == _VContig)
	assInsert (look->ctgAss, assVoid(seg->key), assVoid(iseg)) ;
    }

  aLook = look ;		/* static for cMapAddSeg() */

  iseg = arrayMax (look->segs) ;
  for (i = 0 ; i < keySetMax(kSet) ; ++i)
    { key = arr(kSet, i, KEY) ;
      if (!(obj = bsCreate (key)))
	continue ;

      if (class (key) == _VLocus)
	{ if (bsGetKey (obj, _Positive_clone, &clone) && 
	      (Clone = bsCreate (clone)))
	    { if (bsGetKey (Clone, _pMap, &contig) &&
		  bsGetData (Clone, _bsRight, _Int, &x1) &&
		  bsGetData (Clone, _bsRight, _Int, &x2) &&
		  assFind (look->ctgAss, assVoid(contig), &ictg))
		cMapAddSeg (key, assInt(ictg), x1, x2) ;
	      bsDestroy (Clone) ;
	    }
	}
      else if (class(key) == _VClone)
	{
	  if (bsGetKey (obj, _pMap, &contig) &&
	      bsGetData (obj, _bsRight, _Int, &x1) &&
	      bsGetData (obj, _bsRight, _Int, &x2) &&
	      assFind (look->ctgAss, assVoid(contig), &ictg))
	    { seg = cMapAddSeg (key, assInt(ictg), x1, x2) ;
	      if (bsGetData (obj, _In_Situ, _Int, &x1) &&
		  bsGetData (obj, _bsRight, _Int, &x2))
		{ seg->flag |= FLAG_INSITU ;
		  seg->p1 = x1/100.0 ; /* NB p is a float, x an int */
		  seg->p2 = x2/100.0 ;
		}
	    }
	  else if (gridClusterKey (key, map, 200))
	    for (k = 0 ; k < arrayMax(map) ; ++k)
	      { grid = arrp (map, k, GRIDMAP) ;
		if (!assFind (look->ctgAss, 
			      assVoid(grid->ctg), &ictg))
		  continue ;
		seg = cMapAddSeg (key, assInt(ictg), grid->x1, grid->x2) ;
		seg->flag = FLAG_HYB ;
	      }
	}
      else if (class(key) == _VSequence)
	{
	  units = arrayReCreate (units, 1024, BSunit) ;
	  if (bsFindTag (obj, _Homol) && bsFlatten (obj, 7, units))
	    for (k = 0 ; k < arrayMax (units) ; k += 7)
	      if (lexReClass (arr(units,k,BSunit).k, &clone, _VClone) &&
		  (Clone = bsCreate (clone)))
		{ if (bsGetKey (Clone, _pMap, &contig) &&
		      bsGetData (Clone, _bsRight, _Int, &x1) &&
		      bsGetData (Clone, _bsRight, _Int, &x2) &&
		      assFind (look->ctgAss, assVoid(contig), &ictg))
		    cMapAddSeg (key, assInt(ictg), 
				x1 + arr(units,k+5,BSunit).i * PMAP_FAC,
				x1 + arr(units,k+5,BSunit).i * PMAP_FAC) ;
		  bsDestroy (Clone) ;
		}
	}
      bsDestroy (obj) ;
    }
  
  arraySort (look->segs, cMapOrder) ;
  
  if (!isFullGenome)
    cMapDraw (look, 0) ;
}

static void cMapSubtract (void) 
{
  int i, j;
  SEG *seg, *seg2 ;
  KEYSET kSet = 0 ;
 
  LOOKGET("cMapSubtract") ;

  if (!keySetActive(&kSet, 0))
    { messout("First select a keySet window, thank you.") ;
      return ;
    }
  for (seg = arrp(look->segs,0,SEG), i = arrayMax(look->segs) ; i-- ; seg++)
    if (keySetFind (kSet, seg->key, &j))
      seg->flag |= FLAG_DELETING ;
  
  for (seg = seg2 = arrp(look->segs,0,SEG), i = arrayMax(look->segs) ; i-- ; seg++)
    if (!(seg->flag & FLAG_DELETING))
      *seg2++ = *seg ;
  arrayMax(look->segs) = seg2 - arrp(look->segs,0,SEG) + 1 ;

  cMapDraw (look, 0) ;
}

static void cMapHighlight (void)
{
  int i, j ;
  SEG *seg ;
  KEYSET kSet = 0 ;
 
  LOOKGET("cMapHighlight") ;

  if (!keySetActive(&kSet, 0))
    { messout("First select a keySet window, thank you.") ;
      return ;
    }
  seg = arrp(look->segs,0,SEG) ;
  i = arrayMax(look->segs) ;
  while(i--)
    { seg->flag &= ~FLAG_HIGHLIGHT ;
      if (keySetFind (kSet, seg->key, &j))
	seg->flag |= FLAG_HIGHLIGHT ;
      seg++ ;
    }
  cMapDraw (look, 0) ;
}

/************* special function for cDNA paper display **********/

static void cMapFullGenome (void)
{
  static char* chroms[] = {"I", "II", "III", "IV", "V", "X"} ;
  KEY key ;
  int i, xOff ;
  LOOKGET ("cdnaPaperDisplay") ;

  graphClear () ;
  graphColor (BLACK) ;
  isFullGenome = TRUE ;
  nx = 100 ;
  ny = 100 ;
  look->mag = 1.0 ;
  for (i = 0 ; i < 6 ; ++i)
    {
      lexword2key (chroms[i], &key, _VMap) ;
      look->key = key ;
      if (!cMapConvert (look, FALSE))
	messcrash ("lost a chromosome!") ;
      look->centre = 0.0 ;
      yCentre = 3.5 + 20 * (i / 3) ;
      xOff = -7 + 20 * (i % 3) ;
      xContig += xOff ; xItem += xOff ;
      graphTextHeight (1.5) ;
      graphText (chroms[i], xContig-7, MAP2GRAPH(look,-1)) ;
      graphTextHeight (0.0) ;

      { int ii ;
	COLUMN *col ;

	for (ii = 0 ; ii < arrayMax (look->columns) ; ++ii)
	  { col = arr(look->columns, ii, COLUMN*) ;
	    graphText (col->name, xItem + col->offset, MAP2GRAPH(look,-2)) ;
	  }
      }

      cMapDraw (look, 0) ;
      xContig -= xOff ; xItem -= xOff ;
    }

  isFullGenome = FALSE ;
}

static void cMapPartGenome (void)
{
  static char* chroms[] = {"II", "III", "X"} ;
  KEY key ;
  int i ;
  COLUMN *col ;
  LOOKGET ("cdnaPaperDisplay") ;

  graphClear () ;
  graphColor (BLACK) ;
  isFullGenome = TRUE ;
  nx = 100 ;
  ny = 100 ;
  look->mag = 0.8 ;
  yCentre = 3.5 ;
  xContig -= 6 ; xItem -= 6 ;

  for (i = 0 ; i < arrayMax (look->columns) ; ++i)
    { col = arr(look->columns, i, COLUMN*) ;
      graphText (col->name, xItem + col->offset, 2) ;
    }

  for (i = 0 ; i < 3 ; ++i)
    {
      lexword2key (chroms[i], &key, _VMap) ;
      look->key = key ;
      if (!cMapConvert (look, FALSE))
	messcrash ("lost a chromosome!") ;
      look->centre = 0.0 ;
      graphTextHeight (1.5) ;
      graphText (chroms[i], xContig-7, MAP2GRAPH(look,-1)) ;
      graphTextHeight (0.0) ;
      cMapDraw (look, 0) ;
      yCentre += look->mag * (look->max + 4) ;
    }

  xContig += 6 ; xItem += 6 ;
  isFullGenome = FALSE ;
}

/*******************************************************************/
/*******************************************************************/

static void cMapColDraw (LOOK look, COLUMN *col)
{ 
  int i, ihis ;
  float *px, y1, y2, w, x1, x2 ;
  KEY *pctg ;
  SEG *seg ;
  Array his, seq, binEnd ;
  int oldColour ;
  static Associator ctg2seg = 0 ;

  if (!col->bin || !col->wfac || !arrayMax (col->coords))
    return ;

  ctg2seg = assReCreate (ctg2seg) ;
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs, i, SEG) ;
      if (class(seg->key) == _VContig)
	assInsert (ctg2seg, assVoid(seg->key), seg) ;
    }

  oldColour = graphColor (col->colour) ;

  if (col->isLine)
    { 
      px = arrp (col->coords, 0, float) ;
      pctg = arrp (col->ctgs, 0, KEY) ;
      for (i = 0 ; i < arrayMax (col->coords) ; ++i, ++px, ++pctg)
	if (assFind (ctg2seg, assVoid(*pctg), &seg))
	  { x1 = seg->x1 + PMAP_FAC * (*px - seg->p1) ;
	    y1 = MAP2GRAPH (look, x1) ;
	    w = col->wfac * array(col->count, i, int) ;
	    graphLine (xItem + col->offset, y1, 
		       xItem + col->offset + w, y1) ;
	  }
    }
  else
    {
      his = arrayCreate (32, int) ;
      seq = arrayCreate (32, int) ;
      binEnd = arrayCreate (32, float) ;

      px = arrp (col->coords, 0, float) ;
      pctg = arrp (col->ctgs, 0, KEY) ;
      ihis = -1 ;
      for (i = 0 ; i < arrayMax (col->coords) ; ++i, ++px, ++pctg)
	if (assFind (ctg2seg, assVoid(*pctg), &seg))
	  { w = seg->x1 + PMAP_FAC * (*px - seg->p1) ;
	    if (arrayMax (col->seqlen))
	      {
#ifdef VARIABLE
		if (ihis == -1)	/* start */
		  array(binEnd,++ihis,float) = w ;
		else if (arr(seq, ihis, int) > 1000000 * col->bin)
		  ++ihis ;
		array(binEnd,ihis+1,float) = w ;
#else
		ihis = (int)(w/col->bin) ;
#endif
		array(his, ihis, int) += arr(col->count, i, int) ;
		array(seq, ihis, int) += arr(col->seqlen, i, int) ;
	      }
	    else
	      ++array(his, (int)(w/col->bin), int) ;
	  }
      
      for (i = 0 ; i < arrayMax (his) ; ++i)
	if ((w = col->wfac * arr(his, i, int)))
	  { x1 = i*col->bin ; 
	    x2 = (i+1)*col->bin ;
	    if (arrayMax (col->seqlen))
	      w *= 1000000.0 / arr(seq, i, int) ;
	    else
	      { w /= col->bin ;
		if (x2 > look->max)
		  w *= (x2 - x1) / (look->max - x1) ;
	      }
	    if (x2 > look->max)
	      x2 = look->max ;
	    y1 = MAP2GRAPH (look, x1) ;
	    y2 = MAP2GRAPH (look, x2) ;
	    graphFillRectangle (xItem + col->offset, y1, 
				xItem + col->offset + w, y2) ;
	  }

      arrayDestroy (his) ;
      arrayDestroy (seq) ;
      arrayDestroy (binEnd) ;
    }
  graphColor (oldColour) ;

  assDestroy (ctg2seg) ;
}

/********** little window for controlling the column *******/

static int colGraphHolder ;

static Graph cMapColGraphCreate (LOOK look, COLUMN *col)
{ 
  int i, n = 0 ;

  if (graphActivate (col->graph))
    { graphPop () ;
      return col->graph ;
    }
  col->graph = graphCreate (TEXT_FIT, "column", 0, 0, .3, .2) ;

  for (i = 0 ; i < arrayMax(col->count) ; ++i)
    n += arr(col->count, i, int) ;
  graphText (messprintf ("%d objects", n), 2, 0.5) ;

  graphFloatEditor ("Histogram bin", &col->bin, 1, 2, 0) ;
  graphFloatEditor ("Width factor", &col->wfac, 1, 3, 0) ;
  graphIntEditor ("Colour", &col->colour, 1, 4, 0) ;
  graphIntEditor ("Offset", &col->offset, 1, 5, 0) ;
  graphWordEditor ("Name", col->name, 30, 1, 6, 0) ;
  graphIntEditor ("Is line", &col->isLine, 1, 7, 0) ;

  graphRedraw () ;

  graphAssociate (&colGraphHolder, col) ;

  return col->graph ;
}

/*******************************************/

static void cMapColEdit (void)
{ 
  int i ;
  LOOKGET("cMapColEdit") ;

  for (i = 0 ; i < arrayMax (look->columns) ; ++i)
    cMapColGraphCreate (look, arr(look->columns, i, COLUMN*)) ;
}

/*******************************************/

static void cMapColFillFromKeySet (LOOK look, COLUMN *col)
{ 
  KEY key, contig ;
  KEYSET kSet = 0 ;
  GRIDMAP *grid ;
  static Array map = 0 ;
  int i, k, x1, x2, n = 0 ;
  OBJ obj ;

  if (!keySetActive (&kSet, 0))
    return ;

  col->coords = arrayReCreate (col->coords, 128, float) ;

  map = arrayReCreate (map, 8, GRIDMAP) ;

  for (i = 0 ; i < keySetMax(kSet) ; ++i)
    { key = arr(kSet, i, KEY) ;
      if (!(obj = bsCreate (key)))
	continue ;

      if (class(key) == _VClone)
	{
	  if (bsGetKey (obj, _pMap, &contig) &&
	      bsGetData (obj, _bsRight, _Int, &x1) &&
	      bsGetData (obj, _bsRight, _Int, &x2))
	    { array (col->coords, n, float) = (x1+x2)/2 ;
	      array (col->ctgs, n, KEY) = contig ;
	      array (col->count, n, int) = 1 ;
	      ++n ;
	    }
	  else if (gridClusterKey (key, map, 200))
	    for (k = 0 ; k < arrayMax(map) ; ++k)
	      { grid = arrp (map, k, GRIDMAP) ;
		array (col->coords, n, float) = (grid->x1+grid->x2)/2 ;
		array (col->ctgs, n, KEY) = grid->ctg ;
		array (col->count, n, int) = 1 ;
		++n ;
	      }
	}
    }
}

static void cMapColFillFromFile (LOOK look, COLUMN *col, FILE *fil)
{ 
  int level, n = 0, x1, x2, k ;
  OBJ obj ;
  KEY key, contig ;
  char *word ;

/* expected file format: CLONE [N] [SEQLEN] */

  level = freesetfile (fil, 0) ;
  while (freecard (level))
    { if (!(word = freeword ()) ||
	  !lexword2key (word, &key, _VClone) ||
	  !(obj = bsCreate (key)))
	continue ;
      if (bsGetKey (obj, _pMap, &contig) &&
	  bsGetData (obj, _bsRight, _Int, &x1) &&
	  bsGetData (obj, _bsRight, _Int, &x2))
	{ array (col->coords, n, float) = (x1+x2)/2 ;
	  array (col->ctgs, n, KEY) = contig ;
	  k = 1 ; freeint (&k) ;
	  array (col->count, n, int) = k ;
	  if (freeint (&k))
	    array (col->seqlen, n, int) = k ;
	  ++n ;
	}
      bsDestroy (obj) ;
    }

  filclose (fil) ;
}

/*******************************************/

static void cMapColFinalise (void *p)
{ COLUMN *col = (COLUMN*)p ;
  if (graphActivate (col->graph)) 
    graphDestroy() ;
}

static COLUMN *cMapColCreate (LOOK look)
{ 
  COLUMN *col ;

  col = array (look->columns, arrayMax(look->columns), COLUMN*) =
    (COLUMN*) halloc (sizeof (COLUMN), look->handle) ;
  col->bin = 0 ;
  col->wfac = 0 ;
  col->colour = BLACK ;
  col->coords = arrayHandleCreate (128, float, look->handle) ;
  col->ctgs = arrayHandleCreate (128, KEY, look->handle) ;
  col->seqlen = arrayHandleCreate (128, int, look->handle) ;
  col->count = arrayHandleCreate (128, int, look->handle) ;

  blockSetFinalise (col, cMapColFinalise) ;
  return col ;
}

static void cMapColFromKeySet (void)
{ 
  COLUMN *col ;
  LOOKGET("cMapColCreate") ;
  
  col = cMapColCreate (look) ;
  cMapColFillFromKeySet (look, col) ;
  cMapColGraphCreate (look, col) ;
}

static void cMapColFromFile (void)
{ 
  COLUMN *col ;
  FILE *fil ;
  LOOKGET("cMapColCreate") ;

  if (!(fil = filqueryopen (0, 0, "hisinf", "r", "info for column")))
    return ;

  col = cMapColCreate (look) ;
  cMapColFillFromFile (look, col, fil) ;
  cMapColGraphCreate (look, col) ;
}

static void colMake (LOOK look, char *fname, char *name,
		     float bin, float wfac, int colour, int offset)
{ 
  FILE *fil ;
  COLUMN *col ;

  if (!(fil = filopen (fname, "hisinf", "r")))
    return ;

  col = cMapColCreate (look) ;
  col->bin = bin ; col->wfac = wfac ; 
  col->colour = colour ; col->offset = offset ;
  col->isLine = TRUE ;
  strcpy (col->name, name) ;
  cMapColFillFromFile (look, col, fil) ;
}

static void cMapHardColumns (void)
{ 
  LOOKGET("cMapColCreate") ;

  colMake (look, "/nfs/disk37/rd/worm.data/maps.971209/1", "L",
	   0.1, 1, BLACK, 2) ;
  colMake (look, "/nfs/disk37/rd/worm.data/maps.971209/3", "S",
	   0.1, 1, DARKGREEN, 3) ;
  colMake (look, "/nfs/disk37/rd/worm.data/maps.971209/4", "A",
	   0.1, 1, BLUE, 5) ;
  colMake (look, "/nfs/disk37/rd/worm.data/maps.971209/5", "P",
	   0.1, 1, MAGENTA, 6) ;

#ifdef REPS
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/Q", "REP-Q",
	   1.5, 0.4, DARKGRAY, 5) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/D", "REP-D",
	   1.5, 0.3, DARKGRAY, 20) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/N", "REP-N",
	   1.5, 0.08, DARKGRAY, 35) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/T", "REP-T",
	   1.5, 0.2, DARKRED, 5) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/M", "REP-M",
	   1.5, 0.2, DARKGRAY, 5) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/A", "REP-A",
	   1.5, 0.08, DARKRED, 20) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/A", "REP-A",
	   1.0, 0.02, DARKRED, 0) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/A", "REP-D",
	   1.0, 0.05, DARKRED, 5) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/E", "REP-E",
	   1.0, 0.02, DARKRED, 10) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/F", "REP-F",
	   1.0, 0.05, DARKRED, 15) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/G", "REP-G",
	   1.0, 0.05, DARKRED, 20) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/J", "REP-J",
	   1.0, 0.08, DARKRED, 25) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/L", "REP-L",
	   1.0, 0.10, DARKRED, 30) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/M", "REP-M",
	   1.0, 0.05, DARKRED, 35) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/N", "REP-N",
	   1.0, 0.02, DARKRED, 40) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/P", "REP-P",
	   1.0, 0.15, DARKRED, 45) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/Q", "REP-Q",
	   1.0, 0.15, DARKRED, 50) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/R", "REP-R",
	   1.0, 0.20, DARKRED, 55) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/S", "REP-S",
	   1.0, 0.08, DARKRED, 60) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/T", "REP-T",
	   1.0, 0.07, DARKRED, 65) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/V", "REP-V",
	   1.0, 0.03, DARKRED, 70) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/W", "REP-W",
	   1.0, 0.03, DARKRED, 75) ;
  colMake (look, "/nfs/disk37/rd/worm.data/9510/REP/hisinf/X", "REP-X",
	   1.0, 0.10, DARKRED, 80) ;
#endif

  cMapPartGenome () ;
}

/*******************************************************************/
/*******************************************************************/
 
 
