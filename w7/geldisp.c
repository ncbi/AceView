/*  File: geldisp.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **   Display gels and corellates them to the sequence
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 16 12:10 1998 (edgrif)
 * * Jul 16 12:05 1998 (edgrif): Inserted the new, public fmap.h,
 *     AND - removed vile redeclarations of fmap functions
 *         - changed LOOK to LOOKGEL
 *         - changed topMargin to topMargin_G
 * * Apr 29 12:14 1996 (rd)
 * * Dec 16 13:59 1991 (mieg): follow and color the dna
 * Created: Wed Dec 11 18:03:54 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: geldisp.c,v 1.4 2014/11/30 03:20:46 mieg Exp $ */


#include "acedb.h"
#include "fmap.h"
#include "keyset.h"
#include "display.h"
#include "key.h" 
#include "lex.h"
#include "dna.h"
#include "peptide.h"
#include "restriction.h"
#include "query.h"
#include "bs.h"
#include "systags.h"
#include "classes.h"
#include "tags.h"


typedef struct LOOKSTUFF
  { int   magic;        /* == MAGIC */
    Graph graph ;
    int   activeBox, positionBox, currentEntryBox, enzymeBox ;
    char positionBuffer [10] ;
    Array enzymes ;
    BOOL coor ;         /* show coordinates toggle */
    float centre ;
    float mag ;		/* gel units * mag = window coords */
    Array segs,      	/* array of SEG's from the obj */
          boxIndex ;    /* if >256 a SEG, else something special */
    struct
      { float fac,offset ;
      }   cursor ;
  } *LOOKGEL ;



BOOL gelDisplay (void) ;
static void gelDestroy (void) ;
static void gelPick (int box, double x , double y) ;
static void gelMiddleDown (double x, double y) ;
static void gelKbd (int k) ;
static void gelDraw (LOOKGEL look, KEY key) ;
static void gelResize (void) ;
static void gelComputeGel(LOOKGEL look) ;
static void gelSelect (LOOKGEL look, int box) ;
static void gelFollow (LOOKGEL look, double x , double y) ;
static void wholeChromosome (void) ;

static MENUOPT gelMenu[] =
            { {graphDestroy, "Quit"},
	      {help,"Help"},
	       {graphPrint,"Print"},
		{0, 0}
            } ;


/* This is grim...these variables are all in the global name space ready to  */
/* clash with just about anything...sigh....                                 */
/* If there must be globals they should be contained in a struct with a      */
/* good name...that way the globals are referred to as structname.global...  */
/* For now I have changed topMargin as it clashes with map.h                 */
static int MAGIC = 284657 ;
static int minLiveBox , maxLiveBox, scaleBox = 0 ;
static int nx, ny , topMargin_G ;
static float ny2 ;


#define GEL2GRAPH(look,x) (look->mag * ((x) - look->centre) + ny2 + topMargin_G)

typedef struct
  { Array itvs ;
    char *text ;
    KEY seqKey ;
    int from, length ;
    BOOL observed ;  /* not computed from dna, but read in ace data */
  } ENZ ;

typedef struct
{ int len , from ; float log ;  } INTERVAL ;

typedef struct
  { ENZ* enzp ;
    INTERVAL itv ;
    int color ;
    int flag ;
  } SEG ;

static Graph oldGelGraph = 0 ;

/* These colours come from graph.h and fmap.h                                */
#define MAXGELCOLOR 4
static int GELTINT[] = { TINT_CYAN, TINT_LIGHTGREEN, TINT_MAGENTA, TINT_YELLOW } ;
static int GELCOLOR[] = { BLUE, GREEN, MAGENTA, YELLOW } ;
#define FRIEND_COLOR CYAN


#define LOOKGET(name)     LOOKGEL look ; \
                          if (!graphAssFind (&MAGIC, &look)) \
		            messcrash ("graph not found in %s",name) ; \
                          if (look->magic != MAGIC) \
                            messcrash ("%s received a wrong pointer",name)

/*****************************************/
/*********** action routines ***********/

static void gelEdit (char *junk)
{
  int box , i ;
  
  LOOKGET("gelEdit") ;
  
  box = look->currentEntryBox ;
  i  = (box - look->enzymeBox - 1)/2 ;

    /* One empty entry for a new instruction */ 
  if(i == arrayMax(look->enzymes) - 1)
    {
      ENZ *enzp = arrayp(look->enzymes,i+1,ENZ) ; 
      enzp->text = messalloc(80) ;
      enzp->itvs = arrayCreate(80, INTERVAL) ;
      gelDraw(look, 0) ;
      /* Reactivate box  */
      box =  2 * i + look->enzymeBox + 1 ;
      gelPick(box, (float)(-1), 0) ; 
    }
  /* execute */
  gelComputeGel(look) ;
}

/**********************************/
/**********************************/

BOOL gelDisplay (void)
{
  LOOKGEL oldlook, look=(LOOKGEL)messalloc(sizeof(struct LOOKSTUFF)) ;
  KEY  curr = 0 ;

  look->magic = MAGIC;

  look->centre = 0 ;   	/* centre on parent object if poss */
 
  look->boxIndex = arrayCreate (64,SEG) ;
  look->activeBox = 0 ;

    /* Try to keep same magnification */
  if (graphActivate(oldGelGraph) &&
      graphAssFind (&MAGIC, &oldlook) &&
      oldlook &&
      oldlook->magic == MAGIC )
    look->mag = oldlook->mag ;
  else
    look->mag = 10 ;
  
  if (oldGelGraph)
    { 
      gelDestroy () ;
      graphAssRemove (&MAGIC) ;
    }
  else 
    { if (!displayCreate(DtGel))
	goto abort ;
      
      graphRegister (RESIZE,gelResize) ;
      graphRegister (DESTROY, gelDestroy) ;
      graphRegister (PICK, gelPick) ;
      graphRegister (MIDDLE_DOWN, gelMiddleDown) ;
      graphRegister (KEYBOARD, gelKbd) ;
      graphMenu (gelMenu) ;
    }

  graphAssociate (&MAGIC, look) ;
  look->graph = oldGelGraph = graphActive() ;
  look->enzymes =  arrayCreate(5, ENZ) ;
  look->coor = TRUE ;
    
  { ENZ * enzp = arrayp(look->enzymes,0,ENZ) ; 
    enzp->text = messalloc(80) ;
    enzp->itvs = arrayCreate(80, INTERVAL) ;
  }
 
  gelDraw (look, curr) ;

  return TRUE ;

abort :
  messfree (look) ;
  return FALSE ;
}

/************************************************************/
/***************** Registered rourtines *********************/

static void gelDestroy (void)
{ int i ;
  LOOKGET("gelDestroy") ;

  if (arrayExists(look->enzymes))
    { i  = arrayMax(look->enzymes) ;
      while(i--)
	{ messfree(array(look->enzymes,i,ENZ).text) ;
	  arrayDestroy(array(look->enzymes,i,ENZ).itvs) ;
	}
    }
      
  arrayDestroy(look->enzymes) ;
  arrayDestroy(look->segs) ;
  arrayDestroy(look->boxIndex) ;

  oldGelGraph = 0 ;
  look->magic = 0 ;
  look->graph = 0 ;
  messfree (look) ;
}

/***************************************************************/
/******** coordinated dragging of pos and scale boxes **********/
/********  copied from gmapdisp.c ******************************/

#define xScale 6
static int posBox = 0 ;
static LOOKGEL lookDrag ;

static void scaleBoxDrag (float *x, float *y, BOOL isDone)
{
  static BOOL isDragging = FALSE ;
  static float newCentre, yCentre ;
  
  if (isDragging)
    graphXorBox (posBox,xScale-0.5,yCentre) ;
  else
    isDragging = TRUE ;
    
   *x = 0 ;
 
  newCentre = lookDrag->cursor.offset + 
    		(*y - 3) / lookDrag->cursor.fac ;
  yCentre = GEL2GRAPH(lookDrag,newCentre) ;

  if (isDone)
    { lookDrag->centre = newCentre ;
      gelDraw (lookDrag, 0) ;
      isDragging = FALSE ;
    }
  else
    graphXorBox (posBox,xScale-0.5,yCentre) ;
}

static void posBoxDrag (float *x, float *y, BOOL isDone)
{
  static BOOL isDragging = FALSE ;
  static float newCentre,yScaleBox ;

  if (isDragging)
    graphXorBox (scaleBox,0,yScaleBox) ;
  else
    isDragging = TRUE ;

  *x = xScale-0.5 ;

  newCentre = lookDrag->centre + (*y - ny2 - topMargin_G ) / lookDrag->mag ;
  yScaleBox = 5 + 
    (newCentre - lookDrag->cursor.offset) * lookDrag->cursor.fac ;
  if (isDone)
    { lookDrag->centre = newCentre ;
      gelDraw (lookDrag, 0) ;
      isDragging = FALSE ;
    }
  else
    graphXorBox (scaleBox,0,yScaleBox) ;
}

static void makeDragBoxes (void)
{
  float x1,x2,y1,y2 ;
  int ny21 = ny2 + topMargin_G ;

  if (posBox) graphBoxDraw (posBox, WHITE, WHITE) ;
  posBox = graphBoxStart() ;
  graphLine (xScale-0.5,ny21 ,xScale-0.5,ny21+1) ;
  graphLine (200,ny21,200,ny21+1) ;	/* to make box wide */
  graphBoxEnd () ;
  graphBoxDraw (posBox, BLACK, TRANSPARENT) ;

  if (scaleBox) 
    { graphBoxDraw (scaleBox, WHITE, WHITE) ;
    graphBoxDim (scaleBox, &x1, &y1, &x2, &y2) ;
    scaleBox = graphBoxStart() ;
    graphFillRectangle (5.25, y1, 5.75, y2) ;
    graphLine (0, y1, 0, y2) ;	/* to make box wide */
    graphBoxEnd () ;
    graphBoxDraw (scaleBox, BLACK, TRANSPARENT) ;
    }
}

/*****************************************/
/*****************************************/

static int intervalOrder(const void *a, const void *b)
{
  if ( ((const INTERVAL*)a)->log  > ((const INTERVAL*)b)->log )
    return 1 ;
  else
    if ( ((const INTERVAL*)a)->log  > ((const INTERVAL*)b)->log )
      return - 1 ;
    else
      return  0 ;
}

/********/
 
static void gelGetLengths(LOOKGEL look, Array dna , Array sites, int lane, int from)
{
  Array itvs = array(look->enzymes,lane,ENZ).itvs ;
  int i, j = 0, m, n ;
  float y ;

  arrayDestroy(itvs) ;
  itvs =  array(look->enzymes,lane,ENZ).itvs =
    arrayCreate(20,INTERVAL) ;

  arraySort(sites,siteOrder) ;
  n = arrayMax(sites) ? arr (sites, 0, Site).i : 0 ;
  for (i = 1; i < arrayMax(sites); i++)
    { 
      m = array(sites, i, Site).i ;
     	  
      array(itvs, j, INTERVAL).from = n - from ;
      array(itvs, j, INTERVAL).len = m - n ;
      if (m - n < 10)
	n = m - 10 ;  /* Prevent log errors. */
	     /* 
	        Log relation between lenght and distance of migration:
		figure 6.1 page 6.5 ; Molecular cloning, a laboratory manual
		second ed. Sambrook, Fritsch, Maniatis
		ColdSping Harbor lab. press, 1989
	
                log10 isdeclared in mystdlib.h 
		the limits 10 to 10000 (so 0 <= y < 11.5)
		are recalled in a graphTextMessage.
		The inverse relation is used in drawScale.
	     */
      y = 3.5 * ( 4 - log10((double)(m - n)) ) ; 
      array(itvs, j++,INTERVAL).log = (y > 0 ? y : 0 ) ;
      n = m ;
    }
  arraySort(itvs,intervalOrder) ;
}

static void gelComputeGel(LOOKGEL look)
{ int   lane, from, to, length , origin ;
  KEY seqKey ;
  Array dna,  colors, sites ;
  void *fMapLook ;
  char *tt ;
  Stack s , siteNames = 0 ;
 
  if(fMapActive(&dna,&colors,&seqKey,&fMapLook))
    dnaRepaint(colors) ;
  else
    { messout("First select a dna window, thank you.") ;
      return ;
    }
 
  fMapFindZone (fMapLook, &from, &to, &origin) ;
  length = to - from + 1 ;
  lane = (look->currentEntryBox  - look->enzymeBox - 1)/2 ;
  if(lane >= arrayMax(look->enzymes))
    messcrash("Lane error in gelcomputegel") ;
  s = dnacptAnalyseRestrictionBox(array(look->enzymes,lane,ENZ).text, FALSE, &siteNames) ;
  if(!s)
    { gelDraw(look,0) ;
      return ;
    }
  array(look->enzymes,lane,ENZ).seqKey = seqKey ; 
  array(look->enzymes,lane,ENZ).from = from ;
  array(look->enzymes,lane,ENZ).length = length ;
  
  sites = arrayCreate(20,Site) ;
  tt = pepGetTranslationTable (seqKey, 0) ;
  dnacptMultipleMatch(sites, dna, 0, 0, colors,s, 
		      siteNames,from, length, FALSE,0, 0, tt) ;
  gelGetLengths(look, dna, sites, lane, from) ;
  stackDestroy(s) ;
 
  fMapRegisterSites (fMapLook, sites, siteNames) ; /* sites and sitenames will be destroyed there */
  fMapReDrawDNA(fMapLook) ;
  graphActivate(look->graph) ;
  wholeChromosome() ;
}

/***********************************/

static void gelPick (int box,  double x , double y) 
{
  int i ;
  LOOKGET("gelPick") ;


 		/* a control box */
  if(!box)
    return ;
  else if(box == scaleBox)
   { lookDrag = look ;
     makeDragBoxes () ;
     graphBoxDrag (scaleBox,scaleBoxDrag) ;
     return ;
   }
  else if(box == posBox)
    { lookDrag = look ;
      makeDragBoxes () ;
      graphBoxDrag (posBox,posBoxDrag) ;
      return ;
    }
  
  if (box > look->enzymeBox ) /* a textEntry */
    { 
      i = (box - look->enzymeBox - 1)/ 2 ;
      if (i < arrayMax(look->enzymes))
	{ look->activeBox = box ;
	  look->currentEntryBox = 
	    graphTextEntry 
	      ( array(look->enzymes, i,ENZ).text,
	       0,(float) x,0,0) ;
	}
      return ;
     }

   if (box == look->activeBox)         /* a second hit - follow it */
    gelFollow (look,x,y) ;
  else
    gelSelect (look,box) ;
}

/********* use middle button for cursor **********/

static double oldy , oldDy , oldx;
static BOOL dragFast ;
#define DRAGFASTLIMIT 6

static void gelMiddleDrag (double x, double y) 
{
  if(dragFast)
    { graphXorLine (0, oldy - oldDy, DRAGFASTLIMIT, oldy - oldDy) ;
      graphXorLine (0, oldy + oldDy, DRAGFASTLIMIT, oldy + oldDy) ;
    }
  else
    graphXorLine (DRAGFASTLIMIT, oldy, nx, oldy) ;

  oldy = y ;

  if(dragFast)
    { oldDy *= exp ((x - oldx) / 25.) ;
      oldx = x ;
      graphXorLine (0, y - oldDy, DRAGFASTLIMIT, y - oldDy) ;
      graphXorLine (0, y + oldDy, DRAGFASTLIMIT, y + oldDy) ;
    }
  else
    graphXorLine (DRAGFASTLIMIT, y, nx, y) ;
}

static void gelMiddleUp (double x, double y) 
{  float x1,x2,y1,y2 ;
  LOOKGET("gelMiddleUp") ;

  if(dragFast)
    { graphBoxDim (scaleBox, &x1, &y1, &x2, &y2) ;
      look->mag *= (y2 - y1) / (2. * oldDy) ;
      look->centre = look->cursor.offset + 
	(y - (y2 - y1) * .5 - 3.05 - topMargin_G ) / look->cursor.fac ;
    }
  else
   look->centre = look->centre + (y - .5  - ny2 - topMargin_G ) / look->mag ;
  gelDraw (look, 0) ;
}

static void gelMiddleDown (double x, double y) 
{  float x1,x2,y1,y2 ;
  LOOKGET("gelMiddleDown") ;

  if (!scaleBox) return ;
  graphBoxDim (scaleBox, &x1, &y1, &x2, &y2) ;
  oldDy = (y2 - y1) / 2. ;

  lookDrag = look ;
  dragFast = x < DRAGFASTLIMIT ? TRUE : FALSE ;

  if(dragFast)
    { graphXorLine (0, y - oldDy, DRAGFASTLIMIT, y - oldDy) ;
      graphXorLine (0, y + oldDy, DRAGFASTLIMIT, y + oldDy) ;
    }
  else
    graphXorLine (DRAGFASTLIMIT, y, nx, y) ;
   
  oldx = x ;
  oldy = y ;
  graphRegister (MIDDLE_DRAG, gelMiddleDrag) ;	/* must redo */
  graphRegister (MIDDLE_UP, gelMiddleUp) ;
}

/*****************/

static void gelKbd (int k)
{
  int  box ;
  LOOKGET("gelKbd") ;

  if (!look->activeBox)
    return ;

  box = look->activeBox ;
  switch (k)
    {
    case UP_KEY :
      if (box > minLiveBox)
	--box ;
      break ;
    case DOWN_KEY :
      if (box < maxLiveBox - 1)
	++box ;
      break ;
    }
  if (look->activeBox != box)
    { gelSelect (look,box) ;
      gelFollow (look,0.,0.) ;
    }
  else
    gelSelect (look,box) ;
}

/**********************************/

static void gelResize (void)
{
  LOOKGET("gelResize") ;

  gelDraw (look, 0) ;
}

/**********************************/

static void gelToggleCoordinates (void)
{
  LOOKGET("gelResize") ;

  look->coor = ~look->coor ;
  gelDraw (look, 0) ;
}

/*****************************************************************/
/****************** Drawing routines *****************************/

static void drawScale (LOOKGEL look)
{
  float cutoff = 5 / look->mag ;
  float unit = 0.01 ;
  float subunit = 0.001 ;
 
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

  graphBoxStart() ;
/* pas utile en log

  x = look->centre - (ny2-2)/look->mag ;
  if(x<0)
    x = 0 ;
  x = unit * (((x>=0)?1:0) + (int)(x/unit)) ;
  while ((y = GEL2GRAPH(look, x)) < ny  - 1)
    { graphLine (xScale-1.5,y,xScale-0.5,y) ;
      x1 = pow(10.0, 4.0 -  ((double) x)/3.5) ; 
      graphText (messprintf ("%-5.0f",x1),xScale,y-.5) ;
      x += unit ;
    }

  x = look->centre - (ny2-2)/look->mag ;
  if(x<0)
    x = 0 ;
  x = subunit * (((x>=0)?1:0) + (int)(x/subunit)) ;
  while ((y = GEL2GRAPH(look, x)) < ny  - 1)
    { graphLine (xScale-1.0,y,xScale-0.5,y) ;
      x += subunit ;
    }

 //  graphLine (xScale-0.5,2 + topMargin_G ,xScale-0.5,ny - 0.5 ) ;
  graphBoxEnd() ;  
  array(look->boxIndex,posBox=graphBoxStart(),SEG*) = (SEG*)5 ;
  graphFillRectangle (xScale-1,ny2 + topMargin_G ,xScale,ny2+1+topMargin_G) ;
*/
  graphBoxEnd () ;
  if (posBox) graphBoxDraw (posBox,DARKGREEN,TRANSPARENT) ;
}

static void drawChromosome (LOOKGEL look)
{
  float y,pageful ;
  int i , x0, x1 ;
  Array itvs ;

  x1 = 0 ;
  for (i=0; i<arrayMax(look->enzymes); i++)
    { itvs = array(look->enzymes,i,ENZ).itvs ;
      x0 =  arrayMax(itvs) ?
	arrp(itvs, arrayMax(itvs) - 1, INTERVAL)->len  : 0 ;
      if(x0 > x1)
	x1 = x0 ;
    }
  x0 = 0 ; if(!x1) x1 = 1 ;
  look->cursor.fac = (ny - 6.0 -topMargin_G) / (x1 - x0) ;

  graphBoxStart() ;
         /* Draw complete black line */
  graphFillRectangle (2.25,3+topMargin_G,2.75,ny-3) ;
  graphBoxEnd() ;

        /* Draw Green scaleBox */
  pageful = (ny - topMargin_G - 3) / look->mag ;	
  look->cursor.offset = x0 + 0.5*pageful ;
  y = 3 + topMargin_G + (look->centre - look->cursor.offset) * look->cursor.fac ;

 
  scaleBox=graphBoxStart() ;
  graphRectangle (2, y, 3, y + pageful*look->cursor.fac) ;
  graphBoxEnd () ;
  graphBoxDraw (scaleBox,DARKGREEN,GREEN) ;

       /* Draw oblique line outside of boxes */ 
  graphColor (DARKGRAY) ;
  graphLine (2.5, y, xScale - 1, 2.5 + topMargin_G) ;
  graphLine (2.5, y + pageful*look->cursor.fac, xScale-1, ny-1) ;

  graphColor (BLACK) ;
}

/*************************************************************/
/****** magnification control ******/
  /* JTM 18-6 Recenters the whole chromo correctly */
static void wholeChromosome (void)
{ float  x, xmax = 1 ;
  int   i, j ;
  ENZ *enzp ;
  LOOKGET("wholeChromosome") ; 

  for (i = 0 ; i < arrayMax(look->enzymes) ; ++i)
    { enzp = arrayp(look->enzymes,i, ENZ) ;
      if ((j = arrayMax(enzp->itvs)))
	if (( x = arrp(enzp->itvs,j-1,INTERVAL)->log) > xmax)
	  xmax = x ;
    }
  look->mag = (2*ny2-5)/xmax ;
  look->centre = xmax / 2  ; 
  gelDraw (look, 0) ;
}

static void zoomIn (void)
{ LOOKGET("zoomIn") ; look->mag *= 2 ; gelDraw (look, 0) ;}

static void zoomOut (void)
{ LOOKGET("zoomOut") ; look->mag /= 2 ; gelDraw (look, 0) ;}

static MENUOPT buttonOpts[] = {
  {wholeChromosome, "Whole"}, 
  {zoomIn, "Zoom In"},
  {zoomOut, "Zoom Out"},
  {gelToggleCoordinates, "Coordinates"},
   {0, 0}} ;

/**************************************/

static void gelDraw (LOOKGEL look, KEY dummy)
{ int   i, icol, j, ibox ;
  float x, width, y ;
  ENZ *enzp ;
  SEG   *seg ;
  INTERVAL *itvp ;
 
  if(graphActivate(look->graph))
    graphPop() ;
  else
    messcrash("gelDraw lost its graph") ;
  graphClear () ;
  posBox = scaleBox = 0 ;

  graphFitBounds (&nx,&ny) ;
  
  topMargin_G = 5 + 1.4 * arrayMax(look->enzymes) ; 
   if (ny < 10 || nx < xScale - 15)
    { messout ("Sorry, this window is too small for a  gel") ;
      return ;
    }

  ny2 = 0.5 * (ny - topMargin_G)  ;

   look->activeBox = 0 ;
  arrayMax (look->boxIndex) = 0 ;

  drawScale (look) ;
  drawChromosome (look) ;
  graphText("Position:", 4.0, 4.0) ;
  look->positionBox = graphBoxStart() ;
  graphTextPtr(look->positionBuffer,14,  4, 10) ;
  graphBoxEnd() ;
  graphText("Resolution:", 4.0, 5.0) ;
  graphText("10 to 10000", 8.0, 6.0) ;
  graphText("base pairs", 8.0, 7.0) ;

              /* must come at end of preamble boxes*/
  minLiveBox = graphBoxStart() + 1 ;
  graphBoxEnd () ;

  graphColor (DARKGRAY) ;
  graphFillRectangle (xScale-2.6,2,xScale-2.3,ny-0.5) ;
  graphColor (BLACK) ;

  width = (nx - xScale ) / (1 + arrayMax(look->enzymes));
  if (width < 2 )
    { messout ("Sorry, this window is too small for a so many lanes, make it bigger.") ;
      return ;
    }
 
  y = GEL2GRAPH(look, 0) ;
  graphLine(xScale, y, nx,y) ;
  for (i = 0, icol = 0 ; i < arrayMax(look->enzymes) ; ++i)
    { enzp = arrp(look->enzymes,i,ENZ) ;
      if (enzp->observed && icol > 0) icol-- ;
      x = xScale + 5  + i*width ;
      for(j=0; j < arrayMax(enzp->itvs); j++)
	{
	  y = arrp(enzp->itvs, j, INTERVAL)->log ;
	  y = GEL2GRAPH(look,y) ;

	  if (y >= topMargin_G   && y < ny-1)
	    { 
	      ibox = graphBoxStart() ;
	      seg = arrayp(look->boxIndex,ibox,SEG)  ;
	      seg->enzp = enzp ;
	      seg->color = icol % MAXGELCOLOR ;
	      itvp = arrp(enzp->itvs,j,INTERVAL) ;
	      (seg->itv).log = itvp->log ;
	      (seg->itv).len = itvp->len ;
	      (seg->itv).from = itvp->from ;
	      graphRectangle(x, y, x + width, y + .2) ;
	      graphBoxEnd() ;
	      graphBoxDraw(ibox, BLACK, GELCOLOR[seg->color]) ;
	      if (look->coor)
		graphText(messprintf("%d", itvp->len), x , y - .8) ;
	    }
	}
      if (!enzp->observed) icol++ ;
    }
  
  maxLiveBox = arrayMax(look->boxIndex) ;
  graphText ("Enzymes :", 16.0, 2.5) ;
  graphText ("Example : aatgc ; hindIII; ttkswatg ;",
	     26.0, 2.5) ;
                /* Draw buttons last to be on top of graph */
  /*but before the textEntryBox since i may add one at the end */
  graphButtons (buttonOpts, 5, 0.5, nx) ; 
 
 look->enzymeBox = graphBoxStart() ;
  if(look->enzymes)
    { for(i=0; i<arrayMax(look->enzymes); i++)
	look->currentEntryBox = 
		graphTextEntry 
		  ( array(look->enzymes, i,ENZ).text,
		   36, 26, i*1.4 + 3.9, gelEdit) ;
    }
  graphBoxEnd() ;
  look->activeBox = 0 ;

  graphRedraw () ;
}

/*********************************************************/

static void gelSelect (LOOKGEL look, int box)
{
  SEG *seg ;

  if (box >= maxLiveBox || box < minLiveBox)
    return ;

  if (look->activeBox)
    { seg = arrp(look->boxIndex,look->activeBox,SEG) ;
      if (seg)
	graphBoxDraw (look->activeBox, BLACK, GELCOLOR[seg->color]) ;
    }

  look->activeBox = box ;

  seg = arrp(look->boxIndex,look->activeBox,SEG) ;
  graphBoxDraw (look->activeBox, WHITE, RED) ;
  strncpy(look->positionBuffer, messprintf("%d", 
      arrp( look->boxIndex, box, SEG)->itv.len), 9) ;
  graphBoxDraw(look->positionBox, BLACK, YELLOW) ;
}

   /* color the corresponding bands in the other digestions
      and the corresponding piece of the dna 
      */
static void gelCorrespondingBands(LOOKGEL look, Array boxIndex, ENZ *masterEnzp ,
				  KEY seq, int min, int max)
{
  int i, activeBox = look->activeBox, e, f ;
  SEG *segp = 0 ;

  for (i = minLiveBox; i < maxLiveBox; i++)
    { 
      if (i < arrayMax(boxIndex) )
	segp = arrp(boxIndex, i, SEG) ;
      if (i != activeBox)
	graphBoxDraw(i, BLACK, GELCOLOR[segp->color]) ;
      if (!segp->enzp ||
	  segp->enzp == masterEnzp || 
	  segp->enzp->seqKey != seq  )
	continue ;
      f = segp->itv.from + segp->enzp->from ;
      e = segp->itv.from + segp->itv.len ;
      if (( f<min && e > min) ||
	  ( f >=min && f < max) )
	graphBoxDraw(i, BLACK, FRIEND_COLOR) ;
    }
}

	  
static void gelFollow (LOOKGEL look, double x, double y)
{ 
  SEG *seg = arrp(look->boxIndex,look->activeBox,SEG) ;
  KEY seq = seg->enzp->seqKey ;
  int origin = seg->enzp->from ;
  int from = origin + (seg->itv).from ;
  int len = (seg->itv).len ;
  char col = GELTINT[seg->color] ;
  char *ip ;
  KEY seqKey ;
  Array dna,  colors ;
  void *fMapLook ;
 
  if (seg->enzp->observed) return ;
  gelCorrespondingBands(look, look->boxIndex, seg->enzp, seq, from, from + len) ;
  if(!fMapActive(&dna,&colors,&seqKey, &fMapLook)
     || seq != seqKey || !graphActivate(fMapActiveGraph ()))
    { messout
	( messprintf ("Please select the dna window: %s, thank you.",
		      name(seq)) ) ;
      return ;
    }
  if ( len < 0 || from < 0 || from+len > arrayMax(colors))
    messerror( " Length confusion in gelFollow, sorry") ;
  else
    {
      ip = arrp(colors, from, char) ;
      while (len--)
	*ip++ |= col ;
      if (mapColSetByName ("Summary bar", -1))
	fMapReDrawDNA(fMapLook) ;
      else  /* need to open the summary bar */
	{ mapColSetByName ("Summary bar", TRUE) ;
	  fMapDraw (fMapLook, 0) ;
	}
    }
  graphActivate(look->graph) ;
}

/***************************************************************************/
/***************************************************************************/

void gelComparativeDisplay (KEY clone, KEY link) 
{ int box , i, j, k;
  LOOKGEL look ;
  ENZ *enzp ;
  KEYSET ks = clone ? queryKey (clone," > Gel") : 0 ;
  OBJ obj ;
  KEY key ;
  float x, y ;
  Array itvs ;

  gelDisplay() ; 
   if (!graphAssFind (&MAGIC, &look))
     { messerror ("gelComparative display can't find look") ;
       return ;
     }

  if (!keySetMax (ks)  && (obj = bsUpdate (clone))) 
    { lexaddkey ("HindIII", &key, _VMotif) ;
      bsAddKey (obj, _Gel, key) ;
      lexaddkey ("EcoRV", &key, _VMotif) ;
      bsAddKey (obj, _Gel, key) ;
      lexaddkey ("BamHI", &key, _VMotif) ;
      bsAddKey (obj, _Gel, key) ;
      bsSave (obj) ;
      messout ("Observed restriction data are missing, please update Clone %s", 
	       name(clone)) ;
      keySetDestroy (ks) ;
      ks = queryKey (clone," > Gel") ;
    }

  obj = bsCreate (clone) ; 
  for (i = 0 , j = 0 ; i <  keySetMax (ks) ; i++) 
    { /* compute from assembly */
      enzp = arrayp(look->enzymes, j++ ,ENZ) ; 
      enzp->text = messalloc(80) ;
      strncpy (enzp->text, name(keySet(ks, i)), 79) ;
      enzp->itvs = arrayCreate(80, INTERVAL) ;
      box =  look->enzymeBox  + 2 * (j - 1) + 1 ; 
      gelPick(box, (float)(-1), 0) ; 
      look->currentEntryBox = look->enzymeBox - 1 + 2 * j ;
      gelComputeGel(look) ;
      
      /* experimental values */
      bsGoto (obj, 0) ;
      if (bsFindKey (obj, _Gel, keySet (ks, i)) &&
	  bsPushObj (obj) &&
	  bsGetData (obj, _Bands, _Float, &x))
	{
	  enzp = arrayp(look->enzymes, j++, ENZ) ; 
	  enzp->seqKey = 0 ;
	  enzp->from = 0 ;
	  enzp->length = (enzp - 1)->length ;
	  enzp->observed = TRUE ;
	  
	  enzp->text = messalloc(80) ;
	  strncpy (enzp->text, "observed", 79) ;
	  itvs = enzp->itvs = arrayCreate(80, INTERVAL) ;
	  k = 0 ;
	  do 
	    { array(itvs, k, INTERVAL).from = 0 ;
	      array(itvs, k, INTERVAL).len = x ;
	      y = 3.5 * ( 4 - log10((double)(x)) ) ; 
	      array(itvs, k++,INTERVAL).log = (y > 0 ? y : 0 ) ;
	    } while (bsGetData (obj, _bsRight, _Float, &x)) ;
	  arraySort(itvs,intervalOrder) ;
	}
    }
  bsDestroy (obj) ;

      /* one extra line */
  enzp = arrayp(look->enzymes, j++, ENZ) ; 
  enzp->text = messalloc(80) ;
  enzp->itvs = arrayCreate(80, INTERVAL) ;

  box =  look->enzymeBox  + 2 * (j - 1) + 1 ;
  gelPick(box, (float)(-1), 0) ; 

  gelDraw (look, 0) ;
}

/***************************************************************************/
/***************************************************************************/

