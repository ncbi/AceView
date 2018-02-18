/*  file: pmapdisp.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *      Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: physical map display
 * Exported functions: pMapDisplay()
 * HISTORY:
 * Last edited: Dec 22 11:39 1998 (fw)
 *      - Added the option to toggle the display of loci.
 * * Feb  1 12:51 1997 (rd)
 * * Sept22 1996 (mieg): Removed Neil's code, see in older archives if needed
 * * Jun 23 23:19 1992 (rd): change GRIDMAP to handle intervals
 * * Jan 13 01:30 1992 (mieg): Introduced look->mag, suppressed globals showFacXY
                               Horizontal midButton zooming copied from gMap
 * * Jan  6 18:48 1992 (rd): change probe localization to use gridDisp calls
 * * Jan  6 14:50 1992 (rd): add horizontal scale change to menu
 * * Dec 16 14:32 1991 (rd): changed hybridisation -> map rules
 * * Dec  1 18:10 1991 (mieg): added intersection with keyset and FLAG_HIDE
 * * Oct 22 18:54 1991 (mieg): adjusted size of gMapBox 
 * * Oct 21 13:46 1991 (mieg): factor 20 in follow sequence
 * * Oct 18 20:38 1991 (mieg): made seqbox >=2,
 added a sequence box if the exists seq(gene(clone))
 * * Oct 14 22:02 1991 (rd): draw mini contig with scroll box
 * * Oct 14 18:45 1991 (rd): use TextFit window instead of MapScroll
 * * Oct  7 22:20 1991 (rd): changed cursor to redraw properly
 *-------------------------------------------------------------------
 */

/* $Id: pmapdisp.c,v 1.3 2017/03/18 15:30:57 mieg Exp $ */
/*
#define ARRAY_CHECK
#define MALLOC_CHECK
*/

#include "acedb.h"

#include "display.h"
#include "pmap.h"
#include "pmap_.h"
#include "fingerp.h"		/* for fpDisplay */

#include "key.h"
#include "lex.h"
#include "bs.h"
#include "plot.h"
#include "a.h"
#include "sysclass.h"
#include "systags.h"
#include "classes.h"
#include "tags.h"
#include "session.h"
#include "grid.h"
#include "client.h"
#include "help.h"

#include "freeout.h"

#include "bump.h"
#include "bump_.h"		/* provide access to internal members
				   of BumpStruct - YUCK */

/************************************************************/

static void pMapMenu(KEY k) ;

static magic_t GRAPH2PhysMap_ASSOC = "GRAPH2PhysMap";
/* use find the PhysMap pointer for the current graph */

static magic_t PhysMap_MAGIC = "PhysMap_MAGIC";
/* verify the received pointer */

/*..... ..... ..... .....static routines
*/

static PhysMap currentPhysMap(char *callerFuncname);

static void pMapDraw(PhysMap look1, KEY key) ; 
static void pMapBoxDraw(int box, DrawType drawType, SEG *seg);

static void pMapRetitle (PhysMap look) ;
static void pMapDestroy(void);
static void pMapPick(int box, double x, double y);
static void pMapMiddleDown (double x , double y) ;
static void pMapMiddleDrag (double x, double y); 
static void pMapMiddleUp (double x, double y);

static void findLimits(PhysMap look);

static void mapKbd (int k) ;
static void pMapShowAll(void) ;
static void pMapHide(void) ;

static void setDisplayDepths (void) ;
static void pMapShowSegs (PhysMap look, KEY key) ;
static void pMapSelectBox (PhysMap look, int box) ;
static void pMapFollow(PhysMap look, double x, double y);
static void pMapScroll (float *px, float *py, BOOL isDone) ;
static void pMapHighlight(KEYSET k, KEY key);
static void pMapResize (void);
static void drawMiniContig(PhysMap look, float fac);
static void displayAllBuriedClones(PhysMap look) ;


/*..... ..... ..... .....stray extern
*/
extern void pMapToGMap(KEY contig ,KEY key, int position); /*<<--code for this is in gmapdisp.c*/
static int gMapBox=0;

#define ModelToGraphX(__x) (Xmagn*(__x)-Xoffset)
#define GraphToModelX(__x) (((__x)+Xoffset)/Xmagn)
#define ModelToGraphY(__y) (Ymagn*(__y)-Yoffset)
#define ModelToOviewX(__x, __magn) ((((__magn)!=0.0)? (__magn): oviewXmagn)*(__x)-oviewXoffset)
#define GraphToOviewX(__x) ModelToOviewX(GraphToModelX(__x), 0.0)

#define ScrolledCentreOffset(__magn) (look->centre*(__magn)-look->winWidth/2)

/*..... ..... ..... .....anchor points used by pMapShowSegs - these may be modified in other places now
*/
static int
  anchorRule,
  anchorCueBar,
  anchorProbes,
  anchorYACs,
  anchorClones,
  anchorSequence,
  anchorNavigBar,
  anchorGenes,
  anchorRemarks,
  anchorMiniContig,
  anchorScrollHelp;

static int nReg = 10 ; /*depth of display for cosmids (enlarged to show buried clones)*/

/*..... ..... ..... .....
*/
#define DrawLine(__x1, __y1, __x2, __y2) graphLine(ModelToGraphX(__x1), ModelToGraphY(__y1), ModelToGraphX(__x2), ModelToGraphY(__y2))
#define DrawRectangle(__x1, __y1, __x2, __y2) graphRectangle(ModelToGraphX(__x1), ModelToGraphY(__y1), ModelToGraphX(__x2), ModelToGraphY(__y2))
#define DrawFilledRectangle(__x1, __y1, __x2, __y2) graphFillRectangle(ModelToGraphX(__x1), ModelToGraphY(__y1), ModelToGraphX(__x2), ModelToGraphY(__y2))

#define MENU_QUIT 1
#define MENU_HELP 2
#define MENU_PRINT 3
#define MENU_PRESERVE 4

#define MENU_RECALCULATE 6
#define MENU_SHOW_REMARKS 7
#define MENU_SHOW_WORK_REMARKS 8
#define MENU_SHOW_SELECTED 9
#define MENU_HIGHLIGHT_SELECTED 10
#define MENU_REVERT_DISPLAY 11
#define MENU_SHOW_ALL_BURIED_CLONES 49
#define MENU_SHOW_LOCUS 55
#define MENU_SHOW_SCROLLBAR 56

#define MENU_SET_ZOOM 15
#define MENU_SET_DISPLAY_DEPTHS 67


static float Xmagn=0.4, oviewXmagn=0.05, Ymagn=1.3, Xoffset=0.0, oviewXoffset=0.0, Yoffset=(-1.0);

static void  DrawSegText (SEG *seg, float y)
{  
  /* float 
    xL=seg->x0, xR=SegRightEnd(seg), dx = Xmagn*(xR - xL) ; */ 	
  char *text = 0 ;
  float x = SegMidPt(seg) ;
  
  if (IsSet(seg, FLAG_CANONICAL))
    text=messprintf("%s *", name(seg->key));
  else
    text = name(seg->key) ;
/* RD 961024 - incorrect to look at length of clone
** should look at whether overlaps next entry on this line
** much more complicated - too complicated for me now
  if (dx < strlen(text) - 3)
    text = messprintf("%3s-", name(seg->key)) ;
  if (dx < 3)
    text = messprintf("%c-", *name(seg->key)) ;
  if (dx < 2)
    text = 0 ;
*/
  if (text) 
    graphText (text,  ModelToGraphX(x) - strlen(text) / 2.0 ,   
	       ModelToGraphY (y) - 0.88) ;
}


static double oldy , oldDx, oldx;
static BOOL dragFast;
static PhysMap lookDrag=0;

/*..... ..... ..... .....
*/
static KEYSET pMapPendingKeySet=0;

/*..... ..... ..... .....default physical map menu
*/
static FREEOPT pmap_menu[]=
{
  {15 /*no of entries - **NB update this when change menu*/, "Physical Map"},

  {MENU_QUIT, "Quit"},
  {MENU_HELP, "Help"},
  {MENU_PRINT, "Print"},
  {MENU_PRESERVE, "Preserve"},
    
  {MENU_RECALCULATE, "Recalculate"},
  {MENU_SHOW_REMARKS, "Show/Hide All Remarks"},
  {MENU_SHOW_WORK_REMARKS, "Show/Hide Work Remarks"},
  {MENU_SHOW_SELECTED, "Show Selected Objects"},
  {MENU_HIGHLIGHT_SELECTED, "Highlight Selected Objects"},
  {MENU_REVERT_DISPLAY, "Unhighlight, Revert"},
    
  {MENU_SHOW_ALL_BURIED_CLONES, "Show/Hide Buried Clones"},
  {MENU_SHOW_LOCUS, "Show/Hide Loci"},
  {MENU_SHOW_SCROLLBAR, "Show/Hide Scroll Bar"},
  {MENU_SET_ZOOM, "Set Zoom"},
  {MENU_SET_DISPLAY_DEPTHS, "Set display depths"}

};

/************************************************************/

static PhysMap currentPhysMap (char *caller)
{
  PhysMap pmap;

   if (!(graphAssFind(&GRAPH2PhysMap_ASSOC, &pmap)))
    messcrash("%s() could not find PhysMap on graph", caller);
  if (!pmap)
    messcrash("%s() received NULL PhysMap pointer", caller);
  if (pmap->magic != &PhysMap_MAGIC)
    messcrash("%s() received non-magic PhysMap pointer", caller);
  
  return pmap; 
} /* currentPhysMap */


BOOL pMapDisplay (KEY key, KEY from, BOOL isOldGraph) /*called by display()/display.c as a functional*/
{
  KEY contig = 0, clone = 0, pmap ;
  int pos ;
  OBJ  Clone = 0 ;
  BOOL isProbe = FALSE ;
  Array segs = 0 ;
  PhysMap look = 0 ;
  static int isFirst = TRUE ;

  if (isFirst)
    { isFirst = FALSE ;
      if (getenv ("ACEDB_CONTIG9"))
	++pmap_menu->key ;
    }				/* end of RMD change */

  
   if (externalServer) 
     goto abortClone ; /* pmapconvert is much too slow */

  if (class(key) == _VClone)	/* use the contig */
    { clone = key ;

      if (!(Clone = bsCreate(clone)))
	goto abort ;

      if (!bsGetKey(Clone, _pMap, &contig) &&
	  !CanonicalFromBuriedClone(Clone, clone))
	{
	  if (CloneIsGridded(Clone, key))
	    { 
	      bsDestroy (Clone) ;
	      display (key, clone, GRID) ; /* key is a Grid_Clone */
	      goto abort ;
	    }
	  else
	    goto abortNoMap ;
	}
      if (clone == from)
	goto abortClone ;
      if (!contig)
	{
	  bsDestroy (Clone) ;
	  if (!(Clone = bsCreate(clone))) /*value of clone may have changed*/
	    goto abort ;
	  if(!bsGetKey (Clone,_pMap,&contig))
	    goto abortNoMap ; /*even parent not located on a contig*/
	}
      if (!bsGetData(Clone,_bsRight,_Int,&pos))
	goto abortNoMap ;

      bsDestroy (Clone) ;
      if (isProbe)
	clone = key ;
      key = contig ;
    }

  if (class(key) == _VpMap)
    pmap = key ;
  else if (class(key) == _VContig)
    { /* protection against ill defined contigs */
      OBJ Contig = bsCreate (key) ;
      int x1 = 0, x2 = -1 ;
      if (Contig && bsGetData (Contig, _pMap, _Int, &x1))
	bsGetData (Contig, _bsRight, _Int, &x2) ;
      bsDestroy (Contig) ;
      if ( x1 > x2)
	goto abortNoMap ;
      lexaddkey (name(key), &pmap, _VpMap) ;
    }
  else
    goto abort ;

  if (!(segs=arrayGet(pmap, SEG, segFormat)) && (class(key)!=_VContig || !(segs=pMapConvert(NULL, key, pmap))))
    goto abort ;

  if (isOldGraph && 
      graphAssFind (&GRAPH2PhysMap_ASSOC, &look) &&
      look &&
      look->magic == &PhysMap_MAGIC)
  {
    arrayDestroy(look->segs);
  }
  else 
  {
    if (!displayCreate(PMAP))
      goto abort;

    graphRegister(DESTROY, pMapDestroy) ;
    graphRegister(RESIZE, pMapResize);
    graphRegister(PICK, pMapPick) ;
    graphRegister(MIDDLE_DOWN, pMapMiddleDown) ;
    graphRegister(MIDDLE_DRAG, pMapMiddleDrag) ;
    graphRegister(MIDDLE_UP, pMapMiddleUp) ;
    graphRegister(KEYBOARD, mapKbd) ;

    look = (PhysMap) messalloc(sizeof(struct PhysMapStruct));
    look->magic = &PhysMap_MAGIC;
    look->mainGraph = graphActive() ;
    look->box2seg = arrayCreate(32,SEG*) ;
  }

  look->key = pmap ;
  look->segs = segs ;
  look->menu = pmap_menu ;
      /*always use the default menu when showing a contig for the first time*/
  look->activebox = 0 ;
  look->centre = 0 ;

  graphFreeMenu(pMapMenu, look->menu);

  graphAssociate(&GRAPH2PhysMap_ASSOC, look); /*doesn't matter if already been associated*/

  pMapRetitle(look) ;

  if ((!clone && class(from)==_VClone) || 
      class(from)==_VCalcul)  /* Geometrical positioning */
    clone=from;
  findLimits(look);

  if (isGifDisplay)
    { 
      Set (look, LOOK_NO_SCROLL_BAR) ;
      Set (look, LOOK_NO_LOCUS) ; 
      Set (look, LOOK_NO_REMARKS) ;
    }

  if (pMapPendingKeySet)
    {
      pMapHighlight (pMapPendingKeySet, clone); /* calls pMapDraw */
      pMapPendingKeySet = 0 ;
    }
  else
    pMapDraw(look, clone);
  return TRUE;

abortNoMap:
abortClone:
  bsDestroy(Clone);
  display(key, from, TREE);	/* just show as a TREE */

abort:
/*I don't create the look until later and this messfree shd be redundant now:*/
  messfree (look);
  return FALSE;
}

/************************************************************/

static void pMapDestroy(void)
{
  PhysMap look = currentPhysMap("pMapDestroy");

  arrayDestroy(look->segs);
  arrayDestroy(look->box2seg);
  look->magic=0;
  messfree(look);
}

/*************************************************************/

static void pMapRetitle (PhysMap look)
{ 
  KEY pmap = look->key, contig, map ;
  float start, stop ;
  OBJ Contig ;
  char *cp = 0 ;

  if (!lexReClass(pmap, &contig, _VContig))
    return ;
  if ((Contig = bsCreate(contig)))
    { if ( bsGetKey (Contig, _Map, &map) &&
	   bsPushObj (Contig) &&
	   bsGetData (Contig, _Left, _Float, &start) &&
	   bsGetData (Contig, _Right, _Float, &stop) )
	cp = messprintf ("%s: [%4.2f, %4.2f]", name(map), start, stop) ;
      bsDestroy(Contig) ;
    }

  if (!cp)
    cp = messprintf ("Contig %s (Not localised)", name(pmap)) ;
  graphRetitle (strnew (cp, graphHandle())) ;
}

/*..... ..... ..... ..... ..... ..... ..... ..... ..... ..... ..... .....menu
*/
static void pMapMenu(KEY k)
{
  PhysMap look = currentPhysMap("pMapMenu");

  switch ((int) k)
  {
  /*.....display*/
  case MENU_QUIT:
    graphDestroy () ;
    break ;
  case MENU_HELP: 
    helpOn("Physical-map"); 
    break;
  case MENU_PRINT: 
    graphPrint(); 
    break;
  case MENU_PRESERVE: 
    displayPreserve(); 
    break;
  case MENU_RECALCULATE:
    pMapRecalculateContig(look, TRUE);  
    pMapDraw(look, 0); /*mhmp 07.12.98, avoids a bizare crash */
    break ;
  case MENU_SHOW_REMARKS:
    Toggle(look, LOOK_NO_REMARKS);
    pMapDraw(look, 0);
    break ;
  case MENU_SHOW_WORK_REMARKS:
    Toggle(look, LOOK_WORK_REMARKS);
    pMapDraw(look, 0);
    break ;
  case MENU_SHOW_SELECTED: 
    pMapHide(); 
    break;
  case MENU_HIGHLIGHT_SELECTED: 
    pMapHighlight(NULL, 0); 
    break;
  case MENU_REVERT_DISPLAY: 
    pMapShowAll(); 
    break;
  case MENU_SHOW_ALL_BURIED_CLONES: /*switch*/
    Toggle(look, FLAG_SHOW_BURIED_CLONES);
    if (IsSet(look, FLAG_SHOW_BURIED_CLONES))
      displayAllBuriedClones(look);
    pMapDraw(look, 0);
    break;
  case MENU_SET_ZOOM:
    if (messPrompt("Increase to zoom in, decrease to zoom out:", messprintf("%.2f",Xmagn), "fz"))
    {
      float old = Xmagn ; /* mhmp 04.12.98 pour eviter Xmagn = 0 */
      freefloat(&Xmagn);
      if (!Xmagn)
	Xmagn = old ;
      else
	pMapDraw(look, 0);   pMapDraw(look, 0);
    }
    break;
  case MENU_SET_DISPLAY_DEPTHS:
    setDisplayDepths () ;
    break ;
  case MENU_SHOW_LOCUS:
    Toggle(look, LOOK_NO_LOCUS);
    pMapDraw(look,0);
    break;
  case MENU_SHOW_SCROLLBAR:
    Toggle(look, LOOK_NO_SCROLL_BAR);
    pMapDraw(look,0);
    break;

  }
  return;
}

static void pMapZoomIn (void)
{
  PhysMap look = currentPhysMap("pMapZoomIn") ;
  Xmagn *= 1.5 ;
  pMapDraw(look, 0);
}

static void pMapZoomOut (void)
{
  PhysMap look = currentPhysMap("pMapZoomOut") ;
  Xmagn /= 1.5 ;
  pMapDraw(look, 0);
}


/*********************************************************************************************************************/

static void pMapPick(int box, double x, double y) 
{
  PhysMap look = currentPhysMap("pMapPick");

  if (!box) /*pick on background*/
    pMapSelectBox(look, 0); /*deselect current selection*/
  else if (box==look->scrollBox)
    graphBoxDrag (look->scrollBox, pMapScroll);
  else if (arrayMax(look->box2seg)<box)
    messerror("pMapPick received an invalid box number");
  else if (box==look->activebox)  /*second pick on this box: display the object*/
    pMapFollow(look, x, y);
  else /*select this box and highlight it*/
    pMapSelectBox (look,box) ;

  return;
}

#undef X_ARM
#undef Y_ARM
#undef ShiftLeftSegment
#undef DrawXorCross
#undef DrawXorMark
#undef DragXorCrossTo
#undef DragXorMark
#undef ClearModalDragOp

/***********************************/

static void mapKbd (int k)
{
  int box, table ;
  SEG *seg ;
  PhysMap look = currentPhysMap("mapKbd") ;

  if (!look->activebox)
    return ;

  box = look->activebox ;
  if (!(seg = arr(look->box2seg,box,SEG*)))
    return ;
  table = class(seg->key) ;
  switch (k)
    {
    case LEFT_KEY :
      while ( --box > 0 &&
	     (seg = arr(look->box2seg,box,SEG*)) &&
	     class(seg->key) != table) ;
      break ;
    case RIGHT_KEY :
      while (++box < arrayMax(look->box2seg) &&
	     (seg = arr(look->box2seg,box,SEG*)) &&
	     class(seg->key) != table) ;
      break ;
    }

  if (box < arrayMax(look->box2seg) && arr(look->box2seg, box, SEG*))
    { pMapSelectBox (look, box) ;
      pMapFollow (look,0,0) ;
    }
}

/**********************************/

static void findLimits (PhysMap look)
{
  int i, min=1000000000, max= -1000000000;
  Array segs=look->segs;
  SEG *seg;

  for (i=1; i<arrayMax(segs); ++i)
  {
     seg=arrp(segs, i, SEG);
     if (seg->x0<min) min=seg->x0;
     if (max<seg->x1) max=seg->x1;
  }
  look->min=min;
  look->max=max;
  return;
}

/******************************/

static void pMapDraw (PhysMap look1, KEY key)     /* key becomes activebox, unless reverting*/
{
  float y=0.0;
  PhysMap look = currentPhysMap ("pMapDraw") ;

  if (look != look1) invokeDebugger () ;
  graphClear();
  graphLinewidth(0.10);

  if (!arrayMax(look->segs))
    return;

  graphFitBounds(&look->winWidth, &look->winHeight);
  ++look->winWidth;


  graphLine(look->winWidth/2, y, look->winWidth/2, y+0.9);
  pMapShowSegs (look, key) ;
  graphLine(look->winWidth/2, look->winHeight-2.0, 
	    look->winWidth/2, look->winHeight-2.0+0.9);
  graphRedraw();

  return;
}


/*************************************************************/
/******** pmap show intersect with active keyset *************/
 
static void pMapShowAll(void) /*show pmap without any highlighting*/
{
  int i, n ;
  SEG *seg ;
 
  PhysMap look = currentPhysMap("pMapShowAll") ;

  seg = arrp(look->segs,0,SEG) ;
  i = arrayMax(look->segs) ;
  n = FLAG_HIDE | FLAG_HIGHLIGHT ; 
  while(i--)
    {
      seg->flag&= ~n; /*<<--reset these bits*/
      seg++;
    }

  Clear (look, LOOK_WORK_REMARKS);
    /*make sure extra remarks aren't shown*/
  Set(look, FLAG_REVERT_DISPLAY);
    /*<<--this is just to communicate down the call stack without passing arg: it only happens here*/
  pMapDraw (look, 0);
  Clear(look, FLAG_REVERT_DISPLAY);
}

/**********/

static void pMapHide(void)
{
  int i, j, n ;
  SEG *seg ;
  void *dummy ;
  KEYSET kSet = 0 ;
 
  PhysMap look = currentPhysMap("pmapHide") ;

  if(!keySetActive(&kSet, &dummy))
    { messout("First select a keySet window, thank you.") ;
      return ;
    }
  seg = arrp(look->segs,0,SEG) ;
  i = arrayMax(look->segs) ;
  n = FLAG_HIDE ;
  while(i--)
    { seg->flag &= ~n ;  /* bitwise NOT */
      if (!arrayFind(kSet,&seg->key,&j, keySetOrder) &&
	  !arrayFind(kSet,&seg->parent,&j, keySetOrder))
	seg->flag |= FLAG_HIDE ;
      seg++ ;
    }
  pMapDraw (look, 0);
}


/*********************************************/
/******* local display routines **************/

/************* first a mini package for bumping text *************/

static BOOL bumpMayWrite (BUMP bump, char *text, float x, float yy, BOOL doIt)
{ int y = 0 ;
  int n ; float dy = 1000000000.0 ; /* RD 961024 was 6.0 */

  x =  ModelToGraphX (x) -3.5 ;
  n = bumpText (bump, text, &y, &x, dy, FALSE);
  if (!n || n > 999) return FALSE ;

  if (doIt)
    { 
      char buf[1000] ;

      strncpy(buf, text, n) ; buf[n] = 0 ;
      graphText(buf,x,ModelToGraphY(y + yy)) ;
      /*
	graphText (text, ModelToGraphX(bump->end[min])-3.5, 
	ModelToGraphY(bump->y0 + min)) ;
	*/
    }
  return TRUE ;
}


/********* use middle button for cursor **********/


#define ACCEL (1./25.)

/*..... ..... ..... .....
*/
static void pMapMiddleDown (double x, double y) 
{ 
  PhysMap look = currentPhysMap("pMapMiddleDown"); lookDrag=look;

  oldx=x; oldy=y;

 if (anchorMiniContig < y + 8 && anchorMiniContig > y - 8)
  { float x1, x2, y1, y2;
    dragFast = TRUE ;
    graphBoxDim(look->scrollBox, &x1, &y1, &x2, &y2);
    oldDx=0.5*(x2-x1);
    oviewXoffset=ScrolledCentreOffset(oviewXmagn);
    graphXorLine(oldx-oldDx, anchorMiniContig-10, oldx-oldDx, anchorMiniContig+10);
    graphXorLine(oldx+oldDx, anchorMiniContig-10, oldx+oldDx, anchorMiniContig+10);
      /*
      anchorMiniContig is set in pMapShowSegs when the graph is constructed depending 
      on whether this is extendedPanel or not, so this covers both cases
      */
  }
  else if (IsSet(look, FLAG_SHOW_BURIED_CLONES))
  { dragFast = FALSE ;
    graphXorLine(oldx, anchorMiniContig+10, oldx, anchorMiniContig+50);
  }
  else /*default physical map display*/
  { dragFast = FALSE ;
    graphXorLine(oldx, 0, oldx, anchorMiniContig-10);
  }
  return;
}

/*..... ..... ..... .....
*/
static void pMapMiddleDrag (double x, double y) 
{
  PhysMap look = currentPhysMap("pMapMiddleDrag");
  if (dragFast)
  {
    graphXorLine(oldx-oldDx, anchorMiniContig-10, oldx-oldDx, anchorMiniContig+10);
    graphXorLine(oldx+oldDx, anchorMiniContig-10, oldx+oldDx, anchorMiniContig+10);
  }
  else if (IsSet(look, FLAG_SHOW_BURIED_CLONES))
  {
    graphXorLine(oldx, anchorMiniContig+10, oldx, anchorMiniContig+50);
  }
  else /*default physical map display*/
  {
    graphXorLine(oldx, 0, oldx, anchorMiniContig-10);
  }
  oldx=x;
  if (dragFast)
  {
    oldDx*=exp((oldy-y)*ACCEL);
    graphXorLine(oldx-oldDx, anchorMiniContig-10, oldx-oldDx, anchorMiniContig+10);
    graphXorLine(oldx+oldDx, anchorMiniContig-10, oldx+oldDx, anchorMiniContig+10);
  }
  else if (IsSet(look, FLAG_SHOW_BURIED_CLONES))
  {
    graphXorLine(oldx, anchorMiniContig+10, oldx, anchorMiniContig+50);
  }
  else
  {
    graphXorLine(oldx, 0, oldx, anchorMiniContig-10);
  }
  oldy=y;
  return;
}

/*..... ..... ..... .....
*/
static void pMapMiddleUp (double x, double y) 
{ 
  float x1, x2, y1, y2 ;
  PhysMap look = currentPhysMap ("pMapMiddleUp") ;
  

  if (lookDrag !=look)
    { invokeDebugger(); return ; }

  x += .3 ; /* mieg */
  if (dragFast)
  {
    graphBoxDim(look->scrollBox, &x1, &y1, &x2, &y2);
    look->centre+=(x-look->winWidth*0.5)/oviewXmagn;
    Xmagn*=0.5*(x2-x1)/oldDx;
  }
  else
    look->centre+=(x-look->winWidth*0.5)/Xmagn;
  pMapDraw(look, 0);
  return;
}

/*..... ..... ..... .....
*/
static void pMapBoxDraw(int box, DrawType drawType, SEG *seg)
{
  if (seg)
    switch(drawType)
      {
      case HIGHLIT:
	graphBoxDraw(box, BLACK, MAGENTA); break;
      case SELECTED:
	graphBoxDraw(box, BLACK/*WHITE*/, RED); break;
      case DEBUG:
	graphBoxDraw(box, BLACK, GREEN); break;
      case SISTER:
	graphBoxDraw(box, BLACK, LIGHTRED); break;
      case NORMAL:
	if (class(seg->key) == _VSequence || IsSet(seg, FLAG_SEQUENCED))
	  graphBoxDraw(box, BLACK, YELLOW);
	else
	  graphBoxDraw(box, BLACK, WHITE); 
	break;
      case BURIED_CLONE:
	graphBoxDraw(box, BLACK, LIGHTGRAY); break;
      default:
	break;
      }
  else if (box==gMapBox)
    graphBoxDraw(box, BLACK, drawType==SELECTED? GREEN: LIGHTGREEN);
  
  return;
}

/*..... ..... ..... .....
*/
static void pMapScroll (float *px, float *py, BOOL isDone)
{
  PhysMap look = currentPhysMap("pMapScroll") ;

  *py = anchorMiniContig - 0.2 ;
  if (isDone)
    { look->centre += (*px - look->winWidth*0.45) / oviewXmagn ;
      pMapDraw (look, 0);
    }
}

/*..... ..... ..... .....
*/

static void pMapResize (void)
{
  PhysMap look = currentPhysMap("pMapResize");
  pMapDraw(look, 0);
  return;
}

/*..... ..... ..... ..... ..... ..... ..... ..... ..... ..... ..... .....display
*/
static void drawMiniContig(PhysMap look, float fac)
{
  SEG *seg;
  float x, height;
  int i;
  int kHeight ;
  BUMP bump = bumpCreate(6,1) ;
#define DY 0.75

  Xoffset=ScrolledCentreOffset(Xmagn);
  oviewXoffset=ScrolledCentreOffset(oviewXmagn);

  if(IsntSet(look, LOOK_NO_LOCUS))
    {
      graphTextHeight(0.75);
      for (i=0 ; i < arrayMax(look->segs) ; ++i)
	{
	  seg=arrp(look->segs, i, SEG) ;
	  
	  if (class(seg->key)==_VLocus)
	    {
	      x=ModelToOviewX(SegMidPt(seg), fac) ;
	      if (x<0) continue;
	      if (x>look->winWidth) break;
	      kHeight = 0 ;
	      bumpItem (bump, 1, strlen(name(seg->key)), &kHeight, &x) ;
	      height = 1 + kHeight*DY ;
           /* graphLine(x, anchorMiniContig-height+DY, x, anchorMiniContig); */
	      graphText(name(seg->key), x,  anchorMiniContig-height);
	    }
	}
      graphTextHeight(0 /*back to default*/);
    }

  if(IsntSet(look, LOOK_NO_SCROLL_BAR))
    {
      graphFillRectangle(ModelToOviewX(look->min, fac), anchorMiniContig, ModelToOviewX(look->max, fac), anchorMiniContig+0.2);

      look->scrollBox=graphBoxStart();
      
      graphRectangle(GraphToOviewX(0), anchorMiniContig-0.2, GraphToOviewX(look->winWidth), anchorMiniContig+0.4);
      graphBoxEnd();
      graphBoxDraw(look->scrollBox, BLACK, GREEN);
    }
  bumpDestroy (bump) ;
  return;
}

/*..... ..... ..... .....
*/

static void fingerPrint(int box)  /* mieg janv 94 */
{ SEG* seg ;

  PhysMap look = currentPhysMap("fingerPrint");

   seg = array(look->box2seg, box, SEG*) ;
      if(seg && seg->key)
	fpDisplay(seg->key) ;
}

static MENUOPT cloneBoxMenu[] =  /* mieg janv 94 */
    { { (VoidRoutine)fingerPrint,"Finger-Print", }, 
      { 0, 0 } } ;

/*..... ..... ..... ..... ..... ..... ..... ..... ..... ..... ..... .....
*/
static int nYAC = 5 ;
/*  static int nReg = 10 ;	see up when buried */
static int nProbe = 3 ;

/*..... ..... ..... .....
*/

static void setDisplayDepths (void)
{ int old ; /* mhmp 04.12.98 pour eviter des nb de lignes <= 0 */
  if (messPrompt ("Numbers of lines for probes, virtual clones, fingerprinted clones:",
		  messprintf ("%d %d %d", nProbe, nYAC, nReg),
		  "iiiz"))
    { old = nProbe ;
      freeint (&nProbe) ;
      if (nProbe <= 0)
	nProbe = old ;
      old = nYAC ;
      freeint (&nYAC) ;
      if (nYAC <= 0)
	nYAC = old ;
      old = nReg ;
      freeint (&nReg) ;
      if (nReg <= 0)
	nReg = old ;
    }
}

/*..... ..... ..... .....
*/
static void pMapShowSegs(PhysMap look, KEY key)
{
  int  i, box, keybox=0 ;
  float boxRightEdge ;
  float xL, xR, oldxR, y;
  SEG  *seg ;
  char *text ;
  BUMP remBump,geneBump,probeBump ;
  int  iYac = 0, iReg = 0 ;
  int  nYAC=5;
  char boxInfoBuf[256] ;
  int graph_width ;

  graphFitBounds (&graph_width, 0) ;

  /*..... ..... ..... .....find out what arrangement required*/
  if (IsSet(look, FLAG_SHOW_BURIED_CLONES))
  {
    anchorCueBar=0; /*NONE*/
    anchorMiniContig=anchorCueBar+7;
    anchorScrollHelp=anchorMiniContig+1;
    anchorSequence=anchorScrollHelp+1;
    anchorProbes=anchorSequence+1+nProbe;
    anchorYACs=anchorProbes+nYAC;
    anchorClones=anchorYACs+1+nReg;
    anchorNavigBar=anchorRule=anchorClones+2;
    anchorGenes=anchorRule+1;
    anchorRemarks=anchorGenes+2;
  }
  else /*ordinary pmap display*/
  {
    anchorCueBar=0; /*NONE*/
    anchorProbes = nProbe ;
    anchorYACs=anchorProbes+nYAC ;
    anchorClones=anchorYACs+1+nReg;
    anchorSequence=anchorClones+1;
    anchorNavigBar=anchorSequence+1;
    anchorGenes=anchorNavigBar+1;
    anchorRemarks=anchorGenes+2;
    anchorMiniContig=look->winHeight-2;
    anchorScrollHelp=look->winHeight-1; /*make sure it's visible whatever the window size*/
  }

  /*..... ..... ..... .....*/
  probeBump=bumpCreate(nProbe, 0) ;
  geneBump=bumpCreate(2,0) ;
  remBump=bumpCreate(7, 0) ;
  /*..... ..... ..... .....*/

          /* Direct interpolating move from the gmap */
          /* The position is passed approximately in pmap units */

  if (class(key)==_VCalcul)
  {
    xL=(int) KEYKEY(key)-(int) 0x800000;
    look->centre=xL=Trunc(look->min, xL, look->max);
  }
  else if (key)
  {
    for (i=arrayMax(look->segs); --i; )
    {
      if (arrp(look->segs, i, SEG)->key==key)
      {
        look->centre=SegMidPt(arrp(look->segs, i, SEG));
	break;
      }
    }
  }
  if (!key && look->activebox && 
      (seg=arr(look->box2seg, look->activebox, SEG*)))
    key=seg->key;

  look->activebox=0;
  arrayMax(look->box2seg)=0;

  Xoffset=ScrolledCentreOffset(Xmagn);

  /*..... ..... ..... .....draw mini-contig for default display: takes care of buried clones display*/
    drawMiniContig(look, 0.0); 
    if(IsntSet(look, LOOK_NO_SCROLL_BAR))
      {
	graphText(
		  "Mid-mouse button: touch to recenter, drag vertically to zoom.",
		  12, anchorScrollHelp) ;
	
	graphButton("Zoom in", pMapZoomIn, 1, anchorScrollHelp ) ;
	graphButton("Zoom out", pMapZoomOut,graph_width - 9, anchorScrollHelp ) ;
      }
   /*draw green bar - navigation box to the gMap*/
  {
    array(look->box2seg, gMapBox=graphBoxStart(), SEG*)=0;
    graphRectangle(0, ModelToGraphY(anchorNavigBar+0.2), 
		   look->winWidth, ModelToGraphY(anchorNavigBar+0.8));
    graphTextHeight(0.6);
    for (i=10; i<look->winWidth; i+=40)
    {
      graphText("Pick here to go to the genetic map", 
		i, ModelToGraphY(anchorNavigBar+0.25));
    }
    graphTextHeight(0 /*back to default*/);
    graphBoxEnd();
    graphBoxDraw(gMapBox, BLACK, LIGHTGREEN);
  }
  /*..... ..... ..... .....*/
  oldxR = -10000001.0 ;
  arraySort(look->segs, pMapCompareSeg);

  for (i=1 /*seg[0] is a dummy so don't touch it*/; i<arrayMax(look->segs); ++i)
  {
    box = 0 ;
    seg = arrp(look->segs, i, SEG);

    if (IsSet(seg, FLAG_HIDE)) continue;
    xL=seg->x0; xR=SegRightEnd(seg);
    if (ModelToGraphX(xR)<0)
      /*off the left side of the window, but must ensure that the vertical placing of clones, remarks, genes etc. stays constant*/
    {
      if (class(seg->key) == _VClone)
	{
	  if (IsSet(seg, FLAG_FINGERPRINT))
	    {
	      if (oldxR<xL) iReg=0;
	      if (oldxR<xR) oldxR=xR;
	      if (IsSet(seg, FLAG_DISPLAY) || 
		  IsntSet(seg, FLAG_IS_BURIED_CLONE) || 
		  IsSet(look, FLAG_SHOW_BURIED_CLONES))
		/*FLAG_DISPLAY - display this buried clone*/
		{
		  ++iReg;
		}
	    }
	  else if (IsntSet(seg, FLAG_PROBE))
	    ++iYac ;
	}
      else if (class(seg->key) == _VLocus)
	{ graphTextFormat (BOLD) ;
	  bumpMayWrite (geneBump, name(seg->key), SegMidPt(seg), anchorGenes, FALSE);
	  graphTextFormat (PLAIN) ;
	}
      else if (class(seg->key) == _VText &&
	     !(IsSet(look, LOOK_NO_REMARKS) 
#ifdef ACEDB5
	       && IsSet(seg, FLAG_REMARK)
#endif
	       ) &&
	       !(IsntSet(look, LOOK_WORK_REMARKS) && IsSet(seg, FLAG_WORK_REMARK)))
	{
	  if (IsSet(seg, FLAG_MORE_TEXT))
	    text = messprintf("%s(*)", name(seg->key));
	  else
	    text = name(seg->key);
	  bumpMayWrite(remBump, text, SegMidPt(seg), anchorRemarks, FALSE);
	}
      else 
	continue; /*skip drawing*/
    }
    if (look->winWidth<ModelToGraphX(xL)) continue;	/* no need to be so careful */

    *boxInfoBuf = 0 ;		/* start from scratch */

    if (class(seg->key) == _VClone)
      { 

	if (IsSet(seg, FLAG_IS_BURIED_CLONE) && IsntSet(seg, FLAG_DISPLAY) && 
	    IsntSet(look, FLAG_SHOW_BURIED_CLONES)) ;
	else
	  { 
	    int bgColor = TRANSPARENT ;
	    
	    if (IsSet(seg, FLAG_IS_COSMID))
	      strcat (boxInfoBuf, "Cosmid ") ;
	    else if (IsSet(seg, FLAG_IS_YAC))
	      strcat (boxInfoBuf, "YAC ") ;
	    else if (IsSet(seg, FLAG_IS_CDNA))
	      strcat (boxInfoBuf, "cDNA ") ;
	    else if (IsSet(seg, FLAG_IS_FOSMID))
	      strcat (boxInfoBuf, "Fosmid ") ;

	    if (IsSet(seg, FLAG_FINGERPRINT))
	      {	strcat (boxInfoBuf, "Fingerprinted ") ;
		if (oldxR<xL)
		  { 
		    iReg=0; /* start new islands from bottom */
		    if (-10000000.0<oldxR)
		      {
			DrawLine(xL-7, anchorClones, xL-7, anchorClones+1);
			DrawLine(xL-8, anchorClones, xL-8, anchorClones+1);
		      }
		  }
				/* RD 961024 - start box after || */
		box = graphBoxStart() ;
		graphBoxMenu(box, cloneBoxMenu) ;  /* mieg janv 94 */
		
		if (oldxR<xR) oldxR=xR;
		y = (anchorClones+0.5) - (iReg%nReg) - 0.5*((iReg/nReg)%2) ;
		if (IsSet(seg, FLAG_COSMID_GRID))
		  { strcat (boxInfoBuf, "On_cosmid_grid ") ;
		    DrawFilledRectangle(xL, y, xR, y+0.1);
		  }
		else
		  DrawLine(xL, y, xR, y);
		DrawSegText(seg, y) ;
		if (IsSet(seg, FLAG_SEQUENCED))
		  { strcat (boxInfoBuf, "Sequenced") ;
		    bgColor = YELLOW ;
		  }
		++iReg ;
	      }
	    else if (IsSet(seg, FLAG_PROBE))
	      { box = graphBoxStart() ;
		bumpMayWrite(probeBump, name(seg->key), SegMidPt(seg), 
			     anchorProbes-nProbe, TRUE);
		strcat (boxInfoBuf, "Probe") ;
	      }
	    else /*YACs*/
	      { y = anchorYACs - (iYac%nYAC) + 0.5*((iYac/nYAC)%2) ;
		if (IsSet(seg, FLAG_YAC_GRID))
		  DrawFilledRectangle(xL, y, xR, y+0.1);
		else
		  DrawLine(xL, y, xR, y);
				/* RD 961024 - after line */
		box = graphBoxStart() ;
		text = messprintf ("(%s)", name(seg->key)) ;
		graphText(text, 
			  ModelToGraphX(SegMidPt(seg))-3.5, 
			  ModelToGraphY(y-0.4)); 
		/* DrawSegText (seg, y) ; */
		bgColor = WHITE ;
		++iYac ;
		strcat (boxInfoBuf, "YAC") ;
	      }
	    graphBoxEnd() ;
	    graphBoxDraw (box,BLACK, bgColor) ;
	    graphBoxDim (box, 0, 0, &boxRightEdge, 0) ;
	    if (boxRightEdge > 0)
	      graphBoxInfo (box, seg->key, *boxInfoBuf ? strnew(boxInfoBuf, graphHandle()): 0) ;
	  }
      }
    else if (class(seg->key) == _VSequence)
      {
	box=graphBoxStart();
	DrawRectangle(xL, anchorSequence+0.2, xR, anchorSequence+0.8);
	graphBoxEnd();
	graphBoxInfo (box,seg->key, 0) ; 
      }
    else if (class(seg->key) == _VAllele)
      {
	box=graphBoxStart() ;
	DrawLine(xL, anchorSequence+0.8, xL-3, anchorSequence+1.5) ;
	DrawLine(xL, anchorSequence+0.8, xL+3, anchorSequence+1.5) ;
	DrawLine(xL-3, anchorSequence+1.5, xL+3, anchorSequence+1.5) ;
	graphBoxEnd() ;
	graphBoxDim (box, 0, 0, &boxRightEdge, 0) ;
	if (boxRightEdge > 0)
	  graphBoxInfo (box, seg->key, *boxInfoBuf ? strnew(boxInfoBuf, graphHandle()): 0) ;
      }
    else if (class(seg->key) == _VLocus)
      {
	box=graphBoxStart();
	bumpMayWrite(geneBump, name(seg->key), SegMidPt(seg), anchorGenes, TRUE) ;
	graphBoxEnd();
	graphBoxDim (box, 0, 0, &boxRightEdge, 0) ;
	if (boxRightEdge > 0)
	  graphBoxInfo (box, seg->key, *boxInfoBuf ? strnew(boxInfoBuf, graphHandle()): 0) ;
      }
    else if (class(seg->key) == _VText && 
	     !(IsSet(look, LOOK_NO_REMARKS) 
#ifdef ACEDB5
	       && IsSet(seg, FLAG_REMARK)
#endif
	       ) &&
	     !(IsntSet(look, LOOK_WORK_REMARKS) && IsSet(seg, FLAG_WORK_REMARK)))
      {
	box = graphBoxStart() ;
	if (IsSet(seg, FLAG_MORE_TEXT))
	  text = messprintf("%s(*)", name(seg->key));
	else
	  text = name(seg->key) ;
	bumpMayWrite(remBump, text, SegMidPt(seg), anchorRemarks, TRUE);
	graphBoxEnd();
      }
    if (box) /*something drawn*/
      {
/*for showing clones matched against in region while calculating positionClone:*/
      if (class(seg->key) != _VClone ||  /* unset clones are better drawn transparent */
	  IsSet(seg, FLAG_HIGHLIGHT) )
	  pMapBoxDraw(box, IsSet(seg, FLAG_HIGHLIGHT)? HIGHLIT:NORMAL, seg);
      array(look->box2seg, box, SEG*)=seg;
      if(seg->key==key) keybox=box;
    }
  }

  if (keybox && IsntSet(look, FLAG_REVERT_DISPLAY))
    pMapSelectBox(look, keybox);

  /*..... ..... ..... .....*/
  bumpDestroy(remBump);
  bumpDestroy(geneBump);
  bumpDestroy(probeBump);
  return;
}

/*..... ..... ..... ..... ..... ..... ..... ..... ..... ..... ..... .....*/

static Array sisterBox=0, sisterSeg=0 ; /* disgustingly global - Richard */

static void makeSisterList (PhysMap look, int box)
{
  OBJ   obj ;
  SEG	*seg, *subseg ;
  Array b2s = look->box2seg ;
  int	i, n, pos ;
  static KEYSET sisterKeys = 0 ;
  static Array flat = 0 ;
  static KEY _Positive ;

  if (!_Positive)
    lexaddkey ("Positive", &_Positive, 0) ;
 
  seg = arr(b2s, box, SEG*) ;

  sisterKeys = keySetReCreate (sisterKeys) ;
  n = 0 ;
  keySet(sisterKeys, n++) = seg->parent ;

  if ((obj = bsCreate (seg->parent)))
    { flat = arrayReCreate (flat, 32, BSunit) ;
      if (bsFindTag (obj, _Positive) && bsFlatten (obj, 2, flat))
	for (i = 0 ; i < arrayMax(flat) ; i += 2)
	  keySet(sisterKeys, n++) = arr(flat, i+1, BSunit).k ;
      bsDestroy (obj) ;
    }
  keySetSort (sisterKeys) ;
  keySetCompress (sisterKeys) ;

  sisterBox = arrayReCreate (sisterBox, 8, int) ;
  sisterSeg = arrayReCreate (sisterSeg, 8, SEG*) ;
  n = 0 ;
  for (i = box ; --i > 1 && (subseg = arr(b2s, i, SEG*)) && SegMidPt(seg) - 250 < SegMidPt(subseg) ; )
    if (keySetFind (sisterKeys, subseg->parent, &pos))
      { array(sisterBox, n, int) = i ;
	array(sisterSeg, n++, SEG*) = subseg ;
      }
  for (i = box ; ++i < arrayMax(b2s) && (subseg = arr(b2s, i, SEG*)) && SegMidPt(subseg) < SegMidPt(seg) + 250 ; )
    if (keySetFind (sisterKeys, subseg->parent, &pos))
      { array(sisterBox, n, int) = i ;
	array(sisterSeg, n++, SEG*) = subseg ;
      }
}

/*********************************************************/

static void pMapSelectBox (PhysMap look, int box)
/* switch activebox */
{
  int i ;
  SEG *seg = 0, *sSeg = 0 ;

  if (look->activebox)		/* deselect the current selection */
    { seg = arr(look->box2seg, look->activebox, SEG*) ;
      if (seg && IsSet(seg, FLAG_HIGHLIGHT))
	pMapBoxDraw (look->activebox, HIGHLIT, seg) ;
      else
	pMapBoxDraw (look->activebox, NORMAL, seg) ;
      if (seg)			/* i.e. not for special boxes*/
	{ makeSisterList (look, look->activebox) ;
	  for (i = arrayMax(sisterBox) ; i-- ; )
	    { sSeg = arr(sisterSeg, i, SEG *) ;
	      if (IsSet (sSeg, FLAG_HIGHLIGHT))
		pMapBoxDraw (arr(sisterBox, i, int), HIGHLIT, sSeg) ;
	      else
		pMapBoxDraw (arr(sisterBox, i, int), NORMAL, sSeg) ;
	    }
	}
    }

  if (box)			/* select the given box */
    { seg = arr(look->box2seg, box, SEG*);
      pMapBoxDraw(box, SELECTED, seg);
      if (seg)			/*i.e. not for special boxes */
	{ makeSisterList (look, box) ;
	  for (i = arrayMax(sisterBox) ; i-- ; )
	    { sSeg = arr(sisterSeg, i, SEG *) ;
	      if (IsSet (sSeg, FLAG_IS_BURIED_CLONE))
		pMapBoxDraw (arr(sisterBox, i, int), BURIED_CLONE, sSeg) ;
	      else
		pMapBoxDraw (arr(sisterBox, i, int), SISTER, sSeg) ;
	    }
	}
    }

  look->activebox = box ;
  look->currentSelected= box && seg ? seg->key : 0 ;
}

/*********************************************************/

static void pMapFollow(PhysMap look, double x, double y)
{
  extern int sequenceLength(KEY seq);
  KEY  from ;
  SEG  *seg ;

  if (look->activebox == gMapBox)
    { KEY contig ;

      if (!lexword2key(name(look->key), &contig, _VContig))
	return ;
      pMapToGMap (contig, 0, 
	  (int)(look->centre + (x - look->winWidth/2)/Xmagn)) ;
    }
  else
    { seg = arr(look->box2seg,look->activebox,SEG*) ;
      
      if (class(seg->key) == _VLocus)
	display (seg->key, look->key, GMAP) ;
      else if (class(seg->key) == _VSequence)
	{
	  if (SegLen(seg)!=0)
	    from=KEYMAKE(_VCalcul, sequenceLength(seg->key)*x/(20.0*Xmagn*(0.5*SegLen(seg))/*ie old dx*/));
	  else
	    from=0;
	  display(seg->key, from, 0);
	}
      else if (class(seg->key) == _VClone)
	display (seg->key, look->key, TREE) ;
      else if (class(seg->key) == _VText)
	display (seg->parent, look->key, TREE) ;
      else
	display (seg->key, look->key, 0) ;
    }
}

/*.......... ..... .....
*/

static void pMapHighlight(KEYSET k, KEY key)
{
  int i, j;
  SEG *seg;
  void *dummy;
  KEYSET kSet=0; 
  PhysMap look = currentPhysMap("pMapHighlight");

  if (k!=0 && keySetExists(k) /*I probably ought to check this more thoroughly*/) kSet=k;
  if (kSet==0 && !keySetActive(&kSet, &dummy))
  {
     messout("First select a keySet window, thank you.");
     return;
  }
  seg=arrp(look->segs, 0, SEG);
  i=arrayMax(look->segs);
  while (i--)
  {
    seg->flag&= ~FLAG_HIGHLIGHT; /*clear this bit on all segs, except for selected objects:*/
    if (keySetFind(kSet, seg->key, &j)) seg->flag|=FLAG_HIGHLIGHT;
    seg++;
  }
  pMapDraw(look, (k==0? 0: key));
  return;
}

/*..... ..... ..... .....external service call routines
*/
extern void pMapSetHighlightKeySet(KEYSET k) /* called by gridDisp */
{
  pMapPendingKeySet=k;
  return;
}

/*..... ..... ..... .....
*/
extern BOOL pMapGetCMapKeySet(KEY contig, int *xMin, int *xMax, KEYSET clones)
{
  messout("**pMapGetCMapKeySet not implemented.");
  return FALSE;
}


/****************************************************************/
/******** RMD copied from pmaped.c to make stand-alone **********/

BOOL pMapLiftBuriedClones(Array segs, SEG *canonical, BOOL setDisplayFlag)
{				/* RD - BOOL is Neil crap */
  KEY child ;
  OBJ Parent, Child ;
  float x = SegMidPt(canonical) ;
  static BSMARK mark ;

  if ((Parent = bsCreate(canonical->key)))
    {
      if (bsGetKey (Parent, _Canonical_for, &child)) do
  /*at least one valid buried clone: get all buried clones for this seg*/
	{
	  int nBands, left, right ;
	  SEG *seg = arrayp(segs, arrayMax(segs), SEG) ; /* new seg */

	  mark = bsMark (Parent, mark) ;

	  seg->key = child ;
	  seg->parent = canonical->key ;
	  if (bsGetData (Parent, _bsRight, _Int, &left) &&
	      bsGetData (Parent, _bsRight, _Int, &right))
	    { seg->x0 = canonical->x0 + left ;
	      seg->x1 = canonical->x0 + right ;
	      Set(seg, FLAG_BURIED_IS_POSITIONED);
	    }
	  else
	    { if ((Child = bsCreate(child)))
		{ if (!bsGetData (Child, _Bands, _Int, &nBands) ||
		      !bsGetData (Child, _bsRight, _Int, &nBands))
		    nBands = 20 ;
		  bsDestroy (Child) ;
		}
	      else
		nBands = 20 ;
	      seg->x0 = x - nBands/2 ; 
	      seg->x1 = x + nBands/2 ;
	    }
	  Set(seg, FLAG_IS_BURIED_CLONE);
	  Set(seg, FLAG_FINGERPRINT);
	  if (setDisplayFlag) Set(seg, FLAG_DISPLAY);
	  bsGoto (Parent, mark) ;
	} while (bsGetKey (Parent, _bsDown, &child)) ;

      bsDestroy(Parent) ;
      return TRUE;
    }
  else
    return FALSE ;
}

static void displayAllBuriedClones(PhysMap look)
{
  int i ;
  SEG *pseg ;

  if (IsntSet(look, FLAG_BURIED_CLONES_ATTACHED))
    { graphMessage("Attaching buried clones: PLEASE WAIT...") ;
      
      for (i = 1 ; i < arrayMax (look->segs) ; ++i)
	{ pseg = arrp(look->segs, i, SEG) ;
	  if (IsntSet(pseg, FLAG_DISPLAY_BURIED_CLONES))
    /* if it is, then its buried clones are already loaded and checked */
	    pMapLiftBuriedClones(look->segs, pseg, FALSE) ;
	}

      arraySort(look->segs, pMapCompareSeg);
      Set(look, FLAG_BURIED_CLONES_ATTACHED);
      graphUnMessage() ;
    }
  return;
}

/****************** giface hook *****************/

void pMapDrawGIF (void)
{
  KEY key, from = 0 ;
  int x ;
  int z1, z2, winWidth ;
  char *word ;

  while (freestep ('-'))	/* options */
    if ((word = freeword()))
      {
	if (!strcmp (word, "clone") && freecheck("w"))
	  {
	    if (!lexword2key ((word = freeword()), &key, _VClone))
	      {
		freeOutf ("// gif pmap error: clone %s not known\n", word) ;
		return ;
	      }
	  }
	else if (!strcmp (word, "contig") && freecheck("wi"))
	  { 
	    if (!lexword2key ((word = freeword()), &key, _VContig))
	    {
	      freeOutf ("// gif pmap error: contig %s not known\n", word) ;
	      return ;
	    }
	    freeint (&x) ;
	    x |= 0x800000 ;
	    from = KEYMAKE (_VCalcul,x) ;
	  }
	else
	  goto usage ;
      }
    else
      goto usage ;

  if (freecheck ("w"))
    goto usage ;

  { extern BOOL displayReportGif ;
    displayReportGif = FALSE ;
    display (key, from, PMAP) ;	/* the primary display */
    displayReportGif = TRUE ;
  }

  z1 = Xoffset / Xmagn ;
  graphFitBounds (&winWidth, 0) ;
  z2 = z1 + winWidth / Xmagn ;

  if (class(key) == _VClone)
    { OBJ obj = bsUpdate (key) ;
      if (!obj)
	{ freeOutf ("// gif pmap error: clone %s has no data\n", freeprotect(name(key))) ;
	  return ;
	}
      if (!bsGetKey (obj, _pMap, &key))	/* replace key with contig */
	{ freeOutf ("// gif pmap error: clone %s has no pMap\n", freeprotect(name(key))) ;
	  bsDestroy (obj) ;
	  return ;
	}
      bsDestroy (obj) ;
    }

  freeOutf ("// PMAP %s %d %d\n", freeprotect(name(key)), z1, z2) ;

  return ;

usage:
  freeOut ("// gif pmap error: usage: PMAP [-clone <name>]|[-contig <name> <coord>]\n") ;
}
 
 
 
 
 
