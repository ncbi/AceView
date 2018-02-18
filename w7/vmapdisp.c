/*  File: vmapdisp.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Display of the genetic map
 * Exported functions:
       vMapDisplay
 * HISTORY:
 * Last edited: Jan  8 11:43 1999 (fw)
 * * May 29 22:48 1993 (cgc): split off vMapphys.c for contigs and clones
 * * Mar  2 19:58 1993 (rd): submenus for highlight and gMapdata
 * * Feb 27 23:06 1993 (rd): added Jean's addFather() stuff
 * * Nov 22 00:38 1992 (rd): switched to using map.h
 * * Jun  9 13:58 1992 (mieg): fixed pmap2vMap
 * * Mar  4 03:25 1992 (rd): added simplest form of chromosome banding
 * * Jan  6 13:11 1992 (rd): corrected green box scroll bug
 * * Dec 13 16:48 1991 (mieg): Added topMargin
 * * Dec 11 17:24 1991 (mieg): Added drag gene and private vMaps
 * * Nov 29 14:49 1991 (mieg): Added Locus display
 * Created: Fri Oct 18 20:25:36 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: vmapdisp.c,v 1.4 2008/04/07 22:18:19 mieg Exp $ */
/*
#define ARRAY_CHECK
#define MALLOC_CHECK
*/

#include "acedb.h"
#include "bump.h"

#include "display.h"
#include "pick.h"
#include "session.h"
#include "systags.h"
#include "tags.h"
#include "classes.h"
#include "sysclass.h"
#include "bs.h"
#include "lex.h"
#include "a.h"
#include "vmap_.h"

/************************************************************/

static void vMapDestroy (void) ;
static void vMapRecalculate(void) ;
static void vMapPick (int box, double x, double y) ;
static void vMapKbd (int k) ;
static Array vMapConvert (KEY chrom) ;
static void vMapFollow (VerticalMap look, float x, float y) ;
static void vMapClear(void) ;
static void vMapHighlight(void) ;
static void vMapUnHighlight(void) ;
static void vMapHideHighlit(void) ;
static void vMapExportHighlit(void) ;
static void vMapSaveMap (void) ;
static void vMapSetMap (VerticalMap look) ;
static void vMapGotoCmap (void) ;
static void hideHeaderToggle (void) ;
static void vMapSearchProblems(VerticalMap look) ;

static void vMapDrawVoid (void) ;
static void vMapFlip (void) ;
static void vMapDoFlip (VerticalMap look) ;


/* the MapDrawColFunc's in this file */
static void vMapDrawMiniChromBands (LOOK genericLook, float *offset) ;
static void vMapDrawMainMarkers (LOOK genericLook, float *offset) ;
static void vMapDrawAnyInterval (LOOK genericLook, float *offset) ;
static void vMapDrawAnyIntervalNames (LOOK genericLook, float *offset) ;
static void vMapDrawChromBands (LOOK genericLook, float *offset) ;
static void vMapDrawScale (LOOK genericLook, float *offset) ;
static void vMapDrawAnyLocus (LOOK genericLook, float *offset) ;
static void vMapDrawOrderedGenes (LOOK genericLook, float *offset) ;


static MENUOPT vMapMenu[] = {
  {graphDestroy, "Quit"},
  {help, "Help"},
  {graphPrint, "Print Screen"},
  {mapPrint,"Print Whole Map"},
  {displayPreserve, "Preserve"},
  {vMapRecalculate, "Recalculate"},
  {vMapFlip, "Vertical Flip"},
  {vMapSaveMap, "Save Map"},
  {hideHeaderToggle, "Hide Header"},
  {vMapGotoCmap, "Physical Chromo Map"},
  {mapColControl, "Columns"},
  {0, 0}
} ;

static BOOL showMarginal = FALSE ;
static BOOL setCursorDrag = FALSE ; /* cursor operations a la richard */
static Associator anyIntAss = 0 ;
static Associator allMapsAss = 0 ;

/************************************************************/

magic_t GRAPH2VerticalMap_ASSOC = "VerticalMap";

static magic_t VerticalMap_MAGIC = "VerticalMap";

/************************************************************/

     /* Recalculate all vMaps */
void vMapMakeAll(void)
{
  KEY chromosome = 0 ;	/* 0 primes lexNext() */
  Array segs ;

  if (!isWriteAccess ())
    { messout ("Sorry, you do not have Write Access") ;
      return ;
    }


    while (lexNext (_VMap,&chromosome))
      { if ((segs = vMapConvert (chromosome)))
	  arrayDestroy (segs) ;
	if (messIsInterruptCalled ())
	  break ;
      }
}

/**********************************/

BOOL vMapGetPos(KEY from, KEY *chromo, float *xp)
{ 
  OBJ obj ;
  float x1, x2 ;
  BOOL result = FALSE ; 

  if (!from || !(obj = bsCreate (from)))
    return FALSE ;

  if (*chromo && bsFindKey (obj, _Map, *chromo))
    result = TRUE ;

  if (result || bsGetKey (obj, _Map, chromo))
    if (bsPushObj(obj))
      {
	if (bsGetData (obj, _Position, _Float, xp) ||
	    bsGetData (obj, _Multi_Position, _Float, xp))
	  result = TRUE ;
	else if ( bsGetData (obj, _Left, _Float, &x1) &&
		 bsGetData (obj, _Right, _Float, &x2))
	  { *xp = 0.5*(x1+x2) ;
	    result = TRUE ;
	  }
      }
  bsDestroy (obj) ;
  return result ;
}
  
/**********************************/

VerticalMap currentVerticalMap (char *caller)
{
  VerticalMap vmap;

  if (!graphAssFind (&GRAPH2VerticalMap_ASSOC,&vmap))
    messcrash("%s() could not find VerticalMap on graph", caller);
  if (!vmap)
    messcrash("%s() received NULL VerticalMap pointer", caller);
  if (vmap->magic != &VerticalMap_MAGIC)
    messcrash("%s() received non-magic VerticalMap pointer", caller);
  
  return vmap;
} /* currentVerticalMap */

/**********************************/

BOOL vMapDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  float centre, mag ; 
  VerticalMap  oldlook = 0, look ;
  OBJ   obj ;
  KEY   vMapkey ;
  Array segs ;
  BOOL doflip = FALSE ;

  if (key && class(key) != _VvMap && class(key) != _VMap)
    { KEY tmp = from ;
      from = key ;
      if (class(tmp) == _VvMap || class(tmp) == _VMap)
	key = tmp ;
      else
	key = 0 ;
    }
		/* centre on parent object if poss */
  centre = mag = 99999 ;
  if (from)
    {
      if (class(from) == _VCalcul)
	{ centre = (float)(KEYKEY(from))/1000.0 - 1000 ;
	  from = 0 ;
	}
      else if (!vMapGetPos (from, &key, &centre))
	{ display (from, 0, TREE) ;
	  return FALSE ;
	}
    }

  if (!key)
    return FALSE ;

  vMapkey = 0 ; segs = 0 ;
  if (class(key) == _VvMap)
    {
      vMapkey = key ;
      lexReClass (vMapkey, &key, _VMap) ;
    }
  else if (class(key) == _VMap)
      lexReClass (key, &vMapkey, _VvMap) ; /* to arrayGet if possible */
  

  if (
      (
       ( allMapsAss && assFind(allMapsAss, (char*)(&allMapsAss) + vMapkey, &oldlook) &&
	graphActivate(oldlook->graph) )
       ||
       ( isOldGraph && graphAssFind (&GRAPH2VerticalMap_ASSOC, &oldlook)    )
	) 
      &&
      oldlook->magic == &VerticalMap_MAGIC &&
      oldlook->key == vMapkey )
    { if (centre != 99999)
	{ if (oldlook->map->flip)
	    centre = -centre ;
	  oldlook->map->centre = centre ;
	  if (setCursorDrag)
	    mapCursorSet (oldlook->map, centre) ;
	}
      vMapDraw (oldlook, from) ;
      return TRUE ;
    }

  if ((segs = arrayGet(vMapkey,SEG, segFormat)))
    { int i = arrayMax(segs) ;
      SEG* seg = arrp(segs, 0, SEG) - 1 ;
      unsigned int mask = FLAG_RELATED | FLAG_STRESSED | FLAG_ANTI_RELATED | 
	FLAG_ANTI_STRESSED |  FLAG_HAVE_DATA | FLAG_HIDE |  FLAG_HIGHLIGHT | 
	FLAG_MOVED | FLAG_PROBLEM ;
      while (seg++, i--)
	seg->flag &= ~mask ;
    }
  else if ((segs = vMapConvert (key))) ;
  else
    { display (key,from,TREE) ;
      return FALSE ;
    }
  lexReClass(key, &vMapkey, _VvMap) ; /* if segs newly made */
  if (oldlook)
    assRemove(allMapsAss, (char*)(&allMapsAss) + oldlook->key) ;

  if (isOldGraph)
    { graphRetitle (name (key)) ;
      vMapDestroy () ; /* does a messfree(look) */
    }
  else 
    {
      if (!displayCreate (VMAP))
	return FALSE ;

      if (!allMapsAss)
	allMapsAss = assCreate() ;
      graphRetitle (name(key)) ;
      graphRegister (RESIZE,vMapDrawVoid) ;
      graphRegister (DESTROY, vMapDestroy) ;
      graphRegister (KEYBOARD, vMapKbd) ;
    }

  look = (VerticalMap) messalloc (sizeof (struct VerticalMapStruct)) ;  
  look->magic = &VerticalMap_MAGIC;
  look->flag = 0 ;
  look->key = vMapkey ;
  look->boxIndex = arrayCreate (64,int) ;
  look->activeBox = 0 ;
  look->segs = segs ;
  look->lastTrueSeg = arrayMax(look->segs) ;
  look->errorScale = 10 ;
  if ((obj = bsCreate(key)))
    { bsGetData(obj, _Error_scale, _Float, &look->errorScale) ;
      doflip = bsFindTag (obj, _Flipped) ;
      /* rather, i make an ABOUT button => Map in TREE mode 
	 bsGetKey(obj, _Title, &look->titleKey) ;
	 look->remarkStack = stackCreate(50) ;
	 if (bsGetData(obj, _Remark, _Text, &cp))
	 do { pushText(look->remarkStack, cp) ;
	 } while (bsGetData(obj, _bsDown, _Text, &cp)) ;
	 */
      
      bsDestroy(obj) ;
    }
  
  look->graph = graphActive() ;

  look->map = mapCreate (vMapDrawVoid) ;
  mapInsertCol (look->map, 1, TRUE, "Mini Bands", vMapDrawMiniChromBands) ;
  mapInsertCol (look->map, 2, TRUE, "Marker Genes", vMapDrawMainMarkers) ;
  mapInsertCol (look->map, 3, TRUE, "Locator", mapShowLocator) ;
  mapInsertCol (look->map, 4, TRUE, "AnyInterval", vMapDrawAnyInterval) ;
  mapInsertCol (look->map, 5, TRUE, "AnyIntervalNames", vMapDrawAnyIntervalNames) ;
  mapInsertCol (look->map, 6, TRUE, "Chromosome Bands", vMapDrawChromBands) ;
  mapInsertCol (look->map, 7, TRUE, "Contigs", vMapDrawContigs) ;
  mapInsertCol (look->map, 8, FALSE, "Reversed physical", vMapReversedPhysical) ;
  mapInsertCol (look->map, 9, FALSE, "Physical genes", vMapDrawPhysGenes) ;
  mapInsertCol (look->map, 10, TRUE, "Scale", vMapDrawScale) ;
  mapInsertCol (look->map, 11, TRUE, "2 point", vMapDraw2pt) ;
  mapInsertCol (look->map, 12, TRUE, "Multipoint", vMapDrawMultiPt) ;
  mapInsertCol (look->map, 13, FALSE, "Likelihood", vMapDrawDbn) ;
  mapInsertCol (look->map, 14, FALSE, "Ordered genes", vMapDrawOrderedGenes) ;
  mapInsertCol (look->map, 15, TRUE, "AnyLocus", vMapDrawAnyLocus) ;

  vMapSetMap (look) ;		/* reset look->map params etc. */

  graphAssociate (&GRAPH2VerticalMap_ASSOC, look) ; /* attach look to graph */
  graphAssociate (&MAP2LOOK_ASSOC, look) ; /* attach map to look */
  mapAttachToGraph (look->map) ; /* attach the map itself to the graph */

  assInsert (allMapsAss, (char*)(&allMapsAss) + vMapkey, look) ;

  if (mag != 99999)
    look->map->mag = mag ;
  if (centre != 99999)
    look->map->centre = centre ;
  if (doflip)
    vMapDoFlip(look) ;
    
  if (setCursorDrag)
    mapCursorCreate (look->map, 0.01, look->map->centre) ;
  
  vMapDraw (look, from) ;

  return TRUE ;
} /* vMapDisplay */

/************************************************************/
/***************** Registered routines *********************/

static void vMapDestroy (void)
{
  int i ;
  VerticalMap look = currentVerticalMap("vMapDestroy") ;

  if (isWriteAccess())
    for (i = 0 ; i < arrayMax(look->segs) ; ++i)
      if (arrp(look->segs,i,SEG)->flag & FLAG_MOVED)
	{ if (messQuery ("Do you want to save your changes?"))
	    vMapSaveMap () ;
	  break ;
	}

  if (look->data)
    vMapDataDestroy (look) ;

  arrayDestroy (look->segs) ;
  arrayDestroy (look->boxIndex) ;
  stackDestroy (look->remarkStack) ;

  mapDestroy (look->map) ;	/* also detach map from graph */

  graphAssRemove (&GRAPH2VerticalMap_ASSOC) ; /* detach look from graph */
  graphAssRemove (&MAP2LOOK_ASSOC) ; /* detach map from look */

  assRemove(allMapsAss, (char*)(&allMapsAss) + look->key) ;
  look->magic = 0 ;
  messfree (look) ;

  return;
} /* vMapDestroy */

/*************************************************************/

static void vMapRecalculate(void)
{
  int i ;
  KEY chrom, curr ;
  VerticalMap look = currentVerticalMap("vMapRecalculate") ;
  
  if (lexReClass(look->key,&chrom,_VMap))
    { if (look->activeBox && 
	  (i = arr(look->boxIndex, look->activeBox, int)))
	curr = arrp(look->segs, i, SEG)->key ;
      else 
	curr = 0 ;
      arrayDestroy (look->segs) ;
      if ((look->segs = vMapConvert (chrom)))
	{ look->lastTrueSeg = arrayMax(look->segs) ;
	  vMapDraw (look, curr) ;
	}
      else
	{ messout ("Sorry, I have to abandon this vMap") ;
	  graphDestroy () ;
	}
    }
}

/*************************************/

static void vMapPick (int box,  double x , double y) 
{
  VerticalMap look = currentVerticalMap("vMapPick") ;

  if (!box)
    return ;
  if (setCursorDrag && box == look->map->cursor.pickBox)
    graphBoxDrag (look->map->cursor.pickBox, mapCursorDrag) ;
  else if (box == look->map->thumb.box)
    graphBoxDrag (look->map->thumb.box, mapThumbDrag) ;
  else if (box < arrayMax (look->boxIndex) && arr(look->boxIndex,box,int)) /* a SEG */
    {
      if (box == look->activeBox) /* a second hit - follow it */
	vMapFollow (look, (float)x, (float)y) ;
      else
	vMapSelect (look, box) ;
    }

  return;
} /* vMapPick */

/*****************/

static void vMapKbd (int k)
{ 
  float unit, subunit ;
  int i, j, box ;
  SEG *seg ;
  KEY mainKey ;
  VerticalMap look = currentVerticalMap("vMapKbd") ;

  if (FALSE && !look->activeBox && setCursorDrag && look->map->cursor.box)
    { mapFindScaleUnit (look->map, &unit, &subunit) ;

      switch (k)
	{
	case UP_KEY :
	  look->map->cursor.val -= 5*subunit ;
	  break ;
	case DOWN_KEY :
	  look->map->cursor.val += 5*subunit ;
	  break ;
	default:
	  return ;
	}
      mapCursorShift (look->map) ;
      return ;
    }
 
  box = look->activeBox ;

  if (!box ||
      !(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)))
    return ;
	
  mainKey = seg->key ;

  switch (k)
    {
    case UP_KEY :
      if (box < 2)
	return ;
      j = arr(look->boxIndex, box - 1, int) ;
      if (j >= 0 && j < arrayMax(look->segs) &&
	  (seg = arrp(look->segs,j,SEG)) )
	box-- ;
      else
	return ;
      break ;
    case DOWN_KEY :
      if (box >= arrayMax(look->boxIndex) - 1)
	return ;
      j = arr(look->boxIndex, box + 1, int) ;
      if (j >= 0 && j < arrayMax(look->segs) &&
	  (seg = arrp(look->segs,j,SEG)) )
	box++ ;
      else
	return ;
      break ;
    case LEFT_KEY:
      for (box-- ; box > 0 ; box--)
	{ j = arr(look->boxIndex, box, int) ;
	  if (j >= 0 && j < arrayMax(look->segs) &&
	      (seg = arrp(look->segs,j,SEG)) &&
	      seg->key == mainKey )
	    goto foundIt ;
	}
      return ;
    case RIGHT_KEY:
      for (box++ ; box < arrayMax(look->boxIndex) ; box++)
	{ j = arr(look->boxIndex, box, int) ;
 	  if (j >= 0 && j < arrayMax(look->segs) &&
	      (seg = arrp(look->segs,j,SEG)) &&
	      seg->key == mainKey )
	    goto foundIt ;
	}
      return ;
    default:
      return ;
    }

  foundIt:
  if (box < 1 || box >= arrayMax(look->boxIndex))
    for (box = 1 ; box < arrayMax(look->boxIndex) ; box++)
      { j = arr(look->boxIndex, box, int) ;
 	if (j >= 0 && j < arrayMax(look->segs) &&
	    (seg = arrp(look->segs,j,SEG)))
	  break ;
      }
  if (box > 0 && 
      box < arrayMax(look->boxIndex) &&
      (j = arr(look->boxIndex, box, int) ) &&
      j >= 0 && j < arrayMax(look->segs) &&
      (seg = arrp(look->segs,j,SEG)))
    vMapSelect(look, box) ;
}

/**********************************/

static void vMapGotoCmap (void)
{
  KEY chrom ;
  VerticalMap look = currentVerticalMap("vMapGotoCmap") ;

  if (lexReClass (look->key,&chrom,_VMap))
    display (chrom, 0, CMAP) ;
}

/*************************************************************/
/******** vMap show intersect with active keyset *************/
 
static void vMapClear (void)
{
  int i, mask = FLAG_STRESSED | FLAG_RELATED | FLAG_ANTI_RELATED |
    FLAG_ANTI_STRESSED | FLAG_HIDE | FLAG_HIGHLIGHT ;
  VerticalMap look = currentVerticalMap("vMapClear") ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    arrp(look->segs,i,SEG)->flag &= ~mask ;
  look->activeBox = 0 ;
  vMapDraw (look, 0) ;
}

/**********/
 
static void vMapHighlightAll (void)
{
  int i ;
  VerticalMap look = currentVerticalMap("vMapHighlightAll") ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    arrp(look->segs,i,SEG)->flag |= FLAG_HIGHLIGHT ;
  look->activeBox = 0 ;
  vMapDraw (look, 0) ;
}

/**********/

static void vMapHighlight (void)
{
  int i, j ;
  SEG *seg ;
  void *dummy ;
  KEYSET kSet = 0 ;
  VerticalMap look = currentVerticalMap("vMapHighlight") ;

  if (!keySetActive(&kSet, &dummy))
    { messout("First select a keySet window, thank you.") ;
      return ;
    }
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (keySetFind (kSet, seg->key, &j))
	seg->flag |= FLAG_HIGHLIGHT ;
    }
  for (i = 0 ; i < arrayMax(look->boxIndex) ; ++i)
    if ((j = arr(look->boxIndex, i, int)))
      { seg = arrp(look->segs, j, SEG) ;
	if (seg && ( seg->flag & FLAG_HIGHLIGHT))
	  graphBoxDraw (i, BLACK, MAGENTA) ;
      }
}

/**********/

static void vMapUnHighlight (void)
{
  int i, j;
  SEG *seg ;
  void *dummy ;
  KEYSET kSet = 0 ;
  VerticalMap look = currentVerticalMap("vMapUnHighlight") ;

  if (!keySetActive(&kSet, &dummy))
    { messout("First select a keySet window, thank you.") ;
      return ;
    }
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (keySetFind (kSet, seg->key, &j))
	seg->flag &= ~FLAG_HIGHLIGHT ;
    }
  vMapDraw (look, 0) ;
}

/**********/

static void vMapExportHighlit (void)
{
  static KEYSET kset = 0 ;
  int i, j=0 ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("vMapUnHighlight") ;

  kset = keySetReCreate (kset) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->flag & FLAG_HIGHLIGHT)
	keySet(kset,j++) = seg->key ;
    }
  keySetSort (kset) ;

  displayCreate (DtKeySet) ;	/* new window */
  keySetShow (kset, 0) ;
}

/**********/

static void vMapHideHighlit (void)
{
  int i ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("vMapHide") ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->flag & FLAG_HIGHLIGHT)
	seg->flag |= FLAG_HIDE ;
    }
  vMapDraw (look, 0) ;
}

static void vMapKeepHighLit (void)
{
  int i ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("vMapHide") ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->flag & FLAG_HIGHLIGHT)
	seg->flag &= ~(FLAG_HIDE | FLAG_HIGHLIGHT) ;
      else
	seg->flag |= FLAG_HIDE ;
    }
  vMapDraw (look, 0) ;
}

static void vMapHighlightProblems (void)
{
  int i, iseg ; static void *v ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("vMapHide") ;

  if (!graphAssFind(&v, v))
    vMapSearchProblems(look) ;
  graphAssociate(&v, &v) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->flag & FLAG_PROBLEM)
	seg->flag |= FLAG_HIGHLIGHT ;
      else
	seg->flag &= ~FLAG_HIGHLIGHT ;
    }
  for (i = 0 ; i < arrayMax(look->boxIndex) ; ++i)
    if ((iseg = arr(look->boxIndex, i, int)))
      { seg = arrp(look->segs, iseg, SEG) ;
	if (seg && (seg->flag & FLAG_HIGHLIGHT) )
	  graphBoxDraw (i, BLACK, MAGENTA) ;
      }
}

static void vMapUnHide (void)
{
  int i ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("vMapHide") ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      seg->flag &= ~FLAG_HIDE ;
    }
  vMapDraw (look, 0) ;
}

static void vMapReColor (VerticalMap look)
{ SEG *seg ;
  int i, iseg ;
  unsigned int mask = 
    FLAG_RELATED | FLAG_STRESSED |
      FLAG_ANTI_RELATED | FLAG_ANTI_STRESSED ;
  /* mask to be able to pick after highlight is on ! */

  for (i = 0 ; i < arrayMax(look->boxIndex) ; ++i)
    if ((iseg = arr(look->boxIndex, i, int)))
      { seg = arrp(look->segs, iseg, SEG) ;
	if (!seg)
	  continue ;
	if ((seg->flag & FLAG_HIGHLIGHT) &&
	    !(seg->flag & mask))
	  graphBoxDraw (i, BLACK, MAGENTA) ;
	else
	  {
	    if (class (seg->key) == _V2_point_data ||
		class (seg->key) == _VMulti_pt_data)
	      {
		if (seg->flag  & FLAG_RELATED)
		  graphBoxDraw (i, BLACK, WHITE /*YELLOW*/) ;
	      }
	    else if (class (seg->key) == _VClone ||
		     class (seg->key) == _VContig)
	      graphBoxDraw (i, BLACK, YELLOW) ;
	    else
	      {
		if (seg->flag & FLAG_STRESSED)
		  graphBoxDraw (i, WHITE, RED) ;
		else if (seg->flag & FLAG_ANTI_RELATED)
		  graphBoxDraw (i, BLACK, LIGHTBLUE) ;
		else if (seg->flag & FLAG_ANTI_STRESSED)
		  graphBoxDraw (i, WHITE, RED) ;
		else if (seg->flag & FLAG_RELATED)
		  graphBoxDraw (i, BLACK, GREEN) ;
		else if (seg->flag & FLAG_CLONED)
		  graphBoxDraw (i, BLACK, YELLOW) ;
		else
		  graphBoxDraw (i, BLACK, WHITE) ;
	      }
	  }
      }
}

/***************************************************************/
/****** magnification control etc, uses special tag keys *******/

static void vMapSetMap (VerticalMap look) 
{ 
  int i ;
  SEG *seg, *extentSeg = 0, *centreSeg = 0 ;
  MAP map = look->map ;

  graphFitBounds (&mapGraphWidth,&mapGraphHeight) ;
  halfGraphHeight = 0.5 * (mapGraphHeight - topMargin) ; 

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs, i, SEG) ;
      if (seg->key == _Extent)
	extentSeg = seg ;
      if (seg->key == _Centre)
	centreSeg = seg ;
    }

  if (extentSeg)
    { map->min = extentSeg->x ;
      map->max = extentSeg->x + extentSeg->dx ;
    }
  else
    { map->min = 100000000.0 ;
      map->max = -100000000.0 ;
      for (i = 1 ; i < arrayMax(look->segs) ; ++i) /* 0th seg fake */
	{ seg = arrp(look->segs, i, SEG) ;
	  if (seg->x < map->min)
	    map->min = seg->x ;
	  if (seg->x > map->max)
	    map->max = seg->x ;
	}
    }

  map->centre = 0 ; 
  map->mag = 10 ;
  if (centreSeg)
    { map->centre = centreSeg->x ;
      if (centreSeg->dx)
	map->mag = (2*halfGraphHeight-5) / centreSeg->dx ;
    }
}

/******************************************/

static void highlightSeg (int box)
{
  int i ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("itemName") ;

  if (!(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)))
    return ;
  seg->flag ^= FLAG_HIGHLIGHT ;
  if (seg->flag & FLAG_HIGHLIGHT)
    graphBoxDraw (box, BLACK, MAGENTA) ;
  else if (class(seg->key) == _VClone ||
	   class(seg->key) == _VContig ||
	   (class(seg->key) == _VLocus && (seg->flag & FLAG_CLONED)))
    graphBoxDraw (box, BLACK, YELLOW) ;
  else
    graphBoxDraw (box, BLACK, WHITE) ;
}

static void nameSeg (int box)
{
  int i ;
  SEG *seg ;
  float x1, x2, y1, y2 ;
  VerticalMap look = currentVerticalMap("nameSeg") ;

  if (!(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)))
    return ;
  graphBoxDim (box, &x1, &y1, &x2, &y2) ;
  graphTextFormat (BOLD) ;
  graphText (name(seg->key), x2, y1) ;
  graphTextFormat (PLAIN_FORMAT) ;
  graphRedraw () ;
}

/************ cursor stuff ***********/

static void cursorToTop (int box)
{
  int i ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("cursorToTop") ;

  if (!(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)))
    return ;
  mapCursorSet (look->map, seg->x) ;
}

static void cursorToBottom (int box)
{
  int i ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("cursorToBottom") ;

  if (!(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)))
    return ;
  mapCursorSet (look->map, seg->x + seg->dx) ;
}

static void cursorToPos (int box)
{
  int i ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("cursorToPos") ;

  if (!(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)))
    return ;
  mapCursorSet (look->map, seg->x) ;
}

static void topToCursor (int box)
{
  int i ;
  SEG *seg ;
  float a, b ;
  VerticalMap look = currentVerticalMap("topToCursor") ;

  if (!(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)))
    return ;
  a = mapCursorPos(look->map) ;
  mapCursorSet (look->map, seg->x) ;
  b = seg->x + seg->dx ;
  if (b > a) 
    { seg->x = a ; 
      seg->dx = b-a ; 
    }
  else 
    { seg->x = b ; 
      seg->dx = a-b ; 
    }
  seg->flag |= FLAG_MOVED ;
  vMapDraw (look, 0) ;
}

static void bottomToCursor (int box)
{
  int i ;
  SEG *seg ;
  float a, b ;
  VerticalMap look = currentVerticalMap("bottomToCursor") ;

  if (!(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)))
    return ;
  a = seg->x ;
  b = mapCursorPos(look->map) ;
  mapCursorSet (look->map, seg->x + seg->dx) ;
  if (b > a)
    seg->dx = b-a ;
  else
    { seg->x = b ;
      seg->dx = a-b ;
    }
  seg->flag |= FLAG_MOVED ;
  vMapDraw (look, 0) ;
}

static void posToCursor (int box)
{
  int i ;
  SEG *seg ;
  float a ;
  VerticalMap look = currentVerticalMap("posToCursor") ;

  if (!(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)))
    return ;
  a = mapCursorPos(look->map) ;
  mapCursorSet (look->map, seg->x) ;
  seg->x = a ;
  seg->flag |= FLAG_MOVED ;
  vMapDraw (look, 0) ;
}

static void vMapShowMarginal(void)
{
  showMarginal = !showMarginal ;

}

static void vMapSetCursorDrag(void)
{
  VerticalMap look = currentVerticalMap("vMapSetCursorDrag") ;
  
  setCursorDrag = ! setCursorDrag ;
  if (!look->map->cursor.unit)
    look->map->cursor.unit = .01 ;
  vMapDraw (look,0) ;
}

static void vMapAnalysis (void)
{ return ;
}

/*****************************************************************/
/****************** Drawing routines *****************************/

static MENUOPT buttonOpts[] = {
  {mapWhole, "Whole"}, 
  {mapZoomIn, "Zoom In"},
  {mapZoomOut, "Zoom Out"},
  {vMapHighlight, "Highlight..."},
  {vMapShowData, "Map Data"},
  {vMapDragButton, "Drag..."},
  {vMapAnalysis, "Analysis..."},
  {0, 0 }} ;

static MENUOPT analysisMenu[] = {
  {vMapClear, "Clear"},
  {vMapSetCursorDrag, "Cursor dragging on/off"},
  {vMapHighlightProblems, "Highlight Problems"},
  {vMapShowMarginal, "Marginal errors in pink"},
  {0, 0 }} ;


static MENUOPT clearMenu[] = {
  {vMapClear, "Clear"},
  {vMapHighlightAll, "Highlight all"},
  {vMapHighlight, "Highlight selected keyset"},
  {vMapUnHighlight, "Unhighlight selected keyset"},
  {vMapKeepHighLit, "Keep Only Highlighted objects"},
  {vMapHideHighlit, "Hide highlit items"},
  {vMapUnHide, "Unhide hidden items"},
  {vMapExportHighlit, "Export highlit keyset"},
  {0, 0} } ;

static MENUOPT dragMenu[] = {
  {vMapDragButton, "Drag"},
  {vMapDragUndo, "Undo"},
  {0, 0 }} ;

static BOOL isDrawnOrdered ;

/*******************/

static void aboutButton (void)
{
  KEY key ;
  VerticalMap look = currentVerticalMap("aboutButton") ;

  lexReClass (look->key, &key, _VMap) ;
  display (key, 0, TREE) ;
}

static void hideHeaderToggle (void)
{
  VerticalMap look = currentVerticalMap("hideHeaderToggle") ;

  look->flag ^= FLAG_HIDE_HEADER ;
  vMapDraw (look, 0) ;
}

static void vMapDrawVoid (void) 
{ 
  VerticalMap look = currentVerticalMap("vMapDrawVoid") ; 

  vMapDraw (look, 0) ; 
}

void vMapDraw (VerticalMap look, KEY curr)
{
  int i, iseg, trueMax ;
  SEG *seg ;

  if (!curr && look->activeBox && 
      look->activeBox > 0 &&
       look->activeBox < arrayMax(look->boxIndex) &&
      (i = arr(look->boxIndex,look->activeBox,int)) &&
      i >=0 && i < arrayMax (look->segs))
    curr = arrp(look->segs, i, SEG)->key ;

  graphClear () ;

  graphFitBounds (&mapGraphWidth,&mapGraphHeight) ;
  topMargin = 3 ;
  halfGraphHeight = 0.5 * (mapGraphHeight - topMargin) ; 

  if (mapGraphHeight < 10 + topMargin)
    { messout ("Sorry, this window is too small for a genetic vMap") ;
      return ;
    }

  trueMax = arrayMax(look->segs) ;
  arrayMax(look->segs) = look->lastTrueSeg ;
  arraySort (look->segs, vMapOrder) ;
  arrayMax(look->segs) = trueMax ;

  look->boxIndex = arrayReCreate (look->boxIndex,50,int) ;
  look->activeBox = 0 ;

  isDrawnOrdered = FALSE ;
  mapDrawColumns (look->map) ;

  if (look->flag & FLAG_HIDE_HEADER)
    look->messageBox = 0 ;
  else
    { *look->messageText = 0 ;
      look->messageBox = graphBoxStart () ;
      graphTextPtr (look->messageText, 1, 2, 127) ;
      graphBoxEnd () ;
      graphBoxDraw (look->messageBox, BLACK, CYAN) ;
      
      graphButton ("About", aboutButton , mapGraphWidth - 6.4, 3.2) ;
/*
      if (iskey (look->titleKey))
	{ graphTextFormat (BOLD) ;
	  graphText (name (look->titleKey),3, 3) ;
	  graphTextFormat (PLAIN_FORMAT) ;
	}
      
      stackCursor(look->remarkStack, 0) ;
      i = 3 ;
      while(cp = stackNextText(look->remarkStack))
	graphText(cp, mapGraphWidth - 20, i++) ;
*/
      i = graphButtons (buttonOpts, 5, 0.5, mapGraphWidth) ; 
      mapColMenu (i);			/* "Whole" */
					/* "Zoom In"  NOMENU */
					/* "Zoom Out" NOMENU */
      graphBoxMenu (i+3, clearMenu) ;	/* "Highlight..." */
      graphBoxMenu (i+4, vMapDataMenu); /* "Map Data" */
      graphBoxMenu (i+5, dragMenu);     /* "Drag..." */
      graphBoxMenu (i+6, analysisMenu); /* "Analysis..." */
    }

         /* Re pick */
  if (curr)
    for (i = 0 ; i < arrayMax(look->boxIndex) ; ++i)
      if ((iseg = arr(look->boxIndex, i, int)))
	{ seg = arrp(look->segs, iseg, SEG) ;
	  if (seg && seg->key == curr) 
	    { vMapSelect(look, i);	/* after drawing messageBox! */
	      break ;
	    }
	}
         /* Re color */

  vMapReColor (look) ;
  graphRegister (PICK, vMapPick) ; /* redo because of ColControl */
  graphMenu (vMapMenu) ;
  
  graphRedraw () ;
}

/**********************/
 
static void vMapDoFlip (VerticalMap look)
{
  int i; float y ;
  SEG   *seg ;

  look->map->flip = !look->map->flip ;
  i = arrayMax(look->segs) ;
  seg = arrp(look->segs,0,SEG) - 1 ;
  while (seg++, i--)    
    if ( seg->key == _Centre ||
	(seg->flag & (FLAG_ANY_LOCUS | FLAG_MULTIPLE_LOCUS)))
	seg->x = - seg->x ;
    else if (seg->key == _Extent ||
	     (seg->flag & FLAG_ANY_INTERVAL))
	seg->x = - seg->x - seg->dx ;
  arraySort(look->segs, vMapOrder) ;
  y = look->map->min ;
  look->map->min =  - look->map->max ;
  look->map->max = - y ;
  look->map->centre = - look->map->centre ;
  if (setCursorDrag && look->map->cursor.box)
    { look->map->cursor.val = - look->map->cursor.val ;
      mapCursorShift (look->map) ;
    }
}
 
static void vMapFlip (void)
{
  VerticalMap look = currentVerticalMap("vMapFlip") ; 

  vMapDoFlip(look) ;
  vMapDraw (look, 0) ;
}
/**********************/

static MENUOPT vGeneMenu[] = {
  {(VoidRoutine)cursorToPos, "Set cursor"},
  {(VoidRoutine)cursorToTop, "Cursor to upper limit"},
  {(VoidRoutine)cursorToBottom, "Cursor to lower limit"},
  {(VoidRoutine)posToCursor, "Move to Cursor"},
  {(VoidRoutine)highlightSeg, "Highlight On/Off"},
  {(VoidRoutine)vMapShowDataFromMenu, "Show data"},
  {(VoidRoutine)vMapMakeDbn, "Likelihood dbn"},
  {0, 0 }} ;

static void vMapDrawOrderedGenes (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  BUMP  bump ;
  float y ;
  int   i, ibox, n, x ;
  SEG   *seg ;
  int   hideHeader = (look->flag & FLAG_HIDE_HEADER) ;

  *offset += 1 ;		/* bump bug */
  x = (mapGraphWidth - *offset)/2 ;
  if (x < 12) x = 12 ;
  bump = bumpCreate (x, 0) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->flag & FLAG_HIDE)
	continue ;
      y = MAP2GRAPH(look->map,seg->x) ;
      if (y > topMargin && y < mapGraphHeight-1  &&
	  class(seg->key) == _VLocus && 
	  (seg->flag & FLAG_WELL_ORDERED) &&
	  !(seg->flag & FLAG_PHYS_GENE))
	{ array(look->boxIndex,ibox=graphBoxStart(),int) = i ;
	  x = 0 ;
	  /* was:
	  bumpItem (bump, strlen(name(seg->key))+1, 1, &x, &y) ;
	  graphText (name(seg->key), *offset+x, y-0.5) ;
	  if (seg->flag & FLAG_CLONED)
	    if (hideHeader)
	      { graphLine (*offset+x, y+0.3, 
			   *offset+x+0.8*strlen(name(seg->key)), y+0.3) ;
		graphBoxDraw (ibox, BLACK, WHITE) ;
	      }
	    else
	      graphBoxDraw (ibox, BLACK, YELLOW) ;
	  else
	    graphBoxDraw (ibox, BLACK, WHITE) ;
	    */
	  if ((n = bumpText (bump, name(seg->key), &x, &y, 1, TRUE)))
	    { char *ccp = name(seg->key), *ccq = ccp +n, ccc = *ccq ;
	      *ccq = 0 ; graphText (ccp, *offset+x, y-0.5) ; *ccq = ccc ;
	      if (seg->flag & FLAG_CLONED)
		{
		  if (hideHeader)
		    { graphLine (*offset+x, y+0.3, 
				 *offset+x+0.8*strlen(name(seg->key)), y+0.3) ;
		      graphBoxDraw (ibox, BLACK, WHITE) ;
		    }
		  else
		    graphBoxDraw (ibox, BLACK, YELLOW) ;
		}
	      else
		graphBoxDraw (ibox, BLACK, WHITE) ;
	    }
	  graphBoxEnd() ;
	  if (setCursorDrag)
	    graphBoxMenu (ibox, vGeneMenu) ;
	}
    }

  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;
  isDrawnOrdered = TRUE ;
}


static void vMapDrawAnyLocus (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  float y , errorScale = look->errorScale ;
  int   i, ibox, n, x ;
  SEG   *seg ;
  BUMP  bump ;

  *offset += 1 ;		/* bump bug */
  if (*offset > mapGraphWidth-6)
    return ;
  bump = bumpCreate (mapGraphWidth-*offset, 0) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->flag & FLAG_HIDE)
	continue ;
      if (!(seg->flag & (FLAG_ANY_LOCUS | FLAG_MULTIPLE_LOCUS)))
	continue ;
      y = MAP2GRAPH(look->map,seg->x) ;
      if (y > topMargin && y < mapGraphHeight-1)
	{ 
	  x = errorScale*seg->dx ;
	  if ((n = bumpText (bump, name(seg->symbol), &x, &y, 1, TRUE)))
	    { char *ccp = name(seg->symbol), *ccq = ccp +n, ccc = *ccq ;
	      array(look->boxIndex,ibox=graphBoxStart(),int) = i ;
	      *ccq = 0 ; graphText (ccp, *offset+x, y-0.5) ; *ccq = ccc ; 
	      graphBoxDraw (ibox, BLACK, WHITE) ;
	      graphBoxEnd () ;
	      if (setCursorDrag)
		graphBoxMenu (ibox, vGeneMenu) ;
	    }
	}
    }

  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;
}

/*****************************************/

static void vMapDrawChromBands (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  float y1, y2 ;
  int   i, ibox ;
  SEG   *seg ;
  
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->flag & FLAG_HIDE)
	continue ;
      if (class(seg->key) != _VChrom_Band)
	continue ;
      y1 = MAP2GRAPH(look->map,seg->x) ;
      y2 = MAP2GRAPH(look->map,seg->x + seg->dx) ;
      if (y2 > topMargin && y1 < mapGraphHeight-1)
	{ array(look->boxIndex,ibox=graphBoxStart(),int) = i ;
	  if (y1 < topMargin) 
	    y1 = topMargin ;
	  if (y2 > mapGraphHeight -1)
	    y2 = mapGraphHeight -1 ;
	  graphRectangle (*offset+0.5,y1,*offset+1.5,y2) ;
	  graphBoxEnd() ;
	  if (seg->flag & FLAG_DARK_BAND)
	    graphBoxDraw (ibox, BLACK, DARKGRAY) ;
	  else if (seg->flag & FLAG_NOR)
	    graphBoxDraw (ibox, BLACK, LIGHTGRAY) ;
	  else if (seg->flag & FLAG_COLOUR)
	    graphBoxDraw (ibox, BLACK, (seg->flag >> 28) + WHITE) ;
	  else  
	    graphBoxDraw (ibox, BLACK, WHITE) ;
	}
    }
  *offset += 2.0 ;
}

/***********************/

MENUOPT vMapIntervalMenu[] = {
  {(VoidRoutine)nameSeg, "Show Name"},
  {(VoidRoutine)highlightSeg, "Highlight"},
  {(VoidRoutine)cursorToTop, "Cursor to top"},
  {(VoidRoutine)cursorToBottom, "Cursor to bottom"},
  {(VoidRoutine)topToCursor, "Move top to cursor"},
  {(VoidRoutine)bottomToCursor, "Move bottom to cursor"},
  {(VoidRoutine)vMapShowDataFromMenu, "Show data"},
  {0, 0} } ;

static int hideHeader ;

static float anyIntervalOffset ;

static void vMapDrawAnyInterval (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  float y1, y2, x ;
  int   i, ibox, ix ;
  SEG   *seg ;
  void *v ;
  BUMP  bump = bumpCreate (mapGraphWidth, 0) ;

  anyIntAss = assReCreate(anyIntAss) ;

  hideHeader = look->flag & FLAG_HIDE_HEADER ;

  anyIntervalOffset = *offset += 1.25 ;
  
  graphTextHeight (0.8) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->flag & FLAG_HIDE)
	continue ;
      if (class(seg->key) == _VChrom_Band) /* dedicated display */
	continue ;
      if (class(seg->key) == _VContig) /* dedicated display */
	continue ;
      if (!(seg->flag & FLAG_ANY_INTERVAL))
	continue ; 
      y1 = MAP2GRAPH(look->map,seg->x) ;
      y2 = MAP2GRAPH(look->map,seg->x + seg->dx) ;
      if (y2 > topMargin && y1 < mapGraphHeight-1)
	{ array(look->boxIndex,ibox=graphBoxStart(),int) = i ;
	  ix = 0 ; 
	  if (y1 < topMargin) y1 = topMargin ;
	  if (y2 > mapGraphHeight) y2 = mapGraphHeight ;
	  bumpItem (bump,1,(y2-y1+1)+0.2,&ix,&y1) ;
	  x = *offset + ix ;
	  if (y1 > topMargin) 
	    graphLine (x - 0.25, y1, x + 0.25, y1) ; 
	  else
	    { graphLine (x - 0.25, topMargin + .5, x, topMargin) ; 
	      graphLine (x + 0.25, topMargin + .5, x, topMargin) ; 
	    }
	  if (y2 < mapGraphHeight)
	    graphLine (x - 0.25, y2, x+0.25, y2) ; 
	  else
	    { graphLine (x - 0.25, mapGraphHeight - .5, x, mapGraphHeight) ; 
	      graphLine (x + 0.25, mapGraphHeight - .5, x, mapGraphHeight) ; 
	    }
	  graphLine (x, y1, x, y2) ;
	  graphBoxEnd () ;
	  v = (char*)(&anyIntAss) + ix ;
	  assInsert(anyIntAss, seg, v) ;
	  graphBoxDraw (ibox, BLACK, TRANSPARENT) ;
	  if (setCursorDrag)
	    graphBoxMenu (ibox, vMapIntervalMenu) ;
	}
    }

  graphTextFormat (PLAIN_FORMAT) ;
  graphTextHeight (0.0) ;

  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;
}

static void vMapDrawAnyIntervalNames (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  float y1,  y2, x, n ;
  int   i, ibox, ix , ix1, nmax = 0 ;
  SEG   *seg ;
  void *v ;
  BUMP  bump = bumpCreate (1, 0) ;

  hideHeader = look->flag & FLAG_HIDE_HEADER ;

  *offset += 1.25 ;

  graphTextHeight (0.8) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->flag & FLAG_HIDE)
	continue ;
      if (!(seg->flag & FLAG_ANY_INTERVAL))
	continue ; 
      y1 = MAP2GRAPH(look->map,seg->x) ;
      y2 = MAP2GRAPH(look->map,seg->x + seg->dx) ;
      if (y2 > topMargin && y1 < mapGraphHeight-1)
	{ 
	  n = 0.65*strlen(name(seg->symbol)) ;
	  if (n > nmax)
	    nmax = n ;
	  ix = 0 ; 
	  if (y1 < topMargin) y1 = topMargin ;
	  if (y2 > mapGraphHeight) y2 = mapGraphHeight ;
	  bumpTest (bump,1,1,&ix,&y1) ;
	  if (y1 > y2 + 3) /* too far down, drop it */
	    continue ;
	  bumpRegister (bump,1,1,&ix,&y1) ;
	  array(look->boxIndex,ibox=graphBoxStart(),int) = i ;
	  x = *offset + ix ;
	  if(y1 < mapGraphHeight)
	    graphText(name(seg->symbol), x, y1) ;
	  graphBoxEnd () ;
	  graphBoxDraw (ibox, BLACK, TRANSPARENT) ;
	  if (setCursorDrag)
	    graphBoxMenu (ibox, vMapIntervalMenu) ;
	  if (anyIntAss && assFind (anyIntAss, seg, &v))
	    { ix1 = (char*)v - (char*)(&anyIntAss) ;
	      y1 += .5 ;
	      if (y1 < mapGraphHeight &&
		  y1 < y2)
		graphCircle(anyIntervalOffset + ix1, y1, .4) ;
	    }
	}
    }

  graphTextFormat (PLAIN_FORMAT) ;
  graphTextHeight (0.0) ;

  *offset += nmax + 1 ;
  bumpDestroy (bump) ;
}

/***************************************************/

static void vMapDrawScale (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  float unit = 0.01 ;
  float subunit = 0.001 ;
  float x, xx, y ;
  char *cp = 0 ;
  MAP map = look->map ;

  mapFindScaleUnit (map, &unit, &subunit) ;
  if (unit >= 1)
    look->resolution = 0 ;
  else if (unit >= .1)
    look->resolution = 1 ;
  else if (unit >= .01)
    look->resolution = 2 ;
  else if (unit >= .001)
    look->resolution = 3 ;

  x = GRAPH2MAP(map, topMargin+1) ;
  x = unit * (((x>=0)?1:0) + (int)(x/unit)) ;
  while ((y = MAP2GRAPH(map, x)) < mapGraphHeight - 1)
    { graphLine (*offset+0.5,y,*offset+1.5,y) ;
      xx = look->map->flip ? -x : x ;
      cp = messprintf ("%-4.*f", look->resolution, xx) ;
      graphText (cp, *offset+2, y-0.5) ;
      x += unit ;
    }

  x = GRAPH2MAP(map, topMargin+1) ;
  x = subunit * (((x>=0)?1:0) + (int)(x/subunit)) ;
  while ((y = MAP2GRAPH(map, x)) < mapGraphHeight - 1)
    { graphLine (*offset+1.0,y,*offset+1.5,y) ;
      x += subunit ;
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

  if (setCursorDrag)
    mapCursorDraw (look->map, *offset) ;

  *offset += 2 ;
  if (cp)
    *offset += strlen(cp) ;
}

/***********************************************************/

static void vMapDrawMainMarkers (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  float y, xx, max = 0 ;
  int i ;
  SEG* seg ;

  graphTextHeight (0.75) ;
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->flag & FLAG_MARKER)
	{ y = MAP2WHOLE(look->map, seg->x) ;
	  graphText (name(seg->key),0.5,y-0.25) ;
	  xx = .65 * strlen(name(seg->key)) ;
	  if (xx > max)
	    max = xx ;
	}
    }
  graphTextHeight (0) ;
  *offset += max + .5 ;
}

/***********************************************************/

static void vMapDrawMiniChromBands (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  float y, dy, x = *offset + 1 , l=0, max = 0 , ylast = - 1000 ;
  int i ;
  SEG* seg ;
  BOOL isBands = FALSE ;
  MAP map = look->map ;

  graphTextHeight (0.75) ;
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (class(seg->key) == _VChrom_Band)
	{ isBands = TRUE ;
	  y = MAP2WHOLE(look->map, seg->x + seg->dx/2.) ;
	  dy =  0.5 * seg->dx * map->thumb.fac ;
	  if (seg->flag & FLAG_CENTROMERE)
	    { graphColor (WHITE) ;
	      graphFillRectangle (x, y-dy, x + 2., y+dy) ;
	      graphColor (BLACK) ;
	      graphLine (x, y-dy, x + 2., y+dy) ;
	      graphLine (x + 2., y-dy, x, y+dy) ;
	      graphLine (x, y-dy, x + 2., y-dy) ;
	      graphLine (x + 2., y+dy, x, y+dy) ;
	    }	      
	  else 
	    { if (seg->flag & FLAG_DARK_BAND)
		graphColor (DARKGRAY) ;
	      else if (seg->flag & FLAG_NOR)
		graphColor (LIGHTGRAY) ;
	      else
		graphColor (WHITE) ;
	      if (seg->flag >> 28)
		graphColor((seg->flag  >> 28) + WHITE) ;
	      graphFillRectangle (x, y-dy, x + 2., y+dy) ;
	      graphColor (BLACK) ;
	      graphRectangle (x, y-dy, x + 2., y+dy) ;
	      if (y > ylast + 4)
		{ graphText (name(seg->key),0.5,y-0.25) ;
		  l = .65 * strlen(name(seg->key)) ;
		}
	      ylast = y ;
	      if (l > max)
		max = l ;
	    }
	}
    }
 
  graphTextHeight (0) ;
  if (isBands)
    *offset += 2 + max ;
}

/**************************************************************/
/******************** conversion routines *********************/

int vMapOrder (const void *a, const void *b)  /* for arraySort() call */
{
  const SEG *seg1 = (const SEG*)a, *seg2 = (const SEG*)b ;
  float diff = seg1->x - seg2->x ;
  
  if (diff > 0)
    return 1 ;
  else if (diff < 0)
    return -1 ;
  else if (seg2->key != seg1->key)
    return seg2->key - seg1->key ;
  else
    return seg2->flag - seg1->flag ;
    /* never give a random result please */
}

/**********************************/

   /* the order in which the fathers are added
      implies that downMost father is the best
      and any son is better than his father,
      at least I hope so !
      */

static KEYSET addFather(KEY chrom)
{
  int n = 0, i, max;
  OBJ Chrom ;
  KEYSET a = keySetCreate(), a1 ; KEY key ;

  lexSetStatus(chrom, CALCULSTATUS) ;  /* prevent recursions */

  if ((Chrom = bsCreate(chrom)))
    { if (bsGetKey(Chrom, _Inherits_from, &key))
	do 
	  { if (!(CALCULSTATUS & lexGetStatus(key)))
	      { a1 = addFather(key) ;
		max = keySetMax(a1) ;
		for (i = 0 ; i < max ; i++)
		  keySet(a, n++) = keySet(a1, i++) ;
		keySetDestroy(a1) ;
	      }
	  } while (bsGetKey(Chrom, _bsDown, &key)) ;
      bsDestroy(Chrom) ;
      keySet(a, n++) = chrom ;
    }
  return a ;
}

/**********************************/

static Array vMapConvert (KEY chrom)
{
  static BSMARK mark = 0 ;
  KEY originalChrom = chrom;
  KEY clone, vMap, locus, col, symbol, symbolTag = 0 ; 
  OBJ Chrom, Locus;
  Array segs, loci, class2symbol ;
  KEYSET chKS ; int nCh ;
  SEG *seg = 0 ;
  float x, dx, x1, x2 , xmin = 9999999, xmax = -999999 ;
  int i = 0, iLocus, ns1 = 0, flag ;

  if (!(Chrom = bsCreate (chrom)))
    { messout ("No data for chromosome %s",name(chrom)) ;
      return 0 ;
    }
  if (bsFindTag(Chrom, _Non_graphic))
    { bsDestroy(Chrom) ;
      return 0 ;
    }
  bsDestroy(Chrom) ;

  if (lexReClass(chrom, &vMap,_VvMap) &&
      !lexlock (vMap))
    { messout ("Sorry, %s is locked (being processed elsewhere)",
	       name(vMap)) ; 
      return 0 ;
    }

  i = 256 ;
  while (i--)
    lexClearClassStatus(i, CALCULSTATUS); 
  segs = arrayCreate (128,SEG) ;

  chKS = addFather(chrom) ;
  if (!keySetMax(chKS))
    messout("No data found for this chromosome, sorry") ;

  loci = arrayCreate(100, BSunit) ;

/****** look recusivelly for symbols  *********/
  class2symbol = keySetCreate () ;
  keySet (class2symbol, 256) = 0 ;

  for (nCh = keySetMax(chKS) ; nCh-- ;)
    if ((chrom = keySet(chKS, nCh)) &&
	(Chrom = bsCreate(chrom)))
    { char *cls, *tg ;
      int cl ; KEY t ;
      if (bsFindTag (Chrom,_Symbol) &&
	  bsFlatten(Chrom, 2, loci))
	for (iLocus = 0 ; iLocus < arrayMax(loci) - 1 ; iLocus += 2)
	  { cls = arr(loci, iLocus, BSunit).s ;
	    tg  = arr(loci, iLocus + 1, BSunit).s ;
	    if ((cl = pickWord2Class (cls)) &&
		!keySet (class2symbol, cl) &&
		lexword2key (tg, &t, _VSystem))
	      keySet (class2symbol, cl) = t ;
	  }
      bsDestroy (Chrom) ;
    }
	
  i = 256 ;
  while (i--)
    if (!keySet (class2symbol, i))
      keySet (class2symbol, i) = pickList [class(i)].symbol ;

/*********************************************/
  i = 0 ;
  
  for (nCh = keySetMax(chKS) ; nCh-- ;)
    if ((chrom = keySet(chKS, nCh)) &&
	(Chrom = bsCreate(chrom)))
    { arrayMax(loci) = 0 ;
      if (bsFindTag (Chrom,_Contains) &&
	  bsFlatten(Chrom, 2, loci))
	for (iLocus = 1 ; iLocus < arrayMax(loci) ; iLocus += 2)
	  { locus = arr(loci, iLocus, BSunit).k ;
	      
	    if (!locus ||  (CALCULSTATUS & lexGetStatus(locus)) ||
		!(Locus = bsCreate (locus)))
	      continue ;
	    
	    symbolTag = keySet(class2symbol, class(locus)) ;
	    x = 10000.0 ; dx = 0 ;
	    if (bsFindKey (Locus, _Map, chrom) &&
		bsPushObj(Locus) )
	       { ns1 = i ;
		 if (bsGetData (Locus, _Position, _Float, &x))
		   { mark = bsMark(Locus,mark) ;
		     if (bsPushObj(Locus))
		       { bsGetData (Locus, _Error, _Float, &dx) ;
			 bsGoto(Locus, mark) ;
		       }
		     lexSetStatus(locus, CALCULSTATUS) ;
		     seg = arrayp(segs,i++,SEG) ;
		     seg->key = locus ;
		     seg->x = x ;
		     seg->dx = dx ;
		     seg->flag = FLAG_ANY_LOCUS ;
		     if (bsFindTag (Locus, _Well_ordered) ||
			 bsFindTag (Locus, _Main_Marker))
		       seg->flag |= FLAG_WELL_ORDERED ;

		     if (x < xmin)
		       xmin = x ;
		     if(x > xmax)
		       xmax = x ;
		   }  
		 symbol = 0 ;
		 if (symbolTag)
		   bsGetKey(Locus, symbolTag, &symbol) ;
		
		 if (bsGetData (Locus, _Multi_Position, _Float, &x))
		   do
		     { 
		       mark = bsMark(Locus,mark) ;
		       if (bsPushObj(Locus))
			 { bsGetData (Locus, _Error, _Float, &dx) ;
			   bsGoto(Locus, mark) ;
			 }
		       lexSetStatus(locus, CALCULSTATUS) ;
		       seg = arrayp(segs,i++,SEG) ;
		       seg->key = locus ;
		       if (symbol)
			 seg->symbol = symbol ;
		       else
			 seg->symbol = locus ;
		       seg->x = x ;
		       seg->dx = dx ;
		       seg->flag = FLAG_MULTIPLE_LOCUS ;
		       if (x < xmin)
			 xmin = x ;
		       if(x > xmax)
			 xmax = x ;
		       bsGoto(Locus, 0) ;
		       if (bsFindTag(Locus, _Main_Marker))
			 seg->flag |= FLAG_MARKER ;
		       bsGoto(Locus, mark) ;
		     } while (bsGetData(Locus, _bsDown, _Float, &x)) ;
		 
		 if (bsGetData (Locus, _Left, _Float, &x1) &&
		     bsGetData (Locus, _Right, _Float, &x2))
		   { 
		     lexSetStatus(locus, CALCULSTATUS) ;
		     seg = arrayp(segs,i++,SEG) ;
		     seg->key = locus ;
		     seg->x = x1 ;
		     seg->dx = x2 - x1 ;
		     if (seg->dx < 0)
		       seg->dx = -seg->dx ;
		     seg->flag = FLAG_ANY_INTERVAL ;
		     if (seg->flag >> 28)
		graphColor((seg->flag  >> 28) + WHITE) ;

		     if (x1 < xmin)
		       xmin = x1 ;
		     if(x2 > xmax)
		       xmax = x2 ;
		   }    
	       }
	    bsGoto(Locus,0) ;

	    if (i > ns1)
	      { flag = 0 ;
		if (!bsFindTag (Locus, _Compound))
		  { if (bsFindTag (Locus, _Duplication))
		      flag |= FLAG_DUPLICATION ;
		    if (bsFindTag (Locus, _Deletion))
		      flag |= FLAG_DEFICIENCY ;
		  }
		if (bsFindTag (Locus,_Centromere))
		  flag |= FLAG_CENTROMERE ;
		if (bsFindTag (Locus,_Dark))
		  flag |= FLAG_DARK_BAND ;
		if (bsFindTag (Locus,_NOR))
		  flag |= FLAG_NOR ;
		if (bsFindTag (Locus,_Colour) &&
		    bsGetKeyTags(Locus, _bsRight, &col))
		  { flag |= FLAG_COLOUR ;
		    flag |= ((col - _WHITE) & 15) << 28;
		  }
		if (bsFindKey(Locus, _Main_Marker, chrom))
		  flag |= FLAG_MARKER ;
		if (symbolTag && bsGetKey(Locus, symbolTag, &symbol))
		  seg->symbol = symbol ;
		else
		  seg->symbol = locus ;
		if (bsGetKey (Locus,_Positive_clone,&clone))
		  { seg->flag |= FLAG_CLONED ;
		    seg = arrayp(segs,i++,SEG) ;
		    seg->key = clone ;
		    seg->x = x ;
		    seg->dx = 0 ;
		    seg->flag = 0 ;
		  }
	      
		for (; ns1 < i ; ns1++)
		  array(segs,ns1,SEG).flag |= flag ;
	      }
		
	    bsDestroy (Locus) ;
	  }
      bsDestroy (Chrom) ;
    }

  arrayDestroy(loci) ;
  keySetDestroy(chKS) ;
  keySetDestroy(class2symbol) ;

  vMapAddPhysGenes (segs) ;
	
  if (arrayMax(segs) && (Chrom = bsCreate(originalChrom)))
    { float xx1, xx2, centre = 0, width = 10.0 ;

      seg = arrayp(segs,arrayMax(segs),SEG) ;
      seg->key = _Centre ;
      seg->x = (xmax + xmin)/2.;
      seg->dx = (xmax - xmin)/3. ;
      if (bsGetData (Chrom, _Centre, _Float, &centre))
	{ seg->x = centre ;
	  if (bsGetData (Chrom, _bsRight, _Float, &width))
	    seg->dx = width ;
	}

      seg = arrayp(segs,arrayMax(segs),SEG) ;
      seg->key = _Extent ;
      seg->x = xmin ;
      seg->dx = xmax - xmin ;
      if (bsGetData (Chrom, _Extent, _Float, &xx1))
	{ seg->x = xx1 ;
	  if (bsGetData (Chrom, _bsRight, _Float, &xx2))
	    seg->dx = xx2-xx1 ;
	}

      seg = arrayp(segs,arrayMax(segs),SEG) ;
      seg->key = _Int ;
      seg->x = -1000000 ;	/* to avoid boxIndex == 0 */
    
      arraySort (segs, vMapOrder) ;
      lexaddkey (name(originalChrom), &vMap, _VvMap) ;
      arrayStore (vMap, segs, segFormat) ;
      bsDestroy (Chrom) ;
    }
  else
    arrayDestroy(segs) ;
  
  if (vMap &&
      lexiskeylocked(vMap))
    lexunlock(vMap) ;

  return segs ;
}

/*********************************************************/

#define ZSEGDATA 	  if (seg->flag & FLAG_ANY_INTERVAL) \
	    { z1 = seg->x ; z2 = seg->x + seg->dx ; } \
	  else \
	    if (showMarginal) \
	      { z1 = seg->x - seg->dx ; z2 = seg->x + seg->dx ; } \
	    else \
	      z1 = z2 = seg->x 


void vMapSelect (VerticalMap look, int box)
{
  int i, j , k ;
  unsigned int  mask ;
  float x1, x2, z1, z2 ;
  SEG *seg ;
  KEYSET neighbours, inside, outside ;
  KEYSET contains, nocontains, overlap, nooverlap, a ; 
  KEY mainKey ;
  OBJ obj ;
  double dd = -12.60 ;
  float ff = -12.60, ee ;

/* ee is an estimate of rounding errors 
   when computing with float
*/
  ee = ff - dd ; if (ee < 0) ee = -ee ; ee *= 4 ;

  if (box < 1 ||
      box >= arrayMax(look->boxIndex) ||
      !(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)) ||
      !(a =  bsKeySet(seg->key)))
    { look->activeBox = 0 ;
      return ;
    }

  mainKey = seg->key ;
  look->activeBox = box ;

  if (seg->flag & FLAG_ANY_INTERVAL)
    { x1 = seg->x ; x2 = seg->x + seg->dx ; }
  else
    if (showMarginal)
      { x1 = seg->x - seg->dx ; x2 = seg->x + seg->dx ; }	     
    else
      { x1 = x2 = seg->x ; }

  mask = FLAG_STRESSED | FLAG_RELATED | FLAG_ANTI_RELATED |
    FLAG_ANTI_STRESSED ; 
  mask = ~mask ;
  i = arrayMax(look->segs) ; 
  seg = arrp(look->segs,0,SEG) - 1 ;
  while (seg++, i--)
    seg->flag &= mask ;

   neighbours = keySetCreate() ;

  for(i=0, j=0 ; i<keySetMax(a) ; i++)
    if (class(keySet(a,i)))
      keySet(neighbours,j++) = keySet(a,i) ;
  keySetDestroy(a) ;

  inside = keySetCreate() ;
  outside = keySetCreate() ;
  contains = keySetCreate() ;
  nocontains = keySetCreate() ;
  overlap = keySetCreate() ;
  nooverlap = keySetCreate() ;

  if ((obj = bsCreate(mainKey)))
    { a = arrayCreate(12, BSunit) ;
      if (bsFindTag(obj,_Inside) && bsFlatten(obj, 2, a))
	{ for (i=0; i < arrayMax(a) ; i+= 2 )
	    keySet(inside,i/2) = array(a,i+1,BSunit).k ;
	}
      arrayMax(a) = 0 ;

      if (bsFindTag(obj,_Contains) && bsFlatten(obj, 2, a))
	{ for (i=0; i < arrayMax(a) ; i+= 2 )
	    keySet(contains,i/2) = array(a,i+1,BSunit).k ;
	}
      arrayMax(a) = 0 ;
 
      if (bsFindTag(obj,_Overlaps) && bsFlatten(obj, 2, a))
	{ for (i=0; i < arrayMax(a) ; i+= 2 )
	    keySet(overlap,i/2) = array(a,i+1,BSunit).k ;
	}
      arrayMax(a) = 0 ;
 
      if (bsFindTag(obj,_Does_not_Overlap) && bsFlatten(obj, 2, a))
	{ for (i=0; i < arrayMax(a) ; i+= 2 )
	    keySet(nooverlap,i/2) = array(a,i+1,BSunit).k ;
	}
      arrayMax(a) = 0 ;
 
      if (bsFindTag(obj,_Does_not_Contain) && bsFlatten(obj, 2, a))
	{ for (i=0; i < arrayMax(a) ; i+= 2 )
	    keySet(nocontains,i/2) = array(a,i+1,BSunit).k ;
	}
      arrayMax(a) = 0 ;
 
      if (bsFindTag(obj,_Outside) && bsFlatten(obj, 2, a))
	{ for (i=0; i < arrayMax(a) ; i+= 2 )
	    keySet(outside,i/2) = array(a,i+1,BSunit).k ;
	}
      arrayDestroy(a) ;
      bsDestroy(obj) ;
    }

  keySetSort(inside) ;
  keySetSort(outside) ;
  keySetSort(contains) ;
  keySetSort(nocontains) ;
  keySetSort(overlap) ;
  keySetSort(nooverlap) ;
	
  i = arrayMax(look->boxIndex) ;
  while(--i)
    { k = arr(look->boxIndex, i, int) ;
      if (!k)  /* a control box */
	continue ;
      seg = arrp(look->segs,k,SEG) ;
      if (keySetFind(neighbours, seg->key, &j))
	seg->flag |= FLAG_RELATED ;

  /* x1-x2 are the extremities of the selected seg who contains... z1-z2 segs */
      if (keySetFind(inside, seg->key, &j))
	{
	  ZSEGDATA ;
	  if (x1 >= z1 - ee && x2 <= z2 + ee)
	    seg->flag |= FLAG_RELATED ;
	  else if (x1 > z2 + ee || x2 < z1 - ee)
	      seg->flag |= FLAG_PROBLEM | FLAG_STRESSED ;
	  else 
	    seg->flag |= FLAG_MARGINAL ;
	}		  
      if (keySetFind(outside, seg->key, &j))
	{
	  ZSEGDATA ;
	  if (x1 >= z2 + ee || x2 <= z1 - ee)
	    seg->flag |= FLAG_ANTI_RELATED ;
	  else if (x1 > z1 + ee && x2 <  z2 - ee)
	    seg->flag |= FLAG_PROBLEM | FLAG_ANTI_STRESSED ;
	  else 
	    seg->flag |= FLAG_MARGINAL ;
	}		  
      if (keySetFind(contains, seg->key, &j))
	{
	  ZSEGDATA ;
	  if (x1  <= z1 + ee && x2 >= z2 - ee)
	    seg->flag |= FLAG_RELATED ;
	  else if (x1 > z2 + ee || x2 < z1 - ee)
	      seg->flag |= FLAG_PROBLEM | FLAG_STRESSED ;
	  else 
	    seg->flag |= FLAG_MARGINAL ;
	}		  
      if (keySetFind(nocontains, seg->key, &j))
	{
	  ZSEGDATA ;
	  if (x1 > z2 - ee || x2 < z1 + ee)
	    seg->flag |= FLAG_ANTI_RELATED ;
	  else if (x1  <= z1 - ee && x2 >= z2 + ee)
	    seg->flag |= FLAG_PROBLEM | FLAG_ANTI_STRESSED ;
	  else 
	    seg->flag |= FLAG_MARGINAL ;
	}		  
      if (keySetFind(overlap, seg->key, &j))
	{
	  ZSEGDATA ;
	  if ((x1  >= z1 - ee && x1 <= z2 + ee) ||
	      (x2 >= z1 - ee && x2 <= z2 + ee))
	    seg->flag |= FLAG_RELATED ;
	  else 
	    seg->flag |= FLAG_PROBLEM | FLAG_STRESSED ;
	}		  
      if (keySetFind(nooverlap, seg->key, &j))
	{
	  ZSEGDATA ;
	  if (x1 > z2 - ee  || x2 < z1 + ee)
	    seg->flag |= FLAG_ANTI_RELATED ;
	  else 
	    seg->flag |= FLAG_PROBLEM | FLAG_ANTI_STRESSED ;
	}		  
      
      if(seg->key == mainKey) 
	seg->flag |= FLAG_RELATED ;
    }
  keySetDestroy(neighbours) ;
  keySetDestroy(inside) ;
  keySetDestroy(outside) ;
  keySetDestroy(contains) ;
  keySetDestroy(nocontains) ;
  keySetDestroy(overlap) ;
  keySetDestroy(nooverlap) ;

  vMapReColor (look) ;
  if (class (mainKey) == _V2_point_data ||
      class (mainKey) == _VMulti_pt_data)
    graphBoxDraw(look->activeBox, BLACK, YELLOW) ;
  else
    graphBoxDraw (look->activeBox, BLACK, CYAN) ;

  if (look->messageBox)
    { seg = arrp(look->segs, arr(look->boxIndex,box,int), SEG) ;
      if (class (seg->key) == _VMulti_pt_data)
	vMapReportMultiPt (look, seg) ;
      else  if (class (seg->key) == _V2_point_data)
	  vMapReport2pt (look, seg) ;
      else  if (class (seg->key) == _VClone ||
		class (seg->key) == _VContig)
	{ strncpy (look->messageText, name (seg->key), 100) ;
	  strcat (look->messageText, messprintf(" %.2f", seg->x)) ;
	  if (class(seg->key) == _VContig)
	    strcat (look->messageText, messprintf (" %.2f", 
						   seg->x + seg->dx)) ;
	}
      else  if (class (seg->key) == _VLocus ||
		class (seg->key) == _VInterval)
	{
	  strncpy (look->messageText, name (seg->key), 100) ;
	  strcat (look->messageText, messprintf(" %.2f", seg->x)) ;
	  if (class(seg->key) == _VLocus)
	    {
	      if (seg->flag & FLAG_PHYS_GENE)
		strcat (look->messageText, messprintf (" interpolated")) ;
	      else
		strcat (look->messageText, messprintf (" [%.2f]", 
						       seg->dx)) ;
	    }
	  else
	    strcat (look->messageText, messprintf (" %.2f", 
						   seg->x + seg->dx)) ;
	}
      else
	{
	  strncpy (look->messageText, name (seg->key), 100) ;
	  strcat (look->messageText, messprintf(" %.2f  %.2f", seg->x, seg->dx)) ;
	}
      
      graphBoxDraw (look->messageBox, BLACK, CYAN) ;
    }
}


static void vMapSearchProblems(VerticalMap look) 
{
  int i, j , k , k1 ;
  float x1, x2, z1, z2 ;
  SEG *seg ;
  KEYSET  inside = 0 , outside = 0, contains = 0, a ;
  KEY mainKey ;
  OBJ obj ;

  messStatus("Searching Map Problems") ;
   k1 = arrayMax(look->segs) ;
  while(k1--)
    { seg = arrp(look->segs,k1,SEG) ;
      mainKey = seg->key ;
      if (seg->flag & FLAG_ANY_INTERVAL)
	{ x1 = seg->x ; x2 = seg->x + seg->dx ; }
      else
	{ x1 = seg->x - seg->dx ; x2 = seg->x + seg->dx ; }	    
      
      inside = keySetReCreate(inside) ;
      outside = keySetReCreate(outside) ;
      contains = keySetReCreate(contains) ;
      if ((obj = bsCreate(mainKey)))
	{ a = arrayCreate(12, BSunit) ;
	  if (bsFindTag(obj,_Inside) && bsFlatten(obj, 2, a))
	    { for (i=0; i < arrayMax(a) ; i+= 2 )
		keySet(inside,i/2) = array(a,i+1,BSunit).k ;
	    }
	  arrayMax(a) = 0 ;
	  
	  if (bsFindTag(obj,_Contains) && bsFlatten(obj, 2, a))
	    { for (i=0; i < arrayMax(a) ; i+= 2 )
		keySet(contains,i/2) = array(a,i+1,BSunit).k ;
	    }
	  arrayMax(a) = 0 ;
	  
	  if (bsFindTag(obj,_Outside) && bsFlatten(obj, 2, a))
	    { for (i=0; i < arrayMax(a) ; i+= 2 )
		keySet(outside,i/2) = array(a,i+1,BSunit).k ;
	    }
	  arrayDestroy(a) ;
	  bsDestroy(obj) ;
	}
      
      keySetSort(inside) ;
      keySetSort(outside) ;
      keySetSort(contains) ;
      
      k = arrayMax(look->segs) ;
      while(k--)
	{ seg = arrp(look->segs,k,SEG) ;
	  
	  /* x1-x2 are the extremities of the selected seg who contains... z1-z2 segs */
	  if (keySetFind(inside, seg->key, &j))
	    {
	      if (seg->flag & FLAG_ANY_INTERVAL)
		{ z1 = seg->x ; z2 = seg->x + seg->dx ; }
	      else
		{ z1 = seg->x - seg->dx ; z2 = seg->x + seg->dx ; }	    
	      if (x1 > z2 || x2 < z1)
		seg->flag |= FLAG_PROBLEM ;
	    }		  
	  if (keySetFind(contains, seg->key, &j))
	    {
	      if (seg->flag & FLAG_ANY_INTERVAL)
		{ z1 = seg->x ; z2 = seg->x + seg->dx ; }
	      else
		{ z1 = seg->x - seg->dx ; z2 = seg->x + seg->dx ; }	    
	      if (x1 > z2 || x2 < z1)
		seg->flag |= FLAG_PROBLEM ;
	    }		  
	  if (keySetFind(outside, seg->key, &j))
	    {
	      if (seg->flag & FLAG_ANY_INTERVAL)
		{ z1 = seg->x ; z2 = seg->x + seg->dx ; }
	      else
		{ z1 = seg->x - seg->dx ; z2 = seg->x + seg->dx ; }	    
	      if (x1 > z1 && x2 <  z2)
		seg->flag |= FLAG_PROBLEM ;
	    }		  
	}
    }

  keySetDestroy(contains) ;
  keySetDestroy(inside) ;
  keySetDestroy(outside) ;
}

/*****************************************************************************************/

static void vMapFollow (VerticalMap look, float x, float y)
{ KEY displayKey = 0 ;
  SEG *seg = arrp(look->segs, 
		  arr(look->boxIndex,look->activeBox,int),
		  SEG) ;

   if (class(seg->key) == _VContig)
     vMapToPMap (look, y) ;
   else
     {
       lexword2key (VMAP, &displayKey, _VDisplay) ;
       if (class(seg->key) == _VMap)
	 display (seg->key, 0, 0) ;
       else if (strcmp(name(pickDisplayKey (seg->key)), VMAP))
	 display (seg->key, look->key, 0) ;
       else
	 display (seg->key, look->key, TREE) ;
     }
}

/************* richard's public Save Map ************/

static void vMapSaveMap (void)
{
  KEY map ;
  OBJ obj ;
  int i, oldMax ;
  int nLocus = 0, nInterval = 0, nContig = 0 ;
  SEG *seg ;
  float x ;
  VerticalMap look = currentVerticalMap("vMapSaveMap") ;
  
  if (!isWriteAccess())
    { messout("Sorry, you do not have Write Access");
      return ;
    }

  lexReClass (look->key, &map, _VMap) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if ( (seg->flag & FLAG_MOVED) &&
	  (obj = bsUpdate (seg->key)) )
	{ bsAddKey (obj, _Map, map) ;
	  bsPushObj (obj) ;
	      
	  if (seg->flag & FLAG_ANY_INTERVAL)
	    { bsAddData (obj, _Left, _Float, &seg->x) ;
	      x = seg->x + seg->dx ;
	      bsAddData (obj, _Right, _Float, &x) ;
	      ++nInterval ;
	    }
	  else
	    { bsAddData (obj, _Position, _Float, &seg->x) ;
	      if (bsPushObj(obj))
		bsAddData (obj, _Error, _Float, &seg->dx) ;
	      ++nLocus ;
	    }
	  bsSave (obj) ;
	  seg->flag &= ~FLAG_MOVED ;
	}
    }

  messout ("Saved map and %d loci, %d intervals, %d contigs",
	   nLocus, nInterval, nContig) ;

  oldMax = arrayMax(look->segs) ;
  arrayMax(look->segs) = look->lastTrueSeg ;
  arrayStore (look->key, look->segs, segFormat) ;
  arrayMax(look->segs) = oldMax ;
}

/**************************************************************/
/*********  Private Maps **************************************/
/*
static void vMapSaveDetails(Array segs, KEY chromo)
{ int i = arrayMax(segs) , n = 0 ;
  SEG *seg = arrp(segs, 0, SEG) - 1 ;
  OBJ obj ;
  float x1, x2 ;

  while(seg++, i--)
    if (seg->flag & FLAG_MOVED)
      {if (obj = bsUpdate(seg->key))
	 { bsAddKey(obj, _Map, chromo) ;
	   if (bsPushObj(obj))
	     if (seg->flag & FLAG_ANY_INTERVAL)
	       { x1 = seg->x ;
		 x2 = seg->x + seg->dx ;
		 bsAddData(obj, _Left, _Float, &x1) ;
		 bsAddData(obj, _Right, _Float, &x2) ;
	       }
	     else
	       { x1 = seg->x ;
		 x2 = seg->dx ;
		 bsAddData(obj, _Position, _Float, &x1) ;
		 if (bsPushObj(obj))
		   bsAddData(obj, _Error, _Float, &x2) ;
	       }
	   bsSave(obj) ;
	   seg->flag &= ~FLAG_MOVED ;
	 }
       else
	 messout("Sorry %s is locked elsewhere", name(seg->key)) ;
       n++ ;
     }
  messout("I edited %d objects", n) ;
}
*/ 
/*
static void vMapSaveMapPrivate (void)
{ KEY key , oldCh , vMap ; 
  char *cp , buffer[41] ;
  OBJ obj ;
  VerticalMap look = currentVerticalMap("vMapSaveMap") ;
 
  
  if(!isWriteAccess())
    { messout("Sorry, you do not have Write Access");
      return ;
    }

  if (!messPrompt ("Save as :",name(look->key),"t"))
    return ;
	
  strncpy(buffer, freeword(), 40) ;
  if (lexword2key(buffer, &key,_VMap) &&
	!messQuery ("Overwrite existing map ?"))
      return ;
  
  lexaddkey(buffer, &key,_VMap) ;
  
  if (!key)
    return ;
  
  if (obj = bsCreate(key))
    { if (!bsFindTag(obj, _Author) ||  // as in the official chromo 
	  bsGetData(obj, _Author, _Text, &cp) &&
	  strcmp(cp, thisSession.name))
	{ messout ("This chromo is not yours, save it under some other name, sorry") ;
	  bsDestroy(obj) ;
	  return ;
	}
        bsDestroy(obj) ;
    }

  lexReClass(look->key, &oldCh, _VMap) ;
  if (key != oldCh) 
    { obj = bsUpdate(key) ;
      bsAddKey(obj, _Inherits_from, oldCh) ;
      bsAddData(obj, _Author, _Text, thisSession.name) ;
      bsSave(obj) ;
    }

  vMapSaveDetails(look->segs, key) ;      

  lexaddkey(name(key), &vMap, _VvMap) ;
  if (!lexlock(vMap))
    { messout ("Sorry, vMap %s is locked by someone else", name(vMap)) ;
      return ;
    }

  arrayStore (vMap,look->segs,segFormat) ;
  lexunlock (vMap) ;
  look->key = vMap ;
}
*/

/***************************************************************/

static BOOL pMapToMyVMap (KEY contig, KEY from, int x)
{
  KEY vMap, map ;
  OBJ obj ;
  SEG *seg ;
  float y1 = -1000000, y2 = 1000000, y ;
  int i, p1 = -1000000, p2 = 1000000, x1, x2 ;
  Array segs ;
  BOOL inContig = FALSE ;

  if (!(obj = bsCreate(contig)))
    { /* messout ("Can't open contig object") ; */
      return FALSE ;
    }
  if (!bsGetKey (obj, _Map, &map) || !map)
    { bsDestroy (obj) ;
     /*  messout ("This contig is not assigned to a chromosome") ; */
      return FALSE ;
    }
  bsGetData (obj, _pMap, _Int, &p1) ;
  bsGetData (obj, _bsRight, _Int, &p2) ;
  bsDestroy (obj) ;

		/* next get the segs */
  if (!lexReClass (map, &vMap, _VvMap) ||
      !(segs = arrayGet (vMap, SEG, segFormat)))
    {/* messout ("Please recalculate the genetic map of %s",
	       name(map)) ;
      */
      return FALSE ;
    }
  
  if (from)			/* look for it in the segs */
    for (i = arrayMax(segs) ; i-- ;)
      if (arr(segs,i,SEG).key == from)
	{ display (vMap, from, 0) ;
	  arrayDestroy (segs) ;
	  return TRUE ;
	}

 /* strategy: find 2 vMapped clones around or by x and interpolate */

  for (i = 0 ; i < arrayMax(segs) ; ++i)
    { seg = arrp(segs, i, SEG) ;
      if (seg->key == contig)
	{ inContig = TRUE ;
	  y1 = seg->x ; 
	  y2 = seg->x + seg->dx ;
	}
      else if (inContig &&
	       seg->x > y1 && seg->x < y2 &&
	       class(seg->key) == _VClone &&
	       (obj = bsCreate(seg->key)))
	{ if (bsFindKey (obj, _pMap, contig) &&
	      bsGetData (obj, _bsRight, _Int, &x1) &&
	      bsGetData (obj, _bsRight, _Int, &x2))
	    { x1 = 0.5 * (x1+x2) ;
	      if (x1 < x)
		{ y1 = seg->x ; p1 = x1 ; }
	      else
		{ y2 = seg->x ; p2 = x1 ; }
	    }
	  bsDestroy (obj) ;
	}
    }
  arrayDestroy (segs) ;

  if (!inContig || y1 < -999999 || y2 > 999999)
    { /* messout ("Sorry, I could not find two flanking markers") */
      return FALSE ;
    }								     
  else
    { y = y1 + (y2 - y1) * (x - p1) / (p2 - p1) ;
	from = KEYMAKE(_VCalcul, 1000.0*(y + 1000.0)) ;
      display (vMap, from, 0) ;
      return TRUE ;
    }
}

void pMapToMyGMap (KEY contig, KEY from, int x)
{
  extern void pMapToGMap (KEY contig, KEY from, int x);	/* gmapsys.c */

  if (!pMapToMyVMap (contig,from, x))
    pMapToGMap (contig,from, x) ;
}

/***************************************************************/
/************* end of file ****************/

 
