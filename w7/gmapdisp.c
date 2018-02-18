/*  File: gmapdisp.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Display of the genetic map
 * Exported functions:
 *     gMapDisplay
 *     gMapMakeAll
 *     gMapHighlightKey
 *     gMapUnselect
 *     gMapColourBox
 *     gMapOrder
 *     gMapConvert
 * HISTORY:
 * Last edited: Dec 21 13:56 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: gmapdisp.c,v 1.5 2015/08/12 14:11:59 mieg Exp $ */

#include "gmap.h"
#include "main.h"		/* for checkWriteAccess() */

magic_t GeneticMap_MAGIC = "GeneticMap";

static COLPROTO columnsTab[] = {
  &gMapMainMarkersColumn,
  &gMapMiniChromBandsColumn,
  &mapLocatorColumn,
  &mapScrollerColumn,
  &gMapJTMIntervalColumn,
  &gMapRDIntervalColumn,
  &gMapChromBandsColumn,
  &gMapContigColumn,
  &gMapRevPhysColumn,
  &gMapPhysGenesColumn,
  &mapScaleColumn,
  &gMapTwoPtColumn,
  &gMapMultiPtColumn,
  &gMapLikelihoodColumn,
  &gMapPointColumn,
  &gMapRemarkColumn,
  &spacerColumn,
  &gMapPosNegColumn,
  &gMapSubMapColumn,
  0
};


static COLDEFAULT defaultColumns = { 
  { &gMapMainMarkersColumn, TRUE, "Marker_loci" },
  { &gMapMiniChromBandsColumn, TRUE, "Chromsome_picture" },
  { &mapScrollerColumn, FALSE, "Scroller" },
  { &mapLocatorColumn, TRUE, "Locator" },
  { &gMapJTMIntervalColumn, TRUE, "Intervals" },
  { &gMapContigColumn, TRUE, "Contigs" },
  { &gMapRevPhysColumn, TRUE, "Reversed_physical" },
  { &gMapPhysGenesColumn, TRUE, "Physical_genes" },
  { &mapScaleColumn, TRUE, "Scale" },
  { &gMapTwoPtColumn, FALSE, "Two_point_data" },
  { &gMapMultiPtColumn, FALSE, "Multi_point_data" },
  { &gMapLikelihoodColumn, FALSE, "Likelihood_distribution" },
  { &gMapPointColumn, TRUE, "Loci" },
  { &gMapSubMapColumn, FALSE, "submap" },
  { 0, 0, 0}
};

/************************************************************/


GeneticMap currentGeneticMap (char *caller)
{
  GeneticMap gMap;

  gMap = (GeneticMap)(currentMapControl())->look;

  if (!gMap)
    messcrash("%s() received NULL GeneticMap pointer", caller);
  if (gMap->magic != &GeneticMap_MAGIC)
    messcrash("%s() received non-magic GeneticMap pointer", caller);
  
  return gMap;
} /* currentGeneticMap */



void gMapMakeAll (void)
{ KEY chromosome = 0;
  BOOL isold = FALSE;

  if (!isWriteAccess ())
    { messout ("Sorry, you do not have Write Access") ;
      return ;
    }

  while (lexNext (_VgMap,&chromosome))
    arrayKill(chromosome);
  
  chromosome = 0;

  while (lexNext (_VMap,&chromosome))
    { gMapDisplay(chromosome, 0, isold);
      isold = TRUE;
      if (messIsInterruptCalled ())
	break ;
    }

  if (isold)
    graphDestroy();

}




/********************************/
/* Main menu                    */
/********************************/
static void flip(void)
{ 
  MAPCONTROL map = currentMapControl();
  map->mag = -map->mag;
  controlDraw();
}

static void gMapGotoCmap (void)
{
  MAPCONTROL map = currentMapControl() ;
  display (map->key, 0, CMAP) ;
}


static void gMapRecalculate(void)
{
  MAPCONTROL map = currentMapControl() ;
  KEY key ;
  float cen = map->centre, mag = map->mag ;

  /* If multi map then parentKey will be set, else it won't*/
  if(map->parentKey)
    key = map->parentKey;
  else
    key= map->key;

  /* destroy a cached set of segs */
  gMapCacheKill(key);

  map->key = 0 ;
  gMapDisplay (key, 0, TRUE) ;	/* reuse the graph */

  map = currentMapControl() ;
  map->centre = cen ;
  map->mag = mag ;
  controlDraw() ;
}


static void gMapSaveMap (void)
{
  GeneticMap look = currentGeneticMap("gMapSaveMap");
  
  if (!checkWriteAccess())
    return ;

  gMapSaveDetails (look->segs, look->map->key) ;      
}

static void gMapSaveMapPrivate (void)
{
  KEY key;
  char *cp , buffer[41] ;
  OBJ obj ;
  MAPCONTROL map = currentMapControl();
  GeneticMap look = currentGeneticMap("gMapSaveMapPrivate");
  
  if (!checkWriteAccess())
      return ;

  if (!messPrompt ("Save as :", name(map->key), "t"))
    return ;
	
  strncpy(buffer, freeword(), 40) ;
  if (lexword2key(buffer, &key, _VMap) &&
	!messQuery ("Overwrite existing map ?"))
      return ;
  
  lexaddkey (buffer, &key, _VMap) ;
  
  if (!key)
    return ;
  
  if ((obj = bsCreate(key)))
    { if (!bsFindTag(obj, _Author) ||  /* as in the official chromo */
	  (bsGetData(obj, _Author, _Text, &cp) &&  strcmp(cp, thisSession.name))
	  )
	{ messout ("This map is not yours, save it under some other name, sorry") ;
	  bsDestroy(obj) ;
	  return ;
	}
      bsDestroy(obj) ;
    }

  if (key != map->key) 
    { obj = bsUpdate(key) ;
      bsAddKey(obj, _From_map, map->key) ;
      bsAddData(obj, _Author, _Text, thisSession.name) ;
      bsSave(obj) ;
    }

  gMapSaveDetails (look->segs, key) ;      
}

static void setMagnification (void)
{ 
  MAPCONTROL map = currentMapControl();
  float old = map->mag ; /* avoid zero, mhmp, dec 98 */

  if (graphPrompt ("Give new magnification", messprintf ("%.2f", map->mag), "fz"))
    freefloat (&map->mag) ;
  if (!map->mag)
    map->mag = old ;
  else
    controlDraw () ;
}

static void removeThisMap(void)
{ 
  controlDeleteMap(currentColControl("removeThisMap"),
		   currentMapControl());
}

static void hideHeaderToggle()
{
  MAPCONTROL map = currentMapControl();
  COLCONTROL control = map->control;

  control->hideHeader = !control->hideHeader;
  controlDraw () ;
}

static MENUOPT gMapMenu[] = {
  { graphDestroy, "Quit" },
  { help, "Help" },
  { graphPrint, "Print Screen" },
  { controlPrint,"Print Whole Map" },
  { removeThisMap, "Delete Current Map" },
  { setMagnification, "Set Magnification" },
  { displayPreserve, "Preserve" },
  { gMapRecalculate, "Recalculate" },
  { flip, "Vertical Flip" },
  { gMapSaveMap, "Save Map" },
  { gMapSaveMapPrivate, "Private Save" },  
  { gMapGotoCmap, "Physical Chromo Map" },
  { hideHeaderToggle, "Hide Header" },
  { 0, 0 }
} ;


/*************************************************************/
/******** gMap show intersect with active keyset *************/
/*************************************************************/ 


static void gMapClear (void)
{
  GeneticMap look = currentGeneticMap("gMapClear");

  look->highlight = keySetReCreate(look->highlight);
  look->hidden = keySetReCreate(look->hidden);

  controlDraw(); /* deals with stressed et al */

  return;
} /* gMapClear */

static void gMapHighlightAll (void)
{
  int i ;
  GeneticMap look = currentGeneticMap("gMapHighlightAll");

  look->highlight = keySetReCreate(look->highlight);
  
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    keySet(look->highlight, keySetMax(look->highlight)) =
      arrp(look->segs,i,GMAPSEG)->key;

  keySetSort(look->highlight);
  keySetCompress(look->highlight);
  controlUnselectAll();
}

static void gMapHighlight (void)
{
  void *dummy ;
  KEYSET kSet = 0 ;
  GeneticMap look = currentGeneticMap("gMapHighlight");

  if (!keySetActive(&kSet, &dummy))
    { messout("First select a keySet window, thank you.") ;
      return ;
    }

  keySetDestroy(look->highlight);
  look->highlight = keySetCopy(kSet);
  keySetSort(look->highlight);
  
  controlDraw();
}

static void gMapUnHighlight (void)
{
  void *dummy ;
  KEYSET kSet = 0 ;
  KEYSET tmp;
  GeneticMap look = currentGeneticMap("gMapUnHighlight");

  if (!keySetActive(&kSet, &dummy))
    { messout("First select a keySet window, thank you.") ;
      return ;
    }

  tmp = look->highlight;

  look->highlight = keySetMINUS(tmp, kSet);
  keySetDestroy(tmp);

  controlDraw();
}

static void gMapExportHighlit (void)
{
  KEYSET kset ;
  GeneticMap look = currentGeneticMap("gMapExportHighlit");

  kset = keySetCopy(look->highlight) ;

  displayCreate (DtKeySet) ;	/* new window */
  keySetShow (kset, 0) ;
}

static void gMapHideHighlit (void)
{
  GeneticMap look = currentGeneticMap("gMapHideHighlit");
  KEYSET tmp;

  tmp = look->hidden;
  look->hidden = keySetOR(look->hidden, look->highlight);
  keySetDestroy(tmp);

  controlDraw() ;
}

static void gMapHideUnHighlit (void)
{
  int i ;
  GeneticMap look = currentGeneticMap("gMapHideUnHighlit");
  KEYSET a = keySetCreate (), b, c ;


  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    keySet(a, keySetMax(a)) = arrp(look->segs,i,GMAPSEG)->key;
  keySetSort(a) ;
  keySetCompress(a) ;

  b = keySetMINUS(a,look->highlight) ;

  c = look->hidden ;
  look->hidden = keySetOR(b, c);

  keySetDestroy(a) ;
  keySetDestroy(b) ;
  keySetDestroy(c) ;

  controlDraw() ;

  return;
} /* gMapHideUnHighlit */

static void gMapUnHide (void)
{
  GeneticMap look = currentGeneticMap("gMapUnHide");

  look->hidden = keySetReCreate(look->hidden);
  controlDraw();

  return;
} /* gMapUnHide */



void gMapHighlightKey(MAPCONTROL map, KEY key) 
/* public routine */
{
  GeneticMap look = (GeneticMap)(map->look);
  
  keySet(look->highlight, keySetMax(look->highlight)) = key;
  
  keySetSort(look->highlight);
  keySetCompress(look->highlight);

  controlUnselectAll(); /* So that object is magenta, not cyan after redraw */
  controlDrawControl(map->control); 

  return;
} /* gMapHighlightKey */


static MENUOPT clearMenu[] = {
  { gMapClear, "Clear and show all" },
  { gMapHighlightAll, "Highlight all" },
  { gMapHighlight, "Highlight selected keyset" },
  { gMapUnHighlight, "Unhighlight selected keyset" },
  { gMapHideHighlit, "Hide highlit items" },
  { gMapHideUnHighlit, "Hide unhighlit items" },
  { gMapUnHide, "Unhide hidden items" },
  { gMapExportHighlit, "Export highlit keyset" },
  { 0, 0 }
} ;

static COLOUROPT highlightButton[] = {
  { (ColouredButtonFunc)gMapHighlight, 0, BLACK, WHITE, "Highlight...", 0 }
};

/* Note, as the button to which the above Menu is attached is only
   drawn when it is the current map (it's in the header), we can 
   rely on currentMap() to get us our map.
*/
/**************************************************/
/*     gMap display code                          */
/**************************************************/

static void gMapFinalise(void *p)
{
  int i ;
  MAPCONTROL map = (MAPCONTROL)p;
  GeneticMap look = (GeneticMap)(map->look);

  if (isWriteAccess())
    for (i = 0 ; i < arrayMax(look->segs) ; ++i)
      if (arrp(look->segs,i,GMAPSEG)->flag & FLAG_MOVED)
	{ if (messQuery ("Do you want to save your changes?"))
	    gMapSaveDetails (look->segs, map->key) ;
	  break ;
	}

  /* destroy in order of declaration in gmap.h */
  /* destroyed on map handle 
     keySetDestroy(look->orderedLoci) ;
     assDestroy (look->posAss);
     */
  keySetDestroy(look->highlight);
  keySetDestroy(look->hidden);
  keySetDestroy(look->positive);
  keySetDestroy(look->negative);
  keySetDestroy(look->neighbours);
  keySetDestroy(look->twoPt);
  keySetDestroy(look->multi);
  /* 
     keySetDestroy(look->negativeMarkers); never allocated directly 
     */

  arrayDestroy (look->segs);
  arrayDestroy (look->remarkSegs);

  look->magic = 0 ;
  messfree(look);
}

static void gMapKybd (int k)
{
  MAPCONTROL map = currentMapControl();

  switch (k)
    {
    case UP_KEY :
      if (map->cursor.box)
	{ if (map->mag > 0)
	    --map->cursor.val ;
	  else
	    ++map->cursor.val ;
	  mapControlCursorShift (map) ;
	}
      break ;
    case DOWN_KEY :
      if (map->cursor.box)
	{ if (map->mag > 0)
	    ++map->cursor.val ;
	  else
	    --map->cursor.val ;
	  mapControlCursorShift (map) ;
	}
      break ;
    }
}

static float gMapDrawHeader (MAPCONTROL map)
{
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);
  int x ;
  float y ;
  
  *look->messageText = 0;

  look->aboutBox = graphBoxStart();
  graphTextFormat (BOLD) ;
  if (look->titleKey)
    { x = strlen (name(look->titleKey)) ;
      graphText (name(look->titleKey),1,0.5) ;
    }
  else
    { x = strlen (name(map->key)) ;
      graphText (name(map->key),1,0.5) ;
    }
  graphTextFormat (PLAIN_FORMAT) ;
  graphBoxEnd();
  
  if (map->noButtons) /* may suppress header buttons for WWW */
    y = 1.4;
  else
    y = controlDrawButtons (map, x+2, 0.5, control->graphWidth) ; 

  look->messageBox = graphBoxStart () ;
  graphTextPtr (look->messageText, 1, y + 0.25, MAXMESSAGETEXT-1) ;
  graphBoxEnd () ;
  graphBoxDraw (look->messageBox, BLACK, PALECYAN) ;

  return y + 1.5 ;
}

void gMapAddToHeader(MAPCONTROL map, char *string)
{
  char *mt = ((GeneticMap)map->look)->messageText;

  strncat(mt, string, 
	  strlen(mt)<MAXMESSAGETEXT-2 ? MAXMESSAGETEXT-(2+strlen(mt)) : 0);
}

static void gMapHeaderPick (MAPCONTROL map, int box, double x, double y)
{
  GeneticMap look = (GeneticMap)(map->look);
  
  if (box == look->aboutBox)
    display(map->key, 0, TREE);
}

static BOOL getCentre(KEY from, KEY map, float *ret)
{ if (class(from) == _VCalcul)
    { *ret = (float)(KEYKEY(from))/1000.0 - 1000;
      return TRUE;
    }
  
  if (gMapGetMapObject(from, map, 0, 20, ret, 0, 0, 0))
    return TRUE;
  else
    return FALSE;
}

static MAPCONTROL gMapDoMap(COLCONTROL control, 
			    KEY key, KEY from, KEY view, int nmaps)
{
  MAPCONTROL map = mapControlCreate (control, gMapFinalise) ;
  GeneticMap look ;
  float centre;
  OBJ obj;

  /* look freed in gMapFinalise to make sure it's around when needed */
  look = (GeneticMap) messalloc (sizeof (struct GeneticMapStruct));
  
  map->look = (void *)look;
  look->map = map;  /* how incestuous! */
  
  look->magic = &GeneticMap_MAGIC;
  look->highlight = keySetCreate();
  look->hidden = keySetCreate();
  look->positive =  keySetCreate();
  look->negative  = keySetCreate();
  look->neighbours = keySetCreate();
  look->multi = keySetCreate();
  look->twoPt = keySetCreate();
  /* messalloc automatically initialise to zero 
  look->negativeMarkers = 0;
  look->posAss = 0;
  look->titleKey = 0;
  look->remarkSegs = 0;
  look->segs = 0;
  */

  if ((obj = bsCreate(key)))
    { bsGetKey(obj, _Title, &look->titleKey) ;
      bsDestroy(obj) ;
    }
    map->displayName = GMAP;
  if(nmaps ==1){
    map->name = handleAlloc(0, map->handle, 1+strlen(name(key)));
    strcpy(map->name, name(key));
    map->parentKey = 0;
  }
  else{
    map->name = handleAlloc(0, map->handle, 1+strlen(name(from)));
    strcpy(map->name, name(from));
    map->parentKey = from;
  }  
  map->key = key;
  map->prototypes = columnsTab;
  map->menu = gMapMenu;
  map->drawHeader = gMapDrawHeader;
  map->headerPick = gMapHeaderPick;
  map->keyboard = gMapKybd;
  map->topMargin = 3;
  map->viewTag = str2tag("Gmap");
  controlAddMap(control, map);

  getPos(map, 0, 0); /* init map->min, map->max, 
			map->mag, map->centre + position associator */
  controlConvertPerm(map); /* call all the conversion routines */
  
  /* make the initial set of columns */
  controlMakeColumns(map, defaultColumns, view, nmaps>1);
  
  /* A button which don't really belong to any column,
     so we tack them on the end here */
  
  controlAddButton(map, highlightButton, clearMenu);
  
  if (getCentre(from, map->key,  &centre))
    { map->centre = centre ;
      mapControlCursorSet (map, centre) ;
    }
  else
    if (control->hideHeader)
      /* if hideheader set in the view default to displaying it all */
      { BOOL isNeg;
	map->centre = (map->max + map->min) / 2.0 ;  
	isNeg = (map->mag < 0) ;
	map->mag = (map->control->graphHeight - 2.0) /
    			(map->max - map->min) ;
  
	if (isNeg)
	  map->mag = -map->mag ;
      }

  return map ; 
}

/* display newMap in a new MAP on the current control, if not
   shown already, and centre and size correctly
*/

BOOL gMapFollowMap (COLCONTROL control, MAPCONTROL oldMap, KEY newMapKey)
{ 
  float x, dx, x1, x2 ;
  float y1, y2 ;
  int i ;
  MAPCONTROL newMap = 0 ;

				/* do we have it already? */
  for (i=0; i < arrayMax(control->maps); i++)
    { newMap = arr(control->maps, i, MAPCONTROL);
      if (newMap->key == newMapKey)
	break ;
    }
				/* if not, make a new map */
  if (i == arrayMax(control->maps))
    if (!(newMap = gMapDoMap (control, newMapKey, 0, 0, 2)))
      return FALSE ;		/* 2 forces minimal-view */
  if (!newMap)
    return FALSE ;
				/* get ends on current map */
  if ((gMapGetMapObject (newMapKey, oldMap->key, 0, 20, &x, &dx, 0, 0) 
      		!= FLAG_ANY_INTERVAL) || !dx)
    return FALSE ;
  x1 = x - dx ; x2 = x + dx ;
				/* extent in its own coords */
  y1 = newMap->min ; y2 = newMap->max ;
				/* new centre, mag */
  newMap->centre = y1 + (oldMap->centre - x1) * (y2 - y1) / (x2 - x1) ;
  newMap->mag = oldMap->mag * (x2 - x1) / (y2 - y1) ;
  
  controlDraw () ;
  return TRUE ;
}

static Array allWindows = 0;

static void gMapControlRemove(void *cp)
{ int i;
  COLCONTROL c = *((COLCONTROL *)cp);
  for ( i=0; i<arrayMax(allWindows); i++)
    if (arr(allWindows, i, COLCONTROL) == c)
      arr(allWindows, i, COLCONTROL) = 0;
}

KEY gMapSelectView (KEY v)
{
  static KEY view = 0 ;
  if (class(v)) { view = v ; return 0 ; }
  else { v = view ; view = 0 ; return v ; }
}

BOOL gMapDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  COLCONTROL control;
  int i;
  KEYSET maps = keySetCreate();
  Associator views = assCreate();
  BSMARK mark = 0;
  KEY _VManyMap, _VView, mm;
  BOOL haveControl = FALSE ;

  lexaddkey("ManyMap", &mm, _VMainClasses);
  _VManyMap = KEYKEY(mm);
  lexaddkey("View", &mm, _VMainClasses);
  _VView = KEYKEY(mm);

/* The first stage, generate a set of maps to display, object to highlight */

  if (key && class(key) == _VMap)
   { keySet(maps,keySetMax(maps)) = key; 
     if (from && class(from) == _VView) 
       /* can display from a view - use that */
       assInsert(views, assVoid(key), assVoid(from));
   }
  /* display something from a map, use that */
  else if (from && class(from) == _VMap)
    { keySet(maps,keySetMax(maps)) = from;
      from = key;
    }
  
  /* displaying a mapped object, get the maps it's on. */
  /* this code gets multimaps, manymaps also */
  else
    { KEY k;
      OBJ obj = bsCreate(key);
      if (obj)
	{ if (bsFindTag(obj, str2tag("Display_map")))
	    { if (bsGetKey(obj, str2tag("Display_map"), &k))
		/* display_maps over-rides Maps */
		do {
		  keySet(maps,keySetMax(maps)) = k;
		} while (bsGetKey(obj, _bsDown, &k));
	    }
	  else if (bsGetKey(obj, _Map, &k))
	    do {
	      keySet(maps,keySetMax(maps)) = k;
	      if (class(key) == _VManyMap)
		{ KEY v;
		  mark = bsMark(obj, mark);
		  if (bsGetKey(obj, _bsRight, &v))
		    assInsert(views, assVoid(k), assVoid(v));
		  bsGoto(obj, mark);
		}
	      break ; /*  mieg: many maps must come explicitelly, not by defaults */
	    } while (bsGetKey(obj, _bsDown, &k));
	  bsDestroy (obj) ;
	  bsMarkFree(mark);
	}
      from = key;
    }

  if (keySetMax(maps) == 0 && class (key) == _VSequence) /* infer the map */
    { 
      KEY seq = key, clone = 0 ;
      OBJ obj = 0 ;
      
      /* seach a clone */
      while (seq && !clone)
	{
	  if ((obj = bsCreate (seq)))
	    {
	      seq = 0 ;
	      if (!bsGetKey (obj, _Clone, &clone))
		bsGetKey (obj, _Source, &seq) ;
	      bsDestroy (obj) ;
	    }
	}
      if (clone)
	{
	  KEY k ;
	  float x ;

	  if (gMapPhysClone2map (clone, 0, &k, &x))
	    {
	      keySet(maps,keySetMax(maps)) = k;
	      from = KEYMAKE (_VCalcul, (int)((x + 1000) * 1000)) ;
	    }
	}
    }

  if (keySetMax(maps) == 0) /* failed to find any maps */
    { if (from) 
	display (from, 0, TREE) ;
      keySetDestroy(maps);
      assDestroy(views);
      return FALSE ;
    }

  /* Remove any non-graphic maps from the list and display them as trees */
  keySetSort(maps); /* for KeySetFind */
  for (i=0; i<arrayMax(maps); i++)
    { KEY map = arr(maps, i, KEY);
      OBJ obj = bsCreate(map);
      if (obj && bsFindTag(obj, _Non_graphic))
	{ keySetRemove(maps, map);
	  display(map, from, TREE);
	}
      bsDestroy(obj);
    }
  if (keySetMax(maps) == 0)
    { keySetDestroy(maps);
      assDestroy(views);
      return FALSE;
    }

  if (!allWindows)
    allWindows = arrayCreate(10, COLCONTROL);

	/* rd changed the logic here: only match windows if isOldGraph */
  if (isOldGraph)
    { 
/* Stage two, see if we have any window's open with precisely the right
   maps already displayed, this is inspired by ace3. Note that if we have
   views tied down to any maps, we don't re-use, we re-draw */
  
      for (i=0 ; i<arrayMax(allWindows); i++)
	{ int j;
	  COLCONTROL c = arr(allWindows, i, COLCONTROL);
	  if (!c)
	    continue; /* deleted one */
	  if (arrayMax(c->maps) != keySetMax(maps))
	    continue; /* different number, failed */
	  for (j=0; j<arrayMax(c->maps); j++)
	    { MAPCONTROL m = arr(c->maps, j, MAPCONTROL);
	      if (!keySetFind(maps, m->key, 0))
		break;
	      if (assFind(views, assVoid(m->key), 0))
		break;
	    }
	  if (j == arrayMax(c->maps))
	    break; /* found the match */
	}
  
      if (i != arrayMax(allWindows)) /* found match */
	{ float centre;
	  control = arr(allWindows, i, COLCONTROL);
	  graphActivate(control->graph);
	  for (i=0; i<arrayMax(control->maps); i++)
	    { MAPCONTROL map = arr(control->maps, i, MAPCONTROL);
	      if (getCentre(from, map->key, &centre))
		{ map->centre = centre;
		  mapControlCursorSet(map, centre);
		}
	    }
	  haveControl = TRUE ;
	}
    }

  if (!haveControl)
/* Stage three make a new control, in new or existing window, then display
   the maps in it */
    { Stack buff = stackCreate (127) ;
      COLCONTROL *cp;

      for (i=0; i<arrayMax(maps); i++)
	{ KEY mk = arr(maps, i, KEY);
	  catText (buff, name(mk)) ;
	  catText (buff, " ") ;
	}

      control = colControlCreate(isOldGraph, stackText (buff,0), GMAP);
      stackDestroy (buff) ;
      array(allWindows, arrayMax(allWindows), COLCONTROL) = control;
      /* The next bit ensures that the pointer to the control in allWindows
	 is removed when it goes away */
      cp = (COLCONTROL *) handleAlloc(gMapControlRemove, 
				      control->handle,
				      sizeof(COLCONTROL));
      *cp = control;
      
      for (i=0; i<arrayMax(maps); i++)
	{ KEY mk = arr(maps, i, KEY), view;
	  void *vview;
	  if ((view = gMapSelectView (0)) /* mieg feb 5 select view for gif display */
	      ||( assFind(views, assVoid(mk), &vview) && (view = (KEY)assInt(vview))))
	    gMapDoMap(control, mk, from, view, arrayMax(maps));
	  else
	    gMapDoMap(control, mk, from, 0, arrayMax(maps));
	}
    }
  
  control->from = from;
  controlDraw();
  keySetDestroy(maps);
  assDestroy(views);
  return TRUE;
}

/*********************************************************************/
/* Some routines shared by the locus chrom_band and interval columns */
/*********************************************************************/	

void gMapUnselect(MAPCONTROL map)
/* this avoids every column having to know about all map-global
   highlighting data, even if it's in fact used only between a couple
   of columns. This function unsets everything that's shared in the map */
{
  GeneticMap look = (GeneticMap)map->look;

  /* look->selectKey is always valid */
  look->neighboursInfoValid = FALSE; /* look->neighbours valid when true */
  look->friendsInfoValid = FALSE;    /* look->positive, negative */
  look->dataInfoValid = FALSE;       /* look_multi, 2pt */
  look->coOrdsInvalid = FALSE;
  look-> negativeMarkersValid = FALSE; /* valid when submap contig picked */
}

BOOL gMapIsNeighbour(MAPCONTROL map, KEY key)
{
  GeneticMap look = (GeneticMap)map->look;
  
  if (key == map->control->activeKey)
    return TRUE;
  
  /* make sure our map is in charge */
  if (map->control->activeInstance && 
      map == map->control->activeInstance->map)
    { 
      if (look->neighboursInfoValid && keySetFind(look->neighbours, key, 0))
	return TRUE;
      if (key == look->selectKey)
	return TRUE;
    }
  
  return FALSE;
}

int gMapOverlap(MAPCONTROL map, KEY key, float z1, float z2)
{
  GeneticMap look = (GeneticMap)map->look;
  COLCONTROL control = map->control;
  int flags = 0;
  float x1 = look->x1;
  float x2 = look->x2;
  double dd = -12.60 ;
  float ff = -12.60, ee ;

/* ee is an estimate of rounding errors 
   when computing with float
*/
  ee = ff - dd ; if (ee < 0) ee = -ee ; ee *= 4 ;
  

  /* make sure our map is in charge */
  if (control->activeInstance && map == control->activeInstance->map)
    {  
      /* x1-x2 are the extremities of the selected seg who contains... */
      /* z1-z2 segs */
      
      if (look->friendsInfoValid && keySetFind(look->positive, key, 0))
	{ if (look->coOrdsInvalid)
	    flags |= FLAG_RELATED;
	  else if (x1 > z2 + ee || 
		   z1 - ee > x2)
	    flags |= FLAG_STRESSED;
	  else if ((x1 >= z1 - ee && x2 <= z2 + ee) ||
		   (z1 - ee >= x1 && z2 + ee <= x2))
	    flags |= FLAG_RELATED ;
	  else			/* staggered overlap */
	    flags |= FLAG_STRESSED ;
	}		  
      
      if (look->friendsInfoValid && keySetFind(look->negative, key, 0))
	{ if (look->coOrdsInvalid)
	    flags |= FLAG_ANTI_RELATED;
	  else if (x1 >= z2 + ee || 
		   z1 - ee >= x2)
	    flags |= FLAG_ANTI_RELATED ;
	  else if ((x1 > z1 - ee && x2 < z2 + ee) ||
		   (z1 - ee  > x1 && z2 + ee < x2))
	    flags |= FLAG_ANTI_STRESSED ;
	  else 			/* staggered overlap */
	    flags |= FLAG_ANTI_STRESSED ;
	}		  
    }

  return flags;
}

/* NB this is used by the markers column also */
int gMapOrder (const void *a, const void *b)  /* for arraySort() call */
{
  const GMAPSEG *seg1 = (const GMAPSEG*)a, *seg2 = (const GMAPSEG*)b ;
  float diff, x1, x2;

  if (seg1->flag & FLAG_ANY_INTERVAL)
    x1 = seg1->x - seg1->dx; /* sort on start posn of intervals */
  else
    x1 = seg1->x;

  if (seg2->flag & FLAG_ANY_INTERVAL)
    x2 = seg2->x - seg2->dx;
  else
    x2 = seg2->x;
  
  diff = x1 - x2 ;

  if (diff > 0)
    return 1 ;
  else if (diff < 0)
    return -1 ;
  else
    return seg2->key - seg1->key ; /* never give a random result please */
}


/* Convert Routine for locus, interval, and Chrom_band */
void *gMapConvert (MAPCONTROL map, void *params)
{
#define NPOS 1000
  GeneticMap look = (GeneticMap)map->look;
  KEY chrom, clone, locus, cache;
  OBJ Locus;
  Array segs, loci;
  GMAPSEG *seg;
  float x[NPOS], dx[NPOS];
  int i, j, flags, count;


  /* are they cached */
  if ((cache = gMapCacheKey(map->key, "main")) &&
      (segs = arrayGet(cache, GMAPSEG, segFormat)))
    { look->segs = segs;
      return segs;
    }
  segs = arrayCreate (128,GMAPSEG) ;
  /* This array is explicitly destroyed in gMapFinalise */
  
  loci = gMapMapObjects(map->key) ;

  for (i = 0; i<arrayMax(loci); i++)
    { locus = arr(loci, i, KEYPAIR).obj;
      chrom = arr(loci, i, KEYPAIR).map;
      count = NPOS;
      flags = gMapGetMapObject(locus, map->key, chrom, 
			       20, x, dx, &Locus, &count);
      if (!flags)
	continue;
	  
      if (bsFindTag (Locus, _Well_ordered))
	flags |= FLAG_WELL_ORDERED ;

      if (bsFindTag (Locus, _Centromere))
	flags |= FLAG_CENTROMERE ;
      
      if (!bsFindTag (Locus, _Compound))
	{ if (bsFindTag (Locus, _Duplication))
	    flags |= FLAG_DUPLICATION ;
	  if (bsFindTag (Locus, _Deletion) ||
	      bsFindTag (Locus, str2tag("Deficiency")))
	    flags |= FLAG_DEFICIENCY ;
	}
      
      if (bsGetKey (Locus,_Positive_clone,&clone))
	flags |= FLAG_CLONED ;
      
      
      for (j=0; j<count; j++)
	{ seg = arrayp(segs,arrayMax(segs),GMAPSEG) ;
	  
	  seg->key = locus;
	  seg->x = x[j] ;
	  seg->dx = dx[j] ;
	  seg->flag = flags;
	}
      
      bsDestroy (Locus) ;
    }
  
  arrayDestroy(loci);
  
  arraySort (segs, gMapOrder) ; 

  gMapCacheStore(map->key, "main", segs, segFormat);
  look->segs = segs ; /* Make this display list commonly available. */
  return (void *)segs ;
}

void gMapSaveDetails (Array segs, KEY map)
{
  OBJ obj ;
  int i ;
  GMAPSEG *seg ;
  float x ;
  int nLocus = 0, nInterval = 0 ;

  /* destroy a cached set of segs */
  gMapCacheKill(map);

  for (i = 0 ; i < arrayMax(segs) ; ++i)
    { seg = arrp(segs,i,GMAPSEG) ;
      if (seg->flag & FLAG_MOVED)
	{
	  if (seg->flag & FLAG_ANY_LOCUS)
	    {
	      if ((obj = bsUpdate (seg->key)))
		{
		  bsAddKey (obj, _Map, map) ;
		  bsPushObj (obj) ;
		  bsAddData (obj, _Position, _Float, &seg->x) ;
		  bsPushObj (obj) ;
		  bsAddData (obj, _Error, _Float, &seg->dx) ;
		  bsSave (obj) ;
		  ++nLocus ;
		}
	    }
	  else if (seg->flag & FLAG_ANY_INTERVAL)
	    {
	      if ((obj = bsUpdate (seg->key)))
		{
		  bsAddKey (obj, _Map, map) ;
		  bsPushObj (obj) ;
		  x = seg->x - seg->dx;
		  bsAddData (obj, _Left, _Float, &x) ;
		  x = seg->x + seg->dx ;
		  bsAddData (obj, _Right, _Float, &x) ;
		  bsSave (obj) ;
		  ++nInterval ;
		}
	    }
	  else
	    messout ("%s has moved, but I don't know how to save it",
		     name(seg->key)) ;
	}
      seg->flag &= ~FLAG_MOVED ;
    }

  messout ("Saved %d loci, %d intervals", nLocus, nInterval) ;
}

/************* end of file ****************/



 
 
 
 
 
 
