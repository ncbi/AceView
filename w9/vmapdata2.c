/*  File: vmapdata2.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: basic routines for Muli_pt and 2pt data
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  7 17:47 1999 (fw)
 * Created: Tue Nov 30 18:42:30 1993 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: vmapdata2.c,v 1.5 2015/08/12 14:12:03 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "sysclass.h"
#include "classes.h"
#include "systags.h"
#include "tags.h"
#include "query.h"
#include "lex.h"
#include "a.h"
#include "vmap_.h"

/************************************************************/

#define myMapDraw(_look, _box) vMapDraw(_look, _box)

static BOOL getMulti (KEY) ;
static SEG* insertMulti (VerticalMap look, KEY key) ;
static BOOL logLikeMulti (float *ll) ;
static void multiBoundCheck (KEY locus, float *min, float *max, 
			     KEY *minLoc, KEY *maxLoc) ;
static BOOL get2p (KEY key, KEY *loc1, KEY *loc2, KEY *type) ;
static SEG* insert2p (VerticalMap look, KEY key) ;
static float logLike2pt (float dist, KEY type) ;
static void calcAll2pt (void) ;
static BOOL getPos (KEY key, float *y) ;
static void setPos (VerticalMap look) ;
static void getAllData (VerticalMap look) ;
static void clearMapData (void) ;

static KEYSET badData = 0 ;

KEY dataKey ;
KEYSET newBadData = 0 ; 
Array loci = 0 ;
Array counts = 0 ;

/**************************************************/

static KEYSET orDestroy (KEYSET a, KEYSET b) /* utility */
{
  KEYSET c = keySetOR (a, b) ;
  keySetDestroy (a) ;
  keySetDestroy (b) ;
  return c ;
}

/************************************/

static void checkDataStruct (VerticalMap look)
{
  if (!look->data)
    { look->data = (vMapData) messalloc (sizeof(struct vMapDataStruct)) ;
      look->data->loc2two = assCreate () ;
      look->data->loc2multi = assCreate () ;
      look->data->loc2inside = assCreate () ;
      look->data->loc2outside = assCreate () ;
    }
}

static void getSegData (VerticalMap look, SEG *seg)
{
  KEYSET a, b, c ;

  if (seg->flag & FLAG_HAVE_DATA)
    return ;
  seg->flag |= FLAG_HAVE_DATA ;

  if (!newBadData)		/* must be done somewhere */
    newBadData = keySetCreate () ;

  checkDataStruct (look) ;

  messStatus ("Getting map data") ;

  if (class(seg->key) == _VLocus)
    { a = queryKey (seg->key, ">2_point") ;
      if (keySetMax(a))
	assInsert (look->data->loc2two, assVoid(seg->key), a) ;
      else
	keySetDestroy (a) ;
      a = queryKey (seg->key, ">Multi_point") ;
      if (keySetMax(a))
	assInsert (look->data->loc2multi, assVoid(seg->key), a) ;
      else
	keySetDestroy (a) ;
      a = queryKey (seg->key, ">Df_Dup ; Inside ; >Rearrangement") ;
      b = queryKey (seg->key, ">Inside") ;
      c = orDestroy (a, b) ;
      if (keySetMax(c))
	assInsert (look->data->loc2inside, assVoid(seg->key), c) ;
      else
	keySetDestroy (c) ;
      a = queryKey (seg->key, ">Df_Dup ; Outside ; >Rearrangement") ;
      b = queryKey (seg->key, ">Outside") ;
      c = orDestroy (a, b) ;
      if (keySetMax(c))
	assInsert (look->data->loc2outside, assVoid(seg->key), c) ;
      else
	keySetDestroy (c) ;
    }
  else if (class(seg->key) == _VInterval)
    { a = queryKey (seg->key, ">Df_Dup ; Inside ; >Locus") ;
      b = queryKey (seg->key, ">Inside") ;
      c = orDestroy (a, b) ;
      if (keySetMax(c))
	assInsert (look->data->loc2inside, assVoid(seg->key), c) ;
      else
	keySetDestroy (c) ;
      a = queryKey (seg->key, ">Df_Dup ; Outside ; >Locus") ;
      b = queryKey (seg->key, ">Outside") ;
      c = orDestroy (a, b) ;
      if (keySetMax(c))
	assInsert (look->data->loc2outside, assVoid(seg->key), c) ;
      else
	keySetDestroy (c) ;
    }
}

static void getAllData (VerticalMap look)
{
  int i ;

  if (look->flag & FLAG_ALL_DATA)
    return ;
  look->flag |= FLAG_ALL_DATA ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    getSegData (look, arrp(look->segs, i, SEG)) ;
}

static void regetData (void)
{
  int i ;
  VerticalMap look = currentVerticalMap("regetData") ;

  look->flag &= ~FLAG_ALL_DATA ;
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    arrp(look->segs,i,SEG)->flag &= ~FLAG_HAVE_DATA ;

  getAllData (look) ;
}

/*********************************************/

void vMapDataDestroy (VerticalMap look)
{
  KEYSET kset ;
  void *k ;
  vMapData data = look->data ;

  if (!data)
    return ;

  for (k = 0 ; assNext (data->loc2two, &k, &kset) ; )
    keySetDestroy (kset) ;
  for (k = 0 ; assNext (data->loc2multi, &k, &kset) ; )
    keySetDestroy (kset) ;
  for (k = 0 ; assNext (data->loc2inside, &k, &kset) ; )
    keySetDestroy (kset) ;
  for (k = 0 ; assNext (data->loc2outside, &k, &kset) ; )
    keySetDestroy (kset) ;
  keySetDestroy (data->orderedLoci) ;

  messfree (data) ;
}

/**********************************/ /* fast position lookup */

static Array posArray = 0 ;
static KEY locusBase ; 

static BOOL getPos (KEY key, float *y)
{ 
  key -= locusBase ;
  if (key > arrayMax(posArray) ||
      arr(posArray, key, float) < -999999)
    return FALSE ;
  *y = arr(posArray, key, float) ;
  return TRUE ;
}

static void setPos (VerticalMap look)
{ 
  int i ;
  KEY max ;
  SEG *seg ;
      
  posArray = arrayReCreate (posArray, 1024, float) ;
  max = locusBase = KEYMAKE(_VLocus, 0) ;

  checkDataStruct (look) ;
  look->data->orderedLoci = keySetReCreate (look->data->orderedLoci) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (class(seg->key) == _VLocus && seg->key > max)
	max = seg->key ;
    }
  max -= locusBase ;

  for (i = 0 ; i < max ; ++i)
    array (posArray, i, float) = -1000000 ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)  
    { seg = arrp(look->segs,i,SEG) ;
      if (class(seg->key) == _VLocus)
	array (posArray, seg->key - locusBase, float) = seg->x ;
      if (seg->flag & FLAG_CLONED &&
	  !(seg->flag & FLAG_PHYS_GENE))
	keySetInsert (look->data->orderedLoci, seg->key) ;
    }
}

static void setTestPos (KEY key, float pos)
{
  key -= locusBase ;
  if (key > arrayMax(posArray))
    { int i ;
      for (i = arrayMax(posArray) ; i < key ; ++i)
	array(posArray, i, float) = -1000000 ;
    }
  array(posArray, key, float) = pos ;
}

/*********************************************************/
/****************** display routines *********************/

static void contig2data (VerticalMap look, SEG *seg)
{
  int i, x ; 
  KEYSET  genes ;
  KEY map ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    arrp(look->segs,i,SEG)->flag &= 
      ~(FLAG_STRESSED | FLAG_RELATED | FLAG_ANTI_RELATED) ;

  if (!lexReClass (look->key, &map, _VMap))
    return ;
  genes = queryKey (map, 
     messprintf(">Locus ; >Positive_clone pMap = %s ; "
		">Positive_locus", name(seg->key))) ;
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    if (keySetFind (genes, arrp(look->segs,i,SEG)->key, &x))
      arrp(look->segs,i,SEG)->flag |= FLAG_RELATED ;
  keySetDestroy(genes) ;

  myMapDraw (look, 0) ;
}

/********************/

static void rearrangement2data (VerticalMap look, SEG *seg)
{
  KEYSET genes ;
  int i, x ; 
  void *vs, *va ;
  SEG *seg1 ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    arrp(look->segs,i,SEG)->flag &= 
      ~(FLAG_STRESSED | FLAG_RELATED | 
	FLAG_ANTI_STRESSED | FLAG_ANTI_RELATED) ;

  getSegData (look, seg) ;

  vs = assVoid(seg->key) ; 
  if (assFind (look->data->loc2inside, vs, &va))
    { genes = (KEYSET) va ;
      for (i = 0 ; i < arrayMax(look->segs) ; ++i)
	if (keySetFind (genes, arrp(look->segs,i,SEG)->key, &x))
	  { seg1 = arrp(look->segs,i,SEG) ;
	    if (seg1->x >= seg->x && seg1->x <= seg->x + seg->dx)
	      seg1->flag |= FLAG_RELATED ;
	    else
	      seg1->flag |= FLAG_STRESSED ; 
	  }
    }

  if (assFind (look->data->loc2outside, vs, &va))
    { genes = (KEYSET) va ;
      for (i = 0 ; i < arrayMax(look->segs) ; ++i)
	if (keySetFind (genes, arrp(look->segs,i,SEG)->key, &x))
	  { seg1 = arrp(look->segs,i,SEG) ;
	    if (seg1->x > seg->x && seg1->x < seg->x + seg->dx)
	      seg1->flag |= FLAG_ANTI_STRESSED ;
	    else
	      seg1->flag |= FLAG_ANTI_RELATED ;
	  }
    }

  myMapDraw(look, 0) ;
}

/****************/

static void locus2data (VerticalMap look, SEG *seg)
{
  int i, x ;
  float dummy ;
  void *va, *vs = assVoid(seg->key) ;
  SEG *seg1, *seg0 = seg ;
  KEYSET data, intervals ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    arrp(look->segs,i,SEG)->flag &= 
      ~(FLAG_STRESSED | FLAG_RELATED | 
	FLAG_ANTI_STRESSED | FLAG_ANTI_RELATED) ;

  getSegData (look, seg) ;

  if (assFind (look->data->loc2two, vs, &va))
    { data = (KEYSET) va ;
      for (i = 0 ; i < keySetMax(data) ; ++i)
	if ((seg = insert2p (look, keySet(data, i))))
	  seg->flag |= FLAG_RELATED ;
    }

  if (assFind (look->data->loc2multi, vs, &va))
    { data = (KEYSET) va ;
      for (i = 0 ; i < keySetMax(data) ; ++i)
	if ((seg = insertMulti (look, keySet(data, i))))
	  {
	    if (logLikeMulti (&dummy))
	      seg->flag |= FLAG_RELATED ;
	    else
	      seg->flag |= FLAG_STRESSED ;
	  }
    }

  if (assFind (look->data->loc2inside, vs, &va))
    { intervals = (KEYSET) va ;
      for (i = 0 ; i < arrayMax(look->segs) ; ++i)
	if (keySetFind (intervals, arrp(look->segs,i,SEG)->key, &x))
	  { seg1 = arrp(look->segs,i,SEG) ;
	    if (seg0->x >= seg1->x && seg0->x <= seg1->x + seg1->dx)
	      seg1->flag |= FLAG_RELATED ;
	    else
	      seg1->flag |= FLAG_STRESSED ;
	  }
    }

  if (assFind (look->data->loc2outside, vs, &va))
    { intervals = (KEYSET) va ;
      for (i = 0 ; i < arrayMax(look->segs) ; ++i)
	if (keySetFind (intervals, arrp(look->segs,i,SEG)->key, &x))
	  { seg1 = arrp(look->segs,i,SEG) ;
	    if (seg0->x > seg1->x && seg0->x < seg1->x + seg1->dx)
	      seg1->flag |= FLAG_ANTI_STRESSED ;
	    else
	      seg1->flag |= FLAG_ANTI_RELATED ;
	  }
    }

  myMapDraw (look, 0) ;
}

/********************/

static void seg2data (VerticalMap look, SEG *seg)
{
  if (class(seg->key) == _VContig)
    contig2data (look, seg) ;
  else if (class(seg->key) == _VInterval)
    rearrangement2data (look, seg) ;
  else if (class(seg->key) == _VLocus)
    locus2data (look, seg) ;
  else
    messout ("First pick a gene or rearrangment") ;
}

void vMapShowDataFromMenu (int box)
{
  int i ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("ShowLocusData") ;

  if ((i = arr(look->boxIndex, box, int)) &&
      (seg = arrp(look->segs, i, SEG)))
    { if (box != look->activeBox)
	vMapSelect (look, box) ;
      seg2data (look, seg) ;
    }
}

void vMapShowData (void)
{
  int i ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("gMapShowData") ;

  clearMapData () ;
  if ((i = arr(look->boxIndex, look->activeBox, int)) &&
      (seg = arrp(look->segs, i, SEG)))
    seg2data (look, seg) ;
}

/************************************************/
/************ likelihood stuff ******************/

static float logLikeLocus (KEY locus)
{
  KEYSET data ;
  float ll = 0, bit, pos, y1, y2 ;
  KEY loc1, loc2, type ;
  int i, ip ;
  VerticalMap look = currentVerticalMap("logLikeLocus") ;

  if (!getPos (locus, &pos))
    return ll ;
  
  if (assFind (look->data->loc2multi, assVoid(locus), &data))
    for (i = keySetMax(data) ; i-- ;)
      if (!(badData && keySetFind (badData, keySet(data, i), &ip)) &&
	  getMulti (keySet(data, i)) &&
	  logLikeMulti (&bit))
	ll += bit ;

  if (assFind (look->data->loc2two, assVoid(locus), &data))
    for (i = keySetMax(data) ; i-- ;)
      if (!(badData && keySetFind (badData, keySet(data, i), &ip)) &&
	  get2p (keySet(data, i), &loc1, &loc2, &type) &&
	  getPos (loc1, &y1) && getPos (loc2, &y2))
	{
	  if (y2 > y1)
	    ll += logLike2pt (y2-y1, type) ;
	  else
	    ll += logLike2pt (y1-y2, type) ;
	}
  return ll ;
}

static BOOL boundFind (KEY locus, float *min, float *max,
		       KEY *minLoc, KEY *maxLoc)
{
  KEYSET in, out ;
  float x, y, lx ;
  int i, j, ip ;
  VerticalMap look = currentVerticalMap("boundFind") ;
  
  *min = -1000000 ; *max = 1000000 ;
  *minLoc = 0 ; *maxLoc = 0 ;
  if (!getPos (locus, &x))
    return FALSE ;

  if (assFind (look->data->loc2inside, assVoid(locus), &in))
    for (i = keySetMax(in) ; i-- ;)
      if (!(badData && keySetFind (badData, keySet(in, i), &ip)) &&
	  assFind (look->data->loc2outside, assVoid(keySet(in,i)), &out))
	for (j = keySetMax(out) ; j-- ;)
	  if (getPos (keySet(out,j), &y))
	    {
	      if (y < x && y > *min)
		{ *min = y ; *minLoc = keySet(out, j) ; }
	      else if (y > x && y < *max)
		{ *max = y ; *maxLoc = keySet(out, j) ; }
	    }

  if (assFind (look->data->loc2outside, assVoid(locus), &out))
    for (i = keySetMax(out) ; i-- ;)
      if (!(badData && keySetFind (badData, keySet(out, i), &ip)) &&
	  assFind (look->data->loc2inside, assVoid(keySet(out,i)), &in))
	for (j = keySetMax(in) ; j-- ;)
	  if (getPos (keySet(in,j), &y))
	    {
	      if (y < x && y > *min)
		{ *min = y ; *minLoc = keySet(in, j) ; }
	      else if (y > x && y < *max)
		{ *max = y ; *maxLoc = keySet(in, j) ; }
	    }

  if (assFind (look->data->loc2multi, assVoid(locus), &in))
    for (i = keySetMax(in) ; i-- ;)
      if (!(badData && keySetFind (badData, keySet(in, i), &ip)) &&
	  getMulti (keySet(in, i)))
	multiBoundCheck (locus, min, max, minLoc, maxLoc) ;

  if (keySetFind (look->data->orderedLoci, locus, &ip))
    for (i = keySetMax(look->data->orderedLoci) ; i-- ;)
      if (keySet(look->data->orderedLoci,i) != locus &&
	  getPos(keySet(look->data->orderedLoci,i), &y))
	{
	  if (y < x && y > *min)
	    { *min = y ; *minLoc = keySet(look->data->orderedLoci,i) ; }
	  else if (y > x && y < *max)
	    { *max = y ; *maxLoc = keySet(look->data->orderedLoci,i) ; }
	}
  
  if (*min < -999999)
    { lx = logLikeLocus (locus) - 3 ;
      for (y = 0.1 ; y < 10 ; y *=2)
	{ setTestPos (locus, x-y) ;
	  if (logLikeLocus (locus) < lx)
	    break ;
	}
      if (logLikeLocus (locus))
	*min = x-y ;
      setTestPos (locus, x) ;
    }

  if (*max > 999999)
    { lx = logLikeLocus (locus) - 3 ;
      for (y = 0.1 ; y < 10 ; y *=2)
	{ setTestPos (locus, x+y) ;
	  if (logLikeLocus (locus) < lx)
	    break ;
	}
      if (logLikeLocus (locus))
	*max = x+y ;
      setTestPos (locus, x) ;
    }

  return (*min > -999999 && *max < 999999) ;
}

/**********************************************/

#define NBIN 50
float dbn[NBIN], dbnMin, dbnMax, dMax ; 
SEG *vDbnSeg ;
#define DBN_X(z) ((((z)+1)*dbnMax + (NBIN-(z))*dbnMin)/(NBIN+1))

static BOOL makeDbn (KEY locus)
{
  int i ;
  float oldx ;

  dMax = -1E20 ; ;

  if (!getPos (locus, &oldx))
    return FALSE ;

  for (i = 0 ; i < NBIN ; ++i)
    { setTestPos (locus, DBN_X(i)) ;
      dbn[i] = logLikeLocus (locus) ;
      if (dbn[i] > dMax)
	dMax = dbn[i] ;
    }
  for (i = 0 ; i < NBIN ; ++i)
    dbn[i] = exp(dbn[i] - dMax) ;

  setTestPos (locus, oldx) ;
  return TRUE ;
}

static void getBadData (void)
{
  int i, t ;
  KEY key ;
  KEYSET active ;

  if (keySetExists (badData))
    keySetDestroy (badData) ;
  
  if (lexword2key ("bad_map_data", &key, _VKeySet))
    badData = arrayGet (key, KEY, "k") ;

  if (keySetActive (&active, 0))
    { if (!badData)
	badData = keySetCreate () ;
      for (i = 0 ; i < keySetMax(active) ; ++i)
	{ t = class(arr(active,i,KEY)) ;
	  if (t == _V2_point_data ||
	      t == _VMulti_pt_data ||
	      t == _VInterval)
	    keySetInsert (badData, arr(active,i,KEY)) ;
	}
      if (!keySetMax (badData))
	keySetDestroy (badData) ;
    }
}

void vMapDrawDbn (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  int i ;
  float x ;
  KEY minLoc, maxLoc ;

  if (!vDbnSeg || !getPos (vDbnSeg->key, &x))
   return ;

  getSegData (look, vDbnSeg) ;
  getBadData () ;
  if (!boundFind (vDbnSeg->key, &dbnMin, &dbnMax, &minLoc, &maxLoc))
    { if (dbnMin < -999999)
	dbnMin = x - 10 ;
      if (dbnMax > 999999)
	dbnMax = x + 10 ;
    }
  makeDbn (vDbnSeg->key) ;

  graphColor (RED) ;
  x = MAP2GRAPH(look->map, dbnMin) ;
  graphLine (*offset, x, *offset+6, x) ;
  x = MAP2GRAPH(look->map, dbnMax) ;
  graphLine (*offset, x, *offset+6, x) ;

  graphColor (BLUE) ;
  for (i = 0 ; i < NBIN ; ++i)
    { x = DBN_X(i) ;
      x = MAP2GRAPH(look->map, x) ;
      graphLine (*offset, x, *offset + 6*dbn[i], x) ;
    }
  graphColor (BLACK) ;
  *offset += 6 ;
}

void vMapMakeDbn (int box)
{
  int i ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("makeDbn") ;

  if (!(i = arr(look->boxIndex, box, int)) ||
      !(seg = arrp(look->segs, i, SEG)))
    return ;
  if (box != look->activeBox)
    vMapSelect (look, box) ;
  vDbnSeg = seg ;
  mapColSetByName ("Likelihood", TRUE) ;
  myMapDraw (look, seg->key) ;
  return ;
}

/*****************************************************/

static void gjmOptAll (void)
{
  int i, j, ibestmin, n = 0, ibestmax ;
  SEG *seg ;
  float x, best, gain = 0, rms = 0, diff, worst = 0 ;
  KEY worstKey=0, minLoc, maxLoc ;
  VerticalMap look = currentVerticalMap("gjmOptAll") ;

  getAllData (look) ;
  getBadData () ;
  newBadData = keySetReCreate (newBadData) ;

  for (j = 0 ; j < arrayMax(look->segs) ; ++j)
    { seg = arrp(look->segs,j,SEG) ;
      if (class(seg->key) != _VLocus ||
	  seg->flag & FLAG_HIGHLIGHT ||
	  seg->flag & FLAG_PHYS_GENE ||
	  !boundFind (seg->key, &dbnMin, &dbnMax, &minLoc, &maxLoc) ||
	  !makeDbn (seg->key))
	continue ;

      best = dbn[0] ; ibestmin = ibestmax = 0 ;
      for (i = 1 ; i < NBIN ; ++i)
	{ if (dbn[i] >= best)
	    { if (dbn[i] > best)
		{ best = dbn[i] ; ibestmin = i ; }
	      ibestmax = i ;
	    }
	}

      best = 0.5*(ibestmin+ibestmax) ; /* for if flat dbn */
      x = DBN_X(best) ;
      diff = (x - seg->x) ; if (diff < 0) diff = -diff ;
      rms += diff*diff ;
      gain += dMax - logLikeLocus(seg->key) ;
      if (diff > worst)
	{ worst = diff ;
	  worstKey = seg->key ;
	}
      ++n ;

      if (x != seg->x)
	{ seg->x = x ;
	  seg->flag |= FLAG_MOVED ;
	}
      setTestPos (seg->key, x) ; /* reset test position */

      for (i = ibestmin ; i > 0 ; --i)
	if (dbn[i] < 0.1)
	  break ;
      dbnMin = DBN_X(i) ;
      for (i = ibestmax ; i < NBIN-1 ; ++i)
	if (dbn[i] < 0.1)
	  break ;
      dbnMax = DBN_X(i) ;
      diff = 0.5 * (dbnMax - dbnMin) ;
      if (diff != seg->dx)
	{ seg->dx = diff ;
	  seg->flag |= FLAG_MOVED ;
	}
    }

  if (n)
    messout ("%d loci\n"
	     "RMS moverment %.3f\n"
	     "Average gain %.2f\n"
	     "%s moved %f\n", 
	     n, sqrt(rms/n), gain/n, name(worstKey), worst) ;
  myMapDraw (look, 0) ;

  if (keySetMax (newBadData))
    { displayCreate (DtKeySet) ;	/* new window */
      keySetShow (newBadData, 0) ;
      newBadData = keySetCreate () ; /* make a new one for here */
    }
}

/**************/

static void clearMapData (void)
{
  VerticalMap look = currentVerticalMap("ClearMapData") ;

  arrayMax(look->segs) = look->lastTrueSeg ;
  if (look->data)
    look->data->key2seg = assReCreate(look->data->key2seg) ;
  mapColSetByName ("Likelihood", FALSE) ;
  myMapDraw (look,0) ;
}

static void highlightDataless (void)
{
  int i ;
  SEG *seg ;
  KEY minLoc, maxLoc ;
  VerticalMap look = currentVerticalMap("highlightDataless") ;

  getAllData (look) ;
  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (class(seg->key) == _VLocus &&
	  !boundFind (seg->key, &dbnMin, &dbnMax, &minLoc, &maxLoc))
	seg->flag |= FLAG_HIGHLIGHT ;
    }
}

static void check2point (void)
{
  KEY key = 0 ;
  KEY loc1, loc2, type ;

  while (lexNext (_V2_point_data, &key))
    if (get2p (key, &loc1, &loc2, &type) &&
	!logLike2pt (1, type))
      fprintf (stderr,"Bad 2_point %s\n", name (key)) ;
}

static void positionIntervals (void)
{
  void* v ;
  KEY rearr ;
  Array in, out ;
  int i ;
  float minIn, maxIn, maxBelow, minAbove, x, x1, x2 ;
  SEG *seg ;
  VerticalMap look = currentVerticalMap("positionRearrangements") ;

  getBadData () ;
  newBadData = keySetReCreate (newBadData) ;
  getAllData (look) ;

  v = 0 ;
  while (assNext (look->data->loc2inside, &v, &in))
    if (class(rearr = (KEY)assInt(v)) == _VInterval)
      { for (i = 0 ; i < arrayMax(look->segs) ; ++i)
	  if (arrp(look->segs, i, SEG)->key == rearr)
	    break ;
	if (i >= arrayMax(look->segs))
	  continue ;
	seg = arrayp(look->segs,i,SEG) ; /* new if necessary */
	
	minIn = 1000000 ;
	maxIn = -1000000 ;
	for (i = keySetMax(in) ; i-- ;)
	  if (getPos (keySet(in,i), &x))
	    { if (x < minIn) minIn = x ;
	      if (x > maxIn) maxIn = x ;
	    }
	if (minIn > 999999)
	  continue ;

	maxBelow = -1000000 ;
	minAbove = 1000000 ;
	if (assFind (look->data->loc2outside, assVoid(rearr), &out))
	  for (i = keySetMax(out) ; i-- ;)
	    if (getPos (keySet(out,i), &x))
	      { if (x < minIn)
		  { if (x > maxBelow) maxBelow = x ; }
	      else if (x > maxIn) 
		{ if (x < minAbove) minAbove = x ; }
	      else
		keySetInsert (newBadData, rearr) ;
	      }

	if (maxBelow < -999999) 
	  x1 = minIn - 0.01 ;
	else 
	  x1 = 0.5 * (minIn + maxBelow) ;
	if (minAbove > 999999) 
	  x2 = maxIn + 0.01 ;
	else 
	  x2 = 0.5 * (maxIn + minAbove) ;
	seg->key = rearr ;
	if (seg->x != x1 || seg->dx != x2-x1)
	  { seg->x = x1 ;
	    seg->dx = x2 - x1 ;
	    seg->flag |= FLAG_MOVED ;
	  }
      }

  myMapDraw (look, 0) ;

  if (keySetMax (newBadData))
    { displayCreate (DtKeySet) ; /* new window */
      keySetShow (newBadData, 0) ;
      newBadData = keySetCreate () ; /* make a new one for here */
    }
}

static void geneSummary (void)
{
  int i, j ;
  SEG *seg ;
  float x = 0, best ;
  KEY minLoc, maxLoc, clone ;
  OBJ obj ;
  FILE *fil ;
  VerticalMap look = currentVerticalMap("geneSummary") ;

  if (!(fil = filqueryopen (0, 0, "out", "a", "Output file to add genes to")))
    return ;

  getAllData (look) ;
  getBadData () ;

  for (j = 0 ; j < arrayMax(look->segs) ; ++j)
    { seg = arrp(look->segs,j,SEG) ;
      if (class(seg->key) != _VLocus ||
	  seg->flag & FLAG_HIGHLIGHT ||
	  seg->flag & FLAG_PHYS_GENE ||
	  !boundFind (seg->key, &dbnMin, &dbnMax, &minLoc, &maxLoc) ||
	  !makeDbn (seg->key))
	continue ;

      getPos (seg->key, &x) ;
      fprintf (fil,"%s	%s	%.2f	", name(seg->key), name(look->key), x) ;
      if (minLoc)
	fprintf (fil,"%s	", name(minLoc)) ;
      else
	fprintf (fil,"	") ;
      if (maxLoc)
	fprintf (fil,"%s	", name(maxLoc)) ;
      else
	fprintf (fil,"	") ;

      best = dbn[0] ;
      for (i = 1 ; i < NBIN ; ++i)
	if (dbn[i] > best)
	  best = dbn[i] ;
      best -= LOG_10 ;
      if (minLoc && dbn[0] > best)
	{ getPos (minLoc, &x) ;
          fprintf (fil,"%.2f	", x) ;
	}
      else
	{ for (i = 1 ; i < NBIN ; ++i)
	    if (dbn[i] > best)
	      break ;
	  fprintf (fil,"%.2f*	", DBN_X(i)) ;
	}
      if (maxLoc && dbn[NBIN-1] > best)
	{ getPos (maxLoc, &x) ;
          fprintf (fil,"%.2f	", x) ;
	}
      else
	{ for (i = NBIN-2 ; i >= 0 ; ++i)
	    if (dbn[i] > best)
	      break ;
	  fprintf (fil,"%.2f*	", DBN_X(i)) ;
	}
      if (seg->flag & FLAG_CLONED && (obj = bsCreate (seg->key)))
	{ if (bsGetKey (obj, _Positive_clone, &clone))
	    fprintf (fil, "%s", name(clone)) ;
	  bsDestroy (obj) ;
	}
      fprintf (fil,"\n") ;
    }

  filclose (fil) ;
}

static void rearrSummary (void)
{
  KEY rearr ;
  Array in, out ;
  KEY minLoc=0, maxLoc=0, maxBelowLoc, minAboveLoc ;
  int i, j ;
  float minIn, maxIn, maxBelow, minAbove, x, x1, x2 ;
  SEG *seg ;
  FILE *fil ;
  BOOL isBad ;
  VerticalMap look = currentVerticalMap("positionRearrangements") ;

  if (!(fil = filqueryopen (0, 0, "out", "a", "Output file to add rearr to")))
    return ;

  getAllData (look) ;
  getBadData () ;

  for (j = 0 ; j < arrayMax(look->segs) ; ++j)
    { seg = arrp(look->segs,j,SEG) ;
      if (class(seg->key) != _VInterval ||
	  seg->flag & FLAG_HIGHLIGHT)
	continue ;

      rearr = seg->key ;
      if (!assFind (look->data->loc2inside, assVoid(rearr), &in))
	continue ;

      minIn = 1000000 ;
      maxIn = -1000000 ;
      for (i = keySetMax(in) ; i-- ;)
	if (getPos (keySet(in,i), &x))
	  { if (x < minIn) { minIn = x ; minLoc = keySet(in,i) ; }
	    if (x > maxIn) { maxIn = x ; maxLoc = keySet(in,i) ; }
	  }
      if (minIn > 999999)
	continue ;

      maxBelow = -1000000 ;
      minAbove = 1000000 ;
      maxBelowLoc = 0 ;
      minAboveLoc = 0 ;
      isBad = FALSE ;
      if (assFind (look->data->loc2outside, assVoid(rearr), &out))
	for (i = keySetMax(out) ; i-- ;)
	  if (getPos (keySet(out,i), &x))
	    {
	      if (x < minIn)
		{ if (x > maxBelow) 
		    { maxBelow = x ; maxBelowLoc = keySet(out,i) ; }
		}
	      else if (x > maxIn) 
		{ if (x < minAbove) 
		    { minAbove = x ; minAboveLoc = keySet(out,i) ; }
		}
	      else
		isBad = TRUE ;
	    }

      if (maxBelow < -999999) 
	x1 = minIn ;
      else 
	x1 = 0.5 * (minIn + maxBelow) ;
      if (minAbove > 999999) 
	x2 = maxIn ;
      else 
	x2 = 0.5 * (maxIn + minAbove) ;

      fprintf (fil, "%s	%s	%.2f	%.2f	",
	       name(seg->key), name(look->key), x1, x2) ;
      if (maxBelowLoc) 
	fprintf (fil, "%s", name(maxBelowLoc)) ;
      fprintf (fil, "	%s	%s	", name(minLoc), name(maxLoc)) ;
      if (minAboveLoc) 
	fprintf (fil, "%s", name(minAboveLoc)) ;
      if (isBad)
	fprintf (fil, "	*") ;
/*
	{ fprintf (fil, "	|") ;
	  for (i = keySetMax(out) ; i-- ;)
	    if (getPos (keySet(out,i), &x) && x > minIn && x < maxIn)
	      fprintf (fil, " %8s", name(keySet(out,i))) ;
	}
*/
      fprintf (fil, "\n") ;
    }

  filclose (fil) ;
}

static void shortSummary (void)
{
  int i ;
  float x = 0 ;
  SEG *seg ;
  FILE *fil ;
  VerticalMap look = currentVerticalMap("shortSummary") ;

  if (!(fil = filqueryopen (0, 0, "out", "a", 
			    "Output file to add short to")))
    return ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (class(seg->key) != _VLocus ||
	  seg->flag & FLAG_HIGHLIGHT)
	continue ;

      getPos (seg->key, &x) ;
      if (seg->flag & FLAG_PHYS_GENE)
	fprintf (fil, "%s	%s	%.2f	P\n", 
		 name(seg->key), name(look->key), x) ;
      else
	fprintf (fil, "%s	%s	%.2f\n", 
		 name(seg->key), name(look->key), x) ;
    }

  filclose (fil) ;
}

MENUOPT vMapDataMenu[] = {
  { clearMapData,	"Clear Map Data" },
  { gjmOptAll,		"Global max likelihood" },
  { highlightDataless,	"Highlight unbound" },
  { check2point,	"Check 2 point" },
  { regetData,		"Update data" },
  { positionIntervals,	"Position Rearrangements" },
  { calcAll2pt,		"Calc all 2point" },
  { geneSummary,	"Gene Summary" },
  { rearrSummary,	"Rearr Summary" },
  { shortSummary,	"Short Summary" },
  { 0, 0 }
} ;

/********************************************************/
/************** Multi_pt_data ***************************/

static void addMultiData (OBJ obj)
{
  int i, count ;
  KEY locus ;

  for (i = 0 ; ; ++i)
    { if (!bsGetKey (obj, _bsRight, &locus))
	{ messerror ("No terminating locus in multi results") ;
	  return ;
	}
      if (i >= arrayMax(loci))
	array (loci, i, KEY) = locus ;
      else if (locus != array(loci,i,KEY))
	{ messerror ("Mismatching locus %s in multi results", 
		     name(locus)) ;
	  return ;
	}
      if (!bsGetData (obj, _bsRight, _Int, &count))
	return ;
      if (i >= arrayMax(counts))
	array(counts,i,int) = count ;
      else
	array(counts,i,int) += count ;
    }
}

static BOOL getMulti (KEY key)
{ 
  int index, i ; void *v ;
  OBJ obj ;
  static Associator keyCache = 0 ;
  static Array lociCache, countsCache ;

  dataKey = key ;

  if (!keyCache)
    { keyCache = assCreate () ;
      lociCache = arrayCreate (1024,KEY) ;
      countsCache = arrayCreate (1024,int) ;
    }

  arrayMax(loci) = 0 ;
  arrayMax(counts) = 0 ;

  if (assFind (keyCache, assVoid(key), &v))
    { index = assInt(v) ;
      for (i = 0 ; ; ++i, ++index)
	{ array (loci, i, KEY) = arr (lociCache, index, KEY) ;
	  if (arr(countsCache, index, int) < 0)
	    break ;
	  array (counts, i, int) = arr (countsCache, index, int) ;
	}
      return (arrayMax(loci) > 2) ;
    }

  if (!(obj = bsCreate (key)))
    return FALSE ;

  if (bsFindTag (obj, _Combined))
    addMultiData (obj) ;
  else 
    { if (bsFindTag (obj, _A_non_B))
	addMultiData (obj) ;
      if (bsFindTag (obj, _B_non_A))
	addMultiData (obj) ;
    }

  index = arrayMax (lociCache) ;
  assInsert (keyCache, assVoid(key), assVoid(index)) ;
  for (i = 0 ; ; ++i, ++index)
    { array(lociCache, index, KEY) = arr(loci, i, KEY) ;
      if (i >= arrayMax(counts))
	break ;
      array(countsCache, index, int) = arr(counts, i, int) ;
    }
  array(countsCache, index, int) = -1 ;	/* terminator */

  bsDestroy (obj) ;
  return (arrayMax(loci) > 2) ;
}

static SEG* insertMulti (VerticalMap look, KEY key) 
{ 
  SEG *seg ; 
  float y, min, max ;
  int i ;

  if (!look->data->key2seg)
    look->data->key2seg = assCreate () ;
  if (assFind (look->data->key2seg, assVoid(key), &seg))
    return seg ;

  if (!getMulti (key))
    {
#if !defined(MACINTOSH)
      if (iskey(key) > 1)	/* i.e. not just a name */
	fprintf (stderr, "Bad multi datum %s", name (key)) ;
#endif
      return 0 ;
    }
 
  min = 1000000 ; max = -1000000 ;
  for (i = 0 ; i < arrayMax(loci) ; ++i)
    if (getPos (arr(loci,i,KEY), &y))
      { if (y < min) min = y ;
	if (y > max) max = y ;
      }

  if (min > max)
    return 0 ;

  seg = arrayp(look->segs,arrayMax(look->segs),SEG) ;
  seg->key = key ;
  seg->flag = 0 ;
  seg->x = min  ; 
  seg->dx = max - min ;
  assInsert (look->data->key2seg, assVoid(key), seg) ;
  return seg ;
}

/*******************/

static BOOL logLikeMulti (float *ll)
{
	/* assumes single recombinant in locA, locB interval */
  int i = 0, n = arrayMax(loci) ;
  float y, z, denom, sum, least, most ;

  if (!getPos (arr(loci, 0, KEY), &y) || 
      !getPos (arr(loci, n-1, KEY), &z) ||
      y == z)
    return FALSE ;
  
  *ll = 0 ;
  sum = 0 ;			/* number of recombinants */

  if (z > y)
    { denom = z - y ;
      while (i < n)
	{ least = y ; most = y ;
	  while (++i < n && arr(counts,i-1,int) == 0)
	    if (getPos(arr(loci,i,KEY), &y))
	      {
		if (y < least)
		  least = y ;
		else if (y > most)
		  most = y ;
	      }
	  --i ;

	  if (sum)
	    {
	      if (least <= z)		/* order violates recombinants */
		{ keySetInsert (newBadData, dataKey) ;
		  return FALSE ;
		}
	      else
		{ *ll += sum * log((least - z) / denom) ;
		  sum = 0 ;
		}
	    }
	  z = most ;

	  while (++i < n)
	    { sum += arr(counts,i-1,int) ;
	      if (getPos(arr(loci,i,KEY), &y))
		break ;
	    }
	}
    }
  else
    { denom = y - z ;
      while (i < n)
	{ least = y ; most = y ;
	  while (++i < n && arr(counts,i-1,int) == 0)
	    if (getPos(arr(loci,i,KEY), &y))
	      {
		if (y < least)
		  least = y ;
		else if (y > most)
		  most = y ;
	      }
	  --i ;

	  if (sum)
	    {
	      if (most >= z)	/* order violates recombinants */
		{ keySetInsert (newBadData, dataKey) ;
		  return FALSE ;
		}
	      else
		{ *ll += sum * log((z - most) / denom) ;
		  sum = 0 ;
		}
	    }
	  z = least ;

	  while (++i < n)
	    { sum += arr(counts,i-1,int) ;
	      if (getPos(arr(loci,i,KEY), &y))
		break ;
	    }
	}
    }

  return TRUE ;
}

/*******************/

static void multiBoundCheck (KEY locus, float *min, float *max, 
		      KEY *minLoc, KEY *maxLoc)
{
	/* assumes single recombinant in locA, locB interval */

  int i, n = arrayMax(loci) ;
  float x, y ;

  getPos (locus, &x) ;

  for (i = 0 ; i < n ; ++i)
    if (arr(loci, i, KEY) == locus)
      break ;
  while (i < n)
    if (arr(counts, i++, int))
      break ;
  while (i < n)
    if (getPos (arr(loci, i++, KEY), &y))
      {
	if (y < x && y > *min)
	  { *min = y ; *minLoc = arr(loci, i-1, KEY) ; }
	else if (y > x && y < *max)
	  { *max = y ; *maxLoc = arr(loci, i-1, KEY) ; }
      }

  for (i = n-1 ; i >= 0 ; --i)
    if (arr(loci, i, KEY) == locus)
      break ;
  while (--i >= 0)
    if (arr(counts, i, int))
      break ;
  while (i >= 0)
    if (getPos (arr(loci, i--, KEY), &y))
      {
	if (y < x && y > *min)
	  { *min = y ; *minLoc = arr(loci, i+1, KEY) ; }
	else if (y > x && y < *max)
	  { *max = y ; *maxLoc = arr(loci, i+1, KEY) ; }
      }
}

/*******************/

static void drawMultiItem (VerticalMap look, SEG *seg, float x) 
{
  int i, n , sens = 0 ;
  float y1=0, y2=0, ymin, ymax ;
  
  if (!getMulti (seg->key))
    return ;
  n = arrayMax(loci) ;
 
  graphPointsize (0.8) ;
  if (getPos (arr(loci, 0, KEY), &y1))
    graphPoint (x, MAP2GRAPH(look->map, y1)) ;
  ymin = ymax = y1 ;
  if (getPos (arr(loci, n-1, KEY), &y2))
    graphPoint (x, MAP2GRAPH(look->map, y2)) ;
  
  graphPointsize (0.6) ;
  for (i = 1 ; i < n ; ++i)	/* redo n-1 to get counts line */
    if (getPos (arr(loci,i,KEY), &y2))
      { graphPoint (x, MAP2GRAPH(look->map, y2)) ;
	if (y2 != y1 && arr (counts, i-1, int))
	  { 
	    if (sens)
	      { if ((y2-y1) * sens  < 0)
		  { graphColor (RED) ;
		    graphFillRectangle (x - .2, MAP2GRAPH(look->map, y1), 
			       x + .2, MAP2GRAPH(look->map, y2)) ;
		  }
	      }
	    else
	      sens = y2 - y1  > 0 ? 1 : -1 ;
	    graphLine (x-0.25, MAP2GRAPH(look->map, (y1+y2)/2),
		       x+0.25, MAP2GRAPH(look->map, (y1+y2)/2)) ;
	    graphColor (BLACK) ;
	  }
	y1 = y2 ;
	if (y2 < ymin) ymin = y2 ;
	if (y2 > ymax) ymax = y2 ;
      }
  graphLine (x, MAP2GRAPH(look->map, ymin), 
	     x, MAP2GRAPH(look->map, ymax)) ;

}

void vMapDrawMultiPt (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  int   i, off = *offset ;
  SEG   *seg ;

  loci = arrayReCreate (loci, 8, KEY) ;
  counts = arrayReCreate (counts, 8, int) ;

  setPos (look) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (class(seg->key) == _VMulti_pt_data &&
	  !(seg->flag & FLAG_HIDE))
	{ array(look->boxIndex,graphBoxStart(),int) = i ;
	  drawMultiItem (look, seg, ++off) ;
	  graphBoxEnd() ;
	}
    }
  *offset = off+1 ;
}

void vMapReportMultiPt (VerticalMap look, SEG *seg)
{
  int i ;
  char *buf = look->messageText ;

  *buf = 0 ;
  if (!getMulti (seg->key))
    return ;

  strcpy (look->messageText, "Multi_pt ") ;
  strcat (look->messageText, name(dataKey)) ;
  strcat (look->messageText, ": ") ;
  strcat (look->messageText, name(arr(loci,0,KEY))) ;
  for (i = 1 ; i < arrayMax(loci) ; ++i)
    strcat (look->messageText, messprintf (" %d %s",
	    arr(counts,i-1,int), name(arr(loci,i,KEY)))) ;
}

/********************************************************/
/************************ 2 point ***********************/

static BOOL get2p (KEY key, KEY *loc1, KEY *loc2, KEY *type)
{ 
  OBJ obj ;
  BOOL result ;
  int x, i, index ;
  float mean, err ;
  static Associator keyCache = 0 ;
  static Array cache ;
  void* v, *vk = assVoid(key) ;
  dataKey = key ;

  if (!keyCache)
    { keyCache = assCreate () ;
      cache = arrayCreate (1024,KEY) ;
    }

  arrayMax(counts) = 0 ;
      
  if (assFind (keyCache, vk, &v))
    { index = assInt (v) ;
      *loc1 = arr(cache, index++, KEY) ;
      *loc2 = arr(cache, index++, KEY) ;
      *type = arr(cache, index++, KEY) ;
      for (i = 0 ; arr(cache, index, int) >= 0 ; ++i, ++index)
	array (counts, i, int) = arr (cache, index, int) ;
      return TRUE ;
    }

  if (!(obj = bsCreate (key)))
    return FALSE ;

  result = bsGetKey (obj, _Locus_1, loc1) &&
           bsGetKey (obj, _Locus_2, loc2) ;
  if (result)
    { 
      index = arrayMax(cache) ;
      assInsert (keyCache, assVoid(key), assVoid(index)) ;
      array(cache, index++, KEY) = *loc1 ;
      array(cache, index++, KEY) = *loc2 ;
      if (bsGetKeyTags (obj, _Calculation, type))
	{ array(cache, index++, KEY) = *type ;
	  for (i = 0 ; bsGetData (obj, _bsRight, _Int, &x) ; ++i)
	    array(counts,i,int) = x ;
	  for (i = 0 ; i < arrayMax(counts) ; ++i, ++index)
	    array(cache, index, int) = arr(counts, i, int) ;
	}
      else
	{ *type = 0 ;
	  mean = err = 0 ;
	  bsGetData (obj, _Distance, _Float, &mean) ;
	  bsGetData (obj, _Error, _Float, &err) ;
	  array(cache, index++, KEY) = 0 ;
	  array(counts, 0, KEY) = (KEY)(1000.0 * mean) ;
	  array(counts, 1, KEY) = (KEY)(1000.0 * err) ;
	  array(cache, index++, KEY) = (KEY)(1000.0 * mean) ;
	  array(cache, index++, KEY) = (KEY)(1000.0 * err) ;
	}
      array(cache, index, int) = -1 ; /* terminator */
    }

  bsDestroy (obj) ;
  return result ;
}

static SEG* insert2p (VerticalMap look, KEY key) 
{
  SEG *seg ; 
  float y1, y2 ;
  KEY loc1, loc2, type ;

  if (!look->data->key2seg)
    look->data->key2seg = assCreate () ;
  if (assFind (look->data->key2seg, assVoid(key), &seg))
    return seg ;

  if (!get2p (key, &loc1, &loc2, &type) ||
      !getPos(loc1, &y1) || !getPos(loc2, &y2))
    return 0 ;

  seg = arrayp(look->segs, arrayMax(look->segs), SEG) ;
  seg->key = key ;
  seg->flag = 0 ;
  if (y1 > y2)
    { seg->x = y2 ; seg->dx = y1 - y2 ; }
  else
    { seg->x = y1 ; seg->dx = y2 - y1 ; }
  assInsert (look->data->key2seg, assVoid(key), seg) ;
  return seg ;
}

/****************/

static BOOL check2pCounts (KEY type)	/* check arrayMax(counts) */
{
  if (type ==  _Full ||			/* WT X Y XY */
      type == _Backcross ||		/* WT X Y XY */
      type == _Sex_full)		/* WT X Y XY */
    {
      if (arrayMax(counts) != 4)
	{ messerror ("Bad counts size = %d for 2 point type %s",
		     arrayMax(counts), name (type)) ;
	  return FALSE ;
	}
    }
  else if (type == _Recs_all)		/* X Y ALL */
    {
      if (arrayMax(counts) != 3)
	{ messerror ("Bad counts size = %d for 2 point type %s",
		     arrayMax(counts), name (type)) ;
	  return FALSE ;
	}
    }
  else if
    ( type == 0 ||	/* no type given - counts contains distance, err */
    type == _One_recombinant ||	/* WT X */
    type == _Selected ||		/* X XY */
    type == _One_all ||		/* X ALL */
    type == _One_let ||		/* X ALL */
    type == _Tested ||		/* H X */
    type == _Selected_trans ||	/* X XY */
    type == _Back_one ||		/* WT X */
    type == _Sex_one ||		/* WT X */
    type == _Sex_cis ||		/* X ALL */
    type == _Dom_one ||		/* WT nonWT */
    type == _Dom_selected ||		/* WT X */
    type == _Dom_semi ||		/* XD ALL */
    type == _Dom_let ||		/* WT ALL */
    type == _Direct ||		/* R T */
    type == _Complex_mixed )	/* X ALL */
      {
	if (arrayMax(counts) != 2)
	  { if (type)
	      messerror ("Bad counts size = %d for 2 point type %s",
			 arrayMax(counts), name (type)) ;
	    return FALSE ;
	  }
      }
  
  return TRUE ;
}

/***********/

static float logLike2pt (float dist, KEY type)
{
  /* log probability of observed counts at given distance */

  float p = dist / 100.0 ; /* p(recombination) */
	/* if no interference 0.5 * (1 - exp(-0.02 * dist)) */
  float XY = (1-p) * (1-p) ;
  float W = 2 + XY ;		/* 3 - 2*p + p*p */
  float X = 1 - XY ;
#define N1 (arr(counts,0,int))
#define N2 (arr(counts,1,int))
#define N3 (arr(counts,2,int))
#define N4 (arr(counts,3,int))

  if (!check2pCounts (type))
    return 0 ;
  if (p < 0.0001)
    { p = 0.0001 ;
      XY = (1-p) * (1-p) ;
      W = 2 + XY ;
      X = 1 - XY ;
    }
  if (p > 0.9999)
    { p = .9999 ;
      XY = (1-p) * (1-p) ;
      W = 2 + XY ;
      X = 1 - XY ;
    }

  
/* segregation of recessives from xy/++ hermaphrodites */
    if (type == _Full)		/* WT X Y XY */
      return N1*log(W) + (N2+N3)*log(X) + N4*log(XY) ;
    else if (type == _One_recombinant)	/* WT X */
      return N1*log(W) + N2*log(X) ;
    else if (type == _Selected)		/* X XY */
      return N1*log(X) + N2*log(XY) ;
    else if (type == _One_all)		/* X ALL */
      return N1*log(X) + (N2-N1)*log(W+X+XY) ;
    else if (type == _Recs_all)		/* X Y ALL */
      return (N1+N2)*log(X) + (N3-N1-N2)*log(W+XY) ;
/* segregration from xy+/++z  - z is linked recessive lethal */
    else if (type == _One_let)		/* X ALL */
      return N1*log(X) + (N2-N1)*log(2-X) ;
/* recessive segregation from x+/+y hermaphrodites */
    else if (type == _Tested)		/* H X */
      return N1*log(2*p*(1-p)) + (N2-N1)*2*log(1-p) ;
    else if (type == _Selected_trans)	/* X XY */
      return N1*log(1-p*p) + N2*log(p*p) ;
    else if (type == _Backcross)		/* WT X Y XY */
      return (N1+N4)*log(1-p) + (N2+N3)*log(p) ;
    else if (type == _Back_one)		/* WT X */
      return N1*log(1-p) + N2*log(p) ;
    else if (type == _Sex_full)		/* WT X Y XY */
      return (N1+N4)*log(p) + (N2+N3)*log(1-p) ;
    else if (type == _Sex_one)		/* WT X */
      return N1*log(p) + N2*log(1-p) ;
    else if (type == _Sex_cis)		/* X ALL */
      return N1*log(p) + (N2-N1)*log(2-p) ;
    else if (type == _Dom_one)		/* WT nonWT */
      return N1*log(X) + N2*log(W+X+XY) ;
    else if (type == _Dom_selected)		/* WT X */
      return N1*log(X) + N2*log(XY) ;
    else if (type == _Dom_semi)		/* XD ALL */
      return N1*log(2*p*(1-p)) + (N2-N1)*log(4-2*p*(1-p)) ;
    else if (type == _Dom_let)		/* WT ALL */
      return N1*log(X) + (N2-N1)*log(W) ;
    else if (type == _Direct)		/* R T */
      return N1*log(p) + (N2-N1)*log(1-p) ; 
    else if (type == _Complex_mixed)	/* X ALL */
      return N1*log(X) + (N2-N1)*log(5-X) ;
    else if (type == 0)		/* Gaussian N1 = mu, N2 = err (2sd) */
      {
	if (N2)
	  { p = (N1 - 1000*dist) / (0.5*N2) ;
	    return -p*p/2 ;
	  }
	else
	  return 0 ;
      }
    else
      {
	messerror ("Unknown calculation type %s for 2point",
		   name (type)) ;
	return 0 ;
      }
}

/****************/

static BOOL best2p (KEY type, float *best, float *lo, float *hi)
{
  int n ;
  float p = 0, p22 = 0, p12 = 0 ;
  float testScore, x, x1, x2 ;
#define N1 (arr(counts,0,int))
#define N2 (arr(counts,1,int))
#define N3 (arr(counts,2,int))
#define N4 (arr(counts,3,int))

  if (!check2pCounts (type))
    return 0 ;

/* segregation of recessives from xy/++ hermaphrodites */
    if (type ==  _Full)			/* WT X Y XY */
      p22 = 2.0 * (N2+N3) / (float)(N1+N2+N3+N4) ;
    else if (type ==  _One_recombinant)	/* WT X */
      p22 = 3.0 * N2 / (float)(N1+N2) ;
    else if (type ==  _Selected)		/* X XY */
      p22 = N1 / (float)(N1+N2) ; 
    else if (type ==  _One_all)		/* X ALL */
      p22 = 4 * N1 / (float)N2 ;
    else if (type ==  _Recs_all)		/* X Y ALL */
      p22 = 2 * (N1+N2) / (float)N3 ; 			
/* segregration from xy+/++z  - z is linked recessive lethal */
    else if (type ==  _One_let)		/* X ALL */
      p22 = N1 / (float)N2 ; 				
/* recessive segregation from x+/+y hermaphrodites */
    else if (type ==  _Tested)		/* H X */
      p = N1 / (float)(2*N2 - N1) ; 			
    else if (type ==  _Selected_trans)	/* X XY */
      p = sqrt (N2 / (float)(N1+N2)) ; 			
    else if (type ==  _Backcross)		/* WT X Y XY */
      p = (N2+N3) / (float)(N1+N2+N3+N4) ; 		
    else if (type ==  _Back_one)		/* WT X */
      p = N2 / (float)(N1+N2) ; 			
    else if (type ==  _Sex_full)		/* WT X Y XY */
      p = (N1+N4) / (float)(N1+N2+N3+N4) ; 		
    else if (type ==  _Sex_one)		/* WT X */
      p = N1 / (float)(N1+N2) ; 			
    else if (type ==  _Sex_cis)		/* X ALL */
      p = 2 * N1 / (float)N2 ; 				
    else if (type ==  _Dom_one)		/* WT nonWT */
      p22 = 4 * N1 / (float)(N1 + N2) ;			
    else if (type ==  _Dom_selected)		/* WT X */
      p22 = N1 / (float)(N1 + N2) ;			
    else if (type ==  _Dom_semi)		/* XD ALL */
      p12 = 2 * N1 / (float)N2 ;			
    else if (type ==  _Dom_let)		/* WT ALL */
      p22 = 3 * N1 / (float)N2 ;			 
    else if (type ==  _Direct)		/* R T */
      p = N1 / (float)N2 ;				
    else if (type ==  _Complex_mixed)	/* X ALL */
      p22 = 5 * N1 / (float)N2 ;			
    else if (type ==  0)			/* Gaussian N1 = mu, N2 = err */
      {
	*best = N1 / 1000.0 ;
	*hi = (N1 + N2) / 1000.0 ;
	*lo = (N1 - N2) / 1000.0 ;
	if (*lo < 0) *lo = 0 ;
	return (N2 > 0) ;
      }
    else
      {
	messerror ("Unknown calculation type %s for 2point",
		   name (type)) ;
	return FALSE ;
      }

		/* solve quadratic eqn if necessary */
  if (p22)			/* p22 = 2p-p^2 */
    p = 1 - sqrt(1 - p22) ;
  else if (p12)			/* p12 = p-p^2 */
    p = (1 - sqrt(1 - 4*p12)) / 2 ;

  *best = 100 * p ;
				/* find *lo, *hi by splitting */
  testScore = logLike2pt (*best, type) - LOG_10 ;

  if (logLike2pt (0, type) > testScore)
    *lo = 0 ;
  else
    { x1 = 0 ; x2 = *best ;
      for (n = 6 ; n-- ; )
	{ x = (x1 + x2) / 2 ;
	  if (logLike2pt (x, type) < testScore)
	    x1 = x ;
	  else
	    x2 = x ;
	}
      *lo = (x1 + x2) / 2 ;
    }

  if (!*best)
    { x1 = 0.01 ;
      if (logLike2pt (x1, type) > testScore)      
	{ *hi = x1 ;
	  return TRUE ;
	}
    }
  else
    x1 = *best ;
  x2 = 2*x1 ;
  while (x2 < 100.0 && logLike2pt (x2, type) > testScore) x2 *= 2 ;
  for (n = 6 ; n-- ; )
    { x = (x1 + x2) / 2 ;
      if (x >= 100) x = 99.9 ;
      if (logLike2pt (x, type) < testScore)
	x2 = x ;
      else
	x1 = x ;
    }
  *hi = x ;

  return TRUE ;
}

/******* jean's 2pt draw -- richard's in gmapdata2.c ********/

static void  draw2ptItem(VerticalMap look, SEG *seg,float x) 
{
  KEY loc1, loc2 ;
  KEY type ;
  float y1, y2, c, best, lo, hi ;
 
  if (!get2p (seg->key, &loc1, &loc2, &type) ||
      !getPos (loc1, &y1) || !getPos (loc2, &y2))
    return ;

  if (y1 > y2)
    { float tmp = y1 ; 
      y1 = y2 ;
      y2 = tmp ;
    }
  
     /* join the 2 genes */
  graphLine (x, MAP2GRAPH(look->map, y1), 
	     x, MAP2GRAPH(look->map, y2)) ;


  if (best2p (type, &best, &lo, &hi))
    { c = (y1 + y2) / 2 ;

      /* a red rectangle from true to best positions */
      graphColor(RED) ;
      
      graphFillRectangle(x-.2,MAP2GRAPH(look->map,y1),x+.2,MAP2GRAPH(look->map,c - best/2)) ;
      graphFillRectangle(x-.2,MAP2GRAPH(look->map,y2),x+.2,MAP2GRAPH(look->map,c + best/2)) ;
      
      graphColor(BLACK) ;
      graphRectangle(x-.2,MAP2GRAPH(look->map,y1),x+.2,MAP2GRAPH(look->map,c - best/2)) ;
      graphRectangle(x-.2,MAP2GRAPH(look->map,y2),x+.2,MAP2GRAPH(look->map,c + best/2)) ;


      /* a green error box around best position, hence hiding the red box */
      graphColor(GREEN) ;
      graphFillRectangle(x-.25,MAP2GRAPH(look->map,c - lo/2),
			 x+.25,MAP2GRAPH(look->map,c - hi/2)) ;
      graphFillRectangle(x-.25,MAP2GRAPH(look->map,c + lo/2),
			 x+.25,MAP2GRAPH(look->map,c + hi/2)) ;
      
      graphColor(BLACK) ;
      graphRectangle(x-.25,MAP2GRAPH(look->map,c - lo/2),
		     x+.25,MAP2GRAPH(look->map,c - hi/2)) ;
      graphRectangle(x-.25,MAP2GRAPH(look->map,c + lo/2),
		     x+.25,MAP2GRAPH(look->map,c + hi/2)) ;
    }

    /* point the genes, must come last */
  graphPointsize(.7) ;
  graphPoint(x,MAP2GRAPH(look->map,y1)) ;
  graphPoint(x,MAP2GRAPH(look->map,y2)) ;
}

void vMapDraw2pt (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  int   i, off = *offset ;
  SEG   *seg ;

  counts = arrayReCreate (counts, 4, 0) ;

  setPos (look) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (class(seg->key) == _V2_point_data &&
	  !(seg->flag & FLAG_HIDE))
	{ array(look->boxIndex,graphBoxStart(),int) = i ;
	  draw2ptItem (look, seg, ++off) ;
	  graphBoxEnd() ;
	}
    }
  *offset = off+1 ;
}

void vMapReport2pt (VerticalMap look, SEG *seg)
{
  KEY loc1, loc2, type ;
  int i ;
  float best, lo, hi, y1, y2 ;
  char *buf = look->messageText ;

  *buf = 0 ;
  if (!get2p (seg->key, &loc1, &loc2, &type))
    return ;

  strncpy (buf, messprintf ("2point %s: %s %s: %s",
    name(seg->key), name(loc1), name(loc2), name(type)), 100) ;

  if (type)			/* 27 chars plenty for counts */
    for (i = 0 ; i < arrayMax(counts) ; ++i)
      strcat (buf, messprintf (" %d", arr(counts, i, int))) ;
  else
    strcat (buf, messprintf (" %.2f %.2f",
	     arr(counts,0,int)/1000.0, arr(counts,1,int)/1000.0)) ;

  if (getPos (loc1, &y1) && getPos (loc2, &y2))
    strcat (buf, messprintf (": %.2f", (y1>y2) ? (y1-y2) : (y2-y1))) ;

  if (best2p (type, &best, &lo, &hi))
    strncat (buf, messprintf (" [%.2f, %.2f, %.2f]", 
	      lo, best, hi), 127 - strlen(buf)) ;
}

/******************/

static void calcAll2pt (void)
{
  static KEYSET done = 0 ;
  KEYSET kset ;
  void *k ;
  int i, n = 0 ;
  KEY loc1, loc2, type ;
  float best, lo, hi ;
  FILE *fil ;
  VerticalMap look = currentVerticalMap("calcAll2pt") ;

  if (!(fil = filqueryopen (0, 0, "ace", "a", "Output ace file:")))
    return ;

  getAllData (look) ;
  done = keySetReCreate (done) ;

  for (k = 0 ; assNext (look->data->loc2two, &k, &kset) ;)
    for (i = 0 ; i < keySetMax(kset) ; ++i)
      if (!keySetInsert (done, keySet(kset, i)) &&
	  get2p (keySet(kset, i), &loc1, &loc2, &type) &&
	  best2p (type, &best, &lo, &hi))
	{ fprintf (fil, "2_point_data \"%s\"\n", name (dataKey)) ;
	  fprintf (fil, "Calc_distance %.2f\n", best) ;
	  fprintf (fil, "Calc_lower_conf %.2f\n", lo) ;
	  fprintf (fil, "Calc_upper_conf %.2f\n", hi) ;
	  fprintf (fil, "\n") ;
	  ++n ;
	}

  filclose (fil) ;
}

/******************************************************************************/
/******************************************************************************/

