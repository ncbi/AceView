/*  File: gmapsubmapcol.c
 *  Author: Ian Longden (il@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: Submap display column for the genetic map
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  8 16:26 1999 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: gmapsubmapcol.c,v 1.3 2014/11/30 03:20:46 mieg Exp $ */

#include "gmap.h"
#include "query.h"
#include "bump.h"

/************************************************************/

#define submapmarkerFormat "kffff"

typedef struct {
  KEY key;
  float gpos;        /* global positon i.e. chrom map coordinates*/
  float lpos;        /* local position with regards to submap*/
                     /* to calculate pos use extent + lpos to
			get overall position in contig */
  float newpos;      /* global coors after moved by lcoi midpoint*/
  float posoncontig;   /* global coors with repsect to contig */
} SUBMAPMARKERS;


typedef struct {
  KEY key;
  int reject;                    /* True if markers lie outside the
				    tolerance contig (coors from global) */
  int rejectbyloci;              /* True if markers lie outside the
				    tolerance contig (coors calc from loci) */
  float gstart, gend;            /* global start/end on chrom map*/
  float extentLeft, extentRight; /* extent i.e. size */
  float lociLeft,lociRight;      /* coors calculated from markers global positions. */
  Array submapmarkers;           /* contains the above marker information */
  KEYSET negMarkers;             /* keyset of the markers which are rejected */
  KEYSET negMarkersbyLoci;       /* as above but for when contigs placed by loci */ 
  KEYSET markers;                /* a Keyset of all the markers */
} GMAPSUBMAPSEG;

#define submapFormat "kiiffffff"
typedef struct
{ KEY key;
  int reject;                    /* True if markers lie outside the tolerance contig (coors from global) */
  int rejectbyloci;              /* True if markers lie outside the tolerance contig (coors calc from loci) */
  float gstart, gend;            /* global start/end on chrom map*/
  float extentLeft, extentRight; /* extent i.e. size */
  float lociLeft,lociRight;      /* coors calculated from markers global positions. */
} SUBMAPTEMP;

#define LOCI 0
#define INTERVAL 1
#define LINK 2

typedef struct
{ KEY key;
  int type;
  int index,subindex;
} BOX2KEY;

#define ALL 1
#define REJECTED 2
#define NON_REJECTED 3
 
typedef struct SubPrivStruct {
  char *lociQuery, *intQuery;
  KEYSET lociKeyset, intervalKeyset;
  int show;     /* 1 ALL 2 rejected 3 Non rejected */
  int tolerance;
  BOOL posByLoci;
  BOOL showAllMarkers;
  float scale;
  Array segs;
  Array box2key;
  int currmap;
} *SUBPRIV;

/************************************************************/

int gMapSubMapOrder (const void *a,const  void *b);

/* we know ym as this is the top/bot of the screen Therefore only need to calculate xm */
BOOL getCutoffwithBorder(float x1, float y1, float x2, float y2, float ym, float *xm)
{
  float m;
  float c;
  
  m = (y2-y1)/(x2-x1);
  c = y1-(m*x1);
  if(m)
    *xm = (ym - c)/m;
  else
    return FALSE;
  return TRUE;
} /* getCutoffwithBorder */


BOOL gMapIsNegMarkers(MAPCONTROL map, KEY key)
{
  GeneticMap look = (GeneticMap)(map->look);
  
  /* make sure our map is in charge */
  if (map->control->activeInstance && 
      map == map->control->activeInstance->map)
    { 
      if (look->negativeMarkersValid)
	if(look->negativeMarkers)
	  if( keySetFind(look->negativeMarkers, key, 0))
	    return TRUE;
    }
  
  return FALSE;
} /* gMapIsNegMarkers */


static BOOL submapSetSelect(COLINSTANCE instance, int box, double x, double y)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  SUBPRIV private = instance->private;
  BOX2KEY *seg = arrp(private->box2key, box, BOX2KEY);
  GMAPSUBMAPSEG *seg2;
  GeneticMap look = (GeneticMap)(map->look);

  gMapUnselect(map);
  look->selectKey = seg->key;
  look->friendsInfoValid = TRUE;
  look->neighboursInfoValid = TRUE;
  (void)gMapNeighbours(map, &look->neighbours, seg->key);
  (void)gMapPositive(map, &look->positive, seg->key);
  if(seg->type == INTERVAL)
    { look->negativeMarkersValid = TRUE;
      (void)gMapNegative(map, &look->negative, seg->key);
      seg2 = arrp(private->segs, seg->index, GMAPSUBMAPSEG);
      if(private->posByLoci)
       look->negativeMarkers = seg2->negMarkersbyLoci;
      else
       look->negativeMarkers = seg2->negMarkers;
      private->currmap = seg->key;
    }
  else
    look->negativeMarkers = 0;

  control->activeKey = seg->key;
   if(seg->type == INTERVAL)
     controlDraw();
  return FALSE;
} /* submapSetSelect */


static BOOL submapUnSelect(COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);

  control->activeKey = 0;
  *look->messageText = 0;
  
  if(look->messageBox)
    graphBoxDraw(look->messageBox, -1, -1);
  
  return FALSE;
} /* submapUnSelect */


static BOOL getMarkerLocations(KEY markerkey, KEY subkey, KEY mainkey, float *gpos, float *lpos)
{
  OBJ obj;
  int i;
  Array loci;
  float loc=0.000005,glob=0.000005;

  obj = bsCreate(markerkey);
  if(!obj)
    return FALSE;

  
  *lpos = *gpos = -99.0;
  loci = arrayCreate(200,BSunit);
  if(bsFindTag(obj,_Map) && bsFlatten(obj, 3, loci)){
    if(arr(loci, 1, BSunit).k == str2tag("Position")){
      for(i=0; i+2<=arrayMax(loci); i+=3){
	if(arr(loci, i, BSunit).k == mainkey)
	  glob = arr(loci, i+2, BSunit).f;
	else if(arr(loci, i, BSunit).k == subkey)
	  loc = arr(loci, i+2, BSunit).f;
      }
    }
  }
  if(loc > 0.005 || loc < -0.005)  /* becouse we can get small error number from incorrect*/
    *lpos = loc;                     /* reading of number as what is being read may not be a number */
  if(glob > 0.005 || glob < -0.005)
    *gpos = glob;
  bsDestroy(obj);
  arrayDestroy(loci);
  if((*lpos != -99.0 && *gpos != -99.0) && (*lpos != 0.0 && *gpos != 0.0))
    return TRUE;
  else {                    /* perhaps it's a pac etc so find midpoint */
    obj = bsCreate(markerkey);
    if(!obj)
      return FALSE;
    loci = arrayCreate(200,BSunit);
    if(bsFindTag(obj,_Map) && bsFlatten(obj, 4, loci)){
	for(i=0; i+7<=arrayMax(loci); i+=8){
	  if(arr(loci, i, BSunit).k == mainkey)
	    glob = (arr(loci, i+3, BSunit).f+arr(loci, i+7, BSunit).f)/2.0;
	  else if(arr(loci, i, BSunit).k == subkey)
	    loc = (arr(loci, i+3, BSunit).f+arr(loci, i+7, BSunit).f)/2.0;
	}
      }

    bsDestroy(obj);
    arrayDestroy(loci);

    if(loc > 0.005 || loc < -0.005)  /* becouse we can get small error number from incorrect*/
      *lpos = loc;                     /* reading of number as what is being read may not be a number */
    if(glob > 0.005 || glob < -0.005)
      *gpos = glob;
    if((*lpos != -99.0 && *gpos != -99.0) && (*lpos != 0.0 && *gpos != 0.0))
      return TRUE;
  }      
  return FALSE;
} /* getMarkerLocations */


static void gMapGetSubMapData(KEY submapkey, MAPCONTROL map, Array segs)
{
  OBJ obj;
  SUBMAPMARKERS *marker;
  Array loci;
  GMAPSUBMAPSEG *seg;
  int i=0;
  KEY locus;
  float gpos, lpos, extentLeft, extentRight;
  static BOOL errormess=TRUE;

  obj = bsCreate(submapkey);
  if(!obj)
    return ;

  extentLeft = extentRight =  -10.0;
  
  if (!bsGetData(obj, str2tag("Extent"), _Float, &extentLeft) ||
      !bsGetData(obj, _bsRight, _Float, &extentRight)){
    if(errormess){
      messerror("Map %s is a submap of map %s, but does not have an extent",
		name(submapkey), name(map->key));
      errormess = FALSE;
    }
    return;
  }
  else
    {
      seg = arrayp(segs, arrayMax(segs), GMAPSUBMAPSEG);
      seg->submapmarkers = 0;
      seg->key = submapkey;
    }
  if(seg->extentLeft <= seg->extentRight)
    {seg->extentLeft = extentLeft;
     seg->extentRight = extentRight;
    }
  else
    { seg->extentLeft = extentRight;
      seg->extentRight = extentLeft;
   }
  loci = arrayCreate(200,BSunit);
  if(bsFindTag(obj,_Contains) && bsFlatten(obj, 2, loci)) /* all the contains bits */
    { for(i=1; i<arrayMax(loci); i+=2)
	{ locus = arr(loci, i, BSunit).k;
	  if(getMarkerLocations(locus, submapkey, map->key, &gpos, &lpos))
	    { if(!seg->submapmarkers)
		seg->submapmarkers = arrayCreate(10, SUBMAPMARKERS);
	      marker = arrayp(seg->submapmarkers, arrayMax(seg->submapmarkers), SUBMAPMARKERS);
	      marker->key = locus;
	      marker->lpos = lpos;
	      marker->gpos = gpos;
	    }
	}
    }



  /* need to test Map for the correct map. */


  if(bsFindTag(obj,str2tag("Map")) && bsFlatten(obj, 5, loci))
    { i = 0; 
      while(i+8 < arrayMax(loci)){
	if(arr(loci,i,BSunit).k == map->key){
	  seg->gstart = arr(loci, i+3, BSunit).f;
	  seg->gend = arr(loci, i+8, BSunit).f;
	  break;
	}
	i+=10;
      }
    }
  if(seg->gstart > seg->gend)
    { seg->gstart = arr(loci, i+8, BSunit).f;
      seg->gend = arr(loci, i+3, BSunit).f;
    }

  arrayDestroy(loci);

  bsDestroy(obj);
  return ;
} /* gMapGetSubMapData */


static void calcRejects(COLINSTANCE instance)
{
  int i,j;
  GMAPSUBMAPSEG *seg;
  SUBMAPMARKERS *marker;
  SUBPRIV private = instance->private;
  float percentTol = ((float)private->tolerance*0.01);
  float range;
  Array segs = private->segs;
  MAPCONTROL map = instance->map;
  GeneticMap look= (GeneticMap)(map->look);

  look->negativeMarkers = 0;
  for(i = 0; i < arrayMax(segs); i++)
    { seg = arrayp(segs, i, GMAPSUBMAPSEG);
      seg->rejectbyloci = seg->reject = FALSE;

      /* create the keysets for wether the markers lie within the tolerance */
      if(seg->negMarkers)
	keySetDestroy(seg->negMarkers);
      seg->negMarkers = keySetCreate();

      if(seg->negMarkersbyLoci)
	keySetDestroy(seg->negMarkersbyLoci);
      seg->negMarkersbyLoci = keySetCreate();	  

      if(seg->markers)
	keySetDestroy(seg->markers);
      seg->markers = keySetCreate();	  

      if(seg->key == private->currmap)
	{
	  if(private->posByLoci)
	    look->negativeMarkers = seg->negMarkersbyLoci;
	  else
	    look->negativeMarkers = seg->negMarkers;
	}
      

      if(seg->submapmarkers)
	{ range = seg->gend - seg->gstart;
	  for(j = 0; j < arrayMax(seg->submapmarkers); j++)
	    { 
	      /* check that the marker lies within the tolerance of the contig for the */
	      /*loci calculated positions */
	      marker = arrp(seg->submapmarkers,j,SUBMAPMARKERS);
	      keySet(seg->markers, arrayMax(seg->markers)) = marker->key;
	      /* ONLY REJECT THOSE BASED ON THOSE IN QUERY */
	      if(private->lociKeyset && !keySetFind(private->lociKeyset, marker->key, 0))
		continue;
	      
	      if(marker->gpos < seg->lociLeft)
		{ if((marker->gpos + (percentTol * range)) < seg->lociLeft)
		  { seg->rejectbyloci = TRUE; 
		    keySet(seg->negMarkersbyLoci, arrayMax(seg->negMarkersbyLoci)) = marker->key;
		  }
		}
	      else if( marker->gpos > seg->lociRight)
		{ if((marker->gpos - (percentTol * range)) > seg->lociRight)
		    { seg->rejectbyloci = TRUE;
		      keySet(seg->negMarkersbyLoci, arrayMax(seg->negMarkersbyLoci)) = marker->key;
		    } 
		}
	      /* calculate the same for the standard positions */
	      if(marker->gpos < seg->gstart)
		{ if((marker->gpos + (percentTol * range)) < seg->gstart)
		    { seg->reject = TRUE; 
		      keySet(seg->negMarkers, arrayMax(seg->negMarkers)) = marker->key;
		    }
		}
	      else if( marker->gpos > seg->gend)
		{ if((marker->gpos - (percentTol * range)) > seg->gend)
		    { seg->reject = TRUE;
		      keySet(seg->negMarkers, arrayMax(seg->negMarkers)) = marker->key;
		    }
		}
	    }
	  keySetSort(seg->negMarkersbyLoci);
	  keySetSort(seg->negMarkers);
	  keySetSort(seg->markers);
	}
    }
} /* calcRejects */


static void calcNewPos(COLINSTANCE instance)
/* The position of the interval is calculated from mean position of its markers.
The markers postion (to the interval) is then recalulated.
*/
{ SUBPRIV private = instance->private;
  double offset, scale=0.0;
  int i, j, count;
  GMAPSUBMAPSEG *seg;
  SUBMAPMARKERS *marker;
  float range,temp;
  
  for(i = 0; i < arrayMax(private->segs); i++)
    { seg = arrp(private->segs, i, GMAPSUBMAPSEG);
      if(private->scale!= 0.0)
	range = seg->extentRight - seg->extentLeft;
      else
	range = (seg->gend - seg->gstart);
      if(range){
	if(private->scale!= 0.0){
	  scale =  private->scale;
	  range*=private->scale;
	}
	else
	  scale = (seg->gend - seg->gstart) / (seg->extentRight - seg->extentLeft);
      }
      else 
	continue;
      if(seg->submapmarkers)
	/* POSITION ON ONLY THOSE IN QUERY */
	{ for(j = 0 , offset = 0, count = 0; j < arrayMax(seg->submapmarkers); j++)
	    { marker = arrayp(seg->submapmarkers, j, SUBMAPMARKERS);
	      if(private->lociKeyset)
		{ if(keySetFind(private->lociKeyset, marker->key, 0))
		   { offset += (double)marker->gpos;
		     count++;
		    }
		}  
	      else  
		{ offset += (double)marker->gpos;
		  count++;
		}
	    }  
	    
	  if(count)
	    { offset = offset / (float)count;
	      seg->lociLeft = offset - (range/2);
	      seg->lociRight = offset + (range/2);
	    }
	  else                                        /* ??? MAYBE THESE SHOULD BE IGNORED ???? */
	    { seg->lociLeft = seg->gstart;
	      seg->lociRight = seg->gend;
	    }
	  for(j = 0; j < arrayMax(seg->submapmarkers); j++)
	    { marker = arrp(seg->submapmarkers, j, SUBMAPMARKERS);
	      temp = (marker->lpos - seg->extentLeft) * scale;
	      /*  if(seg->key == private->currmap)
		printf("temp = %f,scale = %f,lpos = %f,locileft =%f,lociright =%f\n",temp,scale,marker->lpos,seg->lociLeft,seg->lociRight);*/
	      marker->newpos = seg->lociLeft + temp;
	      marker->posoncontig =  seg->gstart + temp;
	    }
	}
      else
	{ seg->lociLeft = seg->gstart;
	  seg->lociRight = seg->gend;
	}
    }
  arraySort(private->segs, gMapSubMapOrder);
  calcRejects(instance);

  return;
} /* calcNewPos */


int gMapSubMapOrder (const void *a, const void *b)  /* for arraySort() call */
{
  const GMAPSUBMAPSEG *seg1 = (const GMAPSUBMAPSEG*)a, *seg2 = (const GMAPSUBMAPSEG*)b ;
  float diff;
  
  /*  diff = seg1->gstart - seg2->gstart ;*/
  diff = seg1->lociLeft - seg2->lociLeft;

  if (diff > 0)
    return 1 ;
  else if (diff < 0)
    return -1 ;
  else
    return seg2->key - seg1->key ; /* never give a random result please */
} /* gMapSubMapOrder */


void *submapConvert2(MAPCONTROL map, void *params)
{
  Array segs = 0;
  return segs;
} /* submapConvert2 */


void *submapConvert(MAPCONTROL map)
{
  int i;
  Array loci;
  KEY locus;
  Array segs;
  char str1[200];
  KEY cache;
  GMAPSUBMAPSEG *seg=0;
  Array tempsegs;
  SUBMAPTEMP *tempseg;

  getPos(map, 0, 0);

  /*Are they Cached */
  if ((cache = gMapCacheKey(map->key, "subs")) &&
      (tempsegs = arrayGet(cache, SUBMAPTEMP, submapFormat)))
    { segs = arrayCreate(arrayMax(tempsegs), GMAPSUBMAPSEG);
      for(i = 0; i < arrayMax(segs); i++)
	{ tempseg = arrp(tempsegs, i, SUBMAPTEMP);
	  seg= arrp(segs, i, GMAPSUBMAPSEG);
	  seg->reject = tempseg->reject ;
	  seg->rejectbyloci = tempseg->rejectbyloci ;
	  seg->gstart = tempseg->gstart ;
	  seg->gend = tempseg->gend ;
	  seg->extentLeft = tempseg->extentLeft ;
	  seg->extentRight = tempseg->extentRight ;
	  seg->lociLeft = tempseg->lociLeft ;
	  seg->key = tempseg->key ;
	  seg->lociRight = tempseg->lociRight ;
	}
      arrayDestroy(tempsegs);
      for(i = 0; i < arrayMax(segs); i++)
	{ 
	  sprintf(str1,"subs.marker.%d", i);
	  if ((cache = gMapCacheKey(map->key, str1)))
	    { seg = arrp(segs, i, GMAPSUBMAPSEG);
	      seg->submapmarkers = arrayGet(cache, SUBMAPMARKERS, submapmarkerFormat);
	    }
	  else
	    seg->submapmarkers = 0;
	}
    return segs;
    }

  segs = arrayCreate (128, GMAPSUBMAPSEG) ; /* is destroyed explicitly by subprivdestroy */
  loci = gMapMapObjects(map->key); /* does not seem to get any map objects ???? */

  for (i = 0; i <arrayMax(loci); i++)
    { locus = arr(loci, i, KEYPAIR).obj;
      if(class(locus) == _VMap)
	gMapGetSubMapData(locus, map, segs);      /* the data for the sub map */
    }
  arrayDestroy(loci);
  /*save cache */
  arraySort(segs, gMapSubMapOrder);
  tempsegs = arrayCreate(arrayMax(segs), SUBMAPTEMP);
  for(i = 0; i < arrayMax(segs); i++)
    { tempseg = arrp(tempsegs, i, SUBMAPTEMP);
      seg= arrp(segs, i, GMAPSUBMAPSEG);
      tempseg->reject = seg->reject ;
      tempseg->rejectbyloci = seg->rejectbyloci ;
      tempseg->gstart = seg->gstart ;
      tempseg->gend = seg->gend ;
      tempseg->extentLeft = seg->extentLeft ;
      tempseg->extentRight = seg->extentRight ;
      tempseg->lociLeft = seg->lociLeft ;
      tempseg->key = seg->key ;
      tempseg->lociRight = seg->lociRight ;
    }
  
  gMapCacheStore(map->key, "subs", tempsegs, submapFormat);
  arrayDestroy(tempsegs);

  for(i=0;i<arrayMax(segs);i++)
    { sprintf(str1,"subs.marker.%d", i);
      seg = arrp(segs, i, GMAPSUBMAPSEG);
      if(seg->submapmarkers)
	gMapCacheStore(map->key, str1, seg->submapmarkers, submapmarkerFormat);
    }

  return (void *) segs; 
} /* submapConvert */


static void submapDraw(COLINSTANCE instance, float *offset)
{
  float y1, y2, y3, y4, x, z1, x3 = 0;
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  SUBPRIV private = instance->private;
  GeneticMap look= (GeneticMap)(map->look);
  GMAPSEG *seg2;
  GMAPSUBMAPSEG *seg;
  SUBMAPMARKERS *marker;
  BOX2KEY *box2key;
  int i, box, j=0, ix, cmax = 0,x2 = 0,k;
  BUMP bump;
  float lx1, ly1, lx2, ly2;
  *offset +=1;

  if(private->box2key)
    arrayDestroy(private->box2key);
  private->box2key = arrayCreate(10, BOX2KEY);

  /* a quick dummy run to see the  width (x3) */ 
  bump = bumpCreate(24, 0);
  for(i= 0; i< arrayMax(private->segs); i++)
    { seg = arrp(private->segs,i,GMAPSUBMAPSEG);
      if(private->intervalKeyset && !keySetFind(private->intervalKeyset, seg->key, 0))
	continue;
      if(keySetFind(look->hidden, seg->key, 0))
	continue;

      if(private->posByLoci)
	{ if(private->show == REJECTED && !seg->rejectbyloci)
	    continue;
	  else if(private->show == NON_REJECTED && seg->rejectbyloci)
	    continue;
	  y1 = MAP2GRAPH(map, seg->lociLeft);
	  y2 = MAP2GRAPH(map, seg->lociRight);
	}
      else
	{ if(private->show == REJECTED && !seg->reject)
	    continue;
	  else if(private->show == NON_REJECTED && seg->reject)
	    continue;
          y1 = MAP2GRAPH(map, seg->gstart);
	  y2 = MAP2GRAPH(map, seg->gend);
	}

      if(y2 < control->topMargin || y1 > control->graphHeight-1)
	continue;
      if(y1 < control->topMargin) y1 = control->topMargin;
      if(y2 > control->graphHeight) y2 = control->graphHeight; 
      z1 =y1;
      ix = 0;
      bumpItem(bump, 3, (y2-y1+1)+0.2, &ix, &z1);
      if(cmax < ix)
        cmax = ix;
    }
  x3 =*offset + (cmax) + 6.0;

  bumpDestroy(bump);
  cmax = 0;
  bump = bumpCreate(24, 0);
  for(i = 0; i < arrayMax(private->segs); i++)
    { seg = arrp(private->segs, i, GMAPSUBMAPSEG);

      if(private->intervalKeyset && !keySetFind(private->intervalKeyset, seg->key, 0))
	continue;
      if(keySetFind(look->hidden, seg->key, 0))
	continue;

      if(private->posByLoci)
	{ if(private->show == REJECTED && !seg->rejectbyloci)
	    continue;
	  else if(private->show == NON_REJECTED && seg->rejectbyloci)
	    continue;
	  y1 = MAP2GRAPH(map, seg->lociLeft);
	  y2 = MAP2GRAPH(map, seg->lociRight);
	}
      else
	{ if(private->show == REJECTED && !seg->reject)
	    continue;
	  else if(private->show == NON_REJECTED && seg->reject)
	    continue;
          y1 = MAP2GRAPH(map, seg->gstart);
	  y2 = MAP2GRAPH(map, seg->gend);
	}

      if(y2 < control->topMargin || y1 > control->graphHeight-1)
	continue;
      if(y1 < control->topMargin) y1 = control->topMargin;
      if(y2 > control->graphHeight) y2 = control->graphHeight; 
      z1 =y1;
      ix = 0;
      bumpItem(bump,3,(y2-y1+1)+0.2, &ix, &z1);
      x = *offset + ix;
      if(cmax < ix)
        cmax = ix;
      box = graphBoxStart();
      box2key = arrayp(private->box2key, box, BOX2KEY);
      box2key->key = seg->key;
      box2key->type = INTERVAL;
      box2key->index = i;
      controlRegBox(instance, box, box2key);
      if(y1 > control->topMargin)
	graphLine(x+0.25, y1, x+1.75, y1);
      if(y2 < control->graphHeight)
	graphLine(x+0.25, y2, x+1.75, y2);
      graphLine(x+0.25, y1, x+0.25, y2);
      graphLine(x+1.75, y1, x+1.75, y2);    
      graphBoxEnd();

      if(seg->submapmarkers)
	{ for(j = 0; j < arrayMax(seg->submapmarkers); j++)
	    { marker = arrp(seg->submapmarkers, j, SUBMAPMARKERS);
	    /* ONLY DRAW BAND IF MARKER IS IN THE KEYSET */
	      if(private->lociKeyset && !keySetFind(private->lociKeyset, marker->key, 0))
		continue;
	      if(private->posByLoci)
		y3 = MAP2GRAPH(map, marker->newpos);
	      else
		y3 = MAP2GRAPH(map, marker->posoncontig);
	      y4 = MAP2GRAPH(map, marker->gpos);
	      if((y3 > control->topMargin) && (y3 < control->graphHeight)){
		box = graphBoxStart();
		box2key = arrayp(private->box2key, box, BOX2KEY);
		box2key->key = marker->key;
		box2key->type = LINK;
		box2key->index = i;
		box2key->subindex = j;
		controlRegBox(instance, box, box2key);
		graphRectangle(x+1.75, y3-0.25, x+2.0, y3+0.25);
		graphBoxEnd();
	      }
	     
	      if(private->showAllMarkers || seg->key == private->currmap){
		if(gMapIsNeighbour(map, marker->key)){
		  lx1 = x+2.0;
		  lx2 = x3;
		  ly1 = y3;
		  ly2 = y4;
		  if((y4 < control->topMargin)){
		    getCutoffwithBorder(x+2.0, y3, x3, y4,control->topMargin,&lx2);
		    ly2 = control->topMargin;
		  }
		  if(y4 > control->graphHeight){
		    getCutoffwithBorder(x+2.0, y3, x3, y4,control->graphHeight,&lx2);
		    ly2 = control->graphHeight;
		  }
		  if(y3 < control->topMargin){
		    getCutoffwithBorder(x+2.0, y3, x3, y4,control->topMargin,&lx1);
		    ly1 = control->topMargin;
		  }
		  if(y3 > control->graphHeight){
		    getCutoffwithBorder(x+2.0, y3, x3, y4,control->graphHeight,&lx1);
		    ly1 = control->graphHeight;
		  }
		  graphLine(lx1, ly1, lx2, ly2);
		  /* printf("%s (%d) gpos = %f, lpos = %f, newpos =%f, posoncontig = %f\n",name(marker->key),marker->key,marker->gpos,marker->lpos,marker->newpos,marker->posoncontig);*/
		  
		}
	      }
	    		
	      if((y4 > control->topMargin) && (y4 < control->graphHeight)){
		if(gMapIsNeighbour(map, marker->key)){
		  if(private->showAllMarkers || seg->key == private->currmap){
		    /* graphLine(x+2.0, y3, x3, y4);*/
		    box = graphBoxStart();
		    box2key = arrayp(private->box2key, box, BOX2KEY);
		    box2key->key = marker->key;
		    box2key->type = LINK;
		    box2key->index = i;
		    box2key->subindex = j;
		    controlRegBox(instance, box, box2key);
		    graphRectangle(x3-0.25, y4-0.25, x3+0.25, y4+0.25);
		    graphBoxEnd();
		  }
		}
	      }
	    }
	}
    }
  bumpDestroy(bump);
  *offset = x3 +2.0;

  /* NEED TO DISPLAY ALL MARKERS  IN THE QUERY NOT JUST THOSE WHICH ARE GLOBAL?????? */
  bump = bumpCreate(1,0);
  if(private->showAllMarkers){
    for(i = 0; i< arrayMax(look->segs); i++)
      { seg2 = arrp(look->segs, i, GMAPSEG);
      /* ONLY DRAW BAND IF MARKER IS IN THE KEYSET */
      if(!(seg2->flag & FLAG_ANY_LOCUS))
	continue;
      if(private->lociKeyset && !keySetFind(private->lociKeyset, seg2->key, 0))
	continue;
      y4 = MAP2GRAPH(map,seg2->x);
      x2 = 0;
      if((y4 > control->topMargin +0.5 ) && (y4 < control->graphHeight))
	{ if(bumpText (bump, name(seg2->key), &x2, &y4, 1, TRUE) && 
	     y4 < control->graphHeight-0.5)
	  { box = graphBoxStart();
	    box2key = arrayp(private->box2key, box, BOX2KEY);
	    box2key->key = seg2->key;
	    box2key->type = LOCI;
	    box2key->index = i;
	    box2key->subindex = j;
	    controlRegBox(instance, box, box2key);
	    graphText(name(seg2->key), *offset, y4-0.5);
	    graphBoxEnd();
	  }
	}
    }
  }
  else
    { for(k = 0; k < arrayMax(private->segs); k++)
      { seg =  arrp(private->segs, k, GMAPSUBMAPSEG);
        if(seg->key == private->currmap){
	  for(i = 0; i < arrayMax(look->segs); i++)
	    { seg2 = arrp(look->segs, i, GMAPSEG);
	      if(seg->markers && !keySetFind(seg->markers, seg2->key, 0))
		continue;
	      y4 = MAP2GRAPH(map,seg2->x);
	      
	      x2 = 0;
	      if((y4 > control->topMargin +0.5 ) && (y4 < control->graphHeight))
		{ if(bumpText (bump, name(seg2->key), &x2, &y4, 1, TRUE) && 
		     y4 < control->graphHeight-0.5)
		  { box = graphBoxStart();
		  box2key = arrayp(private->box2key, box, BOX2KEY);
		  box2key->key = seg2->key;
		  box2key->type = LOCI;
		  box2key->index = i;
		  box2key->subindex = j;
		  controlRegBox(instance, box, box2key);
		  graphText(name(seg2->key), *offset, y4-0.5);
		  graphBoxEnd();
		  }
		}
	    }
	}
      }
    }
  bumpDestroy(bump);
  *offset += 12.0;

  return;
} /* submapDraw */


static void submapDoColour(COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  SUBPRIV private = instance->private;
  BOX2KEY *seg = arrp(private->box2key, box, BOX2KEY);
  GeneticMap look= (GeneticMap)(map->look);

  
  if((seg->key == look->selectKey) && (instance == control->activeInstance))
    { if(seg->type != LINK)
        graphBoxDraw(box, BLACK, CYAN);
      else
	graphBoxDraw(box, BLACK, CYAN);
      strncpy(look->messageText, name(seg->key), 100);
      if(look->messageBox)
	graphBoxDraw(look->messageBox, -1, -1);
      control->activeBox = box;
    }
  else if (gMapIsNegMarkers(map, seg->key))
    graphBoxDraw(box, BLACK, RED);
  else if (gMapIsNeighbour(map, seg->key))
    graphBoxDraw(box, BLACK, PALECYAN);
  else
    graphBoxDraw(box, BLACK, WHITE); 

  return;
} /* submapDoColour */


  
/* MODELS for COLUMN SUBMAP

   Submaps Loci_Query UNIQUE Text
           Interval_Query UNIQUE Text
           Show UNIQUE Int
	   Tolerance UNIQUE Int
	   Pos_By_Loci
	   Show_All_Markers
	   Scale UNIQUE float
*/

static void submapSave(COLINSTANCE instance, OBJ init)
{
  SUBPRIV private = instance->private;

  if(private->lociQuery)
    bsAddData(init, str2tag("Loci_Query"), _Text, &private->lociQuery);

  if(private->intQuery)
    bsAddData(init, str2tag("Interval_Query"), _Text, &private->intQuery);

  if(private->show)
    bsAddData(init, str2tag("Show"), _Int, &private->show);

  if(private->tolerance)
    bsGetData(init, str2tag("Tolerance"), _Int, &private->tolerance);
  
  if(private->posByLoci)
     bsAddTag(init, str2tag("pos_By_Loci")); 

  if(private->showAllMarkers)
     bsAddTag(init, str2tag("Show_All_Markers")); 

  if(private->scale)
    bsGetData(init, str2tag("Scale"), _Float, &private->scale);

  return;
} /* submapSave */

struct configLocals
{ char lociQuery[280],intQuery[280];
  int show;     /* 1 ALL 2 rejected 3 Non rejected */
  int tolerance;
  BOOL posByLoci;
  BOOL showAllMarkers;
  float scale;
};


static BOOL submapConfigure(COLINSTANCE instance)
{
  SUBPRIV private = instance->private;
  struct configLocals *cf =(struct configLocals *) messalloc(sizeof(struct configLocals));
  float y;
  int tag;

  if(controlCreateConfig(instance, cf, "Configure submap column",0.8,0.4))
    { /*initialise data */
      if(private->lociQuery)
	strcpy(cf->lociQuery, private->lociQuery);
      if(private->intQuery)
	strcpy(cf->intQuery, private->intQuery);
      cf->show = private->show;
      cf->tolerance = private->tolerance;
      cf->posByLoci = private->posByLoci;

      y = 1.0;

      graphTextEditor("Loci Query:", cf->lociQuery, 280, 4.0, y+=1.2,0);
      graphTextEditor("Interval Query:", cf->intQuery, 280, 4.0, y+=1.2,0);

      tag = graphRadioCreate("Show:", &cf->show, 4.0, y+=1.2);
      graphAddRadioEditor("All",tag, 1, 12.0, y);
      graphAddRadioEditor("Rejected", tag, 2, 20.0, y);
      graphAddRadioEditor("Not Rejected", tag, 3, 31.0, y);
      graphSetRadioEditor(tag, cf->show);

      graphIntEditor("Tolerance:", &cf->tolerance, 4.0, y+=1.2, 0);
      graphToggleEditor("Position Intervals by Loci", &cf->posByLoci, 4.0, y+=1.2);
      graphToggleEditor("Show all Markers", &cf->showAllMarkers, 4.0, y+=1.2);
      graphFloatEditor("Scale:", &cf->scale, 4.0, y+=1.2, 0);
      graphRedraw();
    }

  return FALSE;
} /* submapConfigure */


static void submapPointFinal(COLINSTANCE instance, void *locals, BOOL ok)
{
  struct configLocals *cf = locals;
  SUBPRIV private = instance->private;

  if(ok)
    { 
      private->show = cf->show;

      if(private->tolerance != cf->tolerance)
	{ private->tolerance = cf->tolerance;
	  calcRejects(instance);
	}

      if(private->intQuery)
	messfree(private->intQuery);
      if(private->intervalKeyset)
	keySetDestroy(private->intervalKeyset);
      
      if (strlen(cf->intQuery) != 0)
	{ private->intQuery = handleAlloc(0, 
				   instance->handle, 
				   1+strlen(cf->intQuery));
          strcpy(private->intQuery, cf->intQuery);
	  private->intervalKeyset = queryKey(instance->map->key, private->intQuery);
	}
      else
	{ private->intQuery = 0;
          private->intervalKeyset = 0;
	}

      if(private->lociQuery)
	messfree(private->lociQuery);
      if(private->lociKeyset)
	keySetDestroy(private->lociKeyset);
      
      if (strlen(cf->lociQuery) != 0)
	{ private->lociQuery = handleAlloc(0, 
				   instance->handle, 
				   1+strlen(cf->lociQuery));
          strcpy(private->lociQuery, cf->lociQuery);
	  private->lociKeyset = queryKey(instance->map->key, private->lociQuery);
	}
      else
	{ private->lociQuery = 0;
          private->lociKeyset = 0;
	}
      
      private->posByLoci = cf->posByLoci;
      private->showAllMarkers = cf->showAllMarkers;
      private->scale = cf->scale;

      calcNewPos(instance);
    }
  else
    messfree(cf);

  return;
} /* submapPointFinal */


static void subPrivDestroy (void *p)
     /* BlockFinalisation function for SUBPRIV-type */
{
  SUBPRIV private = (SUBPRIV)p;
  GMAPSUBMAPSEG *seg;
  int i;

  for(i = 0; i < arrayMax(private->segs); i++)
    { seg = arrp(private->segs, i, GMAPSUBMAPSEG);
      if(seg->submapmarkers)
	arrayDestroy(seg->submapmarkers);
      if(seg->negMarkers)
	keySetDestroy(seg->negMarkers);
      if(seg->negMarkersbyLoci)
	keySetDestroy(seg->negMarkersbyLoci);
    }
  arrayDestroy(private->segs);

  return;
} /* subPrivDestroy */


static void submapFollowBox(COLINSTANCE instance, int box, double x, double y)
{ MAPCONTROL map = instance->map;
  SUBPRIV private = instance->private;
  BOX2KEY *seg = arrp(private->box2key, box, BOX2KEY);
 
  BOOL isMap = (class(seg->key) == _VMap);

  if(isMap)
    { displayPreserve();
      display(seg->key, map->key, 0);
    }
  else
    display(seg->key, map->key, TREE);
}
  

static BOOL subMapCreate(COLINSTANCE instance, OBJ init)
{
  char *s1;
  int i;
  float f;
  SUBPRIV private;
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
    
  private = (SUBPRIV)halloc(sizeof(struct SubPrivStruct),
			    instance->handle);
  blockSetFinalise (private, subPrivDestroy);

  instance->draw = submapDraw;
  instance->setSelectBox = submapSetSelect;
  instance->unSelectBox = submapUnSelect;
  instance->save = submapSave;
  instance->doColour = submapDoColour;
  instance->followBox = submapFollowBox;
  instance->configFinal = submapPointFinal;
  instance->configure = submapConfigure;

  instance->private = private;
  private->lociQuery = 0;
  private->intQuery = 0;
  private->lociKeyset = private->intervalKeyset = 0;
  private->show = 1;
  private->tolerance = 7;
  private->posByLoci = TRUE;
  private->showAllMarkers = TRUE;
  private->scale = 0.0;
  private->currmap = 0;

  control->activeInstance = instance;

  private->segs = submapConvert(map);
  
  /*read in query's */
 
  if(init)
    {
      if(bsGetData(init, str2tag("Loci_Query"), _Text, &s1))
	{ private->lociQuery =strnew(s1,instance->handle);
	  private->lociKeyset = queryKey(instance->map->key, private->lociQuery);
       }
      
      if(bsGetData(init, str2tag("Interval_Query"), _Text, &s1))
	{ private->intQuery = strnew(s1,instance->handle);
	  private->intervalKeyset = queryKey(instance->map->key, private->intQuery);
       }
      
      if(bsGetData(init, str2tag("Show"), _Int, &i))
	private->show = i;
      
      if(bsGetData(init, str2tag("Tolerance"), _Int, &i))
	private->tolerance = i;
      
      if(bsFindTag(init, str2tag("pos_By_Loci")))
	private->posByLoci = TRUE;
      else
	private->posByLoci = FALSE;

      if(bsFindTag(init, str2tag("Show_All_Markers")))
	private->showAllMarkers = TRUE;
      else
	private->showAllMarkers = FALSE;

      if(bsGetData(init, str2tag("Scale"), _Float, &f))
	private->scale = f;
      
    }

  calcNewPos(instance);
  return TRUE;
} /* subMapCreate */


struct ProtoStruct gMapSubMapColumn = {
  0,
  subMapCreate,
  0,
  "Submaps",
  0,
  FALSE,
  submapConvert2,
  0,
  "The submap column can display any maps contained in the displayed "
  "object with its global loci. "
  "Global loci are described as having positions in both the submap "
  "and the object being displayed. "
  "The map set is then filtered using the Interval Query in the configuration "
  "to give the set of Map objects to draw. The loci set is filtered "
  "using the Loci Query to give the loci to use and draw.\n"
  "Submap intervals are drawn as vertical boxes. From these boxes are "
  "drawn lines to the loci which are displayed as text\n."
  "From the configuration you can position the submaps from either "
  "the position listed in the object or from the mean position of "
  "the loci (the filtered set) it contains. The tolerance is used to "
  "calculate wether to reject submaps based on wether all the loci "
  "lie with the tolerance of the start and end of the submap.\n"
  "You can choose to select which submap to draw based on this "
  "rejection by selecting Show to be either All, Rejected or NOT "
  "Rejected.\n"
  "Upon clicking a submap the loci corresponding to this item "
  "are highligted blue if they are within the tolerance of the "
  "start and end of the submap else they are highlighted red. "
  "By clicking on a loci the submaps which it is contained in "
  "are highlighted."
};


/************************ eof *******************************/
