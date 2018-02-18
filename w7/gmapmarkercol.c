/*  File: gmapmarkercol.c
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: marker display columns for the genetic map
 * Exported functions:
 *     gMapMainMarkersColumn
 *     gMapMiniChromBandsColumn
 * HISTORY:
 * Last edited: Mar  5 01:46 1999 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: gmapmarkercol.c,v 1.1.1.1 2002/07/19 20:23:00 sienkiew Exp $ */

/* NB uses GMAPSEG format, and gMapOrder */
#include "gmap.h"

static void *gMapMarkerConvert (MAPCONTROL map, void *params)
{
  KEY chrom = map->key;
  KEY locus;
  OBJ Chrom, Locus;
  Array segs, loci;
  GMAPSEG *seg;
  float x, dx;
  int i, flags;

  segs = arrayHandleCreate (20,GMAPSEG, map->handle) ;

/* If we inherit from another map, we look there for markers */
  while ((Chrom = bsCreate(chrom)))
    { if (bsFindTag (Chrom,_Main_Marker))
	goto gotit;
      if (!bsGetKey(Chrom, _From_map, &chrom))
	break; /* can't find markers, can't go back */
      bsDestroy(Chrom);
    }
  bsDestroy(Chrom);
  return segs; /* failed to open map, return empty array */

 gotit:
  loci = arrayCreate(100, BSunit);
      
  if (bsFindTag (Chrom,_Main_Marker) && bsFlatten(Chrom, 2, loci))
    for (i = 1; i<arrayMax(loci); i += 2)
      { locus = arr(loci, i, BSunit).k;
	flags = gMapGetMapObject(locus, chrom, 0, 20, &x, &dx, &Locus, 0);
	if (!flags)
	  continue;
	
	
	if (bsFindTag (Locus,_Centromere))
	  flags |= FLAG_CENTROMERE ;
	if (bsFindTag (Locus,_Dark))
	  flags |= FLAG_DARK_BAND ;
	if (bsFindTag (Locus,_NOR))
	  flags |= FLAG_NOR ;
	if (bsFindTag (Locus,_p_Telomere))
	  flags |= FLAG_P_TELOMERE;
	if (bsFindTag (Locus,_q_Telomere))
	  flags |= FLAG_Q_TELOMERE;
	
	seg = arrayp(segs,arrayMax(segs),GMAPSEG) ;
	seg->key = locus;
	seg->x = x;
	seg->dx = dx;
	seg->flag = flags;
	
	bsDestroy (Locus) ;
      }
  
  bsDestroy (Chrom) ;
  arrayDestroy(loci);
  arraySort (segs, gMapOrder) ; 

  return (void *)segs ;
}

/***************************/
/* marker chromosome bands */
/***************************/
static void miniChromBandsDraw (COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  float y, dy, x = *offset + 0.5 , l, max = 0 , ylast = - 1000 ;
  int i ;
  GMAPSEG* seg ;
  BOOL isBands = FALSE ;
  Array segs = instance->private;
  int firstSeg, lastSeg, countDir;

  if (map->mag > 0)
    { firstSeg = 0;
      lastSeg = arrayMax(segs);
      countDir = 1;
    }
  else
    { firstSeg = arrayMax(segs)-1 ;
      lastSeg = -1;
      countDir = -1;
    }
  
  graphTextHeight (0.75) ;
  for (i = firstSeg; i != lastSeg; i += countDir)
    { seg = arrp(segs,i,GMAPSEG) ;
      if (seg->flag & FLAG_ANY_INTERVAL)
	{ isBands = TRUE ;
	  y = MAP2WHOLE(map, seg->x) ;
	  dy =  seg->dx *map->thumb.fac ;
	  
	  if (seg->flag & FLAG_DARK_BAND)
	    graphColor (DARKGRAY) ;
	  else if (seg->flag & FLAG_NOR)
	    graphColor (LIGHTGRAY) ;
	  else
	    graphColor (WHITE) ;
	  
	  if (seg->flag & FLAG_CENTROMERE)
	    { graphFillRectangle (x, y-dy, x + 2., y+dy) ;
	      graphColor (BLACK) ;
	      graphLine (x, y-dy, x + 2., y+dy) ;
	      graphLine (x + 2., y-dy, x, y+dy) ;
	      graphLine (x, y-dy, x + 2., y-dy) ;
	      graphLine (x + 2., y+dy, x, y+dy) ;
	    }	      
	  else if ((seg->flag & FLAG_Q_TELOMERE || 
		    seg->flag & FLAG_P_TELOMERE) &&
		   (dy>1.0 || dy<(-1.0)))
	    { /* dy<0 if flipped */
	      /* P is at the top unless flipped */
	      BOOL up;
	      if (dy>0)
		{ up = seg->flag & FLAG_P_TELOMERE;
		}
	      else
		{ up = seg->flag & FLAG_Q_TELOMERE;
		  dy = -dy;
		}
	      if (up)
		{ graphFillRectangle(x, y+1-dy, x+2, y+dy);
		  graphFillArc(x+1, y+1-dy, 1, 0, 180);
		  graphColor(BLACK);
		  graphArc(x+1, y+1-dy, 1, 0, 180);
		  graphLine(x + 2., y+dy, x, y+dy);
		  graphLine(x, y+dy, x, y+1-dy);
		  graphLine(x+2., y+dy, x+2., y+1-dy);
		}
	      else /* down */
		{ graphFillRectangle(x, y-dy, x+2, y+dy-1);
		  graphFillArc(x+1, y+dy-1, 1, 180, 180);
		  graphColor(BLACK);
		  graphArc(x+1, y+dy-1, 1, 180, 180);
		  graphLine(x + 2., y-dy, x, y-dy);
		  graphLine(x, y-dy, x, y+dy-1);
		  graphLine(x+2., y-dy, x+2., y+dy-1);
		}
	    }	      
	  else
	    { graphFillRectangle (x, y-dy, x + 2., y+dy) ;
	      graphColor (BLACK) ;
	      graphRectangle (x, y-dy, x + 2., y+dy) ;
	    }

	  if (y > ylast + 1)
	    { graphText (name(seg->key),x+2.5,y-0.25) ;
	      l = .65 * strlen(name(seg->key)) ;
	      if (l > max)
		max = l ;
	      ylast = y ;
	    }
	}
    }
  
  graphTextHeight (0) ;
  if (isBands)
    *offset += 3.5 + max ;
}

static BOOL miniChromBandsCreate(COLINSTANCE instance, OBJ init)
{ 

  instance->draw = miniChromBandsDraw;
  instance->private = *(instance->proto->convertResults);
  return TRUE;
}

struct ProtoStruct gMapMiniChromBandsColumn = {
  0,
  miniChromBandsCreate,
  0,
  "Marker_intervals",
  0,
  FALSE,
  gMapMarkerConvert,
  0,
  "This column displays intervals as chromosome bands. The whole extent of "
    "the map is always shown, making it useful for navigation around the map. "
      "All interval objects two to the right of the tag Main_markers are "
	"drawn. The style of drawing is affected by the presence of the "
	  "following tags in the interval objects.\n\n"
	    "Dark\nNOR\nCentromere\np_telomere\nq_telomere\n"
};

/****************************************************************************/
/*                            Marker Loci                                   */
/****************************************************************************/

static void mainMarkersDraw(COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  float y, xx, max = 0 ;
  int i, box ;
  GMAPSEG* seg ;
  Array segs = instance->private;

  graphTextHeight (0.75) ;

  box = graphBoxStart () ;
  for (i = 0 ; i < arrayMax(segs) ; ++i)
    { seg = arrp(segs,i,GMAPSEG) ;
      if (seg->flag & FLAG_ANY_LOCUS)
	{ y = MAP2WHOLE(map, seg->x) ;
	  graphText (name(seg->key), *offset+0.5, y-0.25) ;
	  xx = (isGifDisplay ? 1 : .65) * strlen(name(seg->key)) ;
	  if (xx > max)
	    max = xx ;
	}
    }
				/* draw lines to ensure box is full size */
  y = MAP2WHOLE(map,map->min) ;
  graphLine (*offset, y, *offset + max/2, y) ;
  y = MAP2WHOLE(map,map->max) ;
  graphLine (*offset, y, *offset + max/2, y) ;
  graphBoxEnd() ;
  array(instance->map->control->boxIndex, box, COLINSTANCE) = instance;

  graphTextHeight (0) ;
  *offset += max + 1.0 ;
}

static void mainMarkerPick(COLINSTANCE instance,
			   int box, double x, double y)
{
  MAPCONTROL map = instance->map ;
  map->centre = WHOLE2MAP(map, y + MAP2WHOLE(map,map->min)) ;
  controlDraw() ;
} /* mainMarkerPick */

static BOOL mainMarkerCreate(COLINSTANCE instance, OBJ init)
{ 
  instance->draw = mainMarkersDraw;
  instance->private = *(instance->proto->convertResults);
  instance->pick = mainMarkerPick ;
  return TRUE;
}

struct ProtoStruct gMapMainMarkersColumn = {
  0,
  mainMarkerCreate,
  0,
  "Marker_points",
  0,
  FALSE,
  gMapMarkerConvert,
  0,
  "This column displays points. The whole extent of "
    "the map is always shown, making it useful for navigation around the map. "
      "All point  objects two to the right of the tag Main_markers are "
	"drawn."
};

