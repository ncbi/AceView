/*  File: vmapphys.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 	vMapAddPhysGenes, vMapDrawPhysGenes, vMapDrawContigs,
	vMapReversedPhysical, vMapToPMap, pMapToVMap
 * HISTORY:
 * Last edited: Jan  7 17:51 1999 (fw)
 * Created: Sat May 29 18:52:31 1993 (cgc)
 *-------------------------------------------------------------------
 */

/* $Id: vmapphys.c,v 1.4 2014/11/30 03:20:49 mieg Exp $ */

#include "acedb.h"
#include "bump.h"

#include "display.h"
#include "systags.h"
#include "tags.h"
#include "classes.h"
#include "sysclass.h"
#include "lex.h"
#include "a.h"
#include "bs.h"
#include "vmap_.h"

/************************************************************/

#define myMapIntervalMenu vMapIntervalMenu

typedef struct
  { KEY key , symbol ;
    float x, dx ;
    unsigned int flag ;
  } MY_VMAPSEG ;
#define MY_vmapFormat "kkffi"
#define MY_VMap _VvMap

typedef struct
  { KEY key;
    KEY parent; /* clone associated with this item */
    float x0, x1;
    unsigned int flag;
  } PSEG ;
#define psegFormat "kkffi"
#define FMINUSINF (-1000000000.0)

/***********************************************/

void vMapAddPhysGenes (Array gsegs)
{
  int	i, j, k ;
  KEY   pmap ;
  Array psegs ;
  PSEG  *pseg=0, *qseg=0 ;
  SEG   *gseg, *seg ;
  float y1, y2, y, x1, x2 ;
  static Associator gene2pos = 0 ;
  union { float f ; void* v ; } convert ;

  gene2pos = assReCreate (gene2pos) ;
  for (i = 0 ; i < arrayMax(gsegs)-1 ; ++i)
    { gseg = arrp(gsegs,i,SEG) ;
      if (gseg->flag & FLAG_CLONED)
	{ convert.f = gseg->x ;
	  if (!gseg->x)
	    convert.f += 0.000001 ;
	  assInsert (gene2pos, assVoid (gseg->key), convert.v) ;
	}
    }

  for (i = 0 ; i < arrayMax(gsegs) ; ++i)
    { gseg = arrp(gsegs,i,SEG) ;
/*
      if (class(gseg->key) == _VContig)
	printf ("%s with flag %x\n", name(gseg->key), gseg->flag) ;
*/
      if (class(gseg->key) == _VContig &&
	  !(gseg->flag & FLAG_HIDE) && 
	  lexReClass (gseg->key, &pmap, _VpMap) &&
	  (psegs = arrayGet (pmap, PSEG, psegFormat)))
	{ for (j = 0 ; j < arrayMax(psegs) ; ++j)
	    { pseg = arrp(psegs,j,PSEG) ;
	      if (pseg->x0 > FMINUSINF)
		break ;
	    }
	  y1 = y2 = gseg->x ;
	  x1 = x2 = pseg->x0 ; x2 += 0.00001 ; /* for safety */
	  for (j = 1 ; j < arrayMax(psegs) ; ++j)
	    { pseg = arrp(psegs, j, PSEG) ;
	      if (class(pseg->key) == _VLocus)
		if (!assFind (gene2pos, assVoid(pseg->key), &y))
		  { if (pseg->x0 >= x2)
		      { y1 = y2 ;
			x1 = x2 ;
			for (k = j ; k + 1 < arrayMax(psegs) ; ++k)
			  { qseg = arrp(psegs, k, PSEG) ;
			    if (class(qseg->key) == _VLocus &&
				qseg->x0 > x1 &&
				assFind (gene2pos, 
					 assVoid(qseg->key), &y2))
			      break ;
			  }
			x2 = qseg->x0 ;
			if (k >= arrayMax(psegs)-1)
			  y2 = gseg->x + gseg->dx ; /* contig end */
		      }
		    seg = arrayp(gsegs,arrayMax(gsegs),SEG) ;
		    seg->key = pseg->key ;
		    if (x2 != x1)
		      seg->x = y1 + (y2-y1) * (pseg->x0-x1) / (x2-x1) ;
		    else
		      seg->x = y1 ;
		    seg->dx = 0 ;
		    seg->flag = FLAG_PHYS_GENE | FLAG_CLONED ;
		    convert.f = seg->x ;
		    assInsert (gene2pos, assVoid(seg->key), convert.v) ;
		  }
	    }
	  arrayDestroy (psegs) ;
	}
    }
  
  return;
} /* vMapAddPhysGenes */

/************************************************/

void vMapDrawPhysGenes (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  BUMP  bump = bumpCreate (12, 0) ;
  float y ;
  int   i, ibox, ix ;
  SEG   *seg ;
  int   hideHeader = (look->flag & FLAG_HIDE_HEADER) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      y = MAP2GRAPH(look->map,seg->x) ;
      if (!(seg->flag & FLAG_PHYS_GENE) ||
	  seg->flag & FLAG_HIDE ||
	  y < topMargin + 1 ||
	  y > mapGraphHeight - 1)
	continue ;
      array (look->boxIndex, ibox=graphBoxStart(), int) = i ;
      ix = 0 ;
      bumpItem (bump, strlen(name(seg->key))+1, 1.0, &ix, &y) ;
      graphText (name(seg->key), *offset + ix, y) ;
      graphBoxEnd () ;
      graphBoxDraw (ibox, BLACK, hideHeader ? WHITE : YELLOW) ;
    }

  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;

  return;
} /* vMapDrawPhysGenes */

/******************************************************/

void vMapDrawContigs (LOOK genericLook, float *offset)
{				/* based on drawRearrangements */
  VerticalMap look = (VerticalMap)genericLook;
  BUMP  bump = bumpCreate (mapGraphWidth, 0) ;
  float y1, y2, x, y ;
  int   i, j, ibox, ix ;
  SEG   *seg ;

  *offset += 1 ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (class(seg->key) != _VContig ||
	  seg->flag & FLAG_HIDE)
	continue ;
      y1 = MAP2GRAPH(look->map,seg->x) ;
      y2 = MAP2GRAPH(look->map,seg->x + seg->dx) ;
      if (y2 > topMargin+1 && y1 < mapGraphHeight-1)
	{ array(look->boxIndex,ibox=graphBoxStart(),int) = i ;
	  ix = 0 ; 
	  bumpItem (bump,1,(y2-y1)+0.2,&ix,&y1) ;
	  x = *offset + ix ;
	  if (y1 < topMargin + 1) 
	    y1 = topMargin + 1;
	  if (y2 > mapGraphHeight)
	    y2 = mapGraphHeight ;
	  graphLine (x-0.5, y1, x-0.5, y2) ;
	  graphLine (x+0.5, y1, x+0.5, y2) ;

	  for (j = i ; j < arrayMax(look->segs) ; ++j)
	    { seg = arrp(look->segs,j,SEG) ;
	      if (class(seg->key) != _VClone)
		continue ;
	      y = MAP2GRAPH(look->map, seg->x) ;
	      if (y < y1 || seg->flag & FLAG_HIDE)
		continue ;
	      if (y > y2)
		break ;
	      array(look->boxIndex,graphBoxStart(),int) = j ;
	      graphRectangle (x-0.5, y-0.1, x+0.5, y+0.1) ;
	      graphBoxEnd() ;
	    }
	  graphBoxEnd () ;
	  graphBoxDraw (ibox, BLACK, YELLOW) ;
	  graphBoxMenu (ibox, myMapIntervalMenu) ;
	}
    }

  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;

  return;
} /* vMapDrawContigs */

/*****************************************************/

void vMapReversedPhysical (LOOK genericLook, float *offset)
{
  VerticalMap look = (VerticalMap)genericLook;
  float y, oldy=0 ;
  int   i, x, oldx ;
  KEY   contig = 0, new ;
  SEG   *seg ;
  OBJ	obj ;

  graphColor (RED) ;

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (class(seg->key) != _VClone ||
	  !(obj = bsCreate (seg->key)))
	continue ;
      y = MAP2GRAPH(look->map,seg->x) ;
      if (bsFindKey (obj, _pMap, contig) &&
	  bsGetData (obj, _bsRight, _Int, &x))
	{ if (x < oldx && y > oldy)
	    graphFillRectangle (*offset, oldy, *offset+0.5, y) ;
	  oldx = x ;
	  oldy = y ;
	}
      else if (bsGetKey (obj, _pMap, &new) &&
	       bsGetData (obj, _bsRight, _Int, &oldx))
	{ oldy = y ;
	  contig = new ;
	}
      bsDestroy (obj) ;
      if (oldy < topMargin + 1)
	oldy = topMargin + 1 ;
      else if (oldy  > mapGraphHeight - 1)
	break ;
    }

  graphColor (BLACK) ;

  *offset += 0.5 ;

  return;
} /* vMapReversedPhysical */

/*****************************************************/

void vMapToPMap (VerticalMap look, float y)
{
  int i, p1, p2, c1, c2 ;
  float x1, x2, y1, y2 ;
  SEG *seg = arrp(look->segs, 
		  arr(look->boxIndex,look->activeBox,int),
		  SEG) ;
  KEY contig = seg->key, from ;
  OBJ obj ;

  graphBoxDim (look->activeBox, &x1, &y1, &x2, &y2) ;
  y = GRAPH2MAP(look->map, y+y1) ;

  if ((obj = bsCreate (contig)) &&
      bsGetData (obj, _pMap, _Int, &p1) &&
      bsGetData (obj, _bsRight, _Int, &p2))
    { y1 = seg->x ;
      y2 = seg->x + seg->dx ;
    }
  else
    { y1 = -1000000 ;
      y2 = 1000000 ;
    }
  bsDestroy (obj) ;		/* OK if 0 */

  for (i = 0 ; i < arrayMax(look->segs) ; ++i)
    { seg = arrp(look->segs,i,SEG) ;
      if (class(seg->key) == _VClone &&
	  seg->x > y1 && seg->x < y2 &&
	  (obj = bsCreate (seg->key)) && 
	  bsFindKey (obj, _pMap, contig) && 
	  bsGetData (obj, _bsRight, _Int, &c1) &&
	  bsGetData (obj, _bsRight, _Int, &c2))
	{
	  if (seg->x < y)
	    { y1 = seg->x ; p1 = 0.5*(c1+c2) ; }
	  else
	    { y2 = seg->x ; p2 = 0.5*(c1+c2) ; }
	}
      bsDestroy (obj) ;
    }

  if (y1 > -999999 && y2 < 999999)
    { i = 0x800000L + p1 + (p2 - p1) * (y - y1) / (y2 - y1) ;
      from = KEYMAKE (_VCalcul,i) ;
      display (contig, from, 0) ;
    }
  else
    messout ("Sorry, I can't find flanking points to "
	     "interpolate from.") ;

  return;
} /* vMapToPMap */

/******************************/

void pMapToVMap (KEY contig, KEY from, int x)
{
  KEY vmap, map ;
  OBJ obj ;
  SEG *seg ;
  float y1=0, y2=0, y ;
  int i, p1 = -1000000, p2 = 1000000, x1, x2 ;
  Array segs ;
  BOOL inContig = FALSE ;

  if (!(obj = bsCreate(contig)))
    { messout ("Can't open contig object") ;
      return ;
    }
  if (!bsGetKey (obj, _Map, &map) || !map)
    { bsDestroy (obj) ;
      messout ("This contig is not assigned to a chromosome") ;
      return ;
    }
  bsGetData (obj, _pMap, _Int, &p1) ;
  bsGetData (obj, _bsRight, _Int, &p2) ;
  bsDestroy (obj) ;

		/* next get the segs */
  if (!lexReClass (map, &vmap, MY_VMap) ||
      !(segs = arrayGet (vmap, SEG,MY_vmapFormat)))
    { messout ("Please recalculate the genetic map of %s",
	       name(map)) ;
      return ;
    }
  
  if (from)			/* look for it in the segs */
    for (i = arrayMax(segs) ; i-- ;)
      if (arr(segs,i,SEG).key == from)
	{ display (vmap, from, 0) ;
	  arrayDestroy (segs) ;
	  return ;
	}

 /* strategy: find 2 vmapped clones around or by x and interpolate */

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
    messout ("Sorry, I could not find two flanking markers") ;
  else
    { y = y1 + (y2 - y1) * (x - p1) / (p2 - p1) ;
	from = KEYMAKE(_VCalcul, 1000.0*(y + 1000.0)) ;
      display (vmap, from, 0) ;
    }

  return;
} /* pMapToVMap */

/************ end of file *************/
