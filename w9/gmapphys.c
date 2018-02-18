/*  File: gmapphys.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 *	pMapToGMap
 *      gMapGetClones
 *      gMapPhysGenesColumn
 *      gMapContigColumn
 *      gMapRevPhysColumn
 * HISTORY:
 * Last edited: Dec 21 14:31 1998 (fw)
 * Created: Sat May 29 18:52:31 1993 (cgc)
 *-------------------------------------------------------------------
 */

/* $Id: gmapphys.c,v 1.5 2014/11/30 03:20:49 mieg Exp $ */

#include "gmap.h"
#include "freeout.h"
#include "bump.h"

/***********************************************/

typedef struct 
{ KEY key;
  float x, dx;
} GMAPPHYSEG;

typedef struct
  { KEY key;
    KEY parent; /* clone associated with this item */
    float x0, x1;
    unsigned int flag;
  } PSEG ;
#define pMapFormat "kkffi"
#define FMINUSINF (-1000000000.0)

/***********************************************/

static int gMapPhysOrder (const void *a, const void *b)  /* for arraySort() call */
{
  const GMAPPHYSEG *seg1 = (const GMAPPHYSEG*)a, *seg2 = (const GMAPPHYSEG*)b ;
  float diff, x1, x2;

  if (class(seg1->key) == _VContig)
    x1 = seg1->x - seg1->dx; /* sort on start posn of intervals */
  else
    x1 = seg1->x;

  if (class(seg2->key) == _VContig)
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
} /* gMapPhysOrder */


static void *gMapPhysConvert (MAPCONTROL map, void *params)
/* This does conversion for the physical genes and contig columns */
{
  GeneticMap look  = (GeneticMap)(map->look);
  int	i, j, k ;
  KEY   pmap ;
  Array psegs = 0 ;
  PSEG  *pseg=0, *qseg=0 ;
  GMAPPHYSEG   *seg ;
  float y1, y2, x1, x2, x, dx;
  Associator done;
  KEY chrom, clone, locus, cache;
  OBJ Locus= 0 ;
  Array segs = 0, loci = 0 ;

  getPos(map, 0, 0); /* make sure look->orderedLoci is valid */

  /* are they cached */
  if ((cache = gMapCacheKey(map->key, "phys")) &&
      (segs = arrayGet(cache, GMAPPHYSEG, "kff")))
    return segs;
  
  segs = arrayHandleCreate (128,GMAPPHYSEG,map->handle) ;
  done = assCreate();
  loci = gMapMapObjects(map->key);

  for (i = 0; i <arrayMax(loci); i++)
    { locus = arr(loci, i, KEYPAIR).obj;
      chrom = arr(loci, i, KEYPAIR).map;
      if (!gMapGetMapObject(locus, map->key, chrom, 20, &x, &dx, &Locus, 0)){
	bsDestroy(Locus);
	continue;
      }
      if (bsGetKey (Locus,_Positive_clone,&clone))
	{ 
	  seg = arrayp(segs, arrayMax(segs), GMAPPHYSEG) ;
	  seg->key = clone ;
	  seg->x = x ;
	  seg->dx = 0 ;
	}
      
      bsDestroy(Locus);
      
      if (class(locus) == _VContig)
	{ seg = arrayp(segs, arrayMax(segs), GMAPPHYSEG);
	  seg->key = locus;
	  seg->x = x;
	  seg->dx = dx;
	  
	  
	  /* now see if there's a pMap with the same name as the contig, and
	     interpolate genes from it */
	  if (lexReClass (locus, &pmap, _VpMap) &&
	      (psegs = arrayGet (pmap, PSEG, pMapFormat)))
	    { for (j = 0 ; j < arrayMax(psegs) ; ++j)
		{ pseg = arrp(psegs,j,PSEG) ;
		  if (pseg->x0 > FMINUSINF)
		    break ;
		}
	      y1 = y2 = x - dx; /* Contig start */
	      x1 = x2 = pseg->x0 ; x2 += 0.00001 ; /* for safety */
	      for (j = 1 ; j < arrayMax(psegs) ; ++j)
		{ pseg = arrp(psegs, j, PSEG) ;
		  if (class(pseg->key) == _VLocus)
		    if (!keySetFind(look->orderedLoci, pseg->key, 0) &&
			!assFind(done, assVoid(pseg->key), 0))
		      { if (pseg->x0 >= x2)
			  { y1 = y2 ;
			    x1 = x2 ;
			    for (k = j ; k + 1 < arrayMax(psegs) ; ++k)
			      { qseg = arrp(psegs, k, PSEG) ;
				if (class(qseg->key) == _VLocus &&
				    qseg->x0 > x1 &&
				    getPos (map, qseg->key, &y2))
				  break ;
			      }
			    x2 = qseg->x0 ;
			    if (k >= arrayMax(psegs)-1)
			      y2 = x + dx ; /* contig end */
			  }
			seg = arrayp(segs,arrayMax(segs), GMAPPHYSEG) ;
			seg->key = pseg->key ;
			if (x2 != x1)
			  seg->x = y1 + (y2-y1) * (pseg->x0-x1) / (x2-x1) ;
			else
			  seg->x = y1 ;
			seg->dx = 0 ;
			assInsert (done, assVoid(seg->key), 0) ;
		      }
		}
	      arrayDestroy (psegs) ;
	    }
	}
    }
  
  
  arrayDestroy(loci);
  assDestroy(done);
  
  arraySort(segs, gMapPhysOrder);
  gMapCacheStore(map->key, "phys", segs, "kff");
  return (void*) segs;
} /* gMapPhysConvert */


/*****************************************************************************/
/*           Physical Genes                                                  */
/*****************************************************************************/
static BOOL physGenesSetSelect(COLINSTANCE instance, 
			    int box,
			    double x,
			    double y)
{ 
  COLCONTROL control = instance->map->control;
  GeneticMap look = (GeneticMap)(instance->map->look);
  GMAPPHYSEG *seg = (GMAPPHYSEG *) arr(control->boxIndex2, box, void *); 

  look->selectKey = seg->key;
  gMapUnselect(instance->map);
  look->neighboursInfoValid = TRUE;
  (void)gMapNeighbours(instance->map, &look->neighbours, seg->key);

  return FALSE;
} /* physGenesSetSelect */


static void physGenesDoColour(COLINSTANCE instance, int box)
     /* NB this used by the rearrangements column also */
{ 
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);
  GMAPPHYSEG *seg = (GMAPPHYSEG *) arr(control->boxIndex2, box, void *); 
  
  if (seg->key == look->selectKey && instance == control->activeInstance)
    { graphBoxDraw(box, BLACK, CYAN);
      control->activeBox = box;
      strncpy(look->messageText, name(seg->key), 100);
      strcat(look->messageText, messprintf(" %.2f", seg->x));
      strcat(look->messageText, messprintf(" interpolated"));
      
      if (look->messageBox) 
	graphBoxDraw(look->messageBox, -1, -1);
    }
  else if (keySetFind(look->highlight, seg->key, 0))
    graphBoxDraw(box, BLACK, MAGENTA);
  else if (gMapIsNeighbour(map, seg->key))
    graphBoxDraw(box, BLACK, PALECYAN);   
  else
    graphBoxDraw(box, BLACK, WHITE);

  return;
} /* physGenesDoColour */


void physGenesSummary (void)
{
  MAPCONTROL map = currentMapControl() ;
  COLCONTROL control = map->control ;
  COLINSTANCE instance = 0 ;
  GMAPPHYSEG   *seg ;
  int i ;
  Array segs ;
  FILE *fil ;


  if (!(fil = filqueryopen (0, 0, "out", "a", "Output file to add physical genes to")))
    return ;

  for (i = 0 ; i < arrayMax (control->instances) ; ++i)
    { instance = arr(control->instances, i, COLINSTANCE) ;
      if (!strcmp (instance->name, "Physical_genes") &&
	  instance->map == map)
	break ;
    }
  if (i == arrayMax(control->instances))
    { messout ("Sorry, can't find a physical genes column!") ;
      return ;
    }

  segs = instance->private ;
  for (i = 0 ; i < arrayMax(segs) ; ++i)
    { seg = arrp(segs,i,GMAPPHYSEG) ;
      fprintf (fil, "%s %s %.2f\n",
	       name(seg->key), name(map->key), seg->x) ;
    }

  filclose (fil) ;

  return;
} /* physGenesSummary */


static void physGenesDraw (COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);
  BUMP  bump = bumpCreate (12, 0) ;
  float y ;
  int   i, ibox, ix , n ;
  GMAPPHYSEG   *seg ;
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
  
  for (i = firstSeg; i != lastSeg ; i += countDir)
    { seg = arrp(segs,i,GMAPPHYSEG) ;
      if (class(seg->key) != _VLocus || keySetFind(look->hidden, seg->key, 0))
	continue;
      y = MAP2GRAPH(map,seg->x) ;
      if (y < control->topMargin + 1 ||
	  y > control->graphHeight - 1)
	continue ;
      ibox = graphBoxStart();
      array (control->boxIndex, ibox, COLINSTANCE) = instance;
      array (control->boxIndex2, ibox, void *) = (void*) seg;
      ix = 0 ;
      /* was bumpItem (bump, strlen(name(seg->key))+1, 1.0, &ix, &y) ; */
      if ((n = bumpText (bump, name(seg->key), &ix, &y, 1, TRUE)) &&
	  y < control->graphHeight-1)
	{ char *ccp = name(seg->key), *ccq = ccp +n, ccc = *ccq ;
	  *ccq = 0 ; graphText (ccp, *offset+ix, y-0.5) ; *ccq = ccc ;
	  remarkRegister(map, seg->key, seg->x);
	}
      graphBoxEnd () ;
      graphBoxInfo (ibox, seg->key,0) ; 
      if (seg->key == control->from)
	control->fromBox = ibox;
    }

  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;

  return;
} /* physGenesDraw */


static BOOL physGenesUnselect(COLINSTANCE instance, 
			      int box)
{ 
  GeneticMap look = (GeneticMap)(instance->map->look);

  
  *look->messageText = 0;
  if (look->messageBox)
    graphBoxDraw(look->messageBox, -1, -1);

  return FALSE;
} /* physGenesUnselect */


static void physGenesFollowBox(COLINSTANCE instance, int box, 
			      double x, double y)
{
  MAPCONTROL map = instance->map;
  GMAPPHYSEG *seg = (GMAPPHYSEG *) arr(map->control->boxIndex2, box, void *); 
  
  display(seg->key, map->key, TREE);

  return;
} /* physGenesFollowBox */


static BOOL physGenesCreate(COLINSTANCE instance, OBJ init)
{ 
  instance->draw = physGenesDraw;
  instance->setSelectBox = physGenesSetSelect;
  instance->unSelectBox = physGenesUnselect;
  instance->doColour = physGenesDoColour;
  instance->followBox = physGenesFollowBox;
  instance->private = *(instance->proto->convertResults);

  return TRUE;
} /* physGenesCreate */

  
struct ProtoStruct gMapPhysGenesColumn = {
  0,
  physGenesCreate,
  0,
  "Physical_genes",
  0,
  FALSE,
  gMapPhysConvert,
  0
};

/****************************************************************************/
/*              Contigs                                                     */
/****************************************************************************/

static BOOL contigSetSelect(COLINSTANCE instance, int box, double x, double y)
{ 
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);
  GMAPPHYSEG *seg = (GMAPPHYSEG *) arr(control->boxIndex2, box, void *);

  look->selectKey = seg->key;
  gMapUnselect(map);

  if (class(seg->key) == _VContig)
    { /* only know neighbours of contigs, not clones */
      keySetDestroy(look->neighbours);
      look->neighbours = queryKey(map->key,
				  messprintf(">Locus ; "
					     ">Positive_clone pMap = %s; "
					     ">Positive_locus", 
					     name(seg->key)));
      look->neighboursInfoValid = TRUE;
    }

  return FALSE;
} /* contigSetSelect */


static BOOL contigUnselect(COLINSTANCE instance, 
			   int box)
{ 
  GeneticMap look = (GeneticMap)(instance->map->look);

  
  *look->messageText = 0;
  if (look->messageBox)
    graphBoxDraw(look->messageBox, -1, -1);

  return FALSE;
} /* contigUnselect */


static void contigDoColour(COLINSTANCE instance, int box)
{ 
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GMAPPHYSEG *seg = (GMAPPHYSEG *) arr(control->boxIndex2, box, void *);
  GeneticMap look = (GeneticMap)(map->look);

  if (look->selectKey == seg->key && instance == control->activeInstance)
    {
      control->activeBox = box;
      graphBoxDraw(box, BLACK, CYAN);

      strncpy(look->messageText, name(seg->key), 100);
      if (class(seg->key) == _VContig)
	strcat(look->messageText, 
	       messprintf(" %.2f %.2f", seg->x - seg->dx, seg->x + seg->dx));
      else
	strcat(look->messageText, 
	       messprintf(" %.2f", seg->x ));
      
      if (look->messageBox)
	graphBoxDraw(look->messageBox, -1, -1);
    }
  else
    {
      graphBoxDraw (box, BLACK, YELLOW) ;
    }

  return;
} /* contigDoColour */


static void contigFollowBox(COLINSTANCE instance, int box, double x, double y)
{
  MAPCONTROL map = instance->map;
  GMAPPHYSEG *seg = (GMAPPHYSEG *) arr(map->control->boxIndex2, box, void *); 
  Array segs = instance->private;

 if (class(seg->key) == _VClone) 
   display(seg->key, map->key, 0);
  
  if (class(seg->key) == _VContig)
    { /* Go to Pmap */
      int i, p1, p2, c1, c2 ;
      float x1, x2, y1, y2 ;
      KEY contig = seg->key, from ;
      OBJ obj ;
      
      graphBoxDim (box, &x1, &y1, &x2, &y2) ;
      y = GRAPH2MAP(map, y+y1) ;

      if ((obj = bsCreate (contig)) &&
	  bsGetData (obj, _pMap, _Int, &p1) &&
	  bsGetData (obj, _bsRight, _Int, &p2))
	{ y1 = seg->x - seg->dx ;
	  y2 = seg->x + seg->dx ;
	}
      else
	{ y1 = -1000000 ;
	  y2 = 1000000 ;
	}
      bsDestroy (obj) ;		/* OK if 0 */
      
      for (i = 0 ; i < arrayMax(segs) ; ++i)
	{ seg = arrp(segs,i,GMAPPHYSEG) ;
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
    }
  
  return;
} /* contigFollowBox */


static void contigDraw (COLINSTANCE instance, float *offset)
{				/* based on drawRearrangements */
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);
  BUMP  bump = bumpCreate (control->graphWidth, 0) ;
  float y1, y2, x, y ;
  int   i, j, ibox, box2, ix ;
  GMAPPHYSEG   *seg, *seg2 ;
  Array segs = instance->private;
  int firstSeg, lastSeg, countDir;

  if (map->mag > 0)
    { firstSeg = 0;
      lastSeg = arrayMax(segs);
      countDir = 1;
    }
  else
    { firstSeg = arrayMax(segs)-1;
      lastSeg = -1;
      countDir = -1;
    }

  *offset += 1 ;

  for (i = firstSeg ; i != lastSeg ; i += countDir)
    { seg = arrp(segs,i,GMAPPHYSEG) ;
      if (class(seg->key) != _VContig)
	continue ;
      if (keySetFind(look->hidden, seg->key, 0))
	continue;
      y1 = MAP2GRAPH(map,seg->x - seg->dx) ;
      y2 = MAP2GRAPH(map,seg->x + seg->dx) ;
      if (y1>y2)
	{ float tmp = y1; y1 = y2; y2 = tmp; }
      if (y2 > control->topMargin+1 && y1 < control->graphHeight-1)
	{ ibox = graphBoxStart();
	  array(control->boxIndex, ibox, COLINSTANCE) = instance ;
	  array(control->boxIndex2, ibox, void *) = (void*) seg;
	  ix = 0 ; 
	  if (y1 < control->topMargin) 
	    y1 = control->topMargin;
	  if (y2 > control->graphHeight)
	    y2 = control->graphHeight;
	  bumpItem (bump,1,(y2-y1)+0.2,&ix,&y1) ; 
	  x = *offset + ix ;
	  graphLine (x-0.5, y1, x-0.5, y2) ;
	  graphLine (x+0.5, y1, x+0.5, y2) ;
	  remarkRegister(map, seg->key, seg->x);

	  for (j = firstSeg ; j != lastSeg ; j += countDir)
	    { seg2 = arrp(segs,j,GMAPPHYSEG) ;
	      if (class(seg2->key) != _VClone)
		continue ;
	      if (keySetFind(look->hidden, seg2->key, 0))
		continue;
	      y = MAP2GRAPH(map, seg2->x) ;
	      if (y < y1)
		continue ;
	      if (y > y2)
		break ;
	      box2 = graphBoxStart();
	      array(control->boxIndex, box2, COLINSTANCE) = instance ;
	      array(control->boxIndex2, box2, void *) = (void *) seg2;
	      graphRectangle (x-0.5, y-0.1, x+0.5, y+0.1) ;
	      graphBoxEnd() ;
	      remarkRegister(map, seg2->key, seg2->x);
	      if (seg2->key == control->from)
		control->fromBox = box2;
	    }
	  graphBoxEnd () ;
	  if (seg->key == control->from)
	    control->fromBox = ibox;
	}
    }

  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;

  return;
} /* contigDraw */


static BOOL contigCreate(COLINSTANCE instance, OBJ init)
{ 
  instance->draw = contigDraw;
  instance->setSelectBox = contigSetSelect;
  instance->unSelectBox = contigUnselect;
  instance->doColour = contigDoColour;
  instance->followBox = contigFollowBox;
  instance->private = *(instance->proto->convertResults);

  return TRUE;
} /* contigCreate */


struct ProtoStruct gMapContigColumn = {
  0,
  contigCreate,
  0,
  "Contigs",
  0,
  FALSE,
  gMapPhysConvert,
  0
};

/***************************************************************************/
/*               Reversed physical                                         */
/***************************************************************************/

static void revPhysDraw(COLINSTANCE instance, float *offset)
{
  MAPCONTROL map =  instance->map;
  COLCONTROL control = map->control;
  float y = 0, oldy = 0 ;
  int   i, x, oldx ;
  KEY   contig = 0, new ;
  GMAPPHYSEG   *seg ;
  OBJ	obj = 0 ;
  Array segs = instance->private;
  int firstSeg, lastSeg, countDir;

  if (map->mag > 0)
    { firstSeg = 0;
      lastSeg = arrayMax(segs);
      countDir = 1;
    }
  else
    { firstSeg = arrayMax(segs)-1;
      lastSeg = -1;
      countDir = -1;
    }

  graphColor (RED) ;

  for (i = firstSeg ; i != lastSeg ; i+= countDir)
    { seg = arrp(segs,i,GMAPPHYSEG) ;
      if (class(seg->key) != _VClone ||
	  !(obj = bsCreate (seg->key)))
	continue ;
      y = MAP2GRAPH(map,seg->x) ;
      if (bsFindKey (obj, _pMap, contig) &&
	  bsGetData (obj, _bsRight, _Int, &x))
	{ if (x < oldx && y > oldy)
	    { graphFillRectangle (*offset, oldy, *offset+0.5, y) ;
	      remarkRegister(map, contig, seg->x);
	    }
	  oldx = x ;
	  oldy = y ;
	}
      else if (bsGetKey (obj, _pMap, &new) &&
	       bsGetData (obj, _bsRight, _Int, &oldx))
	{ oldy = y ;
	  contig = new ;
	}
      bsDestroy (obj) ;
      if (oldy < control->topMargin + 1)
	oldy = control->topMargin + 1 ;
      else if (oldy  > control->graphHeight - 1)
	break ;
    }

  graphColor (BLACK) ;

  *offset += 0.5 ;

  return;
} /* revPhysDraw */


static BOOL revPhysCreate(COLINSTANCE instance, OBJ init)
{ 
  instance->draw = revPhysDraw;
  instance->private = *(instance->proto->convertResults);

  return TRUE;
} /* revPhysCreate */


struct ProtoStruct gMapRevPhysColumn = {
  0,
  revPhysCreate,
  0,
  "Reversed_physical",
  0,
  FALSE,
  gMapPhysConvert,
  0
};

/*****************************************************/

struct Clone 
{ KEY key;
  float x;
};

Array gMapGetClones(KEY key)
{ 
  Array loci = 0, clones = 0;
  KEY chrom, clone, locus;
  OBJ Locus = 0;
  float x;
  int i;
  struct Clone *seg;

  loci = gMapMapObjects(key);
  clones = arrayCreate(100, struct Clone);

  for (i = 0; i <arrayMax(loci); i++)
    { locus = arr(loci, i, KEYPAIR).obj;
      chrom = arr(loci, i, KEYPAIR).map;
      if (!gMapGetMapObject(locus, key, chrom, 20, &x, 0, &Locus, 0))
	continue;
      
      if (bsGetKey (Locus, _Positive_clone, &clone))
	{ 
	  seg = arrayp(clones, arrayMax(clones), struct Clone) ;
	  seg->key = clone ;
	  seg->x = x ;
	}
      
      bsDestroy(Locus);
    }
  arrayDestroy(loci);

  return clones;
} /* gMapGetClones */


static BOOL pMapToGMapDo (KEY contig, KEY from, int x, BOOL show, KEY *mapp, float *posp)
{
  KEY map, map1 ;
  OBJ obj, Map = 0 ;
  struct Clone *seg ;
  float y1, y2, y, dy ;
  int i, p1 = -1000000, p2 = 1000000, x1, x2 ;
  Array clones;

/* First open the contig and find where  it's mapped in both maps */

  if (!(obj = bsCreate(contig)))
    { if (show) messout ("Can't open contig object") ;
      return FALSE ;
    }
  if (!bsGetKey (obj, _Map, &map) || !map)
    { bsDestroy (obj) ;
      if (show) messout ("This contig is not assigned to a chromosome") ;
      return FALSE ;
    }
  map1 = map ;  /* look to give preference is possible to a genetic_map */
  while (map1)
    {
      if ((Map = bsCreate(map1)))
	{
	  if (bsFindTag (Map, str2tag("Genetic_map")))
	    { map = map1 ; break ; }
	  bsDestroy (Map) ;
	}
      map1 = 0 ;
      bsGetKey (obj, _bsDown, &map1) ;
    }
  bsDestroy (Map) ;

  bsGetData (obj, _pMap, _Int, &p1) ;
  bsGetData (obj, _bsRight, _Int, &p2) ;
  bsDestroy (obj) ;
 
		
  if (from)			/* look for from directly */
    { if (gMapGetMapObject(from, map, 0, 20, 0, 0, 0, 0))
	{ if (show) 
	    { display (map, from, 0) ; return TRUE ; }
	  else
	    return FALSE;
	}
    }

 /* strategy: find 2 gMapped clones around or by x and interpolate */

 if (FLAG_ANY_INTERVAL & gMapGetMapObject(contig, map, 0, 20, &y, &dy, 0, 0))
   { y1 = y - dy; /* these form the bounds for searching for clones */
     y2 = y + dy;

     clones = gMapGetClones(map);
     for (i = 0 ; i < arrayMax(clones) ; ++i)
       { seg = arrp(clones, i, struct Clone) ;
	 if (seg->x > y1 && seg->x < y2 &&
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
     arrayDestroy (clones) ;

     if (y1 > -999999 && y2 < 999999)
       { y = y1 + (y2 - y1) * (x - p1) / (p2 - p1) ;
	 from = KEYMAKE(_VCalcul, 1000.0*(y + 1000.0)) ;
	 if (show)  display (map, from, 0) ;
	 if (mapp) *mapp = map ;
	 if (posp) *posp = y ;
	 return TRUE ;
       }
   }
  if (show) 
    messout ("Sorry, I could not find two flanking markers") ;

  return FALSE ;
} /* pMapToGMapDo */


void pMapToGMap (KEY contig, KEY from, int x)
{
  pMapToGMapDo (contig, from, x, TRUE, 0, 0) ;
}

BOOL  gMapPhysClone2map (KEY clone, KEY *seqp, KEY *mapp, float *xp)
{
  OBJ Clone = 0 ; 
  int x1, x2 ;
  KEY seq, contig, map ;
  float x ;
  BOOL ok = FALSE ;
  
  if ((Clone = bsCreate (clone)) &&
      bsGetKey (Clone, _Sequence, &seq) &&
      bsGetKey (Clone, _pMap, &contig) &&
      bsGetData (Clone, _bsRight, _Int, &x1) &&
      bsGetData (Clone, _bsRight, _Int, &x2) &&
      pMapToGMapDo (contig, 0, x1 + x2 / 2, FALSE, &map, &x))
    {
      x1 = 100 * x ;
      if (xp) *xp = (float)x1/100.0 ;
      if (mapp) *mapp = map ; 
      if (seqp) *seqp = seq ;
      ok = TRUE ;
    }
  bsDestroy (Clone) ;

  return ok ;
} /* gMapPhysClone2map */


void gMapPhysNameClones (void)
{
  KEY map, clone = 0, seq ;
  int level ;
  float x ;
  FILE *f = filopen ("clone_full_name","ace","w") ;

  if (!f)
    return ;
  level = freeOutSetFile (f) ;

  while (lexNext (_VClone, &clone))
    {
      if (gMapPhysClone2map (clone, &seq, &map, &x))
	freeOutf ("Sequence %s\nFull_name %s(%s:%g)\n\n", name(clone), name(seq), name(map), x) ;
    }
  freeOutClose (level) ;

  return;
} /* gMapPhysNameClones */

/************ end of file *************/





 
