/*  File: align2.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **  Realigns the various maps in dna unit                    **
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  8 11:32 1999 (fw)
 * * May 16 22:35 1993 (cgc): rd started again from scratch
 * * Dec 11 12:18 1991 (mieg): myText for non interactive messages
 * Created: Fri Nov 29 14:31:10 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: align.c,v 1.4 2015/08/18 23:24:27 mieg Exp $ */

#include "acedb.h"		/* must be outside for mac */

#include "display.h"
#include "sysclass.h"
#include "systags.h"
#include "tags.h"
#include "classes.h"
#include "session.h"
#include "a.h"
#include "bs.h"
#include "dna.h"
#include "lex.h"
#include "plot.h"
#include "query.h"
#include "gmap.h"
#include "vmap.h"

static Graph alignGraph = 0 ;
static int line ;

static void alignDisplay (void) ;

typedef struct
  { KEY key ;
    KEY parent ; /* clone associated with this item */
    float x0, x1 ;
    unsigned int flag ;
  } PMAPSEG ;
#define pMapFormat "kkffi"

/*****************************************************/
/* Routine to get an array of gMap Clone and contigs */
/*****************************************************/

/*********************************************/
/*********************************************/

   /* Graphic or non graphic printouts */

static void textGraphic(char *cp, int x, int y)
{
  line += y ;
  graphText(cp, x, line) ;
}

static void textNonGraphic(char *cp, int x, int y)
{
  if (!y)
    messdump ("  ") ;
  while (y--)
    messdump ("\n") ;
  messdump (cp) ;
}

static void (*myText)(char*,int,int) ;

/********************************************/

typedef struct { float x, y ; } PAIR ;

/* strategy is to write an ace file */

static void alignAllContigs (void)
{
  KEY contig, pmap, map ;
  Associator map2pairs = 0;
  Array psegs, pairs, units  ;
  int i, j, pmin, pmax ;
  PAIR *pair, *pair2 ;
  PMAPSEG *pseg ;
  float f1, f2 ;
  FILE *fil ;
  OBJ  obj ;
  BOOL hasPrinted ;

  alignDisplay () ;

  myText ("Before aligning you need up-to-date pMaps.", 
	  1, 1) ;
  myText ("If necessary, make these first.", 1, 1) ;
  myText ("", 1, 1) ;

  if (!(fil = filqueryopen (0, 0, "ace", "w", "Ace file to write")))
    return ;
  /* process contigs */
  contig = 0 ;
  units = arrayCreate (12, BSunit) ;
  while (lexNext (_VContig, &contig))
    { if (!lexReClass (contig, &pmap, _VpMap) ||
	  !(psegs = arrayGet (pmap, PMAPSEG, pMapFormat)))
	continue ;
      pmin = 1000000 ; 
      pmax = -1000000 ;
      map2pairs = assReCreate (map2pairs) ;
      for (i = 0 ; i < arrayMax(psegs) ; ++i)
	{ pseg = arrp(psegs,i,PMAPSEG) ;
	  if (class (pseg->key) == _VClone)
	    { if (pseg->x0 < pmin)
		pmin = pseg->x0 ;
	      if (pseg->x1 > pmax)
		pmax = pseg->x1 ;
	    }
	  else if (class (pseg->key) == _VLocus &&
		   (obj = bsCreate(pseg->key)))
	    { if (bsFindTag (obj, _Map) && 
		  bsFlatten (obj, 3, units))
		for (j = 0 ; j < arrayMax(units) ; j += 3)
		  if (arr(units, j+1, BSunit).k == _Position)
		    { map = arr(units, j, BSunit).k ;
		      if (!assFind (map2pairs, assVoid(map), &pairs))
			assInsert (map2pairs, assVoid(map),
			   (pairs = arrayCreate(8, PAIR))) ;
		      pair = arrayp(pairs, arrayMax(pairs), PAIR) ;
		      pair->x = pseg->x0 ;
		      pair->y = arr(units, j+2, BSunit).f ;
		    }
	      bsDestroy (obj) ;
	    }
	}
      arrayDestroy (psegs) ;
      hasPrinted = FALSE ;
      for (map = 0 ; assNext (map2pairs, &map, &pairs) ; )
	{
	  if (!hasPrinted)
	    { fprintf (fil, "\nContig %s\n", name(contig)) ;
	      fprintf (fil, "pMap %d %d\n", pmin, pmax) ;
	      hasPrinted = TRUE ;
	    }
	  if (arrayMax(pairs) >= 2)
	    { pair = arrp(pairs, 0, PAIR) ;
	      pair2 = arrp(pairs, arrayMax(pairs)-1, PAIR) ;
	      f1 = pair->y - (pair->x - pmin) * 
		(pair2->y - pair->y) / (pair2->x - pair->x) ;
	      f2 = pair2->y + (pmax - pair2->x) * 
		    (pair2->y - pair->y) / (pair2->x - pair->x) ;
	      fprintf (fil, "Map %s Left %.3g // %d hits\n"
		            "Map %s Right %.3g\n",
		       name(map), f1, arrayMax(pairs), name(map), f2) ;
	    }
	  else
	    { pair = arrp(pairs, 0, PAIR) ;
	      f1 = pair->y - (pair->x - pmin) * 0.004 ;
	      f2 = pair->y + (pmax - pair->x) * 0.004 ;
	      fprintf (fil, "Map %s Position %.3g Error %.3g // 1 hit\n",
		       name(map), pair->y, (pmax - pmin) * 0.004) ;
	    }
	  arrayDestroy (pairs) ;
	}
      if (messIsInterruptCalled ())
	break ;
    }

  arrayDestroy (units) ;
  assDestroy (map2pairs) ;

  filclose (fil) ;
  return ;
}

/*******************************************************************/
/****************** code to make sequence maps *********************/

/* rd 950826: sorry - very worm specific just now */

typedef struct {
  KEY key ;
  KEY map ;
  int leftEnd, rightEnd, leftMap ;
} CINFO ;

int cinfoCompare (const void *A, const void *B) /* for qsort() */
{ 
  CINFO *a = (CINFO*)A, *b = (CINFO*)B ;

  if (a->map > b->map) return 1 ;
  if (a->map < b->map) return -1 ;
  if (a->leftMap > b->leftMap) return 1 ;
  if (a->leftMap < b->leftMap) return -1 ;
  return 0 ;
}

void sMapMake (FILE *fil, BOOL isSeq)		/* make sequence maps */
{
  KEYSET contigs, genCan, genClone ;
  int i, offset = 0, x1, x2, j ;
  float pos ;
  KEY key, map ;
  OBJ obj ;
  CINFO *cinf ;
  Array cInfo ;
  static BSMARK mark = 0 ;
  extern Array pMapConvert (void*, KEY, KEY) ;
#define CONTIG_GAP 50		/* nominal size of gap */

	/* find contigs to put on the maps (from last time) */

  contigs = query (0, "FIND Map Sequence-* ; FOLLOW Contig") ;
  cInfo = arrayCreate (keySetMax(contigs), CINFO) ;
  for (i = 0 ; i < keySetMax(contigs) ; ++i)
    if ((obj = bsCreate (keySet(contigs, i))))
      { if (bsGetKey (obj, _Map, &map)) do
	  { mark = bsMark (obj, mark) ;
	    if (strlen (name(map)) > 9 &&
		!strncmp (name(map),"Sequence-",9) &&
		bsPushObj (obj) &&
		bsGetData (obj, _Left, _Float, &pos))
	      { cinf = arrayp(cInfo, arrayMax(cInfo), CINFO) ;
		cinf->key = keySet(contigs,i) ;
		cinf->map = map ;
		cinf->leftMap = pos + 0.5 ; /* 0.5 for rounding */
		bsGoto (obj, mark) ; 	    /* to get out of subobj */
		if (!bsGetData (obj, _pMap, _Int, &cinf->leftEnd) ||
		    !bsGetData (obj, _bsRight, _Int, &cinf->rightEnd))
		  { Array segs ;
		    lexaddkey (name(cinf->key), &map, _VpMap) ;
		    segs = pMapConvert (0, cinf->key, map) ;
		    arrayDestroy (segs) ;
		    bsDestroy (obj) ;
		    obj = bsCreate (cinf->key) ;
		    if (!bsGetData (obj, _pMap, _Int, &cinf->leftEnd) ||
			!bsGetData (obj, _bsRight, _Int, &cinf->rightEnd))
		      { messout ("Problem with ends of %s", name(cinf->key)) ;
			--arrayMax (cInfo) ;
		      }
		  }
		break ;		/* out of do-while loop */
	      }
	    bsGoto (obj, mark) ;
	  }
        while (bsGetKey (obj, _bsDown, &map)) ;
	bsDestroy (obj) ;
	if (messIsInterruptCalled ())
	  break ;
      }
    else
      messout ("Failed to find object %s", keySet(contigs,i)) ;

	/* sort the contigs into map order, and rewrite the map objects */

  qsort (arrp(cInfo,0,CINFO), arrayMax(cInfo), 
	 			sizeof(CINFO), cinfoCompare) ;
  map = 0 ;
  for (i = 0 ; i < arrayMax(cInfo) ; ++i)
    { cinf = arrp(cInfo, i, CINFO) ;
      if (cinf->map != map)
	{ if (map)
	    { fprintf (fil, "\nMap %s\n", name(map)) ;
	      fprintf (fil, "Extent 0 %d\n", offset) ;
	    }
	  map = cinf->map ;
	  fprintf (fil, "\nMap %s\n", name(map)) ;
	  fprintf (fil, "-D Contig\n") ;
	  fprintf (fil, "-D %s\n", isSeq ? "Sequence" : "Clone") ;
	  offset = 0 ;
	}
      else
	offset += CONTIG_GAP ;
      fprintf (fil, "\nContig %s\n", name(cinf->key)) ;
      cinf->leftMap = offset ;
      fprintf (fil, "Map %s Left %d\n", name(cinf->map), offset) ;
      offset += (cinf->rightEnd - cinf->leftEnd) ;
      fprintf (fil, "Map %s Right %d\n", name(cinf->map), offset) ;
    }
  if (map)
    { fprintf (fil, "\nMap %s\n", name(map)) ;
      fprintf (fil, "Extent 0 %d\n", offset) ;
    }

	/* now put objects on map */

  genCan = query (0, "FIND Genome_Sequence") ;
  genClone = keySetCreate() ;
  for (i = 0 ; i < keySetMax (genCan) ; ++i)
    if (lexword2key (name(keySet(genCan,i)), &key, _VClone))
      keySet(genClone, keySetMax(genClone)) = key ;
  keySetSort (genClone) ;
  keySetDestroy (genCan) ;

  for (i = 0 ; i < keySetMax(genClone) ; ++i)
    if ((obj = bsCreate (keySet (genClone, i))))
      { if (!bsGetKey (obj, _pMap, &map))
	  { 
	    if ( /* just bsGet, if() is for compiler happiness */
		bsGetKey (obj, _Approximate_Match_to, &key) ||
		bsGetKey (obj, _Exact_Match_to, &key) ||
		bsGetKey (obj, _Funny_Match_to, &key))
	      {} ;
	    bsDestroy (obj) ;
	    if (!(obj = bsCreate (key)))
	      continue ;
	    if (!bsGetKey (obj, _pMap, &map))
	      { bsDestroy (obj) ;
		continue ;
	      }
	  }
	if (bsGetData (obj, _bsRight, _Int, &x1) &&
	    bsGetData (obj, _bsRight, _Int, &x2))
	  { for (j = 0 ; j < arrayMax(cInfo) ; ++j)
	      if ((cinf = arrp(cInfo, j, CINFO))->key == map)
		break ;
	    if (j < arrayMax(cInfo))
	      { x1 += cinf->leftMap - cinf->leftEnd ;
		x2 += cinf->leftMap - cinf->leftEnd ;
		fprintf (fil, "\n%s %s\n", 
			 isSeq ? "Sequence" : "Clone",
			 name(keySet(genClone,i))) ;
		fprintf (fil, "Map %s Left %d\n", name(cinf->map), x1) ;
		fprintf (fil, "Map %s Right %d\n", name(cinf->map), x2) ;
	      }
	  }

	bsDestroy (obj) ;
	if (messIsInterruptCalled ())
	  break ;
      }

  arrayDestroy (cInfo) ;
  keySetDestroy (genClone) ;

  fprintf (fil, "\n// end of file\n") ;
}

void sMapMakeAll (void)		/* make sequence maps */
{
  FILE *fil ;
  BOOL isSeq ;

  if (!(fil = filqueryopen (0, 0, "ace", "w", "Ace file for result")))
    return ;

  isSeq = messQuery ("Do you want Sequence objects on the map (Yes),"
		     " or Clone objects (No)?") ;

  sMapMake (fil, isSeq) ;
  filclose (fil) ;
}

/*****************************************/
/************* action routines ***********/
/*****************************************/

static void alignDestroy (void)
{    
  alignGraph = 0 ;
  myText = textNonGraphic ;
}

/**********************************************************/

extern 
  void pMapMakeAll(void), cMapMakeAll(void), gMapMakeAll(void);

static MENUOPT alignMenu[] =
  {
   {graphDestroy, "Quit"},
   {help, "Help"},
   {graphPrint,"Print"},
   {pMapMakeAll,"Make pMaps"},
   {vMapMakeAll,"Make vMaps"},
   {gMapMakeAll,"Make gMaps"},
   {cMapMakeAll,"Make cMaps"},
   {sMapMakeAll,"Make sequence map ace file"},
   {alignAllContigs, "Align contigs - write ace file"},
    {0, 0} 
   } ;

/*************************************************************/

static void alignDisplay (void)
{
  if (!graphActivate (alignGraph))
    messcrash ("alignDisplay lost its graph.") ;
 
  graphPop () ;
  graphClear () ;
  graphButtons (alignMenu,3.,1.,40.) ;

  graphText ("To halt any action hit F4.", 5, 8.85) ;
  
  line = 8 ;
  graphTextBounds (80, line += 1 ) ;
  graphRedraw () ;
}

/************************************************/
/***********  public routines   ****************/

void alignMaps (void)
{
  myText = textGraphic ;

  if (graphActivate (alignGraph))
    { graphPop () ;
      return ;
    }

  if (!(alignGraph = displayCreate(DtAlign)))
    return ;

  graphTextBounds (80,50) ;   /* needed to for text box sizing */
  graphRegister (DESTROY,alignDestroy) ;
  graphMenu (alignMenu) ;  
 
  alignDisplay() ;
}

void alignMapsNonInteractive (void)
{
  void (*oldText)(char*,int,int) = myText ;

  myText = textNonGraphic ;
  messStatus ("Aligning maps") ;
  alignAllContigs () ;
  myText = oldText ;
}

/*************************************************************/
/*************************************************************/
