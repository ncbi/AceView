/*  File: pmapconvert.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: convert stuff from pmap - trying to split it
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 17 17:09 1998 (fw)
 * * Feb 10 19:39 1994 (rd): whole file cleaned up by Richard 940210
 * Created: Tue Nov 30 23:39:24 1993 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: pmapconvert.c,v 1.2 2014/11/30 03:20:46 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "key.h"
#include "pmap.h"
#include "pmap_.h"
#include "lex.h"
#include "bs.h"
#include "a.h"
#include "sysclass.h"
#include "systags.h"
#include "classes.h"
#include "tags.h"
#include "session.h"
#include "grid.h"		/* for GRIDMAP */
#include "client.h"


#define KeyToSeg(__segs, __key, __kseg) \
{ \
  int __i, __nSegs=arrayMax(__segs); \
  (__kseg)=NULL; \
  for (__i=1 /*seg[0] is the class*/; __i<__nSegs; __i++) \
  { \
    SEG *__seg=arrayp((__segs), __i, SEG); \
    if (__seg!=NULL && 0<iskey(__seg->key) && __seg->key==(__key)) \
    { \
      (__kseg)=__seg; break; \
    } \
  } \
}
  /*returns NULL if seg not found*/

void pMapMakeAll (void) 	/* to be called during update, and from alignMap button */
{
  KEY contig = 0, pmap ;	/* 0 primes lexNext() */
  Array segs ;

  if (!isWriteAccess())
    { messout("Sorry, you do not have Write Access.");
      return ;
    }

  while (lexNext (_VContig,&contig))
    { lexaddkey (name (contig), &pmap, _VpMap) ;
      if ((segs = pMapConvert (0, contig, pmap)))
	arrayDestroy (segs) ;
    }
}

/***************************************************************/

int pMapCompareSeg (const void *a, const void *b) /*for call to qsort handled by arraySort*/

/*  1/ seg[0] gets sorted, so any change to the sort primitive must preserve its pre-eminence
    2/ buried clones must come after their parents
*/
{
  const SEG *s = (const SEG *)a, *t = (const SEG *)b ;

  if (!s || !t) 
    messcrash ("pMapCompareSeg found a NULL seg") ;

  if (IsSet(s, FLAG_IS_BURIED_CLONE) &&
      !IsSet(s, FLAG_BURIED_IS_POSITIONED) && 
      s->parent == t->parent)
    return -1 ;
  else if (IsSet(t, FLAG_IS_BURIED_CLONE) && 
	   !IsSet(t, FLAG_BURIED_IS_POSITIONED) && 
	   t->parent == s->parent)
    return 1 ;
  else if (SegMidPt(s) <= SegMidPt(t))
    return -1 ;
  else
    return 1 ;
}

/***********************************/

static Associator probeHash = 0 ;
static BOOL hashProbes = FALSE ;

static void pMapAddSeg(Array segs, int *pj, KEY clone, OBJ Clone, int ctgL)
{
  extern int sequenceLength (KEY sequence) ;
  static KEY poly1 = 0 ;
  static Array flatA = 0 ;
  char *cp ;
  int	i, j = *pj ;
  KEY   locus, sequence, probe ;
  OBJ   Sequence ;
  int   ctgR = ctgL ;
  float midX;
  SEG 	*seg ;
  Array probeArray ;
  GRIDMAP *map ;

			/* initialise statics */
  if (!poly1)
    { int _VGrid ;
      KEY key ;
      if (lexword2key ("Grid", &key, _VMainClasses) &&
	  (_VGrid = KEYKEY(key)) &&
	  lexword2key("POLY1", &poly1, _VGrid)) ; /* poly1 found */
      else
	poly1 = 1 ;		/* nonsense KEY */
    }
  flatA = arrayReCreate (flatA, 8, BSunit) ;

  bsGetData(Clone, _bsRight, _Int, &ctgR) ; /* ctgL otherwise */

  seg = arrayp(segs,j++,SEG) ;
  seg->key = clone ;
  seg->parent = clone ;
  seg->x0 = ctgL ; 
  seg->x1 = ctgR ; 
  midX = 0.5 * (ctgL + ctgR) ;

  seg->flag = 0 ;
  
  if (bsFindTag (Clone,_Bands))	/* FingerPrint no good because of Flag */
    seg->flag |= FLAG_FINGERPRINT ;
  if (bsFindTag (Clone,_Canonical_for))
    seg->flag |= FLAG_CANONICAL ;
  if (bsFindTag (Clone,_Cosmid_grid) || bsFindTag (Clone,_Canon_for_cosmid))
    seg->flag |= FLAG_COSMID_GRID ;
  if (poly1 != 1 && bsFindKey (Clone,_Gridded, poly1))
    seg->flag |= FLAG_YAC_GRID ;

  if (bsFindTag (Clone,str2tag("Cosmid")))
    seg->flag |= FLAG_IS_COSMID ;
  else if (bsFindTag (Clone,str2tag("YAC")))
    seg->flag |= FLAG_IS_YAC ;
  else if (bsFindTag (Clone,str2tag("Fosmid")))
    seg->flag |= FLAG_IS_FOSMID ;
  else if (bsFindTag (Clone,str2tag("cDNA")))
    seg->flag |= FLAG_IS_CDNA ;

  if (bsGetKey (Clone,_Sequence, &sequence) &&
     !externalServer &&                         /* too slow on xclient, suppressed */
      (i = sequenceLength (sequence)))
    { 
      float dx = i / 1200.0 ; /* nominal (inflated) kb per map unit */
				/* should be dependent on a ?Contig field */
      seg->flag |= FLAG_SEQUENCED ;
      seg = arrayp(segs,j++,SEG) ;
      seg->key = sequence ;
      seg->parent = clone ;
      if (dx < 1)
	dx = 1 ;		/* so you can see it */
      seg->x0 = midX - 0.5 * dx ;
      seg->x1 = midX + 0.5 * dx ;

      if ((Sequence = bsCreate (sequence)))
	{ BOOL isTransposon = FALSE ;
	  if (bsFindTag (Sequence, _Allele) && 
	      bsFlatten (Sequence, 4, flatA))
	    for (i = 0 ; i < arrayMax (flatA) ; i += 4)
	      { 
		cp = arr(flatA, i+3, BSunit).s ;
		
		if (cp) while (*cp)
		  if (!strchr ("AGCTagct-", *cp++))
		    { isTransposon = TRUE ; break ;}
	 
		
		if (isTransposon) /* a transposon */
		  { float x = 0.5 * (arr(flatA, i+1, BSunit).i + 
				     arr(flatA, i+2, BSunit).i) / 1200.0 ;
		  seg = arrayp(segs,j++,SEG) ;
		  seg->key = arr(flatA, i, BSunit).k ;
		  seg->parent = clone ;
		  seg->x0 = midX + 0.5 * (x - dx) ;
		  seg->x1 = FMINUSINF ;
		  }
	      }
	  bsDestroy (Sequence) ;
	}
    }

  if (bsFindTag(Clone, _Remark)
      && bsFlatten(Clone,2,flatA))
    for (i = 0 ; i<arrayMax(flatA) ; i+=2 )
      { seg = arrayp(segs,j++,SEG) ;
	seg->key = arr(flatA, i+1, BSunit).k ;
	seg->parent = clone ;
        seg->x0=midX; seg->x1=FMINUSINF;
	seg->flag = 0 ;
#ifdef ACEDB5
        seg->flag = FLAG_REMARK ;
#endif
	if (arr(flatA, i, BSunit).k  != _General_remark)
	  seg->flag |= FLAG_WORK_REMARK ;
      } 

  if (hashProbes && bsGetKey (Clone,_Positive_locus,&locus)) do
    { if (!(assFind(probeHash, assVoid(locus), &probeArray) &&
	    arrayExists(probeArray)))
	{ probeArray = arrayCreate (8, GRIDMAP) ;
	  assInsert (probeHash,  assVoid(locus), probeArray) ; 
	}
      map = arrayp(probeArray, arrayMax(probeArray), GRIDMAP) ;
      map->x1=ctgL;
      map->x2=ctgR;
      map->ctg = 1 ;		/* 0 value is special */
    } while (bsGetKey (Clone,_bsDown,&locus)) ;


  if (hashProbes && bsGetKey (Clone,_Positive_probe,&probe)) do
    { if (!(assFind(probeHash, assVoid(probe), &probeArray) &&
	    arrayExists(probeArray)))
	{ probeArray = arrayCreate (8, GRIDMAP) ;
	  assInsert (probeHash, assVoid(probe), probeArray) ; 
	}
      map = arrayp(probeArray, arrayMax(probeArray), GRIDMAP) ;
      map->x1=ctgL;
      map->x2=ctgR;
      map->ctg = 1 ;		/* 0 value is special */
    } while (bsGetKey (Clone,_bsDown,&probe)) ;
  
  *pj = j ;
}

/********************************/

Array pMapConvert (PhysMap look, KEY contig, KEY pmap) 
     /* extract display information for contig into segs array */
{
  static Array mapArray ;
  OBJ	Clone, Contig ;
  KEY 	clone, jprobe ;
  int 	i, j, m, ctgL;
  SEG   *seg;
  Array segs = 0, probeArray, flagBuriedClonesOf = 0 ;
  void *v;
  
   if (externalServer) 
     { 
       return 0 ; /* much too slow */
       externalServer (contig, 0, 0, TRUE) ;/* get object and neighbours */
     }
  if (!(Contig = bsCreate (contig)))
    { messout ("Can't find raw data for contig %s.",
	       name(contig)) ;
      return 0 ;
    }
  if (!lexlock (pmap))
    { messout ("Structure for %s map is locked by another user.",
	       name (pmap)) ;
      return 0 ;
    }


  if (look && look->segs && IsntSet(look, FLAG_SHOW_BURIED_CLONES))
    { flagBuriedClonesOf = arrayCreate(10, KEY) ; /* save per clone DISPLAY_BURIED information */
      for (i = 1, m = 0 ; i < arrayMax(look->segs) ; i++)
	if (IsSet((seg = arrp(look->segs, i, SEG)), FLAG_DISPLAY_BURIED_CLONES))
	  array(flagBuriedClonesOf, m++, KEY) = seg->key;
    }

  if (look && !(segs = look->segs))
    segs = arrayGet(pmap, SEG, segFormat) ;
  segs = arrayReCreate(segs, 128, SEG) ;

  probeHash = assReCreate (probeHash) ;
  mapArray = arrayReCreate (mapArray, 8, GRIDMAP) ;
  
  j = 1 ;			/* segs are boxes, so start at 1 */
  hashProbes = TRUE ;
  if (bsGetKey(Contig, _Clone, &clone)) do
    if ((Clone=bsCreate(clone)))
      { if (!CloneIsBuriedClone(Clone) && /* else deal with it in pMapLiftBuriedClones */
	    bsFindKey (Clone, _pMap, contig) && 
	    bsGetData(Clone, _bsRight, _Int, &ctgL))
	  pMapAddSeg(segs, &j, clone, Clone, ctgL);
        bsDestroy(Clone);
      } 
    while (bsGetKey(Contig, _bsDown, &clone)) ;
  hashProbes = FALSE ;

				/* now cluster probes, loci */
  v = 0 ; probeArray = 0 ;
  while (assNext (probeHash, &v, &probeArray))
    { clone = assInt(v) ;
      if ((Clone = bsCreate(clone)))
	{ gridCluster (probeArray, mapArray, 200) ;
	  for (i = 0 ; i < arrayMax(mapArray) ; ++i)
	    { float x1 = arrp(mapArray, i, GRIDMAP)->x1 ;
              float x2 = arrp(mapArray, i, GRIDMAP)->x2 ;
              int xx = (x1 + x2) / 2 ;
	      if (class(clone) == _VClone)
		{ jprobe = j ;
		  pMapAddSeg(segs, &j, clone, Clone, xx) ;/* WAS: (int)0.5*(x1+x2)) == 0 ! */
		  seg = arrp(segs, jprobe, SEG) ;
		}
	      else
		{ seg = arrayp(segs,j++,SEG) ;
		  seg->key = clone ;
		  seg->parent = clone ;
		  seg->x0 = xx ;
		  seg->flag = 0 ;
		}
	      seg->flag |= FLAG_PROBE ;
              seg->x1 = FMINUSINF ;
	      /*        if (seg->x0 == 0.0) seg->x0 = ERROR ; */
              if (x2 < x1)
                seg->flag |= FLAG_ERROR ;
	    }
	  bsDestroy(Clone) ;
	}
      arrayDestroy (probeArray) ;
    }

  /*  recompute buried clones, if any are being displayed, to preserve consistency 
      over recalculations of the contig due to editing operations. */

  if (look)
    {
      if (IsSet(look, FLAG_SHOW_BURIED_CLONES))
	{
	  int nSegs = arrayMax(segs) ; /* since will be adding to segs during loop */
	  for (i = 1 ; i < nSegs ; ++i)
	    if ((seg = arrp(segs, i, SEG)))
	      pMapLiftBuriedClones(segs, seg, FALSE) ;
	  Set(look, FLAG_BURIED_CLONES_ATTACHED);
	}
      else if (flagBuriedClonesOf)
	{
	  for (i = 0 ; i < arrayMax(flagBuriedClonesOf) ; ++i)
	    {
	      KeyToSeg(segs, arr(flagBuriedClonesOf, i, KEY), seg);
	      if (seg)
		{ 
		  Set(seg, FLAG_DISPLAY_BURIED_CLONES);
		  pMapLiftBuriedClones(segs, seg, TRUE);
		}
	    }
	  arrayDestroy(flagBuriedClonesOf);
	}
    }

  array(segs, 0, SEG).x0 = array(segs, 0, SEG).x1 = FMINUSINF ; /* so it stays first */
   arraySort(segs, pMapCompareSeg) ;

  bsDestroy (Contig) ;

  if ((Contig = bsUpdate(contig)))
    { int x = arrp(segs, 1, SEG)->x0 ;
      bsAddData (Contig, _pMap, _Int, &x) ;
      i = arrayMax(segs)-1 ;
      do x = arrp(segs, i, SEG)->x1 ;
      while (x == FMINUSINF && --i > 0) ;
      if (i)
	bsAddData (Contig, _bsRight, _Int, &x) ;
      bsSave (Contig) ;
    }

  arrayStore (pmap, segs, segFormat) ;
  lexunlock (pmap) ;

  return segs ;
}

/***********************************************************************/
/* Neil note: if RecalculateContig is called with
   fromScratch==FALSE, any affected segs should be updated manually
   beforehand (see updateSegs in pmapcons.c) 
*/

void pMapRecalculateContig (PhysMap look, BOOL fromScratch)
{
  KEY  contig;

   if (!lexlock(look->key))
     { messout ("Sorry the physical map %s is locked elsewhere",
		name(look->key));
     return ;
     }
   lexunlock(look->key) ;

   if (fromScratch)
     {
       if (lexword2key (name(look->key), &contig, _VContig))
	 look->segs = pMapConvert (look, contig, look->key) ;
     }
   else
     {
       if (lexlock(look->key)) /* already tested, cannot fail */
	 { arraySort(look->segs, pMapCompareSeg) ;
	 arrayStore(look->key, look->segs, segFormat) ;
	 lexunlock(look->key) ;
	 }
     }

   look->scrollBox = 0 ;
   look->activebox = 0 ;
}

/************************* end of file **************************/
 
