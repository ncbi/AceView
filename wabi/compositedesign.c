#include "ac.h"
#include "cdna.h"
#include "makemrna.h"

/**********************************************************************************/
/**********************************************************************************/
/* This strategy is design to handle deep sequencing
 * We construct an oriented graph where each exon is a vertex
 * Each intronis an arc weighted by the number and samples of the supporting tags
 * The introns are ordered by weight
 * For each intron in decreasing order extend it into a path as follows
 *   If the next intron is not yet in a path, extend to the following exon connected by the heaviest intron
 *   If it is attach to the lowest path arriving at this exon from the other side
 *   unless the way i arrive on the exon is much smaller than the way the best path exits the exon
 *   or enough full paths have been constructed
 * Stop when all introns are used or enough paths have been constructed 
 */
typedef struct compositeDesignStruct { Array exons, introns, covering ; BitSet bbIn, bbOut ; KEYSET reads ; int c1, c2 ; AC_HANDLE h ; } DS ;
typedef struct dsVertexStruct { int a1, a2, donor, acceptor, cover, path, type, score, nn ; } DSX ;

/**********************************************************************************/

static void showDs (SC *sc, Array aa)
{
  DSX *up ;
  int ii ;
  
  if (aa)
    {
      for (ii = 0, up = arrp (aa, 0, DSX) ; ii < arrayMax (aa) ; ii++, up++)
	{
	  fprintf( stderr, "%d:: %d %d\t %d %d\t d=%d a=%d c=%d s=%d path=%d\t"
		   , ii, up->a1, up->a2, sc->a1 + up->a1 - 1, sc->a1 + up->a2 - 1
		   , up->donor, up->acceptor, up->cover, up->score, up->path
		   ) ;
	  
	  if (gGene & up->type) fprintf (stderr, "Gene ") ;
	  if (gGap & up->type) fprintf (stderr, "Gap ") ;
	  if (gLink & up->type) fprintf (stderr, "Link ") ;
	  if (gGhost & up->type) fprintf (stderr, "Ghost ") ;
	  if (gFuseToGhost & up->type) fprintf (stderr, "FuseToGhost ") ;
	  if (gDroppedGene & up->type) fprintf (stderr, "DroppedGene ") ;
	  if (gSuspect & up->type) fprintf (stderr, "Suspect ") ;
	  if (gMicro & up->type) fprintf (stderr, "Micro-") ;
	  if (gX & up->type) fprintf (stderr, "Exon ") ;
	  if (gI & up->type) fprintf (stderr, "Intron ") ;
	  if (gJ & up->type) fprintf (stderr, "Fuzzy ") ;
	  if (gS & up->type) fprintf (stderr, "SL ") ;
	  if (gS0 & up->type) fprintf (stderr, "SL0 ") ;
	  if (gReal5p & up->type) fprintf (stderr, "r5p ") ;
	  if (gD & up->type) fprintf (stderr, "Debut ") ;
	  if (gDF & up->type) fprintf (stderr, "DebutFuzzy ") ;
	  if (gCompleteCDS & up->type) fprintf (stderr, "CDS ") ;
	  if (gF & up->type) fprintf (stderr, "Fin ") ;
	  if (gFF & up->type) fprintf (stderr, "FinFuzzy ") ;
	  if (gReal3p & up->type) fprintf (stderr, "r3p ") ;
	  if (gA & up->type) fprintf (stderr, "polyA ") ;
	  if (gB & up->type) fprintf (stderr, "Alter ") ;
	  if (gStolen & up->type) fprintf (stderr, "Stolen ") ;
	  if (gPredicted & up->type) fprintf (stderr, "Predicted ") ;
	  
	  fprintf (stderr, "\n") ;
	}
      showDs (0, 0) ;
    }
  return ;
} /* showDs */

/**********************************************************************************/

static int dsScoreOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->score - b->score ; if (n) return -n ; /* large scores first */
  n = a->type - b->type ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  
  return 0 ;
} /* dsScoreOrder */

/**********************************************************************************/

static int dsA1Order (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  n = a->type - b->type ; if (n) return n ;
  
  return 0 ;
} /* dsA1rder */

/**********************************************************************************/

static void mrnaDesignGetElements (S2M *s2m, SC *sc, DS *ds, Array smrnas)
{
  int ii = 0, jje = 0, jji = 0, j, j2, j2Max, a1, *ip ;
  SMRNA *smrna ;
  Array introns, exons ;
  DSX *dsx, *dsy ;
  HIT *up ;
  Array ks ;
  AC_HANDLE h = ac_new_handle () ;
  BOOL debug = FALSE ;

  exons = ds->exons = arrayHandleCreate (100, DSX, ds->h) ;
  introns = ds->introns = arrayHandleCreate (100, DSX, ds->h) ;
  ds->reads = keySetHandleCreate (ds->h) ;

  if (1)
    {
      /* extract all exon fragments and all introns */
      if (smrnas && arrayMax(smrnas))
	for (ii = 0 ; ii < arrayMax(smrnas) ; ii++)
	  {
	    smrna = arrp (smrnas, ii, SMRNA) ;
	    if (! smrna->clones)
	      continue ;
	    for (j = 0, up = arrp(smrna->hits, 0, HIT) ; j < arrayMax(smrna->hits) ; up++, j++)
	      {
		if ((gX | gI) & up->type) /* exon */
		  {
		    if (gX & up->type) /* exon */
		      dsx = arrayp (exons, jje++, DSX) ;
		    if (gI & up->type) /* exon */
		      dsx = arrayp (introns, jji++, DSX) ;
		    dsx->a1 = up->a1 + smrna->a1 - 1 ;
		    dsx->a2 = up->a2 + smrna->a1 - 1 ;
		    dsx->type = up->type & (gX | gI) ;
		  }
	      }
	  }  
    }
  if (0)  /* stupid idea, i resuse all est from all genes in the current gene */
    {
      BOOL isDown =  (sc->a1 < sc->a2) ? TRUE : FALSE ;
      /* extract all exon fragments and all introns */
      for (j = 0, up = arrp(s2m->plainHits, 0, HIT) ; j < arrayMax(s2m->plainHits) ; up++, j++)
	{
	  if ((gX | gI) & up->type) /* exon */
	    {
	      if (gX & up->type) /* exon */
		dsx = arrayp (exons, jje++, DSX) ;
	      if (gI & up->type) /* exon */
		dsx = arrayp (introns, jji++, DSX) ;
	      if (isDown)
		{
		  dsx->a1 = up->a1 - sc->a1 + 1 ;
		  dsx->a2 = up->a2 - sc->a1 + 1 ;
		}
	      else
		{
		  dsx->a1 = sc->a1 - up->a2 + 1 ;
		  dsx->a2 = sc->a1 - up->a1 + 1 ;
		}
	      dsx->type = up->type & (gX | gI) ;
	    }
	}
    }  

   if (debug)
     {
       fprintf (stderr, "... getElements Z\n") ;
       showDs (sc, exons) ;
     }
   arraySort (exons, dsA1Order) ;
   arrayCompress (exons) ;
   if (debug)
     {
       fprintf (stderr, "... getElements A\n") ;
       showDs (sc, exons) ;
     }
   /* now we split the exons in subparts */ 
   ks = arrayHandleCreate (2 * arrayMax (exons) + 1, int, h) ;
   for (j = j2 = 0, dsx = arrayp (exons, j, DSX) ; j < arrayMax (exons) ; j++, dsx++)
     {
      array (ks, j2++, int) = dsx->a1 ;
      array (ks, j2++, int) = dsx->a2 + 1 ;
     }
   arraySort (ks, intOrder) ;
   arrayCompress (ks) ;
   j2Max = keySetMax (ks) ;
   ds->exons = arrayHandleCreate (100, DSX, ds->h) ;
   for (jje = j = j2 = 0, dsx = arrayp (exons, j, DSX), ip = arrp (ks, 0, int) ; j < arrayMax (exons) ; j++, dsx++)
     {
       a1 = dsx->a1 ;
       while (j2 > 0 && *ip > a1) { j2-- ; ip-- ; }
       while (*ip < dsx->a1) { j2++; ip++ ;}
       while (++j2 < j2Max && *(++ip) <= dsx->a2 + 1)
	 {
	   dsy = arrayp (ds->exons, jje++, DSX) ;
	   dsy->a1 = a1 ; dsy->a2 = *ip - 1 ; dsy->type = gX ; a1 = *ip ;
	 }
     }

   arraySort (ds->exons, dsA1Order) ;
   arrayCompress (ds->exons) ;
   if (debug) 
     {
       fprintf (stderr, "... getElements B\n") ;
       showDs (sc, ds->exons) ;
     }

   arraySort (introns, dsA1Order) ;
   arrayCompress (introns) ;
   if (debug)
     showDs (sc, introns) ; 

   ac_free (h) ;
   return ;
} /* mrnaDesignGetElements */

/**********************************************************************************/
/* get the coounts in introns included in XI XY XW */	
static void mrnaDesignGetIntronSupport (KEY intron, DSX *vp)
{
  int ir ;
  OBJ Intron = bsCreate (intron) ;
  static Array units = 0 ;
  BSunit *uu ;

  units = arrayReCreate (units, 12, BSunit) ;
  if (Intron)
    {
      if (bsFindTag (Intron, _Other) || bsFindTag (Intron, _ct_ac))
	vp->type |= gSuspect ;
      else
	vp->type &= ~gSuspect ;
      if (bsGetArray (Intron, str2tag ("Group_U"), units, 6))
	{
	  for (ir = 0 ; ir < arrayMax (units) ; ir += 6)
	    {
	      uu = arrp (units, ir, BSunit) ;
	      if (uu[2].f > vp->cover) /* seqs */
		vp->cover = uu[2].f ;
	      if (uu[4].f > vp->cover) /* tags */
		vp->cover = uu[4].f ;
	    }
	}
      if (bsGetArray (Intron, str2tag ("Validated_U"), units, 2))
	{
	  for (ir = 0 ; ir < arrayMax (units) ; ir += 2)
	    {
	      uu = arrp (units, ir, BSunit) ;
	      if (uu[1].i > vp->cover)
		vp->cover = uu[1].i ;
	    }
	}
      if (bsGetArray (Intron, str2tag ("RNA_seq"), units, 1))
	{
	  for (ir = 0 ; ir < arrayMax (units) ; ir += 1)
	    {
	      uu = arrp (units, ir, BSunit) ;
	      if (uu[0].i > vp->cover)
		vp->cover = uu[0].i ;
	    }
	}
      if (bsGetArray (Intron, str2tag ("Magic_any_any"), units, 1))
	{
	  for (ir = 0 ; ir < arrayMax (units) ; ir += 1)
	    {
	      uu = arrp (units, ir, BSunit) ;
	      if (uu[0].i > vp->cover)
		vp->cover = uu[0].i ;
	    }
	}
      bsDestroy (Intron) ;
    }

  return ;
} /* mrnaDesignGetIntronSupport */
		    
/**********************************************************************************/
/* get the XI or XD Composite supporting the introns and their counts */
static void  mrnaDesignGetSupport  (DS *ds, S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  HIT *up ;
  DSX *vp, *ssp, *ssp2 ;
  int ii, jj, ir, iss, iss2, a1, a2, b1, b2, ssMax, u1 ;
  BOOL isDown ;
  BOOL debug = FALSE ;
  BOOL hasData = FALSE ;
  Array units = 0 ;
  Array ss = 0, ss2 = 0 ;
  BSunit *uu ;
  OBJ Est = 0 ;
  KEYSET xks, xws = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  ss = arrayHandleCreate (256, DSX, h) ; iss = 0 ;
  xks = keySetHandleCreate (h) ;
  xws = keySetHandleCreate (h) ;

  /* grab the introns in the plainHits EST alignments */
  isDown = (sc->a1 < sc->a2) ? TRUE : FALSE ;
  units = arrayHandleCreate (12, BSunit, h) ;
  for (ii = 0, up = arrp (s2m->plainHits, 0, HIT) ; ii < arrayMax (s2m->plainHits) ; ii++, up++)
    {
      if (! keyFindTag (up->est, _Composite))
	continue ;
      a1 = isDown ? up->a1 - sc->a1 + 1 : sc->a1 - up->a2 + 1 ;
      a2 = isDown ? up->a2 - sc->a1 + 1 : sc->a1 - up->a1 + 1 ;
      /* identify the introns in the intron table */
      if (up->type & gI)
	for (jj = 0, vp = arrp (ds->introns, 0, DSX) ; jj < arrayMax (ds->introns) ; jj++, vp++)
	  {
	    if (vp->a1 != a1 || vp->a2 != a2)
	      continue ;
	    if ((Est = bsCreate (up->est)))
	      {
		if (bsGetArray (Est, _Composite, units, 5))
		  {
		    keySet (ds->reads, keySetMax (ds->reads)) = up->est ;
		    if (! strncmp (name(up->est), "XD_", 3))
		      {
			if (arrayMax (units) >= 5)
			  {
			    uu = arrp (units, 0, BSunit) ;
			    if (vp->donor < uu[1].i) vp->donor = uu[1].i ; 
			    if (vp->acceptor < uu[3].i) vp->acceptor = uu[3].i ; 
			  }
		      }
		    else
		      {
			if (arrayMax (units) >= 1)
			  {
			    int k ;
			    uu = arrp (units, 0, BSunit) ;
			    k =  uu[0].i ; 
			    if (bsFindTag (Est, _Other))
			      k /= 5 ;
			    if (vp->cover < k)
			      vp->cover  = k ;
			  }
		      }
		  }
		if (bsGetArray (Est, _Intron, units, 1) &&
		    arrayMax (units) == 1
		    )
		  {
		    KEY intron = 0 ;
		    if (bsGetKey (Est, _Intron, &intron))
		      mrnaDesignGetIntronSupport (intron, vp) ;
		  }
		bsDestroy (Est) ;
	      }
	    if ((up->type & gSuspect) && ! (vp->type & gJ))
	      { 
		vp->type |= gJ ; /* divide only once */
		vp->cover /= 5 ;
		vp->donor /= 5 ;
		vp->acceptor /= 5 ;
		vp->score /= 5 ;
	      }
	  }
      /* register the support of the exon support */
      if (up->type & gX)
	{
	  if ( ! strncmp (name(up->est), "XC_", 3) || ! strncmp (name(up->est), "XF_", 3))  /* coverons */
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  hasData = TRUE ;
		  if (bsGetArray (Est, _Composite, units, 1) && arrayMax (units) >= 1)
		    for (ir = 0 ; ir < arrayMax (units) ; ir ++)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 + 10 ; ssp->a2 = a2 - 10 ; ssp->cover = uu[0].i ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }
		  bsDestroy (Est) ;
		}	  
	      continue ;
	    }      
	  
	  if ( ! strncmp (name(up->est), "XG_", 3) || ! strncmp (name(up->est), "XH_", 3) || ! strncmp (name(up->est), "XE_", 3))
	    {
	      int iss0 = iss ;
	      if ((Est = bsCreate (up->est)))
		{
		  hasData = TRUE ;
		  if (bsGetArray (Est, _Composite, units, 3) && arrayMax (units) >= 3)
		    for (ir = 0 ; ir < arrayMax (units) ; ir += 3)
		      {
			uu = arrp (units, ir, BSunit) ;
			if (ir) 
			  { 
			    /* oldb1 = b1 ;  */
			    /*  oldc = uu[2].i > uu[-1].i ? uu[-1].i : uu[2].i ; */
			  } 
			b1 = a1 + uu[0].i - 1 ; 	
			b2 = a1 + uu[1].i - 1 ; 
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = b1 ; ssp->a2 = b2 ; ssp->cover = uu[2].i ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }
		  bsDestroy (Est) ;
		}
	      /* shrink the higher expressed segments */
	      if ( ! strncmp (name(up->est), "XG_", 3) || ! strncmp (name(up->est), "XH_", 3))
		{
		  int ii ;
		  for (ii = iss0 ; ii < iss - 1 ; ii++)
		    {
		      DSX *ssq ;
		      ssp = arrayp (ss, ii, DSX) ;
		      ssq = ssp + 1 ;
		      if (ssp->a1 < ssp->a2 && ssp->a2 == ssq->a1)
			{
			  if (ssp->cover < ssq->cover)
			    { ssp->a2 += 10 ; ssq->a1 += 10 ; }
			  if (ssp->cover > ssq->cover)
			    { ssp->a2 -= 10 ; ssq->a1 -= 10 ; }
			}
		      if (ssp->a1 > ssp->a2 && ssp->a2 == ssq->a1)
			{
			  if (ssp->cover < ssq->cover)
			    { ssp->a2 -= 10 ; ssq->a1 -= 10 ; }
			  if (ssp->cover > ssq->cover)
			    { ssp->a2 += 10 ; ssq->a1 += 10 ; }
			}
		    }
		  ssp = arrayp (ss, iss0, DSX) ;
		  ssp->a1 += (ssp->a1 < ssp->a2 ? +10 : -10) ;
		  ssp = arrayp (ss, iss - 1, DSX) ;
		  ssp->a2 += (ssp->a1 < ssp->a2 ? -10 : +10) ;
		  if (debug)
		    fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		}
	      continue ;
	    }	  
      
	  if (! strncmp (name(up->est), "XI_", 3) ||  ! strncmp (name(up->est), "XK_", 3))
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  BOOL j ;
		  hasData = TRUE ;
		  if (bsGetArray (Est, _Composite, units, 1) && arrayMax (units) >= 1)
		    for (ir = 0 ; ir < arrayMax (units) ; ir += 1)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ;
			/* only count as supported the region close to the intron */
			j = keySetInsert (xks, up->est) ^ up->reverse ;
			if (j) /* first exonic apparition */
			  ssp->a1 = a2 - 30 ;
			else
			  ssp->a2 = a1 + 30 ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }	    
		  bsDestroy (Est) ;
		}
	      continue ;
	    }	  
      
	  if (! strncmp (name(up->est), "XW_", 3) || ! strncmp (name(up->est), "XY_", 3))
	    {  /* double introns incoded in the name */
	      if (keySetInsert (xks, up->est) && /* first exonic apparition */
		  (Est = bsCreate (up->est)))
		{
		  int cover = 0 ;
		  const char *ccp = strstr (name(up->est), "__") ;
		  if (! bsFindTag (Est, _Other) &&
		      ccp && bsGetData (Est, _Composite, _Int, &cover) && cover > 0)
		    {
		      int k, u1 = 0, u2 = 0, v1 = 0, v2 = 0 ;
		      char c = 0 ;
		      hasData = TRUE ;
		      k = sscanf (ccp + 2, "%d_%d_%d_%d%c", &u1, &u2, &v1, &v2, &c) ;
		      if (k == 4) /* found 6 numbers, zero terminated */
			{
			  int isDown1 = isDown ? 1 : -1 ;
			  int isDown2 = (u1 < u2 ? 1 : -1) ;
			  int da1 = sc->a1 ;

			  if (! ds->c1)
			    {
			      OBJ Cosmid = bsCreate (s2m->cosmid) ;
			      int c1 = 0 ;
			      KEY map ;

			      if (Cosmid)
				{
				  if (bsGetKey (Cosmid, _IntMap, &map) &&
				      bsGetData (Cosmid, _bsRight, _Int, &c1)
				      )
				    ds->c1 = c1 ;
				  bsDestroy (Cosmid) ;
				}
			    }
			  da1 = sc->a1 + ds->c1 - 1 ;

			  ssp = arrayp (ss, iss++, DSX) ;  /* first exon */
			  ssp->a1 = isDown1 * (u1 - 30 * isDown2 - da1) + 1  ; 
			  if (0) ssp->a1 = isDown1 * (up->a1 - sc->a1) + 1 ;
			  ssp->a2 = isDown1 * (u1 - isDown2   - da1 ) + 1 ;
			  ssp->cover = cover ;
			  if (debug) 
			      fprintf(stderr, "A::%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;

			  ssp = arrayp (ss, iss++, DSX) ;  /* second exon */
			  ssp->a1 =  isDown1 * (u2 + isDown2  - da1) + 1  ; 
			  ssp->a2 =  isDown1 * (v1 -isDown2  - da1) + 1  ; 
			  ssp->cover = cover ;
			  if (debug)
			      fprintf(stderr, "B::%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;

			  ssp = arrayp (ss, iss++, DSX) ;  /* third exon */
			  ssp->a1 =  isDown1 * (v2 + isDown2  - da1) + 1   ; 
			  ssp->a2 =  isDown1 * (v2 + 30 * isDown2  - da1) + 1  ; 
			  if (0) ssp->a2 = isDown1 * (up->a2 - sc->a1) + 1 ;
			  ssp->cover = cover ;
			  if (debug)
			      fprintf(stderr, "C::%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
			}
		    }	    
		  bsDestroy (Est) ;
		}
	      continue ;
	    }	  

	  if (! strncmp (name(up->est), "XD_", 3))
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  BOOL j ;
		  hasData = TRUE ;
		  if (bsGetArray (Est, _Composite, units, 5) && arrayMax (units) >= 5)
		    for (ir = 0 ; ir < arrayMax (units) ; ir += 5)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			j = keySetInsert (xks, up->est) ^ up->reverse ;
			if (j) /* first exonic apparition */
			  { ssp->a1 = a2 - 30 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ; }
			else
			  {  
			    j = keySetInsert (xws, up->est) ;
			    if (j)  /* second exonic apparition = supported central exon */
			      { ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[4].i ; }
			    else
			      { ssp->a1 = a1 ; ssp->a2 = a1 + 30 ; ssp->cover = uu[4].i ; }
			  }
			if (debug)
			  fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }	    
		  bsDestroy (Est) ; 
		}
	      continue ;
	    }	  
	  
	  if (! strncmp (name(up->est), "XA_", 3)) 
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 2) && arrayMax (units) >= 2)
		    for (ir = 0 ; ir < arrayMax (units) ; ir += 2)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = (uu[1].i > uu[0].i ? uu[1].i : uu[0].i) ;
			ssp->type |= (gA | gX) ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d\n", name(up->est), a1, a2, uu[0].i) ;
		      }	    
		  bsDestroy (Est) ;
		}
	      continue ;
	    }
	  if (! strncmp (name(up->est), "Xcds_", 5)) 
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 1) && arrayMax (units) >= 1)
		    for (ir = 0 ; ir < arrayMax (units) ; ir++)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ;
			ssp->type |= (gCompleteCDS | gX) ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d\n", name(up->est), a1, a2, uu[0].i) ;
		      }	    
		  bsDestroy (Est) ;
		}
	      continue ;
	    }
	  if (! strncmp (name(up->est), "Xends_", 6))
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 1) && arrayMax (units) >= 1)
		    for (ir = 0 ; ir < arrayMax (units) ; ir ++)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ;
			 if (! strncmp (name(up->est), "Xends_ELF_", 9) && isDown) 
			   ssp->type |= (gS | gX) ;
			 if (! strncmp (name(up->est), "Xends_ERF_", 9) && isDown) 
			   ssp->type |= (gReal3p | gX) ;
			 if (! strncmp (name(up->est), "Xends_ERR_", 9) && ! isDown) 
			   ssp->type |= (gS | gX) ;
			 if (! strncmp (name(up->est), "Xends_ELR_", 9) && ! isDown) 
			   ssp->type |= (gReal3p | gX) ;
			if (debug)
			  fprintf(stderr, "....%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }
		  bsDestroy (Est) ;
		}
	      continue ;
	    }
	  
	  if (! strncmp (name(up->est), "XSL", 3))
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 2) && arrayMax (units) >= 2)
		    for (ir = 0 ; ir < arrayMax (units) ; ir += 2)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = (uu[1].i > uu[0].i ? uu[1].i : uu[0].i) ;
			ssp->type |= (gS | gX) ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d:%d \n", name(up->est), a1, a2, uu[0].i, uu[4].i) ;
		      }	    
		  bsDestroy (Est) ;
		}
	      continue ;
	    }
	  
	  if (1)
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 1) && arrayMax (units) >= 1)
		    for (ir = 0 ; ir < arrayMax (units) ; ir ++)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }
		  bsDestroy (Est) ;
		}	  
	      continue ;
	    }
	}
    }
  if (! hasData)
    arrayMax (ss) = 0 ;
  if (! arrayMax (ss))
    goto done ;
  /* at each position choose the maximal value */
  arraySort (ss, dsA1Order) ;

  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetSupport phase 1 done\n") ;
      showDs (sc, ss) ;
    }
  if (1)
    {
      float intronBonus = 3.0 ;
      for (jj = 0, vp = arrp (ds->introns, 0, DSX) ; jj < arrayMax (ds->introns) ; jj++, vp++)
	if (! (vp->type & gSuspect))
	  vp->cover *= intronBonus ;
	else
	  vp->cover /= 2 ;
    }
   
  /*  register all boundaries */
  ds->covering = ss2 = arrayHandleCreate (2 * arrayMax (ss), DSX, ds->h) ;
  for (jj = iss = 0, ssp = arrp (ss, 0, DSX) ; iss < arrayMax (ss) ; iss++, ssp++)
    {
      if(ssp->a1 >= ssp->a2) { ssp->cover = 0 ; continue ; }
      ssp2 = arrayp (ss2, jj++, DSX) ;
      ssp2->a1 = ssp2->a2 = ssp->a1 ; ssp2->type = ssp->type & gDebut ; 
      ssp2 = arrayp (ss2, jj++, DSX) ;
      ssp2->a1 = ssp2->a2 = ssp->a2 ; ssp2->type = ssp->type & gFin ; 
      ssp2 = arrayp (ss2, jj++, DSX) ; /* this point is needed to go back to zero at the end of the intervals */
      ssp2->a1 = ssp2->a2 = ssp->a2 + 1 ;
      ssp2 = arrayp (ss2, jj++, DSX) ; /* this point is needed to go back to zero at the end of the intervals */
      ssp2->a1 = ssp2->a2 = ssp->a1 - 1 ;
    }
  arraySort (ss2, dsA1Order) ; 
  for (iss = 0, ssp = arrp (ss, 0, DSX) ; iss < arrayMax (ss) ; iss++, ssp++)
    { /* so that compress will eliminate coords duplicates created with different flags */
      for (iss2 = iss+1, ssp2 = ssp + 1 ;  iss2 < arrayMax (ss) && ssp2->a1 == ssp->a1 && ssp2->a2 == ssp->a2; iss2++, ssp2++)
	{ ssp2->type |= ssp->type ; ssp->type = ssp2->type ; }
    }
  arrayCompress (ss2) ;
  /* register the max values in each zone */ 
  for (jj = iss = iss2 = 0, ssp = arrp (ss, 0, DSX), ssp2 = arrp (ss2, 0, DSX) ; iss < arrayMax (ss) ; iss++, ssp++)
    {
      /* position ssp2 on ssp->a1 */
      while (iss2 > 0 && ssp2->a1 > ssp->a1) { iss2-- ; ssp2-- ; }
      while (iss2 < arrayMax (ss2) && ssp2->a1 < ssp->a1) { iss2++ ; ssp2++ ; }
      /* possibly increment the coverage */
      while (iss2 < arrayMax (ss2) && ssp2->a1 <= ssp->a2) 
	{ 
	  if (ssp2->cover < ssp->cover)
	    ssp2->cover = ssp->cover ;
	  iss2++ ; ssp2++ ; 
	}
    }
  if (debug)
    {
      fprintf (stderr, "max value in each zone\n") ;
      showDs (sc, ss2) ;
    }

  /* look for weak introns donor/acceptor */
   ii = 0 ; ssp = arrp (ss2, 0, DSX) ; ssMax = arrayMax (ss2) ;
   if (1)
     for (jj = 0, vp = arrp (ds->introns, 0, DSX) ; ii < ssMax && jj < arrayMax (ds->introns) ; jj++, vp++)
       {
	 DSX* vp2 ;
	 int jj2 ;
	 long int cover = vp->cover ;

         cover *= 100 ;

	 a1 = vp->a1 ;
	 a2 = vp->a2 ;
	 
	 /* compare the cover to the donor exon (a1 - 1) and the retained intron (a1) */
	 while (ii > 0 && ssp->a1 > a1-5) { ii-- ; ssp-- ; }
	 for ( ; ssp->a1 <= a1+5 && ii < ssMax ; ssp++, ii++)
	   if (ssp->a1 >= a1 - 5 && ssp->a1 <= a1 + 5)
	     if (ssp->cover > cover)
	       { vp->donor = -1 ; break ; }
	 /* compare to other introns starting at same site */
	 for (jj2 = jj - 1, vp2 = vp - 1 ; jj2 >= 0  && vp2->a1 == a1 ; jj2--, vp2--)
	    if (vp2->cover > cover)
	       { vp->donor = -1 ; break ; }
	 for (jj2 = jj + 1, vp2 = vp + 1 ; jj2 <  arrayMax (ds->introns)  && vp2->a1 == a1 ; jj2++, vp2++)
	    if (vp2->cover > cover)
	       { vp->donor = -1 ; break ; }
	 /* compare the cover to the donor exon (a2 + 1) and the retained intron (a2) */
	 while (ii > 0 && ssp->a1 > a2 - 5) { ii-- ; ssp-- ; }
	 for ( ; ssp->a1 <= a2 + 5 && ii < ssMax ; ssp++, ii++)
	   if (ssp->a1 >= a2 - 5  && ssp->a1 <= a2 + 5)
	     if (ssp->cover > cover)
	       { vp->acceptor = -1 ; break ; }
	 for (jj2 = jj + 1, vp2 = vp + 1 ; jj2 <  arrayMax (ds->introns)  && vp2->a1 < a2 ; jj2++, vp2++)
	    if (vp2->a2 == a2 && vp2->cover > cover)
	       { vp->acceptor = -1 ; break ; }
	 if (vp->donor == -1 && vp->acceptor == -1)
	   vp->cover = vp->score = 0 ; 
	 if (ii) { ii-- ; ssp-- ; }
       }
  

  /* compute the average support of each exon fragment */
  ii = 0 ; ssp = arrp (ss2, 0, DSX) ; ssMax = arrayMax (ss2) ;
  for (jj = 0, vp = arrp (ds->exons, 0, DSX) ; jj < arrayMax (ds->exons) ; jj++, vp++)
    {
      long int nn, nnn ;
      int type = 0 ;
      nnn = nn = 0 ;
      a1 = vp->a1 ;
      a2 = vp->a2 ;
      while (ii > 0 && ssp->a1 > a1) { ii-- ; ssp-- ; }
      u1 = a1 ; if (u1 < ssp->a1) u1 = ssp->a1 ; nn = ssp->cover ; 
      type = ssp->type ; 
      for ( ; ii < ssMax && ssp->a1 < a1 ; ii++, ssp++) 
	{
	  nn = ssp->cover ;  type = ssp->type ; 
	}
      if (0) vp->type |= (type & gDebut) ; type = 0 ;
      for ( ; ii < ssMax && ssp->a1 < a2 ; ii++, ssp++)
	{
	  nnn += nn * (ssp->a1 - u1) ;
	  if (ssp->a1 >= u1) u1 = ssp->a1 ; 
	  nn = ssp->cover ;
	}
      if (0 && ii < ssMax && a2 == ssp->a1)
	vp->type |= (ssp->type & gFin) ; 
      if (a2 >= u1) 
	nnn += nn * (a2 - u1 + 1) ;

      vp->cover = nnn / (a2 - a1 + 1) ;
      if (nnn >= 0 && vp->cover < 1) vp->cover = 1 ;
    }

  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetSupport phase 2 done\n") ;
      showDs (sc, ss) ;
    }

  if (1) /* Xcds is used to mask Xends */
    for (jj = 0, vp = arrp (ss, 0, DSX) ; jj < arrayMax (ss) ; jj++, vp++)
      if (vp->type & gCompleteCDS)
	for (iss = 0, ssp = arrp (ss2, 0, DSX) ; iss < arrayMax (ss2) ; iss++, ssp++)
	  {
	    if (
		ssp->a1 > vp->a1 &&
		ssp->a1 <= vp->a2
		)
	      {
		ssp->type |=  gCompleteCDS ;
		ssp->type &= ~gReal3p ;  /* Xcds is used to mask Xends */
		ssp->type &= ~gA ;  /* Xcds is used to mask Xends */
	      }
	  }

  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetSupport phase 3 done\n") ;
      showDs (sc, ss2) ;
    }

 if (1) /* reestablish the polyA flags */
    for (jj = 0, vp = arrp (ds->exons, 0, DSX) ; jj < arrayMax (ds->exons) ; jj++, vp++)
      {
	for (iss = 0, ssp = arrp (ss2, 0, DSX) ; iss < arrayMax (ss2) ; iss++, ssp++)
	  {
	    if ((ssp->type & gCompleteCDS) &&
		ssp->a2 < vp->a2 &&
		ssp->a2 >= vp->a1
		)
	      vp->type |= gCompleteCDS ;
	    if ((ssp->type & gS) &&
		ssp->a2 < vp->a2 &&
		ssp->a2 >= vp->a1
		)
	      vp->type |= gS ;
	    if ((ssp->type & (gA | gReal3p)) &&
		ssp->a1 > vp->a1 &&
		 ssp->a1 <= vp->a2
		)
	      {
		vp->type |= (ssp->type & (gA | gCompleteCDS | gReal3p)) ;
		if (vp->type & gCompleteCDS &  gReal3p) vp->type &= ~gReal3p ;
		if (vp->type & gCompleteCDS &  gA) vp->type &= ~gA ;
	      }
	  }
     }
   if (debug)
     {
       fprintf (stderr, "mrnaDesignGetSupport  done\n") ;
       showDs (sc, ds->exons) ;
     }
 done:
   ac_free (h) ;
  return ;
} /* mrnaDesignGetSupport */

/**********************************************************************************/
/* copy the exon up0 into ds->exons, possibly splitting it */
static void mrnaDesignSplitOneExon (DS *ds, SC *sc, DSX *up0)
{
  DSX *up, *xp ;

  
  if ( up0->type & gX &&      
      (up0->type & (gReal5p | gS | gS0)) &&
      (up0->type & (gReal3p | gA)))
    {
      int ii ;
      int iMax = arrayMax (ds->covering) ;
      int a1 = up0->a1 ;
      int a2 = up0->a2 ;
      int cover = 0 ;
      unsigned int type = 0 ;
      
      for (ii = 0, xp = arrp (ds->covering, 0, DSX) ; ii < iMax ; xp++, ii++)
	{
	  if (xp->a1 >= a1 && xp->a1 <= a2)
	    {
	      if  (ii && 
		   ( 5 * xp->cover < (xp-1)->cover || 
		     xp->cover > 5 * (xp-1)->cover)
		   )
		{
		  if (xp->a1 <= a1) continue ;
		  up = arrayp (ds->exons, arrayMax (ds->exons), DSX) ;
		  up->a1 = a1 ;
		  up->a2 = xp->a1 - 1 ;
		  a1 = xp->a1 ;
		  up->type = type | gX ;
		  up->cover = cover/(up->a2 - up->a1 + 1) ;
		  cover = 0 ;
		  type = 0 ;
		}
	      type |= xp->type ;
	      cover += xp->cover * ((xp+1)->a1 - xp->a1) ;
	    }
	}
      up = arrayp (ds->exons, arrayMax (ds->exons), DSX) ;
      up->a1 = a1 ;
      up->a2 = xp->a1 - 1 ;
      if (up->a2 + 1 > up->a1)
	up->cover = cover/(up->a2 - up->a1 + 1) ;
    }
  else
    {
      up = arrayp (ds->exons, arrayMax (ds->exons), DSX) ;
      *up = *up0 ;
    }
  
  return ;
} /* mrnaDesignSplitOneExon */

/**********************************************************************************/

static void mrnaDesignSplitExons (DS *ds, SC *sc)
{
  BOOL debug = FALSE ;

  if (debug)
    {
      fprintf (stderr, "mrnaDesignSpliExon starts\n") ;
      showDs (sc, ds->covering) ;
      fprintf (stderr, "   .. original exons\n") ;
      showDs (sc, ds->exons) ;
    }

  if (arrayMax (ds->exons))
    {
      int ii ;
      DSX *up ;
      Array aa = arrayCopy (ds->exons) ;

      arrayMax (ds->exons) = 0 ;
      for (ii = 0, up = arrp (aa, 0, DSX) ; ii < arrayMax (aa) ; ii++, up++)
	mrnaDesignSplitOneExon (ds, sc, up) ; 
      arrayDestroy (aa) ;
    }

  if (debug)
    {
      fprintf (stderr, "   .. final exons\n") ;
      showDs (sc, ds->exons) ;
    }

  return ;
} /* mrnaDesignSplitExons  */

/**********************************************************************************/
/* extend the path Up via the best supported route */
static void mrnaDesignExtendDown (DS *ds, Array segs, KEYSET ks, int path, int nn0, int *nIntronp)
{ 
  int a2, ii, nn, score ; 
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;

  up = arrp (segs, nn0, DSX) ;
  if (up->type & (gA | gReal3p)) 
    return ;
  a2 = up->a2 + 1 ;
  for (score = nn = 0, vp = up + 1, ii = nn0 + 1 ; ii < sMax && vp->a1 <= a2 ; ii++, vp++)
    {
      if (vp->a1 != a2) continue ;
      if (score < vp->score) { score = vp->score ; nn = vp->nn ; wp = vp ; }
    }
  if (score)
    {
      /* kill alternative branches  with too low relative score */
      if (0)
	for (vp = up + 1, ii = nn0 + 1 ; ii < sMax && vp->a1 <= a2 ; ii++, vp++)
	  {
	    if (vp->a1 != a2) continue ;
	    if (score > 200 * vp->score) { vp->score = 0 ; }
	  }

      keySet (ks, arrayMax(ks)) = nn ;
      if (! wp->path) wp->path = path ;
      if (up->type & gI)
	(*nIntronp)++ ;

      if (0) fprintf (stderr, "Path=%d  %d:%d >> %d:%d\n", path, up->a1, up->a2, wp->a1, wp->a2) ;
      mrnaDesignExtendDown (ds, segs, ks, path, nn, nIntronp) ;
    }
  return ;
} /* mrnaDesignExtendDown */

/**********************************************************************************/
/* extend the path Up via the best supported route */
static void mrnaDesignExtendUp (DS *ds, Array segs, KEYSET ks, int path, int nn0, int *nIntronp, int *nStartp)
{ 
  int a1, ii, nn, score ; 
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;

  up = arrp (segs, nn0, DSX) ;
  if (*nStartp && (up->type & gI))
    return ;
  if (*nIntronp && (up->type & (gA | gReal3p)))
    return ;
  a1 = up->a1 - 1 ;
  for (score = nn = ii = 0, vp = arrp (segs, 0, DSX) ; ii < sMax && vp->a1 <= a1 ; ii++, vp++)
    {
      if (vp->a2 != a1) continue ;
      if (score < vp->score) { score = vp->score ; nn = vp->nn ; wp = vp ; }
    }
  if (score)
    {
      /* kill alternative branches  with too low relative score */
      if (0)
	for (ii = 0, vp = arrp (segs, 0, DSX) ; ii < sMax && vp->a1 <= a1 ; ii++, vp++)
	  {  
	    if (vp->a2 != a1) continue ;
	    if (score > 200 * vp->score) { vp->score = 0 ; }
	  }

      keySet (ks, arrayMax(ks)) = nn ;
      if (! wp->path) wp->path = path ;
      if (up->type & gI)
	(*nIntronp)++ ;
      if (up->type & gS)
	(*nStartp)++ ;

      if (0) fprintf (stderr, "Path=%d  %d:%d << %d:%d\n", path, up->a1, up->a2, wp->a1, wp->a2) ;
      mrnaDesignExtendUp (ds, segs, ks, path, nn, nIntronp, nStartp) ;
    }
  return ;
} /* mrnaDesignExtendUp */

/**********************************************************************************/

static int mrnaDesignExport (S2M *s2m, DS *ds, Array segs, KEYSET ks, int path, Array smrnas)
{
  int a1 = 0, ii,jj,  nn ;
  DSX *up ;
  SMRNA *smrna ;
  HIT *vp = 0 ;
  BOOL debug = FALSE ;

  smrna = arrayp(smrnas, path - 1, SMRNA) ;
  arrayDestroy (smrna->dnas) ;
  arrayMax (smrnas) = path ;
  if (! smrna->hits)
    smrna->hits = arrayHandleCreate (keySetMax(ks), HIT, s2m->h) ;
  else
    smrna->hits = arrayReCreate (smrna->hits, keySetMax(ks), HIT) ;
  if (! smrna->clones)
    smrna->clones = keySetHandleCreate (s2m->h) ; 
  keySetSort (ks) ; /* now the elements are in a1 order and correctly chained */
  for (ii = jj = 0 ; ii < keySetMax (ks) ; ii++)
    {
      nn = keySet (ks, ii) ;
      up = arrp (segs, nn, DSX) ;
      if (!ii && (up->type & gI))
	continue ;
      if (!vp) a1 = smrna->a1 = up->a1 ; smrna->a2 = up->a2 ;
      if (!vp || up->type != vp->type)
	{	
	  vp = arrayp (smrna->hits, jj++, HIT) ;
	  vp->a1 = up->a1 - a1 + 1 ;
	}
      vp->a2 = up->a2 - a1 + 1 ;
      vp->type = up->type ;
      if (debug) fprintf(stderr, "Path %d :: %d :: %d :: %s %d %d\t%d\td%d:a%d:c%d\t%d%s\ttvp %d:%d\n"
	      , path, ii, nn, up->type & gX ? "exon" : "intron", up->a1, up->a2, up->score, up->donor, up->acceptor, up->cover
	      , up->path, path == up->path ? "*****" : ""
	      , vp->a1, vp->a2
	      ) ;
    }
  return path ;
} /* mrnaDesignExport */

/**********************************************************************************/

static int mrnaDesignFindPaths (S2M *s2m, SC *sc, DS *ds, Array smrnas)
{
  int ii, jj, ii2, maxScore, path, a2Max, nIntron, nStart ;
  int eeMax = arrayMax (ds->exons) ;
  int iiMax = arrayMax (ds->introns) ;
  Array segs, segs2 ;
  DSX *up, *vp, *vp2 ; ;
  KEYSET ks = keySetHandleCreate (ds->h) ;

  /* group exons/introns */
  segs = arrayHandleCreate (eeMax + iiMax +1, DSX, ds->h) ;
  for (ii = jj = 0, vp = arrp (ds->exons, 0, DSX) ; jj < eeMax ; ii++, jj++, vp++)
    {	    
      up =  arrayp (segs,ii, DSX) ;
      *up = *vp ; up->score = vp->cover ; up->nn = ii ;
    }
  for (jj = 0, vp = arrp (ds->introns, 0, DSX) ; jj < iiMax ; ii++, jj++, vp++)
    {	    
      up =  arrayp (segs,ii, DSX) ;
      *up = *vp ; up->score = vp->cover ; up->nn = ii ;
      if (up->score < vp->donor) up->score = vp->donor ;
      if (up->score < vp->acceptor) up->score = vp->acceptor ;
      up->score <<= 1 ;
      if (0 && ! up->score) up->score = 1 ;
    }
  arraySort (segs, dsA1Order) ;
  eeMax = arrayMax (segs) ;
  for (ii = 0, vp = arrp (segs, 0, DSX) ; ii < eeMax ; ii++, vp++)
    vp->nn = ii ;

  if (0)
    showDs (sc, segs) ;

  /* kill the non extremal exon extensions which are below 15 bp and do not hook onto an intron */
  for (ii = a2Max = 0, vp = arrp (segs, 0, DSX) ; ii < eeMax ; ii++, vp++)
    if (a2Max < vp->a2) a2Max = vp->a2 ;
 
 for (ii = 0, vp = arrp (segs, 0, DSX) ; ii < eeMax ; ii++, vp++)
    {
      if ((vp->type & gX) && vp->a2 < vp->a1 + 15)
	{
	  int ok1 = 0, ok2 = 0 ;
	  int a1 = vp->a2 + 1, a2 = vp->a1 - 1 ;
	  if (vp->a1 == 1)
	    ok1++ ;
	  else
	    for (jj = ii - 1, vp2 = vp - 1 ; ok1 == 0 && jj >= 0 ; jj--, vp2--)
	      if (vp2->a2 == a2)
		{
		  if (1 || (vp2->type & gI) ||
		      (vp->score > .8 * vp2->score && vp->score < 1.2 * vp2->score))
		    ok1 = 1 ;
		}
	  if (vp->a2 == a2Max) 
	    ok2 = 1 ;
	  else
	    for (jj = ii + 1, vp2 = vp + 1 ; ok2 == 0 && jj < eeMax && vp2->a1 <= a1 ; jj++, vp2++)
	      {
		if (vp2->a1 == a1)
		  {
		    if (1 || (vp2->type & gI) ||   /* vp2 intron sticcking to vp exon */
			(vp->score > .8 * vp2->score && vp->score < 1.2 * vp2->score)) /* vp2 exon with similar score directly sticking to vp exon */
		      ok2 = 1 ;
		  }
	      }
	  if (ok1 + ok2 < 2) 
	    vp->score = -1 ;
	}
    }


  /* sort a copy by score */
  segs2 = arrayHandleCopy (segs, ds->h) ;
  arraySort (segs2, dsScoreOrder) ;

  for (ii2 = 0, path = 0, vp2 = arrp (segs2, 0, DSX), maxScore = ((vp2->type & gI) ? vp2->score/6 :  vp2->score) ; ii2 < eeMax && (100 * vp2->score > maxScore || vp2->score >= 1) ; ii2++, vp2++)
    {
      /* find highest scoring seg not yet incorporated in a path */
      vp = arrp (segs, vp2->nn, DSX) ;
      if (vp->path || ! vp->score)
	continue ;
      /* do not start on a retained intron */
      if (vp->type & gX)
	{
	  DSX *wp, *wp2 ;
	  int jj2 ;
	  BOOL ok2 = TRUE ;

	  for (jj2 = 0, wp2 = arrp (segs2, 0, DSX) ; ok2 && jj2 < ii2 ; jj2++, wp2++)
	    {
	      wp =  arrp (segs, wp2->nn, DSX) ;
	      if ((wp->type & gI) && wp->a1 <= vp->a1 && wp->a2 >= vp->a2 && wp->score > 3 * vp->score && 4*vp->score < maxScore)
		ok2 = FALSE ;  /* wp is a true  intron better and including  my stupid  vp retained intron (i.e. artefactual exon) */
	    }
	  if (! ok2)
	    continue ;
	}

      path++ ;
      vp->path = path ;
      ks = keySetReCreate (ks) ;
      keySet (ks, 0) = vp2->nn ;
      nIntron = (vp2->type & gI) ? 1 : 0 ;
      nStart = (vp2->type & gS) ? 1 : 0 ;
      mrnaDesignExtendDown (ds, segs, ks, path, vp2->nn, &nIntron) ;
      mrnaDesignExtendUp (ds, segs, ks, path, vp2->nn, &nIntron, &nStart) ;
      mrnaDesignExport (s2m, ds, segs, ks, path, smrnas) ;
    }
  return path ;
} /* mrnaDesignFindPaths */

/**********************************************************************************/
/* flag the reads contributing introns that were actually rejected as void or poor */
static void mrnaDesignCleanUp (S2M *s2m, SC *sc, DS *ds, Array smrnas)
{
  int ii ;
  KEYSET ks ;
  vTXT txt ;
  const char **errors = 0 ;
  AC_DB db = ac_open_db (0, errors) ; /* local database, cannot fail */

  if (! keySetMax (ds->reads)) 
    return ;
  ks = query (ds->reads, "IS XY_* && (ct_ac || other || small_deletion) && (ct_ac > 1 || other > 1 || (ct_ac && other) || ! composite > 50)") ;
  if (keySetMax (ks))
    {
      txt = vtxtHandleCreate (ds->h) ;
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
	vtxtPrintf (txt, "Sequence %s\n-D Is_read\n\n", name(keySet(ks, ii))) ;
      ac_parse (db, vtxtPtr (txt), errors, 0, ds->h) ;
    } 
  keySetDestroy (ks) ;
  
  ks = query (ds->reads, "Intron && ! from_gene  && ! composite > 50 ; >cdna_clone") ;
  if (keySetMax (ks))
    {
      txt = vtxtHandleCreate (ds->h) ;
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
	vtxtPrintf (txt, "cDNA_clone %s\nIgnore_this_clone\n\n", name(keySet(ks, ii))) ;
      ac_parse (db, vtxtPtr (txt), errors, 0, ds->h) ;
    } 
  keySetDestroy (ks) ;
  ac_free (db) ;
  return ;
} /* mrnaDesignCleanUp */

/**********************************************************************************/
/**********************************************************************************/

BOOL mrnaDesignUsingCompositeStrategy (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  AC_HANDLE h = 0 ;
  DS ds ;
  int i ;
  HIT *up ;
  BOOL ok = FALSE ;
  BOOL debug = FALSE ;

  if (0 || arrayMax (smrnas) < 1) 
    return FALSE; 

  /* check if we are using composite reads */
  for (i = 0, up = arrp (s2m->geneHits, 0, HIT), ok = FALSE ; !ok && i < arrayMax (s2m->geneHits) ; i++, up++)
    if (keyFindTag (up->est, _Composite))
      ok = TRUE ;
  if (! ok)
    return FALSE ;

  memset (&ds, 0, sizeof (ds)) ;
  ds.h = h = ac_new_handle () ;

  mrnaDesignGetElements (s2m, sc, &ds, smrnas) ;
  mrnaDesignGetSupport (&ds, s2m, sc, gmrna, smrnas) ;
  mrnaDesignSplitExons (&ds, sc) ;
  if (debug)
    {
      showDs (sc, ds.exons) ;
      showDs (sc, ds.introns) ; 
    }
  smrnas = arrayReCreate (smrnas, 12, SMRNA) ;
  mrnaDesignFindPaths (s2m, sc, &ds, smrnas) ;
  mrnaDesignCleanUp (s2m, sc, &ds, smrnas) ;
  
  ac_free (h) ;
  return ok ;
} /*  mrnaDesignUsingCompositeStrategy */

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
