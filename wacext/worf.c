#include "ac.h"
#include "freeout.h"
#include "utils.h"
#include "vtxt.h"
#include "dict.h"

static BOOL debug = FALSE ;

typedef struct worf_struct { int nGenes
			     , hasGF, hasHmm, hasGFandHMM,  hasGForHMM
			     , hasOnlyGF, hasOnlyHMM
			     , hasOnlyGF_gfNegative, hasOnlyHMM_hmmNegative
			     , hasGForHMM_gfNegative, hasGForHMM_hmmNegative
			     , hasGFandHMM_gfNegative, hasGFandHMM_hmmNegative
			     , hasTG, hasOST, nonCodingTg, nonCodingTgHasOst
			     , isWellCoding /* Danielle's selection for Worf5, june 28 2005 */
			     , type[1000], myType[1000]
			     , codingComplete, codingNonComplete
			     , codingCompleteHasOst,  codingNonCompleteHasOst 
			     , codingCompleteHasOstPartial,  codingNonCompleteHasOstPartial 
			     , codingCompleteHasOstComplete,  codingNonCompleteHasOstComplete
			     , codingCompleteNoOst,  codingNonCompleteNoOst 
			     , codingCompleteNoOstInGene,  codingNonCompleteNoOstInGene
			     , codingCompleteNoOstTestedNegative,  codingNonCompleteNoOstTestedNegative 
			     , codingCompleteHasOstCompleteNoHole
			     , codingNonCompleteHasOstCompleteNoHole
			     ;DICT *dict
                             ; Array gg
			     ; } WORF ;

typedef struct ggStruct { int type, gene ; } GG ;
static BOOL ggAdd (WORF *w, char *type, AC_OBJ gene) ;

/*****************************************************************************/
/*****************************************************************************/
/* Report on success or failure of oligo choice in the worm ORFeome project */

/* uitility returns the IntMap values
   if *mapp == NULL, it is filled
   else gene->IntMap must == *mapp is inforced
*/

static BOOL worfIntMap (AC_OBJ gene, AC_OBJ *mapp, int *m1p, int *m2p, AC_HANDLE h)
{
  AC_TABLE mm = ac_tag_table (gene, "IntMap", h) ;
  int m1, m2 ;
  AC_OBJ mymap = 0 ;

  if (mm)
    {
      mymap = ac_table_obj (mm, 0, 0, h) ;
      m1 = ac_table_int (mm, 0, 1, 0) ;
      m2 = ac_table_int (mm, 0, 2, 0) ;
      
      if (m2 &&
	  ((! (*mapp)) || ac_obj_equal (*mapp, mymap))
	  )
	{
	  if  (! (*mapp)) *mapp = mymap ;
	  *m1p = m1 ; *m2p = m2 ;
	  return TRUE ;
	}
    }
  return FALSE ;
} /* worfIntMap */

/*****************************************************************************/

static int worfSuccess (AC_OBJ gene, AC_OBJ map, int r1, int r2, AC_HANDLE h)
{
  /*   char *types [] = { "Success", "Failed", "Not tested", "Modif", 0 } ; */
  AC_OBJ ost = 0 ;
  AC_ITER iter = 0 ;
  int nn = 2 ; /* not tested */
  int u1, u2 ;

  iter = ac_objquery_iter (gene, ">OST", h) ;

  while (ac_free (ost), ost = ac_iter_obj (iter))
    {
      if (!worfIntMap (ost, &map, &u1, &u2, h))
	continue ;
      
      if (
	   r1 > u1 - 10 && r1 < u1 + 10 &&
	   r2 > u2 - 10 && r2 < u2 + 10 
	   )
	{
	  if (ac_has_tag (ost, "Not_amplified"))
	    { nn = 1 ;  /* FAILED, continue looping */ }
	  if (ac_has_tag (ost, "amplified"))
	    { nn = 0 ; goto done ; /* SUCCESS */ }
	}
    }

 done:
  ac_free (ost) ;
  return nn ;
} /* worfSuccess  */

/*****************************************************************************/

static int worfOneOligoMrna (AC_OBJ gene, AC_OBJ mrna, WORF *w, BOOL showNegative)
{
  AC_HANDLE h = handleCreate () ;
  AC_OBJ prod  ;
  AC_OBJ map = 0 ;
  AC_TABLE prods = ac_tag_table (mrna, "Product", h) ;
  AC_TABLE spls = 0 ;
  int nn = 0 , ir, jr, iProdLength, pLen ;
  int sType = 1 ;
  int a1, a2, x1, x2 ; /* splicing coords of the mRNA */
  int m1, m2 ; /*  IntMap du mRNA */
  int p1, p2 ; /*  position du producit sur le mRNA */
  int q1, q2 ; /*  position de p1, p2 remappe sur le genome */
  int r1, r2 ; /*  intmap des bouts de la protein */

  if (! worfIntMap (mrna, &map, &m1, &m2, h))
    goto done ;

  spls = ac_tag_table (mrna, "Splicing", h) ;
  for (ir = 0 ; spls && prods && ir < prods->rows ; ir++)
    {
      prod = ac_table_obj (prods, ir, 0, h) ;
      if (! ac_has_tag (prod, "Best_product") ||
	  ! ac_has_tag (prod, "COOH_complete")
	  )
	continue ;
      p1 = ac_table_int (prods, ir, 1, -1) ;
      p2 = ac_table_int (prods, ir, 2, -1) ;
      if (p1 < 0 || p2 < 0 || p2 < p1 + 10)
	continue ;
      pLen = (p2 - p1 + 1 - 3)/3 ;
      iProdLength = 0 ;
      if (pLen > 80) iProdLength = 1 ;
      if (pLen > 200) iProdLength = 2 ;
      if (pLen > 300) iProdLength = 3 ;

      /* get the spliced values of p1 p2 */
      q1 = q2 = -1 ;
      for (jr = 0 ; q2 < 0 && jr < spls->rows ; jr++)
	{
	  a1 = ac_table_int (spls, jr, 0, -1) ;
	  a2 = ac_table_int (spls, jr, 1, -1) ;
	  x1 = ac_table_int (spls, jr, 2, -1) ;
	  x2 = ac_table_int (spls, jr, 3, -1) ;
	  if (strstr (ac_table_printable (spls, jr, 4, ""), "xon") &&
	      x2 > 0)
	    {
	      if (q1 < 0 &&
		  x1 <= p1 && x2 >= p1)
		q1 = a1 + p1 - x1 ;
	      if (q2 < 0 &&
		  x1 <= p2 && x2 >= p2)
		q2 = a1 + p2 - x1 ;
	    }
	}
      if (q2 < 0)
	continue ;
      if (m1 < m2) { r1 = m1 + q1 - 1 ; r2 = m1 + q2 - 1 ; }
      else { r1 = m1 - q1 + 1 ; r2 = m1 - q2 + 1 ; }
      sType = worfSuccess (gene, map, r1, r2, h) ;
      
      if (showNegative)
	{ if (sType == 1) nn = 1 ; }
      else
	{
	  w->type[0 + 10*sType + iProdLength]++ ;
	  w->hasTG++ ;  ggAdd (w, "hasTG", gene) ;

	  nn++ ;
	}
    }
  
 done:
  ac_free (h) ;
  return nn ;
} /* worfOneOligoMrna */

/*****************************************************************************/

static int worfOneOligoPg (AC_OBJ gene, AC_OBJ pg, WORF *w, BOOL showNegative)
{
  AC_HANDLE h = handleCreate () ;
  BOOL isHmm = strncasecmp (ac_name (pg), "hmm", 3) ? FALSE : TRUE ;
  int nn, iProdLength  = 0 ; 
  AC_TABLE exons = ac_tag_table (pg, "Source_Exons", h) ;
  AC_OBJ map = 0 ;
  int pLen = 0, x1, x2, ir ;
  int gType = 0, sType = 1 ;
  int m1, m2 ; /*  IntMap du pg */
  int p1, p2 ; /*  position du produit sur le pg */
  int r1, r2 ; /*  intmap des bouts de la protein */

  nn = 0 ;
  
  if (! worfIntMap (pg, &map, &m1, &m2, h))
    goto done ;

  for (p1 = p2 = -1, ir = 0 ; exons && ir < exons->rows ; ir++)
    {
      x1 = ac_table_int (exons, ir, 0, -1) ;
      x2 = ac_table_int (exons, ir, 1, -1) ;
      if (isHmm && ! strstr (ac_table_printable (exons, ir, 2, ""), "xon"))
	continue ;
      if (p1 < 0) p1 = x1 ;
      p2 = x2 ;
      pLen += x2 - x1 + 1 ;      
    }
   
  iProdLength = 0 ;
  pLen /= 3 ; if (! isHmm) pLen-- ; /* exclude the stop in genefinder case */
  if (pLen > 80) iProdLength = 1 ;
  if (pLen > 200) iProdLength = 2 ;
  if (pLen > 300) iProdLength = 3 ;

  if (m1 < m2) { r1 = m1 + p1 - 1 ; r2 = m1 + p2 - 1 ; }
  else { r1 = m1 - p1 + 1 ; r2 = m1 - p2 + 1 ; }
  sType = worfSuccess (gene, map, r1, r2, h) ;


  if (showNegative)
    { if (sType == 1) nn = 1 ; }
  else
    {
      if (isHmm){ gType = 1 ;  w->hasHmm++ ;   ggAdd (w, ">hasHmm", gene) ;}
      else { gType = 2 ; w->hasGF++ ;  ggAdd (w, "hasGF", gene) ; }
      w->type[100 * gType + 10*sType + iProdLength]++ ;
    }
 done:
  ac_free (h) ;
  return nn ;
} /* worfOneOligoPg */

/*****************************************************************************/

static int pgTestedNegative (AC_OBJ gene, AC_KEYSET ks)
{
  int neg = 0 ;
  AC_OBJ pg = 0;
  AC_ITER iter = ac_keyset_iter (ks, TRUE, 0) ;

  while (ac_free (pg), pg = ac_iter_obj (iter))
    if (worfOneOligoPg (gene, pg, 0, TRUE))
      neg++ ;
  
  ac_free (iter) ;
  return neg ;
}  /* pgTestedNegative */

/*****************************************************************************/

static void worfOneOligo (AC_OBJ gene, WORF *w, int type)
{
  AC_HANDLE h = handleCreate () ;
  AC_ITER iter ;
  AC_OBJ obj = 0 ;
  int nn = 0 ;
  
  w->nGenes++ ;
  
  /* deal with the cDNA supported cases */
  iter = ac_objquery_iter (gene, "{Follow transcribed_gene} SETMINUS {Follow transcribed_gene ; >to_be_fused_with  !Intron_boundaries}  ; Follow mRNA", h) ;
  while (ac_free (obj), obj = ac_iter_obj (iter))
    nn += worfOneOligoMrna (gene, obj, w, FALSE) ;
  
  /* deal with the models */
  iter = ac_objquery_iter (gene, "Follow genefinder", h) ;
  while (ac_free (obj), obj = ac_iter_obj (iter))
    nn += worfOneOligoPg (gene, obj, w, FALSE) ;
  
  ac_free (h) ;
} /* worfOneOligo */


/*****************************************************************************/

static void worfShowOligo (WORF *w, int type)
{
  vTXT txt = vtxtCreate () ;
  int iGeneType, iProdLength, iType, nn = 0 ;
  char *geneTypes[] = { "AceView", "hmm", "wormbase", 0 } ;
  char *prodLengths[] = { "<80 AA", "80-200 AA", "200-300 AA", ">300 AA", 0} ;
  char *types [] = { "Success", "Failed", "Not tested", "Modif", 0 } ;

  vtxtPrintf (txt, "%12s%12s", "TYPE", "ORF") ;
  for (iProdLength = 0 ; prodLengths[iProdLength] ; iProdLength++)
    vtxtPrintf (txt, "%12s", prodLengths[iProdLength]) ;
  vtxtPrint (txt, "\n") ;
  for (iGeneType = 0 ; geneTypes [iGeneType] ; iGeneType++)
    {
      vtxtPrint (txt, "\n") ;
      for (iType = 0 ; types[iType] ; iType++)
	{ 
	  vtxtPrintf (txt, "%12s%12s", geneTypes [iGeneType], types[iType]) ;
	  for (iProdLength = 0 ; prodLengths[iProdLength] ; iProdLength++)
	    {
	      nn =  w->type[100*iGeneType + 10*iType + iProdLength] ;
	      vtxtPrintf (txt, "%12d", nn) ;
	    }
	vtxtPrint (txt, "\n") ;
	}
    } 
  vtxtPrint (txt, "\n") ;
  freeOut (vtxtPtr (txt)) ;
  vtxtDestroy (txt) ;
}  /* worfShowOligo */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/



static BOOL worfGetBestMP (AC_OBJ tg, AC_OBJ *bestMrnap, AC_OBJ *bestProductp
			   , int *p1p, int *p2p
			   , int type, int *geneType, int *intronType)
			   
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER mrnas = ac_objquery_iter (tg, ">mRNA ; IS *.a", h) ;
  AC_ITER products = 0 ;
  int ir, nIntron = 0, nIntronOut = 0 ;
  int p1, p2, n ;
  AC_OBJ bestMrna = 0, bestProduct = 0 ;

  if (mrnas)
    bestMrna = ac_iter_obj (mrnas) ;

  if (!bestMrna)
    {
      if (debug) fprintf (stderr, "%s is non coding\n", ac_name(tg)) ;
      return FALSE ;
    }
   
  products = ac_objquery_iter (bestMrna, "> product ; good_product && best_product", 0) ;
  if (products)
    bestProduct = ac_iter_obj (products) ;

  p1 = p2 = 0 ;

  if (bestProduct)
    {
      AC_TABLE tt = ac_tag_table (bestMrna, "Product", 0) ;
      int pLen = 0 ;

      for (ir = 0 ; tt && ir < tt->rows ; ir++)
	if (!strcmp (ac_table_printable (tt, ir, 0, ""), ac_name(bestProduct)))
	  { p1 = ac_table_int (tt, ir, 1, 0) ; p2 = ac_table_int (tt, ir, 2, 0) ; }
      ac_free (tt) ;
      pLen = p2 - p1 + 1 ; pLen /= 3 ;
      tt = ac_tag_table (bestMrna, "Splicing", 0) ;

      for (ir = 0 ; p1 > 0 && p2 > 0 && tt && ir < tt->rows ; ir++)
	if (strstr (ac_table_printable (tt, ir, 4, ""), "tron") &&
	    strstr (ac_table_printable (tt, ir, 5, ""), "_ag"))
	  {
	    nIntron++ ;
	    if (p1 > ac_table_int (tt, ir, 3, 0) ||
		p2 <  ac_table_int (tt, ir, 2, 0))
	      nIntronOut++ ;
	  }
      ac_free (tt) ;

      if ( 0 &&
	   (
	    (nIntronOut >= 3 && pLen < 300) ||
	    (nIntronOut >= 2 && pLen < 200) ||
	    (!nIntron >= 2 && pLen < 80) 
	    )	  
	   )
	{
	  if (debug) 
	    fprintf (stderr, "%s non coding pLen=%d aa, nIntron=%d, nIntronOut = %d\n"
			      , ac_name(tg), pLen, nIntron, nIntronOut) ;
	  ac_free (bestProduct) ;
	}
    }

  if (bestProduct)
    {
      *bestMrnap = bestMrna ; *bestProductp = bestProduct ;
      *p1p = p1 ; *p2p = p2 ;
      n = p2 - p1 ;
      if (n < 240)
	*geneType = 0 ;
      else if (n < 600)
	*geneType = 1 ;
      else if (n < 900)
	*geneType = 2 ;
      else
	*geneType = 3 ;

      if (nIntron > 0 && nIntronOut == 0)
	*intronType = 0 ;
      else if (nIntron >= 2 && nIntronOut == 1)
	*intronType = 1 ;
      else if (nIntron > 2 && nIntronOut == 2)
	*intronType = 2 ;
      else if (nIntron >= 3 && nIntronOut >= 3)
	*intronType = 3 ;
      else if (!nIntron)
	*intronType = 4 ;
      else if (nIntron == nIntronOut )
	*intronType = 3 ;
      return TRUE ;
    }

  ac_free (bestMrna) ;
  if (debug) fprintf (stderr, "%s is non coding\n", ac_name(tg)) ;
  return FALSE ;
}

/*****************************************************************************/

static int worfGetOst (AC_OBJ bestMrna, AC_OBJ bestProduct, int p1, int p2, int type)
{
  AC_HANDLE h = handleCreate () ;
  AC_TABLE clones, coveringReads ;
  AC_TABLE tt = ac_tag_table (bestMrna, "Constructed_from", h) ;
  AC_TABLE reads = 0 ;
  int ir, jr, kr, x1, x2, y1, y2, ok = 0, ok2, oldOk ;
  AC_OBJ read, read2, clo, bestClone = 0 ;
  BOOL hasMgc = FALSE , hasHinv = FALSE ;

  clones = ac_tag_table (bestMrna, "cdna_clone", h) ;

  if (type == 4) /* other */
    {
      for (jr = 0 ; clones && !hasMgc && jr < clones->rows ; jr++)
	{
	  
	  clo = ac_table_obj (clones, jr, 0, h) ;
	  if (ac_has_tag (clo, "MGC"))
	    hasMgc = TRUE ;
	  if (ac_has_tag (clo, "MGC") || ! ac_has_tag (clo, "HINV_libs"))
	    hasHinv = TRUE ;
	}
      if (hasMgc || ! hasHinv)
	{ ok = 0 ; goto done ; }
    }

  if (type == 5 || type == 6)
    coveringReads = ac_tag_table (bestProduct, "Covered_by", h) ;
  if (type == 6) /* other exact, no exact MGC */
    {
      BOOL isCovering = FALSE, mgcCovering = FALSE ;
      for (ir = 0 ; coveringReads && ! mgcCovering && ir < coveringReads->rows ; ir++)
	if (ac_table_int (coveringReads, ir, 2, -1) != 0 ||
	    ac_table_int (coveringReads, ir, 4, -1) != 0)
	  {
	    read = ac_table_obj (coveringReads, ir, 0, h) ;
	    clo = ac_tag_obj (read, "cDNA_clone", h) ;
	    if (ac_has_tag (clo, "MGC"))
	      mgcCovering = TRUE ;
	    if (ac_has_tag (clo, "MGC") || ! ac_has_tag (clo, "HINV_libs"))
	      isCovering = TRUE ;
	  }
      if (!isCovering || mgcCovering)
	{ ok = 0 ; goto done ; }
    }

  for (jr = 0 ; clones && jr < clones->rows ; jr++)
    {
      clo = ac_table_obj (clones, jr, 0, h) ;
      if (0 && type == 1 && ! ac_has_tag (clo, "Worf1") && ! ac_has_tag (clo, "Worf3"))
	continue ;
      if (1 && type == 1 && strncasecmp (ac_name (clo), "ost", 3))
	continue ;
      if (0 && type == 1 &&  ac_has_tag (clo, "suspected_internal_deletion"))
	continue ;
      if (type == 2 && (ac_name (clo))[0] != 'y')
	continue ;
      if (type == 3 && ! ac_has_tag (clo, "MGC"))
	continue ;
      if (type == 4 && ( ac_has_tag (clo, "MGC") || ! ac_has_tag (clo, "HINV_libs")))
	continue ;
      if (type == 5 && ! ac_has_tag (clo, "MGC"))
	continue ;
      if (type == 6 && ( ac_has_tag (clo, "MGC") || ! ac_has_tag (clo, "HINV_libs")))
	continue ;
      reads = ac_tag_table (clo, "Read", h) ;
      
      y1 = 99999999; y2 = 0 ;
      for (ok2 = kr = 0 ; reads && kr < reads->rows ; kr++)
	{
	  read = ac_table_obj (reads, kr, 0, h) ;
	  if (type == 5 || type == 6)
	    {
	      BOOL isCovering = FALSE ;
	      for (ir = 0 ; coveringReads && ir < coveringReads->rows ; ir++)
		{
		  read2 = ac_table_obj (coveringReads, ir, 0, h) ;
		  if (!ac_obj_equal (read, read2))
		    continue ;
		  if (ac_table_int (coveringReads, ir, 2, -1) == 0 &&
		      ac_table_int (coveringReads, ir, 4, -1) == 0)
		    isCovering = TRUE ;
		}
	      if (!isCovering)
		continue ;
	    }
	  for (ir = 0 ; tt && ir < tt->rows ; ir++)
	    {
	      read2 = ac_table_obj (tt, ir, 2, h) ;
	      if (!ac_obj_equal (read, read2))
		continue ;
	      ok2 |= 1 ;
	      x1 = ac_table_int (tt, ir, 0, 0) ;
	      x2 = ac_table_int (tt, ir, 1, 0) ;
	      if (x1 <= p1 + 3 && x2 >= p1)
		{ ok2 |= 2 ; if (x2 > y2) y2 = x2 ; }
	      if (x1 <= p2 && x2 >= p2 - 3)
		{ ok2 |= 4 ; if (x1 < y1) y1 = x1 ; }
	    }
	}
      oldOk = ok ;
      switch (ok2)
	{
	case 1: 
	  ok |= ok2 ; break ;
	case 3: /* covers the ATG */
	  if (oldOk == 3)
	    bestClone = clo ; /* better than a 3p complete */
	case 5:
	  ok |= 3 ; break ; /* partial on any side */
	case 7: 
	  if (y1 < y2) ok2 |= 8 ; /* complete and no hole */
	  ok |= ok2 ; /* complete */
	  break ; 
	}
      if (ok > oldOk)
	bestClone = clo ;
    }
 done:
  if (bestClone && (ok == 3 || ok == 7) && type == 2)
    /* dictAdd (w->dict, ac_name(bestClone), 0)  */ ;
  ac_free (h) ;
  return ok ;
}

/*****************************************************************************/

static void worfOneStats (AC_OBJ gene, WORF *w, int type)
{
  AC_HANDLE h = handleCreate () ;
  AC_OBJ bestMrna = 0, bestProduct = 0 ;
  int neg1, neg2, p1, p2, geneType, intronType ;
  AC_KEYSET ks1, ks2 ;
  AC_OBJ tg = ac_tag_obj (gene, "Transcribed_gene", h) ;

  w->nGenes++ ;  ggAdd (w, "nGenes", gene) ;

  ks1 = ac_objquery_keyset (gene, "! transcribed_gene ; > genefinder ; IS hmm*", h) ;
  ks2 = ac_objquery_keyset (gene, "! transcribed_gene ; > genefinder ; ! IS hmm*", h) ;

  if (ac_keyset_count (ks1))
    { w->hasHmm++ ;  ggAdd (w, "hasHmm", gene) ;}

  if (ac_keyset_count (ks2))
    { w->hasGF++ ;  ggAdd (w, "hasGF", gene) ;}

  neg1 = pgTestedNegative (gene, ks1) ;
  neg2 = pgTestedNegative (gene, ks2) ;
  if (ac_keyset_count (ks1) && ac_keyset_count (ks2))
    {
      w->hasGFandHMM++ ;
      ggAdd (w, "hasGFandHMM", gene) ;

      if (neg1) 
	{
	  w->hasGFandHMM_hmmNegative++ ;
	  ggAdd (w, "hasGFandHMM_hmmNegative", gene) ;
	}
      if (neg2)
	{
	  w->hasGFandHMM_gfNegative++ ;
	  ggAdd (w, "hasGFandHMM_gfNegative", gene) ;
	}
    } 
  if (ac_keyset_count (ks1) && !ac_keyset_count (ks2))
    {
      w->hasOnlyHMM++ ;
      ggAdd (w, "hasOnlyHMM", gene) ;
      
      if (neg1)
	{
	  w->hasOnlyHMM_hmmNegative++ ;
	  ggAdd (w, "hasOnlyHMM_hmmNegative", gene) ;
	}
    }
  if (!ac_keyset_count (ks1) && ac_keyset_count (ks2))
    {
      w->hasOnlyGF++ ;
      ggAdd (w, "hasOnlyGF", gene) ;

      if (neg2)
	{
	  w->hasOnlyGF_gfNegative++ ;
	  ggAdd (w, "hasOnlyGF_gfNegative", gene) ;
	}
    }
  if (ac_keyset_count (ks1) || ac_keyset_count (ks2))
    {
      w->hasGForHMM++ ;
      ggAdd (w, "hasGForHMM", gene) ;

      if (neg1)
	{
	  w->hasGForHMM_hmmNegative++ ;
	  ggAdd (w, "hasGForHMM_hmmNegative", gene) ;
	}
      if (neg2)
	{
	  w->hasGForHMM_gfNegative++ ;
	  ggAdd (w, "hasGForHMM_gfNegative", gene) ;
	}
    }

  if (tg)
    {
      w->hasTG++ ; ggAdd (w, "hasTG", gene) ;
      if (ac_keyset_count (ac_objquery_keyset (tg, ">cdna_clone OST*", h)))
	 {
	   w->hasOST++ ;
	   ggAdd (w, "hasOST", gene) ;
	 }
      if (!worfGetBestMP (tg, &bestMrna, &bestProduct, &p1, &p2, type, &geneType, &intronType))
	{
	  w->nonCodingTg++ ;
	  ggAdd (w, "nonCodingTg", gene) ;

	  if (ac_keyset_count (ac_objquery_keyset (tg, ">cdna_clone ;IS  OST*", h)))
	    {
	      w->nonCodingTgHasOst++ ;
	      ggAdd (w, "nonCodingTgHasOst", gene) ;
	    }
	}
      else
	{
	  BOOL isComplete = ac_has_tag (bestProduct, "Complete") ;
	  int firstAtg = ac_tag_int (bestProduct, "First_ATG", -1) ;

	  if (firstAtg < 1) firstAtg = 1;

	  if (isComplete)
	    { w->codingComplete++ ;  ggAdd (w, "codingComplete", gene) ; }

	  else
	    { w->codingNonComplete++ ; ggAdd (w, "codingNonComplete", gene) ; }

	  switch (worfGetOst (bestMrna, bestProduct, p1 + 3 * (firstAtg - 1)+ 9, p2 - 9, type))
	    {
	    case 0:
	      if (isComplete)
		{
		  { w->codingCompleteNoOst++ ;  ggAdd (w, "codingCompleteNoOst", gene) ;}
		  
		  if (!ac_keyset_count (ac_objquery_keyset (tg, ">cdna_clone OST*", h)))
		    { 
		      w->codingCompleteNoOstInGene++ ; 
		      ggAdd (w, "codingCompleteNoOstInGene", gene) ;
		    }
		  if (worfOneOligoMrna (gene, bestMrna, 0, TRUE))
		   {
		     w->codingCompleteNoOstTestedNegative++ ;  
		     ggAdd (w, "codingCompleteNoOstTestedNegative", gene) ;
		   }
		}
	      else
		{
		  w->codingNonCompleteNoOst++ ;  ggAdd (w, "codingNonCompleteNoOst", gene) ;

		  if (!ac_keyset_count (ac_objquery_keyset (tg, ">cdna_clone OST*", h)))
		    { 
		      w->codingNonCompleteNoOstInGene++ ; 
		      ggAdd (w, "codingNonCompleteNoOstInGene", gene) ;
		    }
		  if (worfOneOligoMrna (gene, bestMrna, 0, TRUE))
		   { 
		     w->codingNonCompleteNoOstTestedNegative++ ;  
		     ggAdd (w, "codingNonCompleteNoOstTestedNegative", gene) ;
		   }
		}
	      break ;
	    case 1:
	    case 3:
	    case 5:
	      if (isComplete)
		{
		  w->codingCompleteHasOst++ ; 
		  ggAdd (w, "codingCompleteHasOst", gene) ;

		  w->codingCompleteHasOstPartial++ ; 
		  ggAdd (w, "codingCompleteHasOstPartial", gene) ;
		}
	      else
		{
		  w->codingNonCompleteHasOst++ ;  
		  ggAdd (w, "codingNonCompleteHasOst", gene) ;

		  w->codingNonCompleteHasOstPartial++ ;
		  ggAdd (w, "codingNonCompleteHasOstPartial", gene) ;
		}
	      break ;
	    case 7:
	      w->myType[10*geneType + intronType]++ ;
	      if (isComplete)
		{
		  w->codingCompleteHasOst++ ;
		  ggAdd (w, "codingCompleteHasOst", gene) ;

		  w->codingCompleteHasOstComplete++ ; 
		  ggAdd (w, "codingCompleteHasOstComplete", gene) ;
		}
	      else
		{
		  w->codingNonCompleteHasOst++ ;
		  ggAdd (w, "codingNonCompleteHasOst", gene) ;

		  w->codingNonCompleteHasOstComplete++ ; 
		  ggAdd (w, "codingNonCompleteHasOstComplete", gene) ;
		}
	      break ;
	    case 15:
	      w->myType[10*geneType + intronType]++ ;
	      if (isComplete)
		{
		  w->codingCompleteHasOst++ ; 
		  ggAdd (w, "codingCompleteHasOst", gene) ;

		  w->codingCompleteHasOstComplete++ ;
		  ggAdd (w, "codingCompleteHasOstComplete", gene) ;

		  w->codingCompleteHasOstCompleteNoHole++ ; 
		  ggAdd (w, "codingCompleteHasOstCompleteNoHole", gene) ;
		}
	      else
		{
		  w->codingNonCompleteHasOst++ ; 
		  ggAdd (w, "codingNonCompleteHasOst", gene) ;
		  
		  w->codingNonCompleteHasOstComplete++ ; 
		  ggAdd (w, "codingNonCompleteHasOstComplete", gene) ;

		  w->codingNonCompleteHasOstCompleteNoHole++ ;
		  ggAdd (w, "codingNonCompleteHasOstCompleteNoHole", gene) ;
		}
	      break ;
	    }
	  w->type[10*geneType + intronType]++ ;
	  if (
	      (geneType == 0 && intronType == 0) ||
	      (geneType == 1 && intronType == 0) ||
	      (geneType == 2 && intronType == 0) ||
	      (geneType == 3 && intronType == 0) ||
	      (geneType == 1 && intronType == 1) ||
	      (geneType == 2 && intronType == 1) ||
	      (geneType == 3 && intronType == 1) ||
	      (geneType == 3 && intronType == 2) ||
	      (geneType == 3 && intronType == 3) ||
	      (geneType == 2 && intronType == 4) ||
	      (geneType == 3 && intronType == 4)
	      )
	    ggAdd (w, "isWellCoding", gene) ;
	}
    }
  ac_free (bestMrna) ;   /* not allocated on h */
  ac_free (bestProduct) ;/* not allocated on h */
  ac_free (h) ;

  /* categories utiles
     
   * 1: has_ost
   *    1a: supported by other clones
   *       1a1: good ATG (modulo first A)
   *          1a11:  good ATG and we know the mrna is complete
   *       1a2: good stop
   *          1a21: good start + stop
   *       1a3: has mutations
   *       1a4: several OST in this cluster

   *    1b: supported by prediction in WS140
   *       1b1: good ATG (modulo first A)
   *       1b2: good stop
   *          1b21: good start + stop
   *       1b3: has mutations
   *       1a4: several OST in this cluster

   *    1c: not matching any clone or any prediction
   *       1c1: uniquely mapped here
   *       1c4: several OST in this cluster

   * 2: no_ost
   *    2a: supported by other clones
   *       2a1: at least 200 AA
   *       2a2: at least one intron
   *    2b: just predicted
   *       2a1: at least 200 AA
   *       2a2: at least one intron
   */

}

/*****************************************************************************/

static void worfShow (WORF *w, int type)
{
  vTXT txt = vtxtCreate () ;
  char *ttt[] = { "OST" , "YK", "MGC", "Other", "MGC Exact", "Other Exact" } ;
  char *tt = ttt[type-1] ;
  int jj = 0 ;

  vtxtPrintf (txt, "\t%d protein coding genes are defined in AceView\n", w->nGenes) ;
  vtxtPrintf (txt, "\t  %d have some cDNA support\n", w->hasTG) ;
  vtxtPrintf (txt, "\t      %d have some OST support\n", w->hasOST) ;
if (0)
  {
    vtxtPrintf (txt, "\t    %d seem non coding: (no intron && < 80aa) \n", w->nonCodingTg) ;
    vtxtPrintf (txt, "\t                    or (nIntronOutCDS>=2 && < 200aa)\n") ;
    vtxtPrintf (txt, "\t                    or (nIntronOutCDS>=3 && < 300aa)\n") ;
  }
else
  {
    vtxtPrintf (txt, "\t    %d seem non coding: no variant .a encoding a best product\n", w->nonCodingTg) ;
  }
  vtxtPrintf (txt, "\t      %d have OST support\n", w->nonCodingTgHasOst) ;
  vtxtPrintf (txt, "\t      %d do not have OST support\n", w->nonCodingTgHasOst - w->nonCodingTgHasOst) ;
  vtxtPrintf (txt, "\t    %d variant .a encodes a  complete product\n", w->codingComplete) ;
  vtxtPrintf (txt, "\t      %d Coding Complete with %s\n", w->codingCompleteHasOst, tt) ;
  vtxtPrintf (txt, "\t        %d Coding Complete with Partial %s\n", w->codingCompleteHasOstPartial, tt) ;
  vtxtPrintf (txt, "\t        %d Coding Complete with Complete %s\n", w->codingCompleteHasOstComplete, tt) ;
  vtxtPrintf (txt, "\t              %d and no hole over the CDS\n", w->codingCompleteHasOstCompleteNoHole) ;
  vtxtPrintf (txt, "\t      %d Coding Complete with no  %s\n", w->codingCompleteNoOst, tt) ;
  vtxtPrintf (txt, "\t        %d Coding Complete with no  %s in gene\n", w->codingCompleteNoOstInGene, tt) ;
  vtxtPrintf (txt, "\t        %d already tested negative\n", w->codingCompleteNoOstTestedNegative) ;
  vtxtPrintf (txt, "\t    %d variant .a encodes a possibly incomplete product\n",  w->codingNonComplete) ;
  vtxtPrintf (txt, "\t      %d Coding incomplete with %s\n", w->codingNonCompleteHasOst, tt) ;
  vtxtPrintf (txt, "\t        %d Coding incomplete with Partial %s\n", w->codingNonCompleteHasOstPartial, tt) ;
  vtxtPrintf (txt, "\t        %d Coding incomplete with Complete %s\n", w->codingNonCompleteHasOstComplete, tt) ;
  vtxtPrintf (txt, "\t              %d and no hole over the CDS\n", w->codingNonCompleteHasOstCompleteNoHole) ;

  vtxtPrintf (txt, "\t      %d Coding incomplete with no %s\n", w->codingNonCompleteNoOst, tt) ;
  vtxtPrintf (txt, "\t        %d already tested negative\n", w->codingNonCompleteNoOstTestedNegative) ;
  vtxtPrintf (txt, "\t        %d Coding incomplete with no %s in gene\n", w->codingNonCompleteNoOstInGene, tt) ;
  vtxtPrintf (txt, "\t  %d are pure ab initio models\n", w->nGenes - w->hasTG) ;

  vtxtPrintf (txt, "\t    %d are predicted by hmm OR by genefinder\n", w->hasGForHMM) ;
  vtxtPrintf (txt, "\t        %d  at least one gf was tested negative\n", w->hasGForHMM_gfNegative) ;
  vtxtPrintf (txt, "\t        %d  at least one hmm was tested negative\n", w->hasGForHMM_hmmNegative) ;
  vtxtPrintf (txt, "\t      %d are predicted only by hmm\n", w->hasOnlyHMM) ;
  vtxtPrintf (txt, "\t        %d  at least one hmm was tested negative\n", w->hasOnlyHMM_hmmNegative) ;
  vtxtPrintf (txt, "\t      %d are predicted only by genefinder\n", w->hasOnlyGF) ;
  vtxtPrintf (txt, "\t        %d  at least one gf was tested negative\n", w->hasOnlyGF_gfNegative) ;
  vtxtPrintf (txt, "\t      %d are predicted by hmm and by genefinder\n", w->hasGFandHMM) ;
  vtxtPrintf (txt, "\t        %d  at least one gf was tested negative\n", w->hasGFandHMM_gfNegative) ;
  vtxtPrintf (txt, "\t        %d  at least one hmm was tested negative\n", w->hasGFandHMM_hmmNegative) ;

  for (jj = 0 ; jj < 2 ; jj++)
  {
    if (jj == 1)
      vtxtPrintf (txt
		, "\nTypes of products supported by an %s clone\n", tt) ;
    else
      vtxtPrintf (txt
		  , "\nTypes of the best product of all genes\n") ;
    {
      int ir, jr, tot[100], total, toth, nn ;
      char *gTypes[] = { "<80 AA", "80-200 AA", "200-300 AA", ">300 AA", 0} ;
      char *iTypes [] = { "All introns", ">=2 introns", ">=3 introns", "all out or", "no intron", 0};
      char *iTypes2 [] = { "inside CDS", "1 out CDS", "2 out CDS", ">=3 out CDS", " ", 0};
      vtxtPrintf(txt, "%014s", " ") ;
      for (ir = 0 ; iTypes[ir] ; ir++)
	{
	  vtxtPrintf(txt, "%014s", iTypes[ir]) ;
	}
      vtxtPrintf (txt, "\n") ; 
      vtxtPrintf(txt, "%014s", "length") ;
      for (ir = 0 ; iTypes[ir] ; ir++)
	{
	  vtxtPrintf(txt, "%014s", iTypes2[ir]) ;
	  tot[ir] = 0 ;
	}
      vtxtPrintf (txt, "%014s\n", "Total") ;
      total = 0 ;
      for (jr = 0 ; gTypes[jr] ; jr++)
	{
	  toth = 0 ;
	  vtxtPrintf(txt, "%014s", gTypes[jr]) ;
	  for (ir = 0 ; iTypes[ir] ; ir++)
	    {
	      nn =  jj ? w->myType[10*jr + ir] : w->type[10*jr + ir] ;
	      vtxtPrintf(txt, "%14d", nn) ;
	      tot [ir] += nn ;
	      toth += nn ;
	      total += nn ;
	    }
	  vtxtPrintf (txt, "%14d\n", toth) ;
	}
      vtxtPrintf(txt, "%014s", "Total") ;
      for (ir = 0 ; iTypes[ir] ; ir++)
	vtxtPrintf(txt, "%14d", tot[ir]) ;
	vtxtPrintf(txt, "%14d", total) ;
      vtxtPrintf (txt, "\n") ;
    }
  }
  freeOut (vtxtPtr (txt)) ;
  vtxtDestroy (txt) ;
} /* worfShow */

/*****************************************************************************/

static void worfShowTerse (WORF *w, int type)
{
  vTXT txt = vtxtCreate () ;
  int jj = 0 ;

  for (jj = 0 ; jj < 2 ; jj++)
  {
    {
      int ir, jr, nn ;
      char *gTypes[] = { "<80 AA", "80-200 AA", "200-300 AA", ">300 AA", 0} ;
      char *iTypes [] = { "All introns", ">=2 introns", ">=3 introns", "all out or", "no intron", 0};

      for (jr = 0 ; gTypes[jr] ; jr++)
	{
	  vtxtPrintf (txt, "%d\t%d", jj, jr) ;
	  for (ir = 0 ; iTypes[ir] ; ir++)
	    {
	      nn =  jj ? w->myType[10*jr + ir] : w->type[10*jr + ir] ;
	      vtxtPrintf (txt, "\t%d", nn) ;
	    }
	  vtxtPrint (txt, "\n") ;
	}
    }
  }
  freeOut (vtxtPtr (txt)) ;
  vtxtDestroy (txt) ;
} /* worfShowTerse */

/*****************************************************************************/
/*****************************************************************************/

static int ggOrder (const void *a, const void *b)
{
  const GG *ga = (const GG*)a, *gb = (const GG*) b ;
  int n ;
  
  n = ga->type - gb->type ; if (n) return n ;
  n = ga->gene - gb->gene ;
  
  return n ;
}

/*****************************************************************************/

static BOOL ggAdd (WORF *w, char *type, AC_OBJ gene)
{
  GG *gg ;

  if (!w->dict)
    return FALSE ;
  if (!type || ! *type)
    messcrash ("Bad call to ggAdd type= %s gene = %s\n"
	       , type && *type ? type : "null type"
	       , ac_name (gene)) ;
  
  gg = arrayp (w->gg, arrayMax (w->gg), GG) ;
  dictAdd (w->dict, type, &(gg->type)) ;
  dictAdd (w->dict, ac_name(gene), &(gg->gene)) ;
  
  return TRUE ;
}  /* ggAdd */

/*****************************************************************************/

static void worfShowKeySets (WORF *w, int type)
{
  int i, oldType = -1 ;
  GG *gg ;

  arraySort (w->gg, ggOrder) ;
  
  for (i = 0, gg = arrp (w->gg, 0, GG) ; i < arrayMax (w->gg) ; i++, gg++)
    {
      if (gg->type != oldType)
	{
	  freeOutf ("\nKEYSET \"worf_%s\"\n", dictName (w->dict, gg->type)) ;
	  oldType = gg->type ;
	}
      freeOutf ("Gene \"%s\"\n", dictName (w->dict, gg->gene)) ;
    }
  freeOutf ("\n") ;
}  /* worfShowTerse */

/*****************************************************************************/
/*****************************************************************************/

static int worfStats (AC_DB db, const char *mask, int type, int terse, int oligo, FILE *f)
{
  AC_HANDLE h = handleCreate () ;
  AC_ITER iter ;
  AC_OBJ gene = 0 ;
  WORF worf ;
  
  memset (&worf, 0, sizeof (worf)) ;

  if (f)
    {
      worf.dict = dictHandleCreate (10000, h) ;
      worf.gg = arrayHandleCreate (30000, GG, h) ;
    }

  iter = ac_dbquery_iter (db, messprintf("Find gene %s ; ! non_protein_coding && SMAP && (transcribed_gene || genefinder)", mask ? mask : ""), h) ;

  while (ac_free (gene), gene = ac_next_obj (iter))
    {
      if (oligo)
	worfOneOligo (gene, &worf, type) ;
      else
	worfOneStats (gene, &worf, type) ;
    }
  if (terse)
    worfShowTerse (&worf, type) ;
  else if (oligo)
    worfShowOligo (&worf, type) ;
  else
    worfShow (&worf, type) ;

  if (f)
    {
      int outlevel = freeOutSetFile (f) ;
   
      worfShowKeySets (&worf, type) ;

      freeOutClose (outlevel) ;
    }

  ac_free (h) ;
  return worf.nGenes ;
}  /* worfStats */

/*****************************************************************************/
/*****************************************************************************/
typedef struct worfDumpStruct { int pos ; int ln ; }WD ;

static int wdOrder (const void *va, const void *vb)
{
  const WD *a = (const WD *)va ;
  const WD *b = (const WD *)vb ;
  int n ;

  n = a->ln - b->ln ;
  if (n) return n > 0 ? -1 : 1 ; /* long first */
  n = a->pos - b->pos ;
  return n ;
}

static int worfDump (AC_DB db, const char *mask)
{
  AC_HANDLE h = handleCreate () ;
  AC_ITER iter ;
  AC_TABLE gPrimers = 0, gProd = 0 ;
  AC_OBJ product = 0, mrna = 0, gene = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  Stack s = stackHandleCreate (100000, h) ;
  const char *ccp, *newname ;
  int a1, a2, ir, iwd = 0 ;
  Array wdArray = arrayCreate (10000, WD) ;
  WD *wd ;
  
  stackTextOnly (s) ;
  freeOutf ("#NewName\tGene\tProduct\tLeft\tPrimer1\tRight\tPrimer2\tbp\tmRNA\n") ;
  ccp = messprintf ("Find keyset %s ; expand; >transcribed_gene; >mrna ; IS *.a ; > product best_product ; Complete || Open_length > 300 ; Primers ", mask ? mask : "xxx") ;
  iter = ac_dbquery_iter (db, ccp, h) ;
  while (ac_free (product), product = ac_iter_obj (iter))
    {
      vtxtClear (txt) ;
      gene = ac_tag_obj (product, "GeneBox", h) ;
      mrna = ac_tag_obj (product, "Mrna", h) ;
      if (! gene || ! mrna) 
	continue ;
      newname = ac_tag_printable (gene, "NewName", ac_name(gene)) ;
      ac_free (gPrimers) ;
      gPrimers = ac_tag_table (product, "Primers", h) ;
      if (gPrimers->rows < 2)
	continue ;

      vtxtPrintf (txt, "%s\t%s\t%s\t", ac_name (gene), newname, ac_name(product)) ;
 
      /* left primer temperature then primer sequence */
      ccp = ac_table_printable (gPrimers, 0, 1, "nnn") ;
      vtxtPrintf (txt, "%s\t", ccp) ;
      ccp = ac_table_printable (gPrimers, 0, 0, "nnn") ;
      vtxtPrintf (txt, "%s\t", ccp) ;

      /* right primer temperature then primer sequence */
      ccp = ac_table_printable (gPrimers, 1,  1, "nnn") ;
      vtxtPrintf (txt, "%s\t", ccp) ;
      ccp = ac_table_printable (gPrimers, 1, 0, "nnn") ;
      vtxtPrintf (txt, "%s\t", ccp) ;

      /* dna length and sequence */
      gProd = ac_tag_table (mrna, "Product", h) ;
      a1 = a2 = -1 ;
      for (ir = 0 ; ir < gProd->rows ; ir++)
	{
	  if (strcmp (ac_name (product), ac_table_printable (gProd, ir, 0, "toto")))
	    continue ;
	  a1 = ac_table_int (gProd, ir, 1, -1) ;
	  a2 = ac_table_int (gProd, ir, 2, -1) ;
	}
      if (a1 > 0 && a2 > 0)
	{
	  ccp = ac_zone_dna (mrna, a1, a2, h) ;
	  vtxtPrintf (txt, "%d\t%s\n", strlen (ccp), ccp) ;

	  /* success */
	  wd = arrayp (wdArray, iwd++, WD) ;
	  wd->pos = stackMark (s) ;
	  wd->ln = strlen (ccp) ;
	  pushText (s, vtxtPtr (txt)) ;
	}

    }
  /* export */
  arraySort (wdArray, wdOrder) ;
  for (ir = 0 ; ir < iwd ; ir++)
    {
      wd = arrp (wdArray, ir, WD) ;
      freeOutf ("%s", stackText (s, wd->pos)) ;
    }

  ac_free (h) ;
  return 0 ;
} /* worfDump */

/*****************************************************************************/
/*****************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: worf -db ACEDB (-oligo | -ost | -yk | -dump) [-o outfile] [-terse] [-g gene_mask]\n") ;
  fprintf (stderr, "// Example:  worf -db /home/mieg/yk -yk -g '2b*'\n") ;
  fprintf (stderr, "//    -oligo: say which oligo succeded, failed, could be re predicted\n") ;
  fprintf (stderr, "//    -ost: say which mrna.a are well supported by an OST clone\n") ;
  fprintf (stderr, "//    -yk: say which mrna.a are well supported by an YK clone\n") ;
  fprintf (stderr, "//    -o: export in outfile the details of all keysets\n") ;
  fprintf (stderr, "//    -terse: only applies to ost/yj mode, exports a tab delimited table without legend\n") ;
  fprintf (stderr, "//    -dump: expport a table suitable for philippe lamesh: -g is a worf_keyset_name\n") ;
  exit (1) ;
}

/*****************************************************************************/

int main (int argc, const char **argv)
{
  int nn = 0 ;
  FILE *f = 0 ;
  char *s = "ok" ;
  const char *outfilename = 0 ;
  const char *dbName = 0 ;
  const char *geneMask = "*" ;
  int outlevel = 0 ;
  int type = 0, terse = 0, oligo = 0 ;
  AC_DB db ;

  /* consume optional args */
  outfilename = getArgV (&argc, argv, "-o") ;
  /* read absolute args */
  dbName = getArgV (&argc, argv, "-db") ;
  geneMask = getArgV (&argc, argv, "-g") ;

  if (getArg (&argc, argv, "-oligo"))
    { oligo = 1 ; type = 1 ; }
  else if (getArg (&argc, argv, "-ost"))
    type = 1 ;
  else if (getArg (&argc, argv, "-yk"))
    type = 2 ;
  else if (getArg (&argc, argv, "-dump"))
    type = 10 ;

  if (getArg (&argc, argv, "-terse"))
    { terse = 1 ; if (oligo) type = 0 ; }

  if (!dbName || ! type)
    usage () ;
  db = ac_open_db (dbName, &s);
  if (!db)
    messcrash ("Failed to open db %s, error %s", dbName, s) ;

  if (outfilename)
    f = filopen (outfilename, 0, "w") ;

  outlevel = freeOutSetFile (stdout) ;	

  switch (type)
    {
    case 10: /* dump */ 
      nn = worfDump (db, geneMask) ;
      break ;
    default:
      nn = worfStats (db, geneMask, type, terse, oligo, f) ;
      break ;
    }

  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;

  ac_db_close (db) ;

  fprintf (stderr, "Analysed %d genes\n", nn) ;
  sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
