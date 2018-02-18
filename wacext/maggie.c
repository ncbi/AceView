/*  File: maggie.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2004
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the COMPARATOR genome database package, written by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This code works as a client against an acedb server
 *  and should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  The present code and the acedb schema for this project
 *  are available from http://www.aceview.org/aligold
 */

#include "../wac/ac.h"
#include "keyset.h"
#include <errno.h>
#include "bitset.h"
#include "freeout.h"
#include "vtxt.h"
#include "dict.h"

static void usage (void) ;
typedef struct mmStruct {
  AC_DB db ; 
  BOOL isIntron, isMrna, isCDS, isMapProbe, analyse ;
  BOOL hlsStatsNm, hlsStatsAceView, hitStatsNm, hitStatsAceView, isMaggie, isNmDb ;
  BOOL hitSupportAceView, flagAmbiguousProbe, flagAmbiguousProbeAceView ;
  BOOL removeNegativeHits, exportG2Nm2hits, exportP2Nm, flagRemoveNtHit ;
  BOOL exportAceViewMrnaAsGenomic, exportNmAsGenomic, exportMaqcMapping, exportMaqcMappingPerProbe ;
  BOOL showExonHisto ;
  BOOL showIntronHisto ;
  const char *template, *keysetName, *nm_file, *mrna_file ;
  Array results ;
  Array completes ;
  Array atgs ;
  Array stops ;
  Array introns ;
  Array intronCount ;
  Array exons ;
  Array exonCount ;
  Array exonHisto ;
  Array intronHisto ;
  Array atgHisto ;
  Array stopHisto ;
  Array mrnaHisto ;
  Array cdsHisto ;
  Array exon1Histo [16] ;
  Array exonLHisto [16] ;
  Array intronLastHisto ;
  Array intron1Histo [16] ;
  Array intronLHisto [16] ;
  int mrna2intronCount ;
  int mrna2exonCount ;
  DICT *dict ;
  Array mg ;
} MM ;

/*************************************************************************************/

/* uitility returns the IntMap values
   if *mapp == NULL, it is filled
   else gene->IntMap must == *mapp is inforced
*/

static BOOL magIntMap (AC_OBJ gene, AC_OBJ *mapp, int *m1p, int *m2p, AC_HANDLE h)
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

/***************************************************************************************/

static int magOrder (const void *a, const void *b)
{
  const int xa = *(const int *)a, xb = *(const int *)b ;
  
  return xa - xb ;
}


/*************************************************************************************/

static int magOneGene (MM *mm, int nS, AC_TABLE gSp, AC_HANDLE h)
{
  int i, j, ir, jr, kr, mx = 0, x, dx, x0, y, y0 ;
  Array xx = arrayCreate (2*nS, int) ;
  Array xy = arrayCreate (2*nS, int) ;
  BitSet bb = bitSetCreate (200, h) ;

  /* record all exons extremities, min==1 and max */
  for (ir = jr = kr = 0 ; ir < nS ; ir++)
    if (strstr (ac_table_printable (gSp, ir, 2, "toto"), "tron") &&
	strstr (ac_table_printable (gSp, ir, 3, "toto"), "_ag") 
	)
      {
	array (xx, 2*jr, int) = ac_table_int (gSp, ir, 0, 0) ;
	x = array (xx, 2*jr + 1, int) = ac_table_int (gSp, ir, 1, 0) ;
	jr++ ;
	if (x > mx) mx = x ; 
      } 
    else if (strstr (ac_table_printable (gSp, ir, 2, "toto"), "xon")
	)
      {
	array (xy, 2*kr, int) = ac_table_int (gSp, ir, 0, 0) ;
	x = array (xy, 2*kr + 1, int) = ac_table_int (gSp, ir, 1, 0) ;
	kr++ ;
	if (x > mx) mx = x ; 
      } 
  /* sort */
  arraySort (xx, magOrder) ;
  arrayCompress (xx) ;
  arraySort (xy, magOrder) ;
  arrayCompress (xy) ;
  /* flag the exon regions */
  bitUnSet (bb, nS) ; /* make room */
  for (ir = 0 ; ir < nS ; ir++)
    if (strstr (ac_table_printable (gSp, ir, 2, "toto"), "xon"))
      {
	x = ac_table_int (gSp, ir, 0, 0) ;
	for (i = 0 ; i < arrayMax (xy) ; i++)
	  if (x == arr (xy, i, int))
	    break ;
	x = ac_table_int (gSp, ir, 1, 0) ;
	for ( i++ ; i < arrayMax (xy) ; i++)
	  {
	    bitSet (bb, i) ;
	    if (x == arr (xy, i, int))
	      break ;
	  }
      }
  /* cumulate the exons regions */
  for (i = dx = 0, x0 = 1 ; i < arrayMax (xy) ; i++)
    {
      x = arr (xy, i, int) ;
      if (! bit (bb, i))
	dx += x - x0 - 1 ;
      x0 = x ;
    }
  mx -= dx ;  /* goto spliced coordinates */
  /* remove the introns regions */
  for (i = 0 ; i < arrayMax (xx) ; i++)
    {
      x = arr (xx, i, int) ; dx = 0 ;
      for (j = 0, y0 = 1 ; j < arrayMax (xy) ; j++)
	{
	  y = arr (xy, j, int) ;
	  if (y > x) y = x ;
	  if (bit (bb, j))
	    { dx += y - y0 + 1 ; }
	  y0 = y ;
	  if (y >= x)
	    break ;
	}
      arr (xx, i, int) = dx ;  /* goto spliced coordinates */
    }
  for (i = 0 ; i < arrayMax (xx) ; i++)
    {
      x = arr (xx, i, int) ;
      x = (100 * x + 1)/mx ;  /* rescale */
      array (mm->results, x, int) += 1 ;
    }
  return 1 ;
}  /* magOneGene */

/*************************************************************************************/

static int magOneMrna (MM *mm, BOOL isComplete, BOOL isCOOH_Complete,
		       int nS, AC_TABLE gSp, int p1, int p2, int *mxp, AC_HANDLE h)
{
  int i, ir, j, jr, mx = 0, a1, a2, x1, x2, x ;
  Array xx = arrayCreate (2*nS, int) ;

  /* record all introns positions */
  for (ir = jr = 0 ; ir < nS ; ir++)
    if (strstr (ac_table_printable (gSp, ir, 4, "toto"), "tron") &&
	strstr (ac_table_printable (gSp, ir, 5, "toto"), "_ag") 
	)
      {
	x1 = ac_table_int (gSp, ir, 2, 0) ;
	x2 = ac_table_int (gSp, ir, 3, 0) ;

	if (mm->isCDS && (x1 < p1 || x2 > p2))
	  continue ;

	array (xx, 2*jr, int) = x1 ;
	x = array (xx, 2*jr + 1, int) = x2 ; 
	jr++ ;
	if (x > mx) mx = x ; 
      } 
    else if (strstr (ac_table_printable (gSp, ir, 4, "toto"), "xon")
	)
      {
	x = ac_table_int (gSp, ir, 3, 0) ;
	if (x > mx) mx = x ; 
      } 
  /* sort */
  arraySort (xx, magOrder) ;
  arrayCompress (xx) ;

  if (mm->isCDS) mx = p2 - p1 + 1 ;
  for (i = 0 ; i < arrayMax (xx) ; i++)
    {
      x = arr (xx, i, int) ;
      if (mm->isCDS) x = x - p1 + 1 ;
      x = (100 * x + 1)/mx ;  /* rescale */
      array (isComplete ? mm->completes : mm->results, x, int) += 1 ;
    }
  *mxp = mx ;
  /* record all introns lengths */
  {
    mm->mrna2intronCount++ ;
    for (ir = jr = 0 ; ! isComplete && ir < nS ; ir++)
      if (strstr (ac_table_printable (gSp, ir, 4, "toto"), "tron") &&
	  strstr (ac_table_printable (gSp, ir, 5, "toto"), "_ag") 
	  )
	{
	  a1 = ac_table_int (gSp, ir, 0, 0) ;
	  a2 = ac_table_int (gSp, ir, 1, 0) ;
	  x1 = ac_table_int (gSp, ir, 2, 0) ;
	  x2 = ac_table_int (gSp, ir, 3, 0) ;
	  
	  if (mm->isCDS && (x1 < p1 || x2 > p2))
	    continue ;
	  
	  x = (100 * x1 + 1)/mx ;  /* rescale */
	  if (x < 0) x = 0 ;
	  array (mm->introns, x, int) += (a2 - a1 + 1) ;
	  array (mm->intronCount, x, int) += 1 ;
	  
	  /* record exon histo */
	  array (mm->intronHisto, a2 - a1 + 1, int)++ ;
	} 
  }
  if (isComplete && p1 >= 0 && mx - p2 >= 0) /* record atg/stp histos */
    {
      array (mm->atgHisto, p1, int)++ ;
      array (mm->stopHisto, mx - p2, int)++ ;
      array (mm->cdsHisto, p2 - p1 + 1, int)++ ;
    }
  if (! isComplete)
    array (mm->mrnaHisto, mx, int)++ ;
    
  /* record all exons lengths */
  if (isCOOH_Complete && !isComplete)
    {
      mm->mrna2exonCount++ ;
      for (ir = jr = 0 ; ir < nS ; ir++)
	if (strstr (ac_table_printable (gSp, ir, 4, "toto"), "xon")
	    )
	  {
	    a1 = ac_table_int (gSp, ir, 0, 0) ;
	    a2 = ac_table_int (gSp, ir, 1, 0) ;
	    x1 = ac_table_int (gSp, ir, 2, 0) ;
	    x2 = ac_table_int (gSp, ir, 3, 0) ;
	    
	    if (mm->isCDS && (x1 < p1 || x2 > p2))
	      continue ;
	    
	    x1 = (1000 * x1 + 1)/mx ;  /* rescale */
	    if (x1 < 0) x1 = 0 ;
	    x2 = (1000 * x2 + 1)/mx ;  /* rescale */
	    if (x2 < 0) x2 = 0 ;
	    /* count touched sections */
	    for (i = j = 0 ; i < 100 ; i++) 
	      if (x1 < 10 * (i+1) && x2 > 10 *i) j++ ;
	    /* record touched sections */
	    for (i = 0 ; j > 0 && i < 100 ; i++) 
	      if (x1 < 10 * (i+1) && x2 > 10 *i) 
		{
		  array (mm->exons, i, int) += (a2 - a1 + 1) * 100 / j ;
		  array (mm->exonCount, i, int) += 100/j ;
		} 
	    /* record exon histo */
	    array (mm->exonHisto, a2 - a1 + 1, int)++ ;
	  }
    }
  return 1 ;
}  /* magOneMrna */

/*************************************************************************************/

static int magGene (MM *mm, AC_OBJ gene)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE gProd ;
  AC_TABLE gSpl = ac_tag_table (gene, "Splicing", h) ;
  int nMag = 0, nS = gSpl ? gSpl->rows : 0, mx ;
  BOOL isCOOH_Complete ;

  if (nS > 1)
    {
      if (mm->isMrna)
	{
	  int x1 = 0, x2 = 0, ir ;
	  AC_OBJ prod ;

	  gProd =  ac_tag_table (gene, "Product", h) ;
	  for (ir = 0 ; gProd && ir < gProd->rows ; ir++)
	    {
	      prod = ac_table_obj (gProd, ir , 0, h) ;
	      isCOOH_Complete = ac_has_tag (prod, "COOH_COMPLETE") ;
		if (ac_has_tag (prod, "Best_product") && !mm->isCDS)
		  nMag = magOneMrna (mm, 0, isCOOH_Complete, nS, gSpl, x1, x2, &mx, h) ;

	      if (ac_has_tag (prod, "Best_product") &&
		  ac_has_tag (prod, "complete"))
		{
		  x1= ac_table_int (gProd, ir, 1, 0) ;
		  if (x1 < 0) x1 = 1 ; /* happens when we steal from */
		  x2 = ac_table_int (gProd, ir, 2, 0) ;

		  nMag = magOneMrna (mm, 1, isCOOH_Complete, nS, gSpl, x1, x2, &mx, h) ;

		  if (mm->isCDS)  { x2 = x2 - x1 + 1 ; x1 = 1 ; }
		  if (x2 > mx) x2 = mx ;
		  x1 = (100 * x1 + 1)/mx ;  /* rescale */
		  array (mm->atgs, x1, int)++ ;
		  x2 = (100 * x2 + 1)/mx ;  /* rescale */
		  array (mm->stops, x2, int)++ ;
		}	      
	    }	  
	}
      else
	nMag = magOneGene (mm, nS, gSpl, h) ;
    }

  ac_free (h) ;
  return nMag ;
} /* magGene */

/*************************************************************************************/

static void magOneIntronExonAbsolute (MM *mm, AC_OBJ mrna)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE spl = ac_tag_table (mrna, "Splicing", h) ;
  int x, ir, nIntrons = 0, nx, ni ;
    
  for (ir = 0 ; ir < spl->rows ; ir++)
    if (strstr (ac_table_printable (spl, ir, 4, ""), "tron"))
      nIntrons++ ;
  for (ir = ni = nx = 0 ; ir < spl->rows ; ir++)
    {
      if (strstr (ac_table_printable (spl, ir, 4, ""), "tron"))
	{
	  x = ac_table_int (spl, ir, 1, 0) - ac_table_int (spl, ir, 0, 0) + 1 ;
	  ni++ ;
	  if (nIntrons >= ni + 1 && ni <= 16)
	    array (mm->intron1Histo[ni-1], x, int) += 1 ;
	  if (nIntrons == ni && ni <= 16)
	    array (mm->intronLHisto[ni-1], x, int) += 1 ;
	  if (nIntrons == ni)
	    array (mm->intronLastHisto, x, int) += 1 ;
	}	
      if (strstr (ac_table_printable (spl, ir, 4, ""), "xon"))
	{
	  x = ac_table_int (spl, ir, 1, 0) - ac_table_int (spl, ir, 0, 0) + 1 ;
	  /* use ni */
	  if (nIntrons >= ni + 2 && ni < 16)
	    array (mm->exon1Histo [ni], x, int) += 1 ;
	  if (nIntrons - ni < 16 &&  ni >= 2)
	    array (mm->exonLHisto [nIntrons - ni], x, int) += 1 ;
	}	
    }
  ac_free (h) ;
  return ;
} /* magIntronExonAbsolute */

/*************************************************************************************/

static void magIntronExonAbsolute (MM *mm, AC_OBJ mrna0)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER iter = ac_objquery_iter (mrna0, "CLASS mrna ; COUNT {>product; best_product ; ! complete } > 0 && ! gap", h) ;
  AC_OBJ mrna = 0 ;
    
  while (ac_free (mrna), mrna = ac_iter_obj (iter))
    magOneIntronExonAbsolute (mm, mrna) ;
  ac_free (h) ;
  return ;
} /* magIntronExonAbsolute */

/*************************************************************************************/

static void histoStats (Array histo, int *totalp, int *avp, int *medianp, int *sigmap)
{
  int x, nn, tot = 0, cumul = 0, tot2 = 0, median = 0 ;
  double s2 = 0, s3 = 0 ;

  for (x = 0 ; x < arrayMax (histo) ; x++)
    {
      nn = arr (histo, x, int) ;
      cumul += x*nn ; tot += nn ;
      s2 += x*x*nn ;
    }
  for (x = 0 ; x < arrayMax (histo) ; x++)
    {
      nn = arr (histo, x, int) ;
      tot2 += nn ;
      if (2*tot2 > tot) { median = x ; break ; }
    }
  *totalp = tot ;
  if (!tot) tot = 1 ;
  *avp = cumul / tot ;
  *medianp = median ;
  s3 = cumul ; s3 = s3/tot ; s3 = s3 * s3 ;
  *sigmap = (int) (sqrt (s2/tot - s3)) ;
} /* histoStats */

/*************************************************************************************/

static void magShowResults (MM *mm)
{
  int i, nn, nn1, nx, nx1 ;
  if (mm->isMrna)
    printf ("%3s\t%5s\t%5s\t%8s\t%5s\t%3s\t%6s\t%8s\t%5s\n"
	    , "  #", "all", "compl", "atg", "stops"
	    , "  #", "int-ln", "#exon", "exon-ln") ;
  else
    printf ("#\tall\n") ;
  for (i = 0 ; i < arrayMax (mm->results) ; i++)
    if (mm->isMrna)
      {
	nn = nn1 = array (mm->intronCount, i, int) ;
	if (nn <= 1) nn1 = 1 ;
	nx = nx1 = array (mm->exonCount, i, int) ;
	if (nx <= 1) nx1 = 1 ;
	printf ("%3d\t%5d\t%5d\t%8d\t%5d\t"
		, i
		, array (mm->results, i, int)
		, array (mm->completes, i, int)
		, array (mm->atgs, i, int)
		, array (mm->stops, i, int)
		) ;
	printf ("%3d\t%5d\t%8d\t%5d\n"
		, i
		, array (mm->introns, i, int)/nn1
		, nx1/100
		, array (mm->exons, i, int)/nx1
		) ;
      }
    else
      printf ("%d\t%d\n"
	      , i
	      , array (mm->results, i, int)
	      ) ;
  

  printf ("\n") ;
  if (1)
    {
      int ii, total, av, median, sigma ;

      histoStats (mm->mrnaHisto, &total, &av, &median, &sigma) ;
      printf ("%8d\tmRNAs average length\t%5d bp, median\t%5d\tsigma\t%5d\n"
	      , total, av, median, sigma) ;
      histoStats (mm->cdsHisto, &total, &av, &median, &sigma) ;
      printf ("%8d\tCDSs average length\t%5d\tbp, median\t%5d\tsigma\t%5d\n"
	      , total, av, median, sigma) ;
      histoStats (mm->exonHisto, &total, &av, &median, &sigma) ;
      printf ("%8d\texons average length\t%5d\tbp, median\t%5d\tsigma\t%5d\n"
	      , total, av, median, sigma) ;
      
      for (ii = 0 ; ii < 16 ; ii++)
	{
	  histoStats (mm->exon1Histo[ii], &total, &av, &median, &sigma) ;
	  printf ("%8d\t%d th exon, average length\t%5d\tbp, median\t%5d\tsigma\t%5d\n"
	      , total, ii+1, av, median, sigma) ;
	}
      for (ii = 15 ; ii >= 0 ; ii--)
	{
	  histoStats (mm->exonLHisto[ii], &total, &av, &median, &sigma) ;
	  printf ("%8d\tlast - %d th exon, average length\t%5d\tbp, median\t%5d\tsigma\t%5d\n"
	      , total, ii, av, median, sigma) ;
	}

      for (ii = 0 ; ii < 16 ; ii++)
	{
	  histoStats (mm->intron1Histo[ii], &total, &av, &median, &sigma) ;
	  printf ("%8d\t%d th intron, average length\t%5d\tbp, median\t%5d\tsigma\t%5d\n"
	      , total, ii + 1, av, median, sigma) ;
	}
      for (ii = 0 ; ii < 16 ; ii++)
	{
	  histoStats (mm->intronLHisto[ii], &total, &av, &median, &sigma) ;
	  printf ("%8d\tlast is %d th intron, average length\t%5d\tbp, median\t%5d\tsigma\t%5d\n"
	      , total, ii+1, av, median, sigma) ;
	}
      histoStats (mm->intronLastHisto, &total, &av, &median, &sigma) ;
      printf ("%8d\tlast is any intron, average length\t%5d\tbp, median\t%5d\tsigma\t%5d\n"
	      , total, av, median, sigma) ;

    }

  if (1)
    {
      int total, av, median, sigma ;

      histoStats (mm->atgs, &total, &av, &median, &sigma) ;
      printf ("Found %8d MET, average position %d%%, median %d%%\n"
	      , total, av, median) ;
    }
  if (1)
    {
      int total, av, median, sigma ;

      histoStats (mm->stops, &total, &av, &median, &sigma) ;
      printf ("Found %8d stops, average position %d%%, median %d%%\n"
	      , total, av, median) ;
    }
  if (0 && mm->showExonHisto)
    {
      printf ("\n\nExonHisto limited at 300\n") ;  
      printf ("%3s%5s%5s%5s%5s%5s%5s%5s\n","bp","all","1","2","3","-3","-2","-1") ;
      for (i = 0 ; i <= 300 && i < arrayMax (mm->exonHisto) ; i++)
	{
	  nn = arr (mm->exon1Histo[0], i, int) ; printf ("%5d", nn) ;
	  printf ("\n") ;
	}
    }
} /* magShowResults */

/*************************************************************************************/

static void magAnalyse (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER iter ;
  AC_OBJ gene = 0 ;
  AC_KEYSET ks = 0 ;
  int ii, nn = 0 ;
  char *doMrna = mm->isMrna ? " ; >mrna  gt_ag > 0 && ! gap " : "" ;

  mm->atgs = arrayHandleCreate (101, int, h) ;
  mm->stops = arrayHandleCreate (101, int, h) ;
  mm->results = arrayHandleCreate (101, int, h) ;
  mm->completes = arrayHandleCreate (101, int, h) ;
  mm->introns = arrayHandleCreate (101, int, h) ;
  mm->intronCount = arrayHandleCreate (101, int, h) ;
  mm->exons = arrayHandleCreate (101, int, h) ;
  mm->exonCount = arrayHandleCreate (101, int, h) ;
  mm->exonHisto = arrayHandleCreate (10000, int, h) ;
  mm->intronHisto = arrayHandleCreate (10000, int, h) ;
  mm->atgHisto = arrayHandleCreate (10000, int, h) ;
  mm->stopHisto = arrayHandleCreate (10000, int, h) ;
  mm->cdsHisto = arrayHandleCreate (10000, int, h) ;
  mm->mrnaHisto = arrayHandleCreate (10000, int, h) ;
  mm->intronLastHisto = arrayHandleCreate (10000, int, h) ;

  for (ii = 0 ; ii < 16 ; ii++)
    {
      mm->exon1Histo [ii] = arrayHandleCreate (10000, int, h) ;
      mm->exonLHisto [ii] = arrayHandleCreate (10000, int, h) ;
      mm->intron1Histo [ii] = arrayHandleCreate (10000, int, h) ;
      mm->intronLHisto [ii] = arrayHandleCreate (10000, int, h) ;
    }

  array (mm->results, 100, int) = 0 ; /* make room */
  if (mm->keysetName)
    {
      ks = ac_read_keyset (mm->db, mm->keysetName, h) ;
      if (!ks || !ac_keyset_count (ks))
	{ 
	  fprintf(stderr, "cannot find keyset %s", mm->keysetName) ;
	  ac_db_close (mm->db) ;
	  usage () ;
	}
      iter = ac_keyset_iter (ks, TRUE, h) ;
    }
  else if (mm->template)
    {
      iter = ac_dbquery_iter (mm->db 
			      , messprintf ("Find transcribed_gene %s ; gt_ag %s "
					    , mm->template, doMrna)			    
			      , h) ;
    }
  else
    {
       iter = ac_dbquery_iter (mm->db 
			      , messprintf ("Find transcribed_gene gt_ag %s"
					    , doMrna)
			      , h) ;
    }
  while (ac_free (gene), gene = ac_next_obj (iter))
    {
      nn += magGene (mm, gene) ;
      magIntronExonAbsolute (mm, gene) ;
    }
  
  arr (mm->results, 100, int) = 0 ;
  arr (mm->completes, 100, int) = 0 ;
  magShowResults (mm) ;
  printf ("// Analysed %d %s\n\n", nn, ! mm->isMrna ? "genes" : "transcripts") ;
  fprintf (stderr, "// Analysed %d %s\n\n", nn, ! mm->isMrna ? "genes" : "transcripts") ;

  ac_free (h) ;
  return ;
} /* magAnalyse */

/*************************************************************************************/
/*************************************************************************************/
/* export the probe mapped on the mrna using 
 * acem
 *   align_probe -mrna
 */
typedef struct mgStruct { int probe, onNm, nHit, nExon, nClonesInMrna, nClonesInGene, chrom, mrna, gene, a1, a2, x1, x2 ; } MG ;

static int mgOrder (const void *a,const   void *b)
{
  const MG *ma = (const MG *)a, *mb = (const MG *)b ;
  int nn ;

  nn = ma->probe - mb->probe ; if (nn) return nn ;
  nn = ma->nHit - mb->nHit ; if (nn) return nn ;
  nn = ma->chrom - mb->chrom ; if (nn) return nn ;
  nn = ma->nExon - mb->nExon ; if (nn) return nn ;
  nn = ma->x1 - mb->x1 ; if (nn) return nn ;
  nn = ma->x2 - mb->x2 ; if (nn) return nn ;
  nn = ma->a1 - mb->a1 ; if (nn) return nn ;
  return 0 ;
} /* mgOrder */

/*************************************************************************************/

static void magMapOneProbe (MM *mm, AC_OBJ probe)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ gMrna = 0, mrna = 0, chrom = 0, gene ;
  AC_TABLE spl = 0, hits = ac_tag_table (probe, "nm_hit", h), mrnas ;
  int ir, jr, v1, v2, u0, u1, u2, a1, a2, x1, x2, y1, y2, b1, b2, m1, m2, da, xm1, xm2 ;
  int mPrim = 0, mChrom = 0, mMrna = 0, mGene = 0, nHit, nExon = 0 ;
  int ii1, ii2 ;
  int nClonesInMrna = 0, nClonesInGene = 0 ;
  MG *mg = 0 ;
  BOOL ok ;
  char * signature[10000] ;

  ii1 = arrayMax (mm->mg) ;

  for (ir = nHit = 0 ; hits && ir < hits->rows && nHit < 4000 ; ir++)
    {
      nClonesInMrna = nClonesInGene = nExon = 0 ;
      gMrna = ac_table_obj (hits, ir, 0, h) ;
      u0 = u1 = ac_table_int (hits, ir, 1, -9999) ;
      u2 = ac_table_int (hits, ir, 2, -9999) ;
      if (! gMrna || u1 < 0 || u2 < 0 || u1 >= u2)
	continue ;
      if (! magIntMap (gMrna, &chrom, &m1, &m2, h))
	continue ;

      v1 = 1 ; v2 = u2 - u1 + 1 ;
      if (ac_has_tag (gMrna, "Exon_coord")) /* we are on an NM */
	{
	  dictAdd (mm->dict, ac_name (probe), &mPrim) ;
	  dictAdd (mm->dict, ac_name (gMrna), &mMrna) ;
	  dictAdd (mm->dict, ac_name (chrom), &mChrom) ;
	  gene = ac_tag_obj (gMrna, "Model_of_gene", h) ;
	  if (gene)
	    dictAdd (mm->dict, ac_name(gene), &mGene) ;

	  spl = ac_tag_table (gMrna, "Exon_coord", h) ;
	  
	  ii2 = arrayMax (mm->mg) ;

	  for (jr = 0 ; spl && jr < spl->rows ; jr++)
	    {
	      a1 = ac_table_int (spl, jr, 0, 0) ;
	      a2 = ac_table_int (spl, jr, 1, 0) ;
	      x1 = ac_table_int (spl, jr, 2, 0) ;
	      x2 = ac_table_int (spl, jr, 3, 0) ;
 	
	      y1 = u1 < x1 ? x1 : u1 ;
              y2 = u2 < x2 ? u2 : x2 ;
	      if (y1 <= y2)
		{
		  nHit++ ; 
		  /* register */
		  nExon++ ;
		  mg = arrayp (mm->mg, arrayMax (mm->mg), MG) ;
		  mg->onNm = TRUE ;
		  mg->probe = mPrim ;
		  mg->chrom = mChrom ;
		  mg->mrna = mMrna ;
		  mg->gene = mGene ;
		  mg->nClonesInGene = ac_tag_int (gene, "Nb_cDNA_clone", 0) ;
		  mg->nClonesInMrna = ac_tag_int (gMrna, "Nb_cDNA_clone", 0) ;
		  b1 = a1 + y1 - x1 ;
		  b2 = a1 + y2 - x1 ;
		  if (m1 < m2)
		    {
		      mg->a1 = m1 + b1 - 1 ;
		      mg->a2 = m1 + b2 - 1 ;
		    }
		  else
		    {
		      mg->a1 = m1 - b1 + 1 ;
		      mg->a2 = m1 - b2 + 1 ;
		    }
		  mg->nHit = nHit ; 
		  mg->nExon = nExon = jr + 1 ;
		  mg->x1 = y1 ; mg->x2 = y2 ;
		}
	    }
	}
      else
	{
	  mrnas = ac_tag_table (gMrna, "mRNAs", h) ;
	  mrna = ac_table_obj (mrnas, 0, 0, h) ;
	  xm1 =  ac_table_int (mrnas, 0, 1, 0) ;
	  xm2 =  ac_table_int (mrnas, 0, 2, 0) ;
	  if (!mrna)
	    continue ;
	  ac_free (chrom) ;
	  if (! magIntMap (mrna, &chrom, &m1, &m2, h))
	    continue ;
	  if (u1 < 1 || u2 > xm2)
	    continue ;
	  gene = ac_tag_obj (mrna, "From_gene", h) ;
	  nClonesInGene = ac_keyset_count (ac_objquery_keyset (gene, ">mrna;>cdna_clone",h)) ;
	  nClonesInMrna = ac_keyset_count (ac_objquery_keyset (mrna, ">cdna_clone",h)) ;
	  dictAdd (mm->dict, ac_name (probe), &mPrim) ;
	  dictAdd (mm->dict, ac_name (chrom), &mChrom) ;
	  dictAdd (mm->dict, ac_name (mrna) + 3, &mMrna) ;
	  dictAdd (mm->dict, ac_name(gene), &mGene) ;
	  if (strstr (ac_name (chrom), "NT"))
	    continue ;
	  spl = ac_tag_table (mrna, "Splicing", h) ;
	  
	  nHit++ ; 
	  ii2 = arrayMax (mm->mg) ;
	}
      /* register */
      for (jr = nExon = da = 0, ok = FALSE ; spl && !ok && jr < spl->rows ; jr++)
	{
	  if (strstr (ac_table_printable (spl, jr, 4, "toto"), "tron")) /* intron */
	    {
	      da += ac_table_int (spl, jr, 7, -9999) ;
	      if (da < 0) ok = TRUE ;
	      continue ;
	    }
	  if (!strstr (ac_table_printable (spl, jr, 4, "toto"), "xon")) /* exon ? */
	    { ok = TRUE ; continue ; } /* discard Gap and other unknown stuff */
	  a1 = da + ac_table_int (spl, jr, 0, -9999) ;   /* exon ok */
	  a2 = da + ac_table_int (spl, jr, 1, -9999) ;
	  x1 = ac_table_int (spl, jr, 2, -9999) ;
	  x2 = ac_table_int (spl, jr, 3, -9999) ;
	  if (x1 < 0 || x2 < 0)
	    continue ;
	  if (x1 <= u1 && x2 >= u1)
	    {
	      nExon++ ;
	      b1 = a1 + u1 - x1 ;
	      if (x1 <= u2 && x2 >= u2)
		{ b2 = a1 + u2 - x1 ; ok = 1 ; v2 = u2 ; }
	      else
		{ b2 = a2 ; v2 = x2 ; }
	      
	      mg = arrayp (mm->mg, arrayMax (mm->mg), MG) ;
	      mg->probe = mPrim ;
	      mg->chrom = mChrom ;
	      mg->mrna = mMrna ;
	      mg->gene = mGene ;
	      mg->nClonesInGene = nClonesInGene ;
	      mg->nClonesInMrna = nClonesInMrna ;
	      if (m1 < m2)
		{
		  mg->a1 = m1 + b1 - 1 ;
		  mg->a2 = m1 + b2 - 1 ;
		}
	      else
		{
		  mg->a1 = m1 - b1 + 1 ;
		  mg->a2 = m1 - b2 + 1 ;
		}
	      mg->nHit = nHit ; 
	      mg->nExon = nExon ;
	      mg->x1 = u1 - u0 + 1 ; mg->x2 = v2 - u0 + 1 ;
	      
	      if (!ok)
		u1 = x2 + 1 ;
	    }
	}
      if (mg && nExon > 0 && nHit > 0)
	{
	  int jj ;
	  BOOL isNew ;
	  
	  /* construct a signature string with all coordinates */
	  signature [nHit - 1] = halloc (100 + 40*nExon, h) ;
	  sprintf (signature [nHit-1], "%d/", mg->chrom) ;
	  for (jj = ii2 ; jj < arrayMax (mm->mg) ; jj++)
	    strcat (signature [nHit - 1], messprintf ("%d/%d/%d/%d", mg->a1, mg->a1, mg->x1, mg->x2)) ;
	  
	  /* compare to previous */
	  for (jj = 0, isNew = TRUE ; isNew && jj < nHit - 1 ; jj++)
	    if (signature [jj] && !strcmp (signature [jj], signature [nHit-1]))
	      {
		int j ;
		MG *mg2 ;
		isNew = FALSE ;
		for (j = 0, mg2 = arrp (mm->mg, 0, MG) ; j < arrayMax (mm->mg) ; mg2++, j++)
		  if (mg2->nHit == jj+1) mg2->nClonesInMrna += nClonesInMrna ;
	      }
		  
	  if (!isNew)   /* if not new, destroy */
	    {
	      nHit-- ;
	      arrayMax (mm->mg) = ii2 ;
	    }
	}
    }
  ac_free (h) ;
} /* magMapOneProbe */

/*************************************************************************************/
/*   Export the table */
static void magMapProbeExportTable (MM *mm)
{
  int ii ;
  MG *mg ;
  DICT *dict = mm->dict ;

  freeOutf  ("#probe\t#Hit\t#Exon\tchrom\ta1\ta2\tx1\tx2\tgene\t#clones\tmrna\t#clones (base 1, not zero)\n") ;
  for (ii = 0, mg = arrp (mm->mg, 0, MG) ; ii < arrayMax (mm->mg) ; ii++, mg++)
    {
      freeOutf ("%s\t%d\t%d\t"
		, dictName(dict, mg->probe), mg->nHit, mg->nExon) ;
      freeOutf ("%s\t%d\t%d\t%d\t%d\t"
		, dictName(dict, mg->chrom)
		, mg->a1, mg->a2, mg->x1, mg->x2
		/* ,  dictName(dict, mg->mrna) */
		) ;
      freeOutf ("%s\t%d\t%s\t%d\n"
		, dictName(dict, mg->gene)
		, mg->nClonesInGene
		, dictName(dict, mg->mrna)
		, mg->nClonesInMrna
		/* ,  dictName(dict, mg->mrna) */
		) ;
    }
  return ;
} /* magMapProbeExportTable */

/*************************************************************************************/

static void magMapProbeShow (MM *mm)
{

  if (!mm->mg || ! arrayMax (mm->mg))
    { freeOutf ("# Empty map\n") ; return ; }
  
  arraySort (mm->mg, mgOrder) ;
  arrayCompress (mm->mg) ;
  
  magMapProbeExportTable (mm) ;

  return ;
} /* magMapProbeShow */

/*************************************************************************************/

static void magMapProbe (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER probes ;
  AC_OBJ probe = 0 ;

  mm->dict = dictHandleCreate (100000, h) ; dictAdd (mm->dict, "avoid zero", 0) ;
  mm->mg = arrayHandleCreate (100000, MG, h) ;

  probes = ac_dbquery_iter (mm->db
			     , messprintf ("Find Probe %s ; Nm_Hit", mm->template ? mm->template : "")
			     , h) ;

  while (ac_free (probe), probe = ac_iter_obj (probes))
    magMapOneProbe (mm, probe) ; 

  magMapProbeShow (mm) ;

  ac_free (h) ;
  return ;
} /* magMapProbe */

/*************************************************************************************/

static void magExportAceViewMrnaAsGenomic (MM *mm)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_ITER genes ;
  AC_ITER mrnas ;
  AC_TABLE table ;
  int ir, nn ;
  AC_KEYSET ks ;
  AC_OBJ gene = 0, mrna = 0 ;
  char *dna ;
  int nnngenes = 0, nnnmrnas = 0 ;

  genes = ac_dbquery_iter (mm->db, "find gene", h) ;
  while (ac_free (gene), gene = ac_iter_obj (genes))
    {
      printf ("\nGene \"%s\"\n", ac_name (gene)) ;
      if (ac_has_tag (gene, "Balise"))
	printf ("Main\n") ;
      else if (ac_has_tag (gene, "Cloud_gene"))
	printf ("Cloud\n") ;
      else
	printf ("Putative\n") ;
      ks = ac_objquery_keyset (gene, ">transcribed_gene ; > cdna_clone", h) ;
      nn = ac_keyset_count (ks) ;
      ac_free (ks) ;
      if (nn)
	printf ("Nb_cDNA_clone %d\n", nn) ;
      table = ac_tag_table (gene, "GeneId", h) ;
      for (ir = 0; table && ir < table->rows ; ir++)
	printf ("GeneId \"%s\"\n", ac_table_printable (table, ir, 0, "toto")) ;
      ac_free (table) ;
      mrnas = ac_objquery_iter (gene, ">transcribed_gene;>mrna", h) ;
      while (ac_free (mrna), mrna = ac_iter_obj (mrnas))
	{
	  h1 = ac_new_handle () ;

	  printf ("\nSequence \"GM_%s\"\n", ac_name(mrna)) ;
	  printf ("Genomic\n") ;
	  printf ("Model_of_gene \"%s\"\n", ac_name(gene)) ;
	  table = ac_tag_table (mrna, "IntMap", h1) ;
	  if (table && table->rows && table->cols >= 3)
	    printf ("IntMap \"%s\" %d %d\n"
		    , ac_table_printable (table, 0, 0, "toto")
		    , ac_table_int (table, 0, 1, 0)
		    , ac_table_int (table, 0, 2, 0)
		    ) ;

	  table = ac_tag_table (mrna, "Splicing", h1) ;
	  for (ir = 0; table && ir < table->rows ; ir++)
	    if (strstr (ac_table_printable (table, ir, 4, "toto"), "xon"))
	      printf ("Exon_coord %d %d %d %d \n"
		    , ac_table_int (table, ir, 0, 0)
		    , ac_table_int (table, ir, 1, 0)
		    , ac_table_int (table, ir, 2, 0)
		    , ac_table_int (table, ir, 3, 0)
		    ) ;

	  ks = ac_objquery_keyset (mrna, ">cdna_clone", h1) ;
	  nn = ac_keyset_count (ks) ;

	  if (nn)
	    printf ("Nb_cDNA_clone %d\n", nn) ;
	  dna = ac_obj_dna (mrna, h1) ;
	  printf ("\nDNA  \"GM_%s\"\n%s\n\n", ac_name(mrna), dna) ;
	  nnnmrnas++ ;
	  ac_free (h1) ;
	}
      ac_free (mrnas) ;
      nnngenes++ ;
    }
  fprintf (stderr, "\n// exported %d genes %d mrnas\n", nnngenes, nnnmrnas) ;
  ac_free (h) ;
  return ;
} /* magExportAceViewMrnaAsGenomic */

/*************************************************************************************/

static void magExportNmAsGenomic (MM *mm)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_ITER nms ;
  AC_ITER genes ;
  AC_TABLE table ;
  int ir ;
  AC_OBJ nm = 0, gene = 0 ;
  char *dna ;
  int nnngenes = 0, nnnNms = 0 ;

  nms = ac_dbquery_iter (mm->db, "find EST Ref_seq", h) ;
  while (ac_free (nm), nm = ac_iter_obj (nms))
    {
      h1 = ac_new_handle () ;

      printf ("Sequence \"%s\"\n", ac_name (nm)) ;
      printf ("Genomic\n") ; 
      table = ac_tag_table (nm, "GeneId", h1) ;
      for (ir = 0; table && ir < table->rows ; ir++)
	printf ("GeneId \"%s\"\n", ac_table_printable (table, ir, 0, "toto")) ;
 
      genes = ac_objquery_iter (nm, ">from_gene ; > Gene", h1) ;
      while (ac_free (gene), gene = ac_iter_obj (genes))
	{
	  printf ("Model_of_gene \"%s\"\n", ac_name(gene)) ;
	  nnngenes++ ;
	}

      dna = ac_obj_dna (nm, h1) ;
      printf ("\nDNA  \"%s\"\n%s\n\n", ac_name(nm), dna) ;
      nnnNms++ ;
      ac_free (h1) ;
    }
  fprintf (stderr, "// exported %d NMs from %d genes\n", nnnNms, nnngenes) ;
  ac_free (h) ;
  return ;
} /* magExportNmAsGenomic */

/*************************************************************************************/
/*************************************************************************************/
/* herman and aamggie projetcs
 *  count how many genes are hit by each type of probe etc
 */
typedef struct hitStatStruct { int nam, iType, p, p_nm, p_1_nm, p_2_nm, p_3_nm, p_4_nm, p_5_nm, p_6_nm, p_7_nm, p_8_nm, p_9_nm, p_10_nm, nm, nm_p, p_1g, p_1gc, p_1ga, g_p1, gid_p1 ; AC_KEYSET nm_ps, g_p1s ; } HSS ;
typedef struct libStatStruct { KEY nm, gene, locus ; int nLib, nLib2 ; int libs[12] ;} HLS ;

static int hlsOrder (const void *a, const void *b)
{
  const HLS *ha = (const HLS *)a, *hb = (const HLS *)b ;
  int n ;
  
  n = ha->nLib - hb->nLib ;
  if (n) return -n ; /* many lib first */

  n = ha->nLib2 - hb->nLib2 ;
  if (n) return n ; /* small lib first */
  
  if (ha->gene != hb->gene) return (ha->gene < hb->gene) ? -1 : 1 ;
  if (ha->nm != hb->nm) return (ha->nm < hb->nm) ? -1 : 1 ;
  return 0 ;
} /* hlsOrder */

static void hlsReport (const char **types, Array hls, char *title)
{
  int nn, n, i, ii, total = 0, ngeneid = -1 ;
  HLS *hh ;
  
  hh = arrayp (hls, 0, HLS) ;
  nn = hh->nLib ;
  for (n = ii = 0, hh = arrayp (hls, 0, HLS) ; ii < arrayMax(hls) ; hh++, ii++)
    if (hh->nLib == nn) n++ ;
    else
      {
	total += n ;
	if (n)
	  printf ("%2d\tplatform test\t%d\t%ss\t%d\t%2.2f%%\t%d\n"
		  , nn, n, title
		  , total, (100.0 * total)/arrayMax (hls)
		  , ngeneid
		  ) ;
	nn = hh->nLib ;
	n = 1 ;
      }
  total += n ;
  printf ("%2d\tplatform test\t%d\t%ss\t%d\t%2.2f%%\n", nn, n, title, total, (100.0 * total)/arrayMax (hls)) ;
  printf ("\ttotal\t%8d\t%ss\n", total, title) ;
  printf ("analysed\t%8d\tlines\n", arrayMax (hls)) ; 

  {
    FILE *fil = filopen (messprintf ("%s_hits_per_number_of_platform", title), "txt", "w") ;
    int n6 ;

    if (fil)
      {
	fprintf (fil, "# %s\tLocus\t%s\t%s\t%s", title, "#platform", "in6", types[0]) ;
	for (i = 1 ; types[i] ; i++) fprintf (fil, "\t%s", types[i]) ;
	fprintf (fil, "\n") ;
	for (n = ii = 0, hh = arrayp (hls, 0, HLS) ; ii < arrayMax(hls) ; hh++, ii++)
	  {
	    for (n6 = 1, i = 0 ; i < 6 ; i++) 
	      if (! hh->libs[i])
		n6 = 0 ;
	    fprintf (fil, "%s\t%s\t%8d\t%3d"
		     , (char *) ac_key_name(hh->nm)
		     , (char *) ac_key_name(hh->locus)
		     , hh->nLib, n6) ;
	    for (i = 0 ; types[i] ; i++) fprintf (fil, "\t%3d", hh->libs[i]) ;
	    fprintf (fil, "\n") ;
	  }
	filclose (fil) ;
      }
  }
}

static void hlsStats (MM *mm, BOOL isAceView)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_KEYSET ks ;
  AC_ITER iter ;
  AC_OBJ nm = 0 ;
  HLS *hh ;
  int ii, n, nn = 0 ;
  Array hls = arrayHandleCreate (25000, HLS, h) ;
  Array hlg = arrayHandleCreate (25000, HLS, h) ;
  /*
    const char *types[] = { "ABI", "AFX", "AGL", "GEH", "ILM",  "NCI", "GEX", "QGN", "TAQ", "EPP", 0} ;
    BOOL hasCluster[] =   { FALSE,  TRUE, FALSE, FALSE, FALSE, FALSE,  FALSE, FALSE,  FALSE, FALSE, FALSE } ;
  */
  const char *types[] = { "ABI", "AFX", "AGL", "GEH", "ILM",  "NCI", "MLT", 0} ;
  BOOL hasCluster[] =   { FALSE,  TRUE, FALSE, FALSE, FALSE, FALSE,FALSE,FALSE } ;

  if (mm->isMaggie) 
    {
      types[0] = "Affy_" ;
      types[1] = "GE_" ;
      types[2] = "Agilent_" ;
      types[3] = 0 ;
      hasCluster[0] = FALSE ;
      hasCluster[1] = FALSE ;
      hasCluster[2] = TRUE ;
      hasCluster[3] = FALSE ;
    }

  nn = 0 ; iter = ac_dbquery_iter (mm->db, messprintf ("Find sequence genomic && %s"
						       , mm->template ? mm->template : "genomic")
				   , h) ;
  while (ac_free (nm), nm = ac_iter_obj (iter))
    {
      hh = arrayp (hls, nn++, HLS) ;
      hh->nm = ac_obj_key (nm) ;
      if (isAceView)
	hh->locus = ac_tag_key (nm, "Model_of_gene", 0) ;
      else
	hh->locus = ac_tag_key (nm, "Locus_GB", 0) ;
      if (ac_has_tag (nm, "probe_nm_hit"))
	for (ii = 0; types[ii] ; ii++)
	  {
	    if (isAceView)
	      ks = ac_objquery_keyset (nm, messprintf (">probe_nm_hit; IS \"%s%s*\" %s ; nm_hit ; COUNT {>Probe_set ; non_ambiguous} > 0 || COUNT {>nm_hit; >Model_of_gene} = 1"
						       , types [ii]
						       , mm->isMaggie ? "" : "_"
						       , hasCluster[ii] ? "&& ! IS *_at*" : ""), 0) ;
	    else
	      ks = ac_objquery_keyset (nm, messprintf (">probe_nm_hit; IS \"%s%s*\" %s ; nm_hit ; COUNT {>nm_hit; >GeneId} = 1"
						       , types [ii]
						       , mm->isMaggie ? "" : "_"
						       , hasCluster[ii] ? "&& ! IS *_at*" : ""), 0) ;
	    hh->libs [ii] = n = ac_keyset_count (ks) ;
	    if (n) {  hh->nLib++ ; hh->nLib2 += ii ; }
	    ac_free (ks) ;
	  }
    }

  nn = 0 ; 
  if (isAceView)
    iter =  ac_dbquery_iter (mm->db, messprintf ("Find sequence genomic && %s ; >Model_of_gene"
						       , mm->template ? mm->template : "genomic")
			     , h) ;
  else
    iter = ac_dbquery_iter (mm->db, "Find geneid", h) ;
  while (ac_free (nm), nm = ac_iter_obj (iter))
    {
      hh = arrayp (hlg, nn++, HLS) ;
      hh->nm = ac_obj_key (nm) ; 
      if (isAceView)
	hh->locus = ac_obj_key (nm) ;
      else
	{ 
	  AC_OBJ nm2 = ac_tag_obj (nm, "Sequence", 0) ; 
	  hh->locus = ac_tag_key (nm2, "Locus_GB", 0) ;
	  ac_free (nm2) ;
	}
      for (ii = 0; types[ii] ; ii++)
	{
	  if (isAceView)
	    ks = ac_objquery_keyset (nm, messprintf (">Genefinder; >probe_nm_hit; IS \"%s%s*\" %s ; nm_hit ; COUNT {>Probe_set ; non_ambiguous} > 0 || COUNT {>nm_hit; >Model_of_gene} = 1"
						     , types [ii]
						     , mm->isMaggie ? "" : "_"
						     , hasCluster[ii] ? "&& ! IS *_at_*" : ""), 0) ;
	  else
	    ks = ac_objquery_keyset (nm, messprintf (">Sequence; >probe_nm_hit; IS \"%s%s*\" %s ; nm_hit ; COUNT {>nm_hit; >GeneId} = 1"
						     , types [ii]
						     , mm->isMaggie ? "" : "_"
						     , hasCluster[ii] ? "&& ! IS *_at_*" : ""), 0) ;
	  hh->libs [ii] = n = ac_keyset_count (ks) ;
	  if (n) {  hh->nLib++ ; hh->nLib2 += ii ; }
	  ac_free (ks) ;
	}
    }
  
  arraySort (hls, hlsOrder) ;
  arraySort (hlg, hlsOrder) ;
  if (isAceView)
    {
      hlsReport (types, hls, "Transcript") ;
      hlsReport (types, hlg, "Gene") ;
    }
  else
    {
      hlsReport (types, hls, "NM") ;
      hlsReport (types, hlg, "GeneId") ;
    }
  
  ac_free (h) ;
}

static void magOneHitStats (MM *mm, HSS *hp, const char *type, BOOL hasCluster, BOOL isNmDb, AC_HANDLE h)
{ 
  AC_HANDLE h1 = ac_new_handle () ;
  AC_KEYSET probes, probesNonAmbiguous, ks, ks2 ;
  char *qq = hprintf (h1, 
		      "IS  \"%s%s*\" %s ; %s"
		      , type
		      , mm->isMaggie ? "" : "_"
		      , hasCluster ? " && NOT IS *_at_* " : ""
		      , mm->template ?  mm->template : ""
		      ) ;

  
  probes = ac_dbquery_keyset (mm->db, messprintf ("Find probe %s ", qq), h1) ;
  hp->p = ac_keyset_count (probes) ;
  
  ks = ac_ksquery_keyset (probes,"nm_hit" , h1);
  hp->p_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  if (isNmDb) /* we use exortNmAsGenomic */
    probesNonAmbiguous = ac_ksquery_keyset (probes," nm_hit && (COUNT {>nm_hit ; >geneid} = 1)" , h1); 
  else /* we use exportAceViewAsGenomic */
      probesNonAmbiguous = ac_ksquery_keyset (probes,"nm_hit ; COUNT {>Probe_set ; non_ambiguous} > 0 || COUNT {>nm_hit ; >Model_of_gene} = 1 " , h1);

  hp->p_1g = ac_keyset_count (probesNonAmbiguous) ;

  ks = ac_ksquery_keyset (probesNonAmbiguous,"COUNT nm_hit = 1" , h1);
  hp->p_1_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  ks = ac_ksquery_keyset (probesNonAmbiguous,"COUNT nm_hit = 2" , h1);
  hp->p_2_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  ks = ac_ksquery_keyset (probesNonAmbiguous,"COUNT nm_hit = 3" , h1);
  hp->p_3_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  ks = ac_ksquery_keyset (probesNonAmbiguous,"COUNT nm_hit = 4" , h1);
  hp->p_4_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  ks = ac_ksquery_keyset (probesNonAmbiguous,"COUNT nm_hit = 5" , h1);
  hp->p_5_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  ks = ac_ksquery_keyset (probesNonAmbiguous,"COUNT nm_hit = 6" , h1);
  hp->p_6_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  ks = ac_ksquery_keyset (probesNonAmbiguous,"COUNT nm_hit = 7" , h1);
  hp->p_7_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  ks = ac_ksquery_keyset (probesNonAmbiguous,"COUNT nm_hit = 8" , h1);
  hp->p_8_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  ks = ac_ksquery_keyset (probesNonAmbiguous,"COUNT nm_hit = 9" , h1);
  hp->p_9_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  ks = ac_ksquery_keyset (probesNonAmbiguous,"COUNT nm_hit = 10" , h1);
  hp->p_10_nm = ac_keyset_count (ks) ;
  ac_free (ks) ;
  
  hp->nm_ps = ac_ksquery_keyset (probesNonAmbiguous,">nm_hit" , h);
  hp->nm_p = ac_keyset_count (hp->nm_ps) ;
  
  if (isNmDb) /* we use exortNmAsGenomic */
    {
      hp->g_p1s = ac_ksquery_keyset (probesNonAmbiguous, ">nm_hit ; >geneid", h) ;
      hp->g_p1 = ac_keyset_count (hp->g_p1s) ;
      
      hp->gid_p1 = hp->g_p1 ;

      ks2 = ac_ksquery_keyset (probesNonAmbiguous,  "COUNT {>nm_hit } = COUNT { >nm_hit ; >geneid; >Sequence probe_nm_hit}", h1) ;
      hp->p_1gc = ac_keyset_count (ks2) ;
      ac_free (ks2) ;
      ks2 = ac_ksquery_keyset (probesNonAmbiguous,  "COUNT {>nm_hit } < COUNT { >nm_hit ; >geneid; >Sequence probe_nm_hit}", h1) ;
      hp->p_1ga = ac_keyset_count (ks2) ;
      ac_free (ks2) ;
    }
  else /* we use exportAceViewAsGenomic */
    {
      hp->g_p1s = ac_ksquery_keyset (probesNonAmbiguous, ">nm_hit ; >Model_of_gene", h) ;
      hp->g_p1 = ac_keyset_count (hp->g_p1s) ;
      
      ks2 = ac_ksquery_keyset (hp->g_p1s, ">geneid", 0) ;
      hp->gid_p1 = ac_keyset_count (ks2) ;
      ac_free (ks2) ;
      
      ks2 = ac_ksquery_keyset (probesNonAmbiguous,  "COUNT {>nm_hit } = COUNT { >nm_hit ; >Model_of_gene; >Genefinder probe_nm_hit}", h1) ;
      hp->p_1gc = ac_keyset_count (ks2) ;
      ac_free (ks2) ;
      ks2 = ac_ksquery_keyset (probesNonAmbiguous,  "COUNT {>nm_hit } < COUNT { >nm_hit ; >Model_of_gene; >Genefinder probe_nm_hit}", h1) ;
      hp->p_1ga = ac_keyset_count (ks2) ;
      ac_free (ks2) ;
    }
  ac_free (h1) ;
}

static void magHitStatsReportDoubles (DICT *dict, Array hss, int pass)
{
  AC_HANDLE h1 = ac_new_handle () ;
  AC_KEYSET ks = 0 ;
  HSS *hp, *hp2 ; 
  int nn, nType, nType2 ;
  
  if (pass == 0) printf ("%s", "m_p1") ;
  if (pass == 1) printf ("%s", "g_p1") ;
  
  for (nType = 0; nType < arrayMax (hss); nType++)
    { 
      hp = arrayp (hss, nType, HSS) ;
      printf ("\t%s", dictName (dict, hp->nam)) ;
    }
  printf ("\n") ;
  
  for (nType = 0; nType < arrayMax (hss); nType++)
    {
      hp = arrayp (hss, nType, HSS) ;
      printf ("%s", dictName (dict, hp->nam)) ;
      for (nType2 = 0; nType2 < arrayMax (hss); nType2++)
	{
	  hp2 = arrayp (hss, nType2, HSS) ;
	  if (pass == 0) 
	    {
	      ks = ac_copy_keyset (hp->nm_ps, h1) ;
	      nn = ac_keyset_and (ks, hp2->nm_ps) ;
	    }
	  else if (pass == 1)
	    {
	      ks = ac_copy_keyset (hp->g_p1s, h1) ;
	      nn = ac_keyset_and (ks, hp2->g_p1s) ;
	    }
	  else
	    nn = 0 ;
	  printf ("\t%d", nn) ;
	  ac_free (ks) ;	  
	}
      printf ("\n") ;
    }
  
  ac_free (h1) ;
}

static void magHitStatsReport (const char **types, Array hss)
{
  HSS *hp ; 
  int nType ;
  
  printf ("%s\t%8s\t%8s\t%8s"
	  , "type", "p", "p_m",  "p_1g"
	  ) ;
   printf ("\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s"
	   , "p1g_1_m", "p1g_2_m", "p1g_3_m",  "p1g_4_m", "p1g_5_m"
	   , "p1g_6_m", "p1g_7_m", "p1g_8_m", "p1g_9_m", "p1g_10m"
	  ) ;
   printf ("\t%8s\t%8s\t%8s\t%8s\t%8s\n" 
	   , "m_p1", "p_1gc", "p_1ga", "g_p1", "gid_p1"
	  ) ;
 
  for (nType = 0, hp = arrayp (hss, nType, HSS); types [nType]; hp++,  nType++)
    printf ("%s\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\n"
	    , types [nType]
	    , hp->p
	    , hp->p_nm
	    , hp->p_1g
	    , hp->p_1_nm
	    , hp->p_2_nm
	    , hp->p_3_nm
	    , hp->p_4_nm
	    , hp->p_5_nm
	    , hp->p_6_nm
	    , hp->p_7_nm
	    , hp->p_8_nm
	    , hp->p_9_nm
	    , hp->p_10_nm
	    , hp->nm_p
	    , hp->p_1gc
	    , hp->p_1ga
	    , hp->g_p1  
	    , hp->gid_p1
	    ) ;
} /* magHitStatsReport */

static void magHitStatsCreateDoubles (DICT *dict, Array hss, int stop, AC_HANDLE h)
{
  int i, j ;
  char *cp ;
  HSS *h1, *h2, *h3 ;

  h3 = arrayp (hss, arrayMax(hss), HSS) ;
  dictAdd (dict, "ALL", &(h3->nam)) ;

  i = 0 ;
  h1 = arrayp (hss, i, HSS) ;
  h3->nm_ps = ac_copy_keyset (h1->nm_ps, h) ;
  h3->g_p1s = ac_copy_keyset (h1->g_p1s, h) ;
  for (j = i + 1 ; j < stop ; j++)
    {
      h2 = arrayp (hss, j, HSS) ;
      ac_keyset_and (h3->nm_ps, h2->nm_ps) ;
      ac_keyset_and (h3->g_p1s, h2->g_p1s) ;
    }

  for (i = 0 ; i < stop ; i++)
    for (j = i + 1 ; j < stop ; j++)
      {
	h3 = arrayp (hss, arrayMax(hss), HSS) ;
	h1 = arrayp (hss, i, HSS) ;
	h2 = arrayp (hss, j, HSS) ;
	cp = messprintf ("%s/%s", dictName (dict, h1->nam), dictName (dict, h2->nam)) ;
	dictAdd (dict, cp, &(h3->nam)) ;
	h3->nm_ps = ac_copy_keyset (h1->nm_ps, h) ;
	ac_keyset_and (h3->nm_ps, h2->nm_ps) ;
	h3->g_p1s = ac_copy_keyset (h1->g_p1s, h) ;
	ac_keyset_and (h3->g_p1s, h2->g_p1s) ;
      }
} /* magHitStatsCreateDoubles */


static void magHitStats (MM *mm, BOOL isNmDb)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *types[] = { "ABI", "AFX", "AGL", "GEH", "ILM",  "NCI", "MLT", "GEX", "QGN", "TAQ", "EPP", 0} ;
  BOOL hasCluster[] =   { FALSE,  TRUE, FALSE, FALSE, FALSE, FALSE,  FALSE,  FALSE, FALSE,  FALSE, FALSE, FALSE } ;
  Array hss = arrayCreate (32, HSS) ;
  HSS *hp ;
  int nType ;
  DICT *dict = dictHandleCreate (60, h) ;
  
  if (mm->isMaggie) 
    {
      types[0] = "Affy_" ;
      types[1] = "GE_" ;
      types[2] = "Agilent_" ;
      types[3] = 0 ;
      hasCluster[0] = FALSE ;
      hasCluster[1] = FALSE ;
      hasCluster[2] = TRUE ;
      hasCluster[3] = FALSE ;
    }
  /* count the types */
  /* count the probes, probes with hits, genes in each type */
  for (nType = 0 ; types [nType]; nType++)
    {
      hp = arrayp (hss, nType, HSS) ;
      dictAdd (dict, types [nType], &(hp->nam)) ;
      magOneHitStats (mm, hp, types [nType], hasCluster[nType], isNmDb, h) ;
    }
  /* report */
  magHitStatsReport (types, hss) ;
  /* report compatibilities */
  magHitStatsCreateDoubles (dict, hss, mm->isMaggie ? 3 : 6, h) ;
  magHitStatsReportDoubles (dict, hss, 0) ;
  magHitStatsReportDoubles (dict, hss, 1) ;

  ac_free (h) ;
}

/*************************************************************************************/
/*************************************************************************************/
/* parse files 
 * tg_with_probe2assembled_from.ace  
 * tg_with_probe2mrna.ace
 *
 * export for each gene the number of EST for each probe
 */

typedef struct hsaStruct { 
  int mrna ;  /* index of the mrna name in dict */
  int tgNam ; /* index of the tg name in dict */
  int estNam ; /* index of the tg name in estDict */
  BOOL isNm ; /* this est is an NM */
  int m1 ;    /* offset of the mRNA in the tg */
  int a1, a2 ; /* est exon support */
} HSA ;

/****************************************/

static BOOL hitSupportAceViewParseOffset (DICT *dict, Array tgs, Array aa)
{
  BOOL ok = FALSE ;
  int mrna, imrna = 0, tgNam, m1, level ;
  int nmrna = 0, nm1 = 0, ntg = 0 ;
  char *cp ;
  HSA *hsa, *tgsa = 0 ;
  FILE *f = filopen ("tg_with_probe2mrna", "ace", "r") ;

  if (!f) 
    {
      messout ("tg_with_probe2mrna.ace should contain the mRNA offset in the tg\n") ;
      return FALSE ;
    }
  level = freesetfile (f, 0) ;
  while (ntg >= 0 && freecard (level))
    {
      cp = freeword () ;
      if (cp && ! strcmp (cp, "Transcribed_gene") && freeword () && (cp = freeword ()))
	{
	  ntg++ ;
	  dictAdd (dict, cp, &tgNam) ;
	  tgsa = arrayp (tgs, tgNam, HSA) ;
	  tgsa->tgNam = tgNam ;
	  tgsa->m1 = imrna ;
	  continue ;
	}
      if (cp && ! strcmp (cp, "mRNA") && (cp = freeword ()))
	{
	  nmrna++ ;
	  dictAdd (dict, messprintf("GM_%s", cp), &mrna) ;
	  hsa = arrayp (aa, imrna++, HSA) ;
	  hsa->mrna = mrna ;
	  hsa->tgNam = tgNam ;
	  if (freeint (&m1))
	    { ok = TRUE ; hsa->m1 = m1 ; nm1++ ;}
	}
    }
  freeclose (level) ; /* will close f */
  printf ("Found %d names %d tg %d mRNA %d offset\n", dictMax (dict), ntg, nmrna, nm1) ;

  return ok ;
} /* hitSupportAceViewParseOffset */

/****************************************/

static void hitSupportAceViewReport (MM *mm, DICT *dict, int iTg2Mrna, int tgNam, Array aa, Array af)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *probeDict = dictHandleCreate (64, h) ;
  Array support = arrayHandleCreate (100, int, h) ;
  Array nmSupport = arrayHandleCreate (100, int, h) ;
  int i, ir, jr, iaa, iaf, p1, p2, q1, q2, b1, b2, a1, a2, x1, x2 ;
  int nEstInGene, nNmInGene, nMrnaInGeneWithHit, nTypesInGene ;
  int nEstSupport, nEstHitByAProbe, nNmSupport, nMrnaHits, nExonProbeExon ;
  BOOL isUp ;
  HSA *hsa, *afp ; 
  AC_TABLE gIds, probes, exons, hits ;
  AC_OBJ mrna = 0 ;
  AC_OBJ gene = 0 ;
  AC_OBJ probe = 0 ;

  /* loop in the gene to count the EST and the NMs */
  support = arrayReCreate (support, 100, int) ;
  nmSupport = arrayReCreate (nmSupport, 100, int) ;
  nEstInGene = nEstHitByAProbe = nNmInGene =  nMrnaInGeneWithHit = nTypesInGene = 0 ;
  for (iaf = 0 ; iaf < arrayMax (af) ; iaf++)
    {
      afp = arrp (af, iaf, HSA) ;
      if (!array (support, afp->estNam, int))
	{
	  nEstInGene ++ ; 
	  array (support, afp->estNam, int) = 1 ;
	  if (afp->isNm)
	    nNmInGene++ ;
	}
    }
 
  /* loop on all mRNAs of this gene */
  for (iaa = iTg2Mrna ; iaa < arrayMax (aa) ; iaa++)
    {
      hsa = arrp (aa, iaa, HSA) ;
      if (hsa->tgNam != tgNam)
	break ;
      mrna = ac_get_obj (mm->db, "Sequence", dictName (dict, hsa->mrna), h) ;
      if (!mrna || ! ac_has_tag (mrna, "probe_nm_hit"))
	continue ;
      gene = ac_tag_obj (mrna, "Model_of_gene", h) ;
      if (!gene)
	continue ;
      exons = ac_tag_table (mrna, "Exon_coord", h) ;
      if (!exons)
	continue ; 
      /* count all mRnas of this gene with a probe hit */
      nMrnaInGeneWithHit = ac_keyset_count (ac_objquery_keyset (gene, ">Genefinder probe_nm_hit", h)) ;
      nTypesInGene = ac_keyset_count (ac_objquery_keyset (gene, ">Genefinder ; >probe_nm_hit ; >Probe_Set ; >library", h)) ;
      gIds = ac_tag_table (gene, "GeneId", h) ;
      
      /* loop on all probes of this mRNA */
      probes = ac_tag_table (mrna, "probe_nm_hit", h) ;
      for (ir = 0 ; probes  && ir < probes->rows ; ir++)
	{
	  probe = ac_table_obj (probes, ir, 0, h) ; 
	  if (dictFind (probeDict, ac_name (probe), 0))
	    continue ;
	  /* avoid ambiguous probes */
	  {
	    AC_KEYSET kks = ac_objquery_keyset (probe
				     , "COUNT {>Probe_set ; non_ambiguous} > 0 || COUNT {>nm_hit; >Model_of_gene} = 1"
				     , h) ;
	    int nkks = ac_keyset_count (kks) ;
	    ac_free (kks) ;
	    if (!nkks)
	      {
		dictAdd (probeDict, ac_name (probe), 0) ;
		continue ;
	      }
	  }
	  hits = ac_tag_table (probe, "nm_hit", h) ;
	  if (! hits)
	    continue ;
	  q1 = q2 = 0 ;
	  nExonProbeExon = 0 ;
	  for (jr = 0 ; jr < hits->rows ; jr++)
	    {
	      if (ac_obj_key (mrna) != ac_table_key (hits, jr, 0, 0))
		continue ;
	      q1 = ac_table_int (hits, jr, 1, 0) ; /* coords on spliced mrna */
	      q2 = ac_table_int (hits, jr, 2, 0) ; 
	      break ;
	    }
	  isUp = FALSE ; 
	  if (q1 > q2) { int q3 = q1 ; q1 = q2 ; q2 = q3 ; isUp = TRUE ; }
	  if (q1 == q2)
	    continue ;
	  support = arrayReCreate (support, 100, int) ;
	  nmSupport = arrayReCreate (nmSupport, 100, int) ;
	  nMrnaHits = ac_keyset_count (ac_objquery_keyset (probe
							   , messprintf (">nm_hit ; Model_of_gene = \"%s\"", ac_name(gene))
							       , h)) ;
	  while (q1 < q2) /* we may hit several exons */
	    {
	      p1 = 1 ; p2 = 0 ;
	      for (jr = 0 ; jr < exons->rows ; jr++)
		{
		  a1 = ac_table_int (exons, jr, 0, 0) ;
		  a2 = ac_table_int (exons, jr, 1, 0) ;
		  x1 = ac_table_int (exons, jr, 2, 0) ;
		  x2 = ac_table_int (exons, jr, 3, 0) ;
		  b1 = x1 < q1 ? q1 : x1 ;
		  b2 = x2 > q2 ? q2 : x2 ;
		  if (b1 <= b2) 
		    {
		      p1 = a1 + q1 - x1 ;
		      p2 = p1 + b2 - b1 ;
		      q1 = q1 + b2 - b1 + 1 ;
		      nExonProbeExon++ ;
		      break ;
		    }
		}
	      if (p1 > p2)
		break ;
	      
	      p1 = p1 + hsa->m1 - 1 ;
	      p2 = p2 + hsa->m1 - 1 ;

	      /* loop in the gene support to find the probe support */
	      for (iaf = 0 ; iaf < arrayMax (af) ; iaf++)
		{
		  afp = arrp (af, iaf, HSA) ;
		  b1 = afp->a1 < p1 ? p1 : afp->a1 ;
		  b2 = afp->a2 > p2 ? p2 : afp->a2 ;
		  if (b1 < b2) /* the hit is supported */
		    {
		      array (support, afp->estNam, int) = 1 ;
		      if (afp->isNm)
			array (nmSupport, afp->estNam, int) = 1 ;
		    }
		}
	    }

	  dictAdd (probeDict, ac_name (probe), 0) ;
	  /* count the support */
	  for (nEstSupport = i = 0 ; i < arrayMax (support) ; i++)
	    if (arr (support, i, int)) nEstSupport++ ;
	  for (nNmSupport = i = 0 ; i < arrayMax (nmSupport) ; i++)
	    if (arr (nmSupport, i, int)) nNmSupport++ ;

	  /* report */
	  printf ("%d\t", nTypesInGene) ;
	  printf ("%s\t", ac_name(gene)) ;
	  /* ; separated geneIds */
	  for (i = 0 ; gIds && i < gIds->rows ; i++)
	    printf ("%s%s"
		    , i ? ";" : ""
		    , ac_table_printable (gIds, i, 0, "")
		    ) ;
	  printf ("\t") ;
	  /* gene support */ 
	  { 
	    int ntg = 52954, ncl = 4018074 ;
	    float magic = (ncl + 1.0)/(ntg + 1.0) ;
	    printf ("%d\t%d\t%2.2f%%\t%d\t"
		    ,  nMrnaInGeneWithHit, nEstInGene
		    , 100.0 * nEstInGene/magic
		    , nNmInGene
		    ) ;
	  }
	  /* probe support */
	  printf ("%s\t", ac_name(probe)) ;
	  printf ("%s\t", isUp ? "-" : "+") ;
	  printf ("%d\t", nExonProbeExon - 1) ;
	  {
	    int v1, v2, v3 ;

	    v1 = nMrnaInGeneWithHit ? nMrnaInGeneWithHit : 1 ;
	    v2 = nEstInGene ? nEstInGene : 1 ; /* nEstHitByAProbe ? nEstHitByAProbe : 1 ; */
	    v3 = nNmInGene ? nNmInGene : 1 ;
	    printf ("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
		    , nMrnaHits
		    , (10000 * nMrnaHits/v1 + 100 * (v1 - 1)/v1)/100
		    , nEstSupport
		    , (10000 * nEstSupport/v2 +  100 * (v2 - 1)/v2)/100
		    , nNmSupport
		    , (10000 * nNmSupport/v3  +  100 * (v3- 1)/v3)/100
		    , nEstHitByAProbe, nEstInGene
		    ) ;
	  }	  
	}      
    }

  ac_free (h) ;
} /* hitSupportAceViewReport  */

/****************************************/

static BOOL hitSupportAceViewAnalyse (MM *mm, DICT *dict, Array tgs, Array aa)
{
  BOOL ok = FALSE ;
  int tgNam, estNam, a1, a2, level, iTg2Mrna = 0 ;
  int ntg = 0, iaf = 0 ;
  char *cp ;
  HSA *afp ;
  Array af = 0 ; /* array of the gene support */
  DICT *estDict = 0 ;
  FILE *f = filopen ("tg_with_probe2assembled_from", "ace", "r") ;

  if (!f) 
    {
      messout ("tg_with_probe2assembled_from should contain the support of the tg\n") ;
      return FALSE ;
    }
  level = freesetfile (f, 0) ;
  while (ntg >= 0 && freecard (level))
    {
      cp = freeword () ;
      if (cp && ! strcmp (cp, "Transcribed_gene") && freeword () && (cp = freeword ()))
	{
	  if (iaf)
	    hitSupportAceViewReport (mm, dict, iTg2Mrna, tgNam, aa, af) ;
	  tgNam = 0 ;
	  ntg++ ;
	  dictFind (dict, cp, &tgNam) ;
	  af = arrayReCreate (af, 1000, HSA) ;
	  iaf = 0 ;
	  dictDestroy (estDict) ;
	  estDict = dictHandleCreate (1000, 0) ;
	  iTg2Mrna = -1 ;
	  if (tgNam >= 0 && tgNam < arrayMax(tgs))
	    iTg2Mrna = arr (tgs, tgNam, HSA).m1 ;
	  continue ;
	}
      if (cp && tgNam && iTg2Mrna >= 0 && ! strcmp (cp, "Assembled_from"))
	{ 
	  if (freeint (&a1) && freeint (&a2) && (cp = freeword()))
	    {
	      dictAdd (estDict, cp, &estNam) ;
	      afp = arrayp (af, iaf++, HSA) ;
	      afp->estNam = estNam ;
	      afp->a1 = a1 ;
	      afp->a2 = a2 ;
	      afp->isNm = (!strncasecmp (cp, "NM_", 3) ? TRUE : FALSE) ;
	    }
	  continue ;
	}
    }
  freeclose (level) ; /* will close f */
  printf ("Found %d names %d tg\n", dictMax (dict), ntg) ;

  arrayDestroy (af) ; 
  dictDestroy (estDict) ;
  return ok ;
} /* hitSupportAceViewAnalyse */

/****************************************/

static void hitSupportAceView (MM *mm, BOOL isNmDb)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *dict = dictHandleCreate (80000, h) ;
  Array aa = arrayHandleCreate (80000, HSS, h) ;
  Array tgs = arrayHandleCreate (80000, HSS, h) ;
  
  dictAdd (dict, "toto", 0) ;
  if (hitSupportAceViewParseOffset (dict, tgs, aa))
    hitSupportAceViewAnalyse (mm, dict, tgs, aa) ;

  ac_free (h) ;
} /* hitSupportAceView */

/*************************************************************************************/
/*************************************************************************************/

static void flagAmbiguousOneProbe_Set (MM *mm, AC_OBJ Probe_set)
{
  AC_HANDLE h = ac_new_handle () ; 
  AC_OBJ gid = 0, probe = 0, mrna = 0 ;
  AC_KEYSET genes, ks ;
  AC_ITER gids, probes ;
  AC_TABLE hits = 0, intmap = 0, exons = 0 ; 
  int ir, jr, nn ;
  int a1, a2, x1, x2, m1, m2, p1, p2, q1, q1Old ;
  KEY chrom, chromOld ;
  BOOL ok ;
  
  /* may be all the genes contain the same NM*/
  genes = ac_objquery_keyset (Probe_set, ">Probe ; >Nm_Hit ; >Model_of_gene", h) ;
  nn = ac_keyset_count (genes) ;
  gids = ac_ksquery_iter (genes, ">geneid", h) ;

  ok = FALSE ;
  while (gid = 0 , !ok && (gid = ac_iter_obj (gids)))
    {
      ks = ac_objquery_keyset (gid, ">Gene", h) ;
      if (nn == ac_keyset_count (ks))
	{
	  ok = TRUE ;  /* this Probe_set is NOT ambiguous */
	  printf ("PROBE_SET %s\nNon_ambiguous Repeated_gene\n\n"
		  , freeprotect (ac_name(Probe_set))
		  ) ;
	  goto done ;
	}
    }

  if (mm->isNmDb)
    goto done ;

  /* may be coordinates are the same */
  probes = ac_objquery_iter (Probe_set, ">probes", h) ;
  while (ac_free (hits), ac_free (probe), !ok && (probe = ac_iter_obj (probes)))
    {
      hits = ac_tag_table (probe, "nm_hit", h) ;
      ok = TRUE ;
      for (q1Old = chromOld = ir = 0 ; hits && ok && ir < hits->rows 
	     ; ac_free (mrna), ac_free (exons), ac_free (intmap),  ir++)
	{
	  mrna = ac_table_obj (hits, ir, 0, h) ;
	  p1 = ac_table_int (hits, ir, 1, 0) ;
	  p2 = ac_table_int (hits, ir, 2, 0) ;
	  if (p1 > p2)
	    { ok = FALSE ; continue ; }
	  intmap = ac_tag_table (mrna, "IntMap", h) ;
	  m1 = m2 = 0 ;
	  if (! intmap)
	    { ok = FALSE ; continue ; }
	  chrom = ac_table_key (intmap, 0, 0, 0) ;
	  m1 = ac_table_int (intmap, 0, 1, 0) ;
	  m2 = ac_table_int (intmap, 0, 2, 0) ;
	  if (!m2)
	    { ok = FALSE ; continue ; }
	  if (! chromOld) chromOld = chrom ;
	  if (chrom != chromOld)
	    { ok = FALSE ; continue ; }
	  exons = ac_tag_table (mrna, "Exon_coord", h) ;
	  if (! exons)
	    { ok = FALSE ; continue ; }
	  for (q1 = jr = 0 ; jr < exons->rows ; jr++)
	    {  
	      a1 = ac_table_int (exons, jr, 0, 0) ;
	      a2 = ac_table_int (exons, jr, 1, 0) ;
	      x1 = ac_table_int (exons, jr, 2, 0) ;
	      x2 = ac_table_int (exons, jr, 3, 0) ;
	      if (p1 >= x1 && p2 <= x2)
		{ q1 = a1 + p1 - x1 ; break ; }
	    }
	  if (m1 < m2) q1 = m1 + q1 - 1 ;
	  else q1 = m1 - q1 + 1 ;
	  if (! q1Old) q1Old = q1 ;
	  if (q1 != q1Old) /* not the same match */
	    { ok = FALSE ; continue ; }
	}
      if (ok)
	printf ("PROBE_SET %s\nNon_ambiguous Shed_gene\n\n"
		, freeprotect (ac_name(Probe_set))
		) ;
      goto done ;
    }

 done:
  ac_free (h) ;
} /*flagAmbiguousOneProbe_Set */

/*************************************************************************************/

static void flagAmbiguousProbe (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ Probe_set = 0 ;
  AC_ITER iter = ac_dbquery_iter (mm->db
				  , messprintf ("Find Probe_set %s ; COUNT {>Probe ; >nm_hit ; >%s} > 1"
						, mm->template ? mm->template  : ""
						, mm->isNmDb ? "geneId" : "model_of_gene"
						)
				  , h) ;

  while (ac_free (Probe_set), (Probe_set = ac_iter_obj (iter)))
    flagAmbiguousOneProbe_Set (mm, Probe_set) ;

  ac_free (h) ;
} /* flagAmbiguousProbe */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static void flagAmbiguousOneProbeAceView (MM *mm, AC_OBJ Probe)
{
  AC_HANDLE h = ac_new_handle () ; 
  AC_OBJ mrna = 0, bestMrna ;
  AC_KEYSET genes, genes2, gids = 0, ks, gidmrnas = 0 ;
  AC_TABLE tbl, tbl2, hits = 0, intmap = 0, exons = 0 ; 
  int ir, jr, nn, n2, n3, myClass = 0 ;
  int a1, a2, x1, x2, m1, m2, bestDm, bestErr, p1, p2, q1, q1Old ;
  KEY chrom, chromOld ;
  BOOL ok, bestMrnaHasGid = FALSE ;
  int probeLength = ac_tag_int (Probe, "Length", 9999) ;

  /* may be all the genes contain the same NM*/
  ok = FALSE ;
  freeOutf ("Probe %s\n-D Type\n", freeprotect (ac_name(Probe))) ;
  genes = ac_objquery_keyset (Probe
			      , "{>Nm_exact_Hit ;>In_mrna} SETOR {>mrna_exact_hit}; >from_gene;>gene"
			      , h) ;
  if (0) /* feb 11 2006, export to MAQC, we no longer wish to salvage based on geneid */
    gids = ac_ksquery_keyset (genes, ">geneid", h) ;
  tbl = ac_keyset_table (genes, 0, -1, FALSE, h) ;
  nn = tbl ? tbl->rows : 0 ; 

  freeOutf ("-D Single_exact_geneid\n") ;
  freeOutf ("-D Exact_gene\n") ;
  for (ir = 0 ; ir < nn ; ir++)
    freeOutf ("Exact_gene %s\n", freeprotect (ac_table_printable (tbl, ir, 0, ""))) ;
  if (nn == 1)
    {
      freeOutf ("Single_exact_gene %s\n", freeprotect (ac_table_printable (tbl, 0, 0, ""))) ;
      freeOutf ("Confirmed_gene %s\n", freeprotect (ac_table_printable (tbl, 0, 0, ""))) ;
      ok = TRUE ;
      myClass = 10 ;
    }

  n2 = gids ? ac_keyset_count (gids) : 0 ;
  if (n2 == 1)
    {
      if (nn == 1)
	{
	  tbl2 = ac_keyset_table (gids, 0, -1, FALSE, h) ;
	  freeOutf ("Single_exact_geneid %s\n"
		    , freeprotect (ac_table_printable (tbl2, 0, 0, ""))) ;
	}
      else
	{
	  genes2 = ac_ksquery_keyset (gids, ">gene", h) ;
	  n3 = ac_keyset_count (genes2) ;
	  if (n3 == nn)
	    {
	      tbl = ac_keyset_table (gids, 0, -1, FALSE, h) ;
	      freeOutf ("Single_exact_geneid %s\n"
			, freeprotect (ac_table_printable (tbl, 0, 0, ""))) ;
	      freeOutf ("Confirmed_gene %s\n"
			, freeprotect (ac_table_printable (tbl, 0, 0, ""))) ;
	      ok = TRUE ; /* this Probe is NOT ambiguous */ 
	      myClass = 110 ;
	    }
	}
    }
  
  
  /* register the coordinates and see if they are the same */
  if (1)
     {
       int nr, ns, p1, p2, q1, q2, g1, g2, errpos ;

       bestErr = 9999 ;
       ac_free (mrna) ; ac_free (exons) ;
       hits = ac_tag_table (Probe, "mrna_hit", h) ;
       for (ir = 0 ; hits && ir < hits->rows 
	      ; ac_free (mrna), ac_free (exons), ir++)
	 {
	   mrna = ac_table_obj (hits, ir, 0, h) ;
	   exons = ac_tag_table (mrna, "Splicing", h) ;
	   if (! exons)
	     continue ; 
	   p1 = ac_table_int (hits, ir, 1, -100) ;
	   p2 = ac_table_int (hits, ir, 2, -100) ;
	   nr = ac_table_int (hits, ir, 3, 0) ;
	   ns = ac_table_int (hits, ir, 4, 0) ;
	   g1 = ac_table_int (hits, ir, 5, 0) ;
	   g2 = ac_table_int (hits, ir, 6, 0) ;
	   errpos = ac_table_int (hits, ir, 7, 0) ;

	   if (nr < bestErr) 
	     {
	       if (ok && ac_keyset_count (genes))
		 {
		   AC_KEYSET myGene = ac_objquery_keyset (mrna, ">from_gene ; >gene", h) ;
		   if (!ac_keyset_and (myGene, genes)) /* this is not one of the exact genes */ 
		     bestErr = nr ;
		   ac_free (myGene) ;
		 }
	       else
		 bestErr = nr ;
	     }
	   if (p2 < 0) 
	     continue ;
	   for (q1 = q2 = -100, jr = 0 ; jr < exons->rows ; jr++)
	     {  
	       a1 = ac_table_int (exons, jr, 0, 0) ;
	       a2 = ac_table_int (exons, jr, 1, 0) ;
	       x1 = ac_table_int (exons, jr, 2, 0) ;
	       x2 = ac_table_int (exons, jr, 3, 0) ;
	       if (p1 >= x1 && p1 <= x2)
		 { q1 = a1 + p1 - x1 ; break ; }
	     }
	   for ( ; jr < exons->rows ; jr++)
	     {  
	       a1 = ac_table_int (exons, jr, 0, 0) ;
	       a2 = ac_table_int (exons, jr, 1, 0) ;
	       x1 = ac_table_int (exons, jr, 2, 0) ;
	       x2 = ac_table_int (exons, jr, 3, 0) ;
	       if (p2 >= x1 && p2 <= x2)
		 { q2 = a1 + p2 - x1 ; break ; }
	     }
	   
	   if (q1 > 0 && q2 > 0)
	     freeOutf ("mrna_hit %s %d %d %d %d %d %d %d\n"
		       , freeprotect (ac_name (mrna))
		       , p1, p2, nr, ns, q1, q2, errpos
		       ) ;
	 }
     }

  /* may be the coordinates are the same */
  bestDm = 0 ; bestMrna = 0 ;
  if (!ok && nn > 1)
    {
      hits = ac_tag_table (Probe, "mrna_exact_hit", h) ;
      ok = TRUE ;
      for (q1Old = chromOld = ir = 0 ; hits && ok && ir < hits->rows 
	     ; ac_free (exons), ac_free (intmap),  ir++)
	{
	  mrna = ac_table_obj (hits, ir, 0, h) ;
	  p1 = ac_table_int (hits, ir, 1, 0) ;
	  p2 = ac_table_int (hits, ir, 2, 0) ;
	  if (p1 > p2)
	    { ok = FALSE ; continue ; }
	  intmap = ac_tag_table (mrna, "IntMap", h) ;
	  m1 = m2 = 0 ;
	  if (! intmap)
	    { ok = FALSE ; continue ; }
	  if (strstr (ac_table_printable (intmap, 0, 0, ""), "|NT"))
	    continue ; /* ignore NT hits */
	  if (strstr (ac_table_printable (intmap, 0, 0, ""), "|Hs"))
	    continue ; /* ignore NT hits */
	  chrom = ac_table_key (intmap, 0, 0, 0) ;
	  m1 = ac_table_int (intmap, 0, 1, 0) ;
	  m2 = ac_table_int (intmap, 0, 2, 0) ;
	  if (!bestMrnaHasGid &&
	      (
	       (m1 < m2 && m2 - m1 > bestDm) ||
	       (m1 > m2 && m1 - m2 > bestDm)
	       )
	      )
	    { bestDm = m1 < m2 ? m2 - m1 : m1 - m2 ; bestMrna = mrna ; }
	  if (!bestMrnaHasGid)
	    {
	      gidmrnas = ac_objquery_keyset (mrna, "COUNT {>from_gene;>gene; geneid} > 0", h) ;
	      if (gidmrnas && ac_keyset_count (gidmrnas))
		{
		  bestMrnaHasGid = TRUE ;
		  bestDm = m1 < m2 ? m2 - m1 : m1 - m2 ; bestMrna = mrna ;
		}
	    }
	  if (!m2)
	    { ok = FALSE ; continue ; }
	  if (! chromOld) chromOld = chrom ;
	  if (chrom != chromOld)
	    { ok = FALSE ; continue ; }
	  exons = ac_tag_table (mrna, "Splicing", h) ;
	  if (! exons)
	    { ok = FALSE ; continue ; }
	  for (q1 = jr = 0 ; jr < exons->rows ; jr++)
	    {  
	      a1 = ac_table_int (exons, jr, 0, 0) ;
	      a2 = ac_table_int (exons, jr, 1, 0) ;
	      x1 = ac_table_int (exons, jr, 2, 0) ;
	      x2 = ac_table_int (exons, jr, 3, 0) ;
	      if (p1 >= x1 && p1 <= x2)
		{ q1 = a1 + p1 - x1 ; break ; }
	    }
	  if (m1 < m2) q1 = m1 + q1 - 1 ;
	  else q1 = m1 - q1 + 1 ;
	  if (! q1Old) q1Old = q1 ;
	  if (q1 != q1Old) /* not the same match */
	    { ok = FALSE ; continue ; }
	}
      if (ok && bestMrna)
	{
	  ks = ac_objquery_keyset (bestMrna, ">from_gene ; > gene", h) ;
	  tbl2 = ac_keyset_table (ks, 0, -1, FALSE, h)  ;
	  if (tbl2 && tbl2->rows == 1)
	    {
	      freeOutf ("Single_exact_gene %s\n"
			, freeprotect (ac_table_printable (tbl2, 0, 0, ""))
			/* 	, g1, g2 */
			) ;
	      freeOutf ("Confirmed_gene %s\n"
			, freeprotect (ac_table_printable (tbl2, 0, 0, ""))) ;
	      myClass = 10 ;
	    }
	}
    }
  /* look for exon_exon_hit accross an intron scar */
  freeOutf ("-D Last_exon_hit\nExon_exon_hit\n") ;
  if (1)
    {
      hits = ac_tag_table (Probe, "mrna_hit", h) ;
      
      for (ir = 0 ; hits && ir < hits->rows 
	     ; ac_free (mrna), ac_free (exons), ir++)
	{
	  mrna = ac_table_obj (hits, ir, 0, h) ;
	  p1 = ac_table_int (hits, ir, 1, 0) ;
	  p2 = ac_table_int (hits, ir, 2, 0) ;
	  if (p1 > p2)
	    continue ; 
	  exons = ac_tag_table (mrna, "Splicing", h) ;
	  if (! exons)
	    continue ; 
	  for (jr = 0 ; jr < exons->rows - 1 ; jr++)
	    {  
	      x1 = ac_table_int (exons, jr, 2, 0) ;
	      x2 = ac_table_int (exons, jr, 3, 0) ;
	      if (p1 >= x1 && p1 <= x1 && p2 > x2)
		{ 
		  freeOutf ("Exon_exon_hit %s\n", freeprotect (ac_name (mrna))) ;
		  break ; 
		}
	    }
	  jr =  exons->rows - 1 ;
	  x1 = ac_table_int (exons, jr, 2, 0) ;
	  x2 = ac_table_int (exons, jr, 3, 0) ;
	  if (p1 >= x1 && p1 <= x1 && p2 > x2)
	    {
	      ks = ac_objquery_keyset (mrna, ">Product ; best_product && Down_stop", h) ;
	      if (ac_keyset_count (ks) == 1)
		freeOutf ("Last_exon_hit %s\n", freeprotect (ac_name (mrna))) ;
	    }
	}
    }

  if (! probeLength) /* probe length often missing  */
    {
      const char *cp = ac_tag_printable (Probe, "Motif", 0) ;
      if (cp && (probeLength = strlen (cp)))
	freeOutf ("Length %d\n", ir) ;
    } 

  genes2 = ac_objquery_keyset (Probe, "{>Nm_Hit ;>In_mrna} SETOR {>mrna_hit}; >from_gene;>gene", h) ;
  n2 = ac_keyset_count (genes2) ;
  if (n2 > 1 || nn > 0)
    {
      AC_KEYSET genes3 = ac_ksquery_keyset (genes2, "NOT IntMap = \"*|NT*\" && NOT IntMap = \"*|Hs*\"", h) ;
      int n3 = ac_keyset_count (genes3) ;
      if (nn > 0 || n3 > 0)
	genes2 = genes3 ;
    }
  ac_keyset_minus (genes2, genes) ;
  tbl = ac_keyset_table (genes2, 0, -1, FALSE, h) ;
  n2 = tbl ? tbl->rows : 0 ;
  freeOutf ("-D Approximate_gene\n") ;
  freeOutf ("-D Single_approximate_gene\n") ;
  for (ir = 0 ; ir < n2 ; ir++)
    freeOutf ("Approximate_gene %s\n", freeprotect (ac_table_printable (tbl, ir, 0, ""))) ;


  if (!ok && nn == 0 && n2 == 1)
    {
      freeOutf ("Single_approximate_gene %s\n", freeprotect (ac_table_printable (tbl, 0, 0, ""))) ;
      freeOutf ("Confirmed_gene %s\n", freeprotect (ac_table_printable (tbl, 0, 0, ""))) ;
      if (20 * bestErr <  probeLength) /* bestErr is the number of errors of bestNonExact */
	myClass = 11 ;
      else
	myClass = 12 ;
    }
  else if (!ok && nn > 1) /* several solutions */
    myClass = 30 ;
  else if (ok && n2 >= 1) /* single best exact but also some approximates */
    {
      if (20 * bestErr <  probeLength) /* bestErr is the number of errors of bestNonExact */
	myClass = 31 ;
      else
	myClass = 32 ;
    }
  else if (!ok && n2 > 1) /* several solutions */
    {
      if (20 * bestErr <  probeLength) /* bestErr is the number of errors of bestNonExact */
	{ 
	  int secondBestErr = bestErr ;
	  if (5 * secondBestErr <  probeLength)
	    myClass = 31 ;
	  else
	    myClass = 32 ; 
	}
      else
	myClass = 33 ;
    }
  if (myClass >=  10 && myClass <= 12 &&
       ! ac_has_tag (Probe, "mRNA_hit") && 
       ac_has_tag (Probe, "NM_exact_hit")
       )
    myClass = 11 ;
  if ((myClass == 10 || myClass == 11) &&
       ! ac_has_tag (Probe, "mRNA_hit") && 
       ac_has_tag (Probe, "NM_hit")
       )
    myClass = 12 ;

  if (! myClass && 
       ac_has_tag (Probe, "mrna_anti_hit"))
    myClass = 50 ;
  
  if ( ! myClass &&
	    ( ac_has_tag (Probe, "Genome_exact_hit") ||
	      ac_has_tag (Probe, "Genome_approximate_hit"))
	    )
    {
      ir = 0 ;
      tbl = ac_tag_table (Probe, "Genome_exact_hit", h) ;
      ir += tbl ? tbl->rows : 0 ;
      tbl = ac_tag_table (Probe, "Genome_approximate_hit", h) ;
      ir += tbl ? tbl->rows : 0 ;
      if (ir > 0)
	myClass = 40 ; 
    }
  if ( ! myClass)
    myClass = 60 ; /* unmapped */

  if (1) /* signal analysis */
    {
      int sa, sc, sd, sb ;
      
      sa = ac_tag_float (Probe, "Signal_A", -99999) ;
      sc = ac_tag_float (Probe, "Signal_C", -99999) ;
      sd = ac_tag_float (Probe, "Signal_D", -99999) ;
      sb = ac_tag_float (Probe, "Signal_B", -99999) ;
      
      if (sa > -90000 && sc > -90000 && sd > -90000 && sb > -90000)
	{
	  if (sa > sc && sc > sd && sd > sb)
	    freeOutf ("ACDB_plus\n") ;
	  if (sa < sc && sc < sd && sd < sb)
	    freeOutf ("ACDB_minus\n") ;
	}
    }
  ir = myClass % 100 ; ir /= 10 ;
  if (myClass > 0)
    freeOutf ("Probe_class %d\nProbe_grade%d\n", myClass, ir) ;

  ac_free (h) ;
  freeOutf ("\n") ;
} /* flagAmbiguousOneProbeAceView */

/*************************************************************************************/

static void flagAmbiguousProbeAceView (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ Probe = 0 ;
  AC_ITER iter = ac_dbquery_iter (mm->db
				  , messprintf ("Find Probe ; IS %s"
						, mm->template ? mm->template  : "*"
						)
				  , h) ;

  while (ac_free (Probe), (Probe = ac_iter_obj (iter)))
    flagAmbiguousOneProbeAceView (mm, Probe) ;

  ac_free (h) ;
} /* flagAmbiguousProbe */

/*************************************************************************************/
/*************************************************************************************/
/* remove the hits to NT contigs if there exist a hit to a non NT contig */
static void flagRemoveOneNtHits (AC_OBJ mrna)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ mrna2 = 0, probe = 0 ;
  AC_TABLE tbl = 0, tbl2 = 0 ;
  int a1, a2, ir, jr, nok = 0, nok2 ;
  BOOL baddy ;
  KEY k1, k2 ;

  k1 = ac_obj_key (mrna) ;
  tbl = ac_tag_table (mrna, "Probe_hit", h) ;

  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      probe = ac_table_obj (tbl, ir, 0, h) ;
      /* count its errors against the NT mrna */
      a1 = a2 = nok = 0 ;
      for (jr = 0 ; tbl && jr < tbl->rows ; jr++)
	{
	  k2 = ac_table_key (tbl, jr, 0, 0) ;
	  if (k1 == k2)
	    {
	      a1 = ac_table_int (tbl, jr, 1, 0) ;
	      a2 = ac_table_int (tbl, jr, 2, 0) ;
	      nok = a2 - a1 + ac_table_int (tbl, jr, 3, 0) ;
	      break ;
	    }
	}
      baddy = FALSE ;
      if (a1 > a2)
	baddy = TRUE ;
      else
	tbl2 = ac_tag_table (probe, "mrna_hit", h) ;
      /* count the error in the other mrnas */
      for (jr = 0 ; ! baddy && tbl2 && jr < tbl2->rows ; jr++)
	{
	  mrna2 = ac_table_obj (tbl2, jr, 0, h) ;
	  if (! strstr (ac_tag_printable (mrna2, "IntMap", "NT"), "NT"))
	    { /* found a good non NT mrna */
	      a1 = ac_table_int (tbl2, jr, 1, 0) ;
	      a2 = ac_table_int (tbl2, jr, 2, 0) ;
	      nok2 = a2 - a1 + ac_table_int (tbl2, jr, 3, 0) ;
	      if (a1 < a2 && nok2 >= nok)
		baddy = TRUE ;
	    }
	}
      ac_free (tbl2) ;
      if (baddy)
	{
	  freeOutf ("Probe %s\n", freeprotect (ac_name(probe))) ;
	  freeOutf ("-D mrna_hit %s\n", freeprotect (ac_name(mrna))) ;
	  freeOutf ("-D mrna_exact_hit %s\n", freeprotect (ac_name(mrna))) ;
	  freeOutf ("\n") ;
	}
    }
  ac_free (h) ;
} /* flagRemoveOneNtHits */

/*************************************************************************************/
/* remove the hits to NT contigs if there exist a hit to a non NT contig */
static void flagRemoveMultiTg2gene (AC_OBJ tg)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = ac_tag_table (tg, "Gene", h) ;
  int ir, bestIr ;
  char *cp, buffer[1000] ;

  strncpy (buffer, ac_name(tg), 999) ;
  cp = buffer + strlen(buffer) - 1 ;
  while (cp > buffer && *cp != '.') cp-- ;
  if (*cp == '.') *cp = 0 ;
  for (ir = 0, bestIr = -1 ; tbl && bestIr < 0 && ir < tbl->rows ; ir++)
    {
      if (! strcmp (buffer, ac_table_printable (tbl, ir, 0, "")))
	bestIr = ir ; 
    }
  if (bestIr == -1) bestIr = 0 ;
  freeOutf ("Transcribed_gene %s\n", freeprotect (ac_name(tg))) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    if (ir != bestIr)
      freeOutf ("-D Gene %s\n", freeprotect (ac_table_printable (tbl, ir, 0, ""))) ;

  freeOutf ("\n") ;
  ac_free (h) ;
} /* flagRemoveMultiTg2gene */

/*************************************************************************************/
/* remove the hits to NT contigs if there exist a hit to a non NT contig */
static void flagReverseNmAntiHit (AC_OBJ probe)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = ac_tag_table (probe, "nm_hit", h) ;
  int ir, a1, a2, nr, ns ;
  const char *cp ;

  for (ir = 0 ; ir < tbl->rows ; ir++)
    {
      a1 = ac_table_int (tbl, ir, 1, 0) ;
      a2 = ac_table_int (tbl, ir, 2, 0) ;
      nr = ac_table_int (tbl, ir, 3, 0) ;
      ns = ac_table_int (tbl, ir, 4, 0) ;
      if (a1 < a2) /* this is not an anti hit */
	continue ;
      cp = ac_table_printable (tbl, ir, 0, 0) ;
      freeOutf ("Probe %s\n", freeprotect (ac_name(probe))) ;
      freeOutf ("-D NM_exact_hit %s\n", freeprotect (cp)) ;
      freeOutf ("-D NM_hit %s\n", freeprotect (cp)) ;
      freeOutf ("NM_anti_hit %s %d %d %d %d\n\n", freeprotect (cp), a1, a2, nr, ns) ;
      freeOutf ("Sequence %s\n",  freeprotect (cp)) ;
      freeOutf ("probe_NM_anti_hit %s %d %d %d %d\n\n", freeprotect (ac_name(probe)), a1, a2, nr, ns) ;
    }

  freeOutf ("\n") ;
  ac_free (h) ;
} /* flagReverseNmAntiHit */

/*************************************************************************************/
/* remove the hits to NT contigs if there exist a hit to a non NT contig */
static void flagRemoveNtHit (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ mrna = 0, tg = 0, probe = 0 ;
  AC_ITER iter ;
  
  if (1)
    {
      iter = ac_dbquery_iter (mm->db, "Find mrna intmap = *nt_*", h) ;
      while (ac_free (mrna), mrna = ac_iter_obj (iter))
	flagRemoveOneNtHits (mrna) ;
    }
  if (1)
    {
      iter = ac_dbquery_iter (mm->db, "Find transcribed_gene COUNT gene > 1", h) ;
      while (ac_free (tg), tg = ac_iter_obj (iter))
	flagRemoveMultiTg2gene (tg) ;
    }
  if (1)
    {
      iter = ac_dbquery_iter (mm->db, "Find probe nm_hit", h) ;
      while (ac_free (probe), probe = ac_iter_obj (iter))
	flagReverseNmAntiHit (probe) ;
    }
  ac_free (h) ;
} /* flagRemoveNtHits */

/*************************************************************************************/
/*************************************************************************************/

static void exportServerDataOneMrna (AC_OBJ mrna)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ tg = 0, gene = 0, product = 0 ;
  AC_TABLE tbl = 0 ;
  int ir, a1, a2, x1, x2 ;
  const char *cp ;

  freeOutf ("\n\nmRNA %s\n", freeprotect (ac_name(mrna))) ;
  tbl = ac_tag_table (mrna, "Splicing", h) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      a1 = ac_table_int (tbl, ir, 0, 0) ;
      a2 = ac_table_int (tbl, ir, 1, 0) ;
      x1 = ac_table_int (tbl, ir, 2, 0) ;
      x2 = ac_table_int (tbl, ir, 3, 0) ; 
      cp = ac_table_printable (tbl, ir, 4, "") ;
      freeOutf ("Splicing %d %d %d %d %s\n"
		, a1, a2, x1, x2, cp) ;
    }
  tbl = ac_tag_table (mrna, "Constructed_from", h) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      cp = ac_table_printable (tbl, ir, 2, "") ;
      if (!strncmp ("NM_", cp, 3))
	{
	  a1 = ac_table_int (tbl, ir, 0, 0) ;
	  a2 = ac_table_int (tbl, ir, 1, 0) ;
	  x1 = ac_table_int (tbl, ir, 3, 0) ;
	  x2 = ac_table_int (tbl, ir, 4, 0) ; 
	  cp = ac_table_printable (tbl, ir, 2, "") ;
	  freeOutf ("Constructed_from %d %d %s %d %d\n"
		  , a1, a2, cp, x1, x2) ;
	}
    }
  tbl = ac_tag_table (mrna, "IntMap", h) ;
  if (tbl && tbl->cols > 2)
    {
      ir = 0 ;
      a1 = ac_table_int (tbl, ir, 1, 0) ;
      a2 = ac_table_int (tbl, ir, 2, 0) ;
      cp = ac_table_printable (tbl, ir, 0, "") ;
      freeOutf ("IntMap %s %d %d\n", cp, a1, a2) ;
    }
  tbl = ac_tag_table (mrna, "Product", h) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      a1 = ac_table_int (tbl, ir, 1, 0) ;
      a2 = ac_table_int (tbl, ir, 2, 0) ;
      product = ac_table_obj (tbl, ir, 0, h) ;
      if (ac_has_tag (product, "best_product") &&
	  ac_has_tag (product, "down_stop"))
	{
	  freeOutf ("Product %s %d %d\n\n"
		    ,  freeprotect (ac_name(product))
		    , a1, a2
		    ) ;
	  freeOutf ("Product %s\nbest_product\ndown_stop\n"
		    ,  freeprotect (ac_name(product))
		    ) ;
	  break ;
	}
    }
  tg = ac_tag_obj (mrna, "From_gene", h) ;
  if (tg)
    {
      gene = ac_tag_obj (tg, "gene", h) ;
      freeOutf ("\nTranscribed_gene %s\n", freeprotect (ac_name(tg))) ;
      if (gene)
	freeOutf ("Gene %s\n", freeprotect (ac_name(gene))) ;
      tbl = ac_tag_table (tg, "mrna", h) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  a1 = ac_table_int (tbl, ir, 1, 0) ;
	  a2 = ac_table_int (tbl, ir, 2, 0) ;
	  cp = ac_table_printable (tbl, ir, 0, "") ;
	  if (cp && !strcmp (cp, ac_name(mrna)))
	    freeOutf ("mrna %s %d %d\n"
		    , freeprotect (cp), a1, a2) ;
	}
      tbl = ac_tag_table (tg, "IntMap", h) ;
      if (tbl && tbl->cols > 2)
	{
	  ir = 0 ;
	  a1 = ac_table_int (tbl, ir, 1, 0) ;
	  a2 = ac_table_int (tbl, ir, 2, 0) ;
	  cp = ac_table_printable (tbl, ir, 0, "") ;
	  freeOutf ("IntMap %s %d %d\n", cp, a1, a2) ;
	}
      freeOutf ("\n") ;
    }
  if (gene && (tbl = ac_tag_table (gene, "GeneId", h)))
    {
      freeOutf ("\nGene %s\n", freeprotect (ac_name(gene))) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  cp = ac_table_printable (tbl, ir, 0, "") ;
	  if (cp)
	    freeOutf ("GeneId %s\n", cp) ;
	}
      freeOutf ("\n") ;
    }
  freeOutf ("\n") ;
  ac_free (h) ;
} /* exportServerDataOneMrna */

/*************************************************************************************/

static void exportServerData (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ mrna = 0 ;
  AC_ITER iter = 0 ;
  AC_KEYSET nms = 0, mrnas = 0, nms2mrnas = 0 ;
 
  nms = ac_read_keyset (mm->db, mm->nm_file, h) ;
  if (nms)
    nms2mrnas = ac_ksquery_keyset (nms, ">in_mrna", h) ;
  mrnas = ac_read_keyset (mm->db, mm->mrna_file, h) ;

  if (mrnas)
    {
      if (nms2mrnas)
	ac_keyset_or (mrnas, nms2mrnas) ;
      iter = ac_keyset_iter (mrnas, TRUE, h) ;
      while (ac_free (mrna), (mrna = ac_iter_obj (iter)))
	exportServerDataOneMrna (mrna) ;
    }
  ac_free (h) ;
} /* exportServerData */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static void removeNegativeHitsOneProbe_Set (MM *mm, AC_OBJ Probe_set)
{
  AC_HANDLE h = ac_new_handle () ; 
  AC_OBJ probe = 0 ;
  AC_TABLE hits, probes ;
  int ir, jr, p1, p2 ;
  BOOL ok, nok ;
  
  ok = nok = FALSE ;
  probes = ac_tag_table (Probe_set, "Probe", 0) ;
  for (ir = 0 ; probes && ir < probes->rows ; ir++)
    {
      probe = ac_table_obj (probes, ir, 0, h) ;
      freeOutf ("Probe %s\n", freeprotect (ac_name(probe))) ;
      hits = ac_tag_table (probe, "nm_hit", h) ;
      for (jr = 0 ; hits && jr < hits->rows ; jr++)
	{
	  p1 = ac_table_int (hits, jr, 1, 0) ;
	  p2 = ac_table_int (hits, jr, 2, 0) ;
	  if (p1 < p2) 
	    ok = TRUE ; /* at least one positive hit */
	  else
	    {
	      if (1) ok = TRUE ; /* kill all negatives */
	      nok = TRUE ;
	      freeOutf ("-D nm_hit %s %d %d\n"
			  , freeprotect (ac_table_printable (hits, jr, 0, "toto"))
			  , p1, p2) ;
	      freeOutf ("-D nm_exact_hit %s %d %d\n"
			  , freeprotect (ac_table_printable (hits, jr, 0, "toto"))
			  , p1, p2) ;
	    }
	}
      freeOutf ("\n") ;
    }

  ac_free (h) ;
} /* removeNegativeHitsOneProbe_Set */

/*************************************************************************************/

static void removeNegativeHits (MM *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER iter = ac_dbquery_iter (mm->db, "Find Probe_set", h) ; 
  AC_OBJ Probe_set = 0 ;

  while (ac_free (Probe_set), (Probe_set = ac_iter_obj (iter)))
    removeNegativeHitsOneProbe_Set (mm, Probe_set) ;

  ac_free (h) ;
} /* removeNegativeHits */

/*************************************************************************************/
/*************************************************************************************/
/* -exportG2Nm2hits : export 3 column table gene/NM(multi)/non-ambiguous-probe(multi) */
static void exportG2Nm2hits (MM *mm)
{
  AC_HANDLE h1, h = ac_new_handle () ;
  AC_ITER probes, iter = ac_dbquery_iter (mm->db, "Find locus ; COUNT {>sequence_gb ; probe_nm_hit} > 0", h) ; 
  AC_OBJ locus = 0, probe = 0 ;
  AC_TABLE nms ;
  int ir ;

  while (ac_free (locus), (locus = ac_iter_obj (iter)))
    {
      h1 = ac_new_handle () ;
      printf ("%s\t", ac_name(locus)) ;
      nms = ac_tag_table (locus, "sequence_gb", h1) ;
      for (ir = 0 ; ir < nms->rows ; ir++)
	printf ("%s%s"
		, ir ? ";" : ""
		, ac_table_printable (nms, ir, 0, "")
		) ;
      printf ("\t") ;

      probes = ac_objquery_iter (locus, ">sequence_gb; >probe_nm_hit", h1) ;
      ir = -1 ;
      while (ac_free (probe), ir++, (probe = ac_iter_obj (probes)))
	printf ("%s%s"
		, ir ? ";" : ""
		, ac_name (probe)
		) ;
      printf ("\n") ;
 
      ac_free (h1) ;
    }

  ac_free (h) ;
} /* exportG2Nm2hits */

/*************************************************************************************/
/*************************************************************************************/
/* -exportMaqcMapping export 2006_02_11 for MAQC */
typedef struct  mmxStruct { KEY nm, gene, mrna, probe, geneid, geneid2, chrom ; int a1, a2, gg1, gg2, nr, ns, ln, probe_class ;} MMX ;
static int mmxOrder (const void *a, const void *b)
{
  const MMX *ma = (const MMX *)a, *mb = (const MMX *)b ;
  int nn = 0 ;

  if (ma->gene && mb->gene)
    nn =  (int) (ma->gene) -  (int) (mb->gene) ;
  else if (ma->gene)
    nn = -1 ; /* aceview gene first */
  else if (mb->gene)
    nn =  1 ; /* aceview gene first */
  else if (ma->geneid && mb->geneid)
    {
      nn =  (int) (ma->geneid) -  (int) (mb->geneid) ;
      if (! nn && ma->nm && mb->nm)
	nn =  (int) (ma->nm) -  (int) (mb->nm) ;
    }
  if (nn)
    return nn ;
  /* now we are in the same gene or at least the same nm */
  nn = ma->a1 - mb->a1 ;
  if (nn) return nn ;
  nn = ma->a2 - mb->a2 ;
  if (nn) return nn ;
  nn = (int) (ma->probe) - (int)(mb->probe) ;
  return nn ;
} /* mmxOrder */

/*************************************************************************************/

static KEY exportMaqcMappingGeneId (AC_OBJ Gene, int mmxA1, KEY chrom, int tg1, int tg2)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ GeneId = 0, Pg = 0 ;
  AC_TABLE tbl2 = 0, tbl3 = 0, tbl4 = 0 ;
  int jr, kr, pg1, pg2, ig = 0, mid ;
  Array pgs = arrayHandleCreate (12, MMX, h) ;
  MMX *mmx2 ;
  KEY geneid = 0 ;

  tbl2 = ac_tag_table (Gene, "GeneId", h) ;
  if (chrom && tbl2 && tbl2->rows > 1)
    {
      for (jr = 0 ; tbl2 && jr < tbl2->rows && jr < 10 ; jr++)
	{ 
	  GeneId = ac_table_obj (tbl2, jr, 0, h) ;
	  if (GeneId)
	    { 
	      geneid = ac_obj_key (GeneId) ;
	      tbl3 = ac_tag_table (GeneId, "Predicted_gene", h) ;
	      for (kr = 0 ; tbl3 && kr < tbl3->rows && kr < 10 ; kr++)
		{
		  Pg = ac_table_obj (tbl3, kr, 0, h) ;
		  tbl4 = ac_tag_table (Pg, "IntMap", h) ;
		  if (chrom != ac_table_key (tbl4, 0, 0, 0))
		    continue ;
		  mmx2 = arrayp (pgs, ig++, MMX) ;
		  mmx2->geneid2 = ac_obj_key (GeneId) ;
		  pg1 = ac_table_int (tbl4, 0, 1, 0) ;
		  pg2 = ac_table_int (tbl4, 0, 2, 0) ;
		  if (tg1 < tg2)
		    { 
		      mmx2->a1 = pg1 - tg1 ;
		      mmx2->a2 = pg2 - tg1 ;
		    }
		  else
		    { 
		      mmx2->a1 = - pg1 + tg1 ;
		      mmx2->a2 = - pg2 + tg1 ;
		    }
		}
	    }
	}
      arraySort (pgs, mmxOrder) ;
      arrayCompress (pgs) ;
      for (ig = 0 ; ig < arrayMax (pgs) - 1 ;ig++)
	{
	  mmx2 = arrayp (pgs, ig, MMX) ;
	  if ((mmx2 + 1)->a2 > mmx2->a2)
	    { 
	      mid = ((mmx2 + 1)->a1 + mmx2->a2)/2 ;
	      (mmx2 + 1)->a1 = mid + 1 ; 
	      mmx2->a2 = mid ;
	    }
	  else
	    *(mmx2 + 1) = *mmx2 ;
	}
      arrayCompress (pgs) ;
      for (ig = 0 ; ig < arrayMax (pgs) ;ig++)
	{
	  mmx2 = arrayp (pgs, ig, MMX) ;
	  if (mmxA1 >= mmx2->a1 && mmxA1 < mmx2->a2)
	    {  geneid = mmx2->geneid2 ; break ; }
	}
    }
  else if (tbl2)
    geneid = ac_table_key (tbl2, 0, 0, 0) ;

  ac_free (h) ;
  return geneid ;
}  /* exportMaqcMappingGeneId */

/*************************************************************************************/
/* register the coordinates and see if they are the same */
static BOOL exportMaqcMapping2Mrna (AC_OBJ Probe, AC_OBJ Mrna, int p1, int p2, int *a1p, int *a2p)
{
  int a1, a2, x1, x2, q1, q2, jr ;
  BOOL ok = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE exons = ac_tag_table (Mrna, "Splicing", h) ;

  exons = ac_tag_table (Mrna, "Splicing", h) ;
  for (q1 = q2 = -100, jr = 0 ; jr < exons->rows ; jr++)
    {  
      a1 = ac_table_int (exons, jr, 0, 0) ;
      a2 = ac_table_int (exons, jr, 1, 0) ;
      x1 = ac_table_int (exons, jr, 2, 0) ;
      x2 = ac_table_int (exons, jr, 3, 0) ;
      if (p1 >= x1 && p1 <= x2)
	{ q1 = a1 + p1 - x1 ; break ; }
    }
  for ( ; jr < exons->rows ; jr++)
    {  
      a1 = ac_table_int (exons, jr, 0, 0) ;
      a2 = ac_table_int (exons, jr, 1, 0) ;
      x1 = ac_table_int (exons, jr, 2, 0) ;
      x2 = ac_table_int (exons, jr, 3, 0) ;
      if (p2 >= x1 && p2 <= x2)
	{ q2 = a1 + p2 - x1 ; break ; }
    }
  
  if (q1 > 0 && q2 > 0)
    { 
      *a1p = q1 ; *a2p = q2 ; ok = TRUE ;
    }
  ac_free (h) ;
  return ok ;
}

/*************************************************************************************/
/* export the gene specific probes */
static void exportMaqcMapping (MM *mm)
{
  AC_HANDLE h2 = 0, h = ac_new_handle () ;
  AC_ITER genes ;
  AC_TABLE tbl = 0, probes, mrnas, imap, geneIds ; 
  AC_OBJ Probe = 0, Nm = 0, Tg = 0, Gene = 0, Mrna = 0, Gf = 0 ;
  Array aa = 0, gfs = 0 ;
  int i, ir, ip, im, ix = 0, m1, iexact, a1, a2, p1, p2, nr, ns, gg1, gg2, gf1, gf2 ;
  MMX *mmx ;
  KEY nm, mrna, geneid, tg, gId, chrom ;
  KEYSET ks = arrayHandleCreate (64, KEY, h) ;
  BOOL ok = FALSE ;  /* used to export a single position per probe */
  char *cp, buff[1000] ;
  char *myHit[] = {"mrna_hit" , "mrna_hit" , "nm_exact_hit", "nm_hit", 0} ;

  genes =  ac_dbquery_iter (mm->db
			    , messprintf("Find gene maqc_probe ; IS %s"
					 , mm->template ? mm->template : "*")
			    , h) ;
  aa = arrayHandleCreate (1000, MMX, h) ; 
  gfs = arrayHandleCreate (10, MMX, h) ;
  while (ac_free (Gene), (Gene = ac_iter_obj (genes)))
    {
      ac_free (h2) ;
      h2 = ac_new_handle () ; 
      aa = arrayReCreate (aa, 1000, MMX) ;
      gfs = arrayReCreate (gfs, 1000, MMX) ;
      probes =  ac_tag_table (Gene, "maqc_probe", h2) ; 
      Tg = ac_tag_obj (Gene, "Transcribed_gene", h2) ; 
      tg = ac_obj_key (Tg) ;
      mrnas = ac_tag_table (Tg, "mRNA", h2) ;
      geneid = ac_tag_key (Gene, "GeneId", 0) ; 
      gId = 0 ;
      geneIds = ac_tag_table (Gene, "GeneId", h2) ; 
      imap = ac_tag_table (Tg, "IntMap", h2) ;
      chrom = ac_table_key (imap, 0, 0, 0) ;
      gg1 = ac_table_int (imap, 0, 1, 0) ;
      gg2 = ac_table_int (imap, 0, 2, 0) ;
      if (geneIds && geneIds->rows > 1)
	{
	  tbl = ac_tag_table (Tg, "Matching_genefinder_gene", h2) ;

	  for (ir = 0 ; tbl && ir < tbl->rows && !ok  ; ir++) 
	    {
	      Gf = ac_table_obj (tbl, ir, 0, h2) ;
	      imap = ac_tag_table (Gf, "IntMap", h2) ;
	      gf1 = ac_table_int (imap, 0, 1, 0) ;
	      gf2 = ac_table_int (imap, 0, 2, 0) ;
	      if (gg1 < gg2) { gf1 = gf1 - gg1 ; gf2 = gf2 - gg1 ; }
	      else { gf1 = gg1 - gf1 ; gf2 = gg1 - gf2 ; }
	      mmx = arrayp (gfs, ir, MMX) ;
	      mmx->a1 = gf1 ; mmx->a2 = gf2 ; 
	      mmx->probe = ac_tag_key (Gf, "GeneId_pg", 0) ; /* to please mmxOrder */
	    }
	  arraySort (gfs, mmxOrder) ;
	  arrayCompress (aa) ;
	  for (ir = 0 ; ir < arrayMax(gfs) ; ir++)
	    {
	      mmx = arrayp (gfs, ir, MMX) ;
	      mmx->geneid = mmx->probe ; mmx->probe = 0 ;
	      if (ir > 0)
		{
		  a1 = ((mmx-1)->a2 + mmx->a1)/2 ;
		  (mmx-1)->a2 = mmx->a1 = a1 ; 
		}
	    }
	}
      for (ip = 0 ; probes && ip < probes->rows ; ip++)
	{
	  Probe = ac_table_obj (probes, ip, 0, h2) ;
	  if (!strncmp (ac_name(Probe), "MLT_", 4))
	    continue ;
	  if (!ac_has_tag (Probe, "Probe_class"))
	    continue ;
	  ks = keySetReCreate (ks) ;
	  ok = FALSE ;
	  for (iexact = 0 ; !ok && iexact < 4 ; iexact++)
	    {
	      tbl = ac_tag_table (Probe, myHit[iexact], h2) ;
	      for (ir = 0 ; !ok && tbl && ir < tbl->rows && !ok  ; ir++)
		{
		  nm = 0 ;
		  mrna = ac_table_key (tbl, ir, 0, 0) ;
		  p1 = ac_table_int (tbl, ir, 1, 0) ;
		  p2 = ac_table_int (tbl, ir, 2, 0) ;
		  nr = ac_table_int (tbl, ir, 3, 0) ; 
		  ns = ac_table_int (tbl, ir, 4, 0) ;
		  /* locate the mrna in the Tg */
		  if (iexact < 2) /* not nm case */
		    {
		      for (im = -1, i = 0 ; im == -1 && mrnas && i < mrnas->rows ; i++)
			if (mrna == ac_table_key (mrnas, i, 0, 0))
			  im = i ;
		      if (im == -1) 
			continue ; 
		      a1 = a2 = 10000000 ;
		      if (iexact == 0)
			{
			  m1 = ac_table_int (mrnas, im, 1, 0) ; /* mRNA start in Tg */
			  a1 = ac_table_int (tbl, ir, 5, 0) ;
			  a2 = ac_table_int (tbl, ir, 6, 0) ;

			  if (!a2 && p2 > 0) 
			    {
			      Mrna = ac_table_obj (tbl, ir, 0, h2) ;
			      if (exportMaqcMapping2Mrna (Probe, Mrna, p1, p2, &a1, &a2))
				{
				  freeOutf ("ZZZZ\tProbe %s", freeprotect (ac_name (Probe))) ;
				  freeOutf ("\tmrna_hit %s %d %d %d %d %d %d\n"
					    , freeprotect (ac_name (Mrna))
					    , p1, p2, nr, ns, a1, a2
					    ) ;
				}
			    }
			  if (!a2)
			    continue ;
			  a1 += m1 ; a2 += m1 ;
			}
		      else
			{
			  m1 = ac_table_int (mrnas, im, 1, 0) ; /* mRNA start in Tg */
			  a1 = m1 + p1 + 10000000 ;
			  a2 = m1 + p2 + 10000000 ;
			}
		      gId = 0 ;
		      if (arrayMax(gfs) > 1)
			for (ir = 0 ; !gId && ir < arrayMax(gfs) ; ir++)
			  {
			    mmx = arrayp (gfs, ir, MMX) ;
			    if (a1 < mmx->a2)
			      gId = mmx->geneid ;
			  }
		    }		  
		  else  /* nm case */
		    { 
		      Nm = ac_table_obj (tbl, ir, 0, h2) ;
		      if (ac_tag_key (Nm, "From_gene", 0) != tg)
			continue ;
		      m1 = -1 ;
		      a1 = 10000000+ac_table_int (tbl, ir, 1, 0) ;
		      a2 = 10000000+ac_table_int (tbl, ir, 2, 0) ;
		      nm = ac_table_key (tbl, ir, 0, 0) ;
		      gId = ac_tag_key (Nm, "GeneId", 0) ;
		    }
		  if (!strncmp (ac_name (Probe), "TAQ_",4) ||  !strncmp (ac_name (Probe), "ABI_",4))
		    a1 = 20000000 ;
		  ok = TRUE ;
		  mmx = arrayp (aa, ix++, MMX) ;
		  mmx->gene = ac_obj_key (Gene) ; 
		  mmx->probe = ac_obj_key (Probe) ; 
		  mmx->ln = ac_tag_int (Probe, "Length", 0) ;
		  mmx->a1 = a1 ;
		  mmx->a2 = a2 ;
		  mmx->nr = nr ;
		  mmx->ns = ns ;
		  mmx->probe_class = ac_tag_int (Probe, "Probe_class", 0) ;
		  mmx->geneid = gId ? gId : geneid ; 
		  mmx->nm = nm ;
		  mmx->chrom = chrom ;
		  mmx->gg1 = gg1 ;
		  mmx->gg2 = gg2 ;
		}
	    }
	}
    

      /* we are now ready to export this gene */
      arraySort (aa, mmxOrder) ;
      arrayCompress (aa) ;
      
      for (ix = 0, mmx = arrp (aa, 0, MMX) ; ix < arrayMax (aa) ; mmx++, ix++)
	{
	  /* gene or NM */
	  if (mmx->gene)
	    freeOutf ("%s", ac_key_name (mmx->gene)) ;
	  else if (mmx->nm)
	    freeOutf ("%s", ac_key_name (mmx->nm)) ;
	  else
	    continue ;
	  if (mmx->geneid)
	    {
	      freeOutf ("\t%s", ac_key_name (mmx->geneid)) ;
	    }
	  else
	    freeOutf ("\t-") ;
	  if (mmx->chrom)
	    freeOutf ("\t%s\t%d\t%d", ac_key_name (mmx->chrom), mmx->gg1, mmx->gg2) ;
	  else
	    freeOutf ("\t-\t-\t-") ;
	  if (mmx->a1 > 0)
	    freeOutf ("\t%012d", mmx->a1) ;
	  else
	    freeOutf ("\t%s", "ERROR missing coordinates") ;
	  if (mmx->probe)
	    {
	      strncpy (buff, ac_key_name (mmx->probe), 4) ;
	      buff[3] = 0 ;
	      freeOutf ("\t%s", buff) ;
	      if (0) /* probe set */
		{
		  if (!strcmp (buff, "AFX_"))
		    {
		      strcpy (buff, ac_key_name (mmx->probe)) ;
		      cp = strstr (buff, "_at") ;
		      if (cp)
			{
			  *cp = 0 ;
			  freeOutf ("\t%s", buff) ;
			}
		      else
			freeOutf ("\t%s", "ERROR") ;
		    }
		  else
		    freeOutf ("\t-") ;
		}
	      freeOutf ("\t%s", ac_key_name (mmx->probe) + 4) ; 
	      if (mmx->nr && mmx->ns == mmx->ln) mmx->ns -= 1 ;
	      freeOutf ("\t%d\t%d\t%d\t%d", mmx->probe_class, mmx->ln, mmx->ns, mmx->nr) ;
	    }
	  else
	    {
	      freeOutf ("\t%s", "ERROR") ;
	      freeOutf ("\t%s", "ERROR") ;
	      freeOutf ("\t%s", "ERROR") ;
	      freeOutf ("\t%d\t%d\t%d\t%d", mmx->probe_class, mmx->ln, mmx->ns, mmx->nr) ;
	    }
	  freeOutf ("\n") ;		
	}
    } 
  ac_free (h) ;
} /* exportMaqcMapping */

/*************************************************************************************/
/*************************************************************************************/

static void exportMaqcMappingPerProbe (MM *mm)
{
  AC_HANDLE h2 = 0, h = ac_new_handle () ;
  AC_ITER probes ;
  AC_TABLE tbl = 0, tbl2 = 0 ;
  AC_OBJ Probe = 0, Nm = 0, Tg = 0, Gene = 0, Mrna = 0 ;
  int ir, p1, p2, m1, m2, tg1, tg2, a1, a2, probe_class, nr, ns, ln, pcl ;
  KEY chrom = 0 ;
  KEY geneid = 0, gene = 0 ;
  BOOL ok = FALSE ;  /* used to export a single position per probe */
  char buff[1000] ;
  char *nmTag[] = { "NM_hit",  "NM_hit", "NM_hit",  "NM_hit",  "NM_anti_hit",  "NM_anti_hit", "" , "" , "" , "" , "" , "" , "" , "" , "" } ;
  char *mrnaTag[] = { "mRNA_hit",  "mRNA_hit", "mRNA_hit",  "mRNA_hit",  "mRNA_anti_hit",  "mRNA_anti_hit", "" , "" , "" , "" , "" , "" , "" , "" , "" } ;

  freeOutf ("# Array\tProbe\tGene\tGeneId\tTranscript\tQuality\tp1\tp2\tprobe-length\tlongest-match\tn-error\n") ;

  probes = ac_dbquery_iter (mm->db, 
		  messprintf ("Find probe probe_class && length && ! IS MLT* ; IS %s"
			      , mm->template ? mm->template : "*")
			    , h) ;
  while (ac_free (Probe), (Probe = ac_iter_obj (probes)))
    {
      ac_free (h2) ;
      h2 = ac_new_handle () ;
      
      strncpy (buff, ac_name (Probe), 4) ;
      buff[3] = 0 ;
      
      /*** probe_class ***/
      pcl = probe_class = ac_tag_int (Probe, "Probe_class", 0) ;
      probe_class =  probe_class%100 ;
      probe_class /= 10 ;
      tbl = ac_tag_table (Probe, nmTag[probe_class], h2) ;
      for (ir = 0, Gene = Nm = 0  ; tbl && ir < tbl->rows && !ok ; ir++, ac_free (Nm), ac_free (Gene))
	{
	  Nm = ac_table_obj (tbl, ir, 0, h2) ;
	  p1 = ac_table_int (tbl, ir, 1, 0) ;	    
	  p2 = ac_table_int (tbl, ir, 2, 0) ;
	  nr = ac_table_int (tbl, ir, 3, 0) ;
	  ns = ac_table_int (tbl, ir, 4, 0) ;
	  ln = ac_tag_int (Probe, "Length", 0) ;
	  if (ns == ln && nr > 0) ns-- ;
	  geneid = ac_tag_key (Nm, "GeneId", 0) ; 
	  gene = ac_tag_key (Nm, "From_gene", 0) ;
	  /* export tne NM hits */
	  freeOutf ("%s", buff) ;
	  freeOutf ("\t%s", ac_name (Probe) + 4) ;
	  if (Gene)
	    freeOutf ("\t%s", ac_key_name (gene)) ;
	  else
	    freeOutf ("\t-") ;
	  if (geneid)
	    freeOutf ("\t%s", ac_key_name (geneid)) ;
	  else
	    freeOutf ("\t-") ;
	  freeOutf ("\t%s", ac_name (Nm)) ;
	  if (1 ||  strcmp (buff, "ILM"))
	    p1 = p2 = 0 ; 
	  freeOutf ("\t%d", pcl) ;
	  ln = ac_tag_int (Probe, "Length", 0) ;
	  if (ns == ln && nr > 0) ns-- ;
	  freeOutf ("\t%d\t%d\t%d\t%d\t%d\n", p1, p2, ln, ns, nr) ;
	}

      tbl = ac_tag_table (Probe, mrnaTag[probe_class], h2) ;
      for (ir = 0, Gene = Mrna = Tg = 0 ; tbl && ir < tbl->rows && !ok ; 
	   ir++, ac_free (Gene), ac_free (Tg), ac_free (Mrna))
	{
	  Mrna = ac_table_obj (tbl, ir, 0, h2) ; 
	  p1 = a1 = ac_table_int (tbl, ir, 1, 0) ;	    
	  p2 = a2 = ac_table_int (tbl, ir, 2, 0) ;
	  nr = ac_table_int (tbl, ir, 3, 0) ;
	  ns = ac_table_int (tbl, ir, 4, 0) ;
	  a1 = ac_table_int (tbl, ir, 5, 0) ;	    
	  a2 = ac_table_int (tbl, ir, 6, 0) ;
	  ln = ac_tag_int (Probe, "Length", 0) ;
	  if (ns == ln && nr > 0) ns-- ;
	 
	  tbl2 = ac_tag_table (Mrna, "IntMap", h2) ;
	  m1 = ac_table_int (tbl2, 0, 1, -1) ;
	  m2 = ac_table_int (tbl2, 0, 2, 0) ;
	  Tg = ac_tag_obj (Mrna, "From_gene", h2) ;
	  ac_free (tbl2) ;
	  tbl2 = ac_tag_table (Tg, "IntMap", h2) ;
	  chrom = ac_table_key (tbl2, 0, 0, 0) ;
	  tg1 = ac_table_int (tbl2, 0, 1, -1) ;
	  tg2 = ac_table_int (tbl2, 0, 2, 0) ;
	  ac_free (tbl2) ;
	  if (tg1 < tg2)
	    {
	      a1 = m1 - tg1 + a1 ;
	      a2 = m1 - tg1 + a2 ;
	    }
	  else
	    {
	      a1 = - m1 + tg1 + a1 ;
	      a2 = - m1 + tg1 + a2 ;
	    }
	  Gene = ac_tag_obj (Tg, "Gene", h2) ;
	  if (Gene)
	    geneid = exportMaqcMappingGeneId (Gene, a1, chrom, tg1, tg2) ;

	  /* export tne aceview hits */
	  freeOutf ("%s", buff) ;
	  freeOutf ("\t%s", ac_name (Probe) + 4) ;
	  if (Gene)
	    freeOutf ("\t%s", ac_name(Gene)) ;
	  else
	    freeOutf ("\t-") ;
	  if (geneid)
	    freeOutf ("\t%s", ac_key_name(geneid)) ;
	  else
	    freeOutf ("\t-") ;
	  freeOutf ("\t%s", ac_name(Mrna)) ;
	  if (1 || strcmp(buff, "ILM"))
	    p1 = p2 = 0 ; 
	  freeOutf ("\t%d", pcl) ;
	  freeOutf ("\t%d\t%d\t%d\t%d\t%d\n", p1, p2, ln, ns, nr) ;
 	}
    }
  ac_free (h2) ; 
  ac_free (h) ;
} /* exportMaqcMappingPerProbe */

/*************************************************************************************/
/*************************************************************************************/
/*  -exportP2Nm : export 2 column table Probe/NM(multi)/Gene(multi) */
static void exportP2Nm (MM *mm)
{
  AC_HANDLE h1, h = ac_new_handle () ;
  AC_ITER genes, iter = ac_dbquery_iter (mm->db, "Find probe nm_hit", h) ; 
  AC_OBJ probe = 0, gene = 0 ;
  AC_TABLE hits ;
  int ir ;

  while (ac_free (probe), (probe = ac_iter_obj (iter)))
    {
      h1 = ac_new_handle () ;
      printf ("%s\t",ac_name(probe)) ;
      hits = ac_tag_table (probe, "nm_hit", h1) ;
      for (ir = 0 ; ir < hits->rows ; ir++)
	printf ("%s%s"
		, ir ? ";" : ""
		, ac_table_printable (hits, ir, 0, "")
		) ;
      printf ("\t") ;

      genes = ac_objquery_iter (probe, ">nm_hit; >Locus_gb", h1) ;
      if (genes)
	{
	  ir = -1 ;
	  while (ac_free (gene), ir++, (gene = ac_iter_obj (genes)))
	    printf ("%s%s"
		    , ir ? ";" : ""
		    , ac_name (gene)
		    ) ;
	}
      printf ("\n") ;
 
      ac_free (h1) ;
    }

  ac_free (h) ;
} /* exportP2Nm */

/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: maggie ACEDB [-o outfile] [-mapProbe] [-mrna] [-cds] [-intron] [-t template] [-k keyset]\n") ;
  fprintf (stderr, "// Example:  maggie locusid 384D8-2 \n") ;
  fprintf (stderr, "//   -mapProbe : export the alignments of the probes mapped on the genome/mRNAs\n") ;
  fprintf (stderr, "//   -analyse : reports something possibly oboslete\n") ;
  fprintf (stderr, "//   -hlsStatsNm : report the number NM hit by somany libraries\n") ;
  fprintf (stderr, "//   -hlsStatsAceView : report the number NM hit by somany libraries\n") ;
  fprintf (stderr, "//   -hitStatsNm : report the number of probes by category\n") ;
  fprintf (stderr, "//   -hitStatsAceView : report the number of probes by category\n") ;
  fprintf (stderr, "//   -hitSupportAceView: for each gene, probe hits + number of supporting ESTs\n") ;
  fprintf (stderr, "//   -exportG2Nm2hits : export 3 column table gene/NM(multi)/non-ambiguous-probe(multi)\n") ;
  fprintf (stderr, "//   -exportP2Nm : export 2 column table Probe/NM(multi)/Gene(multi)\n") ;
  fprintf (stderr, "//   -removeNegativeHits : remove negative hits if at least one positive hit exists\n") ;
  fprintf (stderr, "//   -removeNtHits : remove hits on NT contigs if a better one exists\n") ;
  fprintf (stderr, "//   -flagAmbiguousProbe : export a list of ambiguous probes\n") ;
  fprintf (stderr, "//   -flagAmbiguousProbeAceView : exports the probe type\n") ;
  fprintf (stderr, "//   -exportMaqcMapping : exports the probe mapping for MAQC\n") ;
  fprintf (stderr, "//   -exportMaqcMappingPerProbe : exports the probe mapping for MAQC\n") ;
  fprintf (stderr, "//   -cds : forget the UTR\n") ;
  fprintf (stderr, "//   -isMaggie : applies to h??Stats to select the set of libraries\n") ;
  fprintf (stderr, "//   -isNmDb : applies to flagAmbiguosProbe, to select the schema\n") ;
  fprintf (stderr, "//   -t template : a query on class tg\n") ;
  fprintf (stderr, "//   -k keyset : a keyset.list of class tg\n") ;
  fprintf (stderr, "//   -mrna : after choosing templete/keyset, move from class tg into classs mrna\n") ;
  fprintf (stderr, "//   -exportAceViewMrnaAsGenomic : export from public server all mRNAs as GM_ genomic + detailled data\n") ;
  fprintf (stderr, "//   -exportNmAsGenomic : export from public server all NM as genomic + detailled data\n") ;
  fprintf (stderr, "//   -nm_file ff1 -mrna_file ff2 : ff1 is a .list of nm, ff2 a .list of mrna, exports relevant .ace data about them\n") ;
  fprintf (stderr, "//   -out filename : redirects the output in that file, useful for batch jobs\n") ;

  exit (1) ;
}

/*************************************************************************************/

int main (int argc, const char **argv)
{
  FILE *f = 0 ;
  const char *s = "ok" ;
  const char *outfilename = 0 ;
  const char *dbName ;
  int outlevel = 0 ;
  MM mm ;
  
  memset (&mm, 0, sizeof (MM)) ;
  /* consume optional args */
  getCmdLineOption (&argc, argv, "-out", &outfilename ) ;
  mm.analyse = getCmdLineOption (&argc, argv, "-analyse", 0) ;
  mm.isMrna = getCmdLineOption (&argc, argv, "-mrna", 0) ;
  mm.isMapProbe = getCmdLineOption (&argc, argv, "-mapProbe", 0) ;
  mm.hlsStatsNm = getCmdLineOption (&argc, argv, "-hlsStatsNm", 0) ;
  mm.hlsStatsAceView = getCmdLineOption (&argc, argv, "-hlsStatsAceView", 0) ;
  mm.hitStatsNm = getCmdLineOption (&argc, argv, "-hitStatsNm", 0) ;
  mm.hitStatsAceView = getCmdLineOption (&argc, argv, "-hitStatsAceView", 0) ; 
  mm.hitSupportAceView = getCmdLineOption (&argc, argv, "-hitSupportAceView", 0) ; 
  mm.isMaggie = getCmdLineOption (&argc, argv, "-isMaggie", 0) ;
  mm.isNmDb = getCmdLineOption (&argc, argv, "-isNmDb", 0) ;
  mm.isCDS = getCmdLineOption (&argc, argv, "-cds", 0) ;
  mm.exportAceViewMrnaAsGenomic = getCmdLineOption (&argc, argv, "-exportAceViewMrnaAsGenomic", 0) ;
  mm.exportNmAsGenomic = getCmdLineOption (&argc, argv, "-exportNmAsGenomic", 0) ;
  mm.exportG2Nm2hits = getCmdLineOption (&argc, argv, "-exportG2Nm2hits", 0) ;
  mm.exportP2Nm = getCmdLineOption (&argc, argv, "-exportP2Nm", 0) ;
  mm.exportMaqcMapping = getCmdLineOption (&argc, argv, "-exportMaqcMapping", 0) ;
  mm.exportMaqcMappingPerProbe = getCmdLineOption (&argc, argv, "-exportMaqcMappingPerProbe", 0) ;
  mm.flagAmbiguousProbe = getCmdLineOption (&argc, argv, "-flagAmbiguousProbe", 0) ;
  mm.flagAmbiguousProbeAceView = getCmdLineOption (&argc, argv, "-flagAmbiguousProbeAceView", 0) ;
  mm.flagRemoveNtHit = getCmdLineOption (&argc, argv, "-removeNtHits", 0) ;
  mm.removeNegativeHits = getCmdLineOption (&argc, argv, "-removeNegativeHits", 0) ;
  getCmdLineOption (&argc, argv, "-t", &(mm.template)) ;
  getCmdLineOption (&argc, argv, "-k", &(mm.keysetName)) ;
  getCmdLineOption (&argc, argv, "-nm_file", &(mm.nm_file)) ;
  getCmdLineOption (&argc, argv, "-mrna_file", &(mm.mrna_file)) ;

  /* read absolute args */
  dbName = argc>=2 ? argv[1] : 0 ;
  if (!dbName) 
    usage () ;

  mm.db = ac_open_db (dbName, &s);
  if (!mm.db)
    messcrash ("Failed to open db %s, error %s", dbName, s) ;

  if (outfilename)
    {
      f = filopen (outfilename, 0, "w") ;
      if (f)
	outlevel = freeOutSetFile (f) ;
    }
  if (!outlevel)
    outlevel = freeOutSetFile (stdout) ;	

  if (mm.exportAceViewMrnaAsGenomic)
    magExportAceViewMrnaAsGenomic (&mm) ;
  else if (mm.exportNmAsGenomic)
    magExportNmAsGenomic (&mm) ;
  else if (mm.isMapProbe)
    magMapProbe (&mm) ;
  else if (mm.hitStatsNm)
    magHitStats (&mm, TRUE) ;
  else if (mm.hitStatsAceView)
    magHitStats (&mm, FALSE) ;
  else if (mm.hlsStatsNm)
    hlsStats (&mm, FALSE) ;
  else if (mm.hlsStatsAceView)
    hlsStats (&mm, TRUE) ;
  else if (mm.hitSupportAceView)
    hitSupportAceView (&mm, TRUE) ;
  else if (mm.flagAmbiguousProbe)
    flagAmbiguousProbe (&mm) ;
  else if (mm.flagAmbiguousProbeAceView)
    flagAmbiguousProbeAceView (&mm) ;
  else if (mm.flagRemoveNtHit)
    flagRemoveNtHit (&mm) ;
  else if (mm.removeNegativeHits)
    removeNegativeHits (&mm) ;
  else if (mm.exportMaqcMapping)
    exportMaqcMapping (&mm) ;
  else if (mm.exportMaqcMappingPerProbe)
    exportMaqcMappingPerProbe (&mm) ;
  else if (mm.exportG2Nm2hits)
    exportG2Nm2hits (&mm) ;
  else if (mm.exportP2Nm)
     exportP2Nm (&mm) ;
  else if (mm.analyse)
    magAnalyse (&mm) ;
  else if (mm.nm_file && mm.mrna_file)
    exportServerData (&mm) ;
  else
    usage () ;

  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;

  ac_db_close (mm.db) ;

  sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

