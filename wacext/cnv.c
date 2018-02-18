/*  File: cnv.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2015
 * -------------------------------------------------------------------
 * AceView is free software; you can redistribute it and/or
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
 * This file is part of the AceView project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov
 */

/*
 * Ananlysis of the CNV: copy number variations
 * from the exome data
 */
#define ARRAY_CHECK

#include "ac.h"
#include "channel.h"

typedef struct cnvStruct {
  AC_HANDLE h ;
  Array runs, genes, coverons, ends ;
  Array gene2coverons ;
  DICT *runDict ;
  DICT *chromDict ;
  DICT *geneDict ;
  DICT *coveronDict ;
  Stack info ;
  const char *project ;
  const char *title ;
  const char *chrom ;
  const char *outFileName ;
  const char *runsSortedFileName ;
  const char *runsAllFileName ;
  const char *runsControlFileName ;
  const char *gtitleFileName ;
  const char *geneSpongeFileName ;
  const char *coveronSpongeFileName ;
  const char *coverageFileName ;
  const char *u2uFileName ;
  const char *u2nuFileName ;
  CHAN *todo, *done ;
  int max_threads ;
  int nRuns[2] ;
  int nRunsY[2] ;
  float M, M1, M2 ;
  BOOL gzo ;
  BOOL unique, non_unique ; /* should we parse .u.cnv.sponge files OR .nu.cnv.sponge of .cnv.sponge */
  AC_DB db ;
} CNV ;

typedef struct getCoverageStruct { CNV *cnv ; int chrom ; int run ; } GCV ;

typedef struct runStruct {
  int type ;
  int run ;
  int sex ; /* 0 -> autosome, 1 : Y, 2 : X */
  float cumul, cumulX, cumulY ;
  Array cover, geneCover ;
  int machine, sample, tissue, title, runId, sortingTitle, sortingTitle2, otherTitle ;
  int histo[101] ;
} RC ;

typedef struct geneStruct {
  int coveron, gene ;
  int chrom, a1, a2, ln ;
  BOOL isDown ;
  int type ;
  float cumul[2], median[2], u2u, u2nu ;
  KEYSET coveron2genes ;
} GC ;

/*************************************************************************************/

typedef struct endStruct {
  int gene ;
  int chrom, a1 ;
  int type ;
} EC ;

/*************/

static int endOrder (const void *a, const void *b)
{
  const EC *up = (const EC *)a, *vp = (const EC *)b ;
  int z ;

  z = up->chrom - vp->chrom ; if (z) return z ;
  z = up->a1 - vp->a1 ; if (z) return z ;
  z = up->type - vp->type ; if (z) return z ;
  z = up->gene - vp->gene ; if (z) return z ;

  return 0 ;
} /* endOrder*/

/*************************************************************************************/
/* only runs in these lists will be considered */
static void cnvParseTitles (CNV *cnv)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (cnv->gtitleFileName, FALSE, h) ;
  int run ;
  DICT *dict = cnv->runDict ;
  RC *rc ;
  Stack s = cnv->info ;

  if (ai)
    {
      aceInSpecial (ai, "\n") ;
      while (aceInCard (ai))
	{
	  char cutter ;
	  const char *ccp = aceInWordCut (ai, "\t", &cutter) ;
	  if (! ccp) continue ;
	  ccp = ac_unprotect (ccp, h) ;
	  if ( ! dictFind (dict, ccp, &run))
	    continue ;
	  rc = arrayp (cnv->runs, run, RC) ;

	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->machine = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->sample = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  /* rc->system = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ; */
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->tissue = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->title = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  /* rc->system2 = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ; */
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->runId = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  /* rc->machine2 = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ; */
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->sortingTitle = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->sortingTitle = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->sortingTitle2 = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	  aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ; 
	  rc->otherTitle = stackMark (s) ; pushText (s, ac_unprotect (ccp, h)) ;
	}
    }
  ac_free (h) ;
  return ;
} /* cnvParseTitles */

/*************************************************************************************/
/* only runs in these lists will be considered */
static int cnvParseSpongeFile (CNV *cnv, int type)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, nE = arrayMax (cnv->ends) ;
  int gene ;
  const char *ff = 0 ;
  DICT *dict = 0 ;
  Array genes = 0 ;
  BOOL debug = 0 ;

  switch (type)
    {
    case 1: ff = cnv->geneSpongeFileName ; genes = cnv->genes ; dict = cnv->geneDict ; break ;
    case 2: ff = cnv->coveronSpongeFileName ; genes = cnv->coverons ; dict = cnv->coveronDict ; break ;
    }
  if (ff)
    {
      AC_HANDLE h1 = ac_new_handle () ;
      ACEIN ai= aceInCreate (ff, FALSE, h1) ;
      int gene, chrom, a1, a2 ;
      GC *gc ;
      char cutter ;
      const char *ccp ;
      BOOL isDown ;

      if (ai)
	{
	  aceInSpecial (ai, "\n") ;
	  while (aceInCard (ai))
	    if ((ccp = aceInWord (ai))) /* mrna name ignore */
	    {
	      aceInStep (ai, '\t') ; if (! (ccp = aceInWord (ai))) continue ; /* exon numbering, ignore  */
	      aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ;  /* chrom name */
	      if (cnv->chrom && strcasecmp (ccp, cnv->chrom))
		continue ;
	      dictAdd (cnv->chromDict, ccp, &chrom) ;
	      aceInStep (ai, '\t') ; if (! aceInInt (ai, &a1)) continue ;
	      aceInStep (ai, '\t') ; if (! aceInInt (ai, &a2)) continue ;
	      aceInStep (ai, '\t') ; ccp = aceInWordCut (ai, "\t", &cutter) ; if (! ccp) continue ;  /* gene name */
	      dictAdd (dict, ccp, &gene) ;
	      nn++ ;
	      gc = arrayp (genes, gene, GC) ;
	      if (type == 1) gc->gene = gene ;
	      else gc->coveron = gene ;
	      gc->type = type ;
	      gc->chrom = chrom ;
	      isDown = TRUE ;
	      if (a1 > a2)
		{ int a0 = a1 ; a1 = a2 ; a2 = a0 ; isDown = FALSE ; }
	      gc->a1 = a1 ; gc->a2 = a2 ; gc->isDown = isDown ;
	      /* ln = a2 - a1 + 1 ; */
	    }
	}
      ac_free (h1) ;
    }
  
  /* construct a sorted list of all ends */
  for (gene = 1 ; gene < arrayMax (genes) ; gene++)
    {
      GC *gc = arrayp (genes, gene, GC) ;
      EC *ec = arrayp (cnv->ends, nE++, EC) ;
      ec->chrom = gc->chrom ; ec->gene = gene ;
      ec->a1 = gc->a1 ; ec->type = 2 * type ; 
      
      ec = arrayp (cnv->ends, nE++, EC) ;
      ec->chrom = gc->chrom ; ec->gene = gene ;
      ec->a1 = gc->a2 ; ec->type = 2 * type +  1 ; 
    }

  arraySort (cnv->ends, endOrder) ;
  arrayCompress (cnv->ends) ;

  fprintf (stderr, "// parsed %d lines in sponge file type %d\n", nn, type) ;
  if (debug && type == 1)
    {
      int i ; EC *ec = arrayp (cnv->ends, 0, EC) ;
      nE = arrayMax (cnv->ends) ;
      for (i = 0 ; i < nE ; i++, ec++)
	{
	  EC *ec = arrayp (cnv->ends, i, EC) ;
	  fprintf (stdout, "%d\t%s\t%s\t%d\t%d\n", i
		   , ec->type < 4 ? dictName(cnv->geneDict, ec->gene) :  dictName(cnv->coveronDict, ec->gene) 
		   , dictName(cnv->chromDict, ec->chrom)
		   , ec->a1, ec->type
		   ) ;
	}
    }

  ac_free (h) ;

  return nn ;
} /* cnvParseSpongeFile */

 /*************************************************************************************/
/* only runs in these lists will be considered */
static int cnvParseRunLists (CNV *cnv, int type)
{
  int nn = 0 ;
  const char *ff = 0 ;

  switch (type)
    {
    case 0: ff = cnv->runsSortedFileName ; break ;
    case 1: ff = cnv->runsControlFileName ; break ;
    case 2: ff = cnv->runsAllFileName ; break ;
    }
  if (ff)
    {
      AC_HANDLE h1 = ac_new_handle () ;
      ACEIN ai= aceInCreate (ff, FALSE, h1) ;
      int run = 0 ;
      RC *rc ;
      const char *ccp ;
      DICT *dict = cnv->runDict ;

      if (ai)
	while (aceInCard (ai))
	  if ((ccp = aceInWord (ai)))
	    {
	      dictAdd (dict, ccp, &run) ;
	      rc = arrayp (cnv->runs, run, RC) ;
	      rc->run = run ;
	      rc->type = type ; nn++ ;
	    }
      ac_free (h1) ;
    } 
  fprintf (stderr, "// parsed %d runs of type %d\n", nn, type) ;
  return nn ;
} /* cnvParseRunLists */

/*************************************************************************************/
/*************************************************************************************/
/* measure the unicity coefficient */
static int cnvParseU2u (CNV *cnv)
{
  int pass ;
  int coveron ;
  DICT *dict = dict = cnv->coveronDict ;
  Array coverons = cnv->coverons ;
  GC *gc ;
  const char *ccp ;

  for (pass = 0 ; pass < 2 ; pass++)
    {
      int nn = 0 ;
      AC_HANDLE h = ac_new_handle () ;
      const char *fNam = pass == 0 ? cnv->u2uFileName :  cnv->u2nuFileName ;
      ACEIN ai = aceInCreate (fNam, FALSE, h) ;

      if (! ai)
	messcrash ("Cannot open the u2u file %s", fNam) ;

      aceInSpecial (ai, "\n") ;
      while (aceInCard (ai))
	{
	  int i, j, k ;
	  float z = 0 ;
	  
	  ccp = aceInWord (ai) ;
	  if (! ccp || *ccp == '#')
	    continue ;
	  aceInStep (ai, '\t') ; aceInWord (ai) ;
	  aceInStep (ai, '\t') ;
	  ccp = aceInWord (ai) ; /* coveron name */ 
	  if (! dictFind (dict, ccp, &coveron))
	    continue ;
	  gc = arrayp (coverons, coveron, GC) ;
	  for (j = k = 0, i = 4 ; i < 11 ; k++, i++)
	    if (aceInStep (ai, '\t')) 
	      { aceInWord (ai) ; j++ ; }
	  if (j < k) /* check that we jumped the correct number of columns */
	    continue ;
	  nn++ ;
	  if (aceInStep (ai, '\t') && aceInFloat (ai, &z)) 
	    switch (pass)
	      {
	      case 0 : gc->u2u = z ; break ;
	      case 1 : gc->u2nu = z ; break ;
	      }
	}
      fprintf (stderr, "// parsed %d lines in %s file %s\n"
	       , nn
	       , pass ? "u2nu" : "u2u"
	       , fNam
	       )  ;
      ac_free (h) ;
    }
  
  return 1 ;
} /* cnvParseU2u */

/*************************************************************************************/
/*************************************************************************************/

static void  cnvGeneCoveronRegister (CNV *cnv, int gene, int coveron)
{
  KEYSET g2c ;
  int i ;

  g2c = array (cnv->gene2coverons, gene, KEYSET) ;
  if (! g2c)
    g2c = array (cnv->gene2coverons, gene, KEYSET) = keySetCreate () ;
  for (i = 0 ; i < keySetMax (g2c) ; i++)
    if (keySet (g2c, i) == coveron)
      return ;
  keySet (g2c, i) = coveron ;
  return ;
}

/*************************************************************************************/
/* scan the ends file and associate to each coveron the genes it intersects */

static int cnvCoverons2Genes (CNV *cnv)
{
  EC *ec ;
  GC *gc ;
  int gene = 0, coveron = 0 ;
  int i, iMax = 0, ii, iiMax = arrayMax (cnv->ends),n,  nn = 0, coveronMax ;
  KEYSET genes = keySetCreate () ;
  BOOL ok ;

  coveronMax = arrayMax (cnv->coverons) ; 
  if (! coveronMax)
    return 0 ;

  for (ii = 0, ec = arrp (cnv->ends, ii, EC) ; ii < iiMax ; ii++, ec++)  
    {
      switch (ec->type)
	{
	case 2: /* gene start */
	  ok = FALSE ;
	  gene = ec->gene ;
	  for (i = 0 ; i < iMax && ! ok ; i++)
	    if (keySet (genes, i) == gene)
	      { ok = TRUE ; break ; }
	  for (i = 0 ; i < iMax && ! ok ; i++)
	    if (keySet (genes, i) == 0) 
	      { keySet (genes, i) = gene ; ok = TRUE ; break ; }
	  if (! ok) 
	    { keySet (genes, iMax) =  gene ; iMax++ ; }
	  if (coveron) /* add to the gene list of that coveron */
	    {
	      gc = arrayp (cnv->coverons, coveron, GC) ;
	      if (! gc->gene)  /* first gene of this coveron */
		{ gc->gene = gene ;  nn++ ; }
	      else if (! gc->coveron2genes)       /* second gene */
		{
		  gc->coveron2genes = arrayCreate (4, KEY) ; /* do not hook on a handle, there will be 10^4 of these */
		  keySet (gc->coveron2genes, 0) = gc->gene ;
		  keySet (gc->coveron2genes, 1) = gene ;
		}
	      else       /* append to the list */
		{
		  ok = FALSE ; 
		  for (i = 0 ; i < arrayMax (gc->coveron2genes) && ! ok ; i++)
		    if (keySet (gc->coveron2genes, i) == gene)
		      ok = TRUE ;
		  if (! ok)
		    keySet (gc->coveron2genes, i) = gene ;
		}
	    }
	  break ;
	case 3: /* gene end */
	  ok = FALSE ;
	  gene = ec->gene ;
	  for (i = 0 ; i < iMax && ! ok ; i++)
	    if (keySet (genes, i) == ec->gene) 
	      { keySet (genes, i) = 0 ; ok = TRUE ; break ; }
	  if (ok && i == iMax - 1) 
	    { iMax-- ; }
	  break ;
	case 4: /* coveron start */
	  coveron = ec->gene ;
	  gc = arrayp (cnv->coverons, coveron, GC) ;
	  n = 0 ; gene = 0 ; 
	  for (i = 0 ; i < iMax ; i++)
	    if (keySet (genes, i))
	      { gene = keySet (genes, i) ; n++ ; }
	  if (n) 
	    nn++ ;
	  if (n == 1) 
	    {
	      gc->gene = gene ; 
	    }
	  else if (n > 1)
	    {
	      gc->coveron2genes = arrayCreate (4+iMax, KEY) ; /* do not hook on a handle, there will be 10^4 of these */
	      for (i = n = 0 ; i < iMax  ; i++)
		if (keySet (genes, i))
		  keySet (gc->coveron2genes, n++) = keySet (genes, i) ;  
	      gc->gene = keySet (genes, 0) ;
	    }
	  break ;
	case 5: /* coveron end */
	  coveron = 0 ;
	  break ;
	}
    }
  keySetDestroy (genes) ;  

  /* create the back cross link */

  for (coveron = 1, gc = arrp (cnv->coverons, coveron, GC) ; coveron < coveronMax ; coveron++, gc++)
    if (gc->coveron && (gc->coveron2genes || gc->gene))
      {
	int i, gMax = 0 ;
	if (gc->coveron2genes && (gMax = keySetMax (gc->coveron2genes)))
	  for (i = 0 ; i < gMax ; i++)
	    {
	      gene = keySet (gc->coveron2genes, i) ;
	      cnvGeneCoveronRegister (cnv, gene, coveron) ;
	    }
	else if (gc->gene)
	  cnvGeneCoveronRegister (cnv, gc->gene, coveron) ;
      }
	  
  if (1)
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (cnv->outFileName, ".coveron2genes", FALSE, h) ;
      int n = 0, i, coveron, coveronMax = dictMax (cnv->coveronDict) + 1 ;

      aceOutDate (ao, "##", "Each gene is associated to each coveron it partially overlaps ina non stranded way") ;
      aceOut (ao, "# Ordinal\tCoveron internal number\tNumber of overlapping gene\tCoveron\tChomosome\tBegin of coveron\tEnd\tGenes\n" ) ;
      for (coveron = 1, gc = arrp (cnv->coverons, coveron, GC) ; coveron < coveronMax ; coveron++, gc++)
	if (gc->coveron && (gc->coveron2genes || gc->gene))
	  {
	    iMax = gc->coveron2genes ? keySetMax (gc->coveron2genes) : 1 ;
	    
	    aceOutf (ao, "%d\t%d\t%d\t%s\t%s\t%d\t%d\t"
		     , ++n
		     , coveron
		     , iMax
		     , dictName (cnv->coveronDict, gc->coveron)
		     , dictName (cnv->chromDict, gc->chrom)
		     , gc->a1, gc->a2
		     ) ;
	    if (gc->coveron2genes)
	      {	
		for (i = 0 ; i < iMax ; i++)
		  aceOutf (ao, "%s%s", i ? ", " : ""
			   , dictName (cnv->geneDict, keySet (gc->coveron2genes, i))
			   ) ;
	      }
	    else
	      aceOutf (ao, "%s", dictName (cnv->geneDict, gc->gene)) ;
	    
	    aceOut (ao, "\n") ;
	  }
      ac_free (h) ;
    }
  fprintf (stderr, "// Associated %d coverons to their genes\n", nn) ;
  return nn ;
} /* cnvCoverons2Genes */

/*************************************************************************************/
/*************************************************************************************/

static BOOL cnvDoParseOneCoverage (CNV *cnv, int run, int chrom, ACEIN ai, AC_HANDLE h)
{
  RC *rc = arrp (cnv->runs, run, RC) ;
  Array aa = rc->cover ;
  DICT *dict = cnv->coveronDict ;
  int i, nn = 0, nn1 = 0,  coveron = 0 ;
  char *cp ;
  float z ;
  BOOL debug = FALSE ;

  if (debug) fprintf (stderr, "# cnvDoParseOneCoverage %s", aceInFileName (ai)) ;
  while (aceInCard (ai)) /* closes ai->f and kills access to  aceInFileName (ai) */
    {
      aceInWord (ai) ; aceInStep (ai,'\t') ; aceInWord (ai) ;  /* jump col 1, 2 */
      aceInStep (ai,'\t') ; cp = aceInWord (ai) ;               /* coveron is in col 3 */
      if (!cp || ! dictFind (dict, cp, &coveron)) continue ;
      for (i = 4 ; i<= 10 ; i++) { aceInStep (ai,'\t') ; aceInWord (ai) ; }  /* jump to col 11 */
      aceInStep (ai,'\t') ; if (! aceInFloat (ai, &z)) continue ;
      array (aa, coveron, float) = z ;
      nn++ ; if (z > 0) nn1++ ;
    }
  if (debug) fprintf (stderr, " contains %d/%d positive values\n", nn1, nn) ;
  return nn > 0 ? TRUE : FALSE ;
} /* cnvDoParseOneCoverage  */

/*************************************************************************************/
#include <sys/file.h>
static void cnvParseOneCoverage (void *vp)
{
  CNV *cnv = (CNV*)vp ;
  ACEIN ai = 0 ;
  GCV gcv ;
  char *buffer[1024], *cp ;
  BOOL b = FALSE ;
  int nGet = 0, nPut = 0 ;

  memset (buffer, 0, sizeof(buffer)) ;
  while (channelMultiGet (cnv->todo, &gcv, 1, gcv))
    {
       AC_HANDLE h = ac_new_handle () ;
       nGet++ ;
       cp = hprintf (h, "%s/%s/%s.cnv.%ssponge"
		     , cnv->coverageFileName
		     , dictName (cnv->chromDict, gcv.chrom )
		     , dictName (cnv->runDict, gcv.run)
		     , cnv->unique ? "u." : (cnv->non_unique ?  "nu." : "")
		      ) ; 
       if (! access (cp, R_OK))
	 {
	   ai = aceInCreate (cp, FALSE, h) ;
	   b = ai ? cnvDoParseOneCoverage (cnv, gcv.run, gcv.chrom, ai, h) : FALSE ;
	 }
       else if (! strstr (cp, "|NT") && ! strstr (cp, "RAND")&& ! strstr (cp, "UN") && ! strstr (cp, "Un") && ! strstr (cp, "mito"))
	 fprintf (stderr, "# cannot open %s\n", cp) ;
       ac_free (h) ;
       nPut++ ;
      channelPut (cnv->done, &b, BOOL) ;
    }
  fprintf (stderr, "...cnvParseOneCoverage   nGet = %d  nPut = %d\n", nGet, nPut) ; 
  return ;
} /*  cnvParseOneCoverage */

/*************************************************************************************/

/* loop of the sponge files of all relevant runs and relevant chromosomes */
static int cnvParseCoverage (CNV *cnv)
{
  int k, nToDo = 0, nDone = 0, nDoneOk = 0 ;
  int chrom, run ;
  int cMax = arrayMax (cnv->coverons) ;
  int gMax = arrayMax (cnv->genes) ;
  GCV gcv ;

  cnv->todo = channelCreate (dictMax (cnv->runDict) * dictMax (cnv->chromDict), GCV, cnv->h) ;
  cnv->done = channelCreate (6 + 0 * dictMax (cnv->runDict) * dictMax (cnv->chromDict), BOOL, cnv->h) ;

  if (0)
    {
      BOOL b = TRUE, bb[5] ;
      for (k = 0 ; k < 6 ; k++)
	 channelPut (cnv->done, &b, BOOL) ;

      k = channelMultiGet (cnv->done, bb, 5, BOOL) ;
      k += channelMultiGet (cnv->done, bb, 5, BOOL) ;
      fprintf (stderr, "--- k = %d\n", k) ;
    }      
  /* allocate all the memory in advance */
  for (run = 1 ; run <= dictMax (cnv->runDict) ; run++)
    {
      RC *rc = arrayp (cnv->runs, run, RC) ;
      rc->cover = arrayHandleCreate (cMax, float, cnv->h) ;
      rc->geneCover = arrayHandleCreate (gMax, float, cnv->h) ;
    }
  fprintf (stderr, "--- created rc->cover with cMax = %d\n", cMax) ;
  /* fill the todo channel, it is large enough so that it will not block */
  for (run = 1 ; run <= dictMax (cnv->runDict) ; run++)
    if (arr (cnv->runs, run, RC).type)
      for (chrom = 1 ; chrom <= dictMax (cnv->chromDict) ; chrom++)
	{
	  gcv.cnv = cnv ; gcv.run = run ; gcv.chrom = chrom ; 
	  channelPut (cnv->todo, &gcv, GCV) ;
	  nToDo++ ;
	}
  channelClose (cnv->todo) ;

  fprintf (stderr, "// cnvParseOneCoverage  prepared %d runs * %d chroms = %d requests using %d threads \n"
	   , dictMax (cnv->runDict), dictMax (cnv->chromDict), nToDo,  (cnv->max_threads > 1 ? cnv->max_threads : 1)
	   ) ;
  
  /* launch k concurent threads */ 
  fprintf (stderr, "---- %s : cnvParseCoverage  starts %d threads\n"
	   , timeShowNow()
	   , (cnv->max_threads > 3 ? cnv->max_threads - 2 : 1)
	  ) ;
  for (k = 0 ; k <  (cnv->max_threads > 3 ? cnv->max_threads - 2 : 1)  ; k++)
    wego_go (cnvParseOneCoverage, cnv, CNV) ;

  /* wait for the results */
  {
    BOOL bb[1001] ;
    while (nDone < nToDo)
      {
	int i, j ;
	j = channelMultiGet (cnv->done, bb, 1000, BOOL) ;
	nDone += j ;
	for (i = 0 ; i < j ; i++)
	  if (bb[i]) 	
	    nDoneOk++ ;
	if (0)
	  if (j > 1 || nDone > 12700)
	    fprintf (stderr, ".... j=%d todo=%d done=%d\n", j, nToDo, nDone) ;
      }
  }
  fprintf (stderr, "---- %s : cnvParseCoverage  received %d results, %d success\n"
	   , timeShowNow ()
	   , nDone, nDoneOk
	   ) ;
  return nDoneOk ;
}  /* cnvParseCoverage */

/*************************************************************************************/
/* cumulate the covers per run and per coveron, distinguishing by run type
 * compute the median values 
 */
static void cnvCumulCoverage (CNV *cnv, BOOL byCoveron)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL byGene = byCoveron ? FALSE : TRUE ;
  int run, coveron, type, kCov ;
  Array medians = arrayHandleCreate (dictMax (cnv->runDict) + 1, float, h) ;
  Array kk = arrayHandleCreate (dictMax (cnv->runDict) + 1, float, h) ;
  GC *gc ; 
  RC *rc ;
  int chromX = 0, chromY = 0 ;
  Array myCoverons = byCoveron ? cnv->coverons : cnv->genes ;
  DICT *myDict = byCoveron ? cnv->coveronDict : cnv->geneDict ;
  int runMax =  dictMax (cnv->runDict) + 1 ;
  int *nRuns = cnv->nRuns ;

  dictFind (cnv->chromDict, "X", &chromX) ;
  dictFind (cnv->chromDict, "Y", &chromY) ;

  for (type = 0 ; type < 2 ; type++)
    nRuns[type] = 0 ;
  for (run = 1, rc = arrayp (cnv->runs, run, RC)  ; run < arrayMax (cnv->runs) ; rc++, run++)
    {
      rc->cumul = 0 ;
      if (rc->type > 0 && rc->type < 3) nRuns[rc->type - 1]++ ;
    }
  fprintf (stderr, "// There are %d control runs plus %d spectator runs\n", nRuns[0], nRuns[1]) ;
  if (! byCoveron)
    {
      int gene, geneMax = arrayMax (cnv->genes) ;

      for (gene = 1, gc = arrp (cnv->genes, gene, GC) ; gene < geneMax ; gene++, gc++)
	if (gc->gene && gc->gene < arrayMax (cnv->gene2coverons))
	  {
	    KEYSET g2c = array (cnv->gene2coverons, gc->gene, KEYSET) ;
	    int i, type, iMax = g2c ? keySetMax (g2c) : 0 ;
	    int run, runMax = arrayMax(cnv->runs) ;
	    RC *rc2 ;

	    gc->ln = 0 ; 
	    for (i = 0 ; i < iMax ; i++)
	      {
		int coveron2 = keySet (g2c, i) ;
		GC *gc2 = coveron2 < arrayMax (cnv->coverons) ?  arrp (cnv->coverons, coveron2, GC) : 0 ;
		int ln = gc2 ? (gc2->a1 < gc2->a2 ?  gc2->a2 - gc2->a1 + 1 : gc2->a1 - gc2->a2 + 1) : 0 ;
		gc->ln += ln ;
		gc->u2u += gc2->u2u ;
		gc->u2nu += gc2->u2nu ;
	      }
	    for (run = 1, rc2 = arrayp (cnv->runs, run, RC)  ; run < runMax ; rc2++, run++)
	      for (type = 1 ; type <= 2 ; type++)
		if (rc2->type == type)
		  {
		    for (i = 0 ; i < iMax ; i++)
		      {
			int coveron2 = keySet (g2c, i) ;
			
			if (rc2->cover && coveron2 < arrayMax (rc2->cover))
			  {
			    float z = arr (rc2->cover, coveron2, float) ;
			    float *zp = arrayp (rc2->geneCover, gc->gene, float) ;
			    *zp += z ;
			    gc->cumul[type-1] += z ;
			    if (0 && gene == 1000)
			      fprintf (stderr, "* %s %s %s %f %f\n"
				       , dictName(cnv->geneDict, gene)
				       , dictName(cnv->runDict, run)
				       , dictName(cnv->coveronDict, coveron2)
				       , z
				       ,  array (rc2->geneCover, gc->gene, float)
				       ) ;

			  }
		      }
		    if (0 && gene == 1000)
		      fprintf (stderr, "******** %s %s %f\n"
			       , dictName(cnv->geneDict, gene)
			       , dictName(cnv->runDict, run)
			       ,  array (rc2->geneCover, gc->gene, float)
			       ) ;
		  }
	  } 
    }
  
  /* we must loop by coveron since we wish to compute the median coverage of each coveron */
  /* we loop twice to establish the sex */
  
  for (kCov = 0, coveron = 1, gc = arrp (myCoverons, coveron, GC) ; coveron <= dictMax (myDict); coveron++, gc++)
    if (
	((gc->coveron && byCoveron) || (gc->gene && byGene)) &&
	(gc->chrom == chromX || gc->chrom == chromY)
	)
      for (type = 1 ; type <= 2 ; type++)
	{
	  RC *rc ;
	  Array myCover ;

	  arrayMax (kk) = 0 ;
	  
	  for (run = 1, rc = arrayp (cnv->runs, run, RC)  ; run < runMax ; rc++, run++)
	    { 
	      myCover =  byCoveron ? rc->cover : rc->geneCover ;
	      if (rc->type == type && myCover && coveron < arrayMax (myCover))
		{
		  float z = arr (myCover, coveron, float) ;
		  if (gc->chrom == chromX)
		    rc->cumulX += z ;
		  if (gc->chrom == chromY)
		    rc->cumulY += z ;
		}
	    }
	}
   
  for (run = 1, rc = arrayp (cnv->runs, run, RC)  ; run < runMax ; rc++, run++)
    { 
      if (0 && rc->cumulX > 0)  
	fprintf (stdout, "%s\t%f\t%f\t%f\n"
		 , dictName(cnv->runDict, rc->run)
		 , rc->cumulX, rc->cumulY
		 , 100.0 * rc->cumulY/(1 + rc->cumulX + rc->cumulY)
		 ) ;
      if (rc->cumulX > 0 && rc->cumulX > 3.0 * (rc->cumulX + rc->cumulY))
	rc->sex = 2 ; /* fille */
      else
	{
	  rc->sex = 1 ; /* garcon : on double les X */
	  if (rc->type && rc->type < 3)
	    cnv->nRunsY[rc->type - 1]++ ;
	}
    }
  if (cnv->nRunsY[0] < 1) cnv->nRunsY[0] = 1 ;
  if (cnv->nRunsY[1] < 1) cnv->nRunsY[1] = 1 ;
  
  for (kCov = 0, coveron = 1, gc = arrp (myCoverons, coveron, GC) ; coveron <= dictMax (myDict); coveron++, gc++)
    if ((gc->coveron && byCoveron) || (gc->gene && byGene))
      for (type = 1 ; type <= 2 ; type++)
	{
	  int k = 0, k0 = 0, k1 = 0 ;
	  RC *rc ;
	  Array myCover ;

	  arrayMax (kk) = 0 ;
	  
	  for (run = 1, rc = arrayp (cnv->runs, run, RC)  ; run < runMax ; rc++, run++)
	    { 
	      /* if gc->chrom == Y : only compute the median among boys */
	      if (chromY && gc->chrom == chromY && rc->sex != 1)
		continue ;
	      myCover =  byCoveron ? rc->cover : rc->geneCover ;
	      if (rc->type == type && myCover && coveron < arrayMax (myCover))
		{
		  float z = arr (myCover, coveron, float) ;
		  if ( 1 &&
		       ((rc->sex == 1 && chromX && gc->chrom == chromX) ||  (chromY && gc->chrom == chromY))
		       )
		    { z *= 2 ; arr (myCover, coveron, float) = z ; } /* double while computing norms, half when exporting */
		  array (kk, k++, float) = z ;
		  rc->cumul += z ; k0++ ; if (z > 0) k1++ ;
		  gc->cumul[type-1] += z ; 

		}
	    }
	  if (! keySetMax(kk))
	    continue ;
	  arraySort (kk, floatOrder)  ;
	  k = arrayMax (kk)/2 ;
	  if (gc->cumul[0] > 800 &&
	      (type > 1 || (0 * 3 * k0 < 4 * k1)) /* a coveron should be detected in most wild type runs */
	      )
	    gc->median[type-1] = arr (kk, k, float) ;
	  else
	    gc->median[type-1] = 0 ;
	  if (type == 1 && gc->median[type-1] > 0)
	    array (medians, kCov++, float) = gc->median[type-1] ; 
	}
 
  /* then we compute the median M of the median coverage of each coveron
   * then sigma : the median of the absolute value of the deviations up and down
   * the we remove points at zero or futher than 3 sigma
   */
  if (1)
    {
      float M = 0, M1, M2, sigmam, sigmap ;
      int k, kMax = arrayMax (medians), kp, km, nDiscardedP = 0, nDiscardedM = 0 ;
      Array sortedMedians = arrayHandleCopy (medians, h) ; 
      Array dp = arrayHandleCreate (kMax, float, h) ;
      Array dm = arrayHandleCreate (kMax, float, h) ;
      
      arraySort (sortedMedians, floatOrder)  ;  
      M = array (sortedMedians, kMax/2, float) ;
      
      /* construct the table of deviations around the median */
      for (k = kp = km = 0 ; k < kMax ; k++)
	{
	  float z = array (medians, k, float) - M ;
	  if (z > 0)
	    array (dp, kp++, float) = z ;
	  else
	    array (dm, km++, float) = -z ;
	}
      
      arraySort (dp, floatOrder)  ;  k = arrayMax (dp)/2 ; sigmap = array (dp, k, float) ;
      arraySort (dm, floatOrder)  ;  k = arrayMax (dm)/2 ; sigmam = array (dm, k, float) ;
      
    /* flag all coverons whese median breaks the 3 MAD rule : m outside [M - 3 sm, M + 3 sp] */
      M1 = M - 3 * sigmam ; if (M1 < .1) M1 = .1 ; if (M1 < M/10) M1 = M/10 ;
      M2 = M + 3 * sigmap ; 
      if (byCoveron)  /*  discard 3 MADs */
	for (coveron = 1, gc = arrp (myCoverons, coveron, GC) ; coveron <= dictMax (myDict); coveron++, gc++)
	  {
	    float z =  gc->median[0] ;
	    if (0 && z < M1)
	      nDiscardedM++ ; 
	    if (z > M2)
	      nDiscardedP++ ;
	  }
      cnv->M = M ; cnv->M1 = M1 ; cnv->M2 = M2 ; 
      fprintf (stderr, " Median of Median coverage %.1f, sigmam = %.1f,  sigmap = %.1f, M - 3 sigmaM = %.1f, M + 3 sigmaP = %.1f, discared operons %d + %d, kept %d\n"
	       , M, sigmam, sigmap, M1, M2, nDiscardedM, nDiscardedP, kMax - nDiscardedM - nDiscardedP 
	       ) ;
    }
  
  ac_free (h) ;
  return ;
} /* cnvCumulCoverage  */

/*************************************************************************************/

static void cnvReportCoverage (CNV *cnv, BOOL byCoveron)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL byGene = byCoveron ? FALSE : TRUE ;
  ACEOUT ao = aceOutCreate (cnv->outFileName, byCoveron ? ".cnv_compare_coverons.OLD.txt" : ".cnv_compare_genes.OLD.txt", cnv->gzo, h) ;
  int run, coveron, nn1 = 0, nn2 = 0 ;
  int runMax = dictMax (cnv->runDict) ;
  RC *rc ;
  GC *gc ;
  float  globalCumul1 = 0, globalCumul2 = 0 ;
  int cMax = arrayMax (cnv->coverons) ;


  aceOutDate (ao, "##", cnv->title) ;
  
  
  aceOutf (ao, "# Coveron\tChromosome\ta1\ta2\tCoveron Length (bp)\tGene\tAbsolute average coverage in normal tissues\tAbsolute average coverage in tumor tissues\trelative average coverage : tumor/normal\t\tRun") ;
  for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
    if (rc->type) 
      aceOutf (ao, "\t%s", dictName (cnv->runDict, run)) ;
  aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t\tTitle") ;
  for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
    if (rc->type)
      aceOutf (ao, "\t%s", stackText (cnv->info, rc->title)) ;
  aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t\tSorting Title") ;
  for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
    if (rc->type)
      aceOutf (ao, "\t%s",  stackText (cnv->info, rc->sortingTitle)) ;
  aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t\tSorting Title 2") ;
  for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
    if (rc->type) 
      aceOutf (ao, "\t%s",  stackText (cnv->info, rc->sortingTitle2)) ;
  aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t\tOther Title") ;
  for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
    if (rc->type) 
      aceOutf (ao, "\t%s",  stackText (cnv->info, rc->otherTitle)) ;
  aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t\tRunID") ;
  for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
    if (rc->type)
      aceOutf (ao, "\t%s",  stackText (cnv->info, rc->runId)) ;
  aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t\tSample") ;
  for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
    if (rc->type)
      aceOutf (ao, "\t%s",  stackText (cnv->info, rc->sample)) ;
  
  if (1) /* normalize relative to expected value = sum this row * sum of this column / grand total */
    {
      for (coveron = 1, gc = arrp (cnv->coverons, coveron, GC) ; coveron < cMax ; coveron++, gc++)
	{
	  if ((byCoveron && ! gc->coveron)  || (byGene && gc->coveron)) continue ;
	  if (gc->median[0] <= 0) continue ;
	  globalCumul1 += gc->median[0] ;
	  globalCumul2 += gc->median[1] ;
	}
      if (globalCumul1 == 0) globalCumul1 = 1 ;
      if (globalCumul2 == 0) globalCumul2 = 1 ;

      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	   {
	     if (rc->type >= 1) 
	       for (coveron = 1, gc = arrp (cnv->coverons, coveron, GC) ; coveron < cMax ; coveron++, gc++)
		 if ((byCoveron && gc->coveron)  || (byGene && ! gc->coveron)) 
		   {
		     float z = rc->cumul * gc->median[0] ;  
		     if (gc->median[0] <= 0) continue ;
		     array (rc->cover, coveron, float) = 0 ; /* make room */
		     if (z == 0)  
		       arr (rc->cover, coveron, float) = -1 ;
		     else
		       arr (rc->cover, coveron, float) *= globalCumul1 / z  ;
		   }
	   }
    }
  
  
  for (coveron = 1, gc = arrp (cnv->coverons, coveron, GC) ; coveron < cMax ; coveron++, gc++)
    {
      float z ;
      int len = gc->a2 - gc->a1 + 1 ;
      if (len < 1) len = 1 ;
      
      if ((byCoveron && ! gc->coveron)  || (byGene &&  gc->coveron)) 
	continue ; 
      if (1) /* discard 3 MAD */ 
	{
	  float z =  gc->median[0] ;
	  if (z < cnv->M1)
	    continue ;
	  if (z > cnv->M2)
	    continue ;
	}
      nn1++ ;
      if (byCoveron)
	{
	  aceOutf (ao, "\n%s\t%s\t%d\t%d\t%d", dictName (cnv->coveronDict, coveron), dictName (cnv->chromDict, gc->chrom), gc->a1, gc->a2, len) ;
	  if (gc->coveron2genes && keySetMax (gc->coveron2genes)) 
	    {
	      int i ;
	      nn2++ ;
	      for (i = 0 ; i < keySetMax (gc->coveron2genes) ; i++)
		aceOutf (ao, "%s%s", i ? ", " : "\t",  dictName (cnv->geneDict, keySet (gc->coveron2genes, i))) ;
	    }
	  else if (gc->gene) 
	    {
	      nn2++ ;
	      aceOutf (ao, "\t%s", dictName (cnv->geneDict, gc->gene)) ;
	    }
	  else 
	    aceOutf (ao, "\t") ;
	}
      else
	aceOutf (ao, "\n%s\t%s\t%d\t%d\t%d\t", dictName (cnv->geneDict, gc->gene), dictName (cnv->chromDict, gc->chrom), gc->a1, gc->a2, len) ;

      aceOutf (ao, "\t%.0f", gc->median[0]/len) ;
      aceOutf (ao, "\t%.0f", gc->median[1]/len) ;
      
      z =  gc->median[0] + gc->median[1] ; if (z == 0) z = 1 ;
      aceOutf (ao, "\t%.2f", (gc->median[1]/z) * (globalCumul1 + globalCumul2)/globalCumul2) ;
      
      aceOutf (ao, "\t\t") ;
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	{
	  float z =  array (rc->cover, coveron, float) ;
	  if (rc->type)
	    {
	      if (z >= 0) 
		aceOutf (ao, "\t%g", arr (rc->cover, coveron, float)) ;
	      else
		aceOutf (ao, "\t") ;
	    }	    
	}
    }
  aceOutf (ao, "\n") ;
  fprintf (stderr, "# Exported %d coveron coverages, %d matching a gene, in file %s\n"
	   , nn1, nn2
	   , aceOutFileName (ao)
	   ) ;
  
  for (coveron = 1, gc = arrp (cnv->coverons, coveron, GC) ; coveron < cMax ; coveron++, gc++)
    ac_free (gc->coveron2genes) ;
  ac_free (h) ;

} /* cnvReportCoverage */

/*************************************************************************************/

static void cnvExportCoverage (CNV *cnv, BOOL byCoveron, BOOL normalize)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL byGene = byCoveron ? FALSE : TRUE ;
  AC_TABLE groups = 0, groups2runs = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  ACEOUT ao ;
  ACEOUT aoC = aceOutCreate (cnv->outFileName, byCoveron ? 
			     (normalize ? ".cnv_compare_coverons.scaled.txt" : ".cnv_compare_coverons.raw.txt") : 
			     (normalize ? ".cnv_compare_genes.scaled.txt" : ".cnv_compare_genes.raw.txt") 
, cnv->gzo, h) ;
  ACEOUT aoL = aceOutCreate (cnv->outFileName, byCoveron ? 
			     (normalize ? ".cnv_compare_coverons.lost.scaled.txt" : ".cnv_compare_coverons.lost.raw.txt") : 
			     (normalize ? ".cnv_compare_genes.lost.scaled.txt" : ".cnv_compare_genes.lost.raw.txt")
, cnv->gzo, h) ;
  ACEOUT aoH = aceOutCreate (cnv->outFileName, byCoveron ? 
			     (normalize ? ".cnv_compare_coverons_histo.scaled.txt" : ".cnv_compare_coverons_histo.raw.txt") :
			     (normalize ? ".cnv_compare_genes_histo.scaled.txt" : ".cnv_compare_genes_histo.raw.txt")
			     , cnv->gzo, h) ;
  int iao, i, run, coveron, nn1 = 0, nn2 = 0 ;
  int runMax = dictMax (cnv->runDict) ;
  BOOL isLost ;
  RC *rc ;
  GC *gc ;
  float  globalCumul1 = 0, globalCumul2 = 0 ;
  Array myCoverons = byCoveron ? cnv->coverons : cnv->genes ;
  int cMax = arrayMax (myCoverons) ;
  int chromX = 0, chromY = 0 ;
  KEYSET ksG2R[100] ;
  int iGroup, iGroupMax = 1 ;
  int *nRuns = cnv->nRuns ;
  int *nRunsY = cnv->nRunsY ;
  const char *errors = 0 ;
  Array aGroup = arrayHandleCreate (128, float, h) ;
  Array aY = arrayHandleCreate (1 + arrayMax (cnv->runs), float, h) ;

  memset (ksG2R, 0, sizeof (ksG2R)) ;

  dictFind (cnv->chromDict, "X", &chromX) ;
  dictFind (cnv->chromDict, "Y", &chromY) ;

  if (cnv->db && cnv->project)
    {
      vtxtClear (txt) ;
      vtxtPrintf (txt, "colonne 1 \n class Run \n Condition Runs && CNV && ! Add_counts && project == \"%s\"\n", cnv->project) ;
      vtxtPrintf (txt, "colonne 2 \nFrom 1\nClass Text\nTag Sorting_title\n\n") ;
      groups = ac_tablemaker_table (cnv->db, vtxtPtr (txt), 0, ac_tablemaker_text , 0 , "2+1", &errors, h) ;
 
      vtxtClear (txt) ;
      vtxtPrintf (txt, "colonne 1 \n class Run \n Condition Runs && ! Add_counts && project == \"%s\"\n", cnv->project) ;
      vtxtPrintf (txt, "colonne 2 \nFrom 1\nClass Run\nTag Runs\n\n") ;
      vtxtPrintf (txt, "colonne 3 \nFrom 1\nClass Run\nTag Subgroup\n\n") ;
      groups2runs = ac_tablemaker_table (cnv->db, vtxtPtr (txt), 0, ac_tablemaker_text , 0 , 0, &errors, h) ;

      iGroupMax = 1 + (groups ? groups->rows : 0) ; 
      if (iGroupMax >= 100) iGroupMax = 99 ;

      for (iGroup = 1 ; iGroup < iGroupMax ; iGroup++)
	{
	  int jr, kMax = 0 ;
	  KEY key =  ac_table_key (groups, iGroup - 1, 0, 0) ;
	  KEYSET ks = ksG2R[iGroup] = keySetHandleCreate (h) ;

	  for (jr = 0 ; key && groups2runs && jr < groups2runs->rows ; jr++)
	    {
	      int run1 = 0 ;
	      if (ac_table_key (groups2runs, jr, 0, 0) == key && 
		  dictFind (cnv->runDict, ac_table_printable (groups2runs, jr, 1, "-"), &run1)
		  )
		keySet (ks, kMax++) = run1 ;
	      if (ac_table_key (groups2runs, jr, 0, 0) == key && 
		  dictFind (cnv->runDict, ac_table_printable (groups2runs, jr, 2, "-"), &run1) 
		  )
		keySet (ks, kMax++) = run1 ;
	    }
	  keySetSort (ks) ;
	  keySetCompress (ks) ;
          fprintf (stderr, "Group %d has %d runs\n", iGroup, keySetMax (ks)) ;
	}
    }
  
  
  if (1) /* normalize relative to expected value = sum this row * sum of this column / grand total of control samples */
    {
      for (coveron = 1, gc = arrp (myCoverons, coveron, GC) ; coveron < cMax ; coveron++, gc++)
	{
	  if ((byCoveron && ! gc->coveron)  || (byGene && gc->coveron)) continue ;
	  if (gc->median[0] <= 0) continue ;
	  globalCumul1 += gc->median[0] ;
	  globalCumul2 += gc->median[1] ;
	}
      if (globalCumul1 == 0) globalCumul1 = 1 ;
      if (globalCumul2 == 0) globalCumul2 = 1 ;
      
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	memset (rc->histo, 0, sizeof (rc->histo)) ;
      
      /* renormalize chromY (after we computed the global cumul) */
      if (chromY && gc->chrom == chromY)
	{
	  int t ;
	  Array myCover = byCoveron ? rc->cover : rc->geneCover ;
	  for (t = 1; t < 3 ; t++)  /* 1:control, 2:any */
	    {
	      int n = 0 ;
	      aY = arrayReCreate (aY, 1 + arrayMax (cnv->runs), float) ;
	      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
		if (rc->type == t && rc->sex == 1) 
		  array (aY, n++, float) = arr (myCover, coveron, float) ;
	      arraySort (aY, floatOrder) ;
	      gc->median[t-1] = array (aY, n/2, float) ;
	    }
	}
				     
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	{
	  if (rc->type >= 1) 
	    for (coveron = 1, gc = arrp (myCoverons, coveron, GC) ; coveron < cMax ; coveron++, gc++)
	      {
		Array myCover = byCoveron ? rc->cover : rc->geneCover ;
		float z = rc->cumul * gc->median[0] ;  
		
		if (coveron >= arrayMax (myCover)) break ;
		/* Limit to autosomes (exclude X, Y which are artificially haploid in boys  */
		if (0 && (gc->chrom == chromX || gc->chrom == chromY))
		  continue ;
		
		if ((byCoveron && ! gc->coveron)  || (byGene &&  gc->coveron)) 
		  continue ; 
		isLost = FALSE ;
		if (1) /* discard 3 MAD */ 
		  {
		    float z1 =  gc->median[0] ;
		    if (0 && z1 == 0)
		      continue ;
		    if (0 && z1 < cnv->M1)
		      isLost = TRUE ;
		    if (0 && byCoveron && z1 > cnv->M2)
		      continue ;
		    if ( z1 < cnv->M1)
		      z = rc->cumul * cnv->M1 ;  
		  }
		
		if (isLost)
		  { if (normalize) if (rc->cumul) arr (myCover, coveron, float) /=  rc->cumul ; }
		if (z == 0)  
		  { if (normalize)arr (myCover, coveron, float) = -1 ; }
		else
		  { if (normalize) arr (myCover, coveron, float) *= 2 * globalCumul1 / z  ;} /* 2 because the normal is to be diploid */
		
		z = arr (myCover, coveron, float) ;
		if (0 &&
		    ((rc->sex == 1 && chromX && gc->chrom == chromX) || ( chromY && gc->chrom == chromY))
		    )
		  arr (myCover, coveron, float) = z = z/2 ;
		if (0 && ! isLost)
		  {
		    int zk = 20 * z + .4999 ; 
		    if (zk > 100) zk = 100 ; 
		    if (zk < 0) zk = 0 ; 
		    if (zk >= 0)
		      rc->histo[zk] ++ ;
		}		     
	      }
	}
    }
 
  for (iao = 0 ; iao < 3 ; iao++)
    {
      char groupBuf[9 * iGroupMax + 16] ;

      memset (groupBuf, 0, sizeof(groupBuf)) ;
      memset (groupBuf, '\t', 6 + 8 * iGroupMax) ;

      switch (iao)
	{
	case 0 : ao = aoC ; break ;
	case 1 : ao = aoL ; break ;
	case 2 : ao = aoH ; break ;
	}
      aceOutDate (ao, "##", cnv->title) ;
      aceOutf (ao, "## Genome wide %s gene base normalized relative to the median of the normal samples, excluding the 3 MAD.\n"
	       , byCoveron ? "coveron" : "gene"
	       ) ;
      aceOutf (ao, "## Only paired end aligned over 90 %% of their length with less than 3 missmatches are considered.\n") ;
      if (ao == aoH)
	aceOutf (ao, "## The histograms are limited to the autosomes, i.e. excluding the X and Y chromosomes.\n") ;
      aceOutf (ao, "## Unicity is computed as the ratio 100*u/u+nu where u is the number of bases aligned uniquely and nu non-uniquelly inside the coveron, summed over all runs\n") ;
      
      aceOutf (ao, "# %s\tChromosome\ta1\ta2\t%s Length (bp)\tUnicity\tGene\tCumul in normal tissues\tAverage coverage in normal tissues\tAverage coverage in tumor tissues\tMedian coverage in normal tissues\tMedian coverage in tumor tissues\trelative median coverage : tumor/normal\t"
	       , byCoveron ? "Coveron" : "Gene" , byCoveron ? "Coveron" : "Gene"
	       ) ;
      aceOutf (ao, "\tAverage\tMedian\tSamples\tLost (< 0.1 copy)\tLow (0.1 to 1.5 copy)\tMid (1.5 to 2.5 copies)\tHigh (2.5 to 3.5 copies)\tVery_high (over 3.5 copies)\t\t") ;
      for (iGroup = 1 ; iGroup < iGroupMax ; iGroup++)
	{
	  const char *ccp = ac_table_printable (groups, iGroup - 1, 0, " ")  ;
	  aceOutf (ao, "\tAverage in %s\tMedian in %s\tSamples\tLost in %s\tLow in %s\tMid in %s\tHigh in %s\tVery_high in %s\t"
		   , ccp
		   , ccp
		   , ccp
		   , ccp
		   , ccp
		   , ccp
		   , ccp
		   ) ;
	}
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	if (rc->type) 
	  aceOutf (ao, "\t%s", dictName (cnv->runDict, run)) ;
      aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t%sTitle", groupBuf) ;

      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	if (rc->type)
	  aceOutf (ao, "\t%s", stackText (cnv->info, rc->title)) ;
      aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t%sSorting Title", groupBuf) ;
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	if (rc->type)
	  aceOutf (ao, "\t%s",  stackText (cnv->info, rc->sortingTitle)) ;
      aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t%sSorting Title 2", groupBuf) ;
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	if (rc->type) 
	  aceOutf (ao, "\t%s",  stackText (cnv->info, rc->sortingTitle2)) ;
      aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t%sOther Title", groupBuf) ;
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	if (rc->type) 
	  aceOutf (ao, "\t%s",  stackText (cnv->info, rc->otherTitle)) ;
      aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t%sRunID", groupBuf) ;
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	if (rc->type)
	  aceOutf (ao, "\t%s",  stackText (cnv->info, rc->runId)) ;
      aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t%sSample", groupBuf) ;
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	if (rc->type)
	  aceOutf (ao, "\t%s",  stackText (cnv->info, rc->sample)) ;
      aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t%sCumulated aligned base pairs", groupBuf) ;
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	if (rc->type)
	  aceOutf (ao, "\t%f", rc->cumul) ;
      aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t%sSorting Title", groupBuf) ;
      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
	if (rc->type)
	  aceOutf (ao, "\t%s",  stackText (cnv->info, rc->sortingTitle)) ;
      
  
      if (iao == 2)
	{
	  for (i = 0 ; i <= 100 ; i++)
	    {
	      aceOutf (ao, "\n#\t\t\t\t\t\t\t\t\t\t\t%s%f", groupBuf, 5.0 * i / 100.0) ;
	      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
		if (rc->type)
		  aceOutf (ao, "\t%d", rc->histo[i]) ;
	    }
	}
      else
	{
	  for (coveron = 1, gc = arrp (myCoverons, coveron, GC) ; coveron < cMax ; coveron++, gc++)
	    {
	      float z ;
	      int len = gc->a2 - gc->a1 + 1 ;
	      if (len < 1) len = 1 ;
	      if (! byCoveron && gc->ln) len = gc->ln ;
	      if ((byCoveron && ! gc->coveron)  || (byGene &&  gc->coveron)) 
		continue ; 
	      isLost = FALSE ;
	      if (1) /* discard 3 MAD */ 
		{
		  float z =  gc->median[0] ;
		  if (gc->median[0] + gc->median[1] == 0)
		    continue ;
		  if (1 && z < cnv->M1)
		    isLost = TRUE ;
		  if (1 && byCoveron && z > cnv->M2)
		    continue ;
		}
	      if ((isLost && iao != 1) || (!isLost && iao == 1))
		continue ;
	      nn1++ ;
	      if (byCoveron)
		{
		  aceOutf (ao, "\n%s\t%s\t%d\t%d\t%d\t"
			   , dictName (cnv->coveronDict, coveron)
			   , dictName (cnv->chromDict, gc->chrom)
			   , gc->a1, gc->a2, len
			   ) ;
		  if (gc->u2u + gc->u2nu > 100)
		    aceOutf (ao, "%.1f",  100.0*(gc->u2u)/(gc->u2u+gc->u2nu +.0001)) ;

		  if (gc->coveron2genes && keySetMax (gc->coveron2genes)) 
		    {
		      int i ;
		      nn2++ ;
		      for (i = 0 ; i < keySetMax (gc->coveron2genes) ; i++)
			aceOutf (ao, "%s%s", i ? ", " : "\t",  dictName (cnv->geneDict, keySet (gc->coveron2genes, i))) ;
		    }
		  else if (gc->gene) 
		    {
		      nn2++ ;
		      aceOutf (ao, "\t%s", dictName (cnv->geneDict, gc->gene)) ;
		    }
		  else 
		    aceOutf (ao, "\t") ;
		}
	      else
		{
		  int a1 = gc->isDown ? gc->a1 : gc->a2 ;
		  int a2 = gc->isDown ? gc->a2 : gc->a1 ;
		  aceOutf (ao, "\n%s\t%s\t%d\t%d\t%d\t"
			   , dictName (cnv->geneDict, gc->gene)
			   , dictName (cnv->chromDict, gc->chrom)
			   , a1, a2, len
			   ) ;
		  if (gc->u2u + gc->u2nu > 100)
		    aceOutf (ao, "%.1f",  100.0*(gc->u2u)/(gc->u2u+gc->u2nu+.0001)) ;
		  aceOutf (ao, "\t%s", dictName (cnv->geneDict, gc->gene)) ;
		}
	      
	      aceOutf (ao, "\t%.3f", gc->cumul[0]) ;
	      if (chromY && gc->chrom == chromY)
		{
		  aceOutf (ao, "\t%.3f", gc->cumul[0]/(nRunsY[0]*len)) ;
		  aceOutf (ao, "\t%.3f", gc->cumul[1]/(nRunsY[1]*len)) ;
		}
	      else
		{
		  aceOutf (ao, "\t%.3f", gc->cumul[0]/(nRuns[0]*len)) ;
		  aceOutf (ao, "\t%.3f", gc->cumul[1]/(nRuns[1]*len)) ;
		}
	      aceOutf (ao, "\t%.3f", gc->median[0]/len) ;
	      aceOutf (ao, "\t%.3f", gc->median[1]/len) ;
	      
	      z =  gc->median[0] + gc->median[1] ; if (z == 0) z = 1 ;
	      aceOutf (ao, "\t%.2f", (gc->median[1]/z) * (globalCumul1 + globalCumul2)/globalCumul2) ;
	      
	      /* compute the number of occurence in the 5 categories */
	      for (iGroup = 0 ; iGroup < iGroupMax ; iGroup++)
		{
		  int jj, N = 0 ;
		  float X = 0, limit[] = { 0, .1, 1.5, 2.5, 3.5, 99999} ;

		  aGroup = arrayReCreate (aGroup, 128, float) ;
		  aceOutf (ao, "\t\t") ;
		  for (jj = 0 ; jj <= 5 ; jj++)
		    {
		      int n = 0 ;
		      BOOL ok = FALSE ;
		      float lim1 = jj > 0 ? limit[jj-1] : -99999 ;
		      float lim2 = jj > 0 ? limit[jj] : 999999 ;
		      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
			if (rc->type >= 1)
			  {
			    switch (iGroup)
			      {
			      case 0: 
				ok = TRUE ;  /* no condition */
				break ;
			      default:
				ok = keySetFind (ksG2R[iGroup], run, 0) ;
			      }
			    if (ok)
			      {  
				Array myCover =  byCoveron ? rc->cover : rc->geneCover ;
				float z = coveron < arrayMax (myCover) ?  arr (myCover, coveron, float) : 0 ;
				/* was done before 
				  if ((rc->sex == 1 && chromX && gc->chrom == chromX) || gc->chrom == chromY)
				  z /= 2 ;
				*/
				if (z >= lim1 && z < lim2) 
				  n++ ;
				if (jj == 0) 
				  {
				    X += z ; array (aGroup, N, float) = z ;  N++ ; 
				  }
			      }
			  }
		      if (jj == 0)
			{
			  X = N > 0 ? X/N : 0 ;
			  aceOutf (ao, "\t%.2f", X) ;
			  arraySort (aGroup, floatOrder) ;
			  aceOutf (ao, "\t%.2f", array (aGroup, N/2, float)) ;
			}
		      
		      aceOutf (ao, "\t%d", n) ;
		    }
		}
	      aceOutf (ao, "\t") ;
	      for (run = 1, rc = arrp (cnv->runs, run, RC) ; run <= runMax ; run++, rc++)
		if (rc->type >= 1) 
		  {  
		    Array myCover =  byCoveron ? rc->cover : rc->geneCover ;
		    float z =   coveron < arrayMax (myCover) ? arr (myCover, coveron, float) : 0 ;
		    /* was done before 
		       if ((rc->sex == 1 && chromX && gc->chrom == chromX) || gc->chrom == chromY)
		      z /= 2 ;
		    */
		    
		    if (z >= 0) 
		      aceOutf (ao, "\t%g", z) ;
		    else
		      aceOutf (ao, "\t") ;
		  }	    
	    }
	}
      aceOutf (ao, "\n") ;
    }

  fprintf (stderr, "# Exported %d coveron coverages, %d matching a gene, in file %s\n"
	   , nn1, nn2
	   , aceOutFileName (ao)
	   ) ;
  
  ac_free (h) ;
  return ;
} /* cnvExportCoverage */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: cnv -compare -... \n"
	   "//      try: -h --help --hits_file_caption \n"
	   "// Example:  cnv \n"
	   
	   "// -o fileName : output file prefix, all exported files will start with that prefix\n"
  	   "// -gzo : the output file is gziped\n"
	   "//\n"
	   ) ;

  
  if (argc > 1)
    {
      fprintf (stderr,
	       "//\n// You said: %s\n", commandBuf) ;
      fprintf (stderr,
	       "// ########## ERROR: I do not understand the argument%s ", argc > 2 ? "s" : "") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (stderr, "%s ", argv[i]) ;
      fprintf (stderr, "\n") ;
    }
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  char *cp ;
  CNV cnv ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;
  const char *dbName = 0 ;

  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&cnv, 0, sizeof (CNV)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  cnv.h = h ;

  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, 1, argv) ;

  { 
    int ix ;
    for (ix = 0, cp = commandBuf ;  ix < argc && cp + strlen (argv[ix]) + 1 < commandBuf + 3900 ; cp += strlen (cp), ix++)
      sprintf(cp, "%s ", argv[ix]) ;
  }
  /* consume optional args */
  getCmdLineOption (&argc, argv, "-o", &(cnv.outFileName)) ;
  cnv.gzo = getCmdLineOption (&argc, argv, "--gzo", 0) ;
  getCmdLineOption (&argc, argv, "--project", &(cnv.project)) ;
  getCmdLineOption (&argc, argv, "--title", &(cnv.title)) ;
  getCmdLineOption (&argc, argv, "--chrom", &(cnv.chrom)) ;
  getCmdLineOption (&argc, argv, "--runs", &(cnv.runsSortedFileName)) ;
  getCmdLineOption (&argc, argv, "--runsAll", &(cnv.runsAllFileName)) ;
  getCmdLineOption (&argc, argv, "--runsControl", &(cnv.runsControlFileName)) ;
  getCmdLineOption (&argc, argv, "--gtitle", &(cnv.gtitleFileName)) ;
  getCmdLineOption (&argc, argv, "--geneSponge", &(cnv.geneSpongeFileName)) ;
  getCmdLineOption (&argc, argv, "--coveronSponge", &(cnv.coveronSpongeFileName)) ;
  getCmdLineOption (&argc, argv, "--coverage", &(cnv.coverageFileName)) ;
  cnv.unique = getCmdLineOption (&argc, argv, "--u", 0) ;
  cnv.non_unique = getCmdLineOption (&argc, argv, "--nu", 0) ;
  getCmdLineOption (&argc, argv, "--u2u", &(cnv.u2uFileName)) ;
  getCmdLineOption (&argc, argv, "--u2nu", &(cnv.u2nuFileName)) ;
  getCmdLineOption (&argc, argv, "--db", &dbName) ;

  if (dbName)
    {
      const char *errors = 0 ;
      
      cnv.db = ac_open_db (dbName, &errors);
      if (! cnv.db)
	messcrash ("Failed to open db %s, error %s", dbName, errors) ;
    }

  cnv.runDict = dictHandleCreate (1000, h) ;
  cnv.chromDict = dictHandleCreate (1000, h) ;
  cnv.geneDict = dictHandleCreate (1000, h) ;
  cnv.coveronDict = dictHandleCreate (1000, h) ;
  cnv.runs = arrayHandleCreate (256, RC, h) ;
  cnv.genes = arrayHandleCreate (100000, GC, h) ;
  cnv.coverons = arrayHandleCreate (100000, GC, h) ;
  cnv.ends = arrayHandleCreate (100000, EC, h) ;
  cnv.info = stackHandleCreate (10000, h) ;
  cnv.gene2coverons = arrayHandleCreate (10000, KEYSET, h) ;
  stackTextOnly (cnv.info) ;
  
  if (1)
    {
      cnv.max_threads = 4 ;
      arrayReport (-2) ;
      getCmdLineInt (&argc, argv, "--max_threads", &(cnv.max_threads)) ;
      
      if (cnv.max_threads < 2) 
	cnv.max_threads = 2 ;
      wego_max_threads (cnv.max_threads) ;
    }
  
  if (argc > 1)
    usage ("Unknown argument", argc, argv) ;

  fprintf (stderr, "// %s start : max_threads=%d\n", timeShowNow(),cnv.max_threads) ;
  cnvParseRunLists (&cnv, 0) ;
  cnvParseRunLists (&cnv, 2) ; /* all runs to be compared to the control (may optionally include control runs) */
  cnvParseRunLists (&cnv, 1) ; /* the control groups must be read after */
  
  if (cnv.gtitleFileName) 
    cnvParseTitles (&cnv) ;
 
  if (cnv.geneSpongeFileName) cnvParseSpongeFile (&cnv, 1) ;
  if (cnv.coveronSpongeFileName) cnvParseSpongeFile (&cnv, 2) ;
 
  if (arrayMax (cnv.ends))
    cnvCoverons2Genes (&cnv) ;
  
  if (cnv.u2uFileName && cnv.u2nuFileName)  cnvParseU2u (&cnv) ;

  if (cnv.coverageFileName  && arrayMax (cnv.coverons))
    {
      cnvParseCoverage (&cnv) ;
      cnvCumulCoverage (&cnv, FALSE) ;
      cnvExportCoverage (&cnv, FALSE, FALSE) ; /* by gene */
      cnvExportCoverage (&cnv, FALSE, TRUE) ; /* by gene */
      cnvReportCoverage (&cnv, FALSE) ; /* by gene */
      cnvCumulCoverage (&cnv, TRUE) ;
      cnvExportCoverage (&cnv, TRUE, FALSE) ; /* by coveron */
      cnvExportCoverage (&cnv, TRUE, TRUE) ; /* by coveron */
      cnvReportCoverage (&cnv, TRUE) ; /* by coveron */
    }
  if (cnv.db) 
    ac_db_close (cnv.db) ;
  ac_free (cnv.h) ;
  
  if (1)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// %s done: max memory %d Mb\n", timeShowNow(), mx) ;
     }
  
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
