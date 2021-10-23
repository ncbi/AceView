 /*  File: tricoteur.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2005
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

 * 2017_09_07
 * strategic idea

 * SNP and InDel DETECTION

 * Analyze the TCT alignments, the bar code is in column 20
 * Rationalize away the multiplets: same a1 coordinates and same barcode
 *    while creating the consensus of the match errors and overlaps
 * Count expression per gene and transcripts in all target classes
 *    add the hierarchic counts to detect the missing transcripts of the top hierarchy
 * Construct count and export all variations 
 * Locate all gene fusion mrna1, mrna2, coords, delta-x
 *   export the delta-x histogram
 * Construct the consesnus (bridge) for all SNPs and gene fusions
 * Rextend the partial reads and count the 'validations'

 * 1: Read the hits file, select a chromosomal zone of say 10 Mb
 *    Consider only the extremal coords of each read ali
      Construct the coverage wiggle


 * 2: Delineate the segments a--------X==X==X---------b
 *    such that X is a brk, x==X is shorter than DELTA a---X and X----b are longer
 *    and X is in a place covered >= 10 and #brk/#cover > 10%

 * 3: Read the fastc file of the reads touching one of the segments

 * 4: Loop on all segments, get the relevant ali, sort them
      Anchor at a, and extend base per base
        each time there is minor letter > 5%, biffurcate

 * 5: Merge the left and right consensus if both exist
      Export the consensus either as a del-ins or as an overhang

 *    DONE

 * SNP VALIDATION

 * export 31mers, count them in all the original fasta/fastq/fastc sequencing files

 * SNP REPORT

 * using the database, export nice reports

 */

#define VERSION "1.1"

#define MALLOC_CHECK   
#define ARRAY_CHECK   

#include "ac.h"
#include "channel.h"
#include "query.h"
#include "bql.h"

/* snp chaining distance */

typedef struct tsnpCallerTable {
  AC_HANDLE h ;
  const char *outFileName ;
  BOOL gzi, gzo ;
  ACEIN ai ; 
  ACEOUT ao ;

  const char *target_class ;
  const char *target ;
  int  t1, t2 ;

  BOOL intron_only ;
  const char *referenceGenome ;
  const char *inFileOfFileList ;
  const char *inFileList ;
  const char *project ;
  const char *force ;
  const char *zone ;
  const char *filter ;
  const char *remap2genome ;
  const char *remap2genes ;

 
  int minSnpCover, minSnpCount ;
  float minSnpFrequency ;
  int minIntron ;
  int maxIntron ;

  int max_threads ;
  int maxLanes ;

  int nAna ; /* default 4, number of analyzers launched in parallel wego */
  DICT *snpDict ;
  DICT *geneDict ;
  DICT *runDict ;
  DICT *targetDict ;
  DICT *chromDict ;
  DICT *target_classDict ;
  DICT *selectDict ;
  DICT *varTypeDict ;

  Array runs ;
  Array snps ;
  int snpDetected, snpExported ;
  KEYSET wdelta ;
  BOOL doubleReport ;

  BOOL mergeCounts ;
  BOOL dbReport ;
  BOOL dbTranslate ;
  BOOL dropMonomodal ;
  const char *dbName ;
  const char *wiggleDir ;
  AC_DB db ; 
  AC_TABLE runMetaData ;

  Array target2geneAtlas ;
  Array target2exonAtlas ;

  BOOL makeWords ;

  int pure, high, mid, low, ref, wellCovered ; /* prevalence */
  CHAN *getLaneInChan, *getLaneOutChan, *analyzeChan ;
  CHAN *doneChan ;
} TCT ;


typedef struct snpStruct {
  int snp ; /* self, needed after we sort the table */

  int type ; /* offset in typeDict <<= 3, then use or not split counts in bits 0x3 and slide  */
  int varType ; /* offset in varTypeDict */
  int target, a1, a2 ; /* target and position as read, can be edited by calling the database for the remapping */
  BOOL select ;

  int coverp, coverm, acoverp, acoverm ; /* coverage or mutant or wild type on plus/minus trand of template, using pairs */
  int mp, mm, wp, wm, awp, awm, amp, amm ; /* mutant, wild-type or donor, acceptor. p/m: read in forward or reverse direction */
  int topRun ;
  float minFrequency ; /* initialized at 1000 */
  float maxFrequency ; /* initialized at -1000 */
  float alleleFrequency ;

  int tested, rejected, notMeasurable, wellCovered ; /* prevalence */
  int pure, high, mid, low, ref ;
  int xpure, xhigh, xmid, xlow, xref ;
  int ppure, phigh, pmid, plow, pref ;
  float pureC, highC, midC, lowC, refC ;
  float xpureC, xhighC, xmidC, xlowC, xrefC ;
  float ppureC, phighC, pmidC, plowC, prefC ;
  int wiggleProblemHigh ;
  int wiggleProblemLow ;
  Array counts ;
} SNP ;

typedef struct runStruct {
  int run ; /* self, needed after we sort the table */
  int pure, high, mid, low, ref, tested, rejected, notMeasurable, wellCovered ; /* prevalence */
  int pureC, highC, midC, lowC, refC, testedC, rejectedC, notMeasurableC, wellCoveredC ; /* prevalence */
  int xpure, xhigh, xmid, xlow, xref, xtested, xrejected, xnotMeasurable, xwellCovered ; /* prevalence */
  int ppure, phigh, pmid, plow, pref, ptested, prejected, pnotMeasurable, pwellCovered ; /* prevalence */
  int ppureC, phighC, pmidC, plowC, prefC, ptestedC, prejectedC, pnotMeasurableC, pwellCoveredC ; /* prevalence */
  int *types ;
  int *typesR ;
  int *typesC ;
} RC ;

/*************************************************************************************/
/*************************************************************************************/
 
static int tctSnpA1Order  (const void *va, const void *vb)
{
  const SNP *up = (const SNP *)va, *vp = (const SNP*)vb ;
  int n ;
  
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;

  return 0 ;
} /* tctSnpA1Order */

/*************************************************************************************/
#define MAXTYPE 6
/* TYPE up->counts(run *MAXTYPE + i)  ::  0/1: var-donor/accp, 2/3 :ref d/a, 4/5: wiggle d/a */
static int tctSnpParseOne (TCT *tct, ACEIN ai)
{
  int nn = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  DICT *selectDict = tct->selectDict ;
  DICT *runDict = tct->runDict ;
  DICT *typeDict = dictHandleCreate (8, h) ;

  dictAdd (typeDict, "Var", 0) ;
  dictAdd (typeDict, "VarDonor", 0) ;
  dictAdd (typeDict, "VarAccep", 0) ;

  dictAdd (typeDict, "Ref", 0) ;
  dictAdd (typeDict, "RefDonor", 0) ;
  dictAdd (typeDict, "RefAccep", 0) ;

  while (aceInCard (ai))
    {
      char cutter, *cq, *cp = aceInWordCut (ai, "\t", &cutter) ;
      const char *ccp ;
      int snp, type, target, run, delta, subDelIns = 0 ;
      int a0, a1, a2, a3, da = 0 ;
      SNP *up ;

      if (! cp || *cp == '#' || *cp == '/')
	continue ;
    
      if (0)
	{ /* the tsf file has no zone */
	  if (tct->zone && strcasecmp (cp, tct->zone))
	    continue ;
	  
	  aceInStep (ai,'\t') ;
	  cp = aceInWordCut (ai, "\t", &cutter) ;
	  if (! cp || *cp == '#' || *cp == '/')
	    continue ;
	} 
      if (0)
	{
	  cq = strstr (cp, "__") ;
	  if (!cq)
	    continue ;
	  *cq = 0 ; cq += 2 ;
	  if (! dictFind (typeDict, cp, &type))
	    continue ;
	  cp = cq ;
	  cq = strstr (cp, "___") ;
	  if (cq)
	    { *cq = 0 ; cq += 2 ; }
	}
      if (selectDict && ! dictFind (selectDict, cp, 0))
	continue ;
      if (! strcmp (cp, "NC_045512:9693:DelIns_118_119::"))  /* HACK 2020_06_10, bad variant */
	continue ;

      if (dictAdd (tct->snpDict, cp, &snp))
	nn++ ; /* else it is already detected */
      da = subDelIns == 0 ;

      a1 = a2 = target = 0 ;
      if (1)
	{
	  char *cr, *cs ;
	  cr = strchr (cp, ':') ;
	  if (cr)
	    {
	      *cr = 0 ;
	      if (*cp)
		dictAdd (tct->targetDict, cp, &target) ;
	      *cr++ = ':' ;
	      cs = strchr (cr, ':') ;
	      if (cs)
		{
		  *cs = 0 ;
		  sscanf (cr, "%d:", &a1) ;
		  *cs = ':' ;
		}
	      a2 = a1 + da ; /* must be modified for long deletions, ok for sub and inserts */
	    }
	}

      ccp = dictName (tct->snpDict, snp) ;
      if (strstr (ccp, ":Sub:"))
	{ 
	  subDelIns = 1 ;
	  da = 0 ;
	}
      else if (strstr (ccp, ":Sub_"))  /* bases a1 to a2 = a1 + da are modified */
	{
	  int k = 1 ; char cc = 0 ;
	  const char *ccq, *ccr ;

	  subDelIns = 2 ; da = 1 ; /* wild guess */
	  ccq = strstr (ccp, ":Sub_") + 5 ; 
	  if (sscanf (ccq, "%d%c", &k, &cc) == 2 && cc == ':')
	    da = k - 1 ;
	  else
	    {
	      ccr = strchr(ccq, ':') ; 
	      if (ccr) 
		da = ccr - ccq ; 
	    }
	}
      else if (strstr (ccp, ":Del:"))
	{ subDelIns = 3 ; da = 2 ; }  /* bases a1 and a2 = a1 + da are correct */
      else if (strstr (ccp, ":Del_"))
	{ 
	  int k = 1 ; char cc = 0 ;
	  const char *ccq, *ccr ;

	  subDelIns = 4 ; ccq = strstr (ccp, ":Del_") + 5 ; 
	  if (sscanf (ccq, "%d%c", &k, &cc) == 2 && cc == ':')
	    da = k+1 ;
	  else
	    {
	      ccr = strchr(ccq, ':') ; 
	      if (ccr) 
		da = ccr - ccq + 1 ; 
	      else 
		subDelIns = 7 ; 
	    }
	}
      else if (strstr (ccp, ":Ins:"))
	{ subDelIns = 5 ; da = 1 ; }  /* bases a1 and a2 = a1 + 1 are correct */
      else if (strstr (ccp, ":Ins_"))
	{ subDelIns = 6 ; da = 1 ; }
      else if (strstr (ccp, ":DelIns_"))
	{
	  int k = 1 ; char cc = 0 ;
	  const char *ccq, *ccr ;

	  subDelIns = 7 ; da = 1 ; /* wild guess */
	  ccq = strstr (ccp, ":DelIns_") + 8 ; 
	  if (sscanf (ccq, "%d%c", &k, &cc) == 2 && cc == '_')
	    da = k + 1 ;
	  else
	    {
	      ccr = strchr(cq, '_') ; 
	      if (ccr) 
		da = ccr - ccq + 1 ; 
	    }
	}
      a2 = a1 + da ;

      up = arrayp (tct->snps, snp, SNP) ;
      up->snp = snp ; /* self */
      up->a1 = a1 ;
      up->a2 = a2 ;
      up->target = target ;
      if (! up->counts)
	{
	  up->counts = arrayCreate (MAXTYPE * (1+dictMax (runDict)), KEY) ;

	  up->minFrequency = 1000 ;
	  up->maxFrequency = -1000 ;
	}
      
      aceInStep (ai, '\t') ;
      cp = aceInWordCut (ai, "\t", &cutter) ;
      
      if (! cp || *cp == '#' || *cp == '/')
	continue ;
      if (tct->mergeCounts)
	dictAdd (runDict, cp, &(run)) ;
      if (! dictFind (runDict, cp, &(run)))
	continue ;
      
      aceInStep (ai, '\t') ;  /* the format, expect iiii or ii */
      cp = aceInWordCut (ai, "\t", &cutter) ; /* the format of this tsf file, ignore */
      if (! cp || (strcmp(cp,"10") && strcmp(cp,"10i") && strcmp(cp,"iiiiiiiiii")))
	continue ;
      
      aceInStep (ai,'\t') ;
      aceInInt (ai, &a0) ;   /* forward counts */
      aceInStep (ai,'\t') ;
      aceInInt (ai, &a1) ;  /* reverse counts (direction of the read, not of the template as in pair sequencing */
      /* in the next 2 coulumns, we have access to the stand coverage */
      aceInStep (ai,'\t') ;
      aceInInt (ai, &a2) ;   /* forward counts */
      aceInStep (ai,'\t') ;
      aceInInt (ai, &a3) ;  /* reverse counts (direction of the read, not of the template as in pair sequencing */

      switch (type)
	{
	case 1: /* Variant */
	  up->type |= 0x1 ;   /* split var counts */
	  keySet (up->counts, MAXTYPE * run + 1) = a0 + a1 ; /* double register */
	  up->acoverp += a2 ;
	  up->acoverm += a3 ;
	case 2: /* VarDonor*/
	  delta = 0 ;
	  up->coverp += a2 ;
	  up->coverm += a3 ;
	  break ;
	case 3: /* VarAccep */
	  delta = 1 ;
	  break ;

	case 4: /* Ref */
	  up->type |= 0x2 ;   /* split var counts */
	  keySet (up->counts, MAXTYPE * run + 3) = a0 + a1 ; /* double register */
	  up->acoverp += a2 ;
	  up->acoverm += a3 ;
	case 5: /* RefDonor */
	  delta = 2 ;
	  up->coverp += a2 ;
	  up->coverm += a3 ;
	  break ;
	case 6: /* RefAccep */
	  delta = 3 ;
	  up->acoverp += a2 ;
	  up->acoverm += a3 ;
	  break ;
	}
      up->type = (up->type & 0x3) +  8 * subDelIns ;
      
      keySet (up->counts, MAXTYPE * run + delta) = a0 + a1 ; /* merge the strands */

      switch (type)
	{
	case 1: /* Variant */
	  up->amp += a0 ;
	  up->amm += a1 ;
	case 2: /* VarDonor*/
	  up->mp += a0 ;
	  up->mm += a1 ;
	  break ;
	case 3: /* VarAccep */
	  up->amp += a0 ;
	  up->amm += a1 ;
	  break ;
	case 4: /* Ref */
	  up->awp += a0 ;
	  up->awm += a1 ;
	case 5: /* RefDonor */
	  up->wp += a0 ;
	  up->wm += a1 ;
	  break ;
	case 6: /* RefAccep */
	  up->awp += a0 ;
	  up->awm += a1 ;
	  break ;
	}
    }

  ac_free (h) ;
  return nn ;
} /* tctSnpParseOne */

/*************************************************************************************/
/*************************************************************************************/

/* import the runs metadata from the database */

static int tctGetSelectedSnps (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  tbl = 0 ;
  const char *errors = 0 ;
  const char *ccp ;
  char *cp, *cq, *sep ;
  int ir, nn = 0, target ;
  DICT *targetDict = tct->targetDict ;
  DICT *snpDict = tct->snpDict ;
  vTXT txt = vtxtHandleCreate (h) ;

  if (tct->force)
    {
      cp = strnew (tct->force, h) ;
      sep = " where (" ;
      while (cp)
	{
	  cq = strchr (cp, ',') ;
	  if (cq) *cq = 0 ;
	  vtxtPrintf (txt, "%s s#%s ", sep, cp) ;
	  sep = " or " ;
	  if (cq) 
	    cp = cq + 1 ;
	  else
	    cp = 0 ;
	}
      vtxtPrintf (txt, " ) ") ;
    }

  tbl = ac_bql_table (tct->db, hprintf (h, "select s, c, a1, a2 from s in ?variant %s , c in s->intmap, a1 in c[1], a2 in c[2] where a2", vtxtPtr (txt)), 0, 0, &errors, tct->h) ;

  if (tbl)
    for (ir = 1 ; ir < tbl->rows ; ir++)
      {
	SNP *up ;
	int snp ;
	int a1 = ac_table_int (tbl, ir, 2, 0) ;
	int a2 = ac_table_int (tbl, ir, 3, 0) ;
	ccp = ac_table_printable (tbl, ir, 0, 0) ;
	dictAdd (snpDict, ccp, &snp) ;
	ccp = ac_table_printable (tbl, ir, 1, 0) ;
	dictAdd (targetDict, ccp, &target) ;
	up = arrayp (tct->snps, snp, SNP) ;
	up->snp = snp ;
	up->target = target ;
	up->a1 = a1 ;
	up->a2 = a2 ;
	up->select = TRUE ;
	nn++ ;
      }

  ac_free (h) ;
  return nn ;
} /* tctGetSelectedSnps */

/*************************************************************************************/

static int hisMutIsMyMut (SNP *up, SNP *vp)
{
  if (up == vp)
    return 0 ;
  if (1) /* his sides cover my mut A, so he is +B or ++, and his +B and ++ counts reinforce my + count */
    {
      int a1 = up->a1 ;
      int a2 = up->a2 ;
      int b1 = vp->a1 ;
      int b2 = vp->a2 ;
      int aT = up->type / 8 ;
      int bT = vp->type / 8 ;

      if (aT == 1)
	{
	  switch (bT)
	    {
	    case 2: /* my single sub is inide a muti sub, they are most likely identical multi sub */
	      if (a1 >= b1 && a1 <= b2)
		return 0 ;
	    case 1:
	      if (a1 <= b1 + 14 && a2 >= b2 - 14)
		return 1 ; /* his mut is my ref */
	      break ;
	    case 3:
	    case 4:
	    case 5:
	    case 6:
	    case 7: /* base vp->a1 and vp->a2 are correct */
	      if (bT < 7  || b2 - b1 < 14)
		{
		  if (a2 <= b1 && a2 >= b1 - 14)
		    return 1 ; /* his mut is my ref */
		  if (a1 <= b2 + 14 && a1 >= b2)
		    return 1 ; /* his mut is my ref */
		}
	      break ;
	    }
	}
      if (aT == 2)
	{
	  switch (bT)
	    {
	    case 2: /* my single sub is inide a muti sub, they are most likely identical multi sub */
	    case 1:
	    case 3:
	    case 4:
	    case 5:
	    case 6:
	      if (b2 - b1 < 14 && 
		  (
		   (a1 > b2 && a2 < b2 + 14 ) ||
		   (a2 < b1 && a1 > b1 - 14 ) 
		   )
		  )
		return 1 ; /* his mut is my ref */
	      break ;
	    case 7:
	      break ;
	    }
	}
      if (aT >=3 && aT <= 7 && a2 - a1 < 14)
	{
	  switch (bT)
	    {
	    case 1:
	    case 2:
	    case 5:
	    case 6:
	    case 7: /* base vp->a1 and vp->a2 are correct */
	      if (a1 >= b2 - 14 && a2 <= b1 + 14)
		return 1 ; /* his mut is my ref */
	      break ;
	    case 3:
	    case 4:
	      if (a2 >= b1 - 14 && a2 <= b1)
		return 1 ; /* his mut is my ref */
	      if (a1 >= b2 && a1 <= b2 + 14)
		return 1 ; /* his mut is my ref */
	      break ;
	    }
	}
    }

  return 0 ; /* remainder, for example a complex multisub */
}

/*************************************************************************************/

static void tctSnpParse (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;

  tct->targetDict = dictHandleCreate (64, tct->h) ;
  tct->snpDict = dictHandleCreate (1000, tct->h) ;
  if (tct->inFileOfFileList)
    {
      vTXT txt = 0 ;
      char *sep = "" ;
      ACEIN ai = aceInCreate (tct->inFileOfFileList, 0, h) ;
      if (!ai)
	messcrash ("Cannot open -file_of_file_list %s\n", tct->inFileOfFileList) ;

      txt = vtxtHandleCreate (h) ;
      while (aceInCard (ai))
	{
	  char cutter, *cp ;
	  while ((cp = aceInWordCut (ai, " ,\t", &cutter)))
	    {
	      if (cp && *cp && *cp != '#')
		{
		  vtxtPrintf (txt, "%s%s",sep, cp) ;
		  sep = "," ;
		}
	    }
	}
      if (*sep)
	tct->inFileList = vtxtPtr (txt) ;
    }

  tct->runs = arrayHandleCreate (300, RC, tct->h) ;
  tct->snps = arrayHandleCreate (30000, SNP, tct->h) ;
  if (tct->inFileList)
    {
      ACEIN ai = 0 ;
      int nn = 0 ;
      const char *error = 0 ;

      /* check that all files exist */
      if (! aceInCheckList (tct->inFileList, &error, h))
	messcrash ("Bad parameter -f, %s"
		       ,  error
		       ) ;	
 
      /* parse a validated list of files */
      while ((ai = aceInCreateFromList (ai, &nn, tct->inFileList, tct->gzi, h)))
	{
	  tct->snpDetected += tctSnpParseOne (tct, ai) ;
	}
     }
  if (tct->force) tctGetSelectedSnps (tct) ;
  arraySort (tct->snps, tctSnpA1Order) ;

  /* coallesce the neighbours */
  if (arrayMax (tct->snps) > 1)
    {
      int snp, snpMax = arrayMax (tct->snps) ;
      SNP *up, *vp ;
      int ir, irMax = dictMax (tct->runDict) ;
      
      for (snp = 1, up = arrp (tct->snps, 0, SNP) ; snp < snpMax ; up++, snp++)
	{
	  if (up->counts) /* count problems */
	    for (ir = 1 ; ir <= irMax ;ir++)
	      {
		int s ;
		int a1 = up->a1 ;
		int a2 = up->a2 ;
		KEY *kp = arrayp (up->counts, ir * MAXTYPE, KEY) ;
		int m5 = kp[0] ; /* mutant */
		int m3 = kp[1] ; /* mutant */
		int r5 = kp[2] ; /* ref */
		int r3 = kp[3] ; /* ref */
		
		for (s = snp - 30, vp = up - 30 ; s < snp + 30 && s < snpMax ; s++, vp++)
		  {
		    if (s <= 0)
		      continue ;
		    if (vp->counts && vp->a1 >= a1 - 14 && vp->a1 <= a1 +14)
		      {
			KEY *kp2 = arrayp (vp->counts, ir * MAXTYPE, KEY) ;
			int m52 = kp2[0] ; /* variant, including self */
			int r52 = kp2[2] ; /* variant, including self */
			switch (hisMutIsMyMut (up, vp))
			  {
			  case 2:
			    m5 += m52 ;
			    break ;
			  case 1:
			    if (r5 < m52 + r52) r5 = m52 + r52 ; /* all the counts of my neighbous are my wild type, unless they are my mutant */
			    break ;
			  }
		      }
		    if (vp->counts && vp->a1 >= a2 - 14 && vp->a2 <= a1 +14)
		      {
			KEY *kp2 = arrayp (vp->counts, ir * MAXTYPE, KEY) ;
			int m32 = kp2[1] ; /* variant, including self */
			int r32 = kp2[3] ; /* variant, including self */
			switch (hisMutIsMyMut (up, vp))
			  {
			  case 2:
			    m3 += m32 ;
			    break ; 
			  case 1:
			    if (r3 < m32 + r32) r3 = m32 + r32 ; /* all the counts of my neighbous are my wild type, unless they are my mutant */
			    break ;
			  }
		      }
		  }
		kp[4] = (m5 > r5 ? 2*m5 : 2*r5+1) ;  /* union with future wiggle value */
		kp[5] = (m3 > r3 ? 2*m3 : 2*r3+1) ;  /* if the neighbour is wild, i am wild */
	      }
	}
    }
  /* accept the coallesced values */
  if (arrayMax (tct->snps) > 1)
    {
      int snp, snpMax = arrayMax (tct->snps) ;
      SNP *up ;
      int ir, irMax = dictMax (tct->runDict) ;
      
      for (snp = 1, up = arrp (tct->snps, 0, SNP) ; snp < snpMax ; up++, snp++)
	{
	  if (up->counts) /* count problems */
	    for (ir = 1 ; ir <= irMax ;ir++)
	      {
		KEY *kp = arrayp (up->counts, ir * MAXTYPE, KEY) ;
		int m5 = kp[0] ; /* mutant */
		int m3 = kp[1] ; /* mutant */
		int r5 = kp[2] ; /* ref */
		int r3 = kp[3] ; /* ref */
		int w5 = kp[4] ; /* ref */ if (w5 & 0x1) w5 = - (w5-1)/2 ; else w5 = w5/2 ;
		int w3 = kp[5] ; /* ref */ if (w3 & 0x1) w3 = - (w3-1)/2 ; else w3 = w5/2 ;
		if (w5 > m5) { kp[0] = w5 ;}
		if (-w5 > r5) { kp[2] = -w5 ;}
		if (w3 > m3) { kp[1] = w3 ;}
		if (-w3 > r3) { kp[3] = -w3 ;}
		kp[4] = m3 + r3 ; /* corrected estimation of the donor site coverage */
		kp[5] = m5 + r5 ; /* acceptor = mutant-variant + wild-reference + others */
	      }
	}
    }
		
  ac_free (h) ;

  return  ;
} /* tctSnpParse */

/*************************************************************************************/
/*************************************************************************************/

static void tctWiggleParseOne (TCT *tct, int run, int target)
{
  AC_HANDLE h = ac_new_handle () ;
  char *fNam = hprintf (h, "%s/%s/%s/R.chrom.frns.u.BF.gz", tct->wiggleDir, dictName (tct->runDict, run), dictName (tct->targetDict, target)) ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;
  int line = 0, start = 0, step = 0, x = 0 ;
  float w ;
  Array ww = 0 ;
  char *cp ;

  if (ai)
    { /* parse the whole file, it is simpler and unavoidable anyway. We get a tmp sorted wiggle, we can then access it from the unsorted snp array */
      line++ ;
      ww = arrayHandleCreate (1000000, float, h) ;
      while (aceInCard (ai))
	{
	  if (!step)
	    {
	      if (line > 12) break ; /* badly formatted file */
	      cp = aceInWord (ai) ;
	      if (!cp || strcmp (cp, "fixedStep"))
		continue ;
	      cp = aceInWord (ai) ;
	      if (!cp)
		continue ;
	      cp = strstr (cp,"chrom=") ;
	      if (!cp)
		continue ;
	      cp += 6 ;
	      if (!cp || strcmp (cp, dictName (tct->targetDict, target)))
		continue ;
	      
	      cp = aceInWord (ai) ;
	      if (!cp)
		continue ;
	      cp = strstr (cp,"start=") ;
	      if (!cp)
		continue ;
	      cp += 6 ;
	      if (sscanf (cp, "%d", &start) < 1 || start<= 0)
		continue ;
	      
	      cp = aceInWord (ai) ;
	      if (!cp)
		continue ;
	      cp = strstr (cp,"step=") ;
	      if (!cp)
		continue ;
	      cp += 5 ;
	      if (sscanf (cp, "%d", &step) < 1 || step<= 0)
		continue ;
	      continue ;
	    }
	  x++ ; w = 0 ;
	  if (aceInFloat (ai, &w))
	    array (ww, x, float) = w ;	  
	}
    }

  if (ww && arrayMax (ww)) /* we have a wiggle, with known start and step, we scan the snps and update the wiggle values at a1 and a2 */
    {
      int kk, pass, snp, snpMax = arrayMax (tct->snps) ;
      int zz[2] ;
      int ww0[2] ;
      SNP *up ;
      int stop = start + step * arrayMax (ww) ;
      float wwww[5], w0, dw, dw0, z ;
      if (!step) step = 1 ;
      for (snp = 1, up = arrp (tct->snps, 0, SNP) ; snp < snpMax ; up++, snp++)
	{
	  KEY *kp = 0 ;
	  if (! up->select && tct->minSnpFrequency && (up->maxFrequency < 0 || up->maxFrequency < tct->minSnpFrequency))
	    continue ;
	  if (up->a1 < start+10 || up->a1 > stop - 10)
	    continue ;
	  if (! up->counts)
	    {
	      up->counts = arrayCreate (MAXTYPE * (1+dictMax (tct->runDict)), KEY) ;
	      
	      up->minFrequency = 1000 ;
	      up->maxFrequency = -1000 ;
	    }
	  kp = arrayp (up->counts, run * MAXTYPE, KEY) ; 
	  memset (wwww, 0, sizeof(wwww)) ;
	  ww0[0] = kp[4] ; /* was set as sum of negighbours */
	  ww0[1] = kp[5] ; /* was set as sum of neighbours */
	  for (kk = 0 ; kk < 2 ; kk++) /* donor acceptor */
	    {
	      int i, a1, a01 = kk == 0 ? up->a1 : up->a2 ;
	      if (! kk || up->a1 >= up->a2 + 3 || up->a1 <= up->a2 - 3)  /* reuse the same wwww */
		{
		  memset (wwww, 0, sizeof(wwww)) ;
		  for (pass = 0 ; pass < 5 ; pass++)
		    {
		      switch (pass)
			{
			case 0: a1 = a01 ; break ;
			case 1: a1 = a01 - 1 ; break ;
			case 2: a1 = a01 + 1 ; break ;
			case 3: a1 = a01 - 10 ; break ;
			case 4: a1 = a01 + 10 ; break ;
			}
		      if (a1)
			{
			  int b1 = (a1 - start)/step ;
			  float db = a1 - start - step * b1 ; db /= 1.0 * step ; /* how much we need of next position */
			  if (db >= 1) 
			    messcrash ("a1=%d b1=%d step=%d start=%d db=%f should be < 1, file %s", up->a1, step*b1, step, start, db, fNam) ;
			  if (b1 < arrayMax (ww) - 1)
			    {
			      z = (1 - db) * arr (ww, b1, float) + db * arr (ww, b1+1, float) ;
			      wwww[pass] = z ;
			    }
			}
		    }
		}
	      /* select best estimate, closest to the word counts */
	      if (! up->counts)
		up->counts = arrayCreate (MAXTYPE * (1+dictMax (tct->runDict)), KEY) ;
	      w0 = ww0[kk] ;
	      dw0 = 1000000000 ; /* horrible value */
	      z = 0 ;
	      for (i = 0 ; i < 5 ; i++) /* try direct positions */
		{
		  dw = wwww[i] - w0 ; if (dw < 0) dw = -dw ; if (dw0 > dw) { dw0 = dw ; z = wwww[i] ; }
		}
	      for (i = 1 ; i < 4 ; i += 2) /* try averages */
		{
		  dw = (wwww[i] + wwww[i+1])/2 - w0 ; if (dw < 0) dw = -dw ; if (dw0 > dw) { dw0 = dw ; z = (wwww[i] + wwww[i+1])/2 ; }
		}
	      zz[kk] = z ;
	    }
	  z = zz[0] + zz[1] ;
	  w0 = ww0[0] + ww0[1] ;
	  if (z > 20 && z > 2 * w0) 
	    { zz[0] = - zz[0] ; up->wiggleProblemHigh++ ; }
	  if (w0 > 20 && 2 * z < w0)
	    { zz[0] = w0/2 ; zz[1] = -z/2 ;  up->wiggleProblemLow++ ; }
	  kp[4] = zz[0] ; /* register best estimate and replace count gathering */
	  kp[5] = zz[1] ; /* register best estimate and replace count gathering */
	}
    }
  
  ac_free (h) ;
  return  ;
} /* tctWiggleParseOne */

/*************************************************************************************/

static void tctWiggleParse (TCT *tct)
{
  int run, runMax = dictMax (tct->runDict) ;
  int target, targetMax = dictMax (tct->targetDict) ;

  if (! tct->wiggleDir)
    return ; 

  tct->wdelta = keySetHandleCreate (tct->h) ;
  for (run = 1 ; run <= runMax ; run++)
    for (target = 1 ; target <= targetMax ; target++)
      tctWiggleParseOne (tct, run, target) ;

  if (keySetMax (tct->wdelta))
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (tct->outFileName, ".snp_wiggle_delta.txt", 0, h) ;
      int i, iMax = keySetMax (tct->wdelta) ;

      aceOutDate (ao, "###", "Histo of delta frequency relative to the wiggles") ;
      aceOutf (ao, " ## delta = 100 * w / (v + r), where w is the wiggle coverage and v and r and the reported numbers of 31-mers representing at the same position the reference (r) and the variant (v)\n") ;
      aceOut (ao, "# delta\tNumber of occurences\n") ;
      for (i = 0 ; i <= iMax ; i++)
	aceOutf (ao, "%d\t%d\n", i, keySet (tct->wdelta, i)) ;
      ac_free (h) ;
    }

  return  ;
} /* tctWiggleParse */


/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

typedef struct remapStruct { int target_class, gene, target, x1, x2, chrom, a1, a2 ; } RMP ;

/*****************/

static int atlasOrder1 (const void *va, const void *vb)
{
  const RMP *a = (const RMP *)va, *b = (const RMP *)vb ;
  int n ;
  n = a->target_class - b->target_class ; if (n) return n ;
  n = a->target - b->target ; if (n) return n ;
  n = a->x1 - b->x1 ; if (n) return n ;
  n = a->x2 - b->x2 ; if (n) return n ;
  n = a->chrom - b->chrom ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;

  return 0 ;
} /* atlasOrder1 */

/*****************/

static int atlasOrder2 (const void *va, const void *vb)
{
  const RMP *a = (const RMP *)va, *b = (const RMP *)vb ;
  int n ;
  n = a->chrom - b->chrom ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  n = a->target_class - b->target_class ; if (n) return n ;
  n = a->target - b->target ; if (n) return n ;
  n = a->x1 - b->x1 ; if (n) return n ;
  n = a->x2 - b->x2 ; if (n) return n ;

  return 0 ;
} /* atlasOrder1 */

/*************************************************************************************/
/* a tampon in bp, to allocate the snp to the promotor region of the gene */
#define PROMOTOR_DELTA 500

static int tctCreateAtlas (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai ;
  Array atlas = 0, geneAtlas = 0, map, geneMap ;
  const char *ccp ;
  RMP *up ;
  int nn = 0, x1, x2, a1, a2, target, chrom, gene, target_class ;

  x1 = x2 = a1 = a2 = 0 ;
  if (! tct->target_classDict)
    {
      tct->targetDict = dictHandleCreate (4, tct->h) ;
      tct->target_classDict = dictHandleCreate (4, tct->h) ;
      tct->chromDict = dictHandleCreate (64, tct->h) ;
      tct->target2geneAtlas = arrayHandleCreate (4, Array, tct->h) ;
      tct->target2exonAtlas = arrayHandleCreate (4, Array, tct->h) ;
    }

  if (tct->remap2genome)
    ai = aceInCreate (tct->remap2genome, FALSE, h) ;
  else
    ai = aceInCreate (tct->remap2genes, FALSE, h) ;
  aceInSpecial (ai, "\t\n") ;

  while (ai && aceInCard (ai))
    {
      ccp = aceInWord (ai) ; /* target_class */
      if (! ccp || *ccp == '#')
	continue ;
      if (tct->target_class && strcmp (ccp, (tct->target_class)))
	continue ;
      ccp = "any" ; /* ignore target class for the moment */
      dictAdd (tct->target_classDict, ccp, &target_class) ;
      atlas = array (tct->target2exonAtlas, target_class, Array) ;
      geneAtlas = array (tct->target2geneAtlas, target_class, Array) ;
      if (! atlas)
	{
	  atlas = arrayHandleCreate (100000, Array, tct->h) ;
	  array (tct->target2exonAtlas, target_class, Array) = atlas ;
	}
      if (! geneAtlas)
	{
	  geneAtlas = arrayHandleCreate (10000, Array, tct->h) ;
	  array (tct->target2geneAtlas, target_class, Array) = geneAtlas ;
	}
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* target */
      if (! ccp || *ccp == '#')
	continue ; 
      dictAdd (tct->targetDict, ccp, &target) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x1))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x2))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* chrom */
      if (! ccp || *ccp == '#')
	continue ; 
      dictAdd (tct->chromDict, ccp, &chrom) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a1))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a2))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* chrom */
      if (! ccp || *ccp == '#')
	continue ;
      dictAdd (tct->targetDict, ccp, &gene) ;
      if (1)
	{
	  map = array (atlas,  tct->remap2genome ? target : chrom, Array) ;
	  if (! map)
	    map = array (atlas, tct->remap2genome ? target : chrom, Array) = arrayHandleCreate (8, RMP, tct->h) ;
	  geneMap = array (geneAtlas,  tct->remap2genome ? target : chrom, Array) ;
	  if (! geneMap)
	    geneMap = array (geneAtlas, tct->remap2genome ? gene : chrom, Array) = arrayHandleCreate (8, RMP, tct->h) ;
	  nn++ ;
	  up = arrayp (map, arrayMax (map), RMP) ;
	  up->target = target ;
	  up->chrom = chrom ;
	  if (tct->remap2genome || a1 < a2)
	    {
	      up->x1 = x1 ; up->x2 = x2 ;
	      up->a1 = a1 ; up->a2 = a2 ;
	    }
	  else
	    {
	      up->x1 = x2 ; up->x2 = x1 ;
	      up->a1 = a2 ; up->a2 = a1 ;
	    }

	  up =  arrayMax (geneMap) ? arrayp (geneMap, arrayMax (geneMap) - 1, RMP) : 0 ;
	  if (! up || up->target != gene)
	    up = arrayp (geneMap, arrayMax (geneMap), RMP) ;
	  up->target = gene ;
	  up->chrom = chrom ;
	  if (a1 < a2)
	    {
	      a1 -= PROMOTOR_DELTA ;
	      if (! up->a1 || a1 < up->a1) up->a1 = a1 ; 
	      if (a2 > up->a2) up->a2 = a2 ;
	      up->x1 = 1 ; up->x2 = up->a2 - up->a1 + 1 ;
	    }
	  else
	    {
	      a1 += PROMOTOR_DELTA ;
	      if (! up->a1 || a2 < up->a1) up->a1 = a2 ; 
	      if (a1 > up->a2) up->a2 = a1 ;
	      up->x2 = 1 ; up->x1 = up->a2 - up->a1 + 1 ;
	    }
	}
    }

  for (a1 = 0 ; atlas && a1 < arrayMax (atlas) ; a1++)
    {
      map = array (atlas, a1, Array) ;
      if (map)
	arraySort (map, tct->remap2genome ? atlasOrder1 : atlasOrder2) ;
    }
  for (a1 = 0 ; geneAtlas && a1 < arrayMax (geneAtlas) ; a1++)
    {
      map = array (geneAtlas, a1, Array) ;
      if (map)
	{
	  arraySort (map, tct->remap2genome ? atlasOrder1 : atlasOrder2) ;
	  arrayCompress (map) ;
	}
    }
   
  fprintf (stderr, "tctCreateAtlas found %d exons in file %s\n"
	   , nn
	   , tct->remap2genome ? tct->remap2genome :  tct->remap2genes
	   ) ;

  ac_free (h) ;
  return nn ;
} /* tctCreateAtlas */

/*************************************************************************************/
/* 
 * We import all the relevant data from the ZZ database
 * Remap the transcript variants into genome coordinates
 */

static BOOL tctRemap1Do (TCT *tct, int mrna, int x1, int x2, int *chromp, int *a1p, int *a2p, int *strandp)
{
  int ii, target_class = 0 ;
  RMP *up ;
  Array atlas, map ;
  
  dictAdd (tct->target_classDict, "any", &target_class) ;
  if (tct->target2exonAtlas && 
      target_class < arrayMax (tct->target2exonAtlas) && 
      (atlas =  arr (tct->target2exonAtlas, target_class, Array)) &&
      mrna < arrayMax(atlas) &&
      (map = array (atlas, mrna, Array))
      )
    {
      for (ii = 0, up = arrp (map, 0, RMP) ; ii < arrayMax (map) ; ii++, up++)
	{
	  if (up->x1 <= x1 && up->x2 >= x1)
	    {
	      if (up->a1 < up->a2)
		{ *a1p = up->a1 + x1 - up->x1 ; *a2p = up->a1 + x2 - up->x1 ; *strandp = 1 ; *chromp = up->chrom ; return TRUE ; }
	      else
		{ *a1p = up->a1 - x1 + up->x1 ; *a2p = up->a1 - x2 + up->x1 ; *strandp = - 1 ; *chromp = up->chrom ; return TRUE ; }
	    }
	  else if (x1 < up->x1)
	    break ;
	}
    }
  return FALSE ;
} /* tctRemap1Do */

/*************************************************************************************/
/* scan the VariantDB acedb database
 * add the remap info 
 * Remap the transcript variants into genome coordinates
 */
static int tctRemap1 (TCT *tct)
{
  AC_HANDLE  h1 = 0, h = ac_new_handle () ;
  AC_ITER iter ;
  AC_OBJ variant = 0 ;
  AC_TABLE mrnaTable = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  int x1, x2, a1, a2, strand, chrom, mrna, nn1 = 0, nn2 = 0 ;
  const char *mm ;
  const char *errors = 0 ;
  
  if (tct->db)
    {
      iter = ac_query_iter (tct->db, TRUE, "find variant mRNA && ! IntMap", 0, h) ;
      while (ac_free (variant), ac_free (h1), variant = ac_iter_obj (iter))
	{
	  h1 = ac_new_handle () ;
	  nn1++ ;
	  mrnaTable = ac_tag_table (variant, "mRNA", h1) ;
	  if (mrnaTable)
	    {
	      a1 = a2 = 0 ;
	      x1 = ac_table_int (mrnaTable, 0, 1, 0) ;
	      x2 = ac_table_int (mrnaTable, 0, 2, 0) ;
	      mm = ac_table_printable (mrnaTable, 0, 0, 0) ;
	      if (mm && dictFind (tct->targetDict, mm, &mrna) && tctRemap1Do (tct, mrna, x1, x2, &chrom, &a1, &a2, &strand))
		{
		  nn2++ ;
		  vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name (variant), h1)) ;
		  vtxtPrintf (txt, "IntMap %s %d %d\n\n", dictName (tct->chromDict, chrom), a1, a2) ; 
		}
	    }
	}
      ac_parse (tct->db, vtxtPtr (txt), &errors, 0, h) ; 
    }
  fprintf(stderr, "tctRemap1 found %d variants remapped %d\n", nn1, nn2) ;

  if (errors && *errors) fprintf(stderr, "tctRemap parsing error %s\n", errors) ;
  ac_free (h1) ;
  ac_free (h) ;
  return nn2 ;
} /* tctRemap1 */

/*************************************************************************************/
/* 
 * We import all the relevant data from the ZZ database
 *   Remap the genome variants into transcript coordinates
 *   Remap to the genebox without giving precise coordinates
 */

static BOOL tctRemap2geneBoxDo (TCT *tct, int chrom, int pos, int *JJp, int *geneBoxp, int *x1p, int *strandp)
{
  int ii, target_class = 0 ;
  RMP *up ;
  Array atlas, map ;
  
  dictAdd (tct->target_classDict, "any", &target_class) ;
  if (tct->target2geneAtlas && 
      target_class < arrayMax (tct->target2geneAtlas) && 
      (atlas =  arr (tct->target2geneAtlas, target_class, Array)) &&
      chrom < arrayMax(atlas) &&
      (map = array (atlas, chrom, Array))
      )
    {
      for (ii = *JJp, up = arrp (map, ii, RMP) ; ii < arrayMax (map) ; ii++, up++)
	{
	  if (up->a1 <= pos && up->a2 >= pos)
	    {
	      if (up->x1 < up->x2)
		{ *x1p = up->x1 + pos - up->a1 - PROMOTOR_DELTA ; *strandp = 1 ; *geneBoxp = up->target ; *JJp = ii + 1 ; return TRUE ; }
	      else
		{ *x1p = up->x1 - pos + up->a1  - PROMOTOR_DELTA ; *strandp = - 1 ; *geneBoxp = up->target ; *JJp = ii + 1 ;  return TRUE ; }
	    }
	  else if (pos < up->a1)
	    break ;
	}
    }
  return FALSE ;
} /* tctRemap2Do */

static BOOL tctRemap2Do (TCT *tct, int chrom, int pos, int *JJp, int *mrnap, int *x1p, int *strandp)
{
  int ii, target_class = 0 ;
  RMP *up ;
  Array atlas, map ;
  
  dictAdd (tct->target_classDict, "any", &target_class) ;
  if (tct->target2exonAtlas && 
      target_class < arrayMax (tct->target2exonAtlas) && 
      (atlas =  arr (tct->target2exonAtlas, target_class, Array)) &&
      chrom < arrayMax(atlas) &&
      (map = array (atlas, chrom, Array))
      )
    {
      for (ii = *JJp, up = arrp (map, ii, RMP) ; ii < arrayMax (map) ; ii++, up++)
	{
	  if (up->a1 <= pos && up->a2 >= pos)
	    {
	      if (up->x1 < up->x2)
		{ *x1p = up->x1 + pos - up->a1 ; *strandp = 1 ; *mrnap = up->target ; *JJp = ii + 1 ; return TRUE ; }
	      else
		{ *x1p = up->x1 - pos + up->a1 ; *strandp = - 1 ; *mrnap = up->target ; *JJp = ii + 1 ;  return TRUE ; }
	    }
	  else if (pos < up->a1)
	    break ;
	}
    }
  return FALSE ;
} /* tctRemap2Do */

/*************************************************************************************/
/* scan the VariantDB acedb database
 * add the remap info 
 *   Remap the genome variants into transcript coordinates\n"
 */
static int tctRemap2 (TCT *tct)
{
  AC_HANDLE  h1 = 0, h = ac_new_handle () ;
  AC_ITER iter ;
  AC_OBJ variant = 0 ;
  AC_TABLE intMapTable = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  int pos, pos2, x1 = 0, strand, chrom, mrna = 0, geneBox = 0, nn1 = 0, nn2 = 0, JJ=0 ;
  const char *chromNam ;
  const char *errors = 0 ;
  
  iter = ac_query_iter (tct->db, TRUE, "find variant IntMap && ! geneBox", 0, h) ;
  while (ac_free (variant), ac_free (h1), variant = ac_iter_obj (iter))
    {
      h1 = ac_new_handle () ;
      nn1++ ;
      intMapTable = ac_tag_table (variant, "IntMap", h1) ;
      if (intMapTable)
	{
	  chromNam = ac_table_printable (intMapTable, 0, 0, 0) ;
	  pos = ac_table_int (intMapTable, 0, 1, 0) ;
	  pos2 = ac_table_int (intMapTable, 0, 2, 0) ;
	  if (pos2 == 1 || pos2 > pos) pos2 = 1 ;
	  else if (pos2 == -1 || pos2 < pos) pos2 = 1 ;
	  
	  if (pos && dictFind (tct->chromDict, chromNam, &chrom))
	    {
	      JJ = 0 ;
	      if (! ac_has_tag (variant, "mRNA"))
		while (tctRemap2Do (tct, chrom, pos, &JJ, &mrna, &x1, &strand))
		  {
		    nn2++ ;
		    vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name (variant), h1)) ;
		    vtxtPrintf (txt, "mRNA %s %d %d\n\n", dictName (tct->targetDict, mrna), x1, strand * pos2) ; 
		  }
	      JJ = 0 ;
	      while (tctRemap2geneBoxDo (tct, chrom, pos, &JJ, &geneBox, &x1, &strand))
		{
		  nn2++ ;
		  vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name (variant), h1)) ;
		  vtxtPrintf (txt, "GeneBox %s %d %d\n\n", dictName (tct->targetDict, geneBox), x1, strand * pos2) ; 
		}
	    }
	}
    }
  fprintf(stderr, "tctRemap2 found %d variants remapped %d\n", nn1, nn2) ;
  if (vtxtPtr (txt))
    {
      ACEOUT ao = aceOutCreate ("toto", 0, 0, h) ;
      aceOut (ao,  vtxtPtr (txt)) ;
      ac_free (ao) ;
    }
  ac_parse (tct->db, vtxtPtr (txt), &errors, 0, h) ; 
  if (*errors) fprintf(stderr, "tctRemap parsing error %s\n", errors) ;
  ac_free (h1) ;
  ac_free (h) ;
  return nn2 ;
} /* tctRemap2 */

/*************************************************************************************/
/*************************************************************************************/
/* import the runs metadata from the database */

static int tctGetRunList (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *errors = 0 ;
  AC_TABLE tbl = tct->runMetaData = ac_bql_table (tct->db, hprintf (h, "select r, t, st from p in class project where p == \"%s\", r in p->run where r ISA runs, t in r->title, st in r->Sorting_title", tct->project), 0, "st", &errors, tct->h) ;
  int ir ;
  const char *ccp ;

  if (! tbl)
    messcrash ("No run belong to project:%s in database \n", tct->project ? tct->project : "unspecified") ;
  tct->runDict = dictHandleCreate (tbl->rows + 10, tct->h) ;

  /* ATTENTION: dans la table BQL il ne faut avoir que des data uniques associees au run pour avoir toujours exactement 1 ligne par run, sinon on a une erreur d'attribution des nombres qui eux utilisent runDict pour la memorisation */
  for (ir = 0 ; ir < tbl->rows ; ir++)
    {
      ccp = ac_table_printable (tbl, ir, 0, 0) ;
      if (! ccp || ! dictAdd (tct->runDict, ccp, 0))
	messcrash ("Asynchrony it is necessary to ensure that the dict and the table agree on the numbering of the runs ") ;
      if (strcmp (ccp, dictName(tct->runDict, ir+1)))
	messcrash ("Asynchrony it is necessary to ensure that the dict and the table agree on the numbering of the runs ") ;
    }

  ac_free (h) ;
  return dictMax (tct->runDict) ;
} /* tctGetRunList */

/*************************************************************************************/

static void tctReportFrequency (ACEOUT ao, ACEOUT ao2, TCT *tct, int line, SNP *up, AC_KEYSET ks, BOOL intertwinFrequencyCounts)
{
  DICT *runDict = tct->runDict ;
  int ir, irMax = runDict ? dictMax (runDict) : 0 ;
  AC_TABLE tbl = tct->runMetaData ;
  BOOL hasW = tct->wiggleDir ? TRUE : FALSE ;
  
  if (ao)
    aceOutf (ao, "\t") ;
  if (line) /* titles */
    {
      if (ao)
	{
	  aceOutf (ao, "\tVariant") ;
	  for (ir = 1 ; ir <= irMax ;ir++)
	    {
	      aceOutf (ao, "\t%s", ac_table_printable (tbl, ir - 1, line-1, "")) ;
	      if (intertwinFrequencyCounts)
		aceOutf (ao, "\t%s counting v:variant, r:reference, w:wiggle at 5' and 3' sites if different", ac_table_printable (tbl, ir - 1, line-1, "")) ;
	    }
	}
    }
  else
    {
      if (ao)
	aceOutf (ao, "\t%s", dictName (tct->snpDict, up->snp)) ;
      for (ir = 1 ; ir <= irMax ;ir++)
	{
	  int m5 = 0, m3 = 0, r5 = 0, r3 = 3, w5 = 0, w3 = 0, z3, z5 ;
	  char wiggleProblem = ' ' ;
	  float f5, f3, f ;
	  f3 = f5 = f = -20 ;
	  if (up->counts)
	    {
	      KEY *kp = arrayp (up->counts, MAXTYPE*ir, KEY) ;
	      m5 = kp[0] ;
	      m3 = kp[1] ;
	      r5 = kp[2] ;
	      r3 = kp[3] ;
	      w5 = kp[4] ; if (w5 < 0) {  wiggleProblem = '#' ; w5 = -w5 ; }
	      w3 = kp[5] ; 
	      if (w3 < 0)  /* wiggle too low */
		{ 
		  wiggleProblem = '@' ; 
		  z3 = z5 = w5  ; /* true neighbour coverage */
		  w3 = w5 ;
		} 
	      else
		{
		  z3 = m3 + r3 ; if (w3 > 2 * z3) z3 = w3 ;
		  z5 = m5 + r5 ; if (w5 > 2 * z5) z5 = w5 ;
		}
	      f5 = (z5 >= 20 ?  100.0 * m5 /z5 : -20) ;
	      f3 = (z3 >= 20 ?  100.0 * m3 /z3 : -20) ;
	      f  = (z3 + z5 >= 40 ?  100.0 * (m3 + m5) /(z3 + z5) : -20) ;
	      if (ao2)
		{ /* tsf export of same data */
		  aceOutf (ao2, "%s\t%s" , dictName (tct->snpDict, up->snp), dictName(tct->runDict, ir)) ;
		  aceOutf (ao2, "\t10\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n", m5,m3,r5,r3,w5,w3, z3,z5,f3, f5) ;
		}
	    }
	  if (ao)
	    {
	      aceOut (ao, "\t") ;
	      if (0 && (f3  > f5+10 || f3 < f5 - 10))
		{ aceOutPercent (ao, f5) ; aceOut (ao, " ") ; aceOutPercent (ao, f3) ; }
	      else
		aceOutPercent (ao, f) ;
	      if (intertwinFrequencyCounts)
		{ 
		  { /* ATTENTION duplicated in ReportCounts */
		    if (tct->doubleReport)
		      {
			if (up->type & 0x1)
			  aceOutf (ao, "\t%cv:%d", wiggleProblem, m5) ;
			else
			  aceOutf (ao, "\t%cv:%d_%d", wiggleProblem, m5,m3) ;
			if (up->type & 0x2)
			  aceOutf (ao, " r:%d", r5) ;
			else
			  aceOutf (ao, " r:%d_%d", r5,r3) ;
			if (hasW)
			  {
			    if ((up->type & 0x3)  == 0x3)
			      aceOutf (ao, " w:%d", (w3+w5)/2) ;
			    else
			      aceOutf (ao, " w:%d_%d", w5, w3) ;
			  }
		      }
		    else
		      {
			aceOutf (ao, "\t%cv:%d r:%d", wiggleProblem, (m5 + m3), (r5 + r3)/2) ;
			if (hasW)
			  aceOutf (ao, " w:%d", (w5 + w3)/2) ;
		      }
		    
		  }
		}
	    }
	}
    }
} /* tctReportFrequency */

/*************************************************************************************/

static void tctReportCounts (ACEOUT ao, TCT *tct, int line, SNP *up, AC_KEYSET ks, BOOL best)
{
  DICT *runDict = tct->runDict ;
  int ir, irMax = dictMax (runDict) ;
  AC_TABLE tbl = tct->runMetaData ;
  BOOL hasW = tct->wiggleDir ? TRUE : FALSE ;

  if (1 || ! best) aceOutf (ao, "\t") ;
  if (line) /* titles */
    {
      if (best)
	{
	  aceOut (ao, "\t% variant reads 100*v/(v+ref)");
	  aceOut (ao, "\tLowest frequency\tHighest frequency") ;
	  aceOutf (ao, "\tTop run\t%s", "Top run counts v:variant, r:reference, w:wiggle at 5' and 3' sites if different") ;
	}
      else
	for (ir = 1 ; ir <= irMax ;ir++)
	  aceOutf (ao, "\t%s counting v:variant, r:reference, w:wiggle at 5' and 3' sites if different", ac_table_printable (tbl, ir - 1, line-1, "")) ;
    }
  else
    {
      int topRun = up->topRun ;
      if (best)
	{
	  int w, m ;
	  float f = -20 ;
	  
	  if (0)
	    {
	      if (up->awm > 5 * up->wm)  up->wm = up->awm ;
	      if (up->wm > 5 * up->awm)  up->awm = up->wm ;
	      if (up->awp > 5 * up->wp)  up->wp = up->awp ;
	      if (up->wp > 5 * up->awp)  up->awp = up->wp ;
	    }
	  w = up->wp + up->wm ; m = up->mp + up->mm ; f = (w + m >= 20 ? 100.0 * m / ((float)(m + w)) : -20) ;
	  aceOut (ao, "\t") ;      aceOutPercent (ao, f) ;

	  f = up->minFrequency < 200 ? up->minFrequency : -tct->minSnpCover ; aceOut (ao, "\t") ;      aceOutPercent (ao, f) ;
	  f = up->maxFrequency > -tct->minSnpCover  ? up->maxFrequency : -tct->minSnpCover ; aceOut (ao, "\t") ;      aceOutPercent (ao, f) ;
	}

      if (best && !topRun)
	aceOutf (ao, "\t\t") ;
      else
	for (ir = 1 ; ir <= irMax ;ir++)
	  { 
	    int m5 = 0, m3 = 0, r5 = 0, r3 = 0, w5 = 0, w3 = 0, z3, z5 ;
	    char wiggleProblem = ' ' ;

	    if (best)
	      {
		if (ir != topRun)
		  continue ;
		aceOutf (ao, "\t%s", dictName(tct->runDict, topRun)) ;			
	      }
	    if (up->counts)
	      {
		KEY *kp = arrayp (up->counts, MAXTYPE*ir, KEY) ;
		m5 = kp[0] ;
		m3 = kp[1] ;
		r5 = kp[2] ;
		r3 = kp[3] ;
		w5 = kp[4] ; if (w5 < 0) {  wiggleProblem = '#' ; w5 = -w5 ; }
		w3 = kp[5] ; 
		if (w3 < 0)  /* wiggle too low */
		  { 
		    wiggleProblem = '@' ; 
		    z3 = z5 = w5  ; /* true neighbour coverage */
		    w3 = w5 ;
		  } 
		else
		  {
		    z3 = m3 + r3 ; if (w3 > 2 * z3) z3 = w3 ;
		    z5 = m5 + r5 ; if (w5 > 2 * z5) z5 = w5 ;
		  }
	      }
	    
	    /* ATTENTION duplicated in ReportFrequency */
	    if (tct->doubleReport)
		{
		  if (up->type & 0x1)
		  aceOutf (ao, "\t%cv:%d", wiggleProblem, m5) ;
		  else
		    aceOutf (ao, "\t%cv:%d_%d", wiggleProblem, m5,m3) ;
		  if (up->type & 0x2)
		    aceOutf (ao, " r:%d", r5) ;
		  else
		    aceOutf (ao, " r:%d_%d", r5,r3) ;
		  if (hasW)
		  {
		    if ((up->type & 0x3) == 0x3)
		      aceOutf (ao, " w:%d", (w3+w5)/2) ;
		    else
		      aceOutf (ao, " w:%d_%d", w5, w3) ;
		  }
		}
	      else
		{
		  aceOutf (ao, "\t%cv:%d r:%d", wiggleProblem, (m5 + m3), (r5 + r3)/2) ;
		  if (hasW)
		    aceOutf (ao, " w:%d", (w5 + w3)/2) ;
		}
	  }
    }
} /* tctReportCounts */

/*************************************************************************************/

static void tctGetGlobalCounts (TCT *tct)
{
  DICT *runDict = tct->runDict ;
  int ir, irMax = dictMax (runDict) ;
  int snp, snpMax = arrayMax (tct->snps) ;
  SNP *up ;
  int minSnpCount = tct->minSnpCount ;
  int minSnpCover = tct->minSnpCover ;

  if (snpMax > 1)
    for (snp = 1, up = arrp (tct->snps, snp, SNP) ; snp < snpMax ; snp++, up++)
      {
	float alleleFrequency = 0 ;
	int nAlleleFrequency = 0 ;

	up->wiggleProblemHigh = 0 ;
	up->wiggleProblemLow = 0 ;
	for (ir = 1 ; ir <= irMax ;ir++)
	  {
	    float f = -minSnpCover ;
	    if (up->counts)
	      {
		float m5, m3, r5, r3, w5, w3, f3, f5, z3, z5 ;
		KEY *kp = arrayp (up->counts, MAXTYPE*ir, KEY) ;
		
		m5 = kp[0] ;
		m3 = kp[1] ;
		r5 = kp[2] ;
		r3 = kp[3] ;
		w5 = kp[4] ; if (w5 < 0) { w5 = -w5 ; }
		w3 = kp[5] ; 
		if (w3 < 0)  /* wiggle too low */
		  { 
		    z3 = z5 = w5  ; /* true neighbour coverage */
		    w3 = w5 ;
		  } 
		else
		  {
		    z3 = m3 + r3 ; if (w3 > 2 * z3) z3 = w3 ;
		    z5 = m5 + r5 ; if (w5 > 2 * z5) z5 = w5 ;
		  }
		if (0 && w3 + w5 == 0)
		  {
		    if (r3 > 5 * r5) r5 = kp[2] = r3 ;
		    if (r5 > 5 * r3) r3 = kp[3] = r5 ;
		  }

		if (m3 + m5 + r3 + r5 >= 2 * minSnpCover)
		  {
		    up->wellCovered++ ; 
		    f = f3 = f5 = -1000 ;
		    z3 = m3 + r3 ; if (w3 > 2 * z3) z3 = w3 ;
		    z5 = m5 + r5 ; if (w5 > 2 * z5) z5 = w5 ;
		    if ((m5 >= minSnpCount || m5 == 0) && m5 + r5 >= 20) f5 = 100.0 * m5 /z5 ;
		    if ((m3 >= minSnpCount || m3 == 0) && m3 + r3 >= 20) f3 = 100.0 * m3 /z3 ;
		    
		    /* take the max */
		    f = (m3 + m5 == 0 || m3 + m5 >= minSnpCount ? 100.0 * (m5 + m3) /(z3 + z5) : f) ;
		    if (f < f5) f = f5 ;
		    if (f < f3) f = f3 ;
		    if (f >= 95)        up->pure++ ; 
		    else if (f >= 80)   up->high++ ; 
		    else if (f >= 20)   up->mid++  ; 
		    else if (f >= 5)    up->low++  ; 
		    else                up->ref++  ; 
		    
		    if (f >= 0) { alleleFrequency += 2 * f ; nAlleleFrequency++ ; }
		    
		    if (up->maxFrequency < f)
		      { up->maxFrequency = f ; up->topRun = ir ; }
		    if (f >= 0 && up->minFrequency > f)
		      up->minFrequency = f ; 
		    if (minSnpCount  && up->maxFrequency < -minSnpCount && m3 + m5 > 0 && m3 + m5 < minSnpCount) 
		      { up->maxFrequency = -minSnpCount ;  up->topRun = ir ; }
		  }
	      }
	  }
	if (nAlleleFrequency)
	  up->alleleFrequency = alleleFrequency / nAlleleFrequency ;
      }
} /* tctGetGlobalCounts */

/*************************************************************************************/

static void tctReportGlobalCounts (ACEOUT ao, TCT *tct, int line, SNP *up, AC_KEYSET ks)
{
  aceOutf (ao, "\t") ;
  if (line) /* titles */
    {
      aceOut (ao, "\tReads supporting the variant on strand plus\tReads supporting the variant on strand minus") ;
      aceOut (ao, "\tReads supporting the reference  on strand plus\tReads supporting the reference on strand minus") ;


      aceOut (ao, "\tCoverage on plus strand of template\tCoverage on minus strand of template") ;
    }
  else
    {
      if (up->type & 0x1) aceOutf (ao, "\t%d", up->mp) ; else aceOutf (ao, "\t%d_%d", up->mp,up->amp) ;
      if (up->type & 0x1) aceOutf (ao, "\t%d", up->mm,up->amm) ;  else aceOutf (ao, "\t%d_%d", up->mm,up->amm) ;
      if (up->type & 0x2) aceOutf (ao, "\t%d", up->wp) ; else aceOutf (ao, "\t%d_%d", up->wp,up->awp) ;
      if (up->type & 0x2) aceOutf (ao, "\t%d", up->wm) ; else aceOutf (ao, "\t%d_%d", up->wm,up->awm) ;

      if ((up->type & 0x3) == 0x3) aceOutf (ao, "\t%d", up->coverp) ; else      aceOutf (ao, "\t%d_%d", up->coverp, up->acoverp) ;
      if ((up->type & 0x3)== 0x3) aceOutf (ao, "\t%d", up->coverm) ; else      aceOutf (ao, "\t%d_%d", up->coverm, up->acoverm) ;
    }
} /* tctReportGlobalCounts */

/*************************************************************************************/

static void tctReportPrevalence (ACEOUT ao, TCT *tct, int line, SNP *up, AC_KEYSET ks)
{
  aceOutf (ao, "\t") ;
  if (line == 1) /* titles */
    {
      aceOutf (ao, "\tRuns not measurable (not covered at least 20 times)\tMeasurable (covered at least 20 times)\tReference (<5%%)\tLow (5%%-20%%)\tMid (20%%-80%%)\tHigh (80%%-95%%)\tPure (over 95%%)") ;
      aceOut (ao, "\tAllele frequency in cohort (assuming diploidy)");
    }
  else if (line)
    aceOut (ao, "\t\t\t\t\t\t\t\t") ;
  else
    {
      int runMax = dictMax (tct->runDict) ;
      
      aceOutf (ao, "\t%d", runMax - up->wellCovered) ;
      aceOutf (ao, "\t%d", up->wellCovered) ;
      aceOutf (ao, "\t%d", up->ref) ;
      aceOutf (ao, "\t%d", up->low) ;
      aceOutf (ao, "\t%d", up->mid) ;
      aceOutf (ao, "\t%d", up->high) ;
      aceOutf (ao, "\t%d", up->pure) ;
      aceOutf (ao, "\t%.2f", up->alleleFrequency) ;
    }
} /* tctReportPrevalence */

/*************************************************************************************/

static void tctReportWiggleProblem (ACEOUT ao, TCT *tct, int line, SNP *up, AC_KEYSET ks)
{
  BOOL hasW = tct->wiggleDir ? TRUE : FALSE ;

  if (line == 1) /* titles */
    {
      if (hasW)
	aceOutf (ao, "\tWiggle too high\tWiggle too low") ;
      else
	aceOutf (ao, "\t\t") ;
    }
  else if (line)
    aceOut (ao, "\t\t") ;
  else
    {
      if (hasW) 
	aceOutf (ao, "\t%d\t%d", up->wiggleProblemHigh, up->wiggleProblemLow) ;
      else
	aceOut (ao, "\t\t") ;
    }
} /* tctReportWiggleProblem */

/*************************************************************************************/

static void tctReportCoding (ACEOUT ao, TCT *tct, int line, SNP *up, AC_KEYSET ks)
{
  aceOutf (ao, "\t") ;
  if (line == 1) /* titles */
    {
      aceOutf (ao, "\tHGVS_name but shifted 5' to be compatible with VCF") ;
    }
  else if (line)
    aceOut (ao, "\t") ;
  else 
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *errors = 0 ;
      AC_TABLE tbl = ac_bql_table (tct->db, "select s, map, gNam, rNam, pNam, pNam2, tit from s in @, map in s->IntMap, gNam in s->gName, rNam in s->rName, pNam in s->pName, pNam2 in s->Observed__protein_sequence[3], tit in s->title", ks, 0, &errors, h) ;      
      if (tbl)
	{
	  const char *map = ac_table_printable (tbl, 0,1, "") ;
	  const char *gNam = ac_table_printable (tbl, 0,2, "") ;
	  const char *rNam = ac_table_printable (tbl, 0,3, 0) ;
	  const char *pNam = ac_table_printable (tbl, 0,4, "") ;
	  const char *pNam2 = ac_table_printable (tbl, 0,5, "") ;
	  const char *title = ac_table_printable (tbl, 0,6, 0) ;
	  
	  if (title)
	    aceOutf (ao, "\t%s", title) ;
	  else
	    {
	      aceOutf (ao, "\t%s:", map) ;
	      aceOutf (ao, "%s", gNam) ;
	      if (rNam) aceOutf (ao, ", %s", rNam) ;
	      if (pNam[0] || pNam2[0]) aceOutf (ao, " (%s)[%s]", pNam,pNam2) ;
	    }
	}
      else
	aceOut (ao, "\t") ;

      ac_free (h) ;
    }
} /* tctReportCoding */

/*************************************************************************************/

static void tctReportSnippet (ACEOUT ao, TCT *tct, int line, SNP *up, AC_KEYSET ks)
{
  aceOutf (ao, "\t") ;
  if (line == 1) /* titles */
    {
      aceOutf (ao, "\tGenomic 51mer reference > variant") ;
      aceOutf (ao, "\tProtein 17aa reference > variant") ;
    }
  else if (line)
    aceOut (ao, "\t\t") ;
  else 
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *errors = 0 ;
       AC_TABLE tbl = ac_bql_table (tct->db, "select s, pref, pvar, gref, gvar from s in @, pref in s->Reference_protein_sequence[2], pvar in s->Observed__protein_sequence[2], gref in s->Reference_genomic_sequence, gvar in s->Observed__genomic_sequence", ks, 0, &errors, h) ;      
      if (tbl)
	{
	  const char *pR, *pV, *gR, *gV ;
	  
	  pR = ac_table_printable (tbl, 0,1,"") ;
	  aceOutf (ao, "\t%s", pR) ;
	  pV = ac_table_printable (tbl, 0,2,0) ;
	  if (pR && pV && !strcmp (pR,pV))
	    aceOutf (ao, " =") ;
	  else if (pV)
	    aceOutf (ao, " > %s", pV) ;
	  gR = ac_table_printable (tbl, 0,3,"") ;
	  aceOutf (ao, "\t%s", gR) ;
	  gV = ac_table_printable (tbl, 0,4,0) ;
	  if (gV) aceOutf (ao, " > %s", gV) ;
	}
      else
	aceOut (ao, "\t\t") ;

      ac_free (h) ;
    }
} /* tctReportSnippet */

/*************************************************************************************/

static void tctReportMethod (ACEOUT ao, TCT *tct, int line, SNP *up, AC_KEYSET ks)
{
  if (line == 1) /* titles */
    {
      aceOutf (ao, "\tMethod") ;
    }
  else if (line)
    aceOut (ao, "\t") ;
  else 
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *errors = 0 ;
      AC_TABLE tbl = ac_bql_table (tct->db, "select s, m from s in @, m in s=>method", ks, 0, &errors, h) ;      
      if (tbl)
	{
	  aceOutf (ao, "\t%s", ac_table_printable (tbl, 0,1,"")) ;
	}
      else
	aceOut (ao, "\t") ;

      ac_free (h) ;
    }
} /* tctReportCoding */

/*************************************************************************************/

static void tctReportVCF (ACEOUT ao, TCT *tct, int line, SNP *up, AC_KEYSET ks)
{
  aceOutf (ao, "\t") ;
  if (line == 1) /* titles */
    {
      aceOutf (ao, "\tto\tTarget\tVCF coordinate: substituted base or base before indel\tReference\tVariant\tModified bases") ;
    }
  else if (line)
    aceOut (ao, "\t\t\t\t\t\t") ;
  else if (up)
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *errors = 0 ;

      AC_TABLE tbl = ac_bql_table (tct->db, "select s, map, a1, a2, w1, w2, typ, ms, md, mi, from s in @, map in s->IntMap, v in s#VCF, a1 in map[1], a2 in map[2], w1 in v[2], w2 in v[3], typ in s->typ, ms in s->Multi_substitution, md in s->Multi_deletion, mi in s->Multi_insertion", ks, 0, &errors, h) ;      
      if (tbl)
	{
	  const char *ccp, *typ ;
	  int a1 = ac_table_int (tbl, 0, 2, 0) ;
	  int ms = ac_table_int (tbl, 0, 7, 0) ;
	  int a2 = ac_table_int (tbl, 0, 3, 0) ;
	  /*
	  int md = ac_table_int (tbl, 0, 8, 0) ;
	  int mi = ac_table_int (tbl, 0, 9, 0) ;
	  */
	  aceOutf (ao, "\t%d\t%s", a2, ac_table_printable (tbl, 0,1,"")) ;
	  aceOutf (ao, "\t%d", a1);

	  ccp = ac_table_printable (tbl, 0,4,"") ;  aceOutf (ao, "\t%s", ccp) ;
	  ccp = ac_table_printable (tbl, 0,5,"") ;  aceOutf (ao, "\t%s", ccp) ;
	  typ =  ac_table_printable (tbl, 0, 6, "")   ; /* type */
	  if (typ)
	    {
	      const char *ccp = strchr (typ, '>') ;
	      int nn = ccp ? ccp - typ : 0 ;
	      if (! nn && ms && typ[ms] == '2') nn = ms ;	      
	      if (! nn)
		aceOutf (ao, "\t%s", typ) ;
	      else
		{
		  int i ;
		  aceOutf (ao, "\t%c%d%c", ace_upper(typ[0]), a1, ace_upper(typ[nn+1])) ;
		  for (i = 1 ; i < nn ; i++)
		    if (ace_upper(typ[i]) !=  ace_upper(typ[nn+1+i]))
		      aceOutf (ao, ",%c%d%c", ace_upper(typ[i]), a1+i, ace_upper(typ[nn+1+i])) ;
		}
	    }
	}
      else
	aceOut (ao, "\t\t\t\t\t\t") ;

      ac_free (h) ;
    }
} /* tctReportVCF */

/*************************************************************************************/

static void tctReportLine (ACEOUT ao, ACEOUT ao2, TCT *tct, int line, SNP *up)
{
  AC_HANDLE h = 0 ;
  AC_KEYSET ks = 0 ;
  BOOL intertwinFrequencyCounts = TRUE ;

  if (up)
    {
      h = ac_new_handle () ;
      {
	const char *ccp = dictName (tct->snpDict, up->snp) ;
	if (tct->db)
	  ks = ac_dbquery_keyset (tct->db, hprintf (h, "find Variant IS \"%s\"", ccp), h) ;
      }
    }
 
  if (ao)
    {
      if (1) tctReportMethod (ao, tct, line, up, ks) ;
      if (1) tctReportWiggleProblem (ao, tct, line, up, ks) ;
      if (1) tctReportVCF (ao, tct, line, up, ks) ;
      if (1) tctReportCoding (ao, tct, line, up, ks) ;
      if (1) tctReportSnippet (ao, tct, line, up, ks) ;
      if (1) tctReportPrevalence (ao, tct, line, up, ks) ;
      if (1) tctReportCounts (ao, tct, line, up, ks, TRUE) ;
      if (1) tctReportGlobalCounts (ao, tct, line, up, ks) ;
    }
  if (1) tctReportFrequency (ao, ao2, tct, line, up, ks, intertwinFrequencyCounts);
  if (ao)
    {
      if (! intertwinFrequencyCounts) tctReportCounts (ao, tct, line, up, ks, FALSE) ;
      aceOut (ao, "\tZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ\n") ;
    }
  ac_free (h) ;
} /* tctReportLine */

/*************************************************************************************/
/* the rejected sites count as rejected in the run only if they are measurable
 * and never count as non meaurable
 */
static void tctCountRejectedSnps (TCT *tct, SNP *up)
{
  KEYSET counts = up->counts ;
  int ir, irMax = dictMax (tct->runDict) ;
  int minSnpCover = tct->minSnpCover ;
  if (! counts)
    return ;

  for (ir = 1 ; ir <= irMax ;ir++)
    {
      RC *rc = arrayp (tct->runs, ir, RC) ;
      float m5, m3, r5, r3 ;
      KEY *kp = arrayp (up->counts, MAXTYPE*ir, KEY) ;
      
      m5 = kp[0] ;
      m3 = kp[1] ;
      r5 = kp[2] ;
      r3 = kp[3] ;
     
      if (m3 + m5 + r3 + r5 >= 2 * minSnpCover)
	{
	  int k = up->varType ; 
	  if (k < 0 || k > dictMax(tct->varTypeDict))
	    k = 0 ;
	  rc->rejected++ ;
	  if (! rc->typesR)
	    rc->typesR = halloc (dictMax(tct->varTypeDict)+4 * sizeof(int), tct->h) ;
	  rc->typesR[k]++ ;
	}
    }
} /* tctCountRejected */

/*************************************************************************************/
/* the rejected sites are not rejected !
 * ditinguish well covered, ... up to pure
 */
static void tctCountAcceptedSnps (TCT *tct, SNP *up)
{
  KEYSET counts = up->counts ;
  int ir, irMax = dictMax (tct->runDict) ;
  int minSnpCover = tct->minSnpCover ;
  int minSnpCount = tct->minSnpCount ;
  if (! counts)
    return ;
  for (ir = 1 ; ir <= irMax ;ir++)
    {
      RC *rc = arrayp (tct->runs, ir, RC) ;
      float m5, m3, r5, r3, f, f3, f5, w3, w5, z3, z5 ;
      KEY *kp = arrayp (up->counts, MAXTYPE*ir, KEY) ;
      int k ;

      if (! rc->types)
	rc->types = halloc (dictMax(tct->varTypeDict)+4 * sizeof(int), tct->h) ;
      if (! rc->typesC)
	rc->typesC = halloc (dictMax(tct->varTypeDict)+4 * sizeof(int), tct->h) ;
      if (! rc->typesR)
	rc->typesR = halloc (dictMax(tct->varTypeDict)+4 * sizeof(int), tct->h) ;

      m5 = kp[0] ;
      m3 = kp[1] ;
      r5 = kp[2] ;
      r3 = kp[3] ;
      w5 = kp[4] ;
      w3 = kp[5] ;
     
      rc->tested++ ;
      if (m3 + m5 + r3 + r3 < 2 * minSnpCover)
	{ rc->notMeasurable++ ; continue ; }
      rc->wellCovered++ ;

      f = f3 = f5 = -1000 ;
      z3 = m3 + r3 ; if (w3 > 2 * z3) z3 = w3 ;
      z5 = m5 + r5 ; if (w5 > 2 * z5) z5 = w5 ;
      if ((m5 >= minSnpCount || m5 == 0) && m5 + r5 >= 20) f5 = 100.0 * m5 /z5 ;
      if ((m3 >= minSnpCount || m3 == 0) && m3 + r3 >= 20) f3 = 100.0 * m3 /z3 ;
      
      /* take the max */
        f = (m3 + m5 == 0 || m3 + m5 >= 2 * minSnpCount ? 100.0 * (m5 + m3) /(z3 + z5) : f) ;
      if (f < f5) f = f5 ;
      if (f < f3) f = f3 ;
      if (f >= 95)        rc->pure++ ; 
      else if (f >= 80)   rc->high++ ;
      else if (f >= 20)   rc->mid++ ; 
      else if (f >= 5)    rc->low++ ;
      else                rc->ref++ ;
		    
      k = up->varType ;
      if (k < 0 || k > dictMax(tct->varTypeDict))
	k = 0 ;
      rc->types[k] ++ ;
      rc->typesC[k] += (m3 + m5)/2 ; 
    }
} /* tctCountAccepted */

/*************************************************************************************/

static int tctReport (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, nn = 0, snpMax  = arrayMax (tct->snps) ;
  DICT *snpDict = tct->snpDict ;

  ACEOUT ao = aceOutCreate (tct->outFileName, ".snp_frequency_table.txt", tct->gzo, h) ;
  ACEOUT ao2 = tct->outFileName ? aceOutCreate (tct->outFileName, ".snp_counts.tsf", tct->gzo, h) : 0 ;
  aceOutDate (ao, "####", "SNPs and variations table") ;
  if (ao2)
    {
      aceOutDate (ao2, "####", "SNPs and variations count variants, reference and wiggle cover at donor and acceptor sites :\n") ;
      aceOutf (ao2, "# SNP\tRun\t10\tVar5\tVar3\tRef5\tRef3\tCover5\tCover3\tWiggle5\tWiggel3\tFreq5\tFreq3\n") ;
    }
  tct->doubleReport = TRUE ;

  aceOutf (ao, "# Title") ; tctReportLine (ao, ao2, tct, 2, 0) ;
  aceOutf (ao, "# Sorting_title") ; tctReportLine (ao, ao2, tct, 3, 0) ;
  aceOutf (ao, "# Run") ; tctReportLine (ao, ao2, tct, 1, 0) ;
    
  for (ii = 1 ; ii < snpMax  ; ii++)
    {
      SNP *up = arrayp (tct->snps, ii, SNP) ;

      if (tct->t1 && (up->a1 > tct->t2 || up->a2 < tct->t1 || up->a1 < tct->t1 - 1000 || up->a2 > tct->t2 + 1000)) continue ;

      if (! up->select && tct->minSnpFrequency > 0 && (up->maxFrequency < 0 || up->maxFrequency < tct->minSnpFrequency))
	continue ;
      if (up->coverp + up->coverm + up->acoverp + up->acoverm < 2 * tct->minSnpCover)
	continue ;
      if (up->mp + up->mm + up->amp + up->amm < 2 * tct->minSnpCount)
	continue ;
      if (! up->select && tct->dropMonomodal && up->pure + up->high + up->mid == 0)
	{
	  tctCountRejectedSnps (tct, up) ;
	  continue ;
	}
      tctCountAcceptedSnps (tct, up) ;
      aceOutf (ao, "%s", dictName (snpDict, up->snp)) ;
      tctReportLine (ao, ao2, tct, 0, up) ;
      nn++ ;
    }
  ac_free (h) ;
  return nn ;
} /* tctReport */

/*************************************************************************************/

static int tctMergeExport (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, nn = 0, snpMax  = arrayMax (tct->snps) ;

  ACEOUT ao = aceOutCreate (tct->outFileName, ".snp_frequency.tsf", tct->gzo, h) ;
  aceOutDate (ao, "####", "SNPs counts and frequency table") ;
  aceOutf (ao, "# SNP\tRun\t10\tVar5\tVar3\tRef5\tRef3\tCover3\tCover5\tWiggle5\tWiggel3\tFreq5\tFreq3\n") ;

  for (ii = 1 ; ii < snpMax  ; ii++)
    {
      SNP *up = arrayp (tct->snps, ii, SNP) ;

      if (tct->t1 && (up->a1 > tct->t2 || up->a2 < tct->t1 || up->a1 < tct->t1 - 1000 || up->a2 > tct->t2 + 1000)) continue ;

      if (! up->select && tct->minSnpFrequency > 0 && (up->maxFrequency < 0 || up->maxFrequency < tct->minSnpFrequency))
	continue ;
      if (up->coverp + up->coverm + up->acoverp + up->acoverm < 2 * tct->minSnpCover)
	continue ;
      if (up->mp + up->mm + up->amp + up->amm < 2 * tct->minSnpCount)
	continue ;
      if (! up->select && tct->dropMonomodal && up->pure + up->high + up->mid == 0)
	{
	  continue ;
	}
      tctReportLine (0, ao, tct, 0, up) ;
      nn++ ;
    }
  ac_free (h) ;
  return nn ;
} /* tctMergeExport */

/*************************************************************************************/
/*************************************************************************************/
/* export snp stats per run that will go to the Ali object then to qc_ssummary */
static void tctReportRuns (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tct->outFileName, "snp_runs.ace", tct->gzo, h) ;
  int i, ir, irMax = arrayMax (tct->runs) ;
  DICT *varTypeDict = tct->varTypeDict ;
  int iMax = dictMax (varTypeDict) ;
  RC *rc ;

  if (irMax)
    for (ir = 1, rc = arrp (tct->runs, ir, RC) ; ir < irMax ; rc++, ir++)
      {
	aceOutf (ao, "Ali \"%s\"\n", dictName (tct->runDict, ir)) ;
	aceOut (ao, "-D SNP\n") ;
	
	if (rc->tested)
	  {
	    aceOutf (ao, "Tested_sites %d any\n", rc->tested) ;
	    aceOutf (ao, "Rejected_sites %d any\n", rc->rejected) ;
	    aceOutf (ao, "Not_measurable_sites %d any\n", rc->notMeasurable) ;
	  }
	if (rc->wellCovered)
	  {
	    aceOutf (ao, "Measured_sites Genomic %d any %d reference %d low %d mid %d high %d pure\n"
		     , rc->wellCovered
		     , rc->ref
		     , rc->low
		     , rc->mid
		     , rc->high
		     , rc->pure
		     ) ;	     
	    aceOutf (ao, "Support_by_sites Genomic %%.0f any %%.0f reference %%.0f low %%.0f mid %%.0f high %%.0f pure\n"
		     , rc->wellCoveredC
		     , rc->refC
		     , rc->lowC
		     , rc->midC
		     , rc->highC
		     , rc->pureC
		     ) ;	     
	  }
	
	if (rc->xwellCovered)
	  aceOutf (ao, "Measured_sites Exonic %d any %d reference %d low %d mid %d high %d pure\n"
		   , rc->xwellCovered
		   , rc->xref
		   , rc->xlow
		   , rc->xmid
		   , rc->xhigh
		   , rc->xpure
		   ) ;
	if (rc->pwellCovered)
	  {
	    aceOutf (ao, "Measured_sites Protein_changing %d any %d reference %d low %d mid %d high %d pure\n"
		     , rc->pwellCovered
		     , rc->pref
		     , rc->plow
		     , rc->pmid
		     , rc->phigh
		     , rc->ppure
		     ) ;	    
	    aceOutf (ao, "Support_by_sites Protein_changing %%.0f any %%.0f reference %%.0f low %%.0f mid %d high %%.0f pure\n"
		     , rc->pwellCoveredC
		     , rc->prefC
		     , rc->plowC
		     , rc->pmidC
		     , rc->phighC
		     , rc->ppureC
		     ) ;	    
	  }
	for (i = 1 ; i <= iMax ; i++)
	  {
	    int n = rc->types[i] ;
	    int nR = rc->typesR[i] ;
	    float nS = rc->typesC[i] ;
	    if (n)
	      aceOutf (ao, "SNP_profile any %s %d %d %.0f\n", dictName (varTypeDict, i), n, nR, nS) ;
	  }
	for (i = 0 ; i < 1 ; i++)
	  {
	    int n = rc->types[i] ;
	    int nR = rc->typesR[i] ;
	    float nS = rc->typesC[i] ;
	    if (n)
	      aceOutf (ao, "SNP_profile any %s %d %d %.0f\n", "Other", n, nR, nS) ;
	  }

	aceOut (ao, "\n") ;
      }
  ac_free (h) ;
  return ;
} /* tctReportRuns */

/*************************************************************************************/
/*************************************************************************************/
static DICT *tctMakeVarTypeDict (AC_HANDLE h) ;
static void tctDbReport (TCT *tct)
{
  tct->varTypeDict = tctMakeVarTypeDict (tct->h) ;

  tctSnpParse (tct) ;
  if (tct->snps && arrayMax (tct->snps))
    {
      tctGetGlobalCounts (tct) ;
      tctWiggleParse (tct) ;
      tct->snpExported = tctReport (tct) ;
      tctReportRuns (tct) ;
    }
} /* tctDbReport */
  
/*************************************************************************************/
/*************************************************************************************/
static DICT *tctMakeVarTypeDict (AC_HANDLE h) ;
static void tctMergeCounts (TCT *tct)
{
  tct->runDict = dictHandleCreate (256, tct->h) ;
  tct->varTypeDict = tctMakeVarTypeDict (tct->h) ;
  tct->minSnpCount = tct->minSnpFrequency = 0 ;
  tctSnpParse (tct) ;
  if (tct->snps && arrayMax (tct->snps))
    {
      tctGetGlobalCounts (tct) ;
      tctMergeExport (tct) ;
    }
} /* tctMergeCounts */
  
/*************************************************************************************/
/*************************************************************************************/
/* slide deletions and insertions, speclial treat for position 76 */

static BOOL tctSlide (const char *dna, int dnaLn, int a10, int da, int *dxp, int *slidep)
{
  int i, dx, dy, a1 = a10, xx ; 
  BOOL isDim = FALSE ;
  
  /* slide left */
  dx = 0 ; i = a1 - 1 ;
  while (i > 0 && dna[i] == dna[i + da])
    { i-- ; dx++ ; }
  a1 -= dx ;
  /* slide right */
  dy = 0 ; i = a1 ;
  while (i + da < dnaLn - 1 && dna[i] == dna[i + da])
    { i++ ; dy++ ; }
  if (dy >= da)
    isDim = TRUE ;
  xx = 75 ;
  if (a1 < xx && a1 + dy >= xx && da > 1000)
    { dy = xx - a1 ; a1 += dy ; }

  *slidep = dy ;
  *dxp = a1 - a10 ;
  return isDim ;
} /* tctSlide */

/*************************************************************************************/
/* slide deletions and insertions, speclial treat for position 76 */

static BOOL tctSlideDup (const char *dna, int dnaLn, int a10, int da, int *dxp, int *slidep, const char *insert, char *buf)
{
  int i, j, k, dx, dy, a1 = a10 ;
  BOOL isDup = FALSE ;
  
  /* slide left */
  dx = 0 ; i = a1 - 1 ; j = da - 1 ;
  while (i > 0 && dna[i] == ace_lower(insert [j]))
    { i-- ; dx++ ; j-- ; if (j < 0) j += da ; }
  a1 -= dx ;
  /* copy the rotated buffer */
  for (k = 0 ; k < da ; k++)
    buf[k] = insert[(j+1+k+da) % da] ;
  buf[k] = 0 ;
  /* slide right */
  dy = 0 ; i = a1 ;
  while (i < dnaLn - 1 && dna[i] == ace_lower(insert[j]))
    { i++ ; dy++ ; j = (j+1) % da ; }
  if (dy >= da)
    isDup = TRUE ;

  *slidep = dy ;
  *dxp = a1 - a10 ;
  return isDup ;
} /* tctSlideDup */

/*************************************************************************************/

static const char *Rdna = 0 ;
static Array Gdna = 0 ;
static int RdnaLn = 0 ;
static int myMrna = 0 ;
static int g2m = 0 ;
static int dnaLn = 0 ;
BOOL fromMrna = FALSE ;

static char myRdna (int a) 
{
  int m = a + (fromMrna ? 0 : g2m) ;
  int cc = 'n' ;
  
  if (m >=0 && m < RdnaLn)
    {
      cc = Rdna[m] ;
      if (cc == 't') cc = 'u' ;
    }
  return cc ;
} /* myRdna */

/*************************************************************************************/

static int myTranslate (char *buf, char tBuf1[], char tBuf3[], int m1, int p1, AC_HANDLE h) 
{
  Array aa = arrayHandleCreate (strlen (buf) + 8, char, h) ;
  char * translationTable = pepGetTranslationTable(myMrna, 0) ; 
  char cc, *cp1, *cp2, *cB1, *cB3 ;
  int frame = 0, j, i = 50 + strlen (buf) ;
  char *tbuf = halloc (i, h) ;
  
  array(aa, i, char) = 0 ;
  arrayMax(aa) = i ;
  cp1 = arrp (aa, 0, char) ; cp2 = buf ;
  for (i = 0 ; *cp2 ; cp2++)
    { cc = dnaEncodeChar [(int)*cp2] ; if (cc) { *cp1++ = cc ; i++ ; }}
  *cp1 = 0 ; 
  if (i >= arrayMax(aa)) messcrash ("n importe quoi dans myTranslate") ;
  arrayMax (aa) = i ;
  memset (tbuf, 0, i) ;
  tbuf[0] = 0 ;
  i = m1 ;
  i = i % 3 ;
  i += 9 ;
  i = i % 3 ;
  cB1 = tBuf1 ;   cB3 = tBuf3 ; 
  for (j = 0 ; j < i ; j++)
    *cB3++ = '-' ;
  if (arrayMax (aa) > 2)
    for (cp1 = arrp (aa, i, char)  ; i < arrayMax (aa) - 2 ; i+= 3 , cp1 += 3 )
      {
	int j ;
	char cc = e_codon (cp1, translationTable) ;
	const char *pep = pepShortName[(int)cc] ;
	if (i>=19 && i <=21) cc = ace_lower (cc) ;
	*cB1++ = cc ;
	if (pep[0] == 'X') pep = "---" ;
	if (cc == '*') pep = "Ter" ;
	for (j = 0 ; j < 3 ; j++)
	  *cB3++ = *pep++ ;
      }
  *cB1 = *cB3 = 0 ;
  return frame ;
} /* myTranslate */

/*************************************************************************************/
/* if i kill Met1, then use Leu2_Met124del  (i.e. we now use Met14 as our new Met)
 * stops are caller Ter
 * extension: the stop has a substitution :  Ter110GlyextTer31
 * (Arg123LysfsTer34) 
 */
static void tctSetPName (vTXT txt, KEY product, char *pR1, char *pV1, char *pRef, char *pVar, int m1, int p1, int frame, int fs) 
{
  char *cp, *cq ;
  char * translationTable = pepGetTranslationTable(myMrna, 0) ; 
  char *sep = "" ;
  char buf[6];
  char buf2[6];
  char idem[4] ;
  int i=0, j, k = 0, m ;
  m1 /= 3 ; m1++ ;
  
  if (1) /* sub del ins */
    {
      if (1) vtxtPrintf (txt, "\npName p.%s", name(product)) ;
      cp = pRef ; cq = pVar ;
      while (*cp == '-') { cp++; cq++ ; i++; }
      while (*cp)
	{
	  m = m1 - 7+(i+2)/3 ;
	  if (m == m1)
	    {  
	      for (j = 0 ; j < 3 ; j++)
		idem[j] = cp[j] ;
	      idem[3] = 0 ;
	    }
	  if (strncmp (cp, cq,3))
	    {
	      k++ ;
	      for (j = 0 ; j < 3 ; j++)
		buf[j] = cp[j] ;
	      buf[3] = 0 ;
	      for (j = 0 ; j < 3 ; j++)
		buf2[j] = cq[j] ;
	      buf2[3] = 0 ;
	      if (! strncmp (buf2, "Ter", 3))
		{ 
		  vtxtPrintf (txt, "%s%s%dTer", sep, buf,m,buf2) ;
		  break ;
		}
	      else if (fs == 0)
		vtxtPrintf (txt, "%s%s%d%s", sep, buf,m,buf2) ;
	      else if (fs == -3) /* del in frame */
		vtxtPrintf (txt, "%s%ddel", sep, buf,m) ;
	      else if (fs < 0 && fs % 3 == 0)  /* multi_deletion in frame */
		{
		  vtxtPrintf (txt, "%s%s%ddel", sep, buf,m) ;
		}
	      if (fs)
		{
		  if (fs == -3)
		    vtxtPrintf (txt, "%ddel%s", m, buf) ;
		  else if (fs == 3)
		    vtxtPrintf (txt, "%d_%dins%s", m,m+1, buf) ;
		  else if (fs <0 && (-fs % 3 == 0))
		    vtxtPrintf (txt, "%d_%ddel", m, m - fs/3 -1) ;
		  else if (fs >0 && (fs % 3 == 0))
		    vtxtPrintf (txt, "%d_%dins(%d)", m, m + 1, fs/3) ;
		  else if (fs %3)
		    vtxtPrintf (txt, "%dfs(%d)", m, fs) ;
		  
		  break ; 
		}

	      sep = ";" ;
	     
	      if (m >= 1 && ! strncmp (cp, "Ter", 3))
		{
		  char *cr = cq ;
		  int i1 = i ;
		  while (i1 < 51 && strncmp (cr, "Ter",3))
		    { cr+= 3 ; i1 += 3 ; }
		  m = m1 - 7+(i1+2)/3 ;
		  if (!strncmp (cr, "Ter",3))
		    vtxtPrintf (txt, "extTer%d",m) ;
		  else if (Gdna)
		    { /* the mRNA is too short to find the Ter, look in the genome */
		      const char *ccr ; 
		      int a = m1 - g2m + 2 + (m1 %3) ; /* start of the Ter codon on the genome */  
		      for (ccr = arrp (Gdna, a, char) ; *cr && a < dnaLn - 2 ; ccr += 3, a += 3)
			if (e_codon (ccr, translationTable) == '*')
			  break ;
		      vtxtPrintf (txt, "extTer%d", m1 + (a - m1 + g2m + 2)/3) ;
		    }
		  break ;
		}
	    }
	  cp += 3 ; cq += 3 ; i += 3 ;
	}
      if (k == 0) vtxtPrintf (txt, "%s%d=", idem,m1) ;
      if (0) vtxtPrintf (txt, ")\"") ;
    }
  /* clean up the pR1 pV1 */
  if (m1 > 0)
    {
      cp = strchr (pR1 + 6, '*') ;
      if (cp) *++cp = 0 ;
    }
  if (m1 > 0)
    {
      cp = strchr (pV1 + 6, '*') ;
      if (cp) *++cp = 0 ;
    }
  return ; 
}  /* tctSetPName */

/*************************************************************************************/
/* biologits, like Socrates, never heard of zero */
static int socrate (int x)
{
  return x > 0 ? x : x - 1 ; 
} /* socrate */

/*************************************************************************************/

static int tctSetGName (vTXT txt, TCT *tct, AC_OBJ Snp, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = ac_tag_table (Snp, "IntMap", h) ;
  AC_TABLE Rtbl = ac_tag_table (Snp, "mRNA", h) ;
  KEY seq = ac_table_key (tbl, 0, 0, 0) ;
  KEY mrna = ac_tag_key (Snp, "mRNA", 0) ;
  static KEY oldSeq = 0 ;
  static KEY oldMrna = 0 ;
  static const char *dna = 0 ;
  BOOL ok = FALSE ;
  char bufR[52], bufV[52] ;
  char RbufR[52], RbufV[52] ;
  
  fromMrna = ac_has_tag (Snp, "Found_in_mRNA") ;
  /*
    KEY product = 0 ; 
    int mrnaA1 = 0 ;   // a coord of the first base of the mRNA in the map 
  int productX1 = 0 ; // x coord of the first base of the coding region in the mRNA 
  */
  memset (bufR, 0, sizeof (bufR)) ;
  memset (bufV, 0, sizeof (bufV)) ;
  memset (RbufR, 0, sizeof (bufR)) ;
  memset (RbufV, 0, sizeof (bufV)) ;


  if (! seq)
    {
      if (dna != Rdna) ac_free (dna) ;
      dna = 0 ;
      dnaLn = 0 ;
      oldSeq = 0 ;      
      ac_free (Gdna) ;
    }
  else if (seq != oldSeq && ! fromMrna)
    {
      AC_OBJ Seq = ac_tag_obj (Snp, "Parent_sequence", h) ;

      ac_free (dna) ;
      oldSeq = seq ;
	
      dna = Seq ? ac_obj_dna (Seq, h0) : 0 ;
      dnaLn = dna ? strlen (dna) : 0 ;
      ac_free (Gdna) ;
      Gdna = dnaGet (seq) ;
    }

  if (! mrna)
    {
      if (RdnaLn) ac_free (Rdna) ;
      Rdna = dna ;
      RdnaLn = 0 ;
      oldMrna = 0 ;
      if (!seq) 
	return 0 ;
    }
  else if (mrna != oldMrna)
    {
      AC_OBJ Mrna = ac_tag_obj (Snp, "mRNA", h) ;

      if (RdnaLn) ac_free (Rdna) ;
      oldMrna = mrna ;
      Rdna = Mrna ? ac_obj_dna (Mrna, h0) : 0 ;
      RdnaLn = Rdna ? strlen (Rdna) : 0 ;
      ac_free (Mrna) ;
    }
  myMrna = mrna  ;
  if (fromMrna)
    {
      dna = Rdna ;
      dnaLn = RdnaLn ;
    }

  if (! dna)
    return 0 ;
  if (tbl && tbl->cols >= 3)
    {

      const char *ccp = ac_table_printable (tbl, 0, 0, "xxx") ;
      int i, j ;
      int fs = 0 ; /* frameshift */
      int dda = 0 ; /* use to get the translation frame */
      int m1 = ac_table_int (Rtbl, 0, 1, 0) ;
      int m2 = ac_table_int (Rtbl, 0, 2, 0) ;
      int da = 1 ;
      int a1 = ac_table_int (tbl, 0, 1, 0) ;
      int a2 = ac_table_int (tbl, 0, 2, 0) ;
      char **sub, *subs[] = {"A2T","A2G","A2C","T2A","T2G","T2C","G2A","G2T","G2C","C2A","C2T","C2G",0} ;
      char **del, *dels[] = {"DelA","DelT","DelG","DelC", 0} ;
      char **ins, *inss[] = {"InsA","InsT","InsG","InsC", 0} ;
      g2m = RdnaLn ? m1 - a1 : 0 ;
      
      if (a2 == 1 || a2 == -1)
	{
	  char *cp = strnew (ac_name(Snp), 0), *cq, *cr ;
	  cq = strchr (cp, ':') ;
	  if (cq)
	    {
	      cq++ ; 
	      cr = strchr (cq, '_') ;
	      if (cr)
		{
		  da = 0 ;
		  cr++ ;
		  while (*cr >= '0' && *cr <= '9')
		    { da = 10 * da + (*cr - '0') ; cr++ ; }
		}
	    }	  
	  if (a2 == 1)
	    a2 = a1 + da + 1 ;
	  else
	    a2 = a1 - (da + 1) ;
	  if (m2 == 1)
	    m2 = m1 + da + 1 ;
	  else
	    m2 = m1 - (da + 1) ;
	  ac_free (cp) ;
	}      
      if (strstr (ac_name(Snp), ":Sub"))
	{
	  if (! ac_has_tag (Snp, "Multi_substitution"))
	    {
	      for (sub = subs ; *sub ; sub++)
		{
		  if (ac_has_tag (Snp, *sub))
		    {
		      int am1 = fromMrna ? m1 : a1 ;
		      char buf[4] ;
		      buf[0] = ace_lower (sub[0][0]) ;
		      buf[1] = '>' ;
		      buf[2] = ace_lower(sub[0][2]) ;
		      buf[3] = 0 ;
		      ok = TRUE ;
		      vtxtPrint (txt, "-D Substitution\n") ; /* cleanup */
		      vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Base %d %s is modified\"\n", name (seq), a1, a2, a1+1, *sub) ;
		      vtxtPrintf (txt,"Typ %s\n%s\n", buf, *sub) ; /* reinstate */ 
		      vtxtPrintf (txt, "gName \"g.%d%c>%c\"\n", a1+1, sub[0][0], sub[0][2]) ;
		      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d%c>%c\"\n", name(mrna), socrate(m1 + 1), sub[0][0], sub[0][2]) ;
		      vtxtPrintf (txt, "VCF %d %c %c\n", a1+1, sub[0][0], sub[0][2]) ;
		      dda = 0 ; /* -1 for Socrates, but +1 because the sub is between a1 and a2, at a1+1 */

		      am1 = fromMrna ? m1 : a1 ;
		      for (i = am1 - 21 + dda, j = 0 ; i < am1 + 30 + dda && i < dnaLn ; i++)
			if (i > 0 && j < 51) 
			  {
			    if (i == am1)
			      { 
				bufR[j] = ace_upper (dna[i]) ;
				RbufR[j] = ace_upper (myRdna(i)) ;
				bufV[j] = sub[0][2] ;
				RbufV[j] = sub[0][2] ;
			      }
			    else
			      {
				bufV[j] = bufR[j] = dna[i] ;
				RbufV[j] = RbufR[j] = myRdna(i) ;
			      }
			    j++ ;
			  }
		      RbufR[j] = bufR[j] = 0 ;
		      RbufV[j] = bufV[j] = 0 ;
		    }
		}
	    }
	  else
	    {
	      tbl = ac_tag_table (Snp, "Multi_substitution", h) ;
	      if (tbl->cols >= 3)
		{
		  int k,  da = a2 - a1 ;
		  const char *ccR = ac_table_printable (tbl, 0, 1, "") ; 
		  const char *ccV = ac_table_printable (tbl, 0, 2, "") ; 
		  int am1 = fromMrna ? m1 : a1 ;
		  int am2 = fromMrna ? m2 : a2 ;
		  vtxtPrint (txt, "-D Substitution\n") ; /* cleanup */
		  vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Bases %d to %d (%d bases) are modified\"\n", name (seq), a1, a2, a1, a2-1, da) ;
		  vtxtPrintf (txt, "Multi_substitution %d %s %s\n", da, ccR, ccV) ; /* reinstate */
		  if (ccV && strlen (ccV) == da)
		    {
		      ok = TRUE ;
		      vtxtPrintf (txt, "VCF %d %s %s\n", a1+1, ccR, ccV) ;
		      vtxtPrintf (txt, "gName \"g.%d_%ddelins%s\"\n", a1,a2 -1, ccV) ;
		      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddelins%s\"\n", name(mrna),socrate (g2m+a1),socrate (g2m+a2 -1), ccV) ;
		      vtxtPrintf (txt, "Typ \"%s2%s\"\n", ccR, ccV) ;
		      dda = 0 ; /* -1 for Socrates, but +1 because the sub is between a1 and a2, starting at a1+1 */
		      
		      for (i = am1 - 21 + dda, j = k = 0 ; i < am1 + 30 + dda && i < dnaLn ; i++)
			if (i > 0 && j < 51)
			  {
			    if (i >= am1 -1 && i <= am2 - 2)
			      { 
				bufR[j] = ace_upper (dna[i]) ;
				bufV[j] = ace_upper (ccV[i - a1 + 1]) ;
				RbufR[j] = ace_upper (myRdna(i)) ;
				RbufV[j] = ace_upper (ccV[i - a1 + 1]) ;
				if (bufR[j] == bufV[j])
				  bufR[j] = bufV[j] = dna[i] ; /* lower the common bases */
				if (RbufR[j] == RbufV[j])
				  RbufR[j] = RbufV[j] = myRdna(i) ; /* lower the common bases */
			      }
			    else
			      {
				bufV[j] = bufR[j] = dna[i] ;
				RbufV[j] = RbufR[j] = myRdna(i) ;
			      }
			    j++ ;
			  }
		      RbufR[j] = bufR[j] = 0 ;
		      RbufV[j] = bufV[j] = 0 ;
		    }
		}
	    }
	}	    
      
      else if (ac_has_tag (Snp, "Deletion"))   /* happens because dimA is under DelA and should be under Diminution */
	{
	  if (! ac_has_tag (Snp, "Multi_deletion"))   /* happens because dimA is under DelA and should be under Diminution */
	    {
	      for (del = dels ; *del ; del++)
		if (ac_has_tag (Snp, *del))
		  {
		    int da = 1 ;
		    int dx = 0 ;
		    int slide = 0 ;
		    int am1 = fromMrna ? m1 : a1 ;
		    BOOL isDim = tctSlide (dna, dnaLn, am1, da, &dx, &slide) ;
		    if (a1 < a2) { a1 += dx ; a2 += dx ;}
		    else { a1 -= dx ; a2 -= dx ;}
		    if (fromMrna) { m1 += dx ; m2 += dx ; }
		    am1 = fromMrna ? m1 : a1 ;
		    fs = -1 ; 

		    vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Base %c %d is deleted\"\n", name (seq), a1, a2, (*del)[3], a1+1) ;
		    vtxtPrintf (txt, "-D Sliding\n") ;
		    vtxtPrint (txt, "-D Deletion\n") ; /* cleanup */
		    vtxtPrintf (txt,"%s\n", *del) ; /* reinstate */
		    
		    
		    ok = TRUE ;
		    if (isDim)
		      {
			vtxtPrintf (txt, "gName \"g.%ddim%c\"\n", a1+1, ace_upper(dna[am1])) ;
			if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%ddim%c\"\n", name(mrna), socrate (g2m+a1+1), ace_upper(dna[am1])) ;
			vtxtPrintf (txt, "Typ \"Dim%c\"\n", ace_upper(dna[am1])) ;
			vtxtPrintf (txt, "Diminution\n") ;
		      }
		    else
		      {
			vtxtPrintf (txt, "gName \"g.%ddel%c\"\n", a1+1, ace_upper(dna[am1])) ;
			if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%ddel%c\"\n", name(mrna), socrate (g2m+a1+1), ace_upper(dna[am1])) ;
			vtxtPrintf (txt, "Typ \"Del%c\"\n", ace_upper(dna[am1])) ;
		      }
		    if (slide) 
		      vtxtPrintf (txt, "Sliding %d\n", slide) ;
		    vtxtPrintf (txt, "VCF %d %c%c %c\n", a1,  dna[am1-1], ace_upper(dna[am1]), dna[am1-1]) ;

		    am1 = fromMrna ? m1 : a1 ;
		    for (i = am1 - 21, j = 0 ; i < am1 + 30 + da && i + da < dnaLn ; i++)
		      if (i > 0 && j < 51)
			{
			  if (i < am1)
			    {
			      bufV[j] = bufR[j] = dna[i] ;
			      RbufV[j] = RbufR[j] = myRdna(i) ;
			    }
			  else if (i == a1)
			    {
			      bufR[j] = ace_upper (dna[i]) ;
			      bufV[j] = '_' ; 
			      RbufR[j] = ace_upper (myRdna(i)) ;
			      RbufV[j] = '_' ; 
			    }
			  else if (i > am1)
			    {
			      bufR[j] =  bufV[j] = dna[i] ;
			      RbufV[j] = RbufR[j] = myRdna(i) ;
			    }
			  j++ ;
			}
		    RbufR[j] = bufR[j] = 0 ;
		    RbufV[j] = bufV[j] = 0 ;
		  }
	    }
	  else
	    {
	      tbl = ac_tag_table (Snp, "Multi_deletion", h) ;
	      if (tbl)
		{
		  int da = a2 - a1 - 1 ;
		  int dx = 0, k, kk ;
		  int slide = 0 ;
		  int am1 = fromMrna ? m1 : a1 ;
		  int am2 = fromMrna ? m2 : a2 ;
		  BOOL isDim = tctSlide (dna, dnaLn, am1, da, &dx, &slide) ;
		  char bufN[15] ;
		  if (a1 < a2) { a1 += dx ; a2 += dx ;}
		  else { a1 -= dx ; a2 -= dx ;}
		  if (fromMrna) { m1 += dx ; m2 += dx ; }
		  am1 = fromMrna ? m1 : a1 ;
		  am2 = fromMrna ? m2 : a2 ;
		  fs = -da ;
		  vtxtPrintf (txt, "-D Sliding\n") ;
		  vtxtPrint (txt, "-D Deletion\n") ; /* cleanup */
		  
		  vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Bases %d to %d (%d bases) are deleted\"\n", name (seq), a1, a2, a1+1, a2 - 1, da) ;
		  ok = TRUE ;
		  /* temporarily store the deleted bases in bufV */
		  for (i = am1, j = 0 ; j < 51 && i < am2 - 1 && i < dnaLn ; i++)
		    bufV[j++] = ace_upper(dna[i]) ;
		  bufV[j] = 0 ;
		  
		  vtxtPrintf (txt,"Multi_deletion %d %s\n", da, bufV) ; /* reinstate */ 
		  if (da <= 20)
		    {
		      vtxtPrintf (txt, "VCF %d %c%s %c\n", a1, dna[am1-1], bufV, dna[am1-1]) ;
		      if (isDim)
			{
			  vtxtPrintf (txt, "gName \"g.%d_%ddim%s\"\n", a1+1,a2-1, bufV) ;
			  if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddim%s\"\n", name(mrna), socrate (g2m+a1+1),socrate (g2m+a2-1), bufV) ;
			  vtxtPrintf (txt, "Typ \"Dim%s\"\n", bufV) ;
			  vtxtPrintf (txt, "Diminution\n") ;
			}
		      else 
			{
			  vtxtPrintf (txt, "gName \"g.%d_%ddel%s\"\n", a1+1,a2-1, bufV) ;
			  if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddel%s\"\n", name(mrna), socrate (g2m+a1+1),socrate (g2m+a2-1), bufV) ;
			  vtxtPrintf (txt, "Typ \"Del%s\"\n", bufV) ;
			    }
		      vtxtPrintf (txt, "Multi_deletion %d %s\n", da, bufV) ;
		    }
		  else
		    {
		      vtxtPrintf (txt, "Typ \"Del_%d\"\n", da) ;
		      vtxtPrintf (txt, "VCF %d %c%d %c\n", a1, dna[am1-1], da, dna[am1-1]) ;
		      vtxtPrintf (txt, "gName \"g.%d_%ddel(%d)\"\n", a1+1,a2-1, da) ;
		      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddel(%d)\"\n", name(mrna), socrate (g2m+a1+1),socrate (g2m+a2-1), da) ;
		      vtxtPrintf (txt, "Multi_deletion %d\n", da) ;
		    }
		  if (slide) 
		    vtxtPrintf (txt, "Sliding %d\n", slide) ;
		  bufN[0] = 0 ;
		  if (da ==  1) sprintf (bufN, "_") ; 
		  if (da ==  2) sprintf (bufN, "__") ; 
		  if (da ==  3) sprintf (bufN, "___") ; 
		  if (da ==  4) sprintf (bufN, "____") ; 
		  if (da ==  5) sprintf (bufN, "_____") ; 
		  if (da ==  6) sprintf (bufN, "______") ; 
		  if (da ==  7) sprintf (bufN, "_______") ; 
		  if (da ==  8) sprintf (bufN, "________") ; 
		  if (da ==  9) sprintf (bufN, "_________") ; 
		  if (da == 10) sprintf (bufN, "__________") ; 
		  if (da == 11) sprintf (bufN, "___________") ; 
		  if (da == 12) sprintf (bufN, "____________") ; 
		  if (da > 12) sprintf (bufN, "_%d_", da) ; 
		  
		  am1 = fromMrna ? m1 : a1 ;
		  am2 = fromMrna ? m2 : a2 ;
		  for (i = am1 - 21, k = kk = j = 0 ; i < am1 + 30 && i < dnaLn ; i++)
		    if (i > 0 && k < 51)
		      {
			if (i < am1)
			  {
			    bufV[k] = bufR[j] = dna[i] ;
			    RbufV[k] = RbufR[j] = myRdna(i) ;
			    k++ ; j++ ;
			  }
			else if (i >= am1 && i < am2 - 1)
			  {
			    bufR[j] = ace_upper (dna[i]) ;
			    RbufR[j] = ace_upper (myRdna(i)) ;
			    j++ ;
			    if (bufN[kk]) 
			      { RbufV[k] = bufV[k] = bufN[kk++] ; k++ ; } 
			    else
			      {
				bufV[k]  = dna[i+da-kk] ;		    
				RbufV[k] = myRdna(i+da-kk) ;
				k++ ;
			      }
			  }
			else if (i >= am2 -1)
			  {
			    if (j < 51) { bufR[j] = dna[i] ; RbufR[j] = myRdna(i) ;j++ ; }
			    bufV[k]  = dna[i+da-kk] ;		    
			    RbufV[k] = myRdna(i+da-kk) ;
			    k++ ;
			  }
		      }
		  RbufR[j] = bufR[j] = 0 ;
		  RbufV[k] = bufV[k] = 0 ;
		}
	    }
	}
      else if (ac_has_tag (Snp, "Insertion"))
	{  
	  if (! ac_has_tag (Snp, "Multi_insertion"))   /* happens because DupA is under InsA and should be under Duplication */
	    {
	      for (ins = inss ; *ins ; ins++)
		if (ac_has_tag (Snp, *ins))
		  {
		    int k ;
		    int da = 1 ;
		    int dx = 0 ;
		    BOOL isDup = FALSE ;
		    char buf[da+1] ;
		    int slide = 0 ;
		    int am1 = fromMrna ? m1 : a1 ;
		    isDup = tctSlideDup (dna, dnaLn, am1, da, &dx, &slide, (*ins) + 3, buf) ;
		    if (a1 < a2) { a1 += dx ; a2 += dx ;}
		    else { a1 -= dx ; a2 -= dx ;}
		    if (fromMrna) { m1 += dx ; m2 += dx ; }
		    am1 = fromMrna ? m1 : a1 ;
		    ok = TRUE ;
		    ok = TRUE ;
		    fs = 1 ;

		    vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Base %c is inserted between base %d and %d\"\n", name (seq), a1, a2, (*ins)[3], a1, a2) ;
		    vtxtPrintf (txt, "VCF %d %c %c%c\n", a1, dna[am1-1], dna[am1-1], ins[0][3]) ;
		    vtxtPrintf (txt, "-D Sliding\n") ;
		    vtxtPrint (txt, "-D Insertion\n") ; /* cleanup */
		    vtxtPrintf (txt,"%s\n", *ins) ; /* reinstate */
		    if (isDup)
		      {
			vtxtPrintf (txt, "gName \"g.%d_%ddup%c\"\n", a1, a2, ins[0][3]) ;
			vtxtPrintf (txt, "Typ \"Dup%c\"\n", ins[0][3]) ;
			if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddup%c\"\n", name(mrna), socrate (g2m+a1), socrate (g2m+a2), ins[0][3]) ;
			vtxtPrintf (txt, "Duplication\n") ;
		      }
		    else
		      {
			vtxtPrintf (txt, "gName \"g.%d_%dins%c\"\n", a1, a2, ins[0][3]) ;
			if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%dins%c\"\n", name(mrna), socrate (g2m+a1), socrate (g2m+a2), ins[0][3]) ;
			vtxtPrintf (txt, "Typ \"Ins%c\"\n", ins[0][3]) ;
		      }
		    if (slide) 
		      vtxtPrintf (txt, "Sliding %d\n", slide) ;
		    
		    for (i = am1 - 21, j = k = 0 ; i < am1 + 30 && i < dnaLn ; i++)
		      if (i > 0 && j < 51) 
			{
			  bufR[j] = dna[i] ; RbufR[j] = myRdna(i) ;  j++ ;
			  if (k < 51) { bufV[k] = dna[i] ; RbufV[k] = myRdna(i) ; k++ ; }
			  if (i == am1-1)
			    { 
			      bufR[j] = '^' ; RbufR[j] = '^' ; j++ ;
			      ccp = buf ; while (*ccp) { bufV[k] = ace_upper(*ccp) ; RbufV[k] = ace_upper(*ccp) ; k++ ; ccp++ ; }
			    }
			}
		    RbufR[j] = bufR[j] = 0 ;
		    RbufV[j] = bufV[j] = 0 ;
		  }
	    }
	  else
	    {
	      tbl = ac_tag_table (Snp, "Multi_insertion", h) ;
	      if (tbl)
		{
		  int am1, da = ac_table_int (tbl, 0, 0, 0) ;
		  const char *ccp = ac_table_printable (tbl, 0, 1, 0) ; 
		  if (0 && ! ccp) {
		    ccp = ac_name (Snp) ;
		    if (ccp) ccp = strchr (ccp, ':') ;
		    if (ccp) ccp = strchr (ccp+1, ':') ;
		    if (ccp) ccp = strchr (ccp+1, ':') ;
		    if (ccp) ccp = strchr (ccp+1, ':') ;
		    if (ccp) ccp = ccp+2 ;
		    da = strlen (ccp) ;
		  }
		  if (ccp && da == strlen (ccp) && a2 - a1 == 1)
		    {
		      int k, dx = 0 ;
		      BOOL isDup = FALSE ;
		      char buf[da+1], bufN[12] ;
		      int slide = 0 ;
		      am1 = fromMrna ? m1 : a1 ;
		      isDup = tctSlideDup (dna, dnaLn, am1, da, &dx, &slide, ccp, buf) ;

		      if (a1 < a2) { a1 += dx ; a2 += dx ;}
		      else { a1 -= dx ; a2 -= dx ;}
		      if (fromMrna) { m1 += dx ; m2 += dx ; }
		      fs = da ;
		      am1 = fromMrna ? m1 : a1 ;
		      
		      vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Between base %d and %d %d bases are inserted\"\n", name (seq), a1, a2, a1, a2, da) ;
		      vtxtPrintf (txt, "-D Sliding\n") ;
		      vtxtPrint (txt, "-D Insertion\n") ; /* cleanup */
		      vtxtPrintf (txt,"Multi_insertion %d\n", da) ; /* reinstate */ 
		      
		      ok = TRUE ;
		      
		      if (da < 20)
			{
			  vtxtPrintf (txt, "VCF %d %c %c%s\n", a1, dna[am1-1], dna[am1-1], buf) ;
			  if (isDup)
			    {
			      vtxtPrintf (txt, "Typ \"Dup%s\"\n", buf) ;
			      vtxtPrintf (txt, "gName \"g.%d_%ddup%s\"\n", a1,a2, buf) ;
			      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddup%s\"\n", name(mrna), socrate (g2m+a1),socrate (g2m+a2), buf) ;
			      vtxtPrintf (txt, "Duplication\n") ;
			    }
			  else
			    {
			      vtxtPrintf (txt, "Typ \"Ins%s\"\n", buf) ;
			      vtxtPrintf (txt, "gName \"g.%d_%dins%s\"\n", a1,a2, buf) ;
			      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%dins%s\"\n", name(mrna), socrate (g2m+a1),socrate (g2m+a2), buf) ;
			    }
			  vtxtPrintf (txt,"Multi_insertion %d %s\n", da, buf) ; /* reinstate */ 
			}
		      else
			{
			  vtxtPrintf (txt, "Typ \"Ins%d\"\n", da) ;
			  vtxtPrintf (txt, "VCF %d %c%d %c\n", a1, dna[am1-1], da, dna[am1-1]) ;
			  vtxtPrintf (txt, "gName \"g.%d_%dins(%d)\"\n", a1,a2, da) ;
			  if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%dins(%d)\"\n", name(mrna), socrate (g2m+a1),socrate (g2m+a2), da) ;
			  vtxtPrintf (txt,"Multi_insertion %d\n", da) ; /* reinstate */ 
			}
		      if (slide) 
			vtxtPrintf (txt, "Sliding %d\n", slide) ;
		      if (da == 1) sprintf (bufN, "^") ; 
		      else if (da == 2) sprintf (bufN, "^^") ; 
		      else  sprintf (bufN, "^%d^", da) ; 
		      
		      for (i = am1 - 21, j = k = 0 ; i < am1 + 30 && i < dnaLn ; i++)
			if (i > 0 && j < 51)
			  {
			    if (i < am1)
			      { bufV[k] = bufR[j] = dna[i] ; RbufV[k] = RbufR[j] = dna[i] ; j++ ; k++ ; }
			    else if (i == am1)
			      {
				if (da)
				  {
				    int kk ;
				    for (kk = 0 ; kk < da && kk < 30 ; kk++)
				      { 
					bufV[j] = RbufV[j] = '-'; j++ ;
					bufR[k] = RbufR[k] = ace_upper(buf[kk]) ; k++ ; 
				      }
				  }
			      }
			    else if (i > am1)
			      {
				if (j < 51) { bufR[j] = dna[i] ; RbufR[j] = myRdna(i) ; j++ ; }
				if (k < 51) { bufV[k] = dna[i] ; RbufV[k] = myRdna(i) ; k++ ; }
			      }
			  }
		      RbufR[j] = bufR[j] = 0 ;
		      RbufV[k] = bufV[k] = 0 ;
		    }
		}
	    }
	}
      else if (ac_has_tag (Snp, "DelIns"))
	{
	  int am1 = 0, am2 = 0 ;
	  tbl = ac_tag_table (Snp, "DelIns", h) ;
	  if (tbl)
	    {
	      int dD = ac_table_int (tbl, 0, 0, 0) ;
	      int dI = ac_table_int (tbl, 0, 1, 0) ;
	      const char *ccD = ac_table_printable (tbl, 0, 2, 0) ; 
	      const char *ccI = ac_table_printable (tbl, 0, 3, 0) ; 

	      fs = dI - dD ;

	      vtxtPrintf (txt, "-D IntMap\nIntMap %s %d %d \"Bases %d to %d (%d bases) are replaced by %d bases\"\n", name (seq), a1, a2, a1, a2, dD, dI) ;
	      if (tbl->cols >= 4 && ccI)
		{
		  int i, j, k ;
		  ok = TRUE ;
		  if (dD < 20 && dI < 20)
		    {
		      am1 = fromMrna ? m1 : a1 ;
		      vtxtPrintf (txt, "Typ \"DelIns%s>%s\"\n", ccD, ccI) ;
		      vtxtPrintf (txt, "VCF %d %c%s %c%s\n", a1, dna[am1-1], ccD, dna[am1-1], ccI) ;
		      vtxtPrintf (txt, "gName \"g.%d_%ddelins%s\"\n", a1,a2, ccI) ;
		      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddelins%s\"\n", name(mrna), socrate (g2m+a1),socrate (g2m+a2), ccI) ;
		    }
		  else
		    {
		      vtxtPrintf (txt, "Typ \"DelIns_%d_%d\"\n", dD, dI) ;
		      vtxtPrintf (txt, "VCF %d %c%s %c%s\n", a1, dna[am1-1], ccD, dna[am1-1], ccI) ;
		      vtxtPrintf (txt, "gName \"g.%d_%ddelins(%d,%d)\"\n", a1,a2, dD, dI) ;
		      if (RdnaLn)  vtxtPrintf (txt, "rName \"%s:c.%d_%ddelins(%d,%d)\"\n", name(mrna), socrate (g2m+a1),socrate (g2m+a2), dD, dI) ;
		    }
		  dda = -1 ;
		  am1 = fromMrna ? m1 : a1 ;
		  am2 = fromMrna ? m2 : a2 ;
		  for (i = am1 - 21 + dda, j = k = 0 ; i < am1 + 30 + dda && i < dnaLn ; i++)
		    if (i > 0 && j < 51) 
		      {
			if (i < am1 - 1 || i >= am2-1) 
			  {
			    if (j < 51)  { bufR[j] = dna[i] ; RbufR[j] = myRdna(i) ; j++ ; }
			    if (k < 51)  { bufV[k] = dna[i] ; RbufV[k] = myRdna(i) ; k++ ; }
			  }
			else if (i == am1 - 1) /* copy the 2 buffers */
			  {
			    ccp = ccD ;
			    while (*ccp && j < 51)
			      { bufR[j] = RbufR[j] = ace_upper(*ccp++) ; j++ ; }
			    ccp = ccI ;
			    while (*ccp && k < 51)
			      { bufV[k] = RbufV[k] = ace_upper(*ccp++) ; k++ ; }
			  }
		      }
		  RbufR[j] = bufR[j] = 0 ;
		  RbufV[k] = bufV[k] = 0 ;
		}
	    }
	}
      
      if (! fromMrna)
	{
	  if (bufR[0]) vtxtPrintf (txt, "Reference_genomic_sequence %s\n", bufR) ;
	  if (0 && bufV[0]) vtxtPrintf (txt, "-D Observed_genomic_sequence %s\n", bufV) ; /* old spelling */
	  if (bufV[0]) vtxtPrintf (txt, "Observed__genomic_sequence %s\n", bufV) ;
	}

      if (RdnaLn)
	{
	  AC_OBJ Mrna = ac_tag_obj (Snp, "mRNA", h) ;
	  int p1 = 0, p2 = 0 ;
	  KEY product = 0 ;
	  if (Mrna)
	    {
	      AC_TABLE pp = ac_tag_table (Mrna, "Product", h) ;
	      int ir ;
	      for (ir = 0 ; pp && ir < pp->rows ; ir++)
		{
		  p1 = ac_table_int (pp, ir, 1 ,0) ;
		  p2 = ac_table_int (pp, ir, 2 ,0) ;
		  if (m1 > p1 && m1 < p2)
		    { product = ac_table_key (pp, ir, 0, 0) ; break ; }
		}
	    }
	  if (product && RbufR[0] && RbufV[0])
	    {
	      /* transmit the coord of the conserved hook base written at offset 21 of the snippet */
	      char pR1[60], pV1[60] ;
	      char pR3[60], pV3[60] ;
	      char *cp, *cq ;
	      int frame, k ;
	      vTXT pp = vtxtHandleCreate (h) ;

	      k = strlen (RbufR) ;
	      for (cp = RbufR + k - 1 ; *cp == 'n' && k > 0 ; cp--, k--)
		*cp = 0 ;
	      k = strlen (RbufV) ;
	      for (cp = RbufV + k - 1 ; *cp == 'n' && k > 0 ; cp--, k--)
		*cp = 0 ;

	      k = strlen (RbufR) ;
	      if (k < 52)
		{
		  for (cp = RbufR + k, cq = bufR + k ; *cq && k < 51 ; cp++, cq++, k++)
		    *cp = *cq ;
		  *cp = 0 ;
		}
	      k = strlen (RbufV) ;
	      if (k < 52)
		{
		  for (cp = RbufV + k, cq = bufV + k ; *cq && k < 51 ; cp++, cq++, k++)
		    *cp = *cq ;
		  *cp = 0 ;
		}
	      
	      /* in RNA, use U and possibly lower the common letters */
	      for (cp = RbufR ; *cp ; cp++)
		if (*cp == 'T') *cp = 'U' ;
	      for (cp = RbufV ; *cp ; cp++)
		if (*cp == 'T') *cp = 'U' ;
	      if (0)
		for (cp = RbufR, cq = RbufV ; *cp ; cp++, cq++) 
		  if (ace_lower(*cp) == ace_lower (*cq)) *cp = *cq = ace_lower (*cp) ;
	      
	      frame = myTranslate (RbufR, pR1, pR3, m1, p1, h) ;
	      myTranslate (RbufV, pV1, pV3, m1, p1, h) ;

	      for (k = 0, cp = pR1, cq = pV1 ; *cp && *cq ; cp++, cq++)
		if (*cp != *cq) 
		  { *cp = ace_lower (*cp) ; *cq = ace_lower (*cq) ; }
	      vtxtPrintf (txt, "Reference_RNAexon_sequence %s\n", RbufR) ;
	      vtxtPrintf (txt, "Observed__RNAexon_sequence %s\n", RbufV) ;

	      tctSetPName (pp, product, pR1, pV1,pR3, pV3, m1, p1, frame, fs) ;
	      cp = vtxtPtr (pp) ;
	      if (pR1[0]) { char *cq = pR1; while (*cq == 'X') cq++ ; vtxtPrintf (txt,   "Reference_protein_sequence %s %s\n", pR3, cq) ; }
	      if (pV1[0]) { char *cq = pV1; while (*cq == 'X') cq++ ; vtxtPrintf (txt,   "Observed__protein_sequence %s %s %s\n", pV3, cq, cp ? cp : "") ; }
	    }
	}
      if (0)
	{
	  char PbufR [51] ;
	  char PbufV [51] ;
	  char geneBuf[1000] ; 
	  KEY product = 0, geneBox = 0 ;
	  int r1 = 0 ;

	  memset (PbufR, 0, sizeof (PbufR)) ;
	  memset (PbufV, 0, sizeof (PbufV)) ;
	  memset (geneBuf, 0, sizeof (geneBuf)) ;
	  
	  if (product)
	    {
	      int px1 = 0, pa1 = 0, dx = pa1 % 3 ;
	      sprintf (geneBuf, "%s %s %d", name (geneBox), name (product), px1) ;
	      memcpy (PbufV, RbufV + dx, r1 - dx) ;
	      memcpy (PbufR, RbufR + dx, r1 - dx) ;

	      /* translate PbufV, PbufR */
	    }
	  if (bufR[0]) vtxtPrintf (txt, "Reference_sequence %s %s %s\n", RbufR, PbufR, geneBuf)  ;
	  if (0 && bufV[0]) vtxtPrintf (txt, "-D Observed__sequence \n") ; /* old spelling */
	  if (bufV[0]) vtxtPrintf (txt, "Observed__sequence %s %s %s\n", RbufV, PbufV, geneBuf) ;
	}
    }

  return ok ? 1 : 0 ;
}

/*************************************************************************************/
/* enter all types of modifs in a logical order */
static DICT *tctMakeVarTypeDict (AC_HANDLE h) 
 {
   DICT *dict = dictHandleCreate (256, h) ;
   char *cp, *cq ;
   /*
   char *Types = "Any,Substitution,Transition,Transversion,Insertion,Deletion,Double insertion,Double deletion,Triple insertion,Triple deletion,Other,A>G,T>C,G>A,C>T,A>T,T>A,G>C,C>G,A>C,T>G,G>T,C>A,Ins A,Ins T,Ins G,Ins C,Del A,Del T,Del G,Del C,Ins AA,Ins TT,Ins GG,Ins CC,Ins AG,Ins CT,Ins AC,Ins GT,Ins TG,Ins CA,Ins TC,Ins GA,Ins AT,Ins TA,Ins GC,Ins CG,Del AA,Del TT,Del GG,Del CC,Del AG,Del CT,Del AC,Del GT,Del TG,Del CA,Del TC,Del GA,Del AT,Del TA,Del GC,Del CG,Ins AAA,Ins TTT,Ins GGG,Ins CCC,Ins AAT,Ins ATT,Ins AAG,Ins CTT,Ins AAC,Ins GTT,Ins TTA,Ins TAA,Ins TTG,Ins CAA,Ins TTC,Ins GAA,Ins GGA,Ins TCC,Ins GGT,Ins ACC,Ins GGC,Ins GCC,Ins CCA,Ins TGG,Ins CCT,Ins AGG,Ins CCG,Ins CGG,Ins ATA,Ins TAT,Ins ATG,Ins CAT,Ins ATC,Ins GAT,Ins AGA,Ins TCT,Ins AGT,Ins ACT,Ins AGC,Ins GCT,Ins ACA,Ins TGT,Ins ACG,Ins CGT,Ins TAG,Ins CTA,Ins TAC,Ins GTA,Ins TGA,Ins TCA,Ins TGC,Ins GCA,Ins TCG,Ins CGA,Ins GAG,Ins CTC,Ins GAC,Ins GTC,Ins GTG,Ins CAC,Ins GCG,Ins CGC,Ins CAG,Ins CTG,Del AAA,Del TTT,Del GGG,Del CCC,Del AAT,Del ATT,Del AAG,Del CTT,Del AAC,Del GTT,Del TTA,Del TAA,Del TTG,Del CAA,Del TTC,Del GAA,Del GGA,Del TCC,Del GGT,Del ACC,Del GGC,Del GCC,Del CCA,Del TGG,Del CCT,Del AGG,Del CCG,Del CGG,Del ATA,Del TAT,Del ATG,Del CAT,Del ATC,Del GAT,Del AGA,Del TCT,Del AGT,Del ACT,Del AGC,Del GCT,Del ACA,Del TGT,Del ACG,Del CGT,Del TAG,Del CTA,Del TAC,Del GTA,Del TGA,Del TCA,Del TGC,Del GCA,Del TCG,Del CGA,Del GAG,Del CTC,Del GAC,Del GTC,Del GTG,Del CAC,Del GCG,Del CGC,Del CAG,Del CTG,Ambiguous" ; 
   */
   char *Types = "Substitution,Transition,Transversion,Deletion,Single_deletion,Double_deletion,Triple_deletion,Deletion_4_30,Long_deletion,Insertion,Single_insertion,Double_insertion,Triple_insertion,Insertion_4_30,Long_insertion,a>g,t>c,g>a,c>t,a>t,t>a,g>c,c>g,a>c,t>g,g>t,c>a,+a,*+a,+t,*+t,+g,*+g,+c,*+c,-a,*-a,-t,*-t,-g,*-g,-c,*-c,++aa,++tt,++gg,++cc,++ag,++ct,++ac,++gt,++tg,++ca,++tc,++ga,++at,++ta,++gc,++cg,*++aa,*++tt,*++gg,*++cc,*++ag,*++ct,*++ac,*++gt,*++tg,*++ca,*++tc,*++ga,*++at,*++ta,*++gc,*++cg,--aa,--tt,--gg,--cc,--ag,--ct,--ac,--gt,--tg,--ca,--tc,--ga,--at,--ta,--gc,--cg,*--aa,*--tt,*--gg,*--cc,*--ag,*--ct,*--ac,*--gt,*--tg,*--ca,*--tc,*--ga,*--at,*--ta,*--gc,*--cg,+++aaa,+++ttt,+++ggg,+++ccc,+++aat,+++att,+++aag,+++ctt,+++aac,+++gtt,+++tta,+++taa,+++ttg,+++caa,+++ttc,+++gaa,+++gga,+++tcc,+++ggt,+++acc,+++ggc,+++gcc,+++cca,+++tgg,+++cct,+++agg,+++ccg,+++cgg,+++ata,+++tat,+++atg,+++cat,+++atc,+++gat,+++aga,+++tct,+++agt,+++act,+++agc,+++gct,+++aca,+++tgt,+++acg,+++cgt,+++tag,+++cta,+++tac,+++gta,+++tga,+++tca,+++tgc,+++gca,+++tcg,+++cga,+++gag,+++ctc,+++gac,+++gtc,+++gtg,+++cac,+++gcg,+++cgc,+++cag,+++ctg,*+++aaa,*+++ttt,*+++ggg,*+++ccc,*+++aat,*+++att,*+++aag,*+++ctt,*+++aac,*+++gtt,*+++tta,*+++taa,*+++ttg,*+++caa,*+++ttc,*+++gaa,*+++gga,*+++tcc,*+++ggt,*+++acc,*+++ggc,*+++gcc,*+++cca,*+++tgg,*+++cct,*+++agg,*+++ccg,*+++cgg,*+++ata,*+++tat,*+++atg,*+++cat,*+++atc,*+++gat,*+++aga,*+++tct,*+++agt,*+++act,*+++agc,*+++gct,*+++aca,*+++tgt,*+++acg,*+++cgt,*+++tag,*+++cta,*+++tac,*+++gta,*+++tga,*+++tca,*+++tgc,*+++gca,*+++tcg,*+++cga,*+++gag,*+++ctc,*+++gac,*+++gtc,*+++gtg,*+++cac,*+++gcg,*+++cgc,*+++cag,*+++ctg,---aaa,---ttt,---ggg,---ccc,---aat,---att,---aag,---ctt,---aac,---gtt,---tta,---taa,---ttg,---caa,---ttc,---gaa,---gga,---tcc,---ggt,---acc,---ggc,---gcc,---cca,---tgg,---cct,---agg,---ccg,---cgg,---ata,---tat,---atg,---cat,---atc,---gat,---aga,---tct,---agt,---act,---agc,---gct,---aca,---tgt,---acg,---cgt,---tag,---cta,---tac,---gta,---tga,---tca,---tgc,---gca,---tcg,---cga,---gag,---ctc,---gac,---gtc,---gtg,---cac,---gcg,---cgc,---cag,---ctg,*---aaa,*---ttt,*---ggg,*---ccc,*---aat,*---att,*---aag,*---ctt,*---aac,*---gtt,*---tta,*---taa,*---ttg,*---caa,*---ttc,*---gaa,*---gga,*---tcc,*---ggt,*---acc,*---ggc,*---gcc,*---cca,*---tgg,*---cct,*---agg,*---ccg,*---cgg,*---ata,*---tat,*---atg,*---cat,*---atc,*---gat,*---aga,*---tct,*---agt,*---act,*---agc,*---gct,*---aca,*---tgt,*---acg,*---cgt,*---tag,*---cta,*---tac,*---gta,*---tga,*---tca,*---tgc,*---gca,*---tcg,*---cga,*---gag,*---ctc,*---gac,*---gtc,*---gtg,*---cac,*---gcg,*---cgc,*---cag,*---ctg" ; 

   for (cp = Types, cq = strchr (Types, '.') ; cp ; cp = cq ? cq + 1 : 0, cq = cp ? strchr (cp, ',') : 0) 
     {
       char cc = cq ? *cq : 0 ;
       dictAdd (dict, cp, 0) ;
       if (cq) *cq = cc ;
     }
   return dict ;
 }

/*************************************************************************************/
/* establish VCF name, gName, rNam, pNam */
static void tctDbTranslate (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tct->outFileName, ".translate.ace", tct->gzo, h) ;
  AC_ITER iter = ac_query_iter (tct->db, TRUE, "Find Variant", 0, h) ;
  int nn = 0, nm = 0, nt = 0 ;
  AC_OBJ Snp = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  
  tct->varTypeDict = tctMakeVarTypeDict (h) ;

  /* NOT DONE: split the 729 multi sub and force create the individual pieces */
  while (ac_free (Snp), Snp = ac_iter_obj (iter))
    {
      nn++ ;
      vtxtClear (txt) ;
      vtxtPrintf (txt, "Variant \"%s\"\n", ac_name (Snp)) ;
      tctSetGName (txt, tct, Snp, h) ;
      /*
      tctSetRname (txt, tct, Snp, h) ;
      tctSetRNam (Snp) ;

      */
      vtxtPrint (txt,"\n") ;
      aceOutf (ao, "%s\n", vtxtPtr (txt)) ;
    }

  fprintf (stderr, "Found %d SNPs, remapped %d, translated %d\n", nn, nm, nt) ;

  ac_free (h) ;
} /* tctDbTranslate */
  
/*************************************************************************************/
/*************************************************************************************/

static int tctMakeWords_Any (TCT *tct, ACEOUT ao, BOOL isMrna)
{
  AC_HANDLE h1 =0, h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  KEYSET ksIn = 0 ;
  const char *zone = tct->zone ? tct->zone : "" ;
  const char *errors = 0 ;
  char *bqlQuery = 0 ;
  BQL_ITER *bqlIter = 0 ;
  int nn = 0, nn1 = 0 ;

  if (isMrna)
    {
      qq = hprintf (h, "%s %s %s %s " 
		    , "Find variant ! Found_in_genome && mRNA" 
		    , tct->filter ? " && (" : ""
		    , tct->filter ? tct->filter : ""
		    , tct->filter ? ") " : ""
		    ) ;
      bqlQuery = hprintf (h, "%s %s"
			  , " select v,ref,var,del, mdel, ins, mins   "
			  , "from v in @, ref in v->Reference_RNAexon_sequence , var in v->Observed__RNAexon_sequence where ref and var, del in v->deletion, mdel in v->multi_deletion, ins in v->insertion, mins in v->multi_insertion"
			  ) ;
    }
  else
    {
      qq = hprintf (h, "%s %s %s %s"
		    , "Find variant Found_in_genome" 
		    , tct->filter ? " && (" : ""
		    , tct->filter ? tct->filter : ""
		    , tct->filter ? ") " : ""
		    ) ;

      bqlQuery = hprintf (h, "%s %s "
			  , " select v,ref,var,del, mdel, ins, mins   "
			  , "from v in @, ref in v->Reference_genomic_sequence , var in v->Observed__genomic_sequence where ref and var, del in v->deletion, mdel in v->multi_deletion, ins in v->insertion, mins in v->multi_insertion"
			  ) ;
    }
  ksIn = query (0, qq) ;
  if (keySetMax (ksIn))   bqlIter = bqlIterCreate (bqlQuery, ksIn, &errors, h) ;
  if (bqlIter) while (ac_free (h1), (tbl = bqlIterTable (bqlIter)))
    {
      KEY v = ac_table_key (tbl, 0, 0, 0) ;
      const char *ref = ac_table_printable (tbl, 0, 1, "xxx") ;
      const char *var = ac_table_printable (tbl, 0, 2, "xxx") ;
      const char *del = ac_table_printable (tbl, 0, 3, 0) ;
      const char *ins = ac_table_printable (tbl, 0, 5, 0) ;
      int mdel = ac_table_int (tbl, 0, 4, 0) ;
      int mins = ac_table_int (tbl, 0, 6, 0) ;
      char buf[128] ;

      nn1++ ;
      if (ref && strlen (ref) > 30 && var && strlen (var) > 30)
	{
	  if (ins && mins == 0)
	    { 
	      memcpy (buf, ref + 6, 15) ;
	      buf[15]='@' ;
	      memcpy (buf+15, ref + 22, 16) ;
	    }
	  else if (ins && mins)
	    {
	      const char *cr = ref + strlen(ref) - 1 ;
	      memcpy (buf, ref + 6, 15) ;
	      while (cr > ref + 20 && *cr != '^') 
		cr-- ;
	      memcpy (buf+15, cr+1, 16) ;
	    }
	  else
	    memcpy (buf, ref + 5 + (del ? 1 : 0), 31) ;
	  buf[32] = 0 ;
	  aceOutf (ao, "%s\tRef__%s\t%s\n", buf, name(v), zone) ;
	  if (del && mdel == 0)
	    { 
	      memcpy (buf, var + 6, 15) ;
	      memcpy (buf+15, var + 22, 16) ;
	    }
	  else if (del && mdel)
	    {
	      const char *cr = var + strlen(var) - 1 ;
	      memcpy (buf, var + 6, 15) ;
	      while (cr > var + 20 && *cr != '_') 
		cr-- ;
	      memcpy (buf+15, cr+1, 16) ;
	    }
	  else
	    memcpy (buf, var + (ins ? 6 : 5), 31) ;
	  buf[32] = 0 ;
	  aceOutf (ao, "%s\tVar__%s\t%s\n", buf, name(v), zone) ;
	  nn++ ;
	}
      if (nn < 0) break ;
    }

  ac_free (ksIn) ;
  ac_free (h) ;
  fprintf (stderr, "# tctMakeWords_Sub isMrna %s exported %d/%d snps\n", isMrna ? "TRUE" : "FALSE", nn, nn1) ;
  return nn ;
} /* tctMakeWords_Any */

/*************************************************************************************/

static int tctMakeWords_Sub (TCT *tct, ACEOUT ao, BOOL isMrna)
{
  AC_HANDLE h1 =0, h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  KEYSET ksIn = 0 ;
  const char *zone = tct->zone ? tct->zone : "" ;
  const char *errors = 0 ;
  char *bqlQuery = 0 ;
  BQL_ITER *bqlIter = 0 ;
  int nn = 0, nn1 = 0 ;

  if (isMrna)
    {
      qq = hprintf (h, "%s %s %s %s " 
		    , "Find variant substitution && ! multi_substitution && ! Found_in_genome && mRNA" 
		    , tct->filter ? " && (" : ""
		    , tct->filter ? tct->filter : ""
		    , tct->filter ? ") " : ""
		    ) ;
      bqlQuery = hprintf (h, "%s %s"
			  , " select v,typ,dna,a1    "
			  , "from v in @, typ in v->typ, m in v->mRNA, a1 in m[1], dna in DNA(m, a1-15, a1+15) where dna"
			  ) ;
    }
  else
    {
      qq = hprintf (h, "%s %s %s %s"
		    , "Find variant substitution && ! multi_substitution && Found_in_genome" 
		    , tct->filter ? " && (" : ""
		    , tct->filter ? tct->filter : ""
		    , tct->filter ? ") " : ""
		    ) ;

      bqlQuery = hprintf (h, "%s %s "
			  , " select v,typ,dna    "
			  , "from v in @, typ in v->typ, chr in v->intMap, a1 in chr[1], seq in v->Parent_sequence, dna in DNA(seq, a1-15, a1+15) where dna"
			  ) ;
    }
  ksIn = query (0, qq) ;
  if (keySetMax (ksIn))   bqlIter = bqlIterCreate (bqlQuery, ksIn, &errors, h) ;
  if (bqlIter) while (ac_free (h1), (tbl = bqlIterTable (bqlIter)))
    {
      KEY v = ac_table_key (tbl, 0, 0, 0) ;
      const char *typ = ac_table_printable (tbl, 0, 1, "xxx") ;
      char cc = typ[2] ;
      const char *dna = ac_table_printable (tbl, 0, 2, 0) ;
      char buf[32] ;

      nn1++ ;
      if (dna && strlen (dna) == 31)
	{
	  aceOutf (ao, "%s\tRef__%s\t%s\n", dna, name(v), zone) ;
	  memcpy (buf, dna, 32) ;
	  buf[15] = cc ; 
	  aceOutf (ao, "%s\tVariant__%s\t%s\n", buf, name(v), zone) ;
	  nn++ ;
	}
      if (nn < 0) break ;
    }

  ac_free (ksIn) ;
  ac_free (h) ;
  fprintf (stderr, "# tctMakeWords_Sub isMrna %s exported %d/%d snps\n", isMrna ? "TRUE" : "FALSE", nn, nn1) ;
  return nn ;
} /* tctMakeWords_Sub */

/*************************************************************************************/

static int tctMakeWords_Del (TCT *tct, ACEOUT ao, BOOL isMrna)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  KEYSET ksIn = 0 ;
  const char *zone = tct->zone ? tct->zone : "" ;
  const char *errors = 0 ;
  char *bqlQuery = 0 ;
  BQL_ITER *bqlIter = 0 ;
  int nn = 0, nn1 = 0 ;

  if (isMrna)
    {
      qq = "Find variant deletion && ! Found_in_genome && mRNA" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,dnaRD,dnaRA,dnaVD,dnaVA    "
			  , "from v in @, m in v->mRNA, a1 in m[1], a2 in m[2], dnaRD in DNA(m, a1-15, a1+15), dnaRA in DNA(m, a2-15, a2+15) where dnaRD, dnaVD in DNA(m, a1-14, a1), dnaVA in DNA(m, a2, a2+15)  where dnaVA"
			  , tct->filter ? " and (" : ""
			  , tct->filter ? tct->filter : ""
			  , tct->filter ? ") " : ""
			  ) ;
    }
  else
    {
      qq = "Find variant deletion && Found_in_genome" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,dnaRD,dnaRA,dnaVD,dnaVA    "
			  , "from v in @, chr in v->intMap, a1 in chr[1], a2 in chr[2], seq in v->Parent_sequence, dnaRD in DNA(seq, a1-15, a1+15) where dnaRD, dnaRA in DNA(seq, a2-15, a2+15) , dnaVD in DNA(seq, a1-14, a1), dnaVA in DNA(seq, a2, a2+15)  where dnaVA"
			  , tct->filter ? " and (" : ""
			  , tct->filter ? tct->filter : ""
			  , tct->filter ? ") " : ""
			  ) ;
    }
  ksIn = query (0, qq) ;
  if (keySetMax (ksIn)) bqlIter = bqlIterCreate (bqlQuery, ksIn, &errors, h) ;
  if (bqlIter) while (ac_free (h1), (tbl = bqlIterTable (bqlIter)))
    {
      KEY v = ac_table_key (tbl, 0, 0, 0) ;
      const char *dnaRD = ac_table_printable (tbl, 0, 1, "") ;
      const char *dnaRA = ac_table_printable (tbl, 0, 2, "") ;
      const char *dnaVD = ac_table_printable (tbl, 0, 3, "") ;
      const char *dnaVA = ac_table_printable (tbl, 0, 4, "") ;

      nn1++ ;
      if (dnaRD && strlen (dnaRD) == 31)
	{
	  aceOutf (ao, "%s\tRefDonor__%s\t%s\n", dnaRD, name(v), zone) ;
	  aceOutf (ao, "%s\tRefAccep__%s\t%s\n", dnaRA, name(v), zone) ;
	  aceOutf (ao, "%s%s\tVariant__%s\t%s\n", dnaVD,dnaVA,name(v), zone) ;
	  nn++ ;
	}
      if (nn < 0) break ;
    }

  ac_free (ksIn) ;
  ac_free (h) ;
  fprintf (stderr, "# tctMakeWords_Del isMrna %s exported %d/%d snps\n", isMrna ? "TRUE" : "FALSE", nn, nn1) ;
  return nn ;
} /* tctMakeWords_Del */

/*************************************************************************************/

static int tctMakeWords_Ins (TCT *tct, ACEOUT ao, BOOL isMrna)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  KEYSET ksIn = 0 ;
  const char *zone = tct->zone ? tct->zone : "" ;
  const char *errors = 0 ;
  char *bqlQuery = 0 ;
  BQL_ITER *bqlIter = 0 ;
  int nn = 0, nn1 = 0 ;

  if (isMrna)
    {
      qq = "Find variant Insertion && ! Found_in_genome && mRNA" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,dnaR,dnaVD,dnaVA,v2    "
			  , "from v in @, m in v->mRNA, a1 in m[1], a2 in m[2], dnaR in DNA(m, a1-15, a1+15), dnaRD in DNA(m, a1-15, a1) where dnaRD, dnaRA in DNA(m, a1+1, a1+15), v1 in v->VCF[2], v2 in v1[1] where v2"
			  , tct->filter ? " and (" : ""
			  , tct->filter ? tct->filter : ""
			  , tct->filter ? ") " : ""
			  ) ;
    }
  else
    {
      qq = "Find variant Insertion && Found_in_genome" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,dnaR,dnaVD,dnaVA,v2    "
			  , "from v in @, chr in v->intMap, a1 in chr[1], a2 in chr[2], seq in v->Parent_sequence, dnaR in DNA(seq, a1-15, a1+15), dnaVD in DNA(seq, a1-15, a1) where dnaVD, dnaVA in DNA(seq, a1+1, a1+15), v1 in v->VCF[2], v2 in v1[1] where v2"
			  , tct->filter ? " and (" : ""
			  , tct->filter ? tct->filter : ""
			  , tct->filter ? ") " : ""
			  ) ;
    }
  ksIn = query (0, qq) ;
  if (keySetMax (ksIn))   bqlIter = bqlIterCreate (bqlQuery, ksIn, &errors, h) ;
  if (bqlIter) 
    while (ac_free (h1), (tbl = bqlIterTable (bqlIter)))
      {
	KEY v = ac_table_key (tbl, 0, 0, 0) ;
	const char *dnaR = ac_table_printable (tbl, 0, 1, "") ;
	const char *dnaVD = ac_table_printable (tbl, 0, 2, "") ;
	const char *dnaVA = ac_table_printable (tbl, 0, 3, "") ;
	const char *v2 = ac_table_printable (tbl, 0, 4, "") ;
	char buf[64], *cp ; const char *cq ;
	int i = 0 ;
	
	nn1++ ;
	i = 0 ; sscanf (v2, "%d", &i) ; if (i>0) continue ; /* we scanned a number */ 
	if (v2 && v2[0]) v2++ ; /* remove the hook base */
	if (dnaR && strlen (dnaR) == 31 && v2 && v2[0])
	  {
	    i = 0 ;
	    cp = buf ; cq = dnaVD ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    cq = v2 ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    cq = dnaVA ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; }
	    *cp = 0 ; buf[32] = 0 ;
	    if (strlen(buf)==31)
	      {
		aceOutf (ao, "%s\tRef__%s\t%s\n", dnaR, name(v), zone) ;
		aceOutf (ao, "%s\tVariant__%s\t%s\n", buf,name(v), zone) ;
		nn++ ;
	      }
	  }
	if (nn < 0) break ;
      }
  
  ac_free (ksIn) ;
  ac_free (h) ;
  fprintf (stderr, "# tctMakeWords_Ins isMrna %s exported %d/%d snps\n", isMrna ? "TRUE" : "FALSE", nn, nn1) ;
  return nn ;
} /* tctMakeWords_Ins */

/*************************************************************************************/

static int tctMakeWords_MultiSub (TCT *tct, ACEOUT ao, BOOL isMrna)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  char *qq ;
  KEYSET ksIn = 0 ;
  const char *zone = tct->zone ? tct->zone : "" ;
  const char *errors = 0 ;
  char *bqlQuery = 0 ;
  BQL_ITER *bqlIter = 0 ;
  int nn = 0, nn1 = 0 ;

  if (isMrna)
    {
      qq = "Find variant (Multi_substitution || DelIns) && ! Found_in_genome && mRNA" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,x1,x2,,y,dnaRD,dnaRA,v1,v2    "
			  , "from v in @, x1 in v->DelIns, x2 in x1[1], y in v->Multi_substitution, m in v->mRNA, a1 in m[1], a2 in m[2], dnaR in DNA(m, a1-15, a1+15), dnaRD in DNA(m, a1-15, a1-1) where dnaRD, dnaRA in DNA(m, a2, a2+15), v1 in v->VCF[2], v2 in v1[1] where v2 "
			  , tct->filter ? " and (" : ""
			  , tct->filter ? tct->filter : ""
			  , tct->filter ? ") " : ""
			  ) ;
    }
  else
    {
      qq = "Find variant (Multi_substitution || DelIns)  && Found_in_genome" ;
      bqlQuery = hprintf (h, "%s %s %s %s %s"
			  , " select v,x1,x2,y,dnaRD,dnaRA,v1,v2    "
			  , "from v in @, x1 in v->DelIns, x2 in x1[1], y in v->Multi_substitution, chr in v->intMap, a1 in chr[1], a2 in chr[2], seq in v->Parent_sequence, dnaRD in DNA(seq, a1-15, a1-1), dnaRA in DNA(seq, a2, a2+15) where dnaRA, v1 in v->VCF[2], v2 in v1[1] where v2 "
			  , tct->filter ? " and (" : ""
			  , tct->filter ? tct->filter : ""
			  , tct->filter ? ") " : ""
			  ) ;
    }
  ksIn = query (0, qq) ;
  if (keySetMax (ksIn)) bqlIter = bqlIterCreate (bqlQuery, ksIn, &errors, h) ;
  if (bqlIter) 
    while (ac_free (h1), (tbl = bqlIterTable (bqlIter)))
      {
	KEY v = ac_table_key (tbl, 0, 0, 0) ;
	int x1 = ac_table_int (tbl, 0, 1, 0) ;
	int x2 = ac_table_int (tbl, 0, 2, 0) ;
	int y = ac_table_int (tbl, 0, 3, 0) ;
	
	const char *dnaRD = ac_table_printable (tbl, 0, 4, "") ;
	const char *dnaRA = ac_table_printable (tbl, 0, 5, "") ;
	const char *v1, *v2 ;
	char bufR[64], bufV[64], *cp ; const char *cq ;
	int i ;
	nn1++ ;
	
	v1 = ac_table_printable (tbl, 0, 6, "") ;
	i = 0 ; sscanf (v1, "%d", &i) ; if (i>0) continue ; /* we scanned a number */ 
	if (x1 && strlen (v1) != x1) continue ;
	if (y && strlen (v1) != y) continue ;

	v2 = ac_table_printable (tbl, 0, 7, "") ; /* v1, v2 are not stable, this is a bug of the library that must be fixed */
	i = 0 ; sscanf (v2, "%d", &i) ; if (i>0) continue ; /* we scanned a number */ 
	if (x2 && strlen (v2) != x2) continue ;
	if (y && strlen (v2) != y) continue ;

	if (1) 
	  {
	    i = 0 ;
	    cp = bufR ; cq = dnaRD ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    v1 = ac_table_printable (tbl, 0, 6, "") ;
	    cq = v1 ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    cq = dnaRA ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; }
	    *cp = 0 ; bufR[32] = 0 ;
	    
	    i = 0 ;
	    cp = bufV ; cq = dnaRD ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    v2 = ac_table_printable (tbl, 0, 7, "") ;
	    cq = v2 ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; } 
	    cq = dnaRA ; while (*cq && i < 31) { *cp++ = *cq++ ; i++ ; }
	    *cp = 0 ; bufV[32] = 0 ;
	    if (strlen(bufR)==31 && strlen(bufV)==31)
	      {
		aceOutf (ao, "%s\tRef__%s\t%s\n", bufR, name(v), zone) ;
		aceOutf (ao, "%s\tVariant__%s\t%s\n", bufV,name(v), zone) ;
		nn++ ;
	      }
	  }
	if (nn < 0) break ;
      }
  
  ac_free (ksIn) ;
  ac_free (h) ;
  fprintf (stderr, "# tctMakeWords_Multi isMrna %s exported %d/%d snps\n", isMrna ? "TRUE" : "FALSE", nn, nn1) ;
  return nn ;
} /* tctMakeWords_MultiSub */

/*************************************************************************************/

static void tctMakeWords (TCT *tct)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tct->outFileName, ".w31", tct->gzo, h) ;

  tctMakeWords_Any (tct, ao, TRUE) ;
  tctMakeWords_Any (tct, ao, FALSE) ;
  if (0)
    {
      tctMakeWords_Sub (tct, ao, TRUE) ;
      tctMakeWords_Sub (tct, ao, FALSE) ;
      tctMakeWords_Del (tct, ao, TRUE) ;
      tctMakeWords_Del (tct, ao, FALSE) ;
      tctMakeWords_Ins (tct, ao, TRUE) ; 
      tctMakeWords_Ins (tct, ao, FALSE) ;
      tctMakeWords_MultiSub (tct, ao, TRUE) ;  
      tctMakeWords_MultiSub (tct, ao, FALSE) ;  
    }
  ac_free (h) ;
} /* tctMakeWords */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/
/* run -run NA12878MOD -laneList tmp/TSNP/NA12878MOD/LaneList -t 20 -target_fasta TARGET/CHROMS/hs.chrom_20.fasta.gz -t1 25000001 -t2 30010000 -minSnpFrequency 18 -minSnpCover 10 -minSnpCount 4 -target_class Z_genome -o tata -maxLanes 4
 */

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  if (argc == 0)
  fprintf (stderr,
	   "// Usage: tsnp -db $DB -project $MAGIC -db_report -i[-f] f.val.txt \n"
	   "//      Post processing of .val.txt SNP counts exported by tricoteur -val\n"
	   "//      try: -h --help \n"
	   "// Example\n"
	   "//     tsnp -i f1.val.txt,f2.val.txt -db tmp/TSNP/$zone -p MMM -db_report\n"
	   "// Several sets of parameters must or may be specified\n"
	   "// SNP COUNTS TO BE ANALYZED: \n"
	   "//   -i f[,g[,h]] : coma separated list of files\n"
	   "//   -f file : file contains a list of files, coma or space or tab or line separated\n"
	   "//      ATTENTION all file names must be local or fully qualified starting at /\n"
	   "//      All the input files are of type .val.txt as exported by tricoteur -val\n"
	   "//      The format is very specific and cannot easily be replicated\n"
	   "//      So the present program can only postprocess tricoteur\n"
	   "// TARGET ZONE:\n"
	   "//   -target_class : target class (A_mito ... Z_genome) [default Z_genome]\n"
	   "//   -t targetName :  [optional] often a chromosome name\n"
	   "//   -t1 int -t2 int :  [optional] analyze only section [t1,t2] of the target\n"
	   "//    example:\n"
	   "//        -target_class  Z_genome -t chr7 -t1 2000000 -t2 3000000\n"
	   "//       Typically t2 - t1 = 10 Mbases, or if t1, t2 are not specified\n"
	   "//       the whole chromosome is analyzed .\n"
	   "// SNP FILTERS, DISCOVERY PHASE\n"
	   "//   -method arbitrary_name: the method is echoed in the output file\n" 
  	   "//   -minSnpCover integer : min coverage [default 10] \n"
	   "//   -minSnpCount integer: [default 4]\n"
	   "//   -minSnpFrequency float: [default 18] minimal MAF precentage\n"
	   "//   -intron : [default off] diifferentiate introns from deletions\n"
	   "//   -intron_only : just detect and report introns\n"
	   "//   -min_intron <int> : [default 30] min reported intron length\n"
	   "//   -max_intron <int> : [default 0]  max reported intron length\n"
	   "// Phase 3 actions: Analyse the .snp files\n"
	   "//   The analyses rely on a the existence of an acedb TSNP_DB/$zone database, aware of genes, transcripts and coding structures\n"
	   "//     and containing a copy of the metadata of the runs, originally hand constructed in MetaDB\n"
	   "//     The parameter -db points to this database, we recommend one database per zone, allowing parallelization\n"
	   "//   -db_remap2genome tmp/METADATA/mrnaRemap.gz  -db ACEDB\n"
	   "//      Remap the transcript variants into genome coordinates\n"
	   "//   -db_remap2genes tmp/METADATA/mrnaRemap.gz -db ACEDB\n"
	   "//      Remap the genome variants into transcript coordinates\n"
	   "//   -db_translate -db ACEDB : translate the mRNA variants (or genome variants remapped to mRNAs) if they map to a protein coding exon\n"
	   "//   -db_count -i count_file -db ACEDB  : scan the input file and adjust in the ACEDB database the variant->population and ->strand counts\n"
	   "// GENE FUSION\n"
	   "//   -target_class : target class (KT_RefSeq ET_av...) [default ET_av]\n"
	   "//   -min_GF integer : [default 5]  filter geneFusionFile on min support \n" 
	   "//   -minOverhang integer : [default 15] minimal number of bases\n"
	   "//   -geneFusion fileName: file of genefusions to be analysed\n"
	   "//      mrna1 a1 a2 mrna2 b1 b2 n (n supports for a jump from mrna1[position a2] to m2[b1]\n"
	   "//      Scan the hit file(s) report for each donor/acceptor read count that support\n"
	   "//         the proposed donor and goes to the acceptor\n"
	   "//         OR align locally OR jump locally OR jump elsewhere\n"
	   "// OUTPUT\n"
	   "//     The program exports a vcf file and several histograms\n"
	   "//   -o fileNamePrefix : output file name, equivalent to redirecting stdout\n"
  	   "//   -gzo : the output files should be gzipped\n"
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
  TCT tct ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&tct, 0, sizeof (TCT)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  tct.h = h ;

  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "--help")||
      getCmdLineBool (&argc, argv, "-h")
      )
    usage (commandBuf, 0, argv) ;
  if (getCmdLineBool (&argc, argv, "-v") ||
      getCmdLineBool (&argc, argv, "--version")
      )
    {
      fprintf (stderr, "variant_caller: %s\n", VERSION) ;
      exit (0) ;
    }
 
  {
    int ix ;
    char *cp ;
    for (ix = 0, cp = commandBuf ;  ix < argc && cp + strlen (argv[ix]) + 1 < commandBuf + 3900 ; cp += strlen (cp), ix++)
      sprintf(cp, "%s ", argv[ix]) ;
  }

  /* consume optional args */
  tct.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  tct.gzo = getCmdLineBool (&argc, argv, "-gzo") ;

  getCmdLineOption (&argc, argv, "-o", &(tct.outFileName)) ;
 
  /* TARGET ZONE */
  tct.target_class = 0 ; /* "Z_genome" ; */
  getCmdLineOption (&argc, argv, "-referenceGenome", &(tct.referenceGenome)) ;
  getCmdLineOption (&argc, argv, "-target_class", &(tct.target_class)) ;
  getCmdLineOption (&argc, argv, "-t", &(tct.target)) ;
  getCmdLineInt (&argc, argv, "-t1", &(tct.t1)) ;
  getCmdLineInt (&argc, argv, "-t2", &(tct.t2)) ;

  tct.intron_only = getCmdLineBool (&argc, argv, "-intron_only") ;

  tct.minIntron = 30 ; /* default */
  tct.maxIntron = 0 ; /* default : do not search introns */
  if (getCmdLineBool (&argc, argv, "-intron"))
    tct.maxIntron = 100000 ; /* default if we search introns */
 
  if (tct.t2 < tct.t1)
    {
      fprintf (stderr, "FATAL ERROR: zone limits expeted t1 < t2, got %d > %d, please try\n\tvariant_caller -help\n", tct.t1, tct.t2) ;
      exit (1) ;
    }
  if (tct.t1 < 0)
    {
      fprintf (stderr, "FATAL ERROR: zone limits expeted t1 > 0 got %d  please try\n\tvariant_caller -help\n", tct.t1) ;
      exit (1) ;
    }

  /* SNP FILTERS */
  tct.minSnpCover = 20 ;
  tct.minSnpCount = 4;
  tct.minSnpFrequency = 5 ;

  getCmdLineInt (&argc, argv, "-minSnpCover", &(tct.minSnpCover)) ;
  getCmdLineInt (&argc, argv, "-minSnpCount", &(tct.minSnpCount)) ;
  getCmdLineFloat (&argc, argv, "-minSnpFrequency", &(tct.minSnpFrequency)) ;


  getCmdLineOption (&argc, argv, "-i", &(tct.inFileList)) ;
  getCmdLineOption (&argc, argv, "--inFiles", &(tct.inFileList)) ;
  getCmdLineOption (&argc, argv, "-f", &(tct.inFileOfFileList)) ;
  getCmdLineOption (&argc, argv, "--fileList", &(tct.inFileOfFileList)) ;
  getCmdLineOption (&argc, argv, "--project", &(tct.project)) ;
  getCmdLineOption (&argc, argv, "-p", &(tct.project)) ;
  getCmdLineOption (&argc, argv, "--force", &(tct.force)) ;
  getCmdLineOption (&argc, argv, "--zone", &(tct.zone)) ;
  getCmdLineOption (&argc, argv, "--filter", &(tct.filter)) ;
  getCmdLineOption (&argc, argv, "--", &(tct.filter)) ;
  getCmdLineOption (&argc, argv, "-db_remap2genome", &tct.remap2genome) ;
  getCmdLineOption (&argc, argv, "-db_remap2genes", &tct.remap2genes) ;

  tct.makeWords = getCmdLineBool (&argc, argv, "--makeWords") ;

  tct.nAna = 4 ;
  getCmdLineInt (&argc, argv, "-nAna", &(tct.nAna)) ;
  
  
  /********* REPORT **********/

  getCmdLineOption (&argc, argv, "-db", &(tct.dbName)) ;
  getCmdLineOption (&argc, argv, "-wiggleDir", &(tct.wiggleDir)) ;
  tct.max_threads = 4 ;
  getCmdLineInt (&argc, argv, "-max_threads", &tct.max_threads) ;
  
  if (tct.minSnpCount > tct.minSnpCover)
    tct.minSnpCount =  tct.minSnpCover ;

  tct.mergeCounts = getCmdLineBool (&argc, argv, "-merge") ;
  tct.dbReport = getCmdLineBool (&argc, argv, "-db_report") ;
  tct.dbTranslate = getCmdLineBool (&argc, argv, "-db_translate") ;
  tct.dropMonomodal = getCmdLineBool (&argc, argv, "-dropMonomodal") ;
  if (tct.dbName)
    {
      const char *errors ;
      tct.db = ac_open_db (tct.dbName, &errors);
      if (! tct.db)
	messcrash ("Failed to open db %s, error %s", tct.dbName, errors) ;
    }
  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }
  filAddDir ("./") ;

  if (tct.dbReport || tct.makeWords)
    {
      if (! tct.db)
	{
	  fprintf (stderr, "-db_report requires -db SnpMetaDataDB, sorry, try -help\n") ;
	  exit (1) ;
	}
      if (! tct.project)
	{
	  fprintf (stderr, "-db_report requires -project $MAGIC, sorry, try -help\n") ;
	  exit (1) ;
	}
      if (! tctGetRunList (&tct))
	messcrash ("No run in -db %s belong to project %s", tct.dbName, tct.project) ;
    }
  if (tct.dbTranslate)
    {
      if (! tct.db)
	{
	  fprintf (stderr, "-db_translate requires -db SnpMetaDataDB, sorry, try -help\n") ;
	  exit (1) ;
	}
    }

  /* parallelization */
   if (tct.max_threads < 4)
     tct.max_threads = 4 ;
   if (tct.max_threads < tct.nAna)
     tct.max_threads = tct.nAna ;
 
  /* check the absolute args */

  wego_max_threads (tct.max_threads) ;
  tct.doneChan = channelCreate (12, int, tct.h) ;

  if (tct.makeWords)
    {
      tctMakeWords (&tct) ;
    }
  if (tct.dbReport)
    {
      tctDbReport (&tct) ;
    }
  if (tct.mergeCounts)
    {
      tctMergeCounts (&tct) ;
    }
  if (tct.dbTranslate)
    {
      tctDbTranslate (&tct) ;
    }
  if (tct.remap2genome)
    {
      tctCreateAtlas (&tct) ;
      tctRemap1 (&tct) ; tctRemap2(&tct) ;
    }
  if (tct.remap2genes)
    {
      tctCreateAtlas (&tct) ;
      messcrash ("tctRemap2genes  not programmed") ;
    }
  wego_flush () ;

  ac_free (tct.ao) ;
  if (1)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// minSnpCount %d minSnpCover %d minSnpFrequency %.1f%%\n"
	       , tct.minSnpCount
	       , tct.minSnpCover
	       , tct.minSnpFrequency
	       ) ;
      if (tct.snps)
	fprintf (stderr, "// SNP detected %d SEG reported %d\n"
		 , tct.snpDetected
		 , tct.snpExported
		 ) ;
     }
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderrn and to ensure all pipes are closed*/
  ac_free (tct.h) ;

  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
