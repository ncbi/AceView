/*  File: estscan.c
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
 */

#include "ac.h"
#include "keyset.h"
#include "dna.h"
#include "mytime.h"
#include "freeout.h"
#include "dict.h"
static char B2[256] ;
static char B2r[256] ;
static int nVerif = 0, nEst = 0 ;

typedef struct estScanStruct {
  AC_HANDLE h ;
  const char *estFileName ;
  const char *genomeFileName, *wordFrequencyFile ;
  const unsigned char *dnaGenome ;
  FILE *estFile, *genomeFile ;
  int wordLength, maskLength, minEntropy, genomeStep, minHitNumber, maxWordRepeat;
  int genomeLength ;
  unsigned int mask ;
  Stack genomeStack ;
  Associator ass ;
} ESTSCAN ;

/*************************************************************************************/
/*************************************************************************************/

static void estScanInit ()
{
  int i ; 
  
  i = 256 ;  while (i--) B2[i] = 0 ;
  i = 256 ;  while (i--) B2r[i] = 0 ;
  B2[A_] = 0x0 ; B2[T_] = 0x3 ;     /* ATTENTION copied in dnasubs.c:dnagetWordUsage */
  B2[G_] = 0x1 ; B2[C_] = 0x2 ;     /* you must keep the 2 identical */
  B2r[A_] = 0x3 ; B2r[T_] = 0x0 ;     /* ATTENTION copied in dnasubs.c:dnagetWordUsage */
  B2r[G_] = 0x2 ; B2r[C_] = 0x1 ;     /* you must keep the 2 identical */

  return ;
} /* estScanInit */

/*************************************************************************************/
/*****************************  actual work ******************************************/
/*************************************************************************************/
/* if length > 16 we have to explicitely verify the hit */
static BOOL estScanVerifyEstHit (const unsigned char *dnaGenome, int a1, Array dnaEst, int x1, int wordLength)
{
  const unsigned char *ccp, *ccq ;
  int maxN = 0 ;

  nVerif++ ;
  ccp = dnaGenome + a1 ;
  ccq = arrp (dnaEst, x1, unsigned char) ;

  while (wordLength--)
    {
      if (*ccp != *ccq)
	{
	  if (*ccp & *ccq) maxN-- ;
	  if (maxN < 0)
	    return FALSE ;
	}
      ccp++ ; ccq++ ;
    }
  return TRUE ;
} /* estScanVerifyEstHit */

/*********************************************************************/

static int estScanGetEstHits (ESTSCAN *pp, const char *estName, Array dna)
{
  int ii, pos, dx1, nHits = 0 ;
  int a1 ;
  char *cp ;
  unsigned int oligo, mask ;
  const void *vp;
  Associator ass = pp->ass ;

  mask = pp->mask ;
  dx1 =  pp->maskLength ;

  if (arrayMax (dna) < pp->wordLength)
    return 0 ;
  /* load the mask */
  oligo =  0 ;
  cp = arrp(dna, 0, char) - 1  ;
  ii = dx1 - 1 ;
  while (cp++, ii--)
    { oligo <<= 2 ; oligo |= B2[(int)(*cp)] ; oligo &= mask ; } 
  
  pos = -1 ; ii = arrayMax(dna) - pp->wordLength + 1 ;
  cp = arrp(dna, 0, char) + dx1 - 2  ;
  if (ii > 0) 
    while (pos++, cp++, ii--, ii>0)
      { 
	Array bucket = 0 ;
	int iBucket = 0 ;
	oligo <<= 2 ; oligo |= B2[(int)(*cp)] ; oligo &= mask ;
	
	if (oligo && ass && assFind(ass, assVoid(oligo), &vp))
	  while (assFindNext(ass, assVoid(oligo), &vp, &bucket, &iBucket))
	    {
	      a1 = assInt(vp)-1 ;
	      if (pp->wordLength <= 16 || estScanVerifyEstHit (pp->dnaGenome, a1
				       , dna, pos
				       , pp->wordLength)
		  )
		{ nHits++ ; cp += 20 ; ii -= 20 ; break ; }
	    }
	if (nHits >= pp->minHitNumber)
	  break ;
      }
  
  return nHits ;
} /* estScanGetEstHits */

/*************************************************************************************/

static BOOL estScanSearch (ESTSCAN *pp, const char *estName, Array dna)
{
  int nHits = estScanGetEstHits (pp, estName, dna) ;

  nEst++ ;
  if (nHits < pp->minHitNumber)
    {
      reverseComplement (dna) ;
      nHits += estScanGetEstHits (pp, estName, dna) ;
      reverseComplement (dna) ;
    }
	
  if (nHits >= pp->minHitNumber)
    {
      freeOutf ("Sequence %s\n", freeprotect (estName)) ;
      return TRUE ;
    }
  return FALSE ;
} /* estScanSearch */

/*************************************************************************************/

static BOOL estScanRun (ESTSCAN *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  int n, state, line = 0, nn = 0, seq = 0, nHits = 0 ;
  char *cp, *cq ; 
  char estName [1000] ;
  Stack s = stackHandleCreate (300, h) ;
  Array dna = arrayHandleCreate (100000, unsigned char, h) ;
  ACEIN ai = aceInCreate (pp->estFileName, 0, h) ;

  if (!ai)
    goto done ;
  aceInSpecial (ai, "\n") ;
  state = nn = 0 ;
  
  while (aceInCard (ai)) 
    {
      line++ ;
      cp = aceInWord (ai) ;
      if (cp) switch (state)
	{
	case 1:
	case 2:
	  if (*cp != '>')
	    {
	      if (state == 1)
		{
		  state = 2 ;
		  pushText (s, cp) ;
		}
	      else
		catText (s, cp) ;
	      break ;
	    }
	  else
	    {
	      cq = stackText (s, seq) ;
	      n = strlen (cq) ;
	      array (dna , n, unsigned char) = 0 ; /* make room */
	      arrayMax (dna) = n ;
	      memcpy (arrp (dna, 0, unsigned char), cq, n) ;
	      dnaEncodeArray (dna) ;
	 
	      if (estScanSearch (pp, estName, dna))
		nHits++ ;
	      state = 0 ; /* and fall thru */
	      stackClear (s) ;
	      dna = arrayReCreate (dna,100000, unsigned char) ; 		
	    }
      /* fall thru to new sequence */
	case 0: /* expecting   >target_name */
	  if (*cp != '>' || ! *(cp+1))
	    {
	      fprintf (stderr, "// bad character in file %s line %d, expecting a '>', got %s"
		       , pp->genomeFileName, line, cp) ;
	      return FALSE ;
	    }
	  strcpy (estName, cp + 1) ;
	  state = 1 ;
	  seq = stackMark (s) ;
	  break ;
	}
    }
  if (state > 0)
    {
      cq = stackText (s, seq) ;
      n = strlen (cq) ;
      array (dna , n, unsigned char) = 0 ; /* make room */
      arrayMax (dna) = n ;
      memcpy (arrp (dna, 0, unsigned char), cq, n) ;
      dnaEncodeArray (dna) ;
      if (estScanSearch (pp, estName, dna))
	nHits++ ;
      state = 0 ; /* and fall thru */
      stackClear (s) ;
    }

 done:
  ac_free (h) ;
  return nHits ;
} /* estScanRun */

/*************************************************************************************/
/**********************************************************************/

static int oligoEntropy2 (int nn, int na, int nt, int ng, int nc)
{
  int ss ;
  static int ee[65] ;
  static int oldNn = 0 ;

  if (nn > 16) nn = 16 ;
  if (nn != oldNn)
    {
      double s = 0, log4 = log(4.0) ; /* divide by log4 to be in base equivalent */ 
      int j ;

      oldNn = nn ;
      for (j = 1 ; j <= nn ; j++)
	{ s = j ; s = s/nn ; ee[j] = - (int) (1000.0 * j * log(s)/log4 ) ; }
      ee[0] = 0 ;
    }

  ss = (ee[na] + ee[ng] + ee[nc] + ee[nt]) / 1000 ;
 
  return ss ;
} /* oligoEntropy2  */

/*************************************************************************************/

static int estScanHashGenome (ESTSCAN *pp)
{
  Associator ass0, ass = 0 ;
  int n, i, j, jx = 16, n1 = 0, nw1 = 0, pos, pass ;
  int na, nt, ng, nc ;
  int minEntropy = pp->minEntropy ;
  unsigned int oligo, mask ;
  char *cp ;
  const void *vp ;
  BOOL goodName ;
  unsigned long  wnb = 0, wnw = 0 ;
  Array words = 0 ;

  i = pp->genomeLength / pp->genomeStep ;
  ass = pp->ass = assBigCreate (2 * i) ;
  ass0 = assBigCreate (2 * pp->genomeLength ) ;

  mask = pp->mask ; jx =  pp->maskLength ;

  if (pp->wordFrequencyFile)
    {
      words = dnaGetWordUsage  (0, 1, &wnb, &wnw, FALSE, pp->wordFrequencyFile, 0) ;
      fprintf (stderr, "Found wnd=%ld  wnw=%ldParsed in file %s\n"
	       , wnb, wnw, pp->wordFrequencyFile) ; 
    }
  for (pass = 1 ; pass < 2 ; pass++)
    for (pos = 0 ; pos < pp->genomeLength - pp->wordLength ; pos += pass ? pp->genomeStep : 1)
      {
	cp = stackText (pp->genomeStack, pos) ;
	oligo = 0 ; i = 0 ; goodName = TRUE ;
	
	cp-- ; i = j = 0 ;
	na = nt = ng = nc = 0 ;
	while (j < jx && *++cp)
	  {
	    j++ ;
	    switch (*cp)
	      {
	      case A_:
		na++ ;
		oligo <<= 2 ; oligo |= B2[(int)*cp] ;
		break ;
	      case T_:
		nt++ ;
		oligo <<= 2 ; oligo |= B2[(int)*cp] ;
		break ;
	      case G_:
		ng++ ;
		oligo <<= 2 ; oligo |= B2[(int)*cp] ;
		break ;
	      case C_:
		nc++ ;
		oligo <<= 2 ; oligo |= B2[(int)*cp] ;
		break ;
	      default:
		goodName = FALSE ;
		break ;
	      } 
	  }
	oligo &= mask ;
	
	if (oligo && goodName && oligoEntropy2 (pp->maskLength, na, nt, ng, nc) > minEntropy)
	  {
	    nw1 = 0 ;
	    if (words)
	      {
		unsigned char cc ;
		unsigned int oligo15 = oligo & 0x3fffffff ;
		int pos1 = oligo15 >> 1 ; /* divide by 2 to be in half char */
		int	rest1 = oligo15 & 0x1 ;

		cc =  arr(words, pos1, unsigned char) ;
		nw1 = cc ;
		if (rest1) nw1 >>= 4 ;
		nw1 &= 0x0f ;
	      }

	    if (nw1 < pp->maxWordRepeat)
	      {
		if (pass == 0)
		  {
		    assMultipleInsert (ass0, assVoid(oligo), assVoid(pos + 1)) ;
		  }
		else
		  {
		    Array bucket = 0 ;
		    int iBucket = 0 ;

		    n = 0 ;
		    if (assFind(ass0, assVoid(oligo), &vp))
		      while (assFindNext(ass0, assVoid(oligo), &vp, &bucket, &iBucket))
			n++ ;
		    if (n < pp->maxWordRepeat)
		      {
			if (pp->wordLength > 16)
			  assMultipleInsert (ass, assVoid(oligo), assVoid(pos + 1)) ;
			else
			  assInsert (ass, assVoid(oligo), assVoid(pos + 1)) ;
		      }
		  }
	      }
	    n1++ ;
	  }
      }
  assDestroy (ass0) ;

  fprintf (stderr, "// Created an associator with %d word in genome of length %d, %s\n"
	   , n1, pp->genomeLength, timeShowNow()) ;
  
  return n1 ;
} /* estScanHashGenome */

/*************************************************************************************/
/****************************** utlities *********************************************/

static BOOL estScanOpenFiles (ESTSCAN *pp)
{
  pp->estFile = filopen (pp->estFileName, 0, "r") ;
  pp->genomeFile = filopen (pp->genomeFileName, 0, "r") ;
  if (!pp->estFile || !pp->genomeFile)
    return FALSE ;
  return TRUE ;
} /* estScanOpenFiles */

/*************************************************************************************/

static BOOL estScanGetGenome (ESTSCAN *pp)
{
  int level, line = 0, nn = 0 ;
  char cc, *cp, *cq ;

  if (!pp->genomeFile)
    return FALSE ;

  pp->genomeStack = stackHandleCreate (10000000, pp->h) ;
  level = freesetfile (pp->genomeFile, 0) ;
  freespecial ("\n") ;

  while (freecard (level)) /* will close pp->estFile */
    {
      line++ ;
      cp = freeword () ;
      if (cp && *cp != '>')
	{
	  cq = cp - 1 ;
	  while (*++cq)
	    {
	      nn++ ;
	      cc = dnaEncodeChar[(int)*cq] ;
	      if (cc)
		*cq = cc ;
	      else
		{
		  fprintf (stderr
			   , "// bad character in file %s line %d, expecting atgc...n, got %s\n"
			   , pp->genomeFileName, line, cq) ;
		  return FALSE ;
		}		     
	    }
	  catText (pp->genomeStack, cp) ;
	}
    }
  pp->genomeFile = 0 ; /* closed by freecard () */
  pp->genomeLength = nn ;
  pp->dnaGenome = (unsigned char *) stackText (pp->genomeStack, 0) ;
  return estScanHashGenome (pp) ;
} /* estScanGetGenome */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr, "// Usage: estScan -e est_fasta_file -g genome_fasta_file -... \n") ;
  fprintf (stderr, "// Example:  estScan -e est.fasta -g chromX.fasta -errMax 2\n") ;  
  fprintf (stderr, 
	   "// -e fileName : the name of a EST target fasta file on which we align the ests\n"
	   "//      this file is read as needed, so there is no limit to its size\n"
	   );
  fprintf (stderr,
	   "// -g : the name of a fasta file containing the genome,\n"
	   "//      (tested with 100 Mb stracches, upper limit depends on harware)\n"
	   ) ;
  fprintf (stderr, 
	   "// -o fileName : outfile name, same as redirecting stdout\n"
	   );

  fprintf (stderr,
	   "// -wordLength wordLength : default 16bp, length of exact match,\n"
	   ) ;
  fprintf (stderr,
	   "// -minEntropy minEntropy : default 12bp, min entropy of hash words,\n"
	   ) ;

  fprintf (stderr,
	   "// -wf fileName [-maxWordRepeat maxWordRepeat] : default 10, max seed occurence in genome,\n"
	   ) ;
  fprintf (stderr,
	   "// -genomeStep genomeStep : default 10bp, how often we hash a genome word,\n"
	   ) ;
  fprintf (stderr,
	   "// -minHitNumber minHitNumber : default 3, minimal number of hits per EST,\n"
	   ) ;


  fprintf (stderr, "// You said: %s\n", commandBuf) ;

      fprintf (stderr, "// ########## ERROR: Unknown argument ") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (stderr, "%s ", argv[i]) ;
      fprintf (stderr, "\n") ;

  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  FILE *f = 0 ;
  char *cp ;
  const char *ccp ;
  const char *outFileName = 0 ;
  int nHits = 0, outlevel = 0, ix ;
  ESTSCAN estScan ;
  AC_HANDLE h = 0 ;
  char commandBuf [1000] ;

  freeinit () ; 
  fprintf (stderr, "// start: %s\n", timeShowNow()) ;
  h = ac_new_handle () ;
  memset (&estScan, 0, sizeof (ESTSCAN)) ;
  estScan.h = h ;

  for (ix = 0, cp = commandBuf ; cp < commandBuf + 900 && ix < argc ; cp += strlen (cp), ix++)
    sprintf(cp, "%s ", argv[ix]) ;
  /* consume optional args */
  getCmdLineOption (&argc, argv, "-o", &outFileName) ;
  getCmdLineOption (&argc, argv, "-e", &(estScan.estFileName)) ;
  getCmdLineOption (&argc, argv, "-g", &(estScan.genomeFileName)) ;
  estScan.wordFrequencyFile = 0 ;
  getCmdLineOption (&argc, argv, "-wf", &(estScan.wordFrequencyFile)) ;

  estScan.wordLength = 16 ;
  if (getCmdLineOption (&argc, argv, "-wordLength", &ccp))
    {
      if (sscanf (ccp, "%d", &ix) == 1 && ix >= 8)   
	estScan.wordLength = ix ;
      else 
	usage (commandBuf, argc, argv) ;
    }
  estScan.maskLength = estScan.wordLength <= 16 ? estScan.wordLength : 16 ;
  ix = estScan.maskLength ; estScan.mask = 0 ;
  while (ix--)
    { estScan.mask <<= 2 ; estScan.mask |= 0x3 ; }
  estScan.minEntropy = 12 ;
  if (getCmdLineOption (&argc, argv, "-minEntropy", &ccp))
    {
      if (sscanf (ccp, "%d", &ix) == 1 && ix >= 5)   
	estScan.minEntropy = ix ;
      else 
	usage (commandBuf, argc, argv) ;
    }
  estScan.genomeStep = 10 ;
  if (getCmdLineOption (&argc, argv, "-genomeStep", &ccp))
    {
      if (sscanf (ccp, "%d", &ix) == 1 && ix >= 1)
	estScan.genomeStep = ix ;
      else 
	usage (commandBuf, argc, argv) ;
    }
  estScan.minHitNumber = 3 ;
  if (getCmdLineOption (&argc, argv, "-minHitNumber", &ccp))
    {
      if (sscanf (ccp, "%d", &ix) == 1 && ix >=0)
	estScan.minHitNumber = ix ;
      else 
	usage (commandBuf, argc, argv) ;
    }

  estScan.maxWordRepeat = 10 ;
  if (getCmdLineOption (&argc, argv, "-maxWordRepeat", &ccp))
    {
      if (sscanf (ccp, "%d", &ix) == 1 && ix >=0)
	estScan.maxWordRepeat = ix ;
      else 
	usage (commandBuf, argc, argv) ;
    }

  if (argc != 1)
    {
      usage (commandBuf, argc, argv) ;
    }
  /* check the absolute args */
  if (!estScan.estFileName)
    {
      fprintf (stderr, "// missing argument -e estFileName\n") ;
      exit (1) ;
    }
  if (! estScan.genomeFileName)
    {
      fprintf (stderr, "// missing argument -g genomeFileName\n") ;
      exit (1) ;
    }
  if (!estScanOpenFiles (&estScan))
    exit (1) ; /* estScanOpen should generate the relevant error message */
  
  if (outFileName)
    {
      f = filopen (outFileName, 0, "w") ;
      if (f)
	outlevel = freeOutSetFile (f) ;
    }
  if (!outlevel)
    outlevel = freeOutSetFile (stdout) ;	

  estScanInit () ;
  if (!estScanGetGenome (&estScan))
    {
      fprintf (stderr, 
	       "// Sorry, I do not understand the genome fasta file,\n"
	       "// please check the name and format\n"
	       "// I expect something like\n"
	       ">contig1\n"
	       "atgtgtccagtagctatattagttc\n"
	       ">contig2\n"
	       "CCTGTGCNNCAT\n"
	       "AATCGCTAAGTGTGCA\n"
	       "// upper lower case can be mixed. the 16 letter dna alphabet is accepted\n"
	       "// We use the standard UPAC coding for representing DNA with a single character\n"
	       "// per base:\n"
	       "// \n"
	       "// Exactly known bases\n"
	       "// \n"
	       "// A\n"
	       "// T  (or U for RNA)\n"
	       "// G\n"
	       "// C\n"
	       "// \n"
	       "// Double ambiguity\n"
	       "// \n"
	       "// R	AG	~Y		puRine\n"
	       "// Y	CT	~R		pYrimidine\n"
	       "// M	AC	~K		aMino\n"
	       "// K	GT	~M		Keto\n"
	       "// S	CG	~S		Strong\n"
	       "// W	AT	~w		Weak\n"
	       "// \n"
	       "// Triple ambiguity\n"
	       "// \n"
	       "// H	AGT	~D		not G	\n"
	       "// B	CGT	~V		not A\n"
	       "// V	ACG	~B		not T\n"
	       "// D	AGT	~H		not C\n"
	       "// \n"
	       "// Total ambiguity\n"
	       "// \n"
	       "// N	ACGT	~N		unkNown\n"
	       "// \n"
	       "// Run without parameter to get help\n"
	       ) ;
      usage (commandBuf, argc, argv) ;
    }

  nHits = estScanRun (&estScan) ;

  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;

  if (estScan.estFile)
    filclose (estScan.estFile) ;
  if (estScan.genomeFile)
    filclose (estScan.genomeFile) ;
  fprintf (stderr, "// done: tested %d est found %d Hits, performed %d verif %s\n"
	   , nEst, nHits, nVerif, timeShowNow()) ;
  fprintf (stderr, "// wordLength = %d, minEntropy = %d, maxWordRepeat = %d, genomeStep = %d, minHitNumber = %d\n"
	   , estScan.wordLength
	   , estScan.minEntropy
	   , estScan.maxWordRepeat
	   , estScan.genomeStep
	   , estScan.minHitNumber
	   ) ;
  if (0) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

