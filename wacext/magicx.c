/*  File: magicx.c
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

 * A system to align translated RNA-seq againts a protein database

 * 1: Start from the Hexacode (see for example Conway and Sloane book: spere packing...
 *    Find s syem of words with minimal ambiguity under frame shift by rotating the letters

 * 2: Translate a database of protein under this optimal hexacode
 *    Each amino acid is replaced by a canonical 6-mer

 * 3: translate the RNA-seq in 6 frames (or 3 frames if stranded), encode as hexamers

 * 4: align using magicblast

 * 5: filter so that x == 0 modulo 6, a == 0 modulo 6
 *    we expect very few out of frame hits

 * 6: sore the protein hits in the usual way borrowed from tblastx

 *    DONE

 * Code created  2017_10_13

 */

#define VERSION "1.1"

/* #define ARRAY_CHECK  */
#include "ac.h"
#include "channel.h"
#include "ac.h"
 
/* How to download the protein database attention it has 100M sequences inside
 * blastdbcmd -db nr -entry all | head
 */
const char *hexacode[] = {"AAAAAA", "CGCGCG", "TCTCTC", "GTGTGT", "AATTTT", "TTAATT", "TTTTAA", "AAGGGG", "GGAAGG", "GGGGAA", "AACCCC", "CCAACC", "CCCCAA", "TTGGCC", " GGTTCC", "GGCCTT", "TTCCGG", "CCTTGG", "CCGGTT", "CGATAT", "ATCGAT", "ATATCG", "CGTATA", "TACGTA", "TATACG", "CGGCGC", "GCCGGC", "GCGCCG", "ATTAGC", " TAATGC", "TAGCAT", "ATGCTA", "GCATTA", "GCTAAT", "TCAGAG", "AGTCAG", "AGAGTC", "TCGAGA", "GATCGA", "GAGATC", "TCCTCT", "CTTCCT", "CTCTTC", "AGGACT", " GAAGCT", "GACTAG", "AGCTGA", "CTAGGA", "CTGAAG", "GTACAC", "ACGTAC", "ACACGT", "GTCACA", "CAGTCA", "CACAGT", "GTTGTG", "TGGTTG", "TGTGGT", "ACCATG", " CAACTG", "CATGAC", "ACTGCA", "TGACCA", "TGCAAC", 0} ;

/* his set f 21 codons was obtained by running the code with option -bestHexaCode
 * only 20 off frame 6 mers are recognized, we could obtain probably zero
 */

const char *hCodons[] = {
  "GAAAAA"
  , "ACACAT"
  , "CCCCGA"
  , "AAAAGA"
  , "TTTTGA"
  , "AATTCT"
  , "TTAACT"
  , "GCGCTG"
  , "CGCGTG"
  , "ATATTG"
  , "CGATGT"
  , "ATCGGT"
  , "GAGACC"
  , "AGAGCC"
  , "TCTCCC"
  , "AGTCGG"
, "TCAGGG"
  , "GTGTAT"
  , "TGTGAT"
  , "CACAAT"
  , "TGCAGC"
  , "CATGGC"
  , "CTCTCC"
  , "CCGGCT"
  , "AAGGAG"
  , "TTGGTC"
  , "TAACCC"
  , 0
} ;

typedef struct mgxStruct {
  AC_HANDLE h ;
  ACEOUT ao ;
  BOOL gzi, gzo, debug ;
  BOOL p2x ; 
  BOOL r2x ; 
  BOOL wantTitle ;
  BOOL bestHexaCode ;
  BOOL isRaw ;
  int max_threads ;
  const char *inFileName ;
  const char *outFileName ;
} MGX ;

typedef struct fsStruct {
  int rotation ;
  int nFrameShift ;
} FS ;


/*************************************************************************************/
 
static int fsOrder (const void *va, const void *vb)
{
  const FS *up = (const FS *)va, *vp = (const FS*)vb ;
  int n ;

  n = up->nFrameShift - vp->nFrameShift ; if (n) return n ;
  n = up->rotation - vp->rotation ; if (n) return n ;

  return 0 ;
} /* tctCssOrder */

/*************************************************************************************/

static void mgxDnaEncode (char *word)
{
  char *cp = word - 1, cc[256] ;

  memset (cc, 0, sizeof (cc)) ;
  cc['A'] = cc['a'] = 1 ;
  cc['T'] = cc['t'] = 2 ;
  cc['G'] = cc['g'] = 3 ;
  cc['C'] = cc['c'] = 4 ;

  while (*++cp)
      *cp = cc [(int)*cp] ;
  return ;
} /* mgxDnaEncode */

/*************************************************************************************/

static void mgxDnaDecode (char *word, int nn)
{
  char *cp = word - 1, cc[256] ;

  memset (cc, 0, sizeof (cc)) ;
  cc[1] = 'A' ;
  cc[2] = 'T' ;
  cc[3] = 'G' ;
  cc[4] = 'C' ;

  while (cp++, nn--)
      *cp = cc [(int)*cp] ;
  return ;
} /* mgxDnDencode */

/*************************************************************************************/
static void mgxRotate(char *rwords, const char *words, int nn)
{
  int ii, i ;

  memcpy (rwords, words, 64 * 8) ;

  for (i = 0 ; i < 6 ; i++) /* we rotate ony 5 letters */
    {
      char *cp ;
      int k = (nn >> (2 * i)) & 0x3 ;
      if (k)
	for (ii = 0, cp = rwords + i ; ii < 64 ; ii++, cp += 8)
	  *cp = ((*cp) + k ) % 4 + 1 ;
    }
  return ;
} /* mgxRotate */

/*************************************************************************************/

static int mgxCountFrameshifts(char *words)
{
  int nn = 0, ii, jj, kk, dx ;
  char buf[12], *cp, *cq ;

  for (ii = 0; ii < 3 * 27 ; ii+=3)
    {
      memcpy (buf, words + 8 * (ii % 64), 6) ;
      for (jj = 0; jj < 3 * 27 ; jj += 3)
	{
	  memcpy (buf + 6, words + 8 * (jj % 64), 6) ;
	  for (dx = 1 ; dx < 6 ; dx++)
	    {
	      cp = buf + dx ;
	      for (kk = 0 ; kk < 3 * 27 ; kk += 3)
		{
		  cq = words + 8 * (kk % 64) ;
		  if (! memcmp (cp, cq, 6))
		    { nn++ ;  }
		}
	    }
	}
    }

  return nn ;
} /* mgxCountFrameshifts */

/*************************************************************************************/

static void mgxBestHexaCode (MGX *mgx)
{
  char words [64*8] ;
  const char *ccp ;
  int ii ;
  FS *fs ;
  Array aa = arrayHandleCreate (1024,FS, mgx->h) ;

  /* copy the hexacode into a single string
   * each of the 64 word is aligned modulo 8
   * then encode the string
   */
  memset (words, 0, sizeof (words)) ;
  for (ii = 0 ; ii < 64 ; ii++)
    {
      ccp = hexacode[ii] ;
      memcpy (words + 8 * ii, ccp, 6) ;
      mgxDnaEncode (words + 8 * ii) ;
    }

  /* rotate the letters in the  words in all possible way, 
   * fixing the first letter
   * we have 5^4 = 1024 possibilities
   */
  for (ii = 0 ; ii < 1024 ; ii++)
    {
      char rwords [64 * 8] ;
      mgxRotate (rwords, words, ii) ;
      fs = arrayp (aa, ii, FS) ;
      fs->rotation = ii ;
      fs->nFrameShift = mgxCountFrameshifts (rwords) ;
    }

  arraySort (aa, fsOrder) ;
  fprintf (stdout, "# Rotation\tNumber of frameshifed words \n") ;
  for (ii = 0, fs = arrp (aa, 0, FS) ; ii < 24 ; fs++, ii++)
    fprintf (stdout, "%d\t%d\n"
	     , fs->rotation
	     , fs->nFrameShift
	     ) ;

  /* export best */
  if (1)
    {
      int kk ;
      char *cp, rwords [64 * 8] ;
      
      ii = arr (aa, 0, FS).rotation ;
      mgxRotate (rwords, words, ii) ;
     
      mgxDnaDecode (rwords, 8 * 64) ;
      for (kk = 0 ; kk < 3 * 27 ; kk += 3)
	{
	  cp = rwords + 8 * (kk % 64) ;
	  fprintf (stdout, ", \"%s\"\n", cp) ;
	}      
    }

  return ; 
 } /* mgxCleanUp */

/*************************************************************************************/

static void mgxProtein2Hexa (MGX *mgx)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (mgx->inFileName, mgx->gzi, h) ;
  ACEOUT ao = aceOutCreate (mgx->outFileName, ".protein.hexa.fasta", mgx->gzo, h) ;
  char *cp ;
  int ln = 0 ;
  BOOL wantTitle = mgx->wantTitle ;
  BOOL isRaw = mgx->isRaw ;

  while (aceInCard (ai))
    {
      ln++ ;
      cp = aceInWord (ai) ;
      if (! cp || ! *cp)
	continue ;
      if (isRaw)
	{
	  aceOutf (ao, " >%s\n", cp) ;
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  if (! cp)
	    continue ;
	}
      else
	{
	  if (*cp == '#' || *cp == '>') 
	    {
	      aceOutf (ao, "%s", cp) ;
	      if (wantTitle)
		{
		  cp = aceInPos (ai) ;
		  if (cp)
		    aceOutf (ao, " %s", cp) ;
		}
	      aceOutf (ao, "\n") ;
	      continue ;
	    }
	}
      cp-- ;
      while (*++cp)
	{
	  int ii = pepEncodeChar [(int)*cp] ; /* 0->26 */
	  if (ii < 0)
	    messcrash ("Invalid character %s at line %d in file %d\n", ln, aceInFileName (ai)) ;
	  aceOut (ao, hCodons[ii]) ;
	}
      aceOutf (ao, "\n") ;
    }

  ac_free (ai) ;
  ac_free (ao) ;
  ac_free (h) ;
} /* mgxProtein2Hexa */

/*************************************************************************************/

static void mgxRna2Hexa (MGX *mgx)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (mgx->inFileName, mgx->gzi, h) ;
  ACEOUT ao = aceOutCreate (mgx->outFileName, ".rna.hexa.fasta", mgx->gzo, h) ;
  int ln = 0 ;

  while (aceInCard (ai))
    {
      int frame ;
      char *cp0, *cp ;

      ln++ ;
      cp0 = cp = aceInWord (ai) ;
      if (! cp || ! *cp)
	continue ;
      if (*cp == '#' || *cp == '>') 
	{
	  aceOutf (ao, "%s", cp) ;
	  cp = aceInPos (ai) ;
	  if (cp)
	    aceOutf (ao, " ", cp) ;
	  aceOutf (ao, "\n") ;
	  continue ;
	}

      for (frame = 0 ; frame < 3 ; frame++)
	{
	  cp = cp0 + frame ;

	  while (cp[0] && cp[1] && cp[2])
	    {
	      int ii ;
	      char buf[3] ;

	      for (ii = 0 ; ii < 3 ; ii++)
		buf[ii] = dnaEncodeChar[(int)cp[ii]] ;

	      ii = pepEncodeChar [(int)codon(buf)] ; /* 0->26 */
	      if (ii < 0)
		messcrash ("Invalid character %s at line %d in file %d\n", ln, aceInFileName (ai)) ;
	      if (0)
		aceOutf (ao, "%c", pepDecodeChar[ii]) ;
	      else
		aceOut (ao, hCodons[ii]) ;
	      cp += 3 ;
	    }
	  if(0)
	    aceOutf (ao, "\t") ;
	  else
	    aceOutf (ao, "aaaaaaaaaaaaaaaaaa") ;
	}
      aceOutf (ao, "\n") ;
    }

  ac_free (ai) ;
  ac_free (ao) ;
  ac_free (h) ;
} /* mgxRna2Hexa */

/*************************************************************************************/

static void mgxCleanUp (MGX *mgx)
{
  return ; 
} /* mgxCleanUp */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: magicx ..... \n"
	   "//      try: -h --help \n"
	   "// Example\n"
	   "//     magicx ...\n"
	   "// OUTPUT\n"
	   "//   -o fileNamePrefix : output file name, equivalent to redirecting stdout\n"
  	   "//   -gzo : the output files should be gziped\n"
	   "//   -debug : more diagnostics\n"
	   "//   -p2x [-title] : encode a protein fasta file as hexacode\n"
	   "//   -r2x [-title] : encode a rna/dna fasta file as hexacode\n"
	   "//   -raw : input sequence file is in raw format : name tab sequence\n"
	   "//   -bestHexaCode : select best among all 4^5 rotated hexacodes\n"
	   "//   -max_threads <int> : [default 8] number of requested threads\n"
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
  MGX mgx ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&mgx, 0, sizeof (MGX)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  mgx.h = h ;

  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "--help")||
      getCmdLineBool (&argc, argv, "-h")
      )
    usage (commandBuf, 1, argv) ;
  if (getCmdLineBool (&argc, argv, "-v") ||
      getCmdLineBool (&argc, argv, "--version")
      )
    {
      fprintf (stderr, "magicx: %s\n", VERSION) ;
      exit (0) ;
    }
 
  {
    int ix ;
    char *cp ;
    for (ix = 0, cp = commandBuf ;  ix < argc && cp + strlen (argv[ix]) + 1 < commandBuf + 3900 ; cp += strlen (cp), ix++)
      sprintf(cp, "%s ", argv[ix]) ;
  }

  /* consume optional args */
  mgx.wantTitle = getCmdLineBool (&argc, argv, "-title") ;
  mgx.p2x = getCmdLineBool (&argc, argv, "-p2x") ;
  mgx.r2x = getCmdLineBool (&argc, argv, "-r2x") ;
  mgx.isRaw = getCmdLineBool (&argc, argv, "-raw") ;
  mgx.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  mgx.gzo = getCmdLineBool (&argc, argv, "-gzo") ;
  mgx.debug = getCmdLineBool (&argc, argv, "-debug") ;

  getCmdLineOption (&argc, argv, "-i", &(mgx.inFileName)) ;
  getCmdLineOption (&argc, argv, "-o", &(mgx.outFileName)) ;
 
  mgx.bestHexaCode = getCmdLineBool (&argc, argv, "-bestHexaCode") ;
  
  mgx.max_threads = 8 ;
  getCmdLineInt (&argc, argv, "-max_threads", &mgx.max_threads) ;
  
  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }
  /* parallelization */
   if (mgx.max_threads < 4)
     mgx.max_threads = 4 ;
   wego_max_threads (mgx.max_threads) ;
   aceInWaitAndRetry (0) ;
   channelClose (0) ;

  /* check the absolute args */

   if (mgx.bestHexaCode)
     mgxBestHexaCode (&mgx) ;
   if (mgx.p2x)
     mgxProtein2Hexa (&mgx) ;
   if (mgx.r2x)
     mgxRna2Hexa (&mgx) ;


  mgxCleanUp (&mgx) ;
  wego_flush () ; 

  ac_free (mgx.ao) ;
  ac_free (mgx.h) ;
  if (1)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      /*
	fprintf (stderr, "// minSnpCount %d minSnpCover %d minSnpFrequency %d SNP\n"
	       , mgx.minSnpCount
	       , mgx.minSnpCover
	       , mgx.minSnpFrequency
	       ) ;
      */
      fprintf (stderr, "// %s done, max memory %d Mb\n", timeShowNow(), mx) ;
     }

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderrn and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/

