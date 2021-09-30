/*  File: tfali.c
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

 * 2018_04_01
 * strategic idea

 * Construcct K-mer index for genome sections, K-odd, as bitSets

 */

#define VERSION "1.1"

/* #define ARRAY_CHECK  */
#include "ac.h"
#include "bitset.h"

#ifdef RESTRICT 
#define Restrict restrict
#else
#define Restrict 
#endif

/* snp chaining distance */

typedef struct tfa_struct {
  AC_HANDLE h ;
  const char *outFileName ;
  const char *genomeFileName ;
  const char *inFileList ;
  const char *inFileName ;
  BOOL gzi, gzo, debug ;
  BOOL makeIndex ;
  ACEIN ai ; 
  ACEOUT ao ;
  Array chroms, bbb ;
  Array wordCount ;
  DICT *targetDict ;
  DICT *seqDict ;
  const char *run ;
  int wordSize ;
} TFA ;

/*************************************************************************************/
/*************************************************************************************/
#define LN 1000000

static void tfaIndexOneSeg (TFA *tfa, Array wordCount, char *cp0, int ln, BitSet bb)
{
  int i ;
  int w = tfa->wordSize ;
  int iMax = ln - w ;
  long int ww = 0 ;
  long int mask = (1 << 2*w) - 1 ;
  char *cp ;  

  for (i = 0, cp = cp0 ; i < iMax ; i++, cp++)
    {
      char cc = cp[w/2] ; /* central base */

      ww <<= 2 ;
      switch (*cp)
	{
	case 'a':
	case 'A':
	  break ;
	case 't':
	case 'T':
	  ww |= 0x3 ;
	  break ;
	case 'g':
	case 'G':
	  ww |= 0x1 ;
	  break ;
	case 'c':
	case 'C':
	  ww |= 0x2 ;
	  break ;
	}
      ww &= mask ;

      /* select the strand */
      switch (cc)
	{
	case 'a':
	case 'g':
	  break ;
	default:
	  continue ;
	}
      bitSet (bb, ww) ;
      array (wordCount, ww, int)++ ;
    }
} /* tfaIndexOneSeg */

/*************************************************************************************/

static void tfaMakeIndex (TFA *tfa)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, iMax = arrayMax (tfa->chroms) ;
  int w = tfa->wordSize ;
  int kk = 0 ;
  Array wordCount ;

  wordCount = tfa->wordCount = arrayHandleCreate (1<<(2*w), int, tfa->h) ;
  tfa->bbb = arrayHandleCreate (10000, BitSet, tfa->h) ;

  for (ii = 1 ; ii < iMax ; ii++)
    {
      Stack s = array (tfa->chroms, ii, Stack) ;
      char *cp = stackText (s, 0) ;
      int ln = strlen (cp) ;
      int j, jMax = (ln + LN - 1)/LN ;

      for (j = 0 ; j < jMax -1 ; j++)
	{
	  BitSet bb ;
	  char *cq = cp + LN ;
	  char cc = *cq ;
	  *cq = 0 ;
	  
	  bb = bitSetCreate (1<<(2*w), tfa->h) ;
	  array (tfa->bbb, kk++, BitSet) = bb ;
	  tfaIndexOneSeg (tfa, wordCount, cp, cq - cp, bb) ;

	  *cq = cc ;
	}
    }

  arraySort (wordCount, intOrder) ;
  iMax = arrayMax (wordCount) ;
  for (ii = iMax - 1 ; ii >= 0 ; ii -= iMax/100)
    printf ("%d\t%d\n", ii, arr (wordCount, ii, int)) ;

  ac_free (h) ; 
  return ;
} /* tfaMakeIndex*/

/*************************************************************************************/

static void tfaParseGenome (TFA *tfa)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (tfa->genomeFileName, tfa->gzi, h) ;
  Stack s = 0 ;
  int n = 0, kk = 0 ;
  long int nn = 0 ;

  if (! ai)
    messcrash ("Cannot open the genome fasta file %s", tfa->genomeFileName) ;

  tfa->chroms = arrayHandleCreate (256, Stack, h) ;
  while (aceInCard (ai))
    {
      const char *cp = aceInWord (ai) ;

      kk++ ;
      if (!cp || *cp == '#')
	continue ;
      if (*cp == '>')
	{
	  dictAdd (tfa->targetDict, cp+1, &n) ;
	  s = array (tfa->chroms, n, Stack) = stackHandleCreate (1000000,h) ;
	}
      else
	{
	  catText (s, cp) ;
	  nn += strlen (cp) ;
	}
    }
  ac_free (h) ;

  fprintf (stderr, "// Parsed %ld bases in %d chroms %d lines\n", nn, n - 1, kk) ;
  return ;
} /* tfaParseGenome */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: tfali -i file_list ... \n"
	   "//      try: -h --help -v --version\n"
	   "// Examples\n"
	   "//     tfali -g genome.fasta -o myindex --makeIndex -w 13 --run run \n"
	   "//       split the genome as K blocs of size 4^(w-3)\n"
	   "//       For each blck construct a bit set of existing non-starnded w-mers\n"
	   "//       w should be odd, the word is possibly reversed so that the central base is A or G\n"
	   "//   -o prefix\n"
	   "//     all output files will be called prefix.*\n"
	   "//     otherwise they are exported to stdout\n"
  	   "//   -gzo : the output files should be gzipped\n"
	   "// Optional debugging level\n"
	   "//   -debug : more diagnostics\n"
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
  TFA tfa ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&tfa, 0, sizeof (TFA)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  tfa.h = h ;
 
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
      fprintf (stderr, "tfali: %s\n", VERSION) ;
      exit (0) ;
    }
  tfa.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  tfa.gzo = getCmdLineBool (&argc, argv, "-gzo") ;
  tfa.debug = getCmdLineBool (&argc, argv, "-debug") ;
  tfa.makeIndex = getCmdLineBool (&argc, argv, "--makeIndex") ;
  getCmdLineOption (&argc, argv, "-o", &(tfa.outFileName)) ;
  getCmdLineOption (&argc, argv, "-g", &(tfa.genomeFileName)) ;
  getCmdLineOption (&argc, argv, "-i", &(tfa.inFileList)) ;
  getCmdLineOption (&argc, argv, "--run", &(tfa.run)) ;
  getCmdLineInt (&argc, argv, "-w", &tfa.wordSize) ;

  if (tfa.wordSize %2 == 0)
    messcrash ("The -w parameter (%d) should be an odd number, please try --help", tfa.wordSize) ;

  if (tfa.makeIndex)
    {
      if (! tfa.genomeFileName)
	messcrash ("--makeIndex requires -g genomeFileName, please try --help", tfa.wordSize) ;
      tfaParseGenome (&tfa) ;
      tfaMakeIndex (&tfa) ;
    }
  if (1)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// %s done, max memory %d Mb\n"
	       , timeShowNow()
	       , mx) ;
     }

   ac_free (tfa.h) ;

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderrn and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
