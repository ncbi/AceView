/*  File:  snpbloom.c
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
 *
 * Purpose of this code
 *  Phase 1:
 *   Read a collection of 51 base snippets, hash them using both strands.
 *
 * Phase 2:
 *   Stream a collection of sequences in fasta or fastc format
 *   For each sequence, count the number of exact support of each sequence
 *   so that they can later be collected and treated and merged using bestali.
*/

/* #define ARRAY_CHECK  */

/* 13:5 means take a 13 mer every 6 bases so a 18 mer counts 1, a 19mer counts 2 */
#define WORDLENGTH 15

#define LATESORT 0
#include "ac.h"
#include "bitset.h"
#include "channel.h"
#include "wego.h"

#define IBMAX 10
typedef struct pairStruct { char *buf ; } SHORT_READ ; /* nam:forward seq:reverse seq */
typedef struct seqStruct { 
  int nam ; /* the name of the sequence: defined in the zone dict */
  int ln, x1, x2 ; /* total length of the sequence, the cordicnates of the fragment contained in the zone */
  int pos ; /* offset of the first base of this fragement in the zone stack */
} SEQ ;

typedef enum { FASTA=0, FASTC, FASTQ, CSFASTA, CSFASTC, RAW, GLOBAL } DNAFORMAT ;

typedef struct  pStruct {
  AC_HANDLE h ;
  AC_HANDLE seqH ;
  const char *inFileName, *outFileName, *snpFileName ;
  const char *run ;
  BOOL gzi, gzo, silent, polyA ;
  BOOL telomere ;
  int max_threads, ibMax ;
  DNAFORMAT dnaFormat ;
  CHAN *genomeChan ;
  Associator ass ;
  DICT *snpDict ;
} PP ;

typedef struct telomereStruct { 
  int thId ; 
  int nSeqs, nTags ;
  int ibMax ;
  char *buf ; 
  PP *pp ;
  CHAN *chan, *done ; 
  KEYSET kmerHisto, polyaHisto ; 
  BitSet bb, blooms[IBMAX] ;
} TELO ;

/*************************************************************************************/
/* analyse a sequence
 *   pass 0: creation : fill a Bloom BitSet, return the number words analysed
 *   pass 1 :counting : return the number of hits to the BitSet
 */
typedef struct int4Struct { int xx[IBMAX] ; } I4 ;

#define SHIFT  (2 * (WORDLENGTH - 1))
#define MASK  ((1 << (2 * WORDLENGTH)) - 1)
#define CENTRALMASK (1 << WORDLENGTH)
#define RIGHTMASK (CENTRALMASK - 1)
#define LEFTMASK ((RIGHTMASK << WORDLENGTH) & MASK)

static I4 snpBloomOneSequence (TELO *tp, char *dna)
{
  I4 mynk ;
  register unsigned int f = 0 ; /* the word buffer in the forward direction */
  register unsigned int r = 0 ; /* its reversecomplement r */
  unsigned int g ;     /* the masked compressed selected word */
  register unsigned int x ;  
  register int n = 0 ;
  register const char *cp = dna - 1 ;
  Associator ass = tp->pp->ass ;
  KEYSET kmerHisto = tp->kmerHisto ;
  
  memset (&mynk, 0, sizeof(mynk)) ;

  while (1)
    {
      cp++ ;
      n++ ;
      switch (*cp)
	{
	case 'a': case 'A': x = 0x0 ; break ;
	case 't': case 'T': x = 0x3 ; break ;
	case 'g': case 'G': x = 0x1 ; break ;
	case 'c': case 'C': x = 0x2 ; break ;
	case ' ' : continue ;
	case '\n' : continue ;  /* merge bases across line breaks */
	case 0:  n = 0 ; x = 0x8 ; break ;
        default: x = n = 0 ; break ;
	}
      f = (x | (f << 2)) & MASK ;
      r = (( (x ^ 0x3) << SHIFT) | (r >> 2)) & MASK ;
      g = (f & CENTRALMASK ? r : f) ;                 /* select the strand */
      g = (g & RIGHTMASK) | ((g >> 1) & LEFTMASK) ;   /* gobble the half base used to select the strand */
      
      if (n >= WORDLENGTH)
	{
	  const void *vp ;
	  int snp ;
	  Array bucket = 0 ;
	  int iBucket = 0 ;

	  while (assFindNext (ass, &g, &vp, &bucket, &iBucket))
	    {
	      snp = assInt (vp) ;
	      keySet (kmerHisto, snp)++ ;
	    }
	}
      if (x == 0x8) break ;
    }  
  return mynk ;
} /* snpBloomOneSequence */

/*************************************************************************************/

static void sblDoCountSnps (void *vp)
{
  TELO *tp = (TELO *)vp ;
  int k, nLines = 0, mult = 0, polyA = 0 ;
  char *cp, *cq ;
  BOOL debug = FALSE ;
  DNAFORMAT dnaFormat = tp->pp->dnaFormat ;
  
  while (1)
    {
       if (tp->pp->max_threads)
	 {
	   if (debug) fprintf (stderr, "sblDoCountSnps waits on channelGet (chan-%d) nTags=%d\n", tp->thId, tp->nTags) ;
	   k = channelGive (tp->chan, &k, int) ; /* wait for data to be ready */
	   if (k == -1)    /* work is over */
	     {
	       if (debug) 
		 fprintf (stderr, "sblDoCountSnps channel %d got k == -1 and calls channelPut(done)\n", tp->thId) ;
	       k = - tp->thId - 1 ;
	       channelPut (tp->done, &k, int) ; 
	       return ;
	     }
	 }
      cp = cq = tp->buf ;
      mult = 1 ;
      while (cp && *cp)
	{
	  nLines++ ;
	  switch (dnaFormat)
	    {
	    case GLOBAL:
	      cq = 0 ;  /* accept the whole buffer at once */ 
	      break ; 
	    case FASTQ:
	      cq = strchr (cp, '\n') ;
	      if (cq) *cq++ = 0 ;
	      if ((nLines % 4) != 1) /* only keep the DNA lanes */
		cp = 0 ; 
	      break ;
	    case RAW:
	      cq = strchr (cp, '\n') ;
	      if (cq) *cq++ = 0 ;
	      break ;
	    case FASTC:
	    case FASTA:
	    default:
	      switch ((int)*cp)
		{
		case 0: 
		case '#': 
		  cq = strchr (cp, '\n') ;
		  cp = 0 ; /* jump this line */
		  break ;
		case '>':
		  cq = strchr (cp, '\n') ;
		  if (cq) *cq = 0 ;
		  mult = fastcMultiplicity (cp+1, 0, 0) ;
		  cp = 0 ; /* jump this line */
		  break ;
		default:
		  if (dnaFormat == FASTA) /* merge successive lines */
		    cq = strstr (cp, "\n>") ;
		  else
		    cq = strchr (cp, '\n') ;
		  break ;
		} 
	    }
	  if (cq) { *cq++ = 0 ; }
	  
	  if (cp)
	    {
	      snpBloomOneSequence (tp, cp) ;

	      tp->nSeqs++ ; /* analizing */
	      tp->nTags += mult ; 
	      if (tp->pp->polyA)
		keySet (tp->polyaHisto, polyA) += mult ;
	    }
	  cp = cq ; /* start of next line */
	}
      k = tp->thId ;
      if (tp->pp->max_threads)
	{
	  if (debug) fprintf (stderr, "sblDoCountSnps channel %d and calls channelPut(done)\n", k) ;
	  if (debug) fprintf (stderr, "sblDoCountSnps offset %d , run %s parse %d lines, analyzed %d fragments, %d distinct fragments, %d polyA : %s\n"
			      , tp->thId,  tp->pp->run, nLines, tp->nTags, tp->nSeqs, polyA, timeShowNow()) ;
	  channelPut (tp->done, &k, int) ;  /* this thread is ready for a new block */
	}
      else
	break ;
    }
} /* sblDoCountSnps */

/*************************************************************************************/
/* read a big number of lines, and regularizes to the start of a sequence uing the carryOver */
static int dblFillBuffer (int fd, char *buf, int size, char *carryOverBuf, int carryOver)
{
  int n, n1 = 0, nn = size ;
  char *cp ;
  BOOL debug = FALSE ;

   if (debug) fprintf (stderr, "dblFillBuffercalls read\n") ;
  if (carryOver)
    {
      memcpy (buf, carryOverBuf, carryOver) ;
      nn -= carryOver ;
    } 
   if (debug) fprintf (stderr, "dblFillBuffercalls read\n") ;
  while (nn > 0)   /* read from the file up to size bytes */
    {
      n = read (fd, buf + carryOver + n1, nn) ;
      if (n < 0)
	messcrash ("Error in dblFillBuf:read") ;
      if (debug)  fprintf (stderr, "-----dblFillBuffer requested %d bytes got %d bytes cumul %d bytes\n", nn, n, n1);
      nn -= n ; n1 += n ;
      if (n == 0)
	break ;
    }
   if (debug) fprintf (stderr, "dblFillBuffer got %d bytes\n", n1) ;
  cp = buf + size - nn ;
  cp[0] = 0 ; cp[1] = 0 ;  /* double zero terminate */
  if (nn == 0)
    {
      char *cp0 = cp ;
      while (cp > buf)  /* check for a reminder */
	if (! strncmp (--cp, "\n>", 2))
	  break ;
      if (cp <= buf)   /* this file does not contain a '>' fasta delimiter */
	cp = cp0 ; 
      *cp++ = 0 ;  /* split out the remainder */
      carryOver = strlen (cp) ;
      if (carryOver)
	memcpy (carryOverBuf, cp, carryOver) ;
    }
  else
    {
      cp[0] = 0 ; cp[1] = 0 ; 
      carryOver = 0 ;
    }
  return carryOver ;
}  /* dblFillBuffer */

/*************************************************************************************/
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static void sblParseFastaFile (PP *pp, CHAN *done, int maxTh, int size, Array telos)
{
  char *carryOverBuf = halloc (1 << 21, 0) ; /* 2 Mbytes */
  int k = 0, n, fd, carryOver = 0, nn = 0 ;
  TELO *tt ;
  BOOL debug = FALSE ;

  if (pp->inFileName)
    fd = open (pp->inFileName, O_RDONLY) ;
  else
    fd = 0 ;
  if (fd < 0)
    messcrash ("dblParseFastaFile cannot open file %s", pp->inFileName) ;
  while (1)
    {
      /* listen on the communication to find an available thread */
      if (pp->max_threads)
	{
	  if (nn < maxTh)
	    k = nn++ ;
	  else
	    k = channelGive (done, &k, int) ; /* the number of the channel which is done */
	}
      tt = arrayp (telos, k, TELO) ; 
      if (debug) fprintf (stderr, "dblParseFastaFile calls dblFillBuffer k=%d\n", k) ;
      carryOver = dblFillBuffer (fd, tt->buf, size -  (1 << 21),  carryOverBuf, carryOver) ;
      if (tt->buf[0] == 0)
	break ;  
      n = 1 ;
      if (pp->max_threads) /* signal that data is ready for thread tt */
	{
	  if (debug) fprintf (stderr, "dblParseFastaFile prepared data and calls channelPut(chan-%d)",k) ;
	  channelPut (tt->chan, &n, int) ; /* data ready */
	}
      else  /* direct procedural call */
	sblDoCountSnps (tt) ;
    }

  if (pp->max_threads) /* synchronize */
    {
      for (k = 0 ; k < maxTh ; k++)
	{ 
	  tt = arrayp (telos, k, TELO) ;
	  if (debug) fprintf (stderr, "dblParseFastaFile is over %d tags in chan-%d\n", tt->nTags, k) ;
	  n = -1 ;
	  if (debug) fprintf (stderr, "dblParseFastaFile is over and calls channelPut(chan-%d)",k) ;
	  channelPut (tt->chan, &n, int) ; /* close the thread */
	}
      n = 0 ; k = maxTh ;
      while (k > 0)
	{
	  if (debug) fprintf (stderr, "dblParseFastaFile waits on channelGet(done)") ;
	  n = channelGive (done, &n, int) ;
	  if (n < 0)  /* count the work done signals */
	    {
	      k-- ;    
	      if (debug) fprintf (stderr, "thread %d done %s\n", -n, timeShowNow()) ;
	    }
	}
    }
  
  ac_free (carryOverBuf) ;
  return ;
} /* dblParseFastaFile */

/*************************************************************************************/

static int sblCountSnps (PP *pp)
{
  AC_HANDLE h  = ac_new_handle () ;
  int ibMax = pp->ibMax ;
  int ii, j, nn = 0, maxTh ;
  int nSeqs = 0, nTags = 0 ;
  CHAN *done = channelCreate (24, int, h) ;
  KEYSET kmerHisto, polyaHisto ;
  TELO *tp ;
  Array telos = arrayHandleCreate (pp->max_threads + 1, TELO, h) ;

  int size = (1 << 22) ; /* 4M */ 
  BOOL debug = FALSE ;
  char *title = "xxx" ;

  if (ibMax > IBMAX)
    messcrash ("ibMax = %d > IBMAX = %d", ibMax, IBMAX) ;
  if (debug) channelDebug (done, TRUE, "chan-done") ;
	 			     
  fprintf (stderr, "--- data preparation done : %s\n", timeShowNow ()) ;

  polyaHisto = keySetHandleCreate (h) ;
  kmerHisto = keySetHandleCreate (h) ;

  maxTh = pp->max_threads ;
  if (maxTh < 1) maxTh = 1 ;
  for (ii = 0 ; ii < maxTh ; ii++)
    {  
      /* allocate all elements */
      tp = arrayp (telos, ii, TELO) ;
      tp->kmerHisto = arrayHandleCreate (dictMax (pp->snpDict) * 2, KEY, h) ;
      
      tp->bb =  bitSetCreate (1 << (2 * WORDLENGTH - 1), h) ;
      bitSet (tp->bb, 0) ; /* for compiler happiness */
      tp->buf = halloc (size, h) ;

      if (pp->max_threads)
	{
	  tp->chan = channelCreate (24, int, h) ;
	  if (debug) channelDebug (tp->chan, TRUE, messprintf ("chan-%d", ii)) ;
	}
      /* select the slicing */
      tp->thId = ii ; 
      tp->done = done ;
      tp->pp = pp ;
      tp->ibMax = ibMax ;
      
      /* start the thread */
      if (pp->max_threads)
	wego_go (sblDoCountSnps, tp, TELO) ;
    }
  
  sblParseFastaFile (pp, done, maxTh, size, telos) ;

  /* cumulate the individual histos */
  for (ii = 0 ; ii < maxTh ; ii++)
    { 
      tp = arrayp (telos, ii, TELO) ;

      for (j = keySetMax(tp->kmerHisto) - 1 ; j >= 0 ; j--)
	keySet (kmerHisto, j) +=  keySet(tp->kmerHisto,j) ;
      
      nTags += tp->nTags ;
      nSeqs += tp->nSeqs ;
    }
  fprintf (stderr, "--- histo cumulated %s\n", timeShowNow()) ;
  /* export */
  if (1)
    {
      KEYSET ks ;
      int n, cumulK, cumulpA = 0, jMax, jGlobalMax = 0 ;
      ACEOUT ao = aceOutCreate (pp->outFileName, ".telomere_distribution.txt", pp->gzo, h) ;

      aceOutDate (ao, "##", "Histogram of the number of fragments with n 13-mers matching the canonical telomere motif ") ;
      aceOutf (ao, "## Found %d fragments, %d distinct fragments\n", nTags, nSeqs) ;
      aceOutf (ao, "#Number of motifs\tNumber of fragments with x TERRA 13mers TTAGGGTTAGGGT\tNumber of fragments with x homopolymers of 13 consecutive A\tCumul TERRA\tCumul pA\n") ;

      jGlobalMax = jMax = keySetMax(polyaHisto) ; 
      for (n = j = 0 ; j < jMax ; j++) 
	n += keySet (polyaHisto, j) ;
      cumulpA  = n ;
      
      ks = kmerHisto ;
      jMax = keySetMax(ks) ;
      if (jGlobalMax < jMax) 
	jGlobalMax = jMax ;
      for (n = j = 0 ; j < jMax ; j++) 
	n += keySet (ks, j) ;
      cumulK = n ;

      aceOutf (ao, "\n%s\tA%d", pp->run, WORDLENGTH) ; for (j = 0 ; j < jGlobalMax ; j++) { n =  keySet (polyaHisto, j) ;aceOutf (ao, "\t%d", n) ; }
      aceOutf (ao, "\n%s\tA%d-cumul", pp->run, WORDLENGTH ) ; for (n = j = 0 ; j < jGlobalMax ; j++) { aceOutf (ao, "\t%d",  cumulpA - n) ; n += keySet (polyaHisto, j) ; } 

      ks = kmerHisto ;
      aceOutf (ao, "\n%s\t%s", pp->run, title) ; for (n = j = 0 ; j < jGlobalMax ; j++) { n = keySet (ks, j) ; aceOutf (ao, "\t%d", n) ; }
      aceOutf (ao, "\n%s\t%s-cumul", pp->run, title) ; for (n = j = 0 ; j < jGlobalMax ; j++) { aceOutf (ao, "\t%d", cumulK - n) ; n += keySet (ks, j) ; }
      
      aceOutf (ao, "\n") ;
    }
  ac_free (h) ;	
  return nn ;
} /* sblCountSnps */

/*************************************************************************************/

static BOOL sblParseSnps (PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (pp->snpFileName, FALSE, h) ;
  char *cp ;
  BOOL ok = FALSE ;
  int snp ;

  pp->snpDict = dictHandleCreate (256, pp->h) ;
  pp->ass = assBigHandleCreate (100000, pp->h) ;
  if (! ai)
    goto done ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      cp  = aceInWord (ai) ;
      if (! cp || *cp == '#') continue ;
      dictAdd (pp->snpDict, cp, &snp) ;
      cp  = aceInWord (ai) ;
      if (! cp || *cp == '#') continue ;
      snp <<= 1 ;  /* wild type right bit = 0 */
      assMultipleInsert (pp->ass, cp, assVoid (snp)) ;
      cp  = aceInWord (ai) ;
      if (! cp || *cp == '#') continue ;
      snp++ ; /* mutant right bit = 1 */
      assMultipleInsert (pp->ass, cp, assVoid (snp)) ;
    }

 done:
  ac_free (h) ;
  return ok ;
}  /* sblParseSnps */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: dnabloom  -t target_fasta_file -i probe_fasta_file [-errMax] -... \n"
	   "//      try: -h --help --hits_file_caption \n"
	   "// Example:  clipalign -p tags.fasta -t chromX.fasta -errMax 2\n"
	   "// -t fileName : the name of a target fasta file on which we align the tags\n"
	   "//      this file is read as needed, so there is no limit to its size\n"
	   "//      all named files will be gzip decompressed if they are called *.gz\n"
	   "// -i : the name of a fasta/fastc/fastq file, possibly .gz, containing the tags to be aligned \n"
	   "//      (tested with 10 Million tags, upper limit depends on hardware)\n"
	   "// -fastq33 : the input is in fastq, the quality of the mismatches start at 33=!(NCBI)\n"
	   "// -silent : suppress title lines and status reports from the output, just report the hits\n"
	   "// -gzo : the output file is gziped\n"
	   "// -run runName : the name of the run to be exported in the results\n"
	   "// -telomere : count the polyA and the TERRA motifs unsing exact 13-mers\n"

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
  PP p ;
  AC_HANDLE h = 0 ;
  char commandBuf [1000] ;

  freeinit () ; 
  messErrorInit (argv[0]) ;
  
  h = ac_new_handle () ;
  memset (&p, 0, sizeof (PP)) ;
 
  p.h = h ;
  
  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (commandBuf, 1, argv) ;
  
  /* parse the arguments */
  getCmdLineOption (&argc, argv, "-o", &(p.outFileName)) ;
  getCmdLineOption (&argc, argv, "-i", &(p.inFileName)) ;
  getCmdLineOption (&argc, argv, "-snpFile", &(p.snpFileName)) ;

  if (1)
    {
      p.dnaFormat = FASTC ; /* default */
      const char *ccp = 0 ;
      if (getCmdLineOption (&argc, argv, "-I", &ccp))
	{
	  if (!strcasecmp (ccp, "global"))
	    p.dnaFormat = GLOBAL ;
	  else if (!strcasecmp (ccp, "raw"))
	    p.dnaFormat = RAW ;
	  else if (!strcasecmp (ccp, "FASTC"))
	    p.dnaFormat = FASTC ;
	  else if (!strcasecmp (ccp, "FASTA"))
	    p.dnaFormat = FASTA ;
	  else if (!strcasecmp (ccp, "FASTQ"))
	    p.dnaFormat = FASTQ ;
	  else
	    messcrash ("Sorry, unknown format in option -I %s, try dnabloom -h", ccp) ;
	}
    }

  p.run = "-"  ; /* default */
  getCmdLineOption (&argc, argv, "-run", &(p.run)) ;
  
  p.gzi = getCmdLineOption (&argc, argv, "-gzi", 0);
  p.gzo = getCmdLineOption (&argc, argv, "-gzo", 0);

  p.polyA = getCmdLineOption (&argc, argv, "-polyA", 0);
 
  p.silent = getCmdLineOption (&argc, argv, "-silent", 0);
  p.telomere = getCmdLineOption (&argc, argv, "-telomere", 0);

  p.max_threads = 0 ;
  getCmdLineInt (&argc, argv, "-max_threads", &p.max_threads);
  p.ibMax = 4 ;
  getCmdLineInt (&argc, argv, "-ibMax", &p.ibMax);
 
  if (p.ibMax > IBMAX)
    messcrash ("Please -ibMax %d cannot exceed IBMAX = %d hard defined in the source code", p.ibMax, IBMAX) ;

  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }
  /* check the absolute args */
  if ((WORDLENGTH & 0x1) == 0)
    messcrash ("WORDLENGTH == %d should be odd",  WORDLENGTH) ;
  if (2 *  WORDLENGTH + 1 > 8 * sizeof (int))
    messcrash ("WORDLENGTH = %d should such that, in bits, 2 *  WORDLENGTH + 1 < sizeof(int) = %d"
	       , WORDLENGTH,  8 * sizeof (int)
	       ) ;

  fprintf (stderr, "// %s start\n", timeShowNow()) ;

  if (p.max_threads)
    {
      wego_max_threads (p.max_threads + 2) ;
    }

  if (p.snpFileName)
    {
      sblParseSnps (&p) ;

      sblCountSnps (&p) ;
      goto done ;
    }
 done:
  
  ac_free (p.h) ;
  if (!p.silent)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// %s done: , max memory %d Mb\n", timeShowNow(),  mx) ;
     }

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

