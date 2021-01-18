/*  File:  dnawalk.c
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
 * Purpose of this code: 
 *    walk on short reads, to construct long frequent words
 *  Strategy:
 *   for all words of length WORDLENGTH
 *     count their occurence in a direct table indexed by the word
 *     if non stranded, use an odd length and reverse the read if the central base is T or G
 *   sort the words by their number of occurence
 *   for all unflagged words
 *     If its one-base left of right extension has same prevalence n up to +- 60%, chainandnx flag
 *    Report the sufficiently prevalent sufficiently long chains

 *    Prefixes (bar codes) are found by counting all prefixes of length 1 2 3 4 ... 8
 *     then sorting alphabetically and reporting frequent guys

 *   Look at what happens and Bonne Chance
*/

#define HUIT 8

#include "ac.h"
#include "channel.h"
#include "wego.h"

typedef struct seqStruct { 
} SEQ ;

typedef struct  pStruct PP ;
typedef struct ttStruct TP ;
typedef struct chainStruct CHAIN ;
typedef struct wordCountStruct WC ;

struct wordCountStruct {
  int n, n3p ;
  long unsigned int f ;
} ;

struct chainStruct {
  CHAIN *left, *right ; 
  int n, n3p, flag, chainLength ;
  long unsigned int f ;
} ;

struct ttStruct {
  PP *pp ; 
  int id ;
  char *buffer ;
} ;

struct pStruct {
  AC_HANDLE h ;
  Array tps ;
  Associator ass ;
  Array chain ;
  BigArray wordCount ;
  long unsigned nSeq, nTag ;
  long unsigned int *wordNumberp ; /* pp is shared by value, so only pouinted vaalues are actually shared */
  const char *inFileName, *outFileName ;
  const char *run ;
  BOOL gzi, gzo, fastq, silent, useMult, doPrefix, doWalk, doCount, read2 ;
  int max_threads ;
  int L, dL ;
  int leftClipAt, rightClipAt ; /* options to clip the original reads if wished */
  CHAN *inChan, *outChan, *prefixChanIn, *prefixChanOut ;
  int *prefix[HUIT] ;
  int WORDLENGTH, vary ;
  DICT *dict ;
  Array counts ;
  int min_count ;
}  ;

/*************************************************************************************/

static int wordCountOrder (const void *a, const void *b)
{
  WC *up = (WC *)a, *vp = (WC *)b ;
  
  int n = up->n - vp->n ; if(n) return -n ; /* big first */
  return 0 ;
} /* wordCountOrder  */

/*************************************************************************************/

void dwReportCount (void *vp)
{
  AC_HANDLE h = ac_new_handle () ;
  PP *pp = (PP*) vp ; 
  int i ;
  long int k ;
  int min_count = pp->min_count ;
  WC *wc ;
  ACEOUT ao = aceOutCreate (pp->outFileName,".tc", 0, h) ;
  DICT *dict = pp->dict ;
  Array counts = pp->counts ;
  BigArray aa = bigArrayHandleCreate (arrayMax (counts), WC, h) ;
  const char *run = pp->run ? pp->run : "x" ;

  for (i = 0 ; i < arrayMax (counts) ; i++)
    {
      k = i ;
      wc = bigArrayp (aa, k, WC) ;
      wc->n = arr (counts, i, int) ;
      wc->f = k ;
    }
  bigArraySort (aa, wordCountOrder) ;

  aceOutDate (ao, "#", "Walk") ;
aceOutf (ao, "# Run\tRead\tMultiplicity\tTouching 3p end of read\n") ;

  for (k = 0, wc = bigArrp (aa, 0, WC) ; k < bigArrayMax (aa) ; wc++, k++)
    if (wc->n > min_count)
      aceOutf (ao, "%s\t%s\t%ld\t%d\n", run, dictName(dict, (i = wc->f)), wc->n, wc->n3p) ;
 
  ac_free (h) ;
} /* dwReportCount  */

/*************************************************************************************/

void dwReportWalk (void *vp)
{
  AC_HANDLE h = ac_new_handle () ;
  PP *pp = (PP*) vp ; 
  int i, ii, jj,  n0, n03p  ;
  int SHIFT = 2 *  pp->WORDLENGTH - 2 ;
  long unsigned int mask = 1 ;
  long unsigned int k, dk ;
  WC *up ;
  Array chain ;
  CHAIN *xp ;
  ACEOUT ao = aceOutCreate (pp->outFileName,".walk", 0, h) ;
  Associator ass ;
  int nFoundLeft = 0, nFoundRight = 0 ;
  double vary = (100.0 - pp->vary)/100.0 ;
  mask = (mask << (2 * pp->WORDLENGTH)) - 1 ;
  /* kill the associator word -> wordCount, sort the word count and reconstruct
     an associator word ->chain 
  */
  if (1) ac_free (pp->ass) ;
  bigArraySort (pp->wordCount, wordCountOrder) ;

  aceOutDate (ao, "#", "Walk") ;
  aceOutf (ao, "# Run\tRead\tSeq\tTag\tScore\tTailScore\tLength\tWord_length\tSequence\tProfile..\n") ;
  {
    register long unsigned int wordNumber = (*pp->wordNumberp)/10000 ;
    register long int wordCountMax =  bigArrayMax (pp->wordCount) ;

    fprintf (stderr, "chain :  wordCountMax = %lu wordNumber = %lu\n", wordCountMax, wordNumber) ;

    for (ii = 0, up = bigArrp (pp->wordCount, 0, WC) ; ii < wordCountMax ; ii++, up++)
      if (up->n < wordNumber)
	break ;
    jj = ii + 1 ; 
    if (1) ass = pp->ass = assBigHandleCreate (jj, pp->h) ;
    chain = pp->chain = arrayHandleCreate (jj, CHAIN, pp->h) ;
    array (chain, jj - 1, CHAIN).n = 0 ; /* make room */
    for (ii = 0, up = bigArrp (pp->wordCount, 0, WC), xp = arrp (chain, 0, CHAIN) ; ii < jj ; ii++, up++, xp++)
      {
	xp->n = up->n ; xp->n3p = up->n3p ; xp->f = up->f ;
	if (xp->f)
	  assInsert (ass, assVoid (xp->f), assVoid (ii)) ;
      }
    ac_free (pp->wordCount) ;
  }
  fprintf (stderr, "chain->max = %d\n", arrayMax (pp->chain)) ;
  /* Now that we have rorganized the data, 
   * we have a relatively small chain array, pointed at by the ass
   * we now try to chain these elements 
   */
  for (ii = jj = 0 ; ii < arrayMax (pp->chain) && jj < 20 ; ii++)
    {
      CHAIN *up = arrp (pp->chain, ii, CHAIN) ;
      CHAIN *leftMost = up ;
      int chainLength = 1 ;

      if (up->flag)
	continue ;
      n0 = up->n ; /* prevalence of highest word */ 
      n03p = up->n3p ; 
      jj++ ;
      up->flag = 1 ;

      /* search a left extension */
      while (n0 > 0)
	{
	  CHAIN *bestVp = 0 ;
	  k = up->f >> 2 ;
	  for (i = 0 ; i < 4 ; i++)
	    {
	      void *vvp ;

	      dk = i ; dk <<= SHIFT ; 
	      dk |= k ;
	      if (assFind (ass, assVoid(dk), &vvp))
		{
		  int t = assInt (vvp) ;
		  CHAIN *vp = arrp (chain, t, CHAIN) ;

		  if (vp->flag) continue ;
		  if (
		      (vp->n < n0 && 16 * vp->n > 7 * n0) ||
		      (vp->n > n0 && 7 * vp->n < 16 * n0)
		      )
		    {  /* found a left extension */
		       if (! bestVp || vp->n > bestVp->n)
			bestVp = vp ;
		    }
		}
	    }	
	  
	  if (bestVp && bestVp->n)
	    {
	      CHAIN *vp = bestVp ;
	      up->left = vp ;
	      vp->right = up ;
	      vp->flag = 1 ;
	      n0 = vp->n ; n03p = vp->n3p ; up = vp ;
	      leftMost = vp ; 
	      nFoundLeft++ ;
	      chainLength++ ;
	    }
	  else/* no more left extension */
	    break ;
	}

         /* search a right extension */
      up = arrp (pp->chain, ii, CHAIN) ;
      n0 = up->n ;
      n03p = up->n3p ;
      while (n0 > 0 && n03p > 0)
	{
	  CHAIN *bestVp = 0 ;
	  k = ((up->f << 2) & mask) ;
	  for (i = 0 ; i < 4 ; i++)
	    {
	      void *vvp ;

	      dk = i ; dk = k | dk ;
	      if (assFind (ass, assVoid(dk), &vvp))
		{
		  int t = assInt (vvp) ;
		  CHAIN *vp = arrp (chain, t, CHAIN) ;

		  if (vp->flag) continue ;
		  if (
		      (
		       (vp->n < n0 && (100 + vary) * vp->n > n0) ||
		       (vp->n > n0 && vp->n < (100 + vary) * n0)
		       )
		      && ( n03p < 30 ||   /* this is because in micro RNA the first n-mers may not at all touch the end of the read, so the numbers fluctuate a lot */
			   (vp->n3p < n03p && (100 + vary) * vp->n3p > n03p) ||
			   (vp->n3p > n03p && vp->n3p < (100 + vary) * n03p)
		       )
		      )
		    {  /* found a right extension */
		      if (! bestVp || vp->n > bestVp->n)
			bestVp = vp ;
		    }
		}
	    }
	  if (bestVp)
	    {
	      CHAIN *vp = bestVp ;
	      up->right = vp ;
	      vp->left = up ;
	      vp->flag = 1 ;
	      nFoundRight++ ;
	      n0 = vp->n ; n03p = vp->n3p ; up = vp ; chainLength++ ;
	    }
	  else/* no more left extension */
	    break ;
	}
      leftMost->chainLength =  chainLength ;
    }

  /* report the chains */
  /* report */
  fprintf (stderr, "chains, found %d left and %d right extensions\n", nFoundLeft, nFoundRight) ;
  for (ii = jj = 0 ; ii < arrayMax (pp->chain) ; ii++)
    { 
      char *atgc = "agct" ; int dk = pp->WORDLENGTH - 1 ;
      CHAIN *vp, *up = arrp (pp->chain, ii, CHAIN) ;

      if (up->chainLength && up->right && !up->left)
	{ 
	  int k = 1, kk ;
	  char buf[2500] ;
	  long unsigned int score = up->n, score3p = up->n3p, tailScore = 0, tailScore3p = 0 ;
	  
	  for (k = 0 ; k < pp->WORDLENGTH ; k++)
	    buf[k] = atgc[(up->f >> (2*(dk - k))) & 0x3] ;
	  for (vp = up->right, kk = up->chainLength - 2 ; vp ; vp = vp->right, kk--)
	    {
	      buf[k++] = atgc[vp->f & 0x3] ;
	      score += vp->n ; score3p += vp->n3p ;
	      if (kk < dk)
		{ tailScore  += vp->n ; tailScore3p += vp->n3p ; }
		
	    }
	  buf[k++] = 0 ;
	  aceOutf (ao, "%s\t%d\t%lu\t%lu\t%lu:%lu\t%lu:%lu\t%d\t%d\t%s"
		   , pp->run, pp->read2 ? 2 : 1
		   , pp->nSeq, pp->nTag
		   , score, score3p, tailScore, tailScore3p
		   , pp->WORDLENGTH + up->chainLength - 1, pp->WORDLENGTH
		   , buf
		   ) ;
	  for (vp = up ; vp ; vp = vp->right)
	    aceOutf (ao, " %d/%d", vp->n, vp->n3p) ;
	  aceOutf (ao, "\n") ;
	}
    }
  ac_free (h) ;
} /* dwReportWalk  */

/*************************************************************************************/

void dwReportPrefix (void *vp)
{
  AC_HANDLE h = ac_new_handle () ;
  PP *pp = (PP*) vp ; 
  int **prefix = pp->prefix ;
  int j, jMax, jj, *jp, nT = 0, nT1 ;
  char *atgc = "agct" ;
  /*   ACEOUT ao = aceOutCreate (pp->outFileName,".prefix", 0, h) ; */
  jp = prefix[0] ;  /* single letters */
  for (j = 0 ; j < 4 ; j++) 
    nT += *jp++ ;
  printf (".\t%d\n", nT) ;
  nT1 = nT / 8 ;
  for (jj = 0 ; jj < HUIT ; jj++)
    { 
      nT1 = nT / (jj == 0 ? 2 : 8) ;
      jp = prefix[jj] ;
      jMax = 1 << (2 * jj + 2) ;
      for (j = 0 ; j < jMax ; j++)
	if (jp[j] > nT1)
	  {
	    char buf[24] ;
	    int k ;
	    for (k = 0 ; k <= jj ; k++)
	      buf[k] = atgc[(j >> (2*(jj - k))) & 0x3] ;
	    buf[k++] = '\t' ;
	    sprintf (buf + k, "%d", jp[j]) ;
	    
	    printf ("%s\n", buf) ;
	  }
    }
  ac_free (h) ;
} /* dwReportPrefix */

/*************************************************************************************/
/* this code communicates via channels */

void dwPrefix (void *vp)
{
  TP *tp = 0 ;
  PP *pp = (PP*) vp ;
  BigArray wordCount = pp->wordCount ;
  Associator ass = pp->ass ;
  char *cp ;
  int SHIFT = 8 * sizeof (long unsigned int) - 2 ;
  int jj, n, mult = 1 ;
  register long unsigned int f = 0 ; /* the word buffer in the forward direction */
  register long unsigned int r = 0 ; /* its reversecomplement r */
  register long unsigned int mask = 1 ;
  register long unsigned int x ;
  int WORDLENGTH = pp->WORDLENGTH ;
  int **prefix = pp->prefix ;
  Array counts = pp->counts ;
  DICT *dict = pp->dict ;

  mask = (mask << (2 * WORDLENGTH)) - 1 ;
  while (channelMultiGet (pp->prefixChanIn, &jj, 1, int))
    {
      long unsigned int wordNumber = 0 ;

      tp = arrayp (pp->tps, jj, TP) ;
      n = 0 ;
      cp = tp->buffer ;
      while (*cp)
	{
	  n++ ;
	  switch (*cp++)
	    {
	    case '0': case 'a': case 'A': x = 0x0 ; break ;
	    case '1': case 't': case 'T': x = 0x3 ; break ;
	    case '2': case 'g': case 'G': x = 0x1 ; break ;
	    case '3': case 'c': case 'C': x = 0x2 ; break ;
	    case ' ' : continue ;
	    case '\n' : continue ;  /* merge bases across line breaks */
	    case '#': 
	      mult = 0 ; memcpy (&mult, cp, 4) ; cp += 4 ;  x = n = f = r = 0 ; 
	      if (dict && counts && cp && *cp)
		{
		  int k ;
		  char *cr = strchr (cp, '#') ;
		  if (cr) *cr = 0 ;
		  dictAdd (dict, cp, &k) ;
		  array (counts, k, int) += mult ;
		  if (cr) *cr = '#' ;
		}	      
	      break ;
	    default: x = n = f = r = 0 ; break ;
	    }
	  f = (x | (f << 2)) & mask ;
	  r = (( (x ^ 0x3) << SHIFT) | (r >> 2)) & mask ;
	  
	  if (n == HUIT)
	    {
	      register long unsigned int k ;
	      int *jp ;
	      for (jj = 0 ; jj < HUIT ; jj++)
		{
		  k = f >> (2 * (HUIT - jj -1 )) ;
		  jp = prefix[jj] + k ;
		  (*jp)++ ;
		}	  
	    }
	  if (n < WORDLENGTH)
	    continue ;
	  if (wordCount && f)
	    {
	      long int t ;
	      void *vvp = 0 ;
	      WC *up ;

	      if (! assFind (ass, assVoid (f), &vvp))
		{
		  t = bigArrayMax (wordCount) ;
		  assInsert (ass, assVoid (f), assVoid (t)) ;
		}
	      else
		t = assULong (vvp) ;
	      up = bigArrayp (wordCount, t, WC) ;
	      up->n += mult ; up->f = f ;
	      if (*cp == '#' || *cp == 0) up->n3p += mult ;
	      wordNumber += mult ;
	    }
	} 
      if (pp->wordNumberp) 
	*(pp->wordNumberp) += wordNumber ;
      ac_free (tp->buffer) ;
      channelPut (pp->inChan, &(tp->id), int) ; /* this tp is again available */
    }
  channelPut (pp->prefixChanOut, &(tp->id), int) ;
}  /* dwPrefix */

/*************************************************************************************/

void dwParse (void *vp)
{
  AC_HANDLE h = ac_new_handle () ;
  TP *tp ;
  PP *pp = (PP *)vp ;
  ACEIN ai = 0 ;
  int ii, state = 0, mult = 1, n ; 
  int nMax = 1000000 ;
  int nFree = .9 * nMax ;
  char *bb ;
  int nx = HUIT ;
  int nLine = 0 ;

  if (pp->doCount || pp->doWalk) nx = 8888888 ;
  ii = channelGive (pp->inChan, &ii, int) ; /* will block until things a re ready */
  ai = aceInCreate (pp->inFileName, pp->gzi, h) ;
  tp = arrayp (pp->tps, ii, TP) ;
  tp->buffer = halloc (nMax, 0) ;
  bb = tp->buffer ;
  nFree = .9 * nMax ;

  while (aceInCard (ai))
    {
      char *cq, *cp = aceInWord (ai) ;
      if (!cp)
	continue ;
      nLine++ ;
      if (pp->fastq)
	{
	  if (nLine % 4 != 1) 
	    continue ;
	  else
	    { state = 2 ; mult = 1 ; }
	}
      if (1)
	{
	  if (state == 2 || *cp == '>')
	    {
	      if (10*nFree < nMax)
		{
		  int i4 = 0 ;
		  TP *tp4 = 0 ;
		  i4 = channelGive (pp->inChan, &i4, int) ; 
		  tp4 = arrayp (pp->tps, i4, TP) ;
		  tp4->buffer = tp->buffer ;
		  tp4->pp = pp ;
		  channelPut (pp->prefixChanIn, &(tp4->id), int) ;
		  /* will activate 	      dwPrefix (tp4) ; */
		  bb = tp->buffer = halloc (nMax, 0) ; bb[0] = 0 ;
		  nFree =.9 *  nMax ;
		}
	      
	      *bb++ = '#' ; nFree-- ;
	      state = 1 ;
	      mult = 1 ;
	      if (pp->useMult && (cq = strchr (cp, '#')))
		{
		  cq++ ; mult = 0 ;
		  while (*cq && *cq <= '9' && *cq >= '0')
		    mult = 10 * mult + (*cq++ - '0') ;
		}
	      memcpy (bb, &mult, 4) ; bb += 4 ; nFree -= 4 ;
	      pp->nSeq++ ; pp->nTag += mult ;
	      if (! pp->fastq)
		continue ;
	    }
	}
      switch (state)
	{
	case 0: /* out */
	  continue ;
	  break ;
	case 1:
	  cq = strchr (cp, '>') ;
	  if (cq) { *cq = 0 ; if (pp->read2) { cp = cq + 1 ; if (cp[0] == '<') cp++; }}
	  n = strlen (cp) ;
	  if (n > nFree)
	    { *bb++ = 0 ; state = 0 ; nFree = 0 ; continue ; }	
	  if (n > nx)
	    { cp[nx] = 0 ; n = nx ; }
	  memcpy (bb, cp, n) ; bb += n ; nFree -= n ;
	}
    }    
  ii = channelGive (pp->inChan, &ii, int) ; 
  channelPut (pp->prefixChanIn, &(tp->id), int) ;
  channelClose (pp->prefixChanIn) ;
  channelMultiGet (pp->prefixChanOut, &ii, 1, int) ;  /* dwPrefix is closed */
  

  if (pp->doPrefix) dwReportPrefix (pp) ;
  if (pp->doWalk && *(pp->wordNumberp)) dwReportWalk (pp) ;
  if (pp->doCount) dwReportCount (pp) ;

  /* done: transmit on channel that we are done */
  channelPut (pp->outChan, &(tp->id), int) ;
  fprintf (stderr, "Parsed %d lines from file %s\n", nLine, aceInFileName (ai)) ;
  return ;
} /* dwParse */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: dnawalk -i <fasta_file> -... \n"
	   "//      try: -h --help --hits_file_caption \n"
	   "// Example:  dnawalk -i f.fasta.gz -walk -wordLength 11 -prefix\n"
	   "// -i : the name of a fasta/fastc/fastq file, possibly .gz, containing the reads \n"
	   "// -fastq : the reads are in fastq format\n"
	   "// -gzo : the output file is gziped\n"
	   "// -run runName : the name of the run to be exported in the results\n"
	   "// -mult : use the read multiplicities, default: FALSE : use fastc merge identical sequences\n"
	   "// -prefix : search for up word 1 to 8 bases highly prevalent at 5' end, probably a barcode\n"
	   "// -walk -wordLength <int> -vary <int>: search for the most frequent subsequences\n"
	   "//     Start from a frequent word of given length (default 25) and extend it one base at a time\n"
	   "//     with the constraint that the support of the successive words vary by less than p percent\n"
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
  int jj ;
  char commandBuf [1000] ;
  long unsigned int wordNumber = 0 ;
  BOOL debug = FALSE ;

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

  p.useMult = getCmdLineBool (&argc, argv, "-mult") ;
  p.doPrefix = getCmdLineBool (&argc, argv, "-prefix") ;
  p.doWalk = getCmdLineBool (&argc, argv, "-walk") ;
  p.doCount = getCmdLineBool (&argc, argv, "-count") ;
  p.run = "-"  ; /* default */
  getCmdLineOption (&argc, argv, "-run", &(p.run)) ;
  
  p.gzi = getCmdLineOption (&argc, argv, "-gzi", 0);
  p.gzo = getCmdLineOption (&argc, argv, "-gzo", 0);

  p.silent = getCmdLineOption (&argc, argv, "-silent", 0);
  p.fastq = getCmdLineOption (&argc, argv, "-fastq", 0);

  p.read2 = getCmdLineOption (&argc, argv, "-read2", 0);

  p.max_threads = 0 ;
  getCmdLineInt (&argc, argv, "-max_threads", &p.max_threads);
  p.min_count = 30 ;
  getCmdLineInt (&argc, argv, "-min_count", &p.min_count);
  p.WORDLENGTH = 25 ; p.vary = 30 ;
  getCmdLineInt (&argc, argv, "-wordLength", &p.WORDLENGTH) ;
  getCmdLineInt (&argc, argv, "-vary", &p.vary) ;

  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }
    
  fprintf (stderr, "// %s start\n", timeShowNow()) ;
  if (1)
    {
      /* initialize the prefix table */
      int jj, n = 1 ;
      for (jj = 0 ; jj < HUIT ; jj++)
	{
	  n <<= 2 ;
	  p.prefix[jj] = halloc (n * sizeof (int), h) ;
	}
    }
  if (p.doCount)
    {
      p.dict = dictHandleCreate (1000000, h) ;
      p.counts = arrayHandleCreate (1000000, int, h) ;
    }
  if (p.doWalk)
    {
      p.wordNumberp = &wordNumber ; 
      p.wordCount = bigArrayHandleCreate (100000, WC, p.h) ;
      bigArray (p.wordCount, 0, WC).n = 0 ; /* avoid zero */
      p.ass = assBigHandleCreate (100000, p.h) ;
    }

  wego_max_threads (p.max_threads + 4) ;

  p.inChan = channelCreate (2, int, p.h) ;
  p.outChan = channelCreate (2, int, p.h) ;
  p.prefixChanIn = channelCreate (2, int, p.h) ;
  p.prefixChanOut = channelCreate (2, int, p.h) ;
  p.tps = arrayHandleCreate (4, TP, h) ;
  channelDebug (p.inChan, debug, "inChan") ;
  channelDebug (p.outChan, debug, "outChan") ;
  channelDebug (p.prefixChanIn, debug, "prefixChanIn") ;
  channelDebug (p.prefixChanOut, debug, "prefixChanOut") ;

  wego_go (dwParse, &p, PP) ;
  wego_go (dwPrefix, &p, PP) ;
  

  if (1)
    {
      /* initialize the prefix table */
      int kk ;
      for (kk = 0 ; kk < 2 ; kk++)
	{
	  TP *tp = arrayp (p.tps, kk, TP) ; 
	  tp->pp = &p ;
	  tp->id = kk ;
	  
	  channelPut (p.inChan, &kk, int) ; /* start jumpAlignImportReadsDict */
	}
      }

  channelMultiGet (p.outChan, &jj, 1, int) ;  /* wait on readDict initialisation */
  
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

#ifdef JUNK	  

# Pre-analysis of the Rat miRNAs using the Magic pipeline
#    Danielle and Jean Thierry-Mieg, Shanghai, feb 17m 2016

# The fastq files for 320 miRNA Illuniman runs were provided by 
# the Leming Shi laboratory in the fastq directory
ls Fastq | sed -e 's/\.txt\.gz//' > myList

# The adaptors were identified by the command
mkdir Walk
foreach run (`cat myList`)
   if (! -e Walk/$run.err) then
     echo -n "$run "
     date
     scripts/submit Walk/$run "bin/dnawalk -prefix -walk -wordLength 13 -fastq -i Fastq/$run.txt.gz -o Walk/$run"
   endif
end
# No 5p adaptor was identified
# The 3p adaptors matched the known Illumina adaptor with variable NNNNNN bar code given in capital on the last line
echo ' ' > adaptors
foreach run (`cat myList`)
  if (! -e Walk/$run.walk) continue
  echo -n "$run\t" >> adaptors
  grep tggaattct Walk/$run.walk | head -1 | gawk '{print $5}' >> adaptors
end
tggaattctcgggtgccaaggaactccagtcacgacgacatctc
tggaattctcgggtgccaaggaactccagtcacgtagagatctc 
TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

# Reads where then clipped where their tails soft-matched the adaptor by the command
 mkdir TC
 foreach run (`cat myList`)
   if (! -e TC/$run.tc) then
     echo -n "$run "
     date
     scripts/submit TC/$run "bin/dna2dna -I fastq -i Fastq/$run.txt.gz -O tc -o TC/$run -rightClipOn TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
   endif
  end
# The complete collection of reads was then collated, sorted by decreasing number of occurences, adding the cumulated counts
cat myList | gawk -F '_' '{print $2}' | sort -u > TissueList
foreach tissue (`cat TissueList`)
  echo -n "$tissue "
  date
  cat TC/ *$tissue*.tc | bin/dna2dna -I tc -O tc | sort -k 2nr |  gawk -F '\t' '{n += $2; printf("%s\t%d\t%d\n",$1,$2,n);}' > allReads.$tissue.tc
end
cat allReads.*.tc | gawk -F '\t' '{if($2>=10)print}' | bin/dna2dna -I tc -O tc | sort -k 2nr |  gawk -F '\t' '{n += $2; printf("%s\t%d\t%d\n",$1,$2,n);}' > allReads.tc
# Reads between 17 and 26 bp occuring at least 10000 times were selected
cat allReads.tc | gawk -F '\t' '{n=length($1);if(n >= 17 && n <= 27 && $2 >= 10000) print;}' | dna2dna -I tc -O fastc -o bestMiRNA

# Reads belonging to this collection were then selected in each run, reformated as pseudo gene-counts, and used in further analysis
cat bestMiRNA.fastc | gawk '/^>/{next;}{print}' > bestMiRNA.list

mkdir Fastc
foreach run (`cat myList`)
   set run2=`echo $run | gawk '{split($1,aa,"_");printf ("%s_%s_%s_%s",aa[2],aa[3],aa[4],aa[5]);}'`
   echo $run2
   if (! -e Fastc/$run2/f.1.fastc.gz && -e TC/$run.tc) then
     mkdir Fastc/$run2
     cat bestMiRNA.list ZZZZZ TC/$run.tc | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1]==1)print;}' | bin/dna2dna -I tc -O fastc -count -gzo -o Fastc/$run2/f.1
     echo $run2/f.1 > Fastc/$run2/LaneList
   endif
 end

ls Fastc > MetaDB/$MAGIC/RunList
cat bestMiRNA.list | gawk '{printf ("Gene %s\nTargeted\n\n",$1);}' > targetedGenes.ace

# construct pseudo gene count
mkdir tmp tmp/GENERUNS

\rm tmp/GENERUNS/Rat_miRNA.miRNA.GENE.u.ace
foreach run (`cat MetaDB/$MAGIC/RunList `)
  zcat Fastc/$run/f.1.fastc.gz | dna2dna -I fastc -O tc | gawk -F '\t' '{if ($2 >= 0) printf ("Gene %s\nRun_U %s 0 %d seq %d tags %d base %d reads %d compRead\n\n",$1,run,$2,$2,2*$2,$2,$2);}' run=$run >> tmp/GENERUNS/Rat_miRNA.miRNA.GENE.u.ace
end



#######################################
###### Construct the groups
cat myList | gawk '{split($1,aa,"_");printf ("Run %s_%s_%s_%s\nRunId %s\nGroup %s\nGroup %s\nGroup %s_weeks\nGroup Rat_%s\nProject Rat_miRNA\nRNA\n\n",aa[2],aa[3],aa[4],aa[5],$1,aa[2],aa[3],aa[4],aa[5]) ;}' > runs.ace

cat myList | gawk '{split($1,aa,"_");printf ("mkdir Fastc/%s_%s_%s_%s\nmv Fastc/%s/ * Fastc/%s_%s_%s_%s\n",aa[2],aa[3],aa[4],aa[5],$1,aa[2],aa[3],aa[4],aa[5]) ;}' > _m


cat Counts/ *.tc | gawk '{l=length($2);if(l>=21 && l<=23)n[$2]+=$3;}END{for(k in n)printf("%s\t%d\n",k,n[k]);}' | sort -k 2nr | gawk '{k+=$2;u++;printf("%d\t%s\t%s\t%d\n",u,$1,$2,k);}' | gawk '{if($1 % 1 == 0 && $3 > 1000)print}' > RESULTS/prevalence.txt
#######################################



scripts/geneindex.tcsh snp4 000
cat targetedGenes.ace >> tmp/GENERUNS/Rat_miRNA.miRNA.GENE.info.ace

cat MetaDB/RunList | gawk '{split($1,aa,"_");printf ("Run %s\nGroup Age_%s_%s\n\n",$1,aa[1],aa[3]);printf("Compare Age_%s\nRuns Age_%s_%s\n\n",aa[1],aa[1],aa[3])}' > toto.ace



#endif
