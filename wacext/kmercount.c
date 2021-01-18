/*  File: kmercount.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2014
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

#define LATESORT 0
#include "ac.h"

typedef struct kmerStruct {
  AC_HANDLE h ;
  const char *probeFileName ;
  const char *refFileName ;
  char ccc[256] ;
  BOOL isFastc ;
  BOOL isFastq ;
  BOOL strand ;
  BOOL gzi, gzo ;
  int len, mult ;
  Array aa ;
  DICT *dict ;
  Associator ass ;
  Stack seq ;
  ACEIN ai ; ACEOUT ao ;
} PP ;

/*************************************************************************************/

static long int kmerLongCount (PP *pp, int pass, long int *nrp)
{
  long int nkmers = 0, nr = 0 ;
  int mult = pp->mult ;
  int i, iMax, len = pp->len, n, good ;
  char *cp, *cq, cc, cc1, *cp1 ;
  char buf[pp->len + 1] ;
  Array aa = pp->aa ;
  DICT *dict = pp->dict ;

  memset (buf, 0, sizeof(buf)) ;

  cp = stackText (pp->seq, 0) ;
  vtextLowerCase (cp) ;
 
  iMax = strlen (cp) - len + 1 ;
  cq = cp + len ;
  for (i = good = 0 ; i < iMax ; i++, cp++, cq++)
    {
      cc = *cq ;
      *cq = 0 ;
      cp1 = cp ; cc1 = cp[len/2] ;
      if (! pp->strand && (cc1 == 'a' || cc1 == 'g'))
	{
	  int j ;
	  /* complement the kmer */
	  for (j = 0 ; j < len ; j++)
	    buf[len - j - 1] = pp->ccc[(int)cp[j]] ;
	  cp1 = buf ;
	}
      if (pass == 0)
	{
	  if (dictAdd (dict, cp1, 0))
	    nr += mult ;
	}
      else
	{
	  if (dictFind (dict, cp1, &n))
	    {
	      array (aa, n, int) += mult ;
	      nr += mult ;
	    }
	}
      nkmers += mult ;
      *cq = cc ;
      if (cc == '>')
	{ cp += len + 1 ; cq += len + 1 ; i += len + 1 ; }
    }

  *nrp += nr ;
  return nkmers ;
} /* kmerLongCount */

/*************************************************************************************/

static long int kmerShortCount  (PP *pp, int pass, long int *nrp)
{
  long int nkmers = 0, nr = 0 ;
  int mult = pp->mult ;
  int i, iMax, len = pp->len, lenR, good ;
  int nw = 1 ;
  char *cp, cc ;
  unsigned long int v, vf, vr, w, mask, centralUp ;
  char buf[pp->len + 1] ;
  BOOL strand = pp->strand ;
  Array aa = pp->aa ;
  Associator ass = pp->ass ;

  memset (buf, 0, sizeof(buf)) ;
  cp = stackText (pp->seq, 0) ;
  iMax = strlen (cp) ;

  v = vf = vr = w = mask = 0 ;
  for (i = 0 ; i < len ; i++)
    mask = (mask<<2) | 0x3 ;
  i = len/2 ; /* take the integral part of len/2 */
  lenR = 2 * len - 2 ;
  centralUp = (0x1 << 2 * i) ; /* prefer one strand */

  for (i = good = 0 ; i < iMax ; i++, cp++)
    {
      good++ ;
      switch (*cp)
	{
	case 'a':
	case 'A':
	  cc = 0x0 ;
	  break ;
	case 't':
	case 'T':
	  cc = 0x3 ;
	  break ;
	case 'g':  /* complementary base must differ on their 0x1 bit */
	case 'G':
	  cc = 0x1 ;
	  break ;
	case 'c':
	case 'C':
	  cc = 0x2 ;
	  break ;
	default:
	  good = 0; cc = 0 ;
	  cp += len + 1 ; i += len + 1 ;
	  break ;
	}
      
      if (! strand)
	{
	  /* construct the seed and its complement */
	  vf =  ((vf << 2) & mask) | cc ;
	  vr = ((vr >> 2) | ((((cc&mask) ^ 0x3) & 0x3) << lenR)) & mask ;
	  
	  /* select a the strand for which the central letter == T or C */
	  if ((vf & centralUp) == centralUp) v = vf ;
	  else v = vr ;
	}
      else
	v =  ((v << 2) & mask) | cc ;
      
      if (good < len)
	continue ;
      
      if (pass == 0)
	{
	  if (v && assInsert (ass, (void*) v, assVoid(nw)))
	    {
	      nw++ ; 
	      nr += mult ; 
	    }
	}
      else
	{
	  if (v && assFind (ass, (void*)v, &w))
	    {
	      array (aa, assInt(w), int) += mult ;
	      nr += mult ;
	    }   
	}
      nkmers += mult ;
    }

  *nrp += nr ;
  return nkmers ;
} /* kmerShortCount */

/*************************************************************************************/

static void kmerReport (PP *pp, long int *nrp)
{ 
  Array aa = pp->aa ; 
  int i, *ip, iMax = arrayMax (aa) ;
  int nnn[132] ;
  
  memset (nnn, 0, sizeof(nnn)) ;

  for (i = 1, ip = arrp (aa, 1, int) ; i < iMax ; ip++, i++)
    {
      if (*ip < 32) nnn[*ip]++ ;
      if (*ip >= 32) nnn[32]++ ;
    }
  iMax = *nrp > 1000 ? 32 : 2 ;
  for (i = 0 ; i <= iMax ; i++)
    printf ("%d\t%d\n", i, nnn[i]) ;
} /* kmerReport */

/*************************************************************************************/

static long int kmerParse (PP *pp, int pass, long int *nrp)
{
  long int nkmers = 0 ;
  int line = 0 ;
  int state = 0 ;
  int nseq = 0 ;
  char *cp ;
  char prefix = '>' ;   /* sequence identifier prefix in  fasta and fastc format */
  Stack seq = pp->seq ;
  BOOL isFastc = pp->isFastc ;
  BOOL isFastq = pp->isFastq ;
  ACEIN ai ;
  
  *nrp = 0 ;
  pp->mult = 1 ; 
  ai = aceInCreate (pass ? pp->probeFileName : pp->refFileName, 0, pp->h) ;
  if (!ai)
    messcrash ("cannot read probe file %s", pp->probeFileName) ;
  aceInSpecial (ai, "\n") ;
  
  while (aceInCard (ai))
    {
      line++ ;
      cp = aceInWord (ai) ;
      if (! cp)
	{
	  if (pass && isFastq)
	    {
	      fprintf (stderr, "Sorry, blank line %d in file %s\n", line, aceInFileName (ai)) ;
	      exit (1) ;
	    }
	  else
	    continue ;
	}
      if (pass == 1 && isFastq) /* the reference sequence (pass == 0) must be in fasta format */
	{
	  switch (line % 4)
	    {
	    case 3 : break ; /* quality identifier line, drop it */
	    case 0 : break ; /* quality line, drop it */
	    case 1 : break ; /* read identifier line, drop it */
	    case 2 :  /* actual sequence of the read, work */
	      stackClear (seq) ;
	      pushText (seq, cp) ;
	      if (pp->len > 31)
		nkmers += kmerLongCount (pp, pass, nrp) ;
	      nseq++ ;
	      if (nseq % 1000000 == 0)
		fprintf(stderr, "%s : Processed %d sequences\n", timeShowNow (), nseq) ;
	      break ;
	    }       
	  continue ;
	}
      if (*cp == prefix)
	{
	  if (state == 2) 
	    {
	      if (pp->len > 31)
		nkmers += kmerLongCount (pp, pass, nrp) ;
	      else
		nkmers += kmerShortCount (pp, pass, nrp) ;
	    }
	  
	  state = 2 ; 
	  stackClear (seq) ;
	  pp->mult = isFastc ? fastcMultiplicity (cp, 0, 0) : 1 ;
	  nseq++ ;
	  if (nseq % 1000000 == 0)
	    fprintf(stderr, "%s : Processed %d sequences\n", timeShowNow (), nseq) ;
	  continue ;
	}
      catText (seq, cp) ; /* accumulate the multiline fasta file */
    }
  if (state == 2)
    {
      if (pp->len > 31)
	nkmers += kmerLongCount (pp, pass, nrp) ;
      else
	nkmers += kmerShortCount (pp, pass, nrp) ;
    }

  if (pass)
    kmerReport (pp, nrp) ;


  ac_free (ai) ;
  return nkmers ;
} /* kmerParse */

/*************************************************************************************/
#define DEBUG 1
#ifdef DEBUG4
static long int kmerShortParse (PP *pp, int pass, long int *nrp)
{
  int N = 200000 ;
  long int nkmers = 0, nr = 0 ;
  int nw = 1 ;
  int mult = 0 ;
  int nseq = 0, nnn[128] ;
  int i, iMax, *ip, len = pp->len, lenR, good ;
  char *cp, cc, prefix = '>' ;
#ifdef DEBUG
  char *cq ;
#endif
  BOOL isFastc = pp->isFastc ;
  BOOL isNonStranded = TRUE ;
  ACEIN ai ;
  Array aa = arrayCreate (N, int) ;
  Associator ass = pp->ass ;
  char buf[pp->len + 1] ;
  unsigned long int v, vf, vr, w, mask, centralUp ;

  if (! pass)
    ass = pp->ass = assBigHandleCreate (N, 0) ;

  memset (nnn, 0, sizeof(nnn)) ;
  memset (buf, 0, sizeof(buf)) ;

  v = vf = vr = w = mask = 0 ;
  for (i = 0 ; i < len ; i++)
    mask = (mask<<2) | 0x3 ;
  i = len/2 ; /* take the integral part of len/2 */
  lenR = 2 * len - 2 ;
  centralUp = (0x1 << 2 * i) ; /* prefer one strand */

  *nrp = 0 ;
  ai = aceInCreate (pass ? pp->probeFileName : pp->refFileName, 0, pp->h) ;
  if (!ai)
    messcrash ("cannot read probe file %s", pp->probeFileName) ;
  aceInSpecial (ai, "\n") ;
  
  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (! cp)
	continue ;
      if (cp && *cp == prefix)
	{
	  mult = isFastc ? fastcMultiplicity (cp, 0, 0) : 1 ;
	  nseq++ ;
	  if (nseq % 1000000 == 0)
	    fprintf(stderr, "%s : Processed %d sequences\n", timeShowNow (), nseq) ;
	  continue ;
	}
      iMax = strlen (cp) ;
      for (i = good = 0 ; i < iMax ; i++, cp++)
	{
	  good++ ;
#ifdef DEBUG
	  cq = cp - len ;
#endif
	  switch (*cp)
	    {
	    case 'a':
	    case 'A':
	      cc = 0x0 ;
	      break ;
	    case 't':
	    case 'T':
	      cc = 0x3 ;
	      break ;
	    case 'g':  /* complementary base must differ on their 0x1 bit */
	    case 'G':
	      cc = 0x1 ;
	      break ;
	    case 'c':
	    case 'C':
	      cc = 0x2 ;
	      break ;
	    default:
	      good = 0; cc = 0 ;
	      cp += len + 1 ; i += len + 1 ;
	      break ;
	    }
	  
	  if (isNonStranded)
	    {
	      /* construct the seed and its complement */
	      vf =  ((vf << 2) & mask) | cc ;
	      vr = ((vr >> 2) | ((((cc&mask) ^ 0x3) & 0x3) << lenR)) & mask ;
	      
	      /* select a the strand for which the central letter == T or C */
	      if ((vf & centralUp) == centralUp) v = vf ;
	      else v = vr ;
	    }
	  else
	    v =  ((vf << 2) & mask) | cc ;
	  
	  if (good < len)
	    continue ;


	  if (pass == 0)
	    {
#ifdef DEBUG1
	      fprintf(stderr, "Ref   \t%s\n", cq) ;
#endif
	      if (v && assInsert (ass, (void*) v, assVoid(nw)))
		{
		  nw++ ; 
		  nr += mult ; 
		}
	    }
	  else
	    {
	      if (v && assFind (ass, (void*)v, &w))
		{
		  array (aa, assInt(w), int) += mult ;
		  nr += mult ;	 
#ifdef DEBUG
		  fprintf(stderr,"Found \t%s\n", cq) ;
#endif
		}   
#ifdef DEBUG
	      else
		fprintf(stderr, "Missed\t%s\n", cq) ;
#endif
	    }
	  nkmers += mult ;
	}
    }

  if (pass)
    {
      iMax = arrayMax (aa) ;
      for (i = 1, ip = arrp (aa, 1, int) ; i < iMax ; ip++, i++)
	{
	  if (*ip < 32) nnn[*ip]++ ;
	  if (*ip >= 32) nnn[32]++ ;
	}
      iMax = nr > 1000 ? 32 : 2 ;
      for (i = 0 ; i <= iMax ; i++)
	printf ("%d\t%d\n", i, nnn[i]) ;
    }
  *nrp = nr ;
  return nkmers ;
} /* kmerShortParse */
#endif

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: kmercount -r reference_file -i read_file -k k_mer_length [-strand] [-fastc] [-fastq] \n"
	   "//      try: -h --help  \n"
	   "// Example: kmercount -r RefSeq.fasta.gz -i xyz.fasta.gz -k 31\n"
	   "//\n"
	   "// -r : Reference fasta file, for example the RefSeq or the AceView transcript\n"
	   "//      Sequence identifiers start with '>',\n"
	   "//      The DNA is given on any number of lines of unlimited length\n"
	   "//      k-mers containing letters other than atgc or ATGC are ignored\n"
	   "// -i : The read file, in fasta format\n"
	   "//     -fastc the read file is in fastc format\n"
	   "//     -fastq the read file is in fastq format\n"
	   "// -k : positive integer k giving the length of the k-mer\n"
	   "//      In the [default] non stranded mode, k must be an odd number\n"
	   "//      If k > 31, the code is slower\n"
	   "// -strand : only count the k-mer matching the top strand of the reference file\n"
	   "//           by default, the read and its complement are both tested\n"
	   "//  Read and reference input files called .gz are gzip decompressed on the fly\n"
	   /*
	   "// -o fileName : output file name, equivalent to redirecting stdout\n"
	   "// -gzo : the output file is gziped\n"
	   */
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
  long int nn0 = 0, nn1 = 0, nn2 = 0 ;
  PP p ;
  AC_HANDLE h = 0 ;
  int N = 1000000 ;

  freeinit () ; 

  h = ac_new_handle () ;
  memset (&p, 0, sizeof (PP)) ;
  p.h = h ;

  p.len = 31 ;

  if (argc < 2)
    usage ("Missing arguments", argc, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage ("", 1, argv) ;

  getCmdLineInt (&argc, argv, "-k", &p.len) ;
  p.strand = getCmdLineOption (&argc, argv, "-strand", 0) ;
  if (! p.strand && p.len % 2 == 0)
    {
      fprintf (stderr, "In non stranded protocol, please select an odd k-mer length. Thank you") ;
      exit (1) ;
    }
  
  getCmdLineOption (&argc, argv, "-r", &p.refFileName) ;
  getCmdLineOption (&argc, argv, "-i", &p.probeFileName) ;
  p.isFastc = getCmdLineOption (&argc, argv, "-fastc", 0) ;
  if ( p.probeFileName && strstr ( p.probeFileName, "fastc"))
    p.isFastc = TRUE ;
  if ( p.probeFileName && strstr ( p.probeFileName, "fastq"))
    p.isFastq = TRUE ;

  p.aa = arrayCreate (N, int) ;
  if (p.len < 32)
    p.ass = assBigCreate (N) ;
  else
    p.dict = dictHandleCreate (N, 0) ; ;
  p.seq = stackHandleCreate (100000, 0) ;

  p.ccc['a'] = 't' ;
  p.ccc['t'] = 'a' ;
  p.ccc['g'] = 'c' ;
  p.ccc['c'] = 'g' ;

  fprintf (stderr, "// %s start\n", timeShowNow()) ;
  nn0 = kmerParse (&p, 0, &nn2) ;
  fprintf (stderr, "// %s ref done: found %ld ref %d-mers %ld distinct\n", timeShowNow(), nn0, p.len, nn2) ;

  if (p.probeFileName)
    nn1 = kmerParse (&p, 1, &nn2) ;

  ac_free (p.ao) ;
  ac_free (p.h) ;
  if (1)
    {
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// %s done: found %ld ref %d=mers, %ld read kmers, mapped %ld, %.2f%%  max memory %d Mb\n", timeShowNow(), nn0, p.len, nn1, nn2, 100.0 * nn2/(0.0001+nn1), mx) ;
     }

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
