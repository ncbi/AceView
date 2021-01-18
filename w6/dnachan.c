#include <ac.h>
#include "channel.h"
#include "wego.h"

typedef enum { FASTA=0, FASTC, FASTQ, CSFASTA, CSFASTC, CCFA, CCFAR, RAW, CRAW} DNAFORMAT ;
typedef struct ppStruct {
  BOOL gzi ;
  BOOL gzo ;
  BOOL silent ;
  const char *inFileName ;
  const char *outFileName ;
  const char *suffix ;
  CHAN *ventChans[64] ;
  CHAN *doneChan ;
  int nPrefix ;
  int nVentilators ;
  int splitMb, maxBytes ; /* maxSize of eac output file, unless one sequence is bigger */
  int nWriters, iWriter, nWritersMask ;
  int blockSize ; /* in Mbytes to fractionate the input */
  DNAFORMAT dnaFormat ;
} PP ;

typedef struct vvStruct {
  int ln ;   /* current length */
  int nSeq ; /* current number of sequences */
  int prefix ; /* 1 to 256 */
  int writer ; /* number of writer agent, typically 4 or 16 */
  char *buf ;
  int max_threads ;
} VV ;


/************************************************************/
/************************************************************/
/************************************************************/

/* read a big number of lines, and regularizes to the start of a sequence uing the carryOver */
static int dnaChanFillBuffer (int fd, char *buf, int size, char *carryOverBuf, int carryOver)
{
  int n, n1 = 0, nn = size ;
  char *cp ;
  BOOL debug = FALSE ;

   if (debug) fprintf (stderr, "dnaChanBuffercalls read\n") ;
  if (carryOver)
    {
      memcpy (buf, carryOverBuf, carryOver) ;
      nn -= carryOver ;
    } 
   if (debug) fprintf (stderr, "dnaChanBuffercalls read\n") ;
  while (nn > 0)   /* read from the file up to size bytes */
    {
      n = read (fd, buf + carryOver + n1, nn) ;
      if (n < 0)
	messcrash ("Error in dnaChanBuf:read") ;
      if (debug)  fprintf (stderr, "-----dnaChanBuffer requested %d bytes got %d bytes cumul %d bytes\n", nn, n, n1);
      nn -= n ; n1 += n ;
      if (n == 0)
	break ;
    }
   if (debug) fprintf (stderr, "dnaChanBuffer got %d bytes\n", n1) ;
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
}  /* dnaChanFillBuffer */

/*************************************************************************************/

static void dnaChanFastaWriter (void *vp)
{
  PP *pp = (PP *)vp ;
  AC_HANDLE h = ac_new_handle() ;
  CHAN *ventChan = pp->ventChan [pp->iWriter] ;
  int kk = pp->nVent ;
  int nn [256] ;
  VV *vv ;

  memset (nn, sizeof(nn], 0) ;
  memset (ln, sizeof(ln], 0) ;

  while (kk && channelGet (ventChan, &vv, VV))
    {
      if (!vv || ! vv->ln)
	{ kk-- ; continue ; }
      ao = arrp (aos, vv->prefix, ACEOUT) ;
      if (ao && ln[vv->prefix] + vv->ln > pp->maxBytes)
	ac_free (ao) ; 	  
      if (!a0)
	{
	  ln[prefix] = 0 ;
	  ao = arrp (aos, vv->prefix, ACEOUT) =
	    aceOutCreate (sprintf (aoBuf, "%s.%d.%d.%s"
				   , pp->outFileName
				   , vv->prefix
				   , 1+nn[vv->prefix]++
				   , pp->suffix
				   )
			  0, 0, h) ;
	  aceOutPrint (ao, vv->buf) ;
	  ln[prefix] += vv->ln ; 
	  free (vv->buf) ;
	}
    }

  channelPut (pp->doneChan, pp->iWriter, int) ;
  ac_free (h) ;
  return ;
}

/*************************************************************************************/

static void dnaChanVentilate (void *vp)
{
  PP *pp = (PP *)vp ;
  CHAN *ventChan = pp->ventChans[pp->iWriter] ;
  int kk = pp->nVentilators ;
  int i ;
  int bufSize = 1 << 20 ;
  nPrefix = pp->nPrefix ;
  VV *vv, *vvv[nPrefix]

    /* vv is a buffer that will only hold reads with a given prefix */
  for (i = 0 ; i < nPrefix ; i++)
    {
      vv = vvv[i] = (VV*) malloc (sizeof(VV), 0) ;
      vv->prefix = i ;
      vv->writer = (i >> (8 - pp->nWritersMask)) ;
      vv->ln = 0 ;
      vv->buf = malloc (bufSize) ;
    }
  /* buf is a big buffer containing fasta or fatc reads */
  while (channelGet (bufChan, &buf, VV))
    {
      if (!vent)
	{ kk-- ; continue ; }
      cp = buf ;
      while (cp[0])
	{
	  cq = strchr (cp, '\n') ; /* name ends */
	  if (!cq)
	    break ;
	  nam = cp ;  ln1 = cq - cp + 1 ; cp = cq + 1 ;
	  cq = strchr (cp, '\n') ; /* dna ends */
	  if (!cq)
	    break ;
	  dna = cp ; ln2 = cq - cp + 1 ;  cp = cq + 1 ;
	  memcpy (aa->base, dna, 4) ;
	  dnaEncodeArray (aa) ;
	  cq = arrp (aa, 0, char*) ;
	  /* read the 4 first letters */
	  x = 
	    ((cq[0] & 0x3) << 6) + 
	    ((cq[1] & 0x3) << 4) + 
	    ((cq[2] & 0x3) << 2) + 
	    (cq[3] & 0x3) ;
	  ln = ln1 + ln2 + 1 ;
	  vv = vvv[x] ;
	  if (vv->ln + ln + 1 > bufSize)
	    {
	      if (vv->ln)
		{
		  channelPut (pp->ventChans[vv->writer], vv, VV*) ;
		  vv = vvv[x] = (VV*) malloc (sizeof(VV), 0) ;
		  vv->prefix = x ;
		  vv->writer = (x >> (8 - pp->nWritersMask)) ;
		  vv->ln = 0 ;
		  y = vv->ln + ln + 1 ;
		  if (y < bufSize)
		    y = bufSize ;
		  vv->buf = malloc (y) ;
		}
	    }
	  memcpy (vv->buf + vv->ln, nam, ln + 1) ;
	  vv->ln += ln ;	  
	}
    }
  /* export remnants */
  for (i = 0 ; i < nPrefix ; i++)
    {
      vv = vvv[i] ;
      if (vv->ln)
	channelPut (pp->ventChans[vv->writer], vv, VV*) ;
      else
	free (vv->buf) ;
    }

  /* signal completion to all the writers */
  for (i = 0 ; i < pp->nWriters ; i++)
    channelPut (pp->ventChans[i], 0, VV*) ;
  return ;
}

/*************************************************************************************/

static void dnaChanFastaReader (PP *pp)
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
      if (debug) fprintf (stderr, "dblParseFastaFile calls dnaChanBuffer k=%d\n", k) ;
      carryOver = dnaChanBuffer (fd, tt->buf, size -  (1 << 21),  carryOverBuf, carryOver) ;
      if (tt->buf[0] == 0)
	break ;  
      n = 1 ;
      if (pp->max_threads) /* signal that data is ready for thread tt */
	{
	  if (debug) fprintf (stderr, "dblParseFastaFile prepared data and calls channelPut(chan-%d)",k) ;
	  channelPut (tt->chan, &n, int) ; /* data ready */
	}
      else  /* direct procedural call */
	dblDoCountTelomere (tt) ;
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
	  if (debug) fprintf (stderr, "dblParseFastaFile waits on channelGive(done)") ;
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

static int dnaChanRun (PP *pp)
{
  AC_HANDLE h  = ac_new_handle () ;
  int i, n ;
  CHAN *done = channelCreate (24, int, h) ;
  int size = (1 << 25) ; /* 32M */ 
  BOOL debug = FALSE ;

  if (pp->blockSize * (1 << 20) > size)
    size = pp->blockSize * (1 << 20) ;

  if (debug) channelDebug (done, TRUE, "chan-done") ;
	 			     
  fprintf (stderr, "--- data preparation done : %s\n", timeShowNow ()) ;
      
  pp->doneChan = channelCreate (2*pp->nWriters, int, h) ;

  maxTh = pp->max_threads ;
  if (maxTh < 1) 
    pp->max_threads = maxTh = 1 ;

  /* declare the reader channel, we can only have one
   * but it can prepare several buffers, limited 
   * by the depth of the corresponding channel
   * so we start it directly
   * It export on bufChan
   */
  tp->bufChan = channelCreate (24, int, h) ;

  /* the buffers are analysized by ventialtors
   * declare the ventilation channel
   * each ventilator export blocks with a given signature
   * in 4 channels with different famillies signatures
   */
  for (i = 0 ; i < pp->nWriters ; i++)
    tp->ventChan[i] = channelCreate (24, int, h) ;
  /* the ventialted buffers are appendedon open files 
   * writers 
   */
  for (i = 0 ; i < 1 ; i++)
    wego_go (dnaChanReadBuf, pp, PP) ;
  for (i = 0 ; i < pp->nVentilators ; i++)
    wego_go (dnaChanVentilate, pp, PP) ;
  for (i = 0 ; pp->nWriters < 1 ; i++)
    { /* the pp structure is copied
       * so pp->iWriter is private to the instance 
       */
      pp->iWriter = i ;
      wego_go (dnaChanWriter, pp, PP) ;
    }
  
  dnaChanFillBuffer (pp) ;
  
  i = pp->nWriters ;
  while (i > 0)
    {
      if (debug) 
	fprintf (stderr, "dnaChanRun waits on channelGive(done)") ;
      n = channelGive (pp->doneChan, &n, int) ;
      i-- ;    
      if (debug) 
	fprintf (stderr, "dnaChanRun: writer %d done", n) ;
    }

  ac_free (h) ;	
  return nn ;
} /* dnaChanRun */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: dnachan [-help]... \n"
	   "//      try: -h -help --help \n"
	   "// Example:  dnachan -i f.fastq -I fastq -o out \n"
	   "// -i fileName : the name of a sequence file, possibly .gz, to be analyzed \n"
	   "// -I input_file_format : one of\n"
	   "//    -I fasta: [default] the input is a fasta file\n"
	   "//    -I fastq: the input is a standard fastq file\n"
	   "//    -I fastc: fastc file (contains #multipliers for repeated sequences)\n"
	   "//    -I csfasta: SOLiD color coded transition fasta file\n"
	   "//    -I csfastc: idem (contains #multipliers for repeated sequences)\n"
	   "//    -I raw: the input is a raw file, one sequenc per line, no identifers\n"
	   "// -silent : suppress stderr status report\n"
	   "// -gzi : force gunzipping the input (useful when piping from stdin)\n"
	   "// -gzo : the output file will be gziped\n"
	   "// -max_threads <integer>: [default 1] try 4, 8..., machine dependant\n"
	   "//   -split <int> : split the output files in chunks of specified number of sequences. We recommend 5000000. \n"
	   "//     The multiple output files are named <output_file>.1 .2 .3 ... .n\n"
	   "//   -splitMb <int> : split the output files in chunks of specified number of Mega bases.  We recommend 250. \n"
	   "//     The multiple output files are named <output_file>.1 .2 .3 ... .n\n"
	   "//   -splitByPrefix n : n=1,2,3 or 4, split the output in 4^n files according to the first n letters\n"
	   "//     Note that to merge input file no option is required, pipe them all to stdin, UNIX way\n"

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

  if (1)
    {
      p.dnaFormat = FASTC ; /* default */
      const char *ccp = 0 ;
      if (getCmdLineOption (&argc, argv, "-I", &ccp))
	{
	  if (!strcasecmp (ccp, "raw"))
	    {
	      p.dnaFormat = RAW ;
	      p.suffix = "raw" ;
	    }
	  else if (!strcasecmp (ccp, "FASTC"))
	    {
	      p.dnaFormat = FASTC ;
	      p.suffix = "fastc" ;
	    }
	  else if (!strcasecmp (ccp, "FASTA"))
	    {
	      p.dnaFormat = FASTA ;
	      p.suffix = "fasta" ;
	    }
	  else if (!strcasecmp (ccp, "FASTQ"))
	    {
	      p.dnaFormat = FASTQ ;
	      p.suffix = "fastq" ;
	    }
	  else
	    messcrash ("Sorry, unknown format in option -I %s, try dnabloom -h", ccp) ;
	}
    }

  p.gzi = getCmdLineOption (&argc, argv, "-gzi", 0);
  p.gzo = getCmdLineOption (&argc, argv, "-gzo", 0);

  p.silent = getCmdLineOption (&argc, argv, "-silent", 0);

  getCmdLineInt (&argc, argv, "-split", &sx.split) ;
  getCmdLineInt (&argc, argv, "-splitMb", &sx.splitMb) ;
  getCmdLineInt (&argc, argv, "-splitByPrefix", &sx.splitByPrefix) ;

  p.max_threads = 0 ;
  getCmdLineInt (&argc, argv, "-max_threads", &p.max_threads);
  p.nWriters = 4 ;
  p.blockSize = 32 ;
  
  if (argc != 1)
    {
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf, argc, argv) ;
    }

  fprintf (stderr, "// %s start\n", timeShowNow()) ;
  aceInWaitAndRetry (0) ; /* does nothing, but triggers the linker on LINUX */
  if (p.max_threads)
    {
      wego_max_threads (p.max_threads + 2) ;
    }

  dnaChanRun (&p) ;
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
