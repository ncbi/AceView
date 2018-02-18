/*
 * Authors: Jean Thierry-Mieg, NCBI, 
 * Oct 2009
 * dnastitch paired end fragments
 *   Input: raw, fasta, fastq, solid
 *   Output: raw, fasta, solid
 *   Filters: length, entropy
 *
*/

#define BITSET_MAKE_BITFIELD   
#include "regular.h"
#include <ac.h>
#include <acedna.h>
#include <aceio.h>
#include <dict.h>
#include <vtxt.h>
#include <math.h>

typedef struct gxStruct { 
  AC_HANDLE h ; 
  AC_DB db ;           /* ACEDB database containing the semantics of the experiment */
  const char *title ;
  const char *inFileName ;
  const char *outFileName ;
  const char *runListFileName ;
  const char *run ;
  Array dna1, dna2, dna2R, err ;
  int overlap ;
  int fragment_length ;
  BOOL small ;
  BOOL gzi, gzo ;
} GX ;



/*************************************************************************************/
/*************************************************************************************/

static BOOL gxDoStitch (GX *gx, ACEOUT ao, char *line, char *nameBuf)
{
  int i, j, n, n1, n1a, n2 ;
  char *sep ;
  register char *cp, *cq ;
  Array d1 = gx->dna1, d2 = gx->dna2, d2R = gx->dna2R ;

  sep = strstr (line, "><") ;
  if (sep)
    {
      /* encode the first segments */
      n = strlen (line) ; 
      array (d1, n+4, char) = 0 ; /* make room */
      array (d2, n+4, char) = 0 ;
      *sep = 0 ;
      cq = line ;
      n = strlen (cq) ; 
      arrayMax(d1) = n ;
      for (cp = arrp(d1, 0, char) ; *cq ; cp++, cq++)
	*cp = dnaEncodeChar[(int)*cq] ;
      *++cp = 0 ; *++cp = 0 ;

      /* encode the complement of the second segments */
      cq = sep + 2 ;
      n = strlen (cq) ; 
      arrayMax(d2) = arrayMax(d2R) = n ;
      for (cp = arrp(d2, 0, char) ; *cq ; cp++, cq++)
	*cp = dnaEncodeChar[(int)*cq] ;
      *++cp = 0 ; *++cp = 0 ;
      memcpy (arrp(d2R, 0, char), arrp(d2, 0, char), n+2) ;
      reverseComplement (d2R) ;
      

      n1 = arrayMax (d1) ; 
      n2 = arrayMax (d2) ;
      if (gx->fragment_length)
	{
	  char nBuf[1 + gx->fragment_length] ;
	  
	  memset (nBuf, 'N', gx->fragment_length) ;
	  nBuf[gx->fragment_length] = 0 ;

	  dnaDecodeArray (d1) ;
	  dnaDecodeArray (d2R) ;
	  
	  /* export the name and the first read, the padding, the complement of the second read */
	  nBuf[gx->fragment_length - n1 - n2] = 0 ;
	  aceOutf (ao, "%s\n%s%s%s\n"
		   , nameBuf
		   , arrp(d1, 0, char)
		   , nBuf
		   , arrp(d2R, 0, char)
		   ) ;
	  /* export the missing n */
	  nBuf[gx->fragment_length - n1 - n2] = 'N' ;

	  return TRUE ;
	}

      /* scan segment 2, check if words  are in catalog */ 
      n1a = gx->small ? 0 : n1 - 8 ;
      for (i = 0, cp = arrp(d1,0,char) ; i <= n1a ; cp+=8, i+=8)
	{
	  for (j = 0, cq = arrp(d2R,j,char) ; j < n2 - 8 ; cq++, j++)
	    if (memcmp (cp,cq, 8) == 0)
	      {
		int ok = 0 ;
		int x1 = i + 1, x2 = i + 8 ; 
		int a1 = j + 1, a2 = j + 8 ;
		char cc, *cx1 = cp, *cx2 = cp + 7, *ca1 = cq, *ca2 = cq + 7 ;

		while (x1 >=  1 && a1 >=  1 && *cx1 == *ca1) { x1-- ; a1-- ; cx1-- ; ca1-- ;}
		x1++ ; a1++ ; cx1++ ; ca1++ ;

		while (x2 <= n1 && a2 <= n2 && *cx2 == *ca2) { x2++ ; a2++ ; cx2++ ; ca2++ ;}
		x2-- ; a2-- ; cx2-- ; ca2-- ;

		if (x2 - x1 < gx->overlap) continue ;
		else if (x1 == 1  && a2 == n2) ok = 1 ;
		else if (x2 == n1 && a2 == 1)  ok = 2 ;
		else if (a1 == 1  && a2 == n2) ok = 3 ;
		else if (x1 == 1  && x2 == n1) ok = 4 ;
		else continue ;

		/* decode */
		dnaDecodeArray (d1) ;
		dnaDecodeArray (d2) ;

		/* export the overlap */
		cc = *(cx2+1) ;  *(cx2+1) = 0 ;
		cp = cx1 ;
		while (*cp) { *cp = ace_upper(*cp) ; cp++ ;}
		aceOutf (ao, "%s\t%s\t%d\t%d\t%s", gx->run,nameBuf, fastcMultiplicity (nameBuf,0,0),x2 - x1 + 1, cx1) ;
		*(cx2+1) = cc ;

		/* export the coordinates */
		aceOutf (ao, "\t%d\t%d\t%d\t%d\t"
			 , x1, x2, a1, a2
			 ) ;

		/* export the left exit overlap */
		if ((ok == 1 || ok == 4) && a1 > 1) aceOutf (ao, "%s\t", cx2 + 1) ;
		else aceOutf (ao, "\t") ;

		/* export the right exit overlap */
		if ((ok == 1 || ok == 3) && x2 < n1) aceOutf (ao, "%s\n", arrp (d2R, arrayMax(d2) - a1 + 1, char)) ;
		else aceOutf (ao, "\n") ;

		return TRUE ; 
	      }
	}
    }
  
  return FALSE ;
} /* gxDoStitch */

/*************************************************************************************/

static void gxStitch (GX *gx)
{
  AC_HANDLE h = 0 ;
  ACEIN ai = aceInCreate (gx->inFileName, gx->gzi, h) ;
  ACEOUT ao = aceOutCreate (gx->outFileName, ".fragments", gx->gzo, h) ;
  int state = 0, nn = 0 ;
  char *cp, nameBuf[1024] ;

  while (aceInCard (ai))
    {
      cp = aceInPos (ai) ;
      if (! cp)
	continue ;
      if (*cp == '#')
	{
	  if (0) aceOutf (ao, "%s\n", cp) ;
	  continue ;
	}
      if (state == 0)
	{
	  if (*cp == '>')
	    {
	      strncpy (nameBuf, cp, 1024) ;
	      state = 1 ;
	    }
	  continue ;
	}
      if (state == 1)
	{
	  if (gxDoStitch (gx, ao, cp, nameBuf))
	    nn++ ;
	  state = 0 ;
	  continue ;
	}
    }

  ac_free (h) ;
} /* gxStitch */

/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// dnastitch: stitch paired-end DNA fragments\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, Mars 2014, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Read a DNA file in fatc format\n"
	    "// If the 2 reads overlap significantly, export a merged fragment\n"
	    "// together with handling the overhanging adaptor consensus\n"
	    "// for fragments shorter than the sequencing reads\n"
	    "//\n"
	    "// Syntax:\n"
	    "// dnastitch [options]\n"
	    "// The options are all optional and may be specified in any order\n"
	    "// Expected fragment length\n"
	    "//   -fragment_length <int>\n"
	    "//    in this case, do not stitch, just pad with n and complement the second read\n"
	    "// DNA Input\n"
	    "//   -i input_file: [default: stdin] fastc sequence file to analyze\n"
	    "// DNA Output\n"
	    "//   -o output_file: [default: stdout] processed data\n"
	    "// gzip (recommended)\n"
	    "//   -gzi -gzo : gzi, decompress the input file; gzo, compress the output file\n"
	    "// \n"
	    "// Sequence Input Format: Fastc\n"
	    "//   fastc format is generated by the dna2dna format converter\n"
	    "//   Example\n"
	    "//   >s1#12\n"
	    "//     ATGCTGTC><TTCGACAG\n"
	    "//   Export\n"
	    "//     atGCTGTCgaa or atgctgtcnnnnnnnnnctgtcgaa\n"
	    "//\n"
	    "// Caveat:\n"
	    "//   Lines starting with '#' are considered as comments and dropped out\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  dna2dna -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  GX gx ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&gx, 0, sizeof (GX)) ;
  gx.h = h ;

  /* optional arguments */
  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  /* input output redirections */
  gx.gzi = getCmdLineOption (&argc, argv, "-gzi", 0);
  gx.gzo = getCmdLineOption (&argc, argv, "-gzo", 0);
  getCmdLineOption (&argc, argv, "-i", &gx.inFileName) ;
  getCmdLineOption (&argc, argv, "-o", &gx.outFileName) ;
  gx.run = "AnyRun" ;
  getCmdLineOption (&argc, argv, "-run", &gx.run) ;

  gx.overlap = 16 ;
  getCmdLineInt (&argc, argv, "-overlap", &gx.overlap) ; 
  getCmdLineInt (&argc, argv, "-fragment_length", &gx.fragment_length) ;
  gx.small = getCmdLineOption (&argc, argv, "-small", 0);

  gx.dna1 = arrayHandleCreate (1024, char, h) ;
  gx.dna2 = arrayHandleCreate (1024, char, h) ;
  gx.dna2R = arrayHandleCreate (1024, char, h) ;
  gx.err =  arrayHandleCreate (64, A_ERR, h) ;
  gxStitch (&gx) ;

  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// done: %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

