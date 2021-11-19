#include "ac.h"

typedef enum { FASTA=0, FASTC, FASTQ, CSFASTA, CSFASTC, CCFA, CCFAR, RAW,  CRAW, TAG_COUNT, CTAG_COUNT, COUNT, LETTER, LETTERPLOT, LETTERSHOW } DNAFORMAT ;

#define LNV 25

typedef struct mirStruct { 
  AC_HANDLE h ;
  const char *inFileName, *outFileName, *mirFileName ;
  BOOL gzi, gzo ;
  DNAFORMAT in ; 
  char fileType [16] ;
  Array mirs ;
  int v5[4][4*LNV], v3[4][4*LNV] ;
  DICT *mirDict ;
} MIR ;

typedef struct wStruct { 
  int w, ln, mult, a, t, g, c ;
} WW ;

/*************************************************************************************/
/*************************************************************************************/

static int mirWwOrder (const void *a, const void *b)
{
  const WW *up = (const WW *) a, *vp = (const WW *) b ;
  int n ;
  
  n = up->ln - vp->ln ; if (n) return -n ;
  n = up->w - vp->w ; if (n) return n ;

  return 0 ;
}

/*************************************************************************************/

static int mirMultOrder (const void *a, const void *b)
{
  const WW *up = (const WW *) a, *vp = (const WW *) b ;
  int n ;
  
  n = up->mult - vp->mult ; if (n) return -n ;
  n = up->w - vp->w ; if (n) return n ;

  return 0 ;
}

/*************************************************************************************/

static void mirRegister (MIR *mir, WW *up, char *cp, int mult)
{
  int k = -4 ;
  int kk = 0 ;
  int *nn ;
  
  switch ((int) *cp)
    {
    case 'a' : case 'A': kk = 0 ; up->a += mult ; break ;
    case 't' : case 'T': kk = 1 ; up->t += mult ; break ;
    case 'g' : case 'G': kk = 2 ; up->g += mult ; break ;
    case 'c' : case 'C': kk = 3 ; up->c += mult ; break ;
    }
  
  nn = mir->v5[kk] ;
  cp-- ;
  while (k += 4, *++cp)
    {
      switch ((int) *cp)
	{
	case 'a' : case 'A': nn[k+ 0] += mult ; break ;
	case 't' : case 'T': nn[k+ 1] += mult ; break ;
	case 'g' : case 'G': nn[k+ 2] += mult ; break ;
	case 'c' : case 'C': nn[k+ 3] += mult ; break ;
	}
    }
} /* mirRegister */

/*************************************************************************************/

static void mirExportCounts (MIR *mir)
{
  int i ;
  int k, nn, kk ;
  WW *up ;
  int NNN = 0, NN[4] ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (mir->outFileName, ".mirVector", mir->gzo, h) ; 

  for (kk = 0 ; kk < 4 ; kk++)
    {
      aceOutf (ao, "\n\n\n") ;
      
      int *nn = mir->v5[kk] ;
      NN[kk] =  nn[0] + nn[1] + nn[2] + nn[3] ;
      NNN += NN[kk] ;
    }



  for (kk = 0 ; kk < 4 ; kk++)
    {
      char V5[LNV+2] ;
      if (10*NN[kk] < NNN)
	continue ;
      aceOutf (ao, "\n\n\n") ;
      
      int *nn = mir->v5[kk] ;
      for (i = 0 ; i < LNV ; i++)
	{
	  int N =  nn[4*i+0] + nn[4*i+1] + nn[4*i+2] + nn[4*i+3] ;
	  char cc = 'N' ;

	  if (nn[4*i+0] > N/2) cc = 'A' ;
	  if (nn[4*i+1] > N/2) cc = 'T' ;
	  if (nn[4*i+2] > N/2) cc = 'G' ;
	  if (nn[4*i+3] > N/2) cc = 'C' ;

	  if (N == 0) 
	    N = 1 ;
	  aceOutf (ao, "%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%c\n"
		  , i, N
		  , nn[4*i+0]
		  , nn[4*i+1]
		  , nn[4*i+2]
		  , nn[4*i+3]
		  , (100.0 * nn[4*i+0])/N
		  , (100.0 * nn[4*i+1])/N
		  , (100.0 * nn[4*i+2])/N
		  , (100.0 * nn[4*i+3])/N
		  , cc
		  ) ;
	  V5[i] = cc ;
	}
      V5[i] = 0 ;
      aceOutf (ao, "\n%s\t%d\tVECTOR_5\n\n", V5, NN[kk]) ;
    }

  arraySort (mir->mirs, mirMultOrder) ;
  if (1)
    for (k = 0, up = arrp (mir->mirs, 0, WW) ; k < 128 ; k++, up++)
      aceOutf (ao, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", dictName (mir->mirDict, up->w), up->ln,up->mult, up->a,up->t,up->g,up->c) ;

  ac_free (h) ;
} /* mirExportCounts */

/*************************************************************************************/

static void mirCountOne (MIR *mir, char *word, int mult)
{
  int i ;
  DICT *dict = mir->mirDict ;
  Array mirs = mir->mirs ;
  WW *up ;
  int iMax = arrayMax (mirs) ;
  
  for (i = 0, up = arrp (mirs, 0, WW) ; i < iMax ; up++, i++)
    {
      const char *ccp ;
      char buf [LNV], *cq ;

      ccp = dictName (dict, up->w) ;
      cq = strstr (word, ccp) ;
      if (cq)
	{
	  int k = cq - word ; /* offset */
	  int wln = strlen (word) ;
	  up->mult += mult ;  /* we found this mir */
	  if (k + up->ln < wln)
	    {
	      cq = word + k + up->ln ;
	      if (k + up->ln + LNV < wln)
		word [k + up->ln + LNV] = 0 ;
	      mirRegister (mir, up, cq, mult) ;
	    }
	  if (0 && k > 0)
	    {
	      memcpy (buf, cq + wln, LNV) ;
	      buf[LNV] = 0 ; /* to be sure */ 
	      mirRegister (mir, up, buf, 0) ;
	    }
	  break ;
	}
    }
} /* mirCountOne */

/*************************************************************************************/

static void mirGetCounts (MIR *mir)
{
  AC_HANDLE  h = ac_new_handle () ;  
  ACEIN ai = aceInCreate (mir->inFileName, mir->gzi, h) ;
  aceInSpecial (ai, "\n\t") ;
  char buf[512] ;
  DNAFORMAT in = mir->in ;
  int mult = 1 ;

  buf[511] = 0 ;
  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      if (! cp || *cp == '#')
	continue ;
      switch (in)
	{
	case TAG_COUNT:
	  strncpy (buf, cp, 510) ;
	  aceInStep (ai, '\t') ;
	  if (aceInInt (ai, &mult))
	    mirCountOne (mir, buf, mult) ; ;
	  break ;
	case FASTC:
	  if(cp[0]=='>')
	    mult = fastcMultiplicity (cp, 0, 0) ;
	  else
	    { mirCountOne (mir, cp, mult) ;  mult = 1 ; }
	  break ;
	case FASTA:
	  if(cp[0] != '>')
	   mirCountOne (mir, cp, mult) ;  
	  break ;
	case RAW:
	    mirCountOne (mir, cp, mult) ; 
	  break ;
	default:
	  messcrash ("this -I format is not yet programmed in mir.c, sorry\n") ;
	}
    }
} /* mirGetFrequentMir */

/*************************************************************************************/

static void mirGetFrequentMir (MIR *mir)
{
  AC_HANDLE  h = ac_new_handle () ;  
  ACEIN ai = aceInCreate (mir->mirFileName, mir->gzi, h) ;
  int k, nn = 0 ;
  DICT *dict ; 
  Array aa ;
  WW *up ;

  aa = mir->mirs = arrayHandleCreate (300, WW, mir->h) ;
  mir->mirDict = dict = dictHandleCreate (300, mir->h) ; 
  if (ai)
    while (aceInCard (ai))
      {
	const char *ccp = aceInWord (ai) ;
	if (! ccp || *ccp == '#')
	continue ;
	dictAdd (dict, ccp, &k) ;
	up = arrayp (aa, nn++, WW) ;
	up->w = k ;
	up->ln = strlen (ccp) ;
      }
  
  arraySort (aa, mirWwOrder) ;
  if (0)
    for (k = 0, up = arrp (aa, 0, WW) ; k < nn ; k++, up++)
      printf ("%s\t%d\n", dictName (dict, up->w), up->ln) ;

  ac_free (h) ;
} /* mirGetFrequentMir */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// mir: Manalysis of short (mir) RNAs\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Analyze short RNAs \n"
	    "//\n"
	    "// Syntax:\n"
	    "// mir [options]\n"
	    "// The options are all optional and may be specified in any order\n"
	    "// DNA Input\n"
	    "//   -i input_file: [default: stdin] sequence file to analyze\n"
	    "//      If the file is gzipped : -i f.fasta.gz\n, it will be gunzipped"
	    "//      Any pipe command prefixed by < is accepted:  -i \'< obj_get f.fasta.gz | gunzip\'"
	    "//   -I input_format: [default: fasta] format of the input file, see formats list below\n"
	    "// Formats\n"
	    "//   raw : one DNA sequence per line in 16 letter IUPAC code, no identifiers\n"
	    "//   raw<n> :  Input only, the DNA is in column Number (example: raw1 == raw, raw2, raw3 ...)\n"
	    "//   raw<n1>yn<n2> :  In addition flag in column n2 must be Y\n"
	    "//      example: raw9yn22  imports DNA from column 9 if column 22 says Y\n"
	    "//      in raw or raw1 mode, the identifier of the sequence is assumed to be in column 2\n"
	    "//      in raw<n> mode, n>1, the identifier of the sequence is assumed to be in column 1\n"
	    "//   tc : Tag Count : DNA sequence in 16 letter IUPAC code 'tab' multiplicity\n"
	    "//   fasta : >identifier 'new line' DNA in 16 letter IUPAC code\n"
	    "//     As input, the DNA can be spread over any number of lines\n"
	    "//     As output, the DNA is exported on a single line except if -maxLineLn <number> is specified \n"
	    "//     If the input does not provide identifiers, e.g. in the raw or tc input formats, the option\n"
	    "//          -prefix <txt> serves to construct the identifiers as 'txt.1 txt.2 txt.3 ...\n"
	    "//   fastc : >identifier#tag_count 'new line' DNA in 16 letter IUPAC code\n"
	    "//     Same as fasta, but the identifier includes the tag count multiplicity after the # symbol\n"
	    "//   fastq : @identifier 'new line'  DNA 'new line' +same_identifier 'new line' quality\n"
	    "//     Qualities are specific of the fastq format and dropped when exporting in any other format\n"
	    "//     Hence -O fastq is only possible if -I fastq, for example to split or filter the input file\n"
	    "//     or if we provide \n"
	    "//\n"
	    "// SOLiD formats and options\n"
	    "//     Recall that a SOLiD sequence starts with an Anchor letter followed by transitions, a natural\n"
	    "//     code corresponding to their ligation sequencing protocol\n"
	    "//   csfasta : LifeTech SOLiD color fasta format: anchor (usually T) followed by transitions 0123\n"
	    "//   csfastc : idem, but with a #multiplicity factor in the sequence identifiers\n"
	    "//   ccfa : modified LifeTech SOLiD format, with 0123 replaced by acgt\n"
	    "//   ccfar : modified LifeTech SOLiD format with 0123 replaced by tgca\n"
	    "//     The advantage of the ccfa format is that it is accepted by standard DNA programs such as\n"
	    "//     Blast/Blat or Bowtie. For example one may align a cfa file onto the direct strand of the genome\n"
	    "//     by transforming both sequences in ccfa format: this is much preferable to decoding cfa to\n"
	    "//     fasta then aligning to the fasta genome, because in such a naive scheme, any missmatch\n"
	    "//     irreversibly disrupts the homology between the genome and the SOLiD sequence\n"
	    "//     downstream of the mismatch\n"
	    "//     The ccfar format is needed to align sequences to the complementary strand. A possible\n"
	    "//     protocol to align a SOLiD file on the genome would be to reformat the file in ccfa,\n"
	    "//     to reformat the genome in both ccfa and ccfar and to run the aligner twice.\n"
	    "//     Preferably, one may align just once using the aceview 'align' program or any other\n"
	    "//     aligner which understands the native SOLiD cfa format and implements the SOLiD error\n"
	    "//     correction mechanism.\n"  
	    "//   -n2a : 'n', 'N' or '.' characters in the input sequence file are converted to 'a', and become\n"
	    "//     exportable in ccfa format. This option should be used for converting a genome file, which\n"
	    "//     usually contains stretches of N, to SOLiD ccfa and ccfar color formats.\n"
	    "//   -jumpAnchor : removes the SOLiD Anchor letter. The resulting color fasta file is\n"
	    "//     no longer decodable into fasta format, but now matches exactly the SOLiD color fasta genome.\n"
	    "//     If the anchor letter was not removed in this way, it would not match the corresponding position \n"
	    "//     on the ccfa genome, resulting in a spurious mismatch in naive aligner programs, such as \n"
	    "//     BLAST/BLAT/Bowtie which, as of December 2009, do not process the color format. \n"
	    "//\n"
	    "// ACTIONS\n"
	    "//   -m mirFileName : file of say 250 very frequenct miRs\n"
	    "//      format : one atgc word per line\n"
	    "//      Search these words as exact in the input file and export the overhangs\n"
	    "// Examples:\n"
	    "// Caveat:\n"
	    "//   Lines starting with '#' are considered as comments and dropped out\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  mir --help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/

static BOOL checkFormat (const char *io, DNAFORMAT *ip, const char *ccp, char *ftype)
{
  int i ;
  const char **f ; 
  const char *ff[] = { "fasta", "fastc", "fastq", "csfasta",  "csfastc", "ccfa", "ccfar", "raw", "craw", "tc", "ctc", "count", 0 } ; 

  for (i = 0 , f = ff ; *f ; i++, f++)
    if (! strcmp (*f, ccp))
      { 
	*ip = i ; 
	if (*io == 'O')
	  strcpy (ftype, ff[i]) ;
	return TRUE ;
      }

  usage (messprintf ("Unknown format in option -%s %s", io, ccp)) ;

  return FALSE ;
} /* checkFormat */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  const char *ccp = 0 ;
  MIR mir ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&mir, 0, sizeof (MIR)) ;
  mir.h = h ;

  /* freeOut (0) ;  needed be the linker 2017_02_26 */
  /* optional arguments */

  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  getCmdLineOption (&argc, argv, "-i", &mir.inFileName) ;
  getCmdLineOption (&argc, argv, "-m", &mir.mirFileName) ;
  getCmdLineOption (&argc, argv, "-o", &mir.outFileName) ;

  if (getCmdLineOption (&argc, argv, "-I", &ccp))
    {
      checkFormat ("I", &mir.in, ccp, mir.fileType) ;
    }
  if (argc > 1)
    usage (messprintf("sorry unknown arguments %s", argv[argc-1])) ;

  if (mir.mirFileName)
    {
      mirGetFrequentMir (&mir) ;
      mirGetCounts (&mir) ;
      mirExportCounts (&mir) ;
    }
  

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
