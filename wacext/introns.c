/*
 * Authors: Jean Thierry-Mieg, NCBI, 
 * Apr 2016
 * Statistical analysis of intron distribution
 *  type, banc de Poisson
*/

#define BITSET_MAKE_BITFIELD   
#include <ac.h>
#include <acedna.h>

typedef struct ppStruct { 
  AC_HANDLE h ; 
  const char *inFileName, *outFileName ;
  int n ; 
} PP ;

typedef struct ddStruct { int chrom, a1, a2, n ; BOOL isDown ; char type ;} DD ;

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
    fprintf  (stderr,
	      "// dna2dna: Multilingual DNA parser\n"
	      "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	      "// Purpose\n"
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
  PP pp ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&pp, 0, sizeof (PP)) ;
  pp.h = h ;


  /* optional arguments */

  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  getCmdLineOption (&argc, argv, "-i", &pp.inFileName) ;
  getCmdLineOption (&argc, argv, "-o", &pp.outFileName) ;



  fprintf (stderr , "Processed %d introns\n", 0) ;
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


