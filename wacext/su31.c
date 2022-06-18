#include "ac.h"
#include <complex.h>
#include "matrix.h"

/* Create june 2022
 * consruct the superalgebra sl(3/1)
 */

#include "ac.h"
#include <complex.h>
#include "matrix.h"

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// su21: Construction of su(2/1) representations and Feynman diagrams\n"
	    "// Authors: Jean Thierry-Mieg, NCBI, 2020-, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Construct the matrices of irreducible and indecomposable representations\n"
	    "// Construct the Casimirs, super Casimirs, Gorelik ghost Casimir\n"
	    "// \n"
	    "// Also compute the anomalies and the feynman diagrams supporting my JHEP su(2/1) papers\n"
	    "//\n"
	    "// Syntax:\n"
	    "// su31 [options]\n"
	    "//   [] [-h] [-help] [--help] : this message\n"
	    "// A: Representations\n"
	    "//   su21 -a1 <int> -a2 <int> -b <int> [-N <int>]\n"
	    "//     export the matrices, Casimirs and verifications for the module with \n"
	    "//     Dynkin lables (a1,a2,b), a1,a2 positive integers, b signed integer\n"
	    "//     Number of generations N (N >= 2)\n"
	    "//       In theory, b can be any complex number,\n"
	    "//     for numerical convenience, we restrict here to signed integers\n"
	    "//     but the formulas like the Casimir eigen values are anlytic in b\n"
	    "//       When a or N are large, many outputs are suppressed, try first a<=3, N<=3\n"
	    "//\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  dna2dna --help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;

  freeinit () ;

  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  /*   BOOL SU3 = getCmdLineBool (&argc, argv, "-su3") ; */
  int NN = 0 ;
    
  getCmdLineInt (&argc, argv, "-N", &NN) ; /* Number of generations >= 2 */
  getCmdLineInt (&argc, argv, "-NN", &NN) ; /* synonim */

  int a1 = 0, a2 = 0, b = 0 ;

  getCmdLineInt (&argc, argv, "-a1", &a1) ;
  getCmdLineInt (&argc, argv, "-a2", &a2) ;
  getCmdLineInt (&argc, argv, "-b", &b) ;

  if (a1 < 0)
    usage ("SU(3) Dynkin weigth a should be a positiver integer") ; 
  if (a2 < 0)
    usage ("SU(3) Dynkin weigth a should be a positiver integer") ; 
  if (NN != 0 && NN < 2)
    usage ("The number of generations N should be an integer >= 2") ;

  if (!NN && (a1 + a2 || b))
    {
      /* 2021_03_18 
       * construct the 15 matrices for the generic irreps of su(3/1) with h.w. (a1,a2,b)
       * verify all commutations relations
       * compute the casimir tensors and operators 

       */
#ifdef JUNK
      Kasimirs (1,0,1, FALSE) ;  /* adjopint */
      Kasimirs (1,0,0, FALSE) ;  /* fundamental */
      Kasimirs (b,a1,a2, TRUE) ;  
  if (0) muInit (h) ;   /* init the 4x4 matrices */
  if (0) muInit2 (h) ;  /* init the 2-families 8x8 rotated matrices */
  if (NN >= 2) muInitNMarcu (a,b, NN) ;  /* init the 2-families 8x8 marcu indecomposable matrices */

#endif
    }
  /* always init, otherwise the gcc linker is unhappy */

  ac_free (h) ;   
  return 0 ;
}

