#include "ac.h"
#include "bitset.h"

/* construct the irrep as a sub verma module
 */

typedef struct saStruct
{ int M, N, rank ;
} SA ; /* SuperAlgebra struct */

int main  (int argc, const char **argv)
{
   AC_HANDLE h = handleCreate () ;
   SA sa ;
   freeinit () ;

   memset (&sa, 0, sizeof (sa)) ;
   getCmdLineInt (&argc, argv, "-m", &sa.M) ; 
   getCmdLineInt (&argc, argv, "-n", &sa.N) ; 
   if (sa.M < 0) messcrash ("argument -m m of SU(m/n) must be positive or null") ;
   if (sa.N < 0) messcrash ("argument -n n of SU(m/n) must be positive or null") ;
   if (sa.M + sa.N < 2) messcrash ("In SU(m/n) m+n must be >2") ;
   if (sa.M == sa.N) messcrash ("SU(m/m) must be treated separatelly") ; 
   sa.rank = sa.M + sa.N ;

#ifdef JUNK
  showCartan (g) ;
  constructSimpleRoots (g) ;
  constructAllRoots (g) ;
  constructRho (g) ;
  constructRotate (g) ;
  showRotate (g) ;
  showAllRoots (g, TRUE) ;

  constructPools (g) ;
  showPools (g) ;

  constructTests (g, 2) ;
  showTests (g) ;
  printf ("%d tests constructed\n", g->nTests) ;
  constructDistances (g) ;
  if (0) showDistances (g) ;
#endif
  messfree (h) ;
  printf ("A bientot\n") ;

  return 0 ;
}
