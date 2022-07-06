#include "ac.h"
#include "matrix.h"

/* construct the table of representations of the basic classical Lie superalgebras 
 */

typedef struct saStruct
{
  AC_HANDLE h ;
  
  const char *type ; /* A,B.C,D,F,G */
  int m, n, rank ;
  int D1, D2 ; /* number of even and odd generators of the adjoiunt rep */
  Array dd ;   /* dims of all the submodules */
  MX chi ;
  MX *mu ; /* the matrices */
  MX *Cartan ;
  ww HW ; /* HIGHEST WEIGHT */
  Array www ; /* weights */
} SA ; /* SuperAlgebra struct */

typedef struct wStruct
{
  int x[8] ;
} WW ; /* representation weigth vector */


/******************************************************************************************************/


static void getCartan (SA *sa)
{
  MX Cartan = 0 ;
  int m = sa->m, n = sa->n, r ;
  switch ((int)sa->type[0])
    {
      case 'A':
	if (m<0 || n<0 || m+n<2)
	  messcrash ("Type A(m/n): m+n should be >=2 m=%d n=%d", m,n) ;
	r = sa->rank = m + n - 1 ;
	Cartan = mxCreate (sa->h, "Cartan" , MX_INT, r, r, 0) ;
	{
	  int i, xx [r][r] ;
	  memset (xx, 0, sizeof (xx)) ;
	  for (i = 0 ; i < r ; i++)
	    {
	      xx[i][i] = 2 ;
	      if (i > 0) xx[i-1][i] = -1 ;
	      if (i < r-1) xx[i+1][i] = -1 ;
	    }
	  if (m*n > 0)
	    {
	      xx[m][m] = 0 ;
	      if (m > 0) xx[m-1][m] = -1 ;
	      if (m < r-1) xx[m+1][m] = 1 ;
	    }
	  mxSet (Cartan, xx) ;
	}
	break ;
	
	default:
	  messcrash ("The Cartan matrix of a superalgebra of type %c is not yet programmed, sorry.", sa->type) ;
    }
  mxShow (Cartan) ;

}  /* getCartan */

/******************************************************************************************************/

static voidKY getHighestWeight (SA *sa)
{
  WW *hw ;
  Array wws ;

  wws = SA->wws = ArrayHandleCreate (64, WW, sa->h) ;
  ww = arrayp (wws, 0, WW) ;
  if (sa->wws)
    sscanf (sa->wws, "%d:%d:%d:%d:%d:%d:%d:%d", hw->x[0] , hw->x[1] , hw->x[2] , hw->x[3] , hw->x[4] , hw->x[5] , hw->x[6] , hw->x[7] ) ;
  else
    messcrash ("Missing argument --s, expectirn -wws 1:0:2;...    giving the Dynkin lables of the highest weight ") ;
  
} /* getHighestWeight */

/******************************************************************************************************/

static BOOL demazure (SA *sa)
{
  BOOL new = FALSE ;
  
  
  
  return new ;
} /* demazure */

/******************************************************************************************************/
/******************************************************************************************************/


int main  (int argc, const char **argv)
{
   AC_HANDLE h = handleCreate () ;
   SA sa ;
   freeinit () ;
   sa.h = h ;

   sa.type = "toto" ;
   memset (&sa, 0, sizeof (sa)) ;
   getCmdLineInt (&argc, argv, "-m", &sa.m) ; 
   getCmdLineInt (&argc, argv, "-n", &sa.n) ;
   getCmdLineOption (&argc, argv, "-type", &sa.type) ; 
   if (sa.m < 0) messcrash ("argument -m m of SU(m/n) must be positive or null") ;
   if (sa.n < 0) messcrash ("argument -n n of SU(m/n) must be positive or null") ;
   if (sa.m + sa.n < 2) messcrash ("In SU(m/n) m+n must be >2") ;
   if (sa.m == sa.n) messcrash ("SU(m/m) must be treated separatelly") ; 
   sa.rank = sa.m + sa.n ;

   getCartan (&sa) ;
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
