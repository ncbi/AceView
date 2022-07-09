#include "ac.h"
#include "matrix.h"

/* construct the table of representations of the basic classical Lie superalgebras 
 */

#define RMAX 8 

typedef struct wStruct
{
  int x[RMAX] ;
  BOOL ok[RMAX] ;
  int mult ;
  int k ;
} WW ; /* representation weigth vector */

typedef struct saStruct
{
  AC_HANDLE h ;
  
  BOOL table ;
  const char *type ; /* A,B.C,D,F,G */
  const char *DynkinWeights ; /* 1:0:2:.... */
  int m, n, rank ;
  int D1, D2 ; /* number of even and odd generators of the adjoiunt rep */
  Array dd ;   /* dims of all the submodules */
  MX chi ;
  Array Cartan ;
  WW hw ; /* HIGHEST WEIGHT */
  Array wws ; /* weights */
  DICT *dict ;
} SA ; /* SuperAlgebra struct */


/******************************************************************************************************/

static int wwOrder (const void *a, const void *b)
{
  const WW *up = (WW *)a ;
  const WW *vp = (WW *)b ;
  int n = 0 ;

  n = up->k - vp->k ;
  return n ;  
} /* wwOrder */

/******************************************************************************************************/


static void getCartan (SA *sa)
{
  Array Cartan = 0 ;
  int m = sa->m, n = sa->n, r = 0, rr ;
  switch ((int)sa->type[0])
    {
    case 'A':
      if (m<0 || n<0 || m+n<2)
	messcrash ("Type A(m/n): m+n should be >=2 m=%d n=%d", m,n) ;
      r = sa->rank = m + n - 1 ;
      rr = r*r ;
      Cartan = arrayCreate (rr, int) ;
      {
	int i ;
	for (i = 0 ; i < r ; i++)
	  {
	    array (Cartan, r*i + i, int) = 2 ;
	    if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	    if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	  }
	if (m*n > 0)
	  {
	      array (Cartan, r*m + m , int) = 0 ;
	      if (m > 0) array (Cartan, r*m + m - 1 , int) = 1 ;
	      if (m < r-1) array (Cartan, r*m + m + 1, int) = -1 ;
	  }
      }
      break ;
      
    case 'B':
      if (m < 2)
	messcrash ("Type B(m)): m should be at least 2:  m=%d n=%d", m,n) ;
      r = m ;
      rr = r*r ;
      Cartan = arrayCreate (rr, int) ;
      {
	int i ;
	for (i = 0 ; i < r ; i++)
	  {
	    array (Cartan, r*i + i, int) = 2 ;
	    if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	    if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	  }
      }
      array (Cartan, r*1 + 0, int) = -2 ;
      break ;
      
    case 'C':
      if (m < 2)
	messcrash ("Type C(m)): m should be at least 2:  m=%d n=%d", m,n) ;
      r = m ;
      rr = r*r ;
      Cartan = arrayCreate (rr, int) ;
      {
	int i ;
	for (i = 0 ; i < r ; i++)
	  {
	    array (Cartan, r*i + i, int) = 2 ;
	    if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	    if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	  }
      }
      array (Cartan, r*0 + 1, int) = -2 ;
      break ;
      
    case 'D':
      if (m<3)
	messcrash ("Type D(m): m should be >=3 and even m=%d n=%d", m,n) ;
      r = sa->rank = m ;
      rr = r*r ;
      Cartan = arrayCreate (rr, int) ;
      {
	int i ;
	for (i = 0 ; i < r ; i++)
	  {
	    array (Cartan, r*i + i, int) = 2 ;
	    if (i > 1) array (Cartan, r*i + i-1, int) = -1 ;
	    if (i> 0 && i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	  }
	array (Cartan, r*0 + 2, int) = -1 ;
	array (Cartan, r*2 + 0, int) = -1 ;
      }
      break ;
      
    case 'E':
      if (m < 6 || m > 8)
	messcrash ("Type E(m): m should be 6, 7 or 8 : m=%d n=%d", m,n) ;
      r = sa->rank = m ;
      rr = r*r ;
      Cartan = arrayCreate (rr, int) ;
      {
	int i ;
	for (i = 0 ; i < r ; i++)
	  {
	    array (Cartan, r*i + i, int) = 2 ;
	    if (i > 1) array (Cartan, r*i + i-1, int) = -1 ;
	    if (i > 0 && i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	  }
	array (Cartan, r*0 + 3, int) = -1 ;
	array (Cartan, r*3 + 0, int) = -1 ;
      }
      break ;
      
    case 'F':
      if (m != 4)
	messcrash ("Type F(4): m should be 4 : m=%d n=%d", m,n) ;
      r = 4 ;
      rr = r*r ;
      Cartan = arrayCreate (rr, int) ;

      {
	int i ;
	for (i = 0 ; i < r ; i++)
	  {
	    array (Cartan, r*i + i, int) = 2 ;
	    if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	    if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	  }
	array (Cartan, r*1 + 2, int) = -2 ;
      }
      break ;
      
    case 'G':
      if (m != 2)
	messcrash ("Type G(2)): m should be 2 m=%d n=%d", m,n) ;
      r = 2 ;
      rr = r*r ;
      Cartan = arrayCreate (rr, int) ;
      {
	array (Cartan, r*0 + 0, int) = 2 ;
	array (Cartan, r*0 + 1, int) = -3 ;
	array (Cartan, r*1 + 1, int) = 2 ;
	array (Cartan, r*1 + 0, int) = -1 ;
      }
      break ;
      
    default:
      messcrash ("The Cartan matrix of a superalgebra of type %c is not yet programmed, sorry.", sa->type) ;
    }
  if (1)
    {
      int i, j ;
      printf ("\n### Cartan Matrix type %s m=%d n=%d rank = %d\n", sa->type, m, n, r) ;
      for (i = 0 ; i < r ; i++)
	{
	  for (j = 0 ; j < r ; j++)
	    printf ("\t%d", arr (Cartan, r*i + j, int)) ;
	  printf ("\n") ;
	}
    }
  sa->Cartan = Cartan ;
}  /* getCartan */

/******************************************************************************************************/

static void getHighestWeight (SA *sa)
{
  WW *hw ;
  Array wws ;
  int k = 0 ;
  char buf[1024] ;
  
  wws = sa->wws = arrayHandleCreate (64, WW, sa->h) ;
  hw = arrayp (wws, 1, WW) ;
  if (sa->wws)
    k = sscanf (sa->DynkinWeights, "%d:%d:%d:%d:%d:%d:%d:%d", &hw->x[0] , &hw ->x[1] , &hw ->x[2] , &hw ->x[3] , &hw ->x[4] , &hw ->x[5] , &hw ->x[6] , &hw ->x[7] ) ;
  else
    messcrash ("Missing argument --s, expectirn -wws 1:0:2;...    giving the Dynkin lables of the highest weight ") ;
  if (k != sa->rank)
    messcrash ("HW: r5ank = %d but you provided %d Dynkin weights\n", sa->rank, k) ;
  hw->mult = 1 ;

  sprintf (buf, "%d:%d:%d:%d:%d:%d:%d:%d", hw->x[0] , hw ->x[1] , hw ->x[2] , hw ->x[3] , hw ->x[4] , hw ->x[5] , hw ->x[6] , hw ->x[7] ) ;
  dictAdd (sa->dict, buf, &k) ;

  hw->k = k ;
} /* getHighestWeight */

/******************************************************************************************************/

static BOOL demazure (SA *sa, int r1, int *dimp, BOOL show)
{
  BOOL new = FALSE ;
  Array wws = sa->wws ;
  int i, dim, rank = sa->rank ;
  static int pass = 0 ;


  char buf[1024] ;

  if (++pass == 1) new = TRUE ;
    
  for (i = 1 ; i < arrayMax (wws) ; i++) 
    {
      WW *w = arrayp (wws, i, WW) ;
      int jMax = w->x[r1] ;
      int r, n1, n2, k2 ;

      n1 = w->mult ;
      n2 = 0 ;
      
      if (jMax > 0)
	{
	  int j ;
	  WW *w2 ;
	  for (j = 1 ; j <= jMax ; j++)
	    {
	      for (r = 0 ; r < rank ; r++)
		w->x[r] -= array(sa->Cartan, rank * r1 + r, int) * j ;
	      memset (buf, 0, sizeof (buf)) ;
	      sprintf (buf, "%d:%d:%d:%d:%d:%d:%d:%d", w->x[0] , w ->x[1] , w ->x[2] , w ->x[3] , w ->x[4] , w ->x[5] , w ->x[6] , w ->x[7] ) ;
	      dictAdd (sa->dict, buf, &k2) ;
	      w2 = arrayp (wws, k2, WW) ;
	      w = arrayp (wws, i, WW) ; 
	      if (! w2->k)
		{
		  *w2 = *w ;
		  w2->k = k2 ;
		  w2->mult = 0 ;
		}
	      for (r = 0 ; r < rank ; r++)
		w->x[r] += array(sa->Cartan, rank * r1 + r, int) * j ;
	      n2 = w2->mult ;
	    }
	  if (n2 < n1) /* populate */
	    for (j = 1 ; j <= jMax ; j++)
	      {
		new = TRUE ;
		for (r = 0 ; r < rank ; r++)
		  w->x[r] -= array(sa->Cartan, rank * r1 + r, int) * j ;
		memset (buf, 0, sizeof (buf)) ;
		sprintf (buf, "%d:%d:%d:%d:%d:%d:%d:%d", w->x[0] , w ->x[1] , w ->x[2] , w ->x[3] , w ->x[4] , w ->x[5] , w ->x[6] , w ->x[7] ) ;
		dictAdd (sa->dict, buf, &k2) ;
		w2 = arrayp (wws, k2, WW) ;
		w2->mult += n1 - n2 ;
		for (r = 0 ; r < rank ; r++)
		  w->x[r] += array(sa->Cartan, rank * r1 + r, int) * j ;
	    }
	    
	}
    }


  for (dim = 0, i = 1 ; i < arrayMax (wws) ; i++)
    {
      WW *w = arrp (wws, i, WW) ;
      dim += w->mult ;
    }

  if (0)   arraySort (wws, wwOrder) ;
  if (show && new)
    {
      printf ("................Demazure, pass %d, r=%d dim = %d\n", pass, r1, dim) ;
      for (i = 1 ; i < arrayMax (wws) ; i++)
	{
	  int j ;
	  WW *w = arrp (wws, i, WW) ;
	  for (j = rank -1 ; j >= 0 ; j--)
	    printf (" %d", w->x[j]) ;
	  printf ("\tmult=%d k=%d\n", w->mult, w->k) ;
	}
      printf ("\n") ;
    }

  *dimp = dim ;
  
  return new ;
} /* demazure */

/******************************************************************************************************/
/******************************************************************************************************/


int main  (int argc, const char **argv)
{
   AC_HANDLE h = handleCreate () ;
   SA sa ;
   BOOL show = FALSE ;
   freeinit () ;


   memset (&sa, 0, sizeof (sa)) ;
   sa.h = h ;
   sa.type = "toto" ;
   sa.DynkinWeights = "0" ;
   getCmdLineInt (&argc, argv, "-m", &sa.m) ; 
   getCmdLineInt (&argc, argv, "-n", &sa.n) ;
   show = getCmdLineBool (&argc, argv, "-show") ;
   sa.table = getCmdLineBool (&argc, argv, "-table") ;
   getCmdLineOption (&argc, argv, "-type", &sa.type) ;
   getCmdLineOption (&argc, argv, "-w", &sa.DynkinWeights) ;
   if (sa.m < 0) messcrash ("argument -m m of SU(m/n) must be positive or null") ;
   if (sa.n < 0) messcrash ("argument -n n of SU(m/n) must be positive or null") ;
   if (sa.m + sa.n < 2) messcrash ("In SU(m/n) m+n must be >2") ;
   if (sa.m == sa.n) messcrash ("SU(m/m) must be treated separatelly") ; 
   sa.rank = sa.m + sa.n ;

   sa.dict = dictHandleCreate (32, sa.h) ;
   
   getCartan (&sa) ;
    printf ("## Weights of the representations\n") ;

   if (sa.table)
     {
       char *w0 = "0:0:0:0:0:0:0:0:" ;
       w0[2*sa.rank] = 0 ;
       sa.DynkinWeights = w0 ;
       sa.dict = dictHandleCreate (32, sa.h) ;
     }

   else 
     {
       int r = -1, dim = 0 ;
       BOOL ok = TRUE ;

       getHighestWeight (&sa) ;
      while (ok)
	 {
	   ok = FALSE ;
	   for (r = sa.rank - 1 ; r >= 0 ; r--)
	     ok |= demazure (&sa, r, &dim, show) ;
	 }
       printf ("*************************** %s m=%d n=%d %s dim=%d\n ", sa.type, sa.m, sa.n, sa.DynkinWeights, dim) ;
     }
   
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
