#include "ac.h"
#include "matrix.h"

#define MALLOC_CHECK
#define ARRAY_CHECK

/* construct the table of representations of the basic classical Lie superalgebras 
 */

#define RMAX 8 

typedef struct wStruct
{
  int x[RMAX] ;
  BOOL isBlack ;
  int mult ;
  int k ;
  int layer ;
  int l2 ;
  int oddPair ;
  BOOL odd, hw ;
} WW ; /* representation weigth vector */

typedef struct saStruct
{
  AC_HANDLE h ;
  
  int method ;
  BOOL table ;
  const char *type ; /* A,B.C,D,F,G */
  const char *DynkinWeights ; /* 1:0:2:.... */
  BOOL hasY ;
  int YY[RMAX], YYd ;
  int m, n, alpha, rank, rank1, rank2 ;
  Array scale, lowerMetric, upperMetric ;
  Array Cartan, Cartan1, Cartan2 ;
  int detgg ;
  int detGG ;
  int hasAtypic ;
  BOOL hasOdd ;
  BOOL isBlack ;
  int hasExtended ;
  int nEven, nOdd ;
  BOOL odd[RMAX] ;
  BOOL odd1[RMAX] ;
  BOOL odd2[RMAX] ;
  BOOL extended[RMAX] ;
  Array Kac ; /* Kac crystal */
  WW evenHw, evenHw1, evenHwD2, evenHw2 ;
  WW oddHw ;
  WW extendedHw ;
  Array evenRoots ; /* odd roots */
  Array negativeOddRoots ; /* odd roots */
  Array wws ; /* weights */
  WW hw, rho, rho0, rho1 ;
  DICT *dict ;
  int pass ;
  int D1, D2 ; /* number of even and odd generators of the adjoiunt rep */
  Array dd ;   /* dims of all the submodules */
  Array atypic ;
  Array wwsShifted ;
  MX chi ;
} SA ; /* SuperAlgebra struct */


static int wwScalarProduct (SA *sa, WW *ww1, WW *ww2) ;
static BOOL demazureEvenOdd (SA *sa, Array wws, int r1, BOOL even, BOOL snailTrail, int *dimEvenp, int *dimOddp, BOOL show) ;
static int demazure (SA *sa, int *dimEvenp, int *dimOddp, BOOL nonExtended, BOOL show) ;

/******************************************************************************************************/

static int wwCreationOrder (const void *a, const void *b)
{
  const WW *up = (WW *)a ;
  const WW *vp = (WW *)b ;
  int n = 0 ;

  n = up->k - vp->k ;
  return n ;  
} /* wwCreationOrder */

/******************************************************************************************************/

static int wwLayerOrder (const void *a, const void *b)
{
  const WW *up = (WW *)a ;
  const WW *vp = (WW *)b ;
  int n = 0 ;

  n = up->layer - vp->layer ; if (n) return n ;
  n = up->k - vp->k ;  if (n) return n ;

  n = up->x[0] - vp->x[0] ;  if (n) return -n ;
  n = up->x[1] - vp->x[1] ;  if (n) return -n ;
  n = up->x[2] - vp->x[2] ;  if (n) return -n ;
  n = up->x[3] - vp->x[3] ;  if (n) return -n ;
  n = up->x[4] - vp->x[4] ;  if (n) return -n ;

  return n ;  
} /* wwLayerOrder */

/******************************************************************************************************/

static int locateWeight (SA *sa, WW *w, BOOL create)
{
  int kk = 0 ;
  char buf[1024] ;

  /* locate it to construct the multiplicity */
  memset (buf, 0, sizeof (buf)) ;
  sprintf (buf, "%d:%d:%d:%d:%d:%d:%d:%d", w->x[0] , w->x[1] , w->x[2] , w->x[3] , w->x[4] , w->x[5] , w->x[6] , w->x[7] ) ;
  if (create)
    dictAdd (sa->dict, buf, &kk) ;
  else
    dictFind (sa->dict, buf, &kk) ;
  w->l2 = wwScalarProduct (sa, w, w) ;
  return kk ;
} /* locateWeight */

/******************************************************************************************************/

static void metricRescale (Array metric, int r, int ii, int jj, int scale)
{
  int i, k ;
    for (i = ii ; i <= jj ; i++)
  for (k = 0 ; k < r ; k++)
    {
      int x = arr (metric, r * i + k, int) ;
      arr (metric, r * i + k, int) = scale * x ;
    }
  return ;
} /* metricRescale */

/******************************************************************************************************/
/* merge and complete the Cartan Matrix for a superalgebra */
static void mergeCartan (SA *sa, int r1, int r2, BOOL hasY, BOOL show)
{
  int i, j ;
  int r = r1 + r2 + (hasY ? 1 : 0) ;
  int rr = r * r ;
  int dx = 0 ;
  
  sa->Cartan = arrayHandleCreate (rr, int , sa->h) ;

  array (sa->Cartan, rr-1, int) = 0 ;  

  sa->rank = r ;
  sa->hasY = hasY ;
  sa->rank1 = r1 ;
  sa->rank2 = r2 ;
  if (r1 > 0)
    {
      for (i = 0 ; i < r1 ; i++)
	for (j = 0 ; j < r1 ; j++)
	  {
	    array (sa->Cartan, r * i + j, int) = arr (sa->Cartan1, r1 * i + j, int) ;
	  }
      dx = r1 ;
    }
  if (hasY)
    {
      if (r1) array (sa->Cartan, r * (r1 - 1) + r1, int) = -1 ;
      if (r1) array (sa->Cartan, r * (r1) + r1 - 1, int) = 1 ;
      if (r2) array (sa->Cartan, r * (r1) + r1 + 1, int) = 1 ;
      if (1 && r1 && r2) array (sa->Cartan, r * (r1) + r1 - 1, int) = -1 ;
      if (r2) array (sa->Cartan, r * (r1+1) + r1, int) = -1 ;
      dx++ ;
    }
  if (r2 > 0)
    {
      WW oldhw2 = sa->evenHw2 ;

      if (0 && r1)
	sa->evenHw1.x[r1] = -1 ;
      for (i = 0 ; i < r ; i++)
	sa->evenHw2.x[i] = 0 ;
      for (i = 0 ; i < r2 ; i++)
	{	
	  for (j = 0 ; j < r2 ; j++)
	    {
	      array (sa->Cartan, r * (i+dx) + j + dx, int) = arr (sa->Cartan2, r2 * i + j, int) ;
	    }
	  sa->evenHw2.x[dx+i] = oldhw2.x[i] ;
	}
    }
  if (hasY)
    {
      int i, ro = r1 ;
      int YY[RMAX] ;
      
      memset (YY, 0, sizeof (YY)) ;

      switch ((int)sa->type[0])
	{
	default:
	  messcrash ("Only A and C have a U(1) Y hypercharge") ;
	  break ;
	case 'A':
	case 'D':
	  for (i = 0 ; i < ro ; i++)
	    YY[i] =  -(r2+1)*(i+1) ;
	  for (i = r-1 ; i > ro ; i--)
	    YY[i] =  (r1+1)*(r-i) ;
	  YY[ro]  = (r1+1)*(r2+1) * (r1 <= r2 ? -1 : -1) ;
	  sa->YYd = (r1 == r2 ? 1 : r1 - r2) ;
	  for (i = 0 ; i < r ; i++)
	    sa->YY[i] = YY[i] ;

	  if (show)
	    {
	      printf ("............................................. YYd=%d  YY[]=", sa->YYd) ;
	      for (i = 0 ; i < r ; i++)
		printf (" %d", sa->YY[i]) ;
	      printf ("\n") ;
	    }
	  break ;
	case 'C':
	  break ;
	}
    }
  else
    {
      if (0)
	{
	  if (r1) sa->oddHw.x[r1 - 1] = 1 ;
	  if (r2) sa->oddHw.x[r1] = 1 ;
	}
    }
} /* mergeCartan */

/******************************************************************************************************/
/* prepare the Cartan Matrix for a simple Lie algebra block */
static void getOneCartan (SA *sa, char *type, int r, int Lie, BOOL show)
{
  int i, r2 = r * r ;
  WW hw ;
  Array Cartan = arrayHandleCreate (r2, int, sa->h) ;

    array (Cartan, r2 - 1, int) = 0 ;
  
  memset (&hw, 0, sizeof (hw)) ;
  switch ((int)type[0])
    {
    case 'A':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      
      hw.x[0] = 1 ;
      hw.x[r-1] += 1 ;
      break ;
      
    case 'B':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      if (r>=2) array (Cartan, r*(r-1) + r-2, int) = -2 ;

      if (r>1) hw.x[r-2] = 1 ;
      break ;
      
    case 'C':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      array (Cartan, r*1 + 0, int) = -2 ;

      hw.x[r-1] = 2 ;
      break ;
      
    case 'D':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 1) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i> 0 && i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      array (Cartan, r*0 + 2, int) = -1 ;
      array (Cartan, r*2 + 0, int) = -1 ;

      hw.x[0] = 1 ;
      break ;
      
    case 'E':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 1) array (Cartan, r*i + i-1, int) = -1 ;
	    if (i > 0 && i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      array (Cartan, r*0 + 3, int) = -1 ;
      array (Cartan, r*3 + 0, int) = -1 ;

      hw.x[r-1] = 1 ;
      break ;
      
    case 'F':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      array (Cartan, r*2 + 1, int) = -2 ;

      hw.x[0] = 1 ;
      break ;
      
    case 'G':
      array (Cartan, r*0 + 0, int) = 2 ;
      array (Cartan, r*0 + 1, int) = -1 ;
      array (Cartan, r*1 + 1, int) = 2 ;
      array (Cartan, r*1 + 0, int) = -3 ;
      
      hw.x[0] = 1 ;
      break ;
      
    default:
      messcrash ("The Cartan matrix of a Lie algebra of type %c is not yet programmed, sorry.", type) ;
    }


  switch (Lie)
    {
    case 0:
      sa->Cartan = Cartan ;
      sa->evenHw = hw ;
      break ;
    case 1:
      sa->Cartan1 = Cartan ;
      sa->evenHw1 = hw ;
      break ;
    case 2:
      sa->Cartan2 = Cartan ;
      sa->evenHw2 = hw ;
      break ;
    }
  
  if (show)
    {
      int i, j ;
      printf ("\n### getOneCartan (lower) Matrix type %s rank = %d\n", type, r) ;
      for (i = 0 ; i < r ; i++)
	{
	  for (j = 0 ; j < r ; j++)
	    printf ("\t%d", arr (Cartan, r*i + j, int)) ;
	  printf ("\n") ;
	}
    }
  
} /* getOneCartan */

/******************************************************************************************************/

static void getCartan (SA *sa, BOOL show)
{
  int m = sa->m, n = sa->n, r = 0, r1 = 0, r2 = 0, ro = 0 ;
  
  switch ((int)sa->type[0])
    {
    case 'A':
      sa->rank = r = m + n - 1 ;
      r1 = m - 1 ; r2 = n - 1 ;
      if (m <= 0 || n < 0 || m+n == 0)
	messcrash ("Type Lie A(m,n) bad m=%d n=%d", m,n) ;
      if (n == 0)
	getOneCartan (sa, "A", r, 0, TRUE) ;
      else
	{
	  sa->hasOdd = TRUE ;
	  sa->odd[m-1] = 1 ;
	  sa->odd1[m-1] = TRUE ;
	  if (sa->method == 3)
	    {
	      if (n >= 2) sa->odd1[m] = TRUE ;
	      else if (n>= 1 && m >= 2) sa->odd2[m-2] = 1 ;
	    }
	  ro = r1 = m - 1 ;
	  if (m >= 2) getOneCartan (sa, "A", m-1, 1, TRUE) ;
	  if (n >= 2) getOneCartan (sa, "A", n-1, 2, TRUE) ;
	  mergeCartan (sa, m-1, n-1, TRUE, show) ;
	  
	  if (r1) sa->evenHw1.x[ro] = 1 ;
	  if (r1 && r2) sa->evenHw1.x[ro] = -1 ;
	  if (r2) sa->evenHw2.x[ro] = 1 ;
	  if (r1) sa->oddHw.x[ro - 1] = 1 ;
	  if (r2) sa->oddHw.x[ro + 1] = 1 ;
	}
      break ;
      
    case 'B':
      sa->rank = r = m ;
      r1 = r2 = 0 ;
      if (m <= 1)
	messcrash ("Type Lie B(m) and Kac OSp(1/2m) cannot have m=%d < 2\n", m) ;
      if (n == 0)
	getOneCartan (sa, "B", r, 0, TRUE) ;
      else
	{
	  sa->n = 0 ;
	  sa->extended[0] = TRUE ;
	  sa->hasOdd = TRUE ;
	    sa->isBlack = TRUE ;
	  getOneCartan (sa, "C", m, 0, TRUE) ;
	  
	  sa->odd[0] = 1 ;
	  sa->oddHw.x[m-1] = 1 ;
	  sa->oddHw.x[0] = -1 ;
	}
      break ;
      
    case 'C':
      sa->rank = r = m ;
      if (m < 1)
	messcrash ("Type Lie C(m) and Kac OSp(2/2m) cannot have m=%d < 1\n", m) ;
      if (n > 0)
	messcrash ("Type Kac C(n) = OSp(2/2n) is called in this program D(1/n)\n") ;
      if (n == 0)
	getOneCartan (sa, "C", r, 0, TRUE) ;
      break ;

    case 'D':   /* D(m/n) = OSp(2m/2n):SO(2m)+Sp(2n):D(m)+C(n) */
      sa->rank = r = m + n ;
      r1 = m ; r2 = n ;
      if (m < 0 || n < 0 || m+n < 1)
	messcrash ("Type Lie D(m) and Kac D(m/n)=OSp(2m/2n) cannot have m+n=%d+%d < 1\n", m,n) ;
      else if (n == 0)
	{
	  if (m == 1)
	    getOneCartan (sa, "A", m, 0, TRUE) ;
	  else if (m == 2)
	    {
	      sa->hasOdd = TRUE ;
	      getOneCartan (sa, "D", 1, 1, TRUE) ;
	      getOneCartan (sa, "D", 1, 2, TRUE) ;
	      mergeCartan  (sa, 1, 1, FALSE, show) ;
	    }
	  else
	    getOneCartan (sa, "D", m, 0, TRUE) ;
	}
      else if (m == 0)
	getOneCartan (sa, "C", n, 0, TRUE) ;
      else if (m > 1 && n > 0)  /* D(m/n) excluding D(1/m0== C(m) */
	{
	  if (m == 2)
	    {
	      getOneCartan (sa, "D", 1, 1, TRUE) ;
	      getOneCartan (sa, "D", 1, 2, TRUE) ;
	      mergeCartan  (sa, 1, 1, FALSE, show) ;
	      sa->Cartan1 = sa->Cartan ;
	      getOneCartan (sa, "C", n, 2, TRUE) ;
	      mergeCartan  (sa, 2, n, FALSE, show) ;
	      sa->extended[2] = TRUE ;
	      if (sa->method == 3)
		{
		  sa->odd1[2] = TRUE ;
		  sa->odd2[3] = TRUE ;
		}
	      sa->hasOdd = TRUE ;
	      sa->evenHw1.x[0] = 2 ;
	      sa->evenHwD2.x[1] = 2 ;
	      sa->oddHw.x[0] = 1 ;
	      sa->oddHw.x[1] = 1 ;
	      sa->oddHw.x[2] = -1 ;
	    }
	}
      else if (m == 1 && n > 0) /* D(1/n) = super-C (n) = OSp(2/2n) */
	{
	  sa->rank = n + 1 ;
	  sa->n = 0 ;
	  sa->hasOdd = TRUE ;
	  sa->hasY = TRUE ;
	  getOneCartan (sa, "C", n, 1, TRUE) ; 
	  mergeCartan (sa, n, 0, TRUE, show) ;
	  sa->oddHw.x[n-1] = 1 ;
	  if (sa->method == 3)
	    {
	      sa->odd1[0] = TRUE ;
	      sa->odd2[1] = TRUE ;
	    }
	}
      
      break ;
      
    case 'E':
      sa->rank = r = m ;
      switch (m)
	{
    case 6:  /* Lie algebra E6 */
	   case 7:     /* Lie algebra E7 */
	case 8:     /* Lie algebra E8 */
	  getOneCartan (sa, "E", r, 0, TRUE) ;
	  break ;
	default:
	  messcrash ("Type Lie E(m) should be m=6,7,8  m=%d n=%d", m,n) ;
	  break ;
	}
      break ;
      
     case 'F':
      sa->rank = r = m ;
      switch (m + n)
	{
	default:
	  messcrash ("Type Lie G(2) or Kac G(3): m should be 2 or 3 m=%d n=%d", m,n) ;
	  break ;
	case 4:     /* Lie algebra G2 */
	  getOneCartan (sa, "F", 4, 0, TRUE) ;
	  break ;
	case 5:     /* Lie F(4) */
	  sa->rank = sa->m = r = 4 ; sa->n = 0 ;
	  r1 = 3 ; r2 = 1 ;
	  sa->extended[3] = TRUE ;
	  sa->hasOdd = TRUE ;
	  sa->oddHw.x[2] = 1 ;
	  sa->oddHw.x[3] = -1 ;
	  getOneCartan (sa, "B", 3, 1, TRUE) ;
	  getOneCartan (sa, "A", 1, 2, TRUE) ;
	  mergeCartan  (sa, r1, r2, FALSE, show) ;
	  break ;
	}
      break ;
      
     case 'G':
      sa->rank = r = m + n ;
      if (m != 2 || n*(n-1) != 0)
	messcrash ("Type Lie G(2) or Kac G(3): m should be 2 or 3 m=%d n=%d", m,n) ;
      switch (m+n)
	{
	case 2:     /* Lie algebra G2 */
	  getOneCartan (sa, "G", 2, 0, TRUE) ;
	  break ;
	case 3:     /* Lie superalgebra G(3) */
	  r1 = 2 ; r2 = 1 ;
	  sa->extended[2] = TRUE ;
	  sa->hasOdd = TRUE ;
	  sa->odd[1] = 1 ;
	  if (sa->method == 3)
	    {
	      sa->odd1[1] = TRUE ;
	      sa->odd2[2] = TRUE ;
	    }
	  sa->oddHw.x[1] = 1 ;
	  sa->oddHw.x[2] = -1 ;
	  getOneCartan (sa, "G", 2, 1, TRUE) ;
	  getOneCartan (sa, "A", 1, 2, TRUE) ;
	  mergeCartan  (sa, r1, r2, FALSE, show) ;
	  break ;
	}
      break ;
      
    default:
        messcrash ("The Cartan matrix of a superalgebra of type %c is not yet programmed, sorry.", sa->type) ;
    }
  if (show)
    {
      int i, j ;
      printf ("\n### Cartan (lower) Matrix type %s m=%d n=%d rank = %d\n", sa->type, m, n, r) ;
      for (i = 0 ; i < r ; i++)
	{
	    for (j = 0 ; j < r ; j++)
	      printf ("\t%d", arr (sa->Cartan, r*i + j, int)) ;
	    printf ("\n") ;
	}
    }
}   /* KK */
 /* getCartan */

  /******************************************************************************************************/
/* the lower Cartan metric is the symmetrized form of the Cartan matrix
*  the upper metric is the inverse of the lower metric
* the scale is the metric relative to cartan
*/
static void getMetric (SA *sa, BOOL show)
{
  Array aa ;
  int i, j ;
  int m = sa->m, n = sa->n, r = sa->rank ;

  aa = sa->lowerMetric = arrayHandleCopy (sa->Cartan, sa->h) ;
  sa->upperMetric = arrayHandleCreate (r*r, int , sa->h) ;
  sa->scale = arrayHandleCreate (r, int , sa->h) ;
  array (sa->scale, r - 1 , int) = 1 ;

  for (i = 1 ; i < r ; i++)
    {
      int x = - arr (aa, r * (i-1) + i , int) ;
      int y = - arr (aa, r * (i+0) + i - 1, int) ;

      if (x * y < 0 )
	{
	  y = -y ;
	  metricRescale (aa, r, i, r - 1, -1) ;
	}
      if (x > y)
	  metricRescale (aa, r, i, r - 1, x/y) ;
      else if (x < y)
	metricRescale (aa, r, 0, i - 1, y/x) ;
      }

  for (i = 0 ; i < r ; i++)
    {
      for (j = 0 ; j < r ; j++)
	{
	   int x = - arr (aa, r * i + j , int) ;
	   int y = - arr (sa->Cartan, r * i + j, int) ;
	   if (x * y != 0)
	     {
	       arr (sa->scale, i , int) = x / y ;
	       break ;
	     }
	}
    }

    
  sa->detgg = mxIntDeterminant (arrp(sa->lowerMetric, 0, int), r) ;
  if (show && sa->lowerMetric)
      {
	
	printf ("\n### Lower Metric type %s m=%d n=%d rank = %d det=%d\n", sa->type, m, n, r, sa->detgg) ;
	for (i = 0 ; i < r ; i++)
	  {
	    for (j = 0 ; j < r ; j++)
	      printf ("\t%d", arr (sa->lowerMetric, r*i + j, int)) ;
	    printf ("\n") ;
	  }
      }
    sa->detGG = mxIntInverse (arrp (sa->upperMetric, 0, int),  arrp(sa->lowerMetric, 0, int), sa->rank) ;
    if (show && sa->upperMetric)
      {

	printf ("\n### Upper Metric type %s m=%d n=%d rank = %d, det=%d\n", sa->type, m, n, r, sa->detGG) ;
	for (i = 0 ; i < r ; i++)
	{
	  for (j = 0 ; j < r ; j++)
	    printf ("\t%d", arr (sa->upperMetric, r*i + j, int)) ;
	  printf ("\n") ;
	}
      }

    if (show && sa->scale)
      {

	printf ("\n### Scale = metric/Cartan\n") ;
	for (i = 0 ; i < r ; i++)
	  printf ("\t%d", arr (sa->scale, i, int)) ;
	printf ("\n") ;
      }
}  /* getMetric */

/******************************************************************************************************/  
#ifdef JUNK
static int wwNaturalProduct (SA *sa, WW *ww1, WW *ww2)
{
  int i, x = 0, rank = sa->rank ;
  
  for (i = 0 ; i < rank ; i++)
    x += ww1->x[i] * ww2->x[i] ;
  
  return x ;
} /* wwNaturalProduct */
#endif
/******************************************************************************************************/  

static int wwScalarProduct (SA *sa, WW *ww1, WW *ww2)
{
  Array scale = sa->scale ;
  Array upperMetric = sa->upperMetric ;
  int i, j, x = 0, rank = sa->rank ;

  for (i = 0 ; i < rank ; i++)
    for (j = 0 ; j < rank ; j++)
      x += ww1->x[i] * arr (upperMetric, rank*i+j,int) * ww2->x[j] * arr (scale, i, int) * arr (scale, j, int) ;
  
  return x ;
} /* wwScalarProduct */

/******************************************************************************************************/  

static void wwY (SA *sa, WW *ww, int *y1p, int *y2p)
{
  int *YY = sa->YY ;
  int i, y1 = 0, y2 = 0 ;
  
  for (i = 0 ; i < sa->rank ; i++)
    y1 += YY[i] * ww->x[i] ;
  
  y2 = sa->YYd ;
  if (y2 < 0) {y2 = -y2 ; y1 = -y1 ;}
  if (y1 == 0)
    y2 = 1 ;
  if (y1 == y2)
    y1 = y2 = 1 ;
  if (y1 == -y2)
    { y1 = -1 ; y2 = 1 ; }
  for (i = 2 ; i<= y2 ; i++)
    {
      int z1= y1/i, z2 = y2/i ;
      if (i*z1 == y1 && i*z2 == y2)
	{ y1 = z1 ; y2 = z2 ; i = 1 ; }
    }
  *y1p = y1 ;
  *y2p = y2 ;
  return  ;
} /* wwY */

/******************************************************************************************************/  

static void wwsShow (SA *sa, char *title, int type, Array wws, WW *hw)
{
  if (wws)  
    {
      int dimE = 0, dimO = 0 ;
      int ii, iMax = arrayMax (wws) ;
      int r, rank = sa->rank ;

      printf ("\n##################### %s: t=%d : hw", title, type) ;
      for (r = 0 ; hw && r < rank ; r++)
	printf (" %d", hw->x[r]) ;
      printf ("\n") ;

      for (ii = 1 ; ii < iMax ; ii++)
	{
	  WW *ww = arrp (wws, ii, WW) ;
	  if (ww->mult)
	    {
	      for (r = 0 ; r < rank ; r++)
		printf (" %d", ww->x[r]) ;
	      printf ("\tmult=%d k=%d l=%d %s %s Lsquare=", ww->mult, ww->k, ww->layer, ww->odd ? "Odd" : "", ww->hw ? "*" : "") ;
	      ww->l2 = wwScalarProduct (sa, ww, ww) ;
	      if (sa->detGG && ww->l2 % sa->detGG == 0)
		printf ("%d", ww->l2/sa->detGG) ;
	      else
		printf ("%d/%d", ww->l2, sa->detGG) ;
	      if (sa->YYd)
		{
		  int y1, y2 ;
		  wwY (sa, ww, &y1, &y2) ;
		  if (y2 <= 1) printf ("\tY=%d",y1) ;
		  else printf ("\tY=%d/%d",y1,y2) ;
		}
	      printf ("\n") ;
	      if (ww->odd)
		dimO += ww->mult ;
	      else
		dimE += ww->mult ;
	    }
	}
      printf ("## %s dimEven=%d dimOodd=%d dim=%d sdim=%d\n", title, dimE, dimO, dimE+dimO, dimE - dimO) ;
    }
  return ;
} /* wwsShow */

/******************************************************************************************************/

static WW getHighestWeight (SA *sa, int type, BOOL create, BOOL show)
{
  WW *hw ;
  Array wws ;
  int k = 0 ;
  int rank = sa->rank ;
  int r ;

  /* reinitialize a single hw.w */
  sa->pass = 0 ;

  wws = sa->wws = arrayHandleCreate (64, WW, sa->h) ;
  hw = arrayp (wws, 1, WW) ;
  hw[-1].layer = -1000 ;
  
  switch (type)
    {
    case 0: /* trivial h.w. */
      break ;
    case -1: /* use as h.w. the lowering simple root */
      hw->odd = TRUE ;
      for (r = 0 ; r < rank ; r++)
	hw->x[r] = sa->oddHw.x[r] ;
      break ;
    case -2: /* construct from the parameters */
      if (sa->wws)
	k = sscanf (sa->DynkinWeights, "%d:%d:%d:%d:%d:%d:%d:%d", &hw->x[0] , &hw ->x[1] , &hw ->x[2] , &hw ->x[3] , &hw ->x[4] , &hw ->x[5] , &hw ->x[6] , &hw ->x[7] ) ;
      else
	messcrash ("Missing argument --s, expectirn -wws 1:0:2;...    giving the Dynkin lables of the highest weight ") ;

      if (k != sa->rank)
	messcrash ("HW: rank = %d but you provided %d Dynkin weights\n", sa->rank, k) ;
      break ;
    case -3:
      *hw = sa->evenHw ;
      if (sa->rank1)
	*hw = sa->evenHw1 ;
      if (sa->rank2)
	{
	  WW *ew = arrayp (wws, 2, WW) ;
	  if (sa->type[0]=='D' && sa->m == 2)
	    {
	      *ew = sa->evenHwD2 ;
	      ew->mult = create ? 1 : 0 ;
	      ew->hw = TRUE ;
	      ew->k = locateWeight (sa, ew, TRUE) ;

	      ew = arrayp (wws, 3, WW) ;
	    }
	  *ew = sa->evenHw2 ;
	  ew->mult = create ? 1 : 0 ;
	  ew->hw = TRUE ;
	  ew->k = locateWeight (sa, ew, TRUE) ;
	}
      break ;
    default:
      if (sa->negativeOddRoots && type <= arrayMax (sa->negativeOddRoots))
	{
	  WW *wo = arrp (sa->negativeOddRoots, type, WW) ;
	  for (r = 0 ; r < rank ; r++)
	    hw->x[r] = sa->hw.x[r] - wo->x[r] ;
	}
      else
	messcrash ("getHighestWeight received type=%d out of range of odd roots", type) ;
      break ;
    }
  hw->mult = create ? 1 : 0 ;
  hw->hw = TRUE ;
  hw->k = locateWeight (sa, hw, TRUE) ;

  if (type == -2)
    {
      sa->hw = *hw ;
      if (show)
	{
	  printf ("####### Highest weight ") ;
	  for (r = 0 ; r < rank ; r++)
	    printf (" %d", hw->x[r]) ;
	  printf ("\n") ;
	}
    }
  return *hw ;
} /* getHighestWeight */

/******************************************************************************************************/
/* construct the KacCrystal */
static void getKacCrystal (SA *sa, BOOL show)
{
  Array wws ;
  Array negativeOddRoots = sa->negativeOddRoots ;
  int r, rank = sa->rank ;
  int k, kMax = sa->nOdd ;
  int ii, iiMax = 1 ;
  WW *w, hw ;

  /* reinitialize on the trivial h.w */
  sa->dict = dictHandleCreate (32, sa->h) ;
  getHighestWeight (sa, 0, 0, show) ;
  wws = sa->wws ;
  hw = arr (wws, 1, WW) ;

  for (k = 0 ; k < kMax ; k++)
    iiMax *= 2 ;   /* 2^negativeOddRoots = size of the KacCrystal */
  for (ii = 0 ; ii < iiMax ; ii++) /* all points in the Kac Crystal, h.w. already included */
    {
      int layer = 0 ;
      int x = ii ;
      BOOL ok = TRUE ;
      vTXT txt = vtxtHandleCreate (0) ;
      WW w0 ;

      memset (&w0, 0, sizeof (w0)) ;
      w0.hw = 1 ;

      for (k = 0 ; k < kMax ; k++)
	{
	  int yes  = x % 2 ; /* odd root k is used in ii */
	  x /= 2 ;
	  if (yes == 1)
	    {
	      WW *wodd = arrayp (negativeOddRoots, k + 1, WW) ;
	      wodd->l2 = wwScalarProduct (sa, wodd, wodd) ;

	      if (arr(sa->atypic,k+1,int) >= 0 && (sa->isBlack || wodd->l2 == 0))
		{
		  vtxtPrintf (txt, " %d", k+1) ;
		  layer++ ;
		  for (r = 0 ; r < rank ; r++)
		    w0.x[r] += wodd->x[r] ;
		}
	      else
		ok = FALSE ;
	    }
	}

      if (0 && !sa->hasY && ok && ii && sa->rank2)
	{
	  ok = FALSE ;
	  for (r = 0 ; r < rank ; r++)
	    if (!sa->odd[r] && ! sa->extended[r] && w0.x[r] != 0)
	      ok = TRUE ;
	  /*
	    if (! ok)
	    w0.isBlack = TRUE ;
	    ok = TRUE ;
	  */
	}
      for (r = 0 ; ok && r < rank ; r++)
	if (0 && !sa->odd[r] && ! sa->extended[r] && w0.x[r] < 0)
	  ok = FALSE ;
      if (ok)
	{
	  int k2 = locateWeight (sa, &w0, TRUE) ;
	  w = arrayp (wws, k2, WW) ; /* new node */
	  w->k = locateWeight (sa, &w0, TRUE)  ;

	  w->hw = TRUE  ;
	  w->layer = layer ;
	  w->odd = layer %2 ;
	  w->mult++ ;
	  for (r = 0 ; ok && r < rank ; r++)
	    w->x[r] = w0.x[r] ;
	  if (show)
	    {
	      printf (".....Kac crystal: ") ;
	      for (r = 0 ; r < rank ; r++)
		printf (" %d", w->x[r]) ;
	      printf ("  :: ii=%d  yes:%s :: ", ii, vtxtPtr (txt)) ; 
	      printf ("\tmult=%d k=%d l=%d %s %s\n", w->mult, w->k, w->layer, w->odd ? "Odd" : "", w->hw ? "*" : "" ) ;
	    }
	}
      ac_free (txt) ;
    }
  arraySort (sa->wws, wwLayerOrder) ;
  sa->Kac = sa->wws ;
  sa->wws = 0 ;
  if (show)
    wwsShow (sa, "Kac crystal dim =", 1+arrayMax(sa->Kac), sa->Kac, &hw) ;
  if (0)
    printf ("# Constructed %d Kac weigthst\n", arrayMax (sa->Kac) - 1) ;
} /* getKacCrystal */

/******************************************************************************************************/

static Array getShiftedCrystal (SA *sa, int atypic, WW *hwAp, BOOL show)
{
  int rank = sa->rank ;
  int ii, kk, iMax ;
  Array wws ;
  Array Kac = sa->Kac ;
  int kMax = arrayMax (Kac) ;
  WW hwA ;

  hwA = getHighestWeight (sa, atypic, TRUE, 0) ;
  *hwAp = hwA ;
  /* reinitialize wws and but not the dictionary */
  wws = sa->wws = arrayHandleCreate (256, WW, sa->h) ;
  
  for (kk = 1 ; kk < kMax ; kk++)
    {
      WW *wo = arrp (Kac, kk, WW) ;
      WW w, *ww ;
      int r, k2 ;
	    
      if (! wo->mult) continue ;
      /* set the corrds of the next point of the tensor product */
      memset (&w, 0, sizeof (w)) ;
      w = hwA ;
      for (r = 0 ; r < rank ; r++)
	w.x[r] += wo->x[r] ;
      
      /* locate it to construct the multiplicity */
      k2 = locateWeight (sa, &w, TRUE)  ;
      ww = arrayp (wws, k2, WW) ; /* new node */
      ww->k = k2 ;
      ww->layer = wo->layer ;
      ww->odd = wo->odd ;
      ww->mult += hwA.mult * wo->mult ;
      for (r = 0 ; r < rank ; r++)
	ww->x[r] = w.x[r] ;
    }
  /* remove the non-extended  negative even weights */
  iMax = arrayMax (wws) ;
  if (0) wwsShow (sa, "Shifted Kac Crystal", atypic, wws, &hwA) ;
  for (ii = 1 ; ii < iMax ; ii++)
    {
      int r ;
      WW *ww = arrp (wws, ii, WW) ;
      if (! ww->mult) continue ;
      for (r = 0 ; r < rank ; r++)
	if (! sa->odd[r] && /* ! sa->extended[r] && */ ww->x[r] < 0)
	  ww->mult = 0 ;
    }

  /* symmetrize the extended  even weights */
  iMax = arrayMax (wws) ;
  if (0)
    {
      int r ;
      for (r = 0 ; r < rank ; r++)
	if (sa->extended[r])
	  demazureEvenOdd (sa, wws, r, TRUE, TRUE, 0, 0, 0) ; 
    }

  if (show)
    wwsShow (sa, "Shifted Kac Crystal before extended Demazure", atypic, wws, &hwA) ;

  if (1)
    for (ii = 1 ; ii < iMax ; ii++)
      {
	int r, rank = sa->rank ;
	WW *ww = arrp (wws, ii, WW) ;
	if (! ww->mult) continue ;
	for (r = 0 ; r < rank ; r++)
	  if (sa->extended[r] && ww->x[r] > 0)
	    {
	      int r2, k, dk = ww->x[r] ;
	      WW w2 ;
	      int k2, mult2 = 0, dm = 0 ;

	      memset (&w2, 0, sizeof (w2)) ;
	      for (r2 = 0 ; r2 < rank ; r2++)
		w2.x[r2] = ww->x[r2] ;
	      w2.x[r] -= 2 * dk ;
	      k2 = locateWeight (sa, &w2, TRUE) ;
	      mult2 = w2.mult ;
	      dm = ww->mult - mult2 ;
	      if (dm > 0)
		for (k = 1 ; k <= dk ; k++)
		  {
		    WW *ww2 ;

		    memset (&w2, 0, sizeof (w2)) ;
		    for (r2 = 0 ; r2 < rank ; r2++)
		      w2.x[r2] = ww->x[r2] ;

		    w2.odd = ww->odd ;
		    w2.x[r] -= 2 * k ;
		    k2 = locateWeight (sa, &w2, TRUE) ;
		    w2.k = k2 ;
		    ww2  = arrayp (wws, k2, WW) ;
		    if (!ww2->k)
		      *ww2 = w2 ;
		    if (ww2->x[r] < 0)
		      ww2->mult += dm ;
		  }
	    }
      }

  if (show)
    wwsShow (sa, "Shifted Kac Crystal", atypic, wws, &hwA) ;
  return wws ;
} /* getShiftedCrystal */

/******************************************************************************************************/

static void intersectCrystals (SA *sa, Array wws, Array xxs)
{
  WW *ww, *xx ;
  int j ;
  
  wwsShow (sa, "INTERSECT wws", 0, wws, 0) ;
  wwsShow (sa, "INTERSECT xxs", 1, xxs, 0) ;
  for (j = 0 ; j < arrayMax (wws) ; j++)
    { /* we use a single dict */
      ww = arrp (wws, j, WW) ;
      xx = arrayp (xxs, j, WW) ;
      if (xx->mult < ww->mult)
	ww->mult = xx->mult ;
    }
  wwsShow (sa, "INTERSECT SUBTRACTED wws", 0, wws, 0) ;
  return ;
} /* intersectCrystals */

/******************************************************************************************************/

static void getHwCrystal (SA *sa, int *dimp, int *sdimp,  BOOL show)
{
  WW hwA ;
  Array wws = getShiftedCrystal (sa, -2, &hwA, TRUE) ;
  Array old = sa->wws ;

  if (sa->hasAtypic)
    {
      int ii ;
      
      for (ii = 0 ; 1 && ii < arrayMax (sa->atypic) ; ii++)
	if (array (sa->atypic, ii, int))
	  {
	    Array xxs = getShiftedCrystal (sa, ii, &hwA, show) ;
	    intersectCrystals (sa, old, xxs) ;
	    ac_free (xxs) ;
	    
	    if (show)
	      wwsShow (sa, "Get Subtracted Atypical Tensor Product", ii, old, &hwA) ;
	  }
      wws = old ;
    }
  sa->wws = wws ;
  if (dimp)
    {
      int dimE = 0, dimO = 0 ;
      WW *ww ;
      int ii ;
      for (ii = 0 ; ii < arrayMax (wws) ; ii++)
	{
	  ww = arrp (wws, ii, WW) ;
	  if (ww->mult > 0)
	    {
	      if (ww->odd)
		dimO += ww->mult ;
	      else
		dimE += ww->mult ; 
	    }
	}
      *dimp = dimE + dimO ;
      *sdimp = dimE - dimO ;
    }
  arraySort (wws, wwLayerOrder) ;
  if (1)
    {
      Array old = wws ;
      WW *ww, *ww1 ;
      int ii ;

      wws = sa->wws = arrayHandleCreate (64, WW, sa->h) ;
      /* clean up the dict and redeclare the locators */
      sa->dict = dictHandleCreate (32, sa->h) ;
      for (ii = 0 ; ii < arrayMax (old) ; ii++)
	{
	  ww = arrp (old, ii, WW) ;
	  if (ww->mult > 0)
	    {
	      int k = locateWeight (sa, ww, TRUE) ;
	      ww1 = arrayp (wws, k, WW) ;
	      *ww1 = *ww ;
	      ww1->k = k ;
	    }
	}
    }
  if (show)
    wwsShow (sa, "Final Atypical Tensor Product", 999, wws, 0) ;

  return ;
} /* getHwCrystal */

/******************************************************************************************************/

static BOOL demazureEvenOdd (SA *sa, Array wws, int r1, BOOL even, BOOL snailTrail, int *dimEvenp, int *dimOddp, BOOL show) 
{
  BOOL new = FALSE ;
  BOOL odd = ! even ;
  int i, dimEven, dimOdd, rank = sa->rank ;
  
  if (sa->pass++ == 0) new = TRUE ;
    
  for (i = 1 ; i < arrayMax (wws) ; i++) 
    {
      WW *w = arrayp (wws, i, WW) ;
      int n1 = w->mult ;
      int n2 = 0 ;
      int na = 0 ;
      int dn = 0 ;
      int jMax = w->x[r1] ; /* default even values */
      if (! even)
	jMax = 1 ;
      if (! even)
	n1 -= w->oddPair ;
      if (n1 <= 0)
	continue ;
      if (jMax > 0)  /* create or check existence of the new weights along the sl(2) submodule */
	{
	  int j, r, k2 ;
	  int dj = snailTrail ? 1 : jMax - 1 ;
	  WW *w2 ;
	  
	  for (j = 1 ; j <= jMax ; j += dj)
	    {
	      /* position to the new weight */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] -= array(sa->Cartan, rank * r + r1, int) * j ;
	      
	      /* locate the k2 index of the WW weight structure of the corresponding member of the multiplet */
	      k2 = locateWeight (sa, w, TRUE)  ;
	      
	      /* create w2(k2) it if needed */
	      w2 = arrayp (wws, k2, WW) ;
	      w = arrayp (wws, i, WW) ;  /* needed because the wws array may be relocated in RAM upon extension */ 
	      if (! w2->k)
		{
		  *w2 = *w ;
		  w2->k = k2 ;
		  w2->mult = 0 ;
		  w2->oddPair = 0 ;
		  w2->hw = FALSE ;
		  w2->l2 = wwScalarProduct (sa, w2, w2) ;
		}
	      n2 = na + w2->mult ;

	      /* reposition w to its original location */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] += array(sa->Cartan, rank * r + r1, int) * j ;
	    }
	}

      dn = n1 - n2 ; /* defect multiplicity */
      if (jMax > 0 && dn > 0)   /* populate the new multiplet */
	{
	  int j, r, k2 ;
	  int dj = snailTrail ? 1 : jMax - 1 ;
	  WW *w2 ;
	  
	  new = TRUE ;
	  for (j = 1 ; j <= jMax ; j += dj)
	    {
	      /* position to the new weight */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] -= array(sa->Cartan, rank * r + r1, int) * j ;
	      
	      /* locate the k2 index of the WW weight structure of the corresponding member of the multiplet */
	      k2 = locateWeight (sa, w, TRUE) ;
	      w2 = arrayp (wws, k2, WW) ;
	      
	      /* increase multiplicity of the new multiplet */
	      if (1)
		{
		  if (! w2->hw)
		    w2->mult += dn ;
		}
	      else /* classic demazure , correct for Lie algebras */
		w2->mult += dn ;

	      if (odd)
		{
		  w->oddPair += dn ;
		  w2->oddPair += dn ;
		  w2->odd = ! w->odd;
		  w2->layer = w->layer + 1 ;
		}
	      /* reposition w to its original location */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] += array(sa->Cartan, rank * r + r1, int) * j ;
	    }
	}
      
    }

  for (dimEven = dimOdd = 0, i = 1 ; i < arrayMax (wws) ; i++)
    {
      WW *ww = arrp (wws, i, WW) ;
      if (ww->odd)
	dimOdd += ww->mult ;
      else 
	dimEven += ww->mult ;
    }

  if (1)   arraySort (wws, wwLayerOrder) ;
  if (0 && new)
    {
      printf ("................Demazure, pass %d, r=%d dimEven = %d dimOdd = %d\n", sa->pass, r1, dimEven, dimOdd) ;
      for (i = 1 ; i < arrayMax (wws) ; i++)
	{
	  int j ;
	  WW *w = arrp (wws, i, WW) ;
	  for (j = 0 ; j < rank ; j++)
	    printf (" %d", w->x[j]) ;
	  printf ("\tmult=%d k=%d l=%d %s\n", w->mult, w->k, w->layer, w->hw ? "*" : "" ) ;
	}
      printf ("\n") ;
    }
  if (1)   arraySort (wws, wwCreationOrder) ;
  
  if (dimEvenp) *dimEvenp = dimEven ;
  if (dimOddp) *dimOddp = dimOdd ; 
  
  return new ;
} /* demazureEvenOdd */

/******************************************************************************************************/

static int demazure (SA *sa, int *dimEvenp, int *dimOddp, BOOL nonExtended, BOOL show)
{
  BOOL ok = TRUE ;
  BOOL debug = FALSE ;
  int r, dimEven = 0, dimOdd = 0 ;

  if (debug)
    wwsShow (sa, "Before Demazure", 0, sa->wws, 0) ;

  while (ok)
    {
      ok = FALSE ;
      for (r = 0 ; r < sa->rank ; r++)
	{
	  if (! sa->odd[r] && ! sa->odd1[r] && ! sa->odd2[r] && ! (nonExtended && sa->extended[r]))
	    ok |= demazureEvenOdd (sa, sa->wws, r, TRUE, TRUE, &dimEven, &dimOdd, show) ;
	  if (ok && debug)
	    wwsShow (sa, "Inside Demazure", 0, sa->wws, 0) ;
	}
      if (sa->hasOdd)
	{
	  switch (sa->method)
	    {
	    case 2:
	      for (r = 0 ; r < sa->rank ; r++)
		{
		  if (sa->odd[r])
		    ok |= demazureEvenOdd (sa, sa->wws, r, FALSE, sa->hasY, &dimEven, &dimOdd, show) ;
		  if (ok && debug)
		    wwsShow (sa, "Inside Demazure Odd", 0, sa->wws, 0) ;
		}
	      break ;
	    case 3:
	      for (r = 0 ; r < sa->rank ; r++)
		{
		  if (sa->odd[r])
		    ok |= demazureEvenOdd (sa, sa->wws, r, FALSE, sa->hasY, &dimEven, &dimOdd, show) ;
		  if (ok && debug)
		    wwsShow (sa, "Inside Demazure Odd", 0, sa->wws, 0) ;
		}
	      break ;
	    }
	}
      if (0 && sa->hasExtended)
	{
	  ok = TRUE ;
	  while (ok)
	    {
	      ok = FALSE ;
	      for (r = 0 ; r < sa->rank ; r++)
		if (sa->extended[r])
		  ok |= demazureEvenOdd (sa, sa->wws, r, TRUE, TRUE, &dimEven, &dimOdd, show) ;
	    }
	}
    }
  
  if (dimEvenp) *dimEvenp = dimEven ;  
  if (dimOddp) *dimOddp = dimOdd ; 

  if (show)
    wwsShow (sa, "Demazure", 0, sa->wws, 0) ;
  
  return dimEven ;
} /* demazure */

/******************************************************************************************************/

static void getNegativeOddRoots (SA *sa, BOOL show)
{
  int dimEven = 0 ;
  int dimOdd = 0 ;
  
  sa->dict = dictHandleCreate (32, sa->h) ;
  
  /* use as h.w. the lowering simple root */
  getHighestWeight (sa, -1, 1, 0) ;
  /* construct the first layer using Demazure */
  demazure (sa, &dimEven, &dimOdd, TRUE, show) ;
  sa->nOdd = dimOdd ; 
  sa->negativeOddRoots = sa->wws ;
  sa->wws = 0 ;
  
  if (show)
    wwsShow (sa, "Negative Odd roots", 1, sa->negativeOddRoots, arrp (sa->negativeOddRoots, 1, WW)) ;
  printf ("## Constructed %d odd roots\n", dimOdd) ;

  return ;
} /* getNegativeOddRoots */

/******************************************************************************************************/


static void getAdjoint (SA *sa, BOOL show)
{
  int dimE = 0 ;
  int dimOdd = 0 ;

    sa->dict = dictHandleCreate (32, sa->h) ;

  /* use as h.w. the lowering simple root */
  getHighestWeight (sa, -3, 1, 0) ;
  /* construct the adjoint layer using Demazure */
  sa->nEven = demazure (sa, &dimE, &dimOdd, FALSE, show) ;

  sa->evenRoots = sa->wws ;
  sa->wws = 0 ;

  if (show)
    wwsShow (sa, "Even Adjoint ", 1, sa->evenRoots, 0) ;
  printf ("# Constructed %d adjoint even roots\n", dimE) ;

  return ;
} /* getNegativeOddRoots */

/******************************************************************************************************/

static void getRho (SA *sa, BOOL show)
{
  WW *ww ;
  int ii, r ;
  int rank = sa->rank ;

  memset (&sa->rho0, 0, sizeof (WW)) ;
  /* rho0: sum of the positive even roots */
  if (sa->evenRoots)
    for (ii = 0 ; ii < arrayMax (sa->evenRoots) ; ii++)
      {
	int x, i, j ;
	BOOL isPositive = TRUE ;
	ww = arrp (sa->evenRoots, ii, WW) ;
	for (i = 0 ; i < r ; i++)
	  {
	    x = 0 ;
	    if (! sa->odd[i])
	      for (j = 0 ; j < r ; j++)
		x +=  arr (sa->upperMetric, rank * i + j, int) * arr (sa->scale, 0, int) * ww->x[j] ;
	    if (x * arr (sa->scale, i, int) < 0)
	      isPositive = FALSE ;
	  }
	if (isPositive)
	  for (r = 0 ; r < rank ; r++)
	    sa->rho0.x[r] += ww->x[r] ; /* use + since we deal with the negative odd roots */
      }
  /* rho1: sum of the null positive odd roots */
  memset (&sa->rho1, 0, sizeof (WW)) ;
  if (sa->negativeOddRoots)
    for (ii = 0 ; ii < arrayMax (sa->negativeOddRoots) ; ii++)
      {
	ww = arrp (sa->negativeOddRoots, ii, WW) ;
	for (r = 0 ; r < rank ; r++)
	  sa->rho1.x[r] += ww->x[r] ; /* use + since we deal with the negative odd roots */
      }
  for (r = 0 ; r < rank ; r++)
    sa->rho.x[r] = sa->rho0.x[r]  + sa->rho1.x[r] ;
  
  if (show)
    {
      int i ;
      printf ("#================================= Sum of the even positive roots : 2 rho_0 = ") ;
      for (i = 0 ; i < rank ; i++)
	printf (" %d", sa->rho0.x[i]) ;
      
      if (sa->negativeOddRoots)
	{
	  printf ("\n#================================= Sum of the odd negative roots : 2 rho_1 = ") ;
	  for (i = 0 ; i < rank ; i++)
	    printf (" %d", sa->rho1.x[i]) ;
	}
      printf ("\n") ;
    }
} /* getRho */

/******************************************************************************************************/

static void getAtypic (SA *sa, BOOL show)
{
  WW hwT ;
  int rank = sa->rank ;
  int ii, r ;
  Array oddRoots = sa->negativeOddRoots ;
  
  sa->atypic = arrayHandleCreate (arrayMax (sa->negativeOddRoots), int, sa->h) ;  

  /* use as h.w. the declared h.w. */
  getHighestWeight (sa, -2, 1, TRUE) ;
  hwT = sa->hw ;  /* the highest weight L */ 
  for (r = 0 ; r < rank ; r++)
    hwT.x[r] = 2 * sa->hw.x[r] + sa->rho.x[r] ;     /* hw.w translated 2(L + rho) */
      
  for (ii = 1 ; ii < arrayMax (oddRoots) ; ii++)
    {
      /* x = < L + rho | beta_i > */
      int i, j, x = 0 ;
      WW *ww = arrp (oddRoots, ii, WW) ;

      for (i = 0 ; i < rank ; i++)
	for (j = 0 ; j < rank ; j++)
	  x += hwT.x[i] * arr (sa->upperMetric, rank*i + j, int) * arr (sa->scale, j, int) * ww->x[j] ;
      if (x == 0)
	{
	  array (sa->atypic, ii, int) = 1 ;
	  sa->hasAtypic = TRUE ;
	}
      if (ii == 1 && x != 0 && ! sa->isBlack)
	messcrash ("The trivial representaion is not atypic 1") ;
    }
  
  if (1 || show)
    {
      int i ;
      
      if (sa->hasAtypic)
	{
	  printf ("\n#=================================   Atypical") ;
	  for (i = 0 ; i < arrayMax (sa->atypic) ; i++)
	    if (array (sa->atypic, i, int))
	      printf (" %d", i) ;
	  printf ("\n") ;
	}
      else
	printf ("\n#=================================  Typical\n") ;
    }

  return ;
} /* getAtypic */

/******************************************************************************************************/
/******************************************************************************************************/


int main  (int argc, const char **argv)
{
   AC_HANDLE h = handleCreate () ;
   SA sa ;
   BOOL show = FALSE ;
   BOOL hasOdd ;
   int dimEven = 0, dimOdd = 0 ;
   freeinit () ;


   memset (&sa, 0, sizeof (sa)) ;
   sa.h = h ;
   sa.method = 3 ;   
   sa.type = "toto" ;
   sa.DynkinWeights = "0" ;
   getCmdLineInt (&argc, argv, "-m", &sa.m) ; 
   getCmdLineInt (&argc, argv, "-n", &sa.n) ;
   sa.alpha = 1 ;
   getCmdLineInt (&argc, argv, "-alpha", &sa.alpha) ;
   show = getCmdLineBool (&argc, argv, "-show") ;
   sa.table = getCmdLineBool (&argc, argv, "-table") ;
   getCmdLineOption (&argc, argv, "-type", &sa.type) ;
   getCmdLineOption (&argc, argv, "-w", &sa.DynkinWeights) ;
   if (sa.m < 0) messcrash ("argument -m m of SU(m/n) must be positive or null") ;
   if (sa.n < 0) messcrash ("argument -n n of SU(m/n) must be positive or null") ;
   
   sa.dict = dictHandleCreate (32, sa.h) ;
   
   getCartan (&sa, show) ;
   getMetric (&sa, show) ;
   
   hasOdd = sa.hasOdd ;
   sa.hasOdd = FALSE ;
   if (hasOdd) /* do this first then destroy the dict */
     getNegativeOddRoots (&sa, show) ;
   /* apply Demazure of the even group */
   getAdjoint (&sa, show) ;
   getRho (&sa, show) ;
   sa.hasOdd = hasOdd ;
   
   if (sa.hasOdd)
     {                        /* contruct the kasCrystal */
       getAtypic (&sa, show) ;
       if (sa.method <= 2)
	 getKacCrystal (&sa, show) ; 
     }
   
   /* construct the h.w of the top layer even module */
   sa.dict = dictHandleCreate (32, sa.h) ;
   getHighestWeight (&sa, -2, TRUE, TRUE) ;
   printf ("*************************** %s m=%d n=%d %s \n "
	   , sa.type, sa.m, sa.n, sa.DynkinWeights) ;

   if (!sa.hasOdd)
     demazure (&sa, &dimEven, &dimOdd, TRUE, FALSE) ;
   else
     {
       switch (sa.method)
	 {
	 case 1:
	   getHwCrystal (&sa, &dimEven, &dimOdd, show) ;
	   /* complete the hw Crystal to a full module */
	   if (1) demazure (&sa, &dimEven, &dimOdd, TRUE, FALSE) ;
	   break ;
	 case 2:
	 case 3:
	   demazure (&sa, &dimEven, &dimOdd, TRUE, show) ;
	   break ;
	 }
     }
	   
   printf  ("################################################## Final Representation dimEven/dimOdd %d / %d\n",  dimEven, dimOdd) ;
   arraySort (sa.wws, wwLayerOrder) ;
   if (show) wwsShow (&sa, "Final representation", 1, sa.wws, &sa.hw) ;
   if (sa.hasOdd) getAtypic (&sa, show) ;
   messfree (h) ;
   printf ("A bientot\n") ;
   
   return 0 ;   
} /* main */

/************************************************************************/
/************************************************************************/
/************************************************************************/
/*
Variation:
  i again compute the shifted crystal and set w->mult = 0 is w->mult>x->mult
  this is favorable for F(4) adjoint because it kill the 1:0:0:0 dim 7 at level 0
  but this could be computed by hand, in the exact way
 
   BUG: the h.w. is not recovered by the shifted crystal which is absurd 
run -type D -m 2 -n 1 -alpha 3  -w 0:0:2 -show
*/
