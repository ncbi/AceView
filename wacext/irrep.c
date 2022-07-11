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
  int layer ;
  BOOL odd, hw ;
} WW ; /* representation weigth vector */

typedef struct saStruct
{
  AC_HANDLE h ;
  
  BOOL table ;
  const char *type ; /* A,B.C,D,F,G */
  const char *DynkinWeights ; /* 1:0:2:.... */
  int m, n, rank ;
  Array Cartan ;
  int hasOdd ;
  BOOL odd[RMAX] ;
  Array Kac ; /* Kac crystal */ 
  Array oddRoots ; /* odd roots */
  Array wws ; /* weights */
  DICT *dict ;
  int pass ;
  int D1, D2 ; /* number of even and odd generators of the adjoiunt rep */
  Array dd ;   /* dims of all the submodules */
  MX chi ;
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
	    int o = m - 1 ;
	    sa->odd[o] = TRUE ;
	    array (Cartan, r*o + o , int) = 0 ;
	    if (o > 0) array (Cartan, r*o + o - 1 , int) = -1 ;
	    else array (Cartan, r*o + o + 1, int) = 1 ;
	  }
      }
      break ;
      
    case 'B':
      if (m < 2)
	messcrash ("Type B(m)): m should be at least 2:  m=%d n=%d", m,n) ;
      sa->rank = r = m ;
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
      sa->rank = r = 4 ;
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
      if (m != 2 && m != 3)
	messcrash ("Type Lie G(2) or Kac G(3): m should be 2 or 3 m=%d n=%d", m,n) ;
      sa->rank = r = m ;
      rr = r*r ;
      Cartan = arrayCreate (rr, int) ;
      {
	array (Cartan, r*0 + 0, int) = 2 ;
	array (Cartan, r*0 + 1, int) = -3 ;
	array (Cartan, r*1 + 1, int) = 2 ;
	array (Cartan, r*1 + 0, int) = -1 ;
      }
      if (m == 3)
	{
	  sa->odd[2] = TRUE ;
	  array (Cartan, r*2 + 2 , int) = 0 ;
	  array (Cartan, r*2 + 1, int) = -1 ;
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
	  if (sa->odd[i])
	    sa->hasOdd = i + 1 ;
	  for (j = 0 ; j < r ; j++)
	    printf ("\t%d", arr (Cartan, r*i + j, int)) ;
	  printf ("\n") ;
	}
    }
  sa->Cartan = Cartan ;
}  /* getCartan */

/******************************************************************************************************/
/* construct the KacCrystal */
static void kacCrystal (SA *sa, BOOL show)
{
  Array wws = sa->wws ;
  int jj = arrayMax (wws) ;
  Array oddRoots = sa->oddRoots ;
  int r, rank = sa->rank ;
  int k, kMax = arrayMax (oddRoots) - 1 ;
  int k2, ii, iiMax = 1 ;
  WW *w1, *w ;
  char buf[1024] ;

  if (show)
    printf("\n####### Kac Crystal \n") ;
  for (k = 0 ; k < kMax ; k++)
    iiMax *= 2 ;   /* 2^oddroots = size of the KacCrystal */
  for (ii = 0 ; ii < iiMax ; ii++) /* all points in the Kac Crystal */
    {
      int layer = 0 ;
      int x = ii ;
      BOOL ok = TRUE ;

      w = arrayp (wws, jj++, WW) ; /* new node */
      w1 = arrayp (wws, 1, WW) ; /* the heighest weight */
      *w = *w1 ;
      for (k = 0 ; k < kMax ; k++)
	{
	  int yes  = x % kMax ; /* odd root k is used in ii */
	  x /= kMax ;
	  if (yes)
	    {
	      WW *wodd = arrayp (oddRoots, k + 1, WW) ;
	      layer++ ;
	      for (r = 0 ; r < rank ; r++)
		w->x[r] -= wodd->x[r] ;
	    }
	}
      for (r = 0 ; ok && r < rank ; r++)
	if (w->x[r] < 0)
	  {
	    /* reject */
	    jj-- ;
	    arrayMax (wws) = jj ;
	    ok = FALSE ;
	  }
      if (ok)
	{
	  memset (buf, 0, sizeof (buf)) ;
	  sprintf (buf, "%d:%d:%d:%d:%d:%d:%d:%d", w->x[0] , w ->x[1] , w ->x[2] , w ->x[3] , w ->x[4] , w ->x[5] , w ->x[6] , w ->x[7] ) ;
	  dictAdd (sa->dict, buf, &k2) ;
	  w->k = k2 ;
	  w->hw = TRUE  ;
	  w->layer = layer ;
	  w->odd = layer %2 ;

	  if (show)
	    {
	      for (r = 0 ; r < rank ; r++)
		printf (" %d", w->x[r]) ;
	      printf ("\tmult=%d k=%d l=%d %s %s\n", w->mult, w->k, w->layer, w->odd ? "Odd" : "", w->hw ? "*" : "" ) ;
	    }
	}
    }
} /* kacCrystal */

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

  hw->hw = TRUE ;
  hw->k = k ;
} /* getHighestWeight */

/******************************************************************************************************/

static BOOL demazure (SA *sa, int r1, int *dimp, int *sdimp, BOOL show)
{
  BOOL new = FALSE ;
  Array wws = sa->wws ;
  int i, dim, sdim, rank = sa->rank ;
  BOOL odd = sa->odd[r1] ;
  
  char buf[1024] ;

  if (sa->pass++ == 0) new = TRUE ;
    
  for (i = 1 ; i < arrayMax (wws) ; i++) 
    {
      WW *w = arrayp (wws, i, WW) ;
      int n1 = w->mult ;
      int n2 = 0 ;
      int na = 0 ;
      int dn = 0 ;
      int jMax = w->x[r1] ; /* default even values */
      BOOL oddLayer = w->odd ;

      
      if (odd) /* assess the multiplicity of the ancestor descendent  */
	{
	  int r, k2 = 0 ;
	  int j = 1 ;

	  if (0 && ! w->hw)  /* only consider the Kac crystal */ 
	    continue ;
	  if (0 && i == 1 && w->x[r1] == 0) /* atypical type one */
	    continue ;
	  if (0 && r1 > 0 && i == 2 && w->x[r1] == w->x[r1-1]+ 1) /* atypical type 2 */
	    continue ;
	  if (0 && r1 > 1 && i == 3 && w->x[r1] == w->x[r1-1]+ w->x[r1-2] + 2) /* atypical type 3 */
	    continue ;
	  oddLayer = ! w->odd ;
	  
	  /* position to ancestor */
	  for (r = 0 ; r < rank ; r++)
	    w->x[r] += array(sa->Cartan, rank * r1 + r, int) * j ;
	  memset (buf, 0, sizeof (buf)) ;
	  sprintf (buf, "%d:%d:%d:%d:%d:%d:%d:%d", w->x[0] , w ->x[1] , w ->x[2] , w ->x[3] , w ->x[4] , w ->x[5] , w ->x[6] , w ->x[7] ) ;
	  if (dictFind (sa->dict, buf, &k2))
	    {
	      WW * w2 = arrayp (wws, k2, WW) ;
	      na = w2->mult ;
	    }
	  else
	    na = 0 ; /* no ancestor */
	  if (0 && w->layer && ! na && arr(wws,1,WW).x[r1]== arr(wws,1,WW).x[r1-1] + 0)
	    na = n1 ;
	  if (na >=  n1)  /* status quo */
	    jMax = 0 ;
	  else
	    jMax = 1 ;   /* odd value */

	  /* reposition w to its original location */
	  for (r = 0 ; r < rank ; r++)
	    w->x[r] -= array(sa->Cartan, rank * r1 + r, int) * j ;
	}
	  
      if (jMax > 0)  /* create or check existence of the new weights along the sl(2) submodule */
	{
	  int j, r, k2 ;
	  WW *w2 ;
	  
	  for (j = 1 ; j <= jMax ; j++)
	    {
	      /* position to the new weight */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] -= array(sa->Cartan, rank * r1 + r, int) * j ;
	      
	      /* locate the k2 index of the WW weight structure of the corresponding member of the multiplet */
	      memset (buf, 0, sizeof (buf)) ;
	      sprintf (buf, "%d:%d:%d:%d:%d:%d:%d:%d", w->x[0] , w ->x[1] , w ->x[2] , w ->x[3] , w ->x[4] , w ->x[5] , w ->x[6] , w ->x[7] ) ;
	      dictAdd (sa->dict, buf, &k2) ;
	      
	      /* create w2(k2) it if needed */
	      w2 = arrayp (wws, k2, WW) ;
	      w = arrayp (wws, i, WW) ;  /* needed because the wws aray may be relocated in RAM upon extension */ 
	      if (! w2->k)
		{
		  *w2 = *w ;
		  w2->k = k2 ;
		  w2->mult = 0 ;
		  w2->hw = FALSE ;
		  w2->odd = oddLayer ;
		  if (odd)
		    w2->layer++ ;
		  if (odd)
		    w2->hw = TRUE ;
		}
	      n2 = na + w2->mult ;

	      /* reposition w to its original location */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] += array(sa->Cartan, rank * r1 + r, int) * j ;
	    }
	}

      dn = n1 - n2 ; /* defect multiplicity */
      if (jMax > 0 && dn > 0)   /* populate the new multiplate */
	{
	  int j, r, k2 ;
	  WW *w2 ;
	  
	  new = TRUE ;
	  for (j = 1 ; j <= jMax ; j++)
	    {
	      /* position to the new weight */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] -= array(sa->Cartan, rank * r1 + r, int) * j ;
	      
	      /* locate the k2 index of the WW weight structure of the corresponding member of the multiplet */
	      memset (buf, 0, sizeof (buf)) ;
	      sprintf (buf, "%d:%d:%d:%d:%d:%d:%d:%d", w->x[0] , w ->x[1] , w ->x[2] , w ->x[3] , w ->x[4] , w ->x[5] , w ->x[6] , w ->x[7] ) ;
	      dictFind (sa->dict, buf, &k2) ; /* notice, w2(k2) was created just above */
	      w2 = arrayp (wws, k2, WW) ;
	      
	      /* increase multiplicity of the new multiplet */
	      if (0)
		{
		  if (w2->mult < w->mult + dn  + (w2->hw ? 1 : 0))
		    w2->mult += dn ;
		}
	      else /* classic demazure , correct for Lie algebras */
		w2->mult += dn ;
	      
	      /* reposition w to its original location */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] += array(sa->Cartan, rank * r1 + r, int) * j ;
	    }
	}
      
    }

  for (dim = sdim = 0, i = 1 ; i < arrayMax (wws) ; i++)
    {
      WW *w = arrp (wws, i, WW) ;
      if (w->odd)
	sdim += w->mult ;
      else
	dim += w->mult ;
    }

  if (0)   arraySort (wws, wwOrder) ;
  if (show && new)
    {
      printf ("................Demazure, pass %d, r=%d dim = %d sdim = %d\n", sa->pass, r1, dim + sdim, dim - sdim) ;
      for (i = 1 ; i < arrayMax (wws) ; i++)
	{
	  int j ;
	  WW *w = arrp (wws, i, WW) ;
	  for (j = 0 ; j < rank ; j++)
	    printf (" %d", w->x[j]) ;
	  printf ("\tmult=%d k=%d l=%d %s %s\n", w->mult, w->k, w->layer, w->odd ? "Odd" : "", w->hw ? "*" : "" ) ;
	}
      printf ("\n") ;
    }

  *dimp = dim + sdim ;
  *sdimp = dim - sdim ;
  
  return new ;
} /* demazure */

/******************************************************************************************************/

static void demazureEven (SA *sa, int *dimp, int *sdimp, BOOL show)
{
  BOOL ok = TRUE ;
  int r, dim = 0, sdim = 0 ;

  while (ok)
    {
      ok = FALSE ;
      for (r = 0 ; r < sa->rank ; r++)
	if (! sa->odd[r])
	  ok |= demazure (sa, r, &dim, &sdim, show) ;
    }
  *dimp = dim ; 
  *sdimp = sdim ;
  return ;
} /* demazure */

/******************************************************************************************************/
/******************************************************************************************************/


int main  (int argc, const char **argv)
{
   AC_HANDLE h = handleCreate () ;
   SA sa ;
   BOOL show = FALSE ;
   int dim = 0, sdim = 0 ;

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

   sa.dict = dictHandleCreate (32, sa.h) ;
   
   getCartan (&sa) ;
   getHighestWeight (&sa) ;

   if (sa.hasOdd)
     {
       int rank = sa.rank ;
       int r, r1 = sa.hasOdd - 1 ;

       /* use as h.w. the lowering simple root */
       WW *w = arrp (sa.wws, 1, WW) ;
       for (r = 0 ; r < rank ; r++)
	 w->x[r] = - array(sa.Cartan, rank * r1 + r, int) ;
       /* construct the first layer using Demazure */
       demazureEven (&sa, &dim, &sdim, show) ;
       sa.oddRoots = sa.wws ;
       sa.wws = 0 ;
       /* reinitialize a single hw.w */
       sa.pass = 0 ;
       sa.dict = dictHandleCreate (32, sa.h) ;
       getHighestWeight (&sa) ;
       /* consruct the antisymmetric Kac crystal */
       kacCrystal (&sa, show) ; 
       /* we end up with wws just containing the Kac Crystal */
     }
   
   /* apply Demazure of the even group */
   printf ("## Weights of the representations\n") ;

   demazureEven (&sa, &dim, &sdim, show) ;
   printf ("*************************** %s m=%d n=%d %s dim=%d sdim = %d\n "
	       , sa.type, sa.m, sa.n, sa.DynkinWeights, dim, sdim) ;
  messfree (h) ;
  printf ("A bientot\n") ;

  return 0 ;
}

/************************************************************************
SU(3/2) shifted trivial, layer 2, should be 6S * 1A   + 3'A * 3S  == 9+6 = 15
................Demazure, pass 17, r=2 dim = 12 sdim = 0

layer 1 and first contact to layer 2 viy the odd root
 0 0 10 0	mult=1 k=1 l=0  *

 0 1 10 1	mult=1 k=2 l=1 Odd *
 1 -1 11 1	mult=1 k=3 l=1 Odd *
 -1 0 11 1	mult=1 k=6 l=1 Odd *

 0 1 11 -1	mult=1 k=4 l=1 Odd *
 1 -1 12 -1	mult=1 k=5 l=1 Odd *
 -1 0 12 -1	mult=1 k=7 l=1 Odd *

 1 0 11 2	mult=1 k=8 l=2  *
 0 2 11 0	mult=1 k=9 l=2  *
 1 0 12 0	mult=1 k=10 l=2  *
 -1 1 11 2	mult=1 k=11 l=2  *
 -1 1 12 0	mult=1 k=12 l=2  *


 **** full layer 2, obtained by even iteration
 **** wrong, we do not expect multiplicity 3, we expect 2
 ****  we have 18 states in place of 15, why ?

1 0 11 2	mult=1 k=8 l=2  *
 1 0 12 0	mult=3 k=10 l=2  *
 1 0 13 -2	mult=1 k=16 l=2  

 0 -1 12 2	mult=1 k=14 l=2  *
 0 -1 13 0	mult=3 k=15 l=2  *
 0 -1 14 -2	mult=1 k=18 l=2  

 -1 1 11 2	mult=1 k=11 l=2  *
 -1 1 12 0	mult=3 k=12 l=2  *
 -1 1 13 -2	mult=1 k=17 l=2  



 2 -2 13 0	mult=1 k=13 l=2  
 0 2 11 0	mult=1 k=9 l=2  *
 -2 0 13 0	mult=1 k=19 l=2  


 ************************************************************************

 SU(3) rep 0:2 == dim 6

...............Demazure, pass 3, r=0 dim = 6 sdim = 6
 0 2	mult=1 k=1 l=0  *
 1 0	mult=1 k=2 l=0  *
 2 -2	mult=1 k=3 l=0  
 -1 1	mult=1 k=4 l=0  *
 0 -1	mult=1 k=5 l=0  
 -2 0	mult=1 k=6 l=0  

 **************************
 SU(3) rep 0:1 ................Demazure, pass 3, r=0 dim = 3 sdim = 3
 0 1	mult=1 k=1 l=0  *
 1 -1	mult=1 k=2 l=0  *
 -1 0	mult=1 k=3 l=0  *

 **************************
  SU(3) rep 1:0 ..............Demazure, pass 2, r=1 dim = 3 sdim = 3
 1 0	mult=1 k=1 l=0  *
 -1 1	mult=1 k=2 l=0  *
 0 -1	mult=1 k=3 l=0  *

*/
