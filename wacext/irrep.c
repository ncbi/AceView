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
  n = up->k - vp->k ;
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
  return kk ;
} /* locateWeight */

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

static void getHighestWeight (SA *sa, int type)
{
  WW *hw ;
  Array wws ;
  int k = 0 ;
  int rank = sa->rank ;
  int r, r1 = sa->hasOdd - 1 ;

  /* reinitialize a single hw.w */
  sa->pass = 0 ;
  sa->dict = dictHandleCreate (32, sa->h) ;


  wws = sa->wws = arrayHandleCreate (64, WW, sa->h) ;
  hw = arrayp (wws, 1, WW) ;

  switch (type)
    {
    case 0: /* trivial h.w. */
      break ;
    case 1: /* use as h.w. the lowering simple root */
      for (r = 0 ; r < rank ; r++)
	hw->x[r] = - array(sa->Cartan, rank * r1 + r, int) ;
      break ;
    case 2: /* construct from the parameters */
      if (sa->wws)
	k = sscanf (sa->DynkinWeights, "%d:%d:%d:%d:%d:%d:%d:%d", &hw->x[0] , &hw ->x[1] , &hw ->x[2] , &hw ->x[3] , &hw ->x[4] , &hw ->x[5] , &hw ->x[6] , &hw ->x[7] ) ;
      else
	messcrash ("Missing argument --s, expectirn -wws 1:0:2;...    giving the Dynkin lables of the highest weight ") ;
      if (k != sa->rank)
	messcrash ("HW: r5ank = %d but you provided %d Dynkin weights\n", sa->rank, k) ;
      break ;
    }
  hw->mult = 1 ;
  hw->hw = TRUE ;
  hw->k = locateWeight (sa, hw, TRUE) ;
} /* getHighestWeight */

/******************************************************************************************************/
/* construct the KacCrystal */
static void getKacCrystal (SA *sa, BOOL show)
{
  Array wws ;
  Array oddRoots = sa->oddRoots ;
  int r, rank = sa->rank ;
  int k, kMax = arrayMax (oddRoots) - 1 ;
  int ii, iiMax = 1 ;
  WW *w1, *w ;
  int atypic[rank] ;

  getHighestWeight (sa, 2) ;
  w = arrp (sa->wws, 1, WW) ;

  memset (atypic, 0, sizeof (atypic)) ;
  if (w->x[sa->hasOdd-1] == 0)
    atypic[0] = 1 ;
  
  /* reinitialize on the trivial h.w */
  getHighestWeight (sa, 0) ;
  wws = sa->wws ;
  if (show)
    printf("\n####### Kac Crystal \n") ;
  for (k = 0 ; k < kMax ; k++)
    iiMax *= 2 ;   /* 2^oddroots = size of the KacCrystal */
  for (ii = 1 ; ii < iiMax ; ii++) /* all points in the Kac Crystal, h.w. already included */
    {
      int layer = 0 ;
      int x = ii ;
      BOOL ok = TRUE ;
      vTXT txt = vtxtHandleCreate (0) ;
      WW w0 ;

      memset (&w0, 0, sizeof (w0)) ;

      w1 = arrayp (wws, 1, WW) ; /* the heighest weight */
      w1->hw = 1 ;
      w0 = *w1 ;
      for (k = 0 ; ok && k < kMax ; k++)
	{
	  int yes  = x % 2 ; /* odd root k is used in ii */
	  x /= 2 ;
	  if (yes == 1)
	    {
	      if (atypic[k])
		ok = FALSE ;
	      else
		{
		  WW *wodd = arrayp (oddRoots, k + 1, WW) ;
		  
		  vtxtPrintf (txt, " %d", k) ;
		  layer++ ;
		  for (r = 0 ; r < rank ; r++)
		    w0.x[r] += wodd->x[r] ;
		}
	    }
	}

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
	      printf ("Kac crystal: ii=%d  yes:%s :: ", ii, vtxtPtr (txt)) ; 
	      for (r = 0 ; r < rank ; r++)
		printf (" %d", w->x[r]) ;
	      printf ("\tmult=%d k=%d l=%d %s %s\n", w->mult, w->k, w->layer, w->odd ? "Odd" : "", w->hw ? "*" : "" ) ;
	    }
	}
      ac_free (txt) ;
    }
  sa->Kac = sa->wws ;
  sa->wws = 0 ;
} /* getKacCrystal */

/******************************************************************************************************/


static void getTensorProduct (SA *sa, BOOL show)
{
  int rank = sa->rank ;
  int r, ii, kk ;
  Array wws ;
  Array old = sa->wws ;
  Array Kac = sa->Kac ;
  int kMax = arrayMax (Kac) ;
  int iMax = arrayMax (old) ;
  
  sa->wws = 0 ;
  getHighestWeight (sa, 0) ; /* reinitialize wws and the dictionary */
  wws = sa->wws ;
  if (1)
    { /* avoid overcounting */
      WW *ww = arrp (wws, 1, WW) ;
      ww->mult = 0 ;
    }
  
  for (ii = 1 ; ii < iMax ; ii++)
    {
      WW *w1 = arrp (old, ii, WW) ;
      
      for (kk = 1 ; kk < kMax ; kk++)
	{
	  WW *wo = arrp (Kac, kk, WW) ;
	  WW w, *ww ;
	  int r, k2 ;

	  /* set the corrds of the next point of the tensor product */
	  memset (&w, 0, sizeof (w)) ;
	  w = *w1 ;
	  for (r = 0 ; r < rank ; r++)
	    w.x[r] += wo->x[r] ;
	  
	  /* locate it to construct the multiplicity */
	  k2 = locateWeight (sa, &w, TRUE)  ;
	  ww = arrayp (wws, k2, WW) ; /* new node */
	  ww->layer = wo->layer ;
	  ww->odd = wo->odd ;
	  ww->mult += w1->mult * wo->mult ;
	  for (r = 0 ; r < rank ; r++)
	    ww->x[r] = w.x[r] ;
	  if (1)
	    {
	      printf ("\n---------------------- ii=%d kk=%d k2 = %d m=%d\n", ii, kk, k2, ww->mult) ;
	      for (r = 0 ; r < rank ; r++)
		printf (" %d", w1->x[r]) ;
	      for (r = 0 ; r < rank ; r++)
		printf (" %d", wo->x[r]) ;
	      for (r = 0 ; r < rank ; r++)
		printf (" %d", ww->x[r]) ;
	    }
	}
    }

  if (show)
    {
      int dimE = 0, dimO = 0 ;
      printf ("\nTensorProduct:") ;
      iMax = arrayMax (wws) ;
      for (ii = 1 ; ii < iMax ; ii++)
	{
	  WW *ww = arrp (wws, ii, WW) ;
	  for (r = 0 ; r < rank ; r++)
	    printf (" %d", ww->x[r]) ;
	  printf ("\tmult=%d k=%d l=%d %s %s\n", ww->mult, ww->k, ww->layer, ww->odd ? "Odd" : "", ww->hw ? "*" : "" ) ;
	  if (ww->odd)
	    dimO += ww->mult ;
	  else
	    dimE += ww->mult ;
	}
      printf ("## TensorProduct dimE=%d dimO=%d dim=%d sdim=%d\n", dimE, dimO, dimE+dimO, dimE - dimO) ;    }
  
  return ;
} /* getTensorProduct */

/******************************************************************************************************/

static BOOL demazure (SA *sa, int r1, int *dimp, int *sdimp, BOOL show)
{
  BOOL new = FALSE ;
  Array wws = sa->wws ;
  int i, dim, sdim, rank = sa->rank ;
  BOOL odd = sa->odd[r1] ;
  
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
	  k2 = locateWeight (sa, w, FALSE) ;
	  if (k2) 
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
	      k2 = locateWeight (sa, w, TRUE)  ;
	      
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

  if (1)   arraySort (wws, wwLayerOrder) ;
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
  if (1)   arraySort (wws, wwCreationOrder) ;
  
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
  if (dimp) *dimp = dim ; 
  if (sdimp) *sdimp = sdim ;
  return ;
} /* demazureEven */

/******************************************************************************************************/


static void getOddRoots (SA *sa, BOOL show)
{
  int dim = 0 ;
  
  /* use as h.w. the lowering simple root */
  getHighestWeight (sa, 1) ;
  /* construct the first layer using Demazure */
  demazureEven (sa, &dim, 0, show) ;
  sa->oddRoots = sa->wws ;
  sa->wws = 0 ;

  printf ("# Constructed %d odd roots\n", dim) ;

  return ;
} /* getOddRoots */

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

   if (sa.hasOdd)
     {                        /* contruct the kasCrystal */
       getOddRoots (&sa, show) ;
       getKacCrystal (&sa, show) ; 
     }
   
   /* apply Demazure of the even group */
   getHighestWeight (&sa, 2) ;
   printf ("## Weights of the representations\n") ;
   
   demazureEven (&sa, &dim, &sdim, show) ;
   printf ("*************************** %s m=%d n=%d %s dim=%d sdim = %d\n "
	   , sa.type, sa.m, sa.n, sa.DynkinWeights, dim, sdim) ;

   if (sa.hasOdd)
     { 
       getTensorProduct (&sa, show) ;
     }
   
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
