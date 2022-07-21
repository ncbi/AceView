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
  int hasAtypic ;
  int nEven, nOdd ;
  BOOL odd[RMAX] ;
  BOOL extended[RMAX] ;
  Array Kac ; /* Kac crystal */
  WW evenHw ;
  Array evenRoots ; /* odd roots */
  Array oddRoots ; /* odd roots */
  Array wws ; /* weights */
  WW hw ;
  DICT *dict ;
  int pass ;
  int D1, D2 ; /* number of even and odd generators of the adjoiunt rep */
  Array dd ;   /* dims of all the submodules */
  Array atypic ;
  Array wwsShifted ;
  MX chi ;
} SA ; /* SuperAlgebra struct */


static int demazure (SA *sa, int *dimp, BOOL show) ;

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
	    sa->extended[r] = TRUE ;
	    array (Cartan, r*o + o , int) = 0 ;
	    if (o > 0) array (Cartan, r*o + o - 1, int) = 1 ;
	    if (o < r-1) array (Cartan, r*o + o + 1, int) = 1 ;
	  }
	
	if (n > 1)
	  {
	    sa->evenHw.x[m] += 1 ;
	    sa->evenHw.x[m + n -2] += 1 ;
	  }
      }

      if (m > 0)
	{
	  sa->evenHw.x[0] += 1 ;
	  sa->evenHw.x[m - 2] += 1 ;
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

      sa->evenHw.x[0] = 1 ;
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

      sa->evenHw.x[0] = 1 ;
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

	sa->evenHw.x[m-1] = 1 ;
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

	sa->evenHw.x[0] = 1 ;
      }
      break ;
      
    case 'G':
      if (m != 2 && m != 3)
	messcrash ("Type Lie G(2) or Kac G(3): m should be 2 or 3 m=%d n=%d", m,n) ;
      if (m == 3)
	n = 1 ;
      sa->rank = r = m + n ;
      rr = r*r ;
      Cartan = arrayCreate (rr, int) ;
      {
	array (Cartan, r*0 + 0, int) = 2 ;
	array (Cartan, r*0 + 1, int) = -1 ;
	array (Cartan, r*1 + 1, int) = 2 ;
	array (Cartan, r*1 + 0, int) = -3 ;

	sa->evenHw.x[0] = 1 ;
      }
      if (m == 3)
	{
	  sa->odd[2] = TRUE ;
	  array (Cartan, r*2 + 2 , int) = 0 ;
	  array (Cartan, r*2 + 1, int) = -1 ;
	  array (Cartan, r*3 + 2 , int) = 1 ;
	  array (Cartan, r*3 + 1, int) = 2 ;

	  sa->evenHw.x[3] = 1 ;
	  sa->extended[3] = 1 ;
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

static void wwsShow (SA *sa, char *title, int type, Array wws)
{
  if (wws)  
    {
      int dimE = 0, dimO = 0 ;
      int ii, iMax = arrayMax (wws) ;
      int r, rank = sa->rank ;

      if (0)arraySort (wws, wwLayerOrder) ;
      printf ("\n##################### %s: t=%d : \n", title, type) ;

      for (ii = 1 ; ii < iMax ; ii++)
	{
	  WW *ww = arrp (wws, ii, WW) ;
	  if (ww->mult)
	    {
	      for (r = 0 ; r < rank ; r++)
		printf (" %d", ww->x[r]) ;
	      printf ("\tmult=%d k=%d l=%d %s %s\n", ww->mult, ww->k, ww->layer, ww->odd ? "Odd" : "", ww->hw ? "*" : "" ) ;
	      if (ww->odd)
		dimO += ww->mult ;
	      else
		dimE += ww->mult ;
	    }
	}
      printf ("## %s dimE=%d dimO=%d dim=%d sdim=%d\n", title, dimE, dimO, dimE+dimO, dimE - dimO) ;
      if (0) arraySort (wws, wwCreationOrder) ;
    }
  return ;
} /* wwsShow */

/******************************************************************************************************/

static void getHighestWeight (SA *sa, int type, BOOL create, BOOL show)
{
  WW *hw ;
  Array wws ;
  int k = 0 ;
  int rank = sa->rank ;
  int r, r1 = sa->hasOdd - 1 ;

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
      for (r = 0 ; r < rank ; r++)
	hw->x[r] = - array(sa->Cartan, rank * r + r1, int) ;
      break ;
    case -2: /* construct from the parameters */
      if (sa->wws)
	k = sscanf (sa->DynkinWeights, "%d:%d:%d:%d:%d:%d:%d:%d", &hw->x[0] , &hw ->x[1] , &hw ->x[2] , &hw ->x[3] , &hw ->x[4] , &hw ->x[5] , &hw ->x[6] , &hw ->x[7] ) ;
      else
	messcrash ("Missing argument --s, expectirn -wws 1:0:2;...    giving the Dynkin lables of the highest weight ") ;

      if (k != sa->rank)
	messcrash ("HW: r5ank = %d but you provided %d Dynkin weights\n", sa->rank, k) ;
      break ;
    case -3:
      *hw = sa->evenHw ;
      break ;
    default:
      if (sa->oddRoots && type <= arrayMax (sa->oddRoots))
	{
	  WW *wo = arrp (sa->oddRoots, type, WW) ;
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

  if (1 || type == -2)
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
} /* getHighestWeight */

/******************************************************************************************************/
/* construct the KacCrystal */
static void getKacCrystal (SA *sa, BOOL show)
{
  Array wws ;
  Array oddRoots = sa->oddRoots ;
  int r, rank = sa->rank ;
  int k, kMax = sa->nOdd ;
  int ii, iiMax = 1 ;
  WW *w ;

  /* reinitialize on the trivial h.w */
  sa->dict = dictHandleCreate (32, sa->h) ;
  getHighestWeight (sa, 0, 0, show) ;
  wws = sa->wws ;

  for (k = 0 ; k < kMax ; k++)
    iiMax *= 2 ;   /* 2^oddroots = size of the KacCrystal */
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
	      WW *wodd = arrayp (oddRoots, k + 1, WW) ;
	      
	      vtxtPrintf (txt, " %d", k) ;
	      layer++ ;
	      for (r = 0 ; r < rank ; r++)
		w0.x[r] += wodd->x[r] ;
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
	  if (1)
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
  sa->Kac = sa->wws ;
  sa->wws = 0 ;
  if (show)
    wwsShow (sa, "Kac crystal dim =", 1+arrayMax(sa->Kac), sa->Kac) ;
  if (0)
    printf ("# Constructed %d Kac weigthst\n", arrayMax (sa->Kac) - 1) ;
} /* getKacCrystal */

/******************************************************************************************************/


static Array getOneTensorProduct (SA *sa, int atypic,   Array old, BOOL show)
{
  int rank = sa->rank ;
  int ii, kk ;
  Array wws ;
  Array Kac = sa->Kac ;
  int kMax = arrayMax (Kac) ;
  int iMax = arrayMax (old) ;
  WW deltaW ;

  sa->wws = 0 ;
  memset (&deltaW, 0, sizeof(deltaW)) ;
  /* reinitialize wws and but not the dictionary */
  getHighestWeight (sa, -2, 0, show) ;
  if (atypic > 0)
    deltaW = arr (sa->oddRoots, atypic, WW) ;

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
	    w.x[r] += wo->x[r] - deltaW.x[r] ;
	  
	  /* locate it to construct the multiplicity */
	  k2 = locateWeight (sa, &w, TRUE)  ;
	  ww = arrayp (wws, k2, WW) ; /* new node */
	  ww->k = k2 ;
	  ww->layer = wo->layer ;
	  ww->odd = wo->odd ;
	  ww->mult += w1->mult * wo->mult ;
	  for (r = 0 ; r < rank ; r++)
	    ww->x[r] = w.x[r] ;
	  if (0)
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
    wwsShow (sa, "Tensor Product", atypic, wws) ;
  return wws ;
} /* getOneTensorProduct */


/******************************************************************************************************/

static void intersectTensorProducts (SA *sa, Array wws, Array xxs)
{
  WW *ww, *xx ;
  int j ;
  
  wwsShow (sa, "INTERSECT wws", 0, wws) ;
  wwsShow (sa, "INTERSECT xxs", 1, xxs) ;
  for (j = 0 ; j < arrayMax (wws) ; j++)
    { /* we use a single dict */
      ww = arrp (wws, j, WW) ;
      xx = arrayp (xxs, j, WW) ;
      if (xx->mult < ww->mult)
	ww->mult = xx->mult ;
    }
  wwsShow (sa, "SUBTRACTED wws", 0, wws) ;
  return ;
} /* intersectTensorProducts */

/******************************************************************************************************/

static void getTensorProducts (SA *sa, int *dimp, int *sdimp,  BOOL show)
{
  Array old = sa->wws ;
  Array wws = getOneTensorProduct (sa, -2, old, show) ;
  if (sa->hasAtypic)
    {
      int ii ;
      
      for (ii = 0 ; ii < arrayMax (sa->atypic) ; ii++)
	if (array (sa->atypic, ii, int))
	  {
	    Array top, xxs ;
	    getHighestWeight (sa, -ii, 0, show) ;
	    demazure (sa, 0, show) ;
	    top = sa->wws ;
	    xxs = getOneTensorProduct (sa, ii, top, show) ;
	    intersectTensorProducts (sa, wws, xxs) ;
	    ac_free (xxs) ;
	    ac_free (top) ;
	    
	    if (show)
	      wwsShow (sa, "Subtracted Atypical Tensor Product", ii, wws) ;
	  }
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
  if (show)
    wwsShow (sa, "Final Atypical Tensor Product", 999, wws) ;

  return ;
} /* getTensorProducts */

/******************************************************************************************************/

static BOOL demazureEven (SA *sa, int r1, int *dimp, BOOL show)
{
  BOOL new = FALSE ;
  Array wws = sa->wws ;
  int i, dim, sdim, rank = sa->rank ;
  
  if (sa->pass++ == 0) new = TRUE ;
    
  for (i = 1 ; i < arrayMax (wws) ; i++) 
    {
      WW *w = arrayp (wws, i, WW) ;
      int n1 = w->mult ;
      int n2 = 0 ;
      int na = 0 ;
      int dn = 0 ;
      int jMax = w->x[r1] ; /* default even values */

	  
      if (jMax > 0)  /* create or check existence of the new weights along the sl(2) submodule */
	{
	  int j, r, k2 ;
	  WW *w2 ;
	  
	  for (j = 1 ; j <= jMax ; j++)
	    {
	      /* position to the new weight */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] -= array(sa->Cartan, rank * r + r1, int) * j ;
	      
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
		}
	      n2 = na + w2->mult ;

	      /* reposition w to its original location */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] += array(sa->Cartan, rank * r + r1, int) * j ;
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
	      
	      /* reposition w to its original location */
	      for (r = 0 ; r < rank ; r++)
		w->x[r] += array(sa->Cartan, rank * r + r1, int) * j ;
	    }
	}
      
    }

  for (dim = sdim = 0, i = 1 ; i < arrayMax (wws) ; i++)
    {
      WW *w = arrp (wws, i, WW) ;
      dim += w->mult ;
    }

  if (1)   arraySort (wws, wwLayerOrder) ;
  if (0 && new)
    {
      printf ("................Demazure, pass %d, r=%d dim = %d sdim = %d\n", sa->pass, r1, dim + sdim, dim - sdim) ;
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
  
  *dimp = dim + sdim ;
  
  return new ;
} /* demazureEven */

/******************************************************************************************************/

static int demazure (SA *sa, int *dimp, BOOL show)
{
  BOOL ok = TRUE ;
  int r, dim = 0 ;

  while (ok)
    {
      ok = FALSE ;
      for (r = 0 ; r < sa->rank ; r++)
	if (! sa->odd[r])
	  ok |= demazureEven (sa, r, &dim, show) ;
    }
  if (dimp) *dimp = dim ; 

  if (show)
    wwsShow (sa, "Demazure", 0, sa->wws) ;
  
  return dim ;
} /* demazure */

/******************************************************************************************************/


static void getOddRoots (SA *sa, BOOL show)
{
  int dim = 0 ;

   sa->dict = dictHandleCreate (32, sa->h) ;

   if (show)
    printf ("\n####### Odd roots \n") ;
   /* use as h.w. the lowering simple root */
   getHighestWeight (sa, -1, 1, show) ;
  /* construct the first layer using Demazure */
  sa->nOdd = demazure (sa, &dim, show) ;
  sa->oddRoots = sa->wws ;
  sa->wws = 0 ;

  if (show)
    wwsShow (sa, "Odd roots", 1, sa->oddRoots) ;
  printf ("# Constructed %d odd roots\n", dim) ;



  return ;
} /* getOddRoots */

/******************************************************************************************************/


static void getAdjoint (SA *sa, BOOL show)
{
  int dimE = 0 ;

  sa->dict = dictHandleCreate (32, sa->h) ;
  if (show)
    printf ("\n####### Odd roots \n") ;
  /* use as h.w. the lowering simple root */
  getHighestWeight (sa, -3, 1, show) ;
  /* construct the first layer using Demazure */
  sa->nEven = demazure (sa, &dimE, show) ;
  sa->evenRoots = sa->wws ;
  sa->wws = 0 ;

  if (show)
    wwsShow (sa, "Even Adjoint ", 1, sa->evenRoots) ;
  printf ("# Constructed %d adjoint even roots\n", dimE) ;

  return ;
} /* getOddRoots */

/******************************************************************************************************/

static void getAtypic (SA *sa, BOOL show)
{
  WW hw ;
  int rank = sa->rank ;
  int r, r1 = sa->hasOdd - 1 ;
  int z = 0 ;

  sa->atypic = arrayHandleCreate (arrayMax (sa->oddRoots), int, sa->h) ;  

  /* use as h.w. the declared h.w. */
  getHighestWeight (sa, -2, 1, TRUE) ;
   hw = sa->hw ;
  
  /* construct the sum of the even roots rho0 */
  /* construct the sum of the even roots rho1 */
  /* rho = rho0 - rho1 */
  /* the correct equation uses the half supersum of the roots, but i do not want to divide by 2 my integers */
  /* for each odd root beta[i] compute   <2 * hw + rho | beta[i]> , so i use 2*hw to compensate the 2*rho */
  /* if (zero, atypic[i] = TRUE */
  
  /* for the moment we do something crude for sl(m/1) ; */
  switch ((int)sa->type[0])
    {
    case 'A':
      for (r = z = 0 ; r + r1 < rank ; r++)
	{
	  if (hw.x[r1] == z )
	    {
	      array (sa->atypic, r+1, int) = 1 ;
	      sa->hasAtypic++ ;
	    }
	  z += hw.x[r1 - r - 1] + 1 ;
	}
      break ;
    default :
      break ;
    }
  if (1 || show)
    {
      int i ;
      
      if (sa->hasAtypic)
	{
	printf ("##############  Aypical") ;
	  for (i = 0 ; i < arrayMax (sa->atypic) ; i++)
	    if (array (sa->atypic, i, int))
	      printf (" %d", i) ;
	  printf ("\n") ;
	}
      else
	printf ("############## Typical\n") ;
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

   if (sa.hasOdd) /* do this first then destroy the dict */
     getOddRoots (&sa, show) ;

   /* apply Demazure of the even group */
   getAdjoint (&sa, show) ;
 
   
   if (sa.hasOdd)
     {                        /* contruct the kasCrystal */
       getAtypic (&sa, show) ;
       getKacCrystal (&sa, show) ; 
     }
   
   /* construct the top layer even module */
   sa.dict = dictHandleCreate (32, sa.h) ;
   getHighestWeight (&sa, -2, TRUE, TRUE) ;
   printf ("## Weights of the representations\n") ;

   printf ("*************************** %s m=%d n=%d %s dim=%d sdim = %d\n "
	   , sa.type, sa.m, sa.n, sa.DynkinWeights, dim, sdim) ;

   demazure (&sa, &dim, show) ;
   if (sa.hasOdd)
     { 
       getTensorProducts (&sa, &dim, &sdim, show) ;
     }
   
   printf  ("Final Representation dim=%d sdim=%d\n",  dim, sdim) ;
   messfree (h) ;
   printf ("A bientot\n") ;
   
   return 0 ;   
} /* main */

/************************************************************************/
/************************************************************************/
/************************************************************************/
