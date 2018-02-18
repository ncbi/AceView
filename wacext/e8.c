#include "ac.h"
#include "bitset.h"

/* 
   This program does some Lie algebra theory
   Given the Cartan matrix and the kac Dynkin weigths
   I construct, by the Demazure algorithm, all the roots.

   The next step is to construct the 3 * (rank+1)
   pools of all roots having scalar product 1,0,-1
   with respect to the simple roots and the higest root.

   The next step is to test the distance between the pools.
   and to study, in a way i do not remember if they are ok
   for the problem of Nicolas.
*/

typedef struct rootStruct { 
  /*********** Description of the algebra ********************/
  char title [64] ;  /* name of the algebra, i.e.e A2, E8... */

  /* The system of simple roots */
  int l ;               /* rank of the algebra */
  Array cartan ;        /* The Cartan matrix of dim l times l */
  Array simpleRoots ;   /* The simple roots, e1, e2, ... el i.e. delta(i,j) */
  Array kd ;            /* The highest root, i.e. the Kac-Dynkin labels */
  Array rotate ;        /* Matrix to change basis and use kd, e2, ... el */
  /* All roots */
  int dim ;      /* dimension = total number of roots (counting te origin only once */
  Array roots ;  /* Each root as a linear sum of simple roots (contravariant coordinates) */
  int *rho ;     /* the Weyl half sum of the positive roots */

  /************ Pooling strategy ******************************/
  int nPools ;   /* number of pools contructed */
  Array pools ;
  int nTests ;         /* number of test considered i.e. all singlet + all doublets */
  Array tests ;        /* which roots belong to each test */
  Array signatures ;   /* which pools are positives for each test */
  Array distances ;    /* table of distances between the signatures */

  AC_HANDLE h  ; } *ROOT ;
/***************************************************************/

static void rootFinalize (void *v)
{
  ROOT g = (ROOT) v ; 
  messfree (g->h) ;
}

static BOOL readCartan (ROOT g, char *cartan, char *kd)
{
  int i, j, k, x ;

  g->cartan = arrayHandleCreate (g->l * g->l, int, g->h) ;
  g->kd = arrayHandleCreate (g->l, int, g->h) ;
  
  freeforcecard (cartan) ; /* ignores the \n */
  for (i = k = 0 ; i < g->l ; i++)
    for (j = 0 ; j < g->l ; j++)
      if (freeint (&x))
	array (g->cartan, k++, int) = x ;
      else
	messcrash ("cannot parse the cartan matrix over = = %d, j = %d", i, j) ;
  
  freeforcecard (kd) ; /* ignores the \n */
  for (i = k = 0 ; i < g->l ; i++)
    if (freeint (&x))
      array (g->kd, k++, int) = x ;
    else
      messcrash ("cannot parse the cartan matrix over = = %d, j = %d", i, j) ;
  
  return TRUE ;
}

static BOOL showCartan (ROOT g)
{
  int i, j, k, x ;
  
  printf ("Cartan matrix of %s, rank = %d", g->title, g->l) ;
  for (i = k = 0 ; g && i < g->l ; i++)
    {
      printf ("\n") ;
      for (j = 0 ; j < g->l ; j++)
	{
	  x = array (g->cartan, k++, int) ;
	  printf ("%4d", x) ;
	}
    }
  printf ("\n\n") ;
  return TRUE ;
}

static ROOT  rootCreate (char *algebra, AC_HANDLE h)
{
  ROOT g = (ROOT) halloc (sizeof (struct rootStruct), h) ;
  char *kd, *cartan = 0 ;

  blockSetFinalise (g, rootFinalize) ;

  g->h = handleCreate () ;

  strncpy (g->title, algebra, 63) ;
  g->roots = arrayHandleCreate (1024, int, g->h) ;

  if (!strcmp (algebra, "A2"))
    {
      g->l = 2 ;
      kd = " 1 1 " ;
      cartan =
	" 2 -1 "
	"-1  2 " ;
    }
  if (!strcmp (algebra, "A3"))
    {
      g->l = 3 ;
      kd = " 1 1 1 " ;
      cartan =
	" 2 -1  0 "
	"-1  2 -1 " 
	"0  -1  2 " ;
    }
  if (!strcmp (algebra, "E8"))
    {
      g->l = 8 ;
      kd = "2 3 4 5 6 3 4 2 " ;
      cartan =
	" 2 -1  0  0  0  0  0  0 "
	"-1  2 -1  0  0  0  0  0 " 
	" 0 -1  2 -1  0  0  0  0 "
	" 0  0 -1  2 -1  0  0  0 "
	" 0  0  0 -1  2 -1 -1  0 "
	" 0  0  0  0 -1  2  0  0 "
	" 0  0  0  0 -1  0  2 -1 "
	" 0  0  0  0  0  0 -1  2 " ;
    }
  if (cartan)
    readCartan (g, cartan, kd) ;
  else
    messcrash ("Undefined algebra") ;
  return g ;
}

static BOOL constructSimpleRoots (ROOT g)
{
  int i ;
  g->simpleRoots = arrayHandleCreate (g->l * g->l, int, g->h) ;

  for (i = 0 ; i < g->l ; i++)
    array (g->simpleRoots, i * g->l + i, int) = 1 ;
  return TRUE ;
}

static int sp (ROOT g, int *ip, int *jp) ;

static void showAllRoots (ROOT g, BOOL doShow)
{
  int jj, i2, *hr, *rr, s ;
  
  hr = g->rho ;
  printf ("This algebra %s of rank %d has %d roots", 
	  g->title, g->l, g->dim + g->l - 1) ;
  if (doShow)
    for (jj = 0 ; jj < g->dim ; jj++)
      {
	printf ("\n%4d: ", jj) ;
	for (i2 = 0 ; i2 < g->l ; i2++)
	  printf ("%3d", array (g->roots, jj * g->l + i2, int)) ;
	rr = arrp (g->roots, g->l * jj, int) ;
	s = sp (g, hr, rr) ;
	printf ("   s=%3d", s) ;
      }
 printf ("\n\n") ;
}
 
static int sp (ROOT g, int *ip, int *jp)
{
  int s = 0 ;
  int i, j ;

  for (i = 0 ; i < g->l ; i++)
    for (j = 0 ; j < g->l ; j++)
      s += (*(ip+i)) * 
	arr (g->cartan, i*g->l + j, int) *
	(*(jp+j));
  return s ;
}

static BOOL constructAllRoots (ROOT g)
{
  BOOL isSame, isNew ;
  int ii, j, jj, i2, s, s1, n, n1, rr[100] ;
  g->roots = arrayHandleCreate (g->l * 256, int, g->h) ;

  /* define the highest root */
  for (j = 0 ; j < g->l ; j++)
    array (g->roots, j, int) = array (g->kd, j, int) ;
  /* run demazure on set of roots untill stable */
  n = n1 = 1 ;
  while (n1)
    {
      if (0)
	{
	  printf ("test: n1=%d", n1) ;
	  showAllRoots (g, TRUE) ;
	}
      n  = arrayMax (g->roots)/g->l ; n1 = 0 ;
      for (ii = 0 ; ii < g->l ; ii++) /* for each simple root sr */
	for (jj = 0 ; jj < n ; jj++) /* for each existing root */
	  {
	    /* compute the scalar product of r and sr */
	    s = sp (g, 
		    arrp (g->simpleRoots, ii * g->l, int),
		    arrp (g->roots, jj * g->l, int)) ;
	  
	    for (s1 = 1 ; s1 <= s ; s1++) /* Demazure completion */
	      {
		/* construct the Weyl symmetric root */
		for (i2 = 0 ; i2 < g->l ; i2++)
		  rr[i2] =
		    array (g->roots, jj * g->l + i2, int) 
		    - s1 * array (g->simpleRoots, ii * g->l + i2, int);
		/* if this root does not yet exist it should be created */
		isNew = TRUE ;
		for (j = 0 ; isNew && j < n + n1 ; j++)
		  {
		    isSame = TRUE ;
		    for (i2 = 0 ; isSame && i2 < g->l ; i2++)
		      if (rr[i2] != array (g->roots, j * g->l + i2, int))
			isSame = FALSE ;
		    if (isSame)
		      isNew = FALSE ;
		  }
		if (isNew) /* ok, create it */
		  {
		    for (i2 = 0 ; i2 < g->l ; i2++)
		      array (g->roots, (n+n1) * g->l + i2, int) = rr [i2] ;
		    n1++ ;
		  }		
	      }
	  }
    }
  g->dim = n ;
  return TRUE ;
}

static void constructRho (ROOT g)
{
  int i, ii, isNeg, *rr ;
  g->rho = (int *) halloc (g->l * sizeof (int), g->h) ;
  
  for (ii = 0 ; ii < g->dim ; ii++)
    {
      rr = arrp (g->roots, g->l * ii, int) ;
      for (i = isNeg = 0 ; !isNeg && i < g->l ; i++)
	if (*(rr+i) < 0)
	  isNeg = TRUE ;
      if (isNeg)
	continue ;
      for (i = 0 ; i < g->l ; i++)
	*(g->rho + i) += *(rr + i) ;
    }
  for (i = 0 ; i < g->l ; i++)
    *(g->rho + i) /= 1 ;
}

static Array matMultiply (Array aa, Array bb, int l, AC_HANDLE h)
{
  int i, j, k, x ;
  Array cc = arrayHandleCreate (l * l, int, h) ;
  for (i = 0 ; i < l ; i++)
    for (j = 0 ; j < l ; j++)
      {
	x = 0 ;
	for (k = 0 ; k < l ; k++)
	  x += arr (aa, i * l + k, int) *  arr (bb, k * l + j, int) ;
	array (cc, i *l + j, int) = x ;
      }
  return cc ;
}
/* construct the rotation matrix to use
 *  e2 e3 e4.. kd */
static BOOL showRotate (ROOT g)
{
  int i, j, k, x ;
  
  printf ("Rotation matrix of %s, rank = %d", g->title, g->l) ;
  for (i = k = 0 ; g && i < g->l ; i++)
    {
      printf ("\n") ;
      for (j = 0 ; j < g->l ; j++)
	{
	  x = array (g->rotate, k++, int) ;
	  printf ("%4d", x) ;
	}
    }
  printf ("\n\n") ;
  /* verif that rotate ^ l = identity ??? */
  {
    Array a1 = 0, an = 0 ;
    
    an = arrayCopy (g->rotate) ; /* power 1 */
    for (i = 1 ; i <= g->l ; i++)
      {
	a1 = an ;
	an = matMultiply (a1, g->rotate, g->l, g->h) ;
	arrayDestroy (a1) ;
      }

    printf ("Rotation matrix power l=%d", g->l) ;
    for (i = k = 0 ; g && i < g->l ; i++)
      {
	printf ("\n") ;
	for (j = 0 ; j < g->l ; j++)
	  {
	    x = array (an, k++, int) ;
	    printf ("%4d", x) ;
	  }
      }
    printf ("\n\n") ;
    arrayDestroy (an) ;
  }
  return TRUE ;
}

static void constructRotate (ROOT g)
{
  int ii ;

  g->rotate = arrayHandleCreate (g->l * g->l, int, g->h) ;
  
  array (g->rotate, g->l * g->l - 1 , int) = 0 ; 
  for (ii = 0 ; ii < g->l; ii++)
    {
      if (ii) array (g->rotate, g->l * ii + ii - 1 , int) = 1 ; 
      array (g->rotate, g->l * ii + g->l - 1 , int) = - array (g->kd, ii, int) ; 
    }
}


static void constructPools (ROOT g)
{
  /* for each simple root and the highest root
     i make a pool with every value of the scalar product
  */
  BitSet bb ; 
  int np, ii, jj, s, *rr, *hr ;

  g->pools = arrayHandleCreate (3 * (g->l + 1), BitSet, g->h) ;
  hr = g->rho ; /* the sum of the positive roots */
  for (ii = 0 ; ii <= g->l ; ii++)
    {
      if (ii == g->l)
	{
	  rr = arrp (g->roots, 0, int) ; /* the highest root */
	  rr = g->rho ;
	}
      else
	rr = arrp (g->simpleRoots, ii * g->l, int) ;
      /* loop on all positive roots */
      for (jj = 0 ; jj <= g->dim ; jj++)
	{
	  s = sp (g, hr,  arrp (g->roots, jj * g->l, int)) ;
	  if (s<0) /* keep only the positive roots */
	    continue ;
	  if (0) /* sue scalar products */
	    {
	      s = sp (g, rr,  arrp (g->roots, jj * g->l, int)) ;
	      if (ii == g->l) s /= 2 ;
	    }
	  else /* use components */
	    {
	      if (ii < g->l)
		s = arr (g->roots, jj *g->l + ii, int) ;
	      else
		{ /* rotate the root and take the last component */
		  int k ;
		  for (s = k = 0 ; k < g->l ; k++)
		    s += arr (g->rotate, (g->l-1) * g->l + k, int) * arr (g->roots, jj *g->l + k, int) ;
		}
	    }
	  s= ((s %2) + 2 ) %2 ;
         if (s <= -1 || s > 1)
	    messcrash ("bad scalar product") ;
	  np = 2 * ii + (0 + s) ; 
	  bb = array (g->pools, np, BitSet) ;
	  if (!bb)
	    bb = array (g->pools, np, BitSet) = bitSetCreate (g->dim, 0) ;
	  bitSet (bb, jj) ;
	}
    }
  g->nPools = arrayMax (g->pools) ;
}

static void showPools (ROOT g)
{
  int np, i ;
  BitSet bb ;

  printf ("\nConstructed %d pools", g->nPools) ;
  for (np = 0 ; np < g->nPools ; np++)
    {
      bb = array (g->pools, np, BitSet) ;
      if (!bb) continue ;
      printf ("\n%4d:", np) ;
      for (i=0 ; i < g->dim; i++)
	if (bit (bb, i))
	  printf ("%4d",i) ;
    }
  printf ("\n\n") ;

}

/* return a bit set of all pools positive for roots in pRoots */
static BitSet positivePools (ROOT g, BitSet pRoots)
{
  BitSet bp, pp = bitSetCreate (g->nPools, g->h) ;
  int i, j ;
  
  for (i = 0 ; i < g->nPools ; i++)
    {
      bp = array (g->pools, i, BitSet) ;
      for (j = 0 ; j < g->dim ; j++)
	if (bit (pRoots, j) && bit (bp, j))
	{
	  bitSet (pp, i) ;
	  break ;
	}
    }
  return pp ;
}

static int pDistance (ROOT g, BitSet b1, BitSet b2)
{
  int j, n = 0 ;

  for (j = 0 ; j < g->dim ; j++)
    if (
	(bit (b1, j) && !bit (b2, j)) ||
	(!bit (b1, j) && bit (b2, j))
	)
      n++ ;
  return n ;
}

static int constructTests (ROOT g, int level)
{
  int i, j, jj = 0, *hr, *rr ;
  BitSet bb ;

  g->tests = arrayCreate (100, BitSet) ;
  g->signatures = arrayCreate (100, BitSet) ;
  hr = g->rho ;
   /* create a test for each single root */
  for (i = 0 ; level >= 1 && i < g->dim ; i++)
    {
      rr = arrp (g->roots, g->l * i, int) ; /* the hith root */
      if (sp (g, hr, rr) < 0)
	continue ;
      bb = array (g->tests, jj, BitSet) =
	bitSetCreate (g->dim, g->h) ;
      bitSet (bb, i) ;
      array (g->signatures, jj, BitSet) = positivePools (g, bb) ;
      jj++ ;
   }
  /* create a test for each pair of roots */
  for (i = 0 ; level >= 2 && i < g->dim ; i++)
    {
      rr = arrp (g->roots, g->l * i, int) ; /* the ith root */
      if (sp (g,hr, rr) < 0)
	continue ;
      for (j = i + 1 ; j < g->dim ; j++)
	{
	  rr = arrp (g->roots, g->l * j, int) ; /* the ith root */
	  if (sp (g,hr, rr) < 0)
	    continue ;
	  bb = array (g->tests, jj, BitSet) =
	    bitSetCreate (g->dim, g->h) ;
	  bitSet (bb, i) ;
	  bitSet (bb, j) ;
	  array (g->signatures, jj, BitSet) = positivePools (g, bb) ;
	  jj++ ;
	}
    }
	
  g->nTests = jj ;
  return jj ;
}

static void showTest (ROOT g, int ii)
{
  BitSet bb = array (g->tests, ii, BitSet) ;
  int i ;

  printf ("Test %3d:", ii) ;
  for (i = 0 ; i < g->dim ; i++)
    if (bit (bb, i))
      printf (" %d", i) ;
  printf ("\n") ;

  bb = array (g->signatures, ii, BitSet) ;
  printf ("sig %3d:", ii) ;
  for (i = 0 ; i < g->dim ; i++)
    if (bit (bb, i))
      printf (" %d", i) ;
  printf ("\n") ;
}

static void showTests (ROOT g)
{
  int ii ;
  
  for (ii = 0 ; ii < g->nTests ; ii++)
    showTest (g, ii) ;

}

static int constructDistances (ROOT g)
{
  int dd, imin = 0, jmin = 0, ii, jj, dMin = 999999 ;
 
  g->distances = arrayHandleCreate (g->nTests * g->nTests, int, g->h) ;
  for (ii = 0 ; ii < g->nTests ; ii++)
    for (jj = 0 ; jj < g->nTests ; jj++)
      {
	dd = pDistance (g, array (g->signatures, ii, BitSet),
			array (g->signatures, jj, BitSet)) ;
	array (g->distances, g->nTests * ii + jj, int) = dd ;
	if (ii != jj && dd < dMin) 
	  {
	    imin = ii ; jmin = jj ; dMin = dd ;
	  }
      }
  printf ("distances constructed mini = %d : %d %d\n", dMin, imin, jmin) ;
  showTest (g, imin) ;
  showTest (g, jmin) ; 
 return dMin ;
}

static void showDistances (ROOT g)
{
  int i, ii, jj ;
  BitSet bb ;

  for (ii = 0 ; ii < g->nTests ; ii++)
    {
      bb = array (g->tests, ii, BitSet) ;
      for (i = 0 ; i < g->dim ; i++)
	if (bit (bb, i))
	  printf ("%d ", i) ;
      printf (" :: ") ;
      for (jj = 0 ; jj < g->nTests ; jj++)
	printf ("%3d", array (g->distances, g->nTests * ii + jj, int)) ;
      printf ("\n") ;
    }
}

int main  (int argc, char *argv[])
{
  ROOT g = 0 ;
  AC_HANDLE h = handleCreate () ;

  g = rootCreate (argc > 1 ? argv[1] : "A2", h) ;
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

  messfree (h) ;
  printf ("A bientot\n") ;

  return 0 ;
}

