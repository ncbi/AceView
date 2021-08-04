#include "ac.h"
#include <complex.h>
#include "matrix.h"

/* Create june 2020
 * edited 27 jan 2021
 *
 * My personnal way of computing the renormalisability of gauge theories
 *
 * A diagram is a hooking of subdiagrams, like propagators and vertex
 * It is an ordered product of functions of g_ij, epsilon_ijkl, sigma_i, sigma_bar_i
 * and has dependence of momenta (p,q,r), at the numerator and at the denominator
 * In addition, it contains group matrices and structure constants.
 *
 * The idea is to represent each term as a structure containing such operators
 * We only want to compute the pole part
 * 
 * To reduce the problem to a numeric expression, we estimate the degree in p,q,r
 * and compute the correct number of formal derivatives with repect to p_mu, q_nu, r_rho.
 * Then the loop integral reduces to a polynome in g_ab, p, q, r
 * Separately we know the overall degree in sigma, sigma_bar, multiply by sigma and trace
 * We extract the overall epsilon by multiplying by an additional epsilon.
 * 
 * expand: allows to flatten the products of polynomes
 * contractIndices: sort and reduces the Einstein repeated indices and the epsilons
 * PaukiTrace:  transforms the ordered products of Pauli matrices into a polynome in g_mn and epsilon
 *
 * Z2: Construct the Feynman diagram of a propagator given the field types
 *        by hooking the relevant Feynam rules
 * Z3, Z4: Idem for 3-fields and 4-fields vertices
 * Ward: Construct the sum of all diagrams that are expected to give a Ward identity 
 *
 * Testing: We must check all subcomponents calculations
 *          and check that we recover the Yang-Mill Ward identites
 *          The idea is that the progam can rerun all checks each time it is modified
 *          So each component comes with its internal checks.
 *
 * Purpose: Prove that QAD is renormalisable.

 * A separate part of the program is meant to check the consistency of reps of SU(2/1)
 * The 8 complex matrices are enumerated, and we then verify the commutators and the anomalies
 * we have explicitly the leptons, the quarks, and several indecomposable reps
 */

#define minAbs 0.0000001
#define GMAX 32
typedef struct termStruct {
  int type ;
  complex float z ;/* complex scalar multiplier. If zero, the whole TT is NULL */
  char sigma[GMAX] ; /* sigma     matrices : non-commutative list of index "ab" means sigma_a sigma-bar_b */
  char sigB[GMAX] ;  /* sigma-bar matrices : non-commutative list of index "ab" means sigma-bar_a sigma_b */
  char g[GMAX] ;     /* Lorentz metric */
  char gg[GMAX] ;    /* group metric */
  char eps[GMAX] ; /* espislon anti symmetric set of n times 4 indices */ 
  char mm[4][GMAX] ;      /*(k p q r)_mu momenta :  "1 ab" means the product p_a p_b */
                          /*  "1 a"  "2 b" means the product p_a q_b */
  int  denom[4] ;     /* number of terms of the form 1/k^, 1/(k+p)^2, 1/(k+p+q)^2, 1/(k+p+q+r)^2 */
  char Id2 ; /* Pauli identity matrix, needed its value is 2 when we trace */
  char freeIndex[GMAX] ;
} TT ;


typedef struct polynomeStruct *POLYNOME ;
struct polynomeStruct {
  int id ;
  BOOL isFlat, isSum, isProduct ;
  POLYNOME p1, p2 ;
  TT tt ;
} ;

static POLYNOME newMultiSum (POLYNOME ppp[]) ;
static POLYNOME expand (POLYNOME pp) ;
static POLYNOME newPolynome (void) ;
static void showPol (POLYNOME pp) ;
static POLYNOME squareMomentaCleanUp (POLYNOME pp, char alpha) ;
static POLYNOME dimIntegralDo (POLYNOME pp, int pass) ;
static BOOL freeIndex (POLYNOME pp) ;

/***********************************************************************************************************************************************/
/***************************************************************** Utilities  ******************************************************************/
/***********************************************************************************************************************************************/

static char firstDummyIndex = 'a' ;
static char newDummyIndex (void)
{
  unsigned char cc = firstDummyIndex++ ;
  
  if (cc >= 254)
    messcrash ("too many dummy indices") ;
  return (char)cc  ;
}

/***********************************************************************************************************************************************/
/* look for all indices, 
 * edit = FALSE: count occurences of a,b,c...
 * edit = TRUE: replace by a continuous list liberating the high values
 */
static void reduceIndicesTtDo (char *cp, KEYSET ks, BOOL edit)
{  
  int i ;
  
  for (i = 0 ; *cp && i<GMAX ; cp++, i++)
    {
      if (! edit)
	keySet (ks, *(unsigned char *)cp) ++ ;
      else
	*cp = keySet (ks, *(unsigned char *)cp) ;
    }
  return ;
} /* reduceIndices */

/***********************************************************************************************************************************************/

static void reduceIndicesTt (POLYNOME pp, KEYSET ks, BOOL edit)
{
  reduceIndicesTtDo (pp->tt.g, ks, edit) ;
  reduceIndicesTtDo (pp->tt.gg, ks, edit) ;
  reduceIndicesTtDo (pp->tt.sigma, ks, edit) ;
  reduceIndicesTtDo (pp->tt.sigB, ks, edit) ;
  reduceIndicesTtDo (pp->tt.eps, ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[0], ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[1], ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[2], ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[3], ks, edit) ;
} /* reduceIndicesTt */

/***********************************************************************************************************************************************/

static void reduceIndicesDo (POLYNOME pp, KEYSET ks, BOOL edit)
{
  if (! pp) 
    return ;
  reduceIndicesDo (pp->p1, ks, edit) ;
  reduceIndicesDo (pp->p2, ks, edit) ;
  if (pp->tt.type)
    reduceIndicesTt (pp, ks, edit) ;
} /* reduceIndices */

/***********************************************************************************************************************************************/

static int reduceIndicesIsProduct (POLYNOME pp)
{
  int i, j, n ;
  KEYSET ks = keySetCreate () ;
  KEYSET ks2 = keySetCreate () ;
  KEYSET ks3 = keySetCreate () ;

  /* count the occurence of all indices */
  reduceIndicesDo (pp, ks, FALSE) ;
  /* contract the list by dropping unused values */
  for (i = j = 0 ; 'a' + i < keySetMax (ks) ; i++)
    {
      n = keySet (ks, 'a' + i) ;
      if (n == 1)    /* block and do not rename the free indices present once */
	{
	  keySet (ks2, 'a' + i)  = 'a' + i ;
	  keySet (ks3, 'a' + i)  = 1 ;       /* block it */ 
	  if (firstDummyIndex <= 'a' + i)
	    firstDummyIndex = 'a' + i + 1 ;
	}
      else if (n > 2)
	{
	  showPol (pp) ;
	  messcrash ("no index should be repeated %d > 2 times", n) ;
	}
    }
  for (i = 0 ; 'a' + i < keySetMax (ks) ; i++)
    {
      n = keySet (ks, 'a' + i) ;
      if (n == 2)
	{
	  for (j = 0 ; ; j++)
	    if (keySet (ks3, 'a' + j)  == 0)  /* first unblocked index  */
	      {
		keySet (ks2, 'a' + i)  = 'a' + j ;
		keySet (ks3, 'a' + j)  = 1 ;
		if (firstDummyIndex <= 'a' + j)
		  firstDummyIndex = 'a' + j + 1 ;
		break ;
	      }
	}
    }
  /* rename all indices */
  reduceIndicesDo (pp, ks2, TRUE) ;

  keySetDestroy (ks) ; 
  keySetDestroy (ks2) ; 
  keySetDestroy (ks3) ; 
  return 'a' + j ;
} /* reduceIndicesIsProduct */

/***********************************************************************************************************************************************/

static POLYNOME reduceIndices (POLYNOME pp)
{
  static int level = 0 ;

  level++ ;
  if (0 && level == 1)
    {
      char a = newDummyIndex () ;
      pp = squareMomentaCleanUp (pp, a) ;
    }
  if (pp && pp->isSum)
    {
      if (pp->p1) 
	pp->p1 = reduceIndices (pp->p1) ;
      if (pp->p2) 
	pp->p2 = reduceIndices (pp->p2) ;
    }
  if (pp && ! pp->isSum) reduceIndicesIsProduct (pp) ;
  level-- ;
  if (level == 0 && pp)
    {
      if (! pp->isSum)
	reduceIndicesIsProduct (pp) ;
    }
  return pp ;
} /* reduceIndices */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static int polOrder (const void *va, const void *vb)
{
  const POLYNOME p1 = *(POLYNOME *)va ;
  const POLYNOME p2 = *(POLYNOME *)vb ;
  float complex z1 = p1->tt.z ;
  float complex z2 = p2->tt.z ;
  int s ;
  p1->tt.z = 0 ;       
  p2->tt.z = 0 ;
  memset (p1->tt.freeIndex, 0, GMAX) ;
  memset (p2->tt.freeIndex, 0, GMAX) ;
  s = memcmp (&(p1->tt), &(p2->tt), sizeof(TT)) ;
  p1->tt.z = z1 ;       
  p2->tt.z = z2 ;       

  return s ;
}

/***********************************************************************************************************************************************/

static void sortPolGetSum (POLYNOME pp, Array aa)
{
  if (! pp)
    return ;
  if (pp->isSum)
    {
      sortPolGetSum (pp->p1, aa) ;
      sortPolGetSum (pp->p2, aa) ;
    }
  else if (pp)
    array (aa, arrayMax (aa), POLYNOME) = pp ;
  return ;
}

/***********************************************************************************************************************************************/

static POLYNOME sortReduceSum (Array aa)
{
  POLYNOME pp = 0 ;
  int ii, jj ;

  if (!aa)
    return 0 ;

  if (arrayMax (aa) > 1)
    arraySort (aa, polOrder) ;

  for (ii = 0 ; ii + 1 < arrayMax (aa) ; ii++)
    {	/* add identical terms */				
      int jj ;
      POLYNOME q1 = arr (aa, ii, POLYNOME) ; 
      if (! q1 || ! q1->tt.type)
	continue ;
      for (jj = ii + 1 ; jj < arrayMax (aa) ; jj++)
	{
	  POLYNOME q2 = arr (aa, jj, POLYNOME) ; 
	  if (q2)
	    {
	      int n = polOrder (&q1, &q2) ;
	      if (n == 0)
		{
		  q1->tt.z += q2->tt.z ;
		  arr (aa, jj, POLYNOME)  = 0 ;
		}
	      else if (n < 0)
		break ;
	    }
	}
    }
  for (ii = jj = 0 ; ii < arrayMax (aa) ; ii++)
    {  /* elimimate the killed terms */
      POLYNOME qq = arr (aa, ii, POLYNOME) ; 
      if (qq && (qq->isSum || qq->isProduct || qq->tt.z != 0))
	arr (aa, jj++, POLYNOME) = qq ;
    }
  arrayMax (aa) = jj ;
  
  if (arrayMax (aa) == 0)
    pp = 0 ;
  else if (arrayMax (aa) == 1)
    pp = arr (aa, 0, POLYNOME) ;
  else 
    {
      POLYNOME qq ;
      
      pp =  newPolynome () ;
      for (qq = pp, ii = 0 ; ii + 1 < arrayMax (aa) ; ii++)
	{
	  qq->isSum = TRUE ;
	  qq->p1 = arr (aa, ii, POLYNOME) ;
	  if (ii < arrayMax (aa) - 2)
	    { qq->p2 = newPolynome () ; qq = qq->p2 ; }
	}
      qq->p2 = arr (aa, ii, POLYNOME) ;
    }
  return pp ;
} /* sortReduceSum */
    
/***********************************************************************************************************************************************/

static POLYNOME sortPol (POLYNOME pp)
{
  if (! pp)
    return 0 ;
  else if (pp->isSum)
    {
      int ii = 0, jj = 0 ;
      Array aa = arrayCreate (32, POLYNOME) ;
      
      sortPolGetSum (pp, aa) ;
      for (ii = jj = 0 ; ii < arrayMax (aa) ; ii++)
	{
	  POLYNOME qq = arr (aa, ii, POLYNOME) ; 
	  qq = sortPol (qq) ;
	  if (qq)
	    arr (aa, jj++, POLYNOME) = qq ;
	}
      arrayMax (aa) = jj ;

      pp = sortReduceSum (aa) ;
      arrayDestroy (aa) ;
    }
  else if (pp->isProduct)
    {
      pp->p1 = sortPol (pp->p1) ;
      pp->p2 = sortPol (pp->p2) ;
      if (pp->p1 && pp->p2)
	return pp ;
      return 0 ;
    }
  else if (! pp->tt.type ||  cabs (pp->tt.z) < minAbs)
    return 0 ;

  return pp ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static BOOL ppIsNumber (POLYNOME pp)
{
  TT tt = pp->tt ;
  BOOL ok = tt.type ;
  if (tt.mm[0][0])  ok = FALSE ;
  if (tt.mm[1][0])  ok = FALSE ;
  if (tt.mm[2][0])  ok = FALSE ;
  if (tt.mm[3][0])  ok = FALSE ;
  if (tt.g[0])  ok = FALSE ;
  if (tt.gg[0])  ok = FALSE ;
  if (tt.sigma[0])  ok = FALSE ;
  if (tt.sigB[0])  ok = FALSE ;
  if (tt.eps[0])  ok = FALSE ;
  return ok ;
} /* ppIsNumber */

/***********************************************************************************************************************************************/
/*********************************************************** Display utilitie ******************************************************************/
/***********************************************************************************************************************************************/

static float  nicePrintFraction (const char *prefix, float x, const char *suffix);

static float complex nicePrint (const char *prefix, float complex z)
{
  float a = creal (z) ;
  float b = cimag (z) ;
  if (b == 0)
    a = nicePrintFraction (prefix, a, "") ;
  else if (a == 0)
    b = nicePrintFraction (prefix, b, "i") ;
  else
    {
      a = nicePrintFraction (prefix, a, "") ;
      b = nicePrintFraction (" + ", b, "i") ;
    }
  return a + I * b ;
} /* nicePrint */

/*******************************************************************************************/
static int niceInt (float x)
{
  int k = x >= 0 ? x + .01 : x - .01 ;
  return k ;
} /* niceInt */

/*******************************************************************************************/

static float nicePrintFraction (const char *prefix, float x, const char *suffix)
{
  int i = 10000 ;
  int s = x >= 0 ? 1 : -1 ;
  
  x *= s ;
  if (niceInt (10000*x) == 0)
    { x = 0 ; s = 1 ; }
  printf ("%s", prefix) ;
  if (niceInt (10000*x) == 10000 * niceInt (x))
    printf ("%.0f", s*x) ;
  else
    for (i = 2 ; i <= 720 ; i++)
      if (niceInt (10000*i*x) == 10000 * niceInt (i*x))
	{
	  printf ("%.0f/%d", s*i*x, i) ;
	  x = ((int) (i * x)) / ((float)i) ;
	  i = 10000 ;
	}
  if (i < 10000)
    printf ("%.3f", s*x) ;
  if (x != 0)
    printf ("%s", suffix) ;

  return s * x ;
} /* nicePrintFraction */

/*******************************************************************************************/

static void niceComplexShow (MX a)
{
  int i,j, iMax = a->shapes[0] ;
  const complex float *zc ;

  mxValues (a, 0, 0, &zc) ;
  printf ("## %s", a->name) ;
  for (i = 0 ; i < iMax ; i++)
    {
      printf ("\n") ;
      for (j = 0 ; j < iMax ; j++)
	nicePrint ("\t", zc[iMax*i + j]) ;
    }
  printf ("\n") ;
  return ;
} /* niceShow */

/*******************************************************************************************/

static void niceFloatShow (MX a)
{
  int i,j, iMax = a->shapes[0] ;
  const float *zf ;

  mxValues (a, 0, &zf, 0) ;
  printf ("## %s", a->name) ;
  for (i = 0 ; i < iMax ; i++)
    {
      printf ("\n") ;
      for (j = 0 ; j < iMax ; j++)
	nicePrintFraction ("\t", zf[iMax*i + j],"") ;
    }
  printf ("\n") ;
  return ;
} /* niceFloatShow */

/*******************************************************************************************/
static void niceIntShow (MX a)
{
  int i,j, iMax = a->shapes[0] ;
  const int *zi ;

  mxValues (a, &zi, 0, 0) ;
  printf ("## %s", a->name) ;
  for (i = 0 ; i < iMax ; i++)
    {
      printf ("\n") ;
      for (j = 0 ; j < iMax ; j++)
	printf ("\t%d", zi[iMax*i + j]) ;
    }
  printf ("\n") ;
  return ;
} /* niceIntShow */

/*******************************************************************************************/

static void niceShow (MX a)
{
  switch (a->type)
    {
    case MX_NULL: break ;
    case MX_BOOL:
    case MX_INT: niceIntShow (a) ; break ;
    case MX_FLOAT: niceFloatShow (a) ; break ;
    case MX_COMPLEX: niceComplexShow (a) ; break ;
    }
  return ;
} /* niceShow */

/*******************************************************************************************/

static void showTT (POLYNOME pp)
{
  TT tt ;
  int i ;
  
  if (!pp)
    return ;
  tt = pp->tt ;  
  if (!tt.type)
    return ;
  
  if (! ppIsNumber (pp))
    {
      if (cabs (tt.z + 1) < minAbs)
	printf (" -") ;
      else if (cabs (tt.z - 1) > minAbs)
	{ tt.z = nicePrint ("", tt.z) ; }
    }
  else
    { tt.z = nicePrint ("", tt.z) ; }
  
  if (*tt.g) printf (" g_%s ", tt.g) ;
  if (*tt.gg) printf (" gg_%s ", tt.gg) ;
  if (*tt.eps) printf (" epsilon_%s ", tt.eps) ;

  if (*tt.sigma) printf (" s_%s ", tt.sigma) ;
  if (*tt.sigB) printf (" sB_%s ", tt.sigB) ;

  if (tt.mm[0][0]) printf (" k_%s ", tt.mm[0]) ;
  if (tt.mm[1][0]) printf (" p_%s ", tt.mm[1]) ;
  if (tt.mm[2][0]) printf (" q_%s ", tt.mm[2]) ;
  if (tt.mm[3][0]) printf (" r_%s ", tt.mm[3]) ;

  i = 0 ;
  if      (tt.denom[0]) printf (" %c k^%d ", i++ ? ' ' : '/', 2 * tt.denom[0]) ;
  if      (tt.denom[1]) printf (" %c (k+p)^%d ", i++ ? ' ' : '/', 2 * tt.denom[1]) ;
  if      (tt.denom[2]) printf (" %c (k+p+q)^%d ", i++ ? ' ' : '/', 2 * tt.denom[2]) ;
  if      (tt.denom[3]) printf (" %c (k+p+q+r)^%d ", i++ ? ' ' : '/', 2 * tt.denom[3]) ;

  return ;
} /* showTT */

/*******************************************************************************************/

static void showPol (POLYNOME pp)
{
  POLYNOME p1, p2 ;
  static int nn = 0, level = 0 ;

   if (!pp)
     {
       if (level == 0) printf ("(NULL) #### zero term ###\n") ;
       return ;
     }
   if (! level) nn = 0 ;
   level++ ;
   p1 = pp->p1 ;
   p2 = pp->p2 ;
   if (p1 == pp) messcrash ("pp == pp->p1 in showPol") ;
   if (p2 == pp) messcrash ("pp == pp->p2 in showPol") ;
   if (p1 && p1 == p2) messcrash ("pp->p1 == pp->p2 in showPol") ;

   if (pp->isProduct)
     { printf ("(");  showPol (p1) ; printf (")*(");  showPol(p2) ; printf(")") ; }
   else if (pp->isSum)
     { printf ("(");  showPol (p1) ; printf (")+") ; if (p2 && !p2->isSum) printf("(");  showPol(p2) ; if (p2 && !p2->isSum) printf(")") ; }
   else if (pp->tt.type)
     { showTT (pp) ; nn++ ; }
   level-- ;

   if (! level) printf ("  ### %d term%s ###\n", nn, nn > 1 ? "s" : "") ;
} /* showPol */

/***********************************************************************************************************************************************/
/******************************************* New Polynomes of all types ************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME newPolynome (void) 
{
  static int id = 0 ;
  int n = sizeof (struct polynomeStruct) ;
  POLYNOME pp = (POLYNOME) halloc (n, 0) ;
  memset (pp, 0, n) ;
  pp->id = ++id ;
  return pp ;
}

/*************************************************************************************************/

static POLYNOME copyPolynome (POLYNOME p1)
{
  POLYNOME pp = 0 ; 
 
  if (p1)
    {
      int id ;
      pp = newPolynome () ;
      id = pp->id ;
      memcpy (pp, p1, sizeof (struct polynomeStruct)) ;
      pp->id = id ;
      if (pp->p1)
	pp->p1 =  copyPolynome (pp->p1) ;
      if (pp->p2)
	pp->p2 =  copyPolynome (pp->p2) ;
    }
  return pp ;
}

/*************************************************************************************************/

static POLYNOME newScalar (complex float z)
{
  POLYNOME p = newPolynome () ;
  p->tt.type = 1 ;
  p->tt.z = z ;
  return p ;
}

/*************************************************************************************************/

static POLYNOME newG (char mu, char nu)
{
  POLYNOME p = newPolynome () ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.g[0] = mu ;
  p->tt.g[1] = nu ;
  return p ;
}

/*************************************************************************************************/
/* newAG (... ,  0) antisymmetric link 1/2(ac bd - ad bc). Optionally adding the i epsilon 
 * newAG (... , +1) is the self-dual projector P+ = 1/4 ( ac bd - ad bc) + i/2 epsilon(abcd)
 * newAG (... , -1) is the self-dual projector P- = 1/4 ( ac bd - ad bc) - i/2 epsilon(abcd)
 *     WE have (P+)^2 = (P+),   (P-)^2 = (P-),   (P+)(P-) = (P-)(P+) = 0 
 */
static POLYNOME newEpsilon (char a, char b, char c, char d)
{
  POLYNOME pp ;
  
  pp = newPolynome () ;
  pp->tt.type = 1 ;
  pp->tt.z = 1.0 ;
  pp->tt.eps[0] = a ;
  pp->tt.eps[1] = b ;
  pp->tt.eps[2] = c ;
  pp->tt.eps[3] = d ;
  return pp ;
}

static POLYNOME newAG (char a, char b, char c, char d, int parity)
{
  POLYNOME pp, ppp[4] ;

  pp = ppp[0] = newPolynome () ;
  pp->tt.type = 1 ;
  pp->tt.z = 1.0/2.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = c ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = d ;

  pp = ppp[1] = newPolynome () ;
  pp->tt.type = 1 ;
  pp->tt.z = -1.0/2.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = d ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = c ;

  ppp[2] = ppp[3] = 0 ;

  if (parity)
    {
      ppp[0]->tt.z /= 2 ;
      ppp[1]->tt.z /= 2 ;

      pp = ppp[2] = newEpsilon (a, b, c, d) ;
      pp->tt.z = (I * parity ) / 4.0 ;
      if (parity == 2)
	return pp ;
    }
  
  pp = newMultiSum (ppp) ;
  return pp ; 
} /* newAG */

/*************************************************************************************************/

static POLYNOME newAG6 (char a, char b, char c, char d, char e, char f)
{
  POLYNOME pp, ppp[7] ;

  pp = ppp[0] = newPolynome () ;
  pp->tt.type = 1 ;
  pp->tt.z = 1.0/1.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = d ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = e ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = f ;

  pp = ppp[1] = newPolynome () ;
  pp->tt.type = 1 ;
  pp->tt.z = -1.0/1.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = d ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = f ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = e ;

  pp = ppp[2] = newPolynome () ;
  pp->tt.type = 1 ;
  pp->tt.z = -1.0/1.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = e ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = d ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = f ;

  pp = ppp[3] = newPolynome () ;
  pp->tt.type = 1 ;
  pp->tt.z = 1.0/1.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = e ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = f ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = d ;

  pp = ppp[4] = newPolynome () ;
  pp->tt.type = 1 ;
  pp->tt.z = -1.0/1.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = f ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = d ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = e ;

  pp = ppp[5] = newPolynome () ;
  pp->tt.type = 1 ;
  pp->tt.z = 1.0/1.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = f ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = e ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = d ;

  ppp[6] = 0 ;
  
  pp = newMultiSum (ppp) ;
  return pp ; 
} /* newAG6 */

/*************************************************************************************************/

static POLYNOME newK (char cc)
{
  POLYNOME p = newPolynome () ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[0][0] = cc ;
  return p ;
}

/*************************************************************************************************/

static POLYNOME newP (char cc)
{
  POLYNOME p = newPolynome () ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[1][0] = cc ;
  return p ;
}

/*************************************************************************************************/

static POLYNOME newQ (char cc)
{
  POLYNOME p = newPolynome () ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[2][0] = cc ;
  return p ;
}

/*************************************************************************************************/

static POLYNOME newR (char cc)
{
  POLYNOME p = newPolynome () ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[3][0] = cc ;
  return p ;
}

/*************************************************************************************************/

static POLYNOME newPQR (int pqr, char mu)
{
  POLYNOME p0 = pqr >=0 ? newK (mu) : 0 ;
  POLYNOME p1 = pqr >=1 ? newP (mu) : 0 ;
  POLYNOME p2 = pqr >=2 ? newQ (mu) : 0 ;
  POLYNOME p3 = pqr >=3 ? newR (mu) : 0 ;
  POLYNOME ppp[5] = {p0, p1, p2, p3, 0} ;

  return newMultiSum (ppp) ;
}

/*************************************************************************************************/

static POLYNOME newSigma (char cc)
{
  POLYNOME p = newPolynome () ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.Id2 = 1 ;
  p->tt.sigma[0] = cc ;
  return p ;
}

/*************************************************************************************************/

static POLYNOME newSigB (char cc)
{
  POLYNOME p = newPolynome () ;
  p->tt.type = 1 ;
  p->tt.z = 1 ;
  p->tt.Id2 = 1 ;
  p->tt.sigB[0] = cc ;
  return p ;
}

/*************************************************************************************************/

static POLYNOME newSum (POLYNOME p1, POLYNOME p2)
{
  if (p1 && p2)
    {
      POLYNOME p = newPolynome () ;
      p->p1 = p1 ;
      p->p2 = p2 ; 
      if (p1->tt.type && p2->tt.type && p1->tt.Id2 + p2->tt.Id2 == 1)
	messcrash ("Cannot add aPauli matrix to a number\n") ;;
      p->isSum = TRUE ;
      if (p1 == p2)
	messcrash ("p1 == p2 in newSum") ;
      return p ;
    }
  else if (p1)
    return copyPolynome (p1) ;
  else if (p2)
    return copyPolynome (p2) ;
  return 0 ;
}

/*************************************************************************************************/

static POLYNOME newProduct (POLYNOME p1, POLYNOME p2)
{
  if (p1 && p2)
    {
      POLYNOME p = newPolynome () ;
      p->p1 = p1 ;
      p->p2 = p2 ;
      p->isProduct = TRUE ;
      return p ;
    }
  return 0 ;
}

/*************************************************************************************************/

static POLYNOME newMultiSum (POLYNOME ppp[])
{
  POLYNOME pp, p1, p2 ;
  int i = 0 ;

  while (ppp[i]) i++ ;
  if (i <= 1) return copyPolynome(ppp[0]) ;

  pp = ppp[--i] ;
  while (i > 0)
    {
      p2 = pp ;
      p1 = ppp[--i] ;
      pp = newSum (p1, p2) ;
    }
  return pp ;
}

/*************************************************************************************************/

static POLYNOME newMultiProduct (POLYNOME ppp[])
{
  POLYNOME pp, p1, p2 ;
  int i = 0 ;

  while (ppp[i]) i++ ;
  if (i <= 1) return copyPolynome (ppp[0]) ;

  pp = ppp[--i] ;
  while (i > 0)
    {
      p2 = pp ;
      p1 = ppp[--i] ;
      pp = newProduct (p1, p2) ;
    }
  return pp ;
}

/***********************************************************************************************************************************************/
/**************************************** Quantum Field theory rules ***************************************************************************/
/***********************************************************************************************************************************************/
/* k/sig to g */
static int indexTrace (char *old, char *new, int sign, AC_HANDLE h)
{
  /*
  char *cp, buf[GMAX] ;
  int i, j, k ;
  static int level = 0 ;
  Array a = arrayNandleCreate () ;
  level++ ;
  for (i = 1 ; i < GMAX && old[i] ; i++) // which index shall i contract with index 0 
    {
      for (k = 0, j = 1 ; j < GMAX && old[j])
	if (i != j)
	  buf[k++] = old[j] ;
      new[0] = old[0] ;
      new[1] = old[i] ;
      indexTrace (buf, new, sign) ;
    }
  level-- ;

*/
  return 0 ;
}

/***********************************************************************************************************************************************/

static int indexTtSort (char *cp, int dx, int sign)
{
  int i, j, k, ss = 1, blockS = 1 ;
  BOOL modif = TRUE ;

  for (i = 0 ; i < dx ; i++)
    blockS *= sign ;
  while (modif)
    {   /* sort inside the goups of length dx (i.e. g_ba -> g_ab  eps_acbd -> - eps_abcd */
      modif = FALSE ;
      for (i = 0 ; i < GMAX + dx - 1 && cp[i + dx - 1] ; i += dx)
	for (j = 0 ; j < dx -1 ; j++)
	  if (cp[i+j+1]  && cp[i+j] > cp[i+j+1])
	    { k = cp[i+j] ; cp[i+j] = cp[i+j+1] ; cp[i+j+1] = k ; ss *= sign ; modif = TRUE ; }
    }
  modif = TRUE ;
  while (modif)
    {   /* sort inside the blocks */
      modif = FALSE ;
      for (i = 0 ; i < GMAX - dx && cp[i+dx] ; i += dx)
	if (cp[i+dx]  && cp[i] > cp[i+dx])
	for (j = 0 ; j < dx ; j++)
	    { k = cp[i+j] ; cp[i+j] = cp[i+j+dx] ; cp[i+j+dx] = k ; ss *= blockS ; modif = TRUE ; }
    }
  return ss ;
} /* indextTtSort */

/***********************************************************************************************************************************************/
/* Einstein contraction rules */
static TT contractTtIndices (POLYNOME pp)
{
  complex float zz = 1 ;
  TT tt = pp->tt ;

  /* sort and search repeated pair of indices inside the metric itself */
  if (tt.type)
    {
      int ii, i, j, k = 0 ;
      char *g, *s ;
      BOOL ok = FALSE ;

      while (! ok)
	{
	  ok = TRUE ;
	  /* sort the indices */
	  for (i = 0 ; i < 4 ; i++)
	    tt.z *= indexTtSort (tt.mm[i], 1, 1) ;
	  tt.z *= indexTtSort (tt.g, 2, 1) ;
	  tt.z *= indexTtSort (tt.eps, 4, -1) ;
	  /* simplify repeated k indices */
	  for (i = 0 ; i < GMAX -1 ; i++)
	    if (tt.mm[0][i] && tt.mm[0][i] == tt.mm[0][i+1] && tt.denom[0]) 
	      {
		tt.denom[0]-- ;   /* divide by k^2 top and bottom */
		for (j = 0 ; i+j < GMAX -2 ; j++)
		  tt.mm[0][i+j] = tt.mm[0][i+j+2] ;
		tt.mm[0][GMAX-2] = 0 ; 
		tt.mm[0][GMAX-1] = 0 ;
		i-- ;  /* scan again the same position */
	      }

	  /* search repeated indices in a single metric */ 
	  for (i = 0, g = tt.g ; g[i] && i < GMAX ; i+=2)
	    { 
	      if (g[i] == g[i+1])
		{ zz *= 4 ; ok = FALSE ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
	      if (g[i] > g[i+1]) /* switch : the Lorentz metric is Abelian */
		{ char cc = g[i] ; g[i] = g[i+1] ; g[i+1] = cc ; ok = FALSE ; }
	    }
	  
	  /* search repeated indices in a pair of metrics, do one modif at a time */
	  for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	    for (j = i + 2 ; ok && g[j] && j < GMAX ; j+= 2)
	      { 
		if (g[i] == g[j])
		  { g[i] = g[j+1] ; for (k = j ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
		else if (g[i] == g[j+1])
		  { g[i] = g[j] ; for (k = j ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
		else if (g[i+1] == g[j])
		  { g[i+1] = g[j+1] ; for (k = j ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
		else if (g[i+1] == g[j+1])
		  { g[i+1] = g[j] ; for (k = j ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
	      }
	  
	  /* search repeated indices between a metric and epsilon */
	  for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	    for (j = 0, s = tt.eps ; ok && s[j] && j < GMAX ; j++)
	      { 
		if (g[i] == s[j])
		  { s[j] = g[i+1] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k] = g[k+1] = 0 ; ok = FALSE ; }
		if (g[i+1] == s[j])
		  { s[j] = g[i] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k] = g[k+1] = 0 ; ok = FALSE ; }
	      }
	
	  /* search repeated indices inside an epsilon */
	  for (i = 0, s = tt.eps ; ok && s[i] && i < GMAX ; i+= 4)
	    for (j = 0 ; j <= 2 ; j++)
	      for (k = j+1 ; k <= 3 ; k++)
		{
		  if (s[i+j] == s[i+k])
		    {
		      tt.z = 0 ;
		      return tt ;
		    }
		}
	
	  /* search repeated indices between epsilon and momenta */
	  for (ii = 0 ; ii < 3 ; ii++)
	    for (i = 0, s = tt.eps ; ok && s[i] && i < GMAX ; i+= 4)
	      for (j = 0 ; j <= 3 ; j++)
		{
		  int k, l, m ; char *t, *u ;
		  for (k = 0, t = tt.mm[ii] ; ok && t[k] && k < GMAX ; k++)		
		    if (s[i+j] == t[k])
		      for (l = 0, u = tt.mm[ii] ; ok && s[l] && l < GMAX ; l++)		
			for (m = 0 ; m <= 3 ; m++)
			  if (k != l && s[i+m] == u[l])
			    {
			      tt.z = 0 ;
			      return tt ;
			    }
		}
      
	  /* search repeated indices between a metric and a sigma */
	  for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	    for (j = 0, s = tt.sigma ; ok && s[j] && j < GMAX ; j++)
	      { 
		if (g[i] == s[j])
		  { s[j] = g[i+1] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
		if (g[i+1] == s[j])
		  { s[j] = g[i] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
	      }
	
	  /* search repeated indices between a metric and a sigma-bar */
	  for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	    for (j = 0, s = tt.sigB ; ok && s[j] && j < GMAX ; j++)
	      { 
		if (g[i] == s[j])
		  { s[j] = g[i+1] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
		if (g[i+1] == s[j])
		  { s[j] = g[i] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; ok = FALSE ; }
	      }
	    
	  /* search repeated indices between a metric and a momentum */
	  for (ii = 0 ; ii < 3 ; ii++)
	    for (i = 0, g = tt.g ; ok && g[i] && i < GMAX - 2 ; i+=2)
	      for (j = 0, s = tt.mm[ii] ; ok && s[j] && j < GMAX ; j++)
		{ 
		  if (g[i] == s[j])
		    { s[j] = g[i+1] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k] = g[k+1] = 0 ; ok = FALSE ; }
		  if (g[i+1] == s[j])
		    { s[j] = g[i] ; for (k = i ; k < GMAX - 2 ; k++) g[k] = g[k+2] ; g[k] = g[k+1] = 0 ; ok = FALSE ; }
		}
	  
	  /* sort alphabetically the momenta, they are Abelian */
	  for (ii = 0 ; ii < 3 ; ii++)
	    for (i = 0, s = tt.mm[ii] ; ok && s[i] && i < GMAX ; i++)
	      for (j = i + 1 ; s[j] && j < GMAX ; j++)
		if (s[j] < s[i])
		  { char cc = s[i] ; s[i] = s[j] ; s[j] = cc ; }
	  
	  /* sort internal dummy k-slash indices, rename them */
	  if (0) 
	    {
	      int nd = 0 ; char *u, *v ;
	      for (ii = 0 ; ii < 3 ; ii++)
		for (i = 0, u = tt.mm[ii] ; ok && u[i] && i < GMAX ; i++)
		  {
		    for (j = 0, v = tt.sigma ; ok && v[j] && j < GMAX ; j++)
		      if (u[i] < 'w' && u[i] == v[j])
			{ u[i] = v[j] = 'w' + nd++ ; ok = FALSE ; }
		    for (j = 0, v = tt.sigB ; ok && v[j] && j < GMAX ; j++)
		      if (u[i] < 'w' && u[i] == v[j])
			{ u[i] = v[j] = 'w' + nd++ ; ok = FALSE ; }
		  }
	    }
	  
	  /* search contiguous sigma_a sigB_a = 4, AND  s_a s_b s_a = -2 s_b AND  abEab = 4 E AND s_a sB_b s_c sB_a = 4 g_bc */
	  for (i = 0, s = tt.sigma, g = tt.g ; ok && s[i] && i < GMAX - 3 ; i++)
	    { 
	      if (s[i] == s[i+1])
		{ zz *= 4 ; for (k = i ; k < GMAX - 2 ; k++) s[k] = s[k+2] ; s[k+1] = s[k+2] = 0 ; ok = FALSE ; }
	      else if (s[i] == s[i+2])
		{ zz *= -2 ; s[i] = s[i+1] ;for (k = i + 1 ; k < GMAX - 2 ; k++) s[k] = s[k+2] ; s[k+1] = s[k+2] = 0 ; ok = FALSE ; }
	      else if (s[i] == s[i+3] && s[i+1] == s[i+4]) 
		{ zz *= 4 ; s[i] = s[i+2] ; for (k = i + 1 ; k < GMAX - 4 ; k++) s[k] = s[k+4] ; for (; k < GMAX ; k++) s[k] = 0 ; ok = FALSE ; }
	      else if (s[i] == s[i+4] && s[i+1] == s[i+5]  && s[i+2] == s[i+6])
		{ zz *= -32 ; s[i] = s[i+3] ; for (k = i + 1 ; k < GMAX - 6 ; k++) s[k] = s[k+6] ; for (; k < GMAX ; k++) s[k] = 0 ; ok = FALSE ; }
	      else if (s[i] == s[i+3])
		{ zz *= 4 ; k = strlen (g) ; g[k] = s[i+1]; g[k+1] = s[i+2] ; g[k+2] = 0 ; for (k = i ; k < GMAX - 4 ; k++) s[k] = s[k+4] ; for (; k < GMAX ; k++) s[k] = 0 ; ok = FALSE ; }
	    }
	  /* search contiguous sigB_a sigma_a = 4, AND  sb_a s_b sB_a = -2 s_b AND sB_a s_b sB_c s_a = 4 g_bc */
	  for (i = 0, s = tt.sigB, g = tt.g ; ok && s[i] && i < GMAX - 3 ; i++)
	    { 
	      if (s[i] == s[i+1])
		{ zz *= 4 ; for (k = i ; k < GMAX - 2 ; k++) s[k] = s[k+2] ; s[k+1] = s[k+2] = 0 ; ok = FALSE ; }
	      else if (s[i] == s[i+2])
		{ zz *= -2 ; s[i] = s[i+1] ;for (k = i + 1 ; k < GMAX - 2 ; k++) s[k] = s[k+2] ; s[k+1] = s[k+2] = 0 ; ok = FALSE ; }
	      else if (s[i] == s[i+3] && s[i+1] == s[i+4]) 
		{ zz *= 4 ; s[i] = s[i+2] ; for (k = i + 1 ; k < GMAX - 4 ; k++) s[k] = s[k+4] ; for (; k < GMAX ; k++) s[k] = 0 ; ok = FALSE ; }
	      else if (s[i] == s[i+4] && s[i+1] == s[i+5]  && s[i+2] == s[i+6])
		{ zz *= -32 ; s[i] = s[i+3] ; for (k = i + 1 ; k < GMAX - 6 ; k++) s[k] = s[k+6] ; for (; k < GMAX ; k++) s[k] = 0 ; ok = FALSE ; }
	      else if (s[i] == s[i+3])
		{ zz *= 4 ; k = strlen (g) ; g[k] = s[i+1]; g[k+1] = s[i+2] ; g[k+2] = 0 ; for (k = i ; k < GMAX - 4 ; k++) s[k] = s[k+4] ; for (; k < GMAX ; k++) s[k] = 0 ; ok = FALSE ; }
	    }
	}
    }
  tt.z *= zz ;
  return tt ;
}

static POLYNOME contractIndices (POLYNOME pp)
{
  POLYNOME p1, p2 ;

  if (!pp)
    return 0 ;
  p1 = pp->p1 ;
  p2 = pp->p2 ;
  p1 = contractIndices (p1) ;
  p2 = contractIndices (p2) ;

  if (pp->tt.type)
    {
      pp->tt = contractTtIndices (pp) ;
      if (pp->tt.z == 0)
	return 0 ;
    }
  return pp ;
}

/* in a product of monomes, the list of symbols must merge */
static int contractTTProducts (POLYNOME pp, POLYNOME p1, POLYNOME p2)
{
  char *u, *v, *w ;
  char buf[GMAX] ;
  int ii, i, j ;
  TT tt = pp->tt ;
  TT t1 = p1->tt ;
  TT t2 = p2->tt ;

  /* merge numbers */
  tt.z = t1.z * t2.z ;
  if (cabs (tt.z) < minAbs)
    {
      if (0)
	{
	  ac_free (p1) ;
	  ac_free (p2) ;
	}
      pp->isFlat = FALSE ;
      return 0 ;
    }

  /* merge denoms */
  for (i = 0 ; i < 4 ; i++)
    {
      if (tt.denom[i] != t1.denom[i] + t2.denom[i])
	pp->isFlat = FALSE ;
      tt.denom[i] = t1.denom[i] + t2.denom[i] ;
    }

  /* merge metrics */
  u = t1.g ; v = t2.g ; w = tt.g ;
  memcpy (buf, w, GMAX) ;
  i = strlen (u) ; j = strlen (v) ;
  if (i+j)
    {
      if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
      while ((*w++ = *u++)) ;  
      w-- ; 
      while ((*w++ = *v++)) ;
    }
  if (memcmp (buf, w, GMAX))
    pp->isFlat = FALSE ;

  u = t1.gg ; v = t2.gg ; w = tt.gg ;
  i = strlen (u) ; j = strlen (v) ;
  memcpy (buf, w, GMAX) ;
  if (i+j)
    {
      if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
      while ((*w++ = *u++)) ;
      w-- ; while ((*w++ = *v++)) ;
    }
  if (memcmp (buf, w, GMAX))
    pp->isFlat = FALSE ;

  /* merge epsilon */
  u = t1.eps ; v = t2.eps ; w = tt.eps ;
  memcpy (buf, w, GMAX) ;
  i = strlen (u) ; j = strlen (v) ;
  if (i+j)
    {
      if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
      while ((*w++ = *u++)) ;  
      w-- ; 
      while ((*w++ = *v++)) ;
    }
  if (memcmp (buf, w, GMAX))
    pp->isFlat = FALSE ;

  /* merge momenta */
  for (ii = 0 ; ii < 4 ; ii++)
    {
      memcpy (buf, w, GMAX) ;
      u = t1.mm[ii] ; v = t2.mm[ii] ; w = tt.mm[ii] ;
      i = strlen (u) ; j = strlen (v) ;
      if (i+j)
	{
	  memset (w, 0, GMAX) ;
	  if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
	  while ((*w++ = *u++)) ;
	  w-- ; while ((*w++ = *v++)) ;
	}
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;
    }

  /* merge matrices */
  if (t1.Id2 || t2.Id2)
    tt.Id2 = 1 ;
  if (t1.sigma[0] && t1.sigB[0])
    messcrash ("cannot have a sigma sigma product, should be sigma sigma->bar") ;
  else if (t2.sigma[0] && t2.sigB[0])
    messcrash ("cannot have a sigma sigma product, should be sigma sigma->bar") ;
  else if (! t1.sigma[0] && ! t1.sigB[0])
    { 
      w = tt.sigma ; v = t2.sigma ; 
      memcpy (buf, w, GMAX) ;
      while ((*w++ = *v++)) ; 
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;

      w = tt.sigB ;  v = t2.sigB  ; 
      memcpy (buf, w, GMAX) ;
      while ((*w++ = *v++)) ; 
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;
    }
  else if (! t2.sigma[0] && ! t2.sigB[0])
    {
      w = tt.sigma ; u = t1.sigma ;
      memcpy (buf, w, GMAX) ;
      while ((*w++ = *u++)) ; 
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;

      w = tt.sigB ;  u = t1.sigB  ;
      memcpy (buf, w, GMAX) ;
      while ((*w++ = *u++)) ; 
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;
    }
  else if (t1.sigma[0] && (strlen(t1.sigma) % 2) && t2.sigma[0])
    messcrash ("cannot have a sigma sigma product, should be sigma sigma->bar") ;
  else if (t1.sigma[0] && (strlen(t1.sigma) % 2 == 0) && t2.sigB[0])
    messcrash ("cannot have a sigma sigma product, should be sigma sigma->bar") ;
  else if (t1.sigma[0])
    { 
      u = t1.sigma ; i = strlen (u) ; w = tt.sigma ;
      memcpy (buf, w, GMAX) ;
      if (i % 2) 
	{
	  v = t2.sigB ; j = strlen (v) ;
	  if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
	  while ((*w++ = *u++)) ;  
	  w-- ; 
	  while ((*w++ = *v++)) ;
	}
      else
	{
	  v = t2.sigma; j = strlen (v) ;
	  if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
	  while ((*w++ = *u++)) ;
	  w-- ; while ((*w++ = *v++)) ;
	}
      if (!tt.sigma[0])
	tt.Id2 = 1 ;
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;
    }
  else if (t1.sigB[0])
    {
      u = t1.sigB ; i = strlen (u) ; w = tt.sigB ;
      memcpy (buf, w, GMAX) ;
      if (i % 2 == 0) 
	{
	  v = t2.sigB ; j = strlen (v) ;
	  if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
	  while ((*w++ = *u++)) ;
	  w-- ; while ((*w++ = *v++)) ;
	}
      else 
	{
	  v = t2.sigma; j = strlen (v) ;
	  if (i+j >= GMAX) messcrash ("Merging too long metrics: %s %s", u, v) ;
	  while ((*w++ = *u++)) ;
	  w-- ; while ((*w++ = *v++)) ;
	}
      if (!tt.sigB[0])
	tt.Id2 = 1 ;
      if (memcmp (buf, w, GMAX))
	pp->isFlat = FALSE ;
    }

  tt.type = 1 ;
  pp->tt = tt ;
  return tt.type ;
}

/*******************************************************************************************/
/* transform X eps_abcd eps_cdef into (X * (g_ae g_bf - g_af g_be) */
static POLYNOME contractEpsilon (POLYNOME pp)
{
  char *cp = 0 ;
  int pass, i1, j1, i, j, n, sign = 1 ;
  
 redo:
  contractIndices (pp) ;
  if (!pp)
    return 0 ;
  if (pp->tt.type && cabs (pp->tt.z) < minAbs)
    return 0 ;
  cp = pp->tt.eps ;
  for (pass = 0 ; pass < 2 ; pass++)
    {
      for (i1 = 0 ; cp[i1] ; i1+= 4)
	for (j1 = i1+4 ; cp[j1] ; j1+= 4)
	  {
	    /* compare the pair i1,j1 of epsilon symbols */
	    for (n = 0, i = i1 ; i < i1 + 4 ; i++)
	      for (j = j1 ; j < j1 + 4 ; j++)
		if (cp[i] == cp[j]) 
		  {
		    n++ ; 
		    if ((j-i)% 2) sign = -sign ;
		  }
	    switch (n)
	      {          /* do the easy cases */
	      case 0:
		break ;
	      case 4:    /* eps^2 = -24, minus in minkovski, plus in euclidean */
		pp->tt.z *= -24 * sign ;
		pp->isFlat = FALSE ;
		/* cleanup the epsilons */
		for (i = j1 ; i < GMAX - 4 ; i++)
		  cp[i] = cp[i+4] ;
		for (i = i ; i < GMAX ; i++)
		  cp[i] = 0 ;
		for (i = i1 ; i < GMAX - 4 ; i++)
		  cp[i] = cp[i+4] ;
		for (i = i ; i < GMAX ; i++)
		  cp[i] = 0 ;
		goto redo ;
		break ;
	      case 3:    /*  - 6 * g_mu_nu up to a sign */
		pp->isFlat = FALSE ;
		pp->tt.z *= -6 * sign ;
		/* locate the surviving index */
		for (i = i1 ; i < i1+4 ; i++)
		  {
		    int ok = 0 ;
		    for (j = j1 ; j < j1 + 4 ; j++)
		      if (cp[i] == cp[j])
			ok = 1 ;

		    if (!ok)
		      {
			char *cq = pp->tt.g ;
			cq += strlen (cq) ;
			cq[0] = cp[i] ;
			cq[1] = 0 ;
		      }
		  }
		for (i = j1 ; i < j1+4 ; i++)
		  {
		    int ok = 0 ;
		    for (j = i1 ; j < i1 + 4 ; j++)
		      if (cp[i] == cp[j])
			ok = 1 ;

		    if (!ok)
		      {
			char *cq = pp->tt.g ;
			cq += strlen (cq) ;
			cq[0] = cp[i] ;
			cq[1] = 0 ;
		      }
		  }
		/* cleanup the epsilons */
		for (i = j1 ; i < GMAX - 4 ; i++)
		  cp[i] = cp[i+4] ;
		for (i = i ; i < GMAX ; i++)
		  cp[i] = 0 ;
		for (i = i1 ; i < GMAX - 4 ; i++)
		  cp[i] = cp[i+4] ;
		for (i = i ; i < GMAX ; i++)
		  cp[i] = 0 ;
		goto redo ;
		break ;
	      }
	    if (pass < 1)
	      continue ;
	    switch (n)
	      {
	      case 2:
		/* transform X eps_abcd eps_cdef into (X * (g_ae g_bf - g_af g_be) */
		if (1)
		  {
		    POLYNOME ppp[4] ;
		    char a[4] ;
		    int ia = 0 ;
		    memset (a, 0, sizeof(a)) ;		    
		    /* locate the 4 surviving index */
		    for (i = i1 ; i < i1+4 ; i++)
		      {
			int ok = 1 ;
			for (j = j1 ; j < j1 + 4 ; j++)
			  if (cp[i] == cp[j])
			    ok = 0 ;
			if (ok)
			  a[ia++] = cp[i] ;
		      }
		    for (i = j1 ; i < j1+4 ; i++)
		      {
			int ok = 1 ;
			for (j = i1 ; j < i1 + 4 ; j++)
			  if (cp[i] == cp[j])
			    ok = 0 ;
			if (ok)
			  a[ia++] = cp[i] ;
		      }
		    
		    ppp[0] = pp ;
		    pp->tt.z *= (-4 * sign) ;
		    ppp[1] = newAG (a[0],a[1],a[2],a[3], 0) ;
		    ppp[2] = 0 ;
		    /* cleanup the epsilons */
		    for (i = j1 ; i < GMAX - 4 ; i++)
		      cp[i] = cp[i+4] ;
		    for (i = i ; i < GMAX ; i++)
		      cp[i] = 0 ;
		    for (i = i1 ; i < GMAX - 4 ; i++)
		      cp[i] = cp[i+4] ;
		    for (i = i ; i < GMAX ; i++)
		      cp[i] = 0 ;
		    pp = newMultiProduct (ppp) ;
		    pp->isFlat = FALSE ;
		    pp->p1->isFlat = FALSE ;
		    pp->p2->isFlat = FALSE ;
		    goto redo ;
		  }
		break ;
	      case 1:
		/* transform X eps_abcd eps_adef into (X * (g_bd g_ce g_df ... 3 terns */
		if (1)
		  {
		    POLYNOME ppp[4] ;
		    char a[6] ;
		    int ia = 0 ;
		    memset (a, 0, sizeof(a)) ;		    
		    /* locate the 6 surviving index */
		    for (i = i1 ; i < i1+4 ; i++)
		      {
			int ok = 1 ;
			for (j = j1 ; j < j1 + 4 ; j++)
			  if (cp[i] == cp[j])
			    ok = 0 ;
			if (ok)
			  a[ia++] = cp[i] ;
		      }
		    for (i = j1 ; i < j1+4 ; i++)
		      {
			int ok = 1 ;
			for (j = i1 ; j < i1 + 4 ; j++)
			  if (cp[i] == cp[j])
			    ok = 0 ;
			if (ok)
			  a[ia++] = cp[i] ;
		      }
		    
		    ppp[0] = pp ;
		    pp->tt.z *= (- 0.5 *sign) ;
		    ppp[1] = newAG6 (a[0],a[1],a[2],a[3],a[4],a[5]) ;
		    ppp[2] = 0 ;

		    /* cleanup the epsilons */
		    /* cleanup the epsilons */
		    for (i = j1 ; i < GMAX - 4 ; i++)
		      cp[i] = cp[i+4] ;
		    for (i = i ; i < GMAX ; i++)
		      cp[i] = 0 ;
		    for (i = i1 ; i < GMAX - 4 ; i++)
		      cp[i] = cp[i+4] ;
		    for (i = i ; i < GMAX ; i++)
		      cp[i] = 0 ;
		    pp = newMultiProduct (ppp) ;
		    pp->isFlat = FALSE ;
		    pp->p1->isFlat = FALSE ;
		    pp->p2->isFlat = FALSE ;
		    goto redo ;
		  }
		break ;
	      }
	  }
    }
  return pp ;
} /* contractEpsilon */
  
/*******************************************************************************************/

static POLYNOME contractProducts (POLYNOME pp)
{
  POLYNOME p1, p2 ;
  BOOL debug = FALSE ;
  static int nn= 0 ;
  if (!pp)
    return 0 ;
  nn++ ;
  p1 = pp->p1 ;
  p2 = pp->p2 ;
  p2 = pp->p2 = contractProducts (p2) ;
  p1 = pp->p1 = contractProducts (p1) ;

  if (debug)  showPol(pp) ;
  if (0 && pp->isProduct && p1 && p2 && p1->tt.type && p2->isProduct && p2->p1 && p2->p1->tt.type)
      contractTTProducts (pp, p1, p2->p1) ;


  if (pp->isSum && p1 && !p2)
    {
      *pp = *p1 ;
      pp->isFlat = FALSE ;
    }

  if (pp->isSum && p2 && !p1)
    {
      *pp = *p2 ;
      pp->isFlat = FALSE ;
    }

  if (pp->isSum && p1 && p2 && ! p1->tt.type && p2->tt.type)
    { /* addition is Abelian */
      POLYNOME q = pp->p1 ;
      pp->p1 = pp->p2 ;
      pp->p2 = q ;
      pp->isFlat = FALSE ;
    }
  
  if (pp->isSum && p1 && p2 && p1->tt.type && p2->tt.type)
    {
      float complex z1 = p1->tt.z ;
      float complex z2 = p2->tt.z ;
      int s ;
      s = polOrder (&p1, &p2) ;

      if (s == 0)
	{
	  *pp = *p1 ;
	  pp->isFlat = FALSE ;
	  pp->tt.z = z1 + z2 ;
	}
      else if (s > 0)
	{ /* addition is Abelian */
	  POLYNOME q = pp->p1 ;
	  pp->p1 = pp->p2 ;
	  pp->p2 = q ;
	  pp->isFlat =  FALSE ;
	}
    }
  

  if (pp->isSum && p1 && p2 && p1->tt.type && p2->isSum && p2->p1 && p2->p1->tt.type)
    {
      float complex z1 = p1->tt.z ;
      float complex z2 = p2->p1->tt.z ;
      int s ;
      p1->tt.z = 0 ;       
      p2->p1->tt.z = 0 ;       
      s = memcmp (&(p1->tt), &(p2->p1->tt), sizeof(TT)) ;
      p1->tt.z = z1 ;       
      p2->p1->tt.z = z2 ;       

      if (s == 0)
	{
	  *pp = *p2 ;
	  pp->isFlat = FALSE ;
	  pp->p1->tt.z = z1 + z2 ;
	}
      else if (s > 0)
	{ /* addition is Abelian */
	  POLYNOME q = pp->p1 ;
	  pp->p1 = pp->p2->p1 ;
	  pp->p2->p1 = q ;
	  pp->isFlat =  FALSE ;
	}
    }
  

  if (pp->isProduct && p1 && p2 && p1->tt.type && p2->tt.type)
    {
      pp->isProduct = FALSE ;
      pp->isFlat = TRUE ;
      contractTTProducts (pp, p1, p2) ;
      if (0)
	{
	  ac_free (p1) ;
	  ac_free (p2) ;
	}
      pp->p1 = pp->p2 = 0 ;
      if (debug) showPol(pp) ;
    }
  return pp ;
}
  
/*******************************************************************************************/
/*******************************************************************************************/
/* incomplet, this only works for pairs of sigma, we need the cases4,6,8 ...
 * which create polynomes in gg, not monomes 
 */
static POLYNOME pauliTraceTT (POLYNOME pp)
{
  TT tt = pp->tt ; 
  int i ;
  char *s = tt.sigma ; 
  char *sb = tt.sigB ;
  complex float parity = I ; 

  pp->isFlat = FALSE ;

  if (s[0] && sb[0])
    messcrash ("Computing the trace of a monome wwhere sigma=%s and sigB=% are both present\n", s, sb) ; 
  if (sb[0])
    { s = sb ; parity = -I ; }
  i = strlen (s) ;
  if (i % 2)
    { tt.z = 0 ; return 0 ; }
  if (i == 0)
    {
      if (tt.Id2)
	pp->tt.z *= 2 ;
    }
  else if (i == 2)
    { 
      char *g = tt.g ;
      int k = strlen (g) ;
      g += k ;
      while ((*g++ = *s++)) ;
      memset (tt.sigma , 0, GMAX) ;
      memset (tt.sigB , 0, GMAX) ;
      pp->tt = tt ;

      pp->tt.z *= 2 ; /* trace (identity) = 2 */
    }
  else if (i == 4)
    { 
      int n, N = 4, NN = 3 ;
      char S[N] ;
      memcpy (S, s, N) ;
      pp->tt.z *= 2 ; /* trace (identity) = 2 */
      char *gg = tt.g ;
      int k = strlen (gg) ;
      POLYNOME ppp[NN+1] ;
      char *z[3] = { "abcd", "acbd", "adbc"} ;
      for (n = 0 ; n < NN ; n++)
	{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
	  int i ;
	  ppp[n] = copyPolynome (pp) ;
	  memset (ppp[n]->tt.sigma , 0, GMAX) ;
	  memset (ppp[n]->tt.sigB , 0, GMAX) ;
	  if (n%2) ppp[n]->tt.z *= -1 ;   /* alternate signs */
	  for (i = 0 ; i < N ; i++)
	    {
	      ppp[n]->tt.g[k+i] = S[z[n][i] - 'a'] ;
	    }
	}
      if (0)
	{
	  ppp[n] = copyPolynome (pp) ;
	  memset (ppp[n]->tt.sigma , 0, GMAX) ;
	  memset (ppp[n]->tt.sigB , 0, GMAX) ;
	  memcpy (ppp[n]->tt.eps, S, 4) ; ppp[n]->tt.eps[4] = 0 ; 
	  ppp[n]->tt.z *= parity ;
	  n++ ;
	}
      ppp[n++] = 0 ; /* zero terminate the list */	
      pp = newMultiSum (ppp) ;
    }
  else if (i == 6)
    { 
      int n, N = 6, NN = 15 ;
      char S[N] ;
      memcpy (S, s, N) ;
      pp->tt.z *= 2 ; /* trace (identity) = 2 */
      memset (tt.sigma, 0, GMAX) ;
      memset (tt.sigB, 0, GMAX) ;
      char *gg = tt.g ;
      int k = strlen (gg) ;
      POLYNOME ppp[NN+1] ;
      char *z[15] = { "abcdef","abcedf","abcfde",
		      "acbdef","acbedf","acbfde",
		      "adbcef","adbecf","adbfce",
		      "aebcdf","aebdcf","aebfcd",
		      "afbcde","afbdce","afbecd"
      } ;
      /*
      char *e1[6] = { "abcdef","bcadef","cdabef","deabcf", "fabcde", 0 } ;
      char *e2[6] = { "acbdef","bdacef","ceabdf","dfabce", "eabcdf", 0 } ;
      */
      for (n = 0 ; n < NN ; n++)
	{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
	  int i ;
	  ppp[n] = copyPolynome (pp) ;
	  memset (ppp[n]->tt.sigma , 0, GMAX) ;
	  memset (ppp[n]->tt.sigB , 0, GMAX) ;
	  if (n%2) ppp[n]->tt.z *= -1 ;   /* alternate signs */
	  for (i = 0 ; i < N ; i++)
	    {
	      ppp[n]->tt.g[k+i] = S[z[n][i] - 'a'] ;
	    }
	}
      ppp[n] = 0 ; /* zero terminate the list */	
      pp = newMultiSum (ppp) ;
    }

  else if (i == -8)
    { 
      int n, N = 8, NN = 105 ;
      char S[N] ;
      memcpy (S, s, N) ;
      pp->tt.z *= 2 ; /* trace (identity) = 2 */
      memset (tt.sigma, 0, GMAX) ;
      memset (tt.sigB, 0, GMAX) ;
      char *gg = tt.g ;
      int k = strlen (gg) ;
      POLYNOME ppp[NN+1] ;
      char *z[105] = {
	"abcdefgh","abcdegfh","abcdehfg",
	"abcedfgh","abcedgfh","abcedhfg",
	"abcfdegh","abcfdgeh","abcfdheg",
	"abcgdefh","abcgdfeh","abcgdhef",
	"abchdefg","abchdfeg","abchdgef",
	
	"acbdefgh","acbdegfh","acbdehfg",
	"acbedfgh","acbedgfh","acbedhfg",
	"acbfdegh","acbfdgeh","acbfdheg",
	"acbgdefh","acbgdfeh","acbgdhef",
	"acbhdefg","acbhdfeg","acbhdgef",

 	"adcbefgh","adcbegfh","adcbehfg",
 	"adcebfgh","adcebgfh","adcebhfg",
 	"adcfbegh","adcfbgeh","adcfbheg",
 	"adcgbefh","adcgbfeh","adcgbhef",
 	"adchbefg","adchbfeg","adchbgef"

 	"aecdbfgh","aecdbgfh","aecdbhfg",
 	"aecbdfgh","aecbdgfh","aecbdhfg",
 	"aecfdbgh","aecfdgbh","aecfdhbg",
 	"aecgdbfh","aecgdfbh","aecgdhbf",
 	"aechdbfg","aechdfbg","aechdgbf"

 	"afcdebgh","afcdegbh","afcdehbg",
 	"afcedbgh","afcedgbh","afcedhbg",
 	"afcbdegh","afcbdgeh","afcbdheg",
 	"afcgdebh","afcgdbeh","afcgdheb",
 	"afchdebg","afchdbeg","afchdgeb"

 	"agcdefbh","agcdebfh","agcdehfb",
 	"agcedfbh","agcedbfh","agcedhfb",
 	"agcfdebh","agcfdbeh","agcfdheb",
 	"agcbdefh","agcbdfeh","agcbdhef",
 	"agchdefb","agchdfeb","agchdbef"

 	"ahcdefgb","ahcdegfb","ahcdebfg",
 	"ahcedfgb","ahcedgfb","ahcedbfg",
 	"ahcfdegb","ahcfdgeb","ahcfdbeg",
 	"ahcgdefb","ahcgdfeb","ahcgdbef",
 	"ahcbdefg","ahcbdfeg","ahcbdgef"
       } ;
      for (n = 0 ; n < NN ; n++)
	{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
	  int i ;
	  ppp[n] = copyPolynome (pp) ;
	  memset (ppp[n]->tt.sigma , 0, GMAX) ;
	  memset (ppp[n]->tt.sigB , 0, GMAX) ;
	  if (n%2) ppp[n]->tt.z *= -1 ;   /* alternate signs */
	  for (i = 0 ; i < N ; i++)
	    {
	      ppp[n]->tt.g[k+i] = S[z[n][i] - 'a'] ;
	    }
	}
      ppp[n] = 0 ; /* zero terminate the list */	
      pp = newMultiSum (ppp) ;
    }




  if (i > 6)
    {
      messcrash ("Trace of more than 6 matrices is not yet programmed") ;
    }
  pp->tt.Id2 = 0 ;
  return pp ;
} /* pauliTraceTT */

/*******************************************************************************************/

static POLYNOME pauliTrace (POLYNOME pp)
{
  POLYNOME p1, p2 ;
  if (!pp)
    return 0 ;
  p1 = pp->p1 ;
  p2 = pp->p2 ;
  if (p1) p1 = pp->p1 = pauliTrace (p1) ;
  if (p2) p2 = pp->p2 = pauliTrace (p2) ;

  if (pp->tt.type)
    {
      TT tt = pp->tt ;
      char *s = tt.sigma ; 
      char *sb = tt.sigB ; 
      if (s[0] && sb[0]) messcrash ("Cannot have sigma=%s and sigmaBar=%s in the same monome", s, sb) ;
      if (tt.Id2)
	{
	  pp->isFlat = FALSE ;
	  pp->tt = contractTtIndices (pp) ;
	  pp = pauliTraceTT (pp) ;
	}
      if (pp && pp->tt.type && cabs (pp->tt.z) < minAbs)
	{ pp = 0 ; }
    }
  return pp ;
} /* pauliTrace */
  
/*******************************************************************************************/
static KEYSET polynomeKs = 0 ;
static void checkPolynome (POLYNOME pp) 
{
  static int level = 0 ;
  static int nn ;
  int i ;

  if (! pp)
    return ;
  if (level == 0)
    {
      if (! polynomeKs)
	polynomeKs = keySetCreate () ;
      nn = 0 ;
    }
  level++ ;
  for (i = 0 ; i < nn ; i++)
    if (keySet (polynomeKs, i) == pp->id) 
      messcrash ("Duplicate node in polynome")  ;
  keySet (polynomeKs, nn++) = pp->id ; 
	
  if (pp->p1 == pp) messcrash ("pp == pp->p1 in checkPolynome") ;
  if (pp->p2 == pp) messcrash ("pp == pp->p2 in checkPolynome") ;
  if (pp->p1 && pp->p1 == pp->p2) messcrash ("pp->p1 == pp->p2 in checkPolynome") ;
  if (pp->p1 && pp->p2 && pp->p1->p2 == pp->p2) messcrash ("pp->p1->p2 == pp->p2 in checkPolynome") ;
  if (pp->p1) checkPolynome (pp->p1) ;
  if (pp->p2) checkPolynome (pp->p2) ;

  level-- ;
}

/*******************************************************************************************/
/* Flatten a polynome */
static POLYNOME expandDo (POLYNOME pp, int force)
{
  POLYNOME p1, p2 ;
  BOOL debug = FALSE ;

  if (!pp)
    return 0 ;
  if (force) pp->isFlat = FALSE ;
  p1 = pp->p1 ;
  p2 = pp->p2 ;
  if (p1 == pp) messcrash ("pp == pp->p1 in expandDo") ;
  if (p2 == pp) messcrash ("pp == pp->p2 in expandDo") ;
  if (p1 && p1 == p2) messcrash ("pp->p1 == pp->p2 in expandDo") ;

  if (pp->tt.type && cabs (pp->tt.z) < minAbs)
    return 0 ;
  if (p1 && p1->tt.type && cabs (p1->tt.z) < minAbs)
    p1 = pp->p1 = 0 ;
  if (p2 && p2->tt.type && cabs (p2->tt.z) < minAbs)
    p2 = pp->p2 = 0 ;
  
  if (1 && pp->tt.type && pp->tt.eps[7])
    {
      pp = contractEpsilon (pp) ;
      if (!pp)	return 0 ;
      if (! pp->tt.eps[7]) pp = expand (pp) ;
      if (!pp)	return 0 ;
      p1 = pp->p1 ;
      p2 = pp->p2 ;
    }
  if (pp->tt.type && cabs (pp->tt.z) < minAbs)
    return 0 ;
  if (p1 && p1->tt.type && cabs (p1->tt.z) < minAbs)
    p1 = pp->p1 = 0 ;
  if (p2 && p2->tt.type && cabs (p2->tt.z) < minAbs)
    p2 = pp->p2 = 0 ;
  
  if (! pp)
    return 0 ;

  if (pp->isProduct)
    {
      if (!p1 || !p2)
	return 0 ;
    }
  if (pp->isSum)
    {
      if (!p1 && !p2)
	return 0 ;
      if (!p2)
	return expandDo (pp->p1, force) ;
      if (!p1)
	return expandDo (pp->p2, force) ;
    }

  if (pp->p2 && (force || ! pp->p2->isFlat))
    p2 = pp->p2 = expandDo (pp->p2, force) ;
  if (pp->p1 && (force || ! pp->p1->isFlat))
    p1 = pp->p1 = expandDo (pp->p1, force) ;
  if (debug) checkPolynome (pp) ;

  if (pp->isProduct) /* check again */
    {
      if (!p1 || !p2)
	return 0 ;
    }
  if (pp->isSum)
    {
      if (!p1 && !p2)
	return 0 ;
      if (!p2)
	return expandDo (pp->p1, force) ;
      if (!p1)
	return expandDo (pp->p2, force) ;
    }

  if (pp->isProduct && p1 && p1->isSum)
    {
      POLYNOME q1 = newPolynome () ;
      POLYNOME q2 = newPolynome () ;
      pp->isSum = TRUE ;
      pp->isProduct = FALSE ;
      q1->isProduct = TRUE ;
      q2->isProduct = TRUE ;
      q1->p2 = p2 ;
      q2->p2 = copyPolynome (p2) ;
      q1->p1 = p1->p1 ;
      q2->p1 = p1->p2 ;
      pp->p1 = q1 ;
      pp->p2 = q2 ;
      if (debug) checkPolynome (pp) ;
      if (pp->p2 && (force || ! pp->p2->isFlat))
	p2 = pp->p2 = expandDo (q2, force) ;
      if (pp->p1 && (force || ! pp->p1->isFlat))
	p1 = pp->p1 = expandDo (q1, force) ;
      if (debug) checkPolynome (pp) ;
    }
  if (pp->isProduct && p2 && p2->isSum)
    {
      POLYNOME q1 = newPolynome () ;
      POLYNOME q2 = newPolynome () ;
      pp->isSum = TRUE ;
      pp->isProduct = FALSE ;
      q1->isProduct = TRUE ;
      q2->isProduct = TRUE ;
      q1->p1 = p1 ;
      q2->p1 = copyPolynome (p1) ;
      q1->p2 = p2->p1 ;
      q2->p2 = p2->p2 ;
      if (q2 && q2->p1 && q2->p1 == q2->p2) messcrash ("q2->p1 == q2->p2 in expandDo") ;
      p1 = pp->p1 = q1 ;
      p2 = pp->p2 = q2 ;
      if (debug) checkPolynome (pp) ;
      if (pp->p2 && (force || ! pp->p2->isFlat))
	p2 = pp->p2 = expandDo (q2, force) ;
      if (pp->p1 && (force || ! pp->p1->isFlat))
	p1 = pp->p1 = expandDo (q1, force) ;
      if (debug) checkPolynome (pp) ;
    }
  if (pp->isSum && p1 && p2 && p1->isSum && p1->p2 && p2->isSum)
    {
      if (debug) checkPolynome (pp) ;
      p2 = pp->p2 = expandDo (pp->p2, force) ;
      if (debug) checkPolynome (p2) ;
      p1 = pp->p1 = expandDo (pp->p1, force) ;
      if (debug) checkPolynome (p1) ;
    }

  if (pp->isSum && p1 && p2 && p1->isSum && p1->p2 && p2->isSum)
    {
      pp->p1 = p1->p1 ;
      p1->p1 = p1->p2 ;
      p1->p2 = pp->p2 ;
      pp->p2 = p1 ;
      p1 = pp->p1 ; p2 = pp->p2 ;
      if (debug) checkPolynome (pp) ;
    }

  if (pp->isSum && p1 && p1->tt.type && p2 && p2->isSum && p2->p1 && p2->p1->tt.type)
    {
      POLYNOME p3 = pp->p2->p1 ;
      float complex z1 = p1->tt.z ;
      float complex z3 = p3->tt.z ;
      int s ;
      p1->tt.z = 0 ;       
      p3->tt.z = 0 ;       
      s = memcmp (&(p1->tt), &(p3->tt), sizeof(TT)) ;
      p1->tt.z = z1 ;       
      p3->tt.z = z3 ;       

      if (s == 0)
	{ /* add p1 inside p2, skip p1, return p2 */
	  p3->tt.z = z1 + z3 ;
	  p2->isFlat = FALSE ;
	  return expandDo (pp->p2, force) ;
	}
      else if (s > 0)
	{ /* addition is Abelian */
	  pp->p2->p1 = p1 ;
	  p1 = pp->p1 = p3 ;
	  pp->isFlat = FALSE ;
	  pp->p2->isFlat = FALSE ;
	}
      p1 = pp->p1 ; p2 = pp->p2 ;
    }
  

  if (! pp)
    return 0 ;

  if (pp->tt.type && cabs (pp->tt.z) < minAbs)
    {
      return 0 ;
    }
  if (! pp->tt.type && !p1 && !p2)
    {
      return 0 ;
    }
  pp->isFlat = TRUE ;
  if (debug) checkPolynome (pp) ;
  return pp ;
} /* expandDo */

/*************************************************************************************************/

static POLYNOME expand (POLYNOME pp)
{
  int force = 1 ;

  pp = sortPol (pp) ;
  if (pp)
    {
      int nn = 12 ;
      pp->isFlat = FALSE ;
      while (pp && ! pp->isFlat && nn-- > 0)
	{
	  pp = expandDo (pp, force) ;
	  pp = contractProducts (pp) ;
	  pp = contractIndices (pp) ;
	  pp = sortPol (pp) ;
	  force = FALSE ;
	}
    }
  return pp ;
}

/*************************************************************************************************/

static POLYNOME killMomenta (POLYNOME pp)
{

  if (pp->p1) pp->p1 = killMomenta (pp->p1) ;
  if (pp->p2) pp->p2 = killMomenta (pp->p2) ;
  if (pp->isProduct)
    {
      if (! pp->p1 || ! pp->p2)
	return 0 ;
    }
  else if (pp->isSum)
    {
      if (! pp->p1)
	{ pp->p1 = pp->p2 ; pp->p2 = 0 ;}
      if (! pp->p2)
	pp = pp->p1 ;
    }
  else if (! pp->tt.type)
    return 0 ;
  else 
    {
      int i ;
      for (i = 1 ; i < 4 ; i++) /* do not kill k , just p,p,r */
	{
	  if (pp->tt.mm[i][0])
	    return 0 ;
	  pp->tt.denom[0] += pp->tt.denom[i] ; 
	  pp->tt.denom[i] = 0 ;
	}
      if (cabs (pp->tt.z) < minAbs)
	return 0 ;
    }
  if(pp) pp->isFlat = FALSE ;
  return pp ;
}

/*************************************************************************************************/
/* compute the dervative of (1/(k + p + q)^2*order) relative to (pqr)_mu */ 
static POLYNOME newDeriveDenom (POLYNOME p0, int pqr, int mu)
{
  TT tt = p0->tt ;
  POLYNOME ppp[6], pp = newPolynome () ;
  int j, nn = 0 ;
  BOOL debug = FALSE ;

  memset (ppp, 0, sizeof(ppp)) ;
  
  for (j = pqr ; j < 4 ; j++)
    {
      if (tt.denom[j]) /* ( k+p)^2 */
	{
	  int i, k, n ;
	  POLYNOME w, vv [5] ;
	  memset (vv, 0, sizeof(vv)) ;
	  if (j >=0) vv[0] = newK (mu) ;
	  if (j >=1) vv[1] = newP (mu) ;
	  if (j >=2) vv[2] = newQ (mu) ;
	  if (j >=3) vv[3] = newR (mu) ;
	  vv[4] = 0 ;
	  for (k = 0 ; k < 4 ; k++)
	    if (vv[k])
	      {
		for (i=0 ; i < 4 ; i++)
		  vv[k]->tt.denom[i] = tt.denom[i] ; /* passive factors in the denom */ 
		n = tt.denom[j] ;   /* the factor we derived */
		vv[k]->tt.denom[j] = n + 1 ;
		vv[k]->tt.z *= (-2 * n ) ;
	      }
	  w = newMultiSum (vv) ; /* 1/(k+p+q)^n  -> w = (k+p+q) = derivee du denominateur */
	  ppp[nn++] = w ;
	}
    }
  pp = 0 ;
  if (nn)
    {
      pp = newMultiSum (ppp) ;
      pp->isFlat = FALSE ;
      if (debug) checkPolynome (pp) ;  
    }
  return pp ;  
} /* newDeriveDenom */

/*************************************************************************************************/
/* partial derivative of a polynome with respect to p_mu or q_mu or r_mu */
static POLYNOME derivePdo (POLYNOME pp, int pqr, int mu)
{
  static     POLYNOME empty = 0 ;
  TT tt ;
  BOOL hasDenom = FALSE ;
  BOOL debug = FALSE ;

  if (! pp) return 0 ;
  if (debug) checkPolynome (pp) ;  
  if (pp->isSum)
    { /* linearity */
      pp->p1 = derivePdo (pp->p1, pqr, mu) ;
      pp->p2 = derivePdo (pp->p2, pqr, mu) ;
      pp->isFlat = FALSE ;
      return pp ;
    }
  if (pp->isProduct)
    {
      POLYNOME q1 = copyPolynome (pp->p1) ;
      POLYNOME q2 = copyPolynome (pp->p2) ;

      q1 = derivePdo (q1, pqr, mu) ;
      q2 = derivePdo (q2, pqr, mu) ;
      
      if (!q1 && ! q2)
	pp = 0 ; 
      else if (q1 && ! q2)
	pp->p1 = q1 ; 
      else if (! q1 && q2)
	pp->p2 = q2 ;
      else
	{ /* Leibnitz */
	  POLYNOME r1 = newProduct (q1, pp->p2) ;
	  POLYNOME r2 = newProduct (pp->p1, q2) ;
	  pp->isSum = TRUE ;
	  pp->isProduct = FALSE ;
	  pp->p1 = r1 ; 
	  pp->p2 = r2 ;
	}
      if (pp) pp->isFlat = FALSE ;
      if (debug) checkPolynome (pp) ;  
      return pp ;
    }

  /* free object */
  if (pqr <1 || pqr > 3)
      messcrash ("You can only partial derive with respect to 1:p, 2:q, 3:r, not %d\n", pqr) ;
    
  /* derive the numerator */
  if (! empty)
    empty = newPolynome () ;

  
  tt = pp->tt ;

  if (! hasDenom)
    {
      int i ;
      for (i = pqr ; i < 4 ; i++)
	if (tt.denom[i])
	  hasDenom = TRUE ;
    }
  if (! tt.mm[pqr][0] && ! hasDenom)
    return 0 ;
  else if (tt.mm[pqr][0] && ! hasDenom)
    {
      char *u0 = tt.mm[pqr] ;
      int i, k, iMax = strlen (u0)  ;
      
      if (iMax == 1) /* simplest case, just add a g_munu and suppress the p_mu */
	{
	  pp = copyPolynome (pp) ;
	  char *v = tt.g, *w = tt.mm[pqr] ;
	  v += strlen (v) ;
	  v[0] = mu ; v[1] = tt.mm[pqr][0] ; v[2] = 0 ;
	  memset (w, 0, GMAX) ;

	  pp->tt = tt ;	  pp->isFlat = FALSE ;
	  if (debug) checkPolynome (pp) ;  
	  return pp ;
	}
      else 
	{
	  POLYNOME qq[iMax+1] ;
	  for (k = 0 ; k < iMax ; k++)
	    {
	      char *v,*w = tt.mm[pqr] ;
	      /* copy the original polynome */
	      qq[k] = newPolynome () ;
	      qq[k]->tt.type = 1 ;
	      qq[k]->tt = tt ;
	      /* replace one dependence on p_alpha by g_mu_alpha */
	      v = qq[k]->tt.mm[pqr] ; 
	      for (i = k ; i < iMax ; i++)
		v[i] = v[i+ 1] ; 
	      v = qq[k]->tt.g ;
	      v += strlen (v) ;
	      v[0] = mu ; v[1] = u0[k] ; v[2] = 0 ;
	      memset (w, 0, GMAX) ;
	    }
	  qq[iMax] = 0 ;
	  pp = newMultiSum (qq) ;
	  if (debug) checkPolynome (pp) ;  
	}
    }
  else if (hasDenom)
    { /* construct the product (f'/g + fg'/g^2) */
      POLYNOME q3, q1prime, q1 = newScalar (1) ;
      POLYNOME q4, q2prime, q2 = newScalar(1) ;
      int i ;

      /* contruct f: copy and remove the denom */
      q1->tt = tt ;
      for (i = 0 ; i < 4 ; i++)
	q1->tt.denom[i] = 0 ;
      /* contruct g: copy the denom, but certainly not the matrices */
      for (i = 0 ; i < 4 ; i++)
	q2->tt.denom[i] = tt.denom[i] ;

      /* derive */
      q1prime = derivePdo (q1, pqr, mu) ;
      q2prime = newDeriveDenom(q2, pqr, mu) ;

      /* construct the 2 products (q1 and q1prime both commute) */
      q3 = newProduct (q1, q2prime) ;
      q4 = q1prime ? newProduct (q2, q1prime) : 0 ;
      pp = q4 ? newSum (q3,q4) : q3 ;
      if (debug) checkPolynome (pp) ;  
      /*
      ac_free (q1) ;      
      ac_free (q2) ;
      ac_free (q3) ;
      ac_free (q4) ;
      ac_free (q5) ;
      ac_free (q1prime) ;
      ac_free (q2prime) ;
      */
    }  
  if (debug) checkPolynome (pp) ;  
  if (pp) pp->isFlat = FALSE ;
  return pp ;
} /* derivePdo */

/*************************************************************************************************/

static POLYNOME deriveP (POLYNOME pp, int pqr, int mu)
{
  BOOL debug = FALSE ;

  pp = derivePdo (pp, pqr, mu) ;
  if (0 && pp) pp = killMomenta (pp) ;
  pp = expand (pp) ;
  if (debug) showPol (pp) ;
  pp = contractProducts (pp) ;
  if (debug) showPol (pp) ;
  pp = contractIndices (pp) ;
  if (debug) showPol (pp) ;
  if (pp) pp->isFlat = FALSE ;
  return pp ;
} /* deriveP */

/*************************************************************************************************/
/*************************************************************************************************/
static POLYNOME dimIntegral (POLYNOME pp) ;
static POLYNOME dimIntegrateByPart (POLYNOME pp) ;
/* dimensional integral, only works on a flat expression of total degree 4 */
static POLYNOME dimIntegralMonome (POLYNOME pp, int state, char *kk, int *np, int pass)
{
  int i, j, k ;
  BOOL debug = FALSE ;

  if (! pp) return 0 ;
  if (debug) checkPolynome (pp) ;  
  if (debug)
    showPol(pp) ;
  switch (state)
    {
    case 0:
      if (1)
	{
	  int nd = 0 ;
	  char buf[GMAX] ;
	  
	  memset (buf, 0, sizeof(buf)) ;
	  dimIntegralMonome (pp, 1, buf, &nd, pass) ;
	  if (debug) checkPolynome (pp) ;  
	  k = 4 + strlen (buf) - 2 * nd ;
	  if (k < 0) /* convergent integral */
	    return 0 ;
	  if (! pass)
	    return pp ;
	  else if (k == 1)
	    return dimIntegrateByPart (pp) ;
	  else if (k == 2)
	    {
	      pp = dimIntegrateByPart (pp) ;
	      if (! pp)
		return 0 ;
	      if (pp->isProduct)
		pp->p1->tt.z /= 2 ;
	      else if (pp->isSum)
		{
		  POLYNOME qq = newScalar (.5) ;
		  pp = newProduct (qq, pp) ;
		}
	      else
		pp->tt.z /= 2 ;
	      return pp ;
	    }
	  else if (k)
	    messcrash ("This integral should have been differentiated, it has k-order %d > 0", k) ;

	  dimIntegralMonome (pp, 2, 0, 0, pass) ; /* clean up all denoms and all k dependencies */
	  if (debug) checkPolynome (pp) ;  
	  /* transform the polynome into a sum of g_munu */
	  k = strlen (buf) ;
	  if (k) 
	    {
	      /* eliminate repeated indices    k^2 / k^4 (k+p)^2 = 1/k^2 (k+p)^2 */
	      for (i = 0 ; i < k - 1 ; i++)
		for (j = i + 1 ; j < k ; j++)
		  if (buf[i] == buf[j])
		    buf[i] = buf[j] = 0 ;
	      for (i = j = 0 ; i < k ; i++)
		if (buf[i]) { buf[j] = buf[i] ; j++ ; }
	      buf[j] = 0 ;
	    }
	  k = strlen (buf) ; /* simplified length */
	  if (k == 0)
	    {
	    }
	  else if (k == 2)
	    { 
	      int n, N = 2, NN = 1 ;
	      POLYNOME qqq[NN+1] ;
	      char *cp ;
	      pp->tt.z /= 4 ;
	      for (n = 0 ; n < NN ; n++)
		{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
		  int i ;
		  qqq[n] = copyPolynome (pp) ;
		  cp = qqq[n]->tt.g ;
		  while (*cp) cp++ ;
		  
		  for (i = 0 ; i < N ; i++)
		    *cp++ = buf[i] ;
		  *cp++ = 0 ;
		}
	      qqq[n] = 0 ; /* zero terminate the list */	
	      pp = newMultiSum (qqq) ;
	    }	  
	  else if (k == 4)
	    { 
	      int n, N = 4, NN = 3 ;
	      POLYNOME qqq[NN+1] ;
	      char *cp ;
	      char *z[3] = {
		"abcd","acbd","adbc"
	      } ;
	      pp->tt.z /= 24 ;
	      for (n = 0 ; n < NN ; n++)
		{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
		  int i ;
		  qqq[n] = copyPolynome (pp) ;
		  cp = qqq[n]->tt.g ;
		  while (*cp) cp++ ;

		  for (i = 0 ; i < N ; i++)
		    *cp++ = buf[z[n][i] - 'a'] ;
		  *cp++ = 0 ;
		}
	      qqq[n] = 0 ; /* zero terminate the list */	
	      pp = newMultiSum (qqq) ;
	    }	  
	  else if (k == 6)
	    { 
	      int n, N = 6, NN = 15 ;
	      POLYNOME qqq[NN+1] ;
	      char *cp ;
	      char *z[15] = {
		"abcdef","abcedf","abcfde",
		"acbdef","acbedf","acbfde",
		"adbcef","adbecf","adbfce",
		"aebcdf","aebdcf","aebfcd",
		"afbcde","afbdce","afbecd"
	      } ;
	      pp->tt.z /= 192 ;
	      for (n = 0 ; n < NN ; n++)
		{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
		  int i ;
		  qqq[n] = copyPolynome (pp) ;
		  cp = qqq[n]->tt.g ;
		  while (*cp) cp++ ;

		  for (i = 0 ; i < N ; i++)
		    *cp++ = buf[z[n][i] - 'a'] ;
		  *cp++ = 0 ;
		}
	      qqq[n] = 0 ; /* zero terminate the list */	
	      pp = newMultiSum (qqq) ;
	    }	  
	  else if (k == 8)
	    { 
	      int n, N = 8, NN = 105 ;
	      POLYNOME qqq[NN+1] ;
	      char *cp ;
	      char *z[105] = {
		"abcdefgh","abcdegfh","abcdehfg",
		"abcedfgh","abcedgfh","abcedhfg",
		"abcfdegh","abcfdgeh","abcfdheg",
		"abcgdefh","abcgdfeh","abcgdhef",
		"abchdefg","abchdfeg","abchdgef",
		
		"acbdefgh","acbdegfh","acbdehfg",
		"acbedfgh","acbedgfh","acbedhfg",
		"acbfdegh","acbfdgeh","acbfdheg",
		"acbgdefh","acbgdfeh","acbgdhef",
		"acbhdefg","acbhdfeg","acbhdgef",
		
		"adcbefgh","adcbegfh","adcbehfg",
		"adcebfgh","adcebgfh","adcebhfg",
		"adcfbegh","adcfbgeh","adcfbheg",
		"adcgbefh","adcgbfeh","adcgbhef",
		"adchbefg","adchbfeg","adchbgef",
		
		"aecdbfgh","aecdbgfh","aecdbhfg",
		"aecbdfgh","aecbdgfh","aecbdhfg",
		"aecfdbgh","aecfdgbh","aecfdhbg",
		"aecgdbfh","aecgdfbh","aecgdhbf",
		"aechdbfg","aechdfbg","aechdgbf",
		
		"afcdebgh","afcdegbh","afcdehbg",
		"afcedbgh","afcedgbh","afcedhbg",
		"afcbdegh","afcbdgeh","afcbdheg",
		"afcgdebh","afcgdbeh","afcgdheb",
		"afchdebg","afchdbeg","afchdgeb",
		
		"agcdefbh","agcdebfh","agcdehfb",
		"agcedfbh","agcedbfh","agcedhfb",
		"agcfdebh","agcfdbeh","agcfdheb",
		"agcbdefh","agcbdfeh","agcbdhef",
		"agchdefb","agchdfeb","agchdbef",
		
		"ahcdefgb","ahcdegfb","ahcdebfg",
		"ahcedfgb","ahcedgfb","ahcedbfg",
		"ahcfdegb","ahcfdgeb","ahcfdbeg",
		"ahcgdefb","ahcgdfeb","ahcgdbef",
		"ahcbdefg","ahcbdfeg","ahcbdgef"
	      } ;
	      if (4) pp->tt.z /= 1920 ;
	      for (n = 0 ; n < NN ; n++)
		{                             /* we need N products of type g_ab g_cd g_ef, then we zero terminate the list */
		  int i ;
		  qqq[n] = copyPolynome (pp) ;
		  cp = qqq[n]->tt.g ;
		  while (*cp) cp++ ;

		  for (i = 0 ; i < N ; i++)
		    *cp++ = buf[z[n][i] - 'a'] ;
		  *cp++ = 0 ;
		}
	      qqq[n] = 0 ; /* zero terminate the list */	
	      pp = newMultiSum (qqq) ;
	    }	  
	  else 
	    messcrash ("Sorry, i cannot yet integrate k^10/k^14\n") ;
	}
      if (debug) checkPolynome (pp) ;  
      break ;

    case 1:
      if (pp->isSum)
	messcrash ("dimIntegralMonome found a sum hidden under a product, the polynome is not flat") ;
      if (pp->isProduct)
	{
	  dimIntegralMonome (pp->p1, 1, kk, np, pass) ;
	  dimIntegralMonome (pp->p2, 1, kk, np, pass) ;
	}
      if (pp->tt.type)
	{
	  int i ;
	  char *u, *v ;
	  
	  for (i = 0 ; i < 4 ; i++)
	    *np += pp->tt.denom[i] ;  /* power of k in 1/(k+...)^2 */
	  u = kk ; u += strlen(u) ; 
	  v = pp->tt.mm[0] ;        /* numerator k_mu k_nu ... */
	  while ((*u++ = *v++)) ;
	  if (kk[GMAX-1])
	    messcrash ("Overflow while collecting the k[] vectors ") ;
	}
      break ;
    case 2:  /* clean up the k dependency */
      if (pp->isProduct)
	{
	  dimIntegralMonome (pp->p1, 2, 0, 0, pass) ;
	  dimIntegralMonome (pp->p2, 2, 0, 0, pass) ;
	}
      if (pp->tt.type)
	{
	  for (i = 0 ; i < 4 ; i++)
	    pp->tt.denom[i] = 0 ;   /* power of k in 1/(k+...)^2 */
	  memset (pp->tt.mm[0], 0, GMAX) ; /* numerator k_mu k_nu ... */
	}
      break ;
    }
  if (pp) pp->isFlat = FALSE ;
  if (debug) checkPolynome (pp) ;
  if (debug)
    showPol(pp) ;
  return pp ;
}  /* dimIntegralMonome */
  
/*************************************************************************************************/

static POLYNOME dimIntegrateByPart (POLYNOME pp) 
{
  POLYNOME ppp[4] ;
  int nn = 0 ;
  int i, j ;
  char mu = newDummyIndex () ;
  static int level = 0 ;
  BOOL debug = FALSE ;

  level++ ;
  for (j = 1 ; j <= 3 ; j++)
    {
      POLYNOME p1d = 0, p1 = copyPolynome (pp) ;
      POLYNOME p2d = 0, p2 = newScalar (1) ;
      POLYNOME puu = 0, pvv = 0 ;
 
     if (debug)
       {
	 int i ;
	 for (i = 0 ; i < level ; i++)
	   printf ("###") ;
	 printf ("# level %d before\n", level) ;
	 showPol (pp) ;
       }
      for (i = 0 ; i < 4 ; i++)
	{
	  p2->tt.denom[i] = p1->tt.denom[i] ;
	  p1->tt.denom[i] = 0 ; 
	}
      for (i = 0 ; i < GMAX ; i++)
	{
	  p2->tt.mm[0][i] = p1->tt.mm[0][i] ;
	  p1->tt.mm[0][i] = 0 ;
	}

      if (0) p1d = deriveP (p1, j, mu) ;
      p2d = deriveP (p2, j, mu) ;
      if (p1d)
	{
	  POLYNOME pm = newScalar (1.0) ;
	  pm->tt.mm[j][0] = mu ;
	  puu = dimIntegralDo (p1d, 1) ;
	  if (puu) puu = newProduct (pm, puu) ;
	}
      if (p2d)
	{
	  POLYNOME pm = newScalar (1.0) ;
	  pm->tt.mm[j][0] = mu ;
	  pvv = dimIntegralDo (p2d, 1) ;
	  if (pvv) pvv = newProduct (pm, pvv) ;
	}
      if (puu) puu = newProduct (p2, puu) ;
      if (pvv) pvv = newProduct (p1, pvv) ;
      if (puu && pvv) 
	ppp[nn++] = newSum (puu, pvv) ;
      else if (puu)
	ppp[nn++] = puu ;
      else if (pvv)
	ppp[nn++] = pvv ;
     if (debug)
       {
	 int i ;
	 for (i = 0 ; i < level ; i++)
	   printf ("###") ;
	 printf ("# level %d after\n", level) ;
	 showPol (pvv) ;
       }
    }
  ppp[nn] = 0 ;

  if (nn > 1)
    pp = newMultiSum (ppp) ;
  else if (nn == 1)
    pp = ppp[0] ;
  else
    pp = 0 ;

  level-- ;
  return pp ;
}

/*************************************************************************************************/
/* dimensional integral, only works on a flat expression of total degree 4 */
static POLYNOME dimIntegralDo (POLYNOME pp, int pass)
{  
  BOOL debug = FALSE ;
  static int level = 0 ;
  static int nnn = 0 ;
  
  if (! pp) return 0 ;
  level++ ;
  if (level == -1) firstDummyIndex += 4 ;
  pp = copyPolynome (pp) ;
  pp = expand (pp) ; 
  pp = contractProducts (pp) ;
  contractIndices (pp) ;
  if (! pp) { level-- ; return 0 ; }
 
  if (debug && level <=2) 
    {
      int i = level ;
      while (i--)
	printf ("##") ;
      showPol(pp) ;
    }
  if (pp->isSum)
    { /* linearity */
      int ii = 0, jj = 0 ;
      Array aa = arrayCreate (32, POLYNOME) ;
      
      sortPolGetSum (pp, aa) ;
      for (ii = jj = 0 ; ii < arrayMax (aa) ; ii++)
	{
      	  POLYNOME qq = arr (aa, ii, POLYNOME) ; 
	  nnn++ ;
	  freeIndex (qq) ;
	  qq = dimIntegralDo (qq, pass) ;
	  freeIndex (qq) ;
	  if (qq)
	    arr (aa, jj++, POLYNOME) = qq ;
	}
      arrayMax (aa) = jj ;
      pp = sortReduceSum (aa) ;
   }
  else
    {
      pp = dimIntegralMonome (pp, 0, 0, 0, pass) ;
    }
  if (pp) 
    {
      pp->isFlat = FALSE ;
      pp = expand (pp) ;
    }
  if (debug && level <=2) 
    {
      int i = level ;
      while (i--)
	printf ("XX") ;
      showPol(pp) ;
    }
  level-- ;
  return pp ;
}   /* dimIntegralDo */

/***********************************************************************************************************************************************/
/* given a list of tensor indices, fuse and eliminate repated index */
static BOOL freeIndexFuse (char *top, char **cpp)
{
  int i, ii ;
  char *cp, nn[256] ;
  memset (nn, 0, sizeof (nn)) ;

  for (ii = 0, cp = cpp[ii] ; cp ; cp = cpp[++ii])
    while (*cp)
      nn[(int)*cp++]++ ;
  for (ii = i = 0 ; i < 256 ; i++)
    if (nn[i] == 1)
      top[ii++] = i ;
    else if (nn[i] > 2)
      messcrash ("Index %c is repeated %d times in this tensor", i, nn[i]) ;
  top[ii] = 0 ;
  if (ii >= GMAX)
    messcrash ("too many indicies in this tensor %s", top) ;
  return TRUE ;

} /* freeIndexFuse */

/***********************************************************************************************************************************************/

static BOOL freeIndex (POLYNOME pp)
{
  POLYNOME p1, p2 ;
  
  if ( !pp)
    return TRUE ; /* no problem */ ;
  p1 = pp->p1 ;
  p2 = pp->p2 ;

  memset (pp->tt.freeIndex, 0, GMAX) ;
  if (pp->isProduct)
    {
      freeIndex (p1) ;
      freeIndex (p2) ;
      if (! p1 && ! p2)
	return TRUE ;
      if (p1 && ! p2)
	{
	  strcpy (pp->tt.freeIndex, p1->tt.freeIndex) ;
	  return TRUE ;
	}
      if (! p1 && p2)
	{
	  strcpy (pp->tt.freeIndex, p2->tt.freeIndex) ;
	  return TRUE ;
	}
      if (p1 && p2)
	{
	  char *cp = pp->tt.freeIndex ;
	  char *cp1 = p2->tt.freeIndex ;
	  char *cp2 = p1->tt.freeIndex ;
	  char *cpp[] = {cp1, cp2, 0 } ;

	  freeIndexFuse (cp, cpp) ;
	}
      return TRUE ; /* a product may have any set of free indicies */
    }
  
  else if (pp->isSum)
    {
      freeIndex (p1) ;
      freeIndex (p2) ;
      if (! p1 && ! p2)
	return TRUE ;
      if (p1 && ! p2)
	{
	  strcpy (pp->tt.freeIndex, p1->tt.freeIndex) ;
	  return TRUE ;
	}
      if (! p1 && p2)
	{
	  strcpy (pp->tt.freeIndex, p2->tt.freeIndex) ;
	  return TRUE ;
	}
      if (p1 && p2)
	{
	  char *cp = pp->tt.freeIndex ;
	  char *cp1 = p1->tt.freeIndex ;
	  char *cp2 = p2->tt.freeIndex ;

	  if (strcmp (cp1, cp2))
	    messcrash ("Index not repeated in a sum  #%s#  #%s#", p1->tt.freeIndex, p2->tt.freeIndex) ;
	  strcpy (cp, cp1) ;
	}
      return TRUE ; /* a product may have any set of free indicies */
    }
  else if (pp->tt.type)
    {
      int jj, n = 0 ;
      char *cp = pp->tt.freeIndex ;
      char *cpp[12] ;

      cpp[n++] = pp->tt.g ;
      cpp[n++] = pp->tt.sigma ;
      cpp[n++] = pp->tt.sigB ;
      cpp[n++] = pp->tt.eps ;
      for (jj = 0 ; jj < 4 ; jj++)
	cpp[n++] = pp->tt.mm[jj] ;
      cpp[n++] = 0 ;
      freeIndexFuse (cp, cpp) ;
    }
  
  return TRUE ;
}

/***********************************************************************************************************************************************/

static POLYNOME dimIntegral (POLYNOME p0)
{
  POLYNOME pp ;
  
  pp = copyPolynome (p0) ;
  pp = expand (pp) ;
  freeIndex (pp) ;
  pp = contractProducts (pp) ;
  freeIndex (pp) ;
  contractIndices (pp) ;
  freeIndex (pp) ;
  pp = dimIntegralDo (pp, 0) ;
  freeIndex (pp) ;
  pp = dimIntegralDo (pp, 1) ;
  freeIndex (pp) ;

  return pp ;
} /* dimIntegral */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/* New vertex A BB H and A H BB */
static POLYNOME vertex_A_B_HB (char mu, char a, char b, int mm[4]) /* A_mu B_nu_rho, momentum of the incoming photon */
{
  POLYNOME pp, ppp[6] ; 
  int nn = 0 ; 
  
  if (mm[0]) { pp = newK (b) ; pp->tt.z = mm[0] ; ppp[nn++] = pp ; }
  if (mm[1]) { pp = newP (b) ; pp->tt.z = mm[1] ; ppp[nn++] = pp ; }
  if (mm[2]) { pp = newQ (b) ; pp->tt.z = mm[2] ; ppp[nn++] = pp ; }
  if (mm[3]) { pp = newR (b) ; pp->tt.z = mm[3] ; ppp[nn++] = pp ; }
  ppp[nn++] = 0 ;

  pp = newMultiSum (ppp) ;
  nn = 0 ;
  if (0) ppp[nn++] = newScalar (4) ;
  ppp[nn++] = pp ;
  ppp[nn++] = newG (mu, a) ;
  ppp[nn++] = 0 ;
  
  pp = newMultiProduct (ppp) ;
  return pp ;
}

/***********************************************************************************************************************************************/
/* New vertex A BB H and A H BB */
static POLYNOME vertex_A_H_BB_do (char mu, char a, char b, int mm[4], int z) /* momentum of the incoming mu-photon */
{
 POLYNOME pp, ppp[6] ; 
  int nn = 0 ; 
  
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;

  if (mm[0]) { pp = newK (d) ; pp->tt.z = z * mm[0] ; ppp[nn++] = pp ; }
  if (mm[1]) { pp = newP (d) ; pp->tt.z = z * mm[1] ; ppp[nn++] = pp ; }
  if (mm[2]) { pp = newQ (d) ; pp->tt.z = z * mm[2] ; ppp[nn++] = pp ; }
  if (mm[3]) { pp = newR (d) ; pp->tt.z = z * mm[3] ; ppp[nn++] = pp ; }
  ppp[nn++] = 0 ;

  pp = newMultiSum (ppp) ;
  return pp ;
  
  nn = 0 ;
  if (0) ppp[nn++] = newScalar (4) ;
  ppp[nn++] = pp ;
  ppp[nn++] = newG (mu, c) ;
  ppp[nn++] = newEpsilon (a,b,c,d) ;
  ppp[nn++] = 0 ;
  
  pp = newMultiProduct (ppp) ;
  

  return pp ;
}

static POLYNOME vertex_A_H_BB (char mu, char a, char b, int mm[4]) /* momentum of the incoming mu-photon */
{
  POLYNOME p1 = vertex_A_H_BB_do (mu, a, b, mm, 1) ;
  return p1 ;
  POLYNOME p2 = vertex_A_H_BB_do (mu, b, a, mm, -1) ;
  return newSum (p1, p2) ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/*   A_mu diffusion of B_ab -> BB_cd */
static POLYNOME vertex_A_B_BB (char mu, char a, char b, char c, char d, int  mm1[4], int  mm2[4])   /* B(mm1)->BB(mm2) in the orientation of the BB->B line */
{
  POLYNOME pp, ppp[6] ; 
  int i, nn = 0 ; 

  for (i = 0 ; i < 4 ; i++)
    {
      if (mm1[i])
	{
	  pp = newScalar (4 * mm1[i]) ;
	  pp->tt.mm[i][0] = a ;
	  pp->tt.g[0] = mu ;
	  pp->tt.g[1] = c ;
	  pp->tt.g[2] = b ;
	  pp->tt.g[3] = d ;
	  ppp[nn++] = pp ;	  
	}
      if (mm2[i])
	{
	  pp = newScalar (4 * mm2[i]) ;
	  pp->tt.mm[i][0] = c ;
	  pp->tt.g[0] = mu ;
	  pp->tt.g[1] = a ;
	  pp->tt.g[2] = b ;
	  pp->tt.g[3] = d ;
	  ppp[nn++] = pp ;	  
	}
    }
  ppp[nn++] = 0 ;	  
  pp = newMultiSum (ppp) ;
  return pp ;
} /* vertex_A_B_BB */

/***********************************************************************************************************************************************/
/*   A_mu diffusion of B_ab -> BB_cd */
static POLYNOME vertex_A_A_A (char mu, char nu, char rho, char  mmu[4], char  mnu[4], char  mrho[4])   /* all 3 are incoming momenta */
{
  POLYNOME p = newK (mu) ;
  int i ;
  for (i = 3 ; i < 4 ; i++)
    p->tt.mm[i][0] = mmu[i] ;
  return p ;
}

/**********************************************************************************************************************************************/
/*   A_mu diffusion of H -> HB */
static POLYNOME vertex_A_H_HB (char mu, int mm[4])  /* 2k+p = (2,1,0,0) : sum of the momentum of the H->HB in that direction */
{
  POLYNOME pp, ppp[6] ; 
  int i, nn = 0 ; 

  for (i = 0 ; i < 4 ; i++)
    if (mm[i])
      {
	pp = newScalar (mm[i]) ;
	pp->tt.mm[i][0] = mu ;
	ppp[nn++] = pp ;	  
      }
  ppp[nn++] = 0 ;	  
  pp = newMultiSum (ppp) ;
  return pp ;
} /* vertex_A_H_HB */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/* This vertex is by itelf self dual */
static POLYNOME vertex_B_PsiR_PsiLB (char mu, char nu)
{
  POLYNOME p1 = newSigma (mu) ;
  POLYNOME p2 = newSigma (nu) ;
  p1->tt.sigma[1] = nu ;
  p2->tt.sigma[1] = mu ;
  p1->tt.z = 1.0/2 ;
  p2->tt.z = -1.0/2 ;

  return newSum (p1, p2) ; ;
}

/***********************************************************************************************************************************************/
/* This vertex is by itelf anti self dual */
static POLYNOME vertex_BB_PsiL_PsiRB (char mu, char nu)
{
  POLYNOME p1 = newSigB (mu) ;
  POLYNOME p2 = newSigB(nu) ;
  p1->tt.sigB[1] = nu ;
  p2->tt.sigB[1] = mu ;
  p1->tt.z = 1.0/2 ;
  p2->tt.z = -1.0/2 ;

  return newSum (p1, p2) ; ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME vertex_A_PsiR_PsiRB (char mu)
{
  POLYNOME p = newSigma (mu) ;
  return p ;
}

/***********************************************************************************************************************************************/

static POLYNOME vertex_A_PsiL_PsiLB (char mu)
{
  POLYNOME p = newSigB (mu) ;
  return p ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME vertex_H_PsiR_PsiLB (void)
{
  POLYNOME p = newScalar (1) ;
  return p ;
}

/***********************************************************************************************************************************************/

static POLYNOME vertex_HB_PsiL_PsiRB (void)
{
  POLYNOME p = newScalar (1) ;
  return p ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME squareMomentaCleanUp (POLYNOME pp, char alpha) 
{
  static int level = 0 ;

  if (level == 0 && ! alpha)
    {
      if (0) pp = reduceIndices (pp) ;
      alpha = newDummyIndex() ;
    } 
  level++ ;

  if (! pp)
    return 0 ;

  if (pp->isSum && pp->p1)
    pp->p1 = squareMomentaCleanUp (pp->p1, alpha) ;
  if (pp->isSum && pp->p2)
    pp->p2 = squareMomentaCleanUp (pp->p2, alpha) ;
  if (pp->tt.type == 1)
    {
      if (pp->tt.mm[1][0] && pp->tt.mm[1][0] == pp->tt.mm[1][1])
	pp->tt.mm[1][0] = pp->tt.mm[1][1] = alpha ;
    }
  pp = expand (pp) ;
  level-- ;

  return pp ;
}
/***********************************************************************************************************************************************/

static POLYNOME momentaCleanUp (POLYNOME pp, char alpha) 
{
  int nn = 0 ;
  POLYNOME ppp[4] ;
  POLYNOME p1, p2, p3, p4, p5 ;
  p1 = copyPolynome (pp) ;
  p2 = deriveP (p1, 1, alpha) ;
  p3 = expand (p2) ;
  p4 = contractIndices (p3) ;
  p5 = expand (p4) ;
  if (p5)
    ppp[nn++] = newProduct (newP(alpha), p5) ;

  p1 = copyPolynome (pp) ;
  p2 = deriveP (p1, 2, alpha) ;
  p3 = expand (p2) ;
  p4 = contractIndices (p3) ;
  p5 = expand (p4) ;
  if (p5)
    ppp[nn++] = newProduct (newQ(alpha), p5) ;
  p1 = copyPolynome (pp) ;
  p2 = deriveP (p1, 3, alpha) ;
  p3 = expand (p2) ;
  p4 = contractIndices (p3) ;
  p5 = expand (p4) ;
  if (p5)
    ppp[nn++] = newProduct (newR(alpha), p5) ;
  ppp[nn++] = 0 ;

  return newMultiSum (ppp) ;
}
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME prop_BB_B (char mu, char nu, char rho, char sig, int pqr)
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  int z = 1 ;

  POLYNOME p1 = newAG (mu,nu,a,b, z) ;
  if (0) p1 = newEpsilon (mu, nu, a,b) ;
  POLYNOME p2 = newPQR (pqr, a) ;
  POLYNOME p3 = newPQR (pqr, c) ;
  POLYNOME p4 = newG  (b, d) ;
  POLYNOME p5 = newAG (c,d,rho,sig, -z) ;
  if (0) p5 = newEpsilon (c,d,rho,sig) ;
  POLYNOME pp, ppp[] = {p1,p2,p3, p4, p5, 0} ;

  p4->tt.denom[pqr] = 2 ;
  pp = newMultiProduct (ppp) ;

  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME prop_AL (char mu, char nu, int pqr)
{
  POLYNOME pp = newG (mu,nu) ;
  pp->tt.denom[pqr] = 1 ;
  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME prop_AR (char mu, char nu, int pqr)
{
  POLYNOME pp = newG (mu,nu) ;
  pp->tt.denom[pqr] = 1 ;
  pp->tt.z = -1 ; /* negative norm */
  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME prop_HB_H (int pqr) 
{
  POLYNOME pp = newScalar (1.0) ;
  pp->tt.denom[pqr] = 1 ;
  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME prop_PsiRB_PsiR (int pqr)
{
  char cc = newDummyIndex () ;
  POLYNOME p1 = newSigB (cc) ;
  POLYNOME p2 = newPQR (pqr,cc) ;
  POLYNOME pp = contractIndices(newProduct (p1, p2)) ;
  p1->tt.denom[pqr] = 1 ;
  return contractProducts (pp) ;
}

/***********************************************************************************************************************************************/

static POLYNOME prop_PsiLB_PsiL (int pqr)
{
  char cc = newDummyIndex () ;
  POLYNOME p1 = newSigma (cc) ;
  POLYNOME p2 = newPQR (pqr,cc) ;
  POLYNOME pp = contractIndices(newProduct (p1, p2)) ;
  p1->tt.denom[pqr] = 1 ;
  return contractProducts (pp) ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME pauliCleanUp (POLYNOME pp)
{
  char w = newDummyIndex () ;
  POLYNOME p1 = pp->tt.sigma[0] || (pp->p1 && pp->p1->tt.sigma[0]) ? newSigB (w) : newSigma (w) ;
  POLYNOME p2 = newProduct (p1, pp) ;
  POLYNOME p3 = expand (p2) ;
  POLYNOME p4 = pauliTrace (p3) ;
  POLYNOME p5 = expand (p4) ;
  POLYNOME p6 = contractIndices (p5) ;
  POLYNOME p7 = newProduct (p1, p6) ;
  POLYNOME p8 = expand (p7) ;
  POLYNOME p9 = contractIndices (p8) ;
  POLYNOME p10 = expand (p9) ;
  if (p10) p10->tt.z /= 2 ;
 return p10 ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME Z2_AA_loopH  (const char *title) 
{
  firstDummyIndex = 'a' ;
    
  char mu = newDummyIndex () ;
  char nu = newDummyIndex () ;
  int ppv[4] = {2,1,0,0} ; /* 2k + p : vertex */
  
  POLYNOME p1 = vertex_A_H_HB (mu, ppv) ; /*(2k + p)_mu */
  POLYNOME p2 = prop_HB_H (1) ;   /* (1/(k)^2 */
  POLYNOME p3 = vertex_A_H_HB (nu, ppv) ; /* (2k + p)_mu */prop_PsiRB_PsiR (1) ; /* (1/(k+p)^2 */
  POLYNOME p4 = prop_HB_H (0) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(newMultiProduct (ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  printf ("### Z2 Psi left avec loop B_mu_nu, expect ::  je_sais_pas * p-slash\n") ;
  showPol (pp) ;
  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_AA_loopH */

/***********************************************************************************************************************************************/

static POLYNOME Z2_AA_loopB (const char *title) 
{
  firstDummyIndex = 'a' ;

  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;
  char g = newDummyIndex () ;
  char h = newDummyIndex () ;
  char i = newDummyIndex () ;
  char j = newDummyIndex () ;
  char n = newDummyIndex () ;

  int mm1[4] = {1,1,0,0} ;    /* k+p */
  int mm2[4] = {1,0,0,0} ;    /* k */

  POLYNOME p0 = newScalar (96) ;
  POLYNOME p1 = vertex_A_B_BB (a,c,d,i,j,mm1,mm2) ;
  POLYNOME p2 = prop_BB_B (c,d,e,f,1) ;   /* (1/(k+p)^2 */
  POLYNOME p3 = vertex_A_B_BB (b,g,h,e,f,mm2,mm1) ;
  POLYNOME p4 = prop_BB_B (g,h,i,j,0) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p0, p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(newMultiProduct (ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;
  freeIndex (pp) ;
  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;
  freeIndex (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  pp = squareMomentaCleanUp (pp, n) ;
  showPol(pp) ;
  pp = expand (pp) ;
  pp = expand (pp) ;
  showPol(pp) ;
  printf ("DONE %s\n\n", title) ;

  return pp ;
} /* Z2_AA_loopB */

/***********************************************************************************************************************************************/

static POLYNOME Z2_AA_loopHB (const char *title) 
{    
  firstDummyIndex = 'a' ;

  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;
  char n = newDummyIndex () ;
  int mm1[4] = {0,1,0,0} ; /*  p  */
  int mm2[4] = {0,-1,0,0} ; /* -p  */

  POLYNOME p1 = prop_HB_H (1) ;   /* (1/(k+p)^2 */
  POLYNOME p2 = vertex_A_B_HB (a, c, d, mm1) ; /* (2k + p)_mu */prop_PsiRB_PsiR (1) ; /* (1/(k+p)^2 */
  POLYNOME p3 = prop_BB_B (c,d,e,f, 0) ;   /* (1/(k)^2 */
  POLYNOME p4 = vertex_A_H_BB (b, e, f, mm2) ; /*(2k + p)_mu */

  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(newMultiProduct (ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;

  if (1)
    {
      printf ("Expand\n") ;
      pp = expand (pp) ;
      showPol (pp) ;
    }
  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp, n) ;
  showPol(pp) ;
  pp = expand (pp) ;
  pp = expand (pp) ;
  showPol(pp) ;
  printf ("DONE %s\n\n", title) ;

  return pp ;
} /* Z2_AA_loopHB */

/***********************************************************************************************************************************************/

static POLYNOME Z2_HH_loopAH (const char *title) 
{
  firstDummyIndex = 'a' ;
    
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char n = newDummyIndex () ;
  int ppv[4] = {1,2,0,0} ; /* 2p + k : vertex */
  
  POLYNOME p1 = vertex_A_H_HB (b, ppv) ; /*(2k + p)_mu */
  POLYNOME p2 = prop_HB_H (1) ;   /* (1/(p+k)^2 */
  POLYNOME p3 = vertex_A_H_HB (a, ppv) ; /* (2k + p)_mu */prop_PsiRB_PsiR (1) ; /* (1/(k+p)^2 */
  POLYNOME p4 = prop_AL (a,b,0) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(newMultiProduct (ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp, n) ;
  showPol (pp) ;
  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_HH_loopAH */

/***********************************************************************************************************************************************/

static POLYNOME Z2_HH_loopAB (const char *title) 
{
  firstDummyIndex = 'a' ;
    
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;
  char n = newDummyIndex () ;
  int m1[4] = {-1,-1,0,0} ; /* -p - k : incoming A momentum */
  int m2[4] = {1,1,0,0} ; /* p + k : vertex */
  
  POLYNOME p1 = vertex_A_H_BB (a,c,d,m1) ; /*(2k + p)_mu */
  POLYNOME p2 = prop_BB_B (c,d,e,f,0) ;   /* (1/(p+k)^2 */
  POLYNOME p3 = vertex_A_B_HB (b,e,f, m2) ; /* (2k + p)_mu */prop_PsiRB_PsiR (1) ; /* (1/(k+p)^2 */
  POLYNOME p4 = prop_AL (b,a,1) ;   /* (1/(k)^2 */
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(newMultiProduct (ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;
  freeIndex (pp) ;
  
  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp, n) ;
  showPol (pp) ;
  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_BB_loopAB */

/***********************************************************************************************************************************************/

static POLYNOME Z2_BB_loopAB (const char *title) 
{
  firstDummyIndex = 'a' ;
    
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;
  char g = newDummyIndex () ;
  char h = newDummyIndex () ;
  char i = newDummyIndex () ;
  char j = newDummyIndex () ;
  char u = newDummyIndex () ;
  char v = newDummyIndex () ;
  char w = newDummyIndex () ;
  char x = newDummyIndex () ;
  
  char n = newDummyIndex () ;
  int m1[4] = {0,1,0,0} ; /*  p : external B */
  int m2[4] = {-1,0,0,0} ; /* -k  : looping B ,  the A propagator takes K=p*/


  POLYNOME p0 = newAG (a,b,u,v,0) ;
  POLYNOME p1 = vertex_A_B_BB (i,u,v,e,f,m1, m2) ; /*(2k + p)_mu */
  POLYNOME p2 = prop_BB_B (e,f,g,h,0) ;   /* (1/(p+k)^2 */
  POLYNOME p3 = vertex_A_B_BB (j,g,h,w,x, m2, m1) ; /* (2k + p)_mu */prop_PsiRB_PsiR (1) ; /* (1/(k+p)^2 */
  POLYNOME p4 = prop_AL (j,i,1) ;   /* (1/(k)^2 */
  POLYNOME p5 = newAG (w,x,c,d,0) ;
  POLYNOME ppp[] = {p0,p1,p2,p3,p4,p5,0} ;

  POLYNOME pp = contractIndices(newMultiProduct (ppp)) ;

  printf ("%s\n", title) ;
  showPol (pp) ;
  freeIndex (pp) ;
  
  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp, n) ;
  showPol (pp) ;
  printf ("DONE %s\n", title) ;

  return pp ;
} /* Z2_BB_loopAB */

/***********************************************************************************************************************************************/

static POLYNOME Z2_BB_loopAH (const char *title) 
{
  firstDummyIndex = 'a' ;
    
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;
  char g = newDummyIndex () ;
  char h = newDummyIndex () ;
  char i = newDummyIndex () ;
  char j = newDummyIndex () ;
  char n = newDummyIndex () ;
  int m1[4] = {1,0,0,0} ; /*  k : vertex */
  int m2[4] = {-1,0,0,0} ; /* -k : vertex */


  if (0)
    {
      POLYNOME px = newEpsilon (a,b,c,d) ;
      POLYNOME py = newEpsilon (c,d,e,f) ;
      POLYNOME pz = newEpsilon (e,f,g,h) ;
      POLYNOME pppx[] = {px, py, pz, 0} ;
      
      POLYNOME ppx = newMultiProduct (pppx) ;
      showPol (ppx) ;
      ppx = expand (ppx) ;
      showPol (ppx) ;
      
      showPol (ppx) ;
      ppx = expand (ppx) ;
      showPol (ppx) ;
    }
 
  
  /*  POLYNOME p0 = newAG(a,b,g,h,1) ; */
  POLYNOME p1 = vertex_A_H_BB (e,g,h,m1) ; /*(2k + p)_mu */
  POLYNOME p2 = prop_HB_H (1) ;   /* (1/(p+k)^2 */
  POLYNOME p3 = vertex_A_B_HB (f,i,j, m2) ; /* (2k + p)_mu */prop_PsiRB_PsiR (1) ; /* (1/(k+p)^2 */
  POLYNOME p4 = prop_AL (f,e,0) ;   /* (1/(k)^2 */
  /*   POLYNOME p5 = newAG(c,d,i,j,-1) ; */
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = newMultiProduct (ppp) ;
  showPol (pp) ;
  pp = contractIndices(pp) ;

  printf ("%s\n", title) ;
  showPol (pp) ;
  freeIndex (pp) ;
  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp, n) ;
  showPol (pp) ;

  pp = expand (pp) ;
  showPol (pp) ;
  p1 = pp->p1 ;
  p2=pp->p2 ;
  showPol (p1) ;
  showPol (p2) ;
  printf ("DONE %s\n", title) ;

  POLYNOME p02 = newAG(a,b,g,h,0) ;
  POLYNOME p52 = newAG(i,j,c,d,0) ;
  POLYNOME ppp2[] = {p02,p2,p52,0} ;
  pp = newMultiProduct (ppp2) ;
  showPol (pp) ;
  pp = contractIndices(pp) ;
  showPol (pp) ;
  POLYNOME ppE = expand (pp) ;
  showPol (pp) ;
  pp = expand (pp) ;
  showPol (pp) ;
  pp = expand (pp) ;
  showPol (pp) ;
  pp = expand (pp) ;
  showPol (pp) ;
  exit (0) ;
  
  p02 = newAG(a,b,g,h,-1) ;
  p52 = newAG(i,j,c,d,1) ;
  POLYNOME ppp3[] = {p02,p1,p52,0} ;
  pp = contractIndices(newMultiProduct (ppp3)) ;
  showPol (pp) ;
  POLYNOME ppG = expand (pp) ;
  showPol (ppG) ;

  pp = newSum (ppE, ppG) ;
  exit (0) ;
  
  p02 = newAG(a,b,g,h,2) ;
  p52 = newAG(i,j,c,d,0) ;
  POLYNOME ppp4[] = {p02,p1,p52,0} ;
  pp = contractIndices(newMultiProduct (ppp4)) ;
  showPol (pp) ;
  pp = expand (pp) ;
  showPol (pp) ;

  p02 = newAG(a,b,g,h,0) ;
  p52 = newAG(i,j,c,d,2) ;
  POLYNOME ppp5[] = {p02,p1,p52,0} ;
  pp = contractIndices(newMultiProduct (ppp5)) ;
  showPol (pp) ;

  
  
  return pp ;
} /* Z2_BB_loopAH */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME Z2_PsiL__B_Psi (void)
{
  char mu = newDummyIndex () ;
  char nu = newDummyIndex () ;
  char rho = newDummyIndex () ;
  char sigma = newDummyIndex () ;

  POLYNOME p1 = vertex_BB_PsiL_PsiRB (mu,nu) ;
  POLYNOME p2 = prop_BB_B (mu,nu,rho,sigma, 0) ;   /* (1/(k)^2 */
  POLYNOME p3 = prop_PsiRB_PsiR (1) ; /* (1/(k+p)^2 */
  POLYNOME p4 = vertex_B_PsiR_PsiLB (rho,sigma) ;
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(newMultiProduct (ppp)) ;
  
  printf ("###############Z2 Psi left avec loop B_mu_nu\n") ;
  showPol (pp) ;

  printf ("Expand\n") ;
  pp = expand (pp) ;
  showPol (pp) ;

  printf ("integrate\n") ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  printf ("pauli cleanup\n") ;
  pp = pauliCleanUp (pp) ;
  showPol (pp) ;
  printf ("### Z2 Psi left avec loop B_mu_nu, expect ::  je_sais_pas * p-slash\n") ;
  showPol (pp) ;
  printf ("Z2_  done\n\n") ;

  firstDummyIndex = 'a' ;
  p1 = newSigma (mu) ; 
  p1->tt.sigma[1] = nu ;
  p3 = newSigB (nu) ; 
  p3->tt.sigB[1] = rho ;
  p1->tt.mm[0][0] = mu ;
  p1->tt.mm[0][1] = rho ;
  p1->tt.denom[0] = 2 ;
  p2 = prop_PsiLB_PsiL (1) ; /* (1/(k+p)^2 */
  ppp[0] = p1 ;
  ppp[1] = p2 ;
  ppp[2] = p3 ;
  ppp[3] = 0  ;
  pp = newMultiProduct (ppp) ;
  showPol(pp) ;
  pp = expand (pp) ;
  showPol (pp) ;
  pp = contractIndices (pp) ;
  showPol(pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;
  pp = contractIndices (pp) ;
  showPol(pp) ;

  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME Z2_PsiL__A_Psi (void)
{
  char mu = newDummyIndex () ;
  char nu = newDummyIndex () ;
  BOOL debug = FALSE ;

  POLYNOME p1 = vertex_A_PsiL_PsiLB (mu) ;
  POLYNOME p2 = prop_AL (mu,nu, 0) ;   /* (1/(k)^2 */
  POLYNOME p3 = prop_PsiLB_PsiL (1) ; /* (1/(k+p)^2 */
  POLYNOME p4 = vertex_A_PsiL_PsiLB (nu) ;
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(newMultiProduct (ppp)) ;
  
  printf ("###############Z2 Psi left avec loop A_mu\n") ;
  showPol (pp) ;

  pp = expand (pp) ;
  if (debug) showPol (pp) ;

  pp = dimIntegral (pp) ;
  printf ("### Z2 Psi left avec loop A_mu, expect  p-slash\n") ;
  showPol (pp) ;
  pp = pauliCleanUp (pp) ;
  showPol (pp) ;
  return pp ;
}

/***********************************************************************************************************************************************/

static POLYNOME Z2_PsiL__H_Psi (void)
{
  BOOL debug = FALSE ;
  POLYNOME p1 = vertex_HB_PsiL_PsiRB () ;
  POLYNOME p2 = prop_HB_H (0) ;          /* (1/(k)^2 */ 
  POLYNOME p3 = prop_PsiRB_PsiR (1) ;   /* (1/(k+p)^2 */ 
  POLYNOME p4 = vertex_H_PsiR_PsiLB () ;
  POLYNOME ppp[] = {p1,p2,p3,p4,0} ;

  POLYNOME pp = contractIndices(newMultiProduct (ppp)) ;
  
  printf ("############# Z2 Psi left avec loop Phi\n") ;
  showPol (pp) ;
  
  expand (pp) ;

  if (debug) showPol (pp) ;
  pp = dimIntegral (pp) ;
  printf ("### Z2 Psi left avec loop A_mu, expect  1/2 p-slash\n") ;  
  showPol (pp) ;
  pp = pauliCleanUp (pp) ;
  showPol (pp) ;

  return pp ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME Z2_H__AH (void) 
{
  char mu = newDummyIndex () ;
  char nu = newDummyIndex () ;
  POLYNOME pp, ppp[6] ;
  int mm1[4] = {1,2,0,0} ; /* 2p + k */
  int mm2[4] = {1,2,0,0} ; /* 2p + k */

  ppp[0] = prop_AL (mu, nu, 0) ; /* 1/k^2 */
  ppp[1] = vertex_A_H_HB (mu, mm1) ;
  ppp[2] = prop_HB_H (1) ; /* 1/(p+k)^2 */
  ppp[3] = vertex_A_H_HB (nu, mm2) ;
  ppp[4] = 0 ;

  pp = newMultiProduct (ppp) ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  printf ("### Z2 H with transient A \n") ;  
  showPol (pp) ;
  pp = reduceIndices (pp) ;
  showPol (pp) ;
  pp = squareMomentaCleanUp (pp, mu) ;
  showPol(pp) ;
  
  return pp ;
} /* Z2_H__AH */

/***********************************************************************************************************************************************/
/* check some integrals */
static void Thooft (void) 
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;
  char g = newDummyIndex () ;
  char h = newDummyIndex () ;

  POLYNOME pp, qq, rr ;

  pp = newPolynome () ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  pp->tt.denom[0] = 1 ;
  pp->tt.denom[1] = 1 ;
  printf ("### Int { 1/k^2 (k+p)^2, expect I1\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;

  pp = newPolynome () ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome () ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  strcpy (qq->tt.g, pp->tt.mm[0]) ;
  pp->tt.denom[0] = 1 ;
  pp->tt.denom[1] = 2 ;
  printf ("### Int { k^a k^b/k^2 (k+p)^4, expect g_ab/4\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = newProduct (pp, qq) ;
  rr = expand (rr) ;
  showPol (rr) ;


  pp = newPolynome () ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome () ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  pp->tt.mm[0][2] = c ;
  pp->tt.mm[0][3] = d ;
  strcpy (qq->tt.g, pp->tt.mm[0]) ;
  pp->tt.denom[0] = 2 ;
  pp->tt.denom[1] = 2 ;
  printf ("### Int { k^a k^b k^c k^d/k^4 (k+p)^4, expect g_ab g_cd/24\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = newProduct (pp, qq) ;
  rr = expand (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;

  pp = newPolynome () ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome () ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  pp->tt.mm[0][2] = c ;
  pp->tt.mm[0][3] = d ;
  pp->tt.mm[0][4] = e ;
  pp->tt.mm[0][5] = f ;
  strcpy (qq->tt.g, pp->tt.mm[0]) ;
  pp->tt.denom[0] = 2 ;
  pp->tt.denom[1] = 3 ;
  printf ("### Int { k^a k^b k^c k^d k^e k^f/k^4 (k+p)^6, expect g_ab g_cd/192\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = newProduct (pp, qq) ;
  showPol (rr) ;
  rr = expand (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;

  pp = newPolynome () ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome () ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  pp->tt.mm[0][2] = c ;
  pp->tt.mm[0][3] = d ;
  pp->tt.mm[0][4] = e ;
  pp->tt.mm[0][5] = f ;
  pp->tt.mm[0][6] = g ;
  pp->tt.mm[0][7] = h ;
  strcpy (qq->tt.g, pp->tt.mm[0]) ;
  pp->tt.denom[0] = 3 ;
  pp->tt.denom[1] = 3 ;
  printf ("### Int { k^a k^b k^c k^d k^e k^f k^g k^h/k^6 (k+p)^6, expect g_ab g_cd/1920\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = newProduct (pp, qq) ;
  rr = expand (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;
 

 pp = newPolynome () ;
  pp->tt.z = 1 ;
  pp->tt.type = 1 ;
  qq = newPolynome () ;
  qq->tt.z = 1 ;
  qq->tt.type = 1 ;

  pp->tt.mm[0][0] = a ;
  pp->tt.mm[0][1] = b ;
  pp->tt.mm[0][2] = c ;
  pp->tt.mm[0][3] = d ;
  pp->tt.mm[0][4] = e ;
  strncpy (qq->tt.g, pp->tt.mm[0],4) ;
  pp->tt.denom[0] = 1 ;
  pp->tt.denom[1] = 1 ;
  pp->tt.denom[2] = 1 ;
  pp->tt.denom[3] = 1 ;
  printf ("### Int { k^a k^b k^c k^d k^e /k^2 (k+p)^2 (k+p+q)^2 (k+p+q+r)^2, expect p_a g_bc g_de/ ?\n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol (pp) ;
  rr = newProduct (pp, qq) ;
  showPol (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;
  rr = expand (rr) ;
  showPol (rr) ;

  return ;
} /* tHooft */

/***********************************************************************************************************************************************/
/* check self duality projectors */
static POLYNOME Hodge (void) 
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;
  char g = newDummyIndex () ;
  char h = newDummyIndex () ;
  char i = newDummyIndex () ;
  char j = newDummyIndex () ;
  char mu = newDummyIndex () ;
  char nu = newDummyIndex () ;
  char rho = newDummyIndex () ;
  char sig = newDummyIndex () ;

  POLYNOME pp, ppp[6], pP1, pP2, pM2,pG1, pG2 ;
  /*
    int mm1[4] = {1,1,0,0} ; //
    int mm2[4] = {-1,-1,0,0} ; //  p + k 
  */

  if (1)
    {
      pP1 = newAG (mu,nu,a,b,1) ;
      /*  pM1 = newAG (mu,nu,a,b,-1) ; */
      pP2 = newAG (rho,sig,a,b,1) ;
      pM2 = newAG (rho,sig,a,b,-1) ;
      pG1 = newAG (mu,nu,a,b,0) ;
      pG2 = newAG (rho,sig,a,b,0) ;
    }
  if (0)
    {
      pp = newG (mu,nu) ;
      strcpy (pp->tt.eps, "abcdabcd") ;
      showPol (pp) ;
      pp = contractEpsilon (pp) ;
      showPol (pp) ;
      
      pp = newG (mu,nu) ;
      strcpy (pp->tt.eps, "abcdabdc") ;
      showPol (pp) ;
      pp = contractEpsilon (pp) ;
      showPol (pp) ;
      
      pp = newG (mu,nu) ;
      strcpy (pp->tt.eps, "abcdacbd") ;
      showPol (pp) ;
      pp = contractEpsilon (pp) ;
      showPol (pp) ;
      
      pp = newG (mu,nu) ;
      strcpy (pp->tt.eps, "abcdabce") ;
      showPol (pp) ;
      pp = contractEpsilon (pp) ;
      showPol (pp) ;
      
      pp = newG (mu,nu) ;
      strcpy (pp->tt.eps, "abcdaebc") ;
      showPol (pp) ;
      pp = contractEpsilon (pp) ;
      showPol (pp) ;
      
      
      pp = newG (mu,nu) ;
      strcpy (pp->tt.eps, "abcdabec") ;
      showPol (pp) ;
      pp = contractEpsilon (pp) ;
      showPol (pp) ;
    }
  
  if (0)
    {
      ppp[0] = newG ('a','b') ;
      ppp[1] = 0 ;
      pp = newMultiSum (ppp) ;
      showPol (pp) ;
      pp = sortPol (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      ppp[0] = newG ('a','b') ;
      ppp[1] = newG ('c','d') ;
      ppp[2] = 0 ;
      pp = newMultiSum (ppp) ;
      showPol (pp) ;
      pp = sortPol (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      ppp[0] = newG ('a','b') ;
      ppp[1] = newG ('c','d') ;
      ppp[2] = newG ('a','b') ;
      ppp[3] = 0 ;
      pp = newMultiSum (ppp) ;
      showPol (pp) ;
      pp = sortPol (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

    }
  if (1)
    {
      ppp[0] = newG ('a','b') ;
      ppp[1] = newG ('c','d') ;
      ppp[2] = newG ('a','b') ;
      ppp[3] = newG ('c','d') ;
      ppp[4] = 0 ;
      pp = newMultiSum (ppp) ;
      showPol (pp) ;
      pp = sortPol (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      exit (0) ;

      ppp[0] = newG ('a','b') ;
      ppp[1] = newG ('c','d') ;
      ppp[2] = newG ('a','b') ;
      ppp[3] = newG ('c','d') ;
      ppp[4] = newG ('a','b') ;
      ppp[5] = 0 ;
      pp = newMultiSum (ppp) ;
      showPol (pp) ;
      pp = sortPol (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
    }

  if (1)
    {
      printf ("Is PP a projector: compute PP^2\n") ;
      ppp[0] = pP1 ;
      ppp[1] = pP2 ;
      ppp[2] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;

      showPol (pp) ;
  
      printf ("Is PG a projector: compute PG PG\n") ;
      ppp[0] = pG1 ;
      ppp[1] = pG2 ;
      ppp[2] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = contractProducts (pp) ;
      contractIndices (pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      printf ("Is PP * PM zero: compute PP PM\n") ;
      ppp[0] = pP1 ;
      ppp[1] = pM2 ;
      ppp[2] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = contractProducts (pp) ;
      contractIndices (pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      printf ("Test : compute PP PG\n") ;
      ppp[0] = pP1 ;
      ppp[1] = pG2 ;
      ppp[2] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      
      printf ("Test : compute PG eps\n") ;
      ppp[0] = pG1 ;
      ppp[1] = newEpsilon (rho,sig,a,b) ;
      ppp[2] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute B propagator\n") ;
      ppp[0] = newAG (a,b,e,f,-1) ; 
      ppp[1] = prop_BB_B (e,f,g,h,0) ;
      ppp[2] = newAG (g,h,c,d,1) ; 
      ppp[3] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = reduceIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = reduceIndices (pp) ;
      pp = expand (pp) ;
      pp = reduceIndices (pp) ;
      pp = expand (pp) ;
      pp = sortPol (pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      
      showPol (pp) ;
    }

  if (1)
    {
      printf ("Test : compute eps eps\n") ;
      ppp[0] = newEpsilon (mu,nu,a,b) ;
      ppp[1] = newEpsilon (rho,b,a,sig) ;
      ppp[0]->tt.z = I/4 ;
      ppp[1]->tt.z = I/4 ;
      ppp[2] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;

      showPol (pp) ;
    }
  
  if (1)
    {
      printf ("Test : compute chiral projectors PP PP = PP\n") ;
      ppp[0] = newAG (a,b,c,d,1) ;
      ppp[1] = newAG (c,d,e,f,1) ;
      ppp[2] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute chiral projectors PP PP PP = PP\n") ;
      ppp[0] = newAG (a,b,c,d,1) ;
      ppp[1] = newAG (c,d,e,f,1) ;
      ppp[2] = newAG (e,f,g,h,1) ;
      ppp[3] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute chiral projectors PP PP PP PP = PP\n") ;
      ppp[0] = newAG (a,b,c,d,1) ;
      ppp[1] = newAG (c,d,e,f,1) ;
      ppp[2] = newAG (e,f,g,h,1) ;
      ppp[3] = newAG (g,h,i,j,1) ;
      ppp[4] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute chiral projectors PP PM = 0\n") ;
      ppp[0] = newAG (a,b,c,d,1) ;
      ppp[1] = newAG (c,d,e,f,-1) ;
      ppp[2] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute chiral projectors PM PP = 0\n") ;
      ppp[0] = newAG (a,b,c,d,-1) ;
      ppp[1] = newAG (c,d,e,f,1) ;
      ppp[2] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : compute chiral projectors PP PM PP = PP\n") ;
      ppp[0] = newAG (a,b,c,d,1) ;
      ppp[1] = newAG (c,d,e,f,-1) ;
      ppp[2] = newAG (e,f,g,h,1) ;
      ppp[3] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractIndices (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Test : double epsilon, no common index\n") ;
      pp = newG (a,b) ;
      strcpy (pp->tt.eps, "cdefghij") ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;
      pp = contractEpsilon (pp) ;
      showPol (pp) ;
    }
  exit (0) ;
  return pp ;
} /* Hodge */

/***********************************************************************************************************************************************/

static POLYNOME Z2_H__AB (void) 
{
  firstDummyIndex = 'a' ;
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;

  POLYNOME pp, ppp[6] ;
  int mm1[4] = {1,1,0,0} ; /* p + k */
  int mm2[4] = {-1,-1,0,0} ; /* p + k */

  if (1)
    {
      ppp[0] = prop_AL (a, b, 1) ; /* 1/(p+k^2 */
      ppp[1] = vertex_A_B_HB (a, c, d, mm1) ;
      ppp[2] = prop_BB_B (c,d,e,f,0) ; /* kk/(k)^4 */
      ppp[3] = vertex_A_H_BB (b, e,f, mm2) ;
      ppp[4] = 0 ; 
    }
  else
    {
      ppp[0] = vertex_A_B_HB (a, c, d, mm1) ;
      ppp[1] = prop_BB_B (c,d,e,f,0) ; /* kk/(k)^4 */
      ppp[1] = newAG (c,d,e,f,-1) ;
      ppp[2] = 0 ; 
    }

  printf ("### Z2 H with new vertex ABH \n") ;  
  pp = newMultiProduct (ppp) ;
  showPol (pp) ;
  pp = expand (pp) ;
  pp = reduceIndices (pp) ;
  showPol (pp) ;
  pp = expand (pp) ;
  pp = reduceIndices (pp) ;
  showPol (pp) ;
  if (pp) pp = dimIntegral (pp) ;

  pp = reduceIndices (pp) ;
  pp = expand(pp) ;
  showPol (pp) ;
  if (pp) pp = contractIndices(pp) ;
  showPol (pp) ;
  if (pp) pp = squareMomentaCleanUp (pp, a) ;
  showPol(pp) ;
  pp = reduceIndices (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z2_H__AB */

/***********************************************************************************************************************************************/

static POLYNOME Z2_B__AB (void) 
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;
  char g = newDummyIndex () ;
  char h = newDummyIndex () ;
  char i = newDummyIndex () ;
  char j = newDummyIndex () ;
  char k = newDummyIndex () ;
  char l = newDummyIndex () ;
  char mu = newDummyIndex () ;
  char nu = newDummyIndex () ;

  POLYNOME pp, ppp[12] ;
  int mm1[4] = {0,1,0,0} ;   /* p + k */
  int mm2[4] = {1,0,0,0} ; /* p + k */
  int mm3[4] = {0,1,0,0} ;   /* p */

  if (1)
    { 
      ppp[0] = newAG (a,b,i,j, -1) ;
      ppp[0]->tt.z = 384 ;
      ppp[1] = prop_AL (mu, nu, 1) ; /* 1/(p+k^2 */
      ppp[2] = vertex_A_B_BB (mu, i, j, e, f, mm1, mm2) ;
      ppp[3] = prop_BB_B (e,f,g,h,0) ; /* kk/(k)^4 */
      ppp[4] = vertex_A_B_BB (nu, g, h, k, l, mm2, mm3) ;
      ppp[5] = newAG (k,l,c,d, 1) ;
      ppp[6] = newScalar (384) ;
      ppp[7] = 0 ;

      pp = newMultiProduct (ppp) ;
      printf ("########### 384 Z2 B with transient A TRES FAUX resultat en p^4 au lieu de p^2 \n") ;  
      showPol (pp) ;
      pp = expand (pp) ;
      pp = dimIntegral (pp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      pp = squareMomentaCleanUp (pp, mu) ;
      showPol(pp) ;
      pp = reduceIndices (pp) ;
      showPol(pp) ;
    }  
  return pp ;
} /* Z2_B__AB */

/***********************************************************************************************************************************************/

static POLYNOME Z2_B__AH (void) 
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;
  char g = newDummyIndex () ;
  char h = newDummyIndex () ;
  char mu = newDummyIndex () ;
  char nu = newDummyIndex () ;
 
  POLYNOME pp, ppp[12] ;
  int mm1[4] = {1,0,0,0} ; /* k */
  int mm2[4] = {-1,0,0,0} ; /* k */

  if (0)
    {
      newDummyIndex () ;
      newDummyIndex () ;
    }
  if (1)
    {
      if (1)
	{
	  ppp[0] = newScalar (1) ; /* newAG (a,b,e,f, -1) ; */
	  ppp[1] = prop_AL (mu, nu, 0) ; /* 1/(p+k^2 */
	  ppp[2] = vertex_A_B_HB (mu, a,b, mm1) ;
	  ppp[3] = prop_HB_H (1) ; /* 1/(k+p)^2 */
	  ppp[4] = vertex_A_H_BB (nu, c,d, mm2) ;
	  ppp[5] = 0 ;
	}
      else
	{
	  ppp[0] = vertex_A_H_BB (nu, g,h, mm2) ;
	  ppp[1] = newAG (a,b,c,d, -1) ;
	  ppp[2] =0 ;
	}

      printf ("### 384 Z2 B with new vertex ABH \n") ;  
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol (pp) ;

      printf ("Integrate\n") ;
      pp = dimIntegral (pp) ;
      showPol(pp) ;
      pp = contractIndices (pp) ;
      if (1) pp = reduceIndices (pp) ;
      showPol(pp) ;

      printf ("Project\n") ;
      e = newDummyIndex () ;
      f = newDummyIndex () ;
      g = newDummyIndex () ;
      h = newDummyIndex () ;
      ppp[0] = newAG (e,f,'a','b', 1) ;
      ppp[1] = pp ;
      ppp[2] = newAG (g,h,'c','d', -1) ;
      ppp[3] = 0 ;
      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      if (1) pp = reduceIndices (pp) ;
      pp = contractIndices (pp) ;
      pp = expand (pp) ;
      if (1) pp = reduceIndices (pp) ;
      showPol (pp) ;
    }  
  return pp ;
} /* Z2_B__AH */

/***********************************************************************************************************************************************/

static POLYNOME Z2_A__BH (void) 
{
  firstDummyIndex = 'a' ;
  
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;

  POLYNOME pp, ppp[6] ;
  int mm1[4] = {0,1,0,0} ; /*  p  */
  int mm2[4] = {0,-1,0,0} ; /* -p  */

  if (1)
    {
      ppp[0] = vertex_A_B_HB (a,c,d, mm1) ;
      ppp[1] = prop_HB_H (1) ;  /* 1/(p+k^2 */
      ppp[2] = vertex_A_H_BB (b,e,f, mm2) ;
      ppp[3] = prop_BB_B (c,d,e,f, 0) ; /* kk/(k)^4 */
      ppp[4] = 0 ;

      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = dimIntegral (pp) ;
      printf ("### ** Z2 A with new vertex ABH \n") ;  
      showPol (pp) ;
      pp = expandDo (pp,2) ;
      showPol(pp) ;
      pp = expand(pp) ;
      showPol(pp) ;
      pp = squareMomentaCleanUp (pp, 0) ;
      showPol(pp) ;
      pp = expand (pp) ;
      pp = expand (pp) ;
      showPol(pp) ;
      printf ("DONE Z2_A__BH\n\n") ;
    }  
  return pp ;
} /* Z2_A__BH */

/***********************************************************************************************************************************************/

static POLYNOME Z2_H__PsiL__Psi (void) 
{
  POLYNOME ppp[6], pp ;

  if (1)
    {
      ppp[0] = newScalar (-1) ; /* Fermion loop */
      ppp[1] = prop_PsiRB_PsiR (0) ;   /* (1/(k)^2 */ 
      ppp[2] = vertex_H_PsiR_PsiLB () ;
      ppp[3] = prop_PsiLB_PsiL (1) ;   /* (1/(k+p+q)^2 */ 
      ppp[4] = vertex_HB_PsiL_PsiRB () ;
      ppp[5] = 0 ;

      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = dimIntegral (pp) ;
      printf ("### Z2 H with Psi loop\n") ;  
      showPol (pp) ;
      pp = pauliTrace (pp) ;
      showPol (pp) ;
      pp = expand (pp) ; 
      pp = squareMomentaCleanUp (pp, 0) ;
      pp = expand (pp) ;
      showPol(pp) ;
    }  
  return pp ;
} /* Z2_H__PsiL__Psi */

/***********************************************************************************************************************************************/

static POLYNOME Z2_A__PsiL__Psi (void) 
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  POLYNOME ppp[6], pp ;

  if (1)
    {
      ppp[0] = newScalar (-1) ; /* Fermion loop */   
      ppp[1] = prop_PsiRB_PsiR (0) ; /* 1/k^2 */
      ppp[2] = vertex_A_PsiR_PsiRB (b) ;
      ppp[3] = prop_PsiRB_PsiR (1) ; /* (1/(k+p+q)^2 */ 
      ppp[4] = vertex_A_PsiR_PsiRB (a) ;
      ppp[5] = 0 ;

      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol(pp) ;
      pp = pauliTrace (pp) ;
      showPol(pp) ;
      pp = expand (pp) ;
      showPol(pp) ;

      pp = dimIntegral (pp) ;
      printf ("### Z2 A with Psi loop\n") ;  
      showPol (pp) ;

      showPol (pp) ;
      pp = expand (pp) ; 
      showPol(pp) ;
      pp = squareMomentaCleanUp (pp, 0) ;
      showPol(pp) ;
      pp = expand (pp) ;
      showPol(pp) ;
      pp = sortPol(pp) ;
      showPol(pp) ;
    }  
  return pp ;
} /* Z2_A__PsiL__Psi */

/***********************************************************************************************************************************************/

static POLYNOME Z2_B__PsiL__Psi (void) 
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  POLYNOME ppp[6], pp ;

  if (1)
    {
      ppp[0] = newScalar (-1) ; /* Fermion loop */
      ppp[1] = prop_PsiRB_PsiR (0) ; /* 1/k^2 */
      ppp[2] = vertex_B_PsiR_PsiLB (a,b) ;
      ppp[3] = prop_PsiLB_PsiL (1) ; /* (1/(k+p+q)^2 */ 
      ppp[4] = vertex_BB_PsiL_PsiRB (c,d) ;
      ppp[5] = 0 ;

      pp = newMultiProduct (ppp) ;
      showPol (pp) ;
      pp = expand (pp) ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      pp = expand (pp) ;
      showPol(pp) ;

      pp = pauliTrace (pp) ;
      showPol(pp) ;
      pp = expand (pp) ;
      printf ("### Z2 B with Psi loop\n") ;  
      showPol (pp) ;

      showPol (pp) ;
      pp = expand (pp) ; 
      showPol(pp) ;
      pp = squareMomentaCleanUp (pp, 0) ;
      showPol(pp) ;
      pp = reduceIndices (pp) ;
      showPol(pp) ;
      pp = expand (pp) ;
      showPol(pp) ;
      pp = sortPol(pp) ;
      showPol(pp) ;
    }  

  return pp ;
} /* Z2_B__PsiL__Psi */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME Z3_A_HB_H__psiLoop (void)
{
  char mu = newDummyIndex () ;
  POLYNOME p0 = newP (mu) ;
  /* calcul plus simple, on calcule successivement q=r=0, r=p=0, p=q=0,  et on additionne */

  /* set q == 0 */
  POLYNOME p11 = prop_PsiRB_PsiR (0) ;   /* (1/(k)^2 */ 
  POLYNOME p12 = vertex_A_PsiR_PsiRB (mu) ;
  POLYNOME p13 = prop_PsiRB_PsiR (1) ;   /* (1/(k+p)^2 */ 
  POLYNOME p14 = vertex_H_PsiR_PsiLB () ;
  POLYNOME p15 = prop_PsiLB_PsiL (2) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p16 = vertex_HB_PsiL_PsiRB () ;
  POLYNOME pppP[] = {p11,p12,p13,p14,p15,p16,0} ;

  POLYNOME PP = contractIndices(newMultiProduct (pppP)) ;

  printf ("############# Z3 A HB B with psi loop\n") ;
  showPol (PP) ;

  PP = dimIntegral (PP) ;
  showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;
  PP = pauliTrace (PP) ;
  showPol (PP) ;
  PP = momentaCleanUp (PP, mu) ;
  printf ("### Z3 A HB H with psi loop, expect (-1) Trace[1]*(p+2q)_a = -2p - 4q\n") ;
  showPol (PP) ;


  return p0 ;
} /* Z3_A_HB_H__psiLoop */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_BB_B__psiLoop (void)
{
  char mu = newDummyIndex () ;
  POLYNOME p0 = newP (mu) ;
  return p0 ;
} /* Z3_A_BB_B__psiLoop */

/***********************************************************************************************************************************************/

static POLYNOME Z3_A_A_A__psiLoop (void)
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char n = newDummyIndex () ;

  /* calcul plus simple, on calcule successivement q=r=0, r=p=0, p=q=0,  et on additionne */

  /* set q == 0 */
  POLYNOME p11 = vertex_A_PsiR_PsiRB (a) ;
  POLYNOME p12 = prop_PsiRB_PsiR (0) ;   /* (1/(k)^2 */ 
  POLYNOME p13 = vertex_A_PsiR_PsiRB (c) ;
  POLYNOME p14 = prop_PsiRB_PsiR (2) ;   /* (1/(k+p+q)^2 */
  POLYNOME p15 = vertex_A_PsiR_PsiRB (b) ;
  POLYNOME p16 = prop_PsiRB_PsiR (1) ;   /* (1/(k+p)^2 */ 



  POLYNOME pppP[] = {p11,p12,p13,p14,p15,p16,0} ;

  POLYNOME PP = contractIndices(newMultiProduct (pppP)) ;

  printf ("############# Z3 A A A with psi loop\n") ;
  showPol (PP) ;

  PP = dimIntegral (PP) ;
  showPol (PP) ;
  printf ("# trace the pauli matrix :\n ") ;
  PP = pauliTrace (PP) ;
  showPol (PP) ;
  PP = momentaCleanUp (PP, n) ;
  PP = expand (PP) ;
  printf ("### Z3 A A A with psi loop, expect g_bc (p+2q)_rho + .. + ..\n") ;
  showPol (PP) ;


  return PP ;
} /* Z3_A_A_A__psiLoop */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME ZF_H_PsiR_PsiLB (void)
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;

  int mm[4] = {1,0,0,0} ;

  POLYNOME p1 = vertex_A_PsiR_PsiRB (b) ;
  POLYNOME p2 = prop_PsiRB_PsiR (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_B_PsiR_PsiLB (e,f) ;

  POLYNOME p4 = prop_BB_B (c,d,e,f, 0) ; /* 1/k^2 */
  POLYNOME p5 = prop_AL (a,b,0) ; /* 1/k^2 */
  POLYNOME p6 = vertex_A_B_HB (a,c,d,mm) ;
  POLYNOME ppp[7] = {p1,p2,p3,p4,p5,p6,0} ;
  POLYNOME pp = newMultiProduct (ppp) ;


  printf ("Z3 New vertex A/H/B contrib to vertex phi-psi-psi : \n") ;
  showPol(p6) ;
  showPol(pp) ;
  pp = expand (pp) ;
  showPol(pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* ZF_H_PsiR_PsiLB */

/***********************************************************************************************************************************************/

static POLYNOME Z3_H_PsiR_PsiLB (void)
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;

  POLYNOME p1 = vertex_A_PsiR_PsiRB (a) ;
  POLYNOME p2 = prop_PsiRB_PsiR (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_H_PsiR_PsiLB () ;
  POLYNOME p4 = prop_PsiLB_PsiL (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = vertex_A_PsiL_PsiLB (b) ;
  POLYNOME p6 = prop_AL (a,b,0) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = newMultiProduct (ppp) ;
  printf ("Z3 Classic vertex H psi with A under: expect 4 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;


  return pp ;
} /* Z3_H_PsiR_PsiLB */

/***********************************************************************************************************************************************/

static POLYNOME ZF_B_PsiR_PsiLB (void)
{
  char mu = newDummyIndex () ;
  char nu = newDummyIndex () ;
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;

  int mm[4] = {1,0,0,0} ;

  POLYNOME p1 = vertex_A_PsiR_PsiRB (b) ;
  POLYNOME p2 = prop_PsiRB_PsiR (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_H_PsiR_PsiLB () ;

  POLYNOME p4 = prop_HB_H (0) ; /* 1/k^2 */
  POLYNOME p5 = prop_AL (a,b,0) ; /* 1/k^2 */
  POLYNOME p6 = vertex_A_B_HB (a,mu,nu,mm) ;
  POLYNOME ppp[7] = {p1,p2,p3,p4,p5,p6,0} ;
  POLYNOME pp = newMultiProduct (ppp) ;


  printf ("Z3 New vertex A/H/B contrib to vertex B_ab -psi-psi : \n") ;
  showPol(p6) ;
  showPol(pp) ;
  pp = expand (pp) ;
  showPol(pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;
  return pp ;
} /* ZF_B_PsiR_PsiLB */

/***********************************************************************************************************************************************/

static POLYNOME Z3_B__AB__PsiR_PsiLB (void)
{
  char mu = newDummyIndex () ;
  char nu = newDummyIndex () ;
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;
  char f = newDummyIndex () ;

  int mm1[4] = {0,0,0,0} ;
  int mm2[4] = {1,0,0,0} ;


  POLYNOME p1 = prop_PsiRB_PsiR (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p2 = vertex_A_PsiR_PsiRB (a) ;
  POLYNOME p3 = prop_AL (a,b,0) ; /* 1/k^2 */
  POLYNOME p4 = vertex_A_B_BB (a,mu,nu,c,d,mm1,mm2) ;
  POLYNOME p5 = prop_BB_B (c,d,e,f, 0) ; /* 1/k^2 */
  POLYNOME p6 = vertex_BB_PsiL_PsiRB (e,f) ;
  POLYNOME ppp[7] = {p1,p2,p3,p4,p5,p6,0} ;
  POLYNOME pp = newMultiProduct (ppp) ;


  printf ("Z3 Classic vertex A/B contrib to vertex B_ab -psi-psi : \n") ;
  showPol(p6) ;
  showPol(pp) ;
  pp = expand (pp) ;
  showPol(pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;
  return pp ;
} /* ZF_B_PsiR_PsiLB */

/***********************************************************************************************************************************************/

static POLYNOME Z3_B_PsiR_PsiLB (void)
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;

  POLYNOME p1 = vertex_A_PsiR_PsiRB (c) ;
  POLYNOME p2 = prop_PsiRB_PsiR (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_B_PsiR_PsiLB (a,b) ;
  POLYNOME p4 = prop_PsiLB_PsiL (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = vertex_A_PsiL_PsiLB (d) ;
  POLYNOME p6 = prop_AL (c,d,0) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = newMultiProduct (ppp) ;
  printf ("Z3 Classic vertex B psi with A under: expect 0 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3_H_PsiR_PsiLB */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static POLYNOME Z3H_A_PsiR_PsiLB (void)
{
  char a = newDummyIndex () ;

  POLYNOME p1 = vertex_H_PsiR_PsiLB () ;
  POLYNOME p2 = prop_PsiLB_PsiL (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_A_PsiL_PsiLB (a) ;
  POLYNOME p4 = prop_PsiLB_PsiL (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = vertex_HB_PsiL_PsiRB () ;
  POLYNOME p6 = prop_HB_H (0) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = newMultiProduct (ppp) ;
  printf ("Z3 Classic vertex A psi with H under: expect 0 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3H_A_PsiR_PsiLB */

/***********************************************************************************************************************************************/

static POLYNOME Z3A_A_PsiR_PsiLB (void)
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;

  POLYNOME p1 = vertex_A_PsiR_PsiRB (b) ;
  POLYNOME p2 = prop_PsiRB_PsiR (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_A_PsiR_PsiRB (a) ;
  POLYNOME p4 = prop_PsiRB_PsiR (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = vertex_A_PsiR_PsiRB (c) ;
  POLYNOME p6 = prop_AL (b,c,0) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = newMultiProduct (ppp) ;
  printf ("Z3 Classic vertex A psi with A under: expect 0 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;

  return pp ;
} /* Z3A_A_PsiR_PsiLB */

/***********************************************************************************************************************************************/

static POLYNOME Z3B_A_PsiR_PsiLB (void)
{
  char a = newDummyIndex () ;
  char b = newDummyIndex () ;
  char c = newDummyIndex () ;
  char d = newDummyIndex () ;
  char e = newDummyIndex () ;

  POLYNOME p1 = vertex_B_PsiR_PsiLB (b,c) ;
  POLYNOME p2 = prop_PsiLB_PsiL (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p3 = vertex_A_PsiL_PsiLB (a) ;
  POLYNOME p4 = prop_PsiLB_PsiL (0) ;   /* (1/(k+p+q)^2 */ 
  POLYNOME p5 = vertex_BB_PsiL_PsiRB (d,e) ;
  POLYNOME p6 = prop_BB_B (d,e,b,c,0) ; /* 1/k^2 */
  POLYNOME ppp[7] = {p1, p2, p3, p4, p5, p6, 0} ;

  POLYNOME pp = newMultiProduct (ppp) ;
  printf ("Z3 Classic vertex A psi with B under: expect 0 \n") ;
  showPol (pp) ;
  pp = dimIntegral (pp) ;
  showPol(pp) ;
  pp = expand (pp) ;
  pp = contractIndices (pp) ;
  showPol(pp) ;

  showPol (pp) ;

  printf ("Z3 TEST \n") ;
  p1 = newSigma ('a') ;
  strcpy (p1->tt.sigma, "bdecbde") ;
  showPol(p1) ;
  p1 = contractIndices (p1) ;
  showPol(p1) ;
  p1 = pauliCleanUp (p1) ;
  showPol (p1) ;
 
  return pp ;
} /* Z3B_A_PsiR_PsiLB */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static BOOL polynomeTest (void)
{
  POLYNOME p1, p2, p3, p4, pp, gac, sga, sgb, sbc , sbd ;
  BOOL debug = FALSE ;
  char alpha ;
  sga = newSigma ('a') ;
  sgb = newSigma ('b') ;
  sbc = newSigB ('c') ;
  sbd = newSigB ('d') ;
  gac = newG ('a','c') ;

  p1 = newSum (sga, sgb) ;
  p2 = newSum (sbc, sbd) ;
  p4 = newProduct (p1,p2) ;
  p3 = newProduct (gac, p4) ;

  if (debug) showPol (p3) ;
  p4 = expand (p3) ;
  if (debug) showPol (p4) ;
  contractIndices (p4) ;
  if (debug) showPol (p4) ;
  pp = contractProducts (p4) ;
  if (debug) showPol (p4) ;
  contractIndices (p4) ;
  if (debug) showPol (p4) ;
  pauliTrace (p4) ;
  if (debug) showPol (p4) ;
  contractIndices (p4) ;

  printf ("Test  g_ac (s_ + s_b) * (s_c + s_d)\n") ;
  if (debug) showPol (p4) ;



  printf ("Test  g_ab p_a p_c\n") ;
  p1 = newG ('a','b') ;
  p2 = newSigB ('a') ;
  p3 = newSigB ('c') ;
  p4 = newProduct (p1,p2) ;
  pp = newProduct (p4,p3) ;
  if (debug) showPol (pp) ;
  pp = contractProducts (pp) ;
  contractIndices (pp) ;
  if (debug) showPol (pp) ;

  p1 = newK ('a') ;
  p1->tt.mm[0][1] = 'b' ;
  p1->tt.mm[0][2] = 'c' ;
  p1->tt.mm[0][3] = 'd' ;
  p1->tt.denom[0] = 4 ;
  showPol(p1) ;
  p4 = copyPolynome(p1) ;
  p2 = dimIntegral (p1) ;
  showPol(p2) ;
  p3 = newProduct (newAG('a','b','e','f',0),p2) ;
  showPol(p3) ;
  p3 = expand(p3) ;
  showPol(p3) ;
  p3 = newProduct (newG('e','f'),p3) ;
  showPol(p3) ;
  p3 = expand(p3) ;
  showPol(p3) ;

  showPol(p4) ;
  p3 = newProduct (newG('c','d'),p4) ;
  showPol(p3) ;
  p2 = dimIntegral (p3) ;
  showPol(p2) ;

  printf ("##### Test des integrales  /(k+p)^2 expect 1\n") ;
  p1 = newScalar (1) ;
  p1->tt.mm[0][0] = newDummyIndex () ;
  p1->tt.denom[0] = 1 ;
  p1->tt.denom[1] = 0 ;
  p1->tt.denom[2] = 1 ;
  showPol(p1) ;
  p2 = dimIntegral (p1) ;
  showPol(p2) ;
  alpha = newDummyIndex () ;
  p3 = momentaCleanUp (p2, alpha) ;
  printf ("##### expect -(p+q)/2 \n") ;
  showPol(p3) ;

  printf ("##### test k3/k^6 cases\n") ;
  p1 = newScalar (1) ;
  p1->tt.mm[0][0] = newDummyIndex () ;
  p1->tt.mm[0][1] = newDummyIndex () ;
  p1->tt.mm[0][2] = newDummyIndex () ;
  p1->tt.denom[0] = 1 ;
  p1->tt.denom[1] = 1 ;
  p1->tt.denom[2] = 1 ;


  showPol(p1) ;

  p2 = dimIntegral (p1) ;
  showPol(p2) ;
  alpha = newDummyIndex () ;
  p3 = momentaCleanUp (p2, alpha) ;
  printf ("##### check please \n") ;
  showPol(p3) ;

  printf ("##### pauli traces\n") ;
  p1 = newSigma ('a') ;
  p1->tt.sigma[0] = 'a' ;
  p1->tt.sigma[1] = 'b' ;

  showPol (p1) ;
  p2 = pauliTrace (p1) ;
  showPol (p2) ;

  p1 = newSigma ('a') ;
  p1->tt.sigma[0] = 'a' ;
  p1->tt.sigma[1] = 'b' ;
  p1->tt.sigma[2] = 'c' ;
  p1->tt.sigma[3] = 'd' ;

  showPol (p1) ;
  p2 = pauliTrace (p1) ;
  showPol (p2) ;

  p1 = newSigma ('a') ;
  p1->tt.sigma[0] = 'a' ;
  p1->tt.sigma[1] = 'b' ;
  p1->tt.sigma[2] = 'c' ;
  p1->tt.sigma[3] = 'd' ;
  p1->tt.sigma[4] = 'e' ;
  p1->tt.sigma[5] = 'f' ;

  showPol (p1) ;
  p2 = pauliTrace (p1) ;
  showPol (p2) ;
  

  printf ("#########################################\n") ;
  firstDummyIndex = 'a' ;
  if (1)
    {
      char a, b, c, d ;
      pp = newScalar (1/2.0) ; /* we derive twice in p_ab */
      pp->tt.mm[0][0] = a = newDummyIndex () ;
      pp->tt.mm[0][1] = b = newDummyIndex () ;
      pp->tt.denom[0] = 1 ;
      pp->tt.denom[1] = 1 ;
      if (0)
	{
	  pp->tt.g[0] = c = newDummyIndex () ;
	  pp->tt.g[1] = d = newDummyIndex () ;
	}
      printf ("# test des integrales  (g_ac k_bd /k^2 (k+p)^2)\n") ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      printf ("# Expect 1/3 p_ab - 1/12 g_ab p_cc\n") ;
      showPol(pp) ;
      printf ("# trace it and expect zero\n") ;
      p1 = newScalar (1) ;
      p1->tt.g[0] = a ;
      p1->tt.g[1] = b ;
      p2 = newProduct (p1, pp) ;
      showPol(p1) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      showPol(p2) ;
      p2 = squareMomentaCleanUp (p2, a) ;
      showPol(p2) ;
    }
  
  
  firstDummyIndex = 'a' ;
  if (0)
    {
      char a,b,c,d ;
      pp = newScalar (1) ;
      pp->tt.mm[0][0] = a = newDummyIndex () ;
      pp->tt.mm[0][1] = b = newDummyIndex () ;
      pp->tt.mm[0][2] = c = newDummyIndex () ;
      pp->tt.mm[0][3] = d = newDummyIndex () ;
      pp->tt.denom[0] = 1 ;
      pp->tt.denom[1] = 1 ;
      printf ("# test des integrales k_abcd / k^2 (k+p)2\n") ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      showPol(pp) ;
      printf ("#  FAUX on voit des indicies libres nouveaux\n") ;
    }
  firstDummyIndex = 'a' ;
  if (1)
    {
      char a,b,c,d ;
      pp = newScalar (1) ;
      pp->tt.mm[0][0] = a = newDummyIndex () ;
      pp->tt.mm[0][1] = b = newDummyIndex () ;
      pp->tt.mm[0][2] = c = newDummyIndex () ;
      pp->tt.mm[0][3] = d = newDummyIndex () ;
      pp->tt.denom[0] = 2 ;
      pp->tt.denom[1] = 2 ;
      printf ("# test des integrales k_abcd / k^4 (k+p)4\n") ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      showPol(pp) ;
      printf ("# trace it, expect 1: CORRECT\n") ;
      p1 = newScalar (1) ;
      p1->tt.g[0] = a ;
      p1->tt.g[1] = b ;
      p1->tt.g[2] = c ;
      p1->tt.g[3] = d ;
      p2 = newProduct (p1, pp) ;
      showPol(p1) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      showPol(p2) ;
    }

  firstDummyIndex = 'a' ;
  if (1)
    {
      char a,b,c,d ;
      pp = newScalar (1) ;
      pp->tt.mm[0][0] = a = newDummyIndex () ;
      pp->tt.mm[0][1] = b = newDummyIndex () ;
      pp->tt.mm[0][2] = c = newDummyIndex () ;
      pp->tt.mm[0][3] = d = newDummyIndex () ;
      pp->tt.denom[0] = 2 ;
      pp->tt.denom[1] = 1 ;
      printf ("# test des integrales k_abcd / k^4 (k+p)2\n") ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      showPol(pp) ;
      printf ("# FAUX je veux du p_ab - g_ab p^2 ou un trc du genre , pas du pur p^2\n") ;
      printf ("# trace it, expect 1: CORRECT\n") ;
      p1 = newScalar (1) ;
      p1->tt.g[0] = a ;
      p1->tt.g[1] = b ;
      p1->tt.g[2] = c ;
      p1->tt.g[3] = d ;
      p2 = newProduct (p1, pp) ;
      showPol(p1) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      showPol(p2) ;
    }

  firstDummyIndex = 'a' ;
  if (1)
    {
      char a,b,c,d,e,f ;
      pp = newScalar (1) ;
      pp->tt.mm[0][0] = a = newDummyIndex () ;
      pp->tt.mm[0][1] = b = newDummyIndex () ;
      pp->tt.mm[0][2] = c = newDummyIndex () ;
      pp->tt.mm[0][3] = d = newDummyIndex () ;
      pp->tt.mm[0][4] = e = newDummyIndex () ;
      pp->tt.mm[0][5] = f = newDummyIndex () ;
      pp->tt.denom[0] = 2 ;
      pp->tt.denom[1] = 2 ;
      pp->tt.denom[2] = 1 ;
      printf ("# test des integrales k_abcdefgh / k^4 (k+p)4 (k+p+q)^4\n") ;
      showPol(pp) ;
      pp = contractIndices (pp) ;
      showPol(pp) ;
      pp = dimIntegral (pp) ;
      showPol(pp) ;
      printf ("# trace it, expect 1: CORRECT\n") ;
      p1 = newScalar (1) ;
      p1->tt.g[0] = a ;
      p1->tt.g[1] = b ;
      p1->tt.g[2] = c ;
      p1->tt.g[3] = d ;
      p1->tt.g[4] = e ;
      p1->tt.g[5] = f ;
      showPol(p1) ;
      p1 = contractIndices (p1) ;
      showPol(p1) ;

      p2 = newProduct (p1, pp) ;

      showPol(p2) ;
      p2 = expand(p2) ;
      showPol(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      showPol(p2) ;
    }

  firstDummyIndex = 'a' ;
  if (1)
    {
      char a,b,c,d ;
      pp = newScalar (1) ;
      pp->tt.eps[0] = a = newDummyIndex () ;
      pp->tt.eps[1] = b = newDummyIndex () ;
      pp->tt.eps[2] = c = newDummyIndex () ;
      pp->tt.eps[3] = d = newDummyIndex () ;
      printf ("# test des espilon\n") ;
      showPol(pp) ;
      p1 = newScalar (1) ;
      p1->tt.g[0] = a ;
      p1->tt.g[1] = b ;
      showPol(p1) ;
      p2 = newProduct (p1, pp) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      printf ("# g_ab eps_abcd  expect zero \n") ;
      showPol(p2) ;

      p3 = newScalar (1) ;
      p3->tt.mm[1][0] = c ;
      p3->tt.mm[1][1] = a ;
      showPol(p3) ;
      p2 = newProduct (p3, pp) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      printf ("# (p)_ab eps_abcd  expect zero \n") ;
      showPol(p2) ;

      p3->tt.mm[1][0] = b ;
      p2 = newProduct (p3, p1) ;
      p2 = expand(p2) ;
      p2 = contractIndices(p2) ;
      p2 = expand(p2) ;
      printf ("# (p)_ab g_ab  expect p_aa \n") ;
      showPol(p2) ;
    }
  return 0 ;
}

/*************************************************************************************************/

static BOOL feynmanDiagrams (void)
{
  BOOL failed = FALSE ;
  POLYNOME p ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 H avec loop Psi, expect je sais pas\n") ;
  if (0) Z2_H__PsiL__Psi () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 A avec loop Psi, expect je sais pas\n") ;
  if (0) Z2_A__PsiL__Psi () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 B avec loop Psi, expect je sais pas\n") ;
  if (0) Z2_B__PsiL__Psi () ;


  firstDummyIndex = 'a' ;
  printf ("============Z2 B avec new vertex expect je sais pas\n") ;
  if (0) Z2_B__AH () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 A avec loop new vertex expect zero\n") ;
  if (0) Z2_A__BH () ;

  if (0) failed = polynomeTest () ;

  firstDummyIndex = 'a' ;
  printf ("============Hodge, check that PP is a projector\n") ;
  if (0) Hodge () ;

  firstDummyIndex = 'a' ;
  printf ("============t'Hooft, check some integrals\n") ;
  if (0) Thooft () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 H avec transient A expect je sais pas\n") ;
  if (0) Z2_H__AH () ;


  firstDummyIndex = 'a' ;
  printf ("============Z2 H avec new vertex ABH\n") ;
  if (0) Z2_H__AB () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 B avec transient A expect je sais pas\n") ;
  if (0) Z2_B__AB () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 B avec new vertex ABH\n") ;
  if (0) Z2_B__AH () ;


  
  firstDummyIndex = 'a' ;
  printf ("============Z2 H avec loop new vertex expect je sais pas\n") ;
  if (0) Z2_H__AB () ;


  firstDummyIndex = 'a' ;
  printf ("============Z2 A avec loop new vertex expectzero\n") ;
  if (0) Z2_A__BH () ;


  firstDummyIndex = 'a' ;
  printf ("============Z2 B avec transient A expect je sais pas\n") ;
  if (0) Z2_B__AH () ;


  firstDummyIndex = 'a' ;
  printf ("============Z2 H avec new vertex expect je sais pas\n") ;
  if (0) Z2_H__AB () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 B avec new vertex expect je sais pas\n") ;
  if (0) Z2_B__AH () ;

  firstDummyIndex = 'a' ;
  printf ("===========Z2 Psi left avec loop Phi, expect  1/2\n") ;
  if (0) Z2_PsiL__H_Psi () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 Psi left avec loop A, expect  -1\n") ;
  if (0) Z2_PsiL__A_Psi () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 Psi left avec loop B, expect -3/2\n") ;
  if (0) Z2_PsiL__B_Psi () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 H avec loop Psi, expect je sais pas\n") ;
  if (0) Z2_H__PsiL__Psi () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 A avec loop Psi, expect je sais pas\n") ;
  if (0) Z2_A__PsiL__Psi () ;

  firstDummyIndex = 'a' ;
  printf ("============Z2 B avec loop Psi, expect je sais pas\n") ;
  if (0) Z2_B__PsiL__Psi () ;

  firstDummyIndex = 'a' ;
  printf ("============Z3 A HB B psi loop, expect je sais pas\n") ;
  if (0) Z3_A_HB_H__psiLoop () ;

  firstDummyIndex = 'a' ;
  printf ("============Z3 A A A psi loop, expect je sais pas\n") ;
  if (1) Z3_A_A_A__psiLoop () ;

    exit (0) ;
  


  firstDummyIndex = 'a' ;
  printf ("============Z3 A BB B psi loop, expect je sais pas\n") ;
  if (1) Z3_A_BB_B__psiLoop () ;

  firstDummyIndex = 'a' ;
  printf ("============Z3 H Psi classic vertex with photon under, expect je hsais pas\n") ;
  if (1) Z3_H_PsiR_PsiLB () ;

  firstDummyIndex = 'a' ;
  printf ("============ZF H Psi new vertex, expect je sais pas\n") ;
  if (1) ZF_H_PsiR_PsiLB () ;

  firstDummyIndex = 'a' ;
  printf ("============Z3 B psi classic vertex with photon under, expect NULL\n") ;
  if (1) Z3_B_PsiR_PsiLB () ;

  firstDummyIndex = 'a' ;
  printf ("============Z3 B psi classic vertex with photon over, expect je sais pas\n") ;
  if (1) Z3_B__AB__PsiR_PsiLB () ;

  firstDummyIndex = 'a' ;
  printf ("============ZF B psi new vertex, expect je sais pas\n") ;
  if (1) ZF_B_PsiR_PsiLB () ;

  firstDummyIndex = 'a' ;
  printf ("============Z3 A psi classic vertex H under, expect je sais pas\n") ;
  if (1) Z3H_A_PsiR_PsiLB () ;

  firstDummyIndex = 'a' ;
  printf ("============Z3 A psi classic vertex A under, expect je sais pas\n") ;
  if (1) Z3A_A_PsiR_PsiLB () ;

  firstDummyIndex = 'a' ;
  printf ("============Z3 A psi classic vertex B under, expect je sais pas\n") ;
  if (1) Z3B_A_PsiR_PsiLB () ;

  p = newPolynome () ;
  printf ("feynmanDiagrams done, %d polynomes created \n", p->id) ;
  return failed ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/*****  SU(2/1) representation theory. This is used by, but does not depend on the analysis above of the Feynman diagrams **********************/
/*****  Casimir studies with Pater Jarvis, mars 2021 **********************************************/
/***** Scalar anomaly paper is below **********************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

typedef struct kasStruct { MX *mu, *QQ, gg, GG, kas2, CHI, ccc, CCC, kas3 ; int a, b, d, d1, d2, d3, d4, chi, scale, COQ ; BOOL isOSp, isSU2 ; AC_HANDLE h ; } KAS ;
typedef struct comtpStruct { int a, b, c, n, s ; } LC ;
static MX KasCommut (MX a, MX b, int sign, KAS *kas) ;
  
static MX *KasimirConstructSU2Matrices (KAS *kas)
{
  MX  muE, muF, muH  ; /* the 5 generators of OSp(2/1) in the Chevalley basis */
  MX *mu ;
  int d, d1 = 0, d2 = 0 ;
  int i ;
  int a ;
  AC_HANDLE h = kas->h ;
  mu = kas->mu = (MX *) halloc (10 * sizeof (MX), kas->h) ;

  kas->chi = 1 ;
  kas->isSU2 = TRUE ;
  
  a = kas->a = kas->a - 2000 ;
  d1 = a + 1 ;
  d2 = 0 ;
  d = d1 + d2 ;

  kas->d = d = d1 + d2 ;
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = 0 ;
  kas->d4 = 0 ;
      
  
  muH = mxCreate (h,  "muH", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF", MX_INT, d, d, 0) ;
  
  int xx[d*d] ;

  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = d1 - 1 - 2*i ;
  mxSet (muH, xx) ;
  mxShow (muH) ;
 
   /* even raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    xx[d * i + i - 1] = i * (a  - i + 1) ;
  
  mxSet (muE, xx) ;
  mxShow (muE) ;
  
  /* even lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d - 1 ; i++)
    {
      if (i != d1 - 1) xx[d * i + i + 1] = 1 ;
    }
  mxSet (muF, xx) ;
  mxShow (muF) ;

  mu[0] =0 ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = 0 ; mu[5] = 0 ; mu[6] = 0 ; mu[7] = 0 ;
  mu[8] = 0 ;
  mu[9] = 0 ;

  return mu ;
} /* KasimirConstructSU2Matrices */

static MX *KasimirConstructOSp1_2Matrices (KAS *kas)
{
  MX  muE, muF, muH, muS, muT  ; /* the 5 generators of OSp(2/1) in the Chevalley basis */
  MX *mu ;
  int d, d1 = 0, d2 = 0 ;
  int i, j ;
  int a ;
  AC_HANDLE h = kas->h ;
  mu = kas->mu = (MX *) halloc (10 * sizeof (MX), kas->h) ;
  kas->chi = 1 ;
  kas->isOSp = TRUE ;
  
  a = kas->a = kas->a - 1000 ;
  d1 = a + 1 ;
  d2 = a ;
  d = d1 + d2 ;

  kas->d = d = d1 + d2 ;
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = 0 ;
  kas->d4 = 0 ;
      
  
  muH = mxCreate (h,  "muH", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF", MX_INT, d, d, 0) ;
  muS = mxCreate (h,  "muS", MX_INT, d, d, 0) ;
  muT = mxCreate (h,  "muT", MX_INT, d, d, 0) ;
  
  int xx[d*d] ;

  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = d1 - 1 - 2*i ;
  for (i = 0 ; i < d2 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = d2 - 1 - 2 * i ;
    }
  mxSet (muH, xx) ;
  mxShow (muH) ;
 
   /* even raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    xx[d * i + i - 1] = i * (a  - i + 1) ;
  for (i = 1 ; i < d2 ; i++)
    xx[d * (d1 + i) + d1 + i - 1] =  i * (a - i ) ;
  
  mxSet (muE, xx) ;
  mxShow (muE) ;
  
  /* even lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d - 1 ; i++)
    {
      if (i != d1 - 1) xx[d * i + i + 1] = 1 ;
    }
  mxSet (muF, xx) ;
  mxShow (muF) ;

 /* odd raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (i+1)  + d1 + i] = (i + 1) ;
      xx[d * (d1 + i)  + i] = d2 - i ;
    }
  mxSet (muS, xx) ;
  mxShow (muS) ;
  
 /* odd lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (i)  + d1 + i] = -1 ;
      xx[d * (d1 + i)  + i + 1] = 1 ;
    }
  mxSet (muT, xx) ;
  mxShow (muT) ;
  
  mu[0] =0 ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = muS ; mu[5] = muT ; mu[6] = 0 ; mu[7] = 0 ;
  mu[8] = 0 ;
  mu[9] = 0 ;

  return mu ;
} /* KasimirConstructOSp1_2Matrices */

/***********************************************************************************************************************************************/

static MX *KasimirConstructAtypicMatrices (KAS *kas)
{
  MX muY, muE, muF, muH, muU, muV, muW, muX, muK1, muK2 ; /* the 8 generators of SU(2/1) in the Chevalley basis */
  MX *mu ;
  int d, d1 = 0, d2 = 0, d3 = 0, d4 = 0 ;
  int i, j ;
  int a = kas->a ;
  AC_HANDLE h = kas->h ;
  
  mu = kas->mu = (MX *) halloc (10 * sizeof (MX), kas->h) ;
  kas->chi = 1 ;
  
   /* atypic 1 */
  d1 = a + 1 ;
  d2 = a ;
  
  d = d1 + d2 ;
  kas->d = d ;
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = 0 ;
  kas->d4 = 0 ;
      
  
  muY = mxCreate (h,  "muY", MX_INT, d, d, 0) ;
  muH = mxCreate (h,  "muH", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF", MX_INT, d, d, 0) ;
  muU = mxCreate (h,  "muU", MX_INT, d, d, 0) ;
  muV = mxCreate (h,  "muV", MX_INT, d, d, 0) ;
  muW = mxCreate (h,  "muW", MX_INT, d, d, 0) ;
  muX = mxCreate (h,  "muX", MX_INT, d, d, 0) ;
  muK1 = mxCreate (h,  "muK1", MX_INT, d, d, 0) ;
  muK2 = mxCreate (h,  "muK2", MX_INT, d, d, 0) ;
  
  int xx[d*d] ;

  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = d1 - 1 - 2*i ;
  for (i = 0 ; i < d2 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = d2 - 1 - 2 * i ;
    }
  mxSet (muH, xx) ;
  mxShow (muH) ;
 
  /* Y hypercharge  Y = diag (a,a...a/a+1,...a+1) */
 
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = -a ;
  for (i = 0 ; i < d2 + d3 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = -a - 1 ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 + i ;
      xx[d * j + j] = -a - 2 ;
    }
  mxSet (muY, xx) ;
  mxShow (muY) ;
  
  /* odd Cartan operator K = diag (0,-1,-2,...-a/-1,-2...-a) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    {
      xx[d * i + i] = -i ;
      j = d1 + i - 1 ;
      xx[d * j + j] = -i ;
    }
  mxSet (muK1, xx) ;
  mxShow (muK1) ;
  
  /* odd Cartan operator K = diag (-a,...-2,-1, 0/--a,...-2,-1) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < a ; i++)
    {
      xx[d * i + i] = -a + i ;
      j = d1 + i ;
      xx[d * j + j] = -a + i ;
    }
  mxSet (muK2, xx) ;
  mxShow (muK2) ;
  
  
   /* even raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    xx[d * i + i - 1] = i * (a  - i + 1) ;
  for (i = 1 ; i < d2 ; i++)
    xx[d * (d1 + i) + d1 + i - 1] =  i * (a - i ) ;
  
  mxSet (muE, xx) ;
  mxShow (muE) ;
  
  /* even lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d - 1 ; i++)
    {
      if (i != d1 - 1) xx[d * i + i + 1] = 1 ;
    }
  mxSet (muF, xx) ;
  mxShow (muF) ;

 /* odd raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (d1+i) + i + 1] = -1 ;
    }
  mxSet (muU, xx) ;
  mxShow (muU) ;
  
 /* odd lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (i + 1)  + d1 + i] = i+1 ;
    }
  mxSet (muV, xx) ;
  mxShow (muV) ;


  /* other raising operator */  
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (d1+i) + i] = -a + i ;
    }
  mxSet (muW, xx) ;
  mxShow (muW) ;

 /* other lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      xx[d * (i)  + d1 + i] = 1 ;
    }
  mxSet (muX, xx) ;
  mxShow (muX) ;
  
  mu[0] = muY ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = muW ; mu[5] = muX ;  mu[6] = muU ; mu[7] = muV ; 
  mu[8] = muK1 ;
  mu[9] = muK2 ;
  
  return mu ;
} /* KasimirConstructAtypicMatrices */

/***********************************************************************************************************************************************/
static MX *KasimirConstructAntiMatrices (KAS *kas)
{
  MX muY, muE, muF, muH, muU, muV, muW, muX, muK1, muK2 ; /* the 8 generators of SU(2/1) in the Chevalley basis */
  MX *mu ;
  int d, d1 = 0, d2 = 0, d3 = 0, d4 = 0 ;
  int i, j, k ;
  int a = kas->a ;
  int b = kas->b ;
  AC_HANDLE h = kas->h ;
  mu = kas->mu = (MX *) halloc (10 * sizeof (MX), kas->h) ;

  kas->chi = -1 ;
  
  if (1)
    {
      d1 = a + 1 ;
      d2 = a + 2 ;
    }

  d = d1 + d2 ;
  kas->d = d ;
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = 0 ;
  kas->d4 = 0 ;
      
  
  muY = mxCreate (h,  "muY", MX_INT, d, d, 0) ;
  muH = mxCreate (h,  "muH", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF", MX_INT, d, d, 0) ;
  muU = mxCreate (h,  "muU", MX_INT, d, d, 0) ;
  muV = mxCreate (h,  "muV", MX_INT, d, d, 0) ;
  muW = mxCreate (h,  "muW", MX_INT, d, d, 0) ;
  muX = mxCreate (h,  "muX", MX_INT, d, d, 0) ;
  muK1 = mxCreate (h,  "muK1", MX_INT, d, d, 0) ;
  muK2 = mxCreate (h,  "muK2", MX_INT, d, d, 0) ;
  
  int xx[d*d] ;
 
  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0, k = a ; i < d2 ; k -= 2, i++)
    {
      if ( i < d1)
	xx[d * i + i] = k ;
      j = d1 + i ;
      xx[d * j + j] = k + 1 ;
    }
  mxSet (muH, xx) ; 
  mxShow (muH) ;


  /* Y hypercharge  Y = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = 2*b -a ;
  for (i = 0 ; i < d2 + d3 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = 2*b -a - 1 ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 + i ;
      xx[d * j + j] = 2*b - a - 2 ;
    }
  mxSet (muY, xx) ;
  mxShow (muY) ;

  /* odd Cartan operator K = diag (a,...2,1,/ a,...2,1,0) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      j = d1 + i  ;
      xx[d * j + j] = d1 - i ;
      if (i < d1)
	xx[d * i + i] = d1 - i ;
    }
  mxSet (muK1, xx) ;
  mxShow (muK1) ;

  /* odd Cartan operator K = diag (1,2,...a/0,1,2...a) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d2 ; i++)
    {
      j = d1 + i  ;
      xx[d * j + j] = i ;
      if (i < d1)
	xx[d * i + i] = i + 1;
    }
  mxSet (muK2, xx) ;
  mxShow (muK2) ;
  
  /* even raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    xx[d * i + i - 1] = i * (a + 1 - i ) ;
  for (i = 1 ; i < d2 ; i++)
    xx[d * (d1 + i) + d1 + i - 1] =  i * (a - i + 2 ) ;
  mxSet (muE, xx) ;
  mxShow (muE) ;

  /* even lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d - 1 ; i++)
    {
      if (i != d1 - 1) xx[d * i + i + 1] = 1 ;
    }
  mxSet (muF, xx) ;
  mxShow (muF) ;

  /* odd raising operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {
      xx[d * (d1+i) + i] = 1 ;
    }
  mxSet (muU, xx) ;
  mxShow (muU) ;
  
/* odd lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {
      xx[d * (i )  + d1 + i] = (d1 - i) ;
    }
  mxSet (muV, xx) ;
  mxShow (muV) ;

  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {
      xx[d * (d1+i+1) + i] = -1 - i ;
    }
  mxSet (muW, xx) ;
  mxShow (muW) ;

  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {
      xx[d * (i)  + d1 + i + 1] = -1 ;
    }
  mxSet (muX, xx) ;
  mxShow (muX) ;
  
  mu[0] = muY ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = muW ; mu[5] = muX ;  mu[6] = muU ; mu[7] = muV ; 
  mu[8] = muK1 ;
  mu[9] = muK2 ;

  return mu ;
} /* KasimirConstructAntiMatrices */

/***********************************************************************************************************************************************/

static MX *KasimirConstructTypicMatrices (KAS *kas)
{
  MX muY, muE, muF, muH, muU, muV, muW, muX, muK1, muK2 ; /* the 8 generators of SU(2/1) in the Chevalley basis */
  MX *mu ;
  int d, d1 = 0, d2 = 0, d3 = 0, d4 = 0 ;
  int i, j ;
  int a = kas->a, b = kas->b ;
  int s = 1 ;  /* scaling U V K1 K2 */
  AC_HANDLE h = kas->h ;
  mu = kas->mu = (MX *) halloc (10 * sizeof (MX), kas->h) ;

  kas->chi = -1 ;
  s = a + 1 ; 
  kas->scale = s * s ;
  
  if (1)
    {
      d1 = d4 = a + 1 ;
      d2 = a + 2 ;
      d3 = a  ;
    }

  d = d1 + d2 + d3 + d4 ;
  int xx[d*d] ;
  const int *xx1 = messalloc (d*d*sizeof(int)) ;
  const int *xx2 = messalloc (d*d*sizeof(int)) ;
      
  kas->d = d ;
  kas->d1 = d1 ;
  kas->d2 = d2 ;
  kas->d3 = d3 ;
  kas->d4 = d4 ;
  
  
  muY = mxCreate (h,  "muY: Y Hypercharge", MX_INT, d, d, 0) ;
  muH = mxCreate (h,  "muH: Even SU(2) Cartan operator", MX_INT, d, d, 0) ;
  muE = mxCreate (h,  "muE: Even raising operator", MX_INT, d, d, 0) ;
  muF = mxCreate (h,  "muF: Even lowering operator", MX_INT, d, d, 0) ;
  muU = mxCreate (h,  "muU: Odd raising operator", MX_INT, d, d, 0) ;
  muV = mxCreate (h,  "muV: Odd lowering operator", MX_INT, d, d, 0) ;
  muW = mxCreate (h,  "muW: Other odd raising operator", MX_INT, d, d, 0) ;
  muX = mxCreate (h,  "muX: Other odd lowering operator", MX_INT, d, d, 0) ;
  muK1 = mxCreate (h,  "muK1: K1 = {U,V}", MX_INT, d, d, 0) ;
  muK2 = mxCreate (h,  "muK2: K2 = {W,X}", MX_INT, d, d, 0) ;
  
  /* even Cartan operator H = diag (a, a-2 .... -a in each SU(2) sector */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = d1 - 1 - 2*i ;
  for (i = 0 ; i < d2 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = d2 - 1 - 2 * i ;
    }
  for (i = 0 ; i < d3 + d3 ; i++)
    {
      j = d1 + d2 + i ;
      xx[d * j + j] = d3 - 1 - 2 * i ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 + i ;
      xx[d * j + j] = d4 - 1 - 2 * i ;
    }
  mxSet (muH, xx) ;
  mxShow (muH) ;

  /* Y hypercharge  Y = diag (a+1,a+1...a+1/a,a,....a) */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    xx[d * i + i] = 2*b - a ;
  for (i = 0 ; i < d2 + d3 ; i++)
    {
      j = d1 + i ;
      xx[d * j + j] = 2*b -a - 1 ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 + i ;
      xx[d * j + j] = 2*b -a - 2 ;
    }
  mxSet (muY, xx) ;
  mxShow (muY) ;

  
  /* odd Cartan operator K1 = diag (a,...2,1,/ a,...2,1,0) */
  /* odd Cartan operator K2 = diag (1,2,...a/0,1,2...a) */
  memset (xx, 0, sizeof (xx)) ;
  mxValues (muY, &xx1, 0, 0) ;
  mxValues (muH, &xx2, 0, 0) ;

    
  for (i = 0 ; i < d * d ; i++)
    xx[i] =  (xx1[i] + xx2[i])/2 ;
  mxSet (muK1, xx) ;
  mxShow (muK1) ;
  
  for (i = 0 ; i < d * d ; i++)
    xx[i] =  (xx1[i] - xx2[i])/2 ;
  mxSet (muK2, xx) ;
  mxShow (muK2) ;
  

  /* even raising operator E */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    {
      j = 0 ;
      xx[d * (j + i) + j + i - 1] = i * (d1 - i) ;
    }
  for (i = 1 ; i < d2 ; i++)
    {
      j = d1 ;
      xx[d * (j + i) + j + i - 1] = i * (d2 - i) ;
    }
  for (i = 1 ; i < d3  ; i++)
    {
      j = d1 + d2 ;
      xx[d * (j + i) + j + i - 1] = i * (d3 - i) ;
    }
  for (i = 1 ; i < d4  ; i++)
    {
      j = d1 + d2 + d3 ;
      xx[d * (j + i) + j + i - 1] = i * (d4 - i) ;
    }
  mxSet (muE, xx) ;
  mxShow (muE) ;

  /* even lowering operator F */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 1 ; i < d1 ; i++)
    {
      j = 0 ;
      xx[d * (j + i - 1) + j + i] = 1 ;
    }
  for (i = 1 ; i < d2 ; i++)
    {
      j = d1 ;
      xx[d * (j + i - 1) + j + i] = 1 ;
    }
  for (i = 1 ; i < d3  ; i++)
    {
      j = d1 + d2 ;
      xx[d * (j + i - 1) + j + i] = 1 ;
    }
  for (i = 1 ; i < d4  ; i++)
    {
      j = d1 + d2 + d3 + 0 ;
      xx[d * (j + i - 1) + j + i] = 1 ;
    }
  mxSet (muF, xx) ;
  mxShow (muF) ;
  
  /* odd raising operator */  
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {

      j = d1 ;
      xx[d * (j + i) + i] = s*b ;
    }
  for (i = 0 ; i < d3 ; i++)
    {
      j = d1 + d2 ;
      xx[d * (j + i ) + i + 1 ] = s * ( a+1 -b) ;
      j = d1 + d2 + d3 ;
      xx[d * (j+i) + i+d1+d2] = b ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 + d2 + d3 ;
      xx[d * (j + i ) + d1 + i + 1 ] = b - s ;
    }
  mxSet (muU, xx) ;
  mxShow (muU) ;

  /* odd lowering operator */
  memset (xx, 0, sizeof (xx)) ;
  for (i = 0 ; i < d1 ; i++)
    {
      j = 0 ;
      xx[d * (j+i) + i+d1] = s - i ;
    }
  for (i = 0 ; i < d4 ; i++)
    {
      j = d1 ;
      xx[d * (j+i+1) + i+d1+d2+d3] = s * (i+1) ;
    }
  for (i = 0 ; i < d3 ; i++)
    {
      j = 0 ;
      xx[d * (j+i+1) + i+d1+d2] = - (i+1) ;
      j = d1 + d2 ;
      xx[d * (j+i) + i+d1+d2+d3] = s * ( a - i) ;
    }
  mxSet (muV, xx) ;
  mxShow (muV) ;
  if (0) exit (0) ;

  /* odd other raising operator */
  muW = KasCommut (muE, muU, -1, kas) ;
  muW->name = "muW" ;
  mxShow (muW) ;
  
  /* odd other oweringing operator */
  muX = KasCommut (muV, muF, -1, kas) ;
  mxShow (muX) ;


  mu[0] = muY ; mu[1] = muE ; mu[2] = muF ; mu[3] = muH ;
  mu[4] = muW ; mu[5] = muX ;  mu[6] = muU ; mu[7] = muV ; 
  mu[8] = muK1 ;
  mu[9] = muK2 ;

  return mu ;
} /* KasimirConstructTypicMatrices */

/***********************************************************************************************************************************************/

static void KasimirCheckSuperTrace (KAS *kas)
{
  int i, j, ii ;
  int d = kas->d ;
  int d1 = kas->d1 ;
  int d2 = kas->d2 ;
  int d3 = kas->d3 ;
  int d4 = kas->d4 ;

  for (ii = 0 ; ii < 10 ; ii++)
    {
      int n = 0 ;
      MX mu = kas->mu[ii] ;
      const int *xx ;

      if (! mu)
	continue ;
      mxValues (mu, &xx, 0, 0) ;
      for (i = 0 ; i < d1 ; i++)
	n +=  xx[d*i + i] ;
      for (i = 0 ; i < d2 ; i++)
	{
	  j = d1 + i ;
	  n -=  xx[d*j + j] ;
	}
      for (i = 0 ; i < d3 ; i++)
	{
	  j = d1 + d2 + i ;
	  n -=  xx[d*j + j] ;
	}
      for (i = 0 ; i < d4 ; i++)
	{
	  j = d1 + d2 + d3 + i ;
	  n +=  xx[d*j + j] ;
	}
      n *= kas->chi ;
      printf ("Str(%s) = %d\n", mu->name, n) ;
    }
  return ;
} /* KasimirCheckSuperTrace */

/***********************************************************************************************************************************************/

static MX KasCommut (MX a, MX b, int sign, KAS *kas)
{
  MX p = mxMatMult (a, b, kas->h) ;
  MX q = mxMatMult (b, a, kas->h) ;
  MX r = mxCreate (kas->h, "r", MX_INT, kas->d, kas->d, 0) ;

  r = sign == 1 ? mxAdd (r, p, q, kas->h) : mxSubstract (r, p, q, kas->h) ;
  
  return r ;
}

/***********************************************************************************************************************************************/

static MX KasCheck (LC *up, KAS *kas)
{
  int d = kas->d ;
  int dd = kas->d * kas->d ;
  MX a = kas->mu[up->a] ;
  MX b = kas->mu[up->b] ;
  MX c = kas->mu[up->c] ;
  MX r = mxCreate (kas->h, "r", MX_INT, d,d,0) ;
  MX s = mxCreate (kas->h, "s", MX_INT, d,d,0) ;
  MX t = mxCreate (kas->h, "t", MX_INT, d,d,0) ;
  const int *xx ;
  int yy [dd] ;
  int i, k ;

  MX ab = KasCommut (a, b, up->s, kas) ;
  int scale = (up->s == 1 ? (kas->scale ? kas->scale : 1) : 1) ;

  mxValues (c, &xx, 0, 0) ;
  for (i = 0 ; i < dd ; i++)
    yy[i] = up->n * scale * xx[i] ;
  mxSet (s, yy) ;
  t = mxSubstract (t, ab, s, kas->h) ;
  mxValues (t, &xx, 0, 0) ;
  for (i = k = 0 ; i < dd ; i++)
    k += xx[i] * xx[i] ;
  if (k > 0)
    {
      mxShow (ab) ;
      mxShow (c) ;
      messcrash ("\nKasChect Failed [%s,%s] = %d * %s\n", a->name, b->name, up->n, c->name) ;
    }
  else
    {
      if (up->s == -1 && up->n == 0)
	printf ("[%s,%s] = 0\n",  a->name, b->name) ;
      else if (up->s == 1 && up->n == 0)
	printf ("{%s,%s} = 0\n",  a->name, b->name) ;
      else if (up->s == -1 && up->n == 1)
	printf ("[%s,%s] = %s\n",  a->name, b->name, c->name) ;
      else if (up->s == -1 && up->n == 1)
	printf ("[%s,%s] = - %s\n",  a->name, b->name, c->name) ;
      else if (up->s == -1)
	printf ("[%s,%s] = %d %s\n",  a->name, b->name, up->n, c->name) ;

      else if (up->s == 1 && up->n == 1)
	printf ("{%s,%s} = %s\n",  a->name, b->name, c->name) ;
      else if (up->s == 1 && up->n == 1)
	printf ("{%s,%s} = - %s\n",  a->name, b->name, c->name) ;
      else if (up->s == 1)
	printf ("{%s,%s} = %d %s\n",  a->name, b->name, up->n, c->name) ;
    }
  return r ;
}

/***********************************************************************************************************************************************/

static void KasimirCheckCommutators (KAS *kas)
{
  LC *up, *XXX ;
  LC Su2XXX[] = {
    {1,2,3,1,-1},
    {3,1,1,2,-1},
    {3,2,2,-2,-1},
    {0,0,0,0,0}
  } ;

  LC Su21XXX[] = {
    {1,2,3,1,-1},
    {3,1,1,2,-1},
    {3,2,2,-2,-1},

    {0,1,1,0,-1},
    {0,2,2,0,-1},
    {0,3,3,0,-1},
    
    {0,4,4,1,-1},
    {0,5,5,-1,-1},
    {0,6,6,1,-1},
    {0,7,7,-1,-1},
    
    {3,4,4,1,-1},
    {3,5,5,-1,-1},
    {3,6,6,-1,-1},
    {3,7,7,1,-1},

    {4,4,9,0,1},
    {5,5,9,0,1},
    {6,6,8,0,1},
    {7,7,8,0,1},
    
    {6,7,8,1,1},
    {4,5,9,1,1},
    {6,5,2,-1,1},
    {4,7,1,-1,1},
    {6,4,1,0,1},
    {5,7,2,0,1},
    
    {1,4,4,0,-1},
    {1,6,4,1,-1},
    {2,4,6,1,-1},
    {2,6,6,0,-1},
    {2,7,5,-1,-1},
    {2,5,5,0,-1},
    {1,5,7,-1,-1},
    {1,7,7,0,-1},
    {0,0,0,0,0}
  } ;

  LC OSp21XXX[] = {
    {1,2,3,1,-1},
    {3,1,1,2,-1},
    {3,2,2,-2,-1},

    {3,4,4,1,-1},
    {3,5,5,-1,-1},
    {1,4,4,0,-1},
    {1,5,4,1,-1},
    {2,4,5,1,-1},
    {2,5,5,0,-1},

    {4,4,1,2,1},
    {5,5,2,-2,1},
    {4,5,3,-1,1},

    {0,0,0,0,0}
  } ;

  XXX = kas->isSU2 ? Su2XXX : (kas->isOSp ? OSp21XXX : Su21XXX) ;
  
  for (up = XXX ; up->s ; up++)
    KasCheck (up, kas) ;
  printf ("SUCCESS (a=%d, 0) all comutators have been verified\n", kas->a) ;
  return ;
} /* KasimirCheckCommutators */

/***********************************************************************************************************************************************/

static void  KasimirLowerMetric (KAS *kas)
{
  MX gg ;
  int i, j, k, k1 ;
  float  n, yy[100] ;
  AC_HANDLE h = ac_new_handle () ;
  
  gg = kas->gg = mxCreate (kas->h,  "gg", MX_FLOAT, 10, 10, 0) ;

  printf ("Metric gg:: ") ;
  memset (yy, 0, sizeof (yy)) ;
  for (i = 0 ; i < 8 ; i++)
    for (j = 0 ; j < 8 ; j++)
      {
	int d = kas->d ;
	int d1 = kas->d1 ;
	int d2 = kas->d2 ;
	int d3 = kas->d3 ;
	int s = kas->scale ;
	MX a = kas->mu[i] ;
	MX b = kas->mu[j] ;
	int COQ = kas->COQ ;
	
	if (!a || !b) continue ;
	MX c = mxCreate (h, "c", MX_INT, d, d, 0) ;
	const int *xx ;
	
	c = mxMatMult (a, b, h) ;
	mxValues (c, &xx, 0, 0) ;
	n = 0 ;
	for (k = 0 ; k < d ; k++)
	  {
	    k1 = COQ ? k % (d/COQ) : k ;
	    n += (k1 < d1 || k1 >= d1 + d2 + d3 ? xx[d*k + k] : - xx[d*k + k]) ;
	  }
	n *= kas->chi ;
	if (s > 1 && i>=4 && j>=4)
	  n /= s ;
	if (COQ)
	  n /= COQ ;
	yy [10*i + j] = n/2.0 ;
	if (n != 0)
	  printf (" %d:%d=%.2f ",i,j,n/2.0) ;
      }
  printf ("\n") ;

  mxSet (gg, yy) ;
  ac_free (h) ;
  return  ;
} /* KasimirMetric */

/***********************************************************************************************************************************************/

static void KasimirUpperMetric (KAS *kas)
{
  MX gg, GG ;
  int i, j ;
  const float *xx ;
  float yy[100] ;

  gg = kas->gg ;
  GG = kas->GG = mxCreate (kas->h,  "gg", MX_FLOAT, 10, 10, 0) ;

  mxValues (gg,0, &xx, 0) ;
  memcpy (yy, xx, sizeof (yy)) ;
 
  printf ("Metric GG:: ") ;
  for (i = 0 ; i < 8 ; i++)
    for (j = 0 ; j < 8 ; j++)
      {
	float z = yy[10*i + j] ;
	if (i>=4 && i < 8)
	  z *= -1 ;
	yy[10*i + j] = z ? 1/z : 0 ;
	if (z)
	  printf (" %d:%d=%.2f",i,j,1/z) ;
      }
  mxSet (GG, yy) ;
  
  printf ("\n") ;
  return ;
} /* KasimirMetric */

/***********************************************************************************************************************************************/

static void KasimirOperatorK2 (KAS *kas)
{
  int i, j, k ;
  int d = kas->d ;
  int a = kas->a ;
  int b = kas->b ;
  AC_HANDLE h = ac_new_handle () ;
  const float *xx ;
  const int *yy ;
  float zz [d*d], dz,dz1 ;
  int s = kas->scale ;
  
  memset (zz, 0, sizeof (zz)) ;
  mxValues (kas->GG, 0, &xx, 0) ;
  for (i = 0 ; i < 8 ; i++)
    for (j = 0 ; j < 8 ; j++)
      if (xx[10*i + j])
	{
	  MX a = kas->mu[i] ;
	  MX b = kas->mu[j] ;

	  if (!a || !b)
	    continue ;
	  MX c = mxMatMult (a, b, h) ;
	  float z = xx[10*i + j] ;

	  mxValues (c, &yy, 0, 0) ;
	  if (s > 1 && i>= 4 && j >= 4)
	    z /= s ;
	  for (k = 0 ; k < d*d ; k++)
	    zz[k] += z * yy[k] ;
	}

  MX kas2 = kas->kas2 = mxCreate (kas->h,  "KAS2", MX_FLOAT, d, d, 0) ;
  dz = -4*b * (b - a - 1)/(a+1.0) ;
  dz1 = 4*(2*b -a - 1)*(2*b - 1)/(a+1.0) ;
  if (dz != 0)
    for (i = 0 ; i < d*d ; i++)
      zz[i] /= dz ;
  mxSet (kas2, zz) ;

  printf ("Quadratic super-Casimir operator KAS2 = -4b(b-a-1)/(a+1) * %f 4(2b-a-1)(2b-1)/(a+1) * %f\n", zz[0],zz[d*d/2]*dz/dz1)  ;
  if (a<3) niceShow (kas2) ;

  ac_free (h) ;
  return ;
} /* KasimirOperatorK2 */

/***********************************************************************************************************************************************/

static void GhostKasimirOperatorK2 (KAS *kas)
{
  int i, j, k, l, m1 ;
  int d = kas->d ;
  AC_HANDLE h = ac_new_handle () ;
  const float *xx ;
  const int *yy ;
  float zz [d*d], dz, dz1 ;

  memset (zz, 0, sizeof (zz)) ;
  mxValues (kas->GG, 0, &xx, 0) ;
  for (i = 4 ; i < 8 ; i++)
    for (j = 4 ; j < 8 ; j++)
      if (xx[10*i + j])
	{
	  MX a = kas->mu[i] ;
	  MX b = kas->mu[j] ;

	  if (!a || !b)
	    continue ;
	  MX c = mxMatMult (a, b, h) ;
	  float z = xx[10*i + j] ;
	  if (z>0)z=1;else z=-1;
	  if (kas->scale) z /= kas->scale ;
	  

	  mxValues (c, &yy, 0, 0) ;
	  for (k = 0 ; k < d*d ; k++)
	    zz[k] -= z * yy[k] ;
	}
  for (i = 4 ; i < 8 ; i++)
    for (j = 4 ; j < 8 ; j++)
      for (k = 4 ; k < 8 ; k++)
	for (l = 4 ; l < 8 ; l++)
	  {
	    int jj = (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l);
	    if (jj)
	      {
		MX a = kas->mu[i] ;
		MX b = kas->mu[j] ;
		MX c = kas->mu[k] ;
		MX d1 = kas->mu[l] ;
		float z = 1 ;
		
		if (!a || !b || !c || !d)
		  continue ;
		MX e = mxMatMult (a, b, h) ;
		MX f = mxMatMult (e,c, h) ;
		MX g = mxMatMult (f,d1, h) ;

		if (kas->scale) z /= (kas->scale * kas->scale) ;
		mxValues (g, &yy, 0, 0) ;
		for (m1 = 0 ; m1 < d*d ; m1++)
		  zz[m1] += z * (jj>0  ? yy[m1] : -yy[m1]) ;
	      }
	  }

  MX kas2 = kas->CHI = mxCreate (kas->h,  "Ghost-Casimir2", MX_FLOAT, d, d, 0) ;
  int a = kas->a, b = kas->b ;
  dz = 6*b * (b - a - 1) ;
  dz1 = -6*(2*b -a - 1)*(2*b - 1) ;
  if (dz != 0)
    for (i = 0 ; i < d*d ; i++)
      zz[i] /= dz ;
  mxSet (kas2, zz) ;
  

  printf ("Quadratic Ghost-Casimir operator KAS2 = 6 * b * (b - a - 1) * %.3f, -6(2b-a-1)(2b-1)%f\n", zz[0],zz[d*d/2]*dz/dz1) ;
  if (kas->a<3) niceShow (kas2) ;

  ac_free (h) ;
  return ;
} /* GhostKasimirOperatorK2 */

/***********************************************************************************************************************************************/
/* Casimir proposed by Peter, july 28 */
static void KasimirOperatorK4 (KAS *kas)
{
  int d = kas->d ;
  AC_HANDLE h = ac_new_handle () ;
  MX U, V, W, X ;
  int ii ;
  
  for (ii = 0 ; ii < 2 ; ii++)
    {
      if (ii==0)
	{
	  U = kas->mu[6] ;
	  V = kas->mu[7] ;
	  W = kas->mu[4] ;
	  X = kas->mu[5] ;
	}
      else
	{
	  V = kas->mu[6] ;
	  U = kas->mu[7] ;
	  X = kas->mu[4] ;
	  W = kas->mu[5] ;
	}
      
      MX Y = kas->mu[0] ;
      
      MX WX = mxMatMult (W,X,h) ;
      MX UWX = mxMatMult (U,WX,h) ;
      MX WXU = mxMatMult (WX,U,h) ;
      MX uwx =  mxCreate (h, "[U,WX]", MX_INT, d, d, 0) ;
      uwx = mxSubstract (uwx, UWX, WXU, h) ;
      MX Vuwx = mxMatMult (V, uwx,h) ;
      MX uwxV = mxMatMult (uwx,V,h) ;
      MX vuwx = mxCreate (h, "{V,[U,WX]}", MX_INT, d, d, 0) ;
      vuwx = mxAdd (vuwx, Vuwx,uwxV, h) ;
      
      MX Y2 = mxMatMult (Y,Y,h) ;
      MX Y2WX = mxMatMult (Y2,WX,h) ;
      MX UY2WX = mxMatMult (U,Y2WX,h) ;
      MX Y2WXU = mxMatMult (Y2WX,U,h) ;
      MX uy2wx =  mxCreate (h, "[U,Y2WX]", MX_INT, d, d, 0) ;
      uy2wx = mxSubstract (uy2wx, UY2WX, Y2WXU, h) ;
      MX Vuy2wx = mxMatMult (V, uy2wx,h) ;
      MX uy2wxV = mxMatMult (uy2wx,V,h) ;
      MX vuy2wx = mxCreate (h, "{V,[U,Y2WX]}", MX_INT, d, d, 0) ;
      vuy2wx = mxAdd (vuy2wx, Vuy2wx,uy2wxV, h) ;
        
      niceIntShow (vuwx) ;
      niceIntShow (vuy2wx) ;
    }
  
  ac_free (h) ;
  return ;
} /* KasimirOperatorK4 */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static void  KasimirLower3tensor (KAS *kas)
{
  int i, j, k, i1, scale ;
  float yy[1000] ;
  float zz ;
  AC_HANDLE h = ac_new_handle () ;
  MX ccc = kas->ccc = mxCreate (kas->h,  "ccc", MX_FLOAT, 10, 10, 10, 0) ;
  int mx0 = 0 ;
  int mx1 = 8 ;

  printf ("Lower ccc:: ") ;
  memset (yy, 0, sizeof (yy)) ;
  for (i = mx0 ; i < mx1 ; i++)
    for (j = mx0 ; j < mx1 ; j++)
      for (k = mx0 ; k < mx1 ; k++)
      {
	int d = kas->d ;
	int d1 = kas->d1 ;
	int d2 = kas->d2 ;
	int d3 = kas->d3 ;
	int s ;
	MX a = kas->mu[i] ;
	MX b = kas->mu[j] ;
	MX c = kas->mu[k] ;
	MX u,v,x,y ;
	MX z = mxCreate (h, "z", MX_INT, d, d, 0) ;
	const int *xx ;
	
	u = mxMatMult (a, b, h) ;
	v = mxMatMult (u, c, h) ;
	x = mxMatMult (a, c, h) ;
	y = mxMatMult (x, b, h) ;

	if (j >= 4 && j <= 7 && k >= 4 && k <= 7)
	  s = -1 ;
	else
	  s = 1 ;
	if (i > 40)
	  s = -s ;
	if (s == -1)
	  z = mxSubstract (z, v, y, h) ;
	else
	  z = mxAdd (z, v, y, h) ;
	mxValues (v, &xx, 0, 0) ;
	zz = 0 ;
	for (i1 = 0 ; i1 < d ; i1++)
	  {
	    int COQ = kas->COQ ;
	    int dd2 = COQ ? d/COQ : d ;
	    int i2 = i1 % dd2 ; 
	    zz += (i2 < d1 || i2 >= d1 + d2 + d3 ? xx[d*i1 + i1] : - xx[d*i1 + i1]) ;
	  }
	zz *= kas->chi ;

	scale = (i>=4 || j >= 4 || k >= 4) ? kas->scale : 0 ;
	if (scale != 0)
	  zz /= scale*scale ;
	
	yy [100*i + 10*j + k] = zz ;
	if (zz != 0)
	  printf (" %d:%d:%d=%f",i,j,k,zz) ;
      }
  mxSet (ccc, yy) ;
  ac_free (h) ;
  return  ;
} /* KasimirLower3tensor */

/***********************************************************************************************************************************************/

static void  KasimirUpper3tensor (KAS *kas)
{
  int i, j, k ;
  float yy[1000] ;
  AC_HANDLE h = ac_new_handle () ;
  MX CCC = kas->CCC = mxCreate (kas->h,  "ccc", MX_FLOAT, 10, 10, 10, 0) ;
  const float *GG ;
  const float *ccc ;

  printf ("Upper CCC:: ") ;
  mxValues (kas->GG, 0,  &GG, 0) ;
  mxValues (kas->ccc, 0, &ccc, 0) ;
  memset (yy, 0, sizeof (yy)) ;
  for (i = 0 ; i < 8 ; i++)
    for (j = 0 ; j < 8 ; j++)
      for (k = 0 ; k < 8 ; k++)
	{
	  int a, b, c ; /* dummy indices */
	  float  z = 0 ;
	  for (a = 0 ; a < 8 ; a++)
	    for (b = 0 ; b < 8 ; b++)
	      for (c = 0 ; c < 8 ; c++)
		z += GG[10*i + a] * GG[10*j + b] * GG[10*k + c] * ccc[100*a + 10 * b + c] ; 

	  if (i>4 || j>4 || k>4) z = -z ;
	  if (z != 0)
	    printf (" %d:%d:%d=%.2f",i,j,k,z) ;
	  yy[100*i + 10*j + k] += z ;
	}
  mxSet (CCC, yy) ;
  ac_free (h) ;
  return  ;
} /* KasimirUpper3tensor */

/***********************************************************************************************************************************************/

static void KasimirOperatorK3 (KAS *kas)
{
  int i, j, k, m ;
  int d = kas->d ;
  AC_HANDLE h = ac_new_handle () ;
  const float *xx ;
  const int *yy ;
  float z, zz [d*d] ;
  int mx0 = 0 ;
  int mx1 = 8 ;
  
  memset (zz, 0, sizeof (zz)) ;
  mxValues (kas->CCC, 0, &xx, 0) ;
  for (i = mx0 ; i < mx1 ; i++)
    for (j = mx0 ; j < mx1 ; j++)
      for (k = mx0 ; k < mx1 ; k++)
	if (xx[100*i + 10*j + k])
	  {
	    MX a = kas->mu[i] ;
	    MX b = kas->mu[j] ;
	    MX c = kas->mu[k] ;
	    MX u = mxMatMult (a, b, h) ;
	    MX v = mxMatMult (u, c, h) ;
	    float n = xx[100*i + 10*j + k] ;
	    
	    mxValues (v, &yy, 0, 0) ;
	    for (m = 0 ; m < d*d ; m++)
	      zz[m] += n * yy[m] ;
	  }

  MX kas3 = kas->kas3 = mxCreate (kas->h,  "KAS3", MX_FLOAT, d, d, 0) ;
  mxSet (kas3, zz) ;

  z = kas->b * (kas->b - kas->a - 1) *(kas->a) ;
  if (z == 0) z = 1 ;
  printf ("\nCubic super-Casimir operator KAS3 b(b-a-1)32/9(a+1) * %f\n", 9*zz[0]/(32*z)) ;
  if (kas->a<6) niceShow (kas3) ;

  ac_free (h) ;
  return ;
} /* KasimirOperatorK3 */

/***********************************************************************************************************************************************/

static void Kasimirs (int a, int b)
{
  KAS kas ;
  memset (&kas, 0,sizeof(KAS)) ;
  kas.h = ac_new_handle () ;
  kas.a = a ;    /* Kac Dynkin weights of the heighest weight */
  kas.b = b ;
  kas.isOSp = FALSE ;
  AC_HANDLE h = kas.h ;
  
  if (a>=2000)
    KasimirConstructSU2Matrices (&kas) ;
  else if (a>=1000)
    KasimirConstructOSp1_2Matrices (&kas) ;
  else if (a>0 && b == 0)
    KasimirConstructAtypicMatrices (&kas) ;
  else if (a >= 0 && b == a + 1)
    KasimirConstructAntiMatrices (&kas) ;
  else
    KasimirConstructTypicMatrices (&kas) ;

  KasimirCheckSuperTrace (&kas) ;

  KasimirCheckCommutators (&kas) ;

  KasimirLowerMetric (&kas) ;
  KasimirUpperMetric (&kas) ;
  exit (0) ;
  
  KasimirOperatorK2 (&kas) ;
  GhostKasimirOperatorK2 (&kas) ;
  if (0) KasimirOperatorK4 (&kas) ;


  MX qmuX = kas.mu[6] ;
  int d = kas.d ;
  	printf ("Verify that the casimir commutes with X and Y\n") ;
	MX CKX = mxMatMult (kas.kas2, qmuX, h) ;
	MX CXK = mxMatMult (qmuX, kas.kas2, h) ;
	MX Com =  mxCreate (h, "[casimir,X]", MX_COMPLEX,d,d, 0) ;
	Com = mxSubstract (Com, CKX, CXK, h) ;
	niceShow (Com) ;
	
	printf ("Verify that the S-casimir anticommutes with X and Y\n") ;
	MX SCKX = mxMatMult (kas.CHI, qmuX, h) ;
	MX SCXK = mxMatMult (qmuX, kas.CHI, h) ;
	MX SCom =  mxCreate (h, "{S-casimir,X}", MX_COMPLEX, d,d, 0) ;
	SCom = mxAdd (SCom, SCKX, SCXK, h) ;
	niceShow (SCom) ;
	
	printf ("Compute the square of the S-casimir\n") ;
	MX SC2 = mxMatMult (kas.CHI,kas.CHI, h) ;
	if (0) SC2 = mxLinearCombine (SC2, 1, SC2, -1, kas.kas2, h) ;
	SC2->name = "S-Casimir square" ;
	niceShow (SC2) ;
	
	printf ("Compute the product of the Casimir by the S-casimir Q^3 = Q\n") ;
	MX SC3 = mxMatMult (kas.kas2, kas.CHI, h) ;
	if (0) SC3 = mxLinearCombine (SC3, 1, SC3, -1, kas.CHI, h) ;
	SC3->name = "S-Casimir cube" ;
	niceShow (SC3) ;

  exit (0);


  

  KasimirLower3tensor (&kas) ;
  KasimirUpper3tensor (&kas) ;
  KasimirOperatorK3 (&kas) ;
  
  exit (0) ;
} /* Kasimirs */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/*****  SU(2/1) representation theory. This is used by, but does not depend on the analysis above of the Feynman diagrams **********************/
/*****  Scalar anomaly paper , indecomposable representations submited to Arxiv and JHEP in My 20, 2020 ****************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/


#define NTYPES 12

MX *neq[NTYPES] ;
MX *coq[NTYPES] ;
MX *Coq[NTYPES] ;
MX nchiT[NTYPES] ;
MX nchiS[NTYPES] ;
MX nchiL[NTYPES] ;
MX nchiR[NTYPES] ;

int ss[] = {4,4,4, 8,8,8, 8,8,8, 8,8,8} ;
MX nn[10], ee[10], qq[10], N2[10], E2[10], Q2[10], N2a[10], E2a[10], Q2a[10], N2b[10], E2b[10], Q2b[10] ;
MX nncoq[10], eecoq[10], eeCoq[10], qqcoq[10] ;
MX chiT, chiS, chiL, chiR ;
MX chiT2, chiS2, chiL2, chiR2 ;
MX SG[4], SB[4] ;
float gg[4][4] = {{-1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}} ;
complex float eps[4][4][4][4] ;
complex float PP[4][4][4][4] ;
complex float PM[4][4][4][4] ;

static void niceShow (MX a) ;
static int c3Mask = 0 ;
static int SU3 = 0 ;
static int myType = -1 ;
static MX muHermite (MX m, AC_HANDLE h)
{
  MX m1 = 0 ;
  m1 = mxMatTranspose (m1, m, h) ;
  mxConjugate (m1, m1, h) ;
  
  return m1 ;
} /* muHermite */

/*************************************************************************************************/

typedef struct fabcStruct {int a,b,c,sign; complex float z; const char *title ;} FABC ;
static void muStructure (void)
{
  AC_HANDLE h = ac_new_handle () ;
  int a,b,c, t ;
  MX mm1, mm2, mm3, mm4 ;
  FABC *f ;
  FABC ff[] = {
    {1,2,3,-1,2.0I, "[m_1,m_2] = 2i m_3"},
    {2,3,1,-1,2.0I, "[m_2,m_3] = 2i m_1"},
    {3,1,2,-1,2.0I, "[m_3,m_1] = 2i m_2"},

    {1,4,7,-1, I,"\n  [m_1,m_4] =  i m_7"},
    {1,5,6,-1,-I,    "[m_1,m_5] = -i m_6"},
    {1,6,5,-1, I,    "[m_1,m_6] =  i m_5"},
    {1,7,4,-1,-I,    "[m_1,m_7] = -i m_4"},

    {2,4,6,-1, I,"\n  [m_2,m_4] =  i m_6"},
    {2,5,7,-1, I,    "[m_2,m_5] =  i m_7"},
    {2,6,4,-1,-I,    "[m_2,m_6] = -i m_4"},
    {2,7,5,-1,-I,    "[m_2,m_7] = -i m_5"},

    {3,4,5,-1, I,"\n  [m_3,m_4] =  i m_5"},
    {3,5,4,-1,-I,    "[m_3,m_5] = -i m_4"},
    {3,6,7,-1,-I,    "[m_3,m_6] = -i m_7"},
    {3,7,6,-1, I,    "[m_3,m_7] =  i m_6"},

    {0,4,5,-1, I,"\n  [m_0,m_4] =  i m_5"},
    {0,5,4,-1,-I,    "[m_0,m_5] = -i m_4"},
    {0,6,7,-1, I,    "[m_0,m_6] =  i m_7"},
    {0,7,6,-1,-I,    "[m_0,m_7] = -i m_6"},

    {4,4,9, 1,1,"\n  {m_4,m_4} = -m_9"},
    {5,5,9, 1,1,    "{m_5,m_5} = -m_9"},

    {6,6,8, 1,1,   "{m_6,m_6} = -m_8"},
    {6,7,0, 1, 0,    "{m_6,m_7} = 0  "},
    {7,7,8, 1,1,   "{m_7,m_7} = -m_8"},
    {4,5,0, 1, 0,"\n  {m_4,m_5} = 0  "},
    {6,7,0, 1, 0,    "{m_6,m_7} = 0  "},

    {4,6,1, 1, -1,"\n  {m_4,m_6} = m_1"},
    {5,6,2, 1, -1,    "{m_5,m_6} = m_2"},
    {4,7,2, 1,1,"\n  {m_4,m_7} = -m_2"},
    {5,7,1, 1, -1,    "{m_5,m_7} = m_1"},


    /*

    */
    {-1,0,0,0,0}
  } ;

  for (f = ff ; f->a >=0 ; f++)
    {
      a = f->a ;
      b = f->b ;
      c = f->c ;

      printf ("# %s\t", f->title) ;
      for (t = 0 ; t < NTYPES ; t++)
	{
	  double z ;
	  if (0 && t != 11) continue ;
	  mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
	  mm2  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
	  mm3  = mxCreate (h,  "mm3", MX_COMPLEX, ss[t], ss[t], 0) ;
	  mm4  = mxCreate (h,  "mm4", MX_COMPLEX, ss[t], ss[t], 0) ;

	  mm1 = mxMatMult (neq[t][a], neq[t][b], h) ;
	  mm2 = mxMatMult (neq[t][b], neq[t][a], h) ;
	  mm3 = mxLinearCombine (mm3, 1, mm1, f->sign, mm2, h) ;
	  mm4 = mxLinearCombine (mm4, 1, mm3, f->z,  neq[t][c], h) ;
	  
	  if (0 && a == 2 && b == 7 && t == 7)
	    {
	      printf ("\n# mu(%d) type %d\n", a, t) ;
	      niceShow  (neq[t][a]) ;
	      printf ("\n# mu(%d) type %d\n", b, t) ;
	      niceShow  (neq[t][b]) ;
	      printf ("\n# [a,b] type %d\n", t) ;
	      niceShow (mm3) ;
	      printf ("\n# [a,b] should be equal to neq[%d][%d]\n", t,c) ;
	      niceShow (neq[t][c]) ;
	      printf ("\n# norm");
	    }
	  z = mxFNorm (mm4) ;
	  if (z < .0000001) z = 0 ;
	  printf ("\t%.2g", z) ;
	  }
	printf ("\n") ;
      }
 
  if (1)
    {
      niceShow (neq[6][6]) ;
      niceShow (neq[6][7]) ;
    }
  ac_free (h) ;
  return ;
} /* muStructure */

/*************************************************************************************************/

static void muSigma (AC_HANDLE h)
{
  int i, j, k, l, n ;
  float z ;

  memset (eps, 0, sizeof(eps)) ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  { /* checked, this is correct */
	    n =(i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l) ;
	    if (n > 0)
	      eps[i][j][k][l] = n = 1 ;
	    else if (n < 0)
	      eps[i][j][k][l] = n = -1 ;
	    if (0) if (n) printf ("espison(%d,%d,%d,%d) = %d\n", i,j,k,l,n) ;
	  }


  complex float sg0[] = {1,0,0,1} ;
  complex float sb0[] = {-1,0,0,-1} ;
  complex float sg1[] = {0,1,1,0} ;
  complex float sg2[] = {0,I,-I,0} ;
  complex float sg3[] = {1,0,0,-1} ;

  SG[0] = mxCreate (h, "Sigma_0", MX_COMPLEX, 2, 2, 0) ;
  SG[1] = mxCreate (h, "Sigma_1", MX_COMPLEX, 2, 2, 0) ;
  SG[2] = mxCreate (h, "Sigma_2", MX_COMPLEX, 2, 2, 0) ;
  SG[3] = mxCreate (h, "Sigma_3", MX_COMPLEX, 2, 2, 0) ;

  SB[0] = mxCreate (h, "SB_0", MX_COMPLEX, 2, 2, 0) ;
  SB[1] = mxCreate (h, "SB_1", MX_COMPLEX, 2, 2, 0) ;
  SB[2] = mxCreate (h, "SB_2", MX_COMPLEX, 2, 2, 0) ;
  SB[3] = mxCreate (h, "SB_3", MX_COMPLEX, 2, 2, 0) ;

  mxSet (SG[0], sg0) ;
  mxSet (SB[0], sb0) ;
  mxSet (SG[1], sg1) ;
  mxSet (SB[1], sg1) ;
  mxSet (SG[2], sg2) ;
  mxSet (SB[2], sg2) ;
  mxSet (SG[3], sg3) ;
  mxSet (SB[3], sg3) ;

  printf ("### Verify that the sigma sigma-bar obey the Clifford algebra sg_i sb_j + sg_j sb_i = 2 g_ij Identity[2] \n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	AC_HANDLE h = ac_new_handle () ;
	MX mmm[3] ; 
	MX mm1 =  mxCreate (h, "ij", MX_COMPLEX, 2, 2, 0) ;
	MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm3 =  mxCreate (h, "ij+ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm4 =  mxCreate (h, "zero", MX_COMPLEX, 2, 2, 0) ;

	mm1 = mxMatMult (SG[i], SB[j], h) ;
	mm2 = mxMatMult (SG[j], SB[i], h) ;

	mmm[0] = SG[i] ;
	mmm[1] = SB[j] ;
	mmm[2] = 0 ;
	mm1 = mxMatListMult (h, mmm) ;
	if (0)
	  {
	    printf ("###### sigma sbar : i=%d j=%d \n", i, j) ;
	    niceShow (mm1) ;
	  }

	mmm[0] = SG[j] ;
	mmm[1] = SB[i] ;
	mmm[2] = 0 ;
	mm2 = mxMatListMult (h, mmm) ;
	if (0) niceShow (mm2) ;
	mm3 = mxLinearCombine (mm3, 1,mm1, 1, mm2, h) ; 
	mm4 = mxLinearCombine (mm4, 1,mm3, -2*gg[i][j], SG[0], h) ; 
	if (0)
	  {
	    niceShow (mm1) ;
	    niceShow (mm2) ;	
	    niceShow (mm3) ;	
	    niceShow (mm4) ;
	  }
	z = mxFNorm(mm4) ;
	if (z > .001)
	  printf ("###### sigma sbar : i=%d j=%d {i,j} = 2 g_ij Id :: verif %g\n",i,j, z) ;
	
	ac_free (h) ;	
      }
  if (0) exit (0) ;

  printf ("### Verify that the sigma sigma-bar obey the Clifford algebra sby_i sg_j + sb_j sg_i = 2 g_ij Identity[2] \n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	AC_HANDLE h = ac_new_handle () ;

	MX mm1 =  mxCreate (h, "ij", MX_COMPLEX, 2, 2, 0) ;
	MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm3 =  mxCreate (h, "ij+ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm4 =  mxCreate (h, "zero", MX_COMPLEX, 2, 2, 0) ;

	mm1 = mxMatMult (SB[i], SG[j], h) ;
	mm2 = mxMatMult (SB[j], SG[i], h) ;
	mm3 = mxLinearCombine (mm3, 1,mm1, 1, mm2, h) ; 

	mm4 = mxLinearCombine (mm4, 1,mm3, -2*gg[i][j], SG[0], h) ; 
	if (0)
	  {
	    niceShow (mm1) ;
	    niceShow (mm2) ;	
	    niceShow (mm3) ;	
	    niceShow (mm4) ;
	  }
	z = mxFNorm(mm4) ;
	if (z > .001)
	  printf ("###### sbar sigma : i=%d j=%d {i,j} = 2 g_ij Id :: verif %g\n",i,j, z) ; 
	
	ac_free (h) ;	
      }

  /* check the projectors */
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    PP[i][j][k][l] = (gg[i][k]*gg[j][l] - gg[i][l]*gg[j][k] + I * eps[i][j][k][l])/4.0 ;
	    PM[i][j][k][l] = (gg[i][k]*gg[j][l] - gg[i][l]*gg[j][k] - I * eps[i][j][k][l])/4.0 ;
	  }

  printf("### Verify that PP is a projector PP^2 = PP\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    float z ;
	    complex float z2 = 0, z1 = PP[i][j][k][l] ;
	    int a, b ;
	    for (a = 0 ; a < 4 ; a++)
	      for (b = 0 ; b < 4 ; b++)
		z2 += PP[i][j][a][b] *gg[a][a] * gg[b][b] * PP[a][b][k][l] ;
	    z = cabsf (z2 - z1) ;
	    if (z > minAbs)
	      printf("PP PP - PP not zero ijkl = %d %d %d %d  zz=%g\n", i,j,k,l,z) ;
	  }

  printf("### Verify that PM is a projector PM^2 = PM\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    float z ;
	    complex float z2 = 0, z1 = PM[i][j][k][l] ;
	    int a, b ;
	    for (a = 0 ; a < 4 ; a++)
	      for (b = 0 ; b < 4 ; b++)
		z2 += PM[i][j][a][b] *gg[a][a] * gg[b][b] * PM[a][b][k][l] ;
	    z = cabsf (z2 - z1) ;
	    if (z > minAbs)
	      printf("PM PM - PM not zero ijkl = %d %d %d %d  zz=%g\n", i,j,k,l,z) ;
	  }

  printf("### Verify that PP is a projector PP PM = 0\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    float z, z2 = 0, z1 = 0 ;
	    int a, b ;
	    for (a = 0 ; a < 4 ; a++)
	      for (b = 0 ; b < 4 ; b++)
		z2 += PP[i][j][a][b] *gg[a][a] * gg[b][b] * PM[a][b][k][l] ;
	    z = fabsf (z2 - z1) ;
	    if (z > minAbs)
	      printf("PP PM  not zero ijkl = %d %d %d %d  zz=%g\n", i,j,k,l,z) ;
	  }

  printf("### Verify that SG SB = PP SG SB\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	AC_HANDLE h = ac_new_handle () ;

	MX mm1 =  mxCreate (h, "ij", MX_COMPLEX, 2, 2, 0) ;
	MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm3 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	int a, b ;

	mm2 = mxMatMult (SG[i], SB[j], h) ;
	mm3 = mxMatMult (SG[j], SB[i], h) ;
	mm1 = mxLinearCombine (mm1, 0.5,mm2, -0.5, mm3, h) ;
	if (0)
	  {
	    printf ("## i=%d j=%d ::\n", i, j) ;
	    niceShow (mm1) ;
	  }
	mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;

	for (a = 0 ; a < 4 ; a++)
	  for (b = 0 ; b < 4 ; b++)
	    {
	      MX mm ;
	      mm =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	      mm = mxMatMult (SG[a], SB[b], h) ;
	      mm3 =  mxCreate (h, "mm2", MX_COMPLEX, 2, 2, 0) ;
	      mm3 = mxLinearCombine (mm3, 1,mm2, PP[i][j][a][b], mm, h) ; 
	      mm2 = mm3 ;
	    }
	mm3 =  mxCreate (h, "zero", MX_COMPLEX, 2, 2, 0) ;
	if (0) niceShow (mm2) ;
	mm3 = mxLinearCombine (mm3, 1,mm1, -1*gg[i][i]*gg[j][j], mm2, h) ;
	z = mxFNorm(mm3) ;
	if (z > minAbs)
	  printf ("###### S_i sb_j not equal PP s sb: sbar i=%d j=%d z = %f\n", i,j,z) ;
	
	ac_free (h) ;	
      }


  printf("### Verify that SB SG = PM SB SG\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      {
	AC_HANDLE h = ac_new_handle () ;

	MX mm1 =  mxCreate (h, "ij", MX_COMPLEX, 2, 2, 0) ;
	MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	MX mm3 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	int a, b ;

	mm2 = mxMatMult (SB[i], SG[j], h) ;
	mm3 = mxMatMult (SB[j], SG[i], h) ;
	mm1 = mxLinearCombine (mm1, 0.5,mm2, -0.5, mm3, h) ;
	if (0)
	  {
	    printf ("## i=%d j=%d ::\n", i, j) ;
	    niceShow (mm1) ;
	  }
	mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;

	for (a = 0 ; a < 4 ; a++)
	  for (b = 0 ; b < 4 ; b++)
	    {
	      MX mm ;
	      mm =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	      mm = mxMatMult (SB[a], SG[b], h) ;
	      mm3 =  mxCreate (h, "mm2", MX_COMPLEX, 2, 2, 0) ;
	      mm3 = mxLinearCombine (mm3, 1,mm2, PM[i][j][a][b], mm, h) ; 
	      mm2 = mm3 ;
	    }
	mm3 =  mxCreate (h, "zero", MX_COMPLEX, 2, 2, 0) ;
	if (0) niceShow (mm2) ;
	mm3 = mxLinearCombine (mm3, 1,mm1, -1*gg[i][i]*gg[j][j], mm2, h) ;
	z = mxFNorm(mm3) ;
	if (z > minAbs)
	  printf ("###### S_i sb_j not equal PP s sb: sbar i=%d j=%d z = %f\n", i,j,z) ;
	
	ac_free (h) ;	
      }

  printf("### Verify that Tr(SG_i SB_j SG_k SB_l = 2 * (g_ijg_kl - g_ik_g_jl+gil_gjk + I eps_ijkl)\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    float z ;
	    float complex z1, z2 ;
	    AC_HANDLE h = ac_new_handle () ;
	    MX mmm[5] ;

	    MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 3, 0) ;
	    MX mm3 =  mxCreate (h, "ji", MX_COMPLEX, 3, 5, 0) ;
	    MX mm4 =  0 ;

	    mm4 = mxMatMult (mm2, mm3, h) ;



	    mmm[0] = mm2 ;
	    mmm[1] = mm3 ;
	    mmm[2] = 0 ;
	    mmm[4] = 0 ;


	    mm4 = mxMatListMult (h, mmm) ;               


	    mmm[0] = SG[i] ;
	    mmm[1] = SB[j] ;
	    mmm[2] = SG[k] ;
	    mmm[3] = SB[l] ;
	    mmm[4] = 0 ;


	    mm4 = mxMatListMult (h, mmm) ;               
	    z1 = mxMatTrace (mm4) ;
	    z2 = 2*(gg[i][j]*gg[k][l] - gg[i][k]*gg[j][l] + gg[i][l]*gg[j][k] + I*eps[i][j][k][l]) ;
	    z = cabs (z2 - z1) ;
	    if (z > minAbs)
	      printf ("###### Trace (sigma ijkl) not equal gg - gg + gg + i epsilon: i=%d j=%d k=%d l=%d z = %f\n", i,j,k,l,z) ;
	    
	    ac_free (h) ;	
	  }

  printf("### Verify that Tr(SB_i SG_j SB_k SG_l = 2 * (g_ijg_kl - g_ik_g_jl+gil_gjk - I eps_ijkl)\n") ;
  for (i = 0 ; i < 4 ; i++)
    for (j = 0 ; j < 4 ; j++)
      for (k = 0 ; k < 4 ; k++)
	for (l = 0 ; l < 4 ; l++)
	  {
	    float z ;
	    float complex z1, z2 ;
	    AC_HANDLE h = ac_new_handle () ;

	    MX mm2 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	    MX mm3 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;
	    MX mm4 =  mxCreate (h, "ji", MX_COMPLEX, 2, 2, 0) ;

	    mm2 = mxMatMult (SB[i], SG[j], h) ;
	    mm3 = mxMatMult (mm2, SB[k], h) ;
	    mm4 = mxMatMult (mm3, SG[l], h) ;

	    z1 = mxMatTrace (mm4) ;
	    z2 = 2*(gg[i][j]*gg[k][l] - gg[i][k]*gg[j][l] + gg[i][l]*gg[j][k] - I*eps[i][j][k][l]) ;
	    z = cabs (z2 - z1) ;
	    if (z > minAbs)
	      printf ("###### Trace (sb sg ijkl) not equal 2 *( gg - gg + gg - i epsilon): i=%d j=%d k=%d l=%d z = %f\n", i,j,k,l,z) ;
	    
	    ac_free (h) ;	
	  }



  return ;
} /* muSigma */

/*************************************************************************************************/

static void muInit (AC_HANDLE h)
{
  int t, i ;
  float s2 = sqrt (2) ;
  float s3 = sqrt (3) ;
  MX mm =  mxCreate (h, "mm", MX_COMPLEX, 4, 4, 0) ;
  MX mm2 = 0 ;

  complex float mu1[] = {0,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,0} ;
  complex float mu2[] = {0,0,0,0, 0,0,-I,0, 0,I,0,0, 0,0,0,0} ;
  complex float mu3[] = {0,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,0,0} ;

  complex float mu0n[] = {1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,-1} ;
  complex float mu0e[] = {0,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0, -2} ;
  complex float mu0SU3[] = {0,0,0,0, 0,1/s3,0,0, 0,0,1/s3,0, 0,0,0, -2/s3} ;
  complex float mu0q[] = {4/3.0,0,0,0, 0,1/3.0,0,0, 0,0,1/3.0,0, 0,0,0,-2/3.0} ;

  complex float mu8n[] = {1,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,0,-1} ;
  complex float mu8e[] = {0,0,0,0, 0,0,0,0, 0,0,-2,0, 0,0,0,-2} ;
  complex float mu8q[] = {4/3.0,0,0,0, 0,4/3.0,0,0, 0,0,-2/3.0,0, 0,0,0,-2/3.0} ;

  complex float mu9n[] = {1,0,0,0, 0,-1,0,0, 0,0,1,0, 0,0,0,-1} ;
  complex float mu9e[] = {0,0,0,0, 0,-2,0,0, 0,0,0,0, 0,0,0,-2} ;
  complex float mu9q[] = {4/3.0,0,0,0, 0,-2/3.0,0,0, 0,0,4/3.0,0, 0,0,0,-2/3.0} ;

  complex float mu4n[] = {0,0,-1/s2,0, 0,0,0,1/s2, 1/s2,0,0,0, 0,1/s2,0,0} ;
  complex float mu5n[] = {0,0,I/s2,0, 0,0,0,-I/s2, I/s2,0,0,0, 0,I/s2,0,0} ;
  complex float mu6n[] = {0,1/s2,0,0, -1/s2,0,0,0, 0,0,0,1/s2, 0,0,1/s2,0} ;
  complex float mu7n[] = {0,-I/s2,0,0, -I/s2,0,0,0, 0,0,0,-I/s2, 0,0,I/s2,0} ;

  complex float mu4e[] = {0,0,0,0, 0,0,0,1, 0,0,0,0, 0,1,0,0} ;
  complex float mu5e[] = {0,0,0,0, 0,0,0,-I, 0,0,0,0, 0,I,0,0} ;
  complex float mu6e[] = {0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0} ;
  complex float mu7e[] = {0,0,0,0, 0,0,0,0, 0,0,0,-I, 0,0,I,0} ;

  complex float mu4q[] = {0,0,-s2/s3,0, 0,0,0,1/s3, s2/s3,0,0,0, 0,1/s3,0,0} ;
  complex float mu5q[] = {0,0,I*s2/s3,0, 0,0,0,-I/s3, I*s2/s3,0,0,0, 0,I/s3,0,0} ;
  complex float mu6q[] = {0,s2/s3,0,0, -s2/s3,0,0,0, 0,0,0,1/s3, 0,0,1/s3,0} ;
  complex float mu7q[] = {0,-I*s2/s3,0,0, -I*s2/s3,0,0,0, 0,0,0,-I/s3, 0,0,I/s3,0} ;

  complex float xT[] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1} ;
  complex float xS[] = {-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1} ;
  complex float xL[] = { 0,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,0} ;
  complex float xR[] = { 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1} ;


  complex float coq0n[] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1} ;
  complex float coq4n[] = {0,0,1/s2,0, 0,0,0,1/s2, -1/s2,0,0,0, 0,1/s2,0,0} ;
  complex float coq5n[] = {0,0,I/s2,0, 0,0,0,I/s2, I/s2,0,0,0, 0,-I/s2,0,0} ;
  complex float coq6n[] = {0,-1/s2,0,0, 1/s2,0,0,0, 0,0,0,1/s2, 0,0,1/s2,0} ;
  complex float coq7n[] = {0,-I/s2,0,0, -I/s2,0,0,0, 0,0,0,I/s2, 0,0,-I/s2,0} ;

  
  complex float coq0e[] = {2*s2/3,0,0,0, 0,2*s2/3,0,0, 0,0,2*s2/3,0, 0,0,0,2*s2/3} ;
  complex float coq4e[] = {0,0,0,0, 0,0,0,1, -1,0,0,0, 0,1,0,0} ;
  complex float coq5e[] = {0,0,0,0, 0,0,0,-I, -I,0,0,0, 0,I,0,0} ;
  complex float coq6e[] = {0,0,0,0, 1,0,0,0, 0,0,0,1, 0,0,1,0} ;
  complex float coq7e[] = {0,0,0,0, I,0,0,0, 0,0,0,-I, 0,0,I,0} ;

  complex float Coq4e[] = {0,0,-2,0, 0,0,0,1, 0,0,0,0, 0,1,0,0} ;
  complex float Coq5e[] = {0,0,2*I,0, 0,0,0,-I, 0,0,0,0, 0,I,0,0} ;
  complex float Coq6e[] = {0,2,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0} ;
  complex float Coq7e[] = {0,-2*I,0,0,0,0,0,0, 0,0,0,-I, 0,0,I,0} ;

  /* ERROR in eq appendix H.2 of Scalar paper:
   * in the paper we should replace sqrt(2) by -sqrt(2) in equation H.2
   * there is probably a related error of sign in H.3
   * No conclusion is modified
   */
  
  complex float coq0q[] = {2*s2/3,0,0,0, 0,2*s2/3,0,0, 0,0,2*s2/3,0, 0,0,0,2*s2/3} ;
  complex float coq4q[] = {0,0,1/s3,0, 0,0,0,s2/s3, -1/s3,0,0,0, 0,s2/s3,0,0} ;
  complex float coq5q[] = {0,0,I/s3,0, 0,0,0,I*s2/s3, I/s3,0,0,0, 0,-I*s2/s3,0,0} ;
  complex float coq6q[] = {0,-1/s3,0,0, 1/s3,0,0,0, 0,0,0,s2/s3, 0,0,s2/s3,0} ;
  complex float coq7q[] = {0,-I/s3,0,0, -I/s3,0,0,0, 0,0,0,I*s2/s3, 0,0,-I*s2/s3,0} ;

  
  chiT = mxCreate (h, "chiT", MX_COMPLEX, 4, 4, 0) ;
  chiS = mxCreate (h,  "chi", MX_COMPLEX, 4, 4, 0) ;
  chiL = mxCreate (h, "chiL", MX_COMPLEX, 4, 4, 0) ;
  chiR = mxCreate (h, "chiR", MX_COMPLEX, 4, 4, 0) ;

  mxSet (chiT, xT) ;
  mxSet (chiS, xS) ;
  mxSet (chiL, xL) ;
  mxSet (chiR, xR) ;

  for (t = 0 ; t < 3 ; t++)
    {
      nchiT[t] = chiT ;
      nchiS[t] = chiS ;
      nchiL[t] = chiL ;
      nchiR[t] = chiR ;
    }
  neq[0] = nn ;
  neq[1] = ee ;
  neq[2] = qq ;

  coq[0] = nncoq ;
  coq[1] = eecoq ;
  Coq[1] = eeCoq ;
  coq[2] = qqcoq ;

  for (i = 0 ; i < 10 ; i++)
    {
      nn[i] = mxCreate (h, messprintf ("nn_%d", i), MX_COMPLEX, 4, 4, 0) ;
      ee[i] = mxCreate (h, messprintf ("ee_%d", i), MX_COMPLEX, 4, 4, 0) ;
      qq[i] = mxCreate (h, messprintf ("qq_%d", i), MX_COMPLEX, 4, 4, 0) ;
      coq[0][i] = mxCreate (h, messprintf ("coq_%d", i), MX_COMPLEX, 4, 4, 0) ;
      coq[1][i] = mxCreate (h, messprintf ("coq_%d", i), MX_COMPLEX, 4, 4, 0) ;
      Coq[1][i] = mxCreate (h, messprintf ("Coq_%d", i), MX_COMPLEX, 4, 4, 0) ;
      coq[2][i] = mxCreate (h, messprintf ("coq_%d", i), MX_COMPLEX, 4, 4, 0) ;
    }
  for (t = 0 ; t < 3 ; t++)
    {
      mxSet (neq[t][1], mu1) ;
      mxSet (neq[t][2], mu2) ;
      mxSet (neq[t][3], mu3) ;
    }
  mxSet (nn[0], mu0n) ;
  if(SU3 == 0)
    mxSet (ee[0], mu0e) ;
  else
    mxSet (ee[0], mu0SU3) ;
  mxSet (qq[0], mu0q) ;

  mxSet (nn[8], mu8n) ;
  mxSet (ee[8], mu8e) ;
  mxSet (qq[8], mu8q) ;

  mxSet (nn[9], mu9n) ;
  mxSet (ee[9], mu9e) ;
  mxSet (qq[9], mu9q) ;

  mxSet (nn[4], mu4n) ;
  mxSet (nn[5], mu5n) ;
  mxSet (nn[6], mu6n) ;
  mxSet (nn[7], mu7n) ;

  mxSet (ee[4], mu4e) ;
  mxSet (ee[5], mu5e) ;
  mxSet (ee[6], mu6e) ;
  mxSet (ee[7], mu7e) ;

  mxSet (qq[4], mu4q) ;
  mxSet (qq[5], mu5q) ;
  mxSet (qq[6], mu6q) ;
  mxSet (qq[7], mu7q) ;


  mxSet (coq[0][0], coq0n) ;
  mxSet (coq[0][4], coq4n) ;
  mxSet (coq[0][5], coq5n) ;
  mxSet (coq[0][6], coq6n) ;
  mxSet (coq[0][7], coq7n) ;
    
  mxSet (coq[1][0], coq0e) ;
  mxSet (coq[1][4], coq4e) ;
  mxSet (coq[1][5], coq5e) ;
  mxSet (coq[1][6], coq6e) ;
  mxSet (coq[1][7], coq7e) ;
  mxSet (coq[1][0], coq0e) ;
  mxSet (Coq[1][4], Coq4e) ;
  mxSet (Coq[1][5], Coq5e) ;
  mxSet (Coq[1][6], Coq6e) ;
  mxSet (Coq[1][7], Coq7e) ;
  

    
  mxSet (coq[2][0], coq0q) ;
  mxSet (coq[2][4], coq4q) ;
  mxSet (coq[2][5], coq5q) ;
  mxSet (coq[2][6], coq6q) ;
  mxSet (coq[2][7], coq7q) ;
    
  if (1) 
    {
      niceShow (qq[1]) ;
      niceShow (qq[2]) ;
      niceShow (qq[3]) ;
      
      niceShow (ee[6]) ;
      niceShow (ee[7]) ;
      niceShow (ee[0]) ;
      niceShow (qq[6]) ;
      niceShow (qq[7]) ;
      niceShow (qq[0]) ;

      niceShow (coq[0][0]) ;
      niceShow (coq[0][4]) ;
      niceShow (coq[0][5]) ;
      niceShow (coq[0][6]) ;
      niceShow (coq[0][7]) ;

	    
      mm = mxMatMult (ee[4],ee[4],h) ;
      niceShow (mm) ;
      
      mm2 = mxMatMult (ee[6],ee[6],h) ;
      niceShow (mm2) ;
    }
} /* muInit */

/*************************************************************************************************/

static MX muComposeMatrix (MX mm, MX a00, MX a01, MX a10, MX a11, complex float x00, complex float x01, complex float x10, complex float x11)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, iMax = 4, iiMax = 8 ;
  int di, dj ;
  const complex float *zc ;
  complex float zz[64] ;

  memset (zz, 0, sizeof (zz)) ;
  if (a00)
    {
      mxValues (a00, 0, 0, &zc) ;
      di = 0 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x00 * zc[iMax * i + j] ;
    }
  if (a01)
    {
      mxValues (a01, 0, 0, &zc) ;
      di = 0 ; dj = 4 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x01 * zc[iMax * i + j] ;
    }
  if (a10) 
    {
      mxValues (a10, 0, 0, &zc) ;
      di = 4 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x10 * zc[iMax * i + j] ;
    }
  if (a11)
    {
      mxValues (a11, 0, 0, &zc) ;
      di = 4 ; dj = 4 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x11 * zc[iMax * i + j] ;
    }

  mxSet (mm, zz) ;
  ac_free (h) ;
  
  return mm ;
} /* muComposeMatrix */

/*************************************************************************************************/

static MX muComposeIntMatrix (int d, MX mm, MX a00, MX a01, MX a10, MX a11, int x00, int x01, int x10, int x11)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, iMax = d, iiMax = 2*d ;
  int di, dj ;
  const int *zi ;
  int zz[4*d*d] ;

  memset (zz, 0, sizeof (zz)) ;
  if (a00)
    {
      mxValues (a00, &zi, 0, 0) ;
      di = 0 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x00 * zi[iMax * i + j] ;
    }
  if (a01)
    {
      mxValues (a01, &zi, 0,0) ;
      di = 0 ; dj = d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x01 * zi[iMax * i + j] ;
    }
  if (a10) 
    {
      mxValues (a10, &zi, 0, 0) ;
      di = d ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x10 * zi[iMax * i + j] ;
    }
  if (a11)
    {
      mxValues (a11, &zi, 0, 0) ;
      di = d ; dj = d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = x11 * zi[iMax * i + j] ;
    }

  mxSet (mm, zz) ;
  ac_free (h) ;
  
  return mm ;
} /* muComposeIntMatrix */

/*************************************************************************************************/

static MX muTripleComposeIntMatrix (int d, MX mm, MX a00, MX a01, MX a02, MX a10, MX a11, MX a12, MX a20, MX a21, MX a22)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, iMax = d, iiMax = 3*d ;
  int di, dj ;
  const int *zi ;
  int zz[9*d*d] ;

  memset (zz, 0, sizeof (zz)) ;
  if (a00)
    {
      mxValues (a00, &zi, 0, 0) ;
      di = 0 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a01)
    {
      mxValues (a01, &zi, 0,0) ;
      di = 0 ; dj = d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a02)
    {
      mxValues (a02, &zi, 0,0) ;
      di = 0 ; dj = 2*d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a10) 
    {
      mxValues (a10, &zi, 0, 0) ;
      di = d ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a11)
    {
      mxValues (a11, &zi, 0, 0) ;
      di = d ; dj = d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a12)
    {
      mxValues (a12, &zi, 0, 0) ;
      di = d ; dj = 2*d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }



   if (a20) 
    {
      mxValues (a20, &zi, 0, 0) ;
      di = 2*d ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }
  if (a21)
    {
      mxValues (a21, &zi, 0, 0) ;
      di = 2*d ; dj = d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = 2*zi[iMax * i + j] ;
    }
  if (a22)
    {
      mxValues (a22, &zi, 0, 0) ;
      di = 2*d ; dj = 2*d ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iiMax * (i + di) + (j + dj)] = zi[iMax * i + j] ;
    }

  mxSet (mm, zz) ;
  ac_free (h) ;
  
  return mm ;
} /* muTripleComposeIntMatrix */

/*************************************************************************************************/

static MX muCoqComposeIntMatrix (int COQ, int d, MX mm, MX mu, MX nu)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, iMax = d, iiMax = COQ*d ;
  int coq, di, dj ;
  const int *zi ;
  int zz[COQ*COQ*d*d] ;

  memset (zz, 0, sizeof (zz)) ;
  if (mu) /* install in the white diagonals the triangular block copies of mu */
    {
      int w = 1 ;
      int diag ;
      int diagMax = nu ? COQ : 1 ; /* if nu==0, just populate the block diagonal: good for SU(2) */
      for (diag = 0 ; diag < diagMax ; diag += 2)
	{
	  w = 1 ;
	  if (diag == 4) w = 24 ;
	  for (coq = 0 ; coq < COQ - diag ; coq++)
	    {
	      mxValues (mu, &zi, 0, 0) ;
	      di = d * (diag + coq) ; dj = d * (coq) ;
	      if (diag == 0) w = 1 ;
	      if (diag == 2 && coq) w = 4*w ; /* gamma */
	      for (i = 0 ; i < iMax ; i++)
		for (j = 0 ; j < iMax ; j++)
		  zz[iiMax * (i + di) + (j + dj)] = w * zi[iMax * i + j] ;
	    }
	}
    }
  if (nu) /* install in the black diagonals the triangular block copies of mu */
    {
      int w = 1 ;
      int diag ;
      int diagMax = nu ? COQ : 1 ; /* if nu==0, just populate the block diagonal */
      for (diag = 1 ; diag < diagMax ; diag += 2)
	{
	  w = 1 ;
	  if (diag == 3) w = 4 ;
	  for (coq = 0 ; coq < COQ - diag ; coq++)
	    {
	      mxValues (nu, &zi, 0, 0) ;
	      di = d * (diag + coq) ; dj = d * (coq) ;
	      if (diag == 1 && coq) w = 2 * w ;
	      if (diag == 3 && coq) w = 8 * w ; /* to be determined */
	      for (i = 0 ; i < iMax ; i++)
		for (j = 0 ; j < iMax ; j++)
		  zz[iiMax * (i + di) + (j + dj)] = w * zi[iMax * i + j] ;
	    }
	}
    }
  mxSet (mm, zz) ;
  ac_free (h) ;
  
  return mm ;
} /* muCoqComposeIntMatrix */

/*************************************************************************************************/

/* extract the Hermitian part of a matrix */
static MX muBiHK (MX a, int sign, AC_HANDLE h)
{
  MX mm = mxCreate (h, "muBiHK", MX_COMPLEX, 4, 4, 0) ;
  if (a)
    {
      int i, j, iMax = 4 ;
      MX at = mxMatTranspose (0, a, h)  ;
      const complex float *za ;
      const complex float *zat ;
      complex float zz[16] ;
      
      memset (zz, 0, sizeof (zz)) ;
      mxValues (a, 0, 0, &za) ;
      mxValues (at, 0, 0, &zat) ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  zz[iMax * (i) + (j)] = 0.5 * (za [iMax * i + j] + sign * conj(zat [iMax * i + j])) ; 
      mxSet (mm, zz) ;

    }
  return mm ;
} /* muBiHK */

/*************************************************************************************************/

/* extract the anti-Hermitian part of a matrix */
static MX muBiH (MX a, AC_HANDLE h)
{
  return muBiHK (a, 1, h) ;
}
static MX muBiK (MX a, AC_HANDLE h)
{
  return muBiHK (a, -1, h) ;
}

/*************************************************************************************************/
/* Mixing 2 famillies, using a pair of angles
 * alpha and beta
 * alpha is the Hermitian angle, it concerns the down quarks
 * beta is the anti-Hermitian angle, it conscerns the up quarks
 *
 * if alpha = beta, or in the electron case
 *   this is just a change of variables global to all the right states
 *   and the representation remains decomposable
 * theta = alpha - beta   could hopefully be the cabbibo angle
 *   It describes the misalignment of te up/c qarks relative to the down/s quarks
 *   We verify here that the mix matrices represent SU(2/1)
 *   We need to verify that the representatin is indecomposable
 *   It has 2 highest weights u_R and c_R
 * and seems to share d_R + s_R with a same phase ? 
 */
static MX muBiComposeMatrix (MX mm, MX a00, MX a01, MX a10, MX a11, int x00, int x01, int x10, int x11)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, iMax = 4, iiMax = 8 ;
  int di, dj ;
  float pi = 3.1415926535 ;
  float alpha = 1*pi/6 ;
  float beta = 1*pi/4 ;

  complex float x00K = x00 * cos (beta) ;
  complex float x01K = x01 * sin (beta) ;
  complex float x10K = x10 * sin (beta) ;
  complex float x11K = x11 * cos (beta) ;
  complex float x00H = x00 * cos (alpha) ;
  complex float x01H = x01 * sin (alpha) ;
  complex float x10H = x10 * sin (alpha) ;
  complex float x11H = x11 * cos (alpha) ;

  const complex float *zcH ;
  const complex float *zcK ;
  complex float zz[64] ;
  MX a00H = muBiH (a00, h) ;
  MX a01H = muBiH (a01, h) ;
  MX a10H = muBiH (a10, h) ;
  MX a11H = muBiH (a11, h) ;
  MX a00K = muBiK (a00, h) ;
  MX a01K = muBiK (a01, h) ;
  MX a10K = muBiK (a10, h) ;
  MX a11K = muBiK (a11, h) ;

  memset (zz, 0, sizeof (zz)) ;
  if (a00)
    {
      mxValues (a00H, 0, 0, &zcH) ;
      mxValues (a00K, 0, 0, &zcK) ;
      di = 0 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  {
	    zz[iiMax * (i + di) + (j + dj)] = 
	      x00H * zcH[iMax * i + j] +
	      x00K * zcK[iMax * i + j] ;
	  }
    }
  if (a01)
    {
      mxValues (a01H, 0, 0, &zcH) ;
      mxValues (a01K, 0, 0, &zcK) ;
      di = 0 ; dj = 4 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  {
	    zz[iiMax * (i + di) + (j + dj)] = 
	      x01H * zcH[iMax * i + j] +
	      x01K * zcK[iMax * i + j] ;
	  }
    }
  if (a10) 
    {
      mxValues (a10H, 0, 0, &zcH) ;
      mxValues (a10K, 0, 0, &zcK) ;
      di = 4 ; dj = 0 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  {
	    zz[iiMax * (i + di) + (j + dj)] = 
	      x10H * zcH[iMax * i + j] +
	      x10K * zcK[iMax * i + j] ;
	  }
    }
  if (a11)
    {
      mxValues (a11H, 0, 0, &zcH) ;
      mxValues (a11K, 0, 0, &zcK) ;
      di = 4 ; dj = 4 ;
      for (i = 0 ; i < iMax ; i++)
	for (j = 0 ; j < iMax ; j++)
	  {
	    zz[iiMax * (i + di) + (j + dj)] = 
	      x11H * zcH[iMax * i + j] +
	      x11K * zcK[iMax * i + j] ;
	  }
    }

  mxSet (mm, zz) ;
  ac_free (h) ;
  
  return mm ;
} /* muBiComposeMatrix */


/*************************************************************************************************/
/* construct the rotated 8x8 mattrices */
static void muInit2 (AC_HANDLE h)
{
  int i, t ;

  chiT2 = mxCreate (h, "chiT", MX_COMPLEX, 8, 8, 0) ;
  chiS2 = mxCreate (h, "chi", MX_COMPLEX, 8, 8, 0) ;
  chiL2 = mxCreate (h, "chiL", MX_COMPLEX, 8, 8, 0) ;
  chiR2 = mxCreate (h, "chiR", MX_COMPLEX, 8, 8, 0) ;

  muComposeMatrix (chiT2, chiT, 0, 0, chiT, 1, 0, 0, 1) ;
  muComposeMatrix (chiS2, chiS, 0, 0, chiS, 1, 0, 0, 1) ;
  muComposeMatrix (chiL2, chiL, 0, 0, chiL, 1, 0, 0, 1) ;
  muComposeMatrix (chiR2, chiR, 0, 0, chiR, 1, 0, 0, 1) ;

  for (t = 3 ; t < NTYPES ; t++)
    {
      nchiT[t] = chiT2 ;
      nchiS[t] = chiS2 ;
      nchiL[t] = chiL2 ;
      nchiR[t] = chiR2 ;
    }

  for (i = 0 ; i < 10 ; i++)
    { 
      N2[i] = mxCreate (h, messprintf ("N2_%d", i), MX_COMPLEX, 8, 8, 0) ;
      E2[i] = mxCreate (h, messprintf ("E2_%d", i), MX_COMPLEX, 8, 8, 0) ;
      Q2[i] = mxCreate (h, messprintf ("Q2_%d", i), MX_COMPLEX, 8, 8, 0) ;

      N2a[i] = mxCreate (h, messprintf ("N2a_%d", i), MX_COMPLEX, 8, 8, 0) ;
      E2a[i] = mxCreate (h, messprintf ("E2a_%d", i), MX_COMPLEX, 8, 8, 0) ;
      Q2a[i] = mxCreate (h, messprintf ("Q2a_%d", i), MX_COMPLEX, 8, 8, 0) ;

      N2b[i] = mxCreate (h, messprintf ("N2b_%d", i), MX_COMPLEX, 8, 8, 0) ;
      E2b[i] = mxCreate (h, messprintf ("E2b_%d", i), MX_COMPLEX, 8, 8, 0) ;
      Q2b[i] = mxCreate (h, messprintf ("Q2b_%d", i), MX_COMPLEX, 8, 8, 0) ;
    }
  neq[3] = N2 ;
  neq[4] = E2 ;
  neq[5] = Q2 ;

  neq[6] = N2a ;
  neq[7] = E2a ;
  neq[8] = Q2a ;
  
  neq[9] = N2b ;
  neq[10] = E2b ;
  neq[11] = Q2b ;
  
  /* even matrices, same block diagonal */
  for (t = 0 ; t < 3 ; t++)
    for (i = 0 ; i < 10 ; i++)
      {
	if (i > 3 && i < 8) continue ;
	muComposeMatrix (neq[t+3][i], neq[t][i], 0, 0, neq[t][i], 1, 0, 0, 1) ;
	muComposeMatrix (neq[t+6][i], neq[t][i], 0, 0, neq[t][i], 1, 0, 0, 1) ;
	muComposeMatrix (neq[t+9][i], neq[t][i], 0, 0, neq[t][i], 1, 0, 0, 1) ;
      }

  /* odd matrices block diagonal */
  for (i = 4 ; i < 8 ; i++)
    for (t = 0 ; t < 3 ; t++)
      muComposeMatrix (neq[t+3][i], neq[t][i], 0, 0, neq[t][i], 1, 0, 0, 1) ;

  /* odd matrices block diagonal + bottom corner */
    for (t = 0 ; t < 3 ; t++)
      {
	muComposeMatrix (neq[t+6][4], neq[t][4], 0, neq[t][5], neq[t][4], 1, 0, 1, 1) ;
	muComposeMatrix (neq[t+6][5], neq[t][5], 0, neq[t][4], neq[t][5], 1, 0, -1, 1) ;
	muComposeMatrix (neq[t+6][6], neq[t][6], 0, neq[t][7], neq[t][6], 1, 0, 1, 1) ;
	muComposeMatrix (neq[t+6][7], neq[t][7], 0, neq[t][6], neq[t][7], 1, 0, -1, 1) ;
      }      
  /* odd matrices block diagonal + top corner */
    for (t = 0 ; t < 3 ; t++)
      {
	muBiComposeMatrix (neq[t+9][4], neq[t][4], neq[t][5], neq[t][5], neq[t][4],1,1,1,1) ;
	muBiComposeMatrix (neq[t+9][5], neq[t][5], neq[t][4], neq[t][4], neq[t][5],1,-1,-1,1) ;
	muBiComposeMatrix (neq[t+9][6], neq[t][6], neq[t][7], neq[t][7], neq[t][6],1,1,1,1) ;
	muBiComposeMatrix (neq[t+9][7], neq[t][7], neq[t][6], neq[t][6], neq[t][7],1,-1,-1,1) ;
      }      
  return ; 
} /* muInit2 */

/*************************************************************************************************/
#ifdef JUNK
  

  
  exit (0) ;
    printf("#### extract the OSp(2/1) sub-superalgebbra Lepton Cabibbo \n") ;
    if (1)
      { 
	/* we extract the generators F Y H X E of OSp(2/1) from the generators of SU(2/1)
	 * for SU(2/1) 
	 *         X = (4 + i5)/2, Y = (4-i5)/2, Z = (6+i7)/2, T = (6 - i7)/2
	 * Now we extract the OSp generators by projection of the odd generators of the eightfold way adjoint of SU(2/1) on the SU(2) axis
         *         OX = (X+T)/2    OY=(Z-Y)/2
	 * Finally
	 *         OE= OX OX,   OF = - OY OY, OH = - OX OY - OY OX
	 *  we now write that as a program
	*/
	MX mx =  mxCreate (h, "m1", MX_COMPLEX, 8, 8, 0) ;
	MX my =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX mz =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX mt =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;

	MX OX =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX OY =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX OH =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
			
	mx = mxLinearCombine (mx, 0.5, E2a[4], 0.5I, E2a[5], h) ;
	my = mxLinearCombine (my, 0.5, E2a[4], -0.5I, E2a[5], h) ;
	mz = mxLinearCombine (mz, 0.5, E2a[6], 0.5I, E2a[7], h) ;
	mt = mxLinearCombine (mt, 0.5, E2a[6], -0.5I, E2a[7], h) ;

	OX = mxLinearCombine (OX, 1, mx, 1, mt, h) ;
	OY = mxLinearCombine (OY, 1, mz, -1, my, h) ;
	MX OE = mxMatMult (OX, OX, h) ;
	MX OF = mxMatMult (OY, OY, h) ;
	OF = mxLinearCombine (OF, -1, OF, 0, OF, h) ;
	OH = mxLinearCombine (OH, -1, mxMatMult (OX, OY, h), -1, mxMatMult (OY, OX, h), h) ;

	OX->name = "OX" ;
	OY->name = "OY" ;
	OE->name = "OE" ;
	OF->name = "OF" ;
	OH->name = "OH" ;

	niceShow (OX) ;
	niceShow (OY) ;
	niceShow (OH) ;
	niceShow (OE) ;
	niceShow (OF) ;

	/* Casimir  HH + 2 EF + 2 Fe + 2 XY - 2 YX */
	MX K1 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K2 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K3 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K4 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K5 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX KK =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	
	MX KHH = mxMatMult (OH, OH, h) ;
	MX KEF = mxMatMult (OE, OF, h) ;
	MX KFE = mxMatMult (OF, OE, h) ;
	MX KXY = mxMatMult (OX, OY, h) ;
	MX KYX = mxMatMult (OY, OX, h) ;
	K1 = mxAdd (K1, KEF, KFE, h) ;
	K2 = mxAdd (K2, KHH, K1, h) ;
	K3 = mxAdd (K3, K1, K2, h) ;
	K4 = mxSubstract (K4, KYX, KXY, h) ;
	KK = mxAdd (K5, K3, K4, h) ;

	const complex float *zz4 ;
	complex float zz45[64] ;
	mxValues (K4, 0, 0, &zz4) ;
	memcpy (zz45, zz4, sizeof (zz45)) ;
	for (i = 0 ; i < 8 ; i++)
	  zz45[8*i + i] -= .5 ;
	mxSet (K4, zz45) ;

	KK->name = "Q_Casimir" ;
	KHH->name = "Q_Casimir HH" ;
	K1->name = "Q_Casimir EF-FE" ;
	K2->name = "Q_Casimmir HH + EF+ FE " ;
	K3->name = "Q_Casimmir HH + 2 EF+ FE " ;
	K4->name = "Q_Casimir XY - YX -1/2" ;
	if (0)
	  {
	    niceShow (KHH) ;
	    niceShow (K1) ;
	    niceShow (K2) ;
	    niceShow (K3) ;
	  }
	niceShow (K4) ;
	niceShow (KK) ;

	printf ("Verify that the casimir commutes with X and Y\n") ;
	MX CKX = mxMatMult (KK, OX, h) ;
	MX CXK = mxMatMult (OX, KK, h) ;
	MX Com =  mxCreate (h, "[casimir,X]", MX_COMPLEX, 8, 8, 0) ;
	Com = mxSubstract (Com, CKX, CXK, h) ;
	niceShow (Com) ;
	
	printf ("Verify that the S-simir anticommutes with X and Y\n") ;
	MX SCKX = mxMatMult (K4, OX, h) ;
	MX SCXK = mxMatMult (OX, K4, h) ;
	MX SCom =  mxCreate (h, "{S-casimir,X}", MX_COMPLEX, 8, 8, 0) ;
	SCom = mxAdd (SCom, SCKX, SCXK, h) ;
	niceShow (SCom) ;
	
	printf ("####### OSp(2/1) Lepton-Cabibbo representation done\n") ;
      }
    printf("#### extract the OSp(2/1) sub-superalgebbra Quark Cabibbo \n") ;
    if (1)
      { 
	/* we extract the generators F Y H X E of OSp(2/1) from the generators of SU(2/1)
	 * for SU(2/1) 
	 *        X = (4 + i5)/2, Y = (4-i5)/2, Z = (6+i7)/2, T = (6 - i7)/2
	 * Now we extract the OSp generators by projection of the odd generators of the eightfold way adjoint of SU(2/1) on the SU(2) axis
         *         OX = (X+T)/2    OY=(Z-Y)/2
	 * Finally
	 *         OE= OX OX,   OF = - OY OY, OH = - OX OY - OY OX
	 * we now write that as a program
	*/
	MX mx =  mxCreate (h, "m1", MX_COMPLEX, 8, 8, 0) ;
	MX my =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX mz =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX mt =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;

	MX OX =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX OY =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
	MX OH =  mxCreate (h, "m2", MX_COMPLEX, 8, 8, 0) ;
			
	mx = mxLinearCombine (mx, 0.5, Q2a[4], 0.5I, Q2a[5], h) ;
	my = mxLinearCombine (my, 0.5, Q2a[4], -0.5I, Q2a[5], h) ;
	mz = mxLinearCombine (mz, 0.5, Q2a[6], 0.5I, Q2a[7], h) ;
	mt = mxLinearCombine (mt, 0.5, Q2a[6], -0.5I, Q2a[7], h) ;

	OX = mxLinearCombine (OX, 1, mx, 1, mt, h) ;
	OY = mxLinearCombine (OY, 1, mz, -1, my, h) ;
	MX OE = mxMatMult (OX, OX, h) ;
	MX OF = mxMatMult (OY, OY, h) ;
	OF = mxLinearCombine (OF, -1, OF, 0, OF, h) ;
	OH = mxLinearCombine (OH, -1, mxMatMult (OX, OY, h), -1, mxMatMult (OY, OX, h), h) ;

	OX->name = "QX" ;
	OY->name = "QY" ;
	OE->name = "QE" ;
	OF->name = "QF" ;
	OH->name = "QH" ;

	niceShow (OX) ;
	niceShow (OY) ;
	niceShow (OH) ;
	niceShow (OE) ;
	niceShow (OF) ;
	
	/* Casimir  HH + 2 EF + 2 Fe + 2 XY - 2 YX */
	MX K1 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K2 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K3 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K4 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX K5 =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	MX KK =  mxCreate (h, "KHH", MX_COMPLEX, 8, 8, 0) ;
	
	MX KHH = mxMatMult (OH, OH, h) ;
	MX KEF = mxMatMult (OE, OF, h) ;
	MX KFE = mxMatMult (OF, OE, h) ;
	MX KXY = mxMatMult (OX, OY, h) ;
	MX KYX = mxMatMult (OY, OX, h) ;
	K1 = mxAdd (K1, KEF, KFE, h) ;
	K2 = mxAdd (K2, KHH, K1, h) ;
	K3 = mxAdd (K3, K1, K2, h) ;
	K4 = mxSubstract (K4, KYX, KXY, h) ;
	KK = mxAdd (K5, K3, K4, h) ;

	if (0)
	  {
	    const complex float *zz4 ;
	    complex float zz45[64] ;
	    mxValues (K4, 0, 0, &zz4) ;
	    memcpy (zz45, zz4, sizeof (zz45)) ;
	    for (i = 0 ; i < 8 ; i++)
	      zz45[8*i + i] -= 0.5 ;
	    mxSet (K4, zz45) ;
	  }
	else
	  K4 = mxLinearCombine (K4, 1, K4, -.25, KK, h) ;
	KK->name = "Casimir" ;
	KHH->name = "Q_Casimir HH" ;
	K1->name = "Q_Casimir EF-FE" ;
	K2->name = "Q_Casimmir HH + EF+ FE " ;
	K3->name = "Q_Casimmir HH + 2 EF+ FE " ;
	K4->name = "Q_Casimir XY - YX -1/2" ;
	if (0)
	  {
	    niceShow (KHH) ;
	    niceShow (K1) ;
	    niceShow (K2) ;
	  }
	niceShow (K3) ;

	niceShow (K4) ;
	niceShow (KK) ;

	printf ("Verify that the casimir commutes with X and Y\n") ;
	MX CKX = mxMatMult (KK, OX, h) ;
	MX CXK = mxMatMult (OX, KK, h) ;
	MX Com =  mxCreate (h, "[casimir,X]", MX_COMPLEX, 8, 8, 0) ;
	Com = mxSubstract (Com, CKX, CXK, h) ;
	niceShow (Com) ;
	
	printf ("Verify that the S-casimir anticommutes with X and Y\n") ;
	MX SCKX = mxMatMult (K4, OX, h) ;
	MX SCXK = mxMatMult (OX, K4, h) ;
	MX SCom =  mxCreate (h, "{S-casimir,X}", MX_COMPLEX, 8, 8, 0) ;
	SCom = mxAdd (SCom, SCKX, SCXK, h) ;
	niceShow (SCom) ;
	
	printf ("Compute the square of the S-casimir\n") ;
	MX SC2 = mxMatMult (K4, K4, h) ;
	SC2 = mxLinearCombine (SC2, 8/9.0, SC2, -1, KK, h) ;
	SC2->name = "S-Casimir square" ;
	niceShow (SC2) ;
	
	printf ("Compute the product of the Casimir by the S-casimir Q^3 = Q\n") ;
	MX SC3 = mxMatMult (KK, K4, h) ;
	SC3 = mxLinearCombine (SC3, 1/2.00, SC3, -1, K4, h) ;
	SC3->name = "S-Casimir cube" ;
	niceShow (SC3) ;
	
	printf ("####### OSp(2/1) Quark-Cabibbo representation done\n") ;
      }
#endif

/*************************************************************************************************/
/* construct the triple Coquereaux matrices where the cartan subalgebra is non diagonal */
static void muInitNCoq (int a, int b, int COQ)
{
  KAS kas, kas2, kasQ ;
  AC_HANDLE h = ac_new_handle () ;
  int i, d ;
  MX *mu, *nu, *QQ ;
  MX qmuY, qmuH, qmuE, qmuF, qmuU, qmuV, qmuW, qmuX, qmuK1, qmuK2 ;
  kas.a = a ;
  kas.b = b ;
  kas.h = h ;
  kas2.a = a ;
  kas2.b = b+1 ;
  kas2.h = h ;
  kasQ.h = h ;

  KasimirConstructTypicMatrices (&kas) ;
  KasimirConstructTypicMatrices (&kas2) ;
  kasQ.COQ = COQ ;
  mu = kas.mu ;
  nu = kas2.mu ; 
  QQ = kasQ.mu = (MX *) halloc (10 * sizeof (MX), kas.h) ;

  kasQ.a = kas.a ;
  kasQ.b = kas.b ;
    
  kasQ.d = d = COQ * kas.d ;
  kasQ.d1 = kas.d1 ;
  kasQ.d2 = kas.d2 ;
  kasQ.d3 = kas.d3 ;
  kasQ.d4 = kas.d4 ;
  kasQ.chi = kas.chi ;
	
  kasQ.scale = kas.scale ;
  QQ[0] = qmuY = mxCreate (h,  "qmuY", MX_INT, d, d, 0) ;
  QQ[3] = qmuH = mxCreate (h,  "qmuH", MX_INT, d, d, 0) ;
  QQ[1] = qmuE = mxCreate (h,  "qmuE: E", MX_INT, d, d, 0) ;
  QQ[2] = qmuF = mxCreate (h,  "qmuF", MX_INT, d, d, 0) ;
  QQ[6] = qmuU = mxCreate (h,  "qmuU", MX_INT, d, d, 0) ;
  QQ[7] = qmuV = mxCreate (h,  "qmuV", MX_INT, d, d, 0) ;
  QQ[4] = qmuW = mxCreate (h,  "qmuW", MX_INT, d, d, 0) ;
  QQ[5] = qmuX = mxCreate (h,  "qmuX", MX_INT, d, d, 0) ;
  QQ[8] = qmuK1 = mxCreate (h,  "qmuK1: K1 = {U,V}", MX_INT, d, d, 0) ;
  QQ[9] = qmuK2 = mxCreate (h,  "qmuK2: K2 = {W,X}", MX_INT, d, d, 0) ;
  
  
  /* even and odd matrices, same block diagonal, use nu in the bottom left */
  if (1) /* flip sign of nu[5] */
  {
    nu[5] = mxLinearCombine (nu[5], -1, nu[5], 0, nu[5], h) ;
    nu[7] = mxLinearCombine (nu[7], -1, nu[7], 0, nu[7], h) ;
  }
  for (i = 1 ; i < 4 ; i++) /* block diagonal SU(2) */
  {
    muCoqComposeIntMatrix (COQ, kas.d, QQ[i], mu[i],0) ;
  }
  for (i = 4 ; i < 8 ; i++) /* triangular odd generators */
  {
    muCoqComposeIntMatrix (COQ, kas.d, QQ[i], mu[i],nu[i]) ;
  }

  QQ[8] = qmuK1 = KasCommut(qmuU,qmuV,1,&kasQ) ;
  QQ[9] = qmuK2 = KasCommut(qmuW,qmuX,1,&kasQ) ;
  QQ[0] = qmuY = mxLinearCombine (qmuY, 1, qmuK1, 1, qmuK2, h) ;
  
  if (1) /* rescale */
  {
    int dd = kasQ.d, d2 = dd*dd ;
     int i, yy[dd*dd], s2 = kasQ.scale ;
     const int *xx ;


     mxValues (qmuK1, &xx, 0, 0) ;
    for (i = 0 ; i < d2 ; i++)
      yy[i] = xx[i]/s2 ;
    mxSet (qmuK1, yy) ;
     mxValues (qmuK2, &xx, 0, 0) ;
    for (i = 0 ; i < d2 ; i++)
      yy[i] = xx[i]/s2 ;
    mxSet (qmuK2, yy) ;

    mxValues (qmuY, &xx, 0, 0) ;
    for (i = 0 ; i < d2 ; i++)
      yy[i] = xx[i]/s2 ;
    mxSet (qmuY, yy) ;
  }

  printf ("###### Coquereaux\n") ;
  for (i = 0 ; i < 1 ; i++)
    niceIntShow (QQ[i]) ;

  
  MX zUV = mxMatMult (qmuU, qmuV, h) ;
  MX zVU = mxMatMult (qmuV, qmuU, h) ;
  niceIntShow (qmuY) ;
  niceIntShow (qmuU) ;
  niceIntShow (qmuV) ;
  niceIntShow (zUV) ;
  niceIntShow (zVU) ;
  
    
  KasimirCheckCommutators (&kasQ) ;

  KasimirLowerMetric (&kasQ) ;
  KasimirUpperMetric (&kasQ) ;
  
  KasimirOperatorK2 (&kasQ) ;
  GhostKasimirOperatorK2 (&kasQ) ;

  if (0) KasimirOperatorK4 (&kasQ) ;

	printf ("Verify that the casimir commutes with X and Y\n") ;
	MX CKX = mxMatMult (kasQ.kas2, qmuX, h) ;
	MX CXK = mxMatMult (qmuX, kasQ.kas2, h) ;
	MX Com =  mxCreate (h, "[casimir,X]", MX_COMPLEX, d, d, 0) ;
	Com = mxSubstract (Com, CKX, CXK, h) ;
	niceShow (Com) ;
	
	printf ("Verify that the S-casimir anticommutes with X and Y\n") ;
	MX SCKX = mxMatMult (kasQ.CHI, qmuX, h) ;
	MX SCXK = mxMatMult (qmuX, kasQ.CHI, h) ;
	MX SCom =  mxCreate (h, "{S-casimir,X}", MX_COMPLEX, d, d, 0) ;
	SCom = mxAdd (SCom, SCKX, SCXK, h) ;
	niceShow (SCom) ;
	
	printf ("Compute the square of the S-casimir\n") ;
	MX SC2 = mxMatMult (kasQ.CHI,kasQ.CHI, h) ;
	SC2->name = "S-Casimir square" ;
	niceShow (SC2) ;
	
	printf ("Compute the product of the Casimir by the S-casimir Q^3 = Q\n") ;
	MX SC3 = mxMatMult (kasQ.kas2, kasQ.CHI, h) ;
	SC3->name = "S-Casimir cube" ;
	niceShow (SC3) ;

  KasimirLower3tensor (&kasQ) ;
  KasimirUpper3tensor (&kasQ) ;
  KasimirOperatorK3 (&kasQ) ;

  exit(0) ;
  return ;
} /* muInitNCoq */

/*************************************************************************************************/

void muConjugate (AC_HANDLE h)
{
  int i, t ;
  for (t = 0 ; t < NTYPES ; t++)
    for (i = 0 ; i < 10 ; i++)
      neq[t][i] = muHermite (neq[t][i], h) ;
} /* muConjugate */

/*************************************************************************************************/

static void casimir2 (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  MX mm1 = 0, mm2 = 0, chi = 0, casimir = 0 ;
  int i, j, t ;
  float complex z ;
  float a, b ;

  printf ("%s\n", title) ;
  for (t = 0 ; t < NTYPES ; t++)
    {
      casimir  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
      printf ("# Casimir2 Type %d\n", t) ;
      for (i = 0 ; i < 8 ; i++) 
	for (j = 0 ; j < 8 ; j++)
	{
	  if (i < 4 && j >= 4)
	    continue ;
	  if (i >= 4 && j < 4)
	    continue ;

	  /* compute the coefficient g_ab */
	  chi = nchiS[t] ;
	  mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
	  mm2  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
	  mm1 = mxMatMult (neq[t][i], neq[t][j], h) ;
	  mm2 = mxMatMult (chi, mm1, h) ;
	  z = mxMatTrace (mm2) ;
	  z = z/2.0 ;
	  a = creal (z) ;
	  b = cimag (z) ;
	  
	  if (a*a + b*b < .01)
	    continue ;
	  casimir = mxLinearCombine (casimir, 1, casimir, z, mm1, h) ;
	  if (0 && i == 3 && j == 3)
	    niceShow (casimir) ;
	}
      niceShow (casimir) ;
    }
  printf ("\n") ;
  ac_free (h) ;
} /* casimir2 */

/*************************************************************************************************/

static void casimir3 (const char *title, BOOL isHyper)
{
  AC_HANDLE h = ac_new_handle () ;
  MX mm1 = 0, mm1a = 0, mm1b = 0, mm2 = 0, mm3 = 0, mm4 = 0, chi = 0, casimir = 0 ;
  int i, j, k, t ;
  float complex z ;
  float a, b ;

  printf ("%s\n", title) ;
  for (t = 0 ; t < NTYPES ; t++)
    {
      if (myType != -1 && t != myType)
	continue ;
      casimir  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
      printf ("# Casimir3 Type %d\n", t) ;
      for (i = 0 ; i < 8 ; i++) 
	for (j = 0 ; j < 8 ; j++)
	  for (k = 0 ; k < 8 ; k++)
	    {
	      switch (c3Mask)
		{
		case 333:   /* abc => expect zero */
		  if (i*j*k == 0)
		    continue ;
		  if (i > 3 || j > 3 || k > 3)
		    continue ;
		  break ;
		case 888:  /* 888 */
		  if (i+j+k > 0)
		    continue ;
		  break ;
		case 833:  /* 8aa */
		  if (i*j*k > 0)
		    continue ;
		  if (i+j+k == 0 || i > 3 || j > 3 || k > 3)
		    continue ;
		  break ;
		case 844:   /* 8ij */
		  if (i*j*k > 0)
		    continue ;
		  if (i+j+k ==0 || i*(i-4) < 0 || j*(j-4) < 0 || k*(k-4) < 0 )
		    continue ;
		  break ;
		case 344:   /* aij */
		  if (i*j*k == 0)
		    continue ;
		  if (i < 4 && j < 4 && k < 4)
		    continue ;
		  break ;
		case 444:   /* ijk => expect zero */
		  if (i < 4 || j < 4 || k < 4)
		    continue ;
		  break ;
		}

	      /* compute the coefficient t_abc */
	      chi = SU3 ? nchiT[t] : nchiS[t] ;
	      mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm1a  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm1b  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm2  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm3  = mxCreate (h,  "mm3", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm1a = mxMatMult (neq[t][j], neq[t][k], h) ;
	      mm1b = mxMatMult (neq[t][k], neq[t][j], h) ;
	      z = 1 ;
	      if (SU3 == 0 && ( j>=4 && k >= 4))
		z = -1 ;
	      mm1 =  mxLinearCombine (mm1, 1, mm1a, z, mm1b , h) ;
	      if (0)
		{
		  niceShow (mm1) ;	
		  niceShow (neq[t][i]) ;	
		}
	      mm2 = mxMatMult (neq[t][i], mm1, h) ;
	      mm3 = mxMatMult (chi, mm2, h) ;
	      z = isHyper ? mxMatTrace (mm2) : mxMatTrace (mm3) ;

	      z *= 9 ;
	      z /= 8 ;
	      if (SU3 == 0) z *= 1 ;
	      if (0) 
		niceShow (mm3) ;	

	      if (t%3 == 2) z *= 3 ; /* 3 quark colors */
	      if (SU3 == 0 && i*j*k == 0 && ! isHyper) z = -z ;/* g^00 and (g^00)cube == -1 */
	      a = creal (z) ;
	      b = cimag (z) ;
	      
	      if (0 && a*a + b*b < .000001) 
		continue ;
	      if (0)
		{
		  printf ("# mm1 Type %d [%d,%d,%d] %.2f %.2f\n", t,i,j,k, a,b) ;
		  niceShow (mm2) ;
		}
	      mm4 = casimir ;
	      casimir  = mxCreate (h,  "casimir3", MX_COMPLEX, ss[t], ss[t], 0) ;
	      casimir = mxLinearCombine (casimir, 1, mm4, z, mm2, h) ;

	      if (SU3 == 1 && j==7 && k == 7)
		{
		  printf ("# Casimir3 Type %d [%d,%d,%d]\n", t,i,j,k) ;
		  niceShow (casimir) ;
		}
	    }
      niceShow (casimir) ;
    }
  printf ("\n") ;
  ac_free (h) ;
} /* casimir3 */

/*************************************************************************************************/

static void mu2p (const char *title)
{
  AC_HANDLE h = ac_new_handle () ;
  MX mm1 = 0, mm2 = 0, chi ;
  int i, j, t, pass, ok ;
  BOOL debug = FALSE ;
  float complex zz[NTYPES] ;

  printf ("%s\n", title) ;
  printf ("# Index\t\tN  \te  \tq  \tf \tCheck\tN2   \tE2   \tQ2   \tf2   \tCheck2\tN2a   \tE2a   \tQ2a   \tf2a   \tCheck2a\tN2b   \tE2b   \tQ2b   \tf2b   \tCheck2b") ;
  for (i = 0 ; i < 8 ; i++) 
    for (j = 0 ; j < 8 ; j++)
      for (pass = 0 ; pass < 2 ; pass++)
	{
	  if (pass == 0)
	    ok = 0;
	  if (pass == 1 && ok == 0)
	    continue ;
	  
	  if (pass == 1)
	    printf ("\n(%d,%d)\t", i, j) ;
	  for (t = 0 ; t < NTYPES ; t++)
	    {
	      float complex z = 0 ;
	      float a, b ;

	      if (i < 4)
		chi = nchiS[t] ;
	      else
		chi = nchiL[t] ;

	      if (debug) niceShow (neq[t][i]) ;
	      if (debug) niceShow (neq[t][j]) ;
	      mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
	      mm2  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;

	      mm1 = mxMatMult (chi, neq[t][i], h) ;
	      mm2 = mxMatMult (mm1, neq[t][j], h) ;
	      if (debug) niceShow (mm1) ;
	      if (debug) niceShow (mm2) ;
	      z = mxMatTrace (mm2) ;
	      zz[t] = z ;

	      a = creal (z) ;
	      b = cimag (z) ;
	      if (pass == 0)
		{
		  if (a != 0 || b != 0)
		    ok = 1;
		}
	      else
		{
		  nicePrint ("\t", z) ;
		    if (t%3 == 2)
		      {
			nicePrint ("\t", zz[t-1] + 3 * zz[t]) ;
			nicePrint ("\t", -4 * zz[t-2] + zz[t-1] + 3 * zz[t]) ;
		      }
		}
	    }
	}
  printf ("\n") ;
  ac_free (h) ;
} /* mu2p */

/*************************************************************************************************/

static void mu3p (const char *title, int type)
{
  AC_HANDLE h = ac_new_handle () ;
  MX chi1, chi2 ;
  MX mm1 = 0, mm2 = 0, mm3 = 0 ;
  int i, j, k, t, pass, ok, sign ;
  float complex zz[NTYPES] ;
  BOOL debug = FALSE ;

  printf ("\n%s\n", title) ;
  printf ("# Index\t\tN  \te  \tq  \tf \tCheck\tN2   \tE2   \tQ2   \tf2   \tCheck2\tN2   \tE2   \tQ2   \tf2   \tCheck2a\tN2b   \tE2b   \tQ2b   \tf2b   \tCheck2b") ;
  for (i = 0 ; i < 8 ; i++) 
    for (j = 0 ; j < 8 ; j++)
      for (k = j ; k < 8 ; k++)
	{
	  switch (type)
	    {
	    case 0: /* f-abc */
	    case 1: /* d-abc */
	      if (j < i || i > 3 || j > 3 || k > 3)
		continue ;
	      break ;
	    case 2: /* f-aij vector-scalar */
	    case 20: /* f-aij vector-scalar */
	    case 21: /* f-aij vector-scalar */
	    case 22: /* d-aij vector-scalar */
	    case 23: /* d-aij vector-scalar */
	      if (i > 3 || j <= 3 || k <= 3)
		continue ;
	      break ;
	    case 4: /* f-aij vector-scalar anomaly */
	      if (i > 3 || j <= 3 || k <= 3)
		continue ;
	      break ;
	    case 3: /* f-abi f-ijk should vanish */
	      if (i <= 3 && (j > 3 || k <= 3))
		continue ;
	      if (i > 3 && (j <= 3 || k <= 3))
		continue ;
	      break ;
	    }
	  for (pass = 0 ; pass < 2 ; pass++)
	    {	      
	      if (pass == 0)
		ok = 0;
	      if (pass == 1 && ok == 0)
		continue ;
	      
	      if (pass == 1)
		printf ("\n(%d,%d,%d)\t", i, j, k) ;
	      for (t = 0 ; t < NTYPES ; t++)
		{
		  float complex z = 0 ;
		  float a, b ;
		  if (debug) niceShow (neq[t][i]) ;
		  if (debug) niceShow (neq[t][j]) ;
		  switch (type)
		    {
		    case 0: /* f-abc symmetrize in mu-nu, skew in bc, use trace:  (L+R) (abc - acb) */
		      chi1 = nchiT[t] ;
		      chi2 = nchiT[t] ;
		      sign = -1 ;
		      break ;
		    case 1: /* d-abc skew-symmetrize in mu-nu, sym in bc, use super trace:  (L-R) (abc + acb) */
		      chi1 = nchiS[t] ;
		      chi2 = nchiS[t] ;
		      sign = 1 ;
		      break ;
		    case 2: /* f-aij use Laij - Raji */
		      chi1 = nchiL[t] ;
		      chi2 = nchiR[t] ;
		      sign = -1 ;
		      break ;
		    case 20: /* f-aij use Trace aij - aji, expect zero */
		      chi1 = nchiT[t] ;
		      chi2 = nchiT[t] ;
		      sign = -1 ;
		      break ;
		    case 21: /* f-aij use STrace aij - aji, expect zero */
		      chi1 = nchiS[t] ;
		      chi2 = nchiS[t] ;
		      sign = -1 ;
		      break ;
		    case 22: /* d-aij use Trace aij - aji, expect zero */
		      chi1 = nchiT[t] ;
		      chi2 = nchiT[t] ;
		      sign = +1 ;
		      break ;
		    case 23: /* d-aij use STrace aij - aji, expect zero */
		      chi1 = nchiS[t] ;
		      chi2 = nchiS[t] ;
		      sign = +1 ;
		      break ;
		    case 3: /* should be null */
		      chi1 = nchiT[t] ;
		      chi2 = nchiT[t] ;
		      sign = 0 ;
		      break ;
		    case 4: /* should be null */
		      chi1 = nchiT[t] ;
		      chi2 = nchiT[t] ;
		      sign = -1 ;
		      break ;
		    }

		  mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
		  mm2  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
		  mm3  = mxCreate (h,  "mm3", MX_COMPLEX, ss[t], ss[t], 0) ;

		  mm1 = mxMatMult (chi1, neq[t][i], h) ;
		  mm2 = mxMatMult (mm1, neq[t][j], h) ;
		  mm3 = mxMatMult (mm2, neq[t][k], h) ;
		  if (debug) niceShow (mm1) ;
		  if (debug) niceShow (mm2) ;
		  z = mxMatTrace (mm3) ;
		  
		  mm1 = mxMatMult (chi2, neq[t][i], h) ;
		  mm2 = mxMatMult (mm1, neq[t][k], h) ;
		  mm3 = mxMatMult (mm2, neq[t][j], h) ;
		  z += sign * mxMatTrace (mm3) ;
		  
		  zz[t] = z ; /* memorize, to be able to compute the Family e + 3*q */
		  a = creal (z) ;
		  b = cimag (z) ;
		  if (pass == 0)
		    {
		      if (a != 0 || b != 0)
			ok = 1;
		    }
		  else
		    {
		      nicePrint ("\t", z) ;
		      if (t%3 == 2)  /* compute the family vertex */
			{
			  nicePrint ("\t", zz[t-1] + 3 * zz[t]) ;
			  nicePrint ("\t", -4 * zz[t-2] + zz[t-1] + 3 * zz[t]) ;
			}
		    }
		}
	    }
	}
  printf ("\n\n") ;
  ac_free (h) ;
} /* mu3p */

/*************************************************************************************************/

static float complex tetraTrace (MX chi, int t, int i, int j, int k, int l)
{
  float complex z = 0 ;
  MX mm1, mm2, mm3, mm4 ;
  AC_HANDLE h = ac_new_handle () ;

  mm1  = mxCreate (h,  "mm1", MX_COMPLEX, ss[t], ss[t], 0) ;
  mm2  = mxCreate (h,  "mm2", MX_COMPLEX, ss[t], ss[t], 0) ;
  mm3  = mxCreate (h,  "mm3", MX_COMPLEX, ss[t], ss[t], 0) ;
  mm4  = mxCreate (h,  "mm4", MX_COMPLEX, ss[t], ss[t], 0) ;

  mm1 = mxMatMult (chi, neq[t][i], h) ;
  mm2 = mxMatMult (mm1, neq[t][j], h) ;
  mm3 = mxMatMult (mm2, neq[t][k], h) ;
  mm4 = mxMatMult (mm3, neq[t][l], h) ;
  z = mxMatTrace (mm4) ;
  
  ac_free (h) ;
  return z ;
} /* tetraTrace */

/******************/

static void mu4p (const char *title, int type)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, k, l, t, a, b, c, d, pass, ok, mult ;
  float complex zz[NTYPES] ;
  
  printf ("\n%s\n", title) ;
  printf ("# Index  \t\tN  \te  \tq  \tf \tCheck\tN2   \tE2   \tQ2   \tf2   \tCheck2\tN2   \tE2   \tQ2   \tf2   \tCheck2\tN2a   \tE2a   \tQ2a   \tf2a   \tCheck2a\tN2b   \tE2b   \tQ2b   \tf2b   \tCheck2b") ;
  for (a = 0 ; a < 6 ; a+=1) 
    for (b = 0 ; b < 8 ; b+=1)
      for (c = 0 ; c < 8 ; c+=1)
	for (d = 0 ; d < 8 ; d+=1)
	  {
	    if (0 && (a-b)*(a-c)*(a-d)*(b-c)*(b-d)*(c-d) == 0)
	      continue ;
	    mult = 1 ;
	    switch (type)
	      {
	      case 0: /* K-abcd 4-vectors, [ab] [cd] usual */
	      case 1: /* K-abcd 4-vectors, [ab] {cd} should be zero */
		/* case {ab} {cd} vanishes because of the g-mu,nu symmetries */
	      case 2: /* anomaly, epsilon mu,nu,tho,sigma, use STr and fully anti-sym in [abcd] */
		if (a > 3 || b > 3 || c > 3 || d > 3)
		  continue ;
		if (a > b || c > d)
		  continue ;
		break ;
	      case 3: /* K=abij 2-vectors 2-scalars */
		i = c ; j = d ;
		if (a > 3 || b > 3 || i < 4 || j < 4)
		  continue ;
		if (a > b || i > j)
		  continue ;
		break ;
	      case 4: /* K=ijkl 4-scalars */
		i = a ; j = b ; k = c ; l = d ;
		if (i < 4 || j < 4 || k < 4 || l < 4)
		  continue ;
		if (i>j || k > l)
		  continue ;
		if (i<j)
		  mult *= 2 ;
		if (k<l)
		  mult *= 2 ;
		break ;
	      }
	    for (pass = 0 ; pass < 2 ; pass++)
	      {	      
		if (pass == 0)
		  ok = 0;
		if (pass == 1 && ok == 0)
		  continue ;
		
		if (pass == 1)
		  printf ("\n(%d,%d,%d,%d)\t", a, b, c, d) ;
		for (t = 0 ; t < 4 && t < NTYPES ; t++)
		  {
		    float complex z = 0 ;
		    switch (type)
		      {
			
		      case 0: /* K-abcd 4-vectors, [ab] [cd] usual use trace and skew symmetrize in (ab) and in (cd) */
			z = 0 ;
			z += tetraTrace (nchiT[t], t, a, b, c, d) ;
			z -= tetraTrace (nchiT[t], t, a, b, d, c) ;
			z -= tetraTrace (nchiT[t], t, b, a, c, d) ;
			z += tetraTrace (nchiT[t], t, b, a, d, c) ;
			break ;
		      case 1: /* K-abcd 4-vectors, [ab] {cd}  use trace expect zero */
			/* case {ab} {cd} vanishes because of the g-mu,nu symmetries */
			z = 0 ;
			z += tetraTrace (nchiT[t], t, a, b, c, d) ;
			z += tetraTrace (nchiT[t], t, a, b, d, c) ;
			z -= tetraTrace (nchiT[t], t, b, a, c, d) ;
			z -= tetraTrace (nchiT[t], t, b, a, d, c) ;
			break ;
		      case 2: /* anomaly, epsilon mu,nu,tho,sigma, use STr and fully anti-sym in [abcd] */
			z = 0 ;
			z += tetraTrace (nchiS[t], t, a, b, c, d) ;
			z -= tetraTrace (nchiS[t], t, a, b, d, c) ;
			z -= tetraTrace (nchiS[t], t, a, c, b, d) ;
			z += tetraTrace (nchiS[t], t, a, c, d, b) ;
			z += tetraTrace (nchiS[t], t, a, d, b, c) ;
			z -= tetraTrace (nchiS[t], t, a, d, c, b) ;
			break ;
		      case 3: /* K-abij 2-vectors, 2-scalars terme direct {ab}(Lij+Rji) */
			/* if we use STrace everywhere, we get zero on lepton + quarks */
			z = 0 ;
			z += tetraTrace (nchiL[t], t, a, b, i, j) ;
			z += tetraTrace (nchiL[t], t, b, a, i, j) ;
			z += tetraTrace (nchiR[t], t, a, b, j, i) ;
			z += tetraTrace (nchiR[t], t, b, a, j, i) ;
			
			/* K-abij 2-vectors, 2-scalars terme croise Laibj + Rajbi */
			z += -2 * tetraTrace (nchiL[t], t, a, i, b, j) ;
			z += -2 * tetraTrace (nchiR[t], t, a, j, b, i) ;
			break ;
		      case 4: /* K-ijkl, 4 scalars symmetrize in {kl} : L(ikjl + iljk) */
			      /* use Strace => zero, use nchiR == 1/2 Trace  => Higgs potential */
			z = 0 ;
			z += tetraTrace (nchiS[t], t, i, k, j, l) ;
			z += tetraTrace (nchiS[t], t, i, l, j, k) ;
			z += tetraTrace (nchiS[t], t, j, k, i, l) ;
			z += tetraTrace (nchiS[t], t, j, l, i, k) ;
		      }
		    z = mult * z ;
		    zz[t] = z ; /* memorize, to be able to compute the Family e + 3*q */
		    if (pass == 0)
		      {
			if (creal (z * conj(z)) > .1)
			  ok = 1;
		      }
		    else
		      {
			nicePrint ("\t", z) ;
			if (t == 2 || t == 5 || t == 8)  /* compute the family vertex */
			  {
			    nicePrint ("\t", (zz[1] + 3 * zz[2])*.3/.8) ;
			    nicePrint ("\t", -4 * zz[0] + zz[1] + 3 * zz[2]) ;
			  }
		      }
		  }
	      }
	  }
  
  printf ("\n\n") ;
  ac_free (h) ;
} /* mu4p */

/*************************************************************************************************/
/*************************************************************************************************/
/*************************************************************************************************/

int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;

  getCmdLineInt (&argc, argv, "-c3Mask", &c3Mask) ;
  getCmdLineInt (&argc, argv, "-t", &myType) ;
  BOOL SU3 = getCmdLineBool (&argc, argv, "-su3") ;
  int COQ = 0 ;
  getCmdLineInt (&argc, argv, "-coq", &COQ) ;

  int a = 0, b = 0 ;
  getCmdLineInt (&argc, argv, "-a", &a) ;
  getCmdLineInt (&argc, argv, "-b", &b) ;
  if (!COQ && (a || b))
    {
      /* 2021_03_18 
       * construct the 8 matrices for the generic irreps of su(2/1) with h.w. (a,b)
       * verify all commutations relations
       * compute the casimir tensors and operators 

       */
      Kasimirs (a,b) ;
      exit (0) ;
    }
  /* always init, otherwise the gcc linker is unhappy */
  if (1) muInit (h) ;   /* init the 4x4 matrices */
  if (0) muInit2 (h) ;  /* init the 2-families 8x8 rotated matrices */
  if (COQ >= 2) muInitNCoq (a,b, COQ) ;  /* init the 2-families 8x8 coquereaux indecomposable matrices */
  /* verification numerique directe de traces de matrices de pauli */ 
  if (0) muSigma (h) ;
  exit (0) ;
  if(0) feynmanDiagrams () ;


  if (0)   exit (0) ;
 
  /* Verifications des traces sur la theorie des groupes pour l'article sur les anomalies scalaires */
  if (1)
    {
      muConjugate (h) ;
      
      
      printf ("########## Compute the relevant traces of products of 2,3,4 SU(2/1) matrices\n") ;
      printf ("########## In each case, the trace is computed for the neutral representation (N), then for leptons (e), quarks (q) and family (e+3*q)\n") ;
      printf ("########## The observation is that leptons and quarks have anomalous traces, but they compensate each other\n") ;
      printf ("########## The family trace, one lepton +  quarks, is proportional to the neutral trace\n") ;
      printf ("########## In the last column, we check that S = e + 3*q - 4*n == 0\n") ;
      
      printf ("########## Verify the commutators,   all computed norms should vanish\n");
      muStructure () ;
    }
  
  if (0) mu2p ("######### Metric\n# For the even generators (a,b=0123), compute the Super-Trace: STr(ab)\n# For the odd generators (i=4567), compute the Left trace: LTr(ij)\n We hope to find the SU(2/1) Super-Killing metric") ;

  if (1) casimir2 ("######### Casimir 2\n# 1/2 g^AB mu_A mu_B,   we hope to find a diagonal matrix") ;

  if (1) casimir3 ("######### Super Casimir 3\n# 1/6 d^ABC mu_A mu_B mu_C,   we hope to find a diagonal matrix", FALSE) ;
  if (0) casimir3 ("######### Hyper Casimir 3\n# 1/6 d^ABC mu_A mu_B mu_C,   we hope to find a diagonal matrix", TRUE) ;

  exit (0) ;

  if (1) /* pure gauge theory, no fermions, attempt to verify the vector coupling to scalars and tensors */
    {
      if (0) Z2_AA_loopH ("######### Vector propagator, Scalar loop\n") ;
      if (0) Z2_AA_loopB ("######### Vector propagator, Tensor loop\n") ;
      if (0) Z2_AA_loopHB ("######### Vector propagator, Tensor-Scalar loop\n") ;
      if (0) Z2_HH_loopAH ("######### Scalar propagator, Vector-Scalar loop\n") ;
      if (0) Z2_HH_loopAB ("######### Scalar propagator, Vector-Tensor loop\n") ;
      if (0) Z2_BB_loopAB ("######### Tensor propagator, Vector-Tensor loop\n") ;
      if (1) Z2_BB_loopAH ("######### Tensor propagator, Vector-Scalar loop\n") ;
      
      exit (0) ;
      
#ifdef JUNK


      /* 2021_02_26
       * il y a chaquefois 6 triangles, mais le nombre qui depend du new vertex FBA n'est pas le meme pour le diagramme mixte
       * pour respecter Ward vecoriel, on attend une non linrearite
       * avec deux solutions theta-FBPhi    soit theta = 1 comme pour les Fermions loop
       * soit theta = 0, mais peut etre meme que theta = 0 est impossible 
       * ce serait deja un resultat publiable
       * mais cela demande sans doute que le propagateur du vecteurr soit super-Killing si on utilise les d_aij
       * mais le probleme est sans doute le meme avec des couplages f_aij ett Killing  
       */
      if (1) Z3_AHH_loopAAH ("######### Vector-scalar-scalar vertex, Scalar_below-vector-vector loop\n") ;
      if (1) Z3_AHH_loopAAB ("######### Vector-scalar-scalar vertex, Tensor_below-vector-vector loop\n") ;
      if (1) Z3_AHH_loopHHA ("######### Vector-scalar-scalar vertex, Vector_below-scalar-scalar loop\n") ;
      if (1) Z3_AHH_loopBBA ("######### Vector-scalar-scalar vertex, Vector_below-tensor-tensor loop\n") ;
      if (1) Z3_AHH_loopBHA ("######### Vector-scalar-scalar vertex, Vector_below-tensor-scalar loop\n") ;
      if (1) Z3_AHH_loopHBA ("######### Vector-scalar-scalar vertex, Vector_below-scalar-tensor loop\n") ;

      if (1) Z3_ABB_loopAAB ("######### Vector-tensor-tensor vertex, Tensor_below-vector-vector loop\n") ;
      if (1) Z3_ABB_loopAAH ("######### Vector-tensor-tensor vertex, Scalar_below-vector-vector loop\n") ;
      if (1) Z3_ABB_loopBBA ("######### Vector-tensor-tensor vertex, Vector_below-tensor-tensor loop\n") ;
      if (1) Z3_ABB_loopHHA ("######### Vector-tensor-tensor vertex, Vector_below-scalar-scalar loop\n") ;
      if (1) Z3_ABB_loopBHA ("######### Vector-tensor-tensor vertex, Vector_below-tensor-scalar loop\n") ;
      if (1) Z3_ABB_loopHBA ("######### Vector-tensor-tensor vertex, Vector_below-scalar-tensor loop\n") ;

      if (1) Z3_ABH_loopAAB ("######### Vector-tensor-scalar vertex, Tensor_below-vector-vector loop\n") ;
      if (1) Z3_ABH_loopAAH ("######### Vector-tensor-scalar vertex, Scalar_below-vector-vector loop\n") ;
      if (1) Z3_ABH_loopBBA ("######### Vector-tensor-scalar vertex, Vector_below-tensor-tensor loop\n") ;
      if (1) Z3_ABH_loopHHA ("######### Vector-tensor-scalar vertex, Vector_below-scalar-scalar loop\n") ;
      if (1) Z3_ABH_loopBHA ("######### Vector-tensor-scalar vertex, Vector_below-tensor-scalar loop\n") ;
      if (1) Z3_ABH_loopHBA ("######### Vector-tensor-scalar vertex, Vector_below-scalar-tensor loop\n") ;

#endif
      exit (0) ;
    }

  

  if (0) mu3p ("######### Triple Vector Vertex\n# Lie algebra f-abc vertex,\n# compute the trace anti-symmetrized in bc: Tr(a[bc])\n# we hope to find the Lie algebra f-123 = 4i", 0) ;

  if (0) mu3p ("######### Adler-Bardeen Anomalous Triple Vector Vertex\n# d-abc anomalous vertex\n# compute the super-trace symmetrized in bc: STr(a{bc})\n# The anomaly should vanish", 1) ;
  if (0) mu3p ("######### Vector Scalar Vertex\n# since  i and j are oriented, do not symmetrized in i,j but use LTr(aij)-RTr(aji)\n# We hope to find the super-algebra d-aij\n", 2) ;
  if (1) mu3p ("######### Vector Scalar Vertex\n# use Trace (aij - aji), expect zero in f=famille\n", 20) ;
  if (1) mu3p ("######### Vector Scalar Vertex STr measure\n# use SuperTrace (aij - aji), expect zero in f=famille\n", 21) ;
  if (0) mu3p ("######### Vector Scalar Vertex Tr measure\n# use Trace (aij - aji), expect irregularities\n", 22) ;
  if (1) mu3p ("######### Vector Scalar Vertex STr vertex\n# use STrace (aij + aji), expect universal d_aij\n", 23) ;
  if (0) mu3p ("######### The other types of triple vertices, i.e. f-abi and f-ijk should be zero because they do not conserve the even/odd grading\n", 3) ;
  if (0) mu3p ("######### Vector scalar anomaly, Tr (a [ij]) should vanish\n", 4) ;
  exit (0) ;
  


  printf ("\n######### Four vector vertices\n# The 3 types of (abcd) symmetrisations are implied by the trace on the Pauli matrices of the Fermion loop\n") ;
  if (1) mu4p ("#########  K-abcd 4 vectors\n# [ab] [cd]: standard Lie Algebra vertex g_mn f^m_ab f^n_cd", 0) ;
  if (1) mu4p ("#########  K-abcd 4 vectors\n# [ab] {cd} should vanish", 1) ;
  if (1) mu4p ("#########  K-abcd 4 vectors anomaly\n# [abcd]", 2) ;

  printf ("\n######### Two vectors, 2 scalars vertices\n# The scalars are oriented, so we do not symmetrize on (ij)\n") ;
  if (1) mu4p ("#########  K-abij 2 vectors, 2 scalars\n# abij: Symmetize in {ab}, use Lij+Rji\n# Then add the K-aibj Symmetrize in {ab}, use (-2)(L.i.j+R.j.i)", 3) ;

  printf ("\n######### Four scalars\n# The scalars are oriented,{ij} incoming, {kl} outcoming\n") ;
  if (1) mu4p ("#########  K-ijkl Symmetrize in {ij} and {kl}, use Likjl + Liljk", 4) ;
	   
  return 0 ;
}

#ifdef JUNK

Oct 17, 2020
  Dear editor

  Thank you for the proofs
  the copyright notice is fine thank you

  I saw 3 errors, the most important is error C.


  A: Page 4, lines 3 and 4 of paragraph 2
  my error, readers of the arxiv preprint notice that I have a sign mistake in my definition of the matrix gamma_5. Please edit the signs 
  the current equation have
    \psi_R = \farx{1}{4} (1 - \chi) (1 + \gamma_5) and  on the next line \psi_L = \frac{1}{4} (1 - \chi) (1 - \gamma_5)
  please flip the sign in front of \gamma_5 to obtain
    \psi_R = \farx{1}{4} (1 - \chi) (1 - \gamma_5) and  on the next line \psi_L = \frac{1}{4} (1 - \chi) (1 + \gamma_5)



  B: page 10, 3 lines before the end of section 5
    A grammatical typo, please delete the second occurence of the word can
      Please edit   "can the symmetry can be restored"
      into  "can the symmetry be restrored"


  C: An important error in the references generated on your side in some weird way
     I beleive i reported this error previously but it was not fixed
     Reference 24 is random, i never quoted or read that paper
     it should be
A. Connes and J. Lott,
\textit{Particle models and noncommutative geometry,}
Nucl. Phys. B18 (1990), 29-47.


New edition requested october 21

												

 A: Page 4, lines 3 and 4 of paragraph 2
  my error, readers of the arxiv preprint notice that I have a sign mistake in my definition of the matrix gamma_5. Please edit the signs 
  the current equation have
    \psi_R = \farx{1}{4} (1 - \chi) (1 + \gamma_5) and  on the next line \psi_L = \frac{1}{4} (1 - \chi) (1 - \gamma_5)
  please flip the sign in front of \gamma_5 to obtain
    \psi_R = \farx{1}{4} (1 - \chi) (1 - \gamma_5) and  on the next line \psi_L = \frac{1}{4} (1 - \chi) (1 + \gamma_5)


Dear editor
Thank you so much.
The grammatical error page 10 is gone and the reference 24 has been fixed

I am so sorry, there is still a very important sign mistake which is my fault
because you did apply my instructions in my request, but my instructions  were utterly wrong
please edit once more these signs
 
 A: Page 4, lines 3 and 4 of paragraph 2
   the current equation have
    \psi_R = \farx{1}{4} (1 - \chi) (1 - \gamma_5) and  on the next line \psi_L = \frac{1}{4} (1 - \chi) (1 + \gamma_5) \psi_L
  please flip three signs, i wish 
    \psi_R = \farx{1}{4} (1 - \chi) (1 + \gamma_5) and  on the next line \psi_L = \frac{1}{4} (1 + \chi) (1 - \gamma_5)


												
#endif
