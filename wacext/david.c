#include "ac.h"


/* 
   This program does a simulation on protein distributions
   
   We suppose the existence on N famillies i=0.1.2.3,.. N of proteins
   Each has a proba p[i] of beeing deleterious
   Each is produced in k[i] copies
     K = sum ( k[i] ) the total number of protein produced is assumed fixed: 
     F = sum ( k[i] p[i] ) the total number of bad proteins is assumed bounded F < FMAX
     S = K! / prod (k[i]!) the "entropy" of the produced proteins is assumed maximal

   Question:
     What is the distribution of k[] and p[]
     We hope and expect a log(p) = k (log (k))
*/

static const char magic[] = "magic" ;

typedef struct ddStruct { 
  const char *magic ;
  int nIter ;          /* required number of iterations */
  int max ;            /* number of types of proteins */
  int *k ;             /* number of protein in each type */
  double *p ;          /* probability of misfolding in percent */
  int K ;              /* sum k[i] = total number of proteins */
  double F ;              /* sum (k[i] p[i]) = 100*total number of misfolded proteins */
  double fMax ;
  int t0, t1, t2 ;
  double S, S0, S1, S2 ;       /* log of the entropy */
 } DD ;

static void ddShow (DD *dd) ;

/***************************************************************/
/* utilities */
static void ddDestroy (DD *dd)
{
  if (dd && dd->magic == magic)
    {
      dd->magic = 0 ;
      messfree (dd->k) ;
      messfree (dd->p) ;
      messfree (dd) ;
    }
} /* ddDestroy  */

/***************************************************************/

static DD* ddCreate (int max)
{
  DD *dd = 0 ;
  if (max > 0)
    {
      dd = (DD *) messalloc (sizeof (DD)) ;
      dd->magic = magic ;
      dd->max = max ;
      dd->k = (int *) messalloc (max * sizeof (int)) ;
      dd->p = (double *) messalloc (max * sizeof (double)) ;
    }
  return dd ;
} /* ddCreate */

/***************************************************************/
/******************** actual work ******************************/
/***************************************************************/
/* compute the log of entropy of the config */
static double ddS (DD *dd)
{
  int i, K ;
  double S ;

  for (i = K = 0 ; i < dd->max ; i++)
    K += dd->k[i] ;
  S = K * log((double) K) ;
  for (i = 0 ; i < dd->max ; i++)
    S -= (dd->k[i] > 0 ? dd->k[i] * log((double) dd->k[i]) : 0) ;

  return S ;
} /* ddS */

/***************************************************************/
/* compute the baddies */
static double ddF (DD *dd)
{
  int i;
  double F = 0 ;

  for (i = 0 ; i < dd->max ; i++)
    F += dd->p[i] * dd->k[i] ;

  return F ;
} /* ddS */

/***************************************************************/

static void ddInit (DD *dd)
{
  int i, j = (dd->max -1)/2;
  int k = dd->K/dd->max - j ;
  double x = 2 ;

  for (i = j = 0 ; i < dd->max ;i++)
    dd->k[i] = k + dd->max - i - 1 ;
  for (i = 0, j = dd->K ; i < dd->max ;i++)
    j -= dd->k[i] ;
  if (j > 0 )
    for (i = dd->max - 1 ; i >= 0 && j-- > 0 ; i--)
      dd->k[i]++ ;
  else if (j < 0 )
    for (i = dd->max - 1 ; i >= 0 && j++ < 0 ; i--)
      dd->k[i]-- ;
  for (i = 0, j = 0 ; i < dd->max ;i++)
    j += dd->k[i] ;

  fprintf (stderr, "Init: sum = %d   should be %d\n", j, dd->K) ;
  
  for (i = 0 ; i < dd->max ; x *= 1.3, i++)
    dd->p[i] = x ;
  dd->F = ddF (dd) ;
  dd->fMax = dd->F/1.7 ;
  dd->S = ddS (dd) ;
  dd->S0 = dd->S ;
  return  ;
} /* ddInit */

/***************************************************************/

static void ddIterate (DD *dd)
{
  int nIter = dd->nIter ; /* number of iteration */
  int s, delta, knew, K, K0 = dd->K ;
  int iter, ii, i, j, max = dd->max, *k = dd->k  ;
  double S, oldS ;
  int oldk[dd->max] ;
  double oldF = ddF (dd) ;
  double F, fMax = dd->fMax ;

  oldS = dd->S = ddS (dd) ;
  memcpy (oldk, k, max*sizeof(int)) ;

  for (iter = 0 ; iter < nIter ; iter++)
    {
      for (i = K = 0 ; i < max ; i++)
	K += k[i] ;

      /* change randomly one of the value by +- 1% or +- 1 */
      j = random () ;                     /* random integer */
      ii = ((j % max) + max) % max ;       /* alway positive */
      s = random () & 0x1 ? 1 : -1 ;      /* random sign */
      if (iter < nIter/2)
	delta = k[ii] > 30 ? k[ii]/30 : 1 ;
      else
	delta = k[ii] > 100 ? k[ii]/100 : 1 ;
      knew  = k[ii] + s * delta ;
      if (knew < 0) knew = 0 ;
      if (knew > K0) knew = K0 ;
      delta = knew - k[ii] ;
      if (delta < 0) { delta = -delta ; s = 1 ; }
      else { s = -1 ; }
      k[ii] = knew ; 
      /* change randomly the other values so that the sum K is fixed */
      j = (randint () % max) | 0x1 ;
      i = randint () % max ;
      for (; delta > 0 ;)
	{
	  i += j ;
	  i %= max ;
	  if (i != ii)
	    { delta-- ; k[i] += s ; }
	}
      
      if (ii > 0 && k[ii] > k[ii-1])
	{ i = k[ii-1] ; k[ii-1] = k[ii] ; k[ii] = i ; }
      if (ii < max - 1 && k[ii] < k[ii+1])
	{ i = k[ii+1] ; k[ii+1] = k[ii] ; k[ii] = i ; }

      for (i = K = 0 ; i < max ; i++)
	K += k[i] ;

      F = ddF(dd) ;
      S = ddS (dd) ;
      if ((oldF > fMax && F < oldF) || ( F <= fMax && S > oldS)) /* accept */
	{
	  if (!dd->S1 && oldF <= fMax) { dd->S1 = S ; dd->t1 = iter ; }
	  if (S > oldS) { dd->S2 = S ; dd->t2 = iter ; }
	  oldF = dd->F = F ;
	  oldS = dd->S = S ;
	  memcpy (oldk, k, max*sizeof(int)) ;
	}
      else /* backtrack */
	memcpy (k, oldk, max*sizeof(int)) ;
      if (0) ddShow (dd) ;
    }
  return  ;
} /* ddIterate */

/***************************************************************/

static void ddShow (DD *dd)
{
  int i, K ;
     
  for (i = K = 0 ; i < dd->max ; i++)
    K += dd->k[i] ;

  printf ("Results\t K = %d\t F = %d, fMax = %d\nS(0) = %g, S1(%d) = %g, S(%d) = %g, S(%d)=%g\n"
	  , K, (int)dd->F,  (int)dd->fMax, dd->S0, dd->t1, dd->S1, dd->t2, dd->S2, dd->nIter, dd->S) ;
  printf ("Type\tk\tp\tlog(k)\tlog(p)\n") ;

  for (i = 0 ; i < dd->max ; i++)
    printf ("%d\t%d\t%3.2f\t%3.2f\t%3.2f\n",
	    i+1, dd->k[i], dd->p[i], log (dd->k[i]), log (dd->p[i])) ;
  return  ;
} /* ddShow */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char *err)
{
  fprintf (stderr, "// Usage: david -n 100 \n") ;
  fprintf (stderr, "// Example:  david -n 100 \n") ;  
  fprintf (stderr, "//    export the distrib of proteins famillies\n") ;
  fprintf (stderr, "\n// ERROR %s\n", err) ;
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  int nn = 100 ;
  int nIter = 300000 ;
  DD *dd = 0 ;
  const char *ccp ;

  freeinit () ; 
 
  /* mandatory arguments */
  if (getCmdLineOption (&argc, argv, "-n", &ccp))
    {
      if (sscanf (ccp, "%d", &nn) == 1 &&
	  nn >0) ;
      else
	usage (messprintf ("Bad argument -n %s", ccp)) ;
    }
   if (getCmdLineOption (&argc, argv, "-nIter", &ccp))
    {
      if (sscanf (ccp, "%d", &nIter) == 1 &&
	  nIter >0) ;
      else
	usage (messprintf ("Bad argument -n %s", ccp)) ;
    }
 
    
  dd = ddCreate (nn) ;
  dd->K = 1000 * dd->max * (dd->max -1) ; dd->K /= 2 ;
  dd->nIter = nIter ;
  ddInit (dd) ; 
  if (0)ddShow (dd) ;
  ddIterate (dd) ;
  ddShow (dd) ;
  ddDestroy (dd) ;
  fprintf (stderr, "// done\n" ) ;
  if (0) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

