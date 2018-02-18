#include "ac.h"
#include "freeout.h"

typedef struct zzStruct { Array aa ; int nmax, s ; float fmax ; AC_HANDLE h; } ZZ ;
typedef struct zzzStruct { int n[100], t ; } ZZZ ;

static void usage (void) ;

/*************************************************************************************/

static int zzGetData (ZZ *zz)
{
  int level, ii = 0, t, nn, s ;
  float fmax = 0, f ;
  ZZZ *up ;

  level = freesetfile (stdin, 0) ;

  while (freecard (level))
    {
      if (freeint (&t))
	{
	  up = arrayp (zz->aa, ii++, ZZZ) ; 
	  up->t = t ;
	  for (s = 2 ; s < zz->s ; s++) 
	    if (freeint (&nn) && nn >= 0 && t > 0)
	      {
		if (nn>=0)
		  {
		    up->n[s-2] = nn ;
		    f = nn ? nn / (float)t : 1/(2.0 * t) ;
		    if (f > fmax) fmax = f ;
		  }
	      }
	    else
	      {
		if (zz->s > 2)
		  fprintf (stderr, "error line %d, awaiting s columns (duration,  %d nb-events) tab separated, got : %s\n"
			   , ii + 1, zz->s - 1, freepos()) ;
		else
		  fprintf (stderr, "error line %d, awaiting s columns (duration,  nb-events) tab separated, got : %s\n"
			   , ii + 1, freepos()) ;
		usage () ;
	      }
	}
    }
  if (!zz->fmax) zz->fmax = fmax ;
  fprintf (stderr, "Found %d patients, maximal frequency fmax = %g\n", ii, fmax) ;

  return ii ;
}

/*************************************************************************************/

static int zzExport (ZZ *zz)
{
  int N, nmax = zz->nmax ;
  int ii, j, n, s ;
  double f, m, q, qt,  total2[zz->s], p[nmax+1], pp[(nmax+zz->s)*(nmax+zz->s)], total = 0 ;
  ZZZ *up ;

  for (j = 0 ; j <= nmax ; j++)
    pp[j] = 0 ;

  N = nmax + 1 ; if (N < zz->s-2) N = zz->s ;
  /* for each patient compute the distribution of the f_j */
  for (s = 2 ; s <= zz->s ; s++)
    {
      for (ii = 0 ; ii < arrayMax (zz->aa) ; ii++)
	{
	  up = arrp (zz->aa, ii, ZZZ) ;
	  qt = 0 ;
	  
	  for (j = 0 ; j <= nmax ; j++)
	    {
	      f = j * zz->fmax / nmax ; 
	      m = f * up->t ;
	      
	      q = exp (-m)   ;  /* i.e. p[0] */
	      for (n = 1 ; n <= up->n[s-2] ; n++)
		q *= m / n ; /* i.e. exp(-m) m ^n / n! , computed recursivelly */
	      p[j] = q ;
	      qt += q ;
	    }
	  /* normalize, so that each patient contributes a total of 1 */
	  for (j = 0 ; j <= nmax ; j++)
	    p[j] /= qt ;
	  /* cumulate to find the distrib of the whole population */
	  for (j = 0 ; j <= nmax ; j++)
	    pp[N*j+s-2] += p[j] ;
	}
    }

  /* normalize */ 
  for (s = 2 ; s <= zz->s ; s++)
    {
      total = 0 ;
      for (j = 0 ; j <= nmax ; j++)
	total += pp[N*j+s-2] ; 
      if (total > 0)
	for (j = 0 ; j <= nmax ; j++)
	  pp[N*j+s-2] /= total ; 
    }
  /* export */
  if (zz->s > 2)
    printf ("Frequency\tProbabilities\n") ;
  else
    printf ("Frequency\tProbability\n") ;

  for (s = 2 ; s <= zz->s ; s++)
    total2[s-2] = 0 ;

  for (j = 0 ; j <= nmax ; j++)
    {
      printf ("%g", j * zz->fmax / nmax) ; 
      for (s = 2 ; s <= zz->s ; s++)
	{
	  printf ("\t%g", pp[N*j+s-2]) ;
	  total2[s-2] += pp[N*j+s-2] ; 
	}
      printf ("\n") ;
    }
  printf ("Total") ;
  for (s = 2 ; s <= zz->s ; s++)
    printf ("\t%g", total2[s-2]) ;
  printf ("\n") ;

  return 0 ;
}

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: laurehisto -n number [-fmax number] [-s number]\n") ;
  fprintf (stderr, "// Example:  laurehisto -n 20 -fmax .3 < data.txt > histo.txt \n") ;  
  fprintf (stderr, "//   -n : number of values in the output, probably between 5 and 100 \n") ;
  fprintf (stderr, "//   -s : last column containing a number of events, default=2 \n") ;
  fprintf (stderr, "//   -fmax : optional, plot for f in range [0 fmax], defaults to max observed f \n") ;
  fprintf (stderr, "//   data.txt : 2 columns, tab delimited,  duration\tnb event\n") ;
  fprintf (stderr, "//   export the Bayes histogram of likely even frequency in the population\n") ;
  fprintf (stderr, "//   assuming a Poisson probability of the eevnts for each patient\n") ;
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  int outlevel = 0 ;
  const char *ici ;
  ZZ zz ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&zz, 0, sizeof (ZZ)) ;
  zz.h = h ;

  if (argc < 2)
    usage () ;

  if (getCmdLineOption (&argc, argv, "-n", &ici))
    {
      if (sscanf (ici, "%d", &(zz.nmax)) == 1 &&
	  zz.nmax > 0) ;
      else
	{
	  fprintf (stderr, "bad -n parameter on command line, should be a positive number\n") ;
	  usage () ;
	}
    }
  zz.s = 2 ;
  if (getCmdLineOption (&argc, argv, "-s", &ici))
    {
      if (sscanf (ici, "%d", &(zz.s)) == 1 &&
	  zz.s > 1) ;
      else
	{
	  fprintf (stderr, "bad -s parameter, should number of last column containong the number of events\n") ;
	  usage () ;
	}
      if (zz.s > 100)
	{
	  fprintf (stderr, "bad -s parameter, sorry, we can only treat 100 columns at a time\n") ;
	  usage () ;
	}
    }
  zz.fmax = 0 ;
  if (getCmdLineOption (&argc, argv, "-fmax", &ici))
    {
      if (sscanf (ici, "%g", &(zz.fmax)) == 1 &&
	  zz.fmax > 0) ;
      else
	{
	  fprintf (stderr, "bad -pmax parameter on command line, should be a positive number\n") ;
	  usage () ;
	}
    }

  zz.aa = arrayHandleCreate (20000, ZZZ, zz.h) ;
  outlevel = freeOutSetFile (stdout) ;	

  fprintf (stderr, "// start: %s\n", timeShowNow()) ;

  zzGetData (&zz) ;
  zzExport (&zz) ;

  fprintf (stderr, "// done: %s\n", timeShowNow()) ;

  if (outlevel)
    freeOutClose (outlevel) ;

  ac_free (zz.h) ;
  if (0) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

