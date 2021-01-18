#include    "../wac/ac.h"
#include    "vtxt.h"
#include "freeout.h"


int firstDivisor (int nn)
{
#define NN 10000              /* max size of prime table */
  static int primeList [NN] ; /* lazy table of known primes */
  static int maxKnown = 0 ;
  int i, p, q ;               /* running integers */
  int testing = 1 ;           /* 1: show prime list when new number computed */

  if (!maxKnown) /* initialise */
    {
      memset (primeList, 0, NN * sizeof (int)) ;
      maxKnown = 1 ;
      primeList [0] = 2 ;
    }
  while (1) /* lazy fill of primeList */
    /* NB if p is our largest known prime, there
       is always another prime q < p*p
    */
    {
      p = primeList [maxKnown - 1] ;
      if (p*p > nn) 
	break ;
      for (q = p + 1 ; firstDivisor (q) ; q++) ;
      if (testing)
	printf ("%6d is prime\n", q) ;
      primeList [maxKnown++] = q ;
      if (maxKnown == NN)
	{ 
	  printf ("Overflow in firstDivisor nn=%d, maxKnown=%d q = %d\n", 
		  nn, maxKnown, q) ;
	  exit (1) ;
	}
    }
  /* testing */
  for (i = 0 ; i < maxKnown ; i++)
    {
      p = primeList [i] ;
      if (p*p > nn) 
	return 0 ; /* nn is prime, I return 0 for ease of use */
      if ( nn == p * (nn/p))
	return p ; /* not a prime */
    }
  return 0 ; /* nn is prime, I return 0 for ease of use */
}

int main (int argc, char *argv[])
{
  char *usage = "Usage: asndump a:annie:2000101 chrom_name out_file_name\n" ;
  char *err = 0, *chrom ;
  char *target = "a:annie:2000101" ;
  char *fNam = 0 ;
  int p, k, level = 0 ;
  FILE *ff=0  ;
  AC_DB db = 0 ;
  AC_OBJ oMap = 0, oGene = 0 ;
  vTXT blkp = vtxtCreate () ;

  k = 8 ; p = firstDivisor (k) ;
  printf ( "first divisor of %d is %d\n\n", k, p) ;

  k = 301 ; p = firstDivisor (k) ;
  printf ( "first divisor of %d is %d\n\n", k, p) ;
  k = 1027 ; p = firstDivisor (k) ;
  printf ( "first divisor of %d is %d\n\n", k, p) ;
  exit (0) ;

  if (argc > 1)
    target = argv[1] ;
  if (argc > 2)
    chrom = argv[2] ;
  else
    messcrash (usage) ;
 
  if (argc > 3)
    fNam = argv[3] ;
  db = ac_open_db (target , &err) ;

  freeinit () ;
  freeOutInit() ;


  if (fNam)
    {
      ff = filopen (fNam,0, "w") ;
      if (!ff)
	messcrash ("cannot open output file %s", fNam) ;
    }
  level = freeOutSetFile (ff ? ff : stdout) ;


  if (err)
    printf ("Error message %s", err) ;
  if (!db)
    messcrash ("could not open %s", target) ;

#ifdef 0
  /* the existence of this code forces function iin sworm.c not to be static 
   * if we need this code again, it should eb inside swrm.c
   */

  gmpJumpPointInit () ;
  
  if ((oMap = ac_get_obj (db, "Map", chrom, 0)))
    {
      ficheChromosomeDump (blkp, db, oMap, 'r')  ;    
    }
  else if ((oGene = ac_get_obj (db, "Gene", chrom, 0)))
    ficheGene (blkp, db, oGene, 'x') ;
  else if ((oGene = ac_get_obj (db, "mrna", chrom, 0)))
    {
      ficheMrna (blkp, db, oGene, 'x') ;
      /*
      AC_HANDLE h = handleCreate () ;
      char *dna = ac_obj_dna (oGene, h) ;
      int ln = strlen (dna) ;
      printf("dna of %s length ln\n", ac_name(oGene), ln) ;
      ac_free (h) ;
      */
    }
  else if ((oGene = ac_get_obj (db, "paper", "p1",0)))
    {
      AC_HANDLE h = handleCreate () ;
      AC_OBJ aut ;
      
      aut = ac_tag_obj (oGene, "Title", h) ;

      printf ("paper %s title %s\n", ac_name(oGene), ac_name(aut)) ;
      ac_free (h) ;
    }
  else
    printf ("cannot find map %s\n", chrom) ;
#endif

  if (vtxtPtr (blkp)) 
    {
      if (0) freeOut (vtxtPtr (blkp)) ;
      if (1 && ff) fprintf (ff, "%s", vtxtPtr (blkp)) ;
    }
 laba:
  vtxtDestroy (blkp) ;
  gmpJumpPointDestroy () ;
      

  freeOutClose (level) ;
  filclose (ff) ;

  ac_free (oMap) ;
  ac_free (oGene) ;
  ac_db_close (db) ;

  freeshutdown () ;
  freeOutShutDown () ;
  return 0 ;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

