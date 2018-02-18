#include "../wac/ac.h"
#include "dict.h"
#include <ctype.h>

#define isIndexed(_cc) (isalnum(_cc)||((_cc)=='-'))

typedef struct absStruct { AC_OBJ key ; DICT *dict ; char *text ; int nw ; } ABS ;

/******************************************************************/

static int tokenize (ABS *aa)
{
  int nw = 0 ;
  char cc, *cp = aa->text, *cq ;
 
  while (*cp)
    {
      cq = cp ;
      while (isIndexed ((int)*cq)) cq++ ;
      cc = *cq ; *cq = 0 ;
      if (cq > cp &&
	  dictAdd (aa->dict, cp, 0))
	nw++ ; 

      *cq = cc ;
      if (cp < cq)
	cp = cq ;
      else
	cp++ ;
    }
  return nw ;
} /* tokenize */

/******************************************************************/

static void addAbstract (Array abs, AC_OBJ key)
{
  int nn = arrayMax (abs) ;
  char *cp ;

  ABS *aa = arrayp (abs, nn, ABS) ;

  if ((cp = ac_longtext (key, 0)))
    {
      aa = arrayp (abs, nn, ABS) ;
      aa->key = key ;
      aa->text = cp ;
      aa->dict = dictCreate (256) ;
      dictAdd (aa->dict, ac_name(key), 0) ;
      
      aa->nw = tokenize (aa) ;
      if (nn > 0 && !aa->nw)
	arrayMax(abs) = nn -1 ;
      if (0) printf ("***%s   %d words\n%s\n\n", ac_name(key), aa->nw, aa->text) ;
    }
  return ;
} /* addAbstract */

/******************************************************************/

static int nCommonWords (Array abs, int ii, int jj)
{
  int nw = 0, i ;
  ABS *aa, *bb ;

  aa = arrp (abs, ii, ABS) ;
  bb = arrp (abs, jj, ABS) ;
  for (i = 1 ; i <= dictMax (aa->dict) ; i++)
    if (dictFind (bb->dict, dictName (aa->dict, i), 0))
      nw++ ;
  return nw ;
} /* nCommonWords */

/******************************************************************/

static BOOL absCompare (Array abs, int ii) 
{
   ABS *aa, *bb ;
   int jj ;
  int n11, n12, n21, n22, h1, h2, v1, v2, ntot ;
  double chi2 = 0, th11, th12, th21, th22, ztot ;

  aa = arrp (abs, ii, ABS) ;
  for (jj = ii + 1 ; jj < arrayMax (abs) ; jj++)
    {
      bb = arrp (abs, jj, ABS) ;
      n12 = n21 = nCommonWords (abs, ii, jj) ;
     
      h1 = arr (abs, ii, ABS).nw ;
      n11 = h1 - n12 ;
      h2 = arr (abs, jj, ABS).nw ;
      n22 = h2 - n21 ;
      ntot = h1 + h2 ;
      v1 = n11 + n21 ; v2 = n12 + n22 ;
      if (h1>0&&h2>0&&v1>0&&v2>0)
	{
	  ztot = ntot ;
	  th11 = h1*v1/ztot ; th12 = h1*v2/ztot ;
	  th21 = h2*v1/ztot ; th22 = h2*v2/ztot ;
	  chi2 =
	    (n11 - th11)*(n11 - th11)/th11 + 
	    (n12 - th12)*(n12 - th11)/th12 +
	    (n21 - th21)*(n21 - th21)/th21 +
	    (n22 - th22)*(n22 - th22)/th22 ;
	  chi2 = n11 * n11 /th11 + n22 *n22/th22 ;
	}
      if (chi2 < 1.96)
	{
	  printf ("LongText %s //  %d words\n%s\n\n", dictName (aa->dict, 0), aa->nw, aa->text) ;
	  printf ("LongText %s //  %d words\n%s\n\n", dictName (bb->dict, 0), bb->nw, bb->text) ;
	  
	  printf ("// %6d%6d%6d\n// %6d%6d%6d\n// %6d%6d%6d\t%8.2f\n\n",
		  n11, n12, h1, n21, n22, h2, v1, v2, ntot, chi2) ;
	}
    }
  return chi2 < 2.96 ? TRUE : FALSE ;
} /* absCompare */

/******************************************************************/

int main (int argc, char *argv[])
{
  const char *err = 0 ;
  char *target = "a:annie:2000101" ;
  AC_DB db = 0 ;
  AC_ITER iter ;
  AC_OBJ k ;
  int ii, nn = 0 ;

  Array abs = 0 ;

  if (argc > 1)
    target = argv[1] ;
  db = ac_open_db (target , &err) ;
  if (err)
    printf ("Error message %s", err) ;
  if (!db)
    messcrash ("could not open %s", target) ;

  abs = arrayCreate (10000, ABS) ;
  iter = ac_query_iter (db, 0, "find paper ; > abstract", 0, 0) ;
  while (--nn && (k = ac_next_obj (iter)))
    {
      addAbstract (abs, k) ;
      ac_free (k) ;  
      /* if (ok) break ; */
    }
  ac_free (iter) ;

  for (ii = 0 ; ii < arrayMax (abs) ; ii++)
    absCompare (abs, ii) ;
      

  ac_db_close (db) ;

  return 0 ;
} /* main */

/******************************************************************/
/******************************************************************/
