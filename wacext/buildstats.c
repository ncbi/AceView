/* file buildstats.c
 * Compute the stats for each human build
 
 * created sept 14 2004, build 35c
 */
  
#include "../wac/ac.h"
#include "dict.h"

/*********************************************************************/
/* phase 1: count genes */
static void countgenes (AC_DB db)
{
  AC_KEYSET ggs, gcl, tgs, mrnas, ps ;
  int ng, ncl, nbal, ntg, nmrna, np ;
  AC_HANDLE h = ac_new_handle () ;

  ggs = ac_dbquery_keyset (db, "find gene intmap", h) ; 
  ng = ac_keyset_count (ggs) ;

  gcl = ac_dbquery_keyset (db, "find gene cloud_gene && intmap && transcribed_gene", h) ; 
  ncl= ac_keyset_count (gcl) ;

  gcl = ac_dbquery_keyset (db, "find gene balise && intmap && transcribed_gene", h) ; 
  nbal= ac_keyset_count (gcl) ;

  tgs = ac_dbquery_keyset (db, "find tg cdna_clone", h) ; 
  ntg = ac_keyset_count (tgs) ;

  mrnas = ac_dbquery_keyset (db, "find mrna from_gene", h) ; 
  nmrna = ac_keyset_count (mrnas) ;

  ps = ac_dbquery_keyset (db, "find product mrna", h) ; 
  np = ac_keyset_count (ps) ;


  printf ("#######################################################\n") ;
  printf ("Number of genes, balise genes, cloud genes, mrnas and products\n") ;
  printf ("Type\tNumber\n") ;
  printf ("Gene\t%d\n", ng) ;
  printf ("Cloud\t%d\n", ncl) ;
  printf ("Balise\t%d\n", nbal) ;
  printf ("tg\t%d\n", ntg) ;
  printf ("mRNA\t%d\n", nmrna) ;
  printf ("Product\t%d\n", np) ;

 
  printf ("\n") ;

  ac_free (h) ;
} /* countgenes */

/*********************************************************************/
/* phase 1: count multi alignememnts */

static void countali (AC_DB db)
{
  AC_KEYSET nm, nmali, nmmulti, mrna, mrnaali, mrnamulti, est, estali, estmulti ;
  int nnm, nnmali, nnmmulti=0, nmrna, nmrnaali, nmrnamulti=0, nest, nestali, nestmulti=0 ;
  AC_HANDLE h = ac_new_handle () ;

  nm = ac_dbquery_keyset (db, "find est ref_seq && DNA", h) ; 
  nnm = ac_keyset_count (nm) ;
  nmali = ac_dbquery_keyset (db, "find est ref_seq && from_gene && DNA", h) ; 
  nnmali = ac_keyset_count (nmali) ;
  if (1)
    {
      nmmulti = ac_dbquery_keyset (db, "find est ref_seq && DNA && COUNT from_gene > 1", h) ; 
      nnmmulti = ac_keyset_count (nmmulti) ;
    }

  mrna = ac_dbquery_keyset (db, "find est  ( ref_mrna || is_mrna )   && DNA && !ref_seq", h) ; 
  nmrna = ac_keyset_count (mrna) ;
  mrnaali = ac_dbquery_keyset (db, "find est  ( ref_mrna || is_mrna )   && DNA && !ref_seq  && from_gene", h) ; 
  nmrnaali = ac_keyset_count (mrnaali) ;
  if (1)
    {
      mrnamulti = ac_dbquery_keyset (db, "find est  ( ref_mrna || is_mrna )   && DNA && !ref_seq && COUNT from_gene > 1", h) ; 
      nmrnamulti = ac_keyset_count (mrnamulti) ;
    }

  est = ac_dbquery_keyset (db, "find est  ! ( ref_mrna || is_mrna )   && DNA && !ref_seq", h) ; 
  nest = ac_keyset_count (est) ;
  estali = ac_dbquery_keyset (db, "find est ! ( ref_mrna || is_mrna )   && DNA && !ref_seq && from_gene", h) ; 
  nestali = ac_keyset_count (estali) ;
  if (0)
    {
      estmulti = ac_dbquery_keyset (db, "find est ! ( ref_mrna || is_mrna )   && DNA && !ref_seq && COUNT from_gene > 1", h) ; 
      nestmulti = ac_keyset_count (estmulti) ;
    }

  if (!nnm) nnm = 1 ;
  if (!nmrna) nmrna = 1 ;
  if (!nest) nest = 1 ;
  printf ("#######################################################\n") ;
  printf ("Number of clone aligned once or many times\n") ;
  printf ("Type\tTotal\taligned (%%)\tmultiply aligned (%%)\n") ;
  printf ("NM\t%d\t%d (%2.2f%%)\t%d (%2.2f%%)\n", nnm,  nnmali, (float)nnmali*100.0/nnm, nnmmulti, (float)nnmmulti*100.0/nnm) ;
  printf ("mRNA\t%d\t%d (%2.2f%%)\t%d (%2.2f%%)\n", nmrna,  nmrnaali, (float)nmrnaali*100.0/nmrna, nmrnamulti, (float)nmrnamulti*100.0/nmrna) ;
  printf ("EST\t%d\t%d (%2.2f%%)\t%d (%2.2f%%)\n", nest,  nestali, (float)nestali*100.0/nest, nestmulti, (float)nestmulti*100.0/nest) ;
  printf ("\n") ;

  ac_free (h) ;
} /* countali */

/*********************************************************************/
/* phase 3: coding potential */
static void codingPotential (AC_DB db)
{
  AC_HANDLE h ;
  AC_ITER tgs = ac_dbquery_iter (db, "Find gene balise ; > transcribed_gene", 0) ;
  AC_OBJ tg = 0 ;
  AC_KEYSET ks, kspp, ksm, kskk ;
  int ntg2mi = 0, ntg2mii = 0, ntg2mnoi = 0 ;
  int ntg = 0 , ni = 0, nnoi = 0, jj, n00 = 0, n30=0, n100=0, n300=0 ;
  BOOL hasIntron, hasTax, hasNM ;
  int n, nt[16], nm[16] ;
  int  nmrna = 0,  nmrnaWithProduct = 0, nmrnaWithSeveralProduct = 0 ;
  int ni1clone = 0, ni2clone = 0, ni1cloneM = 0, ni2cloneM = 0, ni2mrna = 0, ni2mrnai = 0, niAntisense = 0 ;
  int npp = 0, nppgood = 0, nppbest = 0, npp2m = 0, npp2tg = 0, nppgood2m = 0, nppgood2tg = 0
    , npp1clone = 0, npp2clone = 0, npp3clone = 0, npp2kantor = 0
    , npp1clone2tg = 0, npp2clone2tg = 0, npp3clone2tg = 0
    , npp2m1clone = 0, npp2m2clone = 0, npp2m3clone = 0
    , npp2m1clone2tg = 0, npp2m2clone2tg = 0, npp2m3clone2tg = 0
    ;

  n = 16 ; while (n--) { nt[n] = nm[n]=0 ; }

  while (ac_free (tg) , (tg = ac_next_obj (tgs)))
    {
      h = ac_new_handle () ;
      ntg++ ;
      if (ac_has_tag (tg, "gt_ag") || ac_has_tag (tg, "gc_ag"))
	{ hasIntron = TRUE ; ni++ ; }
      else
	{ hasIntron = FALSE ; nnoi++ ; }

      jj = 0 ;
      if (!jj)
	{
	  ks = ac_objquery_keyset (tg, ">mrna ; Longest_CDS > 900", h) ;
	  if (ac_keyset_count (ks))
	    { n300++ ; jj = 1 ; }
	}
      if (!jj)
	{
	  ks = ac_objquery_keyset (tg, ">mrna ; Longest_CDS > 300", h) ;
	  if (ac_keyset_count (ks))
	    { n100++ ; jj = 2 ; }
	}
      if (!jj)
	{
	  ks = ac_objquery_keyset (tg, ">mrna ; Longest_CDS > 90", h) ;
	  if (ac_keyset_count (ks))
	    { n30++ ; jj = 3 ; }
	}
      if (!jj)
	{ jj = 4 ; n00++ ; }

      ks = ac_objquery_keyset (tg, ">mrna", h) ;
      if (hasIntron)
	ntg2mi += ac_keyset_count (ks) ;
      else
	ntg2mnoi += ac_keyset_count (ks) ;

      if (hasIntron)
	{
	  ks = ac_objquery_keyset (tg, ">mrna ; gt_ag || gc_ag", h) ;
	  ntg2mii += ac_keyset_count (ks) ;
	}

      ks = ac_objquery_keyset (tg, ">mrna ; >product ; Tax_common_ancestor != *sapiens*", h) ;
      if (ac_keyset_count (ks))
	hasTax = 1 ;
      else
	hasTax = FALSE ;

      ks = ac_objquery_keyset (tg, ">cdna_clone IS NM*", h) ;
      if (ac_keyset_count (ks))
	hasNM = 1 ;
      else
	hasNM = FALSE ;

      if (jj)
	{
	  n = 4 * (jj-1) + (hasIntron ? 0 : 2) + (hasTax ? 0 : 1) ;
	  nt[n]++ ;
	  n = 4 * (jj-1) + (hasIntron ? 0 : 2) + (hasNM ? 0 : 1) ;
	  nm[n]++ ;
	}

      ks = ac_objquery_keyset (tg, ">mrna", h) ;
      nmrna += ac_keyset_count (ks) ;

      ks = ac_objquery_keyset (tg, ">mrna ; product", h) ;
      nmrnaWithProduct += ac_keyset_count (ks) ;

      ks = ac_objquery_keyset (tg, ">mrna ; COUNT product > 1", h) ;
      nmrnaWithSeveralProduct += ac_keyset_count (ks) ;

      if (hasIntron)
	{
	  ks = ac_objquery_keyset (tg, "cdna_clone", h) ;
	  ni1clone += ac_keyset_count (ks) ;
	  
	  ks = ac_objquery_keyset (tg, "cdna_clone ; > mrna", h) ;
	  ni1cloneM += ac_keyset_count (ks) ;
	  
	  ks = ac_objquery_keyset (tg, "COUNT cdna_clone > 1", h) ;
	  ni2clone += ac_keyset_count (ks) ;
	  
	  ks = ac_objquery_keyset (tg, "(COUNT cdna_clone > 1) ;  >mrna", h) ;
	  ni2cloneM += ac_keyset_count (ks) ;
	  
	  ks = ac_objquery_keyset (tg, "COUNT mRNA > 1", h) ;
	  ni2mrna += ac_keyset_count (ks) ;
	  
	  ks = ac_objquery_keyset (tg, "COUNT {>mRNA ; gt_ag || gc_ag} > 1", h) ;
	  ni2mrnai += ac_keyset_count (ks) ;
	  
	  ks = ac_objquery_keyset (tg, "COUNT {>antisens_to ; gt_ag || gc_ag } > 0", h) ;
	  niAntisense += ac_keyset_count (ks) ;
	}

      ks = ac_objquery_keyset (tg, ">mrna ; >product ! coding_gap &&  Coding_length >= 300", h) ;
      npp += ac_keyset_count (ks) ;
      
      kspp = ac_ksquery_keyset (ks, "best_product", h) ;
      nppbest += ac_keyset_count (kspp) ;

      kspp = ac_ksquery_keyset (ks, "Tax_common_ancestor != *sapiens*", h) ;
      nppgood += ac_keyset_count (kspp) ;
      
      ksm = ac_ksquery_keyset (kspp, ">mrna", h) ;
      nppgood2m += ac_keyset_count (ks) ;
      if (ac_keyset_count (ks)) nppgood2tg++ ;

      kskk = ac_ksquery_keyset (ks, ">kantor", h) ;
      npp2kantor += ac_keyset_count (kskk) ;

      ksm = ac_ksquery_keyset (ks, ">mrna", h) ;
      npp2m += ac_keyset_count (ks) ;
      if (ac_keyset_count (ks)) npp2tg++ ;

      ks = ac_ksquery_keyset (ks, "Covered_by", h) ;
      npp1clone += ac_keyset_count (ks) ;
      if (ac_keyset_count (ks)) npp1clone2tg++ ; 
      
      ks = ac_ksquery_keyset (ksm, "COUNT CDS_covered_by = 2", h) ;
      npp2clone += ac_keyset_count (ks) ;
      if (ac_keyset_count (ks)) npp2clone2tg++ ; 
      
      ks = ac_ksquery_keyset (ksm, "COUNT CDS_covered_by > 2", h) ;
      npp3clone += ac_keyset_count (ks) ;
      if (ac_keyset_count (ks)) npp3clone2tg++ ; 
      
      
      ks = ac_ksquery_keyset (ksm, "COUNT mrna_covered_by = 1", h) ;
      npp2m1clone += ac_keyset_count (ks) ;
      if (ac_keyset_count (ks)) npp2m1clone2tg++ ; 
      
      ks = ac_ksquery_keyset (ksm, "COUNT mrna_covered_by = 2", h) ;
      npp2m2clone += ac_keyset_count (ks) ;
      if (ac_keyset_count (ks)) npp2m2clone2tg++ ; 
      
      ks = ac_ksquery_keyset (ksm, "COUNT mrna_covered_by > 2", h) ;
      npp2m3clone += ac_keyset_count (ks) ;
      if (ac_keyset_count (ks)) npp2m3clone2tg++ ; 
      
      ac_free (h) ;
    }


  printf ("########################################################\n") ;
  printf ("We considered %d balise genes.\n", ntg) ;
  printf ("Size AA\t   All\tintron  \tConserved\tNM\n") ;
  printf ("  >=300\t%6d\ty:%6d\ty:%6d\ty:%6d\n", n300, nt[0]+nt[1], nt[0], nm[0]) ;
  printf ("  >=300\t      \t        \tn:%6d\tn:%6d\n", nt[1], nm[1]) ;
  printf ("  >=300\t      \tn:%6d\ty:%6d\ty:%6d\n", nt[2]+nt[3], nt[2], nm[2]) ;
  printf ("  >=300\t      \t        \tn:%6d\tn:%6d\n", nt[3], nm[3]) ;
  printf ("\n") ;

  printf ("100/300\t%6d\ty:%6d\ty:%6d\ty:%6d\n", n100, nt[4]+nt[5],  nt[4], nm[4]) ;
  printf ("100/300\t      \t        \tn:%6d\tn:%6d\n", nt[5], nm[5]) ;
  printf ("100/300\t      \tn:%6d\ty:%6d\ty:%6d\n", nt[6]+nt[7], nt[6], nm[6]) ;
  printf ("100/300\t      \t        \tn:%6d\tn:%6d\n", nt[7], nm[7]) ;
  printf ("\n") ;

  printf (" 30/100\t%6d\ty:%6d\ty:%6d\ty:%6d\n", n30, nt[8]+nt[9],  nt[8], nm[8]) ;
  printf (" 30/100\t      \t        \tn:%6d\tn:%6d\n", nt[9], nm[9]) ;
  printf (" 30/100\t      \tn:%6d\ty:%6d\ty:%6d\n", nt[10]+nt[11], nt[10], nm[10]) ;
  printf (" 30/100\t      \t        \tn:%6d\tn:%6d\n", nt[11], nm[11]) ;
  printf ("\n") ;

  printf ("   0/30\t%6d\ty:%6d\ty:%6d\ty:%6d\n", n00, nt[12]+nt[13],  nt[12], nm[12]) ;
  printf ("   0/30\t      \t        \tn:%6d\tn:%6d\n", nt[13], nm[13]) ;
  printf ("   0/30\t      \tn:%6d\ty:%6d\ty:%6d\n", nt[14]+nt[15], nt[14], nm[14]) ;
  printf ("   0/30\t      \t        \tn:%6d\tn:%6d\n", nt[15], nm[15]) ;

 
  printf ("\n") ;
  printf ("%d genes with introns have %d mRNA, i.e. %2.2f per gene\n"
	  , ni, ntg2mi, (float)ntg2mi/(ni+1)) ;
  printf ("%d genes have %d mRNA with introns, i.e. %2.2f per gene\n"
	  , ni, ntg2mii, (float)ntg2mii/(ni+1)) ;
  printf ("%d genes without introns have %d mRNA, i.e. %2.2f per gene\n"
	  , nnoi, ntg2mnoi, (float)ntg2mnoi/(nnoi+1)) ;
  printf ("\n") ;
  printf ("\n") ;

  printf ("Analysis of the %6d balise genes\n", ntg) ;
  printf ("The %6d balise genes have: %6d mRNAs,\n", ntg, nmrna) ;
  printf ("                                 %6d mRNAs with product,\n",  nmrnaWithProduct) ;
  printf ("                                 %6d mRNAs with several products.\n", nmrnaWithSeveralProduct) ;
  printf ("They contain %6d products with 100 AA, %6d best, %6d kantor \n", npp, nppbest, npp2kantor) ;
  printf ("             %6d of these products from %d mRNAs for %d genes are conserved\n", nppgood, nppgood2m, nppgood2tg) ;
  printf ("             %6d of these products from %d genes are covered by a single read\n", npp1clone, npp1clone2tg) ;
  printf ("             %6d of the best products from %d genes are covered by 2 clones\n", npp2clone, npp2clone2tg) ;
  printf ("             %6d of the best products from %d genes are covered by more than 2 clones\n", npp3clone, npp3clone2tg) ;
  printf ("They contain %6d mRNAs from %d genes corresponding to these %d products\n", npp2m, npp2tg, npp) ;
  printf ("             %6d of these mRNAs from %d genes are covered by a single clone\n", npp2m1clone, npp2m1clone2tg) ;
  printf ("             %6d of these mRNAs from %d genes are covered by 2 clones\n", npp2m2clone, npp2m2clone2tg) ;
  printf ("             %6d of these mRNAs from %d genes are covered by more than 2 clones\n", npp2m3clone, npp2m3clone2tg) ;
  printf ("\n") ;
  printf ("\n") ;

  printf ("Analysis of the %6d genes with standard introns\n", ni) ;
  printf ("  %6d with >= 1 clone lead to %6d mRNAs\n" , ni1clone, ni1cloneM) ;
  printf ("  %6d with >= 2 clone lead to %6d mRNAs\n" , ni2clone, ni2cloneM) ;
  printf ("    %6d have at least 2 mRNAs\n" , ni2mrna) ;
  printf ("    %6d have at least 2 mRNAs with introns\n" , ni2mrnai) ;
  printf ("  %6d are antisense of a gene with standard introns\n", niAntisense) ;
  printf ("\n") ;
  printf ("\n") ;

  ac_free (tgs) ;
} /* codingPotential */

/*********************************************************************/
/* phase 3: coding potential */
static void maxim (AC_DB db)
{
  AC_HANDLE h ;
  AC_ITER tgs = ac_dbquery_iter (db, "Find transcribed_gene gene && read", 0) ;
  AC_ITER rs ;
  AC_OBJ tg = 0, rr = 0 ;
  DICT *dict = 0 ;
  const char *ucp ;
  int n, nn[32], bestn, maxn, ntg = 0, nmaxim = 0 ;

  while (ac_free (tg) , (tg = ac_next_obj (tgs)))
    {
      h = ac_new_handle () ;
      ntg++ ;
      n = 32 ; while (n--) { nn[n] = 0 ; }
      dict = dictHandleCreate (32, h) ;
      rs = ac_objquery_iter (tg, ">read genecard_id", h) ;
      while (ac_free (rr) , (rr = ac_next_obj (rs)))
	{
	  ucp = ac_tag_printable (rr, "genecard_id", 0) ;
	  if (!ucp)
	    continue ;
	  dictAdd (dict, ucp, &n) ;
	  if (n < 32) nn[n]++ ;
	}
      maxn = -1 ; bestn = 0 ;
      for (n = 0 ; n < 32 ; n++)
	if (nn[n] > maxn)
	  { maxn = nn[n] ; bestn = n ; }
      if (maxn > 0)
	{
	  nmaxim++ ;
	  ucp = ac_tag_printable (tg, "Gene", 0) ;
	  if (ucp)
	    {
	      printf ("Gene \"%s\"\n", ucp) ;
	      ucp = dictName (dict, bestn) ;
	      if (ucp)
		printf ("Genecard_id \"%s\"\n\n", ucp) ;
	    }
	}
      ac_free (h) ;
    }


  printf ("########################################################\n") ;
  printf ("We considered %d genes, attributed %d maxims\n", ntg, nmaxim) ;

  ac_free (tgs) ;
} /* maxim */

/*********************************************************************/

static void usage (void)
{
  printf ("Usage:  buildstats a:machine:port test1 [test2 [test3...]]\n") ;
  printf ("   or: sbuildstats $ACEDB test\n") ;
  printf  (" where test = 1: countgenes\n"
	   "              2: countali\n"
	   "              3: codind potential\n"
	   "              4: maxim-genecard\n"
	   "             99: all but maxim\n\n") ;
  exit (1) ;
}

/*********************************************************************/

int main (int argc, char *argv[])
{
  int ir ;
  const char *err = 0 ;
  char *target = "" ;
   AC_DB db = 0 ;

  if (argc > 2)
    target = argv[1] ;
  else
    usage () ;

  db = ac_open_db (target , &err) ;
  if (err)
    printf ("Error message %s", err) ;
  if (!db)
    messcrash ("could not open %s", target) ;

  for (ir = 2 ; ir < argc ; ir++)
    {
      int j = atoi (argv[ir]) ;
      switch (j)
	{
	case 1:  countgenes (db) ; break ;
	case 2:  countali (db) ; break ;
	case 3:  codingPotential (db) ; break ;
	case 4:  maxim (db) ; break ;
	case 99:
	  countgenes (db) ;
	  countali (db) ;
	  codingPotential (db) ;
	  break ;
	default: usage () ;
	  break ;
	}
    }

  ac_db_close (db) ;
  
  return 0 ;
}

/*********************************************************************/
/*********************************************************************/
