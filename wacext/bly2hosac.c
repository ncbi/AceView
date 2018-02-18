#include "../wac/ac.h"
#include <errno.h>

static int bly2hosac (AC_DB db)
{
  AC_ITER iter = ac_query_iter (db, 0, "Find cdna_clone COUNT from_gene = 1 && COUNT clone_group = 1", 0, 0) ;
  int nn, nSame, nBlyIncludesHosac, nHosacIncludesBly, nNeither, nTwiceUnmapped, nUmInBlywithMappedHosac, nOriginal, nHosacOnly ;
  AC_KEYSET ks1, ks2, ks4, ks31, ks32 ;
  int n1, n2, n31, n32, n4 ;
  AC_OBJ clone ;

  nn = nSame = nBlyIncludesHosac = nHosacIncludesBly = nNeither = nTwiceUnmapped = nUmInBlywithMappedHosac = nOriginal= nHosacOnly = 0 ;
  while ((clone = ac_next_obj (iter)))
    {
      nn++ ;
      ks1 = ks2 = ks4 = ks31 = ks32 = 0 ;
      ks1 = ac_dbquery_keyset (db, messprintf("find cdna_clone %s ; >from_gene ; >cdna_clone", ac_name(clone)), 0) ;
      n1 = ac_keyset_count (ks1) ;
      ks2 = ac_dbquery_keyset (db, messprintf("find cdna_clone %s ; >clone_group ; ! IS *um* ; >contains", ac_name(clone)), 0) ;
      n2 = ac_keyset_count (ks2) ;
      ks31 = ac_copy_keyset (ks1, 0) ;
      ks32 = ac_copy_keyset (ks2, 0) ;
      n31 = ac_keyset_minus (ks31, ks2) ;
      n32 = ac_keyset_minus (ks32, ks1) ;
      if (n1 && n2)
	{
	  if (!n31 && !n32) nSame++ ;
	  else if (n1 && !n31) nHosacIncludesBly++ ;
	  else if (n2 && !n32) nBlyIncludesHosac++ ;
	  else if (n31 && n32) nNeither++ ;
	}
      else if (n1 && !n2) /* mapped only by bly */
	{
	  ks4 = ac_ksquery_keyset (ks1, ">clone_group ! IS *um*", 0) ;
	  n4 = ac_keyset_count (ks4) ;
	  if (n4) nUmInBlywithMappedHosac++ ;
	  else nOriginal++ ;
	}
      else if (!n1 && n2) /* mapped only by hosac */
	{
	  nHosacOnly++ ;
	}
      else if (!n1 && !n2)
	{
	  nTwiceUnmapped++ ;
	}
      ac_free (ks1) ;
      ac_free (ks2) ;
      ac_free (ks31) ;
      ac_free (ks32) ;
      ac_free (ks4) ;
    }
  ac_free (iter) ;

  printf ("total number of clones %d\n", nn) ;
  printf ("nSame %d\n", nSame) ;
  printf ("nHosacIncludesBly %d\n", nHosacIncludesBly) ;
  printf ("nBlyIncludesHosac %d\n", nBlyIncludesHosac) ;
  printf ("nNeither %d\n", nNeither) ;
  printf ("nUmInBlywithMappedHosac %d\n", nUmInBlywithMappedHosac) ;
  printf ("nOriginal %d\n", nOriginal) ;
  printf ("nHosacOnly %d\n", nHosacOnly) ;
  printf ("nTwiceUnmapped %d\n", nTwiceUnmapped) ;

  return nn ;
}


int main (int argc, char **argv)
{
  const char *s = "ok" ;
  char *dbName = argc>=2 ? argv[1] : 0 ;
  AC_DB db = dbName ? ac_open_db (dbName, &s) : 0 ;
  if (!dbName)
    {
      printf ("Usage: bly2hosac a:vesta:12345\n compare hosac and bly clusters\n") ;
      exit (1) ;
    }
  if (!db)
    messcrash ("Failed to open db %s, error %s", dbName, s) ;

  bly2hosac (db) ;

  ac_db_close (db) ;
  return 0 ;
}

