/*
 * Given a database containing a number of objects positioned on
 * the genome using the schema
 *     Location ?any_class UNIQUE Int UNIQUE Int
 * cluster all these objects if they overlap, using model
 *     ?Cluster Child class_name object XREF Cluster
 *
 * Program written for Sam Cartinhour, by Jean Thierry-Mieg
 * created 2007_04_25 mieg, mieg@ncbi.nlm.nih.gov
 */


#include "../wac/ac.h"
#include "keyset.h"
#include <errno.h>
#include <math.h>
#include "bitset.h"
#include "keyset.h"
#include "freeout.h"
#include "vtxt.h"
#include "dict.h"
#include "aceio.h"

BITSET_MAKE_BITFIELD 	 /*	 define bitField for bitset.h ops */

typedef struct zoneStruct { AC_DB db ; Array aa ; const char *classes, *tag, *prefix ; AC_HANDLE h ; BOOL mixStrands, frameSensitive ; float dx ;} ZZ ;
typedef struct zzoneStruct { int strand, cluster ; float a1, a2 ; KEY gene, map ; } ZZZ ;

static void zzExport (ZZ *zz) ;


/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static void zzGetData (ZZ *zz)
{
  int nn = 0, ir ;
  float a1, a2 ;
  char *cp0 = zz->classes ? strnew (zz->classes, 0) : 0 ;
  char cc, *cp = cp0, *cq ;
  char *err ;
  ZZZ *up ;
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  AC_KEYSET ks = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;

  while (cp)
    {
      cq = strstr (cp, "|") ;
      if (cq)
	{
	  cc = *cq ;
	  *cq = 0 ;
	}
      else
	cc = 0 ;
      if (nn++) vtxtPrintf (txt , " SETOR ") ;
      vtxtPrintf (txt , "{Find %s %s}", cp, zz->tag) ;
      if (cc) cp = cq + 1;
      else cp = 0 ;
    }
  
  if (vtxtPtr (txt) && (ks = ac_dbquery_keyset (zz->db, vtxtPtr (txt),  0)))
    {
      if (1)
	{
	  cp = hprintf (h, "select g,m,m1,m2 from g in @, m in g->%s, m1 in m[1], m2 in m1[1]"
			, zz->tag) ;
	  tbl = ac_bql_table (zz->db, cp, ks, "2+3+4", &err, h) ;
	}
      for (ir = 0 ; ir < tbl->rows ; ir++)
	{
	  up = arrayp (zz->aa, ir, ZZZ) ;
	  up->gene = ac_table_key (tbl, ir, 0, 0) ;
	  up->map = ac_table_key (tbl, ir, 1, 0) ;
	  if (zz->dx)
	    {
	      a1 = ac_table_float (tbl, ir, 2, 0) ;
	      a2 = ac_table_float (tbl, ir, 3, a1) ;
	    }
	  else
	    {
	      a1 = ac_table_int (tbl, ir, 2, 0) ;
	      a2 = ac_table_int (tbl, ir, 3, a1) ;
	    }


	  if (a1 > a2)
	    { up->strand = -1 ; up->a1 = a2 ; up->a2 = a1 ; }
	  else
	    { up->strand = 1 ; up->a1 = a1 ; up->a2 = a2 ; }
	  if (zz->mixStrands)
	    up->strand = 0 ;
	}
      }
  ac_free (h) ;
} /* zzGetData */

/*************************************************************************************/

static int zzOrder (const void *a, const void *b)
{
  const ZZZ *up = (const ZZZ*)a, *vp = (const ZZZ *)b ;
  int n ;
  float da ;

  n = up->map - vp->map ;
  if (n)
    return n ;

  n = up->strand -  vp->strand ;
  if (n)
    return n ;

  da = up->a1 - vp->a1 ;
  if (da < 0)
    return -1 ;
  if (da > 0)
    return -1 ;
  da = up->a2 - vp->a2 ;
  if (da < 0)
    return -1 ;
  if (da > 0)
    return -1 ;

  n = up->gene - vp->gene ;
  return n ;
} /* zzOrder */

/*************************************************************************************/

static int zzOrderByCluster (const void *a, const void *b)
{
  const ZZZ *up = (const ZZZ*)a, *vp = (const ZZZ *)b ;
  int n ;
  float da ;

  n = up->cluster - vp->cluster ;
  if (n)
    return n ;

  n = up->strand -  vp->strand ;
  if (n)
    return n ;

  da = up->a1 - vp->a1 ;
  if (da < 0)
    return -1 ;
  if (da > 0)
    return -1 ;
  da = up->a2 - vp->a2 ;
  if (da < 0)
    return -1 ;
  if (da > 0)
    return -1 ;

  n = up->gene - vp->gene ;
  return n ;
} /* zzOrderByCluster */

/*************************************************************************************/

static void zzMerge (ZZ *zz)
{
  int ii, kk, cluster = 0 ;
  float a1, a2 ;
  ZZZ *up, *vp ;
  Array aa = zz->aa ;
  
  /* merge contiguous coding segments and equivalent coding/exonic segments 
   *  for the moment there is a single gene per segment so nothing to merge
   *
   * first pass merge coding and exonic annotations 
   */
  for (ii = 0, up = arrp (aa, 0, ZZZ) ; ii < arrayMax (aa) ; up++, ii++)
    {
      if (up->cluster)
	continue ;
      a1 = up->a1 ;
      a2 = up->a2 + zz->dx ;
      up->cluster = ++cluster ;

      /* look forward for equivalent segment */
      for (kk = ii + 1, vp = up+1 ; kk < arrayMax (aa) ; kk++, vp++)
	{
	  if (vp->map != up->map || vp->strand != up->strand || vp->a1 > a2)
	    continue ;
	  if (zz->frameSensitive && (((int)(up->a1 - a1)) % 3)) /* not in same frame */
	    continue ;
	  if (a2 < vp->a2 + zz->dx)
	    a2 = vp->a2 + zz->dx ;
	  vp->cluster = cluster ;
	}
    }
  return ;
} /* zzMerge */

/*************************************************************************************/

static void zzExport (ZZ *zz)
{
  AC_HANDLE  h = ac_new_handle () ;
  KEY old = 0, oldMap = 0 ;
  int ii, clMax = 0, oldStrand = 0 ;
  float a0, a1 = 0, a2 = 0 ;
  ZZZ *up ;

  freeOutf ("// samclustering :: %s\n\n", zz->classes) ;
  for (ii = 0, up = arrp (zz->aa, 0, ZZZ) ; ii < arrayMax (zz->aa) ; up++, ii++)
    {
      if (up->cluster != old)
	{
	  if (old)
	    {
	      if (oldStrand == -1)
		{ a0 = a1 ; a1 = a2 ; a2 = a0 ; }
	      if (!zz->dx)
		freeOutf ("Location %s %d %d\n"
			  , ac_protect (ac_key_name(oldMap), h)
			  , (int)a1, (int)a2) ;
	      else
		freeOutf ("Location %s %g %g\n"
			  , ac_protect (ac_key_name(oldMap), h)
			  , a1, a2) ;
	    }
	  freeOutf ("\nCluster %s%d\n", zz->prefix, up->cluster) ;
	  old = up->cluster ;
	  oldMap = up->map ;
	  oldStrand = up->strand ;
	  a1 = up->a1 ;
	}
      a2 = up->a2 ;
      if (old > clMax) clMax = old ;
      freeOutf ("%s %s\n", ac_key_class (up->gene), ac_protect (ac_key_name(up->gene), h)) ;
    }
  if (old)
    {
      if (oldStrand == -1)
	{ a0 = a1 ; a1 = a2 ; a2 = a0 ; }
	      if (!zz->dx)
		freeOutf ("Location %s %d %d\n"
			  , ac_protect (ac_key_name(oldMap), h)
			  , (int)a1, (int)a2) ;
	      else
		freeOutf ("Location %s %g %g\n"
			  , ac_protect (ac_key_name(oldMap), h)
			  , a1, a2) ;
    }
  ac_free (h) ;
  freeOutf ("\n// Exported %d objetcs in %d clusters\n", ii, clMax) ;
} /* zzExport */

/*************************************************************************************/
/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: samclustering [-db ACEDB] [-mixStrands] \n") ;
  fprintf (stderr, "// Example:  samclustering -db ZZ -mixStrands \n") ;  
  fprintf (stderr, "//   -db ACEDB database\n") ;
  fprintf (stderr, "//   -classes \'gene|mrna|foo|bar\' : list of classes to be clustered\n") ;
  fprintf (stderr, "//   -mixStrands : fuse the 2 strands\n") ;
  fprintf (stderr, "//   -frameSensitive : only fuse a and b in same frame x1a - x1b = 0 modulo 3\n") ;
  fprintf (stderr, "//   -float : the coordinates used are float rather than integers\n") ;
  fprintf (stderr, "//   -tag tagName : default=Location : the tag used as below in the schema\n") ;
  fprintf (stderr, "//   -prefix text : default=c_ : a prefix for the exported cluster names\n") ;
  fprintf (stderr, "//   -dx number : fuse object if 1+gap is smaller than dx\n") ;
  fprintf (stderr, "// Given a database containing a number of objects positioned on\n"
	   "// the genome using the schema\n"
	   "//     Location ?any_class UNIQUE Int UNIQUE Int\n"
	   "// cluster all these objects if they overlap, an export an ace file using the model\n"
	   "//    ?Cluster Child class_name object XREF Cluster\n"
	  ) ; 
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  FILE *f = 0 ;
  int outlevel = 0 ;
  ZZ zz ;
  AC_HANDLE h = 0 ;
  const char *ici ;
  char *s = "ok" ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&zz, 0, sizeof (ZZ)) ;
  zz.h = h ;

  /* optional temple argument */
  zz.mixStrands = getCmdLineOption (&argc, argv, "-mixStrands", 0) ;
  zz.frameSensitive = getCmdLineOption (&argc, argv, "-frameSensitive", 0) ;
  if (zz.frameSensitive &&  zz.mixStrands)
    {
      fprintf (stderr, "You cannot set at the same time frameSensitive and mixStrands\n") ;
      usage () ;
    }
 /* mandatory database descriptor */
  if (getCmdLineOption (&argc, argv, "-db", &ici))
    {
      zz.db = ac_open_db (messprintf("%s",ici), &s);
      if (!zz.db)
	messcrash ("Failed to open db %s, error %s", ici, s) ;
    }
  else 
     usage () ;

  if (getCmdLineOption (&argc, argv, "-classes", &ici))
    {
      zz.classes = strnew (ici, h) ;
    }
  else 
     usage () ;

  if (getCmdLineOption (&argc, argv, "-tag", &ici))
    zz.tag = strnew (ici, h) ;
   else
     zz.tag = "Location" ;
  if (getCmdLineOption (&argc, argv, "-prefix", &ici))
    zz.prefix = strnew (ici, h) ;
  else
    zz.prefix = "c_" ;
  if (getCmdLineOption (&argc, argv, "-dx", &ici))
    {
      if (sscanf (ici, "%g", &(zz.dx)) == 1 &&
	  zz.dx >= 0) ;
      else
	{
	  fprintf (stderr, "bad -dx parameter on command line, should be a positive number\n") ;
	  usage () ;
	}
    }

  zz.aa = arrayHandleCreate (20000, ZZZ, zz.h) ;
  outlevel = freeOutSetFile (stdout) ;	

  fprintf (stderr, "// start: %s\n", timeShowNow()) ;

  if (1)
    {
      zzGetData (&zz) ;
      arraySort (zz.aa, zzOrder) ;
      
      fprintf (stderr, "// data parsed, start merging: %s\n", timeShowNow()) ;
      zzMerge (&zz) ;
      arraySort (zz.aa, zzOrderByCluster) ;
      
      fprintf (stderr, "// merging done: %s\n", timeShowNow()) ;
      zzExport (&zz) ;
    }
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;

  if (outlevel)
    freeOutClose (outlevel) ;
  if (f) filclose (f) ;
  ac_db_close (zz.db) ;
  ac_free (zz.h) ;
  if (0) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

