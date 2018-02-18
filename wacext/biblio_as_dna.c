#include "../wac/ac.h"


static void encodeAbstract (AC_DB db, AC_OBJ gene)
{
  AC_HANDLE h = handleCreate () ;
  Array aa = 0 ;
  char *lt = 0 ;
  int jj = 0 ;
  char *cp ;
  /* transform the gene into a tg, add its nms as clones */
  
  aa = arrayCreate (3000, char) ;

 
  lt = ac_longtext (gene, h) ;
  
  cp = lt ;
  /*
  if ((cp = strstr (lt, "***LongTextEnd***")))
    *cp = 0 ;
  cp = strstr (lt, "LongText") ;
  while (*cp && *cp != '\n') cp++ ; 
  */
  while (*++cp)
    if (pepEncodeChar[(int) ace_upper(*cp)] > -1)
      array (aa, jj++, char) = ace_upper(*cp) ;
  array (aa, jj++, char) = 0 ;

  if (jj > 80)
    {
      printf (">%s\n", ac_name(gene)) ;
      printf ("%s\n", arrp (aa, 0, char)) ;
    }
  arrayDestroy (aa) ;
  messfree (h) ;
}

int main (int argc, char *argv[])
{
  const char *err = 0 ;
  char *target = "a:annie:2000101" ;
  AC_DB db = 0 ;
  AC_ITER iter ;
  AC_OBJ k ;
  int nn = 0 ;
  
  if (argc > 1)
    target = argv[1] ;
  db = ac_open_db (target , &err) ;
  if (err)
    printf ("Error message %s", err) ;
  if (!db)
    messcrash ("could not open %s", target) ;


  iter = ac_query_iter (db, 0, "find paper; >abstract", 0, 0) ;
  while (--nn && (k = ac_next_obj (iter)))
    {
      encodeAbstract (db, k) ;
      ac_free (k) ;      
    }
  ac_free (iter) ;
  ac_db_close (db) ;

  return 0 ;
}
