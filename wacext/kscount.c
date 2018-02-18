#include "../wac/ac.h"

int main (int argc, char *argv[])
{
  const char *err = 0 ;
  char *target = "a:annie:2000101" ;
  AC_DB db = 0 ;
  AC_ITER iter ;
  AC_OBJ k ;
  AC_KEYSET ks ;
  if (argc > 1)
    target = argv[1] ;
  db = ac_open_db (target , &err) ;
  if (err)
    printf ("Error message %s", err) ;
  if (!db)
    messcrash ("could not open %s", target) ;

  {
    AC_KEYSET genes = ac_dbquery_keyset (db, "Find go_c IS \"membrane\" ; > gene", 0) ;
    AC_TABLE tt = ac_tablemaker_table (db, "web2genelist2", genes, ac_tablemaker_db, 0, 0, 0, 0) ;
    int ir ;
    
    for (ir = 0 ; tt && ir < tt->rows ; ir++)
      {
	printf ("%s" ,
		ac_table_printable (tt, ir, 0, "")
		) ;
	printf (" %s", 
		ac_table_printable (tt, ir, 5, "")
		) ;
	printf (" %s", 
		ac_table_printable (tt, ir, 6, "")
		) ;
	printf (" %s\n",
		ac_table_printable (tt, ir, 8, "")
		) ;
      }

    exit (0) ;
  }



  iter = ac_query_iter (db, 0, "find keyset __*", 0, 0) ;
  while ((k = ac_next_obj (iter)))
    {
      ks = ac_objquery_keyset (k, "expand", 0) ;
      printf ("%s\t%d\n", ac_name(k), ac_keyset_count (ks)) ;
      ac_free (ks) ;
      ac_free (k) ;      
    }
  ac_free (iter) ;
  ac_db_close (db) ;
}
