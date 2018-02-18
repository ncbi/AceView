#include "../wac/ac.h"

int main (int argc, char *argv[])
{
  int ir ;
  char *cp, *cq ;
  const char *err = 0 ;
  char *target = "" ;
  AC_DB db = 0 ;
  AC_ITER iter ;
  AC_OBJ k = 0 ;
  AC_TABLE gP ;

  if (argc > 1)
    target = argv[1] ;
  db = ac_open_db (target , &err) ;
  if (err)
    printf ("Error message %s", err) ;
  if (!db)
    messcrash ("could not open %s", target) ;

  iter = ac_query_iter (db, 0, "find kantor acekog && product", 0, 0) ;
  while (ac_free (k), (k = ac_next_obj (iter)))
    {
      gP = ac_tag_table (k, "Product", 0) ;
      for (ir = 0 ; gP && ir < gP->rows ; ir++)
	{
	  if ((cp = ac_tag_ace (k, "AKG", 0)))
	    {
	      printf ("Product \"%s\"\n"
		      , ac_table_printable (gP, ir, 0, "")
		      ) ;
	      cq = strtok (cp, "\n") ;
	      do
		{
		  if (!strstr (cq, "AceKog_Date"))
		    printf ("%s\n", cq) ;
		} while ((cq = strtok (0, "\n"))) ;
	      printf ("\n") ;
	      ac_free (cp)  ;
	    }
	}

      if (0 && (cp = ac_obj_ace (k, 0)))
	printf ("/////\n%s/////\n", cp) ;
      ac_free (cp)  ;
    }
  ac_free (iter) ;
  ac_db_close (db) ;
  
  return 0 ;
}
