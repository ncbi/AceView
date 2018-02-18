#include "../wac/ac.h"
#include <errno.h>



static void showTable (AC_TABLE t, int x1, int x2)
{
  int x, y ;
  const char *ccp;

  printf("show  (%d,%d) (rows,cols)\n", t->rows, t->cols) ;
  if (t) for (x=x1; x < t->rows && x <= x2; x++)
    {
      for (y=0; y < t->cols+ 10; y++)
	{
	  ccp = ac_table_printable(t,x,y,"toto") ;
	  printf("%s", ccp) ;
	  printf("\t");
	}
      printf("\n");
    }
}

static int ontofind (AC_DB db)
{
  AC_TABLE tt = 0 ;
  int x, nn = 0, nm = 0 ;
  const char *newn, *go_m ;
  unsigned const char *s ;
  if (0) tt = ac_bql_table (db, "select nn, go_m, go_b, go_c from go in class \"GO\" , nn in go->nickname where nn, go_m in go->go_m where  go_m, go_b in go->go_b , go_c in go->go_c", 0, 0, 0, 0) ;
  tt = ac_bql_table (db, "select go, nn, go_m from go in class \"GO\" , nn in go->nickname where nn, go_m in go->Molecular_function where go_m", 0, 0, 0, 0) ;

  if (!tt)
    {
      printf ("bql failed") ;
      return 0 ;
    }
  
  printf("ac_bql_done: found (%d,%d) (rows,cols)\n", tt->rows, tt->cols) ;
  for (x=0 ;  x < tt->rows ; x++)
    {
      showTable (tt, x, x) ;
      /* go = ac_table_printable (tt,x,0,"") ; */
      newn = ac_table_printable (tt,x,1,"") ;
      if (!*newn)
	continue ;
      nn++ ; nm++ ;
      s = ac_command (db, messprintf ("webquery extended \"%s\"", newn), 0, 0) ;
      messfree (s) ;
      go_m = ac_table_printable (tt,x,2,"") ;
      if (1)
	s = ac_command (db, messprintf ("edit GO_m_ace %s", go_m), 0, 0) ;
      else
	printf ("edit GO_m_ace %s\n", go_m) ;
      messfree (s) ;
    }
  printf ("Found %d go_m term\n", nm) ;
  ac_free (tt) ;
  return nn ;
}


int main (int argc, char **argv)
{
  const char *errors ;
  char *dbName = argc>=2 ? argv[1] : "." ;
  AC_DB db = ac_open_db (dbName, &errors);
  if (!db)
    messcrash ("Failed to open db %s, error %s", dbName, errors) ;

  ontofind (db) ;

  ac_db_close (db) ;
  return 0 ;
}

