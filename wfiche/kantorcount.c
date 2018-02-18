#include "../wac/ac.h"

int main(int argc,char **argv)
{
  AC_DB db;
  AC_KEYSET ks1, tks;
  const char *err;
  int prod, blastp, taxblast, pfam, expasy, psort;
  
  if (! argv[1] || !argv[2])
    {
      fprintf(stderr,"usage: kantorcount database_specifier query\n");
      exit(1);
    }
  
  db=ac_open_db (argv[1], &err);
  if (!db)
    {
      fprintf(stderr,"cannot open database: %s\n",err);
      exit(1);
    }
  
  ks1 = ac_dbquery_keyset( db, argv[2], NULL );
  
  prod = ac_keyset_count( ks1 );
  
  tks = ac_ksquery_keyset(ks1, "! Blastp_date", NULL );
  
  blastp = ac_keyset_count( tks );
  
  tks = ac_ksquery_keyset(ks1, "Blastp && !Taxblast_date", NULL );
  
  taxblast = ac_keyset_count( tks );
  
  tks = ac_ksquery_keyset(ks1, "!Pfam_date", NULL );
  
  pfam = ac_keyset_count( tks );
  
  tks = ac_ksquery_keyset(ks1, "!PfamNew_date", NULL );
  
  /*   pfamnew = ac_keyset_count( tks ); */
  
  tks = ac_ksquery_keyset(ks1, "!Expasy_date", NULL );
  
  expasy = ac_keyset_count( tks );
  
  tks = ac_ksquery_keyset(ks1, "!Psort_date" , NULL );
  
  psort = ac_keyset_count( tks );
  
  printf("\nNumber of genes to be kantorized#\n\n");
  printf("               Missing       Done      Total\n");
  printf ("Blastp      %10d %10d %10d\n", blastp, 	prod - blastp, 		prod) ;
  printf ("Taxblast    %10d %10d %10d\n", taxblast,	prod - taxblast,	prod) ;
  printf ("Pfam        %10d %10d %10d\n", pfam, 		prod - pfam, 		prod) ;
  printf ("Expasy      %10d %10d %10d\n", expasy, 	prod - expasy, 		prod) ;
  printf ("Psort       %10d %10d %10d\n", psort,		prod - psort, 		prod) ;

  ac_db_close (db) ;
  return 0 ;
}
