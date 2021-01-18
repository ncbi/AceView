#include <errno.h>
#include <ctype.h>

#include "../wac/ac.h"
/* #include "../wac/ac_.h" */
/* #include <wac/acclient_.h> */
#include <wh/dict.h>

AC_HANDLE h = 0 ;


static void ac_dump (char *cp, int n)
{
  int x;
  int addr;
  addr = 0;
  while (n > 0)
    {
      printf("%03x %4d  ",addr,addr);
      for (x=0; x< 16; x++)
        {
          if (x >= n)
            break;
          if (x % 4 == 0)
            printf(" ");
          printf("%02x ",cp[x]);
        }
      printf("\n          ");
      for (x=0; x< 16; x++)
        {
          if (x >= n)
            break;
          if (x % 4 == 0)
            printf(" ");
          printf(" %c ",isprint((int)cp[x]) ? cp[x] : '.' );
        }
      printf("\n");

      n -= 16;
      cp += 16;
      addr += 16;
    }
}


static void table_print(AC_TABLE t)
{
  int x,y ;
  char b[1000];

  printf("table: ");
  if (!t)
    {
      printf("	empty\n");
      return;
    }
  printf("%d rows %d cols\n",t->rows, t->cols);
  for (x=0; x < t->rows; x++)
    {
      for (y=0; y < t->cols; y++)
	{
	  sprintf(b,"-%s-", ac_table_printable(t,x,y,""));
	  printf("%-20s ",b);
	}
      printf("\n");
    }
}

static AC_DB db;
static char *db_server = "a:localhost:12345";


static void open_test()
{
  const char *s;
  extern char *db_server;
  if (!db)
    {
      printf ("test: open_test\n") ; 
      printf("USING DB SERVER %s\n",db_server);
      db = ac_open_db(db_server, &s);
    }
  if (!db)
    {
      printf("open error %s\nerrno %d\n",s,errno);
      exit(0);
    }
}

static void bql_test_1 (void)
{
  AC_TABLE t, t1;
  AC_HANDLE hh = ac_new_handle () ;
  AC_KEYSET aks = ac_dbquery_keyset (db, "Find bql1 ;  IS > a2", hh) ;
  char *cmd = 
    "select a, i, f, t, d, o from a in class bql1, i in a->ivalue, f in a->fvalue, t in a->tvalue, d in a->dvalue, o in a->ovalue";
  char *cmd1 = 
    "select a, i, f, t, d, o from a in @active:1, i in a->ivalue, f in a->fvalue, t in a->tvalue, d in a->dvalue, o in a->ovalue";

  if (1)
    {
      printf("test: bql_test_1\n");
      t = ac_bql_table (db, cmd, NULL, 0, 0, hh);
      table_print(t);
    }
  if (1)
    {
      printf("test: bql_test_1, -active\n");
      t1 = ac_bql_table (db, cmd1, aks, 0, 0, hh);
      table_print(t1);
    }

  ac_free(hh);
}

static void tblmaker_test (void)
{
  AC_TABLE t, t1;
  AC_HANDLE h = ac_new_handle () ;

  if (0)
    {
      printf("test: tblmaker_test\n");
      
      t = ac_tablemaker_table (db, "test.def", 0, ac_tablemaker_file, 0, 0, 0, h) ;
      table_print(t);
      
      t1 = ac_tablemaker_table (db, "test.def", 0, ac_tablemaker_file, "a1", 0, 0, h) ;
      table_print(t1);
    }
  ac_free (h) ;
}

static void it_test_1()
{
  AC_ITER i;
  AC_OBJ o;

  printf("test: it_test_1\n");

  i = ac_query_iter (db, 1, "find bql1", NULL, h );
  while ((o = ac_next_obj(i)))
    {
      printf("\t%s:%s %d %f\n",ac_class(o),ac_name(o),ac_tag_int(o,"ivalue", -1),
		ac_tag_float(o,"fvalue",-1.0));
      ac_free(o);
    }
  if (!h)
    ac_free(i);

  i = ac_query_iter( db, 0, "find bql1", NULL, h );
  while ((o = ac_next_obj(i)))
    {
      printf("\t%s:%s %d %f\n",ac_class(o),ac_name(o),ac_tag_int(o,"ivalue", -1),
		ac_tag_float(o,"fvalue",-1.0));
      ac_free(o);
    }
  if (!h)
    ac_free(i);
}

static void ks_print(char *name, AC_KEYSET ks)
{
  AC_ITER i;
  AC_OBJ o;
  int n;

  printf("Keyset %s: %d elements\n",name, ks ? ac_keyset_count(ks) : 0);
  if (ks)
    {
      i = ac_keyset_iter (ks, 0, h );
      n = 0;
      while ((o = ac_next_obj(i)))
	{
	  printf("\t%s:%s\n",ac_class(o),ac_name(o));
	  ac_free(o);
	  n++;
	}
      if (n != ac_keyset_count(ks))
	printf("YOW! number of objects iterated (%d) not equal to count (%d)\n",n,ac_keyset_count(ks));
      if (!h)
	ac_free(i);
    }
}

static void ks_test_1()
{
  AC_KEYSET ks1, ks2, ks3, ks4;
  AC_OBJ bql1_a1, bql1_a2, bql1_a3, bql1_a4;

  printf("test: ks1\n");

  ks1 = ac_new_keyset(db, h);
  ks2 = ac_new_keyset(db, h);

  printf("new keyset has %d elements\n",ac_keyset_count(ks1));

  bql1_a1 = ac_get_obj( db, "bql1", "a1", h);
  bql1_a2 = ac_get_obj( db, "bql1", "a2", h);
  bql1_a3 = ac_get_obj( db, "bql1", "a3", h);
  bql1_a4 = ac_get_obj( db, "bql1", "a4", h);

  printf("objects found: %s %s %s %s\n",
	 bql1_a1  != NULL? "y" : "n",
	 bql1_a2  != NULL? "y" : "n",
	 bql1_a3  != NULL? "y" : "n",
	 bql1_a4  != NULL  ? "y" : "n" );

  ac_keyset_add( ks1, bql1_a1);
  ac_keyset_add( ks1, bql1_a2);

  ks_print("ks1",ks1);

  ac_keyset_add( ks2, bql1_a3);
  ac_keyset_add( ks2, bql1_a4);

  ks_print("ks2",ks2);

  ks3 = ac_copy_keyset (ks1, h);

  ks_print("ks3",ks3);

  ac_keyset_or(ks1, ks2);

  ks_print("ks1 = ks1 OR ks2", ks1);

  ac_keyset_and(ks1, ks2);

  ks_print("ks1 = ( ks1 OR ks2) AND ks2", ks1);

  ac_keyset_xor(ks1, ks3);

  ks_print("ks1 = ( ( ks1 OR ks2) AND ks2 ) OR ks3", ks1);

  ks4 = ac_copy_keyset (ks1, h);

  ks_print("ks4",ks4);

  ac_keyset_remove(ks4, bql1_a2);

  ks_print("ks4 remove a2", ks4);

  ac_keyset_minus( ks1, ks4 );

  ks_print("ks1 = ks1 - ks4",ks1);

  if (!h)
    {
      printf("freeing\n");
      ac_free(ks1);
      ac_free(ks2);
      ac_free(ks3);
      ac_free(ks4);
      ac_free(bql1_a1);
      ac_free(bql1_a2);
      ac_free(bql1_a3);
      ac_free(bql1_a4);
    }
}

static void obj_test_1()
{
  AC_OBJ obj, obj2;
  AC_TABLE tb1, tb2, tNest1;

  printf ("test: obj_test_1\n");

  obj = ac_get_obj(db, "arf", "a", h);
  printf("obj found: %s\n", obj  != NULL ? "y" : "n");

  printf("%s %s\n",ac_class(obj), ac_name(obj));

  if (ac_has_tag(obj, "arbitrary_tag"))
    printf("has arbitrary_tag\n");
  else
    printf("YOW! does not have arbitrary_tag\n");

  if (ac_tag_type(obj, "ivalue") != ac_type_int)
    printf("YOW! ivalue tag not followed by integer\n");

  printf("ivalue %d\n",ac_tag_int(obj, "ivalue", -1));
  printf("fvalue %g\n",ac_tag_float(obj, "fvalue", -1.0));
  printf("tvalue %s\n",ac_tag_text(obj, "tvalue", "(null)"));
  printf("dvalue %d\n",ac_tag_date(obj, "dvalue", -1));

  printf("printable ivalue %s\n",ac_tag_printable(obj, "ivalue", "yow!"));
  printf("printable fvalue %s\n",ac_tag_printable(obj, "fvalue", "yow!"));
  printf("printable tvalue %s\n",ac_tag_printable(obj, "tvalue", "yow!"));
  printf("printable dvalue %s\n",ac_tag_printable(obj, "dvalue", "yow!"));

  if ((tNest1 = ac_tag_table (obj, "nest1", 0)))
    {
      printf("printable name (testing #constructs, should find one)    %s\n", ac_table_printable(tNest1, 0, 1, "yow!"));
      ac_free (tNest1) ;
    }
  else
    printf("nest1 failed to construct ac_tag_table on # construct\n") ;
      printf("printable name (testing #constructs, should find dude)    %s\n", ac_tag_printable(obj, "name", "yow!"));

  obj2 = ac_tag_obj(obj, "ovalue", h);

  printf("obj2 = %s %s\n",ac_class(obj2), ac_name(obj2));

  printf("obj2 table:\n");
  tb1 = ac_tag_table( obj, NULL, h);
  table_print(tb1);

  if (ac_tag_type(obj2, "ovalue") != ac_type_empty)
    printf("YOW! ovalue in arf:b object\n");

  printf("obj2 table:\n");
  tb2 = ac_tag_table( obj, NULL, h);
  table_print(tb2);

  if (!h)
    {
      ac_free(obj);
      ac_free(obj2);
      ac_free(tb1);
      ac_free(tb2);
    }
}

static void obj_test_2()
{
  AC_OBJ o;
  AC_TABLE tb;

  printf("obj_test_2:\n");
  o = ac_get_obj(db, "arf", "self_ref", h);
  tb = ac_tag_table( o, NULL, h);
  table_print(tb);

  if (!h)
    {
      printf ("obj_test_2 calls ac-free\n") ;
      ac_free(tb);
      ac_free(o);
    }
  else
    printf ("obj_test_2 h!=0\n") ;
}

static void obj_test_3()
{
  AC_OBJ o;
  AC_TABLE tb;
  printf("obj_test_3:\n");
  o = ac_get_obj(db, "arf", "a", NULL);
  tb = ac_tag_table( o, NULL, NULL);
  table_print(tb);
  ac_free(o);
  ac_free(tb);
}

static void kt_test_1()
{
  AC_TABLE t;
  AC_KEYSET ks;
  int x;
  AC_OBJ o;
  
  printf("test: kt1\n");
  
  if (1)
    {
      printf("find bql1\n");
      ks = ac_dbquery_keyset (db, "find bql1", h) ;
      
      t = ac_keyset_table( ks, 0, -1, TRUE, h);
      
      table_print(t);
      
      if (!h)
	{
	  ac_free(t);
	  ac_free(ks);
	}
    }

  ks = ac_dbquery_keyset (db, "find lots", h) ;

  t = 1 ? ac_keyset_table( ks, 0, 2000, TRUE, h) : 0 ;

  printf("find lots:\n");
  if (t) printf("%d rows %d cols\n",t->rows, t->cols);

  for (x=0; t && x < 30; x++)
    {
      o = ac_table_obj(t, x, 0, NULL);
      if (ac_table_type(t,x,1) != ac_type_empty)
	printf("YOW! row %d col 1 not empty\n",x);
      printf("\t%s\n",ac_name(o));
      ac_free(o);
    }

  if (!h)
    {
      ac_free(t);
      ac_free(ks);
    }
}

static struct a
{
  int num;
  char *name;
}
dict_values[] = 
  {
    { 1, "Global" },
    { 2, "Session" },
    { 3, "Voc" },
    {23, "Display" },
    {24, "MainClasses" },
    {25, "Bat" },
    {26, "KeySet" },
    {27, "Calcul" },
    {28, "Class" },
    {29, "Model" },
    {30, "Text" },
    {31, "LongText" },
    {32, "Image" },
    {33, "Table" },
    {34, "TableResult" },
    {35, "Jade" },
    {36, "View" },
    {37, "FicheView" },
    {38, "Comment" },
    {39, "UserSession" },
    {40, "Query" },
    {41, "Constraint" },
    {42, "Peptide" },
    {43, "Sequence" },
    {44, "Protein" },
    {45, "DNA" },
    {46, "Paper" },
    {47, "Method" },
    {48, "Map" },
    {49, "gMap" },
    {50, "vMap" },
    {51, "MultiMap" },
    {52, "Locus" },
    {53, "Gene" },
    {54, "Allele" },
    {55, "Interval" },
    {56, "2_point_data" },
    {57, "Multi_pt_data" },
    {58, "Clone" },
    {59, "Clone_Grid" },
    {60, "Pool" },
    {61, "Contig" },
    {62, "pMap" },
    {63, "Chrom_Band" },
    {64, "Motif" },
    {65, "BaseCall" },
    {66, "BaseQuality" },
    {67, "BasePosition" },
    {68, "OligoRepeat" },
    {69, "Homology_group" },
    {70, "Map_set" },
    {71, "Doc" },
    {72, "Genetic_code" },
    {73, "Person" },
    {74, "Colour" },
    {75, "Keyword" },
    {76, "SourceCode" },
    {77, "Include" },
    {78, "MatchTable" },
    {79, "cMap" },
    {80, "Tag" },
    {81, "Table_definition" },
    {82, "arf" },
    {83, "foo" },
    {84, "bql1" },
    {85, "foobar" },
    {86, "lots" },
    {-1, 0}
  };

static void dict_test_2()
{
  int n,x;
  DICT *dict;

  printf("test: dict 2\n");

  dict = dictHandleCreate( 256, 0 );

  for (n=0; n < 256; n++)
    {
      for (x=0; dict_values[x].num != -1; x++)
	if (dict_values[x].num == n)
	  goto found;
      dictAdd(dict, "", &x);
      continue;
    found:
      dictAdd(dict, dict_values[x].name, &x);
      continue;
    }

  for (x=0; dict_values[x].num != -1; x++)
    {
      if (! dictFind(dict, dict_values[x].name, &n) )
	printf("YOW! dict fails to find %d -%s-\n",x,dict_values[x].name);
    }

  dictDestroy(dict);
  printf("ok\n");
}

static void tag_table_test_1()
{
  AC_OBJ o;
  AC_TABLE t;

  printf("test: tt1\n");
  printf("w\n");

  t = 0;

  o = ac_get_obj(db, "wide", "w", h);
  t = ac_tag_table(o, "arfs", h);
  table_print(t);

  if (!h)
    {
      ac_free(t);
      ac_free(o);
    }
}


static void tag_table_test_2()
{
  AC_OBJ o;
  AC_TABLE t;
  printf("w1\n");
  o = ac_get_obj(db, "wide", "w1", h);
  t = ac_tag_table(o, "arfs", h);
  if (t)
    {
      printf("YOW! should be empty table\n");
      table_print(t);
      if (!h)
	ac_free(t);
    }
  if (!h)
    ac_free(o);
}


static void tag_table_test_3()
{
  AC_OBJ o;
  AC_TABLE t;

  printf("w2\n");
  o = ac_get_obj(db, "wide", "w2", h);
  t = ac_tag_table(o, "arfs", h);
  table_print(t);

  if (!h)
    ac_free(t);
  if (!h)
    ac_free(o);

}

static void make_handle()
{
  h = ac_new_handle ();
}

static void close_handle ()
{
  ac_free (h) ;
}

static void close_db()
{
  ac_db_close (db) ;
}

static void all_tests();
static void all_handle();

typedef struct fn_list
{
  void (*fn)();
  char *name;
  int special;
} ALLFN ;
static ALLFN* allfn ;

static void do_all(int use_handle)
{
  int x;
  for (x=0; allfn[x].fn; x++)
    {
      switch (allfn[x].special)
	{
	case 2:
	  break;
	case 1:
	  if (use_handle)
	    (*allfn[x].fn)();
	  break;
	case 0:
	case 3:
	  (*allfn[x].fn)();
	  break;
	}
    }
}

static void all_tests()
{
  do_all(0);
}

static void all_handle()
{
  do_all(1);
}

static void command_keyset()
{
  AC_KEYSET ks;
  AC_TABLE tbl;
  
  printf("test: ck\n");
  printf("find bql1:\n");
  ks = ac_command_keyset (db, "find bql1", NULL, h);
  tbl = ac_keyset_table( ks, 0, -1, FALSE, h);
  
  table_print(tbl);
  
  if (!h)
    { 
      ac_free(ks);
      ac_free(tbl);
    }
  
  printf("find arf:\n");
  ks = ac_command_keyset (db, "find arf", NULL, h);
  tbl = ac_keyset_table( ks, 0, -1, FALSE, h);
  table_print(tbl);
  
  if (!h)
    { 
      ac_free(ks);
      ac_free(tbl);
    }

}

#include <sys/time.h>

static void transaction_time()
{
  struct timeval tstart, topen, tloop, tclose;
  int x;
  int t;
  int N = 1000;
  
  if (gettimeofday(&tstart, 0) < 0)
    perror("gettimeofday");
  open_test();
  
  if (gettimeofday(&topen, 0) < 0)
    perror("gettimeofday");
  
  for (x=0; x< N; x++)
    {
      unsigned char *s;
      s = ac_command(db, "help", 0, 0);
      messfree(s);
    }
  
  if (gettimeofday(&tloop, 0) < 0)
    perror("gettimeofday");
  
  ac_db_close (db) ;
  
  if (gettimeofday(&tclose, 0) < 0)
    perror("gettimeofday");
  
  t = ( topen.tv_sec - tstart.tv_sec ) * 1000000 + ( topen.tv_usec - tstart.tv_usec );
  printf("open  %2d.%06d\n",t / 1000000, t % 1000000);
  t = ( tloop.tv_sec - topen.tv_sec ) * 1000000 + ( tloop.tv_usec - topen.tv_usec );
  printf("loop  %2d.%06d ( N = %d )\n",t / 1000000, t % 1000000, N);
  t = t / N ; 
  printf("avg   %2d.%06d\n",t / 1000000, t % 1000000);
  t = ( tclose.tv_sec - tloop.tv_sec ) * 1000000 + ( tclose.tv_usec - tloop.tv_usec );
  printf("close %2d.%06d\n",t / 1000000, t % 1000000);
}


static void command_test()
{
#if 0
  unsigned char *s;
  int l,e;
  // s = (db->ac_partial_command)(db, "find arf\nfind wide\nfind bql1\n", &l, &e);
  s = (db->ac_partial_command)(db, "find arf", &l, &e);
  printf("encore %d\n",e);
  printf("response:\n%s\n",s);
  ac_dump(s, l);
#else
  printf("no command test\n");
#endif
}

static void stack_test()
{
  Stack s;
  char *st;
  int n;
  
  s = stackCreate(1);

  pushText(s,"aaa");
  pushText(s,"bbb");
  catText(s,"BBB");
  pushText(s,"ccc");
  
  st = popText(s);
  printf("pop %s\n",st);
  
  pushText(s,"ddd");
  printf("st = %s\n",st);
  
  stackCursor(s, 0);
  while ((st = stackNextText(s)))
    printf("->%s<-\n",st);
  
  stackClear(s);
  
  st = strdup("this is a  test   string");
  stackTokeniseTextOn( s, st, " ");
  
  printf("tokenize:\n");
  
  stackCursor(s, 0);
  while ((st = stackNextText(s)))
    printf("->%s<-\n",st);
  
  printf("catBinary:\n");
  s = stackCreate(1);
  catBinary (s,"\0",1) ;
  catBinary(s,"\0\0\0\0\0\0\0\0",3);
  catBinary(s,"x",1);
  catBinary(s,"y",1);
  n = stackMark(s);
  st = stackText(s, 0);
  ac_dump (st, n);
  
  printf("catBinary 2:\n");
  s = stackCreate(1);
  pushText(s,"a");
  catBinary(s,"\0\0",2);
  catText(s,"b");
  n = stackMark(s);
  st = stackText(s, 0);
  ac_dump(st, n);
  
  printf("catBinary 3:\n");
  s = stackCreate(1);
  catBinary(s,"1234",4);
  n = stackMark(s);
  printf("found %d bytes\n",n);
  st = stackText(s, 0);
  ac_dump(st, n);
  
}

static void dna_test()
{
  AC_OBJ product, sequence ;
  char *dna, *peptide ;
  
  product = ac_get_obj(db, "Protein", "pep1", NULL);
  peptide = ac_obj_peptide (product, NULL) ;
 
  
  printf("Product pep1\n%s\n", peptide ? peptide : "Missing peptide");

  
  sequence = ac_get_obj(db, "Sequence", "dna1", NULL);
  dna = ac_obj_dna  (sequence, NULL) ;
  printf("Sequence dna1 dna:\n%s\n", dna ? dna : "missing dna");
  if (1)
    {
      dna = ac_zone_dna  (sequence, 6, 14, NULL) ;
      printf("Zone Sequence dna1 6->14 dna:\n%s\n", dna ? dna : "missing zone dna");
      dna = ac_zone_dna  (sequence, 14, 6, NULL) ;
      printf("Reverse Zone Sequence dna1 14->6 dna:\n%s\n", dna ? dna : "missing reverse zone dna");
    }
}

static void rawcmd()
{
  unsigned char *s;
  int len;
  printf("raw command test: ");
  s = ac_command( db, "help", &len, h);
  if (strlen((char*)s) != len && strlen((char*)s) != len - 1)
    printf("length of response not equal to specified length\n");
  else
    printf("ok\n");
  if (!h)
    messfree (s) ;
}

static void parse_test()
{
  AC_KEYSET ks1 = 0;
  const char *etext;
  BOOL b;
  
  b = ac_parse(db, "load : a\ntag4 hello 1\n\n", &etext, &ks1, NULL_HANDLE );
  ks_print("parse_test 1", ks1);
  printf("success %d\n",b);
  printf("etext:\n%s\n",etext);
  ac_free(ks1);
  messfree(etext);
  
  b = ac_parse(db, "load : b\ntag4 b 1\n\nload : c\nsnarf t\ntag4 a\n", &etext, &ks1, NULL_HANDLE );
  ks_print("parse_test 2", ks1);
  printf("success %d\n",b);
  printf("etext:\n%s",etext);
  ac_free(ks1);
  messfree(etext);
}


static void allfnInit (void) ;
int main(int argc, char **argv)
{
  int x;
  if (argv[0]) argv++;
  
  if (argv[0])
    {
      db_server = argv[0];
      argv++;
    }
  
  if (! argv[0] )
    {
      printf("use: acctest testname\n");
      printf("\tsee definition of allfn in acctest.c\n");
      exit(1);
    }
  allfnInit () ;
  open_test () ;
  for ( ; argv[0] ; argv++ )
    {
      for (x=0; allfn[x].fn; x++)
	if (strcmp(allfn[x].name,argv[0]) == 0)
	  {
	    (*allfn[x].fn)();
	    break;
	  }
      if (! allfn[x].fn)
	printf("YOW! don't recognize %s\n",argv[0]);
    }
  
  return 0 ;
}


/*
 * function - function to call
 * name - name of test given by user
 * special - 
 *	0 ordinary test
 *	1 only perform if handle
 *	2 not included when running all tests
 */
static void allfnInit (void)
{
  static ALLFN allfn2[] =
  {
    { open_test,		"o",		3 },
    { make_handle,		"h", 		1 },
    { bql_test_1,		"bql1",		0 },
    { tblmaker_test,		"tblmaker",	0 },
    { ks_test_1,		"ks1", 		0 },
    { obj_test_1,		"obj1",		0 },
    { obj_test_2,		"obj2",		3 },
    { obj_test_3,		"obj3",		0 },
    { kt_test_1,		"kt1", 		3 },
    { tag_table_test_1,		"tt1", 		0 },
    { tag_table_test_2,		"tt2", 		0 },
    { tag_table_test_3,		"tt3", 		0 },
    { command_keyset,		"ck1", 		0 },
    { it_test_1,		"it1",		0 },
    { rawcmd,			"rawcmd",	0 },
    { parse_test, 		"parse",	0 },
    { dict_test_2,		"dict2",	0 },
    { dna_test,		        "dna",		0 },
    { stack_test,               "stack",        0 },

    { close_handle,		"ch", 		1 },
    { close_db,			"c", 		3 },

    { all_tests,		"all", 		2 },
    { all_handle,		"allh",		2 },

    { command_test,		"ct",		2 },
    { transaction_time,		"trans",	2 },


    { 0, 0 }
  };
  allfn = allfn2 ;
}
