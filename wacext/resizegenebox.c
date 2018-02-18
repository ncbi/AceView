#include "../wac/ac.h"

typedef struct { char nam[256] ; int x1, x2 ;} GG ;

static void showGG (Array aa)
{
  GG *gg ;
  int ir ;

  printf ("\n") ;
  if (aa && arrayMax (aa))
    for (ir = 0 ; aa && ir < arrayMax (aa) ; ir++)
      {
	gg = arrayp (aa, ir, GG) ;
	printf ("%d\t%d\t%s\n", gg->x1, gg->x2, gg->nam) ;
      }
  printf ("\n") ;
}

static BOOL exportGene (AC_DB db, AC_OBJ gene)
{
  AC_HANDLE h = handleCreate () ;
  AC_TABLE xx, pgs, imap ;
  AC_OBJ parent = 0, chrommap = 0, tg, mrna ;
  int ir, jr, a1, a2, da1, da2, x1, x2 ;
  BOOL isTg, isDown = TRUE, ok = FALSE ;
  Array aa = 0 ;
  GG *gg ;
  BOOL debug = FALSE ;

  /* transform the gene into a tg, add its nms as clones */

      /* find the size of the gene->mrnas or the gene->pgs */
  isTg = ac_has_tag (gene, "Transcribed_gene") ;
  a1 = a2 = -1 ;
  aa = arrayCreate (12, GG) ;
  if (isTg)
    {
      tg = ac_tag_obj (gene, "Transcribed_gene", h) ;
   
      imap = ac_tag_table (tg, "IntMap", 0) ;
      if (imap && imap->cols >= 2)
	{
	  if (!chrommap) 
	    chrommap = ac_table_obj (imap, 0, 0, h) ;
	  gg = arrayp (aa, arrayMax (aa) , GG) ;
	  strcpy (gg->nam, ac_name(tg)) ;
	  gg->x1 = ac_table_int (imap, 0, 1, 0) ;
	  gg->x2 = ac_table_int (imap, 0, 2, 0) ;
	}
      ac_free (imap) ;
    }
  else
    {
      pgs = ac_tag_table (gene, "Genefinder", h) ;
      for (jr = 0 ; pgs && jr < pgs->rows ; jr++)
	{
	  mrna = ac_table_obj (pgs, jr, 0, h) ;
	  imap = ac_tag_table (mrna, "IntMap", 0) ;
	  if (imap && imap->cols >= 2)
	    {
	      if (!chrommap) 
		chrommap = ac_table_obj (imap, 0, 0, h) ;
	      gg = arrayp (aa, arrayMax (aa) , GG) ;
	      strcpy (gg->nam, ac_name(mrna)) ;
	      gg->x1 = ac_table_int (imap, 0, 1, 0) ;
	      gg->x2 = ac_table_int (imap, 0, 2, 0) ;
	    }
	  ac_free (imap) ;
	}
    }

  if (debug) showGG (aa) ;
  if (arrayMax (aa))
    {
      gg = arrayp (aa, 0, GG) ;
      a1 = gg->x1 ; a2 = gg->x2 ;
      if (a1 < a2) isDown = TRUE ;
      else isDown = FALSE ;      

      for (ir = 0 ; aa && ir < arrayMax (aa) ; ir++)
	{
	  gg = arrayp (aa, ir, GG) ;
	  if (isDown)
	    {
	      if (a1 > gg->x1) a1 = gg->x1 ;
	      if (a2 < gg->x2) a2 = gg->x2 ;
	    }
	  else
	    {
	      if (a1 < gg->x1) a1 = gg->x1 ;
	      if (a2 > gg->x2) a2 = gg->x2 ;
	    }
	}
      
      if (debug) printf ("Maxs a1=%d a2 = %d\n\n", a1, a2) ;
      imap = ac_tag_table (gene, "IntMap", 0) ;
      if (imap && imap->cols >= 2)
	{
	  if (!chrommap) 
	    chrommap = ac_table_obj (imap, 0, 0, h) ;
	  gg = arrayp (aa, arrayMax (aa) , GG) ;
	  gg->x1 = ac_table_int (imap, 0, 1, 0) ;
	  gg->x2 = ac_table_int (imap, 0, 2, 0) ;
	  
	  da1 = a1 - gg->x1 ;
	  da2 = a2 - gg->x2 ; 
	  if (debug) printf ("Maxs da1=%d da2 = %d\n\n", da1, da2) ;
	  if (!isDown) {  da1 = - da1 ; da2 = - da2 ; }
	  if (da1 || da2)
	    {
	      /* edit the intmap */
	      printf ("Gene \"%s\"\n", ac_name(gene)) ;
	      printf ("IntMap \"%s\" %d %d // da1=%d da2=%d\n\n"
		      , ac_name(chrommap), a1, a2, da1, da2) ;
	      
	      
	      /* edit the coords in the parent */
	      parent = ac_tag_obj (gene, "Genomic_sequence", h) ;
	      xx = ac_tag_table (parent, "Genes", 0) ;
	      for (jr = 0 ; jr < xx->rows ; jr++)
		{
		  if (!strcmp (ac_name(gene), ac_table_printable (xx, jr, 0, "")))
		    {
		      /* current coords of gene in parent */
		      x1 = ac_table_int (xx, jr, 1, 0) ;
		      x2 = ac_table_int (xx, jr, 2, 0) ;
		      
		      printf ("Sequence \"%s\"\n", ac_name(parent)) ;
		      printf ("-D Genes %s %d %d\n"
			      , ac_name(gene), x1, x2) ;
		      
		      if (x1 < x2) { x1 += da1 ; x2 += da2 ; }
		      else { x1 -= da1 ; x2 -= da2 ; }
		      printf ("Genes %s %d %d\n"
			      , ac_name(gene), x1, x2) ;
		      printf ("\n") ;
		      ok = TRUE ;
		      break ;
		    }
		  
		}
	      ac_free (xx) ;
	    }
	}
      ac_free (imap) ;
    }
  arrayDestroy (aa) ;
  messfree (h) ;

  return ok ;
}


int main (int argc, char *argv[])
{
  const char *err = 0 ;
  char *target = "a:annie:2000101" ;
  AC_DB db = 0 ;
  AC_ITER iter ;
  AC_OBJ k ;
  int nn = -1 ;

  if (argc > 1)
    target = argv[1] ;
  db = ac_open_db (target , &err) ;
  if (err)
    printf ("Error message %s", err) ;
  if (!db)
    messcrash ("could not open %s", target) ;


  iter = ac_query_iter (db, 0, "find gene smap", 0, 0) ;
  while (--nn && (k = ac_next_obj (iter)))
    {
      exportGene (db, k) ;
      ac_free (k) ;  
      /* if (ok) break ; */
    }
  ac_free (iter) ;
  ac_db_close (db) ;

  return 0 ;
}
