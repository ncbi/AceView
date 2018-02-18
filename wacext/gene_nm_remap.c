#include "../wac/ac.h"

typedef struct { char nam[256] ; int x1, x2 ;} GG ;

static void exportGene (AC_DB db, AC_OBJ gene)
{
  AC_HANDLE h = handleCreate () ;
  AC_TABLE nms, mrnas, pgs, imap, exons ;
  AC_OBJ chrom = 0, chrommap = 0, tg, mrna, tmp ;
  int ir, jr, jr2, a1, a2, len ;
  BOOL isTg, isDown = TRUE ;
  Array aa = 0 ;
  GG *gg ;
  const char *ccp ;

  /* transform the gene into a tg, add its nms as clones */
  nms = ac_tag_table (gene, "NM_Id", h) ;
  if (!nms)
    return ;
  
  aa = arrayCreate (32, GG) ;
  printf ("Transcribed_gene \"tg_%s\"\n", ac_name(gene)) ;
  for (ir = 0 ; ir < nms->rows ; ir++)
    printf ("cDNA_clone \"GB:%s\"\n", ac_table_printable (nms, ir, 0, "-")) ;
  printf ("\n") ;

  printf ("Gene \"%s\"\n", ac_name(gene)) ;
  for (ir = 0 ; ir < nms->rows ; ir++)
    printf ("NM_Id \"%s\"\n", ac_table_printable (nms, ir, 0, "-")) ;
  printf ("\n") ;

      /* find the size of the gene->mrnas or the gene->pgs */
  isTg = ac_has_tag (gene, "Transcribed_gene") ;
  a1 = a2 = -1 ;
  if (isTg)
    {
      tg = ac_tag_obj (gene, "Transcribed_gene", h) ;
      mrnas = ac_tag_table (tg, "mRNA", h) ;
      for (jr = 0 ; jr < mrnas->rows ; jr++)
	{
	  mrna = ac_table_obj (mrnas, jr, 0, h) ;
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
	  exons =  ac_tag_table (mrna, "Splicing", 0) ;
	  printf ("Sequence \"%s\"\nIs_predicted_gene\nModel_of_gene \"%s\"\nABI\nMethod curated\n"
		  ,ac_name(mrna), ac_name(gene)) ;
	  for (jr2 = 0 ; exons && jr2 < exons->rows ; jr2++)
	    {
	      if (exons->cols >= 5 &&
		  strstr (ac_table_printable (exons, jr2, 4, "-"), "xon"))
		printf("Source_exons %d %d\n"
		       , ac_table_int (exons, jr2, 0, 0) 
		       , ac_table_int (exons, jr2, 1, 0) 
		       ) ;
	    }
	  ac_free (exons) ;
	  if ((tmp = ac_tag_obj (mrna, "NM_Id", h)))
	    {
	      printf("NM_id \"%s\"\n", ac_name(tmp)) ;
	      ac_free (tmp) ;
	    }
	  printf ("\n") ;
	}
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
	  exons =  ac_tag_table (mrna, "Source_exons", 0) ;
	  printf ("Sequence \"%s\"\nIs_predicted_gene\nModel_of_gene \"%s\"\nCDS\nMethod curated\n"
		  , ac_name(mrna), ac_name (gene)) ;
	  for (jr2 = 0 ; jr2 < exons->rows ; jr2++)
	    printf("Source_exons %d %d\n"
		   , ac_table_int (exons, jr2, 0, 0) 
		   , ac_table_int (exons, jr2, 1, 0) 
		   ) ;
	  ac_free (exons) ;
	  if ((tmp = ac_tag_obj (mrna, "NM_Id", h)))
	    {
	      printf("NM_id \"%s\"\n", ac_name(tmp)) ;
	      ac_free (tmp) ;
	    }
	  printf ("\n") ;
	}
    }

  if (arrayMax (aa))
    {
      gg = arrayp (aa, 0, GG) ;
      a1 = gg->x1 ; a2 = gg->x2 ;
      if (a1 < a2) isDown = TRUE ;
      else isDown = FALSE ;      
    }
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
  len = isDown ? a2 - a1 + 1 : a1 - a2 + 1 ;

  /* create a pseudo sequence .1 kb up and down of the gene */
  printf ("Sequence \"c_%s\"\nGenomic\nTranscribed_gene \"tg_%s\" 101 %d\n"
	  , ac_name(gene)
	  , ac_name(gene)
	  , len + 100) ;

  for (ir = 0 ;  aa && ir < arrayMax (aa) ; ir++)
    {
      gg = arrayp (aa, ir, GG) ;
      /* export the new positions */
      if (isDown)
	printf ("Subsequence %s %d %d\n"
		, gg->nam
		, gg->x1 - a1 + 101 
		, gg->x2 - a1 + 101
		) ; 
      else
	printf ("Subsequence %s %d %d\n"
		, gg->nam
		, a1 - gg->x1 + 101 
		, a1 - gg->x2 + 101
		) ; 
    }
  printf ("\n") ;

  imap = ac_tag_table (gene, "IntMap", 0) ;
  if (imap && imap->cols >= 2)
    {
      if (!chrommap) 
	chrommap = ac_table_obj (imap, 0, 0, h) ;
      gg = arrayp (aa, arrayMax (aa) , GG) ;
      gg->x1 = ac_table_int (imap, 0, 1, 0) ;
      gg->x2 = ac_table_int (imap, 0, 2, 0) ;
      printf ("Sequence \"c_%s\"\n"
	      ,ac_name(gene)) ;
      if (isDown)
	printf ("Genes %s %d %d\n"
		, ac_name(gene)
		, gg->x1 - a1 + 101 
		, gg->x2 - a1 + 101
		) ; 
      else
	printf ("Genes %s %d %d\n"
		, ac_name(gene)
		, a1 - gg->x1 + 101 
		, a1 - gg->x2 + 101
		) ; 
      printf ("\n") ;
      if ((ccp = ac_tag_printable (gene, "Locus", 0)))
	printf ("Gene \"%s\"\nLocus \"%s\"\n\n",
		ac_name(gene), ccp) ;
      if ((ccp = ac_tag_printable (gene, "LocusId", 0)))
	printf ("Gene \"%s\"\nLocusId \"%s\"\n\n",
		ac_name(gene), ccp ) ;
      if ((ccp = ac_tag_printable (gene, "NewName", 0)))
	printf ("Gene \"%s\"\nNewName \"%s\"\n\n",
		ac_name(gene), ccp) ;
    }
  ac_free (imap) ;

  /* construct the DNA and export it */
  if (chrommap && 
      (chrom = ac_get_obj (db, "Sequence", ac_name(chrommap), h)))
    {
      char *dna = 0 ;
      if (isDown)
	dna = ac_zone_dna (chrom, a1 - 100, a2 + 100, h) ;
      else
	dna = ac_zone_dna (chrom, a1 + 100, a2 - 100, h) ;
      if (dna)
	printf ("DNA \"c_%s\"\n%s\n\n", ac_name(gene), dna) ;
      else 
	printf ("\n\n//ERROR no dna for gene %s\n", ac_name(gene)) ;
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
    }
  ac_free (iter) ;
  ac_db_close (db) ;

  return 0 ;
}
