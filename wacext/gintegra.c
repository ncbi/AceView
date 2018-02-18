/*
 * Exportation des alignemants pour la confernce gintegra 2003 
 */

#include "../wac/ac.h"
#include <errno.h>

static int gCodingFrame (AC_TABLE coding, int b1, int r1)
{
  int a1, a2, x1, ii, frame = 0 ;

  for (ii = 0 ; coding && ii < coding->rows ; ii++)
    {
      a1 = ac_table_int (coding, ii, 0, -9999) ;
      a2 = ac_table_int (coding, ii, 1, -9999) ;
      x1 = ac_table_int (coding, ii, 3, -9999) ;
      
      
      if (b1 >= a1 && b1 <= a2 && x1 >= r1)
	{
	  frame = (b1 - a1 + x1 - r1) % 3 ;
	  frame = (3 - frame) % 3 ; /* bp i should jump */
	  break ;
	}
    }
  return frame ;
} /* gCodingFrame */

static int gExportExons (AC_OBJ gene_id, int ngene, AC_OBJ est, AC_OBJ tg, AC_OBJ chrom,
			 int g1, int g2, 
			 AC_TABLE coding, int r1, int p1, int p2, /* coding part of the mrna in gene coords */
			 int len, int nTg,
			 AC_HANDLE h, int gold)
{
  int type = 3 ; /* 1 to export to fujii
		  * 2 to mimick his export to us !
		  * 3 new ali format
		  * 4 hack to re-split the NEDO/JBIRC genes
		  */
  const char *ccp ;
  char *strand, *geneStrand ;
  int ii, ii1, jj, jjTotal, nExon, nExon2, nExon4 = 0, nHack = 0, frame, oldX1 = 0 ;
  int a1, a2, x1, x2 ; /* a single aligned exon in tg coords */
  int b1, b2 ; /* same in chrom coords */
  int q1, q2 ; /* coding part of the exon in gene coords */
  int w1, w2 ; /* coding part in chrom coords */
  AC_TABLE tt = ac_tag_table (tg, "Assembled_from", h) ;
  /* count all exons to number them in the direction of the EST */
  for (ii = jjTotal = 0 ; ii < tt->rows ; ii++)
    {
      ccp = ac_table_printable (tt, ii, 2, "") ;
      if (strcmp (ccp, ac_name (est)))
	continue ;
      jjTotal++ ;
    }

  for (ii1 = jj = 0 ; ii1 < tt->rows ; ii1++)
    {
      /* ii1 -> ii: export in the direction of the genome */
      ii = (type == 3 || type == 4 || g1 < g2) ? ii1 : tt->rows - ii1 - 1 ; 
      ccp = ac_table_printable (tt, ii, 2, "") ;
      if (strcmp (ccp, ac_name (est)))
	continue ;
      a1 = ac_table_int (tt, ii, 0, -9999) ;
      a2 = ac_table_int (tt, ii, 1, -9999) ;
      x1 = ac_table_int (tt, ii, 3, -9999) ;
      x2 = ac_table_int (tt, ii, 4, -9999) ;
      if (a1 == -9999 ||
	  a1 == -9999 ||
	  a1 == -9999 ||
	  a1 == -9999)
	messcrash ("bad values in Assembled_from tg = %s, est = %s, ii = %d", 
		   ac_name(tg), ac_name (est), ii) ;

      /* intersect the exon with the coding region
	 if (q1 > q2)  the exon is non coding */
      q1 = p1 ; q2 = p2 ; frame = 0 ;
      if (q1 < a1) { frame = jj ? 0 : gCodingFrame (coding, a1, r1) ; q1 = a1 + frame ; }
      if (q2 > a2) q2 = a2 ;

      if (g1 < g2)
	{
	  strand = "+" ;
	  b1 = g1 + a1 - 1 ;
	  b2 = g1 + a2 - 1 ;
	  w1 = g1 + q1 - 1 ;
	  w2 = g1 + q2 - 1 ;
	}
      else
	{
	  strand = "-" ;
	  b1 = g1 - a1 + 1 ;
	  b2 = g1 - a2 + 1 ;
	  w1 = g1 - q1 + 1 ;
	  w2 = g1 - q2 + 1 ;
	}


      nExon = ++jj ; 
      if ((b1 < b2 && x1 < x2) ||( b1 > b2 && x1 > x2))
	{
	  strand = "+" ;
	  nExon2 = nExon ;
	}
      else
	{
	  strand = "-" ;
	  nExon2 = jjTotal - nExon + 1 ;
	}
      if (g1 < g2)
	{
	  geneStrand = "+" ;
	}
      else
	{
	  geneStrand = "-" ;
	}
      if (q1 > q2) w1 = w2 = 0 ; /* the exon is non coding */


      if (type == 3)
	{
	  /*we expect, in the direction of the gene
	    gene = $1 ;
	    acc = $2 ;
	    method = $3 ;
	    nexon = $4 ;
	    chrom = $5 ;
	    a1 = $6 ; a2 = $7 ; x1 = $8 ; x2 = $9 ;
	  */
	  printf("%d\t", ngene) ;
	  printf("%s\t", ac_name(est)) ;
	  if (gold)
	    printf ("Gold\t") ;
	  else if (ac_has_tag (tg, "Method")) 
	    {
	      AC_TABLE gMet = ac_tag_table (tg, "Method", 0) ;
	      const char *mm = ac_table_printable (gMet, 0, 0, "AceView") ;
	      
	      printf ("%s\t", mm) ;
	      ac_free (gMet) ;
	    }
	  else
	    printf ("AceView\t") ;
	  printf("%d\t", nExon) ;
	  printf ("%s\t", ac_name(chrom)) ;
	  
	  printf("%d\t%d\t", b1,b2) ;
	  printf("%d\t%d\t", x1,x2) ;

	  printf ("\n") ;
	}
      else if (type == 4)
	{
	  /*we expect, in the direction of the gene
	    gene = $1 ;
	    acc = $2 ;
	    method = $3 ;
	    nexon = $4 ;
	    chrom = $5 ;
	    a1 = $6 ; a2 = $7 ; x1 = $8 ; x2 = $9 ;
	  */
	  if (x1 < x2 && oldX1 &&
	      x1 < oldX1 - 300)
	    { nHack++ ; nExon4 = nExon - 1 ; }
	  if (x1 > x2 && oldX1 &&
	      x1 > oldX1 + 300)
	    { nHack++ ; nExon4 = nExon - 1 ; }
	  oldX1 = x1 ;
	  printf("%s", ac_name (gene_id)) ;
	  if (nHack) printf("_%c", 'A' + nHack - 1) ;
	  printf("\t") ;
	  printf("%s\t", ac_name(est)) ;
	  if (! gold)
	    printf ("AceView\t") ;
	  else
	    printf ("Gold\t") ;
	  printf("%d\t", nExon - nExon4) ;
	  printf ("%s\t", ac_name(chrom)) ;
	  
	  printf("%d\t%d\t", b1,b2) ;
	  printf("%d\t%d\t", x1,x2) ;

	  printf ("\n") ;
	}
      else
	{
	  /* always export in a stupid way */
	  if (x1 > x2)
	    { int xx = x1 ; x1 = x2 ; x2 = xx ; }
	  if (w1 > w2)
	    { int ww = w1 ; w1 = w2 ; w2 = ww ; }
	  if (b1 > b2)
	    { int bb = b1 ; b1 = b2 ; b2 = bb ; }

	  if (type==1) printf ("%s\t", ac_name(gene_id)) ;
	  printf("%s\t", ac_name(est)) ;
	  if (type==2) 
	    {
	      if (! gold)
		printf ("A\t") ;
	      else
		printf ("Gold\t") ;
	    }
	  printf("%s\t", strand) ;
	  printf("%d\t%d\t", b1,b2) ;
	  printf("%d\t%d\t", w1,w2) ;
	  printf("%d\t", nExon2) ;
	  printf("%d\t%d\t%d\t", x1,x2, len) ;
	  printf ("%s\t", ac_name(chrom)) ;
	  if (type==1) printf ("%s\t", geneStrand) ;
	  printf ("\n") ;
	}
    }    
  return 1 ;
} /* gExportExons */

static int gExportEstTgMrna (AC_OBJ est, AC_OBJ tg, AC_TABLE coding, int r1, int p1, int p2, int nTg, AC_HANDLE h, int gold)
{
  AC_OBJ gene = 0, gene_id = 0, chrom = 0 ;
  int g1, g2 ; /* coords of the gene on the chromosome */
  int len ;
  AC_TABLE tt = 0 ;
  static int ngene = 0 ;
  /* get the intmap of the gene */
  tt = ac_tag_table (tg, "IntMap", h) ;
  if (!tt || ! tt->rows || tt->cols < 3)
    messcrash ("missing IntMap in transcribed_gene %s", ac_name (tg)) ;
  chrom = ac_table_obj (tt, 0, 0, h) ;
  g1 = ac_table_int (tt, 0, 1, 0) ;
  g2 = ac_table_int (tt, 0, 2, 0) ;
  

  /* get the est length */
  tt = ac_tag_table (est, "DNA", h) ;
  if (!tt || ! tt->rows || tt->cols < 2)
    messcrash ("missing DNA in est %s", ac_name (est)) ;
  len = ac_table_int (tt, 0, 1, 0) ;

  /* get the gene_id */
  gene = ac_tag_obj (tg, "Gene", h) ;
  gene_id = gene ? ac_tag_obj (gene, "Gene_id", h) : 0 ;
  if (!gene_id) gene_id = tg ;
  /* get the coding region in gene coordinates */
  gExportExons (gene_id, ++ngene, est, tg, chrom, g1, g2, coding, r1, p1, p2, len, nTg, h, gold) ;
  
  return 1 ;
} /* gExportEstTgMrna */

static AC_TABLE gCodingRegion (AC_OBJ mrna, int *r1p, int *p1p, int *p2p, AC_HANDLE h)
{
  AC_TABLE products = ac_tag_table (mrna, "Product", h) ;
  AC_TABLE tt = ac_tag_table (mrna, "Coding", h) ;
  int ii, x1, x2, y1, y2 ; ;
  
  /* get the longest CDS in spliced mrna coords */
  x1 = x2 = 1 ;
  for (ii = 0 ; ii < products->rows ; ii++)
    {
      y1 = ac_table_int (products, ii, 1, 0) ;
      y2 = ac_table_int (products, ii, 2, 0) ;
      if (y2 - y1 > x2 - x1)
	{ x1 = y1 ; x2 = y2 ; }
    }
  /* find the corresponding unspliced coords , they must exist in Coding */
  *p1p = *p2p = 1 ; *r1p = x1 ;
  for (ii = 0 ; ii < tt->rows ; ii++)
    {
      y1 = ac_table_int (tt, ii, 2, -9999) ;
      if (y1 == x1) *p1p = ac_table_int (tt, ii, 0, 1) ;
      y2 = ac_table_int (tt, ii, 3, -9999) ;
      if (y2 == x2) *p2p = ac_table_int (tt, ii, 1, 1) ; 
    }
  
  return tt ;
} /* gCodingRegion */

static int gExportEstTg (AC_OBJ est, AC_OBJ tg, int nTg, AC_HANDLE h, int gold)
{
  int ii, jj ;
  int m1 ; /* start of mrna in tg coords */
  int r1, p1, p2 ; /* coding region in mrna coords */
  AC_OBJ mrna, tgMrna ;
  AC_TABLE coding = 0, tgMrnas = ac_tag_table (tg, "Mrna", h) ;
  AC_TABLE estMrnas = ac_tag_table (est, "In_mrna", h) ;

  if (tgMrnas && estMrnas)
    for (ii = 0 ; estMrnas && ii < estMrnas->rows ; ii++)
      {
	mrna = ac_table_obj (estMrnas, ii, 0, h) ;
	for (jj = 0 ; jj < tgMrnas->rows ; jj++)
	  {
	    tgMrna = ac_table_obj (tgMrnas, jj, 0, h) ;
	    if (!strcmp (ac_name(mrna), ac_name (tgMrna)))
	      {
		m1 = ac_table_int (tgMrnas, jj, 1, 1) ;
		coding = gCodingRegion (mrna, &r1, &p1, &p2, h) ;
		/* move to gene coords */
		return gExportEstTgMrna (est, tg, coding, r1, m1 + p1 - 1, m1 + p2 - 1, nTg, h, gold) ;
	      }
	  }
      } 
  else
    return gExportEstTgMrna (est, tg, 0, 0, 0, 0, nTg, h, gold) ;
    
  return 0 ;
} /* gExportEstTg */


static int gExportEst (AC_OBJ est, char *template, int gold)
{
  int ii, nn = 0 ;
  AC_HANDLE h = handleCreate () ;
  AC_TABLE tgs = ac_tag_table (est, "From_gene", h) ;
  AC_OBJ tg ;

  for (ii = 0 ; ii < tgs->rows ; ii++)
    {
      tg = ac_table_obj (tgs, ii, 0, h) ;

      if (template && ! ac_has_tag (tg, template))
	continue ;
      if (gold && ! ac_has_tag (tg, "gold"))
	continue ;
      if (gold == 1 && ! ac_has_tag (tg, "excellent") && ! ac_has_tag (tg, "good"))
	continue ;
      if (gold == 2&& (ac_has_tag (tg, "excellent") ||  ac_has_tag (tg, "good")))
	continue ;
      nn += gExportEstTg (est, tg, ii, h, gold) ;
    }
  
  ac_free (h) ;
  return nn ;
} /* gExportEst */



static int gExport (AC_DB db, char *template, int gold)
{
  AC_ITER iter = 0 ; 
  AC_OBJ est ;
  int nn = 0, nest = 0 ;

  if (template && !strcasecmp(template,"is_mrna"))
    {
      iter = ac_query_iter (db, TRUE, "find est is_mrna && from_gene", 0, 0) ; /* is_mrna alignments */
      template = 0 ;
    }
  else if (template && !strcasecmp(template,"is_nm"))
    {
      iter = ac_query_iter (db, TRUE, "find est IS NM* && from_gene", 0, 0) ; /* is_mrna alignments */
      template = 0 ;
    }
  else if (template) /* select some method */
    iter = ac_query_iter (db, TRUE, messprintf ("find tg ; %s; >read ;  from_gene", template), 0, 0) ;
  else
    iter = ac_query_iter (db, TRUE, "find est  from_gene", 0, 0) ; /* all alignments */
  while ((est = ac_next_obj (iter)))
    {
      nn += gExportEst (est, template, gold) ;
      ac_free (est) ;
      nest++ ;
    }
  ac_free (iter) ;
  fprintf (stderr, "// Considered %d ests\n", nest) ;

  return nn ;
} /* gExport */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: gintegra ACEDB [-is_mrna || -is_nm || -gold | -gold2 | AceView | ...]\n") ;
  fprintf (stderr, "// Example:  gintegra AceView >! aceview.ali \n") ;
  fprintf (stderr, "//    -is_nm : export NM alignments\n") ;
  fprintf (stderr, "//    -is_mrna : export \"est is_mrna\" alignments\n") ;
  fprintf (stderr, "//    -gold : export excellent and good GOLD alignments\n") ;
  fprintf (stderr, "//    -gold2 : export other GOLD alignments\n") ;
  fprintf (stderr, "//    AceView: export AceView or any other method alignments\n") ;
  exit (1) ;
}

/*************************************************************************************/

int main (int argc, char **argv)
{
  int nn = 0 ;
  const char *s = "ok" ;
  int gold = 0 ;
  char *dbName = argc>=2 ? argv[1] : 0 ; 
  char *template = argc>=3 ? argv[2] : 0 ;
  AC_DB db ;

  if (!dbName)
    usage () ;

  db = ac_open_db (dbName, &s);
  if (!db)
    messcrash ("Failed to open db %s, error %s", dbName, s) ;

  if (template && !strcmp (template, "-gold"))
    {
      gold = 1 ;
      template = argc>=4 ? argv[3] : 0 ;
    }
  if (template && !strcmp (template, "-gold2"))
    {
      gold = 2 ;
      template = argc>=4 ? argv[3] : 0 ;
    }
  nn = gExport (db, template, gold) ;

  ac_db_close (db) ;

  fprintf (stderr, "Exported %d alignments\n", nn) ;
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
