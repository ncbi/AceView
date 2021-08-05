#include "regular.h"
#include "dict.h"
/*GFF (GTF) Reader. Takes_ GFF transcript information from a file and outputs in ace format to std out in new 'mrna' ace schema.
 * This file was fixed and adated for the Encode project
 *
 * NOTE: GTF spec does NOT define order of lines,
 * but AceDB format mrna objects need an origin to define
 * positions of exons and coding regions relative to.
 *
 * 2021: for several years thgis code has been superseeded by dna2dna -gff
 * since it depends on RPC, it is only compiled if one ask
     make rpcace
 * but is no longer implied by
     make all     
 */

/* 
   ENr111  VEGA_Novel_CDS  exon    -9347   -9177   .       +       .       transcript_id "RP11-90M5.1-001"; gene_id "RP11-90M5.1"
   ; gene_alias "."; exon_id "RP11-90M5.1-001-exon1";
   ENr111  VEGA_Novel_CDS  intron  -9176   -8677   .       +       .       transcript_id "RP11-90M5.1-001"; gene_id "RP11-90M5.1"
   ; gene_alias "."; intron_id "RP11-90M5.1-001-intron1"; intron_status "not_tested";
*/
typedef struct gffStruct {
  int cosmid, method, etype, a1, a2, x1, x2, m1, m2, dot, forward, tr_type, tr, gene_type, gene ; 
} GFF ;

DICT *dict = 0 ;

/***************************************************************/

static int gffOrder (const void *va, const void *vb)
{
  const GFF *a = (const GFF *) va;
  const GFF *b = (const GFF *) vb;
  int nn ;

  nn = a->forward - b->forward ;
  if (nn) return nn ;

  nn = a->gene - b->gene ;
  if (nn) return nn ;

  nn = a->tr - b->tr ;
  if (nn) return nn ;

  nn = a->a1 - b->a1 ;
  if (nn) return a->forward ? nn : -nn ;

  nn = a->a2 - b->a2 ;
  if (nn) return a->forward ? nn : -nn ;

  return 0 ;
}

/***************************************************************/

static void gffError (int ln, char *txt)
{
  char *cq = freepos () ;
  
  if (! cq || ! *cq) cq = "no data" ;
  fprintf (stderr, "// Error line %d: %s at %s\n", ln, txt, cq) ;
  exit (1) ;
}

/***************************************************************/

static void gffParse (int level, int *nlp, int *ngp, int *nmp)
{
  int ii = 0, jj, exon, cds ;
  int g1, m1, m2, oldGene = 0, oldTr = 0 ;
  Array aa = arrayCreate (10000, GFF) ;
  GFF *up ;
  char *cp ;

  dictAdd (dict, "toto", &ii) ;
  dictAdd (dict, "exon", &exon) ;
  dictAdd (dict, "cds", &cds) ;

  /* parse the data */
  while (freecard (level))
    {
      (*nlp)++ ;
      cp = freeword () ;
      if (! cp)
	break ;
      up = arrayp (aa, ii++, GFF) ; /* new line */
      dictAdd (dict, cp, &(up->cosmid)) ;

      cp = freeword () ;
      if (! cp)
	gffError (*nlp, "no method") ;
      dictAdd (dict, cp, &(up->method)) ;
      
      cp = freeword () ;
      if (! cp)
	gffError (*nlp, "no etype") ;
      dictAdd (dict, cp, &(up->etype)) ;
      
      if (! freeint (&(up->a1)))
	gffError (*nlp, "no a1") ;
      if (! freeint (&(up->a2)))
	gffError (*nlp, "no a2") ;
      
      cp = freeword () ;
      if (! cp || *cp != '.')
	gffError (*nlp, "no dot") ;

      cp = freeword () ;
      if (! strcmp (cp,"+"))
	up->forward = TRUE ;
      else if (! strcmp (cp,"-"))
	up->forward = FALSE ;
      else
	gffError (*nlp, "unknown strand") ;

      cp = freeword () ;
      if (! cp)
	gffError (*nlp, "no second dot") ;

      cp = freeword () ;
      if (! cp)
	gffError (*nlp, "no transcript_id") ;
      if (strcmp (cp,"transcript_id"))
	gffError (*nlp, "unknown transcript_id") ;

      cp = freeword () ;
      if (! cp)
	gffError (*nlp, "no tr") ;
      dictAdd (dict, cp, &(up->tr)) ;

      cp = freeword () ;
      if (! cp || *cp != ';')
	gffError (*nlp, "no semi column") ;

      cp = freeword () ;
      if (! cp)
	gffError (*nlp, "no gene_id") ;
      if (strcmp (cp,"gene_id"))
	gffError (*nlp, "unknown gene_id") ;

      cp = freeword () ;
      if (! cp)
	gffError (*nlp, "no gene") ;
      dictAdd (dict, cp, &(up->gene)) ;
    }
  
  arraySort (aa, gffOrder) ;
  /* fix strand */
  for (ii = 0, up = arrp (aa, 0, GFF) ; ii < arrayMax (aa) ; ii++, up++)
    if (! up->forward)
      { jj = up->a1 ; up->a1 = up->a2 ; up->a2 = jj ; }

  /* make coords relative */
  for (ii = 0, up = arrp (aa, 0, GFF) ; ii < arrayMax (aa) ; ii++, up++)
    {
      if (up->gene != oldGene)
	{
	  (*ngp)++ ;
	  g1 = up->a1 ; /* g2 = up->a2 ;  */
	  oldGene = up->gene ;
	}
      if (up->tr != oldTr)
	{
	  (*nmp)++ ;
	  m1 = up->a1 ; m2 = up->a2 ; 
	  oldTr = up->tr ;
	}
      if (up->forward)
	{
	  up->x1 = up->a1 - g1 + 1 ;
	  up->x2 = up->a2 - g1 + 1 ;
	  up->m1 = up->a1 - m1 + 1 ;
	  up->m2 = up->a2 - m1 + 1 ;
	}
      else
	{
	  up->x1 = g1 - up->a1 + 1 ;
	  up->x2 = g1 - up->a2 + 1 ;
	  up->m1 = m1 - up->a1 + 1 ;
	  up->m2 = m1 - up->a2 + 1 ;
	}
    }

  /* export */
  arrayp (aa, arrayMax (aa), GFF) ; /* add an empty line */
  oldTr = 0 ;
  for (ii = 0, up = arrp (aa, 0, GFF) ; ii < arrayMax (aa) ; ii++, up++)
    {
      if (up->tr && up->etype != exon)
	continue ;
      if (oldTr != up->tr)
	{
	  if (oldTr)
	    {
	      printf ("\nSequence \"%s\"\nSubsequence \"%s\" %d %d\n\n"
		      , dictName (dict, (up-1)->cosmid)
		      , dictName (dict, oldTr)
		      , m1, m2) ;
	    }
	  if (up->tr)
	    {
	      printf ("Sequence \"%s\"\n"
		      , dictName (dict, up->tr)) ;
	      printf ("CDS\nIs_predicted_gene\nMethod %s\n"
		      , dictName (dict, up->method)
		      ) ;
	    }
	  oldTr = up->tr ;
	  m1 = up->a1 ; m2 = up->a2 ;
	}
      if (up->tr && up->etype == exon)
	{
	  printf ("Source_exons %d %d\n", up->m1, up->m2) ;
	  if (up->forward && m2 < up->a2) m2 = up->a2 ;
	  if (!up->forward && m2 > up->a2) m2 = up->a2 ;
	}
    }

  /* start again and export the CDS stuff */
  oldTr = 0 ;
  for (ii = 0, up = arrp (aa, 0, GFF) ; ii < arrayMax (aa) ; ii++, up++)
    {
      if (oldTr != up->tr)
	{
	  if (oldTr)
	    {
	      printf ("\nSequence \"%s\"\nSubsequence \"CDS_%s\" %d %d\n\n"
		      , dictName (dict, (up-1)->cosmid)
		      , dictName (dict, oldTr)
		      , m1, m2) ;
	    }
	  if (up->tr)
	    {
	      printf ("Sequence \"CDS_%s\"\n"
		      , dictName (dict, up->tr)) ;
	      printf ("CDS\nIs_predicted_gene\nMethod CDS_%s\n"
		      , dictName (dict, up->method)
		      ) ;
	    }
	  oldTr = up->tr ;
	  m1 = up->a1 ; m2 = up->a2 ;
	}
      if (up->tr && up->etype == cds)
	{
	  printf ("Source_exons %d %d\n", up->m1, up->m2) ;
	  if (up->forward && m2 < up->a2) m2 = up->a2 ;
	  if (!up->forward && m2 > up->a2) m2 = up->a2 ;
	}
    }
} /* gffParse */
 
/***************************************************************/
/***************************************************************/

static void usage (void)
{
  char * usage = "gffread file : adapted for Encode project april 2005";

  printf ("Usage: %s\n", usage);
  exit(1);
} /* usage */

/***************************************************************/

int main (int argc, char * argv[])
{
  char *filNam ;
  FILE *fil = 0 ;
  int  level, nlines = 0, nmrnas= 0, ngenes = 0 ;
  
  freeinit () ;
  dict = dictCreate (10000) ;

  if (argc > 1) 
    filNam = argv[1];
  else 
    usage () ;

  if (! (fil = filopen (filNam, 0, "r")))
    usage () ;
  
  level = freesetfile (fil, 0) ;
  gffParse (level, &nlines, &ngenes, &nmrnas) ;
  freeclose (level) ;
  filclose (fil) ;

  printf ("// done, read %d lines, %d genes, %d mrna\n"
	  , nlines, ngenes, nmrnas) ;
  return 0 ;
} /* main */

/***************************************************************/
/***************************************************************/


