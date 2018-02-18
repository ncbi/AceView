/*
 * Export from aceview a separation of the genome
 * into : intergenic/intronic/exonic regions
 *
 * for each region defined on the genome
 * associate a keyset of genes
 */


#include "../wac/ac.h"
#include "keyset.h"
#include <errno.h>
#include <math.h>
#include "keyset.h"
#include "freeout.h"
#include "vtxt.h"
#include "dict.h"
#include "aceio.h"

typedef enum {ZERO, INTERGENIC, INTRON, EXON, CODINGEXON } zTYPE ;
char *zType[] = {"ZERO", "Intergenic", "Premessenger", "Exonic", "Coding" }  ;

typedef struct worfStruct { AC_DB db ; Array aa, bb ; const char *html, *byGene, *ace, *template; DICT *dict ; AC_HANDLE h ;} NW ;
typedef struct worfStructZ { 
  AC_OBJ tg   , mrna   , product   , clo    , est ; 
  KEY tgKey, mrnaKey, productKey, bestProductKey, cloKey ;
  int a1, a2, x1, x2, err, p1, p2, bestP1, bestP2, stop, coding, ecoding, nonbestprotein, inframe,
    contains5pStop, contains3pStop, atg, partial5p, atg2, mrnaLength
    , partial3p , missing3pRead, missing5pRead, paire
    , insert, gap, endsConfirmed, pGap, score, anomaly, rearray ; } NWZ ;

      /* 0 = new peptide
       * 1 = complete protein, fully sequenced,  atg||atg2 > 0 && stop > 0 && !partial5p > 90bp
       * 2 = complete protein, partial sequence,  atg||atg2 && stop && gap  && !partial5p > 90bp
       * 3 = 5' partial protein, fully sequenced    !(atg||atg2 && stop) && !gap
       * 4 = 3' partial protein, fully sequenced    !(atg||atg2 && stop) && !gap
       * 5 = 5' and 3' partial protein, fully sequenced    !(atg||atg2 && stop) && !gap
       * 6 = 5' partial protein, partial sequence,  !(atg||atg2 && stop) && gap
       * 7 = 3' partial protein, partial sequence,  !(atg||atg2 && stop) && gap
       * 8 = 5' and 3' partial protein, partial sequence,  !(atg||atg2 && stop) && gap
       * 9 = 5' complete protein, 3' end  not known,   paire != 0x3 correct corresponding end
       * 10 = 3' complete protein, 5' end  not known,   paire != 0x3 correct corresponding end
       * 11 = 5' partial and 3' end not known,   paire != 0x1
       * 12 = 3' partial and 5' end not known,   paire == 0x2
       * 13 = contains the stop of the protein
       * 14 = contains a stop before the ATG
       */
static  char * scoreTitle[] = {
  "New peptide"
  , "Complete protein, fully sequenced" 
  , "Complete protein, partial sequence"
  , "5' partial protein, fully sequenced"
  , "3' partial protein, fully sequenced"
  , "5' and 3' partial protein, fully sequenced "
  , "5' partial protein, partial sequence"
  , "3' partial protein, partial sequence"
  , "5' and 3' partial protein, partial sequence"
  , "5' complete protein, 3' end  not known"
  , "3' complete protein, 5' end  not known"
  , "5' partial and 3' end not known"
  , "3' partial and 5' end not known"
  , "contains the stop of the protein"
  , "contains a stop before the ATG"
}  ;

/*************************************************************************************/
/* debug function */
static void nwShow (Array aa)
{
  NWZ *up ;
  int ii ;
  
  if (aa)
    for (ii = 0 ; ii < arrayMax(aa) ; ii++)
      {
	up = arrp (aa, ii, NWZ) ;
	fprintf (stderr, "%3d  %s  %s p1/p2=%d/%d\ta1/a2=%d/%d\n"
		 , ii
		 , ac_key_name(up->cloKey), ac_key_name(up->productKey)
		 , up->p1, up->p2, up->a1, up->a2
		 ) ;
      }
  return ;
} /* nwShow  */

/*************************************************************************************/

static int nwCodingLength (AC_OBJ mrna, int a1, AC_HANDLE h) 
{
  char *dna = ac_obj_dna (mrna, h) ;
  char *cp ;
  int i, imax, n = 0 ;

  while (a1 < 1) { a1 += 3,  n+= 3 ; }
  if (!dna) 
    return 0 ;
  imax = strlen (dna) ;

  for (i = a1 - 1, cp = dna + a1 - 1 ; i < imax - 2 ; i += 3 , cp += 3, n += 3)
    if (!strncmp (cp, "taa", 3) || !strncmp (cp, "tag", 3) || !strncmp (cp, "tga", 3))
      break ;
  return n ;
}

/*************************************************************************************/

static void nwReportClo (Array aa, Array bb, AC_HANDLE h) 
{ 
  int ii ;
  NWZ *up, *vp ;
  char *conf[] = {"None", "5'", "3'", "Both"} ;
  char *frame[] = {"No", "1", "2", "3", "NA"} ;
  
  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      up = arrp (aa, ii, NWZ) ;
      if (! up->clo)
	continue ;
      
      
      if (up->product)
	{
	  if (up->p1 > 3 || ac_has_tag (up->product, "NH2_Complete"))
	    up->endsConfirmed |= 0x1 ;
	  if (up->p2 < up->mrnaLength - 3 || ac_has_tag (up->product, "COOH_Complete"))
	    up->endsConfirmed |= 0x2 ;

	  /* did we clip not as the vector cliping or have a manual annotation */
	  /* is the sequence accurate ? doubtful we can do that */
	  
	  if (up->contains5pStop)
	    up->score = 14 ;
	  else if (up->contains3pStop)
	    up->score = 13 ;
	  else if (up->paire == 0x2)
	    {
	      if (up->partial3p)
		up->score = 12 ;
	      else
		up->score = 10 ; 
	    }
	  else if (up->paire == 0x1 && (!up->stop || up->stop == 4))
	    {
	      if (!up->atg && (!up->atg2 || up->partial5p > 90) && up->a1 > up->p1 + 3)
		up->score = 11 ;
	      else
		up->score = 9 ; 
	    }
	  else if (up->partial3p && (!up->atg && (!up->atg2 || up->partial5p > 90)) && up->a1 > up->p1 + 3)
	    {
	      if (up->gap)
		up->score = 8 ;
	      else
		up->score = 5 ; 
	    }
	  else if (up->partial3p)
	    {
	      if (up->gap)
		up->score = 7 ;
	      else
		up->score = 4 ; 
	    }
	  else if (!up->atg && (!up->atg2 || up->partial5p > 90) && up->a1 > up->p1 + 3)
	    {
	      if (up->gap)
		up->score = 6 ;
	      else
		up->score = 3 ; 
	    }
	  else
	    {
	      if (up->gap)
		up->score = 2 ;
	      else
		up->score = 1 ; 
	    }
	}


      vp = arrayp (bb , arrayMax (bb), NWZ) ;
      *vp = *up ;

      if (1)
	{
	  printf ("%s", ac_name(up->clo)) ;
	  printf ("\t%s", scoreTitle[up->score]) ;
	  printf ("\t%s", conf[up->paire]) ;
	  if (up->gap) 
	    printf ("\t~ %d", up->insert) ;
	  else if (up->paire != 0x3) 
	    printf ("\t> %d", up->insert) ;
	  else
	    printf ("\t%d", up->insert) ;
	  printf ("\t%d", up->coding) ;
	  if (up->err >= 99999)
	    printf ("\tNA") ;
	  else
	    printf ("\t%d", up->err) ;
	  if (up->paire == 0x3)
	    printf ("\t%d", up->gap) ;
	  else
	    printf ("\tNA") ;

	  printf ("\t%d\t%d", up->missing5pRead, up->missing3pRead) ;
	  printf ("\t%s", ac_name(up->mrna)) ;
	  
	  if (!  up->contains5pStop && ! (up->paire & 0x1) && (up->atg == 0 || up->atg == 4))
	    printf ("\tNA") ;
	  else
	    printf ("\t%d",  up->contains5pStop) ;
	  if (!  up->contains3pStop && ! (up->paire & 0x2) && (up->stop == 0 || up->stop == 4))
	    printf ("\tNA") ;
	  else
	    printf ("\t%d",  up->contains3pStop) ;

	  printf ("\t%s\t%s", frame[up->atg], frame[up->atg2]) ;

	  if (!(up->paire & 0x1) && ! (up->paire & 0x1) && (up->atg == 0 || up->atg == 4))
	    printf ("\tNA") ;
	  else
	    printf ("\t%d", up->partial5p) ;

	  if (!(up->paire & 0x2) && ! (up->paire & 0x2) && (up->stop == 0 || up->stop == 4))
	    printf ("\tNA") ;
	  else
	    printf ("\t%d", up->partial3p) ;

	  printf ("\t%s", frame[up->stop]) ;

	  printf ("\t%s", conf[up->endsConfirmed] ) ;

	  printf ("\n") ;
	}
    }
} /* nwReportClo */

/*************************************************************************************/

static int cloneOrder (const void *a, const void *b)
{
  int n ;
  const NWZ *up = (const NWZ *)a, *vp = (const NWZ *)b ;

  n = ac_obj_cmp (up->clo, vp->clo) ;
  if (n) return n ;
  n = ac_obj_cmp (up->mrna, vp->mrna) ;
  if (n) return n ;
  n = up->nonbestprotein - up->nonbestprotein ;
  if (n) return -n ; /* best first */
  n = up->productKey - vp->productKey ;
  if (n) return -n ; /* zero-product last */
  n = up->paire - vp->paire ; /* 5' first */
  if (n) return n ;
  n = ac_obj_cmp (up->est, vp->est) ;
  if (n) return n ;
  return 0 ;
} /* nwOrder */

/*************************************************************************************/
/* pass = 0 annotate the best protein */
/* pass = 1 annotate the non best if we touch it */
/* pass = 2 annotate new prots */

static void nwGetMrnaData (NW *nw, AC_OBJ mrna, AC_OBJ clo, int err, Array aa, AC_HANDLE h)
{
  BOOL ok ;
  int frame, jj, ir, jr, a1, a2, b1, b2, x, x1, x2, p1, p2, p21 ;
  const char *ccp ;
  AC_TABLE tbl, tbl1, products ;
  AC_OBJ tg, worf = 0, product = 0, myclo = 0 ;
  NWZ *up = 0 ;
  int mrnaLength = 0 ;
  vTXT anoTxt = vtxtHandleCreate (h) ;
  
  product = 0 ; p1 = p2 = 0 ;
  products = ac_tag_table (mrna, "Product", h) ;
  /* collect the worf reads */
  jj = 0 ;
  tg = ac_tag_obj (mrna, "From_gene", h) ;
  tbl = ac_tag_table (mrna, "Constructed_from", h) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      jj = ac_table_int (tbl, ir, 1, 0) ;
      if (mrnaLength < jj) 
	mrnaLength = jj ;
    }
  jj = arrayMax (aa) ;
  up = arrayp (aa, jj++, NWZ) ;
  up->clo = clo ; up->tg = tg ; up->mrna = mrna ; 
  /* stabilize the data */
  up->tgKey = ac_obj_key (up->tg) ;
  up->mrnaKey = ac_obj_key (up->mrna) ;
  up->cloKey = ac_obj_key (up->clo) ;
  up->mrnaLength = mrnaLength ;
  up->pGap = ac_tag_int (mrna, "Gap_length", 0) ;
  up->err = err ;

  up->a1 = b1 = -9999 ; b2 = 9999999 ;
  /* find the ends of the clo */
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      worf = ac_table_obj (tbl, ir, 2, h) ;
      myclo = ac_tag_obj (worf, "cDNA_clone", h) ;
      if (!ac_obj_equal (clo, myclo))
	continue ;
      a1 = ac_table_int (tbl, ir, 0, 0) ;
      a2 = ac_table_int (tbl, ir, 1, 0) ;
      x1 = ac_table_int (tbl, ir, 3, 0) ;
      x2 = ac_table_int (tbl, ir, 4, 0) ;
      if (x1 < x2)
	up->paire |= 0x1 ;
      else
	up->paire |= 0x2 ;
      if (up->a1 == -9999)
	{ up->a1 = a1 ; up->a2 = a2 ; }
      if (up->a1 > a1) up->a1 = a1 ;
      if (up->a2 < a2) up->a2 = a2 ;
      if (x1 < x2 && b1 < a2) b1 = a2 ;
      if (x1 > x2 && b2 > a1) b2 = a1 ;
    }

  up->insert = up->a2 - up->a1 + 1 ;
  up->coding = nwCodingLength (up->mrna, up->a1 - 1, h) ;
  if (b1 > 0 && b1 < b2 && b2 <  9999999) up->gap = b2 - b1 - 1 ;
  ok = FALSE ; product = 0 ;

  /* search for a product in frame with the 5' read */
  if (up->paire & 0x1)
    for (jr = 0 ; ! ok && products && jr < products->rows ; jr++)
      {
	p1 = ac_table_int (products, jr, 1, 0) ;
	p2 = ac_table_int (products, jr, 2, 0) ;

	frame = (up->a1 + 9999999 - p1 - 1) % 3 ;
	if (!frame && p1 < up->a1 + up->coding && p2 > up->a1)
	  {
	    ok  = TRUE ;
	    up->product = ac_table_obj (products, jr, 0, h) ;
	    up->p1 = p1 ; up->p2 = p2 ; 
	  }
      }
  else /* we just have a 3', the frame is unknown */
    {
      /* look for longest product intersecting the 5' end of the 3' read */
      for (jr = p21 = 0 ; products && jr < products->rows ; jr++)
	{
	  p1 = ac_table_int (products, jr, 1, 0) ;
	  p2 = ac_table_int (products, jr, 2, 0) ;
	  
	  if (p1 < up->a1 + up->coding && p2 > up->a1)
	    if (p2 - p1 > p21) p21 = p2 - p1 ;
	}
      /* if p21 > 0, find again the best */
      for (jr = 0 ; p21 && products && jr < products->rows ; jr++)
	{
	  p1 = ac_table_int (products, jr, 1, 0) ;
	  p2 = ac_table_int (products, jr, 2, 0) ;
	  
	  if (p1 < up->a1 + up->coding && p2 > up->a1 && p2 - p1 == p21)
	    {
	      ok = TRUE ;
	      up->product = ac_table_obj (products, jr, 0, h) ;
	      up->p1 = p1 ; up->p2 = p2 ; 
	    }
	}
    }
  
  /* find the best product */ 
  up->nonbestprotein = 1 ;
  /* may be it is the one we selected */
  if (up->product && ac_has_tag (up->product, "Best_product"))
    {
      up->bestProductKey = ac_obj_key (up->product) ;
      up->nonbestprotein = 0 ; 
      up->bestP1 = up->p1 ; up->bestP2 = up->p2 ;
    }
  /* if not find it */
  for (jr = 0 ; !up->bestProductKey && products && jr < products->rows ; jr++)
    {
      product = ac_table_obj (products, jr, 0, h) ;
      if (ac_has_tag (product, "Best_product"))
	{
	  p1 = ac_table_int (products, jr, 1, 0) ;
	  p2 = ac_table_int (products, jr, 2, 0) ;
	  up->bestP1 = p1 ; up->bestP2 = p2 ; up->bestProductKey = ac_obj_key (product) ;
	}
    }
  if (up->product)
    {
      up->productKey = ac_obj_key (up->product) ;
      up->ecoding = ac_tag_int (up->product, "Coding_length", 0) ;
    }
  
  tbl1 = ac_tag_table (worf, "Anomalous_clone", h) ;
  vtxtClear (anoTxt) ;
  for (jr = 0 ; tbl1 && jr < tbl1->rows ; jr++)
    {
      ccp = ac_table_printable (tbl1, jr, 0, 0) ;
      if (!ccp )
	continue ;
      else if (!strcmp (ccp, "Comment"))
	{
	  if ((ccp = ac_table_printable (tbl1, jr, 1, 0)))
	    {
	      vtxtDot (anoTxt) ; vtxtPrint (anoTxt, ccp) ;
	    }
	}
      else if (!strcasecmp (ccp, "Internal_priming_on_A_rich"))
	continue ;
      else if (!strcasecmp (ccp, "Fuse_to_clone"))
	continue ;
      else
	{
	  vtxtComma (anoTxt) ; vtxtPrint (anoTxt, ccp) ;
	}	    
    }
  ccp = vtxtPtr (anoTxt) ;
  if (ccp)
    dictAdd (nw->dict, ccp, &(up->anomaly)) ;

  tbl1 = ac_tag_table (up->clo, "rearrayed_as", h) ;
  vtxtClear (anoTxt) ;
  for (jr = 0 ; tbl1 && jr < tbl1->rows ; jr++)
    {
      ccp = ac_table_printable (tbl1, jr, 0, 0) ;
      if (!ccp )
	continue ;
     else
	{
	  vtxtComma (anoTxt) ; vtxtPrint (anoTxt, ccp) ;
	}	    
    }
  ccp = vtxtPtr (anoTxt) ;
  if (ccp)
    dictAdd (nw->dict, ccp, &(up->rearray)) ;

  if (1)
    {
      if (up->a1 == up->p1) 
	up->atg = 1 ; /* see the A of the ATG */
      else if (up->a1 == up->p1 + 1) 
	up->atg = 2 ; /* see the T of the ATG */
      else if (up->a1 == up->p1 + 2) 
	up->atg = 3 ; /* see the G of the ATG */
      else if (up->a1 > up->p1 + 3 && !(up->paire & 0x1))
	up->atg = 4 ; /* NA */
    }
  /* check if the ost corresponds to a first ATG downstreeam */
  x = up->product ? 3*ac_tag_int (up->product, "First_ATG", 0)  : 0 ;
  if (x > 0)
    {
      if (up->a1 == up->p1 + x - 3) 
	up->atg2 = 1 ; /* see the A of the ATG */
      else if (up->a1 == up->p1 + x - 3 + 1) 
	up->atg2 = 2 ; /* see the T of the ATG */
      else if (up->a1 == up->p1 + x - 3 + 2) 
	up->atg2 = 3 ; /* see the G of the ATG */
      else if (up->a1 > up->p1 + x - 3 + 3 && !(up->paire & 0x1))
	up->atg2 = 4 ; /* NA */
    }
  
  if (up->product && up->a1 >  up->p1 + 2)
    {
      if (ac_has_tag (up->product, "NH2_Complete"))
	x = up->a1 - up->p1 ;
      else
	x = up->a1 - 1 ;
      if (up->paire & 0x1)
	up->partial5p = x > 0 ? x-1 : 0 ; /* because we should strat in frame 2 */
      else
	up->missing5pRead = x ;
    }

  if (up->a1 < up->p1 && up->product &&
      (x = ac_tag_int (up->product, "Up_stop", 0)) && /* negative value */
      up->a1 < up->p1 + x
      )
    up->contains5pStop = up->p1 - up->a1 ;	  
  if (((up->a1 + 99999 - up->p1 + 1) %3) == 2)
    up->inframe = TRUE ;
  if (!up->product)
    up->coding = nwCodingLength (up->mrna, up->a1 - 1, h) ;
  else
    {
      int u1 = up->a1 > up->p1 ? up->a1 : up->p1 ;
      int u2 = up->a2 < up->p2 ? up->a2 : up->p2 ;
      x = u2 > u1 ? u2 - u1 + 1 : 0 ; x /= 3 ; x = 3 * x ;
      if (ac_has_tag (up->product, "Down_stop") && up->a2 >= up->p2) x-= 3 ;
      up->coding = x ;
    }
  
  if (up->product && ac_has_tag (up->product, "Down_stop"))
    {
      if (up->a2 == up->p2) up->stop = 3 ;     /* contains 3 letters of our end  */
      if (up->a2 == up->p2 - 1) up->stop = 2 ; /* contains 2 letters of our end  */
      if (up->a2 == up->p2 - 2) up->stop = 1 ; /* contains 1 letters of our end  */
      
      if (up->a2 > up->p2) 
	up->contains3pStop = up->a2 - up->p2 ;
    }
  
  if (up->a2 < up->p2 - 2) 
    {
      if (up->paire & 0x2)
	up->partial3p = up->p2 - up->a2 ;
      else if (!up->stop)
	{
	  up->missing3pRead = up->p2 - up->a2 ; 
	  up->stop = 4 ;
	}
    }
} /* nwGetMrnaData */

/*************************************************************************************/

static void nwGetData (NW *nw)
{
  int err, iclo ;
  AC_OBJ mrna = 0 ;
  AC_HANDLE h ;
  AC_TABLE tblClo ;
  AC_OBJ clo = 0 ;
  Array aa = 0 ;
  AC_ITER iter = ac_dbquery_iter (nw->db, messprintf ("Find tg %s ; > mrna", nw->template ? nw->template : ""), 0) ;
  
  printf ("Clone\tCoding potential\tAligned 5' and 3' reads") ;
  printf ("\tInsert length\tCoding length\tDifferences with the genome\tGap between the 2 reads") ;
  printf ("\tMissing 5' sequence\tMissing 3' sequence") ;
  printf ("\tmRNA") ;
  printf ("\tContains 5' stop, 5' UTR length\tContains 3' stop, 3' UTR length") ;
  printf ("\tStarts on our start\tStarts on first ATG") ;
  printf ("\tPartial 5' missing bp") ;
  printf ("\tPartial 3' missing bp") ;

  printf ("\tEnds on validated stop") ;

  printf ("\tconfirmed protein ends") ;

  printf ("\n") ;

  if (iter)
    while (ac_free (mrna) , mrna = ac_iter_obj (iter))
      {
	h = ac_new_handle () ;
	aa = arrayHandleCreate (12, NWZ, h) ;
	tblClo= ac_tag_table (mrna, "cdna_clone", h) ;  
	for (iclo = 0 ; tblClo && iclo < tblClo->rows ; iclo++)
	  {
	    clo = ac_table_obj (tblClo, iclo, 0, h) ;
	    if (! ac_has_tag (clo, "WORF") &&   /* WORF1  and WORF3 */
		strncasecmp (ac_name (clo), "OST_", 4)) /* any other est */
	      continue ;
	    err = ac_table_int (tblClo, iclo, 3, 99999) ;
	    aa = arrayReCreate (aa, 12, NWZ) ;
	    nwGetMrnaData (nw, mrna, clo, err, aa, h) ; 
	    arraySort (aa, cloneOrder) ; 
	    nwReportClo (aa, nw->bb, h) ;
	  }
	ac_free (h) ;
      }

  ac_free (iter) ;
} /* nwGetData */
 
/*************************************************************************************/
/*************************************************************************************/
 
 static void nwExportHtml (NW *nw)
{
  return ;
} /* nwExportHtml */

/*************************************************************************************/

static int byGeneOrder (const void *a, const void *b)
{
  int n ;
  const NWZ *up = (const NWZ *)a, *vp = (const NWZ *)b ;

  n = keySetAlphaOrder (&(up->tgKey), &(vp->tgKey)) ;
  if (n) return n ;
  n = keySetAlphaOrder (&(up->mrnaKey), &(vp->mrnaKey)) ;
  if (n) return n ;
  n = up->score - vp->score ;
  if (n) return n ;
  n = up->err - vp->err ;
  if (n) return n ;
  n = keySetAlphaOrder (&(up->cloKey), &(vp->cloKey)) ;
  if (n) return n ;
  return 0 ;
} /* nwOrder */

/*************************************************************************************/

static void nwExportByGene (NW *nw)
{
  int ii, x ;
  Array bb = nw->bb ;
  NWZ *up ;
  ACEOUT out = aceOutCreateToFile (nw->byGene, "w", nw->h) ;

  if (bb && arrayMax (bb))
    {
      arraySort (bb, byGeneOrder) ;
      aceOutf (out, "Gene") ;
      aceOutf (out, "\tmRNA") ;
      aceOutf (out, "\tProtein AA") ; /* best_product && good_product OR non-coding 
				       * = if both ends confirmed
				       * > otherwise
				 */
      aceOutf (out, "\tOST clone") ; 
      aceOutf (out, "\tRearrayed") ; 
      aceOutf (out, "\tCoding potential") ;
      aceOutf (out, "\tIncomplete sequence") ; /* NA or missing5pRead + gap + missing3pRead */
      aceOutf (out, "\tLength encoded by the OST") ; /* %d STOP if not fully open
						      * NA  if only one read
						      * ~ insert, if gap
						      * = insert,  otherwise
						      */
      aceOutf (out, "\tAnomaly") ; /* clo->anmomalous_clone */
      aceOutf (out, "\n") ;

      for (ii = 0, up = arrp (bb, ii, NWZ) ; ii < arrayMax (bb) ; ii++, up++)
	{
	  aceOutf (out, "%s\t%s", ac_key_name (up->tgKey), ac_key_name (up->mrnaKey)) ;
	  aceOutf (out, "\t%s %d", up->endsConfirmed == 0x3 ? "=" : ">", (up->p2 - up->p1 + 1)/3) ;
	  aceOutf (out, "\t%s", ac_key_name (up->cloKey)) ;
	  aceOutf (out, "\t%s", up->rearray ? dictName (nw->dict, up->rearray) : "") ;
	  aceOutf (out, "\t%s", scoreTitle [up->score]) ;
	  x = up->missing5pRead + up->gap + up->missing3pRead ;
	  if (up->paire == 0x3 || x)
	    aceOutf (out, "\t%d", x) ;
	  else
	    aceOutf (out, "\tNA") ;
	  if (up->contains5pStop) aceOutf (out, "\tSTOP before ATG") ;
	  else if (up->contains3pStop) aceOutf (out, "\t%d STOP", up->p2 - up->a1 + 1) ;
	  else if (up->paire != 0x3) aceOutf (out, "\tNA") ;
	  else if (up->gap)  aceOutf (out, "\t~ %d", up->a2 - up->a1 + 1) ;
	  else aceOutf (out, "\t= %d", up->a2 - up->a1 + 1) ;
	  aceOutf (out, "\t%s", up->anomaly ? dictName (nw->dict, up->anomaly) : "") ;
	  aceOutf (out, "\n") ;
	}
    }
  aceOutDestroy (out) ;
  return ;
} /* nwExportByGene */

/*************************************************************************************/

static void nwExportAce (NW *nw)
{
  int ii, ok ;
  Array bb = nw->bb ;
  NWZ *up ;
  ACEOUT out = aceOutCreateToFile (nw->ace, "w", nw->h) ;
  char buf25[25] ;
  char *pp = 0 ;

  if (bb && arrayMax (bb))
    {
      for (ii = 0, up = arrp (bb, ii, NWZ) ; ii < arrayMax (bb) ; ii++, up++)
	{
	  ac_protect (ac_key_name(up->mrnaKey), nw->h) ;
	  if (up->product)
	    pp = ac_protect (ac_key_name(up->productKey), nw->h) ;
	  else
	    pp = strnew (ac_protect (messprintf("%s.new", ac_key_name(up->mrnaKey)), nw->h), nw->h) ;

	  aceOutf (out, "cDNA_clone %s\n", ac_protect (ac_key_name(up->cloKey), nw->h)) ;
	  aceOutf (out, "-D ORF_quality %s\n", pp) ; /* clean up previous annotations */
	  aceOutf (out, "ORF_quality %s Date %s\n", pp, timeShow (timeNow(), buf25, 25)) ;
	  aceOutf (out, "ORF_quality %s Title %s\n", pp, ac_protect(scoreTitle [up->score], nw->h)) ;
	  if (up->paire == 0x3)
	    aceOutf (out, "ORF_quality %s Both\n", pp) ;
	  else if (up->paire == 0x1)
	    aceOutf (out, "ORF_quality %s Aligned_5p\n", pp) ;
	  else if (up->paire == 0x2)
	    aceOutf (out, "ORF_quality %s Aligned_3p\n", pp) ;
	  if (up->gap && up->paire == 0x3)
	    aceOutf (out, "ORF_quality %s Insert_length  \"~\" %d bp\n", pp, up->insert) ;
	  else if (up->paire != 0x3)
	    aceOutf (out, "ORF_quality %s Insert_length  \">\" %d bp\n", pp, up->insert) ;
	  else
	    aceOutf (out, "ORF_quality %s Insert_length  \"=\" %d bp\n", pp, up->insert) ;
	  if (up->err)
	    aceOutf (out, "ORF_quality %s Mismatch %d\n", pp, up->err ) ;
	  else 
	    aceOutf (out, "ORF_quality %s No_mismatch\n", pp) ;

	  aceOutf (out, "ORF_quality %s Coding_length %d bp\n", pp, up->coding) ;
	  if (up->product)
	    aceOutf (out, "ORF_quality %s Expected_coding_length %d bp\n", pp, up->ecoding) ;
	  if (!up->productKey || up->nonbestprotein)  
	    {         
	      aceOutf (out, "ORF_quality %s Non_best_protein\n", pp) ;
	      if (up->bestProductKey)
		{
		  aceOutf (out, "ORF_quality %s Best_is %s\n", pp
			   , ac_protect (ac_key_name (up->bestProductKey), nw->h)
			   ) ;
		  aceOutf (out, "ORF_quality %s Best_encodes %d bp\n", pp
			   , up->bestP2 - up->bestP1 - 4
			   ) ;
		  if (up->bestP2 < up->a1)
		    aceOutf (out, "ORF_quality %s Best_is_upstream %d bp\n", pp
			     , up->a1 - up->bestP2 
			     ) ;
		  else if (up->bestP1 > up->a2)
		    aceOutf (out, "ORF_quality %s Best_is_downstream %d bp\n", pp
			     , up->bestP1 - up->a2 
			     ) ;
		  else
		    {
		      if (!up->stop && !up->contains3pStop)
			aceOutf (out, "ORF_quality %s Best_is_in_another_frame\n", pp) ;
		      aceOutf (out, "ORF_quality %s Best_intersects_by %d bp\n", pp
			       , (up->bestP2 < up->a2 ? up->bestP2 : up->a2) 
			       - (up->bestP1 > up->a1 ? up->bestP1 : up->a1) + 1
			       ) ;
		    }
		}
	    }
	  if (up->gap)
	    aceOutf (out, "ORF_quality %s Gap %d bp\n", pp,up->gap ) ;
	  if (up->missing5pRead) 
	    aceOutf (out, "ORF_quality %s Missing_5p %d bp\n", pp, up->missing5pRead) ;
	  if (up->missing3pRead) 
	    aceOutf (out, "ORF_quality %s Missing_3p %d bp\n", pp, up->missing3pRead) ;
	  if (up->contains5pStop) 
	    aceOutf (out, "ORF_quality %s Contains_5p_stop %d \"bp 5p UTR length\"\n", pp, up->contains5pStop) ;
	  if (up->contains3pStop) 
	    aceOutf (out, "ORF_quality %s Contains_3p_stop %d bp\n", pp, up->contains3pStop) ;

	  if (up->product && !up->contains5pStop) 
	    {
	      if (up->atg == 2)
		aceOutf (out, "ORF_quality %s start_on_our_start\n", pp) ;
	      else if (up->atg == 1 || up->atg == 3) 
		aceOutf (out, "ORF_quality %s Out_of_frame_Met %d\n", pp, up->atg) ;
	      else if (up->atg2 == 2 && up->partial5p == 0)
		aceOutf (out, "ORF_quality %s start_on_first_Met\n", pp) ;
	      else if (up->atg2 == 2 && up->partial5p < 90)
		aceOutf (out, "ORF_quality %s start_on_first_Met %d \"bp from our start to OST start\"\n", pp, up->partial5p) ;
	      else if ((up->atg2 == 1 || up->atg2 == 3) && up->partial5p < 90)
		aceOutf (out, "ORF_quality %s Out_of_frame_Met %d\n", pp, up->atg2) ;
	      else if ((up->paire & 0x1) && (ok = (up->a1 + 99999 - up->p1 + 1) %3) != 2)
		aceOutf (out, "ORF_quality %s Out_of_frame_partial %d\n", pp, ok ? ok : 3) ;
	      if ((ok = (up->a1 + 99999 - up->p1 + 1) %3) == 2)
		aceOutf (out, "ORF_quality %s In_frame\n", pp) ;
	    }
	  if (up->product && ! up->contains3pStop)
	    {
	      if ((up->stop & 0x3) || ((up->p2 - up->a2 < 3) && (up->a2 - up->p2 < 1)))
		aceOutf (out, "ORF_quality %s Ends_on_our_end\n", pp) ;
	      ok = (up->a2 + 99999 - up->p1 + 1) % 3 ; if (ok == 0) ok = 3 ;
	      if ((up->paire & 0x2) && ok != 2)
		aceOutf (out, "ORF_quality %s Out_of_frame_end %d \"should be 2\"\n", pp, ok) ;
	      else if (((up->paire & 0x2) || (up->stop & 0x3)) && ok == 2)
		aceOutf (out, "ORF_quality %s In_frame_end\n", pp) ;
	    }
	  if (((up->atg | up->atg2) & 0x3) && (up->stop & 0x3))
	    aceOutf (out, "ORF_quality %s Complete\n", pp) ;
	  if (((up->atg | up->atg2) & 0x2) && (up->stop & 0x2) && !up->gap && ! up->err)
	    aceOutf (out, "ORF_quality %s Perfect\n", pp) ;
	  ok = 0 ;
	  if (up->partial5p && ! (up->atg2 &&  up->partial5p < 90))
	    { ok++ ; aceOutf (out, "ORF_quality %s Partial_5p %d bp\n", pp, up->partial5p) ; }
	  if (up->partial3p) 
	    { ok++ ; aceOutf (out, "ORF_quality %s Partial_3p %d bp\n", pp, up->partial3p) ; }
	  if (ok == 2)
	    aceOutf (out, "ORF_quality %s Fragment\n", pp) ;
	  if (up->product && up->pGap)
	    aceOutf (out, "ORF_quality %s P_gap %d bp\n", pp, up->pGap) ;
	  else if (up->endsConfirmed == 0x3) 
	    aceOutf (out, "ORF_quality %s Confirmed_ends\n", pp) ;
	  else if (up->endsConfirmed & 0x1) 
	    aceOutf (out, "ORF_quality %s Confirmed_Met\n", pp) ;
	  else if (up->endsConfirmed & 0x2) 
	    aceOutf (out, "ORF_quality %s Confirmed_stop\n", pp) ;
	  else 
	    aceOutf (out, "ORF_quality %s Unconfirmed_ends\n", pp) ;
	  if (up->anomaly)
	    aceOutf (out, "ORF_quality %s Anomaly %s\n", pp, ac_protect (dictName (nw->dict, up->anomaly), nw->h)) ;
	  if (up->rearray)
	    aceOutf (out, "ORF_quality %s Rearrayed %s\n", pp, ac_protect (dictName (nw->dict, up->rearray), nw->h)) ;
	  aceOutf (out, "\n") ;
	}
    }
  aceOutDestroy (out) ;
  return ;
} /* nwExportAce */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage: nico_worf [-db ACEDB] [-gene template] [-byGene file_name] [-ace file_name] [-html dir] \n") ;
  fprintf (stderr, "// Example:  nico_worf -db ~/yk -gene unc-32 \n") ;  
  fprintf (stderr, "//    export the orfeome clone clone quality and details\n") ;
  fprintf (stderr, "//  -byGene exports a quality summary sorted by gene\n") ;
  fprintf (stderr, "//  -html dir  exports a collection of html files\n") ;
  fprintf (stderr, "//  -ace file_name export a .ace file OST_clone->ORF_quality\n") ;
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  NW nw ;
  AC_HANDLE h = 0 ;
  const char *ici ;
  const char *s = "ok" ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&nw, 0, sizeof (NW)) ;
  nw.h = h ;

  /* optional temple argument */
  getCmdLineOption (&argc, argv, "-gene", &(nw.template)) ;
  getCmdLineOption (&argc, argv, "-html", &(nw.html)) ;
  getCmdLineOption (&argc, argv, "-byGene", &(nw.byGene)) ;
  getCmdLineOption (&argc, argv, "-ace", &(nw.ace)) ;

  /* mandatory database descriptor */
  if (getCmdLineOption (&argc, argv, "-db", &ici))
    {
      nw.db = ac_open_db (messprintf("%s",ici), &s);
      if (!nw.db)
	messcrash ("Failed to open db %s, error %s", ici, s) ;
    }
  else 
     usage () ;

  nw.aa = arrayHandleCreate (20000, NWZ, nw.h) ;
  nw.bb = arrayHandleCreate (20000, NWZ, nw.h) ;
  nw.dict = dictHandleCreate (10000, nw.h) ;
  fprintf (stderr, "// start: %s\n", timeShowNow()) ;
  nwGetData (&nw) ;
  fprintf (stderr, "// analysis done: %s\n", timeShowNow()) ;
  if (nw.byGene)
    nwExportByGene (&nw) ;
  if (nw.ace)
    nwExportAce (&nw) ;
  if (!nw.html)
    nwExportHtml (&nw) ;
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;

  ac_db_close (nw.db) ;
  ac_free (nw.h) ;
  if (0) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  if (0) nwShow (0) ; /* for compiler happiness */

  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

