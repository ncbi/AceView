#include    "../wfiche/biolog.h"
#include    "../wfiche/gtitle.h"
#include "call.h"

#define TABLEBORDER "2"

#define ExonColor0 "<font color=#000000>"
#define ExonColor1 "<font color=#0000ff>"
#define ExonColor2 "<font color=#008888>"

static char *ficheMrnaDNA (vTXT vtxt, AC_DB db, AC_OBJ oMrna, char style, char *param) ;

/******************************************************************************/
/******************************************************************************/
/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  gene ExpressionProfile function
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
	 
static char *ficheGeneExpressioProfile (vTXT vtxt, AC_DB db, AC_OBJ ocDNA_clone, char style)
{
  char orderBy = 'q' ; /* by quality */
  AC_KEYSET clones = ac_objquery_keyset (ocDNA_clone, 0, 0) ;

  vtxtHr (vtxt, 1, 1) ;
  ficheNewCloneParagraph (vtxt, db, clones, style, orderBy, 0) ;
  vtxtHr (vtxt, 1, 1) ;

  ac_free (clones) ;
  return vtxtPtr (vtxt) ;     
}

/******************************************************************************/
/******************************************************************************/
/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  cDNA_clone FICHE function
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
	 
static char *ficheCDNAClone (vTXT vtxt, AC_DB db, AC_OBJ ocDNA_clone, char style)
{
  char orderBy = 'q' ; /* by quality */
  AC_KEYSET clones = ac_objquery_keyset (ocDNA_clone, 0, 0) ;

  vtxtHr (vtxt, 1, 1) ;
  ficheNewCloneParagraph (vtxt, db, clones, style, orderBy, 0) ;
  vtxtHr (vtxt, 1, 1) ;

  ac_free (clones) ;
  return vtxtPtr (vtxt) ;     
}

/******************************************************************************/

static char *ficheCDNACloneGold (vTXT vtxt, AC_DB db, AC_OBJ clone, char style)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ read = ac_tag_obj (clone,  "Read", h) ;
  AC_OBJ tg = read ? ac_tag_obj (read, "From_gene", h) : 0 ;
  GMP *gmp = gmpCreate (db, 0, tg, 0, 0, 0, style, 'o') ;
  vTXT bfr = vtxtHandleCreate (h) ;
  char *cp, *cq ;

  if (gmp) gmp->est = read ;
  vtxtHr (vtxt, 1, 1) ;
  vtxtPrintf (vtxt, "Gold Clone %s", ac_name(clone)) ;
  vtxtBreak (vtxt) ;
  vtxtPrintf (bfr, ac_obj_ace (clone, h)) ;
  cp = vtxtPtr (bfr) ;
  if (vtxt->markUp)
    {
      while (cp && *cp)
	{
	  cq = strstr (cp, "\n") ;
	  *cq = 0 ;
	  vtxtPrint (vtxt, cp) ;
	  vtxtPrint (vtxt, "<br>\n") ;
	  cp = cq+1 ;
	}
    }
  else 
    {
      if (cp)
	vtxtPrint (vtxt, cp) ;
    }

  ficheGoldCloneParagraph (vtxt, gmp) ;
  vtxtBreak (vtxt) ;
  vtxtHr (vtxt, 1, 1) ;

  gmpDestroy (gmp) ;
  ac_free (h) ;

  return vtxtPtr (vtxt) ;     
}

/******************************************************************************/
/******************************************************************************/
/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  RNAi FICHE function
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
	 
static char *ficheRNAi (vTXT vtxt, AC_DB db, AC_OBJ rnai, char style)
{
  vtxtHr (vtxt, 1, 1) ;
  vtxtPrint (vtxt, "Information about the RNAi experiments may be available from the 'Links' tab above this page") ; 
  vtxtHr (vtxt, 1, 1) ;

  return vtxtPtr (vtxt) ;     
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  mRNA CLONES function
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/					
					
static char *ficheProductBlastP (vTXT vtxt, AC_DB db, AC_OBJ oProduct, char style)
{
  vTXT    blkp;
  GMP *gmp = gmpCreate (db, 0, 0, 0, 0, oProduct,  style, 'b') ;
  blkp = vtxtCreate () ;
  if (style == 'x') vtxtMarkup (blkp) ;

  ficheProductBlastPTableParagraph (blkp, gmp) ;

  messfree (gmp) ;
  if (vtxtPtr (blkp))
    vtxtPrint (vtxt,  vtxtPtr (blkp)) ;
  vtxtDestroy (blkp) ;
  return vtxtPtr (vtxt) ;
}

static char *ficheMrnaCdnaClone (vTXT vtxt, AC_DB db, AC_OBJ oMrna, char style)
{
  AC_HANDLE h = ac_new_handle () ;
  GMP *gmp = gmpCreate (db, 0, 0, oMrna, 0, 0, style, 0) ;
  AC_KEYSET clones = ac_objquery_keyset (oMrna, ">cdna_clone", h) ;
  int nn = clones ? ac_keyset_count (clones) : 0 ;

  vtxtHr (vtxt, 1, 1) ;
  gmpSection (vtxt, gmp,"Supporting clones", messprintf ("%d supporting clones for mRNA %s", nn, ac_name (oMrna))) ;

  if (!clones || ! (nn = ac_keyset_count (clones)))
    vtxtPrint (vtxt, "There are no clones supporting this mRNA.") ;
  else
    ficheNewCloneTable (vtxt, gmp, clones, 'q', 0, 0, 0) ;

  vtxtBreak (vtxt) ;
  vtxtHr (vtxt, 1, 1) ;

  gmpDestroy (gmp) ; 
  ac_free (h) ;
  
  return vtxtPtr (vtxt) ;
}

/*
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  _/
  _/  Transcribed_gene CLONES function
  _/
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/					


static char *ficheGeneClone (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT blkp = vtxtHandleCreate (h) ;
  vTXT buf = vtxtHandleCreate (h) ;
  GMP *gmp = gmpCreate (db, oGene, 0, 0, 0, 0, style, 'c') ;
  int nn ;
  AC_KEYSET ksr1, ksr2, ksc1 ;
  vTXT title = vtxtHandleCreate (h) ;
  int numReads = 0, numBuriedReads = 0, numClones = 0 ;

  if (style == 'x'){  vtxtMarkup (blkp) ; vtxtMarkup (buf) ; }
  if (oGene)
    { 
      if (gmp->tg)
	{
	  vtxtPrint (blkp, "Gene ") ;
	  gmpObjLink (blkp, gmp, oGene, 0) ;
	  gmpBlockInit (blkp, gmp, TRUE,10) ;
 	  if (0) vtxtPrint (blkp, " (click for a full description) ") ;
	  gmpChapter (blkp, gmp, "*clone_AceView_Summary", "Summary") ;
	  ficheNewGeneAceViewSummary (blkp, gmp) ;
	  gmpChapterClose (blkp, gmp, "clone_AceView_Summary", TRUE) ; 
	  ficheNewGeneCompactDiagramChapter (blkp, gmp, 4) ;
	}

      ksr1 = ac_objquery_keyset (gmp->gene, ">transcribed_gene; >Read", h) ;
      ksr2 = ksr1 ? ac_ksquery_keyset (ksr1, ">buries;>buried_est;", h) : 0 ;
	  /*
	    if we had access to the buried_reads->cdna_clones we would do:
	    AC_KEYSET ksc1 = ac_objquery_keyset (gmp->gene, ">transcribed_gene; >cdna_clone", h) ;
	    AC_KEYSET ksc2 = ksr2 ? ac_ksquery_keyset (ksr2, ">cdna_clone", h) : 0 ;
	    numClones = ac_keyset_count (ksc1) + (ksc2 ?  ac_keyset_count (ksc2) : 0) ;
	    numReads = ac_keyset_count (ksr1) + (ksr2 ?  ac_keyset_count (ksr2) : 0) ;
	  */
      ksc1 = ac_objquery_keyset (gmp->gene, ">transcribed_gene; >cdna_clone", h) ; 
      numReads = ksr1 ? ac_keyset_count (ksr1) : 0 ;
      numBuriedReads = ksr1 ? ac_keyset_count (ksr2) : 0 ;
      numClones = ksc1 ? ac_keyset_count (ksc1) : 0 ;
      nn = numBuriedReads ? numReads + numBuriedReads - numClones : 0 ;
      vtxtPrintf (title, "%d %s  clones for gene %s"
		  , numClones
		  , nn ? "representative" : "supporting"
		  , ac_name (oGene)
		  ) ;
      gmpChapter (blkp, gmp, "*clone_table", vtxtPtr (title)) ;
      ficheNewCloneTable (blkp, gmp, ksc1, 'q', 0, 0, 0) ;
      if (nn > 0) vtxtPrintf (blkp, "%d other accessions, mainly redundant or of lesser quality, also belong to this gene but are not reported.", nn) ;
      gmpChapterClose (blkp, gmp, "clone_table", TRUE) ; 

      gmpJumpPointShow (blkp, gmp, TRUE, FALSE) ;
    }

  gmpDestroy (gmp) ;
  if (vtxtPtr (blkp))
    vtxtPrint (vtxt,  vtxtPtr (blkp)) ;

  ac_free (h) ;

  return vtxtPtr (vtxt) ;
}

static char *ficheMrnaSupport (vTXT vtxt, AC_DB db, AC_OBJ tg, char style, char *param)
{
  AC_HANDLE h = ac_new_handle () ;

  ac_free (h) ;
  return "ficheMrnaSupport is not yet programmed" ;
}

static char *ficheTranscribed_geneSupport (vTXT vtxt, AC_DB db, AC_OBJ tg, char style, char *param)
{
  AC_HANDLE h = ac_new_handle () ;
  int ir, ir0 = 0, jr, a1, a2, x1 = 0, x2 = 0 ;
  AC_OBJ oGene, est, clo ;
  GMP *gmp = 0 ;
  char *typeNam = "" ;
  AC_KEYSET ks1 ;
  int type = 0 ;
  vTXT title = vtxtHandleCreate (h) ;
  vTXT blkp = vtxtHandleCreate (h) ;
  vTXT buf = vtxtHandleCreate (h) ;
  AC_TABLE tbl ;
  
  if (param && sscanf (param, "-x1 %d -x2 %d", &x1, &x2) >= 2) ;
  else goto done ;

  oGene = ac_tag_obj (tg, "Gene", h) ;
  if (oGene)
    gmp = gmpCreate (db, oGene, 0, 0, 0, 0, style, 'c') ;
  if (!gmp)
    goto done ;
  gmp->tg = tg ;
  
  ks1 = ac_new_keyset (gmp->db, h) ; 
  if (!type) /* try an intron */
    {
      tbl = ac_tag_table (tg, "Intron_boundaries", h) ; 
      for (ir = jr = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  a1 = ac_table_int (tbl, ir, 2, -1) ;
	  a2 = ac_table_int (tbl, ir, 3, -1) ;
	  if (a1 == x1 && a2 == x2)
	    {
	      clo = ac_table_obj (tbl, ir, 4, h) ;
	      ir0 = ir ; jr++ ;
	      ac_keyset_add (ks1, clo) ;
	      type = 1 ; 
	      typeNam = "intron" ;
	    }
	}
    }

  if (!type) /* try an exon */
    {
      tbl = ac_tag_table (tg, "assembled_from", h) ; 
      for (ir = jr = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  a1 = ac_table_int (tbl, ir, 0, -1) ;
	  a2 = ac_table_int (tbl, ir, 1, -1) ;
	  if (a1 == x1 && a2 == x2)
	    {
	      est = ac_table_obj (tbl, ir, 2, h) ;
	      clo = ac_tag_obj (est, "cdna_clone", h) ;
	      ir0 = ir ; jr++ ;
	      ac_keyset_add (ks1, clo) ;
	      type = 2 ;
	      typeNam = "exon" ;
	    }
	}
    }

  if (style == 'x'){  vtxtMarkup (blkp) ; vtxtMarkup (buf) ; }
  

  if (0) vtxtPrint (blkp, " (click for a full description) ") ;
  if (1)
    gmpBlockInit (blkp, gmp, FALSE, 20) ;
  else
    {
      vtxtPrint (blkp, "Gene ") ;
      gmpObjLink (blkp, gmp, oGene, 0) ;
      vtxtPrintf (blkp, ", cDNA suport for the %s %d %d", typeNam, x1, x2) ;
      vtxtBreak (blkp) ;

      gmpBlockInit (blkp, gmp, TRUE, 20) ;
      gmpChapter (blkp, gmp, "*clone_AceView_Summary", "Summary") ;
      ficheNewGeneAceViewSummary (blkp, gmp) ;
      gmpChapterClose (blkp, gmp, "clone_AceView_Summary", TRUE) ;
      ficheNewGeneCompactDiagramChapter (blkp, gmp, 5) ;
    }

  vtxtPrintf (title, "%d accession%s exactly support the %s%s%s from base %d to %d in gene %s"
	      , jr, _multi(jr)
	      , type == 1 ? ac_table_printable (tbl, ir0, 0, "") : "" 
	      , type == 1 && ac_table_printable (tbl, ir0, 0, 0)? " " : ""
	      , typeNam
	      , x1, x2
	      , ac_name (oGene)
	      ) ;

  ficheNewTgSupportChapter (blkp, gmp, ks1, vtxtPtr (title)) ;
  gmpJumpPointShow (blkp, gmp, TRUE, FALSE) ;
  gmpDestroy (gmp) ;
  if (vtxtPtr (blkp))
    vtxtPrint (vtxt,  vtxtPtr (blkp)) ;

 done:
  ac_free (h) ;
  return vtxtPtr (vtxt) ;
}

static char *ficheGeneSummary (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style)
{
  AC_HANDLE h = ac_new_handle () ;
  int ir ;
  char *cp ;
  vTXT blkp = vtxtHandleCreate (h) ;
  
  GMP *gmp = gmpCreate (db, oGene, 0, 0, 0, 0, style, 'c') ;
  AC_TABLE tbl ;

  if (style == 'r') vtxtMarkup (blkp) ;
  if (oGene)
    { 
	  char *species ;

	  switch (gmp->Spc)
	    {
	    case WORM: species = "worm" ; break ;
	    default: species = "human" ; break ;
	    }
      if (style == 'r')
	vtxtPrintf (vtxt, "<h2>Gene <a href=\"https://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/av.cgi?db=%s&q=%s\">%s</a></h2>\n", species, ac_name(oGene), ac_name(oGene)) ;
      else
	vtxtPrintf (vtxt, "Gene %s\n", freeprotect(ac_name(oGene))) ;
      if ((cp = ficheNewGeneAceViewSummary (blkp, gmp)))
	{
	  if (style == 'r')
	    {
	      vtxtBreak (vtxt) ;
	      vtxtPrint (vtxt, "Summary:\n") ;
	      vtxtPrint (vtxt, cp) ; 
	      vtxtBreak (vtxt) ;
	      vtxtBreak (vtxt) ;
	    }
	  else
	    {
	      tbl = ac_tag_table (oGene, "GeneId", h) ;
	      if (tbl && tbl->rows > 0)
		{
		  vtxtPrint (vtxt, "GeneId") ;
		  for (ir = 0 ; ir < tbl->rows ; ir++)
		    vtxtPrintf (vtxt, " %s", ac_table_printable (tbl, ir, 0, "")) ;
		  vtxtPrint (vtxt, "\n") ;
		}
	      vtxtPrint (vtxt, "Summary\n") ;
	      vtxtPrint (vtxt, cp) ; 
	      vtxtBreak (vtxt) ;
	      vtxtPrint (vtxt, "End_Summary ") ;
	    }
	}
       vtxtPrint (vtxt, "\n") ;
    }
  gmpDestroy (gmp) ;
  ac_free (h) ;
  
  return vtxtPtr (vtxt) ;
}

static char *ficheTranscribed_geneClone (vTXT vtxt, AC_DB db, AC_OBJ oTg, char style)
{
  AC_OBJ oGene = ac_tag_obj (oTg, "Gene", 0) ;

  ficheGeneClone (vtxt, db, oTg, style) ;
  ac_free (oGene) ;
  return vtxtPtr (vtxt) ;
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  DNA FICHE function
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/					

static char *ficheGeneDNAZone (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style, char *param)
{
  char *sDna = ac_obj_dna (oGene, 0)  ;

  if (sDna)
    {
      GMP *gmp = gmpCreate (db, 0, 0, 0, 0, 0, style, 'z') ;
      char *cp, *cp0 ;
      int x1, x2 ;

      cp = cp0 = strnew (sDna, 0) ;


      if (sscanf (param, "-p %d %d", &x1, &x2) == 2 &&
	  x1 >= 1 && x2 <= strlen (cp))
	{ *(cp+x2) = 0; cp += x1 - 1 ; }
      else
	{ x1 = 1 ; x2 = strlen (sDna) ; }
      gmpSection (vtxt, gmp, "tg_dna"
		  , messprintf ("Gene %s, DNA Sequence from %d to %d len %d"
				, ac_name (oGene) , x1, x2, x2 - x1 + 1)) ;

      vtxtSequence (vtxt, cp) ;

      ac_free (sDna) ;
      gmpDestroy (gmp) ;
      messfree (cp0) ;
    }
  else
    vtxtPrint (vtxt, "No DNA sequence associated with this object.") ;
  return vtxtPtr (vtxt) ;
}

typedef struct { int ir, len ; } PLEN ;

static int plenOrder(const void *a, const void *b)
{
  const PLEN *u = (const PLEN *)a, *v = (const PLEN *)b ;
  return v->len - u->len ; /* longest ali first */
}
					
static char *ficheGeneDNA (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style, char *param)
{
  AC_OBJ oProd, oRf, oMrna, oModel ;
  AC_TABLE oProds ;
  char *sPep, *sDna ;
  int ir, ix ;
  Array aa = 0 ;
  PLEN *up ;
  GMP *gmp = 0 ;
  
  if (param && *param)
    return ficheGeneDNAZone (vtxt, db, oGene, style, param) ;
  oProds = ac_tag_table (oGene, "Product", 0) ;
  aa = arrayCreate (20, PLEN) ;
  gmp = gmpCreate (db, oGene, 0, 0, 0, 0, style, 'x') ;
  vtxtClear (vtxt) ;
  vtxtPrint (vtxt, 
	     "This page gives, in fasta format, first the sequence of each protein, "
	     "then the sequence of the mRNA as derived from the underlying genome, "
	     "then the sequence of the AM, where each AM sequence is a "
	     "\"golden path\" composite of cDNAs, where we choose, for each segment, "
	     "the clone compatible with the intron structure of the variant "
	     "that best matches the genome.\n<br/>") ;
  vtxtBreak (vtxt) ;
  vtxtPrint (vtxt, 
	     "If you are interested in some particular sequence data, please go"
	     "to the tables of transcripts for the mRNA "
	     "(sequence reconstructed from the genome or from the consensus of the cDNAs" 
	     "guided by the genome AM, premessenger 3'UTR, and the 5 kb upstream of the "
	     "transcript on the genome (probably containing the promotor "
	     "in case of mRNA with complete CDS). Go to the protein table  for the "
	     "sequences of proteins (deduced from the "
	     "genome), and to the introns and exons table for the sequences of "
	     "introns and exons."
	     "The sequences of primers to amplify the CDS"
	     "are given in each mRNA page."
	     ) ;
  vtxtBreak (vtxt) ;
  vtxtPrint (vtxt, 
	     "Please tell us if you would need other sequences."
	     ) ;
  vtxtBreak (vtxt) ;
  for (ir = 0 ; oProds && ir < oProds->rows ; ir++)
    {
      if ((oProd = ac_table_obj (oProds, ir, 0, 0)) &&
	  (sPep = ac_obj_peptide (oProd, 0)))
	{ 
	  up = arrayp (aa, ir, PLEN) ;
	  up->ir = ir ;
	  up->len = strlen (sPep) ;
	  ac_free (sPep) ;
	}
    }
  arraySort (aa, plenOrder) ;

  for (ix = 0 ; oProds && ix < oProds->rows ; ix++)
    {
      ir = arr (aa, ix, PLEN).ir ;
      if ((oProd = ac_table_obj (oProds, ir, 0, gmp->h)) &&
	  (sPep = ac_obj_peptide (oProd, gmp->h)))
	{
	  gmp->product = oProd ;
	  gmp->mrna = ac_tag_obj (oProd, "mRNA", gmp->h) ;
	  if (ac_has_tag (oGene, "Use_AM"))
	    vtxtPrintf (vtxt, ">%s %dAA, protein deduced from the cDNA sequences<br/>\n"
			, ac_name (oProd) , strlen (sPep)) ;
	  else
	    vtxtPrintf (vtxt, ">%s %dAA, protein deduced from the genome sequence<br/>\n"
			, ac_name (oProd) , strlen (sPep)) ;
	  if ((oModel = ac_tag_obj (oProd, "Identical_to", 0)))
	    {
	      vtxtPrintf (vtxt, "// identical to %s<br/>\n", ac_name(oModel)) ;
	      ac_free (oModel) ;
	    }
	  else
	    {
	      BOOL hasStop = ac_has_tag (oProd, "Down_stop") ;
	      
	      gmp->peptide = sPep ;
	      ficheNewMrnaDecoratedPeptide (vtxt, gmp, hasStop) ;
	    }
	}
      ac_free (oProd) ;
    }
   for (ix = 0 ; oProds && ix < oProds->rows ; ix++)
    {
      ir = arr (aa, ix, PLEN).ir ;
      oMrna = 0 ;
      if ((oProd = ac_table_obj (oProds, ir, 0, 0)) &&
	  (oMrna = ac_tag_obj (oProd, "mRNA", 0)) &&
	  (sDna = ac_obj_dna (oMrna, 0)))
	{
	  ficheMrnaDNA (vtxt, gmp->db, oMrna, gmp->style, 0) ;
	  ac_free (sDna) ;
	} 
      ac_free (oMrna) ;
      ac_free (oProd) ;
    }
  for (ix = 0 ; oProds && ix < oProds->rows ; ix++)
    {
      ir = arr (aa, ix, PLEN).ir ;
      oMrna = oRf = 0 ;
      if ((oProd = ac_table_obj (oProds, ir, 0, 0)) &&
	  (oMrna = ac_tag_obj (oProd, "mRNA", 0)) &&
	  (oRf = ac_tag_obj (oMrna,"RefSEqMaker", 0)) &&
	  (sDna = ac_obj_dna (oRf, 0)))
	{
	  vtxtPrintf (vtxt, ">%s.AM sequence %dbp, derived from cDNA clones with best match to genome<br/>", ac_name (oMrna) , strlen (sDna)) ;
	  vtxtSequence (vtxt, sDna) ;
	  ac_free (sDna) ;
	} 
      ac_free (oRf) ;
      ac_free (oMrna) ;
      ac_free (oProd) ;
    }

  arrayDestroy (aa) ;

  if (gmp->markup)
    vtxtPrint (vtxt, "<br/><br/><i>You may copy paste this page or save the frame source and remove the html markups.</i><br/>\n") ;
  
  gmpDestroy (gmp) ;

  return vtxtPtr (vtxt) ;
}

static char *ficheDNA (vTXT blkp, AC_DB db, AC_OBJ oMrna, char style, char *param)
{
  char *sDna, remote[1024] ;
  int x1 = 0, x2, color, nn, n1 = 0, nLIMIT = 40 ;
  static DICT *dict = 0 ;
  static Array count = 0 ;

  if (!dict)
    {
      dict = dictCreate (500) ;
      count = arrayCreate (500, int) ;
    }

  *remote = 0 ;
  if (param && sscanf (param, "-remote %s -x1 %d -x2 %d -color %d", remote, &x1, &x2, &color) >= 2)
    sDna = ac_zone_dna (oMrna, x1, x2, 0) ;
  else if (param && sscanf (param, "-x1 %d -x2 %d -color %d", &x1, &x2, &color) >= 2)
    sDna = ac_zone_dna (oMrna, x1, x2, 0) ;
  else  if (param && sscanf (param, "-remote %s", remote) >= 1)
    sDna = ac_obj_dna (oMrna, 0) ;
  else
    sDna = ac_obj_dna (oMrna, 0) ;
  
  if (*remote && *(ac_name(oMrna)) == 'y')
    {
      dictAdd (dict, remote, &nn) ;
      n1 = array (count, nn, int) + 1 ;
      array (count, nn, int) = n1 ;
    }
    
  color = color % 3 ;
  if (sDna && n1 < nLIMIT)
    {
      char *title ;
      GMP *gmp = gmpCreate (db, 0, 0, 0, 0, 0, style, 'z') ;
      gmpBlockInit (blkp, gmp, FALSE, 50) ;

      if (!x1)
	title = messprintf (">DNA Sequence %s %d bp", ac_name (oMrna) , strlen (sDna)) ;
      else
	title = messprintf (">DNA fragment %d bp", strlen (sDna)) ;

      vtxtBreak (blkp) ;
      if (style == 'x')
	{
	  vtxtPrintf (blkp, "<br/>\n<font color=#007f7f><b>%s </b></font>", title) ;
	  /* we cannot link back to a sequence class object */
	  vtxtPrintf (blkp, "<br/>\n") ;
	}
      else
	vtxtPrintf (blkp, "%s\n", title) ;
    
       if (style == 'x')
	switch (color)
	  {
	  case 0: vtxtPrint (blkp, ExonColor0) ;  break ;
	  case 1: vtxtPrint (blkp, ExonColor1) ;  break ;
	  case 2: vtxtPrint (blkp, ExonColor2) ;  break ;
	  }
      vtxtSequence (blkp, sDna) ;
      if (style == 'x')
	vtxtPrint (blkp, "</font>") ;
      ac_free (sDna) ;


      gmp->view = 'z' ; /* to link the menus */
      gmpJumpPointShow (blkp, gmp, TRUE, FALSE) ;
      gmpDestroy (gmp) ;
    }
  else if (sDna && n1 >= nLIMIT)
    {
      vtxtPrint (blkp
		 , "<br/>Sorry, this looks like a download of private data, please try again tomorrow or contact us<br/>") ;
      callSystem (messprintf ("echo 'yk dna download n=%d from %s' | Mail mieg", n1, remote)) ;
    }
  else
    vtxtPrint (blkp, "No DNA sequence associated with this object.") ;

  return vtxtPtr (blkp) ;
}

char *ficheMrnaDNA (vTXT blkp, AC_DB db, AC_OBJ oMrna, char style, char *param)
{
  char *sDna, *cp ;
  AC_HANDLE h = handleCreate () ;
  int x1 = 0, x2 ;
  extern   void ficheNewMrnaDnaDecoratedSequence (vTXT blkp, GMP *gmp, int x1, int x2) ;

  if (param && (cp = strstr (param, "-x1")) && sscanf (cp, "-x1 %d -x2 %d", &x1, &x2) == 2) ;
  else
    x1 = x2 =0 ;
 
  if (x1 == 99999 && x2 == 99999) /* hack, the web pages wants the colored genome region */
    {
      ac_free (h) ;
      return ficheMrnaGenomeDNA (blkp, db, oMrna, style, FALSE) ;
    }
  sDna = ac_obj_dna (oMrna, h) ;
  if (sDna)
    {
      GMP *gmp = gmpCreate (db, 0, 0, oMrna, 0, 0, style, 'm') ;
      int ln = strlen (sDna) ;
      char *title ;

      gmpBlockInit (blkp, gmp, FALSE, 50) ;

      if (!x1)
	title = messprintf (">%s mRNA Sequence %d bp, derived from the genome", ac_name (oMrna) , ln) ;
      else
	title = messprintf (">%s.mRNA_fragment from %d to %d, length %d bp", ac_name (oMrna), x1, x2, x2 - x1 + 1) ;
      
       if (style == 'x')
	 {
	   vtxtPrintf (blkp, "<br/>\n<font color=#007f7f><b>%s </b></font>", title) ;
	   if (0)
	     { /* BUG MOI est mal defini */
	      gmpObjLink (blkp, 0 /* should be gmp */, oMrna
			 ,"<img SRC='images/arrowup.gif' border='0' width='14' height='14'></a>") ;
	     }
	  /*	    gmpHelp (blkp, gmp, jmp->file, jmp->help) ; */
	  vtxtPrint (blkp, "<br/>\n") ;
	}
       else
	 vtxtPrintf (blkp, "%s\n", title) ;
       ficheNewMrnaDnaDecoratedSequence (blkp, gmp, x1, x2) ;
       gmp->view = 'g' ; /* to link the menus */
       gmpJumpPointShow (blkp, gmp, TRUE, FALSE) ;
       gmpDestroy (gmp) ;
    }
  else
    vtxtPrint (blkp, "Sorry, there is no DNA sequence associated with this object.") ;
  ac_free (h) ;
  return vtxtPtr (blkp) ;
}

static char *fichePEP (vTXT vtxt, AC_DB db, AC_OBJ oProduct, char style, char *param)
{
  char *sPep, *cp ;
  int x1 = 0, x2 ;

  sPep = ac_obj_peptide (oProduct, 0) ;
  if (sPep &&
      param && (cp = strstr (param, "-x1")) && sscanf (cp, "-x1 %d -x2 %d", &x1, &x2) == 2 &&
      x1 > 0 && x2 >= x1 && x2 < strlen (sPep))
    { sPep [x2] = 0 ; sPep += x1 - 1 ; }
  
  if (sPep)
    {
      BOOL hasStop = ac_has_tag (oProduct, "Down_stop") ;
      GMP *gmp = gmpCreate (db, 0, 0, 0, 0, oProduct, style, 'x') ;
      gmp->peptide = sPep ;
      gmpBlockInit (vtxt, gmp, FALSE, 50) ;

      if (!x1)
	gmpSection (vtxt, gmp, "tg_pep"
		    , messprintf (">Peptide Sequence %s %d AA", ac_name (oProduct) , strlen (sPep))) ;
      else
	gmpSection (vtxt, gmp, "tg_pep"
		    , messprintf (">Peptide fragment %d AA", strlen (sPep))) ;

      if (0) vtxtSequence (vtxt, sPep) ;
      else ficheNewMrnaDecoratedPeptide (vtxt, gmp, hasStop) ;

      gmp->view = 'g' ; /* to link the menus */
      gmpJumpPointShow (vtxt, gmp, TRUE, FALSE) ;
      ac_free (sPep) ;
      gmpDestroy (gmp) ;
    }
  else
    vtxtPrint (vtxt, "No DNA sequence associated with this object.") ;
  return vtxtPtr (vtxt) ;
}


/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  MRNA FICHE function
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static char* ficheMrnaContent (vTXT blkp, GMP *gmp, int section)
{ 
  AC_HANDLE h = ac_new_handle () ;

  if (gmp->markup) 
    {
      vtxtPrint (blkp, "<a href='#' name='top'></a>\n") ;
    }
  gmpBlockInit (blkp, gmp, TRUE, 30) ;

  ficheNewMrnaSummaryChapter (blkp, gmp) ;
  ficheNewGeneGenomeDiagramChapter (blkp, gmp) ;
  ficheNewMrnaFlashDiagramChapter (blkp, gmp) ;
  ficheNewMrnaExpressionCloneSupportChapter  (blkp, gmp) ;
  ficheNewMrnaStructureChapter (blkp, gmp) ;
  ficheNewMrnaProteinChapter (blkp, gmp) ;
  if (gmp->nMrna > 1 && ac_has_tag (gmp->gene, "Structure") && ! ac_has_tag (gmp->gene, "Cloud_gene"))
    {
      vtxtEmptyLine (blkp, 0) ;
      ficheListAndMarkUpMrnas (blkp, gmp, 'S', FALSE) ;
    }
  
  gmpJumpPointShow (blkp, gmp, TRUE, TRUE) ;
  
  ac_free (h) ;
  return vtxtPtr (blkp) ;
}  /* ficheMrnaContent */

char *ficheMrnaDo (vTXT blkp, GMP *gmp, int section)
{
  if (gmp->markup) vtxtMarkup (blkp) ;
  switch (gmp->style)
    {
    case 'r':
    case 's':	   
      fAsnGenerateMRNA (blkp, gmp, 0, 1) ;
    case 'x':
    default:
      if (! section)
	ficheMrnaTitleParagraph (blkp, gmp) ;
      ficheMrnaContent (blkp, gmp, section) ;
    }
  return vtxtPtr (blkp) ;
}  /* ficheMrnaDo */

/******************************************************************************/

char *ficheMrnaFiche (vTXT vtxt, AC_DB db, AC_OBJ oMrna, char style, int section)
{
  GMP *gmp ;

  gmp = gmpCreate (db, 0, 0, oMrna, 0, 0, style, 'm') ;
  if (!gmp)
    {
      vtxtPrint (vtxt, "\nSorry, I do not find the details on this gene.\n\n"
		  "It is possible that it was filtered at the last step because its quality was not good enough."
		  "\n") ;
    }
  else
    ficheMrnaDo (vtxt, gmp, section) ;
  gmpDestroy (gmp) ;
  return vtxtPtr (vtxt) ;
}  /* ficheMrnaFiche */

/******************************************************************************/

static char *ficheMrnaSection (vTXT vtxt, AC_DB db, AC_OBJ oMrna, char style, char *param)
{ 
  int s = 0 ;
  if (sscanf (param, "-s %d", &s) == 1 &&  s >= 1)
    ficheMrnaFiche (vtxt, db, oMrna, style, s) ;
  else
    vtxtPrintf (vtxt, "Error in ficheMrnaSection") ;
  return vtxtPtr (vtxt) ;
} /* ficheGeneSection */

/******************************************************************************/

char *ficheMrna (vTXT vtxt, AC_DB db, AC_OBJ oMrna, char style)
{
  return ficheMrnaFiche (vtxt, db, oMrna, style, 0) ;
}  /* ficheMrna */

/******************************************************************************/
/******************************************************************************/

const char *ficheProduct (vTXT vtxt, AC_DB db, AC_OBJ oProduct, char style)
{
  GMP *gmp = gmpCreate (db, 0, 0, 0, 0, oProduct, style, 'm') ;
  if (!gmp)
    {
      vtxtPrint (vtxt, "\nSorry, I do not find the details on this product.\n\n"
		  "It is possible that it was filtered at the last step because its quality was not good enough."
		  "\n") ;
    }
  else
    ficheMrnaDo (vtxt, gmp, 0) ;
  gmpDestroy (gmp) ;
  return vtxtPtr (vtxt) ;
}  /* ficheProduct */

/******************************************************************************/
/******************************************************************************/

static void ficheGeneDo (vTXT blkp, GMP *gmp, int chapter)
{
  AC_HANDLE h = ac_new_handle () ;

  if (gmp->markup) 
    vtxtMarkup (blkp) ;  
   
  if (chapter == 100000)
    {
      vtxtPrintf (blkp, "%6d.<br>", 100001) ; /* essential in main.js file downloader */
      gmp->superDiv = 2 ;   /* so that the chapters are loaded as blocked */
      
      gmpBlockStart (blkp, gmp, "Gene_MOLECULAR") ;
      ficheNewGeneMOLECULESChapter (blkp, gmp) ;
      gmpBlockClose (blkp, gmp, "Gene_MOLECULAR") ;

      gmpBlockStart (blkp, gmp, "Gene_EXPRESSION_REGULATION") ;
      ficheNewGeneExpressionChapter (blkp, gmp) ;
      gmpBlockClose (blkp, gmp, "Gene_EXPRESSION_REGULATION") ;

      gmpBlockStart (blkp, gmp, "Gene_FUNCTION") ;
      ficheNewGenePhenotypeFunctionChapter (blkp, gmp) ;
      ficheNewGeneBiblioChapter (blkp, gmp, 1) ;
      gmpBlockClose (blkp, gmp, "Gene_FUNCTION") ;

      gmpJumpPointShow (blkp, gmp, FALSE, FALSE) ;
    } 
  else
    {
      /* title info */  
      if (gmp->markup) vtxtPrint (blkp, "<a href='#' name='top'></a>\n") ;
      
      ficheNewGeneTitleParagraph (blkp, gmp) ;   
      
      if (chapter == 0) chapter = 1 ;
      gmpBlockInit (blkp, gmp, TRUE, chapter) ;
      gmp->requestedChapter = chapter ;
      switch (chapter)
	{
	case 0:
	  gmpBlockStart (blkp, gmp, "SUMMARY") ;
	  /* fall thru */
	case 1: /* summary */
	  ficheNewGeneSummaryChapter (blkp, gmp) ;
	  if (1)
	    {
	      if (gmp->Spc == HUMAN)
		ficheNewGeneExpressionProfileChapter (blkp, gmp) ;
	      ficheNewGeneGenomeDiagramChapter (blkp, gmp) ;
	      if (0) ficheNewGeneWiggleDiagramChapter (blkp, gmp, 1) ;
	      ficheNewGeneCompactDiagramChapter (blkp, gmp, 1) ;

		  ficheNewGeneSequenceChapter (blkp, gmp) ;
	   
	      if (gmp->Spc == WORM)
		ficheNewGeneLocatorDiagramChapter (blkp, gmp, TRUE) ;
	      else
		ficheNewGeneLocatorDiagramChapter (blkp, gmp, FALSE) ;
	      ficheNewGeneMrnaDiagramChapter (blkp, gmp) ;
	    }
	  ficheNewGeneBiblioChapter (blkp, gmp, 0) ;
	  break ;
	case 2: /* molecules */
	  if (gmp->nMrna > 1) 
	    {
	      ficheNewGeneAltFeatureChapter (blkp, gmp) ;
  	      ficheNewGeneGenomeDiagramChapter (blkp, gmp) ;
	      ficheNewGeneCompactDiagramChapter (blkp, gmp, 2) ;
	      ficheNewGeneMOLECULESChapter (blkp, gmp) ;
	    }
	  break ;
	case 3: /* expression */ 
	  ficheNewGeneCompactDiagramChapter (blkp, gmp, 3) ;
	  ficheNewGeneExpressionChapter (blkp, gmp) ;
	  break ;
	case 4: /* function */
	  ficheNewGenePhenotypeFunctionChapter (blkp, gmp) ;
	  /* never introduce a second chapter, becuase of the call in showJP 
	   *if need be add it inside ficheNewGenePhenotypeFunctionChapter
	   */
	  break ;
	}
      
      if (!chapter)
	{
	  gmpBlockClose (blkp, gmp, "SUMMARY") ;
	  gmpImportRemainder (blkp, gmp, gmp->gene) ;
	}
      else
	gmpJumpPointShow (blkp, gmp, TRUE, TRUE) ;
    }
  ac_free (h) ;
} /* ficheGeneDo */

/******************************************************************************/

static char *ficheGeneFiche (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style, int chapter)
{
  GMP *gmp ;
  
  /* initialise */
  gmp = gmpCreate (db, oGene, 0, 0, 0, 0, style, 'g') ;
 
  if (!gmp)
    {
      vtxtPrint (vtxt, 
		  "\nSorry, I do not find the details on this gene.\n\n"
		  " It is possible that it was filtered at the last step"
		  "because its quality was not good enough.\n") ;
    }
  else
    ficheGeneDo (vtxt, gmp, chapter) ;

  gmpDestroy (gmp) ;
  return vtxtPtr (vtxt) ;
} /* ficheGeneFiche */				

/******************************************************************************/

static char *ficheGeneSection (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style, char *param)
{
  int s = 0 ;
  if (sscanf (param, "-s %d", &s) == 1 &&  s >= 1)
    ficheGeneFiche (vtxt, db, oGene, style, s) ;
  else
    vtxtPrintf (vtxt, "Error in ficheGeneSection") ;
  return vtxtPtr (vtxt) ;
} /* ficheGeneSection */

/******************************************************************************/
/******************************************************************************/

static void ficheGenomeDo (vTXT blkp, AC_DB db, GMP *gmp, int chapter)
{
  AC_HANDLE h = ac_new_handle () ;

  if (gmp->markup) 
    vtxtMarkup (blkp) ;  
   
  if (chapter == 100000)
    {
      vtxtPrintf (blkp, "%6d.<br>", 100001) ; /* essential in main.js file downloader */
      gmp->superDiv = 2 ;   /* so that the chapters are loaded as blocked */
      
      gmpBlockStart (blkp, gmp, "Genome_summary") ;
      ficheGenomeSummaryChapter (blkp, gmp) ;
      gmpBlockClose (blkp, gmp, "Genome_summary") ;

      gmpBlockStart (blkp, gmp, "Genome_plot") ;
      ficheGenomePlotChapter (blkp, db, gmp, TRUE) ;
      ficheGenomePlotChapter (blkp, db, gmp, FALSE) ;
      gmpBlockClose (blkp, gmp, "Genome_plot") ;

      gmpJumpPointShow (blkp, gmp, FALSE, FALSE) ;
    } 
  else
    {
      /* title info */  
      if (gmp->markup) vtxtPrint (blkp, "<a href='#' name='top'></a>\n") ;
      
      if (chapter == 0) chapter = 1 ;

      gmp->requestedChapter = chapter ;
      switch (chapter)
	{
	case 1: /* summary */
	  ficheNewGeneTitleParagraph (blkp, gmp) ;   
	  gmpBlockInit (blkp, gmp, TRUE, chapter) ;
	  ficheGenomeSummaryChapter (blkp, gmp) ;
	  break ;
	case 2: /* plot */
	  ficheGenomePlotChapter (blkp, db, gmp, FALSE) ;	  
	  ficheGenomePlotChapter (blkp, db, gmp, TRUE) ;	  
	  break ;
	}
    }
  ac_free (h) ;
} /* ficheGenomeDo */

/******************************************************************************/

static char *ficheGenomeFiche (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style, int chapter)
{
  GMP *gmp ;
  
  /* initialise */
  gmp = gmpCreate (db, oGene, 0, 0, 0, 0, style, 'g') ;
 
  if (!gmp)
    {
      vtxtPrint (vtxt, 
		  "\nSorry, I do not find the details on this gene.\n\n"
		  " It is possible that it was filtered at the last step"
		  "because its quality was not good enough.\n") ;
    }
  else
    {
      ficheGenomeDo (vtxt, db, gmp, 1) ;
      ficheGenomeDo (vtxt, db, gmp, 2) ;
    }
  gmpDestroy (gmp) ;
  return vtxtPtr (vtxt) ;
} /* ficheGenomeFiche */				

/****************/

static char *ficheGenome (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style)
{
  return ficheGenomeFiche (vtxt, db, oGene, style, 0) ;
} /* ficheGene */

/******************************************************************************/

static char *ficheGene (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style)
{
  return ficheGeneFiche (vtxt, db, oGene, style, 0) ;
} /* ficheGene */

/******************************************************************************/

char *ficheGeneGene (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style)
{
  return ficheGeneFiche (vtxt, db, oGene, style, 1) ;
} /* ficheGeneGene */

/******************************************************************************/

char *ficheGeneMolecules (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style)
{
  AC_HANDLE h = ac_new_handle () ;
  int n = ac_keyset_count (ac_objquery_keyset (oGene, ">Transcribed_gene;>mRNA",h)) ;
  char *cp ;

  if (n > 1) 
    cp = ficheGeneFiche (vtxt, db, oGene, style, 2) ;
  else
    {
      /*       AC_OBJ oMrna = ac_objquery_obj (oGene, "{>transcribed_gene; >MRNA}SETELSE{>genefinder;>predicted_mrna}", h) ; */
      cp = ficheMrna (vtxt, db, oGene, style) ;
    }
  ac_free (h) ;
  return cp ;
} /* ficheGeneMolecules */

/******************************************************************************/

char *ficheGeneExpression (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style)
{
  return ficheGeneFiche (vtxt, db, oGene, style, 3) ;
} /* ficheGeneExpression */

/******************************************************************************/

char *ficheGeneFunction (vTXT vtxt, AC_DB db, AC_OBJ oGene, char style)
{
  return ficheGeneFiche (vtxt, db, oGene, style, 4) ;
} /* ficheGeneFunction */

/******************************************************************************/
/******************************************************************************/

static char *ficheTranscribed_gene (vTXT vtxt, AC_DB db, AC_OBJ oTranscribed_gene, char style)
{
  char *cp ;
  AC_OBJ oGene = ac_tag_obj (oTranscribed_gene, "Gene", 0) ;

  if (style == 's')
    {
      cp = fAsnGenerateGene (db, oGene, 0, 's') ;
      if (cp)
	vtxtPrint (vtxt,  cp) ;
      messfree (cp) ;
    }
  else
    ficheGene (vtxt, db, oGene, style) ;
   
  ac_free (oGene) ; /* may be */
  return vtxtPtr (vtxt) ;
} /* ficheTranscribed_gene */

/******************************************************************************/
/*********************** MULTIFICHES Entry point ******************************/


typedef char* (*ficheFUNCTION) (vTXT vtxt, AC_DB db, AC_OBJ obj, char style, char *param) ;
struct
{
  char *classname, *viewname ;
  char style ;
  ficheFUNCTION   func ;
}ficheFunctions[] =
{

  {"Transcribed_gene", "FICHE", 'x', (ficheFUNCTION) ficheTranscribed_gene}, /* done */
  {"Transcribed_gene", "SUPPORT", 'x', (ficheFUNCTION) ficheTranscribed_geneSupport}, /* done */

  {"Gene", "FICHE", 'x', (ficheFUNCTION) ficheGene},  /* done */
  
  {"Gene", "FGENE", 'x', (ficheFUNCTION) ficheGeneGene},  /* done */
  {"Gene", "FMOL", 'x', (ficheFUNCTION) ficheGeneMolecules},  /* done */
  {"Gene", "FEXP", 'x', (ficheFUNCTION) ficheGeneExpression},  /* done */
  {"Gene", "FFUNC", 'x', (ficheFUNCTION) ficheGeneFunction},  /* done */
  {"Gene", "GENOME", 'x', (ficheFUNCTION) ficheGenome},  /* done */

  {"Gene", "SECTION", 'x', (ficheFUNCTION) ficheGeneSection},  /* done */

  {"Gene", "Summary", 's', (ficheFUNCTION) ficheGeneSummary},  /* done */
  {"Gene", "Summary", 'r', (ficheFUNCTION) ficheGeneSummary},  /* done */

  /*  {"Gene", "Probe", 's', (ficheFUNCTION) ficheGeneProbe},   done */

  {"Transcribed_gene", "CLONES", 'x', (ficheFUNCTION) ficheTranscribed_geneClone},  /* done */

  {"Gene", "CLONES", 'x', (ficheFUNCTION) ficheGeneClone},  /* done */

  {"MRNA", "FICHE", 'x', (ficheFUNCTION) ficheMrna},   /* xml */ /* done */
  {"MRNA", "SUPPORT", 'x', (ficheFUNCTION) ficheMrnaSupport}, /* done */

  {"MRNA", "SECTION", 'x', (ficheFUNCTION) ficheMrnaSection},   /* xml */ /* done */

  {"MRNA", "FICHE", 's', (ficheFUNCTION) ficheMrna},   /* asngb */ /* done */

  {"MRNA", "FICHE", 'r', (ficheFUNCTION) ficheMrna},   /* asnref */ /* done */

  {"MRNA", "CLONES", 'x', (ficheFUNCTION) ficheMrnaCdnaClone},   /* xml*/

  {"Product", "FICHE", 'x', (ficheFUNCTION) ficheProduct},   /* xml*//* done */

  {"Product", "BLASTP", 'x', (ficheFUNCTION) ficheProductBlastP},   /* xml*//* done */

  {"Sequence", "DNA", 'x', (ficheFUNCTION) ficheDNA},   /* xml */ /* done */

  {"MRNA", "DNA", 'x', (ficheFUNCTION) ficheMrnaDNA},   /* xml */ /* done */

  {"Product", "PEP", 'x', (ficheFUNCTION) fichePEP},   /* xml */ /* done */

  {"Gene", "fasta", 'x', (ficheFUNCTION) ficheGeneDNA},   /* xml */ /* done */

  {"Gene", "ExpProfile", 'x', (ficheFUNCTION) ficheGeneExpressioProfile},   /* xml */ /* done */

  {"cDNA_clone", "CLONES", 'x', (ficheFUNCTION) ficheCDNAClone},   /* xml */ /* done */

  {"cDNA_clone", "Gold", 'x', (ficheFUNCTION) ficheCDNACloneGold},   /* xml */ /* done */

  {"Map", "FICHE", 'r', (ficheFUNCTION) ficheChromosomeDump},  /* done */

  {"cDNA_clone", "FICHE", 'x', (ficheFUNCTION) ficheCDNAClone},   /* xml */ /* done */

  {"RNAi", "FICHE", 'x', (ficheFUNCTION) ficheRNAi},   /* xml */ /* done */

  {0, 0, 0, 0}
};

char *swormView (vTXT blkp, char *dbName, char *clsName, char *objName, char style, char *view, char *params)
{
  AC_DB db;
  AC_ITER lTmp = 0 ;
  AC_OBJ obj = 0, clo ;
  int             i;
  const char *err ;
  char    *      resultBufr ;
  static vTXT          myBlkp = 0;

  vtxtDestroy (myBlkp) ;
  if (!blkp)
    { myBlkp = vtxtCreate () ; blkp = myBlkp;}

  if (! (db = ac_open_db (dbName, &err)))
    {
      vtxtPrintf (blkp, "No database connection %s: %s \n", dbName, err) ;
      return vtxtPtr (blkp) ;
    }
 if (! (obj = ac_get_obj (db, clsName, objName, 0)) )
    {
      /*    if (! (lObj = acFetch (db, 0, acFETCHALLTAGS, "query find %s %s ", clsName, objName)) )
	    {*/
      vtxtPrintf (blkp, "%s %s was not found\n", clsName, objName) ;
      ac_db_close (db)  ;
      return vtxtPtr (blkp) ;
    }

  if ((lTmp = ac_query_iter (db, 0, "find Clone Main_clone", 0, 0)) &&
      (clo = ac_next_obj (lTmp)))
    {
      strcpy (genomeRelease, ac_tag_text (clo , "Genome_release", "WS")) ;
      ac_free (clo) ;
    }
  ac_free (lTmp) ;

  /* find fiche function */
  for (i = 0 ; ficheFunctions[i].func ; i++)
    {
      if (strcasecmp (clsName, ficheFunctions[i].classname)) continue;
      if (view && ficheFunctions[i].viewname && strcasecmp (view, ficheFunctions[i].viewname)) continue;
      if (style && ficheFunctions[i].style && (style != ficheFunctions[i].style)) continue;
      break;
    }

  if (ficheFunctions[i].func)
    {/* form the fiches */
      if (obj)
	{
	  vTXT vtxt = vtxtCreate () ;
	  if (style == 'x')
	    {
	      vtxtMarkup (vtxt) ;
	    }
	  gmpJumpPointInit () ;

	  if ((resultBufr = ficheFunctions[i].func (vtxt, db, obj, style, params)) )
	    {
	      { /* it happens vtxtPrintf has a buffer limitation  at sometimes 32M */
		int nn = strlen (resultBufr) , mx = 0xffff ;
		char cc, *text = resultBufr ;
		while (nn > mx)
		  {
		    cc = text[mx] ;
		    text[mx] = 0 ;
		    vtxtPrint (blkp, text)  ;
		    text += mx ;
		    *text = cc ;
		    nn -= mx ;
		  }
		
		vtxtPrint (blkp, text)  ;
		vtxtPrint (blkp, "\n") ;
		/* messfree (resultBufr) ; */
	      }	 
	    }
	  vtxtDestroy (vtxt) ;
	}
    }

  gmpJumpPointDestroy () ;
  
  ac_free (obj) ;
  ac_db_close (db)  ;

  if (!ficheFunctions[i].func)
  
    {
      vtxtPrintf (blkp, "No Fiche is defined for %s %s\n", clsName, objName) ;
      return vtxtPtr (blkp) ;
    }
  return vtxtPtr (blkp) ;
}

/*****************************************************************************/

char *swormGetGeneTitle (vTXT blkp, char *geneName)
{
   AC_DB db =  ac_open_db ("local", 0) ;
   AC_OBJ oGene ;

   if (db)     
     {
       oGene = ac_get_obj (db, "Gene", geneName, 0) ;
       
       if (oGene)
	 {
	   GMP *gmp = gmpCreate (db, oGene, 0, 0, 0, 0, 's', 'g') ;
	   
	   gtGeneTitle (blkp, gmp, FALSE) ;
	   gmpDestroy (gmp) ;
	   ac_free (oGene) ;
	 }       
       ac_db_close (db) ;
     }
   return vtxtPtr (blkp) ;
} /* swormGetGeneTitle */

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  Main functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/
