#include    "../wfiche/biolog.h"
#include    "../wfiche/gtitle.h"
#include    "../wfiche/pfambeautify.h"

#include "a.h"
#include "../whooks/sysclass.h"

int Spc ;
static char *kantorTitleDB = "local";

/**********************************************************************/
/****************** TITLE SUGGESTION **********************************/
/* a few pfam titles are correct but quite horrible */
static char *kantorCleanPfamTitle (const char *title)
{
  PFB *pfb ;
  char *cp, *title2 ;

  for (pfb = pfamBeautify ; pfb->old; pfb++)
    if (!strcasecmp (title, pfb->old))
      return pfb->new ;
  
  title2 = gtCleanUp (title) ;
  cp =  title2 + strlen (title2) - 9 ;
  if (cp < title2) cp =  title2 ;
  if (strstr (cp, "domain"))
    return messprintf ("%s %s", title2, "containing protein") ;

  return title2 ;
} /* kantorCleanPfamTitle */

/*******************/

static BOOL kantorGetBlastPfamTitle (vTXT blkp, AC_OBJ oKantor, char *kantorName)
{
  int ir;
  int isTitle =  0 ;
  AC_TABLE  gOb ;
  AC_OBJ obj ;
  vTXT pfamHintBlk = vtxtCreate () ;
  vTXT blastHintBlk = vtxtCreate () ;

  vtxtClear (blkp) ;

  gOb = ac_tag_table (oKantor, "BlastP", 0) ;
  if (gOb)  /* blast */
    {
      for (ir = 0; ir < gOb->rows;ir++)
	{
	  if (!strncmp (ac_table_printable (gOb, ir, 0, ""), "NP_",3) && /* avoid nematode NP, they are ours */
	      strstr (ac_table_printable (gOb, ir, 7, ""),"aenorhabditis"))
	    continue ;
	  /* problem: briggsae is badlya notated and has too high scores against elegans */
	  if (strstr (ac_table_printable (gOb, ir, 7, ""),"aenorhabditis briggsae"))
	    continue ;
	  vtxtPrintf (blastHintBlk, "%f %s\n", ac_table_float (gOb, ir, 2, 0), ac_table_printable (gOb, ir, 7, "")) ;
	}
    }
  ac_free (gOb) ;

  if ((gOb = ac_tag_table (oKantor, "Pfam", 0)))
    {
      char *cp ;
      const char *ccp ;
      int nn = 0 ;

      for (ir = 0 ;ir < gOb->rows ; ir++)
	{
	  if (ir && /* avoid duplication if there are several hits to the same pfam */
	      !strcmp (ac_table_printable (gOb, ir, 0, ""), 
		       ac_table_printable (gOb, ir - 1, 0, "")
		       )
	      )
	    continue ;
	  obj = ac_table_obj (gOb, ir, 0, 0) ;
	  ccp = ac_tag_text (obj, "Definition", 0) ;
	  cp = ccp ? kantorCleanPfamTitle (ccp) : 0 ;

	  if (cp && /* avoid duplicated definitions */
	      (!vtxtPtr (pfamHintBlk) || !strstr (vtxtPtr (pfamHintBlk), cp)) &&
	      !pickMatch (cp, "*unknown*"))
	    {
	      if (nn++)
		vtxtPrint (pfamHintBlk, " and") ;
	      vtxtPrintf (pfamHintBlk, " %s", cp) ;
	    }
	  ac_free (obj) ;
	}  
    }
  ac_free (gOb) ;

  if (vtxtPtr (blastHintBlk))
    { /* export Blastp_title && TitleHint */
      extern BOOL titleBestSuggestionANALYZER (vTXT blkp, char *blastHints) ;
      isTitle = titleBestSuggestionANALYZER (blkp, vtxtPtr (blastHintBlk)) ;
    }
  if (!isTitle && vtxtPtr (pfamHintBlk))
    {
      vtxtPrintf (blkp, "Blastp_title \"%s\"\n", vtxtPtr (pfamHintBlk)) ;
      isTitle = 1 ;
    }

  if (!isTitle)
     vtxtPrintf (blkp, "Blastp_title \"Protein of unknown function\"\n") ;
  
  vtxtDestroy (blastHintBlk) ;
  vtxtDestroy (pfamHintBlk) ;
  
  return vtxtPtr (blkp) ? TRUE : FALSE ;
} /* kantorGetBlastPfamTitle */

/*****************************************************************************/
/* export the blast pfam title is available, or construct the UFO title */
static void kantorGetTaxBlastTitle (vTXT blkp, AC_OBJ oKantor, int isComplete)
{
  char *ancestor  ;
  char *cp =  0 ;

  if (ac_has_tag (oKantor, "Taxblast_date") &&
      (ancestor =  fogicArguableCommonAncestors (ac_name (oKantor))))
    {
      if (ac_has_tag (oKantor, "Archaea") ||
	  ac_has_tag (oKantor, "Bacteria"))
	cp = "of ancient origin" ;
      else if (pickMatch (ancestor, "*Eukaryota*"))
	cp =  "of eukaryotic origin" ;
      else if (pickMatch (ancestor, "*fungi*"))
	cp =  "of fungal and metazoan origin" ;
      else if (pickMatch (ancestor, "*Metazoa*"))
	cp =  "of metazoan origin" ;
      else if (!strcasecmp (ancestor, "Bilateria"))
	cp =  "of bilaterial origin" ;
      else if (pickMatch (ancestor, "*Vertebrata*"))
	cp =  "of vertebrate origin" ;
      else if (pickMatch (ancestor, "*Teleost*"))
	cp =  "of vertebrate origin" ;
      else if (pickMatch (ancestor, "*Sushi*"))
	cp =  "of vertebrate origin" ;
      else if (pickMatch (ancestor, "*Amphibi*"))
	cp =  "of vertebrate origin" ;
      else if (pickMatch (ancestor, "*Amniota*"))
	cp =  "of vertebrate origin" ;
      else if (pickMatch (ancestor, "*Eutheria*"))
	cp =  "of mammalian origin" ;
      else if (pickMatch (ancestor, "*Mammalia*"))
	cp =  "of mammalian origin" ;
      else if (pickMatch (ancestor, "*Homo sapiens*"))
	cp =  "human specific" ;
      else if (
	       (!strcasecmp (ancestor, "Pseudocoelomata") ||
		!strcasecmp (ancestor, "Caenorhabditis elegans")) &&
	       isComplete)
	cp =  "nematode specific" ; 
    }
  if (cp)
    vtxtPrintf (blkp, "TaxBlast_title \"%s\"\n", cp) ;
  return ;
} /*  kantorGetTaxBlastTitle */

/*****************************************************************************/

static BOOL kantorGetMultipleFamilyTitle (vTXT blkp, AC_OBJ oKantor, char *kantorName)
{
  int nSelf =  0;
  const char *ccp ;
  char *tag =  0 ;
  AC_TABLE tbl = 0 ;
  int n, ir ;

  if (Spc == WORM)
    tag =  "NCaenorhabditis_elegans" ;
  else if (Spc == HUMAN)
    tag =  "NHomo_sapiens" ;
  else if (Spc == ARA)
    tag =  "NArabidopsis_thaliana" ;
  else if (Spc == MOUSE)
    tag =  "NMus_musculus" ;
  else if (Spc == RAT)
    tag =  "NRattus_norvegicus" ;
  else if (Spc == DROSO)
    tag =  "NDrosophila_melanogaster" ;

  if (tag && (tbl =  ac_tag_table (oKantor, tag, 0)))
   nSelf =  ac_table_int (tbl, 0, 0, 0) ;
  ac_free (tbl) ;

  /* restrict to 2 hits among 6 best to self */
  if (nSelf  >=  4)
    {
      if (Spc == WORM)
	tag =  "*[Caenorhabditis elegans]*" ;
      else if (Spc == HUMAN)
	tag =  "*[Homo sapiens]*" ;
      else if (Spc == ARA)
	tag =  "*[Arabidopsis thaliana]*" ;
      else if (Spc == MOUSE)
	tag =  "*[Mus musculus]*" ;
      else if (Spc == RAT)
	tag =  "*[Rattus_norvegicus]*" ;
      else if (Spc == DROSO)
	tag =  "*[Drosophila melanogaster]*" ;

      n =  0 ; 
      if ((tbl =  ac_tag_table (oKantor, "Blastp", 0)))
	for (n = ir = 0; ir < tbl->rows && ir < 6 ; ir++)
	  if ((ccp =  ac_table_printable (tbl, ir, 7, 0)) &&
	      pickMatch (ccp, tag))
	    n++ ;
      if (n < 2)
	nSelf =  0 ;
      ac_free (tbl) ;
    }
  if (nSelf  >=  4)
    { 
      if (Spc == WORM)
	tag =  " C.elegans" ;
      else if (Spc == HUMAN)
	tag =  " human" ;
      else if (Spc == ARA)
	tag =  "n A.thaliana" ;
      else if (Spc == MOUSE)
	tag =  " mouse" ;
      else if (Spc == RAT)
	tag =  " mrat" ;
      else if (Spc == DROSO)
	tag =  " drosophila" ;
      
      vtxtPrintf (blkp, "Family_title ") ;
      vtxtPrintf (blkp, "\"member of a%s multigene family\"\n", tag) ;
      
      return TRUE ;
    }
  
  return  FALSE ;
} /* kantorGetMultipleFamily */

/*****************************************************************************/
/*****************************************************************************/
/********************************** public functions *************************/

char *kantorReportTitles (char *productName)
{
  AC_DB db =  ac_open_db (kantorTitleDB, 0) ;
  AC_OBJ oProduct ;
  char *ptr = 0 ;
  
  if (db)     
    {
      oProduct = ac_get_obj (db, "Product", productName, 0) ;
      if (oProduct)
	{
	  ptr = gtReportTitles (db, oProduct) ;
	  ac_free (oProduct) ;
	}
    }
  
  ac_db_close (db) ;
  return ptr ;
} /* kantorReportTitles */

/*****************************************************************************/

char *kantorGetTitles (vTXT blkp, char *kantorName, int isComplete)
{
  AC_DB db =  ac_open_db (kantorTitleDB, 0) ;
  AC_OBJ oKantor ;
  
  if (db)     
    {
      oKantor = ac_get_obj (db, "Kantor", kantorName, 0) ;
      
      Spc =  ficheDBCreatureID (db) ;
      if (oKantor)
	{
	  kantorGetBlastPfamTitle (blkp, oKantor, kantorName) ;
	  kantorGetTaxBlastTitle (blkp, oKantor, isComplete) ;  
	  kantorGetMultipleFamilyTitle (blkp, oKantor, kantorName) ;
	  
	  ac_free (oKantor) ;
	}       
      ac_db_close (db) ;
    }
  return vtxtPtr (blkp) ;
} /* kantorGetTitles */

/*****************************************************************************/
/*****************************************************************************/

