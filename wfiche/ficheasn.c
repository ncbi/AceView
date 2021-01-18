#include "../wfiche/biolog.h"
#include "../wfiche/gtitle.h"
#include "utils.h"
#include "dna.h"
#include "peptide.h"
#include "parse.h"
#include "session.h"
#define KANS_LIMIT 1024+1

#define TESTZZZ
/*  SHOW_ALL_ASN_GENE was set iin dec 2002 dump, Kim ask it to be removed in mars 2003 */
#define SHOW_ALL_ASN_GENE 0
/* june 9 2003, Taatian asks that i restore the genes onto the NMs */
#define SHOW_GENE_ON_NM 1
#define SHOW_GENE_ON_NP 0

static BOOL clipGenomic256 = FALSE ;
static char *ficheAsnGeneSeqFeat (vTXT blkp, GMP *gmp, int gStart, int gEnd, 
				  const char *seqName, BOOL xref, BOOL pseudo) ;
static char *ficheAsnProductRef (vTXT blkp, GMP *gmp, AC_OBJ oProduct, BOOL isReal) ;
static void ficheAsnChromoAnnot (vTXT blkp, AC_DB db, AC_OBJ oMap, char *testGene, char style) ;
static void ficheDescrWormSource (vTXT blkp) ;
static int iGlimit = 0 ;

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  ASN Supporting functions 
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static const char* niceChromoName (AC_OBJ oMap)
{
  const char *cp = ac_name (oMap) ;
  if (!strncasecmp ("Chromosome_", cp, 11))
    cp += 11 ;
  return cp ;
}

static int isRefSeqDump = 0 ;
static char* ficheAsnId (const char *id, char *endding)
{
  char *nc = 0 ;
  static char buf[1000] ; 

  if (!endding) 
    endding = "" ;

  if (strlen (id) > 900)
    messcrash ("id too long in  ficheId, sorry::", id) ;
 
  switch (isRefSeqDump)
    {
    case 0:
      sprintf (buf, "local str \"%s%s\"", id, endding) ;
      break ;
    case 1:
      if (!strcmp (id, "CHROMOSOME_I"))
	nc = "NC_003279" ;
      else if (!strcmp (id, "CHROMOSOME_II"))
	nc = "NC_003280" ;
      else if (!strcmp (id, "CHROMOSOME_III"))
	nc = "NC_003281" ;
      else if (!strcmp (id, "CHROMOSOME_IV"))
	nc = "NC_003282" ;
      else if (!strcmp (id, "CHROMOSOME_V"))
	nc = "NC_003283" ;
      else if (!strcmp (id, "CHROMOSOME_X"))
	nc = "NC_003284" ;

      if (nc)
	sprintf (buf, "other  { accession \"%s\" }", nc) ;
      else
	{
	  /* replace global by a local id, june 27, 2003 
	   * sprintf (buf, "general { db \"WormGenes\", tag str  \"%s%s\" }"
	   *             , id, endding) ; 
	   */
	  sprintf (buf, "local str \"WormGenes:%s%s\" ", id, endding) ;
	}
      break ;
    }
  return buf ;  
}

/************************************************************************/

static char *ficheAsnDate (vTXT blkp)
{
  time_t ltm ;struct tm *t ;

  vtxtClear (blkp) ;
  time ( &ltm) ;
  t = localtime (&ltm) ;
  
  vtxtPrintf (blkp, 
	    "std {"
	    "  year %i , " 
	    "  month %i , " 
	    "  day %i "
	    "}"
	    , t->tm_year+1900, t->tm_mon+1, t->tm_mday) ;

  return vtxtPtr (blkp) ;
} /* ficheAsnDate */

static void printfAsnAuthor (vTXT  blk, const char * last, const char * first, const char * initial)
{
  vtxtPrintf (blk, 
	    "{ name"
	    "  name {"
	    "    last \"%s\" , "
	    "    first \"%s\" , "
	    "    initials \"%s\""
	    "} }"
	    , last, first, initial) ;
}

static void printfAsnDate (vTXT  blkp)
{
  time_t ltm ;struct tm *t ;
  time ( &ltm) ;
  t = localtime (&ltm) ;
  
  vtxtPrintf (blkp, 
	    "{"
	    "  year %i , " 
	    "  month %i , " 
	    "  day %i "
	    "}"
	    , t->tm_year+1900, t->tm_mon+1, t->tm_mday) ;
}

void ficheTiledFrom (vTXT  blk, AC_OBJ  oMrna)
{
  AC_TABLE gTil ;
  AC_OBJ oSeq ;
  int 	ir ;
  const char * ptr ;
  AC_HANDLE h = handleCreate () ;

  if ((gTil = ac_tag_table (oMrna, "Tiling_path", h)))
    {
      for (ir = 0 ;ir <gTil->rows ;ir++)
	{
	  oSeq = ac_table_obj (gTil, ir, 2, h) ;
	  ptr = ac_name (oSeq) ;
	  if (!strncmp (ptr, "GGB:", 4))ptr++ ;
	  if (!strncmp (ptr, "GB:", 3))ptr+= 3 ;
	  vtxtPrintf (blk, "%s%s", ir == 0 ? "" : ", ", ptr) ;
	}
    }
  ac_free (h) ;
}

/*****************************************************************************/
/*****************************************************************************/

static void ficheAsnMapSubtype (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = handleCreate () ;
  AC_TABLE gMap ; 
  AC_OBJ oMap = ac_tag_obj (gmp->gene, "Intmap", h) ;
  const char *tNam = niceChromoName (oMap) ;
  float mpos = 0 ;

  vtxtPrintf (blkp, 
	    "    {"
	    "      subtype chromosome, "
	    "      name \"%s\""
	    "    }"
	    , tNam) ;
  
  { /* all genetic positions in a single line */
    vTXT vtxtMap = vtxtCreate () ;
    int nmap = 0 ;
    char *ptr ;

    if ((gMap = ac_tag_table (gmp->gene, "InterpolatedMap", h)) &&
	(mpos = ac_table_float (gMap, 0, 1, -9999)) != -9999)
      {
	vtxtPrintf (vtxtMap,
		   "%s%s;%s%.2f cM (interpolated genetic position)"
		    , nmap++ ? "; " : ""
		    , tNam, mpos > 0 ? "+" : "", mpos
		) ;		  
      }
    
    if ((gMap = ac_tag_table (gmp->gene, "Map", h)) &&
	(mpos = ac_table_float (gMap, 0, 2, -9999)) != -9999)
      {
	vtxtPrintf (vtxtMap,
		    "%s%s;%s%.2f cM (measured genetic position)"
		    , nmap++ ? "; " : ""
		    , tNam, mpos > 0 ? "+" : "", mpos
		) ;		  
      }

    if ((ptr = vtxtPtr (vtxtMap)))
       {
	 vtxtPrintf (blkp, 
		     "    , {"
		    "      subtype map, "
		     "      name \"%s\""
		     "    }"
		     , ptr
		     ) ;		  
       }
    vtxtDestroy (vtxtMap) ;
  }
  
  if (gmp->style == 'r' && 
      (gMap = ac_tag_table (gmp->gene, "IntMap", h)))
    {   
      int nn = ac_table_int (gMap, 0, 2, 0) - ac_table_int (gMap, 0, 1, 0) ;
      if (nn < 0) nn = - nn  ;
      nn++ ; /* include the extremities */
      vtxtPrintf (blkp, 
		"    , {"
		"      subtype map, "
		"      name \"%s; covering %d bp, from base %d to %d on genome release %s\""
		"    }"
		, tNam
		, nn
		  , ac_table_int (gMap, 0, 1, 0)
		  , ac_table_int (gMap, 0, 2, 0)
		  , genomeRelease
		) ;		 
    }
  ac_free (h) ;
} /* ficheAsnMapSubtype */

/*****************************************************************************/

static void ficheAsnTilingClones (vTXT blkp, GMP *gmp)
{
  AC_OBJ oRead, oClone ;
  AC_HANDLE h = handleCreate () ;
  AC_TABLE gTiling = ac_tag_table (gmp->mrna, "Mrna_covered_by", h) ;
  int ir ;
  DICT *dict = dictCreate (12) ;

  for (ir = 0 ; gTiling && ir < gTiling->rows ; ir++)
    {
      oRead = ac_table_obj (gTiling, ir, 0, h) ;
      oClone = ac_tag_obj (oRead, "cDNA_clone", h) ;
      if (!oClone || dictFind (dict, ac_name(oClone), 0))
	continue ;
      dictAdd (dict, ac_name(oClone), 0) ;
      vtxtPrintf (blkp, "%s%s",ir > 0 ? ", " : "", gtYkName (ac_name(oClone))) ;
    }
  dictDestroy (dict) ;
  ac_free (h) ;
} /* ficheAsnTilingClones  */

/*****************************************************************************/

static void ficheAsnCloneSubtype (vTXT blkp, GMP *gmp, AC_OBJ oGene, AC_OBJ oMrna, AC_OBJ oGF, char *asnTag)
{  
  char *ptr ;
  vTXT buf ; buf = vtxtCreate () ;

  if (gmp->style == 'r')
    {
      fichePrimersParagraphContent (buf, gmp) ;
      ficheCloneParagraphContent (buf, gmp) ;
    }
  else if (gmp->mrna && gmp->style == 's')
    {
      ficheAsnTilingClones (buf, gmp) ;
    }
  if ((ptr = vtxtPtr (buf)))
    {
      if (asnTag) vtxtPrintf (blkp, "%s", asnTag) ;
      vtxtPrintf (blkp, 
		  "{"
		  "      subtype clone, "
		  "      name \""
		  ) ;
      vtxtPrintf (blkp, "%s", ptr) ;
      vtxtPrintf (blkp, 
		  "\""
		  "}"
		  ) ;
    }
  vtxtDestroy (buf) ;
} /* ficheAsnCloneSubtype */

/*****************************************************************************/

typedef struct clo2libstruct { AC_OBJ clo, lib ; int idx ; } C2LIB ;
static int C2libOrder (const void *a, const void *b)
{
  const C2LIB *va = (const C2LIB *)a , *vb = (const C2LIB *)b ;
  const char *la, *lb ;

  if (va->idx != vb->idx)
    {
      return va->idx - vb->idx ;
    }
  else
    {
      la = ac_name (va->lib) ; lb = ac_name (vb->lib) ;
      return strcmp (la, lb) ;
    }
}

/* declare the lib here in the order you wish them in the output */
static DICT *ficheCloneLibDict (void)
{
  static DICT *dict = 0 ;

  if (dict)
    return dict ;
  
  dict = dictCreate (12) ;
  dictAdd (dict, "yke", 0) ;
  dictAdd (dict, "cee", 0) ;
  dictAdd (dict, "ycL1", 0) ;
  dictAdd (dict, "ycL2", 0) ;
  dictAdd (dict, "ycL4", 0) ;
  
  dictAdd (dict, "yc", 0) ;

  dictAdd (dict, "cm", 0) ;
  dictAdd (dict, "cem", 0) ;
  
  dictAdd (dict, "Barstead", 0) ;
  dictAdd (dict, "Kim", 0) ;
  dictAdd (dict, "Okkema", 0) ;
  dictAdd (dict, "cem", 0) ; 
  
  dictAdd (dict, "ykm", 0) ;
  dictAdd (dict, "ykms", 0) ;
  
  dictAdd (dict, "sperm", 0) ;

  dictAdd (dict, "ovary", 0) ;
  dictAdd (dict, "mv2h", 0) ;
  dictAdd (dict, "rare", 0) ; 
  dictAdd (dict, "yd", 0) ;

  dictAdd (dict, "cadmium", 0) ;

  return dict ;
}

static void ficheAsnCloneLibSubtype (vTXT blkp, GMP *gmp, char *asnTag)
{
  int ir, jr, idx ;
  const char *ccp ;
  AC_TABLE gTmp=0 ;
  AC_OBJ oClone, oLib ;
  C2LIB *c2lib ;
  Array C2lib = 0 ;
  DICT *dict = ficheCloneLibDict () ; /* do not destroy */
  vTXT buf = vtxtCreate () ;
  AC_HANDLE h = handleCreate () ;

  if (0 &&    /* aug 25, 2005, put it back here rather than in the comment */
      gmp->style == 's') /* july 16, 2003, put it rather in the comment */
    return ;

  /* export only if the clones can be attributed to this variant */
  if (!gmp->pg && !gmp->mrna)
    return ;
  if (gmp->gene && gmp->pg && !gmp->tg)
    {
      gTmp = ac_tag_table (gmp->gene, "genefinder", h) ;
      if (gTmp->rows > 1)
	return ;
      gTmp = ac_tag_table (gmp->gene, "Has_cDNA_clone", h) ;
    }
  else if (gmp->mrna)
    gTmp = ac_tag_table (gmp->mrna, "cDNA_clone", h) ;
  if (!gTmp || !gTmp->rows)
    return ;

  if (asnTag) vtxtPrintf (blkp, "%s", asnTag) ;
  vtxtPrintf (blkp, 
	     "{"			 
	     "   subtype clone-lib, "
	     "   name \""
	     ) ;

  if (gTmp && gTmp->rows >= 1)
    {
      C2lib = arrayCreate (gTmp->rows, C2LIB) ;
      for (ir = jr = 0 ; ir < gTmp->rows ; ir++)
	{
	  oClone = ac_table_obj (gTmp, ir, 0, h) ;
	  oLib = oClone ? ac_tag_obj (oClone, "Library", h) : 0 ;
	  
	  ccp = oLib ? ac_name(oLib) : "toto" ;
	  dictAdd (dict, ccp, &idx) ;
	  c2lib = arrayp (C2lib, ir, C2LIB) ;
	  c2lib->clo = oClone ;
	  c2lib->lib = oLib ;
	  c2lib->idx = idx ;
	}

      arraySort (C2lib, C2libOrder) ;

      for (ir = jr = 0 ; ir < arrayMax (C2lib) ; ir++)
	{
	  c2lib = arrp (C2lib, ir, C2LIB) ;
	  oClone = c2lib->clo ;
	  oLib = c2lib->lib ;
	  if (!ir || 
	      (
	       ac_name(oLib) != ac_name((c2lib-1)->lib) &&
	       ( strcasecmp (ac_name(oLib) , "ykms") || strcasecmp (ac_name((c2lib-1)->lib) , "ykm"))
	       )) /* mix up ykm ykms in the export */
	    {
	      jr = 0 ;
	      vtxtPrintf (buf, "%s%s:"
			  , ir > 0 ? ";\n" : ""
			  , oLib ? ac_tag_printable (oLib, "Title", ac_name (oLib)) : "") ;
	    }
	    
	  vtxtPrintf (buf, "%s%s"
		      , jr++ > 0 ? ", " : " "
		      , gtYkName (ac_name (oClone))) ;
	}
    }

  if (vtxtPtr (buf))
    vtxtPrintf (blkp, "%s", vtxtPtr (buf)) ;
  vtxtPrintf (blkp, 
	     "\""
	     "}"
	     ) ;

  vtxtDestroy (buf) ;
  arrayDestroy (C2lib) ;
  ac_free (h) ;

  return ;
} /* ficheAsnCloneLibSubtype */

/*****************************************************************************/

static void ficheAsnGenePublications (vTXT blkp, GMP *gmp, AC_OBJ oGene, char *asnTag)
{
  int 	npms, ir, dummy ;
  char dummychar ;
  AC_OBJ oPap ;
  AC_TABLE gPubs ;
  const char *pmTxt ;
  Array pms = 0 ;
  extern int arrstrcmp(const void *s1, const void *s2) ;
  AC_HANDLE h = handleCreate () ;

  /* 
   * we had at some point
   gPubs = ac_tag_table (oGene, gmp->style == 's' ? "Genbank_reference" : "Reference", h)) ;
   the idea was to subselect the paper to quote in the submission
  */

  if ((gPubs = ac_tag_table (oGene, "Reference", h)))
    {
      npms = 0 ;  pms = arrayCreate (20, char *) ;
      for (ir = gPubs->rows - 1 ; ir >= 0 ; ir--)
	{
	  oPap = ac_table_obj (gPubs, ir, 0, h) ;
	  pmTxt = ac_tag_printable (oPap, "PMID", 0) ;
	  if (!pmTxt && !strncasecmp (ac_name(oPap), "pm", 2) &&
	      sscanf (ac_name(oPap) + 2, "%d%c", &dummy, &dummychar) == 1)
	    pmTxt = ac_name(oPap) + 2 ;	    
	  /* jump direct submission GB* pseudo citations and non pmid papers */
	  if (!pmTxt)
	    continue ;
	  array (pms, npms++, const char *) = pmTxt ;
	}
      arraySort (pms, arrstrcmp) ;
      arrayCompress (pms) ;

      if (arrayMax (pms))
	vtxtPrintf (blkp, asnTag) ;
      for (ir = 0 ; ir < arrayMax (pms) ; ir++)
	{
	  pmTxt =  array (pms, ir, const char *) ;
	  

	  vtxtPrintf (blkp
		      , " %s { pub {  pmid %s  }}"
		      , ir  ? ", pub ": ""
		      , pmTxt
		    ) ;
	}
    }
  ac_free (h) ;
  arrayDestroy (pms) ;
} /* ficheAsnGenePublications */

/*****************************************************************************/

static void ficheGeneDbXref (vTXT blkp, GMP *gmp, AC_OBJ oGene, AC_OBJ oProduct, char *asnTag)
{ /* nuc-prot set dbxref */
  int ir ;
  const char *ptr ;
  AC_TABLE oGenefinder, oLL ;
  AC_OBJ oGf ;
  BOOL allXref = FALSE ;
  AC_HANDLE h = handleCreate () ;

  if (asnTag)
    vtxtPrintf (blkp, "%s {", asnTag) ; 
    
  vtxtPrintf (blkp, 
	      "  { db \"%s\" , tag str \"%s\"  }" 
	      , gmp->Spc == WORM ? "AceView/WormGenes" : "HumanGenes", ac_name (oGene)
	      ) ;
  if ((oLL = ac_tag_table (oGene, "LocusId", h)))
    {
      const char *best = 0 ;
      for (ir = 0 ; ir < oLL->rows ; ir++)
	{
	  if ((ptr = ac_table_printable (oLL, ir, 0, 0)) &&
	      (!best || lexstrcmp (ptr, best) < 0))
	    best = ptr ;
	}
      if (best)
	vtxtPrintf (blkp, 
		    "  , { db \"LocusID\" , tag str \"%s\"  }" 
		    , best
		    ) ;
    }
  
  if (allXref && gmp->style == 'r' && 
      gmp->Spc == WORM  && 
      (oGenefinder = ac_tag_table (oGene, "Genefinder", h)))
    {
      for (ir = 0 ; ir < oGenefinder->rows ; ir++)
	{
	  oGf = ac_table_obj (oGenefinder, ir, 0, h) ;
	  if (strstr (ac_name (oGf), ".ws6"))
	    continue ;
	  vtxtPrintf (blkp, 
		      "   , { db \"WormBase\" , tag str \"%s\"  }" 
		      , ac_name (oGf)
		      ) ;
	}
    }
  
  if (((allXref && gmp->style == 'r')  || gmp->style == 's') &&
      gmp->Spc == WORM &&  
      (ptr = gtCloneGroup (oGene, TRUE)))
    vtxtPrintf (blkp, 
		"  , {db \"NextDB\", tag str \"%s\"  }"
		, ptr
		) ;
  
  if (allXref && gmp->style == 'r' &&
      gmp->Spc == WORM &&  
      (ptr = gtOst (oGene)))
    vtxtPrintf (blkp, 
		"  , {db \"WorfDB\", tag str \"%s\"  }"
		, ptr
		) ;

  if (0 &&  /* this site does not work */
      allXref && gmp->style == 'r' &&
      gmp->Spc == WORM &&  
      (ptr = gtHyman (oGene)))
    vtxtPrintf (blkp, 
		"  , {db \"HymanRNAiDB\", tag str \"%s\"  }"
		, ptr
		) ;
  /* */
  if (asnTag)
    vtxtPrintf (blkp, " }") ;
  ac_free (h) ;
} /* ficheGeneDbXref */

/*****************************************************************************/

static void ficheMrnaDbXref (vTXT blkp, GMP *gmp, AC_OBJ oGene, AC_OBJ oProduct, char *asnTag)
{ /* nuc-prot set dbxref */
  int ir, nn = 0 ;
  AC_TABLE gTag ;
  AC_OBJ oTmp ;
  AC_HANDLE h = 0 ;

  if (1)  /* to please Kim */
    return ;
  h = handleCreate () ;
/*
* This whole function is dead code
*/

  if (1)
    {
      nn++ ;
      vtxtPrintf (blkp, "%s { ", asnTag) ; 
      ficheGeneDbXref (blkp, gmp, oGene, 0, 0) ;
    }

  if (gmp->style == 'r' &&
      (gTag = ac_tag_table (oGene, "COG", h)))
    {
      if (!nn++)
	{ vtxtPrintf (blkp, asnTag) ; vtxtPrintf (blkp, "{"); }

      for (ir = 0 ;ir <gTag->rows ;ir++)
	{
	  oTmp = ac_table_obj (gTag, ir, 0, h) ;
	  vtxtPrintf (blkp, 
		     "  %s {db \"COG\", tag str \"%s\"  }"
		     ,nn++ > 1 ? "," : "" , ac_name (oTmp)
		     ) ;
	}
    }
  
  /* i kill those db XREF bnecause they are repaeted in the misc-feature section */
  if (0 && (gTag = ac_tag_table (oGene, "GO", h)))
    {
      if (!nn++)
	{ vtxtPrintf (blkp, asnTag) ; vtxtPrintf (blkp, "{"); }

      for (ir = 0 ;ir <gTag->rows ;ir++)
	{
	  oTmp = ac_table_obj (gTag, ir, 0, h) ;
	  vtxtPrintf (blkp, 
		     "  %s {db \"GO\", tag str \"%s\"  }"
		     ,nn++ > 1 ? "," : "" , ac_name (oTmp)
		     ) ;
	}
    }
  
  if (0 && (gTag = ac_tag_table (oGene, "Pfam_GO", h)))
    {
      if (!nn++)
	{ vtxtPrintf (blkp, asnTag) ; vtxtPrintf (blkp, "{"); }
     
      for (ir = 0 ;ir <gTag->rows ;ir++)
	{
	  const char * ccp = gTag->cols >= 3 ? ac_table_printable (gTag, ir, 2, "") : 0 ;
	  if (ccp && strlen (ccp)>3 && !strncmp (ccp, "GO:", 3))
	    vtxtPrintf (blkp, 
		       "  %s {db \"GO\", tag str \"%s\"  }"
		       ,nn++ > 1 ? "," : "" , ccp+3
		       ) ;
	}
    }
  
  if (0 && oProduct && (gTag = ac_tag_table (oProduct, "Pfam", h)))
    {
      AC_OBJ oPfam ;
      const char *tmpNam = "" ;
      if (!nn++)
	{ vtxtPrintf (blkp, asnTag) ; vtxtPrintf (blkp, "{"); }
     
      for (ir = 0 ;ir < gTag->rows ;ir++)
	{
	  if (!strcmp (ac_table_printable (gTag, ir, 0, ""), tmpNam))
	    continue ;
	  oPfam = ac_table_obj (gTag, ir, 0, h) ;
	  tmpNam = ac_tag_printable (oPfam, "accession", "") ;
	  if (!*tmpNam) continue ;
	  if (!strcasecmp (ac_table_tag (gTag, ir, 1, ""), "pfam"))
	    vtxtPrintf (blkp, 
		       "  %s {db \"CDD\", tag str \"pfam%s\"  }"
		       ,nn++ > 1 ? "," : "" , tmpNam+2
		       ) ;
	}
    }
  
  if (nn)
    vtxtPrintf (blkp, " }") ;
  ac_free (h) ;
} /* ficheMrnaDbXref */

/*****************************************************************************/

void ficheAlternativeExons (vTXT  blk, AC_OBJ  oMrna)
{
  /*
    vtxtPrintf (blk, "This variant ") ;
    if (numIntrons >0)vtxtPrintf (blk, "%d", numIntrons) ;
    else vtxtPrintf (blk, "no") ;
    vtxtPrintf (blk, " confirmed intron%s", _multi (numIntrons)) ;
    
    numIntrons = ac_tag_int (oTranscribed_gene, "Nb_confirmed_alternative_introns", 0) ;
    if (numIntrons >1)
    vtxtPrintf (blk, ", %d of which %s alternative", numIntrons, _isare (numIntrons)) ;
    
    vtxtPrintf (blk, ". ") ;
    
    
    AC_OBJ   gTil, oSeq ;
    int 	ir ;
    char * ptr ;
    
    if ((gTil = ac_tag_table (oMrna, "Tiling_path", h)))
    {
    for (ir = 0 ;ir <gTil->rows ;ir++)
    {
    oSeq = ac_tag_table (gTil, ir, 2) ;
    ptr = ac_name (oSeq) ;
    if (!strncmp (ptr, "GGB:", 4))ptr++ ;
    vtxtPrintf (blk, "%s%s", ir == 0 ? "" : ", ", ptr) ;
    }
    }*/
}

static void ficheAsnPFAM (vTXT blkp, AC_OBJ oProduct, char *prRealName, char *sPeptide)
{
  const char *ptr ;
  char *cp ;
  int ir, x1, x2, pepL = strlen (sPeptide) ;
  AC_TABLE gPsortTag ;
  AC_OBJ oDom ;
  char buf[256] ;  
  AC_HANDLE h = handleCreate () ;
  
  /*
    please do not remove this code
    we may activate it
    it validate sin asn but produces no
    output in flat gb
    so the same data is now reported inside
    gtProductDescriptor->productBasePfamHomol ()
    
    *****    site ENUMERATED {
    modified (5) -------------------------- diff with genome
    myristoylation (7) --------------------, N_myristoylation_domain
    signal-peptide (23) , ------------------N_terminal_signal_domain but not N_terminal_domain 
    transit-peptide (24) , 
    transmembrane-region (25) , ------------- Transmembrane_domain
    other (255) } ,    ---------------------les pfam 

  **********************/
 
  if ((gPsortTag = ac_tag_table (oProduct, "Pfam", h)))
    {
      float score = 0 ;
      for (ir = 0 ; ir < gPsortTag->rows && gPsortTag->cols >= 8 ; ir++)
	{
	  oDom = ac_table_obj (gPsortTag, ir, 0, h) ;
	  if (ir > 0 &&
	      !strcmp (ac_name(oDom), ac_table_printable (gPsortTag, ir - 1, 0, "")))
	    {
	      int a1, a2, b1, b2 ;

	      a1 = ac_table_int (gPsortTag, ir - 1, 3, 0) ;
	      a2 = ac_table_int (gPsortTag, ir - 1, 4, 0) ;
	      b1 = ac_table_int (gPsortTag, ir, 3, 0) ;
	      b2 = ac_table_int (gPsortTag, ir, 4, 0) ;
	      if (a1 < b1) a1 = b1 ;
	      if (a2 > b2) a2 = b2 ; 
	      if (4 * (a2 - a1) > 3 * (b2 - b1)) continue ; /* ignore this repetion */
	    }
	      
	  x1 = (ac_table_int (gPsortTag, ir, 3, 0) -1)/3 ;
	  x2 = (ac_table_int (gPsortTag, ir, 4, 0) -1)/3 ;
	  if (x1 < 0) x1 = 0 ;
	  if (x2 > pepL - 1) x2 = pepL - 1 ;
	  if (x2 < x1 + 1) continue ;
	  /*
	    #ifdef ___jonathan_recommendation
	    {
	    data
	    region "ABC transporter transmembrane region. This family
	    represents a unit of six transmembrane helices. Many members of the ABC
	    transporter family (pfam00005) have two such regions" ,
	    comment "ABC_membrane" ,
	    
	    } ,
	    current export: 
	    {
	    data
	    imp {
	    key "PFAM" } ,
	    comment "[Pfam/InterPro description] ankyrin repeat" ,
	    
	    } 
	    
	    #endif
	  */
	  vtxtPrintf (blkp, 
		      ", "
		      "{" 
		      "  data region \""
		      ) ;
	  if ((ptr = ac_tag_printable (oDom, "Definition", ac_name(oDom))))
	    vtxtPrintf (blkp, "[Pfam/InterPro description] %s", gtCleanUp(ptr)) ;
	  if (0 &&
	      (
	       ir == 0 ||
	       strcmp (ac_name(oDom), ac_table_printable (gPsortTag, ir - 1, 0, ""))
	       ))
	    {
	      if ((ptr = ac_tag_printable (oDom, "Comment", 0)))
		vtxtPrintf (blkp, ": %s.", gtCleanUp (ptr)) ;
	    }
	  
	  score = ac_table_float (gPsortTag, ir, 2, -1) ; ptr = ac_table_printable (gPsortTag, ir, 7, "") ;
	  if (ptr) ptr = strstr (ptr, "Evalue") ;
	  if (ptr)								
	    vtxtPrintf (blkp, " HMMER score %.2f %s", score, ptr) ;
	  vtxtPrintf (blkp, 
		      "\""
		      );
	  if (0) vtxtPrintf (blkp,", comment \"%s\"",  ac_name(oDom)) ;
	  
	  vtxtPrintf (blkp, 
		      " , location int {" 
		      "    from %i , " 
		      "    to %i , " 
		      "    strand plus , " 
		      "    id %s" 
		      "  }" 
		      , x1, x2
		      , ficheAsnId (prRealName, 0)) ;
	  
	  if ((ptr = ac_tag_printable (oDom, "Accession", 0)) && !strncasecmp(ptr,"PF",2))
	    vtxtPrintf (blkp, ", dbxref {{ db \"CDD\", tag str \"pfam%s\" }} ",  ptr+2) ;
	  
	  vtxtPrintf (blkp, 
		      "}" ) ;
#ifdef __code_giving_nice_PFAM_tags
	  vtxtPrintf (blkp, 
		      ", "
		      "{" 
		      "  data imp {"
		      "    key \"PFAM\" " 
		      "  }, " 
		      "  comment \""
		      ) ;
	  if ((ptr = ac_tag_printable (oDom, "Definition", ac_name(oDom))))
	    vtxtPrintf (blkp, "[Pfam/InterPro] %s", ptr) ;
	  if (ir == 0 ||
	      strcmp (ac_name(oDom), ac_name(ac_table_printable (gPsortTag, ir - 1, 0, "")))
	    {
	      if ((ptr = ac_tag_printable (oDom, "Comment", 0)))
		vtxtPrintf (blkp, ": %s.", gtCleanUp (ptr)) ;
	    }
	  
	  score = ac_table_float (gPsortTag, ir, 2, -1) ; ptr = ac_tag_printable (gPsortTag, ir, 7, 0)) ;
	  if (ptr) ptr = strstr (ptr, "Evalue") ;
	  if (ptr)								
	    vtxtPrintf (blkp, " HMMER score %.2f %s", score, ptr) ;
	  
	  vtxtPrintf (blkp, 
		      "\", "
		      "  location int {" 
				"    from %i , " 
		      "    to %i , " 
		      "    strand plus , " 
		      "    id %s" 
		      "  }" 
		      , x1, x2
		      , ficheAsnId (prRealName, 0)) ;

	  if ((ptr = ac_tag_printable (oDom, "Accession", 0)))
	    vtxtPrintf (blkp, ", dbxref {{ db \"Pfam\", tag str \"%s\" }} ",  ptr) ;
	  
	  vtxtPrintf (blkp, 
		      "}" ) ;
#endif
	}
    }
  /* Psort_domain
   * because they come in 2 blocks: Nucleic_acid_binding && Psort_domain, 
   * i grab them from one layer up : Homol
   */
  if ((gPsortTag = ac_tag_table (oProduct, "Homol", h)))
    {
      const char *oType, *oDom ;
      char *cq ;
      int jr ;
      
      for (ir = jr = 0 ; ir < gPsortTag->rows && gPsortTag->cols >= 8 ; ir++)
	{
	  oType = ac_table_printable (gPsortTag, ir, 0, 0) ;
	  if (!oType ||
	      ( strcmp (oType, "Nucleic_acid_binding") &&
		strcmp (oType, "Psort_domain")))
	    continue ;
	  oDom = ac_table_tag (gPsortTag, ir, 1, "") ;
	  if (jr++)
	    {
	      int a1, a2, b1, b2 ;
	      const char *oDom1 =  ac_table_tag (gPsortTag, ir - 1, 1, "") ;
	      a1 = ac_table_int (gPsortTag, ir - 1, 4, 0) ;
	      a2 = ac_table_int (gPsortTag, ir - 1, 5, 0) ;
	      b1 = ac_table_int (gPsortTag, ir, 4, 0) ;
	      b2 = ac_table_int (gPsortTag, ir, 5, 0) ;
	      if (a1 < b1) a1 = b1 ;
	      if (a2 > b2) a2 = b2 ; 
	      if (!strcmp (oDom, oDom1) &&
		  4 * (a2 - a1) > 3 * (b2 - b1)) continue ; /* ignore this repetion */
	    }
	  /*
	   * here we can eliminate specific domains from the dump
	   * we wanted that up to mars 30 2003,
	   * but taht day we found that the intersect between pfam search and psort search of 
	   * those was poor (2/3 psort found by pfam), so we decided again to export them
	   if (!strcmp (ac_name(oDom), "Ribosomal_protein_domain") ||
	   !strcmp (ac_name(oDom), "Zinc_finger_domain") ||
	   !strcmp (ac_name(oDom), "RNA_binding_domain") 
	   )
	   continue ;
	  */
	  x1 = (ac_table_int (gPsortTag, ir, 4, 0) -1)/3 ;
	  x2 = (ac_table_int (gPsortTag, ir, 5, 0) -1)/3 ;
	  if (x1 < 0) x1 = 0 ;
	  if (x2 > pepL - 1) x2 = pepL - 1 ;
	  if (x2 < x1 + 1) continue ;
	  vtxtPrintf (blkp, 
		      ", "
		      "{" 
		      "  data region \"[PSORT] "
		      ) ;

	  if ((cq = gtCleanUp (oDom)))
	    { /* tmp copy */ 
	      cp = cq - 1 ;
	      while (*++cp) 
		if (*cp == '_') 
		  {
		    if (*(cp-1) == 'N')
		      *cp = '-' ;
		    else if (*(cp+1) == '4')
		      *cp = '-' ;
		    else
		      *cp = ' ' ;
		  }
	      vtxtPrintf (blkp, cq) ;
	    }	  
	  
	  if (x2 - x1 < 250 && x2 >= x1 && x1 >= 1 && x2 < strlen (sPeptide))
	    {
	      strncpy (buf, sPeptide + x1, 255) ;
	      buf[x2 - x1 + 1] = 0 ;
	      vtxtPrintf (blkp, ": %s", buf) ;
	    }
	  vtxtPrintf (blkp, 
		      "\""
		      );
	  if (0) vtxtPrintf (blkp,", comment \"%s\"", oDom) ;
	  
	  vtxtPrintf (blkp, 
		      " , location int {" 
		      "    from %i , " 
		      "    to %i , " 
		      "    strand plus , " 
		      "    id %s" 
		      "  }" 
		      , x1, x2
		      , ficheAsnId (prRealName, 0)) ;
	  
	  
	  vtxtPrintf (blkp, 
		      "}" ) ;
	}
    }

#if 0
  if ((gPsortTag = ac_tag_table (oProduct, "Transmembrane_domain", h)))
    {
      for (ir = 0 ;ir <gPsortTag->rows ;ir++)
	{
	  oDom = ac_table_obj (gPsortTag, ir, 2, h) ;
	  vtxtPrintf (blkp, 
		      ", "
		      "      {"
		      "        data site transmembrane-region, "
		      "        location int {"
		      "      	   from %d , "
		      "          to %d , "
		      "          id %s"
		      "        }"
		      "      }"
		      , (acInt (oDom)-1)/3, (acInt (oDom+1)-1)/3, ficheAsnId (prRealName, 0)) ;
	}
    }

  if ((gPsortTag = ac_tag_table (oProduct, "N_myristoylation_domain", h)))
    {
      for (ir = 0 ;ir <gPsortTag->rows ;ir++)
	{
	  oDom = ac_table_obj (gPsortTag, ir, 2, h) ;
	  vtxtPrintf (blkp, 
		      ", "
		      "      {"
		      "        data site myristoylation, "
		      "        location int {"
		      "      	   from %d , "
		      "          to %d , "
		      "          id %s"
		      "        }"
		      "      }"
		      , (acInt (oDom)-1)/3, (acInt (oDom+1)-1)/3, ficheAsnId (prRealName, 0)) ;
	}
    }
  if ((gPsortTag = ac_tag_table (oProduct, "N_terminal_signal_domain", h)))
    {
      for (ir = 0 ;ir <gPsortTag->rows ;ir++)
	{
	  oDom = ac_table_obj (gPsortTag, ir, 2, h) ;
	  vtxtPrintf (blkp, 
		      ", "
		      "      {"
		      "        data site signal-peptide, "
		      "        location int {"
		      "      	   from %d , "
		      "          to %d , "
		      "          id %s"
		      "        }"
		      "      }"
		      , (acInt (oDom)-1)/3, (acInt (oDom+1)-1)/3, ficheAsnId (prRealName, 0)) ;
	}
    }
#endif
  ac_free (h) ;
}

/******************************************************************/
/*
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  _/
  _/  MRNA ASN Generator functions 
  _/
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static BOOL isBegin (AC_TABLE tbl, char *txt)
{
  int ir, nn = strlen (txt) ;

  for (ir = 0 ; ir < tbl->rows ; ir++)
    if (!strncasecmp (txt, ac_table_printable (tbl, ir, 0, ""), nn))
      return TRUE ;
  return FALSE ;
}

static void ficheSeqSubmitHeader (vTXT blkp, AC_OBJ oGene, AC_OBJ oMrna, AC_TABLE gcDNA_clone, char style)
{
  int	ir, isComma = 0, isRand ;
  const char  *tNam = 0, *geneName = 0 ;
  AC_OBJ oTmp, oFound = 0 ;
  AC_TABLE gAuthors ;
  AC_HANDLE h = handleCreate () ;

  vtxtPrintf (blkp,  
	     "Seq-submit ::= {"
	     ) ;   
  
  /* -== == == == == == == == == == == == == == == == == = - /
   * -=  SUBMISSION INFO                 = - /
   * -== == == == == == == == == == == == == == == == == = - */      
  if (0) return ;

  vtxtPrintf (blkp,  
	     "sub {"
	     ) ;
  
  /* -== == == == == == == == == == == == == == == == == = - /
   * -=  CONTACT INFO                    = - /
   * -== == == == == == == == == == == == == == == == == = - */
  
  vtxtPrintf (blkp,  
	      "contact {"
	      "  contact {"
	      "    name name {"
	      "      last \"Thierry-Mieg\", "
	      "      first \"Danielle\" }, "
	      "    affil std {"
	      "      affil \"National Conter for Biotechnology Information, NIHZ\" , "
	      "      city \"Bethseda\" , "
	      "      sub \"MD\","
	      "      country \"USA\" , "
	      "      street \"Bdg 38A, 8600, Rockville Pike\" , "
	      "      email \"mieg@ncbi.nlm.nih.gov\" , "
	      "      postal-code \"20894\""
	      "    }"
	      "  } "
	      "}, "
	      ) ;
  
  /* -== == == == == == == == == == == == == == == == == = - /
   * -=  CITATION AND AUTHORS INFO       = - /
   * -== == == == == == == == == == == == == == == == == = - */
  
#define doCommap()		 if (isComma)vtxtPrintf (blkp, ", ") ;isComma = 1 ;
  
  vtxtPrintf (blkp,  
	     "cit {"
	     "  authors {"
	     "    names std {"
	     ) ;
  
  if (style == 's')
    {
      /**** get the author names */
      if ((oTmp = ac_tag_obj (oGene, "First_gb_author", h)))
	{
	  doCommap () ;
	  tNam = ac_tag_printable (oTmp, "Surname", "") ;
	  if (tNam)printfAsnAuthor (blkp, tNam, ac_tag_printable (oTmp, "Firstname", ""), 
				  ac_tag_printable (oTmp, "Initial", "")) ;
	  else printfAsnAuthor (blkp, ac_name (oTmp), "", "") ;
	}
      
      if (isBegin (gcDNA_clone, "y"))
	{ /* # y.kohara t.shini */
	  doCommap () ;printfAsnAuthor (blkp, "Kohara", "Yuji", "Y.") ;
	}
      if ((gAuthors = ac_tag_table (oGene, "Specialist", h)))
	{
	  AC_TABLE oFullName ;
	  for (ir = 0 ;ir <gAuthors->rows ;ir++)
	    {
	      oTmp = ac_table_obj (gAuthors, ir, 0, h) ;
	      oFullName =  ac_tag_table (oTmp,  "Full_name", h) ;
	      if (!oFullName || ! oFullName->rows  || oFullName->cols < 4 ||
		  ! ac_table_printable (oFullName, 0, 3, 0))
		continue ;
	
	      doCommap () ;
	      printfAsnAuthor (blkp, strnew (ac_table_printable (oFullName, 0, 1,""), h )
					     , strnew (ac_table_printable (oFullName, 0, 2, ""), h)
			       , strnew (ac_table_printable (oFullName, 0, 3, ""), h)) ;
	    }
	}
      if (isBegin (gcDNA_clone, "yk"))
	{ /* # y.kohara t.shini */
	  doCommap () ;printfAsnAuthor (blkp, "Shin-i", "Tadasu", "T.") ;
	}
      if (isBegin (gcDNA_clone, "mv") ||  /* # m.vidal */
	 isBegin (gcDNA_clone, "dvpl"))
	{ 
	  doCommap () ;printfAsnAuthor (blkp, "Vidal", "Marc", "M.") ;
	}
      if (isBegin (gcDNA_clone, "cm"))
	{ /*  # c martin r.waterston */
	  isRand = randint ()%2 ;
	  if (isRand == 0)
	    {doCommap () ;printfAsnAuthor (blkp, "Martin", "Chris", "C.") ;}
	  else {doCommap () ;printfAsnAuthor (blkp, "Waterston", "Robert", "R.H.") ;}
	}
      if (isBegin (gcDNA_clone, "yd"))
	{ /*  # D.L. Riddle and A.T. Puzyrev */
	  isRand = randint ()%2 ;
	  if (isRand == 0)
	    {doCommap () ;printfAsnAuthor (blkp, "Riddle", "Don", "D.L.") ;}
	  else {doCommap () ;printfAsnAuthor (blkp, "Puzyrev", "Anatoliy", "A.") ;}
	}
      if (isBegin (gcDNA_clone, "yc"))
	{ /*  # s.sugano y.suzuki */
	  isRand = randint ()%2 ;
	  if (isRand == 0)
	    {doCommap () ;printfAsnAuthor (blkp, "Suzuki", "Yutaka", "Y.") ;}
	  else {doCommap () ;printfAsnAuthor (blkp, "Sugano", "Sumio", "S.") ;}
	}
      if (0)if (isBegin (gcDNA_clone, "CE") || /* # public database */
	      isBegin (gcDNA_clone, "GB") || 
	      isBegin (gcDNA_clone, "EMBL"))
	{ 
	  doCommap () ;
	  if (oFound) printfAsnAuthor (blkp, ac_name (oFound), "Data", "") ;
	  /*printfAsnAuthor (blkp, "Genbank", "Data", "") ;*/
	}
      isRand = randint ()%5;
      if (isRand == 0)
	{doCommap () ;printfAsnAuthor (blkp, "Lowe", "Adam", "A.") ;}
      else if (isRand == 1)
	{doCommap () ;printfAsnAuthor (blkp, "Potdevin", "Michel", "M.") ; }
      else if (isRand == 2)
	{doCommap () ;printfAsnAuthor (blkp, "Simonyan", "Vahan", "V.") ; }
      else if (isRand == 3)
	{doCommap () ;printfAsnAuthor (blkp, "Thierry-Mieg", "Yann", "Y.") ;}
      else if (isRand == 4)
	{doCommap () ;printfAsnAuthor (blkp, "Sienkiewicz", "Marc", "M.") ;}
      doCommap () ;printfAsnAuthor (blkp, "Thierry-Mieg", "Danielle", "D.") ;
      doCommap () ;printfAsnAuthor (blkp, "Thierry-Mieg", "Jean", "J.") ;
    }
  
  /**** close brakets for CITATION AND AUTHORS INFO */
  vtxtPrintf (blkp, 
	      "    } , "
	      "affil std {"
	      "  affil \"National Institue of Genetics\" , "
	      "  div \"Genome biology laboratory, Center for genetic resource information\" ,"
	      "  city \"Mishima\" , "
	      "  country \"Japan\" , "
	      "  postal-code \"411-8540\" "
	      "    }" 
	      "  }" 
	      "}, "
	      ) ;
  
  /* -== == == == == == == == == == == == == == == == == = - /
   * -=  TOOLS INFO                      = - /
   * -== == == == == == == == == == == == == == == == == = - */
  geneName = ac_tag_printable(oGene, "NewName", ac_name(oGene)) ;
  vtxtPrintf (blkp,  
	     "subtype new , "
	     "tool \"AcemblySubmission_0.0\", "
	     "user-tag \"%s\""
	     , geneName) ;
  
  
  /**** close brakets for SUBMISSION INFO  */
  vtxtPrintf (blkp, "}, ") ;

  
  /* -== == == == == == == == == == == == == == == == == = - /
   * -=  DATA INFO                       = - /
   * -== == == == == == == == == == == == == == == == == = - */      
  
  vtxtPrintf (blkp, 
	     "data entrys {"
	     ) ;
  ac_free (h) ;
} /* ficheSeqSubmitHeader */

/******************************************************************/

static char * fAsnPostTreatCommentSection(char * src, AC_HANDLE h)
{
  int i,k,j;
  int isep,ipar=0,ilf=0;
  char ch,* dst;
#define addS(vsym)          {dst[k]=(vsym);k++;}
  
  dst = halloc (strlen(src)*2+1, h) ;
  
  isep=0;
  for(k=0,i=0;src[i];i++){
    ch=src[i];
    
    if(ch=='\r')continue;
    
    if(ch==' ' || ch=='\t'){ /* copy only one of the separators */
      if(isep==0){addS(' ');isep=1;}
      continue;
    }
    
    else if(ch==',' || ch==':' || ch==';' || ch=='.'){ 
      if(isep==1)k--;     /* remove blank before _commas, colons and dots */
      addS(ch);
      if(ch!='.' || src[i+1]<'0' || src[i+1]>'9'){addS(' ');isep=1;} /* if not a dot inside of the number */
      continue;
    }
    
    else if(ch=='\n'){
      for(j=1;src[i+j] && (src[i+j]==' ' || src[i+j]=='\t');j++); /* scan while blank */
      if(src[i+j]=='\n' || src[i+j]=='\r'){ /* if blank line separator is met */
	if(0 && ipar==0){addS('~');addS('\n');addS('~');addS('\n');ipar=1;}
	if(ipar==0){addS('\n');addS('\n');ipar=1;}
	i++;
      }
      else if(ilf==0 && ipar==0){addS('\n');ilf=1;}
      continue;
    }
    ipar=0;
    ilf=0;
    isep=0;
    addS(ch);
  }
  dst[k]=0; dst[k+1]=0;
  
  return dst;
}


/******************************************************************/
/* try to grab the NM from the mrna itself
   then try to grab from the gene if not used by another mrna
   or from the mtching genefinders
*/
  
static AC_OBJ fAsnMrna2NM (GMP *gmp, BOOL isCoding, AC_HANDLE hUser) 
{
  AC_OBJ oNM = 0 ;
  AC_TABLE gMrna ;
  int ir ;
  const char *ccp ;
  DICT *dict = 0 ;
  static DICT *geneDict = 0 ;
  AC_HANDLE h = handleCreate () ;

  if (!geneDict)
    geneDict = dictCreate (12) ;
  else if (!dictFind (geneDict, ac_name(gmp->gene), 0))
    { dictDestroy (geneDict) ;  geneDict = dictCreate (12) ; }
  dictAdd (geneDict, ac_name(gmp->gene), 0) ;
  
#define CHECKNM(_ccp) ((!dictFind (geneDict,_ccp,0)) && ((isCoding && !strncmp (_ccp, "NM_", 3))||(!isCoding && !strncmp (_ccp, "NR_", 3))))
  if (gmp->style != 'r')
    return 0 ;
  if (gmp->oNM)
    return gmp->oNM ; /* we do not need to copy, since the user only deallocates his handle */
  if (gmp->tg)
    {
      /* try to grab the NM from the mrna itself */
      if ((oNM = ac_tag_obj (gmp->mrna, "NM_id", h)) && CHECKNM(ac_name(oNM)))
	goto ok ;
      /* collect all NMs used by some other mrna of this gene */
      dict = dictCreate (12) ;
      for (ir = 1 ; ir <= dictMax (geneDict) ; ir++)
	dictAdd (dict, dictName (geneDict, ir), 0) ;
      if (gmp->tg && (gMrna = ac_tag_table (gmp->tg, "mrna", h)))
	for (ir = 0 ; ir < gMrna->rows ; ir++)
	  if ((ccp = ac_tag_printable (ac_table_obj (gMrna, ir, 0, h), "NM_id", 0)) && CHECKNM(ccp))
	    dictAdd (dict, ccp, 0) ;
      /* collect all NM known in this gene */
      if ((gMrna = ac_tag_table (gmp->gene, "NM_id", h)))
	for (ir = 0 ; ir < gMrna->rows ; ir++)
	  {
	    if ((oNM = ac_table_obj (gMrna, ir, 0, h)) && CHECKNM(ac_name(oNM)) &&
		dictAdd (dict, ac_name (oNM), 0))
	      break ;
	    else 
	      oNM = 0 ;
	  }
      /* steal the NM from the predicted genes of this gene */
      if (!oNM && (gMrna = ac_tag_table (gmp->gene, "genefinder", h)))
	for (ir = 0 ; ir < gMrna->rows ; ir++)
	  {
	    if ((oNM = ac_table_obj (gMrna, ir, 0, h)) && CHECKNM(ac_name(oNM)) &&
		dictAdd (dict, ac_name (oNM), 0)) 
	      break ;
	    else 
	      oNM = 0 ;
	  }
    }
  
  else if (gmp->pg)
    {
      /* try to grab the NM from the mrna itself */
      if ((oNM = ac_tag_obj (gmp->pg, "NM_id", h)) && CHECKNM(ac_name(oNM)) )
	goto ok ;
      /* collect all NMs used by some other pg of this gene */
      dict = dictCreate (12) ;
      for (ir = 1 ; ir <= dictMax (geneDict) ; ir++)
	dictAdd (dict, dictName (geneDict, ir), 0) ;
      if ((gMrna = ac_tag_table (gmp->gene, "Genefinder", h)))
	for (ir = 0 ; ir < gMrna->rows ; ir++)
	  if ((ccp = ac_tag_printable (ac_table_obj (gMrna, ir, 0, h), "NM_id", 0)) && CHECKNM(ccp))
	    dictAdd (dict, ccp, 0) ;
      /* collect all NM known in this gene */
      if ((gMrna = ac_tag_table (gmp->gene, "NM_id", h)))
	for (ir = 0 ; ir < gMrna->rows ; ir++)
	  {
	    if ((oNM = ac_table_obj (gMrna, ir, 0, h)) && CHECKNM(ac_name(oNM)) &&
		dictAdd (dict, ac_name (oNM), 0))
	      break ;
	    else
	      oNM = 0 ;
	  }
    }
 
 ok:     
  if (oNM && CHECKNM(ac_name(oNM)))
    { /* success, save the result in the database */
      dictAdd (geneDict, ac_name (oNM), 0) ;
      gmp->oNM = oNM ;

      if (gmp->tg)
	fprintf (stdout, "\n\nmRNA %s\nNM_Id %s\n\n", ac_name(gmp->mrna), ac_name (oNM)) ;
      else if (gmp->pg)
	fprintf (stdout, "Sequence %s\nNM_Id %s\n\n", ac_name(gmp->pg), ac_name (oNM)) ;
    }
  else
    oNM = 0 ; /* temporary non correct NM i.e. pg.ZK637.8a */
  dictDestroy (dict) ;
  oNM = ac_copy_obj (oNM, hUser) ;
  ac_free (h) ;

  return oNM ;
}  /* fAsnMrna2NM */

/******************************************************************/
/* 
 * called 3 times:
    sworm.c:ficheMrna-> (oMrna, 0, 1, style) (style == r | s)
    fichegraph:         (acOrigin(lObj), cpyBufr, 1, style);
    ficheasn.c:         (oTranscript, 0, 0, 'r') ;
*/
char *fAsnGenerateMRNA (vTXT blkp, GMP *gmp, char * ficheComments, int isHeader)
{
  vTXT buf ;
  AC_HANDLE h = handleCreate () ;
  int		ir, isComma = 0 ;
  int		Total_length, Length_5prime_UTR, Length_3prime_UTR ;
  AC_OBJ oNM = 0 ;
  AC_OBJ  oTranscribed_gene, oProduct = 0, oGene, oMrna ;
  AC_TABLE  oTmp, gPolyA_Signal, gSplicing, gcDNA_clone ;
  const char *ccp, *extName, *NewName = 0, *tNam = 0 ;
  char *ptr, *descrCompletness[4] = {"no-ends", "no-right", "no-left", "complete"} ;
  char geneName[1024] ;
  char PolyA_Signal[1024], mrnaName[1024], prRealName[1024] ;
  int  PolyA_variant = 0, PolyA_signal_pos = 0, completeness ;
  char *sDna, *sPeptide, *titleFromFiche, *commentSection, *ficheTreatedComments = 0 ;
 

#define doComma()		 if (isComma)vtxtPrintf (blkp, ", ") ;isComma = 1 ;

  ficheSectionizeFiche(ficheComments, &titleFromFiche, &commentSection);
  
  /* generate asn style comment section */	
  if(commentSection)
    ficheTreatedComments = fAsnPostTreatCommentSection (commentSection, h);
  
  buf = vtxtCreate () ;
  sprintf (mrnaName, "%s.m", ac_name (gmp->mrna)) ;
  extName = strrchr (mrnaName, '.') ;
  if (!extName)extName = mrnaName+strlen (mrnaName) ;
  
  /* -== == == == == == == == == == == == == == == == == = - /
     * -=  GET AND CHECK MUSTEXIST TAGS    = - /
     * -== == == == == == == == == == == == == == == == == = - */		
  oGene = gmp->gene ;
  oProduct = gmp->product ;
  oMrna = gmp->mrna ;

  if ((oTranscribed_gene = ac_tag_obj (gmp->mrna, "From_gene", h)) == 0 || 
      (gmp->Spc == WORM && (NewName = ac_tag_printable (oTranscribed_gene, "NewName", 0)) == 0) ||
      (sDna = ac_obj_dna (gmp->mrna, h)) == 0 )
    {
      ac_free (h) ;
      return 0 ;
    }

  strcpy (prRealName, ac_name (oProduct)) ;
  strcat (prRealName, ".p") ;

  oNM = fAsnMrna2NM (gmp, TRUE, h) ;  /* export here only NM_* identifiers */
  
  sPeptide = oProduct ? ac_obj_peptide (oProduct, h) : 0 ;

  gcDNA_clone = ac_tag_table (oMrna, "cDNA_clone", h) ;
  
  strcpy (geneName, ac_name (oTranscribed_gene)) ;

  Total_length = strlen (sDna) ;
  Length_5prime_UTR = ac_tag_int (oMrna, "Length_5prime_UTR", 0) ;
  Length_3prime_UTR = ac_tag_int (oMrna, "Length_3prime_UTR", 0) ;
  /* Longest_ORF = ac_tag_int (oMrna, "Longest_ORF", 0) ; */
  /*  geneLoc = ac_tag_printable (oTranscribed_gene, "Gene", 0) ; */
  
  gPolyA_Signal = ac_tag_table (oMrna, "PolyA_Signal", h) ;
  if (gPolyA_Signal)
    {
      strcpy (PolyA_Signal, ac_table_printable (gPolyA_Signal, 0, 0, "")) ;
      if (strstr (PolyA_Signal, "AATAA"))
	{
	  PolyA_signal_pos = ac_table_int (gPolyA_Signal, 0, 1, 0) ;
	  PolyA_variant = 0 ;
        }
      else if (strstr (PolyA_Signal, "Variant"))
	{
	  strcpy (PolyA_Signal, ac_table_printable (gPolyA_Signal, 0, 1, "")) ;
	  PolyA_signal_pos = ac_table_int (gPolyA_Signal, 0, 2, 0) ;
	  PolyA_variant = 1 ;
	}
      PolyA_signal_pos = Total_length-PolyA_signal_pos ;
    }
  else 
    PolyA_Signal[0] = 0 ;
  
  /* completness */
  completeness = 0 ;
  if (oProduct && ac_has_tag (oProduct, "Met") && ac_has_tag (oProduct, "NH2_Complete"))
    completeness |= 0x01 ;
  if (oProduct && ac_has_tag (oProduct, "COOH_Complete")) 
    completeness |= 0x02 ;
  if (oProduct && ac_has_tag (oProduct, "Complete")) 
    completeness |= 0x03 ;
  
  /* -== == == == == == == == == == == == == == == == == = - /
   * -=  START INFORMATION OUTPUT        = - /
   * -== == == == == == == == == == == == == == == == == = - */
  if (isHeader) 
    ficheSeqSubmitHeader (blkp, oGene, oMrna, gcDNA_clone, gmp->style) ;
  
  /* -== == == == == == == == == == == == == == == == == = - /
   * -=  SET INFO                        = - /
   * -== == == == == == == == == == == == == == == == == = - */  
  /* SET NUC Prot  */
  vtxtPrintf (blkp, 
	     "set {"
	     "  class nuc-prot, " 
	     ) ;

  /* -== == == == == == == == == == == == == == == == == = - /
   * -=  DESCR INFO                      = - /
   * -== == == == == == == == == == == == == == == == == = - */      
  { /* DESCR INFO */
    vtxtPrintf (blkp, 
	       "descr {"
	       "  create-date std "
	       ) ; 
    printfAsnDate (blkp) ;
    vtxtPrintf (blkp, 
	       ", " 
	       "name "
	       ) ;

    /* we no longer use titleFromFiche
     * we do NOT want to extract an edited title from the fiche
     * it must be edited as kantor->title
     */
    vtxtClear (buf) ;
    ptr = gtMrnaName (buf, gmp) ;

    vtxtPrintf (blkp, "\"%s\"" , ptr) ;

    /* get the title */
    vtxtClear (buf) ;
    ptr = gtMrnaTitle (buf, gmp) ;
    if (0 && /* kans does not like title */
	ptr) 
      { 
	vtxtPrintf (blkp, ", title \"") ;
	vtxtPrintf (blkp, "%s\"", ptr) ;
      }
    vtxtClear (buf) ;
  }
  
  {
    int nerr ;
    BOOL fM = ac_has_tag (oMrna, "From_AM") ;
	
    if (fM && !ac_has_tag (oMrna, "Tiling_error"))
      fM = 0 ;
    nerr = ac_tag_int (oMrna, "Tiling_error", 0) ;
    if (!nerr)
      fM = 0 ;
    vtxtClear (buf) ;

    if (gmp->style == 's' && ficheTreatedComments)
      ptr = ficheTreatedComments ;
    else
      ptr = ficheNewMrnaSubmissionComment (buf, gmp) ; 
    
    /* was until july 8, 2003 in 's' case : ptr = ficheMrnaContent (buf, gmp))) */
    fM = 0 ; /* july 18, 2003, Danielle no longer want the Note .. */
    if (fM || ptr)
      {
	vtxtPrintf (blkp, ", comment \"Summary:") ;
	if (ptr)
	  { 
	    vtxtPrintf (blkp, " %s", ptr) ;
	    
	    vtxtClear (buf) ;
	  }
	
	if (fM)
	  {
	    vtxtBreak (blkp) ;
	    vtxtPrintf (blkp, "Note that this sequence does not match the genome in %d position%s"
			, nerr
			, nerr > 1 ? "s" : "") ;
	  }
	vtxtPrintf (blkp, ".\"") ;
      }
    
    vtxtClear (buf) ;
  }
  /* -== == == == == == == == == == == == == == == == == = - /
   * -=  SOURCE INFO                     = - /
   * -== == == == == == == == == == == == == == == == == = - */  
  { /* SOURCE INFO   */
    AC_TABLE gExpression_profile ;
    
    /*char tPos[256] ;*/
    if (gmp->style == 'r')
      {
	int iPre = 0, iSto = 0, iAlt = 0 ;
	
	if ((gSplicing = ac_tag_table (oMrna, "Splicing", h)))
	  {
	    int iss ;
	    const char * txt ;

	    for (iss = 0 ;iss <gSplicing->rows ;iss++)
	      {
		txt = ac_table_printable (gSplicing, iss, 4, "") ;
		if (!strcasecmp (txt, "Alternative_Exon"))iAlt++ ;
		if (!strcasecmp (txt, "Predicted_Exon"))iPre++ ;
		else if (!strcasecmp (txt, "Stolen_Exon"))iSto++ ;
	      }
	  }
	vtxtPrintf (blkp, 
		   ", user {"
		   "  type str \"RefGeneTracking\" , "
		   "  data {"
		   "    {"
		   "      label str \"Status\" , "
		   "      data str \"%s\" "
		   "    }"
		   , ( iSto || iPre) ? "Provisional" : "Reviewed") ;

	{
	  AC_TABLE gTil, datab ;
	  AC_OBJ oSeq ;
	  int ir, jr, dummy, nn = 0 ;
	  const char *ptr, *ginum = 0, *acc = 0 ;
	  DICT *tDict = dictCreate (50) ;
	  
	  if ((gTil = ac_tag_table (oMrna, "Tiling_path", h)))
	    {
	      for (ir = 0 ;ir <gTil->rows ;ir++)
		{
		  oSeq = ac_table_obj (gTil, ir, 2, h) ;
		  ptr = ac_name (oSeq) ;
		  if ((datab = ac_tag_table (oSeq, "Database", h)))
		    nn++ ;
		}
	    }
	 
	  vtxtPrintf (blkp, 
		     " , "
		     "    {"
		     "      label str \"Assembly\" , "
		     "      data fields {"
		     ) ;
	  if (nn)
	    {
	      if ((gTil = ac_tag_table (oMrna, "Tiling_path", h)))
		{
		  for (jr = ir = 0 ;ir <gTil->rows ;ir++)
		    {
		      oSeq = ac_table_obj (gTil, ir, 2, h) ;
		      ptr = ac_name (oSeq) ;
		      if (!(datab = ac_tag_table (oSeq, "Database", h)))
			continue ;
		      if (!dictAdd (tDict, ptr, 0))
			continue ;
		      if (!strncmp (ptr, "EMBL:", 4))acc = ptr+5 ;
		      if (!strncmp (ptr, "GGB:", 4))ptr++ ;
		      if (!strncmp (ptr, "GB:", 3))acc = ptr+3 ;
		      if ((datab = ac_tag_table (oSeq, "Database", h)))
			{
			  if (!strcmp (ac_table_printable (datab, 0, 0, ""), "GI"))
			    ginum = ac_table_printable (datab, 0, 1, "") ;
			  if (ginum && !sscanf (ginum, "%d", &dummy))
			    ginum = 0 ;
			  if (!strcasecmp (ac_table_printable (datab, 0, 0, ""), "genbank"))
			    acc = ac_table_printable (datab, 0, 1, "") ;
			}
		      if (!acc) acc = ptr ;
		      if (jr++)vtxtPrintf (blkp, ", ") ;/* fAsnGenerateMRNA  */
		      vtxtPrintf (blkp, 
				  "        {"
				  "          label id %d , "
				  "          data fields {"
				  "            {"
				  "              label str \"accession\" , " 
				  "              data str \"%s\" "
				  "            } "
				  , ir, acc) ;
		      if (ginum)
			vtxtPrintf (blkp, 
				    ", "
				    "            {"
				    "              label str \"gi\" , "
				    "              data int %s "
				    "            }"
				    , ginum) ;
		      vtxtPrintf (blkp, 
				  ""
				  "          }" 
				  "        }"
				  ) ;
		      
		    }
		}
	    }
	  else
	    {
	      AC_OBJ oGF = ac_tag_obj (oGene, "Genefinder", h) ;
	      vtxtPrintf (blkp, 
			  "{"
			  "     label id 0 ,"
			  "     data fields {{"
			  "                  label str \"name\" ,"
			  "                  data  str \"WormBase CDS:%s\""
			  "          }}"
			  " } "
			  , oGF ? ac_name(oGF) : ac_name(oGene) 
			  ) ;	      
	    }
	  dictDestroy (tDict) ;	  
	}
	
	vtxtPrintf (blkp, 
		   "      }"
		   "    }"
		   "  }"
		   "}, "
		   ) ;
      }
    else vtxtPrintf (blkp, ",") ;
      
    
    
    /* may be but what about the clones ficheAsnDescrBioSource (blkp, oMap, "source") ; */	     
    
    vtxtPrintf (blkp, 
	       
	       "source {"
	       ) ;
    ficheDescrWormSource (blkp) ;
    
    vtxtPrintf (blkp, ",  subtype {") ;	 /* open subtypes */

    ficheAsnMapSubtype (blkp, gmp) ;    
    {
      AC_OBJ obj = 0 ;

      if (gcDNA_clone || ac_has_tag (oProduct, "Primers"))
	ficheAsnCloneSubtype (blkp, gmp, gmp->gene, gmp->mrna, 0,  ", ") ;
 
      if (
	  (
	   (gExpression_profile = ac_tag_table (oTranscribed_gene, "Expression_profile", h)) && 
	   gExpression_profile->cols > 3
	   ) ||
	  ac_has_tag (oGene, "Pattern") ||
	  ac_has_tag (oGene, "Expr_pattern") ||
	  (
	   (obj = gtMrna2Product (oMrna, 0)) &&
	   ac_has_tag (obj, "Psort")
	   )
	  )
	ficheAsnCloneLibSubtype (blkp, gmp, ", ") ;
      ac_free (obj) ;
    }
    vtxtPrintf (blkp, "  }") ; /* close subtypes */
    vtxtPrintf (blkp, "}") ; /* close source */
  
    {  /* KEYWORD */
      if (gmp->style == 's')
	vtxtPrintf (blkp, 
		    ", genbank { keywords { \"Worm Transcriptome Project\" }}"
		    ) ;

    }
    {  /* PUBLICATION INFO  */
      isComma = 0 ;
      if (gmp->style == 's')
	{
	  


	  vtxtPrintf (blkp, 
		     ", "
		     "pub { "
		      ) ;
	  vtxtPrintf (blkp, 
		     "  pub {"
		     "    gen {"
		     "      cit \"Unpublished\" , "
		     "      authors {"
		     "        names std {"
		     ) ;
	  
	  /**** get the author names */
	  doComma () ;printfAsnAuthor (blkp, "Kohara", "Yuji", "Y.") ;
	  doComma () ;printfAsnAuthor (blkp, "Shin-i", "Tadasu", "T.") ;
	  doComma () ;printfAsnAuthor (blkp, "Suzuki", "Yutaka", "Y.") ;
	  doComma () ;printfAsnAuthor (blkp, "Sugano", "Sumio", "S.") ; 
	  doComma () ;printfAsnAuthor (blkp, "Thierry-Mieg", "Danielle", "D.") ;  
	  doComma () ;printfAsnAuthor (blkp, "Thierry-Mieg", "Jean", "J.") ; 
	  
	  vtxtPrintf (blkp, 
		     "        }"
		     "      }, "
		     "      title \"The Caenorhabditis elegans transcriptome project, a complementary view of the genome.\""
		     "    }"
		     "  }"
		      ) ;
	  
	  vtxtPrintf (blkp, 
		     "}"
		     ) ;
	  ficheAsnGenePublications (blkp, gmp, oGene, ", pub ") ;
	}
      else if (1 || gmp->style == 'r')
	ficheAsnGenePublications (blkp, gmp, oGene, ", pub") ;
    }
    vtxtPrintf (blkp, "}, ") ; /* close descr */
  } /* END SOURCE-INFO */


  /* -== == == == == == == == == == == == == == == == == = - /
   * -=  SEQ-SET INFO  oMRNA                  = - /
   * -== == == == == == == == == == == == == == == == == = - */      

  {   /* SEQ-SET INFO  sequence Block */
    vtxtPrintf (blkp, 
			    "seq-set {"
			    ) ;
    /* -== == == == == == == == == == == == == == == == == = - /
     * -=  SEQ-MRNA INFO                   = - /
     * -== == == == == == == == == == == == == == == == == = - */  
    { /*  SEQ-MRNA */
      vtxtPrintf (blkp, "seq {" ) ;
      if (0) printf ("YYYY %s", mrnaName) ;
      { /*  SEQ-MRNA INFO */
	BOOL isMrnaComplete = ac_has_tag (oMrna, "Complete") ;
	BOOL isFound3p = ac_has_tag (oMrna, "Found3p") ;
	BOOL isFound5p = ac_has_tag (oMrna, "Found5p") ;
	char *completeness = isMrnaComplete ? "complete" : 
	  (isFound5p ? "has-left" :
	   (isFound3p ? "has-right" : "no-ends")) ;
	{
	  vtxtPrintf (blkp, "id {") ; /* open id */
	  vtxtPrintf (blkp, "%s " , ficheAsnId (mrnaName, 0)) ;
	  if (oNM) /* fAsnGenerateMRNA  */
	    vtxtPrintf (blkp, 
		       ", other { accession \"%s\" } "
		       , ac_name (oNM)
		       ) ;
	  vtxtPrintf (blkp, "}") ; /* close id */
	}
 
	
	{
	  vtxtPrintf (blkp, ", descr {") ; /* open descr */
	  vtxtPrintf (blkp, "    molinfo {" 
		     "      biomol mRNA, completeness %s "
		     "    }"
		     , completeness
		     ) ;
	     /* get the title */
	  ptr = 0 ;
	  if (0 && titleFromFiche && *titleFromFiche)
	    {
	      char *cp ;
	      vtxtClear (buf) ;
	      vtxtPrintf (buf, titleFromFiche) ;
	      ptr = vtxtPtr (buf) ;
	      cp = ptr + strlen (ptr) - 1 ;
	      while (cp >= ptr && (*cp == '\n' || *cp == '\r' || *cp == '\t'))
		*cp-- = 0 ;
	    }
	  if (!ptr || !*ptr)
	    { 
	      vtxtClear (buf) ;
	      ptr = gtMrnaTitle (buf, gmp) ;
	    }
	  if (1 && ptr) /* kans does not like title */
	    { 
	      vtxtPrintf (blkp, ", title \"") ;
	      vtxtPrintf (blkp, "%s\"", ptr) ;
	    }
	  vtxtClear (buf) ;
	  
	  vtxtPrintf (blkp, "}") ; /* close descr */
	}
	
	{  /* dna */
	  int nrm = 0 ;
	  vtxtPrintf (blkp, 
		     ", "
		     "inst {"
		     "  repr raw , " 
		     "  mol rna , " 
		     "  length %i , " 
		     "  seq-data" 
		     "  iupacna \""
		     , Total_length) ;
	  
	  vtextUpperCase (sDna) ;
	  vtxtSequence (blkp, sDna) ;
	  ptr = vtxtPtr (blkp) ; ptr += strlen (vtxtPtr (blkp)) - 1 ;
	  while (*ptr == '\n' || *ptr == ' ') { nrm++ ; ptr-- ; }
	  if (nrm) *ptr = '\"' ;
	  else vtxtPrintf (blkp, "\"") ;
	  vtxtPrintf (blkp, " }") ;
	}
	
	{ /* SEQ-MRNA ANNOTATION INFO */
	  BOOL isFirst = TRUE ;
	  vtxtPrintf (blkp, 
		     ", annot {"
		     "  {"
		     "    data "
		     ) ;
	 
	  { /* SEQ-MRNA ftable GENERAL INFO */
	    
	    vtxtPrintf (blkp, 
		       "ftable {"	
		       ) ;
	    
	    if (SHOW_GENE_ON_NM)
	      {
		ficheAsnGeneSeqFeat (blkp, gmp, 0, Total_length - 1, mrnaName, TRUE, FALSE) ;
		 isFirst = FALSE ; 
	      }
	    /* -== == == == == == == == == == == == == == == == == = - /
	     * -=  SEQ-MRNA 5prime UTR INFO        = - /
	     * -== == == == == == == == == == == == == == == == == = - */                      
	    if (Length_5prime_UTR)
	      {
		/*AC_OBJ   oTmp ;
		  char mBuf[128] ;*/
		if (!isFirst) { vtxtPrintf (blkp, ", ") ;} isFirst = FALSE ;
		
		vtxtPrintf (blkp, 
			   "{ "
			   "  data imp {" 
			   "    key \"5\'UTR\""
			   "  }"
			   ) ;
		vtxtPrintf (blkp, ",  comment \"") ;
		ficheMRNA5PrimeParagraphContent (blkp, gmp) ;
		vtxtPrintf (blkp, "\"") ;
		
		vtxtPrintf (blkp, 
			   "  , location int {"
			   "    from 0 , "  
			   "    to %i , " 
			   "    id %s" 
			   "  }, " 
			   "  exp-ev experimental" 
			   "}"
			   , Length_5prime_UTR-1, ficheAsnId (mrnaName, 0)) ;
	      }
	    /* -== == == == == == == == == == == == == == == == == = - /
	     * -=  SEQ-MRNA 3prime UTR INFO        = - /
	     * -== == == == == == == == == == == == == == == == == = - */                      
	    if (Length_3prime_UTR) 
	      {
		if (!isFirst) { vtxtPrintf (blkp, ", ") ;} isFirst = FALSE ;
		
		vtxtPrintf (blkp, 
			   "{"
			   "  data imp {"
			   "    key \"3\'UTR\""
			   "  } ") ;
		vtxtPrintf (blkp, ",  comment \"The 3\' UTR") ;
		ficheMRNA3PrimeParagraphContent (blkp, gmp) ;
		vtxtPrintf (blkp, "\"") ;
		
		vtxtPrintf (blkp, 
			   "  , location int {" 
			   "    from %i , " 
			   "    to %i , " 
			   "    strand plus , " 
			   "    id %s"
			   "  }, "
			   "  exp-ev experimental" 
			   "}" 
			   , Total_length-Length_3prime_UTR, Total_length-1, ficheAsnId (mrnaName, 0)) ;
		
	      }
	    
	    /* -== == == == == == == == == == == == == == == == == = - /
	     * -=  SEQ-UNSURE                      = - /
	     * -== == == == == == == == == == == == == == == == == = - */                      
	    if ((gSplicing = ac_tag_table (oMrna, "Splicing", h)))
	      {
		int isPred, iexon ;
		const char * pgNam = ac_tag_printable (oTranscribed_gene, "Matching_genefinder_gene", "") ;
		
		for (iexon = 0, ir = 0 ;ir <gSplicing->rows ;ir++)
		  {
		    tNam = ac_table_tag (gSplicing, ir, 4, "") ;
		    isPred = 0 ;
		    if (strstr (tNam, "Exon"))
		      {
			if (
			    gmp->style != 's' && 
			    (!strcasecmp (tNam, "Stolen_Exon") || ! (isPred = strcasecmp (tNam, "Predicted_Exon")))
			    )
			  {
			    if (!isFirst) { vtxtPrintf (blkp, ", ") ;} isFirst = FALSE ;
			    
			    vtxtPrintf (blkp, 
				       "{"
				       "  data imp {"
				       "    key \"unsure\""
				       "  }, "
				       "  comment \"exon %d stolen from %s %s\", "
				       "  location int {" 
				       "    from %i , " 
				       "    to %i , " 
				       "    strand plus , " 
				       "    id %s"
				       "  }"
				       "}" 
				       , iexon+1, !isPred ? "predicted gene" : "variant mRNA", !isPred ? pgNam : ""
					, ac_table_int (gSplicing, ir, 2, 0) - 1
					, ac_table_int (gSplicing, ir, 3, 0) - 1
					, ficheAsnId (mrnaName, 0)) ;
			  }
			iexon++ ;
		      }
		  }
	      }
	    /* -== == == == == == == == == == == == == == == == == = - /
	     * -=  SEQ-EXON                        = - /
	     * -== == == == == == == == == == == == == == == == == = - */                      
	    if ((gSplicing = ac_tag_table (oMrna, "Splicing", h)))
	      {
		int iexon, x1, x2 ;

		for (iexon = 0, ir = 0 ;ir <gSplicing->rows ;ir++)
		  {
		    tNam = ac_table_printable (gSplicing, ir, 4, "") ;
		    x1 = ac_table_int (gSplicing, ir, 2, 0) - 1 ;
		    x2 = ac_table_int (gSplicing, ir, 3, 0) - 1 ;
		    if (x1 < 0) x1 = 0 ;
		    if (x2 > Total_length - 1) x2 = Total_length - 1 ;
		    if (strstr (tNam, "Exon"))
		      {
			if (!strstr (tNam, "Stolen") && !strstr (tNam, "Predicted"))
			  {
			    if (!isFirst) { vtxtPrintf (blkp, ", ") ;} isFirst = FALSE ;
			    
			    vtxtPrintf (blkp, 
				       "{"
				       "  data imp {"
				       "    key \"exon\""
				       "  }, "
				       /* "  except TRUE, " */
				       "  comment \"%s %d length %i bp\", "
				       "  location int {" 
				       "    from %i , " 
				       "    to %i , " 
				       "    strand plus , " 
				       "    id %s"
				       "  }, "
				       /*											"  except-text \"defined by genomic alignment\""*/
				       "}" 
				       , tNam, iexon+1, x2 - x1 + 1
				       , x1, x2, ficheAsnId (mrnaName, 0)) ;
			  }
			iexon++ ;
		      }
		  }
	      }
	    /* -== == == == == == == == == == == == == == == == == = - /
	     * -=  SEQ-INTRON                      = - /
	     * -== == == == == == == == == == == == == == == == == = - */                      
	    if ((gSplicing = ac_tag_table (oMrna, "Splicing", h)))
	      {
		
		for (ir = 0 ;ir <gSplicing->rows ;ir++)
		  {
		    tNam = ac_table_tag (gSplicing, ir, 4, "") ;
		    if (0 && 
			!strcmp (ac_table_printable (gSplicing, ir, 5, ""), "gt_ag"))
		      continue ;
		    
		    if (strstr (tNam, "Intron"))
		      {
			if (!isFirst) { vtxtPrintf (blkp, ", ") ;} isFirst = FALSE ;

			vtxtPrintf (blkp, 
				   "{"
				   "  data imp {"
				   "    key \"misc_feature\""
				   "  }, "
				   "  comment \"%s length %d bp, type %s\", "
				   "  location pnt {"
				   "      point %i , " 
				   "      strand plus , " 
				   "      id %s ,"
			           "   fuzz lim tr"
				   "  }"
				   "%s"
				   "}" 
				    , strnew (ac_table_printable (gSplicing, ir, 4, ""), h)
				    , ac_table_int (gSplicing, ir, 7, 0) 
					      , strnew (ac_table_printable (gSplicing, ir, 5, ""), h)
				    , ac_table_int  (gSplicing, ir - 1, 3, 0) - 1
				    , ficheAsnId (mrnaName, 0)
				    , !strcmp (strnew (ac_table_printable (gSplicing, ir, 5, ""),h), "gt_ag")
				    || strstr (tNam, "Stolen") ?
				    "" : ",  exp-ev experimental") ;
		      }
		  }
	      }
	    /* -== == == == == == == == == == == == == == == == == = - /
	     * -=  SEQ-ANTISENS                    = - /
	     * -== == == == == == == == == == == == == == == == == = - */                      
	    if ((ficheTGAntisensParagraphContent (buf, gmp, TRUE)))
	      {
		if (!isFirst) { vtxtPrintf (blkp, ", ") ;} isFirst = FALSE ;
		
		vtxtCleanHtml (buf) ;
		vtxtPrintf (blkp, 
			    "{"
			    "  data imp {"
			    "    key \"misc_binding\""
			    "  }, "
			    "  comment \"%s\", "
			    "  location whole %s ,"
			    "  qual { "
			    "    { qual \"bound_moiety\" , val \"RNA\" }"
                            "       } "
			    "}"
			   , vtxtPtr (buf), ficheAsnId (mrnaName, 0)) ;
		vtxtClear (buf) ;
	      }

	    /* -== == == == == == == == == == == == == == == == == = - /
	     * -=  SEQ-MRNA PolyA INFO             = - /
	     * -== == == == == == == == == == == == == == == == == = - */                      
	    if (ac_has_tag (oMrna, "PolyA_found"))
	      { 
		if (!isFirst) { vtxtPrintf (blkp, ", ") ;} isFirst = FALSE ;
		vtxtPrintf (blkp, 
			   "{"
			   "  data imp {"
			   "    key \"polyA_site\""
			   "  }" 
			   
			   ) ;
		
		if ((oTmp = ac_tag_table (oMrna, "PolyA_found", h) ))
		  {
		    int nPolyA = 0, x1, x2 ;
		    
		    vtxtPrintf (blkp, ", comment \"") ;
		    for (x1 = ir = 0 ;ir <oTmp->rows ;ir++)
		      {
			x2 = ac_table_int (oTmp, ir, 0, 0) ;
			if (x2 > x1 + 10 || x2 < x1 - 10) nPolyA++ ;
			x1 = x2 ;
		      }
		    if (nPolyA > 1)
		      vtxtPrintf (blkp, "%d different polyA addition sites ", nPolyA) ;
		    else
		      vtxtPrintf (blkp, "PolyA ") ;
		    vtxtPrintf (blkp, "visible in %s", 
			       oTmp->rows > 20 ? " (20 shown)" : "") ;
		    for (ir = 0 ;ir <20 && ir <oTmp->rows ;ir++)
		      {
			vtxtPrintf (blkp, "%s%s", (ir ? ", ":""), gtYkName (ac_table_printable (oTmp, ir, 1, ""))) ;
		      }
		    vtxtPrintf (blkp, "\"") ;
		  }
		vtxtPrintf (blkp, 
			   "  , location pnt {" 
			   "    point %i , "
			   "    strand plus , "
			   "    id %s"
			   "  }, " 
			   "  exp-ev experimental"
			   "}"
			   , Total_length - 1, ficheAsnId (mrnaName, 0)) ;
	      }
	    /* -== == == == == == == == == == == == == == == == == = - /
	     * -=  SEQ-MRNA PolyAsignal UTR INFO   = - /
	     * -== == == == == == == == == == == == == == == == == = - */                      
	    {
	      
	      if (PolyA_Signal[0])
		{
		  if (!isFirst) { vtxtPrintf (blkp, ", ") ;} isFirst = FALSE ;
		  vtxtPrintf (blkp, 
			     "{" 
			     "  data imp {"
			     "    key \"polyA_signal\" " 
			     "  }, " 
			     "  comment \"%s %s\", "
			     "  location int {" 
			     "    from %i , " 
			     "    to %i , " 
			     "    strand plus , " 
			     "    id %s" 
			     "  }" 
			     "}"
			     , ( PolyA_variant ? "variant " : "standard " )
			     , PolyA_Signal
			     , PolyA_signal_pos, PolyA_signal_pos+5, ficheAsnId (mrnaName, 0)) ;
		}
	    }

	    /**** close brakets for SET-MRNA ANNOTATION INFO  */
	    vtxtPrintf (blkp, 
		       "    }"
		       "  }"
		       "}"
		       ) ;
	  }
	}
      }
      
      /**** close brakets for SEQ-MRNA INFO  */
      vtxtPrintf (blkp, "}") ;
    }
    /* -== == == == == == == == == == == == == == == == == = - /
       * -=  SEQ-PROTEIN INFO                = - /
       * -== == == == == == == == == == == == == == == == == = - */  
    if (oProduct) 
      {  /* SEQ-PROTEIN INFO  */
	vtxtPrintf (blkp, ", seq {" ) ;
	
	{ /* id */
	  vtxtPrintf (blkp, "  id { ") ;
	  vtxtPrintf (blkp, "%s",ficheAsnId (prRealName, 0)) ;
	  if (oNM &&
	      (ccp = ac_tag_printable (oNM, "NPtext", 0)))
	    vtxtPrintf (blkp, 
		       ", other { accession \"%s\" } "
		       , ccp
		       ) ;/* fAsnGenerateMRNA  */
	  vtxtPrintf (blkp, "}") ; /* close id */

	  vtxtPrintf (blkp, 
		     "  , descr {"  
		     "    molinfo {" 
		     "    biomol peptide , "  
		     "    tech concept-trans , "  
		     "    completeness %s } "
		      /* no title here, this contradicts sequin expectations */
		     "  }"
		     , descrCompletness[completeness]) ;
	}
	
	{ /* AA */
	  char *cp = sPeptide ;
	  if (*cp == 'L' && (completeness & 0x01))
	    *cp = 'M' ;

	  vtxtPrintf (blkp, 
		     ", inst {" 
		     "    repr raw , " 
		     "    mol aa , " 
		     "    length %i , "       
		     "    seq-data" 
		     "    ncbieaa \""
		     , strlen (sPeptide)) ;
	  vtxtSequence (blkp, sPeptide) ;
	  vtxtPrintf (blkp, "\""
		     "}, ") ;
	}
   
	{ /* prot  annots */
	  vtxtPrintf (blkp, 
		     "annot {" 
		     "  {" 
		     "    data ftable {" 
		     ) ;
	  if (SHOW_GENE_ON_NP)
	    {
	      ficheAsnGeneSeqFeat (blkp, gmp, 0, strlen (sPeptide) - 1, prRealName, FALSE, FALSE) ;
	      vtxtPrintf (blkp, " ,") ;
	    }
	  vtxtPrintf (blkp, "{ data prot ") ;
	  ficheAsnProductRef (blkp, gmp, oProduct, TRUE) ;
	  if (! ac_has_tag (oProduct, "Complete"))
	    vtxtPrintf (blkp, ", partial TRUE") ;
	    vtxtPrintf (blkp, 
			"      , location int {" 
			"          from %i , " 
			"          to %i , " 
			"          id %s" 
			, 0
			, strlen (sPeptide) - 1
			, ficheAsnId (prRealName, 0)) ;
	  if (! ac_has_tag (oProduct, "NH2_complete"))
	    vtxtPrintf (blkp, 
			"       , fuzz-from lim lt" ) ;
	  if (! ac_has_tag (oProduct, "COOH_complete"))
	    vtxtPrintf (blkp, 
			"       , fuzz-to lim gt" ) ;
	  vtxtPrintf (blkp,
		      "        }" ) ;
	  vtxtPrintf (blkp,
		      "      }" ) ;

	  ficheAsnPFAM (blkp, oProduct, prRealName, sPeptide) ;
	  vtxtPrintf (blkp, 
		     ""
		     "    }" 
		     "  }"
		     "}"
		     ) ;
	}
	
	/**** close brakets for SEQ-PROTEIN INFO  */
	vtxtPrintf (blkp, "}") ;
      }
    /**** close brakets for SEQ-SET INFO  */
    vtxtPrintf (blkp, "}, ") ;
  }
 
  if (oProduct)
    { /* mrna seq-set  ANNOTATION INFO */
      int p1 = -9999, p2 = -9999 ; /* position of the product in mrna coordinates */
      AC_TABLE gProd = ac_tag_table (gmp->mrna, "Product", h) ;
      
      for (ir = 0 ; p1 == -9999 && ir < gProd->rows ; ir++)
	if (!strcmp (ac_name(oProduct), ac_table_printable (gProd, ir, 0, "")))
	  { p1 = ac_table_int (gProd, ir, 1, -9999); p2 = ac_table_int (gProd, ir, 2, -9999); }
      if (p1 != -9999)
	{
	  BOOL isKozak = FALSE ; 
	  if (ac_has_tag (oProduct, "NH2_Complete") &&
	      ac_has_tag (oProduct, "at_position_1") &&
	      ac_has_tag (oProduct, "First_Kozak") &&
	      ac_tag_int (oProduct, "First_ATG", 1) > 1)
	    isKozak = TRUE ;
	  vtxtPrintf (blkp, "annot {" ) ;
	  {
	    vtxtPrint (blkp, 
		       "  {"
		       "    data ftable {" 
		       "      {" 
		       "        data cdregion {" 
		       "          frame one, " 
		       "          code { id 1 }"
		       ) ;
	    if (isKozak)
	      vtxtPrintf (blkp, 
			  "       , code-break {"
			  "           {"
			  "             loc int {"
			  "               from %d, to %d,"
			  "		      strand plus,"
			  "               id %s"
			  "               },"
			  "             aa ncbieaa 77"
			  "           }"
			  "         }"
			  , p1 - 1, p1 + 1,  ficheAsnId (mrnaName, 0)) ;
	    vtxtPrintf (blkp, 
			"        }, " 
			"        %s"
			"        %s"
			"        product whole %s , " 
			, (completeness >= 0 && completeness  < 3) ? "  partial TRUE , " : ""
			, isKozak ? ", except TRUE," : ""
			, ficheAsnId (prRealName, 0)
			) ;
	    vtxtPrintf (blkp, 
			"        location int {" 
			"          from %i , " 
			"          to %i , " 
			"          id %s" 
			, p1 - 1
			, p2 - 1 
			, ficheAsnId (mrnaName, 0)) ;
	    if (! (completeness & 0x1))
	      vtxtPrintf (blkp, 
			  "       , fuzz-from lim lt" ) ;
	    if (! (completeness & 0x2))
	      vtxtPrintf (blkp, 
			  "       , fuzz-to lim gt" ) ;
	    vtxtPrint (blkp,
		       "        }" 
			) ;
	    if (isKozak)
	      vtxtPrintf (blkp,  
			  " ,except-text \"alternative start codon\""
			  ) ;

	  }
	  
	  ficheMrnaDbXref (blkp, gmp, oGene, oProduct, ", dbxref ") ;
	  
	  vtxtPrintf (blkp, 
		      "      }" 
		      "    }" 
		      "  }" 
		      "}" 
		      ) ;
	}
    } /* mrna seq-set annot info */
  
  /**** close brakets for DATA ->SET INFO  */
  vtxtPrintf (blkp, "}") ;

  /**** close brakets for DATA INFO  */
  if (isHeader)vtxtPrintf (blkp, "}") ;

  /**** close brakets for INFO OUTPUT */
  if (isHeader)vtxtPrintf (blkp, "}") ;


  if (vtxtPtr (blkp)) 
    vtxtPrintf (blkp, "\n\n") ;

  vtxtDestroy (buf) ;

  ac_free (h) ;

  return vtxtPtr (blkp) ;
} /* fAsnGenerateMRNA */


/*
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  _/
  _/  TG ASN Generator functions 
  _/
  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/***************************************************************************/
/***************************************************************************/
/*  GENE_ref tool-box*/

/***********************/
/*  GENE_ref dumper */
/* The ASN schema is
		

 Gene-ref ::= SEQUENCE {
 locus VisibleString OPTIONAL ,        -- Official gene symbol
 allele VisibleString OPTIONAL ,       -- Official allele designation
 desc VisibleString OPTIONAL ,         -- descriptive name
 maploc VisibleString OPTIONAL ,       -- descriptive map location
 pseudo BOOLEAN DEFAULT FALSE ,        -- pseudogene
 db SET OF Dbtag OPTIONAL ,            -- ids in other dbases
 syn SET OF VisibleString OPTIONAL ,   -- synonyms for locus
 locus-tag VisibleString OPTIONAL }    -- systematic gene name (e.g., MI0001, ORF0069)
*/

/***********************/

static BOOL ficheIsPseudogene (AC_OBJ oGene)
{
  BOOL ok = FALSE ;
  AC_OBJ oGF= 0 ;

  if (!ac_has_tag (oGene, "Transcribed_gene") &&
      (oGF = ac_tag_obj (oGene, "Genefinder", 0)) &&
      (strstr (ac_tag_printable (oGF, "Method", ""), "seudogene")))
    ok = TRUE ;

  ac_free (oGF) ;
  return ok ;
} /* ficheIsPseudogene */

static BOOL ficheAsnGeneRef (vTXT blkp, GMP *gmp, BOOL xref)
{
  AC_TABLE gMap ;
  DICT *dict = 0 ;
  int ir ;
  const char *ccp ;
  float mpos ;
  AC_OBJ oGene = gmp->gene ;
  AC_HANDLE h = handleCreate () ;
  vTXT buf ;   buf = vtxtCreate () ;

  vtxtPrintf (blkp, 
	    "{"
	    "  locus \"%s\" "
	    , ac_name (oGene)) ;
  

  if (0 && /* superflu */
      (ccp = gtGeneDescriptor (buf, gmp)))
    {
      vtxtPrintf (blkp, ", "
	      "  desc "
		) ;
      vtxtPrintf (blkp, "\"%s\" ", ccp) ;
    }
  if (1 && /* tatiana asks me to put the gene title here, sept29 2003 */
      (ccp = gtGeneTitle (buf, gmp, TRUE)))
    {
      vtxtPrintf (blkp, ", "
	      "  desc "
		) ;
      vtxtPrintf (blkp, "\"%s\" ", ccp) ;
    }
  
  if (xref)
    {
      vTXT vtxtMap = vtxtCreate () ;
      int nmap = 0 ;
      char *ptr ;

      if ((gMap = ac_tag_table (oGene, "InterpolatedMap", h)) && 
	  (mpos = ac_table_float (gMap, 0, 1, -9999)) != -9999)
	{
	  vtxtPrintf (vtxtMap
		       , "%s%s;%s%3.2f (interpolated genetic position)"
		      , nmap++ ? "; " : ""
		      , ac_table_printable (gMap, 0, 0, ""), mpos > 0 ? "+" : "", mpos) ;
	}
      if ((gMap = ac_tag_table (oGene, "Map", h)) && 
	 (mpos = ac_table_float (gMap, 0, 2, -9999)) != -9999)
	{
	  vtxtPrintf (vtxtMap
		      , "%s%s;%s%3.2f (measured genetic position)" 
		      , nmap++ ? "; " : ""
		      , ac_table_printable (gMap, 0, 0, ""), mpos > 0 ? "+" : "", mpos) ;
	}

      if ((ptr = vtxtPtr (vtxtMap)))
	vtxtPrintf (blkp, ",  maploc \"%s\" ", ptr) ; 

      vtxtDestroy (vtxtMap) ;
    }

  /* pseudogenes */
  if (ficheIsPseudogene (oGene))
    vtxtPrintf (blkp, ", pseudo TRUE") ;
  /* fill the dic with  all the aliases of the gene */

  if (gmp->style == 'r') /* not in submissions */
    {
      dict = gtGeneAliases (gmp,  TRUE) ;
      if (dictMax (dict) > 1)
	{
	  vtxtPrintf (blkp, ", ") ;
	  vtxtPrintf (blkp, "  syn {") ;
	  for (ir = 1 ;ir <= dictMax (dict) ;ir++)
	    {
	      vtxtPrintf (blkp, 
			  "    \"%s\"%s", dictName (dict, ir), ir == dictMax (dict) ? "" :", ") ;
	    }
	  vtxtPrintf (blkp, "  }" ) ;
	}
    }
  
  
  if (1 && /* i suppress that because it comes in an ugly way in the gb flat file as the first note */
      (ccp = ac_tag_printable (oGene, "NewName", 0)))
    vtxtPrintf (blkp, 
	      ", "
	      "  locus-tag \"%s\" "
	      , ccp) ;
  
  vtxtPrintf (blkp, "}") ;
  dictDestroy (dict) ;
  vtxtDestroy (buf) ;
  ac_free (h) ;

  return TRUE ;
} /* ficheAsnGeneRef  */

static char *ficheAsnPubSet (vTXT blkp, AC_OBJ oGene, char *tag)
{ 
  AC_OBJ oPap ;
  AC_TABLE tbl ;
  int nn, ir, nPmid ;
  const char *ptr ;
  AC_HANDLE h = handleCreate () ;

  vtxtClear (blkp) ;
  if ((tbl = ac_tag_table (oGene, tag, h)))
    {
      
      for (nn = ir = 0 ; ir < tbl->rows ; ir++)
	{
	  oPap = ac_table_obj (tbl, ir, 0, h) ;
	  if (ac_has_tag (oPap, "PMID"))
	    nn++ ;
	}
      if (nn)
	{
	  for (nn = ir = 0 ; ir < tbl->rows ; ir++)
	    {
	      oPap = ac_table_obj (tbl, ir, 0, h) ;
	      nPmid = 0 ;
	      if ((ptr = ac_tag_printable (oPap, "PMID", 0)) &&
		  sscanf (ptr, "%d", &nPmid) &&
		  nPmid)
		vtxtPrintf (blkp, "%s pmid %d", 
			  nn++ ? ", ":"", nPmid) ;
	    }
	}
    }
  ac_free (h) ;

  return vtxtPtr (blkp) ;
} /* ficheAsnPubSet */

/* -== == == == == == == == == == == == == == == == == = - /
* -=  Genefinder dump    = - /
* -== == == == == == == == == == == == == == == == == = - */		
static void fichePerGeneFinderTrnaDump (vTXT  blkp, const char *chromName, GMP *gmp, int *isAny)
{
  char gfName[KANS_LIMIT] ;
  AC_HANDLE h = handleCreate () ;
  AC_OBJ oGene = gmp->gene, oGF = gmp->pg ;
  char* sDna = ac_obj_dna (oGF, h) ;
  int Total_length ;
  BOOL pseudo = FALSE ;
  char *ptr, *type = 0 ;
  const char *ccp ;
  vTXT buf ; buf = vtxtCreate () ;
  
  if (!sDna) return ;
  Total_length = strlen (sDna) ;
  sprintf (gfName, "%s.r", gtGfName (oGF)) ;
  
  
  if (isAny) vtxtPrintf (blkp, ", ") ;
  { /* SEQ class predicted tRNA */
    vtxtPrintf (blkp, 
		"seq {"
		) ;
    if (0) printf ("XXXX %s", gfName) ;
    
    {
      AC_OBJ oNM ;

      vtxtPrintf (blkp, "id {") ; /* open id */
      vtxtPrintf (blkp, "%s " , ficheAsnId (gfName, 0)) ;
      if ((oNM = fAsnMrna2NM (gmp, FALSE, h))) /* export here only NR_* identifiers */
	vtxtPrintf (blkp, 
		    ", other { accession \"%s\" } "
		    , ac_name (oNM)
		    ) ;
      vtxtPrintf (blkp, "}") ; /* close id */
    }
    
    if ((ptr = gtPredictedTrnaTitle (buf, gmp, 0, &type, h))) /* get the type */
      vtxtClear (buf) ;
      

    {
      char *type2 ; BOOL poor ;
      vtxtPrintf (blkp, ", descr {") ; /* open descr */ 
      vtxtPrintf (blkp, 
		  "  create-date std "
		  ) ;
      /* 
	 ACEDB                   ASN

	 Unprocessed_mRNA 	 premsg
	 CDS                     mRNA

         tRNA                    tRNA
         rRNA                    rRNA
         snRNA                   snRNA
         scRNA                   scRNA
                                 snoRNA
	 stRNA			 scRNA
         miRNA                   scRNA
         misc_RNA                transcribed-RNA
         Pseudogene              transcribed-RNA
         Transposon              transcribed-RNA
      */
      type2 = type ; poor = FALSE ; pseudo = FALSE ;
      if (!strcasecmp(type, "RNA")) { pseudo = TRUE ; poor = TRUE ;type2 = "transcribed-RNA" ; }
      if (!strcasecmp(type, "rRNA")) type2 = "rRNA" ;
      if (!strcasecmp(type, "stRNA")) type2 = "scRNA" ;
      if (!strcasecmp(type, "miRNA")) type2 = "scRNA" ;
      if (!strcasecmp(type, "misc_RNA"))  { poor = TRUE ; type2 = "transcribed-RNA" ; }
      if (!strcasecmp(type, "Pseudogene")) { pseudo = TRUE ; poor = TRUE ; type2 = "transcribed-RNA" ; }
      if (!strcasecmp(type, "Transposon")) {  poor = TRUE ; type2 = "transcribed-RNA" ; }

      printfAsnDate (blkp) ;
      vtxtPrintf (blkp, ",    molinfo {" 
		  "      biomol %s "
		  "    }"
		  , type2
		  ) ;
      
      if (!strcmp(type, "tRNA") && (ccp = ac_tag_printable (oGF, "Method", 0)))
	vtxtPrintf (blkp, ", comment \"%s annotated by WormBase using %s\"", type, ccp) ; 
   

      vtxtPrintf (blkp, 
		  ", user"
		  "{"
		  "  type str \"RefGeneTracking\" , "
		  "  data {"
		  "         {"
		  "           label str \"Status\" , "
		  "           data str \"%s\" "
		  "         }," 
		  "         { "
		  "           label str \"Assembly\" ,"
		  "           data fields {{"
		  "                          label id 0 ,"
		  "                          data fields {{"
		  "                                         label str \"name\" ,"
		  "                                         data  str \"WormBase %s:%s\""
		  "                                       }}"
		  "                       }}"
                  "         }"
		  "       } "
		  " } "
		  , poor ? "Provisional" : "Reviewed", type, ac_name(oGF) 
		  ) ;
      
      if (1 && /* kans does not like title */
	  (ptr = gtPredictedTrnaTitle (buf, gmp, TRUE, &type, h)))
	{ 
	  vtxtPrintf (blkp, ", title \"") ;
	  vtxtPrintf (blkp, "%s\"", ptr) ;
	}
      vtxtClear (buf) ;
      {
	vtxtPrintf (blkp, ", source {" ) ;	 /* open source */
	ficheDescrWormSource (blkp) ;
	
	{
	  AC_OBJ oPprod = 0 ;

	  vtxtPrintf (blkp, ",  subtype {") ;	 /* open subtypes */

	 
	  ficheAsnMapSubtype (blkp, gmp) ; 
	  
	  if (ac_has_tag (oGF, "Primers"))
	    ficheAsnCloneSubtype (blkp, gmp, oGene, 0, oGF,  ", ") ;

	  if (
	      ac_has_tag (oGene, "Pattern") ||
	      ac_has_tag (oGene, "Expr_pattern") ||
	      (
	       (oPprod = gtPredictedMrna2Product (oGF, 0)) &&
	       ac_has_tag (oPprod, "Psort")
	       )
	      )
	    ficheAsnCloneLibSubtype (blkp, gmp, ", ") ;
	  ac_free (oPprod) ;
	  vtxtPrintf (blkp, "  }") ; /* close subtypes */
	}
	vtxtPrintf (blkp, " }") ; /* close source */
      }
      
      ficheAsnGenePublications (blkp, gmp, oGene, ", pub") ;
      
      vtxtPrintf (blkp, "}") ; /* close descr */
    }
    
    { /* dna */
      vtxtPrintf (blkp, 
		  ", "
		  "inst {"
		  "  repr raw , " 
		  "  mol rna , " 
		  "  length %i , " 
		  "  seq-data" 
		  "  iupacna \""
		  , Total_length
		  ) ;
      
      vtextUpperCase (sDna) ;
      vtxtSequence (blkp, sDna) ;
      vtxtPrintf (blkp, "\" }") ;
    }
    
    if (SHOW_GENE_ON_NM)
      { /*  predicted trna seq  ANNOTATION INFO */
	vtxtPrintf (blkp, ", annot {" ) ; /* open annot */
	{
	  vtxtPrintf (blkp, 
		      "  {"
		      "    data ftable {" 
		      ) ;
	  ficheAsnGeneSeqFeat (blkp, gmp, 0, Total_length - 1, gfName, FALSE, pseudo) ;
	  vtxtPrintf (blkp, 
		      "    }" 
		      "  }"
		      ) ;
	}
      
      vtxtPrintf (blkp, "}") ;  /* close annot */
    } /* predicted trna seq annot info */
    
    vtxtPrintf (blkp, "}") ;  /* close seq predicted tRNA */
  } 
  *isAny = 1 ;
  
  ac_free (h) ;
  vtxtDestroy (buf) ;
} /* fichePerGeneFinderTrnaDump */

/* -== == == == == == == == == == == == == == == == == = - /
* -=  Genefinder dump    = - /
* -== == == == == == == == == == == == == == == == == = - */		

static void fichePerGeneFinderNucProtDump (vTXT  blkp, const char *chromName, GMP *gmp, int *isAny)
{
  char protName[KANS_LIMIT], gfName[KANS_LIMIT] ;
  AC_HANDLE h = handleCreate () ;
  AC_OBJ oGene = gmp->gene ;
  AC_OBJ oGF = gmp->pg ;
  char *sDna = ac_obj_dna (oGF, h) ;
  char *sPeptide = ac_obj_peptide (oGF, h) ;
  AC_OBJ oPredictedMrna =  ac_tag_obj (oGF, "Predicted_mRNA", h) ;
  AC_OBJ oPredictedProduct =  oPredictedMrna ?  gtMrna2Product (oPredictedMrna, h) : 0 ;
  AC_OBJ oNM = fAsnMrna2NM (gmp, TRUE, h) ; /* export here only NM_* identifiers */
  int Total_length ;
  const char *ccp ;
  char *ptr, *descrCompletness[4] = {"no-ends", "no-right", "no-left", "complete"} ;
  int completeness = 3 ;
  vTXT buf ; buf = vtxtCreate () ;
  
  Total_length = strlen (sDna) ;
  sprintf (protName, "%s.p", gtGfName (oGF)) ; 
  sprintf (gfName, "%s.m", gtGfName (oGF)) ;
  
  completeness = 3 ;
  if (ac_has_tag (oGF, "Start_not_found"))
    completeness &= ~(0x2) ;
  if (ac_has_tag (oGF, "End_not_found"))
    completeness &= ~(0x1) ;
  
  if (completeness == 3 && !ac_has_tag (oGF, "has_exact_est"))
    completeness = -1 ;
  
  if (isAny) vtxtPrintf (blkp, ", ") ;
  { /* SET class predicted NUC Prot */
    vtxtPrintf (blkp, 
		"set {"
		"  class nuc-prot, " 
		) ;
    
    { /* DESCR INFO */
      char *supportStatus[5] = { "Inferred", "Provisional",   "Provisional",   "Provisional", "Validated" };
      int supportLevel = 0, ncl = 0 ;

      vtxtPrintf (blkp, 
		  "descr {"
		  "  create-date std "
		  ) ; 
      printfAsnDate (blkp) ;
      if (0) /* tatiana */
	{
	  vtxtPrintf (blkp, 
		      ", " 
		      "name "
		      ) ;
	  vtxtPrintf (blkp, "\"%s\"" 
		      , gtPredictedMrnaName (buf, gmp, oGF)
		      ) ;
	  /* get the title */
	  if (0 && /* kans does not like title */
	      (ptr = gtPredictedMrnaTitle (buf, gmp)))
	    { 
	      vtxtPrintf (blkp, ", title \"") ;
	      vtxtPrintf (blkp, "%s\"", ptr) ;
	    }
	  vtxtClear (buf) ;
	}
      supportLevel = gtPredictedMrnaSupport (oGF, &ncl) ;
      vtxtPrintf (blkp, ", comment \"") ;

      vtxtClear (buf) ;
      if (TRUE && /* show the summary */
	  (ptr = ficheNewMrnaSubmissionComment (buf, gmp)))
	{
	  vtxtPrint (blkp, ptr) ;
	  vtxtClear (buf) ;
	}
      else
	switch (supportLevel)
	  {
	    /* REVIEWED REFSEQ: This record has been curated by NCBI staff. The
	       reference sequence was derived from totoid.
	    */
	  case 4:
	    /* 
	       VALIDATED REFSEQ: This record has undergone preliminary review of
	       the sequence, but has not yet been subject to final NCBI review.
	       The reference sequence was derived from totoid.
	    */
	    /* validated */
	    vtxtPrintf (blkp, "Fully supported by %d cDNA clones%s"
			, ncl > 0 ? ncl : -ncl
			, ncl > 0 ? "" : " mainly from the Worm Transcriptome Project"
			) ;
	    break ;
	  case 3:  
	    /* 
	       PROVISIONAL REFSEQ: This record has not yet been subject to final
	       NCBI review. The reference sequence was derived from
	       WormBase:totoid.    */
	    vtxtPrintf (blkp, "Partial support by %d cDNA clones%s"
			, ncl > 0 ? ncl : -ncl
			, ncl > 0 ? "" : " mainly from the Worm Transcriptome Project"
			) ;
	    break ;
	  case 2:                     /* provisional */
	    vtxtPrintf (blkp, "Partial support by the Vidal ORFeome project") ;
	    break ;
	  case 1:                     /* provisional */
	    vtxtPrintf (blkp, "Supported by expression data") ;
	    break ;
	  case 0:	    /*  MODEL REFSEQ: This record is predicted by genome sequence analysis 
				and is not yet supported by experimental evidence. The reference
				sequence was derived from totoid. */
	    vtxtPrintf (blkp, ", comment \"no cDNA, RT-PCR or expression data") ;
	    break ;
	  }  
      vtxtPrintf (blkp, ".\" ") ;

      vtxtPrintf (blkp, 
		  ", user"
		  "{"
		  "  type str \"RefGeneTracking\" , "
		  "  data {"
		  "         {"
		  "           label str \"Status\" , "
		  "           data str \"%s\" "
		  "         }," 
		  "         { "
		  "           label str \"Assembly\" ,"
		  "           data fields {{"
		  "                          label id 0 ,"
		  "                          data fields {{"
		  "                                         label str \"name\" ,"
		  "                                         data  str \"WormBase CDS:%s\""
		  "                                       }}"
		  "                       }}"
                  "         }"
		  "       } "
		  " } "
		  , supportStatus[supportLevel], ac_name(oGF) 
		  ) ;
      
      {
	vtxtPrintf (blkp, ", source {" ) ;	 /* open source */
	ficheDescrWormSource (blkp) ;

	{
	  AC_OBJ oPprod = 0 ;
	  vtxtPrintf (blkp, ",  subtype {") ;	 /* open subtypes */
	  
	  ficheAsnMapSubtype (blkp, gmp) ; 

	  if (ac_has_tag (oGF, "Primers"))
	    ficheAsnCloneSubtype (blkp, gmp, 0, 0, oGF,  ", ") ;

	  if (
	      ac_has_tag (oGene, "Pattern") ||
	      ac_has_tag (oGene, "Expr_pattern") ||
	      (
	       (oPprod = gtPredictedMrna2Product (oGF, 0)) &&
	       ac_has_tag (oPprod, "Psort")
	       )
	      )
	    ficheAsnCloneLibSubtype (blkp, gmp, ", ") ;
	  ac_free (oPprod) ;
	  vtxtPrintf (blkp, "  }") ; /* close subtypes */
	}
	vtxtPrintf (blkp, " }") ; /* close source */
      }
      
      ficheAsnGenePublications (blkp, gmp, oGene, ", pub") ;
      vtxtPrintf (blkp, "}, ") ; /* close descr */
    }
    
    /* -== == == == == == == == == == == == == == == == == = - /
     * -=  SEQ-SET INFO   oGF                 = - /
     * -== == == == == == == == == == == == == == == == == = - */      
    
    {   /*  predicted SEQ-SET INFO  sequence Block */
      vtxtPrintf (blkp, 
		 "seq-set {"
		 ) ;
      /* -== == == == == == == == == == == == == == == == == = - /
       * -=  SEQ-MRNA INFO                   = - /
       * -== == == == == == == == == == == == == == == == == = - */  
      { /*  predicted SEQ-MRNA */
	vtxtPrintf (blkp, "seq {" ) ;
	if (0) printf ("XXXX %s", gfName) ;

	{
	  vtxtPrintf (blkp, "id {") ; /* open id */
	  vtxtPrintf (blkp, "%s " , ficheAsnId (gfName, 0)) ;
	  if (oNM)
	    vtxtPrintf (blkp, 
			", other { accession \"%s\" } "
			, ac_name (oNM)
			) ;
	  vtxtPrintf (blkp, "}") ; /* close id */
	}
	
	{
	  BOOL isMrnaComplete = ac_has_tag (oGF, "has_exact_est") ;
	  isMrnaComplete = TRUE ; /* KANS, i cannot make that work */
	  vtxtPrintf (blkp, ", descr {") ; /* open descr */
	  vtxtPrintf (blkp, "    molinfo {" 
		     "      biomol mRNA, completeness %s "
		     "    }"
		     , isMrnaComplete ? "complete" : "partial"
		     ) ;
	  
	  if (1 && /* kans does not like title */
	      (ptr = gtPredictedMrnaTitle (buf, gmp)))
	    { 
	      vtxtPrintf (blkp, ", title \"") ;
	      vtxtPrintf (blkp, "%s\"", ptr) ;
	    }
	  vtxtClear (buf) ;
	  
	  vtxtPrintf (blkp, "}") ; /* close descr */
	}
	
	{ /* dna */
	  vtxtPrintf (blkp, 
		     ", "
		     "inst {"
		     "  repr raw , " 
		     "  mol rna , " 
		     "  length %i , " 
		     "  seq-data" 
		     "  iupacna \""
		     , Total_length
		     ) ;
	  vtextUpperCase (sDna) ;
	  vtxtSequence (blkp, sDna) ;
	  vtxtPrintf (blkp, "\" }") ;
	}
	
	
	if (SHOW_GENE_ON_NM)
	  { /*  predicted mrna seq  ANNOTATION INFO */
	    vtxtPrintf (blkp, ", annot {" ) ; /* open annot */
	    {
	      vtxtPrintf (blkp, 
			  "  {"
			  "    data ftable {" 
			  ) ;
	      ficheAsnGeneSeqFeat (blkp, gmp, 0, Total_length - 1, gfName, FALSE, FALSE) ;
	      vtxtPrintf (blkp, 
			  "    }" 
			  "  }"
			  ) ;
	    }
	    
	    vtxtPrintf (blkp, "}") ;  /* close annot */
	  } /* predicted mrna seq annot info */

	vtxtPrintf (blkp, "}") ;
      } /**** SEQ-MRNA INFO  */
      
      {  /*  predicted SEQ-PROTEIN INFO  */
	vtxtPrintf (blkp, ", seq {" ) ;

	{ /* id */
	  vtxtPrintf (blkp, "  id { ") ;
	  vtxtPrintf (blkp, "%s", ficheAsnId (protName, 0)) ;
	  
	  if (gmp->style == 'r')
	    {
	      if (oNM &&
		  (ccp = ac_tag_printable (oNM, "NPtext", 0)))
		vtxtPrintf (blkp, 
			    ", other { accession \"%s\" } "
			    , ccp
			    ) ;
	    }
	  vtxtPrintf (blkp, "}") ; /* close id */
	  vtxtPrintf (blkp, 
		      "  , descr {"  
		      "      molinfo {" 
		      "               biomol peptide , "  
		      "               tech concept-trans"  
		     ) ;
	  if (completeness >= 0 && completeness  <= 3) /* worth reporting */
	    vtxtPrintf (blkp, 
			"               , completeness %s "
			, descrCompletness[completeness]
		     ) ;
	  vtxtPrintf (blkp, 
		      "               }" 
		      /* no title here bad for sequin */
		      "           }"
		     ) ;
	}
	
	{ /* predicted AA */	
	  char *cp = sPeptide ;
	  if (*cp == 'L') *cp = 'M' ;

	  vtxtPrintf (blkp, 
		     ", inst {" 
		     "  repr raw , " 
		     "  mol aa , " 
		     "  length %i , "       
		     "  seq-data" 
		     "  ncbieaa \""
		     , strlen (sPeptide)) ;
	  vtxtSequence (blkp, sPeptide) ;
	  vtxtPrintf (blkp, "\""
		     "}, ") ;
	}
	
	{ /* predicted prot  annots */
	  vtxtPrintf (blkp, 
		     "annot {" 
		     "  {" 
		     "    data ftable {" 
		     ) ;
	  vtxtPrintf (blkp, "{ data prot ") ;
	  ficheAsnProductRef (blkp, gmp, oGF, FALSE) ;
	  vtxtPrintf (blkp, 
		     "        , location whole %s" 
		     "      }" 			
		     , ficheAsnId (protName, 0)) ;
	  if (SHOW_GENE_ON_NP)
	    {
	      vtxtPrintf (blkp, ", ") ;
	      ficheAsnGeneSeqFeat (blkp, gmp, 0, strlen (sPeptide) - 1, protName, FALSE, FALSE) ;
	    }
	  {
	    AC_OBJ pMrna, pProduct ;
	    pMrna = ac_tag_obj (oGF, "Predicted_mrna", h) ;
	    if (pMrna && (pProduct =  gtMrna2Product (pMrna, h)))
	      ficheAsnPFAM (blkp, pProduct, protName, sPeptide) ;
	  }

	  vtxtPrintf (blkp, "} } }") ;
	}
	/**** close brakets for SEQ-PROTEIN INFO  */
	vtxtPrintf (blkp, "}") ;
      }

    /**** close brakets for SEQ-SET INFO  */
      vtxtPrintf (blkp, "}, ") ;
    }
  
    if (1)
      { /*  predicted mrna seq-set  ANNOTATION INFO */
	vtxtPrintf (blkp, "annot {" ) ; /* open annot */
	{
	  vtxtPrintf (blkp, 
		      "  {"
		      "    data ftable {" 
		      "      {" 
		      "        data cdregion {" 
		      "          frame one, " 
		      "          code { id 1 }"
		      "        }, " 
		      "        %s"
		      "        product whole %s , " 
		      , (completeness >= 0 && completeness  < 3) ? "  partial TRUE , " : ""
		      , ficheAsnId (protName, 0)
		      ) ;
	  vtxtPrintf (blkp, 
		      "        location int {" 
		      "          from %i , " 
		      "          to %i , " 
		      "          id %s" 
		      "        }" 
		      , 0, Total_length - 1 , ficheAsnId (gfName, 0)) ;
	  oPredictedMrna =  ac_tag_obj (oGF, "Predicted_mRNA", h) ;
	  oPredictedProduct =   gtMrna2Product (oPredictedMrna, h) ;

	  ficheMrnaDbXref (blkp, gmp, oGene, oPredictedProduct, ", dbxref ") ;
	  vtxtPrintf (blkp, 
		     "      }" 
		     "    }" 
		     "  }"
		     ) ;
	}
	
	vtxtPrintf (blkp, "}") ;  /* close annot */
      } /* predicted mrna seq-set annot info */
 

    vtxtPrintf (blkp, "}") ; /* close set nuc prot*/
  }
  *isAny = 1 ;
  vtxtDestroy (buf) ;
  ac_free (h) ;

  return ;
} /* fichePerGeneFinderNucProtDump */

/* -== == == == == == == == == == == == == == == == == = - /
   * -=  Prodcut DUMP = - /
   * -== == == == == == == == == == == == == == == == == = - */  
/*
  Prot-ref ::= SEQUENCE {
    name SET OF VisibleString OPTIONAL ,      -- protein name
    desc VisibleString OPTIONAL ,      -- description (instead of name)
    ec SET OF VisibleString OPTIONAL , -- E.C. number (s)
    activity SET OF VisibleString OPTIONAL ,  -- activities


*/

static char *ficheAsnProductRef (vTXT blkp, GMP *gmp, AC_OBJ oProduct, BOOL isReal)
{
  char *ptr ;
  BOOL doQuote = TRUE ;
  vTXT buf ; buf = vtxtCreate () ;

  vtxtPrintf (blkp, 
	     "     {" 
	     "          name { " 
	     "            \""
	     ) ;
  vtxtPrintf (blkp, 
	     "%s\""
	     , isReal ? gtProductName (buf, gmp, oProduct) : gtPredictedProductName (buf, gmp, oProduct, TRUE)
	     ) ;
  
  
  vtxtPrintf (blkp, "}") ;
  
  if (0 && /* redundant with the summary */
      (ptr = gtGeneFunctionalDescriptor (buf, gmp, oProduct, doQuote)))
    {
      vtxtPrintf (blkp, ", activity ") ;
      vtxtPrintf (blkp, "{ %s }", ptr) ;
    }
  vtxtPrintf (blkp, 
	     "        }"
	     ) ;
  
  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
}

/* -== == == == == == == == == == == == == == == == == = - /
 * -=  MRNA DUMP = - /
 * -== == == == == == == == == == == == == == == == == = - */  
/*
  RNA-ref ::= SEQUENCE {
  type ENUMERATED {            -- type of RNA feature
  unknown (0) , 
  premsg (1) , 
  mRNA (2) , 
  tRNA (3) , 
  rRNA (4) , 
  snRNA (5) , 
  scRNA (6) , 
  snoRNA (7) , 
  other (255) } , 
  pseudo BOOLEAN OPTIONAL ,  
  ext CHOICE {
  name VisibleString ,        -- for naming "other" type
   tRNA Trna-ext } OPTIONAL }  -- for tRNAs
*/

static BOOL fichaAsnTrnaType (AC_OBJ pg, int *aatypep, int *codontypep)
{
  int i, nn ;
  char buf[4] ;
  const char *text = ac_tag_printable (pg, "tRNA", 0) ;

  if (text && strlen(text) > 4 && text[3] == ' ')
    {
      for (i = nn = 0 ; i < 3 ; i++, text++)
	{
	  nn *= 4 ;
	  switch (*text)
	    { /* Jonathan Kans' coding */
	    case 'T': nn += 0 ; buf[i] = T_ ; break ;
	    case 'C': nn += 1 ; buf[i] = C_ ; break ;
	    case 'A': nn += 2 ; buf[i] = A_ ; break ;
	    case 'G': nn += 3 ; buf[i] = G_ ; break ;
	    default:
	      return FALSE ;
	    }
	}
      buf[i] = 0 ;
    }
  else
    return FALSE ;
  
#ifndef AC_TEST
  {
    /* in the june 2003 release, i had aatype and codon type switched, and nn was 4 times too large ! */
    char *pn = pepShortName[(int) codon(buf)] ;
    if (strstr (text, pn)) /* wormbase used the codon */
      { /* success */
	*aatypep = (int) codon(buf); 
	*codontypep = nn ;
	return TRUE ;
      }
    pn = pepShortName[(int) antiCodon(buf)] ;
    if (strstr (text, pn)) /* wormbase used the anticodon */
      { /* success */
	*aatypep = (int) antiCodon(buf);
	*codontypep = nn ;
	return TRUE ;
      }
  }
#endif
  /* success */  
  return FALSE ;
}


static char *ficheAsnRnaRef (vTXT blkp, GMP *gmp, BOOL isPredicted, char *type)
{
  BOOL pseudo, poor ;
  int aatype, codontype ;
  char *type2 = 0 ;
  vTXT buf ; buf = vtxtCreate () ;

  if (type)
    {
#ifdef DO_NOT_DELETE_THIS_JUNK
      /* in principle, to do the same as in the nuc-prot case we shoud use as above */
      type2 = type ; poor = FALSE ; pseudo = FALSE ;
      if (!strcasecmp(type, "RNA")) { pseudo = TRUE ; poor = TRUE ;type2 = "transcribed-RNA" ; }
      if (!strcasecmp(type, "stRNA")) type2 = "scRNA" ;
      if (!strcasecmp(type, "miRNA")) type2 = "scRNA" ;
      if (!strcasecmp(type, "misc_RNA"))  { poor = TRUE ; type2 = "transcribed-RNA" ; }
      if (!strcasecmp(type, "Pseudogene")) { pseudo = TRUE ; poor = TRUE ; type2 = "transcribed-RNA" ; }
      if (!strcasecmp(type, "Transposon")) {  poor = TRUE ; type2 = "transcribed-RNA" ; }
      type2 = type ; poor = FALSE ; pseudo = FALSE ;
#else
      /* but in real asn life, the schema here is more limited (less predefined types, so we say */
      if (!strcasecmp(type, "RNA")) { pseudo = TRUE ; poor = TRUE ;type2 = "other" ; }
      if (!strcasecmp(type, "tRNA")) type2 = "tRNA" ;
      if (!strcasecmp(type, "rRNA")) type2 = "rRNA" ;
      if (!strcasecmp(type, "snRNA")) type2 = "snRNA" ;
      if (!strcasecmp(type, "stRNA")) type2 = "scRNA" ;
      if (!strcasecmp(type, "miRNA")) type2 = "scRNA" ;
      if (!strcasecmp(type, "misc_RNA"))  { poor = TRUE ; type2 = "other" ; }
      if (!strcasecmp(type, "Pseudogene")) { pseudo = TRUE ; poor = TRUE ; type2 = "other" ; }
      if (!strcasecmp(type, "Transposon")) {  poor = TRUE ; type2 = "other" ; }
#endif
    }
  vtxtPrintf (blkp, 
	    "{"
	    "    type %s"
	      , type2 ? type2 : "mRNA"
	    ) ;
  if (type &&
      !strcasecmp(type, "tRNA") &&
      gmp->pg &&
      fichaAsnTrnaType (gmp->pg, &aatype, &codontype))
    {
/*
 
  ext CHOICE {
        name VisibleString ,        -- for naming "other" type
        tRNA Trna-ext } OPTIONAL }  -- for tRNAs

  example:
                   ftable {
                     {
                       data
                         rna {
                           type tRNA ,
                           ext
                             tRNA {
                               aa
                                 ncbieaa 77 ,
                               codon {
                                 35 } } } ,
                       location


so i convert
ATG  = 16 * A + 4 T + G = 32 + 0 + 3 = 35
and so on for all dna triplets
using your TCAG 0123 ordering


and i convert

M ->   (int)'M' == 77

and so on using any of the one-char classic amino acid code

*/

  vtxtPrintf (blkp, 
	    "    , ext tRNA { aa ncbieaa %d, codon { %d } }"
	      , aatype, codontype
	    ) ;
    }
#if 0
  /* this else does not validate */
  else
    {
/*
	ispredicted	type	printed
	T		T	gtPredictedTrnaTitle (buf, gmp, FALSE, &type)
	T		F	gtPredictedMrnaName (buf, gmp->pg)
	F		X	gtMrnaName (buf, gmp)
*/
      vtxtPrintf (blkp,
		  "    , ext name \"%s\""
		  , isPredicted ? 
		  type ? gtPredictedTrnaTitle (buf, gmp, FALSE, &type, h) : gtPredictedMrnaName (buf, gmp->pg) 
		  : gtMrnaName (buf, gmp)) ; 
    }
#endif
  if (pseudo || poor) pseudo = poor = FALSE ; /* for compiler happiness */
  vtxtPrintf (blkp, "}") ;
  
  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
}

/* -== == == == == == == == == == == == == == == == == = - /
* -=  location    = - /
* -== == == == == == == == == == == == == == == == == = - */		
/* type 1: predicted_gene, 2: mRNA, 3: product*/
static char* ficheAsnMrnaLocation (vTXT  blkp, GMP *gmp, int type)
{
  int ir, isAny, x1, x2, a1, p1, p2, c1, c2, u1, u2 ; 
  AC_TABLE oChrom = 0 ;
  const char *oChromName = 0 ;
  int start = 0 ;
  int end = 0 ;
  BOOL isDown, is5pPartial = FALSE, is3pPartial = FALSE, is5p, is3p ;
  AC_OBJ oTmp = 0 ;
  AC_TABLE gSpli = 0, gProd = 0 ;
  AC_HANDLE h = handleCreate () ;

  p1 = 888888 ; p2 = -888888 ;
  switch (type)
    {
    case 1: /* predicted_gene */
      gSpli = ac_tag_table (oTmp = gmp->pg, "Source_Exons", h) ;
      break ;
    case 2: /* mRNA */
      if (! ac_has_tag (gmp->mrna, "Found5p"))
	is5pPartial = TRUE ;
      if (! ac_has_tag (gmp->mrna, "Found3p"))
	is3pPartial = TRUE ;
      gSpli = ac_tag_table (oTmp = gmp->mrna, "Splicing", h) ;
      break ;
    case 3: /* product */
      if (! ac_has_tag (gmp->product, "NH2_complete"))
	is5pPartial = TRUE ;
      if (! ac_has_tag (gmp->product, "COOH_complete"))
	is3pPartial = TRUE ;
      gSpli = ac_tag_table (oTmp = gmp->mrna, "Splicing", h) ;
      gProd = ac_tag_table (gmp->mrna, "Product", h) ;
      for (ir = 0 ; gProd && ir < gProd->rows ; ir++)
	if (!strcmp (ac_name (gmp->product), ac_table_printable (gProd, ir, 0, "")))
	  {
	    p1 = ac_table_int (gProd, ir, 1, 999999) ;
	    p2 = ac_table_int (gProd, ir, 2, -999999) ;
	    break ;
	  }
      break ;
    default:
      messcrash ("Unknown type %d in ", type) ;
    }      
  oChrom = ac_tag_table (oTmp, "IntMap", h) ;
  oChromName = ac_table_printable (oChrom, 0, 0, "") ;
  start = ac_table_int (oChrom, 0, 1, 0) - 1 ;
  end = ac_table_int (oChrom, 0, 2, 0) - 1 ;
  isDown = start < end ? TRUE : FALSE ;

  vtxtClear (blkp) ;
  if (gSpli)
    for (isAny = 0, ir = 0 ;ir < gSpli->rows ; ir++)
      {
	a1 = c1 = ac_table_int (gSpli, ir, 0, 0) ;
	c2 = ac_table_int (gSpli, ir, 1, 0) ;
	is5p = is3p = FALSE ;

	switch (type)
	  {
	  case 1: /* pg */
	    break ;
	  case 2: /* mrna */  
	    if (!strstr (ac_table_printable (gSpli, ir, 4, ""), "xon"))
	      continue ;
	    if (ir == 0 && is5pPartial) is5p = TRUE ;
	    if (ir == gSpli->rows - 1 && is3pPartial) is3p = TRUE ;
	    break ;
	  case 3: /* product */
	    if (!strstr (ac_table_printable (gSpli, ir, 4, ""), "xon"))
	      continue ;
	    x1 = ac_table_int (gSpli, ir, 2, 0) ;
	    x2 = ac_table_int (gSpli, ir, 3, 0) ;
	    c1 = x1 < p1 ? p1 : x1 ; 
	    c2 = x2 < p2 ? x2 : p2 ; 
	    if (c1 > c2)
	      continue ;
	    if (c1 == p1 && is5pPartial) is5p = TRUE ;
	    if (c2 == p2 && is3pPartial) is3p = TRUE ;
	    c1 = a1 + c1 - x1 ; c2 = a1 + c2 - x1 ;
	    break ;
	  }
	
	if (isDown)
	  {
	    u1 = start + c1 - 1 ;
	    u2 = start + c2 - 1 ;
	  } 
	else 
	  {
	    u2 = start - c1 + 1 ;
	    u1 = start - c2 + 1 ;
	  }	
	vtxtPrintf (blkp, 
		    "%s"
		    "    int {"
		    "      from %d , "
		    "      to %d , "
		    "      strand %s , "
		    "      id %s %s %s"
		    "    }"
		    , isAny++?", ":""
		    , u1, u2
		    , isDown ? "plus" : "minus"
		    , ficheAsnId (oChromName, 0)
		    , ((isDown && is5p) || (!isDown && is3p)) ? ", fuzz-from lim lt" : ""
		    , ((isDown && is3p) || (!isDown && is5p)) ? ", fuzz-to lim gt" : ""
		  ) ;
      }
  ac_free (h) ;
  return vtxtPtr (blkp) ;  
}  /* ficheAsnMrnaLocation */


   /*
   * We export here the geneBox accoring to the following schema
   * then directly called from this function
   * we export the mRNAs and the predictedGene as successive Seq-feat
   *
   Seq-feat ::= SEQUENCE {
    data SeqFeatData ,           -- the specific data  // i.e. mrna
    partial BOOLEAN OPTIONAL ,    -- incomplete in some way?
    except BOOLEAN OPTIONAL ,     -- something funny about this?
    comment VisibleString OPTIONAL , 
    product Seq-loc OPTIONAL ,    -- product of process
    location Seq-loc ,            -- feature made from
    qual SEQUENCE OF Gb-qual OPTIONAL ,  -- qualifiers
    title VisibleString OPTIONAL ,   -- for user defined label
    ext User-object OPTIONAL ,    -- user defined structure extension
    cit Pub-set OPTIONAL ,        -- citations for this feature
    exp-ev ENUMERATED {           -- evidence for existence of feature
        experimental (1) ,        -- any reasonable experimental check
        not-experimental (2) } OPTIONAL , -- similarity, pattern, etc
    xref SET OF SeqFeatXref OPTIONAL ,   -- cite other relevant features
    dbxref SET OF Dbtag OPTIONAL ,  -- support for xref to other databases
    pseudo BOOLEAN OPTIONAL ,     -- annotated on pseudogene?
    except-text VisibleString OPTIONAL } -- explain if except = TRUE
   }
*/
static char *ficheAsnProductGenomeLocation (vTXT blkp, GMP *gmp)
{
  char *ptr ; 
  BOOL isKozak = FALSE ;
  vTXT buf ;   buf = vtxtCreate () ;
  
  vtxtClear (blkp) ;  
  vtxtPrintf (blkp, "{") ; /* open ficheAsnProductGenomeLocation */
  /* Export [Product-ref] */
  /*
    vtxtPrintf (blkp, "data prot ") ;
    ficheAsnProductRef (blkp, oGene, oProduct) ;
  */
  vtxtPrintf (blkp, 
	      "data cdregion { frame 1, code { id 1 }}" 
	      ) ;
  if (ac_has_tag (gmp->product, "NH2_Complete") &&
      ac_has_tag (gmp->product, "First_Kozak") &&
      ac_has_tag (gmp->product, "at_position_1") &&
      ac_tag_int (gmp->product, "First_ATG", 1) > 1)
    isKozak = TRUE ;
 /* partial */
  if (!ac_has_tag (gmp->product, "NH2_complete") || !ac_has_tag (gmp->product, "COOH_Complete"))
    vtxtPrintf (blkp, ", partial TRUE") ;
  if (gmp->useAm || isKozak)
    vtxtPrintf (blkp, ", except TRUE ") ;
   
  /* except BOOLEAN OPTIONAL ,     -- something funny about this? */

  /* product of process */
  vtxtPrintf (blkp, ",  product  whole %s" , ficheAsnId (ac_name (gmp->product), ".p")) ;
  
  /* export location */	
  if ((ptr = ficheAsnMrnaLocation (buf, gmp, 3)))
    vtxtPrintf (blkp, 								   
	      ", location mix {"
	      "  %s"
	      "}"
	      , ptr) ;
  
  /* title */
  if (0 &&
      (ptr = gtProductTitle (buf, gmp)))
    {
       vtxtPrintf (blkp, ", title ") ;
       vtxtPrintf (blkp, "\"%s\"", ptr) ;
    }
 
  /* citations */
  if (0 &&   /* tatiana */
      (ptr = ficheAsnPubSet (buf, gmp->gene, "Reference")))
    {
      vtxtPrintf (blkp, ", cit pub {") ;
      vtxtPrintf (blkp, ptr) ;
      vtxtPrintf (blkp, "}") ;
    }
  /* experimental, always true in class mrna */
  vtxtPrint (blkp, ", exp-ev experimental") ;
  if (isKozak) 
    vtxtPrint (blkp, ", except-text \"alternative start codon\"") ;
  /* the mrna dbxref */
  ficheMrnaDbXref (blkp, gmp, gmp->gene, gmp->product, ", dbxref ") ;
 
  if (gmp->useAm)
    vtxtPrintf (blkp, ", except-text \"unclassified translation discrepancy\" ") ;
  vtxtPrintf (blkp, "}") ; /* close ficheAsnProductGenomeLocation */

  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
}  /* ficheAsnProductGenomeLocation */

/*********************************************************************/

static char *ficheAsnMrnaGenomeLocation (vTXT blkp, GMP *gmp)
{
  char *ptr ; 
  vTXT buf ;   buf = vtxtCreate () ;
  
  vtxtClear (blkp) ;
  /* Export [Rna-ref] */
  vtxtPrintf (blkp, "{  data rna ") ;
  ficheAsnRnaRef (blkp, gmp, 0, 0) ;
  
  /* partial */
  if (!ac_has_tag (gmp->mrna, "Complete"))
    vtxtPrintf (blkp, ", partial TRUE") ;
  if (gmp->useAm)
    vtxtPrintf (blkp, ", except TRUE ") ;

  /* product of process */
  vtxtPrintf (blkp, ",  product  whole %s" , ficheAsnId (ac_name (gmp->mrna), ".m")) ;
  
  
  /* export location */	
  if ((ptr = ficheAsnMrnaLocation (buf, gmp, 2)))
    vtxtPrintf (blkp, 								   
	      ", location mix {"
	      "  %s"
	      "}"
	      , ptr) ;
  
  /* title */
  if (1 && /* kans does not like title */
      (ptr = gtMrnaTitle (buf, gmp)))
    {
       vtxtPrintf (blkp, ", title ") ;
       vtxtPrintf (blkp, "\"%s\"", ptr) ;
    }
 
  /* citations */
  if (0 && /* tatiana */
      (ptr = ficheAsnPubSet (buf, gmp->gene, "Reference")))
    {
      vtxtPrintf (blkp, ", cit pub {") ;
      vtxtPrintf (blkp, ptr) ;
      vtxtPrintf (blkp, "}") ;
    }
  /* experimental, always true in class mrna */
  vtxtPrintf (blkp, ", exp-ev experimental") ;
  
  /* the mrna dbxref */

  ficheMrnaDbXref (blkp, gmp, gmp->gene, 0, ", dbxref ") ;
  if (gmp->useAm)
    vtxtPrintf (blkp, ", except-text \"unclassified transcription discrepancy\" ") ;

  vtxtPrintf (blkp, "}") ;

  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
}  /* ficheAsnMrnaGenomeLocation */

/*********************************************************************/
/*********************************************************************/

static char *ficheAsnPredictedMrnaGenomeLocation (vTXT blkp, GMP *gmp)
{
  char *type=0, *ptr ;
  vTXT buf ;   buf = vtxtCreate () ;
  
  vtxtClear (blkp) ;
  
  if (!ac_has_tag(gmp->pg,"CDS"))
    {
      gtPredictedTrnaTitle (blkp, gmp, 0, &type, 0) ;
      vtxtClear (blkp) ;
    }
  /* Export [Rna-ref] */
  vtxtPrintf (blkp, "{  data rna ") ;
  ficheAsnRnaRef (blkp, gmp, 1, type) ;
  ac_free (type) ;
  /* partial , default=TRUE do not export  */
  if ((0 && ac_has_tag (gmp->pg, "Start_not_found")) || 
      ac_has_tag (gmp->pg, "End_not_found"))
    vtxtPrintf (blkp, ", partial TRUE") ;

  /* vtxtPrintf (blkp, ", partial TRUE") ;  since this is just a CDS */

  /* comment */
  if (0) vtxtPrintf (blkp, 								   
	    " , comment \"ZZZZZZ THIS IS A predicted COMMENT ZZZZZZZZZZZZZZ\""
	    ) ;
  /* product of process */
  vtxtPrintf (blkp, ",  product  whole %s" , ficheAsnId ( gtGfName (gmp->pg), ac_has_tag(gmp->pg,"CDS") ? ".m" : ".r")) ;
  
    /* export location */	
  if ((ptr = ficheAsnMrnaLocation (buf, gmp, 1)))
    vtxtPrintf (blkp, 								   
	      ", location mix {"
	      "  %s"
	      "}"
	      , ptr) ;
  
  /* title */
  if (1 && /* kans does not like title */
      (ptr = gtPredictedMrnaTitle (buf, gmp)))
    {
       vtxtPrintf (blkp, ", title ") ;
       vtxtPrintf (blkp, "\"%s\"", ptr) ;
    }
  /* citations */
  if (0 && /* tatiana */
      (ptr = ficheAsnPubSet (buf, gmp->gene, "Reference")))
    {
      vtxtPrintf (blkp, ", cit pub {") ;
      vtxtPrintf (blkp, ptr) ;
      vtxtPrintf (blkp, "}") ;
    }

  /* experimental, always true in class mrna */
  if (gtIsExperimentalGene (gmp->gene, gmp->pg))
    vtxtPrintf (blkp, ", exp-ev experimental") ;
  
  /* the predictedgene dbxref */
  ficheMrnaDbXref (blkp, gmp, gmp->gene, 0, ", dbxref ") ;
  
  vtxtPrintf (blkp, "}") ; /* close ficheAsnPrdictedProductGenomeLocation */

  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
}  /* ficheAsnPredictedMrnaGenomeLocation */

/*********************************************************************/

static char *ficheAsnPredictedProductGenomeLocation (vTXT blkp, GMP *gmp)
{
  char *ptr ;

  vTXT buf ;   buf = vtxtCreate () ;
  
  vtxtClear (blkp) ;
  
  vtxtPrintf (blkp, "{ data cdregion { frame 1, code { id 1 }}" ) ;
  
  /* partial   */
  if (ac_has_tag (gmp->pg, "Start_not_found") || ac_has_tag (gmp->pg, "End_not_found"))
    vtxtPrintf (blkp, ", partial TRUE") ;

  /* comment */
  if (0) vtxtPrintf (blkp, 								   
	    " , comment \"ZZZZZZ THIS IS A predicted COMMENT ZZZZZZZZZZZZZZ\""
	    ) ;
  /* product of process */
  vtxtPrintf (blkp, ",  product  whole %s" , ficheAsnId ( gtGfName (gmp->pg), ".p")) ;
  
    /* export location */	
  if ((ptr = ficheAsnMrnaLocation (buf, gmp, 1)))
    vtxtPrintf (blkp, 								   
	      ", location mix {"
	      "  %s"
	      "}"
	      , ptr) ;
  
  /* title */
  if (1 && /* kans does not like title */
      (ptr = gtPredictedProductName (buf, gmp, gmp->pg, TRUE)))
    {
       vtxtPrintf (blkp, ", title ") ;
       vtxtPrintf (blkp, "\"%s\"", ptr) ;
    }
  /* citations */
  if (0 && /* tatiana */
      (ptr = ficheAsnPubSet (buf, gmp->gene, "Reference")))
    {
      vtxtPrintf (blkp, ", cit pub {") ;
      vtxtPrintf (blkp, ptr) ;
      vtxtPrintf (blkp, "}") ;
    }

  /* experimental, always true in class mrna */
  if (gtIsExperimentalGene (gmp->gene, gmp->pg))
    vtxtPrintf (blkp, ", exp-ev experimental") ;
  
  /* the predictedgene dbxref */
  ficheMrnaDbXref (blkp, gmp, gmp->gene, gmp->product, ", dbxref ") ;
  
  vtxtPrintf (blkp, "}") ; /* close ficheAsnPrdictedProductGenomeLocation */

  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
}  /* ficheAsnPredictedProductGenomeLocation */

/*********************************************************************/

/*********************************************************************/
/*********************************************************************/

/* -== == == == == == == == == == == == == == == == == = - /
   * -=  GENE DUMP = - /
   * -== == == == == == == == == == == == == == == == == = - *  
   *
   * We export here the geneBox accoring to the following schema
   * then directly called from this function
   * we export the mRNAs and the predictedGene as successive Seq-feat
   *
   Seq-feat ::= SEQUENCE {
   data SeqFeatData ,           -- the specific data // i.e gene
   comment VisibleString OPTIONAL , 
   location Seq-loc ,            -- feature made from
   title VisibleString OPTIONAL ,  == DEFINITION LINE in the gb flat file
   cit Pub-set OPTIONAL ,        -- citations for this feature
   exp-ev ENUMERATED {           -- evidence for existence of feature
   experimental (1) ,        -- any reasonable experimental check
   not-experimental (2) } OPTIONAL , -- similarity, pattern, etc
   dbxref SET OF Dbtag OPTIONAL ,  -- support for xref to other databases
   }
*/
static char *ficheAsnGeneSeqFeat (vTXT blkp, GMP *gmp, int gStart, int gEnd, 
				  const char *seqName, BOOL xref, BOOL pseudo)
{
  char *ptr = 0 ;
  int a1, a2 ;
  char *strand ;
  vTXT buf ;   buf = vtxtCreate () ;
    
  vtxtPrintf (blkp, 								   
	    "{"
	    "  data gene "
	    ) ;

  /* Export [Gene-ref] */
  ficheAsnGeneRef (blkp, gmp, xref) ;
  
   if (0) ptr =  ficheNewMrnaSubmissionComment (buf, gmp) ;
   else if (0) ptr =  gtGeneTitle (buf, gmp, TRUE) ;
   if (ptr)
     {
       vtxtPrintf (blkp, " , comment \"Title: %s\"", ptr) ;
     }
  vtxtClear (buf) ;
  /* export location */
  if (gStart < gEnd) { a1 = gStart ; a2 = gEnd ; strand = "plus" ; }
  else  { a2 = gStart ; a1 = gEnd ; strand = "minus" ; }
  vtxtPrintf (blkp, 								   
	    " , location int {"
	    "    from %d , "
	    "	 to %d , "
	    "    strand %s , "
	    "    id %s"
	    "  }"
	    , a1, a2, strand
	    , ficheAsnId (seqName, 0)) ;
  
  if (1 && /* kans does not like title */
      (ptr = gtGeneTitle (buf, gmp, TRUE)))
    {
      vtxtPrintf (blkp, "  , title ") ;
      vtxtPrintf (blkp, "\"%s\"", ptr) ;
    }
  

  if (0 && /* tatiana */
      (ptr = ficheAsnPubSet (buf, gmp->gene, "Reference")))
    {
      vtxtPrintf (blkp, ", cit pub {") ;
      vtxtPrintf (blkp, ptr) ;
      vtxtPrintf (blkp, "}") ;
    }
  /* experimental */
  
  /* the gene dbxref */
  if (1 && xref)
    ficheGeneDbXref (blkp, gmp, gmp->gene, 0, ", dbxref ") ;
  
  if (pseudo) vtxtPrintf (blkp, ", pseudo TRUE") ;

  /*fichePerGeneDumpSource (blkp, gmp->gene) ; */
  
  vtxtPrintf (blkp, 
	    "}"
	    ) ; 
  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
}  /* ficheAsnGeneSeqFeat */

/************************************************************************/
  /* resize the gene to incorporate all its mRNAs and predictions */

/*
* if there are mrna
*	use the min of the mrna starts and the gene start
*	use the max of the mrna ends and the gene end
* else
*	use the min of the predicted gene starts, ignoring the gene start
*	use the max of the predicted gene ends, ignoring the gene end
*
* if the gene start < gene end, the gene goes the other way; adjust
* your idea of min/max appropriately
* 
*/
static BOOL ficheResizeGeneBox (AC_OBJ oGene, AC_KEYSET mrnas, AC_KEYSET predictedGenes, int *gStartp, int *gEndp)
{
  AC_OBJ oTranscript ;
  AC_TABLE gIntMap, mIntMap ;
  int	ir, gStart, gEnd, mStart, mEnd, nMrna = 0, nPg = 0 ;
  BOOL isDown ;
  AC_HANDLE h = handleCreate () ;

  if (1)
    {
      nPg = predictedGenes ? ac_keyset_count (predictedGenes) : 0 ;
      nMrna = mrnas ? ac_keyset_count (mrnas) : 0 ;
    }

  gIntMap = ac_tag_table (oGene, "IntMap", h) ;
  gStart = ac_table_int (gIntMap, 0, 1, 0) - 1 ;
  gEnd = ac_table_int (gIntMap, 0, 2, 0) - 1 ;
  
  isDown = gStart < gEnd ? TRUE : FALSE ;
  ir = 0 ;
  if (nMrna)
    {
      AC_ITER iter = ac_keyset_iter (mrnas, TRUE, 0) ;
      while ((oTranscript = ac_next_obj (iter)))
	{
	  mIntMap = ac_tag_table (oTranscript, "IntMap", h) ;
	  mStart = ac_table_int (mIntMap, 0, 1, 0) - 1 ;
	  mEnd = ac_table_int (mIntMap, 0, 2, 0) - 1 ;
	  
	  if (0 && ir++== 0)
	    { gStart = mStart ; gEnd = mEnd ; }
	  if (isDown)
	    { 
	      if (mStart < gStart) gStart = mStart ;
	      if (mEnd > gEnd) gEnd = mEnd ;
	    }
	  else
	    { 
	      if (mStart > gStart) gStart = mStart ;
	      if (mEnd < gEnd) gEnd = mEnd ;
	    }
	  ac_free (oTranscript) ;
	}
      ac_free (iter) ;
    }
  
  if (nPg && !nMrna)
    {
      AC_ITER iter = ac_keyset_iter (predictedGenes, TRUE, 0) ;
      while ((oTranscript = ac_next_obj (iter)))
	{
	  mIntMap = ac_tag_table (oTranscript, "IntMap", h) ;
	  mStart = ac_table_int (mIntMap, 0, 1, 0) - 1 ;
	  mEnd = ac_table_int (mIntMap, 0, 2, 0) - 1 ;
	  
	  if (ir++== 0)
	    { gStart = mStart ; gEnd = mEnd ; }
	  if (isDown)
	    { 
	      if (mStart < gStart) gStart = mStart ;
	      if (mEnd > gEnd) gEnd = mEnd ;
	    }
	  else
	    { 
	      if (mStart > gStart) gStart = mStart ;
	      if (mEnd < gEnd) gEnd = mEnd ;
	    }
	  ac_free (oTranscript) ;
	}
      ac_free (iter) ;
    }

  *gStartp = gStart ;
  *gEndp = gEnd ;

  ac_free (h) ;
  return TRUE ;
} /* ficheResizeGeneBox */

/************************************************************************/

static int fichePerGeneSeqSetDump (vTXT  blkp, AC_DB db, AC_OBJ oGene, char *asnTag, char style)
{
  AC_KEYSET mrnas = 0, predictedGenes = 0 ;
  AC_OBJ oTranscript ;  
  AC_HANDLE h = handleCreate () ;
  int nMrna, nPg ;
  int isAny = 0 ;
  const char *chromName ;
  vTXT buf ;   buf = vtxtCreate () ;

  if (!ac_tag_obj (oGene, "Transcribed_gene", h) && !ac_tag_obj (oGene, "Genefinder", h) )
    { ac_free (h) ; return 0 ; }

  chromName = ac_tag_printable (oGene, "IntMap", 0) ;
  nMrna = nPg = 0 ;
  /*
    mrnas = ac_objquery_keyset (oGene, ">Transcribed_gene ; >mRNA ; !Gap_length && (COUNT {>cDNA_clone ; library = GB } > 0 || COUNT {>cDNA_clone ; >Read ; IsMrna } > 0)", h) ;
  */
  mrnas = ac_objquery_keyset (oGene, ">Transcribed_gene ; >mRNA ; product &&  ! bad_quality && ! gap", h) ;
  nMrna = ac_keyset_count (mrnas) ;
  predictedGenes = ac_objquery_keyset (oGene, ">genefinder ; CLASS predicted_gene && ! IS *ws6 && ! bad_quality", h) ;
  nPg = ac_keyset_count (predictedGenes) ;
  if (!nMrna && !nPg)
    goto abort ;

  /* export the transcribed_gene */
  if (nMrna)
    {
      AC_ITER iter = ac_keyset_iter (mrnas, TRUE, 0) ;
      while ((oTranscript = ac_next_obj (iter)))
	{
	  char *result ;
	  GMP *gmp = gmpCreate (db, oGene, 0, oTranscript, 0, 0, style, 'm') ;
	  vTXT buf = vtxtCreate () ;
	  
	  result = fAsnGenerateMRNA (buf, gmp, 0, 0) ;
	  if (result)
	    {
	      vtxtComma (blkp) ;
	      vtxtPrintf (blkp, "%s", result) ;
	    }
	  gmpDestroy (gmp) ;
	  ac_free (oTranscript) ;
	  vtxtDestroy (buf) ;
	}
      ac_free (iter) ;
    }
    
  /* export the predicted genes */
  if (!nMrna && predictedGenes) 
    {
      AC_ITER iter = ac_keyset_iter (predictedGenes, TRUE, 0) ;
      while ((oTranscript = ac_next_obj (iter)))
	{
	  GMP *gmp = gmpCreate (db, 0, 0, 0, oTranscript, 0, style, 'm') ;
	  isAny = 1 ;
	  
	  if (gmp)
	    {
	      if (ac_has_tag (oTranscript, "CDS"))
		{
		  if (gmp->mrna)
		    fichePerGeneFinderNucProtDump (blkp, chromName, gmp, &isAny) ;
		}
	      else
		fichePerGeneFinderTrnaDump (blkp, chromName, gmp, &isAny) ;
	      gmpDestroy (gmp) ;
	    } 
	  ac_free (oTranscript) ;
	}
      ac_free (iter) ;
    }
 abort:  
  vtxtDestroy (buf) ;
  ac_free (h) ;

  return isAny ;
} /* fichePerGeneNucProtDump */

/************************************************************************/

static int fichePerGeneLocationDump (vTXT  blkp, AC_DB db, GMP *gmp, char *asnTag, char style)
{
  AC_KEYSET mrnas = 0, predictedGenes = 0 ;
  AC_OBJ oMrna, oGF, oGene = gmp->gene ;
  int gStart, gEnd, nMrna, nPg ;
  int isAny = 0 ;
  char *ptr ;
  AC_HANDLE h = handleCreate () ;
  AC_TABLE oChrom = ac_tag_table (oGene, "IntMap", h) ;
  const char *chromName ;
  vTXT buf ;   buf = vtxtCreate () ;

  if (!ac_has_tag (oGene, "Transcribed_gene") && !ac_has_tag (oGene, "Genefinder") )
    { ac_free (h) ; return 0 ; }
  
  if (!oChrom)
    messcrash ("Missing IntMap in gene %s", ac_name (oGene)) ;
  chromName = ac_tag_printable (oGene, "IntMap", 0) ;

  nMrna = nPg = 0 ;
  mrnas = ac_objquery_keyset (oGene
			      , ">Transcribed_gene ; >mRNA ; product && ! bad_quality && ! gap"
			      , h) ;
  nMrna = ac_keyset_count (mrnas) ;
  predictedGenes = ac_objquery_keyset (oGene
				       , ">genefinder ;  CLASS predicted_gene && ! IS *ws6 && ! bad_quality"
				       , h) ;
  nPg = ac_keyset_count (predictedGenes) ;
  /*
   * these are reference returns from acFetch
   */
  if (!nMrna && !nPg)
    goto abort ;

  ficheResizeGeneBox (oGene, mrnas, predictedGenes, &gStart, &gEnd) ;

	/*
	* gStart and gEnd are the output of that function.
	* They are the min, max coordinates of anything in the gene.
	*/

  /* export the geneBox */
  
  if (1)
    {
      BOOL pseudo = FALSE ;
      AC_KEYSET ks2 = 0 ;

      if (isAny++)vtxtPrintf (blkp, ", ") ;
      else { if (asnTag) vtxtPrintf (blkp, asnTag) ; }
      if (nPg == 1 && nMrna == 0 && 
	  (ks2 = ac_ksquery_keyset (predictedGenes, "Pseudogene", 0)) &&
	  ac_keyset_count (ks2))
	pseudo = TRUE ;
      ac_free (ks2) ;
      ficheAsnGeneSeqFeat (blkp, gmp, gStart, gEnd, chromName, TRUE, pseudo) ;
    }
  
  /* export the mrnas and products */
  if (nMrna)
    {
      AC_ITER iter = ac_keyset_iter (mrnas, TRUE, 0) ;
      while ((oMrna = ac_next_obj (iter)))
	{
	  GMP *gmp2 = gmpCreate (gmp->db, 0, 0, oMrna, 0, 0, style, 'm') ;
	  gmp2->product = gtMrna2Product (oMrna, h) ;
	  if ((ptr = ficheAsnMrnaGenomeLocation (buf, gmp2)))
	    {  
	      if (isAny++)vtxtPrintf (blkp, ", ") ;
	      else { if (asnTag) vtxtPrintf (blkp, asnTag) ; } /* <- this case impossible */
	      vtxtPrintf (blkp, "%s", ptr ) ;
	    }
	  if (gmp2->product)
	    if ((ptr = ficheAsnProductGenomeLocation (buf, gmp2)))
	      {  
		if (isAny++)vtxtPrintf (blkp, ", ") ;
		else { if (asnTag) vtxtPrintf (blkp, asnTag) ; }
		vtxtPrintf (blkp, "%s", ptr ) ;
	      }
	  ac_free (gmp2->product) ;
	  gmpDestroy (gmp2) ;
	  ac_free (oMrna) ;
	}
      ac_free (iter) ;
    }
    
  /* export the predicted genes */
  if (!nMrna && predictedGenes)
    {
      AC_ITER iter = ac_keyset_iter (predictedGenes, TRUE, 0) ;
      while ((oGF = ac_next_obj (iter)))
	{
	  GMP *gmp = gmpCreate (db, 0, 0, 0, oGF, 0, style, 'm') ;
	  
	  if (gmp && 
	      (!ac_has_tag (oGF, "CDS") || gmp->mrna) &&
	      (ptr = ficheAsnPredictedMrnaGenomeLocation (buf, gmp)))
	    {
	      if (isAny++)vtxtPrintf (blkp, ", ") ;
	      else { if (asnTag) vtxtPrintf (blkp, asnTag) ; }
	      vtxtPrintf (blkp, "%s", ptr ) ;
	    }
	  if (gmp &&
	      ac_has_tag (oGF, "CDS") && gmp->mrna && 
	      (ptr = ficheAsnPredictedProductGenomeLocation (buf, gmp)))
	    {
	      if (isAny++)vtxtPrintf (blkp, ", ") ;
	      else { if (asnTag) vtxtPrintf (blkp, asnTag) ; }
	      vtxtPrintf (blkp, "%s", ptr ) ;
	    }
	  gmpDestroy (gmp) ; 
	  ac_free (oGF) ;
	} 
      ac_free (iter) ;
    }
 abort:  
  vtxtDestroy (buf) ;

  ac_free (h) ;
  return isAny ;
} /* fichePerGeneLocationDump */

/************************************************************************/

static void ficheAsnSeqSetAllGenes (vTXT blkp, AC_DB db, AC_OBJ oMap, char *testGene, char style)
{
  int iG, nn = 0 ;
  AC_OBJ oGene ; 
  AC_HANDLE h = handleCreate () ;
  AC_TABLE gGene = ac_tag_table (oMap, "Gene_i", h) ;

  for (iG = 0 ; iG < gGene->rows && (nn < iGlimit || !iGlimit || testGene) ; nn++, iG += (!testGene && iGlimit) ? 50 : 1)
    {
      oGene = ac_table_obj (gGene, iG, 0, h) ;
      if (testGene && strcasecmp (testGene, ac_name (oGene))) continue ;
      if (0 && strcasecmp ("5M", ac_name (oGene)) > 0) continue ;
      /* ccp = ac_name (oGene) ;  for debugging */
      /* if (*cp+3 == '-' && strstr (cp+3, ".")) continue ; */
      fichePerGeneSeqSetDump (blkp, db, oGene, 0, style) ;
    }
  ac_free (h) ;
} /* ficheAsnSeqSetAllGenes  */

/************************************************************************/
/* this fucntion is the public interface used by the fiche_multiline_editor */
char *fAsnGenerateGene (AC_DB db, AC_OBJ oGene, int isFeature, char style)
{
  vTXT blk ; blk = vtxtCreate () ;
  
  if (oGene)
    {  
      GMP *gmp = gmpCreate (db, oGene, 0, 0, 0, 0, style, 'g') ;
      if (isFeature)
	fichePerGeneSeqSetDump (blk, db, gmp->gene, 0, style) ;
      else
	fichePerGeneLocationDump (blk, db, gmp, 0, style) ;
      gmpDestroy (gmp) ;
    }
  return vtxtPtr (blk) ;
}

/************************************************************************/
/************************************************************************/

static int ficheGetChromoDna (AC_DB db, AC_OBJ oMap, char **cpp, AC_HANDLE h)
{
  int level, len = 0 ;
  Stack s = 0 ;
  char *ptr = 0, *cp ;
  FILE *f = 0 ;
  
  if (!clipGenomic256)
    {
      /*
	vdirGetCurrent (fileName, sizeof (fileName)-1) ;
	sprintf (fileName+strlen (fileName), "/ASN/genomicDna/%s.dna", ac_name (oMap)) ;
	
	sprintf(fileName,"%s.dna",ac_name(oMap));
	
	
      */
      
      f = filopen (messprintf("ASN/genomicDna/%s", ac_name (oMap)), "dna", "r") ;
      if (!f)
	{
	  ac_command (db
		      , messprintf ("query find Sequence IS \"%s\"", ac_name (oMap))
		      , 0, h) ;
	  ac_command (db
		      , messprintf ("dna ASN/genomicDna/%s.dna", ac_name (oMap))
		      , 0, h) ;
	  f = filopen (messprintf("ASN/genomicDna/%s", ac_name (oMap)), "dna", "r") ;
	}
      if (!f)
	messcrash ("Could neither read nor create the dna file ASN/genomicDna/*.dna") ;
      
      s = stackHandleCreate (20000000,h) ;
      level = freesetfile (f, 0) ;
      if (freecard (level)) /* jump title line */
	while (freecard (level))
	  {
	    cp = freepos () ;
	    vtextUpperCase (cp) ;
	    catText (s, cp) ;
	  }
      freeclose (level) ;
      ptr = stackText (s, 0) ;
      len = strlen (ptr) ;
    }
  else
    { ptr = strnew ("ATGCATGCATGCATGCATGCATGCATGC", h) ; len = strlen (ptr) ; }
	/* count the length and uppercase */

  *cpp = ptr ;

  return len ;
} /* ficheGetChromoDna */

/************************************************************************/

static void ficheAsnChromoComment (vTXT blkp)
{
  vtxtPrintf (blkp, "\"") ;
  vtxtPrintf (blkp, 
#ifdef TEST
	     "ZZZZZZZZZZ Chromosome comment ZZZZZZZZZZ"
#else
	     "The annotations are based on data from the C.elegans Genome Sequencing Consortium, "
	     "the Worm Transcriptome Project, the Worm ORFeome project and "
	     "other publicly available sources\n\n"
	     "The genomic sequence, imported from Wormbase release WS130 (Aug 15 2004), "
	     "has been produced by the C.elegans Genome Sequencing Consortium, "
	     "led by John Sulston and Bob Waterston.\n"
	     "Their coding sequences were predicted by "
	    "integration and manual review of the following data : computer analysis "
	    "using the program Genefinder (P. Green and L. Hillier, personal communication), "
	    "the large scale EST project of Yuji Kohara "
	    " (http://www.ddbj.nig.ac.jp/c-elegans/html/CE_INDEX.html) "
	    "the C.elegans ORFeome cloning project led by Marc Vidal (http://worfdb.dfci.harvard.edu/), "
	    "as well as similarity to other proteins using BlastX analyses (http://blast.wustl.edu/), "
	    "sequence conservation with C. briggsae using Jim Kent's WABA alignment program "
	    " (Genome Research 10:1115-1125, 2000), individual C. elegans GenBank submissions, "
	    "and personal communications with C. elegans researchers. tRNAs were predicted "
	    "using the program tRNAscan-SE (Lowe, T.M. and Eddy, S.R., 1997, Nucl. Acids. "
	    "Res., 25, 955-964.). The mitochondria was sequenced by David Wolstenholme "
	    " (Genbank X54252).\n"
            "NCBI annotations consisted in replacing, whenever possible, the predicted "
	     "CDS by the complete mRNAs available in Genbank and in enhancing the analysis of "
            "all proteins, predicted or confirmed by the Kohara Worm Transcriptome Project. "
	    "The analysis uses BLASTP, TaxBlast, PFAM (Bateman et al., Nucleic "
	    "Acids Res  2002 Jan 1 ;30 (1):276-80), PSORT2 (Nakai, K. and Kanehisa, M., "
	    "Genomics 14, 897-911, 1992 ; Nakai, K., Advances in protein chemistry 54, "
	    "277-344, 2000) and Expasy (pI and molecular weight). COG were provided for "
	    "2841 predicted CDS by Tatusov et al. (Nucleic Acids Res 29, 22-28, 2001). "
	    "More information about those genes is available from LocusLink, AceView/WormGenes and Wormbase."
#endif
	    ) ;

  vtxtPrintf (blkp, "\"") ;
} /* ficheAsnChromoComment */

/************************************************************************/
/*
  Bioseq ::= SEQUENCE {
  id SET OF Seq-id ,            -- equivalent identifiers
  descr Seq-descr OPTIONAL , -- descriptors
  inst Seq-inst ,            -- the sequence data
  annot SET OF Seq-annot OPTIONAL }
*/

static void ficheAsnSeqChromo (vTXT blkp, AC_DB db, AC_OBJ oMap, char *testGene, char style)
{
  char *iupac ;
  AC_HANDLE h = handleCreate () ;  
  int len = ficheGetChromoDna (db, oMap, &iupac, h) ;

  vtxtPrintf (blkp, 
	    "seq {\n"
	    "  id { %s } , "
	    "  descr {"
	    "    molinfo {"
	    "      biomol genomic }"
	     , ficheAsnId (ac_name (oMap), 0)
	    ) ; 

  vtxtPrintf (blkp, ", name ") ;
  vtxtPrintf (blkp, "\"Caenorhabditis elegans chromosome %s, complete sequence.\""
	    , niceChromoName (oMap)
	    ) ; 

  vtxtPrintf (blkp, ", title ") ;
  vtxtPrintf (blkp, "\"Caenorhabditis elegans chromosome %s, complete sequence.\""
	    , niceChromoName (oMap)
	    ) ; 

  vtxtPrintf (blkp, ", comment ") ;
  ficheAsnChromoComment (blkp) ;
  vtxtPrintf (blkp, 
	      ", user {"
	      "    type str \"RefGeneTracking\","
	      "    data {"
	      "           {"
	      "             label str \"Status\","
	      "             data str \"Reviewed\""
	      "            },"
	      "           {"
	      "             label str \"Assembly\","
	      "             data fields {"
	      "	                          {"
	      "                             label id 0,"
	      "	                             data fields {"
	      "                                           {"
	      "                                              label str \"name\","
	      "                                              data str \"WormBase genome release WS97\""
	      "                                            }"
	      "	                                          }"
	      "                            }"
	      "                          }"
	      "              }"
	      "            }"
	      " 	}"
	      );
  

  vtxtPrintf (blkp, 
	      ", pub {"
	      "pub {"
	      "article {"
	      "title {"
	      "name \"Genome sequence of the nematode C. elegans: a platform for "
	      "investigating biology.\" } ,"
	      "authors {"
	      "names "
              "std {"
	      "{  name  consortium \"The C. elegans Sequencing Consortium\" } } ,"
	      " affil std { affil "
              " \"The Washington University Genome Sequencing Center, Box "
	      "8501, 4444 Forest Park Parkway, St. Louis, MO 63108, USA. "
	      "worm@watson.wustl.edu\" }} ,"
	      "from journal {   title { iso-jta \"Science\" } ,"
              "imp { date   std {  year 1998 ,    month 12 , day 11 } ,"
	      " volume \"282\" ,"
	      " issue \"5396\" ,   pages \"2012-2018\" } } } ,"
	      "pmid 9851916 ,"
	      "muid 99069613 }}"
	      ) ;

  vtxtPrintf (blkp, 
	      /* 	    ", pub { pub { pmid 9851916 } }" */
	    "    } , "
	    "  inst {"
	    "    repr raw , "
	    "    mol dna , "
	    "    length %d , "
	    "    seq-data"
	    "      iupacna \""
	    , len
	    ) ;
  /* the actual dna */  
  vtxtPrintf (blkp, iupac) ;  

  ac_free (h) ;
  vtxtPrintf (blkp, " \"   }") ;/* closing inst */
  ficheAsnChromoAnnot (blkp, db, oMap, testGene, style) ; /* directly in blkp */
  vtxtPrintf (blkp, " }") ;/* closing seq */
} /* ficheAsnSeqChromo  */

/************************************************************************/
/************************************************************************/

static void ficheAsnChromoSeqSet (vTXT blkp, AC_DB db, AC_OBJ oMap, char *testGene, char style)
{
  /*directly in blkp :: do not  vtxtClear (blkp) ; */
  vtxtPrintf (blkp, ",  seq-set {") ; /* open seq-set */

  ficheAsnSeqChromo (blkp, db, oMap, testGene, style) ;
/*--*/
  ficheAsnSeqSetAllGenes (blkp, db, oMap, testGene, style) ;

  vtxtPrintf (blkp, "}") ; /* close  seq-set */
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
/*
Seq-feat ::= SEQUENCE {
    id Feat-id OPTIONAL , 
    data SeqFeatData ,           -- the specific data
    comment VisibleString OPTIONAL , 
    product Seq-loc OPTIONAL ,    -- product of process
    location Seq-loc ,            -- feature made from
    cit Pub-set OPTIONAL ,        -- citations for this feature
    exp-ev ENUMERATED {           -- evidence for existence of feature
        experimental (1) ,        -- any reasonable experimental check
        not-experimental (2) } OPTIONAL , -- similarity, pattern, etc
*/

static void ficheAsnSeqFeaturesChromo (vTXT blkp, AC_DB db, AC_OBJ oMap)
{
  vTXT buf ; buf = vtxtCreate () ;
  
  vtxtPrintf (blkp, "{") ;
  vtxtPrintf (blkp, "data region \"Whole chromosome\"") ;
  vtxtPrintf (blkp, ", comment ") ;
  ficheAsnChromoComment (blkp) ;
  vtxtPrintf (blkp, ", product whole %s", ficheAsnId (ac_name (oMap), ".seq")) ;
  vtxtPrintf (blkp, ", location id %s", ficheAsnId (ac_name (oMap), 0)) ;
  vtxtPrintf (blkp, 
	      ", cit  pub {"
	      "pub {"
	      "article {"
	      "title {"
	      "name \"Genome sequence of the nematode C. elegans: a platform for"
	      "investigating biology.\" } ,"
	      "authors {"
	      "names "
              "std {"
	      "{  name  consortium \"The C. elegans Sequencing Consortium\" } } ,"
	      " affil std { affil "
              "  \"The Washington University Genome Sequencing Center, Box"
	      "8501, 4444 Forest Park Parkway, St. Louis, MO 63108, USA."
	      "worm@watson.wustl.edu\" }} ,"
	      "from journal {   title { iso-jta \"Science\" } ,"
              "imp { date   std {  year 1998 ,    month 12 , day 11 } ,"
	      " volume \"282\" ,"
	      " issue \"5396\" ,   pages \"2012-2018\" } } } ,"
	      "pmid 9851916 ,"
	      "muid 99069613 }}"
	      ) ;
 
  vtxtPrintf (blkp, ", exp-ev experimental") ;
  vtxtPrintf (blkp, "}") ; /* close  SET of seq feat */
  vtxtDestroy (buf) ;
}  /* ficheAsnSeqFeaturesChromo */

/************************************************************************/
/*
Seq-annot ::= SEQUENCE {
   data  ftable SET OF Seq-feat }
*/
static void ficheAsnAnnotOne (vTXT blkp, AC_DB db, AC_OBJ obj, int type, char style)
{
  vTXT buf ; buf = vtxtCreate () ;

  switch (type)
    {
    case 1: /* chromo */
      vtxtPrintf (blkp, "{") ; /* open SEQUENCE */
      vtxtPrintf (blkp, "data ftable {") ; /* open SET of seq feat */
      ficheAsnSeqFeaturesChromo (blkp, db, obj) ;
      vtxtPrintf (blkp, "}") ; /* close  SET of seq feat */
      vtxtPrintf (blkp, "}") ; /* close  Sequence et of Seq-annot */
      break ;
    case 2:  /* gene */
      {
	GMP *gmp = gmpCreate (db, obj, 0, 0, 0, 0, style, 'g') ;
	if (fichePerGeneLocationDump (blkp, db, gmp, "{data ftable {", style))
	  {
	    vtxtPrintf (blkp, "}") ; /* close  SET of seq feat */
	    vtxtPrintf (blkp, "}") ; /* close  Sequence set of Seq-annot */
	  }
	gmpDestroy (gmp) ;
      }
      break ;
    }
 
  vtxtDestroy (buf) ;

} /* ficheAsnAnnotOne  */

/************************************************************************/

static void ficheAsnAnnotAllGenes (vTXT blkp, AC_DB db, AC_OBJ oMap, char *testGene, int isAny, char style)
{
  int iG, nn = 0  ;
  AC_HANDLE h = handleCreate () ;
  AC_OBJ oGene ;
  AC_TABLE gGene = ac_tag_table (oMap, "Gene_i", h) ;

  for (iG = 0 ; iG < gGene->rows && (nn < iGlimit || !iGlimit || testGene) ; iG += (!testGene && iGlimit) ? 50 : 1)
    {
      oGene = ac_table_obj (gGene, iG, 0, h) ;
      if (testGene && strcasecmp (testGene, ac_name (oGene)))continue ; 
      if (isAny++) vtxtPrintf (blkp, ", ") ;
      ficheAsnAnnotOne (blkp, db, oGene, 2, style) ;
      nn++ ;
    }
  ac_free (h) ;
	
} /* ficheAsnAnnotAllGenes  */

/************************************************************************/

static void ficheAsnAnnotChromo (vTXT blkp, AC_DB db, AC_OBJ oMap, char style)
{ 
  ficheAsnAnnotOne (blkp, db, oMap, 1, style) ;
} /* ficheAsnAnnotChromo */

/************************************************************************/
/************************************************************************/
/* annot SET OF Seq-annot  */
static void ficheAsnChromoAnnot (vTXT blkp, AC_DB db, AC_OBJ oMap, char *testGene, char style)
{
  int isAny = 0 ;
  /*directly in blkp :: do not  vtxtClear (blkp) ; */
  vtxtPrintf (blkp, ", annot {") ; /* open Set of Seq-annot */

  if (0) { ficheAsnAnnotChromo (blkp, db, oMap, style) ; isAny = 1 ; }
  ficheAsnAnnotAllGenes (blkp, db, oMap, testGene, isAny, style) ;

  vtxtPrintf (blkp, "}") ; /* close  Set of Seq-annot */

} /* ficheAsnChromoAnnot  */

/************************************************************************/
/************************************************************************/
/************************************************************************/

static void ficheDescrWormSource (vTXT blkp)
{
#ifdef JUNK
  this may help for human
 , gmp->Spc == HUMAN ? taxname " Homo sapiens"
	       , gmp->Spc == HUMAN ? id "9606" : "6239"
	       , gmp->Spc == HUMAN ? "Homo" :"Caenorhabditis"
	       , gmp->Spc == HUMAN ? "sapiens" :"elegans"
	       , gmp->Spc == HUMAN ? "" :",          div \"INV\""
#endif
  vtxtPrintf (blkp, 
	     "      genome genomic , " 
	     "      origin natural , "
	     "      org {"
	     "           taxname \"Caenorhabditis elegans\" , "
	     "           common \"worm\" , "
	     "           db {"
	     "                {"
	     "                  db \"taxon\" , "
	     "                  tag"
	     "                  id 6239 "
	     "               } } , "
	     "           orgname { "
	     "                     name"
	     "                     binomial {"
	     "                        genus \"Caenorhabditis\" , "
	     "                        species \"elegans\" } , "
	     "                     lineage \"Eukaryota ; Metazoa ; Nematoda ; Chromadorea ; Rhabditida ; Rhabditoidea ; Rhabditidae ; Peloderinae ; Caenorhabditis\" , "
	     "                     gcode 1 , "
	     "                     mgcode 5 , "
	     "                     div \"INV\""
	     "                    }"
	     "           } " 
	     ) ;
}



static  char *ficheAsnChromoDescrBioSource (vTXT blkp, AC_OBJ oMap)
{
  vtxtClear (blkp) ;

  vtxtPrint (blkp, "{" ) ;
  ficheDescrWormSource (blkp) ;
  vtxtPrintf (blkp, 
	     ", "
	     "subtype {"
	     "           { subtype chromosome, name \"%s\" } "
	     "       }" 
	     , niceChromoName (oMap)
	     ) ;
  vtxtPrint (blkp, 
	     "}"
	     ) ;

  return vtxtPtr (blkp) ;
} /* ficheAsnChromoDescrBioSource  */

/************************************************************************/

static char *ficheAsnChromoPub (vTXT blkp)
{
  char *ptr ;
  vTXT buf ; buf = vtxtCreate () ;
  ptr = ficheAsnDate (buf) ;
  /* submission citation and Consortium paper */
  vtxtPrintf (blkp, 
	      "{" /* open a sequence of Pub-equiv */
	      "      pub {"
	      "        sub {"
	      "          authors {"
	      "            names"
	      "            std {"
	      "              {"
	      "                name"
	      "                name {"
	      "                  last \"The NCBI REFSEQ group.\" } } } , "
	      "            affil"
	      "            std {"
	      "              affil \"National Center for Biotechnology Information, NIH\","
	      "              city \"Bethseda\" , "
	      "              sub \"MD\","
	      "              country \"USA\" , "
	      "              street \"Bdg 38A, 8600, Rockville Pike\" , "
	      "              email \"mieg@ncbi.nlm.nih.gov\" , "
	      "              postal-code \"20894\""
	      "                 } } , "
	      "          medium other, "
	      "          date %s "	
	      "        }" /*, close submission */
	      "        "  /* pmid 9851916 THE c.elegans paper Science 1998 */
	      "       }" /* close  pub (submission type) */
	      " }"   /* close */
	      , ptr) ;
  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
} /* ficheAsnChromoPub */

/***************************************
Seqdesc ::= CHOICE {
    name VisibleString ,         -- a name for this sequence
    title VisibleString ,        -- a title for this sequence
  OBSOLETE  org Org-ref ,                -- if all from one organism
    comment VisibleString ,      -- a more extensive comment
    pub Pubdesc ,                -- a reference to the publication
    create-date Date ,           -- date entry first created/released
    update-date Date ,           -- date of last update
    source BioSource }           -- source of materials, includes Org-ref

Org-ref ::= SEQUENCE {
    taxname VisibleString OPTIONAL ,   -- preferred formal name
    common VisibleString OPTIONAL ,    -- common name
    mod SET OF VisibleString OPTIONAL , -- unstructured modifiers
    db SET OF Dbtag OPTIONAL ,         -- ids in taxonomic or culture dbases
    syn SET OF VisibleString OPTIONAL ,  -- synonyms for taxname or common
    orgname OrgName OPTIONAL }

***************************/

static char *ficheAsnChromoDescr (vTXT blkp, AC_OBJ oMap)
{
  char *ptr ;
  vTXT buf ; buf = vtxtCreate () ;
  vtxtClear (blkp) ;
  
  /* name title org */
  vtxtPrintf (blkp, 
	     "{"   /* open descr */
	      "     title \"C.elegans chromosome %s complete annotated sequence.\""
	      , niceChromoName (oMap) 
	     ) ;

   /* skip the comment of the chromosome in the header part */
 
   /* pub */
  if ((ptr = ficheAsnChromoPub (buf)))
    vtxtPrintf (blkp, 
	       " , pub %s" /* expect a sequence of pub Pub-equiv */
	       , ptr) ;
   /* dates */
  if ((ptr = ficheAsnDate (buf)))
    vtxtPrintf (blkp, 
	       " , create-date %s"
	       , ptr) ;

   /* source */
  if ((ptr = ficheAsnChromoDescrBioSource (buf, oMap)))
    vtxtPrintf (blkp, 
	       " , source%s"
	       , ptr) ;
	  
  vtxtPrintf (blkp, "}") ; /* close descr */
  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
} /* ficheAsnChromoDescr */
 
/************************************************
 * Seq-entry ::= CHOICE {
 *      seq Bioseq , 
 *       set Bioseq-set } 
 *
 *
 * Bioseq-set ::= SEQUENCE {      -- just a collection
    class gen-prod-set , 
    release VisibleString OPTIONAL , 
    date Date OPTIONAL , 
    descr Seq-descr OPTIONAL , 
    seq-set SEQUENCE OF Seq-entry , 
    annot SET OF Seq-annot OPTIONAL }


 *************************************************/

static void ficheAsnSeqEntry (vTXT blkp, AC_DB db, AC_OBJ oMap, char *testGene, char style)
{ 
  char *ptr ;
  vTXT buf ; buf = vtxtCreate () ;
  
  vtxtPrintf (blkp, 
	      "Seq-entry ::= set { "
	      "class gen-prod-set"
	      ) ;
  /* info in the the descr is echoed to all the nuc-prot of the submission */
  if ((ptr = ficheAsnDate (buf)))
    vtxtPrintf (blkp, 
		", date %s"
		, ptr) ;
  
  if ((ptr = ficheAsnChromoDescr (buf, oMap)))
    {
      vtxtPrintf (blkp, ", descr ") ;
      vtxtPrintf (blkp, "%s", ptr) ;
    }  
  ficheAsnChromoSeqSet (blkp, db, oMap, testGene, style) ; /*directly in blkp */
   
  vtxtPrintf (blkp, "}") ; /* close the seq-entry */
  
  vtxtDestroy (buf) ;
} /* ficheAsnSeqEntry  */

/************************************************************************/
/************************************************************************/

char *testGene = 0 ;
/*

static void fiche_magic (char *s)
{
  printf("fiche magic:\n");
  if (! s)
    return;
  printf("\t%s\n",s);
  if (*s == '-')
    {
      clipGenomic256 = atoi(s) ;
      return;
    }
  if (*s == '+')
    {
      iGlimit = atoi(s);
      return;
    }
  if (*s == '/')
    {
    }
  testGene = strdup(s);
}

*/

char *ficheChromosomeDump (vTXT vtxt, AC_DB db, AC_OBJ oMap, char style)
{
  int nn ;
  AC_KEYSET aks = 0 ;
  vTXT blk ; blk = vtxtCreate () ;

  {
    AC_ITER lTmp ;
    AC_OBJ clo ;

    if ((lTmp = ac_query_iter (db, 0, "find Clone Main_clone", 0, 0)) &&
	(clo = ac_next_obj (lTmp)))
      {
	strcpy (genomeRelease, ac_tag_text (clo , "Genome_release", "WS")) ;
	ac_free (clo) ;
      }
    ac_free (lTmp) ;
  }

  if (0) testGene = "let-2" ;   /* non coding */
  if (0) testGene = "3a784" ;   /* has_ost */
  if (0) testGene = "3k842" ;   /* has_exact_est */
  if (0) testGene = "1a661" ;   /* */
  if (0) testGene = "unc-11" ;   /*   on chromosome I */
  if (0) testGene = "ceh-20C" ;   /* ceh-10  on chromosome III */
  if (0) testGene = "3f903" ;   /* ceh-10  on chromosome III */
  if (0) testGene = "pyc-1" ;   /* mir-35  rop-1  on chromosome V */
  if (0) testGene = "him-4" ;   /* ace-1 him-4 on X */
  if (0) testGene = "cca-1" ;   /* pme-6 ace-1  on chromosome X */
  if (0) testGene = "XE582" ;   /* with antisens */
  if (0) testGene = "ubq-2" ;   /*mab-21 == AM_gene, XF264 2G1 pme-6 ace-1  on chromosome X */
  if (0) testGene = "daf-5" ;
  if (0) testGene = "XB553" ;   /* very simple gene on X */
  if (0) testGene = "dyn-1" ;   /* very simple gene on X */
  if (0) testGene = "XA939" ;   /* 5' partial single exon on X */
  if (0) testGene = "XA693" ;   /* 5' partial single exon on X */
  if (0) testGene = "XF584" ;   /* asn bug gene inconsistant */
  if (0) testGene = "3K90" ;    /* starts on L */
  if (0) iGlimit = 200 ;
  if (0) clipGenomic256 = TRUE ;

  isRefSeqDump  = 1 ; /* horrible static in this file, controls asn identifier style */

  nn = ac_keyset_count (aks = ac_objquery_keyset (oMap, messprintf (">Gene_i ; IS \"%s\"", testGene ? testGene : "*"), 0)) ;
  ac_free (aks) ;
  if (nn)
    ficheAsnSeqEntry (blk, db, oMap, testGene, style) ;
  
  {   /* horrible hack to remove some trailing comas */
    char *cq, *cp = vtxtPtr (blk) ;
    if (cp) while (*cp)
      {
	if (*cp == ',')
	  {
	    cq = cp ;
	    while (*++cq)
	      switch (*cq)
		{
		case ',':
		case '}':
		  *cp = ' ' ; goto done ; break ;
		case '\n':
		case ' ':
		  break ; /* iterate on cq */
		case 0:
		default:
		  goto done ; break ;
		}
	  }
      done:
	cp++ ;
      }
  }
  if (vtxtPtr (blk))
    vtxtPrintf (vtxt, "%s", vtxtPtr (blk)) ;
  vtxtDestroy (blk) ;
  return vtxtPtr (vtxt) ;
}

/************************************************************************/
/************************************************************************/
/************************************************************************/




