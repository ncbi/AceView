/*  File: webquery.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg. 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEMBLY extension to 
        the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
     A dedicated query system to drive the AceView web site
 * Exported functions:
 * HISTORY:
 * Created: Mars 2001 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */

#define ARRAY_CHECK
#define MALLOC_CHECK

#include <ctype.h>
#include "acedb.h"
#include "bitset.h"
/* for webquery */
#include "lex.h"
#include "a.h"
#include "session.h"
/* for dna */
#include "dna.h"
#include "cdna.h"
#include "query.h"
#include "freeout.h"
#include "longtext.h"
#include "../whooks/systags.h"
#include "../whooks/tags.h"
#include "../whooks/sysclass.h"
#include "../whooks/classes.h"
#include "pick.h"
#include "dict.h"
#include "session.h"
#include "vtxt.h"


#define OWQMAX 20000 /* this seems to lose only 10% of the triplets */

/* very ad hoc code to run as fast as possible a general search
taylored for the aceView interface to the worm database
I heavilly use the fact that names are consistent, avoiding
the follow clause and all sorts of disgusting hacks
*/

typedef struct owbStruct { int flag, extend, loop ; char *cl, *follow ; } OWB ;
/* flag = 0 if we search this class directly (c=..) but not generaly (c=0 || c=extended) 
 * 0x1 search these name in a general search
 * 0x2 take as is
 * 0x4 find product $1 that is, this is a tag in the object
 * 0x8 filter

 * extend
 * 0x1 search also *name
 * 0x2 search also name*
 * 0x3 search also *name*
 * 0x4 first cut away the .number version endding 
 * 0x8 first cut away the .alpha version endding 
 * 0x10 search also name.*
 * loop  Search asm* should find locus asmr AND gene asm.1
 * 0x0 stop before this class if classes above gave a result
 * 0x1 also search that class even if the one above gave results
*/
static OWB owbMenu [] =
{
  /* only these classes are used in the first blocking search
     ideally, we should be blocking only on 'informatics identifiers'
     where we iterate first on exact name, then on name* then *name* 
     according to what is allowed in the extension flag
  */
  { 9, 2, 0, "Locus",           "{>Gene} $/{>locusLink;>Gene} $/ {>Ref_id;>model_of_gene} $/ {>Sequence;>model_of_gene}"},
  { 9, 2, 1, "LocusLink",       "{>Gene} $/{>geneid;>Gene} $/{>geneid;>Gene} $/ {>sequence; >from_gene ; > gene} $/ {>geneid;>predicted_gene;>model_of_gene} $/ {>geneid;>predicted_gene;>Matching_transcribed_gene;>Gene}" } ,
  { 9, 2, 1, "NewName",  ">Gene"} ,
  { 3, 16, 1, "Gene", 0} ,
  { 9, 2, 1, "Gene_title",  ">Gene"} ,
  { 9, 0, 1, "PFAM" ,           ">Accession; >Quoted_in; > Gene"} , /* because pfam is also a gene name */
  { 9, 0, 0, "LocusId",         "{>Gene} $/ {>Locuslink;>Gene } $/ {>Predicted_gene;>model_of_gene} $/ {>Sequence;>model_of_gene} $/{>predicted_gene;>Matching_transcribed_gene;>Gene}  $/ {>Predicted_gene;>from_gene ; >gene}" } ,
  { 9, 0, 1, "GeneId",          "{>Gene} $/ {>Locuslink;>Gene } $/ {>Predicted_gene;>model_of_gene} $/ {>Sequence;>model_of_gene} $/{>predicted_gene;>Matching_transcribed_gene;>Gene}  $/ {>Predicted_gene;>from_gene ; >gene}" } ,
  { 9, 0, 1, "WbId",          ">Gene" } ,
  { 9,10, 0, "mRNA",              "{>From_gene ; >Gene} $/ {>From_prediction;>Model_of_gene}"} ,
  /* NB we extend predicted_gene in a special way in aceviewmain */
  { 9, 2, 0, "Predicted_gene",  "{>Model_of_gene} $/ {>Matching_transcribed_gene ; >Gene} $/ {geneid; > predicted_gene;>from_gene ; >gene}"} ,
  { 9, 4, 0, "GenBank",         "{>gene_nm} $/ {>Predicted_gene_NM ; >Model_of_gene}" } , /* nm better as above */
  { 9, 0, 1, "Clone_group",   "{>Contains}" },
  { 9, 0, 0, "Unigene",   ">Gene" },
  { 3, 5, 0, "cDNA_clone",   0 },
  { 9, 4, 0, "Sequence", "{>cdna_clone} $/ {IS XM* ; >locusLink;>sequence;>from_gene;>gene} " }, /* "Follow  From_gene ; >gene"} , */

  /* non blocking search, used only in complete search */
  { 8, 0, 0, "NewName",  " >GeneOld"} ,
  { 8, 2, 0, "Transcribed_gene", "{>Gene}"} ,
  { 8, 3, 0, "Locus_description",       "{>Gene } $| {>Gene_descriptor} $| {>Omim ; >Gene}"} ,
  { 8, 2, 0, "Product",         ">  GeneBox"} , /* psort is actually a tag directly defined in the product */
  /*  { 8, 0, 0, "PFAM" ,           "{> Gene} $/ { >Product ; >  GeneBox}"} , */
  { 8, 3, 0, "Mesh" ,           "{IS *} $| {>alias_of} $| {>hidden_alias} $| {>hidden_alias_of} ; {IS *} $| {>reference ; COUNT Gene < 10} $| {>Extern; !BioCarta } $| {>OMIM_title} $| {>GAD_title} $| {>KEGG_pathway KEGG_disease}  ; >Gene"} ,
  { 8, 3, 0, "Extern" ,           "!BioCarta ; {>Gene} $| {>geneid ; > gene} $|  {GeneId_mol; >gene}"} ,
  { 8, 3, 0, "Go_b" ,           ">gene_ace"} ,
  { 8, 3, 0, "Go_c" ,           "{>gene_ace} $| {>gene_psort}"} ,
  { 8, 3, 0, "Go_m" ,           " >gene_ace "} ,
  { 4, 3, 0, "Motif" ,        0} ,
  { 4, 3, 0, "Psort" ,        0} ,
  { 8, 0, 0, "Expr_pattern" ,        "> Gene"} ,
  { 8, 0, 0, "Strain" ,        "> Gene"} ,
  { 8, 0, 0, "RNAi" ,        "> Gene"} ,
  { 8, 0, 0, "COG" ,        "> Gene"} ,
  { 8, 0, 0, "OMIM" ,        "> Gene"} ,
  { 8, 0, 0, "Paper",       " COUNT Gene < 100 ; > Gene"} ,  /* note the 100 limit */
  { 8, 3, 0, "Author",      "> Gene"} ,  
  { 0, 0, 0, 0 }
} ;

/***************************************************************/

static KEYSET owqDefault (void)
{
  KEY key ;
  KEYSET ks ;
  
  lexaddkey ("__OMIM_genes_4", &key, _VKeySet) ;
  ks = arrayGet (key, KEY, "k") ;
  
  if (!ks)
    {
      ks = query (0, "{Find Gene LocusLink && pastille_disease && ! IS *and*} $/ {Find gene IS *-* && smap } ") ;
      if (keySetMax(ks) && sessionGainWriteAccess () )
	arrayStore (key, ks, "k") ;
    }
  return ks ;
}

/***************************************************************/
/*************** lazy evaluation *******************************/
/***************************************************************/

KEYSET OWQGene2NewName (KEYSET ks, BOOL force)
{
  int i, j, n = ks ? keySetMax (ks) : 0 ;
  KEYSET ks1 = keySetCreate () ;
  KEY gene, newname ;
  BOOL touched = FALSE ;
  static KEYSET g2nn = 0 ;
  static KEY g2nnKey = 0 ;

  if (!g2nn || force)
    {
      keySetDestroy (g2nn) ;
      lexaddkey ("OWQ_g2nn", &g2nnKey, _VKeySet) ;
      g2nn = arrayGet (g2nnKey, KEY, "k") ;
      if (!g2nn)
	g2nn = arrayCreate (lexMax (_VGene), KEY) ;
    }

  for (i = j = 0 ; i < n ; i++)
    {
      gene = keySet (ks, i) ;
      if (class (gene) == _VGene)
	{
	  newname = keySet (g2nn, KEYKEY (gene)) ;
	  if (!newname)
	    {
	      newname = keyGetKey (gene, _NewName) ;
	      if (!newname) 
		newname = 1 ;
	      keySet (g2nn, KEYKEY (gene)) = newname ;
	      touched = TRUE ;
	    }
	  if (class (newname))
	    keySet (ks1, j++) = newname ;
	}
    }

  if (touched && isWriteAccess())
    arrayStore (g2nnKey, g2nn, "k") ;
  keySetSort (ks1) ;
  keySetCompress (ks1) ;
  return ks1 ;
} /* OWQGene2NewName */

/***************************************************************/

KEYSET OWQProduct2Gene (KEYSET ks, BOOL force)
{
  int i, j, n = ks ? keySetMax (ks) : 0 ;
  KEYSET ks1 = keySetCreate () ;
  KEY gene, product ;
  BOOL touched = FALSE ;
  static KEYSET p2g = 0 ;
  KEY p2gKey = 0 ;
  KEY _GeneBox = str2tag ("GeneBox") ;

  if (!p2g || force)
    {
      keySetDestroy (p2g) ;
      lexaddkey ("OWQ_p2g", &p2gKey, _VKeySet) ;
      p2g = arrayGet (p2gKey, KEY, "k") ;
      if (!p2g)
	p2g = arrayCreate (lexMax (_VmProduct), KEY) ;
    }

  for (i = j = 0 ; i < n ; i++)
    {
      product = keySet (ks, i) ;
      if (class (product) == _VmProduct)
	{
	  gene = keySet (p2g, KEYKEY (product)) ;
	  if (!gene)
	    {
	      gene = keyGetKey (product, _GeneBox) ;
	      if (!gene) 
		gene = 1 ;
	      keySet (p2g, KEYKEY (product)) = gene ;
	      touched = TRUE ;
	    }
	  if (class (gene))
	    keySet (ks1, j++) = gene ;
	}
    }

  if (p2gKey && touched && isWriteAccess())
    arrayStore (p2gKey, p2g, "k") ;
  keySetSort (ks1) ;
  keySetCompress (ks1) ;
  return ks1 ;
} /* OWQProduct2Gene */

/***************************************************************/

static void owqStoreCannedList (KEYSET ks, char *cp, int type, BOOL isSaturated, FILE *fil) 
{ 
  KEY key ;
  int i ;
  char *suffix ;
  if (ks && isSaturated)
    {
      lexaddkey ("Text", &key, _VKeySet) ;
      keySetMax (ks) = 1 ;
      keySet (ks, 0) = key ;
    }
  switch (type)
    {
    case 1:  suffix = ".ext" ; break ; 
    case 3:  suffix = ".3" ; break ;
    case 5:  suffix = ".5" ; break ;
    default: suffix = "" ; break ;
    }
  if (cp && *cp && strlen(cp) >= 3)
    {
      if (fil)
	{
	  fprintf (fil, "KEYSET %s__%s%s\n"
		   , ks && keySetMax (ks) ? "" : "_"
		   , cp
		   , suffix
		   ) ;
	  for (i = 0 ; ks && i < keySetMax (ks) ; i++)
	    {
	      key = keySet (ks, i) ;
	      fprintf (fil, "%s:%s\n", className(key), freeprotect(name(key))) ;
	    }
	  fprintf (fil, "\n") ;
	}
      else if (sessionGainWriteAccess ())
	{
	  if (ks && keySetMax (ks))
	    {
	      lexaddkey (messprintf("__%s%s", cp, suffix), &key, _VKeySet) ;
	      arrayStore (key, ks, "k") ;
	    }
	  else
	    lexaddkey (messprintf("___%s%s", cp, suffix), &key, _VKeySet) ;
	}
    }
}

/***************************************************************/
/* set found = TRUE even if the set is known to be empty */
static KEYSET owqGetCannedList (char *cp, int *foundp, int type)
{ 
  KEY key ;
  KEYSET ks = 0 ;
  char *suffix ;

  *foundp = 0 ;
  switch (type)
    {
    case 1:  suffix = ".ext" ; break ; 
    case 3:  suffix = ".3" ; break ;
    case 5:  suffix = ".5" ; break ;
    default: suffix = "" ; break ;
    }

  if (cp && *cp)
    {
      if (lexword2key (messprintf("__%s%s", cp, suffix)
		       , &key, _VKeySet))
	{
	  if ((ks = arrayGet (key, KEY, "k")))
	    *foundp = 2 ;
	}
      else if (lexword2key (messprintf ("___%s%s", cp, suffix)
			    , &key, _VKeySet) &&
	       !(lexGetStatus (key) & EMPTYSTATUS))
	*foundp = 2 ;
    }
  if (ks && keySetMax (ks) == 1 && class (keySet (ks, 0)) == _VKeySet)
    *foundp = 1 ;
  return ks ;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

static KEYSET owqFilter (KEYSET ks, OWB *owb)
{
  KEYSET ks1 = 0, ks2 = 0 ;
  KEY key ;
  
  if (owb->flag & 0x2)   /* direct search */
    ks1 = query (ks, messprintf ("CLASS %s", owb->cl)) ;
  else if (owb->flag & 0x4)  /* search those as tags in genes and products */
    {
      int i, j ;
      char *cp ;
      KEY tag = 0 ;
      vTXT txt = vtxtCreate () ;

      vtxtPrintf (txt, "                 ") ; /* space for find ... */
      for (i = j = 0 ; i < keySetMax (ks); i++)
	{
	  key = keySet (ks, i) ;
	  if (strcasecmp (className (key), owb->cl) ||
	      !lexword2key (name(key), &tag, 0))
	    continue ;
	  if (j++)
	    vtxtPrintf (txt, " OR ") ;
	  vtxtPrintf (txt, name(tag)) ;
	}
      if (j)
	{
	  cp = vtxtPtr (txt) ;
	  memcpy (cp, "Find Gene", 9) ;   /* ATTENTION we are writing in the reserved blank space */
	  ks1 = query (0, cp) ;
	  if (!keySetMax (ks1))
	    {
	      keySetDestroy (ks1) ;

	      memcpy (cp, "Find Product", 12) ; /* ATTENTION word product is a longer than word gene */
	      ks2 = query (0, cp) ;
	      ks1 = OWQProduct2Gene (ks2, FALSE) ;
	      keySetDestroy (ks2) ;
	    }
	}
      vtxtDestroy (txt) ;
    }
  else if (owb->flag & 0x8) /* follow filter */
    {     
      ks2 = query (ks, messprintf ("CLASS %s", owb->cl)) ;
      ks1 = query (ks2, owb->follow) ;
    }
  keySetDestroy (ks2) ;
  return ks1 ;
}

/***************************************************************/
/* does the eventual *name* extension and follows the filters */
static KEYSET owqSearch (OWB *owb, char *qName, int ii, KEYSET ksCandidate, char *cl)
{
  KEYSET ks = 0, ks2 = 0 ;
  char cc = 0, *cp = 0, *qq = 0, *f_cl, *sep ;

  if (!cl && /* if cl, do not tamper with the ennding */
      (owb->extend & (0x4 | 0x8)))
    {
        /* clean up NM version number */
      cp = qName + strlen (qName) - 1 ;
      while (cp > qName  && *cp != '.') cp-- ;
      if (*cp == '.' && strlen (cp) < 3)
	{
	  if (cp > qName + 3 &&
	      (owb->extend & 0x4) &&   /* accession.version, kill .version */
	      isdigit((int)*(cp-1)) && isdigit((int)*(cp-2)) && isdigit((int)*(cp-3))	&&
	      isdigit((int)*(cp+1)))
	    { cc = *cp ; *cp = 0 ; }
	  else if ((owb->extend & 0x8) && /* mrna end as .a .b */
		   isalpha ((int)*(cp+1)))
	    { cc = *cp ; *cp = 0 ; }
	}

    }

  f_cl = ksCandidate ? "CLASS" : "Find" ;
  sep = ksCandidate ? "&&" : "" ;
  switch (ii)
    {
    case -1:
    case 0:
      if (!cl &&  (owb->extend & 0x10))
	qq = messprintf ("%s %s %s (IS \"%s\" OR IS \"%s.*\")", f_cl, owb->cl, sep, qName, qName) ;
      else
	qq = messprintf ("%s %s %s IS \"%s\"", f_cl, owb->cl, sep, qName) ;
      break ;
    case 1: /* right extend */
      if (owb->extend & 0x2)
	qq = messprintf ("%s %s %s IS \"%s*\"", f_cl, owb->cl, sep, qName) ;
      break ;
    case 2: /* left or double  extend */
      if ((owb->extend & 0x3)== 0x3)
	qq = messprintf ("%s %s %s IS \"*%s*\"", f_cl, owb->cl, sep, qName) ;
      else if (owb->extend & 0x1)
	qq = messprintf ("%s %s %s IS \"*%s\"", f_cl, owb->cl, sep, qName) ;
      break ;
    }
  
  if (qq)   /* should we search */
    ks2 = query (ksCandidate, qq) ;
  if (ks2 && keySetMax (ks2)) 
    ks = owqFilter (ks2, owb) ; 
  if (cc) /* restore the query after filtering */
    *cp = cc ;
  keySetDestroy (ks2) ; 
  if (ks && !keySetMax (ks))
    keySetDestroy (ks) ;
  return ks ;
}  /* owqSearch */

/***************************************************************/
typedef struct { KEY gene, map ; int x ; } OWQMAP ;

static int owqOrder (const void *va, const void *vb)
{
  const  OWQMAP *up = (const OWQMAP *)va, *vp = (const OWQMAP *)vb ;

  if (up->map && !vp->map)
    return -1 ;
  if (!up->map && vp->map)
    return 1 ;
  if (!up->map && !vp->map)
    return lexstrcmp (name(up->gene), name(vp->gene)) ;

  if (up->map != vp->map)
    return lexstrcmp (name(up->map), name(vp->map)) ;
  return up->x - vp->x ;
}

/***************************************************************/
/* search each indexed class untill first hit */
static KEYSET owqOneWord (char *cl, char *qNameExact, char *qName, KEYSET ksCandidate, BOOL isExact, BOOL debug)
{
  KEYSET ks = 0, ks1 = 0, ks2 = 0 ;
  int ii ;
  OWB *owb ;
  
  /* loop at -1, for classes with complex name: "cataract, marner like" */
  for (ii = qNameExact && strcasecmp (qNameExact, qName) ? -1 : 0 ; 
       !ks && ii < (isExact ? 1 : 3) ; ii++)   /* loop on name only, not on name* *name* */
    for (owb = owbMenu ; (!ks || owb->loop) && owb->cl ; owb++)
      if ( (!cl && (owb->flag & 0x1)) || (cl && !strcasecmp (cl, owb->cl)) )
	{
	  ks1 = ks ;
	  ks2 = owqSearch (owb, ii == -1 ? qNameExact : qName, ii,  ksCandidate, cl) ; 
	  if (ks1 && ks2)
	    {
	      ks = keySetOR (ks1, ks2) ;
	      keySetDestroy (ks1) ;
	      keySetDestroy (ks2) ;
	    }
	  else if (ks1)
	    ks = ks1 ;
	  else
	    ks = ks2 ;
	}
  return ks ;
} /* owqOneWord */

/**************/

static void owqSortByChromo (KEYSET ks)
{
  Array hits = 0 ;
  int i, nn = keySetMax(ks) ;
  OWQMAP *hh ;
  OBJ Gene = 0 ;

  if (!ks || keySetMax(ks) < 2)
    return ;

  hits = arrayCreate (nn, OWQMAP) ;
  i = nn ;
  while (i--)
    {
      hh = arrayp (hits, i, OWQMAP) ;
      hh->gene = keySet (ks, i) ;
      if ((Gene = bsCreate (hh->gene)))
	{
	  if (bsGetKey (Gene, _IntMap, &hh->map))
	    bsGetData (Gene, _bsRight, _Int, &(hh->x)) ;
	  bsDestroy (Gene) ;
	}
    }
  arraySort (hits, owqOrder) ;
  for (i = 0 ; i < nn ; i++)
    { 
      hh = arrayp (hits, i, OWQMAP) ;
      keySet (ks, i) = hh->gene ;
    }
  arrayDestroy (hits) ;
}

/***************************************************************/

static KEYSET owqChromSearch (char *mapName)
{
  KEYSET ks = 0 ;
  char *cp = 0, *locus1 = 0, *locus2 = 0, *longMapName ;
  int x1 = 0, x2 = 0, mode = 0, a1, a2, ii, jj ;
  KEY map = 0, gene ;
  OBJ Gene = 0 ;

  while ((cp = freeword()))
    {
      if (!strcmp (cp, "-coord"))
	{
	  if (freeint (&x1) &&
	      freeint (&x2))
	    mode |= 1 ;
	  else
	    x1 = x2 = 0 ;
	}
      else if (!strcmp (cp, "-locus1"))
	{
	  if ((cp = freeword()) && *cp != '-')
	    {
	      locus1 = strnew (cp, 0) ;
	      mode |= 2 ;
	    }
	}
      else if (!strcmp (cp, "-locus2"))
	{
	  if ((cp = freeword()) && *cp != '-')
	    {
	      locus2 = strnew (cp, 0) ;
	      mode |= 4 ;
	    }
	}
    }
  
  if (x1 < x2) { int xtmp = x1 ; x1 = x2 ; x2 = xtmp ; }

  /* note that genes are smaller and faster to search than tg */
  if (lexword2key (mapName, &map, _VMap))
    ks = query (0, messprintf ("Find Map IS %s OR IS *_%s ; >Gene_i", mapName, mapName)) ;
  
  if (keySetMax (ks) && (mode & 1))
    {
      longMapName = strnew (messprintf ("CHROMOSOME_%s",mapName), 0) ;
      for (ii = jj = 0 ; ii < keySetMax (ks) ; ii++)
	{
	  gene = keySet (ks, ii) ;
	  if ((Gene = bsCreate (gene)))
	    {
	      if (bsGetKey (Gene, _IntMap, &map) &&
		  (!strcmp (name(map), mapName) ||!strcmp (name(map), longMapName)) && 
		  bsGetData (Gene, _bsRight, _Int, &a1) &&
		  bsGetData (Gene, _bsRight, _Int, &a2) &&
		  (x1 <= a1 || x1 <= a2) &&
		  (x2 >= a1 || x1 >= a2) )
		{
		  if (ii != jj)
		    keySet (ks, jj) = gene ;
		  jj++ ;
		}
	      bsDestroy (Gene) ;
	    }
	}
      keySetMax (ks) = jj ;
      messfree (longMapName) ;
    }
 
  if (!keySetMax (ks)) 
    keySetDestroy (ks) ;
  else
    {
      keySetSort (ks) ;
      keySetCompress (ks) ;
      if (0) owqSortByChromo (ks) ;
    }
  messfree (locus1) ;
  messfree (locus2) ;
     
  return ks ;
}  /* owqChromSearch */

/***************************************************************/
/***************************************************************/
#ifdef JUNK
static KEYSET ficheViewGrep (char *text)
{
  KEYSET ks = query (0, "Find ficheview") ;
  int i, j = 0 ;
  KEY key ;

  if (1)
    for (i = j = 0 ; i < keySetMax (ks) ; i++)
      {
	key = keySet (ks, i) ;
	if (longTextGrep (key, text, 0))
	  keySet (ks, j++) = key ;
      }
  keySetMax (ks) = j ;
  
  return ks ;
}

/***************************************************************/

static int ficheViewFilter (KEYSET ks, Stack s, BOOL reMap) 
{
  KEY key, key1 ;
  int i, j ;
  
  for (i = j = 0 ; i < keySetMax (ks) ; i++)
    {
      key = keySet (ks, i) ;
      if (longTextGrep (key, 0, s))
	{
	  if (1 || /* remap is not correct for ficheview */
	      !reMap || !lexReClass(key, &key1, _VGene))
	    key1 = key ;
	  keySet (ks, j++) = key1 ;
	}
    }
  keySetMax (ks) = j ;
  return j ;
}
#endif

/***************************************************************/
  /* collect is autoKs the autoXrefed classes */

static KEYSET owqExtGrep (char *buffer, BOOL debug)
{
  KEYSET autoKs = keySetCreate() ; 
  int n = 0, t = 256 ;
  BOOL doInterrupt = FALSE ;
  KEY k ;
  char *cp ;

  while(!doInterrupt && t--)
    {
      if (pickXref(t))  /* pick all cross referenced objects */
	{ 
	  for (k = 0 ; lexNext (t,&k) ;)
	    { 
	      cp = name(k) ;

	      if (messIsInterruptCalled())
		{ doInterrupt = TRUE ;
		  break ;
		}

	      if (pickMatch (cp, buffer))
		keySet(autoKs,n++) = k ;
	    }
	}
    }
  if (n > 1)
    {
      keySetSort(autoKs) ;
      keySetCompress(autoKs) ;
    }
  return autoKs ;
} /* owqExtGrep */

/***************************************************************/
/* text must be 8 char long */
static KEYSET owqExt3 (char *text, int create, int *foundp, int type, BOOL debug, FILE *fil)
{
  KEYSET ks = 0, ks1, ks2, ks3, ks4, ks5 ;
  OWB *owb ;
  BOOL isSaturated = FALSE ;
  char kname[6] ;
  
  switch (type)
    {
    case 3:
      kname[0] = text[1] ;
      kname[1] = text[2] ;
      kname[2] = text[3] ;
      kname[3] = 0 ;
      break ;
    case 5:
      kname[0] = text[1] ;
      kname[1] = '_' ;
      kname[2] = text[3] ;
      kname[3] = '_' ;
      kname[4] = text[5] ;
      kname[5] = 0 ;
      break ;
    }

  ks = owqGetCannedList (kname, foundp, type) ;
  if (!fil && *foundp)
    {
      if (debug)
	printf ("owqExt3 Old %s %d entries %s\n"
		, kname, ks ?keySetMax(ks) : 0
		, *foundp == 1 ? "saturated" : ""
		) ;
      return ks ;
    }
  if (!create)
    return 0 ;

  for (owb = owbMenu ;  ; owb++)   /* search all classes  */
    {
      if (owb->cl)
	{
	  /* index all classes which should be searched by name
	   * either directly owb->flag & 0x1
	   * but of course also text like  owb->flag & 0x8
	   */
	  ks2 = query (0, messprintf("Find %s IS \"%s\"", owb->cl, text)) ;
	}
      else /* last round, do the Text search */
	{
	  if (debug) printf ("owqExt3 Grep %s ", text + 1) ;
	  ks1 = owqExtGrep (text, debug) ; 
	  if (debug) printf ("%d\n", ks1 ? keySetMax(ks1) : 0) ;
	  if (debug) printf ("owqExt3 LongGrep %s ", text + 1) ;
	  ks3 = (create == 2 ? longGrepPlain (text) : keySetCreate ()) ;
	  if (debug) printf ("%d\n", ks3 ? keySetMax(ks3) : 0) ;
	  ks4 = keySetOR (ks1, ks3) ; 
	  keySetDestroy (ks1) ; keySetDestroy (ks3) ;
	  ks5 = keySetCreate () ; /* ficheViewGrep (text) ; */
	  ks2 = keySetOR (ks4, ks5) ; 
	  keySetDestroy (ks4) ; keySetDestroy (ks5) ;
	}
      if (ks)
	{
	  ks1 = ks ;
	  ks = keySetOR (ks1, ks2) ;
	  keySetDestroy (ks1) ; keySetDestroy (ks2) ;
	}
      else  /* first time we initialise ks */
	{ ks = ks2 ; ks2 = 0 ; }
      if (!owb->cl)
	break ;
      if (ks && keySetMax (ks) >= 2 * OWQMAX)
	{ 
	  keySetSort (ks) ;
	  keySetCompress (ks) ;
	  if (keySetMax (ks) >= OWQMAX)
	    { isSaturated = TRUE ; goto foundMax ; }
	}
      if (!owb->cl) /* break here to write the above code only once */
	break ;
    }
  
  keySetSort (ks) ;
  keySetCompress (ks) ;
 foundMax:
  *foundp = 2 ; 
  if (keySetMax (ks) >= OWQMAX)
    {
      isSaturated = TRUE ; 
      *foundp = 1 ;  /* saturated */
    }
  owqStoreCannedList (ks, kname, type, isSaturated, fil) ;
  if (debug) printf ("owqExt3 New %3s %d entries %s\n"
		     , kname, ks ?keySetMax(ks) : 0
		     , *foundp == 1 ? "saturated" : "") ;
  return ks ;
} /* owqExt3  */

/***************************************************************/
#define isIndexed(_cc,_type) (_type==3?(isalnum((int)_cc)||((_cc)=='-')):(isalpha((int)_cc)||((_cc)=='-')))

/* hash exactly 3 bytes down to 6 bits == 64 possibilities */
static unsigned int OWQHash (char *cp)
{
  register int n = 6, i ;
  register unsigned int j, x = 0 ;
  register int rotate =  13 ;
  register int leftover = 8*sizeof(int) - rotate ;

  i = 3 ; while (i--)
    x = ace_upper (*cp++) ^ (( x >> leftover) | (x << rotate)) ; 
				/* compress down to n bits */
  for (j = x, i = n ; i < 8*sizeof(int) ; i += n)
    j ^= (x >> i) ;
  j &= (1 << n) - 1 ;
 
  return j ;
}  /*  OWQHash */

static int owqPrecompute (int mask, FILE *fil, DICT *dict)
{
  int type, nn, n2, nfound, nfound2, nsaturated ;
  int c1, c2, c3 ;
  char text[8], text2[3] ;
  KEYSET ks2 = 0 ;
  int found ;
  BOOL debug = TRUE ;
  
  memset (text, 0, sizeof(text)) ;
  memset (text2, 0, sizeof(text2)) ;
  mask &= 0x3f ; /* we hard impose 6 bit hashing == 64 cases */
  if (!fil)
    return 0 ;
  nn = n2 = nfound = nfound2 = nsaturated = 0 ;

  for (type = 3 ; type < 6 ; type += 2)
    for (c1=0;c1<256;c1++)
      {
	if (c1 != ace_lower(c1) ||
	    !isIndexed(c1,type))
	  continue ;
	
	for (c2=0;c2<256;c2++)
	  {
	    if (c2 != ace_lower(c2) ||
		!isIndexed(c2,type))
	      continue ;
	    
	    for (c3=0;c3<256;c3++)
	      {
		if (c3 != ace_lower(c3) ||
		    !isIndexed(c3,type))
		  continue ;
				
		switch (type)
		  {
		  case 3:
		    text[0] = '*' ;
		    text[1] = text2[0] = c1 ;
		    text[2] = text2[1] = c2 ;
		    text[3] = text2[2] = c3 ;
		    text[4] = '*' ;
		    text[5] = 0 ;
		    break ;
		  case 5:
		    text[0] = '*' ;
		    text[1] = text2[0] = c1 ;
		    text[2] = '?' ;
		    text[3] = text2[1] = c2 ;
		    text[4] = '?' ;
		    text[5] = text2[2] = c3 ;
		    text[6] = '*' ;
		    text[7] = 0 ;
		    break ;
		  }
				
		if (mask != OWQHash(text2)) 
		  continue ;

		if (0 && strcmp(text,"*pro*")) continue ;
		fprintf (fil ? fil : stderr, "// txt=%s\tm=%d\n", text,  OWQHash(text2)) ;

		if (dictFind (dict, text, 0))
		  {
		    nfound++ ;
		    continue ;
		  }
		found = 0 ;
		if (1) /* 0 means do not test, just recompute all */
		  ks2 = owqExt3 (text, 2, &found, type, debug, fil) ;
		dictAdd (dict, text, 0) ;
		switch (found)
		  { 
		  case 2: /* success */
		    nfound2 += (ks2 ? keySetMax (ks2) : 0) ;
		    keySetDestroy (ks2) ;
		    nfound++ ;
		    continue ;
		  case 1: /* saturation */
		    nsaturated++ ;
		    keySetDestroy (ks2) ;
		    nfound++ ;
		    continue ;
		  case 0: /* go on and compute */
		    break ;
		  }
		nn++ ;
		ks2 = owqExt3 (text, 2, &found, type, debug, fil) ;
		if (ks2)
		  {
		    n2 += keySetMax (ks2) ;
		    keySetDestroy (ks2) ;
		  }
	      }
	  }
      }
      
  messout ("// PreCompute mask=%d, %d (->%d keys)  new triplets, %d (->%d keys) were known before, %d saturated(>%d)"
	   ,mask , nn, n2, nfound, nfound2, nsaturated, OWQMAX) ;
  return nn ;
} /* owqPrecompute */

/***************************************************************/
/* keep classes relevant to this search */
static void owqClassFilterCandidate (char *word, KEYSET ks, BOOL isExtended)
{
  static char classFast [256], classLong [256], isFirst = 1 ;
  KEY key ;
  int ii, jj ;
  char *starWord = 0 ;

  if (isFirst)
    {
      OWB *owb, *owb1 ;
      KEY dummy ;

      isFirst = 0 ;
      ii = 256 ;
      while (ii--) classFast [ii] = classLong [ii] = 0 ;
      
      /* remove from the menu the classes unknown in this database */
      for (owb = owb1 = owbMenu ;  ; owb++)
	{
	  if (owb->cl && ! lexword2key (owb->cl, &dummy, _VClass))
	    continue ;
	  if (owb1 < owb) *owb1 = *owb ;
	  owb1++ ;
	  if (!owb->cl) break ;
	}

      for (owb = owbMenu ; owb->cl ; owb++)
	if ((ii = pickWord2Class (owb->cl)))
	  {
	    if (owb->flag & 0x1) classFast [ii] = 1 ;
	    if ((owb->extend & 0x3) == 0x3) classLong [ii] = 1 ;
	  }
      if ((ii = pickWord2Class ("Text")))
	classLong [ii] = 1 ;
      if ((ii = pickWord2Class ("LongText")))
	classLong [ii] = 2 ;
      if ((ii = pickWord2Class ("FicheView")))
	classLong [ii] = 2 ;
    }
  
  starWord = word ? messprintf ("*%s*", word) : 0 ;
  for (ii = jj = 0 ; ii < keySetMax (ks) ; ii++)
    {
      key = keySet (ks, ii) ;
      if (
	  (!isExtended && classFast[class(key)]) ||
	  (
	   isExtended && 
	   /* classLong[class(key)] && march 27, in case of extended search extend the cass names */
	   (!word || classLong[class(key)] == 2 || pickMatch (name(key), starWord))
	   ) 
	  )
	{
	  if (ii > jj) keySet (ks, jj) = keySet (ks, ii) ;
	  jj++ ;
	}
    }
  keySetMax (ks) = jj ;
} /* owqClassFilterCandidate */

/***************************************************************/
/* find the intersect of ks0 with all the 3 letter index present in "text" */
static KEYSET owqCandidate (char *word, KEYSET ks0, BOOL isExtended, BOOL *isSaturatedp, BOOL debug)
{
  char *cp, text[8] ;
  KEYSET ks1, ks2 ;
  KEYSET ks = ks0 ? keySetCopy (ks0) : 0 ;
  int found, type ;
  int create = 0 ;
  int ii, jj ;
  static BitSet bb = 0 ;
  BOOL needFilter = FALSE ;

  memset (text, 0, sizeof(text)) ;

  if (ks) /* may be filtering will kill ksInit */
    {
      owqClassFilterCandidate (word, ks, isExtended) ;
      if (!ks)
	return 0 ;
    }

  *isSaturatedp = TRUE ;
  bb = bitSetReCreate (bb, 256) ;
  /* if nameSearch == 1, we do not create the 3-index */
  for (create = 0 ; create < 2 ; create++)
    {	  
      if (create && !isExtended)  /* should we create ? */
	break ;   /* only in extended case when we search for Text etc */

      for (type = 3 ; type < 6 ; type += 2)
	for (ii = 0, cp = word ; *cp ; ii++, cp++)
	  {
	    switch (type)
	      {
	      case 3:
		if (strlen(cp) < type ||
		    !isIndexed(*cp,type) ||
		    !isIndexed(*(cp+1),type) ||
		    !isIndexed(*(cp+2),type))
		  continue ;
		text[0] = '*' ;
		text[1] = *cp ;
		text[2] = *(cp+1) ;
		text[3] = *(cp+2) ;
		text[4] = '*' ;
		text[5] = 0 ;
		break ;
	      case 5:
		if (strlen(cp) < type ||
		    !isIndexed(*cp,type) ||
		    !isIndexed(*(cp+1),3) ||
		    !isIndexed(*(cp+2),type) ||
		    !isIndexed(*(cp+3),3) ||
		    !isIndexed(*(cp+4),type))
		  continue ;
		text[0] = '*' ;
		text[1] = *cp ;
		text[2] = '?' ;
		text[3] =  *(cp+2) ;
		text[4] = '?' ;
		text[5] =  *(cp+4) ;
		text[6] = '*' ;
		text[7] = 0 ;
	      }
	  
	    /* bb is used to do a single call although we loop on create */
	    jj = 2*ii + (type==5?1:0) ;
	    if (bitt (bb, jj)) continue ;
	    
	    if (ks)
	      {
		if (needFilter)
		  owqClassFilterCandidate (create ? word : 0, ks, isExtended) ;
		needFilter = FALSE ;
		if (create &&  /* should we create ? */
		    (
		     type == 5 ||
		     keySetMax(ks) < 300  /* if small ks, it is useless, our ks is a bit too large 
					   * but we no longer care */
		    ))
		  goto basta ;
	      }
	    ks2 = owqExt3 (text, create, &found, type, debug, 0) ;
	    switch (found)
	      {
	      case 0: /* missing */
		continue ;
	      case 1: /* saturated */ 
		bitSet (bb, jj) ;/* mark but ignore saturated results */
		keySetDestroy (ks2) ;
		continue ;
	      case 2: /* success */
		bitSet (bb, jj) ;
		break ;
	      }
	    *isSaturatedp = FALSE ;
	    if (ks && ks2)
	      {
		ks1 = ks ;
		ks = keySetAND (ks1, ks2) ;
		keySetDestroy (ks1) ; keySetDestroy (ks2) ;
	      }
	    else  /* first time we initialise ks */
	      { ks = ks2 ; ks2 = 0 ; needFilter = TRUE ; }
	    if (!ks)
	      return 0 ;
	  }
    }
  
 basta:
  if (ks && needFilter)
    owqClassFilterCandidate (create ? word : 0, ks, isExtended) ;
  if (ks && !keySetMax(ks))
    keySetDestroy (ks) ;

  return ks ;
} /* owqCandidate */

/***************************************************************/

/* loop on all sentences, search for every long word */
static KEYSET  owqExtMultiWordCandidate (char *qName, Stack s, BOOL *isSaturatedp, BOOL debug)
{  
  char *cp, *cq, *word ;
  KEYSET ksCandidate = 0, ks1 = 0, ks2 = 0 ;
  BOOL isSaturated  ;    /* each word */

  *isSaturatedp = TRUE ; /* global for the whole sentence */
  word = messalloc (strlen(qName)+ 1) ;
  for (cp = qName, cq = word ;  ; cp++)
    {
      *cq = *cp ;
      if (!*cp || *cp == ' ')
	{
	  *cq = 0 ;
	  if (strlen(word) >= 3)
	    {
	      ks2 = owqCandidate (word,  ksCandidate, /*isExtended=*/TRUE, &isSaturated, debug) ;
	      if (!isSaturated)
		{
		  *isSaturatedp = FALSE ;
		  if (ksCandidate && ks2)
		    {
		      ks1 = ksCandidate ;
		      ksCandidate = keySetAND (ks1, ks2) ;
		      keySetDestroy (ks1) ; keySetDestroy (ks2) ;
		    }
		  else  /* first time we initialise ksCandidate */
		    { ksCandidate = ks2 ; ks2 = 0 ; } 
		  if (!ksCandidate || !keySetMax(ksCandidate))
		    {
		      keySetDestroy (ksCandidate) ;
		      break ;
		    }
		}
	    }
	  pushText (s, messprintf("*%s*",word)) ;
	  cq = word ;
	}
      else 
	cq++ ;
      if (!*cp)
	break ;
    }
  messfree (word) ;
  return ksCandidate ;
} /* owqExtMultiWordCandidate */

/***************************************************************/
/* filters for the existence of 1 and 2 letters things */
static void owqExtMultiWordFilter (KEYSET ks, Stack s)
{
  int ii, jj ;
  char *cp ;
  KEY key ;
  
  for (ii = jj = 0 ; ii < keySetMax(ks) ; ii++)
    {
      key = keySet (ks, ii) ; 
      stackCursor(s, 0) ; 
      while ((cp = stackNextText(s)))
	if (!pickMatch (name(key), messprintf("*%s*", cp)))
	  break ;
      if (!cp) /* end of stack reached: success */
	{
	  if (ii > jj) keySet (ks, jj) = key ;
	  jj++ ;
	}
    }
  keySetMax (ks) = jj ;
  return ;
} /* owqExtMultiWordFilter */

/***************************************************************/
/* loop on all sentences, search for every long word */
static KEYSET  owqExtended (char *qName, BOOL *isSaturatedp)
{  
  char *cp ;
  KEYSET ks = 0, ks0, ks1, ks2, ks3, ks4 ;
  BOOL debug = TRUE ;
  OWB *owb ;
  int i, j, x3, x5, n3, n5 ;
  BOOL isError = FALSE ;
  Stack s = stackCreate (128) ;
  char timeBuf[25] ;

  *isSaturatedp = FALSE ;
  for (x3 = x5 = n3 = n5 = 0, cp = qName ; !n5 && *cp ; cp++)
    {
      if (isIndexed(*cp,3)) 
	x3++ ;
      else
	x3 = 0 ;
      if (isIndexed(*cp,5)) 
	x5++ ;
      else
	x5 = 0 ;
      if (x3 >= 3) n3++ ;
      if (x5 >= 5) n5++ ;
    }

  if (n3 > 0) /* at least one word of 3 meaningful char */
    ks = owqExtMultiWordCandidate (qName, s, isSaturatedp, debug) ;

  if (ks  && arrayMax(ks) <= 25000)
    {
      ks0 = query (ks, "NOT CLASS LongText") ;
      ks1 = query (ks, "CLASS LongText") ;
      if (! n5 && keySetMax (ks1) > 3000)
	keySetMax (ks1) = 1 ;  /* do not search long text for short frequent words */
      if (debug) printf("// %s owqExtMulti %d candidate %d longText "
			, timeShow (timeNow(), timeBuf, 25), keySetMax(ks), keySetMax(ks1)) ;
      keySetDestroy (ks) ;

      owqExtMultiWordFilter (ks0, s) ; /* verify name(key) for short words */
      longTextFilter (ks1, s, TRUE) ;  /* remap the abstracts into papers */

      ks = query (ks0, "CLASS Text") ;
      ks2 = query (ks, ">quoted_in") ;
      ks3 = query (ks0, "! CLASS Text") ;
      keySetDestroy (ks) ;
      keySetDestroy (ks0) ;

      ks4 =  keySetOR (ks1, ks2) ;
      ks0 =  keySetOR (ks3, ks4) ;
      keySetDestroy (ks1) ; keySetDestroy (ks2) ;
      keySetDestroy (ks3) ; keySetDestroy (ks4) ;

      if (debug) printf("// %d to be filtered", keySetMax(ks0)) ;
      
      ks = keySetCreate () ;
      for (j = 0, owb = owbMenu ; owb->cl ; owb++)
	{
	  ks1 = owqFilter (ks0, owb) ; 
	  for (i = 0 ; ks1 && i < keySetMax (ks1) ; i++)
	    if (class (keySet (ks1, i)) == _VGene) /* some queries return clones */
	      keySet (ks, j++) = keySet (ks1, i) ;
	  keySetDestroy (ks1) ;
	}
      if (debug) printf("// %s ---> %d found ", timeShow (timeNow(), timeBuf, 25), keySetMax (ks)) ;
      keySetDestroy (ks0) ;
      keySetSort (ks) ;
      keySetCompress (ks) ;
      if (debug) printf(" ==> %d compressed %s\n", keySetMax (ks), qName) ;
    }
  else if (ks)
    *isSaturatedp = TRUE ;

  if (!ks || !keySetMax(ks) || *isSaturatedp)
    {
      KEY key ;
      char *cp ;
      
      if (*isSaturatedp)
	{ 
	  cp =  "Sorry, too many hits. Please try a longer or more specific sentence" ;
	}
      else if (!n3)
	cp = "Sorry, this word is too short, enter at least 3 letters or digits" ;
      else
	cp = "Sorry, locus not found." ;
	
      keySetDestroy (ks) ;
      lexaddkey (cp, &key, _VKeySet) ;
      ks = keySetCreate () ;
      keySet (ks, 0) = key ;
      arrayStore (key, ks, "k") ; /* needed in case we destroyed all keyset
				     which sets empty flag and impeaches list */
      isError = TRUE ; /* do not can the result */
    }

  stackDestroy (s) ;
  
  /* store general result */
  if (!isError)
    owqStoreCannedList (ks, qName, 1, 0, 0) ;

  return ks ;
} /* owqExtended */

/***************************************************************/

static char* owqCleanQuery (char *cl, char *qName, int *nwordp, int *ncharp)
{
  char *cp, *cq = qName - 1, *cq2 ; 
  int i = 0, nw, nchar, isSpace, isStar ;
  KEY key ;

  while (*++cq) 
    switch (*cq)
      {
      case '%':  /* web way to protect spaces etc */
	if (*(cq+1) >= '0' && *(cq+1)<='9'&& *(cq+2) >= '0' && *(cq+2)<='9')
	  { *cq = ' ' ; cq2 = cq + 1 ; while (*cq2) { *cq2 = *(cq2+2) ; if (*cq2) cq2++ ; } }
	else if (*(cq + 1) == '%') /* acedb hates % in queries, replace it by a joker */
	  { *cq = '?' ; cq2 = cq + 1 ; while (*cq2) { *cq2 = *(cq2+1) ; if (*cq2) cq2++ ; } }
	break ;
      case '*':
      case '?':
      case ' ':
	break ;
      case '\n':
      case '\t':
	*cq = ' ' ;
	break ;
      default:
	i++ ;
	*cq = ace_lower (*cq) ; 
	break ;
      }
  if ((i < 1) )
    *qName = 0 ;
  
  cp = qName + strlen (qName) ;
  
  *cp++ = ' ' ;  *cp = 0 ; /* add a space */
  cp = qName - 1; nw = 0 ; isSpace = 1 ; isStar = 0 ;
  while (*++cp)
    switch (*cp)
      {
      case '?':
      case '*':
	isStar = 1 ;
	break ;
      case ',': *cp=' ';
      case ';': *cp=' ';
      case '(':
      case ')':
      case ' ': isSpace = 1 ; 
	break ;
      default:
	/* remove silly words which are always saturated */
	if (isSpace)
	  {
	    if (!strncasecmp (cp, "caenorhabditis ", 15))
	      memset (cp, ' ', 14) ;
	    else if (!strncasecmp (cp, "elegans ", 8))
	      memset (cp, ' ', 7) ;
	    else if (!strncasecmp (cp, "homo ", 5))
	      memset (cp, ' ', 4) ;
	    else if (!strncasecmp (cp, "sapiens ", 8))
	      memset (cp, ' ', 7) ;
	    else if (!strncasecmp (cp, "gene ", 5))
	      memset (cp, ' ', 4) ;
	    else if (!strncasecmp (cp, "and ", 4))
	      memset (cp, ' ', 3) ;
	    else if (!strncasecmp (cp, "or ", 3))
	      memset (cp, ' ', 2) ;
	    else if (!strncasecmp (cp, "xor ", 4))
	      memset (cp, ' ', 3) ;
	    else
	      nw++ ;
	  }
	isSpace = 0 ;
	break ;
      }

  /* clean up leading, multiple middle spaces  and trailing spaces */
  cp = cq = qName ;
  while (*cp == ' ') cp++ ; /* clean leading spaces */
  while (*cp)
    {
      if (*cp == ' ' && cq > qName && *(cq-1) == ' ') ; /* avoid multiple spaces */
      else { if (cq != cp) *cq = *cp ; cq++ ; }
      cp++ ;
    }
  *cq = 0 ;
  /* clean up multiple * */
  cp = cq = qName ;
  while (*cp)
    {
      if (*cp == '*' && cq > qName && *(cq-1) == '*') ; /* avoid multiple spaces */
      else { if (cq != cp) *cq = *cp ; cq++ ; }
      cp++ ;
    }
  *cq = 0 ;

 /* clean terminal space and dots */
  while (--cq >= qName  && 
	 (*cq == ' ' ||
	  *cq == '.' ||
	  *cq == '_'))
    *cq = 0 ;

  /* horrible hack because of gene famillies like actin actin.1 actin.2 */
  if (!cl && !isStar && lexword2key (cp, &key, _VGene) &&
      lexword2key (messprintf("%s.1",cp), &key, _VGene))
    {
      cp = qName + strlen (qName) ;
      *cp++ = '*' ;  *cp = 0 ; /* add a star to get the whole familly */
    }

  for (nchar = 0, cp = qName ; *cp ; cp++)
    {
      if (isIndexed(*cp,3)) 
	nchar++ ;
    }

  *nwordp = nw ;
  *ncharp = nchar ;
  return qName ;
} /* owqCleanQuery  */

/***************************************************************/

static void owqPrecomputeCommand (void)
{
  FILE *fil =  0;
  DICT *dict = dictCreate (10000) ;
  int mask = 0 ;
  char *cp ;

  freeint (&mask) ;

  if (mask == -1)
    {
      if (sessionGainWriteAccess())
	{
	  KEYSET ks, ks1 ;
	  
	  freeOutf ("// owqPrecomputeCommand -1: gene2newname\n") ;
	  ks = query (0,"Find Gene") ;
	  ks1 = OWQGene2NewName (ks, TRUE) ;
	  keySetDestroy (ks) ;      
	  keySetDestroy (ks1) ;
	  
	  freeOutf ("// owqPrecomputeCommand -1: product2gene\n") ;
	  ks = query (0,"Find Product") ;
	  ks1 = OWQProduct2Gene (ks, TRUE) ;
	  keySetDestroy (ks) ;      
	  keySetDestroy (ks1) ;
	}
      else
	freeOutf ("// owqPrecomputeCommand failed, no write access\n") ;
    }
  else
    {
      if ((cp = freeword ())) 
	fil = filopen (cp, 0, "w") ;
      if (mask == 0)
	{
	  int i ;
	  for (i = 1 ; i <= 64 ; i++)
	    owqPrecompute (i, fil, dict) ;
	}
      else
	owqPrecompute (mask, fil, dict) ;
      if (fil) filclose (fil) ;
    }
  dictDestroy (dict) ;
  return ;
} /* owqPrecomputeCommand */

/***************************************************************/

static KEYSET owqDoQuery (KEYSET ks0, char *qName, int nword, int nchar, BOOL *isSaturatedp)
{
  KEYSET ks = 0, ksCandidate = 0 ;
  BOOL debug = FALSE ;
  int found ;

  /* if single word first search as a direct name */
  if (nword == 1) 
    {
      ks = owqGetCannedList (qName, &found, 1) ;
      switch (found)
	{
	case 0: /* search needed */
	  break ;
	case 1: /* saturated */
	  *isSaturatedp = TRUE ;
	case 2:	      
	  goto done ;
	}
      ksCandidate = owqCandidate (qName, 0, 0, isSaturatedp, debug) ;
      ks = owqOneWord (0, 0, qName, ksCandidate, nchar < 3 ? TRUE : FALSE, debug) ;
      keySetDestroy (ksCandidate) ;
    }
  
  if (!ks) /* now search extended */
    { 
      ks = owqGetCannedList (qName, &found, 1) ; 
      switch (found)
	{
	case 0: /* search needed */
	  break ;
	case 1: /* saturated */ 
	  *isSaturatedp = TRUE ;
	case 2:	
	  *isSaturatedp = FALSE ;
	  goto done ;
	}
      ks = owqExtended (qName, isSaturatedp) ; 
    }

 done:
  if (ks && ks0 && ! *isSaturatedp)
    {
      KEYSET ks2 = ks ;
      ks = keySetAND (ks0, ks2) ;
      keySetDestroy (ks2) ;
    }
  return ks ;
}  /* owqDoQuery */

/***************************************************************/
/******************* PUBLIC ************************************/
/***************************************************************/

/* to be used from query language */
KEYSET owqQuery (KEYSET ks0, char *qName)
{
  int nword = 0, nchar = 0 ;
  BOOL isSaturatedDummy = FALSE ;

  if (qName)
    owqCleanQuery (0, qName, &nword, &nchar) ;
  
  return owqDoQuery (ks0, qName, nword, nchar, &isSaturatedDummy) ;
} /* owqQuery */

/***************************************************************/

KEYSET cdnaOptimisedWebQuery (KEYSET ksActive)
{
  KEYSET ks = 0 ;
  char *cp, *qName = 0 ;
  char cl[1024], qBuf[4096], qBufExact[4096]  ;
  int nword = 0, nchar ;
  BOOL debug = FALSE ;
  BOOL isCtx = FALSE ;
  BOOL isSaturated = FALSE ;

  cDNAAlignInit () ;

  cl[0] = 0 ;
  while ((cp = freeword ()))
    {
      if (!strcasecmp (cp, "-ctx"))
	isCtx = TRUE ; 
      else if (!strcasecmp (cp, "Precompute"))
	{ strncpy (cl, cp, 1023) ; qName = "" ; break ; }
      else if (!strcasecmp (cp, "Gene2NewName"))
	{ strncpy (cl, cp, 1023) ; qName = "" ; break ; }
      else 
	{
	  if (cp && strcasecmp (cp, "0") && strcasecmp (cp, "extended"))
	    strncpy (cl, cp, 1023) ;

	  cp = freeword () ;   
	  if (0 && /* this hack is stupid, if anything it should go to geneid */
	      cp && !strncmp (cp, "LOC",3)) /* a Hack */
	    {
	      strcpy (cl, "GeneId") ; /* strcpy (cl, "LocusId") ; */
	      cp = cp+3 ;
	    }
	  if (! strcasecmp (cl, "LocusId"))
	    strcpy (cl, "GeneId") ;
	  if (cp && *cp)
	    { 
	      qName = qBuf ; strncpy (qBuf, cp, 4095) ; strncpy (qBufExact, cp, 4095) ;
	      owqCleanQuery (cl, qName, &nword, &nchar) ;
	    }
	  break ;
	}
    }
 
  if (!qName || !strcmp (qName, "*") || !strcasecmp (qName, "disease") || !strcasecmp (qName, "disease genes"))
    ks = isCtx ? keySetCopy (ksActive) : owqDefault () ; 
  else if (!*cl)  /* 0 or extended class */
    ks = owqDoQuery (isCtx ? ksActive : 0, qName, nword, nchar, &isSaturated) ;
  else if (!strcasecmp (cl, "Precompute"))
    owqPrecomputeCommand () ;
  else if (!strcasecmp (cl, "Gene2NewName"))
    ks = OWQGene2NewName (ksActive, FALSE) ;
  else if (!strcasecmp(cl,"chrom"))
    ks = owqChromSearch (qName) ;
  else   /* user asks for a particular class/qName */
    {
      ks = owqOneWord (cl, qBufExact, qName, 0, TRUE, debug) ;
      /* june 9, 2005, we get hits from tatiana using the newname
       * because the genetic name is too recent 
       */
      if ((!ks || !keySetMax (ks)) && !strcasecmp (cl, "Gene"))
	{
	  keySetDestroy (ks) ;
	  ks = owqOneWord ("NewName", qBufExact, qName, 0, TRUE, debug) ;
	}
    }
  return ks ;
} /* cdnaOptimisedWebQuery */

/***************************************************************/
/***************************************************************/
/***************************************************************/
