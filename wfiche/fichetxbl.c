#include "../wfiche/taxtree.h"
#include "../wfiche/biolog.h"
#include "../wfiche/gtitle.h"
#include "../wac/ac.h"
#include "vtxt.h"


/*******************************************************/
/*******************************************************/
/*******************************************************/
/* TAXBLAST output analysis */

typedef struct taxstruct { int level, elevel, hits, cumul, mark, bestHit ; BOOL isSpecies ; BOOL isPrint ; } TAXTREE ;

static void showTT (Array a, Stack s)
{
  TAXTREE *tp ;
  int i ;

  if (arrayExists (a))
    for (i = 0 ; i < arrayMax (a) ; i++)
      {
	tp = arrp (a, i, TAXTREE) ;
	printf ("%3d level=%d\thits=%d\tcumul=%d\tname=%s\n"
		, i, tp->level, tp->hits, tp->cumul, stackText (s, tp->mark)) ;
      }
}

static void ficheTAXTreeFillArray (AC_OBJ oProduct, Array a, Stack s)
{
  int i, ii, level ;
  char *cp, *cq, buf[1000], buf2[1000] ;
  const char *ccp ;
  AC_TABLE nTag ;
  TAXTREE *tp ;
  AC_HANDLE h = handleCreate () ;

  level = freesettext (myTaxTree, 0) ; ii = 0 ;
  while (freecard (level))
    {
      cp = freepos () ;
      if (!*cp) /* this happens becuase mytaxTree ends on \n\000
		   so we get a card for the last empty string \000
		   but we do not want to change mytaxTree becuase it
		   is also used by kantormegaparse
		*/
	continue ;
      tp = arrayp (a, ii++, TAXTREE) ;
      /* count dots to find the level */
      i = 0 ;
      while (*cp++ == '.') i++ ;
      tp->level = i/2 ;

      /* split the tag name */
      buf[0] = 'N' ;
      cp-- ; cq = buf + 1 ;
      while (*cp != '.') *cq++ = *cp++ ;
      *(cq - 1) = 0 ; /* remove terminal space */
      tp->mark = pushText (s, buf + 1) ;
      
      if (strlen(buf) < 2)
	{
	  arrayMax (a)-- ;
	  continue ;
	}
      /* species are flagged '-1' */
      if (strstr (cp, " -1"))
	tp->isSpecies = TRUE ; 

      /* now look for the tag in this protein */
      cq = buf ;
      while (*++cq) if (*cq == ' ') *cq = '_' ;
      if ((nTag = ac_tag_table (oProduct, buf, h)) &&
	  ac_table_int (nTag, 0, 0, -1) != -1)
	{
	  tp->hits = ac_table_int (nTag, 0, 0, -1) ;
	  if ((ccp=ac_table_printable (nTag, 0, 2, 0)))
	    {
	      /* looks like:  sp|P54706|COFI_DICDI */
	      ccp-- ;
	      while (*++ccp) 
		if (*ccp == '|')
		  { ccp++ ; break ; }
	      strcpy (buf2, ccp) ;
	      cp = buf2 ; cp-- ;
	      while (*++cp)
		if (*cp == '|')
		  { *cp = 0 ; break ; }
	      if (strlen (buf2))
		tp->bestHit =  pushText (s, buf2) ;
	    }
	}
    }
  ac_free (h) ;
} /* ficheTAXTreeFillArray  */

/*******************************************************/

static void ficheTAXTreeCumulate (Array a, Stack s)
{
  TAXTREE *tp, *tq ;
  int ii, j ;
  
  /* cumulate at the end each node contains the sum of the lower + his hits */
  for (ii = arrayMax (a) - 1, tp = arrp (a, ii, TAXTREE) ; ii >= 0 ; tp--, ii--)
    {
      tp->cumul += tp->hits ;
      if (tp->cumul)
	for (j = ii - 1, tq = tp -1 ; j >= 0 ; j--, tq--)
	  if (tq->level < tp->level)
	    {
	      tq->cumul += tp->cumul ; 
	      break ;
	    }
    }
}/* ficheTAXTreeCumulate */

/*******************************************************/
static void spacer (vTXT vtxt, GMP *gmp, int n)
{
  char buf[1000] ;
 
  vtxtBreak (vtxt) ;
	  
  if (n > 0)
    {
      memset (buf, '-', 2*n ) ;
      buf [2 * n] = 0 ;
      vtxtPrintf (vtxt, buf) ;
    }
}

static void ficheTAXTreeShowArray (vTXT vtxt, GMP *gmp, Array a, Stack s)
{
  TAXTREE *tp, *tq ;
  int ii, j ;
  
  for (ii = 0 ; ii < arrayMax (a) ; ii++)
    {
      tp = arrp (a, ii, TAXTREE) ;

      tp->isPrint = FALSE ; /* we print if */
      if (!tp->cumul) /* i have a hit */
	goto others ;
      if (tp->level == 0 || tp->isSpecies || tp->hits)
	tp->isPrint = TRUE ; /* keep level zero */
      for (tq = tp + 1, j = ii+1 ; !tp->isPrint && j < arrayMax (a) ; tq++, j++)
	{           /* or i have a rich brother */
	  if (tq->cumul)
	    {
	      if (tq->level == tp->level)
		tp->isPrint = TRUE ; /* found brother */
	      if (tq->level <= tp->level)
		break ; /* back to parent */
	    }
	}
      for (tq = tp - 1, j = ii-1 ; !tp->isPrint && j >= 0 ; tq--, j--)
	{           /* or i have a brother */
	  if (tq->cumul)
	    {
	      if (tq->level == tp->level)
		tp->isPrint = TRUE ; /* found brother */
	      if (tq->level <= tp->level)
		break ; /* back to parent */
	    }
	}
      /* compute the effective level
       * .e. do not put spacers for skipped lines
       */

      if (tp->level)
	for (tq = tp - 1, j = ii-1 ; j >= 0 ; tq--, j--)
	  if (tq->isPrint)
	    {
	      if (tq->level == tp->level)
		{ tp->elevel = tq->elevel ; break ; }
	      else if (tq->level < tp->level)
		{ tp->elevel = tq->elevel + 1 ; break ; }
	    }

      if (!tp->isPrint)  /* used to show others, so we need its elevel */
	goto others ;

      /* spacer */
      spacer (vtxt, gmp, tp->elevel) ;
     /* data */
      if (tp->isSpecies)
	vtxtItalic (vtxt, stackText (s, tp->mark)) ;
      else
	vtxtPrintf (vtxt, "%s", stackText (s, tp->mark)) ;
      vtxtBold (vtxt, messprintf (" %d hit%s", tp->cumul, tp->cumul > 1 ? "s" : "")) ;
      if (tp->bestHit)
	{
	  vtxtPrintf (vtxt, ", best hit: ") ;
	  gmpURL (vtxt, gmp
		  , messprintf ("https://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=%s"
				, stackText (s, tp->bestHit))
		  , stackText (s, tp->bestHit)) ;
	}
      /* other of same kind */
    others:
      if (tp->level > (tp+1)->level)
	{
	  int n = 0 ;
	  TAXTREE *tp2 = tp ;

	  for (j = ii , tq = tp ; j >= 0 ; j--, tq--)
	    {
	      if (tq->isPrint && tq->level == tp2->level)
		{ n += tq->hits ; continue ; }
	      if (tq->level < (tp+1)->level)
		break ; /* tq section is not finished */
	      if (tq->isPrint && tq->level < tp2->level)
		{
		  n = tq->hits ;
		  if (n > 0 && tq->mark)
		    {
		      spacer (vtxt, gmp, tq->elevel + 1) ;
		      
		      vtxtPrintf (vtxt, "Other %s %d hit%s"
				  , stackText (s, tq->mark)
				  , n , n>1 ? "s" : ""
				  ) ;
		      tq->mark = 0 ;
		    }
		  break ;
		}
	      if (tq->level < tp2->level)
		tp2 = tq ;
	    }
	}
    }	

  /* now we recursivelly close the last line */
  tp = arrp (a, arrayMax (a) - 1 , TAXTREE) ;

  for (j = arrayMax (a) - 2 ; j >= 0 ; j--)
    {
      tq = arrp (a, j, TAXTREE) ;
      if (tq->level >= tp->level || 
	  tp->mark == tq->mark)
	continue ;
      if (tq->hits && tq->mark)
	{
	  spacer (vtxt, gmp, tp->elevel) ;
	  vtxtPrintf (vtxt, "Other %s %d hit%s"
		      , stackText (s, tq->mark)
		      , tq->hits , tq->hits > 1 ? "s" : ""
		      ) ; 
	  tq->mark = 0 ;
	}
      tp = tq ;
    }

  /* now report the remaining others */
}/* ficheTAXTreeShowArray */
 
/*******************************************************/

static char *ficheTAXTreeShowClosestAncestor (Array a, Stack s)
{
  TAXTREE *tp ;
  int n0, ii, cmax=0 ;

  ficheTAXTreeCumulate (a, s) ;

  /* first treat level 0 ancestors */
  for (n0 = ii = 0, tp = arrp (a, 0, TAXTREE) ; ii < arrayMax (a) ; ii++, tp++)
    if (tp->level == 0 && tp->cumul > 0) n0++ ;
  if (n0 > 1)
    {
      vTXT vtxt = vtxtCreate () ;
      char *cp ;

      for (n0 = ii = 0, tp = arrp (a, 0, TAXTREE) ; ii < arrayMax (a) ; ii++, tp++)
	if (tp->level == 0 && tp->cumul > 0) 
	  { 
	    if (n0++) vtxtPrintf (vtxt, " and") ; 
	    vtxtPrintf (vtxt, " %s", stackText (s, tp->mark)) ;
	  }
      cp = strnew (vtxtPtr (vtxt), 0) ;
      vtxtDestroy (vtxt) ;
      return cp ;
    }

   
   /* closest ancestor is the lowest guy richer than his  sons */
   /* find max */
   for (ii = 0, tp = arrp (a, 0, TAXTREE) ; ii < arrayMax (a) ; ii++, tp++)
    if (tp->cumul > cmax) cmax = tp->cumul ;
  /* find lowest as rich as max */
  if (cmax)
    for (ii = arrayMax (a) - 1, tp = arrp (a, ii, TAXTREE) ; ii >= 0 ; ii--, tp--)
      if (tp->cumul == cmax)
	return strnew (stackText (s, tp->mark), 0) ;

  return 0 ;
}/* ficheTAXTreeShowClosestAncestor */

/*******************************************************/
/* just count the hits */ 
static int ficheTAXTreeCountHits (GMP *gmp)
{ 
  int ir, nHits = 0 ;
  AC_OBJ oTmp = gmp->kantor ? gmp->kantor : gmp->product ;
  AC_TABLE gTag = oTmp ? ac_tag_table (oTmp, "Tax_count", 0) : 0 ;

  for (ir = 0 ; gTag && ir < gTag->rows ; ir++)
    nHits += ac_table_int (gTag, ir, 1, 0) ;

  ac_free (gTag) ;
  return nHits ;
} /* ficheTAXTreeCountHits */

/*******************************************************/
/* return a developed marked pruned tree */	
static char *ficheTAXTreeStatement (vTXT vtxt, GMP *gmp)
{
  AC_OBJ oTmp = gmp->kantor ? gmp->kantor : gmp->product ;
  BOOL debug = FALSE ;

  if (myTaxTree && ac_has_tag (oTmp, "Tax_count"))
    {
      Stack s = stackCreate (4000) ;
      Array a = arrayCreate (64, TAXTREE) ; 

      ficheTAXTreeFillArray (oTmp, a, s) ;

      ficheTAXTreeCumulate (a, s) ;
      if (debug) showTT (a, s) ;
      ficheTAXTreeShowArray (vtxt, gmp, a, s) ;
      
      arrayDestroy (a) ;
      stackDestroy (s) ;
    }    

  return vtxtPtr (vtxt) ;
}  /* ficheTAXTreeStatement  */

/*******************************************************/
/* returns the name of the closest common ancestor */
char *ficheTAXClosestAncestorStatement (GMP *gmp)
{
  char *ptr = 0 ;
  AC_OBJ oTmp = gmp->kantor ? gmp->kantor : gmp->product ;

  if (ac_has_tag (oTmp, "Tax_count"))
    {
      Stack s = stackCreate (4000) ;
      Array a = arrayCreate (64, TAXTREE) ; 
      
      ficheTAXTreeFillArray (oTmp, a, s) ;
      ptr = ficheTAXTreeShowClosestAncestor (a, s) ;

      arrayDestroy (a) ;
      stackDestroy (s) ;
    }

  return ptr ;
} /* ficheTAXClosestAncestorStatement */

/*******************************************************/
/*******************************************************/
/* -===================================- /
 * -=  TAX TREE                  =- /
 * -===================================- */	
char *fichePRODUCTBlastPParagraphContent (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE gBlastp = 0 ; 
  char	dateBuffer[256] ;
  int  allHits = 0 ; 
  
  if (!ac_has_tag (gmp->kantor, "Blastp_Date") &&  !ac_has_tag (gmp->kantor, "Taxblast_Date"))
    return 0 ;
  
  timeShowFormat (ac_tag_date (gmp->kantor, "Blastp_Date", timeNow ()), "%b %d, %Y", 
		 dateBuffer, sizeof (dateBuffer)) ; 
  
  vtxtPrintf (blkp, "BlastP analysis was run at ") ; 
  gmpURL (blkp, gmp, "https://www.ncbi.nlm.nih.gov/BLAST", "NCBI") ; 
  vtxtPrintf (blkp, " against the nr database on %s" , dateBuffer) ; 
  gBlastp = ac_tag_table (gmp->product, "BlastP", h) ;
  if ((allHits = ficheTAXTreeCountHits (gmp)))
    vtxtPrintf (blkp, ", and identified %d hit%s with E value better than 0.001"
		, allHits, _multi (allHits)) ; 
  else if (gBlastp &&
	   (allHits = gBlastp->rows))
    vtxtPrintf (blkp, ", and identified %d hit%s"
		, allHits, _multi (allHits)) ; 
  else
    vtxtPrintf (blkp, ", yielding no significant hit at threshold 0.001") ;
  if (gBlastp)
    {
      vtxtPrint (blkp, ". The best hits are") ;
      vtxtBreak (blkp) ;
      allHits = gBlastp->rows ;
      if (allHits > 1 && gmp->markup)
	{  
	  if (allHits > 60)
	    {
	      char *cp = hprintf (h, "<a href=\"javascript:openAceViewAction ('Product', '%s', 'blastp')\"><i>more</i></a><br>\n", ac_name (gmp->product)) ; 
	      ficheProductBlastPTableContent (blkp, gmp, 18, cp) ;

	    }
	  else
	    ficheProductBlastPTableContent (blkp, gmp, 0, 0) ;
	}
    }
  ac_free (h) ;
  return vtxtPtr (blkp) ;
}

char *fichePRODUCTTaxTreeParagraphContent (vTXT blkp, GMP *gmp)
{
  vTXT buf = vtxtCreate () ;
  char *ptr ;
  
  if (gmp->markup) 
    vtxtMarkup (buf) ;

  if (!gmp->markup)
    fichePRODUCTBlastPParagraphContent (blkp, gmp) ;

  if ((ptr = ficheTAXTreeStatement (buf, gmp)))
    {
      vtxtDot (blkp) ;
      vtxtPrintf (blkp, "Based on TaxBlast") ;
      vtxtPrintf (blkp, " (BlastP E < .001) related proteins "
		  "are found in the following organism(s): ") ;
      vtxtBreak (blkp) ;
      vtxtPrintf (blkp, "%s", ptr) ;
    } 

  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
} /* fichePRODUCTTaxTreeParagraphContent */

/***************************************************************************/

void fichePRODUCTTaxTreeParagraph (vTXT blkp, GMP *gmp)
{
  vTXT buf = vtxtCreate () ;
  if (gmp->markup)
    vtxtMarkup (buf) ;

  ficheNewAceKogStatement (buf, gmp, FALSE) ;
  vtxtBreak (buf) ;
  fichePRODUCTTaxTreeParagraphContent (buf, gmp) ; 
  if (vtxtPtr (buf))
    {
      gmpSection (blkp, gmp, "Taxonomy", "Closest homologues in other species") ;
      vtxtPrintf (blkp, "%s", vtxtPtr (buf)) ; 
    }
  vtxtDestroy (buf) ;
} /* fichePRODUCTTaxTreeParagraph */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  CALLS FROM ACEDB
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/* -===================================- /
* -=  COMMON ANCESTORS                =- /
* -===================================- */	

char *fogicArguableCommonAncestors (const char *namKantor)
{
  int  spc ; 
  AC_DB db ; 
  AC_OBJ oKantor = 0 ;
  char *whatToReturn=0 ; 
  GMP *gmp = (GMP *) messalloc (sizeof (GMP)) ; 
  
  if ((db = ac_open_db ("local", 0)) && 
      (oKantor = ac_get_obj (db, "kantor", namKantor, 0)) )
    {
      if (ac_has_tag (oKantor, "Tax_tree"))
	{
	  gmp->kantor = oKantor  ; 
	  whatToReturn = ficheTAXClosestAncestorStatement (gmp) ; 
	}
      else 
	{
	  if (ac_has_tag (oKantor, "Blastp_date") &&
	      (ac_has_tag (oKantor, "Taxblast_date") || !ac_has_tag (oKantor, "Blastp")))
	    {
	      spc=ficheDBCreatureID (db) ; 
	      whatToReturn = strnew (SpcI[spc].speciesName, 0) ;
	    }
	  else whatToReturn=0 ; 
	}
    }
  ac_free (oKantor) ;  
  ac_db_close (db) ; 
  messfree (gmp) ; 
  return  whatToReturn ;
}


/* -===================================- /
 * -=  BEST NAME                       =- /
 * -===================================- */	

char *fogicBestNameTranscribed_Gene (char *nameTG)
{
  AC_OBJ tg = 0, gene = 0 ;
  vTXT blk = 0 ; 
  AC_DB db ; 
  GMP *gmp ;
  char *cp ;

  if ((db = ac_open_db ("local", 0)) && 
      (tg = ac_get_obj (db, "transcribed_gene", nameTG, 0)) &&
      (gene = ac_tag_obj (tg, "Gene", 0)))
    {
      blk = vtxtCreate () ;
      gmp = gmpCreate (db, gene, 0, 0, 0, 0, 0, 'g') ;
      if (gmp->markup) vtxtMarkup (blk) ;
      ficheTGBestNameStatement (blk, gmp, gmp->tg, 0) ; 
      gmpDestroy (gmp) ;
    } 
  ac_free (tg) ;
  ac_free (gene) ;
  ac_db_close (db) ; 
  cp = strnew (vtxtPtr (blk), 0) ;
  vtxtDestroy (blk) ;

  return cp ;
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
