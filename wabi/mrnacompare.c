#include "regular.h"
#include "acedb.h"
#include "cdna.h"
#include "../wac/ac.h"
#include "vtxt.h"
#include "parse.h"

typedef struct mi_struct { AC_OBJ mrna ; int a1, a2, x1, x2 ; const char * type ; 
  int isExon, isIntron, isStandard, isFuzzy, isFirst, isNH2Complete, isExactFirst, isPossibleFirst, isCentral, isLast, isRedundant;} MI ;

/*********************************************************************/
/*********************************************************************/

static void mcShowMis (Array mis) 
{
  int ii ;
  MI *mi ;

  if (arrayExists(mis) && arrayMax(mis))
    for (ii = 0 ; ii < arrayMax(mis) ; ii++)
      {
	mi = arrayp (mis, ii, MI) ;
	printf ("%d::\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s %s %s %s %s %s%s %s\n", 
		ii, ac_name(mi->mrna), 
		mi->a1, mi->a2, mi->x1, mi->x2, 
		mi->type,
		mi->isFirst && mi->isNH2Complete ? "complete" : "",
		mi->isFirst ? "first" : "",
		mi->isCentral ? "central" : "",
		mi->isLast ? "last" : "",
		mi->isExon ? "exon" : "",
		mi->isFuzzy ? "fuzzy " : "",
		mi->isStandard ? " standard" : "",
		mi->isIntron ? "intron" : "",
		mi->isRedundant ? "redundant" : ""
		) ;
      }
} /* mcShowMis */

/*********************************************************************/

static int mcValidatePossibleFirstExons (Array mis)
{
  MI *mi, *mj ;
  int ii, jj, nn = 0 ;
 
  if (arrayMax (mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isPossibleFirst)
	  continue ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii != jj &&
		mj->isExon &&
		mi->x1 >= mj->x1 && mi->x1 <= mj->x2)
	      { 
		mi->isPossibleFirst = FALSE ;
		break ;
	      }
	  }
	if (mi->isPossibleFirst)
	  mi->isFirst = TRUE ;
	else
	  mi->isCentral = FALSE ;    /* eaten up */
	mi->isPossibleFirst = TRUE ; /* needed in ovelappingExons */
      }
  return nn;
} /* mcValidatePossibleFirstExons */

/*********************************************************************/

static Array mcGetMis (AC_DB db, KEY gene, AC_TABLE mrnas, AC_HANDLE h, BOOL isCoding, int *nMrnap, int *nMrnaUsplp)
{
  Array mis = arrayHandleCreate (200, MI, h) ;
  int i, ii, jj, a1, a2, x1, x2, y1, y2, ymax = 0, c1, c2 ;
  AC_OBJ mrna = 0, product = 0 ;
  MI *mi;
  AC_TABLE spl = 0 ;
  KEY key ;

  *nMrnaUsplp = *nMrnap = 0 ;
  for (ii = jj = 0 ; ii < mrnas->rows ; ii++)
    {
      ac_free (spl) ;
      ac_free (mrna) ;
      ac_free (product) ;
      mrna = ac_table_obj (mrnas, ii, 0, h) ;
      if (! ac_has_tag (mrna, "gt_ag") && ! ac_has_tag (mrna, "gc_ag"))
	{ (*nMrnaUsplp)++ ; continue ; }
      key = ac_tag_key (mrna, "Gene", 0) ;
      if (key != gene)
	continue ;
      (*nMrnap)++ ;
      c1 = ac_tag_int (mrna, "Length_5prime_UTR", 0) ;
      c2 = ac_tag_int (mrna, "Length_3prime_UTR", 0) ;
      product = ac_tag_obj (mrna, "Product", h) ;
      a1 = ac_table_int (mrnas, ii, 1, -9999) ;
      a2 = ac_table_int (mrnas, ii, 2, -9999) ;
      spl = ac_tag_table (mrna, "Splicing", h) ;
      /*
      spl =  ac_aql_table (db, messprintf("select  a1, a2, b1, b2, t, fuz from m in class \"mrna\" where m like \"%s\", a1 in m->Splicing, a2 in a1[1],  b1 in a1[2], b2 in a1[3], t in a1[4], fuz in a1[5]", ac_name(mrna)), NULL, h) ;
      */
      if (spl->rows) 
	ymax = ac_table_int (spl, spl->rows - 1, 3, -9999) ;
      for (i = 0 ; i < spl->rows ; i++)
	{
	  if (0 &&
	      strstr(ac_table_tag (spl, i, 4, ""), "tolen") &&
	      ! (i == 0 && ac_tag_table (mrna, "Transpliced_to", 0)))
	    continue ;
	  if (0 &&
	      strstr(ac_table_tag (spl, i, 4, ""), "redicted"))
	    continue ;

	  x1 = ac_table_int (spl, i, 0, -9999) ; /* genome coord */
	  x2 = ac_table_int (spl, i, 1, -9999) ; /* genome coord */
	  y1 = ac_table_int (spl, i, 2, -9999) ; /* spliced coord */
	  y2 = ac_table_int (spl, i, 3, -9999) ; /* spliced coord */
	  if (isCoding && y2 < c1) continue ;
	  if (isCoding && y1 < c1) { x1 += c1 - y1 ; y1 = c1 ; }
	  if (isCoding && y1 > ymax - c2) continue ;
	  if (isCoding && y2 > ymax - c2) { x2 -= (y2 - ymax + c2) ; y2 = ymax - c2 ; }
	  mi = arrayp (mis, jj++, MI) ;
	  mi->mrna = mrna ; mi->a1 = a1 ; mi->a2 = a2 ;
	  mi->x1 = a1 + x1 ;
	  mi->x2 = a1 + x2 ;

	  mi->type = ac_table_tag (spl, i, 4, "") ;
	  if (strstr (mi->type, "xon"))  mi->isExon = TRUE ;
	  if (strstr (mi->type, "ron"))
	    mi->isIntron = TRUE ;
	  if (strstr (mi->type, "ron") && 
	      strstr(ac_table_text (spl, i, 5, "toto"), "uzz"))
	    mi->isFuzzy = TRUE ;
	  if (strstr (mi->type, "ron") && 
	      (
	       strstr(ac_table_text (spl, i, 5, "toto"), "gt_ag") ||
	       strstr(ac_table_text (spl, i, 5, "toto"), "gc_ag")
	       ))
	    mi->isStandard = TRUE ;
	  if (i == 0 && mi->isExon)
	    {	
	      AC_TABLE tt1 = 0 ;
	      int x ;
	      const char *sl ;

	      mi->isPossibleFirst = TRUE ;
	      if (ac_has_tag (product, "NH2_Complete"))
		{ mi->isFirst = TRUE ; mi->isNH2Complete = TRUE ; }
	      if ((tt1 = ac_tag_table (mrna, "Transpliced_to", 0)))
		{
		  for (x = 0 ; x < tt1->rows ; x++)
		    {
		      sl = ac_table_printable (tt1, x, 0, 0) ;
		      if (sl &&
			  !strncasecmp (sl, "SL", 2) &&
			  *(sl+2) > '0')
			mi->isExactFirst = TRUE ;
		    }
		  ac_free (tt1) ;
		}
	    }
	  if (i == spl->rows - 1 &&  mi->isExon && ac_has_tag (product, "COOH_Complete"))
	    mi->isLast = TRUE ;
	  if (/* i > 0 && i < spl->rows - 1 && */  mi->isExon)
	    mi->isCentral = TRUE ;
	} 
      /*
	a1 = 999999;
	for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
	{
	if (mi->x1 < a1) { a1 = mi->x1 ; i = ii ; }
	}
	arrp (mis, 0, MI)->isFirst = TRUE ;
      */

    }
  ac_free (spl) ;
  ac_free (mrna) ; 
  ac_free (product) ;
  mcValidatePossibleFirstExons (mis) ;
  return mis ;
} /* mcGetMis */

/*********************************************************************/

static int mcCountSharedIntrons (Array mis)
{
  MI *mi, *mj ;
  int ii, jj, n, nn = 0 ;

  if (arrayMax(mis))
    {
      for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
	mi->isRedundant = FALSE ;
      for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
	{
	  if (!mi->isIntron)
	    continue ;
	  if (ii < arrayMax(mis) - 1)
	    for (n = 0, jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	      {
		if (mj->isIntron &&
		    mi->x1 == mj->x1 && mi->x2 == mj->x2)
		  { 
		    if (n) nn++ ;
		    mj->isRedundant = TRUE ; continue ;
		  }
	      }
	}
    }
  return nn;
}

/*********************************************************************/

static int mcCountSharedFirstExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, n, nn = 0 ;

  /* remove the first exon label from redundant cases */
  /* first treat all cases where an exact can eat another exact or not first */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isExactFirst || !mi->isFirst)
	  continue ;
	for (n= 0, jj = 0, mj = arrp (mis, jj, MI) ; mi->isFirst && jj < arrayMax(mis) ; mj++, jj++)
	  { 
	    if (ii == jj || !mj->isFirst)
	      continue ;
	    if (mi->x2 == mj->x2)
	      {
		if (mi->x1 == mj->x1)
		  { 
		    mj->isFirst = FALSE ; 
		    if (!n++) nn++ ;
		    continue ;
		    
		  }
		if (!mj->isExactFirst && 
		    mi->x1 <= mj->x1)
		  { 
		    mj->isFirst = FALSE ; 
		    if (!n++) nn++ ;
		    continue ;		
		  }
	      }
	  }
      }
	    
  /* now treat the general case
    we need 2 rounds because a long and short first may merge
    together but not both with an exact, hence we do not
    have transitivity
  */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isFirst)
	  continue ;
	if (ii < arrayMax(mis) - 1)
	  for (n= 0, jj = ii+1, mj = arrp (mis, jj, MI) ; mi->isFirst && jj < arrayMax(mis) ; mj++, jj++)
	    { 
	      if (!mj->isFirst)
		continue ;
	      if (mi->x2 == mj->x2)
		{
		  if (mi->isExactFirst && mj->isExactFirst && 
		      mi->x1 == mj->x1)
		    { 
		      mj->isFirst = FALSE ; 
		      if (!n++) nn++ ;
		      continue ;
		    }
		  if (mi->isExactFirst && !mj->isExactFirst && 
		      mi->x1 <= mj->x1)
		    { 
		      mj->isFirst = FALSE ; 
		      if (!n++) nn++ ;
		      messcrash ("should not happen") ;
		    }
		  if (!mi->isExactFirst && mj->isExactFirst && 
		      mi->x1 >= mj->x1)
		    { 
		      mi->isFirst = FALSE ; 
		      if (!n++) nn++ ;
		      messcrash ("should not happen") ;
		    }
		  if (!mi->isExactFirst && !mj->isExactFirst)
		    { 
		      mj->isFirst = FALSE ; 
		      if (!n++) nn++ ;
		      continue ;		
		    }
		}
	    }
      }
  return nn ;
} /* mcCountSharedFirstExons */

/*********************************************************************/

static int mcCountAltFirstExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;
  BOOL foundFirst = FALSE ;

  /* see if my first exon is alternative */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isFirst) continue ;
	foundFirst = TRUE ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii == jj) continue ;
	    if (
		(mi->isFirst && mi->x2 > mj->x2) || /* somebody to my left */
		
		(mj->isFirst && mi->x2 <= mj->x2)   /* a first to my right, i count on shared beeing unflagged  */
		)
	      {  nn++ ; break ; } /* loop on next ii */
	  }
      }
  if (!nn && foundFirst) nn = 1 ;
  return nn ;
} /* mcCountAltFirstExon */

/*********************************************************************/

static int mcCountNonOverlappingAltFirstExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;
  BOOL foundFirst = FALSE ;

  /* see if my first exon is alternative */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      { 
	if (!mi->isFirst) continue ;
	foundFirst = TRUE ;
	nn = 0 ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii == jj) continue ;
	    if (
		(mi->isFirst && mi->x2 > mj->x2) || /* somebody to my left */
		
		(mj->isFirst && mi->x2 <= mj->x2)   /* a first to my right, i count on shared beeing unflagged  */
		)
	      {  nn++ ; break ; } /* loop on next ii */
	  }
	if (nn) mi->isFirst = 3 ;
      } 

  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      { 
	if (mi->isFirst < 2) continue ;
	if (ii < arrayMax(mis) - 1)
	  for (jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	    {
	      if (ii == jj) continue ;
	      if (mj->isFirst < 3) continue ;
	      if (mi->x1 < mj->x2 && mi->x2 > mj->x1) /* contact */
		mj->isFirst = 2 ; 
	    }
      } 
  
  nn = 0 ;
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      { 
	switch (mi->isFirst)
	  {
	case 3: nn++ ; mi->isFirst = 1 ; break ;
	  case 2: mi->isFirst = 1 ; break ;
	  }
      } 
  if (!nn && foundFirst) nn = 1 ;
  return nn ;
} /* mcCountNonOverlappingAltFirstExons */

/*********************************************************************/

static int mcCountPossiblePromotors (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;
  BOOL foundFirst = FALSE ;

  /* see if my first exon is alternative */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      { 
	if (!mi->isFirst) continue ;
	foundFirst = TRUE ;
	nn = 0 ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii == jj) continue ;
	    if (
		(mi->isFirst && mi->isNH2Complete && mi->x2 > mj->x2) || /* somebody to my left */
		
		(mj->isFirst && mj->isNH2Complete && mi->x2 <= mj->x2)   /* a first to my right, i count on shared beeing unflagged  */
		)
	      {  nn++ ; break ; } /* loop on next ii */
	  }
	if (nn) mi->isFirst = 3 ;
      } 

  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      { 
	if (mi->isFirst < 2) continue ;
	if (ii < arrayMax(mis) - 1)
	  for (jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	    {
	      if (ii == jj) continue ;
	      if (mj->isFirst < 3) continue ;
	      if (mi->x1 < mj->x2 && mi->x2 > mj->x1) /* contact */
		mj->isFirst = 2 ; 
	    }
      } 
  
  nn = 0 ; 
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      { 
	switch (mi->isFirst)
	  {
	  case 3: nn++ ; mi->isFirst = 1 ; break ;
	  case 2: mi->isFirst = 1 ; break ;
	  }
      } 
  if (!nn && foundFirst) nn = 1 ;
  return nn ;
} /* mcCountPossiblePromotors */

/*********************************************************************/

static int mcCountOverlappingCentralExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;
  
  /* remove the last exon label from redundant cases */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isCentral || mi->isFirst || mi->isLast || mi->isPossibleFirst)
	  continue ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii == jj ||
		mi->x2 < mj->x1 || mi->x1 > mj->x2 ||
		!mj->isCentral || mj->isFirst || mj->isLast || mj->isPossibleFirst)
	      continue ;
	    {  nn++ ; break ; } /* loop on next ii */
	  }
      }
  return nn ;
} /* mcCountAltFirstExon */

/*********************************************************************/

static int mcCountSkippedCentralExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;
  
  /* remove the last exon label from redundant cases */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isCentral || mi->isFirst || mi->isLast)
	  continue ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii == jj ||
		mi->x1 < mj->x1  || mi->x2 > mj->x2 ||
		!mj->isIntron)
	      continue ;
	    {  nn++ ; break ; } /* loop on next ii */
	  }
      }
  return nn ;
} /* mcCountSkippedCentralExons */

/*********************************************************************/

static int mcCount5pSkippedCentralExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;
  
  /* remove the last exon label from redundant cases */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isCentral || mi->isFirst || mi->isLast)
	  continue ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii == jj ||
		mi->x2 > mj->x1 ||
		!mj->isFirst)
	      continue ;
	    {  nn++ ; break ; } /* loop on next ii */
	  }
      }
  return nn ;
} /* mcCount5pSkippedCentralExons */

/*********************************************************************/

static int mcCount3pSkippedCentralExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;
  
  /* remove the last exon label from redundant cases */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isCentral || mi->isFirst || mi->isLast)
	  continue ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii == jj ||
		mi->x1 < mj->x2 ||
		! mj->isLast)
	      continue ;
	    {  nn++ ; break ; } /* loop on next ii */
	  }
      }
  return nn ;
} /* mcCount3pSkippedCentralExons */

/*********************************************************************/

static int mcCountSharedCentralExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, n, nn = 0 ;
  int realCentral ;
  
  /* remove the last exon label from redundant cases */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isCentral)
	  continue ;
	realCentral = 0 ;
	if (!mi->isFirst && !mi->isLast) realCentral = 1 ;
	if (ii < arrayMax(mis) - 1)
	  for (n = 0, jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	    {
	      if (!mj->isCentral)
		continue ;
	      if (mi->x1 == mj->x1 &&
		  mi->x2 == mj->x2)
		{ 
		  mj->isCentral = FALSE ; 
		  if (!mj->isFirst && !mj->isLast) realCentral++ ;
		  if (realCentral >= 2 && !n++) nn++ ;
		}
	    }
      }
  
  return nn ;
} /* mcCountSharedCentralExons */

/*********************************************************************/

static int mcCountSharedLastExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, n, nn = 0 ;

  /* remove the last exon label from redundant cases */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isLast)
	  continue ;
	if (ii < arrayMax(mis) - 1)
	  for (n = 0, jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	    {
	      if (!mj->isLast)
		continue ;
	      if (mi->x1 == mj->x1)
		{ 
		  if (!n++) nn++ ;
		  mj->isLast = FALSE ; 
		  mj->isCentral = FALSE ;
		}
	    }
      }
  return nn ;
} /* mcCountSharedLastExons */

/*********************************************************************/

static int mcCountAltLastExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;
  BOOL foundLast = FALSE ;

  /* see if my last exon is alternative */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isLast)
	  continue ;
	foundLast = TRUE ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; mi->isLast && jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii == jj) continue ;
	    if (
		(mi->x1 < mj->x1) || /* somebody to my right */
		
		(mj->isLast)  /* a last to my left */
		)
	      {  nn++ ; break ; } /* loop on next ii */
	  }
      }
  if (!nn && foundLast) nn = 1 ;
  return nn ;
} /* mcCountAltLastExon */

/*********************************************************************/

static int mcCountNonOverlappingAltLastExons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;
  BOOL foundLast = FALSE ;

  /* see if my first exon is alternative */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      { 
	if (!mi->isLast)
	  continue ;
	foundLast = TRUE ;
	nn = 0 ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii == jj) continue ;
	    if (
		(mi->x1 < mj->x1) || /* somebody to my right */
		
		(mj->isLast)  /* a last to my left */
		)
	      { nn++ ; break ; } /* loop on next ii */
	  }
	if (nn) mi->isLast = 3 ;
      } 
  
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      { 
	if (mi->isLast < 2) continue ;
	if (ii < arrayMax(mis) - 1)
	  for (jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	    {
	      if (ii == jj) continue ;
	      if (mj->isLast < 3) continue ;
	      if (mi->x1 < mj->x2 && mi->x2 > mj->x1) /* contact */
		mj->isLast = 2 ; 
	    }
      } 
  
  nn = 0 ;
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      { 
	switch (mi->isLast)
	  {
	  case 3: nn++ ; mi->isLast = 1 ; break ;
	  case 2: mi->isLast = 1 ; break ;
	  }
      } 
  if (!nn && foundLast) nn = 1 ;
  return nn ;
} /* mcCountNonOverlappingAltLastExons */

/*********************************************************************/

static int mcCountExons (Array mis)
{
  MI *mi ;
  int ii, nn = 0 ;
  
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isExon)
	  continue;
	if (mi->isFirst || mi->isLast || mi->isCentral)
	  nn++ ;
      }
  return nn ;
} /* mcCountExons */

/*********************************************************************/

static int mcCountIntrons (Array mis)
{
  MI *mi ;
  int ii, nn = 0 ;

  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isIntron || mi->isRedundant)
	  continue;
	nn++ ;
      }
  return nn ;
} /* mcCountIntrons */

/*********************************************************************/

static int mcCountStandardIntrons (Array mis)
{
  MI *mi ;
  int ii, nn = 0 ;

  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isIntron || mi->isRedundant || !mi->isStandard)
	  continue;
	nn++ ;
      }
  return nn ;
} /* mcCountStandardIntrons */

/*********************************************************************/

static int  mcCountRetainedIntrons (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;

  /* see if my last exon is alternative */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isIntron || mi->isRedundant)
	  continue ;
	/* do total square count because one is intron other is exon */
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii != jj &&
		mj->isExon && !mj->isRedundant &&
		mi->x1 > mj->x1 && mi->x2 < mj->x2)
	      { nn++ ; break ; }
	  }
      }
  return nn ;
} /* mcCountRetainedIntrons */

/*********************************************************************/

static int  mcCountMutExclusives (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, n, nn = 0, nk = 0, xx2 ;

  /* see if 2 sucessive mi intron end like 2 succesive mj introns */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) - 3 ; mi++, ii++)
      {
	if (!mi->isIntron || !(mi+1)->isExon ||
	    mi->isFuzzy || (mi+2)->isFuzzy ||
	    !(mi+2)->isIntron || (mi+2)->mrna != mi->mrna)
	  continue ;
	xx2 = (mi+2)->x2 ;
	if (ii < arrayMax(mis) - 1)
	  for (n = 0, jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) - 3 ; mj++, jj++)
	    { 
	      if (mj->isIntron && !mj->isFuzzy && !(mj+2)->isFuzzy &&
		  mi->x1 == mj->x1 &&
		  (mj+1)->isExon &&
		  (mj+2)->mrna == mj->mrna &&
		  (mj+2)->isIntron &&
		  (mj+2)->x2 == xx2 &&
		  (mi->x2 > (mj+2)->x1 || mj->x2 > (mi+2)->x1) /* middle exons do not overlap */
		  )
		{ mi->isIntron = 2 ; if (!n++) nn++ ; break ; }
	    }
      } 
  if (nn)
    {
      for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
	{
	  if (mi->isIntron != 2)
	    continue ; 
	  mi->isIntron = 1 ; 
	  nk++ ;
	  if (ii < arrayMax(mis) - 1)
	    for (jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) - 2 ; mj++, jj++)
	      { 
		if (mj->isIntron != 2)
		  continue ; 
		if (mi->x1 == mj->x1 && mi->x2 == mj->x2)
		  mj->isIntron = 1 ;
	      }
	}
    }
  
  return nk ;
} /* mcCountMutExclusives */

/*********************************************************************/

static int  mcCountInternalDonors (Array mis, int limit)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0, dx ;
  
  /* look for diff in donor site */ 
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isIntron || mi->isFuzzy || mi->isRedundant)
	  continue ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  { /* look on all, do not filter redundant mj, since we break anayway and
	       since they may have a longer preceding exon */
	    if (ii != jj &&
		mj->isIntron &&  !mj->isFuzzy &&
		mi->x1 < mj->x1 &&
		mi->x2 == mj->x2 &&
		(mj-1)->mrna == mj->mrna &&
		(mj-1)->isExon &&
		(mj-1)->x1 < mi->x1 &&
		(mj-1)->x2 > mi->x1
		)
	      {
		dx =  mj->x1 - mi->x1 ;
		if (
		    (limit == 0) ||
		    (dx == limit)
		    )
		{ nn++ ; break ; }
	      }
	  }
      }
  return nn ;
} /* mcCountInternalDonorSite */

/*********************************************************************/

static int  mcCountInternalAcceptors (Array mis, int limit)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0, dx ;
  
  /* look for diff in acceptor site */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isIntron || mi->isFuzzy || mi->isRedundant)
	  continue ;
	for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	  {
	    if (ii != jj &&
		mj->isIntron && !mj->isFuzzy &&
		mi->x2 < mj->x2 &&
		mi->x1 == mj->x1 &&
		(mi+1)->mrna == mi->mrna &&
		(mi+1)->isExon &&
		(mi+1)->x1 < mj->x2 &&
		(mi+1)->x2 > mj->x2 
		)
	      {
		dx = mi->x2 > mj->x2 ? mi->x2 - mj->x2 :  mj->x2 - mi->x2 ;
		if (
		    (limit == 0) ||
		    (dx == limit)
		    )
		  { nn++ ; break ; }
	      }
	  }
      }
  return nn ;
} /* mcCountInternalAcceptorSite */

/*********************************************************************/
 
static int  mcCountCassettes (Array mis, int limit, int *nkp, int *nkp2)
{
  MI *mi, *mj, *mk ;
  int ii, jj, kk, n, nn = 0, nk = 0, nexon ;
  AC_OBJ mrna ;

  /* see if 2 (mj mk) successive (or immediataly successive) introns match mi intron */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->isIntron  || mi->isFuzzy || mi->isRedundant)
	  continue ; 
	for (n = jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) - 2 ; mj++, jj++)
	  {
	    if (ii == jj) continue ;
	    if (mj->isIntron && mi->x1 == mj->x1)
	      {
		mrna = mj->mrna ; nexon = 0 ;
		for (kk = jj+1, mk = arrp (mis, kk, MI) ; kk < arrayMax(mis) ; mk++, kk++)
		  {
		    if (mrna != mk->mrna)
		      break ;
		    if (mk->isExon) nexon++ ;
		    if (limit && nexon > limit)
		      break ;
		    if (mk->isIntron && mk->x2 == mi->x2)
		      {
			if (!limit && nexon < 2) break ;
			if (!n++) nn++ ; 
			while (mk-- > mj) if (mk->isExon) mk->isExon |= 0x6 ;
			goto nextjj ;
		      }
		  }
	      } 
	  nextjj:		  
	    continue ;
	  }
      }
  if (nn)
    {
      for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
	{
	  if (!(mi->isExon & 0x2))
	    continue ; 
	  mi->isExon &= ~0x2 ;
	  nk++ ;
	  if (ii < arrayMax(mis) - 1)
	    for (jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) - 2 ; mj++, jj++)
	      { 
		if (!(mj->isExon & 0x2))
		  continue ; 
		if (mi->x1 == mj->x1 && mi->x2 == mj->x2)
		  mj->isExon &= ~0x2 ;
	      }
	}
    }
  *nkp = nk ;
  if (!limit) /* second pass */
    {
      nk = 0 ; 
      if (arrayMax(mis))
	for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
	  {
	    if (!(mi->isExon & 0x4))
	      continue ; 
	    mi->isExon &= ~0x4 ;
	    nk++ ;
	    if (ii < arrayMax(mis) - 1)
	      for (jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) - 2 ; mj++, jj++)
		{ 
		  if (!(mj->isExon & 0x4))
		    continue ; 
		  if (mi->x1 == mj->x1 && mi->x2 == mj->x2)
		    mj->isExon &= ~0x4 ;
		}
	  }
      *nkp2 = nk ;
    }
  
  return nn ;
} /* mcCountCassettes */

/*********************************************************************/

static int  mcCountCitroenIntrons (Array mis)
{
  MI *mi, *mj, *mk  ;
  int ii, jj, kk, nn = 0 ;
  AC_OBJ mrna ;
  BOOL foundIntron ;

  /* see 2 feet of mi intron go in 2 exons of shared mj mrna */
 if (arrayMax(mis))
   for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
     {
       if (!mi->isIntron  || mi->isFuzzy || mi->isRedundant)
	 continue ;
       for (jj = 0, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) - 2 ; mj++, jj++)
	 {
	   if (ii != jj &&
	       mj->isExon &&
	       (mj->x1 < mi->x1 && mj->x2 > mi->x1))
	     {
	       mrna = mj->mrna ; foundIntron = FALSE ;
	       for (kk = jj+1, mk = arrp (mis, kk, MI) ; kk < arrayMax(mis) ; mk++, kk++)
		 {
		   if (mrna != mk->mrna)
		     break ;
		   if (mk->isIntron && !mk->isFuzzy) /* avoid 2 exons separated by a gap */
		     foundIntron = TRUE ;
		   if (
		       foundIntron &&
		       mk->isExon &&
		       mk->mrna == mrna &&
		       mk->x1 < mi->x2 && mk->x2 > mi->x2
		       ) 
		     { nn++ ; goto nextii ; }
		 }
	     }
	 }
     nextii:
       continue ;
     }
 return nn ;
} /* mcCountSuspectedDeletion */

/*********************************************************************/
#ifdef JUNK
 prototype for similar code
static int  mcCount (Array mis)
{
  MI *mi, *mj  ;
  int ii, jj, nn = 0 ;

  /* see */
  if (arrayMax(mis))
    for (ii = 0, mi = arrp (mis, 0, MI)  ; ii < arrayMax(mis) ; mi++, ii++)
      {
	if (!mi->mrna)
	  continue ;
	if (ii < arrayMax(mis) - 1)
	  for (jj = ii+1, mj = arrp (mis, jj, MI) ; jj < arrayMax(mis) ; mj++, jj++)
	    {
	    }
      }
  return nn ;
} /* mcCount */
#endif

/*********************************************************************/

typedef struct res_struct { /* PURE unsigned int, NEEDED becuase of mcCumul */
  unsigned int ntg, nmrna, nMrnaUspl, nExons, nIntrons, nStandardIntrons ;
  unsigned int nSharedIntrons ;
  unsigned int nSharedFirstExons ;
  unsigned int nAltFirstExons ;
  unsigned int nNonOverlappingAltFirstExons ;
  unsigned int nPossiblePromotors ;
  unsigned int nNonOverlappingAltLastExons ;
  unsigned int nSharedCentralExons ;
  unsigned int nOverlappingCentralExons ;
  unsigned int nSkippedCentralExons ;
  unsigned int n5pSkippedCentralExons ;
  unsigned int n3pSkippedCentralExons ;
  unsigned int nSharedLastExons ;
  unsigned int nAltLastExons ;
  unsigned int nRetainedIntrons ;
  unsigned int nMutExclusives ;
  unsigned int nInternalDonorSite ;
  unsigned int nInternalDonorSite3 ;
  unsigned int nInternalDonorSite6 ;
  unsigned int nInternalDonorSite9 ;
  unsigned int nInternalDonorSite12 ;
  unsigned int nInternalDonorSite15 ;
  unsigned int nInternalDonorSite18 ;
  unsigned int nInternalAcceptorSite ;
  unsigned int nInternalAcceptorSite3 ;
  unsigned int nInternalAcceptorSite6 ;
  unsigned int nInternalAcceptorSite9 ;
  unsigned int nInternalAcceptorSite12 ;
  unsigned int nInternalAcceptorSite15 ;
  unsigned int nInternalAcceptorSite18 ;
  unsigned int nCassettedExons ;
  unsigned int nCassettingIntrons ;
  unsigned int nMultiCassettingIntrons ;
  unsigned int nMultiCassettedExons ;
  unsigned int nAnyCassettedExons ;
  unsigned int nCitroenIntrons ;
} RES ;

static void mcAnalyseMis (Array mis, RES *r)
{
  int nn1=0, nn2=0, nn3=0 ;
  r->nSharedIntrons += mcCountSharedIntrons (mis) ;
  r->nSharedCentralExons += mcCountSharedCentralExons (mis) ; /* must come before first/last */

  r->nOverlappingCentralExons += mcCountOverlappingCentralExons  (mis) ;
  r->nSkippedCentralExons  += mcCountSkippedCentralExons  (mis) ;
  r->n5pSkippedCentralExons  += mcCount5pSkippedCentralExons  (mis) ;
  r->n3pSkippedCentralExons  += mcCount3pSkippedCentralExons  (mis) ;

  r->nSharedFirstExons += mcCountSharedFirstExons (mis) ;
  r->nAltFirstExons += mcCountAltFirstExons (mis) ;
  r->nNonOverlappingAltFirstExons += mcCountNonOverlappingAltFirstExons (mis) ;
  r->nPossiblePromotors += mcCountPossiblePromotors (mis) ;
  r->nSharedLastExons += mcCountSharedLastExons (mis) ;
  r->nAltLastExons += mcCountAltLastExons (mis) ;
  r->nNonOverlappingAltLastExons += mcCountNonOverlappingAltLastExons (mis) ;
  r->nExons += mcCountExons (mis) ; /* must come after shared first last */
  r->nIntrons += mcCountIntrons (mis) ;
  r->nStandardIntrons += mcCountStandardIntrons (mis) ;

  r->nRetainedIntrons += mcCountRetainedIntrons (mis) ;
  r->nMutExclusives += mcCountMutExclusives (mis) ;
  r->nInternalDonorSite += mcCountInternalDonors (mis, 0) ;
  r->nInternalDonorSite3 += mcCountInternalDonors (mis, 3) ;
  r->nInternalDonorSite6 += mcCountInternalDonors (mis, 6) ;
  r->nInternalDonorSite9 += mcCountInternalDonors (mis, 9) ;
  r->nInternalDonorSite12 += mcCountInternalDonors (mis, 12) ;
  r->nInternalDonorSite15 += mcCountInternalDonors (mis, 15) ;
  r->nInternalDonorSite18 += mcCountInternalDonors (mis, 18) ;
  r->nInternalAcceptorSite += mcCountInternalAcceptors (mis, 0) ;
  r->nInternalAcceptorSite3 += mcCountInternalAcceptors (mis, 3) ;
  r->nInternalAcceptorSite6 += mcCountInternalAcceptors (mis, 6) ;
  r->nInternalAcceptorSite9 += mcCountInternalAcceptors (mis, 9) ;
  r->nInternalAcceptorSite12 += mcCountInternalAcceptors (mis, 12) ;
  r->nInternalAcceptorSite15 += mcCountInternalAcceptors (mis, 15) ;
  r->nInternalAcceptorSite18 += mcCountInternalAcceptors (mis, 18) ;
  r->nCassettingIntrons += mcCountCassettes (mis, 1, &nn1, &nn2) ;
  r->nMultiCassettingIntrons += mcCountCassettes (mis, 0, &nn3, &nn2) ;
  r->nCassettedExons += nn1 ;
  r->nAnyCassettedExons += nn2 ;
  r->nMultiCassettedExons += nn3 ;

  r->nCitroenIntrons += mcCountCitroenIntrons (mis) ;
} /* analyseMis */

static void mcShowRes (RES *r)
{
  printf ("\n\nComparisons\n\n") ;
  printf ("ntg=%d,  nSpliced_mRNA=%d, nUnspliced_mRNA=%d,\n  nExons=%d  nIntrons=%d  nStandardIntrons=%d\n\n",	
	  r->ntg, 
	  r->nmrna, r->nMrnaUspl,
	  r->nExons, r->nIntrons, r->nStandardIntrons
	  ) ;

  printf("AltFirstExons=%d\nnNonOverlappingAltFirstExons=%d\nnPossiblePromotors=%d\n\n"
	 ,r->nAltFirstExons
	 ,r->nNonOverlappingAltFirstExons 
	 ,r->nPossiblePromotors
	 ) ;

  printf("n5pSkippedCentralExons=%d\nnSkippedCentralExons=%d\n\tnAnyCassettedExons=%d\n"
	 ,r->n5pSkippedCentralExons
	 ,r->nSkippedCentralExons
	 ,    r->nAnyCassettedExons
	 ) ;

  printf("\t\tnCassettedExons=%d\n\t\tnCassettingIntrons=%d\n\t\tnMultiCassettedExons=%d\n\t\tnMultiCassettingIntrons=%d\n\t\tnMutExclusives=%d\n\t\tnCitroenIntrons=%d\n"
	 ,        r->nCassettedExons
	 ,        r->nCassettingIntrons
	 ,        r->nMultiCassettedExons
	 ,        r->nMultiCassettingIntrons
	 ,        r->nMutExclusives 
	 ,        r->nCitroenIntrons
	 ) ;

  printf("nOverlappingCentralExons=%d\n\tnInternalDonorSite=%d\n\t\tnInternalDonorSite369_18=%d %d %d %d %d %d\n"
	 ,r->nOverlappingCentralExons
	 ,    r->nInternalDonorSite
	 ,        r->nInternalDonorSite3
	 ,        r->nInternalDonorSite6
         ,        r->nInternalDonorSite9
         ,        r->nInternalDonorSite12
 	 ,        r->nInternalDonorSite15
	 ,        r->nInternalDonorSite18
	 ) ;

  printf("\tnInternalAcceptorSite=%d\n\t\tnInternalAcceptorSite369_18 %d %d %d %d %d %d\n"
	 ,    r->nInternalAcceptorSite
	 ,        r->nInternalAcceptorSite3
         ,	  r->nInternalAcceptorSite6
	 ,        r->nInternalAcceptorSite9
         ,	  r->nInternalAcceptorSite12
	 ,        r->nInternalAcceptorSite15
	 ,	  r->nInternalAcceptorSite18
	 ) ;

  printf("nRetainedIntrons=%d\n"
	 ,    r->nRetainedIntrons
	 /* ,        r->nOrfBreakingRetainedIntrons */
	 ) ;

  

  printf("n3pSkippedCentralExons=%d\n\nnAltLastExons=%d\nnNonOverlappingAltLastExons=%d\n\n"
	 ,r->n3pSkippedCentralExons
	 ,r->nAltLastExons 
	 ,r->nNonOverlappingAltLastExons 
	 ) ;


  /*
    r->nSharedIntrons,
    r->nSharedFirstExons,
    r->nSharedCentralExons,
    r->nSharedLastExons,
  */
}

/*********************************************************************/
/* actually create an ace file to stote the data in the gene */
static BOOL mcRes2ace (KEY gene, RES *r, vTXT blkp)
{
  vtxtPrintf (blkp, "Gene %s\n-D Structure\n", name (gene)) ;
  /*nAltCentralExons */
  if (r->nIntrons) vtxtPrintf (blkp, "%s %d\n", "nIntrons" , r->nIntrons ) ;
  if (r->nStandardIntrons) vtxtPrintf (blkp, "%s %d\n", "nStandardIntrons" , r->nStandardIntrons ) ;
  if (r->nExons) vtxtPrintf (blkp, "%s %d\n", "nExons" , r->nExons ) ;
  if (0 && r->nSharedIntrons) vtxtPrintf (blkp, "%s %d\n", "nSharedIntrons" , r->nSharedIntrons ) ;
  if (0 && r->nSharedFirstExons) vtxtPrintf (blkp, "%s %d\n", "nSharedFirstExons" , r->nSharedFirstExons ) ;
  if (0 && r->nSharedCentralExons) vtxtPrintf (blkp, "%s %d\n", "nSharedCentralExons" , r->nSharedCentralExons ) ;
  if (r->nOverlappingCentralExons) vtxtPrintf (blkp, "%s %d\n", "nOverlappingCentralExons" , r->nOverlappingCentralExons ) ;
  if (r->nSkippedCentralExons) vtxtPrintf (blkp, "%s %d\n", "nSkippedCentralExons" , r->nSkippedCentralExons ) ;
  if (r->n5pSkippedCentralExons) vtxtPrintf (blkp, "%s %d\n", "n5pSkippedCentralExons" , r->n5pSkippedCentralExons ) ;
  if (r->n3pSkippedCentralExons) vtxtPrintf (blkp, "%s %d\n", "n3pSkippedCentralExons" , r->n3pSkippedCentralExons ) ;
  if (0 && r->nSharedLastExons) vtxtPrintf (blkp, "%s %d\n", "nSharedLastExons" , r->nSharedLastExons ) ;
  if (r->nAltFirstExons) vtxtPrintf (blkp, "%s %d\n", "nAltFirstExons" , r->nAltFirstExons ) ;
  if (r->nmrna) vtxtPrintf (blkp, "%s %d\n", "nSpliced_mRNA" , r->nmrna ) ;
  if (r->nMrnaUspl) vtxtPrintf (blkp, "%s %d\n", "nUnspliced_mRNA" , r->nMrnaUspl ) ;
  if (r->nNonOverlappingAltFirstExons) vtxtPrintf (blkp, "%s %d\n", "nNonOverlappingAltFirstExons" , r->nNonOverlappingAltFirstExons) ; 
  if (r->nPossiblePromotors) vtxtPrintf (blkp, "%s %d\n", "nPossiblePromotors" , r->nPossiblePromotors) ; 
  if (r->nNonOverlappingAltLastExons) vtxtPrintf (blkp, "%s %d\n", "nNonOverlappingAltLastExons" , r->nNonOverlappingAltLastExons) ;
  if (r->nAltLastExons) vtxtPrintf (blkp, "%s %d\n", "nAltLastExons" , r->nAltLastExons ) ;
  if (r->nRetainedIntrons) vtxtPrintf (blkp, "%s %d\n", "nRetainedIntrons" , r->nRetainedIntrons ) ;
  if (r->nMutExclusives) vtxtPrintf (blkp, "%s %d\n", "nMutExclusives" , r->nMutExclusives ) ;
  if (r->nInternalDonorSite) vtxtPrintf (blkp, "%s %d\n", "nInternalDonorSite" , r->nInternalDonorSite ) ;
  if (r->nInternalDonorSite3 
      + r->nInternalDonorSite6
      + r->nInternalDonorSite9
      + r->nInternalDonorSite12
      + r->nInternalDonorSite15
      + r->nInternalDonorSite18) 
    vtxtPrintf (blkp, "%s %d %d %d %d %d %d\n"
		, "nInternalDonorSite369_18" 
		, r->nInternalDonorSite3
		, r->nInternalDonorSite6
		, r->nInternalDonorSite9
		, r->nInternalDonorSite12
		, r->nInternalDonorSite15
		, r->nInternalDonorSite18
		) ;
  if (r->nInternalAcceptorSite) vtxtPrintf (blkp, "%s %d\n", "nInternalAcceptorSite" , r->nInternalAcceptorSite ) ;
  if (r->nInternalAcceptorSite3 
      + r->nInternalAcceptorSite6
      + r->nInternalAcceptorSite9
      + r->nInternalAcceptorSite12
      + r->nInternalAcceptorSite15
      + r->nInternalAcceptorSite18) 
    vtxtPrintf (blkp, "%s %d %d %d %d %d %d\n"
		, "nInternalAcceptorSite369_18" 
		, r->nInternalAcceptorSite3
		, r->nInternalAcceptorSite6
		, r->nInternalAcceptorSite9
		, r->nInternalAcceptorSite12
		, r->nInternalAcceptorSite15
		, r->nInternalAcceptorSite18
		) ;
  if (r->nCassettedExons) vtxtPrintf (blkp, "%s %d\n", "nCassettedExons" , r->nCassettedExons ) ;
  if (r->nCassettingIntrons) vtxtPrintf (blkp, "%s %d\n", "nCassettingIntrons" , r->nCassettingIntrons ) ;
  if (r->nMultiCassettedExons) vtxtPrintf (blkp, "%s %d\n", "nMultiCassettedExons" , r->nMultiCassettedExons ) ;
  if (r->nMultiCassettingIntrons) vtxtPrintf (blkp, "%s %d\n", "nMultiCassettingIntrons" , r->nMultiCassettingIntrons ) ;
  if (r->nAnyCassettedExons) vtxtPrintf (blkp, "%s %d\n", "nAnyCassettedExons" , r->nAnyCassettedExons ) ;
  if (r->nCitroenIntrons) vtxtPrintf (blkp, "%s %d\n", "nCitroenIntrons" , r->nCitroenIntrons ) ;
  vtxtPrintf (blkp, "\n") ;
  
  return TRUE ;
}

/*********************************************************************/
#include "query.h"
#include "bs.h"

static BOOL mcAddAlterSpliceDetails (AC_DB db, KEY gene, RES *res, vTXT blkp, BOOL isCoding)
{
  AC_HANDLE h;
  AC_TABLE mrnas = 0 ;
  AC_OBJ Tg = 0 ;
  AC_ITER iter ;
  int nn, nu, nMrna = 0 ;
  BOOL ok = FALSE ;
  Array mis = 0 ;
  OBJ Gene = 0 ;
  KEY tg = 0 ;
  
  h = handleCreate () ;
  /* 2007_02_13: give priority to the main tg over the shed tg */
  iter = gene ? ac_dbquery_iter (db, messprintf ("Find gene %s ; {>transcribed_gene gt_ag || gc_ag} SETELSE {>Transcribed_gene}", freeprotect (name(gene))), h) : 0 ;

  while (ac_free (Tg), iter && (Tg = ac_iter_obj (iter)))
    {
      mrnas = ac_tag_table (Tg, "mrna", h) ;
      /*
	mrnas = ac_aql_table (db, messprintf("select m, a1, a2 from tg in class \"transcribed_gene\" where tg like \"%s\", m in tg->mrna, a1 in m[1], a2 in m[2]", name(tg)), NULL, h) ;
      */
      
      if (mrnas && mrnas->rows)
	{  
	  res->ntg++ ;
	  mis = mcGetMis (db, gene, mrnas, h, isCoding, &nn, &nu) ; 
	  nMrna += nn + nu ;
	  res->nmrna += nn ;
	  res->nMrnaUspl += nu ;
	  if (0) mcShowMis (mis) ;
	  mcAnalyseMis (mis, res) ;
	  if (0)
	    {
	      printf("%s ", name(tg)) ;
	      mcShowRes (res) ;
	    }
	}
    }
  if (gene && nMrna >= 1)
    ok = mcRes2ace (gene, res, blkp) ;

  ac_free (h) ;
  return ok ;
} /* mcAddAlterSpliceDetails */

/*********************************************************************/

static int mcTest (KEYSET ks)
{
  AC_DB db ;
  AC_OBJ tg, mrna ;
  AC_ITER iter ;
  AC_TABLE t ;
  int i, x1, x2, a1, a2 ;
  AC_HANDLE h = handleCreate () ;
  
  if ((db = ac_open_db ("local", 0)))
    {
      iter = ac_query_iter (db, TRUE, "Find tg 1b*", NULL, h) ;
      while ((tg = ac_next_obj (iter)))
	{ 
	  printf ("tg %s ", ac_name (tg)) ;
	  t = ac_tag_table (tg, "assembled_from", h) ;
	  for (i = 0 ; t && i < t->rows; i++)
	    {
	      mrna = ac_table_obj (t, i, 2, 0) ;
	      a1 = ac_table_int (t, i, 0, 999999) ;
	      a2 = ac_table_int (t, i, 1, 999999) ;
	      x1 = ac_table_int (t, i, 3, 999999) ;
	      x2 = ac_table_int (t, i, 4, 999999) ;
	      printf ("tg %s %d %d mrna %s %d %d", ac_name (tg), a1, a2, ac_name(mrna), x1, x2) ;
	      
	      printf ("\n") ;
	    } 
	  printf ("\n") ;
	}
      printf ("\n") ;
    }
  messfree (h) ;
  ac_free (db) ;
  return 1 ;
}

/*********************************************************************/
static void mcCumul (RES *res2, RES *res)
{
  unsigned int *jp = (unsigned int *)res2, *ip = (unsigned int *)res ;
  int i = sizeof(RES)/sizeof(int) ;
  while (i--) *jp++ += *ip++ ;

}
/* called with ks from the acembly menu, with key from makemrna.c */
int mcAddKeysetAlterSpliceDetails (KEYSET ks, KEY key, BOOL isCoding) 
{
  int nn = 0, ii ;
  KEY *kp ;
  AC_DB db ;
  vTXT bfr = 0 ;
  RES *res, *res2 ;

  if (0)
    return mcTest (ks) ;

  bfr = vtxtCreate () ;
  res2 = (RES*) halloc (sizeof (RES), 0) ;
  if ((db = ac_open_db ("local", 0)))
    {
      if (ks && keySetMax(ks))
	{
	  kp = arrp (ks, 0, KEY) - 1 ;
	  ii = keySetMax(ks) ;
	  while (kp++, ii--)
	    {
	      res = (RES*) halloc (sizeof (RES), 0) ;
	      
	      if (mcAddAlterSpliceDetails (db, *kp, res, bfr, isCoding))
		nn++ ;
	      mcCumul (res2, res) ;
	      memset (res, 0, sizeof (RES)) ;
	    }
	}
      else if (key)
	mcAddAlterSpliceDetails (db, key, res2, bfr, isCoding) ;
      if (vtxtPtr (bfr))
	parseBuffer (vtxtPtr (bfr), 0) ;
      if (ks)
	mcShowRes (res2) ;
	
      ac_db_close (db) ;
    }

  messfree (res2) ;

  vtxtDestroy (bfr) ;
  
  return nn ;  
} /* mcAddKeysetAlterSpliceDetails */

/*********************************************************************/
/*********************************************************************/
