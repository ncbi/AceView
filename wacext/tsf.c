/*  File: tsf.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2005
 * -------------------------------------------------------------------
 * AceView is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the AceView project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.acedb.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov

 * 2018_04_01
 * tsf   'tag sample format' 
 * strategic idea
 * Define a format semantically equivalent to a 2 dimensional table format
 * and interchangeable with a table format such that each line in the file
 * is providing the value of a single cell, or a small set of contiguous cells.
 * The order of the lines is immaterial, so they can be produced by the source
 * program emitting the data in random order.
 * When 2 values are associated to the same cell, the program applies
 * an action: update (take last value), add (sum the values), 
 * juxtapose (a then b produces a:b), merge (a then a produces a)
 * compute (a then b produces 7*sin(a) + cos(b)), only in case of 
 * parsing 2 files

 * The format has 3 + n fields
 *     1: T tag  (understood as equivalent to the caption of a column table)
 *     2: S sample  (understood as equivalent to the caption of a column lime)
 *     3: F format (descibes the content of a cell or cell group)
 *   ...: v values conforming to the format

 * There is a similarity with the vcf format (variant cll fromat) used in genetics
 * but the objective and the data are manipulated are different.
 
 * For a human readble presentation, use tsf -I tsf -O table, to export as table
 */

#define VERSION "1.1"

/*
#define ARRAY_CHECK  
#define MALLOC_CHECK  
*/
#include "ac.h"
#include "topology.h"

/* snp chaining distance */

typedef struct tsf_struct {
  AC_HANDLE h ;
  
  const char *inFileList ;
  const char *outFileName ;
  const char *separator ;
  const char *title ;
  const char *sample ;
  const char *sortCommand ;
  const char *tagSelectList ;
  const char *tagRenameList ;
  const char *sampleSelectList ;
  const char *compute ;

  const char *tagSeparator ;    /* default : */
  const char *sampleSeparator ; /* default : */
  const char *valueSeparator ;  /* default , */
  const char *fieldSeparator ;  /* default \t */
  
  DICT *tagDict ;
  DICT *sampleDict ;
  DICT *tagOldNameDict ;
  DICT *tagNewNameDict ;
  DICT *tagSelectDict ;
  DICT *sampleSelectDict ;
  DICT *tableCaptionDict ;
  DICT *columnCaptionDict ;
  DICT *valDict ;
  KEYSET tagRenameKs ;
  KEYSET skipKs ;
  KEYSET tagDepth ;
  AC_TABLE results ;

  BOOL gzi, gzo, transpose, debug ;
  BOOL isTableIn, isTableOut ;
  BOOL merge ;
  BOOL sum, min, max, replace ; 
  BOOL tagSelectOnly ;
  BOOL sampleSelectOnly ;
  BOOL chronoOrder ;
  int skip1, skip2, skip3 ;
  int sampleDepth ;

  ACEIN ai ; 
  ACEOUT ao ;
  BigArray hits ;
  int  nInputFiles ;
} TSF ;

#define MAXTAGDEPTH 24
#define MAXCOL 64 

typedef struct hit_struct {
  int tag, sample, n ;
  char types[MAXTAGDEPTH] ;
  long int x[MAXTAGDEPTH] ;
  float z[MAXTAGDEPTH] ;
} HIT ;


#ifdef RESTRICT 
#define Restrict restrict
#else
#define Restrict 
#endif

/*************************************************************************************/

static int tsfSampleOrder (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va ;
  const HIT *vp = (const HIT *)vb ;
  int n ;

  n = up->sample - vp->sample ; if (n) return n ;
  n = up->tag - vp->tag ; if (n) return n ;

  return 0 ;
} /* tsfSampleOrder */

/*************************************************************************************/

static int tsfTagOrder (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va ;
  const HIT *vp = (const HIT *)vb ;
  int n ;

  n = up->tag - vp->tag ; if (n) return n ;
  n = up->sample - vp->sample ; if (n) return n ;

  return 0 ;
} /* tsfTagOrder */

/*************************************************************************************/
/*************************************************************************************/

static void tsfSampleSelectInit (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle() ;
  char *cq, *cp = strnew (tsf->sampleSelectList, h) ;
 
  tsf->sampleSelectOnly = TRUE ;
  tsf->sampleSelectDict = dictHandleCreate (256, tsf->h) ;
  while (cp)
    {
      cq = strchr (cp, ',') ;
      if (cq) *cq++ = 0 ;
      if (!cp[0])
	{
	  fprintf (stderr, "parameter --sampleSelect cannot be parsed, try --help,  expecting comma delimited list a,b,c and so on, received %s"
		 , tsf->sampleSelectList) ;
	  exit (1) ;
	}

      if (strchr (cp, '*'))
	tsf->sampleSelectOnly = FALSE ;
      else
	{
	  dictAdd (tsf->sampleSelectDict, cp, 0) ;
	  dictAdd (tsf->sampleDict, cp, 0) ;
	}
      cp = cq ;
    }
  ac_free (h) ;

} /* tsfSampleSelectInit */

/*************************************************************************************/

static void tsfTagSelectInit (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle() ;
  char *cq, *cp = strnew (tsf->tagSelectList, h) ;
 
  tsf->tagSelectOnly = TRUE ;
  tsf->tagSelectDict = dictHandleCreate (256, tsf->h) ;
  while (cp)
    {
      cq = strchr (cp, ',') ;
      if (cq) *cq++ = 0 ;
      if (!cp[0])
	{
	  fprintf (stderr, "// --tagSelect ERROR, try --help,  expecting comma delimited list a,b,c and so on, received %s"
		 , tsf->tagSelectList) ;
	  exit (1) ;
	}

      if (strchr (cp, '*'))
	tsf->tagSelectOnly = FALSE ;
      else
	{
	  dictAdd (tsf->tagSelectDict, cp, 0) ;
	  dictAdd (tsf->tagDict, cp, 0) ;
	}
      cp = cq ;
    }
  ac_free (h) ;

} /* tsfTagSelectInit */

/*************************************************************************************/

static void tsfTagRenameInit (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle() ;
  char *cr, *cq, *cp = strnew (tsf->tagRenameList, h) ;
 
  tsf->tagOldNameDict = dictHandleCreate (256, tsf->h) ;
  tsf->tagNewNameDict = dictHandleCreate (256, tsf->h) ;
  tsf->tagRenameKs = keySetHandleCreate (tsf->h) ;
  while (cp)
    {
      int old, new ;

      cq = strchr (cp, ',') ;
      if (cq) *cq++ = 0 ;
      
      cr = strchr (cp, ':') ;
      if (!cr || ! cr[1] || cp[0] == ':')
	{
	  fprintf (stderr, "// --tagRename ERROR, try --help,  expecting a1:b1,a2:b2 found %s in string %s"
		 , cp, tsf->tagRenameList) ;
	  exit (1) ;
	}
      *cr++ = 0 ;
      if (! strncmp (cr, "_P_", 3)) /* do not parse percentages */
	{
	  fprintf (stderr, "// --tagRename ERROR, _P_ is reserved, you cannot rename %s to %s", cp, cr)  ;
	  exit (1) ;
	}

      dictAdd (tsf->tagOldNameDict, cp, &old) ;
      dictAdd (tsf->tagNewNameDict, cr, &new) ;
      keySet (tsf->tagRenameKs, old) = new ;

      cp = cq ;
    }
  ac_free (h) ;
} /* tsfTagSelectInit */

/*************************************************************************************/

static void tsfInit (TSF *tsf)
{
  tsf->tagDict = dictHandleCreate (256, tsf->h) ;
  tsf->sampleDict = dictHandleCreate (10000, tsf->h) ;
  tsf->hits = bigArrayHandleCreate (100000, HIT, tsf->h) ;
  tsf->tableCaptionDict = dictHandleCreate (256, tsf->h) ;
  tsf->columnCaptionDict = dictHandleCreate (256, tsf->h) ;
  tsf->valDict = dictHandleCreate (256, tsf->h) ;
  tsf->tagDepth = keySetHandleCreate (tsf->h) ;
  if (tsf->tagRenameList)
    tsfTagRenameInit (tsf)  ;
  if (tsf->tagSelectList)
    tsfTagSelectInit (tsf)  ;
  if (tsf->sampleSelectList)
    tsfSampleSelectInit (tsf)  ;
  return  ;
} /* tsfInit */

/*************************************************************************************/

static long int tsfParseTable (TSF *tsf, ACEIN ai)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, kk, sample = 0 ;
  const char *ccp, *ccq ;
  char cutter ;
  HIT *hit ;
  /*
  KEYSET tagDepth = tsf->tagDepth ; 
  DICT *tagDict = tsf->tagDict ;
  */
  DICT *sampleDict = tsf->sampleDict ;
  BigArray hits = tsf->hits ;
  long int iMax = bigArrayMax (hits) ;
  int sampleDepth = 1 ;
  vTXT sampleTxt = vtxtHandleCreate (h) ;
  KEYSET tagsKs = keySetHandleCreate (h) ;
  KEYSET nKs = keySetHandleCreate (h) ;
  KEYSET posKs = keySetHandleCreate (h) ;

  if (tsf->sample)
    {
      ccp = tsf->sample ;
      ccq = strstr (ccp, "..") ;
      sampleDepth = 1 ;
      while (ccq)
	{
	  sampleDepth++ ;
	  ccq = strstr (ccp + 2, "..") ;
	}
      tsf->sampleDepth = sampleDepth ;	
      dictAdd (sampleDict, ccp, &sample) ;
    }

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      ccp = aceInWordCut (ai, "\t", &cutter) ;
      if (! ccp)
	continue ;
      if (!strncmp (ccp, "###", 3))
	continue ;
      if (!strncmp (ccp, "##", 3))
	{
	  dictAdd (tsf->tableCaptionDict, aceInPos (ai), 0) ;
	  continue ;
	}
      if (!strncmp (ccp, "#", 1)) /* establish tag names */
	{

	  continue ;
	}

      /* parse the samples */ 
      vtxtClear (sampleTxt) ;
      for (i = 0 ; i < sampleDepth ; i++)
	{
	  if (i>0)
	    vtxtPrintf (sampleTxt, "..") ;
	  vtxtPrintf (sampleTxt, "%s", ccp) ;
	  ccp = aceInWordCut (ai, "\t", &cutter) ;
	}
      
      if (tsf->sampleSelectOnly && tsf->sampleSelectDict 
	  && ! dictFind (tsf->sampleSelectDict, ccp, 0)
	  )
	continue ;
      dictAdd (sampleDict, vtxtPtr (sampleTxt), &sample) ;
      
      kk = 0 ;
      while (aceInStep (ai, '\t')) 
	{
	  kk++ ;
	  ccp = aceInWordCut (ai, "\t", &cutter) ;
	  if (ccp)
	    {
	      char cc ;
	      long int x ;
	      if (sscanf (ccp, "%ld%c", &x, &cc) == 1)
		{
		  int pos = keySet (posKs, kk) ;
		  
		  hit = bigArrayp (hits, iMax++, HIT) ;
		  hit->tag = keySet (tagsKs, kk) ;
		  hit->sample = sample ;
		  hit->n = keySet (nKs, kk) ;
		  hit->x[pos] = x ;
		}
	    }
	}
    }
  ac_free (h) ;
  return 0 ;
}/*  tsfParseTable */

/*************************************************************************************/

static long int tsfParseTsf (TSF *tsf, ACEIN ai)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, n, tag, sample = 0, x, line = 0 ;
  float z ;
  const char *ccp, *ccq ;
  HIT *hit ;
  KEYSET tagDepth = tsf->tagDepth ; 
  DICT *valDict = tsf->valDict ;
  DICT *tagDict = tsf->tagDict ;
  DICT *sampleDict = tsf->sampleDict ;
  BigArray hits = tsf->hits ;
  long int iMax = bigArrayMax (hits) ;
  int sampleDepth ;

  if (tsf->sample)
    {
      ccp = tsf->sample ;
      ccq = strstr (ccp, "..") ;
      sampleDepth = 1 ;
      while (ccq)
	{
	  sampleDepth++ ;
	  ccq = strstr (ccp + 2, "..") ;
	}
      tsf->sampleDepth = sampleDepth ;	
      dictAdd (sampleDict, ccp, &sample) ;
    }

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      line++ ;
      if (! ccp)
	continue ;
      if (!strncmp (ccp, "###", 3))
	continue ;
      if (!strncmp (ccp, "##", 3))
	{
	  dictAdd (tsf->tableCaptionDict, aceInPos (ai), 0) ;
	  continue ;
	}
      if (!strncmp (ccp, "#", 1))
	{
	  dictAdd (tsf->columnCaptionDict, aceInPos (ai), 0) ;
	  continue ;
	}
      if (! strncmp (ccp, "_P_", 3)) /* do not parse percentages */
	continue ; 
      if (tsf->tagOldNameDict)
	{
	  int old, new ;
	  if (dictFind (tsf->tagOldNameDict, ccp, &old))
	    {
	      new = keySet (tsf->tagRenameKs, old) ;
	      ccp = dictName (tsf->tagNewNameDict, new) ;
	    }
	}
      if (tsf->tagSelectOnly && tsf->tagSelectDict 
	  && ! dictFind (tsf->tagSelectDict, ccp, 0)
	  )
	continue ;
      if (1)
	{
	  int mx = 0 ;
	  messAllocMaxStatus (&mx) ; 
	  if (mx > 50000)
	    invokeDebugger () ;	    
	}
     dictAdd (tagDict, ccp, &tag) ;
     aceInStep (ai, '\t') ;
     ccp = aceInWord (ai) ;
     if (! ccp || *ccp == '#')
       continue ;
     if (tsf->sampleSelectOnly && tsf->sampleSelectDict 
	 && ! dictFind (tsf->sampleSelectDict, ccp, 0)
	 )
       continue ;
     if (! tsf->sample)
       {
	 ccq = strstr (ccp, "..") ;
	 sampleDepth = 1 ;
	  while (ccq)
	    {
	      sampleDepth++ ;
	      ccq = strstr (ccq + 2, "..") ;
	    }
	  if (sampleDepth > tsf->sampleDepth)
	    tsf->sampleDepth = sampleDepth ;
	  if (sampleDict)  /* sample_name , ignore in merge case */
	    dictAdd (sampleDict, ccp, &sample) ;
       }
     else
       dictAdd (sampleDict, tsf->sample, &sample) ;
     aceInStep (ai, '\t') ;
     ccp = aceInWord (ai) ;
     if (! ccp)
       continue ;
     n = strlen (ccp) ;
     if (n < 0 || n >= MAXTAGDEPTH)
       messcrash ("Overflow, number of fields = %d > %d at line %d of file %s.\nWe do not expect more than %d value following a tag, please edit the code or reformat  the data"
		  , MAXTAGDEPTH
		  , aceInStreamLine (ai)
		  , aceInFileName (ai)
		  , n) ;
     hit = bigArrayp (hits, iMax++, HIT) ;
     hit->tag = tag ;
     hit->sample = sample ;
     memcpy (hit->types, ccp, n+1) ;
     
     for (i = 0 ; i < n ; i++)
       {
	 if (! aceInStep (ai, '\t'))
	   break ;
	  hit->n = i + 1 ;
	  switch (hit->types[i])
	    {
	    case 'i':
	      if (aceInInt (ai, &x))
		hit->x[i] = x ;
	      else
		messcrash ("Non integer number in line %d column %d file %s\n"
			   , 4 + i
			   , aceInStreamLine (ai)
			   , aceInFileName (ai)
			   ) ;
	      break ;
	    case 'f':
	      if (aceInFloat (ai, &z))
		hit->z[i] = z ;
	      else
		messcrash ("Non float number in line %d column %d file %s\n"
			   , 4 + i
			   , aceInStreamLine (ai)
			   , aceInFileName (ai)
			   ) ;
	      break ;
	    default :
	      ccp = aceInWord (ai) ;
	      if (ccp)
		{
		  int k = 0 ;
		  dictAdd (valDict, ccp, &k) ; 
		  hit->x[i] = k ; 
		}
	      else
		messcrash ("Non word in line %d column %d file %s\n"
			   , 4 + i
			   , aceInStreamLine (ai)
			   , aceInFileName (ai)
			   ) ;
	      break ;
	    }
       }
     if (i > keySet (tagDepth, tag))
       keySet (tagDepth, tag) = i ;
     
    }
  ac_free (h) ;

  return bigArrayMax (hits) ;
} /*  tsfParseTsf */

/*************************************************************************************/

static void tsfParse (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle () ;
  
  if (tsf->inFileList)
    {
      ACEIN ai = 0 ;
      int nn = 0 ;
      const char *error = 0 ;

      /* check that all files exist */
      if (! aceInCheckList (tsf->inFileList, &error, h))
	messcrash ("Bad parameter -i, %s"
		       ,  error
		       ) ;	
 
      /* parse a validated list of files */
      while ((ai = aceInCreateFromList (ai, &nn, tsf->inFileList, tsf->gzi, h)))
	{
	  if (tsf->isTableIn)
	    tsfParseTable (tsf, ai) ;
	  else
	    tsfParseTsf (tsf, ai) ;
	}
     }
  else
    {
      ACEIN ai = aceInCreate (0, tsf->gzi, h) ;
      if (tsf->isTableIn)
	tsfParseTable (tsf, ai) ;
      else
	tsfParseTsf (tsf, ai) ;
    }

  ac_free (h) ;

  return  ;
} /* tsfParse */

/*************************************************************************************/
static KEYSET tsfMyLines = 0 ;
static KEYSET tsfMyCols = 0 ;

static int tsfMyOrder (const void *va, const void *vb)
{
  const HIT *up = (const HIT *)va ;
  const HIT *vp = (const HIT *)vb ;
  int n ;

  n = keySet (tsfMyLines, up->tag -1) - keySet (tsfMyLines, vp->tag -1) ; if (n) return n ;
  n = keySet (tsfMyCols, up->sample -1) - keySet (tsfMyCols, vp->sample -1) ; if (n) return n ;

  return 0 ;
} /* tsfTagOrder */

/*************************************************************************************/

static void tsfChronoOrder (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle () ;
  HIT *hit ;
  int ii, iiMax = bigArrayMax (tsf->hits) ;
  int i, tMax = dictMax (tsf->tagDict) ;
  int j, sMax = dictMax (tsf->sampleDict) ;
  KEYSET cols = keySetHandleCreate (h) ;
  KEYSET lines = keySetHandleCreate (h) ;
  Array bb = 0 ;
  BOOL debug = FALSE ;

  for (i = 0 ; i < tMax ; i++)
    keySet (lines, i) = i ;
  for (j = 0 ; j < sMax ; j++)
    keySet (cols, j) = j ;
  
  if (tsf->chronoOrder)
    {
      bb = arrayHandleCreate (tMax * sMax, double, h) ; 
      array (bb, tMax * sMax - 1, double) = 0 ;
      for (ii = 0, hit = bigArrp (tsf->hits, 0, HIT) ; ii < iiMax ; hit++, ii++)
	{
	  double z  ;
	  
	  i = hit->tag - 1 ;
	  j = hit->sample - 1 ;
	  z = hit->z[0] + hit->x[0] ;
	  array (bb, i * sMax + j, double) = z ;
	}
      topoChronoOrder (bb, lines, cols) ;
    }

  tsfMyLines = arrayHandleCopy (lines, 0) ;
  tsfMyCols = arrayHandleCopy (cols, 0) ;

  for (i = 0 ; i < tMax ; i++)
    keySet (tsfMyLines, keySet (lines, i)) = i ;
  for (j = 0 ; j < sMax ; j++)
    keySet (tsfMyCols, keySet (cols, j)) = j ;

  if (debug)
    for (ii = 0, hit = bigArrp (tsf->hits, 0, HIT) ; ii < iiMax ; hit++, ii++)
      {
	float z  ;
	
	i = hit->tag ;
	j = hit->sample ;
	z = hit->z[0] + hit->x[0] ;
	printf("%d\t%d\t%.1f\n",  i, j, z) ; 
      }

  if (tsf->chronoOrder)
    bigArraySort (tsf->hits, tsfMyOrder) ;

  if (debug)
    for (ii = 0, hit = bigArrp (tsf->hits, 0, HIT) ; ii < iiMax ; hit++, ii++)
      {
	float z  ;
	
	i = hit->tag ;
	j = hit->sample ;
	z = hit->z[0] + hit->x[0] ;
	printf("%d\t%d\t%.1f\n",  i, j, z) ; 
      }



  ac_free (h) ;
  return ;
}

/*************************************************************************************/
/*************************************************************************************/

static long int tsfMerge (TSF *tsf)
{
  BigArray hits = tsf->hits ;
  long int ii, jj, iMax = bigArrayMax (hits) ;
  HIT *hit, *hit2 ;

  /* sort and cumulate the values */
  bigArraySort (hits, tsfTagOrder) ;
  for (ii = 0, hit = bigArrp (hits, 0, HIT) ; ii < iMax ; ii++, hit++)
    {
      if ( !hit->tag)
	continue ;

      for (jj = ii + 1, hit2 = hit+1 ; jj < iMax && hit2->tag == hit->tag && hit2->sample == hit->sample ; jj++, hit2++)
	{
	  int n, n1, n2, i ;
	  
	  n1 = hit->n ;
	  n2 = hit->n ;
	  n = n1 < n2 ? n1 : n2 ;
	  if (strncmp (hit->types, hit2->types, n))
	    break ;
	  for (i = 0 ; i < n ; i++)
	    if (hit->types[i] != i && hit->x[i] && hit2->x[i] && hit->x[i] != hit2->x[i])
	      break ;
	  /* ok, all text fileds are equal: merge */	  
	  if (n1 < n2)
	    memcpy (hit->types, hit2->types, n2) ;
	  hit2->tag = 0 ;
	  n = n1 > n2 ? n1 : n2 ;
	  hit->n = n ;
	  for (i = 0 ; i < n ; i++)
	    switch (hit->types[i])
	      {
	      case 'i':
		hit->x[i] += hit2->x[i] ;
		break ;
	      default:
		if (! hit->x[i])
		  hit->x[i] = hit2->x[i] ;
	      }
	}
    } 
  /* clean up */
  for (ii = jj = 0, hit2 = hit = bigArrp (hits, ii, HIT) ; ii < iMax ; ii++, hit++)
    {
      if ( !hit->tag)
	continue ;
      if (hit != hit2)
	*hit2 = *hit ;
      jj++ ; hit2++ ; 
    } 
  bigArrayMax (hits) = jj ;
  return bigArrayMax (hits) ;
} /* tsfMerge */

/*************************************************************************************/
/*************************************************************************************/

static void tsfExportCaption (TSF *tsf, ACEOUT ao)
{
  int i, iMax = dictMax (tsf->tableCaptionDict) ;
  
  aceOutDate (ao, "###"
	      ,  tsf->title ? tsf->title : "No title"
	      ) ;
  
  for (i = 1 ; i <= iMax ; i++)
    aceOutf (ao, "## %s\n", dictName (tsf->tableCaptionDict, i)) ;

return ;
} /* tsfExportCaption */

/*************************************************************************************/

static void tsfExportTsf (TSF *tsf, ACEOUT ao)
{
  long int ii, iMax = bigArrayMax (tsf->hits) ;
  HIT * Restrict hit = iMax ? bigArrp (tsf->hits, 0, HIT) : 0 ;
  DICT * Restrict valDict = tsf->valDict ;
  DICT * Restrict tagDict = tsf->tagDict ;
  DICT * Restrict sampleDict = tsf->sampleDict ;
  int k, hasP = 0 ;

  for (ii = 0 ; ii < iMax ; ii++)
    {
      int i, n ;
      HIT hh = hit[ii] ;
      const char *tagName = dictName (tagDict, hh.tag) ;

      aceOutf (ao, "%s", tagName) ;
      aceOutf (ao, "\t%s", dictName (sampleDict, hh.sample)) ;
      n = hh.n ;
      aceOutf (ao, "\t%s", hh.types) ; 
      for (i = 0 ; i < n ; i++)
	switch (hit->types[i])
	  {
	  case 'i':
	    aceOutf (ao, "\t%ld", hh.x[i]) ;
	    break ; 
	  case 'f':
	    aceOutf (ao, "\t%f", hh.z[i]) ;
	    break ; 
	  default:
	    k = hh.x[i] ;
	    aceOutf (ao, "\t%s", dictName (valDict, k)) ;
	    break ;
	  }
      aceOutf (ao, "\n") ;

      if (strncmp (tagName, "P_", 2))
	hasP = ii + 1 ;
    }
  
  /* export the percentages */
  for (ii = 0 ; ii < hasP ; ii++)
    {
      int i, n ;
      double z ;
      HIT hh = hit[ii] ;
      const char *tagName = dictName (tagDict, hh.tag) ;

      if (strncmp (tagName, "P_", 2))
	continue ;

      /* Export the percentages */
      aceOutf (ao, "_%s", tagName) ;
      aceOutf (ao, "\t%s", dictName (sampleDict, hh.sample)) ;
      n = hh.n ;
      aceOutf (ao, "\t%s\t100%%", hh.types) ; 
      z = hh.x[0] ;
      if (z == 0) z = 100 ;
      z = 100.0/z ;
      for (i = 1 ; i < n ; i++)
	switch (hit->types[i])
	  {
	  case 'i':
	    aceOut (ao, "\t") ;
	    aceOutPercent (ao, z * hh.x[i]) ;
	    aceOut (ao, "%") ;
	    break ; 
	  case 'f':
	    aceOut (ao, "\t") ;
	    aceOutPercent (ao, z * hh.z[i]) ;
	    aceOut (ao, "%") ;
	    break ; 
	  default:
	    k = hh.x[i] ;
	    aceOutf (ao, "\t%s", dictName (valDict, k)) ;
	    break ;
	  }
      aceOutf (ao, "\n") ;	
    }

  return ;
} /* tsfExportTsf */

/*************************************************************************************/
/* check if the kk first fileds are identical with previous line */
static int tsfSame (int same, int kk, vTXT txt, vTXT oldTxt)
{
  if (same + kk < 0 ) /* we need to study the line */
    {
      const char *cp = vtxtPtr (txt) ;
      const char *cq = vtxtPtr (oldTxt) ;
      
      if (cp && cq)
	{
	  int n = strlen (cp) ;
	  if (strncmp (cp, cq, n))
	    return kk ;
	}
    }
  return same ;
} /* tsfSame */

/*************************************************************************************/

static void tsfExportTable (TSF *tsf, ACEOUT ao)
{
  AC_HANDLE h = ac_new_handle () ;
  long int ii, iiMax = bigArrayMax (tsf->hits) ;
  int i, jj, kk, k1 ;
  HIT * Restrict hit = iiMax ? bigArrp (tsf->hits, 0, HIT) : 0 ;
  DICT * Restrict tagDict = tsf->tagDict ;
  DICT * Restrict valDict = tsf->valDict ;
  DICT * Restrict sampleDict = tsf->sampleDict ;
  int oldTag = 0 ;
  int sampleMax = dictMax (sampleDict) ;
  int pass ;
  KEYSET skip = keySetHandleCreate (h) ;
  
  for (i = 0 ; i < tsf->skip1 ; i++)
    keySet (skip, i)++ ;
  for (i = 0 ; i < tsf->skip2 ; i++)
    keySet (skip, i) += 2 ;
  for (i = 0 ; i < tsf->skip3 ; i++)
    keySet (skip, i) += 3 ;
  
  aceOutf (ao, "# Tag") ;
  for (i = 1 ; i <= sampleMax ; i++)
    {
      int i1 = keySet (tsfMyCols, i-1) + 1 ;
      aceOutf (ao, "\t%s", dictName (sampleDict, i1)) ;
    }
  
  for (pass = 0 ; pass < 2 ; pass++)
    for (ii = 0 ; ii < iiMax ; ii++)
      {
	int tag = hit[ii].tag ;
	const char *tagName ;
	
	if (!tag || tag == oldTag)
	  continue ;
	oldTag = tag ;
	tagName = dictName (tagDict, tag) ;
	if (pass == 0 && ! strncmp (tagName, "P_", 2))
	  continue ;
	if (pass == 1 && strncmp (tagName, "P_", 2))
	  continue ;
	
	aceOutf (ao, "\n%s", dictName (tagDict, tag)) ;
	
	for (kk = 1 ; kk <= sampleMax ; kk++)
	  {
	    BOOL ok = FALSE ;
	    
	    k1 = keySet (tsfMyCols, kk-1) + 1 ;
	    for (jj = ii ; jj < iiMax && hit[jj].tag == tag ; jj++)
	      if (hit[jj].sample == k1)
		{
		  int i, k, n = hit->n ;
		  char* sep = "\t" ;
		  for (i = 0 ; i < n ; i++)
		    {
		      aceOut (ao, sep) ;
		      switch (hit->types[i])
			{
			case 'i':
			  aceOutf (ao, "%ld", hit[jj].x[i]) ;
			  break ; 
			case 'f':
			  aceOutf (ao, "%0f", hit[jj].z[i]) ;
			  break ; 
			default:
			  k = hit[jj].x[i] ;
			  aceOutf (ao, "%s", dictName (valDict, k)) ;
			  break ;
			}
		      sep = "," ;
		    }
		  ok = TRUE ;
		  break ;
		}
	    if (! ok)
	      {
		aceOut (ao, "\t") ;
	      }
	  }
      }
# ifdef JUNK
  /* to be used for pass == 1 */
  /* export the tag percentages */
  for (i = 1 ; i < hasP ; i++)
    {  
      int jMax = keySet (tsf->tagDepth, i) ;
      const char *tagName = dictName (tagDict, i) ;
      
      if (strncmp (tagName, "P_", 2))
	continue ;
      
      for (j = 0 ; j < jMax ; j++)
	{
	  kMax++ ;
	  aceOutf (ao, "\t%%%s", tagName + 2) ;
	  if (jMax > 1)
	    aceOutf (ao, ":%d", j+1) ;
	}
    }

  // la suite je sais pas a quoi cela sert (2020-08-07 */

  for (ii = kk = sample = 0 ; ii < iMax ; ii++)
    {
      int i ;
      int same = - keySetMax (skip) ;
      HIT hh = hit[ii] ;
      
      { /* switch */
	vTXT v = txt ; txt = oldTxt ; oldTxt = v ; 
      } 
      vtxtClear (txt) ;
      if (hh.sample != sample)
	{
	  if (sample)  /* complete the previous line */
	    for ( ; kk < kMax ; kk++)
	      vtxtPrintf (txt, "\t") ;
	  kk = 0 ;
	  sample = hh.sample ;
	}
    
      if (tsf->sampleDepth < 2)
	{
	  vtxtPrintf (txt, "\n%s", dictName (sampleDict, sample)) ;
	  kk++ ;
	  same = tsfSame (same, kk, txt, oldTxt) ;
	}
      else
	{
	  char *cq, *cp = strnew (dictName (sampleDict, sample), 0) ;
	  char *sep = "\n" ;
	  
	  for (i = 0 ; i < tsf->sampleDepth ; i++)
	    {
	      vtxtPrint(txt, sep) ;
	      if (cp)
		{
		  cq = strstr (cp, "..") ;
		  if (cq) 
		    *cq = 0 ;
		  vtxtPrintf (txt, "%s", cp) ;
		  kk++ ;
		  same = tsfSame (same, kk, txt, oldTxt) ;
		}
	      sep = "\t" ;
	      if (cq) cp = cq + 2 ;
	      else cp = 0 ;
	    }
	}
      
      
      for (i = 0 ; i <= tagMax ; i++)
	{
	  int k, jMax = keySet (tsf->tagDepth, i) ;
	  for (j = 0 ; j < jMax ; j++)
	    switch (hit->types[i])
	      {
	      }
	  same = tsfSame (same, kk, txt, oldTxt) ;
	}
      
      /* export the percentages */
      for (i = 1 ; i < hasP ; i++)
	{  
	  double z ;
	  int jMax = keySet (tsf->tagDepth, i) ;
	  const char *tagName = dictName (tagDict, i) ;
	  
	  if (strncmp (tagName, "P_", 2))
	    continue ;
	  
	  kk++ ; 
	  vtxtPrintf (txt, "\t100%%") ; 
	  z = (hit->types[0] == 'i' ? hit->x[0] : hit->z[0]) ;
	  if (z == 0) z = 100 ;
	  z = 100.0/z ;
	  same = tsfSame (same, kk, txt, oldTxt) ;
	  
	  for (j = 1 ; j < jMax ; j++)
	    {
	      float zj = (hit->types[i] == 'i' ? hit->x[i] : hit->z[i]) ;
	      kk++ ; 
	      vtxtPrintf (txt, "\t") ;
	      vtxtPercent (txt, z * zj) ;
	      vtxtPrintf (txt, "%%") ;
	      same = tsfSame (same, kk, txt, oldTxt) ;
	    }
	}
      
      if (same > 0)
	for (i = 0 ; i < keySet (skip, same) ; i++)
	  aceOut (ao, vtxtPtr (txt)) ;
      if (vtxtPtr (txt))
	aceOut (ao, vtxtPtr (txt)) ;
    }
#endif
  aceOut (ao, "\n") ;
  ac_free (h) ;

  return ;
} /* tsfExportTable */

/*************************************************************************************/

static void tsfTransposeTable (TSF *tsf, ACEOUT ao)
{
  if (! tsf->transpose)
    return  tsfExportTable (tsf, ao) ;
  else
    { 
      AC_HANDLE h = ac_new_handle () ;
      Stack s = stackHandleCreate (10000, h) ;
      Stack s2 = stackHandleCreate (10000, h) ;
      ACEOUT bo = aceOutCreateToStack (s, h) ;
      ACEIN bi ;
      int ii, jj, iMax, jMax ;
      Array kss = arrayHandleCreate (10000, KEYSET, h) ;
      KEYSET ks ;
      char *sep, *word, cutter ;

      /* export to bo the direct table */
      tsfExportTable (tsf, bo) ;
      ac_free (bo) ;

      /* parse the formated table*/
      bi =  aceInCreateFromText (stackText (s, 0), 0, h) ;
      aceInSpecial (bi, "\n") ;
      iMax = jMax = 0 ;
      while (aceInCard (bi))
	{
	  ks = array (kss, iMax++, KEYSET) = keySetHandleCreate (h) ; ; 
	  jj = -1 ;

	  while (1)
	    {
	      int k = stackMark (s2) ;
	      cutter = 0 ;
	      aceInStep (bi, '\t') ;
	      word = aceInWordCut (bi, "\t", &cutter) ;
	      jj++ ;
	      if (0) fprintf (stderr, "%s :: ii=%d jj=%d k=%d\n", word, iMax,jj,k) ;
	      keySet (ks, jj) = k ;
	      if (word)
		pushText (s2, word) ;
	      if (jj + 1 > jMax)
		jMax = jj + 1 ;
	      if (cutter != '\t')
		break ;
	    }
	}
      ac_free (bi) ;

      /* export a transposed table */
      sep = "" ;
      for (jj = 0 ; jj < jMax ; jj++)
	{
	  for (ii = 0 ; ii < iMax ; ii++)
	    {
	      ks = array (kss, ii, KEYSET) ; 
	      if (ks && jj < keySetMax (ks))
		word = stackText (s2, keySet (ks, jj)) ;
	      else
		word = "" ;
	      aceOutf (ao, "%s%s", sep, word) ;
	      sep = "\t" ;
	    }
	  sep = "\n" ;
	}
      aceOutf (ao, "%s", sep) ;
  
      ac_free (h) ;
    }
  return ;
}  /* tsfTransposeTable */

/*************************************************************************************/

static void tsfExport (TSF *tsf)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (tsf->outFileName, tsf->isTableOut ? ".txt" : ".tsf", tsf->gzo, h) ;

  if (ao)
    {
      tsfExportCaption (tsf, ao) ;
      if (tsf->isTableOut)
	{
	  bigArraySort (tsf->hits, tsfTagOrder) ;
	  tsfChronoOrder (tsf) ; /* used to position tsfMyLines/Cols */
	  tsfTransposeTable (tsf, ao) ;
	}
      else
	{
	  bigArraySort (tsf->hits, tsfTagOrder) ;
	  tsfExportTsf (tsf, ao) ;
	}
    }
  ac_free (h) ;
} /* tsfExport */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/
/* sample -sample NA12878MOD -laneList tmp/TSNP/NA12878MOD/LaneList -t 20 -target_fasta TARGET/CHROMS/hs.chrom_20.fasta.gz -t1 25000001 -t2 30010000 -minSnpFrequency 18 -minSnpCover 10 -minSnpCount 4 -target_class Z_genome -o tata -maxLanes 4
 */

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// tsf: A generic tool to manipulate numeric data and tables\n"
	   "//   Contact: Jean  Thierry-Mieg, NCBI  mieg@ncbi.nlm.nih.gov\n"
	   "//   This code is public and distributed without restrictions or guarantees\n"
	   "// Usage: tsf [parameters] [-i <file_list>] [-o <output>]\n"
	   "//      try: -h --help -v --version\n"
	   "// Examples:\n"
	   "//     tsf -i f1.tsf,f2.tsf\n"
	   "//       Merge several .tsf files\n"
	   "//     tsf -i f1.tsf,f2.tsf --group my_group_of_samples\n"
	   "//       For each tag, cumul counts of all samples\n"
	   "//     tsf -i f1.tsf -O table -o my_table\n"
	   "//       Import a tsf file, export a table called my_table.txt\n"
	   "//     tsf -I table -i f1.txt,f2.txt\n"
	   "//       Merge several tables into a probably larger table\n"
	   "//       adding up the values with corresponding captions\n"
	   "//     tsf -I table --transpose --lineSelect a,b,c --colSelect A,B\n"
	   "//       Manipulate a generic tab delimited table\n"
	   "// INPUT:\n"
	   "//   -i <file_list> : short form of\n"
	   "//   --in_file_list <file_list>\n"
	   "//     A comma delimited list of input files, for example\n"
	   "//       -i f1,f2,f3\n"
	   "//     otherwise, read the data from stdin\n"
	   "//   --gzi\n"
	   "//     gunzip the input\n"
	   "//       This is automatic for all input files called .gz\n"
	   "//       The option is useful when reading stdin, or a file not named .gz\n"
	   "// INPUT FORMAT:\n"
	   "//   -I [tsf|table]: short form of\n"
	   "//   --In_format [tsf|table] : default tsf\n"
	   "//     The input files are provided in tsf or in table format (defined below)\n"
	   "// OUPUT FORMAT:\n"
	   "//   -O [tsf|table]: short form of\n"
	   "//   --Out_format [tsf|table] : default tsf\n"
	   "//     The output files are exported in tsf or in table format (defined below)\n"
	   "// OUPUT:\n"
	   "//   -o <file_name>: short form of\n"
	   "//   --out <file_name>\n"
  	   "//   --gzo : the output files will be gzipped\n"
	   "//     Depending on -O and -gzo, export to file_name[.tsf|.txt][.gz]\n"
	   "//     If there is no -o paramater, export to stdout\n"
	   "// TSF FORMAT:  -I tsf -O [tsf|table]\n"
	   "//     The tsf format is a very simple autocumulative format\n"
	   "//        tag sample format values\n"
	   "//     Where the format announces the number/types of values. Example\n"
	   "//        tag sample 5 x1 x2 x3 ... xn\n"
	   "//     Where n (n<=24) is the number of xi fields on this line and all\n"
	   "//     xi are long signed integral numbers: ...,-2,-1,0,1,2,...\n"
	   "//     The formats are i (integer), f (float), t (text)\n"
	   "//     using a capital (IFT) allows mutivalued fields\n"
	   "//     a number can be used to repeat a field: 5t3i == tttttiii\n"
	   "//     a number with no letter defaults to integer: 3 == 3i == iii\n"
	   "//   -i file_list\n"
	   "//     Merge .tsf files, for example tsf files exported by this program.\n"
	   "//     For example, if a large experiment is split in data batches and\n"
	   "//     each batch produces a .tsf file, merging the files will provide\n"
	   "//     the global counts of the experiment. Example:"
	   "//     Input:\n"
	   "//        Tag1 sample_1 2  200 124\n"
	   "//        Tag1 sample_1 2  60  41\n"
	   "//     Command: tsf \n"
	   "//     Output:\n"
	   "//        Counts Tag1 sample_1 2 260 165\n"
	   "// ACTION\n"
	   "//       When several lines in single file or in separate files refer to\n"
	   "//     the same tag and same sample, one must decide how to merge the values\n"
	   "//       t : single valued text, the last value wins.\n"
	   "//           a  a => a\n"
	   "//           a  b => b\n"
	   "//       T : muti valued text: join removing duplicates\n"
	   "//           a  a => a\n"
	   "//           a  b => a,b\n"
	   "//       i : [--sum: default] add corresponding values.\n"
	   "//           17  17 => 34\n"
	   "//           17  31 => 48\n"
	   "//       i : [--replace] last value wins\n"
	   "//           17  17 => 17\n"
	   "//           17  31 => 31\n"
	   "//       I : muti valued integer: join removing duplicates\n"
	   "//           4,5  4 => 4,5\n"
	   "//           4,7  5 => 4,5,7\n"
	   "//     For single valued {i,f} numbers, one can require a calculation\n"
	   "//   --sum  [default] (add up single valued integers in multiple lines of one or several files)\n"
	   "//          4 9 => 13\n"
	   "//   --replace : do not add single valued integers\n"
	   "//          4 9 => 9\n"
	   "//   --min : take the minimum of correspomding single valued integers\n"
	   "//          4 9 => 4\n"
	   "//   --max : take the maximum of correspomding single valued integers\n"
	   "//          4 9 => 9\n"
	   "//   --compute \"a + 3*b - c + 7\"  (only for multiple input files)\n"
	   "//          4 9 5 => 4 + 3*9 - 5 + 7 = 33\n"
	   "//   The --compute is only valid when parsing a list of files and the\n"
	   "//   variables a,b,c ... in the equation apply to file 1,2,3 ...\n"
	   "//   The calculation is applied to each single valued number of the files.\n"
	   "//   Absent values default as zero. Multi valued numbers {I,F} are merged\n"
	   "//   as in the examples above, but not computed\n"
	   "//   -g group_name] : short form of\n"
	   "//   --group group_name]\n"
	   "//     This parameters can be used to merge sets of sample into a group\n"
	   "//     and cumulate the counts in each tag.  Example:\n"
	   "//     Input:\n"
	   "//        Tag1 sample_1 2  10  25\n"
	   "//        Tag2 sample_1 2  20  41\n"
	   "//        Tag1 sample_2 2  60  30\n"
	   "//     Command: tsf --group my_group\n"
	   "//     Output:\n"
	   "//        Tag1 my_group 2 70 55\n"
	   "//        Tag2 my_group 2 20  41\n"
	   "//     Note that the order of the lines is irrelevant, they can come\n"
	   "//     from a single file or from differentt files of the -i file list\n"
	   "//   P_ and _P_ percentage tags\n"
	   "//     Tags called P_ have a special meaning\n"
	   "//     They follow the standard format and are cumulated as usual\n"
	   "//     But the first field is expected to be the denominator of the other\n"
	   "//     columns and whenever a tag P_ is encountered, the program also\n"
	   "//     exports a _P_ line with the corresponding percentages.\n"
	   "//     Please notice that on input the _P_ tags are ignored,\n"
	   "//     and recomputed from the final values in the corresponding P_ lines\n"
	   "//   -O tsf : default, the data are exported in tsf format\n"
	   "//   -O table : Transforms a tsf file into a table format\n"
	   "//     When a tag is mutivalued, the columns are called Tag:1, Tag:2, Tag:3\n"
	   "//     or by subtags, if the tag name contains : characters. For example\n"
	   "//        A:B:C name 5 110 11 22 33 44\n"
	   "//     will be exported as a table with column names A B C:1 C:2 C:3\n"
	   "//     If the input tag is called P_ the percentages are also exported. e.g.\n"
	   "//        P_A:B:C name 5 110 11 22 33 44\n"
	   "//     will be exported as columns %%A %%B %%C:1 %%C:2 %%C:3\n"
	   "//     with values 100%%  10%%  20%% 30%% 40%%\n"
	   "// SEPARATORS:\n"
	   "//  --TS <sep>: short form of\n"
	   "//  --tagSeparator : default : (column character)\n"
	   "//  --SS <sep>: short form of\n"
	   "//  --sampleSeparator : default : (column character)\n"
	   "//     allows to name each subfield of a mutivalued tag\n"
	   "//     and to mname the sample hierarchically. Example\n"
	   "//       Top_speed:Mileage Ford:model_T:1911  fi 43.3  20\n"
	   "//     The fi format announces 2 value (float, int)\n"
	   "//     namely the speed and mileage of a Ford, model_T from 1911\n"
	   "//  --VS <sep>: short form of\n"
	   "//  --valueSeparator : default , (comma character)\n"
	   "//     separates multivalued values. Example\n"
	   "//       Color Ford:Mustang:1980 T red,blue,yellow\n"
	   "//     The T format announces a possibly mutivalued text\n"
	   "//  --FS <sep>: short form of\n"
	   "//  --fieldSeparator : default \\t (tab character)\n"
	   "//      separates columns in tables and T,S,F,value fields in TSF format\n"
	   "// TABLE FORMAT: -I table -O [tsf|table]\n"
	   "//   -t  <title> : short form of\n"
	   "//   --title <title> : provide a title for the table\n"
	   "//   ### Automatically set file identifier: filename, title and date.\n"
	   "//     Lines starting with ### are ignored when reading a file\n"
	   "//   ##  Table caption. Lines starting with ## are treated as caption lines,\n"
	   "//     they are reexported as is, in the same order, ignoring duplicated lines.\n"
	   "//   ## SAMPLE A:B provide the title for the sample column(s)\n"
	   "//   #   Column caption. The most recent line starting with # is treated as a\n"
	   "//     list of tags used to interpret the following lines. This provides an easy \n"
	   "//     way to merge tables with different sets of columns, into a larger table\n"
	   "//   All other lines start with a \"sample\" name, followed by values\n"
	   "// TAG SELECTION and RENAMING:\n"
	   "//   --tagRename <list>: example  --tagRename 'c1,d1;c2,d2;c3:d3'\n"
	   "//     Renames tag c1 as d1, tag c2 as d2 an so on. This applies to the tag names\n"
	   "//     in tsf format and to the column names in table format.\n"
	   "//     Tag renaming is always applied before all other operations.\n"
	   "//   --tagSelect <list> : \n"
	   "//     example --tagSelect 't1,t2'\n"
	   "//       export only tags/columns t1,t2 in that order.\n"
	   "//     example --tagSelect 't1,t2,t3,*'\n"
	   "//       export tags/columns t1,t2,t3 then all other columns.\n"
	   "//     Otherwise, the tags/columns are exported in order of first occurence\n"
	   "// SAMPLE SELECTION and ORDERING\n"
	   "//   --sampleSelect <list> : \n"
	   "//     example --sampleSelect 's1,s2,s3'\n"
	   "//       export only lines/samples s1,s2,s3 in that order.\n"
	   "//     example --sampleSelect 's1,s2,s3,*'\n"
	   "//       export lines/samples s1,s2,s3 then all other samples.\n"
	   "//     Otherwise, the samples/lines are exported in order of first occurence\n"
	   "//   --sort \"parameters\"\n"
	   "//     Calls UNIX sort, to sort the lines of a table.\n"
	   "//     example   --sort \"-k 1,1 -k 3,3n\"\n"
	   "//     sort the lines on column 1, then numerically on column 3\n"
	   "//     This supersedes the order implied by --sampleSelect\n"
	   "// PRESENTATION: TRANSPOSITION and  line skipping\n"
	   "//   --transpose\n"
	   "//       Transpose a table, exchanging lines and columns, i.e. tags and samples.\n"
	   "//   --skip[123] <column number> : skip lines\n"
	   "//         example:  --skip2 4\n"
	   "//     export 2 blank lines each time there is a new value in columns 1 to 4\n"
	   "// HELP\n"
	   "//    -h : short form of \n"
	   "//    --help : this on line help\n"
	   ) ;

  
  if (argc > 1)
    {
      fprintf (stderr,
	       "//\n// You said: %s\n", commandBuf) ;
      fprintf (stderr,
	       "// ########## ERROR: Sorry, I do not understand the argument%s ", argc > 2 ? "s" : "") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (stderr, "%s ", argv[i]) ;
      fprintf (stderr, "// Please try : tsf -h or tsf --help\n") ;
      fprintf (stderr, "\n") ;
    }
  exit (1) ;
} /* usage */

/*************************************************************************************/

int main (int argc, const char **argv)
{
  TSF tsf ;
  AC_HANDLE h = 0 ;
  char commandBuf [4000] ;
  const char *ccp ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&tsf, 0, sizeof (TSF)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  tsf.h = h ;
 
  if (argc < 2)
    usage (commandBuf, argc, argv) ;
  if (getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "--help")||
      getCmdLineBool (&argc, argv, "-h")
      )
    usage (commandBuf, 1, argv) ;
  if (getCmdLineBool (&argc, argv, "-v") ||
      getCmdLineBool (&argc, argv, "--version")
      )
    {
      fprintf (stderr, "tsf: %s\n", VERSION) ;
      exit (0) ;
    }
  tsf.gzi = getCmdLineBool (&argc, argv, "--gzi") ;
  tsf.gzo = getCmdLineBool (&argc, argv, "--gzo") ;
  tsf.transpose = getCmdLineBool (&argc, argv, "--transpose") ;
  tsf.debug = getCmdLineBool (&argc, argv, "--debug") ;

  getCmdLineOption (&argc, argv, "-i", &(tsf.inFileList)) ;
  getCmdLineOption (&argc, argv, "--in_file_list", &(tsf.inFileList)) ;
  getCmdLineOption (&argc, argv, "-o", &(tsf.outFileName)) ;
  getCmdLineOption (&argc, argv, "--out", &(tsf.outFileName)) ;
  if (getCmdLineOption (&argc, argv, "-O", &(ccp)) || getCmdLineOption (&argc, argv, "--Out_format", &(ccp)))
    {
      if (! strcmp (ccp, "tsf")) ;
      else if (! strcmp (ccp, "table"))
	tsf.isTableOut = TRUE ;
      else
	{
	  fprintf (stderr, "Bad parameter -O %s, please try tsf --help\n", ccp) ;
	  exit (1) ;
	}
    }
  if (getCmdLineOption (&argc, argv, "-I", &(ccp)) || getCmdLineOption (&argc, argv, "--In_format", &(ccp)))
    {
      if (! strcmp (ccp, "tsf")) ;
      else if (! strcmp (ccp, "table"))
	tsf.isTableIn = TRUE ;
      else
	{
	  fprintf (stderr, "Bad parameter -I %s, please try tsf --help\n", ccp) ;
	  exit (1) ;
	}
    }

  getCmdLineOption (&argc, argv, "--tagRename", &(tsf.tagRenameList)) ;
  getCmdLineOption (&argc, argv, "--tagSelect", &(tsf.tagSelectList)) ;
  getCmdLineOption (&argc, argv, "--sampleSelect", &(tsf.sampleSelectList)) ;
  getCmdLineOption (&argc, argv, "--sort", &(tsf.sortCommand)) ;

  tsf.merge = getCmdLineBool (&argc, argv, "-m") || getCmdLineBool (&argc, argv, "--merge") ;
  tsf.replace = getCmdLineBool (&argc, argv, "--replace") ;
  tsf.min = getCmdLineBool (&argc, argv, "--min") ;
  tsf.max = getCmdLineBool (&argc, argv, "--max") ;
  tsf.sum = getCmdLineBool (&argc, argv, "--sum") ;
  tsf.chronoOrder = getCmdLineBool (&argc, argv, "--chronoOrder") ;
  getCmdLineOption (&argc, argv, "--compute", &(tsf.compute)) ;

  getCmdLineOption (&argc, argv, "--group", &(tsf.sample)) ;
  getCmdLineOption (&argc, argv, "-g", &(tsf.sample)) ;
  getCmdLineOption (&argc, argv, "-t", &(tsf.title)) ;
  getCmdLineOption (&argc, argv, "--title", &(tsf.title)) ;

  tsf.tagSeparator    = ":" ;
  tsf.sampleSeparator = ":" ;
  tsf.valueSeparator  = "," ;
  tsf.fieldSeparator  = "\t" ;
  getCmdLineOption (&argc, argv, "--TS", &(tsf.tagSeparator)) ;
  getCmdLineOption (&argc, argv, "--SS", &(tsf.sampleSeparator)) ;
  getCmdLineOption (&argc, argv, "--VS", &(tsf.valueSeparator)) ;
  getCmdLineOption (&argc, argv, "--FS", &(tsf.fieldSeparator)) ;
  getCmdLineOption (&argc, argv, "--tagSeparator", &(tsf.tagSeparator)) ;
  getCmdLineOption (&argc, argv, "--sampleSeparator", &(tsf.sampleSeparator)) ;
  getCmdLineOption (&argc, argv, "--valueSeparator", &(tsf.valueSeparator)) ;
  getCmdLineOption (&argc, argv, "--fieldSeparator", &(tsf.fieldSeparator)) ;

  getCmdLineInt (&argc, argv, "--skip1", &tsf.skip1) ;
  getCmdLineInt (&argc, argv, "--skip2", &tsf.skip2) ;
  getCmdLineInt (&argc, argv, "--skip3", &tsf.skip3) ;

  if (argc > 1)
    {
      fprintf (stderr, "Unknown parameter %s, please try tsf --help\n", argv[1]) ;
      exit (1) ;
    }
  fprintf (stderr, "//%s : Start \n"
	   , timeShowNow()
	   ) ;
  
  tsfInit (&tsf) ;
  tsfParse (&tsf) ;
  tsfMerge (&tsf) ;
  tsfExport (&tsf) ;

  if (tsf.tagDict)
    { 
      int mx ;
      messAllocMaxStatus (&mx) ; 
      fprintf (stderr, "// %s done, %d files %d tags %d sample analyzed, max memory %d Mb\n"
	       , timeShowNow()
	       , tsf.nInputFiles
	       , dictMax (tsf.tagDict), dictMax (tsf.sampleDict)
	       , mx) ;
     }
  ac_free (tsf.h) ;

  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderrn and to ensure all pipes are closed*/
  return 0 ;
} /* main */
 
/*************************************************************************************/
/*************************************************************************************/
