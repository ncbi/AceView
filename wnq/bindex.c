/*  file: bindex.c
 *  Author: Richard Durbin et Jean Thierry-Mieg  
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
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
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Optimizer for acedb lookups
 *              maintains a lookup table for the presence of tags
 *              in objects.
 *
 * Since we know the paths to tags, we can be sure some tag is absent
 * if it is not in the model or if its parent in not.
 * The call is 
 *   bIndexFind(key,tag)
 *     return : 0, tag is absent from key
 *              1, not clear
 *              2, present
 *
 * Algorithm.
 *
 *    ttArray[tag] -> tt ( a renumbered tag, with only the useful tags)
 * this may look useless bust costs little, and can be used for subclasses
 *
 *    ctArray[class] -> a[ct]  (array per class, noted a in the code)
 *    a[ct][tt]  -> nBix
 *    bixArray[nBix] -> bix
 *
 *     bix->tag = tag, optimist result = 2
 *        ->tag =  a parent of tag result, optimistic = 1
 *     then we look for bix->->bitset(KEYKEY(key))
 *     and return 0 or optimistic.
 *
 *
 * The idea is to use this system inside bsFindTag in case we are not
 * in a subobj and in the query package BEFORE opening the obj.
 *
 * The bitsets etc will be stored from references in the Class class.
 * Exported functions:
 *    int  bIndexVersion () ;    // to know if index is uptodate
 *    int  bIndexFind(key,tag) ; // is tag in obj, real question
 *    void bIndexObject (KEY key, OBJ obj) ; // reindex a given object
 *    void bIndexInit (int) ;
 *    void bIndexNewModels() ; // when models change
 *    int  bIndexStatus (int *nTablep, int *nKbp) ; // status report
 * HISTORY:
 * Last edited: Nov 12 17:01 2001 (edgrif)
 * Created: August 1 1997, jtm/rd
 * CVS info:   $Id: bindex.c,v 1.12 2015/09/25 17:29:10 mieg Exp $
 *-------------------------------------------------------------------
 */

/* $Id: bindex.c,v 1.12 2015/09/25 17:29:10 mieg Exp $ */

#define DEFINE_OBJ
typedef struct sobj *OBJ ;

#include "acedb.h"
#include "bitset.h"
#include "bs_.h"
#include "bs.h"
#include "a.h"
#include "query.h"
#include "query_.h"
#include "model.h"
#include  "bindex.h"

#define CURRENT_VERSION 4

static KEY _VTag = 0, _Index = 0 ;
static Array ctArray = 0 ;
static int lastTt = 3 ;
static Array ttArray = 0 ;      /* all tags for which some indexing is available */
static Array bixArray = 0 ;  /* array of BIX */
typedef struct { KEY tag ; KEY filter ; int classe ; KEYSET path ; BitSet bb ; BOOL dirty ; int state ; } BIX ;
typedef struct { KEY tag ; int classe ; KEY bbKey ; } BIX1 ;
static int nIndexedNodes = 0, nIndexedObjs = 0 ;
static BOOL indicesMayBeReadFromDisk = FALSE ;
static BOOL indicesAreUptodate = FALSE; 
static BOOL indicesDirty = FALSE; 
  /* set to true if parent session has established the indices 
     or if i ran the indexAllDataBase button 
     */

static BOOL debug = FALSE ;

/*************************************************************************************/

/*
 * Uses index to find out whether a tag is present in an object, if there is
 * some problem with the index, it opens the object and checks for the tag.
 */
BOOL bIndexTag(KEY key, KEY tag)
{ 
  BOOL result = FALSE ;

  if (key && tag && tag >= _Date && class(tag) == 0)
    {
      switch (bIndexFind(key, tag))
	{
	case BINDEX_TAG_ABSENT:
	  result = FALSE ;
	  break ;
	case BINDEX_TAG_UNCLEAR:
	  {
	    OBJ obj = NULL ;
	    
	    if ((obj = bsCreate(key)))
	      {
		result = bsFindTag (obj, tag) ;
		bsDestroy (obj) ;
	      }
	    break ;
	  }
	case BINDEX_TAG_PRESENT:
	  result = TRUE ;
	  break ;
	}
    }

  return result ;
}


/*************************************************************************************/

/*
 * Uses index to find out whether a tag is present in an object, if there is
 * some problem with the index or the tag is definitely there, then opens
 * the object and then tries to get the key from the tag.
 */
BOOL bIndexGetKey(KEY key, KEY tag, KEY *key_out)
{ 
  BOOL result = FALSE ;

  if (key && tag && key_out && tag >= _Date && class(tag) == 0)
    {
      OBJ obj = NULL ;

      if (bIndexFind(key, tag) && (obj = bsCreate(key)))
	{
	  result = bsGetKey(obj, tag, key_out) ;	    /* key_out not changed on failure. */
	  bsDestroy (obj) ;
	}
    }

  return result ;
}

/*******************************************************************/

/*
 * As for bIndexGetKey() but adapted to tag2 system (e.g. SMAP).
 */
BOOL bIndexGetTag2Key(KEY key, KEY tag, KEY *key_out)
{
  BOOL result = FALSE ;

  if (key && tag && key_out && tag >= _Date && class(tag) == 0)
    {
      OBJ obj = NULL ;
      KEY tag2 = 0 ;

      if (bIndexFind(key, tag) && (obj = bsCreate(key)))
	{
	  if (bsGetKeyTags(obj, tag, &tag2))
	    {
	      BSMARK mark = NULL ;

	      do
		{
		  mark = bsMark (obj, mark) ;
		  if (bsGetKey (obj, _bsRight, key_out))    /* key_out not changed on failure. */
		    {
		      result = TRUE ;
		      break ;
		    }
		  bsGoto (obj, mark) ;
		} while (bsGetKeyTags (obj, _bsDown, &tag2)) ;
	    }
	  bsDestroy (obj) ;
	}
    }

  return result ;
}

/*************************************************************************************/

int bIndexVersion (int parent_version)
{
  if (parent_version == CURRENT_VERSION)
    {
      indicesMayBeReadFromDisk = TRUE ;
      return CURRENT_VERSION ;
    }
  else if (parent_version == -1 && /* question at end of session */
	   indicesAreUptodate)
    return CURRENT_VERSION ;

  return 0 ;
}

/*************************************************************************************/
/* status report */
int bIndexStatus (int *nTablep, int *nKbp) 
{
  int i, nn ; 
  BIX* bix ;
   
  if (bixArray)
    {
      *nTablep = arrayMax(bixArray) ;
      
      for (nn = 0, i = 0 ; i < arrayMax(bixArray) ; i++)
	{
	  bix = arrp(bixArray,i,BIX) ; 
	  if (bix->bb)
	    nn += arrayMax(bix->bb) ;
	}
      *nKbp = (nn * sizeof (int)) >> 10 ; /* in kilobytes */
    }

  return bIndexVersion (-1) ;
}

/*************************************************************************************/

BindexFindResult bIndexFind (KEY originalKey, KEY tag)
{
  int result = BINDEX_TAG_ABSENT ;
  int tt, cc = class(originalKey) ;
  Array a ;
  BIX* bix ;
  int nBix, optimisticResult = BINDEX_TAG_PRESENT ;
  KEY key, kk ;
  
  key = lexAliasOf (originalKey) ;	/* follow aliases */

  if (!bsIsTagInClass (class(key), tag))
    return BINDEX_TAG_ABSENT ;

  if (!indicesAreUptodate || tag < _Date)
    return BINDEX_TAG_UNCLEAR ; 

  if (!bixArray)
    return BINDEX_TAG_UNCLEAR ;  /* don't know */

  if (pickList[cc].protected)
    return BINDEX_TAG_UNCLEAR ;                         /* don't mess with protected classes */

  if (pickList[cc].type != 'B')        /* no tag can exist in non-B classes */
    return BINDEX_TAG_ABSENT ;  /* gets called from a query on mixed keysets */

  if (cc >= arrayMax(ctArray))
    return BINDEX_TAG_UNCLEAR ; /* don't know */

  a = arr(ctArray,cc,Array) ;
  if (!arrayExists(a)) 
    return BINDEX_TAG_UNCLEAR ; /* don't know */
  /* this happens when a protected class (not indexed)
     is unprotected in a different kernel version
     it is a weird idea to unprotect but was done for UserSession
     */

  if (tag >= arrayMax(ttArray))
    return BINDEX_TAG_ABSENT ;

  tt = arr(ttArray,tag,int) ;
  if (!tt ||             /* tag absent from every model */
      tt >= arrayMax(a))  /* tag absent from model */
    return BINDEX_TAG_ABSENT ;

  nBix = arr(a,tt,short) ;
  if (!nBix)   /* i.e. array a was never filled */
    return BINDEX_TAG_UNCLEAR ;  /* don't know */
  /* really tag is absent from model and i should return BINDEX_TAG_ABSENT
     but this is not yest tested
     */

  if (nBix >= arrayMax(bixArray))
    return BINDEX_TAG_UNCLEAR ; /* don't know */

  bix = arrp(bixArray,nBix,BIX) ;
  optimisticResult = bix->tag == tag ? BINDEX_TAG_PRESENT : BINDEX_TAG_UNCLEAR ;
  kk = KEYKEY(key) ;

  if (!bix->bb)
    return BINDEX_TAG_ABSENT ;   /* never seen this tag in this class */

  result = bitt(bix->bb,kk) ? optimisticResult : BINDEX_TAG_ABSENT ;

  return result ;
}

/**************************************************************************/

unsigned long int bIndexTagCount (int cc, KEY tag)
{
  int result = 0 ;
  int tt ;
  Array a ;
  BIX* bix ;
  int nBix ;
  
  if (!bsIsTagInClass (cc, tag))
    return -1 ;

  if (!indicesAreUptodate || tag < _Date)
    return -1 ;

  if (!bixArray)
    return -1 ;

  if (pickList[cc].protected)
    return -1 ;

  if (pickList[cc].type != 'B')        /* no tag can exist in non-B classes */
    return -1  ;  /* gets called from a query on mixed keysets */

  if (cc > arrayMax(ctArray))
    return -1  ; /* don't know */

  a = arr(ctArray,cc,Array) ;
  if (!arrayExists(a)) 
    return -1 ; /* don't know */
  /* this happens when a protected class (not indexed)
     is unprotected in a different kernel version
     it is a weird idea to unprotect but was done for UserSession
     */

  if (tag >= arrayMax(ttArray))
    return -1 ;

  tt = arr(ttArray,tag,int) ;
  if (!tt ||             /* tag absent from every model */
      tt >= arrayMax(a))  /* tag absent from model */
    return BINDEX_TAG_ABSENT ;

  nBix = arr(a,tt,short) ;
  if (!nBix)   /* i.e. array a was never filled */
    return -1 ;  /* don't know */
  /* really tag is absent from model and i should return BINDEX_TAG_ABSENT
     but this is not yest tested
     */

  if (nBix >= arrayMax(bixArray))
    return -1 ; /* don't know */

  bix = arrp(bixArray,nBix,BIX) ;
  if (bix->tag != tag)
    return -1 ;

  if (!bix->bb)
    return 0 ;   /* never seen this tag in this class */

  if (bit (bix->bb, 0)) /* model was flagged by old acedb version */
    {
      bix->dirty = TRUE ; 
      bitUnSet (bix->bb, 0) ;
    }
  result = bitSetCount (bix->bb) ;

  return result ;
}

/**************************************************************************/

void bIndexObject (KEY key, OBJ obj)
{
  static BitSet bb = 0 ;
  int n = 0 ;
  BIX *bix ;
  int kk, cc ;
  Array a = 0 ;
  unsigned char mask ;

  /* allow obj = 0, this is used to zero the index of a destroyed object 
     and also to treat named based subclassesd 
     */
  if (obj && !key) key = obj->key ;
  kk = KEYKEY(key), cc = class(key) ;
  
  if (pickList[cc].protected) return ; /* don't mess with protected classes */ 
  if (pickList[cc].type != 'B')        /* no tag can exist in non-B classes */
    messcrash("bIndexObject called on non B class %s", className(key)) ;   
  if (!kk)
    return ; /* do not index the ?class_model object itself  */
  a = pickList[class(key)].conditions ;
  if (arrayExists(a))
    { 
      bb = bitSetReCreate (bb, arrayMax(a)) ;
      mask = queryIsA (obj, key, a, bb) ; /* will create/destroy obj if needed and zero */
      lexSetIsA (key, mask) ;
    }
  else
    bb = bitSetReCreate (bb, 32) ;

  if (!bixArray) return ;
  n = arrayMax(bixArray) ;
  bix = n ? arrp(bixArray,0,BIX) : 0 ;
  bix-- ;
  while (bix++, n--)
    if(bix->classe == cc)
      {
	if (bix->filter)
	  { 
	    if (bitt(bb, bix->filter))
	      {
		if (!bix->bb)
		  bix->bb = bitSetCreate(lexMax(cc),0);
		if (!bitt(bix->bb,kk))
		  { bix->dirty = TRUE ; 
		  indicesDirty = TRUE ;
		  bitSet(bix->bb,kk) ;
		  }  
		nIndexedNodes++ ;
	      }
	    else
	      {
		if (bix->bb && bitt(bix->bb,kk))
		  { 
		    bix->dirty = TRUE ;  
		    indicesDirty = TRUE ;
		    bitUnSet(bix->bb,kk) ;
		  } 
	      }
	  }
	else if (bix->tag)
	  {
	    if (obj && bsFindTag(obj,bix->tag))
	      {
		if (!bix->bb)
		  bix->bb = bitSetCreate(lexMax(cc),0);
		if (!bitt(bix->bb,kk))
		  { bix->dirty = TRUE ; 
		  indicesDirty = TRUE ;
		  bitSet(bix->bb,kk) ;
		  }  
		nIndexedNodes++ ;
	      }
	    else
	      {
		if (bix->bb && bitt(bix->bb,kk))
		  { 
		    bix->dirty = TRUE ;  
		    indicesDirty = TRUE ;
		    bitUnSet(bix->bb,kk) ;
		  } 
	      }
	  }
      }

  nIndexedObjs++ ;
}

/*********/

static void bIndexClass (int classe)
{
  int n = 0 ;
  OBJ obj ;
  KEY key ; 

  if (pickList[classe].type != 'B'  ||
      pickList[classe].protected) 
    return ;

  n = lexMax(classe) ; /* costly, loads the lexique */
  while (n--)
    {
      key = KEYMAKE(classe,n) ;
      obj = bsCreate(key) ;
      bIndexObject(key, obj) ;
      if (obj) bsDestroy(obj) ;
      }
}

/*************************************************************************************/

static void bIndexDatabase (void)
{
  int i = 256 ;
  int nIndexedClasses;

  messStatus ("Indexing (this may take several minutes)") ;
  freeOut ("// Indexing (this may take several minutes)\n") ;
  nIndexedObjs = 0 ;
  nIndexedNodes = 0 ;
  nIndexedClasses = 0;

  while (i--)
    if (pickList[i].type == 'B'  &&
	!pickList[i].protected) 
      {
	bIndexClass(i) ;
	++nIndexedClasses;
      }
  indicesAreUptodate = TRUE ;

  freeOutf ("// Indexed %d tags in %d objects totalling %d nodes \n", 
	    lastTt, nIndexedObjs, nIndexedNodes) ;

  messdump ("Indexed %d tags in %d classes containing %d objects "
	    "totalling %d nodes", 
	    lastTt, nIndexedClasses, nIndexedObjs, nIndexedNodes) ;

  return;
} /* bIndexDatabase */

/**********************************************************************/
/* tools to initialise the indexing */

static void bIndexCreate (int cc, KEY tag, KEY targetTag, int filter)
{
  int tt = 0, targetTt = 0, nBix ;
  Array a = 0 ;
  BIX *bix ;

  if (pickList[cc].protected)
    return ;                    /* don't mess with protected classes */
  if (pickList[cc].type != 'B')
    messcrash("Non B class i bIndexCreate") ;

  if (debug) printf("Create an Index  bix for class %s tag %s targetTag %s filter %d\n",
	   className(KEYMAKE(cc,0)), name(tag), name(targetTag), filter) ; 

  /* is this a new tag ? */
  if (filter)
    {
      lexaddkey (messprintf("_bix_%d_%d",cc, filter), &targetTag, _VSystem) ;
      tag = targetTag ;
    }

  tt = array(ttArray,tag,int) ;
  if (!tt)
    array(ttArray,tag,int) =  tt = ++lastTt ;
  
  targetTt = array(ttArray,targetTag,int) ;
  if (!targetTt)
    array(ttArray,targetTag,int) =  targetTt =  ++lastTt ;


  /* is this a new class */
  a = array(ctArray,cc,Array) ;
  if (!a)
    { 
      a = array(ctArray,cc,Array)  = arrayCreate(2*tt,short) ;
      indicesDirty = TRUE ;
    }

  /* locate target Bix */
  
  nBix = array(a,tt,short) ;
  if (nBix && 
      (bix = arrayp(bixArray,nBix, BIX)) &&
      bix->tag == tag)
    { bix->state |= 1 ; return ;} /* this tag is directly indexed */

  /* locate target Bix */

  nBix = array(a,targetTt,short) ;
  if (!nBix)
    nBix = array(a,targetTt,short) = arrayMax(bixArray) ;

  bix = arrayp(bixArray,nBix, BIX) ;
  if (bix->tag != targetTag)  /* happens when a tag becomes indexable */
    { 
      if (bix->tag)
	{ nBix = array(a,targetTt,short) = arrayMax(bixArray) ;
	  bix = arrayp(bixArray,nBix, BIX) ;
	}
      bix->dirty = TRUE ;
      bix->state |= 1 ;
      if (bix->tag) /* tag was known before, must be reindexed */
	bix->state |= 2 ;
      bix->tag = targetTag ;
      bix->classe = cc ;
      bix->path = bsGetPath (bix->classe, bix->tag) ;
      bitSetDestroy (bix->bb) ;
      array(a,targetTt,short) = nBix ; /* point back to new bix */
    }
  else
    bix->state |= 1 ;
  bix->filter = filter ;
  if (bix->classe != cc)
    messcrash("inconsistency in bIndexCreate, sorry") ;
  array(a,tt,short) = nBix ; /* point to same target */
  if (debug)
    printf("Created bix %d for class %s tag %s (tt=%d) targetTag %s (targetTt=%d)\n",
	   nBix, className(KEYMAKE(cc,0)), name(tag),tt, name(targetTag),targetTt) ; 
}


/*************************************************************************************/

static void bIndexInitStatics (void) 
{
  KEY dummy ;

  lexaddkey ("Tag", &dummy, _VMainClasses) ; _VTag = KEYKEY(dummy) ;
  lexaddkey("Index",&_Index, _VSystem) ;
}

/*************************************************************************************/
/* recursive exploration of the models
 * a masterTag may occur anyway down the recursion on user request
 */

static void bIndexRegisterSecondaryTag(OBJ obj, int cc, BS bs, KEY masterTag)
{
  if(!bs) return ;
  if (bs->key == _UNIQUE)
    { bIndexRegisterSecondaryTag(obj, cc, bs->right, masterTag) ; return ; }
  if (bs->key == _XREF && bs->right)
    { bIndexRegisterSecondaryTag(obj, cc, bs->right->right, masterTag) ; return ;}
  if (class(bs->key) || 
      bs->key < _Date)
    return ;    /* break on values */
  if (bs->down)
    bIndexRegisterSecondaryTag(obj, cc, bs->down, masterTag) ;
  /* no condition, just index every tag */
  masterTag = bs->key ;
#ifdef ACEDB4
  {  /* done to convince pre code before march.2000 to index everything */ 
    KEY tag = 0 ;
    lexaddkey(name(masterTag),&tag,_VTag) ;
    if (tag) bsAddKey(obj,_Index,tag) ; 
  }
#endif
  if (bs->right)
    bIndexRegisterSecondaryTag(obj, cc, bs->right, masterTag) ;
  if (debug) printf("bIndexRegisterSecTag class %s, tag %s masterTag %s\n",
	 className(KEYMAKE(cc,0)), name(bs->key), name(masterTag)) ; 
  if (masterTag && bs->key) 
    bIndexCreate(cc,bs->key,masterTag,0) ;
}

/*************************************************************************************/

static void bIndexRegisterTags (void)
{
  KEY cc ; unsigned char mask ;
  BS bs ;
  KEY classKey = 0, model ;
  OBJ obj = 0, Model = 0 ;

  while (lexNext(_VClass,&classKey))
    if ((obj = bsUpdate (classKey)))
      {
	cc = classKey ; mask = 0 ;
	pickIsA (&cc, &mask) ;
	if (cc && !mask &&
	    pickList[cc].type == 'B' &&
	    !pickList[cc].protected)
  	    if ((model =  pickList[cc].model) &&
		((Model = bsCreate(model))))
	      {
		bs = Model->root->right ;
		bIndexRegisterSecondaryTag(obj, cc, bs, 0) ;
		bsDestroy (Model) ;
	      }
	bsSave (obj) ;
      }
}

/*************************************************************************************/

static void bIndexRegisterFilters (void)
{
  int i, cc = 256 ;
  Array a = 0 ;

  while (cc--)
    if (pickList[cc].type == 'B' && !pickList[cc].protected)
      {
	a = pickList[cc].conditions ;
	if (arrayExists(a) && arrayMax(a) > 8)
	  for (i = 8 ; i < arrayMax(a) ; i++)
	    bIndexCreate(cc,0, 0, i) ;
      }
}

/*************************************************************************************/

static BOOL bIdenticalPath (KEYSET p1, KEYSET p2)
{
  int i, j ;
  KEY *kp1, *kp2 ;
  
  if (!p1 || !p2)
    return FALSE ;

  i = keySetMax(p1) ; j = keySetMax(p2) ; 
  /* do not rely on max of a path, look for a terminating zero */
  kp1 = arrp(p1, 0, KEY) ; 
  kp2 = arrp(p2, 0, KEY) ; 
  while (*kp1 && i-- && j--) /* for some reason arrayMax(path) is meaningless */
    if (*kp1++ != *kp2++) 
      return FALSE ;
  if ((i >= 0 && *kp1) || (j >= 0 && *kp2))
    return FALSE ;
  return TRUE ;
}

/*************************************************************************************/
/* 
   When new models are read, one should adapt the indexation
   system, but the objects already indexed should be ok 

   there is a problem for those tags that are indexed and are moved
   because bIndex will reply 2 (present) but bsFindTag will not find them
   so it is necessary to check that the whole path did not change
*/

void bIndexNewModels(void)
{
  BIX *bix ;
  int i, nBix, nIndexedClasses ;
  KEYSET path = 0 ;
  Array a = 0 ;

  if (!indicesAreUptodate)
    {
      bIndexInit(BINDEX_FORCE) ;
      return ;
    }

  bIndexInitStatics() ;
  messStatus ("Indexing") ;
  freeOut ("// Indexing \n");
  nBix = arrayMax(bixArray) ;
  bix = arrayp(bixArray,0, BIX) - 1 ;
  while (bix++, nBix--) 
    bix->state = 0 ;

  bIndexRegisterTags () ; 
  nBix = arrayMax(bixArray) ;
  bix = arrayp(bixArray,0, BIX) - 1 ;
  a = arrayCreate (256,int) ;
  while (bix++, nBix--) 
    switch (bix->state)
      {
      case 0: /* this bix is useless, destroy it */    
	if (0 && bix->tag)
	  messdump ("Detroying the index for class %s targetTag %s\n",
		    className(KEYMAKE(bix->classe,0)), name(bix->tag)) ; 
	if (bix->bb)
	  bitSetDestroy (bix->bb) ;
	keySetDestroy (bix->path) ;
	bix->tag = 0 ;
	bix->classe = 0 ;
	bix->dirty = TRUE ; indicesDirty = TRUE ;
	break;
      case 1:   /* this bix is reused, check it */
	path = bix->path ;
	bix->path = bsGetPath (bix->classe, bix->tag) ;
	if (bix->bb && !bIdenticalPath (path, bix->path))
	  {
	    if (path) 
	      messdump ("In model ?%s, tag %s has moved, reindexing\n",
			className(KEYMAKE(bix->classe,0)), name(bix->tag)) ;
	    array(a,bix->classe,int) = 1 ; 	
	    bix->dirty = TRUE ; indicesDirty = TRUE ;
	  }
	keySetDestroy (path) ;
	break ; 
      case 3:	/* htis tag became indexable */
	keySetDestroy (bix->path) ;
	bix->path = bsGetPath (bix->classe, bix->tag) ;
	array(a,bix->classe,int) = 1 ; 
	bix->dirty = TRUE ;  indicesDirty = TRUE ;
	break ;
      }
  i = arrayMax(a) ;  
  nIndexedObjs = 0 ;
  nIndexedClasses = 0 ;
  nIndexedNodes = 0 ;

  while(i--)
    if (arr(a,i,int))
      { bIndexClass(i); nIndexedClasses++ ; }
  arrayDestroy (a) ;
  freeOutf ("// Reindexed %d tags in %d classes containing %d objects totalling %d nodes \n", 
	    lastTt, nIndexedClasses, nIndexedObjs, nIndexedNodes) ;

  messdump ("Indexed %d tags in %d classes containing %d objects "
	    "totalling %d nodes", 
	    lastTt, nIndexedClasses, nIndexedObjs, nIndexedNodes) ;

  return;
} /* bIndexNewModels */

/*********************************************************************/

BOOL  bIndexSave(void)
{
  int i, j ;
  KEY key ;
  BIX *bix ; BIX1* bix1 ;
  Array cta = 0, bixa = 0 ;

  if (!indicesDirty || !indicesAreUptodate)
    return FALSE ;
  
  /* ghost copy ttArray */
  lexaddkey("__ttArray",&key,_VCalcul) ;
  arrayStore (key,ttArray,"k") ;

  /* ghost copy ctArray */
  i = arrayMax(ctArray) ;
  cta = arrayCreate (i, KEY) ;
  while (i--)
    { 
      Array a1, a = array (ctArray, i, Array) ;
      if (!a)
	continue ;
      lexaddkey(messprintf("__cta_%d",i),&key,_VCalcul) ;
      j = arrayMax(a) ;
      a1 = arrayCreate (j, int) ;
      while(j--)
	array(a1,j,int) = array(a,j,short) ;
      arrayStore (key,a1,"i") ;
      arrayDestroy (a1) ;
      array(cta,i,KEY) = key ;
    }
  lexaddkey("__ctArray",&key,_VCalcul) ;
  arrayStore (key,cta,"k") ;
  arrayDestroy (cta) ;


  /* ghost copy bixArray */
  i = arrayMax(bixArray) ;
  bixa = arrayCreate (i, BIX1) ;
  while (i--)
    { 
      BitSet a ;
      bix = arrayp(bixArray, i, BIX) ;
      bix1 = arrayp(bixa, i, BIX1) ;
      bix1->tag = bix->tag ;
      bix1->classe = bix->classe ;
      a = bix->bb ;
      key = 0 ;
      if (!a)
	{  /* look for old stuff and destroy it */
	  lexword2key(messprintf("__bixa_%d",i),&key,_VCalcul) ;
	  if (key && iskey(key) == 2) /* happens from bIndexNewModels */
	    {
	      arrayKill(key) ;  
	      bix->dirty = 0 ;
	    }
	  continue ;
	}
      lexaddkey(messprintf("__bixa_%d",i),&key,_VCalcul) ;
      if (bix->dirty)
	bigArrayStore (key,a,"i") ;
      bix1->bbKey = key ;
      bix->dirty = 0 ;
    }
  lexaddkey("__bixArray",&key,_VCalcul) ;
  arrayStore (key,bixa,"kik") ;
  arrayDestroy (bixa) ;

  indicesDirty = FALSE ;
  return TRUE ;
}

/*************************************************************************************/

static BOOL bIndexInitFromDisk(void)
{
  int i, j, *ip ;
  KEY key ;
  BIX *bix = 0 ; BIX1* bix1 ;
  Array tta = 0, cta = 0, bixa = 0, a = 0 , a1 = 0 ;

  if (!indicesMayBeReadFromDisk)
    return FALSE ;

  bIndexInitStatics () ;

  if (!lexword2key("__ttArray",&key,_VCalcul) ||
      !(tta = arrayGet(key,KEY,"k")))
    goto abort ;
  
  if (!lexword2key("__ctArray",&key,_VCalcul) ||
      !(cta = arrayGet(key,KEY,"k")))
    goto abort ;
  
  if (!lexword2key("__bixArray",&key,_VCalcul) ||
      !(bixa = arrayGet(key,BIX1,"kik")))
    goto abort ;
  
  /* reestablish ttArray */
  ttArray = tta ; tta = 0 ;
  lastTt = 3 ;
  ip = arrp(ttArray,0,int) - 1 ;
  i = arrayMax(ttArray) ;
  while(ip++, i--)
    if (*ip > lastTt) lastTt = *ip ;

  /* reestablish ctArray */
  i = arrayMax(cta) ;
  ctArray = arrayCreate (i, Array) ;
  while (i--)
    { 
      key = array(cta,i,KEY) ;
      if (key)
	{
	  a1 = arrayGet(key, int, "i") ;
	  if (a1)
	    { 
	      j = arrayMax(a1) ;
	      a = arrayCreate (j, short) ;
	      while(j--)
		array(a,j,short) = array(a1,j,int) ;
	      arrayDestroy (a1) ;
	      array(ctArray,i,Array) = a ; a = 0 ;
	    }
	  else
	    goto abort ;
	}
    }
  arrayDestroy (cta) ;

  /* reestablish bixArray */
  i = arrayMax(bixa) ;
  bixArray = arrayCreate (i, BIX) ;
  while (i--)
    { 
      bix = arrayp(bixArray, i, BIX) ;
      bix1 = arrayp(bixa, i, BIX1) ;
      bix->tag = bix1->tag ;
      bix->classe = bix1->classe ;
      bix->path = bsGetPath (bix->classe, bix->tag) ;
      key = bix1->bbKey ;
      if (key)
	{
	  BigArray a = bigArrayGet(key, unsigned int, "i") ;
	  if (a)
	    { bix->bb = a ; a = 0 ; }
	  else
	    goto abort ;
	}
      bix->dirty = 0 ;
    }
  arrayDestroy (bixa) ;

  indicesAreUptodate = TRUE ; 
  return TRUE ;

abort:
  arrayDestroy (tta) ;
  arrayDestroy (bixa) ;
  arrayDestroy (a) ;
  arrayDestroy (ttArray) ;
  if (ctArray)
    {
      i = arrayMax(ctArray) ;
      while (i--)
	arrayDestroy (array(ctArray, i, Array)) ;
      arrayDestroy (ctArray) ;
    }
  if (bixArray)
    {
      i = arrayMax(bixArray) ;
      while (i--)
	{ 
	  bitSetDestroy(array(bixArray, i, BIX).bb) ;
	  keySetDestroy (bix->path) ;
	}
      arrayDestroy (bixArray) ;
    }
  
  return FALSE ;
}

/*************************************************************************************/

static void bIndexDestroy (void)
{
  int i ;
  
  if (bixArray)
    for (i = 0 ; i < arrayMax(bixArray) ; i++)
      { 
	arrayDestroy (arrp(bixArray, i, BIX)->path) ;
	bitSetDestroy (arrp(bixArray, i, BIX)->bb) ;
      }
  arrayDestroy (bixArray) ;
  arrayDestroy (ttArray) ;
  if (ctArray)
    for (i = 0 ; i < arrayMax(ctArray) ; i++)
      arrayDestroy (array(ctArray, i, Array)) ;
  arrayDestroy (ctArray) ;
}

/*************************************************************************************/

/* bIndexInit is called in several different ways, see bindex.h for a        */
/* description of the various values of  askbIndexInit.                      */
/*                                                                           */
/* I think the flow control in this routine is not good, sometime it could   */
/* do with tidying up to make it clearer.                                    */
/*                                                                           */
void bIndexInit(BindexInitType askbIndexInit)
{
  int minSession = 2 ; /* 2, to automatically index a new database
                        * 0, not to
			*/

  if (askbIndexInit == BINDEX_AFTER_READMODELS)
    { 
      indicesAreUptodate = FALSE ; /* will force total reindexing  at the end of readmodels */
      bIndexDestroy () ;
      readModels () ;
      return ;
    }
  if (askbIndexInit == BINDEX_FORCE)
    goto force ;

  if (getenv("ACEDB_NO_INDEX"))  /* if true, indices wont be used, even if available */
    {
      indicesAreUptodate = FALSE ; 
      return ;
    }

  if (thisSession.session > -1 && (indicesAreUptodate  || bIndexInitFromDisk()))
    return ;

  if (getenv("ACEDB_INDEX")) 
    minSession = 2 ; /* if true, a new database will be automatically indexed */

  if (thisSession.session > minSession &&
      !messQuery("Do you want to index this database, this may be slow, "
		 "but will accelerate future access"))
    return ;

force:

  if (!sessionGainWriteAccess())
    {
      if (askbIndexInit == BINDEX_FORCE)
	messout ("Sorry, you cannot gain write access, "
		 "I cannot index this database") ;
      return ;
    }

  bIndexInitStatics () ;

  indicesDirty = TRUE ;

  bIndexDestroy () ;
  bixArray = arrayCreate (512,BIX) ;
  array(bixArray,0,BIX).tag = 0 ; /* avoid bix zero */
  ctArray = arrayCreate(256,Array) ;
  ttArray = arrayCreate(512,int) ;

  lastTt = 3 ;

  bIndexRegisterTags () ;
  bIndexRegisterFilters () ;
  bIndexDatabase() ;

  return ;
}

/*************************************************************************************/
/* clean up all memory created in bIndexInit */
void bIndexShutDown (void)
{
  /*  ctArray ttArray  bixArray */ 
  int n ;
  BIX *bix ;
  Array *ap ;
 
  n = bixArray ? arrayMax(bixArray) : 0 ;
  bix = n ? arrp(bixArray,0,BIX) - 1 : 0 ;
  while (bix++, n--)
    {
      bitSetDestroy (bix->bb) ;
      arrayDestroy (bix->path) ;
    }
  arrayDestroy (bixArray) ;
  
  n = ctArray ? arrayMax (ctArray) : 0 ;
  ap = n ? arrp(ctArray,0,Array) - 1 : 0 ;
  while (ap++, n--)
    arrayDestroy (*ap) ;
  arrayDestroy (ctArray) ;
  arrayDestroy (ttArray) ;
}

/*************************************************************************************/
/*************************************************************************************/
