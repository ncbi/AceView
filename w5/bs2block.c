/*  File: bs2block.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: bs2block.c,v 1.8 2014/05/31 01:30:52 mieg Exp $ 
 * Description:
 * Exported functions:
 * HISTORY:
    mieg: jan 97 
          change double recursion to down loop +rigth recursion and
          tested the code
 * Last edited: Dec  4 12:58 1998 (fw)
 * Created: Wed Apr 12 18:28:47 1995 (rd)
 *-------------------------------------------------------------------
 */

#include "acedb.h"

#include "../wh/byteswap.h"
#include "block.h"
#include "bs_.h"
#include "systags.h"
#include "adisk.h"

/*
   strategy: 
     make an Array bData of BSdata, with node keys and data,
        and a parallel Array cData of char, with text and tree structure
     when storing, first find the sizes of these arrays, so no tests
        needed when converting
     use new arrayStore/Get2 functions to store/restore these
*/

/**********************************************/

/* static Array cData ; */
static char *cp ;
static BSdata *bp = 0 ;
static int nc, nb ;
static AC_HANDLE myHandle ;
static int MM = 2 ; /* 1: double recursion, 2: loop/recursion */
static BS nodeUnPack (void)
{ 
  BS bs0 = BSalloc(), bs = bs0 ;
  char c ;

plusbas:
  bs->key = maybeSwapKEY(bp->key) ; bp++ ;
  if (bs->key <= _LastC)
    { char *s ; int i = 0 ;
      bs->bt = BTalloc() ;
      while(*cp--) i++ ; cp += i + 1 ;
      s = bs->bt->cp = halloc (i + 1, myHandle) ; 
      while ((*s++ = *cp--)) ;
    }
  else if (bs->key <= _LastN)
    {bs->n.key =  maybeSwapKEY(bp->key) ; bp++ ;}
  else
    bs->key = lexAliasOf(bs->key) ;

  c = *cp-- ;
  if (c & 0x10)
    { 
      bs->right = nodeUnPack () ;
      bs->right->up = bs ;
    }
  if (c & 0x01)
    {
      if (MM == 1)
	{
	  bs->down = nodeUnPack () ;
	  bs->down->up = bs ;
	}
      else
	{
	  bs->down = BSalloc() ;
	  bs->down->up = bs ;
	  bs = bs->down ;
	  goto plusbas ;
	}
      
    }

  return bs0 ;
}

static void nodePack (BS bs)
{
plusbas:
  (bp++)->key =  maybeSwapKEY(bs->key) ;
  if (bs->key <= _LastC)
    {
      if (bs->bt && bs->bt->cp)
	{ char *s = bs->bt->cp ;
	while ((*cp-- = *s++)) ;
	}
      else
	*cp-- = 0 ;
    }
  else if (bs->key <= _LastN)
    (bp++)->key = maybeSwapKEY(bs->n.key) ;

  *cp = 0 ;
  if (bs->right)
    *cp |= 0x10 ;
  if (bs->down)
    *cp |= 0x01 ;
  cp-- ;

  if (bs->right)
    nodePack (bs->right) ;
  if (bs->down)
    {
      if (MM == 1)
	nodePack (bs->down) ;
      else
	{
	  bs = bs->down ;
	  goto plusbas ;
	}
    }
}

static void nodeSize (BS bs)
{ 
 plusbas:
  nb++ ; nc++ ;
  if (bs->key <= _LastC)
    { bs->size = bs->bt && bs->bt->cp ? strlen(bs->bt->cp) + 1 : 1 ;
      nc += bs->size ;
    }
  else if (bs->key <= _LastN)
    nb++ ;

  if (bs->right)
    nodeSize (bs->right) ;
  if (bs->down)
    {
      if (MM==1)
	nodeSize (bs->down) ;
      else
	{ bs = bs->down ; goto plusbas ; }
    }
}

/**********************************************/
/************ external functions **************/
#define BSDATASIZE sizeof(BSdata)
#include "a.h"
#include "pick.h"

BS newBsTreeGet (KEY key, AC_HANDLE handle) 
{ 
  BS bs ;
  int nt ;
  Array bData = 0 ;		/* of BSdata */

  key = lexAliasOf(key) ;
  if (!key || (pickType(key) != 'B'))
    messcrash ("bad call to newBsTreeGet : %s:%s", className(key), name(key)) ;
  
  if (!aDiskArrayGet (key, 0, &bData, 0) ||
      !bData || !(nt = arrayMax(bData))) 
    {
      bs=BSalloc() ;
      bs->key=key;
    }
  else
    {
      bp = arrp (bData, 0, BSdata) ;
      myHandle = handle ;
      cp = (char*)(arrp (bData, nt - 1, BSdata)) + BSDATASIZE - 1 ;
      bs = nodeUnPack () ;
    }

  arrayDestroy (bData) ;
  return bs ;
}

static Array  waitingSet = 0 ;
typedef struct wwStruct { KEY key; Array a ; DISK d ; } WW ;

void newBsTreeStore (BS bs, BOOL dowait) /*, Array* bdp)  , Array *cdp) */
{ 
  int nt ;
  KEY key = bs->key ;
  extern void testAdiskArrayStore (KEY key, Array a) ;
  Array bData = 0 ;		/* of BSdata */

  key = lexAliasOf(key) ;
  if (!key || (pickType(key) != 'B') || !bs)
    messcrash ("bad call to newBsTreeGet : %s:%s", className(key), name(key)) ;

  nb = nc = 0 ;
  nodeSize (bs) ;		/* get BSdata and char array sizes */

  if (nb & 0x01)		/* 8 byte boundary good for disk? */
    ++nb ;
  nt = nb + nc/BSDATASIZE + 1 ;
  bData = arrayCreate (nt, BSdata) ; 
  array (bData, nt-1, BSdata).key = 0 ;

  bp = arrp (bData, 0, BSdata) ;
  cp = (char*)(arrp (bData, nt - 1, BSdata)) + BSDATASIZE - 1 ;
  nodePack (bs) ;

  if (dowait)
    {
      WW *ww ;

      if (!waitingSet) 
	  waitingSet = arrayCreate (128, WW) ;

      nt = arrayMax (waitingSet) ;
      ww = arrayp(waitingSet, nt, WW) ;
      ww->key = key ; ww->a = bData ; 
      ww->d = aDiskAssign (key, arrayMax(bData), sizeof(BSdata), 0, 0) ;  
    }
  else
    {
      aDiskArrayStore (key, 0, bData, 0) ;
      arrayDestroy (bData) ;
    }
}

static int wwOrder (const void *a, const void *b)
{ return ((const WW*)a)->d - ((const WW*)b)->d ; } 

void newBsTreeFlush (void)
{
  WW *ww ; int i ;

  if (waitingSet && arrayMax(waitingSet))
    {
      arraySort (waitingSet, wwOrder) ;
      i = arrayMax (waitingSet) ;
      ww = arrayp(waitingSet, 0, WW) - 1;
      while (ww++, i--)
	{
	  aDiskArrayStore (ww->key, 0, ww->a, 0) ;
	  arrayDestroy (ww->a) ;
	}
    }
  arrayDestroy (waitingSet) ;
}
/**************** end of file *****************/
