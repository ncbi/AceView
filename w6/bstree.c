/*  File: bstree.c
 *  Author: R.Durbin & J.Thierry-Mieg.
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
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
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Transforms back and forth the B blocks into BS trees 
 * Exported functions: 
 * HISTORY:
 * Last edited: Nov  9 13:53 1999 (fw)
 * * Aug 26 17:01 1999 (fw): added this header
 * Created: Thu Aug 26 16:59:46 1999 (fw)
 * CVS info:   $Id: bstree.c,v 1.11 2014/11/10 18:05:09 mieg Exp $
 *-------------------------------------------------------------------
 */

/* @(#)bstree.c	1.13 4/15/97 */


#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct bblock *BBLOCK ;
typedef BBLOCK BP ;
#define DEFINE_OBJ
typedef struct sobj *OBJ ;

#include "acedb.h"

#include "bs_.h"
#include "b_.h"
#include "block.h"
#include "byteswap.h"

#include <ctype.h>

#if defined(MACINTOSH)
#include "memory.h"
#endif

static BS     bsTreeRead     (BP bp,NODE nn,BS b0,BOOL doUnbreakComment, AC_HANDLE handle);
static void   bsTreeReadCont  (BS *bsp, AC_HANDLE handle);
static BP     bsTreeFindBlock  (BS  bs,KEY *kp, AC_HANDLE handle) ;
static void   bsTreeStore2   (BP bp, KEY kk, BS  bs);
static NODE   bsTreeSt2      (BS  bs,BP bp);
static void   bsTreeStoreContinuation  (KEY k,BP *bpp, BS  bs) ;
static int    bsTreeCellSize     (BS  bs, AC_HANDLE handle);
static int    bsTreeBreakComment(BS bs, AC_HANDLE handle);
static void   bsTreeUnbreakComment(BS bs, BS bs1, AC_HANDLE handle);
static void   bsPackTimeStamps (BS bs, KEY currStamp) ;
static void   bsUnpackTimeStamps (BS bs, KEY currStamp) ;

/* extern int BSTEST; */
int bsTreeTEST = FALSE ;
static  char BSbuffer[BLOC_SIZE];  /*to allow for the \t*/

/**************************************************************/
/**************************************************************/
                 /* cannot fail */
BS bsTreeGet(KEY key, AC_HANDLE handle)
{
 NODE nn;
 BP bp;
 BS bs;

 chrono("bsTreeGet") ;

 if(!iskeyold(key))
   { bs=BSalloc() ;
     bs->key=key;
   }
 else
   {
     Bget(&bp,&nn,key) ;  /* *bp,*nn becomes the origin of the branch*/
     /*  unalias on the fly does not work because lexAliasOf
	 may trigger a diskRead in the middle of the current disk read
	 creating havoc
	 it would work if we had read all lexiques in advance, but 
	 this would terribly slow down the start up
	 may be doing the lexAliasOf always in bsAdd() is not so costly
	 also it will remove the __lexHasNewAlias horrible static
     */
     bs = bsTreeRead (bp,nn,(BS)0,TRUE, handle) ;

     bsUnpackTimeStamps (bs->right, 0) ;

     blockunpinn(key);
   }
 
 chronoReturn(); 
 return bs ;
}

/********************************************************************/
                 /*store => free*/
        /* because the splitting of big objects destroys the tree */

/* XQX
* This function gets called on objects that need to be saved.  It
* really just splits up large objects into smaller trees that will
* each fit in a database block, then uses bsTreeStore2() to do the
* actual storing.
*/
void bsTreeStore (BS  bs, AC_HANDLE handle)
{ static Array warn = 0 ;

  KEY key=bs->key,kk ;
  int
    nc=0,         /*number of continuation blocks*/
    ncells, nl=0, nr=0,  nmax ;
  BS bs2, bsl,bsr, bsu, bsn, bsc;
  BP bp , bpc;

  chrono("bsTreeStore") ;

  /* XQX
  * Get the timestamps (which are really _VUserSession object keys)
  * off each BS node and insert them as nodes of their own in the 
  * object tree.
  */
  bsPackTimeStamps (bs->right, 0) ;

           /* set bs -> size, no more recursive sizing needed later */
  ncells=bsTreeCellSize(bs, handle);
           /* set bp to a proper top block */
  if (ncells > 100000)
    { BS bb = bs ;
      int n ;

      while (bb->up) bb = bb->up ;
      n = class(bb->key) ;
      if (!warn) warn = arrayCreate(256, int) ;
      array(warn,n,int)++ ;
      if (array(warn,n,int) == 1)
	messdump ("Class %s, object %s has %d > 100000 cells.\n"
		  "This is just a warning, acedb has no hard limits on the mumber of cells per object,\n"
		  "but the performance degrades on very large objects, and it possible you should \n"
		  "reconsider your models design.  In particular, you may be overusing XREF or ?Text.\n"
		  "This message is issued once per offending class and per code run\n",
		  className(bb->key), name(bb->key), ncells) ;
    }

  /* XQX
  * bpc, bp contain an extremely large data structure.
  */
  bpc = bp = bsTreeFindBlock(bs,&kk, handle) ;
  nmax = BNODEMAX-2 ;

  /* XQX
  * BNODEMAX is how many cells fit into a database block.  I guess
  * this loop breaks the tree into pieces somehow.  If you have a
  * small object, the loop body never executes.
  */
  while ((ncells=bs->size) >= BNODEMAX)
    {
      bs2=bs;
      nl = nr = 0 ;
      while (ncells > nmax)
        { bsl = bs2->down ; bsr = bs2->right ;
          nl = bsl ? bsl->size : 0 ;
          nr = bsr ? bsr->size : 0 ;
	  if (nl > nr)
	    { ncells = nl ; bs2 = bsl ; }
	  else
	    { ncells = nr ; bs2 = bsr ; }
        }

      if (!bs2 || ncells < 12)
        messcrash (" Error in bsTreeStore , check bsTreeCellSize() :\n%s",
		   name(bs->key));

      /* bs2 is unhooked from bsu and hooked to the new bsn */

      nc++ ;
      bsn = BSalloc() ;
      bsn->key = key ;
      bsn->right = bs2 ;
      bsu = bs2->up ; bs2->up = bsn ;
      bsTreeStoreContinuation (key, &bpc, bsn) ;
      		   /* bsc will now replace the bs2 tree under bsu */
      bsc = BSalloc () ;
      bsc->up = bsu ;
      if (nl > nr)
        bsu->down = bsc;
      else
        bsu->right = bsc;
      bsc->key = _continuationKey  ;
      bsc->n.i = nc ;
      bsc->size = 2;
      ncells -= 2 ; /* The tree is pruned change the size above bsc*/
      while(bsu)
        {
          bsu->size -= ncells;
          bsu = bsu->up;
        }
    }

				/* save top of tree in original block */
  bsTreeStore2 (bp,kk,bs) ;
				/* discard anything further down */
  blockSetEnd (bpc) ;
  blockrewrite (kk) ;

  chronoReturn () ;
}

/******************************************************************/
                 /* recursive pice,  calls bsTreeBreakComment */
static BS bsxxx ;
static int sizexxx ;
static AC_HANDLE hxxx ;
static void  bsTreeCellSize2   ()

{
  if(!bsxxx)
    return ;

  bsxxx->size = 0;

  if(bsxxx->down)
    {
     bsxxx = bsxxx->down ;
     bsTreeCellSize2() ;
     bsxxx = bsxxx->up ;
     bsxxx->size += sizexxx ;
   }

  if(bsxxx->right)
    {
     bsxxx = bsxxx->right ;
     bsTreeCellSize2() ;
     bsxxx = bsxxx->up ;
     bsxxx->size+= sizexxx ;
   }

  if  (bsxxx->key<=_LastC)
    {
      if(bsxxx->bt && bsxxx->bt->cp)
	{
	  register int i = strlen(bsxxx->bt->cp) ;
	  if(i > NODEX * BNODEMAX /2)
	    {
	      sizexxx = bsTreeBreakComment(bsxxx, hxxx);
	      return ;
	    }
	  else
	    { if (i == 0)
		bsxxx->size += 2; /* SRK - empty string takes two also */
	      else
		bsxxx->size+=1+(i+NODEX-1)/NODEX;
	    }
	}
      else   /* a null comment occupies 2 keys, sorry */
	  bsxxx->size+=2;
    }

  else if(bsxxx->key<=_LastN)
    bsxxx->size += 2 ;

  else
    bsxxx->size += 1 ;  /* generic case of a KEY */

  sizexxx = bsxxx->size ;
  return ;
}

/******************************************************************/
                 /* calls bsTreeBreakComment */
static int    bsTreeCellSize    (BS bs, AC_HANDLE handle)
{
  if(!bs)return(0);

  chrono("bsTreeCellSize") ;
  bsxxx = bs ;
  hxxx = handle ;
  bsTreeCellSize2() ;

  chronoReturn(); return(bs->size);
}
/**************************************************************/
static int  bsTreeBreakComment(BS bs, AC_HANDLE handle)
{
  char * old = bs->bt->cp, * line ;
  BS bs1, bsu;

 int n ;

 chrono("bsTreeBreakComments") ;

 uLinesText(old, (BNODEMAX * NODEX)/ 2);
 line = uPopLine(old) ;

 while(TRUE)
   {
     if (!bs->bt)
       bs->bt = BTalloc() ;
     bs->bt->cp = strnew (line, handle) ;
     n = 2+(strlen(bs->bt->cp)+NODEX-1)/NODEX;  /* michel  was 1+.. */
     if (bs->down)
       n += bs->down->size ;
     if (bs->right)
       n += bs->right->size ;
     bs->size = n ;

     line = uPopLine(old);
     if (!line)
       break;
     bs1 = BSalloc();

     bs1 -> key = bs->key ;
     bs -> key = _NextC;

     bs1->up = bs->up;
     bs->up = bs1;
     bs1->down = bs;       /* move the right hook to the top */
     bs1->right = bs->right ;
     if(bs1->right)
       {
	 bs1->right->up = bs1 ;
	 bs->size -=  bs->right->size ;
       }
     bs->right = 0 ;

     bsu = bs1->up;
     if (bsu)
       {
	 if (bsu->down == bs) bsu->down = bs1;
	 else
	   if (bsu->right == bs) bsu->right = bs1;
	   else
	     messcrash("link missing in bsTreeBreakComment");
       }
     bs = bs1;
   }
 messfree(old);
 bsxxx = bs ;
 chronoReturn(); return n;
}

/**************************************************************/
static void bsTreeUnbreakComment (BS bs, BS bs1, AC_HANDLE handle)
{

 char * u = bs->bt ? bs->bt->cp : 0 ,
      * v = bs1->bt ? bs1->bt->cp : 0 ,
      * cq ;
 int n = ( u ? strlen(u) : 0 )  + ( v ? strlen(v) : 0 ) ;

 chrono("bsTreeUnbreakComments") ;

 cq = halloc(n + 1, handle);
 *cq = 0;
 if(u)
   {
     strcat(cq,u);
     bs->bt->cp = 0 ;
   }
 if(v)
   {
     strcat(cq,v);
     bs1->bt->cp = 0 ;
   }


 if (!bs->bt)
   bs->bt = BTalloc() ;
 bs->bt->cp = cq;
                 /* The right hook is already on bs */
 if ((bs->down = bs1->down)) /* i do mean =, not == */
   bs->down->up = bs ;
 bs1->down = bs1->right =  bs1->up = 0;
 BSfree(bs1);
 chronoReturn() ;

}

/**************************************************************/
      /*note that if b0=0 we do not read bs->down */

static BS  bsTreeRead (BP bp,NODE nn,BS b0, BOOL doUnbreakComment, AC_HANDLE handle)
{
  NODEP r=bp->n+nn ;
  BS bs=BSalloc() ;
  BS bs1 ;

  chrono("bsTreeRead") ;

  bs->up=b0;
  bs->key=maybeSwapKEY(r->ck.key) ;

  /* down hook */
  bs->down = (b0 && r->d1) ?
    bsTreeRead (bp,r->d1,bs,doUnbreakComment, handle) : NULL;

  if (bs->key<=_LastC)
    {
      if(r->d2)
	{
	  Bnode2str(BSbuffer,bp,&r);
	  bs->bt = BTalloc();
	  bs->bt->cp = strnew(BSbuffer, handle) ;

	  /* special right hook for compatibility with earlier version*/
	  bs->right = r->d1  ?
	    bsTreeRead (bp,r->d1,bs,doUnbreakComment, handle) : NULL;
	}
      if(    doUnbreakComment
	     && (bs1 = bs->down)
	     && (bs1->key == _NextC) )
        bsTreeUnbreakComment(bs, bs1, handle) ;
      chronoReturn(); return(bs);
      /* this call may have changed r1 */
    }

  else  if(bs->key<=_LastN)
    {           /* In a number a BS cell corrresponds to 2 nodes */
      nn=r->d2;
      r = bp->n + nn ; /* Thus the right hook
			  will start from the second cell */
      bs->n.i  = maybeSwapInt(r->ck.x.i) ;
      if(bs->key == _continuationKey)
	{
	  bsTreeReadCont(&bs, handle) ; 
	  chronoReturn(); return(bs);
	}
    }
                                                       /* right hook */
    bs->right = r->d2  ? bsTreeRead (bp,r->d2,bs,doUnbreakComment, handle) : NULL;

    chronoReturn(); return(bs);
}

/**************************************************************/
/* note that if b0=0 we do not read bs->down 
 * new code: loop down, recurse right
 * 2014_11 this code is not use, i do not know why, probably it was not fully tested
 * but it would be faster
 */
#ifdef JUNK
static BS  bsTreeReadNew (BP bp,NODE nn,BS b0, BOOL doUnbreakComment, AC_HANDLE handle)
{
  NODEP r ;
  BS bs, bs0 ;
  BS bs1 ;

  r=bp->n+nn ;
  bs = bs0 = BSalloc() ;
  while (bs)
    {
      bs->up = b0 ; 
      bs->key=maybeSwapKEY(r->ck.key) ;
      
      if (bs->key<=_LastC)
	{
	  if(r->d2)
	    {
	      Bnode2str(BSbuffer,bp,&r);
	      bs->bt = BTalloc();
	      bs->bt->cp = strnew(BSbuffer, handle) ;
	      
	      /* special right hook for compatibility with earlier version*/
	      bs->right = r->d1  ?
		bsTreeRead (bp,r->d1,bs,doUnbreakComment, handle) : NULL;
	    }
	  if(    doUnbreakComment
		 && (bs1 = bs->down)
		 && (bs1->key == _NextC) )
	    bsTreeUnbreakComment(bs, bs1, handle) ;
	  /* this call may have changed r1 */
	}
      
      else  if(bs->key<=_LastN)
	{           /* In a number a BS cell corrresponds to 2 nodes */
	  nn=r->d2;
	  r = bp->n + nn ; /* Thus the right hook
			      will start from the second cell */
	  bs->n.i  = maybeSwapInt(r->ck.x.i) ;
	  if(0 && bs->key == _continuationKey)
	    {
	      bsTreeReadCont(&bs, handle) ; 
	    }
	}
      /* right hook */
      if (r->d2)
	bs->right = bsTreeRead (bp,r->d2,bs,doUnbreakComment, handle) ;
      else
	bs->right = NULL ;
      /* down hook */
      if (b0 && r->d1)
	{
	  b0 = bs ;
	  bs->down = BSalloc () ;
	  r=bp->n+r->d1 ;
	  bs = bs->down ;
	}
      else
	bs = 0 ;
    }
  return(bs0);
}
#endif
/**************************************************************/
                 /* replaces *bsp by the top
                    of the corresponding tree
                    freeing the relevant cells
                    */

static void  bsTreeReadCont(BS *bsp, AC_HANDLE handle)
{
  NODE nn;
  BP bp;
  BS bs = *bsp,
     up = bs->up,
     bs1 = 0,
     bs2 = 0;

  chrono("bsTreeReadCont") ;

  bp = blockGetContinuation(bs->n.i) ;
  nn= bp->hat.top ;
  bs1=bsTreeRead (bp,nn,(BS)0,TRUE, handle);

  if(bs1)
    bs2 = bs1->right;
  if(!bs2)
       messcrash("Anomaly  bs2 =0 in bsTreeReadCont " );
  if(!up || bs->down || bs->right )
       messcrash("Anomaly in bs in bsTreeReadCont ");
  bs2->up = up;
  BSfree(bs);
  BSfree(bs1);
  *bsp = bs2;
  chronoReturn();
}

/**************************************************************/

static BP  bsTreeFindBlock  (BS  bs, KEY *kkp, AC_HANDLE handle)
{
  KEY key = bs->key, kk, k2 ;
  int
    freecells, n1,
    ncells   = bs->size,	/* number of cells needed to store bs */
    isold    = iskeyold (key),
    ispruned = FALSE;
  BP bp;
  BS bs2;
  NODE nn ;
  NODEP r;

  chrono("bsTreeFindBlock") ;

  if (isold || ncells > BNODEMAX/2)
    kk = key ;
  else				/* find an existing block */
    kk = blockfriend (key) ;

  if (isold || iskeyold(kk))
    Bget (&bp, &nn, kk) ; /* *bp,*nn becomes the origin of the branch */
  else
    { blockpinn (kk, &bp) ;
      nn = 0 ;
      Binit (bp) ;
    }
				/* bp is pinned now */
  if (isold && nn)
    { Bprune (bp, nn, 0) ;
      ispruned = TRUE ;
    }

	/**** 940802: IMPORTANT CHANGE HERE (I think) *****/
  if (bp->hat.top)		/* something left in the block */
    { freecells = Bfreelen (bp) ;
      if (ncells > freecells)
	{ if (ncells < BNODEMAX/2) /* stay shared */
	    { r = bp->n + bp->hat.top; /* find the top object */
	      k2 = kk ;

	      kk = maybeSwapKEY(r->ck.key) ;

	      blockunpinn (k2) ;
	      Bget (&bp, &nn, kk) ;
	      nn = bp->hat.top ;
	      n1 = Bbranchlen (bp, nn) ;
	      while (ncells + n1 > BNODEMAX/2)
		{ r = bp->n + nn ; nn = r->d1 ;
		  if (!nn) break;
		  n1 = Bbranchlen (bp, nn) ;
		}
	      if (nn)
		{ bs2 = bs ;
		  while (bs2->down) bs2 = bs2->down ;
		  bs2->down = bsTreeRead (bp, nn, bs2, FALSE, handle);
		  Bprune (bp, nn, TRUE) ;
		  ispruned = TRUE ;
		}
	    }

	  if (ispruned)
	    blockrewrite (kk) ;
	  else
	    blockunpinn (kk) ;

	  kk = key ;
	  blockreallocate (kk, &bp) ;
	  Binit (bp) ;
	}
    }

  *kkp = kk ;
  chronoReturn () ;
  return bp ;
}

/**************************************************************/
     /* Removes key from cache1, hence from disk. Needed by bsFuseObjects */
void bsTreeKill  (KEY key)
{
  BP bp;
  NODE nn ;

  chrono("bsTreeKill") ;

  if (!iskeyold(key))
    return ;

  Bget(&bp,&nn,key);   /* *bp,*nn becomes the origin of the branch*/
  if(nn)
    Bprune(bp,nn,0);

  if (bp->hat.top)   /* There is something else on the block */
    { NODEP r ; KEY kk ;

      blockunpinn(key);
      r=bp->n+bp->hat.top;   /*Find the top object*/
#ifdef ACEDB4
      kk=maybeSwapKEY(r->ck.key);
#else
      kk=r->ck.key;
#endif
      Bget(&bp,&nn,kk);

      blockSetEnd(bp) ;    /* discard anything further down */
      blockrewrite(kk) ;
    }
  else
    blockSetEmpty(bp) ;
 lexkill (key) ;           /* discard reference from lex1 */
 chronoReturn() ;
}

/**************************************************************/

/*
* this function recurses through a tree, copying it into a
* database block buffer.  
*
* Each call to Bnewnode() gets the address of the next NODEP in
* the buffer.  You call Bnewnode() and then fill it in with
* the value you want - you never see any call that looks like
* "write out the node".
*
* A NODEP is:
*	a short named "d1"
*		( d1 is the number of cells in the tree )
*	a short named "d2" 
*		( d2 is the number of cells deeper in the tree )
*	a union that contains either a KEY, a BSdata,
*		or just a character array.
* It happens that the short is 2 bytes and the union is 4 bytes,
* for a total of 8 bytes on disk per NODEP.
*
* Each NODEP contains one tag from the object tree.  But if
* the tag is one of the special tags that contains data (_Int, 
* etc) then there is a second NODEP to contain the data value.
* The d1, d2 fields of the second NODEP have no meaning, but
* may contain arbitrary data.
* 
*/

static void bsTreeStore2 (BP bp, KEY kk, BS bs)
{
  NODE nn,nn2;
  NODEP r;

  chrono("bsTreeStore2") ;
  if (bs->size > Bfreelen(bp))
    { bsTreeDump(bs) ;
      messcrash ("Overflow in bsTreeStore2 (%s): bs->size = %d, Bfreelen = %d\n",
		 name(kk), bs->size, Bfreelen(bp)) ;
    }

  nn=bsTreeSt2(bs,bp);
  nn2=nn;
  while(r=bp->n+nn2,r->d1) nn2=r->d1;
  r->d1=bp->hat.top;
  bp->hat.top=nn;

  bsTreePrune(bs);
  chronoReturn() ;
}

/**************************************************************/

static void bsTreeStoreContinuation  (KEY key, BP *bpp, BS  bs)
{ NODE nn ;
  chrono("bsTreeStoreContinuation") ;

  blockSetNext(bpp) ;
		/******* IMPORTANT CHANGE HERE ********/
  if ((*bpp)->hat.h.disk)	/* reused, not grabbed */
    { BLOCKHEADER h = (*bpp)->hat.h ; /* save disk header info */
      Binit (*bpp) ;
      (*bpp)->hat.h = h ; /* restore disk linking */
      (*bpp)->hat.top = 0 ;
      (*bpp)->hat.free = 1 ;
    }
  else
    Binit (*bpp) ;

  if (bs->size > Bfreelen(*bpp))
    { bsTreeDump(bs) ;
      messcrash ("Overflow in bsTreeStoreContinuation (%s): bs->size = %d, Bfreelen = %d\n",
		 name(key), bs->size, Bfreelen(*bpp)) ;
    }

  nn = bsTreeSt2 (bs, *bpp) ; /* 2 instructions needed si *bpp change */
  (*bpp)->hat.top = nn ;

  bsTreePrune (bs) ;
  chronoReturn () ;
}

/************************************************************************/

static NODE bsTreeSt2 (BS bs,BP bp)
{
  NODEP r = Bnewnode(bp) ;

  r->ck.key=maybeSwapKEY(bs->key);

  if(bs->down) 
    r->d1=bsTreeSt2(bs->down,bp);

  if((bs->key)<=_LastC)
    { 
      /* XQX
       * keys <= _LastC are type A objects
       *
       * In practice, this means any tag that can have more than 4 bytes
       * of data.  This means TEXT, comments, and several types of
       * dna or protein.  We just grab as many NODEP as we need to store
       * the data.
       */
      NODEP rnew = r ;  /* Bstr2node will grab additional nodes */
      Bstr2node(bs,bp,&rnew);
      /* comments  support right hooks in funny way */
      if(bs->right)
	rnew->d1 = bsTreeSt2(bs->right,bp) ;
    }
  else if((bs->key)<=_LastN)
    {	     /* grab a second node to write the Union */
      /* XQX
       * keys > _LastC and <= _LastN are data elements
       * We write one node to contain the tag (_Int, etc)
       * and a second node to contain the actual data.
       * ( TEXT fields don't come here - they are actually stored
       *   as type A objects.  All you see here is an object key.)
       */
      NODEP rnew = Bnewnode(bp) ;

      r->ck.key = maybeSwapKEY(bs->key) ;
      rnew->ck.x.i = maybeSwapInt(bs->n.i) ;

      r->d2 = (NODE) (rnew  - bp->n) ;
      if(bs->right)
	rnew->d2=bsTreeSt2(bs->right,bp);
    }
  else /* standard key */
    {
      /* XQX
       * All other keys that we might see here are tag nodes.
       */
      if(bs->right) 
	r->d2=bsTreeSt2(bs->right,bp);
    }

  return (NODE) ( r - bp->n) ;
}

/************************************************************************/
/* Next two routines recursively unpack/pack the timestamp information. */
/* As usual, they loop over a column of bs nodes, and recurse right.    */
/* In the new bs -> array packing, could be done as data are packed.    */
/************************************************************************/

/* XQX
* Timestamps are weird.  They are stored in the struct BS actually ON the
* data item to be flagged, but the disk representation does not
* contain a timestamp in each object.  Instead, the timestamp is
* inserted into the object tree at various places.  When you load
* the object from the disk, the timestamp nodes are removed and the
* timestamps are propagated to every node in the object tree.
* 
* timestamps are NOT a representation of the modification time.  They
* are a key of a _VUserSession session object.
*/

/* XQX
* unpack happens when loading from disk - The object tree has already
* been loaded and stored as a linked group of BS.  This walks the
* tree and eliminates _VUserSession nodes.  A _VUserSession node
* contains a timestamp, which is copied to "nearby nodes".
*
*/

static void bsUnpackTimeStamps (BS bs, KEY currStamp)
{
  while (bs)
    if (class(bs->key) == _VUserSession)
      { currStamp = bs->key ;
	if (bs->up->right == bs)
	  bs->up->right = bs->down ;
	else
	  bs->up->down = bs->down ;
	if (bs->down)
	  bs->down->up = bs->up ;
	{ BS condemned = bs ; bs = bs->down ; BSfree (condemned) ; }
      }
    else
      { bs->timeStamp = currStamp ;
	if (bs->right)				/* recursion */
	  bsUnpackTimeStamps (bs->right, currStamp) ;
	bs = bs->down ;
      }
}


/* XQX
* unpack happens when saving to disk - this function walks the
* object tree finding timestamps and inserts new nodes into
* the tree wherever the traversal finds a timestamp different
* from the thing above it.
*
* currStamp is a _VUserSession key.
*/

static void bsPackTimeStamps (BS bs, KEY currStamp)
{
  while (bs)
    { if (bs->timeStamp != currStamp) /* add in a timestamp */
	{ BS new = BSalloc() ;
	  currStamp = bs->timeStamp ;
	  new->key = currStamp ;
	  new->up = bs->up ;
	  if (bs->up->right == bs)
	    bs->up->right = new ;
	  else
	    bs->up->down = new ;
	  bs->up = new ;
	  new->down = bs ;
	}
      if (bs->right)
	bsPackTimeStamps (bs->right, currStamp) ;
      bs = bs->down ;
    }
}

/************************************************************************/
/************************************************************************/

