/*  Last edited: Dec  4 13:33 1998 (fw) */

/* $Id: bstools.c,v 1.3 2007/03/21 21:02:49 mieg Exp $ */

      /***************************************************************/
      /***************************************************************/
      /**  File bstools.c :                                         **/
      /**  Manipulates the bs objects.                              **/
      /***************************************************************/
      /***************************************************************/
      /*                                                             */
      /*        routines are public :                                */
      /*                                                             */
      /*         R.Durbin & J.Thierry-Mieg.                          */
      /*                    last modified  20/9/91 by JTM.           */
      /*                                                             */
      /***************************************************************/

#include "acedb.h"

#include "bs_.h"
#include "dump.h"
/* #define CHRONO */
#include "chrono.h"
#include "parse.h"

#include <ctype.h>

extern int MAXKNOWNCACHE ;   /* sizeof cache2, defined in blocksub*/


/**************************************************************/

/* Manage BS allocation
   hope they speed things up/prevent fragmentation
   in any case, they will help monitoring memory usage
   NB they must be paired.
*/

static AC_HANDLE BShandle = 0 ;
static Stack freeBSstack = 0 ;
static int nBSused = 0, nBSalloc = 0 ;        /* useful for debug */
BOOL BS_Cache_is_Full = FALSE ; /* used by w5/objcache.c */
int BS_Cache_limit = 0 ;
#define BSHOWSIZE sizeof(struct bshow)
BS BSalloc (void)       /* self managed calloc */
{ static int blocSize = 2048 ;
  BS p ;
  int n ;
 
  chrono ("bsAlloc") ;
  
  if (!freeBSstack)
    { 
      BS_Cache_limit = (MAXKNOWNCACHE << 10) / sizeof(struct bshow) ;
      BShandle = handleCreate () ;
      freeBSstack = stackHandleCreate (4*blocSize, BShandle) ;
    }

  if (stackEmpty (freeBSstack))
    { 
      if (nBSalloc + blocSize < BS_Cache_limit)
	{ n = blocSize ; 
	if (blocSize < BS_Cache_limit/4) blocSize *= 2 ;
	}
      else
	n = BS_Cache_limit / 8 ;
      nBSalloc += n ;
      p = (BS) halloc (n * BSHOWSIZE, BShandle) - 1 ;
      while (p++, n--)
	push (freeBSstack,p,BS) ;
    }

  p = pop (freeBSstack,BS) ;
  memset (p, 0, BSHOWSIZE) ;

  if (++nBSused > BS_Cache_limit)
    BS_Cache_is_Full = TRUE ;
  chronoReturn() ;
  return p ;
}

void BSfree (BS bs)
{ chrono("BSfree") ;
  if (bs->bt)
    BTfree(bs->bt) ;
  push (freeBSstack,bs,BS) ;

  if (--nBSused < BS_Cache_limit)
    BS_Cache_is_Full = FALSE ;
  chronoReturn() ;
}

void BSstatus (int *used, int *alloc, int *memp)
{ *used = nBSused ; *alloc = nBSalloc ; *memp = nBSalloc * BSHOWSIZE ;
}

void BSshutdown (void)
{
  messfree (BShandle) ;
}
/**************************************************************/
/**************************************************************/

         /* idem for BT */
AC_HANDLE BThandle = 0 ;
static Stack freeBTstack = 0 ;
static int nBTused = 0 , nBTalloc = 0 ;    /* useful for debug */

BT BTalloc (void)       /* self managed calloc */
{
  static int blocSize = 2048 ;
  BT p ;
  int i ;

  if (!freeBTstack)
    {
      BThandle = handleCreate () ;
      freeBTstack = stackHandleCreate (4*blocSize, BThandle) ;
    }
  if (stackEmpty (freeBTstack))
    { p = (BT) halloc (blocSize * sizeof (struct btext), BThandle) ;
      for (i = blocSize ; i-- ; ++p)
        push (freeBTstack,p,BT) ;
      nBTalloc += blocSize ;
      blocSize *= 2 ;
    }
  p = pop (freeBTstack,BT) ;
  memset (p, 0, sizeof (struct btext)) ;
  ++nBTused ;
  return p ;
}

void BTfree (BT bt)
{
  if(bt->cp)
    bt->cp = 0 ; /* do not messfree, it belongs on 
		  * the obj->handle and the user may be referencing it
		  */
  push (freeBTstack,bt,BT) ;
  --nBTused ;
}


void BTstatus (int *used, int *alloc, int *memp)
{ *used = nBTused ; *alloc = nBTalloc ; *memp = nBTalloc * sizeof (struct btext) ;
}

void BTshutdown (void)
{
  messfree (BThandle) ;
}

/******************************************************************/
/******************************************************************/
static void   bsTreePrune2(BS bs)
{
/* old code- double recursion
  if(!bs)return;
  bsTreePrune2 (bs->right);
  bsTreePrune2 (bs->down);
  BSfree (bs);          
  bs=0;
************
New code, down loop, right recursion: much faster
*/
  BS bs1 ;

  while (bs)
    { if (bs->right && bs != bs->right) /* happens on purpose in models */      
	bsTreePrune2 (bs->right) ;
      bs1 = bs ; bs = bs->down ;
      BSfree(bs1); 
    }
}

/*******************/

void   bsTreePrune(BS bs)
{
  if(!bs)return;
  
  chrono("bsTreefree") ;

  if(bs->up)          /* unhook */
    {
      if (bs->up->down == bs)
	bs->up->down = 0 ;
      else if (bs->up->right == bs)
	bs->up->right = 0 ;
      else 
	messcrash("double link error in bsTreePrune") ;
      bs->up = 0 ;
    }

  if (bs) bsTreePrune2(bs);

  chronoReturn() ;

}

/************************************************************************/
/************************************************************************/
       /* BS Copying system */
/*************************************************************************/

static BS bsDoTreeCopy (BS bs, BS bsCopyUp, AC_HANDLE h)
{ BS bsCopy = BSalloc(), bsCopyTop = bsCopy ;
/* I BSalloc out of loop to singularise bsCopyTop */

  while (bs)
    { *bsCopy = *bs ;
      bsCopy->up = bsCopyUp ;
                                       /* local work */
      if(bs->bt)
	{ bsCopy->bt = BTalloc() ;
	  *bsCopy->bt = *bs->bt ;
	  if(bs->bt->cp) 
	    bsCopy->bt->cp = strnew (bs->bt->cp, h) ;
	}
                                      /* right recursion */
      if (bs->right)
	bsCopy->right = bsDoTreeCopy (bs->right, bsCopy, h) ; 

                                      /* down loop */
      if (!(bs = bs->down))
	break ;

      bsCopyUp = bsCopy ;
      bsCopy->down = BSalloc () ;
      bsCopy = bsCopy->down ;
    }
      
  return bsCopyTop ;
}

/**************************/
    /* copies a branch of a tree */
BS bsTreeCopy(BS bs, AC_HANDLE h)
{ return bs ? bsDoTreeCopy (bs, 0, h) : 0 ;
}

/******************************************************************/
   /* copies the whole b object */
BOOL bsCarbonCopy (KEY new, KEY old)
{ static Stack s = 0, s1 = 0 ;
  char *cp ;

#ifdef READONLY
  return FALSE ;
#else

  s = stackReCreate(s, 300) ;
  s1 = stackReCreate(s1, 300) ;

  dumpKey(old, 0, s) ;

          /* clean up new object */
  pushText(s1, messprintf("-D %s : \"%s\"", className(new), name(new))) ;
  catText(s1, "\n\n") ;
          /* fill it back */
  catText(s1, messprintf("%s : \"%s\"", className(new), name(new))) ;
  catText(s1, "\n") ;
  
  cp = stackText(s, 0) ;
  while (*cp++ != '\n') ;    /* skip name(key) */
  catText(s1, cp) ;       
  catText(s1, "\n") ;
  catText(s1, "\n") ;

  return  
    parseBuffer(stackText(s1,0), 0) ;
    
#endif
}

/*********************************/
/* called from objcache status function */
int  bsTreeSize (void *vbs)	  /* counts the number of in the branch */
{ 
  int n = 0 ;
  BS bs = (BS) vbs ;
  if (!bs) return 0 ;

  while (bs) /* down loop, right recursion */
    { if (bs->right) n +=  bsTreeSize (bs->right) ;
      n++ ;
      bs = bs->down ;
    }
  return n ;
}

/*****************************************************/

void bsTreeDump (BS bs)
{
  static int depth = 0 ;
  int i ;

  while (bs)
    { 
      for (i = 2*depth ; i-- ;)
	fputc (' ', stderr) ;
      /* keep stderr here since it is for debugging */
      fprintf (stderr,"%8lx : k.%d(%s)   up.%lx down.%lx right.%lx \n", 
	       (unsigned long) bs,
                               bs->key,
	                       name(bs->key), 
	       (unsigned long) bs->up, 
	       (unsigned long) bs->down, 
	       (unsigned long) bs->right) ;
      
      if (bs->right)
	{ ++depth ;
	bsTreeDump (bs->right) ;
	--depth ;
	}
      /* was:  if (bs->down)    bsTreeDump (bs->down) ; */
      bs = bs->down ; /* mieg: replaced by the while(bs) loop */
    }
}

/************************************************************************/
/************************************************************************/


