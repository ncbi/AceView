/*  File: blocksub.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  9 15:11 1998 (fw)
 *	-	everywhere a "if(!(*q)->addr)" test is performed within
 *		a CHECKBLOCK #define, use   *q = KEY2LEX(key) to ensure
 *		that "q" is not a dangling pointer reference (see objcache.c too)
 * * Feb 27 13:19 1996 (srk)
 * Last edited: Feb 27 13:19 1996 (srk)
 * * Jan 22 10:24 1992 (mieg): 
    Session object should not be shared.
    Session and bat obj should be modified en place.
 * Created: Wed Jan 22 10:24:47 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: blocksub.c,v 1.7 2014/02/20 00:10:55 mieg Exp $ */


      /***************************************************************/
      /**  File blocksubs.c :                                       **/
      /**  Handles the cache of the ACeDB program.                  **/
      /***************************************************************/
      /***************************************************************/
      /*                                                             */
      /*  Twelve routines are public :  avail,  rewrite, reallocate, */
      /*   show, mark&save, init,  get/write, pinn/unpinn.           */
      /*              blockfriend, blockNext.                        */
      /*      The cache is an area which contains BLOCKMAX blocks.   */
      /*                                                             */
      /* The cache has two piles controlled by get and pinn/unpinn : */
      /*   1 : the "usedblocks", a LILO stack, on top of which each  */
      /* block is put when  "get" (implying a load the first time    */
      /* or a pop, or a new). These blocks are lost, LILO way, when  */
      /* the cache is full, automatically.                           */
      /*   2 : the "pinnedblocks". Those blocks stay in the cache.   */
      /* This pile is controlled explicitely by blockpinn and unpinn */
      /* Pinns implies a get but unpinn just moves back the block to */
      /* the "usedblocks" list.                                      */
      /*                                                             */
      /* If an object is too large to fit in a single block it       */
      /* is written in a set of blocks linked to the rigth an        */
      /* handled as a whole by the stacks .                          */
      /*                                                             */
      /* blockmark specifies that a block has been modified.         */
      /* blockwrite forces a copy of block to disk, if marked.       */
      /* It has  no action on the cache.                             */
      /* It can be called explicitely. It is invoked by save.        */
      /* blockSave saves the marked block to disk. It is called      */
      /* from the main menu.                                         */
      /*                                                             */
      /* Blockshow is for debugging. Blockavail gives the status.    */
      /*                                                             */
      /* The actual disk handling is done in the separate file       */
      /* disksubs.c invoked only by the static routines load/unload. */
      /*                                                             */
      /*         R.Durbin & J.Thierry-Mieg.                          */
      /*                    last modified    7/3/1991 by JTM.        */
      /*                                                             */
      /***************************************************************/

#undef CHRONO
/* n.b. The TURBOC compiler warns on each use of KEY2LEX(k) that there is
a possible loss of significant digits, which I dont understand since
an LEXP lxi=KEYTOLEX(k) is by itself a huge pointer which can hold
a long, KEY2LEX is defined in ACEDB.h */

#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct block *BLOCKP ;
typedef BLOCKP BP ;

#include "acedb.h"
#include "lex_bl_.h"
#include "disk_.h"   /* used by block load/unload*, disk_ used by blockOrder */
#include "block.h"
#include "disk.h"    /* for function prototypes */

    /*The following are default values. 
      They are overridden by cachesize.wrm in getCacheSize() below 
      */
static int BLOCKMAX =  1000 ; /* default number of blocks in the cache */
int MAXKNOWNCACHE =  1000 ;    /* default number of cacheEntries in the cache, used by objcache.c */
extern int INITIALDISKSIZE ;

static void   blockload   (LEXP *q, ALLOCP v, KEY key);
static ALLOCP blockgrab   (KEY key) ;
static void   blockfree   (ALLOCP v);
static void   blockunload (ALLOCP v);
static BOOL   blockisload (ALLOCP v, KEY k, LEXP *q);
static void   blockisunld (ALLOCP v);
static int    blocknewdsk (BP p);
static BOOL   blockNextKey (BP p, KEY *kp) ;
static ALLOCP freeblocks = 0 ,
              pinnedblocks = 0 ,
              usedbltop = 0 ,
              usedblend = 0 ;

/*
#define TESTBLOCK
#define CHECKBLOCK
*/

extern void invokeDebugger(void) ;
static int npinned = 0, nfree = 0, nused = 0 ;

/**************************************************************/
static int  nc1r = 0, nc1w = 0, ncc1r = 0, ncc1w = 0 ;

void blockStatus(int *nc1rp, int* nc1wp, int *ncc1rp, int* ncc1wp)
{
  *nc1rp = nc1r ; *nc1wp = nc1w ; *ncc1rp = ncc1r ; *ncc1wp = ncc1w ;
}

/**************************************************************/

BOOL blockCheck (void)
{ char *vp, *vp0 = 0 ;
  ALLOCP u,v ;
  Associator table = 0 ;
  int m, n = 0 , nu = 0, nf = 0, np = 0 ;
  int list, badlist ;
  static char* listNames[] = {"","FREE","PINNED","USED"} ;
  BOOL result = TRUE ;

  table = assCreate () ;
  
  list = 1 ; m = n ;
  for (v = freeblocks ; v ; v = v->next)
    { ++n ; nf++;
      if (!assInsert (table, v, vp0 + list))
	goto loop_found ;
      if (v->next && v->next->up != v)
	{ messout (
 "blockCheck: bad link at %d in %s list 0x%x -> 0x%x back to 0x%x",
		   n-m, listNames[list], v, v->next, v->next->up) ;
	  result = FALSE ;
	}
    }
  list = 2 ; m = n ;
  for (u = pinnedblocks ; u ; u = u->next)
    for (v = u ; v ; v = v->right)
      { ++n ; np++;
	if (!assInsert (table, v, vp0 + list))
	  goto loop_found ;
	if (v->next && v->next->up != v)
	{ messout (
 "blockCheck: bad link at %d in %s list 0x%x -> 0x%x back to 0x%x",
		   n-m, listNames[list], v, v->next, v->next->up) ;
	  result = FALSE ;
	}
      }
  list = 3 ; m = n ;
  for (u = usedbltop ; u ; u = u->next)
    for (v = u ; v ; v = v->right)
      { ++n ; nu++ ;
	if (!assInsert (table, v, vp0 + list))
	  goto loop_found ;
	if (v->next && v->next->up != v)
	{ messout (
 "blockCheck: bad link at %d in %s list 0x%x -> 0x%x back to 0x%x",
		   n-m, listNames[list], v, v->next, v->next->up) ;
	  result = FALSE ;
	}
      }

  if (n != BLOCKMAX)
    { messout ("blockCheck found %d blocks instead of %d\n"
	       "%d free, %d used %d pinned",
	       n, BLOCKMAX, nf, nu,np) ;
      result = FALSE ;
    }

  assDestroy (table) ;
  return result ;

 loop_found:
  assFind (table, v, &vp) ;
  badlist = vp - vp0 ;
  messout (
"blockCheck: loop found from 0x%x in %s list to %s list - %d'th item",
	   v, listNames[list], listNames[badlist], n) ;
  assDestroy (table) ;
  return FALSE ;
}

int blockMax (void)  /* for status */
{   
  return (BLOCKMAX * sizeof(BLOCK)) >> 10 ;           ;
}

/*************************************************************/
   /* import cache size or used defaults */

static FREEOPT cacheSizeOptions[] =
{
  {3, "Cache"},
  {'1', "CACHE1"},
  {'2', "CACHE2"},
  {'3', "DISK"},
  } ;

static void getCacheSize()
{
  int n  ; KEY option ;
  char *cp = sessionFilName ("wspec/cachesize","wrm","r") ;
  FILE *fil = cp ? filopen(cp, 0, "r") : 0 ;
    
  if(!fil)
    return ;

  freespecial ("\n\t\"\\") ;

  while(freeread(fil))
    if (freekey(&option, cacheSizeOptions))
      switch (option) 
	{
	  case '1':
	  freenext() ; freestep ('='); freenext() ;
	  if (freeint(&n))
	    BLOCKMAX = n ;
	  else
	    messerror("In semkey.wrm, I cannot read the number after CACHE1") ;
	  break ;
	  case '2':
	  freenext() ; freestep ('='); freenext() ;
	  if (freeint(&n))
	    MAXKNOWNCACHE = n ;
	  else
	    messerror("In semkey.wrm, I cannot read the number after CACHE2") ;
	  break ;
	  case '3':
	  freenext() ; freestep ('='); freenext() ;
	  if (freeint(&n))
	    INITIALDISKSIZE = n;
	  else
	    messerror("In semkey.wrm, I cannot read the number after CACHE2") ;
	  break ;
	}
  filclose(fil) ;
}

/********************************************************************/
                      /* constructs the Cache */
                      /*Allocates m<=n blocks and their control area*/
                      /* returns the number m really allocated*/
                      /*Called at least once on entering the program*/
                         /*     Public;  Calls nothing */
static AC_HANDLE blockHandle = 0 ;

void blockDoInit(BOOL extend)
{ int n ;
  void *v,*w;
  ALLOCP p;
  BP q;
  register unsigned int i;
  register int m;
  static BOOL firstpass = TRUE ;
  static int blocmax ;

  if (!extend && !firstpass)
    return ;
  
  if (firstpass)
    { freeblocks=
	pinnedblocks=usedbltop=usedblend=(ALLOCP )NULL;
      getCacheSize() ;
      firstpass = FALSE ;
      blocmax = BLOCKMAX ;
      blockHandle = handleCreate () ;
    }
  else
    BLOCKMAX += blocmax ;
  
/*  printf("blockDoInit now: %d\n", BLOCKMAX) ; */
  n = blocmax ;
  i=(unsigned int)n*(unsigned int)sizeof(struct alloc);
  /* v =  messalloc(i) ;  handleAlloc(0, blockHandle, i) ; */
  v = halloc (i, blockHandle) ; 

  i=(unsigned int)n*(unsigned int)sizeof(BLOCK);
  /* w = messalloc(i) ; */
  w = halloc (i, blockHandle) ; 

  memset((char *)w,0,(int)i);
  m=n;p=(ALLOCP)v; q=(BP)w;
 
  while(m--) 
    { p->next=freeblocks; 
      p->up=p->right = 0;
      p->p=q;
      p->ismodified=0;
      p->ispinned=0;
      if (freeblocks) freeblocks->up=p;
      freeblocks=p;
      p++;
      q++;
      nfree++;
    }
}

void blockInit(void)
{ blockDoInit (FALSE) ;
}

static Array blockSaveArray = 0 ;
void blockShutDown (void)
{
  arrayDestroy (blockSaveArray) ;
  messfree (blockHandle) ;
}

/**************************************************************/

static int blockOrder(const void *a, const void*b)
{
 return 
 (int)(*(const DISK *)
           (((const ALLOCP)a)->p) - *(const DISK *)(((const ALLOCP)b)->p) ) ;
}

/**************************************************************/

                 /* Called when saving a session 
		  * or when there is no more freeblocks.
		  * Gathers all or a fraction of flagged blocks
		  * and saves them in order, to try better the
		  * disk performances. 
                  * Returns the number of blocks written
		  *  Public;   Calls blockunload 
		  * Called by saveAll and blockGrab
		  */

void blockSave(BOOL whole)
{
  Array a = blockSaveArray ;
  register int i = 0 , m = 0, n = 0, max ;
  ALLOCP v ;
  
 /* check on calls that isWriteAccess == TRUE */
  
  chrono("blockSave");
  
  v = usedblend ;
  if (!v)
    return ;
  
  if (!whole)
    {
      m = 0 ;
      while(m++, (m < 100) && (v = v->up))
	if (v->ismodified) n++ ;
      m /= 2 ;  
    }
  if (n || whole)
    {
      blockSaveArray= a = arrayReCreate(a, 200, ALLOCP) ;
      v = usedblend ; n = 0 ; i = 0 ;
      while (v)
	{
	  if(v->ismodified)
	    array(a,i++,ALLOCP) = v ;
	  v = v->up ;
	  if(!whole && n++ > m)
	    break ;
	}
      
      arraySort(a,blockOrder); 
      
      max = arrayMax(a) ;
      for(i=0;i<max;i++)
	blockunload(arr(a,i,ALLOCP)) ;
    }
  chronoReturn() ;  
}

/**************************************************************/

/* Called to accelerate parsefile */
typedef struct preLoadStruct {KEY key; ALLOCP v;} PRELOAD ;

static int preLoadOrder (const void *a, const void *b)
{
  const PRELOAD *pa = (const PRELOAD *)a, *pb = (const PRELOAD *)b ;
  ALLOCP va = pa->v, vb = pb->v ;
   return 
     (int)(*(DISK *)(va->p) - *(DISK *)(vb->p)) ;
}

#include "bs.h"
void blockPreLoad (KEYSET ks)
{
  static Array a = 0 ;
  register int i = 0 , j = 0, max = 0 ;
  ALLOCP v ;
  KEY *kp, key ;
  LEXP q ; 
  PRELOAD *pl ;
 /* check on calls that isWriteAccess == TRUE */
  
  chrono("blockPreLoad");

  max = arrayMax(a) ;
  a = arrayReCreate(a, max, PRELOAD) ;

  for(i=0, j = 0, kp = arrp(ks, 0, KEY) ; i<max; i++, kp++)
    {
      key = *kp ;
      q=KEY2LEX(key) ;
      v = q->addr ;
      if (v)
	{
	  pl = arrayp (a, j++, PRELOAD) ;
	  pl->key = key ; pl->v = v ;
	}
    }
  arraySort(ks,preLoadOrder); 
  if (j)
    {
      OBJ obj = 0 ;

      pl = arrp (a, 0, PRELOAD) - 1 ;
      while (pl++, j--)
	if ((obj = bsCreate (pl->key))) bsDestroy (obj) ;
    }
  chronoReturn() ;  
}

/******************************************************************/
                          /* Gives the cache status*/
                          /*    Public;  Calls nothing*/

void blockavail(int *used,int *pinned,int *free,int *modif)
{
  register ALLOCP v, w;
  register int i;

  i=0;v=usedbltop;
    while((w=v))
      { while ((w=w->right)) i++;
	v=v->next; i++;}
    *used=i;
  i=0;v=pinnedblocks;
    while((w=v))
      { while ((w=w->right)) i++;
	v=v->next; i++;}
    *pinned=i;
  i=0;v=freeblocks;
    while(v) {v=v->next; i++;}
    *free=i;
  i=0;v = usedblend ;
    while((w=v))
      { if(v->ismodified)
	  { while ((w=w->right)) i++;
	    i++ ;
	  }
	v=v->up;
      }
    *modif=i;
}

/****************************************************************/
    /* For debugging, prints out the content of the cache
       Public;  Calls nothing
    */
#ifdef TESTBLOCK
#include "graph.h"
/**************************************************************/
             /*blockshow_keys, static, called by blockshow*/
static void blockshwks(ALLOCP v,int *line)
{ int i=1;
 KEY key=0;
 char *cp ;

 if(!(v->ismodified))
 while(blockNextKey(v->p,&key))
   {if(i>=80)                
      { (*line)++;i=1; }
   cp = messprintf (" : %s ",name(key)) ;
   graphText(cp,i,*line);
   i += strlen(cp) ;
  }
 (*line)++;
}
/**************************************************************/
             /*blockshow_keys, static, called by blockshow*/
static void blockshmod(ALLOCP v,int *line)
{  int i=1;
 KEY key=0;
 char *cp ;

 if(v->ismodified)
 while(blockNextKey(v->p,&key))
  {if(i>80)                   
     { (*line)++;i=1;}
   cp = messprintf (" : %s ",name(key)) ;
   graphText(cp,i,*line);
   i += strlen(cp) ;
  }
 (*line)++;
 }

#endif

/**************************************************************/

void blockshow(void)
{
#ifdef TESTBLOCK
  register ALLOCP v ;
  int line = 1 ,bused,bfree,bpinned,bmodif;
  int nc1r, nc1w, ncc1r, ncc1w ;

  graphClear() ;
  blockStatus (&nc1r, &nc1w, &ncc1r, &ncc1w) ;
  graphText(messprintf(
  "Cache1 usage : %d obj served, %d saved, %d read from disk, %d saved to disk",
               nc1r, nc1w, ncc1r, ncc1w),
            1,line++);
  line++ ;
  blockavail(&bused,&bpinned,&bfree,&bmodif);
  graphText(messprintf(
  "Cache1 present content: %d used, %d modified, %d pinned, %d freeblocks",
               bused,bmodif,bpinned,bfree),
            1,line++);

  graphText("Modified blocks  ",1,line) ; line++ ;
  for (v = usedblend ; v ; v = v->up)
    if (v->ismodified)
      blockshmod(v,&line) ;

  line += 2 ; graphText("Pinned blocks  ",1,line) ; line++ ;
  for(v = pinnedblocks ; v ; v = v->next)
    blockshwks(v,&line);

  line += 2 ; graphText("Used blocks  ",1,line) ; line++ ;
  for(v = usedbltop ; v ; v = v->next)
    blockshwks(v,&line);

  graphTextBounds (100, 5+line);
  graphRedraw() ;
#endif
}

/*****************************************************************/
      /*Returns a pinned block to the top of the usedblock line*/
                   /*      Public;   Calls nothing.  */
void blockunpinn(KEY kk)
{
 LEXP q=KEY2LEX(kk);
 ALLOCP u,v,w;

 chrono("blockunpinn");
 v=q->addr;
 if (v->ismodified) nc1w++ ; /* statistics */
 if(!v->ispinned)
   messcrash(
	     "Unbalanced pinning found in blockunpinn \n%s",
	     name(kk));
             
 if(--(v->ispinned))
   {     chronoReturn() ; return;}
             /* the block is pinned in some other way*/
 if(v->ismodified & 2)
   blockunload(v); /*a blockwrite is pending*/
   
                                    u=v->up;
                                    w=v->next;
                                    if(!u) pinnedblocks=w;
                                       else u->next=w;
                                    if (w) w->up=u;
                                    v->next=usedbltop;
                                    v->up=(ALLOCP)NULL;
                                    if (usedbltop) usedbltop->up=v;
                                    usedbltop=v;
                                    if(!usedblend) usedblend=v;

   npinned-- ; nused++ ;
  if (usedblend && usedblend->ispinned) invokeDebugger() ;

 chronoReturn() ; return ;
}
/*****************************************************************/
                      /* Clears the disk address of an object */
                      /* and reallocate a new empty block for it */
                      /* This sub does not touch the former block */
                      /*  Public ;  Called by BSstore     */
                      /*            Calls blockpinn */
void blockreallocate(KEY kk,BP *p)
{
 LEXP q=KEY2LEX(kk);
 
 chrono("blockreallocate");

 q->addr=0;
 q->dlx.dk=0;
 blockpinn(kk,p) ;

 chronoReturn() ;
}

/********************************************************************/
               /*pinns a block, it must then be unpinned explicitely*/
                      /*before the memory can be recovered*/
                      /* return the number of blocks used */
                      /*    Public;   Calls blockget */
static ALLOCP activeBlock = 0 ;
int blockpinn(KEY kk,BP *p)
{
  LEXP q=KEY2LEX(kk);
  ALLOCP u,v=q->addr,w;
  int n = 0;
  chrono("blockpinn");
  nc1r++ ;
  if (!v)
    { v = blockgrab(kk) ; /* get a new block */
/* NOTE: blockgrab calls blockfree which calls blockisunld which calls */
/* blockNextKey which calls KEY2LEX, which can extend the lexi2 array and */
/* leave q as a dangling pointer. Let this be a dire warning to all those  */
/* who have sold their soul to the devil whose name is FORTRAN. Arrays and */
/* pointers make a volatile mix. Here we re-call KEY2LEX to fix the problem */
      q = KEY2LEX(kk);
    }
  else                /* get it from used list */
    { if (v->ispinned)
	messcrash("Double pinning of %s",name(kk));
      u=v->up;
      w=v->next;
      if(!u) 
	usedbltop=w;
      else
	u->next=w;
      if(w)
	w->up=u;
      else
	usedblend=u;

      if (usedblend && usedblend->ispinned) invokeDebugger() ;
      nused-- ;

    }
  v->ispinned++ ;

  npinned++ ;

  v->next = pinnedblocks;
  v->up = 0 ;
  if (pinnedblocks)
    pinnedblocks->up = v;
  pinnedblocks = v;
 
  if(!q->addr)
    {
      if(q->dlx.dk)
	{
	  blockload(&q,v,kk) ;
	  if (q->addr != v)
	    messcrash("blockload failed key %d %s",
		      kk,name(kk)) ;
	}
      else
	{
	  v->ismodified=1;
	  q->addr = v;
	  v->p->h.key = kk ;
	}
    }

  *p=v->p;
  
  activeBlock = v ;
  n = 1;        /* counts the cache blocks used by this object */
  while ((v = v->right)) n++ ;
  
  chronoReturn() ;
  return n ;
}

/*********************************************************************/

BOOL blockNext(BP *bpp)
{ 
 if (activeBlock->p != *bpp)
   messcrash("blockNext called out of context");
 if ((activeBlock = activeBlock->right))
   {
     *bpp = activeBlock->p ;
     chrono ("blockNextFound") ;
     chronoReturn() ;     return TRUE ;
   }
 return FALSE ;
}

/*********************************************************************/

BP blockGetContinuation(int nc)
{ register ALLOCP v =  activeBlock ;
  while(v->up && (v->up->right == v))
    v = v->up ;

  while(nc--)
    { 
      if(v->right)
	v = v->right ;
      else
	messcrash("Cant find continuation %d of block %s",
		  nc, name(v->p->h.key));
    }
  activeBlock = v ;
  return v->p ;
}

/*********************************************************************/

void blockSetNext(BP *bpp)
{
  if (activeBlock->p != *bpp)
    messcrash("blockSetNext called out of context");
  if (!(activeBlock->right))
    activeBlock->right = blockgrab(0) ;
  activeBlock->right->up = activeBlock ;
  activeBlock = activeBlock->right ;
  *bpp = activeBlock->p ;
}

/*********************************************************************/

void blockSetEnd(BP bp)
{ 
  ALLOCP a;
  
  if (!activeBlock || activeBlock->p != bp)
    messcrash("blockSetEnd called out of context");
  a = activeBlock->right ; /* so a not on used or pinned list */
  activeBlock->right = 0 ;
  while(a)
    { diskfree(a->p) ;
      a->up = 0;
      a->ispinned = FALSE ;
      a->ismodified = FALSE ;
      a->next = freeblocks ;
      if (freeblocks)
	freeblocks->up = a;
      freeblocks = a;
      a = a->right ;
      freeblocks->right = 0 ; /* must be final */
    }
}

/*********************************************************************/
   /* Release the whole disk space */
void blockSetEmpty(BP bp)
{ 
  ALLOCP a,u,w ;
  
  if (!activeBlock || activeBlock->p != bp)
    messcrash("blockSetEmpty called out of context");
  a = activeBlock ;
  if (!a->ispinned)
    messcrash ("block must be pinned to call blockSetEmpty") ;

  u = a->up ;		/* unhook from pinned list */
  w = a->next;
  if (!u) pinnedblocks = w ;
  else u->next = w ;
  if (w) w->up = u ;
  
  while(a)
    { diskfree(a->p) ;
      a->up = 0 ;
      a->ispinned = FALSE ;
      a->ismodified = FALSE ;
      a->next = freeblocks ;
      if (freeblocks)
	freeblocks->up = a ;
      freeblocks = a ;
      a = a->right ;
      freeblocks->right = 0 ;	/* must be final */
    }
}

/*******************************************************************/
              /* to load a block from disk   : return 0 if success */
              /*  Static; called by blockget.                      */
              /*          Calls diskblockread.      */

static void blockload(LEXP *q, ALLOCP v, KEY key)
{ 
  DISK d=(*q)->dlx.dk;
  int nn = 0 ;
  
  chrono("blockload");
  ncc1r++ ;

  /* Read the first block from q->dlx.dkess */
  diskblockread(v->p,d) ;
  /* it may be shared */
  if (!blockisload(v,key,q))
    messcrash ("blockload loaded a wrong block: %x = %s:%s in place of %x = %s:%s\n",
	       v->p->h.key, className(v->p->h.key), name(v->p->h.key), 
	       key, className(key), name(key)) ;

  /* Read the continuations directly */
#ifdef CHECKBLOCK
  *q = KEY2LEX(key) ;
  if(!(*q)->addr)
    messcrash("key %d %s key2lex= %d", key, name(key), KEY2LEX(key)) ;
#endif
  
  while ((d = v->p->h.nextdisk))
    {
      v->right = blockgrab(key) ;
      v->right->up = v ;
      v = v->right ;
      nn++ ;
#ifdef CHECKBLOCK
      *q = KEY2LEX(key) ;
      if(!(*q)->addr)
	messcrash("key %d %s key2lex= %d", key, name(key), KEY2LEX(key)) ;
#endif
      /* they are private */
      diskblockread(v->p,d) ;
      if (key != v->p->h.key)
	messcrash(" blockload loaded a wrong block : %s\n",
		  name(key)) ;
    }
#ifdef CHECKBLOCK
  *q = KEY2LEX(key) ;
  if(!(*q)->addr)
    messcrash("key %d %s key2lex= %d", key, name(key), KEY2LEX(key)) ;
#endif
 chronoReturn() ; 
}

/*****************************************************************/
                     /*returns a free allocp address or NULL*/
                     /*  Static; called byblockload */
                     /*          calls blockfree */
static ALLOCP blockgrab(KEY kk)
{
  ALLOCP v;
  int recursive = 0 ;
  static KEY key ;
  chrono("blockgrab");
  
  if(recursive)
    messcrash("recursive call of blockgrab");
  recursive++ ;
  if (kk) key = kk ;
 blfree:
  if(!freeblocks && usedblend) 
    { if (isWriteAccess())
	blockSave(FALSE) ;
      blockfree(usedblend);
    }
  
  if(!freeblocks) 
    {
      if (isWriteAccess() || pickList[class(key)].protected)
	{ if (recursive++ < 2)
	  { blockDoInit (TRUE) ;
	  goto blfree ;
	  }
	else
	  messcrash("Error : blockgrab failure(%d).\n"
		    "All blocks are pinned. "
		    "Balance Pinn/Unpinn or Increase CACHE1 which is "
		    "defined in the configuration file cachesize.wrm.",
		    BLOCKMAX
		    ) ;
	}
      else
	{
	  if(!isWriteAccess())
	    { 
	      if (messQuery("The cache1 is full\n"
			    "Do you want write access ?"))
		{
		  if(sessionGainWriteAccess()) 
		    goto blfree ;	/* could get write access */
		}
	    }
	  /* still no write access */
	  
	  if(recursive++ < 2)
	    {
	      blockDoInit (TRUE) ;
	      goto blfree ;
	    }
	  else
	    messExit("Block grab : the cache1 is full "
		     "(see wspec/cachesize.wrm) "
		     "and write access is denied.") ;
	}
    }
  
  v=freeblocks;
  freeblocks=v->next;

  nfree-- ;
  if (v->ispinned || v->ismodified) invokeDebugger() ;
  if (freeblocks && freeblocks->ispinned) invokeDebugger() ;

  if (freeblocks) 
    freeblocks->up=(ALLOCP)NULL;
  v->right = v->next= v->up = (ALLOCP)NULL;
  v->ispinned = v->ismodified = recursive = FALSE ;
  memset(v->p,0,sizeof(BLOCK)) ;
  chronoReturn() ; return(v);
}

/*******************************************************************/
                        /*frees the relevant block*/
                        /*  Static; called by blockgrab */
                        /*          calls blockunload */
static void blockfree(ALLOCP v)
{                       /*allways unpinn before free*/
  ALLOCP u,w;


  if(!isWriteAccess())
    while(v && v->ismodified)
      v = v->up ;   /* buffer the modified blocks */
  
  if(!v)
    return ;

  chrono("blockfree");
  blockunload(v) ;
  blockisunld(v);
  
  u=v->up;
  w=v->next;
  if(!u) usedbltop=w;
  else u->next=w;
  if(w) w->up=u;
  else usedblend=u;
  v->next=freeblocks;
  v->up=(ALLOCP)NULL;
  if (freeblocks) freeblocks->up=v;
  freeblocks=v;

  nused-- ; nfree++ ;
  if (usedblend && usedblend->ispinned) invokeDebugger() ;
  if (freeblocks->ispinned) invokeDebugger() ;


  while ((w=v->right))
    { 
      freeblocks = w ;

      if (freeblocks->ispinned) invokeDebugger() ;
      nfree++ ;

      v->up = w;
      v->right = 0 ;
      w->up = 0;
      w->next = v;
      v = w;
    }
 chronoReturn() ;
}
/**************************************************************/
  /* to mark a key   as modified : returns 0 */
                      /*  Public ; calls nothing */
void blockmark(KEY k)
{
  LEXP q=KEY2LEX(k);
  ALLOCP v=q->addr ;

  chrono("blockmark");
  if(!v) 
    messcrash(" blockmarking an unloaded key :%s\n",
                            name(k));
  v->ismodified |= 1;
  if(!(v->ispinned))
    messcrash("WARNING : blockmarking an unpinned block %s",name(k));

 chronoReturn() ;
}

/**************************************************************/
               /* to rewrite a block sharing the addresses*/
                      /*  Public ; calls nothing */
                      /* called by BSstore */
void blockrewrite(KEY key)
{ LEXP q = KEY2LEX(key) ;
 ALLOCP v = q->addr;

 chrono("blockrewrite");

 blockisload(v,key, &q);               /*to share memory address*/
 v->ismodified |= 1 ; /* == blockmark ,  to ensure rewriting*/
 blockunpinn(key);

 chronoReturn() ;
}

/**************************************************************/
  /* to force-write a key to disk if it is modified : return 0 if success */
                      /*  Public ;  Calls blockunload*/
void blockwrite(KEY k)
{
  LEXP q=KEY2LEX(k);
  ALLOCP v=q->addr;

  chrono("blockwrite");
  if (v)      /*else key is not in memory*/
    {
      if(v->ispinned)
	{
	  messout("%s%s",
		  "Sorry : blockwrite invoking a pinned block,",
		  "\n this is not an error, just a delay");
	  v->ismodified |=2; /*will unload as soon as unpinned*/
	}
      else
	blockunload(v) ;
    }

 chronoReturn() ;
}


/**************************************************************/
             /* to unload a block to disk   : return 0 if success */
             /*        Static, called by blockfree write and save */
             /*                calls diskalloc and diskblockwrite */
static void blockunload(ALLOCP v)
{
  ALLOCP w ;
  BP p=v->p, q;
  KEY key = 0 ; /* to init blockNextKey */

  
  if(!v->ismodified) 
     return ;
  
  if (!isWriteAccess())
    return ;

  v->ismodified=0;

  chrono("blockunload");
  ncc1w++ ;

  if (!blockNextKey(v->p,&key))
    {    /* i.e. the block is really empty */
      diskfree(v->p) ;
      chronoReturn() ;
      return ;
    }
  v->p->h.key = key ;
  p = v->p ;
 
  diskalloc(p) ;
  blocknewdsk(p) ;
      
  while(v)
    { p = v->p ;
      if ((w = v->right))
	{ q = w->p ;
	  q->h.key = key ;
	  diskalloc(q) ;
	  p->h.nextdisk = q->h.disk ;
	}
      else
	p->h.nextdisk = 0 ;
  
      diskblockwrite(p) ;

      v = v->right ;
    }
  chronoReturn() ;
}

/**************************************************************/
/**************************************************************/
              /* block is loaded static, called by blockload*/
static BOOL blockisload (ALLOCP v, KEY k, LEXP *q)
{ LEXP q1 ;
  KEY key=0;
  BOOL found = FALSE ;

  while(blockNextKey(v->p,&key))
	{ q1 = KEY2LEX(key) ;
	  q1->addr = v ;
	  if(key==k)
	    found = TRUE ;  /*k has been found as hoped */
	}
  *q = KEY2LEX(k) ;  /* may move because of arrayExtend(lex) */
  return found ;
}
/**************************************************************/
              /* block is unloaded static, called by blockfree*/
static void blockisunld (ALLOCP v)
{ LEXP q ;
  KEY key=0;
  BP p = v->p ;

  while(blockNextKey(p,&key))
   { q = KEY2LEX(key) ;
     if (q->addr == v)
       { q->addr = 0 ; lex2clear(key) ; }
   }
}

/**************************************************************/
              /* block new disk address static, called by blockload*/
              /* returns the number of objects in the block */
static int  blocknewdsk (BP p)
{
 KEY key=0;
 register int i=0;
 register LEXP q;
 DISK d = p->h.disk;
/*
 if(p->h.key == _continuationKey)
   return 1 ;
*/
 while(blockNextKey(p,&key))
               {
		 i++;
                 q=KEY2LEX(key);
                 if(q->dlx.dk!=d)
                      {q->dlx.dk=d;
                       lexmark(class(key));
                       }
                 }
 return i;
}
/**************************************************************/
        /* Returns a loaded unpinned key of same class, or key.
	 * Looks only for keys of same class, already touched
	 * during the same session and still in the cache.
         * Called from bstree.c : bsTreeFindBlock 
	 *  only when key has never been saved before. 
	 */

KEY blockfriend(KEY key)
{
  register  int t = class(key);
  register ALLOCP v = usedbltop ; 
  KEY k = 0;  
  int i = 100 ; 
  if(!v)
    return key ;
    
   /* Session object should not be shared.
    * Suppose you did, 
    * may be in session 8 you share session 5 block, hence 
    * session 5 would be relocated in a block allocated by session 8.
    * If you then destroyed session 8 before session 5, the 
    * session 5 obj would be lost.
    *
    * Caveat: this also holds for Bat obj if at some point shared
    * arrays were introduced.
    */
  if(class(key) == _VSession  || class(key) == _VBat)
    return key ;

  chrono("blockfriend") ;

  while (i-- && v)
    {
      if(
	 !v->right &&
	 (   v->ismodified 
	  || (v->p->h.session == thisSession.session)
	  || (v->p->h.session == 0)
          )
	 && blockNextKey(v->p,&k) && (class(k) == t)   
	 )
	{ key = k ; break ; } 
      k = 0 ;
      v = v->next ;
    }
  chronoReturn() ; return(key) ;
}
/**************************************************************/
extern int Bnextkey(BP p, KEY *key) ; /* in bsubs */

              /* block is loaded static, called by blockload*/
static BOOL blockNextKey (BP p, KEY *kp)
{
  static BP bp ;
 
  if(*kp)
    { if(bp != p)
	messcrash("WARNING : bad call to blockNextKey \n %s",
		 name(*kp));
   }
  else
    bp = p ;

  switch( bp -> h.type )
    {
    case 'B'  :
      return Bnextkey(p,kp) ;
    default:
      if (*kp)
	return FALSE ;   /* single object by block */
      return (*kp= bp-> h.key) ? TRUE : FALSE ;
    }
}

/**************************************************************/
/**************************************************************/
