/*  Last edited: Apr 16 18:17 1993 (mieg) */

/* $Id: block.h,v 1.2 2015/08/17 22:53:54 mieg Exp $ */

 
                   /*block.h                         */
                   /*public functions of blocksubs.c */
                   /* handling the cache             */
 
#ifndef DEFINE_BLOCK_H
#define DEFINE_BLOCK_H

#ifndef DEF_BP
#define DEF_BP
typedef void* BP ;
#endif

 void blockInit(void);
                 /*Allocates BLOCKMAX blocks and their control area*/

 void blockwrite(KEY k);     /*force write the relevant block*/
 void blockrewrite(KEY key);
 
 int blockpinn(KEY k, BP *p); /*reads and pinns the relevant block*/
 void blockreallocate(KEY kk,BP *p);
 BOOL blockNext(BP *bpp); /* After blockGet, gives the continuations */
 BP blockGetContinuation(int nc) ; /* gives the nc continuation */
 void  blockSetNext(BP *bpp);
 void blockSetEnd(BP bp) ;
 void blockSetEmpty(BP bp) ;
 void blockunpinn(KEY k);     /*unpinns the relevant block*/
 
 void blockmark(KEY k);      /*marks  the relevant block as modified*/

 int blockMax (void)   ;
 void blockSave(BOOL whole);       /*writes everything (half) back to disk*/
 void blockavail(int *used,int *pinned,int *free,int *modif) ;
 void blockStatus(int *nc1rp, int* nc1wp, int *ncc1rp, int* ncc1wp) ;

void blockshow(void);       /*scrmess the cache content*/
KEY blockfriend(KEY key) ;  /* gives a loaded key of same class */ 

void  blockShutDown (void) ;
#endif 



