/*  Last edited: Apr 25 16:14 1995 (srk) */

/* $Id: b_.h,v 1.2 2010/06/15 17:54:49 mieg Exp $ */

                    /* b_.h                            */
                    /* definition of the B block       */
                    /* private to bsubs and b2bs.c     */

/***************************************************************************/
#ifndef DEFINE_B__h
#define DEFINE_B__h 1

#ifndef DEF_BP 
#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct bblock *BBLOCK ;
typedef BBLOCK BP ;
#endif

#include "disk_.h"        /* defines BLOCKHEADER */

typedef short NODE;
 
typedef struct bhat *BHAT;
struct bhat { BLOCKHEADER h ;
              NODE top, free;
            }
              ;  
 /* h.type should be 'B' = Bblock  */
 /* h.key is the key of the first object in the block */

#include "bs_.h"    /* defines the BSdata Union */ 

#define NODEX sizeof(BSdata)  /* >= sizeof(KEY) I hope */
union char_key {
                 KEY key ;   
                 BSdata x ;
                 char c[NODEX] ;
               }
                 ;
typedef struct bnode {
               NODE d1, d2;
               union char_key ck;
             } 
                *NODEP ;

#define BNODEMAX ((BLOC_SIZE-sizeof(struct bhat))/sizeof(struct bnode))

struct bblock {
  struct bhat hat;
  struct bnode n[BNODEMAX];
}
                        ;
 


 
/***************************************************************************/
 
 void   Binit            (BP p);
 void   Bget             (BP *bp,NODE *nn,KEY key);
 void   Bdump            (KEY key);
 void   Bnode2str        (char *cq, BP bp, NODEP *r);
 void   Bstr2node        (BS bs, BP bp, NODEP *r);
 void   Bprune           (BP bp, NODE nn,BOOL bothsides);
 int    Bbranchlen       (BP bp, NODE nn);
 int    Bfreelen         (BP bp);
 NODEP  Bnewnode         (BP bp);
 BOOL   Bnextkey         (BP p,KEY *key);
 void   Balias           (BP p,KEY key, KEY newKey);
#endif
 


