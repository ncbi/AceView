/*  Last edited: Apr 25 15:56 1995 (srk) */


 
             /*********************************************/

/* $Id: disk__.h,v 1.2 2010/06/15 17:54:49 mieg Exp $ */
             /* disk__.h                                  */
             /* diskblock generic structure               */
             /*********************************************/
 
/***************************************************************************/

#ifndef DEF_DISK___h
#define DEF_DISK___h 

#ifndef DEF_BP
#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct block *BLOCKP ;
typedef BLOCKP BP ;
#endif

#include "disk_.h"   /* defines BLOCKHEADER */


#define BLKMX (BLOC_SIZE - sizeof(BLOCKHEADER))
typedef struct block
  {
    BLOCKHEADER  h ;
    char c[BLKMX];
  }
    BLOCK;   /*the transfer unit between disk and cache*/



#endif

