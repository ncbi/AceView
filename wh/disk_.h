/*  Last edited: Apr 25 15:14 1995 (srk) */

/* $Id: disk_.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

 
             /*********************************************/
             /* disk_.h                                   */
             /* disk block universal structure            */
             /* replicated in all classes                 */
             /*********************************************/
 
#ifndef DEF_DISK__H 
#define DEF_DISK__H

#ifndef DEF_DISK
#include "disk.h"
#endif

#define BLOC_SIZE 1024      /*number of char in a block*/
typedef struct  blockheader
  { DISK disk,                 /* where */
         nextdisk ;
    int session ;              /* when */
    KEY key ;                  /* who */
    char type ;                /* how */
  } BLOCKHEADER ;


#endif


