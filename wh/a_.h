/*  Last edited: Dec  4 14:51 1998 (fw) */

/* $Id: a_.h,v 1.2 2010/06/15 17:54:48 mieg Exp $ */

 
 
                    /* a_.h                            */
                    /* definition of the A block       */
                    /* private to asubs.c              */

 
/***************************************************************************/

#ifndef DEFINE_A__h
#define DEFINE_A__h
 
#ifndef DEF_BP
#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct ablock *ABLOCK ;
typedef ABLOCK BP ;
#endif
 
#include "disk_.h"   /* defines BLOCKHEADER */
#include "regular.h"   /* for Array definition */

typedef struct ahat *AHAT;
struct ahat { BLOCKHEADER h ;
	      int size, max, localmax ;
#ifdef ACEDB4
	      KEY parent ;	/* B object that "owns" me */
#endif
            }
              ;  
 /* h.type should be 'A' = permanent Array  */
 /* h.key is the key of the only object in the block */
 /* size, max matches the array */
 /* localmax is the number of elements stored in this block */

#define AMAX      (BLOC_SIZE-sizeof(struct ahat))
struct ablock {
  struct ahat hat;
  char  n[AMAX];
}
                        ;

typedef Array OBJ_HANDLE ;
#define DEF_OBJ_HANDLE

#endif

/*****/
/*****/



