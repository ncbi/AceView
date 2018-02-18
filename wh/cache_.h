/*  Last edited: Nov 18 18:23 1998 (fw) */

/* $Id: cache_.h,v 1.5 2016/11/16 23:43:59 mieg Exp $ */
 
   /* wh/cache_.h , 
    * Private to cachesubs.c who handles the second cache.
    */

#ifndef DEFINE_CACHE__h
#define DEFINE_CACHE__h

#ifndef DEF_CACHE_HANDLE
#define DEF_CACHE_HANDLE
  typedef void* CACHE_HANDLE ;
#endif

#define DEF_CACHE
#define CACHEMAGIC 987634

typedef struct scache *CACHE ;
struct scache { KEY key;
		int magic ;
		char type ;
		int lock ;
		int hasBeenModified;
		int refCount ;
		CACHE  up, next;
		CACHE_HANDLE x ;
                Associator tagAss ; /* accelerates accessing tag pointing to many keys */
  AC_HANDLE handle, assHandle ;
	      } ;

#define SIZE_OF_CACHE  sizeof(struct scache) 
 
  /* Note that a straight lock of an cache opened read only is
     forbidden because of possible concurrent access */

/******************************************************************/

/* functions private to the cache package, 
   that allow objcachedisp.c access to cache system variables
   within this file */

CACHE cacheGetKnownTop (void);
CACHE cacheGetLockedTop (void);

#endif
/**************************** eof ************************************/
