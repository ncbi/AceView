/*  Last edited: Jul 12 12:26 1996 (srk) */

/* $Id: cache.h,v 1.4 2007/03/21 21:03:14 mieg Exp $ */
 
   /* wh/cache.h , 
    * Prototypes of public functions of cache package 
    */

#ifndef DEFINE_CACHE_h
#define DEFINE_CACHE_h

#ifndef DEF_CACHE
#define DEF_CACHE
  typedef void* CACHE ;
#endif

#ifndef DEF_CACHE_HANDLE
  typedef void* CACHE_HANDLE ;
#endif

  /* Note that a straight lock of an cache opened read only is
     forbidden because of possible concurrent access */

int cacheStatus (int *used, int *locked, int*known, int *modified, 
		 int *nc2rp, int *nc2wp, int *nccrp, int *nccwp) ;
CACHE cacheCreate(KEY k, CACHE_HANDLE *bsp, BOOL newCopy) ;
CACHE cacheUpdate(KEY k, CACHE_HANDLE *bsp) ;
void cacheClone(CACHE source, CACHE dest, CACHE_HANDLE *bsp) ;
void cacheDestroy(CACHE v) ;
void cacheSave(CACHE v) ;
void cacheCheck(KEY key, CACHE v) ;
AC_HANDLE cacheStoreHandle (CACHE v) ; /* make it accessible where needed */
void cacheList (int type) ;

#endif
/******************************************************************/
 
