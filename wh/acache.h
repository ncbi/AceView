
/* $Id: acache.h,v 1.2 2015/08/17 22:53:54 mieg Exp $ */

/* prototypes for the adisk routines
 *  new way to store data in test for acedb5
 */


#ifndef DEFINE_ACACHE_H
#define DEFINE_ACACHE_H

void aCacheKill (KEY key) ;
BOOL aCacheGet (KEY key, KEY *pp, Array *ap, char *format) ;
void aCacheStore (KEY key, KEY parent, Array a, char *format) ;
void aCacheStoreDestroy (KEY key, KEY parent, Array a, char *format) ;

#endif
