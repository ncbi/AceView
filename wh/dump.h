/*  Last edited: Aug 19 18:30 1998 (rd) */

/* $Id: dump.h,v 1.2 2007/10/29 20:36:32 mieg Exp $ */

#ifndef DEF_DUMP_H
#define DEF_DUMP_H

#include "keyset.h"
#include "query.h"
#include "bs.h"

extern BOOL dumpTimeStamps;
extern BOOL dumpComments;

void dumpAll (void) ;
void dumpAllNonInteractive (char *directory, BOOL split) ;

BOOL dumpClassPossible (int classe, char style) ; /* 'a' for ace style */
BOOL dumpKey (KEY key, FILE* ff, Stack ss) ;

    /* following dump on freeOut */
BOOL dumpKeyOut (KEY key) ; 
BOOL dumpKeyBeautifully (KEY key, char style, COND cond)  ;
void niceDump (OBJ objet, char style) ;

BOOL dnaDumpKeyCstyle (KEY key) ;
BOOL peptideDumpKeyCstyle (KEY key) ;

BOOL aceBinaryDump (FILE* f, Stack buf, KEY k)  ;
BOOL aceBinaryParse (int level, KEY key) ;
void aceBinaryDoDump (unsigned const char *ccp, int nn) ;


#endif
 
