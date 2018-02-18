/*  Created feb 12, 1997 mieg */

/* $Id: client.h,v 1.1.1.1 2002/07/19 20:23:18 sienkiew Exp $ */

#ifndef DEFINE_CLIENT_H
#define DEFINE_CLIENT_H

#include "table.h"
/***************************************************************/
/*      communications used by the xace and tace clients       */
/***************************************************************/
/* These func pointers are allocated in session.c
 * and initialised in xclient.c
 */

/* 0: close, 1:open, -1:ask state of connection to server */
extern BOOL   (*externalServerState) (int nn) ;
/* Triggers a query on the server side, parses the returned ace file */
extern KEYSET (*externalServer) (KEY key, char *query, char *grep, BOOL getNeighbours) ;
/* A set of calls to export an ace diff of any edition
 * could be used more generally to produce a journal
 */
extern void (*externalSaver) (KEY key, int action) ;
extern void (*externalAlias) (KEY key, char *oldname, BOOL keepOldName) ;
/* sends the journal of ace diffs back to the server */
extern void (*externalCommit) (void) ;
extern KEYSET oldSetForServer ; /* declared in queryexe.c */
extern TABLE* (*externalTableMaker) (char *quer) ;
extern char*  (*externalServerName) (void) ; 

#endif
/****************************************************************/


