/*  File: spread.h
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: public header for spreadsheet operations
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 17:21 1998 (fw)
 * Created: Thu Aug  6 13:40:24 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: spread.h,v 1.7 2007/03/21 21:03:14 mieg Exp $ */
  
#ifndef DEF_SPREAD_H  
#define DEF_SPREAD_H

#include "acedb.h"

#include "table.h"

/************************************************************/

typedef struct SpreadStruct *SPREAD; /* opaque forward declaration */

/************************************************************/

/* sprdop.c */
SPREAD spreadCreate (void) ;
void uSpreadDestroy (SPREAD spread) ;
#define spreadDestroy(spread) (uSpreadDestroy(spread), (spread) = 0)

/* clear the DNA and other long arrays already exported */
void spreadCleanUp (SPREAD spread) ;

TABLE *spreadToTable (SPREAD spread, AC_HANDLE h) ;
KEYSET spreadFilterKeySet (SPREAD spread, KEYSET ks);
KEYSET spreadGetKeyset (SPREAD spread) ;
KEYSET spreadGetPrecalculatedKeySet (SPREAD spread, KEY key, char *cr);
KEY spreadRecomputeKeySet (SPREAD spread, KEYSET ks, int last, int minx, int maxx) ;
int  spreadDoDump (SPREAD spread, int sep0, char style, BOOL beginTable) ;
BOOL spreadDumpLegend (SPREAD spread, char style) ;

/* sprddef.c */
BOOL spreadDoSaveInObj (SPREAD spread, KEY tableKey) ;
BOOL spreadDoReadDefinitions (SPREAD spread, KEY key, FILE *f,
			      Stack s, const char *parms, BOOL noParms) ;

#endif /* DEF_SPREAD_H */
 
