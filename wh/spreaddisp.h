/*  File: spreaddisp.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: public header for 
 *              tablemaker (spreadsheet) display functions
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 17:11 1998 (fw)
 * Created: Tue Dec 22 16:25:35 1998 (fw)
 *-------------------------------------------------------------------
 */


#ifndef _SPREADDISP_H
#define _SPREADDISP_H

#include "acedb.h"

#include "spread.h"

void spreadDispCreate (BOOL oldGraph);
void spreadTableDisplay (const char *tableName, const char *parms) ;

#endif /* _SPREADDISP_H */
