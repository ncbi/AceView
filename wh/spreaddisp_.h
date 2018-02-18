/*  File: spreaddisp.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: private header for modules within the package
 *              that deals with table maker spreadsheet displays
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 17:33 1998 (fw)
 * Created: Tue Dec 22 16:25:35 1998 (fw)
 *-------------------------------------------------------------------
 */


#ifndef _SPREADDISP__H
#define _SPREADDISP__H

#include "spread_.h"		/* include private spread-ops header,
				   this spreaddisp package extends the
				   spread ops package. */
#include "spreaddisp.h"		/* include public header for spread disp
				   this header completes the interface
				   for spread display ops for members
				   of the spreadDisp-package*/

extern magic_t GRAPH2SPREAD_ASSOC;

/* sprddata.c */
SPREAD currentSpread (char *caller);
void spreadImportKeySet (void);
void spreadSwitchColonnes (void);
void spreadDisplayData(SPREAD spread);
void spreadSelectFromMap(SPREAD spread, int line, int colonne);

/* sprdctrl.c */
void spreadDefineColonne (BOOL force) ;
void spreadShow (void);
void spreadInitColonne (SPREAD spread);

/* sprdmap.c */
void spreadMapCreate (void);


#endif /* _SPREADDISP__H */
