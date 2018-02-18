/*  File: update.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: module to add official data updates to the database
 * Exported functions:
 *              updateData
 *              updateDoAction
 * HISTORY:
 * Last edited: Dec  9 14:59 1998 (fw)
 * Created: Wed Dec  9 14:50:47 1998 (fw)
 *-------------------------------------------------------------------
 */

#ifndef _UPDATE_H
#define _UPDATE_H

#include "acedb.h"

void updateData(void);		  /* open update control window */

void updateDoAction (BOOL doAll); /* do the update, 
				     either one-by-one (doAll == FALSE)
				     or all of them (doAll == TRUE) */


#endif /* _UPDATE_H */
