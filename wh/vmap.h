/*  File: vmap.h
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: header file for public vMap operations
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  8 13:17 1999 (fw)
 * Created: Thu Nov 26 05:15:10 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: vmap.h,v 1.1.1.1 2002/07/19 20:23:20 sienkiew Exp $ */

#ifndef _VMAP_H
#define _VMAP_H

#include "acedb.h"

/* public DisplayFunc for display type "VMAP" */
BOOL vMapDisplay (KEY key, KEY from, BOOL isOldGraph);


/* public function to recalculate all maps */
void vMapMakeAll(void);

#endif /* _VMAP_H */

/********** end of file **********/
 
