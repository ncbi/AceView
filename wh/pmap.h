/*  File: pmap.h
 *  Author: Neil Laister (nl1@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 11 15:03 1998 (fw)
 * Created: Fri Jan 17 15:43:30 1992 (nl1)
 *-------------------------------------------------------------------
 */

/* $Id: pmap.h,v 1.1.1.1 2002/07/19 20:23:20 sienkiew Exp $ */

#ifndef _PMAP_H
#define _PMAP_H

#include "acedb.h"

void pMapSetHighlightKeySet (KEYSET k); /*argument NULL causes reset*/
BOOL pMapDisplay (KEY key, KEY from, BOOL isOldGraph); /* DisplayFunc */


#endif /* _PMAP_H */
