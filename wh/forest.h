/*  File: forest.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: public header for the forest display module
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 11 10:14 1998 (fw)
 * Created: Fri Dec 11 10:12:25 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: forest.h,v 1.1.1.1 2002/07/19 20:23:20 sienkiew Exp $ */

#ifndef _FOREST_H
#define _FOREST_H

BOOL forestDisplayKeySet (char *title, KEYSET fSet, BOOL isOldGraph);

BOOL forestDisplay (KEY key, KEY from, BOOL isOldGraph);

#endif /* _FOREST_H */
