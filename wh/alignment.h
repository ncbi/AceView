/*  File: alignment.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: public header for alignment.c
 * Exported functions:
 *              alignDumpKeySet
 *              alignDumpKey
 * HISTORY:
 * Last edited: Nov 19 15:01 1998 (fw)
 * Created: Thu Nov 19 14:59:51 1998 (fw)
 *-------------------------------------------------------------------
 */

#ifndef _ALIGNMENT_H
#define _ALIGNMENT_H

int alignDumpKeySet (KEYSET kSet, FILE *fil);
BOOL alignDumpKey (KEY key, FILE *fil);

#endif /* _ALIGNMENT_H */
