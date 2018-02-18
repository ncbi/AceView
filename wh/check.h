/*  File: check.h
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: May 21 16:18 1994 (rd)
 * Created: Fri Aug 20 12:21:54 1993 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: check.h,v 1.1.1.1 2002/07/19 20:23:18 sienkiew Exp $ */

#ifndef DEF_CHECK_H
#define DEF_CHECK_H

BOOL checkRegExpSyntax (char *cp) ;
BOOL checkText (char *text, KEY cst) ;
BOOL checkKey (KEY new, KEY cst) ;
BOOL checkData (void *xp, KEY type, KEY cst) ;

#endif
