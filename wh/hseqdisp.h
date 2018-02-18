/*  File: hseqdisp.h
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: header file for public hseqdisp operations
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  8 13:17 1999 (fw)
 * Created: Thu Nov 26 05:15:10 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: hseqdisp.h,v 1.4 2007/07/03 17:00:58 mieg Exp $ */

#ifndef _HSEQDISP_H
#define _HSEQDISP_H

#include "acedb.h"

/* public DisplayFunc for display type "HSEQDISP" */
BOOL hseqDisplay (KEY key, KEY from, BOOL isOldGraph) ;
BOOL glocDisplay (KEY key, KEY from, BOOL isOldGraph) ;
BOOL glocBigDisplay (KEY key, KEY from, BOOL isOldGraph) ;
BOOL htileDisplay (KEY key, KEY from, BOOL isOldGraph) ;

#endif /* _HSEQDISP_H */

/********** end of file **********/
 
