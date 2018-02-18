/*  File: banner.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: public header for banner.c
 * Exported functions:
 *              bannerMainStrings 
 *              bannerDataInfoStrings
 *              bannerWrite
 * HISTORY:
 * Last edited: Dec  4 14:51 1998 (fw)
 * Created: Thu Nov 19 15:38:21 1998 (fw)
 *-------------------------------------------------------------------
 */


#ifndef _BANNER_H
#define _BANNER_H

#include "regular.h"	   /* libutil header for Array and BOOL types */

Array bannerMainStrings (char *program,
			 BOOL isGraphical,
			 BOOL isAcembly) ;
Array bannerDataInfoStrings (void) ;
void bannerWrite (Array strings) ;

#endif /* _BANNER_H */
