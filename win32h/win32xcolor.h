
/*-------------------------------------------------------------------
 *  File:  Win32XColor.h
 *  Author: Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *          adapted from code written by Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *  This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 *	Also contains X Window colormap emulation types, et al.
 *
 *  HISTORY:
 *  Last edited: Jun 7 15:42 1996 (rbrusk): WIN32 to 4.3
 *		-	now only used in graphramp.c
 * * Feb 18 22:22 1996 (rbrusk): WIN32 to 4.2 
 *  Created: Dec 22 03:17 1995 (rbrusk)
 *-------------------------------------------------------------------
 */

/* $Id: win32xcolor.h,v 1.1.1.1 2002/07/19 20:23:27 sienkiew Exp $ */

/* Partial emulation of the X Windows colormap data structures? */

typedef struct
{
	unsigned long 	pixel;  /* Pixel value after allocation 	*/
	unsigned short 	red,	/* Red intensity (0 - 65,535) 		*/
					green,	/* Green intensity (0 - 65,535) 	*/
	  				blue;	/* Blue intensity (0 - 65,535) 		*/
	char			flags;	/* Used when storing colors			*/
	char			pad;	/* Just so structure size is even	*/
} XColor ;

/***** end of file ********/
 
