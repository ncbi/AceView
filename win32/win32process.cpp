/*  File: win32process.cpp
 *  Author: Richard M. Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) Richard M. Bruskiewich (1996)
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Win32 version of external process management & communication
 * Exported functions: callMosaic. 
 *
 * HISTORY:
 * Last edited: Jun 5 14:00 1996 (rbrusk): WIN32 to 4.3 
 * Created: Jan  1 14:33 1996 (rmb)
 *-------------------------------------------------------------------
 */

// $Id: win32process.cpp,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $

#include <stdafx.h>

extern "C"
{
#include "regular.h"
}

// process.h: This is the Microsoft, not the Kocab & Senger version!
// Must ensure that the K & S one is not on any include path during
// WIN32 compilation of this module!
#include <process.h>

static int mosaicPID = -1;

extern "C" BOOL callMosaic(char *URL)
{
	if(mosaicPID > 0) // is the mosaic process already open?
	{
		messout("Type ALT+Tab to go to open Web Browser to access ACEDB HTML help") ;  
		return TRUE ;
	}

	char *mosaic = getenv("WEB_BROWSER") ;
	if(!mosaic || !*mosaic) mosaic = "mosaic" ;

	if(!*URL)
		mosaicPID = _spawnlp(	_P_NOWAIT, mosaic, NULL ) ;
	else
		mosaicPID = _spawnlp(	_P_NOWAIT, mosaic, "-home", URL, NULL ) ;
	
	if(mosaicPID == -1)
		return FALSE ;
	else
		return TRUE ;
}
 
