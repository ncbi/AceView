/*  File: ACEDialogs.c
 *  Author:	Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *          using graphxt.c code developed by Richard Durbin.
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: graphxt.c globals ported to Win32 version of graph package
 * Exported functions: SysBeep(), graphOut(), graphPrompt() & graphQuery() 
 * HISTORY:
 * Created: Feb 24 17:12 1996 (rbrusk): WinAce port to Visual C++ Rel. 4.0
 *		-	Procedural indirection introduced for practical reasons,
 *			to the  graph*() functions exported above.  win32Graph*()
 *			functions found in \acedb\win32\windialogs.cpp
 *		-	isGraphics global defined here for practical reasons...
 *-------------------------------------------------------------------
 */

/* $Id: acedialogs.c,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $ */

#include "regular.h"

extern int isGraphics ;

extern void win32SysBeep(int n) ;
void SysBeep(int n) { win32SysBeep(n) ; }

extern void win32GraphOut(char *message) ;
void graphOut(char *message) { win32GraphOut( message ) ; }

extern char *win32GraphPrompt(char *prompt, char *dfault ) ;
BOOL graphPrompt(char *prompt, char *dfault, char *fmt)
{
	static char *answer ;

	while(TRUE)
	{
		answer = win32GraphPrompt( prompt, dfault ) ;

		if ( answer != NULL )
		{
			freeforcecard ( answer ) ;
			/* Verify result */
			if ( freecheck ( fmt ) )
				return TRUE ;
			else	/* Complain then repeat prompt */
				win32GraphOut("Sorry, invalid response. Try again or cancel" ) ;
		}
		else
			return FALSE ;/* Cancel! */
	}
}

extern BOOL win32GraphQuery (char *text) ;
BOOL graphQuery (char *text) { return win32GraphQuery ( text ) ; }
	
 
