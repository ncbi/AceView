
/*  File: win32util.cpp
 *  Author:	Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *          using graphxt.c code developed by Richard Durbin.
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 *	General utility (generally C encapsulated) functions under WIN32
 * Exported functions:
 *  dbgPos() - used for debug tracing: file/line position reporting
 *	WinTrace() - a TRACE function
 *  AceASSERT() - an ASSERT function
 *	NoMemoryTracking() -  disables the post-mortem MFC memory dump
 * HISTORY:
 * Last edited: Jan 12 11:15 1997 (rbrusk)
 * Created: Jan 11 12:20 1997 (rbrusk)
 *-------------------------------------------------------------------
 */
/* $Id: win32util.cpp,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */

#include "stdafx.h"
#include "win32.h"

#if defined(_DEBUG)

#if 1 // use simpler dbgPos for now
extern "C" const char *dbgPos( const char *caller, int lineno, const char *called )
{
	return called ;
}
#else // !1 == VERBOSE dbgPos
// Will be automatically destroyed at the end of the program?
static CStringList dbgHeap ;

// This procedure creates a debug "position" name composed of a "caller"
// function name and line no, prefixed to the name of the "called" function
extern "C" const char *dbgPos( const char *caller, int lineno, const char *called )
{
	char lineBuf[20] ; int pos ;
	CString a = caller, n = itoa(lineno, lineBuf, 10 ), c = called ;
	// just get the filename out of a pathname
	if((pos = c.ReverseFind('\\')) >= 0) 
		c = c.Mid(pos+1) ;
	dbgHeap.AddHead( (const char*)(a + "(" + n + ")|" + c) ) ;
	return ((const char *)(dbgHeap.GetHead()) ) ;
}
#endif // #if 1

/////////////////////////////////////////////////////////////////////////
// WIN32 ACEDB Trace Function
extern "C" void WinTrace(char *prompt, int id)
{
	TRACE(prompt,id) ;
}

extern "C" void AceASSERT(int condition)
{
	ASSERT(condition) ;
}

// To prevent the post-mortem MFC memory dump!
extern "C" void NoMemoryTracking()
{
	AfxEnableMemoryTracking(FALSE) ;
}

#endif // defined(_DEBUG)

 
