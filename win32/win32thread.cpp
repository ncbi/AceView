
/*  File: win32thread.cpp
 *  Author: Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 *	WIN32 multithreading functionality for winace et al.
 *
 * Exported functions:
 *	hACETHREAD initParseThread(FILE *fil)
 *	void getParseThreadLine(hACETHREAD ptHandle, int isSkip, KEYSET ks) ;
 *	void closeParseThread(hACETHREAD handle) ;
 *
 * HISTORY:
 * Last edited:Mar 30 14:30 1997 (rbrusk)
 * Created: Mar 30 14:30 1997 (rbrusk)
 *-------------------------------------------------------------------
 */

/* $Id: win32thread.cpp,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */

/* WIN32 multithreading support */

#include "stdafx.h"

#if 0

extern "C"
{
#include "regular.h"
#include "keyset.h"
#include "session.h"
}

#include "parsedata.h"

#define WIN32THREAD_CPP
#include "win32thread.h"


IMPLEMENT_DYNAMIC(CParseData, CObject)

CParseData::CParseData(FILE *pFile)
{
	m_pParseThread = 0 ;
	FILE *m_pFile = pFile ;
	int m_isSkip = FALSE ;
	KEYSET m_ks = 0 ;
}

extern void parseLine (int isSkip, KEYSET ks) ;

UINT AFX_CDECL AceParseThreadProc(LPVOID pParam)
{
	ASSERT( pParam != 0 ) ; 
	CParseData *pParseData = (CParseData *)pParam ;
	while ( TRUE /* waiting for data? */)
	{
		/* Lock the pParseData object?*/
		parseLine(pParseData->m_isSkip, pParseData->m_ks) ;
		/* Release the pParseData object? */
	}
	return 0 ;
}

extern "C" hACETHREAD initParseThread(FILE *pFile)
{
	ASSERT( pFile != 0 ) ;
	CParseData *pParseData = new CParseData( pFile ) ;
	pParseData->m_pParseThread = AfxBeginThread( &AceParseThreadProc, (LPVOID)pParseData ) ;
	return ((hACETHREAD)pParseData );
}

extern "C" void getParseThreadLine(hACETHREAD ptHandle, int isSkip, KEYSET ks)
{
	ASSERT( ptHandle != 0 ) ;
	CParseData *pParseData = (CParseData *)ptHandle ;

	if( TRUE /* pParseData object is not locked */ )
	{
		/* Lock the pParseData object?*/
		pParseData->m_isSkip = isSkip ;
		pParseData->m_ks = ks;
		/* Release the pParseData object? */
	}
}

extern "C" void closeParseThread(hACETHREAD ptHandle)
{
	ASSERT( ptHandle != 0 ) ;
	CParseData *pParseData = (CParseData *)ptHandle ;
	delete pParseData->m_pParseThread ;
	delete pParseData ;
}

#include "graph.h"

extern "C" void parseGraphBoxDraw (int k, int fcol, int bcol)
{
	/*	Need to convert this to some form of
		threadsafe communication? */
	graphBoxDraw (k, fcol, bcol) ;
}

static CWinThread *pAceThread = 0 ;
struct CAceThreadData
{
	int *argcp ; 
	char **argv ;
} ;

UINT AFX_CDECL AceThreadProc(LPVOID pParam)
{
	ASSERT( pParam != 0 ) ; 
	CAceThreadData *pAceThreadData = (CAceThreadData *)pParam ;

	sessionInit (pAceThreadData->argcp, pAceThreadData->argv) ;

	/* What do I do here? AskOrders()? */

	return 0 ;
}

extern "C" void winAceInit(int *argcp, char **argv)
{
	CAceThreadData *pAceThreadData = new  CAceThreadData;
	pAceThreadData->argcp = argcp ;
	pAceThreadData->argv = argv ;

	pAceThread = AfxBeginThread( &AceThreadProc, (LPVOID)pAceThreadData ) ;

	/* Wait for completion of sessionInit() before proceeding? */

}
#else

typedef void *hACETHREAD ;
typedef void *KEYSET ;

/*
extern "C" hACETHREAD initParseThread(FILE *pFile)
{
	// stub 
	return ((hACETHREAD)0 );
}
*/
extern "C" void getParseThreadLine(hACETHREAD ptHandle, int isSkip, KEYSET ks)
{
	/* stub */
}

extern "C" void closeParseThread(hACETHREAD ptHandle)
{
	/* stub */
}

extern "C" void parseGraphBoxDraw (int k, int fcol, int bcol)
{
	/* stub */
}

extern "C" void winAceInit(int *argcp, char **argv)
{
	/* stub */
}
#endif

 
