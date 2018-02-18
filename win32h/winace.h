/*  WinAce.h : main header file for the WINACE application
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * HISTORY:
 * Last edited: Jun 10 10:15 1996 (rbrusk):
 *		-	added m_Cursor to CWinAceApp (see also CMainFrame and CGraphView)
 * Created: Jul 21 18:30 1995 (rbrusk)
 *-------------------------------------------------------------------
 */
/* $Id: winace.h,v 1.1.1.1 2002/07/19 20:23:27 sienkiew Exp $ */

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

// To trick the resource.h file to define _APS_NEXT_COMMAND_VALUE
#define APSTUDIO_INVOKED
#include "resource.h"

extern "C"
{
#include "graph_.h"
}

#include "graphloop.h"

/////////////////////////////////////////////////////////////////////////////
// CWinAceApp:
// See WinAce.cpp for the implementation of this class
//

class CWinAceApp : public CWinApp
{

private:

	// reentrant graphLoop() management:
	// This ObList functions as a global stack variable which
	// keeps track of blocking and non-blocking calls
	// to graphLoop() and graphLoopReturn()
	// Can't specify CGraphLoop here because of #include file order? 
	static CTypedPtrList<CObList,CGraphLoop*> m_graphLoops ;

public:
	CWinAceApp();					
	virtual ~CWinAceApp() ;

	HCURSOR m_Cursor ; // Current active application mouse cursor

// Attributes

	// IsGraphLoop() returns TRUE if graphLoop() is executing, ie. non empty stack
	static BOOL IsGraphLoop() { return !CWinAceApp::m_graphLoops.IsEmpty() ; }

	// CurrentGraphLoop() returns the current CGraphLoop  context object
	// which is executing without removing it
	static CGraphLoop *CurrentGraphLoop()
		{ return CWinAceApp::m_graphLoops.GetHead() ; }

	// pushGraphLoop(  ) and popGraphLoop( )
	// manage the CGraphLoop context object stack
	static BOOL pushGraphLoop( CGraphLoop *pGL ) ;
	static BOOL popGraphLoop( CGraphLoop *pGL ) ; 
	static void CWinAceApp::SplashBox() ;

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CWinAceApp)
	public:
	virtual BOOL InitInstance();
	virtual int  ExitInstance();
	virtual void AddToRecentFileList(LPCTSTR lpszPathName);
	virtual int  Run();
	virtual BOOL PumpMessage();
	virtual BOOL OnIdle(LONG lCount);
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CWinAceApp)
	afx_msg void OnAppAbout();
	afx_msg void OnUpdateRecentFileMenu(CCmdUI* pCmdUI);
	afx_msg	BOOL OnOpenRecentFile(UINT nID) ;
	afx_msg void OnAppExit();
	afx_msg void OnGraphPreferences();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

extern CWinAceApp theApp ;

/////////////////////////////////////////////////////////////////////////////
 
 
 
