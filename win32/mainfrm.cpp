/*  File: mainfrm.cpp
 *  Author: Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *          adapted from code written by Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: implementation of the CMainFrame class
 * HISTORY:
 * Last edited:May 30 18:20 1996 (rbrusk): WIN32 to 4.3
 *		-	removed intFlag and OnInterrupt() ;
 *			Use WIN32 "Alt+Ctrl+Del" task interrupt for crises interruption!
 * * May 20 13:05 1996 (rbrusk): changed #include's
 * Created: Jul 4 12:20 1995 (rbrusk): ACEDB port to Win32
 *-------------------------------------------------------------------
 */

// $Id: mainfrm.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "mainfrm.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CMainFrame

IMPLEMENT_DYNAMIC(CMainFrame, CMDIFrameWnd)

BEGIN_MESSAGE_MAP(CMainFrame, CMDIFrameWnd)
	//{{AFX_MSG_MAP(CMainFrame)
	ON_WM_CREATE()
	ON_COMMAND(ID_VIEW_MAINWINDOW, OnViewMainWindow)
	ON_WM_PALETTECHANGED()
	ON_WM_QUERYNEWPALETTE()
	ON_WM_SYSCOMMAND()
	ON_COMMAND(ID_HELP, OnHelp)
	ON_COMMAND(ID_GRAPH_PREFERENCES, OnGraphPreferences)
	//}}AFX_MSG_MAP
	// Global help commands
	ON_COMMAND(ID_HELP_INDEX, CMDIFrameWnd::OnHelpIndex)
	ON_COMMAND(ID_HELP_USING, CMDIFrameWnd::OnHelpUsing)
	ON_COMMAND(ID_HELP, CMDIFrameWnd::OnHelp)
	ON_COMMAND(ID_CONTEXT_HELP, CMDIFrameWnd::OnContextHelp)
	ON_COMMAND(ID_DEFAULT_HELP, CMDIFrameWnd::OnHelpIndex)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// arrays of IDs used to initialize control bars
	
// toolbar buttons - IDs are command buttons
static UINT BASED_CODE buttons[] =
{
	// same order as in the bitmap 'toolbar.bmp'
	ID_FILE_NEW,
	ID_FILE_OPEN,
	ID_FILE_SAVE,
		ID_SEPARATOR,
	ID_EDIT_CUT,
	ID_EDIT_COPY,
	ID_EDIT_PASTE,
		ID_SEPARATOR,
	ID_FILE_PRINT,
	ID_APP_ABOUT,
	ID_CONTEXT_HELP,
};

static UINT BASED_CODE indicators[] =
{
	ID_SEPARATOR,           // status line indicator
	ID_INDICATOR_CAPS,
	ID_INDICATOR_NUM,
	ID_INDICATOR_SCRL,
};

/////////////////////////////////////////////////////////////////////////////
// CMainFrame construction/destruction

CMainFrame::CMainFrame()
{
	m_pACEDB_Graph = 0 ; // Main Window Graph hook 
}

CMainFrame::~CMainFrame()
{
}

int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CMDIFrameWnd::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	if (!m_wndToolBar.Create(this) ||
		!m_wndToolBar.LoadBitmap(IDR_MAINFRAME) ||
		!m_wndToolBar.SetButtons(buttons,
		  sizeof(buttons)/sizeof(UINT)))
	{
		TRACE0("Failed to create toolbar\n");
		return -1;      // fail to create
	}

	if (!m_wndStatusBar.Create(this) ||
		!m_wndStatusBar.SetIndicators(indicators,
		  sizeof(indicators)/sizeof(UINT)))
	{
		TRACE0("Failed to create status bar\n");
		return -1;      // fail to create
	}

	/*if( !m_mainDialogBar.Create(this, IDD_MAIN_DIALOG_BAR, CBRS_TOP, IDC_MAIN_ACE_DIALOG) )
	{
		TRACE0("Failed to create ACEDB Main Window dialog control bar\n");
		return -1;      // fail to create
	}
	m_mainDialogBar.EnableDocking(CBRS_ALIGN_ANY);
	DockControlBar(&m_mainDialogBar); */
	EnableDocking(CBRS_ALIGN_ANY);

	// TODO: Delete these three lines if you don't want the toolbar to
	//  be dockable
	m_wndToolBar.EnableDocking(CBRS_ALIGN_ANY);
	DockControlBar(&m_wndToolBar);

	// TODO: Remove this if you don't want tool tips
	m_wndToolBar.SetBarStyle(m_wndToolBar.GetBarStyle() |
		CBRS_TOOLTIPS | CBRS_FLYBY);

	return 0;
}

/////////////////////////////////////////////////////////////////////////////
// CMainFrame diagnostics

#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
	CMDIFrameWnd::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
	CMDIFrameWnd::Dump(dc);
}

#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CMainFrame message handlers

void CMainFrame::OnViewMainWindow() 
{
/*  This has to be implemented differently now?

	static BOOL MainAceWindowIsOpen = FALSE ;

/*
	if( !MainAceWindowIsOpen )
	{
	   // Set the check mark in main menu display
	   // Then unhide the main ace window?
	   m_MainAceWindow->ShowWindow(SW_SHOW) ;
	   m_MainAceWindow->BringWindowToTop() ;
	   MainAceWindowIsOpen = TRUE ;	
	}
	else
	{
	   // Reset the check mark in main menu display
	   // Hide Main Window dialog
	   m_MainAceWindow->ShowWindow(SW_HIDE) ;	
	   MainAceWindowIsOpen = FALSE ;
	}
	*/
}

BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs) 
{
	cs.x	=	0 ;	
	cs.y	=	0 ;
	cs.cx	=	CGraphView::gScreenX ;
	cs.cy	=	CGraphView::gScreenY ;  
	
	return CMDIFrameWnd::PreCreateWindow(cs);
}

extern CPalette gPalette ;

void CMainFrame::OnPaletteChanged(CWnd* pFocusWnd) 
{
	// IF "this" Window initiated the palette change OR
	// the graphics, hence gPalette, are not yet initialized, THEN return.
	if(pFocusWnd == this || !isGraphics ) return ;

	// else...
	OnQueryNewPalette() ;

	// CMDIFrameWnd::OnPaletteChanged(pFocusWnd); // base class call not needed? 
}

BOOL CMainFrame::OnQueryNewPalette() 
{	
	// graphics, hence gPalette not yet initialized?
	if( !isGraphics ) return TRUE ;  // it's logical palette == System palette
    
    CDC *dc = GetDC();
	ASSERT_VALID(dc) ;

    CPalette *hOldPal = dc->SelectPalette( &gPalette, FALSE);
	ASSERT_VALID( hOldPal ) ;

    UINT i = dc->RealizePalette( );
    if (i != 0) {
        InvalidateRect( NULL, TRUE); 
    }

    dc->SelectPalette( hOldPal, TRUE);
    dc->RealizePalette( );
    ReleaseDC( dc );
	return i ;
	// return CMDIFrameWnd::OnQueryNewPalette(); // base class call not needed?
}

extern "C" void wormClose() ; /* normal exit */

void CMainFrame::OnSysCommand(UINT nID, LPARAM lParam) 
{
	UINT sysCode = nID & 0xFFF0 ;
	if( sysCode == SC_CLOSE )
		wormClose() ;
	else
		CMDIFrameWnd::OnSysCommand(nID, lParam);
}

extern "C" void helpOn(char *text) ; 

void CMainFrame::OnHelp() 
{
	helpOn(0) ;	
}

void CMainFrame::OnGraphPreferences() 
{
	// TODO: Add your command handler code here
	
}
 
 
