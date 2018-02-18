
/* WinAce.cpp : implementation file
 *  Author: Richard Bruskiewich	(rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Defines the class behaviors for the WIN32 application.  
 * Exported:  class CGraphLoop
 * HISTORY:
 * Last edited: Jan 14 14:56 1996 (rmb) 
 * Created: Jun 25 11:56 1995 (rmb)
 *-------------------------------------------------------------------
 */
/* $Id: winace.cpp,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */

#include "stdafx.h"
#include "win32.h"

extern "C" {
#include "acedb.h"
#include "session.h"
#include "array.h"
}

#include "WinAce.h"
#include "SplashBox.h"

#include "mainfrm.h"

#include "CGraph.h"

#include "AceFileDoc.h"
#include "AceFileView.h"
#include "ACEDBProfile.h"
#include "ACEDBProfileView.h"
#include "PlainView.h"
#include "FullScrollView.h"
#include "TextScrollView.h"
#include "Fitview.h"
#include "TextFitView.h"
#include "MapScrollView.h"
#include "PixelScrollView.h"
#include "PixelFitView.h"
// #include "FullEditView.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

#define new DEBUG_NEW

/////////////////////////////////////////////////////////////////////////////
// CWinAceApp

BEGIN_MESSAGE_MAP(CWinAceApp, CWinApp)
	//{{AFX_MSG_MAP(CWinAceApp)
	ON_COMMAND(ID_APP_ABOUT, OnAppAbout)
	ON_UPDATE_COMMAND_UI(ID_FILE_MRU_FILE1, OnUpdateRecentFileMenu)
	ON_COMMAND_EX_RANGE(ID_FILE_MRU_FILE1, ID_FILE_MRU_FILE16, OnOpenRecentFile)
	ON_COMMAND(ID_APP_EXIT, OnAppExit)
	ON_COMMAND(ID_GRAPH_PREFERENCES, OnGraphPreferences)
	//}}AFX_MSG_MAP
	// Standard file based document commands
	ON_COMMAND(ID_FILE_NEW, CWinApp::OnFileNew)
	ON_COMMAND(ID_FILE_OPEN, CWinApp::OnFileOpen)
	// Standard print setup command
	ON_COMMAND(ID_FILE_PRINT_SETUP, CWinApp::OnFilePrintSetup)
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CWinAceApp construction

CWinAceApp::CWinAceApp()
{
	// Place all significant initialization in InitInstance
}

CWinAceApp::~CWinAceApp()
{
}

/////////////////////////////////////////////////////////////////////////////
// A pointer to the one and only CWinAceApp object

CWinAceApp theApp ;

/////////////////////////////////////////////////////////////////////////////
// CWinAceApp initialization

extern char **CommandLineToArgv(LPSTR lpCmdLine, int * pArgc) ;
extern void CleanUpApp() ;
static int argc ;
static char **argv ;

void openACEDB() ;
BOOL closeACEDB() ;

BOOL CWinAceApp::InitInstance()
{
	// Standard initialization
	// If you are not using these features and wish to reduce the size
	//  of your final executable, you should remove from the following
	//  the specific initialization routines you do not need.

	Enable3dControls();

	SetRegistryKey( "ACEDB4" ) ; // Use the system registry to store profile data

	LoadStdProfileSettings(4);  // Load standard profile options

	CGraphView::GetGraphDevAttributes() ;

	// Register the application's document templates.  Document templates
	//  serve as the connection between documents, frame windows and views.

	CMultiDocTemplate* pDocTemplate;

	// ACEDB profile file type
	pDocTemplate = new CMultiDocTemplate(
		IDR_MAINFRAME,
		RUNTIME_CLASS(CACEDBProfile),
		RUNTIME_CLASS(CGraphWindow),          // standard MDI child frame
		RUNTIME_CLASS(CACEDBProfileView));
	AddDocTemplate(pDocTemplate);

	// Ace files are treated as a separate document type, with a CEdit view
	pDocTemplate = new CMultiDocTemplate(
		IDR_ACE_FILE_TYPE,
		RUNTIME_CLASS(CAceFileDoc),
		RUNTIME_CLASS(CGraphWindow),          // standard MDI child frame
		RUNTIME_CLASS(CAceFileView));
	AddDocTemplate(pDocTemplate);

	// create main MDI Frame window
	CMainFrame* pMainFrame = new CMainFrame;
	if (!pMainFrame->LoadFrame(IDR_MAINFRAME, WS_MAXIMIZE |
	 						   WS_OVERLAPPEDWINDOW |
	 						   WS_VSCROLL | WS_HSCROLL ))
		return FALSE;
	m_pMainWnd = pMainFrame ;
	HICON icon = LoadIcon(IDR_MAINFRAME) ;
	pMainFrame->SetIcon(icon, TRUE) ;

	// Default cursor used unless graphWaitCursor() is called
	m_Cursor = NULL ;

	// Enable DDE Execute open
	EnableShellOpen();
	RegisterShellFileTypes(TRUE /* == WIN95! */);

	// Enable drag/drop open
	m_pMainWnd->DragAcceptFiles();

	pMainFrame->ShowWindow(m_nCmdShow);
	pMainFrame->UpdateWindow();

	/* Framework encapsulations of ACEDB GraphTypes:
	        PLAIN, TEXT_SCROLL, TEXT_FIT, MAP_SCROLL, 
		    PIXEL_SCROLL, TEXT_FULL_SCROLL, PIXEL_FIT, TEXT_FULL_EDIT
		These are internal WIN32 framework document types, which
		should be invisible to the user from the File Manager, hence
		I add the Doc templates here...
	*/
	pDocTemplate = new CMultiDocTemplate(
		IDR_PLAIN_GRAPH_TYPE,
		RUNTIME_CLASS(CGraph),
		RUNTIME_CLASS(CGraphWindow),          // standard MDI child frame
		RUNTIME_CLASS(CPlainGraphView));
	AddDocTemplate(pDocTemplate);

	pDocTemplate = new CMultiDocTemplate(
		IDR_TEXT_SCROLL_GRAPH_TYPE,
		RUNTIME_CLASS(CGraph),
		RUNTIME_CLASS(CGraphWindow),          // standard MDI child frame
		RUNTIME_CLASS(CTextScrollGraphView));
	AddDocTemplate(pDocTemplate);

	pDocTemplate = new CMultiDocTemplate(
		IDR_TEXT_FIT_GRAPH_TYPE,
		RUNTIME_CLASS(CGraph),
		RUNTIME_CLASS(CGraphWindow),          // standard MDI child frame
		RUNTIME_CLASS(CTextFitGraphView));
	AddDocTemplate(pDocTemplate);

	pDocTemplate = new CMultiDocTemplate(
		IDR_MAP_SCROLL_GRAPH_TYPE,
		RUNTIME_CLASS(CGraph),
		RUNTIME_CLASS(CGraphWindow),          // standard MDI child frame
		RUNTIME_CLASS(CMapScrollGraphView));
	AddDocTemplate(pDocTemplate);

	pDocTemplate = new CMultiDocTemplate(
		IDR_PIXEL_SCROLL_GRAPH_TYPE,
		RUNTIME_CLASS(CGraph),
		RUNTIME_CLASS(CGraphWindow),          // standard MDI child frame
		RUNTIME_CLASS(CPixelScrollGraphView));
	AddDocTemplate(pDocTemplate);

	pDocTemplate = new CMultiDocTemplate(
		IDR_TEXT_FULL_SCROLL_GRAPH_TYPE,
		RUNTIME_CLASS(CGraph),
		RUNTIME_CLASS(CGraphWindow),          // standard MDI child frame
		RUNTIME_CLASS(CTextFullScrollGraphView));
	AddDocTemplate(pDocTemplate);

	pDocTemplate = new CMultiDocTemplate(
		IDR_PIXEL_FIT_GRAPH_TYPE,
		RUNTIME_CLASS(CGraph),
		RUNTIME_CLASS(CGraphWindow),          // standard MDI child frame
		RUNTIME_CLASS(CPixelFitGraphView));
	AddDocTemplate(pDocTemplate);

/**** Not yet supported?
	pDocTemplate = new CMultiDocTemplate(
		IDR_TEXT_FULL_EDIT_GRAPH_TYPE,
		RUNTIME_CLASS(CGraph),
		RUNTIME_CLASS(CGraphWindow),          // standard MDI child frame
		RUNTIME_CLASS(CRichTextEditView));
	AddDocTemplate(pDocTemplate);
*/

	// Display ACEDB Banner
	
	if (!getenv("ACEDB_NO_BANNER"))
		CWinAceApp::SplashBox() ;

	// Initialize ACEDB
	
	argv = CommandLineToArgv( m_lpCmdLine, &argc) ;

	openACEDB() ; // open the current ACEDB database

	return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
private:
	Array m_banner, m_dataInfo ;

public:
	CAboutDlg();
	virtual ~CAboutDlg();

// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	CString	m_aceBanner;
	//}}AFX_DATA

// Implementation
protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//{{AFX_MSG(CAboutDlg)
	afx_msg void OnAcedbBanner();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

// AceIsAlive is modified in openACEDB() and closeACEDB() below
BOOL AceIsAlive = FALSE ;

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	m_aceBanner = _T("");
	//}}AFX_DATA_INIT

	m_banner = bannerMainStrings("WinAce", TRUE ) ;
	int i ;

	if(m_banner)
		for (i = 0 ; i < arrayMax(m_banner) ; ++i)
		{
			m_aceBanner += arr(m_banner, i, char*) ;
			m_aceBanner += "\r\n" ;
		}
	else // m_banner == 0 
		m_aceBanner += "\t\tWarning: ACEDB Banner Not Found?\r\n" ;


	if(AceIsAlive)
	{
		m_dataInfo = bannerDataInfoStrings() ;
		if(m_dataInfo)
		{	m_aceBanner += "\r\nCurrent Database Information:\r\n\r\n" ;
			for (i = 0 ; i < arrayMax(m_dataInfo) ; ++i)
			{
				m_aceBanner += arr(m_dataInfo, i, char*) ;
				m_aceBanner += "\r\n" ;
			}
		}
	} else m_dataInfo = 0 ;
}
CAboutDlg::~CAboutDlg()
{
	arrayDestroy(m_banner) ;
	arrayDestroy(m_dataInfo) ;
}

	void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	DDX_Text(pDX, IDC_ACEDB_BANNER, m_aceBanner);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
	ON_EN_SETFOCUS(IDC_ACEDB_BANNER, OnAcedbBanner)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

// App command to run the dialog
void CWinAceApp::OnAppAbout()
{
    CAboutDlg aboutDlg ;
	aboutDlg.DoModal();
}

// App command to run the dialog
void CWinAceApp::SplashBox()
{
    CSplashBox splashBox ;
	splashBox.DoModal();
}

/////////////////////////////////////////////////////////////////////////////
// CAboutDlg commands

void CAboutDlg::OnAcedbBanner() 
{
	// TODO: Add your control notification handler code here
	UpdateData(FALSE) ; 
}

/////////////////////////////////////////////////////////////////////////////
// CWinAceApp commands

extern "C" void *Win32MainWindow()
{
  return (void *) theApp.m_pMainWnd ;
}

static int ACEDB_ExitCode = 0;
extern void saveSessionParameters() ;

int CWinAceApp::ExitInstance() 
{
	saveSessionParameters() ;

	int WIN32_ExitCode = CWinApp::ExitInstance();
	return( (ACEDB_ExitCode > WIN32_ExitCode) ? ACEDB_ExitCode : WIN32_ExitCode ) ; 
}

/////////////////////////////////////////////////////////////////////////////
// MRU file list overridden to be disabled

void CWinAceApp::AddToRecentFileList(LPCTSTR lpszPathName) 
{
}

void CWinAceApp::OnUpdateRecentFileMenu(CCmdUI* pCmdUI) 
{
}

BOOL CWinAceApp::OnOpenRecentFile(UINT nID)
{
	return FALSE;
}

// Memory leak checking
//#define _DUMP_MEM_STATS  // disabled for now
#if defined(_DEBUG) && defined(_DUMP_MEM_STATS)
CMemoryState before, after, diff ;
#endif

extern "C" void wormInit( int argc, char **argv ) ;
extern "C" void askOrders(void) ;

void openACEDB()
{
// Memory leak checking
#if defined(_DEBUG) && defined(_DUMP_MEM_STATS)
	before.Checkpoint() ;
#endif

    wormInit ( argc, argv ) ;
	askOrders();
	AceIsAlive = TRUE ;
}

/*************************************************/

/* a WIN32 variant of xacemain wormClose() 
   which does NOT call wormClose2() yet... */

extern "C" BOOL DoIQuit (void) ;

extern "C" void wormClose2() ;
int MainGraphID = 1; // Assume that graph #1 is the main one?

BOOL closeACEDB()
{
	
	/* Return if not an error crash
	   and user is NOT ready to exit ACEDB */
	if( !( ACEDB_ExitCode || DoIQuit() ) ) return FALSE;
	
	// Otherwise, continue exiting but tidy up down ACEDB first...
	if( !AceIsAlive )
	{
		AfxMessageBox("UNRECOVERABLE WINACE ERROR!\nmesscrash() invoked during wormClose2()?",
					  MB_OK | MB_ICONSTOP ) ;
		theApp.m_pMainWnd->SendMessage(WM_CLOSE) ;
		AfxThrowUserException() ;
	}
	else
		AceIsAlive = FALSE ;  // Avoids infinite recursion in messcrash()?

	// therefore, close down current ACEDB database
	// ...first, activate the MainGraphWindow
	graphActivate( MainGraphID ) ;
	wormClose2();

	return TRUE ;
}

extern "C" void winExit( int exitcode = 0 ) /* Windows implementation of "exit()" */
{
	ACEDB_ExitCode = exitcode ;
	static BOOL AceClosed = FALSE ;

	// Return if not error condition and User is not ready to close current ACEDB database
	if( !AceClosed )
	   if( !( AceClosed = closeACEDB() ) ) return ;  // Return if user not ready to leave ACEDB

	/*// Return if not error condition and User is not ready to exit WinAce shell
	if( !( 	ACEDB_ExitCode || 
	   		messQuery ("Do you really want to quit the WinAce shell?") ) ) return ;	*/

	if (!getenv("ACEDB_NO_BANNER"))
    	messout ("À bientot!... \n\n") ;

	// Memory leak checking
#if defined(_DEBUG) && defined(_DUMP_MEM_STATS)
	after.Checkpoint() ;
	if(diff.Difference( before, after ) )
	{
		TRACE("\nMemory Leaks!!!\n") ;
		diff.DumpStatistics() ;
		before.DumpAllObjectsSince() ;
	}
	AfxEnableMemoryTracking(FALSE) ; // To prevent the post-mortem dump?
#endif

	// Otherwise, exit the WinAce shell too?
	theApp.m_pMainWnd->SendMessage(WM_CLOSE) ;

#ifdef _DEBUG

#endif
}

extern "C" void wormClose(void) /* normal exit */
{
	winExit(0) ;
}

extern "C" void winCrash(void) /* unconditional error exit */
{
	ASSERT(FALSE) ;
	winExit(1) ;
}

void CWinAceApp::OnAppExit() 
{
	wormClose() ;	
}

//****** graphLoop() management, including main loop? ********************************
// Modified from CWinThread to handle graphLoop()

extern void AFXAPI _AfxTraceMsg(LPCTSTR lpszPrefix, const MSG* pMsg) ;
extern BOOL InGraphLoop() ;
extern int	ExitGraphLoop() ;
extern BOOL LoopIsDone() ;

BOOL CWinAceApp::PumpMessage()
{
	ASSERT_VALID(this);

#ifdef _DEBUG
	if (m_nDisablePumpCount != 0)
	{
		TRACE0("Error: CWinAceApp::PumpMessage called when not permitted.\n");
		ASSERT(FALSE);
	}
#endif

	if (!::GetMessage(&m_msgCur, NULL, NULL, NULL))  // !Message => WM_QUIT
	{
#ifdef _DEBUG
		if (afxTraceFlags & traceAppMsg)
			TRACE0("CWinAceApp::PumpMessage - Received WM_QUIT.\n");

		// application must die if not in recursive graphLoop()
		// Note: prevents calling message loop things in 'ExitInstance'
		// will never be decremented
		if( !InGraphLoop() ) m_nDisablePumpCount++;
#endif
		return FALSE; 
	}

#ifdef _DEBUG
	if (afxTraceFlags & traceAppMsg)
		_AfxTraceMsg(_T("PumpMessage"), &m_msgCur);
#endif

	// process this message
	if (m_msgCur.message != WM_KICKIDLE && !PreTranslateMessage(&m_msgCur))
	{
		::TranslateMessage(&m_msgCur) ;
		::DispatchMessage(&m_msgCur);
	}
	return TRUE;
}

// Main running routine until application exits
// Modified from CWinThread && CWinApp code
// to handle recursive graphLoop() with or without blocking
int CWinAceApp::Run()
{
	ASSERT_VALID(this);

	if (m_pMainWnd == NULL && AfxOleGetUserCtrl())
	{
		// Not launched /Embedding or /Automation, but has no main window!
		TRACE0("Warning: m_pMainWnd is NULL in CWinApp::Run - quitting application.\n");
		AfxPostQuitMessage(0);
	}

	// for tracking the idle time state
	BOOL bIdle = TRUE;
	LONG lIdleCount = 0;

	// acquire and dispatch messages until a WM_QUIT message is received.
	for (;;)
	{
		// phase1: check to see if we can do idle work
		while (bIdle &&
			!::PeekMessage(&m_msgCur, NULL, NULL, NULL, PM_NOREMOVE))
		{
			// call OnIdle while in bIdle state
			if (!OnIdle(lIdleCount++))
				bIdle = FALSE; // assume "no idle" state
		}

		// phase2: pump messages while available
		do
		{
			// pump message, but quit on WM_QUIT ... 
			if ( !PumpMessage() )
				return ExitInstance();
			// or if within a current graphLoop() which isDone
			else if( InGraphLoop() && LoopIsDone() )
				// Then, exit the most recent theApp->Run() in graphLoop()
				// returning the graphLoop()->retval
				return ExitGraphLoop() ;

			// reset "no idle" state after pumping "normal" message
			if (IsIdleMessage(&m_msgCur))
			{
				bIdle = TRUE;
				lIdleCount = 0;
			}

		} while (::PeekMessage(&m_msgCur, NULL, NULL, NULL, PM_NOREMOVE));
	}

	ASSERT(FALSE);  // not reachable
}

BOOL CWinAceApp::OnIdle(LONG lCount) 
{
	// Set main class pick window back to "Ready"
	// the first time OnIdle() is called each time?
#ifdef ACEDB
	if(lCount == 1L) mainActivity(0) ;
#endif
	
	return CWinApp::OnIdle(lCount);
}

#include "Preferences.h"

void CWinAceApp::OnGraphPreferences() 
{
    TRACE("\n*** Entering CWinAceApp::OnGraphPreferences():\n") ;
	CPreferences Preferences("WinAce Preferences") ; // With all pages

	if( Preferences.DoModal() ==  IDOK)
	{
		TRACE("\n*** Returning IDOK from Preferences Dialog in "
			  "CWinAceApp::OnGraphPreferences()\n" ) ;
		if( Preferences.ReloadACEDB() )
		{
			graphOut("Database profile changed. Exiting WinAce. Reload program to continue!") ;
			ACEDB_ExitCode = 1 ; // Trick closeACEDB into closing without prompting
			closeACEDB() ;		 // Close the existing database profile
			ACEDB_ExitCode = 0 ; // but, there are really no errors...yet
			m_pMainWnd->PostMessage(WM_CLOSE) ;
		}
	}
}
 
 
