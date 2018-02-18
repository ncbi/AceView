/* $Id: mainfrm.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// mainfrm.h : interface of the CMainFrame class
//
/////////////////////////////////////////////////////////////////////////////

class CMainFrame : public CMDIFrameWnd
{
	DECLARE_DYNAMIC(CMainFrame)
public:
	CMainFrame();

protected:
	static int m_screenx, m_screeny ; // VGA pixel screen

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMainFrame)
	protected:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CMainFrame();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:  // control bar embedded members
	CStatusBar  m_wndStatusBar;
	CToolBar    m_wndToolBar;
	CDialogBar  m_mainDialogBar;

// Attributes
	Graph_ m_pACEDB_Graph ; // pointer to ACEDB Graph_ data structure

public:
    CDialogBar *GetMainDialogBar() { return &m_mainDialogBar; }
	void SetMainGraph( Graph_ g ) { m_pACEDB_Graph = g; }

// Generated message map functions
protected:
	//{{AFX_MSG(CMainFrame)
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnViewMainWindow();
	afx_msg void OnPaletteChanged(CWnd* pFocusWnd);
	afx_msg BOOL OnQueryNewPalette();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnHelp();
	afx_msg void OnGraphPreferences();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
 
 
