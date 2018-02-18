/* $Id: splashbox.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// SplashBox.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CSplashBox dialog

class CSplashBox : public CDialog
{
// Construction
public:
	CSplashBox(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CSplashBox)
	enum { IDD = IDD_SPLASH_BOX };
		// NOTE: the ClassWizard will add data members here
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSplashBox)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CSplashBox)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
 
