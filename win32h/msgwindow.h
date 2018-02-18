/* $Id: msgwindow.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// msgwindow.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CMsgWindow dialog

class CGraphWindow ;

class CMsgWindow : public CDialog
{
private:
	CWnd *m_pOwner ;
// Construction
public:
	CMsgWindow( const CString& msgTitle, const CString& msgText, CWnd* pOwner );   // standard constructor

// Dialog Data
	//{{AFX_DATA(CMsgWindow)
	enum { IDD = IDD_MSGWINDOW };
	CEdit	m_MsgWindow;
	CString	m_MsgText;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMsgWindow)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual void PostNcDestroy();
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CMsgWindow)
	afx_msg void OnClose();
	afx_msg void OnFocus();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
 
 
