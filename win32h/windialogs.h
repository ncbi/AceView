/* $Id: windialogs.h,v 1.1.1.1 2002/07/19 20:23:27 sienkiew Exp $ */
// graphprompt.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CGraphPrompt dialog

class CGraphPrompt : public CDialog
{
// Construction
public:
	CGraphPrompt(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CGraphPrompt)
	enum { IDD = IDD_ACEDB_PROMPT };
	CString	m_userdata;
	CString	m_prompt;
	//}}AFX_DATA

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CGraphPrompt)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation

protected:

	// Generated message map functions
	//{{AFX_MSG(CGraphPrompt)
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
 
