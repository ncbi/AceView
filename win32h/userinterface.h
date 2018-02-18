/* $Id: userinterface.h,v 1.1.1.1 2002/07/19 20:23:27 sienkiew Exp $ */
// UserInterface.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CUserInterfacePage dialog

class CUserInterfacePage : public CPropertyPage
{
	DECLARE_DYNCREATE(CUserInterfacePage)

// Construction
public:
	CUserInterfacePage();
	~CUserInterfacePage();

// Dialog Data
	//{{AFX_DATA(CUserInterfacePage)
	enum { IDD = IDD_INTERFACE_LOOK };
	CEdit	m_WEB_BROWSER;
	CEdit	m_EMAIL_ADDRESS;
	BOOL	m_ACEDB_NO_BANNER;
	CString	m_EmailAddress;
	CString	m_WebBrowser;
	//}}AFX_DATA


// Overrides
	// ClassWizard generate virtual function overrides
	//{{AFX_VIRTUAL(CUserInterfacePage)
	public:
	virtual void OnOK();
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	// Generated message map functions
	//{{AFX_MSG(CUserInterfacePage)
	virtual BOOL OnInitDialog();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

};
 
