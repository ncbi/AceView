/* $Id: dbprofile.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// DBProfile.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDBProfilePage dialog

class CDBProfilePage : public CPropertyPage
{
	DECLARE_DYNCREATE(CDBProfilePage)

// Construction
public:
	CDBProfilePage();
	~CDBProfilePage();

// Attributes
	BOOL m_ReloadACEDB ;

// Dialog Data
	//{{AFX_DATA(CDBProfilePage)
	enum { IDD = IDD_DB_PROFILE };
	CEdit	m_ACEDB_DATA;
	CEdit	m_ACEDB_COMMON;
	CEdit	m_ACEDB;
	CString	m_ACEDB_Directory;
	CString	m_ACEDB_DATA_Directory;
	CString	m_ACEDB_COMMON_Directory;
	//}}AFX_DATA


// Overrides
	// ClassWizard generate virtual function overrides
	//{{AFX_VIRTUAL(CDBProfilePage)
	public:
	virtual void OnOK();
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	// Generated message map functions
	//{{AFX_MSG(CDBProfilePage)
	virtual BOOL OnInitDialog();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

};
 
