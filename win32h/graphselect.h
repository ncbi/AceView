/* $Id: graphselect.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// graphselect.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CSelector dialog

class CSelector : public CDialog
{
// Attributes
private:
	KEY m_selectedKey ;
	FREEOPT *m_options ;

// Construction
public:
	CSelector(FREEOPT *options, CWnd* pParent = NULL) ;		// standard constructor
	KEY GetSelectedKey() const { return m_selectedKey ; }

// Dialog Data
	//{{AFX_DATA(CSelector)
	enum { IDD = IDD_SELECTOR };
	CComboBox	m_keyList;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSelector)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CSelector)
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	afx_msg void OnSetFocusKeyList();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
 
