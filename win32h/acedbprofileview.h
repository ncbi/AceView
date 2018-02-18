/* $Id: acedbprofileview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// ACEDBProfileView.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CACEDBProfileView form view

#ifndef __AFXEXT_H__
#include <afxext.h>
#endif

class CACEDBProfileView : public CFormView
{
protected:
	CACEDBProfileView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CACEDBProfileView)

// Form Data
public:
	//{{AFX_DATA(CACEDBProfileView)
	enum { IDD = IDD_ACEDB_PROFILE };
		// NOTE: the ClassWizard will add data members here
	//}}AFX_DATA

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CACEDBProfileView)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CACEDBProfileView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CACEDBProfileView)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
 
