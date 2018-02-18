/* $Id: textfitview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// TextFitView.h : header file
//

#include "fitview.h"

/////////////////////////////////////////////////////////////////////////////
// CTextFitGraphView view

class CTextFitGraphView : public CFitGraphView
{
protected:
	DECLARE_DYNCREATE(CTextFitGraphView)

// Attributes
public:
	CGraph* GetDocument();

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTextFitGraphView)
	protected:
	virtual void OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint);
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	//}}AFX_VIRTUAL

public:
// Implementation
	CTextFitGraphView();
	virtual ~CTextFitGraphView();

protected:

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CTextFitGraphView)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in TextFitGraphView.cpp
inline CGraph* CTextFitGraphView::GetDocument()
   { return (CGraph*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
 
 
