/* $Id: pixelfitview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// PixelFitView.h : header file
//

#include "fitview.h"

/////////////////////////////////////////////////////////////////////////////
// CPixelFitGraphView view

class CPixelFitGraphView : public CFitGraphView
{
protected:
	DECLARE_DYNCREATE(CPixelFitGraphView)

// Attributes
public:
	CGraph* GetDocument();

// Operations

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPixelFitGraphView)
	protected:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	virtual void OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint);
	//}}AFX_VIRTUAL

public:
// Implementation
	CPixelFitGraphView() ;
	virtual ~CPixelFitGraphView();

protected:

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CPixelFitGraphView)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in PixelFitGraphView.cpp
inline CGraph* CPixelFitGraphView::GetDocument()
   { return (CGraph*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
 
 
