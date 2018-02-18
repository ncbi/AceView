/* $Id: fitview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// FitView.h : header file
//

#ifndef _FITVIEW_H_
#define _FITVIEW_H_

/////////////////////////////////////////////////////////////////////////////
// CFitGraphView view

class CFitGraphView : public CGraphView
{
protected:
	CFitGraphView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CFitGraphView)

// Attributes
public:
	CGraph* GetDocument();

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFitGraphView)
	public:
	virtual void OnInitialUpdate();
	protected:
	virtual void OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint);
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CFitGraphView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CFitGraphView)
	afx_msg void OnSize(UINT nType, int cx, int cy);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in FitGraphView.cpp
inline CGraph* CFitGraphView::GetDocument()
   { return (CGraph*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
 
#endif //NOT defined _FITVIEW_H_
 
