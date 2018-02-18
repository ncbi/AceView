/* $Id: fullscrollview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// FullScrollView.h : interface of the CTextFullScrollGraphView class
//
/////////////////////////////////////////////////////////////////////////////

class CTextFullScrollGraphView : public CGraphView
{
protected:
	DECLARE_DYNCREATE(CTextFullScrollGraphView)

public:
// Attributes
	CGraph* GetDocument();

// Operations

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTextFullScrollGraphView)
	protected:
	virtual BOOL OnScrollBy(CSize sizeScroll, BOOL bDoScroll = TRUE);
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	//}}AFX_VIRTUAL

public:
// Implementation
	CTextFullScrollGraphView();
	virtual ~CTextFullScrollGraphView();

protected:

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

// Generated message map functions
	//{{AFX_MSG(CTextFullScrollGraphView)
	afx_msg void OnLineDown();
	afx_msg void OnLineUp();
	afx_msg void OnPageDown();
	afx_msg void OnPageUp();
	afx_msg void OnPageTop();
	afx_msg void OnPageBottom();
	afx_msg void OnScrollHome();
	afx_msg void OnScrollEnd();
	afx_msg void OnScrollLeft();
	afx_msg void OnScrollRight();
	afx_msg void OnScrollSOL();
	afx_msg void OnScrollEOL();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in TextFullScrollGraphView.cpp
inline CGraph* CTextFullScrollGraphView::GetDocument()
   { return (CGraph*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
 
 
