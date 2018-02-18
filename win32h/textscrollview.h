/* $Id: textscrollview.h,v 1.1.1.1 2002/07/19 20:23:27 sienkiew Exp $ */
// TextScrollView.h : interface of the CTextScrollGraphView class
//
/////////////////////////////////////////////////////////////////////////////

class CTextScrollGraphView : public CTextFullScrollGraphView
{
protected:
	DECLARE_DYNCREATE(CTextScrollGraphView)

// Attributes
public:
	CGraph* GetDocument();

// Operations

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTextScrollGraphView)
	protected:
	virtual BOOL OnScrollBy(CSize sizeScroll, BOOL bDoScroll = TRUE);
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	//}}AFX_VIRTUAL

// Implementation
public:
	CTextScrollGraphView();
	virtual ~CTextScrollGraphView();

protected:

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

// Generated message map functions
	//{{AFX_MSG(CTextScrollGraphView)
	afx_msg void OnScrollHome();
	afx_msg void OnScrollEnd();
	afx_msg void OnScrollLeft();
	afx_msg void OnScrollRight();
	afx_msg void OnScrollSOL();
	afx_msg void OnScrollEOL();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in TextScrollGraphView.cpp
inline CGraph* CTextScrollGraphView::GetDocument()
   { return (CGraph*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
 
 
