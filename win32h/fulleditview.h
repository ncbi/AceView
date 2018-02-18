/* $Id: fulleditview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// FullEditView.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CTextFullEditGraphView view

class CTextFullEditGraphView : public CEditView
{
protected:
	CTextFullEditGraphView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CTextFullEditGraphView)

// Attributes
public:
	CGraph* GetDocument();

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTextFullEditGraphView)
	protected:
	virtual void OnDraw(CDC* pDC);      // overridden to draw this view
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CTextFullEditGraphView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CTextFullEditGraphView)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in FullEditGraphView.cpp
inline CGraph* CTextFullEditGraphView::GetDocument()
   { return (CGraph*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
 
 
