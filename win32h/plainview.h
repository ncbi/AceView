/* $Id: plainview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// PlainView.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CPlainGraphView view

class CPlainGraphView : public CGraphView
{
protected:
	DECLARE_DYNCREATE(CPlainGraphView)

// Attributes
public:
	CGraph* GetDocument();

	// Operations

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPlainGraphView)
	//}}AFX_VIRTUAL

// Implementation
	CPlainGraphView(); 
	virtual ~CPlainGraphView();

protected:

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CPlainGraphView)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in PlainGraphView.cpp
inline CGraph* CPlainGraphView::GetDocument()
   { return (CGraph*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
 
 
