/* $Id: mapscrollview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// MapScrollView.h : interface of the CMapScrollGraphView class
//
/////////////////////////////////////////////////////////////////////////////

class CMapScrollGraphView : public CGraphView
{
protected:
	DECLARE_DYNCREATE(CMapScrollGraphView)

// Attributes
public:
	CGraph* GetDocument();

// Operations

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMapScrollGraphView)
	//}}AFX_VIRTUAL

// Implementation
	CMapScrollGraphView();
	virtual ~CMapScrollGraphView();

protected:

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

// Generated message map functions
	//{{AFX_MSG(CMapScrollGraphView)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in MapScrollGraphView.cpp
inline CGraph* CMapScrollGraphView::GetDocument()
   { return (CGraph*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
 
 
