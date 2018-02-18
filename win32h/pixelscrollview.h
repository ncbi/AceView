/* $Id: pixelscrollview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// PixelScrollView.h : interface of the CPixelScrollGraphView class
//
/////////////////////////////////////////////////////////////////////////////

class CPixelScrollGraphView : public CGraphView
{
protected:
	DECLARE_DYNCREATE(CPixelScrollGraphView)

// Attributes
public:
	CGraph* GetDocument();

// Operations

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPixelScrollGraphView)
	//}}AFX_VIRTUAL

// Implementation
	CPixelScrollGraphView();
	virtual ~CPixelScrollGraphView();

protected:

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif


// Generated message map functions
	//{{AFX_MSG(CPixelScrollGraphView)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in PixelScrollGraphView.cpp
inline CGraph* CPixelScrollGraphView::GetDocument()
   { return (CGraph*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
 
 
