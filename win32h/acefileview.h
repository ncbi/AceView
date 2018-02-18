/* $Id: acefileview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// AceFileView.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CAceFileView view

class CAceFileView : public CEditView
{
protected:
	CAceFileView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CAceFileView)

// Attributes
public:
	CAceFileDoc* GetDocument();

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAceFileView)
	protected:
	virtual void OnDraw(CDC* pDC);      // overridden to draw this view
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CAceFileView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CAceFileView)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in AceFileView.cpp
inline CAceFileDoc* CAceFileView::GetDocument()
   { return (CAceFileDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////
 
