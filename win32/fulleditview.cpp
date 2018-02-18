// $Id: fulleditview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
// FullEditView.cpp : implementation file
//

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "FullEditView.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CTextFullEditGraphView

IMPLEMENT_DYNCREATE(CTextFullEditGraphView, CRichEditView)

CTextFullEditGraphView::CTextFullEditGraphView()
{
}

CTextFullEditGraphView::~CTextFullEditGraphView()
{
}


BEGIN_MESSAGE_MAP(CTextFullEditGraphView, CRichEditView)
	//{{AFX_MSG_MAP(CTextFullEditGraphView)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CTextFullEditGraphView drawing

void CTextFullEditGraphView::OnDraw(CDC* pDC)
{
	CDocument* pDoc = GetDocument();
	// TODO: add draw code here
	CRichEditView::OnDraw(pDC) ;
}

/////////////////////////////////////////////////////////////////////////////
// CTextFullEditGraphView diagnostics

#ifdef _DEBUG
void CTextFullEditGraphView::AssertValid() const
{
	CRichEditView::AssertValid();
}

void CTextFullEditGraphView::Dump(CDumpContext& dc) const
{
	CRichEditView::Dump(dc);
}

CTextFullEditGraphDoc* CTextFullEditGraphView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CTextFullEditGraphDoc)));
	return (CTextFullEditGraphDoc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CTextFullEditGraphView message handlers
 
 
