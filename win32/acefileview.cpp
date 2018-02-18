// $Id: acefileview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
// AceFileView.cpp : implementation file
//

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "AceFileDoc.h"
#include "AceFileView.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CAceFileView

IMPLEMENT_DYNCREATE(CAceFileView, CEditView)

CAceFileView::CAceFileView()
{
}

CAceFileView::~CAceFileView()
{
}


BEGIN_MESSAGE_MAP(CAceFileView, CEditView)
	//{{AFX_MSG_MAP(CAceFileView)
	// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CAceFileView drawing

void CAceFileView::OnDraw(CDC* pDC)
{
	CEditView::OnDraw(pDC) ;
}

/////////////////////////////////////////////////////////////////////////////
// CAceFileView diagnostics

#ifdef _DEBUG
void CAceFileView::AssertValid() const
{
	CEditView::AssertValid();
}

void CAceFileView::Dump(CDumpContext& dc) const
{
	CEditView::Dump(dc);
}

CAceFileDoc* CAceFileView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CAceFileDoc)));
	return (CAceFileDoc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CAceFileView message handlers

BOOL CAceFileView::OnPreparePrinting(CPrintInfo* pInfo) 
{
	// TODO: Add your specialized code here and/or call the base class
	
	return CEditView::OnPreparePrinting(pInfo);
}
 
