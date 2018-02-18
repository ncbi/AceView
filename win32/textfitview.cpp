/* FitView.cpp : implementation of the CTextFitGraphView class
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: TextFitView graph implementation
 * HISTORY:
 * Last edited: 
 * Created: Jul 21 18:40 1995 (rbrusk)
 *-------------------------------------------------------------------
 */
// $Id: textfitview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
//

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "TextFitView.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CTextFitGraphView

IMPLEMENT_DYNCREATE(CTextFitGraphView, CFitGraphView)

CTextFitGraphView::CTextFitGraphView()
{
}

CTextFitGraphView::~CTextFitGraphView()
{
}


BEGIN_MESSAGE_MAP(CTextFitGraphView, CFitGraphView)
	//{{AFX_MSG_MAP(CTextFitGraphView)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CFitGraphView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CFitGraphView::OnFilePrintPreview)
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CTextFitGraphView diagnostics

#ifdef _DEBUG
void CTextFitGraphView::AssertValid() const
{
	CFitGraphView::AssertValid();
}

void CTextFitGraphView::Dump(CDumpContext& dc) const
{
	CFitGraphView::Dump(dc);
}

CGraph* CTextFitGraphView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CGraph)));
	return (CGraph*)m_pDocument;
}

#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CTextFitGraphView printing


/////////////////////////////////////////////////////////////////////////////
// CTextFitGraphView message handlers

// This OnUpdate() overrides the CFitGraphView() version
void CTextFitGraphView::OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint) 
{
/*
	CSize graphSize( uToXrel( gActive->uw ),uToYrel( gActive->uh ) ) ;
	SetScrollSizes(MM_TEXT, graphSize ) ;
	ResizeParentToFit() ;
*/
	if(GraphReSize)
	{
		SetScrollBounds(this) ;
		GraphReSize = FALSE ;
	}
	Invalidate( EraseGraph ) ;
	EraseGraph = FALSE ;
}

BOOL CTextFitGraphView::PreCreateWindow(CREATESTRUCT& cs) 
{
	// "FIT" windows are not scrollable
	cs.style &= ~(WS_HSCROLL  |  WS_VSCROLL ) ;

	return CFitGraphView::PreCreateWindow(cs);
}
 
 
