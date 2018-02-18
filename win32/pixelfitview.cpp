/* PixelFitView.cpp : implementation file
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: PixelFitView graph implementation
 * HISTORY:
 * Last edited: 
 * Created: Jul 21 18:40 1995 (rbrusk)
 *-------------------------------------------------------------------
 */
// $Id: pixelfitview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
//

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "PixelFitView.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CPixelFitGraphView

IMPLEMENT_DYNCREATE(CPixelFitGraphView, CFitGraphView)

CPixelFitGraphView::CPixelFitGraphView()
{
}

CPixelFitGraphView::~CPixelFitGraphView()
{
}


BEGIN_MESSAGE_MAP(CPixelFitGraphView, CFitGraphView)
	//{{AFX_MSG_MAP(CPixelFitGraphView)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CPixelFitGraphView diagnostics

#ifdef _DEBUG
void CPixelFitGraphView::AssertValid() const
{
	CFitGraphView::AssertValid();
}

void CPixelFitGraphView::Dump(CDumpContext& dc) const
{
	CFitGraphView::Dump(dc);
}

CGraph* CPixelFitGraphView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CGraph)));
	return (CGraph*)m_pDocument;
}

#endif // _DEBUG

/////////////////////////////////////////////////////////////////////////////
// CPixelFitGraphView message handlers


BOOL CPixelFitGraphView::PreCreateWindow(CREATESTRUCT& cs) 
{
	// "FIT" windows are not scrollable
	cs.style &= ~(WS_HSCROLL  |  WS_VSCROLL ) ;

	return CFitGraphView::PreCreateWindow(cs);
}

void CPixelFitGraphView::OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint) 
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
 
 
