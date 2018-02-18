/* FitView.cpp : implementation file:
 * 
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: an abstract super class for fit graphs
 * HISTORY:
 * Last edited: 
 * Created: Jun 11 11:15 1996 (rbrusk): adapted from pixelfitgraphview.cpp
 *-------------------------------------------------------------------
 */
/* $Id: fitview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $ */


#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "FitView.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CFitGraphView

IMPLEMENT_DYNCREATE(CFitGraphView, CGraphView)

CFitGraphView::CFitGraphView()
{
}

CFitGraphView::~CFitGraphView()
{
}


BEGIN_MESSAGE_MAP(CFitGraphView, CGraphView)
	//{{AFX_MSG_MAP(CFitGraphView)
	ON_WM_SIZE()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CFitGraphView diagnostics

#ifdef _DEBUG
void CFitGraphView::AssertValid() const
{
	CGraphView::AssertValid();
}

void CFitGraphView::Dump(CDumpContext& dc) const
{
	CGraphView::Dump(dc);
}

CGraph* CFitGraphView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CGraph)));
	return (CGraph*)m_pDocument;
}

#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CFitGraphView message handlers


void CFitGraphView::OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint) 
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
 
 

void CFitGraphView::OnInitialUpdate() 
{
	CGraphView::OnInitialUpdate();
	
	// TODO: Add your specialized code here and/or call the base class
	
}

void CFitGraphView::OnSize(UINT nType, int cx, int cy) 
{
	CGraphView::OnSize(nType, cx, cy);
	
	// TODO: Add your message handler code here
	
}
 
