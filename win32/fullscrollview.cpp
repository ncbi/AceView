/* FullScrollView.cpp : implementation of the CTextFullScrollGraphView class
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: FullScrollView graph implementation
 * HISTORY:
 * Last edited: 
 * Created: Jul 21 18:40 1995 (rbrusk)
 *-------------------------------------------------------------------
 */
// $Id: fullscrollview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
//

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "FullScrollView.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CTextFullScrollGraphView

IMPLEMENT_DYNCREATE(CTextFullScrollGraphView, CGraphView)

BEGIN_MESSAGE_MAP(CTextFullScrollGraphView, CGraphView)
	//{{AFX_MSG_MAP(CTextFullScrollGraphView)
	ON_COMMAND(ID_LINE_DOWN, OnLineDown)
	ON_COMMAND(ID_LINE_UP, OnLineUp)
	ON_COMMAND(ID_PAGE_DOWN, OnPageDown)
	ON_COMMAND(ID_PAGE_UP, OnPageUp)
	ON_COMMAND(ID_PAGE_TOP, OnPageTop)
	ON_COMMAND(ID_PAGE_BOTTOM, OnPageBottom)
	ON_COMMAND(ID_SCROLL_HOME, OnScrollHome)
	ON_COMMAND(ID_SCROLL_END, OnScrollEnd)
	ON_COMMAND(ID_SCROLL_LEFT, OnScrollLeft)
	ON_COMMAND(ID_SCROLL_RIGHT, OnScrollRight)
	ON_COMMAND(ID_SCROLL_SOL, OnScrollSOL)
	ON_COMMAND(ID_SCROLL_EOL, OnScrollEOL)
	//}}AFX_MSG_MAP
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CGraphView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CGraphView::OnFilePrintPreview)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CTextFullScrollGraphView construction/destruction

CTextFullScrollGraphView::CTextFullScrollGraphView()
{
	// TODO: add construction code here
}

CTextFullScrollGraphView::~CTextFullScrollGraphView()
{
}

/////////////////////////////////////////////////////////////////////////////
// CTextFullScrollGraphView printing


/////////////////////////////////////////////////////////////////////////////
// CTextFullScrollGraphView diagnostics

#ifdef _DEBUG
void CTextFullScrollGraphView::AssertValid() const
{
	CGraphView::AssertValid();
}

void CTextFullScrollGraphView::Dump(CDumpContext& dc) const
{
	CGraphView::Dump(dc);
}

CGraph* CTextFullScrollGraphView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CGraph)));
	return (CGraph*)m_pDocument;
}
#endif / /_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CTextFullScrollGraphView message handlers

void CTextFullScrollGraphView::OnPageDown() 
{
	CSize sizeScroll(0,pageSize.cy) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnPageUp() 
{
	CSize sizeScroll(0,-pageSize.cy) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnLineDown() 
{
	CSize sizeScroll(0,lineSize.cy) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnLineUp() 
{
	CSize sizeScroll(0,-lineSize.cy) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnScrollHome() 
{
	CSize sizeScroll(-pageSize.cx, 0) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnScrollEnd() 
{
	CSize sizeScroll(pageSize.cx, 0) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnScrollLeft() 
{
	CSize sizeScroll(-lineSize.cx, 0) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnScrollRight() 
{
	CSize sizeScroll(lineSize.cx, 0) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnPageTop() 
{
	CPoint spos = GetScrollPosition() ;
	CSize sizeScroll(0,-spos.y) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnPageBottom() 
{
	CPoint spos = GetScrollPosition() ;
	CSize sizeScroll(0,graphSize.cy-spos.y) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnScrollSOL() 
{
	CPoint spos = GetScrollPosition() ;
	CSize sizeScroll(-spos.x,0) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

void CTextFullScrollGraphView::OnScrollEOL() 
{
	CPoint spos = GetScrollPosition() ;
	CSize sizeScroll(graphSize.cx-spos.x,0) ;
	OnScrollBy(sizeScroll, TRUE) ;
}

BOOL CTextFullScrollGraphView::OnScrollBy(CSize sizeScroll, BOOL bDoScroll) 
{
	if(bDoScroll)
	{
		// Initial conditions
		CPoint spos = GetScrollPosition() ; 
		CRect cr ;
		GetClientRect(&cr) ;
		cr.OffsetRect(spos) ;  // Initial Viewport

		// Calculate new scroll position
		spos.x += sizeScroll.cx ;
		spos.y += sizeScroll.cy ;

		// Clip new scroll position to bounds
		if(spos.x < 0) spos.x = 0 ;
		if(spos.x > graphSize.cx) spos.x = graphSize.cx ;
		if(spos.y < 0) spos.y = 0 ;
		if(spos.y > graphSize.cy ) spos.y = graphSize.cy ;

		// Calculate update drawing rectangle
		if( sizeScroll.cx ) 
		{
			if( sizeScroll.cx > 0 )
			{
				cr.left = cr.right ;
				cr.right += sizeScroll.cx ;
			}
			else /* sizeScroll.cx < 0  */
			{
				cr.right = cr.left ;
				cr.left += sizeScroll.cx ;
			}
		}
		if( sizeScroll.cy ) 
		{
			if( sizeScroll.cy > 0 )
			{
				cr.top = cr.bottom ;
				cr.bottom += sizeScroll.cy ;
			}
			else /* sizeScroll.cy < 0  */
			{
				cr.bottom = cr.top ;
				cr.top += sizeScroll.cy ;
			}
		}
		ScrollToPosition(spos) ;
		VIEWPTR->InvalidateRect(&cr,FALSE) ;
		return TRUE ;
	}
	return FALSE ;
}

BOOL CTextFullScrollGraphView::PreCreateWindow(CREATESTRUCT& cs) 
{
	// Full scrolling...
	cs.style |= (WS_HSCROLL  |  WS_VSCROLL ) ;

	return CGraphView::PreCreateWindow(cs);
}
 
 
