/* TextScrollView.cpp : implementation of the CTextScrollGraphView class
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: TextScrollView graph implementation
 * HISTORY:
 * Last edited: 
 * Created: Jul 21 18:40 1995 (rbrusk)
 *-------------------------------------------------------------------
 */
// $Id: textscrollview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
//

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "FullScrollView.h"
#include "TextScrollView.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CTextScrollGraphView

IMPLEMENT_DYNCREATE(CTextScrollGraphView, CTextFullScrollGraphView)

BEGIN_MESSAGE_MAP(CTextScrollGraphView, CTextFullScrollGraphView)
	//{{AFX_MSG_MAP(CTextScrollGraphView)
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
// CTextScrollGraphView construction/destruction

CTextScrollGraphView::CTextScrollGraphView()
{
	// TODO: add construction code here

}

CTextScrollGraphView::~CTextScrollGraphView()
{
}

/////////////////////////////////////////////////////////////////////////////
// CTextScrollGraphView printing


/////////////////////////////////////////////////////////////////////////////
// CTextScrollGraphView diagnostics

#ifdef _DEBUG
void CTextScrollGraphView::AssertValid() const
{
	CGraphView::AssertValid();
}

void CTextScrollGraphView::Dump(CDumpContext& dc) const
{
	CGraphView::Dump(dc);
}

CGraph* CTextScrollGraphView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CGraph)));
	return (CGraph*)m_pDocument;
}

#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CTextScrollGraphView message handlers

void CTextScrollGraphView::OnScrollHome() 
{
	// DO NOTHING implementation overriding CTextFullScroll
	return ;
}

void CTextScrollGraphView::OnScrollEnd() 
{
	// DO NOTHING implementation overriding CTextFullScroll
	return ;
}

void CTextScrollGraphView::OnScrollLeft() 
{
	// DO NOTHING implementation overriding CTextFullScroll
	return ;
}

void CTextScrollGraphView::OnScrollRight() 
{
	// DO NOTHING implementation overriding CTextFullScroll
	return ;
}

void CTextScrollGraphView::OnScrollSOL() 
{
	// DO NOTHING implementation overriding CTextFullScroll
	return ;
}

void CTextScrollGraphView::OnScrollEOL() 
{
	// DO NOTHING implementation overriding CTextFullScroll
	return ;
}

BOOL CTextScrollGraphView::OnScrollBy(CSize sizeScroll, BOOL bDoScroll) 
{
	if(bDoScroll)
	{
		// Initial conditions
		CPoint spos = GetScrollPosition() ; 
		CRect cr ;
		GetClientRect(&cr) ;
		cr.OffsetRect(spos) ;  // Initial Viewport

		// Calculate new scroll position
		spos.y += sizeScroll.cy ;

		// Clip new scroll position to bounds
		if(spos.y < 0) spos.y = 0 ;
		if( spos.y > graphSize.cy ) spos.y = graphSize.cy ;

		// Calculate update drawing rectangle
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

BOOL CTextScrollGraphView::PreCreateWindow(CREATESTRUCT& cs) 
{
	// "TEXTSCROLL" graphs can only be vertically scrolled by user
	cs.style &= ~WS_HSCROLL ;
	cs.style |=  WS_VSCROLL ;
		
	return CTextFullScrollGraphView::PreCreateWindow(cs);
}
 
 
