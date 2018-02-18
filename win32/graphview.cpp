/* GraphView.cpp : implementation file
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: ACE Graph associated SubDev == "CScrollView" window implementation
 * HISTORY:
 * Last edited: Jun 21 06:15 1996 (rbrusk):
 *	-	Use gDC->SetBkMode( WIN32_TRANSPARENT ) for all graphics; fix OnPrepareDC()...
 *		(COLOR_SQUARES doesn't work with OPAQUE)
 *	-	Fix SetScrollBounds() switch for PIXEL_SCROLL (setting graphSize properly)
 *	-	Fix greyMap initialization and PIXELS_RAW display (almost)
 * * Jun 10 10:15 1996 (rbrusk):
 *		-	code reorganized a bit and a new architecture underway
 *		-	implemented graphWaitCursor() (the XWindows way; may not work well in WIN32?)
 * * Jun 9 10:15 1996 (rbrusk):
 *		-	merged graphwin32.cpp with graphview.cpp:
 * * Jun 8 14:00 1996 (rbrusk):
 *		-	merged graphwin32lib.cpp with graphview.cpp:
 *		-	The following variables and methods were absorbed into GraphView:
 *			-	oldfcol, oldbcol
			-	drawBox(),  setColors() 
 * * May 1 9:46 1996 (rbrusk): ACEKey in OnChar() and OnKeyDown()
 * * Apr 30 02:32 1996 (rbrusk): graphIsBlocked() for blocked loops
 * * Mar  23 14:40 1996 (rbrusk): Ace 4.2 upgrade
 * * Jan  2 02:39 1996 (rbrusk): Ace 4.1 upgrade
 * Created: Jul 21 18:40 1995 (rbrusk)
 *-------------------------------------------------------------------
 */
/*  Merged File: graphwin32lib.cpp
 *
 * Exported functions:  graphDevActivate, 
 			graphBoxDraw, graphRedraw, graphClipDraw
			graphXorLine, graphXorBox,
			gFontInfo, printFontInit, printFontFinish

 * HISTORY:
 * Last edited: Jun 7 16:12 1996 (rbrusk): define getVisual() here
 * * Jun 1-2 23:40 (rbrusk):
 * *	-	JIT Fonts, Pens & Brushes; pens with variable linewidths
 *			createAcePen() replaced with GetAcePen()
 * * May 31 10:05 1996 (rbrusk): WIN32 to 4.3
 * *	-	COLOR_SQUARES functionality updated
 * *	-	SetForeground() replaced by setColors()
 * * Apr 30 04:05 1996 (rbrusk): Initialize loopInfo in graphInit()
 * * Apr 21 00:35 1996 (rbrusk): Printing enhancements...
 * * Feb 20 02:39 1996 (rbrusk): General enhancements; WIN32 to 4.2
 * * Feb 3 23:45 1996 (rbrusk): Enhanced font management & display 
 * * Jul 9 22:00 1995 (rbrusk): converted to .cpp file; starting to
 *                                     fill in functions with Framework and CDC code
 * Recreated: Jul 4 00:20 1995 (rbrusk): Port to Win32 using graphxtlib.c
 *                                       as a starting point;
 *                                       all graphDev*() moved to WinAce.cpp
 ******************* from graphxlib.c **************************************
 * * Feb 19 12:21 1996 (daz)
 * * Jan 10 16:43 1996 (rd)
 * * Nov 14 17:41 1992 (rd): FILL_ARC, BackgroundImage (8/11/92)
 * * Aug 14 17:19 1992 (rd): COLOR_SQUARES
 * * Jul 25 12:19 1992 (mieg): TEXT_PTR_PTR
 * * Dec  4 13:55 1991 (mieg): ACEDB_COLOR to force isMono = FALSE
 * * Oct 14 15:56 1991 (rd): added clipping of all vectors, except
	only test the centre (top left) of circles, points, text.
 *-------------------------------------------------------------------
 *
 *  File: graphwin32.cpp (merged with graphview.cpp: 9/6/96)
 *
 * Description: API interface level of Win32 version of graph package
 * Exported functions: many
 * HISTORY:
 * Last Edited: Jun 3 13:30 1996 (rbrusk): WIN32 to 4.3
 *		-	graphCopy(), graphPaste(), graphPostBuffer(), graphPasteBuffer()
 * * May 9 22:10 1996 (rbrusk): Various updates; merged all device code
 * Recreated: Jul 10 22:04 1995 (rbrusk): Win32 (C++) Replacement for graphxt.c
 *-------------------------------------------------------------------
 */

/* $Id: graphview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $ */

#include "stdafx.h"
#include "winace.h"
#include "mainfrm.h"
#include "cgraph.h"

extern "C"
{
#include "key.h"
}

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif
#ifndef VERBOSE_DEBUG
#define VERBOSE_DEBUG
#endif
#define new DEBUG_NEW

/////////////////////////////////////////////////////////////////////////////
// CGraphView

IMPLEMENT_DYNCREATE(CGraphView, CScrollView)

#if defined(NEW_WIN32_GRAPHS)

CGraphView *CGraphView::gViewPort = 0 ;

// Default SubGraph Size at start of program is whole screen?
int CGraphView::m_def_left = 0 ;
int CGraphView::m_def_top = 0 ;
int CGraphView::m_def_width = CGraphView::gScreenX ;
int CGraphView::m_def_height = CGraphView::gScreenY ;

void CGraphView::InitSubDevStruct(Graph_ graph, int x, int y, int w, int h)
{
	// Set Graph (Child Window) Position and Size within MainFrame window

    m_def_left = (x>=0)?x:0;
	m_def_top  = (y>=0)?y:0;
	m_def_width  = graph->w	= (w>=0)?w:0 ;
	m_def_height = graph->h	= (h>=0)?h:0 ;

	// Graph type specific actions
	switch( graph->type )
	{
		// Set undefined uw and uh window size parameters
		// now that you have device info handy
		case PIXEL_FIT:
		case TEXT_FIT:
		case TEXT_SCROLL:
		case TEXT_FULL_SCROLL:
		case TEXT_FULL_EDIT:
			graph->uw = XtoUrel(graph->w) ; 
			graph->uh = YtoUrel(graph->h) ;
			break ;

		case PLAIN:
		case MAP_SCROLL:
		case PIXEL_SCROLL:
			break ;  // uw, uh already set elsewhere?

		default: 
			messcrash ("Invalid graph type %d requested in graphSubDevCreate()",graph->type) ;
	}
}

#endif // defined(NEW_WIN32_GRAPHS)

#if defined(NEW_WIN32_GRAPHS)
int CGraphView::nextControlID()
{
	return CGraphWindow::nextControlID() ;
}
#endif

CGraphView::CGraphView()
{
#if defined(NEW_WIN32_GRAPHS)
	// gActive is assumed to be set properly prior to CGraphView creation!?
	m_pACEDB_Graph = (CGraphStruct *)gActive ;

	m_MessageWndID = nextControlID() ;
	m_pMessageWnd = NULL ;
#endif

	EraseGraph = GraphReSize = TRUE ;
	ScrollBoundsSet = FALSE ;
	pageSize = lineSize = CSize(0,0) ;
    m_currentBrush = NULL ;
	m_myCursorOn = FALSE ;
	setMyCursor (Arrow) ;
	m_otherCursor = NULL ;
	MIDBUTTONFLAG = MK_SHIFT ;  // Shift + Left Button == Middle Button emulated?
	MBMode = DESTROY ; // Nonsense mouse button drag mode?
}

//////////////////////////////////////////////////////////////////////////////
// ACEDB/WIN32 box interface function

void CGraphView::DeleteImages()
{
	POSITION mapPos ;
	char *data ;
	LPIMAGESET pImageSet ;
	LPPIXELIMAGE pLast, pNext ;

	mapPos = GetImages().GetStartPosition() ;
	while( mapPos != NULL )
	{
		GetImages().GetNextAssoc (mapPos, data, pImageSet) ;

#if defined(VERBOSE_DEBUG)
		TRACE("Destroying Image(data == %p, pImageSet == %p)\n", data, pImageSet) ;
#endif

		pLast = pNext = pImageSet->next ;
		ASSERT(pNext) ; // always at least one bitmap?
		do {
			pNext = pNext->next ;
			delete pLast->pBitMap ;
			delete pLast ;
		} while(pLast = pNext) ;
		// assume that the pixel image memory block was malloc'd
		// by the caller with the expectation that I free it here
		free( data ) ;
		delete pImageSet ;
	}
}

#if defined(NEW_WIN32_GRAPHS)

void CSubGraphWindow::DeleteBoxViews()
{
	// Destroy all boxViews
	int maxbox = m_boxViews.GetUpperBound() ;
	for( int i = 0; i <= maxbox; ++i )
	{
		CGraphBox *GBox = m_boxViews.GetAt(i) ;
		delete GBox ;
	}
	m_boxViews.RemoveAll() ;
}

#endif // defined(NEW_WIN32_GRAPHS)

CGraphView::~CGraphView()
{
	DeleteImages() ;

#if defined(NEW_WIN32_GRAPHS)
	if( m_pMessageWnd )
	// in case there is a MESSAGE_DESTROY callback 
		CloseGraphMessage() ;

	DeleteBoxViews() ;
#endif // defined(NEW_WIN32_GRAPHS)
}

BEGIN_MESSAGE_MAP(CGraphView, CScrollView)
	//{{AFX_MSG_MAP(CGraphView)
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_RBUTTONDOWN()
	ON_WM_MENUSELECT()
	ON_WM_SIZE()
	ON_WM_LBUTTONDBLCLK()
	ON_WM_KEYDOWN()
	ON_WM_CHAR()
	ON_WM_CREATE()
	ON_COMMAND(ID_EDIT_COPY, OnEditCopy)
	ON_COMMAND(ID_EDIT_PASTE, OnEditPaste)
	ON_WM_MBUTTONDOWN()
	ON_WM_MBUTTONUP()
	ON_WM_NCHITTEST()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CGraphView class statics initialized

// Actual machine dependent WIN32 display characteristics
// gScreenX, gScreenY are the display dimensions in pixels
int 	CGraphView::gScreenX = 320,
	 	CGraphView::gScreenY = 200,
	 	CGraphView::gNumPixelBits = 8,
		CGraphView::gNumGraphPlanes = 1,
		CGraphView::gColourRes  ;

float CGraphView::gXAspect = 1.0 ;
float CGraphView::gYAspect = 1.0 ;
CGraphView::FONTPALETTE *CGraphView::gFont = &screenFont ;
CDC *CGraphView::gDC = NULL ;
HCURSOR CGraphView::Arrow = NULL , CGraphView::HourGlass = NULL ;

/////////////////////////////////////////////////////////////////////////////
// CGraphView characteristics

void CGraphView::GetGraphDevAttributes() 
{
	CDC IC ;
	if(!IC.CreateIC("DISPLAY","","CON:",NULL) )
		messcrash("GetGraphDevAttributes():Can't access display.drv characteristics?") ;

	if( (gNumPixelBits = IC.GetDeviceCaps(BITSPIXEL)) < 8 )
		messcrash("WinAce needs minimum 256 colour graphics mode to run!") ;

	gScreenX = IC.GetDeviceCaps(HORZRES) ;	
	gScreenY = IC.GetDeviceCaps(VERTRES) ;	
	gNumGraphPlanes = IC.GetDeviceCaps(PLANES) ;
	TRACE("\n\n*** GraphDevAttributes: %d x %d resolution, %d bits/pixel in %d planes\n",
			gScreenX,gScreenY,gNumPixelBits,gNumGraphPlanes) ;	
}

/////////////////////////////////////////////////////////////////////////////
// CGraphView drawing

CGraphWindow *CGraphView::GetGraphWindow()
{
	CGraphWindow *pGraphDev = (CGraphWindow *)GetParent() ;
	ASSERT_VALID( pGraphDev ) ;
	return( pGraphDev ) ;
}

#if defined(NEW_WIN32_GRAPHS)
CGraph *CGraphView::GetGraphDev()
{
	CGraph *pGraphDev = (CGraph *)GetDocument();
	ASSERT_VALID( pGraphDev ) ;
	return( pGraphDev ) ;
}
#endif //  defined(NEW_WIN32_GRAPHS)

#if defined(OLD_WIN32_GRAPHS)
CGraphStruct *CGraphView::GetGraph( BOOL flg ) // flg is a dummy parameter)
{
	return( GetGraphDev()->GetGraph() ) ;
}
#endif //  defined(OLD_WIN32_GRAPHS)

BOOL CGraphView::SetActiveGraph( BOOL ReDraw )
{
	CGraphStruct *pGraph = GetGraph() ;

	// FROZEN flag postpones transfer
	// of control from a graphDestroy() dying graph to another graph
	if( !FROZEN && pGraph )
	{
		gActive = (Graph_)pGraph ;

		gDev = gActive->dev ;  // is_a CGraph object
		gSubDev = gActive->subdev ;
		gStk = gActive->stack ;
		gBox = gActive->nbox ? gBoxGet (gActive->currbox) : 0 ;

#if defined(NEW_WIN32_GRAPHS)
		SETVIEWPORT 
#endif
		if( ReDraw ) graphRedraw() ;  // repaint the newly activated graph if indicated
		return TRUE ;
	}
	else
		return FALSE ; // if active graph not changed
}

void CGraphView::Activate(BOOL flg) // flg is a dummy parameter
{
	// Activate the parent frame window...
	CGraphWindow *pParentGraph =
		(CGraphWindow *)gActive->parent->subdev ;
	pParentGraph->MDIActivate() ;
	if( pParentGraph->IsIconic() )
		pParentGraph->MDIRestore() ;

	// ...but set this subgraph to be the active window?
	SetActiveWindow() ;
	
	SETVIEWPORT // Just in case...

}
/*********************** Main Drawing Routines *****************************/

extern BOOL blockRedraw ;  // see windialogs.cpp
extern BOOL AceIsAlive ;

// This routine is reentrant in the sense that previous
// gActive and gDC's are saved in block scope variables (on the stack)
void CGraphView::OnDraw(CDC* pDC)
{
	// Don't try to draw anything if ACEDB is closing down
	// or an ACEDB graphOut() or similar query is displayed
	if( blockRedraw || !AceIsAlive ) return ;

	ASSERT_VALID(pDC) ;

#if defined(REALLY_VERBOSE)
	TRACE("\n\t*** Entering OnDraw(graphId[%d])...\n",gActive->id) ;
#endif

	// Else, proceed with draw
	// First, save ACEDB graph state
	Graph_	oldGraph 		= gActive ;
	Dev		oldDev   		= gDev ;
	SubDev	oldSubDev   	= gSubDev ;

#if defined(NEW_WIN32_GRAPHS)
	// Currently active viewport; set whenever gDev and gSubDev are set?
	CGraphView *oldViewPort = gViewPort ; 
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this?
#endif

	Box     oldBox 			= gBox ;
	Stack   oldStk 			= gStk ;
	CDC *	oldDC			= gDC ;		// likely to be NULL?

	// and reset graph state to that associated
	// with this CGraphView, without redraw
	// and check if anything is on the drawing stack
	if( !GetGraphWindow()->SetActiveGraph(FALSE) || !gStk )
	{
		TRACE("CGraphView::OnDraw(): graph not redrawable") ;
	}
	else
	// All ACEDB graph state variables should be valid for drawing
	// at this point, so draw anything within the clipping region
	{
		CRect rc ;
		gDC = pDC ;

		int clipRetVal = gDC->GetClipBox( &rc ) ;
		if(clipRetVal == ERROR || clipRetVal == NULLREGION)
			goto OnDrawExit ;

#if defined(REALLY_VERBOSE)
		TRACE("\n\t\tLogical Clipping box = (%d,%d,%d,%d)\n",rc.left,rc.top,rc.right,rc.bottom ) ;
#endif

		// If the clipping region is visible? then, draw the
		// graph contents lying within defined clipping region on gDC
		if( rc.right > rc.left && rc.bottom > rc.top )
		{
			Box box0 = gBoxGet(0) ;
			InitDrawing( box0 ) ;
			drawBox( box0, rc ) ;  // the main show!
		}
	}

OnDrawExit:
	// Restore ACEDB graph state
	gDC		=	oldDC ;
	gActive	=	oldGraph ;
	gDev   	=	oldDev ;
	gSubDev	=	oldSubDev ;

#if defined(NEW_WIN32_GRAPHS)
	gViewPort = oldViewPort; 
#endif

	gBox 	=	oldBox;
	gStk 	=	oldStk;
}


/////////////////////////////////////////////////////////////////////////////
// CGraphView diagnostics

#ifdef _DEBUG
void CGraphView::AssertValid() const
{
	CScrollView::AssertValid();
}

void CGraphView::Dump(CDumpContext& dc) const
{
	CScrollView::Dump(dc) ;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CGraphView message handlers

extern CPalette gPalette ;

void CGraphView::OnActivateView(BOOL bActivate, CView* pActivateView, CView* pDeactiveView) 
{
	if( bActivate )
	{
		if( pActivateView == pDeactiveView )
		{
#if !defined(NEW_WIN32_GRAPHS)
			VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif
			// Application window is being reactivated with this same active view
			// Select/RealizePalette() here?
			CClientDC pDC( this ) ;
			ASSERT_VALID( &pDC ) ;

			pDC.SelectPalette(&gPalette, FALSE) ;
			pDC.RealizePalette() ;
		}
	} 
	
	CScrollView::OnActivateView(bActivate, pActivateView, pDeactiveView);
}


typedef void (*MouseFunc)(double x, double y) ;

void CGraphView::mouseCall( GraphEvent GEType, CPoint point )
{
	Graph_ pGraph = GetGraph() ;
	if( pGraph && pGraph->func[ GEType ] != NULL)
	{
		double ux = (double)(graphEventX = (float)XtoUabs( point.x ) ) ; 
		double uy = (double)(graphEventY = (float)YtoUabs( point.y ) ) ;
		( (MouseFunc)( pGraph->func[ GEType ] ) )(ux,uy) ;
	}
}

extern "C" void gLeftDown (float x, float y) ;
extern "C" void gMiddleDown (float x, float y) ;
extern "C" BOOL graphIsBlocked() ;

void CGraphView::OnLButtonDown(UINT nFlags, CPoint point) 
{
	if( graphIsBlocked() ) return ;

#if !defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	CClientDC dc(this) ;
	OnPrepareDC(&dc) ;
	dc.DPtoLP(&point) ;

	int old = graphActive() ;  // guard against change in active graph
	graphEventX = (float)XtoUabs( point.x ) ; 
	graphEventY = (float)YtoUabs( point.y ) ;

	if( nFlags & MIDBUTTONFLAG )
	{
		MBMode = MIDDLE_DOWN; 
		gMiddleDown( graphEventX, graphEventY ) ;
	}
	else
	{
		MBMode = LEFT_DOWN ;
		gLeftDown( graphEventX, graphEventY ) ;
	}
	if( !graphExists( old ) )
		return ;  // old may be destroyed at this point?  
	else
		SetCapture() ; // Get mouse movement, is case you are dragging

	// CScrollView::OnLButtonDown(nFlags, point);
}

void CGraphView::OnLButtonUp(UINT nFlags, CPoint point) 
{
	if( graphIsBlocked() ) return ;

#if !defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	CClientDC dc(this) ;
	OnPrepareDC(&dc) ;
	dc.DPtoLP(&point) ;

	// Check first whether or not this view captured the mouse
	if( GetCapture() == this )
	{
		ReleaseCapture() ;  // release mouse drag event capture

		// Don't bother checking for MIDBUTTONFLAG
		// but assume MBMode set in OnLButtonDown()
		switch( MBMode )
		{
			case LEFT_DOWN:
				mouseCall( LEFT_UP, point ) ;
				break ;
			case MIDDLE_DOWN: // Emulated middle button?
				mouseCall( MIDDLE_UP, point ) ;
				break ;
			default:  // do nothing?
				TRACE0("**** CGraphView::OnLButtonUp() non-fatal error condition? ****\n") ;
		}	
		MBMode = DESTROY ; // reset to nonsense mode
	}
	// else CScrollView::OnLButtonUp(nFlags, point) ; 	  // is this necessary?
}

void CGraphView::OnMButtonDown(UINT nFlags, CPoint point) 
{
	if( graphIsBlocked() ) return ;

#if !defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	CClientDC dc(this) ;
	OnPrepareDC(&dc) ;
	dc.DPtoLP(&point) ;

	MBMode = MIDDLE_DOWN; 

	graphEventX = (float)XtoUabs( point.x ) ; 
	graphEventY = (float)YtoUabs( point.y ) ;

	// graphsub.c implementation of LEFT_DOWN/PICK?

	int old = graphActive() ;  // guard against change in active graph
	gMiddleDown( graphEventX, graphEventY ) ;
	if( !graphExists( old ) ) return ;  // old may be destroyed at this point?  

	SetCapture() ; // Get mouse movement, is case you are dragging

	// CScrollView::OnMButtonDown(nFlags, point);
}

void CGraphView::OnMButtonUp(UINT nFlags, CPoint point) 
{
	if( graphIsBlocked() ) return ;

#if !defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	CClientDC dc(this) ;
	OnPrepareDC(&dc) ;
	dc.DPtoLP(&point) ;

	// Check first whether or not this view captured the mouse
	if( GetCapture() == this )
	{
		ReleaseCapture() ;  // release mouse drag event capture
		mouseCall( MIDDLE_UP, point ) ;
		MBMode = DESTROY ; // reset to nonsense mode
	}
	// else CScrollView::OnMButtonUp(nFlags, point) ; 	  // is this necessary?
}
void CGraphView::OnMouseMove(UINT nFlags, CPoint point) 
{
	if( graphIsBlocked() ) return ;

#if !defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	CClientDC dc(this) ;
	OnPrepareDC(&dc) ;
	dc.DPtoLP(&point) ;

	if( GetCapture() == this )
	{
		switch( MBMode )
		{
			case LEFT_DOWN:
				mouseCall( LEFT_DRAG, point ) ;
				break ;
			case MIDDLE_DOWN:
				mouseCall( MIDDLE_DRAG, point ) ;
				break ;
			default:  // do nothing?
				TRACE0("**** CGraphView::OnMouseMove() non-fatal error condition? ****\n") ;
		}
	}	
	// else CScrollView::OnMouseMove(nFlags, point); 	  // is this necessary?
}

#include "key.h"

//*******************************Keyboard related functions *******************************


static BOOL ShiftOn = FALSE ;
static BOOL ControlOn = FALSE ;
BOOL VK2AceKey(UINT &nChar )
{
	switch( nChar )
	{
		case VK_LEFT:
			nChar = LEFT_KEY;
			break ;
		case VK_RIGHT:
			nChar = RIGHT_KEY;
			break ;
		case VK_UP:
			nChar = UP_KEY;
			break ;
		case VK_DOWN:
			nChar = DOWN_KEY;
			break ;
		case VK_PRIOR:
			nChar = PAGE_UP_KEY;
			break ;
		case VK_NEXT:
			nChar = PAGE_DOWN_KEY;
			break ;
		case VK_HOME:
			nChar = HOME_KEY;
			break ;
		case VK_END:
			nChar = END_KEY;
			break ;
		case VK_INSERT:
			nChar = INSERT_KEY;
			break ;
		case VK_DELETE:
			nChar = DELETE_KEY;
			break ;
		default: return FALSE ;
	}
	return TRUE ;
}

// Catching
void CGraphView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
	if( graphIsBlocked() ) return ;
	
#if !defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	Graph_ pGraph = GetGraph() ;
	if( !pGraph ||								// if not valid graph
		pGraph->func[ KEYBOARD ] == NULL ) 	    // with no KEYBOARD callback registered
			return ;

#if defined(VERBOSE_DEBUG)
	TRACE("CGraphView::OnKeyDown(	nChar == %.4x | nFlags == %.4x)\n", nChar, nFlags) ;
	TRACE("\tnRepCnt == ") ;
#endif

	while( nRepCnt-- )
	{
#if defined(VERBOSE_DEBUG)
		TRACE("%.4x ",nRepCnt) ;
#endif
		// if is an extended key (not an alphanumeric character)
		// then send to keyboard callback,
		// otherwise, leave to base class to handle
		if( VK2AceKey( nChar ) )
			((KeyboardFunc)(pGraph->func[KEYBOARD]))((int)nChar) ;
		else
			CScrollView::OnKeyDown(nChar, nRepCnt, nFlags);
	}
#if defined(VERBOSE_DEBUG)
	TRACE("\n*** Exiting CGraphView::OnKeyDown()\n\n") ;
#endif
}

void Ext2AceKey(UINT &nChar )
{
#if defined(VERBOSE_DEBUG)
	TRACE("\n*** Entering Ext2AceKey(nChar == %.4x)\n",nChar) ;
#endif
	switch( nChar )
	{
		case 75:
			nChar = LEFT_KEY;
			break ;
		case 77:
			nChar = RIGHT_KEY;
			break ;
		case 72:
			nChar = UP_KEY;
			break ;
		case 80:
			nChar = DOWN_KEY;
			break ;
		case 73:
			nChar = PAGE_UP_KEY;
			break ;
		case 81:
			nChar = PAGE_DOWN_KEY;
			break ;
		case 71:
			nChar = HOME_KEY;
			break ;
		case 79:
			nChar = END_KEY;
			break ;
		case 82:
			nChar = INSERT_KEY;
			break ;
		case 83:
			nChar = DELETE_KEY;
			break ;
		default: /* do nothing */ ;
	}
}

void CGraphView::OnChar(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
	if( graphIsBlocked() ) return ;

#if !defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	Graph_ pGraph = GetGraph() ;
	if( !pGraph ||								// if not valid graph
		pGraph->func[ KEYBOARD ] == NULL ) 	    // with no KEYBOARD callback registered
			return ;

#if defined(VERBOSE_DEBUG)
	TRACE("CGraphView::OnChar(	nChar == %.4x | nFlags == %.4x)\n", nChar, nFlags) ;
	TRACE("\tnRepCnt == ") ;
#endif

	while( nRepCnt-- )
	{
#if defined(VERBOSE_DEBUG)
		TRACE("%.4x ",nRepCnt) ;
#endif
		if(  nFlags & KF_EXTENDED )
			Ext2AceKey( nChar ) ; // convert extended keys?
			 
		// Then send ACEDB the key
		( (KeyboardFunc)( pGraph->func[ KEYBOARD ] ) )( (int)nChar ) ;
	}
#if defined(VERBOSE_DEBUG)
	TRACE("\n*** Exiting CGraphView::OnChar()\n\n") ;
#endif
}

extern void gRightDown( CPoint point, CPoint origin ) ;

void CGraphView::OnRButtonDown(UINT nFlags, CPoint point) 
{
	if( graphIsBlocked() ) return ;

#if defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	CScrollView::OnRButtonDown(nFlags, point);
	CPoint origin = GetScrollPosition() ;
	gRightDown( point , origin ) ;  
}

extern int currentMenuBox ;
UINT selectedMenuItem = 0 ;

void CGraphView::OnMenuSelect(UINT nItemID, UINT nFlags, HMENU hSysMenu) 
{
	if( graphIsBlocked() ) return ;

	if( nFlags == 0xFFFF && !hSysMenu ) return ; // menu completion

	if( !(	nFlags &	MF_GRAYED 	 ||
			nFlags &	MF_DISABLED  ||
			nFlags &	MF_SEPARATOR ||
			nFlags &	MF_POPUP ) )
		selectedMenuItem = nItemID ;
	else
		selectedMenuItem = 0 ;  
}


void CGraphView::OnSize(UINT nType, int cx, int cy) 
{
	if( graphIsBlocked() ) return ;

	if( !gActive || !gDev || !gSubDev || FROZEN || gActive->isBlocked ) return ;

#if defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	if( nType == SIZE_MAXIMIZED || nType == SIZE_RESTORED )
	{
		GraphReSize = TRUE ;
		m_NewX = cx ; m_NewY = cy ;

        Graph_ pGraph = GetGraph() ;
        if( pGraph &&
                ( pGraph->type == TEXT_FIT  ||
                  pGraph->type == PIXEL_FIT   ) )
        {
                pGraph->w = cx ; // new width in pixels?
                pGraph->h = cy ; // new height in pixels?
                pGraph->uw = XtoUrel(cx) ;
                pGraph->uh = YtoUrel(cy) ;
                gUpdateBox0() ;

                if(pGraph->func[ RESIZE ] != NULL)
                                ( pGraph->func[ RESIZE ] )() ;
        } 
    }
    CScrollView::OnSize(nType, cx, cy);
}


void CGraphView::OnPrepareDC(CDC* pDC, CPrintInfo* pInfo) 
{

#if defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	CScrollView::OnPrepareDC(pDC, pInfo);
	pDC->SetBkColor( PALETTERGB(0xFF,0xFF,0xFF) ) ;
	pDC->SetBkMode( WIN32_TRANSPARENT ) ;	// because OPAQUE does not work?*/
	pDC->SetROP2( R2_COPYPEN ) ;
}

void CGraphView::OnInitialUpdate()
{

#if defined(NEW_WIN32_GRAPHS)
	VIEWPTR = this ;  // isa CGraphView under this model, so set to this
#endif

	GraphReSize = TRUE ;
	OnUpdate(NULL,0,NULL) ;
}


void CGraphView::ReSizeGraph()
{
	Graph_ pGraph = GetGraph() ;
	if( !FROZEN && pGraph &&
		( pGraph->type == TEXT_FIT  || 
		  pGraph->type == PIXEL_FIT   ) )
	{
		pGraph->w = m_NewX ; // new width in pixels?
		pGraph->h = m_NewY ; // new height in pixels?
		pGraph->uw = XtoUrel(m_NewX) ;
		pGraph->uh = YtoUrel(m_NewY) ;
		gUpdateBox0() ;

#if defined(VERBOSE_DEBUG)
		Box box = gBoxGet(0) ;
		TRACE(	"\n\t*** CGraphView::OnSize: Box[0](%4.2g,%4.2g,%4.2g,%4.2g)\n",
				box->x1,box->x2,box->y1,box->y2) ;
#endif

		if(pGraph->func[ RESIZE ] != NULL)
				( pGraph->func[ RESIZE ] )() ;
	}
}
		
void CGraphView::OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint) 
{
	if(GraphReSize)
	{
		GraphReSize = FALSE ;
		ReSizeGraph() ;
		SetScrollBounds(this) ;
	}
	else
		Invalidate( EraseGraph ) ;
	EraseGraph = FALSE ;
}

void CGraphView::OnLButtonDblClk(UINT nFlags, CPoint point) 
{
	if( graphIsBlocked() ) return ;

	// Simulate a double left button pick
	OnLButtonDown( nFlags,  point) ; 
	OnLButtonUp( nFlags,  point) ;
	OnLButtonDown( nFlags,  point) ; 
	OnLButtonUp( nFlags,  point) ;
	 
	CScrollView::OnLButtonDblClk(nFlags, point);
}

/******************** Printing Implementation ******************/
/* Note: The printing process is currently non-reentrant due to the 
   manner in which fonts are created and deleted */

#include "CAcePrintPage.h"

BOOL CGraphView::OnPreparePrinting(CPrintInfo* pInfo) 
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

extern void printInit() ;
extern void printFinish() ;

void CGraphView::OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo) 
{
	// Extract the particulars about the printing options
	LPDEVMODE printer = pInfo->m_pPD->GetDevMode() ;
	CAcePrintPage *aPP = new CAcePrintPage(printer) ;

	// The MS MFC documentation says to GlobalFree this but 
	// this crashes the program... Q: If I ignore this requirement
	// will I get an intolerable memory leak?
	// ASSERT( GlobalFree((HGLOBAL)printer) == NULL ) ;

	pInfo->m_lpUserData = (LPVOID)aPP ;

	// Print GDI as a "fit-to-page" graphic using g?Aspect scaling
	// For future enhancements such as headers & footers, need to use
	// pInfo->m_rectDraw instead of aPP.Page* variables, which
	// should rather be the initial full size of m_rectDraw
	CGraphView::gXAspect = (float)aPP->PageWidth() / (float)gActive->w ;
	CGraphView::gYAspect = (float)aPP->PageLength() / (float)gActive->h ;

	// pDC->SetMapMode (MM_ANISOTROPIC) ;
	// pDC->SetWindowExt ( aPP->PageWidth(), aPP->PageLength() ) ;
	// pDC->SetViewportExt ( gActive->w, gActive->h ) ;
	// graph package printInit() must be called
	// AFTER gXAspect and gYAspect are set
	printInit() ;

	CScrollView::OnBeginPrinting(pDC, pInfo);
}

void CGraphView::OnPrint(CDC* pDC, CPrintInfo* pInfo) 
{
	ASSERT_VALID(pDC) ;
	CGraphView::OnDraw(pDC) ;
}

void CGraphView::OnEndPrinting(CDC* pDC, CPrintInfo* pInfo) 
{
	CAcePrintPage *aPP = (CAcePrintPage *)pInfo->m_lpUserData ;
	delete aPP ;

	printFinish() ;
	// pDC->SetMapMode (MM_TEXT) ;

	CGraphView::gXAspect = 1.0 ;
	CGraphView::gYAspect = 1.0 ;

	CScrollView::OnEndPrinting(pDC, pInfo);
}

int CGraphView::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CScrollView::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	m_NewX = lpCreateStruct->cx ;
	m_NewY = lpCreateStruct->cy ;

	return 0;
}

/****************** 8/6/96:  Start of merged graphwin32lib.cpp code *************/


// View aspect may differ between screen and printing (to fit)
// These macros provide the appropriate adjustments to logical coordinates
#define dw(x) (int)((x)*CGraphView::gXAspect)
#define dh(y) (int)((y)*CGraphView::gYAspect)

//***************************** Fonts... ****************************************************
/* font policy: use simple fonts so that they should all be there
   load them all in at beginning
   fontInit assumes default font is first in list
   note: this is a display specific object - so fontInit needs the
     display to have been opened.
*/

#define FONT(format,height) (*(CGraphView::gFont))[FNTIDX(format,height)]

// Entry name strings for Font descriptions in the INI registration database?
static CString fontFaceName[NUM_FONTFACES] = {
	"Plain","Fixed Width","Bold","Italic","Greek"
} ;

static CString fontSizeName[NUM_FONTSIZES] = {
	"Default","Tiny","Small","Medium","Large","Huge"
} ;

// Default Font faces
static CString defaultFontFace[NUM_FONTFACES] = {
	"Courier New",	/* Plain Font */
	"Courier New",	/* Fixed Width Font */
	"Courier New",	/* Bold Font */
	"Courier New",	/* Italic Font */
	"Symbol",		/* Greek Font */
} ;

// Default Font sizes
static int defaultFontSize[NUM_FONTSIZES] = // font heights only, widths are inferred
{
	12,			/* Default Font */
	8,			/* Tiny Font */
	10,			/* Small Font */
	12,			/* Medium Font */
	14,			/* Large Font */
	20			/* Huge Font */
} ;

CString	CGraphView::fontFace[NUM_FONTFACES] ; // Actual font face name strings
int		CGraphView::fontSize[NUM_FONTSIZES] ; // Actual font heights in logical units

CGraphView::FONTPALETTE CGraphView::screenFont;		// The "palette" of fonts available for screen display
CGraphView::FONTPALETTE CGraphView::printFont ;	// The "palette" of fonts available for printing on current page

static int	currFontdy = 0, currFontWidth = 0, currFontHeight= 0 ; // current font dimensions

static FONTFACE format2font (int format)
{
	static int oldFormat = PLAIN_FORMAT ;
	if (format >= 0)
		oldFormat = format ;
	else if(format > FIXED_WIDTH)
		format = PLAIN_FORMAT ;
	else
		format = oldFormat ;
	switch (format)
	{
		case PLAIN_FORMAT:
			return Plain ;
		case FIXED_WIDTH:
			return FixedWidth ;
		case ITALIC:
			return Italic ;
		case BOLD:
			return Bold ;
		case GREEK:
			return Greek ;
		default:
			TRACE("Invalid format code encountered in format2font()?\n") ;
			ASSERT(FALSE) ;
			return Plain ;
	}
}

static FONTSIZE height2font (int height)
{
	static float oldHeight = 1.;
	if (height >= 0)
		oldHeight = height;
	else
		height = oldHeight ;
	switch (height)
	{
		case 0:	 /* default font */
			return DefaultSize ;
		case 1: case 2: case 3: case 4: case 5: case 6:
			return Tiny ;
		case 7: case 8:
			return Small ;
		case 9: case 10:
			return Medium ;
		case 11: case 12:
			return Large ;
		case 13: case 14:
			return DefaultSize ; // ?? As per X Windows version
		default:
			return Huge ;
	}
	return DefaultSize ;
}

void CGraphView::createAceFont( CFont **font, int iFace, int iSize )
{
	// Destroy any previous font at that index
	if( font[FNTIDX(iFace,iSize)] != NULL) delete font[FNTIDX(iFace,iSize)] ;

	font[FNTIDX(iFace,iSize)] = new CFont() ;
	int h = dh(fontSize[iSize]) ;

	switch(iFace)
	{
		case Plain:
			if( font[FNTIDX(Plain,iSize)]->CreateFont( 
                        h,0,0,0,FW_NORMAL,0,0,0,
						ANSI_CHARSET,
						OUT_TT_PRECIS,
						CLIP_DEFAULT_PRECIS,
						PROOF_QUALITY,
						FIXED_PITCH | TMPF_TRUETYPE | FF_MODERN,
						(const char *)fontFace[Plain]) )
			return ;
			else break ;
		case FixedWidth:
			if( font[FNTIDX(FixedWidth,iSize)]->CreateFont(
                        h,0,0,0,FW_NORMAL,0,0,0,
						ANSI_CHARSET,
						OUT_TT_PRECIS,
						CLIP_DEFAULT_PRECIS,
						PROOF_QUALITY,
						FIXED_PITCH | TMPF_TRUETYPE | FF_MODERN,
						(const char *)fontFace[FixedWidth]) )
				return ;
			else break ;
		case Bold:
			if( font[FNTIDX(Bold,iSize)]->CreateFont(
						h,0,0,0,FW_BOLD,0,0,0,
						ANSI_CHARSET,
						OUT_TT_PRECIS,
						CLIP_DEFAULT_PRECIS,
						PROOF_QUALITY,
						FIXED_PITCH | TMPF_TRUETYPE | FF_MODERN,
						(const char *)fontFace[Bold]) )
				return ;
			else break ;
		case Italic:
			if( font[FNTIDX(Italic,iSize)]->CreateFont( 
                        h,0,0,0,FW_NORMAL,TRUE,0,0,
						ANSI_CHARSET,
						OUT_TT_PRECIS,
						CLIP_DEFAULT_PRECIS,
						PROOF_QUALITY,
						FIXED_PITCH | TMPF_TRUETYPE | FF_MODERN,
						(const char *)fontFace[Italic]) )
				return ;
			else break ;
		case Greek:
			if( font[FNTIDX(Greek,iSize)]->CreateFont(
                        h,0,0,0,FW_NORMAL,0,0,0,
						SYMBOL_CHARSET,
						OUT_TT_PRECIS,
						CLIP_DEFAULT_PRECIS,
						PROOF_QUALITY,
						FIXED_PITCH | TMPF_TRUETYPE | FF_DONTCARE,
						(const char *)fontFace[Greek]) )
				return ;
			else break ;
	
			/************TRYING ANOTHER FONT APPROACH *************

		case Plain:
			if( font[FNTIDX(Plain,iSize)]->CreatePointFont( 
                        fontSize[iSize]*10,
						(const char *)fontFace[Plain], gDC) ) 
				return ;
			else break ;
		case FixedWidth:
			if( font[FNTIDX(FixedWidth,iSize)]->CreatePointFont(
						fontSize[iSize]*10,
						(const char *)fontFace[FixedWidth], gDC) )
				return ;
			else break ;
		case Bold:
			if( font[FNTIDX(Bold,iSize)]->CreatePointFont(
						fontSize[iSize]*10,
						(const char *)fontFace[Bold], gDC) )
				return ;
			else break ;
		case Italic:
			if( font[FNTIDX(Italic,iSize)]->CreatePointFont( 
						fontSize[iSize]*10,
						(const char *)fontFace[Italic], gDC) )
				return ;
			else break ;
		case Greek:
			if( font[FNTIDX(Greek,iSize)]->CreatePointFont(
						fontSize[iSize]*10,
						(const char *)fontFace[Greek], gDC) )
				return ;
			else break ;

		***********************************************************/

		default:
			TRACE("Invalid font face passed to createAceFont()?\n") ;
			ASSERT(FALSE) ;
			iFace = iSize = 0 ;
	}
	// If i fall through here, this is a problem...
	messcrash(messprintf("Could not initialize %s %s font to %s of height %d\n",
							fontSizeName[iSize],fontFaceName[iFace],
							fontFace[iFace],fontSize[iSize])) ;
}

// The master font generator:
// Working on a "Just-in-Time" principal, this function
// returns an existing font, if it already exists or
// if the requested font is not yet created, this
// function creates it...
CFont *CGraphView::AceFont (FONTFACE face, FONTSIZE size )
{
	CFont *fnt ;
	if( (fnt = FONT( face, size )) == NULL )
	{
		TRACE("*** AceFont(): %s %s font created!\n",
				fontFaceName[face],fontSizeName[size]) ;
		createAceFont( screenFont, face, size ) ;
	}
	return ( FONT(face,size) ) ;
}

BOOL gFontInfo (int height, int* w, int* h)
{
	CFont *fnt = CGraphView::AceFont( format2font(-1), height2font(height) ) ;
	TEXTMETRIC fontMetric ;
	BOOL Success ;

	if( !CGraphView::gDC )
	{
		CWnd *pWnd ;
		if( gDev && VIEWPTR )
		{
			// Get the active Graph device?
			pWnd = (CWnd *)VIEWPTR ;
		}
		else
		{
			// Get the Application frame window
			pWnd = theApp.m_pMainWnd ;
		}
		TRY
		{
			CWindowDC *pGraphDC = new CWindowDC( pWnd ) ;
			CFont *oldFont = pGraphDC->SelectObject(fnt) ;
			Success = pGraphDC->GetTextMetrics( &fontMetric ) ; 
			pGraphDC->SelectObject(oldFont) ;
			delete pGraphDC ;
		}
		CATCH( CResourceException , e )
		{
			messcrash("WinAce Resource Exception: no device context for gFontInfo() ?") ;
		}
		END_CATCH

	}
	else
	{
		CFont *oldFont = CGraphView::gDC->SelectObject(fnt) ;
		Success = CGraphView::gDC->GetTextMetrics( &fontMetric ) ;
		CGraphView::gDC->SelectObject(oldFont) ;
	}


	if( Success )
	{
		// split the difference between Ave and Max character width
		if (w) *w = (int)( fontMetric.tmAveCharWidth );
		if (h) *h = (int)( fontMetric.tmAscent + fontMetric.tmDescent );
		return TRUE ;
	}
	else
		return FALSE ;
}

static int lastFontHeight ;

static CFont *fontSet (int height)
{
	ASSERT(CGraphView::gDC != NULL) ;
	TEXTMETRIC fontMetric ;
	
	CFont *fnt = CGraphView::AceFont( format2font(-1), height2font(height) ) ;
	if (!fnt) fnt = CGraphView::AceFont( Plain, DefaultSize ) ;

	// TRACE1("Calling fontSet(int %d) in graphWin32lib\n",height) ;
	CGraphView::gDC->SelectObject( fnt ) ;
	if( CGraphView::gDC->GetTextMetrics( &fontMetric ) )
	{ 
		currFontdy = fontMetric.tmAscent ;
		currFontWidth = fontMetric.tmAveCharWidth ;
		currFontHeight = fontMetric.tmAscent + fontMetric.tmDescent ;
	}
	else
		currFontdy = currFontWidth = currFontHeight = 0 ; 

	lastFontHeight = height ;

	return fnt ;
}

void CGraphView::graphFontInit()
{
	static int isDone = FALSE ; CString fontSizeStr ;

	if (isDone) return ;

	// Load the font descriptions: face and sizes
	for (int iFace = Plain ; iFace < NUM_FONTFACES ; ++iFace)
	{
		fontFace[iFace] =
			theApp.GetProfileString("WinAceFonts", 
					fontFaceName[iFace], defaultFontFace[iFace]) ;
	}

	for (int iSize = DefaultSize ; iSize < NUM_FONTSIZES ; ++iSize)
	{
		fontSize[iSize] =
			theApp.GetProfileInt("WinAceFonts", 
					fontSizeName[iSize], defaultFontSize[iSize]) ;
	}

	// JIT Fonts will be created as they are requested;
	// Just initialize font pointers to NULL for now...
	for ( iFace = Plain ; iFace < NUM_FONTFACES ; ++iFace)
		for ( iSize = DefaultSize ; iSize < NUM_FONTSIZES ; ++iSize)
		{
			screenFont[FNTIDX(iFace,iSize)] = NULL ; // To prevent deletion of an uninitialized font
			printFont[FNTIDX(iFace,iSize)]  = NULL ; // To prevent deletion of an uninitialized font
		}

	isDone = TRUE ;
}
			   
void CGraphView::graphFontFinish()
{
	for (int iFace = Plain ; iFace < NUM_FONTFACES ; ++iFace)
		if (!theApp.WriteProfileString("WinAceFonts", fontFaceName[iFace], fontFace[iFace] ) )
			TRACE(	"*** Warning: Font Face Description \"%s\" for the \"%s\" Font "
					" could not be written to program database?\n",
					fontFace[iFace], fontFaceName[iFace]) ;

	for (int iSize = DefaultSize ; iSize < NUM_FONTSIZES ; ++iSize)
		if (!theApp.WriteProfileInt("WinAceFonts", fontSizeName[iSize], fontSize[iSize] ) )
			TRACE(	"*** Warning: Font Size Descriptor \"%d\" for the \"%s\" Font "
					" could not be written to program database?\n",
					fontSize[iSize], fontSizeName[iSize]) ;

	for ( iFace = Plain ; iFace < NUM_FONTFACES ; ++iFace)
		for ( iSize = DefaultSize ; iSize < NUM_FONTSIZES ; ++iSize)
		{
			if( screenFont[FNTIDX(iFace,iSize)] != NULL)
				delete screenFont[FNTIDX(iFace,iSize)] ;
			if( printFont[FNTIDX(iFace,iSize)] != NULL) // not generally true...
				delete printFont[FNTIDX(iFace,iSize)] ;
		}
}			   

void CGraphView::printInit() // called by CGraphView::OnBeginPrint()
{
	for (int iFace = Plain ; iFace < NUM_FONTFACES ; ++iFace)
		for (int iSize = DefaultSize ; iSize < NUM_FONTSIZES ; ++iSize)
		{
			// Assume that you only need to create printer fonts
			// corresponding to JIT created (non-NULL) screen fonts?
			if( screenFont[FNTIDX(iFace,iSize)] != NULL )
			{
				// To prevent deletion of an uninitialized font
				printFont[FNTIDX(iFace,iSize)] = NULL ;

				// Initialize printer fonts scaled
				// to current printer paper size
				createAceFont( printFont, iFace, iSize ) ;
			}
		}
	gFont = &printFont ;
}			   

// Only recreates non-NULL fonts, that is, ones previously used
void CGraphView::reCreateFont( int iFace, int iSize )
{
	if(	screenFont[FNTIDX(iFace,iSize)] != NULL )
		createAceFont( screenFont, iFace, iSize ) ;
}

void CGraphView::printFinish() // called by CGraphView::OnEndPrint()
{
	gFont = &screenFont ;
	
	for (int iFace = Plain ; iFace < NUM_FONTFACES ; ++iFace)
		for (int iSize = DefaultSize ; iSize < NUM_FONTSIZES ; ++iSize)
			if(printFont[FNTIDX(iFace,iSize)] != NULL)
			{
				delete printFont[FNTIDX(iFace,iSize)] ;
				printFont[FNTIDX(iFace,iSize)] = NULL ;
			}
}			   

static CFont *fontSetFormat (int format) 
{
	ASSERT(CGraphView::gDC != NULL) ;
	TEXTMETRIC fontMetric ;
	CFont *fnt = CGraphView::AceFont( format2font(format),height2font(-1) ) ;
	if (!fnt) fnt = CGraphView::AceFont( Plain, DefaultSize ) ;

	switch (format)
	{
		case PLAIN_FORMAT:
		case FIXED_WIDTH:
		case ITALIC:
		case BOLD:
		case GREEK:
			CGraphView::gDC->SelectObject( fnt ) ;
			break ;
		default:
			messcrash("Invalid font format specified in drawBox()") ;
	}
	if( CGraphView::gDC->GetTextMetrics( &fontMetric ) )
	{ 
		currFontdy = fontMetric.tmAscent ;
		currFontWidth = fontMetric.tmAveCharWidth ;
		currFontHeight = fontMetric.tmAscent + fontMetric.tmDescent ;
	}
	else
		currFontdy = currFontWidth = currFontHeight = 0 ; 

	return fnt ;
}

/********** next color control - screen specific ***********/
#include "win32xcolor.h" // X Windows partial emulation types

static UCHAR greyMap[256] ;   

/* win32 specific:

typedef struct tagPALETTEENTRY {
    BYTE        peRed;
    BYTE        peGreen;
    BYTE        peBlue;
    BYTE        peFlags;
} PALETTEENTRY ;

enum Colour...

 */
#define PALETTESIZE	NUM_TRUECOLORS + /*  ACEDB Hard coded maximum # of .gif colours == */ 128 

PALETTEENTRY aceColorMap[PALETTESIZE] =
{  /* Initialize with the TRUECOLORS for now...*/
	{0xFF, 0xFF, 0xFF, 0x00}, // WHITE 
	{0x00, 0x00, 0x00, 0x00}, // BLACK
	{0xC0, 0xC0, 0xC0, 0x00}, // LIGHTGRAY
	{0x80, 0x80, 0x80, 0x80}, // DARKGRAY
	{0xFF, 0x00, 0x00, 0x00}, // RED
	{0x00, 0xFF, 0x00, 0x00}, // GREEN
	{0x00, 0x00, 0xFF, 0x00}, // BLUE
	{0xFF, 0xFF, 0x00, 0x00}, // YELLOW
	{0x00, 0xFF, 0xFF, 0x00}, // CYAN
	{0xFF, 0x00, 0xFF, 0x00}, // MAGENTA
	{0xFF, 0x80, 0x80, 0x00}, // LIGHTRED
	{0x80, 0xFF, 0x80, 0x00}, // LIGHTGREEN
	{0x80, 0x80, 0xFF, 0x00}, // LIGHTBLUE
	{0x80, 0x00, 0x00, 0x00}, // DARKRED
	{0x00, 0x80, 0x00, 0x00}, // DARKGREEN
	{0x00, 0x00, 0x80, 0x00}, // DARKBLUE
	{0xFF, 0xC0, 0xC0, 0x00}, // PALERED
	{0xC0, 0xFF, 0xC0, 0x00}, // PALEGREEN
	{0xC0, 0xC0, 0xFF, 0x00}, // PALEBLUE
	{0xFF, 0xFF, 0xC0, 0x00}, // PALEYELLOW
	{0xC0, 0xFF, 0xFF, 0x00}, // PALECYAN
	{0xFF, 0xC0, 0xFF, 0x00}, // PALEMAGENTA
	{0x80, 0x80, 0x00, 0x00}, // BROWN
	{0xFF, 0x80, 0x00, 0x00}, // ORANGE
	{0xFF, 0xC0, 0x00, 0x00}, // PALEORANGE
	{0xC0, 0x00, 0xC0, 0x00}, // PURPLE
	{0x80, 0x00, 0xC0, 0x00}, // VIOLET
	{0xC0, 0x00, 0xFF, 0x00}, // PALEVIOLET
	{0xA0, 0xA0, 0xA4, 0x00}, // GRAY
	{0xD0, 0xD0, 0xD8, 0x00}, // PALEGRAY
	{0xC0, 0x00, 0x00, 0x00}, // CERISE
	{0xA6, 0xCA, 0xF0, 0x00}  // MIDBLUE
} ; 

//Initialized in graphColorInit()...
// External for Win32 Active Color Palette
CPalette gPalette ;  

static char* colorName[] = 
  	{
  		"white","black","lightgray","darkgray",
   		"red","green","blue",
   		"yellow","cyan","magenta",
		"lightred","lightgreen","lightblue",
		"darkred","darkgreen","darkblue",
		"palered","palegreen","paleblue",
		"paleyellow","palecyan","palemagenta",
		"brown"," orange","paleorange",
		"purple","violet","paleviolet",
		"gray","palegray","cerise","midblue"
	} ;

extern "C" BOOL gWriteColors = FALSE ;

// Current colormap with values [0..255]
static unsigned char 	red[PALETTESIZE],
	 					green[PALETTESIZE],
	 					blue[PALETTESIZE] ;

// Normalized intensities of color components are [0..1],
// mapped from [0..255] colormap above
static float	nred[PALETTESIZE],
	 			ngreen[PALETTESIZE],
	 			nblue[PALETTESIZE] ;

// ACEDB JIT brush palette
static CBrush *aceBrush[PALETTESIZE] ;

// I assume here that CBrush colors/fill patterns don't change
// once the JIT CBrush is created... 
CBrush *CGraphView::GetAceBrush (int colour)
{
	if( aceBrush[colour] == NULL )
	{
		TRY // to create one
		{
			aceBrush[colour] = new CBrush(  PALETTEINDEX( colour ) );
		}
		CATCH( CResourceException , e )
		{
			messcrash("createAceBrush() could not create a CBrush?") ;
		}
		END_CATCH
	}
	return aceBrush[colour] ;
}

struct penBucket {
	CPen *pen ;
	int lineWidth ;
} ;

typedef penBucket *LPPENBUCKET ;
typedef CTypedPtrList<CPtrList,LPPENBUCKET> PENBUCKETBRIGADE ;
static PENBUCKETBRIGADE acePen[PALETTESIZE] ;

// First draft JIT colour pen manager creates  a new pen if NULL pen for 
// given colour OR recreates the pen if the pen is the wrong width
// (Note: a second version should keep track of variable sized pens)
CPen *CGraphView::GetAcePen (int colour)
{
	PENBUCKETBRIGADE &penSet = acePen[colour] ;
	LPPENBUCKET currPen = NULL , nextPen ;

	// look for a pre-existing pen of the given linewidth?
	if( !penSet.IsEmpty() )
	{
		POSITION pos = penSet.GetHeadPosition() ;
		while( pos != NULL )
		{
			nextPen = penSet.GetNext(pos) ;

			// a pen was found, but is is of the right linewidth?
			if( nextPen->lineWidth == getLineWidth() )
				{ currPen = nextPen ; break ; } // found
		}
	}

	if( currPen == NULL ) // pen not found, then create
	{
		TRACE(	"*** GetAcePen(): Creating an acePen of colour %d and width %d\n",
				colour, getLineWidth() ) ;
		TRY
		{
			// create the pen record...
			currPen = new penBucket ;
			currPen->lineWidth = getLineWidth() ;
			currPen->pen = new CPen( PS_SOLID, currPen->lineWidth, PALETTEINDEX( colour ) ) ;
			// ... then add it to the bucket brigade
			penSet.AddHead(currPen) ;

		}
		CATCH( CResourceException , e )
		{
			delete currPen ;
			messcrash("GetAcePen() could not create a CPen?") ;
		}
		END_CATCH
	}
	return currPen->pen ;
}


void CGraphView::graphColorInit ()
{
	if(gWriteColors) return ; // Don't reinitialize twice (i.e. on reloads)

	/**** First, initialize the palette *******/

	LOGPALETTE pal ;
	int i;
	float fac = 1.0 / 0xffff ;

	// aceColorMap is an array of PALETTEENTRY colour definitions defined in WinColor.h
	LPPALETTEENTRY pColours = aceColorMap ;
	for(i = 0; i <  PALETTESIZE; ++i )
	{
		if( i >= NUM_TRUECOLORS )  // Initialize upper palette with greyramp?
		{	int j = i - NUM_TRUECOLORS ;
			// Default to a greyscale ramp for upper palette
			pColours[i].peRed = pColours[i].peGreen = pColours[i].peBlue = 2*j ;
			pColours[i].peFlags = 0x00 ;
			greyMap[2*j] = greyMap[2*j+1] = i ;  // Assumes a linear ramp in intensities
		}
		nred[i] = fac * (red[i] = pColours[i].peRed) ;
		ngreen[i] = fac * (green[i] = pColours[i].peGreen) ;
		nblue[i] = fac * (blue[i] = pColours[i].peBlue) ;
	} 	
	pal.palVersion = 0x300 ;
	pal.palNumEntries = 1;
	if( !( gPalette.CreatePalette( &pal ) &&
	       gPalette.ResizePalette( PALETTESIZE )  &&
	       gPalette.SetPaletteEntries( 0, PALETTESIZE, pColours ) ) )
	    messcrash("Could not initialize the Windows color palette in graphColorInit()!") ;
	
	// Initialize the array of JIT CPens and CBrushes
	for( int colour = 0; colour < PALETTESIZE; ++ colour )
	{
		aceBrush[colour] = NULL ;	// First, initialize to NULL...
		// Note: acePen[colour] is part of a static
		// CTypePtrList object already initialized to empty
	}
	gWriteColors = TRUE ; // space in palette for GDI colors now available...
}

extern "C" BOOL getVisual(void)
{
	CGraphView::graphColorInit() ;	// if this doesn't messcrash... 
	return TRUE ;		//...then colours must be available!
}

void CGraphView::graphColorFinish ()
{
	// Initialize the array of CPens and CBrushes
	for( int colour = 0; colour < PALETTESIZE; ++ colour )
	{
		if(aceBrush[colour] != NULL )
			delete aceBrush[colour] ;
		
		PENBUCKETBRIGADE &penSet = acePen[colour] ;
		if(  !penSet.IsEmpty() )
		{
			// Destroy all pens in the penSet
			POSITION pos = penSet.GetHeadPosition() ;
			while( pos != NULL )
			{
				LPPENBUCKET deadPen = penSet.GetNext(pos) ;
				delete deadPen->pen ;
				delete deadPen ;
			}
		}
	}
}

/**************** pixel images ****************/

extern "C" BOOL graphRawMaps (unsigned char *forward, int *reverse) 
{
  int i ;

  if (!gWriteColors)
    return FALSE ;

  if (reverse)
    for (i = 0 ; i < 256 ; ++i)
      reverse[i] = 0 ;
  for (i = 0 ; i < 256 ; ++i)
    { if (forward)
	forward[i] = greyMap[i] ;
      if (reverse)
	reverse[greyMap[i]] = i ;
    }
  return TRUE ; 
}

extern "C" BOOL graphWritableColors (void)
{
    // Should be TRUE if graphColorInit() has been called
    return gWriteColors ;
}

extern "C" void graphColorMaps (float *red, float *green, float *blue) /* get */
{	 
	int i ;
	for (i = 0 ; i < 256 ; ++i)
	{
		red[i]		= nred[i] ;
		green[i]	= ngreen[i] ;
		blue[i]		= nblue[i] ;
	}
}

extern "C" BOOL graphSetColorMaps (unsigned char *redIntensity,
			unsigned char *greenIntensity,
			unsigned char *blueIntensity)
{
	LPPALETTEENTRY pColours = aceColorMap;
	int palColour, colour ;
	float fac = 1.0 / 0xffff ;

	for(colour = 0, palColour = NUM_TRUECOLORS; palColour < PALETTESIZE ; ++colour, ++palColour )
	{
		pColours[palColour].peRed 	= red[colour]	= redIntensity[colour] ;
		pColours[palColour].peGreen	= green[colour] = greenIntensity[colour] ;
		pColours[palColour].peBlue 	= blue[colour]	= blueIntensity[colour] ;
		pColours[palColour].peFlags	= 0x00 ;

		// Map onto normalized intensities of colour components
		nred[colour] = fac * red[colour] ;
		ngreen[colour] = fac * green[colour]  ;
		nblue[colour] = fac * blue[colour] ;
	} 	
	if( !( gPalette.SetPaletteEntries(NUM_TRUECOLORS , /*  ACEDB Hard coded maximum # of .gif colours == */ 128, pColours ) ) )
	{
	    	messout("Could not reset .gif colours in the Windows color palette in graphSetColorMaps()!") ;
		return  FALSE ;
	}

	return  TRUE ;
}

// This function provides a level of procedural indirection for palette setting
extern "C" void graphStoreColors (XColor *colors, int numColors) 
{	 
	 int i ;

	 // Coerce numColors into acceptable range
	 numColors = (numColors <= /*  ACEDB Hard coded maximum # of .gif colours == */ 128 )
	 		? numColors : /*  ACEDB Hard coded maximum # of .gif colours == */ 128 ;

	 // Separate out and set the color components
	 for(i = 0; i < numColors; ++i )
	 {
	 	red[i] = colors[i].red ;
	 	green[i] = colors[i].green ;
	 	blue[i] = colors[i].blue ;
	 }
	 graphSetColorMaps( red, green, blue ) ;
}

/***** WIN32 Colordialog box...
	CColorDialog cdlg;
	if( cdlg.DoModal() == IDOK)
	{
		COLORREF pixel = cdlg.GetColor() ;
		TRACE1("COLORREF == %d",pixel );
	}

*/
				
CBitmap * CGraphView::rawPixelImage(char *data, int w, int h, int len)
{
	LPIMAGESET pImageSet ;

#if defined(VERBOSE_DEBUG)
	TRACE("Entering rawPixelImage(data = %p, w == %d, h == %d, len == %d)\n",data, w, h, len) ;
#endif

	if( !GetImages().Lookup (data, pImageSet) )	// need to create a new image?
	{
		if (len % 4) 
			messcrash ("len for PixelsRaw = %d is not a multiple of 4", w) ;

		/* xim = XCreateImage (display, visual, 8, ZPixmap, 
		  0, data, 
		  (unsigned int) w, (unsigned int) h,
		  32, len) ;
		*/

		pImageSet = new IMAGESET ;
		pImageSet->data = data ;
		pImageSet->next = 0 ;
		(GetImages())[data] = pImageSet ;
	}

	LPPIXELIMAGE pBM = pImageSet->next ;
	while( pBM )
	{
		if( pBM->w == w &&
			pBM->h == h &&
			pBM->len == len )
			   break ;
		pBM = pBM->next ;
	}
	
	if(!pBM)
	{
		pBM = new PIXELIMAGE ;
		pBM->w = w ;
		pBM->h = h ; 
		pBM->len = len ;
		pBM->next = pImageSet->next ;
		pImageSet->next = pBM ; // chain new onto the head
		pBM->pBitMap = new CBitmap ;
		if( !pBM->pBitMap->CreateBitmap (w,h, gNumGraphPlanes, gNumPixelBits,NULL) )
		{
			messout("Sorrry... CGraphView::rawPixelImage() could not CreateBitmap()?") ;
			delete pBM->pBitMap; delete pBM ;
			return NULL ;
		}
	}
	// (re)set the bits in the image as required
	pBM->pBitMap->SetBitmapBits (w*h,(const void FAR *)data) ;
	return pBM->pBitMap ;
}

//********************* WIN32 Drawing Tools *******************************/

// ACEDB Rel 4.3 code: GraphView::setColors() replaces SetForeground() in a manner 
// similar to the use of setXcolors() w/XSet[Fore/Back]ground() in graphxlib.c 
// oldfcol and oldbcol are variables for use in certain functions

// Trick setColors() into unconditional reinitialization of pens and brushes 
// to the existing colours, by setting the "old" colours to complements of "new"
void CGraphView::setColors()
{
	unsigned char fcol = getFCol() , bcol = getBCol() ; // get existing colours...
	setFCol (~fcol) ; setBCol(~bcol) ;  // set to complement colours
	setColors ( fcol, bcol ) ;			// then reset colours
}

void CGraphView::setColors ( unsigned char fcol, unsigned char bcol )
{
	// Coerce colour values into proper range if necessary
	fcol = fcol < NUM_TRUECOLORS ? fcol: NUM_TRUECOLORS-1 ;
	bcol = bcol < NUM_TRUECOLORS ? bcol: NUM_TRUECOLORS-1 ;

	// Remember the current foreground and background colours as "old"
	oldfcol = getFCol() ; oldbcol = getBCol() ;

	// Obviously, only act when changed
	if( fcol != oldfcol ) 
	{
		setFCol (fcol) ;

		CPen *currPen = GetAcePen (fcol) ;
		VERIFY( ( gDC->SelectObject (currPen) ) != NULL) ;

		CBrush *currBrush = GetAceBrush (fcol) ;
		VERIFY( ( gDC->SelectObject (currBrush) ) != NULL) ;
        setBrush( currBrush ) ;

		gDC->SetTextColor ( PALETTEINDEX(fcol) ) ;
	}	
	if( bcol != oldbcol )
	{
		setBCol (bcol) ;
		gDC->SetBkColor ( PALETTEINDEX(bcol) ) ;
	}
}

/******* now draw boxes **********/

static int psize, psize2 ;

void CGraphView::InitDrawing(Box box)
{
	if ( !box ) return ;

	ASSERT_VALID( gDC ) ;
	gDC->SelectPalette(&gPalette, FALSE) ;
	gDC->RealizePalette() ;

	// Set the pen line width first...
	setLineWidth ( box->linewidth ) ;
	// then the drawing colours...
	setFCol (box->fcol) ; setBCol(box->bcol) ; setColors () ;

	fontSet (uToYrel(box->textheight)) ;
	fontSetFormat (PLAIN_FORMAT) ;
	psize = uToXrel(box->pointsize) ;
	psize2 = psize ? psize/2 : 1 ;
}

// Since the "FILL_RECTANGLE" type of functionality
// is reused several times in this code...
void CGraphView::drawFillRect(CRect *rect)
{
	// Since one presumes that the interest here is on the filling
	// and not on the drawing of the boundary of the rectangle,
	// you will likely always want to set a 1 pixel line width 
	int oldLW = setLineWidth() ;

	// Draw rectangle filled with current brush
	// Brush may be a weird one, like HOLLOW_BRUSH or XOR'd too?
	gDC->Rectangle( rect ) ;

	// then, reset line width...
	setLineWidth (oldLW) ;
}


#define xclip(z)  if (z < clip.left) z = clip.left ; \
		  else if (z > clip.right) z = clip.right
#define yclip(z)  if (z < clip.top) z = clip.top ; \
		  else if (z > clip.bottom) z = clip.bottom
#define xsafe(z)  if (z < -30000) z = -30000 ; \
		  else if (z > 30000) z = 30000
#define ysafe(z)  if (z < -30000) z = -30000 ; \
		  else if (z > 30000) z = 30000

void CGraphView::drawBox ( Box box, CRect &clip )
{
	float  	t ;
	int    	r, s ;
	char   	*text ;
	int    	action ;
	unsigned char fcol ;
#ifdef DEBUG_RECURSION
	static int recursionLevel = 0 ;
	++recursionLevel ;
#endif


	int	x1 = uToXabs(box->x1) ,
		y1 = uToYabs(box->y1) ,
		x2 = uToXabs(box->x2) ,
		y2 = uToYabs(box->y2) ;

	CRect boxRect(x1,y1, x2,y2) ;

	// If box is outside the clipping region, then just skip over
	// current box contents on stack, until stackAtEnd() or BOX_END

	if (box->x1 > box->x2 || box->y1 > box->y2 ||
		boxRect.left > clip.right ||
			boxRect.top > clip.bottom ||
				boxRect.right < clip.left ||
					boxRect.bottom  < clip.top )
	{
		int nDeep = 1 ;
		stackCursor (gStk,box->mark) ; /* sets position to mark */
		while (!stackAtEnd (gStk))
			switch (action = stackNext (gStk,int))
			{
				case BOX_END:
					if (!--nDeep)
#ifdef DEBUG_RECURSION
					TRACE("Level %d return non-draw box-end\n",
					 recursionLevel--) ;
#endif
					return ;                        /* exit point */
					break ;
				case BOX_START:
					r = stackNext (gStk, int) ;
					++nDeep ;
					break ;
				case COLOR: case TEXT_FORMAT:
					r = stackNext (gStk,int) ; 
					break ;
				case LINE_WIDTH: case TEXT_HEIGHT: case POINT_SIZE:
					t = stackNext (gStk,float) ;
					break ;
				case LINE: case RECTANGLE: case FILL_RECTANGLE:
					t = stackNext (gStk,float) ;
					t = stackNext (gStk,float) ;
					t = stackNext (gStk,float) ;
					t = stackNext (gStk,float) ;
					break ;
				case PIXELS: case PIXELS_RAW:
					t = stackNext (gStk,float) ;
					t = stackNext (gStk,float) ;
					text = stackNext (gStk,char*) ;
					r = stackNext (gStk, int) ;
					r = stackNext (gStk, int) ;
					r = stackNext (gStk, int) ;
					if (action == PIXELS)
					{
						t = stackNext (gStk,float) ;
						t = stackNext (gStk,float) ;
					}
					break ;
				case POLYGON :  case LINE_SEGS:
					r = stackNext (gStk, int) ;
					while (r--)
					{
						t = stackNext (gStk,float) ;
						t = stackNext (gStk,float) ;
					}
					break ;
				case CIRCLE: case POINT: 
				case TEXT: case TEXT_UP: case TEXT_PTR: case TEXT_PTR_PTR: 
				case COLOR_SQUARES: case FILL_ARC: case ARC:
					t = stackNext (gStk,float) ;
					t = stackNext (gStk,float) ;
					switch (action)
					{
						case CIRCLE:
							t = stackNext (gStk,float) ;
							break ;
						case FILL_ARC: case ARC:
							t = stackNext (gStk,float) ;
							t = stackNext (gStk,float) ;
							t = stackNext (gStk,float) ;
							break ;
						case POINT:
							break ;
						case TEXT: case TEXT_UP:
							text = stackNextText (gStk) ;
							break ;
						case TEXT_PTR:
							text = stackNext (gStk,char*) ;
							break ;
						case TEXT_PTR_PTR:
							text = *stackNext (gStk,char**) ;
							break ;
						case COLOR_SQUARES:
							text = stackNext (gStk,char*) ;
							r = stackNext (gStk, int) ;
							r = stackNext (gStk, int) ;
							text = (char*) stackNext (gStk,int*) ;
							break ;
						break ;
					}
					break ;
				case IMAGE:
					{
						void* gim = stackNext (gStk, void*) ;
					}
					break ;
				default:
					messout ("Invalid draw action %d received",action) ;
			}
#ifdef DEBUG_RECURSION
			TRACE("Level %d return bottom of stack non-draw\n", recursionLevel--) ;
#endif
			return ;  // return from invisible box
	}

	// gDC->SetBkMode( WIN32_TRANSPARENT ) ;	// because OPAQUE does not work..*/

	if (box->bcol != TRANSPARENT)
	{
		xclip(x1) ; xclip(x2) ; yclip(y1) ; yclip(y2) ;

		//  setXcolors(gc, box->bcol, box->bcol);
		setColors( box->bcol, box->bcol);
		CRect rect( dw(x1), dh(y1), dw(x2), dh(y2) ) ;

		// XFillRectangle (DWG, x1, y1, (x2-x1)+1, (y2-y1)+1) ;
		drawFillRect(&rect) ;
	}

	// color = colorInd[box->fcol & 0xf] ;
	// setXcolors(gc, fcol, box->bcol);
	fcol = box->fcol ;
	setColors( fcol, box->bcol);
	
	stackCursor (gStk,box->mark) ; /* sets position to mark */

	// Box assumed visible: render on device context!

	while (!stackAtEnd (gStk))
		switch (action = stackNext (gStk,int))
		{
			case BOX_END:
#ifdef DEBUG_RECURSION
				TRACE("Level %d return\n", recursionLevel--) ;
#endif
				return ;                        /* exit point */
			case BOX_START:
				r = stackNext (gStk, int) ;
#ifdef DEBUG_RECURSION
				TRACE("Level %d mark %d calls %d\n", 
					 recursionLevel, box->mark, r) ;
#endif
				drawBox ( gBoxGet (r), clip ) ;    /* recursion */
				// setXcolors(gc, fcol, box->bcol);
				setColors( fcol, box->bcol);
				break ;
			case COLOR:
				x1  = stackNext (gStk,int) ; 
				if (x1 == FORECOLOR) 
					fcol = box->fcol; 
				else if (x1 == BACKCOLOR)
					fcol = box->bcol;
				else
					fcol = x1;
				// setXcolors(gc, fcol, box->bcol);
				setColors( fcol, box->bcol);
				break ;
			case LINE_WIDTH:
				t = stackNext (gStk,float) ;
				setLineWidth( t ) ;
				break ;
			case TEXT_HEIGHT:
				t = stackNext (gStk,float) ; x1 = uToYrel(t) ;
				fontSet (x1) ;
				break ;
			case TEXT_FORMAT:
				x1 = stackNext (gStk,int) ; 
				fontSetFormat (x1) ;
				break ;
			case POINT_SIZE:
				t = stackNext (gStk,float) ;
				psize = dw(uToXrel(t)) ; psize2 = psize/2 ; 
				if (!psize2) psize2 = 1 ;
				break ;
			case LINE: case RECTANGLE: case FILL_RECTANGLE:
				t = stackNext (gStk,float) ; x1 = uToXabs(t) ; xsafe(x1) ; x1 = dw(x1) ;
				t = stackNext (gStk,float) ; y1 = uToYabs(t) ; ysafe(y1) ; y1 = dh(y1) ;
				t = stackNext (gStk,float) ; x2 = uToXabs(t) ; xsafe(x2) ; x2 = dw(x2) ;
				t = stackNext (gStk,float) ; y2 = uToYabs(t) ; ysafe(y2) ; y2 = dh(y2) ;
				switch (action)
				{
				case LINE:
					{
					CPoint start( x1, y1 ) ;
					CPoint finish( x2, y2 ) ;
					gDC->MoveTo ( start );
					gDC->LineTo ( finish );  // using gPen
					}
					break ;
				case RECTANGLE:
					if (x2 == x1) x2 = x1+1 ;
					if (y2 == y1) y2 = y1+1 ;
					if (x2 < x1) { r = x1 ; x1 = x2 ; x2 = r ; }
					if (y2 < y1) { r = y1 ; y1 = y2 ; y2 = r ; }
					{
						// XRectangle (DWG, x1, y1, (x2-x1), (y2-y1));
						CRect rect(x1,y1,x2,y2) ;

						// This older implementation works but is not completely to ACEDB spec
                        // in that printed versions of rectangles should have linewidth borders
						gDC->FrameRect ( &rect, getBrush() ) ; 

						/*  // This newer implementation does not work because both
						    // HOLLOW_BRUSH and NULL_BRUSH are not TRANSPARENT
						    CBrush *oldBrush =
                                (CBrush *)gDC->SelectStockObject( HOLLOW_BRUSH or NULL_BRUSH ) ;

						    // If you are NOT printing then use a pen linewidth
						    // of 1 pixel size to draw the rectangle
						    if(!gDC->IsPrinting())
							    drawFillRect (&rect) ;
						    else
							    gDC->Rectangle ( &rect ) ;

						    gDC-> SelectObject(oldBrush) ; // restore the old brush
                        */
					}
					break ;
				case FILL_RECTANGLE:
					if (x2 == x1) x2 = x1+1 ;
					if (x2 < x1) { r = x1 ; x1 = x2 ; x2 = r ; }
					if (y2 == y1) y2 = y1+1 ;
					if (y2 < y1) { r = y1 ; y1 = y2 ; y2 = r ; }

					// XFillRectangle (DWG, x1, y1, (x2-x1), (y2-y1));
					{
						CRect rect( x1, y1, x2, y2 ) ;
						drawFillRect (&rect) ;
					}
					break ;
				}
				break ;
			case PIXELS: case PIXELS_RAW:
				{ 
					int xbase, ybase, w, h, len ;
					char *pixels ;

					t = stackNext (gStk,float) ; x1 = uToXabs(t) ; xbase = dw(x1) ;
					t = stackNext (gStk,float) ; y1 = uToYabs(t) ; ybase = dh(y1) ;
					pixels = stackNext (gStk,char*) ;
					w = stackNext (gStk, int) ;
					h = stackNext (gStk, int) ;
					len = stackNext (gStk, int) ;
					if (action == PIXELS)
					{
						t = stackNext (gStk,float) ; x2 = uToXabs(t) ;
						t = stackNext (gStk,float) ; y2 = uToYabs(t) ;
						/* suz, xsafe in place of xclip in this case */
						xsafe(x1) ; xsafe(x2) ; ysafe(y1) ; ysafe(y2) ;
					}
					else
					{
						x2 = x1 + w ;
						y2 = y1 + h ;
						xclip(x1) ; xclip(x2) ; yclip(y1) ; yclip(y2) ;
					}
					if (x1 < x2 && y1 < y2)
					{ //if (action == PIXELS)
					// xim = pixelsImage (pixels, w, h, len, 
					//		   (x2-x1), (y2-y1)) ;
					//  else

					  if (action == PIXELS_RAW)
						{
							x1 = dw(x1) ; y1 = dh(y1) ;
							x2 = dw(x2) ; y2 = dh(y2) ;
							
							CBitmap *bm = rawPixelImage (pixels, w, h, len) ; 
							if (bm)
							{
					//			XPutImage (DWG, xim, x1-xbase, y1-ybase, 
					//	   					x1, y1, 
					//	   					(x2 < x1+w) ? x2-x1+1 : x2-x1, 
					//	   					(y2 < y1+h) ? y2-y1+1 : y2-y1) ;
					//
					//			This is my best guess on how this should work in WIN32?
					
								CDC mDC ;
								if(!( mDC.CreateCompatibleDC( gDC ) &&
									  mDC.SelectObject(bm) ) )
									 messout("drawBox() could not create and select PIXELS_RAW image?") ;
								else
								{
									if( CGraphView::gXAspect == 1.0 && CGraphView::gYAspect == 1.0 )
									{
										if(!( gDC->BitBlt(xbase, ybase,
											(x2 < x1+w) ? x2-x1+1 : x2-x1,
											(y2 < y1+h) ? y2-y1+1 : y2-y1,
											 &mDC, 0, 0, SRCCOPY) ) )
											messout("drawBox() could not BitBlt PIXELS_RAW image?") ;
									}
									else // need to stretch pattern to fit page scale
									{
										if(!( gDC->StretchBlt(xbase, ybase,
											(x2 < x1+w) ? x2-x1+1 : x2-x1,
											(y2 < y1+h) ? y2-y1+1 : y2-y1,
											 &mDC, 0, 0, w, h, SRCCOPY) ) )
											messout("drawBox() could not StretchBlt PIXELS_RAW image?") ;
									}
								}
							}
						}
					}
				}
				break ;
			case POLYGON :  case LINE_SEGS:
				{				   
					int numVertices = stackNext (gStk, int) ;
					{
						// XPoint *XPts = (XPoint *) messalloc (numVertices * sizeof (XPoint)) ;
						// XPoint *XPt = XPts ;
						LP_WIN_PT polyPts = (LP_WIN_PT) messalloc (numVertices * WIN_PT_SIZE) ;
						LP_WIN_PT polyPt  = polyPts ;

						for (r = 0 ; r < numVertices ; r++)
						{
							t = stackNext(gStk,float) ; 
							x1 = uToXabs(t) ; xsafe(x1) ;
							polyPt->x = dw(x1) ;

							t = stackNext(gStk,float) ; 
							y1 = uToYabs(t) ; ysafe(y1) ;
							polyPt->y = dh(y1) ;
							polyPt++;

						}
	    				if (action == POLYGON)
		    				/* NOTE this assumes the polygons are of the simplest
							 * possible type, with no intersecting lines */
							// XFillPolygon (DWG, XPts, numVertices, 
							//		Convex, CoordModeOrigin) ;
							gDC->Polygon( polyPts, numVertices );
						else
							gDC->Polyline( polyPts, numVertices );

						// messfree (XPts) ;
						messfree( polyPts ) ;

					}
				}
				break ;
			case CIRCLE: case POINT: 
			case TEXT: case TEXT_UP: case TEXT_PTR: case TEXT_PTR_PTR: 
			case COLOR_SQUARES: case FILL_ARC: case ARC:
				t = stackNext (gStk,float) ; x1 = uToXabs(t) ;
				t = stackNext (gStk,float) ; y1 = uToYabs(t) ;
				switch (action)
				{
					case CIRCLE:      /* given center x1, y1, radius r*/
						t = stackNext (gStk,float) ; r = uToXrel(t) ;
						if (!r) r = 1 ;
						if (x1 - r < clip.right && x1 + r > clip.left &&
							y1 - r < clip.bottom && y1 + r > clip.top)
						{
							// XDrawArc (DWG, x1-r, y1-r, 2*r, 2*r, 0, 64*360)
							CRect rect( dw(x1 - r), dh(y1 - r), dw(x1 + r), dh(y1 + r) ) ;
							gDC->Ellipse( &rect ) ;
						}
						break ;
					case ARC:
					case FILL_ARC:
						t = stackNext (gStk,float) ; r = uToXrel(t) ;
						{
							// Convert degrees to radians...
							// (2 * pi) / 360 = 0.0174533
							float ang =  0.0174533 * stackNext (gStk,float) ;
							float angDiff =  0.0174533 * stackNext (gStk,float) ;
							if (!r) r = 1 ;
							if (x1 - r < clip.right && x1 + r > clip.left &&
								y1 - r < clip.bottom && y1 + r > clip.top)
							{
								int xs, ys, xe, ye ; // Start and end points
								xs = x1 + (int)((double)r * cos( ang )) ;
								ys = y1 + (int)((double)r * sin( ang )) ;
								xe = x1 + (int)((double)r * cos( ang+angDiff )) ; 
								ye = y1 + (int)((double)r * sin( ang+angDiff )) ;
								if (action == ARC )
									// XDrawArc (DWG, x1-r, y1-r, 2*r, 2*r, 
									//	  (int)(ang*64), (int)(angDiff*64))
									gDC->Arc( dw(x1-r), dh(y1-r),
											  dw(x1+r), dh(y1+r),
											  dw(xs), dh(ys), 
											  dw(xe), dh(ye) ) ;
								else
									//XFillArc ( DWG, x1-r, y1-r, 2*r, 2*r, 
									//          (int)(ang*64), (int)(angDiff*64))
								   /* NB upper left corner and 64ths of a degree */
									gDC->Chord( dw(x1-r), dh(y1-r),
											    dw(x1+r), dh(y1+r),
											    dw(xs), dh(ys), 
											    dw(xe), dh(ye) ) ;
							}
						}
						break ;
					case POINT:
						if (x1 - psize2 < clip.right && x1 + psize2 > clip.left &&
							y1 - psize2 < clip.bottom && y1 + psize2 > clip.top)
						{
							//XFillRectangle (DWG, x1-psize2, y1-psize2, psize, psize)
							CRect rect(  dw(x1)-psize2, dh(y1)-psize2, dw(x1)+psize2, dh(y1)+psize2 ); 
							drawFillRect (&rect) ;
						}
						break ;
					case TEXT:
						text = stackNextText (gStk) ;
						s = strlen(text) ;
						if (x1 < clip.right && x1 + s*currFontWidth > clip.left &&
							y1 < clip.bottom && y1 + currFontHeight > clip.top)
							//XDrawString(DWG, x1, y1+currFontdy, text, s)
							gDC->TextOut(dw(x1), dh(y1), text, s ) ;
						break ;
				    case TEXT_UP:
				        text = stackNextText (gStk) ;
					    s = strlen(text) ;
					    { int oldHeight = lastFontHeight ;
					      int i ;
					      CFont *fnt = fontSet (oldHeight*0.6) ;
					      /* next lines remove line-spacing allowance */
					      /* fontHeight = fnt->max_bounds.ascent + fnt->max_bounds.descent ; */
	      
					      if (x1 < clip.right && x1 + currFontWidth > clip.left &&
						  y1 - s*currFontHeight < clip.bottom && y1 > clip.top)
						for (i = 0 ; i < s ; ++i)
						  //XDrawString (DWG, x1, (y1-(s-i-1)*fontHeight)-uToYabs(0.2),
						  //	       &text[i], 1)
						  gDC->TextOut( dw(x1), (dh(y1)-(s-i-1)*currFontHeight)-dh(uToYabs(0.2)),
						  		&text[i], 1)	;
					      fontSet (oldHeight) ;
					    }
				        break ;
					case TEXT_PTR:
						text = stackNext (gStk,char*) ;
						s = strlen(text) ;
						if (x1 < clip.right && x1 + s*currFontWidth > clip.left &&
							y1 < clip.bottom && y1 + currFontHeight > clip.top)
							//XDrawString(DWG, x1, y1+currFontdy, text, s)
							gDC->TextOut(dw(x1), dh(y1), text, s ) ;
						break ;
					case TEXT_PTR_PTR:
						if (!(text = *stackNext (gStk,char**)))
							break ;
						s = strlen(text) ;
						if (x1 < clip.right && x1 + s*currFontWidth > clip.left &&
							y1 < clip.bottom && y1 + currFontHeight > clip.top)
							//XDrawString(DWG, x1, y1+currFontdy, text, s)
							gDC->TextOut(dw(x1), dh(y1), text, s ) ;
						break ;
					case COLOR_SQUARES: {
						int *tints ;
						text = stackNext (gStk,char*) ;
						r = stackNext (gStk, int) ;
						s = stackNext (gStk, int) ;
						tints = stackNext (gStk, int*) ;
						x2 = 0 ;
						y2 = uToYrel (1.0) ;
						if ( x1 >= clip.right || x1 + r*uToXrel(1.0)  <= clip.left ||
							 y1 >= clip.bottom || y1 + y2 <= clip.top )
							break ;
						int oldColor = fcol ;
						while (r--)
						{
							int newColor;
							switch (*text^((*text)-1)) /* trick to find first on bit */
							{ 
								case -1:   newColor = WHITE; break;
								case 0x01: newColor = tints[0]; break;
								case 0x03: newColor = tints[1]; break;
								case 0x07: newColor = tints[2]; break;
								case 0x0f: newColor = tints[3]; break;
								case 0x1f: newColor = tints[4]; break;
								case 0x3f: newColor = tints[5]; break;
								case 0x7f: newColor = tints[6]; break;
								case 0xff: newColor = tints[7]; break;
							}

							if ( newColor != oldColor )
							{
								if (x2 > 0 )
								{
									//setXcolors(gc, oldColor, box->bcol) ;
									setColors( oldColor, box->bcol) ;

									//XFillRectangle (DWG, x1, y1, x2, y2) ;
									CRect rect( dw(x1), dh(y1), dw(x1 + x2), dh(y1 + y2) ) ;
									drawFillRect (&rect) ;
									x1 += x2; x2 = 0; 
								}
								oldColor = newColor ;
							}
							x2 += uToXrel(1.0) ;
							text += s ;
						}
						if (x2 > 0)
						{
							// setXcolors(gc, old, box->bcol) ;
							setColors( oldColor, box->bcol) ;

							CRect rect( dw(x1), dh(y1), dw(x1 + x2), dh(y1 + y2) ) ;
							
							//XFillRectangle(DWG, x1, y1, x2, y2)
							drawFillRect (&rect) ;
						} 
						// setXcolors(gc, fcol, box->bcol);
						setColors( fcol, box->bcol);
					} break ;
				}
				break ;
			case IMAGE:
				{
					void* gim = stackNext (gStk, void*) ;
					// What do I do here with an image?
				}
				break ;
			default:
				messout ("Invalid action %d received in drawBox()",action) ;
		}
#ifdef DEBUG_RECURSION
	TRACE("Level %d return at bottom of drawBox()\n", recursionLevel--) ;
#endif
}

extern "C" void graphWhiteOut (void)
{
	if ( !VIEWPTR ) return ;

	// XClearWindow (display, window) ; /* RMD 28/5/92 */
	// XFlush (display) ;

	VIEWPTR->EraseGraph = TRUE ;
}

// graphClipDraw() sets up an externally specified
// clipping rectangle for drawing box 0 of the graph
extern "C" void graphClipDraw (int x1, int y1, int x2, int y2) 
{
	if (!gActive->stack || !gDev || !VIEWPTR ) return ;

//#if defined(VERBOSE_DEBUG)
//	TRACE("\t\t*** Entering graphClipDraw(%d,%d,%d,%d)...\n",x1,y1,x2,y2) ;
//#endif
	
	// Total Virtual ACEDB screen area is arbitrarily set at
	// 30,000 x 30,000 logical units (pixels) squared
	// with origin (0,0)
	CRect rc(0,0,30000,30000) ;

	// Adjust clipping boundaries
	if (rc.left   < x1 )   rc.left = x1 ;
	if (rc.right  > x2 )  rc.right = x2 ;
	if (rc.top    < y1 )    rc.top = y1 ;
	if (rc.bottom > y2 ) rc.bottom = y2 ;

	// Need to convert from virtual graph (logical, virtual)
	// to CWnd client area (CDC Device, CView Client) coordinates?
	CClientDC dc(VIEWPTR) ;
	VIEWPTR->OnPrepareDC(&dc) ;
	dc.LPtoDP(&rc) ;
	rc.InflateRect( 1, 1 ) ;

	// Force the next WM_PAINT to draw this rectangle
	VIEWPTR->InvalidateRect(&rc, VIEWPTR->EraseGraph) ;
	VIEWPTR->EraseGraph = FALSE ;
}

extern "C" void graphBoxDraw (int k, int fcol, int bcol)
{
	if (!gActive->stack || !gDev || !VIEWPTR ) return ;

	Box box = gBoxGet (k) ;

	if (fcol >= 0)
		box->fcol = fcol ;
	if (bcol >= 0)
		box->bcol = bcol ;

	// Set up the box clipping rectangle
	// Total Virtual ACEDB screen area is arbitrarily set at
	// 30,000 x 30,000 logical units (pixels) squared
	// with origin (0,0)
	int x1 = uToXabs(box->x1)-1 ,
		y1 = uToYabs(box->y1)-1 ,
		x2 = uToXabs(box->x2)+1 ,
		y2 = uToYabs(box->y2)+1 ;
	CRect rc(0,0,30000,30000) ;

	// Adjust clipping boundaries
	if (rc.left   < x1 )   rc.left = x1 ;
	if (rc.right  > x2 )  rc.right = x2 ;
	if (rc.top    < y1 )    rc.top = y1 ;
	if (rc.bottom > y2 ) rc.bottom = y2 ;

	// Need to convert from virtual graph (logical, virtual)
	// to CWnd client area (CDC Device, CView Client) coordinates?
	CClientDC dc(VIEWPTR) ;
	VIEWPTR->OnPrepareDC(&dc) ;
	dc.LPtoDP(&rc) ;
	rc.InflateRect( 1, 1 ) ;

	// Force the next WM_PAINT to draw this rectangle
	VIEWPTR->InvalidateRect(&rc, VIEWPTR->EraseGraph) ;
	VIEWPTR->EraseGraph = FALSE ;

	/*
	// Force the next WM_PAINT to draw this rectangle
    if(VIEWPTR->EraseGraph)
        VIEWPTR->RedrawWindow(&rc,NULL,
            RDW_ERASE|RDW_INVALIDATE|RDW_UPDATENOW|RDW_NOCHILDREN ) ;
    else
        VIEWPTR->RedrawWindow(&rc,NULL,
            RDW_INVALIDATE|RDW_UPDATENOW|RDW_NOCHILDREN ) ;

	// VIEWPTR->InvalidateRect(&rc, VIEWPTR->EraseGraph) ;
	VIEWPTR->EraseGraph = FALSE ;
	*/
}

extern "C" void graphRedraw (void)
{
	if (!gActive || !gActive->stack || !gDev || !VIEWPTR )
		return ;

	gActive->isClear = FALSE ;
	VIEWPTR->Invalidate() ;
}

inline void drawXorLine( CPoint pt0, CPoint pt1 )
{
	unsigned char fcol = BLACK, bcol = BLACK ;

	int oldMode = CGraphView::gDC->SetROP2( R2_XORPEN ) ;

	// Note: setColors() returns the old fcol, bcol
	// in VIEWPTR->oldfcol and VIEWPTR->oldbcol ;
	VIEWPTR->setColors( fcol, bcol ) ;

	CGraphView::gDC->MoveTo( pt0 ) ;  // from point (x0,y0)
	CGraphView::gDC->LineTo( pt1 ) ;  // draw line to (x1,y1)
	
	CGraphView::gDC->SetROP2( oldMode ) ;

	// ...so I can reset them here
	VIEWPTR->setColors( VIEWPTR->oldfcol, VIEWPTR->oldbcol ) ;
}

extern "C" void graphXorLine (float x0, float y0, float x1, float y1 )
{
	if(!gDev || !VIEWPTR ) return ;

	CPoint pt0( uToXabs( x0 ), uToYabs( y0 ) ) ;
	CPoint pt1( uToXabs( x1 ), uToYabs( y1 ) ) ;

	if( CGraphView::gDC )
		drawXorLine(pt0,pt1) ;   // draw with valid CGraphView::gDC?

	else // get a CDC to draw in?
	{
		ASSERT_VALID(VIEWPTR) ;
	
		TRY
		{
			CClientDC pDC(VIEWPTR) ;
			CGraphView::gDC = &pDC ;

			drawXorLine( pt0, pt1 ) ;
			
			CGraphView::gDC = 0 ;

		}
		CATCH(CResourceException,e)
		{
			messcrash("graphXorLine() could not obtain a device context to draw in?") ;
		}
		END_CATCH
	}
}

inline void drawXorBox(Box box,float x,float y)
{
	int x1,y1,x2,y2 ;

	int oldMode = CGraphView::gDC->SetROP2( R2_XORPEN ) ;
	unsigned char fcol = box->fcol, bcol = box->bcol ;

	// Note: setColors() returns the old fcol, bcol
	// in VIEWPTR->oldfcol and VIEWPTR->oldbcol ;
	VIEWPTR->setColors( fcol, bcol ) ;
	
	x1 = uToXabs( x ) ;
	y1 = uToYabs( y ) ;
	x2 = uToXabs(x + (box->x2 - box->x1)) ;
	y2 = uToYabs(y + (box->y2 - box->y1)) ;

	CRect rect(x1, y1, x2, y2) ;
	VIEWPTR->drawFillRect (&rect) ;  // draw XOR Filled rectangle

	CGraphView::gDC->SetROP2( oldMode ) ;

	// ...so I can reset them here
	VIEWPTR->setColors( VIEWPTR->oldfcol, VIEWPTR->oldbcol ) ;
}

extern "C" void graphXorBox(int k, float x, float y)
{
	if(!gDev || !VIEWPTR ) return ;

	Box box = gBoxGet(k) ;

	if( CGraphView::gDC )
		drawXorBox(box,x,y) ;   // draw with valid CGraphView::gDC?

	else // get a CDC to draw in?
	{
		TRY
		{
			ASSERT_VALID(VIEWPTR) ;
	
			CClientDC pDC(VIEWPTR) ;
			CGraphView::gDC = &pDC ;

			drawXorBox(box,x,y) ;

			CGraphView::gDC = 0 ;
		}
		CATCH(CResourceException,e)
		{
			messcrash("graphXorLine() could not obtain a device context to draw in?") ;
		}
		END_CATCH
	}
}


/********************** Graph initialization ************************/

//****** Graphic Device initialization/cleanup/activation ***********/

extern "C" void graphDevActivate (int flag)
{
	if( !gActive || !gDev || !gSubDev) return ;

	CGraphWindow *pGraphWindow ;

	if( flag ) // activating Window?
	{

#if defined(OLD_WIN32_GRAPHS)
		pGraphWindow = DEVPTR ;
		pGraphWindow->MDIActivate() ;
#endif

#if defined(NEW_WIN32_GRAPHS)
			SUBDEVPTR( Activate, TRUE )
#endif
	}
	else // deactivating current devices and subdevices
	{
		// let ACEDB handle deactivation of the graph itself now
		// gActive = 0 ;
		gDev = 0 ;
		gSubDev = 0 ; 

#if defined(NEW_WIN32_GRAPHS)
		CGraphView::gViewPort = 0 ;
#endif
	}
}

/****************** 9/6/96:  Start of merged graphwin32.cpp code *************/

//****** ACEDB Specific Graphics Device Control Management Functions *****

UINT Ace2VKey( int aceKey )
{
	UINT virtKey;

	switch( aceKey )
	{
		// Mouse events?

		case LEFT_DRAG:
		case MIDDLE_DRAG:
		case RIGHT_DRAG:
			virtKey = 0 ;  // Drag events not handled yet
			break ;
		case LEFT_DOWN: 	// == 0
			virtKey = WM_LBUTTONDOWN ;
			break ;
		case LEFT_UP:	 	// == 2
			virtKey = WM_LBUTTONUP ;
			break ;
		case MIDDLE_DOWN:	// == 3
			virtKey = WM_MBUTTONDOWN ;
			break ;
		case MIDDLE_UP:		// == 5
			virtKey = WM_MBUTTONUP ;
			break ;
		case RIGHT_DOWN:	// == 6
			virtKey = WM_RBUTTONDOWN ;
			break ;
		case RIGHT_UP: 		// == 8
			virtKey = WM_RBUTTONUP ;
			break ;

		// Keyboard events?

		case LEFT_KEY:
			virtKey = VK_LEFT;
			break ;
		case RIGHT_KEY:
			virtKey = VK_RIGHT;
			break ;
		case UP_KEY:
			virtKey = VK_UP;
			break ;
		case DOWN_KEY:
			virtKey = VK_DOWN;
			break ;
		case PAGE_UP_KEY:
			virtKey = VK_PRIOR;
			break ;
		case PAGE_DOWN_KEY:
			virtKey = VK_NEXT;
			break ;
		case HOME_KEY:
			virtKey = VK_HOME;
			break ;
		case END_KEY:
			virtKey = VK_END;
			break ;
		case INSERT_KEY:
			virtKey = VK_INSERT;
			break ;
		case DELETE_KEY:
			virtKey = VK_DELETE;
			break ;
		case SPACE_KEY:
			virtKey = VK_SPACE;
			break ;
		case ESCAPE_KEY:
			virtKey = VK_ESCAPE;
			break ;
		case RETURN_KEY:
			virtKey = VK_RETURN;
			break ;
		// TAB_KEY is same as "PICK" event but I don't think
		// that this should be a problem in this context?
		case TAB_KEY:
			virtKey = VK_TAB;
			break ;

		// BACKSPACE_KEY cannot be used because
		// it is same as "RIGHT_UP"

		// case BACKSPACE_KEY:   
		// 	virtKey = VK_BACK;
		//	break ;

		default:
			if( aceKey > 0x20 && 	// > SPACE key
				aceKey < 127 )  	// < DEL key 
			// Then. just send whatever printable ASCII you were given?
			virtKey = aceKey ;

			else  // cannot translate?
				virtKey = 0 ;		
	}
	return virtKey ;
}

//**************************** graphxxxxBounds() functions *************************


extern "C" void graphMapBounds (float ux, float uy, float uw, float uh, float aspect)
{
	graphPlainBounds ( ux,  uy,  uw,  uh,  aspect) ; // temporary implementation?
}

void SetScrollBounds( CGraphView *pGraphView )
{
	// Need the current gActive for scroll values...
	if ( !gActive ) return ;
	
	// Need some valid CGraphView to specify scroll bounds for...
	if ( !pGraphView )
		if( !( gDev && VIEWPTR ) )
			return ;
		else
		{
			ASSERT_VALID(VIEWPTR) ;
			pGraphView = VIEWPTR ;
		}

	if ( !pGraphView->ScrollBoundsSet )
		pGraphView->ScrollBoundsSet = TRUE ;

	int dx, dy, xls, yls, xps, yps ;
	float cw, ch ;
	BOOL shrinkOnly = TRUE ;

	dx = dy = xls = yls = xps = yps = 0 ;
	pGraphView->graphSize.cx = gActive->w ;
	pGraphView->graphSize.cy = gActive->h ;

	switch(gActive->type)
	{
		case PLAIN:
			pGraphView->SetScaleToFitSize( pGraphView->graphSize ) ;
			return ;

		case TEXT_FIT:
		case PIXEL_FIT:
			shrinkOnly = FALSE ;
			break ;

		case MAP_SCROLL:
		case PIXEL_SCROLL:
			dx = dy = 1 ;
			xls = yls = 10 ;
			xps = yps = 100 ;
			break ;
		case TEXT_FULL_SCROLL:
		case TEXT_FULL_EDIT:
			xps = 10 ; xls = 1; // enable horizontal scrolling allowed
			// deliberate fall through... others are yls = yps = 0
		case TEXT_SCROLL:
			yps= 10 ;  yls = 1; // enable vertical scrolling allowed
			graphTextInfo(&dx,&dy,&cw,&ch) ;
			break ;

		default:
			ASSERT(FALSE) ;
	}
	pGraphView->pageSize.cx = dx*xps ;
	pGraphView->pageSize.cy = dy*yps ;
	pGraphView->lineSize.cx = dx*xls ;
	pGraphView->lineSize.cy = dy*yls ;
	pGraphView->SetScrollSizes(MM_TEXT, pGraphView->graphSize,
								pGraphView->pageSize, pGraphView->lineSize ) ;
	if(!IsSubGraph(gActive))
		pGraphView->ResizeParentToFit(shrinkOnly) ;
}

extern "C" void graphTextBounds (int nx, int ny)
{
	if (gActive->type != TEXT_SCROLL && gActive->type != TEXT_FULL_SCROLL)
		messcrash ("textBounds called on invalid graph type %d",gActive->type) ;

	/* reset graph device dimensions */
	gActive->ux = 0 ;			
	gActive->uy = 0 ;
	gActive->uw = nx ;			
	gActive->uh = ny ;
	gActive->w = uToXrel( nx ) ;
	gActive->h = uToYrel( ny ) ;

	gUpdateBox0() ;

	// Forces the bounds change to be reset in the CGraph?
	SetScrollBounds() ;
}

extern "C" float graphFakeBounds (float ny)
// SRK added this to make it compile, may need to call window system
{ float old = gActive->uh ;
  gActive->uh = ny ;
  gActive->h = ny * gActive->yFac ;

  return old ;
}

extern "C" void graphPixelBounds (int nx, int ny)
{
 	// reset graph device dimensions
 	// # of pixels, device coordinates?
	gActive->ux = 0 ;			
	gActive->uy = 0 ;
	gActive->uw = gActive->w = nx ;
	gActive->uh = gActive->h = ny ;
	
	gUpdateBox0() ;

	// Forces the bounds change to be reset in the CGraph?
	SetScrollBounds() ;
}


/******* process all outstanding events in the queue *******/

extern "C" void graphProcessEvents (void)
{
	// WIN32 handles this automatically?
}

//////////////////////////////////////////////////////////////////////////////////////

// If the graph is scrollable, graphGoto()
// tries to center x,y in visible region 
extern "C" void graphGoto (float x, float y)
{
	if (!gActive || !gDev || 
		!VIEWPTR || !(VIEWPTR->ScrollBoundsSet) ) return ;

	switch(gActive->type)
	{
		case PLAIN:
		case TEXT_FIT:
		case PIXEL_FIT:
			TRACE1("Graph #%d is not scrollable in graphGoto()?\n", gActive->id ) ;
			break ;

		case TEXT_SCROLL:
		case TEXT_FULL_SCROLL:
		case PIXEL_SCROLL:
		case MAP_SCROLL:
		{
			CPoint centre( uToXabs(x),uToXabs(y) ) ;
			VIEWPTR->ScrollToPosition(centre) ;
		}
			break ;

		default:
			ASSERT(FALSE) ;
	}
}

extern "C" void graphWhere (float *x1, float *y1, float *x2, float *y2)
{
	CPoint pt = VIEWPTR->GetScrollPosition() ;
	CSize  sz = VIEWPTR->GetTotalSize() ;
	// Beware of NULL pointers!
	if(x1) *x1 = XtoUabs( pt.x );
	if(y1) *y1 = YtoUabs( pt.y ) ;
	if(x2) *x2 = XtoUabs( pt.x + sz.cx ) ;
	if(y2) *y2 = YtoUabs( pt.y + sz.cy ) ;
}

static void sendMouseEvent( UINT MsgID, float x, float y )
{
	if( !gDev || !VIEWPTR ) return ;
	
	static LPARAM lParam ;
	CPoint pt( (int)uToXabs(x),(int)uToYabs(y) ) ;

	CClientDC dc(VIEWPTR) ;
	VIEWPTR->OnPrepareDC(&dc) ;
	dc.LPtoDP(&pt) ;

	lParam = (LPARAM)MAKELONG( (WORD)pt.x, (WORD)pt.y ) ;
	TRACE3("graphEvent() mouse action# == %d, Device Position = (%d,%d)\n",
			MsgID, LOWORD(lParam), HIWORD(lParam) ) ;
	theApp.m_pMainWnd->PostMessage(MsgID, 0, lParam) ;  // no fwkeys == wparam set
	return ;
}

static void sendChar(UINT virtKey) // for printable ASCII characters
{
	static WPARAM wParam ;
	wParam = (WPARAM)virtKey ;			
	UINT scanCode = MapVirtualKey(virtKey,0) ;  // May not work for some keys?
	static LPARAM lParam ;
	lParam = (scanCode << 16) & 0x0F00 ;

	TRACE2("graphEvent() virtual character key sent == %c, Scan Code = %c\n",
			 virtKey, LOBYTE(HIWORD(lParam)) ) ;
	theApp.m_pMainWnd->PostMessage(WM_CHAR, wParam, lParam) ;
}

static void sendKey(UINT virtKey)
{
	static WPARAM wParam ;
	wParam = (WPARAM)virtKey ;			
	UINT scanCode = MapVirtualKey(virtKey,0) ;  // May not work for some keys?
	static LPARAM lParam ;
	lParam = (scanCode << 16) & 0x0F00 ;

	TRACE2("graphEvent() virtual character key sent == %c, Scan Code = %c\n",
			 virtKey, LOBYTE(HIWORD(lParam)) ) ;
	theApp.m_pMainWnd->PostMessage(WM_KEYDOWN, wParam, lParam) ;
	theApp.m_pMainWnd->PostMessage(WM_KEYUP, wParam, lParam) ;
}

extern "C" void graphEvent (int action, float x, float y)
/* sends an event over the "wire" to active window
   currently can only handle key events defined in keys.h, 
   printable ascii events and mouse events */
{
	UINT aceKey ;
	
	if( !gDev || !VIEWPTR ) return ;
		return ;

	if(action < 0 || !(aceKey = Ace2VKey(action)) )
	{
		messcrash ("Invalid action number specified to graphEvent()\n");
		return ;
	}

	// ASSERT(valid action)
	// Set the global graphEvent (user) coordinates
	graphEventX = x ;
	graphEventY = y ;

	if (action <= RIGHT_UP) /* means a mouse event; aceKey is the WIN32 MsgID */
		// ASSERT( action == { 0,2,3,5,6,8 } )
		sendMouseEvent( aceKey, x, y ) ;
	
	else if( action <= 127 ) /* keyboard events */
	 	// ASSERT( 8 (== RIGHT_UP) < action <= 127 AND printable )
		sendChar( aceKey ) ;		// Post a WIN32 WM_CHAR event
	else
		sendKey( aceKey ) ;			// send WM_KEYDOWN/KEYUP events for keys defined in key.h
}

extern "C" void graphRetitle (char *title)
{
	char graphTitle[80] ;
	sprintf( graphTitle, "[%u]:",gActive->id ) ;
	strcat( graphTitle, title ) ; 
	GRAPHPTR->SetTitle( graphTitle ) ; // reset window titlebar to graph->id:title only
}


extern "C" void graphPop (void)
{
	if( !gDev || !VIEWPTR ) return ;

#if !defined(NEW_WIN32_GRAPHS)
	CGraphWindow *pGraphWnd = (CGraphWindow *)(VIEWPTR->GetParentFrame()) ;
	pGraphWnd->MDIActivate() ;
	if(pGraphWnd->IsIconic())
		pGraphWnd->MDIRestore() ;
#else
	SUBDEVPTR( Activate, TRUE )
#endif
}

void CGraphView::OpenGraphMessage( const CString& text )
{
	// ignore additional calls to OpenGraphMessage()
	// generated by redundant graphMessage() calls
	if(	m_pMessageWnd ) return ;

	CString msgTitle = GetGraph()->name ;
	msgTitle += " Message" ;

	m_pMessageWnd = new CMsgWindow( msgTitle, text, this) ;

}

void CGraphView::CloseGraphMessage( BOOL flg )  // flg is a dummy parameter)
{
	// block recursive calls: window already destroyed?
	if( !m_pMessageWnd ) return ;

	// for blocking recursive calls...
	CMsgWindow *pMW = (CMsgWindow *)m_pMessageWnd ;
	m_pMessageWnd = 0 ;
	
	// Note: func[MESSAGE_DESTROY]() sometimes recursively calls
	//       CloseGraphMessage() (via graphUnMessage() ) hence the block above
	if( gActive->func[MESSAGE_DESTROY] )
		(gActive->func[MESSAGE_DESTROY])() ;
		
	// Send this message to the base class of CMsgWindow
	pMW->DestroyWindow() ;
}
 
extern "C" void graphMessage (char *text)
{
#if !defined(NEW_WIN32_GRAPHS)
	DEVPTR->OpenGraphMessage( text ) ;
#else
	SUBDEVPTR( OpenGraphMessage, text )
#endif
}

extern "C" void graphUnMessage (void) 
{
#if !defined(NEW_WIN32_GRAPHS)
	DEVPTR->CloseGraphMessage( ) ;
#else
	SUBDEVPTR( CloseGraphMessage, TRUE )
#endif
}

/*****************************************************************/

/* callbacks, associators and all that */

float graphEventX, graphEventY ;
	/* for users to get at event position, in user coords */

// This routine monitors graph interrupts (a.k.a. ESCAPE )
// set by the accelerator key monitoring in CGraphWindow
// This interrupt is distinct from the "System Break" key
// which has more severe (messcrash()) consequences...
BOOL graphInterruptCalled(void)
{ 
	MSG msgCur;

	if( ::PeekMessage(&msgCur, NULL, WM_KEYFIRST, WM_KEYLAST, PM_REMOVE) )
	{
		TRACE("\n*** graphInterruptCalled( message = %.4x, key == %.4x)\n",
				msgCur.message,	msgCur.wParam  ) ;
		if( msgCur.message == WM_KEYDOWN &&
			msgCur.wParam == VK_ESCAPE )
		{
			do ::PeekMessage(&msgCur, NULL, WM_KEYFIRST, WM_KEYLAST, PM_REMOVE) ;
			while(!(msgCur.message == WM_KEYUP && msgCur.wParam == VK_ESCAPE) ) ; 
			return TRUE ;
		}
	}
	return FALSE ;
}

void graphScreenSize (float *sx, float *sy, float *fx, float *fy, int *px, int *py)
{
  int fw, fh ;  /* Font width and height */
  gFontInfo (0, &fw, &fh) ;

  /* gScreenX and gScreenY are determined at runtime during initialization
     and are normalized to SCREENX_ASPECT and SCREENY_ASPECT
	 device independent screen size */ 
  if (sx) *sx = SCREENX_ASPECT ;  
  if (sy) *sy = SCREENY_ASPECT ;

  if (fx) *fx = (float)CGraphView::gScreenX/fw ;
  if (fy) *fy = (float)CGraphView::gScreenY/fh ;

  if (px) *px = CGraphView::gScreenX ;
  if (py) *py = CGraphView::gScreenY ;
}

////////////////////////////////////////////////////////////////////////////
// ACEDB Cut/Copy/Paste functionality

static CString postBuffer, pasteBuffer ;
static HGLOBAL pasteBufferHandle = NULL, copyBufferHandle = NULL ;

/* write to internal ACEDB screen cut/paste buffer */
extern "C" void  graphPostBuffer (char *text)
{
	postBuffer.Empty() ;
	postBuffer = text ;
}

extern CString &dos2unix(const CString &source) ;
extern CString &unix2dos(const CString &source) ;

/* returns Clipboard cut/paste buffer CF_TEXT contents */
extern "C" char *graphPasteBuffer (void)
{
	union jester { const char *yin ; char *yang ; } fool ;
	char *pText ;

	pasteBuffer.Empty() ;

	// Open clipboard
	if(!OpenClipboard(NULL))
	{
		graphOut("graphPasteBuffer(): Could not open Clipboard?") ;
	}
	else // read clipboard
	{
		// Get the Global memory object with the clipboard data
		if( ((pasteBufferHandle = GetClipboardData(CF_TEXT) ) == NULL)
			|| ((pText = (char *)GlobalLock(pasteBufferHandle)) == NULL) )
			graphOut("graphPasteBuffer(): Could not access Clipboard data?") ;
		else // Valid text data found in clipboard
		{
			// Translate clipboard data from
			// DOS format into UNIX/ACEDB format
			pasteBuffer = dos2unix(pText) ;

			// Release the global lock
			GlobalUnlock(pasteBufferHandle) ;
		}
	}

	// return the pastebuffer (as an ordinary char* pointer)
	fool.yin = (LPCTSTR)pasteBuffer ;
	return fool.yang ;
}

// User initiated copying of of internal ACEDB screen cut/paste buffer data
// buffer data buffer dataonto the Windows clipboard (for subsequent graphPaste)
void graphCopy()
{
	CString copyBuffer ;

	// Translate ACEDB postBuffer data from UNIX/ACEDB to DOS format
	copyBuffer = unix2dos(postBuffer) ;

	// Open clipboard
	if(!OpenClipboard(NULL))
		graphOut("Sorrry: Could not open Clipboard?") ;
	else // Write to clipboard
	{
		// Allocate a moveable global memory object
		if(copyBufferHandle == NULL)
			copyBufferHandle =
				GlobalAlloc(GMEM_MOVEABLE | GMEM_DDESHARE | GMEM_ZEROINIT,
							copyBuffer.GetLength()+1) ;
		else			
			copyBufferHandle =
				GlobalReAlloc( copyBufferHandle,
							   GMEM_MOVEABLE | GMEM_ZEROINIT,
							   copyBuffer.GetLength()+1) ;
		if( copyBufferHandle == NULL)
			graphOut("graphCopy() memory handle error: Cannot copy data to Clipboard?") ;
		else
		{
			// Global lock it
			char *pText ;
			if((pText = (char *)GlobalLock(copyBufferHandle)) == NULL )
				graphOut("graphCopy() memory lock error: Cannot copy data to Clipboard?") ;
			else
			{
				// Copy text into the global memory object
				strcpy( pText, (LPCTSTR)copyBuffer ) ;
				// Global unlock it
				GlobalUnlock( copyBufferHandle ) ; // do I care if this fails?
				// Write it to the Clipboard
				SetClipboardData(CF_TEXT,copyBufferHandle) ;
			}
		}
	}
}

// User initiated pasting of clipboard contents into ACEDB
void graphPaste()
{
	// Gets Clipboard CF_TEXT data
	char *text = graphPasteBuffer() ;
	//... then sends it to ACEDB as graphEvent()
	while(*text) graphEvent(*text++,0,0) ;
}

void CGraphView::OnEditCopy() 
{
	graphCopy() ;
}

void CGraphView::OnEditPaste() 
{
	graphPaste() ;
}

/********* graphWaitCursor() management ************/

void CGraphView::setMyCursorOn()
{
	// remember the current cursor for restoration later
	m_otherCursor = SetCursor(m_myCursor) ;
	m_myCursorOn = TRUE ;
}

void CGraphView::setMyCursorOff()
{
	// reset to the old cursor, if any exists
	if( m_otherCursor != NULL )
	{
		SetCursor (m_otherCursor) ;
		m_otherCursor = NULL ; // then, forget it until next time
	}
	m_myCursorOn = FALSE ;
}

UINT CGraphView::OnNcHitTest(CPoint point) 
{
	UINT hitType = CScrollView::OnNcHitTest(point) ;
	if(hitType == HTCLIENT)
	{
#if defined(REALLY_VERBOSE)
		TRACE("OnNcHitTest(HTCLIENT): My cursor is %s the client window\n",
				myCursorOn()?"within":"entering") ;
#endif
		if( !myCursorOn() ) setMyCursorOn() ; 
	}

	return hitType ;
}

void graphWaitCursor (BOOL on)
{
  if(!gActive || !gDev || !VIEWPTR) return ;

  // Just select the current cursor for display
  // but let the current viewport control its display
  if (on)
	  VIEWPTR->setMyCursor (CGraphView::Arrow) ;
  else
	  VIEWPTR->setMyCursor (CGraphView::HourGlass) ;
}

/********** End of File ********/
 

 
 
 
