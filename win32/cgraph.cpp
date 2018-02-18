/*  CGraph.cpp : implementation file
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: ACE Graph associated SubDev == "CScrollView" window implementation
 * HISTORY:
 * Last edited: Jun 10 10:15 1996 (rbrusk):
 *		-	implemented graphWaitCursor() (the XWindows way; may not work well in WIN32?)
 * Created: Jul 21 18:40 1995 (rbrusk)
 *-------------------------------------------------------------------
 */
 /* $Id: cgraph.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $ */

#include "stdafx.h"
#include "win32.h"
#include "winace.h"
#include "cgraph.h"
#include "graphwin32_.h" 

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CGraph

IMPLEMENT_DYNCREATE(CGraph, CDocument)

#define new DEBUG_NEW

//////////////////////////////////////////////////////////////////////////
// Initialize CGraph class static variables


////////////////////////////////////////////////////////////////////////////
// CGraph Constructor/Destructor
 
CGraph::CGraph()
{
#if !defined(NEW_WIN32_GRAPHS)
	m_pACEDB_Graph = 0 ;
#endif
}

CGraph::~CGraph()
{
}


BEGIN_MESSAGE_MAP(CGraph, CDocument)
	//{{AFX_MSG_MAP(CGraph)
	ON_COMMAND(ID_FILE_SEND_MAIL, OnFileSendMail)
	ON_UPDATE_COMMAND_UI(ID_FILE_SEND_MAIL, OnUpdateFileSendMail)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()



/////////////////////////////////////////////////////////////////////////////
// CGraph diagnostics

#ifdef _DEBUG
void CGraph::AssertValid() const
{
	CDocument::AssertValid();
}

void CGraph::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);

	dc	<< "CGraph specific data members: \n"

#if !defined(NEW_WIN32_GRAPHS)
		<< "\tGraph_ ptr == " <<  m_pACEDB_Graph 
#endif

		<< "\n\tCGraph::Frozen semaphore: " << (Frozen ? "True" : "False")
		<< "\n\t***** End of CGraph Dump? *****\n\n" ;

}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CGraph serialization

void CGraph::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// CGraph commands

void CGraph::InitGraph()
{
}

BOOL CGraph::OnOpenDocument(LPCTSTR lpszPathName) 
{
	// if (!CDocument::OnOpenDocument(lpszPathName))
	//	return FALSE;
	
	InitGraph();	
	return TRUE;
}

BOOL CGraph::OnNewDocument()
{
	// if (!CDocument::OnNewDocument())
	//	 return FALSE;
	InitGraph();
	return TRUE;
}

// This flag postpones deactivation of a dying gActive and associated
// gDev, gStk, gBox variables during a graphDestroy() triggered
// OnCloseDocument() call;  note: this mechanism is not reentrant
BOOL FROZEN = FALSE ;

void CGraph::OnCloseDocument() 
{
	FROZEN = TRUE ;
	CDocument::OnCloseDocument();
	FROZEN = FALSE ;
}

/******* graphloop management data structures ***************************/
static Graph_ blockGraph = 0 ;
static int loopLevel = 0 ;

typedef struct LoopInfoStruct {
  Graph_ g;
  BOOL isBlocking ;
  BOOL isDone ;
  int retval ;
} LOOPINFO ;

// static Associator assw2g = 0 ;
static Array loopInfo = 0 ;

//****** Graphic Device Initialization/Finalization ******************************

extern void graphMenuInit() ;
extern void graphMenuFinish() ;

extern "C" void graphDevInit (int *argcptr, char **argv)
{
	// assw2g = assCreate();
	loopInfo = arrayCreate (4, LOOPINFO) ;

	// I should probably initialize ACEDB CGraph Windows
	// parameters and resources here, such as the ACEDB
	// colour palette, brushes, pens etc.?

	VERIFY( CGraphView::Arrow = theApp.LoadStandardCursor (IDC_ARROW) );
	VERIFY( CGraphView::HourGlass = theApp.LoadStandardCursor (IDC_WAIT) );

	graphMenuInit() ;
	CGraphView::graphColorInit() ;
	CGraphView::graphFontInit() ;
}

extern "C" void graphDevStart (Graph_ g)
{
	// Don't need to do anything, since the graphics event loop to
	// capture user events for each CGraph-associated CView is
	// automatically setup by the Windows framework?
}


extern "C" void graphDevFinish(void)
{
	// Don't need to do too much here since most of the 
	// graphics system shutdown is taken care of by WinApp
	// Just release resources allocated in graphDevInit()?

	graphMenuFinish() ;
	CGraphView::graphColorFinish() ;
	CGraphView::graphFontFinish() ;
}

/******* graphloop management********************************************/
/******* can be recursive and either blocking or non-blocking ***********/
extern "C" BOOL graphIsBlocked()
{
	if (blockGraph && !(gActive->notToBeBlocked) && gActive != blockGraph)
	  return TRUE;
	else
	  return FALSE ;
}

int graphLoop (BOOL isBlock)
{
  LOOPINFO *info ;
  int retval ;
#ifdef ACEDB
  extern mainActivity (char*) ;
#endif

  if (!gDev)
    return 0 ;

  info = arrayp(loopInfo, ++loopLevel, LOOPINFO) ;
  info->g = gActive ;
  info->isBlocking = isBlock ;
  info->retval = 0 ;
  info->isDone = FALSE ;

  if (isBlock)
    blockGraph = gActive ;

  // while (!info->isDone) 
  retval = theApp.Run() ;

  info = arrp(loopInfo, --loopLevel, LOOPINFO) ;
  if (info->isBlocking)
    blockGraph = info->g ;
  else
    blockGraph = 0 ;

  return retval ;
}

#ifdef _DEBUG
static void loopStatus (void)
{
  int i ;
  LOOPINFO *loop = arrp(loopInfo, 0, LOOPINFO) ;

  for (i = 1, ++loop ; i <= loopLevel ; ++i, ++loop)
    printf (" level %d, isDone %d, isBlocking %d, ret %d, graph %x\n",
	    i, loop->isDone, loop->isBlocking, loop->retval, loop->g) ;
}
#endif

BOOL graphLoopReturn (int retval)
{
  int i = loopLevel ;

#ifdef _DEBUG
  printf ("LoopReturn %d on graph %x\n", retval, gActive) ;
  loopStatus() ;
#endif

  while (i)
    if (arrp(loopInfo, i, LOOPINFO)->g == gActive &&
	!arrp(loopInfo, i, LOOPINFO)->isDone)
      break ;
    else
      --i ;
  if (i)
    { arrp(loopInfo, i, LOOPINFO)->retval = retval ;
      while (i <= loopLevel) arrp(loopInfo, i++, LOOPINFO)->isDone = TRUE ;
#ifdef _DEBUG
      printf ("Return TRUE\n") ; loopStatus() ;
#endif
      return TRUE ;
    }
#ifdef _DEBUG
  printf ("Return FALSE\n") ; loopStatus() ;
#endif
  return FALSE ;
}

BOOL InGraphLoop()
{
	return( loopLevel ? TRUE : FALSE ) ;
}

// returns current graphLoop() return value
// (assumes that graphLoop() itself is terminated subsequently)
// Returns error "-1" if no unterminated graphLoop() is running
int ExitGraphLoop()
{
	return( !loopLevel ? -1 :arrp(loopInfo, loopLevel, LOOPINFO)->retval ) ;
}

// Returns "isDone" status of current graphLoop()
// Assumes that if no graphLoop is running, that top level
// WIN32 event loop is terminating, hence it returns "TRUE"
BOOL LoopIsDone()
{
	return( arrp(loopInfo, loopLevel, LOOPINFO)->isDone ) ;
}

/********************OLD IMPLEMENTATION (Pre-June 11/96) ************
   Within Win32, graphDevCreate() and graphDevDestroy() are the key routines
   for forging the link between the ACEDB graphics world and the Windows Framework
   world of CDocuments, CMDI{Frame/Child}Windows and CViews et al.
   
   The "Document" is a Windows C++ encapsulation of the "Graph_" construct,
   of a type specific to the category of graph.  The view is one of a set of
   views based upon graph type (and to some extent, what the code knows about
   the dataset creating the graph display, i.e. is_a main window or a tree display
   => I can perhaps used special Windows control tools to format the view???
   
   The contents of Dev in the Graph_ graph structure is a pointer (cast to void) to
   the C++ document encapsulating the Graph_ structure and which points to all
   the relevent graphical entities (Views, Doc template, et al) for the C++
   implementation of the graph display window.

 ************** NEW IMPLEMENTATION (POST-Sept 1/96) *************

	The objective of the new implementation of the underlying WIN32 ACEDB graphics
	device is both to rationalize the graphics code and to properly implement subgraphs.

	MFC Document/Frame/View Framework paradigm of the old implementation is the foundationm of this
	new implementation, with the ACEDB "graph->dev" dereferencing the CGraph (document) object
	for which  a child frame window is automatically created.

	The global DEVPTR complements GRAPHPTR and points to the active Dev (CMDIChildWindow) window.
	The global SUBDEVPTR points to the currently active Subdev.    
	The fundamental difference between Dev's and Subdev's is that all Subdev's are
	child CWnd's of a parent Dev window.

	However, this new implementation differs from the original implementation in that
	SubDev refers to either the MFC CMDIChildWindow of the ACEDB graph object or the
	CGraphView window of a ACEDB subgraph object.

 ****************************************************************************/


#if defined(NEW_WIN32_GRAPHS) 

#include "PlainView.h"
#include "FullScrollView.h"
#include "TextScrollView.h"
#include "Fitview.h"
#include "TextFitView.h"
#include "MapScrollView.h"
#include "PixelScrollView.h"
#include "PixelFitView.h"
// #include "FullEditView.h"

static CCreateContext *GetGraphContext(int type)
{
	static CCreateContext pCC ;
	pCC.m_pCurrentDoc		= NULL ;
	pCC.m_pNewDocTemplate	= NULL ;
	pCC.m_pLastView			= NULL ;
	pCC.m_pCurrentFrame		= NULL ;
	switch(type)
	{
		case PIXEL_FIT:
			pCC.m_pNewViewClass =
				RUNTIME_CLASS (CPixelFitGraphView) ;
			break ;  
		case TEXT_FIT:
			pCC.m_pNewViewClass =
				RUNTIME_CLASS (CTextFitGraphView) ;
			break ;  
		case TEXT_SCROLL:
			pCC.m_pNewViewClass =
				RUNTIME_CLASS (CTextScrollGraphView) ;
			break ;  
		case TEXT_FULL_SCROLL:
			pCC.m_pNewViewClass =
				RUNTIME_CLASS (CTextFullScrollGraphView) ;
			break ;  
/*		case TEXT_FULL_EDIT: // Not yet implemented?
			pCC.m_pNewViewClass =
				RUNTIME_CLASS (CTextFullEditGraphView) ;
			break ; */  
		case PLAIN:
			pCC.m_pNewViewClass =
				RUNTIME_CLASS (CPlainGraphView) ;
			break ;  
		case MAP_SCROLL:
			pCC.m_pNewViewClass =
				RUNTIME_CLASS (CMapScrollGraphView) ;
			break ;  
		case PIXEL_SCROLL:
			pCC.m_pNewViewClass =
				RUNTIME_CLASS (CPixelScrollGraphView) ;
			break ;  
		default: 
			messcrash ("Invalid graph type %d in GetGraphContext()", type) ;
	}
	 return &pCC ;
}

#endif // defined(NEW_WIN32_GRAPHS) 

// Note: coordinates are device independent x(0.0..1.3), y(0.0..1.0)
extern "C" void graphDevCreate (Graph_ graph, float x, float y, float w, float h)
{
	ASSERT( graph != NULL ) ;
	TRACE("\n***graphDevCreate( graph[%d]== %s | Type: %d )\n",
			graph->id,graph->name,graph->type ) ;
	TRACE("\t\tOrigin(%4.2f, %4.2f), Size(%4.2f, %4.2f)\n",x,y,w,h) ;

	ASSERT( graph->dev == NULL) ;
	ASSERT( graph->type >= 0 && graph->type < NUMGRAPHTYPES ) ;

	if(graph->type == TEXT_FULL_EDIT )
	{
		messout("The TEXT_FULL_EDIT graph type is not implemented yet!") ;
		return ;
	}
	
	// Set the graph as the active one
	graph->subdev = graph->dev = 0 ; // not opened yet!
  	graphActivate (graph->id) ;

	// Set Window Relative Graph (Child Window) Position and Size
	// in Windows CDC display units, within MainFrame window

	int wx	= (int) ( (x/SCREENX_ASPECT)*(float)CGraphView::gScreenX ) ; 
	int wy	= (int) ( (y/SCREENY_ASPECT)*(float)CGraphView::gScreenY ) ;
	int ww	= (int) ( (w/SCREENX_ASPECT)*(float)CGraphView::gScreenX ) ;
	int wh	= (int) ( (h/SCREENY_ASPECT)*(float)CGraphView::gScreenY ) ;

	CGraphWindow::InitDevStruct(graph, wx, wy, ww, wh) ;

	// Construct a CGraph document pseudonym for the CDocTemplate 
	char docPseudonym[MAXPATHLEN] ;
	char graphTitle[80] ;
	sprintf( graphTitle, "[%u]:",graph->id ) ;
	strcat( graphTitle, graph->name ) ; 

	strcpy( docPseudonym, graphTitle ) ;
	strcat( docPseudonym, AceGraphDesc[graph->type].FNExt ) ;

	// then, open the device
	CGraph *pGraphDev ;

	// MUST use "VERIFY" rather than "ASSERT" here to ensure that
	// the expression is evaluated in the Release version of the program!!! (rmb) Jan 24/96

#if defined(NEW_WIN32_GRAPHS)
	FROZEN = TRUE ; // To postpone window activation event until after creation?
#endif

	VERIFY( (pGraphDev = (CGraph *)theApp.OpenDocumentFile( docPseudonym ) ) != NULL ) ;

#if defined(NEW_WIN32_GRAPHS)
	FROZEN = FALSE ; 
#endif

	pGraphDev->SetTitle( graphTitle ) ; // reset window titlebar to "id:graphName" only

	// Set the current ACEDB device and subdevice to the graph child frame window
	// alias CGraph document & its associated first view window
	gDev = graph->dev = (void *)pGraphDev ;

#if defined(OLD_WIN32_GRAPHS)
	gSubDev = graph->subdev =(void *)(pGraphDev->GetGraphView()) ;
	pGraphDev->SetGraph(graph) ; // current graph hook in WIN32 window object
#else
	gSubDev = graph->subdev =(void *)pGraphDev->GetGraphWindow() ;	// Gets the new frame window?

	// Gets the new document view? Sets Viewptr
	SETVIEWPORT
#endif
}  

extern "C" void graphDevDestroy (void)
{
	  // exit all unterminated graphLoops() for this graph
	  while (graphLoopReturn (0)) ;

	// Note: Another previously active window is automatically
	// reactivated by OnCloseDocument(), hence
	// gActive and gDev reset to that other CGraph==Graph_
	if( !gActive || !gActive->dev ) return ;

#if defined(OLD_WIN32_GRAPHS)
		TRACE1("Entering graphDevDestroy() with Graph id == %d \n", GRAPHPTR->GetGraph()->id ) ; 
#else
		Graph_ theGraph ;
		USESUBDEVPTR(theGraph, GetGraph, TRUE )
		TRACE1("Entering graphDevDestroy() with Graph id == %d \n", theGraph->id ) ; 
#endif
		GRAPHPTR->OnCloseDocument() ;
		// gDev is not valid at this point, so set it to NULL
		// but gActive is still set to gDev's former owner ACEDB graph?
		gActive->dev = gDev = 0 ;
		gActive->subdev = gSubDev = 0 ;

#if defined(NEW_WIN32_GRAPHS)
		CGraphView::gViewPort = 0 ;		
#endif

}

static CGraphView *GetNewSubGraphView(int type)
{
	CGraphView *pGV ;
	switch(type)
	{
		case PIXEL_FIT:
			pGV = new CPixelFitGraphView ;
			break ;  
		case TEXT_FIT:
			pGV = new CTextFitGraphView ;
			break ;  
		case TEXT_SCROLL:
			pGV = new CTextScrollGraphView ;
			break ;  
		case TEXT_FULL_SCROLL:
			pGV = new CTextFullScrollGraphView ;
			break ;  
/*		case TEXT_FULL_EDIT: // Not yet implemented?
			pGV = new CTextFullEditGraphView ;
			break ; */  
		case PLAIN:
			pGV = new CPlainGraphView ;
			break ;  
		case MAP_SCROLL:
			pGV = new CMapScrollGraphView ;
			break ;  
		case PIXEL_SCROLL:
			pGV = new CPixelScrollGraphView ;
			break ;  
		default: 
			messcrash ("Invalid graph type %d in GetNewSubGraphView()", type) ;
	}
	 return pGV ;
}

// note: coordinates are window/display relative pixel dimensions, not device independent
extern "C" void graphSubDevCreate (Graph_ graph, float x, float y, float w, float h)
{
#if !defined(NEW_WIN32_GRAPHS)
	TRACE("\n\t***** Executing graphSubDevCreate() stub for graph->id == %d *****\n", graph->id) ;
#else // defined(NEW_WIN32_GRAPHS)
	ASSERT( graph != NULL ) ;
	TRACE("\n***graphSubDevCreate( graph[%d]== %s | Type: %d )\n",
			graph->id,graph->name,graph->type ) ;
	TRACE("\t\tOrigin(%4.2f, %4.2f), Size(%4.2f, %4.2f)\n",x,y,w,h) ;

	ASSERT( graph->dev != NULL) ; // contrary to graphs...
	ASSERT( graph->type >= 0 && graph->type < NUMGRAPHTYPES ) ;
	
	// Set the graph as the active one
	CWnd *pParent = (CWnd *)(graph->parent->subdev) ;
	graph->subdev = 0 ; // not opened yet!
  	graphActivate (graph->id) ; // but set as gActive

	CSubGraphWindow::InitSubDevStruct(graph, (int)x, (int)y, (int)w, (int)h) ;

	// Creates a subgraph window device
	CSubGraphWindow *pGraphSubDev = GetNewSubGraphView(graph->type) ;
	gSubDev = graph->subdev = (void *)pGraphSubDev ;
	SetScrollBounds(pGraphSubDev ) ;

	// Does the CCreateContext trigger CGraphView creation as a child window?

	graph->isBlocked = TRUE ; // Need to set this to block an indirect infinite recursion loop

	pGraphSubDev->Create (  NULL, NULL, WS_CHILD|WS_VISIBLE,
							pGraphSubDev->GetGraphRect(),
							pParent /*Parent CWnd */,
						    graph->id /* reused as the Child ID */,  NULL ) ;
	graph->isBlocked = FALSE ;

	SETVIEWPORT
	ASSERT( VIEWPTR ) ;  // Is a CGraphView created by this point?

#endif // defined(NEW_WIN32_GRAPHS)
}

extern "C" void graphSubDevDestroy (void) 
{
#if !defined(NEW_WIN32_GRAPHS)
	TRACE("\n\t***** Executing graphSubDevDestroy() stub *****\n") ;
#else // defined(NEW_WIN32_GRAPHS)
	if( gActive && gActive->dev &&
		gActive->subdev &&
		IsSubGraph( gActive ) )  // only destroy true subgraphs here... 
	{
		Graph_ theGraph = gActive ;
		CSubGraphWindow *theSubGraph = (CSubGraphWindow *)gActive->subdev ;
		theSubGraph->DestroyWindow() ;
		gActive = theGraph ; // restore gActive, just in case
		gActive->dev = 0 ;
		gActive->subdev = 0 ;  // set the subgraph graph to 0?
		gSubDev = gActive->parent->subdev ;  // reset to parent's subdev
		SUBDEVPTR( SetActiveGraph, TRUE )
	}
#endif // defined(NEW_WIN32_GRAPHS)
}

 
 
