/*  CGraph.h : header file
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: ACE Graph associated Dev == "CDocument" header file
 * HISTORY:
 * Last edited: Apr 30 02:32 1996 (rbrusk): graphBlocked() for blocked loops
 * * Mar  23 14:40 1996 (rbrusk): Ace 4.2 upgrade
 * * Jan  2 02:39 1996 (rbrusk): Ace 4.1 upgrade
 * Created: Jul 21 18:40 1995 (rbrusk)
 *-------------------------------------------------------------------
 */

/* $Id: cgraph.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */

#ifndef CGRAPH_H
#define CGRAPH_H
#define NEW_WIN32_GRAPHS
#if !defined(NEW_WIN32_GRAPHS)
#define OLD_WIN32_GRAPHS
#endif

#if !defined(VERBOSE_DEBUG)
//#define VERBOSE_DEBUG
//#define REALLY_VERBOSE
#endif

#include "CMultiAss.h"
#include "win32menus.h"
#include "graphbox.h" 
#include "msgwindow.h"

typedef  struct GraphStruct CGraphStruct; // to permit some transparent code changes for new design
class CGraph;  // a forward reference ?

// WIN32 mirror of ACEDB graph boxes...
// Found in CGraphWindow's and CSubGraphWindow's
typedef CTypedPtrArray<CObArray,CGraphBox*> BOXVIEWARRAY ;

#include "graphview.h"
#include "graphwin.h"

/////////////////////////////////////////////////////////////////////////////
// CGraph document

class CGraph : public CDocument
{
protected:
	CGraph();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CGraph)

#if !defined(NEW_WIN32_GRAPHS)  // old implementation...
// Attributes
	Graph_ m_pACEDB_Graph ; // pointer to ACEDB Graph_ data structure
#endif

public:
	static BOOL Frozen ;  	// a public semaphore to postpone gActive deactivation
							// by CGraphWindow::OnActivateView() during graphDevDestroy()
							// invoked CGraph::OnCloseDocument() process 

// Operations
public:
#if defined(OLD_WIN32_GRAPHS)  // old implementation...
	void SetGraph( Graph_ g ) { m_pACEDB_Graph = g ; }
	Graph_ GetGraph() const { return m_pACEDB_Graph; }
#endif

	CGraphView *GetGraphView() const
	{
		POSITION pos = GetFirstViewPosition() ;
		return( (CGraphView *)GetNextView( pos ) );
	}

	CGraphWindow *GetGraphWindow() const
	{
		CGraphView *pGraphView = GetGraphView() ;
		return( (CGraphWindow *)(pGraphView->GetParentFrame()) );
	}

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CGraph)
	public:
	virtual void OnCloseDocument();
	virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
	protected:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

protected:
	void InitGraph() ;


// Implementation
public:
	virtual ~CGraph();
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CGraph)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
} ;

// WIN32 gDev is a CGraph object, with members borrowed from xtsubs.c

#define GRAPHPTR ((CGraph *)gDev)

#if defined(OLD_WIN32_GRAPHS) 
#define DEVPTR   ((VIEWPTR)?((CGraphWindow *)(VIEWPTR->GetParentFrame())):0)	// Is_a CGraphWindow
#define SUBDEVPTR   ((GRAPHPTR)?((CGraphView *)(GRAPHPTR->GetGraphView())):(CGraphView *)0)					// Is_a CGraphView
#endif // defined(OLD_WIN32_GRAPHS) 
					
#if defined(NEW_WIN32_GRAPHS)

#define DEVPTR		GRAPHPTR   /* Is_a CGraph document == gDev */

// SUBDEVPTR checks if gSubDev is NULL
// then returns either a CGraphWindow or a CSubGraphWindow typed pointer
#define IsSubGraph(g) ( (g)->subdev && !((CObject *)(g)->subdev)->IsKindOf(RUNTIME_CLASS(CGraphWindow)) )? \
							TRUE  : /* if subdev is not NULL and is not a CGraphWindow, else */ FALSE 

#define SUBDEVPTR(fn,arg) if(gSubDev) \
							if( ((CObject *)gSubDev)->IsKindOf(RUNTIME_CLASS(CGraphWindow)) ) \
								((CGraphWindow *)gSubDev)->fn(arg) ; /* Isa Graph? */  \
							else \
								((CSubGraphWindow *)gSubDev)->fn(arg) ; /* isa Subgraph? */							
						
#define USESUBDEVPTR(s,fn,arg) if(gSubDev) { \
								if( ((CObject *)gSubDev)->IsKindOf(RUNTIME_CLASS(CGraphWindow)) ) \
									(s) = ((CGraphWindow *)gSubDev)->fn(arg) ; /* Isa Graph? */  \
								else \
									(s) = ((CSubGraphWindow *)gSubDev)->fn(arg) ; /* isa Subgraph? */ \
							   } else (s) = 0 ;
						
/* VIEWPTR should likely be set once somewhere, then dereferenced */
#define SETVIEWPORT		if(gSubDev) { \
							if( ((CObject *)gSubDev)->IsKindOf(RUNTIME_CLASS(CGraphWindow)) ) \
								VIEWPTR = (CGraphView *)((CGraphWindow *)gSubDev)->GetActiveView() ; /* Isa Graph? */  \
							else \
								VIEWPTR = (CGraphView *)gSubDev ; /* isa Subgraph? Use gSubDev directly? */ \
						} else VIEWPTR = (CGraphView *)0 ;
						
#endif  // defined(NEW_WIN32_GRAPHS)

#define FROZEN CGraph::Frozen 

// Resolution independent screen size is 1.3 x 1.0
#define SCREENX_ASPECT 1.3
#define SCREENY_ASPECT 1.0

extern "C" void setTextAspect( Graph_ g ) ;  // see typeInitialise() in graphcon.c
extern void	SetScrollBounds( CGraphView *pGraphView = NULL ) ;

#if defined(OLD_WIN32_GRAPHS)  
#define BOXVIEWS (DEVPTR->GetBoxViews())
#define BOXVIEW(i) (DEVPTR->BoxView(i))
#endif

#endif // defined(CGRAPH_H)

//********************** End of File ***************************************

 
 
 
