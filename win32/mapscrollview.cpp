/* MapScrollView.cpp : implementation of the CMapScrollGraphView class
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
// $Id: mapscrollview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
//

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "MapScrollView.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CMapScrollGraphView

IMPLEMENT_DYNCREATE(CMapScrollGraphView, CGraphView)

BEGIN_MESSAGE_MAP(CMapScrollGraphView, CGraphView)
	//{{AFX_MSG_MAP(CMapScrollGraphView)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CGraphView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CGraphView::OnFilePrintPreview)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CMapScrollGraphView construction/destruction

CMapScrollGraphView::CMapScrollGraphView()
{
	// TODO: add construction code here

}

CMapScrollGraphView::~CMapScrollGraphView()
{
}

/////////////////////////////////////////////////////////////////////////////
// CMapScrollGraphView printing


/////////////////////////////////////////////////////////////////////////////
// CMapScrollGraphView diagnostics

#ifdef _DEBUG
void CMapScrollGraphView::AssertValid() const
{
	CGraphView::AssertValid();
}

void CMapScrollGraphView::Dump(CDumpContext& dc) const
{
	CGraphView::Dump(dc);
}

CGraph* CMapScrollGraphView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CGraph)));
	return (CGraph*)m_pDocument;
}

#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CMapScrollGraphView message handlers
 
 
