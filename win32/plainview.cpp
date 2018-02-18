/* PlainView.cpp : implementation file
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
// $Id: plainview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
//

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "PlainView.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CPlainGraphView

IMPLEMENT_DYNCREATE(CPlainGraphView, CGraphView)

CPlainGraphView::CPlainGraphView()
{
}

CPlainGraphView::~CPlainGraphView()
{
}


BEGIN_MESSAGE_MAP(CPlainGraphView, CGraphView)
	//{{AFX_MSG_MAP(CPlainGraphView)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CPlainGraphView diagnostics

#ifdef _DEBUG
void CPlainGraphView::AssertValid() const
{
	CGraphView::AssertValid();
}

void CPlainGraphView::Dump(CDumpContext& dc) const
{
	CGraphView::Dump(dc);
}

CGraph* CPlainGraphView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CGraph)));
	return (CGraph*)m_pDocument;
}

#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CPlainGraphView message handlers

 
 
