/* ACEDBProfileView.cpp : implementation file
 *  Author: Richard Bruskiewich	(rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Property page dialog for database configuration.  
 * Exported:  class CACEDBProfile
 * HISTORY:
 * Last edited: Jun 10 10:00 1996 (rbrusk): WIN32 to 4.3
 * Created: May 1996 (rbrusk)
 *-------------------------------------------------------------------
 */
// $Id: acedbprofileview.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $


#include "stdafx.h"
#include "WinAce.h"
#include "ACEDBProfileView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CACEDBProfileView

IMPLEMENT_DYNCREATE(CACEDBProfileView, CFormView)

CACEDBProfileView::CACEDBProfileView()
	: CFormView(CACEDBProfileView::IDD)
{
	//{{AFX_DATA_INIT(CACEDBProfileView)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
}

CACEDBProfileView::~CACEDBProfileView()
{
}

void CACEDBProfileView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CACEDBProfileView)
		// NOTE: the ClassWizard will add DDX and DDV calls here
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CACEDBProfileView, CFormView)
	//{{AFX_MSG_MAP(CACEDBProfileView)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CACEDBProfileView diagnostics

#ifdef _DEBUG
void CACEDBProfileView::AssertValid() const
{
	CFormView::AssertValid();
}

void CACEDBProfileView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CACEDBProfileView message handlers
 
