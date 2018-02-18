/* ACEDBProfile.cpp : implementation file
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
// $Id: acedbprofile.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $

#include "stdafx.h"
#include "WinAce.h"
#include "ACEDBProfile.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CACEDBProfile

IMPLEMENT_DYNCREATE(CACEDBProfile, CDocument)

CACEDBProfile::CACEDBProfile()
{
	EnableAutomation();
}

BOOL CACEDBProfile::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

CACEDBProfile::~CACEDBProfile()
{
}

void CACEDBProfile::OnFinalRelease()
{
	// When the last reference for an automation object is released
	// OnFinalRelease is called.  The base class will automatically
	// deletes the object.  Add additional cleanup required for your
	// object before calling the base class.

	CDocument::OnFinalRelease();
}


BEGIN_MESSAGE_MAP(CACEDBProfile, CDocument)
	//{{AFX_MSG_MAP(CACEDBProfile)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

BEGIN_DISPATCH_MAP(CACEDBProfile, CDocument)
	//{{AFX_DISPATCH_MAP(CACEDBProfile)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_DISPATCH_MAP
END_DISPATCH_MAP()

// Note: we add support for IID_IACEDBProfile to support typesafe binding
//  from VBA.  This IID must match the GUID that is attached to the 
//  dispinterface in the .ODL file.

// {219D81C0-B2BC-11CF-A8DB-444553540000}
static const IID IID_IACEDBProfile =
{ 0x219d81c0, 0xb2bc, 0x11cf, { 0xa8, 0xdb, 0x44, 0x45, 0x53, 0x54, 0x0, 0x0 } };

BEGIN_INTERFACE_MAP(CACEDBProfile, CDocument)
	INTERFACE_PART(CACEDBProfile, IID_IACEDBProfile, Dispatch)
END_INTERFACE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CACEDBProfile diagnostics

#ifdef _DEBUG
void CACEDBProfile::AssertValid() const
{
	CDocument::AssertValid();
}

void CACEDBProfile::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CACEDBProfile serialization

void CACEDBProfile::Serialize(CArchive& ar)
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
// CACEDBProfile commands
 
