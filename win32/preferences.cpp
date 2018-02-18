// $Id: preferences.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
// Preferences.cpp : implementation file
//

#include "stdafx.h"
#include "WinAce.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include "Preferences.h"

/////////////////////////////////////////////////////////////////////////////
// CPreferences

IMPLEMENT_DYNAMIC(CPreferences, CPropertySheet)

CPreferences::CPreferences(UINT nIDCaption, UINT specs, CWnd* pParentWnd, UINT iSelectPage)
	:CPropertySheet(nIDCaption, pParentWnd, iSelectPage)
{
	InitSheet(specs) ;
}

CPreferences::CPreferences(LPCTSTR pszCaption, UINT specs, CWnd* pParentWnd, UINT iSelectPage)
	:CPropertySheet(pszCaption, pParentWnd, iSelectPage)
{
	InitSheet(specs) ;
}

CPreferences::~CPreferences()
{
	if(m_DBProfilePage != NULL) delete m_DBProfilePage ;
	if(m_UserInterfacePage != NULL) delete m_UserInterfacePage ;
	if(m_FontPreferencePage != NULL) delete m_FontPreferencePage ;
}


BEGIN_MESSAGE_MAP(CPreferences, CPropertySheet)
	//{{AFX_MSG_MAP(CPreferences)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPreferences message handlers

void CPreferences::InitSheet(UINT specs) 
{
	if(specs & DBPROFILEPAGE)
	{
		m_DBProfilePage = new CDBProfilePage ;
		AddPage(m_DBProfilePage) ;
	}
	else m_DBProfilePage = NULL ;

	if(specs & USRINTFCPAGE)
	{
		m_UserInterfacePage = new CUserInterfacePage ;
		AddPage(m_UserInterfacePage) ;
	}
	else m_UserInterfacePage = NULL ;


	if(specs & FONTPREFPAGE)
	{
	    m_FontPreferencePage = new CFontPreferencePage ;
	    AddPage(m_FontPreferencePage) ;
	}
	else m_FontPreferencePage = NULL ;

}

BOOL CPreferences::ReloadACEDB()
{
	BOOL reload = FALSE ;

	TRACE("\n*** Entering CPreferences::ReloadACEDB()\n") ;

	if( m_DBProfilePage != NULL && m_DBProfilePage->m_ReloadACEDB )
			reload = TRUE ;

	// if( m_UserInterfacePage != NULL ) { }
	
	//if(m_FontPreferencePage != NULL) { }

	return reload ;
}
 
