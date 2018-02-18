// $Id: graphbox.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
// graphbox.cpp : implementation file
//

#include "stdafx.h"
#include "win32.h"

extern "C"
{
#include "regular.h"
}

#include "WinAce.h"
#include "CGraph.h"

#if defined(VERBOSE_DEBUG)
#undef VERBOSE_DEBUG
#endif

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CGraphBox

CGraphBox::CGraphBox( CView *pParent, int k, Box box )
{
	m_parentView = pParent ;
	m_boxMenu = NULL ;
	m_boxMenuType = NO_MENU ;
	m_boxData = box ;
	m_boxNo = k ;
}

CGraphBox::~CGraphBox()
{
	// delete any existing "Free Menu" associated with this graphBox
	// Other menus automatically deleted upon MDIChildFrame closure?
	m_parentView = NULL ;
	m_boxData = NULL ;
	m_boxNo = -1 ;
}

void *CGraphBox::SetMenu( BoxMenuType boxMenuType, void *boxMenu )  // (re)set box menu
{
#if defined(VERBOSE_DEBUG)
	TRACE("Entering CGraphBox[%d]::SetMenu( boxMenuType == %d, CObject *boxMenu == %p )\n",
		   m_boxNo, boxMenuType, boxMenu ) ;
	TRACE("\tm_boxMenuType == %d, m_boxMenu == %p\n",m_boxMenuType, m_boxMenu) ;
#endif
	if(	m_boxMenu == boxMenu )
	{
//		TRACE("\n\n**** WARNING: Duplicate SetMenu( boxMenu == %p ) ignored?\n\n", boxMenu ) ;
		return m_boxMenu ;
	}
	
	m_boxMenuType = boxMenuType ;
	return (m_boxMenu = boxMenu) ;
 }

/////////////////////////////////////////////////////////////////////////////
// CGraphBox message handlers

 
