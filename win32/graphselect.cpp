/*  File: graphselect.cpp : implementation file
 *  Author:	Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *			adapted from code written by Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 *  WIN32 version of graphSelect written using the Microsoft Windows MFC.
 * Exported functions:
          BOOL graphSelect
 * HISTORY:
 * Last edited: Jan 4 00:10 1996 (rmb): Ported graphSelect() to WIN32
 * UNIX Version last edited: Jun 30 11:49 1993 (mieg)
 * Created: Mon Jul 20 22:25:49 1992 (mieg)
 *-------------------------------------------------------------------
 * gha  04/30/93  change spacing, bold
 * gha  05/03/93  change title, layout
 */

// $Id: graphselect.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "graphselect.h"

extern "C"
{
#include "acedb.h"
#include "keyset.h"
}

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

#define new DEBUG_NEW

/////////////////////////////////////////////////////////////////////////////
// CSelector dialog


CSelector::CSelector(FREEOPT *options, CWnd* pParent /*= NULL*/)
	: CDialog(CSelector::IDD, pParent)
{
	//{{AFX_DATA_INIT(CSelector)
	//}}AFX_DATA_INIT
	m_selectedKey = 0 ;
	m_options = options ;
}

void CSelector::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CSelector)
	DDX_Control(pDX, IDC_KEYLIST, m_keyList);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CSelector, CDialog)
	//{{AFX_MSG_MAP(CSelector)
	ON_CBN_SETFOCUS(IDC_KEYLIST, OnSetFocusKeyList)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CSelector message handlers

static char *title ;

BOOL CSelector::OnInitDialog() 
{
	title = m_options->text ;
	if( !(title && *title) )
		title = "Selector:" ;
	
	CDialog::OnInitDialog();

	SetWindowText( (LPCSTR)title ) ;

	for(int i = 1 ; i <= m_options->key; ++i )
	{
		int keyIdx = m_keyList.AddString( (LPCSTR) (m_options[i].text) ) ;
		if( keyIdx == CB_ERR )
			messcrash("Selector key creation error?") ;
		else if( keyIdx == CB_ERRSPACE ) 
			messcrash("Selector key creation error: Out of String Space?") ;
		else m_keyList.SetItemData( keyIdx, (DWORD) (m_options[i].key) ) ;
	}
	m_keyList.SetCurSel(0) ;

	return TRUE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}

void CSelector::OnOK() 
{
	int keyIdx = m_keyList.GetCurSel() ;
	if( keyIdx == CB_ERR )
		m_selectedKey = 0 ;
	else
	{
		DWORD key = m_keyList.GetItemData( keyIdx ) ;
		if( key == CB_ERR )
			m_selectedKey = 0 ;
		else 
			m_selectedKey = (KEY)key ;
	}

	CDialog::OnOK();
}

BOOL graphSelect (KEY *kpt, FREEOPT *options)
{ 
	CSelector Selector(options) ;
	int retCode = Selector.DoModal() ;
	if( retCode == -1 || retCode == IDABORT )
	{
		ASSERT(FALSE) ;  // if could not succeed in opening Selector?
		return FALSE;
	}
	*kpt = Selector.GetSelectedKey() ;
	if( retCode  == IDOK && *kpt)
		return TRUE ;
	else
		return FALSE ;
}


void CSelector::OnSetFocusKeyList() 
{
		m_keyList.ShowDropDown() ;		
}
 
