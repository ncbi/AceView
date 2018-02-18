/* windialogs.cpp : implementation file
 *  Author: Richard Bruskiewich	(rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: WIN32 implementations of user graph*() input I/O functions.  
 * Exported:  graphPrompt(), SysBeep(), graphOut(), graphQuery() 
 * HISTORY:
 * Last edited:  
 * Created: Aug 25 11:56 1995 (rmb)
 *-------------------------------------------------------------------
 */
// $Id: windialogs.cpp,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $

#include "stdafx.h"
#include "win32.h"

extern "C"
{
#include "regular.h" /* for freesubs */
}
#include <string.h>

#include "WinAce.h"
#include "windialogs.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CGraphPrompt dialog


CGraphPrompt::CGraphPrompt(CWnd* pParent /*=NULL*/)
	: CDialog(CGraphPrompt::IDD, pParent)
{
	//{{AFX_DATA_INIT(CGraphPrompt)
	m_userdata = _T("");
	m_prompt = _T("");
	//}}AFX_DATA_INIT
}


void CGraphPrompt::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CGraphPrompt)
	DDX_Text(pDX, IDC_ACEDB_PROMPT_USERDATA, m_userdata);
	DDX_Text(pDX, IDC_ACEPROMPT, m_prompt);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CGraphPrompt, CDialog)
	//{{AFX_MSG_MAP(CGraphPrompt)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// ACEDB interface functions

extern "C" void win32SysBeep(int n) { MessageBeep(MB_OK); }

extern BOOL blockRedraw = FALSE ;
extern "C" void win32GraphOut(char *message)
{
	blockRedraw = TRUE ;
	MessageBox(NULL,message,"ACEDB Message",MB_OK | MB_APPLMODAL );
	blockRedraw = FALSE ;
}

extern "C" char *win32GraphPrompt( char *prompt, char *dfault )
{
   static char data[80] ;  // Hopefully enough space?

	CGraphPrompt dlg ;
	dlg.m_userdata = dfault ;
	dlg.m_prompt = prompt ;

	blockRedraw = TRUE ;
	int Result = dlg.DoModal() ;
	blockRedraw = FALSE ;

	if( Result == IDOK )
	{
		strcpy(data,dlg.m_userdata) ;
		return &data[0] ;
	}
	else	/* Cancel! */
		return NULL ;
}

extern "C" BOOL win32GraphQuery (char *text)
{
	blockRedraw = TRUE ;
	int Result = MessageBox(NULL,text,"ACEDB Message",MB_YESNO) ;
	blockRedraw = FALSE ;

	if( Result == IDYES )
		return TRUE ;
	else	/* Cancel! */
		return FALSE ;
}


/////////////////////////////////////////////////////////////////////////////
// CGraphPrompt message handlers

void CGraphPrompt::OnOK() 
{
	UpdateData(TRUE) ;
	if ( m_userdata.IsEmpty() )
		CDialog::OnCancel(); /* blank data => Cancel? */
	else
		CDialog::OnOK();
}
 
