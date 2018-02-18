/* msgwindow.cpp : implementation file
 *  Author: Richard Bruskiewich	(rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:   Component of WIN32 emulation of popup graphMessage() window
 * Exported:  class CMsgWindow
 * HISTORY:
 * Last edited: Sep 5 10:17 1996 (rbrusk): m_pOwner
 * *	Jan 14 14:45 1996 (rbrusk)
 * Created: Aug 25 11:56:08 1995 (rmb)
 *-------------------------------------------------------------------
 */
// $Id: msgwindow.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $


#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CMsgWindow dialog

// This value should be adjusted as necessary
// to suit CMsgWindow dialog edit control width
// Leave it about 10 characters shorter than the width
// to allow for wrap around completion to next word boundary
#define MAX_LINE_LENGTH 50
#define NEWLINE "\r\n"

CMsgWindow::CMsgWindow( const CString& msgTitle, const CString& msgText, CWnd* pOwner /*=NULL*/)
	: CDialog(CMsgWindow::IDD, pOwner )
{
	m_pOwner = pOwner ; // Who owns the message window?
	/* convert msgText from Unix to MSDOS format */
	CString DOSText ;
	int nc = msgText.GetLength() ;
	int i,	 /* current source character position */
		lpos /* current output character position */ ;

	for( lpos = i = 0; i< nc; )
	{
		// If source text EOL...
		if( msgText[i] == '\n' )
		{
			while(	i < nc && msgText[i] == '\n')
			// Then, convert all consecutive
			// source newlines to carriage return+newline?
			{ DOSText += NEWLINE ; ++i ; }
			// then reset output line character position to start of line
			lpos = 0 ;
		}
		else
			// Wrap the text around whenever the
			// line gets too long?  Break at next whitespace or EOL 
			if( !(++lpos % MAX_LINE_LENGTH) )
			{
				char c ; // shorthand for msgText[i]...
				// Find next word boundary break
				while(	i < nc 	// while characters still available   
						&& 		// and not EOL or whitespace or hyphen...
						!((c = msgText[i]) == '\n' ||
						 c == ' ' ||
						 c == '-'  ||
						 c == '\t'    ) )
					{ DOSText += c ;  ++i ; } // ...copy character

				// If NOT EOL
				if(	i< nc )
					switch( c )
					{
						case '-': 
							// If a hyphen was encountered, append it to output
							DOSText += '-' ;
							// Then, terminate line	with NEWLINE and skip to next character
							DOSText += NEWLINE ; ++i ;
							break ;
						case '\n':
							// If source text EOL encountered
							while(	i < nc && msgText[i] == '\n')
							// Then, convert all consecutive
							// source newlines to NEWLINE?
							{ DOSText += NEWLINE ; ++i ; }
							break ;
						default: /* White space? */
							// Then, just terminate line with NEWLINE and skip to next character
							DOSText += NEWLINE ;
					}

				// Skip leading whitespace to end of message or next word
				while(	i < nc 	// while characters still available   
						&& 		// and whitespace seen
						((c = msgText[i]) == ' ' || c == '\t') )
					++i ; // ...skip whitespace
				
				// then, reset line character position to start of line
				lpos = 0 ;
			}
			else
				{ DOSText += msgText[i] ; ++i ; } // else, just copy current character and skip
	}
	
	//{{AFX_DATA_INIT(CMsgWindow)
	m_MsgText = _T(DOSText);
	//}}AFX_DATA_INIT

	if( !Create( CMsgWindow::IDD, pOwner ) )
		messcrash("Could not create graphMessage?") ;
	SetWindowText( msgTitle ) ;
}


void CMsgWindow::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CMsgWindow)
	DDX_Control(pDX, IDC_MSGTEXT, m_MsgWindow);
	DDX_Text(pDX, IDC_MSGTEXT, m_MsgText);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CMsgWindow, CDialog)
	//{{AFX_MSG_MAP(CMsgWindow)
	ON_WM_CLOSE()
	ON_EN_SETFOCUS(IDC_MSGTEXT, OnFocus)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CMsgWindow message handlers

void CMsgWindow::PostNcDestroy() 
{																	   
	delete this ; 
}

void CMsgWindow::OnClose() 
{
    if( m_pOwner->IsKindOf(RUNTIME_CLASS(CGraphWindow)) )
		((CGraphWindow *)m_pOwner)->CloseGraphMessage() ;
	else
		((CSubGraphWindow *)m_pOwner)->CloseGraphMessage() ;

}

void CMsgWindow::OnFocus() 
{
	m_MsgWindow.SetSel(-1,0,FALSE) ;
}
 
 
