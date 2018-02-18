// $Id: splashbox.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
// SplashBox.cpp : implementation file
//

#include "stdafx.h"
#include "WinAce.h"
#include "SplashBox.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CSplashBox dialog


CSplashBox::CSplashBox(CWnd* pParent /*=NULL*/)
	: CDialog(CSplashBox::IDD, pParent)
{
	//{{AFX_DATA_INIT(CSplashBox)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
}


void CSplashBox::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CSplashBox)
		// NOTE: the ClassWizard will add DDX and DDV calls here
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CSplashBox, CDialog)
	//{{AFX_MSG_MAP(CSplashBox)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSplashBox message handlers
 
