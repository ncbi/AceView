// $Id: userinterface.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
// UserInterface.cpp : implementation file
//

#include "stdafx.h"
#include "WinAce.h"
#include "UserInterface.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CUserInterfacePage property page

IMPLEMENT_DYNCREATE(CUserInterfacePage, CPropertyPage)

CUserInterfacePage::CUserInterfacePage() : CPropertyPage(CUserInterfacePage::IDD)
{
	//{{AFX_DATA_INIT(CUserInterfacePage)
	m_ACEDB_NO_BANNER = FALSE;
	m_EmailAddress = _T("");
	m_WebBrowser = _T("");
	//}}AFX_DATA_INIT
}

CUserInterfacePage::~CUserInterfacePage()
{
}

void CUserInterfacePage::DoDataExchange(CDataExchange* pDX)
{
	CPropertyPage::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CUserInterfacePage)
	DDX_Control(pDX, IDC_WEB_BROWSER, m_WEB_BROWSER);
	DDX_Control(pDX, IDC_EMAIL_ADDRESS, m_EMAIL_ADDRESS);
	DDX_Check(pDX, IDC_ACEDB_NO_BANNER, m_ACEDB_NO_BANNER);
	DDX_Text(pDX, IDC_EMAIL_ADDRESS, m_EmailAddress);
	DDX_Text(pDX, IDC_WEB_BROWSER, m_WebBrowser);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CUserInterfacePage, CPropertyPage)
	//{{AFX_MSG_MAP(CUserInterfacePage)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CUserInterfacePage message handlers

BOOL CUserInterfacePage::OnInitDialog() 
{
	CPropertyPage::OnInitDialog();
	
	m_ACEDB_NO_BANNER = (getProgParameter("ACEDB_NO_BANNER") == NULL) ? FALSE : TRUE ;
	m_EmailAddress = getProgParameter("EMAIL_ADDRESS") ;
	m_WebBrowser = getProgParameter("WEB_BROWSER") ;

	return FALSE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}

void CUserInterfacePage::OnOK() 
{
	TRACE("\tUserInterfacePage::OnOK()\n") ;
	if( m_WEB_BROWSER.GetModify() )
	{
		TRACE("\t\tSetting Web Browser Program\t== %s\n",m_WebBrowser) ;
		setProgParameter((LPCTSTR)("WEB_BROWSER="+m_WebBrowser) ) ;
	}

	if( m_EMAIL_ADDRESS.GetModify() )
	{
		TRACE("\t\tSetting EMail Address\t\t== %s\n",m_EmailAddress) ;
		setProgParameter((LPCTSTR)("EMAIL_ADDRESS="+m_EmailAddress) ) ;
	}

	TRACE("\t\tSetting ACEDB_NO_BANNER\t\t== %s\n",
		m_ACEDB_NO_BANNER? "TRUE":"FALSE") ;
	if(m_ACEDB_NO_BANNER)
		setProgParameter("ACEDB_NO_BANNER=1" ) ;
	else
		setProgParameter("ACEDB_NO_BANNER="/* NULL */ ) ;

	CPropertyPage::OnOK();
}
 
