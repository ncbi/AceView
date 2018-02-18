// $Id: dbprofile.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
// DBProfile.cpp : implementation file
//

#include "stdafx.h"
#include "WinAce.h"

extern "C" {
#include "mystdlib.h"
}

#include "DBProfile.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDBProfilePage property page

IMPLEMENT_DYNCREATE(CDBProfilePage, CPropertyPage)

CDBProfilePage::CDBProfilePage() : CPropertyPage(CDBProfilePage::IDD)
{
	//{{AFX_DATA_INIT(CDBProfilePage)
	m_ACEDB_Directory = _T("");
	m_ACEDB_DATA_Directory = _T("");
	m_ACEDB_COMMON_Directory = _T("");
	//}}AFX_DATA_INIT
}

CDBProfilePage::~CDBProfilePage()
{
}

void CDBProfilePage::DoDataExchange(CDataExchange* pDX)
{
	CPropertyPage::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDBProfilePage)
	DDX_Control(pDX, IDC_ACEDB_DATA, m_ACEDB_DATA);
	DDX_Control(pDX, IDC_ACEDB_COMMON, m_ACEDB_COMMON);
	DDX_Control(pDX, IDC_ACEDB, m_ACEDB);
	DDX_Text(pDX, IDC_ACEDB, m_ACEDB_Directory);
	DDX_Text(pDX, IDC_ACEDB_DATA, m_ACEDB_DATA_Directory);
	DDX_Text(pDX, IDC_ACEDB_COMMON, m_ACEDB_COMMON_Directory);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDBProfilePage, CPropertyPage)
	//{{AFX_MSG_MAP(CDBProfilePage)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDBProfilePage message handlers

BOOL CDBProfilePage::OnInitDialog() 
{
	CPropertyPage::OnInitDialog();

	m_ReloadACEDB = FALSE ;
	m_ACEDB_Directory = getProgParameter("ACEDB") ;
	m_ACEDB_DATA_Directory = getProgParameter("ACEDB_DATA") ;
	m_ACEDB_COMMON_Directory = getProgParameter("ACEDB_COMMON") ;

	return FALSE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}


void CDBProfilePage::OnOK() 
{
	TRACE("\tDBProfilePage::OnOK()\n") ;

	if( m_ACEDB.GetModify() )
	{
		TRACE("\t\tSetting ACEDB\t\t\t\t== %s\n",m_ACEDB_Directory) ;
		setProgParameter((LPCTSTR)("ACEDB="+m_ACEDB_Directory) ) ;
		m_ReloadACEDB = TRUE ;
	}

	if( m_ACEDB_DATA.GetModify() )
	{
		TRACE("\t\tSetting ACEDB_DATA\t\t\t== %s\n",m_ACEDB_DATA_Directory) ;
		setProgParameter((LPCTSTR)("ACEDB_DATA="+m_ACEDB_DATA_Directory) ) ;
		m_ReloadACEDB = TRUE ;
	}

	if( m_ACEDB_COMMON.GetModify() )
	{
		TRACE("\t\tSetting ACEDB_COMMON\t\t== %s\n",m_ACEDB_COMMON_Directory) ;
		setProgParameter((LPCTSTR)("ACEDB_COMMON"+m_ACEDB_COMMON_Directory) ) ;
		m_ReloadACEDB = TRUE ;
	}
	
	CPropertyPage::OnOK();
}
 
