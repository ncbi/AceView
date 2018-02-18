// $Id: acefiledoc.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
// AceFileDoc.cpp : implementation file
//

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "AceFileDoc.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CAceFileDoc

IMPLEMENT_DYNCREATE(CAceFileDoc, CDocument)

CAceFileDoc::CAceFileDoc()
{
}

BOOL CAceFileDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

CAceFileDoc::~CAceFileDoc()
{
}


BEGIN_MESSAGE_MAP(CAceFileDoc, CDocument)
	//{{AFX_MSG_MAP(CAceFileDoc)
	ON_COMMAND(ID_FILE_SEND_MAIL, OnFileSendMail)
	ON_UPDATE_COMMAND_UI(ID_FILE_SEND_MAIL, OnUpdateFileSendMail)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()



/////////////////////////////////////////////////////////////////////////////
// CAceFileDoc diagnostics

#ifdef _DEBUG
void CAceFileDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CAceFileDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CAceFileDoc serialization

void CAceFileDoc::Serialize(CArchive& ar)
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
// CAceFileDoc commands
 
