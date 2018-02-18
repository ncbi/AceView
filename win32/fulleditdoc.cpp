// $Id: fulleditdoc.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
// FullEditDoc.cpp : implementation of the CTextFullEditGraphDoc class
//

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h"
#include "FullEditDoc.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CTextFullEditGraphDoc

IMPLEMENT_DYNCREATE(CTextFullEditGraphDoc, CGraph)

BEGIN_MESSAGE_MAP(CTextFullEditGraphDoc, CGraph)
	//{{AFX_MSG_MAP(CTextFullEditGraphDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CTextFullEditGraphDoc construction/destruction

CTextFullEditGraphDoc::CTextFullEditGraphDoc()
{
	// TODO: add one-time construction code here

}

CTextFullEditGraphDoc::~CTextFullEditGraphDoc()
{
}

BOOL CTextFullEditGraphDoc::OnNewDocument()
{
	if (!CGraph::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
// CTextFullEditGraphDoc serialization

void CTextFullEditGraphDoc::Serialize(CArchive& ar)
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
// CTextFullEditGraphDoc diagnostics

#ifdef _DEBUG
void CTextFullEditGraphDoc::AssertValid() const
{
	CGraph::AssertValid();
}

void CTextFullEditGraphDoc::Dump(CDumpContext& dc) const
{
	CGraph::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CTextFullEditGraphDoc commands
 
