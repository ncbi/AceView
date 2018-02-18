/* $Id: preferences.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// Preferences.h : header file
//

#include "DBProfile.h"
#include "UserInterface.h"
#include "CGraph.h"
#include "FontPreference.h"

/////////////////////////////////////////////////////////////////////////////
// CPreferences

#define DBPROFILEPAGE 0x00000001 /* Environment Parameters Page */
#define USRINTFCPAGE  0x00000002 /* User Interface Parameters Page */
#define FONTPREFPAGE  0x00000004 /* Font Parameters Page */
#define ALLPAGES (DBPROFILEPAGE | USRINTFCPAGE | FONTPREFPAGE)

class CPreferences : public CPropertySheet
{
	DECLARE_DYNAMIC(CPreferences)

// Construction
protected:
	virtual void InitSheet(UINT specs) ; // Common code for construction
public:
	CPreferences(UINT nIDCaption, UINT specs = ALLPAGES, CWnd* pParentWnd = NULL, UINT iSelectPage = 0);
	CPreferences(LPCTSTR pszCaption, UINT specs = ALLPAGES, CWnd* pParentWnd = NULL, UINT iSelectPage = 0);

// Attributes
protected:
	CDBProfilePage *m_DBProfilePage ;
	CUserInterfacePage *m_UserInterfacePage ;
	CFontPreferencePage *m_FontPreferencePage ;

public:

// Operations
public:
	// Returns true if ACEDB needs to be reloaded due to parameter changes
	virtual BOOL ReloadACEDB() ; 

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPreferences)
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CPreferences();

	// Generated message map functions
protected:
	//{{AFX_MSG(CPreferences)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
 
