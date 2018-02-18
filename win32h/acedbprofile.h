/* $Id: acedbprofile.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// ACEDBProfile.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CACEDBProfile document

class CACEDBProfile : public CDocument
{
protected:
	CACEDBProfile();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CACEDBProfile)

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CACEDBProfile)
	public:
	virtual void OnFinalRelease();
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	protected:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CACEDBProfile();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CACEDBProfile)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
	// Generated OLE dispatch map functions
	//{{AFX_DISPATCH(CACEDBProfile)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_DISPATCH
	DECLARE_DISPATCH_MAP()
	DECLARE_INTERFACE_MAP()
};
 
