/* $Id: acefiledoc.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// AceFileDoc.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CAceFileDoc document

class CAceFileDoc : public CDocument
{
protected:
	CAceFileDoc();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CAceFileDoc)

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAceFileDoc)
	protected:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CAceFileDoc();
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CAceFileDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
 
