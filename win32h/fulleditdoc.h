/* $Id: fulleditdoc.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// FullEditDoc.h : interface of the CTextFullEditGraphDoc class
//
/////////////////////////////////////////////////////////////////////////////

class CTextFullEditGraphDoc : public CGraph
{
protected: // create from serialization only
	CTextFullEditGraphDoc();
	DECLARE_DYNCREATE(CTextFullEditGraphDoc)

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTextFullEditGraphDoc)
	public:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CTextFullEditGraphDoc();
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CTextFullEditGraphDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
 
