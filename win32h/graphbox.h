/* $Id: graphbox.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// graphbox.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CGraphBox window

class CGraphBox : public CObject
{
// Construction
public:
	CGraphBox ( CView *pParent, int k, Box box ) ;

// Attributes
protected:
	CView *m_parentView ;
	BoxMenuType m_boxMenuType ;
	void *m_boxMenu ; // Is a MENUBITS structure; see ACEMenu.h
	Box m_boxData ;
	int m_boxNo ;

// Operations
public:

	BoxMenuType GetMenuType() const { return m_boxMenuType ; }

	void *SetMenu( BoxMenuType boxMenuType, void *boxMenu = NULL ) ;  // (re)set box menu
	void *GetMenu() { return( m_boxMenu ) ; }

// Implementation
public:
	virtual ~CGraphBox();

};

/////////////////////////////////////////////////////////////////////////////
 
