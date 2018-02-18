/*  GraphWin.h : header file
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
*/
/* $Id: graphwin.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */

#ifndef GRAPHWIN_H
#define GRAPHWIN_H

/////////////////////////////////////////////////////////////////////////////
// CGraphWindow frame

class CGraphWindow : public CMDIChildWnd
{
	DECLARE_DYNCREATE(CGraphWindow)

private:
	// left and top are screen coordinates?
	static int 	m_def_left, m_def_top, m_def_width, m_def_height ;

	// WIN32 instantiation of boxes
	BOXVIEWARRAY m_boxViews ;

#if defined(NEW_WIN32_GRAPHS) 
	CGraphStruct *m_pACEDB_Graph ;
#endif

	UINT m_MessageWndID ;
	CMsgWindow *m_pMessageWnd ; // graphMessage(), graphUnMessage() item

public:
	CGraphWindow();		// protected constructor used by dynamic creation

// Operations
public:
#if defined(NEW_WIN32_GRAPHS)  
	static void CreateSubGraph(Graph_ graph, float x, float y, float w, float h) ;
#endif

#if defined(NEW_WIN32_GRAPHS)
	void Activate(BOOL flg = TRUE) ;
#endif

	BOOL SetActiveGraph( BOOL ReDraw = TRUE ) ;
	CGraphStruct *GetGraph( BOOL flg = TRUE ) ;	// Gets ACEDB graph associated with this frame window
												// flg is a dummy parameter
#if !defined(NEW_WIN32_GRAPHS)
	BOXVIEWARRAY &GetBoxViews() { return m_boxViews ; }
#else
	BOXVIEWARRAY *GetBoxViews(BOOL flg = TRUE) { return &m_boxViews ; } // flg is a dummy parameter
#endif

	CGraphBox *BoxView(int i) { return m_boxViews[i] ; }

	static void InitDevStruct(Graph_ graph, int x, int y, int w, int h) ;

	static CSize GetGraphSize()  
		{  CSize winSize(m_def_width,m_def_height); return winSize; }

	static CPoint GetGraphOrigin()  
		{  CPoint winOrigin(m_def_left,m_def_top); return winOrigin ; }

	static CRect GetGraphRect()  
	{
		CRect winRect ( m_def_left, m_def_top,
						m_def_left+m_def_width,	
						m_def_top+m_def_height );
		return winRect ;
	}

	static int nextControlID() ;

	void OpenGraphMessage( const CString& text ) ;
	void CloseGraphMessage( BOOL flg = TRUE ) ;  // flg is a dummy parameter

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CGraphWindow)
	public:
	protected:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	virtual void PostNcDestroy();
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CGraphWindow();
	void DeleteBoxViews() ;	// used in destructor

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CGraphWindow)
	afx_msg BOOL OnQueryNewPalette();
	afx_msg void OnPaletteChanged(CWnd* pFocusWnd);
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnChildActivate();
	afx_msg void OnDownArrow();
	afx_msg void OnLeftArrow();
	afx_msg void OnRightArrow();
	afx_msg void OnUpArrow();
	afx_msg void OnActivate(UINT nState, CWnd* pWndOther, BOOL bMinimized);
	afx_msg UINT OnNcHitTest(CPoint point);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#endif // #ifndef GRAPHWIN_H


 
 
 
