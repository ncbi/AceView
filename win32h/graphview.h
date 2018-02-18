/* graphview.h : header file for graphview.cpp
 *  Author: Richard Bruskiewich, rbrusk@octogene.medgen.ubc.ca
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * HISTORY:
 * Last edited: Jun 8 14:00 1996 (rbrusk):
 *		-	merged graphwin32lib.cpp with graphview.cpp:
 * Created: Jul 21 18:40 1995 (rbrusk)
 *********************************************************************************/

/* $Id: graphview.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */

#ifndef GRAPHVIEW_H
#define GRAPHVIEW_H

//
//******************** Font Constant Definitions *******************

// Enumeration of ACEDB abstract font faces
enum FONTFACE {	Plain,FixedWidth,Bold,Italic,Greek,NUM_FONTFACES} ;

// Enumeration of ACEDB abstract font sizes
enum FONTSIZE {	DefaultSize,Tiny,Small,Medium,Large,Huge,NUM_FONTSIZES} ;

#define NUM_FONTS NUM_FONTFACES*NUM_FONTSIZES
#define FNTIDX(format,height) format*NUM_FONTSIZES+height

class CGraphWindow ; // forward reference
/////////////////////////////////////////////////////////////////////////////
// CGraphView view

class CGraphView : public CScrollView
{
protected:
	void InitDrawing(Box box) ;

	CGraphView() ;           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CGraphView)

// Attributes
private:

#if defined(NEW_WIN32_GRAPHS)  // new alternate implementation...
// Attributes
	// left and top are screen coordinates?
	static int 	m_def_left, m_def_top, m_def_width, m_def_height ;

	// WIN32 instantiation of boxes
	BOXVIEWARRAY m_boxViews ;
	CGraphStruct *m_pACEDB_Graph ;

	UINT m_MessageWndID ;
	CMsgWindow *m_pMessageWnd ; // graphMessage(), graphUnMessage() item

#endif

	// ACEDB linewidth converted to GDI logical coordinates
	int m_currentLineWidth ;

	// Current ACEDB foreground and background colours
	unsigned char	m_currentFCol, m_currentBCol ;
    CBrush *m_currentBrush ;

	struct pixelImage
	{
		int w, h, len ;
		CBitmap *pBitMap ;
		pixelImage *next ;
	} ;
	struct imageSet
	{
 		char *data ;  // master data for bitmaps
		struct pixelImage *next ;
	} ;
	typedef struct imageSet IMAGESET, *LPIMAGESET ;
	typedef struct pixelImage   PIXELIMAGE,   *LPPIXELIMAGE ;
	typedef CTypedPtrMap<CMapPtrToPtr, char*, LPIMAGESET> IMAGEMAP;
	IMAGEMAP m_images ;

protected:
	// The font palettes are shared by all graphs...
	typedef CFont *FONTPALETTE[NUM_FONTS] ;
	// The "palette" of fonts available for screen display
	static FONTPALETTE screenFont ;
	// The "palette" of fonts available for printing on current page
	static FONTPALETTE printFont ;
	// A pointer to the currently drawable palette of CFont pointers: screen or printer

public:
#if defined(NEW_WIN32_GRAPHS)  // new alternate implementation...
// Operations

	static void InitSubDevStruct(Graph_ graph, int x, int y, int w, int h) ;

	static CGraphView *GetNewSubGraphView(int type) ; // returns a CGraphView object of given subclass type

	static CRect GetGraphRect()  
		{
			CRect winRect ( m_def_left, m_def_top,
							m_def_left+m_def_width,	
							m_def_top+m_def_height );
			return winRect ;
		}

	void SetGraph( Graph_ g ) { m_pACEDB_Graph = g ; }
	BOOL SetActiveGraph( BOOL ReDraw = TRUE ) ;
	void Activate(BOOL flg = TRUE) ;  // flg is a dummy parameter

	// Gets ACEDB graph associated with this frame window
	// flg is a dummy parameter
	CGraphStruct *GetGraph( BOOL flg = TRUE ) const { return m_pACEDB_Graph; } 

	// Currently active viewport; set whenever gDev and gSubDev are set?
	static CGraphView *gViewPort ; 

	BOXVIEWARRAY *GetBoxViews(BOOL flg = TRUE) { return &m_boxViews ; } // flg is a dummy parameter
	CGraphBox *BoxView(int i) { return m_boxViews[i] ; }

	static int nextControlID() ;
	void OpenGraphMessage( const CString& text ) ;
	void CloseGraphMessage( BOOL flg = TRUE )  ; // flg is a dummy parameter) ;

#else
	CGraphStruct *GetGraph( BOOL flg = TRUE ) ;		// returns associated ACEDB Graph_ object
													// flg is a dummy parameter 
#endif // defined(NEW_WIN32_GRAPHS)

	// Actual machine dependent WIN32 display characteristics
	// gScreenX, gScreenY are the display dimensions in pixels
	static int 	gScreenX, gScreenY,
	 			gNumPixelBits, gNumGraphPlanes, gColourRes  ;
	static void GetGraphDevAttributes() ;

	// Global class variable for currently active GDI device context
	static CDC *gDC ;

	// Previously set ACEDB foreground and background colours
	// Not guaranteed to be initialized unless setColors() has been called
	unsigned char oldfcol, oldbcol ;

	// Win32 Global variables and macros for scaling drawBox() GDI
	// to screen or to printer device dimensions
	static float gXAspect, gYAspect ;

	static CString	fontFace[NUM_FONTFACES] ; // Actual font face name strings
	static int		fontSize[NUM_FONTSIZES] ; // Actual font heights in logical units

	static CFont *AceFont( FONTFACE face, FONTSIZE size ) ;
	static FONTPALETTE *gFont ;

	CSize 	graphSize,
			pageSize,
			lineSize ;
	static HCURSOR	Arrow ,	HourGlass ;

	// This variable is used to signal a pseudo "middle button" mouse event
	// based on usage of the control or shift key...
	// This value should be user configurable?
	UINT MIDBUTTONFLAG ;

	// This variable keeps track of GraphEvent initiating a mouse move process 
	GraphEvent MBMode ; // Nonsense mouse button drag mode?

// Operations
protected:
	void DeleteImages() ;	// used in destructor

	static void createAceFont( CFont **font, int iFace, int iSize ) ;
	CPen *GetAcePen( int colour ) ;
	CBrush *GetAceBrush( int colour ) ;
	HCURSOR m_myCursor, m_otherCursor ;
	BOOL m_myCursorOn ;
	IMAGEMAP& GetImages() { return m_images ; }

public:
	static void graphColorInit() ;
	static void graphColorFinish() ;
	static void graphFontInit() ;
	static void graphFontFinish() ;
	
	// Reinitializes non-NULL fonts of given face and size
	static void CGraphView::reCreateFont( int iFace, int iSize ) ;

	// Overloaded functions to set and return the current colours
	// Resets the current GDI foreground and background colours...
	void setColors() ; // reinitializes the current FCol and BCol...
	void setColors(unsigned char fcol, unsigned char bcol) ;

	// Overloaded functions to set and return the current linewidth
	// 1.	Integer linewidth are assumed to be in GDI logical coordinates and used directly
	int setLineWidth( int linewidth = 0 ) // defaults to 0 == 1 pixel thick lines
	{
		int oldLineWidth = m_currentLineWidth ;
		m_currentLineWidth = linewidth ;
		setColors() ; // Need to reset the current pens to the specified linewidth
		return oldLineWidth ;
	}
	// 2.	Float linewidths are assumed to be ACEDB user coordinates for convertion into GDI logical coordinates
	// The dw(uToXrel(t)) formula is calculated explicitly to maintain floating point accuracy
	int setLineWidth( float lineWidth ) // defaults to current gActive linewidth
	{
		return( setLineWidth ((int)(lineWidth*gActive->xFac*gXAspect)) ) ;
	}

	int getLineWidth() const { return m_currentLineWidth ; }

	void setFCol ( unsigned char fcol = gActive->color  ) { m_currentFCol = fcol ; }
	unsigned char getFCol () const { return m_currentFCol ; }
	void setBCol ( unsigned char bcol = gActive->color ) { m_currentBCol = bcol ; }
	unsigned char getBCol () const { return m_currentBCol ; }

	void setBrush ( CBrush *brush ) { m_currentBrush = brush ; }
	CBrush *getBrush () const { return m_currentBrush ; }

	
	void setMyCursor(HCURSOR cursor) { m_myCursor = cursor ; } ;

	BOOL myCursorOn() const { return m_myCursorOn ; } ;
	void setMyCursorOn() ;
	void setMyCursorOff() ;

	// Since the "FILL_RECTANGLE" type of functionality
	// is reused several times in the code...
	void drawFillRect (CRect *rect) ;

	// The master drawing function!
	void drawBox ( Box box, CRect &clip ) ;

	BOOL GraphReSize, ScrollBoundsSet, EraseGraph ;

	CGraphWindow *GetGraphWindow() ; // Parent CGraphWindow (CFrameWnd) of ViewPort
	CGraph *GetGraphDev() ;			// returns associated CGraph "device" object

	CBitmap * rawPixelImage(char *data, int w, int h, int len) ;

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CGraphView)
	public:
	virtual void OnPrepareDC(CDC* pDC, CPrintInfo* pInfo = NULL);
	protected:
	virtual void OnDraw(CDC* pDC);      // overridden to draw this view
	virtual void OnInitialUpdate();     // first time after construct
	virtual void OnActivateView(BOOL bActivate, CView* pActivateView, CView* pDeactiveView);
	virtual void OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint);
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnPrint(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);
	//}}AFX_VIRTUAL

// Implementation
protected:
#if defined(NEW_WIN32_GRAPHS)
	void DeleteBoxViews() ;	// used in destructor
#endif

	virtual ~CGraphView();

	void mouseCall( GraphEvent GEType, CPoint point ) ;

	int m_NewX, m_NewY ;

	void ReSizeGraph() ;

	static void printInit() ;		// called by CGraphView::OnBeginPrint()
	static void printFinish() ;	// called by CGraphView::OnEndPrint()

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CGraphView)
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMenuSelect(UINT nItemID, UINT nFlags, HMENU hSysMenu);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnLButtonDblClk(UINT nFlags, CPoint point);
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnChar(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnEditCopy();
	afx_msg void OnEditPaste();
	afx_msg void OnMButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMButtonUp(UINT nFlags, CPoint point);
	afx_msg UINT OnNcHitTest(CPoint point);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#if defined(NEW_WIN32_GRAPHS)
// a subgraph is a CGraphView window?
#define CSubGraphWindow CGraphView
#define VIEWPTR  (CGraphView::gViewPort)  
/* VIEWPORT defined in cgraph.h for NEW_WIN32_GRAPHS? */
#else
#define VIEWPTR  ((CGraphView *)gSubDev)
#define VIEWPORT gSubDev
#endif

#endif // !defined GRAPHVIEW_H
 
