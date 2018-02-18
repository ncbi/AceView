/* $Id: fontpreference.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
// FontPreference.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CFontPreferencePage dialog

class CFontPreferencePage : public CPropertyPage
{
	DECLARE_DYNCREATE(CFontPreferencePage)

private:
	CString ShowFontDlg(int iFace) ;

	BOOL	faceChanged[NUM_FONTFACES],
			sizeChanged[NUM_FONTSIZES] ;

// Construction
public:
	CFontPreferencePage();
	~CFontPreferencePage();

// Dialog Data
	//{{AFX_DATA(CFontPreferencePage)
	enum { IDD = IDD_FONTS };
	CString	m_Bold_Face;
	CString	m_Greek_Face;
	CString	m_Italic_Face;
	CString	m_Fixed_Width_Face;
	CString	m_Plain_Face;
	UINT	m_Default_Font_Size;
	UINT	m_Huge_Font_Size;
	UINT	m_Large_Font_Size;
	UINT	m_Medium_Font_Size;
	UINT	m_Small_Font_Size;
	UINT	m_Tiny_Font_Size;
	//}}AFX_DATA


// Overrides
	// ClassWizard generate virtual function overrides
	//{{AFX_VIRTUAL(CFontPreferencePage)
	public:
	virtual void OnOK();
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	// Generated message map functions
	//{{AFX_MSG(CFontPreferencePage)
	virtual BOOL OnInitDialog();
	afx_msg void OnNewPlainFont();
	afx_msg void OnNewFixedWidthFont();
	afx_msg void OnNewBoldFont();
	afx_msg void OnNewItalicFont();
	afx_msg void OnNewGreekFont();
	afx_msg void OnChangeDefaultFontSize();
	afx_msg void OnChangeHugeFontSize();
	afx_msg void OnChangeLargeFontSize();
	afx_msg void OnChangeMediumFontSize();
	afx_msg void OnChangeSmallFontSize();
	afx_msg void OnChangeTinyFontSize();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

};
 
