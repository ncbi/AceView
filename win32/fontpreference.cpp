// $Id: fontpreference.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $
// FontPreference.cpp : implementation file
//

#include "stdafx.h"
#include "WinAce.h"
#include "cgraph.h"
#include "FontPreference.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CFontPreferencePage property page

IMPLEMENT_DYNCREATE(CFontPreferencePage, CPropertyPage)

CFontPreferencePage::CFontPreferencePage() : CPropertyPage(CFontPreferencePage::IDD)
{
	//{{AFX_DATA_INIT(CFontPreferencePage)
	m_Bold_Face = _T("");
	m_Greek_Face = _T("");
	m_Italic_Face = _T("");
	m_Fixed_Width_Face = _T("");
	m_Plain_Face = _T("");
	m_Default_Font_Size = 0;
	m_Huge_Font_Size = 0;
	m_Large_Font_Size = 0;
	m_Medium_Font_Size = 0;
	m_Small_Font_Size = 0;
	m_Tiny_Font_Size = 0;
	//}}AFX_DATA_INIT
}

CFontPreferencePage::~CFontPreferencePage()
{
}

void CFontPreferencePage::DoDataExchange(CDataExchange* pDX)
{
	CPropertyPage::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CFontPreferencePage)
	DDX_Text(pDX, IDC_BOLD_FONT, m_Bold_Face);
	DDX_Text(pDX, IDC_GREEK_FONT, m_Greek_Face);
	DDX_Text(pDX, IDC_ITALIC_FONT, m_Italic_Face);
	DDX_Text(pDX, IDC_FIXED_WIDTH_FONT, m_Fixed_Width_Face);
	DDX_Text(pDX, IDC_PLAIN_FONT, m_Plain_Face);
	DDX_Text(pDX, IDC_DEFAULT_FONT_HEIGHT, m_Default_Font_Size);
	DDV_MinMaxUInt(pDX, m_Default_Font_Size, 0, 128);
	DDX_Text(pDX, IDC_HUGE_FONT_HEIGHT, m_Huge_Font_Size);
	DDV_MinMaxUInt(pDX, m_Huge_Font_Size, 0, 128);
	DDX_Text(pDX, IDC_LARGE_FONT_HEIGHT, m_Large_Font_Size);
	DDV_MinMaxUInt(pDX, m_Large_Font_Size, 0, 128);
	DDX_Text(pDX, IDC_MEDIUM_FONT_HEIGHT, m_Medium_Font_Size);
	DDV_MinMaxUInt(pDX, m_Medium_Font_Size, 0, 128);
	DDX_Text(pDX, IDC_SMALL_FONT_HEIGHT, m_Small_Font_Size);
	DDV_MinMaxUInt(pDX, m_Small_Font_Size, 0, 128);
	DDX_Text(pDX, IDC_TINY_FONT_HEIGHT, m_Tiny_Font_Size);
	DDV_MinMaxUInt(pDX, m_Tiny_Font_Size, 0, 128);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CFontPreferencePage, CPropertyPage)
	//{{AFX_MSG_MAP(CFontPreferencePage)
	ON_BN_CLICKED(IDC_NEW_PLAIN_FONT, OnNewPlainFont)
	ON_BN_CLICKED(IDC_NEW_FIXED_WIDTH_FONT, OnNewFixedWidthFont)
	ON_BN_CLICKED(IDC_NEW_BOLD_FONT, OnNewBoldFont)
	ON_BN_CLICKED(IDC_NEW_ITALIC_FONT, OnNewItalicFont)
	ON_BN_CLICKED(IDC_NEW_GREEK_FONT, OnNewGreekFont)
	ON_EN_CHANGE(IDC_DEFAULT_FONT_HEIGHT, OnChangeDefaultFontSize)
	ON_EN_CHANGE(IDC_HUGE_FONT_HEIGHT, OnChangeHugeFontSize)
	ON_EN_CHANGE(IDC_LARGE_FONT_HEIGHT, OnChangeLargeFontSize)
	ON_EN_CHANGE(IDC_MEDIUM_FONT_HEIGHT, OnChangeMediumFontSize)
	ON_EN_CHANGE(IDC_SMALL_FONT_HEIGHT, OnChangeSmallFontSize)
	ON_EN_CHANGE(IDC_TINY_FONT_HEIGHT, OnChangeTinyFontSize)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFontPreferencePage message handlers

BOOL CFontPreferencePage::OnInitDialog() 
{
	CPropertyPage::OnInitDialog();
	
	m_Plain_Face		= CGraphView::fontFace[Plain] ;
	m_Fixed_Width_Face	= CGraphView::fontFace[FixedWidth] ;
	m_Bold_Face			= CGraphView::fontFace[Bold] ;
	m_Italic_Face		= CGraphView::fontFace[Italic] ;
	m_Greek_Face		= CGraphView::fontFace[Greek] ;
	m_Default_Font_Size	= CGraphView::fontSize[DefaultSize] ;
	m_Tiny_Font_Size	= CGraphView::fontSize[Tiny] ;
	m_Small_Font_Size	= CGraphView::fontSize[Small] ;
	m_Medium_Font_Size	= CGraphView::fontSize[Medium] ;
	m_Large_Font_Size	= CGraphView::fontSize[Large] ;
	m_Huge_Font_Size	= CGraphView::fontSize[Huge] ;

	// No font faces or sizes changed yet...
	for(int iFace = Plain; iFace < NUM_FONTFACES ; iFace++ )
		faceChanged[iFace] =  FALSE ;

	for(int iSize = DefaultSize; iSize < NUM_FONTSIZES ; iSize++ )
		sizeChanged[iSize] =  FALSE ;

	return FALSE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}

void CFontPreferencePage::OnOK() 
{
	TRACE("\tCFontPreferencePage::OnOK()\n") ;

	CPropertyPage::OnOK();

	if( sizeChanged[DefaultSize] )
	{
		TRACE("\t\tResetting default font size to\t== %d\n",
				m_Default_Font_Size) ;
		CGraphView::fontSize[DefaultSize] = m_Default_Font_Size ;
	}
	if( sizeChanged[Tiny] )
	{
		TRACE("\t\tResetting tiny font size to\t== %d\n",
				m_Tiny_Font_Size) ;
		CGraphView::fontSize[Tiny] = m_Tiny_Font_Size ;
	}
	if( sizeChanged[Small] )
	{
		TRACE("\t\tResetting Small font size to\t== %d\n",
				m_Small_Font_Size) ;
		CGraphView::fontSize[Small] = m_Small_Font_Size ;
	}
	if( sizeChanged[Medium] )
	{
		TRACE("\t\tResetting Medium font size to\t== %d\n",
				m_Medium_Font_Size) ;
		CGraphView::fontSize[Medium] = m_Medium_Font_Size ;
	}
	if( sizeChanged[Large] )
	{
		TRACE("\t\tResetting Large font size to\t== %d\n",
				m_Large_Font_Size) ;
		CGraphView::fontSize[Large] = m_Large_Font_Size ;
	}
	if( sizeChanged[Huge] )
	{
		TRACE("\t\tResetting Huge font size to\t== %d\n",
				m_Huge_Font_Size) ;
		CGraphView::fontSize[Huge] = m_Huge_Font_Size ;
	}
	
	// For all faces and sizes of fonts
	for(int iFace = Plain; iFace < NUM_FONTFACES ; iFace++ )
		for(int iSize = DefaultSize; iSize < NUM_FONTSIZES ; iSize++ )
			// If the font specification has changed...
			if(	faceChanged[iFace] || sizeChanged[iSize] )
				CGraphView::reCreateFont( iFace, iSize ) ; // ...then, recreate the font
	
	CPropertyPage::OnOK();
}

#define FNTDLGMASK fntDlg.m_cf.Flags =  \
		((fntDlg.m_cf.Flags | CF_FIXEDPITCHONLY | CF_SCREENFONTS) & \
		 ~(CF_SHOWHELP | CF_APPLY | CF_EFFECTS | CF_NOSTYLESEL | CF_NOSIZESEL))

// Simplistic implementation of the
// user specification of a new font:
// Just look for the font face returned?
CString CFontPreferencePage::ShowFontDlg(int iFace) 
{
	graphOut(	"Note: Any font size or styles specified is \n"
				"are currently ignored in the font dialog.." ) ;

	CFontDialog fntDlg ;
	FNTDLGMASK ;

	if(fntDlg.DoModal() == IDOK)
	{
		// Assume that the font face has changed
		CGraphView::fontFace[iFace] = fntDlg.GetFaceName() ;
		faceChanged[iFace] = TRUE ;
	}
	return CGraphView::fontFace[iFace] ;
}

void CFontPreferencePage::OnNewPlainFont() 
{
	m_Plain_Face = ShowFontDlg(Plain) ;
	UpdateData(FALSE) ;
}

void CFontPreferencePage::OnNewFixedWidthFont() 
{
	m_Fixed_Width_Face = ShowFontDlg(FixedWidth) ;
	UpdateData(FALSE) ;
}

void CFontPreferencePage::OnNewBoldFont() 
{
	m_Bold_Face = ShowFontDlg(Bold) ;
	UpdateData(FALSE) ;
}

void CFontPreferencePage::OnNewItalicFont() 
{
	m_Italic_Face = ShowFontDlg(Italic) ;
	UpdateData(FALSE) ;
}

void CFontPreferencePage::OnNewGreekFont() 
{
	m_Greek_Face = ShowFontDlg(Greek) ;
	UpdateData(FALSE) ;
}

void CFontPreferencePage::OnChangeDefaultFontSize() 
{
	sizeChanged[DefaultSize] = TRUE ;
}

void CFontPreferencePage::OnChangeTinyFontSize() 
{
	sizeChanged[Tiny] = TRUE ;
}

void CFontPreferencePage::OnChangeSmallFontSize() 
{
	sizeChanged[Small] = TRUE ;
}

void CFontPreferencePage::OnChangeMediumFontSize() 
{
	sizeChanged[Medium] = TRUE ;
}

void CFontPreferencePage::OnChangeLargeFontSize() 
{
	sizeChanged[Large] = TRUE ;
}

void CFontPreferencePage::OnChangeHugeFontSize() 
{
	sizeChanged[Huge] = TRUE ;
}
 
