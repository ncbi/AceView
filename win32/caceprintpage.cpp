
/*  File: CAcePrintPage.cpp
 *  Author: Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *          adapted from code written by Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Win32 Printing Page Data Structure
 * Exported functions:  CAcePrintPage() constructor,
						pageLength(), pageWidth(), isMetric() ;
 * HISTORY:
 * Created: Apr 20 21:00 1996 (rbrusk)
 *-------------------------------------------------------------------
 */

/* $Id: caceprintpage.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $ */
#include "stdafx.h"
#include "caceprintpage.h"

// This function calculates pageLength and pageWidth with the 
// expected page sizes normalized into printer dots (from dpi, via dots per 0.1 mm)
// if pageSize is NULL, then the pageLength and pageWidth arguments are returned (scaled to dots)
// function returns "TRUE" for a "metric" page size, FALSE for imperial sizes
BOOL CAcePrintPage::calcPageSize(int pageSize, float YDPI, float XDPI, 
								 int *pageWidth,int *pageLength)
{
	if( YDPI > 0 ) // a real dpi
		YDPI *= (1.0 / 254.0) ; // rescale dpi to dots per 0.1 mm
	else
		YDPI = 1.0 ;  // unscaled (or should I use "default" scalings?)

	if( XDPI > 0 ) // a real dpi
		XDPI *= (1.0 / 254.0) ; // rescale dpi to dots per 0.1 mm
	else
		XDPI = YDPI ;  // YDPI is isotropic resolution dpi

	if(!pageSize)
	// use pageWidth & pageLength (converted from 1/10 mm to dot values)
	{
		*pageLength = (int)(*pageLength*YDPI) ; 
		*pageWidth = (int)(*pageWidth*XDPI) ; 
		return TRUE ;
	}

	switch( pageSize ) // as defined for DEVMODE dmPaperSize
	{
		case DMPAPER_LEGAL:		// Legal, 8 1/2 by 14 inches
			*pageLength	= (int)(14*254*YDPI) ;
			*pageWidth	= (int)(8.5*254*XDPI) ;
			break ;
		case DMPAPER_CSHEET:	// C Sheet, 17 by 22 inches
			*pageLength	= (int)(22*254*YDPI) ;
			*pageWidth	= (int)(17*254*XDPI) ;
			break ;
		case DMPAPER_DSHEET:	// D Sheet, 22 by 34 inches
			*pageLength	= (int)(34*254*YDPI) ;
			*pageWidth	= (int)(22*254*XDPI) ;
			break ;
		case DMPAPER_ESHEET:	// E Sheet, 34 by 44 inches
			*pageLength	= (int)(44*254*YDPI) ;
			*pageWidth	= (int)(34*254*XDPI) ;
			break ;
		case DMPAPER_TABLOID:	// Tabloid, 11 by 17 inches
		case DMPAPER_11X17:		// 11-by-17-inch sheet
			*pageLength	= (int)(17*254*YDPI) ;
			*pageWidth	= (int)(11*254*XDPI) ;
			break ;
		case DMPAPER_LEDGER:	// Ledger, 17 by 11 inches
			*pageLength	= (int)(11*254*YDPI) ;
			*pageWidth	= (int)(17*254*XDPI) ;
			break ;
		case DMPAPER_STATEMENT:	// Statement, 5 1/2 by 8 1/2 inches
			*pageLength	= (int)(8.5*254*YDPI) ;
			*pageWidth	= (int)(5.5*254*XDPI) ;
			break ;
		case DMPAPER_EXECUTIVE:	// Executive, 7 1/4 by 10 1/2 inches
			*pageLength	= (int)(10.5*254*YDPI) ;
			*pageWidth	= (int)(7.25*254*XDPI) ;
			break ;
		case DMPAPER_A3:		// A3 sheet, 297 by 420 millimeters
			*pageLength	= (int)(4200*YDPI) ;
			*pageWidth	= (int)(2970*XDPI) ;
			break ;
		case DMPAPER_A4:		// A4 Sheet, 210 by 297 millimeters
		case DMPAPER_A4SMALL:	// A4 small sheet, 210 by 297 millimeters
			*pageLength	= (int)(2970*YDPI) ;
			*pageWidth	= (int)(2100*XDPI) ;
			break ;
		case DMPAPER_A5:		// A5 sheet, 148 by 210 millimeters
			*pageLength	= (int)(2100*YDPI) ;
			*pageWidth	= (int)(1480*XDPI) ;
			break ;
		case DMPAPER_B4:		// B4 sheet, 250 by 354 millimeters
			*pageLength	= (int)(3540*YDPI) ;
			*pageWidth	= (int)(2500*XDPI) ;
			break ;
		case DMPAPER_B5:		// B5 sheet, 182-by-257-millimeter paper
			*pageLength	= (int)(2570*YDPI) ;
			*pageWidth	= (int)(1820*XDPI) ;
			break ;
		case DMPAPER_QUARTO:	// Quarto, 215-by-275-millimeter paper
			*pageLength	= (int)(2750*YDPI) ;
			*pageWidth	= (int)(2150*XDPI) ;
			break ;
		case DMPAPER_10X14:		// 10-by-14-inch sheet
			*pageLength	= (int)(14*254*YDPI) ;
			*pageWidth	= (int)(10*254*XDPI) ;
			break ;
		case DMPAPER_ENV_9:		// #9 Envelope, 3 7/8 by 8 7/8 inches
			*pageLength	= (int)(8.875*254*YDPI) ;
			*pageWidth	= (int)(3.875*254*XDPI) ;
			break ;
		case DMPAPER_ENV_10:	// #10 Envelope, 4 1/8 by 9 1/2 inches
			*pageLength	= (int)(9.5*254*YDPI) ;
			*pageWidth	= (int)(4.125*254*XDPI) ;
			break ;
		case DMPAPER_ENV_11:	// #11 Envelope, 4 1/2 by 10 3/8 inches
			*pageLength	= (int)(10.375*254*YDPI) ;
			*pageWidth	= (int)(4.5*254*XDPI) ;
			break ;
		case DMPAPER_ENV_12:	// #12 Envelope, 4 3/4 by 11 inches
			*pageLength	= (int)(11*254*YDPI) ;
			*pageWidth	= (int)(4.75*254*XDPI) ;
			break ;
		case DMPAPER_ENV_14:	// #14 Envelope, 5 by 11 1/2 inches
			*pageLength	= (int)(11.5*254*YDPI) ;
			*pageWidth	= (int)(5*254*XDPI) ;
			break ;
		case DMPAPER_ENV_DL:	// DL Envelope, 110 by 220 millimeters
			*pageLength	= (int)(2200*YDPI) ;
			*pageWidth	= (int)(1100*XDPI) ;
			break ;
		case DMPAPER_ENV_C5:	// C5 Envelope, 162 by 229 millimeters
			*pageLength	= (int)(2290*YDPI) ;
			*pageWidth	= (int)(1620*XDPI) ;
			break ;
		case DMPAPER_ENV_C3:	// C3 Envelope,  324 by 458 millimeters
			*pageLength	= (int)(4580*YDPI) ;
			*pageWidth	= (int)(3240*XDPI) ;
			break ;
		case DMPAPER_ENV_C4:	// C4 Envelope,  229 by 324 millimeters
			*pageLength	= (int)(3240*YDPI) ;
			*pageWidth	= (int)(2290*XDPI) ;
			break ;
		case DMPAPER_ENV_C6:	// C6 Envelope,  114 by 162 millimeters
			*pageLength	= (int)(1620*YDPI) ;
			*pageWidth	= (int)(1140*XDPI) ;
			break ;
		case DMPAPER_ENV_C65:	// C65 Envelope, 114 by 229 millimeters
			*pageLength	= (int)(2290*YDPI) ;
			*pageWidth	= (int)(1140*XDPI) ;
			break ;
		case DMPAPER_ENV_B4:	// B4 Envelope,  250 by 353 millimeters
			*pageLength	= (int)(3530*YDPI) ;
			*pageWidth	= (int)(2500*XDPI) ;
			break ;
		case DMPAPER_ENV_B5:	// B5 Envelope,  176 by 250 millimeters
			*pageLength	= (int)(2500*YDPI) ;
			*pageWidth	= (int)(1760*XDPI) ;
			break ;
		case DMPAPER_ENV_B6:	// B6 Envelope,  176 by 125 millimeters
			*pageLength	= (int)(1250*YDPI) ;
			*pageWidth	= (int)(1760*XDPI) ;
			break ;
		case DMPAPER_ENV_ITALY: // Italy Envelope, 110 by 230 millimeters
			*pageLength	= (int)(2300*YDPI) ;
			*pageWidth	= (int)(1100*XDPI) ;
			break ;
		case DMPAPER_ENV_MONARCH:	// Monarch Envelope, 3 7/8 by 7 1/2 inches
			*pageLength	= (int)(7.5*254*YDPI) ;
			*pageWidth	= (int)(3.875*254*XDPI) ;
			break ;
		case DMPAPER_ENV_PERSONAL:	// 6 3/4 Envelope, 3 5/8 by 6 1/2 inches
			*pageLength	= (int)(6.5*254*YDPI) ;
			*pageWidth	= (int)(3.625*254*XDPI) ;
			break ;
		case DMPAPER_FANFOLD_US:	// US Std Fanfold, 14 7/8 by 11 inches
			*pageLength	= (int)(14.875*254*YDPI) ;
			*pageWidth	= (int)(11*254*XDPI) ;
			break ;
		case DMPAPER_FANFOLD_STD_GERMAN:	// German Std Fanfold, 8 1/2 by 12 inches
			*pageLength	= (int)(12*254*YDPI) ;
			*pageWidth	= (int)(8.5*254*XDPI) ;
			break ;
		case DMPAPER_FOLIO:		// Folio, 8-1/2-by-13-inch paper
		case DMPAPER_FANFOLD_LGL_GERMAN:	// German Legal Fanfold, 8 1/2 by 13 inches
			*pageLength	= (int)(13*254*YDPI) ;
			*pageWidth	= (int)(8.5*254*XDPI) ;
			break ;
		case DMPAPER_LETTER:	// Letter, 8 1/2 by 11 inches
		case DMPAPER_LETTERSMALL:	// Letter Small, 8 1/2 by 11 inches
		case DMPAPER_NOTE:		// Note, 8 1/2 by 11 inches
		default: 	// DMPAPER_LETTER:Letter, 8 1/2 by 11 inches
			*pageLength	= (int)(11*254*YDPI) ;
			*pageWidth	= (int)(8.5*254*XDPI) ;
		}
		switch( pageSize ) // as defined for DEVMODE dmPaperSize
	{
		case DMPAPER_A4:		// A4 Sheet, 210 by 297 millimeters
		case DMPAPER_A3:		// A3 sheet, 297 by 420 millimeters
		case DMPAPER_A4SMALL:	// A4 small sheet, 210 by 297 millimeters
		case DMPAPER_A5:		// A5 sheet, 148 by 210 millimeters
		case DMPAPER_B4:		// B4 sheet, 250 by 354 millimeters
		case DMPAPER_B5:		// B5 sheet, 182-by-257-millimeter paper
		case DMPAPER_QUARTO:	// Quarto, 215-by-275-millimeter paper
		case DMPAPER_ENV_DL:	// DL Envelope, 110 by 220 millimeters
		case DMPAPER_ENV_C5:	// C5 Envelope, 162 by 229 millimeters
		case DMPAPER_ENV_C3:	// C3 Envelope,  324 by 458 millimeters
		case DMPAPER_ENV_C4:	// C4 Envelope,  229 by 324 millimeters
		case DMPAPER_ENV_C6:	// C6 Envelope,  114 by 162 millimeters
		case DMPAPER_ENV_C65:	// C65 Envelope, 114 by 229 millimeters
		case DMPAPER_ENV_B4:	// B4 Envelope,  250 by 353 millimeters
		case DMPAPER_ENV_B5:	// B5 Envelope,  176 by 250 millimeters
		case DMPAPER_ENV_B6:	// B6 Envelope,  176 by 125 millimeters
		case DMPAPER_ENV_ITALY: // Italy Envelope, 110 by 230 millimeters
			return TRUE ; // metric page size

		case DMPAPER_LETTER:	// Letter, 8 1/2 by 11 inches
		case DMPAPER_LEGAL:		// Legal, 8 1/2 by 14 inches
		case DMPAPER_CSHEET:	// C Sheet, 17 by 22 inches
		case DMPAPER_DSHEET:	// D Sheet, 22 by 34 inches
		case DMPAPER_ESHEET:	// E Sheet, 34 by 44 inches
		case DMPAPER_LETTERSMALL:	// Letter Small, 8 1/2 by 11 inches
		case DMPAPER_TABLOID:	// Tabloid, 11 by 17 inches
		case DMPAPER_LEDGER:	// Ledger, 17 by 11 inches
		case DMPAPER_STATEMENT:	// Statement, 5 1/2 by 8 1/2 inches
		case DMPAPER_EXECUTIVE:	// Executive, 7 1/4 by 10 1/2 inches
		case DMPAPER_FOLIO:		// Folio, 8-1/2-by-13-inch paper
		case DMPAPER_10X14:		// 10-by-14-inch sheet
		case DMPAPER_11X17:		// 11-by-17-inch sheet
		case DMPAPER_NOTE:		// Note, 8 1/2 by 11 inches
		case DMPAPER_ENV_9:		// #9 Envelope, 3 7/8 by 8 7/8 inches
		case DMPAPER_ENV_10:	// #10 Envelope, 4 1/8 by 9 1/2 inches
		case DMPAPER_ENV_11:	// #11 Envelope, 4 1/2 by 10 3/8 inches
		case DMPAPER_ENV_12:	// #12 Envelope, 4 3/4 by 11 inches
		case DMPAPER_ENV_14:	// #14 Envelope, 5 by 11 1/2 inches
		case DMPAPER_ENV_MONARCH:	// Monarch Envelope, 3 7/8 by 7 1/2 inches
		case DMPAPER_ENV_PERSONAL:	// 6 3/4 Envelope, 3 5/8 by 6 1/2 inches
		case DMPAPER_FANFOLD_US:	// US Std Fanfold, 14 7/8 by 11 inches
		case DMPAPER_FANFOLD_STD_GERMAN:	// German Std Fanfold, 8 1/2 by 12 inches
		case DMPAPER_FANFOLD_LGL_GERMAN:	// German Legal Fanfold, 8 1/2 by 13 inches
		default: 	// DMPAPER_LETTER:Letter, 8 1/2 by 11 inches
			return FALSE ; // Imperial page size
	}
}

// 8 1/2" x 11" sheet default page 
#define DEFAULTRESOLUTION 300   // dpi
#define DEFAULTWIDTH (int)(8.5*254)
#define DEFAULTLENGTH (int)(11*254)

CAcePrintPage::CAcePrintPage( LPDEVMODE printer )
{
	int paperSize =
		(printer->dmFields & DM_PAPERSIZE)?
			printer->dmPaperSize : 0 ;

	m_XRes =
		(printer->dmFields & DM_PRINTQUALITY)?
			printer->dmPrintQuality : DEFAULTRESOLUTION ;

	m_YRes =
		(printer->dmFields & DM_YRESOLUTION)?
			printer->dmYResolution : 0 ;

	m_pageWidth =
		(printer->dmFields & DM_PAPERWIDTH)?
			printer->dmPaperWidth : DEFAULTWIDTH ;

	m_pageLength =
		(printer->dmFields & DM_PAPERLENGTH)?
			printer->dmPaperLength : DEFAULTLENGTH ;

	m_isMetric = calcPageSize(paperSize, (float)m_XRes,(float)m_YRes, &m_pageWidth, &m_pageLength) ;

	if(	(printer->dmFields & DM_ORIENTATION) && 
		printer->dmOrientation == DMORIENT_LANDSCAPE)
	{
		int temp = m_pageLength ;
		m_pageLength = m_pageWidth ; 
		m_pageWidth = temp ;
	}
}
 
 
