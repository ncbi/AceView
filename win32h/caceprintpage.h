/* $Id: caceprintpage.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */
/*  File: CAcePrintPage.h
 *  Author: Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *          adapted from code written by Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Win32 Printing Page Data Structure header file
 * Exported functions:  CAcePrintPage() constructor,
						getPageLength(), getPageWidth(), getMetric() ;
 * HISTORY:
 * Created: Apr 21 21:00 1996 (rbrusk)
 *-------------------------------------------------------------------
 */
class CAcePrintPage : public CObject
{
private:
	int m_pageLength, m_pageWidth ;
	int m_XRes, m_YRes ;  // made float in calcPageSize() for metric conversion
	BOOL m_isMetric ;
	BOOL calcPageSize(	int pageSize, float XDPI,float YDPI,
						int *pageWidth, int *pageLength) ;

public:
	CAcePrintPage( LPDEVMODE printer ) ;
	int PageLength() const { return m_pageLength ; } ;
	int PageWidth() const { return m_pageWidth ; } ;
	BOOL isMetric() const  { return m_isMetric ; } ;
} ;

 
