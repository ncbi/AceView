
/*  File: parsedata.h
 *  Author: Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 *	Ancillary class for WIN32 multithreading functionality for winace et al. 
 *  See win32thread.cpp for implementation specifics...
 *
 * HISTORY:
 * Last edited:Apr 4 03:30 1997 (rbrusk)
 * Created: Apr 4 03:30 1997 (rbrusk)
  *-------------------------------------------------------------------
 */

/* $Id: parsedata.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */

class CParseData : public CObject
{
	DECLARE_DYNAMIC(CParseData)

public:

// Attributes

	CWinThread *m_pParseThread ;
	FILE *m_pFile ;
	int m_isSkip ;
	KEYSET m_ks ;
	CMutex pHold ;

// Operations
	CParseData(FILE *pFile) ;
} ;


 
