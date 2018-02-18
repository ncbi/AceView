/*  File: winfilquery.cpp
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 *  WIN32 version of filquery.c module
 *
 * Description: filqueryopen()
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 6 01:40 1996 (rbrusk): bug fix for null filname's
 * *  Jan 6 20:00 1996 (rmb): enhanced dialog parameters 
 * Created: Tue Jul  7 2:25 1995 (rmb)
 *-------------------------------------------------------------------
 */

/* $Id: winfilquery.cpp,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */

#include "stdafx.h"
#include "win32.h"

extern "C"
{
#include "regular.h"
}

extern "C" FILE *filqueryopen ( char *dirname, char *filname, char *ending, char *spec, char *title )
{
	static CString CurrDir ;
	char Filter[25], filter1[] = " files", filter2[] = "*.", *pFilter = Filter ;
	char fnameBuf[MAXPATHLEN], *pfnameBuf = fnameBuf ;
	
	// First, get a filename via the standard WIN32 file dialog
	
   	CFileDialog dlg(TRUE,ending,dirname) ;  // does anything else need to be specified?

 lao:
	if(!dirname || !*dirname)
		dlg.m_ofn.lpstrInitialDir	= (!CurrDir.IsEmpty() ? (const char *)CurrDir : NULL) ;
	else
		dlg.m_ofn.lpstrInitialDir	= dirname ;

   	if(filname)
		strcpy(pfnameBuf,filname) ;
	else
		fnameBuf[0] = '\0' ;
	dlg.m_ofn.lpstrFile 		= pfnameBuf ;
   	dlg.m_ofn.lpstrDefExt 		= ending ;
   	dlg.m_ofn.lpstrTitle 		= title ;

	if( !ending || !*ending )
		pFilter = "ace files\0*.ace\0\0" ;
	else
	{
		strcpy(pFilter,ending) ;
		strcat(pFilter,filter1) ;
		pFilter += strlen(pFilter)+1 ;  // Advance just beyond first string
		strcpy(pFilter,filter2) ;
		strcat(pFilter,ending) ;
		pFilter += strlen(pFilter)+1 ;
		*pFilter = '\0' ; // Attach second null character to filter string
	}	
   	dlg.m_ofn.lpstrFilter = Filter ;

		switch( spec[0] )
	{
		case 'w': dlg.m_ofn.Flags |= (OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST); break ;
		case 'r': dlg.m_ofn.Flags |= (OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST);  break ;
		default: ;/* do nothing */
	}

   	if( dlg.DoModal() == IDOK )
	{

		// Otherwise, get the full pathname

   		CString pn = dlg.GetPathName();
		int pathLen = dlg.m_ofn.nFileOffset ;

		// Save the current directory for future reference
		if(pathLen > 3)
			CurrDir = pn.Left(pathLen-1) ; // Remove one backslash?
		else if(pathLen == 3)  // else if root directory, keep backslash
			CurrDir = pn.Left(pathLen) ;
		/* else, something weird is going on; ignore? */

		// Then, open the file if possible

  		FILE*	fil = 0 ;

		if (  spec[0] == 'w' && (fil = fopen (pn, "r")) )
		{
			if ( fil != stdin && fil != stdout && fil != stderr)
				fclose (fil) ; 
			fil = 0 ;
			if (messQuery (messprintf ("Overwrite %s?",pn ) ) )
			{ 
				if (fil = fopen (pn, spec) )
					goto bravo ;
				else
					messout ("Sorry, can't open file %s for writing", pn ) ;
			}
			goto lao ;
		}
		else if (!(fil = fopen ( pn, spec)))
		{
			messout ("Sorry, can't open file %s", pn ) ;
			goto lao ;
		}
	bravo:
		return fil ;
	}
	else
		return 0 ; // failure

}

/********************* end of file ********************/
 
 
