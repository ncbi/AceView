
/*  File: win32lib.cpp
 *  Author: Richard Bruskiewich (rbrusk@octogene.medgen.ubc.ca)
 *-------------------------------------------------------------------
 * This file is part of the Windows NT/'95 port of the
 *  ACEDB genome database package
 *
 *  ACEDB was largely written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:	a general library of useful Win32 routines
 *   
 * Exported functions: CommandLineToArgv(), CleanUpApp(), 
 			DosToPosix(), getINIParameter(),
 *			LoadWinAceINIParameters() and SaveWinAceINIParameters(),
 *			dbgPos()
 * HISTORY:
 * Last edited: Jun 12 17:08 1996 (rbrusk):
 * * JIT saving of session parameters to the registry
 * * Jun 12 17:08 1996 (rbrusk): dos2unix(), unix2dos() fix
 * * Jun 3 17:48 1996 (rbrusk): dos2unix(), unix2dos()
 * * May 20 23:15 1996 (rbrusk): setSessionParameters()
 * * May 9 10:05 1996 (rbrusk): acedb_editor() stub implementation 
 * * April 17 11:32 1996 (rbrusk): dbgPos() implemented
 * * Jan 4 18:20 1996 (rmb): Converted to .cpp
 * 		-	Imported or created DosToPosix(), getINIParameter(),
 *			LoadWinAceINIParameters() and SaveWinAceINIParameters() 
 * * Jun 22 02:24 1995 (rmb): CommandLineToArgv() 
 * Created: Jun 8 21:30 1995 (rmb)
 *-------------------------------------------------------------------
 */

/* $Id: win32lib.cpp,v 1.1.1.1 2002/07/19 20:23:25 sienkiew Exp $ */

#include "stdafx.h"
#include "win32.h"

#include "WinAce.h" /* Also graph_.h=>graph.h=>array.h=>regular.h etc. */

/* #undef getenv() macro WIN32 #define'd in mystdlib.h
   so that you can invoke the real O/S call down below */
#undef getenv

#include <tchar.h>
#include <wtypes.h>
#include <malloc.h>

#include "Preferences.h"

//*************************************************************************
// In UNIX ACEDB, "getenv" is used extensively to obtain user set program
// parameters.  For WIN32, an "initialization" file is also provided;
// In addition, command line switches may provide additional information.
// Finally, it is conceivable to provide a "Preferences" setting command option
// during the execution of the program, to give the user power to set and reset
// program options. It seems reasonable to merge all these parameter setting
// mechanisms under one umbrella... This module.
//*************************************************************************
static CMapStringToString iniMap, envMap ;

static BOOL envParameter(CString param)
{
	int idx ;
	if( (idx = param.Find('=')) <= 0 )
		return FALSE ;
	TRY
	{
		// 	else, treat as an environmental parameter spec
		CString label = param.SpanExcluding("=") ;
		CString value =	param.Mid(++idx) ;
		envMap.SetAt(label,value) ;
	}
	CATCH(CMemoryException,e)
	{
	   messcrash("envParameter() memory exception?");
	}
	END_CATCH
	return TRUE ;
}

#define FIRST_ARG 1  /* if invoking program name not argv[0],
                        then this symbol should be set to 1! */
static char defaultProgName[] = "WinAce" ;

static char ** _argvlocal = 0 ;  // Keep track of argv messalloc() for memory management

char **CommandLineToArgv(LPSTR lpCmdLine, int * pArgc)
{
	int i;
	char ** argv;
	char *argvbuf[80]; /* That's a linefull of args!! */
	char * cp;

	/* Check for null arguments  */
	if(	!(cp = (char *)lpCmdLine) )
	{
		messcrash("*** Fatal Error: NULL lpCmdLine ptr in WinAceMain:CommandLineToArgv()? ***");
	}

	if( FIRST_ARG ) argvbuf[0] = defaultProgName ;

	/* Scan the command line string */
    argvbuf[*pArgc = FIRST_ARG] = cp;
	while( *cp )
	{
       /* Skip over leading spaces */
       while( *cp == ' ') ++cp;
	   if(*cp == '\0') break ; /* EOL? */
	   /* 	Else, was not a blank so increment argument count
	   		and set argv to *argument */
	   argvbuf[(*pArgc)++]=cp;

       /* Skip over argument to next blank or to EOL delimiter */
	   while( *cp != ' ' /* Blank */ && *cp != '\0' ) ++cp;
	   if(*cp == '\0') break ;  /* EOL encountered? */
	   /*else, space found? Terminate string and continue */
	   *cp++ = '\0'	;
    }
 	argvbuf[*pArgc] = NULL ; // NULL ptr

	// Now, scan arguments for local environment parameter settings
	// The only criterion is the presence of an "equal" sign embedded
	// in the argument string.  Such strings are extracted and removed
	// from the argv list, if they correspond to valid parameter settings
	int s, d ; d = s = 1;
	while( s < *pArgc ) 
	{
		argvbuf[d] = argvbuf[s] ;
		if( !envParameter(argvbuf[s++]) ) d++ ;
	}
	// reset pArgc
	argvbuf[*pArgc = d] = NULL ;

	// Can't use messalloc() here because system not set up at this point?
	if( (_argvlocal = argv = (char **)messalloc((*pArgc)*(sizeof(char *))) ) == NULL )
	  {
		   messcrash("WinAce()::CommandLineToArgv() malloc() error?");
	  }
    else
       for( i=0; i < *pArgc ; ++i) argv[i] = argvbuf[i];  // Build list of pointers

    return argv;
}

// saveSessionParameters() saves current program parameters, as recorded in the iniMap,
// into the external WIN32 O/S program registry. NOTE: THIS IS A SELECTIVE SAVE
// Certain parameters, (e.g."username") are not saved, for practical reasons
// Such parameters are filtered out by the locally defined "Saveable() function
static BOOL Saveable(CString &key)
{
	if( !key.CompareNoCase("USERNAME") )
		return FALSE ;
	else
		return TRUE ;
}

void saveSessionParameters()
{
	CString key,value ;
	POSITION pos = iniMap.GetStartPosition() ;
	while( pos != NULL )
	{
		iniMap.GetNextAssoc(pos,key,value) ;
		if( Saveable(key) )
			theApp.WriteProfileString( "WinAce", key, value) ;
	}				              		
}

extern "C" void setProgParameter(const char *eStr )
{
	ASSERT(eStr != NULL) ;

	CString src , key, value ;
	TRY
	{
		src		= eStr ;
		key		= src.SpanExcluding("=")  ;
		value	= src.Mid(src.Find('=')+1) ;

		if(value.IsEmpty()) 
		{
			iniMap.RemoveKey(key) ;
			// Need to remove old keys in the O/S Registry
			// to ensure that old registry database value
			// isn't seen in a subsequent getProgParameter() access?
			theApp.WriteProfileString( "WinAce", key, NULL) ;
			return ;
		}

		iniMap.SetAt(key,value) ;
		if( Saveable(key) )
			theApp.WriteProfileString( "WinAce", key, value) ;
	}
	CATCH(CMemoryException,e)
	{
		src.Empty() ;
		key.Empty() ;
		value.Empty() ;
		messcrash("WinAce()::setProgParameter() memory exception?");
	}
	END_CATCH
}

// This routine retrieves the string value associated with a
// program parameter identifier (such as ACEDB_DATA)
// The routine does a sequential search in a number of places...
// First, it searches the current local key map of parameters 
// previously accessed or changed (this list starts out empty)
// Then, looks for commandline overrides of parameters
// Next, looks for set environment variables
// Finally, looks in the program registry for the variables
// A side effect is the changing of the .ini file (via iniMap)
extern "C" char *getProgParameter(const char *key)
{
	TRACE("\n*** Entering getINIParameter(key == %s) *** ",key) ;
	static CString INIStr ;
	CString keyStr = key ;
	union jester { const char *yin ; char *yang ; } fool ;
	
	INIStr.Empty() ;
	keyStr.MakeUpper() ; // search not case sensitive? 

	// First, look in the iniMap containing modified
	// or previously accessed .ini values
	if( iniMap.Lookup(keyStr,INIStr) )
		goto Success ;

	// Then, look in the commandline override envMap values
	if( envMap.Lookup(keyStr,INIStr) )
		goto Success ;

	//  Next, look for a set OS environment variable
	if( !(INIStr = getenv((const char *)keyStr)).IsEmpty())
		goto Success ;

	// Finally, look in the .ini file
	INIStr = theApp.GetProfileString( "WinAce", keyStr ) ;
	if( INIStr.IsEmpty() )
	{
		TRACE("keyStr(%s) parameter is NULL?\n",
				(const char *)keyStr ) ;
		return NULL ;
	}

Success:

	// Record the current setting in the internal registry as a 
	// side effect of accessing parameter, for future reference.
	// This will result in the environment parameter being reset in the
	// external WIN32 O/S Registry Database, at session termination
	TRACE("Setting keyStr(%s) == INIStr(%s)\n",
		   (const char *)keyStr, (const char *)INIStr ) ;
	setProgParameter( (LPCTSTR)(keyStr+"="+INIStr) ) ;
	fool.yin = (const char *)INIStr ;
	return fool.yang ;

}
extern "C" void getSessionParameters() /* see session.c */
{
    TRACE("\n*** Entering getSessionParameters():\n") ;
	// Use only the database profile page of the Preferences dialog?
	CPreferences Preferences("Database Profile",DBPROFILEPAGE) ;

	if( Preferences.DoModal() ==  IDOK)
	{
		TRACE("\n*** Returning IDOK from Preferences Dialog in getSessionParameters()\n") ;
		// Since getSessionParameters() is called from within
		// sessionInit() before the database is fully initialized
		// it should be safe to simply return if ReloadACEDB is indicated
		// and let sessionInit() Retry its sessionFilname() initialization
		if( Preferences.ReloadACEDB() ) return ;
	}
	// If the user did not update the database profile
	// then I am no better off than before, so I crash...
	messcrash("You must tell me where the database is. Until then, I die!") ;
}

extern "C" void acedb_editor(char *tmpfile)
{
	messout( "WIN32 acedb_editor not yet implemented") ;
}

// These two functions interconvert strings between
// multi-line text strings with null termination and
// EOL conventions for dos and unix respectively
// that is, <cr><lf> (DOS) to/from <lf> only (UNIX)
static CStringList duHeap ;
CString &dos2unix(const CString &source)
{
	CString *target;
	TRY
	{
		target = new CString ;

		BOOL crFlag = FALSE ; // flag for deferred copying of \r's
		for(int i; i < source.GetLength(); i++ )
		{
			if(crFlag)
			{
				if( source[i] != '\n' )
					*target += '\r' ;
				crFlag = FALSE ;
			}
			if( source[i] == '\r' )
			{
				crFlag = TRUE ;
				continue ;
			}
			*target += source[i] ;
		}
		duHeap.AddHead( (LPCTSTR)(*target) ) ;
	}
	CATCH(CMemoryException,e)
	{
		messcrash("dos2unix(): Memory exception") ;
	}
	END_CATCH
	return *target ;
}

CString &unix2dos(const CString &source)
{
	CString *target;
	TRY
	{
		target = new CString ;

		for(int i; i < source.GetLength(); i++ )
			if( source[i] == '\n' )
				*target += "\r\n" ;
			else
				*target += source[i] ;
		duHeap.AddHead( (LPCTSTR)(*target) ) ;
	}
	CATCH(CMemoryException,e)
	{
		messcrash("unix2dos(): Memory exception") ;
	}
	END_CATCH
	return *target ;
}

/**************** End of File ******************/

 
 
 
 
 
