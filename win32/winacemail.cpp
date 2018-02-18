
/*  File: winacemail.cpp
 *  Author: R. Bruskiewich
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: Implements a CMC compliant email interface for WinAce
 * HISTORY:
 * Last edited: 
 * Created: Oct 13 11:20 1996 (rbrusk)
 *-------------------------------------------------------------------
 */
/* $Id: winacemail.cpp,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */

#include "stdafx.h"
#include "win32.h"

extern "C" {
#include "regular.h"
}

#include <xcmc.h>

extern "C" void win32Mail(char *args)
{
	TRACE("Entering win32Mail with args == %s\n", args ) ;
	messout( "WIN32 email functionality not yet implemented") ;
}

 
