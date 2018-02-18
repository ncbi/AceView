/*  File: aqlerror.c
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: %W% %G%
 * Description: error handling within AQL
 * Exported functions:
 * HISTORY:
 * Last edited: Sep  2 14:20 2004 (edgrif)
 * * Aug  5 10:58 1998 (fw): only aqlError left, which now reports
                             errors to the strings inside the AQL object
 * * Jul 14 12:00 1998 (fw): error message reporting now chops
                             up the message string line by line
 * Created: Wed Jul 22 13:01:56 1998 (fw)
 *-------------------------------------------------------------------
 */
/**********************************************************************/

#include <wh/acedb.h>
#include <wh/aceio.h>
#include <waql/aql_.h>

/***************************************************************/
/************************ error reporting **********************/
/***************************************************************/

void aqlError (AQL aql, int errNo, int errPos, char *format, ...)
/* 
   Called from anywhere in the aql query processing code
   for FATAL errors, i.e. when the processing can't proceed.
   It jumps out of the code and into the error reporting.
   A unique error number is used for every invalid ocurrence,
   the position is usually node->pos, and the message is
   a format-string with arguments just as in printf() 

   NOTE, this function could do nasty array-bounds-overwrites with
         the error strings, if the allocated space isn't big enough.
*/
{
  va_list args ;
  char buf[22] ;
  int start, end ;
  char *cp, *startLine;

  if (!aql)
    messcrash("aqlerror.c:aqlError() - called with NULL pointer");

  if (errNo == 0)		/* execution interrupted */
    {
      /* jump back to the calling code */
      longjmp (*(aql->errorJmpEnv), errNo) ;
    }

  /* report error to calling code  */
  aql->IsError = TRUE;
  aql->errorNumber = errNo;
  aql->errorMessage = (char*)halloc (1000, aql->handle);
  aql->errorReport = (char*)halloc (2000, aql->handle);

  va_start (args,format) ;
  vsprintf (aql->errorMessage, format, args) ;
  va_end (args) ;

  if (errPos)
    { 
      /* pretty error reporting, if the error position is known (i.e. non-ZERO) */
      start = errPos - 9 ; 
      if (start < 0) 
	start = 0 ;

      end = errPos + 9 ; 
      if (end >= strlen (aql->query->text))
	end = strlen(aql->query->text) - 1 ;

      strncpy (buf, aql->query->text + start, end-start+1) ;
      buf[end-start+1] = '\0' ;


      sprintf (aql->errorReport,
	       "// AQL error %3d around: '%s'\n"
	       "//                        %*s\n", errNo, buf, errPos-start,"^") ;
    }
  else
    {
      sprintf (aql->errorReport,
	       "// AQL error %d: \n", errNo) ;
    }

  /* chop up the message string nicely into lines starting with *** */
  cp = startLine = aql->errorMessage;

  while (*cp)
    {
      if (*cp == '\n')
	{
	  strcat  (aql->errorReport, "// ");
	  strncat (aql->errorReport, startLine, cp-startLine);
	  strcat  (aql->errorReport, "\n");

	  startLine = cp+1;
	}
      cp++;
    }
  /* append last line when cp reached '\0' */
  strcat (aql->errorReport, "// ");
  strcat (aql->errorReport, startLine);
  strcat (aql->errorReport, "\n");

  /* jump back to the calling code */
  longjmp (*(aql->errorJmpEnv), errNo) ;


  /* We should never get here. */
  return ;
} /* aqlError */

/************************** EOF **********************************/

