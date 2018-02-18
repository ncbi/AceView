/*  File: version.c
 *  Author: Ed Griffiths (edgrif@mrc-lmba.cam.ac.uk)
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
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm1.cnusc.fr
 *
 * Description: This code reports the version of the ACEDB code, this
 *              has the following implications:
 *
 *              - The code will only read databases that report the
 *                same major version from the super block held in the
 *                database.
 *              - The code will write databases that have in their
 *                super block the major/minor/update levels from this
 *                file.
 *              - The version is the overall version of the acedb code
 *                as a whole, not the version of any one library.
 *              - This is the version that most users will talk about
 *                and understand.
 *
 *              See aceversion.h for more description of the interface,
 *              and wdoc/version.html for further information on
 *              version numbers.
 *
 * Exported functions: See aceversion.h
 *
 * HISTORY:
 * Last edited: Aug 10 10:28 2001 (edgrif)
 * * Mar  5 09:36 1999 (edgrif): Removed "Version" from the version string,
 *              Jean and Richard don't like it.
 * * Jan 12 16:15 1999 (fw): version 4_8 is now used for development trunk after 
 *              release of 4_7b to the users. This means that executables
 *              carrying the 4_8 stamp are known to be our development code,
 *              whereas 4_7b code is a user-release.
 * * Dec 10 13:29 1998 (edgrif): Add comments to allow automatic extraction
 *              of version data from this file by a perl script.
 * * Dec  3 10:42 1998 (edgrif): Change around some calls/names to fit in with
 *              acelib naming convention & support new version reporting.
 * * Oct 27 15:26 1998 (edgrif): Replace malloc with messalloc, redo error messages,
 *              they replicated the messcrash macro.
 * *  (done by il at some time): added i to string_len for the /0 character and 
 *              sizeof(char) so that no memory errors on SGI/SUN. 
 * Created: Wed Apr 29 13:49:20 1998 (edgrif)
 *-------------------------------------------------------------------
 */

#include <string.h>
#include <stdio.h>
#include <wh/regular.h>
#include <wh/version.h>					    /* For put string macros etc. */
#include <wh/aceversion.h>




/* **************  EDIT THIS BIT TO CHANGE RELEASE/UPDATE etc.   ************************  */
/*                                                                                         */
/* The official release numbers/letters for the WHOLE of acedb, please see the file        */
/* wdoc/version.html for how these numbers should be changed.                              */
/*                                                                                         */
#if defined(ACEDB5)

/* These are not parsed by the perl script.                                  */
#define ACEDB_VERSION          5
#define ACEDB_RELEASE          0
#define ACEDB_UPDATE           ""
#define ACEDB_RELEASE_DIR "RELEASE.EXPERIMENTAL"

#else

/* THIS SECTION IS PARSED BY A wtools/aceGetVersion                          */
/* DO NOT ALTER THE BELOW COMMENT TAGS WITHOUT ALTERING THIS SCRIPT.         */
/*                                                                           */
/*          you can edit these bits to change the numbers but don't          */
/*                             change the layout.                            */
/*                                         |                                 */
/*                                         V                                 */
/* ACE_VERSION_START                                                         */
#define ACEDB_VERSION                      4
#define ACEDB_RELEASE                      7
#define ACEDB_UPDATE                      "ncbi"
/* ACE_VERSION_END                                                           */
/*                                                                           */
/* END OF PERL PARSED SECTION.                                               */


/* THIS SECTION IS PARSED BY A wtools/aceSetReleaseDir                       */
/* DO NOT ALTER THE BELOW COMMENT TAGS WITHOUT ALTERING THIS SCRIPT.         */
/*                                                                           */
/* ACE_RELEASE_DIR_START                                                     */
#define ACEDB_RELEASE_DIR "RELEASE.EXPERIMENTAL"
/* ACE_RELEASE_DIR_END                                                       */
/*                                                                           */
/* END OF PERL PARSED SECTION.                                               */

#endif
/*                                                                                         */
/* **************************************************************************************  */



/* Used in constructing version string.                                                    */
#define ACEDB_TITLE           "ACEDB"
#define ACEDB_LINKDATE_WORD   UT_COMPILE_PHRASE
#define ACEDB_SEPARATOR       "_"


/* this will only be the link date if this module replaces "linkdate.c" as the module      */
/* that is compiled EVERY time a relink is done of executables that need to display the    */
/* link date. __DATE__ is in the form  "Mon dd yyyy", __TIME__ is in the form  "hh:mm:ss"  */
/* both are ANSI standards.                                                                */
#define ACEDB_LINK_DATE  __DATE__ " " __TIME__


/* Any Unix executable that has this module linked in will display the text strings below  */
/* when anyone uses the 'what' command on the executable, but I don't know what the        */
/* equivalent is on macs and windows...at least on those machines the below stuff will     */
/* be embedded in the executable.                                                          */
/*                                                                                         */
char *copyright =
"@(#) \n"
"@(#) ------------------------------------------------------------------------\n"
"@(#) "ACEDB_TITLE" "UT_MAKESTRING(ACEDB_VERSION)"."UT_MAKESTRING(ACEDB_RELEASE)ACEDB_UPDATE",  "ACEDB_LINKDATE_WORD" "ACEDB_LINK_DATE" \n"
UT_COPYRIGHT()
"@(#) ------------------------------------------------------------------------\n"
"@(#) \n" ;


/* Interface functions:                                                                    */

/* Returns the version number.                                                             */
int aceGetVersion(void)
  {
  return(ACEDB_VERSION) ;
  }


/* Returns the release number.                                                             */
int aceGetRelease(void)
  {
  return(ACEDB_RELEASE) ;
  }


/* Returns the update letter.                                                              */
char *aceGetUpdate(void)
  {
  return(ACEDB_UPDATE) ;
  }

/* Returns the directory where the current code was built from. Normally this is undefined */
/* when developers do builds. When we do the monthly builds however, the     */
/* wtools/aceSetReleaseDir script is run against the monthly build copy of   */
/* this file to insert the build directory name in to ACEDB_RELEASE_DIR.     */
/* This enables us to find out exactly which monthly build a user is         */
/* compiling code from.                                                      */
char *aceGetReleaseDir(void)
  {
  return(ACEDB_RELEASE_DIR) ;
  }


/* Returns the link date.                                                                  */
char *aceGetLinkDate(void)
  {
  return(ACEDB_LINK_DATE) ;
  }


/* Returns a standard string with the link date.                                           */
char *aceGetLinkDateString(void)
  {
  static char *linkdate_string = ACEDB_LINKDATE_WORD " " ACEDB_LINK_DATE ;

  return(linkdate_string) ;
  }


/* Returns the new format version string for modules which need to show the acedb version. */
char *aceGetVersionString(void)
  {
  static char *version_string = 
    ACEDB_TITLE " "
    UT_MAKESTRING(ACEDB_VERSION) ACEDB_SEPARATOR UT_MAKESTRING(ACEDB_RELEASE) ACEDB_UPDATE ;

  return(version_string) ;
  }


/*****************************/
/* Checks whether user specified the "-version" command line option, and     */
/* prints out the version to stdout and exits. Note that we use a raw        */
/* printf here because messout may have been pointed to somewhere graphical  */
/* by xace etc. which we don't want.                                         */
/*                                                                           */
void checkForCmdlineVersion(int *argcp, const char **argv)
{
  if (getCmdLineOption (argcp, argv, "-version", NULL))
    {
      printf("%s,  build dir: %s, %s\n",
	     aceGetVersionString(), aceGetReleaseDir(), aceGetLinkDateString()) ;
      exit(EXIT_SUCCESS) ;
    }

  return ;
}
