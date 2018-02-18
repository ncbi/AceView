/*  File: dbpath.h
 *  Author: Simon Kelley (srk@sanger.ac.uk)
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
 *      Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *      Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 *-------------------------------------------------------------------
 */
#ifndef ACEDB_DBPATH_H
#define ACEDB_DBPATH_H

BOOL dbPathInit(char *dbpath_from_command_line);

/* Return file in either $ACEDB/dir or $ACEDB_COMMON/dir,
   checks for existance of file.
   It doesn't make sense to have a "w" spec to this call. */
char *dbPathFilName (char *dir, char *name, char *ending, char *spec, 
		     AC_HANDLE handle);

/* Return file in $ACEDB/dir, checks for read access ("r" spec) or 
   a writable directory ("w" spec). */
char *dbPathStrictFilName (char *dir, char *name, char *ending, char *spec, 
			   AC_HANDLE handle);

/* Make fileName in $ACEDB/dir, no checking is done, the filename 
   is made by simple string concatentation of $ACEDB/dir/name.ending */
char *dbPathMakeFilName (char *dir, char *name, char *ending, 
			 AC_HANDLE handle);

/* All these functions return memory alloced on the handle: it is the 
   callers responsibility to free it */

/* The scripts to run exteral programs are different under windows, and I want
   want the same database to be cross-mountable to both windows and unix
   so the scripts under windows are in "winscripts"

   Set a macro here to determine the name of the script directory.
   there should be no direct references to wscripts in the code */

#ifdef __CYGWIN__
#define WSCRIPTS "winscripts"
#else
#define WSCRIPTS "wscripts"
#endif
 


#endif /* ACEDB_DBPATH_H */ 
