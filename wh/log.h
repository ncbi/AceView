/*  File: log.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 2001
 *-------------------------------------------------------------------
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
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Interface to acedb logging interface.
 * HISTORY:
 * Last edited: Mar 25 15:31 2001 (edgrif)
 * Created: Mon Feb 26 13:17:06 2001 (edgrif)
 * CVS info:   $Id: log.h,v 1.1.1.1 2002/07/19 20:23:21 sienkiew Exp $
 *-------------------------------------------------------------------
 */
#ifndef ACEDB_LOG_H
#define ACEDB_LOG_H

typedef enum _aceLogExitType {ACELOG_NORMAL_EXIT = 0, ACELOG_ABNORMAL_EXIT,
			      ACELOG_CRASH_EXIT, ACELOG_CLOSE_LOG} aceLogExitType ;

BOOL openAceLogFile(char *logfile_name, ACEOUT *logfile_out, BOOL new_session) ;
BOOL appendAceLogFile(ACEOUT logfile_out, char *text) ;
BOOL closeAceLogFile(ACEOUT logfile_out, aceLogExitType exit_type) ;
ACEOUT startNewAceLogFile(ACEOUT current_log, char *logfile_name,
			  char *log_msg, BOOL new_session) ;

#endif	/* ACEDB_LOG_H */
