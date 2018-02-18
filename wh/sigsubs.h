/*  File: sigsubs.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1999
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
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: public header for signal handler module sigsubs.c
 * Exported functions:
 * HISTORY:
 * Last edited: Mar  7 10:46 2001 (edgrif)
 * Created: Fri Apr 30 13:41:09 1999 (fw)
 *-------------------------------------------------------------------
 */

#ifndef ACEDB_SIGSUBS_H
#define ACEDB_SIGSUBS_H

/* Cmd line option to turn off signal handling.                              */
#define NOSIGCATCH_OPT "-nosigcatch"

/* Standard typedef for signal handler function, why isn't this in the       */
/* POSIX header ??                                                           */
typedef void (*SignalHandler)(int) ;


/* Functions to control acedb signal handling.                               */

/* Initialise and turn on/off signal catching in the process.                */
void signalCatchInit(BOOL initCtrlC, BOOL no_catch) ;	    /* in some cases you may not
							       want to register the standard
							       Ctrl-C sig handler */

void signalCatchOn(BOOL initCtrlC) ;
void signalCatchOff(BOOL initCtrlC) ;



/*                 Utility functions.                                        */

/* Initialise a POSIX signal handler struct to empty, except for the signal  */
/* function.                                                                 */
void signalInitSigAction(struct sigaction *sigact, SignalHandler sig_func) ;

/* Set action for a signal back to the system default so that the signal is  */
/* not caught by the process.                                                */
BOOL signalSetNoCatch(int sig) ;



#endif	/* ACEDB_SIGSUBS_H */
