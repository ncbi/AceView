/*  File: command.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
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
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: header file for functions in command.c
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 24 22:22 2001 (edgrif)
 * * Nov 23 09:46 1999 (edgrif): Insert new choixmenu values for admin cmds.
 * Created: Thu Aug 26 17:06:53 1999 (fw)
 * CVS info:   $Id: command.h,v 1.8 2009/02/25 15:19:22 mieg Exp $
 *-------------------------------------------------------------------
 */

#ifndef ACEDB_COMMAND_H
#define ACEDB_COMMAND_H

#include "acedb.h"
#include "query.h"		/* for COND struct */
#include "spread.h"		/* for SPREAD struct */
#include "dict.h"


/* Opaque reference to an instance of an ace command.                        */
typedef struct _AceCommandStruct *AceCommand ; 

/************************************************************/

AceCommand aceCommandDoCreate (int choix, FILE *fil, Stack s);
void aceCommandSwitchActiveSet (AceCommand look, KEYSET ks, KEYSET ks2) ;
AceCommand aceCommandCreate (int choix, int toto) ;
void aceCommandDestroy (AceCommand look);

void commandExecute (int level, BOOL noSubshell, 
		     BOOL localIsInteractive, 
		     FILE *fil, Stack s, int choix, KEYSET activeKs) ;

KEY aceCommandDoExecute (AceCommand look, int level, 
			 int localIsInteractive, 
			 KEY option, int maxChar);

Stack commandStackExecute (AceCommand look, const char *command) ; /* no fioritures */

void aceCommandNoWrite (AceCommand look); /* drop write access */

/******************/

void tStatus (void) ;

extern void set_user_prompt (const char *prompt);


#endif
