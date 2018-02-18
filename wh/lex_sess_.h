/*  File: lex_sess_.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: header for functions shared between
 *              the session manager and the lexer
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 23 18:02 1998 (fw)
 * Created: Mon Nov 23 17:59:01 1998 (fw)
 *-------------------------------------------------------------------
 */

#ifndef _LEX_SESS__H
#define _LEX_SESS__H

void lexSessionStart (void);
void lexSessionEnd (void);	/* track new/touched keys */
void lexOverLoad(KEY key, DISK disk);
void lexReadGlobalTables(void);
void lexClear(void);

#endif /*  _LEX_SESS__H */

/***************************** eof ****************************/
