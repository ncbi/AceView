/* $Id: main.h,v 1.3 2015/08/17 16:12:51 mieg Exp $ */
/*  File: main.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: include for mainpick.c
     only for graphical builds of acedb
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 11 15:13 1998 (fw)
 * Created: Wed Nov 18 16:52:27 1998 (fw)
 *-------------------------------------------------------------------
 */

#ifndef _MAIN_H
#define _MAIN_H

#include "graph.h"		/* for Graph type */

void helpOnClass (void);
int getCurrentClass ();

BOOL  checkWriteAccess(void);

void pickGetArgs(int *argcp, char **argv);
void pickDraw (void);	      /* also called when write access changes */
void pickReReadClasses (void);/* called when models change */
Graph pickCreate (void);
void pickPopMain(void);

void mainPickReport (int n) ;	/* used by mainpick.c & xclient.c */

void messStatusDisplayInWindow(const char * text);

void quitFunctionRegister (VoidRoutine func);


void classListCreate (Array *classList);

int ksetClassComplete (char *text, int len, int classe);
void mainKeySetComplete (KEYSET ks, char *text, int len);

#endif /* _MAIN_H */
