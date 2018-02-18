/*  File: parse.h
 *  Author: Jonathan Hodgkin (cgc@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm1.cnusc.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  9 15:40 1998 (fw)
 * Created: Mon Nov 22 14:03:14 1993 (cgc)
 *-------------------------------------------------------------------
 */

/* $Id: parse.h,v 1.5 2005/09/12 03:17:37 mieg Exp $ */

#include "acedb.h"

	/* keyset ks is filled with list of parsed objects if non-zero */

BOOL parseFile (FILE *fil, int lineBox, KEYSET ks) ;
BOOL parseNamedFile (char *fileName, FILE *fileHandle, int newLineBox, KEYSET ks) ;
BOOL parsePipe (FILE *fil, int lineBox, KEYSET ks) ; /* will pclose */
void parseOneFile (FILE *fil, KEYSET ks) ;
BOOL parseBuffer (char *text, KEYSET ks) ;
/* convenience call private to AceC */
BOOL parseAceCBuffer (const char *text, Stack s, void *ace_out, KEYSET ks) ;
BOOL parseLevel (int level, KEYSET ks) ;
int parseErrorsLastParse;     /* hack to sneak a number out of deep nested calls in parse code */


extern char parseLineText[64] ;
extern BOOL overRideParseProtection;
extern BOOL parseKeepGoing;

/********* end of file ***********/
