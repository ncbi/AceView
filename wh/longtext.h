/*  File: longtext.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: header for longtext.c module
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  2 17:03 1998 (fw)
 * Created: Wed Dec  2 17:01:49 1998 (fw)
 *-------------------------------------------------------------------
 */

#ifndef _LONGTEXT_H
#define _LONGTEXT_H

#include "acedb.h"

BOOL longTextDump (FILE *f, Stack buffer, KEY k);
BOOL javaDumpLongText (KEY key);
BOOL longTextParse (int level, KEY key);
KEYSET longGrep (char *template); /* remap to paper */
KEYSET longGrepPlain (char *template) ; /* returns longtexts */
int longTextFilter  (KEYSET ks, Stack s, BOOL reMap) ; /* filter a keyset of abstarcts */
BOOL bsMinilibDumpLongText (FILE *f, Stack buffer, KEY k);
BOOL longTextGrep (KEY key, char *template, Stack sWhat) ;


#endif /* _LONGTEXT_H */
