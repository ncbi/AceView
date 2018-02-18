/*  File: flag.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1997
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: flag.h,v 1.2 2007/03/16 18:27:52 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr  3 12:12 1997 (rd)
 * Created: Sat Feb  1 12:53:00 1997 (rd)
 *-------------------------------------------------------------------
 */

#include "dict.h"

typedef unsigned int FLAGSET ;

#define NON_FLAG (~0)		/* all bits - same as in table.h */
#define MAXFLAG (sizeof(FLAGSET)*8)

FLAGSET flag (char *set, char *s) ;
char *flagName (char *set, FLAGSET flag) ;
DICT *flagSetDict (const char *set) ;
char *flagNameFromDict (DICT *dict, FLAGSET flag) ;

/******************** end of file ***********************/
 
 
 
