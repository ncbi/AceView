/*  File: restriction.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  4 14:47 1998 (fw)
 * * Jul 23 14:59 1998 (edgrif): First I added this header block, then
 *      I removed the redeclaration of an fmap function in fmap.h
 * Created: Thu Jul 23 14:59:08 1998 (edgrif)
 *-------------------------------------------------------------------
 */

/* $Id: restriction.h,v 1.3 2007/11/29 05:31:47 mieg Exp $ */

#ifndef DEF_RESTRICTION
#define DEF_RESTRICTION

#include "regular.h"

typedef struct { int i ; int mark, type ;} Site ; /* see fmapsequence.c */
int siteOrder (const void *a, const void *b) ;

void dnacptMultipleMatch (Array sites, 
			  Array dna, Array protein, int frame,
			  Array colors, Stack s, Stack sname, 
			  int from, int length, 
			  BOOL amino, int maxError, int maxN, char *translationTable) ;
     
Stack dnacptAnalyseRestrictionBox (char *text, BOOL amino, Stack *snp) ;
void dnaRepaint (Array colors) ;

void dnacptFingerPrintCompute (Array dna, int from, int to,
			       Array colors, Array hind3, Array sau3a,
			       Array fp);

#endif
