/*  File: dna.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
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
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jun 12 13:26 2003 (edgrif)
 * * May 24 10:26 1999 (edgrif): Add some missing extern func decs.
 * * Jul 23 14:41 1998 (edgrif): Remove redeclarations of fmap functions.
 * * Jul 23 09:06 1998 (edgrif): I have added this header, which was
 *      missing (that's why it shows me as author). I have removed the
 *      declaration of fMapFindDNA which is now in fmaps public header.
 * Created: Thu Jul 23 09:06:15 1998 (edgrif)
 *-------------------------------------------------------------------
 */

/* $Id: dna.h,v 1.12 2017/07/27 00:04:14 mieg Exp $ */

#ifndef DEF_DNA_H
#define DEF_DNA_H

#include "acedna.h"
#include "keyset.h"
/*
On disk, we store the dna as 4 bases per byte, if there is no ambiguity, or
2 bases per byte otherwise.
*/
Array dnaGet (KEY key) ;
Array dnaHandleGet (KEY key, AC_HANDLE h) ;
Array dnaGetWithErrors (KEY key) ;

int dnaAdd (Array a, KEY seq, int start, int stop, BOOL noMismatch) ; /* 0: mismatch, 1: absent, 2: success */
BOOL dnaStoreDestroy (KEY key, Array dna) ; 
			/* beware - sets length in sequence object, and allways destroys the dna array  */
void dnaExonsSort (Array a) ;	/* a is of BSunit from flattening Source_exons */
void dnaExonsSort3 (Array a, int nn) ; /*a variant of above */
BOOL dnaSubClass (KEY seq, KEY *dnaKeyp) ;   /* creates DNA name xxx:yyy, i.e. className(seq):name(seq) except for ?Sequence */
BOOL dnaReClass (KEY dna, KEY *seqp) ;  /* gets the parent class of dna, by default a Sequence */
BOOL dnaInClass (int classe) ; /* does this class supports dna ? */

				/* Handles for dnacpt package */
void dnaAnalyse (void) ;

                    /* export the usage tables of all 12 mers */
Array dnaGetWordUsage (Array kSet, int step, unsigned long *nbp, unsigned long *nwp, 
		       BOOL create, const char *fName, int comb) ;

				/* FastA dumping routines */
BOOL dnaDumpFastAKey (KEY key, FILE *fil, Stack s, char style) ;
BOOL dnaDumpFastAKeyWithErrors (KEY key, FILE *fil, Stack s, char style) ;
BOOL dnaZoneDumpFastAKey (KEY key, FILE *fil, Stack s, char style, int x1, int x2, BOOL withMissmatch, BOOL noClassNam) ;
int dnaDumpFastAKeySet (KEYSET kSet, FILE *fil, Stack s) ;
				/* if fil==0 these call dnaFileOpen */
FILE *dnaFileOpen (void) ;
void saucisseTest (KEYSET ks) ;

/*                  NO CODE AFTER THIS                                       */
#endif
