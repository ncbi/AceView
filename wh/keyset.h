/*  File: keyset.h
 *  Author: R Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: public header for keyset operations.
 *              This file is part of acedb.h and NOT to be included
 *              by other source files.
 *              The KEYSET operations are built upon the Array ops
 *              provided by the utilities library libfree.a
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 11 09:44 1998 (fw)
 * Created: Fri Dec 11 09:42:41 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: keyset.h,v 1.7 2015/03/11 20:20:06 mieg Exp $ */

#ifndef DEFINE_KEYSET_H
#define DEFINE_KEYSET_H

#include "regular.h"		/* header for libfree.a */

/***************************************************************/
/*  a KEYSET is an ordered array of KEYs.                         */
/***************************************************************/

typedef Array KEYSET ;  /* really KEYSET = array(,,KEY) always ordered */
#define keySetCreate()		arrayCreate(32,KEY)
#define keySetHandleCreate(h)	arrayHandleCreate(32,KEY,h)
#define keySetReCreate(s)	arrayReCreate(s,32,KEY)
#define keySet(s,i)		array(s,i,KEY)
#define keySetDestroy(s)	arrayDestroy(s)

#define keySetInsert(s,k)	arrayInsert(s,&(k),keySetOrder)
#define keySetRemove(s,k)	arrayRemove(s,&(k),keySetOrder)
#define keySetSort(s)		arraySort((s),keySetOrder) 
#define keySetCompress(s)	arrayCompress(s)
#define keySetFind(s,k,ip)	arrayFind ((s),&(k),(ip),keySetOrder)
#define keySetMax(s)		arrayMax(s)
#define keySetExists(s)		(arrayExists(s) && (s)->size == sizeof(KEY))
#define keySetCopy(s)		arrayCopy(s)
#define keySetHandleCopy(s,h)		arrayHandleCopy(s,h)

unsigned int keySetNCompress (KEYSET ks, unsigned int min) ; /* keep keys repeated min times */
KEYSET  keySetAND (KEYSET x, KEYSET y) ;
KEYSET  keySetOR (KEYSET x, KEYSET y) ;
KEYSET  keySetXOR (KEYSET x, KEYSET y) ;
KEYSET  keySetMINUS (KEYSET x, KEYSET y) ;
int     keySetOrder (const void *a, const void*b) ;
int     keySetAlphaOrder (const void *a, const void*b) ;
KEYSET  keySetHeap (KEYSET source, unsigned int nn, int (*order)(KEY *, KEY *)) ;
KEYSET  keySetNeighbours (KEYSET ks) ;
KEYSET  keySetAlphaHeap (KEYSET ks, unsigned int nn) ;    /* jumps aliases/deletes */
KEYSET  keySetAlphaHeapAll (KEYSET ks, unsigned int nn) ; /* do not jump aliases */
unsigned int	keySetCountVisible (KEYSET ks) ;
unsigned int     keySetPartition (Array kss) ;  /* input: kss is an Array of KEYSET
					    *  output: all pairwise intersects are empty 
					    */

/**************************************************************/

BOOL keySetActive (KEYSET *setp, void** lookp) ;
void keySetSelect () ;
BOOL keySetDump(FILE *f, Stack buffer, KEYSET s);
KEYSET keySetRead (int level) ; /* construct a keyset with recognised keys */
void keySetKill (KEYSET ks) ; /* kill and lexkill all members of the keyset */

#endif
/*************************************************************/


