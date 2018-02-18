/*  File: fmap.h
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Public header file for fmap = DNA sequence package
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  8 11:27 1999 (edgrif)
 * * Jan  8 11:27 1999 (edgrif): Added missing func dec for fMapGifAlign.
 * * Nov 19 09:11 1998 (edgrif): TINT_<colour> defs should be unsigned.
 * * Jul 27 10:38 1998 (edgrif): Add final declarations of external
 *      fmap routines.
 * * Jul 23 09:09 1998 (edgrif): The original fmap has mostly now gone
 *      into the private header fmap_.h. This header now includes only
 *      fmap functions that are used externally to fmap.
 * Created: Sat Jul 25 20:28:37 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: fmap.h,v 1.5 2017/02/15 20:39:28 mieg Exp $ */

#ifndef ACEDB_FMAP_H
#define ACEDB_FMAP_H
/*                  NO CODE BEFORE THIS                                      */

#include "acedb.h"				  /* necessary include files */
#include "graph.h"
#include "bs.h"
#include "map.h"

/* Colours used for display by fMap, placed here for use by packages         */
/* cooperating with fMap.                                                    */
/* (internal note: tint values must match fMapTints array in fmapsequence.c) */
/* You will probably assign these to a char, make sure its declared as       */
/* unsigned.                                                                 */
#define TINT_HIGHLIGHT1 0x01U			  /* highest priority, dna highlighting.*/
#define TINT_HIGHLIGHT2 0x02U			  /* highlight friends.*/
#define TINT_RED        0x04U 
#define TINT_LIGHTGRAY  0x08U 
#define TINT_MAGENTA    0x10U 
#define TINT_CYAN       0x20U 
#define TINT_LIGHTGREEN 0x40U 
#define TINT_YELLOW     0x80U 

/************************************************************/
/* global LOOK when no graph for gif package */

extern LOOK fMapGifLook;

/************************************************************/

/* Public fmap functions.                                                    */
/*                                                                           */
int fmapView (KEY view, KEY tag) ; /* controls the fmap display, mieg, 2004 */
BOOL fMapDisplay (KEY key, KEY from, BOOL isOldGraph) ;
Graph fMapActiveGraph(void) ;
BOOL fMapActive(Array *dnap, Array *dnaColp, KEY *seqKeyp, void** lookp) ;
BOOL fMapFindZone(void *vv, int *min, int *max, int *origin) ;
void fMapDraw(LOOK look, KEY from) ;
void fMapReDrawDNA(void *look) ;
Array fMapFindDNA(KEY seq, int *start, int *stop) ;

/* Used by gifacemain currently.                                             */
void* fMapGifGet (void*) ;	/* returns a handle used by other fMapGif*() */
void fMapGifDisplay (void*) ;
void fMapGifFeatures (void*) ;
void fMapGifActions (void*) ;
void fMapGifDNA (void*) ;
void fMapGifColumns (void*) ;
void fMapGifDestroy (void*) ;
void fMapGifAlign (void* opaqueLook) ;

int fMapQueryColor (KEY key) ;

/* Used by gfcode.c & hexcode.c                                              */
void fMapAddGfSite (int type, int pos, float score, BOOL comp) ;
void fMapAddGfCodingSeg (int pos1, int pos2, float score, BOOL comp) ;
void fMapAddGfSeg (int type, KEY key, int x1, int x2, float score) ;

/* Used by map stuff in w7 at least.                                         */
void fMapShow3FrameTranslation (LOOK look, float *offset) ;
void fMapShowGeneTranslation (LOOK look, float *offset) ;
/* mhmp 09.09.98 */
void fMapShowUpGeneTranslation (LOOK look, float *offset) ;/*mhmp 17.06.98*/
void fMapShowDownGeneTranslation (LOOK look, float *offset) ;
void fMapShowORF (LOOK look, float *offset) ;

/* Used in a number of fmap routines + in geldisp.c                          */
void fMapRegisterSites(void *v, Array sites, Stack sname) ;

/* Used in plotseq.c, trace.c and other wabi programs...                     */
void fMapClearDNA (LOOK look) ;
void fMapPleaseRecompute (void *look) ;
void fMapDrawFromTrace (LOOK look, KEY from, int xx, int type) ;
void fMapBoxInfo (LOOK look, int box, void *seg) ;

/*                  NO CODE AFTER THIS                                       */
#endif /* ACEDB_FMAP_H */
 
