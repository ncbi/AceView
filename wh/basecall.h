/*  File: basecall.h
 *  Author: Jean Thierry-Mieg
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
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 15 14:38 1999 (fw)
 * Created: Mon Dec 14 15:40:24 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: basecall.h,v 1.5 2017/08/07 20:09:55 mieg Exp $ */

#ifndef _BASECALL_H
#define _BASECALL_H

#include "graph.h"
#include "map.h"

#include "../wstaden/Read.h"   /* Staden's */

#include "dna.h"
#include "../whooks/classes.h"
#include "bs.h"

#include "../whooks/tags.h"

#include "bs.h"
#include "../whooks/systags.h"
#include "../whooks/sysclass.h"

#include "lex.h"
#include "query.h"
#include "dnaalign.h"
#include "dna.h"
#include "cdna.h"


#define A_COLOR GREEN
#define T_COLOR RED
#define G_COLOR DARKGRAY /* lightgray in base call explicit */
#define C_COLOR CYAN
#define N_COLOR YELLOW

#define BC_LOW 0x10
#define BC_HAND 0x20
#define BC_ACE 0x40
#define BC_TAG 0x80

#define NF 32
/************************************************************/

typedef struct LaneStruct
  { KEY key, dnaKey ;
    Read*   seq ;
    int x1, x2, x3;  /* coordonnes de base clipTop, end, extend wiz wDna */
    int   clipTop, clipEnd, clipExtend ;      /* clipping in base coordinates */
    int   previousEnd, nextTop ;  /* for cdna splicing */
    int   vectorTop, vectorEnd, qualityEnd ;
    int   t1, t2 ;		/* zone shown in trace coordinates */
    int   seqClipTop, seqClipEnd ; /* clipping in base in unedited sequence */
    int   xClipTop, xClipEnd, xClipExtend; /* clipping in basePos[] coordinates */
    int   handClipTop, handClipEnd;   /* clipping in basePos[] coordinates */
    int scf ;  /* 0: closed, 1: dna failed, 2: dna ok (traceGetLane)
		  3: trace failed, 4: scf trace, 5: ctftrace (baseCallGetSeq) */
    int hide, lastSpliceHint, getNewSpliceHint, showTraceBeforeClipTop ;
    int jackPotOriginLong, jackPotOriginShort, jackPotJump, jackPotSens, jackPotMatchLength ;
    BOOL upSequence , isTouched, isAligned, isShown ;
    int color ;
    int  hasTag ; /* 0: unknown, 1: absent, 2: present */
    int laneMin, laneMax ;  /* trace height range */
    float dy, ddy ; /* local offset, relative mag */
    Array errArray ; 
    int nerr ; /* number of errors shown in traceEdLanes */
    Array dna; 
    float xHideButton ;
    int minBaseBox, minBase, maxBase ; /* in curve region */
    int edLaneBaseBox, laneBaseCallBox, greenBox ; 
    int boxMoreClip, boxLessClip,boxMoreClipTop, boxLessClipTop, boxName, boxClose, boxNewSpliceHint ;
    float edLaneBaseBoxOffSet, laneBaseCallBoxOffSet ;
    int laneBaseCallBoxMin, laneBaseCallBoxMax ;
    Array base, basePos, baseQuality ; /* the edited ABI base call in ace format */
    KEY baseCallKey ; /* the stored base call + pos in BASECALL format */
    KEY basePositionKey ; /* the stored position in char (dx) format */
    KEY baseQualityKey ; /* the stored quality */
    Array baseCall ; /* the acedb base call in BASECALL */
    int maxPos ;
    Array tags ; /* BSunits depth 4 */
    KEY clone ;
    int favorDelete ;
    float traceRenormalize ;
    int traceShift ;
    BOOL isPoor ;
  } LANE ;

/************** LOOK structure **************************/

struct LookStruct
{
    void*   magic;		/* == &TRACELOOK_MAGIC */
    Graph graph ;
    MAP   map ;
            /* concensus dna */
    int graphHeight, graphWidth, topMargin ; /* geometry of the graph */
    KEY key, dnaKey, link ; /* key is the contig */
    Array dna, dnaR ;
    int mode ; /* of enum AcemblyMode if > 0 dont edit the concensus */
    int   wmin, wc, wmax ;	/* min, centre  and max wDNA pos on screen */
    int sens ; /* orientation du contig dans le link */
    int edMin, edMax ;  /* min max in the editor space (left) */
    float edMag ;   /* relative mag of the editor */
    int hide, hideUp ; /* what we want to show */
    int next ; /* nextwhat */
    char *summary ;
    int
      activeBox, minLiveBox,
      minWildBase, minWildBox, maxWildBox, summaryBox ;

    Array box2seg ;		/* SEG if >minLiveBox, <maxLiveBox */
    Associator traceAss ;
    Array 
      lanes , /* developped traces */
      baseBoxes,
      laneShown,
      hints ;
    Array hits ; /* cdna splicing */
    KEY gene ;   /* cdna splicing */
    LANE *activeLane ;
    int traceWidth, origin ;
    KEYSET modifiedTraces ;
    Graph fMapGraph ; void* fMapLook ;
    BOOL modified, showClip, consensus ;
    BOOL coord ; /* to show the coord of the individual traces */
    int geneMin, geneMax, geneBox,  hideDots ;
    BOOL select ;
    KEYSET selectedKs ; /* to select just a few traces */
  }  ;

/* struct segStruct { KEY key ;
		   int iShort, iLong ;
		 } seg ;
*/
typedef struct hintStructure 
 { KEY key ; LANE *lane ; int iShort, iLong, isSpliceHint ; float y ; } HINT ;

typedef struct baseCallStruct
 { int x ; /* the position in trace coordinates */
   short y ; /* the height of the trace there */
   char t ; /* 0,1,2,3 for A,G,C,T */
   char flag ;
 } BASECALL ;


int seqEnergyOfDerivee (Read* seq, int min, int max) ;
typedef struct Extremum { int x ; short y ; } Extremum ; 
void baseCallStore (LANE *lane) ;
BOOL baseCallGet (LANE *lane) ;
BOOL baseCallGetSeq (LANE *lane) ;
BOOL baseCallFakePos (LANE *lane) ; /* to use trace.c when traces not available */

BOOL baseCorrel(Read* seq1, int x1, BOOL direct,
	       Read* seq2, int x2, int ll, int nstep, int step, int *bestdxp) ;

int baseCallPatchLane (Array dnaDirect, Array dnaReverse, LANE *lane) ;
BOOL baseCallPatchContig (KEY contig, int *nnp) ; /* nnp += result, you must init *nnp to zero */
BOOL baseCallUnclipLane (LANE *lane, KEY type, int *dxp) ;
int baseCallUnclipKeySet (KEY link, KEYSET ks, KEY type, int *dxp) ;
int baseCallClipContig2Max (KEY link, int max, int *dxp) ;

void laneEditSave (LOOK look, LANE *lane, int pos, int n) ;
void  laneMakeErrArray (LOOK look, LANE *lane) ;
void laneDestroy (LANE *) ;
void findXclipping (LANE *lane) ;

BOOL traceGetLane(LOOK look, LANE *lane) ;
void trackContig (LOOK look, int z1, int z2, BOOL whole) ;
void monDnaForget(void *look, KEY key) ;
void lanePlotDerivee (LOOK look, LANE *lane) ;
void lanePlotPeriodicite (LOOK look, LANE *lane) ;
void baseCallRedoLaneBaseCall (LANE *lane, int *np) ;
void nnLaneTrain (LOOK look, LANE *lane) ;
void nnAssemblyTrain (LOOK look) ;
void nnQualityTrain (LANE *lane, int pos, int rate) ;
char nnBaseCall (LANE *lane, int pos) ;
BOOL findBaseCall (LANE *lane) ;
char* baseCallNameGuess (KEY key, int type) ;
BOOL baseCallFlagRnaEditing (KEY tg, int *nreadp, int *nagp, int *nagrp) ;
BOOL baseCallRedoBaseCall (KEY key, int *np) ;
#endif /* _BASECALL_H */
