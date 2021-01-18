/*  File: acembly.h
 *  Author: Danielle et jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@frmop11.bitnet
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 13 23:04 1996 (mieg)
 * Created: Fri Oct 20 20:18:19 1995 (mieg)
 *-------------------------------------------------------------------
 */
/* $Id: acembly.h,v 1.4 2017/07/27 00:04:14 mieg Exp $ */
 
             /*********************************************/
             /* Acembly.h                                   */
             /* type definitions and size limits          */
             /*********************************************/
 
#ifndef DEF_Acembly_h
#define DEF_Acembly_h
 
#include "acedb.h"
#include "query.h"

enum  AcemblyMode { SHOTGUN = 0, MUTANT, CDNA} ;
extern int acemblyMode ;

extern KEY _Clips,
 _Later_part_of, _Hand_Clipping, 
_OriginalBaseCall, _Align_to_SCF, _SCF_Position,_Quality ;

Array blyDnaGet (KEY link, KEY contig, KEY mykey) ;
Array trackErrors (Array  dna1, int pos1, int *pp1, 
		   Array dna2, int pos2, int *pp2, int *NNp, Array err, int mode) ;
Array baseCallCptErreur (Array dnaD, Array dnaR, Array dnaCourt, BOOL isUp,
			 int x1, int x2, int x3, 
			 int *clipTopp, int *cEndp, int *cExp, int mode) ;

void autoAnnotateGene (KEY gene, KEY annot) ;

void kantorProduct (char * keyname) ;
void ficheGraph (KEY gene) ;
int ficheAceRun (char * keyname,void * funccall,void * callback,int params) ;
extern WEBQUERYFUNC owqQuery ;
BOOL queryRegisterWebQuery (WEBQUERYFUNC f) ;

void traceGraphDestroy (void) ; /* defined in trace.c */
void defComputeTace (int level, KEYSET ks) ;
int abiFixDoubleContig (KEY link, KEY contig, int *ip);
void abiFixFinish (void) ;
void abiFixDoFinish (KEY link) ;
int abiFixDoubleContig (KEY link, KEY contig, int *ip);
int abiFixDouble (KEY link, KEYSET ks, int *ip) ;
int abiFixDoubleContig (KEY link, KEY contig, int *ip) ;
void abiFixExtendContig (KEY link, KEY contig) ;
void abiFixExtend (KEY link, KEYSET ks) ;

void defCptOpen (KEY link) ;
BOOL baseCallPatchContig (KEY contig, int *nnp);
int baseCallUnclipContig (KEY contig, KEY type, int *dxp) ;
int baseCallTileContig (KEY contig, KEY type, int *dxp) ;
KEY lastAssembly (void) ;
int baseCallUnclipKeySet (KEY link, KEYSET ks, KEY type, int *dxp) ;
int baseCallClipContig2Max (KEY link, int max, int *dxp) ;
void baseCallMakeSubclones (KEYSET ks) ;
BOOL baseCallRedoBaseCall (KEY key, int *np) ;int baseCallTileContigs (KEYSET ks, KEY type, int *dxp) ;
void statisticsCountGroup (KEYSET recu) ;

void nnContigTrain (KEY contig) ;
void doAssembleAllTraces (KEY target, KEYSET ks, char fonc) ;

void cDNAAlignInit (void) ;
int trackBadQuality (KEYSET ks) ;
void htilePrecompute (KEYSET ks0);

int geneZoneExport (KEYSET ks) ;

#endif
 

