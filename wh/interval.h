/*  File: interval.h
 *  Author: Ulrich Sauvage (ulrich@kaa.crbm.cnrs-mop.fr) 
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 11 17:31 1998 (fw)
 * Created: Wed Oct 13 12:56:32 1993 (ulrich)
 *-------------------------------------------------------------------
 */

/* $Id: interval.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

#ifndef DEFINE_interval_h
#define DEFINE_interval_h

#ifndef NON_GRAPHIC
#include "graph.h"
#endif /* !NON_GRAPHIC */

#define F_YES 0x01000000       /* (1 << 24) donnees experimentales */
#define F_NO 0x02000000        /* (1 << 25) */
#define F_PERE 0x01000000
#define F_MERE 0x02000000
#define F_PEME 0x03000000      /* (3 << 24) */
#define F_WHO 0x00ffffff       /* ((1 << 24) - 1) */
#define F_NO_FLAG 0x03ffffff   /* ((1 << 26) - 1) nettoyage de tous les flags meme F_ZERO */
#define F_FLAG 0xff000000      /* (255 << 24) */
#define F_V_FLAG 0xf8000000    /* (248 << 24) donne les flags autre que F_ZERO */
#define F_ZV_FLAG 0xfc000000   /* (252 << 24) donne les flags y compris F_ZERO */
#define F_ZERO 0x04000000      /* (1 << 26) Filter */
#define F_ISOLE 0x08000000     /* (1 << 27) */
#define F_SINGLE 0x10000000    /* (1 << 28) */
#define F_NON_DIAG 0x20000000  /* (1 << 29) attention par donnee et non par ligne comme les autres */
#define F_SUP_OVER 0x40000000  /* (1 << 30) */

#define F_G_ZERO 8                /* flag sur look->whatDis */
#define F_E_ZFLAG 16
#define F_E_TFLAG 32
#define F_E_FLAG 48
#define F_E_YN 64
#define F_UNK_T 7

#define F_SHOWM 8                 /* flag sur look->mapStatus */
#define F_SOMAR 16
#define F_HIST 32
#define F_SORT 64
#define F_ASSSEG 128
#define F_BOUT 7
#define F_IS_FIXED 256
#define F_NO_CONTIG 512

#define M_ASSEMB 1                /* choix de la methode dans look->method */
#define M_FINGPR 2
#define M_DEFDUP 3
#define M_EDWARD 4
#define M_FIN 5

#define INTTREEMAG 940627105

typedef struct { KEYSET ks, x ;} TREE_DEF ;
typedef struct INTTREESTUFF { int magic ;
			      int prev, yTree ;
			      BOOL display ;
			      KEY map ;
			      KEYSET def, genes ;
			      int nd, nm ;
			      Array cOrder, mind, dinm, segment, tree, tabledis ;
			      Associator distances ;
			      TREE_DEF trd ;
			    } *INTTREE ;

typedef struct DEFCPTSTUFF { int magic, id ;
			     int nm, nd, dlimit, taux ; /* nb ligne, nb colone, distance limite */
			     int Line ; /* for displays */
			     int mapStatus, method, whatDis, nbContig ;
			     int nboligBox, disBox, choixBox, paramBox ;
			     int step, tour ;
			     char nboligo[16], distance[16], choix[16], param[30] ;
			     BOOL manEntry, display, assembling ;
			     Array linOrder, colOrder, maillon ;
			     Array defInMar, marInDef ;
			     Array knownPairs, knownOrder, dataArray, gmap ;
			     KEYSET def, mar, rejected, actif ;
			     KEY selectedMap, link ;
			     Associator assDnaGet ;
                             KEYSET taceActif ; /* do not destroy this one */
#ifndef NON_GRAPHIC
			     Graph defMapGraph, defTreeGraph, defMapCtlGraph ;
#else
                             int defMapGraph, defTreeGraph, defMapCtlGraph ;
#endif /* !NON_GRAPHIC */
			   } *DEFCPT ;

/************************************************************/

void intrinsicTreeDestroy(void) ;
void intCptSupFlag(Array donnee, KEY flag) ; /* Supprime les flags d'un array de KeySet */
int intCptTree(DEFCPT look) ; /* RETURNS THE NUMBER OF SEGMENTS */
void defCptAddSeqIn(KEY link, KEY key) ;
void plotSeqDestroy(void) ;
KEY defCptExecuteCommand (KEY link, KEYSET def, KEYSET actif, char *com, char *param) ;
BOOL defCptOrderKeySet (KEYSET aTrier, KEY tag, KEY link, KEY sub) ;
DEFCPT defCptGetLook (KEY link) ;
void defCptForget (KEY link, KEY key) ;
void defCptDestroyLook (KEY link) ;
void defCptChangeLook (DEFCPT look, KEY oldLink, KEY newLink) ;

/* de abifix.c */
int trackVector (KEYSET ks, KEY key, BOOL force) ;

#endif /* DEFINE_interval_h */
