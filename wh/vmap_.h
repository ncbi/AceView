/*  File: vmap_.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1999
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: private header for the VMAP (Vertical Map) package
 *              to be included only by members of that package (vmapXXX.c)
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  8 13:18 1999 (fw)
 * Created: Thu Jan  7 16:53:30 1999 (fw)
 *-------------------------------------------------------------------
 */


/* $Id: vmap_.h,v 1.2 2008/03/23 20:43:09 mieg Exp $ */

#ifndef _VMAP__H
#define _VMAP__H

#include "acedb.h"

#include "vmap.h"
#include "map.h"
#include "graph.h"

#if defined(applec)
#include "Math.h"
#endif

/************************************************************/

typedef struct vMapDataStruct {
  Associator loc2two ;	/* all loc2 Ass's map key -> KEYSET */
  Associator loc2multi ;
  Associator loc2inside ;
  Associator loc2outside ;
  KEYSET orderedLoci ;
  Associator key2seg ;
} *vMapData ;


typedef struct VerticalMapStruct {
  magic_t *magic;        /* == VerticalMap_MAGIC */
  KEY   key, titleKey ;
  int   activeBox ;
  int	  messageBox ;
  char  messageText[128] ;
  int   resolution ;
  float errorScale ;  /* to scale the error on gene locations */
  unsigned int flag ;
  Array segs;      	/* array of SEG's from the obj */
  Array boxIndex ;    /* box -> int (index into look->segs) */
  Stack remarkStack ;
  int lastTrueSeg ;	/* arrayMax(look->segs) after Convert */
  vMapData data ;      /* private data for vMapdata.c */
  MAP map ;
  Graph graph ;
} *VerticalMap;

extern magic_t GRAPH2VerticalMap_ASSOC;

#define FLAG_ALL_DATA		0x00000001
#define FLAG_HIDE_HEADER	0x00000002

typedef struct
  { KEY key, symbol ;
    float x, dx ;
    unsigned int flag ;
  } SEG ;
#define segFormat "kkffi"

  /* seg->flag flags  i use the leftmost octet to store colour*/
#define FLAG_CLONED		0x00000001
#define FLAG_MARKER		0x00000002
#define FLAG_RELATED		0x00000004
#define FLAG_STRESSED		0x00000008
#define FLAG_ANTI_RELATED	0x00000010
#define FLAG_PHYS_GENE		0x00000020
#define FLAG_ANTI_STRESSED      0x00000040
#define FLAG_HAVE_DATA		0x00000080
#define FLAG_WELL_ORDERED	0x00000100
#define FLAG_HIDE		0x00000200
#define FLAG_HIGHLIGHT		0x00000400
#define FLAG_CENTROMERE		0x00000800
#define FLAG_DARK_BAND		0x00001000
#define FLAG_NOR		0x00002000
#define FLAG_MOVED		0x00004000
#define FLAG_DEFICIENCY		0x00008000
#define FLAG_DUPLICATION	0x00010000
#define FLAG_BALANCER   	0x00020000
#define FLAG_CHROM_BAND    	0x00040000
#define FLAG_COLOUR      	0x00080000 /* colour is coded as << 28 */
#define FLAG_ANY_LOCUS   	0x00100000
#define FLAG_ANY_INTERVAL   	0x00200000
#define FLAG_MULTIPLE_LOCUS   	0x00400000
#define FLAG_PROBLEM    	0x00800000

#define FLAG_MARGINAL (FLAG_RELATED | FLAG_STRESSED)

extern KEY dataKey ;
extern KEYSET newBadData ;
extern Array loci ;     /* workspace for recombination data */
extern Array counts ;   /* workspace for recombination data */

#define LOG_10 	  2.302585	/* log(10) */

extern MENUOPT vMapDataMenu[] ;
extern MENUOPT vMapIntervalMenu[] ;

/************************************************************/


	/* in vmapdisp.c */
VerticalMap currentVerticalMap (char *caller);
int  vMapOrder (const void *a, const void *b) ; /* for arraySort() call */ 

	/* in vmapdata2.c */
void vMapDraw (VerticalMap look, KEY key) ;

void vMapShowData (void) ;
void vMapUpArrow (MAP m, float x, float y) ;
void vMapDownArrow (MAP m, float x, float y) ;
void vMapSelect (VerticalMap look, int box) ;

/* the MapDrawColFunc's */
void vMapDraw2pt (LOOK genericLook, float *offset) ;
void vMapDrawMultiPt (LOOK genericLook, float *offset) ;


      /* in vMapdrag.c */
void dragFindAllErrors(void) ; 
void vMapDragButton (void) ;
void vMapDragGene(int box) ;
void vMapDragUndo (void);


void vMapAddPhysGenes (Array segs) ;
void vMapToPMap (VerticalMap look, float y) ;
void pMapToVMap (KEY contig, KEY from, int x) ;

void vMapMakeDbn (int box) ;
void vMapShowDataFromMenu (int box) ;
void vMapReport2pt (VerticalMap look, SEG *seg) ;
void vMapReportMultiPt (VerticalMap look, SEG *seg) ;
void vMapDataDestroy (VerticalMap look) ;

/* the MapDrawColFunc's */
void vMapDrawContigs (LOOK genericLook, float *offset) ;
void vMapDrawPhysGenes (LOOK genericLook, float *offset) ;
void vMapReversedPhysical (LOOK genericLook, float *offset) ;
void vMapDrawDbn (LOOK genericLook, float *offset) ;


#endif /* _VMAP__H */

/********** end of file **********/
 
