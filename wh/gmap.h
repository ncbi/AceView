/*  File: gmap.h
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: header file for gmap operations
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  8 11:38 1999 (edgrif)
 * * Jan  8 11:37 1999 (edgrif): Add missing prototype for gMapDrawGIF.
 * * Dec 10 16:01 1998 (edgrif): Add prototype for gMapPhysClone2map.
 * Created: Thu Nov 26 05:15:10 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: gmap.h,v 1.4 2008/03/23 20:43:50 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "graphAcedbInterface.h"			    /* graph -> ace header. */
#include "colcontrol.h"
#include "lex.h"
#include "bs.h"
#include "a.h"
#include "systags.h"
#include "classes.h"
#include "sysclass.h"
#include "tags.h"
#include "session.h"
#include "query.h"
/*#include "bump.h"*/

/*#include "graphAcedbInterface.h"*/

#if defined(applec)
#include "Math.h"
#endif

typedef struct TwoPtData {
  KEY loc1;
  KEY loc2;
  KEY type;
  float distance;
  float error;
  float y1, y2;
  int n1, n2, n3, n4; /* must be consecutive */
  int count;
} *TWOPTDATA;

typedef struct MultiPtData {
  Array counts;
  Array loci;
  float max, min;
} *MULTIPTDATA;

typedef struct {
  KEY obj;
  KEY map;
} KEYPAIR;

 
#define MAXMESSAGETEXT 200

typedef struct GeneticMapStruct
{
  magic_t *magic;        /* == GeneticMap_MAGIC */
  MAPCONTROL map ;
  int	  messageBox ;
  char  messageText[MAXMESSAGETEXT] ;
  int   aboutBox;
  KEY   titleKey;
  /* Global data about our displayed keyset. These are set up by getPos */
  KEYSET orderedLoci; /* the keys of all the ordered loci */
  Associator posAss; /* associator for looking up locus positions */
  /* columns to display map data in */
  COLINSTANCE twoPtCurrentColumn;
  COLINSTANCE multiPtCurrentColumn;
  COLINSTANCE dbnCurrentColumn;
  
  KEYSET highlight; /* always valid*/
  KEYSET hidden;
  /* here we have the total state needed to do selection display */
  KEY selectKey;     /* key of object currently selected */
  KEYSET positive;  /* only valid when active box is */
  KEYSET negative;
  KEYSET neighbours;
  KEYSET twoPt;
  KEYSET multi;
  KEYSET negativeMarkers;
  float x1, x2; /* co-ords of seg the keysets relate to */
  BOOL friendsInfoValid; /* set when above keysets are valid */
  BOOL neighboursInfoValid;
  BOOL dataInfoValid;
  BOOL coOrdsInvalid;
  BOOL negativeMarkersValid;
  /* The segs which represent the Contains items are ubiquitous, 
     so we store them here. */
  Array segs ;
  /* The remark column is atypical */
  Array remarkSegs;
} *GeneticMap;

typedef struct
  { KEY key;
    float x, dx ;
    unsigned int flag ;
  } GMAPSEG ;
#define segFormat "kffi"


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
#define FLAG_P_TELOMERE		0x00000200
#define FLAG_Q_TELOMERE	        0x00000400
#define FLAG_CENTROMERE		0x00000800
#define FLAG_DARK_BAND		0x00001000
#define FLAG_NOR		0x00002000
#define FLAG_MOVED		0x00004000
#define FLAG_DEFICIENCY		0x00008000
#define FLAG_DUPLICATION	0x00010000
#define FLAG_BALANCER   	0x00020000
#define FLAG_CHROM_BAND   	0x00040000
#define FLAG_COLOUR      	0x00080000 /* colour is coded as << 28 */
#define FLAG_ANY_LOCUS   	0x00100000
#define FLAG_ANY_INTERVAL   	0x00200000
#define FLAG_MULTIPLE_LOCUS   	0x00400000
#define FLAG_PROBLEM    	0x00800000
#define FLAG_NEIGHBOUR		0x01000000


/* exports of gmapdisp.c */
GeneticMap currentGeneticMap (char *callerFuncName);
BOOL gMapDisplay (KEY key, KEY from, BOOL isOldGraph);
void gMapMakeAll (void);
void gMapHighlightKey(MAPCONTROL map, KEY key);
void gMapUnselect(MAPCONTROL map);
int gMapOverlap(MAPCONTROL map, KEY key, float z1, float z2);
BOOL gMapIsNeighbour(MAPCONTROL map, KEY key);
int gMapOrder (const void *a, const void *b);
void *gMapConvert (MAPCONTROL map, void *params);
void gMapSaveDetails (Array segs, KEY map);
void gMapAddToHeader(MAPCONTROL map, char *string);
BOOL gMapFollowMap (COLCONTROL control, MAPCONTROL oldMap, KEY newMapKey) ;

/* exports of gmapconvert.c */
KEY gMapCacheKey(KEY map, char *suffix);
void gMapCacheStore(KEY map, char *suffix, Array a, char *format);
void gMapCacheKill(KEY map);
BOOL gMapNeighbours(MAPCONTROL map, KEYSET *keyset, KEY key);
BOOL gMapPositive(MAPCONTROL map, KEYSET *ks, KEY key);
BOOL gMapNegative(MAPCONTROL map, KEYSET *ks, KEY key);
BOOL gMapMultiPt(MAPCONTROL map, KEYSET *ks, KEY key);
BOOL gMap2Pt(MAPCONTROL map, KEYSET *ks, KEY key);
BOOL gMapGetMultiPtData(MAPCONTROL map, 
			KEY multi, 
			AC_HANDLE handle, 
			MULTIPTDATA *ret);
BOOL gMapGet2PtData(MAPCONTROL map, 
		    KEY key, 
		    AC_HANDLE handle, 
		    TWOPTDATA *ret);
BOOL getPos (MAPCONTROL map, KEY key, float *y);
void setTestPos (MAPCONTROL map, KEY key, float pos);
int gMapGetMapObject
(KEY key, KEY map, KEY hint, int ind, float *xret, float *dxret, OBJ *objp, int *count);
Array  gMapMapObjects(KEY chrom);
KEY gMapKeyOnTag(KEY key, KEY tag);

void gMapDrawGIF (void) ;


/* exports of gmapmarkercol.c */
extern struct ProtoStruct gMapMainMarkersColumn;
extern struct ProtoStruct gMapMiniChromBandsColumn;

/* exports of gmapintervalcol.c */
extern struct ProtoStruct gMapChromBandsColumn;
extern struct ProtoStruct gMapJTMIntervalColumn;
extern struct ProtoStruct gMapRDIntervalColumn;

/* exports of gmaplocuscol.c */ 
extern struct ProtoStruct gMapPointColumn;

/* exports of gmapsubmapcol.c */ 
extern struct ProtoStruct gMapSubMapColumn;

/* exports of gmapdatacol.c */
void multiPtAddKey(COLINSTANCE instance, KEY key);
extern struct ProtoStruct gMapMultiPtColumn;
void twoPtAddKey(COLINSTANCE instance, KEY key);
extern struct ProtoStruct gMapTwoPtColumn;
void dbnAddKey (COLINSTANCE instance, KEY locus);
extern struct ProtoStruct gMapLikelihoodColumn;

/* exports of gmapdata.c */
BOOL logLikeMulti (MAPCONTROL map, MULTIPTDATA data, KEY dataKey, float *ll);
void multiBoundCheck (MAPCONTROL map, KEY locus, KEY multikey, 
		      float *min, float *max, KEY *minLoc, KEY *maxLoc);
float logLike2pt (TWOPTDATA data, float dist);
BOOL best2p (TWOPTDATA data, 
	     float *best, 
	     float *lo, 
	     float *hi);
float logLikeLocus (MAPCONTROL map, KEY locus);
BOOL boundFind (MAPCONTROL map, KEY locus, float *min, float *max,
		       KEY *minLoc, KEY *maxLoc);
void gMapGetBadData(void);
extern MENUOPT gjmMenu[] ;

/* exports of gmapphys.c */
extern struct ProtoStruct gMapPhysGenesColumn;
extern struct ProtoStruct gMapContigColumn;
extern struct ProtoStruct gMapRevPhysColumn;
void pMapToGMap (KEY contig, KEY from, int x) ;
Array gMapGetClones(KEY key);
BOOL  gMapPhysClone2map (KEY clone, KEY *seqp, KEY *mapp, float *xp) ;


/* exports of gmapremarkcol.c */
void remarkRegister(MAPCONTROL map, KEY key, float coord);
extern struct ProtoStruct gMapRemarkColumn;

/* exports og gmapposnegcol.c */
extern struct ProtoStruct gMapPosNegColumn;

#define LOG_10 	  2.302585	/* log(10) */


/********** end of file **********/
 
 
 
