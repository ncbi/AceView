/*  File: pepdisp.h
 *  Author: Clive Brown (cgb@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 21 14:28 1998 (fw)
 * Created: Tue Nov  7 14:40:38 1995 (cgb)
 *-------------------------------------------------------------------
 */

/* #define ARRAY_CHECK */

#include "acedb.h"

#include "display.h"
#include "colcontrol.h"
#include "lex.h"
#include "classes.h"
#include "sysclass.h"
#include "systags.h"
#include "tags.h"
#include "bs.h"
#include "peptide.h"
#include "a.h"
#include "method.h" 
#include <assert.h>

#define BIG_FLOAT 10000000
#define SMALL 0.0000001
#define PEP_DEBUG
#define DEF_WRAP 7
#define DEF_RES ""
#define MAXMESSAGETEXT 200
#define MAXACTIVETEXT 12
#define activecol CYAN
#define friendcol PALECYAN
#define keyfriendcol PALECYAN


/* the LOOK structure for peptide display */

typedef struct {
  MAPCONTROL map ;		/* the corresponding map */
  Array pep ;/* actual peptide sequence - coded 0..20 */
  int messageBox,aboutBox;
  char messageText[MAXMESSAGETEXT];
  int activeregionBox,activeStart,activeEnd;
  char activeText[MAXACTIVETEXT];
  KEY titleKey;
  int startcol;
  int endcol;
  BOOL block;
} PEPLOOK ;

extern void colBox(Array *colMap);
extern BOOL pepSequenceCreate (COLINSTANCE, OBJ );
extern BOOL hydrophobCreate (COLINSTANCE, OBJ );
extern BOOL homolCreate (COLINSTANCE, OBJ );
extern BOOL homolNameCreate (COLINSTANCE, OBJ );
extern BOOL pepFeatureCreate (COLINSTANCE, OBJ );
extern BOOL activeZoneCreate (COLINSTANCE, OBJ );
extern void * homolConvert (MAPCONTROL, void *);
extern void * FeatureConvert (MAPCONTROL, void *);

/***************** configuration data ******************/
 
 
 
