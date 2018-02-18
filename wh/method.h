/*  File: method.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: methods header file for fmap & pepmap = sequence displays
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  4 14:48 1998 (fw)
 * * Jul 27 09:34 1998 (edgrif): Removed #define METHOD_HOMOL, it's not
 *      found anywhere else....
 * * Jun 24 13:40 1998 (edgrif): Removed reference to methodInitialise,
         this is now an internal routine.
 * Created: Sat Jul 25 20:28:37 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: method.h,v 1.2 2005/08/24 22:59:45 mieg Exp $ */

#include "regular.h"

typedef struct 
{ unsigned int flags ;
  KEY col, upColor ;		/* mieg */
  float min, max ;
  float width ;
  char symbol ;
  float priority ;
  float gaFac ;			/* for gene assembly */
  float histBase ;
} METHOD ;

extern Array methodInfo ;	/* of METHOD* indexed by KEYKEY(method) */

#define METHOD_BLASTN		0x00000001
#define METHOD_DONE		0x00000002
#define METHOD_FRAME_SENSITIVE	0x00000004
#define METHOD_STRAND_SENSITIVE	0x00000008
#define METHOD_SCORE_BY_OFFSET  0x00000010
#define METHOD_SCORE_BY_WIDTH   0x00000020
#define METHOD_BLIXEM_X		0x00000040
#define METHOD_BLIXEM_N		0x00000080

/* There is no such routine as fMapSetMethod so I'm not sure if this is      */
/* even needed now...                                                        */
#define METHOD_INTERNAL		0x00000100 /* from fMapSetMethod */

#define METHOD_BUMPABLE		0x00000200
#define METHOD_PERCENT          0x00000400
#define METHOD_CALCULATED	0x00000800 /* used in addOldSegs */
#define METHOD_EMBL_DUMP	0x00001000
#define METHOD_SCORE_BY_HIST	0x00002000

#define METHOD_SCORE		0x00002030 /* NB combination */

#define METHOD_BELVU    	0x00004000 /* esr, for PEPMAP */
#define METHOD_BLIXEM_P		0x00008000

#define METHOD_SHOW_UP_STRAND	0x00010000 /* mieg */

#define METHOD_FEATURE		0x20000000

typedef struct
{ KEY method ;
  int x1, x2 ;
  float score ;
} HOMOLINFO ;

/**************** functions *******************************/

int methodAdd (KEY view, KEY key) ;
void methodSet (char *name, int col, unsigned int flags, float priority,
		float width, char symbol, float min, float max) ;
METHOD *method (KEY view, KEY key) ;

/**********************************************************/
 
 
 
 
 
 
 
