/*  File: acemblyhook.c
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
 * Last edited: Dec  9 16:10 1998 (fw)
 * Created: Fri Oct 20 20:21:47 1995 (mieg)
 *-------------------------------------------------------------------
 */

/* @(#)acemblyhook.c	1.3 4/22/96 */

#include "acembly.h"
#include "classes.h"
#include "display.h"
#include "pick.h"
#include "lex.h"
#include "bs.h"			/* for isTimeStamps */
#include "query.h"

KEY   
_Clips,
_Later_part_of,
_OriginalBaseCall,
_Align_to_SCF,
_Quality,
_SCF_Position ;

static void parseArrayInit (void)
{
  extern BOOL baseCallDump (FILE*, Stack s , KEY k) ;
  extern BOOL baseCallParse (int level, KEY key) ;
  extern ParseFunc parseFunc[] ;
  extern DumpFunc dumpFunc[] ;

  parseFunc[_VBaseCall]  = baseCallParse ;
  dumpFunc[_VBaseCall]   = baseCallDump ;

                 /*   name            type       display     protected  case sensitive     */
		 /*    |                | visible   |  private   |       |    update names */
	         /*    |                |    |      |     |      |       |      |          */
  pickSetClassOption ("BasePosition",  'A', 'h',   ZERO, FALSE, FALSE, FALSE, FALSE ) ;
  pickSetClassOption ("BaseQuality",   'A', 'h',   ZERO, FALSE, FALSE, FALSE, FALSE ) ;
  pickSetClassOption ("BaseCall",      'A', 'h',   ZERO, FALSE, FALSE, FALSE, FALSE ) ;
}

void dnaAlignInit (void) ;


void acemblyInit (void)
{   
  isTimeStamps = FALSE ;
  parseArrayInit () ;
  dnaAlignInit () ;

  lexaddkey ("Clips", &_Clips, 0) ; 
  lexaddkey ("Later_part_of", &_Later_part_of, 0) ;
  lexaddkey ("OriginalBaseCall", &_OriginalBaseCall, 0) ;
  lexaddkey ("Align_to_SCF", &_Align_to_SCF, 0) ;
  lexaddkey ("Quality", &_Quality, 0) ;
  lexaddkey ("SCF_Position", &_SCF_Position, 0) ;

  queryRegisterWebQuery (owqQuery) ;
}



