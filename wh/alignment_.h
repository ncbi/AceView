/*  File: alignment_.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: private header for alignment package
 *              comprising of alignment.c (generic)
 *              and alignmentdisp.c (graphical display)
 * HISTORY:
 * Last edited: Nov 19 11:48 1998 (fw)
 * Created: Thu Nov 19 11:46:21 1998 (fw)
 *-------------------------------------------------------------------
 */
#ifndef _ALIGNMENT__H
#define _ALIGNMENT__H


#define GAP 21

typedef struct {		/* Alignment component structure */
  KEY key ;
  int start, end ;
  int flag ;
  char *seq ;
  int *pos ;
} ALIGN_COMP ;

typedef struct {		/* Alignment structure */
  KEY key ;
  Array comp ;			/* of ALIGN_COMP */
  AC_HANDLE handle ;
} ALIGNMENT ;

ALIGNMENT *alignGet (KEY key);

#endif /* _ALIGNMENT__H */
