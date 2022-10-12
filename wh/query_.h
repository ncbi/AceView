/*  File: query_.h 1.2
 *  Author: Gary Aochi LBL
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 11 10:05 1998 (fw)
 * 	-	changed fn() => fn(void) function prototypes to please compiler
 * * May 21 17:49 1994 (rd)
 * Created: Wed Apr 28 11:39:17 1993 (mieg)
 *-------------------------------------------------------------------
 * gha 04/28/93  removed line "#define graphTextScroll..." (jtm had added it)
 * gha 04/29/93  removed #define CLASSLIST_MAX, classlist static to *.c files
 * gha 06/14/93  split querybuild.c from querydisp.c
 */

/* $Id: query_.h,v 1.4 2015/09/25 17:29:06 mieg Exp $ */


#ifndef DEF_QUERY__H
#define DEF_QUERY__H

#include "acedb.h"
#include "graph.h"

#define BUFFER_SIZE 81   /* size of char buffer (strings) +1 for NULL. */
#define QBUFF_MULT  5    /* factor by buffer_size, for query command buffer */


unsigned char queryIsA(OBJ obj, KEY key, Array a, BitSet bb) ;

/* in querybuild.c */
int qbuild_selected_class ;  /* last entered/used class in a query */

/* in querydisp.c */
extern char resbuffer[QBUFF_MULT*BUFFER_SIZE]; /* for forming query commands */
extern KEYSET query_undo_set;

extern Graph querGraph ;
extern Graph qbeGraph ;
extern Graph qbuildGraph ;
     
#define ARR2STRING(a, i)  arrayp(array(a, i, Array), 0, char)
#include "pick.h"       /* for pickClass2Word(t) and pickVocList[] */


/* following defined in querydisp.c */
void queryCreate (void);
void queryWindowShowQuery(char *buffer) ;
void queryCommandEdit(void) ;
void queryCommandUndo(void) ;

/* following defined in querybuild.c */
void qbuildCreate (void);
int textAlphaOrder (const void *a, const void *b) ;
void destroyMenu(FREEOPT *menup, int max);
void qbuildNewClass(int classe);

/* following defined in qbedisp.c */
void qbeCreate (void);

#endif /* DEF_QUERY__H */

