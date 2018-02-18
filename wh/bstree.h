/*  File: bstree.h
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Packing unpacking kernel routines
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 19 10:01 1998 (fw)
 * Created: Tue Mar 10 03:40:07 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: bstree.h,v 1.3 2007/03/21 21:03:14 mieg Exp $ */

#ifndef DEFINE_BSTREE_h
#define DEFINE_BSTREE_h
 
void bsTreePrune      (CACHE_HANDLE bs) ; /* prunes bs and everything beyond from tree */
CACHE_HANDLE   bsTreeCopy       (CACHE_HANDLE bs, AC_HANDLE handle) ;
void bsTreeDump       (CACHE_HANDLE bs) ;
CACHE_HANDLE   bsTreeGet        (KEY key, AC_HANDLE handle) ;
void bsTreeStore      (CACHE_HANDLE bs, AC_HANDLE handle) ;  /* implies destruction since fools the tree */
void bsTreeKill (KEY key) ;  /* removes key from first cache */
CACHE_HANDLE   newBsTreeGet        (KEY key, AC_HANDLE handle) ;
void newBsTreeStore      (CACHE_HANDLE bs, BOOL wait) ;  /* non destructive */
void newBsTreeFlush (void) ; /* flush the waiting objects */
int  bsTreeSize (void *vbs);  /* counts the number of in the branch */

/* the void is a can't include bs_.h to define BS here - utter crap involved */
#endif
