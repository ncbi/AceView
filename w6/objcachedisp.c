/*  File: objcachedisp.c
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: graphical module that shows cache details
 * Exported functions:
 *              void cacheShow (void)
 * HISTORY:
 * Last edited: Nov 18 18:27 1998 (fw)
 * Created: Wed Nov 18 18:14:04 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: objcachedisp.c,v 1.1.1.1 2002/07/19 20:22:58 sienkiew Exp $ */

/************************************************************/
    /* For debugging, prints out the content of the cache
       Public;  Calls nothing; needs some work
    */

#include "acedb.h"
#include "graph.h"
#include "chrono.h"

#include "cache_.h"
#include "cache.h"
#include "bstree.h"
#include "bs.h"

/************************************************************/

void cacheShow(void)
{ CACHE v ;
  int nc, nt, ntt = 0, line = 1, ballocated, blocked, bknown, bmodif;
  int bsused, bsalloc, bsmem, ncr, ncw, nccr, nccw ;
  chrono("cacheShow") ;
  graphClear () ;

  graphText 
    ("Cache 2", 1, line++) ;
  graphText 
    ("contains the most recently used B objects",
     4, line++) ;
  line ++ ;
  cacheStatus (&ballocated, &blocked, &bknown, &bmodif, &ncr, &ncw, &nccr, &nccw) ;
  graphText(messprintf(
   "Size: %d allocated cache entries", ballocated), 1, line++) ;
 
  graphText(messprintf(
   "%d known, %d locked, %d modified Objects",
               bknown, blocked, bmodif),
            4, line++);
  graphText(messprintf(
   "Cache 2 access: %d read, %d save;  needing %d cache1-read, %d cache1-write",
               ncr, ncw, nccr, nccw),
            4, line++);
   BSstatus (&bsused, &bsalloc,&bsmem) ;
   graphText(messprintf(
    "referencing %d/%d BS, cells used/allocated",
    bsused,bsalloc), 4, line++) ;

  line ++ ; nt = 0 ;
  graphText("Update mode (still inside bsUpdate/bsSave)",1,line) ; line++ ;
  for(v = cacheGetLockedTop() ; v ; v = v->next)
    { nc = bsTreeSize (v->x) ; nt += nc ;
      graphText (messprintf("%3d:(%5d cells) %s:%s",
			    v->refCount, nc, 
			    className(v->key), name(v->key)), 3, line++) ;
    }
  graphText(messprintf("Total: %5d cells", nt), 1, line++) ; ntt += nt ;

  line ++ ; nt = 0 ;
  graphText("Locked mode (still inside bsCreate/bsDestroy)",1,line) ; line++ ;
  for (v = cacheGetKnownTop() ; v ; v = v->next)
    if (v->refCount)
    { nc = bsTreeSize (v->x) ; nt += nc ;
      graphText (messprintf("%3d:(%5d cells) %s:%s",
			    v->refCount, nc, 
			    className(v->key), name(v->key)), 3, line++) ;
    }
  graphText(messprintf("Total: %5d cells", nt), 1, line++) ; ntt += nt ;

  line ++ ; 
  graphText("Modified objects",1,line) ; line++ ; nt = 0 ;
  for (v = cacheGetKnownTop() ; v ; v = v->next)
    if (!v->refCount && v->hasBeenModified)
    { nc = bsTreeSize (v->x) ; nt+= nc ;
      graphText (messprintf("%3d:(%5d cells) %s:%s",
			    v->refCount, nc, 
			    className(v->key), name(v->key)), 3, line++) ;
    }
  graphText(messprintf("Total: %5d cells", nt), 1, line++) ; ntt += nt ;

  line ++ ; graphText("Known objects",1,line) ; line++ ; nt = 0 ;
  for(v = cacheGetKnownTop() ; v ; v = v->next)
    if (!v->refCount && !v->hasBeenModified)
    { nc = bsTreeSize (v->x) ; nt += nc ;
      graphText (messprintf("%3d:(%5d cells) %s:%s",
			    v->refCount, nc, 
			    className(v->key), name(v->key)), 3, line++) ;
    }
  graphText(messprintf("Total: %5d cells", nt), 1, line++) ; ntt += nt ;
  graphText(messprintf("Cumul: %5d cells", ntt), 1, line++) ; 
  /* mhmp line 8 au lieu de 7 07.10.97 (27.05.97)*/
  graphText (messprintf ("%5d cells referenced, %d lines", ntt, line), 15,8) ; 
  graphTextBounds (100,line + 5) ;
  graphRedraw() ;
  chronoReturn() ;
}

/************************* eof *****************************************/
