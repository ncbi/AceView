/*  File: bqldisplay.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 2017
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: @(#)bqldisplay.c
 * Description:
 **  Constructs BQL queries and execute them in a graphic interface
 * Exported functions:
 **  bqlDisplayCreate
 * HISTORY:
 * Created : Jan 16, 2017 (mieg)
 *-------------------------------------------------------------------
 */

#define MALLOC_CHECK
#define ARRAY_CHECK


#include "acedb.h"
#include "querydisp.h"
#include "graph.h"

#define BBQLMAGIC 63866231

typedef struct BBQLSTUFF {
  AC_HANDLE h ;
  int magic ;
  Graph graph ;
  KEYSET ksActive ;
} BBQL ;

static Graph bbqlGraph = 0 ;

void  bqlDisplayCreate (void)
{ 
  BBQL *bbql ;
  KEYSET oldSet = 0 ;
  void   *look ;

  if(graphActivate(bbqlGraph))
    {
      graphPop() ;
      return ;
    }
  
  bbqlGraph = displayCreate(DtBqlDisplay) ;
  if (!bbqlGraph)
    return ;
  bbql = (BBQL *) messalloc(sizeof(struct BBQLSTUFF));
  bbql->h = ac_new_handle () ;
  bbql->magic = BBQLMAGIC;
  bbql->graph = bbqlGraph ; /* provision for multi windows */
  
  
  keySetActive(&oldSet,&look);
  if (oldSet)       /* if keyset exists, initialiaze the taglist */
    bbql->ksActive = oldSet ;
  
  return ;
} /* bqlDisplayCreate */

  
