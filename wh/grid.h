/*  File: grid.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: header file for griddisp functions
 * Exported functions: gridCluster(), GRIDMAP structure
 *              (DisplayFunc) gridDisplay
 *              gridCluster
 *              gridClusterKey
 *              the GRIDMAP struct (also used by pmapconvert.c/cmapdisp.c)
 * HISTORY:
 * Last edited: Dec 11 15:08 1998 (fw)
 * Created: Mon Jan  6 18:58:11 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: grid.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

#ifndef _GRID_H
#define _GRID_H

#include "acedb.h"

typedef struct {
  KEY ctg ;
  float x1 ;
  float x2 ;
  int clump ;	   /* pos->clump filled with map index */
} GRIDMAP ;

BOOL gridDisplay (KEY key, KEY from, BOOL isOldGraph) ;
/* the public DisplayFunc */

void gridCluster (Array pos, Array map, float range) ; 
BOOL gridClusterKey (KEY key, Array map, float range) ;
/* range is the maximum gap size tolerated in maintaining a cluster */

#endif /* _GRID_H */
