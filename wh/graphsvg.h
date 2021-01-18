/**********************************************************************
 * File: graph.h 
 * Authors: Richard Durbin plus 
 *		 Jean Thierry-Mieg
 *		 and Christopher Lee
 * Copyright (C) J Thierry-Mieg and R Durbin, 1991-97
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 * Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Public header file for the graph package, this package
 *              is independent of the acedb database code but requires
 *		the array and messubs packages.
 *
 * Exported functions:
 * HISTORY:
 * Last edited: May 12 13:47 1999 (edgrif)
 * * May 12 13:46 1999 (edgrif): Added graphSetBlockMode call.
 * * Apr 16 13:34 1999 (edgrif): Add prototype for graphBoxSetMenu
 * * Jan  8 11:31 1999 (edgrif): Add correct prototype for graphPS.
 * * Dec 16 15:10 1998 (edgrif): Removed waitCursor calls, now done internally.
 * * Nov 19 15:02 1998 (edgrif): : Fixed callback func. dec. for
 *              graphTextScrollEditor.
 * * Oct 22 14:15 1998 (edgrif): Added message output functions that are
 *              independent of device level (X Windows etc), e.g. graphOut.
 *              Removed isGraphics flag.
 * * Oct 13 14:22 1998 (edgrif): Removed some acedb specific functions that
 *              do not belong here.
 * * May 14 08:00 1997 (rbrusk):
 *              - introduced GRAPH_FUNC_DCL symbol; see mystdlib.h
 *              - added graphSysBeep(), extern int menuBox
 *              - added an extra void* parameter to graphInit() to help
 *                provide for other program initialization data (in WIN32)
 *              - help() & helpOn() moved from acedb.h to graph.h w/ GRAPH_FUNC_DCL
 * * Oct 21 16:54 1996 (il) 
 * Created: Jan  1992 (rd)
 **********************************************************************/
/*  $Id: graphsvg.h,v 1.1 2017/10/24 15:52:03 mieg Exp $ */

#ifndef DEF_GRAPHSVG_H
#define DEF_GRAPHSVG_H

#include "regular.h"
#include "graph.h"
#include "aceio.h" 

BOOL svgGraphExport (Graph gId, ACEOUT out, BOOL do_size) ;

#endif
