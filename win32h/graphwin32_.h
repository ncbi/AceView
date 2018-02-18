/* graphwin32_.h : header file
 *  Author: Richard Bruskiewich
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 *   version of graph.h for inclusion within graphics package in
 *   files that use Win32 GDI calls.  #include "graph_.h before this file
 * HISTORY:
 * Last edited: Jun 11 03:15 1996 (rbrusk):
 *		-	code reorganized a bit and a new architecture underway
 * Created: Jul 21 18:40 1995 (rbrusk)
 *-------------------------------------------------------------------
 */

/* $Id: graphwin32_.h,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $ */

struct GRAPH_DOC_DATA
{
  int IDR ;
  char * FNExt ;
} ;

extern GRAPH_DOC_DATA AceGraphDesc[] ;


				


 
 
