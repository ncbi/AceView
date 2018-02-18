/*  file: win32print.c
 *  Author: Richard Bruskiewich, adapted from code written by
 *  Danielle et Jean Thierry-Mieg (mieg@kaa.cnrs-mop.fr)
 *  Copyright (C) R. Bruskiewich, J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:   WIN32 adaptation of graphprint.c
 * Exported functions:  graphPrint(), graphBoundsPrint()
 * HISTORY:
 * Last edited: Aug 27 18:00 1995 (rbrusk)
 * Created: Tue Mar  1 11:56:08 1994 (mieg)
 *-------------------------------------------------------------------
 */

// $Id: win32print.cpp,v 1.1.1.1 2002/07/19 20:23:26 sienkiew Exp $

#include "stdafx.h"
#include "win32.h"
#include "WinAce.h"
#include "CGraph.h" // gDev instantiation

//extern "C"
//{
//#include "acedb.h"
//}

extern "C" void graphPrint (void) 
{
  Graph targetGraph = graphActive () ;

  if (!graphExists (targetGraph) )
    return ;

  VIEWPTR->SendMessage(WM_COMMAND,ID_FILE_PRINT,NULL) ;
}

extern "C" void graphBoundsPrint (float uw, float uh, GraphFunc drawFunc )
{
  float olduw = gActive->uw ;
  int oldw = gActive->w ;
  float olduh = gActive->uh ;
  int oldh = gActive->h ;

  gActive->uw = uw ;
  gActive->w = (int)(gActive->uw * gActive->xFac) ;
  gActive->uh = uh ;
  gActive->h = (int)(gActive->uh * gActive->yFac) ;

  if (drawFunc)
    (*drawFunc) () ;

  graphPrint () ;

  gActive->uw = olduw ;
  gActive->w = oldw ;
  gActive->uh = olduh ;
  gActive->h = oldh ;
}

/*********************** end of file ***************************/
 
