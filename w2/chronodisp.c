/*  File: chrono.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Graphic interface to chrone
 *              This package is self contained since any function         
 *              called from here cannot be timed.
 *              If you compile without flag CHRONO all calls to the     
 *              chrono package just disappear out of the code.          
 *              See wh/chrono.h w1/chronoexe.c for how this works.
 *                                                             *
 *  1 routines are public :                                    *
 *     chronoShow                                              *
 *                                                             *
 *     If you compile without flag CHRONO all calls to the     *
 *     chrono package just disappear out of the code.          *

 * Exported functions: see chrono.h
 * HISTORY:
 * Created: Mon Jun 15 14:42:32 1990 (mieg)
 *-------------------------------------------------------------------
 */
/* $Id: chronodisp.c,v 1.1.1.1 2002/07/19 20:22:55 sienkiew Exp $ */

#include "acedb.h"

#include "array.h"
#include "chrono_.h"
#include "chrono.h"
#include "graph_.h" 
#include "freeout.h"

static int running = FALSE ;
static Graph chronoGraph = 0;
void chronoShow(void) ;

/******************************************/
/********* Graphic Routines ***************/
/******************************************/

static void chronoDoStart (void)
{
  if (!running && chronoStart ())
    {
      running = TRUE ;
      chronoShow () ;
    }
}

static void chronoDoStop (void)
{
  if (running == 1)
    { chronoStop () ; chronoShow () ; running = 0 ;}
}


static MENUOPT chronoMenu[]={
        {graphDestroy,"Quit"},
        {help,"Help"},
	{graphPrint,"Print"},
        {chronoDoStart,"Start"},
        {chronoDoStop,"Stop"},
        {chronoShow,"Show"},
        {0,0} };

/******************************************/

static void localDestroy(void)
{
  chronoStop() ;
  chronoGraph = 0 ;
}

/******************************************/

void chronoShow(void)
{
  int ll = 4, box = 0 ;
  
  if(! graphActivate(chronoGraph))
    {
      if (getGraphAcedbDisplayCreate() != NULL)
	chronoGraph =  (getGraphAcedbDisplayCreate())(getGraphAcedbChronoName()) ;
      else
	chronoGraph = graphCreate (TEXT_SCROLL, "", 0.2, 0.1, 0.5, 0.7) ;
      
      graphRegister(DESTROY, localDestroy) ;
    }
  else
    graphPop() ;
  
  graphClear();
  graphTextFormat(FIXED_WIDTH) ;
  
  graphColor (BLACK) ;
  
  switch (running)
    {
    case 0:
      graphText ("The chronometer is not yet running,",4,ll) ;
      graphText (" please press the Start button",4,ll + 1.5) ;
      ll++ ;
      break ;
    case 1:
      {
	Stack s = stackCreate (50) ;
	int level = freeOutSetStack (s) ;
	char cc, *cp, *cq ;

	chronoReport () ;
	freeOutClose (level) ;
	
	cp = stackText (s, 0) ;
	while (*cp == '/') cp++ ;
	while (*cp)
	  {
	    cq = cp ;
	    while (*cq && *cq != '\n') cq++ ;
	    cc = *cq ; *cq = 0 ;
	    box = graphBoxStart () ;
	    graphText (cp, 8, ll++) ;
	    graphBoxEnd () ;
	    if (!(ll%5)) graphBoxDraw (box, BLACK, PALEBLUE) ;
	    if (cc)
	      cp = cq + 1 ;
	    else
	      cp = cq ;
	  }
      }
      break ;
    }

  if (ll == 4)
    {
      graphText ("Nothing to report, ", 8, ll++) ;
      graphText ("no chrono enabled subroutine has yet been used", 8, ll + .5) ;
    }
  ll += 2 ;

  box = graphButtons(chronoMenu, 2, 2, 60) ;
  switch (running)
    {
    case 0:
      graphBoxDraw (box + 4, BLACK, LIGHTBLUE) ;
      break ;
    case 1:
      graphBoxDraw (box + 3, BLACK, LIGHTBLUE) ;
      break ;
    }
  graphTextBounds (100, ll + 3) ;
  graphRedraw() ;
  graphMenu (chronoMenu);
}
 

/*****************************************/
/*****************************************/
 
/*****************************************/
