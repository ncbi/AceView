/*  File: gifsubs.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: gifsubs.c,v 1.3 2015/08/11 22:14:46 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 09:32 1998 (edgrif)
 * * Dec 22 09:32 1998 (edgrif): Added dummy graphBusyCursor routine.
 * * Oct 22 10:56 1998 (edgrif): Changed graphOut etc. to graphWinOut etc.
 *              to match new device independent layer in graph routines,
 *              hopefully a temporary fix until graphis properly fitted up.
 * Created: Thu Feb 29 13:29:57 1996 (rd)
 *-------------------------------------------------------------------
 */

#include "regular.h"

#include "graph_.h"

/*********** globals required elsewhere **********/

float graphEventX, graphEventY ;
int menuBox ;
BOOL gWriteColors = TRUE ;
float  gPixAspect = 1.0 ;

/*********** routines doing nothing **************/

void graphDevInit (int *argcptr, char **argv) {}
void graphDevFinish (void) {}
void graphProcessEvents (void) {}
void graphGoto (float x, float y)  {}
void graphEvent (int action, float x, float y) {}
void graphPop (void) {}
void graphMessage (char *text) {}
void graphUnMessage (void) {}
BOOL graphInterruptCalled (void) { return FALSE ; }
void graphDevActivate (int show)  {}
void graphClipDraw (int x1, int y1, int x2, int y2) {}
void graphWhiteOut (void) {} 
void graphXorLine (float x1, float y1, float x2, float y2) {}
void graphXorBox(int k, float x, float y) {}
void graphBusyCursor (BOOL on) {}
void graphBoxMenu (int box, MENUOPT *menu) {}
void graphBoxFreeMenu (int box, FreeMenuFunction proc, FREEOPT *opts) {}
void graphNewBoxMenu (int box, MENU menu) {}
void graphRedraw (void) {}
void devMenuDestroy (void) {}
BOOL graphWindowSize (float *wx, float *wy, float *ww, float *wh) { return FALSE ; }

/*************************************************/
static int loopValue = 0 ;

BOOL graphLoopReturn (int retval) 
{ loopValue = retval ;
  return TRUE ;
}

int graphLoop (BOOL isBlock)
{ return loopValue ; }

/*************************************************/

/* screenx, screeny are for positioning and sizing created graphs
   SUN and other common screens are 1152x900 (1.3x900 = 1170)
   could/should get these from DISPLAY server in graphDevInit()
*/
static float screenx = 900 ;
static float screeny = 900 ;

static void fakeResize (Graph_ graph)
{ 
  switch (graph->type)
    {
    case PLAIN : case MAP_SCROLL:
      graphFacMake () ;
      break ;
    case TEXT_FIT: case PIXEL_FIT:
      graph->uw = graph->w / graph->xFac ;
      graph->uh = graph->h / graph->yFac ;
      break ;
    case TEXT_SCROLL:
    case PIXEL_SCROLL:
    case TEXT_FULL_EDIT:
    case TEXT_FULL_SCROLL:
      break ;
    }
}

static char* defaultName (Graph_ graph)
{ 
  static int gnum = 0 ;
  return strnew (messprintf ("graph-%d", ++gnum), graph->handle) ;
}

void graphDevCreate (Graph_ graph, float x, float y, float w, float h)
{
  graph->dev = graph->subdev = defaultName (graph) ;
  graph->w = w*screenx ;
  graph->h = h*screeny ;
  fakeResize (graph) ;		/* set up size state in graph structure */
}  

void graphSubDevCreate (Graph_ graph, float x, float y, float w, float h)
{
  graph = gActive ;
  graph->subdev = defaultName (graph) ;
  graph->w = w ;
  graph->h = h ;
  fakeResize (graph) ;		/* set up size state in graph structure */
}  

void graphDevDestroy (void) { gActive->dev = gDev = 0 ;	}
void graphSubDevDestroy (void) { gActive->subdev = gSubDev = 0 ; }

void graphRetitle (const char *title)
{ /*if (gSubDev)
    messfree (gSubDev) ; DOn;t free part of gActivate->handle ?? */
  gSubDev = strnew (title, gActive->handle) ;
}

/*************************************************/

void graphMapBounds (float ux, float uy, float uw, float uh, float aspect)
{
  if (gActive->type != MAP_SCROLL)
    messcrash ("MapBounds called on invalid graph type %d",
	       gActive->type) ;

  gActive->ux = ux ;
  gActive->uy = uy ;
  gActive->uw = uw ;
  gActive->uh = uh ;
  gActive->aspect = aspect ;

  gActive->w = gActive->h * gActive->uw / 
    (gActive->uh * gActive->aspect * gPixAspect) ;
  graphFacMake () ;
}

void graphScreenSize (float *sx, float *sy, float *fx, float *fy, int *px, int *py)
{
	/* no screen in this case: arbitrary limit 10000 each coord */
  int w, h ;
  gFontInfo (0, &w, &h) ;

  if (sx) *sx = (float)10000/screenx ;
  if (sy) *sy = (float)10000/screeny ;

  if (fx) *fx = (float)10000/w ;
  if (fy) *fy = (float)10000/h ;

  if (px) *px = 10000 ;
  if (py) *py = 10000 ;
}

void graphTextBounds (int nx, int ny)
{
  if (gActive->type != TEXT_SCROLL && gActive->type != TEXT_FULL_SCROLL)
    messcrash ("textBounds called on invalid graph type %d",
	       gActive->type) ;

  gActive->uw = nx ;
  gActive->uh = ny ;

  gActive->h = ny * gActive->yFac ;
  if (gActive->type == TEXT_FULL_SCROLL)
    gActive->w = nx * gActive->xFac ;
}

/* fake modif, used in forestdisp.c to prepare a longer printout
 * and restore old at the end
 */
float graphFakeBounds (float ny)
{
  float old = gActive->uh ;
  gActive->uh = ny ;
  gActive->h = ny * gActive->yFac ;

  return old ;
}

void graphPixelBounds (int nx, int ny)
{
  if (gActive->type != PIXEL_SCROLL)
    messcrash ("pixelBounds called on invalid graph type %d",
	       gActive->type) ;

  gActive->uw = nx ;
  gActive->uh = ny ;
  gActive->h = ny ;
  gActive->w = nx ;
}

/*************************************************/

void graphWhere (float *x1, float *y1, float *x2, float *y2)
{
  *x1 = gActive->ux ;
  *x2 = gActive->ux + gActive->uw ;
  *y1 = gActive->uy ;
  *y2 = gActive->uy + gActive->uh ;
}

BOOL gFontInfo (int height, int* w, int* h)
	/* hack - OK I think */
{
  if (w) *w = 8 ;
  if (h) *h = 13 ;
  return TRUE ;
}

/*************************************************/

BOOL graphRawMaps (unsigned char *forward, int *reverse) 
{
  int i ;
  if (reverse)
    for (i = 0 ; i < 256 ; ++i)
      reverse[i] = i ;
  if (forward)
    for (i = 0 ; i < 256 ; ++i)
      forward[i] = i ;
  return TRUE ;
}

BOOL graphWritableColors (void) { return TRUE ; }

static struct { float red, green, blue ; } colors[256] ;

void graphColorMaps (float *red, float *green, float *blue)
{
  int i ;

  for (i = 0 ; i < 256 ; ++i)
    { red[i] = colors[i/2].red ;
      green[i] = colors[i/2].green ;
      blue[i] = colors[i/2].blue ;
    }
}

BOOL graphSetColorMaps (unsigned char *red,
			unsigned char *green,
			unsigned char *blue)
{
  int i ;
  for (i = 0 ; i < 128 ; ++i)
    { colors[i*2].red = red[i] ;
      colors[i*2].green = green[i] ;
      colors[i*2].blue = blue[i] ;
    }
  return TRUE ;
}

/*************************************************/

void graphBoxDraw (int k, int fcol, int bcol)
{
  Box box = gBoxGet (k) ;

  if (fcol >= 0)
    box->fcol = fcol ;
  if (bcol >= 0)
    box->bcol = bcol ;
}

/*************************************************/

void graphWinOut(const char *text, char *ignored) { printf ("// %s\n", text) ; }
BOOL graphWinQuery(const char *text) { return freequery (text) ; }
BOOL graphWinPrompt(const char *prompt, char *dfault, char *fmt) { return freeprompt(prompt, dfault, fmt) ; }

/********* cut and paste stuff ************/

static char *cutbuf = 0 ;

char *graphPasteBuffer (void) { return cutbuf ? cutbuf : "" ; }

void graphPostBuffer (char *text)
{ if (cutbuf) messfree (cutbuf) ;
  cutbuf = strnew (text, 0) ;
}
/* dummy files */
BOOL graphWebBrowser (const char *link){
  return FALSE;
}

/********** end of file **********/
 
 
 
 
