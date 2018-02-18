/* $Id: palette.c,v 1.10 2015/08/18 23:23:57 mieg Exp $ */
/*  Last edited: Nov 19 15:53 1998 (fw) */

#include "acedb.h"
#include "lex.h"
#include "pick.h"
#include "whooks/systags.h"
#include "whooks/sysclass.h"
#include "../w2/graphcolour.h"
#include "pref.h"
#include "session.h"
#include "bs.h"
#include "query.h"
#include "colours.h"	

#include <ctype.h>

#define BOOLEAN 1
#define FLOAT 2
#define INT 3
#define COLOUR 4

#include "graph.h"




typedef struct { Graph graph ; int box, color ; } PP ;
static PP *globalPp = 0 ;
static void paletteDraw (PP *pp) ;
static AC_HANDLE h = 0 ;

static MENUOPT  paletteMenu[] = {
  {graphDestroy, "Quit"},
  {graphPrint, "Print"},  
  {0, 0} };


static void paletteChooser (void *arg)
{
  /*
    int col = assInt (arg) ;
  */
  if (globalPp)
    {
      globalPp->color =  assInt (arg) ;
      paletteDraw  (globalPp) ;
    }
}

static void paletteDraw (PP *pp)
{
  int i, ii, color, height, width ;
  float line ;
  Array buttons = arrayCreate (33, COLOUROPT) ;
  COLOUROPT *cc, *oldcc = 0 ;
  char buf[64][13] ;
  KEY myCol[64] =
  {
 _PALEGRAY 
, _LIGHTGRAY
, _GRAY
,  _DARKGRAY

, _PALEBLUE 
, _LIGHTBLUE
, _LAVANDER
, _DARKBLUE

, _PALECYAN
, _LIGHTCYAN
, _CYAN 
, _DARKCYAN

, _PALEGREEN
, 0
, _LIGHTGREEN 
, _DARKGREEN 

, _PALEYELLOW 
, _WHITE 
, _YELLOW 
, _BLACK

, _PALEORANGE 
, _LIGHTORANGE
, _ORANGE 
, _BROWN 

, _PALERED
, _LIGHTRED
, _RED
, _DARKRED

, _PALEMAGENTA 
, _LIGHTMAGENTA
, _MAGENTA
, _PURPLE

, _PALEVIOLET 
, _LIGHTVIOLET 
, _VIOLET 
, _DARKVIOLET

, _BLUE1, _BLUE2, _BLUE3, _BLUE4
, _BLUE5, _BLUE6, _BLUE7, _BLUE8

, _GREEN1, _GREEN2, _GREEN3, _GREEN4
, _GREEN5, _GREEN6, _GREEN7, _GREEN8

, _RED1, _RED2, _RED3, _RED4
, _RED5, _RED6, _RED7, _RED8

, _MIDBLUE
, _BLUE 
, _GREEN
, _CERISE

  } ;

  graphClear () ;
  if (0)
    {
      int box = graphBoxStart () ;
      graphRectangle (10, 10, 300, 300) ;
      graphBoxEnd () ;
      graphBoxDraw (box, BLACK, RED) ;
      goto done ;
    }

  for (ii = 63 ; ii >= 0 ; ii--)
    {
      color = myCol[ii] ;

      cc = arrayp (buttons, ii, COLOUROPT) ;
      cc->f = paletteChooser ;
      cc->arg = assVoid (color) ;
      cc->fg = BLACK ;
      cc->bg = color - _WHITE + WHITE ;
      sprintf (buf[ii], "%s", color ? name(color) : "_") ;
      for (i = strlen(buf[ii]) ; i < 12 ; i++)
	buf[ii][i] = '_' ;
      buf[ii][i] = 0 ;
      cc->text = buf[ii] ;
      cc->next = oldcc ;
    }
      
  line = 1 ;
  graphText ("Selected colour", 1, 1) ;
  pp->box = graphBoxStart () ;
  graphRectangle (18, 1, 24, 3) ;
  graphBoxEnd() ;
  graphBoxDraw (pp->box, BLACK, pp->color  - _WHITE + WHITE) ;

  graphColor (pp->color  - _WHITE + WHITE) ;
  graphTextFormat (BOLD) ;
  if (pp->color) graphText (name(pp->color), 1, 2) ;
  graphTextFormat (PLAIN_FORMAT) ;
  graphColor (BLACK) ;

  line = 5 ;
  {
    int jj, red,green,blue,grey ;
    
    jj = 3 * (pp->color  - _WHITE) ;
    
    if (jj > 3 * NUM_TRUECOLORS) jj = 0 ;
    graphGetColourAsInt (pp->color  - _WHITE + WHITE, &red, &green, &blue, &grey) ;
    
    graphText (messprintf ("  Red: %3d", red), 27, 1) ;
    graphText (messprintf ("Green: %3d", green), 27, 2) ;
    graphText (messprintf (" Blue: %3d", blue), 27, 3) ;
    graphText (messprintf (" => Grey: %3d", grey), 42, 2) ;
  }
  
  line += 2 ;
  graphText ("You can change the shades, but not the names, so be reasonable", 
			2,line) ;
  line += 2 ;
  graphFitBounds (&width, &height) ;
  graphColouredButtons (arrp(buttons, 0, COLOUROPT), 4, &line, width) ;
  
 done:
  graphRedraw () ;
}

void paletteDisplay (void)
{
  if (globalPp && graphActivate (globalPp->graph))
    {
      graphPop () ;
      graphClear();
    }
  else
    {
      Graph paletteGraph = graphCreate (TEXT_SCROLL, "Palette",
				  .2, .2, 0.6, 0.52) ;
      graphMenu (paletteMenu);
      globalPp = (PP *) messalloc(sizeof(PP)) ;
      globalPp->graph = paletteGraph ;
      globalPp->color = _PALECYAN ;
      graphRegister (RESIZE, paletteDisplay) ;
      h = handleCreate () ;
    }
  
  paletteDraw (globalPp) ;
}

