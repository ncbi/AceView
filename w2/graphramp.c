/*  File: graphxlib.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 * this file written by F Wobus, with contributions from
 * R Durbin and E Sonnhammer.
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: greyRampTool for easy real-time manipulation
                on the colors used in graphPixel-images
 * Exported functions:  graphRampTool

 * HISTORY:
 * Last edited: Dec  4 11:50 1998 (fw)
 * * Jun  7 15:33 1996 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: graphramp.c,v 1.1.1.1 2002/07/19 20:22:55 sienkiew Exp $ */

#include "regular.h"

#include "graph_.h"

/**************************************************/
/**************************************************/
/*  ColorMap-Tool for the ACeDB-GraphPackage */
/*  written by Friedemann Wobus (fw@sanger.ac.uk) */
/*  includes Richard's function graphGreyRamp, which interacts with
    the graphRampTool function in both directions
    -23.03 1994 added resize feature 
    now the display will be resized when window size changes
    -24.03 1994 added interaction between graphGreyRamp and graphRampTool
    a call of graphGreyRamp from outside graphRampTool will alter
    the display of the rampTool window
    changed undo to be a undo/redo function
    -06.04 1994 made rampMin and rampMax global to the library
 */
#if (defined SGI || defined SUN)
  extern int atoi (const char*) ;
#endif
static void colorSetDisplay (void) ;

int rampMin=0, rampMax=255 ;	/* global within graphlibrary */
Graph rampGraph = 0 ;
int greyChanged ;


static int oldmin, oldmax ;
static BOOL SolidLines = TRUE ;
static UCHAR greyMap[256] ;

void graphGreyRamp (int min, int max)
{
  int i ;
  float fac ;
  Graph oldgraph;
  unsigned char ramp[128] ;

  min = min <= 0 ? 0 : min ;
  max = max >=255 ? 255 : max ;
  if (min == max)
    max = min+1 ;

  if (min < max)
    { for (i = 0 ; i < min ; ++i)
	ramp[i/2] = 0x0 ;
      fac = 0xff / (float)(max-min) ;
      for ( ; i < max && i < 256 ; ++i)
	ramp[i/2] =  fac * (i - min) ; 
      for ( ; i < 256 ; ++i)
	ramp[i/2] = 0xff ;
    }
  else
    { for (i = 0 ; i < max ; ++i)
	ramp[i/2] = 0xff ;
      fac = 0xff / (float)(min-max) ;
      for ( ; i < min && i < 256 ; ++i)
	ramp[i/2] = fac * (min - i) ; 
      for ( ; i < 256 ; ++i)
	ramp[i/2] = 0x0 ;
    }

  if (!graphSetColorMaps (ramp, ramp, ramp))
    return ;

  rampMin = min ; rampMax = max ;

  oldgraph = graphActive();
    
  if (graphExists (rampGraph))
    { if (min != oldmin && max != oldmax) colorSetDisplay () ; }

  graphActivate(oldgraph);

}

#define MAXY 20
#define MINY 110
#define BORDER 20
#define SIZEX(x) ((int)((x) * ((float)winx/CWINX)))
#define SIZEY(x) ((int)((x) * ((float)winy/CWINY)))
#define UNSIZEX(x) ((int)((x) * (CWINX/(float)winx)))
#define UNSIZEY(x) ((int)((x) * (CWINY/(float)winy)))
#define CORRX(x) (((x) % 4 == 0) ? (x) : (((int)((x)/4))+1)*4 )

static int defMin, defMax, oldx ;
static int sliderBox, maxTextBox, minTextBox ;
static int maxBox, minBox, thresBox ;
static char maxTx[4] ;
static char minTx[4] ;
static UCHAR *wedge = 0 ;
static float CWINX = 335 ;
static float CWINY = 2*BORDER+(MINY-MAXY) ;
static int winx, winy ;

static void colorDrawWindow (void) ;
static void colorPickSlider (int) ;
static void colorDragSlider (double,double) ;
static void colorUpSlider (double,double) ;
static void colorDragThres (double,double) ;
static void colorUndo (void) ;
static void colorXorLines (void) ;
static void colorSolidLines (void) ;
static void colorEndEntry (char *) ;
static void colorSwap (void) ;

static MENUOPT colorFuncMenu[] = {
  {graphDestroy,	"Quit"},
  {colorUndo,		"Undo"},
  {colorSwap,		"Swap"},
  {0,0} };

static void colorDestroy (void)
{
  /*messfree (wedge) ;*/ /* No - done by XDestroyImage in graphDestroy */
}

void graphRampTool (void)
{
  /* I don't know why this dotter fix sent by Christian Iseli (Bug SANgc02086) has */
  /* this bit of code in it which is never included....note the bug in that the    */
  /* the static notfirst is not initialised.                                       */
#if 0  /* CI */
  static int notfirst;

  if (!graphWritableColors())
    { 
      if (notfirst)
	messout ("Pixel values are readonly: colormap is being "
		 "shared with other GRAPH package applications "
		 "if possible") ;
      else
	notfirst = 1;
      return ;
    }
#endif /* 0 CI */


  if (graphActivate (rampGraph)) /* checks graphExists() first */
    { graphPop () ;
      return;
    }

  rampGraph = graphCreate (PIXEL_FIT, "Greyramp Tool", 0, 0,
			    CWINX/900.0, CWINY/900.0) ;
  graphRawMaps (greyMap, 0) ;	/* to fill local greyMap copy */
  colorDrawWindow () ;
  defMin =  rampMin ; defMax =  rampMax ;
}

static void colorDrawWindow (void)
{
  int i,j;
  float x1, y1 ;

  rampMin = rampMin <= 0 ? 0 : rampMin ;
  rampMax = rampMax >=255 ? 255 : rampMax ;

  graphActivate (rampGraph) ;
  graphClear () ;

  graphRegister (PICK,		colorPickSlider) ;
  graphRegister (RESIZE,	colorDrawWindow) ;
  graphRegister (DESTROY,	colorDestroy) ;
  graphMenu (colorFuncMenu); 

  graphWhere (0, 0, &x1, &y1) ; winx = x1 ; winy = y1 ;

  /* if (wedge) messfree (wedge) ;*/ /* No - done by XDestroyImage in graphDestroy */


  wedge = (UCHAR *) malloc /* Not messalloc since freed by XDestroyImage in graphDestroy */
    (CORRX(SIZEX(256)) * SIZEY((MINY-MAXY-20)) * sizeof(UCHAR)) ;

  graphBoxStart();
  graphPixelsRaw ((char*)wedge, SIZEX(256), SIZEY(MINY-MAXY-20),
		  CORRX(SIZEX(256)), SIZEX(BORDER), SIZEY(MAXY+10)) ;
  for (j = 0 ; j < SIZEY(MINY-MAXY-20) ; j++)
    for (i = 0 ; i < CORRX(SIZEX(256)) ; i++)
      wedge[j*CORRX(SIZEX(256))+i] = greyMap[UNSIZEX(i)] ;
  graphBoxEnd();

  graphColor (BLACK);

  graphRectangle (SIZEX(BORDER-9), SIZEY(MAXY-10),
		  SIZEX(BORDER+255+9), SIZEY(MINY+10));

  maxBox = graphBoxStart () ;
  graphLine (SIZEX(BORDER+rampMax-5), SIZEY(MAXY-5),
	     SIZEX(BORDER+rampMax), SIZEY(MAXY+5));
  graphLine (SIZEX(BORDER+rampMax), SIZEY(MAXY+5),
	     SIZEX(BORDER+rampMax+5), SIZEY(MAXY-5));
  graphLine (SIZEX(BORDER+rampMax-5), SIZEY(MAXY-5),
	     SIZEX(BORDER+rampMax+5), SIZEY(MAXY-5));
  graphBoxEnd ();

  minBox = graphBoxStart () ;
  graphLine (SIZEX(BORDER+rampMin-5), SIZEY(MINY+5),
	     SIZEX(BORDER+rampMin), SIZEY(MINY-5));
  graphLine (SIZEX(BORDER+rampMin), SIZEY(MINY-5),
	     SIZEX(BORDER+rampMin+5), SIZEY(MINY+5));
  graphLine (SIZEX(BORDER+rampMin-5), SIZEY(MINY+5),
	     SIZEX(BORDER+rampMin+5), SIZEY(MINY+5));
  graphBoxEnd ();
  
  thresBox = graphBoxStart () ;
  graphRectangle (SIZEX(BORDER+((rampMin+rampMax)/2)-5),
		  SIZEY(((MINY-MAXY)/2+MAXY)-5),
		  SIZEX(BORDER+((rampMin+rampMax)/2)+5),
		  SIZEY(((MINY-MAXY)/2+MAXY)+5) ) ;
  graphBoxEnd () ;
  graphBoxDraw (thresBox, RED, TRANSPARENT) ;

  graphRectangle (SIZEX(BORDER+255+15),    SIZEY(MAXY-7),
		  SIZEX(BORDER+255+15)+34, SIZEY(MAXY-7)+15);
  graphRectangle (SIZEX(BORDER+255+15),    SIZEY(MINY-7),
		  SIZEX(BORDER+255+15)+34, SIZEY(MINY-7)+15);

  graphButton ("Quit", graphDestroy, SIZEX(BORDER+255+15), SIZEY(MAXY+15)) ;
  graphButton ("Undo", colorUndo, SIZEX(BORDER+255+15), SIZEY(MAXY+15+22)) ;
  graphButton ("Swap", colorSwap, SIZEX(BORDER+255+15), SIZEY(MAXY+15+44)) ;

  sprintf (maxTx,"%d", rampMax);
  maxTextBox = graphTextEntry (maxTx,4,
			       SIZEX(BORDER+255+15)+1,SIZEY(MAXY-7)+1,
			       colorEndEntry);
  sprintf (minTx,"%d", rampMin);
  minTextBox = graphTextEntry (minTx,4,
			       SIZEX(BORDER+255+15)+1,SIZEY(MINY-7)+1,
			       colorEndEntry); 

  graphEntryDisable ();
  graphRedraw ();

  if (SolidLines) {
      oldmin = rampMin; oldmax = rampMax;
      colorSolidLines();
  }

  /* Disable first Xorline */
  oldmin = -100 ; oldmax = -100 ;
}

static void colorPickSlider (int box)
{
    /* if (!graphActivate(rampGraph)) return; */
    SolidLines = FALSE;
    colorDrawWindow();
    SolidLines = TRUE;

    /* This should in make XorLines at LEFT_DOWN, but for some
     * obscure reason it doesn't work @!#*&^@#&^%@!#@!&^%@# (esr)
    oldmin = rampMin; oldmax = rampMax;
    colorXorLines();
    */

  if (box == maxBox || box == minBox)
    {
      defMin =  rampMin ; defMax =  rampMax ;
      sliderBox = box;
      graphRegister (LEFT_DRAG,	colorDragSlider) ;
      graphRegister (LEFT_UP,	colorUpSlider) ;
    }
  if (box == thresBox)
    {
      defMin =  rampMin ; defMax =  rampMax ;
      sliderBox = box ;		/* for colorUpSlider */
      graphRegister (LEFT_DRAG, colorDragThres) ;
      graphRegister (LEFT_UP,   colorUpSlider) ;
    }
  else if (box == maxTextBox)
    {
      defMin =  rampMin ; defMax =  rampMax ;
      graphTextEntry (maxTx,0,0,0,0) ;
    }
  else if (box == minTextBox)
    {
      defMin =  rampMin ; defMax =  rampMax ;
      graphTextEntry (minTx,0,0,0,0) ;
    }

}

static void colorDragSlider (double x,double y)
{
  x = UNSIZEX(x) ;  y = UNSIZEY(y) ;

  if (x < BORDER) x = BORDER; 
  if (x > BORDER+255) x = BORDER+255;

  graphBoxDraw (sliderBox,GREEN,TRANSPARENT);
  if (sliderBox == maxBox)
    {
      graphBoxShift (sliderBox, SIZEX(x-5), SIZEY(MAXY-5));
      rampMax= x-BORDER ;
      graphBoxShift (thresBox,
		     SIZEX(BORDER+((rampMin+rampMax)/2)-5),
		     SIZEY(((MINY-MAXY)/2+MAXY)-5) );
      if (rampMin == rampMax)
	{ 
	  if (x<oldx)
	    x-- ;
	  else  
	    x++ ;
	  if (x < BORDER) x = BORDER+1; 
	  if (x > BORDER+255) x = BORDER+255-1;
	  rampMax = x-BORDER ;
	  graphBoxShift (sliderBox, SIZEX(x-5), SIZEY(MAXY-5));  
	  graphBoxShift (thresBox,
			 SIZEX(BORDER+((rampMin+rampMax)/2)-5),
			 SIZEY(((MINY-MAXY)/2+MAXY)-5) );
	}
      sprintf (maxTx,"%d", rampMax);
      graphBoxDraw (maxTextBox, -1, -1) ;
    }
  else
    {
      graphBoxShift (sliderBox, SIZEX(x-5), SIZEY(MINY-5)) ;
      rampMin = x-BORDER ;
      graphBoxShift (thresBox,
		     SIZEX(BORDER+((rampMin+rampMax)/2)-5),
		     SIZEY(((MINY-MAXY)/2+MAXY)-5) );
      if (rampMin == rampMax)
	{ 
	  if (x < oldx)
	    x-- ;
	  else  
	    x++ ;
	  if (x < BORDER) x = BORDER+1; 
	  if (x > BORDER+255) x = BORDER+255-1;
	  rampMin = x-BORDER ;
	  graphBoxShift (sliderBox,SIZEX(x-5),SIZEY(MINY-5)) ;
	  graphBoxShift (thresBox,
			 SIZEX(BORDER+((rampMin+rampMax)/2)-5),
			 SIZEY(((MINY-MAXY)/2+MAXY)-5) ) ;
	}
      sprintf (minTx,"%d", rampMin) ;
      graphBoxDraw (minTextBox, -1, -1) ;
    }

  colorXorLines () ;
  oldx = x; oldmin = rampMin; oldmax = rampMax;
  colorXorLines () ;

  graphGreyRamp (rampMin, rampMax) ;
}

static void colorUpSlider (double x,double y)
{
  x = UNSIZEX(x) ;  y = UNSIZEY(y) ;

  if (x < BORDER) x = BORDER; 
  if (x > BORDER+255) x = BORDER+255;

  if (rampMin == rampMax)
    { if (x < BORDER) x = BORDER+1; 
      if (x > BORDER+255) x = BORDER+255-1; }

  if (sliderBox == thresBox)
    graphBoxDraw (sliderBox,RED,TRANSPARENT) ;
  else
    graphBoxDraw (sliderBox,BLACK,TRANSPARENT) ;

  graphEntryDisable () ;
  graphRegister (LEFT_DRAG, 0) ;
  graphRegister (LEFT_UP,   0) ;

   
  if (!graphWritableColors())
    {
/*
   Erik's dotter specific stuff - we can't have that
   we need to register a stati callingGraph or so.

      graphActivate (dotGraph);
*/
      greyChanged = 1;
      graphRedraw ();
    }
  colorDrawWindow();
}

static void colorDragThres (double x,double y)
{
  int d = abs(rampMin-rampMax)/2 ;

  x = UNSIZEX(x) ;  y = UNSIZEY(y) ;
  x -= BORDER ;
  if (d > x) x = d ;
  if (255-d < x) x = 255 - d ;

  if (rampMin < rampMax)
    { rampMin = x-d ; rampMax = x+d ;
      if (rampMin == rampMax) rampMax = rampMin+1 ;
    }
 else
    { rampMin = x+d ; rampMax = x-d ;
      if (rampMin == rampMax) rampMin = rampMax+1 ;
    } 
  if (rampMin == 256 || rampMax == 256) 
    { --rampMin ; --rampMax ;}
  graphBoxDraw (thresBox,GREEN,TRANSPARENT) ;
  graphBoxShift (sliderBox,
		 SIZEX(BORDER+((rampMin+rampMax)/2)-5),
		 SIZEY(((MINY-MAXY)/2+MAXY)-5)) ;

  colorSetDisplay ();
}


static void colorUndo (void) /* restore old values */
{
  int m ;
  m = rampMin; rampMin = defMin; defMin = m;
  m = rampMax; rampMax = defMax; defMax = m;
  colorDrawWindow ();
  colorSetDisplay ();
}
static void colorSwap (void)
{
  int mem;
  mem = rampMin; rampMin = rampMax; rampMax = mem;
  colorDrawWindow ();
  colorSetDisplay ();
}


static void colorXorLines (void)
{
  /* The middle slanted lines */
  graphXorLine (SIZEX(oldmin+BORDER),
		SIZEY(MINY-7),
		SIZEX(BORDER+((oldmin+oldmax)/2)-
		((float)(oldmax-oldmin)/255.0)*20.0),
		SIZEY(((MINY-MAXY)/2+MAXY)+7)) ;
  graphXorLine (SIZEX(oldmax+BORDER),
		SIZEY(MAXY+7),
		SIZEX(BORDER+((oldmin+oldmax)/2)+
		((float)(oldmax-oldmin)/255.0)*20.0),
		SIZEY(((MINY-MAXY)/2+MAXY)-7)) ;


  /* The horizontal lines */
  if (oldmax > oldmin)
    {
      graphXorLine (SIZEX(BORDER), SIZEY(MINY-7),
		    SIZEX(oldmin+BORDER-1), SIZEY(MINY-7));
      graphXorLine (SIZEX(BORDER+oldmax+1), SIZEY(MAXY+7),
		    SIZEX(BORDER+255), SIZEY(MAXY+7));
    }
  else if (oldmax < oldmin)
    {
      graphXorLine (SIZEX(BORDER), SIZEY(MAXY+7),
		    SIZEX(oldmax+BORDER-1), SIZEY(MAXY+7)) ;
      graphXorLine (SIZEX(BORDER+oldmin+1), SIZEY(MINY-7),
		    SIZEX(BORDER+255), SIZEY(MINY-7)) ;
    }
}

static void colorSolidLines (void)
{
  /* The middle slanted lines */
  graphLine (SIZEX(oldmin+BORDER),
		SIZEY(MINY-7),
		SIZEX(BORDER+((oldmin+oldmax)/2)-
		((float)(oldmax-oldmin)/255.0)*20.0),
		SIZEY(((MINY-MAXY)/2+MAXY)+7)) ;
  graphLine (SIZEX(oldmax+BORDER),
		SIZEY(MAXY+7),
		SIZEX(BORDER+((oldmin+oldmax)/2)+
		((float)(oldmax-oldmin)/255.0)*20.0),
		SIZEY(((MINY-MAXY)/2+MAXY)-7)) ;


  /* The horizontal lines */
  if (oldmax > oldmin)
    {
      graphLine (SIZEX(BORDER), SIZEY(MINY-7),
		    SIZEX(oldmin+BORDER-1), SIZEY(MINY-7));
      graphLine (SIZEX(BORDER+oldmax+1), SIZEY(MAXY+7),
		    SIZEX(BORDER+255), SIZEY(MAXY+7));
    }
  else if (oldmax < oldmin)
    {
      graphLine (SIZEX(BORDER), SIZEY(MAXY+7),
		    SIZEX(oldmax+BORDER-1), SIZEY(MAXY+7)) ;
      graphLine (SIZEX(BORDER+oldmin+1), SIZEY(MINY-7),
		    SIZEX(BORDER+255), SIZEY(MINY-7)) ;
    }
}

static void colorEndEntry (char* text)
{
  int memMin = rampMin ,memMax = rampMax ;
  BOOL changed = FALSE ;  /* FALSE = max, TRUE = min */

  rampMin = atoi (minTx);
  rampMax = atoi (maxTx);
  if (rampMin != memMin)
    changed = TRUE;
  else if (rampMax != memMax) 
    changed = FALSE;

/* checking, whether the new value is allowed or not */

  if (rampMin == rampMax ||
      (rampMin<0 || rampMin>255) || (rampMax<0 || rampMax>255))
    {
      if (changed)
	{
	  rampMin = memMin ;
	  sprintf (minTx,"%d", rampMin) ;
	  graphBoxDraw (minTextBox, -1 , -1) ;
	}
      else if (!changed)
	{
	  rampMax = memMax ;
	  sprintf (maxTx,"%d", rampMax) ;
	  graphBoxDraw (maxTextBox, -1 , -1) ;
	}
    }
  else
    colorSetDisplay ();
}


static void colorSetDisplay (void)
{
  graphActivate (rampGraph) ;

  graphBoxShift (maxBox, SIZEX(rampMax+BORDER-5), SIZEY(MAXY-5)) ;
  graphBoxShift (minBox, SIZEX(rampMin+BORDER-5), SIZEY(MINY-5)) ;
  graphBoxShift (thresBox,
		 SIZEX(BORDER+((rampMin+rampMax)/2)-5),
		 SIZEY(((MINY-MAXY)/2+MAXY)-5) );

  sprintf (maxTx,"%d", rampMax) ;
  graphBoxDraw (maxTextBox, -1, -1) ;
  sprintf (minTx,"%d", rampMin) ;
  graphBoxDraw (minTextBox, -1, -1) ;


  colorXorLines () ;
  oldmin = rampMin ; oldmax = rampMax ;
  colorXorLines () ;
  graphGreyRamp (rampMin, rampMax);
}
/******************/
/******************/
 
