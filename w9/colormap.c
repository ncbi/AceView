/*  Last edited: Dec  4 14:20 1998 (fw) */

/* $Id: colormap.c,v 1.2 2012/01/20 01:21:35 mieg Exp $ */

#include "regular.h"

#include "graph.h"
#include "key.h"

#define MAXY 20
#define MINY 110
#define BORDER 20

static int min = 0 ;
static int max = 255 ;
static int oldmin, oldmax, defMin, defMax, oldx;
static int sliderBox, maxTextBox, minTextBox ;
static int maxBox, minBox;
static char maxTx[3];
static char minTx[3];
static Graph colorGraph = 0 ;

static void colorPickSlider (int);
static void colorDragSlider (double,double);
static void colorUpSlider (double,double);
static void colorUndo (void);
static void colorXorLines (void);
static void colorEndEntry (char*);
static void colorSwap (void);
static void colorSetDisplay (void);

static MENUOPT colorFuncButtons[] = {
  {graphDestroy,"Quit"},
  {colorUndo,"Undo"},
  {colorSwap,"Swap"},
  {0,0} };

static unsigned char rawMap[256];
static char wedge[MINY-MAXY-20][256] ;

void colorMap (void)
{
  int i,j;

  if (!graphWritableColors())
    { messout ("Pixel values are readonly: colormap is being "
	       "shared with other GRAPH package applications "
	       "if possible") ;
      return ;
    }

  if (graphActivate (colorGraph)) /* checks graphExists() first */
    { graphPop () ;
      return;
    }

  colorGraph = graphCreate (PIXEL_SCROLL, "ColorMap", 0, 0, 0.38,
			    (float)(2*BORDER+(MINY-MAXY))/900) ;
  graphRegister ( PICK, colorPickSlider );
  graphMenu (colorFuncButtons); 

  graphBoxStart();
  graphRawMaps (rawMap, 0) ;
  graphPixelsRaw ((char*)wedge, 256, MINY-MAXY-20, 256, BORDER, MAXY+10) ;
  for (j = 0 ; j < MINY-MAXY-20 ; j++)
    for (i = 0 ; i < 256 ; i++)
      wedge[j][i] = rawMap[i] ;
  graphBoxEnd();

  graphColor (BLACK);

  graphRectangle (BORDER-9,MAXY-10,BORDER+255+9,MINY+10);

  maxBox = graphBoxStart () ;
  graphLine (BORDER+max-5,MAXY-5,BORDER+max,MAXY+5);
  graphLine (BORDER+max,MAXY+5,BORDER+max+5,MAXY-5);
  graphLine (BORDER+max-5,MAXY-5,BORDER+max+5,MAXY-5);
  graphBoxEnd ();

  minBox = graphBoxStart () ;
  graphLine (BORDER+min-5,MINY+5,BORDER+min,MINY-5);
  graphLine (BORDER+min,MINY-5,BORDER+min+5,MINY+5);
  graphLine (BORDER+min-5,MINY+5,BORDER+min+5,MINY+5);
  graphBoxEnd ();
  

  graphRectangle (BORDER+255+15,MAXY-7,BORDER+255+15+34,MAXY+8);
  graphRectangle (BORDER+255+15,MINY-7,BORDER+255+15+34,MINY+8);

  graphButtons (colorFuncButtons,BORDER+255+15,MAXY+15,1);
  
  sprintf (maxTx,"%d",max);
  maxTextBox = graphTextEntry (maxTx,4,BORDER+255+16,MAXY-6, colorEndEntry);
  sprintf (minTx,"%d",min);
  minTextBox = graphTextEntry (minTx,4,BORDER+255+16,MINY-6, colorEndEntry); 

  graphEntryDisable ();
  graphRedraw ();

  defMin =  min ; defMax =  max ;
  oldmin = -100 ; oldmax = -100 ;
}

static void colorPickSlider (int box)
{
  if (box == maxBox || box == minBox)
    {
      sliderBox = box;
      graphRegister ( LEFT_DRAG, colorDragSlider );
      graphRegister ( LEFT_UP, colorUpSlider );
    }
  else if (box == maxTextBox)
    graphTextEntry (maxTx,0,0,0,0) ;

  else if (box == minTextBox)
    graphTextEntry (minTx,0,0,0,0) ;

}

static void colorDragSlider (double x,double y)
{
  if (x < BORDER) x = BORDER; if (x > BORDER+255) x = BORDER+255;

  graphBoxDraw (sliderBox,GREEN,TRANSPARENT);
  if (sliderBox == maxBox)
    {
      graphBoxShift (sliderBox,x-5,MAXY-5);
      max= x-BORDER ;
      if (min == max)
	{ if (x<oldx)
	    x--;
	else  
	  x++;
	  if (x < BORDER) x = BORDER+1; if (x > BORDER+255) x = BORDER+255-1;
	  max= x-BORDER ;
	  graphBoxShift (sliderBox,x-5,MAXY-5);  
	}
      sprintf (maxTx,"%d",max);
      graphBoxDraw (maxTextBox, -1, -1) ;
    }
  else
    {
      graphBoxShift (sliderBox,x-5,MINY-5);
      min= x-BORDER;
      if (min == max)
	{ if (x<oldx)
	    x--;
	else  
	  x++;
	  if (x < BORDER) x = BORDER+1; if (x > BORDER+255) x = BORDER+255-1;
	  min= x-BORDER;
	  graphBoxShift (sliderBox,x-5,MINY-5);  
	}
      sprintf (minTx,"%d",min);
      graphBoxDraw (minTextBox, -1, -1) ;
    }

  colorXorLines ();
  oldx = x; oldmin = min; oldmax = max;
  colorXorLines ();

  graphGreyRamp (min, max);

}

static void colorUpSlider (double x,double y)
{
  if (x < BORDER) x = BORDER; if (x > BORDER+255) x = BORDER+255;
  if (min == max)
    { if (x < BORDER) x = BORDER+1; if (x > BORDER+255) x = BORDER+255-1; }

  graphBoxDraw (sliderBox,BLACK,TRANSPARENT);
  graphEntryDisable ();
  graphRegister ( LEFT_DRAG, 0 );
  graphRegister ( LEFT_UP,   0 );
}

static void colorUndo (void)
{
  min = defMin; /* restore old values */
  max = defMax;
  colorSetDisplay ();
}

static void colorXorLines (void)
{
  graphXorLine (oldmin+BORDER,MINY-7,oldmax+BORDER,MAXY+7) ;
  if (oldmax > oldmin)
    {
      graphXorLine (BORDER,MINY-7,oldmin+BORDER-1,MINY-7);
      graphXorLine (BORDER+oldmax+1,MAXY+7,BORDER+255,MAXY+7);
    }
  else if (oldmax < oldmin)
    {
      graphXorLine (BORDER,MAXY+7,oldmax+BORDER-1,MAXY+7);
      graphXorLine (BORDER+oldmin+1,MINY-7,BORDER+255,MINY-7);
    }
}

static void colorEndEntry (char* text)
{
  int memMin=min ,memMax=max;
  BOOL changed = FALSE ;  /* FALSE = max, TRUE = min */

  min = atoi (minTx);
  max = atoi (maxTx);
  if (min != memMin)
    changed = TRUE;
  else if (max != memMax) 
    changed = FALSE;

/* checking, whether the new value is allowed or not */

  if (min == max || (min<0 || min>255) || (max<0 || max>255))
    {
      if (changed)
	{
	  min = memMin;
	  sprintf (minTx,"%d",min);
	  graphBoxDraw (minTextBox, -1 , -1 );
	}
      else if (!changed)
	{
	  max = memMax;
	  sprintf (maxTx,"%d",max);
	  graphBoxDraw (maxTextBox, -1 , -1 );
	}
    }
  else
    colorSetDisplay ();
}

static void colorSwap (void)
{
  int mem;
  mem = min; min = max; max = mem;
  colorSetDisplay ();
}

static void colorSetDisplay (void)
{
  graphBoxShift (maxBox,max+BORDER-5,MAXY-5);  
  graphBoxShift (minBox,min+BORDER-5,MINY-5);  

  sprintf (maxTx,"%d",max);
  graphBoxDraw (maxTextBox, -1, -1) ;
  sprintf (minTx,"%d",min);
  graphBoxDraw (minTextBox, -1, -1) ;

  colorXorLines ();
  oldmin = min; oldmax = max;
  colorXorLines ();

  graphGreyRamp (min, max);
}
