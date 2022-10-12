/*  file: graphxlib.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: X Low level stuff for graph package.
 *		Based on code from Chris Lee.
 * Exported functions:  graphDevActivate, 
 *			graphBoxDraw, graphRedraw, graphClipDraw
 *			graphXorLine, graphXorBox,
 *			gFontInfo
 *
 * HISTORY:
 * Last edited: May 12 13:56 1999 (edgrif)
 * * May 12 13:56 1999 (edgrif): Fixed small bug in cursor setting in graphClipDraw.
 * * Dec 15 23:23 1998 (edgrif): Added busy cursor on/off to graphClipDraw.
 * * Oct 13 10:42 1998 (edgrif): Fix initialisation bug for cursors, should
 *              be initialised to 'None', replaced #if ACEDB with function
 *              call. Fix some boring stuff about incorrect casts in rawImage.
 * This is a hacked version.
 * * Mar  30      1998 Christian Iseli: Fixed problems with MSB host versus
	LSB displays (and vice versa...)
 * * Mar  3       1998 Christian Iseli: more hacks to add the grey ramp
 	tool.
 * * Mar  2       1998 Christian Iseli: hacked rawImage to try to cope
 	with non 8-bit displays.  The problem of freeing stuff is
 	not handled yet.  Also, the grey ramp tool does not work, so
 	a best spread approach is taken...
 * * Apr 24 23:34 1996 (rd)
 * * Feb 19 12:21 1996 (daz)
 * * Jan 10 16:43 1996 (rd)
 * * Nov 14 17:41 1992 (rd): FILL_ARC, BackgroundImage (8/11/92)
 * * Aug 14 17:19 1992 (rd): COLOR_SQUARES
 * * Jul 25 12:19 1992 (mieg): TEXT_PTR_PTR
 * * Dec  4 13:55 1991 (mieg): ACEDB_COLOR to force isMono = FALSE
 * * Oct 14 15:56 1991 (rd): added clipping of all vectors, except
	only test the centre (top left) of circles, points, text.
 *-------------------------------------------------------------------
 */

/* $Id: graphxlib.c,v 1.15 2017/05/24 16:32:57 mieg Exp $ */

#define  XLIB_DISPLAY
#include "regular.h"

#include "graphxt_.h"

#include <X11/Xaw/Scrollbar.h>
#include <X11/cursorfont.h>

extern char *acefont ;		/* from graphcon.c */

Colormap getGlobalColorMap (void) ; /* for graphxt.c */

static void fontInit (Display * display) ;
static void colorInit (Display * display, Screen * screen) ;
static void stippleInit (Window window) ;
Visual *getVisual (void);

/******* statics for drawing - set in graphDevActivate() *********/

static Display	*display ;
static Window   window ;
static GC	gc, gcInvert, gcStipple ;

#define DWG display,window,gc  /* used in all drawing calls */
#define DG  display,gc         /* used in all GC change calls */

static unsigned long gcValuesMask ;
static XGCValues gcValues;

static BOOL isMono;

/* daz stuff for clean colormaps */
BOOL cmapDebug = FALSE ;
int preferredVisual = -1;  /* Can be set to select a non-default visual */

Colormap globalColorMap = 0 ;	/* also used in graphramp.c */
static	int	visualDepth; /* Depth in planes of the visual that is used */

static  int	bits_per_pixel; /* Pixmap info.  */
static	int	byte_order; /* More pixmap info.  */
static unsigned long greyPixels[256];
extern int greyChanged;


extern int supressXError;
extern BOOL installOwnColorMap;

/* stuff from Theron Friedman (tzf@cs.arizona.edu) to prevent
   program quitting in Openlook/MacX
*/
Atom myExitAtom ;		/* used in graphLoop */
static Atom exitHintAtom ;

Pixmap GreyPixmap; 
unsigned char grey_data[] = {0x55, 0xaa, 0x55, 0xaa, 0x55, 0xaa, 0x55, 0xaa};


/*******************************************************/

static void xlibInit (void)
{
  int screen_num = DefaultScreen(display) ;
  int myclass ;
  int mydepth;

  gc = XCreateGC (display, RootWindow(display,screen_num), 0, 0) ;
  XSetLineAttributes (display, gc, 0, LineSolid, CapRound, JoinRound) ;
  gcInvert = XCreateGC (display, RootWindow(display,screen_num), 0, 0) ;
  XSetFunction (display, gcInvert, GXinvert) ;
  fontInit (display) ;

  mydepth = DefaultDepth (display, screen_num) ; 
  myclass = getVisual()->class; 

  /* for subgraph borders */
  GreyPixmap = XCreatePixmapFromBitmapData(display, 
					   RootWindow(display, screen_num),
					   (char *) grey_data, 4, 4,
					   1, 0, mydepth);


  /* NOTE the 'implied' acedb dependency here...                             */
  isMono = (getenv ("ACEDB_MONO") != 0 || visualDepth == 1 || myclass < 2) ;

  if(getenv("ACEDB_COLOUR") || getenv("ACEDB_COLOR")) isMono = FALSE ;

  if (isMono)
    { stippleInit (RootWindow (display, screen_num)) ;
      gcStipple = XCreateGC (display, RootWindow(display,screen_num), 0, 0) ;
      XSetFillStyle (display, gcStipple, FillOpaqueStippled) ;
    }

	/* tzf quit-safety stuff */
  exitHintAtom = XInternAtom (display, "WM_PROTOCOLS", FALSE) ;
  myExitAtom = XInternAtom (display, "WM_DELETE_WINDOW", FALSE) ;

}

void graphDevActivate (int show) 
{
  display = XtDisplay (gDev->popup) ;
  window = XtWindow (gSubDev->base) ;
  if (!gc)
    xlibInit () ;
  XChangeProperty (display, XtWindow(gDev->popup), 
		   exitHintAtom, XA_ATOM, 32, 
		   PropModeReplace, (UCHAR*)&myExitAtom, 1) ;
}


/********* next the fonts - needs work to clean up *********/
/* Stuff like -acefont needs to go...just use the normal X font flag.        */

/* font policy: use simple fonts so that they should all be there
   load them all in at beginning
   fontInit assumes default font is first in list
   note: this is a display specific object - so fontInit needs the
     display to have been opened.
*/

enum FONTTYPE {F8x13, F5x8, F6x9, F6x10, F6x12, F9x15, 
		 F8x13I, F5x8I, F6x9I, F6x10I, F6x12I, F9x15I, 
		 F8x13B, F5x8B, F6x9B, F6x10B, F6x12B, F9x15B, 
                 G13, G6, G8, G10, G12, G15,
		 NUM_FONTS} ;
static char fontName[NUM_FONTS][120] ;

/*      {"8x13","5x8","6x9","6x10","6x12","9x15"} ; */

#define Fdefault F8x13
static XFontStruct *fns[NUM_FONTS] ;
static int fontdy ;  /* is simply the ascent used to offset base
			line to make drawtext start at upper left */
static int fontWidth, fontHeight ;



static void getFontNames(void)
  {
  char *cp ;
  int  i = NUM_FONTS, curr_level, level ;

  freeinit() ;						    /* Just to make sure. */

  /* ACEDB-GRAPH INTERFACE: acedb can set a font file of its own to be the   */
  /* input for this routine, if it succeeds it alters the freeNNN streamlevel*/
  /* so we check to see if this has happened. This seems a bit of a hack,    */
  /* but this makes it clear that the function does change the stream level. */
  level = curr_level = freeCurrLevel() ;

  if (getGraphAcedbXFont() != NULL) (getGraphAcedbXFont())(&level) ;

  /* This is the default action, this uses a set of hard coded fonts.        */
  if (level == curr_level)
    {
    level = freesettext("8x13 5x8 6x9 6x10 6x12 9x15 "
			"8x13 5x8 6x9 6x10 6x12 9x15 "
			"8x13bold 5x8 6x9 6x10 6x12 9x15bold "
			"8x13 5x8 6x9 6x10 6x12 9x15 "
			, "") ;
    }


  while (i--)
    *fontName[i] = 0 ;

  i = 0 ;
  while (freecard(level))
    while ((cp = freeword()))
      strncpy(fontName[i++], cp,100) ;


  /* This is not strictly an ACEDB dependency, it's a pity about the name    */
  /* of the font, it would have been better to stick to the standard X way   */
  /* of naming fonts.                                                        */
  if (acefont) strncpy(fontName[0], acefont, 100) ;

  return ;
  }


static void fontInit (Display* display)
  {
  int i ;
  static int isDone = FALSE ;

  if (isDone == FALSE)
    {
    isDone = TRUE ;

    getFontNames() ;

    if (!(fns[0] = XLoadQueryFont (display, fontName[0])))
      {
      strncpy(fontName[0], "6x13",100) ;
      if (!(fns[0] = XLoadQueryFont (display, fontName[0])))
	    messcrash ("Can't load default font %s",fontName[0]) ;
      }
    for (i = 1 ; i < NUM_FONTS ; ++i)
      {
      if (!(fns[i] = XLoadQueryFont (display, fontName[i])))
	{
	fprintf(stderr,"Font %s does not exist - using default\n", fontName[i]) ;
	fns[i] = fns[0] ;
	}
      }

    }

  return ;
  }


static int height2font (int height, int format)
{
  static float oldHeight = 1.;
  static int oldFormat = PLAIN_FORMAT ;
  if (height >= 0)
    oldHeight = height;
  else
    height = oldHeight ;
  if (format >= 0)
    oldFormat = format ;
  else
    format = oldFormat ;
  switch (height)
    {
    case 9999:			/* old default font */
      switch (format)
	{
	case PLAIN_FORMAT:
	case FIXED_WIDTH:
	  return F8x13;
	case ITALIC:
	  return F8x13I;
	case BOLD:
	  return F8x13B;
	case GREEK:
	  return G13;
	}
    case 0:			/* default font */
      switch (format)
	{
	case PLAIN_FORMAT:
	case FIXED_WIDTH:
	  return F9x15 ;
	case ITALIC:
	  return F9x15I ;
	case BOLD:
	  return F9x15B ;
	case GREEK:
	  return G15;
	}
       break ;
    case 1: case 2: case 3: case 4: case 5: case 6:
      switch (format)
	{
	case PLAIN_FORMAT:
	case FIXED_WIDTH:
	  return F6x10;
	case ITALIC:
	  return F6x10I ;
	case BOLD:
	  return F6x10B ;
	case GREEK:
	  return G6;
	}
    case 7: case 8:
      switch (format)
	{
	case PLAIN_FORMAT:
	case FIXED_WIDTH:
	  return F5x8 ;
	case ITALIC:
	  return F5x8I ;
	case BOLD:
	  return F5x8B ;
	case GREEK:
	  return G8;
	}
    case 9: case 10:
      switch (format)
	{
	case PLAIN_FORMAT:
	case FIXED_WIDTH:
	  return F6x10;
	case ITALIC:
	  return F6x10I ;
	case BOLD:
	  return F6x10B ;
	case GREEK:
	  return G10;
	}
      return F6x10;
    case 11: case 12:
      switch (format)
	{
	case PLAIN_FORMAT:
	case FIXED_WIDTH:
	  return F6x12;
	case ITALIC:
	  return F6x12I ;
	case BOLD:
	  return F6x12B ;
	case GREEK:
	  return G12;
	}
    case 13: case 14:
      switch (format)
	{
	case PLAIN_FORMAT:
	case FIXED_WIDTH:
	  return F8x13;
	case ITALIC:
	  return F8x13I ;
	case BOLD:
	  return F8x13B ;
	case GREEK:
	  return G13;
	}
    default:
      switch (format)
	{
	case PLAIN_FORMAT:
	case FIXED_WIDTH:
	  return F9x15 ;
	case ITALIC:
	  return F9x15I ;
	case BOLD:
	  return F9x15B ;
	case GREEK:
	  return G15;
	}
    }
  return F8x13;
}


BOOL gFontInfo (int height, int* w, int* h)
	/* note that this depends on an integer height, and so
	   is only reliable in the long term for height 0 */
  {
  XFontStruct* fnt = fns[height2font (height,0)] ;

  if (!fnt)
    {
    fontInit (XtDisplay (root_widget)) ;
    if (!(fnt = fns[height2font (height,0)]))
      return FALSE ;
    }

  if (w) *w = fnt->max_bounds.width;
  if (h) *h = fnt->ascent + fnt->descent;

  return TRUE ;
  }



static int lastFontHeight ;


static XFontStruct *fontSet (int height)
{
  XFontStruct *fnt = fns[height2font (height,-1)] ;

  if (!fnt)
    fnt = fns[height2font(0,0)] ; 

  XSetFont (DG, fnt->fid) ;
  fontdy = fnt->ascent ;
  fontWidth = fnt->max_bounds.width ;
  fontHeight = fnt->ascent + fnt->descent ;

  lastFontHeight = height ;
  return fnt ;
}

static XFontStruct *fontSetFormat (int format) /* unfinished */
{
  XFontStruct *fnt = fns[height2font (-1, format)] ;

  if (!fnt)
    fnt = fns[height2font(0,0)] ;

  XSetFont (DG, fnt->fid) ;
  fontdy = fnt->ascent ;
  fontWidth = fnt->max_bounds.width ;
  fontHeight = fnt->ascent + fnt->descent ;

  return fnt ;
}

XFontStruct * getDefaultFont(void)
/* this little function is needed to set fonts in text editor. see graphxt.c */
{ return fns[Fdefault];
}

/********** next color control - screen specific ***********/

static unsigned long colorInd[256] ; /* really: a power of 2 larger than NUM_TRUECOLORS] */
/* old version - X11 colours
   static char* colorName[] = 
   {"white","black","lightgray","dimgray",
   "red","green","blue",
   "yellow","cyan","magenta",
   "pink","yellowgreen","skyblue",
   "violet","darkgreen","cadetblue"} ;
*/

static BOOL nearestColor (Colormap cmap, XColor *col)
{ 
  int i, imin, dis, dismin, x ;
  int r = col->red >> 8 ;
  int g = col->green >> 8 ;
  int b = col->blue >> 8 ;
  static XColor *cols = 0 ;
  static int npix = 256 ;	/* tosh1 has 4 bit graphics! */

  return 0 ;
  if (!cols)
    { if (DefaultDepth (display, DefaultScreen(display)) < 8)
	{ npix = 1 << DefaultDepth (display, DefaultScreen(display)) ;
	  fprintf (stderr, "screen depth problem: only %d pixel values available\n", npix) ;
	}
      cols = (XColor*) messalloc (npix*sizeof(XColor)) ;
      for (i = 0 ; i < npix ; ++i)
	cols[i].pixel = i ;
    }

  XQueryColors (display, cmap, cols, npix) ;
  dismin = 1 << 20 ;
  for (i = imin = 0 ; i < npix ; ++i)
    { x = r - (cols[i].red >> 8) ; dis = x*x ;
      x = g - (cols[i].green >> 8) ; dis += x*x ;
      x = b - (cols[i].blue >> 8) ; dis += x*x ;
      if (dis < dismin)
	{ imin = i ;
	  dismin = dis ;
	}
    }
		/* try to alloc closest - to get a read lock */
  col->red = cols[imin].red ;
  col->green = cols[imin].green ;
  col->blue = cols[imin].blue ;
  col->pixel = imin;

  return TRUE ;
}


static void colorInit (Display *display, Screen *screen)
{
  int i, red,green,blue,grey ;
  XColor color ;
  Colormap cmap = globalColorMap ;

  for (i = 0 ; i < NUM_TRUECOLORS ; ++i)
    { 
      graphGetColourAsInt (i, &red, &green, &blue, &grey) ;
      color.red = red << 8 ;
      color.green = green << 8 ;
      color.blue = blue << 8 ;
      color.flags = DoRed | DoGreen | DoBlue ;
      if (XAllocColor (display, cmap, &color) ||
	      nearestColor (cmap, &color))
	    colorInd[i] = color.pixel ;
	  else			/* messout here creates a crash */
	    { fprintf (stderr, "Can't allocate color %d %d %d\n", 
		       red, green, blue) ;
	      colorInd[i] = WhitePixelOfScreen(screen) ;
	    }
	}
}


/************ stipples for monochrome systems ***********/

static unsigned short whiteMask[] = {0xffff,0xffff,0xffff,0xffff} ;
static unsigned short lightMask[] = {0x7f7f,0xfbfb,0xdfdf,0xfefe,
				     0xf7f7,0xbfbf,0xfdfd,0xefef} ;
static unsigned short mixMask[]   = {0x7777,0xdddd,0xbbbb,0xeeee} ;
static unsigned short halfMask[]  = {0xaaaa,0x5555,0xaaaa,0x5555} ;
static unsigned short darkMask[]  = {0x8888,0x2222,0x4444,0x1111} ;
static unsigned short blackMask[] = {0x0000,0x0000,0x0000,0x0000} ;
static Pixmap  whiteStipple, lightStipple, mixStipple,
  	       halfStipple,  darkStipple,  blackStipple ;
static Pixmap* stipple[] =    /* 32 of them : should be a power of 2 larger than NUM_TRUECOLORS] */
  { &whiteStipple, &blackStipple, &lightStipple, &darkStipple,
    &halfStipple,  &halfStipple,  &halfStipple,
    &mixStipple,   &mixStipple,   &mixStipple,
    &lightStipple, &lightStipple, &lightStipple,
    &darkStipple,  &darkStipple,  &darkStipple,
    &lightStipple, &lightStipple, &lightStipple,
    &lightStipple, &lightStipple, &lightStipple,
    &darkStipple,  &halfStipple,  &lightStipple,
    &darkStipple,  &halfStipple,  &lightStipple,
    &halfStipple,  &lightStipple, &halfStipple,  &halfStipple
  } ;

static void stippleInit (Window window)
{
  if (getenv("ACEDB_GRAY") || 
      getenv("ACEDB_GREY") )  /* Some screens seem inverted */ 
  {
    blackStipple = XCreateBitmapFromData (display, window,
					  (char*) whiteMask, 16, 4) ;
    darkStipple = XCreateBitmapFromData (display, window,
					  (char*) lightMask, 16,8) ;
    halfStipple   = XCreateBitmapFromData (display, window,
					  (char*) mixMask, 16, 4) ;
    mixStipple  = XCreateBitmapFromData (display, window,
					  (char*) halfMask, 16, 4) ;
    lightStipple  = XCreateBitmapFromData (display, window,
					  (char*) darkMask, 16, 4) ;
    whiteStipple = XCreateBitmapFromData (display, window,
					  (char*) blackMask, 16, 4) ;
  }
  else
  {
    whiteStipple = XCreateBitmapFromData (display, window,
					  (char*) whiteMask, 16, 4) ;
    lightStipple = XCreateBitmapFromData (display, window,
					  (char*) lightMask, 16,8) ;
    mixStipple   = XCreateBitmapFromData (display, window,
					  (char*) mixMask, 16, 4) ;
    halfStipple  = XCreateBitmapFromData (display, window,
					  (char*) halfMask, 16, 4) ;
    darkStipple  = XCreateBitmapFromData (display, window,
					  (char*) darkMask, 16, 4) ;
    blackStipple = XCreateBitmapFromData (display, window,
					  (char*) blackMask, 16, 4) ;
  }

						       
  if (!whiteStipple || !lightStipple || !mixStipple || 
      !halfStipple || !darkStipple || !blackStipple)
    messcrash ("Could not make stipples for monochrome X graphics") ;
}

/**************** pixel images ****************/

static UCHAR greyMap[256] ;
static XColor colors[256] ;

BOOL gWriteColors = FALSE ;

/*
 * Set up the pixel values for the greyramp
 * Be careful to avoid the ones we need for our standard colors.
 */

static
void insertGreyRamp(unsigned long pixels[])
{
    int i,j, red,green,blue,grey ;
    BOOL pixelTaken[256];
    XColor color;
    static XColor *cols = 0 ;

  /*
   * First go through the standard colors.  Find their nearest
   * analogue in the default colourmap, and use that pixel value to
   * allocate their color in this colormap.  
   */
  Colormap defCMap = DefaultColormap(display,DefaultScreen(display));

  /*
   * Do a blanket copy of the colors in the default colormap into our
   * own colormap.  This way we minimize disruption to the default
   * colormap when our colormap is installed.
   */

  cols = (XColor*) messalloc (256*sizeof(XColor)) ;
  for (i = 0 ; i < 256 ; ++i) {
    cols[i].pixel = i ;
  }
  XQueryColors (display, defCMap, cols, 256) ;
  XStoreColors(display,globalColorMap,cols,256);

  for(j=0;j<256;j++) {
    pixelTaken[j] = FALSE;
  }
  for (i = 0 ; i < NUM_TRUECOLORS ; ++i)
    { 
      graphGetColourAsInt (i, &red, &green, &blue, &grey) ;
      color.red = red << 8 ;
      color.green = green << 8 ;
      color.blue = blue << 8 ;
      color.flags = DoRed | DoGreen | DoBlue ;
      if (nearestColor (defCMap, &color)) {
	colorInd[i] = color.pixel ; /* Use the matching pixel value */
	pixelTaken[color.pixel] = TRUE;

	/* Set the color anyway to what we wanted, since we have a R/W
	 * colormap anyway.  Making sure the pixel value is the same as one
	 * with a color similar in the default colormap just minimises
	 * distress when the colormap changes.
	 */
	graphGetColourAsInt (i, &red, &green, &blue, &grey) ;
	color.red = red << 8 ;
	color.green = green << 8 ;
	color.blue = blue << 8 ;
	XStoreColor(display,globalColorMap,&color);
      }
      else			/* messout here creates a crash */
	{ fprintf (stderr, "Can't allocate color %d %d %d\n", 
		   red, green, blue) ;
	  colorInd[i] = WhitePixelOfScreen(XtScreen (root_widget)) ;
	}
    }

    if (cmapDebug) printf("**\n");
    for (i = 0, j = 255; i < 256; i++, j--) 
      { 
	while (j >=0 && pixelTaken[j]) {
	  j--;
	}
	/* This is the grey value that we want to allocate */
	colors[i].red = colors[i].green = colors[i].blue = i << 9;
	colors[i].flags = DoRed | DoGreen | DoBlue;

	if (j >= 0 && i < 128)
	  {
	    colors[i].pixel = pixels[j];
	    greyMap[2 * i] = colors[i].pixel;
	    greyMap[2 * i + 1] = colors[i].pixel;
	    if (cmapDebug) printf(" %d ",(int)pixels[j]);
	  }
      }
    if (cmapDebug) printf("\n**\n");
}

/*
 * Handles the case where we are sharing the default colormap
 */

static void
shareColorMap (Colormap cmap)
{
    int i;
      unsigned long planes[1] ;
  unsigned long pixels[256];

    if (cmapDebug) printf("** trying to share default colormap \n");
    if (XAllocColorCells (display, cmap, FALSE, planes, 0, pixels, 256)) {
	/* allocate our own read-write cells */
	if (cmapDebug) printf("** have allocated 256 r/w cells !\n");
	for (i = 0 ; i < 256 ; ++i) { 
		colors[i].red = colors[i].green = colors[i].blue = i << 9 ;
	      colors[i].flags = DoRed | DoGreen | DoBlue ;
	      colors[i].pixel = pixels[i] ;
	      if (i < 128)
		greyMap[2*i] = greyMap[2*i+1] = pixels[i] ;
        }
	  XStoreColors (display, cmap, colors, 256) ;
	  gWriteColors = TRUE ;
    } else {
	/* We can't have RW cells, so grab readonly cells */
	if (cmapDebug) printf("** am forced to use R/O cells \n");
	for (i = 0 ; i < 256 ; ++i) { 
		colors[i].red = colors[i].green = colors[i].blue = i << 9 ;
	      colors[i].flags = DoRed | DoGreen | DoBlue ;

	      XAllocColor(display,cmap,colors+i);
	      if (i < 128)
		greyMap[2*i] = greyMap[2*i+1] = colors[i].pixel ;
        }
	  gWriteColors = FALSE ;
    }
    colorInit (display, XtScreen (root_widget));
}

Visual *getVisual (void)	/* also allocates colors */
{
  Colormap cmap ;
  XColor dummyC ;
  unsigned long dummyULong ;
  int dummyInt, *data ;
  unsigned long bytesAfter ;
  static BOOL isFirst = TRUE ;
  static Visual *visual ;
  Window rootWin = RootWindow(display,DefaultScreen(display));
  BOOL cmapAllocated = FALSE;
  BOOL cmapCorrupt = FALSE;
  BOOL cmapPreAllocated = FALSE;
  BOOL cmapChangeable = FALSE;
  int screen_num;
  unsigned long pixels[256], planes[1] ;
  Atom colorMapAtom;
  Atom dummyType;

  screen_num = DefaultScreen(display) ;
  /* mydepth = DefaultDepth(display, screen_num) ; */

  if (isFirst)
    {
      if (!display)
	messcrash ("Color maps requested before a window was created") ;

      isFirst = FALSE ;

	  if (preferredVisual != -1) { /* Select a specific visual */
		XVisualInfo	vi,*vip;
		int n;
		vi.visualid = preferredVisual;
		vip = XGetVisualInfo(display,VisualIDMask,&vi,&n);
		if (vip == NULL) {
			fprintf(stderr,"Can't find visual with id = %d",preferredVisual);
			exit(1);
		}
		visual = vip->visual;
		visualDepth = vip->depth;
		XFree(vip);
	  } else {
		visual = DefaultVisual (display, screen_num) ;
		visualDepth = DefaultDepth (display, screen_num) ;
	  }
	  if (cmapDebug) {
		printf("** visual id    = %x \n", assInt (visual->visualid));
		printf("** visual class = %d \n", visual->class);
		printf("** visualDepth  = %d \n", visualDepth);
	  }

	  switch(visual->class) {
		case DirectColor:
			cmapChangeable = TRUE; /* Can we change the colormap */
			break;
		case TrueColor:
			cmapChangeable = FALSE;
			break;
		case PseudoColor:
			cmapChangeable = TRUE;
			break;
		case GrayScale:
			cmapChangeable = TRUE;
			break;
		case StaticColor:
			cmapChangeable = FALSE;
			break;
		case StaticGray:
			cmapChangeable = FALSE;
			break;
		default:
			fprintf(stderr,"Unknown visual type = %d",visual->class);
			exit(1);
	  }

	  /* If we can't change the colormap, then there is no point doing
	   * anything fancy with colormaps, so we can allocate our colors in
	   * the default colormap and then return the visual
	   */
	  if (!cmapChangeable) {
	    int i, nb;
	    XPixmapFormatValues *pixinfo;
	    if (cmapDebug) {
	      fprintf(stderr,"** cmap not changeable, no point allocating or "
		      " playing with the colormap ** \n");
	    }
	    /* Still initialise the standard acedb colors*/
	    globalColorMap = DefaultColormap(display,DefaultScreen(display));
	    colorInit(display, XtScreen (root_widget)); 
	    gWriteColors = FALSE;
	    /* Get Pixmap info.  */
	    pixinfo = XListPixmapFormats(display, &nb);
	    for (i = 0; i < nb; i++)
	      if (bits_per_pixel < visualDepth
		  || (pixinfo[i].bits_per_pixel < bits_per_pixel
		      && pixinfo[i].bits_per_pixel >= visualDepth))
		bits_per_pixel = pixinfo[i].bits_per_pixel;
	    byte_order = ImageByteOrder(display);
	    if (cmapDebug)
	      fprintf(stderr,"** Selected %d bit pixmaps, ByteOrder = %d\n",
		      bits_per_pixel, byte_order);
	    for (i = 0; i < 256; i++) 
	      {
		colors[i].flags = DoRed | DoGreen | DoBlue ;
		colors[i].red = colors[i].green = colors[i].blue = i << 9;

		if (i < 128)
		  greyMap[2 * i] = greyMap[2 * i + 1] = i;
		if (!XAllocColor (display, globalColorMap, &colors[i]))
		fprintf (stderr, "Can't allocate grey %d\n", i);
		greyPixels[i] = colors[i].pixel;
	      }
	    return visual;
	  }

	  /* cmap is changeable, but does it have enough planes ? */
	  if (visualDepth<8) {
		/* There are not enough planes to do a greyramp and so
		 * we should just init the colormap and continue from here 
		 */
		if (cmapDebug) {
			printf("** too few planes for greyramp - just doing basic colors "
				"** \n");
		}
		globalColorMap = DefaultColormap(display,DefaultScreen(display));
		colorInit(display, XtScreen (root_widget));
		gWriteColors = FALSE;
		return visual;
      }	

		/*
		 * Loop around until we create the colormap correctly.
		 */

      do {
	if (cmapDebug) printf("** Checking for existing CMAP atom\n"); 
	colorMapAtom = XInternAtom (display, "ACEDB_ColorMap", 1) ;

	if (colorMapAtom == None || cmapCorrupt) 
	  { if (cmapDebug) printf("** no atom, or corrupt atom. \n"); 
	    /* We will have to create our own colormap */

	    /* Create the colormap atom */
	    if (!cmapCorrupt) 
	      { if (cmapDebug) printf("** create new CMAP atom\n");
		colorMapAtom = XInternAtom (display, "ACEDB_ColorMap",0);
		if (colorMapAtom == 0) 
		  messcrash ("Failed to create X server atom 'ACEDB_ColorMap'");
	      }

			/*
			 * Create a color map
			 */
	    if (cmapDebug) printf("** create colormap\n");
	    if (!installOwnColorMap) {
		/*
		 * If we have got to the stage where we would have to 
		 * create our own colormap, don't do it - break out and
		 * use the default colormap instead.
		 */
		 break;
	    }
	    cmap = XCreateColormap (display, 
				    rootWin,
				    visual,AllocNone);

	    if (cmap == 0)
	      messcrash ("Failed to create X color map") ;
	    if (cmapDebug) 
	      printf ("** CMAP %x created. Storing atom. \n", assInt(cmap)) ;

	    	     /* Store the colormap value in the colormap atom */
	    XChangeProperty (display,
				 rootWin,
				 colorMapAtom,
				 colorMapAtom,
				 8,
				 PropModeReplace, 
				 (unsigned char *)(&cmap),
				 sizeof(cmap));
	    cmapAllocated = TRUE;
	    if (cmapDebug) printf("** CMAP atom allocated.\n");
	  } 
	else 
	  {		/* Using an existing colormap */
	    if (cmapDebug) printf("** retrieving existing colormap atom\n");
	    if (Success != 
		XGetWindowProperty (display,
				    rootWin,
				    colorMapAtom,
				    0, /* offset 0 */
				    sizeof(cmap)/4,
				    0, /* don't delete */
				    AnyPropertyType,
				    &dummyType, /* Actual type */
				    &dummyInt, /* actual format */
				    &dummyULong, /* n items */
				    &dummyULong, /* bytes after */
				    (unsigned char **)&data)) 
	      {	data = 0;
		if (cmapDebug) printf("** Atom corrupt.\n");
	      }
	    if (!data)
	      { /* Data is corrupt - destroy atom and try again */
		if (cmapDebug) printf("** CMAP atom data corrupt.\n");
		XFlush(display);
		supressXError = 1;
		XDeleteProperty(display,rootWin,colorMapAtom);
		XFlush(display);
		supressXError = 0;
		cmapCorrupt = TRUE;
	      } 
	    else 
	      {	
		/* data points to 4 bytes containing the colormap id */
		cmap = *data; 
		if (cmapDebug) printf("** existing colormap atom = %x\n",assInt(cmap)) ;
		/* Check that colormap is ok */
		XFlush(display);
		supressXError = 1;
		if (!XLookupColor(display,cmap,"white",&dummyC,&dummyC))
		  { cmapCorrupt = TRUE;
		    if (cmapDebug) printf("** CMAP corrupt.\n");
		  } 
		else 
		  { cmapAllocated = TRUE;
		    if (cmapDebug) printf("** CMAP reused.\n");
		    cmapPreAllocated = TRUE;
		  }
		XFlush(display);
		supressXError = 0;
	      }
	    if (cmapCorrupt) { 
			if (cmapDebug) printf("** Corrupt colormap is freed.\n");
	    }
	    if (data) XFree((caddr_t)data);
	  }
      } while (!cmapAllocated) ;

      if (!installOwnColorMap) {
	/* We are sharing the default colormap with others */
	globalColorMap = DefaultColormap(display,DefaultScreen(display));
	shareColorMap(globalColorMap);
	return visual;
      } else {
	  globalColorMap = cmap;
      }

      /* By now we have an allocated colormap - either a preexisting
       * one, or a new one that we have created.
       */

      if (cmapDebug) printf("** CMAP =%x\n",assInt(cmap));
      cmapCorrupt = FALSE;

      if (cmapPreAllocated) 
	{ if (Success != XGetWindowProperty(display,
					    rootWin,
					    colorMapAtom,
					    sizeof(cmap)/4, /* offset */
					    sizeof(pixels[0]) * 256 / 4,
					    0, /* don't delete */
					    AnyPropertyType,
					    &dummyType, /* Actual type */
					    &dummyInt, /* actual format */
					    &dummyULong, /* n items */
					    &bytesAfter, /* bytes after */
					    (unsigned char **)&data))
	    { if (cmapDebug) printf("** Pixels corrupt.\n");
	      cmapCorrupt = TRUE;
	    } 
	  else 
	    { if (bytesAfter != 0) messcrash ("Colormap atom corrupted") ;
	      if (data) 
		{ if (cmapDebug) printf("** Pixels retrieved ok.\n") ;
		  memcpy (pixels, data, sizeof(pixels[0]) * 256) ;
	      
		  insertGreyRamp(pixels);
		      /* if (cmapDebug) printf(" %d ",pixels[i]); */
		}
		  /* if (cmapDebug) printf("\n"); */
	      else
		{ if (cmapDebug) 
		    printf("** Pixels retrieved corrupt.\n");
		  cmapCorrupt = TRUE;
		}
	    }
	}

      if (!cmapPreAllocated || cmapCorrupt)
	{ 
  	  XAllocColorCells (display, cmap, FALSE,planes,0,pixels,256) ;
	  insertGreyRamp(pixels);
	  /* if (cmapDebug) printf("\n"); */
	  XStoreColors (display, globalColorMap, colors, 256);

	  /* Append pixel values to colormap atom */

	  if (cmapDebug) printf("*** Appending pixel values to property\n");
	  XChangeProperty (display,
			       rootWin,
			       colorMapAtom,
			       colorMapAtom,
			       8,
			       PropModeAppend, 
			       (unsigned char *)(pixels),
			       sizeof(pixels[0])*256);
	}
      gWriteColors = TRUE ;
    }
  return visual ;
}

Colormap getGlobalColorMap (void)
{ 
  if (!globalColorMap)
    { display = XtDisplay (root_widget) ;
      getVisual () ;
    }

  return globalColorMap ;
}

BOOL graphRawMaps (unsigned char *forward, int *reverse) 
{
  int i ;

  if (!getVisual())
    return FALSE ;

  if (reverse)
    for (i = 0 ; i < 256 ; ++i)
      reverse[i] = 0 ;
  for (i = 0 ; i < 256 ; ++i)
    { if (forward)
	forward[i] = greyMap[i] ;
      if (reverse)
	reverse[greyMap[i]] = i ;
    }
  return TRUE ; 
}

BOOL graphWritableColors (void)
{
  getVisual () ;
  return gWriteColors ;
}

void graphColorMaps (float *red, float *green, float *blue) /* get */
{
  int i ;
  float fac = 1.0 / 0xffff ;

  for (i = 0 ; i < 256 ; ++i)
    { red[i] = fac * colors[i/2].red ;
      green[i] = fac * colors[i/2].green ;
      blue[i] = fac * colors[i/2].blue ;
    }
}

BOOL graphSetColorMaps (unsigned char *red,
			unsigned char *green,
			unsigned char *blue)
{
  int i ;

  if (visualDepth > 8) {
    if (red == green && red == blue)
      for (i = 0 ; i < 256 ; ++i)
	colors[i].pixel = greyPixels[red[i] >> 1];
    else
      for (i = 0 ; i < 256 ; ++i) {
	colors[i].red = red[i] << 8 ;
	colors[i].green = green[i] << 8 ;
	colors[i].blue = blue[i] << 8 ;
	if (!XAllocColor (display, globalColorMap, &colors[i]))
	  fprintf (stderr, "Can't allocate grey %d\n", i);
      }
    return TRUE ;
  }

  if (!getVisual() || !gWriteColors) /* () were missing, mieg 22 9 04 */
    return FALSE ;
  
  for (i = 0 ; i < 256 ; ++i)
    { colors[i].red = red[i] << 8 ;
      colors[i].green = green[i] << 8 ;
      colors[i].blue = blue[i] << 8 ;
    }

  XStoreColors(display, globalColorMap, colors, 256) ;
  return TRUE ;
}

/**************** function from suzi ****************/

static void pixelRubber (char *pixels, UCHAR *data, int padImage,
			 int wPixels, int hPixels, int lenPixels,
			 int wImage, int hImage)
{
  int i, j, jRow ;
  UCHAR *row ;
  float wFac, hFac ;
  int *sample = 0 ;

  hFac = (hPixels-1) / (float) (hImage-1) ;
  if (wPixels != wImage)
    { wFac = (wPixels-1) / (float)(wImage-1) ;
      sample = (int*) messalloc (sizeof(int) * wImage) ;
      for (i = 0 ; i < wImage ; ++i)
	sample[i] = 0.5 + i * wFac ;
    }
  for (j = 0 ; j < hImage ; ++j) 
    { jRow = 0.5 + j * hFac ;
      row = (UCHAR *)(pixels + lenPixels*jRow) ;
      if (wPixels == wImage)
	for (i = 0 ; i < wImage ; i++)
	  *data++ = greyMap[*row++] ;
      else
	for (i = 0 ; i < wImage ; i++)
	  *data++ = greyMap[row[sample[i]]] ;
      data += padImage ;
    }
  messfree (sample) ;
}

/***********************/

static XImage *pixelsImage (char *pixels,
			    int wPixels, int hPixels, int lenPixels,
			    int wImage, int hImage)
{
  int lenImage, padImage ;
  Visual *visual ;
  char *data ;
  XImage *xim ;

  if (!(visual = getVisual ()) || /* checks screen depth, visual class */
      hImage < 2 || wImage < 2)
    return 0 ;

  if (!gSubDev->images)
    gSubDev->images = assCreate () ;

  if (assFind (gSubDev->images, pixels, &xim)) 
    {
      if (xim->width == wImage && 
	  xim->height == hImage)
	return xim ;
      else
	{
	  XDestroyImage (xim) ;
	  /* DAZ */
	  if (!assRemove(gSubDev->images,pixels))
	    {
	      messcrash("Expected to be able to remove image from array");
	    }
	}
    }

  lenImage = wImage ; 
  while (lenImage % 4)
    ++lenImage ;
  padImage = lenImage - wImage ;
  data = (char*) malloc (lenImage*hImage) ;  /* Not messalloc since freed by XDestroyImage in graphDestroy */
  xim = XCreateImage (display, visual, 8, ZPixmap, 
		      0, data, 
		      (unsigned int) wImage, (unsigned int) hImage,
		      32, lenImage) ;
  if (!xim)
    { free (data) ;
      return 0 ;
    }

  pixelRubber (pixels, (UCHAR*)data, padImage,
	       wPixels, hPixels, lenPixels, wImage, hImage) ;

  assMultipleInsert (gSubDev->images, pixels, xim);
  
  return xim ;
}



/* I'll bet good money that 'char *' is used to point to image data all over the place */
/* when it should be 'unsigned char *' */ 
static XImage *rawImage (unsigned char *data, int w, int h, int len)
{
  Visual *visual ;
  XImage *xim = NULL;

  if (!(visual = getVisual()))
    return 0 ;
  
  if (!gSubDev->images)
    gSubDev->images = assCreate () ;

  if (bits_per_pixel <= 8)
    greyChanged = 0;

  if (assFind(gSubDev->images, data, &xim))
    {
      const void *vp = xim ;
      Array bucket = 0 ;
      int iBucket = 0 ;
      while (assFindNext (gSubDev->images, data, &vp, &bucket, &iBucket))
	/* CI if (xim->data == data && xim->width == w && xim->height == h) */
	{
	  xim = (void *)vp ;
	  if (xim->width == w && xim->height == h)
	    {
	      if (!greyChanged)
		return xim ;
	    } 
	  else
	    greyChanged = 0;
	}
    }

  if (xim == NULL)
    greyChanged = 0;
  
  if (len % 4) 
    messcrash ("len for PixelsRaw = %d is not a multiple of 4", w) ;

  if (bits_per_pixel > 8) {
    int i, tot = len * h;
    unsigned char *ptr = data;

    /* Hugh, seems we have to take care about the LSB-MSB brain damage...  */
    ptr = data;
    if (bits_per_pixel == 16) {
      unsigned short *ldata;
      unsigned char *lptr;

      if (greyChanged)
	{
	ldata = (unsigned short *) xim->data;
	lptr = (unsigned char *)xim->data;
	}
      else
	{
	ldata = (unsigned short *) malloc(tot << 1);
	lptr = (unsigned char *) malloc(tot << 1);
	}

      if (byte_order == 0)
	/* Little endian...  */
	for (i = 0; i < tot; i++) {
	  unsigned short tem = colors[*ptr++].pixel;
	  *lptr++ = tem;
	  *lptr++ = tem >> 8;
	}
      else
	/* Big endian...  */
	for (i = 0; i < tot; i++) {
	  unsigned short tem = colors[*ptr++].pixel;
	  *lptr++ = tem >> 8;
	  *lptr++ = tem;
	}
      if (!greyChanged)
	xim = XCreateImage (display, visual, visualDepth, ZPixmap, 
			    0, (char *) ldata, 
			    (unsigned int) w, (unsigned int) h,
			    32, len * 2) ;
    } else if (bits_per_pixel == 24) {
      unsigned long *ldata;
      unsigned char *lptr;

      if (greyChanged)
	{
	ldata = (unsigned long *) xim->data;
	lptr = (unsigned char *) xim->data;
	}
      else
	{
	ldata = (unsigned long *) malloc(tot << 2);
	lptr = (unsigned char *) malloc(tot << 2);
	}

      if (byte_order == 0)
	/* Little endian...  */
	for (i = 0; i < tot; i++) {
	  unsigned long tem = colors[*ptr++].pixel;
	  *lptr++ = tem;
	  *lptr++ = tem >> 8;
	  *lptr++ = tem >> 16;
	}
      else
	/* Big endian...  */
	for (i = 0; i < tot; i++) {
	  unsigned long tem = colors[*ptr++].pixel;
	  *lptr++ = tem >> 16;
	  *lptr++ = tem >> 8;
	  *lptr++ = tem;
	}
      if (!greyChanged)
	xim = XCreateImage (display, visual, visualDepth, ZPixmap, 
			    0, (char *) ldata, 
			    (unsigned int) w, (unsigned int) h,
			    32, len * 3) ;
    } else {
      /* Must be 32...  */
      unsigned long *ldata;
      unsigned char *lptr;

      if (greyChanged)
	{
	ldata = (unsigned long *) xim->data;
	lptr = (unsigned char *) xim->data;
	}
      else
	{
	ldata = (unsigned long *)malloc(tot << 2);
	lptr = (unsigned char *)malloc(tot << 2);
	}

      if (byte_order == 0)
	/* Little endian...  */
	for (i = 0; i < tot; i++) {
	  unsigned long tem = colors[*ptr++].pixel;
	  *lptr++ = tem;
	  *lptr++ = tem >> 8;
	  *lptr++ = tem >> 16;
	  *lptr++ = tem >> 24;
	}
      else
	/* Big endian...  */
	for (i = 0; i < tot; i++) {
	  unsigned long tem = colors[*ptr++].pixel;
	  *lptr++ = tem >> 24;
	  *lptr++ = tem >> 16;
	  *lptr++ = tem >> 8;
	  *lptr++ = tem;
	}
      if (!greyChanged)
	xim = XCreateImage (display, visual, visualDepth, ZPixmap, 
			    0, (char *) ldata, 
			    (unsigned int) w, (unsigned int) h,
			    32, len * 4) ;
    }
  } else {
    xim = XCreateImage (display, visual, 8, ZPixmap, 
			0, (char *)data, 
			(unsigned int) w, (unsigned int) h,
			32, len) ;
  }

  if (xim && !greyChanged)
    assMultipleInsert (gSubDev->images, data, xim) ;

  greyChanged = 0;

  return xim ;
}

/* Marat introduced this to centralize setting colours.  He sets a stipple
   for the foreground colour.  Is this really better?
*/

static void setXcolors(GC gc, int fcol, int bcol)
{
  if (isMono)
    {
      XSetBackground (DG, BlackPixelOfScreen(XtScreen(root_widget)));
      XSetForeground (DG, WhitePixelOfScreen(XtScreen(root_widget)));
      XSetStipple (display, gc, *stipple[fcol & 0x1f]) ;  /* mask possibly smaller than NUM_TRUECOLORS] */
      XSetFillStyle( display, gc, FillOpaqueStippled);
    } 
  else 
    {
      XSetForeground (DG, colorInd[fcol & 0xff]) ;
      XSetBackground (DG, colorInd[bcol & 0xff]) ;
    }
}

/******* now draw boxes **********/

static int psize, psize2 ;
static int clipx1, clipx2, clipy1, clipy2 ;

#define xclip(z)  if (z < clipx1) z = clipx1 ; \
		  else if (z > clipx2) z = clipx2
#define yclip(z)  if (z < clipy1) z = clipy1 ; \
		  else if (z > clipy2) z = clipy2
#define xsafe(z)  if (z < -30000) z = -30000 ; \
		  else if (z > 30000) z = 30000
#define ysafe(z)  if (z < -30000) z = -30000 ; \
		  else if (z > 30000) z = 30000

static void drawBox (Box box)
{
  float  t ;
  int    x1, x2, y1, y2, r, s, old, fcol ;
  char   *text ;
  int    action ;
#ifdef DEBUG_RECURSION
  static int recursionLevel = 0 ;
#endif
  if (!stackExists (gStk)) /* does happen when playing fast with many windows
			    * during expose events on uncomplete windows
			    * this is probably due to fancy X11 event queing 
			    */
    return ;
#ifdef DEBUG_RECURSION
  ++recursionLevel ;
#endif
  if (box->x1 > box->x2 || box->y1 > box->y2 ||
      (x1 = uToXabs(box->x1)) > clipx2 ||
      (y1 = uToYabs(box->y1)) > clipy2 ||
      (x2 = uToXabs(box->x2)) < clipx1 ||
      (y2 = uToYabs(box->y2)) < clipy1)
    { int nDeep = 1 ;
      stackCursor (gStk,box->mark) ; /* sets position to mark */
      while (!stackAtEnd (gStk))
	switch (action = stackNext (gStk,int))
	  {
	  case BOX_END:
	    if (!--nDeep)
	      {
#ifdef DEBUG_RECURSION
		fprintf (stderr,"  Level %d return non-draw box-end\n",
			 recursionLevel--) ;
#endif
		return ;                        /* exit point */
	      }
	    break ;
	  case BOX_START:
	    r = stackNext (gStk, int) ;
	    ++nDeep ;
	    break ;
	  case COLOR: case TEXT_FORMAT:
	    r = stackNext (gStk,int) ; 
	    break ;
	  case LINE_WIDTH: case TEXT_HEIGHT: case POINT_SIZE:
	    t = stackNext (gStk,float) ;
	    break ;
	  case LINE: case RECTANGLE: case FILL_RECTANGLE:
	    t = stackNext (gStk,float) ;
	    t = stackNext (gStk,float) ;
	    t = stackNext (gStk,float) ;
	    t = stackNext (gStk,float) ;
	    break ;
	  case PIXELS: case PIXELS_RAW:
	    t = stackNext (gStk,float) ;
	    t = stackNext (gStk,float) ;
	    text = stackNext (gStk,char*) ;
            r = stackNext (gStk, int) ;
            r = stackNext (gStk, int) ;
            r = stackNext (gStk, int) ;
	    if (action == PIXELS)
	      { t = stackNext (gStk,float) ;
		t = stackNext (gStk,float) ;
	      }
	    break ;
          case POLYGON : case LINE_SEGS:
            r = stackNext (gStk, int) ;
            while (r--)
	      { t = stackNext (gStk,float) ;
                t = stackNext (gStk,float) ;
	      }
            break ;
	  case CIRCLE: case POINT: 
	  case TEXT: case TEXT_UP: case TEXT_PTR: case TEXT_PTR_PTR: case TEXT_EXTERNAL: 
	  case COLOR_SQUARES: case FILL_ARC: case ARC:
	    t = stackNext (gStk,float) ;
	    t = stackNext (gStk,float) ;
	    switch (action)
	      {
	      case CIRCLE:
		t = stackNext (gStk,float) ;
		break ;
	      case FILL_ARC: case ARC:
		t = stackNext (gStk,float) ;
		t = stackNext (gStk,float) ;
		t = stackNext (gStk,float) ;
		break ;
	      case POINT:
		break ;
	      case TEXT: case TEXT_UP:
		text = stackNextText (gStk) ;
		break ;
	      case TEXT_PTR:
		text = stackNext (gStk,char*) ;
		break ;
          case TEXT_PTR_PTR:
		text = *stackNext (gStk,char**) ;
		break ;
          case TEXT_EXTERNAL:
		text = *stackNext (gStk,char**) ;
		break ;
	      case COLOR_SQUARES:
		text = stackNext (gStk,char*) ;
		r = stackNext (gStk, int) ;
		r = stackNext (gStk, int) ;
		text = (char*) stackNext (gStk,int*) ;
		break ;
	      }
	    break ;
	  case IMAGE:
	    (void)stackNext (gStk, void*) ; 
	    break ;
	  default:
	    messout ("Invalid draw action %d received",action) ;
	  }
#ifdef DEBUG_RECURSION
      fprintf (stderr, "  Level %d return bottom of stack nondraw\n", 
	       recursionLevel--) ;
#endif
      return ;
    }

  if (box->bcol != TRANSPARENT)
    { xclip(x1) ; xclip(x2) ; yclip(y1) ; yclip(y2) ;
      setXcolors(gc, box->bcol, box->bcol);
      XFillRectangle (DWG, x1, y1, (x2-x1)+1, (y2-y1)+1) ;
    }

  fcol = box->fcol;
  setXcolors(gc, fcol, box->bcol);
  
  stackCursor (gStk,box->mark) ; /* sets position to mark */

  while (!stackAtEnd (gStk))
    switch (action = stackNext (gStk,int))
      {
      case BOX_END:
#ifdef DEBUG_RECURSION
	fprintf (stderr,"  Level %d return\n", recursionLevel--) ;
#endif
        return ;                        /* exit point */
      case BOX_START:
        r = stackNext (gStk, int) ;
#ifdef DEBUG_RECURSION
	fprintf (stderr, "Level %d mark %d calls %d\n", 
		 recursionLevel, box->mark, r) ;
#endif
        drawBox (gBoxGet (r)) ;             /* recursion */
	setXcolors(gc, fcol, box->bcol);
        break ;
      case COLOR:
        x1 = stackNext (gStk,int) ; 
	if (x1 == FORECOLOR) 
	  fcol = box->fcol; 
	else if (x1 == BACKCOLOR)
	  fcol = box->bcol;
	else fcol = x1;
	setXcolors(gc, fcol, box->bcol);
        break ;
      case LINE_WIDTH:
        t = stackNext (gStk,float) ;
	gcValues.line_width = uToXrel(t) ;
	gcValuesMask = GCLineWidth; 
        XChangeGC (DG, gcValuesMask, &gcValues);
        break ;
      case TEXT_HEIGHT:
        t = stackNext (gStk,float) ; x1 = uToYrel(t) ;
        fontSet (x1) ;
        break ;
      case TEXT_FORMAT:
        x1 = stackNext (gStk,int) ; 
        fontSetFormat (x1) ;
        break ;
      case POINT_SIZE:
        t = stackNext (gStk,float) ;
        psize = uToXrel(t) ; psize2 = psize/2 ; 
	if (!psize2) psize2 = 1 ;
        break ;
      case LINE: case RECTANGLE: case FILL_RECTANGLE:
        t = stackNext (gStk,float) ; x1 = uToXabs(t) ; xsafe(x1) ;
        t = stackNext (gStk,float) ; y1 = uToYabs(t) ; ysafe(y1) ;
        t = stackNext (gStk,float) ; x2 = uToXabs(t) ; xsafe(x2) ;
        t = stackNext (gStk,float) ; y2 = uToYabs(t) ; ysafe(y2) ;
   	switch (action)
          {
          case LINE:
	    if (fcol == TRANSPARENT) break ; /* mieg 2007_09_30 */
            if (!((x1 < clipx1 && x2 < clipx1) ||
		  (x1 > clipx2 && x2 > clipx2) ||
		  (y1 < clipy1 && y2 < clipy1) ||	/* not perfect */
		  (y1 > clipy2 && y2 > clipy2) )) /* excludes both ends left,right,above or below */
	      XDrawLine (DWG, x1, y1, x2, y2);
            break ;
          case RECTANGLE:
	    if (x2 == x1) x2 = x1+1 ;
	    if (y2 == y1) y2 = y1+1 ;
	    if (x2 < x1) { r = x1 ; x1 = x2 ; x2 = r ; }
	    if (y2 < y1) { r = y1 ; y1 = y2 ; y2 = r ; }
            if (x1 <= clipx2 && y1 <= clipy2 && x2 >= clipx1 && y2 >= clipy1)
	      XDrawRectangle (DWG, x1, y1, (x2-x1), (y2-y1));
            break ;
          case FILL_RECTANGLE:
      	    if (x2 == x1) x2 = x1+1 ;
	    if (x2 < x1) { r = x1 ; x1 = x2 ; x2 = r ; }
	    if (y2 == y1) y2 = y1+1 ;
	    if (y2 < y1) { r = y1 ; y1 = y2 ; y2 = r ; }
            if (x1 <= clipx2 && y1 <= clipy2 && x2 >= clipx1 && y2 >= clipy1)
	      XFillRectangle (DWG, x1, y1, (x2-x1), (y2-y1));
            break ;
	  }
	break ;
      case PIXELS: case PIXELS_RAW:
	{ 
	  int xbase, ybase, w, h, len ;
	  char *pixels ;
	  XImage *xim ;

	  t = stackNext (gStk,float) ; x1 = uToXabs(t) ; xbase = x1 ;
	  t = stackNext (gStk,float) ; y1 = uToYabs(t) ; ybase = y1 ;
	  pixels = stackNext (gStk,char*) ;
	  w = stackNext (gStk, int) ;
	  h = stackNext (gStk, int) ;
	  len = stackNext (gStk, int) ;
	  if (action == PIXELS)
	    { t = stackNext (gStk,float) ; x2 = uToXabs(t) ;
	      t = stackNext (gStk,float) ; y2 = uToYabs(t) ;
		  /* suz, xsafe in place of xclip in this case */
	      xsafe(x1) ; xsafe(x2) ; ysafe(y1) ; ysafe(y2) ;
	    }
	  else
	    { x2 = x1 + w ;
	      y2 = y1 + h ;
	      xclip(x1) ; xclip(x2) ; yclip(y1) ; yclip(y2) ;
	    }
	  if (x1 < x2 && y1 < y2)
	    { if (action == PIXELS)
		xim = pixelsImage (pixels, w, h, len, 
				   (x2-x1), (y2-y1)) ;
	      else
		xim = rawImage ((unsigned char *)pixels, w, h, len) ; 
	      if (xim)
		XPutImage (DWG, xim, x1-xbase, y1-ybase, 
			   x1, y1, 
			   (x2 < x1+w) ? x2-x1+1 : x2-x1, 
			   (y2 < y1+h) ? y2-y1+1 : y2-y1) ;
	    }
	}
	break ;
      case POLYGON : case LINE_SEGS :
	{ int numVertices = stackNext (gStk, int) ;
	  XPoint *XPts = (XPoint *) messalloc (numVertices * 
					       sizeof (XPoint)) ;
	  XPoint *XPt = XPts ;

	  for (r = 0 ; r < numVertices ; r++)
	    { t = stackNext(gStk,float) ; 
	      x1 = uToXabs(t) ; xsafe(x1) ; XPt->x = x1 ;
	      t = stackNext(gStk,float) ; 
	      y1 = uToYabs(t) ; ysafe(y1) ; XPt->y = y1 ;
	      XPt++;
	    }
        /* NOTE this assumes the polygons are of the simplest
         * possible type, with no intersecting lines */
	  if (action == POLYGON)
	    XFillPolygon (DWG, XPts, numVertices, 
			  Convex, CoordModeOrigin) ;
	  else
	    XDrawLines (DWG, XPts, numVertices, CoordModeOrigin);
	  messfree (XPts) ;
	}
        break ;
      case CIRCLE: case POINT: 
      case TEXT: case TEXT_UP: case TEXT_PTR: case TEXT_PTR_PTR: case TEXT_EXTERNAL:
      case COLOR_SQUARES: case FILL_ARC: case ARC:
        t = stackNext (gStk,float) ; x1 = uToXabs(t) ;
	t = stackNext (gStk,float) ; y1 = uToYabs(t) ;
        switch (action)
          {
          case CIRCLE:      /* given center x1, y1, radius r*/
            t = stackNext (gStk,float) ; r = uToXrel(t) ;
	    if (!r) r = 1 ;
	    if (x1 - r < clipx2 && x1 + r > clipx1 &&
		y1 - r < clipy2 && y1 + r > clipy1)
	      XDrawArc (DWG, x1-r, y1-r, 2*r, 2*r, 0, 64*360);
		    /* NB upper left corner and 64ths of a degree */
            break ;
	  case ARC:
            t = stackNext (gStk,float) ; r = uToXrel(t) ;
	    { float ang = stackNext (gStk,float) ;
	      float angDiff = stackNext (gStk,float) ;
	      if (!r) r = 1 ;
	    if (x1 - r < clipx2 && x1 + r > clipx1 &&
		y1 - r < clipy2 && y1 + r > clipy1)
		XDrawArc (DWG, x1-r, y1-r, 2*r, 2*r, 
			  (int)(ang*64), (int)(angDiff*64)) ;
	      /* NB upper left corner and 64ths of a degree */
	    }
	    break ;
	  case FILL_ARC:
            t = stackNext (gStk,float) ; r = uToXrel(t) ;
	    { float ang = stackNext (gStk,float) ;
	      float angDiff = stackNext (gStk,float) ;
	      if (!r) r = 1 ;
	    if (x1 - r < clipx2 && x1 + r > clipx1 &&
		y1 - r < clipy2 && y1 + r > clipy1)
		XFillArc (DWG, x1-r, y1-r, 2*r, 2*r, 
			  (int)(ang*64), (int)(angDiff*64)) ;
	      /* NB upper left corner and 64ths of a degree */
	    }
	    break ;
          case POINT:
	    if (x1 - psize2 < clipx2 && x1 + psize2 > clipx1 &&
		y1 - psize2 < clipy2 && y1 + psize2 > clipy1)
	      XFillRectangle (DWG, x1-psize2, y1-psize2, psize, psize) ;
            break ;
          case TEXT:
            text = stackNextText (gStk) ;
	    s = strlen(text) ;
	    if (x1 < clipx2 && x1 + s*fontWidth > clipx1 &&
		y1 < clipy2 && y1 + fontHeight > clipy1)
	      XDrawString(DWG, x1, y1+fontdy, text, s) ; 
            break ;
	  case TEXT_UP:
            text = stackNextText (gStk) ;
	    s = strlen(text) ;
	    { int oldHeight = lastFontHeight ;
	      int i ;
	      XFontStruct *fnt = fontSet (oldHeight*0.6) ;
			/* next lines remove line-spacing allowance */
	      /* fontHeight = fnt->max_bounds.ascent + fnt->max_bounds.descent ; */
	      fontHeight = uToYabs(0.65);
	      
	      fontdy = fnt->max_bounds.ascent ;
	      if (x1 < clipx2 && x1 + fontWidth > clipx1 &&
		  y1 - s*fontHeight < clipy2 && y1 > clipy1)
		for (i = 0 ; i < s ; ++i)
		  XDrawString (DWG, x1, (y1-(s-i-1)*fontHeight)-uToYabs(0.2),
			       &text[i], 1) ;
	      fontSet (oldHeight) ;
	    }
            break ;
	  case TEXT_PTR:
	    text = stackNext (gStk,char*) ;
	    s = strlen(text) ;
	    if (x1 < clipx2 && x1 + s*fontWidth > clipx1 &&
		y1 < clipy2 && y1 + fontHeight > clipy1)
	      XDrawString(DWG, x1, y1+fontdy, text, s) ;
	    break ;
	  case TEXT_PTR_PTR:
	    if (!(text = *stackNext (gStk,char**)))
	      break ;
	    s = strlen(text) ;
	    if (x1 < clipx2 && x1 + s*fontWidth > clipx1 &&
		y1 < clipy2 && y1 + fontHeight > clipy1)
	      XDrawString(DWG, x1, y1+fontdy, text, s) ; 
	    break ;
	  case TEXT_EXTERNAL:
	    {
	      char **cpp = stackNext (gStk,char**) ;
	      s = *stackNext (gStk,int *);
	      if (cpp && !(text = *cpp))
		{
		  if (x1 < clipx2 && x1 + s*fontWidth > clipx1 &&
		      y1 < clipy2 && y1 + fontHeight > clipy1)
		    XDrawString(DWG, x1, y1+fontdy, text, s) ; 
		}
	    }
	    break ;
	  case COLOR_SQUARES:
	    { int *tints ;
	      text = stackNext (gStk,char*) ;
	      r = stackNext (gStk, int) ;
	      s = stackNext (gStk, int) ;
	      tints = stackNext (gStk, int*) ;
	      x2 = 0;
	      y2 = uToYrel (1.0) ;
	      if (x1 >= clipx2 || x1 + r*uToXrel(1.0) <= clipx1 ||
		  y1 >= clipy2 || y1 + y2 <= clipy1)
		break ;
	      old = fcol ;
	      while (r--)
		{ int newColor=0;
		  switch (*text^((*text)-1)) /* trick to find first on bit */
		    { 
		    case -1:   newColor = WHITE; break;
		    case 0x01: newColor = tints[0]; break;
		    case 0x03: newColor = tints[1]; break;
		    case 0x07: newColor = tints[2]; break;
		    case 0x0f: newColor = tints[3]; break;
		    case 0x1f: newColor = tints[4]; break;
		    case 0x3f: newColor = tints[5]; break;
		    case 0x7f: newColor = tints[6]; break;
		    case 0xff: newColor = tints[7]; break;
		    }
		  
		  if (newColor != old)
		    { if (x2 > 0)
			{ setXcolors(gc, old, box->bcol) ;
			  XFillRectangle (DWG, x1, y1, x2, y2) ;
			  x1 += x2; x2 = 0; 
			}
		      old = newColor ;
		    }
		  x2 += uToXrel(1.0) ;
		  text += s ;
		}
	      if (x2 > 0)
		{ setXcolors(gc, old, box->bcol) ;
		  XFillRectangle(DWG, x1, y1, x2, y2);
		} 
	      setXcolors(gc, fcol, box->bcol);
	      break ;
	    }
          }
	break ;
      case IMAGE:
	(void)stackNext (gStk, void*) ;
	break ;
      default:
	messout ("Invalid action %d received in drawBox()",action) ;
      }
#ifdef DEBUG_RECURSION
  fprintf (stderr, "  Level %d return at botom of func\n", 
	   recursionLevel--) ;
#endif
}

static void setBoxDefaults (Box box)
{
  if (!box) return ;

  gcValues.line_width = uToXrel (box->linewidth) ;
  gcValuesMask = GCLineWidth ; 
  XChangeGC (DG, gcValuesMask, &gcValues) ;

  fontSet (uToYrel(box->textheight)) ;
  fontSetFormat (PLAIN_FORMAT) ;
  psize = uToXrel(box->pointsize) ;
  psize2 = psize/2 ;
  if (!psize2)
    psize2 = 1 ;
}

void graphBoxDraw (int k, int fcol, int bcol)
{
  Box box = gBoxGet (k) ;

  if (fcol >= 0)
    box->fcol = fcol ;
  if (bcol >= 0)
    box->bcol = bcol ;

  if (!gDev || !gSubDev || !gSubDev->isExposed || gActive->isClear)
    return ;

  clipx1 = uToXabs(box->x1) ; if (clipx1 < 0) clipx1 = 0 ;
  clipx2 = uToXabs(box->x2) ; if (clipx2 > 30000) clipx2 = 30000 ;
  clipy1 = uToYabs(box->y1) ; if (clipy1 < 0) clipy1 = 0 ;
  clipy2 = uToYabs(box->y2) ; if (clipy2 > 30000) clipy2 = 0 ;
  if (clipx1 > clipx2 || clipy1 > clipy2)
    return ;

  setBoxDefaults (box) ;
  XSetClipMask (DG, None) ;
  if (isMono)
    XSetClipMask (display, gcStipple, None) ;
  /*  fprintf(stderr, "\n\ngraphBoxDraw(%d) calling drawBox\n", k) ; */
  drawBox (box) ;
  XFlush (display) ;
}


/* This routine is the main routine for filling windows with their contents. */
/* It can be called because a window is new, because its contents have       */
/* changed, because the window has been resized or exposed.                  */
/*                                                                           */
void graphClipDraw (int x1, int y1, int x2, int y2) /* expose action */
{
  static XRectangle clipRect[1] ;
  Box box = gBoxGet(0) ;

  /* We only do any drawing if there is an active window and subwindow, if   */
  /* there is anything on the stack to draw and the window has not been      */
  /* cleared of everything.                                                  */
  /* Note, as drawing can take a while, we turn on a busy cursor and then     */
  /* turn it off once finished.                                               */
  if (gDev && gSubDev && gActive->stack && !gActive->isClear)
    {
      gSubDev->isExposed = TRUE ;
      
      setBoxDefaults (box) ;

      clipx1 = (x1 > 0) ? x1 : 0 ; 
      clipx2 = (x2 < 30000) ? x2 : 30000 ; 
      clipy1 = (y1 > 0) ? y1 : 0 ; 
      clipy2 = (y2 < 30000) ? y2 : 30000 ;

      if (clipx1 <= clipx2 && clipy1 <= clipy2)
	{
	  graphBusyCursorAll (TRUE);
	  
	  clipRect->x = clipx1 ; clipRect->y = clipy1 ;
	  clipRect->width = clipx2 - clipx1 + 1 ; 
	  clipRect->height = clipy2 - clipy1 + 1 ;
	  XSetClipRectangles (DG, 0, 0, clipRect, 1, Unsorted) ;
	  if (isMono)
	    XSetClipRectangles (display, gcStipple, 
				0, 0, clipRect, 1, Unsorted) ;
	  /*  fprintf(stderr, "\n\ngraphClipDraw calling drawBox\n") ; */
	  drawBox (box) ;
	  
	  XFlush (display) ;

	  graphBusyCursorAll (FALSE) ;
	}
    }

  return ;
}


void graphRedraw (void)	/* now do this by creating expose event */
{
  unsigned short w, h ;	/* true width, height */
  Widget sbar ;
  float  pos ;
  XExposeEvent ev ;

  if (!gActive)
    return ;
  gActive->isClear = FALSE ;
  if (!gActive->stack || !gDev || !gSubDev || !gSubDev->isExposed)
    return ;

  ev.type = Expose ;
  ev.display = XtDisplay (gDev->popup) ;
  ev.window = XtWindow (gSubDev->base) ;
  ev.x = ev.y = 0 ;
  XtVaGetValues (gSubDev->viewport, XtNwidth, &w, XtNheight, &h, NULL) ;
  ev.width = w ; ev.height = h ;
  if ((sbar = XtNameToWidget (gSubDev->viewport,"horizontal")))
    { XtVaGetValues (sbar, XtNtopOfThumb, &pos, NULL) ;
      ev.x = pos*gActive->w ;
    }
  if ((sbar = XtNameToWidget (gSubDev->viewport,"vertical")))
    { XtVaGetValues (sbar, XtNtopOfThumb, &pos, NULL) ;
      ev.y = pos*gActive->h ;
    }
  XSendEvent (ev.display, ev.window, True, ExposureMask, (XEvent*)&ev) ;
  XFlush (display) ;
}

void graphWhiteOut (void)
{
  if (!gDev || !gSubDev)
    return ;

  if (gSubDev->boxMenus)		/* clear up menus */
    assClear (gSubDev->boxMenus) ;

  XClearWindow (display, window) ; /* RMD 28/5/92 */
  XFlush (display) ;
}

void graphXorLine (float x1, float y1, float x2, float y2)
{
  int ix1 = uToXabs(x1) ;
  int iy1 = uToYabs(y1) ;
  int ix2 = uToXabs(x2) ;
  int iy2 = uToYabs(y2) ;
  
  if (gDev && gSubDev)
    XDrawLine (display, window, gcInvert, ix1, iy1, ix2, iy2);
}

void graphXorBox(int k, float x, float y)
{
  Box box = gBoxGet (k) ;
  int x1 = uToXabs(x) ;
  int y1 = uToYabs(y) ;
  int x2 = x1 + uToXrel(box->x2 - box->x1) ;
  int y2 = y1 + uToYrel(box->y2 - box->y1) ;

  if (x2 == x1 || y2 == y1)
    graphXorLine (x, y, x + (x2-x1), y + (y2-y1)) ;
  else if (gDev && gSubDev)
    XDrawRectangle (display, window, gcInvert, 
		    x1, y1, (x2-x1), (y2-y1));
}

void graphBusyCursor (BOOL on)
{
  static Cursor cursor = None, LeftArrow = None, Watch = None ;

  if (cursor == None)
    { LeftArrow  = XCreateFontCursor(display, XC_top_left_arrow);
      Watch      = XCreateFontCursor(display, XC_watch);
      cursor = LeftArrow ;
    }    

  /* protection : if this was called by the alarm handler on exit
     of the program, the window won't exist anymore, which results
     in a core-dump, so better exit before that */
  if (!gActive->subdev)
    return ;

  if (on) cursor = Watch ;
  else cursor = LeftArrow ;
  
  XDefineCursor(display, XtWindow(gActive->subdev->viewport), 
		cursor) ;

  XFlush(XtDisplay(gActive->subdev->viewport)) ;

  return ;
}

/************* end of file ****************/
