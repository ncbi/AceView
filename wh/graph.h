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
/*  $Id: graph.h,v 1.12 2015/08/12 14:12:11 mieg Exp $ */

#ifndef DEF_GRAPH_H
#define DEF_GRAPH_H

#include "regular.h"		/* libutil header */

#include "colours.h"		/* for enum Colour{..} */

				/* library EXPORT/IMPORT symbols */
#if defined (WIN32)
#include "win32libspec.h"
#else
#define GRAPH_FUNC_DCL
#define GRAPH_VAR_DCL	extern
#define GRAPH_FUNC_DEF
#define GRAPH_VAR_DEF
#endif

typedef int Graph ;             /* really graph->id, hide  data structure */

UTIL_VAR_DCL KEY isGifDisplay ; /* used to unset useless buttons etc, now equal to the 'view' */
GRAPH_VAR_DCL int graphContrastLookup[];
#define graphContrast(x) (graphContrastLookup[(x)])

enum TextFormat { PLAIN_FORMAT, BOLD, ITALIC, GREEK, 
		   FIXED_WIDTH} ; /* The implementation is machine dependant */

#define BACK_COLOR WHITE
#define FORE_COLOR BLACK
	/* max 256 colours (box->fcol, bcol are unsigned char) */

/* GraphEvent note the mouse events are required to be continuous in
   this way. Can't touch the order.
   DESTROY must be last event and no values must be assigned to the
   enums so that they are contiguous.
*/
enum GraphEvent {LEFT_DOWN, LEFT_DRAG, LEFT_UP,
		 MIDDLE_DOWN, MIDDLE_DRAG, MIDDLE_UP,
		 RIGHT_DOWN, RIGHT_DRAG, RIGHT_UP,
		 PICK, KEYBOARD, IDLE, RESIZE, 
		 MESSAGE_DESTROY,
		 DESTROY
		} ;



typedef void (*GraphFunc)(void) ;
typedef void (*MouseFunc)(double x, double y) ;
typedef BOOL (*IdleFunc)(int box, double x, double y) ;
typedef void (*PickFunc)(int box, double x, double y) ;
typedef void (*EntryFunc)(char*) ;
typedef BOOL (*TimerFunc)(int timerIdentifier,void * param);

typedef BOOL (*ToggleFunc)(BOOL) ;
typedef void (*ColouredButtonFunc)(void *arg) ;
typedef int  (*GraphCompletionFunc)(char *text, int len) ; /* return # */

   /* If you change the enumeration of graphTypes, 
    * change also the function pickGetGraphTypes
    */
enum GraphType  { PLAIN, TEXT_SCROLL, TEXT_FIT, MAP_SCROLL, 
		  PIXEL_SCROLL, TEXT_FULL_SCROLL, PIXEL_FIT,
		  TEXT_FULL_EDIT,
		NUMGRAPHTYPES } ; /* always last */

 

GRAPH_VAR_DCL float    graphEventX, graphEventY ;
GRAPH_VAR_DCL FREEOPT  graphColors[] ;         /* defined in graphDevSubs.c */


/**********************************************************/
/************ functions used only by the Kernel ***********/
/**********************************************************/ 
GRAPH_FUNC_DCL void    graphInit (int *argcptr, char **argv) ;
GRAPH_FUNC_DCL void    graphFinish (void) ;
GRAPH_FUNC_DCL void    graphProcessEvents (void) ; /* process events in queue */
				/* NB should use graphLoop instead */



/**** the following blocking functions force a user reply ****/
/**** they should be accessed via messXXX()               ****/
GRAPH_FUNC_DCL void graphOut (const char *text) ;
GRAPH_FUNC_DCL void graphError (const char *text) ;
GRAPH_FUNC_DCL void graphEnd (const char *text) ;
GRAPH_FUNC_DCL BOOL graphQuery (const char *text) ;
GRAPH_FUNC_DCL BOOL graphPrompt (const char *prompt, char *dfault, char *fmt) ;


/* graphSelect can either show selections in a scrolled subwindow (default)  */
/* or in plain window, graphSetSelectDisplayMode() is used to set the mode.  */
typedef enum _GraphSelectDisplayMode {GRAPH_SELECT_SCROLLED,
				      GRAPH_SELECT_PLAIN} GraphSelectDisplayMode ;
GRAPH_FUNC_DCL void graphSetSelectDisplayMode(GraphSelectDisplayMode mode) ;
GRAPH_FUNC_DCL BOOL graphSelect (KEY *kpt, FREEOPT *options) ;





/**********************************************************/
/********** functions used by the applications ************/
/**********************************************************/
 
GRAPH_FUNC_DCL Graph   graphCreate (int type, 		/* of enum graphType */
		     const char *name, 	/* used on title bar */
		     float x, float y, float w, float h) ;
                /* x,y,w,h are in units of screen height = 1.0 */
GRAPH_FUNC_DCL BOOL    graphActivate (Graph graph) ;	/* returns graphExists() */
GRAPH_FUNC_DCL BOOL    graphActivateChild (void) ; /* retruns graphActivate gActive->children */
GRAPH_FUNC_DCL Graph   graphActive (void) ;            /* returns active graph */
GRAPH_FUNC_DCL BOOL    graphExists (Graph g) ;		/* is g a valid graph? */


/* Graphs can be of three types for purposes of user interactions:           */
/*       GRAPH_BLOCKABLE - normal window (default), input blocked when a     */
/*                         modal dialog is posted.                           */
/*   GRAPH_NON_BLOCKABLE - normal window but input is not blocked when a     */
/*                         modal dialog is posted.                           */
/*        GRAPH_BLOCKING - modal dialog window, cannot itself be blocked,    */
/*                         there can only be one of these at a time.         */
/*                                                                           */
typedef enum _GraphWindowBlockMode {GRAPH_BLOCKABLE, GRAPH_NON_BLOCKABLE,
				    GRAPH_BLOCKING} GraphWindowBlockMode ;
GRAPH_FUNC_DCL void graphSetBlockMode(GraphWindowBlockMode mode) ;

/************** timer function - added by Vahan ***************/
GRAPH_FUNC_DCL BOOL graphTimerSet(int msecs,int timerIdentifier,TimerFunc func,void * param);
/************** subgraphs ***************/

GRAPH_FUNC_DCL Graph graphSubCreate (int type, float ux, float uy, float uw, float uh) ;
GRAPH_FUNC_DCL Graph graphParent (void) ;	/* 0 if not a subgraph */
GRAPH_FUNC_DCL BOOL  graphContainedIn (Graph parent) ;

/***************************************************************/ 
/******** everything from here on acts on active graph *********/

GRAPH_FUNC_DCL void    graphDestroy (void) ;
GRAPH_FUNC_DCL int     graphLoop (BOOL isBlock) ; /* X event loop - now reentrant */
GRAPH_FUNC_DCL BOOL	graphLoopReturn (int retval) ; /* returns */
				/* old names for these functions */
#define graphStart(g)	graphLoop(FALSE)
#define graphBlock()	graphLoop(TRUE)
#define graphUnBlock(x)	graphLoopReturn(x)

GRAPH_FUNC_DCL AC_HANDLE graphHandle (void) ; /* gets freed when graph is freed */
GRAPH_FUNC_DCL AC_HANDLE graphClearHandle (void) ; /* freed when graph is cleared */

GRAPH_FUNC_DCL void    graphCleanUp (void) ;	/* kills all but the active graph */

GRAPH_FUNC_DCL BOOL    uGraphAssociate (void* in, void* out) ;
#define graphAssociate(xin,xout)   uGraphAssociate((void*)(xin),(void*)(xout))
GRAPH_FUNC_DCL BOOL    uGraphAssFind (void *in, void* *out) ;
#define graphAssFind(xin,pxout)	uGraphAssFind(xin,(void**)(pxout))
GRAPH_FUNC_DCL BOOL    graphAssRemove (void *in) ;

GRAPH_FUNC_DCL void    graphPlainBounds (float ux, float uy, float uw, float uh, float aspect) ;
GRAPH_FUNC_DCL void    graphMapBounds (float ux, float uy, float uw, float uh, float aspect) ;
GRAPH_FUNC_DCL void    graphTextBounds (int nx, int ny) ;
GRAPH_FUNC_DCL void    graphPixelBounds (int nx, int ny) ;
GRAPH_FUNC_DCL float   graphFakeBounds (float ny) ; /* in lines, to prepare a long print on a virtual window */
GRAPH_FUNC_DCL void    graphFitBounds (int *nx, int *ny) ;  /* no scroll - in text coords */
GRAPH_FUNC_DCL void    graphTextInfo (int *dx, int *dy, float *cw, float *ch) ;
GRAPH_FUNC_DCL void    graphScreenSize (float *sx, float *sy, float *fx, float *fy, int *px, int *py) ;
		/* screen size in graphCreate units, TEXT units, pixels */
GRAPH_FUNC_DCL BOOL    graphWindowSize (float *xp, float *yp, float *wp, float *hp) ; 
                /* window position and size in graphCreate units */
GRAPH_FUNC_DCL void    graphGoto (float x, float y) ;  /* tries to center (x,y) in screen */
GRAPH_FUNC_DCL void    graphWhere (float *x1, float *y1, float *x2, float *y2) ;
		/* coords of visible section (when scrollbars) */
GRAPH_FUNC_DCL void    graphRedraw (void) ;            /* complete redraw */
GRAPH_FUNC_DCL void    graphClear (void) ;
GRAPH_FUNC_DCL void    graphRetitle (const char *name) ;
GRAPH_FUNC_DCL void    graphPop (void) ;

GRAPH_FUNC_DCL void    graphMessage (char *text) ;     /* Non blocking message */
GRAPH_FUNC_DCL void    graphUnMessage (void) ;         /* Applications generally use displayBlock */

GRAPH_FUNC_DCL void graphPS(char *myfilname, char *mail, char *print, 
			    char *title, BOOL doColor, BOOL isRotated,
			    float scaleFactor, int numPages) ; /* Write a PostScript image to a file. */
void graphPSdefaults (BOOL *isRotated, float *scaleFactor, int *numPages) ;
							    /* Get defaults for PostScript image. */
GRAPH_FUNC_DCL void    graphPrint(void);	/* JTM, to simplify the application codes */
GRAPH_FUNC_DCL void    graphBoundsPrint (float uw, float uh, void (*draw)(void)) ;
				/* to cheat the graph size for multipage prints */
GRAPH_FUNC_DCL void    graphEvent (int action, float x, float y) ;     /* posts an event */

GRAPH_FUNC_DCL int     graphLastBox(void); 
GRAPH_FUNC_DCL int     graphBoxStart(void) ;           /* returns box number */
GRAPH_FUNC_DCL void    graphBoxEnd (void) ;
GRAPH_FUNC_DCL BOOL    graphBoxClear (int box) ;
GRAPH_FUNC_DCL BOOL    graphBoxMarkAsClear (int box) ;
GRAPH_FUNC_DCL void    graphBoxDim (int box, float *x1, float *y1, float *x2, float *y2) ;
GRAPH_FUNC_DCL int     graphBoxAt (float x, float y, float *rx, float *ry);
GRAPH_FUNC_DCL void    graphBoxDraw (int box, int fcol, int bcol) ;
GRAPH_FUNC_DCL void    graphClipDraw (int x1, int y1, int x2, int y2);
				      /* expose action */

GRAPH_FUNC_DCL void    graphBoxDrag (int box, void (*clientRoutine)(float*,float*,BOOL)) ;
GRAPH_FUNC_DCL void    graphBoxShift (int box, float xbase, float ybase) ;
GRAPH_FUNC_DCL void    graphBoxSetPick (int box, BOOL pick) ;
GRAPH_FUNC_DCL void    graphBoxSetMenu (int k, BOOL menu) ;

GRAPH_FUNC_DCL void    graphLine (float x0, float y0, float x1, float y1) ;
GRAPH_FUNC_DCL void    graphRectangle (float x0, float y0, float x1, float y1) ;
GRAPH_FUNC_DCL void    graphFillRectangle (float x0, float y0, float x1, float y1) ;
GRAPH_FUNC_DCL void    graphPolygon (Array pts) ;	/* filled polygon */
GRAPH_FUNC_DCL void    graphLineSegs (Array pts) ;	/* line segment */
GRAPH_FUNC_DCL void    graphCircle (float xcen, float ycen, float rad) ;
GRAPH_FUNC_DCL void    graphArc (float xcen, float yxen, float rad,
		  float ang, float angDiff) ; 
		/* degrees anticlockwise from 3 o'clock */
GRAPH_FUNC_DCL void    graphFillArc (float xcen, float yxen, float rad,
		      float ang, float angDiff) ; 
		/* degrees anticlockwise from 3 o'clock */
GRAPH_FUNC_DCL void    graphPoint (float x, float y) ;
GRAPH_FUNC_DCL void    graphText (const char *text, float x0, float y0) ;
GRAPH_FUNC_DCL void    graphTextUp (const char *text, float x0, float y0) ;
     	   /* updatable text - length is length to extend the box */
GRAPH_FUNC_DCL void    graphTextPtr (const char *text, float x0, float y0, int length) ; 
GRAPH_FUNC_DCL void    graphTextPtrPtr (const char **text, float x0, float y0, int length) ;
GRAPH_FUNC_DCL void    graphTextExternal (const char **text, float x0, float y0, int * length) ;
GRAPH_FUNC_DCL void    graphColorSquares (char *cols, float x0, float y0, int length, 
			   int skip, int *tints) ;

GRAPH_FUNC_DCL void	graphPixels (char *pixels, int w, int h, int lineWidth,
		     float x0, float y0, float x1, float y1) ;
GRAPH_FUNC_DCL void	graphPixelsRaw (char *pixels, int w, int h, int lineWidth,
		        float x0, float y0) ;
GRAPH_FUNC_DCL BOOL    graphRawMaps (unsigned char *forward, int *reverse) ;
   /* forward[] maps [0,255] -> PixelsRaw values, 
      reverse[] inverts, giving max value so can test for 255 
      return FALSE if pixel images do not exist */
GRAPH_FUNC_DCL void    graphColorMaps (float *red, float *green, float *blue) ;
   /* colormap values are real in [0,1], and map from [0,255] */

GRAPH_FUNC_DCL BOOL    graphSetColorMaps (unsigned char *red,
			   unsigned char *green,
			   unsigned char *blue) ;
   /* set colors with the given maps */
   /* arrays are 128 bytes long      */
GRAPH_VAR_DCL  Graph rampGraph ;
GRAPH_FUNC_DCL void    graphGreyRamp (int min, int max) ;
GRAPH_FUNC_DCL void    graphRampTool (void) ;
GRAPH_FUNC_DCL BOOL    graphWritableColors (void) ;
 
   /* the next routines does not add to the box : used for cursors */
GRAPH_FUNC_DCL void    graphXorLine (float x0, float y0, float x1, float y1) ;
GRAPH_FUNC_DCL void    graphXorBox (int kbox, float x, float y) ;

   /* The following routines all return the old value, and just return
      the current value if given a negative argument (except 
      graphRegister). */
GRAPH_FUNC_DCL float   graphPointsize (float x) ;
GRAPH_FUNC_DCL float   graphLinewidth (float x) ;
GRAPH_FUNC_DCL float   graphTextHeight (float height) ; /* if 0 then standard text */
GRAPH_FUNC_DCL int     graphTextFormat(int format) ; 

#if defined(WIN32)
/* Specifies the location of the alignment origin for x,y coordinates of text.
   LEFT, CENTRE or RIGHT, and TOP or BOTTOM are permitted alignment values.
   Note: Logical OR ("|") can be used to combine flags:
   LEFT, CENTRE or RIGHT can be combined with either TOP or BOTTOM
   The default alignment is "TOP|LEFT" */
enum TextAlignment { LEFT = 0, CENTRE = 1, RIGHT = 2, TOP = 4, BOTTOM = 8 } ; 
GRAPH_FUNC_DCL int graphTextAlign (int textAlignment) ;
#endif

GRAPH_FUNC_DCL int     graphColor (int color) ;         /* color is a GraphColor */
GRAPH_FUNC_DCL void    graphBoxColor (int k, int fcol, int bcol) ; /* col==-1: no change */

#define  graphRegister(k,f) uGraphRegister(k, (GraphFunc)(f))
GRAPH_FUNC_DCL GraphFunc uGraphRegister (int k, GraphFunc) ;   /* k is a GraphEvent */
GRAPH_FUNC_DCL char*     graphHelp (char* item) ;

GRAPH_FUNC_DCL BOOL graphHtmlViewer (const char* helpFilename);
GRAPH_FUNC_DCL BOOL graphWebBrowser (const char* link); /* remote-controlled netscape via X-atoms */

typedef struct colouropt
{ void (*f)(void *arg);
  void *arg;
  int fg;
  int bg;
  const char *text;
  struct colouropt *next;
} COLOUROPT;

typedef struct menuopt
  { void (*f)(void);
    const char *text;
  } MENUOPT ;

/* typecast for MenuFunction's => struct menuopt x->f VoidRoutine's */
#define MENU_FUNC_CAST (void(*)(void))

typedef void (*MenuRoutine)(KEY) ;
typedef void (*MenuFunction)(int) ;
typedef void (*FreeMenuFunction)(KEY, int) ;

GRAPH_VAR_DCL int menuBox ;
GRAPH_FUNC_DCL void    graphBoxFreeMenu (int box, FreeMenuFunction proc, FREEOPT *options) ;
#define   graphFreeMenu(x,y)  graphBoxFreeMenu(0,(FreeMenuFunction)(x),(y))
GRAPH_FUNC_DCL void	graphBoxMenu (int box, MENUOPT *menu) ;
#define   graphMenu(x)  graphBoxMenu(0,(x))
GRAPH_FUNC_DCL int  graphMenuTriangle (BOOL filled, float x, float y) ;
	/* graphics routine for drawing a pulldown symbol */

GRAPH_FUNC_DCL void    graphBubbleInfo (int box, KEY key, char *ficheName, char *bubbleName) ; /* offset in buble stack */


GRAPH_FUNC_DCL void    graphBoxInfo (int box, KEY key, char *text) ;
GRAPH_FUNC_DCL void    graphBoxInfoFile (FILE *fil) ;
GRAPH_FUNC_DCL void    graphBoxInfoWrite (void) ;
 
GRAPH_FUNC_DCL int     graphButton (const char* text, VoidRoutine func, float x, float y) ;
GRAPH_FUNC_DCL int     graphColouredButton (const char *text,
			     ColouredButtonFunc func,
			     void *arg,
			     int fg,
			     int bg,
			     float x,
			     float y);
GRAPH_FUNC_DCL int     graphButtons (MENUOPT *buttons, float x0, float y0, float xmax) ;
GRAPH_FUNC_DCL int     graphColouredButtons (COLOUROPT *buttons, float x0, float *y0p, float xmax) ;

/**************** new menus ***********/

#ifndef MENU_DEFINED
#define MENU_DEFINED
typedef void *MENU, *MENUITEM ;	/* public handles */
typedef void (*MENUFUNCTION)(MENUITEM) ;
#endif

GRAPH_FUNC_DCL void	graphNewBoxMenu (int box, MENU menu) ;
#define graphNewMenu(x) graphNewBoxMenu(0,(x))



#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
/* These routines appear not to be coded or referenced anywhere....          */
/*                                                                           */
GRAPH_FUNC_DCL int	graphNewButton (MENUITEM item, float x, float y) ;
GRAPH_FUNC_DCL int	graphNewButtons (MENU menu, float x0, float y0, float xmax) ;

/* graphics routine for drawing a pulldown symbol - moved from menu.c/h */
GRAPH_FUNC_DCL int  menuDownSelTriangle (BOOL filled, float x, float y) ;

#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */


/************** text entry ***************/

typedef void (*GraphEntryFunc)(char*) ;

GRAPH_FUNC_DCL int   graphTextEntry (char* cp, int len, float x, float y, GraphEntryFunc fn) ;
GRAPH_FUNC_DCL int   graphTextScrollEntry ( char* text, int len, int wlen, 
	 				  float x, float y, GraphEntryFunc fn) ;
GRAPH_FUNC_DCL int   graphCompletionEntry (GraphCompletionFunc f, char* text, int len, 
		 			 float x, float y, GraphEntryFunc fn) ;
GRAPH_FUNC_DCL int   graphCompScrollEntry (GraphCompletionFunc f, char* text, int len, int wlen, 
					 float x, float y, GraphEntryFunc fn) ;
GRAPH_FUNC_DCL void  graphEntryDisable (void) ; /* disables TextEntries */

#define TEXT_ENTRY_FUNC_CAST (void (*)(char*))

/******** interaction with rest of workspace *******/

GRAPH_FUNC_DCL BOOL  graphInterruptCalled(void) ;
GRAPH_FUNC_DCL void  graphPostBuffer (char *text) ; /* write to screen cut/paste buffer */
GRAPH_FUNC_DCL char* graphPasteBuffer (void) ;	/* get screen cut/paste buffer contents */

#if defined(WIN32)
GRAPH_FUNC_DCL void  graphCopy(void) ;		  /* User initiated copy */
GRAPH_FUNC_DCL void  graphPaste(void) ;		  /* User initiated paste */
#endif

/******** editors - single item labelled textEntry boxes *********/
GRAPH_FUNC_DCL void editorToggleChange(int box);
GRAPH_FUNC_DCL void graphToggleEditor (char *label, BOOL *p, float x, float y);
GRAPH_FUNC_DCL void graphAddRadioEditor(char *text, int tag, int index, float x, float y);
GRAPH_FUNC_DCL void graphSetRadioEditor(int tag,int index);
GRAPH_FUNC_DCL int  graphRadioCreate(char *text,int *p,float x, float y);
GRAPH_FUNC_DCL void graphColourEditor(char *label, char *text, int *p,float x, float y);


GRAPH_FUNC_DCL int graphIntEditor (char *label, int *p, float x, float y, 
		     BOOL (*checkFunc)(int)) ;
GRAPH_FUNC_DCL int graphFloatEditor (char *label, float *p, float x, float y,
		     BOOL (*checkFunc)(float)) ;
GRAPH_FUNC_DCL int graphWordEditor (char *label, char *text, int len, float x, float y,
		    BOOL (*checkFunc)(char*)) ;
GRAPH_FUNC_DCL int graphTextEditor (char *label, char *text, int len, float x, float y,
		     BOOL (*checkFunc)(char*)) ;
GRAPH_FUNC_DCL int graphTextScrollEditor (char *label, char *text, int len, int wlen, float x, float y,
		    BOOL (*checkFunc)(char* text, int box)) ;
GRAPH_FUNC_DCL int graphMultiLineEditor (char *label, char *text, int len, int width, int height, 
					 float x, float y, 
					 BOOL (*checkFunc)(char* text, int box)) ;
GRAPH_FUNC_DCL BOOL graphCheckEditors (Graph graph, BOOL callCheckFunctions) ;
GRAPH_FUNC_DCL BOOL graphUpdateEditors (Graph graph) ;
void graphWhiteOut (void) ;	/* used to clear screen */

/********* FLASH stuff, using swfc tool box  *********/
BOOL swfGraphExport (Graph gId, ACEOUT out, BOOL do_size) ;
BOOL swfGraphResize (float width, float height) ; /* resize a flash graph (in graphFitBounds units) */

/********* GIF stuff, using Tim Boutell's gd package (c) Cold Spring Harbor *********/
GRAPH_FUNC_DCL BOOL graphGIFfile (FILE *fil) ;
GRAPH_FUNC_DCL BOOL graphGIFname (char *name, BOOL dumpBox) ;
int graphGIF (Graph gId, ACEOUT out, BOOL do_size) ;
GRAPH_FUNC_DCL char* graphGIFPtr (int* size) ;
GRAPH_FUNC_DCL int  graphGIFread (FILE *fil, float x0, float y0) ;
GRAPH_FUNC_DCL void gifLeftDown (int x, int y) ;
BOOL graphDumpBubble (char *fileName) ;
KEY gMapSelectView (KEY v) ;


/********* help function prototypes  ******/

/* this void-void function is used in menustructures or for graphButtons
   to bring up the help-topic currently registered for this window */
/* all other help-system stuff is in help.h (not graph-specific) */
GRAPH_FUNC_DCL void  help (void) ;

/********* Some WIN32 specific stuff ******/
#if defined(WIN32)
/*
 * This function is a user interface item so I put it in the 
 * graphics library...
 * Generates a system default standard beep sound 
 */
GRAPH_FUNC_DCL void graphSysBeep(void) ;

#endif /* defined(WIN32) */


#endif /*  !defined DEF_GRAPH_H */

 
 
