/*  File: graph_.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: version of graph.h for inclusion within graphics package
 *              contains device independent private information and then 
 *              includes graph.h (public information).
 *              use graphx_.h, graphsun_.h etc to get at device dependent stuff
 *
 * Exported functions:
 * HISTORY:
 * Last edited: May 12 14:00 1999 (edgrif)
 * * May 12 13:47 1999 (edgrif): Added state to graphstruct for cursor setting,
 *              plus some strings for error messages.
 * * Apr 16 13:34 1999 (edgrif): Move graphBoxSetMenu to graph.h
 * * Mar 24 14:36 1999 (edgrif): Fix SANgc03641, see NONMENU flag.
 * * Jan  5 13:48 1999 (edgrif): Move graphPS to public header.
 * * Dec 16 15:10 1998 (edgrif): Moved waitCursor calls to internal busyCursor calls.
 * * Oct 22 14:19 1998 (edgrif): Added message interface for use in providing
 *              device dependent routines: graphWinOut etc.
 * * Oct 13 14:19 1998 (edgrif): Added func decs for graph/acedb interface.
 *              Mostly these are access functions to get 'callback' type
 *              functions. Plus added this standard format header.
 * Created: Tue Oct 13 14:17:04 1998 (edgrif)
 *-------------------------------------------------------------------
 */

/* $Id: graph_.h,v 1.6 2015/08/12 14:12:11 mieg Exp $ */

#ifndef DEF_GRAPH__H
#define DEF_GRAPH__H

#include "regular.h"
#include "dict.h"
#include "graphAcedbInterface.h" /* includes graph.h */
#include "w2/graphcolour.h"

#define GRAPH_MAGIC     17864

#ifdef DEV_DEFINED
typedef struct DevStruct *Dev ;
typedef struct SubDevStruct *SubDev;
#define DEV_SIZE        sizeof(struct DevStruct)
#define SUBDEV_SIZE     sizeof(struct SubDevStruct)
#else
typedef void *Dev ;
typedef void *SubDev ;
#endif


/* See GraphEvent enums in graph.h, this is the number of enums defined      */
/* there, in GraphEvent the last event is always DESTROY.                    */
enum {NEVENTS = DESTROY + 1} ;




typedef struct GraphStruct
  { int           magic ;         /* check for is a graph */
    char          *name ;
    int           id ;            /* unique graph number for process */
    int           type ;          /* from enum GraphType */

    BOOL          window_being_constructed ;		    /* Used for busy cursor setting. */
    GraphWindowBlockMode
                  block_mode ;				    /* Used for busy cursor currently. */


    float         ux,uy,uw,uh ;   /* user boundaries */
    float         aspect ;        /* x units per inch / y units per inch */
    int           xWin,yWin ;     /* userToWin offsets */
    float         xFac,yFac ;     /*  and scale factors */
    int           textX,textY ;   /* text size in pixels */
    int           h,w ;           /* device coords (0,0) to (w-1,h-1) */
    void          (*func[NEVENTS])() ; /* Used in graphxt with either 0 or 2 args */
    float         pointsize ;
    float         linewidth ;
    float         textheight ;
    int           color ;
    int           textFormat ;
    Array         boxes ;
    Stack         stack ;
    int           nbox ;
    int           currbox ;
    BOOL          isClear ; 

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
    BOOL          isBlocked ; /* used in mapPrint() */
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */


    Stack         boxstack ;
    Associator    assoc ;
    char*         help ;
    Associator    buttonAss ;
    Array	  boxInfo ;
    Array	  editors ;
    AC_HANDLE  handle ;	/* freed on graphDestroy */
    AC_HANDLE  clearHandle ;	/* freed every time graphClear() called */

				/* windows system specific stuff */
    SubDev        subdev ;	/* private subdev contains viewport */
    Dev           dev ;		/* all graphs point to the Dev of their window */

    struct GraphStruct *children ;       /* linked list of all subgraphs  */
    struct GraphStruct *parent ;         /* graph in which I am contained */
    struct GraphStruct *nxt ;            /* next graph  */

    Array bubbleInfo ;
    DICT *bubbleDict ;
  } *Graph_ ;

#define GRAPH_SIZE      sizeof(struct GraphStruct)

typedef struct
        { float     x1, y1, x2, y2 ;
          int       mark ;
          float     linewidth ;
          float     textheight ;
          float     pointsize ;
          unsigned char   fcol, bcol ;
	  unsigned char   flag, format;
#if defined(BETTER_BOX_SHIFT)	/* from rbrusk */
	  int	    id, parentid ;
#endif
        } BoxStruct, *Box ;

typedef struct
	{ KEY key ;
	  char *text ;
	} BOXINFO ;


typedef struct
	{ int box ;
          KEY key ;
	  int fName ;
	  int bName ;
	} BUBBLEINFO ;

/********* Dec 19 16:00 1995 (rbrusk): the GraphAction POINT enum member
 * re#define'd to acePOINT to avoid conflict with a Win32 API typedef;
 * Also need a typedef workaround for use of the 
 * WIN32 POINT types in our ACEDB code 
 ************************************************************/
#if defined(WIN32)
#ifdef _WINDEF_
typedef POINT 	WIN_PT ;
typedef LPPOINT LP_WIN_PT ;
static int WIN_PT_SIZE = sizeof (POINT) ;
#endif
#define POINT acePOINT
#endif

enum GraphAction
        { LINE, RECTANGLE, FILL_RECTANGLE, CIRCLE, POINT, 
	  TEXT, TEXT_UP, TEXT_PTR, TEXT_PTR_PTR, TEXT_EXTERNAL,
          COLOR, TEXT_FORMAT, TEXT_HEIGHT, LINE_WIDTH, POINT_SIZE,
          BOX_START, BOX_END, PIXELS, PIXELS_RAW, COLOR_SQUARES, 
	  POLYGON, FILL_ARC, IMAGE, ARC, LINE_SEGS
        } ;

/* NOPICK & NOMENU exclude the box from being included in the search for     */
/* which box is to receive button clicks.                                    */
#define GRAPH_BOX_NOPICK_FLAG	0x01
#define GRAPH_BOX_MENU_FLAG	0x02
#define GRAPH_BOX_BUTTON_FLAG	0x04
#define GRAPH_BOX_INFO_FLAG	0x08
#define GRAPH_BOX_ENTRY_FLAG	0x10
#define GRAPH_BOX_TOGGLE_FLAG	0x20
#define GRAPH_BOX_NOMENU_FLAG	0x40

/* Some useful defines for messages from Graph.                              */
#define GRAPH_INTERNAL_LOGIC_ERROR_MSG "Graph: Internal logic error, \
please report to acedb development, error was - "
#define GRAPH_RESTART_ERROR_MSG "You should save your work and restart \
the application as soon as possible."



/* externals only for use within graph package: gXxxx */
extern Graph_   gActive ;       /* active graph */
extern Dev      gDev ;          /* gActive->dev */
extern SubDev   gSubDev ;       /* gActive->subdev */
extern Box      gBox ;          /* active box from gActive */
extern Stack    gStk ;          /* gActive->stack */
extern float    gPixAspect ;    /* device x:y pixel ratio */
extern BOOL     gWriteColors ;  /* writable colours */
extern int      gIsCreating ;   /* for Activate control during Creation */

extern int rampMin, rampMax;	/* current greyramp values */

#define uToXrel(x) (int)(gActive->xFac*(x))
#define uToYrel(x) (int)(gActive->yFac*(x))
#define uToXabs(x) (-gActive->xWin + uToXrel(x))
#define uToYabs(x) (-gActive->yWin + uToYrel(x))
#define XtoUrel(x) ((x)/gActive->xFac)
#define YtoUrel(x) ((x)/gActive->yFac)
#define XtoUabs(x) XtoUrel((x)+gActive->xWin)
#define YtoUabs(x) YtoUrel((x)+gActive->yWin)

#define UtextX (gActive->textX/gActive->xFac)
#define UtextY (gActive->textY/gActive->yFac)

void graphFacMake (void) ;
void graphDevInit (int *argcptr, char **argv) ;
void graphDevFinish (void) ;
void graphDevStart (Graph_ g) ;          /* starts notifier loop system */
void graphDevActivate (int flag) ;
void graphDevCreate (Graph_ g, float x, float y, float w, float h) ;
void graphDevDestroy (void) ;
				/* Subdev's -- maratb */
void graphSubDevCreate (Graph_ graph, float x, float y, float w, float h);
void graphSubDevDestroy (void);

void graphDeleteContents (Graph_ theGraph); /* kills boxes/stacks etc. */

void devMenuDestroy (void) ;

void graphASCII (char *myfilname, char *mail, char *print, char *title) ;

Box  gBoxGet (int k) ;             /* gets Box pointer for index k */
void gLeftDown (float x, float y) ;
void gMiddleDown (float x, float y) ;
BOOL gIdle (float x, float y) ;
void gBoxClear (Box box) ;
BOOL gFontInfo (int height, int* w, int* h) ;
void gUpdateBox0 (void) ;


/* These have been hidden as to use them requires some understanding of how  */
/* the graph package works internally + we now have a mechanism for showing  */
/* busy cursors.                                                             */
void graphBusyCursor (BOOL) ;
void graphBusyCursorAll (BOOL) ;
void graphSetBusyCursorInactive(BOOL inactive) ;


void graphDevTextFacSet (void) ;   /* sets fac according to text */

/********* in filquery.c, to be registered in graphInit() **********/

FILE *graphQueryOpen (char *dname, char *fname, char *end, 
		      const char *spec, const char *title) ;


/********************  Graph message facility        *************************/
/* These routines will one day be the device dependant routines that are     */
/* called by the abstract public graph routines.                             */
/*                                                                           */
/* Text displayed on button for graph blocking messages.                     */
#define GRAPH_OUT_MSG     "Continue"
#define GRAPH_ERROR_MSG   "Continue"
#define GRAPH_END_MSG     "Exit"

/* The device dependent routines.                                            */
void graphWinOut (const char *text, char *label) ;
BOOL graphWinQuery(const char *text) ;
BOOL graphWinPrompt (const char *prompt, char *dfault, char *fmt) ;



/********************  Graph  <-->  Acedb interface  *************************/
/* Do not use this interface for anything other than setting up graph for    */
/* use with Acedb. It is here to enable the two to cooperate, not for others */
/* to make use of...you have been warned.                                    */
/*                                                                           */
GRAPH_FUNC_DEF void initOverloadFuncs() ;

GRAPH_FUNC_DEF VoidCharRoutine  getGraphAcedbMainActivity(void) ;

GRAPH_FUNC_DEF GraphCharRoutine getGraphAcedbDisplayCreate(void) ;
GRAPH_FUNC_DEF char *getGraphAcedbChronoName(void) ;
GRAPH_FUNC_DEF char *getGraphAcedbViewName(void) ;
GRAPH_FUNC_DEF char *getGraphAcedbFilqueryName(void) ;

GRAPH_FUNC_DEF MENUOPT *getGraphAcedbPDMenu(void) ;
GRAPH_FUNC_DEF CharVoidRoutine getGraphAcedbStyle(void) ;
GRAPH_FUNC_DEF char *getGraphAcedbLogin(void) ;

GRAPH_FUNC_DEF char *getGraphAcedbSessionName(void) ;
GRAPH_FUNC_DEF VoidFILEStackIntRoutine getGraphAcedbGetFonts(void) ;

GRAPH_FUNC_DEF VoidFILEKEYRoutine getGraphAcedbClassPrint(void) ;

GRAPH_FUNC_DEF VoidIntRoutine getGraphAcedbXFont(void) ;

GRAPH_FUNC_DEF VoidDisplayWindowRoutine getGraphAcedbWcommInit(void) ;
GRAPH_FUNC_DEF VoidVoidRoutine getGraphAcedbWcommFinish(void) ;
GRAPH_FUNC_DEF BoolXeventRoutine getGraphAcedbWprop(void) ;


GRAPH_FUNC_DEF VoidVIEWCOLCONTROLRoutine getGraphAcedbViewAnon(void) ;
GRAPH_FUNC_DEF BOOLMENURoutine getGraphAcedbSaveView(void) ;
GRAPH_FUNC_DEF VoidMAPCONTROLVIEWRoutine getGraphAcedbSaveViewMsg(void) ;
GRAPH_FUNC_DEF VoidMENUSPECRoutine getGraphAcedbViewMenu(void) ;

GRAPH_FUNC_DEF BOOLMAPCONTROLKEYRoutine getGraphAcedbGetSetView(void) ;
GRAPH_FUNC_DEF VoidKEYRoutine getGraphAcedbResetKey(void) ;
GRAPH_FUNC_DEF VoidCOLOBJRoutine getGraphAcedbResetSave(void) ;
GRAPH_FUNC_DEF VoidIntFREEOPT getGraphAcedbFreemenu(void) ;
GRAPH_FUNC_DEF VoidOBJFUNCMAPCONTROL getGraphAcedbMapLocate(void) ;
GRAPH_FUNC_DEF VoidOBJFUNCMAPCONTROL getGraphAcedbScale(void) ;
GRAPH_FUNC_DEF VoidOBJFUNCSPACERPRIV getGraphAcedbSpace(void) ;
GRAPH_FUNC_DEF VoidCharKEY getGraphAcedbAddMap(void) ;

/*****************************************************************************/



#endif /* DEF_GRAPH__H */

/***** end of file ********/
 
