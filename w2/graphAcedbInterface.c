/*  File: graphAcedbInterface.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: The graph package allows the registering of functions to
 *              'overload' some of it's functions. This is necessary for 
 *              ACEDB applications so that they can load users preferences
 *              about where windows should be placed etc.
 *              This file allows the graph package to be written without
 *              the need for #ifdef ACEDB macros within the graph library
 *              to change its behaviour.
 * Exported functions:
 * HISTORY:
 * Last edited: Oct  8 16:36 1998 (edgrif)
 * Created: Wed Sep 16 11:58:55 1998 (edgrif)
 *-------------------------------------------------------------------
 */

#include "graph_.h"

KEY isGifDisplay ; /* used to unset useless buttons etc */

/* Acedb function to show message in acedb main window.                      */
static VoidCharRoutine  graphAcedbMainActivityFunc = NULL ;


/* Strings used to search for various display types                          */
static GraphCharRoutine graphAcedbDisplayCreateFunc = NULL ;
static char *graphAcedbChronoName = NULL ;
static char *graphAcedbViewName = NULL ;
static char *graphAcedbFilqueryName = NULL ;


/* The printer menu that appears when you click on the print window.         */
static MENUOPT *graphAcedbPDMenu = NULL ;
static CharVoidRoutine graphAcedbStyleFunc = NULL ;
static char *graphAcedbLogin = NULL ;


/* Postscript file printing.                                                 */
static char *graphAcedbSessionName = NULL ;
static VoidFILEStackIntRoutine graphAcedbGetFonts = NULL ;


/* BoxInfo routine to print extra acedb specific information.                */
static VoidFILEKEYRoutine graphAcedbClassPrintFunc = NULL ;


/* graph to Xlib specific stuff.                                             */
static VoidIntRoutine graphAcedbXFontFunc = NULL ;

/* graph to Xt specific stuff (mostly to do with wcomm functions).           */
static VoidDisplayWindowRoutine graphAcedbWcommInitFunc = NULL ;
static VoidVoidRoutine graphAcedbWcommFinishFunc = NULL ;
static BoolXeventRoutine graphAcedbWPropFunc = NULL; 



/* viewedit routines specific to acedb.                                      */
static VoidVIEWCOLCONTROLRoutine graphAcedbViewAnonFunc = NULL ;
static BOOLMENURoutine graphAcedbSaveViewFunc = NULL ;
static VoidMAPCONTROLVIEWRoutine graphAcedbSaveViewMsgFunc = NULL ;
static VoidMENUSPECRoutine graphAcedbViewMenuFunc = NULL ;


/* colcontrol routines specific to acedb.                                    */
static BOOLMAPCONTROLKEYRoutine graphAcedbSetViewFunc = NULL ;
static VoidKEYRoutine graphAcedbResetKeyFunc = NULL ;
static VoidCOLOBJRoutine graphAcedbResetSaveFunc = NULL ;
static VoidIntFREEOPT graphAcedbFreemenuFunc = NULL ;
static VoidOBJFUNCMAPCONTROL graphAcedbMapLocateFunc = NULL ;
static VoidOBJFUNCMAPCONTROL graphAcedbScaleFunc = NULL ;
static VoidOBJFUNCSPACERPRIV graphAcedbSpaceFunc = NULL ;
static VoidCharKEY graphAcedbAddMapFunc = NULL ;


/* Initialise the overloaded function pointers.                              */
/*                                                                           */
GRAPH_FUNC_DEF void initOverloadFuncs()
  {

  /* No general overload functions currently...                              */

  return ;
  }


/* Set/Get routines to access the pointers to the overloadable functions.    */
/* NOTE, set routines are for use externally to the graph package,           */
/*       get routines are for use only in the graph package.                 */

GRAPH_FUNC_DEF VoidCharRoutine setGraphAcedbMainActivity(VoidCharRoutine funcptr)
  {
  VoidCharRoutine old = graphAcedbMainActivityFunc ;
  graphAcedbMainActivityFunc = funcptr ;
  return old ;
  }

GRAPH_FUNC_DEF VoidCharRoutine getGraphAcedbMainActivity(void) { return graphAcedbMainActivityFunc ; }


/* Get/set the names for various displays, names are used by acedb code      */
/* to set up display types via the registered acedb display function.        */
/*                                                                           */
GRAPH_FUNC_DEF GraphCharRoutine setGraphAcedbDisplayCreate(GraphCharRoutine funcptr,
							   char *chrono, char *view, char *filquery)
  {
  GraphCharRoutine old = graphAcedbDisplayCreateFunc ;
  graphAcedbDisplayCreateFunc = funcptr ;

  graphAcedbChronoName = chrono ;
  graphAcedbViewName = view ;
  graphAcedbFilqueryName = filquery ;

  return old ;
  }

GRAPH_FUNC_DEF GraphCharRoutine getGraphAcedbDisplayCreate(void) { return graphAcedbDisplayCreateFunc ; }
GRAPH_FUNC_DEF char *getGraphAcedbChronoName(void) { return graphAcedbChronoName ; }
GRAPH_FUNC_DEF char *getGraphAcedbViewName(void) { return graphAcedbViewName ; }
GRAPH_FUNC_DEF char *getGraphAcedbFilqueryName(void) { return graphAcedbFilqueryName ; }


/* Set/get the graphprint stuff.                                             */
GRAPH_FUNC_DEF void setGraphAcedbPDMenu(MENUOPT *new_menu, CharVoidRoutine stylefunc, char *login)
 {
 graphAcedbPDMenu = new_menu ;
 graphAcedbStyleFunc = stylefunc ;
 graphAcedbLogin = login ;

 return ;
 }
GRAPH_FUNC_DEF MENUOPT *getGraphAcedbPDMenu(void) { return graphAcedbPDMenu ; }
GRAPH_FUNC_DEF CharVoidRoutine getGraphAcedbStyle(void) { return graphAcedbStyleFunc ; }
GRAPH_FUNC_DEF char *getGraphAcedbLogin(void) { return graphAcedbLogin ; }


/* Set/get the graphps stuff.                                                */
GRAPH_FUNC_DEF void setGraphAcedbPS(char *session_name, VoidFILEStackIntRoutine fontfunc)
 {
 graphAcedbSessionName = session_name ;
 graphAcedbGetFonts = fontfunc ;

 return ;
 }
GRAPH_FUNC_DEF char *getGraphAcedbSessionName(void) { return graphAcedbSessionName ; }
GRAPH_FUNC_DEF VoidFILEStackIntRoutine getGraphAcedbGetFonts(void) { return graphAcedbGetFonts ; }



/* Set/get boxInfo Class print routine.                                      */
/*                                                                           */
GRAPH_FUNC_DEF VoidFILEKEYRoutine setGraphAcedbClassPrint(VoidFILEKEYRoutine classprintfunc)
  {
  VoidFILEKEYRoutine old = graphAcedbClassPrintFunc ;

  graphAcedbClassPrintFunc = classprintfunc ;

  return old ;
  }
GRAPH_FUNC_DEF VoidFILEKEYRoutine getGraphAcedbClassPrint(void) { return graphAcedbClassPrintFunc ; }




/* Set/get X font routine                                                    */
/*                                                                           */
GRAPH_FUNC_DEF VoidIntRoutine setGraphAcedbXFont(VoidIntRoutine xfontfunc)
  {
  VoidIntRoutine old = graphAcedbXFontFunc ;

  graphAcedbXFontFunc = xfontfunc ;

  return old ;
  }
GRAPH_FUNC_DEF VoidIntRoutine getGraphAcedbXFont(void) { return graphAcedbXFontFunc ; }



/* Set/get for Xt routines, mostly to do with setting up wcomm.              */
/*                                                                           */
GRAPH_FUNC_DEF void setGraphAcedbXWcomm(VoidDisplayWindowRoutine wcomminitfunc,
					VoidVoidRoutine wcommfinishfunc,
					BoolXeventRoutine wpropfunc)
  {
  graphAcedbWcommInitFunc = wcomminitfunc ;
  graphAcedbWcommFinishFunc = wcommfinishfunc ;
  graphAcedbWPropFunc = wpropfunc ;

  return ;
  }
GRAPH_FUNC_DEF VoidDisplayWindowRoutine getGraphAcedbWcommInit(void) { return graphAcedbWcommInitFunc ; } 
GRAPH_FUNC_DEF VoidVoidRoutine getGraphAcedbWcommFinish(void) { return graphAcedbWcommFinishFunc ; } 
GRAPH_FUNC_DEF BoolXeventRoutine getGraphAcedbWprop(void) { return graphAcedbWPropFunc ; }




/* Set/Get routines for viewedit functions.                                  */
/*                                                                           */
GRAPH_FUNC_DEF void setGraphAcedbView(VoidVIEWCOLCONTROLRoutine anonfunc,
				      BOOLMENURoutine saveviewfunc,
				      VoidMAPCONTROLVIEWRoutine viewsavemsgfunc,
				      VoidMENUSPECRoutine viewmenufunc)
  {
  graphAcedbViewAnonFunc = anonfunc ;
  graphAcedbSaveViewFunc = saveviewfunc ;
  graphAcedbSaveViewMsgFunc = viewsavemsgfunc ;
  graphAcedbViewMenuFunc = viewmenufunc ;

  return ;
  }
GRAPH_FUNC_DEF VoidVIEWCOLCONTROLRoutine getGraphAcedbViewAnon(void) { return graphAcedbViewAnonFunc ; } 
GRAPH_FUNC_DEF BOOLMENURoutine getGraphAcedbSaveView(void) { return graphAcedbSaveViewFunc ; }
GRAPH_FUNC_DEF VoidMAPCONTROLVIEWRoutine getGraphAcedbSaveViewMsg(void)
                                                             { return graphAcedbSaveViewMsgFunc ; }

GRAPH_FUNC_DEF VoidMENUSPECRoutine getGraphAcedbViewMenu(void) { return graphAcedbViewMenuFunc ; }


/* Set/Get routines for colcontrol functions.                                */
/*                                                                           */
GRAPH_FUNC_DEF void setGraphAcedbColControl(BOOLMAPCONTROLKEYRoutine setviewfunc,
					    VoidKEYRoutine resetkeyfunc,
					    VoidCOLOBJRoutine resetsavefunc,
					    VoidIntFREEOPT freemenufunc,
					    VoidOBJFUNCMAPCONTROL maplocatefunc,
					    VoidOBJFUNCMAPCONTROL scalefunc,
					    VoidOBJFUNCSPACERPRIV spacefunc,
					    VoidCharKEY addmap)
  {
  graphAcedbSetViewFunc = setviewfunc ;
  graphAcedbResetKeyFunc = resetkeyfunc ;
  graphAcedbResetSaveFunc = resetsavefunc ;
  graphAcedbFreemenuFunc = freemenufunc ;
  graphAcedbMapLocateFunc = maplocatefunc ;
  graphAcedbScaleFunc = scalefunc ;
  graphAcedbSpaceFunc = spacefunc ;
  graphAcedbAddMapFunc = addmap ;
  return ;
  }
GRAPH_FUNC_DEF BOOLMAPCONTROLKEYRoutine getGraphAcedbGetSetView(void) { return graphAcedbSetViewFunc ; }
GRAPH_FUNC_DEF VoidKEYRoutine getGraphAcedbResetKey(void) { return graphAcedbResetKeyFunc ; }
GRAPH_FUNC_DEF VoidCOLOBJRoutine getGraphAcedbResetSave(void) { return graphAcedbResetSaveFunc ; }
GRAPH_FUNC_DEF VoidIntFREEOPT getGraphAcedbFreemenu(void) { return graphAcedbFreemenuFunc ; }
GRAPH_FUNC_DEF VoidOBJFUNCMAPCONTROL getGraphAcedbMapLocate(void) { return graphAcedbMapLocateFunc ; }
GRAPH_FUNC_DEF VoidOBJFUNCMAPCONTROL getGraphAcedbScale(void) { return graphAcedbScaleFunc ; }
GRAPH_FUNC_DEF VoidOBJFUNCSPACERPRIV getGraphAcedbSpace(void) { return graphAcedbSpaceFunc ; }
GRAPH_FUNC_DEF VoidCharKEY getGraphAcedbAddMap(void) { return graphAcedbAddMapFunc ; }
