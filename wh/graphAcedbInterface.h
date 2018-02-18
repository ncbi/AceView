/*  File: graphAcedbInterface.h
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Defines functions to set up the interface between the
 *              graph package and acedb which allows acedb to alter
 *              the behaviour of the graph package.
 *              
 *              NOTE, this is not a general interface, it works for
 *              acedb but that is no guarantee it will work for
 *              anything else...you have been warned...
 *              
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 13 11:42 1998 (edgrif)
 * * Oct  8 10:01 1998 (edgrif): 
 * Created: Thu Sep 24 13:07:45 1998 (edgrif)
 *-------------------------------------------------------------------
 */


#ifndef DEF_GRAPHACEDBINTERFACE_H
#define DEF_GRAPHACEDBINTERFACE_H

#include <X11/Xlib.h>					    /* for Display/Window */

#include "regular.h"
#include "graph.h"
#include "colcontrol_.h"
#include "../wh/menu.h"



/********************  Graph  <-->  Acedb interface  *************************/
/* Do not use this interface for anything other than setting up graph for    */
/* use with Acedb. It is here to enable the two to cooperate, not for others */
/* to make use of...you have been warned.                                    */
/*                                                                           */


/*---------------------------------------------------------------------------*/
/* Intialisation of the  Graph  <-->  Acedb  interface                       */
/*                                                                           */
							    /* Function types:                   */
typedef void  (*VoidCharRoutine)(char*) ;		    /* void funcname(char*) */
typedef Graph (*GraphCharRoutine)(char*) ;		    /* Graph funcname(char*)*/
typedef char  (*CharVoidRoutine)(void) ;		    /* char funcname(void) */
typedef void  (*VoidFILEStackIntRoutine)(FILE**, Stack, int*) ;
							    /* void funcname(FILE**, Stack, int*) */
typedef void  (*VoidFILEKEYRoutine)(FILE *, KEY) ;	    /* void funcname(FILE*, KEY) */
typedef void  (*VoidIntRoutine)(int *level) ;		    /* void funcname(int*) */
typedef void  (*VoidDisplayWindowRoutine)(Display*, Window) ; /* void funcname(Display, Window) */
typedef void  (*VoidVoidRoutine)(void) ;		    /* void funcname(void) */
typedef BOOL  (*BoolXeventRoutine)(XEvent*) ;		    /* BOOL funcname(XEvent*) */
typedef void  (*VoidVIEWCOLCONTROLRoutine)(VIEWWINDOW, COLCONTROL) ;
							    /* void funcname(VIEWWINDOW, COLCONTROL) */
typedef BOOL  (*BOOLMENURoutine)(MENU) ;

typedef BOOL  (*BOOLMAPCONTROLKEYRoutine)(MAPCONTROL, KEY) ;

typedef void  (*VoidMAPCONTROLVIEWRoutine)(MAPCONTROL, VIEWWINDOW) ;

typedef void  (*VoidMENUSPECRoutine)(MENUSPEC**) ;

typedef void  (*VoidKEYRoutine)(KEY*) ;

typedef void  (*VoidCOLINSTANCEOBJ)(COLINSTANCE, OBJ) ;	    /* Note the linked typedef. */
typedef void  (*VoidCOLOBJRoutine)(VoidCOLINSTANCEOBJ*) ;

typedef void  (*VoidIntFREEOPT)(int, FREEOPT*) ;

typedef void  (*VoidOBJFUNCMAPCONTROL)(OBJ, VoidCOLINSTANCEOBJ*, MAPCONTROL) ;

typedef void  (*VoidOBJFUNCSPACERPRIV)(OBJ, VoidCOLINSTANCEOBJ*, SPACERPRIV) ;

typedef void  (*VoidCharKEY)(char*, KEY*) ;


GRAPH_FUNC_DEF VoidCharRoutine  setGraphAcedbMainActivity(VoidCharRoutine mainactivityfunc) ;

GRAPH_FUNC_DEF GraphCharRoutine setGraphAcedbDisplayCreate(GraphCharRoutine displayfunc,
							   char *chrono, char *view, char *filquery) ;

GRAPH_FUNC_DEF void setGraphAcedbPDMenu(MENUOPT *new_menu, CharVoidRoutine stylefunc, char *login) ;
void pdExecute (void) ;					    /* Not good, this is really internal
							       to the graph package in graphprint.c */


GRAPH_FUNC_DEF void setGraphAcedbPS(char *session_name, VoidFILEStackIntRoutine fontfunc) ;


GRAPH_FUNC_DEF VoidFILEKEYRoutine setGraphAcedbClassPrint(VoidFILEKEYRoutine classprintfunc) ;

GRAPH_FUNC_DEF VoidIntRoutine setGraphAcedbXFont(VoidIntRoutine xfontfunc) ;

GRAPH_FUNC_DEF void setGraphAcedbXWcomm(VoidDisplayWindowRoutine wcomminitfunc,
					VoidVoidRoutine wcommfinishfunc,
					BoolXeventRoutine wpropfunc) ;

GRAPH_FUNC_DEF void setGraphAcedbView(VoidVIEWCOLCONTROLRoutine anonfunc,
				      BOOLMENURoutine saveviewfunc,
				      VoidMAPCONTROLVIEWRoutine viewsavemsgfunc,
				      VoidMENUSPECRoutine viewmenufunc) ;


GRAPH_FUNC_DEF void setGraphAcedbColControl(BOOLMAPCONTROLKEYRoutine setviewfunc,
					    VoidKEYRoutine resetkeyfunc,
					    VoidCOLOBJRoutine resetsavefunc,
					    VoidIntFREEOPT freemenufunc,
					    VoidOBJFUNCMAPCONTROL maplocatefunc,
					    VoidOBJFUNCMAPCONTROL scalefunc,
					    VoidOBJFUNCSPACERPRIV spacefunc,
					    VoidCharKEY addmap) ;


/*---------------------------------------------------------------------------*/
/* Graph routines specific to the   Graph  <-->  Acedb  interface            */
/*                                                                           */

void controlMakeColumns(MAPCONTROL map, COLDEFAULT colInit, KEY view, BOOL) ;
int controlGetColour(OBJ obj);
void controlSetColour(OBJ obj, int colour);



/*****************************************************************************/
#endif /* DEF_GRAPHACEDBINTERFACE_H */
