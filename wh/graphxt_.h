/*  Last edited: Jun 10 15:00 1996 (rd) */

/* $Id: graphxt_.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */
/* graphxt_.h
   version of graph.h for inclusion within graphics package in
     files that use X11 (Athena Widgets and Intrinsics) calls.
   contains device dependent private information and then 
     includes graph_.h, which declares device independent private 
     information and includes graph.h (public information).
*/

#include<X11/IntrinsicP.h> /* need a private header so that we can change 
			      the geometry manager for the Composite - maratb*/
#include<X11/StringDefs.h>
#include<X11/Shell.h>

struct SubDevStruct
{
  BOOL   isExposed ;    /* do not draw until exposed */

  Widget viewport;      /* viewport onto base for scrolling */
  Widget base;          /* base widget */

  Widget menu0 ;
  Associator boxMenus ;
  Associator images ;
  void*  message ;	/* actually a Graph */
};

struct DevStruct
{
  Widget popup ;		/* outer parent */
  struct SubDevStruct *subdev ;	/* first subdev */
} ;

extern int menuBox ;		/* used to transmit ID of menu box */

extern XtAppContext app_con;	/* application context */
extern Widget root_widget;	/* parent of all */
extern void promptReturn (Widget w, XEvent* ev, String *params, Cardinal *num_params) ;

/* from graphxlib.c */
extern Pixmap GreyPixmap ;
Colormap getGlobalColorMap (void) ;

#include <X11/Xatom.h>
extern Atom myExitAtom ;

#define DEV_DEFINED

#include "graph_.h"
 
/****** end of file ******/
 
