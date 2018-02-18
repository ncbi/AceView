/*  File: xtsubs.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *      with help from Christopher Lee in this module.
 *
 * Description: popUp's for X version of graph package
 * Exported functions: graphOut, Query, Prompt, Menu and FreeMenu
 * HISTORY:
 * Last edited: May 13 10:37 1999 (edgrif)
 * * May 13 10:36 1999 (edgrif): Make graphWinPrompt blocking for busy cursor.
 * * Mar 10 10:25 1999 (edgrif): Removed XRaiseWindow/XMapWindow from localLoop
 *              and added proper visibility check (SANgc03382).
 * * Dec 17 09:54 1998 (edgrif): Added busy cursor to menu routines.
 * * Oct 26 18:17 1998 (fw): thicker outline for background menus
 * * Oct 22 14:24 1998 (edgrif): Added a series of message display calls
 *              that are called by the device independent layer of graph,
 *              this layering is unfortunately all broken.
 * * Oct 21 13:54 1998 (edgrif): Removed isGraphics, superceeded by
 *              function callback mechanism for messubs (see mieg below).
 * * Jun 10 00:24 1992 (mieg): graphOut = printf if !isGraphics
 * * Mar 31 19:00 1992 (mieg): added castings to please IBM
 * Created: Thu Jan  2 02:11:04 1992 (rd)
 *-------------------------------------------------------------------
 */

#include "regular.h"

#include "menu_.h"
#include "graphxt_.h"

#include <X11/Xaw/Dialog.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/SimpleMenu.h>
#include <X11/Xaw/SmeBSB.h>
#include <X11/Xaw/SmeLine.h>

static Arg args[10];		/* convenience variable */

static Widget popshell, dialog ;
	/* I tried keeping these constant and mapping/unmapping
	   components, but the calls are too variable, so the
	   geometry got into a tangle.  Domage */

static BOOL isLocalLoop = FALSE ;
static BOOL result ;


static void makeVisible(Widget w, XtPointer client_data, XEvent *event, Boolean *propagate) ;


/****** basic dialog box for everything but menus ******/

static void makeDialog (const char *text)
{
  Position x,y ;
  int n ;
  char *text2 = strnew (text, 0) ;;

  if (gDev && gSubDev)
    { XtTranslateCoords (gSubDev->viewport, 
			 (Position) 25, (Position) 20,
			 &x, &y) ;
      if (x < 25) x = 25 ;
      if (x > 900) x = 900 ;
      if (y < 20) y = 20 ;
      if (y > 700) y = 700 ;
    }
  else
    { x = 500 ; y = 400 ; }

  n = 0;
  XtSetArg (args[n], XtNx, x) ;	n++ ;
  XtSetArg (args[n], XtNy, y) ;	n++ ;

  popshell = XtCreatePopupShell ("Please Reply", 
				 transientShellWidgetClass,
				 root_widget, args, n) ;
  dialog = XtCreateManagedWidget ("dialog", dialogWidgetClass,
				  popshell, NULL, 0) ;
  XtVaSetValues (dialog, XtNlabel, uBrokenText (text2, 40), NULL) ;

  /* This event handler will keep the dialog visible to the user.            */
  XtAddEventHandler(popshell, VisibilityChangeMask, FALSE, makeVisible, NULL) ;

  XBell (XtDisplay(popshell), 0) ;
  messfree (text2) ;
}


/* Keep the message dialog on top, this is important because the dialog may  */
/* be blocking and if the dialog gets hidden under other windows the user    */
/* will be left with an xace that does not respond.....                      */
static void makeVisible(Widget w, XtPointer client_data, XEvent *event, Boolean *propagate)
  {
  XVisibilityEvent *visible = (XVisibilityEvent *)event ;

  if (visible->state == VisibilityPartiallyObscured || visible->state == VisibilityFullyObscured)
    {
    XRaiseWindow(XtDisplay(w), XtWindow(w)) ;
    }

  return ;
  }


static void endAction (Widget w, XtPointer client_data, XtPointer call_data)
{
  XtPopdown (popshell) ;
  isLocalLoop = FALSE ;
  result = (client_data != 0) ;
}


static void localLoop (void)
{ 
  XEvent ev ;
  XtPopup (popshell, XtGrabExclusive) ;
  /*mhmp 24.11.97*/
  /* Widget w = popshell ; */
  isLocalLoop = TRUE ;
  while (isLocalLoop)
    {
      XtAppNextEvent (app_con, &ev) ;
      XtDispatchEvent (&ev) ;
    }
  XtDestroyWidget (popshell) ;
}
/* mhmp*/



/* This set of routines provide the device dependent layer that the high     */
/* level functions of graph expect to be able to call, e.g. graphOut calls   */
/* graphWinOut etc.                                                          */
/* These routines must _NOT_ be called by application code, this breaks the  */
/* encapsulation of the graph package.                                       */
/*                                                                           */
void graphWinOut (const char *text, char *label)
{
  Widget button ;

  makeDialog (text) ;

  if (label == NULL) label = "Continue" ;		    /* Crummy default. */

  button = XtCreateManagedWidget (label, commandWidgetClass, dialog, NULL, 0) ;
  XtAddCallback (button, XtNcallback, endAction, (XtPointer) FALSE) ;
      
  localLoop () ;

  return ;
}

BOOL graphWinQuery(const char *text)
{
  Widget yes, no ;

  makeDialog (text) ;

  yes = XtCreateManagedWidget ("Yes", commandWidgetClass, dialog, NULL, 0) ;
  XtAddCallback (yes, XtNcallback, endAction, (XtPointer) TRUE) ;

  no = XtCreateManagedWidget ("No", commandWidgetClass, dialog, NULL, 0) ;
  XtAddCallback (no, XtNcallback, endAction, (XtPointer) FALSE) ;

  localLoop () ;

  return result ;
}


void promptReturn (Widget w, XEvent* ev, String *params, Cardinal *num_params)
{
  isLocalLoop = FALSE ;
  result = TRUE ;
}



#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
BOOL graphBadPrompt (char *prompt, char *dfault, char *fmt)
{
  Widget ok, cancel ;

  makeDialog (prompt) ;

  ok = XtCreateManagedWidget ("OK", commandWidgetClass, dialog, NULL, 0) ;
  XtAddCallback (ok, XtNcallback, endAction, (XtPointer) TRUE) ;

  cancel = XtCreateManagedWidget ("Cancel", commandWidgetClass, dialog, NULL, 0) ;
  XtAddCallback (cancel, XtNcallback, endAction, (XtPointer) FALSE) ;

  while (TRUE)	/* iterate on format mismatch */
    {
      XtVaSetValues (dialog, XtNvalue, dfault, NULL) ;
      
      localLoop () ;
      
      if (!result)
	return FALSE ;
      freeforcecard (dfault = XawDialogGetValueString (dialog)) ;
      if (freecheck (fmt))
	return TRUE ;
      XtVaSetValues (dialog, XtNlabel,
       "Sorry, invalid response. Try again or cancel", NULL) ;
    }
}
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */


/*** another version of graphPrompt based on graph package ***/

static void okReturn (void) { graphUnBlock(TRUE) ; }
static void cancelReturn (void) { graphUnBlock (FALSE) ;}
static void entryProcess (char *text) { graphUnBlock (TRUE) ; }

BOOL graphWinPrompt (const char *prompt, char *dfault, char *fmt)
     /* version written in graph package because X doesn't work! */
{
  char buf[4096], *cp, *pp ;
  float line = 0.5 ;
  Graph gsave = gActive ? gActive->id : 0 ;
  static MENUOPT buttons[] = { {okReturn, "OK"},
			       {cancelReturn, "Cancel"},
				{0, 0}
			     } ;

  graphCreate (TEXT_SCROLL, "Please Reply",
			 0.5,0.4,0.45,0.25) ;

  graphSetBlockMode(GRAPH_BLOCKING) ;

  graphHelp ("Prompt_Window") ; /* mieg 18-8-94 */

  pp = strnew (prompt, 0) ;
  uLinesText (pp, 40) ;
  while ((cp = uNextLine (pp)))
    graphText (cp, 1, line++) ;
  messfree (pp) ;
  memset (buf,0,4096) ;/* JTM otherwise by typing in you collect garbage */
  strncpy (buf, dfault, 4095) ;
  line += 0.5 ;

	/*	graphTextScrollEntry (buf, 1024, 40, 2, line++, entryProcess) ;*/
	graphTextScrollEntry (buf, 4095, 40, 2, line++, entryProcess) ;

  line += 0.5 ;
   graphButtons (buttons, 4, line++, 40) ; 
  line += 0.5 ;

  graphRedraw() ;

/* TRUE && ... circumvents a dec alpha optimiser bug ! */
  while (TRUE && (result = graphBlock())) /* iterate on format mismatch */
    {
      if (!result)
	{ graphActivate (gsave) ; return FALSE ; }
      freeforcecard (buf) ;
      if (freecheck (fmt))
	break ;
      graphText ("Sorry, invalid response. Try again or cancel", 1, line) ;
      graphRedraw() ;
    }

  graphDestroy () ;
  graphActivate (gsave) ;
  return result ;
}

/******* menus ********/

static Associator oldMenus, freeMenus, menuFunc ;

void devMenuDestroy(void) {}

static void registerMenu (int k, Widget menu)
{
  Box box = gBoxGet (k) ;
  union { int i;
	  void *v;
	} u;

  u.v = 0;
  if (k)
    { if (!gSubDev->boxMenus)
	gSubDev->boxMenus = assCreate() ;
      u.i = 4*k;
      assRemove (gSubDev->boxMenus, u.v) ;
      assInsert (gSubDev->boxMenus, u.v, menu) ;
    }
  else
    gSubDev->menu0 = menu ;

  box->flag |= GRAPH_BOX_MENU_FLAG ;
}

static void unregisterMenu (int k)
{
  Box box = gBoxGet (k) ;
  union { int i;
	  void *v;
	} u;

  u.v = 0;
  u.i = 4*k;

  if (!k)
    gSubDev->menu0 = 0 ;
  else if (gSubDev->boxMenus)
    {
      BOOL isPreviousMenu = assRemove (gSubDev->boxMenus, u.v) ;

      if (isPreviousMenu)
	/* There was a menu attached to this box before, so we
	   have removed it and re-instate the original state
	   of the box, as if no menu was ever attached to it.
	   The box will return to the default menu behaviour.
	   The background menu will be used again for this box */
	box->flag &= ~GRAPH_BOX_MENU_FLAG ;
      else
	/* There was NO menu attached to this box before.
	   This was probably caused by calling graphBoxMenu(box, 0)
	   on the box just after creating it. This box is therefor
	   marked to have no menu at all and will not behave
	   according to default menu behaviour (which displays the
	   background menu by default */
	box->flag |= GRAPH_BOX_MENU_FLAG ;
    }
}

/**********************/

static void freeMenuSelect (Widget w, XtPointer data, XtPointer junk)
{
  FreeMenuFunction func ;

  /* If there is a function attached to a menu item call it.                 */
  if (assFind (menuFunc, (void *)w, &func))
    {
    graphBusyCursorAll (TRUE);

    (*func)((KEY) assInt(data), menuBox) ;

    graphBusyCursorAll (FALSE) ;
    }

  return ;
}



static void oldMenuSelect (Widget w, XtPointer data, XtPointer junk)
  {
  graphBusyCursorAll (TRUE);

  (*(MenuFunction)data)(menuBox) ;

  graphBusyCursorAll (FALSE) ;

  return ;
  }



/*********************************/


/* test code......                                                           */
#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
/* These are the original translations lifted direct from the SimpleMenu code*/
static char defaultTranslations[] =
    "<EnterWindow>:     highlight()             \n\
     <LeaveWindow>:     unhighlight()           \n\
     <BtnMotion>:       highlight()             \n\
     <BtnUp>:           MenuPopdown() notify() unhighlight()"; 
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

static char newTranslations[] =
    "Lock<EnterWindow>:     highlight()             \n\
     Lock<LeaveWindow>:     unhighlight()           \n\
     Lock<BtnMotion>:       highlight()             \n\
     Lock<BtnUp>:           MenuPopdown() notify() unhighlight() \n\
     <EnterWindow>:     highlight()             \n\
     <LeaveWindow>:     unhighlight()           \n\
     <BtnMotion>:       highlight()             \n\
     <BtnUp>:           MenuPopdown() notify() unhighlight()"; 




void graphBoxMenu (int box, MENUOPT *options)
{
  struct menubits 
    { Widget widg;
      MENUOPT *opts;
    } *menu;
  MENUOPT *o;
  int i, n=0, iBucket = 0 ;
  Array bucket = 0 ;
  Widget entry ;

  if (!gDev || !gSubDev)
    return ;
  if (!options)
    { unregisterMenu (box) ;
      return ;
    }

  if (!oldMenus)
    oldMenus = assCreate() ;

  iBucket = 0 ; bucket = 0 ;
  if (assFind (oldMenus, (void *)(options->f), &menu))
    { while (assFindNext(oldMenus, (void *)(options->f), (const void **)&menu, &bucket, &iBucket))
	{ o = menu->opts;
	  for (i=0; o[i].f && options[i].f; i++)
	    if (o[i].f != options[i].f ||
		strcmp(o[i].text, options[i].text)) 
	      break ;
	  if (!o[i].f && !options[i].f)
	    {	  /* Found an existing menu here, resuse it */
	      registerMenu (box, menu->widg);
	      return;
	    }
	  continue;
	}
    }

  /*  make a new one here. */

  for (i=0; options[i].f; i++);
  o = (MENUOPT *) messalloc((i+1)*sizeof(MENUOPT));

  for (i=0; options[i].f; i++)
    { o[i].f = options[i].f;
      o[i].text =  options[i].text ? strnew (options[i].text, 0) : "" ;
    }
  o[i].f =  0;
  o[i].text = 0;
  menu = (struct menubits *) messalloc(sizeof(struct menubits));
  menu->opts = o;

  n = 0;
  if (box == 0)
    {
      /* give the background menu a thicker outline,
	 it therefor stands out from the boxmenus */
      XtSetArg (args[n], XtNborderWidth, 2) ;	n++ ;
    }
  if (1) {   XtSetArg(args[n], XtNwidth, 40) ; n++ ; }
  /* XtNheight */

  menu->widg = XtCreatePopupShell ("menu", simpleMenuWidgetClass, 
				   root_widget, args, n) ;

  XtVaSetValues(menu->widg,
		XtNtranslations,XtParseTranslationTable(newTranslations),
		NULL) ;

  assMultipleInsert (oldMenus, (void *)(options->f), menu) ;
  for (i=0; o[i].f; i++)
    { 
      if (o[i].f == menuSpacer &&
	  (strcmp(o[i].text, "") == 0))
	{
	  /* separator */
	  entry = XtCreateManagedWidget (o[i].text,
					 smeLineObjectClass,
					 menu->widg, NULL, 0) ;
	  continue;
	}

      entry = XtCreateManagedWidget (o[i].text,
				     smeBSBObjectClass,
				     menu->widg, NULL, 0) ;
      XtAddCallback (entry, XtNcallback,
		     oldMenuSelect, (XtPointer) o[i].f) ;
    }
  
  registerMenu (box, menu->widg) ;
}


void graphBoxFreeMenu (int box, FreeMenuFunction proc, FREEOPT *options)
{
  struct menubits 
    { Widget widg;
      FREEOPT *opts;
      FreeMenuFunction proc;
    } *menu;
  FREEOPT *o;
  int i, iBucket ;
  Array bucket = 0 ; 
  Widget entry ;
  int n;

  if (!gDev || !gSubDev)
    return ;
  if (!proc || !options)
    { unregisterMenu (box) ;
      return ;
    }

  if (!freeMenus)
    freeMenus = assCreate() ;

  iBucket = 0 ;
  bucket = 0 ;
  if (assFind (freeMenus, (void *)proc, &menu))
    { while (assFindNext(freeMenus, (void *)proc, (const void **)&menu, &bucket, &iBucket))
	{ o = menu->opts;
	  if (menu->proc != proc) goto notfound;
	  if (o->key != options->key ) goto notfound;
	  for (i=1; i <= o->key; i++)
	    { if (o[i].key != options[i].key) goto notfound;
	      if ( 0 != strcmp(o[i].text, options[i].text)) goto notfound;
	    }
	  /* Found an existing menu here, resuse it */
	  registerMenu (box, menu->widg);
	  return;
	notfound: continue;
	}
    }
  /*  make a new one here. */
  o = (FREEOPT *) messalloc(((options->key)+1)*sizeof(FREEOPT));

  for (i=0; i<= options->key; i++)
    { o[i].key = options[i].key;
      o[i].text = options[i].text ? strnew (options[i].text, 0) : "" ;
    }
  menu = (struct menubits *) messalloc(sizeof(struct menubits));
  menu->opts = o;
  menu->proc = proc;

  n = 0;
  if (box == 0)
    {
      /* give the background menu a thicker outline,
	 it therefor stands out from the boxmenus */
      XtSetArg (args[n], XtNborderWidth, 2) ;	n++ ;
    }

  menu->widg = XtCreatePopupShell ("menu", simpleMenuWidgetClass, 
				   root_widget, args, n) ;
  XtVaSetValues(menu->widg,
		XtNtranslations,XtParseTranslationTable(newTranslations),
		NULL) ;

  assMultipleInsert (freeMenus, (void *)proc, menu) ;

  messalloccheck();
  if (!menuFunc) menuFunc = assCreate() ;
  
  for (i=1; i <= o->key; i++)
    { union { Widget w;
	      void *v;
	    } u;
      u.v = 0;
      if (o[i].key == 0 && (strcmp(o[i].text, "") == 0))
	{
	  /* menu seperator */
	  u.w = entry = XtCreateManagedWidget (o[i].text,
					       smeLineObjectClass,
					       menu->widg, NULL, 0) ;
	  continue;
	}

      /* normal menu item is a text widget with a callback register */
      u.w = entry = XtCreateManagedWidget (o[i].text,
					   smeBSBObjectClass,
					   menu->widg, NULL, 0) ;

      XtAddCallback (entry, XtNcallback,
		     freeMenuSelect, (XtPointer) assVoid(o[i].key)) ;
      assInsert(menuFunc, u.v, (void *)proc) ;
    }
  
  registerMenu (box, menu->widg) ;
}

/******************************* new menus *************************/

static Associator newMenus ;

static void newMenuSelect (Widget w, XtPointer data, XtPointer junk)
{
  menuSelectItem ((MENUITEM)data) ;
}

void graphNewBoxMenu (int box, MENU menu)
{
  struct menuPlus 
    { MENU m ;
      Widget w ;
    } *plus ;
  MENUITEM item, item2 ;
  Widget entry ;
  char label[256] ;
  int n, iBucket = 0 ;
  Array bucket ;

  if (!gDev || !gSubDev) return ;
  if (!menu) { unregisterMenu (box) ; return ; }
  if (!newMenus) newMenus = assCreate() ;

  iBucket = 0 ;
  bucket = 0 ;
  if (assFind (newMenus, menu, &plus))
    { while (assFindNext (newMenus, menu, (const void **)&plus, &bucket, &iBucket))
	{ if (strcmp (menu->title, plus->m->title))
	    continue ;
	  for (item = menu->items, item2 = plus->m->items ; 
	       item && item2 ;
	       item = item->down, item2 = item2->down)
	    if (item->label != item2->label ||
		item->func != item2->func ||
		item->flags != item2->flags ||
		item->call != item2->call ||
		item->value != item2->value ||
		item->ptr != item2->ptr ||
		item->submenu != item2->submenu)
	      break ;
	  if (!item && !item2)	/* found an existing menu - reuse it */
	    { registerMenu (box, plus->w) ;
	      return ;
	    }
	}
    }
				/* make a new menu */
  plus = (struct menuPlus*) messalloc (sizeof (struct menuPlus)) ;
  plus->m = menuCopy (menu) ;

  n = 0;
  if (box == 0)
    {
      /* give the background menu a thicker outline,
	 it therefor stands out from the boxmenus */
      XtSetArg (args[n], XtNborderWidth, 2) ;	n++ ;
    }
  plus->w = XtCreatePopupShell ("menu", simpleMenuWidgetClass,
				root_widget, args, n) ;

  XtVaSetValues(plus->w,
		XtNtranslations,XtParseTranslationTable(newTranslations),
		NULL) ;

  if (plus->m->title)
    XtVaSetValues (plus->w, XtNlabel, plus->m->title, NULL) ;
  assMultipleInsert (newMenus, menu, plus) ;
  for (item = plus->m->items ; item ; item = item->down)
				/* important to use plus->m, not menu */
    { unsigned int flags = item->flags ;
      if (flags & MENUFLAG_HIDE) 
	continue ;
      if (flags & MENUFLAG_SPACER)
	{ entry = XtCreateManagedWidget (item->label,
					 smeLineObjectClass,
					 plus->w, NULL, 0) ;
	  continue ;
	}

      if (flags & MENUFLAG_TOGGLE_STATE)
	strcpy(label, "*") ;
      else
	strcpy(label, " ") ;
      strncat(label, item->label, 253) ;
      label[254] = 0 ;
      strcat(label, " ") ;

      entry = XtCreateManagedWidget (label,
				     smeBSBObjectClass,
				     plus->w, NULL, 0) ;
      if (flags & MENUFLAG_DISABLED)
	XtVaSetValues (entry, XtNsensitive, False, NULL) ;
      else
	XtAddCallback (entry, XtNcallback,
		       newMenuSelect, (XtPointer) item) ;
    }

  registerMenu (box, plus->w) ;
}

/************ end of file *************/
