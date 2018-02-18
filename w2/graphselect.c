/*  File: graphselect.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 **  New version of graphSelect written using the graphAce toolBox.
 * Exported functions:
          BOOL graphSelect
 * HISTORY:
 * Last edited: May 12 13:51 1999 (edgrif)
 * * May 12 13:45 1999 (edgrif): Added graphSetBlockMode for blocked menus.
 * * Dec 17 09:36 1998 (edgrif): change var. called 'main' to mainGraph (!)
 * * Sep 29 14:21 1998 (edgrif): Removed #ifdef ACEDB, instead the
 *              caller can now select two different display modes for
 *              graphSelect (one is the old '#ifdef ACEDB' stuff).
 * Created: Mon Jul 20 22:25:49 1992 (mieg)
 *-------------------------------------------------------------------
 * gha  04/30/93  change spacing, bold
 * gha  05/03/93  change title, layout
 */

/* $Id: graphselect.c,v 1.1.1.1 2002/07/19 20:22:55 sienkiew Exp $ */

#include "acedb.h"
#include "graph_.h"
#include "keyset.h"

static void popPick(int box) ;
static void popCancel(void) ;
static void popPlainDraw(FREEOPT *f) ;
static void popScrolledDraw(FREEOPT *f) ;


static GraphSelectDisplayMode display_mode = GRAPH_SELECT_SCROLLED ;
static KEYSET box2key = 0 ;
static KEY popChoice ;
static int minBox ;
static Graph mainGraph = 0 ;
static MENUOPT popMenu[] = {{popCancel, "Cancel"}, {0, 0}} ;


/*---------------------------------------------------------------------------*/
/* External routines                                                         */
/*                                                                           */


/* Set mode to be either scrolled or non-scrolled selections.                */
/*                                                                           */
GRAPH_FUNC_DCL void graphSetSelectDisplayMode(GraphSelectDisplayMode mode)
  {
  display_mode = mode ;

  return ;
  }



GRAPH_FUNC_DCL BOOL graphSelect (KEY *kpt, FREEOPT *f)
  { 
  Graph oldGraph = graphActive() ;
  BOOL result = FALSE ;


  /* ACEDB-GRAPH INTERFACE: if acedb has registered a message function then  */
  /* call it.                                                                */
  /* The string below should be in the acedb code...                         */
  if (getGraphAcedbMainActivity() != NULL)
            (getGraphAcedbMainActivity())("Please answer Selector") ;


  /* Select style of display.                                                */
  if (display_mode == GRAPH_SELECT_SCROLLED) popScrolledDraw(f) ;
  else popPlainDraw(f) ;
  
  graphActivate(mainGraph);

  result = graphBlock() ; 

  if (graphActivate(mainGraph)) graphDestroy();

  mainGraph = 0 ;
  graphActivate(oldGraph) ;

  if (result) *kpt = popChoice ;

  return result ;
  }


/*---------------------------------------------------------------------------*/
/* Internal routines.                                                        */
/*                                                                           */


static void popPick(int box)
{ 
  if (box <= minBox)
    return ;
  popChoice = keySet(box2key, box);
  graphActivate(mainGraph);
  graphUnBlock((int) popChoice ) ; 

  return;
}


static void popCancel(void)
{
  popChoice = 0 ;
  graphActivate(mainGraph);
  graphUnBlock(0) ; 
}


/* Display the selections in their own scrolled subwindow.                   */
/*                                                                           */
static void popScrolledDraw(FREEOPT *f)
{
  int i, max = f->key, maxLen = 30, px, py ;
  KEY key ;

  box2key = keySetReCreate(box2key) ;
  mainGraph = graphCreate(TEXT_FIT,
	      messprintf("Selector : %s", f->text), .3, .3, .25, .3) ;

  graphSetBlockMode(GRAPH_BLOCKING) ;

  graphButton("CANCEL", (GraphFunc) popCancel, 1, 0.5) ;
  graphTextFormat(BOLD);           /* gha debug 4/30/93 */
  graphText("Single-click to select", 1, 2.2) ;
  graphFitBounds(&px, &py);
  graphMenu (popMenu) ;

  graphSubCreate(TEXT_SCROLL, 1, 4, px - 2, py - 4 );
  graphRegister (PICK, popPick) ;
  graphMenu (popMenu) ;
  graphTextFormat(PLAIN_FORMAT);   

  minBox = graphBoxStart() ;
  for (i = 0, f++ ; i  < max ; f++, i++)
    { keySet(box2key, graphBoxStart()) = key = f->key ;
      graphText(f->text, 1, i) ;
      if (strlen(f->text) > maxLen)
	maxLen = strlen(f->text) ;
      graphBoxEnd() ;
    }
  graphBoxEnd() ;
  graphBoxMenu (minBox, popMenu);
  graphTextBounds(maxLen + 3, max + 4) ;
  graphRedraw() ;

  graphActivate(mainGraph);
  graphRedraw();
}


/* Display the selections in a single plain window the whole of which can    */
/* be scrolled.                                                              */
/*                                                                           */
static void popPlainDraw(FREEOPT *f)
{
  int i, max = f->key, maxLen = 30 ;
  KEY key ;

  box2key = keySetReCreate(box2key) ;
  mainGraph = graphCreate(TEXT_SCROLL,
	      messprintf("Selector : %s", f->text), .3, .3, .25, .3) ;

  graphSetBlockMode(GRAPH_BLOCKING) ;

  graphRegister (PICK, popPick) ;
  graphMenu (popMenu) ;

  graphTextFormat(BOLD);           /* gha debug 4/30/93 */
  graphText("Select or ", 1, 0.1) ;
  graphTextFormat(PLAIN_FORMAT);   
  graphButton("CANCEL", (GraphFunc) popCancel, 12, 0.1) ;

  minBox = graphBoxStart() ;
  for (i = 0, f++ ; i < max ; f++, i++)
    { keySet(box2key, graphBoxStart()) = key = f->key ;

      graphText(f->text, 1, i+3) ;

      if (strlen(f->text) > maxLen)
	maxLen = strlen(f->text) ;
      graphBoxEnd() ;
    }
  graphBoxEnd() ;
  graphBoxMenu(minBox, popMenu);

  graphTextBounds (maxLen + 3, max + 3) ;

  graphRedraw() ;

}



/****************************************************/
/****************************************************/
