/*  File: oxhomlist
 *  Written by Jo Dicks (jld@bioch.ox.ac.uk)
 *-------------------------------------------------------------------
 * This file is part of the ACeDB comparative mapping package, 
 * written by Jo Dicks and Michelle Kirby
 * The original ACeDB system was written by :
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *      Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *
 * Description:  to display lists of homology data from an Oxford Grid cell
 * HISTORY:
 * Created: January 1994 
 * Last edited: August 1994
 *-------------------------------------------------------------------------
 */

/* $Id: oxhomlist.c,v 1.2 2017/01/19 21:53:21 mieg Exp $ */

#include "acedb.h"
#include "keyset.h"
#include "graph.h"
#include "key.h"
#include "lex.h"
#include "bs.h"
#include "systags.h"
#include "classes.h"
#include "sysclass.h"
#include "tags.h"
#include "disptype.h"
#include "display.h"
#include "query.h"
#include "session.h"
#include "oxgrid.h"
#include "pmap.h"

static void hlDraw (void) ;
static void hlPick (int box) ;
static void hlSelect (OX ox, int box) ;
static void hlFollow (OX ox) ;
static void showTree (void) ;
static void hlKeyboard (int k) ;
static void hlExportAll (void) ;
static void hlExportM1 (void) ;
static void hlExportM2 (void) ;

enum BoxNames { BACKGROUND=0, EXPORT_BOX, MIN_BOX } ;

static MENUOPT homlistMenu[] =
            { {graphDestroy, "Quit"},
	      {help, "Help"},
              {showTree, "Show as Text"},
	      {graphPrint, "Print"},
	      {displayPreserve, "Preserve"},
	       {0, 0}
            } ;

static MENUOPT exportMenu[] =
            { {hlExportAll, "Export all loci"},
              {hlExportM1, "Export map 1 loci"},
              {hlExportM2, "Export map 2 loci"},
	       {0, 0}
            } ;

/*********************************/

void oxhomlist (void)
{
  int i, box, count ;
  SEG *seg ;

  OXGET ("oxhomlist") ;
  
  count = 0 ;
  ox->list = 0 ;
  box = ox->stayBox = ox->listBox = ox->gridBox ;

  for (i = 0 ; i <  arrayMax (ox->segs) ; i++)
    {seg = arrayp (ox->segs, i, SEG) ;
     if ((seg->box == box) && !(seg->flag & FLAG_HIDE))
     count++ ;}

/* Check if the cell is empty - if not, start a graph */

  if (count == 0)
    {messout ("No homologies within this cell") ;
     if (graphExists (ox->hlGraph)) 
       {graphActivate (ox->hlGraph) ;
        graphDestroy () ;} 
     }
  else
    {if (graphExists (ox->hlGraph))
       {graphActivate (ox->hlGraph) ;
        graphPop () ;}
     else 
       {if (!(ox->hlGraph = graphCreate (TEXT_SCROLL, "Homology List", .5, .15, .5, .7)))
        return ;
        ox->hlGraph = graphActive () ;

        graphAssociate (&OX_MAGIC, ox) ;      
        graphRegister (MESSAGE_DESTROY, displayUnBlock) ;
        graphRegister (PICK, (GraphFunc) hlPick) ;
        graphRegister (KEYBOARD, (GraphFunc) hlKeyboard) ;
        graphMenu (homlistMenu) ;
        graphTextBounds (40, 100) ;
       }

      hlDraw () ;
     }

}

static void hlDraw (void)
{
  int i, h, y, box, count ;
  int wmax, wmax1 = 0 , wmax2 = 0 ;
  KEY m1, m2 ;
  SEG *seg ;
  OXGET("hlDraw") ;

  graphClear () ;

  m1 = m2 = 0 ;
  box = ox->stayBox ;

  if (graphButton ("Export...", hlExportAll, 3, 3) != EXPORT_BOX)
    messcrash ("box screwup at %d in hlDraw", EXPORT_BOX) ;   

  graphBoxMenu (EXPORT_BOX, exportMenu) ;

  for (i = 0 ; i <  arrayMax (ox->segs) ; i++) 
    {seg = arrayp (ox->segs, i, SEG) ;
     seg->hbox = 0 ;
     if (seg->box == box && !(seg->flag & FLAG_HIDE))
       {wmax = strlen (name (seg->gene1)) ;
        if (wmax > wmax1) wmax1 = wmax ;
        wmax = strlen (name (seg->gene2)) ;
        if (wmax > wmax2) wmax2 = wmax ;
        m1 = seg->map1 ;
        m2 = seg->map2 ;
       }
     }

  wmax = strlen (name (ox->s1)) + 4 ;
  if (wmax > wmax1) wmax1 = wmax ;
  wmax = strlen (name (m1)) ;
  if (wmax > wmax1) wmax1 = wmax ;
  wmax = strlen (name (ox->s2)) + 4 ;
  if (wmax > wmax2) wmax2 = wmax ;
  wmax = strlen (name (m2)) ;
  if (wmax > wmax2) wmax2 = wmax ;
  wmax1 += 4 ; wmax2 += 4 ;
  
/* Draw the homology boxes */

  count = 0 ; y = 11 ; h = MIN_BOX ;
  for (i = 0 ; i < arrayMax (ox->segs) ; i++)
    {seg = arrayp (ox->segs, i, SEG) ;
     if (seg->box == box && !(seg->flag & FLAG_HIDE)) 
      {if (h != graphBoxStart ())
         messcrash ("problems with homlist %d", h) ;
       graphTextFormat (BOLD) ;
       graphText (name (seg->gene1), 5, y);
       graphBoxEnd () ;
       graphBoxDraw (h, BLACK, WHITE) ;
       if ((h+1) != graphBoxStart ())
         messcrash ("problems with homlist") ;
       graphTextFormat (BOLD) ;
       graphText (name (seg->gene2), 5 + wmax1, y);
       graphBoxEnd () ;
       graphBoxDraw (h + 1, BLACK, WHITE) ;
       seg->hbox = h ;
       h += 2 ;
       y += 1 ;
       count += 1 ;
       m1 = seg->map1 ;
       m2 = seg->map2 ;
       }
     }

  graphTextFormat(PLAIN_FORMAT);

/* Draw the headings */

  ox->list = count ;
  if (count == 1) graphText ("1 homology found", 3, 1) ;
  else {graphText (messprintf ("%d homologies found", count), 3, 1) ;}
  graphText (messprintf ("%s map", name (ox->s1)), 5, 7) ;
  graphText (messprintf ("%s map", name (ox->s2)), 5+wmax1, 7) ;
  graphText (messprintf ("%s", name(m1)), 5, 8) ;
  graphText (messprintf ("%s", name(m2)), 5+wmax1, 8) ;

/* Draw a box around the data */
  
  graphColor (BLUE) ;
  graphLine (3, 6, 3+wmax1+wmax2, 6) ;
  graphLine (3, 10, 3+wmax1+wmax2, 10) ;
  graphLine (3+wmax1, 6, 3+wmax1, y+1) ;
  graphLine (3, y+1, 3+wmax1+wmax2, y+1) ;
  graphLine (3, 6, 3, y+1) ;
  graphLine (3+wmax1+wmax2, 6, 3+wmax1+wmax2, y+1) ;
  graphColor (BLACK) ;

  ox->listBox = 0 ;

  graphRedraw () ;

}

static void hlKeyboard (int k)
{
  int box = 0, box2 ;

  OXGET ("hlKeyboard") ;

  if (!ox->listBox) return ;

  if (k == RETURN_KEY && (ox->listBox >= MIN_BOX && ox->listBox <= MIN_BOX + (ox->list)*2))
    {hlFollow (ox) ;
     return ;
    }

  if (ox->listBox >= MIN_BOX && ox->listBox <= MIN_BOX + (ox->list)*2)
    {
      switch (k)
	{
	case LEFT_KEY:
	  if (ox->listBox >= MIN_BOX-1 + 2) box = ox->listBox - 1 ;
	  else box = ox->listBox ;
	  break ;
	case RIGHT_KEY:
	  if (ox->listBox <= (MIN_BOX-1 + (ox->list)*2 - 1)) box = ox->listBox + 1 ;
	  else box = ox->listBox ;
	  break ;
	case UP_KEY:
	  if (ox->listBox >= MIN_BOX-1 + 3) box = ox->listBox - 2 ;
	  else box = ox->listBox ;
	  break ;
	case DOWN_KEY:
	  if (ox->listBox <= (MIN_BOX-1 + (ox->list)*2 - 2)) box = ox->listBox + 2 ;
	  else box = ox->listBox ;
	  break ;
	}
      box2 = ox->listBox ;
      ox->listBox = box ;
      hlFollow (ox) ;
      graphActivate (ox->hlGraph) ;
      if (ox->listBox) graphBoxDraw (box2, BLACK, WHITE) ;
      if (box)
	graphBoxDraw (box, WHITE, MAGENTA) ;
    }
}

static void hlPick (int box) 
{ 
  OXGET ("hlPick") ;

  /* Usual Pick routine */

      if (box >= 2 &&
	  box <= 2*arrayMax (ox->segs))           /* a seg */
	{ if (box == ox->listBox)
	    hlFollow (ox) ;
          graphActivate (ox->hlGraph) ;
          hlSelect (ox, box) ;
	}
          
}

static void hlSelect (OX ox, int box)
{
  /* Usual Select routine */

  if (ox->listBox == box) return ;

  if (ox->listBox) graphBoxDraw (ox->listBox, BLACK, WHITE) ;

  ox->listBox = box ;

  graphBoxDraw (ox->listBox, WHITE, MAGENTA) ;

}

static void hlFollow (OX ox)
{
  int i ;
  KEY hit = 0 ;
  SEG *seg ;

  for (i=0 ; i < arrayMax (ox->segs) ; i++)
    {
      seg = arrayp (ox->segs, i, SEG) ;
      if (seg->box == ox->stayBox && seg->hbox == ox->listBox)
	hit = seg->gene1 ;
      if (seg->box == ox->stayBox && !(seg->hbox == 0) && seg->hbox + 1 == ox->listBox )
	hit = seg->gene2 ;
    }
  if (hit)
    display (hit, 0, 0) ;
}

static void showTree (void)
{
  int i ;
  KEY hit = 0 ;
  SEG *seg ;
  OXGET ("showTree") ;

  /* Display a box as Text */
  
  if (!(ox->listBox))
    messout ("You have not selected a cell") ;
  else
    {
      for (i = 0 ; i < arrayMax (ox->segs) ; i++)
	{seg = arrayp (ox->segs, i, SEG) ;
	if (seg->box == ox->stayBox && seg->hbox == ox->listBox)
	  hit = seg->gene1 ;
       if (seg->box == ox->stayBox && seg->hbox + 1 == ox->listBox)
         hit = seg->gene2 ;
	}
      if (hit)
	display (hit, 0, TREE) ;
    }
}

static void hlExportAll (void) 
{ 
  int i, j = 0 ;
  SEG *seg ;
  static KEYSET kset = 0 ;
  
  OXGET ("hlExportAll");

  kset = keySetReCreate (kset) ;

  for (i = 0 ; i < arrayMax (ox->segs) ; ++i)
    { seg = arrp (ox->segs, i, SEG) ;
      if (seg->box == ox->stayBox && !(seg->flag & FLAG_HIDE))
        {keySet (kset, j++) = seg->gene1 ;
         keySet (kset, j++) = seg->gene2 ;
        }
     }
  keySetSort (kset) ;
  keySetCompress (kset) ;
  keySetNewDisplay (kset, "Ox") ;
}

static void hlExportM1 (void) 
{
  int i, j = 0 ;
  SEG *seg ;
  static KEYSET kset = 0 ;

  OXGET ("hlExportM1") ;

  kset = keySetReCreate (kset) ;

  for (i = 0 ; i < arrayMax (ox->segs) ; ++i)
    { seg = arrp (ox->segs, i, SEG) ;
      if (seg->box == ox->stayBox && !(seg->flag & FLAG_HIDE))
        keySet (kset, j++) = seg->gene1 ;
     }
  keySetSort (kset) ;
  keySetCompress (kset) ;
  keySetNewDisplay (kset, "Ox") ;
}

static void hlExportM2 (void) 
{
  int i, j = 0 ;
  SEG *seg ;
  static KEYSET kset = 0 ;

  OXGET ("hlExportM2") ;

  kset = keySetReCreate (kset) ;

  for (i = 0 ; i < arrayMax (ox->segs) ; ++i)
    { seg = arrp (ox->segs, i, SEG) ;
      if (seg->box == ox->stayBox && !(seg->flag & FLAG_HIDE))
        keySet (kset, j++) = seg->gene2 ;
     }
  keySetSort (kset) ;
  keySetCompress (kset) ;
  keySetNewDisplay (kset, "Ox") ;
}

    
  


 
 
 
