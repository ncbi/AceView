/*  File: pairmapdisp.c
 *  Written by Jo Dicks (jld@bioch.ox.ac.uk)
 *-------------------------------------------------------------------
 * This file is part of the ACeDB comparative mapping package, 
 * written by Jo Dicks and Michelle Kirby
 * The original ACeDB system was written by :
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *      Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *
 * Description:  to display single cells of Oxford Grids in detail
 * Exported functions: pairMapCreate() only
 * HISTORY:
 * Last edited: Feb 19 16:49 1997 (srk)
 *	- change NULL's to 0's to please MS Visual C++ compiler
 * * Dec  2 10:19 1996 (srk)
 * Created: December 1993 
 * New models : May 1995
 *-------------------------------------------------------------------
 */

/* $Id: pairmapdisp.c,v 1.2 2017/02/20 01:28:09 mieg Exp $ */

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

BOOL pairMapDisplay (void) ;
int pairMapConvert (OX ox) ;
static void pairMapPick (int box) ; 
static void pairMapSelect (OX ox, int box) ;
static void pairMapFollow (OX ox) ;    
static void pairMapDraw (void) ;
static void pairMapDestroy (void) ;
static void drawyChromosome (OX ox, float xlen) ;
static void drawxChromosome (OX ox, float ylen) ;
static void drawShade (OX ox) ;
static void drawSegs (OX ox) ;
static void pairMapFlip (void) ;
/* static void pairMapShowAll (void) ; */
static void toggleChroms (void) ;
static void toggleBars (void) ;
static void toggleShade (void) ;
static void drawyScale (void) ;
static void drawxScale (void) ;
static void noChroms (void) ;
static void pairClear (void) ;
static void pairLightAll (void) ;
static void pairLight (void) ;
static void pairUnLight (void) ;
static void pairHide (void) ;
static void pairHideLight (void) ;
static void pairUnHide (void) ;
static void pairExport (void) ;
static void pairLightHom (void) ;
static void pairUnLightHom (void) ;
static void pairFindHom (char *string) ;
static void pairFindLoc (char *string) ;
static int pairhomolCompletion (char *cp, int len) ;
static int pairlocusCompletion (char *cp, int len) ;
static void pairZoomIn (void) ;
static void pairZoomOut (void) ;

enum BoxNames { BACKGROUND=0,
                CHROM_BOX, BARS_BOX, FLIP_BOX, 
                HIDE_BOX, SHADE_BOX, 
                FINDHOM_BOX, FINDHOM2_BOX,
                FINDLOC_BOX, FINDLOC2_BOX,
                ZOOM_IN, ZOOM_OUT,
                MESS_BOX, MIN_BOX } ; /* MIN_BOX must be last */

static MENUOPT pairMapMenu[] =
            { {graphDestroy, "Quit"},
	      {help, "Help"},
	      {graphPrint, "Print"},
	      {displayPreserve, "Preserve"},
	       {0, 0}
            } ;

static MENUOPT HighlightMenu[] =
            { {pairClear, "Clear and show all"},
              {pairLightAll, "Highlight all"},
              {pairLightHom, "Highlight homology"},
              {pairUnLightHom, "Unhighlight homology"},
              {pairLight, "Highlight selected keyset"},
              {pairUnLight,"Unhighlight selected keyset"},
              {pairHide, "Restrict to keyset"},
              {pairHideLight, "Hide highlit items"},
              {pairUnHide, "Unhide hidden items"},
              {pairExport, "Export highlit keyset"},
	       {0, 0}
            } ;

static float xOff, yOff ;

/*********************************/

void pairMapCreate (void)
{ 
  int box ;
  BOX *bx ;

  OXGET ("pairMapCreate") ;

/* Find out if the user has selected a cell to look at */

  if (!(ox->gridBox))
    {messout ("You have not selected a cell") ;
     return ;}
  else
    {box = ox->gridBox ;
     if (box >= ox->mlb + ox->map_labels)
      {bx = arrayp (ox->boxes, (box - ox->mlb + 1), BOX) ;
       ox->c1 = bx->map1 ;
       ox->c2 = bx->map2 ;
       if (!pairMapDisplay ()) return ;
      }
    }

}


BOOL pairMapDisplay (void)
{
  int h ;

  OXGET ("pairMapDisplay") ;

  ox->c1min = 0 ; ox->c1max = 0 ; 
  ox->c2min = 0 ; ox->c2max = 0 ;

/* Set up the display as per usual */

/* If no points in the cell have an x AND a y position then give up */

  h = pairMapConvert (ox) ;

  if (!(h > 0)) 
    {messout ("Not enough data to draw a PairMap") ;
     return FALSE ;}
  
  if (graphExists (ox->pairGraph))
    {graphActivate (ox->pairGraph) ;
     graphPop() ;}
  else 
    {if (!displayCreate (PAIRMAP))
       return FALSE ;
     ox->pairGraph = graphActive () ;
    
     graphAssociate (&OX_MAGIC,ox) ;     
     graphRegister (DESTROY, pairMapDestroy) ; 
     graphRegister (MESSAGE_DESTROY, displayUnBlock) ;
     graphRegister (PICK, (GraphFunc) pairMapPick) ;
     graphMenu (pairMapMenu) ;
    }

  pairMapDraw () ;

  return TRUE ;


}




/***********************************************************/
/***************** Registered routines *********************/


/*************************************************/
/************* conversion and drawing ************/

static void pairMapDestroy (void)
{
  OXGET ("pairMapDestroy") ;
  
  arrayDestroy (ox->points) ;
  arrayDestroy (ox->c1segs) ;
  arrayDestroy (ox->c2segs) ;
  ox->c1 = 0 ;
  ox->c2 = 0 ;
  ox->s3 = 0 ;
  ox->s4 = 0 ;

}

int pairMapConvert (OX ox)
{  
  
  KEY band, col, key ;
  KEY map1, map2 ;
  float x1, x2, xmin = 0, xmax = 0 ;
  float ymin = 0, ymax = 0 ;
  float min1 = 0, max1 = 0 ;
  float min2 = 0, max2 = 0 ;
  float f1, f2, f3, f4 ;
  int h, i, box ;
  OBJ yChrom, yBand, xChrom, xBand ;
  OBJ loc1, loc2 ;
  SEG *tseg ;
  CHROMS *c1seg ;
  CHROMS *c2seg ;
  HOMOL *seg ;
  char *sym ;
  Array bands ;
  static BSMARK mark1 = 0 ;
  static BSMARK mark2 = 0 ;

/* Initialise the flags */

  ox->flag |= FLAG_PAIR_CHROMS ;
  ox->flag |= FLAG_SEGS ;
  ox->flag |= FLAG_YCEN ;
  ox->flag |= FLAG_XCEN ;
  ox->flag |= FLAG_BARS ;
  ox->flag |= FLAG_SHADE ;

  ox->pairwide = 60 ;
  ox->pairhigh = 40 ;

/* Check we can make the Map objects and if we can display the Maps */

  if (!(xChrom = bsCreate (ox->c1)))
    { return 0 ;}

  if(bsFindTag (xChrom, _Non_graphic))
    { ox->flag |= ~FLAG_PAIR_CHROMS ;}

  bsDestroy (xChrom) ;

 if (!(yChrom = bsCreate (ox->c2)))
    { return 0 ;}

  if(bsFindTag (yChrom, _Non_graphic))
    {ox->flag |= ~FLAG_PAIR_CHROMS ;}

  bsDestroy (yChrom) ;


/* Get band data on the first chromosome, if there is any */

  i = 0 ;
  bands = arrayCreate (8, BSunit) ;
  ox->c2segs = arrayCreate (128, CHROMS) ;
  ox->c1segs = arrayCreate (128, CHROMS) ;
  ox->points = arrayCreate (8, HOMOL) ;
 
  if ((xChrom = bsCreate (ox->c1)))
   {if (bsFindTag (xChrom, _Chrom_Band) && bsFlatten (xChrom, 1 , bands))
     {for (i = 0 ; i < arrayMax (bands) ; i++ )
       {band = array (bands, i, BSunit).k ; 
        xBand = bsCreate (band) ;
        if (bsFindTag (xBand, _Map) && bsGetKey (xBand, _bsRight, &key)) 
          if (bsPushObj (xBand))
            if (bsGetData (xBand, _Left, _Float, &x1) &&
                bsGetData (xBand, _Right, _Float, &x2))
              {c1seg = arrayp (ox->c1segs, i, CHROMS) ;
               c1seg->key = band ;
               c1seg->x = 0.5*(x1 + x2) ;
               c1seg->dx = x2 - x1 ;
	       c1seg->flag = 0 ;
               ox->flag |= FLAG_XBANDS ;}
        bsDestroy (xBand) ;
        if (x1 < xmin)
          xmin = x1 ;
        if (x2 > xmax)
	  xmax = x2 ;
        }

      for (i = 0 ; i < arrayMax (bands) ; i++ )
       {band = array (bands, i, BSunit).k ; 
        xBand = bsCreate (band) ;
        c1seg = arrayp (ox->c1segs, i, CHROMS) ;
        if (bsFindTag (xBand, _Centromere))
	  {c1seg->flag |= FLAG_CENTROMERE ;
           ox->flag ^= FLAG_XCEN ;}
        if (bsFindTag (xBand, _Dark))
	   c1seg->flag |= FLAG_DARK_BAND ;
	if (bsFindTag (xBand, _NOR))
	   c1seg->flag |= FLAG_NOR ;
	if (bsFindTag (xBand, _Colour) && bsGetKey (xBand, _bsRight, &col))
	   c1seg->flag |= ((col - _WHITE) & 15  ) << 28 ;
        if (bsFindTag (xBand, _Symbol)) 
          {bsGetData (xBand, _bsRight, _Text, &sym) ;
           c1seg->sym = sym ;
           c1seg->flag |= FLAG_SYM ;}
        bsDestroy (xBand) ;
       }
            
      ox->c1min = xmin ;
      ox->c1max = xmax ; 
      bsDestroy (xChrom) ;
     }
          
  }      

/* Get band data on the second chromosome, if there is any */  
 
  i = 0 ;
  if ((yChrom = bsCreate (ox->c2)))
   {if (bsFindTag (yChrom, _Chrom_Band) && bsFlatten (yChrom, 1 , bands))
     {for (i = 0 ; i < arrayMax (bands) ; i++ )
       {band = array (bands, i, BSunit).k ; 
        yBand = bsCreate (band) ;
        if (bsFindTag (yBand, _Map) &&
            bsGetKey (yBand, _bsRight, &key)) 
          if (bsPushObj (yBand))
            if (bsGetData (yBand, _Left, _Float, &x1) &&
               bsGetData (yBand, _Right, _Float, &x2))
             {c2seg = arrayp (ox->c2segs, i, CHROMS) ;
              c2seg->key = band ; 
              c2seg->x = 0.5*(x1+ x2) ;
              c2seg->dx = x2 - x1 ;
	      c2seg->flag = 0 ;
              ox->flag |= FLAG_YBANDS ;}
        bsDestroy (yBand) ;
        if (x1 < ymin)
          ymin = x1 ;
        if (x2 > ymax)
	  ymax = x2 ;
       }

      for (i = 0 ; i < arrayMax (bands) ; i++ )
       {band = array (bands, i, BSunit).k ; 
        yBand = bsCreate (band) ;
        c2seg = arrayp (ox->c2segs, i, CHROMS) ;
        if (bsFindTag (yBand, _Centromere))
	  {c2seg->flag |= FLAG_CENTROMERE ;
           ox->flag ^= FLAG_YCEN ;}
        if (bsFindTag (yBand, _Dark))
	   c2seg->flag |= FLAG_DARK_BAND ;
	if (bsFindTag (yBand, _NOR))
	   c2seg->flag |= FLAG_NOR ;
	if (bsFindTag (yBand, _Colour) &&
	   bsGetKey (yBand, _bsRight, &col))
	   c2seg->flag |= ((col - _WHITE) & 15  ) << 28;
        if (bsFindTag (yBand, _Symbol)) 
          {bsGetData (yBand, _bsRight, _Text, &sym) ;
           c2seg->sym = sym ;
           c2seg->flag |= FLAG_SYM ;}
        bsDestroy (yBand) ;
       }
   
      ox->c2min = ymin ;
      ox->c2max = ymax ;
      bsDestroy (yChrom) ; 
     }
          
   }  


/* Get homology data, Oxford Grid data plus positions on the two chromosomes */


  i = 0 ; h = 0 ;
  box = ox->gridBox ;

  for (i = 0 ; i < arrayMax (ox->segs) ; i++ )
    {tseg = arrayp (ox->segs, i, SEG) ;
     if (tseg->box == box)
      {seg = arrayp (ox->points, h, HOMOL) ;
       seg->hom = tseg->key ;
       seg->key1 = tseg->gene1 ;
       seg->key2 = tseg->gene2 ;
       seg->map1 = tseg->map1 ;
       seg->map2 = tseg->map2 ;
       seg->flag = 0 ;
       if (tseg->flag & FLAG_HIDE) seg->flag |= FLAG_HOMHIDE ;
       seg->x1 = 0 ; seg->dx1 = 0 ; 
       seg->x2 = 0 ; seg->dx2 = 0 ;

       if ((loc1 = bsCreate (seg->key1)))
         {if (bsGetKey (loc1, _Map, &map1)) do
           {mark1 = bsMark (loc1, mark1) ;
            if (chkSpecies (map1, seg->map1))
              {if (bsPushObj (loc1) &&bsGetData (loc1, _Position, _Float, &f1)) 
                {seg->x1 = f1 ;
                 seg->flag |= FLAG_XPOS ;
                 if (bsPushObj (loc1) && bsGetData (loc1, _Error, _Float, &f2))
                   seg->dx1 = f2 ;
                 else seg->dx1 = 0 ;
                 if ((seg->x1-seg->dx1) < min1) 
                   min1 = (seg->x1-seg->dx1) ;
                 if ((seg->x1+seg->dx1) > max1) 
                   max1 = (seg->x1+seg->dx1) ;
                 }
               }
             }while (bsGoto (loc1, mark1), bsGetKey (loc1, _bsDown, &map1)) ;
           bsDestroy (loc1) ;
          }

       if ((loc2 = bsCreate (seg->key2)))
         {if (bsGetKey (loc2, _Map, &map2)) do
           {mark2 = bsMark (loc2, mark2) ;
            if (chkSpecies (map2, seg->map2))
              {if (bsPushObj (loc2) &&bsGetData (loc2, _Position, _Float, &f3)) 
                {seg->x2 = f3 ;
                 seg->flag |= FLAG_YPOS ;
                 if (bsPushObj (loc2) && bsGetData (loc2, _Error, _Float, &f4))
                   seg->dx2 = f4 ;
                 else seg->dx2 = 0 ;
                 if ((seg->x2-seg->dx2) < min2) 
                   min2 = (seg->x2-seg->dx2) ;
                 if ((seg->x2+seg->dx2) > max2) 
                   max2 = (seg->x2+seg->dx2) ;
                 }
               }
             }while (bsGoto (loc2, mark2), bsGetKey (loc2, _bsDown, &map2)) ;
           bsDestroy (loc2) ;
          }


       if (seg->flag & FLAG_XPOS && seg->flag & FLAG_YPOS) h++ ;
      }
    }

  if (!(ox->flag & FLAG_XBANDS && ox->flag & FLAG_YBANDS))
    {ox->c1min = min1-5 ;
     ox->c1max = max1+5 ;
     ox->c2min = min2-5 ;
     ox->c2max = max2+5 ;
    }
     
  arrayDestroy (bands) ;
  return h ;
 

 }


/***************************************/


static void pairMapDraw (void)
{
  int box, i, xmax, ymax ;
  float xl, yl, xlen = 0, ylen = 0 ;
  CHROMS *c1seg, *c2seg ;

  OXGET ("pairMapDraw") ;

/* Draw all the nice buttons, some toggle, some don't */

  graphClear () ;

  if (ox->flag & FLAG_XBANDS && ox->flag & FLAG_YBANDS)
    {if (ox->flag & FLAG_PAIR_CHROMS)
      {if (graphButton ("Hide Chromosomes", toggleChroms, 2, 3) != CHROM_BOX)
       messcrash ("box screwup at %d in pairMapDraw", CHROM_BOX) ;}
    else
      {if (graphButton ("Draw Chromosomes", toggleChroms, 2, 3) != CHROM_BOX)
       messcrash ("box screwup at %d in pairMapDraw", CHROM_BOX) ;}
     }
  else
    {if (graphButton ("Draw Chromosomes", noChroms, 2, 3) != CHROM_BOX)
      messcrash ("box screwup at %d in pairMapDraw", CHROM_BOX) ;}

  if (ox->flag & FLAG_BARS)
      {if (graphButton ("Hide Error Bars", toggleBars, 20, 3) != BARS_BOX)
      messcrash ("box screwup at %d in pairMapDraw", BARS_BOX) ;}
    else
      {if (graphButton ("Draw Error Bars", toggleBars, 20, 3) != BARS_BOX)
      messcrash ("box screwup at %d in pairMapDraw", BARS_BOX) ;}

  if (graphButton ("Flip", pairMapFlip, 37, 3) != FLIP_BOX)
    messcrash ("box screwup at %d in pairMapDraw", FLIP_BOX) ;

  if (graphButton ("Highlight...", pairClear, 43, 3) != HIDE_BOX)
    messcrash ("box screwup at %d in pairMapDraw", HIDE_BOX) ;

  graphBoxMenu (HIDE_BOX, HighlightMenu) ;
  
  if (ox->flag & FLAG_SHADE)
      {if (graphButton (" Shade ", toggleShade, 57, 3) != SHADE_BOX)
      messcrash ("box screwup at %d in pairMapDraw", SHADE_BOX) ;}
    else
      {if (graphButton ("Unshade", toggleShade, 57, 3) != SHADE_BOX)
      messcrash ("box screwup at %d in pairMapDraw", SHADE_BOX) ;}

  graphText ("Find homology :", 2, 5) ;
  *ox->pairhom = 0 ;
  if (graphCompScrollEntry (pairhomolCompletion, ox->pairhom, 20, 20, 18, 5, pairFindHom) != FINDHOM_BOX)
    messcrash ("box screwup at %d in pairMapDraw", FINDHOM_BOX) ;

  graphEntryDisable () ;

  graphText ("Find locus :", 43, 5) ;
  *ox->findloc = 0 ;
  if (graphCompScrollEntry (pairlocusCompletion, ox->pairloc, 20, 20, 56, 5, pairFindLoc) != FINDLOC_BOX)
    messcrash ("box screwup at %d in pairMapDraw", FINDLOC_BOX) ;

  graphEntryDisable () ;

  if (graphButton ("Zoom in", pairZoomIn, 66, 3) != ZOOM_IN)
    messcrash ("box screwup at %d in pairMapDraw", ZOOM_IN) ;

  if (graphButton ("Zoom out", pairZoomOut, 75, 3) != ZOOM_OUT)
    messcrash ("box screwup at %d in pairMapDraw", ZOOM_OUT) ;

/* Draw the header box to name a selected homology on screen */
  
  ox->messageBox = MESS_BOX ;
 *ox->messageText = 0 ;
  ox->messageBox = graphBoxStart () ;
  graphTextPtr (ox->messageText, 2, 7, 80) ;
  graphBoxEnd () ;
  graphBoxDraw (ox->messageBox, BLACK, LIGHTBLUE) ;
      
/* Now draw the titles */
 
  xOff = 0 ; yOff = 0 ;

  graphText (messprintf ("PairMap for map %s and map %s", name (ox->c1), name (ox->c2)), 2, 1); 
  graphLine (0, 9, 100, 9) ;

  xOff += 2 + strlen (name (ox->c2)) ;
  yOff += 10 ;
  
/* If we have bands for both chromsomes we shall draw them, otherwise scales */

  if (ox->flag & FLAG_XBANDS && ox->flag & FLAG_YBANDS)
   {for (i = 0 ; i < arrayMax (ox->c2segs) ; i++)
      {c2seg = arrayp (ox->c2segs, i, CHROMS) ;
       xl = strlen (name (c2seg->key)) ;
       if (xl > xlen) xlen = xl ;
       }
    xOff += 8 + xlen ;
    for (i = 0 ; i < arrayMax (ox->c1segs) ; i++)
      {c1seg = arrayp (ox->c1segs, i, CHROMS) ;
       yl = strlen (name (c1seg->key)) ;
       if (yl > ylen) ylen = yl ;
      }
    ylen *= 0.65 ;
    yOff += 5 + ylen ;
    }
   else
    {xOff += 10 ;
     yOff += 5  ;
    } 

  xmax = arrayMax (ox->c1segs) ;
  ymax = arrayMax (ox->c2segs) ;

  box = MIN_BOX ;
  if (box != graphBoxStart ())
    messcrash ("minLiveBox wrong in pairMapDraw () box %d", box) ;
  graphTextFormat (BOLD) ;
  graphText (name (ox->c1), xOff + (ox->pairwide)/2, 10) ;
  graphBoxEnd () ;
  graphBoxDraw (box, BLACK, WHITE) ;
  box++ ;
  if (box != graphBoxStart ())
    messcrash ("minLiveBox wrong in pairMapDraw () box %d", box) ;
  graphTextFormat (BOLD) ;
  graphText (name (ox->c2), 2, yOff + (ox->pairhigh)/2) ;
  graphBoxEnd () ;
  graphBoxDraw (box, BLACK, WHITE) ;
  graphTextFormat(PLAIN_FORMAT);
  
  if (ox->flag & FLAG_XBANDS && ox->flag & FLAG_YBANDS)
    {if (ox->flag & FLAG_PAIR_CHROMS)
      {ox->pair_labels = 2 + xmax + ymax ;
       drawyChromosome (ox, xlen) ;
       drawxChromosome (ox, ylen) ;
       }
     else
      {ox->pair_labels = 2 ;
       drawyScale () ;
       drawxScale () ;
       }
     }
   else
    {ox->pair_labels = 2 ;
     drawyScale () ;
     drawxScale () ;
     }

  
  if (!(ox->flag & FLAG_SHADE)) drawShade (ox) ;
     
  
/* And draw the box */

  graphLine (xOff, yOff, xOff+ox->pairwide, yOff) ;
  graphLine (xOff+ox->pairwide, yOff, xOff+ox->pairwide, yOff+ox->pairhigh) ;
  graphLine (xOff+ox->pairwide, yOff+ox->pairhigh, xOff, yOff+ox->pairhigh) ;
  graphLine (xOff, yOff+ox->pairhigh, xOff, yOff) ;

  drawSegs (ox) ;
  
  ox->pairBox = 0 ;

  graphRedraw () ;
}

/*************************************************************/


static void drawyChromosome (OX ox, float xlen)
{
  float y0, dy, y, yend ;
  int i, box ;
  CHROMS *c2seg ;

  y0 = -(ox->c2min) ;
  yend = (ox->c2max - ox->c2min) / ox->pairhigh ;

/* Draw the y chromosome alongside the box */

  box = 2 + MIN_BOX ;
  for (i = 0 ; i < arrayMax (ox->c2segs) ; ++i)
    {c2seg = arrayp (ox->c2segs, i, CHROMS) ;
     if (class (c2seg->key) == _VChrom_Band)
       {
	y = yOff + ( y0 + c2seg->x) / yend ;
	dy = ( 0.5 * c2seg->dx) / yend ;
	if (c2seg->flag & FLAG_CENTROMERE)
	  {graphColor (WHITE) ;
	   graphFillRectangle (xOff-4, y-dy, xOff, y+dy) ;
	   graphColor (BLACK) ;
	   graphLine (xOff-3, y-dy, xOff-1, y+dy) ;
	   graphLine (xOff-1, y-dy, xOff-3, y+dy) ;
	   graphLine (xOff-3, y-dy, xOff-1, y-dy) ;
	   graphLine (xOff-1, y+dy, xOff-3, y+dy) ;
           if (box != graphBoxStart ())
           messcrash ("minLiveBox wrong in pairMapDraw () box %d", box) ;
           graphTextHeight (0.5) ;
           graphText (name (c2seg->key), xOff-6-xlen, y-0.25) ;
           graphBoxEnd () ;
           graphBoxDraw (box, BLUE, WHITE) ;
           box++ ; 
           graphTextHeight (1) ;
	   }	      
	else 
	  {if (c2seg->flag & FLAG_DARK_BAND)
	     graphColor (DARKGRAY) ;
	   else if (c2seg->flag & FLAG_NOR)
	     graphColor (LIGHTGRAY) ;
	   else
	     graphColor (WHITE) ;
	   if (c2seg->flag >> 28)
	     graphColor((c2seg->flag  >> 28) + WHITE) ;
	   graphFillRectangle (xOff-4, y-dy, xOff, y+dy) ;
	   graphColor (BLACK) ;
	   graphRectangle (xOff-4, y-dy, xOff, y+dy) ;
           if (box != graphBoxStart ())
           messcrash ("minLiveBox wrong in pairMapDraw () box %d", box) ;
           graphTextHeight (0.5) ;
           graphText (name (c2seg->key), xOff-6-xlen, y-0.25) ;
           graphBoxEnd () ;
           graphBoxDraw (box, BLUE, WHITE) ;
           box++ ; 
           graphTextHeight (1) ;
	   }
        }
     }
}



  static void drawxChromosome (OX ox, float ylen)
{
  float x0, dx, x, xend ;
  int i, box ;
  CHROMS *c1seg ;

  x0 = -(ox->c1min) ;
  xend = (ox->c1max - ox->c1min) / ox->pairwide ;

/* Draw the x chromosome alongside the box */

  box = 2 + MIN_BOX + arrayMax (ox->c2segs) ;
  for (i = 0 ; i < arrayMax (ox->c1segs) ; ++i)
    {c1seg = arrayp (ox->c1segs, i, CHROMS) ;
     if (class(c1seg->key) == _VChrom_Band)
       {
       x = xOff + ox->pairwide - ( x0 + c1seg->x) / xend ; /* bacwards ?? */
       /* x = xOff + ( x0 + c1seg->x) / xend ; */
	dx = ( 0.5 * c1seg->dx) / xend ;
	if (c1seg->flag & FLAG_CENTROMERE)
	  {graphColor (WHITE) ;
	   graphFillRectangle (x-dx, yOff-2.4, x+dx, yOff) ;
	   graphColor (BLACK) ;
	   graphLine (x-dx, yOff-1.4, x+dx, yOff-1) ;
	   graphLine (x-dx, yOff-1, x+dx, yOff-1.4) ;
	   graphLine (x-dx, yOff-1.4, x-dx, yOff-1) ;
	   graphLine (x+dx, yOff-1, x+dx, yOff-1.4) ;
           if (box != graphBoxStart ())
           messcrash ("minLiveBox wrong in pairMapDraw () box %d", box) ;
           graphTextHeight (0.5) ;
           japanese (name (c1seg->key), x-0.2, yOff-3-ylen) ; 
           graphBoxEnd () ;
           graphBoxDraw (box, BLUE, WHITE) ;
           box++ ; 
           graphTextHeight (1) ;
           }	      
        else 
	  {if (c1seg->flag & FLAG_DARK_BAND)
	     graphColor (DARKGRAY) ;
	   else if (c1seg->flag & FLAG_NOR)
	     graphColor (LIGHTGRAY) ;
	   else
             graphColor (WHITE) ;
	   if (c1seg->flag >> 28)
	     graphColor ((c1seg->flag  >> 28) + WHITE) ;
	   graphFillRectangle (x-dx, yOff-2.4, x+dx, yOff) ;
	   graphColor (BLACK) ;
	   graphRectangle (x-dx, yOff-2.4, x+dx, yOff) ;
           if (box != graphBoxStart ())
           messcrash ("minLiveBox wrong in pairMapDraw () box %d", box) ;
           graphTextHeight (0.5) ;
           japanese (name (c1seg->key), x-0.2, yOff-3-ylen) ; 
           graphBoxEnd () ; 
           graphBoxDraw (box, BLUE, WHITE) ;
           box++ ; 
           graphTextHeight (1) ;
     	   }
        }
     } 
}

static void drawyScale (void) 
{
  int i, j, k ;
  float unit , y0 ;

  OXGET ("drawyScale") ;

/* If no banding data or if the chromosomes are hidden, draw a scale instead */

  graphTextHeight (0.5) ;
  y0 = yOff + ox->pairhigh * (1-(ox->c2max / (ox->c2max-ox->c2min))) ;
  unit = ox->pairhigh / (ox->c2max - ox->c2min) ;
  for (i = 0 ; i < (ox->c2max + 1) ; i++)
    graphLine (xOff-3, y0+i*unit, xOff-3.9, y0+i*unit) ;
  for (k = 0 ; k < (ox->c2max + 1) ; k+=10)
   {graphLine (xOff-3, y0+k*unit, xOff-4.5, y0+k*unit) ;
    graphText (messprintf ("%3d",k), xOff - 8, y0+k*unit-0.5) ;
    }
  for (j = 0 ; j < -(ox->c2min - 1) ; j++)
    graphLine (xOff-3, y0-j*unit, xOff-3.9, y0-j*unit) ;
  for (k = 0 ; k < -(ox->c2min - 1) ; k+=10)
   {graphLine (xOff-3, y0-k*unit, xOff-4.5, y0-k*unit) ;
    graphText (messprintf ("%3d",-k), xOff - 8, y0-k*unit-0.5) ;
    }  
  graphLine (xOff-3, y0-(j-1)*unit, xOff-3, y0+(i-1)*unit) ;
  graphTextHeight (1) ;
    

}

static void drawxScale (void) 
{
  int i, j, k ;
  float unit, x0 ;

  OXGET ("drawxScale") ;

/* If no banding data or if the chromosomes are hidden, draw a scale instead */

  graphTextHeight (0.5) ;
  x0 = xOff + ox->pairwide * (ox->c1max / (ox->c1max - ox->c1min)) ;
  unit = ox->pairwide / (ox->c1max - ox->c1min) ;
  for (i = 0 ; i < (ox->c1max + 1) ; i++)
    graphLine (x0 - i*unit, yOff - 1, x0 - i*unit, yOff - 1.6) ;
  for (k = 0 ; k < (ox->c1max + 1) ; k += 10)
   {graphLine (x0 - k*unit, yOff - 1, x0 - k*unit, yOff - 2) ;
    graphText (messprintf ("%3d", k), x0 - k*unit - 2, yOff - 3) ;
    }
  for (j = 0 ; j < -(ox->c1min - 1) ; j++)
    graphLine (x0 + j*unit, yOff - 1, x0 + j*unit, yOff - 1.6) ;
  for (k = 0 ; k < -(ox->c1min - 1) ; k += 10)
   {graphLine (x0 + k*unit, yOff - 1, x0 + k*unit, yOff - 2) ;
    graphText (messprintf ("%3d", -k), x0 + k*unit - 2, yOff - 3) ;
    } 
  graphLine (x0 - (i - 1)*unit, yOff - 1, x0 + (j - 1)*unit, yOff - 1) ;
  graphTextHeight (1) ;

}

static void drawShade (OX ox)
{
  int i, j ;
  float x0, dx, x, xend ;
  float y0, dy, y, yend ;
  CHROMS *c1seg, *c2seg ;


  x0 = -(ox->c1min) ;
  xend = (ox->c1max - ox->c1min) / ox->pairwide ;
  y0 = -(ox->c2min) ;
  yend = (ox->c2max - ox->c2min) / ox->pairhigh ;

  for (i = 0 ; i < arrayMax (ox->c2segs) ; ++i)
    {c2seg = arrayp (ox->c2segs, i, CHROMS) ;
     if (class (c2seg->key) == _VChrom_Band)
       {
        y = yOff + ( y0 + c2seg->x) / yend ;
	dy = ( 0.5 * c2seg->dx) / yend ;
        for (j = 0 ; j < arrayMax (ox->c1segs) ; ++j)
         {c1seg = arrayp (ox->c1segs, j, CHROMS) ;
          if (class (c1seg->key) == _VChrom_Band)
            {
             x = xOff + ox->pairwide - ( x0 + c1seg->x) / xend ;
	     dx = ( 0.5 * c1seg->dx) / xend ;
             if ((c1seg->flag & FLAG_DARK_BAND) && (c2seg->flag & FLAG_DARK_BAND)) 
               {graphColor (DARKGRAY) ;
                graphFillRectangle (x-dx,y-dy,x+dx,y+dy) ;
                graphColor (BLACK) ;}
             else if (((c1seg->flag & FLAG_DARK_BAND) && !(c2seg->flag & FLAG_DARK_BAND)) || (!(c1seg->flag & FLAG_DARK_BAND) && (c2seg->flag & FLAG_DARK_BAND)))
               {graphColor (LIGHTGRAY) ;
                graphFillRectangle (x-dx,y-dy,x+dx,y+dy) ;
                graphColor (BLACK) ;}
             }
           }
         }
       }
   graphLine (xOff, yOff, xOff+ox->pairwide, yOff) ;
   graphLine (xOff, yOff, xOff, yOff+ox->pairhigh) ;
   

}

static void drawSegs (OX ox)
{
  int i, p ;
  float x, dx, y, dy ;
  float x0, xend, y0, yend ;
  HOMOL *seg ;
 
  y0 = -(ox->c2min) ;
  yend = (ox->c2max - ox->c2min) / ox->pairhigh ;
  x0 = -(ox->c1min) ;
  xend = (ox->c1max - ox->c1min) / ox->pairwide ;

/* Draw the homologies, with error bars and with Homology_group names */

  p = MIN_BOX + ox->pair_labels ;
  for (i = 0 ; i < arrayMax(ox->points) ; i++)
    {seg = arrayp (ox->points, i, HOMOL) ;
     if (seg->flag & FLAG_XPOS && seg->flag & FLAG_YPOS)
       if (!(seg->flag & FLAG_HOMHIDE))
         {x = xOff + ox->pairwide - (x0 + seg->x1) / xend ; dx = (seg->dx1) / xend ;
          y = yOff + (y0 + seg->x2) / yend ; dy = (seg->dx2) / yend ;
          if (p != graphBoxStart ())
            messcrash ("problems with pairMap point %d", p) ;
          graphRectangle (x-0.4, y-0.4, x+0.4, y+0.4) ;
          graphBoxEnd () ;
          if (seg->flag & FLAG_HIGHLIGHT) graphBoxDraw (p, WHITE, MAGENTA) ;
          else graphBoxDraw (p, WHITE, BLUE) ;
          seg->pbox = p ;
          p++ ;
          if (ox->flag & FLAG_BARS)
            {graphLine (x, y-dy, x, y+dy) ;
             graphLine (x-dx, y, x+dx, y) ;
             graphLine (x-0.25, y-dy, x+0.25, y-dy) ;
             graphLine (x-0.25, y+dy, x+0.25, y+dy) ;
             graphLine (x-dx, y-0.25, x-dx, y+0.25) ;
             graphLine (x+dx, y-0.25, x+dx, y+0.25) ;
             }
          }
     }
}

extern int ksetClassComplete (char *text, int len, int classe) ;

static int pairhomolCompletion (char *cp, int len)
{
  return ksetClassComplete (cp, len, _VHomology_group) ;
}

static int pairlocusCompletion (char *cp, int len)
{
  return ksetClassComplete (cp, len, _VLocus) ;
}

static void pairFindHom (char* string)
{
  int i ;
  KEY key ;
  HOMOL *seg ;
  BOOL found = FALSE ;

  OXGET ("pairFindHom") ;

  if (!(lexword2key (string, &key, _VHomology_group)))
    {messout ("Sorry, %s is not a homology name", string) ;
     return ;
    }

  for (i = 0 ; i < arrayMax (ox->points) ; i++)
    {seg = arrayp (ox->points, i, HOMOL) ;
     if (chkSpecies (seg->hom, key)) 
       {found = TRUE ;
        pairMapSelect (ox, seg->pbox) ;
       }
     }

  if (found == FALSE) 
   {messout ("Sorry, this homology is not in this grid") ;
    return ;
   }

}

static void pairFindLoc (char* string)
{
  int i ;
  KEY key ;
  HOMOL *seg ;
  BOOL found = FALSE ;

  OXGET ("pairFindLoc") ;

  if (!(lexword2key (string, &key, _VLocus)))
    {messout ("Sorry, %s is not a locus name", string) ;
     return ;
    }

  for (i = 0 ; i < arrayMax (ox->points) ; i++)
    {seg = arrayp (ox->points, i, HOMOL) ;
     if (chkSpecies (seg->key1, key) || chkSpecies (seg->key2, key)) 
       {found = TRUE ;
        pairMapSelect (ox, seg->pbox) ;
       }
     }

  if (found == FALSE) 
   {messout ("Sorry, this locus is not in this grid") ;
    return ;
   }

}

static void pairZoomOut (void) 
{
  OXGET ("pairZoomOut") ;

  if (ox->pairwide >= 30 && ox->pairhigh >= 20)
    {ox->pairwide -= 15 ;
     ox->pairhigh -= 10 ;

     pairMapDraw () ;
    }

}

static void pairZoomIn (void) 
{
  OXGET ("pairZoomIn") ;

  if (ox->pairwide <= 120 && ox->pairhigh <= 80)
    {ox->pairwide += 15 ;
     ox->pairhigh += 10 ;

     pairMapDraw () ;
    }


}


static void pairMapPick (int box) 
{ 
  int ymax ;

  OXGET ("pairMapPick") ;

  /* int  xmax = arrayMax (ox->c1segs) ; */
  ymax = arrayMax (ox->c2segs) ;


  /* The usual - check if the box number is legal, select it if 
     a single click and follow it if a double */

  if (box >= ox->pair_labels + MIN_BOX && 
      box <= (arrayMax (ox->points) + ox->pair_labels + MIN_BOX))/*a seg*/
   {if (box == ox->pairBox)
      pairMapFollow (ox) ;
    graphActivate (ox->pairGraph) ;
    pairMapSelect (ox, box) ; 
   }
  else if (box == MIN_BOX)
    {if (box == ox->pairBox)
       display (ox->c1, 0, 0) ;
     graphActivate (ox->pairGraph) ;
     pairMapSelect (ox, box) ;
    }
  else if (box == MIN_BOX + 1)
    {if (box == ox->pairBox)
       display (ox->c2, 0, 0) ;
     graphActivate (ox->pairGraph) ;
     pairMapSelect (ox, box) ;
    }
  else if (box >= MIN_BOX + 2 && box < MIN_BOX + 2 + ymax)
    {if (box == ox->pairBox)
       display (arrayp (ox->c2segs, box - MIN_BOX - 2, CHROMS)->key, 0, 0) ;
     graphActivate (ox->pairGraph) ;
     pairMapSelect (ox, box) ;
    }
  else if (box >= MIN_BOX + 2 + ymax && box < ox->pair_labels + MIN_BOX)
    {if (box == ox->pairBox)
       display (arrayp (ox->c1segs, box-MIN_BOX-2-ymax, CHROMS)->key, 0, 0) ;
     graphActivate (ox->pairGraph) ;
     pairMapSelect (ox, box) ;
    }
  else if (box == FINDHOM_BOX)
    graphCompScrollEntry (0, ox->pairhom, 0, 0, 0, 0, 0) ;
  else if (box == FINDLOC_BOX)
    graphCompScrollEntry (0, ox->pairloc, 0, 0, 0, 0, 0) ;
        
}

static void pairMapSelect (OX ox, int box)
{
  int i ;
  KEY hit=0, key1=0, key2=0 ;
  HOMOL *seg ;
  BOOL HID = FALSE ;

  /* The usual - return the old box to normal and colour in the new */

  if (ox->pairBox == box) return ;

  for (i = 0 ; i < arrayMax (ox->points) ; ++i)
    {seg = arrp (ox->points, i, HOMOL) ;
     if (seg->pbox == ox->pairBox)
       if (seg->flag & FLAG_HIGHLIGHT) HID = TRUE ;
    }
      
  if (ox->pairBox) 
   {if (ox->pairBox >= ox->pair_labels + MIN_BOX)
     {if (HID == TRUE) graphBoxDraw (ox->pairBox, WHITE, MAGENTA) ;
      else graphBoxDraw (ox->pairBox, WHITE, BLUE) ;
      }
    else if (ox->pairBox == MIN_BOX || ox->pairBox == MIN_BOX + 1)
      graphBoxDraw (ox->pairBox, BLACK, WHITE) ;
    else if (ox->pairBox >= MIN_BOX + 2 && ox->pairBox < MIN_BOX + ox->pair_labels)
      graphBoxDraw (ox->pairBox, BLUE, WHITE) ;
    }
    

  ox->pairBox = box ;

  if (ox->pairBox >= ox->pair_labels + MIN_BOX)
    graphBoxDraw (ox->pairBox, WHITE, GREEN) ;
  else if (ox->pairBox == MIN_BOX || ox->pairBox == MIN_BOX + 1)
    graphBoxDraw (ox->pairBox, WHITE, BLACK) ;
  else if (ox->pairBox >= MIN_BOX+2 && ox->pairBox < MIN_BOX + ox->pair_labels)
    graphBoxDraw (ox->pairBox, WHITE, BLUE) ;

  graphCompScrollEntry (0, ox->pairhom, 0, 0, 0, 0, 0) ; 
  graphEntryDisable () ;
  graphCompScrollEntry (0, ox->pairloc, 0, 0, 0, 0, 0) ;
  graphEntryDisable () ;

  if (ox->messageBox)
    {
      if (ox->pairBox >= ox->pair_labels + MIN_BOX)
	{
	  for (i = 0 ; i < arrayMax (ox->points) ; i++)
	    {
	      seg = arrayp (ox->points, i, HOMOL) ;
	      if (seg->pbox == ox->pairBox)
		{hit = seg->hom ;
		key1 = seg->key1 ;
		key2 = seg->key2 ;}
	    }
	  if (hit)
	    {
	      strncpy (ox->messageText, messprintf ("Homology "), 80) ; 
	      strcat (ox->messageText, messprintf ("%s ", name (hit))) ;
	      strcat (ox->messageText, messprintf ("between %s and %s", name (key1), name (key2))) ;
	    }
	}
    }
  else
    *ox->messageText = 0 ;
 

  graphBoxDraw (ox->messageBox, BLACK, LIGHTBLUE) ;

}

static void pairMapFollow (OX ox)
{
  int i ;
  KEY hit = 0 ;
  HOMOL *seg ;

  /* Similar to usual - seg->pbox needed to make the keyset option work */

  for (i=0 ; i < arrayMax (ox->points) ; i++)
    {seg = arrayp (ox->points, i, HOMOL) ;
     if (seg->pbox == ox->pairBox)
       hit = seg->hom ;
    }
  if (hit)
    display (hit, 0, 0) ;
}


static void noChroms (void)
{
  OXGET ("noChroms") ;

/* Inactivate the Hide chromosomes button if there isn't any band data */

  messout ("Sorry, there's no band data") ;
  return ;

}

static void toggleChroms (void)
{
  OXGET ("toggleChroms") ;

/* Toggle for the Hide chromosomes button */

  if (!(ox->flag & FLAG_XBANDS && ox->flag & FLAG_YBANDS))
    {messout ("Sorry, there's no band data") ;
     return ;}
  else
    {ox->flag ^= FLAG_PAIR_CHROMS ;
     pairMapDraw () ;}

}

static void toggleBars (void)
{ 
 OXGET("toggleBars") ;

/* Toggle for the Hide error bars button */

 ox->flag ^= FLAG_BARS ;
 pairMapDraw () ;

}

static void toggleShade (void)
{
  OXGET ("toggleShade") ;

/* Toggle for the Shade button */

  if (!(ox->flag & FLAG_XBANDS && ox->flag & FLAG_YBANDS))
    {messout ("Sorry, there's no band data") ;
     return ;}
  else
    {ox->flag ^= FLAG_SHADE ;
     pairMapDraw () ;}

}

#ifdef JUNK
static void pairMapShowAll (void)
{
  int i ;
  HOMOL *seg ;
 
  OXGET ("pairMapShowAll") ;

/* Show all homologies, whether in or out of the keyset */

  seg = arrayp (ox->points,0,HOMOL) ;
  i = arrayMax (ox->points) ;
  while (i--)
    {seg->flag &= ~FLAG_HOMHIDE ;  /* bitwise NOT */
     seg++ ;
    }

  pairMapDraw () ;
}
#endif

static void pairMapFlip (void)
{
  int i ;
  float dx, x, min, max ;
  KEY chrom, key, spec ;
  HOMOL *seg ;
  Array temp ;

  OXGET ("pairMapFlip") ;

  ox->flag ^= FLAG_PAIR_FLIP ;

/* Swap all x and y data so that the box will flip over when redrawn */

  spec = ox->s3 ;
  ox->s3 = ox->s4 ;
  ox->s4 = spec ;

  chrom = ox->c1 ;
  ox->c1 = ox->c2 ;
  ox->c2 = chrom ;

  temp = arrayCopy (ox->c1segs) ;
  ox->c1segs = arrayCopy (ox->c2segs) ;
  ox->c2segs = arrayCopy (temp) ;
  arrayDestroy (temp) ;

  min = ox->c1min ;
  ox->c1min = ox->c2min ;
  ox->c2min = min ;

  max = ox->c1max ;
  ox->c1max = ox->c2max ;
  ox->c2max = max ;

/*  if (!(ox->flag & FLAG_YCEN)) cen = TRUE ;
  else cen = FALSE ;
  if (!(ox->flag & FLAG_XCEN)) ox->flag |= FLAG_YCEN ;
  else ox->flag |= ~FLAG_YCEN ;
  if (cen) ox->flag |= FLAG_XCEN ;
  else ox->flag |= ~FLAG_XCEN ; */

  for (i = 0 ; i < arrayMax (ox->points) ; i++)
   {seg = arrayp (ox->points, i, HOMOL) ;
    key = seg->key1 ; seg->key1 = seg->key2 ; seg->key2 = key ;
    x = seg->x1 ; seg->x1 = seg->x2 ; seg->x2 = x ;
    dx = seg->dx1 ; seg->dx1 = seg->dx2 ; seg->dx2 = dx ;}

  pairMapDraw () ; 

}

static void pairClear (void) 
{
  int i ;
  HOMOL *seg ;
 
  OXGET ("pairClear") ;

  seg = arrayp (ox->points, 0, HOMOL) ;
  i = arrayMax (ox->points) ;
  while(i--)
    {seg->flag &= ~FLAG_HOMHIDE ;  /* bitwise NOT */
     seg->flag &= ~FLAG_HIGHLIGHT ;
     seg++ ;
    }

  pairMapDraw () ;

}

static void pairLightAll (void) 
{
  int i ;

  OXGET ("pairLightAll") ;

  for (i = 0 ; i < arrayMax (ox->points) ; ++i)
    arrp (ox->points, i, HOMOL)->flag |= FLAG_HIGHLIGHT ;
    
  ox->pairBox = 0 ;

  pairMapDraw () ;

}

static void pairLight (void) 
{
  int i, j ;
  HOMOL *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;

  OXGET ("pairLight") ;

  if (!keySetActive (&keySet, &dummy))
    { messout ("First select a keySet window, thank you.") ;
      return ;
    }
  for (i = 0 ; i < arrayMax (ox->points) ; ++i)
    { seg = arrp (ox->points, i, HOMOL) ;
      if (keySetFind (keySet, seg->hom, &j) ||
          keySetFind (keySet, seg->key1, &j) ||
          keySetFind (keySet, seg->key2, &j)) 
	seg->flag |= FLAG_HIGHLIGHT ;
    }

  pairMapDraw () ;

}

static void pairUnLight (void) 
{
  int i, j ;
  HOMOL *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;

  OXGET ("pairUnLight") ;

  if (!keySetActive (&keySet, &dummy))
    { messout ("First select a keySet window, thank you.") ;
      return ;
    }
  for (i = 0 ; i < arrayMax (ox->points) ; ++i)
    { seg = arrp (ox->points, i, HOMOL) ;
      if (keySetFind (keySet, seg->hom, &j) ||
          keySetFind (keySet, seg->key1, &j) ||
          keySetFind (keySet, seg->key2, &j)) 
	seg->flag &= ~FLAG_HIGHLIGHT ;
    }

  pairMapDraw () ;

}

static void pairLightHom (void) 
{
  int i ;
  HOMOL *seg ;

  OXGET ("pairLightHom") ;

  if (!(ox->pairBox))
    {messout ("You have not selected a homology") ;
     return ;}

  for (i = 0 ; i < arrayMax (ox->points) ; i++)
    {seg = arrayp (ox->points, i, HOMOL) ;
     if (seg->pbox == ox->pairBox) seg->flag |= FLAG_HIGHLIGHT ;
    }

  pairMapDraw () ;

}

static void pairUnLightHom (void) 
{
  int i ;
  HOMOL *seg ;

  OXGET ("pairUnLightHom") ;

  if (!(ox->pairBox))
    {messout ("You have not selected a homology") ;
     return ;}

  for (i = 0 ; i < arrayMax (ox->points) ; i++)
    {seg = arrayp (ox->points, i, HOMOL) ;
     if (seg->pbox == ox->pairBox) seg->flag &= ~FLAG_HIGHLIGHT ;
    }

  pairMapDraw () ;

}

static void pairHide (void)
{
  int i, j ;
  HOMOL *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;
 
  OXGET ("pairHide") ;

  if (!keySetActive (&keySet, &dummy))
    {messout ("First select a keySet window, thank you.") ;
     return ;
    }

  seg = arrp (ox->points, 0, HOMOL) ;
  i = arrayMax (ox->points) ;
  while(i--)
    {seg->flag &= ~FLAG_HOMHIDE ; 
     if (!keySetFind (keySet, seg->hom, &j) &&
         !keySetFind (keySet, seg->key1, &j) &&
         !keySetFind (keySet, seg->key2, &j)) 
       seg->flag |= FLAG_HOMHIDE ;
     seg++ ;
    }

  pairMapDraw () ;

}

static void pairHideLight (void) 
{
  int i ;
  HOMOL *seg ;

  OXGET ("pairHideLight") ;

  for (i = 0 ; i < arrayMax (ox->points) ; ++i)
    { seg = arrp (ox->points, i, HOMOL) ;
      if (seg->flag & FLAG_HIGHLIGHT)
	seg->flag |= FLAG_HOMHIDE ;
    }
  pairMapDraw () ;

}

static void pairUnHide (void) 
{
  int i ;
  HOMOL *seg ;

  OXGET ("pairUnHide") ;

  for (i = 0 ; i < arrayMax (ox->points) ; ++i)
    { seg = arrp (ox->points, i, HOMOL) ;
      seg->flag &= ~FLAG_HOMHIDE ;
    }
  pairMapDraw () ;


}

static void pairExport (void) 
{
  int i, j, k = 0 ;
  HOMOL *seg ;
  static KEYSET kset = 0 ;
  BOOL found = FALSE ;

  OXGET ("pairExport") ;

  kset = keySetReCreate (kset) ;

  for (i = 0 ; i < arrayMax (ox->points) ; ++i)
    { found = FALSE ;
      seg = arrp (ox->points, i, HOMOL) ;
      if (seg->flag & FLAG_HIGHLIGHT)
       {for (j = 0 ; j < keySetMax (kset) ; j++)
         {if (chkSpecies (seg->hom, keySet (kset, j))) found = TRUE ;}
	if (found == FALSE) keySet (kset, k++) = seg->hom ;
       }
    }
  keySetSort (kset) ;
  keySetCompress (kset) ;
  keySetNewDisplay (kset, "Ox") ;
}









/************** end of file ****************/
 
 
 
 
 
 
