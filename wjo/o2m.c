/*  File: o2m.c
 *  Written by Jo Dicks (jld@bioch.ox.ac.uk)
 *-------------------------------------------------------------------
 * This file is part of the ACeDB comparative mapping package, 
 * written by Jo Dicks and Michelle Kirby
 * The original ACeDB system was written by :
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *      Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *
 * Description:  to display One-2-Many chromosome maps
 * Exported functions: o2mDisplay () only
 * HISTORY:
 * Last edited: Feb 19 16:49 1997 (srk)
 *	- change NULL's to 0's to please MS Visual C++ compiler
 * * Dec  2 10:13 1996 (srk)
 * Created: August 1994
 *-------------------------------------------------------------------
 */

/* $Id: o2m.c,v 1.4 2017/01/19 21:53:21 mieg Exp $ */

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
#include "pmap.h"
#include "oxgrid.h"

BOOL o2mDisplay (void) ;
int o2mConvert (OX ox) ;   
static void o2mDraw (OX ox) ;
static void o2mDestroy (void) ;
static void drawChromosome (OX ox, float xlen) ;
static void drawPoints (OX ox, float xscale) ;
static void drawLines (OX ox, float xscale) ;
static void o2mPick (int box) ;
static void o2mSelect (OX ox, int box) ;
static void o2mFollow (OX ox) ;
static void drawScale (OX ox) ;
static void toggleCells (void) ;
static void toggleChroms (void) ;
static void noChroms (void) ;
static void o2mClear (void) ;
static void o2mLightAll (void) ;
static void o2mLight (void) ;
static void o2mUnLight (void) ;
static void o2mHide (void) ;
static void o2mHideLight (void) ;
static void o2mUnHide (void) ;
static void o2mExport (void) ;
static void o2mLightHom (void) ;
static void o2mUnLightHom (void) ;
static void o2mFindHom (char *string) ;
static void o2mFindLoc (char *string) ;
static int o2mhomolCompletion (char *cp, int len) ;
static int o2mlocusCompletion (char *cp, int len) ;
static void o2mZoomIn (void) ;
static void o2mZoomOut (void) ;


BOOL COLUMN = FALSE ;
BOOL ROW = FALSE ;

KEY _Prefix ;
int _VPrefix ;

enum BoxNames { BACKGROUND=0, 
                O2MEQ_BOX, CHROM_BOX, HIDE_BOX,
                FINDHOM_BOX, FINDHOM2_BOX,
                FINDLOC_BOX, FINDLOC2_BOX,
                ZOOM_IN, ZOOM_OUT,
		MESS2_BOX, MIN_BOX } ;  /* MIN_BOX must be last */

static MENUOPT o2mMenu[] =
            { {graphDestroy, "Quit"},
	      {help, "Help"},
	      {graphPrint, "Print"},
	      {displayPreserve, "Preserve"},
	      {0, 0}
            } ;

static MENUOPT Highlight2Menu[] =
            { {o2mClear, "Clear and show all"},
              {o2mLightAll, "Highlight all"},
              {o2mLightHom, "Highlight homology"},
              {o2mUnLightHom, "Unhighlight homology"},
              {o2mLight, "Highlight selected keyset"},
              {o2mUnLight,"Unhighlight selected keyset"},
              {o2mHide, "Restrict to keyset"},
              {o2mHideLight, "Hide highlit items"},
              {o2mUnHide, "Unhide hidden items"},
              {o2mExport, "Export highlit keyset"},
              {0, 0}
            } ;

static float xOff, yOff ;

/*********************************/

void o2mCreate (void)
{ 
  BOX *bx ;
  int box ;
  OXGET("o2mCreate") ;

  ox->size3 = arrayCreate (64, int) ;
  ox->axis3 = arrayCreate (64, KEY) ;
  COLUMN = FALSE ; ROW = FALSE ; 

  if (!(ox->gridBox))
    {messout ("You have not selected a cell");
     return;}
  else
    {box = ox->gridBox ;
     if (box >= ox->mlb + ox->map_labels)
      {bx = arrayp (ox->boxes, (box - ox->mlb + 1), BOX) ;
       if (messQuery ("Column (y) or Row (n)"))
         {ox->c3 = bx->map1 ;
          ox->s5 = ox->s2 ;
          COLUMN = TRUE ;
          ox->size3 = arrayCopy (ox->size2) ;
          ox->axis3 = arrayCopy (ox->axis2) ;
         }
       else
         {ox->c3 = bx->map2 ;
          ox->s5 = ox->s1 ;
          ROW = TRUE ;
          ox->size3 = arrayCopy (ox->size1) ;
          ox->axis3 = arrayCopy (ox->axis1) ;
         }
          
       if (!o2mDisplay ()) return ; ;
       }
     }
 
}


BOOL o2mDisplay (void)
{
  int h ;

  OXGET("o2mDisplay") ;

  h = o2mConvert (ox) ;
  if (!(h > 0)) 
    {messout ("Not enough data to draw a O2MMap") ;
     return FALSE ;}

  if (graphExists (ox->manyGraph))
    { graphActivate (ox->manyGraph) ;
      graphPop () ;
    }
  else 
    { if (!displayCreate (O2MMAP))
	return FALSE ;
      ox->manyGraph = graphActive () ;
    
      graphAssociate (&OX_MAGIC, ox) ;
      graphRegister (DESTROY, o2mDestroy) ;      
      graphRegister (MESSAGE_DESTROY, displayUnBlock) ;
      graphRegister (PICK, (GraphFunc) o2mPick) ;
      graphMenu (o2mMenu) ;
      graphTextBounds (180, 120) ;
    }

  
  o2mDraw (ox) ;

  return TRUE ;


}

/************************************************************/
/***************** Registered routines *********************/

/*************************************************/
/************* conversion and drawing ************/

static void o2mDestroy (void)
{
  OXGET ("o2mDestroy") ;

  arrayDestroy (ox->many)   ;
  arrayDestroy (ox->csegs)  ;
  arrayDestroy (ox->c3segs) ;
  arrayDestroy (ox->axis3)  ;
  arrayDestroy (ox->size3)  ;
  ox->s5 = 0 ;
  ox->c3 = 0 ;

}


int o2mConvert (OX ox)
{  
  KEY band, col, chr, key ;
  KEY map1, map2 ;
  OBJ Band, Chr, yChrom, yBand ;
  OBJ loc1, loc2 ;
  float x1, x2, ymin = 0, ymax = 0;
  float min2 = 0, max2 = 0 ;
  float f1, f2, f3, f4 ;
  int c, h, i, j, m ;
  int cmin=0, cmax=0, null3=0 ;
  CHROMS *c3seg ;
  BANDS *cseg = 0;
  char *sym ;
  SEG *mseg ;
  HOMOL *seg ;
  Array bands, sizes ;
  static BSMARK mark1 = 0 ;
  static BSMARK mark2 = 0 ;

  ox->flag |= FLAG_O2M_CHROM ;
  ox->flag |= FLAG_O2M_EQUAL ;

  ox->o2mwide = 95 ;
  ox->o2mhigh = 50 ;

  if (!(yChrom = bsCreate (ox->c3)))
    { messout ("No data for chromosome %s", name (ox->c3)) ;
      return 0 ;
    }
  if(bsFindTag (yChrom, _Non_graphic))
    { bsDestroy (yChrom) ;
      ox->flag |= ~FLAG_O2M_CHROM ;
      return 0 ;
    }
  bsDestroy (yChrom) ;

  bands = arrayCreate (8, BSunit);
  sizes = arrayCreate (8, KEY) ;
  ox->c3segs = arrayCreate (8, CHROMS) ;
  ox->csegs = arrayCreate (8, BANDS) ;

  null3 = 0 ;
  for (i = 1 ; i < arrayMax (ox->size3) ; i++)
    {if (array (ox->size3, i, int) == 0) null3++ ;}
  
  if (null3 > 0) 
   {ox->flag |= FLAG_O2M_SIZE ;}

/* BANDS data not yet utilised but can be used to draw points */
/* at their relative position along the y axis as in the TRANSGRID */

  for (j = 1 ; j < arrayMax (ox->axis3) ; j++) 
     {c = 0 ;
      chr = array (ox->axis3, j, KEY) ;
      if ((Chr = bsCreate (chr)))
       {if (bsFindTag (Chr, _Chrom_Band) && bsFlatten (Chr, 1 , bands))
         {for (i = 0 ; i < arrayMax (bands) ; ++i )
           {band = array (bands, i, BSunit).k ;
            Band = bsCreate (band) ; 
            if (bsFindTag (Band, _Map) && bsGetKey (Band, _bsRight, &key))
              if (bsPushObj (Band))
                if (bsGetData (Band, _Left, _Float, &x1) &&
                    bsGetData (Band, _Right, _Float, &x2))
                  {cseg = arrayp (ox->csegs, j, BANDS) ;
                   cseg->chrom = j ;
                   cseg->key = key ;
                   cseg->bands = arrayCreate (8, BSunit) ;
                   array (cseg->bands, c, BSunit).k = band ;
                   array (cseg->bands, c+1, BSunit).f = x1 ;
                   array (cseg->bands, c+2, BSunit).f = x2 ;
                   c += 3 ;}
             bsDestroy (Band) ;
             if (x1 < cmin)
               cmin = x1 ;
             if (x2 > cmax)
             cmax = x2 ;
             }
          cseg->xmin = cmin ;
          cseg->xmax = cmax ;
          }
        }
      bsDestroy (Chr) ;
     }

  arrayDestroy (bands) ;
  bands = arrayCreate (8, BSunit) ;

  if ((yChrom = bsCreate(ox->c3)) )
   {if (bsFindTag (yChrom, _Chrom_Band) && bsFlatten (yChrom, 1 , bands))
     {for (i = 0 ; i < arrayMax (bands) ; i++ )
       {band = array (bands, i, BSunit).k; 
        yBand = bsCreate (band);
        if (bsFindTag (yBand, _Map) &&
            bsGetKey (yBand, _bsRight, &key)) 
          if (bsPushObj(yBand))
            if (bsGetData (yBand, _Left, _Float, &x1) &&
                bsGetData (yBand, _Right, _Float, &x2))
           {c3seg = arrayp (ox->c3segs, i, CHROMS) ;
            c3seg->key = band ;
            c3seg->x = 0.5*(x1+ x2) ;
            c3seg->dx = x2 - x1 ;
	    c3seg->flag = 0 ;
            ox->flag |= FLAG_O2M_BANDS ;}
        bsDestroy (yBand) ;
        if (x1 < ymin)
          ymin = x1 ;
        if (x2 > ymax)
	  ymax = x2 ;
       }

      for (i = 0 ; i < arrayMax (bands) ; i++ )
       {band = array (bands, i, BSunit).k ;
        yBand = bsCreate (band);
        c3seg = arrayp (ox->c3segs, i, CHROMS) ;
        if (bsFindTag (yBand,_Centromere))
	    c3seg->flag |= FLAG_CENTROMERE ;
        if (bsFindTag (yBand,_Dark))
	    c3seg->flag |= FLAG_DARK_BAND ;
	if (bsFindTag (yBand,_NOR))
	    c3seg->flag |= FLAG_NOR ;
	if (bsFindTag (yBand,_Colour) &&
	    bsGetKey (yBand, _bsRight, &col))
	    c3seg->flag |= ((col - _WHITE) & 15  ) << 28;
        if (bsFindTag (yBand,_Symbol)) 
           {bsGetData (yBand, _bsRight, _Text, &sym);
            c3seg->sym = sym ;}
        bsDestroy (yBand) ;
       }   
    
      ox->c3min = ymin;
      ox->c3max = ymax; 
      bsDestroy (yChrom) ;
     }
 

  }       

  i = 0 ; m = 0 ; h = 0 ;

  ox->many = arrayCreate (8, HOMOL) ;

  for (i = 0 ; i < arrayMax (ox->segs) ; i++ )
    {mseg = arrayp (ox->segs, i, SEG) ;
     if (COLUMN)
      {if (mseg->map1 == ox->c3)
        {seg = arrayp (ox->many, m, HOMOL) ;
         seg->hom = mseg->key ;
         seg->key1 = mseg->gene1 ;
         seg->key2 = mseg->gene2 ;
         seg->flag = 0 ;
         if (mseg->flag & FLAG_HIDE) seg->flag |= FLAG_HOMHIDE ;
         seg->y = mseg->y ;
         seg->x1 = 0 ; seg->dx1 = 0 ; 
         seg->x2 = 0 ; seg->dx2 = 0 ;
         m++ ;

         if ((loc1 = bsCreate (seg->key1)))
         {if (bsGetKey (loc1, _Map, &map1)) do
           {mark1 = bsMark (loc1, mark1) ;
            if (chkSpecies (map1, ox->c3))
              {if (bsPushObj (loc1) &&bsGetData (loc1, _Position, _Float, &f1)) 
                {seg->x2 = f1 ;
                 seg->flag |= FLAG_YPOS ;
                 if (bsPushObj (loc1) && bsGetData (loc1, _Error, _Float, &f2))
                   seg->dx2 = f2 ;
                 else seg->dx2 = 0 ;
                 if ((seg->x2-seg->dx2) < min2) 
                   min2 = (seg->x2-seg->dx2) ;
                 if ((seg->x2+seg->dx2) > max2) 
                   max2 = (seg->x2+seg->dx2) ;
                 }
               }
             }while (bsGoto (loc1, mark1), bsGetKey (loc1, _bsDown, &map1)) ;
           }


         if ((loc2 = bsCreate (seg->key2)))
         {if (bsGetKey (loc2, _Map, &map2)) do
           {mark2 = bsMark (loc2, mark2) ;
            if (chkSpecies (map2, ox->c3))
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
           }

         if (seg->flag & FLAG_YPOS) h++ ;
         }
       }
     else if (ROW)
      {if (mseg->map2 == ox->c3)     
        {seg = arrayp (ox->many, m, HOMOL) ;
         seg->hom = mseg->key ;
         seg->key1 = mseg->gene1 ;
         seg->key2 = mseg->gene2 ;
         seg->flag = 0 ;
         if (mseg->flag & FLAG_HIDE) seg->flag |= FLAG_HOMHIDE ;
         seg->y = mseg->x ;
         seg->x1 = 0 ; seg->dx1 = 0 ; 
         seg->x2 = 0 ; seg->dx2 = 0 ;
         m++ ;

         if ((loc1 = bsCreate (seg->key1)))
          {if (bsGetKey (loc1, _Map, &map1)) do
           {mark1 = bsMark (loc1, mark1) ;
            if (chkSpecies (map1, ox->c3))
              {if (bsPushObj (loc1) &&bsGetData (loc1, _Position, _Float, &f1)) 
                {seg->x2 = f1 ;
                 seg->flag |= FLAG_YPOS ;
                 if (bsPushObj (loc1) && bsGetData (loc1, _Error, _Float, &f2))
                   seg->dx2 = f2 ;
                 else seg->dx2 = 0 ;
                 if ((seg->x2-seg->dx2) < min2) 
                   min2 = (seg->x2-seg->dx2) ;
                 if ((seg->x2+seg->dx2) > max2) 
                   max2 = (seg->x2+seg->dx2) ;
                 }
               }
             }while (bsGoto (loc1, mark1), bsGetKey (loc1, _bsDown, &map1)) ;
           }


         if ((loc2 = bsCreate (seg->key2)))
          {if (bsGetKey (loc2, _Map, &map2)) do
           {mark2 = bsMark (loc2, mark2) ;
            if (chkSpecies (map2, ox->c3))
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
           }

         if (seg->flag & FLAG_YPOS) h++ ;      
         }
       }
     }

  if (!(ox->flag & FLAG_O2M_BANDS))
    {ox->c3min = min2-5 ;
     ox->c3max = max2+5 ;
    }

arrayDestroy (bands) ;
arrayDestroy (sizes) ;
return h ;
 
}

static void o2mDraw (OX ox)
{
  int   xmax = 0, ymax = 0 ;
  int   box, i, j, k ;
  float xl, yl, xlen = 0, ylen = 0 ;
  float w1, w2, wtot, h1, h2, htot ;
  float xscale, yscale ;
  CHROMS *c3seg ;

  graphClear () ;

  if (!(ox->flag & FLAG_O2M_EQUAL))
   {if (graphButton (" Equal  sizes ", toggleCells, 1, 3) != O2MEQ_BOX)
    messcrash ("box screwup at %d in o2mDraw", O2MEQ_BOX) ;}
  else
   {if (graphButton ("Relative sizes", toggleCells, 1, 3) != O2MEQ_BOX)
    messcrash ("box screwup at %d in o2mDraw", O2MEQ_BOX) ;}

  if (ox->flag & FLAG_O2M_BANDS)
    {if (ox->flag & FLAG_O2M_CHROM)
      {if (graphButton ("Hide Chromosomes", toggleChroms, 17, 3) != CHROM_BOX)
       messcrash ("box screwup at %d in o2mDraw", CHROM_BOX) ;}
    else
      {if (graphButton ("Draw Chromosomes", toggleChroms, 17, 3) != CHROM_BOX)
       messcrash ("box screwup at %d in o2mDraw", CHROM_BOX) ;}
     }
  else
    {if (graphButton ("Draw Chromosomes", noChroms, 17, 3) != CHROM_BOX)
      messcrash ("box screwup at %d in o2mDraw", CHROM_BOX) ;}

  if (graphButton ("Highlight...", o2mClear, 35, 3) != HIDE_BOX)
    messcrash ("box screwup at %d in o2mDraw", HIDE_BOX) ;

  graphBoxMenu (HIDE_BOX, Highlight2Menu) ;  

  graphText ("Find homology :", 1, 5) ;
  *ox->o2mhom = 0 ;
  if (graphCompScrollEntry (o2mhomolCompletion, ox->o2mhom, 20, 20, 17, 5, o2mFindHom) != FINDHOM_BOX)
    messcrash ("box screwup at %d in o2mDraw", FINDHOM_BOX) ;

  graphEntryDisable () ;

  graphText ("Find locus :", 41, 5) ;
  *ox->o2mloc = 0 ;
  if (graphCompScrollEntry (o2mlocusCompletion, ox->o2mloc, 20, 20, 54, 5, o2mFindLoc) != FINDLOC_BOX)
    messcrash ("box screwup at %d in o2mDraw", FINDLOC_BOX) ;

  graphEntryDisable () ;

  if (graphButton ("Zoom in", o2mZoomIn, 49, 3) != ZOOM_IN)
    messcrash ("box screwup at %d in o2mDraw", ZOOM_IN) ;

  if (graphButton ("Zoom out", o2mZoomOut, 58, 3) != ZOOM_OUT)
    messcrash ("box screwup at %d in o2mDraw", ZOOM_OUT) ;

/* Draw the header box to name a selected homology on screen */   

  ox->message2Box = MESS2_BOX ;
 *ox->message2Text = 0 ;
  ox->message2Box = graphBoxStart () ;
  graphTextPtr (ox->message2Text, 1, 7, 112) ;
  graphBoxEnd () ;
  graphBoxDraw (ox->message2Box, BLACK, LIGHTBLUE) ;

/* Now draw the titles */

  xOff = 0 ; yOff = 0 ;

  graphText (messprintf ("One-to-Many Map for map %s and map_set %s", name (ox->c3), name (ox->s5)), 1, 1);
  graphLine (0, 9.5, 200, 9.5) ; 

  xOff += 2 + strlen (name (ox->c3)) ;
  yOff += 10 ;

  xmax = arrayMax (ox->axis3) - 1 ;
  ymax = arrayMax (ox->c3segs) - 1 ;
  wtot = 0 ; htot = 0 ;
  if (!(ox->flag & FLAG_O2M_EQUAL))
   {for ( j = 1 ; j <= xmax ; j++) {wtot += array (ox->size3, j, int) ;}
    }
  else
   {wtot = xmax + 1;
    }
  htot = ox->c3max - ox->c3min ;
  xscale = ox->o2mwide/wtot ; 
  yscale = ox->o2mhigh/htot ;                                     
  graphColor (BLACK) ;


  for ( i = 1 ; i < arrayMax (ox->axis3) ; i++)
    {yl = strlen (name (array (ox->axis3, i, KEY))) ;
     if (yl > ylen) ylen = yl ;
     }
  ylen *= 0.65 ;
  yOff += 5 + ylen ;

  if (ox->flag & FLAG_O2M_BANDS)
   {for (i = 0 ; i < arrayMax (ox->c3segs) ; i++)
      {c3seg = arrayp (ox->c3segs, i, CHROMS) ;
       xl = strlen (name (c3seg->key)) ;
       if (xl > xlen) xlen = xl ;
       }
    xOff += 8 + xlen ;
    }
  else
    xOff += 10 ;

  box = MIN_BOX ;
  if (box != graphBoxStart ())
    messcrash ("minLiveBox wrong in o2mDraw ()") ;
  graphTextFormat (BOLD) ;
  graphText (name (ox->s5), xOff + (ox->o2mwide)/2, 10) ;
  graphBoxEnd () ;
  graphBoxDraw (box, BLACK, WHITE) ;
  box++ ;
  if (box != graphBoxStart ())
    messcrash ("minLiveBox wrong in o2mDraw ()") ;
  graphTextFormat (BOLD) ;
  graphText (name (ox->c3), 1, yOff + (ox->o2mhigh)/2) ;
  graphBoxEnd () ;
  graphBoxDraw (box, BLACK, WHITE) ;
  box ++ ;
  graphTextFormat(PLAIN_FORMAT);
  
  if (ox->flag & FLAG_O2M_BANDS)
    {for (i = 1 ; i <= xmax ; i++)
      {for ( j = 0 ; j <= ymax ; j++)
        {w1 = 0 ; h1 = 0 ; w2 = 0 ; h2 = 0 ;
         if (!(ox->flag & FLAG_O2M_EQUAL))
           {for ( k = 0 ; k <= i ; k++) {w2 += array (ox->size3, k , int) ;}
            for ( k = 1 ; k <= i ; k++) {w1 += array (ox->size3, k-1 , int) ;}
           }
         else
           {w2 = i ;
            w1 = (i - 1) ;
           }
         for ( k = 0 ; k <= j ; k++) {h2 += (arrayp (ox->c3segs, k, CHROMS)->dx) ;}
         for ( k = 1 ; k <= j ; k++) {h1 += (arrayp (ox->c3segs, k-1, CHROMS)->dx) ;}
         graphRectangle (xOff+(w1*xscale), yOff+(h1*yscale), xOff+(w2*xscale), yOff+(h2*yscale)) ;
         }
       } 
     }
  else
    {for (i = 1 ; i <= xmax ; i++)
      {w1 = 0 ; h1 = 0 ; w2 = 0 ; h2 = 0 ;
       if (!(ox->flag & FLAG_O2M_EQUAL))
           {for ( k = 0 ; k <= i ; k++) {w2 += array (ox->size3, k , int) ;}
            for ( k = 1 ; k <= i ; k++) {w1 += array (ox->size3, k-1 , int) ;}
           }
       else
           {w2 = i ;
            w1 = (i - 1) ;
           }
       graphRectangle (xOff+(w1*xscale), yOff+ox->o2mhigh, xOff+(w2*xscale), yOff) ;
      }
    }

  if (ox->flag & FLAG_O2M_BANDS)
    {for (i = 1 ; i <= xmax ; i++)
       {w1 = 0 ; h1 = 0 ; w2 = 0 ; h2 = 0 ;
        if (!(ox->flag & FLAG_O2M_EQUAL))
          {for ( k = 0 ; k <= i ; k++) {w2 += array (ox->size3, k , int) ;}
           for ( k = 1 ; k <= i ; k++) {w1 += array (ox->size3, k-1 , int) ;}
          }
        else
          {w2 = i ;
           w1 = (i - 1) ;
          }
        if (box != graphBoxStart ())
          messcrash ("minLiveBox wrong in o2mDraw ()") ;
        graphTextHeight (0.5) ;
        japanese (name (array (ox->axis3, i, KEY)), xOff+(w1+w2)*0.5*xscale, yOff-3-ylen) ;
        graphBoxEnd () ;
        graphBoxDraw (box, BLUE, WHITE) ;
        box ++ ;
       } 
     }
  else
    {for (i = 1 ; i <= xmax ; i++)
      {w1 = 0 ; h1 = 0 ; w2 = 0 ; h2 = 0 ;
       if (!(ox->flag & FLAG_O2M_EQUAL))
           {for ( k = 0 ; k <= i ; k++) {w2 += array (ox->size3, k , int) ;}
            for ( k = 1 ; k <= i ; k++) {w1 += array (ox->size3, k-1 , int) ;}
           }
         else
           {w2 = i ;
            w1 = (i - 1) ;
           }
       if (box != graphBoxStart ())
         messcrash ("minLiveBox wrong in o2mDraw ()") ;
       graphTextHeight (0.5) ;
       japanese (name (array (ox->axis3, i, KEY)), xOff+(w1+w2)*0.5*xscale, yOff-3-ylen) ;
       graphBoxEnd () ;
       graphBoxDraw (box, BLUE, WHITE) ;
       box ++ ;
      }
    }

  graphTextHeight (1) ;
  
  if (ox->flag & FLAG_O2M_BANDS)
    {if (ox->flag & FLAG_O2M_CHROM)
       {ox->o2m_labels = 2 + xmax + ymax + 1 ;
        drawChromosome (ox, xlen) ;
       }
     else
       {ox->o2m_labels = 2 + xmax ;
        drawScale (ox) ;
       }
    }
  else
    {ox->o2m_labels = 2 + xmax ;
     drawScale (ox) ;
    }

  drawPoints (ox, xscale); 
  drawLines (ox, xscale);

  ox->manyBox = 0 ;

  graphRedraw () ;
}

static void drawChromosome (OX ox, float xlen)
{
  float y0, dy, y, yend ;
  int i, box, xmax ;
  CHROMS *c3seg ;

  y0 = -(ox->c3min) ;
  yend = (ox->c3max - ox->c3min) / ox->o2mhigh ;
  xmax = arrayMax (ox->axis3) - 1 ;

/* Draw the chromosome alongside the grid */

  box = MIN_BOX + 2 + xmax ;
  for (i = 0 ; i < arrayMax (ox->c3segs) ; ++i)
    {c3seg = arrayp (ox->c3segs, i, CHROMS) ;
     if (class (c3seg->key) == _VChrom_Band)
       {
	y = yOff + ( y0 + c3seg->x) / yend ;
	dy = ( 0.5 * c3seg->dx) / yend ;
	if (c3seg->flag & FLAG_CENTROMERE)
	  { graphColor (WHITE) ;
	    graphFillRectangle (xOff-4, y-dy, xOff, y+dy) ;
	    graphColor (BLACK) ;
	    graphLine (xOff-4, y-dy, xOff, y+dy) ;
	    graphLine (xOff, y-dy, xOff-4, y+dy) ;
	    graphLine (xOff-4, y-dy, xOff, y-dy) ;
	    graphLine (xOff, y+dy, xOff-4, y+dy) ;
            if (box != graphBoxStart ())
            messcrash ("minLiveBox wrong in o2mDraw () box %d", box) ;
            graphTextHeight (0.5) ;
            graphText (name (c3seg->key), xOff-6-xlen, y-0.25) ;
            graphBoxEnd () ;
            graphBoxDraw (box, BLUE, WHITE) ;
            box++ ; 
            graphTextHeight (1) ;
           }	      
	else 
	  { if (c3seg->flag & FLAG_DARK_BAND)
             graphColor (DARKGRAY) ;
	    else if (c3seg->flag & FLAG_NOR)
	     graphColor (LIGHTGRAY) ;
	    else
	     graphColor (WHITE) ;
	    if (c3seg->flag >> 28)
	     graphColor ((c3seg->flag  >> 28) + WHITE) ;
	    graphFillRectangle (xOff-4, y-dy, xOff, y+dy) ;
	    graphColor (BLACK) ;
	    graphRectangle (xOff-4, y-dy, xOff, y+dy) ;
            if (box != graphBoxStart ())
            messcrash ("minLiveBox wrong in o2mDraw () box %d", box) ;
            graphTextHeight (0.5) ;
            graphText (name (c3seg->key), xOff-6-xlen, y-0.25) ;
            graphBoxEnd () ;
            graphBoxDraw (box, BLUE, WHITE) ;
            box++ ; 
            graphTextHeight (1) ;
	  }
       }
    }

  graphTextHeight (0) ;

}

static void drawScale (OX ox) 
{
  int i, j, k ;
  float unit, y0 ;

/* If no banding data or if the chromosomes are hidden, draw a scale instead */

  graphTextHeight (0.5) ;
  y0 = yOff + ox->o2mhigh * (1-(ox->c3max / (ox->c3max-ox->c3min))) ;
  unit = ox->o2mhigh / (ox->c3max - ox->c3min) ;
  for (i = 0 ; i < (ox->c3max + 1) ; i++)
    graphLine (xOff-3, y0+i*unit, xOff-3.9, y0+i*unit) ;
  for (k = 0 ; k < (ox->c3max + 1) ; k+=10)
   {graphLine (xOff-3, y0+k*unit, xOff-4.5, y0+k*unit) ;
    graphText (messprintf ("%3d",k), xOff - 8, y0+k*unit-0.5) ;
    }
  for (j = 0 ; j < -(ox->c3min - 1) ; j++)
    graphLine (xOff-3, y0-j*unit, xOff-3.9, y0-j*unit) ;
  for (k = 0 ; k < -(ox->c3min - 1) ; k+=10)
   {graphLine (xOff-3, y0-k*unit, xOff-4.5, y0-k*unit) ;
    graphText (messprintf ("%3d",-k), xOff - 8, y0-k*unit-0.5) ;
    }  
  graphLine (xOff-3, y0-(j-1)*unit, xOff-3, y0+(i-1)*unit) ;
  graphTextHeight (1) ;
 
}


static void drawPoints (OX ox, float xscale)
{
  int c, i, k, m, w1, w2 ;
  float x, y, y0, yend ;
  HOMOL *seg ;
   
  y0 = -(ox->c3min) ;
  yend = (ox->c3max - ox->c3min) / ox->o2mhigh ;
  m = MIN_BOX + ox->o2m_labels ;
  for (i = 0 ; i < arrayMax (ox->many) ; i++)
    {seg = arrayp (ox->many, i, HOMOL) ;
     if (seg->flag & FLAG_YPOS)
      {if (!(seg->flag & FLAG_HOMHIDE))
        {w1 = 0 ; w2 = 0 ;
         c = seg->y ; 
         if (!(ox->flag & FLAG_O2M_EQUAL))
           {for ( k = 0 ; k <= c ; k++) {w2 += array (ox->size3, k , int) ;}
            for ( k = 1 ; k <= c ; k++) {w1 += array (ox->size3, k - 1 , int) ;}
           }
         else
           {w2 = c ;
            w1 = (c - 1) ;
           }
         x = xOff + (0.5*(w2+w1)*xscale) ;
         y = yOff + (y0 + seg->x2) / yend ;
         if (m != graphBoxStart ())
           messcrash ("problems with One-to-Many-Map point %d", m) ;
         graphRectangle (x-0.4, y-0.4, x+0.4, y+0.4) ;
         graphBoxEnd () ;
         if (seg->flag & FLAG_HIGHLIGHT) graphBoxDraw (m, WHITE, MAGENTA) ;
         else graphBoxDraw (m, WHITE, BLUE) ;
         seg->pbox = m ;
         m++ ;}  
        }    
     }
}

int pOrder (const void *a, const void *b)
{ 
  const POINT *point = (const POINT*)a ;
  const POINT *point2 = (const POINT*)b ;
  float x ;

  x = point->y - point2->y ;
  if (x == 0) return 0 ;
  else if (x > 0) return 1 ;
  return -1 ;
}

static void drawLines (OX ox, float xscale)
{
  int i, j, k, l, w1, w2 ;
  float x, y0, y1, y2, yend ;
  HOMOL *seg ;
  POINT *point, *point2 ;
  Array tpoints ;

  ox->segsize = 25 ;

  graphColor (MAGENTA) ;
  y0 = -(ox->c3min) ;
  yend = (ox->c3max - ox->c3min) / ox->o2mhigh ;
  for (i = 0 ; i < arrayMax (ox->axis3) ; i++)
    {tpoints = arrayCreate (64, POINT) ;
     k = 0 ;
     for (j = 0 ; j < arrayMax (ox->many) ; j++)
       {seg = arrayp (ox->many, j, HOMOL) ;
        if (seg->y == i && seg->x2 != 0)
	  {point = arrayp (tpoints, k, POINT) ;
           point->y = seg->x2 ;
           k++ ;
          }
       }   

     arraySort (tpoints, pOrder) ;
     for (k = 0 ; k < arrayMax (tpoints) - 1 ; k++)
       {point = arrayp (tpoints, k, POINT) ;
        point2 = arrayp (tpoints, k+1, POINT) ;  
        if ((point2->y - point->y) < ox->segsize) 
          {if (!(ox->flag & FLAG_O2M_EQUAL))
             {w1 = 0 ; w2 = 0 ;
              for ( l = 0 ; l <= i ; l++) {w2 += array (ox->size3, l , int) ;}
              for ( l = 1 ; l <= i ; l++) {w1 += array (ox->size3, l - 1 , int) ;}
             }
           else
             {w2 = i ;
              w1 = (i - 1) ;
             }
           x = xOff + (0.5*(w2+w1)*xscale) ;
           y1 = yOff + (y0 + point->y) / yend ;
           y2 = yOff + (y0 + point2->y) / yend ;
           graphLine (x, y1, x, y2) ;
	  }
       }
     arrayDestroy (tpoints) ;
    }
  graphColor (BLACK) ;
}

extern int ksetClassComplete (char *text, int len, int classe) ;

static int o2mhomolCompletion (char *cp, int len)
{
  return ksetClassComplete (cp, len, _VHomology_group) ;
}

static int o2mlocusCompletion (char *cp, int len)
{
  return ksetClassComplete (cp, len, _VLocus) ;
}

static void o2mFindHom (char* string)
{
  int i ;
  KEY key ;
  HOMOL *seg ;
  BOOL found = FALSE ;

  OXGET ("o2mFindHom") ;

  if (!(lexword2key (string, &key, _VHomology_group)))
    {messout ("Sorry, %s is not a homology name", string) ;
     return ;
    }

  for (i = 0 ; i < arrayMax (ox->many) ; i++)
    {seg = arrayp (ox->many, i, HOMOL) ;
     if (chkSpecies (seg->hom, key)) 
       {found = TRUE ;
        o2mSelect (ox, seg->pbox) ;
       }
     }

  if (found == FALSE) 
   {messout ("Sorry, this homology is not in this grid") ;
    return ;
   }

}

static void o2mFindLoc (char* string)
{
  int i ;
  KEY key ;
  HOMOL *seg ;
  BOOL found = FALSE ;

  OXGET ("o2mFindLoc") ;

  if (!(lexword2key (string, &key, _VLocus)))
    {messout ("Sorry, %s is not a locus name", string) ;
     return ;
    }

  for (i = 0 ; i < arrayMax (ox->many) ; i++)
    {seg = arrayp (ox->many, i, HOMOL) ;
     if (chkSpecies (seg->key1, key) || chkSpecies (seg->key2, key)) 
       {found = TRUE ;
        o2mSelect (ox, seg->pbox) ;
       }
     }

  if (found == FALSE) 
   {messout ("Sorry, this locus is not in this grid") ;
    return ;
   }

}

static void o2mZoomOut (void) 
{
  OXGET ("o2mZoomOut") ;

  if (ox->o2mwide >= 30 && ox->o2mhigh >= 20)
    {ox->o2mwide -= 15 ;
     ox->o2mhigh -= 10 ;

     o2mDraw (ox) ;
    }

}

static void o2mZoomIn (void) 
{
  OXGET ("o2mZoomIn") ;

  if (ox->o2mwide <= 120 && ox->o2mhigh <= 80)
    {ox->o2mwide += 15 ;
     ox->o2mhigh += 10 ;

     o2mDraw (ox) ;
    }

}

static void o2mPick (int box) 
{ 
  int xmax ; 

  OXGET ("o2mPick") ;

  xmax = arrayMax (ox->axis3) - 1 ;

  /* The usual - check if the box number is legal, select it if 
     a single click and follow it if a double */

  if (box >= MIN_BOX + ox->o2m_labels && 
      box <= (arrayMax (ox->many) + MIN_BOX + ox->o2m_labels)) /* a seg */
   {if (box == ox->manyBox)
      o2mFollow (ox) ;
    graphActivate (ox->manyGraph) ;
    o2mSelect (ox, box) ; 
   }
  else if (box == MIN_BOX)
    {if (box == ox->manyBox)
       display (ox->s5, 0, 0) ;
     graphActivate (ox->manyGraph) ;
     o2mSelect (ox, box) ;
    }
  else if (box == MIN_BOX + 1)
    {if (box == ox->manyBox)
       display (ox->c3, 0, 0) ;
     graphActivate (ox->manyGraph) ;
     o2mSelect (ox, box) ;
    }
  else if (box >= MIN_BOX + 2 && box < MIN_BOX + 2 + xmax)
    {if (box == ox->manyBox)
       display (array (ox->axis3, box - MIN_BOX - 1, KEY), 0, 0) ;
     graphActivate (ox->manyGraph) ;
     o2mSelect (ox, box) ;
    }
  else if (box >= MIN_BOX + 2 + xmax && box < ox->o2m_labels + MIN_BOX)
    {if (box == ox->manyBox)
       display (arrayp (ox->c3segs, box-MIN_BOX-2-xmax, CHROMS)->key, 0, 0) ;
     graphActivate (ox->manyGraph) ;
     o2mSelect (ox, box) ;
    }
  else if (box == FINDHOM_BOX)
    graphCompScrollEntry (0, ox->o2mhom, 0, 0, 0, 0, 0) ;
  else if (box == FINDLOC_BOX)
    graphCompScrollEntry (0, ox->o2mloc, 0, 0, 0, 0, 0) ;
}

static void o2mSelect (OX ox, int box)
{
  int i, h ;
  KEY group = 0, hit=0, key1 = 0, key2 = 0 ;
  HOMOL *seg ;
  HIDBOX *temp ;
  BOOL HID = FALSE ;
  Array oldhomols = 0 ;
  
  /* The usual - return the old box to normal and colour in the new */
  
  if (ox->manyBox == box) return ;
  
  if (ox->manyBox)
    {oldhomols = arrayCreate (8, HIDBOX) ;
    oldhomols = arrayCopy (ox->ohomols) ;
    ox->ohomols = arrayReCreate (ox->ohomols, 8, HIDBOX) ;
    }
  else 
    ox->ohomols = arrayCreate (8, HIDBOX) ;
  
  for (i = 0 ; i < arrayMax (ox->many) ; ++i)
    {seg = arrp (ox->many, i, HOMOL) ;
    if (seg->pbox == ox->manyBox)
      if (seg->flag & FLAG_HIGHLIGHT) HID = TRUE ;
    if (seg->pbox == box) group = seg->hom ;
    }
  
  h = 0 ;
  for (i = 0 ; i < arrayMax (ox->many) ; ++i)
    {seg = arrp (ox->many, i, HOMOL) ;
    if (chkSpecies (group, seg->hom) && seg->pbox != box)
      {temp = arrayp (ox->ohomols, h, HIDBOX) ;
      temp->box = seg->pbox ;
      if (seg->flag & FLAG_HIGHLIGHT) temp->hid = TRUE ;
      else temp->hid = FALSE ;
      h++ ;
      }
    }
  
  if (ox->manyBox) 
    {if (ox->manyBox >= MIN_BOX + ox->o2m_labels)
      {if (HID == TRUE) graphBoxDraw (ox->manyBox, WHITE, MAGENTA) ;
      else graphBoxDraw (ox->manyBox, WHITE, BLUE) ;
      for (h = 0 ; h < arrayMax (oldhomols) ; h++)
        {temp = arrayp (oldhomols, h, HIDBOX) ;
	if (temp->hid == TRUE) graphBoxDraw (temp->box, WHITE, MAGENTA) ;
	else graphBoxDraw (temp->box, WHITE, BLUE) ;
	}
      }
    else if (ox->manyBox == MIN_BOX || ox->manyBox == MIN_BOX + 1)
      graphBoxDraw (ox->manyBox, BLACK, WHITE) ;
    else if (ox->manyBox >= MIN_BOX + 2 && ox->manyBox < MIN_BOX + ox->o2m_labels)
      graphBoxDraw (ox->manyBox, BLUE, WHITE) ;
    }
  
  ox->manyBox = box ;
  
  if (ox->manyBox >= MIN_BOX + ox->o2m_labels)
    {graphBoxDraw (ox->manyBox, WHITE, GREEN) ;
    for (h = 0 ; h < arrayMax (ox->ohomols) ; h++)
      graphBoxDraw (arrayp (ox->ohomols, h, HIDBOX)->box, WHITE, RED) ;
    }
  else if (ox->manyBox == MIN_BOX || ox->manyBox == MIN_BOX + 1)
    graphBoxDraw (ox->manyBox, WHITE, BLACK) ;
  else if (ox->manyBox >= MIN_BOX + 2 && ox->manyBox < MIN_BOX + ox->o2m_labels)
    graphBoxDraw (ox->manyBox, WHITE, BLUE) ;
  
  graphCompScrollEntry (0, ox->o2mhom, 0, 0, 0, 0, 0) ; 
  graphCompScrollEntry (0, ox->o2mloc, 0, 0, 0, 0, 0) ;
  graphEntryDisable () ;
  
  if (ox->message2Box)
    {
      if (ox->manyBox >= MIN_BOX + ox->o2m_labels)
      {
	for (i=0 ; i < arrayMax (ox->many) ; i++)
	  {
	    seg = arrayp (ox->many, i, HOMOL) ;
	    hit = 0 ;
	    if (seg->pbox == ox->manyBox)
	      {
		hit = seg->hom ;
		key1 = seg->key1 ;
		key2 = seg->key2 ;
	      }
	    
	    if (hit)
	      {
		strncpy (ox->message2Text, messprintf ("Homology "), 80) ; 
		strcat (ox->message2Text, messprintf ("%s ", name (hit))) ;
		strcat (ox->message2Text, 
			messprintf ("between %s and %s", 
				    name (key1), name (key2))) ;
	      }
	  }
      }
    else 
      *ox->message2Text = 0 ;
    }
  graphBoxDraw (ox->message2Box, BLACK, LIGHTBLUE) ;
  arrayDestroy (oldhomols) ;
}
  
  static void o2mFollow (OX ox)
{
  int i ;
  KEY hit = 0 ;
  HOMOL *seg ;

  /* Similar to usual - seg->pbox needed to make the keyset option work */

  for (i=0 ; i < arrayMax (ox->many) ; i++)
    {seg = arrayp (ox->many, i, HOMOL) ;
     if (seg->pbox == ox->manyBox)
       hit = seg->hom ;
    }
  if (hit)
    display (hit, 0, 0) ;
}


static void toggleCells (void) 
{
  OXGET ("toggleCells") ;

  if (ox->flag & FLAG_O2M_SIZE)
    {messout ("Maps have not been assigned relative sizes") ;
     return ;}
  else
    {ox->flag ^= FLAG_O2M_EQUAL ;
     o2mDraw (ox) ;}


}

static void toggleChroms (void) 
{
  OXGET ("toggleChroms") ;
  
  if (!(ox->flag & FLAG_O2M_BANDS))
    {messout ("Sorry, there's no band data") ;
     return ;}
  else
    {ox->flag ^= FLAG_O2M_CHROM ;
     o2mDraw (ox) ;}

}

static void noChroms (void) 
{
  OXGET ("noChroms") ;

  messout ("Sorry, there's no band data") ;
  return ;

}

static void o2mClear (void) 
{
  int i ;
  HOMOL *seg ;
 
  OXGET ("o2mClear") ;

  seg = arrayp (ox->many, 0, HOMOL) ;
  i = arrayMax (ox->many) ;
  while(i--)
    {seg->flag &= ~FLAG_HOMHIDE ;  /* bitwise NOT */
     seg->flag &= ~FLAG_HIGHLIGHT ;
     seg++ ;
    }

  o2mDraw (ox) ;

}

static void o2mLightAll (void) 
{
  int i ;

  OXGET ("o2mLightAll") ;

  for (i = 0 ; i < arrayMax (ox->many) ; ++i)
    arrp (ox->many, i, HOMOL)->flag |= FLAG_HIGHLIGHT ;
    
  ox->manyBox = 0 ;

  o2mDraw (ox) ;

}

static void o2mLight (void) 
{
  int i, j ;
  HOMOL *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;

  OXGET ("o2mLight") ;

  if (!keySetActive (&keySet, &dummy))
    { messout("First select a keySet window, thank you.") ;
      return ;
    }
  for (i = 0 ; i < arrayMax (ox->many) ; ++i)
    { seg = arrp (ox->many, i, HOMOL) ;
      if (keySetFind (keySet, seg->hom, &j) ||
          keySetFind (keySet, seg->key1, &j) ||
          keySetFind (keySet, seg->key2, &j)) 
	seg->flag |= FLAG_HIGHLIGHT ;
    }

  o2mDraw (ox) ;

}

static void o2mUnLight (void) 
{
  int i, j ;
  HOMOL *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;

  OXGET ("o2mUnLight") ;

  if (!keySetActive (&keySet, &dummy))
    { messout ("First select a keySet window, thank you.") ;
      return ;
    }
  for (i = 0 ; i < arrayMax (ox->many) ; ++i)
    { seg = arrp (ox->many, i, HOMOL) ;
      if (keySetFind (keySet, seg->hom, &j) ||
          keySetFind (keySet, seg->key1, &j) ||
          keySetFind (keySet, seg->key2, &j)) 
	seg->flag &= ~FLAG_HIGHLIGHT ;
    }

  o2mDraw (ox) ;

}

static void o2mLightHom (void) 
{
  int i ;
  HOMOL *seg ;

  OXGET ("o2mLightHom") ;

  if (!(ox->manyBox))
    {messout ("You have not selected a homology") ;
     return ;}

  for (i = 0 ; i < arrayMax (ox->many) ; i++)
    {seg = arrayp (ox->many, i, HOMOL) ;
     if (seg->pbox == ox->manyBox) seg->flag |= FLAG_HIGHLIGHT ;
    }

  o2mDraw (ox) ;

}

static void o2mUnLightHom (void) 
{
  int i ;
  HOMOL *seg ;

  OXGET ("o2mUnLightHom") ;

  if (!(ox->manyBox))
    {messout ("You have not selected a homology") ;
     return ;}

  for (i = 0 ; i < arrayMax (ox->many) ; i++)
    {seg = arrayp (ox->many, i, HOMOL) ;
     if (seg->pbox == ox->manyBox) seg->flag &= ~FLAG_HIGHLIGHT ;
    }

  o2mDraw (ox) ;

}



static void o2mHide (void)
{
  int i, j ;
  HOMOL *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;
 
  OXGET ("o2mHide") ;

  if (!keySetActive (&keySet, &dummy))
    {messout ("First select a keySet window, thank you.") ;
     return ;
    }

  seg = arrp (ox->many, 0, HOMOL) ;
  i = arrayMax (ox->many) ;
  while(i--)
    {seg->flag &= ~FLAG_HOMHIDE ; 
     if (!keySetFind (keySet, seg->hom, &j) &&
         !keySetFind (keySet, seg->key1, &j) &&
         !keySetFind (keySet, seg->key2, &j)) 
       seg->flag |= FLAG_HOMHIDE ;
     seg++ ;
    }

  o2mDraw (ox) ;

}

static void o2mHideLight (void) 
{
  int i ;
  HOMOL *seg ;

  OXGET ("o2mHideLight") ;

  for (i = 0 ; i < arrayMax (ox->many) ; ++i)
    { seg = arrp (ox->many, i, HOMOL) ;
      if (seg->flag & FLAG_HIGHLIGHT)
	seg->flag |= FLAG_HOMHIDE ;
    }

  o2mDraw (ox) ;

}

static void o2mUnHide (void) 
{
  int i ;
  HOMOL *seg ;

  OXGET ("o2mUnHide") ;

  for (i = 0 ; i < arrayMax (ox->many) ; ++i)
    { seg = arrp (ox->many, i, HOMOL) ;
      seg->flag &= ~FLAG_HOMHIDE ;
    }

  o2mDraw (ox) ;


}

static void o2mExport (void) 
{
  int i, j, k = 0 ;
  HOMOL *seg ;
  static KEYSET kset = 0 ;
  BOOL found = FALSE ;

  OXGET ("o2mExport") ;

  kset = keySetReCreate (kset) ;

  for (i = 0 ; i < arrayMax (ox->many) ; ++i)
    { found = FALSE ;
      seg = arrp (ox->many, i, HOMOL) ;
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
 
 
 
 
 
 
 
