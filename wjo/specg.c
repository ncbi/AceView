/*  File: specg.c
 *  Written by Jo Dicks (jld@bioch.ox.ac.uk)
 *-------------------------------------------------------------------
 * This file is part of the ACeDB comparative mapping package, 
 * written by Jo Dicks and Michelle Kirby
 * The original ACeDB system was written by :
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *      Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *
 * Description:  to display Species Grids
 * Exported functions: oxDisplay () only
 * HISTORY:
 * Last edited: Feb 19 16:52 1997 (srk)
 *	- change NULL's to 0's to please MS Visual C++ compiler
 * Created: January 1995
 *-------------------------------------------------------------------
 */

/* $Id: specg.c,v 1.4 2017/01/19 21:53:21 mieg Exp $ */

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
 
int specConvert (OX ox) ; 
int pointOrder (const void *a, const void *b) ;
static void specPick (int box) ; 
static void specSelect (OX ox, int box) ;
static void specFollow (OX ox) ; 
static void specDestroy (void) ;   
static void specDraw (OX ox) ;
static void drawChromosome (OX ox, float xlen) ;
static void drawScale (OX ox) ;
static void drawSpots (OX ox, float xscale) ;
static void drawLines (OX ox, float xscale) ;
static void toggleChroms (void) ;
static void noChroms (void) ;
static void specZoomIn (void) ;
static void specZoomOut (void) ;
static void specClear (void) ;
static void specLightAll (void) ;
static void specLightHom (void) ;
static void specUnLightHom (void) ;
static void specLight (void) ;
static void specUnLight (void) ;
static void specHide (void) ;
static void specHideLight (void) ;
static void specUnHide (void) ;
static void specExport (void) ;
static void specFindHom (char *string) ;
static void specFindLoc (char *string) ;
static int spechomolCompletion (char *cp, int len) ;
static int speclocusCompletion (char *cp, int len) ;

enum BoxNames { BACKGROUND=0, 
                CHROM_BOX, HIDE_BOX, 
                ZOOM_IN, ZOOM_OUT,
                FINDHOM_BOX, FINDHOM2_BOX,
                FINDLOC_BOX, FINDLOC2_BOX,
                MESS_BOX, MIN_BOX } ;  /* MIN_BOX must be last */


static MENUOPT specMenu[] =
            { {graphDestroy, "Quit"},
 	      {help, "Help"},
	      {graphPrint, "Print"},
	      {displayPreserve, "Preserve"},
	       {0, 0}
            } ;

static MENUOPT Highlight3Menu[] =
            { {specClear, "Clear and show all"},
              {specLightAll, "Highlight all"},
              {specLightHom, "Highlight homology"},
              {specUnLightHom, "Unhighlight homology"},
              {specLight, "Highlight selected keyset"},
              {specUnLight,"Unhighlight selected keyset"},
              {specHide, "Restrict to keyset"},
              {specHideLight, "Hide highlit items"},
              {specUnHide, "Unhide hidden items"},
              {specExport, "Export highlit keyset"},
	       {0, 0}
            } ;

static float xOff, yOff ;

void specDisplay (void)
{
  int spec ;

  OXGET ("specDisplay") ;
  ox->specflag = 0 ;
  spec = specConvert (ox) ;

  if (spec == -1)
    {messout ("Sorry, definition is bad") ;
     return ;
    }
  else if (spec == 0)
    {messout ("No homologies for this Map") ;
     return ;
    }
  else
    {if (graphExists (ox->specGraph))
      {graphActivate (ox->specGraph) ;
       graphPop () ;}
     else
      {if (!displayCreate (SPECGRID))
       return ;
       graphAssociate (&OX_MAGIC, ox) ;
       ox->specGraph = graphActive () ;
       graphRegister (DESTROY, specDestroy) ;
       graphRegister (MESSAGE_DESTROY, displayUnBlock) ;
       graphRegister (PICK, (GraphFunc) specPick) ;
       graphMenu (specMenu) ;
       graphTextBounds (180, 120) ;}

     specDraw (ox) ;
    }
}


/************************************************************/
/***************** Registered rourtines *********************/

static void specDestroy (void)
{
  OXGET ("specDestroy") ;
 
  arrayDestroy (ox->specsegs) ;
  arrayDestroy (ox->specs) ;
  arrayDestroy (ox->c4segs) ;

  ox->c4 = 0 ;

  ox->specflag = 0 ;

  if (graphAssFind (&OX_MAGIC, &ox)
      && ox
      && (ox->magic == &OX_MAGIC))
     {if (graphExists (ox->specGraph))
	{ graphActivate (ox->specGraph) ;
	  graphDestroy () ;
	}
     }

  ox->chooseGraph = graphActive () ;

}
 
/**********************************/


static void specSelect (OX ox, int box)
{
  int i, h ;
  KEY group = 0, hit = 0, key1 = 0, key2 = 0;
  SEG *seg ;
  HIDBOX *temp ;
  BOOL HID = FALSE ;
  Array oldhomols = 0 ;

  if (ox->specBox == box) return ;

  if (ox->specBox)
    {oldhomols = arrayCreate (8, HIDBOX) ;
     oldhomols = arrayCopy (ox->homols) ;
     ox->homols = arrayReCreate (ox->homols, 8, HIDBOX) ;
    }
  else 
    {ox->homols = arrayCreate (8, HIDBOX) ;
    }
       
  for (i = 0 ; i < arrayMax (ox->specsegs) ; ++i)
    {seg = arrp (ox->specsegs, i, SEG) ;
     if (seg->box == ox->specBox)
       if (seg->flag & FLAG_LIGHT) HID = TRUE ;
     if (seg->box == box) group = seg->key ;
    }

  h = 0 ;
  for (i = 0 ; i < arrayMax (ox->specsegs) ; ++i)
    {seg = arrp (ox->specsegs, i, SEG) ;
     if (chkSpecies (group, seg->key) && seg->box != box)
       {temp = arrayp (ox->homols, h, HIDBOX) ;
        temp->box = seg->box ;
        if (seg->flag & FLAG_LIGHT) temp->hid = TRUE ;
        else temp->hid = FALSE ;
        h++ ;
       }
    }


  if (ox->specBox) 
   {if (ox->specBox >= MIN_BOX + ox->spec_labels)
     {if (HID == TRUE) graphBoxDraw (ox->specBox, WHITE, MAGENTA) ;
      else graphBoxDraw (ox->specBox, WHITE, BLUE) ;
      for (h = 0 ; h < arrayMax (oldhomols) ; h++)
        {temp = arrayp (oldhomols, h, HIDBOX) ;
         if (temp->hid == TRUE) graphBoxDraw (temp->box, WHITE, MAGENTA) ;
         else graphBoxDraw (temp->box, WHITE, BLUE) ;
	}
      }
    else if (ox->specBox == MIN_BOX)
      graphBoxDraw (ox->specBox, BLACK, WHITE) ;
    else if (ox->specBox >= MIN_BOX + 1 && ox->specBox < MIN_BOX + ox->spec_labels)
      graphBoxDraw (ox->specBox, BLUE, WHITE) ;
    }
 
  ox->specBox = box ;

  if (ox->specBox >= MIN_BOX + ox->spec_labels)
    {graphBoxDraw (ox->specBox, WHITE, GREEN) ;
     for (h = 0 ; h < arrayMax (ox->homols) ; h++)
       graphBoxDraw (arrayp(ox->homols, h, HIDBOX)->box, WHITE, RED) ;
    }
  else if (ox->specBox == MIN_BOX)
    graphBoxDraw (ox->specBox, WHITE, BLACK) ;
  else if (ox->specBox >= MIN_BOX+1 && ox->specBox < MIN_BOX + ox->spec_labels)
    graphBoxDraw (ox->specBox, WHITE, BLUE) ;

  graphCompScrollEntry (0, ox->spechom, 0, 0, 0, 0, 0) ; 
  graphEntryDisable () ;
  graphCompScrollEntry (0, ox->specloc, 0, 0, 0, 0, 0) ;
  graphEntryDisable () ;

  if (ox->message3Box)
   {
     if (ox->specBox >= MIN_BOX + ox->spec_labels)
       {for (i = 0 ; i < arrayMax (ox->specsegs) ; i++)
	 {seg = arrayp (ox->specsegs, i, SEG) ;
	 if (seg->box == ox->specBox)
	   {hit = seg->key ;
	   key1 = seg->gene1 ;
	   key2 = seg->gene2 ;
	   }
	 }
       strncpy (ox->message3Text, messprintf ("Homology "), 80) ; 
       strcat (ox->message3Text, messprintf ("%s ", name (hit))) ;
       strcat (ox->message3Text, messprintf ("between %s and %s", name (key1), name (key2))) ;
       }
     else
       *ox->message3Text = 0 ;
   }
  graphBoxDraw (ox->message3Box, BLACK, LIGHTBLUE) ;  
  arrayDestroy (oldhomols) ;

}


static void specPick (int box) 
{
  int specs ;

  OXGET ("specPick") ;

  specs = arrayMax (ox->specs) ;

  if (box >= MIN_BOX + ox->spec_labels &&
      box < (arrayMax (ox->specsegs) + MIN_BOX + ox->spec_labels)) 
    {if (box == ox->specBox)
       specFollow (ox) ;
     specSelect (ox, box) ;
    }
  else if (box == MIN_BOX)
    {if (box == ox->specBox)
       display (ox->c4, 0, 0) ;
     graphActivate (ox->specGraph) ;
     specSelect (ox, box) ;
    }
  else if (box >= MIN_BOX + 1 && box < MIN_BOX + 1 + specs)
    {if (box == ox->specBox)
       display (arrayp (ox->specs, box-MIN_BOX-1, SPEC)->key, 0, 0) ;
     graphActivate (ox->specGraph) ;
     specSelect (ox, box) ;
    }
  else if (box >= MIN_BOX + 1 + specs && box < MIN_BOX + ox->spec_labels)
    {if (box == ox->specBox)
       display (arrayp (ox->c4segs, box-MIN_BOX-1-specs, CHROMS)->key, 0, 0) ;
     graphActivate (ox->specGraph) ;
     specSelect (ox, box) ;
    }
  else if (box == FINDHOM_BOX)
    graphCompScrollEntry (0, ox->spechom, 0, 0, 0, 0, 0) ;
  else if (box == FINDLOC_BOX)
    graphCompScrollEntry (0, ox->specloc, 0, 0, 0, 0, 0) ;
        
 
          
}


static void specFollow (OX ox)
{
  int i ;
  KEY hit = 0 ;
  SEG *seg ;

  for (i = 0 ; i < arrayMax (ox->specsegs) ; i++)
    {seg = arrayp (ox->specsegs, i, SEG) ;
     if (seg->box == ox->specBox)
       hit = seg->key ;
    }
  if (hit)
    display (hit, 0, 0) ;
}


/**********************************/

int specConvert (OX ox)
{ 
  int found, h, i, j, k, l, m, p ;
  int h1, h2 ;
  float x1, x2, x3, x4 ;
  float bmin = 0, bmax = 0 ;
  float min1 = 0, max1 = 0 ;
  char *sym ;
  KEY band, col, group, key ;
  KEY m1, m2, map1, map2 ;
  KEY key1, key2 ;
  OBJ Band, Group, Locus1, Locus2 ;
  POS *p1, *p2 ;
  SEG *seg, *seg2 ;
  SPEC *mapset, *mapset2 ;
  OBJ chr1 = 0 ;
  CHROMS *c4seg ;
  Array bands, homol, homol2, posn1, posn2, tempspecs;
  KEYSET specs, maps, fuzzy, pairwise ;
  static BSMARK mark1 = 0 ;
  static BSMARK mark2 = 0 ;
  BOOL pos1, pos2, found2 ;

  if (!(chr1 = bsCreate (ox->specmap)))
    return -1 ;
  else
    {ox->c4 = ox->specmap ;
    }
  bsDestroy (chr1) ;

  ox->specflag |= FLAG_SPEC_CHROM ;
  ox->specflag &= ~FLAG_SPECBANDS ;
  ox->specwide = 75;
  ox->spechigh = 50 ;

  bands = arrayCreate (8, BSunit) ;
  homol = arrayCreate (32, BSunit) ;
  homol2 = arrayCreate (32, BSunit) ;
  ox->specsegs  = arrayCreate (8, SEG) ;
  ox->c4segs = arrayCreate (8, CHROMS) ;
  ox->specs = arrayCreate (8, SPEC) ;
  tempspecs = arrayCreate (8, SPEC) ;

  k = 0 ;
  specs = keySetCreate () ;
  specs = query(0, ">?Map_set ; Species") ;
  for (i = 0 ; i < keySetMax (specs) ; i++) 
    {found = 0 ;
     maps = keySetCreate () ;
     maps = query (0, messprintf (">?Map_set \"%s\" ; >Map", name (keySet (specs, i))));
     for (j = 0 ; j < keySetMax (maps) ; j++)
       if (chkSpecies (ox->c4, keySet (maps, j))) found++ ;
     if (found == 0)
      {mapset = arrayp (tempspecs, k++, SPEC) ;
       mapset->key = keySet (specs, i) ;
       mapset->homol = 0 ;
       mapset->maps = arrayCreate (8, KEY) ;
       for (l = 0 ; l < keySetMax (maps) ; l++)
         array (mapset->maps, l, KEY) = keySet (maps, l) ;
       }
     arrayDestroy (maps) ;
    }

 
  if ((chr1 = bsCreate (ox->c4)))
    {if (bsFindTag (chr1, _Chrom_Band) && bsFlatten (chr1, 1 , bands))
      {for (i = 0 ; i < arrayMax (bands) ; ++i )
        {band = array (bands, i, BSunit).k ;
         Band = bsCreate (band) ;
         if (bsFindTag (Band, _Map) && bsGetKey (Band, _bsRight, &key))
           if (bsPushObj (Band))
             if (bsGetData (Band, _Left, _Float, &x1) &&
                 bsGetData (Band, _Right, _Float, &x2))
                {c4seg = arrayp (ox->c4segs, i, CHROMS) ;
                 c4seg->key = band ;
                 c4seg->x = 0.5*(x1+x2) ;
                 c4seg->dx = x2 - x1 ;
	         c4seg->flag = 0 ;
                 ox->specflag |= FLAG_SPECBANDS ;
                 }
         bsDestroy (Band) ;
         if (x1 < bmin)
           bmin = x1 ;
         if (x2 > bmax)
	   bmax = x2 ;
         }

       for (i = 0 ; i < arrayMax (bands) ; ++i )
        {band = array (bands, i, BSunit).k ; 
         Band = bsCreate (band) ;
         c4seg = arrayp (ox->c4segs, i, CHROMS) ;
         if (bsFindTag (Band, _Centromere))
	    c4seg->flag |= FLAG_CENTROMERE ;
         if (bsFindTag (Band, _Dark))
	    c4seg->flag |= FLAG_DARK_BAND ;
	 if (bsFindTag (Band, _NOR))
	    c4seg->flag |= FLAG_NOR ;
	 if (bsFindTag (Band, _Colour) && bsGetKey (Band, _bsRight, &col))
	    c4seg->flag |= ((col - _WHITE) & 15  ) << 28 ;
         if (bsFindTag (Band, _Symbol)) 
           {bsGetData (Band, _bsRight, _Text, &sym) ;
            c4seg->sym = sym ;
           }
         bsDestroy (Band) ;
         }   
      }
     ox->c4min = bmin ;
     ox->c4max = bmax ; 
     bsDestroy (chr1) ;
    }   

  h = 0 ; 
  pairwise = keySetCreate () ;
  pairwise = query (0, messprintf (">?Map \"%s\" ; >Locus ; >Pairwise", name (ox->c4))) ;
    for (i = 0 ; i < keySetMax (pairwise) ; i++)
      {group = array (pairwise, i, KEY) ;
       Group = bsCreate (group) ;
       if (bsGetKey (Group, _Pairwise, &key1)) do
         {mark1 = bsMark (Group, mark1) ;
          bsFlatten (Group, 1, homol) ;
          Locus1 = bsCreate (key1) ;
          pos1 = FALSE ;
          if (bsFindTag (Locus1, _Map) && bsGetKey (Locus1, _bsRight, &map1))
           {if (bsPushObj (Locus1) && bsGetData (Locus1, _Position, _Float, &x1))
              {pos1 = TRUE ; x2 = 0 ;
               if (bsPushObj (Locus1)) bsGetData (Locus1, _Error, _Float, &x2) ;
              }
            for (j = 0 ; j < arrayMax (homol) ; j++)
              {key2 = array (homol, j, BSunit).k ;
               Locus2 = bsCreate (key2) ;
               pos2 = FALSE ;
               if (bsFindTag (Locus2, _Map) && bsGetKey (Locus2, _bsRight, &map2))
                 {if (bsPushObj (Locus2) && bsGetData (Locus2, _Position, _Float, &x3))
                   {pos2 = TRUE ; x4 = 0 ;
                    if (bsPushObj (Locus2)) bsGetData (Locus2, _Error, _Float, &x4) ;
                   }
                  for (k = 0 ; k < arrayMax (tempspecs) ; k ++)
                    {mapset = arrayp (tempspecs, k, SPEC) ;
                     for (l = 0 ; l < arrayMax (mapset->maps) ; l++)
                       {if (chkSpecies (map1, ox->c4) && chkSpecies (map2, array (mapset->maps, l, KEY)) && pos1 == TRUE)
                          {found2 = FALSE ;
                           for (p = 0 ; p < arrayMax (ox->specsegs) ; p++)
                             {seg2 = arrayp (ox->specsegs, p, SEG) ;
                              if (chkSpecies (key1, seg2->gene1) && chkSpecies (key2, seg2->gene2)) found2 = TRUE ;
                             }
                           if (found2 == FALSE) 
                             {seg = arrayp (ox->specsegs, h , SEG) ;
                              seg->key = group ;
                              seg->flag = 0 ;
                              seg->spec2 = mapset->key ;
                              seg->x1 = x1 ;
                              seg->dx1 = x2 ;
                              seg->gene1 = key1 ;
                              seg->gene2 = key2 ;
                              seg->map1 = map1 ;
                              seg->map2 = map2 ;
                              (mapset->homol)++ ;       
                              if (seg->x1 < min1) min1 = seg->x1 ;
                              if (seg->x1 > max1) max1 = seg->x1 ;
                              h++ ;
                             }
			  }
                        else if (chkSpecies (map2, ox->c4) && chkSpecies (map1, array (mapset->maps, l, KEY)) && pos2 == TRUE)
                          {found2 = FALSE ;
                           for (p = 0 ; p < arrayMax (ox->specsegs) ; p++)
                             {seg2 = arrayp (ox->specsegs, p, SEG) ;
                              if (chkSpecies (key1, seg2->gene1) && chkSpecies (key2, seg2->gene2)) found2 = TRUE ;
                             }
                           if (found2 == FALSE) 
                             {seg = arrayp (ox->specsegs, h , SEG) ;
                              seg->key = group ;
                              seg->flag = 0 ;
                              seg->spec2 = mapset->key ;
                              seg->x1 = x3 ;
                              seg->dx1 = x4 ;
                              seg->gene1 = key2 ;
                              seg->gene2 = key1 ;
                              seg->map1 = map2 ;
                              seg->map2 = map1 ;
                              (mapset->homol)++ ;
                              if (seg->x1 < min1) min1 = seg->x1 ;
                              if (seg->x1 > max1) max1 = seg->x1 ;
                              h++ ;
                             }
			  }     
                       }
		    }
		 }
	      }
	   }
         bsDestroy (Locus1) ;
        }while (bsGoto (Group, mark1), bsGetKey (Group, _bsDown, &key1)) ;
      bsDestroy (Group) ;
     }

  fuzzy = keySetCreate () ;
  fuzzy = query (0, messprintf (">?Map \"%s\" ; >Locus ; >Fuzzy", name (ox->c4))) ;
    for (i = 0 ; i < keySetMax (fuzzy) ; i++)
    {group = array (fuzzy, i, KEY) ;
     Group = bsCreate (group) ;
       if (bsFindTag (Group, _Fuzzy)) 
         {bsFlatten (Group, 1, homol2) ;
          for (h1 = 0 ; h1 < arrayMax (homol2) - 1 ; h1++)
          for (h2 = h1 + 1 ; h2 < arrayMax (homol2) ; h2++)
	    {key1 = array (homol2, h1, BSunit).k ;
             key2 = array (homol2, h2, BSunit).k ;
             posn1 = arrayCreate (8, POS) ; 
             j = 0 ;
             Locus1 = bsCreate (key1) ;
             if (bsGetKey (Locus1, _Map, &m1)) do
               {mark1 = bsMark (Locus1, mark1) ;
                p1 = arrayp (posn1, j, POS) ;
                p1->pos = FALSE ;
                if (bsPushObj (Locus1) && bsGetData (Locus1,_Position,_Float, &x1))
                  {p1->map = m1 ;
                   p1->x1 = x1 ;
                   p1->pos = TRUE ;
                   p1->dx1 = 0 ;
                   if (bsPushObj (Locus1) && bsGetData (Locus1,_Error, _Float, &x2))
                      p1->dx1 = x2 ;
                  }
                j++ ;
               }while (bsGoto (Locus1, mark1), bsGetKey (Locus1, _bsDown, &m1)) ; 
             posn2 = arrayCreate (8, POS) ;
             j = 0 ;
             Locus2 = bsCreate (key2) ;
             if (bsGetKey (Locus2, _Map, &m2)) do
               {mark2 = bsMark (Locus2, mark2) ;
                p2 = arrayp (posn2, j, POS) ;
                p2->pos = FALSE ;
                if (bsPushObj (Locus2) && bsGetData (Locus2,_Position,_Float, &x3))
                  {p2->map = m2 ;
                   p2->x1 = x3 ;
                   p2->pos = TRUE ;
                   p2->dx1 = 0 ;
                   if (bsPushObj (Locus2) && bsGetData (Locus2,_Error, _Float, &x4))
                      p2->dx1 = x4 ;
                  }
               }while (bsGoto (Locus2, mark2), bsGetKey (Locus2, _bsDown, &m2)) ;

             for (j = 0 ; j < arrayMax (posn1) ; j++)
               for (k = 0 ; k < arrayMax (posn2) ; k++)
                 {map1 = arrayp (posn1, j, POS)->map ;
                  map2 = arrayp (posn2, k, POS)->map ;
                  for (l = 0 ; l < arrayMax (tempspecs) ; l++)
                    {mapset = arrayp (tempspecs, l, SPEC) ;
                     for (m = 0 ; m < arrayMax (mapset->maps) ; m++)
                       {if (chkSpecies (map1, ox->c4) && chkSpecies (map2, array (mapset->maps, m, KEY)) && !(chkSpecies (key1, key2)) && arrayp (posn1, j, POS)->pos == TRUE)
                          {found2 = FALSE ;
                           for (p = 0 ; p < arrayMax (ox->specsegs) ; p++)
                             {seg2 = arrayp (ox->specsegs, p, SEG) ;
                              if (chkSpecies (key1, seg2->gene1) && chkSpecies (key2, seg2->gene2)) found2 = TRUE ;
                             }
                           if (found2 == FALSE) 
                             {seg = arrayp (ox->specsegs, h , SEG) ;
                              seg->key = group ;
                              seg->flag = 0 ;
                              seg->spec2 = mapset->key ;
                              seg->x1 = arrayp (posn1, j, POS)->x1 ;
                              seg->dx1 = arrayp (posn1, j, POS)->dx1 ;
                              seg->gene1 = key1 ;
                              seg->gene2 = key2 ;
                              seg->map1 = map1 ;
                              seg->map2 = map2 ;
                              (mapset->homol)++ ;
                              if (seg->x1 < min1) min1 = seg->x1 ;
                              if (seg->x1 > max1) max1 = seg->x1 ;
                              h++ ;
                             }
			  }
                        else if (chkSpecies (map2, ox->c4) && chkSpecies (map1, array (mapset->maps, m, KEY)) && !(chkSpecies (key1, key2)) && arrayp (posn2, k, POS)->pos == TRUE)
                          {found2 = FALSE ;
                           for (p = 0 ; p < arrayMax (ox->specsegs) ; p++)
                             {seg2 = arrayp (ox->specsegs, p, SEG) ;
                              if (chkSpecies (key1, seg2->gene1) && chkSpecies (key2, seg2->gene2)) found2 = TRUE ;
                             }
                           if (found2 == FALSE) 
                             {seg = arrayp (ox->specsegs, h , SEG) ;
                              seg->key = group ;
                              seg->flag = 0 ;
                              seg->spec2 = mapset->key ;
                              seg->x1 = arrayp (posn2, k, POS)->x1 ;
                              seg->dx1 = arrayp (posn2, k, POS)->dx1 ;
                              seg->gene1 = key2 ;
                              seg->gene2 = key1 ;
                              seg->map1 = map2 ;
                              seg->map2 = map1 ;
                              (mapset->homol)++ ;
                              if (seg->x1 < min1) min1 = seg->x1 ;
                              if (seg->x1 > max1) max1 = seg->x1 ;
                              h++ ;
                             }
                          }   
		       }
		    } 
		 } 
             bsDestroy (Locus1) ;
             bsDestroy (Locus2) ;
             arrayDestroy (posn1) ;
             arrayDestroy (posn2) ;
            }
         }
     bsDestroy (Group) ;
    }

  if (!(ox->specflag & FLAG_SPECBANDS))
    {ox->c4min = min1 - 5 ;
     ox->c4max = max1 + 5 ;
     }

  /*  ox->specs = arrayCopy (tempspecs) ;
  for (k = 0 ; k < arrayMax (ox->specs) ; k ++)
    {mapset = arrayp (ox->specs, k, SPEC) ;
     messout ("%d homologies found", mapset->homol) ;
    }*/
  
  
  k = 0 ;
  for (i = 0 ; i < arrayMax (tempspecs) ; i++)
    {mapset = arrayp (tempspecs, i, SPEC) ;
     if (mapset->homol > 0) 
       {mapset2 = arrayp(ox->specs, k, SPEC) ;
        mapset2->key = mapset->key ;
        mapset2->homol = mapset->homol ;
        mapset2->maps = arrayCopy(mapset->maps) ;
        k++ ;
       }
    }
    

  arrayDestroy (specs) ;
  arrayDestroy (fuzzy) ;
  arrayDestroy (pairwise) ;
  arrayDestroy (bands) ;
  arrayDestroy (homol) ;
  arrayDestroy (homol2) ;
  arrayDestroy (tempspecs) ;
  ox->s = h ;
  return h ;

}

static void specDraw (OX ox)
{
  int box, i, j, k, bands, specs, xmax, ymax ;
  float w1, w2, h1, h2 ;
  float xl, yl, xlen = 0, ylen = 0 ;
  float htot, wtot, xscale, yscale ;
  CHROMS *c4seg ;
  
  graphClear () ;
  graphText (messprintf ("Species Grid of map %s with all other entered map_sets : %d homologies found" , name (ox->c4), ox->s), 1, 1) ;
  
  if (ox->specflag & FLAG_SPECBANDS)
    {if (ox->specflag & FLAG_SPEC_CHROM)
      {if (graphButton ("Hide Chromosomes", toggleChroms, 1, 3) != CHROM_BOX)
       messcrash ("box screwup at %d in specDraw", CHROM_BOX) ;}
    else
      {if (graphButton ("Draw Chromosomes", toggleChroms, 1, 3) != CHROM_BOX)
       messcrash ("box screwup at %d in specDraw", CHROM_BOX) ;}
     }
  else
    {if (graphButton ("Draw Chromosomes", noChroms, 1, 3) != CHROM_BOX)
      messcrash ("box screwup at %d in specDraw", CHROM_BOX) ;}
  
  if (graphButton ("Highlight...", specClear, 19, 3) != HIDE_BOX)
    messcrash ("box screwup at %d in specDraw", HIDE_BOX) ;

  graphBoxMenu (HIDE_BOX, Highlight3Menu) ;  

  if (graphButton ("Zoom in", specZoomIn, 33, 3) != ZOOM_IN)
    messcrash ("box screwup at %d in specDraw", ZOOM_IN) ;

  if (graphButton ("Zoom out", specZoomOut, 42, 3) != ZOOM_OUT)
    messcrash ("box screwup at %d in specDraw", ZOOM_OUT) ;
  
  graphText ("Find homology :", 1, 5) ;
  *ox->spechom = 0 ;
  if (graphCompScrollEntry (spechomolCompletion, ox->spechom, 20, 20, 17, 5, specFindHom) != FINDHOM_BOX)
    messcrash ("box screwup at %d in specDraw", FINDHOM_BOX) ;

  graphEntryDisable () ;

  graphText ("Find locus :", 42, 5) ;
  *ox->specloc = 0 ;
  if (graphCompScrollEntry (speclocusCompletion, ox->specloc, 20, 20, 55, 5, specFindLoc) != FINDLOC_BOX)
    messcrash ("box screwup at %d in specDraw", FINDLOC_BOX) ;

  graphEntryDisable () ;
  
  ox->message3Box = MESS_BOX ;
 *ox->message3Text = 0 ;
  ox->message3Box = graphBoxStart () ;
  graphTextPtr (ox->message3Text, 1, 7, 150) ;
  graphBoxEnd () ;
  graphBoxDraw (ox->message3Box, BLACK, LIGHTBLUE) ;

  xOff = 0 ; yOff = 0 ;
  xOff += 2 + strlen (name (ox->c4)) ;
  graphLine (0, 9.5, 200, 9.5) ; 
  yOff += 10 ;

  xmax = arrayMax (ox->specs) - 1 ; 
  specs = xmax + 1 ;
  ymax = arrayMax (ox->c4segs) - 1 ;
  bands = ymax + 1 ;
  wtot = specs;
  htot = ox->c4max - ox->c4min ; 
  yscale = ox->spechigh/htot ;
  xscale = ox->specwide/wtot ;
  /*messout ("xmax = %d, specs = %d, ymax = %d, bands = %d", xmax, specs, ymax, bands);
  messout ("wtot = %f, htot = %f, yscale = %f, xscale = %f", wtot, htot, yscale, xscale);*/                            
  graphColor (BLACK) ; 

  for (i = 0 ; i < keySetMax (ox->specs) ; i++)
   {yl = strlen (name (keySet (ox->specs, i))) ;
    if (yl > ylen) ylen = yl ;
    }
  ylen *= 0.65 ;
  yOff += 5 + ylen ;

  if (ox->specflag & FLAG_SPECBANDS)
   {for (i = 0 ; i < arrayMax (ox->c4segs) ; i++)
      {c4seg = arrayp (ox->c4segs, i, CHROMS) ;
       xl = strlen (name (c4seg->key)) ;
       if (xl > xlen) xlen = xl ;
       }
    xOff += 8 + xlen ;
    }
  else
    xOff += 10 ;
  
  box = MIN_BOX ;
  if (box != graphBoxStart ())
    messcrash ("minLiveBox wrong in specDraw ()") ;
  graphTextFormat (BOLD) ;
  graphText (name (ox->c4), 1, yOff + (ox->spechigh)/2) ;
  graphBoxEnd () ;
  graphBoxDraw (box, BLACK, WHITE) ;
  box++ ;
  graphTextFormat(PLAIN_FORMAT);

  for (i = 0 ; i <= xmax ; i++)
   {w1 = 0 ; w2 = 0 ;
    w2 = (i+1) ;
    w1 = i ;
    if (box != graphBoxStart ())
      messcrash ("minLiveBox wrong in specDraw ()") ;
    graphTextHeight (0.5) ; 
    japanese (name (arrayp (ox->specs, i, SPEC)->key), xOff+ (w1+w2)*0.5*xscale, yOff-3-ylen) ; 
    graphBoxEnd () ;
    graphBoxDraw (box, BLUE, WHITE) ;
    box ++ ;    
   }

  graphTextHeight (1) ; 
  
  if (ox->specflag & FLAG_SPECBANDS)
   {for (i = 0 ; i <= xmax ; i++)
     {for ( j = 0 ; j <= ymax ; j++)
       {w1 = 0 ; h1 = 0 ; w2 = 0 ; h2 = 0 ;
        w2 = (i+1) ;
        w1 = i ;
        for ( k = 0 ; k <= j ; k++) {h2 += (arrayp(ox->c4segs, k, CHROMS)->dx) ;}
        for ( k = 1 ; k <= j ; k++) {h1 += (arrayp (ox->c4segs, k - 1, CHROMS)->dx) ;}
        
        graphRectangle (xOff+(w1*xscale), yOff+(h1*yscale), xOff+(w2*xscale), yOff+(h2*yscale)) ;
        }
      } 
    }
  else
   {for (i = 0 ; i <= xmax ; i++)
     {w1 = 0 ; w2 = 0 ;
      w2 = (i+1) ;
      w1 = i ;
      graphRectangle (xOff+(w1*xscale), yOff+ox->spechigh, xOff+(w2*xscale), yOff) ;
      }
    }
  
  if (ox->specflag & FLAG_SPECBANDS)
   {if (ox->specflag & FLAG_SPEC_CHROM)
     {ox->spec_labels = 1 + specs + bands ;
      drawChromosome (ox, xlen) ;
     }
    else
     {ox->spec_labels = 1 + specs ;
      drawScale (ox) ;
     }
   }
  else
    {ox->spec_labels = 1 + specs ;
     drawScale (ox) ;
    }
  
  drawSpots (ox, xscale) ;
  drawLines (ox, xscale) ;
  
  ox->specBox = 0 ;
  graphRedraw () ;

}

static void drawChromosome (OX ox, float xlen)
{
  float y0, dy, y, yend ;
  int i, box, specs ;
  CHROMS *c4seg ;
 
  y0 = -(ox->c4min) ;
  yend = (ox->c4max - ox->c4min) / ox->spechigh ;
  specs = arrayMax (ox->specs) ;
  
  box = MIN_BOX + 1 + specs ;
  for (i = 0 ; i < arrayMax (ox->c4segs) ; ++i)
    {c4seg = arrayp (ox->c4segs, i, CHROMS) ;
     if (class (c4seg->key) == _VChrom_Band)
       {
	y = yOff + ( y0 + c4seg->x) / yend ;
	dy = ( 0.5 * c4seg->dx) / yend ;
	if (c4seg->flag & FLAG_CENTROMERE)
	  {graphColor (WHITE) ;
	   graphFillRectangle (xOff-4, y-dy, xOff, y+dy) ;
	   graphColor (BLACK) ;
	   graphLine (xOff-4, y-dy, xOff, y+dy) ;
	   graphLine (xOff, y-dy, xOff-4, y+dy) ;
	   graphLine (xOff-4, y-dy, xOff, y-dy) ;
	   graphLine (xOff, y+dy, xOff-4, y+dy) ;
           if (box != graphBoxStart ())
           messcrash ("minLiveBox wrong in specDraw () box %d", box) ;
           graphTextHeight (0.5) ;
           graphText (name (c4seg->key), xOff-6-xlen, y-0.25) ;
           graphBoxEnd () ;
           graphBoxDraw (box, BLUE, WHITE) ;
           box++ ; 
           graphTextHeight (1) ;
     	   }	      
	else 
	  {if (c4seg->flag & FLAG_DARK_BAND)
	     graphColor (DARKGRAY) ;
	   else if (c4seg->flag & FLAG_NOR)
	     graphColor (LIGHTGRAY) ;
	   else
	     graphColor (WHITE) ;
	   if (c4seg->flag >> 28)
	     graphColor((c4seg->flag  >> 28) + WHITE) ;
	   graphFillRectangle (xOff-4, y-dy, xOff, y+dy) ;
	   graphColor (BLACK) ;
	   graphRectangle (xOff-4, y-dy, xOff, y+dy) ;
           if (box != graphBoxStart ())
           messcrash ("minLiveBox wrong in specDraw () box %d", box) ;
           graphTextHeight (0.5) ;
           graphText (name (c4seg->key), xOff-6-xlen, y-0.25) ;
           graphBoxEnd () ;
           graphBoxDraw (box, BLUE, WHITE) ;
           box++ ; 
           graphTextHeight (1) ;
	   }
        }
     }
}

static void drawScale (OX ox) 
{
  int i, j, k ;
  float unit, y0 ;

/* If no banding data or if the chromosomes are hidden, draw a scale instead */

  graphTextHeight (0.5) ;
  y0 = yOff + ox->spechigh * (1-(ox->c4max / (ox->c4max-ox->c4min))) ;
  unit = ox->spechigh / (ox->c4max - ox->c4min) ;
  for (i = 0 ; i < (ox->c4max + 1) ; i++)
    graphLine (xOff-3, y0+i*unit, xOff-3.9, y0+i*unit) ;
  for (k = 0 ; k < (ox->c4max + 1) ; k+=10)
   {graphLine (xOff-3, y0+k*unit, xOff-4.5, y0+k*unit) ;
    graphText (messprintf ("%3d",k), xOff - 8, y0+k*unit-0.5) ;
    }
  for (j = 0 ; j < -(ox->c4min - 1) ; j++)
    graphLine (xOff-3, y0-j*unit, xOff-3.9, y0-j*unit) ;
  for (k = 0 ; k < -(ox->c4min - 1) ; k+=10)
   {graphLine (xOff-3, y0-k*unit, xOff-4.5, y0-k*unit) ;
    graphText (messprintf ("%3d",-k), xOff - 8, y0-k*unit-0.5) ;
    }  
  graphLine (xOff-3, y0-(j-1)*unit, xOff-3, y0+(i-1)*unit) ;
  graphTextHeight (1) ;
 
}


static void drawSpots (OX ox, float xscale)
{
  int i, j, s ;
  float x=0, y, y0, yend ;
  float w1, w2 ;
  SEG *seg ;

  y0 = -(ox->c4min) ;
  yend = (ox->c4max - ox->c4min) / ox->spechigh ;
  graphTextHeight (0.5) ;
  s = MIN_BOX + ox->spec_labels ;
  for (i = 0 ; i < arrayMax (ox->specsegs) ; i++)
    {seg = arrayp (ox->specsegs, i, SEG) ;
     if (!(seg->flag & FLAG_HIDE))
      {for (j = 0 ; j < arrayMax (ox->specs) ; j++)
         {if (chkSpecies (seg->spec2, arrayp (ox->specs, j, SPEC)->key))
          x = j ;
          }
       w2 = x+1 ;
       w1 = x ;
       x = xOff + (0.5*(w2+w1)*xscale) ;
       y = yOff + (y0 + seg->x1) / yend ; 
       if (s != graphBoxStart ())
         messcrash ("problems with pairMap point %d", s) ;
       graphRectangle (x-0.4, y-0.4, x+0.4, y+0.4) ;
       graphBoxEnd () ;
       if (seg->flag & FLAG_LIGHT) graphBoxDraw (s, WHITE, MAGENTA) ;
         else graphBoxDraw (s, WHITE, BLUE) ;
       graphText (messprintf ("%s", name (seg->map2)), x + 1, y - 0.4) ;
       seg->box = s ;
       s++ ;
      }
    }
  graphTextHeight (1) ;
}


int pointOrder (const void *a, const void *b)
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
  int i, j, k ;
  float x, y0, yend ;
  SPEC *spec ;
  SEG *seg ;
  POINT *point, *point2 ;
  Array tpoints ;

  graphColor (MAGENTA) ;
  y0 = -(ox->c4min) ;
  yend = (ox->c4max - ox->c4min) / ox->spechigh ;
  for (i = 0 ; i < arrayMax (ox->specs) ; i++)
    {spec = arrayp (ox->specs, i, SPEC) ;
     tpoints = arrayCreate (64, POINT) ;
     k = 0 ;
     for (j = 0 ; j < arrayMax (ox->specsegs) ; j++)
       {seg = arrayp (ox->specsegs, j, SEG) ;
        if (chkSpecies (seg->spec2, spec->key))
	    {point = arrayp (tpoints, k, POINT) ;
             point->y = yOff + (y0 + seg->x1) / yend ;
             point->map = seg->map2 ;
             k++ ;
             }
	   
       }
     arraySort (tpoints, pointOrder) ;
     for (k = 0 ; k < arrayMax (tpoints) - 1 ; k++)
       {point = arrayp (tpoints, k, POINT) ;
        point2 = arrayp (tpoints, k+1, POINT) ;  
        /*messout ("Column %d : map1 = %s, pos1 = %f and map2 = %s, pos2 = %f", i+1, name(point->map), point->y, name(point2->map), point2->y) ; */    
        x = xOff + (0.5*(2*i + 1)*xscale) ;
        if (chkSpecies (point->map, point2->map))
          graphLine (x, point->y, x, point2->y) ;
       }
     arrayDestroy (tpoints) ;
    }
  graphColor (BLACK) ;
}


extern int ksetClassComplete (char *text, int len, int classe) ;

static int spechomolCompletion (char *cp, int len)
{
  return ksetClassComplete (cp, len, _VHomology_group) ;
}

static int speclocusCompletion (char *cp, int len)
{
  return ksetClassComplete (cp, len, _VLocus) ;
}

static void specFindHom (char* string)
{
  int i ;
  KEY key ;
  SEG *seg ;
  BOOL found = FALSE ;

  OXGET ("specFindHom") ;

  if (!(lexword2key (string, &key, _VHomology_group)))
    {messout ("Sorry, %s is not a homology name", string) ;
     return ;
    }

  for (i = 0 ; i < arrayMax (ox->specsegs) ; i++)
    {seg = arrayp (ox->specsegs, i, SEG) ;
     if (chkSpecies (seg->key, key)) 
       {found = TRUE ;
        specSelect (ox, seg->box) ;
       }
     }

  if (found == FALSE) 
   {messout ("Sorry, this homology is not in this grid") ;
    return ;
   }

}

static void specFindLoc (char* string)
{
  int i ;
  KEY key ;
  SEG *seg ;
  BOOL found = FALSE ;

  OXGET ("specFindLoc") ;

  if (!(lexword2key (string, &key, _VLocus)))
    {messout ("Sorry, %s is not a locus name", string) ;
     return ;
    }

  for (i = 0 ; i < arrayMax (ox->specsegs) ; i++)
    {seg = arrayp (ox->specsegs, i, SEG) ;
     if (chkSpecies (seg->gene1, key) || chkSpecies (seg->gene2, key)) 
       {found = TRUE ;
        specSelect (ox, seg->box) ;
       }
     }

  if (found == FALSE) 
   {messout ("Sorry, this locus is not in this grid") ;
    return ;
   }

}

static void toggleChroms (void) 
{
  OXGET ("toggleChroms") ;
  
  if (!(ox->specflag & FLAG_SPECBANDS))
    {messout ("Sorry, there's no band data") ;
     return ;}
  else
    {ox->specflag ^= FLAG_SPEC_CHROM ;
     specDraw (ox) ;}

}

static void noChroms (void) 
{
  OXGET ("noChroms") ;

  messout ("Sorry, there's no band data") ;
  return ;

}

static void specClear (void) 
{
  int i ;
  SEG *seg ;

  OXGET ("specClear") ;

  seg = arrayp (ox->specsegs, 0, SEG) ;
  i = arrayMax (ox->specsegs) ;
  while(i--)
    {seg->flag &= ~FLAG_HIDE ;  /* bitwise NOT */
     seg->flag &= ~FLAG_LIGHT ;
     seg++ ;
    }

  specDraw (ox) ;

}

static void specLightAll (void) 
{
  int i ;

  OXGET ("specLightAll") ;

  for (i = 0 ; i < arrayMax (ox->specsegs) ; ++i)
    arrp (ox->specsegs, i, SEG)->flag |= FLAG_LIGHT ;
    
  ox->gridBox = 0 ;

  specDraw (ox) ;


}

static void specLightHom (void) 
{
  int i ;
  SEG *seg ;

  OXGET ("specLightHom") ;

  if (!(ox->specBox))
    {messout ("You have not selected an homology") ;
     return ;}

  for (i = 0 ; i < arrayMax (ox->specsegs) ; i++)
    {seg = arrayp (ox->specsegs, i, SEG) ;
     if (seg->box == ox->specBox) seg->flag |= FLAG_LIGHT ;
    }

  specDraw (ox) ;

}

static void specUnLightHom (void) 
{
  int i ;
  SEG *seg ;

  OXGET ("specUnLightHom") ;

  if (!(ox->specBox))
    {messout ("You have not selected an homology") ;
     return ;}

  for (i = 0 ; i < arrayMax (ox->specsegs) ; i++)
    {seg = arrayp (ox->specsegs, i, SEG) ;
     if (seg->box == ox->specBox) seg->flag &= ~FLAG_LIGHT ;
    }

  specDraw (ox) ;

}

static void specLight (void) 
{
  int i, j ;
  SEG *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;

  OXGET ("specLight") ;

  if (!keySetActive (&keySet, &dummy))
    { messout ("First select a keySet window, thank you.") ;
      return ;
    }
  for (i = 0 ; i < arrayMax (ox->specsegs) ; ++i)
    { seg = arrp (ox->specsegs, i, SEG) ;
      if (keySetFind (keySet, seg->key, &j) ||
          keySetFind (keySet, seg->gene1, &j) ||
          keySetFind (keySet, seg->gene2, &j)) 
	seg->flag |= FLAG_LIGHT ;
    }

  specDraw (ox) ;


}

static void specUnLight (void) 
{
  int i, j ;
  SEG *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;

  OXGET ("specUnLight") ;

  if (!keySetActive (&keySet, &dummy))
    { messout("First select a keySet window, thank you.") ;
      return ;
    }
  for (i = 0 ; i < arrayMax (ox->specsegs) ; ++i)
    { seg = arrp (ox->specsegs, i, SEG) ;
      if (keySetFind (keySet, seg->key, &j) ||
          keySetFind (keySet, seg->gene1, &j) ||
          keySetFind (keySet, seg->gene2, &j)) 
	seg->flag &= ~FLAG_LIGHT ;
    }

  specDraw (ox) ;

}

static void specHide (void) 
{
  int i, j ;
  SEG *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;

  OXGET ("specHide") ;

  if (!keySetActive (&keySet, &dummy))
    {messout ("First select a keySet window, thank you.") ;
     return ;
    }

  seg = arrp (ox->specsegs, 0, SEG) ;
  i = arrayMax (ox->specsegs) ;
  while(i--)
    {seg->flag &= ~FLAG_HIDE ; 
     if (!keySetFind (keySet, seg->key, &j) &&
         !keySetFind (keySet, seg->gene1, &j) &&
         !keySetFind (keySet, seg->gene2, &j)) 
       seg->flag |= FLAG_HIDE ;
     seg++ ;
    }

  specDraw (ox) ;


}

static void specHideLight (void) 
{
  int i ;
  SEG *seg ;

  OXGET ("specHideLight") ;

  for (i = 0 ; i < arrayMax (ox->specsegs) ; ++i)
    { seg = arrp (ox->specsegs, i, SEG) ;
      if (seg->flag & FLAG_LIGHT)
	seg->flag |= FLAG_HIDE ;
    }
  specDraw (ox) ;


}

static void specUnHide (void) 
{
  int i ;
  SEG *seg ;

  OXGET ("specUnHide") ;

  for (i = 0 ; i < arrayMax (ox->specsegs) ; ++i)
    { seg = arrp (ox->specsegs, i, SEG) ;
      seg->flag &= ~FLAG_HIDE ;
    }
  specDraw (ox) ;

}

static void specExport (void) 
{
  int i, j, k = 0 ;
  SEG *seg ;
  static KEYSET kset = 0 ;
  BOOL found = FALSE ;

  OXGET ("specExport") ;

  kset = keySetReCreate (kset) ;

  for (i = 0 ; i < arrayMax (ox->specsegs) ; ++i)
    { found = FALSE ;
      seg = arrp (ox->specsegs, i, SEG) ;
      if (seg->flag & FLAG_LIGHT)
        {for (j = 0 ; j < keySetMax (kset) ; j++)
          {if (chkSpecies (seg->key, keySet (kset, j))) found = TRUE ;}
         if (found == FALSE) keySet (kset, k++) = seg->key ;
        }
     }
  keySetSort (kset) ;
  keySetCompress (kset) ;
  keySetNewDisplay (kset, "Ox") ;
}

static void specZoomOut (void) 
{
  OXGET ("specZoomOut") ;

  if (ox->specwide >= 30 && ox->spechigh >= 20)
    {ox->specwide -= 15 ;
     ox->spechigh -= 10 ;

     specDraw (ox) ;
    }

}

static void specZoomIn (void) 
{
  OXGET ("specZoomIn") ;

  if (ox->specwide <= 120 && ox->spechigh <= 80)
    {ox->specwide += 15 ;
     ox->spechigh += 10 ;

     specDraw (ox) ;
    }

}


/************** end of file ****************/




 
 
 
 
 
