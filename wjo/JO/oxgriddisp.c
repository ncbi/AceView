/*  File: oxgriddisp.c
 *  Written by Jo Dicks (jld@bioch.ox.ac.uk)
 *-------------------------------------------------------------------
 * This file is part of the ACeDB comparative mapping package, 
 * written by Jo Dicks and Michelle Kirby
 * The original ACeDB system was written by :
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *      Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *
 * Description:  to display Oxford Grids
 * Exported functions: oxgridCreate() and oxgridDisplay() 
 * HISTORY:
 * Last edited: Mar 19 16:59 1997 (srk)
 *	- change NULL's to 0's to please MS Visual C++ compiler
 * Created: July 1993 
 * New version using Homology_group and Map_set models
 * Last edited: Dec  5 15:23 1996 (srk)
 *-------------------------------------------------------------------
 */

/* $Id: oxgriddisp.c,v 1.9 1997/11/10 09:15:13 mieg Exp $ */

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
#include "pick.h"
#include "bindex.h"

int OX_MAGIC = 9999 ;
int phenotype_cat = 0 ;

static int oxgridConvert (OX ox) ; 
void oxgridDefine (void) ;
int getMapsetList (void) ;
void getMapList (void) ;
static int homolCompletion (char *cp, int len) ;
static int locusCompletion (char *cp, int len) ;
static void oxgridDestroy (void) ;   
static void oxgridDraw (OX ox) ;
static void oxgridPick (int box) ;
static void oxgridKeyboard (int k) ;
static void oxgridSelect (OX ox, int box) ;
static void oxgridFlip (void) ;
static void oxgridDump(void) ;
static void oxgridReorg (void) ;
static void toggleCells (void) ;

static void oxClear (void) ;
static void oxLightAll (void) ;
static void oxLight (void) ;
static void oxUnLight (void) ;
static void oxHide (void) ;
static void oxHideLight (void) ;
static void oxUnHide (void) ;
static void oxExport (void) ;

static void oxClear2 (void) ;
static void oxEye1 (void) ;
static void oxEye2 (void) ;
static void oxEye3 (void) ;
static void oxEye4 (void) ;
static void oxEye5 (void) ;
static void oxEye6 (void) ;
static void oxEye7 (void) ;
static void oxEye8 (void) ;
static void oxEye9 (void) ;
static void oxEye10 (void) ;
static void oxEye11 (void) ;

static void choosePick (int box) ;
static void oxDefineSpecies1 (KEY k, int box) ;
static void oxDefineSpecies2 (KEY k, int box) ;
static void oxDefineSearch (KEY k, int box) ;
static void oxDefineMap (KEY k, int box) ;
static void oxLightCell (void) ;
static void oxUnLightCell (void) ;
static void oxFindHom (char *string) ;
static void oxFindLoc (char *string) ;
static void oxgridZoomIn (void) ;
static void oxgridZoomOut (void) ;
static void oxgridFish (void) ;
static void oxgridRandom (void) ;
static FILE* f;

#define graphBoxBox(_box) { \
	       float _x1, _y1, _x2, _y2 ; \
	       graphBoxDim (_box, &_x1, &_y1, &_x2, &_y2) ; \
	       graphRectangle (_x1 - 0.4, _y1 - 0.1, _x2 + 0.4, _y2 + 0.1) ; \
	      }

enum BoxNames { BACKGROUND=0, 
                EQUAL_BOX, FLIP_BOX, HIDE_BOX, 
                REORG_BOX, FINDHOM_BOX, FINDHOM2_BOX,
                FINDLOC_BOX, FINDLOC2_BOX,
                FISH_BOX, RANDOM_BOX,
                PHENOTYPE_BOX,
                ZOOM_IN, ZOOM_OUT,
                MIN_LIVE_BOX } ;               /* MIN_LIVE_BOX must be last */


static MENUOPT oxgridMenu[] =
            { {graphDestroy, "Quit"},
              {pairMapCreate, "Pairmap"},
              {o2mCreate, "One-to-Many"},
              {oxLightCell, "Highlight cell"},
              {oxUnLightCell, "Unhighlight cell"},
	      {help, "Help"},
	      {graphPrint, "Print"},
	      {displayPreserve, "Preserve"},
	      {oxgridDump, "Dump"},
	       {0, 0}
            } ;

static MENUOPT LightMenu[] =
            { {oxClear, "Clear and show all"},
              {oxLightAll, "Highlight all"},
              {oxLight, "Highlight selected keyset"},
              {oxUnLight,"Unhighlight selected keyset"},
              {oxHide, "Restrict to keyset"},
              {oxHideLight, "Hide highlit items"},
              {oxUnHide, "Unhide hidden items"},
              {oxExport, "Export highlit keyset"},
              {0, 0}
            } ;

static MENUOPT PhenotypeMenu[] =
            { {oxClear2, "Clear"},
              {oxEye1, "Eye"},
              {oxEye2, "Retina"},
              {oxEye3, "RP"},
              {oxEye4, "Lens"},
              {oxEye5, "Cataract"},
              {oxEye6, "Coloboma"},
              {oxEye7, "Global malformation"},
              {oxEye8, "Localised to eye"},
              {oxEye9, "Generalised"},
              {oxEye10, "Anophthalmos, unilateral"},
              {oxEye11, "Anophthalmos"},
              {0, 0}
            } ;

static MENUOPT StatsMenu[] =
            { {oxgridReorg, "Compute q"},
	       {0, 0}
            } ;


static MENUOPT chooseMenu[] =
            { {graphDestroy, "Quit"},
              {help, "Help"},
	      {graphPrint, "Print"},
	      {displayPreserve, "Preserve"},
              {oxgridDisplay, "OxGrid"},
              {specDisplay, "SpecGrid"},
	       {0, 0}
            } ;

FREEOPT searchChoice[] =
{
  {2, "Normal"},
  {'s', "Include self-homology"},
   {'n', "Normal"}
} ;

static FREEOPT *mapsetList = 0 ;
static Array mapsetListArray = 0 ;
static FREEOPT *mapList = 0 ;
static Array mapListArray = 0 ;



/*********************************/

BOOL oxgridPossible (void)
{
  int n = pickWord2Class("Map_set") ;
  return n && lexMax(n) > 2 ;
}

void oxgridCreate (void)
{ 
 int m ;
 OX ox ;

/* Initialise the OX structure and create the Comparative Map Chooser */

  ox = (OX) messalloc (sizeof (struct OXSTUFF)) ;
  ox->magic = &OX_MAGIC ;
  ox->flag  = 0 ;
  ox->flag &= ~FLAG_SELF ;
 
  m = getMapsetList () ;
  if (m == 0) {
    messout ("Sorry, there are no Map_sets in this database. Two or more are needed for comparative maps.");
    return;
  }
  else if (m == 1) {
    messout ("There is only one Map_set in this database. Two or more are needed for comparative maps.");
    return;
  }    

  if (!displayCreate (OxgridChoose))
    return ;

  graphRegister (PICK, choosePick) ;
  ox->chooseGraph = graphActive() ;
  graphAssociate (&OX_MAGIC, ox) ;

  getMapList () ;
  oxgridDefine () ;

}

int getMapsetList (void)
{
  int i, m ;
  KEY key ;
  KEYSET mapsets ;

/* Create a list of all ?Map_sets currently held in the database */

  mapsets = keySetCreate () ;
  mapsets = query (0, ">?Map_set") ;
  mapsetListArray = arrayCreate (keySetMax (mapsets) + 1 , FREEOPT) ;
  for (i = 0 ; i < keySetMax (mapsets) ; i++)
    {key = keySet (mapsets, i) ;
     arrayp (mapsetListArray, i+1, FREEOPT)->key = key ;
     arrayp (mapsetListArray, i+1, FREEOPT)->text = name (key) ;
     }
  m = keySetMax (mapsets) ;
  keySetDestroy (mapsets) ;
  arrayp (mapsetListArray, 0, FREEOPT)->key = i ;
  arrayp (mapsetListArray, 0, FREEOPT)->text = "Map_sets ..." ;

  mapsetList = arrp (mapsetListArray, 0, FREEOPT) ;
  return m ;
}

void getMapList (void)
{
  int i ;
  KEY key ;
  KEYSET maps ;

/* Create a list of all ?Maps currently held in the database */

  maps = keySetCreate () ;
  maps = query (0, ">?Map") ;
  mapListArray = arrayCreate (keySetMax (maps) + 1 , FREEOPT) ;
    for (i = 0 ; i < keySetMax (maps) ; i++)
    {key = keySet (maps, i) ;
     arrayp (mapListArray, i+1, FREEOPT)->key = key ;
     arrayp (mapListArray, i+1, FREEOPT)->text = name (key) ;
     }
  keySetDestroy (maps) ;
  arrayp (mapListArray, 0, FREEOPT)->key = i ;
  arrayp (mapListArray, 0, FREEOPT)->text = "Maps ..." ;
  
  mapList = arrp (mapListArray, 0, FREEOPT) ;


}


static void choosePick (int box)
{
  KEY key1, key2, key3, key4 ;
  OXGET ("choosePick") ;

  if (!box)
    return ;
  if (box == ox->sp1Box)
   {if (graphSelect (&key1, mapsetList))
      {ox->mapset1 = key1 ;
       ox->mapset1p = freekey2text (ox->mapset1, mapsetList) ;
       }
    graphBoxDraw (box, BLACK, LIGHTBLUE) ;
    
    }
  else if (box == ox->sp2Box)
    {if (graphSelect (&key2, mapsetList))
       {ox->mapset2 = key2 ;
        ox->mapset2p = freekey2text (ox->mapset2, mapsetList) ;
       }
     graphBoxDraw (box, BLACK, LIGHTBLUE) ;
    }
  else if (box == ox->searchBox)
    {if (graphSelect (&key3, searchChoice))
       {ox->search = key3 ;
        ox->searchp = freekey2text (ox->search, searchChoice) ;
       }
     graphBoxDraw (box, BLACK, PALEMAGENTA) ;
    }
  else if (box == ox->mapBox)
    {if (graphSelect (&key4, mapList))
       {ox->specmap = key4 ;
        ox->specmapp = freekey2text (ox->specmap, mapList) ;
       }
     graphBoxDraw (box, BLACK, LIGHTBLUE) ;
    }


}

void oxgridDefine (void)
{
  int box ;
  OXGET ("oxgridDefine") ;

/* Draw the pull-down menu buttons in the Comparative Map Chooser window */

  graphClear() ;
  graphMenu (chooseMenu) ;
  graphButtons (chooseMenu, 1, 1, 55) ;

  graphColor (BLUE) ;
  graphRectangle (3, 3, 55, 17.5) ;
  graphRectangle (3, 18.5, 55, 23) ;
  graphColor (BLACK) ;
 
  graphText ("Map_set 1 :", 7, 5) ;
  box = ox->sp1Box = graphBoxStart () ;
  graphTextPtrPtr (&ox->mapset1p, 20, 5, 30) ;
  ox->mapset1p = freekey2text (ox->mapset1, mapsetList) ;
  graphBoxEnd () ;
  graphBoxBox (box) ;
  if (mapsetList->key > 0) 
    graphBoxFreeMenu (box, (FreeMenuFunction) oxDefineSpecies1, mapsetList) ;
     
  graphText ("Map_set 2 :", 7, 10) ;
  box = ox->sp2Box = graphBoxStart () ;
  graphTextPtrPtr (&ox->mapset2p, 20, 10, 30) ;
  ox->mapset2p = freekey2text (ox->mapset2, mapsetList) ;
  graphBoxEnd () ;
  graphBoxBox (box) ;
  if (mapsetList->key > 0) 
    graphBoxFreeMenu (box, (FreeMenuFunction) oxDefineSpecies2, mapsetList) ;

  graphText ("Homology Search :", 7, 15) ;
  box = ox->searchBox = graphBoxStart () ;
  graphTextPtrPtr (&ox->searchp, 26, 15, 24) ;
  ox->searchp = freekey2text (ox->search, searchChoice) ;
  graphBoxEnd () ;
  graphBoxBox (box) ;
  graphBoxFreeMenu (box, (FreeMenuFunction) oxDefineSearch, searchChoice) ;

  graphText ("Species Grid Map :", 7, 20) ;
  box = ox->mapBox = graphBoxStart () ;
  graphTextPtrPtr (&ox->specmapp, 27, 20, 23) ;
  ox->specmapp = freekey2text (ox->specmap, mapList) ;
  graphBoxEnd () ;
  graphBoxBox (box) ;
  graphBoxFreeMenu (box, (FreeMenuFunction) oxDefineMap, mapList) ;

  graphRedraw() ;

}

static void oxDefineSpecies1 (KEY k, int box)
{
  OXGET ("oxDefineSpecies1") ;

  ox->mapset1 = k ;
  ox->mapset1p = freekey2text (k, mapsetList) ;
  graphBoxDraw (box, BLACK, LIGHTBLUE) ;

}

static void oxDefineSpecies2 (KEY k, int box)
{
  OXGET ("oxDefineSpecies2") ;

  ox->mapset2 = k ;
  ox->mapset2p = freekey2text (k, mapsetList) ;
  graphBoxDraw (box, BLACK, LIGHTBLUE) ;

}

static void oxDefineSearch (KEY k, int box)
{
  OXGET ("oxDefineSearch") ;

  ox->search = k ;
  ox->searchp = freekey2text (k, searchChoice) ;
  if (k == 's') ox->flag |= FLAG_SELF ;
  else if (k == 'n') ox->flag &= ~FLAG_SELF ;
  graphBoxDraw (box, BLACK, PALEMAGENTA) ;

}

static void oxDefineMap (KEY k, int box)
{
  OXGET ("oxDefineMap") ;

  ox->specmap = k ;
  ox->specmapp = freekey2text (k, mapList) ;
  graphBoxDraw (box, BLACK, LIGHTBLUE) ;

}



/***********************************************************************/
/*************************** oxgrid functions **************************/
/***********************************************************************/


void oxgridDisplay (void)
{ 
  int homol ;

  OXGET("oxgridDisplay") ;
  
  /* If homologies are found for a good definition, draw an Oxford Grid */
  
  homol = oxgridConvert (ox) ;
  
  if (homol == -1)
    {messout ("Sorry, definition is bad") ;
    return ;
    }
  else if (homol == 0)
    {messout ("No homologies for these Map_sets") ;
    return ;
    }
  else
    {if (graphExists (ox->graph))
      {graphActivate (ox->graph) ;
      graphPop () ;}
    else
      {if (!displayCreate (OXGRID))
       return ;
      graphAssociate (&OX_MAGIC, ox) ;
      ox->graph = graphActive () ;
      graphRegister (DESTROY, oxgridDestroy) ;
      graphRegister (PICK, (GraphFunc) oxgridPick) ;
      graphRegister (KEYBOARD, (GraphFunc) oxgridKeyboard) ;
      graphRegister (MESSAGE_DESTROY, displayUnBlock) ;
      graphMenu (oxgridMenu) ;
      
      oxgridDraw (ox) ;
      }
    }
}

BOOL chkSpecies (KEY key1, KEY key2)
{ 

/* function used everywhere to check whether two keys have the same name */

  if (strcmp (name (key1), name (key2)) == 0)
    return (TRUE) ;
  else return (FALSE) ;

}


static void oxgridDestroy (void)
{
  Graph g1, g2, g3, g4 ;
  OXGET("oxgridDestroy") ;

/* Destroy all graphs except the Comparative Map Chooser, all arrays and */
/* all keys held in the OX structure */

  arrayDestroy (ox->segs)    ;
  arrayDestroy (ox->boxes)   ;
  arrayDestroy (ox->fishes)  ;
  arrayDestroy (ox->axis1)   ;
  arrayDestroy (ox->axis2)   ;
  arrayDestroy (ox->extent1) ;
  arrayDestroy (ox->extent2) ;

  arrayDestroy (ox->points)  ;
  arrayDestroy (ox->c1segs)  ;
  arrayDestroy (ox->c2segs)  ;

  arrayDestroy (ox->many)    ;
  arrayDestroy (ox->csegs)   ;
  arrayDestroy (ox->c3segs)  ;
  arrayDestroy (ox->axis3)   ;
  arrayDestroy (ox->extent3) ;

  ox->xmax = 0 ;
  ox->ymax = 0 ;

  ox->c1 = 0 ;
  ox->c2 = 0 ;
  ox->c3 = 0 ;

  ox->s1 = 0 ;
  ox->s2 = 0 ;
  ox->s3 = 0 ;
  ox->s4 = 0 ;
  ox->s5 = 0 ;

  if (graphAssFind (&OX_MAGIC, &ox)
      && ox
      && (ox->magic == &OX_MAGIC))
     {g1 = ox->hlGraph ;
      g2 = ox->manyGraph ;
      g3 = ox->pairGraph ;
      g4 = ox->graph ;
      if (graphExists (g1))
	{ graphActivate (g1) ;
	  graphDestroy () ;
	}
      if (graphExists (g2))
	{ graphActivate (g2) ;
	  graphDestroy () ;
	}
      if (graphExists (g3))
	{ graphActivate (g3) ;
	  graphDestroy () ;
	}
      if (graphExists (g4))
	{ graphActivate (g4) ;
	  graphDestroy () ;
	}
   
      }
 
  ox->chooseGraph = graphActive () ;

}




/**********************************/


static void oxgridSelect (OX ox, int box)
{
  if (ox->gridBox == box) return ;

  if (ox->gridBox)
   {if (ox->gridBox >= MIN_LIVE_BOX + ox->xmax + ox->ymax + 2 && ox->gridBox <= MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->xmax*ox->ymax + 1) 
      graphBoxDraw (ox->gridBox, BLACK, WHITE) ;
 
    else if (ox->gridBox >= MIN_LIVE_BOX + 2 && ox->gridBox <= MIN_LIVE_BOX + ox->xmax + ox->ymax + 1)
      graphBoxDraw (ox->gridBox, BLUE, WHITE) ;

    else if (ox->gridBox == MIN_LIVE_BOX || ox->gridBox == MIN_LIVE_BOX + 1)
     {graphTextFormat (BOLD) ;
      graphBoxDraw (ox->gridBox, BLACK, WHITE) ;
      graphTextFormat(PLAIN_FORMAT);
     }
   }
     
  ox->gridBox = box ;
 
  if (ox->gridBox >= MIN_LIVE_BOX + ox->xmax + ox->ymax + 2 && ox->gridBox <= MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->xmax*ox->ymax + 1)
    graphBoxDraw (ox->gridBox, BLACK, GREEN) ;

  else if (ox->gridBox == MIN_LIVE_BOX || ox->gridBox == MIN_LIVE_BOX + 1)
   {graphTextFormat (BOLD) ;
    graphBoxDraw (ox->gridBox, WHITE, BLACK) ;
    graphTextFormat(PLAIN_FORMAT);
   }

  else if (ox->gridBox >= MIN_LIVE_BOX + 2 && ox->gridBox <= MIN_LIVE_BOX + ox->xmax + ox->ymax + 1)
    graphBoxDraw (ox->gridBox, WHITE, BLUE) ;

  graphCompScrollEntry (0, ox->findhom, 0, 0, 0, 0, 0) ; 
  graphEntryDisable () ;
  graphCompScrollEntry (0, ox->findloc, 0, 0, 0, 0, 0) ;
  graphEntryDisable () ;

}

/***********************************/

static void oxgridPick (int box) 
{

  OXGET ("oxgridPick") ;
    
  if (box >= MIN_LIVE_BOX + ox->xmax + ox->ymax + 2 && box <= MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->xmax*ox->ymax + 1)    
    {if (box == ox->gridBox)
       oxhomlist () ;
     graphActivate (ox->graph) ; 
     oxgridSelect (ox, box) ;
    }

  else if (box >= MIN_LIVE_BOX && box <= MIN_LIVE_BOX + ox->xmax + ox->ymax + 1)
   {if (box == ox->gridBox)
      display (arrayp (ox->boxes, box - MIN_LIVE_BOX, BOX)->key, 0, 0) ;
    graphActivate (ox->graph) ;
    oxgridSelect (ox, box) ;
   }

  else if (box == FINDHOM_BOX)
    graphCompScrollEntry (0, ox->findhom, 0, 0, 0, 0, 0) ;

  else if (box == FINDLOC_BOX)
    graphCompScrollEntry (0, ox->findloc, 0, 0, 0, 0, 0) ;
        
}


/*************************************************/

static void oxgridKeyboard (int k)
{
  int box ;
  int x, y ;
  BOX *bx ;
  OXGET("oxgridKeyboard") ;

  if (k == 'f' || k == 'F')
    { oxgridFlip () ;
      return ;
    }

  if ((!ox->gridBox) || (ox->h == 0))
    return ;

/*  messout ("box = %d", ox->gridBox) ;*/

  if (k == RETURN_KEY && (ox->gridBox >= MIN_LIVE_BOX + ox->xmax + ox->ymax + 2 && ox->gridBox <= MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->xmax*ox->ymax + 1))
    {oxhomlist () ;
     return ;
    }
  else if (k == RETURN_KEY && (ox->gridBox >= MIN_LIVE_BOX && ox->gridBox <= MIN_LIVE_BOX + ox->xmax + ox->ymax + 1))
   {display (arrayp (ox->boxes, ox->gridBox-MIN_LIVE_BOX, BOX)->key, 0, 0) ;
    return ;
    }

  if (ox->gridBox >= MIN_LIVE_BOX + ox->xmax + ox->ymax + 2)
   {box = ox->gridBox - MIN_LIVE_BOX ; 
    bx = arrayp (ox->boxes, box, BOX) ;
    x = bx->x ;  y = bx->y ;
   
    switch (k)
   {
    case LEFT_KEY:
      if (x >= 1)
	{x -= 1 ; box = MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->ymax*(x - 1) + y + 1 ;}
      else if (x == 1)
	{box = MIN_LIVE_BOX + ox->xmax + ox->ymax + y + 1 ;}
      break ;
    case RIGHT_KEY:
      if (x < ox->xmax - 1)
	{x += 1 ; box = MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->ymax*(x - 1) + y + 1 ;}
      else if (x == ox->xmax - 1)
	{box = MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->ymax*(ox->xmax - 1) + y + 1 ;}
      break ;
    case UP_KEY:
      if (y > 1)
	{y -= 1 ; box = MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->ymax*(x - 1) + y + 1 ;}
      else if (y == 1)
	{box = MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->ymax*(x - 1) ; }
      break ;
    case DOWN_KEY:
      if (y < ox->ymax - 1)
	{y += 1 ; box = MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->ymax*(x - 1) + y + 1 ;}
      else if (y == ox->ymax - 1)
	{box = MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->ymax*(x - 1) + ox->ymax + 1 ;}
      break ;
    }
    oxgridSelect (ox, box) ;
    }

}

/*************************************************/
/************* conversion and drawing ************/

static void oxgridFlip (void)
{
  Array axis, extent, tmp ;
  KEY gene, map, spec ;
  SEG *seg ;
  int i, j, itmp ;
  float z ;

  OXGET ("oxgridFlip") ;

  ox->flag ^= FLAG_OX_FLIP ;

  spec = ox->s1 ;
  ox->s1 = ox->s2 ;
  ox->s2 = spec ;

  itmp = ox->xmax ;
  ox->xmax = ox->ymax ;
  ox->ymax = itmp ;

/*
  tmp = ox->axis1 ;
  ox->axis1 = ox->axis2 ;
  ox->axis2 = tmp ;

  tmp = ox->extent1 ;
  ox->extent1 = ox->extent2 ;
  ox->extent2 = tmp ;
*/

  axis = arrayCreate (64, KEY) ;
  extent = arrayCreate (64, float) ;

  axis = arrayCopy (ox->axis1) ;
  ox->axis1 = arrayCopy (ox->axis2) ;
  ox->axis2 = arrayCopy (axis) ;
  arrayDestroy (axis) ;

  extent = arrayCopy (ox->extent1) ;
  ox->extent1 = arrayCopy (ox->extent2) ;
  ox->extent2 = arrayCopy (extent) ;
  arrayDestroy (extent) ;


  for (i = 0 ; i < arrayMax (ox->segs) ; i++)
    {seg = arrayp (ox->segs, i, SEG) ;
     j = seg->x ;
     seg->x = seg->y ;
     seg->y = j ;
     seg->box = seg->x*ox->ymax + seg->y + (ox->xmax + ox->ymax + 2) + MIN_LIVE_BOX ;
     gene = seg->gene1 ;
     seg->gene1 = seg->gene2 ;
     seg->gene2 = gene ;
     map = seg->map1 ;
     seg->map1 = seg->map2 ;
     seg->map2 = map ;
     z = seg->z1 ;
     seg->z1 = seg->z2 ;
     seg->z2 = z ;
    }

  for (i = 0 ; i < arrayMax (ox->fishes) ; i++)
    {seg = arrayp (ox->fishes, i, SEG) ;
     j = seg->x ;
     seg->x = seg->y ;
     seg->y = j ;
     seg->box = seg->x*ox->ymax + seg->y + (ox->xmax + ox->ymax + 2) + MIN_LIVE_BOX ;
     map = seg->map1 ;
     seg->map1 = seg->map2 ;
     seg->map2 = map ;
     z = seg->z1 ;
     seg->z1 = seg->z2 ;
     seg->z2 = z ;
     z = seg->z3 ;
     seg->z3 = seg->z4 ;
     seg->z4 = z ;
    }
  
  oxgridDraw (ox) ;  
}


/*******************************CONVERT***********************************/

static int oxgridConvert (OX ox)
{ 
  int i, j, k, l, s ;
  int f, h, h1, h2, h3 ;
  int m1, m2 ;
  int phe1, phe2 ;
  int int1, null1, null2 ;
  OBJ sp1, sp2, Group ; 
  OBJ Locus1, Locus2, Map ;
  KEY group, key, key1, key2 ;
  KEY map1, map2 ;
  SEG *seg, *seg2 ;
  Array maps, maps1, maps2 ;
  Array homol1, homol2, homol3 ;
  KEYSET pairwise = 0, fuzzy = 0, sfuzzy = 0 ;
  BSMARK mark = 0, mark1 = 0, mark2 = 0 ;
  BOOL found ;
  float z1 = 0, z2 = 0 ;

  ox->flag |= FLAG_OX_FLIP ;
  ox->flag |= FLAG_PAIR_FLIP ;
  ox->flag &= FLAG_EQUAL ;

  ox->wide = 75 ;
  ox->high = 50 ;

  if (!(sp1 = bsCreate (ox->mapset1)) || !(sp2 = bsCreate (ox->mapset2)))
    return -1 ;
  else
    {ox->s1 = ox->mapset1 ;
     ox->s2 = ox->mapset2 ;
    }
  bsDestroy (sp1) ;
  bsDestroy (sp2) ;

  ox->xmax = 0 ;
  ox->ymax = 0 ;
  ox->segs = arrayCreate (128, SEG) ;
  ox->fishes = arrayCreate (128, SEG) ;
  ox->extent1 = arrayCreate (64, float) ;
  ox->axis1 = arrayCreate (64, KEY) ;
  ox->extent2 = arrayCreate (64, float) ;
  ox->axis2 = arrayCreate (64, KEY) ;
  maps = arrayCreate (32, BSunit) ;
  maps1 = arrayCreate (32, BSunit) ;
  maps2 = arrayCreate (32, BSunit) ;
  homol1 = arrayCreate (32, BSunit) ;
  homol2 = arrayCreate (32, BSunit) ;
  homol3 = arrayCreate (32, BSunit) ;
  
  null1  = 0 ;
   if ((sp1 = bsCreate (ox->s1)))
     {if (bsFindTag (sp1, _Map) && bsFlatten (sp1, 1, maps))    
        {
	 for (i = 0 ; i < arrayMax (maps) ; i++)
	   {
	     key = array (maps, i, BSunit).k ;
	     array (ox->axis1, i, KEY) = key ;     
	     array (ox->extent1, i*2, float) = 0 ;
	     array (ox->extent1, i*2 + 1, float) = 0 ;
	     if (bIndexFind (key, _Extent) &&
		 (Map = bsCreate (key)))
	       {
		 if (bsGetData (Map,  _Extent, _Float, &z1) &&
		    bsGetData (Map,  _bsRight, _Float, &z2))
		  {
		    array (ox->extent1, i*2, float) = z1 ;
		   array (ox->extent1, i*2 + 1, float) = z2 ;
                   if (!(abs(z1-z2) > 0)) null1++ ;
		  }
		bsDestroy (Map) ;
	       }
            }
         }
      }

  ox->xmax = arrayMax(maps) ;
  arrayDestroy (maps) ;
  maps = arrayCreate (32, BSunit) ;

  null2 = 0 ;
   if ((sp2 = bsCreate (ox->s2))  )
     {if (bsFindTag (sp2, _Map) && bsFlatten (sp2, 1, maps)) 
        {
	 for (i = 0 ; i < arrayMax (maps) ; i++)
	   {
	     key = array (maps, i, BSunit).k ;
	     array (ox->axis2, i, KEY) = key ;
	     array (ox->extent2, i*2, float) = 0 ;
	     array (ox->extent2, i*2 + 1, float) = 0 ;
	     if (bIndexFind (key, _Extent) && (Map = bsCreate (key)))
	       {
		 if (bsGetData (Map,  _Extent, _Float, &z1) &&
		     bsGetData (Map,  _bsRight, _Float, &z2))
		   {
		     array (ox->extent2, i*2, float) = z1 ;
		     array (ox->extent2, i*2 + 1, float) = z2 ;
		     if (!(abs(z1-z2) > 0)) null2++ ;
 		   }
		 bsDestroy (Map) ;
	       }
	    }
         }
      }

  ox->ymax = arrayMax(maps) ;
  arrayDestroy (maps) ; 
  bsDestroy (sp1) ;
  bsDestroy (sp2) ;

  printf ("null1 = %d and null2 = %d\n", null1, null2) ;
  if (null1 > 0 || null2 > 0) 
    {ox->flag |= FLAG_NOSIZE ;
     ox->flag |= FLAG_EQUAL ;
    }
  if (chkSpecies (ox->s1, ox->s2)) ox->flag |= FLAG_PARALOGY ;
  else ox->flag &= ~FLAG_PARALOGY ;
  printf ("ox->flag = %d\n", ox->flag) ;

  f = 0 ;

  for (m1 = 0 ; m1 < ox->xmax ; m1++)
    for (m2 = 0 ; m2 < ox->ymax ; m2++)
      {
	KEYSET fishSet = 0 ;
	int k, ok ;
	KEY map1, map2, fish ;
	float z1, z2, z3, z4 ;
	OBJ Fish = 0 ;
	
	map1 =  keySet (ox->axis1, m1) ;
	map2 =  keySet (ox->axis2, m2) ;
	fishSet =  query (0, messprintf("{Find Map %s ; > Fish} SETAND {Find Map %s ; > Fish}", name(map1), name(map2))) ;

	 
	if (keySetMax(fishSet)) printf("fishing for %s %s got %d fishes \n", name(map1), name(map2), keySetMax(fishSet)) ;
	  
	for (k = 0 ; k < keySetMax(fishSet) ; k++)
	  {
	    fish = keySet(fishSet, k) ;
	    ok = 0 ;
	    if ((Fish = bsCreate (fish)))
	      {
		if (bsFindKey (Fish, _Map, map1) &&
		    bsPushObj (Fish) &&
		    bsGetData (Fish, _Left, _Float, &z1) &&
		    bsGetData (Fish, _Right, _Float, &z3))
		  ok++ ;
		bsGoto (Fish, 0) ;
		if (bsFindKey (Fish, _Map, map2) &&
		    bsPushObj (Fish) &&
		    bsGetData (Fish, _Left, _Float, &z2) &&
		    bsGetData (Fish, _Right, _Float, &z4))
		  ok++ ;
		bsDestroy (Fish) ;
	      }
	    if (ok == 2)
	      {
		seg = arrayp (ox->fishes, f++, SEG) ;
		seg->flag |= FLAG_FISH ;
		seg->box = m1*ox->ymax+m2 + (ox->xmax+ox->ymax+2) + MIN_LIVE_BOX;
		seg->map1 = map1 ;
		seg->map2 = map2 ; 
		seg->x = m1 ;
		seg->y = m2 ;
		seg->z1 = z1 ;
		seg->z2 = z2 ;
		seg->z3 = z3 ;
		seg->z4 = z4 ;
	      }
	  }
	keySetDestroy (fishSet) ;

      }

  printf ("%d fishes found\n", f) ;

  h = 0 ;

  pairwise = keySetCreate () ;
  pairwise = query (0, messprintf ("{>?Map_set \"%s\" ;>Map ; >Locus ; >Pairwise} SETAND {>?Map_set \"%s\" ;>Map ; >Locus ; >Pairwise}", name (ox->mapset1), name (ox->mapset2))) ;
  for (i = 0 ; i < keySetMax (pairwise) ; i++)
    {
      group = array (pairwise, i, KEY) ;
      Group = bsCreate (group) ;
      if (bsGetKey (Group, _Pairwise, &key1)) do
	{
	  mark = bsMark (Group, mark) ;
	  bsFlatten (Group, 1, homol1) ;
	  Locus1 = bsCreate (key1) ;

	  if (bsFindTag (Locus1, _Eye)) phe1 = 1 ;
	  else if (bsFindTag (Locus1, _Retina)) phe1 = 2 ;
	  else if (bsFindTag (Locus1, _RP)) phe1 = 3 ;
	  else if (bsFindTag (Locus1, _Lens)) phe1 = 4 ;
	  else if (bsFindTag (Locus1, _Cataract)) phe1 = 5 ;
	  else if (bsFindTag (Locus1, _Coloboma)) phe1 = 6 ;
	  else if (bsFindTag (Locus1, _Global_malformation)) phe1 = 7 ;   
	  else if (bsFindTag (Locus1, _Localised_to_eye)) phe1 = 8 ;
	  else if (bsFindTag (Locus1, _Generalised)) phe1 = 9 ;
	  else if (bsFindTag (Locus1, _Anophthalmos_unilateral)) phe1 = 10 ;
	  else if (bsFindTag (Locus1, _Anophthalmos)) phe1 = 11 ;
	  else phe1 = 0 ;

	  if (bsGetKey (Locus1, _Map, &map1))
	    do
	      {	j = keySetMax (ox->axis1) ;
		while (j--)
		  if (map1 == keySet(ox->axis1, j))
		    break ;
		if (j < 0)
		  continue ;

		z1 = -9999 ;
		mark1 = bsMark (Locus1, mark1) ;
		if ( bsPushObj (Locus1))
		  bsGetData(Locus1, _Position, _Float, &z1) ;
		bsGoto (Locus1, mark1) ;
		
		
		for (k = 0 ; k < arrayMax (homol1) ; k++)
		  {
		    key2 = array (homol1, k, BSunit).k ;
		    Locus2 = bsCreate (key2) ;

		    if (bsFindTag (Locus2, _Eye)) phe2 = 1 ;
		    else if (bsFindTag (Locus2, _Retina)) phe2 = 2 ;
		    else if (bsFindTag (Locus2, _RP)) phe2 = 3 ;
		    else if (bsFindTag (Locus2, _Lens)) phe2 = 4 ;
		    else if (bsFindTag (Locus2, _Cataract)) phe2 = 5 ;
		    else if (bsFindTag (Locus2, _Coloboma)) phe2 = 6 ;
		    else if (bsFindTag (Locus2, _Global_malformation)) phe2 = 7 ;   
		    else if (bsFindTag (Locus2, _Localised_to_eye)) phe2 = 8 ;
		    else if (bsFindTag (Locus2, _Generalised)) phe2 = 9 ;
		    else if (bsFindTag (Locus2, _Anophthalmos_unilateral)) phe2 = 10 ;
		    else if (bsFindTag (Locus2, _Anophthalmos)) phe2 = 11 ;
		    else phe2 = 0 ;

		    if (bsGetKey (Locus2, _Map, &map2))
		      do
			{
			  l = keySetMax (ox->axis2) ;
			  while (l--)
			    if (map2 == keySet(ox->axis2, l))
			      break ;
			  if (l < 0)
			    continue ; 

			  z2 = -9999 ;	
			  mark2 = bsMark (Locus2, mark2) ;
			  if ( bsPushObj (Locus2))
			    bsGetData(Locus2, _Position, _Float, &z2) ;
			  bsGoto (Locus2, mark2) ;
			  
			  for (m1 = 0 ; m1 < ox->xmax ; m1++)
			    for (m2 = 0 ; m2 < ox->ymax ; m2++)
			      {
				if (chkSpecies (map1, array (ox->axis1, m1, KEY)) && chkSpecies (map2, array (ox->axis2, m2, KEY)))
				  {
				    found = FALSE ;
				    for (s = 0 ; s < arrayMax (ox->segs) ; s++)
				      {
					seg2 = arrayp (ox->segs, s, SEG) ;
					if (chkSpecies (key1, seg2->gene1) && chkSpecies (key2, seg2->gene2)) found = TRUE ;
				      }
				    if (found == FALSE) 
				      {
					seg = arrayp (ox->segs, h , SEG) ;
					seg->key = group ;
					seg->flag = 0 ;
					seg->col = 0 ;
                                        seg->phe1 = phe1 ;
                                        seg->phe2 = phe2 ;
					seg->box = m1*ox->ymax+m2 + (ox->xmax+ox->ymax+2) + MIN_LIVE_BOX;
					seg->x = m1 ;
					seg->y = m2 ;
					seg->gene1 = key1 ;
					seg->gene2 = key2 ;
					seg->map1 = map1 ;
					seg->map2 = map2 ;
					seg->z1 = z1 ;
					seg->z2 = z2 ;
					h++ ;
                                        printf ("homology %d: %s between map %s and map %s\n\n", h, name (group), name(map1), name (map2)) ;
				      }
				  }
				else if (chkSpecies (map2, array (ox->axis1, m1, KEY)) && chkSpecies (map1, array (ox->axis2, m2, KEY)))
				  {
				    found = FALSE ;
				    for (s = 0 ; s < arrayMax (ox->segs) ; s++)
				      {
					seg2 = arrayp (ox->segs, s, SEG) ;
					if (chkSpecies (key1, seg2->gene1) && chkSpecies (key2, seg2->gene2)) found = TRUE ;
				      }
				    if (found == FALSE) 
				      {
					seg = arrayp (ox->segs, h , SEG) ;
					seg->key = group ;
					seg->flag = 0 ;
					seg->col = 0 ;
					seg->phe1 = phe2 ;
                                        seg->phe2 = phe1 ;
					seg->box = m1*ox->ymax+m2 + (ox->xmax+ox->ymax+2) + MIN_LIVE_BOX;
					seg->x = m1 ;
					seg->y = m2 ;
					seg->gene1 = key2 ;
					seg->gene2 = key1 ;
					seg->map1 = map2 ;
					seg->map2 = map1 ;  
					seg->z1 = z1 ;
					seg->z2 = z2 ;
					h++ ;
					printf ("homology %d: %s between map %s and map %s\n\n", h, name (group), name(map1), name (map2)) ;
				      }
				  }     
			      }
			} while (bsGetKey (Locus2, _bsDown, &map2)) ;
		    bsDestroy (Locus2) ;
		  } 
	      } while (bsGetKey (Locus1, _bsDown, &map1)) ;
	  bsDestroy (Locus1) ;
	} while (bsGoto (Group, mark), bsGetKey (Group, _bsDown, &key1)) ;
      bsDestroy (Group) ;
    }
  bsMarkFree (mark) ;
  bsMarkFree (mark1) ;
  bsMarkFree (mark2) ;

/*  
  if (!(ox->flag & FLAG_SELF))
   {fuzzy = keySetCreate () ;
    fuzzy = query (0, messprintf (">?Map_set \"%s\" ; >Map ; >Locus ; >Fuzzy", name (ox->mapset1))) ;
    for (i = 0 ; i < keySetMax (fuzzy) ; i++)
     {group = array (fuzzy, i, KEY) ;
      Group = bsCreate (group) ;
      if (bsFindTag (Group, _Fuzzy)) 
       {bsFlatten (Group, 1, homol2) ;
        for (h1 = 0 ; h1 < arrayMax (homol2) - 1 ; h1++)
          for (h2 = h1 + 1 ; h2 < arrayMax (homol2) ; h2++)
            {key1 = array (homol2, h1, BSunit).k ;
             key2 = array (homol2, h2, BSunit).k ;
             Locus1 = bsCreate (key1) ;
             Locus2 = bsCreate (key2) ;
             if (bsFindTag (Locus1, _Map) && bsFindTag (Locus2, _Map)) 
              {bsGetKey (Locus1, _Map, &map1) ;
               bsGetKey (Locus2, _Map, &map2) ;
               for (m1 = 0 ; m1 < arrayMax (ox->axis1) ; m1++)
                 for (m2 = 0 ; m2 < arrayMax (ox->axis2) ; m2++)
                   {if (chkSpecies (map1, array (ox->axis1, m1, KEY)) && chkSpecies (map2, array (ox->axis2, m2, KEY)) && !(chkSpecies (key1, key2)))
                     {seg = arrayp (ox->segs, h , SEG) ;
                      seg->key = group ;
                      seg->flag = 0 ;
                      seg->box = (m1-1)*ox->ymax+m2-1 + (ox->xmax+ox->ymax+2) + MIN_LIVE_BOX ;
                      seg->x = m1 ;
                      seg->y = m2 ;
                      seg->gene1 = key1 ;
                      seg->gene2 = key2 ;
                      seg->map1 = map1 ;
                      seg->map2 = map2 ;
                      h++ ;
                     }
                    else if (chkSpecies (map2, array (ox->axis1, m1, KEY)) && chkSpecies (map1, array (ox->axis2, m2, KEY)) && !(chkSpecies (key1, key2)))
                     {seg = arrayp (ox->segs, h , SEG) ;
                      seg->key = group ;
                      seg->flag = 0 ;
                      seg->box = (m1-1)*ox->ymax+m2-1 + (ox->xmax+ox->ymax+2) + MIN_LIVE_BOX ;
                      seg->x = m1 ;
                      seg->y = m2 ;
                      seg->gene1 = key2 ;
                      seg->gene2 = key1 ;
                      seg->map1 = map2 ;
                      seg->map2 = map1 ;
                      h++ ;
                      }   
                    } 
                 } 
              bsDestroy (Locus1) ;
              bsDestroy (Locus2) ;
             }
         }
     bsDestroy (Group) ;
     }
   }
*/

/*
    fuzzy = keySetCreate () ;
    fuzzy = query (0, messprintf (">?Map_set \"%s\" ; >Map ; >Locus ; >Fuzzy", name (ox->mapset1))) ;
    for (i = 0 ; i < keySetMax (fuzzy) ; i++)
      {group = array (fuzzy, i, KEY) ;
       Group = bsCreate (group) ;
       if (bsFindTag (Group, _Fuzzy)) 
         {bsFlatten (Group, 1, homol2) ;
          for (h1 = 0 ; h1 < arrayMax (homol2) - 1 ; h1++)
             for (h2 = h1 + 1 ; h2 < arrayMax (homol2) ; h2++)
               {key1 = array (homol2, h1, BSunit).k ;
                key2 = array (homol2, h2, BSunit).k ;
                Locus1 = bsCreate (key1) ;
                Locus2 = bsCreate (key2) ;
                if (bsFindTag (Locus1, _Map) && bsFindTag (Locus2, _Map)) 
                  {bsFlatten (Locus1, 1, maps1) ;
                   bsFlatten (Locus2, 1, maps2) ;
                   for (j = 0 ; j < arrayMax (maps1) ; j++)
                      for (k = 0 ; k < arrayMax (maps2) ; k++)
                          {map1 = array (maps1, j, BSunit).k ;
                           map2 = array (maps2, k, BSunit).k ;
                           for (m1 = 0 ; m1 < arrayMax (ox->axis1) ; m1++)
                               for (m2 = 0 ; m2 < arrayMax(ox->axis2) ; m2++)
                                  {if (chkSpecies (map1, array (ox->axis1, m1, KEY)) && chkSpecies (map2, array (ox->axis2, m2, KEY)) && !(chkSpecies (key1, key2)))
                                      {found = FALSE ;
                                       for (s = 0 ; s < arrayMax (ox->segs) ; s++)
                                          {seg2 = arrayp (ox->segs, s, SEG) ;
                                           if (chkSpecies (key1, seg2->gene1) && chkSpecies (key2, seg2->gene2)) found = TRUE ;
                                          }
                                       if (found == FALSE) 
                                          {seg = arrayp (ox->segs, h , SEG) ;
                                           seg->key = group ;
                                           seg->flag = 0 ;
					   seg->col = 0 ;
                                           seg->box = m1*ox->ymax+m2 + (ox->xmax+ox->ymax+2) + MIN_LIVE_BOX;
                                           seg->x = m1 ;
                                           seg->y = m2 ;
                                           seg->gene1 = key1 ;
                                           seg->gene2 = key2 ;
                                           seg->map1 = map1 ;
                                           seg->map2 = map2 ;
                                           h++ ;
                                          }
                                      }
                                   else if (chkSpecies (map2, array (ox->axis1, m1, KEY)) && chkSpecies (map1, array (ox->axis2, m2, KEY)) && !(chkSpecies (key1, key2)))
                                      {found = FALSE ;
                                       for (s = 0 ; s < arrayMax (ox->segs) ; s++)
                                          {seg2 = arrayp (ox->segs, s, SEG) ;
                                           if (chkSpecies (key2, seg2->gene1) && chkSpecies (key1, seg2->gene2)) found = TRUE ;
                                          }
                                       if (found == FALSE) 
                                          {seg = arrayp (ox->segs, h , SEG) ;
                                           seg->key = group ;
                                           seg->flag = 0 ;
					   seg->col = 0 ;
                                           seg->box = m1*ox->ymax+m2 + (ox->xmax+ox->ymax+2) + MIN_LIVE_BOX;
                                           seg->x = m1 ;
                                           seg->y = m2 ;
                                           seg->gene1 = key2 ;
                                           seg->gene2 = key1 ;
                                           seg->map1 = map2 ;
                                           seg->map2 = map1 ;
                                           h++ ;
                                          }
                                      }   
                                  } 
                          }
                  } 
                bsDestroy (Locus1) ;
                bsDestroy (Locus2) ;
              }
         }
       bsDestroy (Group) ;
    }
*/

/*
  if (ox->flag & FLAG_SELF)
   {arrayDestroy (maps2) ;
    maps2 = arrayCreate (32, BSunit) ;
    sfuzzy = keySetCreate () ;
    sfuzzy = keySetCopy (fuzzy) ;
    for (i = 0 ; i < keySetMax (sfuzzy) ; i++)
      {group = array (sfuzzy, i, KEY) ;
       Group = bsCreate (group) ;
       if (bsFindTag (Group, _Fuzzy)) 
         {bsFlatten (Group, 1, homol3) ;
          for (h3 = 0 ; h3 < arrayMax (homol3) ; h3++)
            {key = array (homol3, h3, BSunit).k ;
             Locus2 = bsCreate (key) ;
             if (bsFindTag (Locus2, _Map)) 
              {bsFlatten (Locus2, 1, maps2) ;
               for (h1 = 0 ; h1 < arrayMax (maps2) - 1 ; h1++)
                for (h2 = h1+1 ; h2 < arrayMax (maps2) ; h2++) 
                 {map1 = array(maps2, h1, BSunit).k ;
                  map2 = array(maps2, h2, BSunit).k ;
                  for (m1 = 0 ; m1 < arrayMax (ox->axis1) ; m1++)
                    for (m2 = 0 ; m2 < arrayMax (ox->axis2) ; m2++)
                      {if (chkSpecies (map1, array (ox->axis1, m1, KEY)) && chkSpecies (map2, array (ox->axis2, m2, KEY)))
                       {found = FALSE ;
                        for (s = 0 ; s < arrayMax (ox->segs) ; s++)
                         {seg2 = arrayp (ox->segs, s, SEG) ;
                          if (chkSpecies (key, seg2->gene1) && chkSpecies (key, seg2->gene2)) found = TRUE ;
                         }
                        if (found == FALSE) 
                         {seg = arrayp (ox->segs, h , SEG) ;
                          seg->key = group ;
                          seg->flag = 0 ;
			  seg->col = 0 ;
                          seg->box =(m1-1)*ox->ymax+m2-1+(ox->xmax+ox->ymax+2)+MIN_LIVE_BOX;
                          seg->x = m1 ;
                          seg->y = m2 ;
                          seg->gene1 = key ;
                          seg->gene2 = key ;
                          seg->map1 = map1 ;
                          seg->map2 = map2 ;
                          h++ ;
                         }
                        if (found == FALSE && (ox->flag & FLAG_PARALOGY))
                         {seg = arrayp (ox->segs, h , SEG) ;
                          seg->key = group ;
                          seg->flag = 0 ;
			  seg ->col = 0 ;
                          seg->box =(m2-1)*ox->ymax+m1-1+(ox->xmax+ox->ymax+2)+MIN_LIVE_BOX;
                          seg->x = m2 ;
                          seg->y = m1 ;
                          seg->gene1 = key ;
                          seg->gene2 = key ;
                          seg->map1 = map2 ;
                          seg->map2 = map1 ;
                          h++ ;
                         }                       
		       }
                       else if (chkSpecies (map2, array (ox->axis1, m1, KEY)) && chkSpecies (map1, array (ox->axis2, m2, KEY)))
                        {found = FALSE ;
                         for (s = 0 ; s < arrayMax (ox->segs) ; s++)
                          {seg2 = arrayp (ox->segs, s, SEG) ;
                           if (chkSpecies (key, seg2->gene1) && chkSpecies (key, seg2->gene2)) found = TRUE ;
                          }
                         if (found == FALSE) 
                          {seg = arrayp (ox->segs, h , SEG) ;
                           seg->key = group ;
                           seg->flag = 0 ;
			   seg->col = 0 ;
                           seg->box = (m1-1)*ox->ymax+m2-1+(ox->xmax+ox->ymax+2)+MIN_LIVE_BOX;
                           seg->x = m1 ;
                           seg->y = m2 ;
                           seg->gene1 = key ;
                           seg->gene2 = key ;
                           seg->map1 = map2 ;
                           seg->map2 = map1 ;
                           h++ ;
                          }
                         if (found == FALSE && (ox->flag & FLAG_PARALOGY))
                          {seg = arrayp (ox->segs, h , SEG) ;
                           seg->key = group ;
                           seg->flag = 0 ;
			   seg->col = 0 ;
                           seg->box = (m2-1)*ox->ymax+m1-1+(ox->xmax+ox->ymax+2)+MIN_LIVE_BOX;
                           seg->x = m2 ;
                           seg->y = m1 ;
                           seg->gene1 = key ;
                           seg->gene2 = key ;
                           seg->map1 = map1 ;
                           seg->map2 = map2 ;
                           h++ ;
                          }
                        }
                      }
                    }
                  }
                bsDestroy (Locus2) ;
               }
            }
          bsDestroy (Group) ;
         }   
      }
*/

  printf ("%d homologies found\n", h) ;

  arrayDestroy (pairwise) ;
  arrayDestroy (fuzzy) ;
  arrayDestroy (sfuzzy) ;
  arrayDestroy (maps1) ;
  arrayDestroy (maps2) ;
  arrayDestroy (homol1) ;
  arrayDestroy (homol2) ;
  arrayDestroy (homol3) ;
  ox->h = h ;
  return h ;
}

/***************************************/


static void oxgridDraw (OX ox)
{
  int   b, d, i, j, k, l, x ;
  int   box, count, dots, ldots ;
  float width1[50], width2[50], height1[50], height2[50] ;
  float diff, xlen, ylen, wtot=0, htot=0 ;
  float xscale=0, yscale=0 ;
  float t1, t2, t3, t4, spx, spy ;
  float fx1, fx2, fy1, fy2 ;
  SEG   *seg ;
  BOX   *bx ;
  float yOff = 0.5 ;
  float xOff = 2 ;

#ifndef RAND_MAX
#define RAND_MAX (32767)
#endif

  /* set up the menus at the top of the grid */
 
  graphClear () ;
  
  graphText (messprintf ("%s-%s Oxford Grid :", name (ox->s1), name (ox->s2)), 1, 1) ; 
  x = strlen (name (ox->s1)) + strlen (name (ox->s2)) + 17 ;
  xOff += 2 + strlen (name (ox->s2)) ;

  if (!(ox->flag & FLAG_EQUAL))
   {if (graphButton (" Equal  sizes ",toggleCells, 1, 3) != EQUAL_BOX)
    messcrash ("box screwup at %d in oxgridDraw", EQUAL_BOX) ;}
  else
   {if (graphButton ("Relative sizes",toggleCells, 1, 3) != EQUAL_BOX)
    messcrash ("box screwup at %d in oxgridDraw", EQUAL_BOX) ;}

  if (graphButton ("Flip", oxgridFlip, 17, 3) != FLIP_BOX)
    messcrash ("box screwup at %d in oxgridDraw", FLIP_BOX) ;  

  if (graphButton ("Highlight...", oxClear, 23, 3) != HIDE_BOX)
    messcrash ("box screwup at %d in oxgridDraw", HIDE_BOX) ;   

  graphBoxMenu (HIDE_BOX, LightMenu) ; 

  if (graphButton ("Statistics...", oxgridReorg, 37, 3) != REORG_BOX)
    messcrash ("box screwup at %d in oxgridDraw", REORG_BOX) ;

  graphBoxMenu (REORG_BOX, StatsMenu) ;

  graphText ("Find homology :", 1, 5) ;
  *ox->findhom = 0 ;
  if (graphCompScrollEntry (homolCompletion, ox->findhom, 20, 20, 17, 5, oxFindHom) != FINDHOM_BOX)
    messcrash ("box screwup at %d in oxgridDraw", FINDHOM_BOX) ;

  graphEntryDisable () ;

  graphText ("Find locus :", 41, 5) ;
  *ox->findloc = 0 ;
  if (graphCompScrollEntry (locusCompletion, ox->findloc, 20, 20, 54, 5, oxFindLoc) != FINDLOC_BOX)
    messcrash ("box screwup at %d in oxgridDraw", FINDLOC_BOX) ;

  graphEntryDisable () ;

  if (ox->hideFish)
    {if (graphButton (" Fish  ", oxgridFish, 52, 3) != FISH_BOX)
       messcrash ("box screwup at %d in oxgridDraw", FISH_BOX) ;}
  else
    {if (graphButton ("No fish", oxgridFish, 52, 3) != FISH_BOX)
       messcrash ("box screwup at %d in oxgridDraw", FISH_BOX) ;}

  if (!(ox->showRandom))
    {if (graphButton ("Random", oxgridRandom, 61, 3) != RANDOM_BOX)
       messcrash ("box screwup at %d in oxgridDraw", RANDOM_BOX) ;}
  else
    {if (graphButton ("Scaled", oxgridRandom, 61, 3) != RANDOM_BOX)
       messcrash ("box screwup at %d in oxgridDraw", RANDOM_BOX) ;}

  if (graphButton ("Eye phenotype...", oxClear2, 69, 3) != PHENOTYPE_BOX)
    messcrash ("box screwup at %d in oxgridDraw", PHENOTYPE_BOX) ;   

  graphBoxMenu (PHENOTYPE_BOX, PhenotypeMenu) ; 

  if (graphButton ("Zoom in", oxgridZoomIn, 87, 3) != ZOOM_IN)
    messcrash ("box screwup at %d in oxgridDraw", ZOOM_IN) ;

  if (graphButton ("Zoom out", oxgridZoomOut, 96, 3) != ZOOM_OUT)
    messcrash ("box screwup at %d in oxgridDraw", ZOOM_OUT) ;


  graphLine (0, 7, 200, 7) ; 
  yOff += 7.5 ;

  /* calculate cell boundaries and scaling factors for drawing the grid */

  for (i = 0 ; i < ox->xmax ; i++)
    {width1[i] = width2[i] = 0 ;}
  for (j = 0 ; j < ox->ymax ; j++)
    {height1[j] = height2[j] = 0 ;}
  wtot = 0 ; htot = 0 ;

  if (!(ox->flag & FLAG_EQUAL))
    {
      width1[0] = 0 ;
      t1 = array(ox->extent1, 0, float) ;
      t2 = array(ox->extent1, 1, float) ;
      width2[0] = t2 > t1 ? t2 - t1 : t1 - t2 ;

      for (j = 1 ; j < ox->xmax ; j++) 
	{ width1[j] = width2[j-1] ;
	  t1 = array(ox->extent1, 2*j, float) ;
	  t2 = array(ox->extent1, 2*j+1, float) ;
	  diff = t2 > t1 ? t2 - t1 : t1 - t2 ;
          width2[j] = width2[j-1] + diff ;
	}

      height1[0] = 0 ;
      t1 = array(ox->extent2, 0, float) ;
      t2 = array(ox->extent2, 1, float) ;
      height2[0] = t2 > t1 ? t2 - t1 : t1 - t2 ;

      for (j = 1 ; j < ox->ymax ; j++) 
	{ height1[j] = height2[j-1] ;
	  t1 = array(ox->extent2, 2*j, float) ;
	  t2 = array(ox->extent2, 2*j+1 , float) ;
	  diff = t2 > t1 ? t2 - t1 : t1 - t2 ;
          height2[j] = height2[j-1] + diff ;
	}
    }
  else
    { 
      for (i = 0 ; i < ox->xmax ; i++)
        {width1[i] = i ;
         width2[i] = i+1 ;
        }
      for (j = 0 ; j < ox->ymax ; j++)
        {height1[j] = j ;
         height2[j] = j+1 ; 
        }
      
    }

  wtot = width2[ox->xmax - 1] ;
  htot = height2[ox->ymax - 1] ;
  xscale = ox->wide/wtot ;
  yscale = ox->high/htot ;

  /* add the length of the map names to the offset values */

  xlen = 0 ;
  for (i = 0 ; i < arrayMax (ox->axis2) ; i++)
    {
      l = strlen (name (array (ox->axis2, i, KEY))) ;
      if (l > xlen) xlen = l ;
    }
  xOff += xlen ;

  ylen = 0 ;
  for ( i = 0 ; i < arrayMax (ox->axis1) ; i++)
    {l = strlen (name (array (ox->axis1, i, KEY))) ;
     if (l > ylen) ylen = l ;
     }
  ylen *= 0.65 ;
  yOff += ylen + 4 ;

 
  /*  set up the clickable boxes for the two species names */
  /*  box values MIN_LIVE_BOX to MIN_LIVE_BOX + 1          */
  /*  b values 0 to 1                                      */
  
  ox->boxes = arrayCreate (128, BOX) ;
  b = 0 ;
  box = MIN_LIVE_BOX ;
  bx = arrayp (ox->boxes, b, BOX) ;
  bx->box = box ;
  bx->key = ox->s1 ;
  bx->x = 0 ;
  bx->y = 0 ;
  if (box != graphBoxStart ())
    messcrash ("minLiveBox wrong in oxgridDraw () for box %d", box) ;
  graphTextFormat (BOLD) ;
  graphText (name (ox->s1), xOff + (ox->wide)/2, 8) ;
  graphBoxEnd () ;
  graphBoxDraw (box, BLACK, WHITE) ;
  b++ ;
  box++ ;
  bx = arrayp (ox->boxes, b, BOX) ;
  bx->box = box ;
  bx->key = ox->s2 ;
  bx->x = 0 ;
  bx->y = 0 ;
    if (box != graphBoxStart ())
    messcrash ("minLiveBox wrong in oxgridDraw () for box %d", box) ;
  graphTextFormat (BOLD) ;
  graphText (name (ox->s2), 2, yOff + (ox->high)/2) ;
  graphBoxEnd () ;
  graphBoxDraw (box, BLACK, WHITE) ;
  b++ ;
  graphTextFormat (PLAIN_FORMAT);
  graphTextHeight (0.5) ;                                       
  graphColor (BLACK) ; 

  /*  set up the clickable boxes for the maps on the x-axis */
  /*  box values MIN_LIVE_BOX + 2 to                        */
  /*             MIN_LIVE_BOX + ox->xmax + 1                */
  /*  b values 2 to (ox->max + 1)                           */
  
  for (i = 0 ; i < ox->xmax ; i++)
    {box = MIN_LIVE_BOX + i + 2 ;
     bx = arrayp (ox->boxes, b , BOX) ;
     bx->box = box ;
     bx->key = array (ox->axis1, i, KEY) ;
     bx->x = 0 ;
     bx->y = 0 ;
     b++ ; 
     if (box != graphBoxStart ())
       messcrash ("minLiveBox wrong in oxgridDraw () for box %d, key = %s", box, name(bx->key)) ;
     japanese (name (array (ox->axis1, i, KEY)), xOff+0.5*(width1[i]+width2[i])*xscale-0.5, yOff-2-ylen) ;
     graphBoxEnd () ;
     graphBoxDraw (box, BLUE, WHITE) ;
    }


  /*  set up the clickable boxes for the maps on the y-axis */
  /*  box values MIN_LIVE_BOX + ox->xmax + 2 to             */
  /*             MIN_LIVE_BOX + ox->xmax + ox->ymax + 1     */
  /*  b values (ox->xmax + 2) to (ox->xmax + ox->ymax + 1)  */

  for (j = 0 ; j < ox->ymax ; j++)
    {box = MIN_LIVE_BOX + ox->xmax + j + 2 ;
     bx = arrayp (ox->boxes, b , BOX) ;
     bx->box = box ;
     bx->key = array(ox->axis2, j, KEY) ;
     bx->x = 0 ;
     bx->y = 0 ;
     b++ ; 
     if (box != graphBoxStart ())
       messcrash ("minLiveBox wrong in oxgridDraw () for box %d, key = %s", box, name(bx->key)) ;
     graphText (messprintf ("%s",name (array (ox->axis2, j, KEY))), xOff-xlen-0.5, yOff + (height1[j]+height2[j])*0.5*yscale - 0.5) ;
     graphBoxEnd () ;
     graphBoxDraw (box, BLUE, WHITE) ;
    }


  /*  set up the clickable boxes for grid cells                             */
  /*  box values MIN_LIVE_BOX + ox->xmax + ox->ymax + 2 to                  */
  /*             MIN_LIVE_BOX + ox->xmax + ox->ymax + ox->ymax*ox->xmax + 1 */
  /*  b values (ox->xmax + ox->ymax + 2) to                                 */
  /*           (ox->xmax + ox->ymax + ox->ymax*ox->xmax + 1)                */

  ox->map_labels = ox->xmax + ox->ymax + 2;
  ox->mlb = MIN_LIVE_BOX ; 
   
  count = 0 ;
  for (i = 0 ; i < ox->xmax ; i++)
    {for ( j = 0 ; j < ox->ymax ; j++)
       {box = MIN_LIVE_BOX + ox->map_labels + ox->ymax*i + j ;
        bx = arrayp (ox->boxes, b , BOX) ;

        dots = 0 ; ldots = 0 ;
        for (d = 0 ; d < arrayMax (ox->segs) ; d++) 
          {seg = arrayp (ox->segs, d, SEG) ;
           if (seg->box == box && !(seg->flag & FLAG_HIDE)) 
             {if (seg->flag & FLAG_LIGHT) ldots++ ;
              else dots++ ;
             }
          }
        bx->box = box ;
        bx->map1 = array (ox->axis1, i, KEY) ;
        bx->map2 = array (ox->axis2, j, KEY) ;
        bx->x = i ;
        bx->y = j ;
        count += dots ;
        count += ldots ;
        b++ ; 

        if (box != graphBoxStart ())
          messcrash ("minLiveBox wrong in oxgridDraw () box %d", box) ;
        graphRectangle (xOff+(width1[i]*xscale), yOff+(height1[j]*yscale), xOff+(width2[i]*xscale), yOff+(height2[j]*yscale)) ;

	if (!ox->hideFish)
	  {graphColor (YELLOW) ;
	   for (d = 0 ; d < arrayMax (ox->fishes) ; d++) 
	     {seg = arrayp (ox->fishes, d, SEG) ;
	      if (seg->box == box)
		{t1 = seg->z1 ; t2 = seg->z2 ; t3 = seg->z3; t4 = seg->z4 ;
		 fx1 = (xOff + (width1[i]*xscale) + (t1*(width2[i]-width1[i])*xscale)) ;
		 fx2 = (xOff + (width1[i]*xscale) + (t3*(width2[i]-width1[i])*xscale)) ;
		 fy1 = (yOff + (height1[j]*yscale) + (t2*(height2[j]-height1[j])*yscale)) ; 
		 fy2 = (yOff + (height1[j]*yscale) + (t4*(height2[j]-height1[j])*yscale)) ; 
		 graphFillRectangle (fx1+0.15, fy1+0.10, fx2-0.15, fy2-0.10) ;
	        }
	     }
	  }
	graphColor (BLACK) ;

	if (ox->showRandom)
	  {graphColor(BLACK) ;
	   for (d = 1 ; d <= dots ; d++)
	     {t1 = rand() ; t3 = (t1/RAND_MAX) ;
	      t2 = rand() ; t4 = (t2/RAND_MAX) ;
	      spx = (xOff+(width1[i]*xscale) + (t3*(width2[i]-width1[i])*xscale)) ;
	      spy = (yOff+(height1[j]*yscale) + (t4*(height2[j]-height1[j])*yscale)) ; 
	      if (t3 < 0.5) spx += 0.15 ;
	      if (t3 > 0.5) spx -= 0.15 ;
	      if (t4 < 0.5) spy += 0.15 ;
	      if (t4 > 0.5) spy -= 0.15 ;
	      graphFillRectangle (spx-0.15, spy-0.15, spx+0.15, spy+0.15) ; 
	     }
	   graphColor (RED) ;
	   for (d = 1 ; d <= ldots ; d++)
	     {t1 = rand() ; t3 = (t1/RAND_MAX) ;
	      t2 = rand() ; t4 = (t2/RAND_MAX) ;
	      spx = (xOff+(width1[i]*xscale) + (t3*(width2[i]-width1[i])*xscale)) ;
	      spy = (yOff+(height1[j]*yscale) + (t4*(height2[j]-height1[j])*yscale)) ; 
	      if (t3 < 0.5) spx += 0.15 ;
	      if (t3 > 0.5) spx -= 0.15 ;
	      if (t4 < 0.5) spy += 0.15 ;
	      if (t4 > 0.5) spy -= 0.15 ;
	      graphFillRectangle (spx-0.15, spy-0.15, spx+0.15, spy+0.15) ; 
	    }
	 }
	else
	  {for (d = 0 ; d < arrayMax (ox->segs) ; d++) 
	     {seg = arrayp (ox->segs, d, SEG) ;
	      if (seg->box == box && !(seg->flag & FLAG_HIDE))
		{t3 = seg->z1 ; t4 = seg->z2 ;

		 if (ox->flag & FLAG_EQUAL)
		   {t1 = array(ox->extent1, 2*i, float) ;
		    t2 = array(ox->extent1, 2*i+1, float) ; 
		    spx = (xOff + (width1[i]*xscale) + (t3/(t2-t1)*xscale)) ;
		    t1 = array(ox->extent2, 2*j, float) ;
		    t2 = array(ox->extent2, 2*j+1, float) ; 
		    spy = (yOff + (height1[j]*yscale) + (t4/(t2-t1)*yscale)) ; 
		  }
		 else 
		   {spx = (xOff + (width1[i]*xscale) + (t3*xscale)) ;
		    spy = (yOff + (height1[j]*yscale) + (t4*yscale)) ; 
		   }

		 if (phenotype_cat == 0 || (seg->phe1 == 0 && seg->phe2 == 0) || (!(seg->phe1 == phenotype_cat) && !(seg->phe2 == phenotype_cat))) 
		    {graphColor (BLACK) ;
		     graphFillRectangle (spx-0.15, spy-0.15, spx+0.15, spy+0.15) ; 
		    }

		 else if (seg->phe1 == phenotype_cat && seg->phe2 == phenotype_cat)
		   {graphColor (RED) ;
		    graphFillRectangle (spx-0.30, spy-0.30, spx+0.30, spy+0.30) ; 
		   }

		 else if (seg->phe1 == phenotype_cat && !(seg->phe2 == phenotype_cat))
		   {graphColor (GREEN) ;
		    graphFillRectangle (spx-0.30, spy-0.30, spx+0.30, spy+0.30) ; 
		   }

		 else if (!(seg->phe1 == phenotype_cat) && seg->phe2 == phenotype_cat)
		   {graphColor (BLUE) ;
		    graphFillRectangle (spx-0.30, spy-0.30, spx+0.30, spy+0.30) ; 
		   }

		 else printf ("Error in Oxford Grid - please report\n") ;
	        }
	     }
	  }
        graphColor (BLACK) ;
        graphBoxEnd () ;
	graphBoxDraw (box, BLACK, TRANSPARENT) ;
      }
   }

  graphTextHeight (1) ;
  graphText (messprintf ("%d homologies found",count), x, 1) ;       

  ox->gridBox = 0 ;

  {
    int 
      maxX = 50 + xOff+(width2[ox->xmax]*xscale),
      maxY = 50 + yOff+(height2[ox->ymax]*xscale) ;
    
    graphTextBounds (maxX, maxY) ;
  }
  
  graphRedraw () ;
}

/*************************************************************/

void japanese (char *text, float x, float y)
{
  int i, n = strlen(text) ;
  char a[2] ;

  a[1] = 0 ;
  for (i = 0 ; i < n ; ++i)
    { a[0] = *text++ ;
      graphText (a, x, y) ;
      y += 0.65 ;
    }
}

static void oxgridZoomOut (void) 
{
  OXGET ("oxgridZoomOut") ;

  if (ox->wide >= 30 && ox->high >= 20)
    {ox->wide /= 1.414 ;
     ox->high /= 1.414 ;

     oxgridDraw (ox) ;
    }

}

static void oxgridZoomIn (void) 
{
  OXGET ("oxgridZoomIn") ;

  if (ox->wide <= 1200 && ox->high <= 800)
    {ox->wide *= 1.414 ;
     ox->high *= 1.414 ;

     oxgridDraw (ox) ;
    }


}

static void oxgridReorg (void)
{
  int i, j ;
  int SAB, SA, SB ;
  float q = 0, q1, q2, q3 ;
  SEG *seg1, *seg2 ;

  OXGET ("oxgridReorg") ;

  SAB = 0 ; SA = 0 ; SB = 0 ;
  for ( i = 0 ; i < (arrayMax (ox->segs) - 1) ; i++)
    {
    for (j = 1 ; j < arrayMax (ox->segs) ; j++)
      {seg1 = arrayp (ox->segs, i, SEG) ;
       seg2 = arrayp (ox->segs, j, SEG) ;
       if (!(seg1->flag & FLAG_HIDE) &&
	   !(seg2->flag & FLAG_HIDE) &&
	   !(seg2->flag & FLAG_FISH) &&
	   !(seg2->flag & FLAG_FISH) )
        {if ((seg1->x == seg2->x) && (seg1->y == seg2->y)) 
           SAB++ ;
         if ((seg1->x == seg2->x) && !(seg1->y == seg2->y))
           SA++ ;
         if (!(seg1->x == seg2->x) && (seg1->y == seg2->y))
           SB++ ;
         }
       }
     }

  SA += SAB ; SB += SAB ;
  q1 = sqrt(SA) ;
  q2 = sqrt(SB) ;
  q3 = q1*q2 ;
  q = 1 - (SAB / q3) ;
  messout ("genome reorganisation measure q = %4.3f", q) ;

}

extern int ksetClassComplete (char *text, int len, int classe) ;

static int homolCompletion (char *cp, int len)
{
  return ksetClassComplete (cp, len, _VHomology_group) ;
}

static int locusCompletion (char *cp, int len)
{
  return ksetClassComplete (cp, len, _VLocus) ;
}


static void oxFindHom (char* string)
{
  int i ;
  KEY key ;
  SEG *seg ;
  BOOL found = FALSE ;

  OXGET ("oxFindHom") ;

  if (!(lexword2key (string, &key, _VHomology_group)))
    {messout ("Sorry, %s is not a homology name", string) ;
     return ;
    }

  for (i = 0 ; i < arrayMax (ox->segs) ; i++)
    {seg = arrayp (ox->segs, i, SEG) ;
     if (chkSpecies (seg->key, key)) 
       {found = TRUE ;
        oxgridSelect (ox, seg->box) ;
       }
     }

  if (found == FALSE) 
   {messout ("Sorry, this homology is not in this grid") ;
    return ;
   }

}

static void oxFindLoc (char* string)
{
  int i ;
  KEY key ;
  SEG *seg ;
  BOOL found = FALSE ;

  OXGET ("oxFindLoc") ;

  if (!(lexword2key (string, &key, _VLocus)))
    {messout ("Sorry, %s is not a locus name", string) ;
     return ;
    }

  for (i = 0 ; i < arrayMax (ox->segs) ; i++)
    {seg = arrayp (ox->segs, i, SEG) ;
     if (chkSpecies (seg->gene1, key) || chkSpecies (seg->gene2, key)) 
       {found = TRUE ;
        oxgridSelect (ox, seg->box) ;
       }
     }

  if (found == FALSE) 
   {messout ("Sorry, this locus is not in this grid") ;
    return ;
   }

}
 
static void toggleCells (void)
{

  OXGET ("toggleCells") ;

  if (ox->flag & FLAG_NOSIZE)
    {messout ("Maps have not been assigned relative sizes") ;
     return ;}
  else
    {ox->flag ^= FLAG_EQUAL ;
     oxgridDraw (ox) ;}

}


/**********************************/

static void oxHide (void)
{
  int i, j ;
  SEG *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;
 
  OXGET ("oxHide") ;

  if (!keySetActive (&keySet, &dummy))
    {messout ("First select a keySet window, thank you.") ;
     return ;
    }

  seg = arrp (ox->segs, 0, SEG) ;
  i = arrayMax (ox->segs) ;
  while(i--)
    {seg->flag &= ~FLAG_HIDE ; 
     if (!keySetFind (keySet, seg->key, &j) &&
         !keySetFind (keySet, seg->gene1, &j) &&
         !keySetFind (keySet, seg->gene2, &j)) 
       seg->flag |= FLAG_HIDE ;
     seg++ ;
    }

  oxgridDraw (ox) ;

}

static void oxClear (void) 
{
  int i ;
  SEG *seg ;
 
  OXGET ("oxClear") ;

  seg = arrayp (ox->segs, 0, SEG) ;
  i = arrayMax (ox->segs) ;
  while(i--)
    {seg->flag &= ~FLAG_HIDE ;  /* bitwise NOT */
     seg->flag &= ~FLAG_LIGHT ;
     seg++ ;
    }

  oxgridDraw (ox) ;

}

static void oxLightCell (void) 
{
  int i ;
  SEG *seg ;
  
  OXGET ("oxLightCell") ;

  if (!(ox->gridBox))
    {messout ("You have not selected a cell") ;
     return ;}

  for (i = 0 ; i < arrayMax (ox->segs) ; i++)
    {seg = arrayp (ox->segs, i, SEG) ;
     if (seg->box == ox->gridBox) seg->flag |= FLAG_LIGHT ;
    }

  oxgridDraw (ox) ;
}

static void oxUnLightCell (void)
{
  int i ;
  SEG *seg ;
  
  OXGET ("oxUnLightCell") ;

  if (!(ox->gridBox))
    {messout ("You have not selected a cell") ;
     return ;}

  for (i = 0 ; i < arrayMax (ox->segs) ; i++)
    {seg = arrayp (ox->segs, i, SEG) ;
     if (seg->box == ox->gridBox) seg->flag &= ~FLAG_LIGHT ;
    }

  oxgridDraw (ox) ;
}


static void oxLightAll (void) 
{
  int i ;

  OXGET ("oxLightAll") ;

  for (i = 0 ; i < arrayMax (ox->segs) ; ++i)
    arrp (ox->segs, i, SEG)->flag |= FLAG_LIGHT ;
    
  ox->gridBox = 0 ;

  oxgridDraw (ox) ;

}

static void oxLight (void) 
{
  int i, j ;
  SEG *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;

  OXGET ("oxLight") ;

  if (!keySetActive (&keySet, &dummy))
    { messout ("First select a keySet window, thank you.") ;
      return ;
    }
  for (i = 0 ; i < arrayMax (ox->segs) ; ++i)
    { seg = arrp (ox->segs, i, SEG) ;
      if (keySetFind (keySet, seg->key, &j) ||
          keySetFind (keySet, seg->gene1, &j) ||
          keySetFind (keySet, seg->gene2, &j)) 
	seg->flag |= FLAG_LIGHT ;
    }

  oxgridDraw (ox) ;

}

static void oxUnLight (void) 
{
  int i, j ;
  SEG *seg ;
  void *dummy ;
  KEYSET keySet = 0 ;

  OXGET ("oxUnLight") ;

  if (!keySetActive (&keySet, &dummy))
    { messout ("First select a keySet window, thank you.") ;
      return ;
    }
  for (i = 0 ; i < arrayMax (ox->segs) ; ++i)
    { seg = arrp (ox->segs, i, SEG) ;
      if (keySetFind (keySet, seg->key, &j) ||
          keySetFind (keySet, seg->gene1, &j) ||
          keySetFind (keySet, seg->gene2, &j)) 
	seg->flag &= ~FLAG_LIGHT ;
    }

  oxgridDraw (ox) ;

}

static void oxHideLight (void) 
{
  int i ;
  SEG *seg ;

  OXGET ("oxHideLight") ;

  for (i = 0 ; i < arrayMax (ox->segs) ; ++i)
    { seg = arrp (ox->segs, i, SEG) ;
      if (seg->flag & FLAG_LIGHT)
	seg->flag |= FLAG_HIDE ;
    }
  oxgridDraw (ox) ;

}

static void oxUnHide (void) 
{
  int i ;
  SEG *seg ;

  OXGET ("oxUnHide") ;

  for (i = 0 ; i < arrayMax (ox->segs) ; ++i)
    { seg = arrp (ox->segs, i, SEG) ;
      seg->flag &= ~FLAG_HIDE ;
    }
  oxgridDraw (ox) ;
}

static void oxExport (void) 
{
  int i, j, k = 0 ;
  SEG *seg ;
  static KEYSET kset = 0 ;
  BOOL found = FALSE ;

  OXGET ("oxExport") ;

  kset = keySetReCreate (kset) ;

  for (i = 0 ; i < arrayMax (ox->segs) ; ++i)
    { found = FALSE ;
      seg = arrp (ox->segs, i, SEG) ;
      if (seg->flag & FLAG_LIGHT)
        {for (j = 0 ; j < keySetMax (kset) ; j++)
          {if (chkSpecies (seg->key, keySet (kset, j))) found = TRUE ;}
         if (found == FALSE) keySet (kset, k++) = seg->key ;
        }
     }
  keySetSort (kset) ;
  keySetCompress (kset) ;
  keySetNewDisplay(kset,"Ox") ;
}

/**********************************/

static void oxClear2 (void) 
{
 
  OXGET ("oxClear2") ;

  phenotype_cat = 0 ;

  oxgridDraw (ox) ;

}

static void oxEye1 (void)
{

  OXGET ("oxEye1") ;

  phenotype_cat = 1 ;

  oxgridDraw (ox) ;

}

static void oxEye2 (void)
{

  OXGET ("oxEye2") ;

  phenotype_cat = 2 ;

  oxgridDraw (ox) ;

}

static void oxEye3 (void)
{

  OXGET ("oxEye3") ;

  phenotype_cat = 3 ;

  oxgridDraw (ox) ;

}

static void oxEye4 (void)
{

  OXGET ("oxEye4") ;

  phenotype_cat = 4 ;

  oxgridDraw (ox) ;

}

static void oxEye5 (void)
{

  OXGET ("oxEye5") ;

  phenotype_cat = 5 ;

  oxgridDraw (ox) ;

}

static void oxEye6 (void)
{

  OXGET ("oxEye6") ;

  phenotype_cat = 6 ;

  oxgridDraw (ox) ;

}

static void oxEye7 (void)
{

  OXGET ("oxEye7") ;

  phenotype_cat = 7 ;

  oxgridDraw (ox) ;

}

static void oxEye8 (void)
{

  OXGET ("oxEye8") ;

  phenotype_cat = 8 ;

  oxgridDraw (ox) ;

}

static void oxEye9 (void)
{

  OXGET ("oxEye9") ;

  phenotype_cat = 9 ;

  oxgridDraw (ox) ;

}

static void oxEye10 (void)
{

  OXGET ("oxEye10") ;

  phenotype_cat = 10 ;

  oxgridDraw (ox) ;

}

static void oxEye11 (void)
{

  OXGET ("oxEye11") ;

  phenotype_cat = 11 ;

  oxgridDraw (ox) ;

}

/*******************************************************************/


static void oxgridDump (void)	/* for debugging */
{
  int i ;
  SEG *seg ;

  OXGET ("oxgridDump") ;

  if (!(f = filopen ("oxgrid", "txt", "w")))
    { messout ("failed to open file oxgrid.txt") ;
      return ;
    }

  fprintf (f, "Dump of %s-%s Oxford Grid\n",name (ox->s1), name (ox->s2)) ;
  fprintf (f, "  Flag is %d\n", ox->flag) ;
  fprintf (f, "  Segs:\n") ;
  fprintf (f, "       Group                 Map1       Locus1") ;
  fprintf (f, "                 Map2       Locus2   Box \n\n") ;
  for (i = MIN_LIVE_BOX + ox->map_labels ; i < arrayMax (ox->segs) ; ++i)
    {seg = arrp (ox->segs, i, SEG) ;
     fprintf (f, "  %10s %20s %12s %20s %12s %5d \n", name (seg->key), name (seg->map1), name (seg->gene1), name (seg->map2), name (seg->gene2), seg->box) ;
    }

  filclose (f) ;

}

/*******************************************************************/

static void oxgridFish (void)	/* for debugging */
{
  OXGET ("oxgridDump") ;

  ox->hideFish = ! ox->hideFish ;
  oxgridDraw (ox) ;
}

/*******************************************************************/

static void oxgridRandom (void)	/* for debugging */
{
  OXGET ("oxgridDump") ;

  ox->showRandom = ! ox->showRandom ;
  oxgridDraw (ox) ;
}

/*******************************************************************/
/*******************************************************************/
