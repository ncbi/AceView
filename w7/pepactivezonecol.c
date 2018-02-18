/*  File: pepfeaturecol.c
 *  Author: Ian Longden (il@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 *      Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *      Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: feature display for fmap package
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 11 10:36 1996 (srk)
 * Created: Tue 24 16:23 1996 (il)
 *-------------------------------------------------------------------
 */

/* $Id: pepactivezonecol.c,v 1.1.1.1 2002/07/19 20:23:00 sienkiew Exp $ */
/* model additions
under ?Colomns put 		pepActiveZone	
*/

#include "mystdlib.h"
#include "pepdisp.h"
BOOL alignDisplay (KEY key, KEY from, BOOL isOldGraph);

static void activeZoneDraw (COLINSTANCE instance, float *offset) 
{ 
  int box ;
  float y1, y2 ;
  PEPLOOK *look = (PEPLOOK*) instance->map->look ;
  COLCONTROL control = instance->map->control ;
      
  y1 = MAP2GRAPH(look->map, look->activeStart) ;
  y2 = MAP2GRAPH(look->map, look->activeEnd) ;
  if (y1 > control->graphHeight -0.5 || y2 < control->topMargin+0.5)
    return ;
  
  if (y1 < control->topMargin+0.5)
    y1 = control->topMargin+0.5 ;
  
  if (y2 > control->graphHeight -0.5)
    y2 = control->graphHeight -0.5 ;
  
  box = graphBoxStart() ;
  graphRectangle (*offset , y1, *offset+1.0, y2) ;
 
  graphBoxEnd() ;
  controlRegBox (instance, box, (void *)1) ;
  graphBoxDraw (box, BLACK, BLUE) ;
  *offset += 3 ;
}

/*******************************************************************************************************/
static BOOL activeZoneUnSetSelect(COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;

  control->activeKey = 0;
  look->startcol = 0;
  look->endcol = 0;
  
  *look->messageText = 0;
  
  if(look->messageBox)
    graphBoxDraw(look->messageBox,-1,-1);

  return FALSE;
}
/****************************************************************************************************/
static BOOL activeZoneSetSelect(COLINSTANCE instance, int box, double x, double y)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;

  control->activeKey = 0;
  control->activeBox = box;
  look->startcol = look->activeStart;
  look->endcol = look->activeEnd;

  return FALSE;

}
/****************************************************************************************************/
static void activeZoneDoColour(COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;

  if(instance == control->activeInstance){
    graphBoxDraw(box, BLACK,activecol);
  }
  else
      graphBoxDraw(box,BLACK,BLUE);
}

/****************************************************************************************************/

extern  BOOL activeZoneCreate (COLINSTANCE instance, OBJ init) 
{
  instance->draw = activeZoneDraw ;
  instance->doColour = activeZoneDoColour; 
  instance->unSelectBox = activeZoneUnSetSelect; 
  instance->setSelectBox =  activeZoneSetSelect;
/*  instance->save = activeZoneSave;*/

  return TRUE;
}
 
