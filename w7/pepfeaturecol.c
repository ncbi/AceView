/* $Id: pepfeaturecol.c,v 1.5 2017/03/18 15:30:57 mieg Exp $ */
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
 * Last edited: Dec 21 13:49 1998 (fw)
 * Created: Tue 24 16:23 1996 (il)
 *-------------------------------------------------------------------
 */

/* model additions
under ?Colomns put 		pepFeature	FEA_bump  // il 27/9/96 
                                                FEA_Query Text
*/

#include "pepdisp.h"
#include "bump.h"
#include "query.h"

/************************************************************/

typedef struct featureColPriv {
  BOOL showNames;
  Array segs;
  BOOL bump;
  char *query;
  COND localQuery;
} *FEATUREPRIV;

typedef struct featuredata {
  KEY key;     /* key of the method */
  int x1,x2;   /* range of the feature */
  float score; 
} featureData;

/************************************************************/

static void featureDraw (COLINSTANCE instance, float *offset) 
{ 
  int i, box, xoff ;
  featureData *seg ;
  float y1, y2, dx, fmax = 0, wid = 0 ;
  METHOD *meth ;
  BUMP bump = 0 ;
  PEPLOOK *look = (PEPLOOK*) instance->map->look ;
  FEATUREPRIV private = instance->private;
  COLCONTROL control = instance->map->control ;
  Array segs;
  OBJ matchObj;

  currentMapControl () ; 
  
/*  segs = (Array) (FeatureConvert (currentMap, NULL)); */
  segs =  * (instance->proto->convertResults);
  if (private->bump)
    bump = bumpCreate (1000, 0) ;
  for (i = 0 ; i < arrayMax (segs) ; ++i)
    { seg = arrp (segs, i, featureData) ;
      matchObj = 0;
      if ((private->localQuery) && 
	  !queryFind3 (private->localQuery, &matchObj, seg->key))
	{ bsDestroy (matchObj);
	  continue;
	}
      bsDestroy (matchObj);
      meth = method (0, seg->key) ;
      
      /* checks to prevent arithmetic crash */
      if (meth->flags & METHOD_SCORE_BY_OFFSET && !meth->min)
	meth->min = 1 ;
      if (meth->flags & METHOD_SCORE_BY_WIDTH && meth->max == meth->min)
	meth->max = meth->min + 1 ;
      
      if (meth->width)
	wid = meth->width/2 ;
      else
	wid = 1 ;
      
      y1 = MAP2GRAPH (look->map, seg->x1) ;
      y2 = MAP2GRAPH (look->map, seg->x2) ;
      if (y1 > control->graphHeight -0.5 || y2 < control->topMargin+0.5)
        continue ;
      
      if (y1 < control->topMargin+0.5)
        y1 = control->topMargin+0.5 ;

      if (y2 > control->graphHeight -0.5)
        y2 = control->graphHeight -0.5 ;

      box = graphBoxStart () ;
      if (meth->flags & METHOD_SCORE_BY_OFFSET)
        { dx = 6 - 4*log10 ((double) (seg->score/meth->min)) ;
          if (dx < 0)
            dx = 0 ;
	  xoff = dx;
	  if (bump)
            {
	      if (look->map->mag < 0)
		{ y1 = -y1 ; y2 = -y2 ;
                bumpItem (bump, 2, (y2-y1), &xoff, &y1) ;
                y1 = -y1 ; y2 = -y2 ;
		}
	      else
		bumpItem (bump, 2, (y2-y1), &xoff, &y1) ;
	    }
          graphRectangle (*offset + wid* (xoff+dx), y1, 
                          *offset + wid* (xoff+dx + 0.85), y2) ;
          if (fmax <  wid* (xoff+dx + 0.85))
            fmax =  wid* (xoff+dx + 0.85);
        }
      else 
        { if (meth->flags & METHOD_SCORE_BY_WIDTH)
            { dx = 0.25 + 0.75* (seg->score - meth->min) / 
                                        (meth->max - meth->min) ;
              if (dx < 0.25) dx = 0.25 ;
              if (dx > 1) dx = 1 ;
            }
          else
            dx = 0.75 ;
          xoff = 1 ;
          if (bump)
            {
	      if (look->map->mag < 0)
		{ y1 = -y1 ; y2 = -y2 ;
                bumpItem (bump, 2, (y2-y1), &xoff, &y1) ;
                y1 = -y1 ; y2 = -y2 ;
		}
	      else
		bumpItem (bump, 2, (y2-y1), &xoff, &y1) ;
	    }
          graphRectangle (*offset + wid* (xoff - dx), y1, 
                          *offset + wid* (xoff + dx), y2) ;
          if (fmax < xoff) fmax = xoff ;
        }

      graphBoxEnd () ;
      graphBoxDraw (box, BLACK, meth->col) ;
      controlRegBox (instance, box, seg) ;
    }
  
  if (bump) bumpDestroy (bump) ;
  if (fmax) *offset += wid * (fmax + 1) ;
}
static void featurePrivDestroy (void * p)
{
 FEATUREPRIV private = (FEATUREPRIV)p;

 if (arrayExists (private->segs)) arrayDestroy (private->segs);
}

static void pepFeatureSave (COLINSTANCE instance, OBJ init) 
{
  FEATUREPRIV private = instance->private;
  
  if (private->bump)
    bsAddTag (init,str2tag ("FEA_bump"));
  if (private->query)
    if (!bsAddData (init,str2tag ("FEA_Query"), _Text, private->query))
      printf ("Error saving query %s\n",private->query);

}
/******************************************************************************************************/
struct configLocals{
  BOOL bump;
  char query[280];
};

static BOOL featureConfigure (COLINSTANCE instance)     
{
  struct  configLocals *cf = (struct configLocals *) messalloc (sizeof (struct configLocals));
  FEATUREPRIV private = instance->private;

  if (controlCreateConfig (instance,cf,"Configure Feature column",0.5,0.15)){

    cf->bump = private->bump;
    if (private->query)
      sprintf (cf->query,"%s",private->query);

    graphToggleEditor ("Bump",&cf->bump,4,8);
    graphTextEditor ("Query for display:",cf->query,280,4,6,0);
    graphRedraw ();
  }
  return FALSE;
}

static void featureConfigFinal (COLINSTANCE instance, void *locals, BOOL ok)
{ struct configLocals *cf = locals;
  FEATUREPRIV private = instance->private;
  
  if (ok)
    { private->bump = cf->bump;
      private->query = 0;
      private->localQuery = 0;
      if (strlen (cf->query) != 0){
	 COND c;
	 private->query = strnew (cf->query, instance->handle);
	 if (condConstruct (private->query, &c))
	   private->localQuery = c;
       }
    }
  else
    messfree (cf);
}

/************************************pephomolcol.c:539: warning: unused variable `i'
pephomolcol.c: In function `homolDoColour':
*******************************************************************/
static BOOL featureUnSetSelect (COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;

  control->activeKey = 0;
  look->startcol = 0;
  look->endcol = 0;
  
  *look->messageText = 0;
  
  if (look->messageBox)
    graphBoxDraw (look->messageBox,-1,-1);

  return FALSE;
}
/****************************************************************************************************/
static BOOL featureSetSelect (COLINSTANCE instance, int box, double x, double y)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;
  featureData *feature = (featureData *) arr (control->boxIndex2,box,void *);

  control->activeKey = feature->key;
  control->activeBox = box;
  look->startcol = feature->x1;
  look->endcol = feature->x2;

  return FALSE;

}
/****************************************************************************************************/
static void featureDoColour (COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;
  featureData *feature = (featureData *) arr (control->boxIndex2,box,void *);
  METHOD *meth = method (0, feature->key) ;

  if (control->activeKey == feature->key && 
      look->startcol == feature->x1 && look->endcol == feature->x2)
    {
      graphBoxDraw (box, BLACK,activecol);
      strncpy (look->messageText, name (feature->key),100);
      strcat (look->messageText, messprintf (" %d  %d ",feature->x1,feature->x2));
      strcat (look->messageText,messprintf (" %4.2f ",feature->score));
      if (look->messageBox)
	graphBoxDraw (look->messageBox,-1,-1);
      control->activeBox = box;
    }
  else
    if (look->endcol ==0 && 
       feature->x1 <= look->startcol && feature->x2 >= look->startcol)
      graphBoxDraw (box,BLACK,friendcol);
    else
      graphBoxDraw (box, BLACK,meth->col);
}

/****************************************************************************************************/

void *FeatureConvert (MAPCONTROL map, void *params)
{
  Array flatA;
  featureData *feature;
  Array segs;
  OBJ obj;
  int i;

  segs = arrayHandleCreate (20,featureData,map->handle);
  flatA = arrayCreate (4, BSunit);
  obj = bsCreate (map->key);
  if (bsFindTag (obj, str2tag ("Feature")) && bsFlatten (obj,4,flatA)){
    for (i=0;i<arrayMax (flatA); i += 4){
      feature = arrayp (segs, arrayMax (segs),featureData);
      feature->key = arr (flatA,i,BSunit).k;
      feature->x1 = arr (flatA,i+1,BSunit).i;
      feature->x2 = arr (flatA,i+2,BSunit).i;
      feature->score = arr (flatA,i+3,BSunit).f;
      methodAdd (0, feature->key);
    }
  }
  arrayDestroy (flatA);
  bsDestroy (obj);
  return segs;
}

static void featureFollowBox (COLINSTANCE instance, int box,double x, double y)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;
  featureData *feature = (featureData *) arr (control->boxIndex2,box,void *);

  if (look->block){
    look->activeStart = feature->x1;
    look->activeEnd = feature->x2;
    displayUnBlock ();
    look->block = FALSE;
    controlDraw ();
  }
}

extern  BOOL pepFeatureCreate (COLINSTANCE instance, OBJ init) 
{
  FEATUREPRIV private;
  char *s1;
  COND c;
  
  instance->draw = featureDraw ;
  instance->configure = featureConfigure;
  instance->doColour = featureDoColour; 
  instance->unSelectBox = featureUnSetSelect; 
  instance->setSelectBox =  featureSetSelect;
  instance->save = pepFeatureSave;
  instance->followBox = featureFollowBox;
  instance->configFinal = featureConfigFinal;

  private = (FEATUREPRIV) handleAlloc (featurePrivDestroy, instance->handle , sizeof (struct featureColPriv));
  instance->private = private;

  private->bump = FALSE;

  if (init){
    if (bsFindTag (init,str2tag ("FEA_bump")))
      private->bump = TRUE;
    
    
    if (bsGetData (init,str2tag ("FEA_Query"),_Text,&s1)){
      private->query = strnew (s1,instance->handle);
      if (condConstruct (private->query, &c))
	private->localQuery = c;
    }
  }
  return TRUE;
}
 
 
 
 
 
 
