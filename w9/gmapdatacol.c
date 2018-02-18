 /*  File: gmapdatacol.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Display columns for maping data
 * Exported functions:
 *      multiPtAddKey
 *      gMapRDMultiPtColumn
 *      gMapJTMMultiPtColumn
 *      twoPtAddKey
 *      gMapRDTwoPtColumn
 *      gMapJTMTwoPtColumn
 *      dbnAddKey
 *      gMapLikelihoodColumn
 * HISTORY:
 * Last edited: Dec 22 09:46 1998 (fw)
 * Created: Tue Nov 30 18:42:30 1993 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: gmapdatacol.c,v 1.1.1.1 2002/07/19 20:23:03 sienkiew Exp $ */

#include "gmap.h"

/********************************************************/
/************** Multi_pt_data ***************************/
typedef struct {
  KEY key;
  float x, dx;
  MULTIPTDATA data;
} MULTIPTSEG;

typedef struct MultiPrivStruct {
  Array segs;
} *MULTIPRIV;

void multiPtAddKey(COLINSTANCE instance, KEY key) 
{ MULTIPRIV private;
  MULTIPTSEG *seg;
  MULTIPTDATA data; 
  int i;
  
  private = (MULTIPRIV) instance->private;

  if (!key) /* clear */
    { arrayMax(private->segs) = 0;
      return;
    }

  /* If we've already got it, don't do it again */
  for (i=0; i<arrayMax(private->segs); i++)
    { seg = arrp(private->segs, i, MULTIPTSEG);
      if (seg->key == key)
	{
	  if (!instance->displayed)
	    instance->displayed = TRUE;
	  return;
	}
    }  

  if (!gMapGetMultiPtData(instance->map, key, instance->handle, &data))
    return;
    
  if (!instance->displayed)
    instance->displayed = TRUE;


  seg = arrayp(private->segs, arrayMax(private->segs), MULTIPTSEG);
  seg->key = key;
  seg->data = data;
  seg->x = data->min; 
  seg->dx = data->max - data->min;
  return;
}

/*******************/


static BOOL multiPtSetSelect (COLINSTANCE instance,
			      int box,
			      double x,
			      double y)
{
  COLCONTROL control = instance->map->control;
  MULTIPTSEG *seg = (MULTIPTSEG *)arr(control->boxIndex2, box, void *);
  GeneticMap look = (GeneticMap)(instance->map->look);

  look->selectKey = seg->key;
  look->multiPtCurrentColumn = instance;  

  /* We can only have neighbours */
  gMapUnselect(instance->map);
  (void)gMapNeighbours(instance->map, &look->neighbours, seg->key);
  look->neighboursInfoValid = TRUE;
  
  return FALSE;

}

static BOOL multiPtUnselect(COLINSTANCE instance, int box)
{
  GeneticMap look = (GeneticMap)(instance->map->look);

  *look->messageText = 0;
  if (look->messageBox)
    graphBoxDraw(look->messageBox, -1, -1);

  return FALSE;
}

static void multiPtDoColours(COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  MULTIPTSEG *seg = (MULTIPTSEG *) arr(control->boxIndex2, box, void *); 
  GeneticMap look = (GeneticMap)(instance->map->look);
  int i;
  float dummy;
  
  if (seg->key == look->selectKey && 
      instance == control->activeInstance) 
    { 
      graphBoxDraw(box, BLACK, CYAN);
      control->activeBox = box;  
      *look->messageText = 0 ;
      strcpy (look->messageText, "Multi_pt ") ;
      strcat (look->messageText, name(seg->key)) ;
      strcat (look->messageText, ": ") ;
      strcat (look->messageText, name(arr(seg->data->loci,0,KEY))) ;
      for (i = 1 ; i < arrayMax(seg->data->loci) ; ++i)
	strcat (look->messageText, 
		messprintf (" %d %s", arr(seg->data->counts,i-1,int), 
			    name(arr(seg->data->loci,i,KEY)))) ;
      if (look->messageBox)
	graphBoxDraw(look->messageBox, -1, -1);
  
    }
  else
    {
      if (control->activeInstance &&
	   map == control->activeInstance->map &&
	   look->dataInfoValid &&
	   keySetFind(look->multi, seg->key, 0))
	{ 
	  BOOL good = logLikeMulti(map, seg->data, seg->key, &dummy) ;
	  if (!instance->map->cambridgeOptions) /* more detailled colouring */
	    graphBoxDraw(box, BLACK, good ? YELLOW : PALERED) ;
	  else
	    graphBoxDraw(box, BLACK, good ? GREEN : DARKGREEN);
	}
      else if (keySetFind(look->highlight, seg->key, 0))
	graphBoxDraw(box, BLACK, MAGENTA);
      else
	graphBoxDraw(box, BLACK, WHITE);
    }
}


static void multiPtRDDraw(COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);
  MULTIPRIV private = instance->private;
  int   i, j, n;
  float y1, y2, x, gy1, gy2;
  int ibox;
  MULTIPTSEG *seg ;
  MULTIPTDATA data;


  for (i = 0 ; i < arrayMax(private->segs) ; ++i)
    { seg = arrp(private->segs,i,MULTIPTSEG) ;
      if (!keySetFind(look->hidden, seg->key, 0))
	{ ibox = graphBoxStart();
	  array(control->boxIndex, ibox, COLINSTANCE) = instance;
	  array(control->boxIndex2, ibox, void *) = (void *) seg;
	  data = seg->data;
	  n = arrayMax(data->loci);
	  y1=0, y2=0 ;
	  x = (*offset)++;
 
	  graphPointsize (0.8) ;
	  if (getPos (map, arr(data->loci, 0, KEY), &y1))
	    if (controlCheckY(map, y1, 0.4, &gy1))
	      graphPoint (x, gy1) ;
	  if (getPos (map, arr(data->loci, n-1, KEY), &y2))
	      if (controlCheckY(map, y2, 0.4, &gy2))
		graphPoint (x, gy2) ;
	  if (gy1 != gy2)
	    { graphLine (x, gy1, x, gy2);
	      remarkRegister(map, seg->key, y1);
	    }

	  graphPointsize (0.6) ;
	  for (j = 1 ; j < n ; ++j)	/* redo n-1 to get counts line */
	    if (getPos (map, arr(data->loci,j,KEY), &y2))
	      { 
		if (controlCheckY(map, y2, 0.3, &gy2))
		  graphPoint (x, gy2) ;
		if (arr (data->counts, j-1, int))
		  if (controlCheckY(map, (y1+y2)/2, 0, &gy1))
		    graphLine (x-0.25, gy1, x+0.25, gy1) ;
		y1 = y2 ;
	      }
	  graphBoxEnd() ;

	  if (seg->key == control->from)
	    control->fromBox = ibox;
	}
    }
  (*offset)++;
}

static void multiPtJTMDraw(COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);
  MULTIPRIV private = instance->private;
  int   i, j, n, sens;
  float y1, y2, ymin, ymax, x, gy1, gy2;
  MULTIPTSEG *seg ;
  MULTIPTDATA data;
  int ibox;

  for (i = 0 ; i < arrayMax(private->segs) ; ++i)
    { seg = arrp(private->segs,i,MULTIPTSEG) ;
      if (!keySetFind(look->hidden, seg->key, 0))
	{ ibox = graphBoxStart();
	  array(control->boxIndex, ibox, COLINSTANCE) = instance;
	  array(control->boxIndex2, ibox, void *) = (void *) seg;
	  data = seg->data;
	  n = arrayMax(data->loci);
	  y1=0, y2=0 ;
	  sens = 0;
	  x = (*offset)++;
	  
	  graphPointsize (0.8) ;
	  if (getPos (map, arr(data->loci, 0, KEY), &y1))
	    if (controlCheckY(map, y1, 0.4, &gy1))
	      graphPoint (x, gy1) ;
	  ymin = ymax = y1 ;
	  if (getPos (map, arr(data->loci, n-1, KEY), &y2))
	    if (controlCheckY(map, y2, 0.4, &gy2))
	      graphPoint (x, gy2) ;
	  
	  graphPointsize (0.6) ;
	  for (j = 1 ; j < n ; ++j)	/* redo n-1 to get counts line */
	    if (getPos (map, arr(data->loci,j,KEY), &y2))
	      { if (controlCheckY(map, y2, 0.3, &gy2))
		  graphPoint (x, gy2) ;
		if (y2 != y1 && gy2 != gy1 && arr (data->counts, j-1, int))
		  { 
		    if (sens)
		      { if ((y2-y1) * sens  < 0)
			  { graphColor (RED) ;
			    graphFillRectangle (x - .2, gy1, x + .2, gy2) ;
			  }
		      }
		    else
		      sens = y2 - y1  > 0 ? 1 : -1 ;
		    if (controlCheckY(map, (y1+y2)/2, 0, &gy1))
		      graphLine (x-0.25, gy1, x+0.25, gy1) ;
		    graphColor (BLACK) ;
		  }
		y1 = y2 ;
		controlCheckY(map, y1, 0, &gy1);
		if (y2 < ymin) ymin = y2 ;
		if (y2 > ymax) ymax = y2 ;
	      }
	  controlCheckY(map, ymin, 0, &gy1);
	  controlCheckY(map, ymax, 0, &gy2);
	  if (gy1 != gy2)
	    { graphLine (x, gy1, x, gy2);
	      remarkRegister(map, seg->key, ymin);
	    }
	  graphBoxEnd() ;

	  if (seg->key == control->from)
	    control->fromBox = ibox;
	}
    }
  (*offset)++;
}

static void multiPtFollowBox(COLINSTANCE instance, int box, double x, double y)
{
  MULTIPTSEG *seg = (MULTIPTSEG *) arr(instance->map->control->boxIndex2, 
				       box, void *); 
  
  display(seg->key, instance->map->key, TREE);
}

static void multiPtDraw(COLINSTANCE instance, float *offset)
{ 
  if (instance->map->cambridgeOptions)
    multiPtRDDraw(instance, offset);
  else
    multiPtJTMDraw(instance, offset);
}

static BOOL multiPtCreate(COLINSTANCE instance, OBJ init)
{
  GeneticMap look = (GeneticMap)(instance->map->look);
  MULTIPRIV private;

  private = (MULTIPRIV)halloc(sizeof(struct MultiPrivStruct),
			      instance->handle);

  instance->private = private;
  private->segs = arrayHandleCreate(50, MULTIPTSEG, instance->handle);

  instance->draw = multiPtDraw;
  instance->setSelectBox = multiPtSetSelect;
  instance->unSelectBox = multiPtUnselect;
  instance->doColour = multiPtDoColours;
  instance->followBox = multiPtFollowBox;

  look->multiPtCurrentColumn = instance;
  return TRUE;
}

static void multiPtDestroy(void *p)
{
  COLINSTANCE instance = (COLINSTANCE)p;
  GeneticMap look = (GeneticMap)(instance->map->look);
  
  if (look->multiPtCurrentColumn == instance)
    look->multiPtCurrentColumn = 0;
}

struct ProtoStruct gMapMultiPtColumn = {
  0,
  multiPtCreate,
  multiPtDestroy,
  "Multi_point",
  0,
  TRUE,
  0,0,0,0,
  0
};


/********************************************************/
/************************ 2 point ***********************/
typedef struct TwoPrivStruct {
  Array segs; /* private array of things to draw */
} *TWOPRIV;

typedef struct {
  float x, dx;
  KEY key;
  TWOPTDATA data;
} TWOPTSEG;


/***********/



static BOOL twoPtSetSelect (COLINSTANCE instance, int box, double x, double y)
{
  COLCONTROL control = instance->map->control;
  TWOPTSEG *seg = (TWOPTSEG *)arr(control->boxIndex2, box, void *);
  GeneticMap look = (GeneticMap)(instance->map->look);

  look->selectKey = seg->key;
  look->twoPtCurrentColumn = instance;


  gMapUnselect(instance->map);
  /* We can only have neighbours */
  (void)gMapNeighbours(instance->map, &look->neighbours, seg->key);
  look->neighboursInfoValid = TRUE;

  return FALSE;
}

static BOOL twoPtUnselect(COLINSTANCE instance, int box)
{
  GeneticMap look = (GeneticMap)(instance->map->look);

  *look->messageText = 0;
  if (look->messageBox)
    graphBoxDraw(look->messageBox, -1, -1);

  return FALSE;
}

static void twoPtFollowBox(COLINSTANCE instance, int box, double x, double y)
{
  TWOPTSEG *seg = (TWOPTSEG *) arr(instance->map->control->boxIndex2, 
				 box, void *); 
  
  display(seg->key, instance->map->key, TREE);
}

static void twoPtDoColour(COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  TWOPTSEG *seg = (TWOPTSEG *)arr(control->boxIndex2, box, void *);
  TWOPTDATA data;
  GeneticMap look = (GeneticMap)(map->look);
  int i, *p;
  float best, lo, hi, y1, y2 ;
  char *buf = look->messageText ;

  if (seg->key == look->selectKey &&
      instance == control->activeInstance)
    { graphBoxDraw(box, BLACK, CYAN);
      control->activeBox = box;
      data = seg->data;
      *buf = 0 ;
      
      strncpy (buf, messprintf ("2point %s: %s %s: %s",
				name(seg->key), 
				name(data->loc1), 
				name(data->loc2), 
				data->type ? name(data->type) : 
				"Simple distance" ), 100) ; /* TODO OK? */

      if (data->type)		     /* 27 chars plenty for counts */
	for (p = &data->n1, i = 0 ; i < data->count ; ++i, p++)
	  strcat (buf, messprintf (" %d", *p)) ;
      else
	strcat (buf, messprintf (" %.2f %.2f",
				 data->distance, 
				 data->error)) ;
      
      if (getPos (map, data->loc1, &y1) && getPos (map, data->loc2, &y2))
	strcat (buf, messprintf (": %.2f", (y1>y2) ? (y1-y2) : (y2-y1))) ;
      
      if (best2p (data, &best, &lo, &hi))
	strncat (buf, messprintf (" [%.2f, %.2f, %.2f]", 
				  lo, best, hi), 127 - strlen(buf)) ;
      
      if (look->messageBox)
 	graphBoxDraw(look->messageBox, -1, -1);
      
    }
  else
    { 
      if (control->activeInstance &&  
	  map == control->activeInstance->map &&
	  look->dataInfoValid &&
	  keySetFind(look->twoPt, seg->key, 0))
	{ if (instance->map->cambridgeOptions)
	    graphBoxDraw(box, BLACK, GREEN);
	  else
	    graphBoxDraw(box, BLACK, YELLOW);
	}
      else if (keySetFind(look->highlight, seg->key, 0))
	graphBoxDraw(box, BLACK, MAGENTA);
      else
	graphBoxDraw(box, BLACK, WHITE);
	
    }
}

static void twoPtRDDraw (COLINSTANCE instance, float *offset)
{ 
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  TWOPRIV private = instance->private;
  GeneticMap look = (GeneticMap)(map->look);
  int   i;
  TWOPTSEG   *seg ;
  int ibox;
  float x;
  float y1, y2, best, lo, hi, cen;
  float gy1, gy2;

  for (i = 0 ; i < arrayMax(private->segs) ; ++i)
    { seg = arrp(private->segs,i,TWOPTSEG) ;
      if (!keySetFind(look->hidden, seg->key, 0))
	{ ibox = graphBoxStart();
	  array(control->boxIndex, ibox, COLINSTANCE) = instance ;
	  array(control->boxIndex2, ibox, void *) = (void *)seg;
	  x = (*offset)++;

	  if (!getPos (map, seg->data->loc1, &y1) || 
	      !getPos (map, seg->data->loc2, &y2))
	    return ;

	  graphPointsize (0.7) ;
	  if (controlCheckY(map, y1, 0.35, &gy1))
	    graphPoint (x, gy1) ;
	  if (controlCheckY(map, y2, 0.35, &gy2))
	    graphPoint (x, gy2) ;
	  if (gy1 != gy2)
	    { graphLine (x, gy1, x, gy2); 
	      remarkRegister(map, seg->key, y1);
	    }

	  if (best2p (seg->data, &best, &lo, &hi))
	    { cen = (y1 + y2) / 2 ;
	      if (controlCheckY(map, cen-best/2, 0, &gy1))
		graphLine (x-0.3, gy1, x+0.3, gy1);
	      if (controlCheckY(map, cen+best/2, 0, &gy1))
		graphLine (x-0.3, gy1, x+0.3, gy1);

	      if (controlCheckY(map, cen-hi/2, 0, &gy1))
		graphLine (x-0.3, gy1, x+0.3, gy1);
	      if (controlCheckY(map, cen-lo/2, 0, &gy2))
		graphLine (x-0.3, gy2, x+0.3, gy2);
	      if (gy1 != gy2)
		{ graphLine(x-0.3, gy1, x-0.3, gy2);
		  graphLine(x+0.3, gy1, x+0.3, gy2);
		}
	      
	      if (controlCheckY(map, cen+hi/2, 0, &gy1))
		graphLine (x-0.3, gy1, x+0.3, gy1);
	      if (controlCheckY(map, cen+lo/2, 0, &gy2))
		graphLine (x-0.3, gy2, x+0.3, gy2);
	      if (gy1 != gy2)
		{ graphLine(x-0.3, gy1, x-0.3, gy2);
		  graphLine(x+0.3, gy1, x+0.3, gy2);
		}

	    }
	  
	  graphBoxEnd() ;
 	  if (seg->key == control->from)
	    control->fromBox = ibox;
	}
    }
  (*offset)++ ;
}

static void twoPtJTMDraw (COLINSTANCE instance, float *offset)
{ 
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  TWOPRIV private = instance->private;
  GeneticMap look = (GeneticMap)(map->look);
  int   i;
  TWOPTSEG   *seg ;
  int ibox;
  float x;
  float y1, y2, best, lo, hi, cen;
  float gy1, gy2, gy3, gy4;
  BOOL t1, t2, t3, t4;

  for (i = 0 ; i < arrayMax(private->segs) ; ++i)
    { seg = arrp(private->segs,i,TWOPTSEG) ;
      if (!keySetFind(look->hidden, seg->key, 0))
	{ ibox = graphBoxStart();
	  array(control->boxIndex, ibox, COLINSTANCE) = instance ;
	  array(control->boxIndex2, ibox, void *) = (void *)seg;
	  x = (*offset)++;
	  
	  if (!getPos (map, seg->data->loc1, &y1) || 
	      !getPos (map, seg->data->loc2, &y2))
	    return ;
	  
	  if (y1 > y2)
	    { float tmp = y1 ; 
	      y1 = y2 ;
	      y2 = tmp ;
	    }
	  
	  /* join the 2 genes */
	  t1 = controlCheckY(map, y1, 0, &gy1);
	  t2 = controlCheckY(map, y2, 0, &gy2);
	  if (gy1 != gy2)
	    { graphLine (x, gy1, x, gy2);
	      remarkRegister(map, seg->key, y1);
	    }
	  
	  if (best2p (seg->data, &best, &lo, &hi))
	    { cen = (y1 + y2) / 2 ;
	      
	      /* a red rectangle from true to best positions */
	      
	      t3 = controlCheckY(map, cen-best/2, 0, &gy3);
	      if (gy1 != gy3)
		{ graphColor(RED);
		  graphFillRectangle(x-.2, gy1, x+.2, gy3) ;
		  graphColor(BLACK);
		  if (t1) 
		    graphLine(x-.2, gy1, x+.2, gy1) ;
		  if (t3)
		    graphLine(x-.2, gy3, x+.2, gy3) ;
		  graphLine(x-.2, gy1, x-.2, gy3) ;
		  graphLine(x+.2, gy1, x+.2, gy3) ;
		}
	      t4 = controlCheckY(map, cen+best/2, 0, &gy4);
	      if (gy2 != gy4)
		{ graphColor(RED);
		  graphFillRectangle(x-.2, gy2, x+.2, gy4) ;
		  graphColor(BLACK);
		  if (t2) 
		    graphLine(x-.2, gy2, x+.2, gy2) ;
		  if (t4)
		    graphLine(x-.2, gy4, x+.2, gy4) ;
		  graphLine(x-.2, gy2, x-.2, gy4) ;
		  graphLine(x+.2, gy2, x+.2, gy4) ;
		}

	      /* a green error box around best position, hence hiding the red box */
	     
	      t3 = controlCheckY(map, cen - lo/2, 0, &gy3);
	      t4 = controlCheckY(map, cen - hi/2, 0, &gy4);
	      if (gy3 != gy4)
		{ graphColor(GREEN) ;
		  graphFillRectangle(x-.25, gy3, x+.25, gy4);
		  graphColor(BLACK);
		  if (t3) 
		    graphLine(x-.2, gy3, x+.2, gy3) ;
		  if (t4)
		    graphLine(x-.2, gy4, x+.2, gy4) ;
		  graphLine(x-.2, gy3, x-.2, gy4) ;
		  graphLine(x+.2, gy3, x+.2, gy4) ;
		}
	      t3 = controlCheckY(map, cen + lo/2, 0, &gy3);
	      t4 = controlCheckY(map, cen + hi/2, 0, &gy4);
	      if (gy3 != gy4)
		{ graphColor(GREEN) ;
		  graphFillRectangle(x-.25, gy3, x+.25, gy4);
		  graphColor(BLACK);
		  if (t3) 
		    graphLine(x-.2, gy3, x+.2, gy3) ;
		  if (t4)
		    graphLine(x-.2, gy4, x+.2, gy4) ;
		  graphLine(x-.2, gy3, x-.2, gy4) ;
		  graphLine(x+.2, gy3, x+.2, gy4) ;
		} 
	      
	      /* a tick at best position */
	      cen = (y1 + y2) / 2 ;
	      if (controlCheckY(map, cen-best/2, 0, &gy1))
		graphLine (x-0.3, gy1, x+0.3, gy1);
	      if (controlCheckY(map, cen+best/2, 0, &gy1))
		graphLine (x-0.3, gy1, x+0.3, gy1);
	    }
	  
	  /* point the genes, must come last */
	  graphPointsize(.7) ;
	  if (controlCheckY(map, y1, 0.35, &gy1))
	    graphPoint(x, gy1) ;
	  if (controlCheckY(map, y2, 0.35, &gy2))
	    graphPoint(x, gy2) ;
	  
	  graphBoxEnd() ;
 	  if (seg->key == control->from)
	    control->fromBox = ibox;
	}
    }
  
  (*offset)++ ;
}

void twoPtAddKey(COLINSTANCE instance, KEY key)
/* Display the the piece of twoPt data whose key is key in this
   instance of the twopoint column. KEY zero is special, and means clear. */
{ TWOPRIV private;
  TWOPTSEG *seg;
  int i;
  TWOPTDATA data;
  
  private = instance->private;
  
  if (!key) /* clear */
    { arrayMax(private->segs) = 0;
      return;
    }

  if (!instance->displayed)
    instance->displayed = TRUE; /* pop it */

/* If we've already got it, don't do it again */
  for (i=0; i<arrayMax(private->segs); i++)
    { seg = arrp(private->segs, i, TWOPTSEG);
      if (seg->key == key)
	return;
    }
  /* now get a new one */
  
  if (!gMapGet2PtData(instance->map, key, instance->handle, &data))
    return;
 
  seg = arrayp(private->segs, arrayMax(private->segs), TWOPTSEG);
  seg->key = key;
  seg->data = data;
  if (data->y1 > data->y2)
    { seg->x = data->y2; seg->dx = data->y1 - data->y2; }
  else
    { seg->x = data->y1; seg->dx = data->y2 - data->y1; }
   
  return;

}  
  
static void twoPtDraw(COLINSTANCE instance, float *offset)
{ 
  if (instance->map->cambridgeOptions)
    twoPtRDDraw(instance, offset);
  else
    twoPtJTMDraw(instance, offset);
}

static BOOL twoPtCreate(COLINSTANCE instance, OBJ init)
{
  GeneticMap look = (GeneticMap)(instance->map->look);
  TWOPRIV private;

  private = (TWOPRIV)halloc(sizeof(struct TwoPrivStruct),
			    instance->handle);

  instance->private = private;
  private->segs = arrayHandleCreate(50, TWOPTSEG, instance->handle);
  
  instance->draw = twoPtDraw;
  instance->setSelectBox = twoPtSetSelect;
  instance->unSelectBox = twoPtUnselect; 
  instance->doColour = twoPtDoColour;
  instance->followBox = twoPtFollowBox;

  look->twoPtCurrentColumn = instance; 
  return TRUE;
}

static void twoPtDestroy(void *p)
{
  COLINSTANCE instance = (COLINSTANCE)p;
  GeneticMap look = (GeneticMap)(instance->map->look);
  
  if (look->twoPtCurrentColumn == instance)
    look->twoPtCurrentColumn = 0;
}

struct ProtoStruct gMapTwoPtColumn = {
  0,
  twoPtCreate,
  twoPtDestroy,
  "Two_point",
  0,
  TRUE,
  0,0,0,0,
  0
};


/******************/
/*     dbn column */
/******************/
#define NBIN 50
#define DBN_X(seg, z) ((((z)+1)*(seg)->dbnMax + \
			(NBIN-(z))*(seg)->dbnMin)/(NBIN+1))

typedef struct DbnSegStruct {
  KEY key;
  float dbn[NBIN], dbnMin, dbnMax;
  float x;
  COLINSTANCE instance;
} *DBNSEG;

typedef struct DbnPrivStruct {
  Array segs;
} *DBNPRIV;

void dbnAddKey (COLINSTANCE instance, KEY locus)
{ DBNPRIV private;
  DBNSEG seg;
  MAPCONTROL map;
  int i ;
  float x, dMax;
  KEY minLoc, maxLoc ;

  private = instance->private;
  map = instance->map;

  if (!locus) /* clear */
    { arrayMax(private->segs) = 0;
      return;
    }
 
  if (!getPos (map, locus, &x))
    return;
  
  gMapGetBadData() ;
  
  if (!instance->displayed)
    instance->displayed = TRUE; 
  
  seg = 0; /* look for existing one */
  for(i=0; i<arrayMax(private->segs); i++)
    if (arr(private->segs, i, DBNSEG) && 
	arr(private->segs, i, DBNSEG)->key == locus)
      seg = array(private->segs, i, DBNSEG);

  if (!seg) /* get new one */
    {
      seg = (DBNSEG)halloc(sizeof(struct DbnSegStruct),
			   instance->handle);
      array(private->segs, arrayMax(private->segs), DBNSEG) = seg;
    }

  seg->x = x;
  seg->instance = instance;
  seg->key = locus;

  if (!boundFind (map, locus, 
		  &seg->dbnMin, &seg->dbnMax, &minLoc, &maxLoc))
    { if (seg->dbnMin < -999999)
	seg->dbnMin = x - 10 ;
      if (seg->dbnMax > 999999)
	seg->dbnMax = x + 10 ;
    }

  dMax = -1E20 ; ;
  for (i = 0 ; i < NBIN ; ++i)
    { setTestPos (map, locus, DBN_X(seg, i)) ;
      seg->dbn[i] = logLikeLocus (map, locus) ;
      if (seg->dbn[i] > dMax)
	dMax = seg->dbn[i] ;
    }
  for (i = 0 ; i < NBIN ; ++i)
    seg->dbn[i] = exp(seg->dbn[i] - dMax) ;

  setTestPos (map, locus, x) ;
  return;
}

static void dbnErase(void *p)
{ DBNSEG seg = (DBNSEG)p;
  DBNPRIV private = seg->instance->private;
  Array segs  = private->segs;
  int i;

  for (i=0; i<arrayMax(segs); i++)
    { if (arr(segs, i, DBNSEG) == seg)
	arr(segs, i, DBNSEG) = 0;
    }

  controlDraw();
}
  

static void dbnDraw (COLINSTANCE instance, float *offset) 
{
  int i, j;
  DBNPRIV private = instance->private;
  DBNSEG seg;
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  float top, bottom, x;
  char buff[7];

  for (i = 0; i<arrayMax(private->segs); i++)
    { seg = arr(private->segs, i, DBNSEG);
      if (!seg) 
	continue; /* deleted one */
      
      /* see if the locus has moved, if so recalculate. */
      /* dbnAddKey will put new data in same seg, hence seg is still valid */
      if (getPos(map, seg->key, &x) && seg->x != x)
	dbnAddKey(instance, seg->key);

	
      bottom = MAP2GRAPH(map, seg->dbnMax);
      top = MAP2GRAPH(map, seg->dbnMin);
      if ((bottom < control->topMargin+2) ||
	  (top> control->graphHeight))
	continue;
	  
      buff[0] = 0;
      strncat(buff, name(seg->key), 6);
      if (strlen(buff) != strlen(name(seg->key)))
	{ buff[0] = 0;
	  strncat(buff, name(seg->key), 4);
	  strcat(buff, "..");
	}
      graphColouredButton(buff, dbnErase, seg, BLACK, WHITE, 
			  *offset+0.5, control->topMargin+0.5); 

      graphColor (RED) ;
        if (top > control->topMargin+2)
	  graphLine (*offset, top, *offset+7, top) ;
      if (bottom < control->graphHeight)
	graphLine (*offset, bottom, *offset+7, bottom) ;

      graphColor (BLUE) ;
      for (j = 0 ; j < NBIN ; ++j)
	{ x = MAP2GRAPH(map, DBN_X(seg, j)) ;
	  if (x > control->topMargin+2 && x < control->graphHeight)
	    graphLine (*offset, x, *offset + 7*seg->dbn[j], x) ;
	}
      graphColor (BLACK) ;
      *offset += 8 ;
    }
}


static BOOL dbnCreate(COLINSTANCE instance, OBJ init)
{
  GeneticMap look = (GeneticMap)(instance->map->look);
  DBNPRIV private;

  private = (DBNPRIV)halloc(sizeof(struct DbnPrivStruct),
			    instance->handle);

  instance->private = private;
  private->segs = arrayHandleCreate(50, DBNSEG, instance->handle);
  instance->draw =  dbnDraw;


  look->dbnCurrentColumn = instance;
  
  return TRUE;
}

static void dbnDestroy(void *p)
{
  COLINSTANCE instance = (COLINSTANCE)p;
  GeneticMap look = (GeneticMap)(instance->map->look);
  
  if (look->dbnCurrentColumn == instance)
    look->dbnCurrentColumn = 0;
}


struct ProtoStruct gMapLikelihoodColumn = {
  0,
  dbnCreate,
  dbnDestroy,
  "Likelihood",
  0,
  TRUE,
  0,0,0,0,
  0
};

