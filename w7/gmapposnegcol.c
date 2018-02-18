/*  File: gmapposnegcol.c
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: positive/negative data column for genetic map
 * Exported functions:
 *     gMapPosNegColumn
 * HISTORY:
 * Last edited: Dec 21 11:59 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: gmapposnegcol.c,v 1.2 2014/11/30 03:20:46 mieg Exp $ */


#include "gmap.h"

/************************************************************/

typedef struct PosNegPrivStruct {
  KEY pointKey, intervalKey;
  Array segs;
  int posColour, negColour, contColour;
  BOOL showAll;
  KEYSET keyset;
  char *query;
  float spacing;
} *POSNEGPRIV;

typedef struct PerSpotDR {
  int parentBox;
  KEY intervalKey, pointKey, dataKey, tag;
} PSDR;

typedef struct DataStruct {
  KEY point;
  KEY tag;
  KEY data;
} DATA;

/************************************************************/

static BOOL setSelect(COLINSTANCE instance, int box, double x, double y)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  POSNEGPRIV private = instance->private;
  GeneticMap look = (GeneticMap)(map->look);
  Array segs = private->segs;
  int i = assInt (controlBoxRegd(instance, box)) ;
  KEY key;
  PSDR *seg = arrp(segs, i, PSDR);

  gMapUnselect(map);

  private->pointKey = seg->pointKey;
  private->intervalKey = seg->intervalKey;
  
  if (seg->pointKey) /* selecting a point */
    key = seg->pointKey;
  else
    key = seg->intervalKey;
  
  look->selectKey = control->activeKey = key;
  look->friendsInfoValid = TRUE;
  look->neighboursInfoValid = TRUE;
  if (seg->pointKey)
    { float x;
      getPos(map, key, &x);
      look->x1 = look->x2 = x;
    }
  else
    { look->coOrdsInvalid = TRUE;
    }
  (void)gMapNeighbours(map, &look->neighbours, 
		       seg->dataKey ? seg->dataKey : key);
  (void)gMapPositive(map, &look->positive, key);
  (void)gMapNegative(map, &look->negative, key);
  
  return FALSE;
} /* setSelect */


static BOOL unSelect(COLINSTANCE instance, int box)
{
  GeneticMap look = (GeneticMap)(instance->map->look);
 
  *look->messageText = 0;
  if (look->messageBox)
    graphBoxDraw(look->messageBox, -1, -1);

  return FALSE;
} /* unSelect */


static void followBox(COLINSTANCE instance, int box, double x, double y)
{
  MAPCONTROL map = instance->map;
  POSNEGPRIV private = instance->private;
  int n = assInt (controlBoxRegd(instance, box)) ;
  PSDR *seg = arrp(private->segs, n, PSDR);

  if (seg->dataKey)
    display(seg->dataKey, map->key, TREE);
  else if (seg->pointKey)
    display(seg->pointKey, map->key, TREE);
  else
    display(seg->intervalKey, map->key, TREE);

  return;
} /* followBox */


static void doColour(COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);
  POSNEGPRIV private = instance->private;
  Array segs = private->segs;
  int i = assInt (controlBoxRegd(instance, box)) ;
  PSDR *seg = arrp(segs, i, PSDR);

  if (seg->pointKey == private->pointKey &&
      seg->intervalKey == private->intervalKey &&
      instance == control->activeInstance)
    { graphBoxDraw(box, -1, CYAN);
      *look->messageText = 0;
      if (seg->pointKey)
	{ gMapAddToHeader(map, name(seg->pointKey));
	  gMapAddToHeader(map, " is ");
	  gMapAddToHeader(map, name(seg->tag));
	  gMapAddToHeader(map, " in ");
	  gMapAddToHeader(map, name(seg->intervalKey));
	  if (seg->dataKey)
	    { gMapAddToHeader(map, " (according to ");
	      gMapAddToHeader(map, name(seg->dataKey));
	      gMapAddToHeader(map, ")");
	    }
	}
      else
	gMapAddToHeader(map, name(seg->intervalKey));

      if (look->messageBox)
	graphBoxDraw(look->messageBox, -1, -1);
      control->activeBox = box;
    }
  else
    { if (seg->pointKey)
        { if (gMapIsNeighbour(map, seg->pointKey))
	    graphBoxDraw(box, -1, PALECYAN); /* dot box */
	  else if (seg->dataKey && gMapIsNeighbour(map, seg->dataKey))
	    graphBoxDraw(box, -1, PALECYAN);
	  else
	    graphBoxDraw(box, -1, TRANSPARENT); /* dot box */
	}
      else
	{ if (gMapIsNeighbour(map, seg->intervalKey))
	    graphBoxDraw(box, -1, PALECYAN); /* interval box*/
	  else
	    graphBoxDraw(box, -1, WHITE); /* intervalbox */
	}    
    }

  if (seg->parentBox)
    graphBoxDraw(seg->parentBox, -1, -1);

  return;
} /* doColour */
	 
 
static int dataOrder(const void *a, const void *b)
{
  KEY k1 = ((const DATA *)a)->point;
  KEY k2 = ((const DATA *)b)->point;

  return k1-k2;
} /* dataOrder */


static void getData(KEY key, Array *a, KEY typeTag)
{
  OBJ obj, comp;
  int i;
  KEY pnData;
  Array loci = arrayCreate(50, BSunit);

  *a = arrayReCreate(*a, 100, DATA);
  
  if ((obj = bsCreate(key)))
    { if (bsFindTag(obj, typeTag) && bsFlatten(obj, 2, loci))
	for (i = 0; i<arrayMax(loci); i+= 2)
	  { KEY k = arr(loci, i+1, BSunit).k;
	    int n = arrayMax(*a);
	    if (k)
	      { array(*a, n, DATA).point = k;
		array(*a, n, DATA).tag = arr(loci, i, BSunit).k;
		array(*a, n, DATA).data = 0;
	      }
	  }

      if (bsGetKey(obj, str2tag("Pos_neg_data"), &pnData))
	do 
	  { KEY mo1, t1, mo2, t2;
	    int n = arrayMax(*a);
	    comp = bsCreate(pnData);
	    if (!bsFindTag(comp, typeTag))
	      { bsDestroy(comp);
		continue;
	      }
	
	    if (bsFindTag(comp, str2tag("Item_1")) &&
		bsGetKeyTags(comp, _bsRight, &t1) &&
		bsGetKey(comp, _bsRight, &mo1) &&
		bsFindTag(comp, str2tag("Item_2")) &&
		bsGetKeyTags(comp, _bsRight, &t2) &&
		bsGetKey(comp, _bsRight, &mo2))
	      {
		if (mo1 == key)
		  { array(*a, n, DATA).point = mo2;
		    array(*a, n, DATA).tag = t2;
		    array(*a, n, DATA).data = pnData;
		  }
		if (mo2 == key)
		  { array(*a, n, DATA).point = mo1;
		    array(*a, n, DATA).tag = t1;
		    array(*a, n, DATA).data = pnData;
		  }
	      }
	    bsDestroy(comp);
	  }
	while (bsGetKey(obj, _bsDown, &pnData));
      
    }  
 
  arraySort(*a, dataOrder); 
  bsDestroy(obj);
  arrayDestroy(loci);

  return;
} /* getData */
  
  

static void posNegDraw (COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  POSNEGPRIV private = instance->private;
  float x  = *offset + 0.5;
  int i ;
  Array segs;
  Array pos = 0, neg = 0;
  

  if (!private->keyset)
    return;

  graphTextHeight(0.8);
  segs = private->segs = arrayReCreate(private->segs, 400, PSDR);

  for (i=0; i < keySetMax(private->keyset); i++)
    {
      enum { NINTS=200 };
      float centre[NINTS], hwidth[NINTS];
      KEY key = keySet(private->keyset, i);
      float endOfName = control->topMargin + 0.65 * strlen(name(key));
      float mt = GRAPH2MAP(map, endOfName);
      float mb = GRAPH2MAP(map, control->graphHeight);
      int j, nints, flag, ibox;

      if (mt > mb)
	{ float tmp = mt; mt = mb; mb  = tmp; }
      nints = NINTS;
    
      /* look for intervals on the current display */
      flag = gMapGetMapObject(key, map->key, 0, 50, centre, hwidth, 0, &nints);
      getData(key, &pos, str2tag("Positive")); /* get mapping info */
      getData(key, &neg, str2tag("Negative"));
      
      if (private->showAll)
	goto found; /* config option */
      
      /* see if we have on interval on the screen */
      if (flag & FLAG_ANY_INTERVAL)
	for (j = 0; j < nints; j++)
	  if (centre[j] - hwidth[j] < mb &&
	      centre[j] + hwidth[j] > mt)
	    goto found;
      
      /* if no good above, look for map points on the screen */
      for (j = 0; j < arrayMax(pos); j++)
	{ float p;
	  if (getPos(map, arr(pos, j, DATA).point, &p))
	    if (p > mt && p < mb)
	      goto found;
	}
      for (j = 0; j < arrayMax(neg); j++)
	{ float p;
	  if (getPos(map, arr(neg, j, DATA).point, &p))
	    if (p > mt && p < mb)
	      goto found;
	}
      
      /* no mapped points or intervals on screen, try next */
      continue;

    found:
      /* print the name and gray line */
      { int intSeg = arrayMax(segs);
	char a[2];
	char *nam = name(key);
	int i, n = strlen(nam);
	float y = control->topMargin;

	ibox = graphBoxStart();
	controlRegBox(instance, ibox, assVoid (intSeg));

	array(segs, intSeg, PSDR).parentBox = 0; 
	array(segs, intSeg, PSDR).intervalKey = key;
	array(segs, intSeg, PSDR).pointKey = 0;
	array(segs, intSeg, PSDR).dataKey = 0;

	a[1] = 0;
	graphColor(BLACK);
	for (i = 0; i < n; i++)
	  { a[0] = *nam++;
	    graphText(a, x, y);
	    y += 0.65;
	  }
	graphColor(GRAY);
	graphLine(x+0.5, endOfName, x+0.5, control->graphHeight-0.1);
      }
      
      /* now do any black lines for intervals */
      for(j=0; j<nints; j++)
	{ float g1, g2;
	  controlCheckY(map, centre[j] - hwidth[j], 0, &g1);
	  controlCheckY(map, centre[j] + hwidth[j], 0, &g2);
	  if (g1 != g2)
	    { if (g1 > g2)
		{ float tmp = g1; g1 = g2; g2 = tmp; }
	      if (g1 < endOfName + 0.2)
		g1 = endOfName + 0.2;
	      if (g2 > endOfName + 0.2)
		{ graphColor(BLACK);
		  graphFillRectangle(x+0.4, g1, x+0.6, g2);
		}
	    }
	}
      
      /* finally, do mapped points */
      for (j = 0; j < arrayMax(pos); j++)
	{ float coord;
	  if (getPos(map, arr(pos, j, DATA).point, &coord))
	    { float gcoord = MAP2GRAPH(map, coord);
	      int n = arrayMax(segs);
	      if ((gcoord > endOfName + 0.5) &&
		  gcoord < control->graphHeight - 0.5)
		{ int box = graphBoxStart();
		  controlRegBox(instance, box, assVoid (n));
		  array(segs, n, PSDR).parentBox = 0;
		  if (arrayFind(neg, &(arr(pos, j, DATA).point), 0, dataOrder))
		    graphColor(private->contColour);
		  else  
		    graphColor(private->posColour);
		  graphFillArc(x+0.5, gcoord, 0.4, 0, 360);
		  array(segs, n, PSDR).pointKey = array(pos, j, DATA).point;
		  array(segs, n, PSDR).intervalKey = key;
		  array(segs, n, PSDR).dataKey = array(pos, j, DATA).data;
		  array(segs, n, PSDR).tag = array(pos, j, DATA).tag;
		  graphBoxEnd();
		  remarkRegister(map, array(pos, j, DATA).point, coord);
		  if (array(pos, j, DATA).data)
		    remarkRegister(map, array(pos, j, DATA).data, coord);
		}
	    }
	}
      graphColor(private->negColour);
      for (j = 0; j < arrayMax(neg); j++)
	{ float coord;
	  if (!arrayFind(pos, &(arr(neg, j, DATA).point), 0, dataOrder) &&
	      getPos(map, arr(neg, j, DATA).point, &coord))
	    { float gcoord = MAP2GRAPH(map, coord);
	      int n = arrayMax(segs);
	      int box = graphBoxStart();
	      controlRegBox(instance, box, assVoid(n));
	      if ((gcoord > endOfName + 0.5) &&
		  gcoord < control->graphHeight - 0.5)
		graphFillArc(x+0.5, gcoord, 0.4, 0, 360);
	      array(segs, n, PSDR).parentBox = 0;
	      array(segs, n, PSDR).pointKey = arr(neg, j, DATA).point;
	      array(segs, n, PSDR).intervalKey = key;
	      array(segs, n, PSDR).dataKey = arr(neg, j, DATA).data;
	      array(segs, n, PSDR).tag = arr(neg, j, DATA).tag;
	      graphBoxEnd();
	      remarkRegister(map, array(neg, j, DATA).point, coord);
		  if (array(neg, j, DATA).data)
		    remarkRegister(map, array(neg, j, DATA).data, coord);
	    }
	}
      
      graphBoxEnd();
      array(segs, arrayMax(segs)-1, PSDR).parentBox = ibox;
      x += private->spacing;
      /* then loop to the next */
    }
  
  arrayDestroy(pos);
  arrayDestroy(neg);
  *offset = x + 0.5;
  graphTextHeight(0.0);

  return;
} /* posNegDraw */

static void posNegPrivDestroy(void *p)
{ POSNEGPRIV private = (POSNEGPRIV)p;
  
  if (private->keyset) 
    keySetDestroy(private->keyset);
}

static void posNegSave(COLINSTANCE instance, OBJ obj)
{
  POSNEGPRIV private = instance->private;

  if (private->query)
    bsAddData(obj, str2tag("RH_query"), _Text, private->query);
  
  if (private->showAll)
    bsAddTag(obj, str2tag("RH_show_all"));

  bsAddData(obj, str2tag("RH_spacing"), _Float, &private->spacing);
  
  bsAddTag(obj, str2tag("RH_positive")); 
  controlSetColour(obj, private->posColour);
  
  bsAddTag(obj, str2tag("RH_negative")); 
  controlSetColour(obj, private->negColour);

  bsAddTag(obj, str2tag("RH_contradictory"));
  controlSetColour(obj, private->contColour);
  
}
struct configLocals {
  char query[280];
  float spacing;
  int colour[3];
  BOOL showAll;
};

static BOOL posNegConfigure(COLINSTANCE instance)
{ POSNEGPRIV private = instance->private;
  struct  configLocals *cf = (struct configLocals *) messalloc(sizeof(struct configLocals));

  if(controlCreateConfig(instance,cf,"Configure pos/neg column",0.8,0.3)){
    
  /*initialise data*/
    if(private->query)
      strcpy(cf->query,private->query);
    cf->showAll = private->showAll;
    cf->colour[0] = private->posColour;
    cf->colour[1] = private->negColour;
    cf->colour[2] = private->contColour;
    cf->spacing = private->spacing;

    graphTextEditor("Query for display:",cf->query,280,4,2,0);
    graphFloatEditor("Interval spacing:",&cf->spacing,4.0,4.0,0);

    graphToggleEditor("Show dataless intervals",&cf->showAll,4.0,11.0);

    graphColourEditor("Positive data"," ",&cf->colour[0],36.0,10.0);
    graphColourEditor("Negative data"," ",&cf->colour[1],36.0,11.1);
    graphColourEditor("Contradictory data"," ",&cf->colour[2],36.0,12.2);
    graphRedraw();
  }
  return FALSE;
}

static void gmapPosNegFinal(COLINSTANCE instance, void *locals, BOOL ok)
{ POSNEGPRIV private = instance->private;
  struct configLocals *cf = locals;
  
  if(ok){
    private->showAll = cf->showAll;
    private->posColour = cf->colour[0];
    private->negColour = cf->colour[1];
    private->contColour = cf->colour[2];
    private->spacing = cf->spacing;

    if (private->keyset)
      keySetDestroy(private->keyset);
    if (private->query)
      messfree(private->query);
    
    if (strlen(cf->query) != 0)
      { private->query = handleAlloc(0, 
				     instance->handle, 
				     1+strlen(cf->query));
	  strcpy(private->query, cf->query);
	private->keyset = queryKey(instance->map->key, private->query);
      }
    else
      { private->keyset = 0;
	private->query = 0;
      }
  }  
  else
    messfree(cf);  
}

static BOOL posNegCreate(COLINSTANCE instance, OBJ init)
{ POSNEGPRIV private;
  char *s1;
  float f;

  instance->draw = posNegDraw;
  instance->setSelectBox = setSelect;
  instance->unSelectBox = unSelect;
  instance->doColour = doColour;
  instance->followBox = followBox;
  instance->configure = posNegConfigure;
  instance->save = posNegSave;
  instance->configFinal = gmapPosNegFinal;
  private = (POSNEGPRIV) handleAlloc(posNegPrivDestroy, instance->handle,
				     sizeof(struct PosNegPrivStruct));
  instance->private = private;

  private->posColour = DARKGREEN;
  private->negColour = BLUE;
  private->contColour = RED;
  private->keyset = 0;
  private->query = 0;
  private->showAll = FALSE;
  private->spacing = 1.0;

  if (init)
    { if (bsGetData(init, str2tag("RH_query"), _Text, &s1))
	{ private->query = handleAlloc(0, instance->handle, 1+strlen(s1));
	  strcpy(private->query, s1);
	  private->keyset = queryKey(instance->map->key, private->query);
	}

      if (bsGetData(init, str2tag("RH_spacing"), _Float, &f))
	private->spacing = f;

       if (bsFindTag(init, str2tag("RH_positive")))
	private->posColour = controlGetColour(init);
      if (bsFindTag(init, str2tag("RH_negative")))
	private->negColour = controlGetColour(init);
      if (bsFindTag(init, str2tag("RH_contradictory")))
	private->contColour = controlGetColour(init);

      if (bsFindTag(init, str2tag("RH_show_all")))
	private->showAll = TRUE;
    }

  return TRUE;
}



struct ProtoStruct gMapPosNegColumn = {
  0,
  posNegCreate,
  0,
  "RH_data",
  0,
  FALSE,
  0,
  0,
  "The RH_data column shows a set of gray vertical lines, one for each "
    "object in the displayed keyset. These are intended to be interval-type "
      "entities, but need not have co-ordinates as such on the map. If they "
	"do have co-ordinates, their extent is shown by black lines drawn "
	  "atop the gray ones. \n\n"
	    "For any point-type objects which have "
	      "co-ordinates in the map and are positive or negative for "
		"the displayed intervals, coloured spots are drawn on the "
		  "interval lines at the mapped position of the point. "
		    "If a point is both postive and negative for an interval, "
		      "the \"contradictory colour\" is used to draw it.\n\n"
  "Unless \"show dataless intervals\" is selected, the drawing of intervals "
    "is suppressed when they have neither co-ordinates or mapped points "
      "on the screen.\n\n"
	"Note that, for this column, the query to display is mandatory; "
	  "without it, nothing is shown.\n"
};







 
 
 
