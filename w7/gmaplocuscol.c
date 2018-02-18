/*  File: gmaplocuscol.c
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Locus display column for the genetic map
 * Exported functions:
 *     gMapPointColumn  
 * HISTORY:
 * Last edited: Dec 21 14:29 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: gmaplocuscol.c,v 1.3 2017/03/18 15:30:56 mieg Exp $ */

#include "gmap.h"
#include "bump.h"

/************************************************************/

BOOL gMapIsNegMarkers(MAPCONTROL map, KEY key);

typedef struct tags {
  char *text;
} TAG;

typedef struct clonedQuerys{
  char *query;
  KEYSET keyset;
  int colour;
} CLONEDQUERY;

typedef struct GenePrivStruct {
  char *query;
  Array clonedQuerys;
  KEYSET keyset;
  int maxWidth;
  float errorScale;
  int pne, pe, nne, ne;
  BOOL showOrdered;
  BOOL showMarginal;
  Array symbolTags;
  Array segs;
  /* may show same key more than once, use these to distinguish */
  float selx, seldx; 
} *GENESPRIV;

/* called when data button pressed, something drawn by us is highlighted.
   add the map data about it to the display. To do this we use the
   rountine twpPtAddKey and multiPtAddKey exported by the data columns.
*/

static void genesShowData(MAPCONTROL map)
{  
  GeneticMap look = (GeneticMap)map->look;
  int i;
  KEYSET k = 0;
  
  if (gMap2Pt(map, &k, look->selectKey))
    { if (look->twoPtCurrentColumn)
	for (i = 0 ; i < keySetMax(k) ; ++i)
	  twoPtAddKey(look->twoPtCurrentColumn, keySet(k, i));
      else
	messout("Cannot display Two Point Data, there is no suitable column "
		"in the current map display.");
    }
  
  if (gMapMultiPt(map, &k, look->selectKey))
    { if (look->multiPtCurrentColumn)
	for (i = 0 ; i < keySetMax(k) ; ++i)
	  multiPtAddKey(look->multiPtCurrentColumn, keySet(k, i));
      else
	messout("Cannot display Multi Point Data, there is no suitable column "
		"in the current map display.");
    }

  keySetDestroy(k);
  controlDrawControl(map->control) ;
}

static void genesShowDataFromButton(void *p)
{
  COLPROTO buttonProto = p;
  MAPCONTROL map = currentMapControl();

  /* pressing this button acts on the currently selected thing, it is 
     therefore imperative to make sure that this belongs to us, so we
     can trust the stored select information. To do this we make sure that
     there is an activeInstance (this is always valid), and that it's
     prototype is the same as the one which placed the button. */
  if (map->control->activeInstance &&
      map->control->activeInstance->proto == buttonProto)
    genesShowData(map);
}


static void genesShowDataFromMenu(int box)
{
  /* As above, but select the item first */
  controlSelectBox(box);
  /* current is OK, we've selected our box */
  genesShowData(currentMapControl()); 
}

  
static void genesShowDbnFromMenu(int box)
{
  MAPCONTROL map;
  GeneticMap look;
  
  controlSelectBox(box);
  map = currentMapControl();
  look = (GeneticMap)map->look;
  if (look->dbnCurrentColumn)
    dbnAddKey(look->dbnCurrentColumn, look->selectKey);
  else
    messout("Cannot display Liklihood Distribution, there is no suitable "
	    "column in the current map display.");
  controlDrawControl(map->control) ;
}


static void genesPosToCursor (int box)
{
  COLCONTROL control = currentColControl("genesPosToCursor");
  COLINSTANCE instance = arr(control->boxIndex, box, COLINSTANCE);
  MAPCONTROL map = instance->map;
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *);
  GENESPRIV private = instance->private;  
  float a; 
  GMAPSEG save;
  int i;

  save = *seg;
  
  a = mapControlCursorPos(map) ;
  mapControlCursorSet (map, save.x) ;

  save.flag |= FLAG_MOVED ;
  setTestPos(map, save.key, a);
  save.x = a ;

  arrayRemove(private->segs, seg, gMapOrder); /* to maintain order */
  arrayInsert(private->segs, &save, gMapOrder); /* to maintain order */
  arrayFind(private->segs, &save, &i, gMapOrder);
  /* must update this, in order for controlSelectBox to set correct values */
  arr(control->boxIndex2, box, void *) = arrp(private->segs, i, GMAPSEG);

  controlSelectBox(box);
  controlDraw();

  return;
} /* genesPosToCursor */


static void genesCursorToPos (int box)
{
  COLCONTROL control = currentColControl("genesCursorToPos");
  MAPCONTROL map = arr(control->boxIndex, box, COLINSTANCE)->map;
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *);
  
  mapControlCursorSet (map, seg->x) ;

  return;
} /* genesCursorToPos */


static void genesHighlightSeg(int box)
{
  COLCONTROL control = currentColControl("genesHighlightSeg");
  MAPCONTROL map = arr(control->boxIndex, box, COLINSTANCE)->map;
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *);

  gMapHighlightKey(map, seg->key);
  
  return;
} /* genesHighlightSeg */



static void genesDoColour(COLINSTANCE instance, int box)
/* NB this used by the rearrangements column also */
{ 
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *); 
  GeneticMap look = (GeneticMap)instance->map->look;
  GENESPRIV private = instance->private;
  BOOL okay;
  int i;
  
  if (seg->key == look->selectKey && instance == control->activeInstance &&
      seg->x == private->selx && seg->dx == private->seldx)
    { graphBoxDraw(box, BLACK, CYAN);
      strncpy(look->messageText, name(seg->key), 100);
      strcat(look->messageText, messprintf(" %.2f", seg->x));
      strcat(look->messageText, messprintf(" [%.2f]", seg->dx));
  
      if (look->messageBox) 
	graphBoxDraw(look->messageBox, -1, -1);
      control->activeBox = box;
    }
  else
    { float x1 = private->showMarginal ? seg->x - seg->dx : seg->x;
      float x2 = private->showMarginal ? seg->x + seg->dx : seg->x;
      int flags = gMapOverlap(map, seg->key, x1, x2);
      if (keySetFind(look->highlight, seg->key, 0))
	graphBoxDraw(box, BLACK, MAGENTA);
      else if (flags & FLAG_STRESSED)
	graphBoxDraw(box, graphContrast(private->pe), private->pe);
      else if (flags & FLAG_ANTI_STRESSED)
	graphBoxDraw(box, graphContrast(private->ne), private->ne);
      else if (flags & FLAG_ANTI_RELATED)
	graphBoxDraw(box, graphContrast(private->nne), private->nne);
      else if (flags & FLAG_RELATED)
	 graphBoxDraw(box, graphContrast(private->pne), private->pne);
      else if (gMapIsNegMarkers(map,seg->key))
	graphBoxDraw(box, BLACK, RED);
      else if (gMapIsNeighbour(map, seg->key))
	graphBoxDraw(box, BLACK, PALECYAN);
      else
	{ okay = FALSE;
	  for (i = 0 ; i < arrayMax(private->clonedQuerys) ; i++)
	    { CLONEDQUERY *q = arrp(private->clonedQuerys,i,CLONEDQUERY) ;
	      if (q->keyset && keySetFind(q->keyset, seg->key, 0)
		  /* && !instance->map->control->hideHeader apb fix*/)
		{ graphBoxDraw(box, BLACK, q->colour);
		  okay = TRUE;
		  break;
		}
	    }
	  if(!okay)
	    graphBoxDraw(box, BLACK, WHITE);
	}
    }
} /* genesDoColour */

static BOOL genesSetSelect(COLINSTANCE instance, int box, double x, double y)
{ MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GMAPSEG *seg = (GMAPSEG *) arr(control->boxIndex2, box, void *); 
  GeneticMap look = (GeneticMap)map->look;
  GENESPRIV private = instance->private;

  look->selectKey = seg->key;
  private->selx = seg->x;
  private->seldx = seg->dx;
  gMapUnselect(map);
  look->neighboursInfoValid = TRUE;
  look->friendsInfoValid = TRUE;
  look->dataInfoValid = TRUE; /* used by data columns */
  look->x1 = private->showMarginal ? seg->x - seg->dx : seg->x;
  look->x2 = private->showMarginal ? seg->x + seg->dx : seg->x; 
  (void)gMapNeighbours(map, &look->neighbours, seg->key);
  (void)gMapPositive(map, &look->positive, seg->key);
  (void)gMapNegative(map, &look->negative, seg->key);
  (void)gMapMultiPt(map, &look->multi, seg->key);  
  (void)gMap2Pt(map, &look->twoPt, seg->key); 
  control->activeKey = seg->key; /* inter-map communication */  

  return FALSE; /* never need to redraw */
}

static BOOL genesUnselect(COLINSTANCE instance, int box)
{ 
  GeneticMap look = (GeneticMap)instance->map->look;

  *look->messageText = 0;

  if (look->messageBox)
    graphBoxDraw(look->messageBox, -1, -1);
  
  return FALSE;
}

static void genesFollowBox(COLINSTANCE instance, int box, double x, double y)
{
  GMAPSEG *seg = (GMAPSEG *) arr(instance->map->control->boxIndex2, 
				 box, void *); 
  
  display(seg->key, instance->map->key, TREE);
}

static char * gMapKeyOnTags(KEY key, Array tags, COLINSTANCE instance)
{ OBJ obj;
  TAG *tag;
  int i;
  int ok;
  KEY key2;
  BSMARK mark=0;
  char text[180];
  char *showName=0;

  if(!arrayMax(tags))
    { showName = strnew(name(key),instance->handle);
      return showName;
    }
  if(!(obj = bsCreate(key)))
    { showName = strnew(name(key),instance->handle);
      return showName;
    }

  mark = bsMark(obj,mark);

  for(i=0;i<arrayMax(tags);i++)
    { tag = arrp(tags,i,TAG);
      bsGoto(obj,mark);
      if(bsGetKey(obj,str2tag(tag->text),&key2))
	{  showName = strnew(name(key2),instance->handle);
	   ok =TRUE;
	 }
      else
	{ bsGoto(obj,mark);
	  ok = bsFindTag(obj,str2tag(tag->text));
	  if(ok){
	    /* what type follows */
	    if(bsType(obj,_bsRight)==_Text)
	      { ok = bsGetData(obj, _bsRight, _Text, &text);
		showName = strnew(text,instance->handle);
	      }
	  }
	}
      if (ok)
	{ bsMarkFree(mark);
	  bsDestroy(obj);
	  return showName;
	}
    }
  bsMarkFree(mark);
  bsDestroy(obj);
  showName = strnew(name(key),instance->handle);
  return showName;

}

static void genesDraw (COLINSTANCE instance,
			      float *offset)
{
  static MENUOPT geneMenu[] = {
    {(VoidRoutine)genesCursorToPos, "Set cursor"},
    {(VoidRoutine)genesHighlightSeg, "Highlight"},
    {(VoidRoutine)genesPosToCursor, "Move to Cursor"},
    {(VoidRoutine)genesShowDataFromMenu, "Show data"},
    {(VoidRoutine)genesShowDbnFromMenu, "Likelihood dbn"},
    {0, 0}
  } ;
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)map->look;
  BUMP  bump ;
  float y ;
  int   i, ibox, n, x ;
  GMAPSEG   *seg ;
  GENESPRIV private = instance->private;
  float orderedWidth;
  int firstSeg, lastSeg, countDir;
  int nameLength;
  char *showName = 0;
  char buf[280];

  if (map->mag > 0)
    { firstSeg = 0;
      lastSeg = arrayMax(private->segs);
      countDir = 1;
    }
  else
    { firstSeg = arrayMax(private->segs)-1;
      lastSeg = -1;
      countDir = -1;
    }


  *offset += 1 ;		/* bump bug */
  if (!private->showOrdered)
    orderedWidth = 0.0;
  else
    {
       if (private->maxWidth)
	 x = private->maxWidth/2;
       else
	 x = (control->graphWidth - *offset)/2 ;
      
       if (x < 12) x = 12 ;
 
      bump = bumpCreate (x, 0) ;
      
      for (i = firstSeg ; i != lastSeg ; i += countDir)
	{ seg = arrp(private->segs,i,GMAPSEG) ;
	  y = MAP2GRAPH(map, seg->x) ;
	  
	  if (!(seg->flag & FLAG_ANY_LOCUS) || 
	      !(seg->flag & FLAG_WELL_ORDERED))
	    continue;
	  if (y < (control->topMargin+0.5) || y > control->graphHeight-0.5)
	    continue;
	  if (private->keyset && !keySetFind(private->keyset, seg->key, 0))
	    continue;
	  if (keySetFind(look->hidden, seg->key, 0))
	    continue ;
	  
	  ibox = graphBoxStart();
	  array(control->boxIndex, ibox, COLINSTANCE) = instance;
	  array(control->boxIndex2, ibox, void *) = seg ;
	  graphBoxInfo(ibox, seg->key, 0) ; /* WWW */
	  x = 0;
	  showName = gMapKeyOnTags(seg->key, private->symbolTags, instance);
	  nameLength = strlen(showName);

	  if ((n = bumpText (bump, showName, &x, &y, 1, TRUE)) &&
	      y < control->graphHeight-0.5)
	    { if(nameLength < n)
	        n = nameLength;
	      strncpy(buf, showName, n);
	      buf[n] = '\0';
	      graphText (buf, *offset+x, y-0.5) ;
	      remarkRegister(map, seg->key, seg->x);
	      /*if (control->hideHeader && 
		  private->clonedKeyset &&
		  keySetFind(private->clonedKeyset, seg->key, 0))*/
		/* use underlining for hard copy */
/*		{ graphLine (*offset+x, y+0.3, 
			     *offset+x+nameLength*0.6, y+0.3) ;
			     }*/
	    }
	  if (seg->key == control->from)
	    control->fromBox = ibox;
	  
	  graphBoxEnd() ;
	  if (map->submenus)
	    graphBoxMenu (ibox, geneMenu) ;
	}
  
      orderedWidth = bumpMax (bump) ;
      bumpDestroy (bump) ;
    }
  
  
  if (private->maxWidth)
    x = private->maxWidth - orderedWidth;
  else
    x = control->graphWidth - (*offset + orderedWidth); 
  /* width to fill the screen */
  
  if (x >= 6 ) 
    {  
      bump = bumpCreate (x, 0) ;
  
      for (i = firstSeg; i != lastSeg; i += countDir)
	{ seg = arrp(private->segs,i,GMAPSEG) ;
	  y = MAP2GRAPH(map, seg->x) ;
	  
	  if (!(seg->flag & FLAG_ANY_LOCUS))
	    continue;
	  if (y < control->topMargin+0.5 || y > control->graphHeight-0.5)
	    continue;
	  if (private->showOrdered && (seg->flag & FLAG_WELL_ORDERED))
	    continue;
	  if (private->keyset && !keySetFind(private->keyset, seg->key, 0))
	    continue;
	  if (keySetFind(look->hidden, seg->key, 0))
	    continue ;
	  
	  ibox = graphBoxStart();
	  array(control->boxIndex, ibox, COLINSTANCE) = instance;
	  array(control->boxIndex2, ibox, void *) = seg ;
	  graphBoxInfo(ibox, seg->key, 0) ; /* WWW */
	  x = private->errorScale * seg->dx;
	  showName = gMapKeyOnTags(seg->key, private->symbolTags, instance);
	  nameLength = strlen(showName);
	  
	  if ((n = bumpText (bump, showName, &x, &y, 1, TRUE)) &&
	      y < control->graphHeight-0.5)
	    { if(nameLength < n)
	        n = nameLength;
	      strncpy(buf, showName, n);
	      buf[n] = '\0';
	      graphText (buf, *offset+orderedWidth+x, y-0.5) ;
	      remarkRegister(map, seg->key, seg->x);
	      /*if (control->hideHeader && 
		  private->clonedKeyset &&
		  keySetFind(private->clonedKeyset, seg->key, 0))*/
		/* use underlining for hard copy */
		/*{ graphLine (*offset+x, y+0.3, 
			     *offset+x+nameLength*0.6, y+0.3) ;
		}*/
	    }	  
	  if (seg->key == control->from)
	    control->fromBox = ibox;
	  
	  graphBoxEnd() ;
	  if (map->submenus)
	    graphBoxMenu (ibox, geneMenu) ;
	}
      
      orderedWidth += (bumpMax (bump));
      bumpDestroy (bump) ;
    }

  *offset += orderedWidth;
  if(showName)
    messfree(showName);
}

static void genesSave(COLINSTANCE instance, OBJ obj)
{
  int i;
  GENESPRIV private = instance->private;

  if (private->query)
    bsAddData(obj, str2tag("Point_query"), _Text, private->query);

  for (i = 0 ; i < arrayMax(private->clonedQuerys) ; i++)
    { CLONEDQUERY *q = arrp(private->clonedQuerys,i,CLONEDQUERY) ;
      bsAddData(obj, str2tag("Point_colour"), _Text, q->query) ;
      controlSetColour (obj, q->colour) ;
    }

  if (private->maxWidth)
    bsAddData(obj, str2tag("Point_width"), _Int, &private->maxWidth);
  if (private->errorScale != 0.0)
    bsAddData(obj, str2tag("Point_error_scale"), _Float, &private->errorScale);
  if (private->showOrdered)
    bsAddTag(obj, str2tag("Point_segregate_ordered"));
  if (private->showMarginal)
    bsAddTag(obj, str2tag("Point_show_marginal"));
  bsAddTag(obj, str2tag("Point_pne")); controlSetColour(obj, private->pne);
  bsAddTag(obj, str2tag("Point_pe")); controlSetColour(obj, private->pe);
  bsAddTag(obj, str2tag("Point_nne")); controlSetColour(obj, private->nne);
  bsAddTag(obj, str2tag("Point_ne")); controlSetColour(obj, private->ne);
  for(i = 0; i < arrayMax(private->symbolTags); i++)
    bsAddData(obj, str2tag("Point_symbol"), _Text, arrp(private->symbolTags,i,TAG)->text);    
}

typedef struct text
{ char query[280];
  int colour;
} TEXTQUERY;

typedef struct text2
{ char text[280];
} TEXTTAG;

struct configLocals {
  char query[280];
  Array colQuerys; /* query for back ground colours */
  Array tags;     /* tags to get the display name */
  int colWidth;
  float errorScale;
  BOOL segregate,showMultMap;
  int colour[4];
};

static BOOL genesConfigure(COLINSTANCE instance)
{ GENESPRIV private = instance->private;
  int i,maxy;
  float y,height;
  struct  configLocals *cf = (struct configLocals *) messalloc(sizeof(struct configLocals));

  if(arrayMax(private->symbolTags) > arrayMax(private->clonedQuerys))
    height = (50.0 + ((arrayMax(private->symbolTags)+3) * 1.2)) * 0.007 ;
  else
    height = (50.0 + ((arrayMax(private->clonedQuerys)+3) * 1.2)) * 0.007 ;

  if(controlCreateConfig(instance,cf,"Configure point column",0.8,height)){
    
  /*initialise data*/
    if(private->query)
      strcpy(cf->query,private->query);

    cf->colQuerys = arrayCreate(4,TEXTQUERY);

    for (i = 0 ; i < arrayMax(private->clonedQuerys) ; i++)
      { CLONEDQUERY *q = arrp(private->clonedQuerys,i,CLONEDQUERY) ;
	TEXTQUERY *t = arrayp(cf->colQuerys,i,TEXTQUERY) ;
	sprintf (t->query, "%s", q->query);
	t->colour = q->colour;
      }
    strcpy (arrayp(cf->colQuerys, arrayMax(cf->colQuerys), TEXTQUERY)->query, "") ;
    strcpy (arrayp(cf->colQuerys, arrayMax(cf->colQuerys), TEXTQUERY)->query, "") ;
    strcpy (arrayp(cf->colQuerys, arrayMax(cf->colQuerys), TEXTQUERY)->query, "") ;

    cf->tags = arrayCreate(4,TEXTTAG);
    for(i = 0 ; i < arrayMax(private->symbolTags) ; i++)
      { TAG *t = arrp(private->symbolTags,i,TAG);
	TEXTTAG *q = arrayp(cf->tags,i,TEXTTAG);
	sprintf(q->text,"%s",t->text);
      }
    strcpy (arrayp(cf->tags, arrayMax(cf->tags), TEXTTAG)->text, "") ;
    strcpy (arrayp(cf->tags, arrayMax(cf->tags), TEXTTAG)->text, "") ;
    strcpy (arrayp(cf->tags, arrayMax(cf->tags), TEXTTAG)->text, "") ;
    
    cf->colWidth = private->maxWidth;
    cf->errorScale = private->errorScale;
    
    cf->colour[0] = private->pne;
    cf->colour[1] = private->pe;
    cf->colour[2] = private->nne;
    cf->colour[3] = private->ne;
 
    cf->segregate= private->showOrdered;
    cf->showMultMap=private->showMarginal;
    
    y=0.0;

    graphTextEditor("Query for display:",cf->query,280,4,y+=2.0,0);
    for (i=0 ; i < arrayMax(cf->colQuerys) ; i++)
      { y+=1.2;
	if(i)
	  graphTextEditor("                     ",arrp(cf->colQuerys,i,TEXTQUERY)->query,280,4,y,0);
	else
	  graphTextEditor("Query for background:",arrp(cf->colQuerys,i,TEXTQUERY)->query,280,4,y,0); 
	graphColourEditor(" "," ",&arrp(cf->colQuerys,i,TEXTQUERY)->colour,46.0,y);
      }
    maxy = y;

    y=2.0;
    for(i = 0 ; i < arrayMax(cf->tags) ; i++)
      { y+=1.2;
	if(i)
	  graphTextEditor("            ",arrp(cf->tags,i,TEXTTAG)->text,30,51.0,y,0);
	else 
	  graphTextEditor("Symbol Tags:",arrp(cf->tags,i,TEXTTAG)->text,30,51.0,y,0);      
      }
    if(y < maxy)
      y=maxy;

    y+=1.0;
    graphIntEditor("Column width Restriction:",&cf->colWidth,4.0,y+=1.0,0);
    graphFloatEditor("Error Scale:",&cf->errorScale,4.0,y+=1.0,0);

    graphToggleEditor("Segregate ordered genes.",&cf->segregate,4.0,y+=1.0);
    graphToggleEditor("Marginal map data as error.",&cf->showMultMap,4.0,y+=1.0);

    graphColourEditor("Positive data, no error."," ",&cf->colour[0],36.0,y+=-2);
    graphColourEditor("Positive data, error in co-ordinates."," ",&cf->colour[1],36.0,y+=1.2);
    graphColourEditor("Negative data, no error."," ",&cf->colour[2],36.0,y+=1.2);
    graphColourEditor("Negative data, error in coordinates."," ",&cf->colour[3],36.0,y+=1.2);

    graphRedraw();
  }
  return FALSE;
}

static void gmapPointFinal(COLINSTANCE instance, void *locals, BOOL ok)
{ struct configLocals *cf = locals;
  GENESPRIV private = instance->private;
  int i;

  if(ok){
    private->showOrdered = cf->segregate;
    private->showMarginal = cf->showMultMap;
    private->pne = cf->colour[0];
    private->pe = cf->colour[1];
    private->nne = cf->colour[2];
    private->ne = cf->colour[3];
    private->maxWidth = cf->colWidth;
    private->errorScale = cf->errorScale;

    if (private->query)
      messfree(private->query);
    if (private->keyset)
      keySetDestroy(private->keyset);

    if (strlen(cf->query) != 0)
      { private->query = handleAlloc(0, 
				     instance->handle, 
				     1+strlen(cf->query));
	strcpy(private->query, cf->query);
	private->keyset = queryKey(instance->map->key, private->query);
      }
    else
      { private->query = 0;
	private->keyset = 0;
      }

    for (i = 0 ; i < arrayMax(private->clonedQuerys) ; ++i)
      keySetDestroy (arrp(private->clonedQuerys, i, CLONEDQUERY)->keyset) ;
    private->clonedQuerys = arrayReCreate (private->clonedQuerys, 4, CLONEDQUERY) ;

    for (i = 0 ; i < arrayMax (cf->colQuerys) ; i++)
      { TEXTQUERY *t = arrp(cf->colQuerys,i,TEXTQUERY) ;
	char *cp = t->query ;
	while (*cp && *cp == ' ') ++cp ;
	if (*cp)
	  { CLONEDQUERY *q = arrayp(private->clonedQuerys, 
				    arrayMax(private->clonedQuerys),
				    CLONEDQUERY) ;
	    q->query = strnew (cp, instance->handle) ;
	    q->colour = t->colour ;
	    q->keyset = queryKey (instance->map->key, q->query);
	  }
      }

    private->symbolTags = arrayReCreate (private->symbolTags, 4, TAG) ;

    for(i=0;i<arrayMax(cf->tags);i++)
      { TEXTTAG *t = arrp(cf->tags,i,TEXTTAG) ;
	char *cp = t->text ;
	while (*cp && *cp == ' ') ++cp ;		
	if(*cp)
	  { TAG *q = arrayp(private->symbolTags,arrayMax(private->symbolTags),TAG);
	    q->text =  strnew(cp, instance->handle);
	  }
      } 
  }
  else
    { arrayDestroy(cf->colQuerys);
      arrayDestroy(cf->tags);
      messfree(cf);
    }
}

static void genesPrivDestroy(void *p)
{ GENESPRIV private = (GENESPRIV)p;
  
  if (private->keyset)
    keySetDestroy(private->keyset);
}

static BOOL genesCreate(COLINSTANCE instance, OBJ init)
{
  char *s1;
  int i ;
  float f;
  Array flatA;
  GENESPRIV private;

  private = (GENESPRIV)halloc(sizeof(struct GenePrivStruct),
			      instance->handle);
  blockSetFinalise (private, genesPrivDestroy);

  instance->draw = genesDraw;
  instance->setSelectBox = genesSetSelect;
  instance->unSelectBox = genesUnselect;
  instance->save = genesSave;
  instance->doColour = genesDoColour;
  instance->followBox = genesFollowBox;
  instance->configure = genesConfigure;
  instance->private = private;
  instance->configFinal = gmapPointFinal;
  private->query = 0;
  private->keyset = 0;
  private->clonedQuerys = 0;
  private->maxWidth = 0;
  private->errorScale = 0.0;
  private->showOrdered = FALSE;
  private->showMarginal = FALSE;
  private->symbolTags = 0;
  private->pne = GREEN;
  private->pe = RED;
  private->nne = LIGHTBLUE;
  private->ne = RED;
  private->segs = *(instance->proto->convertResults);

  private->clonedQuerys = arrayHandleCreate (4, CLONEDQUERY, instance->handle);
  private->symbolTags = arrayHandleCreate (4, TAG, instance->handle);

  if (init)			/* why should this be empty? */
    {

      if (bsGetData (init, str2tag("Point_query"), _Text, &s1))
	{ private->query = strnew (s1, instance->handle) ;
	  private->keyset = queryKey(instance->map->key, private->query);
	}
      
      if (bsGetData(init, str2tag("Point_yellow"), _Text, &s1))
	{ CLONEDQUERY *q = arrayp (private->clonedQuerys, 
				   arrayMax(private->clonedQuerys),
				   CLONEDQUERY) ;
	  q->query = strnew (s1, instance->handle) ;
	  q->keyset = queryKey (instance->map->key, q->query) ;
	  q->colour = _YELLOW - _WHITE ;
	}

      flatA = arrayCreate(20, BSunit);
      if (bsFindTag (init, str2tag("Point_colour")) && 
	  bsFlatten(init,2,flatA))
        for (i=0 ; i < arrayMax(flatA); i+=2)
	  { CLONEDQUERY *q = arrayp (private->clonedQuerys, 
				     arrayMax(private->clonedQuerys),
				     CLONEDQUERY);
	    q->query = strnew (arr(flatA,i,BSunit).s, instance->handle) ;
	    if (arr(flatA,i+1,BSunit).i)
	      q->colour = arr(flatA,i+1,BSunit).i - _WHITE ;
	    else
	      q->colour = _WHITE ;
	    q->keyset = queryKey(instance->map->key, q->query);
	  }
      arrayDestroy(flatA);

      if (bsGetData(init, str2tag("Point_symbol"), _Text, &s1)) do
	{ TAG *t = arrayp(private->symbolTags,arrayMax(private->symbolTags),TAG) ;
	  t->text = strnew (s1, instance->handle) ;
	} while (bsGetData (init, _bsDown, _Text, &s1)) ;

      if (bsGetData(init, str2tag("Point_error_scale"), _Float, &f))
	private->errorScale = f;
      
      if (bsGetData(init, str2tag("Point_width"), _Int, &i))
	private->maxWidth = i;
  
      if (bsFindTag(init, str2tag("Point_segregate_ordered")))
	private->showOrdered = TRUE;
  
      if (bsFindTag(init, str2tag("Point_show_marginal")))
	private->showMarginal = TRUE;
      
      if (bsFindTag(init, str2tag("Point_pne")))
	private->pne = controlGetColour(init);
      if (bsFindTag(init, str2tag("Point_pe")))
	private->pe = controlGetColour(init);
      if (bsFindTag(init, str2tag("Point_nne")))
	private->nne = controlGetColour(init);
      if (bsFindTag(init, str2tag("Point_ne")))
	private->ne = controlGetColour(init);
	  
    }
  return TRUE;
} /* genesCreate */

static void genesInit(COLPROTO proto, MAPCONTROL map)
{ 
  COLOUROPT *button = (COLOUROPT *) messalloc(sizeof(COLOUROPT));
  
  button->text = "Gmap data...";
  button->f = genesShowDataFromButton,
  button->arg = proto;
  button->fg = BLACK;
  button->bg = WHITE;
  button->next = 0;

  controlAddButton(map, button, gjmMenu); /* TODO tidy up gjmMenu */
  
  messfree(button);

  return;
} /* genesInit */
  

struct ProtoStruct gMapPointColumn = {
  genesInit,
  genesCreate,
  0,
  "Points",
  0,
  FALSE,
  gMapConvert,
  0,
  "The point column can display any point object (ie one with a Position or "
  "Multi_Position tag in the #Map_postion). The set of such objects in a map "
  "is found by scanning all objects two to the right of the Contains tag. This "
  "set is then filtered using the query in the configuration to give the set "
  "of objects to draw. A second query gives a set of objects to draw with "
  "yellow background; this is typically used to show cloned genes.\n\n"
  "Objects are moved rightwards depending on the value of their error. "
  "The configuration parameter \"Error scale\" controls the extent to which "
  "this happens; a value of zero eliminates it. If the \"Segregate ordered "
  "genes\" toggle is set, objects with the Well_ordered tag set are drawn on "
  "the left, others on the right.\n\n"
  "The relationship of point objects to other objects on the map, described "
  "by the Positive and Negative tags, is shown using four colours that can "
  "be picked in the configuration panel. \"Positive data, no error\" is used "
  "when a point is recorded to be inside an interval object and their "
  "repspective co-ordinates agree with this. If the co-ordinates do not "
  "agree, the colour for \"Positive data, error in co-ordinates\" is used. "
  "Similarly, if a point is recorded as being outside an interval, one of the "
  "Negative data colours is used, depending on the co-ordinates. "
  "If the \"Marginal map data as error\" toggle is set, then for positive "
  "data, "
  "to avoid an error, the whole region between a point's error bars must "
  "lie in the containing interval, similarly the whole region between the "
  "error bars must lie outside the specified interval for negative data.\n\n"
  "Note that if no column width limit is set, this column will expand to "
  "fill all available space."
};
 
 
 
