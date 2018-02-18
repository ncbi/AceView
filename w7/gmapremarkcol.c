/*  File: gmapremarkcol.c
 *  Author: Simon Kelley (srk@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: derived column for the genetic map
 * Exported functions:
 *     gMapRemarkColumn
 *     remarkRegister
 * HISTORY:
 * Last edited: Jan  8 16:19 1999 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: gmapremarkcol.c,v 1.3 2014/11/30 03:20:46 mieg Exp $ */

#include "gmap.h"
#include "bump.h"

/************************************************************/

typedef struct {
  KEY key;
  float x;
} REMARKSEG;

typedef struct {
  KEY key;
  KEY parent;
  KEY tag;
  float x;
  int next;
} REMARK;

typedef struct RemarkPrivStruct {
  Array remarks;
  /* identify selected box with these */
  KEY selParent;
  KEY selKey;
  KEY selTag;
  float selx;
  /* configuration */
  int maxWidth;
  char *query;
  Array tags;
  Array tagBool;
  KEYSET keyset;
  COND localQuery;
  BOOL duplicates;
  BOOL showParentsNeighbours;
  BOOL showNeighbours;
  BOOL followParent;
  BOOL hide ;
  BOOL parentQuery;
} *REMARKPRIV;

/************************************************************/

static BOOL remarkSetSelect(COLINSTANCE instance, 
			    int box, double xc, double yc)
{ 
  REMARKPRIV private = instance->private;
  MAPCONTROL map = instance->map;
  GeneticMap look = (GeneticMap)(map->look); 
  void *vr = controlBoxRegd(instance, box);
  int r = assInt (vr) ;
  KEY parent = arr(private->remarks, r, REMARK).parent;
  KEY key = arr(private->remarks, r, REMARK).key;
  KEY tag = arr(private->remarks, r, REMARK).tag;
  float x = arr(private->remarks, r, REMARK).x;

  gMapUnselect(map);
  look->selectKey = key;
  look->neighboursInfoValid = TRUE;
  look->neighbours = keySetReCreate(look->neighbours);
  while (r != -1)
    { keySet(look->neighbours, keySetMax(look->neighbours)) = 
	arr(private->remarks, r, REMARK).parent;
      if (private->showParentsNeighbours)
	{ KEYSET a = bsKeySet(arr(private->remarks, r, REMARK).parent);
	  int i;
	  if (a) 
	    { for ( i = 0; i < keySetMax(a); i++)
		if (class(keySet(a,i)))
		  keySet(look->neighbours, keySetMax(look->neighbours)) = 
		    keySet(a,i);
	      keySetDestroy(a);
	    }
	}
      r = arr(private->remarks, r, REMARK).next;
    }
  if (private->showNeighbours)
    { KEYSET a  = bsKeySet(key);
      int i;
      if (a) 
	{ for ( i = 0; i < keySetMax(a); i++)
	    if (class(keySet(a,i)))
	      keySet(look->neighbours, keySetMax(look->neighbours)) = 
		keySet(a,i);
	  keySetDestroy(a);
	}
    } 
  keySetSort(look->neighbours);
  keySetCompress(look->neighbours);
  map->control->activeKey = key;
  private->selKey = key;
  private->selParent = parent;
  private->selTag = tag;
  private->selx = x;

  return FALSE;
} /* remarkSetSelect */


static BOOL remarkUnselect(COLINSTANCE instance, int box)
{ 
  GeneticMap look = (GeneticMap)(instance->map->look);

  
  *look->messageText = 0;

  if (look->messageBox)
    graphBoxDraw(look->messageBox, -1, -1);
  
  return FALSE;
} /* remarkUnselect */


static void remarkFollowBox(COLINSTANCE instance, int box, double x, double y)
{
  void *vr = controlBoxRegd(instance, box);
  int r = assInt (vr) ;
  REMARKPRIV private = instance->private;
  KEY parent = arr(private->remarks, r, REMARK).parent;
  KEY symb = arr(private->remarks, r, REMARK).key;
  
  if (private->followParent)
    display(parent, instance->map->key, TREE);
  else
    display(symb, instance->map->key, TREE);

  return;
} /* remarkFollowBox */


static void remarkDoColour(COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);
  REMARKPRIV private = instance->private;
  void *vr = controlBoxRegd(instance, box);
  int r = assInt (vr) ;
  KEY parent = arr(private->remarks, r, REMARK).parent;
  KEY key = arr(private->remarks, r, REMARK).key;
  KEY tag = arr(private->remarks, r, REMARK).tag;
  float x = arr(private->remarks, r, REMARK).x;
  char *mt = look->messageText;
  
  if (instance == control->activeInstance &&
      tag == private->selTag &&
      parent == private->selParent &&
      key == private->selKey &&
      x == private->selx)
    { graphBoxDraw(box, BLACK, CYAN);
      *mt = 0;
      gMapAddToHeader(map, name(tag));
      if (private->followParent)
	{ gMapAddToHeader(map, " in ");
	  gMapAddToHeader(map, name(parent));
	}
      gMapAddToHeader(map, ": ");
      gMapAddToHeader(map, name(key));
      if (look->messageBox)
	graphBoxDraw(look->messageBox, -1, -1);
      control->activeBox = box;
    }
  else
    { 
      if (control->activeInstance && 
	  instance->map == control->activeInstance->map &&
	  look->neighboursInfoValid && 
	  keySetFind(look->neighbours, key, 0))
	graphBoxDraw(box, BLACK, PALECYAN);
      else
	graphBoxDraw(box, BLACK, WHITE);
    }
} /* remarkDoColour */
 

void remarkRegister(MAPCONTROL map, KEY key, float coord)
{
  GeneticMap look = (GeneticMap)(map->look);
  Array segs = look->remarkSegs;
  REMARKSEG *seg;
  if (segs) /* may be no remark column */
    { seg = arrayp(segs, arrayMax(segs), REMARKSEG);
      seg->key = key;
      seg->x = coord;
    }

  return;
} /* remarkRegister */


static int remarkOrder(const void *a, const void *b)
{
  float x1 = ((const REMARKSEG *)a)->x;
  float x2 = ((const REMARKSEG *)b)->x;
  KEY k1 = ((const REMARKSEG *)a)->key;
  KEY k2 = ((const REMARKSEG *)b)->key;

  if (x1 > x2) return 1;
  if (x1 < x2) return -1;
  
  /* co-ord equal, order on key */
  if (k1 > k2) return 1;
  if (k1 < k2) return -1;

  return 0;
} /* remarkOrder */


static void remarkDraw(COLINSTANCE instance, float *offset)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  GeneticMap look = (GeneticMap)(map->look);
  REMARKPRIV private = instance->private;
  Array segs = look->remarkSegs;
  BUMP bump;
  KEY lastKey;
  float lastx=0;
  int width;
  int usedWidth;
  int i, j, ibox;
  int firstSeg = 0, lastSeg = 0, countDir;
  BOOL isRemark = FALSE;
  Associator dupStomp=0;

  if (!arrayExists(segs) || !arrayMax(segs))
    return ;

  if (!private->duplicates)
    dupStomp = assCreate();

  if (segs)
    {
      if (map->mag > 0)
	{
	  firstSeg = 0;
	  lastSeg = arrayMax(segs);
	  countDir = 1;
	}
      else
	{
	  firstSeg = arrayMax(segs)-1;
	  lastSeg = -1;
	  countDir = -1;
	}
      
      arraySort(segs, remarkOrder);
    }

  if (private->maxWidth)
    width = private->maxWidth;
  else
    width = 1+ control->graphWidth - *offset;

  bump = bumpCreate(width, 0);
  usedWidth = 0;

  private->remarks = arrayReCreate(private->remarks, 50, REMARK);
  lastKey = 0; /* to remove duplicates */

   for (i = firstSeg ; i != lastSeg ; i += countDir)
    { OBJ obj;
      REMARKSEG *seg = arrp(segs, i, REMARKSEG);
      if (seg->key == lastKey && seg->x == lastx)
	continue; /* ignore duplicates */
      lastKey = seg->key;
      lastx = seg->x;
      if (private->keyset && !keySetFind(private->keyset, seg->key, 0))
	continue; 
      obj = bsCreate(seg->key);
      if (!obj)
	continue;
      for (j = 0; j<arrayMax(private->tags); j++)
	{ KEY tag = arr(private->tags, j, KEY);
	  KEY symb;
	  if (!arr(private->tagBool, j, BOOL))
	    continue;
	  if (bsGetKey(obj, tag, &symb))
	    do
	      { int r, x = 0;
		int nameLength;
		float y;

		/* mieg
                   in case !private->parentQuery, run the query locally 
		   cannot use queryFind3 because the query may be
                   Follow, (and not just a filter) 
		   */
		if (private->query && !private->parentQuery)
		  {
		    KEYSET kstest = queryKey (symb, private->query) ;
		    int n = keySetMax (kstest) ;
		    keySetDestroy (kstest) ;
		    if (!n) continue ;
		  }

		r = arrayMax(private->remarks);
		array(private->remarks, r, REMARK).key = symb;
		array(private->remarks, r, REMARK).parent = seg->key;
		array(private->remarks, r, REMARK).tag = tag;
		array(private->remarks, r, REMARK).x = seg->x;
		array(private->remarks, r, REMARK).next = -1;
		
		if (!private->duplicates)
		  { int n;
		    void *vs;
		    if (assFind(dupStomp, assVoid(symb), &vs))
		      { n = assInt (vs) ;
			while (arr(private->remarks, n, REMARK).next != -1)
			  n = arr(private->remarks, n, REMARK).next;
			arr(private->remarks, n, REMARK).next = r;
			continue;
		      }
		    assInsert(dupStomp, assVoid(symb), assVoid(r));
		  }
		
		y = MAP2GRAPH(map, seg->x);
		if (y < control->topMargin + 0.5)
		  y = control->topMargin + 0.5;
		nameLength = strlen(name(symb));
		bumpItem(bump, nameLength+1, 1, &x, &y);
		if ( y<control->graphHeight - 1.0)
		  { char buff[1000];
		    remarkRegister(map, symb, seg->x); /* tricky, watch this */
		    /* Don't use seg now, it may be invalid */ 
		    if (!private->hide)
		      { controlRegBox(instance, ibox = graphBoxStart(), 
				      assVoid(r));
			graphBoxInfo (ibox, symb, 0); /* WWW */
			buff[0] = 0;
			strncat(buff, name(symb), width - x);
			graphText(buff, *offset+x, y);
			graphBoxEnd();
			isRemark = TRUE;
			if (nameLength + x > usedWidth)
			  { if (nameLength + x > width)
			      usedWidth = width;
			  else
			    usedWidth = nameLength + x;
			  }
		      }
		  }
	      }
	    while (bsGetKey(obj, _bsDown, &symb));
	}
      bsDestroy(obj);
    }

  if (isRemark)
    *offset += usedWidth + 1;
  bumpDestroy(bump);
  if (!private->duplicates)
    assDestroy(dupStomp);

  return;
} /* remarkDraw */


static void remarkDrawInit(COLINSTANCE instance)
     /* this is a COLINSTANCE->drawInit() function */
{
  GeneticMap look = (GeneticMap)(instance->map->look);

  look->remarkSegs = arrayReCreate(look->remarkSegs, 50, REMARKSEG);
} /* remarkDrawInit */


struct configLocals {
  char query[280];
  int width;
  BOOL duplicates,parentQuery,showNeighbours,showParentsNeighbours,followParent,hide;
  char *temp[5];
  BOOL toggle[5];
};
  
static BOOL remarkConfigure(COLINSTANCE instance)
{
  REMARKPRIV private = instance->private;
  int i;
  struct  configLocals *cf = (struct configLocals *) messalloc(sizeof(struct configLocals));
  float line = 6.0;

  if(controlCreateConfig(instance,cf,"Configure remark column",0.8,0.3)){
    
  /*initialise data*/
    if(private->query)
      strcpy(cf->query,private->query);
    for(i = 0; i!=5; i++){
      cf->temp[i] = (char *) messalloc(sizeof(char)*80);
      if(i < arrayMax(private->tags)){
	strcpy(cf->temp[i],name(arr(private->tags, i, KEY))); 
	cf->toggle[i] =arr(private->tagBool, i, BOOL); 
      }
    }
    cf->duplicates = private->duplicates;
    cf->parentQuery = private->parentQuery;
    cf->showNeighbours = private->showNeighbours;
    cf->showParentsNeighbours = private->showParentsNeighbours;
    cf->followParent = private->followParent;
    cf->hide = private->hide;
    cf->width = private->maxWidth;

    graphTextEditor("Query for display:",cf->query,280,4,2,0);
    graphIntEditor("Column width Restriction:",&cf->width,4.0,4.0,0);
    graphText("Tags to show:", 1, line++);
    for(i = 0; i!=5; i++)
      { graphTextEditor(" ",cf->temp[i],80,4.0,line,0);
	graphToggleEditor("",&cf->toggle[i],2.0,line++);
	line++;
      }
    line = 6.0;
    graphToggleEditor("Query applies to parent",&cf->parentQuery,40.0,line++);
    graphToggleEditor("Display duplicates",&cf->duplicates,40.0,line++);
    graphToggleEditor("On double-click, show parent",&cf->followParent,40.0,line++);
    graphToggleEditor("Hide this column",&cf->hide,40.0,line++);
    line++;
    graphText("On selecting, highlight:", 39.0,line++);
    graphToggleEditor("Object's neighbours",&cf->showNeighbours,40.0,line++);
    graphToggleEditor("Parent's neighbours",&cf->showParentsNeighbours,40.0,line++);
    
    graphRedraw();
  }

  return FALSE;
} /* remarkConfigure */


static void gmapRemarkFinal(COLINSTANCE instance, void *locals, BOOL ok)
{
  struct configLocals *cf = locals;
  REMARKPRIV private = instance->private;
  int i;

  if(ok){
    private->duplicates=  cf->duplicates;
    private->parentQuery = cf->parentQuery ;
    private->showNeighbours =cf->showNeighbours ;
    private->showParentsNeighbours = cf->showParentsNeighbours ;
    private->followParent = cf->followParent ;
    private->hide = cf->hide;
    private->maxWidth = cf->width;
    
    if (private->keyset)
      keySetDestroy(private->keyset);
    private->keyset = 0;
    if (private->query)
      messfree(private->query);
    private->query = 0;
    
    if (strlen(cf->query) != 0)
      { private->query = strnew(cf->query, instance->handle);
	if (private->parentQuery)
	  private->keyset = queryKey(instance->map->key, private->query);
      }
    
    private->tags = arrayReCreate(private->tags, 10, KEY);
    private->tagBool = arrayReCreate(private->tagBool, 10, BOOL);
    for (i=0; i!=5; i++)
      { if (strlen(cf->temp[i]) != 0)
	  { KEY t;
	    lexaddkey(cf->temp[i], &t, 0);
	    array(private->tags, arrayMax(private->tags), KEY) = t;
	    array(private->tagBool, arrayMax(private->tagBool), BOOL) = cf->toggle[i];
	  }
      }
  }
  else
    messfree(cf);

  return;
} /* gmapRemarkFinal */

      
static void remarkSave(COLINSTANCE instance, OBJ obj)
{ 
  REMARKPRIV private = instance->private;
  int i, d;
  
  if (private->query)
    bsAddData(obj, str2tag("DT_query"), _Text, private->query);
  if (private->maxWidth)
    bsAddData(obj, str2tag("DT_width"), _Int, &private->maxWidth);
  if (!private->duplicates)
    bsAddTag(obj, str2tag("DT_no_duplicates"));
  if (private->showNeighbours)
    bsAddTag(obj, str2tag("DT_neighbours"));
  if (private->showParentsNeighbours)
    bsAddTag(obj, str2tag("DT_parents"));
  if (private->followParent)
    bsAddTag(obj, str2tag("DT_follow_parent"));
  if (private->hide)
    bsAddTag(obj, str2tag("DT_hide"));
  if (!private->parentQuery)
    bsAddTag(obj, str2tag("DT_symbol_query"));
  for (i=0; i<arrayMax(private->tags); i++)
    { bsAddData(obj, str2tag("DT_tag"), _Text, 
		name(arr(private->tags, i, KEY)));
      d = arr(private->tagBool, i, BOOL) ? 0 : 1;
      /* zero means display, for backwards compatibility */
      if(bsType(obj,_bsRight)==_Int) /* check to make sure Int can be saved */
	bsAddData(obj, _bsRight, _Int, &d);
    }

  return;
} /* remarkSave */

static void remarkPrivDestroy(void *p)
     /* blockFinalisation function for REMARKPRIV-type */
{
  REMARKPRIV private = (REMARKPRIV)p;

  if (private->keyset)
    keySetDestroy(private->keyset);
} /* remarkPrivDestroy */


static BOOL remarkCreate(COLINSTANCE instance, OBJ init)
{
  REMARKPRIV private;
  char *s1;
  int i;
  Array tags = arrayCreate(50, BSunit);

  private = (REMARKPRIV)halloc(sizeof(struct RemarkPrivStruct),
			       instance->handle);
  blockSetFinalise (private, remarkPrivDestroy);

  instance->setSelectBox = remarkSetSelect;
  instance->unSelectBox = remarkUnselect;
  instance->followBox = remarkFollowBox;
  instance->doColour = remarkDoColour;
  instance->draw = remarkDraw;
  instance->drawInit = remarkDrawInit;
  instance->save = remarkSave;
  instance->configure = remarkConfigure;
  instance->configFinal = gmapRemarkFinal;

  instance->private = private;
  private->remarks = arrayHandleCreate (50, REMARK, instance->handle);
  private->tags = arrayHandleCreate(10, KEY, instance->handle);
  private->tagBool = arrayHandleCreate(10, BOOL, instance->handle);
  private->maxWidth = 0;
  private->duplicates = FALSE;
  private->showParentsNeighbours = FALSE;
  private->showNeighbours = FALSE;
  private->followParent = FALSE;
  private->hide = FALSE;
  private->parentQuery = TRUE;

  private->keyset = 0;
  private->query = 0;

  if (init)
    { if (bsGetData(init, str2tag("DT_width"), _Int, &i))
	private->maxWidth = i;

      if (bsFindTag(init, str2tag("DT_no_duplicates")))
	private->duplicates = FALSE;
      else
	private->duplicates = TRUE;
      
      if (bsFindTag(init, str2tag("DT_neighbours")))
	private->showNeighbours = TRUE;
	  
      if (bsFindTag(init, str2tag("DT_parents")))
	private->showParentsNeighbours = TRUE;

      if (bsFindTag(init, str2tag("DT_follow_parent")))
	private->followParent = TRUE;

      if (bsFindTag(init, str2tag("DT_hide")))
	private->hide = TRUE;

      if (bsFindTag(init, str2tag("DT_symbol_query")))
	private->parentQuery = FALSE;
      
      if (bsGetData(init, str2tag("DT_query"), _Text, &s1))
	{ private->query = strnew(s1, instance->handle);
	  if (private->parentQuery)
	    private->keyset = queryKey(instance->map->key, private->query);
	}

      if (bsFindTag(init, str2tag("DT_tag")) &&
	  bsFlatten(init, 2, tags))
	for (i = 0; i < arrayMax(tags); i += 2 )
	  { KEY t;
	    lexaddkey(arr(tags, i, BSunit).s, &t, 0);
	    array(private->tags, arrayMax(private->tags), KEY) = t;
	    array(private->tagBool, arrayMax(private->tagBool), BOOL) =
	      ((arr(tags, i+1, BSunit).i == 0) ? TRUE : FALSE);
	  } 
    }

  arrayDestroy(tags);

  return TRUE;
} /* remarkCreate */


struct ProtoStruct gMapRemarkColumn = { 
  0,
  remarkCreate,
  0,
  "Derived_tags",
  0,
  FALSE,
  0,
  0,
  "This column displays information derived from the set of all "
  "objects drawn in columns to the left of it. "
  "This set is first filtered against the keyset resulting from the execution "
  "of the display query, which starts from the map object. These objects are "
  "then searched for the presence of the tags listed in the configuration. "
  "When a tag is found the name of the object to the right of it is displayed "
  "in the derived tags column, at the same position (modulo bumping) as its "
  "parent.\n\n"
  "If Query applies to parent is not set, a different procedure occurs: in "
  "this case the set of objects to be displayed is derived as above, but "
  "without the filtering step, and then each of these objects is matched "
  "against the query.\n\n"
  "If \"display duplicates\" is set, the same object name may appear many "
  "times in different positions if it appears in many parent objects drawn "
  "in other columns. Without \"display duplicates\", second and subsequent "
  "drawings of object names are suppressed and only the one nearest the top "
  "of the display remains.\n\n"
  "When an object in a derived tag column is selected (By clicking on it.) "
  "its parent(s) are coloured pale cyan, optionally it neighbours (All those "
  "objects which it contains pointers to.) are also coloured, and also, "
  "optionally the neighbours of its parents. These options are selected in "
  "the configuration window.\n\n" 
  "The behaviour of objects when double clicked is also controlled in the "
  "configuration window. Either the object itself may be displayed, or its "
  "parent. Note that it does not make much sense to display the parent unless "
  "Display Duplicates is set, as without this, an object may not have a "
  "unique parent.\n\n"
  "Hide this column stops the column from being drawn, but allows the objects "
  "which it would have drawn to be available as parents for other columns "
  "drawn further right.\n"

}; 

/**************************** eof ***************************/
