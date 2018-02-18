/*  File: qbedisp.c
 *  Author: Gary Aochi (aochi@genome.lbl.gov)
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm1.cnusc.fr
 *
 * Description:
 **  Constructs general queries.                              **
     Gary Aochi modified version of qbedisp.c
     with on line help to preform complex queries.
 *     _DDtQuery   -g TEXT_SCROLL -t "Query" -w .73  -height 0.35 -help Query
 *
 * Exported functions:  qbeCreate
 * HISTORY:
 * Last edited: Nov 19 15:03 1998 (edgrif)
 * * Nov 19 15:02 1998 (edgrif): Fixed callback func. dec. for
 *              graphTextScrollEditor.
 * Last edited:Jun 6 11:00 1996 (rbrusk):
 *	-	graphTextScrollEntry() callback function (GraphFunc) casts removed
 * * Jun 4 15:49 1996 (rbrusk):
 *	-	changed fn() => fn(void) function declarations to please compiler
 * * Jun  3 12:11 1996 (rd)
 * * Jan 18 12:55 1993 (mieg): Changed MAGIC and added some static def
 * Created: Feb 09 15:30 1993 (gha): split code from querydisp.c
 *-------------------------------------------------------------------
 * gha 03/29/93  changed querEdit() to queryCommandEdit(), fixed execution
 * gha 04/02/93  redid tags as pickable, reverse-video, but need duplicate ctrl
 * gha 04/05/93  fix mistake for tag = * (boolean), added qbepick-qbeform call.
 * gha 04/08/93  improved AND/OR option and better query command formatting
 * gha 04/12/93  Added & NEXT capability with Int/Float checking for "
 * gha 04/29/93  jean cleaned up stuff, fixed classlist, adding item name box
 *               added code to unhighlight some buttons
 * gha 05/25/93  jean changes plus my own for formation and classlist ;
 * gha 05/28/93  added metadata function
 * gha 06/08/93  changed "Find" to "Search"
 * gha 06/10/93  added feature to search from active keyset
 * gha 06/14/93  split querydisp.c and querybuild.c from querydisp.c
 * gha 06/16/93  cleaned up #hash (not shown at all now), rather than
 *               badly displayed as a box.  Should try to properly fix later.
 * gha 6/23/93   moved metadata stuff to metadisp.c
 * kliu 9/2/93   added triangle box for class selection
 */

/* $Id: qbedisp.c,v 1.7 2017/03/18 15:30:51 mieg Exp $ */
/*
#define ARRAY_CHECK
#define MALLOC_CHECK
*/

#define DEFINE_OBJ
typedef struct sobj *OBJ ;

#include "acedb.h"
#include "bs_.h"        /* 10/12/92 gha */
#include "query_.h"

extern BOOL metaGetData(int metaclass, KEY classkey, BS tag) ;

static int  QBE_MAGIC =  29986 ;
static int nameBox ;
/* Query by example  typedef */
typedef struct QBESTUFF
{ int   magic ;        /* == QBE_MAGIC */
  KEY   key ;          /* key for model-class */
  OBJ   objet ;        /* data stored in structure */
  Graph graph ;        /* window graph id for QBE window */
  int   activebox ;    /* current active box in window */
  Array tagval ;       /* storage of tags and such for QBE interface */
  char* text ;	       /* for typing in things */
  int classe ;          /* for completion mechanism */
} *QBE ;

Graph qbeGraph = 0 ;

#include "display.h"	/* must come after QUER definition */

void qbeCreate (void); 
static void qbeDestroy (void) ;
static char conjunction_op = '&';     
    

typedef struct TAGVALUESTUFF
{
  int  box;                /* current box id */
  char name[BUFFER_SIZE];  /* tag name or entry value */
  int  textformat;         /* for textformat of output (bold/plain) */
  BOOL tag_active;         /* used during query formation */
  BOOL tag_defined;        /* is the tag properly set up; valid? */
  char tv_type;            /* 'char' to save bits; TV_TAG, TV_TEXT_ENTRY.. */
  struct TAGVALUESTUFF *tag;  /* for value :  backwards pointer to tag */
  int  has_value ;         /* for tags...true if has value */
  int  value_number ;      /* for value enumeration */
  BS   bs ;                /* for metadata */
} *TV;


/* types for tv_type field of tagvaluestuff */
#define TV_TAG           0
#define TV_TEXT_ENTRY    1
#define TV_NUMERIC_ENTRY 2
#define TV_DATE_ENTRY 3

/* background and foreground colors */
#define NORMAL_FCOL BLACK
#define NORMAL_BCOL WHITE
#define SELECT_FCOL WHITE
#define SELECT_BCOL BLACK


static TV lasttag ;         /* for qbeDrawBS */
static int lasttag_idx;

static BOOL xref_tag ;      /* for qbeDrawBS */

/* Item Name box */
static int itemBox, itemEntryBox;
static BOOL item_enable ;
static char item_buffer[BUFFER_SIZE] ;

static BOOL from_keyset ;   /* for searching from active keyset */

static int boxMax;

static FREEOPT *classesMenu;

static int qbe_class;
static int qbe_classBox, qbe_classEntryBox, qbe_classTriBox;
static char qbe_class_buffer[BUFFER_SIZE];

static void qbeDisplay(void);
static void qbePick(int box);
static void qbeReKey(void);
static void qbeClear(void) ;
static void qbeMenuExec(void);  /* for QBE_Menu */
static void qbeMenuQuit(void);
static void qbeMenuDraw(void);
static void qbeMenuAnd(void);
static void qbeMenuOr(void);
/*
static void qbeMenuMetaLong(void);
static void qbeMenuMetaFull(void);
*/

static void qbeButtonUndo(void) {
  queryCommandUndo();
  qbeDisplay();  /* gha added to unhighlight undo button 4/29/93 */
}

/* call up query-commander */
static void queryCallup(void) {
  queryCreate();
  qbeDisplay();  /* unhighlight button */
}

/* call up query builder */
static void qbuildCallup(void) {
  qbuildCreate();
  qbeDisplay();  /* unhighlight button */
}

/* execute query */
static void qbeButtonExec(void) {
  qbeMenuExec();
  qbeDisplay();  /* unhighlight button */
}

static void menu_spacer(void) { return; }

static void myNewKeySet (void)
{ newKeySet ("QBE Answer") ;
}

static MENUOPT QBE_Menu[]=
{
  {qbeMenuQuit, "Quit"},
  {help, "Help"},
  {graphPrint, "Print"},
  {qbeMenuDraw, "Refresh"},
  {qbeClear, "Clear"},
  {menu_spacer,""},
  {myNewKeySet, "New KeySet"},
  {qbuildCreate, "Query Builder"},
  {queryCreate, "Query Command Editor"},
  {menu_spacer,""},
  {qbeMenuAnd, "AND"},
  {qbeMenuOr, "OR"},
  {menu_spacer, ""},
/*
  {qbeMenuMetaFull, "Show Full Description"},
  {qbeMenuMetaLong, "Show Long Description"},
*/
  {menu_spacer, ""},
  {qbeMenuExec, "Search"},
  {queryCommandUndo, "Undo"},
   {0,0}
  } ;

static MENUOPT QBE_Buttons[]=
{
  {qbeMenuQuit, "Quit"},
  {help, "Help"},
  {qbeClear, "Clear"},
  {myNewKeySet, "New KeySet"},
  {qbuildCallup, "Builder"},
  {queryCallup, "Command Edit"},
  {qbeButtonExec, "Search"},
  {qbeButtonUndo, "Undo"},
   {0,0}
  } ;

#define QBEGET(name)     QBE qbelook ; \
  if (!graphActivate(qbeGraph)) \
    messcrash("%s can't find qbeGraph",  name); \
  if (!graphAssFind (&QBE_MAGIC,&qbelook)) \
    messcrash ("%s can't find graph",name) ; \
  if (!qbelook) \
    messcrash ("%s received a null pointer",name) ; \
  if (qbelook->magic != QBE_MAGIC) \
    messcrash ("%s received a wrong pointer",name)

/**************************************************************/


static Array classlist = 0;

/* retrieve a class name from db list, note that classlist is different
 * in qbedisp.c and querydisp.c
 */
#define CLASSLIST_NAME(x)  ARR2STRING(classlist, x)

static int graphDownSelTriangle(BOOL filled, float x, float y)
{
  int i = 3, box;
  Array temp ;

  box = graphBoxStart();

  graphLine(x, y+0.25, x+1, y+0.25);
  graphLine(x+1, y+0.25, x+0.5, y+0.75);
  graphLine(x+0.5, y+0.75, x, y+0.25);

  if (filled) {
    temp = arrayCreate(2*i, float) ;
    array(temp, 0, float) = x;
    array(temp, 1, float) = y+0.25;
    array(temp, 2, float) = x+1;
    array(temp, 3, float) = y+0.25;
    array(temp, 4, float) = x+0.5;
    array(temp, 5, float) = y+0.75;
    graphPolygon(temp);
    arrayDestroy(temp);
  }

  graphBoxEnd();

  return box;
}

/* find the total number of listed classes in the menu, if no 
 * list exists, initialize it and return the total number.
 */
static int qbeClassMax(void) 
{
  int n, i, k ;
  char *v, *v0 = 0 ;

  if (graphAssFind(&qbe_class, &v))
    return v - v0 ;

  n = 0;
  classlist = arrayCreate(30, Array);   /* init classlist */
  
  for (i = 1 ; i <= 256 ; i++) {
    if (pickList[i].name  && pickList[i].type == 'B' &&
	pickList[i].visible != 'h') {
      array(classlist, n, Array) = arrayCreate(BUFFER_SIZE, char);
      k = KEYMAKE(i, 0) ;
      strcpy(ARR2STRING(classlist, n), className(k));
      n++;
    }
  }

  arraySort(classlist, textAlphaOrder);

  v = v0 + n ;
  graphAssociate(&qbe_class, v);
  return n;
}


static void destroyClassesMenu(void)
{ int max;
  if (graphAssFind(&classesMenu, &max))  {
    destroyMenu(classesMenu, max);
    graphAssRemove(&classesMenu);
  }
}

/* initialize the classes menu */
static void classesMenuSet(void)
{ char *v, *v0 = 0 ;
  int i, max;
  
  max = qbeClassMax();
  
  /* max should be one larger for header-entry */
  classesMenu = (FREEOPT *)halloc((max + 1) *sizeof(FREEOPT), graphHandle() );
  classesMenu->key = max;
  classesMenu->text = halloc(BUFFER_SIZE, graphHandle()) ;
  
  strcpy(classesMenu->text, "CLASSES");
  
  for (i = 1; i <= max; i++)  {
    (classesMenu + i)->key = i;
    (classesMenu + i)->text = halloc(BUFFER_SIZE, graphHandle()) ;
    strcpy((classesMenu + i)->text, CLASSLIST_NAME(i - 1));
  }
  
  graphAssociate(&classesMenu, v = v0 + max);
}


/* handle the "RIGHT-MOUSE" (menu) selection of the class menu box */
static void classesMenuResult(KEY k, int box)
{
  int temp_class;

  temp_class = pickWord2Class((classesMenu+k)->text);

  if (pickType(KEYMAKE(temp_class, 0)) == 'B') {
    qbe_class = temp_class;
    strcpy(qbe_class_buffer, (classesMenu+k)->text);
    /*    graphTextScrollEntry(qbe_class_buffer, 0, 0, 0, 0, 0); */
    boxMax = 0 ; /* mhmp 17.09.98 */
    qbeReKey();
  }
  else {
    messout("Sorry, class %s has no model to display.", (classesMenu+k)->text);
    qbeDisplay();
  }
}     


static char *currentClassName(void) 
{
  return className(KEYMAKE(qbe_class,0));
}

static BOOL qbe_classEntryHandler(char *text, int unused)
{
  int i, temp_class;
  Graph old = graphActive () ;
  for (i = 1; i <= qbeClassMax(); i++)
    if (!strcasecmp((classesMenu+i)->text, text))
      break;

  if (i > qbeClassMax()) {
    messout("Invalid class name %s.", text);
    strcpy(text, pickClass2Word(qbe_class));
    /*    graphTextScrollEntry(text, 0, 0, 0, 0, 0); */
    graphActivate (old) ;
    return FALSE ;
  }

  temp_class = pickWord2Class((classesMenu+i)->text);
  if (pickType(KEYMAKE(temp_class, 0)) == 'B') {
    qbe_class = temp_class;
    qbeReKey();
  }
  else {
    messout("Chosen class %s has no model.", (classesMenu+i)->text);
    strcpy(text, pickClass2Word(qbe_class));
    qbeDisplay();
  }
  graphActivate (old) ;
  return TRUE ;
}

static FREEOPT fromMenu[] = 
{
  {2, "QUERY ORIGINS"},
  {1, "Active Keyset"},
  {2, "Specified Class"}
};

static void qbeShowQueryCommand(void) ;

static void fromMenuResult(KEY k, int box)
{
  if (k == 1)
    from_keyset = TRUE;
  else 
    from_keyset = FALSE;
  qbeShowQueryCommand();
  qbeDisplay();
}

static char *currentFromName(void) 
{
  if (from_keyset)
    return "KeySet";
  else
    return "Class";
}

static void fromButtonResult(void) 
{
  if (from_keyset)
    from_keyset = FALSE;
  else
    from_keyset = TRUE;
  qbeShowQueryCommand();
  qbeDisplay();
}     


/************************************************************/
static void qbeClear(void) /* mhmp 18.09.98 */
{
  QBEGET("qbeClear") ; 
  if (qbelook && qbelook->key)
    qbelook->key = 0 ;
  qbeReKey() ;
}

/* reset the display for a new class/model */
static void qbeReKey(void)
{ int oldKey = 0 ;
  QBEGET("qbeReKey") ; 
  if (qbelook && qbelook->key)
    oldKey = qbelook->key  ;
  qbelook->key = KEYMAKE(qbe_class, 0);
  if (oldKey == qbelook->key)
    return ;
  boxMax = 0 ;
  bsDestroy (qbelook->objet);
  
  qbelook->tagval = arrayReCreate(qbelook->tagval, 1, TV); /*mhmp 05.11.98 */
  array(qbelook->tagval, 0, TV) = NULL;

  if (pickType (qbelook->key) != 'B' 
      || 
      !(qbelook->objet = bsCreate(qbelook->key)))  {
    if(graphActivate(qbeGraph))  {
      graphDestroy () ;
      qbeGraph = 0 ;
    }
    messout("Bad model for class %s.\n  Removing QBE Window.", 
	    pickClass2Word(qbe_class));
    return;
  }
  graphRetitle (messprintf("Query By Example: %s", 
			   className(qbelook->key))) ;
  
  array(qbelook->tagval, 0, TV) = NULL;
  lasttag = NULL;
  lasttag_idx = 0;
  xref_tag = FALSE;
  /*  boxMax = 0; mhmp 16.09*/
  
  qbeDisplay();
}  


/* Destroys qbe data structure, but not qbe graph  */
static void qbeDestroy (void)
{
  QBE qbelook;
  int i;
  
  if (graphAssFind (&QBE_MAGIC,&qbelook)
      && qbelook
      && qbelook->magic == QBE_MAGIC)
    {
      bsDestroy (qbelook->objet);
          
      if (!graphExists(qbuildGraph) && !graphExists(qbuildGraph))
	arrayDestroy(query_undo_set);
      arrayDestroy(qbelook->tagval);
      qbelook->magic = 0 ;
      
      if (querGraph == 0 || !graphExists(querGraph)) 
	{
	  for (i = 0; i < qbeClassMax(); i++)
	    arrayDestroy(arr(classlist, i, Array));
	  arrayDestroy(classlist);
	}
      
      if (0) /* 2013_01_16   crashes if we init the db and immediatly quit query-by-example */ 
	destroyClassesMenu();
      i = 1 ;
      messfree (qbelook) ;
    }
}

/*************************************************************/

/* QBE Menu responses */

static char *qbeFormQuery(char *buffer)
{
  int i;
  TV tagval, oldtag;
  BOOL notfirst = FALSE, tagnotvalue = TRUE;
  QBEGET("qbeFormQuery") ;

  if (from_keyset)
    strcpy(buffer, "");
  else
    strcpy(buffer, messprintf("Find %s",pickClass2Word(qbe_class)));

  if (item_enable) {
    strcat(buffer, messprintf(" %s", item_buffer));
    notfirst = TRUE ;
  }

  tagval = NULL;

  for (i = 0; i < arrayMax(qbelook->tagval); i++)  {
    oldtag = tagval;
    tagval = array(qbelook->tagval, i, TV);

    /* check if is an entry (rather than a tag (if tag, tagval->tag = NULL)) */
    if (tagval && oldtag != tagval) 
      {
	if (tagval->tag_active) {   /* tag is active */
	  tagnotvalue = TRUE;
	  strcat(buffer,
		 (notfirst ?
		  messprintf(" %c %s", conjunction_op, tagval->name) :
		  messprintf(" %s", tagval->name)));
	  notfirst = TRUE;
	}
	else if (tagval->tag && tagval->tag->tag_active)
	  {  /* value part */
	    if (!strcmp(tagval->name, "*") && !tagnotvalue) {      /* just add a tag (boolean) */
	      strcat(buffer, messprintf(" %c NEXT", conjunction_op));
	      tagnotvalue = FALSE;
	    }
	    else if (strcmp(tagval->name,""))    /* add an  'tag equals value' */
	      { char *cq = tagval->name ; /* mieg allow comparators in qbe */

		while (*cq == ' ') cq++ ;
		if (*cq == '>' || *cq == '<' || *cq == '=') { /* jtm */
		  strcat(buffer, ((tagnotvalue) ? messprintf(" %s", tagval->name) :
				  messprintf(" %c NEXT %s", conjunction_op, tagval->name)));
		  tagnotvalue = FALSE;
		}
		else {
		  strcat(buffer, ((tagnotvalue) ? 
				  ((tagval->tv_type == TV_TEXT_ENTRY) ?
				   messprintf(" = \"%s\"", tagval->name) :
				   messprintf(" = %s", tagval->name)) :
				  ((tagval->tv_type == TV_TEXT_ENTRY) ?
				   messprintf(" %c NEXT = \"%s\"", conjunction_op, 
					      tagval->name) :
				   messprintf(" %c NEXT = %s", conjunction_op, 
					      tagval->name))));
		  tagnotvalue = FALSE;
		}
	      }
	    else if (tagval->value_number < tagval->tag->has_value) {  /* filler */
	      strcat(buffer, messprintf(" %c NEXT", conjunction_op));
	      tagnotvalue = FALSE;
	    }
	  }
      }
  }
  return buffer;
}

static void qbeMenuAnd(void) {
  conjunction_op = '&';
  qbeShowQueryCommand();
}

static void qbeMenuOr(void) {
  conjunction_op = '|';
  qbeShowQueryCommand();
}


static void qbeMenuExec(void)
{
  KEYSET oldSet, newSet ;
  void *look ;
  QBEGET("qbeMenuExec");

  graphCheckEditors(qbeGraph, TRUE) ;
  /* form on query interface window */
  if (graphExists(querGraph)) {
    messStatus("Forming Query Command");
    if(!graphActivate(querGraph))
      messcrash("Cannot reactivate querGraph in qbeMenuExec.\n");
      
    queryWindowShowQuery(qbeFormQuery(resbuffer)); /* calls QUERGET inside*/
    queryCommandEdit();
    if(!graphActivate(qbeGraph))
      messcrash("Cannot reactivate qbeGraph in qbeMenuExec.");
    return;
  }

  /* Execute query */
  if(!keySetActive(&oldSet, &look))
    newKeySet("QBE Answer") ;
  if(keySetActive(&oldSet, &look)) { 
    messStatus("Searching") ;
    newSet = query(oldSet, qbeFormQuery(resbuffer)) ;
    if (newSet != oldSet )  {
      if(keySetExists(query_undo_set) &&
	 query_undo_set != oldSet )
	keySetDestroy(query_undo_set) ;
      query_undo_set = oldSet ;
    }
    keySetShow(newSet,look) ;
  }
  else
    messout("First select a keySet, thank you") ;
}


static void qbeMenuQuit(void)
{
  QBEGET("qbeMenuQuit") ; 
  if(graphActivate(qbeGraph))
    graphDestroy () ;
  qbeGraph = 0;
  qbeDestroy();
}

static void qbeMenuDraw(void)
{
  QBEGET("qbeMenuDraw") ; 
  
  qbeDisplay () ;
  qbePick(qbelook->activebox);
}


/*
 * qbeMetadata:
 * using current bs node in display, traverse paternity and
 * form metadata query of type ">search_class paternity_label"
 */

/*
static void qbeMetadata(int search_class) {
  BS temp ;
  TV tagval ;
  QBEGET("qbeMetadata");

  if (!qbelook->activebox || qbelook->activebox <= itemEntryBox)
    temp = qbelook->objet->root;
  else {
    tagval = array(qbelook->tagval, qbelook->activebox, TV) ;
    while (tagval->tag) {
      tagval = tagval->tag;
    }
    temp = tagval->bs;
  }

  metaGetData(search_class, qbelook->key, temp);
}


static void qbeMenuMetaFull(void) {
  qbeMetadata(_VElement);
}

static void qbeMenuMetaLong(void) {
  qbeMetadata(_VLongText);
}
*/

/*************** drawing package ***********************/


static QBE drawLook ;
static int xmax = 0 ;
static BOOL isModel ;

/* show query command on query command line */
static void qbeShowQueryCommand(void) {
  if (graphExists(querGraph)) 
    { Graph old = graphActive() ;
      graphActivate(querGraph) ;
      queryWindowShowQuery(qbeFormQuery(resbuffer)); /*calls QUERGET inside */

      graphActivate(old) ;
    }
}  

static BOOL notNumeric (char *text, char *c)
{
  char *t = text;
  while (*t) {
    if ((*t > '9' || *t < '0') && *t != '.' && *t != '-') {
      *c = *t;
      return TRUE ;
    }
    t++;
  }
  *c = 0 ;
  return FALSE ;
}

static BOOL itemEntryHandler(char *text, int unused)
{
  Graph old = graphActive() ;
  QBEGET("itemEntryHandler");
  if (strcmp(text, ""))
    item_enable = TRUE ;
  else
    item_enable = FALSE;
  /*   qbeDisplay () ; mhmp 26.01.99 */
  /*    qbePick(qbelook->activebox); mhmp 16.09*/
  /*  qbeShowQueryCommand(); mhmp 26.01.99*/
  graphActivate(old) ;
  return TRUE ;
}


static BOOL  qbeEntryHandler (char *text, int unused)
{
  TV tagval ;
  char c ;
  Graph old = graphActive() ;
  QBEGET("qbeEntryHandler");
  
  tagval = array(qbelook->tagval, nameBox /* qbelook->activebox */ , TV) ;
  if (tagval->tv_type == TV_NUMERIC_ENTRY && notNumeric(text, &c)) {
    messout ("Warning: Non-numeric character %c in numeric entry for tag (%s).",
	     c, tagval->tag ? tagval->tag->name : "no tag");
  }
  /*  qbeDisplay () ;
       qbePick(qbelook->activebox);*/
  qbePick (nameBox) ; /* mhmp 17.09.98 */
  /*  qbeShowQueryCommand(); mhmp 26.01.99*/
  graphActivate (old) ;
  return TRUE ;
}

extern void bsSubFuseType(OBJ obj, BS bs) ;  /* gha debug 6/15/93 */

/*
 * This recursive routine traces thru BS structure and draws the 
 * tree display for a class model, with the value fields and boolean
 * tag-attributes having entry boxes.
 */
static float qbeDrawBS (BS bs, int x, float y, int prev_x, float prev_y)
{
  float yMe = 0 ;
  int xPlus = 0 ;
  char *text = 0 ;
  int box;
  int oldstyle, newstyle;
  int length = 0 ;
  TV tagval = 0 ;
  BOOL new_tag ;
  char timeBuf[25] ;

  yMe = y ;
  new_tag = FALSE ;
  
  if (isModel)
    text =  name (bs->key) ; /* IS A MODEL */
  if (bs->key <= _LastC)
    text = bsText(bs) ; 
  else if (bs->key == _Int)  /* gha here */
    text = messprintf ("%d", bs->n.i) ; 
  else if (bs->key == _Float)
    text = messprintf ("%f", bs->n.f) ; 
  else if (bs->key == _DateType)
    text = bs->n.time ?  timeShow (bs->n.time, timeBuf, 25) : "" ;
  else
    text = name (bs->key) ;   /* gha: need this...tag names */
  
  if (iskey (bs->key) == 2 && class(bs->key) != _VText)  {
    newstyle = BOLD;
    oldstyle = graphTextFormat(newstyle) ;
  }
  else {
    newstyle = PLAIN_FORMAT;
    oldstyle = graphTextFormat(newstyle);
  }
  
  if (text)  { 
    if (bsIsTag(bs))  
      {  /* is a tag  see systags.wrm */
	if (!class(bs->key) && (bs->key) >= 50 && !xref_tag) 
	  {  
	    box = graphBoxStart() ;
	    if (box > boxMax) 
	      {
		boxMax = box;
		new_tag = TRUE;
	      }
	    
	    tagval = array(drawLook->tagval, box, TV) ;
	    /* mhmp 16.09*/
	    if (!tagval ||/* !tagval->tag ||*/ !tagval->tag_defined) {
	      tagval = (TV) halloc(sizeof(struct TAGVALUESTUFF),graphHandle());
	      new_tag = TRUE;
	    }
	    
	    if (new_tag) {
	      tagval->box = box;
	      strcpy(tagval->name, text);
	      tagval->textformat = -1;  /* like a NULL for TextFormat */
	      tagval->tag_active = FALSE;
	      tagval->tag_defined = TRUE;
	      tagval->tv_type = TV_TAG;
	      tagval->tag = 0 ;
	      tagval->has_value = 0;
	      tagval->value_number = 0;
	      tagval->bs = bs;
	      arr(drawLook->tagval, box, TV) = tagval;
	    }
	    
	    lasttag = tagval;
	    lasttag_idx = 0;
	    
	    graphText (tagval->name, x, yMe) ;
	    
	    graphBoxEnd() ;
	    length = strlen(tagval->name);  
	    if (tagval->tag_active == TRUE)
	      graphBoxDraw(tagval->box, SELECT_FCOL, SELECT_BCOL);
	    
	  }  /* end of is tag within range */
	else
	  {
	    if (bs->key == _XREF)
	      xref_tag = TRUE;
	    else if (xref_tag == TRUE)
	      xref_tag = FALSE;
	    length = 0;  /* system and XREF tags ignored */
	  }
      }
    else {   /* is not a tag, i.e. is a value field  */

      /* gha debug 6/15/93:  eventually should add # to QBE here */
      /* if (class(bs->key) && bs->key == KEYMAKE(class(bs->key), 1))
       *	bsSubFuseType(drawLook->objet, bs);
       * else { 
       */

      /* do not show #subtypes at all */
      if (!class(bs->key) 
	  ||
	  (class(bs->key) && bs->key != KEYMAKE(class(bs->key), 1))) 
	{
	  length = 12;
	  box = graphBoxStart() ;
	  tagval = array(drawLook->tagval, box, TV); 
	  if (!tagval || !tagval->tag || !tagval->tag_defined) 
	    {
	      tagval = (TV) halloc(sizeof(struct TAGVALUESTUFF),graphHandle());
	      
	      tagval->box = box;
	      strcpy(tagval->name, "");
	      tagval->textformat = newstyle;
	      tagval->tag_active = FALSE;
	      tagval->tag_defined = TRUE;
	      if (bs->key == _Int || bs->key == _Float)        /* NUMERIC VALUE */
		tagval->tv_type = TV_NUMERIC_ENTRY;
	      else if (bs->key == _DateType)        /* DATE VALUE */
		tagval->tv_type = TV_DATE_ENTRY;
	      else
		tagval->tv_type = TV_TEXT_ENTRY;    /* else assume TEXT */
	      
	      /* tagval->tag = lasttag->tag ;mhmp 17.09.98*/
	      tagval->tag = lasttag;
	      tagval->has_value = 0 ;
	      tagval->value_number = ++lasttag_idx ;
	      tagval->bs = bs;
	      array(drawLook->tagval, box, TV) = tagval;
	    }
	  else if (tagval->tag) 
	    {
	      if (strcmp(tagval->name,"")) 
		tagval->tag->has_value = tagval->value_number;
	      else if (lasttag_idx == 1)
		tagval->tag->has_value = 0;
	    }

	if (x + length > xmax)
	  xmax = x + length ;
      
	/* mieg -> editor
	box = graphTextScrollEntry (tagval->name, BUFFER_SIZE -1,
				    length, x, yMe, qbeEntryHandler) ;
	*/
	nameBox =  graphTextScrollEditor ("",tagval->name,BUFFER_SIZE -1, length, x, yMe, qbeEntryHandler) ;
	box = nameBox ;
	graphAssociate (tagval->name, assVoid(box)) ;

	array(drawLook->tagval, box-1, TV) = tagval;
	array(drawLook->tagval, box, TV) = tagval; 
	array(drawLook->tagval, box+1, TV) = tagval;
	if (box+1 > boxMax)  /* add 1 for cursor */
	  boxMax = box+1;
	
	drawLook->classe = class(bs->key) ; 
	graphBoxEnd() ;
      }   /* end of class stuff */
    } /* end of else (is not a tag) */
    
    yMe += 1.2; /*mhmp 27.10.98   1.5 --> 1.2*/
    
    if (length > xPlus)
      xPlus = length;
  }
  
  graphTextFormat(oldstyle) ;
  
  xPlus += x;
  if (xPlus > xmax)
    xmax = xPlus ;
  
  xPlus = xPlus + 6 - xPlus%4 ;   /* "Tabify" the next position to the right */

  if (bs->right) /* to the right at same y */
    y = qbeDrawBS (bs->right, xPlus, y, xPlus, y) ; 

  if (yMe > y) {
    y = yMe ;
    yMe -= 1.2; /* mhmp 27.10.98 1.5 --> 1.2 */
  }

  if (bs->down)      /* go to previous x position */
    y = qbeDrawBS (bs->down, x, y, xPlus, yMe) ;

  return y ;
}

#define BUTTON_OFFSET 4.3  /* should be at least 4 */
#define CHOICE_OFFSET 2  /* should be about 2 */
#define ITEM_OFFSET   2.5    /* should be about 2 */
static void qbeDisplay (void)
{
  int ymax = 1 ;
  float x1=0, y1=0, x2=0, y2=0;
  int old;
  int fromBox ;
  QBEGET("qbeDisplay");
  
  graphClear () ;
  
  graphButtons(QBE_Buttons, 1, 0.5, 62) ;

  graphText("From:", 1, BUTTON_OFFSET + 0.2);
  old = graphTextFormat(BOLD);
  fromBox = graphButton(currentFromName(), 
			 fromButtonResult, 
			 8, BUTTON_OFFSET);
  graphTextFormat(old);
  graphBoxFreeMenu(fromBox,(FreeMenuFunction) fromMenuResult, fromMenu);
  graphBoxDim(fromBox, &x1, &y1, &x2, &y2);
  graphRectangle(x1-0.4, y1-0.2, x2+0.4, y2+0.3);
  
  graphText("Class:", 18, BUTTON_OFFSET + 0.2);
  qbe_classBox = graphBoxStart();
  graphRectangle(25-0.4, BUTTON_OFFSET-0.1, 47+0.4, BUTTON_OFFSET+1.5);
  graphBoxEnd();
  graphBoxFreeMenu(qbe_classBox, (FreeMenuFunction)classesMenuResult,
		   classesMenu);
  strcpy(qbe_class_buffer, currentClassName());
  qbe_classEntryBox = graphTextScrollEditor("",qbe_class_buffer, BUFFER_SIZE - 1,
					   20, 25, BUTTON_OFFSET + 0.2,
					   qbe_classEntryHandler);
  qbe_classTriBox = graphDownSelTriangle(TRUE, 45.8, BUTTON_OFFSET + 0.2);

  isModel = (KEYKEY(qbelook->key) == 0) ;

  if (qbelook->objet->root && qbelook->objet->root->right)  { 
    drawLook = qbelook ;

    itemBox = graphBoxStart() ;
    graphText("ITEM NAME :", 1, BUTTON_OFFSET + CHOICE_OFFSET );
    graphBoxEnd();

    itemEntryBox = 
      graphTextScrollEditor("",item_buffer,
			   BUFFER_SIZE-1, 30, 13,
			   BUTTON_OFFSET + CHOICE_OFFSET, 
			   itemEntryHandler);
    if (itemEntryBox > boxMax)
      boxMax = itemEntryBox ;

    graphLine(1, BUTTON_OFFSET + CHOICE_OFFSET + 1.3,
	      59, BUTTON_OFFSET + CHOICE_OFFSET + 1.3);

    if (item_enable)
      graphBoxDraw(itemBox, SELECT_FCOL, SELECT_BCOL);
    else
      graphBoxDraw(itemBox, NORMAL_FCOL, NORMAL_BCOL);

    if (!drawLook->text)
      drawLook->text = (char *) halloc (256, graphHandle()) ;
    else
      memset(drawLook->text, 0, (mysize_t) 256) ;
    
    xmax = 0 ;
    ymax = (int) qbeDrawBS (qbelook->objet->root->right, 
			    2, BUTTON_OFFSET + CHOICE_OFFSET + ITEM_OFFSET, 
			    2, BUTTON_OFFSET + CHOICE_OFFSET + ITEM_OFFSET ) ;
  }
  
  graphTextBounds (xmax, ymax) ;
  graphRedraw() ;
  
  if (qbelook->activebox >= boxMax)
    qbelook->activebox = 0 ;
  
  isModel = FALSE ;
}

/************************************************************/

static void qbePick (int box)
{
  TV tagval=0;
  int oldstyle;
  int i;
  KEY key;

  QBEGET("qbePick") ;

  if (box <= 0) {
      graphBoxDraw(qbelook->activebox, NORMAL_FCOL, NORMAL_BCOL);
      qbelook->activebox = 0;
    return ;
  }
  else
    if (box > boxMax) { 
      messerror("qbePick received a box number, %d, beyond range.\n",box);
      return ;
    }
  
  qbelook->activebox = box;
  
  if (qbelook->activebox == qbe_classEntryBox)
    /*    boxMax = 0 ;  mhmp 17.09.98 */
    /* graphTextScrollEntry(qbe_class_buffer, 0, 0, 0, 0, 0) */  ;
  else if (qbelook->activebox == qbe_classTriBox) 
    {
      /* graphTextScrollEntry(qbe_class_buffer, 0, 0, 0, 0, 0); */
      if (graphSelect(&key, classesMenu))
	classesMenuResult(key, 0);
    }
  else if (qbelook->activebox == qbe_classBox)
    {
      /*      boxMax = 0 ;  mhmp 17.09.98 */
      return;
    }
  else if (qbelook->activebox == itemBox) 
    {
      if (item_enable) {
	graphBoxDraw(qbelook->activebox, NORMAL_FCOL, NORMAL_BCOL);
      item_enable = FALSE;
      }
      else {
	graphBoxDraw(qbelook->activebox, SELECT_FCOL, SELECT_BCOL);
	item_enable = TRUE;
      }
    }
  else if (qbelook->activebox == itemEntryBox)
    /* graphTextScrollEntry(item_buffer, 0, 0, 0, 0, 0)  */ ;
  else 
    {
      /* every text entry box is 3 boxes, outer, shading, cursor */
      for (i = qbelook->activebox; i > qbelook->activebox - 3; i--)  
	{
	  tagval = arr(qbelook->tagval, i, TV);
	  if (tagval) break;
	}
      if (tagval && tagval->tag_defined)  
	{
	  if (tagval->tv_type == TV_TEXT_ENTRY ||
	      tagval->tv_type == TV_NUMERIC_ENTRY ||
	      tagval->tv_type == TV_DATE_ENTRY
	      ) 
	    {         /* is a value */
	      oldstyle = graphTextFormat(tagval->textformat);
	      /* mieg->editor 	graphTextScrollEntry(tagval->name, 0, 0, 0, 0, 0); */
	      graphTextFormat(oldstyle);
	      
	      if (tagval->tag && strcmp(tagval->name, "")) {
		graphBoxDraw(tagval->tag->box, SELECT_FCOL, SELECT_BCOL);
		tagval->tag->tag_active = TRUE;
	      }
	    }
	  else {                   /* is a tag */
	    if (tagval->tag_active) {
	      graphBoxDraw(qbelook->activebox, NORMAL_FCOL, NORMAL_BCOL);
	      tagval->tag_active = FALSE;
	    }
	    else {
	      graphBoxDraw(qbelook->activebox, SELECT_FCOL, SELECT_BCOL);
	      tagval->tag_active = TRUE;
	    }	
	    qbeShowQueryCommand();
	  }
	}
      else
	messerror("Bad Item Picked--not a tag or entry box %d\n", 
		  qbelook->activebox);
    }
}


/*********************************************************************/
/********************  public routines   *****************************/

static void qbeCreate2 (void)
{
  QBE newQbeLook;
  KEYSET oldSet = 0 ;
  void *look;
  int i;
  
  if (graphActivate(qbeGraph)) { 
      graphPop() ;
      return ;
    }
  
  qbeGraph = displayCreate (DtQueryByExample); 
  
  if (!qbeGraph)
    return;
  
  newQbeLook=(QBE)messalloc(sizeof(struct QBESTUFF));
  
  newQbeLook->magic = QBE_MAGIC;
  
  if (graphExists(qbuildGraph) && (qbuildGraph != 0)
      && (pickType(KEYMAKE(qbuild_selected_class, 0)) == 'B')) 
    qbe_class = qbuild_selected_class;
  else  
    {
      keySetActive(&oldSet,&look);
      if (!oldSet ||
	  (!(qbe_class = class(keySet(oldSet, 0))) 
	   || pickType(KEYMAKE(qbe_class, 0)) != 'B')) 
	{
	  for (i = 1; i <= qbeClassMax(); i++) 
	    {
	      qbe_class = pickWord2Class(CLASSLIST_NAME(i));
	      if ( pickType(KEYMAKE(qbe_class, 0)) == 'B' &&
		  iskey(KEYMAKE(qbe_class, 0)) == 2) 
		break ;
	    }
	}
    }

  newQbeLook->key = KEYMAKE(qbe_class, 0);

  if (pickType (newQbeLook->key) != 'B' 
      || 
      !(newQbeLook->objet = bsCreate(newQbeLook->key))) {

    if(graphActivate(qbeGraph))
      graphDestroy () ;
    messout("Bad model for class %s.\n Removing QBE Window.", 
	    pickClass2Word(qbe_class));
    return;
  }
  
  newQbeLook->graph = qbeGraph;
  
  /* init classlist */
  qbeClassMax() ;

  graphRetitle (messprintf("Query By Example:   %s", 
			   className(newQbeLook->key), name(newQbeLook->key)));
  graphRegister (DESTROY, qbeDestroy) ;
  graphRegister (PICK,(GraphFunc) qbePick) ;
  graphMenu (QBE_Menu) ;
  graphHelp ("Query_By_Example") ;
  
  graphAssociate (&QBE_MAGIC, newQbeLook) ;
  newQbeLook->activebox = 0 ;
  newQbeLook->tagval = arrayCreate(1, TV);
  array(newQbeLook->tagval, 0, TV) = NULL;
  
  classesMenuSet();
  lasttag = NULL;
  lasttag_idx = 0;
  xref_tag = FALSE;
  boxMax = 0;
  item_enable = FALSE;
}

void qbeCreate (void)
{ qbeCreate2 () ;
  from_keyset = FALSE ;
  qbeDisplay () ;
}

void qbeCreateFromKeySet (KEYSET ksAlpha, int i0)
{ KEY *kp ;
  int i = keySetMax (ksAlpha) - i0 ;

  from_keyset = TRUE ;
  if (i0 >= 0 && i > 0)
    { kp = arrp (ksAlpha, i0, KEY) - 1 ;
      while (kp++, i--)
	if (pickType (*kp) == 'B')
	  { qbe_class = class (*kp) ;
	    break ;
	  }
    }
  qbeCreate2 () ;
  from_keyset = TRUE ;
  qbeDisplay () ;
}

/*****  end of qbedisp.c *****/


 
