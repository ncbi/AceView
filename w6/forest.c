/*  File: forest.c
 *  Author: Michel Potdevin et Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: View multiple trees in a user configurable way
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 17 03:17 1998 (rd)
 * Created: Sept 20 1996 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */

#define DEFINE_OBJ
typedef struct sobj *OBJ ;
#include "acedb.h"

#include "bitset.h"
#include "bs_.h"
#include "biblio.h"
#include "map.h"
#include "client.h"
#include "bindex.h"

/************************************************************/

#define graphBoxBox(_box) { \
               float _x1, _y1, _x2, _y2 ; \
               graphBoxDim (_box, &_x1, &_y1, &_x2, &_y2) ; \
               graphRectangle (_x1 - 0.4, _y1 - 0.1, _x2 + 0.4, _y2 + 0.1) ; \
                }

#define LINE_LENGTH 108
#define MODELFLAG      0x0002
#define MODELNEXTFLAG  0x0004
#define MODELNEXTTAG   0x0008
#define DISCARDFLAG    0x0010

/************************************************************/

		/* addresses of next two used as unique ids */
static int GRAPH2FORESTLOOK_ASSOC ; /* used to find the forest LOOK 
				       from the active graph */
static int FORESTLOOK_MAGIC ;	/* used to verify the forest LOOK */

#define FORESTGET(name) FOREST forest ; \
	if (!graphAssFind (&GRAPH2FORESTLOOK_ASSOC, &forest)) \
          messcrash ("forest-graph not found in %s",name) ; \
        if (!forest) \
          messcrash ("%s received a null forest-pointer",name) ; \
        if (forest->magic != &FORESTLOOK_MAGIC) \
          messcrash ("%s received a non-magic forest-pointer",name)


/************************************************************/
/*  keep destroy and definition alike */
#define SEGMAGIC 6542987
typedef struct segStruct
{
  KEY key ;
  int color, flag ;
  int n ;
  BOOL yellow ;
  int magic ;
} SEG ;

/****************************************************************/

typedef struct ForestStruct
{ 
  void *magic ;			/* == &FORESTLOOK_MAGIC */
  char mot[1024], *title, cq [20] ;
  KEY model ; int classe ;
  Graph graph ;
  MAP map ;

  BOOL showModel, hide, isDownRed, isPageDownRed, 
                  isUpRed, isPageUpRed, isSelect ;

  int  
   firstBox, isGoto,
   max, maxDiscard, nBottom, nTop, nUp,nBegin, nLast, nSelected, 
   lineBegin, lineBeginOld ; 
/*mhmp 27.11.97 lineBeginOld: 
  dessin identique pour un middleclick trop bas a la fin */
  Stack s ;
  KEYSET fSet, alphaSet, selectSet, greenTags, yellowTags, selectedTags,
         tagState, modelTags ;
  BitSet  discard ;
  Array box2segs, segs, lineNumber, size, topLines ;
  Array topLinesIndex, unDiscard ;
  BitSet show, unSelect, verified ;
} *FOREST ;

/****************************************************************/

static int gotoBox, lightBox,upBox, downBox, closeAllBox,  
           pageUpBox, pageDownBox, discardBox ;

static BOOL isRecentre = FALSE ;
static  int  halfHeight ;
static double yOld ;

static void forestUp (void) ;
static void forestDown (void) ;
static void forestPageUp (void) ;
static void forestPageDown (void) ;
static void forestUnverify (FOREST forest) ;
static BOOL forestMayDisplay (FOREST forest, int n) ;
static void forestAddTags (FOREST forest, int nn)  ;
static void forestKbd (int k) ;
static void forestSelect (void) ;
static void forestColor (FOREST  forest) ;

/************************************************************/

static void forestDestroy(void)
{ 
  FORESTGET ("forestDestroy") ;

  keySetDestroy (forest->fSet) ;
  keySetDestroy (forest->alphaSet) ;
  keySetDestroy (forest->selectSet) ;
  keySetDestroy (forest->greenTags) ;
  keySetDestroy (forest->yellowTags) ;
  arrayDestroy (forest->topLines) ;
  arrayDestroy (forest->topLinesIndex) ;/*mhmp 03.12.97*/
  stackDestroy (forest->s) ;
  arrayDestroy (forest->segs) ;
  arrayDestroy (forest->box2segs) ;
  arrayDestroy (forest->modelTags) ;
  arrayDestroy (forest->tagState) ;
  bitSetDestroy (forest->verified) ;
  bitSetDestroy (forest->discard) ;
  bitSetDestroy (forest->show) ;
  bitSetDestroy (forest->unSelect) ;
  arrayDestroy (forest->unDiscard) ;
  arrayDestroy (forest->lineNumber) ;
  arrayDestroy (forest->size) ;
  forest->magic = 0 ;
  messfree (forest->title) ;
  mapDestroy (forest->map) ;
  messfree (forest) ;

  return;
} /* forestDestroy */

/***********************************************************************/

static int fState(FOREST forest, KEY key)
{
  if (key > 0 && key < arrayMax(forest->tagState))
    return array (forest->tagState, key, char) ;
  return 0 ;
}

/***********************************************************************/
/********************** actions  ***************************************/
static void forestDoRedraw (FOREST forest) ;
static void forestDoRedraw2 (FOREST forest) ;
static void forestDraw (FOREST forest) ;
static void forestRedraw (void) ;

/*******************/

static void forestFollow (FOREST forest, int box)
{ 
  int i ;
  SEG *seg;
  BOOL isDisp = FALSE ;

  if (box < 1 ||  box >= arrayMax(forest->box2segs))
    return ;

  seg = arr(forest->box2segs,box,SEG*) ;
  seg->magic = SEGMAGIC ;
  if (seg->n && stackExists (forest->s) && seg->n < stackMark (forest->s))
    {
      graphPostBuffer (stackText(forest->s, seg->n)) ;
      if (seg->yellow)
	{
	  isDisp = TRUE ; /* mhmp 15.06.98 */
	  if (pickType(seg->key) == 'B')
	    display (seg->key, 0, TREE) ;
	  else
	    display (seg->key, 0, 0) ;
	}
      else 
	for (i = forest->firstBox ; i < arrayMax(forest->box2segs) ; i++)
	  {
	    SEG *seg1 = arr(forest->box2segs,i,SEG*) ;
	    
	    if (seg1 && !(seg1->flag & MODELFLAG))
	      {
		if (seg1->yellow)
		  graphBoxDraw (i, BLACK, WHITE) ;
		seg1->yellow = FALSE ;
		if (seg1->n < stackMark (forest->s) &&
		    !lexstrcmp(stackText(forest->s, seg1->n),
			       stackText(forest->s, seg->n)))
		  {
		    graphBoxDraw (i, BLACK, YELLOW) ;
		    seg1->yellow = TRUE ;
		  }
	      }
	  }
    }
  else if (seg->key == forest->model)
    {
      display (seg->key, 0, TREE) ;
      isDisp = TRUE ; /* mhmp 05.05.98 */
    }
  if (!isDisp) 
    {
      /* mhmp 28.05.98 necessaire apres un mot "trop long" (winPos)*/
      memset(forest->mot, 0, 1024) ;
      strncpy (forest->mot,"", 1023) ;
      graphTextEntry (forest->mot, 0, 0, 0, 0) ;
      graphTextEntry (forest->mot, 0, 0, 0, 0) ;
      /**********/
      memset(forest->mot, 0, 1024) ;
      strncpy (forest->mot, stackText(forest->s, seg->n), 1023) ;
      graphTextEntry (forest->mot, 0, 0, 0, 0) ;
      graphEntryDisable() ; 
      graphRegister (KEYBOARD, forestKbd) ;    
      /* mhmp 07.01.98 */
      forestColor(forest) ;
    }
}

static void forestColor (FOREST  forest)
{
  int i ;
  SEG *seg ;
  char buf [1025] ;
  
  if (forest && arrayExists(forest->box2segs))
    {
      if (forest && *forest->mot)
	{
	  memset(buf, 0, 1025) ;
	  i = strlen (forest->mot) ;
	  if (i > 1000) 
	    invokeDebugger() ;
	  strcpy(buf, "*") ;
	  strcat (buf, forest->mot) ;
	  strcat (buf, "*") ;
	}
      for (i = forest->firstBox, seg = arr(forest->box2segs,i,SEG*) ; 
	   i < arrayMax(forest->box2segs) ; seg++, i++) 
	{
	  seg = arr(forest->box2segs,i,SEG*) ;
	  if (!seg || seg->magic != SEGMAGIC || !seg->n)
	    continue ;
	  seg->yellow = FALSE ;
	  if (*buf && 
	      (!(seg->flag & MODELFLAG)) &&
	      seg->n > 0 &&
	      seg->n < stackMark (forest->s) &&
	      pickMatch (stackText(forest->s, seg->n), buf)
	      )
	    { graphBoxDraw (i, BLACK, YELLOW) ; seg->yellow = TRUE ; }
	}
    }
  return;
} /* forestColor */

/*******************/


static void callBiblio (void)
{ 
  KEYSET ks ;
  char *cp ;
  FORESTGET ("callBiblio") ;

  cp = "From forest" ;
  if (messPrompt("Title","","t"))
    cp = freeword() ;     
  ks = keySetCopy (forest->selectSet) ;
  keySetSort (ks) ;
  keySetCompress(ks) ;

  biblioKeySet (cp, ks) ;
  keySetDestroy (ks) ;

  return;
} /* callBiblio */


static void forestPick (int box, double x, double y) 
{ 
  KEY tag ;
  SEG *seg, *seg1 ;
  FORESTGET ("forestPick") ;
  
  if (box == gotoBox)
    { graphTextEntry (forest->cq, 0, 0, 0, 0) ;
      forest->isGoto = TRUE ;
    }
  else  if (box == lightBox)
    graphTextEntry (forest->mot, 0, 0, 0, 0) ;
  else if (box == upBox) 
    forestUp () ;
  else if (box == downBox) 
    forestDown () ;
  else if (box == pageUpBox) 
    forestPageUp () ;
  else if (box == pageDownBox) 
    forestPageDown () ;
  else if (box == closeAllBox) 
    {
      int ii, jj = 1 ;
      
      for (ii = 0, seg1 = arrayp (forest->segs, ii, SEG) ;
	   ii < arrayMax (forest->segs) ; ii++, seg1++)
	if (seg1 && (seg1->flag & MODELFLAG))
	  if (array (forest->tagState, seg1->key, char) == 1)
	    { jj = 4 ; break ; }
      for (ii = 0, seg1 = arrayp (forest->segs, ii, SEG) ;
	   ii < arrayMax (forest->segs) ; ii++, seg1++)
	if (seg1 && (seg1->flag & MODELFLAG))
	  array (forest->tagState, seg1->key, char) = jj ; 
      
      forestUnverify(forest) ;
      forestDoRedraw (forest) ;
    }
  else  if (box < 1 ||                /* check valid box */
      box >= arrayMax(forest->box2segs) ||
      !(seg = arr(forest->box2segs, box, SEG*))
      )
    return ;
  else if (seg->flag & MODELFLAG)
    { tag = seg->key ;

      switch (fState (forest, tag))
	{
	case 0:
	  break ;
	case 1:
	  if (seg->flag & MODELNEXTFLAG)
	    array (forest->tagState, tag, char) = 2 ; 
	  else if (seg->flag & MODELNEXTTAG)
	    array (forest->tagState, tag, char) = 3 ;
	  else
	    array (forest->tagState, tag, char) = 4 ;
	  break ;
	case 2: case 3:
	  array (forest->tagState, tag, char) = 4 ;
	  break ;
	case 4: 
	  array (forest->tagState, tag, char) = 1 ;
	  break ;
	}
      forestUnverify(forest) ;
      forestDoRedraw (forest) ;
    }
  else if (seg->flag == DISCARDFLAG)
    { bitSet (forest->discard, seg->key) ;
      array (forest->unDiscard, forest->maxDiscard, int) = seg->key ;
      forest->maxDiscard ++ ;
   /* mhmp 27.11.97 pour eviter un ecran vide en discardant le dernier */
      if (forest->isDownRed && !forest->isUpRed)
	{
	  forest->lineBeginOld = 0 ; /* mhmp 17.12 */
	  forestUp () ;
	}
      else
	forestDoRedraw (forest) ;
    }
  else if (box > forest->firstBox) 
    forestFollow (forest, box) ;

  return;
} /* forestPick */

/***********************************************************************/

static void forestMiddleDragFast  (double x, double y)
{ 
  FORESTGET ("forestMiddleDragFast") ;

  graphXorLine (0, yOld, forest->map->thumb.x, yOld) ;
  yOld = y ;
  graphXorLine (0, yOld, forest->map->thumb.x, yOld) ;

  return;
} /* forestMiddleDragFast */


static void forestMiddleUpFast  (double x, double y)
{ 
  int n ;
  int nn = 0 ;
  KEY selectKey ;
  FORESTGET ("forestMiddleUpFast") ;

  n = WHOLE2MAP (forest->map, y) ;
  if (n < 0) n = 0 ;
  if (n >= forest->nSelected)
    n =  forest->nSelected -1 ; 

 /* mhmp 17.12.97 pb avec keyset melimelo Taich */
  selectKey = keySet (forest->selectSet, n) ;
  arrayFind (forest->alphaSet, &selectKey, &nn, keySetAlphaOrder) ;
  forestMayDisplay (forest,nn) ;
  forestDoRedraw2 (forest) ;
  /********************************
  mhmp 01.12.97 attention pour forestMayDisplay n index dans forest->alphaSet

  while (n > 0 && !forestMayDisplay (forest,n))
    n-- ;
  selectKey = keySet (forest->selectSet, n) ;
  arrayFind (forest->alphaSet, &selectKey, &n, keySetAlphaOrder) ;
  
  forest->nTop = n ;
  forestDoRedraw (forest) ;

 ***********************************/ 

  /* mhmp 02.12.97 pour avoir le "vrai" forest->nTop */
  while (n > 0)
    { 
      selectKey = keySet (forest->selectSet, n) ;
      arrayFind (forest->alphaSet, &selectKey, &nn, keySetAlphaOrder) ;
      if (forestMayDisplay (forest,nn))
	break ;
      n-- ;
    }
  if (n == 0)
    while (n < forest->nSelected)
      { 
	selectKey = keySet (forest->selectSet, n) ;
	arrayFind (forest->alphaSet, &selectKey, &nn, keySetAlphaOrder) ;
	if (forestMayDisplay (forest,nn))
	  break ;
	n++ ;
      }
  forest->nTop = nn ; 
  forestDoRedraw (forest) ;

  return;
} /* forestMiddleUpFast */


static void forestMiddleDrag  (double x, double y)
{ 
  FORESTGET ("forestMiddleDrag") ;

  graphXorLine (forest->map->thumb.x, yOld, mapGraphWidth, yOld) ;
  yOld = y ;
  graphXorLine (forest->map->thumb.x, y, mapGraphWidth, y) ;

  return;
} /* forestMiddleDrag */


static void forestMiddleUp  (double x, double y)
{
  int i, n, nn1 ;
  int exTop ;
  FORESTGET ("forestMiddleUp") ;

  if (!forest->isPageDownRed || !( y > halfGraphHeight))
    {
      n = y - halfGraphHeight ; 
      if (!arrayMax(forest->topLines)) return ;
      nn1 =  0 ; 
      /* mhmp 04.12.97 -1 a cause du decalage de 1 vers la gauche */
      for (i = 0 ; i + 1 < arrayMax(forest->topLines) ; i++)
	if (array (forest->topLines,i, int) > n)
	  break ; 
	else
	  nn1 = i ;
      forest->lineBegin =  n - array(forest->topLines,nn1, int) ;
      if (forest->lineBegin < 0)
	forest->lineBegin = 0 ;
      forest->nTop = array(forest->topLinesIndex, nn1,  int) ;
    }
  else
    forest->lineBegin = forest->lineBeginOld ;
  isRecentre = TRUE ;
  forest->showModel = FALSE ; 
  forest->lineBeginOld = forest->lineBegin ;
  exTop = forest->nTop ;
  forestDoRedraw2 (forest) ;
  /*mhmp bricolage 08.12.97 voir keyset bug */
  forest->nTop = exTop ;
  forest->lineBegin = forest->lineBeginOld ;
  isRecentre = TRUE ;
  forestDoRedraw (forest) ; /*mhmp 08.12.97*/

  return;
} /* forestMiddleUp */


static void forestMiddleDown (double x, double y)
{ 
  FORESTGET ("forestMiddleDown") ;

  if (!forest->nSelected)
    { graphRegister (MIDDLE_DRAG, 0) ;
      graphRegister (MIDDLE_UP, 0) ;
    }
  else if (x > forest->map->thumb.x)
    { 
      yOld = y ;
      graphXorLine (forest->map->thumb.x, y, mapGraphWidth, y) ;
      graphRegister (MIDDLE_DRAG, forestMiddleDrag) ;
      graphRegister (MIDDLE_UP, forestMiddleUp) ;
    }
  else
    { 
      yOld = y ;
      graphXorLine (forest->map->thumb.x, y, mapGraphWidth, y) ;
      graphRegister (MIDDLE_DRAG, forestMiddleDragFast) ;
      graphRegister (MIDDLE_UP, forestMiddleUpFast) ;
    }

  return;
} /* forestMiddleDown */
 
/***********************************************************************/

static void forestPrint(void)
{ 
  FORESTGET ("forestPrint") ;

  if (!forest->alphaSet || !keySetMax(forest->alphaSet))
    { messout ("Nothing to print, this forest is empty, sorry") ;
      return ;
    }
  forestDraw(forest);

  return;
} /* forestPrint */

/********************** menu actions  **********************************/

static void forestUnDiscard (void)
{
  int i ;
  FORESTGET ("forestUnDiscard") ;
  
  if (forest->maxDiscard)
    {
      i = array (forest->unDiscard, forest->maxDiscard - 1, int) ;
      bitUnSet (forest->discard, i) ;
      forest->maxDiscard -- ;
/*mhmp commencer par le dernier undiscarde 02.12.97*/
      forest->nTop = i ;
    }
  forestDoRedraw (forest) ;

  return;
}

/********************** menu actions  **********************************/

static void forestSaveAs (void)
{ 
  KEY key ;
  char *cp ;
  KEYSET ks ;
  int ii,j, nn ;
  FORESTGET ("forestSaveAs") ;

  if (!forest->nSelected)
    messout ("Nothing to save") ;
  else 
    { 
      if (!checkWriteAccess())
	return ;

      if (!messPrompt("Save the current set in the KeySet class\n"
		      "under the name:","","t")) 
	return ;
      cp = freeword() ;
      if (!*cp)
	return ;
      if (!lexaddkey (cp,&key,_VKeySet) &&
	  !messQuery ("This keyset already exists, "
		      "should I overwrite it ?"))
	  return ;
      if (!key)
	return ;
      
      if (!lexlock(key))
	{ messout ("Sorry, this key is locked elsewhere") ;
	  return ;
	}
      ks = keySetCreate () ;
      ii = keySetMax(forest->alphaSet) ;
      nn = 0 ; j = -1 ;
      while (j++,ii--)
	if (!bit (forest->discard, j) && !bit(forest->unSelect, j))
	  keySet(ks, nn++) = keySet(forest->alphaSet,j) ;
      keySetSort (ks) ;
      keySetCompress(ks) ;
      arrayStore (key,ks,"k") ;
      lexunlock (key) ;

      keySetNewDisplay (ks,name(key)) ;
    }

  return;
} /* forestSaveAs */

/***********************************************************************/
static void forestShowModel (void)
{
  FORESTGET ("forestShowModel") ;
  forest->showModel = !forest->showModel ;
  forestDoRedraw (forest) ;
}

/***********************************************************************/

static void forestDoHide (BOOL hh)
{
  int i ;
  char *cp ;
  FORESTGET ("forestDoHide") ;

  if ((forest->hide = hh))
    {
      i = arrayMax(forest->tagState) ;
      cp = arrp(forest->tagState, 0, char) - 1 ;
      while(cp++, i--)
	if (*cp) *cp = 4 ;
    }
  else
    {
       i = arrayMax(forest->tagState) ;
      cp = arrp(forest->tagState, 0, char) - 1 ;
      while(cp++, i--)
	if (*cp) *cp = 1 ;
    }
  forestUnverify(forest) ;
  forestDoRedraw (forest) ;

  return;
} /* forestDoHide */

static void forestTagsHide (void)
{ forestDoHide(TRUE) ; }

static void forestTagsUnHide (void)
{ forestDoHide(FALSE) ; }

/***********************************************************************/

static void forestUp (void)
{ 
  int ii ;  
  FORESTGET ("forestUp") ;

  if (forest->isUpRed) return ;
  ii = forest->nTop ;
  while (ii--)
    if (forestMayDisplay (forest,ii))
      {
	forest->nTop = ii ;
	break ;
      }
  forestDoRedraw(forest) ;
} /* forestUp */

/***********************************************************************/

static void forestDown (void)
{ 
  int ii ;
  FORESTGET ("forestDown") ;
  
  if (forest->isDownRed) return ;
  ii = forest->nTop ;
  while (++ii < forest->max)
    if (forestMayDisplay (forest,ii))
      { forest->nTop = ii ; break ; }  
  forestDoRedraw(forest) ;
}

/***********************************************************************/

static void forestGoto (char *txt)
{ 
  char *cc ;
  KEY key ;
  int i1, i2, i3, jmax, cmp ;
  FORESTGET ("forestGoto") ;
 
  jmax = keySetMax (forest->alphaSet) ;

  /* mhmp 17.12.97 solution la moins mauvaise */
  /* lexstrcmp pour tomber sur ceux qui existent Bond7 Bond68 Bond007*/
    i1 = 0 ;
  i3 = jmax -1 ;
  
  while (TRUE)
    {
      i2 = (i1 + i3) / 2 ;
      key = keySet (forest->alphaSet, i2) ;
      cc = name (key) ;
      cmp = lexstrcmp (cc, forest->cq) ;
      if (!cmp)
	goto found ;
      else
      if (i2 == i1 || i2 == i3)
	break ;
      else
      if (cmp < 0)
	i1 = i2 ;
      else
	i3 = i2 ;

    }
  /* strcasecmp pour B03 B04 B041 B054  */
  i1 = 0 ;
  i3 = jmax -1 ;

  while (TRUE)
    {
      i2 = (i1 + i3) / 2 ;
      key = keySet (forest->alphaSet, i2) ;
      cc = name (key) ;
      cmp = strcasecmp (cc, forest->cq) ;
      if (!cmp)
	goto found ;
      else
      if (i2 == i1 || i2 == i3)
	break ;
      else
      if (cmp < 0)
	i1 = i2 ;
      else
	i3 = i2 ;

    }

    /* strncasecmp  pour passer au suivant */
    i2++ ;
    i3 = forest->nTop + 1 ;
    while (i3 < jmax && !forestMayDisplay (forest,i3))
      i3++ ;    
    if (i3 < jmax)  
      {
	key = keySet (forest->alphaSet, i3) ;
	cc = name (key) ; 
	if (!strncasecmp(cc,txt,strlen(txt)))
	  i2 = i3 ;
      }
found:
    forestMayDisplay (forest,i2) ;
    forestDoRedraw2 (forest) ; /*mhmp 18.12  pb avec keyset melimelo Taich*/
    while (i2> 0 && !forestMayDisplay (forest,i2))
      i2-- ;
    forest->nTop = i2 ;
    forest->isGoto = 2 ;  /* survives one redraw */
    forestDoRedraw (forest) ;
} /* forestGoto */


static void forestLight (char *txt)
{ 
  FORESTGET ("forestLight") ;

  forestDoRedraw (forest) ;

  return;
} /* forestLight */

/***********************************************************************/

static void forestPageUp (void)
{
  FORESTGET ("forestPageUp") ;

  if (forest->isPageUpRed)
    return ;

  forestMiddleUp (mapGraphWidth, topMargin + 2) ;/* mhmp 04.12.97*/

  return;
} /* forestPageUp */

/***********************************************************************/

static void forestPageDown (void)
{
  FORESTGET ("forestPageDown") ;
  
  if (forest->isPageDownRed)
    return ;

  forestMiddleUp (mapGraphWidth, mapGraphHeight) ;
  
  return;
} /* forestPageDown */

/***********************************************************************/
static void forestKbd (int k) /*mh*/
{
  FORESTGET ("forestKbd") ;

  switch (k)
    {
    case UP_KEY :
      forestUp () ;
      break ;
    case DOWN_KEY :
      forestDown () ;
      break ;
    case PAGE_UP_KEY :
    case LEFT_KEY :
      forestPageUp () ;
      break ;
    case PAGE_DOWN_KEY :
    case RIGHT_KEY :
      forestPageDown () ;
      break ;
    default:
      return ;
    }
} /* forestKbd */

/***********************************************************************/
/*******************  query/select system ******************************/

/***********************************************************************/
/* i order the green tags by depth, this will optimize bIndexSearch
  by hitting the top tags first
  */
static void forestTagsSort (int classe, KEYSET ks)
{
  KEY *kp, tag ;
  KEYSET path ;
  int i, n = keySetMax(ks), p ;

  if (n<2) return ;
  i = n ;
  kp = arrp(ks, 0, KEY) - 1 ;
  while (kp++, i--)
    {
      tag = *kp ;
      path = bsGetPath (classe, tag) ;
      p = path ? keySetMax(path) : 0 ; /* p is the length of the path to the tag */
      if (p > 250) p = 250 ;
      *kp = KEYMAKE (p, tag) ;
      keySetDestroy (path) ;
    }
  keySetSort (ks) ; /* sorts on class then tag value */
  i = n ;
  kp = arrp(ks, 0, KEY) - 1 ;
  while (kp++, i--)
    *kp = KEYKEY (*kp) ;

  return;
} /* forestTagsSort */

/***********************************************************************/

static void forestGetYellowTags (FOREST forest)
{
  KEY *kp ;
  int i , ii = 0, classe ;
  KEYSET ks = keySetCreate();
  KEYSET tmp;
  
  keySetDestroy (forest->yellowTags) ;
  classe = class(keySet(forest->alphaSet, 0)) ; /* should be current model */
  tmp = bsTagsInClass (classe) ;
  if (!tmp)
    return ;
  ii = 0 ; i = keySetMax(tmp) ; kp = arrp(tmp, i - 1, KEY) + 1 ;  
  while (kp--, i--)
    switch  (fState(forest,*kp))
      {
      case 1: case 2:
	keySet(ks, ii++) = *kp ;
	break ;
      default:
	break ;
      }
  forestTagsSort(classe, ks) ;
  forest->yellowTags = ks ;
  keySetDestroy (tmp) ;

  return;
} /* forestGetYellowTags */

/***********************************************************************/

static BOOL forestQueryKey (KEYSET tags, KEY key)
{
  KEY *gtp ;
  int i ;
  OBJ obj = 0 ;
  int found = 1 ;
  
  if (!key || pickType(key) != 'B')
    return 2 ;

  if (bIndexVersion (-1))
    {
      found = 0 ;
      i = keySetMax(tags) ;
      gtp = arrp(tags, 0, KEY) - 1 ;
      while (gtp++, found < 2 && i--)
	{
	  switch (bIndexFind (key, *gtp))
	    {
	    case 0: /* absent */
	      break ;
	    case 1:
	      found = 1 ;
	      break ;
	    case 2:  /* found */
	      found = 2 ;
	      return TRUE ;
	    }
	}
    }
  
  if (!found)
    return 0 ;

  if ((obj = bsCreate (key)))
    {
      i = keySetMax(tags) ;
      gtp = arrp(tags, 0, KEY) - 1 ;
      found = 0 ;
      while (gtp++, found < 2 && i--)
	if (bsFindTag (obj, *gtp))
	  { bsDestroy (obj) ; return TRUE ; }
      bsDestroy (obj) ;
    }
  return 0 ;
} /* forestQueryKey */


static void forestUnverify (FOREST forest) 
{
  int ii ;
 
  ii = keySetMax(forest->alphaSet) ; 
  while (ii--)
    bitUnSet(forest->verified, ii) ;

  return;
} /* forestUnverify */


static BOOL forestQueryN (FOREST forest, int n)
{
  if (!keySetMax(forest->yellowTags))
    { forestAddTags (forest, n) ;
      return TRUE ;
    }

  /*    if (!bit(forest->verified,n)) mhmp 05.12.97*/
    {
      if (forestQueryKey (forest->yellowTags, keySet(forest->alphaSet,n)))
	bitSet (forest->show,n) ;
      else
	bitUnSet (forest->show,n) ;
      forestAddTags (forest, n) ;
    }
  bitSet (forest->verified, n) ;

  return bit (forest->show,n) ;
} /* forestQueryN */


static BOOL forestMayDisplay (FOREST forest, int n)
{
  BOOL isOK;

  if (n < 0 || n >= keySetMax(forest->alphaSet))
    messcrash ("Bad n in forestMayDisplay") ;

  isOK = !bit (forest->unSelect, n) &&
    !bit (forest->discard, n) &&
    forestQueryN (forest, n);

  return isOK;
} /* forestMayDisplay */

/***********************************************************************/
/* looping */
static int forestDoSelect (FOREST forest, BOOL doSelect)
{
  int ii, j, nn ;
  KEY *kp ;

  ii = keySetMax(forest->alphaSet) ; 
  while (ii--)
    bitUnSet(forest->unSelect, ii) ;
  if (!doSelect) 
    {
      keySetDestroy (forest->selectSet) ;
      forest->selectSet = keySetCopy (forest->alphaSet) ;
    }
  else
    {
      forest->selectSet = keySetReCreate(forest->selectSet) ;
      keySetSort (forest->greenTags) ;
      keySetCompress (forest->greenTags) ;
      forestTagsSort(forest->classe, forest->greenTags) ;
      ii = keySetMax(forest->alphaSet) ; nn = 0 ; j = -1 ;
      kp = arrp(forest->alphaSet, 0, KEY) - 1 ;
      while (kp++, j++, ii--)
	{
	  if (!forestQueryKey (forest->greenTags, *kp))
	    bitSet(forest->unSelect, j) ;     /*mhmp ii --> j 28.11.97*/
	  else  /* This is it ! */
	    keySet(forest->selectSet, nn++) = keySet(forest->alphaSet,j) ;
	  if (messIsInterruptCalled())
	    break ;
	}
    }
  keySetDestroy (forest->selectedTags) ; 
  forest->selectedTags = keySetCopy (forest->greenTags) ;
  forest->nSelected = keySetMax(forest->selectSet) ; 
  forest->map->min = 0 ;
  forest->map->max = forest->nSelected ;

  return forest->nSelected ;
} /* forestDoSelect */

/***********************************************************************/
/*  graphic interface */

static void forestSelect (void)
{
  FORESTGET ("forestSelect") ;
  
  if (!forest->showModel)
   {
     forestDoRedraw (forest) ; /* pour garder le bouton Select bleu mhmp*/
     return ;
   }
  forest->isSelect = !forest->isSelect ;
  if (forest->isSelect)
    { 
      if (!keySetMax(forest->greenTags))
	{ 
	  messout ("Please first select some tags from the model") ;
	  forest->isSelect = FALSE ;
	  forestDoSelect (forest, FALSE) ;
	}
      else
	{
	  messout("%s%s%s\n%s",
		  "With the Select mode, only the objects with at ",
		  "least one green (or yellow) tag are selected. ",
		  "Objects with only red tags are hidden. ",
		  "F4 to interrupt" ) ;
	  forestDoSelect (forest, TRUE) ;
	}
    }
  else
    forestDoSelect (forest, FALSE) ;
  /*  graphUnMessage () ; mhmp 26.11.97*/

  forestDoRedraw(forest) ;

  return;
} /* forestSelect */

/***********************************************************************/

static void forestExport (void) /* mhmp 20.02.98 */
{
  char *cp ;
  KEYSET ks ;
  int ii,j, nn ;
  FORESTGET ("forestExport") ;

  ks = keySetCreate () ;
  ii = keySetMax(forest->alphaSet) ;
  nn = 0 ; j = -1 ;
  while (j++,ii--)
    if (!bit (forest->discard, j) && !bit(forest->unSelect, j))
      keySet(ks, nn++) = keySet(forest->alphaSet,j) ;
  keySetSort (ks) ;
  keySetCompress(ks) ;
  if (keySetMax(ks) == 0)
    messout ("The active forest does not contain keys") ;
  else
    {
      cp = "From forest" ;
      if (messPrompt("Title","","t"))
	{
	  cp = freeword() ;
	  keySetNewDisplay (ks,cp) ;
	}
    }

  return;
} /* forestExport */

/***********************************************************************/

static MENUOPT forestMenu[]={
  { graphDestroy,"Quit"},
  { help, "Help"},
  { forestPrint, "Print"},
  { displayPreserve, "Preserve"},
  { 0,0 } };

/***********************************************************************/

static MENUOPT forestMenuButt1[]={
  { graphDestroy,"Quit"},
  { help, "Help"},
  { forestPrint, "Print"},
  { forestUnDiscard, "UnDiscard"},  /* attention a discardBox = 4th button */
  { forestExport, "Export"},
  { forestSaveAs, "Save as"},
  { 0,0 } };

/***********************************************************************/

static MENUOPT forestMenuButt2[]={
  { forestShowModel, "Model show/hide"},
  { forestTagsHide, "Tags hide"},
  { forestTagsUnHide, "Tags show"},
  { forestSelect, "Select"},
  { 0,0 } };

/***********************************************************************/
/* drawing */
/***********************************************************************/
static float oldOffset ;

static void forestDisplayMark (LOOK look, float *offset)
     /* used as a MapColDrawFunc for mapInsertCol() */
{
  float oldTextHeight ;
  char *cc, dd[8], ee[8] ;
  int i, j ;
  float y, yOld = 3 ;
  BOOL idem ;
  KEY key ;
  FOREST forest = (FOREST)look;
  
  if (!forest->nSelected)
    return ;

  memset (ee, 0, sizeof(ee));
  memset (dd, 0, sizeof(dd));
  oldOffset = *offset ;
  *offset += 1 ;                       /*mhmp -6 ?????????*/
  forest->map->thumb.fac =
    (mapGraphHeight - topMargin - 6) / (forest->map->max);

  oldTextHeight = graphTextHeight (0.75) ;
  idem = TRUE ;
  for (j = 0 ; j < forest->nSelected ; j+= 1)  
    { 
      key = keySet(forest->selectSet, j) ;
      cc = name(key) ;
      memset(dd,0,8) ;
      for (i = 0 ; i < 5 ; i+= 1)
	{
	  dd [i] = cc [i] ;
	  if (dd[i] != ee [i])
	    {
	      idem = FALSE ;
	      ee[i] = dd[i] ;
	    }
	}
      if (idem)
	continue ;
      else
	{ idem = TRUE ;
	  y = MAP2WHOLE (forest->map, j) ;
	  if (y > mapGraphHeight)
	    break ;
	  if (y - yOld < 2)
	    continue ;
	  else
	    {
	      yOld = y ;
	      graphLine (0, y, *offset + 2, y) ;
	      graphText (messprintf("%s", dd), 0, y + 0.2) ;
	    }
	}
    }
  graphTextHeight (oldTextHeight) ;
  *offset += 2 ;

  return;
} /* forestDisplayMark */

/***********************************************************************/

static void forestAddTags (FOREST forest, int nn) 
{
  KEY key, *kp1 ;
  int i ;
  KEYSET ks1 ;

  key = keySet (forest->alphaSet, nn) ;
  if (pickType(key) != 'B')
    return;
  ks1 = bsKeySet (key) ;
  if (!ks1) return ;
  i = keySetMax(ks1) ;
  kp1= arrp(ks1, 0, KEY) - 1 ;
  while
    (kp1++, i--)
      if (!class(*kp1) && *kp1 >= _Date && /* a real tag */
	  !array(forest->tagState, *kp1, char))
	array(forest->tagState, *kp1, char) = 1 ; 
  keySetDestroy (ks1) ;

  return;
} /* forestAddTags */

/***********************************************************************/

static char* forestBStext (BS bs)
{
  static char timeBuf[25] ;

  if (bs->key <= _LastC)
    return bsText(bs) ;
  else
    if (bs->key ==_Int)
      return messprintf ("%d", bs->n.i) ;
    else if (bs->key ==_Float)
      return messprintf ("%g", bs->n.f) ;
    else
      if (bs->key ==_DateType)
	return timeShow (bs->n.time, timeBuf, 25) ;
      else
	return name (bs->key);
} /* forestBStext */



static int countSiblings (BS bs)
{
  int n = 0 ;

  while (bs)
    { if (class(bs->key) != _VUserSession)
        n++ ;
      bs = bs->down ;
    }
  return n ;
} /* countSiblings */



static BOOL drawTriangle(BS bs, int x, int y, BOOL isModel)
{
  BOOL isArrow ;
  int n = countSiblings(bs->right) ;

  if (n > 0) 
    {
      graphColor(BLUE) ;
      if (isModel)
	graphText(messprintf("  -----> COUNT"), x, y) ;
      else
	graphText(messprintf("  -----> %d ",n), x, y) ;
      graphColor(BLACK) ;
      isArrow = TRUE ;
    }
  else
    isArrow = FALSE ;
  return isArrow ;
} /* drawTriangle */

static FOREST myForest = 0 ;
static BOOL showTimeStamps = FALSE ; /* could go in menu choice */
static int topLine = 0 ;
/* recursive routine */

/* forestBS (graphic) forestBBS (printer) */

static int forestBS (BS bs, BS bsm, int x, int y, int *yy, 
		     BOOL isRoot, BOOL isModel)	
{
  int yMe = y , y3 ;
  int yyMe = *yy ;
  int xPlus = 0 ;
  BOOL isArrow = FALSE ;
  char *text, *cp ;
  int box;
  int oldTextFormat ;
  SEG *seg ;
  KEY tag ;

  graphFitBounds (&mapGraphWidth,&mapGraphHeight) ;

  if (y > mapGraphHeight + myForest->lineBegin)
    return y ;

  if (isModel && isRoot)
    {
      box = closeAllBox = graphBoxStart() ;
      graphText ("Close All", x, y++) ;
      graphBoxEnd() ;
      graphBoxDraw (box, BLACK, RED) ;
    }

  switch (fState(myForest, bs->key))
    {
    case 2:
      isArrow = TRUE ;
      break ;
    case 3: 
      if (!isModel)
	goto right ;
    case 4:
      if (!isModel)
	goto down ;
    };
  
  bsModelMatch (bs, &bsm) ; /* will fail on Time stamp, ignore returned value */

  if (!showTimeStamps && class(bs->key) == _VUserSession)
    goto right;
  
  text = forestBStext (bs) ;
  if (!text) text = "" ; 

  yMe = y ; /* must do this here again: y can change */
  yyMe = *yy ;
  seg = arrayp (myForest->segs, arrayMax(myForest->segs), SEG) ;
  seg->magic = SEGMAGIC ;
  seg->key = bs->key ;
  seg->color = WHITE ;
  seg->flag = 0 ;
  seg->n = stackMark (myForest->s) ;
  pushText (myForest->s, text) ;
  if (isModel)
    { seg->color = BLUE ;
      tag = seg->key ;
      if (!keySetFind (myForest->modelTags, tag, 0)) 
	keySetInsert (myForest->modelTags, tag) ;
      if (tag < arrayMax(myForest->tagState) &&
	  array (myForest->tagState,tag, char))
	{ 
	  seg->flag = MODELFLAG ;
	  if (bs->right && bs->right->key != _UNIQUE) 
	    seg->flag |= bsIsTag (bs->right) ? MODELNEXTTAG : MODELNEXTFLAG ;
	  switch (fState(myForest, tag))
	    {
	    case 1:
	      seg->color = GREEN ;
	      /* mhmp 10.07.98 */
	      if (myForest->isSelect && myForest->selectedTags &&  
		  keySetFind (myForest->selectedTags, tag, 0))
		seg->color = LIGHTBLUE ;
	      if (isModel)
		keySet(myForest->greenTags, keySetMax(myForest->greenTags)) = tag ;
	      break ;
	    case 2:
	      seg->color = YELLOW ;
	      break ;
	    case 3:
	      seg->color = PALEORANGE ;
	      break ;
	    case 4:
	      seg->color = LIGHTRED ;
	      break ;
	    }
	}
    }
  if (seg->color == BLUE)
    goto right ;

  { 
    oldTextFormat = PLAIN_FORMAT; /* init for compiler happiness */

    box = graphBoxStart() ;
    if ((iskey (bs->key) == 2 &&
	 class(bs->key) != _VText)) 
      oldTextFormat = graphTextFormat(BOLD) ;

    if (bs->key == _Greek)
      oldTextFormat = graphTextFormat(GREEK) ;

    uLinesText (text,LINE_LENGTH) ;
    if ((cp = uNextLine(text)))          
      { 
	y3 = yMe++ - myForest->lineBegin ;
	if (y3 >= topLine)
	  graphText (cp,x,y3) ;	/* write out this node */
	if (strlen(cp) > xPlus)
	  xPlus = strlen(cp) ;
	while ((cp = uNextLine(text)))	/* indent following lines */
	  { 
	    y3 = yMe++ - myForest->lineBegin ;
	    if (y3 >= topLine )
	      graphText (cp,x+2, y3) ;
	    if (strlen(cp)+2 >xPlus)
	      xPlus = strlen(cp)+2 ;
	  }
      }
    if ( bs->right && !isRoot && isArrow)
      /*    if (isDecal)
	{ yyMe++ ;
	  if (*yy >= myForest->lineBegin)
	    isDecal = FALSE ;
	}
      else  mhmp 05.05.98 */
      {
	y3 =  yMe-1- myForest->lineBegin ;
	if  (y3 >= topLine)
	  isArrow = drawTriangle (bs, x+xPlus, y3, isModel) ;
      }
    /*isArrow = drawTriangle (bs, x+xPlus, yMe-1, isModel) ;*/
    if ((iskey (bs->key) ==2 && class(bs->key) != _VText) ||
 	( (bs->size & ON_FLAG) && bs->right ) ||
 	bs->key == _Greek )
      graphTextFormat(oldTextFormat) ;
    graphBoxEnd() ;
    array(myForest->box2segs, box, SEG*) = seg ;

    if (isModel) graphBoxDraw(box, BLACK, seg->color) ;
  }
 right:
  xPlus += x ;
  xPlus = xPlus + 6 - xPlus%4 ;
  if (isRoot)
    { xPlus = x + 2 ;
      y += 1 ;
      *yy +=1 ;
    }
  if (bs->right &&  fState(myForest, bs->key) != 4 && bs->right->key !=_XREF &&  ( !isArrow || isRoot))
    y = forestBS (bs->right, bsm ? bsModelRight(bsm) : 0, xPlus, y, yy, 
		  FALSE, isModel) ; /* to the right at same y */

  if (yMe > y) y  = yMe ;
  if (yyMe > *yy) *yy = yyMe ;
 down:
  if (bs->down) /* && !(bs->size & ATTACH_FLAG)) attach only rightwards */
    y = forestBS (bs->down, bsm, x, y, yy, FALSE, isModel) ;		/* below at new y location */

  return y ;
} /* forestBS */


static int forestBBS (BS bs, BS bsm, int x, int y, 
		      BOOL isRoot) /* recursive routine */
{
  int yMe = y ;
  int xPlus = 0 ;
  BOOL isArrow = FALSE ;
  char *text, *cp ;
  int box;
  int oldTextFormat ;
  SEG *seg ;
  
  switch (fState(myForest, bs->key))
    {
    case 2:
      isArrow = TRUE ;
      break ;
    case 3: 
      goto right ;
    case 4:
      goto down ;
    }
  
  bsModelMatch (bs, &bsm) ; /* will fail on Time stamp, ignore returned value */

  if (!showTimeStamps && class(bs->key) == _VUserSession)
    goto right;
  
  text = forestBStext (bs) ;

  yMe = y ; /* must do this here again: y can change */

  seg = arrayp (myForest->segs, arrayMax(myForest->segs), SEG) ;
  seg->magic = SEGMAGIC ;
  seg->key = bs->key ;
  seg->color = WHITE ;
  seg->flag = 0 ;

  {
    oldTextFormat = PLAIN_FORMAT; /* init for compiler happiness */

    box = graphBoxStart() ;
    if ((iskey (bs->key) == 2 &&
	 class(bs->key) != _VText)) 
      oldTextFormat = graphTextFormat(BOLD) ;

    if (bs->key == _Greek)
      oldTextFormat = graphTextFormat(GREEK) ;

    uLinesText (text,LINE_LENGTH) ;
    if ((cp = uNextLine(text)))
      { graphText (cp,x,yMe++) ;	/* write out this node */
	if (strlen(cp) > xPlus)
	  xPlus = strlen(cp) ;
	while ((cp = uNextLine(text)))	/* indent following lines */
	  { graphText (cp,x+2,yMe++) ;
	    if (strlen(cp)+2 >xPlus)
	      xPlus = strlen(cp)+2 ;
	  }
      }
    if ( bs->right && !isRoot && isArrow)
      isArrow = drawTriangle (bs, x+xPlus, yMe-1, FALSE) ;
    if ((iskey (bs->key) ==2 && class(bs->key) != _VText) ||
 	( (bs->size & ON_FLAG) && bs->right ) ||
 	bs->key == _Greek )
      graphTextFormat(oldTextFormat) ;
    graphBoxEnd() ;
    array(myForest->box2segs, box, SEG*) = seg ;
  }
 right: 
  xPlus += x ;
  xPlus = xPlus + 6 - xPlus%4 ;
  if (isRoot)
    { xPlus = x + 2 ;
      y += 1 ;
    }
  if (bs->right &&  fState(myForest, bs->key) != 4 &&  (!isArrow || isRoot))
    y = forestBBS (bs->right, bsm ? bsModelRight(bsm) : 0, xPlus, y, 
		   FALSE) ; /* to the right at same y */

  if (yMe > y) y = yMe ;
 down:
  if (bs->down) /* && !(bs->size & ATTACH_FLAG)) attach only rightwards */
    y = forestBBS (bs->down, bsm, x, y, FALSE) ;		/* below at new y location */

  return y ;
} /* forestBBS */



static void forestKey (FOREST forest, int nn, int *xp, int *yp)
{
  KEY key ;
  OBJ obj = 0 ;
  int box ;
  SEG *seg ;

  key = keySet (forest->alphaSet, nn) ;

  seg = arrayp (forest->segs, arrayMax(forest->segs), SEG) ;
  seg->magic = SEGMAGIC ;
  array(forest->box2segs, box = graphBoxStart(), SEG*) = seg ;
  seg->key = key ;
  seg->color = WHITE ;
  graphBoxEnd () ;

  if (pickType(key) != 'B' || ! (obj = bsCreate(key)))
    return ;
  myForest = forest ;
  *yp += 1 ;
  *yp = forestBBS (obj->root, bsModelRoot(obj), *xp, *yp, TRUE) ;
  bsDestroy (obj) ;

  return;
} /* forestKey */



static void forestDrawModel (FOREST forest, int *xp, int *yp)
{
 KEY key ;
  OBJ obj = 0 ;
  int box, old, yy ;
  SEG *seg ;

  forest->greenTags = keySetReCreate(forest->greenTags) ;
  key = forest->model ;

  seg = arrayp (forest->segs, arrayMax(forest->segs), SEG) ;
  seg->magic = SEGMAGIC ;
  seg->key = key ;
  seg->color = WHITE ;

  array(forest->box2segs, box = graphBoxStart(), SEG*) = seg ;
  old = graphTextFormat (BOLD) ;
  graphText (messprintf("%s", name(key)+ 1), *xp, *yp) ;
  graphTextFormat (old) ;
  graphBoxEnd () ;

  graphText ("Click tag to: show->green, count->yellow, hide->red",*xp + 10, *yp) ;

  if ((obj = bsCreate(key))) /* garanteed in fact */
    { 
      myForest = forest ;
      yy = 1 ;
      *yp = forestBS (obj->root, bsModelRoot(obj), *xp, *yp, &yy, TRUE, TRUE) ;
      bsDestroy (obj) ;
    }

  return;
} /* forestDrawModel */

/****************************************************/

static void forestDrawKeyUp (FOREST forest, int nn, int *xp, int *yp)
{ 
  KEY key ;
  OBJ obj = 0 ;
  int yy ;

  if (!forestMayDisplay(forest, nn))
    return ;

  key = keySet (forest->alphaSet, nn) ;
 
  array(forest->topLines,arrayMax(forest->topLines), int) =
    *yp - forest->lineBegin ;
  array(forest->topLinesIndex,arrayMax(forest->topLinesIndex), int) = nn ;  
  *yp += 1 ;
  if ((obj = bsCreate(key)))
    { 
      myForest = forest ;
      yy = 1 ;
      *yp = forestBS (obj->root, bsModelRoot(obj), *xp, *yp, &yy, TRUE, FALSE) ;
      bsDestroy (obj);
    }

  return;
} /* forestDrawKeyUp */

/****************************************************/

static void forestDrawKey (FOREST forest,
			   int nn, int *xp, int *yp, BOOL doDraw)
{ 
  KEY key ;
  OBJ obj = 0 ;
  int box, yy, yDis ;
  float hh ;
  SEG *seg ;

  if (!forestMayDisplay(forest, nn))
    return ;

  key = keySet (forest->alphaSet, nn) ;
  
  if (doDraw)
    {	
      array(forest->topLines, arrayMax(forest->topLines), int) =
	*yp -  forest->lineBegin ;
      array(forest->topLinesIndex, arrayMax(forest->topLinesIndex), int) = nn ;
      seg = arrayp (forest->segs, arrayMax(forest->segs), SEG) ;
      seg->magic = SEGMAGIC ;
      seg->key = key ;
      seg->color = WHITE ;
      array(forest->box2segs, box = graphBoxStart(), SEG*) = seg ;
      graphBoxEnd () ;
      
      seg = arrayp (forest->segs, arrayMax(forest->segs), SEG) ;
      seg->magic = SEGMAGIC ;
      array(forest->box2segs, box = graphBoxStart(), SEG*) = seg ;
      seg->key = nn ;
      seg->flag = DISCARDFLAG ;
      hh = graphTextHeight (0.9) ;
      /* mhmp 26.11.97*/
      yDis = *yp -  forest->lineBegin > topMargin + 1 ? 
	*yp + 1 -  forest->lineBegin : topMargin + 2 ;
	graphText ("Discard", 63, yDis - .7) ;

      graphBoxBox (box) ;
      graphBoxEnd () ;
      graphTextHeight (hh) ;
    }
      
   *yp += 1 ;
   if ((obj = bsCreate(key)))
     { 
       myForest = forest ;
       yy = 1 ;
       *yp = forestBS (obj->root, bsModelRoot(obj), *xp, *yp, &yy, TRUE, FALSE) ;
       bsDestroy (obj);
     }

   return;
} /* forestDrawKey */

/****************************************/

static void forestDraw (FOREST forest)
{
  int n = 0, x, y = 0, limit = 1250 ;
  float oldy, yy ;
  Graph old = graphActive() ;

  graphClear () ;

  graphBoxStart() ;
  graphText(forest->title, 0, y) ;
  forestGetYellowTags (forest) ;

  for (x = 0, n= 0, y = 1 ;
       n < forest->max ; n++)
    { 
      if (!forestMayDisplay (forest,n))
	continue ;

      forestKey (forest, n, &x, &y) ;
      if (y > limit)
	{ 
	  if (!messQuery("%d objects already use %d lines, "
			 "do you want to continue ?",  n, y))
	    break ;
	  limit *= 4 ;
	}
    }
  graphBoxEnd() ;

  yy = y ; /* cast to float */
  oldy = graphFakeBounds (yy) ;
  graphRedraw () ; /* michel: should come last */
  graphPrint () ;
  graphActivate (old);
  graphFakeBounds (oldy) ;
  forestRedraw();

  return;
} /* forestDraw */


static void forestDoRedraw2 (FOREST forest)
{ 
  graphFitBounds (&mapGraphWidth,&mapGraphHeight) ; 

  topMargin = 5 ;
  stackClear (forest->s) ;
  pushText (forest->s, "") ;
  graphClear () ;
  forest->segs = arrayReCreate (forest->segs, 256, SEG) ;
  forest->box2segs = arrayReCreate (forest->box2segs, 256, SEG*) ; 
  forest->firstBox = forest->nSelected = 0 ;
  mapDrawColumns (forest->map) ;

  return;
} /* forestDoRedraw2 */


static void forestDoRedraw (FOREST forest)
{ 
  graphFitBounds (&mapGraphWidth, &mapGraphHeight);

  topMargin = 5 ;
  stackClear (forest->s) ;
  pushText (forest->s, "") ;
  graphClear () ;
  forest->segs = arrayReCreate (forest->segs, 256, SEG) ;
  forest->box2segs = arrayReCreate (forest->box2segs, 256, SEG*) ;
  forest->firstBox = forest->nSelected = 0 ;
  mapDrawColumns (forest->map) ;
  graphRegister (PICK, forestPick) ; /* redo because of ColControl */
  graphMenu (forestMenu) ;
  forestDisplayMark ((LOOK)forest, &oldOffset) ; 
  graphRedraw () ;
  forestColor (forest) ; /* mhmp 24.03.98*/

  return;
} /* forestDoRedraw */



/* mhmp 04.12.97 *********************************/
/* calcul des forest->linesTop au-dessus de nTop */

static void forestDisplayTextUp (FOREST forest, float *offset)
{ 
  int nbUp ;
  int n = 0, line = 0, x, y ;
  forest->nUp = forest->nTop - 1 ;
  nbUp = 0 ;
  while ((nbUp < halfGraphHeight) && (forest->nUp >= 0))
    {
      if (forestMayDisplay (forest, forest->nUp))
	nbUp++ ;
      forest->nUp-- ;
    }
  if (forest->nUp < 0) /* mhmp 05.12.97 A voir*/
    forest->nUp = 0 ;

  /*  if (nbUp ==0)
    forest->nUp = forest->nTop ;*/
  /* on decale de -10000 */
  line = topMargin + 1 - 100000 ;
  if (forest->isSelect) line++ ;   
  y = line ;
  topLine = y + 1; /* needed to not block drawing model */
  closeAllBox = 0 ;
  if (forest->showModel)
    { x = *offset ;
    forestDrawModel (forest, &x, &y) ;
    }
  topLine = y ;
  if (keySetMax(forest->yellowTags) && forest->nUp >= 0)
    for (x = *offset, n = forest->nUp ;
	 n <= forest->nTop ; n++)
      forestDrawKeyUp (forest, n, &x, &y) ;

  return;
} /* forestDisplayTextUp */



static void forestTopLines (FOREST forest)
/*********** mhmp 04.12.97 rearrangement des topLines *****/
{ 
  int i, ii = 0 ;

/* ii: rang du 2eme nTop */
  for (i=0 ; i < arrayMax(forest->topLines) ; i++)
    if (array (forest->topLines,i, int) > -1000)
      { 
	ii = i ;
	break ;
      }

/* calcul des topLines au-dessus de nTop en egalant les 2 nTop */ 
  for (i=0 ; i < ii ; i++)
    array (forest->topLines,i, int) = array (forest->topLines,i, int) -
      array (forest->topLines,ii-1, int) + array (forest->topLines,ii, int) ;

/* decalage de 1 a gauche a partir du 1er nTop */
  if (ii > 0)
    for (i=ii-1; i + 1 < arrayMax(forest->topLines) ; i++)
      {
	array (forest->topLines,i, int) = 
	  array (forest->topLines,i+1, int) ;
	array (forest->topLinesIndex,i, int) =
	  array (forest->topLinesIndex,i+1, int) ;
      }

  return;
} /* forestTopLines */

/*******************************************/

static void forestDisplayText (LOOK look, float *offset)
     /* used as a MapColDrawFunc for mapInsertCol() */
{
  int box, n = 0, line = 0 , x, y , ii ;
  int first, last ;
  int selectKey ;
  char *cp ;
  float x1, y1, x2, y2, dx, dx2, dx4, dy, dy2, dy4, dec;
  float oldCentre, delta ;
  KEY key, key1 ;
  FOREST forest = (FOREST)look;

  oldOffset = *offset ;
  *offset += 6 ;
 
  first = last = 999999999 ;
  forestGetYellowTags(forest) ;
  forest->topLines = arrayReCreate (forest->topLines,30, int) ;
  forest->topLinesIndex = arrayReCreate (forest->topLinesIndex,30, int) ;
  if (!isRecentre)
    forest->lineBegin = 0 ;
  isRecentre = FALSE ;
  line = topMargin + 1 ;

  /*************************************/
  /* dry run, search for non empty objects, select the relevant model */
    
  y = line ;
  box = graphBoxStart() ;

  if (forest->isGoto > 0)
    forest->isGoto-- ;
  /* mhmp 02.12.97 pour avoir le bon forest->nTop */
  ii = forest->nTop ;
  while (ii < forest->max)
      {
	if (forestMayDisplay (forest,ii))
	  { 
	    forest->nTop = ii ;
	    break ;
	  }
	ii++ ;
      }

  forest->map->centre = forest->nTop + halfGraphHeight - 1 ; 
  if (forest->map->centre < halfHeight)
    forest->map->centre = halfHeight ;

  forest->isDownRed = TRUE ; forest->isPageDownRed = FALSE ;
  if (keySetMax(forest->yellowTags))
    {
      ii = forest->nTop ;
      while (++ii < forest->max)
	if (forestMayDisplay (forest,ii))
	  { forest->isDownRed = FALSE ; break ; }
    }

  key = keySet(forest->alphaSet, forest->nTop) ;
  if (key)
    {
      lexword2key ( messprintf("?%s",className(key)), &key1, _VModel) ;
      forest->model = key1 ;  
      forest->classe = class (key) ;
    }

  for (x = *offset, n = forest->nTop ;
       y < mapGraphHeight + forest->lineBegin &&  n < forest->max ; n++)
    { 
      if (!forestMayDisplay (forest,n))
	continue ;
      if (first == 999999999)
	first = n ;
      last = n ; 
      if (keySetMax(forest->yellowTags))
	forestDrawKey (forest, n, &x, &y, FALSE) ;
      else
	y+= 2 ;
    }
  forest->nBottom = n ;
  graphBoxEnd() ;
  graphBoxClear (box) ;
    /*************************************/
  /* calcul des forest->topLines  au-dessus de forest->nTop mhmp 04.12.97 */

  forestDisplayTextUp (forest, &oldOffset) ;

  /*************************************/
  /* real drawing */

  box = graphBoxStart() ;
  forest->firstBox = box ;
  line = topMargin + 1 ;
  if (forest->isSelect) line++ ; 
  if (forest->lineBegin != 0) line++ ; /*mhmp 01.12.97*/
  y = line ;
  topLine = y + 1; /* needed to not block drawing model */
  closeAllBox = 0 ;
  if (forest->showModel)
    { x = *offset ;
    forestDrawModel (forest, &x, &y) ;
    }
  topLine = y ;
  if (keySetMax(forest->yellowTags))
    { 
      for (x = *offset, n = forest->nTop ;
	   y < mapGraphHeight + forest->lineBegin &&  n < forest->max ; n++)
      {
	forestDrawKey (forest, n, &x, &y,TRUE) ;
      } 
/********** rearrangement des forest->linesTop *********/

      forestTopLines (forest) ;

/*******************************************/
      /*mhmp 27.11.97*/
      if (/*forest->isDownRed &&*/ y < mapGraphHeight + forest->lineBegin)
	forest->isPageDownRed = TRUE ;
    }
  graphBoxEnd() ; 

 /*************************************/

  /* reports */
 /*mhmp 28.11.97 forest->max-1  --> forest->max */
  cp = messprintf("%d items, %d discarded",
		  forest->max, forest->maxDiscard) ;
  graphText (cp,*offset, topMargin) ;
  ii = strlen(cp) + 2 ;
  if (forest->isSelect)
    graphText(messprintf(", %d selected", forest->nSelected), 
	      *offset + ii, topMargin) ;
  
  /* buttons */
  box = graphButtons (forestMenuButt1, *offset, 0.5, 1000) ;
  discardBox = box + 3 ;  /* 4th box */
  if (forest->maxDiscard)
    graphBoxDraw (discardBox, BLACK, LIGHTBLUE) ;

  box = graphButtons (forestMenuButt2, *offset, 2.2, 1000) ;
  if (forest->showModel)
    graphBoxDraw (box, BLACK, LIGHTBLUE) ;
  if (forest->hide)
    graphBoxDraw (box + 1, BLACK, LIGHTBLUE) ;  
  else
    graphBoxDraw (box + 2, BLACK, LIGHTBLUE) ;  
  if (forest->isSelect)
    graphBoxDraw (box + 3, BLACK, LIGHTBLUE) ;
  if (biblioPossible (forest->fSet))
    graphButton ("Biblio", callBiblio, 51.0, 0.5) ;

  /* mouvements up/down*/
  /*  x1 = 50 ;*/
  /*  y1 = 4 ;*/
  x1 = 52 ;
  y1 = 2 ;
  dx = 2 ;
  dy = 1.2 ;
  dec = 0.05 ;
  dx2 = dx/2 ;
  dx4 = dx2/2 ;
  dy2 = dy/2 ;
  dy4 = dy2/2 ;
  x2 = x1 + dx ;
  y2 = y1 + dy ;

  pageDownBox  = graphBoxStart () ;
  graphRectangle (x1, y1, x2, y2) ;
  graphLine (x1 + dx4, y1 + dec, x1 + dx2, y1 + dy2 + dec) ;
  graphLine (x1 + dx2, y1 + dy2 + dec, x2 - dx4, y1 + dec) ;
  graphLine (x1 + dx4, y1 + dy2 - dec, x1 + dx2, y2 - dec) ;
  graphLine (x1 + dx2, y2 - dec, x2 - dx4, y1 + dy2 - dec) ;
  graphBoxEnd() ;
  graphBoxDraw (pageDownBox, BLACK, GREEN) ;
  downBox  = graphBoxStart () ;
  x1 = x2 + dx2 ;
  x2 = x1 + dx ;
  graphRectangle (x1, y1, x2, y2) ;
  graphLine (x1 + dx4, y1 + dy4, x1 + dx2, y2 - dy4) ;
  graphLine (x1 + dx2, y2 - dy4, x2 - dx4, y1 + dy4) ;
  graphBoxEnd() ;
  graphBoxDraw (downBox, BLACK, GREEN) ;
  x1 = x2 + dx2 ;
  x2 = x1 + dx ;
  upBox  = graphBoxStart () ;
  graphRectangle (x1, y1, x2, y2) ;
  graphLine (x1 + dx4, y2 - dy4, x1 + dx2, y1 + dy4) ;
  graphLine (x1 + dx2, y1 + dy4, x2 - dx4, y2 - dy4) ;
  graphBoxEnd() ;
  graphBoxDraw (upBox, BLACK, GREEN) ;
  x1 = x2 + dx2 ;
  x2 = x1 + dx ;
  
  pageUpBox  = graphBoxStart () ;

  graphRectangle (x1, y1, x2, y2) ;
  graphLine (x1 + dx4, y1 + dy2 + dec, x1 + dx2, y1 + dec) ;
  graphLine (x1 + dx2, y1 + dec, x2 - dx4, y1 + dy2 + dec) ;
  graphLine (x1 + dx4, y2 - dec, x1 + dx2, y1 + dy2 - dec) ;
  graphLine (x1 + dx2, y1 + dy2 - dec, x2 - dx4, y2 - dec) ;

  graphBoxEnd() ;  
  graphBoxDraw (pageUpBox, BLACK, GREEN) ;
  graphText ("Light:", 34, 4) ; /* mhmp 05.05.98 20-->50*/
  lightBox = graphTextScrollEntry (forest->mot, 1023, 50, 41, 4, forestLight) ;
  graphText ("Go to:", 5, 4) ;
  gotoBox = graphTextEntry (forest->cq, 20, 12, 4, forestGoto) ;
  if (!forest->isGoto)
    { 
      graphEntryDisable() ; graphRegister (KEYBOARD, forestKbd) ;
    }
  
  forest->isUpRed = TRUE ; 
  forest->isPageUpRed = FALSE ;
  if (keySetMax(forest->yellowTags))
    {
      ii = forest->nTop ;
      while (ii--)
	if (forestMayDisplay (forest,ii))
	  { forest->isUpRed = FALSE ; break ; }
    }

  if (forest->isUpRed)
    {
      graphBoxDraw (upBox, BLACK, LIGHTRED) ;
      if (!forest->lineBegin)
	{ 
	  forest->isPageUpRed = TRUE ;
	  graphBoxDraw (pageUpBox, BLACK, LIGHTRED) ;
	}
    }

  if (forest->isDownRed)
    graphBoxDraw (downBox, BLACK, LIGHTRED) ;
  if (forest->isPageDownRed)
    graphBoxDraw (pageDownBox, BLACK, LIGHTRED) ; 

  if (forest->nLast)
    graphBoxDraw (pageDownBox, BLACK, LIGHTRED) ;
  forest->nLast = 0 ;
  *offset -= 2  ;
  if (first != 999999999)
    { selectKey = keySet (forest->alphaSet, first) ;
      arrayFind (forest->selectSet, &selectKey, &first, 
		 keySetAlphaOrder) ;  
    }
  if (last != 999999999)
    { selectKey = keySet (forest->alphaSet, last) ;
      arrayFind (forest->selectSet, &selectKey, &last, 
		 keySetAlphaOrder) ;  
    }
  oldCentre = forest->map->centre ;
  if (first == 999999999) {  forest->map->centre = 0 ; first = 999999998 ; last = 999999998 ; }
  else forest->map->centre = (first + last + 1)/2.0 ;
  delta = forest->map->max - forest->map->min ;
  if (!delta)
    delta = 1 ;
  forest->map->thumb.fac = (mapGraphHeight - 6 - topMargin) / delta ;
  forest->map->thumb.halfwidth = .5 * (last - first + 1) *
    (MAP2WHOLE(forest->map, forest->map->max + 1) 
     - MAP2WHOLE(forest->map, forest->map->min))/(1+keySetMax (forest->selectSet)) ;
  *offset -= 1 ;
  mapShowLocator ((void*)forest, offset) ;
  forest->map->centre = oldCentre ; 
  graphRegister (MIDDLE_DOWN, forestMiddleDown) ;
  /*  forest->lineBegin = 0 ;*/

  return;
} /* forestDisplayText */

static void forestRedraw (void)
{
  FORESTGET ("forestRedraw") ;

  graphFitBounds (&mapGraphWidth, &mapGraphHeight) ;

  halfHeight =  (mapGraphHeight - topMargin) / 2 ;
  forestDoRedraw (forest) ;

  return;
} /* forestRedraw */

/***********************************************************************/
/* global setup */
/***********************************************************************/

static Graph forestDisplayCreate(void)
{
  static Array a = 0 ;
  char *title[42], *help[32] ;
  Graph g ;
  float x, y, width, height ;
  int i, type ;

  type = 2 ;
  *title = "" ;
  x = 0.05 ;
  y = 0.05 ;
  width = 1.2 ;
  height = 0.92 ;
  *help = "forest" ;
 messerror (" Don't worry ! I have an implicit geometry for you") ;
 /***** A trick to prevent window superposition ***/
  if (!a)
    a = arrayCreate(20, Graph) ;
  for (i = 0 ; i < arrayMax(a) ; i++)
    if (!graphExists(arr(a,i,Graph)))
      break ;

  x = x + i/30.0 ;
  while (x > 1.0) x -= 1.0 ; /* prevent overflow HJC */

  y = y + i/30.0 ;
  while (y > 1.0) y -= 1.0 ;
  
  /* This is useful if the window manager does not
   * prompt the user, and is called by displayCreate
   */
  
  g =  graphCreate (type, *title, x, y, width, height) ;

  if (g)
    graphHelp(*help) ;
 
  array(a,i,Graph) = g ;

  return g ;
} /* forestDisplayCreate */


static BOOL forestAnticipate (KEYSET fSet)
{
  KEY key ;
  int i , j ;
  BOOL ok = TRUE ;

  if (externalServer)
    { 
      oldSetForServer =  keySetCreate () ;
      for (i = 0, j = 0 ; i < keySetMax(fSet) ; i++)
	{
	  key = keySet (fSet, i) ;
	  if (KEYKEY(key) && pickType(key) == 'B')
	    keySet(oldSetForServer,j++) = keySet(fSet,i) ;
	}
      if (j > 500 &&
	  !messQuery("You are requiring %d objects from the "
		     "remote server, this may be long, "
		     "do you wish to proceed", j))
	ok = FALSE ;
      else
	externalServer (-1, 0, 0, TRUE) ;
      keySetDestroy (oldSetForServer) ;
    }

  return ok ;
} /* forestAnticipate */


BOOL forestDisplayKeySet (char *title, KEYSET fSet, BOOL isOldGraph)
{ 
  FOREST forest = 0 ;
  int i, j, ii ;
  KEY key, *kp , *kp1;
  KEYSET alphaSet = 0 ;

  if (!title)
    title = "" ;
  if (!keySetExists (fSet))
    messcrash ("Bad call to forestDisplay") ;
  ii = 0 ; i = keySetMax(fSet) ; kp = arrayp(fSet,0,KEY) - 1 ;
  while (kp++, i--) 
    {
      key = lexAliasOf(*kp) ;
      if (key!= *kp) 
	{ *kp = key ; ii = 1 ; }
    }
  if (ii) { keySetSort(fSet) ; keySetCompress(fSet) ; }
  if (!forestAnticipate(fSet))
    goto abort ;

  kp = kp1 = arrp(fSet,0,KEY) ; 
  kp-- ; i = keySetMax(fSet) ; j = 0 ;
  while (kp++, i--)
    if (pickType (*kp) == 'B' &&
	iskey (*kp) ==  2)	
      { *kp1++ = *kp ; j++ ; }
  keySetMax(fSet) = j ;

  if (!keySetMax (fSet) || 
      !(alphaSet = keySetAlphaHeap(fSet, keySetMax(fSet))) ||
      !keySetMax (alphaSet))
    { messout ("Sorry, no associated forest") ;
      goto abort ;
    }

  if (isOldGraph && !graphAssFind (&GRAPH2FORESTLOOK_ASSOC, &forest))
    isOldGraph = FALSE ;

  if (isOldGraph) /* called by display.c, activeGraph is ok */
    {
      forestDestroy () ;	/* clears memory of the forest-LOOK 
				   on the active graph */
      myForest = 0 ;
      graphClear () ;
      graphGoto (0,0) ;
    }
  else
    {
      if (!displayCreate ("DtForest")) 
        /* goto abort ; */
	forestDisplayCreate () ;

      graphRegister (DESTROY, forestDestroy) ;
      graphRegister (PICK, forestPick) ;
      graphRegister (RESIZE, forestRedraw) ;
      graphRegister (MIDDLE_DOWN, forestMiddleDown) ;
      graphMenu (forestMenu) ;
    }

  graphRetitle (messprintf("More info: %s", title)) ;
  graphHelp ("Forest") ;

  forest = (FOREST) messalloc (sizeof(struct ForestStruct)) ;
  myForest = forest ;

  graphAssociate (&GRAPH2FORESTLOOK_ASSOC, forest) ;
  graphAssociate (&MAP2LOOK_ASSOC, forest) ; /* make sure the map routine find this look */

  forest->magic = &FORESTLOOK_MAGIC ;
  forest->graph = graphActive() ;

  forest->s = stackCreate (3000) ;

  forest->max = arrayMax(alphaSet) ;
  forest->fSet = fSet ;
  forest->alphaSet = alphaSet ;
  forest->maxDiscard = 0 ;
  forest->lineBegin = 0 ;
  forest->isSelect = FALSE ;
  forest->hide = FALSE ;
  forest->showModel = TRUE ;
  forest->nTop = 0 ; /* start on first object */
  forest->verified = bitSetCreate (arrayMax(forest->fSet),0) ;
  forest->discard = bitSetCreate (arrayMax(forest->fSet),0) ;
  forest->show = bitSetCreate (arrayMax(forest->fSet), 0) ;
  forest->unSelect = bitSetCreate (arrayMax(forest->fSet), 0) ;
  forest->unDiscard = arrayCreate (arrayMax(forest->fSet), int) ;
  forest->lineNumber = arrayCreate (forest->max, int) ;
  forest->size = arrayCreate (forest->max, int) ;
  forest->segs = arrayCreate (256, SEG) ;
  forest->box2segs = arrayCreate (256, SEG*) ;
  forest->firstBox = forest->nSelected = 0 ;
  forest->tagState = arrayCreate (arrayMax(fSet), char) ;
  forest->modelTags = keySetCreate() ;
  forest->title = strnew(title, 0) ;
  forest->map = mapCreate (forestRedraw) ;

  mapInsertCol (forest->map, 1, FALSE, "Mark", forestDisplayMark),
  mapInsertCol (forest->map, 2, TRUE,  "Text", forestDisplayText) ;

  graphRegister (MIDDLE_DOWN, forestMiddleDown) ;
  forest->map->mag = 1  ;

  /*
  graphFitBounds (&graphWidth,&graphHeight) ;
  halfHeight =  (graphHeight - 2) / 2 ;
  forest->map->centre = halfHeight ;
  */
  ii = keySetMax(forest->alphaSet) ; 
  while (ii--)
    bitUnSet(forest->discard, ii) ;
  forestDoSelect (forest, FALSE) ;
  mapNoMag() ; /* no mag with forest */
  forestDoRedraw2 (forest) ;
  forestDoRedraw2 (forest) ;
  forestDoRedraw (forest) ; /* probleme avec pickModel */
  return TRUE ;

abort:
  keySetDestroy (fSet) ; /* manque dans biblio.c */ 
  keySetDestroy (alphaSet) ; 
  return FALSE ;
} /* forestDisplayKeySet */

/***********************************************************************/

/* Standard DisplayFunc prototype for acedb display */
/* return true means a forest graph was constructed */
BOOL forestDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  KEYSET ks = 0 ;

  if (iskey(key) != 2) return FALSE ;
  if (class(key) == _VKeySet)
    { /* get ks from key */
    }
  
  if (!ks) ks = keySetCreate () ;
  keySet(ks,0) = key ;
  return forestDisplayKeySet (messprintf("Forest %s", name(key)), ks, isOldGraph) ;
} /* forestDisplay */

/***********************************************************************/
/*************************** eof ***************************************/
