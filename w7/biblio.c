/*  File: biblio.c
 *  Author: Richard Durbin
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: Displays any bibliographical information
 *               for current object.
 * Exported functions:
 *              biblioKey()
 *              biblioKeySet()
 * HISTORY:
 * Last edited: Dec 22 11:37 1998 (fw)
 * Created: Tue Dec 15 16:55:20 1998 (fw)
 *-------------------------------------------------------------------
 */

/* $Id: biblio.c,v 1.13 2014/05/31 01:31:02 mieg Exp $ */

#include "acedb.h"

#include <ctype.h>
#include "freeout.h"
#include "lex.h"
#include "a.h"
#include "sysclass.h"
#include "classes.h"
#include "systags.h"
#include "tags.h"
#include "query.h"
#include "bs.h"
#include "biblio.h"
#include "client.h"
#include "pick.h"

#ifndef NON_GRAPHIC
#include "graph.h"
#include "display.h"

/************************************************************/

static magic_t GRAPH2BIBLIO_ASSOC = "GRAPH2BIBLIO";
static magic_t BIBLIO_MAGIC = "BIBLIO_MAGIC";

typedef struct BiblioStruct {
  void *magic;
	char mot[1024] ;
  KEYSET bibSet ;
  Stack s ;
  char *cp ; 
  int zero, curr, width, numberRef, nbDiscard, format;
	int topLine, line2, bottomLine, flip, flop, graphWidth, graphHeight ;
  int discardUp, nbLight, nbLightBis, fini, search, dep, searchDep ;
	int firstBox, mmax ;
	int motBox, searchBox ;
  BOOL showAbstracts, showNewNumbers ;
  Array discard, abstract, abstExists, newNumber, sortNumber, author, date, ref, unDiscard ;
} *BIBLIO ;

static void biblioReDraw (void) ;
static void biblioDoReDraw (BIBLIO biblio) ;
static void biblioDraw80 (BIBLIO biblio) ;
static void biblioPrepare (BIBLIO biblio) ;
static void biblioFormatDump (BIBLIO biblio) ;
static void biblioDoTitleGrep (BIBLIO biblio) ;
static void formatNormal (void) ;
static void formatCell (void) ;
static void formatGenDev (void) ;
static void biblioPageUp (void) ;
static void biblioPageDown (void) ;
static void biblioKbd (int k) ;
static int nbBox = 6 ; /* nombre de boxes par reference */
static BOOL discarding = FALSE ;
static float top = 7 ;
static double yOld ;
#define graphBoxBox(_box) { \
               float _x1, _y1, _x2, _y2 ; \
               graphBoxDim (_box, &_x1, &_y1, &_x2, &_y2) ; \
               graphRectangle (_x1 - 0.4, _y1 - 0.1, _x2 + 0.4, _y2 + 0.1) ; \
                }


#define BIBLIOGET(name) BIBLIO biblio ; \
                      if (!graphAssFind (&GRAPH2BIBLIO_ASSOC, &biblio)) \
		        messcrash ("%s() could not find BIBLIO on graph",name) ; \
		      if (!biblio) \
                        messcrash ("%s() received NULL BIBLIO pointer",name) ; \
                      if (biblio->magic != &BIBLIO_MAGIC) \
                        messcrash ("%s received non-magic BIBLIO pointer",name)

/*************************************************************************/

static void biblioDestroy(void)
{
  BIBLIOGET ("biblioDestroy") ;

  keySetDestroy (biblio->bibSet) ;
  stackDestroy (biblio->s) ;
  arrayDestroy (biblio->discard) ;
  arrayDestroy (biblio->abstract) ;
  arrayDestroy (biblio->abstExists) ;
  arrayDestroy (biblio->newNumber) ;
  arrayDestroy (biblio->sortNumber) ;
  arrayDestroy (biblio->author) ;
	arrayDestroy (biblio->date) ;
	arrayDestroy (biblio->ref) ;
  arrayDestroy (biblio->unDiscard) ;
  messfree (biblio->cp) ;
  biblio->magic = 0 ;
  messfree (biblio) ;

  return;
} /* biblioDestroy */

/*************************************************************************/

static void biblioAbstracts (void)
{
  int i ;
  BIBLIOGET ("BiblioAbstracts") ;

  if (!graphCheckEditors (graphActive(), 0))
      return ;  
  biblio->showAbstracts = !biblio->showAbstracts ; /* toggle */
  if (biblio->showAbstracts)
    for (i = 0; i < arrayMax (biblio->bibSet); i++)
      array (biblio->abstract, i, char) = 1 ;
  else
    for (i = 0; i < arrayMax (biblio->bibSet); i++)
      array (biblio->abstract, i, char) = 0 ;  
	biblio->topLine = top ;
	strcpy(biblio->mot, "\000") ;
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ;
	  return;
} /* biblioAbstracts */

/*************************************************************************/

static void biblioPrint80(void)
{
  int max ;
  Graph old = graphActive();
  BIBLIOGET ("BiblioPrint80") ;

  if (!graphCheckEditors (graphActive(), 0))
    return ;

  if (!biblio->bibSet || !(max = keySetMax(biblio->bibSet)))
    { messout ("Nothing to export, this biblio is empty, , sorry") ;
      return ;
    }
  biblioDraw80(biblio);
  graphActivate (old);  
  biblioReDraw();

  return;
} /* biblioPrint80 */

/*************************************************************************/


/* mhmp 06 mai 97 */
static void biblioPostBuffer (int ii)
{
  int i = 0 ;
  char *cp, *cq, cc ;
  BIBLIOGET ("biblioPostBuffer") ;

  cp = cq = stackText(biblio->s, 0) ;
  while (i <= ii)
    { while (*cq && *cq != '\n') cq++ ;
      cp++ ;
      if (*cp == '_')
	i++ ;
      cp = ++cq ;
    }
  cp ++ ;
  while (TRUE)
    { while (*cq && *cq != '\n') cq++ ;
      cq++ ;
      if (*cq == '_' || !*cq)
	{ cc = *cq ; *cq = 0 ;
	  graphPostBuffer (cp) ;
	  *cq++ = cc ; 
	  break ;
	}
    }

  return;
} /* biblioPostBuffer */

/*************************************************************************/

static void biblioPick (int k)
{
  int b, bb, i, j, ii = -1 ;

  BIBLIOGET ("BiblioPick") ;

  if (!graphCheckEditors (graphActive(), 0))
    return ;
	if (k == biblio->motBox) {
		graphTextEntry (biblio->mot, 0, 0, 0, 0) ;
		biblio->nbLight = 0 ;
		return ;
	}
	if (k == biblio->searchBox) {
		graphTextEntry (biblio->mot, 0, 0, 0, 0) ;

		if (!biblio->fini)
			biblio->nbLight = 0 ;
		biblio->search = 1 ;
		biblioDoReDraw (biblio) ;

		if (!strlen(biblio->mot)){
			messout ("Please specify a search string");
			return ;
		}
		if (!biblio->nbLight)
			messout ("No %s!...", biblio->mot) ;
		return ;
	}
	/*if (k < biblio->firstBox - 1)
		return ;*/
	if (k < biblio->firstBox - 1)
		k = biblio->firstBox - 1 - nbBox ;
  bb = k - biblio->zero ;
  b = bb / nbBox ;
  if (bb >= 0 && b < arrayMax (biblio->bibSet)) 
    { j = -1 ; 
      i = ii = -1 ;
      while (j != b)
	{ i = i + 1 ; ii = array (biblio->sortNumber, i, int) -1 ;
	  if (ii != -1 && !array (biblio->discard, ii, char))
	    j = j + 1 ;
	}
      if (!(bb % nbBox))
	{
	  biblioPostBuffer (b) ;
	  if (bb == biblio->curr && bb >= 0 && ii != -1)
	    display (keySet(biblio->bibSet, ii), 0, 0) ;
	  else 
	    { if (biblio->curr > -1)
		{ graphBoxDraw (biblio->zero + biblio->curr, BLACK, WHITE) ;
		  graphBoxDraw (biblio->zero + biblio->curr + 2, BLACK, WHITE) ;
		}
	      biblio->curr = bb ;
				graphBoxDraw (biblio->zero + biblio->curr, BLACK, PALEYELLOW) ;
	      if (ii != -1 && !array (biblio->abstExists, ii, char))
		graphBoxDraw (biblio->zero + biblio->curr + 2, BLACK, PALEYELLOW) ;    
	    }
	}
      else
	if (!((bb - 1) % nbBox))
	  { array (biblio->discard, ii, char) = 1 ;
	    array (biblio->unDiscard, biblio->nbDiscard, int) = ii ;
	    biblio->nbDiscard ++ ;
	    discarding = TRUE ;
			biblio->dep = 1 ;
			biblio->searchDep = 1 ;
	    biblio->numberRef = biblio->numberRef - 1 ;
	    if (bb < biblio->curr)
	      {
		biblio->curr = biblio->curr - nbBox ;
		biblio->discardUp ++ ;
	      }
	    biblioPrepare (biblio) ;
	    biblioDoReDraw (biblio) ;
	  }
	else
	  if (!((bb - 2) % nbBox))
	    if (array (biblio->abstExists, ii, char))
	      { array (biblio->abstract, ii, char) = !array (biblio->abstract, ii, char) ;
		biblio->dep = 1 ;
		biblio->searchDep = 1 ;		
		biblioPrepare (biblio) ;
		biblioDoReDraw (biblio) ;
	      }
	
    }
  else
    { if (biblio->curr > -1)
	{ graphBoxDraw (biblio->zero + biblio->curr, BLACK, WHITE) ;
	  graphBoxDraw (biblio->zero + biblio->curr + 2, BLACK, WHITE) ;
	}
      biblio->curr = -1 ;
    }

  return;
} /* biblioPick */

/*************************************************************************/

static void biblioSort (void)
{
  int i, j, rmin ;
	Array sort, new ;
  float min ;
  BIBLIOGET ("biblioSort") ;

  if (!graphCheckEditors (graphActive(), 0))
      return ;
	sort = arrayCreate (arrayMax(biblio->bibSet), int) ;
	new = arrayCreate (arrayMax(biblio->bibSet), int) ;
  for (i = 0; i < arrayMax (biblio->bibSet) ; i++) {
		min = array (biblio->newNumber, 0, float) ;
		rmin = 0 ;
		for (j = 0; j < arrayMax (biblio->bibSet) ; j++)
			if (min > array (biblio->newNumber, j, float)) {
				min = array (biblio->newNumber, j, float) ;
				rmin = j ;
			}
		array (biblio->newNumber, rmin, float) = 100000 ;
		array (new, i, float) = min ;
		array (sort, i, int) = array (biblio->sortNumber, rmin, int) ;
	}
  biblio->curr = -1 ;
	for (i = 0; i < arrayMax (biblio->bibSet) ; i++) {
		array (biblio->sortNumber, i, int) =  array (sort, i, int) ;
		array (biblio->newNumber, i, int) =  array (new, i, int) ;	
	}		 
	arrayDestroy (sort) ;
	arrayDestroy (new) ;
	strcpy(biblio->mot, "\000") ;
	biblio->topLine = top ;
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ;

  return;
} /* biblioSort */

/*************************************************************************/

static void biblioSortAuthor (void)
{
  int i, j, rmin ;
  char *a1, *a2 ;
	Array sort, new, auth ;

  BIBLIOGET ("biblioSortAuthor") ;

  if (!graphCheckEditors (graphActive(), 0))
      return ;
	sort = arrayCreate (arrayMax(biblio->bibSet), int) ;
	new = arrayCreate (arrayMax(biblio->bibSet), int) ;
	auth = arrayCreate (arrayMax(biblio->bibSet), char *) ;
  for (i = 0; i < arrayMax (biblio->bibSet) ; i++) {
		a1 = array (biblio->author, 0, char *) ;
		rmin = 0 ;
		for (j = 0 ; j < arrayMax (biblio->bibSet); j++) {
			a2 = array (biblio->author, j, char *) ;
			if (strcasecmp (a1, a2)>0)
				{ a1 = a2 ;
	      rmin = j ;
				}
		}
		array (biblio->author, rmin, char *) = "zzzzzz" ;
		array (auth, i, char *) = a1 ;
		array (new, i, int) = array (biblio->newNumber, rmin, int) ;
		array (sort, i, int) = array (biblio->sortNumber, rmin, int) ;
	}
  biblio->curr = -1 ;
	for (i = 0; i < arrayMax (biblio->bibSet) ; i++) {
		array (biblio->sortNumber, i, int) =  array (sort, i, int) ;
		array (biblio->newNumber, i, int) =  array (new, i, int) ;	
	}		 
	arrayDestroy (sort) ;
	arrayDestroy (new) ;
	arrayDestroy (auth) ;
	strcpy(biblio->mot, "\000") ;
	biblio->topLine = top ;
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ;

  return;
} /* biblioSortAuthor */
/*************************************************************************/

static void biblioSortDate (void)
{
  int i, j, rmax ;
  int a1, a2 ;
	Array sort, new, date ;

  BIBLIOGET ("biblioSortDate") ;

  if (!graphCheckEditors (graphActive(), 0))
      return ;
	sort = arrayCreate (arrayMax(biblio->bibSet), int) ;
	new = arrayCreate (arrayMax(biblio->bibSet), int) ;
	date = arrayCreate (arrayMax(biblio->bibSet), int) ;
  for (i = 0; i < arrayMax (biblio->bibSet) ; i++) {
		a1 = array (biblio->date, 0, int) ;
		rmax = 0 ;
		for (j = 0; j < arrayMax (biblio->bibSet); j++) {
			a2 = array (biblio->date, j, int) ;
			if (a1 < a2) {
				a1 = a2 ;
	      rmax = j ;
			}
		}
		array (biblio->date, rmax, int)  = -1 ;
		array (date, i, int) = a1 ;
		array (new, i, int) = array (biblio->newNumber, rmax, int) ;
		array (sort, i, int) = array (biblio->sortNumber, rmax, int) ;
	}
  biblio->curr = -1 ;
	for (i = 0; i < arrayMax (biblio->bibSet) ; i++) {
		array (biblio->sortNumber, i, int) =  array (sort, i, int) ;
		array (biblio->newNumber, i, int) =  array (new, i, int) ;	
	}		 
	arrayDestroy (sort) ;
	arrayDestroy (new) ;
	arrayDestroy (date) ;
	biblio->topLine = top ;
	strcpy(biblio->mot, "\000") ;
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ;

  return;
} /* biblioSortDate */
/*************************************************************************/

static void biblioSortRef (void)
{
  int i, j, rmin ;
  int a1, a2 ;
	Array sort, new, ref ;

  BIBLIOGET ("biblioSortRef") ;

  if (!graphCheckEditors (graphActive(), 0))
      return ;
	sort = arrayCreate (arrayMax(biblio->bibSet), int) ;
	new = arrayCreate (arrayMax(biblio->bibSet), int) ;
	ref =  arrayCreate (arrayMax(biblio->bibSet), int) ;
  for (i = 0; i < arrayMax (biblio->bibSet); i++) {
    a1 = array (biblio->ref, 0, int) ;
		rmin = 0 ;
		for (j = 0; j < arrayMax (biblio->bibSet); j++) {
			a2 = array (biblio->ref, j, int) ;
			if (a1 > a2) {
				a1 = a2 ;
	      rmin = j ;
	    }
		}
		array (biblio->ref, rmin, int)  = 10000 ;
		array (ref, i, int) = a1 ;
		array (new, i, int) = array (biblio->newNumber, rmin, int) ;
		array (sort, i, int) = array (biblio->sortNumber, rmin, int) ;
	}
  biblio->curr = -1 ;
	for (i = 0; i < arrayMax (biblio->bibSet) ; i++) {
		array (biblio->sortNumber, i, int) =  array (sort, i, int) ;
		array (biblio->newNumber, i, int) =  array (new, i, int) ;	
	}		 
	arrayDestroy (sort) ;
	arrayDestroy (new) ;
	arrayDestroy (ref) ;
	strcpy(biblio->mot, "\000") ;
	biblio->topLine = top ;
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ;

  return;
} /* biblioSortRef */

/*************************************************************************/

static void biblioUnDiscard (void)
{
  int i, ii ;
  BIBLIOGET ("biblioUnDiscard") ;

  if (!graphCheckEditors (graphActive(), 0))
      return ;
  if (biblio->nbDiscard)
    { i = array (biblio->unDiscard, biblio->nbDiscard - 1, int) ;
      array (biblio->discard, i, char) = 0 ;
      biblio->nbDiscard -- ;
      biblio->numberRef ++ ;
      ii = 0 ;
      while (i + 1 != array(biblio->sortNumber, ii, int))
	{ ii++ ;
	}
      if (( ii-biblio->discardUp) * nbBox <= biblio->curr)
	 {
	   biblio->curr = biblio->curr + nbBox ;
	   biblio->discardUp -- ;
	 }
    }
	biblio->dep = 1 ;
	biblio->searchDep = 1 ;
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ;

  return;
} /* biblioUnDiscard */
/**********************************************************/
static void biblioDown (void)
{
	BIBLIOGET ("biblioDown") ;

	if (biblio->topLine < top + 5)
		biblio->topLine = top + 5 ;
	else
		if (biblio->topLine < biblio->bottomLine - 1)
			biblio->topLine += 1 ; 
	biblio->dep = 1 ;
	biblio->searchDep = 1 ;
	biblioDoReDraw(biblio) ;
}	
/**********************************************************/
static void biblioUp (void)
{
	BIBLIOGET ("biblioUp") ;

	biblio->topLine -= 1 ;
	if (biblio->topLine < top + 5)
		biblio->topLine = top + 4 ;
	biblio->dep = 1 ;
	biblio->searchDep = 1 ;
	biblioDoReDraw(biblio) ;
}	
/**********************************************************************/
static void biblioPageDown (void)
{
	BIBLIOGET ("biblioPageDown") ;

	if (biblio->line2 >= biblio->graphHeight){
		if (biblio->topLine == top)
			biblio->topLine += biblio->graphHeight - top  -1;
		else
			biblio->topLine += biblio->graphHeight - top - 5 ; 
	}
	biblio->dep = 1 ;
	biblio->searchDep = 1 ;
	biblioDoReDraw(biblio) ;
}	
/**********************************************************/
static void biblioPageUp (void)
{
	BIBLIOGET ("biblioPageUp") ;
	if (biblio->topLine < biblio->graphHeight)
		biblio->topLine -= biblio->graphHeight - top- 1 ;
	else
		biblio->topLine -= biblio->graphHeight - top - 5 ;
	if (biblio->topLine < top)
		biblio->topLine = top ;
	biblio->dep = 1 ;
	biblio->searchDep = 1 ;
	biblioDoReDraw(biblio) ;
}	
/**********************************************************************/
static void biblioKbd (int k)
{
  BIBLIOGET ("biblioKbd") ;

  switch (k)
    {
		case UP_KEY :
      biblioUp () ;
      break ;
    case DOWN_KEY :
      biblioDown () ;
      break ;
    case PAGE_UP_KEY :
    case LEFT_KEY :
      biblioPageUp () ;
      break ;
    case PAGE_DOWN_KEY :
    case RIGHT_KEY :
      biblioPageDown () ;
      break ;
    default:
      return ;
    }
} /* biblioKbd */

/***********************************************************************/

static void biblioNewNumbers (void)
{
  BIBLIOGET ("biblioNewNumbers") ;

  biblio->showNewNumbers = !biblio->showNewNumbers ;
  if (!graphCheckEditors (graphActive(), 0))
    return ;
	biblio->dep = 1 ;
	biblio->searchDep = 1 ;
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ;

  return;
} /* biblioNewNumbers */
  
/**********************************************************************/

static void biblioRenumber (void)
{
  int i, ii, jj ;
  BIBLIOGET ("biblioRenumber") ;

  jj = 0 ;
  for (i = 0; i < arrayMax (biblio->bibSet); i++) 
    { ii = array(biblio->sortNumber, i, int) - 1 ;
      if(!array(biblio->discard, ii, char))
	jj = jj+ 1 ;
      array (biblio->newNumber, i, float) = jj ;
    }
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ;

  return;
} /* biblioRenumber */

/************************************************************************/


static MENUOPT formatMenu[] = {
  { formatNormal,	"Normal" },
  { formatCell,		"Cell" },
  { formatGenDev,	"Gen.& Dev."},
  { 0,0}
};

/***********************************************************************/

static void biblioFormat (void)
{
  BIBLIOGET ("biblioFormat") ;

  biblio->format = -1 ;
  biblioDoReDraw (biblio) ; 
  graphMenu(formatMenu) ;
  graphButtons (formatMenu, 8.5, 4.3, biblio->width) ;

  return;
} /* biblioFormat */

/***********************************************************************/

static void biblioTitleGrep (void)
{
  BIBLIOGET ("biblioTitleGrep") ;

  biblioDoTitleGrep (biblio) ;
  return;
} /* biblioFormat */

/***********************************************************************/

static MENUOPT biblioMenu1[] = {
  { graphDestroy,"Quit"},
  { help,"Help"},
  { biblioPrint80, "Print"},
  { biblioAbstracts, "Abstract on/off"}, 
  { biblioSortAuthor, "Sort by author"},
	{ biblioSortDate, "Sort by date"},
	{ biblioSortRef, "Sort by ref."},
	{0,0}
};
/***********************************************************************/

static MENUOPT biblioMenu2[] = {

  { biblioRenumber, "Renumber 1..n"},
  { biblioSort, "Sort by number"},
  { biblioUnDiscard, "UnDiscard"},
  { biblioNewNumbers, "New numbers"},
	{ biblioTitleGrep, "TitleGrep"},
	{ biblioPageUp, "Page Up"},
	{ biblioPageDown, "Page Down"},
  {0,0}
};
/***********************************************************************/

static MENUOPT biblioMenu3[] = {

	{ biblioFormat, "Format"},
  {0,0}
};
/***********************************************************************/

static void formatNormal (void)
{ 
  BIBLIOGET ("formatNormal") ;

  biblio->format = 0 ;
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ; 

  return;
} /* formatNormal */

/**********************************************************************/

static void formatCell (void)
{ 
  BIBLIOGET ("formatCell") ;

  biblio->format = 1 ;
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ; 

  return;
} /* formatCell */

/*************************************************************************/  

static void formatGenDev (void)
{ 
  BIBLIOGET ("formatGenDev") ;

  biblio->format = 2 ;
  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ; 

  return;
} /* formatGenDev */

/*************************************************************************/  

static void biblioDraw80 (BIBLIO biblio)
{
  int box, n = 0, line = 0 ;
  char *cp, *cq, cc ;
  int level ;

  biblio->s = stackReCreate (biblio->s, 3000) ;
  level = freeOutSetStack (biblio->s) ;
  biblio->width =  80;
  biblioFormatDump (biblio) ;
  freeOutClose (level) ;

  graphClear () ;
  graphText("Topic : ",1.,4.3) ;
  graphTextEntry(biblio->cp, 68, 9.,4.3,0) ;
  cp = cq = stackText(biblio->s, 0) ; n = 0 ;
  graphBoxStart () ;
  line = 6 ;
  while (TRUE)
    { while (*cq && *cq != '\n') cq++ ;
      cc = *cq ; *cq = 0 ;
      if (*cp == '_')
	{ cp++ ;
	  graphBoxEnd() ;
	  box = graphBoxStart () ; 
	  if (!n++)
	    biblio->zero = box ; 
	}
      graphText(cp, 4, line++) ;
      if (!cc)
	break ;
      *cq++ = cc ; cp = cq ; 
    }
  graphBoxEnd () ; 
  graphRedraw () ;

  graphPrint () ;
  
  return;
} /* biblioDraw80 */

/*********************************************/

static void biblioReDraw (void)
{
  BIBLIOGET ("BiblioReDraw") ;
	biblio->topLine = top ;
  biblioDoReDraw (biblio) ;
}
/**********************************************************/
static void biblioMiddleUp (double x, double y)
{ 
	float a, b ;

	BIBLIOGET ("biblioMiddleUp") ;
	
	if (y < top)
		y = top ;
	
	a = (biblio->bottomLine - top) / (biblio->graphHeight - top) ;
	b = top * (biblio->graphHeight - biblio->bottomLine) /
		(biblio->graphHeight - top) ;
	
	biblio->topLine = a * y + b ;
	if (biblio->topLine <= 11)
		biblio->topLine = 11 ;
	if (biblio->topLine > biblio->bottomLine)
		biblio->topLine = biblio->bottomLine - 1 ;
	biblioDoReDraw (biblio) ;

}
/**********************************************************/
static void biblioMiddleDrag  (double x, double y)
{ 
  BIBLIOGET ("biblioMiddleDrag") ;

  graphXorLine (0, yOld, 5, yOld) ;
  yOld = y ;
  graphXorLine (0, y, 5, y) ;

  return;
} /* biblioMiddleDrag */
/**********************************************/
static void biblioMiddleDragCentre  (double x, double y)
{ 
  BIBLIOGET ("biblioMiddleDragCentre") ;

  graphXorLine (5, yOld, biblio->graphWidth, yOld) ;
  yOld = y ;
  graphXorLine (5, y, biblio->graphWidth, y) ;

  return;
} /* biblioMiddleDragCentre */
/**********************************************************/
static void biblioMiddleUpCentre (double x, double y)
{ 
	BIBLIOGET ("biblioMiddleUpCentre") ;
	
	biblio->topLine -= (biblio->graphHeight + top) * 0.5 - y  + 2 ;
	if (biblio->topLine < top)
		biblio->topLine = top ;
	if (biblio->topLine > biblio->bottomLine)
		biblio->topLine = biblio->bottomLine - 1 ;
	biblioDoReDraw (biblio) ;

}
/**********************************************************/
static void biblioMiddleDown (double x, double y) 
{
	BIBLIOGET ("biblioMiddleDown") ;

	if (x < 5) {
		graphXorLine (0, y, 5, y) ;
		graphRegister (MIDDLE_DRAG, biblioMiddleDrag) ;
		graphRegister (MIDDLE_UP, biblioMiddleUp) ;
	}
	else {
		graphXorLine (5, y, biblio->graphWidth, y) ;
		graphRegister (MIDDLE_DRAG, biblioMiddleDragCentre) ;
		graphRegister (MIDDLE_UP, biblioMiddleUpCentre) ;
	}
	biblio->dep = 1 ;
	biblio->searchDep = 1 ;
	return ;
}

/************************************************/
static void biblioScrollBar(void)
{
	float a, b, u, v ;
	int box ;

	BIBLIOGET("biblioScrollBar");

	graphLinewidth(0.4);
	graphLine (2,top,2, biblio->graphHeight - 1) ;
	graphLinewidth(0);
	a = (biblio->graphHeight - 1 - top) /  (biblio->bottomLine - top);
	b = top * (1 - a) ;
	u = biblio->topLine * a + b ;
	v =(biblio->topLine +  biblio->graphHeight) * a + b ;
	if (v > biblio->graphHeight - 1)		v = biblio->graphHeight - 1 ;
	box = graphBoxStart();
	graphRectangle (1.6, u, 2.4, v) ;
	graphBoxEnd();
	graphBoxDraw(box,BLACK,GREEN);
	
}
/**********************************************************/
static void biblioTopic (char *txt)
{ 
  BIBLIOGET ("biblioTopic") ;

	biblioDoReDraw (biblio) ;
	return ;
}
/**********************************************************/

static void biblioSearch (char *txt)
{ 
  BIBLIOGET ("biblioSearch") ;

	biblioDoReDraw (biblio) ;

	if (!strlen(biblio->mot)){
		messout ("Please specify a search string");
		return ;
	}
	if (!biblio->nbLight){
		messout ("No %s!...", biblio->mot) ;
		return ;
	}
	messout ("Search with the Search Button") ;
	biblio->searchDep = 0 ;
	biblio->dep = 0 ;
	biblio->search = 1 ;
  return;
} /* biblioSearch */

/**********************************************************/
  
static void biblioDoReDraw (BIBLIO biblio)
{
  int box, i, ii, j, b, lightBox, nbFound = 0, pos = 0, pos2 = 0, n = 0, line = 0, isLight = 0, prems = 1 ;
	int lineOk, tip, greenLine = 0 ;
  char *cp, *cq, cc ;
  int graph_width, graph_height, length = 0 ;
	char mot [200], buf [200] ;

  graphFitBounds (&graph_width,&graph_height) ;
	if (biblio->nbLightBis && !biblio->searchDep) {
		biblio->nbLight = biblio->nbLightBis ;
		biblio->nbLightBis = 0 ;
	}
	biblio->graphWidth = graph_width ;
	biblio->graphHeight = graph_height ;
	memset(buf, 0, 200) ;
	memset(mot, 0, 200) ;
	biblio->fini = 0 ;
	biblio->firstBox = 0 ;
	line = top ;
	lineOk = 0 ;
	tip = 0 ;
	if (line >= top && line >= biblio->topLine) {
		lineOk = 1 ;
		biblio->line2 = line - biblio->topLine + top ;	
	}
	else
		tip = 4 ;
	biblio->line2 += tip ;
	if (biblio->mot && *biblio->mot)
		{
			length = strlen (biblio->mot) ;
			strncpy (buf, "*", 1);
			strcat (buf, biblio->mot) ;
      strcat (buf, "*") ;
    }
  if (biblio->width !=  graph_width - 5)
    biblioPrepare (biblio) ;
  /* max = keySetMax(biblio->bibSet) ; */
	graphClear () ;

  cp = cq = stackText(biblio->s, 0) ; n = 0 ;
  box = graphBoxStart () ;
  while (TRUE)
    { while (*cq && *cq != '\n') cq++ ;
      cc = *cq ; *cq = 0 ;

      if (*cp == '_')
				{	
					if ((!isLight) && !prems) {
						lightBox = graphBoxStart() ;
						graphBoxEnd () ;
					}
					isLight = 0 ;
					prems = 0 ;
					cp++ ;
				  graphBoxEnd() ;
					box = graphBoxStart () ; 
					if (!n++)
						biblio->zero = box ; 

					box = graphBoxStart () ;
					if (lineOk) {
						graphText ("Discard", 30, biblio->line2) ;
						if (!biblio->firstBox)
							biblio->firstBox = box ;
						graphBoxBox (box) ;
					}
					graphBoxEnd () ;

					box = graphBoxStart () ;
					b = (box - biblio->zero) / nbBox ;
					j = -1 ;
					i = ii = -1 ;
					while (j != b)
						{ i = i + 1 ; ii = array(biblio->sortNumber, i, int) - 1 ;
	            if (!array (biblio->discard, ii, char))
								j = j + 1 ;
						}

					if ( array (biblio->abstExists, ii, char)) {
						if (ii != -1 && !array (biblio->abstract, ii, char)) {
							if (lineOk){
								graphText ("Abstract on ", 40, biblio->line2) ;
								graphBoxBox (box) ;
							}
						}
						else {
							if (lineOk) {
								graphText ("Abstract off", 40, biblio->line2) ;
								graphBoxBox (box) ;
							}
						}
					}
					else
						if (lineOk)
							graphText ("               ", 40, biblio->line2) ;
					graphBoxEnd (); 
					
					if ((biblio->showNewNumbers) && (lineOk))
						graphFloatEditor ("", arrp(biblio->newNumber, i, float), 55, biblio->line2, 0) ;
					else {
						graphBoxStart() ;
						graphBoxEnd() ;
						graphBoxStart() ;
						graphBoxEnd() ; 
					}
				}
			if (!pos && *buf) {
				pos = pickMatch(cp, buf) ;
				if (biblio->searchDep && line < biblio->topLine) 
					pos = 0 ;
				if (pos) {
					if (nbFound >= biblio->nbLight) {
						if (lineOk)
							graphText(cp, 4, biblio->line2) ;
						pos2 = pos ;
						while (pos2--) cp++ ;
						cp-- ;
						strncpy(mot, cp, length) ;
						mot[length] = 0 ;
						lightBox = graphBoxStart () ;
						line++ ;
						if (lineOk)
							graphText(mot, pos + 3, biblio->line2++) ;
						else
							if (line >= top && line >= biblio->topLine) {
								lineOk = 1 ;
								biblio->line2 = line - biblio->topLine + top + tip ; 
							} 
						graphBoxEnd () ;
						if (biblio->searchDep)
							graphBoxDraw(lightBox, BLACK, PALEGREEN) ; 
						else
							graphBoxDraw(lightBox, BLACK, LIGHTGREEN) ;
						greenLine = line - biblio->graphHeight * 0.5 + 5 ;
						isLight = 1 ;
						biblio->fini = 1 ;
						biblio->nbLight++ ;
						if (line >= top && line >= biblio->topLine) {
							lineOk = 1 ;
							biblio->line2 = line - biblio->topLine + top + tip ;
						} 
					}
					else {
						nbFound++ ;
						pos = 0 ;
						line++ ;
						if (lineOk)
							graphText(cp, 4, biblio->line2++) ;
						else
							if (line >= top && line >= biblio->topLine) {
								lineOk = 1 ;
								biblio->line2 = line - biblio->topLine + top + tip ; 
							}
					}
				}
				else {
					line++ ;
					if (lineOk)
						graphText(cp, 4, biblio->line2++) ;
					else
						if (line >= top && line >= biblio->topLine) {
							lineOk = 1 ;
							biblio->line2 = line - biblio->topLine + top + tip ;   
						}
				}
			}
			else {
				line++ ;
				if (lineOk)
					graphText(cp, 4, biblio->line2++) ;
				else
					if (line >= top && line >= biblio->topLine) {
						lineOk = 1 ;
						biblio->line2 = line - biblio->topLine + top + tip ;   
					}
			}
      if (!cc)
				break ;
      *cq++ = cc ; cp = cq ; /* after the \n */
			}
  graphBoxEnd () ;
	if (!discarding)
    discarding = FALSE ; 
	biblio->searchBox = graphBoxStart () ;
	graphText ("Search:", 43, 4.5) ;
	graphBoxBox (biblio->searchBox) ;
	graphBoxEnd () ;
	biblio->motBox = graphTextScrollEntry (biblio->mot, 1023, 50, 51, 4.5, biblioSearch) ;
	
	graphText("Topic : ",1.,5.7) ;
	graphTextScrollEntry (biblio->cp, 69, 68, 9., 5.7, biblioTopic) ;
	graphText (messprintf("Found %d ref", biblio->mmax), 5, 8) ;
	graphTextEntry(biblio->cp, 0, 0, 0, 0) ;
	graphTextEntry(biblio->mot, 0, 0, 0, 0) ;
	graphEntryDisable() ;
	graphRegister (KEYBOARD, biblioKbd) ;
	biblio->bottomLine = line ;
	biblioScrollBar () ;
	graphRedraw () ;
	if (biblio->line2 <= 11 && biblio->mmax) {
		if (!biblio->flip) {
			biblio->flip = 1 ;
			biblioPageUp() ;
			return ;
		}
		else
			biblio->flip = 0 ;
	}
	if (!biblio->flop && pos && !biblio->searchDep) {
		biblio->flop = 1 ;
		biblio->topLine = greenLine ;
		if (biblio->topLine < top)
			biblio->topLine = top ;
		biblio->nbLight-- ;
		biblio->searchDep = 0 ;
		biblioDoReDraw (biblio) ;
		return ;
	}
	else
		biblio->flop = 0 ;
	if (biblio->curr < 0 || 
      biblio->curr >= arrayMax(biblio->bibSet) * nbBox ||
      biblio->zero + biblio->curr > n * nbBox)
    biblio->curr = -1 ;
  else
    { graphBoxDraw (biblio->zero + biblio->curr, BLACK, PALEYELLOW) ;
      b = (biblio->curr) / nbBox ;
      j = -1 ;
      i = ii = -1 ;
      while (j != b)
	{ i = i + 1 ;  ii = array (biblio->sortNumber, i, int) -1 ;
	  if (!array (biblio->discard, ii, char))
	    j = j + 1 ;
	}
      if (ii != -1 && !array (biblio->abstExists, ii, char))
	graphBoxDraw (biblio->zero + biblio->curr + 2, BLACK, PALEYELLOW) ;
    }
	box = graphButtons (biblioMenu1, 1.0, 1.0, 200) ;
	graphButtons (biblioMenu2, 1.0, 2.65, 200) ;
	graphButtons (biblioMenu3, 1.0, 4.3, 200) ;
	
  if (biblio->showAbstracts)
    graphBoxDraw (box + 3, BLACK, PALEBLUE) ;
  if (biblio->nbDiscard)
    graphBoxDraw (box + 9, BLACK, PALEBLUE) ;
  if (biblio->showNewNumbers)
    graphBoxDraw (box + 10, BLACK, PALEBLUE) ;
  switch  (biblio->format)
    {
    case -1:
      graphText (" ", 17.5, 4.4) ;
      break ;
    case 0:
      graphText ("Normal", 10.0, 4.4) ;
      break ;
    case 1:
      graphText ("Cell", 10.0, 4.4) ;
      break ;
    case 2:
      graphText ("Genes & Development", 10.0, 4.4) ;
      break ;
    }
	if (biblio->mot && biblio->dep) {
		if (biblio->search)
			biblio->nbLightBis = biblio->nbLight -1 ;
		biblio->search = 0 ;
		biblio->searchDep = 1 ;
		biblio->dep = 0 ;
		biblio->nbLight = 0 ;
		biblioDoReDraw (biblio) ;
		biblio->searchDep = 0 ;
		return ;
	}
}

/***********************************************************************/

static void biblioPrepare (BIBLIO biblio)
{
  int level ;
  int graph_width ;

  graphFitBounds (&graph_width, 0) ;

  biblio->s = stackReCreate (biblio->s, 3000) ;
  level = freeOutSetStack (biblio->s) ;
	biblio->width =  graph_width - 5 ;
  biblioFormatDump (biblio) ;
  freeOutClose (level) ;

} /* biblioPrepare */

/***********************************************************************/

static void cellAuthor (Array aut, int j, char *author) 
{ 
  char *cc ;
  int i, lastSpace, k = 0, nbDot = 0, nbMinus = 0 ;
  cc = name (keySet (aut, j)) ;
  lastSpace = strlen(cc) - 1 ;

/* cherche le blanc derriere le nom*/
  for (i=0; i<strlen(cc); i++)
    if (cc[i] == ' ')
      lastSpace = i ; 

/* compte . et - dans les initiales*/
  for (i=lastSpace + 1; i<strlen(cc); i++)
    if (cc[i] == '-')
      nbMinus++ ;
    else if (cc[i] == '.')
      nbDot++ ;
  memset (author, 0, 200) ;

/* initiales ? */
  if (strlen (cc) - lastSpace > 4 + nbDot + nbMinus)
    lastSpace = strlen(cc) - 1 ;

/* le nom */
  for (i=0; i<=lastSpace; i++)
    author[k++] = cc[i] ;
  if (author[k-1] == ' ')
    k = k -1 ;
  if (lastSpace != strlen(cc) - 1)
    {
      author[k++] = ',' ;
      author[k++] = ' ' ;
    }

/* les initiales */
  for (i=lastSpace + 1; i<strlen(cc); i++)
    { 
      if (cc[i] == '.') /* initiales L.S. des le depart*/
	continue ;
      author[k++] = cc[i] ;
      if (!nbMinus)
	author[k++] = '.' ; /* . derriere chaque initiale*/
    }
  author[k] = '\0' ;
}
/*************************************************************************/
static void biblioDumpCell (KEYSET bibSet, BOOL showAbstracts, int width)
{
  static Array aut = 0 ;
  int i, j, max, mmax, line, ix, ii, nbLines, nbEmpty, iabst ;
  char *cp, abst ;
  char author [200] ;
  KEY text, kk;
  OBJ Ref ;
  KEY ref ;
  Stack abs, buf = stackCreate (1000) ;
#ifndef NON_GRAPHIC
  BIBLIOGET ("biblioDumpCell");
#endif
 
 if (width <= 0) width = 1 << 24 ;
  if (!keySetExists(bibSet) || !(max = keySetMax(bibSet)))
    { freeOut ("Sorry, no related bibliography") ;
      return ;
    }
  mmax = max ;
#ifndef NON_GRAPHIC
  mmax = biblio->numberRef ;
	biblio->mmax = mmax ;
#endif
	freeOutxy (messprintf("Found %d ref", mmax),1,1);
  
  aut = arrayReCreate (aut, 20, BSunit) ;
  for (i=0, line=3 ; i < max ; i++)
    { 
      abst = 0 ;
      ii = i ;
#ifndef NON_GRAPHIC
      ii = array (biblio->sortNumber, i, int) - 1 ;
			array( biblio->author, i, char *) = "zzzzz" ;  
			array(biblio->date, i, int) = 0 ;
			array(biblio->ref, i, int) = 1000 ;
	/* mhmp 11.02 adapter les zzzzz 0 1000 aux zzzzzz -1 et 10000 des sort */
			if (array (biblio->discard, ii, char)) continue ;

#endif
      line++ ;
      ref=arr(bibSet,ii,KEY)  ;
      j = i ;

#ifndef NON_GRAPHIC
      j = array (biblio->newNumber, i, float) - 1;
#endif

      freeOutxy (messprintf("_%d) ", j + 1), 0, line++) ;

      freeOut(name(ref)) ;
#ifndef NON_GRAPHIC
			if (!strncasecmp(name(ref), "[cgc", 4))
				array( biblio->ref, i, int) = 2 ;
			else
				if (!strncasecmp(name(ref), "[pm", 3))
					array( biblio->ref, i, int) = 1 ;
				else
					if (!strncasecmp(name(ref), "[w", 2))
						array( biblio->ref, i, int) = 4 ;
					else
						array( biblio->ref, i, int) = 3 ;
#endif
      Ref = bsCreate(ref) ;
      if (!Ref) continue ;
      stackClear (buf) ;
      pushText (buf, "") ;
      if (bsFindTag (Ref, _Author) &&
	  bsFlatten (Ref, 1, aut))
	{ stackClear (buf) ;
	  for (j=0 ; j < arrayMax(aut) ; j++)
	    { 
	      cellAuthor (aut, j, author) ;
	      catText  (buf, messprintf ("%s", author)) ;
	      if (j ==  arrayMax(aut) - 2)
		catText (buf, ", and ") ;
	      else
		if (j != arrayMax (aut) - 1)
		    catText (buf, ", ") ;

#ifndef NON_GRAPHIC
		if (!j)
		  array( biblio->author, i, char *) = name (keySet (aut, j)) ; 
#endif
	    }
	}
      if (bsGetData (Ref,_Year,_Int,&ix)) {
				catText (buf, messprintf (" (%d). ", ix)) ;
#ifndef NON_GRAPHIC
					array( biblio->date, i, int) = ix ; 
#endif	
			}
      if(bsGetKey (Ref,_Title,&text))
	catText (buf, messprintf ("%s ", name(text))) ;

      if (bsGetKey (Ref, _Journal, &kk))
	catText (buf, name (kk)) ;
 
      if (bsGetData (Ref, _Volume, _Int, &ix))
	{ 
	  catText (buf, messprintf (" %d",ix)) ;
	  while (bsGetData (Ref,_bsRight,_Text,&cp))
	    catText (buf, messprintf (" %s,", cp)) ;
	}
      else if (bsGetData (Ref, _Volume, _Text, &cp))
	{ 
	  catText (buf, messprintf (" %s,", cp)) ;
	  while (bsGetData (Ref,_bsRight,_Text,&cp))
	    catText (buf, messprintf (" %s", cp)) ;
	}
     
      if (bsGetData (Ref, _Page, _Int, &ix))
	{ 
	  catText (buf, messprintf (" %d",ix)) ;
	  while (bsGetData (Ref,_bsRight,_Int,&ix))
	    catText (buf, messprintf ("-%d", ix)) ;
	  catText (buf, ".") ;
	}
      else if (bsGetData (Ref,_Page, _Text, &cp))
	{ 
	  catText (buf, messprintf (" %s", cp)) ;
	  while (bsGetData (Ref,_bsRight,_Text,&cp))
	    catText (buf, messprintf ("-%s",cp)) ;
	  catText (buf, ".") ;
	}
      uLinesText (stackText(buf,0),width - 4) ;
      while ((cp = uNextLine(stackText(buf,0))))
	freeOutxy (cp,4,line++);
#ifndef NON_GRAPHIC
      if (showAbstracts || array (biblio->abstExists, ii, char))
	{ if (!array (biblio->abstract, ii, char)) 
	  { bsDestroy (Ref) ;
	    continue ;
	  }
	  abst = 1;
	}
#endif
      if (showAbstracts || abst)
	{ if (bsGetKey (Ref, _Abstract, &text) &&
	      (abs = stackGet (text)))
	    { stackClear (buf) ;
	      stackCursor (abs, 0) ;
	      nbLines = 0 ;
	      nbEmpty = 0 ;
	      while ((cp = stackNextText (abs)))
		{ catText (buf, cp) ;
		  catText (buf, "\n") ;	
		  nbLines++ ;
		  if (strlen(cp)) {
		    nbEmpty = 0 ;
		    /* 12.03.02 mhmp */
		    if (strlen(cp) > width - 6 )
			nbLines++ ;
		  }
		  else
		    nbEmpty++ ;
		}
	      stackDestroy (abs) ;
	      uLinesText (stackText(buf,0),width - 6) ;
	      for ( iabst = 1 ; iabst <= nbLines - nbEmpty ; iabst++)
		{ cp = uNextLine(stackText(buf,0)) ;
			freeOutxy (cp,6,line++) ;
		}
	    }
	}
      bsDestroy(Ref) ;
    }
  arrayDestroy (aut);
  stackDestroy (buf) ;
}

/****************************************/
static void genDevAuthor (Array aut, int j, char *author) 
{ 
  char *cc ;
  int i, lastSpace , k = 0, nbDot = 0, nbMinus = 0 ;
  cc = name (keySet (aut, j)) ;
  lastSpace = strlen(cc) - 1 ;

  /* cherche le blanc derriere le nom */
  for (i=0; i<strlen(cc); i++)
    if (cc[i] == ' ')
      lastSpace = i ;

  /* compte . et - dans les initiales*/
  for (i=lastSpace + 1; i<strlen(cc); i++)
    if (cc[i] == '-')
      nbMinus++ ;
    else if (cc[i] == '.')
      nbDot++ ;

  /* Jesus Christ mais pas Capet L.S. ou Guevara Che*/
  if (strlen (cc) - lastSpace > 4 + nbDot + nbMinus)
    lastSpace = strlen(cc) - 1 ;
  
  memset (author, 0, 200) ;

  /* initiales*/
  for (i=lastSpace + 1; i<strlen(cc); i++)
    { 
      if (cc[i] == '.')
	continue ;
      author[k++] = cc[i] ;
      if (!nbMinus)
	author[k++] = '.' ; /* . derriere chaque initiale*/      
    }
  if (lastSpace != strlen(cc) -1)/* initiales ? */
      author[k++] = ' ' ;

  /* le nom */
  for (i=0; i<=lastSpace; i++)
    author[k++] = cc[i] ;
  if (author[k-1] == ' ')
    k = k -1 ;
  author[k] = '\0' ;
}
/*************************************************************************/
static void biblioDumpGenDev (KEYSET bibSet, BOOL showAbstracts, int width)
{
  static Array aut = 0 ;
  int i, j, max, mmax, line, ix, ii, nbLines, nbEmpty, iabst ;
  char *cp, abst ;
  char author [200] ;
  KEY text, kk;
  OBJ Ref ;
  KEY ref ;
  Stack abs, buf = stackCreate (1000) ;
#ifndef NON_GRAPHIC
  BIBLIOGET ("biblioDumpGenDev");
#endif

  if (width <= 0) width = 1 << 24 ;
  if (!keySetExists(bibSet) || !(max = keySetMax(bibSet)))
    { freeOut ("Sorry, no related bibliography") ;
      return ;
    }
  mmax = max ;
#ifndef NON_GRAPHIC
  mmax = biblio->numberRef ;
	biblio->mmax = mmax ;
#endif
	freeOutxy (messprintf("Found %d ref", mmax),1,1);
  
  aut = arrayReCreate (aut, 20, BSunit) ;
  for (i=0, line=3 ; i < max ; i++)
    { 
      abst = 0 ;
      ii = i ;
#ifndef NON_GRAPHIC
      ii = array (biblio->sortNumber, i, int) - 1 ;
			array( biblio->author, i, char *) = "zzzzz" ;  
			array(biblio->date, i, int) = 0 ;
			array(biblio->ref, i, int) = 1000 ;
	/* mhmp 11.02 adapter les zzzzz 0 1000 aux zzzzzz -1 et 10000 des sort */
      if (array (biblio->discard, ii, char)) continue ;

#endif
      line++ ;
      ref=arr(bibSet,ii,KEY)  ;
      j = i ;

#ifndef NON_GRAPHIC
      j = array (biblio->newNumber, i, float) - 1;
#endif

      freeOutxy (messprintf("_%d) ", j + 1), 0, line++) ;

      freeOut(name(ref)) ;
#ifndef NON_GRAPHIC
			if (!strncasecmp(name(ref), "[cgc", 4))
				array( biblio->ref, i, int) = 2 ;
			else
				if (!strncasecmp(name(ref), "[pm", 3))
					array( biblio->ref, i, int) = 1 ;
				else
					if (!strncasecmp(name(ref), "[w", 2))
						array( biblio->ref, i, int) = 4 ;
					else
						array( biblio->ref, i, int) = 3 ;
#endif
      Ref = bsCreate(ref) ;
      if (!Ref) continue ;
      stackClear (buf) ;
      pushText (buf, "") ;
      if (bsFindTag (Ref, _Author) &&
	  bsFlatten (Ref, 1, aut))
	{ stackClear (buf) ;
	  for (j=0 ; j < arrayMax(aut) ; j++)
	    { 
	      if (!j)
		cellAuthor (aut, j, author) ;
	      else
		genDevAuthor (aut, j, author) ;
	      catText  (buf, messprintf ("%s", author)) ;
	      if (j ==  arrayMax(aut) - 2)
		{
		  if (arrayMax (aut) > 2)
		    catText (buf, ", and ") ;
		  else
		    catText (buf, " and ") ;
		}
	      else
		{
		  if (j != arrayMax (aut) - 1)
		    catText (buf, ", ") ;
		}
	      
#ifndef NON_GRAPHIC
		if (!j)
		  array( biblio->author, i, char *) = name (keySet (aut, j)) ; 
#endif /* !NON_GRAPHIC */
	    }
	  if (arrayMax(aut) > 1)
	     catText (buf, ".") ;
	}
      if (bsGetData (Ref,_Year,_Int,&ix)) {
				catText (buf, messprintf (" %d. ", ix)) ;
#ifndef NON_GRAPHIC
					array( biblio->date, i, int) = ix ; 
#endif	
			}
      if(bsGetKey (Ref,_Title,&text))
	catText (buf, messprintf ("%s ", name(text))) ;

      if (bsGetKey (Ref, _Journal, &kk))
	catText (buf, name (kk)) ;
 
      if (bsGetData (Ref, _Volume, _Int, &ix))
	{ 
	  catText (buf, messprintf (" %d",ix)) ;
	  while (bsGetData (Ref,_bsRight,_Text,&cp))
	    catText (buf, messprintf (" %s", cp)) ;
	}
      else if (bsGetData (Ref, _Volume, _Text, &cp))
	{ 
	  catText (buf, messprintf (" %s:", cp)) ;
	  while (bsGetData (Ref,_bsRight,_Text,&cp))
	    catText (buf, messprintf (" %s", cp)) ;
	}
     
      if (bsGetData (Ref, _Page, _Int, &ix))
	{ 
	  catText (buf, messprintf (" %d:",ix)) ;
	  while (bsGetData (Ref,_bsRight,_Int,&ix))
	    catText (buf, messprintf ("-%d", ix)) ;
	  catText (buf, ".") ;
	}
      else if (bsGetData (Ref,_Page, _Text, &cp))
	{ 
	  catText (buf, messprintf (" %s", cp)) ;
	  while (bsGetData (Ref,_bsRight,_Text,&cp))
	    catText (buf, messprintf ("-%s",cp)) ;
	  catText (buf, ".") ;
	}
      uLinesText (stackText(buf,0),width - 4) ;
      while ((cp = uNextLine(stackText(buf,0))))
	freeOutxy (cp,4,line++);
#ifndef NON_GRAPHIC
      if (showAbstracts || array (biblio->abstExists, ii, char))
	{ if (!array (biblio->abstract, ii, char)) 
	  { bsDestroy (Ref) ;
	    continue ;
	  }
	  abst = 1;
	}
#endif
      if (showAbstracts || abst)
	{ if (bsGetKey (Ref, _Abstract, &text) &&
	      (abs = stackGet (text)))
	    { stackClear (buf) ;
	      stackCursor (abs, 0) ;
	      nbLines = 0 ;
	      nbEmpty = 0 ;
	      while ((cp = stackNextText (abs)))
		{ catText (buf, cp) ;
		  catText (buf, "\n") ;	
		  nbLines++ ;
		  if (strlen(cp)) {
		    nbEmpty = 0 ;
		    /* 12.03.02 mhmp */
		    if (strlen(cp) > width - 6 )
		      nbLines++ ;
		  }
		  else
		    nbEmpty++ ;
		}
	      stackDestroy (abs) ;
	      uLinesText (stackText(buf,0),width - 6) ;
	      for ( iabst = 1 ; iabst <= nbLines - nbEmpty ; iabst++)
		{ cp = uNextLine(stackText(buf,0)) ;
	          freeOutxy (cp,6,line++);
		}
	    }
	}
      bsDestroy(Ref) ;
    }
  arrayDestroy (aut);
  stackDestroy (buf) ;
}

/*************************************************************************/

static void biblioDisplay (char *title, KEYSET bibSet)
{
  BIBLIO biblio ;
  int i ;
  OBJ obj ;
  KEY text ;

  if (!title || !keySetExists (bibSet))
    messcrash ("Bad call to biblioCreate") ;
  if (!keySetMax (bibSet))
    { messout ("Sorry, no associated biblio") ;
      return ;
    }

  biblio = (BIBLIO) messalloc (sizeof(struct BiblioStruct)) ;
  biblio->magic = &BIBLIO_MAGIC ;
  biblio->bibSet = bibSet ;
  biblio->s = stackCreate (3000) ;
  biblio->curr = -1 ;
  biblio->nbDiscard = 0 ;
  biblio->discardUp = 0 ;
	biblio->topLine = top ;
	biblio->flip = 0 ;
	biblio->flop = 0 ;
	biblio->dep = 0 ;
	biblio->search = 0 ;
	biblio->searchDep = 0 ;
  biblio->format = -1 ;
	biblio->nbLightBis = 0 ;
  biblio->numberRef = arrayMax(biblio->bibSet) ;
  biblio->cp = messalloc(70) ;
  biblio->discard = arrayCreate (arrayMax(biblio->bibSet), char) ;
  biblio->abstract = arrayCreate (arrayMax(biblio->bibSet), char) ;
  biblio->abstExists = arrayCreate (arrayMax(biblio->bibSet), char) ;
  biblio->newNumber = arrayCreate (arrayMax(biblio->bibSet), float) ;
  biblio->sortNumber = arrayCreate (arrayMax(biblio->bibSet), int) ;
  biblio->author = arrayCreate (arrayMax(biblio->bibSet), char *) ;
	biblio->date = arrayCreate (arrayMax(biblio->bibSet), int) ;
	biblio->ref = arrayCreate (arrayMax(biblio->bibSet), int) ;
  biblio->unDiscard = arrayCreate (arrayMax(biblio->bibSet), int) ;
  for (i = 0; i< arrayMax(biblio->bibSet); i++) 
       { array (biblio->newNumber, i, float) = i + 1 ;
	 array (biblio->sortNumber, i, int) = i + 1 ;
	 if ((obj = bsCreate (arr (bibSet, i, KEY))))
	   { if (bsGetKey (obj, _Abstract, &text) && iskey(text) == 2)
	     array (biblio->abstExists, i, char) = 1 ;
	   bsDestroy (obj) ;
	   }
       }
  strncpy(biblio->cp, title,69); /* *(cp+70) == 0 by messalloc */
  
  displayCreate(DtBiblio) ;
  graphRetitle (title) ;
  graphHelp ("Biblio") ;
  graphRegister (DESTROY, biblioDestroy) ;
  graphRegister (PICK, biblioPick) ;
  graphRegister (RESIZE, biblioReDraw) ;
	graphRegister (MIDDLE_DOWN, biblioMiddleDown) ;
	graphRegister (KEYBOARD, biblioKbd) ;
  graphAssociate (&GRAPH2BIBLIO_ASSOC, biblio) ;

  biblioPrepare (biblio) ;
  biblioDoReDraw (biblio) ;  

  return;
} /* biblioDisplay */

/***********************************************************************/

static void biblioFormatDump (BIBLIO biblio)
{
  switch (biblio->format)
    {
    case -1 :
    case 0 :
      biblioDump (biblio->bibSet, biblio->showAbstracts, biblio->width) ;
      break ;

    case 1 :
      biblioDumpCell (biblio->bibSet, biblio->showAbstracts, biblio->width) ;
      break ;

    case 2 :
      biblioDumpGenDev (biblio->bibSet, biblio->showAbstracts, biblio->width) ;
      break ;
    }

 return; 
} /* biblioFormatDump */

/***********************************************************************/
/***************  Public Graphic display functions *********************/

void biblioKeySet (char *title, KEYSET s)
{ char *cr = title ? strnew (title, 0) : "Biblio" ; /* copy, if comes from messprintf */
  biblioDisplay (cr, biblioFollow (s)) ;
  if (title) messfree (cr) ;
} /* biblioKeySet */

/***********/

void biblioKey(KEY key)
{
  KEYSET set = keySetCreate() ;

  keySet(set,0) = key ;
  biblioKeySet (messprintf("Bibliography attached to %s",
			 name(key)), set) ;
} /* biblioKey */

#endif /* NON_GRAPHIC */

/***********************************************************************/
/***************  Public elementary functions  *************************/

static BOOL biblioKeyPossible2 (KEY key, BOOL doOpen)
{ 
  KEYSET ks = 0 ;
  int i ;
  KEY *kp ;
  BOOL found = FALSE ;

  if (!key || !class(key) || (pickType(key) != 'B'))
    return FALSE ;
  else if (class(key) == _VPaper || class(key) == _VLongText)
    return TRUE ;
  else if (!bsIsClassInClass (class(key), _VPaper))
    return FALSE ;
  if (!doOpen) return TRUE ;
  ks =  bsKeySet(key) ;
  if (!ks) return FALSE ;
  kp = arrp(ks, 0, KEY) - 1 ;
  i = keySetMax(ks) ;
  while (kp++, i--)
    if (class(*kp) == _VPaper)
      { found = TRUE ; break ; }
  keySetDestroy (ks) ;

  return found ;
} /* biblioKeyPossible2 */

BOOL biblioKeyPossible (KEY key)
{ return  biblioKeyPossible2 (key, TRUE) ; }

BOOL biblioPossible (KEYSET s)
{ 
  KEY key = 0 ;
  int i , dx, max ;

  if (!keySetExists (s) || !(max = keySetMax (s)))
    return FALSE ;
  /* try just on class names, very cheap */
  dx = max/60 ; if (!dx) dx = 1 ;
  for (i = 0; i < max ; i += dx)
    {
      key = keySet(s, i) ;
			if (class(key) == _VPaper || class(key) == _VLongText || class(key) == _VPerson)
	return TRUE ;
    }
  /* try again, opening just the models */
  dx = max/10 ; if (!dx) dx = 1 ;
  for (i = 0; i < max ; i += dx)
    {
      key = keySet(s, i) ;
      if (biblioKeyPossible2 (key, FALSE))
	goto maybe ;
    }
  return FALSE ;
 
maybe:
   if ((TRUE || externalServer)
       && max > 1 ) return TRUE ; /* do not wait on server, or probably on any machine */
 /* try again, possibly opening objects */
  dx = max/10 ; if (!dx) dx = 1 ;
  for (i = 0; i < max ; i += dx)
    {
      key = keySet(s, i) ;
      if (biblioKeyPossible2 (key, TRUE))
	return TRUE ;
    }
  return FALSE ;
} /* biblioPossible */


KEYSET biblioFollow (KEYSET s)
{
  KEYSET bib1, bib2, bib3, biba ;

  if (s && keySetMax(s))
    {
      bib1 = query (s, "{NEIGHBOURS} $| { CLASS Person;  >Author; >Paper } $| {>Quoted_in}") ;
                  /* was {>Paper} $| {>Reference} $| {>Quoted_in}") ; */
      bib2 = keySetOR (s, bib1) ;
      bib3 = query (bib2, "CLASS Paper") ; 
      biba = keySetAlphaHeap(bib3, keySetMax(bib3)) ;
      
      keySetDestroy(bib1);
      keySetDestroy(bib2);
      keySetDestroy(bib3);
    }
  else
    biba = keySetCreate () ;

  return biba ;
} /* biblioFollow */


/***********************************************************************/

void biblioDump (KEYSET bibSet, BOOL showAbstracts, int width)
{
  static Array aut = 0 ;
  int i, j, max, mmax, line, ix, ii, nbLines, nbEmpty, iabst ;
  char *cp, abst ;
  KEY text, kk;
  OBJ Ref ;
  KEY ref ;
  Stack abs, buf = stackCreate (3000) ;

#ifndef NON_GRAPHIC
  /* do not use BIBLIOGET here, cause it is used by gifaceserver */
	BIBLIO biblio = 0 ;

  if (graphActive())
    { 
      if (!graphAssFind (&GRAPH2BIBLIO_ASSOC, &biblio))
	biblio = 0 ;
      if (biblio && biblio->magic != &BIBLIO_MAGIC)
	biblio = 0 ;
    }
#endif

  if (width <= 0) width = 1 << 24 ;
  if (!keySetExists(bibSet) || !(max = keySetMax(bibSet)))
    { freeOut ("Sorry, no related bibliography") ;
      return ;
    }
  mmax = max ;
#ifndef NON_GRAPHIC
  mmax = biblio ? biblio->numberRef : max ;
	biblio->mmax = mmax ;
#endif
	freeOutxy (messprintf("Found %d ref", mmax),1,1);
  
  aut = arrayReCreate (aut, 20, BSunit) ;
  for (i=0, line=3 ; i < max ; i++)
    { 
      abst = 0 ;
      ii = i ;
#ifndef NON_GRAPHIC
      if(biblio)
	{ ii = array (biblio->sortNumber, i, int) - 1 ;
	array( biblio->author, i, char *) = "zzzzz" ;  
	array(biblio->date, i, int) = 0 ;
	array(biblio->ref, i, int) = 1000 ;
	/* mhmp 11.02 adapter les zzzzz 0 1000 aux zzzzzz -1 et 10000 des sort */
	if (array (biblio->discard, ii, char)) continue ;
	}

#endif
      line++ ;
      ref=arr(bibSet,ii,KEY)  ;
      j = i ;

#ifndef NON_GRAPHIC
      if (biblio)
	j = array (biblio->newNumber, i, float) - 1;
#endif

      freeOutxy (messprintf("_%d) ", j + 1), 0, line++) ;

      freeOut(name(ref)) ;
#ifndef NON_GRAPHIC
			if (!strncasecmp(name(ref), "[cgc", 4))
				array( biblio->ref, i, int) = 2 ;
			else
				if (!strncasecmp(name(ref), "[pm", 3))
					array( biblio->ref, i, int) = 1 ;
				else
					if (!strncasecmp(name(ref), "[w", 2))
						array( biblio->ref, i, int) = 4 ;
					else
						array( biblio->ref, i, int) = 3 ;
#endif
      Ref = bsCreate(ref) ;
      if (!Ref) continue ;
      
      stackClear (buf) ;
      if(bsGetKey (Ref,_Title,&text))
	{ pushText (buf, name(text)) ;
	uLinesText (stackText(buf,0),width - 4) ;
	while ((cp = uNextLine(stackText(buf,0))))
	  freeOutxy (cp,4,line++);
	}

      if (bsFindTag (Ref, _Author) &&
	  bsFlatten (Ref, 1, aut))
	{ stackClear (buf) ;
	  for (j=0 ; j < arrayMax(aut) ; j++)
	    { if (j)
		catText (buf, messprintf (", %s", name(keySet(aut,j)))) ;
	    else
	      { catText (buf, name(keySet(aut,j))) ;
#ifndef NON_GRAPHIC
		if (biblio)
		  array( biblio->author, i, char *) = name (keySet (aut, j)) ; 
#endif
	      }
	    }
	  catText (buf, ".") ;
	  uLinesText (stackText(buf,0),width - 4) ;
	  while ((cp = uNextLine(stackText(buf,0))))
	    freeOutxy (cp,4,line++);
	}

      j = 0 ;
      if (bsGetKey (Ref, _Journal, &kk))
	{ j = 1 ; freeOutxy (name(kk), 4, line++) ; }
 
      if (bsGetData (Ref, _Volume, _Int, &ix))
	{ if (!j) 
	    { j = 1 ; freeOutxy ("  ", 4, line++) ; }
	    freeOutf (" %d",ix) ;
	  while (bsGetData (Ref,_bsRight,_Text,&cp))
	    freeOutf (" %s", cp) ;
	}
      else if (bsGetData (Ref, _Volume, _Text, &cp))
	{ if (!j) 
	    { j = 1 ; freeOutxy ("  ", 4, line++) ; }
	  freeOut (" ") ; freeOut (cp) ;
	  while (bsGetData (Ref,_bsRight,_Text,&cp))
	    freeOutf (" %s", cp) ;
	}
      
      if (bsGetData (Ref, _Page, _Int, &ix))
	{ if (!j) 
	    { j = 1 ; freeOutxy ("  ", 4, line++) ; }
	  freeOutf (": %d",ix) ;
	  while (bsGetData (Ref,_bsRight,_Int,&ix))
	    freeOutf ("-%d", ix) ;
	}
      else if (bsGetData (Ref,_Page, _Text, &cp))
	{ if (!j) 
	    { j = 1 ; freeOutxy ("  ", 4, line++) ; }
	  freeOutf (", %s", cp) ;
	  while (bsGetData (Ref,_bsRight,_Text,&cp))
	    freeOutf ("-%s",cp) ;
	}
      if (bsGetData (Ref,_Year,_Int,&ix))
	{ if (!j) 
	    { j = 1 ; freeOutxy ("  ", 4, line++) ; }
	  freeOutf (" (%d)",ix) ;
#ifndef NON_GRAPHIC
			array( biblio->date, i, int) = ix ;
#endif	
	}
#ifndef NON_GRAPHIC
      if (biblio &&(showAbstracts || array (biblio->abstExists, ii, char)))
	{ if (!array (biblio->abstract, ii, char)) 
	  { bsDestroy  (Ref) ;
	    continue ;
	  }
	abst = 1;
	}
#endif
      if (showAbstracts || abst)
	{ if (bsGetKey (Ref, _Abstract, &text) &&
	      (abs = stackGet (text)))
	    { stackClear (buf) ;
	      stackCursor (abs, 0) ;
	      nbLines = 0 ;/* mhmp 07/05/97 */
	      nbEmpty = 0 ;
	      while ((cp = stackNextText (abs)))
		{
		  catText (buf, cp) ;
		  catText (buf, "\n") ;
		  nbLines++ ;
		  if (strlen(cp)) {
		    nbEmpty = 0 ;
		    /* 12.03.02 mhmp */
		    if (strlen(cp) > width - 6 )
		      nbLines++ ;
		  }
		  else
		    nbEmpty++ ;
		}
	      stackDestroy (abs) ;
	      uLinesText (stackText(buf,0),width - 6) ;
	      for ( iabst = 1 ; iabst <= nbLines - nbEmpty ; iabst++)
		{ cp = uNextLine(stackText(buf,0)) ;
			freeOutxy (cp,6,line++) ;
		}
	    }
	}	
      bsDestroy(Ref) ;
    }
  arrayDestroy (aut);
  stackDestroy (buf) ;

  return;
} /* biblioDump */

/***********************************************************************/
/**************************************************************/
/* search paper ptitles for gene names */
#ifndef NON_GRAPHIC
static void biblioDoTitleGrep (BIBLIO biblio)
{
  KEY title, paper, *classp, key, *kp ;
  char *cp, cutter;
  int ii, level ,  nX = 0 , xn, allXsn = 0 ;
  KEYSET xs = 0 , ks=0, classList=0, allXs = 0 ;
  static char dName[DIR_BUFFER_SIZE], filName[FIL_BUFFER_SIZE] ;
  FILE *g, *f ; 

  if (!messQuery("This is a titleAnalyser on all paper titles,\n"
		 "You ll be asked for on output file name.\n"
		 "do you want to proceed ?"))
    return ;

  ks = query(biblio->bibSet,"CLASS Paper ; Title") ;
  if (!keySetMax(ks))
    { messout ("First select a keyset containing paper with titles") ;
      keySetDestroy (ks) ;
      return ;
    }

  classList = keySetCreate () ;
  keySet (classList, 0) = _VLocus;
  xs = keySetCreate () ;
  allXs = keySetCreate () ;

  f = filqueryopen(dName, filName, "ace", "w","I need a place to export the results") ;
  g = filopen(filName, "bad","w") ;
  if (!f || !g)
    { 
      messout("Got a cancel file, i quit") ;
      goto abort ;
    }
  
  for (ii = 0 ; ii < keySetMax(ks) ; ii++)
    {
      paper = keySet (ks, ii) ;
      title = keyGetKey (paper, _Title) ;
      xn = 0 ; keySetMax (xs) = 0 ;

      level = freesettext(name(title),"") ;
	  	  
      freespecial("") ;
      freecard(level) ;
	  	  
      while((cp = freewordcut(" .,/;:[]{}()!@#$%^&*+=\\|\"\'~\n",&cutter)) || cutter)
	if (cp)
	  {
	    int n = keySetMax(classList) ;
	    classp = arrp(classList, 0, KEY) ;
	    while (n--)
	      if (lexword2key(cp, &key, *classp++) &&
		  iskey(key) == 2)
		{ 
		  keySet(xs, xn++) = key ;
		  keySet (allXs, allXsn++) = key ;
		  break ;
		}
	    if (islower((int)*cp) &&
		islower((int)*(cp+1)) &&
		islower((int)*(cp+2))&&
		(*(cp+3) == '-') &&
		isdigit((int)*(cp+4)))
	      fprintf(g, "Unknown Gene : %s\n", cp) ;
	  }

      if (xn)
	{ 
	  keySetSort(xs) ;
	  keySetCompress(xs) ;
	  
	  fprintf(f, "\nPaper : %s\nTitle \"%s\"\n", name(paper), name(title)) ;
	  xn = keySetMax(xs) ;
	  nX += xn ;
	  kp = arrp(xs, 0, KEY) - 1 ;
	  while (kp++, xn--)
	    fprintf(f, "%s %s \n", className(*kp), name(*kp)) ;
	  
	}
    }
  keySetSort(allXs) ;
  keySetCompress(allXs) ;
  messout ("Found %d XREF to %d loci in %d titles\n", nX, keySetMax(allXs), keySetMax(ks)) ;

 abort:
  if (f) filclose(f) ;
  if (g) filclose(g) ;

  keySetDestroy(ks) ;
  keySetDestroy(xs) ;  
  keySetDestroy(allXs) ;
  keySetDestroy(classList) ;
}
#endif 
/***********************************************************************/
/***********************************************************************/
