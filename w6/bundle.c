/*  File: dnabundle.c
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

#include "acedb.h"
#include "../wac/ac.h"
#include "vtxt.h"

#include "graph.h"
#include "key.h"		/* for keyboard codes */
#include "a.h"
#include "bs.h"
#include "freeout.h"
#include "dna.h"
#include "peptide.h"
#include "pick.h"
#include "systags.h"
#include "tags.h"		/* for E_mail */
#include "sysclass.h"
#include "session.h"
#include "display.h"
#include "query.h"
#include "help.h"		/* for helpOn() fw-981103 */
#include "bundle.h"	
#include "saucisse.h"

/************************************************************/

typedef struct bundleStruct
{ 
  int   *magic ;        /* == MAGIC */
  AC_HANDLE h ;
  Graph graph ;
  KEYSET ks, ks2, boxKs ;
  Stack s ;
	char mot [128] ;
	int searchBox, activeBox ;
	int graphHeight, graphWidth, topLine, bottomLine, line2 ;
} *BUNDLE ;


static int BUNDLE_MAGIC = 6316298;
static double yOld ;
static int wl0 = 6, th0 = 100 ;
static float top = 6 ;
static float marge = 30 ; /* A REGLER */
#define BUNDLEGET(name)     BUNDLE bundle ; \
     			  if (!graphAssFind (&BUNDLE_MAGIC, &bundle)) \
		            messcrash ("%s can't find graph",name) ; \
                          if (bundle->magic != &BUNDLE_MAGIC) \
                            messcrash ("%s received a wrong pointer",name)

static void bundleSaucisse (void) ;
static void bundleDraw (BUNDLE bundle) ;
static void bundlePageUp (void) ;
static void bundlePageDown (void) ;
static void bundleExport (void) ;
static MENUOPT bundleMenu[] = 
{
  { graphDestroy, "Quit"},
  { help, "Help me"},
  { graphPrint, "Print"}, 
  { bundleExport, "Export"},
  { bundleSaucisse, "Saucisse"},
  { bundlePageUp, "Page Up"},
  { bundlePageDown, "Page Down"},
  {  0, 0}
} ;

/**********************************************************/

static void bundlePageDown (void)
{
  BUNDLEGET ("bundlePageDown") ;
  
  if (bundle->line2 >= bundle->graphHeight)
    bundle->topLine += bundle->graphHeight - top -1 ;
  
  bundleDraw(bundle) ;
}	

/**********************************************************/
static void bundlePageUp (void)
{
	BUNDLEGET ("bundlePageUp") ;
	bundle->topLine -= bundle->graphHeight - top- 1 ;
	if (bundle->topLine < top)
		bundle->topLine = top ;
	bundleDraw(bundle) ;
}	
/**********************************************************/
static void bundleMiddleDrag  (double x, double y)
{ 
  BUNDLEGET ("bundleMiddleDrag") ;

  graphXorLine (0, yOld, 5, yOld) ;
  yOld = y ;
  graphXorLine (0, y, 5, y) ;

  return;
} /* bundleMiddleDrag */
/**********************************************************/
static void bundleMiddleDragCentre  (double x, double y)
{ 
  BUNDLEGET ("bundleMiddleDragCentre") ;

  graphXorLine (5, yOld, bundle->graphWidth, yOld) ;
  yOld = y ;
  graphXorLine (5, y, bundle->graphWidth, y) ;

  return;
} /* bundleMiddleDragCentre */
/**********************************************************/
static void bundleMiddleUp (double x, double y)
{ 
	float a, b ;

	BUNDLEGET ("bundleMiddleUp") ;
	
	if (y < top)
		y = top ;
	
	a = (bundle->bottomLine - top) / (bundle->graphHeight - top) ;
	b = top * (bundle->graphHeight - bundle->bottomLine) /
		(bundle->graphHeight - top) ;
	
	bundle->topLine = a * y + b ;
	if (bundle->topLine < top)
		bundle->topLine = top ;
	if (bundle->topLine > bundle->bottomLine)
		bundle->topLine = bundle->bottomLine - 1 ;
	bundleDraw (bundle) ;

}
/**********************************************************/
static void bundleMiddleUpCentre (double x, double y)
{ 
	BUNDLEGET ("bundleMiddleUpCentre") ;
	
  bundle->topLine -= (bundle->graphHeight + top) * 0.5 - y ;
	if (bundle->topLine < top)
		bundle->topLine = top ;
	if (bundle->topLine > bundle->bottomLine)
		bundle->topLine = bundle->bottomLine - 1 ;
	bundleDraw (bundle) ;
}
/**********************************************************/
static void bundleMiddleDown (double x, double y)
{ 
  BUNDLEGET ("bundleMiddleDown") ;

	yOld = y ;
	if (x < 5) {
		graphXorLine (0, y, 5, y) ;
		graphRegister (MIDDLE_DRAG, bundleMiddleDrag) ;
		graphRegister (MIDDLE_UP, bundleMiddleUp) ;
	}
	else {
		graphXorLine (5, y, bundle->graphWidth, y) ;
		graphRegister (MIDDLE_DRAG, bundleMiddleDragCentre) ;
		graphRegister (MIDDLE_UP, bundleMiddleUpCentre) ;
	}
  return;
} /* bundleMiddleDown */


/**********************************************************/
static BOOL bundleDestroy (void)
{
  BUNDLEGET ("bundleDestroy") ;

  ac_free (bundle->h) ;
  return TRUE ;
}

static void bundleSearch (char *txt)
{
  BUNDLEGET ("BundleSearch") ;
  bundle->topLine = top ;
  bundleDraw(bundle) ;
  return ;
}
/**********************************************************/

static void bundlePick (int k) 
{ 
  int key ;
  
  BUNDLEGET ("BundlePick") ;
  
  if (!graphCheckEditors (graphActive(), 0))
    return ;
  
  if (k == bundle->searchBox)
    {
      graphTextEntry (bundle->mot, 0, 0, 0, 0) ;
      bundleDraw(bundle) ;
      return ;
    }
  
  if ((key = keySet(bundle->boxKs, k))) 
    {
      if (k == bundle->activeBox)
	display (key, 0, TREE) ;
      else 
	{
	  if (bundle->activeBox)
	    graphBoxDraw (bundle->activeBox, BLACK, WHITE) ;
	  bundle->activeBox = k ;
	  graphPostBuffer (name(key)) ;
	  graphBoxDraw (bundle->activeBox, BLACK, RED) ;
	}
    }
}

/**********************************************************/
static BOOL bundleTest (void)
{
  char *cp ;
  
  BUNDLEGET ("BundleTest") ;
  
  return TRUE ;
  for (cp = bundle->mot ; *cp ; cp++) 
    if (!((*cp == 'a') || (*cp == 't') ||  (*cp == 'g') ||
	  (*cp == 'c') || (*cp == 'n') || (*cp == '-')))
      {
	messout("Only a, t, g, c, n and - are authorized") ;
	return FALSE ;
      }
  for (cp = bundle->mot ; *cp ; cp++) 
    if (!(*cp == 'n') && !(*cp == '-'))
      return TRUE ;
  messout("At least one  a, t, g or c is necessary") ;
  return FALSE ;
}

/**********************************************************/
static void bundleScrollBar(void)
{
	float a, b, u, v ;
	int box ;

	BUNDLEGET("bundleScrollBar");

	graphLinewidth(0.4);
	graphLine (2,top,2, bundle->graphHeight - 1) ;
	graphLinewidth(0);
	a = (bundle->graphHeight - 1 - top) /  (bundle->bottomLine - top);
	b = top * (1 - a) ;
	u = bundle->topLine * a + b ;
	v =(bundle->topLine +  bundle->graphHeight) * a + b ;
	if (v > bundle->graphHeight - 1)
		v = bundle->graphHeight - 1 ;
	box = graphBoxStart();
	graphRectangle (1.6, u, 2.4, v) ;
	graphBoxEnd();
	graphBoxDraw(box,BLACK,GREEN);
	
}
/**********************************************************/
static void bundleDraw (BUNDLE bundle)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j, lightBox, pos = 0, pos2 = 0, aligned = 0, gap = 0, line = top, oldFormat = 0, lim = 0 ;
  char *cp, *cq, *cr ;
  int lineOk, key, seqBox ;
  int n = arrayMax (bundle->ks) ;
  BOOL isSmall = FALSE ;
  int graphWidth, graphHeight ;
  int decal, dec, len, lencp, lencp2, width, width1, width2, nbLight = 0, length = 0, prems = 0 ;
  vTXT buf = vtxtHandleCreate (h) ;
  char small[130], mot[128], petit[130] ;
  
  graphRegister (PICK, bundlePick) ;
  graphFitBounds (&graphWidth,&graphHeight) ;
  bundle->graphHeight = graphHeight ;
  bundle->graphWidth = graphWidth ;
  graphClear () ;
  lineOk = 0 ;
  bundle->activeBox = 0 ;
  if (line >= top && line >= bundle->topLine && line <= bundle->topLine + bundle->graphHeight)
    
    lineOk = 1 ;
  bundle->line2 = line - bundle->topLine + top ;
  width1 = (bundle->graphWidth - marge) / 3 ;
  /* securite */
  if (width1 > 127)
    width1 = 127 ;
  if (width1 < 3)
    width1 = 3 ;
  width2 = 2 * width1 + marge ;
  vtxtClear (buf) ;
  memset(mot, 0, 128) ;
  
  if (*bundle->mot)
    {
      length = strlen (bundle->mot) ;

      for (cp = bundle->mot ; *cp ; cp++)
	switch (*cp)
	  {
	  case '?':
	  case 'n':
	    vtxtPrint (buf, ".") ;
	    break ;
	  default:
	    vtxtPrintf (buf, "%c", *cp) ;
	  }
      if (queryRegExpMatch (0, vtxtPtr(buf), 0) < 0)
	{
	  strcpy(bundle->mot, "\000") ;
	  goto done ;	  
	}
    }
  bundle->boxKs = arrayReCreate (bundle->boxKs, 1000, KEY) ;
  graphButtons (bundleMenu, 2, 2, bundle->graphWidth) ;
  graphText (messprintf("%d  active DNA ", n), 3, 3.5) ;
  graphText ("Search:", width2 - 4, 4.5) ;
  bundle->searchBox = graphTextScrollEntry (bundle->mot, 127, 50, width2 + 4, 4.5, bundleSearch) ;
  aligned = 0 ;
  for (i = 0 ; i < keySetMax(bundle->ks) ; i++)
    {
      memset(small, 0, 130) ;
      isSmall = FALSE ;
      width = (bundle->graphWidth - marge) / 3 ;
      cp = cq = cr = stackText(bundle->s, keySet(bundle->ks,i)) ;
      memset(petit, 0, 130) ;
      strncpy(petit, stackText(bundle->s, keySet(bundle->ks,i)), 127) ;
      if (strlen(stackText(bundle->s, keySet(bundle->ks,i))) > 127)
	strcat(petit, "..") ;
      lencp = strlen (cp) ;
      lencp2 = lencp ;
      if (width < lencp)
	isSmall = TRUE ;
      key = keySet(bundle->ks2, i) ;
      decal = dec = gap = 0 ;
      nbLight = 0 ;
      prems = 0 ;
      graphBoxStart () ;			
      while (TRUE)
	{
	  if (vtxtPtr (buf))
	    pos = queryRegExpMatch (cp, vtxtPtr (buf), TRUE) ;
	  else
	    pos = 0 ;
	  if (!lineOk)
	    {
	      if (pos)
		nbLight++ ;
	      break ;
	    }
	  if (pos) 
	    {
	      pos2 = pos ;
	      if (!nbLight)
		gap = width1 + 1 - pos - decal ;
	      nbLight++ ;
	      len = length ;
	      while (pos2--)
		{
		  cp++ ;
		  dec++ ;
		}
	      cp-- ;
	      dec-- ;
	      strncpy(mot, cp, length) ;
	      mot[length] = 0 ;
	      if (!prems) {
		if (isSmall && gap <0)
		  {
		    for (j = gap; j <= 0; j++) 
		      cr++ ;
		    graphText("..", width1 + marge + 3, bundle->line2) ;
		    if (lencp > 128)
		      lencp = 127 ;
		    strncpy (small, cr, lencp) ;
		    if (lencp2 > 127)
		      strcat(small, "..") ;
		    graphText(small, width1 + marge + 5, bundle->line2) ;
		    lim = 128 ;
		    memset(small, 0, 130) ;  /*utile ???  hachement!!!*/
		  }
		else
		  {
		    graphText(petit, width1 + marge + 4 + gap, bundle->line2) ;
		    lim = gap + 127 ;
		  }
		prems++ ;
	      }
	      if (pos + decal + gap <= lim)
		{
		  lightBox = graphBoxStart () ;
		  graphText(mot, pos + decal + width1 + marge + 3 + gap, bundle->line2) ;
		  graphBoxEnd () ;
		  graphBoxDraw(lightBox, BLACK, LIGHTGREEN) ; 
		}
	      if (pos + decal + gap == lim)
		graphText(messprintf(".."), pos + decal + width1 + marge + 3 + gap + length, bundle->line2) ;
	      while (len--)
		{
		  cp++ ;
		  dec++ ;
		}
	      decal = dec ;
	      pos = 0 ;
	    }
	  else 
	    {
	      if (vtxtPtr (buf))
		{
		  if (nbLight) 
		    {
		      graphText(messprintf("%d",nbLight), 5, bundle->line2) ;
		      seqBox = graphBoxStart () ;
		      oldFormat = graphTextFormat(BOLD) ;
		      graphText(messprintf("%s",name(key)), 10, bundle->line2) ;
		      graphTextFormat(oldFormat) ;
		      keySet(bundle->boxKs, seqBox) = key ;
		      graphBoxEnd () ;
		      if (isSmall) 
			{
			  strncpy (small, cq, width1 - 2) ;
			  strcat (small, "..") ;
			  graphText(small, marge, bundle->line2) ;
			}
		      else 
			{
			  graphText(petit, marge, bundle->line2) ;
			}
		    }
		}
	      else
		{
		  graphText(petit, marge, bundle->line2) ;
		  seqBox = graphBoxStart () ;
		  oldFormat = graphTextFormat(BOLD) ;
		  graphText(messprintf("%s",name(key)), 10, bundle->line2) ;
		  graphTextFormat(oldFormat) ;
		  keySet(bundle->boxKs, seqBox) = key ;
		  graphBoxEnd () ;
		}
	      break ;
	    }
	} 
      if (nbLight)
	{
	  line++ ;
	  aligned++ ;
	}
      if (! vtxtPtr (buf))
	line++ ;
      lineOk = 0 ;
      if (line >= top && line >= bundle->topLine && line <= bundle->topLine + bundle->graphHeight)
	lineOk = 1 ;
      bundle->line2 = line - bundle->topLine + top ;
      graphBoxEnd () ;
    }
  if (vtxtPtr (buf))
    graphText (messprintf("%d aligned", aligned), 25, 3.5) ;
  bundle->bottomLine = line ;
  bundleScrollBar () ;
 done:
  graphRedraw () ;
  ac_free (h) ;
}

/**********************************************************/

static void bundleExport (void)
{
  int pos = 0 ;
  int i, j = 0 ;
  char *cp ;
  KEYSET ks ;
  int  key ;
  char buf[130] ;
  Graph  g = graphActive () ;
  
  BUNDLEGET ("BundleExport") ;
  
  memset(buf, 0, 130) ;
  
  if (*bundle->mot){
    strcpy(buf, "*") ;
    strcat (buf, bundle->mot) ;
    strcat (buf, "*") ; 
    for (cp = buf ; *cp ; cp++) if (*cp == 'n') *cp = '?' ;
    if (!bundleTest())
      {
	strcpy(bundle->mot, "\000") ;
	strcpy(buf, "\000") ;
      }
  }
  ks = keySetCreate () ;
  for (i = 0 ; i < keySetMax(bundle->ks) ; i++)
    {
      cp = stackText(bundle->s, keySet(bundle->ks,i)) ;
      key = keySet(bundle->ks2, i) ;
      pos = pickMatch(cp, buf) ;
      if (pos)
	keySet(ks, j++)  = key ;
      
      if (!(*buf))
	keySet(ks, j++)  = key ;
    }
  keySetSort(ks) ;
  keySetCompress(ks) ;
  displayCreate(DtKeySet) ;
  graphRetitle("Exported keyset") ;
  keySetShow (ks,0) ; ks = 0 ; /* will destroy ks */
  keySetSelect () ;
  graphActivate (g) ;
} /* bundleExport */

/**********************************************************/

static void bundleDrawVoid (void)
{
  BUNDLEGET ("bundleDrawVoid") ;
	
  bundleDraw (bundle) ;
} /* bundleDrawVoid */

/**********************************************************/

static BOOL bundleKbd (KEY k)
{
  BUNDLEGET ("bundleKbd") ;

  if (!graphCheckEditors (graphActive(), 0))
    return FALSE ;
  graphTextEntry (bundle->mot, 0, 0, 0, 0) ;
  bundle->topLine = top ;
  bundleDraw(bundle);
  return TRUE ;
} /* bundleKbd */

/**********************************************************/

BOOL bundleDisplay (KEYSET ks, KEYSET ks2, Stack s)
{
  BUNDLE bundle ;
  AC_HANDLE h ;
  if (!displayCreate ("DtBundle"))
    return FALSE ;
  
  graphRetitle ("Bundle") ;
  graphRegister (RESIZE, bundleDrawVoid) ;
  graphRegister (DESTROY, bundleDestroy) ;
  graphRegister (KEYBOARD, bundleKbd) ;
  graphRegister (MIDDLE_DOWN, bundleMiddleDown) ;
  graphHelp ("Bundle") ;

  graphMenu (bundleMenu) ;
  h = ac_new_handle () ;
  bundle = (BUNDLE) halloc (sizeof (struct bundleStruct), h) ;
  bundle->h = h ;
  bundle->magic = &BUNDLE_MAGIC ;
  graphAssociate (&BUNDLE_MAGIC, bundle) ;
  
  bundle ->graph = graphActive() ;
  bundle->ks = arrayHandleCopy (ks, h)  ;
  bundle->ks2 = arrayHandleCopy (ks2, h)  ;
  bundle->boxKs = arrayHandleCreate (1000, KEY, h) ; 
  bundle->s = stackCopy (s, h)  ;
  bundle->topLine = top ;
  bundleDraw (bundle) ;
  return TRUE ;
}

/************************************************************/
/**********************************************************/

static void bundleSaucisse (void)
{
  int i, j, n, level, line, wl = wl0, th = th0 ;
  char *cp, *cq, cc ;
  Saucisse saucisse = 0 ;
  Array a ;
  Stack s = 0 ;
	int graphWidth, graphHeight;

  BUNDLEGET ("bundleSaucisse") ;
  graphRegister (PICK, 0) ;

	graphFitBounds (&graphWidth,&graphHeight) ;
	graphClear() ;
  n = arrayMax(bundle->ks) ;
  a = arrayCreate (120, char) ;
      
	graphButtons (bundleMenu, 2, 2, graphWidth) ;

  graphIntEditor ("wordLength", &wl0, 1, 5, 0) ;
  graphIntEditor ("threshold", &th0, 20, 5, 0) ;

  if (wl < 1) wl = 1 ;
  if (th < 1) th = 1 ;
  saucisse = saucisseCreate (n, wl, th) ;
  for (i = 0 ; i < keySetMax(bundle->ks) ; i++)
    {
      cp = stackText(bundle->s, keySet(bundle->ks,i)) ;
      for (j = 0, cp = stackText(bundle->s, keySet(bundle->ks,i)) ;
	   *cp ; j++, cp++)
	array(a, j, char) = *cp ;
      arrayMax(a) = j ;
      saucisseFill (saucisse, a) ;
    }
  s = stackCreate (6000) ;

  level = freeOutSetStack (s) ;
  saucisseShow (saucisse, 1) ;
  freeOutClose (level) ;
  cp = stackText (s,0) ; line = top ;
  while (*cp)
    {
      cq = cp ; while (*cq && *cq != '\n') cq++ ;
      cc = *cq ; *cq = 0 ;
      graphText (cp, 3, line++) ;
      if (!cc) break ;
      cp = cq + 1 ;
    }
  saucisseDestroy (saucisse) ;
  arrayDestroy (a) ;
  stackDestroy (s) ;
  graphRedraw () ;
}

/************************************************************/
/************************************************************/
/************************************************************/
