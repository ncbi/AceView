/*  file: graphprint.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@kaa.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: May 12 13:50 1999 (edgrif)
 * * May 12 13:50 1999 (edgrif): Added new graphSetBlockMode call.
 * * Sep 25 11:42 1998 (edgrif): Replaced #ifdef ACEDB with new graph/acedb
 *              interface calls. Sadly had to make pdExecute external to
 *              avoid contorted code.
 * * Jul  9 17:31 1998 (fw): for ACEDB I changed the order of the items in the
                             menu - the ACE-typical Quit is now top of the
                             list and does a "Cancel", "OK" will do print
 * Created: Tue Mar  1 11:56:08 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: graphprint.c,v 1.13 2015/08/18 23:23:52 mieg Exp $ */

#include <string.h>
#include "acedb.h"
#include "aceio.h"
#include "graph_.h"  /* to access graph->type */
#include "session.h"
#include "mydirent.h" /* for getwd below */
#include "freeout.h"

#define BUTTON_BACKGROUND PALEBLUE

typedef struct PdStruct { int magic ;
                          AC_HANDLE h ;
			  char status ;
			  char style ;
			  BOOL doCopy, doMail, doPrint ;
			  Graph target ;
			  Graph graph ;
			  char printerBuffer[80] ;
			  BOOL ACEDB_LPR ;
			  char dirBuffer[MAXPATHLEN] ;
			  char filBuffer[FIL_BUFFER_SIZE] ;
			  char copyBuffer[MAXPATHLEN+FIL_BUFFER_SIZE+5] ; /* +5 for extension */
			  char mailerBuffer[80] ;
			  char titleBuffer[80] ;
			  Array printers ;
			  int printerNum ;
			  int titleBox, mailerBox, printerBox, fileBox,
			      firstBox ;
			  int blueBox ;
			  BOOL isRotate ;
			  char scaleText[10] ;
			  float scale ;
			  int firstScaleBox ;
			  int pages ;
			  char pageText[5] ;
			  int pageBox ;
			} *PD;
static Graph pdGraph = 0 ;

#define PD_MAGIC 975632

#define PDGET(name) \
  PD pd ; \
  if (!graphAssFind ((void*)PD_MAGIC, &pd)) \
    messcrash ("graph not found in %s", name) ; \
  if (pd->magic != PD_MAGIC) \
    messcrash ("%s received a wrong pointer", name)

#define pageWidth 558.0
#define pageHeight 756.0

static PD lastPd = 0 ;
static int blueBoxLine ;

static void pdDraw (void) ;
static void pdDrawBlueBox (void) ;
static void copyNameEnd (char *in) ;

/*************************************************************************/

static void pdDestroy (void)
{ PDGET ("pdDestroy") ;

  if (pdGraph)
    pdGraph = 0 ;
}

/*************************************************************************/

static void pdPick (int box)
{ PDGET ("pdPick") ;

  if (!box)
    return ;

  if (pd->doPrint && arrayMax(pd->printers) && 
      box >= pd->firstBox && box < pd->firstBox + arrayMax(pd->printers))
    {
      /* restore old printerbox */
      graphBoxDraw (pd->firstBox + pd->printerNum, BLACK, WHITE) ;

      if (pd->ACEDB_LPR && box == pd->firstBox)
	strcpy (pd->printerBuffer, arr(pd->printers, 0, char*)) ;
      else
	{
	  sprintf (pd->printerBuffer, "lpr -P%s",
		   arr(pd->printers, box - pd->firstBox, char*)) ;
	}
      /* update the printcommand text */
      graphBoxDraw (pd->printerBox, -1, -1) ;

      /* activate new printer box */
      graphBoxDraw (box, BLACK, BUTTON_BACKGROUND) ;

      /* save the number of our current printer for redrawing */
      pd->printerNum = box - pd->firstBox ;
    }

  if (box >= pd->firstScaleBox && box <= pd->firstScaleBox+6)
    {
      int n = box - pd->firstScaleBox ;
      switch (n)
	{
	case 0:
	  pd->scale += 1 ;
	  break ;
	case 1:
	  if (pd->scale >= 1)
	    pd->scale -= 1 ;
	  break ;
	case 2:
	  pd->scale += 0.1 ;
	  break ;
	case 3:
	  if (pd->scale >= 0.1)
	    pd->scale -= 0.1 ;
	  break ;
	case 4:
	  pd->scale += 0.01 ;
	  break ;
	case 5:
	  if (pd->scale >= 0.02)
	    pd->scale -= 0.01 ;
	  break ;
	}
      if (pd->scale <= 0.01) pd->scale = 0.01 ;
      memset (pd->scaleText, 0, 10) ;
      if (pd->scale < 10.0) strcat (pd->scaleText, " ") ;
      strcat (pd->scaleText, messprintf ("%5.3f", pd->scale)) ;
      pd->scaleText[5] = 0 ;
      graphTextEntry (pd->scaleText, 0, 0, 0, 0) ;

      pdDrawBlueBox () ;
    }
}

/*************************************************************************/

void pdExecute (void)
{ 
  BOOL isColor = FALSE ;
  PDGET ("pdExecute") ;
  
  graphActivate (pd->target) ;
  
  switch (pd->style)
    {
    case 'a':
      graphASCII (pd->doCopy ? pd->copyBuffer : 0 ,
		  pd->doMail ? pd->mailerBuffer : 0,
		  pd->doPrint ? pd->printerBuffer : 0,
		  pd->titleBuffer) ;
      break ;
    case 'c':
      isColor = TRUE ; 
      /* fall thru */
    case 'p': /* postscript */
      graphPS (pd->doCopy ? pd->copyBuffer : 0 ,
	       pd->doMail ? pd->mailerBuffer : 0,
	       pd->doPrint ? pd->printerBuffer : 0, 
	       pd->titleBuffer, isColor,
	       pd->isRotate, pd->scale, pd->pages) ;
      break ;
    case 'f': /* pdf */
	if (pd->copyBuffer &&
	    *pd->copyBuffer)
	  {
	    int i = 0 ;
	    char *cp, *tmpBuffer = 0 ;

	    cp = pd->copyBuffer + strlen (pd->copyBuffer) - 1 ;
	    while (*cp != '.' && i < 5) { i++, cp--; }
	    if (cp > pd->copyBuffer && *cp == '.') ;
	    else cp =  pd->copyBuffer + strlen (pd->copyBuffer) ;
	    strcpy (cp, ".tmp") ;
	    tmpBuffer = strnew (pd->copyBuffer, 0) ;
	    strcpy (cp, ".pdf") ;
	    graphPS (tmpBuffer,
		     0,
		     0,
		     pd->titleBuffer, TRUE,
		     pd->isRotate, pd->scale, pd->pages) ;
	    system (messprintf("ps2pdfwr %s %s ; rm %s", tmpBuffer, pd->copyBuffer, tmpBuffer)) ;
	    messfree (tmpBuffer) ;
	  }
      break ;
    case 'g': 
      {
	char cc = 0, *cp ;
	int i = 0 ;
	ACEOUT fo = 0 ;
	
	if (pd->copyBuffer &&
	    *pd->copyBuffer)
	  {
	    cp = pd->copyBuffer + strlen (pd->copyBuffer) - 1 ;
	    while (*cp != '.' && i < 4) { i++, cp--; }
	    if (cp > pd->copyBuffer && *cp == '.') ;
	    else cp =  pd->copyBuffer + strlen (pd->copyBuffer) ;
	    strcpy (cp, ".gif") ;
	    if ((fo = aceOutCreateToFile (pd->copyBuffer, "w", 0)))
	      {
		graphGIF(pd->target, fo, FALSE); 
		aceOutDestroy(fo) ;
	      }
	    if (cc)
	      *cp = cc ;
	  }
      }
      break ;
      
    case 's':  /* SWF embedded drawing */ 
      {
	char *cp ;
	int i = 0 ;
	ACEOUT fo = 0 ;
	
	if (pd->copyBuffer &&
	    *pd->copyBuffer)
	  {
	    cp = pd->copyBuffer + strlen (pd->copyBuffer) - 1 ;
	    while (*cp != '.' && i < 4) { i++, cp--; }
	    if (cp > pd->copyBuffer && *cp == '.') ;
	    else cp =  pd->copyBuffer + strlen (pd->copyBuffer) ;
	    strcpy (cp, ".sc") ;
	    if ((fo = aceOutCreateToFile (pd->copyBuffer, "w", 0)))
	      {
		swfGraphExport (graphActive(), fo, 0) ; 
		aceOutDestroy (fo) ;
	      }
	  }
      }
      break ;
    }
  if (graphActivate (pd->graph))
    graphDestroy ();
  graphActivate (pd->target) ;
}

/*************************************************************************/

static void psButton (void)
{ PDGET ("psButton") ;
  
  pd->style = 'p' ;
  pdDraw() ;
}

/*****************/

static void cpsButton (void)
{ PDGET ("cpsButton") ;
  
  pd->style = 'c' ;
  pdDraw() ;
}

/*****************/

static void asciiButton (void)
{ PDGET ("asciiButton") ;
  
  pd->style = 'a' ;
  pdDraw() ;
}

/*************************************************************************/

static void pdfButton (void)
{ PDGET ("pdfButton") ;
  
  pd->style = 'f' ; 
  pd->doMail = 0 ; /* exclusive */ 
  pd->doPrint = 0 ;
  pdDraw() ;
}

/*************************************************************************/

static void swfButton (void)
{ PDGET ("swfButton") ;
  
  pd->style = 's' ; 
  pd->doMail = 0 ; /* exclusive */ 
  pd->doPrint = 0 ;
  pdDraw() ;
}

/******************************/

static void gifButton (void)
{ PDGET ("gifButton") ;
  
  pd->style = 'g' ; 
  pd->doMail = 0 ; /* exclusive */ 
  pd->doPrint = 0 ;
  pdDraw() ;
}

/*****************/

static void printerButton (void)
{ PDGET ("printerButton") ;
  
  pd->doPrint ^= 0x1 ;
  if (pd->doMail && pd->doPrint) pd->doMail = 0 ; /* exclusive */
  if (!pd->doMail && pd->doPrint)  /* shift to postscript mode mhmp */
    pd->style = 'p';
  else
    if (pd->doMail && !pd->doPrint)  /* shift to ascii mode mhmp */
      pd->style = 'a';
  pdDraw() ;
  if (pd->doPrint) /* mieg */
      graphBoxDraw(pd->printerBox, BLACK, WHITE) ;
/*   graphTextScrollEntry (pd->printerBuffer, 0, 0, 0, 0, 0) ;  */
}

/*****************/

static void mailerButton (void)
{ PDGET ("mailerButton") ;
  
  pd->doMail ^= 0x1 ;
  if (pd->doMail && pd->doPrint) pd->doPrint = 0 ; /* exclusive */
  if (pd->doMail && !pd->doPrint)  /* shift to ascii mode mhmp */
    pd->style = 'a';
  else
    if (!pd->doMail && pd->doPrint)  /* shift to postscript mode mhmp */
      pd->style = 'p';

  pdDraw() ;
  if (pd->doMail)
    graphTextScrollEntry (pd->mailerBuffer, 0, 0, 0, 0, 0) ;
}

/*****************/

static void copyButton (void)
{ PDGET ("copyButton") ;
      
  pd->doCopy ^= 0x1 ;
  pdDraw() ;
  if (pd->doCopy)
    graphTextScrollEntry(pd->copyBuffer, 0, 0, 0, 0, 0) ;
}

/*****************/

static void fileButton (void)
{ FILE *fil = 0 ;
  char *ending ;
  PDGET ("fileButton") ;
  
  copyNameEnd (pd->copyBuffer) ; /* to clean up dir- and filBuffer */
  ending = (pd->style != 'a') ? "ps" : "txt" ;
  fil = filqueryopen (pd->dirBuffer, pd->filBuffer, ending,
		      "w","Print File") ;
  sprintf (pd->copyBuffer, "%s/%s.%s", 
	   pd->dirBuffer, pd->filBuffer, ending) ;
  if (fil) { fclose (fil) ; filremove (pd->copyBuffer, "") ; }
  pdDraw() ;
}

static void landScapeButton (void)
{ PDGET ("landScapeButton") ;
      
  pd->isRotate = TRUE ;
  pdDraw() ;
}

static void portraitButton (void)
{ PDGET ("portraitButton") ;
      
  pd->isRotate = FALSE ;
  pdDraw() ;
}


static void copyNameEnd (char *input)
{
  char *cp, *c1, *c2 ;
  PDGET ("copyNameEnd") ;

  strcpy (pd->dirBuffer, input) ;
  for (cp = pd->dirBuffer; *cp ; ++cp) ;
  --cp ;
  while (*cp != '/' && cp > pd->dirBuffer) --cp ;

  strcpy (pd->filBuffer, ++cp) ;	/* save the name */
  *--cp = 0 ;	     /* cut off the '/' at the end of "dir" */

/* if dir-name contains ../ - handle this as a cd .. inside */
  cp = pd->dirBuffer ;
  while ((c1 = strstr (cp,"/../"))) /* /../ is important to handle ../../ */
    {
      c2 = c1 ; --c2 ;		/* start one char before "/" */
      /* search backwards for next "/" */
      while (*c2 != '/' && c2 > pd->dirBuffer)
	--c2 ;
      *c2 = 0 ;			/* mark this position as the end */
      c1 += 3 ;			/* skip "/.." in append string */
      strcat (c2, c1) ;		/* put both together */
    }
/* now remove /. from dir-name */
  cp = pd->dirBuffer ;
  while ((c1 = strstr (cp, "/./")))
    {
      c2 = c1 ; *c1 = 0 ;
      c2 += 2 ;
      strcat (c1, c2) ;
    }
  if (pd->dirBuffer[0] == '.' && pd->dirBuffer[1] == '/')
    { c1 = pd->dirBuffer ; c1 += 2 ;
      strcpy (pd->dirBuffer, c1) ;
    }
/* remove double slashes from the dir-name */
  cp = pd->dirBuffer ;
  while ((c1 = strstr (cp, "//")))
    {
      c2 = c1 ; *c1 = 0 ;
      c2 += 1 ;
      strcat (c1, c2) ;
    }

/* now remove the ending from the filname */
  for (cp = pd->filBuffer; *cp ; ++cp){} ; --cp ;
  while (*cp != '.' && cp > pd->filBuffer) --cp ;
  if (cp > pd->filBuffer)
    {
      ++cp ;
      if (strcmp(cp, "ps") == 0 || strcmp(cp, "txt") == 0)
	*--cp = 0 ;
    }
  if (!strcmp(pd->filBuffer, ".ps") ||
      !strcmp(pd->filBuffer, ".txt"))
    *pd->filBuffer = 0 ;
  
/* assemble the copyBuffer-path */
/* avoid a double slash at the end */
  if (pd->dirBuffer[strlen(pd->dirBuffer)-1] == '/')
    sprintf (pd->copyBuffer, "%s%s.%s", 
	     pd->dirBuffer, pd->filBuffer,
	     (pd->style != 'a') ? "ps" : "txt") ;
  else
    sprintf (pd->copyBuffer, "%s/%s.%s", 
	     pd->dirBuffer, pd->filBuffer,
	     (pd->style != 'a') ? "ps" : "txt") ;
  
  graphTextScrollEntry(pd->copyBuffer, 0, 0, 0, 0, 0) ;
}

static void scaleTextEnd (char *input)
{
  PDGET ("scaleTextEnd") ;

  pd->scale = atof (pd->scaleText) ;
  if (pd->scale <= 0) pd->scale = 0.01 ;

  pdDrawBlueBox () ;
}

/*************************************************************************/


#define MAXC 1000
static Array getPrinters (AC_HANDLE h)
{
  FILE *printcap ;
  char buf[MAXC], *cp ;
  int n ;
  Array pl = 0 ;
  int c, i = 0 ;


  pl = arrayHandleCreate (10, char*, h) ;
  n = 0 ;
  if ((cp = getenv ("ACEDB_LPR")))
    {
      array (pl, n++, char*) = strnew (cp, h) ;
    }

  if (
      (filName("/etc/printcap", 0, "r") &&
       (printcap = filopen("/etc/printcap", 0, "r")))
      ||
      (filName("/etc/printers.conf", 0, "r") &&
       (printcap = filopen("/etc/printers.conf", 0, "r")))
      )
    { while ((c = getc(printcap)) != EOF)
	if (c == '\\') 
	  while ((c == '\\' || c == ' ' || c == '\n' || c == '\t')
		 && c != EOF)
	    c = getc(printcap);
	else if (c == ':' || c == '|' || c == '#' )
	  { if ( c == ':' || c == '|')
	      { 
		buf[i++] = 0 ;
		array (pl, n++, char*) = strnew (buf, h) ; 
		i = 0;
	      }
	    while (c != '\n' && c != EOF)
	      if ((c = getc(printcap)) == '\\')
		while ((c == '\\' || c == ' ' || c == '\n' || c == '\t') 
		       && c != EOF)
		  c = getc(printcap);
	  }
	else if (i != MAXC && c != '\n' && c != '\t' && c != ' ')
	  buf[i++] = c;
      
      filclose (printcap);
    }
  return pl ;
}

/*************************************************************************/

static void pdDrawBlueBox (void)
{ 
  int line=blueBoxLine ;
  float x1, y1, x2, y2, w, h ;
  PDGET ("pdDrawBlueBox") ;

  graphActivate (pd->target) ;

  if (pd->isRotate)
    {
      x1 = 2 + 0.2 ;
      y1 = line+11-8 + 0.1 ;
      w = gActive->uw / gActive->aspect * pd->scale ;
      h = gActive->uh * pd->scale ;
      x2 = x1 + (35.7/pageHeight * w) ;
      y2 = y1 + (16.3/pageWidth * h) ;

      pd->pages = 1 + gActive->uh * pd->scale/pageWidth ;
     }
  else
    {
      x1 = 2+18-13 + 0.1 ;
      y1 = line + 0.1 ;
      w = gActive->uw / gActive->aspect * pd->scale ;
      h = gActive->uh * pd->scale ;
      x2 = x1 + (26.4/pageWidth * w) ;
      y2 = y1 + (22/pageHeight* h) ;

      pd->pages = 1 + gActive->uh * pd->scale/pageHeight ;
    } 

  graphActivate (pd->graph) ;

  if (pd->blueBox)
    {
      graphBoxClear (pd->blueBox) ;
      graphBoxClear (pd->blueBox+1) ;
      graphBoxClear (pd->blueBox+2) ;
      graphBoxClear (pd->blueBox+3) ;
    }

  pd->blueBox = graphBoxStart () ;
  graphLine (x1, y1, x2, y1) ;
  graphBoxEnd () ;
  graphBoxStart () ;
  graphLine (x2, y1, x2, y2) ;
  graphBoxEnd () ;
  graphBoxStart () ;
  graphLine (x1, y2, x2, y2) ;
  graphBoxEnd () ;
  graphBoxStart () ;
  graphLine (x1, y1, x1, y2) ;
  graphBoxEnd () ;
  graphBoxDraw (pd->blueBox, BLUE, TRANSPARENT) ;
  graphBoxDraw (pd->blueBox+1, BLUE, TRANSPARENT) ;
  graphBoxDraw (pd->blueBox+2, BLUE, TRANSPARENT) ;
  graphBoxDraw (pd->blueBox+3, BLUE, TRANSPARENT) ;

  sprintf (pd->pageText, "%d", pd->pages) ;

  graphBoxDraw (pd->pageBox, -1, -1) ;

  graphRedraw () ;
}



  /* normal style menu for non-ACEDB programs lists menu
     choices in order in which the buttons appear 
     (fw likes that, 980709) */
static MENUOPT pdMenu[] = {
  {pdExecute, " OK "},
  {graphDestroy, "Cancel"},
  {0, 0} } ;


static MENUOPT* pdGetMenu (void)
{
  /* ACEDB-GRAPH INTERFACE: if acedb has registered a different print menu   */
  /* use that, otherwise use default internal one.                           */
  if (getGraphAcedbPDMenu() != NULL) 
    return (getGraphAcedbPDMenu()) ;
  else 
    return (pdMenu) ;
} /* pdGetMenu */

static void pdDraw(void)
{
  int graphWidth, graphHeight ;
  int boxp, boxm, boxgif, boxSwf,boxPdf,  boxc, boxl, box = 0, b ;
  int line ;
  float dy = .1 ;
  double x  ;
  double y  ;
  int n=0, len, ix, iy ;
  float ux, uy, x0, y0, xmax ;
  PDGET ("pdDraw") ;
  
  graphActivate (pd->graph) ;
  graphPop () ;
  
  graphClear() ;
  graphFitBounds (&graphWidth,&graphHeight) ;

  graphMenu (pdGetMenu());

  line = 1 ;

  if (pd->doCopy || pd->doMail || pd->doPrint)
    {
      graphButton (" OK ", pdExecute, 2, line) ;
      graphButton ("Cancel", graphDestroy, 8, line) ;
    }
  else
    {
      graphButton ("Cancel", graphDestroy, 8, line) ;
    }
  line += 2 ;

  graphLine (0, line, graphWidth, line) ;

  line += 1 ;

  graphText ("Document Title :", 2, line) ;
  pd->titleBox = graphTextScrollEntry (pd->titleBuffer, 79, 40, 20,line,0) ;
  line += 2 ;
/* mhmp ex-ligne des formats */

  pd->printerBox = 0 ;
  pd->mailerBox = 0 ;
  pd->blueBox = 0 ;
  pd->pageBox = 0 ;

  boxp = graphButton ("Print", printerButton, 2, line - dy) ;
  if (pd->doPrint)
    {
      graphBoxDraw (boxp, BLACK, BUTTON_BACKGROUND) ;
      graphText("Print Command:", 10, line) ;

	  x = x0 = 25; y = y0 = line - dy ;
      if (arrayMax(pd->printers))
	{
	  xmax = graphWidth ;
	  
	  graphTextInfo (&ix, &iy, &ux, &uy) ;
	  uy += YtoUrel(7) ; 
	  n = 0 ;
	  box = 0 ;
	  while (n < arrayMax(pd->printers))
	    { 
	      len = strlen (pd->ACEDB_LPR && n == 0 ? "ACE default" :
			    arr(pd->printers, n, char*)) ;
	      if (x + ux*len + XtoUrel(15) > xmax && x != x0)
		{ x = x0 ; y += uy ;}
	      box = graphBoxStart () ;
	      graphText (pd->ACEDB_LPR && n == 0 ? "ACE default" :
			 arr(pd->printers, n, char*),
			 x + XtoUrel(3), y + YtoUrel(2)) ;
	      graphRectangle (x, y, gBox->x2 + XtoUrel(3), gBox->y2) ;
	      graphBoxEnd () ;
	      x += ux*len + XtoUrel(15) ; 
	      ++n ;
	    }
	  pd->firstBox = box + 1 - n ;
	  graphBoxDraw (pd->firstBox + pd->printerNum, 
			BLACK, BUTTON_BACKGROUND) ;
	}
      line = y + dy ;
      line += 2 ;
      pd->printerBox = graphBoxStart() ;
      graphTextPtr (pd->printerBuffer, 10, line, 50) ;
      graphBoxEnd () ;
    }
  else
    {
      graphBoxDraw (boxp, BLACK, WHITE) ;
      b = graphBoxStart () ;

      graphText("Print Command:", 10, line) ;

      x = x0 = 25; y = y0 = line - dy ;
      if (arrayMax(pd->printers))
	{
	  xmax = graphWidth ;
	  
	  graphTextInfo (&ix, &iy, &ux, &uy) ;
	  uy += YtoUrel(7) ; 
	  n = 0 ;
	  while (n < arrayMax(pd->printers))
	    { 
	      len = strlen (pd->ACEDB_LPR && n == 0 ? "ACE default" :
			    arr(pd->printers, n, char*)) ;
	      if (x + ux*len + XtoUrel(15) > xmax && x != x0)
		{ x = x0 ; y += uy ;}
	      box = graphBoxStart () ;
	      graphText (pd->ACEDB_LPR && n == 0 ? "ACE default" :
			 arr(pd->printers, n, char*),
			 x + XtoUrel(3), y + YtoUrel(2)) ;
	      graphRectangle (x, y, gBox->x2 + XtoUrel(3), gBox->y2) ;
	      graphBoxEnd () ;
	      graphBoxDraw (box, GRAY, WHITE) ;
	      x += ux*len + XtoUrel(15) ; 
	      ++n ;
	    }
	  pd->firstBox = box + 1 - n ;
	}
      line = y + dy ;
      line += 2 ;
      graphText (pd->printerBuffer, 10, line) ;

      graphBoxEnd () ;
      graphBoxDraw (b, GRAY, WHITE) ;
      graphBoxMenu (b, pdGetMenu());
    }


  line += 2 ;
  boxm = graphButton ("Mail", mailerButton, 2, line - dy) ;
  if (pd->doMail)
    { 
      graphBoxDraw (boxm, BLACK, BUTTON_BACKGROUND) ;
      graphText("Address:", 10, line) ;
      line += 2 ;
      pd->mailerBox = graphTextScrollEntry(pd->mailerBuffer,
					   60, 50, 10, line, 0) ;
    }
  else
    {
      graphBoxDraw (boxm, BLACK, WHITE) ;
      b = graphBoxStart () ;
      graphText("Address:", 10, line) ;
      line += 2 ;
      graphText(pd->mailerBuffer, 10, line) ;
      graphBoxEnd () ;
      graphBoxDraw (b, GRAY, WHITE) ;
      graphBoxMenu (b, pdGetMenu());
    }
  
  line += 2 ;
  boxc = graphButton ("Copy", copyButton, 2, line - dy) ;

  if (pd->doCopy)
    {
      graphBoxDraw (boxc, BLACK, BUTTON_BACKGROUND) ;
      graphText("keep a copy in file:", 10, line) ;
      graphButton("File chooser", fileButton, 32, line - dy) ;
      line += 2 ;
      sprintf (pd->copyBuffer, "%s/%s.%s", 
	       pd->dirBuffer, pd->filBuffer,
	       pd->style != 'a' ? "ps" : "txt") ;
      pd->fileBox = graphTextScrollEntry(pd->copyBuffer,
					 128, 50, 10, line, copyNameEnd) ;
    }
  else
    {
      graphBoxDraw (boxc, BLACK, WHITE) ;
      b = graphBoxStart () ;
      graphText("keep a copy in file:", 10, line) ;
      graphText ("File chooser", 32, line - dy) ;
      line += 2 ;
      sprintf (pd->copyBuffer, "%s/%s.%s", 
	       pd->dirBuffer, pd->filBuffer,
	       pd->style != 'a' ? "ps" : "txt") ;
      graphText (pd->copyBuffer,10, line) ;
      graphBoxEnd () ;
      graphBoxDraw (b, GRAY, WHITE) ;
      graphBoxMenu (b, pdGetMenu());
    }

  graphTextEntry (pd->titleBuffer, 0, 0, 0, 0) ;

  line += 2 ;
/* mhmp 28/04/97 */
  graphText ("Format:", 2, line - dy) ;
  boxp = graphButton ("PostScript", psButton, 11, line - dy) ;
  boxc = graphButton ("Color PS", cpsButton, 23, line - dy) ;
  boxgif = graphButton ("Gif", gifButton, 33, line - dy) ;
  boxPdf = graphButton ("PDF", pdfButton, 38, line - dy) ;
  boxSwf = graphButton ("SWF", swfButton, 43, line - dy) ;
  boxm = graphButton ("Text Only", asciiButton, 49, line - dy) ;

  if (pd->style == 'p')
    graphBoxDraw (boxp, BLACK, BUTTON_BACKGROUND) ;
  else if (pd->style == 'c')
    graphBoxDraw (boxc, BLACK, BUTTON_BACKGROUND) ;
  else if (pd->style == 'a')
    graphBoxDraw (boxm, BLACK, BUTTON_BACKGROUND) ;
  else if (pd->style == 'g')
    graphBoxDraw (boxgif, BLACK, BUTTON_BACKGROUND) ;
  else if (pd->style == 'f')
    graphBoxDraw (boxPdf, BLACK, BUTTON_BACKGROUND) ;
  else if (pd->style == 's')
    graphBoxDraw (boxSwf, BLACK, BUTTON_BACKGROUND) ;

  line += 2 ;
  graphLine (0, line, graphWidth, line) ;
  line += 1.5 ;
  graphColor (GRAY) ;
  graphFillRectangle (2, line, 2+36, line+22) ;

  if (pd->isRotate)
    {
      /* landscape */
      graphColor (WHITE) ;
      graphFillRectangle (2, line+11-8, 2+36, line+11+8) ;
      graphColor (BLACK) ;
      graphRectangle (2, line+11-8, 2+36, line+11+8) ;
    }
  else
    {
      /* portrait */
      graphColor (WHITE) ;
      graphFillRectangle (2+18-13, line, 2+18+13, line+22) ;
      graphColor (BLACK) ;
      graphRectangle (2+18-13, line, 2+18+13, line+22) ;
    }
  
  boxp = graphButton ("Portrait", portraitButton, 44, line+5) ;
  boxl = graphButton ("Landscape", landScapeButton, 44, line+7) ;
  graphBoxDraw (boxp, BLACK, pd->isRotate ? WHITE : BUTTON_BACKGROUND) ;
  graphBoxDraw (boxl, BLACK, pd->isRotate ? BUTTON_BACKGROUND : WHITE) ;

  graphColor (BLACK) ;
  graphRectangle (2, line, 2+36, line+22) ;

  graphText ("Scale :", 44, line+10) ;
  memset (pd->scaleText, 0, 10) ;
  if (pd->scale < 10.0) strcat (pd->scaleText, " ") ;
  strcat (pd->scaleText, messprintf ("%5.3f", pd->scale)) ;
  pd->scaleText[5] = 0 ;
  graphTextEntry (pd->scaleText, 6, 44+8, line+10, scaleTextEnd) ;
  graphText ("Pages : ", 44, line+13) ;



  pd->firstScaleBox = graphBoxStart () ;
  graphRectangle (44+8, line+9, 44+10, line+10) ;
  graphLine (44+8, line+10, 44+9, line+9) ;
  graphLine (44+9, line+9, 44+10, line+10) ;
  graphBoxEnd () ;
  graphBoxStart () ;
  graphRectangle (44+8, line+11, 44+10, line+12) ;
  graphLine (44+8, line+11, 44+9, line+12) ;
  graphLine (44+9, line+12, 44+10, line+11) ;
  graphBoxEnd () ;

  graphBoxStart () ;
  graphRectangle (44+11, line+9, 44+12, line+10) ;
  graphLine (44+11, line+10, 44+11.5, line+9) ;
  graphLine (44+11.5, line+9, 44+12, line+10) ;
  graphBoxEnd () ;
  graphBoxStart () ;
  graphRectangle (44+11, line+11, 44+12, line+12) ;
  graphLine (44+11, line+11, 44+11.5, line+12) ;
  graphLine (44+11.5, line+12, 44+12, line+11) ;
  graphBoxEnd () ;

  graphBoxStart () ;
  graphRectangle (44+12, line+9, 44+13, line+10) ;
  graphLine (44+12, line+10, 44+12.5, line+9) ;
  graphLine (44+12.5, line+9, 44+13, line+10) ;
  graphBoxEnd () ;
  graphBoxStart () ;
  graphRectangle (44+12, line+11, 44+13, line+12) ;
  graphLine (44+12, line+11, 44+12.5, line+12) ;
  graphLine (44+12.5, line+12, 44+13, line+11) ;
  graphBoxEnd () ;


  sprintf (pd->pageText, "%d", pd->pages) ;
  pd->pageBox = graphBoxStart () ;
  graphTextPtr (pd->pageText, 44+10, line+13, 3) ;
  graphBoxEnd () ;

  blueBoxLine = line ;

  pdDrawBlueBox () ;		/* does a redraw */
}


/*************************************************************************/
static PD pdInit(void)
{ 
  char *cp ;
  AC_HANDLE h ;

  if (lastPd)
    return lastPd ;

  h = handleCreate () ;
  lastPd = (PD) halloc (sizeof(struct PdStruct), h) ;
  lastPd->magic = PD_MAGIC ;
  lastPd->h = h ;

  /* ACEDB-GRAPH INTERFACE:                                                  */
  /* Acedb provides colour printing styles and user login info., otherwise   */
  /* this stuff is left out.                                                 */
  if (getGraphAcedbStyle() != NULL)
    {
    lastPd->style = (getGraphAcedbStyle())() ; 
    strcpy (lastPd->mailerBuffer, getGraphAcedbLogin());
    }
  else
    {
    lastPd->style = 'p' ;
    strcpy (lastPd->mailerBuffer, "");
    }

  if ((cp = getenv ("ACEDB_LPR")))
    lastPd->ACEDB_LPR = TRUE ;
  else
    lastPd->ACEDB_LPR = FALSE ;

  lastPd->printers = getPrinters(lastPd->h) ;
  if (arrayMax(lastPd->printers))
    lastPd->doPrint = TRUE ;
  else
    lastPd->doMail = TRUE ; /* without attached printers: mail=default */

  if (lastPd->ACEDB_LPR)
    strcpy (lastPd->printerBuffer, arr(lastPd->printers, 0, char*)) ;
  else if (arrayMax(lastPd->printers))
    sprintf (lastPd->printerBuffer, "lpr -P%s",
	     arr(lastPd->printers, 0, char*)) ;
  else				/* default if everything fails */
    strcpy (lastPd->printerBuffer, "lpr") ;
  lastPd->printerNum = 0 ;

  if (filName ("PS", "", "wd"))
    strcpy (lastPd->dirBuffer, filName ("PS", "", "wd")) ;
  else if (getenv("PWD"))
    strcpy (lastPd->dirBuffer, getenv("PWD")) ;
  else if (!getcwd (lastPd->dirBuffer, MAXPATHLEN))
    strcpy (lastPd->dirBuffer, ".") ;

  return lastPd ;
}

void graphPrint (void) 
{
  PD pd =0 ;
  char *cp, fb[FIL_BUFFER_SIZE] ;
  Graph targetGraph ;
  char *targetGraphName ;
  static int nn = 1 ;
  
  if (lastPd && pdGraph)
    { 
      messout ("Print selector is already open on another graph, sorry") ;
      return ;
    }
  
  strcpy (fb, (cp = graphHelp(0)) ? cp : "acedb" ) ;
  
  targetGraph = gActive->id ;
  targetGraphName = gActive->name ;
  
  if (!graphExists (targetGraph))
    return ;
  
  if (!graphExists(pdGraph))
    { 
      /* calculate the size of the window, depending on the printers
	 that are all drawn to fit on the screen, we load the printers
	 and increase the window size by 20 pixels for each additional
	 (more than 1) line of printer buttons, the printer buttons
	 are drawn into a window 196 pixels wide (from character 
	 position 25 to 62 on the text window) */
      AC_HANDLE h = handleCreate () ;
      Array tmpPr  = getPrinters(h) ;
      float line = 570, x=0, xmax=296;
      int n=0, len ;
      
      while (n < arrayMax(tmpPr))
	{ 
	  len = strlen (getenv ("ACEDB_LPR") && n == 0 ?
			"ACE default" :	arr(tmpPr, n, char*)) ;
	  if (x + 8*len + 15 > xmax && x != 0)
	    { x = 0 ; line += 20 ; }
	  x += 8*len + 15 ; 
	  ++n ;
	}
      messfree (h) ;
      
      pdGraph = graphCreate (TEXT_FIT,
			     targetGraphName == 0 ? "Print/Mail Selection" :
			     messprintf ("Print/Mail Selection  ( %s )", targetGraphName),
			     .4, .25, 500.0/900.0, (line)/900.0) ;
      if (!pdGraph)
	return ;

      graphSetBlockMode(GRAPH_BLOCKING) ;		    /* User must answer print dialog. */
      
      graphRegister (DESTROY, pdDestroy) ;
      graphRegister (RESIZE, pdDraw) ;
      graphRegister (PICK, pdPick) ;
      
      /* ACEDB-GRAPH INTERFACE: if acedb has registered a different print    */
      /* menu use that, otherwise use default internal one.                  */
      if (getGraphAcedbPDMenu() != NULL) graphMenu(getGraphAcedbPDMenu()) ;
      else graphMenu (pdMenu) ;
      
      graphHelp ("Printer_Selection") ;
      
      pd = pdInit () ;
      pd->graph = pdGraph ;
      graphAssociate ((void*)PD_MAGIC, pd) ;
    }
  else
    {
      graphActivate (pdGraph) ;
      graphPop () ;
    }
  
  strcat (fb, messprintf(".%d", nn++)) ;
  strcpy (pd->filBuffer, fb) ;
  sprintf (pd->copyBuffer, "%s/%s.%s", 
	   pd->dirBuffer, pd->filBuffer,
	   (pd->style != 'a') ? "ps" : "txt") ;
  pd->target = targetGraph ;
  memset (pd->titleBuffer, 0, 80) ;
  
  graphActivate (pd->target) ;
  graphPSdefaults (&pd->isRotate, &pd->scale, &pd->pages) ;
  graphActivate (pd->graph) ;
  
  pdDraw () ;
  
  
  /* ACEDB-GRAPH INTERFACE: If an acedb user message function is registered  */
  /* then call it.                                                           */
  if (getGraphAcedbMainActivity() != NULL)
    (getGraphAcedbMainActivity())("Please Finish Print") ; /* mhmp 12.11.97*/
  
  
  graphLoop (TRUE) ;		/* blocking */
  
}

void graphBoundsPrint (float uw, float uh, void (*drawFunc)(void))
{
  int olduw = gActive->uw ;
  int oldw = gActive->w ;
  int olduh = gActive->uh ;
  int oldh = gActive->h ;
  
  gActive->uw = uw ; gActive->w = gActive->uw * gActive->xFac ;
  gActive->uh = uh ; gActive->h = gActive->uh * gActive->yFac ;
  
  if (drawFunc)
    (*drawFunc) () ;
  
  graphPrint () ;
  
  gActive->uw = olduw ;
  gActive->w = oldw ;
  gActive->uh = olduh ;
  gActive->h = oldh ;
}

/*********************** end of file ***************************/
 
 
 
 
 
