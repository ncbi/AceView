/*  File: plot2d.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
     plot2d(char *title, Array a)
     where a is an Array of pairs of float.
 * HISTORY:
 * Last edited: Apr  5 16:19 1993 (mieg)
 * Created: Mon Apr  5 16:00:10 1993 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: plot2d.c,v 1.10 2017/03/18 15:30:51 mieg Exp $ */

#include <acedb.h>
#include <lex.h>
#include <regular.h>
#include <graph.h>
#include <plot.h>

static magic_t Plot2d_MAGIC = "Plot2d";

/*************************************************************************/

typedef struct p2dstuff { 
  magic_t *magic ;
  AC_HANDLE handle ;
  char title[FIL_BUFFER_SIZE], subtitleX[24], subtitleY[24] ;
  int titleBox, subtitleXbox, subtitleYbox ;
  Graph graph ;
  Array xy ;
  float xMin, xMax, yMin, yMax ;
  float x, y ;
  float xTrip, yTrip ;
  int tripBox ;
  KEY kTrip ;
  char xTripBuf [16], yTripBuf [16], kTripBuf [400] ;
  float xCentre, yCentre, xMag, yMag, oldXmag, oldYmag, startXmag ;
  int nx, ny ;
  int leftMargin, topMargin, bottomMargin, rightMargin ;
  int retBox, xyBox, regBox ;
  float xStart, yStart ;
  BOOL reticule, regress, isPrinting ;
  char xBuffer [24], yBuffer [24] ;
  double a, b, r, w ;
  float axisShift ;
} *P2D ;

static void plot2dZoomIn (void) ;
static void plot2dZoomOut (void) ;
static void plot2dWhole (void) ;
static void plot2dDraw (P2D p2d) ;
static void plot2dInit (P2D p2d) ;
static void plot2dReticule (int box, double x, double y) ;
static void plot2dCentre (double x, double y) ;
static void plot2dRegression (void) ;
static void plot2dPrint (void) ;
/*************************************************************************/

#define XGRAPH2MAP(x)   (p2d->xCentre  +  p2d->xMag * ((float)(x) - p2d->leftMargin - p2d->nx/2.))
#define XMAP2GRAPH(x)  (((x) - p2d->xCentre) / p2d->xMag + p2d->leftMargin + p2d->nx/2.)

#define YGRAPH2MAP(y) (p2d->yCentre  -  p2d->yMag * ((float)(y) - p2d->topMargin - p2d->ny/2.))
#define YMAP2GRAPH(y)  ((-(y) + p2d->yCentre) / p2d->yMag + p2d->topMargin + p2d->ny/2.)

/*************************************************************************/
/*************************************************************************/

static MENUOPT plot2dMenu[] =
  { 
    {graphDestroy, "Quit"},
    {help,"Help"},
    {plot2dPrint,"Print"},
    {plot2dZoomIn, "ZoomIn"},
    {plot2dZoomOut, "ZoomOut"},
    {plot2dWhole, "Whole"},
    {plot2dRegression, "Lin.Reg."},
    {0, 0}
  } ;

/*************************************************************************/

static P2D currentP2d (char *caller)
{
  P2D p2d = 0 ;

  if (!graphAssFind (&Plot2d_MAGIC,&p2d))
    messcrash("%s() could not find Plot2d on graph", caller);
  if (!p2d)
    messcrash("%s() received NULL  Plot2d pointer", caller);
  if (p2d->magic != &Plot2d_MAGIC)
    messcrash("%s() received non-magic pointer", caller);
  
  return p2d ;
} /* currentVerticalMap */

/*************************************************************************/
/*********** action routines *************/
/*************************************************************************/

static void plot2dResize (void)
{
  P2D p2d = currentP2d ("plot2dResize") ;

  plot2dDraw (p2d) ;
}

/***************************************/

static void plot2dDestroy(void)
{ 
  P2D p2d = currentP2d ("plot2dDestroy") ;
  
  if (p2d) 
    { 
      p2d->magic = 0 ; 
      arrayDestroy (p2d->xy) ;
      messfree (p2d->handle) ; 
    }
}
/******************************************************************/
/*mhmp 11.05.01  voir aussi plot.h */
/******************************************************************/
    /* returns 1, 2, 5, 10, 20, 50, 100 etc */
int regular(int p)
{
 register int i=1,j=1;
 if (!p)
   return 0 ;
 if(p<0) {j = -1 ; p = -p ; }
 while(i<p) i=10*i;
 i /= 10 ;
 if(j>0) {
   if(5*i>=p) return i*j;
   return 5*i*j;
 }
 else {
   if (5*i>=p) return 5*i*j ;
   return 10*i*j ;
 }
}

/********* start off with some utility routines ********/

static BOOL getNxNy(P2D p2d)
{
  int nx, ny ;
  graphFitBounds (&nx, &ny) ;
  nx = nx - p2d->leftMargin - p2d->rightMargin ;
  ny = ny - p2d->topMargin - p2d->bottomMargin ;
  if (nx > 0 && ny > 0)
    {
      p2d->nx = nx ; p2d->ny = ny ;
      return TRUE ;
    }
  p2d->nx = p2d->ny = 3 ;  /* arbitrary non stupid values */
  return FALSE ;
}

/***************************************/

static void plot2dWhole (void)
{ 
  BOOL regress ;

  P2D p2d = currentP2d ("plot2dWhole") ; 

  regress = p2d->regress ;
  plot2dInit (p2d) ;
  p2d->regress = regress ;
  plot2dDraw (p2d) ;
}

/***************************************/

static void plot2dZoomOut (void)
{ 
  P2D p2d = currentP2d ("plot2dZoomIn") ; 

  p2d->xMag = p2d->oldXmag ;
  p2d->yMag = p2d->oldYmag ;

  if (p2d->xMag < 4.* p2d->startXmag) {
    p2d->xMag *= 2 ;
    p2d->yMag *= 2 ;
  }
  
  plot2dDraw (p2d) ;
}

/***************************************/

static void plot2dZoomIn (void)
{ 
  P2D p2d = currentP2d ("plot2dZoomIn") ; 

  p2d->xMag = p2d->oldXmag ;
  p2d->yMag = p2d->oldYmag ;
  if (p2d->xMag > 1.E-5) {
    p2d->xMag /= 2 ;
    p2d->yMag /= 2 ;
  }
  plot2dDraw (p2d) ;
}
/***************************************/

static void plot2dPrint (void)
{ 
  P2D p2d = currentP2d ("plot2dPrint") ; 

  Graph old = graphActive();
  p2d->isPrinting = TRUE ;

  plot2dDraw (p2d) ;
  graphPrint () ;
  graphActivate (old);  
  p2d->isPrinting = FALSE ;
  plot2dDraw (p2d) ;
}

/***************************************/
static void plot2dReticule (int box, double x, double y)
{
  float x1, y1, x2, y2, xx, yy ;
  float dx, dy, dMin, d ;
  int i ;
  int ddx, ddy ;
  float width, height ;
  POINT2D *pp ;
  P2D p2d = currentP2d ("plot2dReticule") ;
  static KEY oldkey = 0 ;

  if (box == p2d->titleBox)
    { graphTextEntry (p2d->title, 0, 0, 0, 0) ;
      return ;
    }
  if (box == p2d->subtitleXbox)
    { graphTextEntry (p2d->subtitleX, 0, 0, 0, 0) ;
      return ;
    }
  if (box == p2d->subtitleYbox)
    { graphTextEntry (p2d->subtitleY, 0, 0, 0, 0) ;
      return ;
    }

  if (box != p2d->retBox){
    if (p2d->reticule) {
      p2d->xMag = p2d->oldXmag ;
      p2d->yMag = p2d->oldYmag ;
      p2d->reticule = FALSE ;
      plot2dDraw (p2d) ;
    }
    return ;
  }
  graphBoxDim (p2d->retBox, &x1, &y1, &x2, &y2) ;
  p2d->x = XGRAPH2MAP(x + x1) ;
  p2d->y = YGRAPH2MAP(y + y1) ;
  graphTextInfo (&ddx, &ddy, &width, &height) ;
  sprintf (p2d->xBuffer, "x = %g", p2d->x) ;
  sprintf (p2d->yBuffer, "y = %g", p2d->y) ;

  /* recherche du point le plus proche */
  pp  = arrp(p2d->xy, 0, POINT2D) ;
  xx = XMAP2GRAPH(pp->x) -x1;
  yy = YMAP2GRAPH(pp->y) -y1;
  dx = (xx - x) * ddx ;
  dy = (yy - y) * ddy ;
  dMin = dx * dx + dy * dy ;
  p2d->xTrip = pp->x ;
  p2d->yTrip = pp->y ;
  p2d->kTrip = pp->k ;
  pp++ ;
  for (i=1; i< arrayMax(p2d->xy) ; i++, pp++)
    {
      xx = XMAP2GRAPH(pp->x) -x1;
      yy = YMAP2GRAPH(pp->y) -y1;
      dx = (xx - x) * ddx ;
      dy = (yy - y) * ddy ;
      d = dx * dx + dy * dy ;
      if (d < dMin){
	dMin = d ;
	p2d->xTrip = pp->x ;
	p2d->yTrip = pp->y ;
	p2d->kTrip = pp->k ;
      }
    }
  sprintf (p2d->xTripBuf, "x = %g", p2d->xTrip) ;
  sprintf (p2d->yTripBuf, "y = %g", p2d->yTrip) ;
  if (iskey(p2d->kTrip) == 2)
    {
      sprintf (p2d->kTripBuf, "%s", name(p2d->kTrip)) ;
      if (p2d->kTrip != oldkey)
	oldkey = 0 ;
    }
  p2d->reticule = TRUE ;
  p2d->xMag = p2d->oldXmag ;
  p2d->yMag = p2d->oldYmag ;
  plot2dDraw (p2d) ; 
  if (oldkey)
    display (oldkey, 0, TREE) ;
  oldkey = p2d->kTrip ;
}  
/***************************************/
static void plot2dCentre (double x, double y)
{
  float x1, y1, x2, y2 ;
  P2D p2d = currentP2d ("plot2dCentre") ;

  graphBoxDim (p2d->retBox, &x1, &y1, &x2, &y2) ;
  if ( x >= x1 && x <= x2 && y >= y1 && y <= y2) {
    p2d->xCentre = XGRAPH2MAP(x) ;
    p2d->yCentre = YGRAPH2MAP(y) ;
    p2d->xMag = p2d->oldXmag ;
    p2d->yMag = p2d->oldYmag ;
    plot2dDraw (p2d) ;
  }
  else return ;
} 
/**************************************************/
/**************************************************/

static void plot2dXscale (P2D p2d)
{

  float x, xp1, start, end, xMag, sdx ;
  float xPos, yPos, fin, shift ;
  double xx, iDeb, iFin, xxp1 ;
  int xDiv, dx, ddx, i, j, oldi, nx, ny ;
  int xxx ;
  int xMul, isScale, scale= 1 ;
  yPos = p2d->topMargin + p2d->ny ;
  xPos = (p2d->leftMargin + p2d->nx) * 0.75 ;

  start = XGRAPH2MAP(p2d->leftMargin) ;
  end = XGRAPH2MAP(p2d->leftMargin + p2d->nx) ;
  p2d->xStart = start ;
  
  graphFitBounds (&nx, &ny) ;

  xMul = 1. ;
  if (start == end) end = start + 1 ;
  while ((end - start) * xMul < 1.)
    xMul *= 10 ;
  start *= xMul ;
  end *= xMul ;
  if (fabs(end) > fabs(start))
    i = fabs(end) ;
  else 
    i = fabs(start) ;
  oldi = i ;
  j = 1 ;
  while (i /= 10) j++ ; /* count chars */
  xMag = (nx - p2d->leftMargin) /  (end - start);
  dx = 5  * j / xMag ;
  dx = utMainRoundPart(dx) ;
  if(dx <= 0) dx = 1;

  xDiv = 1 ;
  ddx = dx ;
   while (ddx && !(ddx % 10))
    {
      xDiv *= 10 ;
      ddx /= 10 ;
    }
  if (xDiv > 1)
    xDiv /= 10 ;

  if (xDiv > 1 )
    {
      i = oldi/ xDiv ;
      j = 1 ;
      while (i /= 10) j++ ; 
      dx = 5  * j / xMag ;
      dx = utMainRoundPart(dx) ;
      if(dx <= 0) dx = 1;
    }
  /* mhmp 04.05.98 */
  if (dx * xMag < 2.5)
    dx *= 2 ;
  isScale = 0 ;
  if (xDiv > 1 || xMul > 1){
    if (xDiv > xMul) {
      scale = xDiv / xMul ;
      isScale = 1 ;
    }
    else {
      scale = xMul / xDiv ;
      isScale = -1 ;
    }
  }

  if (isScale >= 0) /* mhmp 06.07.99 > ---> >= */
    graphText (messprintf ("scale x: 1/%d", scale), xPos, yPos + 3.0) ;
  else
    graphText (messprintf ("scale x: 1*%d", scale), xPos, yPos + 3.0) ;

  iDeb = regular(utArrondi(start)) ;
  iFin = end ;
/*recalage a gauche */
  if (end - iDeb > dx * 100000) 
    iDeb = start ;
    else
      for (xx = iDeb ; xx < end ; xx += dx)
	if (xx > start) {
	  iDeb =  xx - dx ;
	  iFin = 2*p2d->xCentre * xMul - iDeb ;
	  break ;
	}
  shift = p2d->axisShift ;
  fin = p2d->leftMargin + p2d->nx + shift ;
  p2d->oldXmag = p2d->xMag ;
  p2d->xMag = (p2d->xCentre - (float)iDeb/xMul) / (p2d->nx/2.) ;
  if (p2d->xMag < 1.E-5)
    return ;
  for (xx = iDeb ; xx <= iFin ; xx += dx)
    { 
      x = XMAP2GRAPH((float)xx /xMul) ;
      if (x >= p2d->leftMargin - 0.0001) {
	xxx = utArrondi((float)xx / xDiv + 0.0001);
	graphLine (x, yPos + 1.5, x, yPos + 0.5) ;
	graphText (messprintf ("%d",xxx), x, yPos + 2.0) ;
	/*tirets*/
	xxp1 = xx + dx ;
	xp1 = XMAP2GRAPH((float)xxp1 /xMul ) ;
	sdx = (xp1 - x) / 5. ;
	for (i = 0 ; i < 4 ; i++) {
	  x += sdx ;
	  if (x <= fin + 0.0001)
	    graphLine (x, yPos + 1., x, yPos + 0.5) ;
	}
      }
    }
  graphLine (p2d->leftMargin - shift/2, yPos + shift, 
	     fin, yPos + shift) ;

}

/**************************************************/

static void plot2dYscale (P2D p2d)
{

  float y, yp1, start, end, yMag, sdy ;
  float xPos, yPos, fin, shift ;
  double yy, iDeb, iFin, yyp1 ;
  int yDiv, dy, ddy, nx, ny, i ;
  int yyy ;
  int yMul, isScale, scale = 1 ;
  char *cp ;
  xPos = p2d->leftMargin ;
  yPos = p2d->topMargin + p2d->ny ;

  start = YGRAPH2MAP(p2d->topMargin + p2d->ny) ;
  p2d->yStart = start ;
  end = YGRAPH2MAP(p2d->topMargin) ;

  graphFitBounds (&nx, &ny) ;

  yMul = 1. ;
  if (start == end) end = start + 1 ;
  while ((end - start) * yMul < 1.)
    yMul *= 10 ;
  start *= yMul ;
  end *= yMul ;

  yMag = (ny - p2d->topMargin - p2d->bottomMargin) / (end - start) ;
  dy = 5. / yMag ;
  dy = utMainRoundPart (dy) ;
  if(dy <= 0) dy = 1;
  if (dy * yMag < 2.5)
    dy *= 2 ;
  yDiv = 1 ;
  ddy = dy ;
   while (ddy && !(ddy % 10))
    {
      yDiv *= 10 ;
      ddy /= 10 ;
    }
   if (yDiv > 1)
     yDiv /= 10 ;

  isScale = 0 ;
  if (yDiv > 1 || yMul > 1){
    if (yDiv > yMul) {
      scale = yDiv / yMul ;
      isScale = 1 ;
    }
    else {
      scale = yMul / yDiv ;
      isScale = -1 ;
    }
  }
  if (isScale >= 0) /* mhmp 06.07.99 > ---> >= */
    graphText (messprintf ("scale y: 1/%d", scale), 2,yPos + 3) ;
  else 
    graphText (messprintf ("scale y: 1*%d", scale), 2,yPos + 3) ;    

  iDeb = regular(utArrondi(start)) ;
  iFin = end ;
  /* recalage en bas */
  if (end - iDeb > dy * 100000)
    iDeb = start ;
    else
      for (yy = iDeb ; yy < end ; yy += dy)
	if (yy > start) {
	  iDeb = yy - dy ;
	  iFin = 2*p2d->yCentre * yMul - iDeb ;
	  break ;
	}
  shift = p2d->axisShift ;
  fin = p2d->topMargin - shift ;
  p2d->oldYmag = p2d->yMag ;
  p2d->yMag = (p2d->yCentre - (float)iDeb/yMul) / (p2d->ny/2.) ;
  if (p2d->yMag < 1.E-5)
    return ;
  iFin += dy ;
  for (yy = iDeb ; yy <= iFin ; yy += dy) 
    { 
      y = YMAP2GRAPH((float)yy / yMul) ;
      if ((y <= p2d->topMargin + p2d->ny + 0.0001) &&
	  (y >= fin)) {
	yyy = utArrondi((float)yy / yDiv + 0.0001);
	graphLine (xPos - 1.5, y, xPos - 0.5, y) ;
	cp = messprintf ("%d", yyy) ;
	graphText (cp, p2d->leftMargin - strlen(cp)-1.5,y-.5) ;
	/*tirets*/
	yyp1 = yy + dy ;
	yp1 = YMAP2GRAPH((float)yyp1 /yMul ) ;
	sdy = (yp1 - y) / 5. ;
	for (i = 0 ; i < 4 ; i++) {
	  y += sdy ;
	    if (y >= fin - 0.0001)
	    graphLine (xPos - 1., y, xPos - 0.5, y) ;
	}  
      }
    }
  graphLine (xPos - shift, fin, 
	     xPos - shift, p2d->topMargin + p2d->ny + shift/2) ;
}

/***************************************/
/* intersection de la droite de regression avec le cadre */
/* y = ax + b  topMargin (+ ny) leftMargin (+ nx) */

static void regressInit (P2D p2d, float *x1, float *x2, float *y1, float *y2, 
			 int *nbInterp)
{
  float x, y ;
  int nbInter ;

  nbInter = 0 ;
  x = XMAP2GRAPH((YGRAPH2MAP(p2d->topMargin) - p2d->b) / p2d->a );
  if (x >= p2d->leftMargin && x <= p2d->leftMargin + p2d->nx) {
    nbInter++ ;
    *x1 = x ;
    *y1 = p2d->topMargin ;
  }
  x = XMAP2GRAPH((YGRAPH2MAP(p2d->topMargin + p2d->ny) - p2d->b) / p2d->a) ;
  if (x >= p2d->leftMargin && x <= p2d->leftMargin + p2d->nx) {
    nbInter++ ;
    if (nbInter < 2) {
      *x1 = x ;
      *y1 = p2d->topMargin + p2d->ny ;
    }	
    else {
      *x2 = x ;
      *y2 = p2d->topMargin + p2d->ny;
      *nbInterp = nbInter ;
      return ;
    }
  }
  y = YMAP2GRAPH(p2d->a * XGRAPH2MAP(p2d->leftMargin) + p2d->b) ;
  if ( y >= p2d->topMargin && y <= p2d->topMargin + p2d->ny) {
    nbInter++ ;
    if (nbInter < 2) {
      *x1 = p2d->leftMargin ;
      *y1 = y ;
    }	
    else {
      *x2 = p2d->leftMargin ;
      *y2 = y ;
      *nbInterp = nbInter ;
      return ;
    }
  }   
  y = YMAP2GRAPH(p2d->a * XGRAPH2MAP(p2d->leftMargin + p2d->nx) + p2d->b) ;
  if ( y >= p2d->topMargin && y <= p2d->topMargin + p2d->ny) {
    nbInter++ ;
    *x2 = p2d->leftMargin + p2d->nx ;
    *y2 = y ;
    *nbInterp = nbInter ;
  }
}
/***************************************/

static void plot2dDraw (P2D p2d)
{
  int i, ii ; 
  POINT2D *pp ;
  int nbInter = 0 ;
  float x = 0, y = 0, eps = 0.1, shift;
  float xPos = 0, yPos = 0, x1, x2, y1 = 0, y2 = 0 ;
  graphClear () ;
  if (!getNxNy (p2d))
    { graphRedraw() ; return ; }
  xPos = (p2d->leftMargin + p2d->nx) * 0.75 ;
  yPos = p2d->topMargin + p2d->ny ;
  p2d->titleBox = graphTextEntry(p2d->title, 47, 2,  0.1, 0) ;
  /* mhmp 06.07.99 xPos <--> 2 Xbox <-> Ybox */
  p2d->subtitleXbox = graphTextEntry(p2d->subtitleX, 23, xPos,  yPos + 4., 0) ;
  p2d->subtitleYbox = graphTextEntry(p2d->subtitleY, 23, 2,  yPos + 4., 0) ;
  if (!p2d->isPrinting)
    graphButtons (plot2dMenu, 2, 1.8, p2d->nx) ;
  graphBoxStart () ;
  graphBoxEnd() ;
 
 /* les axes */
  plot2dXscale (p2d) ;
  if (p2d->xMag < 1.E-5)
    goto fin ;

  plot2dYscale (p2d) ;
  if (p2d->yMag < 1.E-5)
      goto fin ;

  /* le cadre */
  x1 = p2d->leftMargin ;
  x2 = p2d->leftMargin + p2d->nx ;
  p2d->retBox = graphBoxStart() ;
  shift = p2d->axisShift ;
  graphLine (x1 - shift/2, p2d->topMargin - shift, 
	     x2 + shift, p2d->topMargin - shift) ;
  graphLine (x2 + shift, p2d->topMargin - shift, 
	     x2 + shift, p2d->topMargin + p2d->ny + shift/2) ;

  /* les points */
  for (i = 0 , pp = arrp (p2d->xy, 0, POINT2D) ; 
       i < arrayMax (p2d->xy) ; i++, pp++)
    {
      x = XMAP2GRAPH(pp->x) ;
      y = YMAP2GRAPH(pp->y) ;
      if (x >= p2d->leftMargin - eps && x <= p2d->leftMargin + p2d->nx + eps &&
	  y >= p2d->topMargin - eps && y <= p2d->topMargin + p2d->ny + eps){
	if ((p2d->reticule) && (pp->x == p2d->xTrip) && 
	    (pp->y == p2d->yTrip)){
	  ii = graphBoxStart() ;
	  graphCircle (x, y, 0.02) ;
	  graphBoxEnd () ;
	  graphBoxDraw(ii,RED,TRANSPARENT) ;
	}
      else
	graphCircle (x, y, 0.02) ;
      }
    }
  graphBoxEnd() ;
  /* droite de regression mhmp 25.05.99 */
  if (p2d->regress) {
    regressInit (p2d, &x1, &x2, &y1, &y2, &nbInter) ;
    if (nbInter == 2) {
      graphColor (BLUE) ;
      graphLine (x1, y1, x2, y2) ;
      graphColor (BLACK) ;
    }
    p2d->regBox = graphBoxStart() ;
    graphColor (BLUE) ;
    graphText (messprintf("y = %g * x + %g", p2d->a, p2d->b),
	       2, 5) ;
    graphText (messprintf ("r = %g w = %g", p2d->r, p2d->w), 38, 5) ;
    graphColor (BLACK) ;
    graphBoxEnd() ;
  }

/* triplets */
  if (p2d->reticule) {
    p2d->tripBox = graphBoxStart () ;
    graphTextPtr (p2d->xTripBuf, 2, 4, 8) ;
    graphTextPtr (p2d->yTripBuf, 20, 4, 8) ;
    graphTextPtr (p2d->kTripBuf, 44, 4, 8) ;
    graphBoxEnd () ;
    graphBoxDraw(p2d->tripBox,RED,TRANSPARENT) ;
  }

  /* le reticule */
  if (p2d->reticule) {
    p2d->xyBox = graphBoxStart() ;
    graphTextPtr (p2d->xBuffer, 2,  3 , 8) ;
    graphTextPtr (p2d->yBuffer, 20,  3 , 8) ;
    graphBoxEnd () ;
    graphBoxMenu(p2d->xyBox, plot2dMenu);
    x = XMAP2GRAPH(p2d->x);
    y = YMAP2GRAPH(p2d->y) ;
    if (x >= p2d->leftMargin - eps && x <= p2d->leftMargin + p2d->nx + eps &&
	y >= p2d->topMargin - eps && y <= p2d->topMargin + p2d->ny + eps){
      graphLine (x - 1, y, x + 1, y) ;
      graphLine (x, y - 1, x, y + 1) ;
    }
  }

  graphRedraw() ;
  return ;
fin:
    { graphText ("You are a little zoom-zoom!", 2, 5) ;
      graphRedraw () ;
    }

}

/**************************************************************/

static void plot2dInit (P2D p2d)
{
  float l, xMin, xMax , yMin, yMax, scaley ;
  int i ; 
  POINT2D *pp ; 
  if (0) { /*test petites valeurs */
    pp = arrp(p2d->xy, 0, POINT2D) ;
    for (i=0; i< arrayMax(p2d->xy) ; i++, pp++)
      { 
	pp->x = 35 + pp->x / 10000  ;
	pp->y *= 100000 ;
      }
  }

  /* avoid numbers so latge that they cannot be converted to integers */
  scaley = 1 ;
 lao:
  pp  = arrp(p2d->xy, 0, POINT2D) ;
  xMin = xMax = pp->x ;
  yMin = yMax = pp->y/scaley ;
  
  for (i=0; i< arrayMax(p2d->xy) ; i++, pp++)
    { 
      pp->y /= scaley ;
      if (pp->x > xMax)
	xMax = pp->x ;
      if (pp->x < xMin)
	xMin = pp->x ;
      if (pp->y > yMax)
	yMax = pp->y ;
      if (pp->y < yMin)
	yMin = pp->y ;
    }
  if (yMax - yMin > 10000000)
    { scaley *= 1000000 ; goto lao ; }
  p2d->xMin = xMin ;
  p2d->xMax = xMax ;
  p2d->yMin = yMin ;
  p2d->yMax = yMax ;

  getNxNy(p2d) ;
    
  l = p2d->xMax - p2d->xMin ;
  p2d->xCentre = p2d->xMin + l/2.0 ;
  /*  p2d->xMag = l/ (.9 * p2d->nx) ;mhmp 29.04*/
  p2d->xMag = l/ (1.0 * p2d->nx) ;
  p2d->oldXmag = p2d->xMag ;
  p2d->startXmag = p2d->xMag ;
  l = p2d->yMax - p2d->yMin ;
  p2d->yCentre = p2d->yMin + l/2.0 ;
  /*  p2d->yMag = l/ (.9 * p2d->ny) ;mhmp 29.04*/
  p2d->yMag = l/ (1.0 * p2d->ny) ;
  p2d->oldYmag = p2d->yMag ;
  p2d->reticule = FALSE ;
  p2d->regress = FALSE ;
  p2d->isPrinting = FALSE ;
  p2d->axisShift = 0.5 ;
}

/*****************************************/

void plot2d (char *title, char *subtitleX, char *subtitleY, Array xy)
{
  P2D p2d = 0 ;
  AC_HANDLE hh = 0 ;

  if (!arrayExists (xy) || !arrayMax(xy) ||
      xy->size != sizeof (POINT2D))
    return ;
  hh = handleCreate () ;

  p2d = (P2D) halloc (sizeof (struct p2dstuff), hh) ;
  p2d->handle = hh ;
	     
  if (!title || !*title)
    title = "Plot" ;
  sprintf(p2d->title, "%s", title) ;
  sprintf(p2d->subtitleX, "%s", subtitleX) ;
  sprintf(p2d->subtitleY, "%s", subtitleY) ;

  p2d->xy = xy ;
  p2d->topMargin = 7 ;/* mhmp 17.03.99 was 5*/
  p2d->bottomMargin = 5 ;/* mhmp 21.05.99 was 4*/
  p2d->leftMargin = 10 ; /* mhmp 25.05.99 was 4 */
  p2d->rightMargin = 2 ;

  if (graphExists(p2d->graph))
    { graphActivate(p2d->graph) ;
      graphPop() ;
    }
  else 
    { 
      if (!graphCreate(TEXT_FIT,title,.2,.2,.6,.55))
	return ;
      graphRetitle (title) ;
      p2d->graph = graphActive() ;
      p2d->magic =  &Plot2d_MAGIC ;
      graphAssociate (&Plot2d_MAGIC,p2d) ;
      graphRegister (DESTROY,plot2dDestroy) ;
      graphRegister (RESIZE,plot2dResize) ;
      
      graphRegister (PICK, plot2dReticule) ;
      graphRegister (MIDDLE_DOWN, plot2dCentre) ;
      /*     graphRegister (MIDDLE_DRAG,  plot2dDrag) ;
      graphRegister (MIDDLE_UP,  plot2dUp) ;
      */
      graphMenu (plot2dMenu) ;
    }
  plot2dInit (p2d) ;  /* sets ->centre, ->mag */
  linearRegression (p2d->xy, &p2d->a, &p2d->b, &p2d->r, &p2d->w) ;
  plot2dDraw (p2d) ;	
}

static void plot2dRegression (void)
{		 
  P2D p2d = currentP2d ("plot2dRegression") ; 

  p2d->regress = !p2d->regress ;
  p2d->xMag = p2d->oldXmag ;
  p2d->yMag = p2d->oldYmag ;
  plot2dDraw (p2d) ;
}
/*************************************************************************/
/*************************************************************************/
