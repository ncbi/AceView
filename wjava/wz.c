/*  File: wz.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@ncbi.nlm.nih.gov
 *
 * Description: Draws an acedb image for the web by invoking the 
 *              Walter Zorn java script drawing library 
 *              distributed as GPL and included in this directory
 *              
 * Exported functions: 
 *      wzGraph2html : active graph exported as a full fledged html document
 * HISTORY:
 * Created Sept 3 2004 (mieg)
 *
 * CVS info:  
 *-------------------------------------------------------------------
 */
#include "regular.h"

#include "graph_.h"
#include "freeout.h" 

#include "../wjava/wz.h"

static float XFac, YFac ;
#define MAXCOL 16
static char* wzCol[MAXCOL + 1] ;
static BOOL DEBUG = TRUE ;

typedef struct { int x, y ; } wzPOINT ;

/*******************************************************************************/

void wzColorInit (void)
{
  int i ;

  for (i = 0 ; i <= MAXCOL ; i++)
    wzCol [i] = "000000" ; /* black */
  wzCol [0] = "ffffff" ; /* white */
} /* wzColorInit */


/*******************************************************************************/

static void wzSetColour (int colour)
{
  static int lastColour = 0 ;

  if (colour < 0 || colour > MAXCOL) colour = 1 ;

  if (colour != lastColour)
    freeOutf ("  jg.setColor(\'#%s\');\n", wzCol [colour]) ; /* example 0000ff */
  lastColour = colour ;
}

/*******************************************************************************/

static void wzSetLineWidth (int new)
{
  static int old = -1 ;
  if (old != new)
    freeOutf ("  jg.setStroke (%d);\n", new) ;
  old = new ;
}

/*******************************************************************************/

static void wzSetFormat (int new)
{
  static int old = -1 ;
  if (0 && old != new)
    freeOutf ("  jg.setStroke (%d);\n", new) ;
  old = new ;
}

/*******************************************************************************/

static void wzFillBox (Box box)
{
  return ;
}

/*******************************************************************************/
/* exports recursivelly the internal boxes
 */
static void wzDrawBox (Box box)
{
  int    i1 ;
  int  x1,x2,y1,y2 ;
  int  colour,  lineWidth, format = 0 ;
  int    action ;
  Box	 nextbox ;
  
  lineWidth = -1;
  colour = box->fcol ;
  wzSetColour (colour) ;

  stackCursor (gStk,box->mark) ;

  while (!stackAtEnd (gStk))
    switch (action = stackNext (gStk, int))
      {
      case BOX_END:           /* exit point */
        return ;                        
      case BOX_START:          /* recursion */
        i1 = stackNext (gStk, int) ;
	if(DEBUG) fprintf(stderr,"BOX_START %d\n",i1);
	nextbox = gBoxGet (i1) ;
	wzFillBox (nextbox) ;
	wzDrawBox (nextbox) ;           

	/* restore current box settings */
	wzSetColour (colour) ;
	wzSetLineWidth (lineWidth) ;
	wzSetFormat (format) ;
	break ;
      case COLOR:
        i1 = stackNext (gStk,int) ;
	colour = i1;
	wzSetColour (colour) ;
        break ;
      case TEXT_FORMAT:
        format = stackNext (gStk,int) ;
	wzSetFormat (format) ;
	/*********/
	if(DEBUG) fprintf(stderr,"TEXT_FORMAT %d\n", format);
        break ;
      case LINE_WIDTH:
	lineWidth = (int)(YFac * stackNext(gStk, float)) ;
	wzSetLineWidth (lineWidth) ;
        break ;
      case TEXT_HEIGHT:
	stackSkip (gStk) ;  /* float currHeight */
	/*********/
        break ;
      case POINT_SIZE:
        {
	  float psize = stackNext (gStk,float) ;
	  /*********/
	  if(DEBUG) fprintf(stderr,"POINT_SIZE %f\n",psize);
	}
        break ;
      case LINE: case RECTANGLE: case FILL_RECTANGLE:
        x1 = (int)(XFac * stackNext (gStk,float));
        y1 = (int)(YFac * stackNext (gStk,float));
        x2 = (int)(XFac * stackNext (gStk,float));
        y2 = (int)(YFac * stackNext (gStk,float));
	switch (action)
	  {
	  case LINE:
	    if (x1 == x2 || y1 == y2)
	      freeOutf ("  jg.drawLine (%d, %d, %d, %d);\n", x1,y1,x2,y2);
	    break ;
	  case RECTANGLE:
	    freeOutf ("  jg.drawRect (%d, %d, %d, %d);\n", x1,y1,x2 - x1,y2 - y1);
	    break ;
	  case FILL_RECTANGLE:
	    freeOutf ("  jg.fillRect (%d, %d, %d, %d);\n", x1,y1,x2 - x1,y2 - y1);
	    break ;
	  }
	break ;
      case PIXELS: case PIXELS_RAW:
	{ 
	  int w ;

	  x1 = (int)stackNext (gStk,float) ;
	  y1 = (int)stackNext (gStk,float) ;
	  stackSkip (gStk) ;
	  w = stackNext (gStk, int) ;
	  stackSkip (gStk) ;
	  stackSkip (gStk) ;
	  if (action == PIXELS)
	    { x2 = (int)stackNext (gStk,float) ;
	      y2 = (int)stackNext (gStk,float) ;
	    }
	  else
	    { x2 = x1 + XtoUrel(w) ;
	      y2 = x2 + YtoUrel(w) ;
	    }
	  
	}
	break ;
      case POLYGON: case LINE_SEGS:   /* consume the polygon */
	{ 
	  int n = stackNext (gStk, int) ;
	  if (n > 2) {
	    wzPOINT *points, *p ;
	    int i;

	    points = (wzPOINT*) messalloc((n - 1) * sizeof(wzPOINT));
	    for(i=0,p=points;i<n-1;i++,p++)
	      {
		p->x = (int)(XFac * stackNext (gStk, float));
		p->y = (int)(YFac * stackNext (gStk, float));
	      }
	    /*
	      if (action == POLYGON)
	      gdImageFilledPolygon(im,points,n-1,alloca_colours[colour]);
	      gdImagePolygon(im,points,n-1,alloca_colours[colour]);
	    */
	    messfree(points);
	  }
	  break ;
	}
      case CIRCLE: case POINT: 
      case TEXT: case TEXT_UP: case TEXT_PTR: case TEXT_PTR_PTR: case TEXT_EXTERNAL:
      case COLOR_SQUARES: case FILL_ARC: case ARC:
        x1 = (int)(XFac * stackNext (gStk,float));
        y1 = (int)(YFac * stackNext (gStk,float));
        switch (action) {
	  int r, a1, a2, dep ;

	case CIRCLE:
	  r = (int)(XFac * stackNext(gStk,float));
	    freeOutf ("  jg.drawEllipse (%d, %d, %d, %d);\n", x1 - r, y1 - r, 2*r, 2*r) ;
	  break ;
	case FILL_ARC: case ARC:
	  r = (int)(2.0 * XFac * stackNext(gStk,float));
	  a1 = (int)stackNext(gStk,float);
	  a2 = (int)stackNext(gStk,float);
	  dep = a1 + a2;
	  if(dep < 0) dep += 360;
	  /*
	    int arr = a1;
	    gdImageArc(im,x1,y1,r,r,dep,arr,alloca_colours[colour]);
	  */
	  break ;
	case POINT:
	  {
	    int psize2 = 1 ;
	    int pszx = (int)(XFac * psize2);
	    int pszy = (int)(YFac * psize2);

	    freeOutf ("  jg.drawRect (%d, %d, %d, %d);\n", x1-pszx,y1-pszy, 2*pszx, 2*pszy) ;
	  }
	  break ;
	case TEXT:
	  {
	    char *text = stackNextText(gStk);
	    char *text2 = messalloc (2 * strlen (text) + 1) ;
	    char *cp = text - 1, *cq = text2 ;

	    while (*++cp) /* javaprotect ? */
	      switch (*cp)
		{
		case '\"':
		case '\'':
		case '\\':
		  *cq++ = '\\' ; *cq++ = *cp ; break ;
		default: *cq++ = *cp ; break ; 
	      }
	    *cq = 0 ;

	    freeOutf ("  jg.drawString (\"%s\", %d, %d);\n", text2, x1, y1) ;
	    messfree (text2) ;
	  }
	  break;
	case TEXT_UP:
	  { 
	    if (1)
	      stackSkip (gStk) ; /* remove the skip if you make this work */
#ifdef JUNK
	    else
	      {
		char *text = stackNextText(gStk);
		int n = strlen(text) ;
		char buf[2] ;
		
		buf[1] = 0 ;
		while (n--)
		  { 
		    buf[0] = *text++ ;
		    /*  im is not declared, who is this ?
		      gdImageString(im,gdFont6x9,x1+1,y1-9*(n+1),buf,
		      alloca_colours[colour]);
		    */
		  }
	      }
#endif
	  }
	  break;
	case TEXT_PTR:
	  { char *text = stackNext(gStk, char *);
	    if(DEBUG) fprintf(stderr,"TEXT_PTR %d %d %s\n",x1,y1,text);
	    /* gdImageString(im,gdFont8x13,x1,y1,text,alloca_colours[colour]); */
	  }
	  break;
	case TEXT_PTR_PTR:
	  { char *text = *stackNext(gStk, char **);
	    if (text && *text)
	      {
		if(DEBUG) fprintf(stderr,"TEXT_PTR_PTR %d %d %s\n",x1,y1,text);
		/* gdImageString(im,gdFont8x13,x1,y1,text,alloca_colours[colour]); */
	      }
	  }
	  break;
	case TEXT_EXTERNAL:
	  { char *text = *stackNext(gStk, char **);
            int len = *stackNext(gStk, int *);
	    if (text && *text)
	      {
		if(DEBUG) fprintf(stderr,"TEXT_EXTERNAL %d %d %s %d\n",x1,y1,text,len);
		/* gdImageString(im,gdFont8x13,x1,y1,text,alloca_colours[colour]); */
	      }
	  }
	  break;
	case COLOR_SQUARES:
	  {
	    char *text = stackNext (gStk, char*) ;
	    int n = stackNext (gStk, int) ;
	    int iskip = stackNext (gStk, int) ;
	    int *tints = stackNext (gStk,int*) ;
	    int col ;
	    if(DEBUG) fprintf(stderr,"COLOR_SQUARES %d %s\n",n,text);
	    x2 = x1 + XFac ; y2 = y1 + YFac ;
	    while (n--) 
	      { switch(*text^((*text)-1))
		    {
		    case 0x01: col = tints[0]; break;
		    case 0x03: col = tints[1]; break;
		    case 0x07: col = tints[2]; break;
		    case 0x0f: col = tints[3]; break;
		    case 0x1f: col = tints[4]; break;
		    case 0x3f: col = tints[5]; break;
		    case 0x7f: col = tints[6]; break;
		    case 0xff: col = tints[7]; break;
		    case -1:
		    default:col = WHITE; break;
		    }
	      /* gdImageFilledRectangle (im,x1,y1,x2,y2,alloca_colours[col]); */
		text += iskip ;
		x1 += XFac ; x2 += XFac ; 
		if (col) x1 += 0 ; /*to please the compiler */
	      }
	  }
	  break ;
	default:
	  if(DEBUG) fprintf(stderr,"Invalid action %d received in drawPSBox",action) ;
	  messout ("Invalid action %d received in wzDrawBox",action) ;
	  break ;     
	}
      }
} /* wzDrawBox */
 
/*******************************************************************************/

static void wzDoBox (int k, int fcol, int bcol)
{
  Box box = gBoxGet (k) ;

  /*
    currFormat = box->format ;
    currHeight = uToYrel(box->textheight) ;
    fontSet () ;
    psize = gActive->pointsize ;
    psize2 = psize/2 ;
  */
  if (fcol >= 0)
    box->fcol = fcol ;
  if (bcol >= 0)
    box->bcol = bcol ;

  wzColorInit () ;
  wzFillBox (box) ;	/* does background */
  wzDrawBox (box) ;
} /* wzDoBox (int box) */

/*******************************************************************************/
/* exports the outer box as a single java script
 */
static void wzOuterBox (int box)
{
  float x1, y1, x2, y2, w, h ; 

  XFac = gActive->xFac ;
  YFac = gActive->yFac ;

  graphBoxDim (0, &x1, &y1, &x2, &y2) ;
  w = XFac * (x2 - x1) ; h = YFac * (y2 - y1) ;
  if (gActive->h < h) h = gActive-> h ;
  if (gActive->w < w) w = gActive-> w ;

  freeOut ("<script language = \'JavaScript\'>\n"
	   "<!--\n"
	   "function myDrawFunction()\n"
	   "{\n"
	   ) ;

  wzDoBox (0, -1, -1) ;
  
  freeOut ("jg.paint(); // draws, in this case, directly into the document\n"
	   "}\n"
	   "\n"
	   "var jg = new jsGraphics(); // draw directly into document\n"
	   "myDrawFunction();\n"
	   "//-->\n"
	   "</script>\n"
	   ) ;
  return ;
} /* wzOuterBox */

/*******************************************************************************/
/* exports the active window using freeOut 
 * as a full fledged html document
 * return the number of bytes written
 */
BOOL wzGraph2html ()
{

  freeOut ("<html>\n"
	   "<head>\n"
	   "<script type=\'text/javascript\' src='jscript/wz_jsgraphics.js'></script>"
	   "</head>\n"
	   "<body>\n"
	   ) ;

  wzOuterBox (0) ;
	   
  freeOut ("</body></html>\n") ;

  return TRUE ;
} /* wzActiveWindow */
 
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

