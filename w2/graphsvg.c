/*  File: svg.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@ncbi.nlm.nih.gov
 *
 * Description: Draws an acedb image for the web in svg format
 *              
 * Exported functions: 
 *      svgGraphExport : active graph exported as a full fledged html document
 * HISTORY:
 * Created Sept 3 2004 (mieg)
 *
 * CVS info:  
 *-------------------------------------------------------------------
 */
#include "regular.h"

#include "graph_.h"
#include "aceio.h" 
#include "whooks/systags.h"
#include "colorRGB.h"


static float XFac, YFac ;
#define MAXCOL 16

static BOOL DEBUG = FALSE ;

typedef struct { int x, y ; } svgPOINT ;

static char *svgColour (char *buf, int colour) ;

/*******************************************************************************/

static char *svgColour (char *buf, int colour)
{
  if (colour >= 0 && colour < NUM_TRUECOLORS)
    sprintf (buf, "rgb(%d,%d,%d)"
	     , colorRGB [3 * colour + 0]
	     , colorRGB [3 * colour + 1]
	     , colorRGB [3 * colour + 2]
	     ) ;
  else
    strcpy (buf, "svg(255, 255, 255)") ;
  
  return buf ;
} /* svgColour */

/*******************************************************************************/

static int svgSetLineWidth (ACEOUT out, int new)
{
  static int old = -1 ;
  if (0 && old != new)
    aceOutf (out, " ???????? \n", new) ;
  old = new ;
  
  return old ;
} /* svgSetLineWidth */

/*******************************************************************************/

static int svgSetTextHeight (ACEOUT out, int new)
{
  static int old = -1 ;
  if (0 && old != new)
    aceOutf (out, " ???????? \n", new) ;
  old = new ;
  
  return old ;
} /* svgSetTextHeight */

/*******************************************************************************/

static void svgSetFormat (ACEOUT out, int new)
{
  static int old = -1 ;
  if (0 && old != new)
    aceOutf (out, " ????????? %d \n", new) ;
  old = new ;
}

/*******************************************************************************/

static void svgFillBox (ACEOUT out, Box box)
{
  int  x1,x2,y1,y2 ;
  char bufBC[128], bufC[128] ;
  
  if (box->bcol && box->bcol < NUM_TRUECOLORS)
    {
      x1 = XFac * box->x1 ;
      x2 = XFac * box->x2 ;
      y1 = YFac * box->y1 ;
      y2 = YFac * box->y2 ;
      aceOutf (out, "<rect x=\"%d\" y=\"%d\" height=\"%d\" width=\"%d\" style=\"stroke:%s;fill:%s;stroke-width:%d\"/>\n"
	       , x1, y1, y2 - y1, x2 - x1
	       , svgColour (bufC, box->fcol)
	       , svgColour (bufBC, box->bcol)
	       , 0
	       ) ;
    }
  return ;
}

/*******************************************************************************/

static BOOL svgAction (ACEOUT out, Box box, int iBox)
{
  BUBBLEINFO *bubble = gActive->bubbleInfo && gActive->bubbleDict ? arrayp(gActive->bubbleInfo, iBox, BUBBLEINFO) : 0 ;  
  KEY key = bubble ? bubble->key : 0 ;
  BOOL hasAction = FALSE ;

  if (key)
    {
      if (0 && key == 1)
	aceOutf (out, " onclick=\"_root.OpenAnchor(0, \"%s\")\n\""
			      , dictName (gActive->bubbleDict, bubble->fName)
			      ) ;
      else if (0)
	aceOutf (out, " onclick=\"alert('%s')\" ", name(key)) ;
      else if (1)
	{
	  aceOutf (out, " onclick=\"javascript:openAceViewLink('%s','%s'); \"  "
		   , className(key), name(key) /* dictName (gActive->bubbleDict, bubble->fName) */
		   ) ;
	  hasAction = TRUE ;
	}
      else 
	{
	  aceOutf (out, " ><a href=\"javascript:openAceViewLink('%s','%s'); \"  "
		   , className(key), name(key) /* dictName (gActive->bubbleDict, bubble->fName) */
		   ) ;
	  hasAction = TRUE ;
	}
    }

  return hasAction ;
}  /* svgAction */

/*******************************************************************************/

static void svgBubble (ACEOUT out, Box box, int iBox)
{
  BUBBLEINFO *bubble = gActive->bubbleInfo && gActive->bubbleDict ? arrayp(gActive->bubbleInfo, iBox, BUBBLEINFO) : 0 ;  
  const char *bubbleInfo = bubble && bubble->bName ? dictName (gActive->bubbleDict, bubble->bName) : 0 ;
  Box box0 = gBoxGet (0) ;
  int xMax = XFac * box0->x2 - 4 ;
  static int kk = 0 ;

  if (bubbleInfo)
    {
      int x1 = XFac * box->x1 ;
      int x2 = XFac * box->x2 ;
      int y1 = YFac * box->y1 ;
      int y2 = YFac * box->y2 ;

      if (x1 > x2) { int dummy = x1 ; x1 = x2 ; x2 = dummy ; }
      if (y1 > y2) { int dummy = y1 ; y1 = y2 ; y2 = dummy ; }
      aceOutf (out, "<rect x=\"%d\" y=\"%d\" height=\"%d\" width=\"%d\" opacity=0 \" />\n"
	       , x1, y1, y2 - y1, x2 - x1
	       ) ;
    }
  if (bubbleInfo)
    {
      int dx = 2, dy = -1 ; 
      int x1 = XFac * ((box->x1 + box->x2)/2.0 + dx) ;
      int y1 = YFac * ((box->y1 + box->y2)/2.0 + dy) ;
      int lineHeight = .8 * YFac ;
      int ddx = 8 ;
      int ddy = 10 ;
      int width = 0 ;
      int height = YFac + 2 ;
      int nLines =  0 ;
      char *buf = strnew (bubbleInfo, 0) ;
      char *cq, *cp = buf ;

      if (0) buf[0] += (kk++) % 6 ;
      while (cp)
	{         /*  count the lines in the bubble */
	  int w ;
	  cq = strchr (cp, '\n') ;
	  if (cq)
	    { 
	      if (cq > cp)
		{
		  nLines++ ;
		  cq++ ;
		  y1 -= lineHeight ;
		  height += .94 * lineHeight ;
		  w = cq - cp ;
		}
	      cq++ ;
	    }
	  else
	    w = strlen(cp) ;
	  if (w > width)
	    width = w ;
	  cp = cq ;
	}

      if (y1 < height + ddy + 2)
	y1 = height + ddy + 2 ;

      width = XFac * 0.42 * width + 16 ;
      x1 -= width/2 + 0.4 * 3 * XFac ;

      if (x1 < 10)
	x1 = 10 ;
      if (x1 + width > xMax)
	x1 -= (x1 + width - xMax) ;
      aceOut (out, "<g  class='bubble'>\n") ;
      aceOutf (out, "<rect class='bubbleFrame' x=\"%d\" y=\"%d\" rx=\"3\" ry=\"3\"  height=\"%d\" width=\"%d\"  style=\"fill:#ffffef\"/>"
	       , x1 - ddx
	       , y1 - ddy 
	       , height 
	       , width + ddx
	       ) ;

      nLines =  0 ; cp = buf ;
      while (cp)
	{         /*  export each line */
	  cq = strchr (cp, '\n') ;
	  if (cq)
	    *cq++ = 0 ;
	  if (cp && *cp) /* !cq || cq > cp */
	    {
	      aceOutf (out, "<text  class='bubbleText' x='%d' y='%d'  ;>\n%s\n</text>\n"
		       , x1, y1 + lineHeight * nLines
		       , cp
		       ) ;				
	      nLines++ ;
	    }
	  cp = cq ;
	}

      aceOut (out, "</g>") ;
      messfree (buf) ;
    }
  return ;
}

/*******************************************************************************/
/* exports recursivelly the internal boxes
 */
static void svgDrawBox (ACEOUT out, Box box, int iBox, int draw)
{
  int    i1 ;
  int  x1,x2,y1,y2 ;
  int  colour = 0,  bgcolour = 0, lineWidth, format, currHeight ;
  int    action ;
  Box	 nextbox ;
  char bufBC[128], bufC[128] ;
 
  lineWidth = 1;
  colour = box->fcol ;
  bgcolour = box->bcol ;

  
  stackCursor (gStk, box->mark) ;
  
  if (iBox)
    {
      aceOutf (out, "<g class='bubbly' iBox='%d' ", iBox) ;
      if (!draw)
	svgAction (out, box, iBox) ;
      aceOut (out, ">\n") ;
    }
  while (!stackAtEnd (gStk))
    switch (action = stackNext (gStk, int))
      {
      case BOX_END:           /* exit point */
	if (! draw) svgBubble (out, box, iBox) ;
	aceOutf (out, "</g iBox='%d'>\n", iBox) ;
        return ;                        
      case BOX_START:          /* recursion */
        i1 = stackNext (gStk, int) ;
	if(DEBUG) fprintf(stderr,"BOX_START %d\n",i1);
	nextbox = gBoxGet (i1) ;
	if (draw) svgFillBox (out, nextbox) ;
	svgDrawBox (out, nextbox, i1, draw) ;           
	
	/* restore current box settings */
	svgSetLineWidth (out, lineWidth) ;
	svgSetFormat (out, format) ;
	break ;
      case COLOR:
        i1 = stackNext (gStk,int) ;
	colour = i1;
        break ;
      case TEXT_FORMAT:
        format = stackNext (gStk,int) ;
	svgSetFormat (out, format) ;
	/*********/
	if(DEBUG) fprintf(stderr,"TEXT_FORMAT %d\n", format);
        break ;
      case LINE_WIDTH:
	lineWidth = (int)(YFac * stackNext(gStk, float)) ;
	svgSetLineWidth (out, lineWidth) ;
        break ;
      case TEXT_HEIGHT:
	currHeight = stackNext (gStk,float) ;
	svgSetTextHeight (out, currHeight) ;
	/*********/
        break ;
      case POINT_SIZE:
        {
	  float psize = stackNext (gStk,float) ;
	  /*********/
	  if(DEBUG) fprintf(stderr,"POINT_SIZE %f\n",psize);
	}
        break ;
      case LINE: 
        x1 = (int)(XFac * stackNext (gStk,float));
        y1 = (int)(YFac * stackNext (gStk,float));
        x2 = (int)(XFac * stackNext (gStk,float));
        y2 = (int)(YFac * stackNext (gStk,float));
	if (draw)
	  aceOutf (out, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" style=\"stroke-width:%d;stroke:%s;\"/>\n"
		 , x1,y1,x2,y2
		 , lineWidth
		 , svgColour (bufC, colour)
		 );
	break ;
      case RECTANGLE: case FILL_RECTANGLE:
        x1 = (int)(XFac * stackNext (gStk,float));
        y1 = (int)(YFac * stackNext (gStk,float));
        x2 = (int)(XFac * stackNext (gStk,float));
        y2 = (int)(YFac * stackNext (gStk,float));
	if (x1 > x2) { int dummy = x1 ; x1 = x2 ; x2 = dummy ; }
	if (y1 > y2) { int dummy = y1 ; y1 = y2 ; y2 = dummy ; }
	if (draw)
	  switch (action)   /* fill:rgb(255,0,0) */
	    {
	    case RECTANGLE:
	      aceOutf (out, "<rect x=\"%d\" y=\"%d\" height=\"%d\" width=\"%d\" style=\"stroke-width:%d;stroke:%s;fill:%s; \" />\n"
		       , x1, y1, y2 - y1, x2 - x1
		       , lineWidth
		       , svgColour(bufC, colour)
		       , svgColour (bufBC, bgcolour)
		       ) ;
	      break ;
	    case FILL_RECTANGLE:
	      svgColour (bufC, colour) ;
	      aceOutf (out, "<rect x=\"%d\" y=\"%d\" height=\"%d\" width=\"%d\" style=\"stroke-width:%d;stroke:%s;fill:%s; \" />\n"
		       , x1, y1, y2 - y1, x2 - x1
		       , lineWidth
		       , bufC
		       , bufC
		       ) ;
	      break ;
	    }
	break ;
      case PIXELS: case PIXELS_RAW:
	{
	  unsigned char *pixels ;   /* consume the pixels */
	  int w, h, line;
	  
	  x1 = (int)stackNext (gStk,float) ;
	  y1 = (int)stackNext (gStk,float) ;
	  pixels = stackNext (gStk, unsigned char*) ;
	  if (0) fprintf (stderr, "%s", pixels) ; /* for compiler happiness */
	  w = stackNext (gStk, int) ;
	  h = stackNext (gStk, int) ;
	  line = stackNext (gStk, int) ; 
	  x2 = x2 + w + h - w - h + line - line  ;
	  if (action == PIXELS)
	    {
	      x2 = (int)stackNext (gStk,float) ;
	      y2 = (int)stackNext (gStk,float) ;
	    }
	  else
	    { 
	      x2 = x1 + XtoUrel(w) ;
	      y2 = x2 + YtoUrel(w) ;
	    }
	  
	}
	break ;
      case POLYGON: case LINE_SEGS:   /* consume the polygon */
	{ 
	  int n = stackNext (gStk, int) ;
	  if (n > 2) 
	    {
	      int i;
	      
	      if (draw)
		aceOutf (out, "<polyline style=\"stroke-width:%d;stroke:%s;fill:%s;\" points=\""
			 , lineWidth
			 , svgColour (bufC, colour) 
			 , svgColour (bufBC, bgcolour) 
			 ) ;
	      for(i=0;i<n;i++)
		{
		  x1 = (int)(XFac * stackNext (gStk, float));
		  y1 = (int)(YFac * stackNext (gStk, float));
		  if (draw)
		    aceOutf (out, " %d , %d ", x1, y1) ;
		}
	      if (draw)
		aceOutf (out, "\"/>\n") ;
	    }
	  break ;
	}
      case CIRCLE: case POINT: 
      case TEXT: case TEXT_UP: case TEXT_PTR: case TEXT_PTR_PTR: case TEXT_EXTERNAL:
      case COLOR_SQUARES: case FILL_ARC: case ARC:
        x1 = (int)(XFac * stackNext (gStk,float));
        y1 = (int)(YFac * stackNext (gStk,float));
        switch (action)
	  {
	    int r, a1, a2, dep ;
	    
	  case CIRCLE:
	    r = (int)(XFac * stackNext(gStk,float));
	    if (draw)
	      aceOutf (out, "<circle style=\"stroke-width:%d;stroke:%s;fill:%s\" cx=\"%d\" cy=\"%d\" r=\"%d\"/>"
		       , lineWidth
		       , svgColour (bufC, colour)  
		       , svgColour (bufBC, bgcolour)  
		       , x1,y1, r
		       ) ;
	    break ;
	  case FILL_ARC: case ARC:
	    r = (int)(2.0 * XFac * stackNext(gStk,float));
	    a1 = (int)stackNext(gStk,float);
	    a2 = (int)stackNext(gStk,float);
	    dep = a1 + a2;
	    if(dep < 0) dep += 360;
	    /* int arr = a1;
	       if (draw)
	       svgImageArc(out, im,x1,y1,r,r,dep,arr,alloca_colours[colour]);
	    */
	    break ;
	  case POINT:
	    { 
	      int psize2 = 1 ;
	      int pszx = (int)(XFac * psize2);
	      int pszy = (int)(YFac * psize2);
	      
	      if (draw)
		aceOutf (out, " ????? (%d, %d, %d, %d);\n", x1-pszx,y1-pszy, 2*pszx, 2*pszy) ;
	    }
	    break ;
	  case TEXT:
	  case TEXT_PTR:
	  case TEXT_PTR_PTR: 
	  case TEXT_EXTERNAL: 
	  case TEXT_UP:
	    { 
	      char *text ;
	      int len ;

	      switch (action)
		{ 
		case TEXT_UP: 
		  text = stackNextText (gStk) ; 
		  text = 0 ; /* i do not yet know how to display that */
		  break ;
		case TEXT: text = stackNextText (gStk) ; break ;
		case TEXT_PTR: text = stackNext (gStk, char *) ; break ;
		case TEXT_PTR_PTR: text = *stackNext(gStk, char **) ; break ;
		case TEXT_EXTERNAL: 
		  text = *stackNext(gStk, char **) ; 
		  len = *stackNext(gStk, int *) ; /* consume len */
		  x1 = x1 + len - len ;
		  break ;
		}
	      if (draw && text && *text)
		aceOutf (out, "<text x=\"%d\" y=\"%d\" fill=\"%s\">%s</text>\n", x1, y1+12, svgColour (bufC, colour), text) ;
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
		{ 
		  switch(*text^((*text)-1))
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
		  /* 
		     svgImageFilledRectangle (out, im,x1,y1,x2,y2,alloca_colours[col]); 
		     svgColour (bufC, colour) ;
		  */
		  text += iskip ;
		  x1 += XFac ; x2 += XFac + col - col ; /* compiler happiness */
		}
	    }
	    break ;
	  default:
	    if(DEBUG) fprintf(stderr,"Invalid action %d received in drawPSBox",action) ;
	    messout ("Invalid action %d received in svgDrawBox",action) ;
	    break ;     
	  }
      }

} /* svgDrawBox */

/*******************************************************************************/

static void svgDoBox (ACEOUT out, int iBox, int fcol, int bcol)
{
  Box box = gBoxGet (iBox) ;

  if (fcol >= 0)
    box->fcol = fcol ;
  if (bcol >= 0)
    box->bcol = bcol ;

  svgFillBox (out, box) ;	/* does background */
  if (1) svgDrawBox (out, box, iBox, 1) ;
  svgDrawBox (out, box, iBox, 0) ;
} /* svgDoBox (int box) */

/*******************************************************************************/
/* exports the outer box as a single java script
 */
static void svgOuterBox (ACEOUT out, int box)
{
  float x1, y1, x2, y2, w, h ; 
  
  XFac = gActive->xFac ;
  YFac = gActive->yFac ;
  
  graphBoxDim (0, &x1, &y1, &x2, &y2) ;
  w = XFac * (x2 - x1) ; h = YFac * (y2 - y1 + 4) ;
  if (gActive->h < h) h = gActive-> h ;
  if (gActive->w < w) w = gActive-> w ;
  
  aceOutf (out, "<svg class='ace_wide' viewBox='%f %f %f %f'  margin-top: 2px; >\n"
	   , x1, y1 - .05 * (y2 - y1), 1.1 * w , 1.1 * h ) ;

  svgDoBox (out, 0, -1, -1) ;  

  aceOut (out, "</svg>\n") ;
  return ;
} /* svgOuterBox */

/*******************************************************************************/
/* resize the SVG window before exportation */
BOOL svgGraphResize (float width, float height)
{
  if (1)
    {
      gActive->w = gActive->xFac * width ;
      gActive->h = gActive->yFac * height ;
    }
  gActive->uw = width ;
  gActive->uh = height ;
  graphFacMake() ;
  gUpdateBox0 () ;

  return TRUE ; 
}

/*******************************************************************************/
/* exports the active window using aceOut 
 * as a full fledged html document
 * return the number of bytes written
 */
BOOL svgGraphExport (Graph gId, ACEOUT out, BOOL do_size)
{
  aceOut (out, "<defs>\n"
	  "  <style type='text/css'>\n"
	  "    .ace_wide {width:100% ; height:auto ;}\n"
	  "    .bubbleText { font-size:8 ; stroke:none ; fill: black; text-align:center; position:relative;overflow-wrap: break-word;}\n"
	  "    .bubble {pointer-events: none; stroke:black ; opacity:0; }\n"
	  "    .bubbly {pointer-events: unset;  }\n"
	  "    .bubbly:hover .bubble { opacity:1;}  \n"
	  
	  "  </style>\n"
	  " </defs>\n"
	  ) ;

  svgOuterBox (out, 0) ;
  
  return TRUE ;
} /* svgActiveWindow */

/*
"<div class="speech-bubble">
    <p><strong>Demo speech bubble</strong></p>
    <p>This is a simple CSS speech bubble.</p>
</div>
*/
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
