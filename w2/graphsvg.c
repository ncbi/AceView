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
 *      svgGraph2html : active graph exported as a full fledged html document
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

static float XFac, YFac ;
#define MAXCOL 16
static char* svgCol[MAXCOL + 1] ;
static BOOL DEBUG = FALSE ;

typedef struct { int x, y ; } svgPOINT ;

/*******************************************************************************/

void svgColorInit (void)
{
  int i ;

  for (i = 0 ; i <= MAXCOL ; i++)
    svgCol [i] = "000000" ; /* black */
  svgCol [0] = "ffffff" ; /* white */
} /* svgColorInit */

/*******************************************************************************/

static int  svgSetColour (ACEOUT out, int colour)
{
  static int lastColour ;
  int col ;

  if (colour < 0 || colour > MAXCOL) colour = 1 ;

  if (0)
    {aceOutf (out, "<stroke=\"%s\">", name(_WHITE + colour - WHITE)) ; /* example 0000ff */
    col = lastColour ;
    lastColour = colour ;
    }
  return col ;
}

/*******************************************************************************/

static int svgUnSetColour (ACEOUT out)
{
  if (0) aceOutf (out, "</stroke>") ;
  return 0 ;
}

/*******************************************************************************/

static void svgSetLineWidth (ACEOUT out, int new)
{
  static int old = -1 ;
  if (0 && old != new)
    aceOutf (out, " ???????? \n", new) ;
  old = new ;
}

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
  if (0) aceOutf (out, " ????????? \n") ;
  return ;
}

/*******************************************************************************/
/* exports recursivelly the internal boxes
 */
static void svgDrawBox (ACEOUT out, Box box)
{
  int    i1 ;
  int  x1,x2,y1,y2 ;
  int  colour = 0,  lineWidth, format, currHeight ;
  int    action ;
  Box	 nextbox ;
  
  lineWidth = 1;
  colour = box->fcol ;
  svgSetColour (out, colour) ;
  
  stackCursor (gStk, box->mark) ;
  
  while (!stackAtEnd (gStk))
    switch (action = stackNext (gStk, int))
      {
      case BOX_END:           /* exit point */
	aceOutf (out, "</g>\n") ;
        return ;                        
      case BOX_START:          /* recursion */
	aceOutf (out, "<g>\n") ;
        i1 = stackNext (gStk, int) ;
	if(DEBUG) fprintf(stderr,"BOX_START %d\n",i1);
	nextbox = gBoxGet (i1) ;
	svgFillBox (out, nextbox) ;
	svgDrawBox (out, nextbox) ;           
	
	/* restore current box settings */
	svgSetColour (out, colour) ;
	svgSetLineWidth (out, lineWidth) ;
	svgSetFormat (out, format) ;
	break ;
      case COLOR:
        i1 = stackNext (gStk,int) ;
	colour = i1;
	svgSetColour (out, colour) ;
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
	    aceOutf (out, "<svg:line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" style=\"stroke:rgb(99,99,99);stroke-width:%d\"/>\n"
		      , x1,y1,x2,y2, lineWidth);
	    break ;
	  case RECTANGLE:
	    aceOutf (out, "<svg:rect  stroke=\"%s\" x=\"%d\"  y=\"%d\" height=\"%d\" width=\"%d\" style=\"stroke-width:%d\"/>\n"
		     , name (colour - BLACK + _BLACK)
		     , x1, y1, y2 - y1, x2 - x1, lineWidth
		     ) ;
	    break ;
	  case FILL_RECTANGLE:
	    aceOutf (out, "<svg:rect fill=\"%s\"  stroke=\"%s\" x=\"%d\" y=\"%d\" height=\"%d\" width=\"%d\" style=\"stroke-width:%d\"/>\n"
		     , name (colour - BLACK + _BLACK)
		     , name (box->bcol - BLACK + _BLACK)
		     , x1, y1, y2 - y1, x2 - x1, lineWidth
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
	  w = stackNext (gStk, int) ;
	  h = stackNext (gStk, int) ;
	  line = stackNext (gStk, int) ;
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
	  if (n > 2) 
	    {
	      int x,y,i;
	      
	      aceOutf (out, "<svg:polyline fill=\"%s\" stroke=\"%s\" stroke-width=\"%d\" points=\""
		       , action == POLYGON ? name (box->bcol - BLACK + _BLACK) : "none"
		       , name (colour - BLACK + _BLACK)
		       , lineWidth
		       ) ;
	      for(i=0;i<n-1;i++)
		{
		  x1 = (int)(XFac * stackNext (gStk, float));
		  y1 = (int)(YFac * stackNext (gStk, float));
		  aceOutf (out, " %d , %d ", x1, y1) ;
		}
	      aceOutf (out, "\"/>\n") ;
	      
	      /*
		<rect fill="red" stroke="black" x="15" y="15" width="100" height="50"/>
		<rect fill="blue" stroke="black" x="150" y="15" width="100" height="50" rx="12" ry="18"/>
		<circle fill="yellow" stroke="black" cx="62" cy="135" r="20"/>
		<ellipse fill="green" stroke="black" cx="200" cy="135" rx="50" ry="20"/>
		<polyline fill="none" stroke="red" stroke-width="2" points="160,200 180,230 200,210 234,220"/>
		<polygon fill="yellow" stroke="black" stroke-width="1" points="49,272 57,297 "/>
	    */
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
	    int r, a1, a2, dep, arr;
	    
	  case CIRCLE:
	    r = (int)(XFac * stackNext(gStk,float));
	    aceOutf (out, "<svg:circle stroke=\"%s\" cx=\"%d\" cy=\"%d\" r=\"%d\"/>"
		     , name (colour - BLACK + _BLACK)
		     , x1,y1, r
		     ) ;
	    break ;
	  case FILL_ARC: case ARC:
	    r = (int)(2.0 * XFac * stackNext(gStk,float));
	    a1 = (int)stackNext(gStk,float);
	    a2 = (int)stackNext(gStk,float);
	    dep = a1 + a2;
	    arr = a1;
	    if(dep < 0) dep += 360;
	    /*
	      svgImageArc(out, im,x1,y1,r,r,dep,arr,alloca_colours[colour]);
	    */
	    break ;
	  case POINT:
	    {
	      int psize2 = 1 ;
	      int pszx = (int)(XFac * psize2);
	      int pszy = (int)(YFac * psize2);
	      
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
		  text = 0 ; /* i do not yet know how to diaply that */
		  break ;
		case TEXT: text = stackNextText (gStk) ; break ;
		case TEXT_PTR: text = stackNext (gStk, char *) ; break ;
		case TEXT_PTR_PTR: text = *stackNext(gStk, char **) ; break ;
		case TEXT_EXTERNAL: 
		  text = *stackNext(gStk, char **) ; 
		  len = *stackNext(gStk, int *) ; /* consume len */
		  break ;
		}
	      if (text && *text)
		{ 
		  char *cp, *cq, *text2 ;
		  
		  text2 = messalloc (2 * strlen (text) + 1) ;
		  cp = text - 1; cq = text2 ;
		  
		  while (*++cp) /* may be this is the same as the function javaprotect ? */
		    switch (*cp)
		      {
		      case '&': break ;   /* ?????? kills svg */
		      case '\"':
		      case '\'':
		      case '\\':
			*cq++ = '\\' ; *cq++ = *cp ; break ;
		      default: *cq++ = *cp ; break ; 
		      }
		  *cq = 0 ;
		  
		  aceOutf (out, "<svg:text x=\"%d\" y=\"%d\">%s</svg:text>\n", x1, y1, text2) ;
		  messfree (text2) ;
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
		  /* svgImageFilledRectangle (out, im,x1,y1,x2,y2,alloca_colours[col]); */
		  text += iskip ;
		  x1 += XFac ; x2 += XFac ;
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

static void svgDoBox (ACEOUT out, int k, int fcol, int bcol)
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

  if (0) svgColorInit () ;
  if (0) svgFillBox (out, box) ;	/* does background */
  svgDrawBox (out, box) ;
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
  w = XFac * (x2 - x1) ; h = YFac * (y2 - y1) ;
  if (gActive->h < h) h = gActive-> h ;
  if (gActive->w < w) w = gActive-> w ;

  svgDoBox (out, 0, -1, -1) ;  

  return ;
} /* svgOuterBox */

/*******************************************************************************/
/* exports the active window using aceOut 
 * as a full fledged html document
 * return the number of bytes written
 */
BOOL svgGraphExport (Graph gId, ACEOUT out, BOOL do_size)
{

  aceOut (out, "<?xml version=\"1.0\" encoding=\"iso-8859-1\" standalone=\"no\"?>\n"
	   "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
	  "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"

	   "<svg xml:space=\"preserve\" width=\"2000\" height=\"3000\"\n"
	   "xmlns=\"http://www.w3.org/2000/svg\">\n"
	  "<style type=\"text/css\">"
	  ".redbox{fill:#FF0000;}"  
	  ".whitewords{font-family:Times-Bold;font-size:36;fill:#FFFFFF;}"
	  "</style>"


	  "<g>"
	  "  <rect class=\"redbox\" x=\"100\" y=\"100\" width=\"460\" height=\"50\" /> \n"
	  "  <text class=\"whitewords\" x=\"120\" y=\"120\" >\n"
	  "    This site is powered by AceView/SVG.</text>"
	  "</g>"
	  ) ;
  /* 
	  "<svg width=\"100%\" height=\"100%\" version=\"1.1\"\n"
	  "xmlns=\"http://www.w3.org/2000/svg\">\n"

	   "<svg xml:space=\"preserve\" width=\"2000\" height=\"3000\"\n"
	   "xmlns=\"http://www.w3.org/2000/svg\">\n"
	  "<style type=\"text/css\">"
	  ".redbox{fill:#FF0000;}"  
	  ".whitewords{font-family:Times-Bold;font-size:36;fill:#FFFFFF;}"
	  "</style>"


	  "<g>"
	  "  <rect class=\"redbox\" x=\"10\" y=\"0\" width=\"460\" height=\"50\" /> \n"
	  "  <text class=\"whitewords\" x=\"20\" y=\"40\" >\n"
	  "    This site is powered by SVG.</text>"
	  "</g>"

 */

  svgOuterBox (out, 0) ;

  aceOut (out, "</svg>\n") ;

  return TRUE ;
} /* svgActiveWindow */
 
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

#ifdef JUNK

example of a full document

<?xml version="1.0" encoding="iso-8859-1"?>
2 <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 12August 1999//EN" "http://www.w3.org/Graphics/SVG/SVG-19990812.dtd">
3 <svg xml:space="preserve" width="1000" height="1000" >
4 <style type="text/css">
5 .redbox{fill:#FF0000;}   ///  define style 1
6 .whitewords{font-family:Times-Bold;font-size:36;fill:#FFFFFF;} ///  define style 2
7 </style>
8 <g>
9 <rect class="redbox" x="10" y="0" width="460" height="50" />   /// use style redbox
10 <text class="whitewords" x="20" y="40" >This site is powered by SVG.</text>
11 </g>
12 </svg>


OTHER example 
<?xml version="1.0" encoding="iso-8859-1" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/SVG/DTD/svg10.dtd">
<svg viewBox="0 0 270 400" width="100%" height="100%" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
	<g id="mainlayer">
		<rect fill="red" stroke="black" x="15" y="15" width="100" height="50"/>
		<rect fill="blue" stroke="black" x="150" y="15" width="100" height="50" rx="12" ry="18"/>
		<circle fill="yellow" stroke="black" cx="62" cy="135" r="20"/>
		<ellipse fill="green" stroke="black" cx="200" cy="135" rx="50" ry="20"/>
		<polyline fill="none" stroke="red" stroke-width="2" points="160,200 180,230 200,210 234,220"/>
		<polygon fill="yellow" stroke="black" stroke-width="1" points="49,272 57,297 114,282 63,314 71,339 49,324 27,339 35,314 13,297 40,297"/>

		<path stroke="black" fill="none" stroke-width="3" d="M150 280l19,10 -22,33 40,3c12,43 44,-83 83,20"/>
		<g stroke="green" fill="none">
			<line x1="15" y1="240" x2="30" y2="200" stroke-width="2"/>
			<line x1="30" y1="240" x2="45" y2="200" stroke-width="4"/>
			<line x1="45" y1="240" x2="60" y2="200" stroke-width="8"/>
			<line x1="60" y1="240" x2="75" y2="200" stroke-width="10"/>
			<line x1="75" y1="240" x2="90" y2="200" stroke-width="12"/>
		</g>
		<g font-size="20px">

			<text x="44" y="88">rect</text>
			<text x="140" y="88">rect (rounded)</text>
			<text x="36" y="180">circle</text>
			<text x="170" y="180">ellipse</text>
			<text x="36" y="263">line</text>
			<text x="156" y="255">polyline</text>

			<text x="16" y="363">polygon</text>
			<text x="140" y="363">path:<tspan x="140" y="383">simple + bezier</tspan></text>
		</g>
	</g>
</svg

#endif

