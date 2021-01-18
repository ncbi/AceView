/*  File: graphps.c
 *  Author: Richard Durbin (rd@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: to provide postscript output for the graph package
 * Exported functions: graphPS(), graphPrint()
 * SCCS: $Id: graphps.c,v 1.8 2020/05/30 16:50:30 mieg Exp $
 * HISTORY:
 * Last edited: Jan  8 11:42 1999 (edgrif)
 * * Jan  5 14:11 1999 (edgrif): Make graphPStitle() internal as PStitle,
 *              it's not used anywhere. Correct code for checking string
 *              pointers in several places (i.e. NULL ptr check AND null
 *              string check required).
 * * Sep 25 16:52 1998 (edgrif): Replace #ifdef ACEDB with graph/acedb,
 *              two references to session and one to acedb PS fonts file.
 * * Oct 1997 (michel): fixed the date printing
 * * Apr 10 (mieg): set for color laser printer (add COLOUR in psfont.wrm)
 * * Aug 14 18:05 1992 (rd): COLOR_SQUARES
 * * Jul 25 12:28 1992 (mieg): TEXT_PTR_PTR
 * Created: Wed May 20 08:30:44 1992 (rd)
 *-------------------------------------------------------------------
 */

#include "regular.h"
#include "call.h"
#include "graph_.h"
#include "mytime.h"		/* for timeShow() / Parse */


static float Grey[] =           /* must correspond to enum Colours */
  {
    1.0, 0.0, 0.7,0.3,		/* old values arranged to be (falsely) distinct */
    0.4,0.5,0.6, 
    0.9,0.7,0.6, 
    0.75,0.80,0.850, 
    0.15, 0.2, 0.25,

    0.9, 0.9, 0.9,
    0.95, 0.95, 0.95,
    0.4, 0.6, 0.8,
    0.5, 0.7, 0.9,
    0.55, 0.85, 0.5, 0.5
  } ;

static float Red[256], Green[256], Blue[256] ;

static FILE *fil ;
static float psize,psize2 ;
static BOOL isColor = FALSE ;

/*****************************************/
static void colorInit (void)
{
  int i, ncol = NUM_TRUECOLORS ;
  float grey ;
  static BOOL firstPass = TRUE ;
  
  if (!firstPass)
    return ;
  if (ncol > 255)
    messcrash ("too many color in graphgd.c:colorInit ") ;
  for (i = 0 ; i < ncol ; i++)
    { 
      graphGetColourAsFloat (i, &(Red[i]), &(Green[i]), &(Blue[i]), &grey) ;
    }
  firstPass = FALSE ;
}

/*****************************************/

static float fontHeight ;

static float currHeight ;
static int currFormat ;

static void fontSet ()
{
  switch (currFormat)
    {
    case ITALIC:
      if (currHeight)
	fprintf (fil, "italicFont %.4f scalefont setfont\n", 
		 currHeight) ;
      else
	fprintf (fil, "italicFont setfont\n") ;
      break ;
    case BOLD:
      if (currHeight)
	fprintf (fil, "boldFont %.4f scalefont setfont\n", 
		 currHeight) ;
      else
	fprintf (fil, "boldFont setfont\n") ;
      break ;
    case FIXED_WIDTH:
      if (currHeight)
	fprintf (fil, "fixedFont %.4f scalefont setfont\n", 
		 currHeight) ;
      else
	fprintf (fil, "fixedFont setfont\n") ;
      break ;
    case PLAIN_FORMAT:
    default:
      if (currHeight)
	fprintf (fil, "defaultFont %.4f scalefont setfont\n", 
		 currHeight) ;
      else
	fprintf (fil, "defaultFont setfont\n") ;
      break ;
    }
}

/***********************************************/

static void fillBox (Box box)   /* does background */

{ 
  if (box->x1 > 999999.0)
    return ;                    /* box still empty */

 if (isColor)
   { if (box->bcol != TRANSPARENT)
       { fprintf (fil, "%.2f %.2f %.2f setrgbcolor\n", 
		  Red[box->bcol], Green[box->bcol], Blue[box->bcol]) ;
	 fprintf (fil, "%.4f %.4f %.4f %.4f rectangle fill\n",
		  box->x1,box->y1,box->x2,box->y2) ;
       }
     fprintf (fil, "%.2f %.2f %.2f setrgbcolor\n", 
		  Red[box->fcol], Green[box->fcol], Blue[box->fcol]) ;
   }
  else
   { if (box->bcol != TRANSPARENT)
       { fprintf (fil, "%.2f setgray\n", Grey[box->bcol]) ;
	 fprintf (fil, "%.4f %.4f %.4f %.4f rectangle fill\n",
		  box->x1,box->y1,box->x2,box->y2) ;
       }
     fprintf (fil, "%.2f setgray\n",Grey[box->fcol]) ;
   }
}

/***********************************************/

static void drawBox (Box box)
{
  int    i1 ;
  float  x1,x2,y1,y2 ;
  int    action ;
  Box	 nextbox ;
  float  colour , red, green, blue ;

  colour = Grey[box->fcol] ;
  red = Red[box->fcol] ;
  green = Green[box->fcol] ;
  blue = Blue[box->fcol] ;

  stackCursor (gStk,box->mark) ;

  while (!stackAtEnd (gStk))
    switch (action = stackNext (gStk, int))
      {
      case BOX_END:
        return ;                        /* exit point */
      case BOX_START:
        i1 = stackNext (gStk, int) ;
	nextbox = gBoxGet (i1) ;
	fillBox (nextbox) ;
        drawBox (nextbox) ;                  /* recursion */
	if (isColor)
	  fprintf (fil, "%.2f %.2f %.2f setrgbcolor\n", 
		   red, green, blue ) ;
	else
	  fprintf (fil, "%.2f setgray\n",colour) ;
	break ;
      case COLOR:
        i1 = stackNext (gStk,int) ;
	if (i1==FORECOLOR) 
	  i1=box->fcol; 
	else if (i1==BACKCOLOR)
	  i1=box->bcol;
	colour = Grey[i1] ;
	red = Red[i1] ;
	green = Green[i1] ;
	blue = Blue[i1] ;
	if (isColor)
	  fprintf (fil, "%.2f %.2f %.2f setrgbcolor\n", red, green, blue) ;
	else
	  fprintf (fil, "%.2f setgray\n", colour) ;
        break ;
      case TEXT_FORMAT:
        currFormat = stackNext (gStk,int) ;
	fontSet () ;
        break ;
      case LINE_WIDTH:
        fprintf (fil,"%f setlinewidth\n",stackNext (gStk,float)) ;
        break ;
      case TEXT_HEIGHT:
	currHeight = stackNext (gStk,float) ;
        fontSet () ;
        break ;
      case POINT_SIZE:
        psize = stackNext (gStk,float) ; psize2 = psize/2 ;
        break ;
      case LINE: case RECTANGLE: case FILL_RECTANGLE:
        x1 = stackNext (gStk,float) ;
        y1 = stackNext (gStk,float) ;
        x2 = stackNext (gStk,float) ;
        y2 = stackNext (gStk,float) ;
        switch (action)
          {
          case LINE:
            fprintf (fil,
                "newpath  %.4f %.4f moveto  %.4f %.4f lineto  stroke\n",
                     x1,y1,x2,y2) ;
            break ;
          case RECTANGLE:
            fprintf (fil,"%.4f %.4f %.4f %.4f rectangle stroke\n",
                     x1,y1,x2,y2) ;
            break ;
          case FILL_RECTANGLE:
            fprintf (fil,"%.4f %.4f %.4f %.4f rectangle fill\n",
                     x1,y1,x2,y2) ;
            break ;
          }
	break ;
      case PIXELS: case PIXELS_RAW:
	{ unsigned char *ip, *pixels, hexVal[256] ; 
	  int w, h, line, i, j, reverse[256] ;
	  float fac ;

	  x1 = stackNext (gStk, float) ;
	  y1 = stackNext (gStk, float) ;
	  pixels = stackNext (gStk, unsigned char*) ;
	  w = stackNext (gStk, int) ;
	  h = stackNext (gStk, int) ;
	  line = stackNext (gStk, int) ;
	  if (action == PIXELS)
	    { x2 = stackNext (gStk, float) ;
	      y2 = stackNext (gStk, float) ;
	      for (i = 0 ; i < 128 ; ++i)
		reverse[i] = i ;
	    }
	  else
	    { x2 = x1 + XtoUrel(w) ;
	      y2 = y1 + YtoUrel(h) ;
	      graphRawMaps (0, reverse) ;
	    }

	  if (rampMin < rampMax)
	    fac = 0xff / (float) (rampMax - rampMin) ;
	  else
	    fac = 0xff / (float) (rampMin - rampMax) ;

          for (i = 0 ; i < 256 ; ++i)
	    { j = reverse[i] ;
	      if (rampMin < rampMax)
		{
		  if (j < rampMin)
		    hexVal[i] = 0 ;
		  else if (j >= rampMax)
		    hexVal[i] = 0xff ;
		  else
		    hexVal[i] = fac * (j - rampMin) ;
		}
	      else
		{
		  if (j < rampMax)
		    hexVal[i] = 0xff ;
		  else if (j >= rampMin)
		    hexVal[i] = 0 ;
		  else
		    hexVal[i] = fac * (rampMin - j) ;
		}
	    }
	  
	  fprintf (fil, "/picstring %d string def\n", w) ;
	  fprintf (fil, "  gsave") ;
	  fprintf (fil, " %.4f %.4f translate ", x1, y1) ;
	  fprintf (fil, " %.4f %.4f scale\n", x2-x1, y2-y1) ;
	  fprintf (fil, "  %d %d 8", w, h) ;
	  fprintf (fil, " [ %d 0 0 %d 0 0 ]", w, h) ;
	  fprintf (fil, " { currentfile picstring"
		   " readhexstring pop } image") ;
	  ip = pixels ;
	  for (j = 0 ; j < h ; j++)
	    { fprintf (fil, "\n  ") ;
	      for (i = 0 ; i < w ; )
		{ fprintf (fil, "%02x", hexVal[*ip++]) ;
		  if (!(++i % 35))
		    fprintf (fil, "\n    ") ;
		}
	      if (i < line)
		ip += (line-i) ;
	    }
	  fprintf (fil, "  grestore\n") ;
	}
	break ;
      case POLYGON: case LINE_SEGS:
	{ int n = stackNext (gStk, int) ;
	  if (n-- > 2)
	    { 
	      float x, y ;
	      x = stackNext (gStk, float) ;
	      y = stackNext (gStk, float) ;
	      fprintf (fil,   "newpath %.4f %.4f moveto\n", x, y) ;
	      while (n--)
		{ 
		  x = stackNext (gStk, float) ;
		  y = stackNext (gStk, float) ;
		  fprintf (fil, "        %.4f %.4f lineto\n", x, y) ;
		} 
	      if (action == POLYGON)
		fprintf (fil,   "        fill\n") ;
	      else
		fprintf (fil,   "        stroke\n") ;
	    }
	  break ;
	}
      case CIRCLE: case POINT: 
      case TEXT: case TEXT_UP: case TEXT_PTR: case TEXT_PTR_PTR: case TEXT_EXTERNAL:
      case COLOR_SQUARES: case FILL_ARC: case ARC:
        x1 = stackNext (gStk,float) ;
        y1 = stackNext (gStk,float) ;
        switch (action)

          {
          case CIRCLE:
            fprintf (fil,"%.4f %.4f moveto\n", x1, y1) ;
	    fprintf (fil,"newpath %.4f %.4f %.4f 0 360 arc stroke\n",
                     x1, y1, stackNext (gStk,float)) ;
            break ;
	  case ARC:
	    { float r = stackNext (gStk, float) ;
	      float ang = stackNext (gStk, float) ;
	      float angdiff = stackNext (gStk, float) ;
	      fprintf (fil,"newpath %.4f %.4f %.4f %.4f %.4f arc stroke\n",
		       x1, y1, r, -ang-angdiff, -ang) ;
	    }
            break ;
	  case FILL_ARC:
	    { float r = stackNext (gStk, float) ;
	      float ang = stackNext (gStk, float) ;
	      float angdiff = stackNext (gStk, float) ;
	      fprintf (fil,"newpath %.4f %.4f moveto\n", x1, y1) ;
	      fprintf (fil,"%.4f %.4f %.4f %.4f %.4f arc fill\n",
		       x1, y1, r, -ang-angdiff, -ang) ;
	    }
            break ;
          case POINT:
            fprintf (fil,"%.4f %.4f %.4f %.4f rectangle fill\n",
                     x1-psize2,y1-psize2,x1+psize2,y1+psize2) ;
            break ;
          case TEXT:
	    fprintf (fil,"%.4f %.4f moveto (",
		     x1,y1) ;
            { char *cp = stackNextText(gStk) ;
	      if (cp) while(*cp)
		{ switch (*cp)
		    { case '(': 
		      case ')': 
		      case '\\':
			fputc('\\',fil) ;
			break ;
		      }
		  fputc(*cp++,fil) ;
		}
	    }
            fprintf (fil,") show\n") ;
            break ;
          case TEXT_UP:
	    fprintf (fil,"gsave %f %f translate -90 rotate 0.6 1.0 scale"
		     " 0.2 0.0 moveto (",
		     x1, y1) ;
            { char *cp = stackNextText(gStk) ;
	      if (cp) while(*cp)
		{ switch (*cp)
		    { case '(': 
		      case ')': 
		      case '\\':
			fputc('\\',fil) ;
			break ;
		      }
		  fputc(*cp++,fil) ;
		}
	    }
            fprintf (fil,") show grestore\n") ;
            break ;
	  case TEXT_PTR:
	    fprintf (fil,"%.4f %.4f moveto (",
		     x1,y1) ;
            { char *cp = stackNext(gStk,char*) ;
	      while(*cp)
		{ switch (*cp)
		    { case '(': 
		      case ')': 
		      case '\\':
			fputc('\\',fil) ;
			break ;
		      }
		  fputc(*cp++,fil) ;
		}
	    }
            fprintf (fil,") show\n") ;
	    break ;
	  case TEXT_PTR_PTR:
	    fprintf (fil,"%.4f %.4f moveto (",
		     x1,y1) ;
            { char *cp = *stackNext(gStk,char**) ;
	      while(cp && *cp)
		{ switch (*cp)
		    { case '(': 
		      case ')': 
		      case '\\':
			fputc('\\',fil) ;
			break ;
		      }
		  fputc(*cp++,fil) ;
		}
	    }
            fprintf (fil,") show\n") ;
	    break ;
	  case TEXT_EXTERNAL:
	    fprintf (fil,"%.4f %.4f moveto (",
		     x1,y1) ;
            { char *cp = *stackNext(gStk,char**) ;
              int len=*stackNext(gStk,char*) ,i;  
	      for(i=0;cp && *cp && i<len ;i++)
		{ switch (*cp)
		    { case '(': 
		      case ')': 
		      case '\\':
			fputc('\\',fil) ;
			break ;
		      }
		  fputc(*cp++,fil) ;
		}
	    }
            fprintf (fil,") show\n") ;
	    break ;        
	  case COLOR_SQUARES:
	    { float old = colour ;
	      float oldRed = red, oldGreen = green, oldBlue = blue;
	      char *text = stackNext (gStk, char*) ;
	      int n = stackNext (gStk, int) ;
	      int s = stackNext (gStk, int) ;
	      int *tints = stackNext (gStk, int*) ;
	      int x2 = 0;
	      int new = WHITE ;
	      while (n--)
		{ switch(*text^((*text)-1))
		    {
		    case -1:   new = WHITE; break;
		    case 0x01: new = tints[0]; break;
		    case 0x03: new = tints[1]; break;
		    case 0x07: new = tints[2]; break;
		    case 0x0f: new = tints[3]; break;
		    case 0x1f: new = tints[4]; break;
		    case 0x3f: new = tints[5]; break;
		    case 0x7f: new = tints[6]; break;
		    case 0xff: new = tints[7]; break;
		    }
		  
		  if ((!isColor && Grey[new] != old) ||
		      (isColor && (Red[new] != oldRed || Green[new] != oldGreen || Blue[new] != oldBlue)))
		    { if (x2 > 0)
			{ if (isColor)
			    fprintf (fil, "%.2f %.2f %.2f setrgbcolor\n", 
				     oldRed, oldGreen, oldBlue) ;
			  else
			    fprintf (fil, "%.2f setgray\n", old) ;

			  fprintf (fil,"%.4f %.4f %.4f %.4f rectangle fill\n",
				   x1, y1, x1+x2, y1+1) ;
			  x1 += x2; x2 = 0;
			}
		      old = Grey[new];
		      oldRed = Red[new];
		      oldGreen = Green[new];
		      oldBlue = Blue[new];
		    }
		  x2++;
		  text += s ;
		}
	      if (x2 > 0)
		{ if (isColor)
		    fprintf (fil, "%.2f %.2f %.2f setrgbcolor\n", 
			     oldRed, oldGreen, oldBlue) ;
		  else
		    fprintf (fil, "%.2f setgray\n", old) ;
		  
		  fprintf (fil,"%.4f %.4f %.4f %.4f rectangle fill\n",
			   x1, y1, x1+x2, y1+1) ;
		}
	      
	      if (isColor)
		fprintf (fil, "%.2f %.2f %.2f setrgbcolor\n", red, green, blue) ;
	      else
		fprintf (fil, "%.2f setgray\n", colour) ;
	      break ;
	    }
	  case IMAGE:
	    stackSkip (gStk) ;
	    break ;
	  default:
	    messout ("Invalid action %d received in drawPSBox",action) ;
	    break ;     
	  }
      }
}

static void graphPSBox (int k, int fcol, int bcol)
{
  Box box = gBoxGet (k) ;

  currFormat = box->format ;
  currHeight = box->textheight ;
  fontSet () ;
  psize = gActive->pointsize ;
  psize2 = psize/2 ;

  if (fcol >= 0)
    box->fcol = fcol ;
  if (bcol >= 0)
    box->bcol = bcol ;

  colorInit () ;
  fillBox (box) ;	/* does background */
  drawBox (box) ;
}


static void PStitle (char *title, char *theDate)
  {
  char *firstLine = "" ;


  /* ACEDB-GRAPH INTERFACE: if acedb has registered a session name, output it*/
  if (getGraphAcedbSessionName() != NULL)
    firstLine = messprintf("%s - %s", getGraphAcedbSessionName(), theDate);

  fprintf (fil,"%%%%title page begin\n") ;
  fprintf (fil,"gsave\n") ;
  fprintf (fil,"20 756 translate 7 -12 scale\n") ;
  fprintf (fil,"0.0 setgray\n") ;
  fprintf (fil,"0 titleFont setfont\n") ;
  fprintf (fil,"0 -5 moveto\n") ;
  fprintf (fil,"( %s) show\n",firstLine) ;

  if (title && *title)
    { 
    fprintf (fil,"0 -4 moveto\n") ;
    fprintf (fil,"( %s) show\n",title) ;
    }

  fprintf (fil,"grestore\n") ;
  fprintf (fil,"%%%%title page end\n") ;

  }

/****************************************************************/

#define pageWidth 558.0
#define pageHeight 756.0

void graphPSdefaults (BOOL *pRot, float *pFac, int *pPage)
{
  float aspect = gActive->aspect ;
  float h = gActive->uh ;
  float w = gActive->uw / aspect ;
  float xfac, yfac, fac ;
  BOOL isRotate = FALSE ;
  
  switch (gActive->type)
    { 
    case TEXT_SCROLL:  
    case TEXT_FULL_SCROLL:
    case TEXT_FIT:
      if (12 * w < pageWidth) /* reduce, but never blow up */
	fac = 12.0 ;
      else if (10 * w < pageWidth)
	fac = 10.0 ;
      else 
	{ if (w > h)
	   { float x = w ; w = h ; h = x ; 
	     isRotate = TRUE ; 
	     if (12 * w < pageWidth) /* reduce, but never blow up */
	       fac = 12.0 ;
	     else if (10 * w < pageWidth)
	       fac = 10.0 ;
	     else
	       fac = .95 * pageWidth / w;
	   }
	  else
	    fac = .95 * pageWidth / w;
	}
      break ;
    case MAP_SCROLL: case PIXEL_SCROLL: case PIXEL_FIT:
    case PLAIN:
      if (w > h)
	 { float x = w ; w = h ; h = x ; isRotate = TRUE ; }
      xfac = .996 * pageWidth / w ;
      yfac = .996 * pageHeight / h ;
      fac = (xfac < yfac) ? xfac : yfac ;


      /* What is this....currently, the '0' means this is never included...  */
      /*                                                                     */
#if 0
      { char *dfault = strnew (messprintf ("%.2f", fac), 0) ;
	if (!graphPrompt (messprintf ("%s mode.  Default magnification gives "
				      "%.2f x %.2f pages.  Confirm or change "
				      "magnification (width > 1 page will "
				      "lose material off right hand side):",
				      isRotate ? "LANDSCAPE" : "PORTRAIT", 
				      w*fac/pageWidth, h*fac/pageHeight),
			dfault, "fz"))
	  { messfree (dfault) ;
	    return ;
	  }
	messfree (dfault) ;
	freefloat (&fac) ;
      }
#endif

      break ;
    default:
      messcrash ("Unknown graph type in graphPS") ;
      return ;
    }
  
  *pRot = isRotate ;
  *pFac = fac ;
  *pPage = 1 + h * fac/pageHeight ;
}

void graphPS (char *myfilname, char *mail, char *print, 
	      char *title, BOOL doColor, BOOL isRotate,
	      float fac, int pages)
{ Stack pc = 0 ;
/* char  name[128] ; */
  char *localfilname = 0 ;
  float aspect = gActive->aspect ;
  int i ;
  static BOOL firstPass = TRUE ;
  static Stack fontNameStack ;
  int fontx, fonty ;
  static char theDate[24];

  timeShowFormat (timeNow(), "%y/%m/%d %H:%M:%S", theDate, 20) ;

  if (!gFontInfo (gActive->textheight, &fontx, &fonty))
    messcrash ("Can't get font info for default font") ;
  fontHeight = fonty ;

  if(firstPass)
    { int nFonts = 0 ;
    firstPass = FALSE ;
    fontNameStack = stackCreate(50) ;
    fil = 0 ;


    /* ACEDB-GRAPH INTERFACE: Get acedb ps fonts if any registered.          */
    if (getGraphAcedbGetFonts() != NULL)
      (getGraphAcedbGetFonts())(&fil, fontNameStack, &nFonts) ;


    if (nFonts < 5)
      { stackClear (fontNameStack) ;
      pushText (fontNameStack,"Helvetica") ;           /* normal */
      pushText (fontNameStack,"Helvetica-Oblique") ;   /* italic */
      pushText (fontNameStack,"Helvetica-Bold") ;      /* bold   */
      pushText (fontNameStack,"Courier") ;             /* fixed  */
      pushText (fontNameStack,"Times-Roman") ;         /* title  */
      }
    }

  isColor = doColor ;
  if (myfilname && *myfilname)
    { 
    localfilname = myfilname ;
    if (!(fil = fopen (localfilname,"w")))
      { messout ("failed to open ps file %s", localfilname) ;
      return ;
      }
    }
  else
    { 
    fil = filtmpopen(&localfilname, "w");
    if (!fil) 
      {
      messout ("failed to open ps file %s", localfilname) ;
      return ;
      }
    }

  stackCursor(fontNameStack,0) ; /* RD adopted Jean's font fix */

  fprintf (fil,"%%!PS-Adobe-1.0    EPSF-1.2\n") ;

  if (title && *title) fprintf (fil,"%%%%Title: %s\n",title);
  else fprintf (fil,"%%%%Title: acedb generated document\n");

  fprintf (fil,"%%%%CreationDate: %s\n",theDate) ;

  /* ACEDB-GRAPH INTERFACE: Get acedb ps fonts if any registered.            */
  if (getGraphAcedbSessionName() != NULL)
    fprintf (fil,"%%%%Creator: %s\n", getGraphAcedbSessionName()) ;


  if (pages == 1)
    { float bx1, bx2, by1, by2 ;
      float ux1, ux2, uy1, uy2 ;

      graphBoxDim (0, &ux1, &uy1, &ux2, &uy2) ;
      ux1 -= gActive->ux ; ux2 -= gActive->ux ;
      uy1 -= gActive->uy ; uy2 -= gActive->uy ;
      ux1 *= fac/aspect ; ux2 *= fac/aspect ;
      uy1 *= -fac ; uy2 *= -fac ;
      if (isRotate)
	{ bx1 = uy1 ; bx2 = uy2 ;
	  by1 = -ux1 ; by2 = -ux2 ;
	  bx1 += pageWidth ; bx2 += pageWidth ;
	  by1 += pageHeight ; by2 += pageHeight ;
	}
      else
	{ bx1 = 18 + ux1 ; bx2 = 18 + ux2 ;
	  by1 = pageHeight + uy2 ;
	  by2 = pageHeight + uy1 ;
	}
      fprintf (fil,"%%%%BoundingBox: %.1f %.1f %.1f %.1f\n",
	       bx1, by1, bx2, by2) ;
    }
  fprintf (fil,"%%%%EndComments\n");

  fprintf (fil,"/rectangle    %% x1 y1 x2 y2 ==> _\n") ;
  fprintf (fil," {/y2 exch def  /x2 exch def /y1 exch def /x1 exch def\n") ;
  fprintf (fil,"  newpath\n") ;
  fprintf (fil,"    x1 y1 moveto x2 y1 lineto x2 y2 lineto x1 y2 lineto\n") ;
  fprintf (fil,"  closepath\n") ;
  fprintf (fil," } def\n\n") ;

  fprintf (fil,"/defaultFont	%% --\n") ;
  fprintf (fil,"  /%s findfont 0 -0.7 matrix translate makefont\n", 
	   stackNextText(fontNameStack)) ;
  fprintf (fil,"  %f %f matrix scale makefont def\n",
	   XtoUrel(fontHeight), XtoUrel(-fontHeight/aspect)) ;

  fprintf (fil,"/scaledFont	%% scale ==> _\n") ;
  fprintf (fil,"  { defaultFont exch scalefont } def\n\n") ;

  fprintf (fil,"/italicFont	%% --\n") ;
    fprintf (fil,"  /%s findfont 0 -0.7 matrix translate makefont\n", 
	   stackNextText(fontNameStack)) ;
  fprintf (fil,"  %f %f matrix scale makefont def\n",
	   XtoUrel(fontHeight), XtoUrel(-fontHeight/aspect)) ;

  fprintf (fil,"/boldFont	%% --\n") ;
  fprintf (fil,"  /%s findfont 0 -0.7 matrix translate makefont\n", 
	   stackNextText(fontNameStack)) ;
  fprintf (fil,"  %f %f matrix scale makefont def\n",
	   XtoUrel(fontHeight), XtoUrel(-fontHeight/aspect)) ;

  fprintf (fil,"/fixedFont	%% --\n") ;
  fprintf (fil,"  /%s findfont 0 -0.7 matrix translate makefont\n", 
	   stackNextText(fontNameStack)) ;
  fprintf (fil,"  %f %f matrix scale makefont def\n",
	   XtoUrel(fontHeight*1.026), XtoUrel(-fontHeight/aspect)) ;  
  /* changed 1.026 from height to width by esr 12/94 */

  fprintf (fil,"/titleFont      %% --\n") ;
  fprintf (fil,"  /%s findfont 0 -0.7 matrix translate makefont\n",
           stackNextText(fontNameStack)) ;
  fprintf (fil,"  %f %f matrix scale makefont def\n",
           XtoUrel(fontHeight)/UtextX, XtoUrel(-fontHeight/aspect)/UtextY) ;
/*                       mhmp /UtextX et Y 03.10.97 */

  if (title && *title) PStitle (title,theDate); 

  if (isRotate)
    fprintf (fil,"gsave\n%f %f translate -90 rotate  ",pageWidth,pageHeight) ;
  else
    fprintf (fil,"gsave\n18 %f translate  ",pageHeight) ;
  fprintf (fil,"%.4f %.4f scale  ",fac/aspect,-fac) ;
  fprintf (fil,"%.4f %.4f translate\n",-gActive->ux,-gActive->uy) ;
  fprintf (fil,"%.6f setlinewidth\n",1.0/fac) ;	       /* default 1/72 inch */

  for (i = 0 ; i < pages ; i++)
    { graphPSBox (0,-1,-1) ;

/* here we will print the page Title */
/*  if (title && *title && (i == 0))  */
/*     PStitle (title,theDate) ; mhmp 06.10.97*/

      fprintf (fil,"\ngsave showpage grestore\n\n") ;

      if (isRotate)
	fprintf (fil,"%f 0 translate\n",-aspect*pageHeight/fac) ;
      else
	fprintf (fil,"0 %f translate\n",-pageHeight/fac) ;
    }

  fprintf (fil,"grestore\n") ;
  fprintf (fil,"%%%%Trailer\n");
  fclose (fil) ;


#if !defined(MACINTOSH)
  if (print && *print &&
      (pages < 3 || messQuery (messprintf ("The document will contains %d pages.  "
					   "Should I print it ?", 
					   pages))))
    {
    pc = stackCreate(100) ;
    pushText (pc, print) ;
    catText (pc," ") ;
    catText (pc, localfilname) ;
    if (callSystem (stackText(pc, 0)))
      messout ("Sorry, for some reason print command %s failed", 
	       stackText(pc, 0)) ;
    stackDestroy (pc) ;
    }

  if (!title || !*title) title = "acedb_mail" ;

  if (mail && *mail)
    {
    if (callSystem (messprintf ("Mail -s \"%s\" %s < %s", 
			    title, mail, localfilname)))
      messout ("Sorry, for some reason mail command '%s' failed",
	       messprintf ("Mail -s \"%s\" %s < %s", title, mail, localfilname)) ;
    }
#endif /* MACINTOSH */

  if (localfilname && (!myfilname || !*myfilname))
    filremove (localfilname,"") ;
}

/******** Public function *************/
/***** end of file *****/


 
 
 
 
