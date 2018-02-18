/*  File: graphgif.c
 *  Author:  Richard Durbin (rd@sanger.cam.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 * -------------------------------------------------------------------
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Draws a gif image using the drawing stack for the
 *              active window.
 *              
 * Code for ACEDB graph package by Guy Decoux, decoux@moulon.inra.fr
 * based on the gd package by Tom Boutell, boutell@netcom.com
 * adapted by Richard Durbin and Jaime Prilusky, ACEDB'95
 *              
 * Exported functions: graphGIF(), graphGIFname()
 * HISTORY:
 * Last edited: May  2 14:13 2003 (edgrif)
 * * May 10 09:57 1999 (edgrif): Fix very old bug in TEXT_PTR_PTR drawing.
 * * Feb  4 20:12 1999 (edgrif): Insert Lincolns fix for no arcs being drawn
 *              for gif only.
 * * Jan 21 14:31 1999 (edgrif): graphGIFfile now static. New graphGIF
 *              function to mirror graphPS. (old graphGIF function withdrawn,
 *              not used anywhere)
 * Created: Thu Jan 21 14:29:33 1999 (edgrif)
 * CVS info:   $Id: graphgd.c,v 1.8 2015/01/25 04:55:27 mieg Exp $
 *-------------------------------------------------------------------
 */

#include <wh/regular.h>
#include <wh/aceio.h>
#include <wh/graph_.h>
#include <wgd/gd.h>
#include <w2/graphcolour.h>
#include <wh/key.h>

#ifndef ESUCCESS
#define ESUCCESS 0
#endif /* !ESUCCESS */

/* Globals...                                                                */
extern gdFontPtr gdFont6x9, gdFont8x13, gdFont8x13bold ;

typedef struct {
  int alloca_colours[NUM_TRUECOLORS];
  gdImagePtr im;
  Graph_ graph;
} gifContext;


static int color2gdColor(gifContext *con, int col)
{
  int r, g, b;
  
  if (col >= NUM_TRUECOLORS)
    col = 0;

  if (con->alloca_colours[col])
    return con->alloca_colours[col];
 
  graphGetColourAsInt(col, &r, &g, &b, NULL);
  con->alloca_colours[col] = gdImageColorAllocate(con->im, r, g, b);
  
  return con->alloca_colours[col];
}

static gdFontPtr GIFfont(gifContext *cont, float textHeight, int textFormat)
{
  int height = cont->graph->yFac * textHeight ;

  if ((textHeight != 0) && (height < 11))
    return gdFont6x9 ;
  else if (textFormat == BOLD)
    return gdFont8x13bold ;
  else
    return gdFont8x13 ;
}

/***********************************************/

static void GIFdrawBox (gifContext *gif, Box box)
{
  int  x1,x2,y1,y2 ;
  float lineWidth = box->linewidth;
  /*  GraphLineStyle lineStyle = box->line_style ; */
  float pointSize = box->pointsize;
  float textHeight = box->textheight;
  int textFormat = box->format;
  int colour = box->fcol;
  int    action ;
  
  if ((box->bcol != TRANSPARENT) && (box->x1 <=  999999.0)) 
    gdImageFilledRectangle(gif->im,
			   (int)(gif->graph->xFac * box->x1),
			   (int)(gif->graph->yFac * box->y1),
			   (int)(gif->graph->xFac * box->x2),
			   (int)(gif->graph->yFac * box->y2),
			   color2gdColor(gif, box->bcol));
  
  stackCursor (gif->graph->stack,box->mark) ;

  while (!stackAtEnd (gif->graph->stack))
    switch (action = stackNext (gif->graph->stack, int))
      {
      case BOX_END:
        return ;                        /* exit point */
      case BOX_START:
        GIFdrawBox (gif, arrp (gif->graph->boxes, 
			       stackNext(gif->graph->stack, int),
			       BoxStruct));
	break ;
      case COLOR:
        colour = stackNext (gif->graph->stack,int) ;
	break ;
      case TEXT_FORMAT:
        textFormat = stackNext (gif->graph->stack,int) ;
	break;
      case LINE_WIDTH:
	lineWidth = stackNext(gif->graph->stack, float);
	break ;
	/*
      case LINE_STYLE:
	lineStyle = stackNext(gif->graph->stack, int);
	break ;
	*/
      case TEXT_HEIGHT:
	textHeight = stackNext (gif->graph->stack,float) ;
	break ;
      case POINT_SIZE:
        pointSize = stackNext (gif->graph->stack,float) ;
        break ;
      case LINE: case RECTANGLE: case FILL_RECTANGLE:
        x1 = (int)(gif->graph->xFac * stackNext (gif->graph->stack,float));
        y1 = (int)(gif->graph->yFac * stackNext (gif->graph->stack,float));
        x2 = (int)(gif->graph->xFac * stackNext (gif->graph->stack,float));
        y2 = (int)(gif->graph->yFac * stackNext (gif->graph->stack,float));
	switch (action)
	  {
	  case LINE:
	    /* This is horrible, we need to upgrade to a level of libgd which
	     * will draw lines with specified width/dashes etc. etc. */
	    /*
	      if (lineStyle == GRAPH_LINE_DASHED)
	      gdImageDashedLine(gif->im,x1,y1,x2,y2,color2gdColor(gif,colour));
	    */
	    if(x1 == x2) 
	      {
		int ext = (int)(gif->graph->xFac * lineWidth);
		gdImageFilledRectangle(gif->im,x1,y1,x2+ext,y2,
				       color2gdColor(gif, colour));
	      } 
	    else if(y1 == y2) 
	      {
		int ext = (int)(gif->graph->yFac * lineWidth);
		gdImageFilledRectangle(gif->im,x1,y1,x2,y2+ext,
				       color2gdColor(gif, colour));
	      } 
	    else
	      gdImageLine(gif->im,x1,y1,x2,y2,color2gdColor(gif,colour));
	    break ;
	  case RECTANGLE:
	    gdImageRectangle(gif->im,x1,y1,x2,y2,color2gdColor(gif, colour));
	    break ;
	  case FILL_RECTANGLE:
	    gdImageFilledRectangle(gif->im,x1,y1,x2,y2,
				   color2gdColor(gif,colour));
	    break ;
	  }
	break ;
      case PIXELS: 
	/* Not implemented */
	(void)stackNext (gif->graph->stack,float) ;
	(void)stackNext (gif->graph->stack,float) ;
	/* Fall through */
      case PIXELS_RAW:
	(void)stackNext (gif->graph->stack,float) ;
	(void)stackNext (gif->graph->stack,float) ;
	break ;
      case POLYGON: case LINE_SEGS:
	{ 
	  int n = stackNext (gif->graph->stack, int) ;
	  if (n > 2) 
	    {	      
	      gdPointPtr points, p;
	      int i;
	      
	      points = (gdPointPtr)messalloc((n - 1) * sizeof(gdPoint));
	      for(i=0,p=points;i<n-1;i++,p++) 
		{
		  p->x = (int)(gif->graph->xFac * 
			       stackNext (gif->graph->stack, float));
		  p->y = (int)(gif->graph->yFac * 
			       stackNext (gif->graph->stack, float));
		}
	      if (action == POLYGON)
		gdImageFilledPolygon(gif->im, points, n-1,
				     color2gdColor(gif, colour));
	      else
		gdImagePolygon(gif->im, points, n-1, 
			       color2gdColor(gif,colour));
	      messfree(points);
	    }
	  break ;
	}
      case CIRCLE: case POINT: 
      case TEXT: case TEXT_UP: case TEXT_PTR: case TEXT_PTR_PTR: 
      case COLOR_SQUARES: case FILL_ARC: case ARC:
	x1 = (int)(gif->graph->xFac * stackNext (gif->graph->stack,float));
	y1 = (int)(gif->graph->yFac * stackNext (gif->graph->stack,float));
	switch (action) 
	  {
	  case CIRCLE:
	    {
	      int r = (int)(gif->graph->xFac * 
			    stackNext(gif->graph->stack,float));
	      gdImageArc(gif->im,x1,y1,2*r,2*r,0,360,
			 color2gdColor(gif, colour));
	      break ;
	    }
	  case FILL_ARC: case ARC:
	    {
	      int r = (int)(2.0 * gif->graph->xFac * 
			    stackNext(gif->graph->stack,float));
	      int a1 = (int)stackNext(gif->graph->stack,float);
	      int a2 = (int)stackNext(gif->graph->stack,float);
	      int dep = a1 + a2;
	      int arr = a1;
	      if(dep < 0) dep += 360;
	      gdImageArc(gif->im,x1,y1,r,r,dep,arr,
			 color2gdColor(gif,colour));
	      break ;
	    }
	  case POINT:
	    {
	      int pszx = (int)(gif->graph->xFac * pointSize / 2);
	      int pszy = (int)(gif->graph->yFac * pointSize / 2);
	      gdImageFilledRectangle(gif->im,x1-pszx,y1-pszy,x1+pszx,y1+pszy,
				     color2gdColor(gif, colour)) ;
	      
	      break ;
	    }
	  case TEXT:
	    { 
	      char *text = stackNextText(gif->graph->stack);
	      gdFontPtr font = GIFfont(gif, textHeight, textFormat);
	      gdImageString(gif->im, font, x1, y1,text,
			    color2gdColor(gif, colour));
	      break;
	    }
	  case TEXT_UP:
	    { 
	      char *text = stackNextText(gif->graph->stack);
	      int n = strlen(text) ;
	      gdFontPtr font = GIFfont(gif, textHeight * 0.6 , textFormat);
	      int textHeightPixels = gif->graph->yFac * textHeight * 0.6;
	      char buf[2] ;
	      buf[1] = 0 ;
	      while (n--)
		{ 
		  buf[0] = *text++ ;
		  gdImageString(gif->im,font,x1+1,y1-textHeightPixels*(n+1),buf,
				color2gdColor(gif, colour));
		}
	      break;
	    }
	  case TEXT_PTR:
	    { 
	      char *text = stackNext(gif->graph->stack, char *);
	      gdFontPtr font = GIFfont(gif, textHeight, textFormat);
	      gdImageString(gif->im,font,x1,y1,text,
			    color2gdColor(gif, colour));
	      break;
	    }
	  case TEXT_PTR_PTR:
	    {
	      char *text = *stackNext(gif->graph->stack, char **);
	      gdFontPtr font = GIFfont(gif, textHeight, textFormat);
	      if (text && *text)
		gdImageString(gif->im,font,x1,y1,text,
			      color2gdColor(gif, colour));
	      break;
	    }
	  case COLOR_SQUARES:
	    {
	      char *text = stackNext (gif->graph->stack, char*) ;
	      int n = stackNext (gif->graph->stack, int) ;
	      int iskip = stackNext (gif->graph->stack, int) ;
	      int *tints = stackNext (gif->graph->stack, int*) ;
	      int col = 0 ;
	      x2 = x1 + gif->graph->xFac ; y2 = y1 + gif->graph->yFac ;
	      while (n--) 
		{ 
		  switch(*text^((*text)-1))
		    {
		    case -1:   col = WHITE; break;
		    case 0x01: col = tints[0]; break;
		    case 0x03: col = tints[1]; break;
		    case 0x07: col = tints[2]; break;
		    case 0x0f: col = tints[3]; break;
		    case 0x1f: col = tints[4]; break;
		    case 0x3f: col = tints[5]; break;
		    case 0x7f: col = tints[6]; break;
		    case 0xff: col = tints[7]; break;
		    }
		  gdImageFilledRectangle (gif->im,x1,y1,x2,y2,
					  color2gdColor(gif, col));
		  text += iskip ;
		  x1 += gif->graph->xFac ; x2 += gif->graph->xFac ;
		}
	      break ;
	    }
	  default:
	    messout ("Invalid action %d received in drawPSBox",action) ;
	    break ;     
	  }
      }
}

/****************** public routines ******************/

int graphGIF (Graph gId, ACEOUT out, BOOL do_size)
{
  int i, size, err;
  void* data;
  gifContext gif;

  gif.graph = gActive ; /* gGetStruct(gId); */
  if (!gif.graph)
    return EINVAL;
  
  for (i = 0 ; i < NUM_TRUECOLORS ; ++i) /* clear colour table!! */
    gif.alloca_colours[i] = 0 ;
  
  if (!(gif.im = gdImageCreate (gif.graph->w, gif.graph->h)))
    return FALSE ;
  
  GIFdrawBox (&gif, arrp (gif.graph->boxes, 0, BoxStruct)) ;
  
  gdImageInterlace (gif.im, 1) ;

  /* making the white transparent speeds up loading over the web LS*/
  gdImageColorTransparent(gif.im,0);

  data = gdImageGifPtr (gif.im, &size) ;
  gdImageDestroy (gif.im) ;

  if (!data)
    return FALSE;
    
  if (do_size)
    err = aceOutf (out, "// %d bytes\n", size);
  else
    err = ESUCCESS;
  
  if (err == ESUCCESS)
    aceOutBinary (out, data, size);
      
  free(data);

  return err;
}

/************************************************/

int graphGIFread (FILE *fil, float x0, float y0)
{ 
  unsigned int *colors;
  char *pixels, *p ;
  int i, j, ncol, x, y ;
  gdImagePtr im;

  if (!(im = gdImageCreateFromGif (fil)))
    return FALSE ;
  
  printf ("Got .gif image: %d x %d with %d colors\n",
	  gdImageSX(im), gdImageSY(im), gdImageColorsTotal(im)) ;

  ncol = gdImageColorsTotal(im) ;

  colors = halloc(ncol*sizeof(unsigned int), graphClearHandle());
  for (i = 0 ; i < ncol ; ++i)
    colors[i] = 
      gdImageRed(im, i)<<16 |
      gdImageGreen(im, i)<<8 |
      gdImageBlue(im, i) ;

  x = gdImageSX(im) ; while (x % 8) ++x ;
  y = gdImageSY(im) ;

  pixels = (char*) halloc(x*y*sizeof(unsigned char), graphClearHandle());

  for (j = 0 ; j < y ; ++j)
    for (p = pixels + j*x, i = 0 ; i < gdImageSX(im) ; ++i)
      *p++ = (char)im->pixels[i][j] ;
  
  i = graphBoxStart () ;
       
  /*
    new ace9 code
    void	graphPixels (unsigned char *pixels, int w, int h, 
    int lineWidth,
    float x0, float y0, 
    unsigned int *colors, int ncols) ;

    graphPixels(pixels, gdImageSX(im), y, x, x0, y0, colors, ncol) ;
  */
  /* current code
     GRAPH_FUNC_DCL void	graphPixels (char *pixels, int w, int h, int lineWidth,
     float x0, float y0, float x1, float y1) ;
  */
  graphPixelsRaw (pixels, gdImageSX(im), y, x, x0, y0) ;
 
  graphBoxEnd () ;
  graphBoxDraw (i, -1, -1) ;

  gdImageDestroy (im) ;

  return i ;
}


void graphGIFLeftDown (int x, int y)
 {
   gLeftDown(x/gActive->xFac, y/gActive->yFac /* , NULL_MODKEY */);

   return ;
 }




/***** end of file *****/
