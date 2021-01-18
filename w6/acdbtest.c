/*  Last edited: Aug 27 11:52 1996 (srk) */

/* @(#)acdbtest.c	1.4    8/27/96 */
#include "acedb.h"
#include "array.h"
#include "lex.h"
#include "classes.h"
#include "query.h" 
#include "graph.h"
extern void verifyModels(void);

void acedbtest(void)
{
  FILE *fil = 0 ;
  Graph g = graphCreate (PIXEL_SCROLL, "GifTest",.2, .2, .4, .6) ;
  float x0=0, y0=0 ;
  static char dir[DIR_BUFFER_SIZE], nom[FIL_BUFFER_SIZE] ;

  if (!*dir) strcpy (dir,"/users/www/ace_images/elegans") ;
  if (g && (fil = filqueryopen(dir, nom, "gif", "r","Choose a gif file")))
    {
      graphClear() ;  
      graphPixelBounds(300,300) ;
      if (graphGIFread (fil, x0, y0))
	{
	  graphBoxStart () ;
	  graphColor (RED) ;
	  graphText ("TOTO", 12, 12) ;
	  graphBoxEnd () ;
	}
      else
	{
	  graphClear() ;
	  graphBoxStart () ;
	  graphColor (RED) ;
	  graphText ("GifReadFailed", 12, 12) ;
	  graphBoxEnd () ;
	}
      graphRedraw() ;
 	
    }
  else
    if (g) graphDestroy () ;
  
  return ;
  


  if (messQuery ("Do you want to compute the clones genetic position by interpolation ?"))
    {
      /* gMapPhysNameClones () ; */
      messout("See $ACEDB/clone_full_name.ace") ;
    }
  /*
    verifyModels();// fails on #types 
    messout("Model check complete.");
    */
}

#ifdef JUNK

The Image program code is available under acedb cvs in directories
im2src, im3src, imhelp.  Here is a relevant bit of code.  Note this
is for monochrome images.


/* this function is called after each zoom operation
 * it will redisplay the image in the image window
 * This function is generic for step #1 and #2
 * To re-zoom the image, just pass im_zoomed as NULL,
 * and a new zoomed bitmap will be created, this is useful
 * if the zoomX/Y has changed since last time */
void redrawImageGraph (Graph graph,
                       float old_zoomX, float old_zoomY,
                       float new_zoomX, float new_zoomY,
                       int x_pos,
                       IMAGEINFO *im, BOOL reZoom)
{
  int nx = 0, ny = 0;           /* dimension of the zoomed image */
  
  graphActivate (graph);

  graphClear();

  nx = (int) (im->width * new_zoomX);
  ny = (int) (im->height * new_zoomY);
 
  if (reZoom)
    /* need to create a new zoomed image */
    {
      graphBusyCursorAll();     /* it may take a while to zoom */

      messfree (im->zoomed);

      /* if we zoom the image bigger than the orignal, we use
       * smoothResize */
      if (nx > im->width || ny > im->height)
        im->zoomed = smoothResize (im->rawdata,
                                   im->width, im->height,
                                   CORRX (nx), ny,
                                   0);
      else
        /* if it is zoomed to be smaller than the original
         * a blockResize will do just fine */
        im->zoomed = BlockResize (im->rawdata,
                                  im->width, im->height,
                                  CORRX (nx), ny,
                                  0);
    }
  
  graphPixelsRaw (im->zoomed,
                  (float) nx, (float) ny, (float) (CORRX (nx)),
                  x_pos, 0);
  
  
  graphPixelBounds (nx + x_pos, ny);
  
  graphLinewidth (1);
  graphColor (BLACK);
  graphRectangle (0, 0, nx + x_pos, ny);
  
  /* we need to call graphRedraw(), once all the elements of the 
   * display have been drawn, but we might want to draw stuff over
   * the top of the image (e.g. lanemap), so we don't call it here yet */

  return;
} /* redrawImageGraph */

#endif
