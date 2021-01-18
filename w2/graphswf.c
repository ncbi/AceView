/*  File: swf.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@ncbi.nlm.nih.gov
 *
 * Description: Draws an acedb image for the web in swf format
 *              
 * Exported functions: 
 *      swfGraph2html : active graph exported as a full fledged html document
 * HISTORY:
 * Created Sept 3 2004 (mieg)
 *
 * CVS info:  
 *-------------------------------------------------------------------
 */
#include "regular.h"
#include "ac.h"
#include "graph_.h"
#include "dict.h" 
#include "aceio.h" 
#include "whooks/systags.h"

static float XFac, YFac ;

static BOOL DEBUG = FALSE ;

typedef struct swfDumpStruct { float x1, x2, y1, y2 ; DICT *colorDict ; Array backgrounds, shapes, filledShapes ; vTXT draw, bb, bb2 ; AC_HANDLE h ; } SDS ;
typedef struct swfPointStruct { int x, y ; } swfPOINT ;
typedef struct swfShapeStruct { int nn, col, bcol, lineWidth ; vTXT txt ; int oldx, oldy ; } SWSH ;

/*******************************************************************************/

static char *swfColName (int col)
{
  char *cp = name(col - WHITE + _WHITE) ;
  if (*cp == '_')
    cp = "BLACK" ;
  return cp ;
  /*
  static char buf[100] ;
  char *cp ;
  int nn ;

  strcpy (buf, name(col - WHITE + _WHITE)) ;
  cp = buf ;
  while ((nn = *cp))
    { *cp++  = ace_lower (nn) ; }
  if (!strncmp (buf, "dark", 4)) 
    memset (buf, ' ', 4) ;
  if (!strncmp (buf, "mid", 3)) 
    memset (buf, ' ', 3) ;
  return buf ;
  */
} /* swfColName */

/*******************************************************************************/
/* exports recursivelly the internal boxes
 */
static void swfDrawBox (SDS *sds, Box box, int iBox)
{
  int  i, i1 ;
  int  x1,x2,y1,y2 ;
  float zx1,zx2,zy1,zy2 ;
  int r = 1, a1, a2, dep ;
  int  fcol, bcol, xcol, format, lineWidth ;
  float textheight ;
  int    action ;
  Box	 nextbox ;
  static int nGlyph = 0 ;
  BUBBLEINFO *bubble = gActive->bubbleInfo && gActive->bubbleDict ? arrayp(gActive->bubbleInfo, iBox, BUBBLEINFO) : 0 ;  
  const char *ccp = bubble && bubble->bName ? dictName (gActive->bubbleDict, bubble->bName) : 0 ;
  SWSH *swsh ;

  x1 = x2 = y1 = y2 = 0 ;
  lineWidth = box->linewidth ;
  if (lineWidth == 0) lineWidth = 1 ;
  /* pointsize = box->pointsize ; */
  textheight = box->textheight ;
  format = box->format ;
  fcol = box->fcol ;  
  bcol = box->bcol ;
  
  stackCursor (gStk, box->mark) ;
  
  while (!stackAtEnd (gStk))
    switch (action = stackNext (gStk, int))
      {
      case BOX_END:           /* exit point */ 
	graphBoxDim (iBox, &zx1, &zy1, &zx2, &zy2) ;
	x1 = zx1* XFac ; x2 = zx2 * XFac ; y1 = zy1 * YFac ; y2 = zy2 * YFac ;
	if (bcol && bcol != TRANSPARENT) /* draw a pseudo rectangle with zero width border */
	  {
	    dictAdd (sds->colorDict,  swfColName (bcol), &i) ;
	    swsh = arrayp (sds->backgrounds, i, SWSH) ;
	    if (!swsh->txt)
	      {
		swsh->txt = vtxtHandleCreate (sds->h) ;
		swsh->bcol = bcol ;
		swsh->oldx = -1 ;
		vtxtPrintf (swsh->txt, ".outline bckg_%s:\n    "
			    , swfColName (bcol)
			    ) ;
	      }
	    if (swsh->oldx != x1 || swsh->oldy != y1)
	      {
		vtxtPrintf (swsh->txt, "M %d %d ", x1, y1) ;
		(swsh->nn)++ ;
	      }
	    /* forth only to delimit a rectangle */
	    vtxtPrintf (swsh->txt, "L %d %d ", x2, y1) ;
	    vtxtPrintf (swsh->txt, "L %d %d ", x2, y2) ;
	    vtxtPrintf (swsh->txt, "L %d %d ", x1, y2) ;
	    vtxtPrintf (swsh->txt, "L %d %d ", x1, y1) ;
	    swsh->oldx = x1 ; swsh->oldy = y1 ;
	    (swsh->nn) += 4 ;
	    if (! (swsh->nn %  25))
	      vtxtPrintf (swsh->txt, "\n    ") ;
	  }

	if (ccp)  /* bubbles method2011, 64 bits machine */
	  {
	    int bbbdx, bbbdx1, bbbdy, bbbdx2 ;
	    int bbx, bby ;
	    char *cp1 ; const char *cq1 ;
	    char buf[2000] ;
	    
	    /* find the geometry of the multine bubble and protect the \n */
	    bbbdy = 1 ;
	    bbbdx = 0 ; 
	    bbbdx1 = 0 ;
	    bbbdx2 = 1998 ; /* protect the buffer length */
	    cq1 = ccp - 1 ;
	    cp1 = buf ;

	    while (bbbdx2--)
	      switch (*++cq1)
		{
		case 0: 
		  bbbdx2 = 0 ; 
		  break ;
		case '\n': 
		  *cp1++ = '\\' ; /* double the backslash */
		  *cp1++ = 'n' ; 
		  bbbdy++ ;
		  if (bbbdx1 > bbbdx) bbbdx = bbbdx1 ;
		  bbbdx1 = 0 ;
		
		  break ;
		case ',':
		case ' ':
		  if (bbbdx1 > 50)
		    {
		      *cp1++ = '\\' ;
		      *cp1++ = 'n' ; 
		      if (bbbdx1 > bbbdx) bbbdx = bbbdx1 ;
		      bbbdx1 = 0 ; bbbdy++ ;
		      break ;
		    }
		  /* else fall thru */
		default:
		  bbbdx1++ ;
		  *cp1++ = *cq1 ;
		  break ;
		}
	    if (bbbdx1 > bbbdx) bbbdx = bbbdx1 ;	 
	    *cp1 = 0 ;
	    /* prepare a sentive area of the correct size */
	    vtxtPrintf (sds->bb, ".box  bc%d  %d %d fill=WHITE\n"
			, nGlyph, (x2 > x1 ? x2 - x1 : 1) , (y2 > y1 ? y2 - y1 : 1)
			) ;
	    /* prepare a border for the bubble */
	    bbbdx = .585 * XFac * (1 + bbbdx) ; bbbdy *=  .7 * YFac ;
	    vtxtPrintf (sds->bb, ".outline su%d:\n   M %d %d L %d %d L %d %d L %d %d L %d %d \n.end\n"
			, nGlyph
			, -10, -15
			, -10, bbbdy
			, bbbdx , bbbdy 
			, bbbdx , -15
			, -10, -15
			) ;
	    /* prepare a filled rectangle delimited by the bubble area */
	    vtxtPrintf (sds->bb, ".filled sv%d outline=su%d color=BLACK fill=BUBBLEBACKGROUND line=1\n"
			,nGlyph , nGlyph 
			) ;
	    /* prepare the text of the bubble */
	    vtxtPrintf (sds->bb, ".text  bt%d font=myfont \"%s\" color=%s\n"
			, nGlyph, buf, "BLACK"
			) ;
	    /* prepare a sprite composed of the filled rectanle and the text */
	    vtxtPrintf (sds->bb, ".sprite s%d\n  .put sv%d\n  .put bt%d scale=15%%\n.end\n"
			, nGlyph, nGlyph, nGlyph
			) ;
	    
	    /* position the bubble in an optimal way */
	    bbx =  (x1+x2)/2 - bbbdx/2 ; bby = (y1 + y2) / 2 - bbbdy / 2 ;
	    if (y2 - y1 < bbbdy + 15) /* flat box. go a bit above */
	      bby -= bbbdy/2 + 7 ;
	      
	    /* do not go out of the window */
	    if (bbx + bbbdx > XFac * sds->x2) bbx = XFac * sds->x2 - bbbdx ;
	    if (bbx <  XFac * sds->x1) bbx = XFac * sds->x1 ;
	    if (bby + bbbdy > YFac * sds->y2) bby = YFac * sds->y2 - bbbdy ;
	    if (bby <  YFac * sds->y1) bby = YFac * sds->y1 ;
	    vtxtPrintf (sds->bb, ".button bb%d\n", nGlyph) ;
	    vtxtPrintf (sds->bb, "  .show bc%d as=area x=%d y=%d\n", nGlyph, x1, y1) ;
	    vtxtPrintf (sds->bb, "  .show s%d as=hover x=%d y=%d\n", nGlyph, bbx, bby) ;
	    
	    if (bubble->fName && bubble->key)
	      {  
		vtxtPrint (sds->bb, "  .on_press:\n") ;
		if (bubble->key == 1)
		  vtxtPrintf (sds->bb, "    _root.OpenAnchor(0, \"%s\");\n"
			      , dictName (gActive->bubbleDict, bubble->fName)
			      ) ;
		else if (bubble->key)
		  vtxtPrintf (sds->bb, "    _root.OpenAceViewLink(\"%s\",\"%s\");\n"
			      , dictName (gActive->bubbleDict, bubble->fName), name(bubble->key)
			      ) ;
		vtxtPrint (sds->bb, "  .end\n") ;
	      }
	    vtxtPrint (sds->bb, ".end\n") ;
	    vtxtPrintf (sds->bb, ".put bb%d\n", nGlyph) ;
	    nGlyph++ ;
	  }
        return ;                        
      case BOX_START:          /* recursion */
        i1 = stackNext (gStk, int) ;
	if(DEBUG) fprintf(stderr,"BOX_START %d\n",i1);
	nextbox = gBoxGet (i1) ;
	swfDrawBox (sds, nextbox, i1) ;
	break ;
      case COLOR:
        i1 = stackNext (gStk,int) ;
	fcol = i1;
        break ;
      case TEXT_FORMAT:
        format = stackNext (gStk,int) ;
	/*********/
	if(DEBUG) fprintf(stderr,"TEXT_FORMAT %d\n", format);
        break ;
      case LINE_WIDTH:
	lineWidth = (int)(YFac * stackNext(gStk, float)) ;
        break ;
      case TEXT_HEIGHT:
	textheight = stackNext (gStk,float) ;
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
	if (x1 < 1) x1 = 1 ;
	if (x2 < 1) x2 = 1 ;
	if (y1 < 1) y1 = 1 ;
	if (y2 < 1) y2 = 1 ;
	if (fcol != TRANSPARENT) /* transparent fcol in acedb is used to create reactive boxes
				  * there is no need yo draw them in flash
				  */
	  switch (action)
	    {
	    case LINE:
	      dictAdd (sds->colorDict, messprintf ("%s_%d", swfColName (fcol), lineWidth), &i) ;
	      swsh = arrayp (sds->shapes, i, SWSH) ;
	      if (!swsh->txt)
		{
		  swsh->txt = vtxtHandleCreate (sds->h) ;
		  swsh->col = fcol ;
		  swsh->lineWidth = lineWidth ;
		  swsh->oldx = -1 ;
		  vtxtPrintf (swsh->txt, ".outline shape_%s_%d:\n    ", swfColName (fcol), lineWidth) ;
		}
	      if (swsh->oldx != x1 || swsh->oldy != y1)
		{
		  vtxtPrintf (swsh->txt, "M %d %d ", x1, y1) ;
		  (swsh->nn)++ ;
		}
	      vtxtPrintf (swsh->txt, "L %d %d ", x2, y2) ;
	      vtxtPrintf (swsh->txt, "L %d %d ", x1, y1) ;
	      swsh->oldx = x1 ; swsh->oldy = y1 ;
	      (swsh->nn)++ ;
	      if (! (swsh->nn %  25))
		vtxtPrintf (swsh->txt, "\n    ") ;
	      break ;
	    case RECTANGLE:
	    case FILL_RECTANGLE:
	      xcol = action == FILL_RECTANGLE ? fcol : bcol ;
	      if (xcol == TRANSPARENT) /* transparent: draw a set of lines */
		{
		  dictAdd (sds->colorDict, messprintf ("%s_%d", swfColName (fcol), lineWidth), &i) ;
		  swsh = arrayp (sds->shapes, i, SWSH) ;
		  if (!swsh->txt)
		    {
		      swsh->txt = vtxtHandleCreate (sds->h) ;
		      swsh->col = fcol ;
		      swsh->bcol = xcol ;
		      swsh->lineWidth = lineWidth ;
		      swsh->oldx = -1 ;
		      vtxtPrintf (swsh->txt, ".outline shape_%s_%d:\n    "
				  , swfColName (swsh->col),  swsh->lineWidth
				  ) ;
		    }
		  if (swsh->oldx != x1 || swsh->oldy != y1)
		    {
		      vtxtPrintf (swsh->txt, "M %d %d ", x1, y1) ;
		      (swsh->nn)++ ;
		    }
		  /* back and forth to obtain a transparent result */
		  vtxtPrintf (swsh->txt, "L %d %d ", x2, y1) ;
		  vtxtPrintf (swsh->txt, "L %d %d ", x2, y2) ;
		  vtxtPrintf (swsh->txt, "L %d %d ", x1, y2) ;
		  vtxtPrintf (swsh->txt, "L %d %d ", x1, y1) ;
		  vtxtPrintf (swsh->txt, "L %d %d ", x1, y2) ;
		  vtxtPrintf (swsh->txt, "L %d %d ", x2, y2) ;
		  vtxtPrintf (swsh->txt, "L %d %d ", x2, y1) ;
		  vtxtPrintf (swsh->txt, "L %d %d ", x1, y1) ;
		  swsh->oldx = x1 ; swsh->oldy = y1 ;
		  (swsh->nn) += 8 ;
		}
	      else
		{
		  dictAdd (sds->colorDict, messprintf ("%s_%s_%d",  swfColName (fcol), swfColName (xcol), lineWidth), &i) ;
		  swsh = arrayp (sds->filledShapes, i, SWSH) ;
		  if (!swsh->txt)
		    {
		      swsh->txt = vtxtHandleCreate (sds->h) ;
		      swsh->col = fcol ;
		      swsh->bcol = xcol ;
		      swsh->lineWidth = lineWidth ;
		      swsh->oldx = -1 ;
		      vtxtPrintf (swsh->txt, ".outline filledshape_%s_%s_%d:\n    "
				  , swfColName (swsh->col), swfColName (xcol), swsh->lineWidth
				  ) ;
		    }
		  if (swsh->oldx != x1 || swsh->oldy != y1)
		    {
		      vtxtPrintf (swsh->txt, "M %d %d ", x1, y1) ;
		      (swsh->nn)++ ;
		    }
		  vtxtPrintf (swsh->txt, "L %d %d ", x2, y1) ;
		  vtxtPrintf (swsh->txt, "L %d %d ", x2, y2) ;
		  vtxtPrintf (swsh->txt, "L %d %d ", x1, y2) ;
		  vtxtPrintf (swsh->txt, "L %d %d ", x1, y1) ;
		  swsh->oldx = x1 ; swsh->oldy = y1 ;
		  (swsh->nn) += 4 ;
		}
	      if (! (swsh->nn %  25))
		vtxtPrintf (swsh->txt, "\n    ") ;
	      break ;
	    }
	break ;
      case PIXELS: case PIXELS_RAW:
	{
	  int w ;
	  
	  x1 = (int)stackNext (gStk,float) ;
	  y1 = (int)stackNext (gStk,float) ;
	  stackSkip (gStk) ; /* consume the pixels */
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
	  int n = stackNext (gStk, int), x0, y0 ;

	  if (n > 2) 
	    {
	      int i;
	      	
	      dictAdd (sds->colorDict, messprintf ("%s_%d", swfColName (fcol), lineWidth), &i) ;
	      swsh = arrayp (sds->shapes, i, SWSH) ;
	      if (!swsh->txt)
		{
		  swsh->txt = vtxtHandleCreate (sds->h) ;
		  swsh->col = fcol ;
		  swsh->lineWidth = lineWidth ;
		  swsh->oldx = -1 ;
		  vtxtPrintf (swsh->txt, ".outline shape_%s_%d:\n    ", swfColName (fcol), lineWidth) ;
		}
	      x0 = (int)(XFac * stackNext (gStk, float));
	      y0 = (int)(YFac * stackNext (gStk, float)); 
	      if (swsh->oldx != x0 || swsh->oldy != y0)
		{ 
		  (swsh->nn)++ ;
		  vtxtPrintf (swsh->txt, "M %d %d ", x0, y0) ;
		}
	      for(i=0;i<n-1;i++)
		{
		  x1 = (int)(XFac * stackNext (gStk, float));
		  y1 = (int)(YFac * stackNext (gStk, float)); 
		  vtxtPrintf (swsh->txt, "L %d %d ", x1, y1) ;
		  (swsh->nn)++ ;
		  if (! (swsh->nn %  25))
		    vtxtPrintf (swsh->txt, "\n    ") ;
		}
	      swsh->oldx = x1 ; swsh->oldy = y1 ;
	      if (i && action == POLYGON) /* close up */
		{
		  vtxtPrintf (swsh->txt, "L %d %d ", x0, y0) ;
		  swsh->oldx = x0 ; swsh->oldy = y0 ;
		}
	    }
	  break ;
	}
      case CIRCLE: case POINT: 
      case TEXT: case TEXT_UP: case TEXT_PTR: case TEXT_PTR_PTR: case TEXT_EXTERNAL:
      case COLOR_SQUARES: case FILL_ARC: case ARC:
	r = 1 ;	  
	x1 = (int)(XFac * stackNext (gStk,float));
	y1 = (int)(YFac * stackNext (gStk,float));
	switch (action)
	  {
	  case CIRCLE:
	    r = (int)(XFac * stackNext(gStk,float));
	    /* fall on POINT */
	  case POINT:
	    vtxtPrintf (sds->draw, ".circle b%d r=%d color=%s line=1 fill=%s\n"
		     , nGlyph, r
		     , swfColName (fcol)
		     , swfColName (box->bcol)
		     ) ;
	    vtxtPrintf (sds->draw, ".put b%d  x=%d y=%d\n"
		     , nGlyph, x1 - r, y1
		     ) ;
	    nGlyph++ ;
	    break ;
	  case FILL_ARC: case ARC:
	    r = (int)(2.0 * XFac * stackNext(gStk,float));
	    a1 = (int)stackNext(gStk,float);
	    a2 = (int)stackNext(gStk,float);
	    dep = a1 + a2;
	    if(dep < 0) dep += 360;
	    /*
	      int arr = a1;
	      swfImageArc(sds, im,x1,y1,r,r,dep,arr,alloca_fcols[fcol]);
	    */
	    break ;
	  case TEXT:
	  case TEXT_PTR:
	  case TEXT_PTR_PTR: 
	  case TEXT_EXTERNAL: 
	  case TEXT_UP:
	    { 
	      char *text = 0 ;
	      
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
		  stackSkip (gStk) ; /* consume len */
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
		      case '&': break ;   /* ?????? kills swf */
		      case '\"':
			/* case '\'': */
		      case '\\':
			*cq++ = '\\' ; *cq++ = *cp ; break ;
		      default: *cq++ = *cp ; break ; 
		      }
		  *cq = 0 ;
		  vtxtPrintf (sds->draw, ".text t%d font=myfont \"%s\" color=%s size=%d%% \n"
			      , nGlyph, text2
			      , swfColName (fcol)
			      , textheight ? (int)(16 * textheight) : 16
			   ) ;
		  vtxtPrintf (sds->draw, ".put t%d x=%d y=%d\n", nGlyph++, x1, (int)(y1 + .8*YFac)) ;
		  
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
		  /* swfImageFilledRectangle (sds, im,x1,y1,x2,y2,alloca_fcols[col]); */
		  text += iskip ;
		  x1 += XFac ; x2 += XFac ; if (col) x2 += 0 ; /* for compiler happiness */
		}
	    }
	    break ;
	  default:
	    if(DEBUG) fprintf(stderr,"Invalid action %d received in drawPSBox",action) ;
	    messout ("Invalid action %d received in swfDrawBox",action) ;
	    break ;     
	  }
      }
} /* swfDrawBox */

/*******************************************************************************/

static void swfDoBox (SDS *sds, int iBox, int fcol, int bcol)
{
  Box box = gBoxGet (iBox) ;

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

  swfDrawBox (sds, box, iBox) ;
} /* swfDoBox (int box) */

/*******************************************************************************/
/* exports the outer box as a single java script
 */
static void swfOuterBox (ACEOUT out, int box, BOOL dumpBubble)
{
  float x1, y1, x2, y2, w, h ; 
  int ii, iw, ih ;
  char *filNam = "box.swf", *fontFile ;
  SDS sds ;
  SWSH *swsh ;

  sds.h = ac_new_handle () ;
  sds.colorDict = dictHandleCreate (128, sds.h) ; 
  sds.draw = vtxtHandleCreate (sds.h) ;
  sds.bb = vtxtHandleCreate (sds.h) ;
  sds.bb2 = vtxtHandleCreate (sds.h) ; /* bubble active zone */
  sds.backgrounds = arrayHandleCreate (128, SWSH, sds.h) ; /* one per bcol */
  sds.shapes = arrayHandleCreate (128, SWSH, sds.h) ; /* one per col/bcol */
  sds.filledShapes = arrayHandleCreate (128, SWSH, sds.h) ; /* one per col/bcol/linewidth */
  XFac = gActive->xFac ;
  YFac = gActive->yFac ;

  graphBoxDim (0, &x1, &y1, &x2, &y2) ;
  w = XFac * (x2 - x1) ; h = YFac * (y2 - y1) ;
  if (gActive->h < h) h = gActive-> h ;
  if (gActive->w < w) w = gActive-> w ;
  sds.x1 = x1 ; sds.x2 = x2 ;
  sds.y1 = y1 ; sds.y2 = y2 ;
  
  iw = w ; ih = h ;
  aceOutf (out, ".flash bbox=%dx%d filename=\"%s\" compress\n", iw, ih, filNam) ;

  if ((fontFile = filName("font/Tahoma","ttf","r")))
    aceOutf (out, ".font myfont filename=%s\n", fontFile) ;
  else 
    messerror ("Cannot find the font file font/Tahoma.ttf in graphswf.c") ;
  aceOutf (out, ".box fond width=%d height=%d color=gray fill=white line=2\n", iw, ih) ; 
  aceOut (out, ".put fond\n") ;
  swfDoBox (&sds, 0, -1, -1) ;  

  /* export the different sorts of glyps */
  /* backgrounds */
  for (ii = 0 ; ii < arrayMax (sds.backgrounds) ; ii++)
    {
      swsh = arrp (sds.backgrounds, ii, SWSH) ;
      if (swsh->txt)
	{
	  aceOut (out, vtxtPtr (swsh->txt)) ;
	  aceOut (out, "\n.end\n") ;
	  aceOutf (out, ".filled bbckg_%s  outline=bckg_%s color=%s fill=%s line=0\n"
		   , swfColName (swsh->bcol)
		   , swfColName (swsh->bcol)
		   , swfColName (swsh->bcol)
		   , swfColName (swsh->bcol)
		   ) ;
	  aceOutf (out, ".put  bbckg_%s\n"
		   , swfColName (swsh->bcol)
		   ) ;
	}
    }
  /* filledShapes */
  for (ii = 0 ; ii < arrayMax (sds.filledShapes) ; ii++)
    {
      swsh = arrp (sds.filledShapes, ii, SWSH) ;
      if (swsh->txt)
	{
	  aceOut (out, vtxtPtr (swsh->txt)) ;
	  aceOut (out, "\n.end\n") ;
	  aceOutf (out, ".filled bfilledshape_%s_%s_%d outline=filledshape_%s_%s_%d color=%s fill=%s line=%d\n"
		   , swfColName (swsh->col), swfColName (swsh->bcol), swsh->lineWidth
		   , swfColName (swsh->col), swfColName (swsh->bcol), swsh->lineWidth
		   , swfColName (swsh->col), swfColName (swsh->bcol)
		   , swsh->lineWidth
		   ) ;
	  aceOutf (out, ".put bfilledshape_%s_%s_%d\n"
		   , swfColName (swsh->col), swfColName (swsh->bcol), swsh->lineWidth
		   ) ;
	}
    }
  /* texts */
  aceOut (out, vtxtPtr (sds.draw)) ;
  /* lines */
  for (ii = 0 ; ii < arrayMax (sds.shapes) ; ii++)
    {
      swsh = arrp (sds.shapes, ii, SWSH) ;
      if (swsh->txt)
	{
	  aceOut (out, vtxtPtr (swsh->txt)) ;
	  aceOut (out, "\n.end\n") ;
	  aceOutf (out, ".filled bshape_%s_%d outline=shape_%s_%d color=%s fill=none line=%d\n"
		   , swfColName (swsh->col), swsh->lineWidth
		   , swfColName (swsh->col), swsh->lineWidth
		   , swfColName (swsh->col)
		   , swsh->lineWidth
		   ) ;
	  aceOutf (out, ".put bshape_%s_%d\n"
		   , swfColName (swsh->col), swsh->lineWidth
		   ) ;
	}
    }
  /* bubbles: exported as flash buttons */
  if (dumpBubble)
    {
      if (1) /* old method */
	{
	  if (vtxtPtr (sds.bb))
	    aceOut (out, vtxtPtr (sds.bb)) ;
	  if (DEBUG) fprintf(stderr, "%s\n", vtxtPtr (sds.bb)) ;
	}
      else
	{
	  if (vtxtPtr (sds.bb2))
	    {
	      aceOut (out, ".outline bubbleZone:\n") ;
	      aceOut (out, vtxtPtr (sds.bb2)) ;
	      aceOut (out, ".end\n") ;
	      aceOut (out, ".button bubbleButton\n") ;
	      aceOut (out, "  .show bubbleZone  as=area\n") ;
	      aceOut (out, "  .show bbb as=hover x=200 y=200\n") ;
	      aceOut (out, "  .on_move_in:\n") ;
	      aceOut (out, "     _root.bubble_text=\"Danou est la plus belle\" ;\n") ;
	      aceOut (out, "  .end\n") ;
	      aceOut (out, "  .on_press:\n") ;
	      aceOut (out, "    _root.OpenAceViewLink(\"hello\",\"world\") ;;\n") ;
	      aceOut (out, "  .end\n") ;
	      aceOut (out, ".end\n") ;
	      aceOut (out, ".put bubbleButton\n") ;
	    }
	}
    }
  
  ac_free (sds.h) ;

  return ;
} /* swfOuterBox */

/*******************************************************************************/
/* resize a flash window before exportation */
BOOL swfGraphResize (float width, float height)
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
 * as a .swf document that will be compiled into a flash document
 * return the number of bytes written
 */
BOOL swfGraphExport (Graph gId, ACEOUT out, BOOL dumpBubble)
{
  swfOuterBox (out, 0, dumpBubble) ;
 
  if (1)
    {
      aceOut (out, ".action:\n") ;
      aceOut (out, "  function OpenAceViewLink (oclass, onam) {\n") ;

      aceOut (out, "  this.geturl(\'javascript:openAceViewLink (\\\'\'+oclass+\'\\\',\\\'\'+ onam+\'\\\')\' );\n");

      aceOut (out, "  }\n") ;
      aceOut (out, "  function OpenAnchor (oclass, onam) {\n") ;

      aceOut (out, "  this.geturl(\'javascript:openAnchor (\\\'\'+oclass+\'\\\',\\\'\'+ onam+\'\\\')\' );\n");

      aceOut (out, "  }\n") ;
      aceOut (out, ".end\n") ;
    }
  aceOut (out, ".end\n") ;
  aceOut (out, ".exit\n") ;

  return TRUE ;
} /* swfActiveWindow */
 
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

#ifdef JUNK

example of a full document

.flash filename="box.swf"
    .box b1 100 100 color=yellow fill=red
    .put b1 pin=center scale=0%
    .frame 100
    .change b1 pin=center scale=100%
    .frame 200
    .change b1 pin=center scale=0%
.end

#endif

