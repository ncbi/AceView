/*  File: fmapstatus.c
 *  Author: Mike Holman
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: sequence display and manipulation for fmap
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 16 14:42 1998 (fw)
 * * Jul 16 10:10 1998 (edgrif): Introduce private header fmap_.h
 * Created: Sat Jul 25 22:36:47 1992 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: fmapstatus.c,v 1.1.1.1 2002/07/19 20:22:59 sienkiew Exp $ */

#include "fmap_.h"
#include "dna.h"
#include "systags.h"
#include "display.h"
#include "bump.h"

/************************************************/ 

static Graph colorkeyGraph = 0 ;
static BOOL colorKeyWanted = FALSE ;

/* cosmid_color  - used to generate color related to each cosmid status. */

#define CSM (3.8)         /* (C)osmid (S)tatus (M)ax right side for drawing. */

static enum Colour cosmid_color (KEY val, float * x, char * buf)
{
  enum Colour tcolor = BLACK ;
  KEY _Shotgun_complete;
  KEY _Library_construction;
  KEY _Received;
  KEY _Dna_made;
    
    lexaddkey("Shotgun_complete", &_Shotgun_complete, 0) ;
    lexaddkey("Library_construction", &_Library_construction, 0) ;
    lexaddkey("Received", &_Received, 0) ;
    lexaddkey("Dna_made", &_Dna_made, 0) ;
    
    if (val == _Archived)
      { if (*x < CSM) {*x = CSM; tcolor = BLACK; strcpy(buf,"a");} }
    else if (val == _Submitted)
      { if (*x < CSM-0.3) {*x = CSM-0.3;tcolor = DARKRED; strcpy(buf,"sb");} }
    else if (val == _Finished)
      { if (*x < CSM-0.6) {*x = CSM-0.6;tcolor = RED; strcpy(buf,"f");} }
    else if (val == _Contiguous)
      { if (*x < CSM-0.9) {*x = CSM-0.9;tcolor = LIGHTRED; strcpy(buf,"c");} }
    else if (val == _Shotgun_complete)
      { if (*x < CSM-1.2) {*x = CSM-1.2;tcolor = YELLOW; strcpy(buf,"sc");} }
    else if (val == _Shotgun)
      { if (*x < CSM-1.5) {*x = CSM-1.5;tcolor = GREEN; strcpy(buf,"s");} }
    else if (val == _Library_construction)
      { if (*x < CSM-1.8) {*x = CSM-1.8;tcolor = BLUE; strcpy(buf,"l");} }
    else if (val == _Dna_made)
      { if (*x < CSM-2.1) {*x = CSM-2.1;tcolor = LIGHTBLUE; strcpy(buf,"d");} }
    else if (val == _Received)
      { if (*x < CSM-2.4) {*x = CSM-2.4;tcolor = LIGHTGRAY;strcpy(buf,"r");} }
    else 
      {	if (*x < 0.0) {*x = 0.0;tcolor = WHITE; strcpy(buf,"-");} }
 
    /* what happens if tcolor is not set, hmmmm??*/
    return tcolor ;
  }      

/************************************************/ 

static void colorKeyDestroy (void) /* may be called indirectly */
{    
  colorkeyGraph = 0 ;
  graphDestroy() ;
}

static void colorKeyDismiss (void) /* explicit dismiss button */
{    
  colorKeyWanted = FALSE ;  /* this is a static bool */
  graphDestroy() ;
}

/************************************************/ 

static void fMapShowColorKey (void)
{ float x1=1.0, x2=14.0, i=1.0 ;
  Graph old = graphActive () ;

  if(graphActivate(colorkeyGraph))
    { 
    graphPop() ;
    graphActivate (old) ;
    return ;
    }

  colorkeyGraph = graphCreate (TEXT_FIT,"Sequencing status color code",0.0,0.0,0.32,0.23) ;
  /*  graphRectangle (0.5, 0.5, 75, 40) ;*/
  graphRegister (DESTROY,colorKeyDestroy) ;

  /* Display Cosmid Group Colors and Cosmid Status Colors */

  graphTextHeight (0.9) ;
  
  graphColor (BLACK) ;
  graphText ("Avery Group", x1, i+.1) ;
  /* line */
  graphLine (x2-1.5, i+.5, x2-.5,  i+.5) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.33) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.66) ;
  graphColor (BLACK) ;
  graphFillRectangle (x2, i, x2+4.0, i+1.0) ;
  graphColor (BLACK) ;
  graphText ("Archived", x2+6.5, i+.1) ;
  /* line */
  graphLine (x2+4.5, i+.5, x2+5.5, i+.5) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.33) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.66) ;

  i += 1.1 ;
  
  graphColor (BLACK) ;
  graphText ("Crick Group", x1, i+.1) ;
  /* line */
  graphLine (x2-1.5, i+.5, x2-.5,  i+.5) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.33) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.66) ;
  /* rectangle */
  graphColor (DARKRED) ;
  graphFillRectangle (x2, i, x2+4.0, i+1.0) ;
  graphColor (BLACK) ;
  graphText ("Submitted", x2+6.5, i+.1) ;
  /* line */
  graphLine (x2+4.5, i+.5, x2+5.5, i+.5) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.33) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.66) ;

  i += 1.1 ;
  
  graphColor (RED) ;
  graphFillRectangle (x2, i, x2+4.0, i+1.0) ;
  graphColor (BLACK) ;
  graphText ("Finished", x2+6.5, i+.1) ;
  /* line */
  graphLine (x2+4.5, i+.5, x2+5.5, i+.5) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.33) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.66) ;

  i += 1.1 ;
  
  graphColor (LIGHTRED) ;
  graphFillRectangle (x2, i, x2+4.0, i+1.0) ;
  graphColor (BLACK) ;
  graphText ("Contiguous", x2+6.5, i+.1) ;
  /* line */
  graphLine (x2+4.5, i+.5, x2+5.5, i+.5) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.33) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.66) ;

  i += 1.1 ;
  
  graphColor (YELLOW) ;
  graphFillRectangle (x2, i, x2+4.0, i+1.0) ;
  graphColor (BLACK) ;
  graphText ("Shotgun complete", x2+6.5, i+.1) ;
  /* line */
  graphLine (x2+4.5, i+.5, x2+5.5, i+.5) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.33) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.66) ;

  i += 1.1 ;
  
  graphColor (BLACK) ;
  graphText ("Franklin Group", x1, i+.1) ;
  /* line */
  graphLine (x2-1.5, i+.5, x2-.5,  i+.5) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.33) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.66) ;
  /* rectangle */
  graphColor (GREEN) ;
  graphFillRectangle (x2, i, x2+4.0, i+1.0) ;
  graphColor (BLACK) ;
  graphText ("Shotgun", x2+6.5, i+.1) ;
  /* line */
  graphLine (x2+4.5, i+.5, x2+5.5, i+.5) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.33) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.66) ;

  i += 1.1 ;
  
  graphColor (BLACK) ;
  graphText ("Watson Group", x1, i+.1) ;
  /* line */
  graphLine (x2-1.5, i+.5, x2-.5,  i+.5) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.33) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.66) ;
  /* rectangle */
  graphColor (BLUE) ;
  graphFillRectangle (x2, i, x2+4.0, i+1.0) ;
  graphColor (BLACK) ;
  graphText ("Library construction", x2+6.5, i+.1) ;
  /* line */
  graphLine (x2+4.5, i+.5, x2+5.5, i+.5) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.33) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.66) ;

  i += 1.1 ;

/* new mdh */
  graphColor (BLACK) ;
  graphText ("Cold Spring", x1, i+.1) ;
  /* line */
  graphLine (x2-1.5, i+.5, x2-.5,  i+.5) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.33) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.66) ;
/* new mdh */
  graphColor (LIGHTBLUE) ;
  graphFillRectangle (x2, i, x2+4.0, i+1.0) ;
  graphColor (BLACK) ;
  graphText ("Dna made", x2+6.5, i+.1) ;
  /* line */
  graphLine (x2+4.5, i+.5, x2+5.5, i+.5) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.33) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.66) ;

  i += 1.1 ;
  
/* new mdh */
  graphColor (BLACK) ;
  graphText ("ABI", x1, i+.1) ;
  /* line */
  graphLine (x2-1.5, i+.5, x2-.5,  i+.5) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.33) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.66) ;
/* new mdh */
  graphColor (LIGHTGRAY) ;
  graphFillRectangle (x2, i, x2+4.0, i+1.0) ;
  graphColor (BLACK) ;
  graphText ("Received", x2+6.5, i+.1) ;
  /* line */
  graphLine (x2+4.5, i+.5, x2+5.5, i+.5) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.33) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.66) ;

  i += 1.1 ;
  
  graphColor (BLACK) ;
  graphText ("Not Assigned", x1, i+.1) ;
  /* line */
  graphLine (x2-1.5, i+.5, x2-.5,  i+.5) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.33) ;
  graphLine (x2-.5,  i+.5, x2-1.0, i+.66) ;
  /* rectangle */
  graphRectangle (x2, i, x2+4.0, i+1.0) ;
  graphText ("Gap", x2+6.5, i+.1) ;
  /* line */
  graphLine (x2+4.5, i+.5, x2+5.5, i+.5) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.33) ;
  graphLine (x2+4.5, i+.5, x2+5.0, i+.66) ;

  i += 2.5 ;
  graphTextHeight (0.0) ;
  graphButton ("Dismiss", colorKeyDismiss, 13.5, i) ;

  graphRedraw();
  /* graphStart (colorkeyGraph); 
     NO, dont want to do that, i think, please mail me about it, jean */
  graphActivate (old) ;
}


/* code for symbolicStatus. */


static void symbolicStatus (LOOK look, float offset)
  {
  int i, j, click, dx, dy ;
  float y1 = MAP2GRAPH(look->map, look->min) ;
  float y2 = MAP2GRAPH(look->map, look->max) ;
  float x1, x2, cw, ch, tmp ;
  Array a = arrayCreate(20,BSunit) ;
  Array b = arrayCreate (8, KEYZONE) ;
  enum Colour tcolor = BLACK ;
  char text[25], buf[25] ;
  SEG *seg ;
  OBJ obj = 0 ;
  KEY val, _Nbr_gel_readings ;
  
  lexaddkey("Nbr_gel_readings", &_Nbr_gel_readings, 0) ;

  /* display max box. */
  graphColor (BLACK) ;
  graphRectangle (offset+1.0, y1, offset+CSM, y2) ;

  strcpy(text, "");
  for (i = 1 ; i < arrayMax(look->segs) ; ++i)
    {
    seg = arrp(look->segs,i,SEG) ;

    if (seg->x1 >= look->max || seg->x2 < look->min || 
	!(seg->type == SEQUENCE || seg->type == SEQUENCE_UP) ||
	arrp(look->seqInfo, seg->data.i, SEQINFO)->flags != 0)
      continue ;
    
    obj = bsCreate(seg->key) ;
    
    if (seg->x1 < look->min)
      y1 = MAP2GRAPH(look->map, look->min) ;
    else
      y1 = MAP2GRAPH(look->map, seg->x1) ;
    
    if (seg->x2 > look->max)
      y2 = MAP2GRAPH(look->map, look->max) ;
    else
      y2 = MAP2GRAPH(look->map, seg->x2) ;
    
    if (y1 < 2 + topMargin) y1 = 2 + topMargin ;
    
    if (bsFindTag (obj, _Status) && bsFlatten(obj, 3, a))
      {
      /* Find # of reads in database if available. */
      if (bsFindTag (obj, _Nbr_gel_readings))
	{
	click = 0;
	if (bsGetData (obj, _bsRight, _Int, &click))
	  sprintf (buf, "%d", click);
	else
	  strcpy (buf, "none");
	}
      else
        strcpy (buf, "none");
      
      x1 = 1.0;
      x2 = 0.0;
      for (j=0; j<arrayMax(a); j+=3)
	{
	val = arr(a, j, BSunit).k;
	tcolor = cosmid_color (val, &x2, text) ;
	}
      
      array(look->boxIndex, graphBoxStart(), int) = i;

      graphColor (tcolor) ;

      graphFillRectangle (offset+x1, y1, offset+x2, y2) ;
      graphBoxEnd () ;
      keyZoneAdd (b, seg->key, y1, y2, i) ;

      /* Text output beside rectangles. */
      graphTextHeight (.2) ;
      graphTextInfo (&dx, &dy, &cw, &ch) ;
      graphColor (BLACK) ;

      if ((tmp = offset+x2-strlen(text)*cw) < offset+x1+CSM)
        tmp = offset+x1+CSM+0.2 ;
      if (tcolor == cosmid_color(_Shotgun, &cw, text)) /* cw is dummy. */
        graphText (buf, tmp, .5*(y1+y2-ch)) ;
      graphTextHeight (0.0) ;
      }
    bsDestroy (obj) ;
    }
  arrayDestroy (a) ;
  }

/* show COSMID STATUS column */

void prepit (float y, float xoff, char *str, BOOL top) 
{
  char * ptr = str;
  char foo[2];
  size_t n = strlen(str);

  /* assign bottom location */
  float yoff  = -0.85 * (float)n;  /* y offset */
  float chunk = 0.8; /* increment for each char in y dim */
  float bord_ht = -0.9 * (float)n;

  /* assign top location if indicated. */
  if (top) {
    yoff = 0.2;
    chunk = 0.8;
    bord_ht *= -1.0;
  }

  /* draw the cosmid left/right guide posts. */
  graphColor (WHITE) ;
  graphFillRectangle (xoff, y, xoff+1.4, y+bord_ht) ;
  graphColor (BLACK) ;
  graphRectangle (xoff, y, xoff+1.4, y+bord_ht) ;
  graphColor (BLUE) ;
  graphTextHeight (0.2) ;

  /* loop through and display chars. */
  foo[1] = '\0';
  while (*ptr) {
    foo[0] = *ptr++ ;
    graphText (foo, xoff+0.26, y+yoff) ;
    yoff += chunk ;
  }
}


#define FIN_COL_W (9.0)

void fMapShowStatus (LOOK look, float *offset)
{
  SEG *seg ;
  int n = 0 ;
  float y1 = MAP2GRAPH(look->map, look->min) ;
  float y2 = MAP2GRAPH(look->map, look->max) ;

  /* fudge distance to allow for chromosome directional indicators. */
  *offset += 0.5 ;

  seg = arrp(look->segs,0,SEG) ;

  if (strstr(name(seg->key),"elegan")) {
    prepit(y1, *offset, "RIGHT", TRUE);
    prepit(y2, *offset, "LEFT", FALSE);
  } else {
    prepit(y1, *offset, "LEFT", TRUE);
    prepit(y2, *offset, "RIGHT", FALSE);
  }
  
/*
  graphColor (BLACK) ;
  graphLine (*offset+0.6, y1+3.5, *offset+0.6, y1+5.5) ;
  graphLine (*offset+0.1, y1+5.0, *offset+0.6, y1+5.5) ;
  graphLine (*offset+0.6, y1+5.5, *offset+1.1, y1+5.0) ;
*/
  
  *offset += 1.0 ;
  symbolicStatus (look, *offset) ;

  strcpy (look->segTextBuf, 
	  messprintf ("%d bases, %d analysed", 
		      look->max - look->min, n)) ;
  graphBoxDraw (look->segBox, -1, -1) ;
  graphColor (BLACK) ;

  *offset += FIN_COL_W ;
}

/* MDH cos status */
/***********************************************************************/
/***********************************************************************/

/* MDH show status , was in mdh:fmapfeatures.c*/

#define EMPTY_COLOR WHITE

static enum Colour group_color (char * buf) 
  { if ('W' == *buf || 'w' == *buf)
      return BLUE ;
    else if ('C' == *buf || 'c' == *buf)
      return DARKRED ;
    else if ('F' == *buf || 'f' == *buf)
      return GREEN ;
    else if ('A' == *buf || 'a' == *buf)
      return BLACK ;
    else if ('S' == *buf || 's' == *buf)
      return LIGHTBLUE ;
    else if ('B' == *buf || 'b' == *buf)
      return LIGHTGRAY ;
    else
      return EMPTY_COLOR ;
  }

/******************** text features *****************************/

static int bumpOff ;

static void bumpWrite (BUMP bump, char *text, int x, float y)
{
  bumpItem (bump, strlen(text)+1, 1, &x, &y) ;
  graphText (text, bumpOff+x, y-0.5) ;
}

/***********************************************************************/

void fMapShowOrigin (LOOK look, float *offset)
{
  int x, i ;
  float y1, y2, y ;
  enum Colour fcolor, pcolor ;
  OBJ  obj = 0 ;
  SEG *seg ;
  BUMP bump = 0 ;
  Array a = 0 ;
  KEY grp ;
  float myTopMargin = 2 + topMargin ;

  /* add keys needed for processing. */
  KEY _Production_Group,
      _Finishing_Group ;

  if (arrayMax(look->segs) < 2) return ;
  /*  if (!getEnv("ACEDB_SHOW_STATUS")) return ; */
  lexaddkey("Production_Group", &_Production_Group, 0) ;
  lexaddkey("Finishing_Group", &_Finishing_Group, 0) ;

  /*   if (!colorKeyWanted)  always show the button ? */
    graphButton ("Status Info", fMapShowColorKey, *offset, topMargin + .2) ;
  if (colorKeyWanted)
    fMapShowColorKey () ;
  bump = bumpCreate (3, 0) ; 
  a = arrayCreate (8, KEYZONE) ;
  /* draw the cosmid bars */
  *offset += 1 ;
  i = (look->flag & FLAG_REVERSE) ? arrayMax(look->segs)-1 : 1 ;
  while (i && i < arrayMax(look->segs))
    { seg = arrp(look->segs,i,SEG) ;
      if (seg->x1 > look->max || seg->x2 < look->min)
	goto loop1 ;
      y1 = MAP2GRAPH(look->map,seg->x1) ;
      y2 = MAP2GRAPH(look->map,seg->x2+1) ;	/* to cover full base */
      if (y1 >= mapGraphHeight) continue ;
      if (y2 <= myTopMargin) continue ;
      if (y1 < myTopMargin) y1 = myTopMargin ;
      if (y1 > mapGraphHeight) y1 = mapGraphHeight ;
      if (y2 < myTopMargin) y2 = myTopMargin ;
      if (y2 > mapGraphHeight) y2 = mapGraphHeight ;

      if ((seg->type == SEQUENCE || seg->type == SEQUENCE_UP) &&
	  !(arrp (look->seqInfo, seg->data.i, SEQINFO)->flags & SEQ_CANONICAL) && 
	  !strstr(name(seg->key),"SUPER"))
	{ x = 0 ;
	  y = seg->x1 ;		/* bump in DNA coords so OK reversed */
	  bumpItem (bump, 1, seg->x2-seg->x1+1, &x, &y) ;

	  obj = bsCreate(seg->key) ;

/*
	  if (bsFindTag (obj,_Production_Group) &&
	      bsGetData (obj, _bsRight, _Text, &buf))
	    pcolor = group_color(buf) ;
*/
	  pcolor = EMPTY_COLOR ;
	  if (bsFindTag (obj,_Production_Group) &&
              bsGetKey (obj, _Production_Group, &grp))
	    pcolor = group_color(name(grp)) ;

/*
	  if (bsFindTag (obj,_Finishing_Group) &&
	      bsGetData (obj, _bsRight, _Text, &buf))
 	    fcolor = group_color(buf) ;
*/	  
	  fcolor = EMPTY_COLOR ;
	  if (bsFindTag (obj,_Finishing_Group) &&
	      bsGetKey (obj, _Finishing_Group, &grp))
 	    fcolor = group_color(name(grp)) ;
	  
          array(look->boxIndex, graphBoxStart(), int) = i;

          /* Show who did Production on this cosmid. */
          graphColor(pcolor) ;
          if (pcolor == EMPTY_COLOR) {  /* then empty rectangle */
	    graphColor(BLACK) ;
            graphRectangle (*offset+x, y1, *offset+x+0.4, y2) ;
	    }
          else
            graphFillRectangle (*offset+x, y1, *offset+x+0.4, y2) ;

          /* Show who did Finishing on this cosmid. */
          graphColor(fcolor) ;
          if (fcolor == EMPTY_COLOR) {  /* then empty rectangle */
	    graphColor(BLACK) ;
            graphRectangle (*offset+x+0.4, y1, *offset+x+0.9, y2) ;
	    }
          else
            graphFillRectangle (*offset+x+0.4, y1, *offset+x+0.9, y2) ;
          graphBoxEnd() ;
          keyZoneAdd (a, seg->key, y1, y2, i) ;

	  graphColor(BLACK) ;
	  bsDestroy (obj) ;
	}
      loop1: if (look->flag & FLAG_REVERSE) --i ; else ++i ;
    }
  *offset += bumpMax(bump) ;
  bumpDestroy (bump) ;

  /* then the text */
  bump = bumpCreate (12, 0) ;
  bumpOff = *offset ;
  arraySort (a, keyZoneOrder) ;
  for (i = 0 ; i < arrayMax(a) ; ++i)
    { KEYZONE *z = arrp(a, i, KEYZONE) ;
      array(look->boxIndex,graphBoxStart(),int) = z->iseg ;
      bumpWrite (bump, name(z->key), 0, 0.5*(z->y1 + z->y2)) ;
      graphBoxEnd() ;
    }
  *offset += bumpMax (bump) ;
  bumpDestroy (bump) ;

  arrayDestroy (a) ;
}

/* MDH show status */
/***********************************************************************/
/***********************************************************************/

 
