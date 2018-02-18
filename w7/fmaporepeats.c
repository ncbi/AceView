
#include "acedb.h"
#include "graph.h"
#include "fmap_.h"
#include "dna.h"
#include "a.h"
#include "lex.h"
#include "../whooks/sysclass.h"
#include "../whooks/classes.h"


/*********************************************************************/

void fMapcDNAShowDavidRepeats (LOOK look, float *offset)
{
  KEY rKey ;
  SEG *seg ;

  /* internal cuisine */
  unsigned char *ip ;
  int n, n1, n2, i0, ii, j, dx ;
  float x = *offset, y = 0, dy = 0 ;
  Array aa = 0 ; 

  /* choose a step with a resolution which is not too high */

  dy = MAP2GRAPH (look->map, 1) - MAP2GRAPH (look->map, 0) ;
  if (dy > 0) dx = 1/dy ;
  else return ;

  if (dx < 1) dx = 1 ;
  if (dx > 100) return ;

  y = MAP2GRAPH (look->map, look->min) ;
  dy = MAP2GRAPH (look->map, dx) - MAP2GRAPH (look->map, 0) ;
  
  /* find the offset and load the array */
  for (i0 = -1, ii = 1, seg = arrp(look->segs,ii,SEG)  ; ii < arrayMax(look->segs) ; seg++, ii++)
    { 
      if (seg->type == SEQUENCE &&
	  seg->x1 <= look->min && seg->x2 >= look->max &&
	  lexword2key (messprintf("%s.i100", name(seg->key)), &rKey, _VOligoRepeat) &&
	  (aa = arrayGet (rKey, unsigned char, "c")))
	{ i0 = seg->x1 ; break ; }
    }
  if (!aa)
    return ;

  /* draw the column */
  for (ii = look->min ; ii < look->max ; ii += dx, y += dy)
    {
      n = n1 = n2 = 0 ;
      if (ii >= i0 && ii + dx -i0 < arrayMax(aa))
	{
	  ip = arrp (aa, ii - i0, unsigned char) ;
	  for (j = 0 ; j < dx ; ip++, j++) 
	    { n = *ip ; if (n1 < n) n1 = n ; n2 += n ; }
	  graphText (messprintf ("%3d %5d", n1, n2), x, y) ;
	}
    }
  arrayDestroy (aa) ;

  *offset += 10 ;
}

/*********************************************************************/

void fMapcDNAShowOligoRepeats (LOOK look, float *offset)
{
  /* parameters controlling the drawing */
  int saturation = 8 ; /* number of repeats at which we saturate */
  float scale = 4 ;   /* maximal horizontal elongation in char units */
  float yStep = .5 ;  /* max vertical density in char units */
  int oColorMax = 4, oColor [] = { BLUE, GREEN, ORANGE, RED } ; 

  /* internal cuisine */
  int ii, iiMax, iol, iolMax, ol, dx, phase = 0 ;
  int jj, nj, xj ; /* counters to group values at high zoom */
  float x = *offset, xx, y = 0, dy = 0, maxHeight = 0 ;
  Array aa = 0 ; int iaa ;   /* the lineArray */
  unsigned char *pp ;
  
  if (1)
    {
      fMapcDNAShowDavidRepeats (look, offset) ;
      return ;
    }
  if (! look->oligoRepeats || 
      ! arrayExists (look->oligoRepeats->oLengths) || 
      ! arrayMax (look->oligoRepeats->oLengths) ||
      arrayMax (look->oligoRepeats->profile) < look->min ||
      *look->map->activeColName == '-')
    return ;

  iolMax = arrayMax (look->oligoRepeats->oLengths) ;

  for (iol = 0 ; iol < iolMax ; iol++)
    {
      ol = arr (look->oligoRepeats->oLengths, iol, int) ;
      
      /* phase is needed to select points for a given oligo lenth */
      phase = ((look->oligoRepeats->phase - look->min) % iolMax + iolMax ) %iolMax;

      /* choose a step with a resolution which is not too high */
      dx = 0 ;
      dy = MAP2GRAPH (look->map, 1) - MAP2GRAPH (look->map, 0) ;
      if (dy > 0) dx = yStep/dy ;
      dx = dx/iolMax ; dx = dx * iolMax ;
      if (dx < iolMax)
	dx = iolMax ;
      nj = dx/iolMax ;
      y = MAP2GRAPH (look->map, look->min + phase + ol/2) ;
      dy = MAP2GRAPH (look->map, dx) - MAP2GRAPH (look->map, 0) ;

      aa = arrayReCreate (aa, 2*(look->max - look->min)/dx + 8, float) ;
      array (aa, 2*(look->max - look->min)/dx + 7, float) = 0 ; /* make room */
      iiMax = arrayMax (look->oligoRepeats->profile) ;
      xj = 0 ;
      for (iaa = 0, jj = 0, ii = look->min + phase + iol, pp = arrp (look->oligoRepeats->profile, ii, unsigned char) ;
	   ii < look->max && ii < iiMax ; ii += iolMax, pp += iolMax)
	{
	  if (*pp > xj) xj = *pp ;
	  if (++jj == nj)
	    {
	      jj = 0 ;
	      if (xj > saturation) xj = saturation ;
	      xx = xj * scale / saturation ;
	      if (xx > maxHeight)  maxHeight = xx ;
	      array (aa, iaa++, float) = x + xx ;
	      array (aa, iaa++, float) = y ;
	      y += dy ; xj = 0 ;
	    }
	}
      if (ii < look->max)
	{
	  array (aa, iaa++, float) = x ;
	  array (aa, iaa++, float) = y ;
	  array (aa, iaa++, float) = x ;
	  array (aa, iaa++, float) = MAP2GRAPH (look->map, look->max) ;
	}
      if (iaa > 2*(look->max - look->min)/dx + 8)
	messcrash ("Jean has to eat his hat") ;
      
      arrayMax (aa) = iaa ;
      graphColor (oColor[iol % oColorMax]) ;
      graphLineSegs (aa) ;
    }
  graphColor (BLACK) ;
  arrayDestroy (aa) ;
  *offset += maxHeight + 1 ;
}

/*********************************************************************/
/*********************************************************************/


