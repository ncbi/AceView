/*  File: pepgraphcol.c
 *  Author: Clive Brown (cgb@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: pepgraphcol.c,v 1.1.1.1 2002/07/19 20:23:00 sienkiew Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 17 17:08 1998 (fw)
 * Created: Tue Nov  7 14:41:12 1995 (cgb)
 *-------------------------------------------------------------------
 */

#include "pepdisp.h" 

typedef struct HydroPrivStruct {
  KEY pointKey, intervalKey;
  int win, colWidth;
  int mainBox;
  int scaleBox;
  BOOL showZero,fixScale;
  BOOL refresh;
  Array fitVals;
  float minval,maxval,zeropos;
  } *HYDROPRIV;

/********************************************** hydrophobicity column **************************************************/

static void calcHydroph (Array hydroVals, Array pep, int win)
{
  float sum = 0 ; /* mieg : was not initialised */
  int i,j;
  static float hydrophob[] = { 1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
	                       1.9, -3.5, -1.6,-3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3, -0.49 } ;

  assert(win > 1) ;
  assert(pep) ; assert (hydroVals) ;

  for(i=0; i< arrayMax(pep) + win; i++)
    {
      for(j = i ; j < i+win; j++)
	if (j < arrayMax(pep))
	  sum += hydrophob[(int)array(pep, j, char)];
     
      if (win > 0)  /* mieg */
	array(hydroVals, i+(win/2) ,float) = sum/win;
      sum =0;
    }
    
}
/***********************************************************************************************************************/
static float xposFit (Array inVals, float xmin, float xmax,float * plotmin, float * plotmax)
{
 int i;
 float minval= BIG_FLOAT,maxval=SMALL,unit;
 
 assert(xmin <= xmax);
 
 for(i=0;i< arrayMax(inVals);i++)
   {
     if(minval > array(inVals,i,float)) minval = array(inVals,i,float);
     if(maxval < array(inVals,i,float)) maxval = array(inVals,i,float);
   }
 assert (maxval >= minval);
 *plotmin = minval;
 *plotmax = maxval;

 unit = (xmax-xmin) / (((maxval - minval)==0)? SMALL : maxval-minval);
 
 /* normalise */
 
 for(i=0;i< arrayMax(inVals);i++)
   array(inVals,i,float) -= minval;
 
 for(i=0; i< arrayMax(inVals);i++)
   array(inVals,i,float) = xmin + (unit * (array(inVals,i,float)));
 
 return(xmin + (unit * (0 - minval)));
}

/***********************************************************************************************************************/
static float xposFitFixed (Array inVals, float xmin, float xmax,float * plotmin, float * plotmax)
{
 int i;
 float minval= BIG_FLOAT,maxval=SMALL,unit;
 
 assert(xmin <= xmax);
 
 minval = -2;
 maxval = 2;
 assert (maxval >= minval);
 *plotmin = minval;
 *plotmax = maxval;

 unit = (xmax-xmin) / (((maxval - minval)==0)? SMALL : maxval-minval);

 for(i=0;i< arrayMax(inVals);i++)
   array(inVals,i,float) -= minval;
 
 for(i=0; i< arrayMax(inVals);i++)
   array(inVals,i,float) = xmin + (unit * (array(inVals,i,float)));
 
 return(xmin + (unit * (0 - minval)));
}

/**********************************************************************************************************************/
static int hydroResetVals (COLINSTANCE instance, float offset)
{
  PEPLOOK *look = (PEPLOOK*) instance->map->look ;
  register HYDROPRIV private = instance->private;

  private->fitVals = arrayReCreate(private->fitVals, 1, float);  
  calcHydroph (private->fitVals, look->pep, private->win); 
  if (private->fixScale)
    private->zeropos = xposFitFixed(private->fitVals, offset, offset+ private->colWidth, &private->minval, &private->maxval);
  else
    private->zeropos = xposFit(private->fitVals, offset, offset+ private->colWidth, &private->minval, &private->maxval);
  private->refresh = FALSE;

  return 0 ;
}
/***********************************************************************************************************************/
static void hydrophobDraw (register COLINSTANCE instance, float *offset)
{
  PEPLOOK *look = (PEPLOOK*) instance->map->look ;
  COLCONTROL control = instance->map->control ;
  HYDROPRIV private = instance->private;
  register int i;
  float x, y;
  register float oldx , oldy ;
  float top;
  int colStart = *offset+3, colStop = *offset+3+private->colWidth;
  char numBuff[10];
  
  if (!look->pep)
    return ;
  
  hydroResetVals (instance , colStart);
  
  oldx = private->zeropos;
  
  top = MAP2GRAPH(instance->map, control->topMargin+0.5);
  if(top < control->topMargin+0.5) top = control->topMargin+0.5;
  oldy = top;
  graphColor(BLACK);
  graphLine(colStart,  top, colStart, MAP2GRAPH(instance->map, arrayMax(look->pep))) ;
  graphColor(BLACK);
  graphLine(colStop,  top , colStop, MAP2GRAPH(instance->map, arrayMax(look->pep))) ;
  
  graphTextHeight(0.1);
  
  if(private->showZero)
    {
      graphText ("0",private->zeropos, top);
      graphColor(GREEN);
      graphLine(private->zeropos, top, private->zeropos, MAP2GRAPH(instance->map, arrayMax(look->pep))) ;
      graphColor (BLACK);
    }

  if (sprintf(numBuff, "%.1f", private->minval))
    graphText (numBuff, colStart-2, top-1);
  if (sprintf(numBuff, "%.1f", private->maxval))
    graphText (numBuff, colStop-1, top-1);
  graphTextHeight(0);
  
  graphColor(RED);
  for (i = 0; i < arrayMax(look->pep); i++)
    { 
      x = array(private->fitVals, i, float);
      y = MAP2GRAPH(instance->map, i)+0.5 ;  /*Hmm Not quite accurate enough */
      if(y < top) y = top;
      graphLine ( oldx, oldy, x, y) ;
      oldx = x ; 
      oldy = y ;
    }
  graphColor(BLACK);
  
  *offset += private->colWidth + 5;
  
  
}
/***********************************************************************************************************************/

struct configLocals {
  BOOL showZeroBar;
  BOOL fixScale;
  int win;
  int colWidth;
};

static BOOL pepGraphConfig(COLINSTANCE instance)
{
  HYDROPRIV private = instance->private;
  struct  configLocals *cf = (struct configLocals *) messalloc(sizeof(struct configLocals));
  float line = 2.0;

  if(controlCreateConfig(instance,cf,"Hydrophobicity Configure",0.5,0.15)){
    
    cf->showZeroBar = private->showZero;
    cf->fixScale = private->fixScale;
    cf->win = private->win;
    cf->colWidth = private->colWidth;
    
    graphToggleEditor("Show Zero Bar",&cf->showZeroBar,4.0,line++);
    graphToggleEditor("Fixed Scaling",&cf->fixScale,4.0,line++);
    graphIntEditor("Hydrophobicity Calculation Window :",&cf->win,4.0,line++,0);
    graphIntEditor("Hydrophobicity Display Width :",&cf->colWidth,4.0,line++,0);

    graphRedraw();

    return FALSE; /* i do not want to redraw the pepdisplay screen YET */
  }
  else
    return FALSE;
}

static void pepConfigFinal(COLINSTANCE instance, void *locals,BOOL ok)
{ struct configLocals *cf = locals;
  HYDROPRIV private = instance->private;  
  
  if (ok)
    { 
      if (cf->win <= 1) /* mhmp 11.12.98 */
	cf->win = private->win;
      if (cf->colWidth < 0)
	cf->colWidth = private->colWidth ;
      private->showZero = cf->showZeroBar;
      private->fixScale = cf->fixScale;
      private->win = cf->win;
      private->colWidth = cf->colWidth;
    }
  else
    messfree(cf);
}

/******************************************************************************************************************/

static void hydroPrivDestroy(void *p)
/* Block finalisation function */
{
 HYDROPRIV private = (HYDROPRIV)p;

 if arrayExists(private->fitVals) arrayDestroy(private->fitVals);

}
/******************************************************************************************************************/
static void pepGraphSave(COLINSTANCE instance, OBJ init)
{
  HYDROPRIV private = instance->private;

  if(private->showZero)
    bsAddTag(init,str2tag("HP_Show_Zero_bar"));
  if(private->fixScale)
    bsAddTag(init,str2tag("HP_Fixed_Scaling"));
  if(private->win)
    bsAddData(init,str2tag("HP_Calculation_window"),_Int,&private->win);
  if(private->colWidth)
    bsAddData(init,str2tag("HP_Display_width"),_Int,&private->colWidth);
}
/******************************************************************************************************************/
extern  BOOL hydrophobCreate (COLINSTANCE instance, OBJ init)
{
  HYDROPRIV private;
  int i;
  
  instance->draw = hydrophobDraw ; 
  instance->configure = /*hydroConfigure*/pepGraphConfig;
  instance->save = pepGraphSave;
  instance->configFinal = pepConfigFinal;
 
  private = (HYDROPRIV)halloc(sizeof(struct HydroPrivStruct), instance->handle);
  blockSetFinalise (private, hydroPrivDestroy);
  
  private->fitVals = arrayCreate(10,float);
  private->win = 10; 
  private->showZero = FALSE;
  private->zeropos = 0;
  private->refresh = TRUE;
  private->colWidth = 4; 
  private->fixScale = FALSE;
  instance->private = private;

  if(init){
    if(bsFindTag(init,str2tag("HP_Show_Zero_bar")))
      private->showZero = TRUE;
    if(bsFindTag(init,str2tag("HP_Fixed_Scaling")))
      private->fixScale = TRUE;
    if(bsGetData(init, str2tag("HP_Calculation_window"), _Int,&i))
      private->win = i;
    if(bsGetData(init, str2tag("HP_Display_width"), _Int,&i))
      private->colWidth = i;

  }

  return TRUE ;
}


 
 
