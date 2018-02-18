/*  File: pepseqcol.c
 *  Author: Clive Brown (cgb@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * $Id: pepseqcol.c,v 1.4 2016/09/03 00:30:45 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 21 13:48 1998 (fw)
 * Created: Tue Nov  7 14:40:06 1995 (cgb)
 *-------------------------------------------------------------------
 */

#include "pepdisp.h" 

/************************************************************/

typedef struct selectedRes {  /* this now per pepchar rather per each seq save space and time */
  int resColour;
} RESSELECT;

typedef struct pepSeqColPriv {
  BOOL showRes;
  Array colMap;  /* is an array of two arrays . 1) the characters (pepDecodeChar)
		                                2) Colours for each of these. */
  int wrapNum;
  /* drag- active section*/
  int firstbox,lastbox,startbox,finalbox;
  float x1,x2,y1,y2;
  /* end drag bits */
  } *PEPSEQPRIV;

static void seqDragActive(COLINSTANCE instance, int box) ;

int AMINOdefaultcol[]={YELLOW,LIGHTBLUE,LIGHTRED,LIGHTRED,YELLOW,CYAN,GREEN,YELLOW,GREEN,YELLOW,YELLOW,
LIGHTRED,CYAN,LIGHTRED,GREEN,LIGHTRED,LIGHTRED,YELLOW,YELLOW,YELLOW,GRAY};
/***********************************************************************************************************************/
struct configLocals {
  BOOL showRes;
  int wrapNum;
  int colour[22];
};

static void setColourMap (char * ress,COLINSTANCE instance)
{
  int i;
  PEPSEQPRIV private = instance->private;
  
  private->colMap = arrayHandleCreate(2,Array,instance->handle);
  array(private->colMap,0,Array)  = arrayHandleCreate(strlen(ress),char,instance->handle); /* residue names */
  array(private->colMap,1,Array) = arrayHandleCreate(strlen(ress),int,instance->handle);  /* residue colour */
  for(i=0;  ress[i];i++){
    array(arr(private->colMap,0,Array),i,char) = ress[i];
    array(arr(private->colMap,1,Array),i,int) = AMINOdefaultcol[i]; 
  }
  
}


static void setColourErik()
{
  extern void graphRedoColourBoxes(void);
  int i;
  struct configLocals *cf;

  graphAssFind(&GRAPH2COLCONTROL_LOCALS_ASSOC, &cf); 

  for(i=0;i < 20;i++)
    {
      cf->colour[i] = AMINOdefaultcol[i]; 
    }
  graphRedoColourBoxes();

  return;
} /* setColourErik */



static void setColourNull()
{
  extern void graphRedoColourBoxes(void);
  int i;
  struct configLocals *cf;

  graphAssFind(&GRAPH2COLCONTROL_LOCALS_ASSOC, &cf); 
  for(i=0;i < 21;i++)
    {
      cf->colour[i] = WHITE; 
    }
  graphRedoColourBoxes();

  return;
} /* setColourNull */



static void setColourToby()
{
  extern void graphRedoColourBoxes(void);
  struct configLocals *cf;

  graphAssFind(&GRAPH2COLCONTROL_LOCALS_ASSOC, &cf); 

  cf->colour[0] =cf->colour[7] =cf->colour[10] =cf->colour[4] =cf->colour[18] 
    =cf->colour[9] =cf->colour[17] = cf->colour[1] =  MIDBLUE;
  
  cf->colour[11] =cf->colour[13] =cf->colour[15] =cf->colour[16] = GREEN;
  
  cf->colour[14] = cf->colour[8] = RED;
  
  cf->colour[3] =  cf->colour[2] = PURPLE;
  cf->colour[19] =   LIGHTBLUE;
  cf->colour[5] =  ORANGE;
  cf->colour[12] = YELLOW;
  cf->colour[6] =  LIGHTRED;
  graphRedoColourBoxes();

  return;
} /* setColourToby */

/***********************************************************************************************************************/
 
static void seqsDoColour(COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  PEPSEQPRIV private = instance->private;
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;
  int index = assInt (controlBoxRegd(instance, box)) ;
  char cp;

  if(index == look->startcol && instance == control->activeInstance){
    cp = pepDecodeChar[(int)(arr(look->pep,index, char))] ;
    graphBoxDraw(box,BLACK,activecol);
    strncpy(look->messageText, messprintf("%c",cp) ,100);
    strcat(look->messageText, messprintf(" [%d] ",index));
    if(look->messageBox)
      graphBoxDraw(look->messageBox,-1,-1);
    control->activeBox = box;
  }
  else if(instance == control->activeInstance){
    if(!private->showRes) graphBoxDraw(box,BLACK,WHITE);
    else graphBoxDraw(box,graphContrastLookup[array(arr(private->colMap,1,Array),arr(look->pep,index,char),int)]
		      ,array(arr(private->colMap,1,Array),arr(look->pep,index,char),int));
  }
  else{
    if(look->endcol == 0){
      if(!private->showRes) graphBoxDraw(box,BLACK,WHITE);
      else graphBoxDraw(box,graphContrastLookup[array(arr(private->colMap,1,Array),arr(look->pep,index,char),int)]
			,array(arr(private->colMap,1,Array),arr(look->pep,index,char),int));
    }
    else if( index >= look->startcol &&  index <= look->endcol)
      graphBoxDraw(box,BLACK,friendcol);
    else{
      if(!private->showRes) graphBoxDraw(box,BLACK,WHITE);
      else graphBoxDraw(box,graphContrastLookup[array(arr(private->colMap,1,Array),arr(look->pep,index,char),int)]
			,array(arr(private->colMap,1,Array),arr(look->pep,index,char),int));
    }
  }      
}


static BOOL seqUnSelect(COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  PEPLOOK *look = map->look;
  COLCONTROL control = map->control;

  control->activeKey = 0;
  look->startcol = 0;
  look->endcol = 0;

  *look->messageText = 0;
  
  if(look->messageBox)
    graphBoxDraw(look->messageBox,-1,-1);
  
  return FALSE;
}

static BOOL seqSelect(COLINSTANCE instance, int box, double x, double y)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;

  if(box!=0)
    seqDragActive(instance, box);
  look->startcol = assInt (controlBoxRegd(instance, box)) ;
  look->endcol = 0;
  control->activeKey = 0;
  control->activeBox = box;
  return FALSE;
}
/********************************* ??? Drag-Active ??? ******************************************************/


static COLINSTANCE activeInstance;

static int findNewBox(int startbox,float x,float y)
{
  PEPSEQPRIV private = activeInstance->private;
  float b1=0.0,b2=0.0,b3=0.0,b4=0.0,bold=0.0;
  int temp;
  BOOL okay = TRUE;

  graphBoxDim(startbox,&private->x1,&private->y1,&private->x2,&private->y2);
  if(y > private->y1 && y < private->y2 && x < private->x2 && x> private->x1)
    return startbox;
  temp = startbox;
  if(y > private->y2){
    while(okay){
      temp += private->wrapNum;
      if(temp <= private->finalbox && temp >= private->firstbox)
	graphBoxDim(temp,&b1,&b2,&b3,&b4);
      else
	okay = FALSE;
      if(b2 > y)
	okay = FALSE;
    }
    temp -=private->wrapNum;
  }
  else if(y < private->y1){
    while(okay){
      temp -= private->wrapNum;
      if(temp <= private->finalbox && temp >= private->firstbox)
	graphBoxDim(temp,&b1,&b2,&b3,&b4);
      else
	okay = FALSE;
      if(b4 < y)
	okay = FALSE;
    }
    temp +=private->wrapNum;
  }
  if(x > private->x2){
    okay = TRUE;
    graphBoxDim(temp,&b1,&b2,&b3,&b4);
    bold = b2;
    while(okay){
      temp++;
      if(temp <= private->finalbox && temp >= private->firstbox)
	graphBoxDim(temp,&b1,&b2,&b3,&b4);
      else
	okay = FALSE;
      if(b1 > x || b2 != bold)
	okay = FALSE;
      else
	bold = b2;
    }
    temp--;
  }
  else if(x < private->x1){
    okay = TRUE;
    graphBoxDim(temp,&b1,&b2,&b3,&b4);
    bold= b2;
    while(okay){
      temp--;
      if(temp <= private->finalbox && temp >= private->firstbox)
	graphBoxDim(temp,&b1,&b2,&b3,&b4);
      else
	okay = FALSE;
      if(b3 < x || b2 != bold)
	okay = FALSE;
      else
	bold = b2;
    }
    temp++;
  }
  return temp;
}

static void controlLeftDrag(double x, double y)
{
  MAPCONTROL map = activeInstance->map;
  PEPLOOK *look = map->look;
  PEPSEQPRIV private = activeInstance->private;
  int incr,colour,newbox,i,index;
  BOOL inlimits = TRUE;

  newbox = findNewBox(private->lastbox,x,y); 

  if(newbox != private->lastbox){
    i = private->lastbox;
    if(private->lastbox < newbox )
      incr = 1;
    else
      incr = -1;
    while(inlimits){
      if((i>=private->startbox && i<= newbox) || (i <= private->startbox && i>= newbox)){
	colour = LIGHTGRAY;
      }
      else{
	index = assInt (controlBoxRegd(activeInstance, i)) ;
	if(private->showRes)
	  colour= array(arr(private->colMap,1,Array),arr(look->pep,index,char),int);
	else
	  colour = WHITE;
      }
      graphBoxDraw(i,BLACK,colour);
      if(i==newbox)
	inlimits= FALSE;
      else
	i+=incr;
      if( i > private->finalbox || i <private->firstbox)
	inlimits= FALSE;
    }
    private->lastbox = i;
  }
}


static void controlLeftUp(double x, double y){
  PEPSEQPRIV private = activeInstance->private;
  PEPLOOK *look = (PEPLOOK*) activeInstance->map->look;

  graphRegister(LEFT_DRAG, 0);
  graphRegister(LEFT_UP, 0);
  if(private->startbox!=private->lastbox){
    if(private->startbox < private->lastbox){
      look->activeStart = assInt (controlBoxRegd(activeInstance,private->startbox)) ;
      look->activeEnd = assInt (controlBoxRegd(activeInstance, private->lastbox)) ;
    }
    else{
      look->activeStart = assInt (controlBoxRegd(activeInstance,private->lastbox)) ;
      look->activeEnd = assInt (controlBoxRegd(activeInstance, private->startbox)) ;
    }
    strncpy(look->activeText,messprintf("%d %d",look->activeStart,look->activeEnd),11);

    if(look->activeregionBox)
      graphBoxDraw(look->activeregionBox,-1,-1);
    controlDraw();
  }
}

static void seqDragActive(COLINSTANCE instance, int box)
{
  PEPSEQPRIV private = instance->private;

  activeInstance = instance;
  private->startbox = private->lastbox = box;
  graphRegister(LEFT_DRAG, controlLeftDrag);
  graphRegister(LEFT_UP,controlLeftUp);
}
/************************************ End active Zone ***************************************************/

/***********************************************************************************************************************/
static BOOL setBoxes(COLINSTANCE instance, int i, int j, float *offset, int *offFact, float y)
{
 PEPLOOK *look = (PEPLOOK*) instance->map->look;
 char cp;

 cp = pepDecodeChar[(int)(arr(look->pep, i, char))] ;
    
 if(instance->map->mag < 2.6)
   graphText (messprintf("%c",cp), j, y - 0.6) ;
 else
   graphText (messprintf("%c",cp), *offset, y - 0.5) ;

 return(TRUE);
}

/*****************************************************************************************************/

static void pepSequenceDraw (COLINSTANCE instance, float *offset)
{
  PEPLOOK *look = (PEPLOOK*) instance->map->look ;
  COLCONTROL control = instance->map->control ;	/*  */
  PEPSEQPRIV private = instance->private; 
  int i, j, offFact; 
  int row,box;
  float start;
  BOOL first = TRUE;

  if (!look->pep)
    return ;


  /* Calculate the number of seq per line. Auto shrink facility */
  i=1;
  start = MAP2GRAPH(instance->map, i);
  while(MAP2GRAPH(instance->map,i)-start <  2 )
    i++;
  offFact = i-2;
  if(offFact < 1)
    offFact = 1;
  private->wrapNum = (offFact+1)/2;

  /* offFact now contains the number of items per line allowed */

  for (i = 0 ; i < arrayMax(look->pep);)
    { 
      if (offFact > 0) /* mieg, otherwise you may loop */
	{ 
	  for(j= *offset, row=i; j< *offset + offFact  &&  i < arrayMax(look->pep); j+=2, i++)
	    {
	      box =  graphBoxStart();
	      if(MAP2GRAPH(instance->map, i) >= control->topMargin+0.5 && MAP2GRAPH(instance->map,i) < control->graphHeight-0.5)  
		{  
		  if(first){ /* first and last box needed for the drag select */
		    private->firstbox = box;
		    first = FALSE;
		  }
		  private->finalbox = box;
		  setBoxes(instance, i, j ,offset, &offFact, 
			   (MAP2GRAPH(instance->map, row))); /*helix taken out as requisted */
		  controlRegBox (instance, box, assVoid (i)) ;
		  if(private->showRes)
		    graphBoxDraw(box, BLACK, array(arr(private->colMap,1,Array),arr(look->pep,i,char),int));
		  else
		   graphBoxDraw(box, BLACK,WHITE); 
		}
	      graphBoxEnd();
	    }	  
	} 
    }
  
   *offset += offFact+2;
  
}
/***********************************************************************************************************************/

static BOOL pepSeqConfig(COLINSTANCE instance)
{
  PEPSEQPRIV private = instance->private;
  struct  configLocals *cf = (struct configLocals *) messalloc(sizeof(struct configLocals));
  float line = 2.0,x,y;
  int i,maxX,maxY;
  static MENUOPT pepseqcolourMenu[] = {
    { setColourErik,"Erik's" },
    { setColourToby,"Toby's" },
    { setColourNull,"Null"},
    { 0, 0 }
  } ;
  

  if(controlCreateConfig(instance,cf,"Peptide Sequence Configure",0.8,0.3)){
    
    /* add the default colout button */
    i = graphButton("Default Colours...",setColourErik,35.0,0.5);
    graphBoxMenu (i,pepseqcolourMenu);


    /* initailise the data structure */
    cf->showRes = private->showRes;
    cf->wrapNum = private->wrapNum;
    for(i=0;i<arrayMax(arr(private->colMap,1,Array))-1;i++){
      cf->colour[i] = array(arr(private->colMap,1,Array),i,int);
    }

    /* draw the configuration */
    graphToggleEditor("Colour Residues",&cf->showRes,4.0,line++);
    graphIntEditor("Residues per wrap:",&cf->wrapNum,4.0,line++,0);
    graphFitBounds(&maxX,&maxY);
    maxY -= 4;
    y = 0;x= 33; 
    for(i=0;i</*arrayMax(arr(private->colMap,0,Array))-1*/21;i++){
      y+=2;
      if(y>maxY){
	y=2.0;
	x+= 10;
      }
      graphColourEditor(pepShortName[(int)array(arr(private->colMap,0,Array),i,char)]
			,messprintf("%c",array(arr(private->colMap,0,Array),i,char)),&cf->colour[i],x,y);
    }
    graphRedraw();
    
    return FALSE; /* i do not want to redraw the pepdisplay screen YET */
  }
  else
    return FALSE;
}

static void seqConfigFinal(COLINSTANCE instance, void *locals,BOOL ok)
{ struct configLocals *cf = locals;
  PEPSEQPRIV private = instance->private;
  MAPCONTROL map = instance->map;
  int i;

  if(ok){
    if (cf->wrapNum <= 0) /* mhmp 10.12.98 */
      cf->wrapNum = private->wrapNum;    
    private->showRes = cf->showRes;
     map->mag *=(float)((float)private->wrapNum/(float)cf->wrapNum);
     private->wrapNum = cf->wrapNum; 
     for(i=0;i<arrayMax(arr(private->colMap,1,Array))-1;i++){
       array(arr(private->colMap,1,Array),i,int) = cf->colour[i] ;
     }
   }
  else
    messfree(cf);
}

/****************************************************************************************************************/
static void pepSeqSave(COLINSTANCE instance, OBJ init)
{
  PEPSEQPRIV private = instance->private;
  int i;

  if(private->showRes)
    bsAddTag(init,str2tag("PS_Highlight_residue"));

  if(private->wrapNum)
    bsAddData(init,str2tag("PS_Residues_per_wrap"),_Int,&private->wrapNum);

  for(i=0;i<arrayMax(arr(private->colMap,0,Array))-1;i++){
    if(AMINOdefaultcol[i] !=array(arr(private->colMap,1,Array),i,int)){
      if (bsAddData (init, str2tag("PS_Colours"), _Text, pepShortName[(int)(array(arr(private->colMap,0,Array),i,char))])){
	controlSetColour(init,array(arr(private->colMap,1,Array),i,int));
      }
      else{
	 messerror ("Invalid pepSequence model - needs PS_Colours Text #Colour") ;
	 printf("%d %s\n",i,pepShortName[(int)(array(arr(private->colMap,0,Array),i,char))]);
	 break ;
       }
    }
  }
}
/************************************************************************************************************************/
static int findindex(COLINSTANCE instance,char *str1){
  PEPSEQPRIV private = instance->private;
  char str2[4];
  int i;

  for(i=0;i<strlen(pepDecodeChar);i++){
    sprintf(str2,"%s",pepShortName[(int)array(arr(private->colMap,0,Array),i,char)]);
    if(strcmp(str1,str2)==0)
      return i;
  }
  return -1;
}

extern BOOL pepSequenceCreate (COLINSTANCE instance, OBJ init)
{ 
  PEPSEQPRIV private;
  PEPLOOK *look = (PEPLOOK *)instance->map->look;
  int i,k=0;
  Array flatA;
  char *temp;
		    
  instance->draw = pepSequenceDraw ; 
  instance->configure = pepSeqConfig;
  instance->setSelectBox = seqSelect;
  instance->unSelectBox =seqUnSelect;
  instance->doColour = seqsDoColour;
  instance->save = pepSeqSave;
  instance->configFinal = seqConfigFinal;

  private = (PEPSEQPRIV)halloc(sizeof(struct pepSeqColPriv), instance->handle);
  private->wrapNum = DEF_WRAP;
  private->showRes = FALSE;
  instance->private = private;

  setColourMap (pepDecodeChar,instance);
 
  look->endcol = 0;

  instance->map->hasProjectionLines = FALSE;

  if(init){
    if(bsFindTag(init,str2tag("PS_Highlight_residue")))
      private->showRes = TRUE;
    if(bsGetData(init, str2tag("PS_Residues_per_wrap"), _Int,&i))
      private->wrapNum = i;
        
    flatA = arrayCreate(3, BSunit);
    if(bsFindTag(init, str2tag("PS_Colours")) && bsFlatten(init,3,flatA)){
      for(i=0;i<arrayMax(flatA);i+=3){
	temp = arr(flatA,i,BSunit).s;
	k = findindex(instance, temp);
	if(k >= 0)
	  array(arr(private->colMap,1,Array),k,int) = arr(flatA,i+1,BSunit).i-_WHITE;
	else
	  messerror ("could not find %s in Amino Acid Table's. Please check.",temp) ;
      }
    }
    arrayDestroy(flatA);
  }
    
  
  return TRUE ;
}


/***********************************************************************************************************************/








 
 
 
 
 
 
 
 
