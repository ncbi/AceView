/*  File: pephomolcol.c
 *  Author: Clive Brown (cgb@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1995
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: pephomolcol.c,v 1.7 2017/02/15 20:36:41 mieg Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 21 13:41 1998 (fw)
 *		-	MENU_FUNC_CAST instead of (GraphFunc) in belvuMenu[] MENUOPTS
 * * May 20 18:48 1996 (esr)
 * Created: Tue Nov  7 14:40:54 1995 (cgb)
 *-------------------------------------------------------------------
 */

#include "pepdisp.h"
#include "call.h"
#include "bump.h"

Stack BelvuMatchStack ;

typedef struct homolColPriv {
  Array segs;
  int homolBox; /* to diff between homol of the same key */
  float x;
  BOOL bump;
  int width;
  } *HOMOLPRIV;

BOOL alignDisplay (KEY key, KEY from, BOOL isOldGraph);

/*********** Configure code for the pep cols.. ****************/

struct configLocals{
  BOOL bump;
};

static BOOL homolConfigure (COLINSTANCE instance)
     
{
  HOMOLPRIV private = instance->private;
  struct configLocals *cf = (struct configLocals *) messalloc (sizeof (struct configLocals));
  float line = 2.0;

  if (controlCreateConfig (instance,cf,"Configure pepdisp name column",0.5,0.15)){
    /* initialise data */
    cf->bump = private->bump;

    graphToggleEditor ("Bump",&cf->bump,4.0,line);
  
    graphRedraw ();
  }
  return FALSE;

}

static void pepConfigFinal (COLINSTANCE instance, void *locals, BOOL ok){
  struct configLocals *cf = locals;
  HOMOLPRIV private = instance->private;

  if (ok){
    private->bump = cf->bump;
  }
  else
    messfree (cf);
}


/************************ homol columns ******************************/

typedef struct { 
  KEY key ;		/* subject, i.e. target for match */
  KEY meth ;		/* type of match, e.g. BLASTP - controls display etc. */
  float score ;
  int qstart ;		/* query start and end, i.e. in self */
  int qend ;
  int sstart ;		/* subject start and end */
  int send ;
} HOMOL ;

static int homolOrder (const void *va, const void *vb)
{
  const HOMOL *a = (const HOMOL*) va, *b = (const HOMOL*) vb ;
  return ((a->qstart + a->qend)/2 > (b->qstart + b->qend)/2) ? 1 : -1 ;
}

extern void * homolConvert (MAPCONTROL map, void *params)
{ 
  OBJ obj ;
  Array flatA, homols ;
  HOMOL *homol ;
  int i ;
 
  homols = arrayHandleCreate (8, HOMOL, map->handle) ;
  
  flatA = arrayCreate (8, BSunit) ;
  obj = bsCreate (map->key);
  
  if (bsFindTag (obj, _Homol) && bsFlatten (obj, 8, flatA))   /* flatten tree in steps of 8 */
    for (i = 0 ; i < arrayMax (flatA) ; i += 8) 
      {
	if (arr (flatA, i+2, BSunit).k != 0){/* will need to something else ??? temp fix */
	  homol = arrayp (homols, arrayMax (homols), HOMOL) ;
	  homol->key = arr (flatA, i+1, BSunit).k ;
	  homol->meth = arr (flatA, i+2, BSunit).k ;
	  homol->score  = arr (flatA, i+3, BSunit).f ;
	  homol->qstart = arr (flatA, i+4, BSunit).i ;
	  homol->qend   = arr (flatA, i+5, BSunit).i ;
	  homol->sstart = arr (flatA, i+6, BSunit).i ;
	  homol->send   = arr (flatA, i+7, BSunit).i ;
	  methodAdd (0, homol->meth) ;
	}
      } 

  bsDestroy (obj) ;
  arrayDestroy (flatA) ;
  arraySort (homols, homolOrder); /* just trying it */
  return homols ;
}

static void pepMapShowText (int box) 
{
  MAPCONTROL map;
  HOMOL *homol;

  map = currentMapControl ();
  homol = arr (map->control->boxIndex2, box, void *);

  display (homol->key,map->key,TREE);
}

/****************************************************************************************************************/

static void bsPushIntData (OBJ obj, int queryInt)
{
  int searchInt;

  bsGetData (obj, _bsRight, _Int, &searchInt);
  while (searchInt != queryInt && bsGetData (obj, _bsDown, _Int, &searchInt));
  if (searchInt != queryInt) messcrash ("Could not find Int %d", queryInt);
  bsPushObj (obj);
}


/****************************************************************************************************************/
static void pepMapPfamWWW (int box) 
{
  MAPCONTROL map;
  HOMOL *homol;

  map = currentMapControl ();
  homol = arr (map->control->boxIndex2, box, void *);

  callScript ("PfamWWW_script", name (homol->key)) ;
}


/****************************************************************************************************************/
static void pepMapShowBelvu (int box) 
{
  MAPCONTROL map;
  PEPLOOK *look;
  HOMOL *homol;
  OBJ obj;
  KEY key, _Align_homol, _Alignment, _Segs;
  char *cp, *cpep, *pep;
  int i, searchInt;
  float searchFloat;

  map = currentMapControl ();
  look = (PEPLOOK*)map->look;
  homol = arr (map->control->boxIndex2, box, void *);
  obj = bsCreate (map->key);

  BelvuMatchStack = stackReCreate (BelvuMatchStack, 10);

  /* Pass the stuff below on to Belvu */
  pushText (BelvuMatchStack, 
	   messprintf ("%s/%d-%d %.1f", name (map->key), homol->qstart, homol->qend, homol->score));

  /* The sequence segment */
  pep = messalloc (homol->qend+1);
  cpep = pep;
  cp = arrp (look->pep, 0, char);
  for (i = homol->qstart; i <= homol->qend; i++, cpep++)
      *cpep = pepDecodeChar[ (int)* (cp+i-1)];
  pushText (BelvuMatchStack, pep);
  messfree (pep);
  
  /* The Segs */
  lexaddkey ("Alignment", &_Alignment, 0);
  lexaddkey ("Align_homol", &_Align_homol, 0);
  lexaddkey ("Segs", &_Segs, 0);

  bsGetKey (obj, _Align_homol, &key);
  while (key != homol->key && bsGetKey (obj, _bsDown, &key));
  if (key != homol->key) messcrash ("Could not find key to %s", name (homol->key));
  bsPushObj (obj);

  bsGetKey (obj, _bsRight, &key);
  while (key != homol->meth && bsGetKey (obj, _bsDown, &key));
  if (key != homol->meth) messcrash ("Could not find key to %s", name (homol->meth));
  bsPushObj (obj);

  bsGetData (obj, _bsRight, _Float, &searchFloat);
  while (searchFloat != homol->score && bsGetData (obj, _bsDown, _Float, &searchFloat));
  if (searchFloat != homol->score) messcrash ("Could not find score %f", homol->score);
  bsPushObj (obj);

  bsPushIntData (obj, homol->qstart);
  bsPushIntData (obj, homol->qend);
  bsPushIntData (obj, homol->sstart);
  bsPushIntData (obj, homol->send);

  /* Cannot use bsFlatten or bsGetArray, which uses bsFlatten, since I then have to give a
     maximum number of segs.  Instead iterate till all are found */

  bsFindTag (obj, _Segs);
  bsPushObj (obj);

  while (bsGetData (obj, _bsRight, _Int, &searchInt)) {
    pushText (BelvuMatchStack, messprintf ("%d ", searchInt));
    bsPushObj (obj);
  }

  alignDisplay (homol->key, 0, 0);
}

/****************************************************************************************************************/
static void pepMapShowBlixemPext (int box) 
{
#if !defined (NO_POPEN)
  HOMOL *homol ;
  Array homols;
  MAPCONTROL currentMap = currentMapControl ();
  PEPLOOK *look;
  METHOD *meth ;
  int i;
  FILE *pipe ;
  char *cp, *script, *qseq ;
  char blixem_script[] = "wscripts/blixem_script" ;
  COLCONTROL control = currentMap->control ;
  COLINSTANCE instance = array (control->boxIndex,box,COLINSTANCE);

  homols = * (instance->proto->convertResults);

  /* Open pipe to Blixem */
  /* Can't call callScriptPipe since we want write access and buildCommand 
     is not an exported function in the graph library */
  if ((cp = filName (blixem_script, 0, "x")))
    script = cp ;
  else {
    /* script = blixem_script ;*/
    messout ("Could not find an executable %s", blixem_script);
    return;
  }
  printf ("Calling \"%s\"\n", script);
  fflush (stdout);
  pipe = (FILE *)popen (script, "w");

  look = (PEPLOOK*) currentMap->look;
  qseq = messalloc (arrayMax (look->pep)+1);
  for (cp = qseq, i = 0; i < arrayMax (look->pep);)
      *cp++ = pepDecodeChar[ (int)arr (look->pep, i++, char)] ;

  fprintf (pipe, "%s\n%s\n", currentMap->name, qseq);
  fprintf (pipe, "# exblx\n# BLASTP\n");
  for (i = 0 ; i < arrayMax (homols) ; i++)
  { 
    homol = arrp (homols, i, HOMOL) ;
    meth = method (0, homol->meth) ;
    if (meth->flags & METHOD_BLIXEM_P)
      fprintf (pipe, "%.0f (+1) %d %d %d %d %s\n",
	       (meth->flags & METHOD_BELVU ? -3 : homol->score),  
	       homol->qstart, 
	       homol->qend,   
	       homol->sstart, 
	       homol->send,   
	       (meth->flags & METHOD_BELVU ? "GREEN" : name (homol->key)));
  }

  fflush (stdout);
  fprintf (pipe, "%c\n", EOF) ; /* To close the pipe, sigh */
  fflush (pipe);
  /* pclose is no good - waits till Blixem is finished. */

  messfree (qseq);
#endif  /* not NO_POPEN */
}

/****************************************************************************************************************/
#include "blxview.h"
static void pepMapShowBlixemPint (int box) 
{
  HOMOL *homol ;
  Array homols, a;
  MAPCONTROL currentMap = currentMapControl ();
  PEPLOOK *look;
  METHOD *meth ;
  MSP MSPlist, *msp, *msp2;
  KEY key; 
  COLCONTROL control = currentMap->control ; 
  COLINSTANCE instance = array (control->boxIndex,box,COLINSTANCE);

  int i;
  char *cp, *qseq, opts[] = "P";

  homols = * (instance->proto->convertResults);

  look = (PEPLOOK*) currentMap->look;
  qseq = messalloc (arrayMax (look->pep)+1);
  for (cp = qseq, i = 0; i < arrayMax (look->pep);)
      *cp++ = pepDecodeChar[ (int)arr (look->pep, i++, char)] ;

  MSPlist.next = 0;
  msp = &MSPlist ; 
  for (i = 0 ; i < arrayMax (homols) ; i++)
  { 
    homol = arrp (homols, i, HOMOL) ;
    meth = method (0, homol->meth) ;

    if (meth->flags & METHOD_BLIXEM_P) {
      msp->next = (MSP*)messalloc (sizeof (MSP));
      msp = msp->next;

      msp->score = homol->score;
      msp->qstart = homol->qstart;
      msp->qend = homol->qend;
      msp->sstart = homol->sstart;
      msp->send = homol->send;
      strcpy (msp->frame, " (+1)");
      if (meth->flags & METHOD_BELVU) {
	  msp->score = -3;
	  msp->id = GREEN;
      }
      else strncpy (msp->sname, name (homol->key), FULLNAMESIZE);
    }
  }

  /* Get sequences */
  for (msp = MSPlist.next; msp ; msp = msp->next) {
    if (msp->sseq) continue;
      
    if (lexword2key (msp->sname, &key, _VPeptide) && (a = arrayGet (key, char, "c"))) { 
      pepDecodeArray (a) ;
      msp->sseq = arrp (a, 0, char) ;
      a->base = 0 ;
      arrayDestroy (a) ;

      /* Avoid duplication of sequence */
      for (msp2 = msp->next; msp2 ; msp2 = msp2->next)
	  if (!strcmp (msp->sname, msp2->sname)) 
	      msp2->sseq = msp->sseq;
    }
  }

  blxview (qseq, currentMap->name, 1, 0, MSPlist.next, opts) ;

}

static MENUOPT belvuMenu[] = {
  { MENU_FUNC_CAST pepMapShowText, "Show as Text" },
  { MENU_FUNC_CAST pepMapShowBelvu, "Analyse in Belvu" },
  { MENU_FUNC_CAST pepMapPfamWWW, "Pfam WWW page" },
  { MENU_FUNC_CAST pepMapShowBlixemPint, "Analyse in Blixem, internal" },
  { MENU_FUNC_CAST pepMapShowBlixemPext, "Analyse in Blixem, external" },
  { 0, 0 }
} ;

static MENUOPT blixemPMenu[] = {
  { MENU_FUNC_CAST pepMapShowText, "Show as Text"},
  { MENU_FUNC_CAST pepMapShowBlixemPint, "Analyse in Blixem, internal"},
  { MENU_FUNC_CAST pepMapShowBlixemPext, "Analyse in Blixem, external"},
  { 0, 0 }
} ;


static int mystrcasestr (char *s1, char *s2) /* mieg: strcasestr is a different standard call */
{
  char *cp;
  int s1len = strlen (s1), 
  s2len = strlen (s2);

  for (cp = s1; cp <= s1+s1len-s2len; cp++)
    if (!strncasecmp (cp, s2, s2len)) return 1;

  return 0;
  /* Does this function really not exist ? There must be a better way */
}

/****************************************************************************************************************/
static void homolDraw (COLINSTANCE instance, float *offset) 
{ 
  float y1, y2, fmax = 0;
  PEPLOOK *look = (PEPLOOK*) instance->map->look ;
  COLCONTROL control = instance->map->control ;
  int i, box ,xoff;
  BUMP bump = 0 ;
  HOMOL *homol ;
  Array homols;
  METHOD *meth ;
  HOMOLPRIV private = instance->private;
  
  currentMapControl ();
  
  homols = * (instance->proto->convertResults);
       
  if (private->bump)
    bump = bumpCreate (1000, 0) ;
  for (i = 0 ; i < arrayMax (homols) ; i++)
    { 
      homol = arrp (homols, i, HOMOL) ;

      /* Make sure it's the right method for this column */
      if (!mystrcasestr (instance->name, name (homol->meth))) continue;

      meth = method (0, homol->meth) ;
	  
      /* checks to prevent arithmetic crash */
    
      if (meth->flags & METHOD_SCORE_BY_OFFSET && !meth->min)
	meth->min = 1 ;
      if (meth->flags & METHOD_SCORE_BY_WIDTH && meth->max == meth->min)
	meth->max = meth->min + 1 ;

      y1 = MAP2GRAPH (instance->map, homol->qstart);
      y2 = MAP2GRAPH (instance->map, homol->qend);

      if (y2 > control->graphHeight -0.5){
	if (y1 > control->graphHeight -0.5)
	  continue;
	else
	  y2 = control->graphHeight -0.5;
      }
      if (y1 < control->topMargin+0.5){
	if (y2 < control->topMargin + 0.5)
	  continue; 
	else
	  y1 = control->topMargin+0.5;
      }
      box = graphBoxStart () ;
      if (meth->flags & METHOD_SCORE_BY_OFFSET)
	{ float dx = 6 - 4*log10 ((double) (homol->score/meth->min)) ;
	  if (dx < 0)
	    dx = 0 ;
	  graphRectangle (*offset + dx, y1, *offset + dx + .85, y2) ;
	  if (fmax < dx)
	    fmax = dx ;
	}
      else if (meth->flags & METHOD_SCORE_BY_WIDTH)
	{ float dx = 0.5 + 0.5 * (homol->score - meth->min) / 
	    				 (meth->max - meth->min) ;
	  if (dx < 0.5) dx = 0.5 ;
	  if (dx > 1) dx = 1 ;
          xoff = 1 ;
          if (bump)
	    {
	      if (look->map->mag < 0)
		{
		  y1 = -y1 ; y2 = -y2 ;
		  bumpItem (bump, 2, (y2-y1), &xoff, &y1) ;
		  y1 = -y1 ; y2 = -y2 ;
		}
	      else
		bumpItem (bump, 2, (y2-y1), &xoff, &y1) ;
	    }
          graphRectangle (*offset + (xoff - dx), y1, 
                          *offset + (xoff + dx), y2) ;
          if (fmax < xoff+dx) fmax = xoff+dx ;

	}
      else
	{ graphRectangle (*offset + 0.25, y1, *offset + 1.75, y2) ;
	  if (!fmax) fmax = 1;
	}
      graphBoxEnd () ;

      graphBoxDraw (box, BLACK, meth->col) ;
      controlRegBox (instance, box, homol) ;
      if (meth->flags & METHOD_BELVU)
	graphBoxMenu (box, belvuMenu) ;
      else if (meth->flags & METHOD_BLIXEM_P)
	graphBoxMenu (box, blixemPMenu) ;

    }
 
  if (fmax) *offset += fmax + 2 ;

  if (bump)
    bumpDestroy (bump);
}

/**************************************************************************************************/

static void homolPrivDestroy (void * p)
{
 HOMOLPRIV private = (HOMOLPRIV)p;

 if (arrayExists (private->segs)) arrayDestroy (private->segs);
}

/**************************************************************************************************/

static BOOL homolUnSetSelect (COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;

  control->activeKey = 0;
  look->startcol = 0;
  look->endcol = 0;
  
  *look->messageText = 0;
  
  if (look->messageBox)
    graphBoxDraw (look->messageBox,-1,-1);

  return FALSE;
}

static BOOL homolSetSelect (COLINSTANCE instance, int box, double x, double y)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  HOMOLPRIV private = instance->private;
  PEPLOOK *look = map->look;
  HOMOL *homol = (HOMOL *) arr (control->boxIndex2,box,void *);

  control->activeKey = homol->key;
  control->activeBox = box;
  look->startcol = homol->qstart;
  look->endcol = homol->qend;

  graphBoxDim (box,&private->x,0,0,0);
  return FALSE;

}

static void homolDoColour (COLINSTANCE instance, int box)
{
  MAPCONTROL map = instance->map;
  HOMOLPRIV private = instance->private;
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;
  HOMOL *homol = (HOMOL *) arr (control->boxIndex2,box,void *);
  METHOD *meth ;
  float x1;

  meth = method (0, homol->meth) ;
  graphBoxDim (box,&x1,0,0,0);
  if (homol->key == control->activeKey && instance == control->activeInstance){
    if ( x1 == private->x && 
     look->startcol == homol->qstart && look->endcol == homol->qend){
      graphBoxDraw (box, BLACK,activecol);
      strncpy (look->messageText, name (homol->key),100);
      strcat (look->messageText, messprintf (" %d  %d ",homol->sstart,homol->send));
      strcat (look->messageText, name (homol->meth));
      strcat (look->messageText,messprintf (" %4.2f ",homol->score));
      strcat (look->messageText, messprintf (" (%d - %d) ",homol->qstart,homol->qend));
      if (look->messageBox)
	graphBoxDraw (look->messageBox,-1,-1);
      control->activeBox = box;
    }
    else
      graphBoxDraw (box, BLACK,keyfriendcol);
  }
  else if (homol->key == control->activeKey )
    graphBoxDraw (box, BLACK,friendcol);
  else{
    if (look->startcol <= homol->qstart && look->endcol >= homol->qend)
      graphBoxDraw (box, BLACK, friendcol);
    else if ((look->endcol == 0)&& (look->startcol >= homol->qstart && look->startcol <= homol->qend))
      graphBoxDraw (box, BLACK, friendcol);
    else
      graphBoxDraw (box, BLACK, meth->col);
  }
}
/*************************************************************************************************/
static void homolFollowBox (COLINSTANCE instance, int box,double x, double y)
{
  MAPCONTROL map = instance->map;
  COLCONTROL control = map->control;
  HOMOL *homol = (HOMOL *) arr (control->boxIndex2,box,void *);
  PEPLOOK *look = map->look;

  if (look->block){
    look->activeStart = homol->qstart;
    look->activeEnd = homol->qend;
  }
  display (homol->key, map->key, 0);
}
/*************************************************************************************************/
static void pepHomolSave (COLINSTANCE instance, OBJ init) 
{
  HOMOLPRIV private = instance->private;
  
  if (private->bump)
    bsAddTag (init,str2tag ("HOM_bump"));
}

/*************************************************************************************************/
extern  BOOL homolCreate (COLINSTANCE instance, OBJ init) 
{
  HOMOLPRIV private;
  
  instance->draw = homolDraw ;
  instance->configure = homolConfigure;
  instance->doColour = homolDoColour; 
  instance->unSelectBox = homolUnSetSelect; 
  instance->setSelectBox =  homolSetSelect;
  instance->followBox = homolFollowBox;
  instance->save = pepHomolSave;
  instance->configFinal = pepConfigFinal;

  private = (HOMOLPRIV) handleAlloc (homolPrivDestroy, instance->handle , sizeof (struct homolColPriv));
  private->bump = FALSE;
  instance->private = private;

  if (init){
    if (bsFindTag (init,str2tag ("HOM_bump")))
      private->bump = TRUE;
  }

  return TRUE ;
}

/*********************************** END **********************************************************/
/* put the names colomn stuff here to share the homol data */

static void homolNameDraw (COLINSTANCE instance, float *offset) 
{ 
  float y1, y2, fmax = 0, namemax = 0 ,y3;
  PEPLOOK *look = (PEPLOOK*) instance->map->look ;
  COLCONTROL control = instance->map->control ;
  int i, box ,xoff;
  BUMP bump = 0 ;
  HOMOL *homol ;
  Array homols;
  HOMOLPRIV private = instance->private;
  
  currentMapControl ();
  
  homols = * (instance->proto->convertResults);
 
  if (private->bump)
    bump = bumpCreate (private->width, 0) ;
  for (i = 0 ; i < arrayMax (homols) ; i++)
    { 
      homol = arrp (homols, i, HOMOL) ;
      
      if (!mystrcasestr (instance->name, name (homol->meth)))
	continue;
      method (0, homol->meth) ;
      
      y1 = MAP2GRAPH (instance->map, homol->qstart);
      y2 = MAP2GRAPH (instance->map, homol->qend);
      
      if (y2 > control->graphHeight -0.5){
	if (y1 > control->graphHeight -0.5)
	  continue;
	else
	  y2 = control->graphHeight -0.5;
      }
      if (y1 < control->topMargin+0.5){
	if (y2 < control->topMargin + 0.5)
	  continue; 
	else
	  y1 = control->topMargin+0.5;
      }
      /* Show name of homologous object - either here or in text-area to the right */
      xoff = 0 ;
      fmax = 0;
      y3 = (y1+y2)/2;
      if (bump)
	{
	  if (look->map->mag < 0)
	    {
	      y1 = -y1 ; y2 = -y2 ;
	      bumpItem (bump,strlen (name (homol->key))+1,1, &xoff, &y3) ;
	      y1 = -y1 ; y2 = -y2 ;
	    }
	  else
	    bumpItem (bump,strlen (name (homol->key))+1,1, &xoff, &y3) ; 
	}
      box = graphBoxStart ();
      graphText (name (homol->key), *offset + xoff, y3);
      if (strlen (name (homol->key))+2 > namemax)
	namemax = strlen (name (homol->key))+2; 
      if (fmax < xoff + namemax)
	fmax = xoff +namemax;
      graphBoxEnd ();
      controlRegBox (instance, box, homol) ;
    }
  *offset += fmax+4;
  if (private->bump)
    bumpDestroy (bump);
  
}

static void pepHomolNameSave (COLINSTANCE instance, OBJ init) 
{
  HOMOLPRIV private = instance->private;
  
  if (private->bump)
    bsAddTag (init,str2tag ("HOM_NAME_bump"));
  if (private->width)
    bsAddData (init,str2tag ("HOM_NAME_width"),_Int,&private->width);
}

struct configLocalsName{
  int width;
  BOOL bump;
};

static BOOL homolNameConfigure (COLINSTANCE instance)
     
{
  HOMOLPRIV private = instance->private;
  struct configLocalsName *cf = (struct configLocalsName *) messalloc (sizeof (struct configLocalsName));
  float line = 2.0;

  if (controlCreateConfig (instance,cf,"Configure pepdisp name column",0.5,0.15)){
    /* initialise data */
    cf->width = private->width;
    cf->bump = private->bump;

    graphIntEditor ("HomolName Display Width :",&cf->width,4.0,line++,0);
    graphToggleEditor ("Bump",&cf->bump,4.0,line++);
  
    graphRedraw ();
  }

  return FALSE;
}

static void pepNameConfigFinal (COLINSTANCE instance, void *locals, BOOL ok){
  struct configLocalsName *cf = locals;
  HOMOLPRIV private = instance->private;

  if (ok){
    private->bump = cf->bump;
    private->width = cf->width;
  }
  else
    messfree (cf);
}

extern  BOOL homolNameCreate (COLINSTANCE instance, OBJ init) 
{
  HOMOLPRIV private;
  int x;
  
  instance->draw = homolNameDraw ;
  instance->configure = homolNameConfigure;
  instance->doColour = homolDoColour; 
  instance->unSelectBox = homolUnSetSelect; 
  instance->setSelectBox =  homolSetSelect;
  instance->followBox = homolFollowBox;
  instance->save = pepHomolNameSave;
  instance->configFinal = pepNameConfigFinal;

  private = (HOMOLPRIV) handleAlloc (homolPrivDestroy, instance->handle , sizeof (struct homolColPriv));
  private->bump = FALSE;
  private->width = 30;
  instance->private = private;

  if (init){
    if (bsFindTag (init,str2tag ("HOM_NAME_bump")))
      private->bump = TRUE;
    if (bsGetData (init,str2tag ("HOM_NAME_width"),_Int,&x))
      private->width = x;
  }
  
  return TRUE ;
}




 
 
 
 
 
 
 
 
 
 
 
 
