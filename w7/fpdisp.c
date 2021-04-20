/*  File: fpdisp.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 17 03:18 1998 (rd)
 * * Jul 23 15:36 1998 (edgrif): Add fmap.h public header for fmap
 *      function declares.
 * * Jun  3 11:41 1996 (rd)
 * Created: Thu Jan  6 11:12:36 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: fpdisp.c,v 1.5 2016/11/23 19:15:12 mieg Exp $ */

#include "acedb.h"

#include "display.h"
#include "lex.h"
#include "sysclass.h"
#include "classes.h"
#include "tags.h"
#include "systags.h"
#include "bs.h"
#include "plot.h"
#include "session.h"
#include "fingerp.h"
#include "pick.h"
#include "dna.h"
#include "restriction.h"
#include "fmap.h"

/************************************************************/

struct LookStruct {
  /* Look not used, we just complete the type with nothing,
     so we can use the type LookStruct to pass it between functions */
  int emptyStruct;
} ;

typedef struct GelStruct {
  KEY clone, clef ;
  Array fp, bp ;
  float maxFp, maxBp ;
} GEL ;

static BOOL isBp = TRUE ;
static float maxBp, maxFp ;
static Array bands = 0 ;
static float rounding = 1 ;
static Array fpSelectedBands = 0 , fpSelectedClefs = 0 , fpAllClefs = 0 , fpAllClones = 0;

static Graph fpGraph = 0 ;
static MAP fpMap = 0 ;
static Array fpGels = 0 ;
static Array fpBoxIndex = 0 ;

static int fpActiveBox = 0, minLiveBox = 0 , minClefBox = 0 ;
static int fpActiveGel = 0 ;
static char fpRoundingBuffer[25] ;

/************************************************************/

static BOOL  fingerPrintData (Array *fpp, Array* bpp,
			      KEY clone, KEY *clef) ;
static BOOL  fpRead (void) ;	/* returns FALSE if unsuccessful */
static Array fpLength2bp(Array fp) ;
static Array fpBp2Length(Array bp) ;

static void  cptAllPredictedBands (void);

static void fpDraw(void) ;
static void fpPick (int box,  double x , double y) ;

/************************************************************/

static void fpHisto(void)
{
  Array histo ;
  int i ;
  
  if (!fpRead())
    return;

  if(!bands) 
    return ;
  
  histo  = arrayCreate(5000,int);
  i = arrayMax(bands);
  while(i--)
    array(histo,arr(bands,i,short),int)++;

  plotHisto ("Bands distribution", histo);

  return;
} /* fpHisto */
 
/******************************************/
 
static BOOL fpRead(void)
{
  static char fileName[FIL_BUFFER_SIZE] ;
  static char dirName[DIR_BUFFER_SIZE] ;
  int i = 0 , n ;
  void *vp = &n;
  unsigned int nb = 20000, nr, size = sizeof(short);
  FILE *fpfile, *corfile;
   
  if (arrayExists(bands))
    return TRUE;

  if (!dirName[0])
    strcpy(dirName, sessionFilName(0,0,0)) ;

  if (!(corfile = filqueryopen (dirName, fileName, "cor", "r", 
				"Finger-print rawdata file")))
    return FALSE;
  filclose (corfile) ;
  if (strlen(fileName) > 4 && fileName[strlen(fileName) - 4] == '.') 
    fileName[strlen(fileName) - 4] = '\0' ;
  bands = arrayCreate(nb,short);
  if ((fpfile = filopen (messprintf("%s/%s", dirName, fileName),"bnd","rb")))
    { do
	{ array(bands,(++i)*nb,short) = 0; /* to create space */
	  vp = arrp(bands,nb*(i-1),short); /* possible relocation */
	  fprintf(stderr,"I read %d bands\n",i*nb);
	}
      while ((nr=fread(vp,size,nb,fpfile)) == nb) ;
      
      arrayMax(bands) -= nb - nr; /* artificial space removed */
    }
  else
    { messout ("I reformat the ascii file %s.cor into a binary file %s.bnd",
	       fileName, fileName);
      if (!(corfile= filopen(messprintf("%s/%s", dirName, fileName),"cor","r")) ||
	  !(fpfile= filopen(messprintf("%s/%s", dirName, fileName),"bnd","wb")))
	return FALSE ;
      
      array(bands,(++i)*nb,short) = 0; /* to create space */
      freespecial ("\n\t\"\\") ;
      while (freeread(corfile))
	if(freeint(&n))
	  array(bands,i++,short) = n ;
	else
	  { messout("Error parsing line %d of cor file : %s",
		    i, freeword());
	    break;
	  }
      
      filclose(corfile) ;
      fwrite(arrp(bands,0,short),size,arrayMax(bands),fpfile);
    }
  
  messout("I read %d gel bands", arrayMax(bands)) ;
  fpHisto() ;
  filclose(fpfile);
  return TRUE;
} /* fpRead */

/*********************************************************************/
/********************  public routines   *****************************/
/*********************************************************************/
  /* set *clefp or start with *clefp == 0, and repeat to get them all */
/***********************************************************/

static BOOL fingerPrintData (Array *fpp, Array* bpp, KEY clone, KEY *clefp) 
{
  static int mm = 0 ;
  int from, i, n ;
  Array  fp = 0 ;
  float x , *xp ;
  KEY direc ;
  OBJ Clone ;
  BOOL isGotLength = FALSE ;

  if (!(Clone = bsCreate(clone)))
    return 0 ;

  if (*clefp)
    { if (!bsFindKey (Clone, _Gel, *clefp) || !bsGetKeyTags (Clone, _bsDown , clefp))
	return 0 ;
    }
  else
    if (!bsGetKeyTags (Clone, _Gel, clefp)) /* generic format */
      {
	if (!fp && class(clone) == _VClone &&  /* fall back on nematode specific format */
	    bsGetData(Clone,_Bands,_Int,&from) &&
	    bsGetData(Clone,_bsRight,_Int,&n))
	  { 
	    if(!bands)
	      {
		if (!fpRead())
		  return 0;
	      }
	    
	    isGotLength = TRUE ;
	    if (bands)
	      {
		if(from+n <= arrayMax(bands))
		  { 
		    KEY fpKey ;
		    
		    lexaddkey ("Observed", &fpKey, _VMotif) ;
		    *clefp = fpKey ;
		    fp = arrayCreate(n, float) ;
		    while(n--)
		      array(fp,n,float) = (float) (array(bands,from+n,short)) ;
		    if (arrayMax(fp)) /* && isWriteAccess())*/
		      { bsDestroy(Clone) ;
			if ((Clone = bsUpdate(clone)))
			  { bsAddTag (Clone, _Gel) ;
			    bsAddKey (Clone, _bsRight, fpKey) ;
			    bsPushObj (Clone) ;
			    bsAddTag (Clone, _Band_Lengths) ;
			    /* first value allready stored */
			    for (n = 0 ; n < arrayMax(fp) ; n++)
			      { x = arr(fp, n, float) ;
				bsAddData(Clone, _bsRight, _Float, &x) ;
			      }
			    bsSave(Clone) ;
			    Clone = 0 ;
			  }		
		      }
		  }
		else
		  messout
		    ("Clone %s\n Band %d > Max %d of data file",
		     name(clone), from+n, arrayMax(bands)) ;
	      }
	  }
	goto laba ;
      }

  if (bsPushObj (Clone))
      { direc = 0 ;
	i = 0 ;
	if ( bsGetData (Clone, _Bands, _Float, &x))
	  direc = _bsRight ;
	else if ( bsGetData (Clone, _Band_Lengths, _Float, &x))
	  { direc = _bsRight ;
	    isGotLength = TRUE ;
	  }
	else if ( bsGetData (Clone, _Band, _Float, &x))
	  direc = _bsDown ;
	fp = arrayCreate (8, float) ;
	if (direc)
	  do { array (fp, i++, float) = x ;
	     } while (bsGetData (Clone, direc, _Float, &x)) ;
	
      }
  
  
 laba:
  if (!fp && mm++ < 3)
    messout("Sorry, no finger print data known about clone %s",
	    name(clone)) ;

  if (Clone)
    bsDestroy(Clone) ;
  
  if (rounding == 0 || rounding <= 0.001) rounding = .001 ;
  if (fp)
    for (n = 0 , xp = arrp(fp, 0, float); n < arrayMax(fp) ; n++, xp++)
      { i = *xp / rounding ; 
	*xp = i * rounding ;
      }
  if (!fp || !arrayMax(fp))
    { arrayDestroy (fp) ;
      return FALSE ;
    }
  else
    arraySort (fp, floatOrder) ;
    
  if (isGotLength)
    { isBp = FALSE ;
      *fpp = fp ;
      *bpp = fpLength2bp(fp) ;
    }
  else
    { *bpp = fp ;
      *fpp = fpBp2Length (fp) ;
    }
  return TRUE ;
} /* fingerPrintData */

/*************************************************************/
/*************************************************************/

  /* donnees is a boolean matrix clones/lengths */
  /* donnees must be created as 2 arrays of KEYSET */
void pmapCptGetData (Array marInDef, Array defInMar, KEYSET clones, KEYSET bands)
{ KEY key = 0 , kb, clef ;
  Array aa,  fp, bp, xp ;
  Associator assB ;
  void *vp, *dummy ;
  int  i, j, nB = 0 , nC, nBmax = 0 , nCmax ;
  KEYSET *ksp ;
  OBJ obj ;

  if (!keySetActive(&aa, &dummy))
    { messout("First select a keyset containing clones") ;
      return ;
    }
  j = 0 ;
  for (i = 0 ; i < keySetMax(aa) ; i++)
    { key = keySet(aa, i) ;
      if (class(key) == _VClone && (obj = bsCreate(key)))
	{ if (bsFindTag(obj, _Bands))
	    keySet(clones, j++) = key ;
	  bsDestroy(obj) ;
	}
    }
  nCmax = keySetMax(clones) ;
  if (!nCmax)
    { messout("First2 select a keyset containing clones") ;
      return ;
    }
  assB = assCreate() ;
  i = nCmax ;
  ksp = arrayp(marInDef, i - 1, KEYSET) ;
  while(i--)
    *ksp-- = keySetCreate() ;
  for (nC = 0 ; nC < nCmax ; nC++) 
    { clef = 0 ;
      while (fingerPrintData(&fp, &bp, keySet(clones, nC), &clef))
	{ if (keySetExists(fpSelectedClefs) &&
	      ! keySetFind (fpSelectedClefs, clef, &j))
	    continue ;
	  xp = isBp ? bp : fp ;
	  j = arrayMax(xp) ;
	  while (j--)
	    { lexaddkey(messprintf("%s.%g", name(clef), arr(xp, j, float)), 
			&kb, _VCalcul) ;
	      if (!assFind (assB, assVoid(kb), &vp))
		{ vp = assVoid(nBmax) ;
		  assInsert(assB, assVoid(kb), vp) ;
		  nB = nBmax ;
		  array(defInMar, nB, KEYSET) = keySetCreate() ;
		  keySet(bands, nBmax++) = kb ;
		}
	      else
		nB = assInt(vp) ;
	      ksp = arrp(marInDef, nC, KEYSET) ;
	      keySet(*ksp, keySetMax(*ksp)) = nB ;
	      ksp = arrp(defInMar, nB, KEYSET) ;
	      keySet(*ksp, keySetMax(*ksp)) = nC ;
	    }
	  arrayDestroy(fp) ;
	  arrayDestroy(bp) ;
	}
    }
  assDestroy(assB) ;
  i = arrayMax(marInDef) ;
  ksp = arrp(marInDef, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    keySetSort(*ksp) ;
  i = arrayMax(defInMar) ;
  ksp = arrp(defInMar, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    keySetSort(*ksp) ;

  return;
} /* pmapCptGetData */

/***************************************************************/
/***************************************************************/


static void fpDestroy(void)
{
  int i ;

  if (graphExists(fpGraph))
    {
      fpGraph = 0 ;
      i = arrayMax (fpGels) ;
      while(i--)
	{ arrayDestroy(arrp(fpGels, i, GEL)->fp) ;
	  arrayDestroy(arrp(fpGels, i, GEL)->bp) ;
	}
      arrayDestroy(fpGels) ;
      arrayDestroy(fpBoxIndex) ;
      arrayDestroy(fpSelectedBands) ;
      arrayDestroy(fpSelectedClefs) ;
      arrayDestroy(fpAllClefs) ;
      arrayDestroy(fpAllClones) ;

      mapDestroy(fpMap) ;
      graphAssRemove (&MAP2LOOK_ASSOC) ;
    }

  return;
} /* fpDestroy */

/**************************************************************/

static void fpRmClone(void)
{ int n , i , j ;
  KEY clone ;

  if (!graphActivate(fpGraph))
    return ;
  n = fpActiveGel ;
  if (n == 100000)
    return ;
  clone = arr (fpGels, n, GEL).clone ;
  i = arrayMax(fpGels) ; j = 0 ;
  while (i--)
    if (clone == arr (fpGels, i, GEL).clone)
     { arrayDestroy(arr(fpGels, i, GEL).fp) ;
       arrayDestroy(arr(fpGels, i, GEL).bp) ;
       for (j = i ; j < keySetMax(fpGels) - 1 ; j++)
	 array(fpGels, j, GEL) = array(fpGels, j + 1, GEL) ; 
       arrayMax(fpGels)-- ;
     }
  
  fpDraw() ;

  return;
} /* fpRmClone */

static void fpUnBlock (void)
{
  displayUnBlock() ;
  if (graphActivate(fpGraph))
    graphRegister(PICK, fpPick) ;

  return;
} /* fpUnBlock */

/*************************************/

static BOOL fpInterpret (int *np)
{ int i, j, dummy, n = *np ;
  KEY clef ;

  n /= 2 ;
  
  for (i = 0 , j = 0 ; i < arrayMax(fpGels) ; i++)
    { clef = arr(fpGels, i, GEL).clef ;
      if (keySetFind (fpSelectedClefs, clef, &dummy))
	{ if (j == n)
	    break ;
	  j++ ;
	}
    }
  if (i >= arrayMax(fpGels))
    return FALSE ;
  *np = i ;
  return TRUE ;
} /* fpInterpret */

/*************************************/
/*
static void fpSort(void)
{ int n = arrayMax(fpGels) ;
  KEYSET ks = keySetCreate() ;

  while (n--)
    keySet(ks, n) = arrp(fpGels, n, GEL)->clone ;
  
  messout("No direct coding yet, use mapping from main menu, sorry") ;
}
*/
/*************************************/

static void fpSelectAND (Array x, Array y)
{
  float *xp, *yp, *zp ;
  register int i = arrayMax(x), j = arrayMax(y) , k = 0 ;

  xp = zp = arrp (x, 0, float) ;
  yp = arrp (y, 0, float) ;

  while(i && j)
    { 
                      /*success, skip further in the 3 index */
      if(*xp == *yp)
	{ *zp = *xp ;
	  i-- ; j-- ;
	  xp++ ; yp++ ; zp++ ;
	  k++ ;
	}
      else      /*recall that every index is in increasing order*/
	if (*xp < *yp)
	  { i--, xp++ ; } 
	else
	  { j--; yp++ ;}
    }
  arrayMax (x) = k ;

  return;
} /* fpSelectAND */

/*************************************/

static void fpDoSelectBands (int box)
{
  int n ;
  Array fp ;
  Array old ;

  old = fpSelectedBands ;
  fpSelectedBands = 0 ;
  if (old && !arrayMax(old))
    arrayDestroy(old) ;
  
  if (box > minLiveBox)
    { n = box - minLiveBox - 1 ;
      if (!fpInterpret (&n))
	return ;
      fpActiveGel = n ;
      fp = isBp ? array(fpGels, n, GEL).bp 
	: array(fpGels, n, GEL).fp ;
	 
      if (old)
	fpSelectAND (old, fp) ;
      else
	old = arrayCopy (fp) ;

      fpSelectedBands  = old ;
      fpDraw() ;
    }
  
  return;
} /* fpDoSelectBands */

/**************************************************************/

static void fpSelectBands (void)
{
  int  n = 0 , max ;
   
  n = fpActiveGel ;
  max = arrayMax(fpGels) ;
  if (n >= max)  
    return ;
  
  graphRegister (MESSAGE_DESTROY, fpUnBlock) ;
  graphRegister (PICK, fpDoSelectBands) ;
  displayBlock ((BlockFunc)fpUnBlock,
		"Pick in succession clones or lanes, their intersect will turn thick. "
		"This will continue until you remove this message or click outside the gel.") ;
} /* fpSelectBands */

/*************************************/

static void fpPickSwitch (int box)
{
  int n, a, b;
  GEL gel ;

  if (box > minLiveBox)
    { n = box - minLiveBox - 1 ;
      if (!fpInterpret (&n))
	return ;

      if (fpActiveGel == 100000)
	{ fpActiveGel = n ;
	  fpDraw() ;
	}
      else
	{ a = fpActiveGel ; b = n ;
	  if (a == b)
	    return ;
	  gel = array(fpGels, a, GEL) ;
	  array(fpGels, a, GEL) = array(fpGels, b, GEL) ;
	  array(fpGels, b, GEL) = gel ;
	
	  fpActiveGel = 100000 ;
          fpDraw() ;
	}
    }

  return;
} /* fpPickSwitch */

static void fpSwitch (void)
{ int  n = 0 , max ;
   
  n = fpActiveGel ;
  max = arrayMax(fpGels) ;

  if (n < 0 || n >= max)  
    return ;
  if (n == max - 1)
    {
      if (n>0) n-- ; 
      else return ;
    }

  if (max >= 2)
    { graphRegister (MESSAGE_DESTROY, fpUnBlock) ;
      graphRegister (PICK, fpPickSwitch) ;
      displayBlock ((BlockFunc)fpUnBlock,
		    "Pick in succession a lane and a target."
		    "This will continue until you remove this message or click outside the gel.") ;
      fpActiveGel = 100000 ; /* to start parity in switch routine */
      return ;
    }
/* difficille a cause du clef select 
  gel = array(fpGels, n , GEL) ;
  array(fpGels, n , GEL) =  array(fpGels, n + 1, GEL) ;
  array(fpGels, n + 1, GEL) = gel ;
  fpDraw() ;
*/  

  return;
} /* fpSwitch */

/*************************************/

static void fpDoAddClone(KEY clone)
{ 
  fpDisplay (clone) ;
  displayRepeatBlock () ;
}

static void fpAddClone(void)
{ 
  if (!graphActivate(fpGraph))
    return ;
  graphRegister (MESSAGE_DESTROY, fpUnBlock) ;
  displayBlock (fpDoAddClone,
	"Double-clicking a clone will show its finger-printing. "
	"This will continue until you remove this message") ;

  return;
} /* fpAddClone */

static void fpAddKeySet(void)
{
  KEYSET aa ;
  void *dummy ;
  int i ;

  if (!graphActivate(fpGraph))
    return ;
  if (!keySetActive(&aa, &dummy))
    { messout("First select a keyset containing clones") ;
      return ;
    }

  fpClearDisplay() ;
  for (i=0; i < keySetMax(aa) ; i++)
    fpDisplay(keySet(aa,i)) ;

  return;
} /* fpAddKeySet */

/*************************************/

static Array fpLength2bp(Array fp)
{
  Array bp;
  double x ;
  int i = arrayMax(fp) ;
	   
  bp = arrayCreate(arrayMax(fp), float) ;

  while(i--)
    { x = arr(fp, i, float) ;
      /* excellent Fit mathematica sur 79 data d'etalonage  */
      x = log(x) ;
      x = 133.164 - 54.4586 * x + 7.9092 * x * x - .389621 * x * x * x ;
      x = exp(x) ;
      array(bp, i, float) = x ;
    }
  arraySort (bp, floatOrder) ;
  return bp ;
} /* fpLength2bp */

/*************************************/

static Array fpBp2Length(Array bp)
{ Array fp ;
  double x; float y ;
  int i ;

  i = arrayMax(bp) ;
  fp = arrayCreate(arrayMax(bp), float) ;
  while(i--)
    { x = arr(bp, i, float) ;
      if (x < 1)
	x = 1 ;
	     /* 
	        Log relation between lenght and distance of migration:
		figure 6.1 page 6.5 ; Molecular cloning, a laboratory manual
		second ed. Sambrook, Fritsch, Maniatis
		ColdSping Harbor lab. press, 1989
	
                log10 isdeclared in mystdlib.h 
		the limits 10 to 10000 (so 0 <= y < 11.5)
		are recalled in a graphTextMessage.
		The inverse relation is used in drawScale.
	     */
      y = 3.5 * ( 4 - log10(x)) ;
      array(fp, i, float) = y ;
    }
  arraySort (fp, floatOrder) ;
  return fp ;
} /* fpBp2Length */

/*************************************/

static void fpRounding(char *cp) 
{ int n, i, level ;
  float x ;
  GEL *gelp ;
  KEYSET ks = keySetCreate () ;

  level = freesettext(cp,"") ;
  freecard (level) ;
  if (freefloat(&x) && (x>0) && rounding != x )
    { rounding = x ;
      freeclose(level) ;
    }
  else
    { freeclose(level) ;
      return ;
    }

  sprintf(fpRoundingBuffer, "%g", rounding) ;

  i = arrayMax(fpGels) ;
  for (i = 0 ; i < keySetMax (fpGels) ; i++)
    { gelp = arrp(fpGels, i, GEL) ;
      keySet (ks, i) = gelp->clone ;
    }

  n = arrayMax(fpGels) ;
  while(n--)
    { arrayDestroy(arr(fpGels, n, GEL).fp) ;
      arrayDestroy(arr(fpGels, n, GEL).bp) ;
    }  
  fpGels = arrayReCreate (fpGels, 12, GEL) ;
  fpAllClefs = keySetReCreate (fpAllClefs) ;
  fpAllClones = keySetReCreate (fpAllClones) ;

  for (i = 0 ; i < keySetMax (fpGels) ; i++)
    fpDisplay (keySet(ks,i)) ;
  keySetDestroy (ks) ;
	       
  fpDraw() ;
} /* fpRounding */

/**************************************************************/

static void fpTransform(void)
{
  float mx, x1;
  KEY clef ;
  int i, dummy ;
  Array bp ;

  isBp = !isBp ;
  mx = -1 ;
lao:
  if (!keySetExists(fpSelectedClefs))
    fpSelectedClefs = keySetCopy (fpAllClefs) ;
  keySetSort (fpSelectedClefs) ;
  keySetCompress (fpSelectedClefs) ;

 
  for (i = 0 ; i < arrayMax(fpGels) ; i++)
    { clef = array(fpGels, i, GEL).clef ;
      if (keySetFind (fpSelectedClefs, clef, &dummy))
	{ bp = isBp ? array(fpGels, i, GEL).bp : array(fpGels, i, GEL).fp ;
	  x1 = arr (bp, arrayMax(bp) - 1, float) ;
	  if (x1 > mx) mx = x1 ;
	}
    }
  if (mx == -1)
    { mx = 0 ;
      keySetDestroy (fpSelectedClefs) ;
      goto lao ;
    }

  fpMap->max = mx ;
  fpMap->centre = mx / 2 ;
  fpMap->mag = (mapGraphHeight - topMargin - 5) /
    			(1.05 * mx) ;
  fpDraw() ;

  return;
} /* fpTransform */

/***********************************************************/

static void fpMainLines (LOOK look, float *offset)
     /* used as a MapColDrawFunc for mapInsertCol() */
{
  float y , mx = (isBp ? maxBp : maxFp) ;
  int i , di = mx /10 , j ;

  if (di<=1) di = 2 ;
  j = 1 ;while (j < di) j*= 10 ;
  j /= 10 ; 
  if (j > 1)
    di = di - (di%j) ;
  i = di/j ;
  if (i>5)
    di = 5*j ;

  graphTextHeight (0.75) ;
  for (i = 0 ; i < mx ; i += di)
    { y = MAP2WHOLE(fpMap, i) ;
      graphText (messprintf("%d", i),*offset + 0.5,y-0.25) ;
    }
  graphTextHeight (0) ;
  *offset += 4 ;

  return;
} /* fpMainLines */

/*************************************/

static void fpPick (int box,  double x , double y) 
{
  int n , nn ;
  KEY clone, clef ;

  if (!box)
    return ;
  if (box == fpMap->cursor.pickBox)
    {
      if (!isGifDisplay)
	graphBoxDrag (fpMap->cursor.box, mapCursorDrag) ;
    }
  else if (box == fpMap->thumb.box)
    {  /* coordinates are realtive to thumb.box */
      if (!isGifDisplay)
	graphBoxDrag (fpMap->thumb.box, mapThumbDrag) ;
      else
	{
	  fpMap->centre = WHOLE2MAP (fpMap, y + MAP2WHOLE(fpMap,topMargin + 1)) ;
	  fpDraw () ;
	}
    }
  else if (box > minLiveBox)
    { n = box - minLiveBox - 1 ;
      nn = n & 1 ;
      if (!fpInterpret (&n))
	return ;
      clone = array(fpGels, n, GEL).clone ;
      if (nn)
	display(clone, (KEY)0, 
		!strcmp(name(pickDisplayKey(clone)), DtPmapFingerprint) ? TREE : 0) ;
      if (fpActiveGel != n)
	{ /* graphBoxDraw(minLiveBox + 2 + 2*fpActiveGel,
		       BLACK, WHITE) ;
	  graphBoxDraw(minLiveBox + 2 + 2*fpActiveGel,
		       BLACK, WHITE) ;
	  */
	  fpActiveGel = n ;
	  fpDraw() ;
	}
    }
  else if (box > minClefBox)
    { int dummy  = box - minClefBox - 1 ;

      if (dummy <0 || dummy >= keySetMax(fpAllClefs))
	return ;
      clef = keySet (fpAllClefs, dummy) ;
      
      if (keySetFind (fpSelectedClefs, clef, &dummy))
	keySetRemove (fpSelectedClefs, clef) ;
      else
	keySetInsert (fpSelectedClefs, clef) ;
	
      fpDraw() ;
    }

  return;
} /* fpPick */

/**************************************************************/

static void fpDrawGels(LOOK look, float *offset)
     /* used as a MapColDrawFunc for mapInsertCol() */
{ 
  float *xb, x = *offset, y, dx;
  float oldwidth = graphLinewidth(-1.0) ;
  Array fp, bp ;
  KEY clef ;
  int i, j, box, max, max1, dummy ; 

  max = arrayMax(fpGels) ;
  if (!max)
    return ;
 

  minClefBox= graphBoxStart() ;
  graphText ("Select:", 2 , 2) ;
  graphBoxEnd() ;
  dx = 8 ;
  for (i = 0 ; i < keySetMax (fpAllClefs) ; i++)
    { clef = keySet(fpAllClefs, i) ;
      box = graphBoxStart() ;
      graphText (name(clef), dx , 2) ;
      graphBoxEnd () ;
      
      dx += strlen(name(array(fpGels, i, GEL).clef)) + 1 ;

      if (keySetFind (fpSelectedClefs, clef, &dummy))
	graphBoxDraw (box, BLACK, GREEN) ;
      else
	graphBoxDraw (box, BLACK, WHITE) ;
    }
 
  for (max1 = 0, i = 0 ; i < arrayMax(fpGels) ; i++)
    { clef = array(fpGels, i, GEL).clef ;
      if (keySetFind (fpSelectedClefs, clef, &dummy))
	max1++ ;
    }

  if (!max1) max1++ ;
  dx = (mapGraphWidth - x - 6)/(1.2 * max1) ;
  minLiveBox= graphBoxStart() ;
  graphBoxEnd() ;
  for (i = 0 ; i < arrayMax(fpGels) ; i++)
    { clef = array(fpGels, i, GEL).clef ;
      if (!keySetFind (fpSelectedClefs, clef, &dummy))
	continue ;
      graphBoxStart() ;
      graphText (name(array(fpGels, i, GEL).clef), x , topMargin + 1) ;

      fp = array(fpGels, i, GEL).fp ;
      bp = array(fpGels, i, GEL).bp ;
      j = arrayMax(bp) ;
      while(j--)
	{ y = array(bp, j, float) ;
	  xb = isBp ? arrayp(bp, j, float) : arrayp(fp, j, float) ;
	  if (fpSelectedBands &&
	      arrayFind(fpSelectedBands, xb, &dummy, floatOrder))
	    { graphColor(RED) ;
	      graphLinewidth(.04) ;
	    }
	  y = MAP2GRAPH(fpMap, *xb) ;
	  if (y > topMargin + 2)
	    graphLine(x , y , x + dx, y) ;
	  graphLinewidth(oldwidth) ;
	  graphColor(BLACK) ;
	}
      graphBoxEnd() ;
      box = graphBoxStart() ;
      graphText (name(array(fpGels, i, GEL).clone), x , topMargin) ;
      graphBoxEnd() ;
      if (i == fpActiveGel)
	graphBoxDraw(box, BLACK, LIGHTBLUE) ;
      x += dx * 1.2 ;
    }
  graphText(isBp ? "Base Pairs" : "mm/10", *offset + 1 , mapGraphHeight - 2) ;
  *offset = x ;

  return;
} /* fpDrawGels */

/**************************************************************/

static MENUOPT fpMenu[] = {
  {graphDestroy,	"Quit"},
  {help,		"Help"},
  {graphPrint,		"Print Screen"},
  {mapPrint,		"Print Whole Map"},
  {fpClearDisplay,	"Clear"},
  {mapColControl,	"Columns"},
  {0, 0}
} ;

static MENUOPT buttonOpts[] = {
  {mapWhole,		"Whole"}, 
  {mapZoomIn,		"Zoom In"},
  {mapZoomOut,		"Zoom Out"},
  {fpAddClone,		"Add Clone"},
  {fpAddKeySet,		"Add KeySet"},
  {fpRmClone,		"Remove Clone"},
  {fpClearDisplay,	"Clear"},
  {fpSwitch,		"Switch"},
  {fpSelectBands,	"Intersect"},
  {fpTransform,		"millimetre/10 <-> base pairs"},
  {cptAllPredictedBands, "predicted bands"},
  {0, 0}
} ;

/**************************************************************/

static void fpDraw(void)
{
  float x1, x2, y1, y2 ;
 
  if (!graphActivate(fpGraph))
    return ;
  
  graphClear () ;
  graphFitBounds (&mapGraphWidth, &mapGraphHeight) ;
  
  halfGraphHeight = 0.5 * (mapGraphHeight - topMargin) ; 

  if (mapGraphHeight < 10 + topMargin)
    { messout ("Sorry, this window is too small for a finger print data") ;
      return ;
    }

 
  fpBoxIndex = arrayReCreate (fpBoxIndex,50,int) ;
  fpActiveBox = 0 ;

  graphButtons (buttonOpts, 2, 3.5, mapGraphWidth * .9) ;
  graphText ("Rounding :", 10, 1) ;
  graphTextEntry (fpRoundingBuffer, 4, 22, 1, fpRounding) ;
  graphBoxDim (0, &x1, &y1, &x2, &y2) ;
  topMargin = y2 + 1 ;
  mapDrawColumns (fpMap) ;
  graphMenu (fpMenu) ;
  
  graphRedraw () ;

  return;
} /* fpDraw */
 
/**************************************************************/

void fpClearDisplay (void)
{
  int n ;
    
  if (!graphActivate(fpGraph))
    return ;

  n = arrayMax(fpGels) ;
  while(n--)
    { arrayDestroy(arr(fpGels, n, GEL).fp) ;
      arrayDestroy(arr(fpGels, n, GEL).bp) ;
    }  
  fpGels = arrayReCreate (fpGels, 12, GEL) ;
  fpAllClefs = keySetReCreate (fpAllClefs) ;
  keySetDestroy (fpSelectedClefs) ; /* don t recreate here */
  fpAllClones = keySetReCreate (fpAllClones) ;
  fpDraw() ;

  return;
} /* fpClearDisplay */
 
/**************************************************************/

void fpDisplay (KEY clone)
{
  int i, n, dummy ;
  Array fp, bp ; 
  KEY clef = 0 ;
  float mx, x1 ;

  if (!fingerPrintData(&fp, &bp, clone, &clef) )
    return ;
  
  if (!keySetExists(fpAllClefs))
    fpAllClefs = keySetCreate () ;
  if (!keySetExists(fpAllClones))
    fpAllClones = keySetCreate () ;
  if (!graphExists(fpGraph))
    { fpGraph = graphCreate (TEXT_FIT,"Finger Printing", .2, .2, .3, .8) ;
      if (!fpGraph)
	return ;
      isBp = TRUE ;
      fpMap = mapCreate (fpDraw) ;
      mapInsertCol (fpMap, 1, TRUE, "Main Bands", fpMainLines) ;
      mapInsertCol (fpMap, 2, TRUE, "Locator", mapShowLocator) ;
      mapInsertCol (fpMap, 3, TRUE, "Scale", mapShowScale) ;
      mapInsertCol (fpMap, 4, TRUE, "Bands", fpDrawGels) ;
      mapCursorCreate (fpMap, 1, 0) ; 
      graphRegister (DESTROY, fpDestroy) ;
      graphRegister (RESIZE, fpDraw) ;
      graphRegister (PICK, fpPick) ;

      /* don't associate a the maps on this graph with any look -
	 This fp-display is not a proper look in the sense that fMap
	 is a proper look, but we still use the mapPackage. */
      graphAssociate (&MAP2LOOK_ASSOC, 0) ;

      graphMenu(fpMenu) ;
      fpGels = arrayReCreate(fpGels, 12, GEL) ;
      array(fpGels, 0, GEL).clone = clone ;
      fpActiveGel = 0 ;
      array(fpGels, 0, GEL).fp = fp ;
      arrayp (fpGels, 0, GEL)->clef = clef ;
      keySetInsert (fpAllClefs, clef) ;
      array (fpGels, 0, GEL).bp = bp ;
      fpMap->mag = 1 ;
      fpMap->centre = 500 ;
      sprintf(fpRoundingBuffer, "%g", rounding) ;
    }
  else
    { graphActivate (fpGraph) ;
      graphPop() ;
      n = keySetMax(fpGels) ;
      while(n--)
	if (arrp(fpGels, n, GEL)->clone == clone)
	  return ;
      n = keySetMax(fpGels) ;	
      fpActiveGel = n ;  
      arrayp (fpGels, n, GEL)->clone = clone ;
      arrayp (fpGels, n, GEL)->clef = clef ;
      keySetInsert (fpAllClefs, clef) ;
      arrp (fpGels, n, GEL)->fp = fp ;
      arrp (fpGels, n, GEL)->bp = bp ;
    }
  
  fpMap->min = 0 ;
  i = arrayMax(fp) ;
  while(i--)
    { if (arr(fp, i, float) > maxFp)
	maxFp = arr(fp, i, float) ;
      if (arr(fp, i, float) > maxBp)
	maxBp = arr(bp, i, float) ;
    }
  keySet(fpAllClones, keySetMax(fpAllClones)) = clone;
  while (fingerPrintData(&fp, &bp, clone, &clef))
    { n = keySetMax(fpGels) ;	
      fpActiveGel = n ;  
      arrayp(fpGels, n, GEL)->clone = clone ;
      arrayp(fpGels, n, GEL)->clef = clef ;
      keySetInsert (fpAllClefs, clef) ;
      arrp(fpGels, n, GEL)->fp = fp ;
      arrp(fpGels, n, GEL)->bp = bp ;
      i = arrayMax(fp) ;
      while(i--)
	{ if (arr(fp, i, float) > maxFp)
	    maxFp = arr(fp, i, float) ;
	  if (arr(fp, i, float) > maxBp)
	    maxBp = arr(bp, i, float) ;
	}
    }
 
  keySetSort (fpAllClefs) ;
  keySetCompress (fpAllClefs) ;
  mx = -1 ;
lao:
  if (!keySetExists(fpSelectedClefs))
    fpSelectedClefs = keySetCopy (fpAllClefs) ;
  keySetSort (fpSelectedClefs) ;
  keySetCompress (fpSelectedClefs) ;

 
  for (i = 0 ; i < arrayMax(fpGels) ; i++)
    { clef = array(fpGels, i, GEL).clef ;
      if (keySetFind (fpSelectedClefs, clef, &dummy))
	{ bp = isBp ? array(fpGels, i, GEL).bp : array(fpGels, i, GEL).fp ;
	  x1 = arr (bp, arrayMax(bp) - 1, float) ;
	  if (x1 > mx) mx = x1 ;
	}
    }
  if (mx == -1)
    { mx = 0 ;
      keySetDestroy (fpSelectedClefs) ;
      goto lao ;
    }


  graphFitBounds (&mapGraphWidth, &mapGraphHeight) ;
  halfGraphHeight = 0.5 * (mapGraphHeight - topMargin) ;

  fpMap->centre = mx / 2 ;
  fpMap->mag = (mapGraphHeight - topMargin - 5) /
    			(1.05 * mx) ;
  fpMap->min = 0 ;
  fpMap->max = mx ;

  fpDraw() ;

  return;
} /* fpDisplay */
  
/**************************************************************/

BOOL fingerPrintDisplay (KEY clone, KEY from, BOOL reuse)
{ fpDisplay (clone) ;
  return graphExists (fpGraph) ;
}

/**************************************************************/

static void cptPredictedBands (KEY clone)
{
  OBJ Clone = bsCreate(clone);
  KEY seqKey ;
  Array dna;
  Array hind3, sau3a , fp ; 
  int from = 0, to = 0, n;
  int x;
  float f;
  KEY predicted;

  if (!Clone) return;

  if (bsGetKey(Clone, _Sequence, &seqKey) &&
      (dna = fMapFindDNA (seqKey, &from, &to)))
    {

      hind3 = arrayCreate(30,int) ;
      sau3a = arrayCreate(30,int) ;
      fp = arrayCreate(30,int) ;
      
      dnacptFingerPrintCompute(dna, from, to, 0, hind3, sau3a, fp);
      
      /* attach to clone */
      lexaddkey("Predicted from dna", &predicted, _VMotif);
      
      bsDestroy(Clone) ;
      if ((Clone = bsUpdate(clone)))
	{
	  bsAddTag (Clone, _Gel) ;
	  bsAddKey (Clone, _bsRight, predicted) ;
	  bsPushObj (Clone) ;
	  bsAddTag (Clone, _Bands) ;
	  /* first value allready stored */
	  for (n = 0 ; n < arrayMax(fp) ; n++)
	    { x = arr(fp, n, int) ; f = x;
	    bsAddData(Clone, _bsRight, _Float, &f) ;
	    }
	  bsSave(Clone) ;
	  Clone = 0 ;
	}		
      
      arrayDestroy(dna);
      arrayDestroy(hind3);
      arrayDestroy(sau3a);
      arrayDestroy(fp);
  }
  bsDestroy(Clone);

  return;
} /* cptPredictedBands */

/**************************************************************/

static void cptAllPredictedBands (void)
{
  int i;
  KEYSET aa = keySetCopy(fpAllClones);
  if (!aa) return;
  fpClearDisplay() ;

  /* for each displayed clone */
  for (i=0; i<keySetMax(aa); i++) {
    cptPredictedBands (keySet(aa, i));
    fpDisplay(keySet(aa,i)) ;
  }
  keySetDestroy(aa);

  return;
} /* cptAllPredictedBands */

/**************************************************************/
/************************ eof *********************************/
 
