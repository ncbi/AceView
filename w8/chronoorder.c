/*  File: chronoorder.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@frmop11.bitnet
 *
 * Description:
   This well known algorithm tries to reorder a matrix containing 
   TRUE/FALSE descriptors untill all ther TRUE are as much as possible
  in the central diagonal.
  In this way, one reorders concurently the lines and columns.

  As far as I know, this algorithm is very fast and very robust but its convergence is not
  garanteed on random data.

* Implementation details:

   ATTENTION: What is called Ligne in the code appears as column in the display. 
   mat(i,j) is given as input and the matrix itself is never touched
   il, ic are arrays of int that get reordered and point on the lines
   and col of mat(i,j)
   wl, wc are weights attached to the lines and col of mat()
   wl[2] is always attached to the line mat(2,*) but its
   evaluation depens on ic(*)

   I treat differently the positives results 2 or *
   from the unknown 0 or . and negatives 1 or -
   I clean up solitary 2 recursively, 

   Also, the repulse parameter expels the - from the vertical
   columns on display called Lignes in the code.
 * Exported functions:
 * HISTORY:
 * Last edited: Mar  6 13:47 1995 (mieg)
 * Created: Fri Oct 30 19:10:15 1992 (mieg)
 *-------------------------------------------------------------------
 */

 * Last edited: Mar  6 12:57 1995 (mieg)

/*  $Id: chronoorder.c,v 1.2 2003/03/17 05:04:34 mieg Exp $ */

#include "acedb.h"
#include "keyset.h"
#include "../wh/menu.h"
#include "graph.h"
#include "topology.h"
#include "interval.h"

static int nCol = 1, nLignes = 1 ;
#define mat(i,j) (arr(m,  (i)*nCol + (j), unsigned char))
static Array wc = 0 , wl = 0 , ic = 0, il = 0, m , boxing = 0 , cost = 0 , nameIndex ;
static Array colO, linO ; /* Tableaux de depart */
static Stack nameStack = 0 ;

static float t = 0 , t1 = 0 , t2 = 0 ;
static int b0, b1, b2 ;
static char tBuffer[8] , t1Buffer[8] , t2Buffer[8] ;
static void rotate(void) ;
static void boxingOut(void) ;
static void showMat(void) ;
static void refine(void) ;
static void export(void) ;
static void shouldWeDestroy(void) ;
static void interToDef(void) ;
static int matColOrder(void *a, void *b) ;

static void simplify(void)
{ int i, j, n ;

  if (!messPrompt("Restrict to the top square of size ","30","i"))
    return ;
  freeint(&n) ;
  for (i = n ; i < nLignes ; i++)
    for (j = n ; j < nCol ; j++)
      mat(i,j) |= F_SIMPLIFY ;
}

static MENUOPT matMenu[] =
            {   shouldWeDestroy, "Quit",
		help,"Help",
		graphPrint,"Print",

		rotate, "Rotate",
		refine, "Refine",
		simplify, "Simplify",
		interToDef, "Save Order",  
		export, "AceDump",
		0, 0
            } ;

static float minWeight = 0 , maxWeight = 0 ;

static void weightCol(int x1, int y1, int x2, int y2)
{ int i, n, j = y2 , i1, j1, k1, k2 ;
  float p = 0 , p1 ;
   
  while (j-- > y1)
    { n = 0 ;
      p = 0 ;
      j1 = arr(ic, j, int) ;
      i = x2 ;
      
      while (i-- > x1)
	{ i1 = arr(il, i, int) ;
	  if (mat(i1,j1) == F_YES)
	    { n += i + 1 + t1*
		 (arr(wl,i1,float)  - minWeight) ;
	      p++ ; 
	    }
	}

      if(t2)
	{
	  i = x2 ; k2 = x1 ;
	  while(i-- > x1)
	    { i1 = arr(il, i, int) ;
	      if (mat(i1,j1) == F_YES)
		{ k2 = i ;
		  break ;
		}
	    }
	  for (k1 = x2 , i = x1 ; i < x2 ; i++)
	    { i1 = arr(il, i, int) ;
	      if (mat(i1,j1) == F_YES)
		{ k1 = i ;
		  break ;
		}
	    }
	  for (i = k1 + 1 ; i < k2 ; i++)
	    { i1 = arr(il, i, int) ;
	      if (mat(i1,j1) == F_NO)
		{ n += t2 *
		    (2*i > k1 + k2 ? 1 : -1 ) ;
		  
		}
	    }
	}
	  
      p1 = arr(wc, j1, float) = minWeight +
	(p>0 ? n/p + randfloat()*t  : 0) ;
      if (p1 > maxWeight)
	maxWeight = p1 ;
    }
}

static void weightLignes(int x1, int y1,  int x2, int y2)
{ int  i = x2, i1, n,  j , j1 ;
  float p = 0 , p1 ;
  char c = 2 ;

   while (i-- > x1)
    { n = 0 ;
      p = 0 ;
      i1 = arr(il, i, int) ;
      j = y2 ;
      
      while (j-- > y1)
	{ j1 = arr(ic, j, int) ;
	  if (mat(i1,j1) == c)
	    { n += j  + 1 + t1*
		(arr(wc,j1,float) - minWeight) ; 
	      p++ ; 
	    }
	}
      p1 = arr(wl, i1, float) = minWeight +
	( p>0 ? n/p + randfloat()*t : 0 ) ;
      if (p1 > maxWeight)
	maxWeight = p1 ;
    }
}

static void findCost(int x1, int y1,  int x2, int y2)
{ int i, n, j = y2 , i1, j1 , old , k1, k2 ;
 
  while (j-- > y1)
    { n = 0 ;
      old = 9877  ;
      j1 = arr(ic, j, int) ;
      i = x2 ; 
 
      i = x2 ; k2 = x1 ;
      while(i-- > x1)
	{ i1 = arr(il, i, int) ;
	  if (mat(i1,j1) != F_NO)
	    break ;
	}
      i = x2 ; k2 = x1 ;
      while(i-- > x1)
	{ i1 = arr(il, i, int) ;
	  if (mat(i1,j1) == F_YES)
	    { k2 = i ;
	      break ;
	    }
	}
      for (k1 = x2 , i = x1 ; i < x2 ; i++)
	{ i1 = arr(il, i, int) ;
	  if (mat(i1,j1) == F_YES)
	    { k1 = i ;
	      break ;
	    }
	}
      for (i = k1 + 1 ; i < k2 ; i++)
	    { i1 = arr(il, i, int) ;
	      if (mat(i1,j1) == F_NO)
		{ n += t2 *
		    (2*i > k1 + k2 ? 1 : -1 ) ;
		  
		}
	    }
	     while(i-- > x1)
	{ i1 = arr(il, i, int) ;
	  if ((mat(i1,j1) == F_NO ||
	       mat(i1,j1) == F_YES) &&
	      mat(i1,j1) != old)
	    { n++ ;
	      old = mat(i1,j1) ;
	    }
	}
      array(cost, j1, int) = n ;
    }
}

static int insertCost(int k1, int j2, int x1, int y1,  int x2, int y2)
{ int i, n, j = j2 , i1, j1 , old ;
 
  while (j-- > y1)
    { n = 0 ;
      old = 9877  ;
      j1 = arr(ic, j, int) ;
      i = x2 ; 
      while(i-- > x1)
	{ i1 = arr(il, i, int) ;
	  if ((mat(i1,j1) == F_YES ||
	       mat(i1,j1) == F_NO ) &&
	      mat(i1,j1) != old)
	    { n++ ;
	      old = mat(i1,j1) ;
	    }
	}
      array(cost, j1, int) = n ;
    }
  return n ;
}


static void expell(int x1, int y1,  int x2, int y2)
{
  int i, n, j = y2 , i1, j1 ;
 
  while (j-- > y1)
    { n = 0 ;
      j1 = arr(ic, j, int) ;
      if (arr(cost, j1, int) <= 2)
	continue ;
      i = x2 ; 
      while(i-- > x1)
	{ i1 = arr(il, i, int) ;
	  mat(i1,j1) |= F_EXPELL ;
	}
    }
}

static void doInsert(int j2, int k2)
{ int j , j1 ;
  j = nCol ;
  while(j--)
    { 
      j1 = arr(ic, j, int) ;
      arr(wc, j1, float) = j ;
    }
  j1 = arr(ic, j2, int) ;
  arr(wc, j1, float) = k2 - .5 ;
  arraySort(ic, matColOrder) ;
}

static void reinstall(int x1, int y1,  int x2, int y2)
{
  int  n, j = y2 , j1 , max = -2 , max2, k1, k2 ;
 
 lao:
  max = 2 ;
  while (j-- > y1)
    { 
      j1 = arr(ic, j, int) ;
      if (arr(cost, j1, int) > max)
	{ k1 = j ;
	  max = arr(cost, j1, int) ;
	}
    }
  if (max = 2)
    return ;
  max2 = 100000;
  j = y2 ; k2 = y1 ;
  while(j-- > y1)
    { n = insertCost(k1, j, x1, y1, x2, y2) ;
      if (n > max2)
	{ max2 = n ; 
	  k2 = j ;
	}
    }
  doInsert(j, k2) ;
  goto lao ;  
}


static void refine(void)
{ int i, i1, i2, j1, j2 ;
  
  cost = arrayReCreate(cost, nCol, int) ;
  
  for (i = 0 ; i < arrayMax(boxing) ; )
    { minWeight = maxWeight + 1 ;
      i1 = arr(boxing, i++, int) ;
      j1 = arr(boxing, i++, int) ;
      i2 = 1 + arr(boxing, i++, int) ;
      j2 = 1 + arr(boxing, i++, int) ;
	
      findCost(i1, j1, i2, j2) ;
      expell(i1, j1, i2, j2) ;
      reinstall(i1, j1, i2, j2) ;
    }
  showMat() ;
}

static void showMat(void)
{ int i, i1, j, j1 , i2, j2 , nn = 7 ;

  graphClear() ;
  i = nLignes ; j = nCol ;
  graphButtons(matMenu, 3,4, 40) ;
  while(i--)
    { i1 = arr(il, i, int) ;
      j = nCol ;
      while (j--)
	{ j1 = arr(ic, j, int) ;
	  switch (mat(i1, j1))
	    { 
	    case 0:
	      break ;
	    case F_NO:
	      graphText("|", 5 + 2*i, nn + j) ;
	      break ;
	    case F_YES:
	      graphText("*", 5 + 2*i, nn + j) ;
	      break ;
	    default:
	      if(mat(i1,j1) & F_YES_NO)
		graphText(".", 5 + 2*i, nn + j) ;
	      break ;
	    }
	}
    }

  boxingOut() ;
  graphColor(RED) ;
  for (i = 0 ; i < arrayMax(boxing) ; )
    { i1 = arr(boxing, i++, int) ;
      j1 = arr(boxing, i++, int) ;
      i2 = arr(boxing, i++, int) ;
      j2 = arr(boxing, i++, int) ;
      graphRectangle(5 + 2*i1, nn + j1, 6 + 2*i2, nn+1 + j2) ;
    }
  graphColor(BLACK) ;

  t = t - (int)t ;
  if (t<0) t *= -1 ;
  sprintf(tBuffer,"%d",(int)(100*t + .1)) ;
  graphText("Temperature:", 1,1) ;
  b0 = graphTextEntry(tBuffer,6,15,1,0) ;

  if (t1<0) t1 *= -1 ;
  sprintf(t1Buffer,"%d",(int)(100*t1 + .1)) ;
  graphText("Drag:", 23,1) ;
  b1 = graphTextEntry(t1Buffer,6,29,1,0) ;

  if (t2<0) t2 *= -1 ;
  sprintf(t2Buffer,"%d",(int)(100*t2 + .1)) ;
  graphText("Repulse:",38,1) ;
  b2 = graphTextEntry(t2Buffer,6,47,1,0) ;

  for (i=0; i < nLignes; i++)
    { i1 = arr(il, i, int) ;
      graphText(messprintf("%d:%4.1f",
			 i1,
			 arr(wl, i1, float)), 
	      8*i, 2) ;
    }
 
 for (i=0; i < nCol; i++)
    { i1 = arr(ic, i, int) ;
      graphText(messprintf("%d:%4.1f",
			 i1, 
			 arr(wc, i1, float)), 
	      8*i, 3) ;
    }

  if(nameStack)
    for (j=0; j < nCol; j++)
      { j1 = arr(ic, j, int) ;
	graphText(stackText(nameStack, array(nameIndex,j1,int)), 
		  8 + 2*nLignes, nn + j) ;
      }

  graphTextBounds(8 + 2*nLignes, nn + 3 + nCol) ;
  graphRedraw() ;
}
 
static void export(void)
{ FILE *f ; 
  int x, i, i1, j, j1 ;
  char dirname[DIR_BUFFER_SIZE], filname[FIL_BUFFER_SIZE] ;

  if(!nameStack)
    { messout("Names are missing, I can t export") ;
      return ;
    }

  strcpy(filname,"chrono") ;
  f = filqueryopen(dirname, filname,"ace","w","File to export the results") ;
  if (!f)
    return ;
  
  
  for (j=0; j < nCol; j++)
    { j1 = arr(ic, j, int) ;
      for (i=0; i < nLignes ; i++)
	{ i1 = arr(il, i, int) ;
	  if (mat(i1,j1) == F_YES) 
	    break ;
	}
      if (i < nLignes)
	fprintf(f, "Locus %s\n Map %s Position %d\n\n", 
		stackText(nameStack, array(nameIndex,j1,int)), filname, j) ;
    }
  
  for (i=0 ; i < nLignes ; i++)
    { i1 = arr(il, i, int) ;
      
      for (j=0; j < nCol; j++)
	{ j1 = arr(ic, j, int) ;
	  if (mat(i1,j1) == F_YES) 
	    break ;
	}
      if (j < nCol)
	{ x = j ;
	  fprintf(f, "Interval %s\n Map %s Left %d\n", 
		  stackText(nameStack, array(nameIndex,nCol + i1,int)), 
		  filname, x) ;
	  for (j=x; j < nCol; j++)
	    { j1 = arr(ic, j, int) ;
	      if (mat(i1,j1) == F_YES) 
		x = j ;
	    }
	  fprintf(f, "Map %s Right %d\n\n", filname, x) ;
	}
    }  
  filclose(f) ;
}

static int matLineOrder(void *a, void *b)
{ int i = *(int*)a, j = *(int*)b ;
  float x = arr(wl, i, float), y = arr(wl,j,float) ;
  return
    x < y ? -1
      : (x == y ? 0 : 1) ;
}

static int matColOrder(void *a, void *b)
{ int i = *(int*)a, j = *(int*)b ;
  float x = arr(wc, i, float), y = arr(wc,j,float) ;
  return
    x < y ? -1
      : (x == y ? 0 : 1) ;
}

static void zeroWeight (BOOL c) 
{ int i ;
  if (c)
    { i = nLignes ;
      while (i--)
	arr(wl, i, float) = 0. ;
    }
  else
    { i = nCol ;
      while (i--)
	arr(wc, i, float) = 0. ;
    }
}

static void jetteLesZeros (BOOL c) 
{ int i ;
  if (c)
    { i = nLignes ;
      while (i--)
	if (!arr(wl, i, float))
	  arr(wl, i, float) = maxWeight ;
    }
  else
    { i = nCol ;
      while (i--)
	if (!arr(wc, i, float))
	  arr(wc, i, float) = maxWeight ;
    }
}

static BOOL jetteLigne(int x1, int y1,  int x2, int y2)
{ int  i = x2, i1, i2,  j , j1 ;
  BOOL found = FALSE ;

  while (i-- > x1 + 1)
    {
      i1 = arr(il, i, int) ;
      i2 = arr(il, i - 1, int) ;
      j = y2 ;
      
      while (j-- > y1)
	{ j1 = arr(ic, j, int) ;
	  if (mat(i1,j1) &&
	      mat(i1,j1) != mat(i2,j1))
	    goto different ;
	}
      arr(wl, i1, float) = maxWeight ;
      j = y2 ;
       while (j-- > y1)
	{ j1 = arr(ic, j, int) ;
	  mat(i1,j1) |= F_JETTE_LIGNE ;
	}
     found = TRUE ;
    different:
      continue ;
    }
  return found ;
}

static BOOL jetteCol(int x1, int y1,  int x2, int y2)
{ int  i , i1,  j = y2 , j1 , j2;
  BOOL found = FALSE ;

  while (j-- > y1 + 1)
    {
      j1 = arr(ic, j, int) ;
      j2 = arr(ic, j - 1, int) ;
      i = x2 ;
      
      while (i-- > x1)
	{ i1 = arr(il, i, int) ;
	  if (mat(i1,j1) &&
	      mat(i1,j1) != mat(i1,j2))
	    goto different ;
	}
      arr(wc, j1, float) = maxWeight ;
      i = x2 ;
       while (i-- > x1)
	{ i1 = arr(il, i, int) ;
	  mat(i1,j1) |= F_JETTE_COL ;
	}
     found = TRUE ;
    different:
      continue ;
    }
  return found ;
}

static void jetteLesDoubles (BOOL c) 
{ BOOL foundDouble = FALSE ;
  int i , i1, j1, i2, j2 ;
  for (i = 0 ; i < arrayMax(boxing) ; )
    { minWeight = maxWeight + 1 ;
      i1 = arr(boxing, i++, int) ;
      j1 = arr(boxing, i++, int) ;
      i2 = 1 + arr(boxing, i++, int) ;
      j2 = 1 + arr(boxing, i++, int) ;
      
      if (c)
	foundDouble |= jetteLigne(i1, j1, i2, j2) ;
      else
	foundDouble |= jetteCol(i1, j1, i2, j2) ; 
    }

  if (foundDouble)
    { 
      if (c)
	arraySort(il, matLineOrder) ;
      else
	arraySort(ic, matColOrder) ;
    }
}

static void rotate(void)
{ int tt , tt1, tt2 ;
  int i, i1,  j1 , i2, j2 ;
  float t0 ;
  static BOOL c = TRUE ;

  freeforcecard(t1Buffer) ;
  if (freeint(&tt1) && tt1 >= 0 && tt1 <= 1000)
    t1 = tt1/100.0 ;

  freeforcecard(t2Buffer) ;
  if (freeint(&tt2) && tt2 >= 0 && tt2 <= 1000)
    t2 = tt2/100.0 ;

  freeforcecard(tBuffer) ;
  if (freeint(&tt) && tt >= 0 && tt <= 100)
    t = tt/100.0 ;

  t0 = t ;
  tt = (int)(100*t) % 100 ;
  do
    { boxingOut() ;
      c = !c ;
      minWeight = maxWeight = 0 ;
      zeroWeight (c) ;
      for (i = 0 ; i < arrayMax(boxing) ; )
	{ minWeight = maxWeight + 1 ;
	  i1 = arr(boxing, i++, int) ;
	  j1 = arr(boxing, i++, int) ;
	  i2 = 1 + arr(boxing, i++, int) ;
	  j2 = 1 + arr(boxing, i++, int) ;
	  
	  if (c)
	    weightLignes(i1, j1, i2, j2) ;
	  else
	    weightCol(i1, j1, i2, j2) ;
	}
      maxWeight++ ;
      jetteLesZeros(c) ;
      if (c)
	arraySort(il, matLineOrder) ;
      else
	arraySort(ic, matColOrder) ;
      jetteLesDoubles(c) ;
      tt /= 2 ;
      t = tt/100.0 ;
    }  while(tt > 0) ;    /* do while construct */
      
  t = t0 ;  /* reestablish for next round  */
  showMat() ;
}

static void  localPick(int box)
{
 if (box == b0)
   graphTextEntry(tBuffer,0,0,0,0) ;
 else if (box == b1)
   graphTextEntry(t1Buffer,0,0,0,0) ;
 else if (box == b2)
   graphTextEntry(t2Buffer,0,0,0,0) ;
}



static BOOL cleanCol(void)
{ int  i, i1, j = nCol, nn = 0 , r = 0 ;
  
  
   while (j--)
     { i = nLignes ;
       nn = 0;

       while (i--)
	 if (mat(i,j))
	   { i1 = i ;
	     nn++ ;
	   }
       
       if (nn == 1 &&
	   mat(i1,j) == 2 )
	 { r++; mat(i1,j) |= F_COL_ISOLE ;}
     }
  return r > 0 ;
}


static BOOL cleanLignes(void)
{ int  i = nLignes, j, j1, nn = 0 , r = 0 ;
/*  
  i = nLignes ; ip = arrp(il, 0,int) - 1 ;
  while(ip++, i--)
  { j = nCol ;
  jp = arrp(ic,i,int) -1 ;
  nn = 0 ;
  while(jp++, j--)
  if(mat(*ip,*jp))
  nn++ ;
  if(nn < 2)
  { push(trivial, 1) ;
  push(trivial, *ip) ;
  }
  */     
  
  while (i--)
     { j = nCol ;
       nn = 0;

       while (j--)
	 if (mat(i,j) > 0)
	   { j1 = j ;
	     nn++ ;
	   }
       
       if (nn == 1 &&
	   mat(i,j1) == 2 )
	 { r++ ; mat(i,j1) |= F_L_ISOLEE ; }
     }
  return r > 0 ;
}

static void boxingOut(void)
{ int i, j, i1, j1, i2, j2, n = 0 ;
  BOOL found ;
  
  boxing = arrayReCreate(boxing, 12, int) ;
  i1 = j1 = 0 ;
 
lao:      /* look for top corner of box */
      /* Start at top corner , scan lines and try to increase i1 */
  found = FALSE ;
  for (i = i1 ; i < nLignes ; i++)
    { found = FALSE ;
      for (j = j1 ; j < nCol ; j++)
	if (mat(arr(il,i,int), arr(ic,j,int)) == F_YES)
	  { found = TRUE ;
	    i2 = i ;
	    break ;
	  }
      if (found)
	  break ;
    }
  if (!found)
    return ;
  i1 = i ;
    
  for (j = j1 ; j < nCol ; j++)
    { found = FALSE ;
      for (i = i1 ; i < nLignes ; i++)
	if (mat(arr(il,i,int), arr(ic,j,int)) == F_YES)
	  { found = TRUE ;
	    j2 = j ;
	    break ;
	  }
      if (found)
	  break ;
    }
  j1 = j ;

  array(boxing, n++, int) = i1 ;
  array(boxing, n++, int) = j1 ;


     /* Now look bottom corner of box */
  found = TRUE ;
  while(found)
    { found = FALSE ;
      for (i = i1 ; i <= i2 ; i++)
	{ 
	  for (j = j1 ; j < nCol ; j++)
	    if (mat(arr(il,i,int), arr(ic,j,int)) == F_YES && j > j2)
	      { found = TRUE ;
		j2 = j ;
	      }	
	}
      for (j = j1 ; j <= j2 ; j++)
	{ 
	  for (i = i1 ; i < nLignes ; i++)
	    if (mat(arr(il,i,int), arr(ic,j,int)) == F_YES && i > i2)
	      { found = TRUE ;
		i2 = i ;
	      }	  
	}
    }

  array(boxing, n++, int) = i2 ;
  array(boxing, n++, int) = j2 ;

  i1 = i2 + 1 ; j1 = j2 + 1 ;
  goto lao ;
}

/**********************************************/

static void interToDef(void)
{
  int i ;
  
      i = nLignes ;
      while (i--)
	arr(linO, i, int) = arr(il, i, int) ;
      i = nCol ;
      while(i--)
	arr(colO, i, int) = arr(ic, i, int) ;
}

/**********************************************/

static void localDestroy(void)
{ 
  arrayDestroy(nameIndex) ;
  stackDestroy(nameStack) ;
}

/**********************************************/

static void shouldWeDestroy(void)
{ if ( messQuery("Do you really want to quit chrono"))
    graphDestroy() ;
}
	  
/**********************************************/

void chronologicalOrdering(Array mm, Array cocolO, Array lilinO, Stack names)
{ int i = 0, nc, nl ;

  linO = lilinO ;
  colO = cocolO ;
  nl = arrayMax(linO) ;
  nc = arrayMax(colO) ;
  wl = arrayReCreate(wl, nl, float) ;
  wc = arrayReCreate(wc, nc, float) ;
  if (arrayExists(il))
    arrayDestroy(il) ;
  if (arrayExists(ic))
    arrayDestroy(ic) ;
  il = arrayCopy(linO) ;
  ic = arrayCopy(colO) ;
  nCol = nc ;
  nLignes = nl ;
  nameStack = names ;
  m = mm ;
  
  nameIndex = arrayCreate(nl, int) ;

  if(names)
    { stackCursor(nameStack,0) ; i = 0 ;
      do { array(nameIndex,i++,int) = stackPos(nameStack) ;
	 } while(stackNextText(nameStack)) ;
    }

  array(wl, nLignes - 1, float) = 0 ;  /* make room */
  array(wc, nCol - 1, float) = 0 ;

  graphCreate(TEXT_SCROLL,"Interval Map", .1, .4, .6, .7) ;
  graphRegister(DESTROY, localDestroy) ;
  graphRegister(PICK, localPick) ;
  graphMenu(matMenu) ;

  while (cleanCol() && cleanLignes());
  showMat() ;
}
