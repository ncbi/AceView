/*  File: saucisse.c
 *  Author: mieg
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * $Id: saucisse.c,v 1.5 2007/03/17 01:38:15 mieg Exp $
 * Description: 
     Saucisse algorithm

     Works on DNA or any unsigned char Array
     Counts the number of exact occurence of all words of length wordfLength
     Displays the statistics for length <= wordLength

 * HISTORY:
 * Created: March 1998 (mieg)
 *-------------------------------------------------------------------
 */

#include "acedb.h"
#include "freeout.h"
#include "peptide.h"

#define MAGIC 32491

typedef struct saucisseStruct {int magic, size, wordLength ; Array tree, histo ; int isDna ; } *Saucisse ;

#define  _SAUCISSE_DEF
#include "saucisse.h"

/**********************************************************/
/* construction d'un arbre de mots de longueurs -> 16 */

#define NDIM 24
#define MAXLENGTH 16 

static int code[256] ;

typedef struct { int val, count, length, up, nn[NDIM] ; } NODE ;

/**********************************************************/

static int saucisseAddNode (Array tree, int zcurr, unsigned char cc, Array histo)
{
  int z = zcurr, n = code[cc] ;
  NODE *y, *x = arrp (tree, zcurr, NODE) ;
  
  array(histo, cc, int) ++ ;

  if (n < NDIM)
    {
      z = x->nn[n] ;
      if (!z)
	z = x->nn[n] = arrayMax(tree) ;
      y = arrayp (tree, z, NODE) ;
      x = arrp (tree, zcurr, NODE) ; /* may have moved */
      y->val = cc ; y->length = x->length + 1 ; y->count++ ;
      y->up = zcurr ;
    }

  return z ;
}

/**********************************************************/

typedef struct { int z, n ;} RES ;
static int resOrder (const void *a, const void *b)
{
  const RES *r1 = (const RES *)a, *r2 = (const RES *)b ;
  return r2->n - r1->n ;
}

/**********************************************************/

static void makeCode (Array a, int isDna)
{
  unsigned char *cp = arrp (a, 0, unsigned char) - 1 ;
  int i = arrayMax(a) ;
  Array res = arrayCreate (256, RES) ;
  RES *rp ;

  switch (isDna)
    {
    case 1:
      code[(int)'a']++ ;  code[(int)'t']++ ;  code[(int)'g']++ ;  code[(int)'c']++ ;
      break ;
    case 2:
      i = 128 ;
      while (i--)
	if (i >= 'A' && i <= 'Z' && (pepEncodeChar[i]) != -1)
	  code[i]++ ; 
      break ;
    default:
      while (cp++, i--) code[(int)*cp]++ ;
      break ;
    }

  i = 256; 
  while (i--) 
    {
      rp = arrayp(res, i, RES) ;
      rp->z = i ; rp->n = code[i] ;
    }
  arraySort (res, resOrder) ;
  for (i = 0 ; i < 256 ; i++)
    {
      rp = arrayp(res, i, RES) ;
      code[rp->z] = i ; 
    }
}

/**********************************************************/

void  saucisseFill (Saucisse s, Array a)
{
  int i, j, zc ; unsigned char *cp ;
  NODE *x ;
  static BOOL firstPass = TRUE ;
  
  if (s->wordLength < 1 || s->wordLength > 16) return ;
  if (a->size != 1)
    {
      messerror ("Saucisse only deals with unsigned char arrays") ;
      return ; 
    }

  if (firstPass)
    { 
      firstPass = FALSE ; makeCode(a, s->isDna) ;
    }
    
   for (j = 0 ; j < arrayMax(a) ; j++)
     {
       if (arr(a,j,unsigned char) == 128)  /* forget strings starting on 128 */
	 continue ;
       x = arrayp(s->tree, 0, NODE) ; 
       x->count++ ;
       for (zc = 0, i=j, cp = arrp(a,i,unsigned char) ;
	    i < arrayMax(a) && i < j + s->wordLength ; i++, cp += s->size)
	 if (i > j && code[(int)*cp] >= NDIM)
	   break ;
	 else 
	   zc = saucisseAddNode (s->tree, zc, *cp, s->histo) ;
    }
}

/**********************************************************/
/* export the nn most frequent words */

static void saucisseDoShow (Array tree, Array histo, int isDna, int wordLength, int threshHold)
{
  NODE *x ;
  int i, j, z, mini, nr, n1, nlet ;
  KEYSET ks = keySetCreate () ;
  Array res = arrayCreate (200, RES) ;
  RES *rp ;

  if (!tree || !arrayMax(tree))
    { freeOut ("Empty tree\n") ; goto abort ; }

  for (i = nlet = 0 ; i < arrayMax(histo) ; i++)
    {
      rp = arrayp(res, i, RES) ;
      rp->z = i ; rp->n = arr (histo, i, int) ; nlet += rp->n ;
    }
  arraySort (res, resOrder) ;
  x = arrp (tree, 0, NODE) ; 
  freeOutf ("\n%d words of length %d analysed\n", x->count, wordLength) ;
  freeOutf ("\n %d letters ::\n", nlet) ;
  for (i=j=0 ; i < arrayMax(res) ; i++)
    {
      rp = arrayp(res, i, RES) ;
      if (code[rp->z] >= NDIM) continue ;
      if (rp->n) 
	{
	  switch (isDna)
	    {
	    case 1:
	    case 2:
	      freeOutf ("%d:%c ", rp->n, rp->z) ;
	      break ;
	    default:
	      freeOutf ("%d:%d ", rp->n, rp->z - 128) ;
	      break ;
	    }
	}
      if (j++%10 == 9) freeOutf("\n") ;
    }
  freeOut ("\n Loners (letters not counted above)\n") ;
  for (i=j=0 ; i < arrayMax(res) ; i++)
    {
      rp = arrayp(res, i, RES) ;
      if (code[rp->z] < NDIM) continue ;
      if (rp->n) 
	{
	 switch (isDna)
	    {
	    case 1:
	    case 2:
	      freeOutf ("%d:%c ", rp->n, rp->z) ;
	      break ;
	    default:
	      freeOutf ("%d:%d ", rp->n, rp->z - 128) ;
	      break ;
	    }
	} 
      if (j++%10 == 9) freeOutf("\n") ;
    }

  x = arrp (tree, 0, NODE) ; 
  res = arrayReCreate(res, 200, RES) ;
  arrayMax(res) = 0 ;  nr = 0 ;
  mini = threshHold ;

  freeOutf ("\nTotal count %d, threshold %d\n", x->count, mini) ;
  /* accumulate upwards */
  for (i = arrayMax(tree) - 1 ; i >= 0 ; i--)
    {
      x = arrp (tree, i, NODE) ; 
      n1 = x->count ;
      if (x->up && (n1 < mini || (isDna == 2 && wordLength > x->length )))
	{ x->count = 0 ; }
    }

  for (i = 0, x = arrp (tree, 0, NODE) ; i < arrayMax(tree) ; x++, i++)
    if (x->count)
      {
	rp = arrayp(res, nr++, RES) ;
	rp->z = i ; rp->n = x->count ;
      }
	
  arraySort (res, resOrder) ;

  freeOutf ("\nFrequent multi letter words, cutout = %d\n", mini) ;
  i = arrayMax (res) ;
  rp = arrp (res, 0, RES) - 1 ;
  while (rp++, i--)
    {
      z = rp->z ; j = arrp (tree, z, NODE)->length ; keySetMax(ks) = 0 ;
      if (isDna == 2 && j < wordLength)
	continue ;
      while (z)
	{ x = arrp (tree, z, NODE) ; keySet(ks, j--) = z ; z = x->up ; }
      freeOutf ("%d :: ", rp->n) ;
      for (j=1 ; j < keySetMax (ks) ; j++)
	{
	  x = arrp (tree, keySet (ks, j), NODE) ; 
	  if (isDna)
	    freeOutf (" %c ", x->val) ;
	  else
	  freeOutf (" %d ", x->val - 128) ;
	}
      freeOut ("\n") ;
      if (isDna == 2 && i < arrayMax (res) - 20)
	break ;
    }

 abort:
  keySetDestroy (ks) ;
  arrayDestroy (res) ;
}


/**********************************************************/

void saucisseShow (Saucisse s, int threshHold)
{
  if (s && s->magic == MAGIC)
    saucisseDoShow (s->tree, s->histo, s->isDna, s->wordLength, threshHold) ;
}

Saucisse saucisseCreate (int dim, int wordLength, int isDna)
{
  Saucisse s = messalloc (sizeof (struct saucisseStruct)) ;

  s->magic = MAGIC ;
  s->histo = arrayCreate (256, int) ;
  s->tree = arrayCreate (dim, NODE) ;
  s->size = 1 ;  /* saucisse only handles unsigned char or DNA */ 
  if (wordLength < 1 || wordLength > 16) 
    {
      messout ("wordLength %d not in range [1-16], adjusted to 5", wordLength) ;
      wordLength = 5 ;
    }
  s->wordLength = wordLength ;
  s->isDna = isDna ;
 
  return s ;
} 

void  uSaucisseDestroy (Saucisse s) 
{
  if (s && s->magic == MAGIC)
    {
      arrayDestroy (s->tree) ;
      arrayDestroy (s->histo) ;
      
      messfree (s) ;
    }
}

/**********************************************************/
/**********************************************************/
#ifdef JUNK

void saucisseTest (KEYSET ks)
{
}

#else

#include "dna.h" /* just needed in saucisseTest */
#include "peptide.h"

void saucisseTest (KEYSET ks)
{
  int i, threshHold = 0 ; char *cp ;
  Array test = 0 ;
  Saucisse saucisse = 0 ;
  int mode = 1, wordLength = 0 ; /* 0: scf trace, 1, 2 dna, 3 pep */

  if (!messPrompt("mode[0-3] (0: scf trace, 1: give one dna, 2: dna active keyset, 3: pep active set)   wordLength[1-16]  threshold[1-1000]","1 5 5","iii"))
    return ;
  freeint(&mode) ;  freeint(&wordLength) ; freeint(&threshHold) ;
  if (wordLength < 1 || wordLength > 16) 
    {
      messout ("wordLength must be between 1 and 16 included") ;
      return ;
    }
  if (threshHold < 0 || threshHold > 1000)
    {
      messout ("threshHold must be between 0 and 1000 included") ;
      return ;
    }
  saucisse =  saucisseCreate (2000, wordLength, TRUE) ;
  switch (mode)
    {
    case 1: /* enter your own dna in interactive mode */

      while (messPrompt ("sequence", "atatgctagtcta","w"))
	{
	  cp = freeword() ;
	  if (cp)
	    {
	      i = strlen(cp) ;
	      test = arrayReCreate (test,i, unsigned char) ;
	      arrayMax(test) =  strlen(cp) ;
	      memcpy (test->base, cp, i) ;
	      saucisseFill (saucisse, test) ;
	    }
	}
      arrayDestroy (test) ;
      break ;
    case 2: /* work on DNA keyset */
      if (!keySetMax(ks))
	{ 
	  messout ("Please select a keyset of sequence") ;
	  return ;
	}
  
      saucisse =  saucisseCreate (2000, wordLength, TRUE) ;
      for (i = 0 ; i < keySetMax(ks) ; i++)
	if ((test = dnaGet (keySet(ks,i))))
	  {
	    dnaDecodeArray (test) ;
	     saucisseFill (saucisse, test) ;
	    arrayDestroy (test) ;
	  }
      break ;
    case 3: /* work on Peptide keyset */
      if (!keySetMax(ks))
	{ 
	  messout ("Please select a keyset of sequence") ;
	  return ;
	}
  
      saucisse =  saucisseCreate (2000, wordLength, 2) ;
      for (i = 0 ; i < keySetMax(ks) ; i++)
	if ((test = peptideGet (keySet(ks,i))))
	  {
	    pepDecodeArray (test) ;
	     saucisseFill (saucisse, test) ;
	    arrayDestroy (test) ;
	  }
      break ;
    }  
  saucisseShow (saucisse, threshHold) ;
  saucisseDestroy (saucisse) ;
}

#endif

/**********************************************************/
/**********************************************************/
/*

Exemple (from scf trace files)

histogram of values

115778:8 58736:120 43364:7 42680:9 24328:10 23568:6 17336:134 9280:11
8864:5 4096:4 3096:12 1888:133 1568:135 1280:13 1264:132 896:3 816:131
744:14 504:119 496:130 480:2 424:15 360:129 352:1 320:128 232:16
216:118 208:127 200:17 168:18 168:0 152:125 144:126 114:19 88:117
80:123 80:124 56:116 56:122 48:20 40:108 40:27 35:21 32:36 32:114
32:113 32:104 32:23 32:109 32:33 32:34 25:22 24:111 24:110 24:107
24:101 24:52 24:115 24:96 24:67 24:29 24:26 24:86 16:25 16:53 16:65
16:51 16:56 16:45 16:69 16:39 16:28 16:89 16:38 16:74 16:103 16:37
16:85 16:78 16:121 16:102 16:81 16:97 16:106 8:66 8:63 8:30 8:80 8:88
8:35 8:90 8:71 8:60 8:41 8:94 8:59 8:58 8:75 8:32 8:57 8:100 8:55 8:42
8:72 8:54 8:105 8:93 8:50 8:40 8:48 8:47 8:46 8:112 8:79 8:43 8:83
0:99 0:98 0:95 0:92 0:91 0:87 0:84 0:77 0:76 0:73 0:70 0:64 0:62 0:24
0:61 0:49 0:44 0:82 0:31 0:68

most frequent words
123945 words  :: 
 11377 ::  4  3  5  6  10  12 
11197 ::  7  7  8  6  9  5  10 
10988 ::  8  8  10  7  10  8  8 
9005 ::  12  2  1  5  10  8  8 
6658 ::  10  8  14  6  10  8  8 
6247 ::  8  12  8  9  8  8  8 
 // + 12 = 12 seconds

*/

/**********************************************************/
/**********************************************************/
