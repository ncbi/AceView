/*  File: basepad.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  4 13:31 1998 (fw)
 * Created: Thu Feb  1 23:49:22 1996 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: basepad.c,v 1.4 2014/11/10 18:05:09 mieg Exp $ */

/* a project can have one and only one padded assembly
   there should be: 

   ?Clone Padded_assembly UNIQUE ?Sequence

   where Sequence Is_assembly
   All this should be maintained by these functions
*/

#include "acedb.h"

#include "dna.h"
#include "a.h"

#ifdef ACEMBLY  
#include "dnaalign.h"		/*  for alignToolsAdjustLink() */
#endif /* ACEMBLY */
 
static KEY _Padded, _Unpadded, _Quality, _SCF_Position ;
static KEY _Align_to_SCF ;

/* 
static int idSize = 0 ;
static int *identity = 0 ;		// used for identity index2 
static void idCheck (int max)
{ 
  int i ;

  ++max ;
  if (max > idSize)
    { idSize = 2*max ;
      messfree (identity) ;		// OK to messfree (0) 
      identity = (int*) messalloc (idSize*sizeof(int)) ;
      for (i = 0 ; i < idSize ; ++i)
	identity[i] = i ; 
    }
}
*/

static void initialise (void)
{ 
  static BOOL done = FALSE ;

  if (!done)
    { lexaddkey ("Padded", &_Padded, 0) ;
      lexaddkey ("Unpadded", &_Unpadded, 0) ;
      lexaddkey ("Clipping", &_Clipping, 0) ;
      lexaddkey ("Quality", &_Quality, 0) ;
      lexaddkey ("SCF_Position", &_SCF_Position, 0) ;
      lexaddkey ("Assembled_from", &_Assembled_from, 0) ;
      lexaddkey ("Subsequence", &_Subsequence, 0) ;
      lexaddkey ("Align_to_SCF", &_Align_to_SCF, 0) ;
      lexaddkey ("DNA", &_DNA, 0) ;

      done = TRUE ;
    }
}

/*******************************************************/

/* Gobble up non indexed characters in a char array 

   throughout this, index[] runs over indices 1..n as used in the
   sequence object, not 0..(n-1) as usual in C
*/

void fixArray (OBJ obj, KEY tag, int *index)
{ 
  register int j ;
  KEY key ;
  Array old, new ;

  if (bsGetKey (obj, tag, &key) && 
      (old = arrayGet (key, char, "c")))
    { j = arrayMax(old)-1 ;
      new = arrayCreate(index[j+1], char) ;
      array(new,index[j+1]-1,char) = arr(old,j,char) ;
      while (j--)   /* must run backwards to overwrite pad values */
	arr(new,index[j+1]-1,char) = arr(old,j,char) ;
      arrayStore (key, new, "c") ;
      { int len = arrayMax(new) ;
	bsAddData (obj, _bsRight, _Int, &len) ;
      }
      arrayDestroy (old) ; arrayDestroy (new) ;
    }
}


/* Nearly identical, but dnaStoreDestroy is special */

void fixDNA (KEY seq, Array s, int *index)
{ 
  register int j ;
  KEY key ;
  OBJ obj = 0 ;
  Array snew ;

  if ((obj = bsCreate (seq)) && bsGetKey (obj, _DNA, &key))
    { snew = arrayCreate(index[arrayMax(s)], char) ;
      j = arrayMax(s)-1 ; 
      array(snew,index[j+1]-1,char) = arr(s,j,char) ;
      while (j--) 
	arr(snew,index[j+1]-1,char) = arr(s,j,char) ;

      dnaStoreDestroy (key, snew) ;
    }

  bsDestroy (obj) ;
  arrayDestroy (s) ;
}

/***********************************************/
#ifdef JUNK
static int *depadSequence (KEY seq, BOOL defragment)
{ 
  Array s = dnaGet (seq) ;
  OBJ obj = bsUpdate (seq) ;
  int *index = 0 ;
  KEY key ;
  int len ;
  register int i, j ;
  
  { KEYSET testKs = queryKey (seq, ">Assembled_from  Padded") ;
    i = keySetMax (testKs) ;
    keySetDestroy (testKs) ;
  }

  if (!obj)
    goto done ;

  if (!s  || (!i && !bsFindTag (obj, _Padded)))
    goto do_children ;

				/* build index */
  index = (int*) messalloc ((1 + arrayMax(s)) * sizeof(int)) ;
  for (i = 0, j = 0 ; i < arrayMax(s) ; ++i)
    { if (arr(s,i,char)) /* jump character zero (the pad) */
	++j ; 
      index[i+1] = j ;
    }

  if (j < i)  /* pads removed from this object */
    { bsSave (obj) ;
      fixDNA (seq, s, index) ; s = 0 ;
      obj = bsUpdate (seq) ;
      fixArray (obj, _Quality, index) ;
      fixArray (obj, _SCF_Position, index) ;
      bsCoordIndex (obj, index) ;
    }

  /* now do Assembled_from and Align_to_SCF, i.e. relationship entries */

  { Array oldA = arrayCreate (8192, BSunit) ;
    Array newA = arrayCreate (32768, BSunit) ;
    int *index2, i1, i2 ;

				/* do Assembled_from - the hard bit */
    if (bsGetArray (obj, _Assembled_from, oldA, 5))
      { for (i = 0, j = 0 ; i < arrayMax(oldA) ; i += 5)
	  { key = arr(oldA,i,BSunit).k ;
	    index2 = depadSequence (key, defragment) ;	/* RECURSION */
	    if (!index2) 
	      { idCheck (arr(oldA,i+4,BSunit).i) ;
		index2 = identity ;
	      }
	    i1 = arr(oldA,i+1,BSunit).i ;
	    i2 = arr(oldA,i+3,BSunit).i ;
	    if (arr(oldA,i+2,BSunit).i > i1) /* forward alignment */
	      { len = arr(oldA,i+2,BSunit).i - i1 + 1 ;
		while (len > 0 && /* skip initial padded region */
		       (index[i1] == index[i1-1] ||
			index2[i2] == index2[i2-1]))
		  { --len ; ++i1 ; ++i2 ; }
		while (len > 0)
		  { if (j && 
			(index[i1] == arr(newA,j-3,BSunit).i + 1) &&
			(index2[i2] == arr(newA,j-1,BSunit).i + 1))
		      j -= 5 ;
		    else
		      { array(newA,j+4,BSunit).i = 0 ; /* make space */
			arr(newA,j,BSunit).k = key ;
			arr(newA,j+1,BSunit).i = index[i1] ;
			arr(newA,j+3,BSunit).i = index2[i2] ;
		      }
				/* find aligned segment */
		    while (len > 0 && 
			   (index[i1+1]-index[i1] == index2[i2+1]-index2[i2]))
		      { --len ; ++i1 ; ++i2 ; }
		    arr(newA,j+2,BSunit).i = index[i1] ;
		    arr(newA,j+4,BSunit).i = index2[i2] ;
		    j += 5 ;
				/* skip unaligned segment */
		    do { --len ; ++i1 ; ++i2 ; } 
		    while (len > 0 && 
			   (index[i1] == index[i1-1] ||
			    index2[i2] == index2[i2-1])) ;
		  }
	      }
	    else		/* reverse alignment */
	      { len = i1 - arr(oldA,i+2,BSunit).i + 1 ;
		while (len > 0 && /* skip initial padded region */
		       (index[i1] == index[i1-1] ||
			index2[i2] == index2[i2-1]))
		  { --len ; --i1 ; ++i2 ; }
		while (len > 0)
		  { if (j && 
			(index[i1] == arr(newA,j-3,BSunit).i - 1) &&
			(index2[i2] == arr(newA,j-1,BSunit).i + 1))
		      j -= 5 ;
		    else
		      { array(newA,j+4,BSunit).i = 0 ; /* make space */
			arr(newA,j,BSunit).k = key ;
			arr(newA,j+1,BSunit).i = index[i1] ;
			arr(newA,j+3,BSunit).i = index2[i2] ;
		      }
				/* find aligned segment */
		    while (len > 0 && 
			   (index[i1-1]-index[i1-2] == index2[i2+1]-index2[i2]))
		      { --len ; --i1 ; ++i2 ; }
		    arr(newA,j+2,BSunit).i = index[i1] ;
		    arr(newA,j+4,BSunit).i = index2[i2] ;
		    j += 5 ;
				/* skip unaligned segment */
		    do { --len ; --i1 ; ++i2 ; }
		    while (len > 0 && 
			   (index[i1] == index[i1-1] ||
			    index2[i2] == index2[i2-1])) ;
		  }
	      }
	    if (index2 != identity)
	      messfree (index2) ;
	  }
        if (defragment)  /* needed in acembly, cDna not yet supported */
	  {  
	    KEY old = 0 ; int a1, b1, a2, b2 , j1, j2 ;
	    a1 = a2 = b1 = b2 = 0 ;
	    arrayDestroy (oldA) ; oldA = newA ;
	    newA = arrayCreate (arrayMax (oldA), BSunit) ;
	    for (key = 0, i = 0, j = 0; i < arrayMax(oldA) ; i += 5)
	      { 
		key = arr(oldA,i,BSunit).k ;
		i1 = arr(oldA,i+1,BSunit).i ;
		i2 = arr(oldA,i+3,BSunit).i ;
		j1 = arr(oldA,i+2,BSunit).i ;
		j2 = arr(oldA,i+4,BSunit).i ;
		if (key != old)
		  { 
		    if (old)
		      {
			array(newA,j+4,BSunit).i = 0 ; /* make space */
			arr(newA,j,BSunit).k = old ;
			arr(newA,j+1,BSunit).i = a1 ;
			arr(newA,j+3,BSunit).i = a2 ;
			arr(newA,j+2,BSunit).i = b1 ;
			arr(newA,j+4,BSunit).i = b2 ;
			j+= 5 ;
		      }
		    a1 = i1 ; a2 = i2 ; b1 = j1 ; b2 = j2 ; 
		  }
		else
		  {
		    b1 = j1 ; b2 = j2 ; 
		  }
		old = key ;
	      }
	    if (key) /* register the last one */
		  { 
		    array(newA,j+4,BSunit).i = 0 ; /* make space */
		    arr(newA,j,BSunit).k = old ;
		    arr(newA,j+1,BSunit).i = a1 ;
		    arr(newA,j+3,BSunit).i = a2 ;
		    arr(newA,j+2,BSunit).i = b1 ;
		    arr(newA,j+4,BSunit).i = b2 ;
		    a1 = i1 ; a2 = i2 ; b1 = j1 ; b2 = j2 ; 
		  }
	  }
	bsAddArray (obj, _Assembled_from, newA, 5) ;
      }
    oldA = arrayReCreate (oldA, 8192, BSunit) ;
    newA = arrayReCreate (newA,32768, BSunit) ;
				/* and Align_to_SCF */
    if (bsGetArray (obj, _Align_to_SCF, oldA, 4))
      { for (i = 0, j = 0 ; i < arrayMax(oldA) ; i += 4)
	  { i1 = arr(oldA,i,BSunit).i ;
	    i2 = arr(oldA,i+2,BSunit).i ;
	    len = arr(oldA,i+1,BSunit).i - i1 + 1 ;
	    while (len > 0 && (index[i1] == index[i1-1]))
	      { --len ; ++i1 ; ++i2 ; }	/* skip initial padded region */
	    while (len > 0)
	      { if (j && 
		    (index[i1] == arr(newA,j-2,BSunit).i + 1) &&
		    (i2 == arr(newA,j-1,BSunit).i + 1))
		  j -= 4 ;
		else
		  { array(newA,j+3,BSunit).i = 0 ; /* make space */
		    arr(newA,j,BSunit).i = index[i1] ;
		    arr(newA,j+2,BSunit).i = i2 ;
		  }
				/* find aligned segment */
		while (len && (index[i1+1]-index[i1] == 1))
		  { --len ; ++i1 ; ++i2 ; }
		arr(newA,j+1,BSunit).i = index[i1] ;
		arr(newA,j+3,BSunit).i = i2 ;
		j += 4 ;
				/* skip unaligned segment */
		do { --len ; ++i1 ; ++i2 ; }
		while (len > 0 && (index[i1] == index[i1-1])) ;
	      }
	  }
	bsAddArray (obj, _Align_to_SCF, newA, 4) ;
      }

    arrayDestroy (newA) ;
    arrayDestroy (oldA) ;
  }

 do_children:
				/* process subsequences */
				/* so can depad an assembly */
  if (bsGetKey (obj, _Subsequence, &key)) do
    { int *index2 = depadSequence (key, defragment) ;	/* RECURSION */
      messfree (index2) ;
    } while (bsGetKey (obj, _bsDown, &key)) ;

				/* done - now mark as unpadded */
  if (bsFindTag (obj, _Padded))
    bsAddTag (obj, _Unpadded) ; 

 done:
  arrayDestroy (s) ;
  bsSave (obj) ;
  return index ;
}
#endif

static void fixPhrapClips (KEY contig)
{ 
  int x1, x2, y1, y2, c1, c2, i ;
  KEY key ;
  OBJ Contig, Read ;
  Array units ;


  Contig = bsUpdate (contig) ;  
  if (!Contig)
    return ;

  units = arrayCreate(50, BSunit) ;
  bsGetArray (Contig, _Assembled_from, units, 5) ;
  bsFindTag (Contig, _Assembled_from) ;
  bsRemove (Contig) ;

  for (i = 0 ; i < arrayMax(units) ; i += 5)
    { key = arr (units, i, BSunit).k ;
      x1 = arr (units, i + 1, BSunit).i ;
      x2 = arr (units, i + 2, BSunit).i ;
      y1 = arr (units, i + 3, BSunit).i ;
      y2 = arr (units, i + 4, BSunit).i ;
      if ((Read = bsCreate (key)))
	{ 
	  if (bsGetData (Read, _Clipping, _Int, &c1) &&
	      bsGetData (Read, _bsRight, _Int, &c2))
	    { 
	      if (x1 < x2) 
		{ 
		  if (c1 > y1) { x1 += c1 - y1 ; y1 = c1 ;}
		  if (c2 < y2) { x2 += c2 - y2 ; y2 = c2 ;}
		}
	      else
		{ 
		  if (c1 > y1) { x1 -= c1 - y1 ; y1 = c1 ;}
		  if (c2 < y2) { x2 -= c2 - y2 ; y2 = c2 ;}
		}
	    }
	  bsDestroy (Read) ;
	}
      arr (units, i + 1, BSunit).i = x1 ;
      arr (units, i + 2, BSunit).i = x2 ;
      arr (units, i + 3, BSunit).i = y1 ;
      arr (units, i + 4, BSunit).i = y2 ;
    }

  bsAddArray (Contig, _Assembled_from, units, 5) ;
  bsSave (Contig) ;
  arrayDestroy (units) ;
}


BOOL depad (char *name) /* , BOOL defragment) */
{ KEYSET contigs ;
  KEY key, contig ;
  int i ;

  initialise () ;

  if (!name)
    { KEYSET ks = 
	query (0, "FIND Clone ; Main_clone ; FOLLOW Assembly ; Padded") ;
      if (!keySetMax(ks))
	{ keySetDestroy (ks) ; return FALSE ; }
      if (keySetMax(ks) > 1)
	{ messerror ("%d padded assemblies, aborting depad", 
		     keySetMax(ks)) ;
	  keySetDestroy (ks) ;
	  return FALSE ;
	}
      key = keySet(ks, 0) ;
      keySetDestroy (ks) ;
    }
  else if (!lexword2key (name, &key, _VSequence))
    { messout ("Can't find sequence %s", name) ;
      return FALSE ;
    }

  contigs = queryKey (key, ">Subsequence") ;
  i = keySetMax(contigs) ;
  while (i--)
  { contig = keySet(contigs,i) ;
    fixPhrapClips (contig) ;
   /* index = depadSequence (contig, defragment) ;
    messfree (index) ;
    */ 
  }
#ifdef ACEMBLY  
  alignToolsAdjustLink (key, contigs, 0);
#endif /* ACEMBLY */
  keySetDestroy (contigs) ;
  return TRUE ;
}

/*******************************************************************/
/********************** now for the padding ************************/

void padRead (KEY seq, int from, int to, int *index)
{ 
  OBJ obj = 0 ;
  Array s ;

  if (!obj)
    return ;

  if ((s = dnaGet(seq)))
    fixDNA (seq, s, index) ;
  obj = bsUpdate (seq) ;
  fixArray (obj, _Quality, index) ;
  fixArray (obj, _SCF_Position, index) ;
  bsCoordIndex (obj, index) ;

		/* Align_to_SCF */
  
  bsAddTag (obj, _Padded) ;

  bsSave (obj) ;
}

typedef struct { BSunit key, c1, c2, r1, r2 ; } AF ;

static int afOrder (const void* x, const void *y)
{ const AF *a = (const AF*) x, *b = (const AF*) y ;
  if (a->key.k > b->key.k) return 1 ;
  if (a->key.k < b->key.k) return -1 ;
  if (a->r1.i > b->r1.i) return 1 ;
  if (a->r1.i < b->r1.i) return -1 ;
  return 0 ;
}

void padContig (KEY seq)
{ 
  OBJ obj = bsCreate (seq) ;
  KEY oldKey = 0 ;
  register int i, j ;
  int c1, oldR2, gap, n ;
  int *ctgPads = 0, *index = 0 ;
  Array s = dnaGet (seq) ;
  Array oldA = 0 ;
  Array afA = 0 ;
  AF *af ;
  
  if (!s || !obj || bsFindTag (obj, _Padded))
    goto done ;

  /* Two pass strategy.
     First pass: make ctgPads.  
     For now rely on good ordering of _Assembled_from lines -- 
     should sort, but that breaks Array encapsulation
  */

  n = 1 + arrayMax(s) ;
  ctgPads = (int*) messalloc (n * sizeof(int)) ;
			/* rely on clearing to 0 by messalloc */
  oldA = arrayCreate (8192, BSunit) ; oldR2 = 0 ;
  if (bsGetArray (obj, _Assembled_from, oldA, 5))
    { afA = arrayCreate (arrayMax (oldA)/5, AF*) ;
      for (i = 0, j = 0 ; i < arrayMax(oldA) ; i += 5)
	array (afA,j++,AF*) = (AF*) arrp (oldA, i, BSunit) ;
      arraySort (afA, afOrder) ;
      for (i = 0, j = 0 ; i < arrayMax(afA) ; ++i)
	{ af = arr(afA,i,AF*) ;
	  c1 = af->c1.i ; if (af->c2.i < c1) ++c1 ; /* reverse */
	  gap = af->r1.i - oldR2 - 1 ;
	  if (af->key.k == oldKey && 
	      gap > ctgPads[c1]) /* pads needed */
	    ctgPads[c1] = gap ;
	  oldKey = af->key.k ; 
	  oldR2 = af->r2.i ;
	}
    }

  /* Now make index, and if nontrivial pad DNA, Quality and Position, 
     and fix internal coordinates
  */

  index = (int*) messalloc (n * sizeof(int)) ;
  for (i = 1, j = 1 ; i < n ; ++i, ++j)
    { j += ctgPads[i] ;
      index[i] = j ;
    }
  
  if (j > n)			/* pads added */
    {				/* pad sequence and store back */
      bsDestroy (obj) ;
      fixDNA (seq, s, index) ; s = 0 ;
      obj = bsUpdate (seq) ;
      fixArray (obj, _Quality, index) ;
      fixArray (obj, _SCF_Position, index) ;
      bsCoordIndex (obj, index) ;
    }

  /* Second pass: build index for each read and pad that, make new
     Assembled_from lines.
  */

  if (arrayMax(oldA)) for (i = 0 ; i < arrayMax(afA) ; ++i)
    { messcrash ("code not yet written, rd owes jean a beer") ;
    /* index2 = arrayReCreate (index2, 1000, int) ; */
      
    }

 done:
  arrayDestroy (s) ;
  arrayDestroy (afA) ; 
  arrayDestroy (oldA) ;
  bsSave (obj) ;
  return ;
}

BOOL pad (char *name)
{ 
  OBJ obj ;
  KEY key ;

  initialise () ;

  if (!lexword2key (name, &key, _VSequence))
    return FALSE ;

  if (!(obj = bsCreate (key)))
    return FALSE ;

  if (bsFindTag (obj, _Padded))
    { bsDestroy (obj) ;
      return TRUE ;
    }

  depad (0) ;		/* clear any currently padded assembly */

  bsDestroy (obj) ;
  return TRUE ;
}

/**************************** end of file ***************************/
/********************************************************************/

