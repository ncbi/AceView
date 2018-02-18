
#include "acedb.h"
#include "a.h"
#include "bs.h"
#include "orepeats.h"
#include "freeout.h"
#include "../whooks/classes.h"
#include "../whooks/systags.h"
#include "lex.h"

/*********************************************************************/

void oligoRepeatDoDestroy (OREPEAT *or)
{
  if (or)
    {
      arrayDestroy (or->oLengths) ;
      arrayDestroy (or->profile) ;
    }
}

/*********************************************************************/

void oligoRepeatRC (OREPEAT *or, int max)
{
  unsigned char cc, *cp, *cq ;

  if (or)
    {
      cc = array (or->profile, max - 1, unsigned char) ; /* make room */
      cp = arrp (or->profile, 0, unsigned char) ;
      cq = arrp (or->profile, arrayMax(or->profile) - 1, unsigned char) ;
      while (cp < cq)
	{
	  cc = *cp ; *cp = *cq ; *cq = cc ;
	  cp++ ; cq-- ;
	}
      or->phase = or->phase ; /* i think this is correct */
    }
}

/*********************************************************************/

OREPEAT *oligoRepeatRegister (OREPEAT *oligoRepeats, int x1, int x2, KEY key, AC_HANDLE h ) 
{
  OREPEAT *or = oligoRepeats ;
  Array a = arrayGet (key, unsigned char, "c") ;
  int i, n, n1, n2, x10, phase ;
  unsigned char *bc ;
  BOOL isUp = x1 < x2 ? FALSE : TRUE ;
  BOOL isOld = TRUE ;

  arrayDestroy (a) ;

  if (!a || arrayMax(a) < 3)
    goto abort ;

  n = arrayMax (a) ; n1 = 0 ;
  bc = arrp (a, 0, unsigned char) ;
  while (*bc && n > 0)
    { n1++ ; n-- ; bc++ ;}
  n2 = n - 1 ;
  x10 = x1 ;
  if (x1 < 0) { n2 += x1 ; x1 = 0 ; }
  if (n2 < 0)
    goto abort ;

  if (!or)
    {
      isOld = FALSE ;
      or = (OREPEAT *) halloc (sizeof (OREPEAT), h) ;
    }
  
 
  if  (!arrayExists (or->profile))
    or->profile = arrayHandleCreate (x1 + n2, unsigned char, h) ;
  array (or->profile, x1 + n2 - 1, unsigned char) = 0  ; /* make room */

  if (x10 > 0) x10 = 0 ;
  memcpy (arrp (or->profile, x1, unsigned char), arrp (a, n1 + 1 + x10, unsigned char), n2 - x10) ;
  for (i = 0 ; i < n1 ; i++)
    {
      phase = ((x10 - x1 + i % n1) + n1) % n1 ; 
      if (!isOld || phase == or->phase)
	break ;
    }
  if (isOld  && phase != or->phase)
    messerror ("inconsistent phases in oligoRepeatRegister") ;
  or->phase = phase ;
  if (isUp)
    {

    }
 
  arrayDestroy (or->oLengths) ;  /* should be merged */
  or->oLengths = arrayHandleCreate (n, int, h) ;
  for (i = 0 ; i < n1 ; i++)
    array (or->oLengths, i, int) = arr (a, i, unsigned char) ;
  
 abort:
  arrayDestroy (a);
  return or ;
}

/*********************************************************************/

BOOL oligoRepeatDump (FILE* f, Stack s, KEY k) 
{
 int level = 0 ;
  Array a ;
  int n, x, i ; 
  unsigned char *bc ;
  int ntypes = 0 ;

  a = arrayGet (k, unsigned char, "c") ;
  if (!a || !arrayMax(a))
    { arrayDestroy (a) ;
      return FALSE ;
    }

  if (f)
    level = freeOutSetFile (f) ;
  else if (s)
    level = freeOutSetStack (s) ;

  n = arrayMax (a) ;
  bc = arrp (a, 0, unsigned char) ;
  while (*bc && n > 0)
    { ntypes++ ; n-- ; bc++ ;}
  freeOutf ("%d ", ntypes) ;
  freeOutf ("%d \n", n - ntypes - 1) ;
  for (i = 0 ; i < ntypes ; i++)
    freeOutf ("%d ", arr (a, i, unsigned char)) ;
  freeOutf ("\n") ;
  while (n)
    { for (i = 0 ; i < 50 && n ; i++, --n)
	{ x = *bc++ ;
	  freeOutf ("%d ", x) ;
	}
      freeOut ("\n") ;
    }
  freeOut ("\n") ;

  arrayDestroy (a) ;
  if (level)
    freeOutClose (level) ;
  return TRUE ;
}

/*********************************************************************/

BOOL oligoRepeatParse (int level, KEY key) 
{
  Array a = arrayCreate (1000, unsigned char) ;
  int x, n = 0, n1 = 0, n2 = 0, line = 0  ;
  char *cp = 0 ;

  if (class(key) != _VOligoRepeat)
    messcrash ("oligorepeatParse called on a non-Oligorepeat key") ;

  while (freecard (level) && freeint (&x))
    {
      line++ ;
      do 
	{ 
	  if (x < 0)
	    {
	      messerror ("Junk at line %d while rerdeaing OligoRepeat %s: %s",
			 freestreamline(level), name(key), cp ? cp : "") ;
	      goto abort ;
	    }
	  if (n == 0) n1 = x ; /* number of oLength */
	  else if (n == 1) n2 = x ; /* number of values */
	  else
	    {
	      if (x > 255)
		x = 255 ;
	      array(a, n - 2, unsigned char) = x ; 
	    }
	  n++ ;
	} while (freeint (&x)) ;
      if (line == 2)  array(a, n++ - 2, unsigned char) = 0; 
      if ((cp = freeword()) )
	{
	  messerror ("Junk at line %d while rerdeaing OligoRepeat %s: %s",
		     freestreamline(level), name(key), cp) ;
	  goto abort ;
	}
    }
  if ((cp = freeword()))
    {
      messerror ("Junk at line %d while rerdeaing OligoRepeat %s: %s",
		 freestreamline(level), name(key), cp) ;
      goto abort2 ;
    }

  if (n < 3)
    {
      messerror ("Error parsing Oligorepeat %s at line %d (n != nlengths + npoints + 1)", 
		 name(key), freestreamline(level)) ;
      goto abort2 ;
    }

  if (n != n1 + n2 + 3)
    {
      messerror ("Error parsing Oligorepeat %s at line %d (n != nlengths + npoints + 1)", 
	     name(key), freestreamline(level)) ;
      goto abort2 ;
    }
  printf ("orepeat read %s  %d\n", name(key), n2) ;
  if (FALSE && arr (a, n1, unsigned char) != 0)
    { 
      messerror ("Error parsing Oligorepeat %s at line %d (nlengths list not zero terminated)", 
		 name(key), freestreamline(level)) ;
      goto abort2 ;
    }

  if (arrayMax(a))
    { KEY seq ;
      OBJ obj ;

      lexaddkey (name(key), &seq, _VSequence) ;
      if ((obj = bsUpdate (seq)))
	{ 
	  bsAddKey (obj, str2tag("OligoRepeat"), key) ;
	  bsAddData (obj, _bsRight, _Int, &n) ;
	  bsSave (obj) ;
	}

      arrayStore (key, a, "c") ;
    }

  arrayDestroy (a) ;
  return TRUE ;

 abort:
  messerror ("Error parsing Oligorepeat %s at line %d (not an int 0-255)", 
	     name(key), freestreamline(level)) ;
 abort2:
  arrayDestroy (a) ;
  return FALSE ;
}

/*********************************************************************/
/*********************************************************************/
