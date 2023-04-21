#include "ac.h"

typedef enum  {ZERO=0, NUM, PLUS, MULT, E, F, H, K } OP ;
 
typedef struct polynomeStruct PP ;
typedef struct polynomeStruct { 
  OP op ;
  PP *p1, *p2 ;
  int n ;
  double z ;
} ;

typedef struct cartanStruct {
  int rank ;
  int *Cartan ;
  int *HW ;
} CC ;

static PP* ppNew (OP op, int n,  AC_HANDLE h) ;
static PP* ppMult (PP *p1, PP *p2,  AC_HANDLE h) ;
static PP* ppPlus (PP *p1, PP *p2,  AC_HANDLE h) ;
static BOOL ppEqual (PP *p1, PP *p2) ;

/***********************************************************************/

static CC *newCartan (char *algebra, int rank0, int rank1, AC_HANDLE h)
{
  int i, rank, rk ;
  CC *cc = halloc (sizeof (CC), h) ;
  rank = rank0 + rank1 ;

  if (! strcmp (algebra, "A"))
    {
      rank1 = 0 ; rank = rank0 + rank1 ; rk = rank  ;
      cc->Cartan = halloc (rk * rk * sizeof(int), h) ;
      cc->HW = halloc (rk * sizeof(int), h) ;
      cc->rank = rank ;
      for (i = 0 ; i < rank ; i++)
	{
	  cc->Cartan[rk*i + i] = 2 ;
	  if (i > 0) cc->Cartan[rk*i + i - 1] = -1 ;
	  if (i < rank - 1) cc->Cartan[rk*i + i + 1] = -1 ;
	}
    }
  else if (! strcmp (algebra, "Z"))
    {

    }
  else
    messcrash ("Unknown algebra class %s, i only know classes: A", "algebra") ;
  return cc ;
} /* newCartan */

/***********************************************************************/

static void showPp (CC *cc, PP *pp, int level)
{
  BOOL debug = FALSE ;

  switch (pp->op)
    {
    case ZERO:
      fprintf (stderr, " 0 ") ;
      if (debug) fprintf (stderr, "[%p]  ",pp) ;
      break ;
    case NUM:
      fprintf (stderr, " %d ", pp->n) ;
      if (debug) fprintf (stderr, "[%p]  ",pp) ;
      break ;
    case E:
      fprintf (stderr, " E%d ", pp->n) ;
      if (debug) fprintf (stderr, "[%p]  ",pp) ;
      break ;
    case F:
      fprintf (stderr, " F%d ", pp->n) ;
      if (debug) fprintf (stderr, "[%p]  ",pp) ;
      break ;
    case H:
      fprintf (stderr, " H%d ", pp->n) ;
      if (debug) fprintf (stderr, "[%p]  ",pp) ;
      break ;
     case K:
       {
	 int i ;
	 fprintf (stderr, " |K") ;
	 for (i = 0 ; i < cc->rank ; i++)
	   fprintf (stderr, " %d", cc->HW[i]) ;
	 fprintf (stderr, "> ") ;
       }
      break ;
    case PLUS:
      showPp (cc ,pp->p1, level + 1) ;
      fprintf (stderr, " + ") ;
      if (debug) fprintf (stderr, "[%p]  ",pp) ;
      showPp (cc ,pp->p2, level + 1) ;
      break ;
    case MULT:
      if (pp->p1 && pp->p1->op == PLUS)
	fprintf (stderr, "(") ;
      showPp (cc ,pp->p1, level + 1) ;
      if (pp->p1  && pp->p1->op == PLUS)
	fprintf (stderr, ")") ;
      fprintf (stderr, " * ") ;
      if (debug) fprintf (stderr, "[%p]  ",pp) ;
      if (pp->p2 && pp->p2->op == PLUS) 
	fprintf (stderr, "(") ;
      showPp (cc ,pp->p2, level + 1) ;
      if (pp->p2  && pp->p2->op == PLUS)  
	fprintf (stderr, ")") ;
      break ;
    }
  if (level <= 0) 
    fprintf (stderr, "\n") ;
} /* showPp */

/***********************************************************************/
/* perform a * (b + c) -> a*b + a*c
 * retur TRUE if something happened
 * false if no expansion was performed
 */
static BOOL ppExpand (CC *cc, PP *pp, AC_HANDLE h)
{
  BOOL modif = FALSE ;
  PP *p1 = pp ? pp->p1 : 0 ;
  PP *p2 = pp ? pp->p2 : 0 ;
  BOOL debug = FALSE ;

  if (! pp || pp->op == NUM || pp->op == ZERO || pp->op == E || pp->op == F || pp->op == H || pp->op == K) return FALSE ;
  ppExpand (cc, p1, h) ;
  ppExpand (cc, p2, h) ;

  modif = TRUE ;
  while (modif)
    { 
      if (! pp || pp->op == NUM || pp->op == ZERO || pp->op == E || pp->op == F || pp->op == H || pp->op == K)
	break ;
      p1 = pp->p1 ;
      p2 = pp->p2 ;
      modif = FALSE ;
      if (! p1 && ! p2)
	{ pp->op = ZERO ; pp->n = pp->z = 0 ; modif = TRUE ; continue ; }
      else if (! p1)
	{
	  *pp = *p2 ;  
	  if (debug) fprintf (stderr, "--- ppExpand  : free p2 %p\n", p2) ;
	  if (0) ac_free (p2) ;
	  modif = TRUE ; continue ;
	}
      else if (! p2)
	{
	  *pp = *p1 ;
	  if (debug) fprintf (stderr, "--- ppExpand  : free p1 %p\n", p1) ;
	  if (0) ac_free (p1) ; 
	  modif = TRUE ; continue ;
	}
      else if (p1->op  == MULT && p2->op == K)
	{
	  PP *p11 = p1->p1 ;
	  PP *p12 = p1->p2 ;
	  
	  pp->p1 = p11 ;
	  pp->p2 = ppMult (p12, p2, h) ;
	  ppExpand (cc, pp->p2, 0) ;
	  modif = TRUE ; continue ;
	} 
     else if (p1->op  == pp->op && (pp->op == PLUS || pp->op == MULT))
	{
	  PP *p11 = p1->p1 ;
	  PP *p12 = p1->p2 ;
	  
	  pp->p1 = p11 ;
	  pp->p2 = p1 ;
	  p1->p1 = p12 ;
	  p1->p2 = p2 ;
	  modif = TRUE ; continue ;
	} 
      else if (p1->op == NUM && p2->op == NUM)
	{
	  if (pp->op == PLUS)
	    { 
	      pp->n = p1->n + p2->n ;
	      pp->z = p1->z + p2->z ;
	    }
	  else if (pp->op == MULT)
	    { 
	      pp->n = p1->n * p2->n ;
	      pp->z = p1->z * p2->z ;
	    } 
	  pp->op = NUM ; modif = FALSE ; 
	  if (pp->z == 0) { pp->op = ZERO ; modif = TRUE ; }
	  if (debug) fprintf (stderr, "--- ppExpand  : free p1/p2 %p/%p\n", p1, p2) ;
	  continue ;
	}
      else if (p2->op == NUM) /* numbers are abelian and move p1 */
	{
	  pp->p1 = p2 ; pp->p2 = p1 ;
	  p1 =  pp->p1 ; p2 =  pp->p2 ;
	  modif = TRUE ; continue ;
	}
      else if (pp->op == MULT && (p1->op == ZERO ||  p2->op == ZERO))
	{
	  pp->op = ZERO ; pp->n = pp->z = 0 ;
	  p1 =  pp->p1 = 0 ; p2 =  pp->p2 = 0 ;
	  modif = FALSE ; continue ;
	}
      else if (pp->op == PLUS && p1->op == ZERO && p2->op == ZERO)
	{
	  pp->op = ZERO ; pp->n = pp->z = 0 ;
	  p1 =  pp->p1 = 0 ; p2 =  pp->p2 = 0 ;
	  modif = FALSE ; continue ;
	}
       else if (pp->op == PLUS && p1->op == ZERO)
	{
	  *pp = *p2 ;
	  modif = TRUE ; continue ;
	}
      else if (pp->op == PLUS && p2->op == ZERO)
	{
	  *pp = *p1 ;
	  modif = TRUE ; continue ;
	}
      else if (pp->op == PLUS && p2 == 0)
	{
	  *pp = *p1 ;
	  modif = TRUE ; continue ;
	}
      else if (1 && pp->op == PLUS && p1 && p2 && p1->op == MULT && p2->op == MULT &&
		p1->p1 && p1->p1->op == NUM &&
		p2->p1 && p2->p1->op == NUM &&
	       p1->p2 && 
	       ppEqual (p1->p2, p2->p2)
	       )
	{
	  *pp = *p1 ; p1->p1->n +=  p2->p1->n ; p1->p1->z +=  p2->p1->z ;
	  if ( p1->p1->z == 0) p1->p1->op = ZERO ;
	  modif = TRUE ; continue ;
	}
      else if (pp->op == PLUS && p2->op == PLUS && p2->p1->op == ZERO)
	{
	  pp->p2 = p2->p2 ;
	  modif = TRUE ; continue ;
	}
      else if (p1->op == NUM) 
	{
	  if (pp->op == PLUS && p2->op == PLUS && p2->p1 && p2->p1->op == NUM)
	    { 
	      pp->p2 = p2->p2 ;
	      pp->p1 = p2->p1 ;
	      p2->p1->n += p1->n ;
	      p2->p1->z += p1->z ;
	      if (p2->p1->z == 0) p2->p1->op = ZERO ;
	      if (debug) fprintf (stderr, "--- ppExpand  : num free p1 %p\n", p1) ;
	      if (0) ac_free (p1) ;  
	      p1 =  pp->p1 ; p2 =  pp->p2 ;
	      modif = TRUE ; continue ;
	    }
	  else if (pp->op == MULT && p2->op == MULT && p2->p1 && p2->p1->op == NUM)
	    { 
	      pp->p2 = p2->p2 ;
	      pp->p1 = p2->p1 ;
	      p2->p1->n *= p1->n ;
	      p2->p1->z *= p1->z ;
	      if (p2->p1->z == 0) p2->p1->op = ZERO ;
	      if (debug) fprintf (stderr, "--- ppExpand  : mult free p1 %p\n", p1) ;
	      if (0) ac_free (p1) ;  
	      p1 =  pp->p1 ; p2 =  pp->p2 ;
	      modif = TRUE ; continue ;
	    }
	}
   else if (pp->op == MULT && p1 && p2 && p2->op == PLUS)
      {
	pp->op = PLUS ;
	pp->p1 = ppMult (p1, p2->p1, h) ;
	pp->p2 = ppMult (p1, p2->p2, h) ;
	if (debug) fprintf (stderr, "--- ppExpand  : mult plus free p2 %p\n", p2) ;
	if (0) ac_free (p2) ;
	if(debug) showPp (cc ,pp, 0) ;
	ppExpand (cc, pp->p1, h) ;
	if(debug) showPp (cc ,pp, 0) ;
	ppExpand (cc, pp->p2, h) ;
	if(debug) showPp (cc ,pp, 0) ;
	modif = TRUE ; continue ;
      } 
     else if (pp->op == MULT && p1 && p2 && p1->op == PLUS)
	{
	  pp->op = PLUS ;
	  if (debug) fprintf (stderr, "--- ppExpand  : mult plus p1=%p p1->p1=%p p1->p2=%p p2=%p\n", p1,p1->p1,p1->p2,p2) ;
	  pp->p1 = ppMult (p1->p1, p2, h) ;
	  if (debug) fprintf (stderr, "--- ppExpand  : mult plus p1=%p p1->p1=%p p1->p2=%p p2=%p\n", p1,p1->p1,p1->p2,p2) ;
	  pp->p2 = ppMult (p1->p2, p2, h) ;
	  if(debug) showPp (cc ,pp, 0) ;
	  ppExpand (cc, pp->p1, h) ;
	  if(debug) showPp (cc ,pp, 0) ;
	  ppExpand (cc, pp->p2, h) ;
	  if(debug) showPp (cc ,pp, 0) ;
	  if (debug) fprintf (stderr, "--- ppExpand  : mult plus free p1 %p\n", p1) ;
	  if (0) ac_free (p1) ;
	  modif = TRUE ; continue ;
	}
    else if (1 && pp->op == MULT && p1 && p1->op == E && p2->op == F)
      {
	if (p1->n == p2->n)
	  {
	    pp->op = PLUS ;
	    pp->p1 = ppNew (H, p1->n, h) ;
	    pp->p2 = ppMult (p2, p1, h) ;
	    modif = TRUE ; continue ;
	  }
	else
	  {
	    pp->p1 = p2 ;
	    pp->p2 = p1 ; 
	    modif = TRUE ; continue ;
	  }
      }
    else if (pp->op == MULT && p1 && p1->op != NUM && p2->op == MULT && p2->p1 && p2->p1->op == NUM)
      {
	
	    PP *p11 = pp->p1 ; 
	    PP *p22 = p2->p2 ; 
	    PP *p21 = ppNew (ZERO, 0, h) ;

	    *p21 = *p2->p1 ; p21->p2 = 0 ;
	    pp->p1 = p21 ;
	    pp->p2 = ppMult (p11, p22, h) ;
	    ppExpand (cc, pp->p2, 0) ;
	    modif = TRUE ; continue ;

      }
      else if (pp->op == MULT && p1 && p1->op == E && p2->op == MULT && p2->p1 && p2->p1->op == F)
      {
	if (p1->n == p2->p1->n)
	  {
	    PP *p11 = pp->p1 ; 
	    PP *p22 = p2->p2 ; 
	    PP *p21 = ppNew (ZERO, 0, h) ;

	    *p21 = *p2->p1 ; p21->p2 = 0 ;
	    pp->op = PLUS ;
	    pp->p1 = ppMult (ppNew (H, p1->n, h), p22, h) ;
	    pp->p2 = ppMult (p21, ppMult (p11, p22, h), h) ;
	    ppExpand (cc, pp->p1, 0) ;
	    ppExpand (cc, pp->p2, 0) ;
	    modif = TRUE ; continue ;
	  }
	else
	  {
	    PP *p11 = pp->p1 ; 
	    PP *p22 = p2->p2 ; 
	    PP *p21 = ppNew (ZERO, 0, h) ;

	    *p21 = *p2->p1 ; p21->p2 = 0 ;
	    pp->p1 = p21 ;
	    pp->p2 = ppMult (p11, p22, h) ;
	    ppExpand (cc, pp->p2, 0) ;
	    modif = TRUE ; continue ;
	  }
      }
      else if (pp->op == MULT && p1 && p1->op == H && p2->op == MULT && p2->p1 && p2->p1->op == F)
      {
	int x = cc->Cartan [cc->rank * (p1->n - 1) + p2->p1->n - 1] ; /* cartan (p1->n, p2->p1->n) ; */
	if (x)
	  {
	    PP *p11 = pp->p1 ; 
	    PP *p22 = p2->p2 ; 
	    PP *p21 = ppNew (ZERO, 0, h) ;

	    *p21 = *p2->p1 ; p21->p2 = 0 ;
	    pp->op = PLUS ;
	    pp->p1 = ppMult (ppNew (NUM, -x, h),  ppMult (p21,p22,h), h) ;
	    pp->p2 = ppMult (p21, ppMult (p11, p22, h), h) ;
	    ppExpand (cc, pp->p1, 0) ;
	    ppExpand (cc, pp->p2, 0) ;
	    modif = TRUE ; continue ;
	  }
	else
	  {
	    PP *p11 = pp->p1 ; 
	    PP *p22 = p2->p2 ; 
	    PP *p21 = ppNew (ZERO, 0, h) ;

	    *p21 = *p2->p1 ; p21->p2 = 0 ;
	    pp->p1 = p21 ;
	    pp->p2 = ppMult (p11, p22, h) ;
	    ppExpand (cc, pp->p2, 0) ;
	    modif = TRUE ; continue ;
	  }
      }

     else if (1 && pp->op == MULT && p1 && p1->op == E && p2->op == K)
      {
	pp->op = ZERO ; pp->p1 = pp->p2 = 0 ; pp->n = pp->z = 0 ;
	modif = FALSE ; continue ;
      }
     else if (1 && pp->op == MULT && p1 && p1->op == H && p2->op == K)
      {
	p1->op = NUM ;
	p1->n = p1->z = cc->HW[p1->n - 1] ;
	modif = TRUE ; continue ;
      }
    }
  
  return modif ;
} /* ppExpand */

/***********************************************************************/

static PP* ppNew (OP op, int n,  AC_HANDLE h)
{
  PP *pp = (PP *) halloc (sizeof (PP), h) ;
  BOOL debug = FALSE ;

  if (debug) fprintf (stderr, "--- ppNew : alloc pp %p\n", pp) ;
  pp->op = op ;
  pp->n = pp->z = n ;
  return pp ;
} /* ppNew */

/***********************************************************************/

static PP* ppCopy (PP *p,  AC_HANDLE h)
{
  if (p)
    {
      PP *pp = (PP *) halloc (sizeof (PP), h) ;
      
      *pp = *p ;
      if (p->p1) pp->p1 = ppCopy (p->p1, h) ;
      if (p->p1) pp->p2 = ppCopy (p->p2, h) ;
      
      return pp ;
    }
  return 0 ;
} /* ppCopy */

/***********************************************************************/

static BOOL ppEqual (PP *p1, PP *p2)
{
  if (! p1 && !p2) return TRUE ;
  if (! p1 && p2) return FALSE ;
  if (! p2 && p1) return FALSE ;
  if (p1->op != p2->op) return FALSE ;
  if (p1->n != p2->n) return FALSE ;
  if (p1->z != p2->z) return FALSE ;
  if (! ppEqual (p1->p1, p2->p1)) return FALSE ;
  if (! ppEqual (p1->p2, p2->p2)) return FALSE ;
  return TRUE ;
} /* ppEqual */

/***********************************************************************/

static PP* ppPlus (PP *p1, PP *p2,  AC_HANDLE h)
{
    PP *p11, *pp = 0 ;
  BOOL debug =  FALSE ;

  if (debug) fprintf (stderr, "--- ppPlus : alloc pp %p, add %p to %p\n", pp, p1, p2) ;
  p1 = ppCopy (p1, h) ;
  p2 = ppCopy (p2, h) ;
  p11 = p1 ;
  if (p11->op == PLUS)
    {
      while (p11->p2->op == PLUS) p11 = p11->p2 ;
      if (! p11->p2) { p11->p2 = p2 ; return p1 ; }
      pp = (PP *) halloc (sizeof (PP), h) ;
      pp->p1 = p11->p2 ; pp->p2 = p2 ; p11->p2 = pp ; pp->op = PLUS ;
      return p1 ;
    }
  
  pp = (PP *) halloc (sizeof (PP), h) ;
  pp->op = PLUS ;
  pp->p1 = p1 ;
  pp->p2 = p2 ;
  return pp ;
} /* ppPlus */

/***********************************************************************/

static PP* ppMult (PP *p1, PP *p2,  AC_HANDLE h)
{
  PP *p11, *pp = 0 ;
  BOOL debug =  FALSE ;

  if (debug) fprintf (stderr, "--- ppMult : alloc pp %p, multiply %p by %p\n", pp, p1, p2) ;
  p1 = ppCopy (p1, h) ;
  p2 = ppCopy (p2, h) ;
  p11 = p1 ;
  if (p11->op == MULT)
    {
      while (p11->p2->op == MULT) p11 = p11->p2 ;
      if (! p11->p2) { p11->p2 = p2 ; return p1 ; }
      pp = (PP *) halloc (sizeof (PP), h) ;
      pp->p1 = p11->p2 ; pp->p2 = p2 ; p11->p2 = pp ; pp->op = MULT ;
      return p1 ;
    }
  
  pp = (PP *) halloc (sizeof (PP), h) ;
  pp->op = MULT ;
  pp->p1 = p1 ;
  pp->p2 = p2 ;
  return pp ;
} /* ppMult  */

/***********************************************************************/

static PP* ppMultiOp (CC *cc, OP op, const char *txt,  AC_HANDLE h)
{
  int i, n = strlen (txt) ;
  PP *pp = 0 ;
  const char *cp = txt ;

  if (n < 1)
    pp = ppNew (ZERO, 0, h) ;
  else if (n >= 1)
    {
      i = *cp - '0' ;
      if (i > cc->rank + 1)
	messcrash ("ppMulti %s invokes number %d > rank = %d", txt, i, cc->rank) ;
      pp = ppNew (op, i, h) ;
      for (cp = txt + 1 ; *cp ; cp++)
	{
	  i = *cp - '0' ;
	  if (i > cc->rank + 1)
	    messcrash ("ppMulti %s invokes number %d > rank = %d", txt, i, cc->rank) ;
	  pp = ppMult (pp, ppNew (op, i, h), h) ;
	}
    }

  return pp ;
} /* ppMultOp  */

/***********************************************************************/
/***********************************************************************/

#define VMAX 110
static double popEval (double *popm, double *popf, double *dpm, double *dpf, double *mm, double *fm, BOOL doShow, double *errmp, double *errfp) 
{
  int i, j ;
  double z, errm = 0, errf = 0 ;
  double vm[VMAX] ; /* vie moyenne restante at age x */
  double vf[VMAX] ; /* vie moyenne restante at age x */
  double vmf[VMAX] ; /* vie moyenne restante at age x */
  double popmf[VMAX] ; 

  popm[0] = 1000000 ;
  for (i = 1 ; i < VMAX ; i++)
    popm[i] = popm[i-1] * (1 - dpm[i-1]) ;
  popf[0] = 1000000 ;
  for (i = 1 ; i < VMAX ; i++)
    popf[i] = popf[i-1] * (1 - dpf[i-1]) ;
  popmf[0] = 1000000 ;
  for (i = 1 ; i < VMAX ; i++)
    popmf[i] = popmf[i-1] * (1 - dpf[i-1]) ;

  for (i = 0 ; i < VMAX ; i++)
    {  /* compute life expectancy at age i */
      z = 0 ;
      for (j = i+1 ; j < VMAX ; j++)
	z +=  popm[j] ;
      z /= popm[i] ;
      vm[i] = z ;
 
      z = 0 ;
      for (j = i+1 ; j < VMAX ; j++)
	z +=  popf[j] ;
      z /= popf[i] ;
      vf[i] = z ;
 
      z = 0 ;
      for (j = i+1 ; j < VMAX ; j++)
	z +=  1 - (1 - popf[j]/popf[i]) *  (1 - popm[j]/popm[i]) ;
      vmf[i] = z ;

      if (i >= 60 && i <= 99)
	errm += (vm[i] - mm[i]) *  (vm[i] - mm[i]) ;
      if (i >= 60 && i <= 99)
	errf += (vf[i] - fm[i]) *  (vf[i] - fm[i]) ;
    }
  if (doShow)
    {
      printf ("#Age\tPopulation Homme\tPopulation Femme\tProba homme\tProba femme\tEsperance de vie homme\tEsperance de vie femme\tObservation chez les hommes\tObservation chez les femmes\tAttendu pour un couple\n") ;
      for (i = 60 ; i <= 99 ; i++)
	printf ("%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", i, 1000000*popm[i]/popm[60], 1000000*popf[i]/popf[60], dpm[i], dpf[i],mm[i], fm[i], vmf[i]) ;
    }
  *errmp = errm ; *errfp = errf ;
  return errm ;
}

static void viager (void)
{
  /* */
  int i, n ; 
  double em, ef ;
  double z, popm[VMAX], popf[VMAX] ;
  double dpm[VMAX] ; /* probab de disparaitre dans l'annee */
  double dpf[VMAX] ; /* probab de disparaitre dans l'annee */
  double mm[VMAX] ; /* vie moyenne restante at age x */
  double fm[VMAX] ; /* vie moyenne restante at age x */

  for (i = 0 ; i < 60 ; i++) fm[i] = mm[i] = 20 ;
  fm[60] =  26.95 ; mm[60] = 23.74 ;
  fm[61] =  26.02 ; mm[61] = 22.83 ;
  fm[62] =  25.09 ; mm[62] = 21.93 ; 
  fm[63] =  24.16 ; mm[63] = 21.04 ; 
  fm[64] =  23.23 ; mm[64] = 20.15 ; 
  fm[65] =  22.30 ; mm[65] = 19.27 ; 
  fm[66] =  21.37 ; mm[66] = 18.39 ; 
  fm[67] =  20.45 ; mm[67] = 17.52 ; 
  fm[68] =  19.53 ; mm[68] = 16.67 ; 
  fm[69] =  18.62 ; mm[69] = 15.82 ; 

  fm[70] =  17.72 ; mm[70] = 14.99 ; 
  fm[71] =  16.83 ; mm[71] = 14.18 ; 
  fm[72] =  15.95 ; mm[72] = 13.39 ; 
  fm[73] =  15.08 ; mm[73] = 12.61 ; 
  fm[74] =  14.23 ; mm[74] = 11.86 ; 
  fm[75] =  13.39 ; mm[75] = 11.13 ; 
  fm[76] =  12.57 ; mm[76] = 10.42 ; 
  fm[77] =  11.77 ; mm[77] = 9.74 ; 
  fm[78] =  11.00 ; mm[78] = 9.08 ; 
  fm[79] =  10.26 ; mm[79] = 8.46 ; 

  fm[80] =  9.54 ; mm[80] = 7.86 ; 
  fm[81] =  8.86 ; mm[81] = 7.30 ; 
  fm[82] =  8.21 ; mm[82] = 6.76 ; 
  fm[83] =  7.60 ; mm[83] = 6.25 ; 
  fm[84] =  7.02 ; mm[84] = 5.77 ; 
  fm[85] =  6.48 ; mm[85] = 5.33 ; 
  fm[86] =  5.97 ; mm[86] = 4.92 ; 
  fm[87] =  5.51 ; mm[87] = 4.54 ; 
  fm[88] =  5.07 ; mm[88] = 4.18 ; 
  fm[89] =  4.68 ; mm[89] = 3.85 ; 

  fm[90] =  4.31 ; mm[90] = 3.55 ; 
  fm[91] =  3.98 ; mm[91] = 3.27 ; 
  fm[92] =  3.67 ; mm[92] = 3.03 ; 
  fm[93] =  3.39 ; mm[93] = 2.82 ; 
  fm[94] =  3.14 ; mm[94] = 2.62 ; 
  fm[95] =  2.90 ; mm[95] = 2.45 ; 
  fm[96] =  2.68 ; mm[96] = 2.28 ; 
  fm[97] =  2.48 ; mm[97] = 2.13 ; 
  fm[98] =  2.29 ; mm[98] = 1.99 ; 
  fm[99] =  2.12 ; mm[99] = 1.86 ; 





  z = .013 ;
  for (i = 0 ; i < VMAX ; i++)
    {
      dpm[i] = z ;
      dpf[i] = z ;
      if (i > 80) z = .05 ;
      if (i > 90) z = .2 ;
    }
  
  if (0) for (i = 0 ; i < VMAX ; i++)
    {
      z = .05 ;
      dpm[i] = z ;
      dpf[i] = z ;
    }
   
   for (n = 0 ; n < 1000 ; n++)
    {
      double dp2m[VMAX] ;
      double dp2f[VMAX] ;
      double em0, ef0 ;
      
      popEval (popm, popf, dpm, dpf, mm, fm, FALSE, &em0, &ef0) ;

      for (i = 60 ; i < VMAX ; i++)
	{
	  double am = dpm[i] ;
	  double af = dpf[i] ;
	  if (n < 20) { dpm[i] *= 1.1 ; dpf[i] *= 1.1 ; }
	  else { dpm[i] *= 1.01 ; dpf[i] *= 1.01 ; }
	  popEval (popm, popf, dpm, dpf, mm, fm, FALSE, &em, &ef) ;
	  if (em < em0) 
	    {
	      dp2m[i] = dpm[i] ;
	    }
	  else
	    {
	      if (n < 20)
		dp2m[i] = .90 * am ; 
	      else if (n < 50)
		dp2m[i] = .99 * am ; 
	      else if (n < 500)
		dp2m[i] = .99 * am ; 
	      else
		dp2m[i] = .9999 * am ; 
	    }
	  if (ef < ef0) 
	    {
	      dp2f[i] = dpf[i] ;
	    }
	  else
	    {
	      if (n < 20)
		dp2f[i] = .90 * af ; 
	      else if (n < 50)
		dp2f[i] = .99 * af ; 
	      else if (n < 500)
		dp2f[i] = .999 * af ; 
	      else
		dp2f[i] = .9999 * af ; 
	    }
	  dpm[i] = am ;
	  dpf[i] = af ;
	}
      for (i = 60 ; i < VMAX ; i++)
	{
	  dpm[i] = dp2m[i] ;
	  dpf[i] = dp2f[i] ;
	}
      fprintf (stderr, "%d\t%g\t%g\n", n, em0, ef0) ;
    }

  popEval (popm, popf, dpm, dpf, mm, fm, TRUE, &em, &ef) ;
}

/***********************************************************************/
/***********************************************************************/

int main (int argc, char *argv[])
{
  AC_HANDLE h = ac_new_handle () ;
  PP *p1, *p2, *p3, *p23, *p4, *p5, *pk ;
  CC *cc ;

  messErrorInit (argv[0]) ;

  if (0) 
    {
      viager () ;
      exit (0) ;
    }
  cc = newCartan ("A", 2, 0, h) ;
  cc->HW[0] = cc->HW[1] = 1 ; /* adjointe */
  
  p1 = ppPlus (ppNew (E, 1, h), ppNew (E, 2, h), h) ;
  p2 = ppNew (F, 1, h) ;
  p3 = ppNew (F, 2, h) ;
  p23 = ppPlus (p2, p3, h) ;
  p4 = ppMult (p1, p23, h) ;
  pk = ppNew (K, 3, h) ;
  p5 = ppMult (p4, pk, h) ;
  if (1)
    {
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }

  if (1)
    {
      p5 = ppMult (p1, pk, h) ;
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }

  if (1)
    {
      p5 = ppMult (ppNew (E, 1, h), ppNew (F, 1, h), h) ;
      p5 = ppMult (p5, pk, h) ;
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }

   if (1)
    {
      p1 = ppMult (ppNew (E, 1, h), ppNew (E, 1, h), h) ;
      p2 = ppMult (ppNew (F, 1, h), ppNew (F, 1, h), h) ;
      p5 = ppMult (p1, p2, h) ;
      p5 = ppMult (p5, pk, h) ;
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }
   if (1)
    {
      p1 = ppMult (ppMult (ppNew (E, 1, h), ppNew (E, 1, h), h), ppNew (E, 1, h), h) ;
      p2 =  ppMult (ppMult (ppNew (F, 1, h), ppNew (F, 1, h), h), ppNew (F, 1, h), h) ; ;
      p5 = ppMult (p1, p2, h) ;
      p5 = ppMult (p5, pk, h) ;
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }

  if (1)
    {
      p5 = ppPlus (ppNew (NUM, 1, h), ppNew (NUM, -1, h), h) ;
      p5 = ppMult (p5, pk, h) ;
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }
  if (1)
    {
      p5 = ppPlus (ppMult (ppNew (NUM, 1, h), pk, h), ppMult (ppNew (NUM, -1, h), pk, h), h) ;
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }
  if (1)
    {
      p5 = ppPlus (ppNew (ZERO, 0, h), pk, h) ;
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }
  if (1)
    {
      p5 = ppMult (ppNew (ZERO, 0, h), pk, h) ;
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }

  if (1)
    {
      p1 = ppMultiOp (cc, E, "1221", h) ;
      p2 = ppMultiOp (cc, F, "1221", h) ;
      p5 = ppMult (p1, ppMult (p2, pk, h), h) ;
      p5 = ppMult (ppMult (p1,p2,h), pk,h) ;
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }
  if (1)
    {
      p1 = ppMultiOp (cc, E, "122211", h) ;
      p2 = ppMultiOp (cc, F, "112221", h) ;
      p5 = ppMult (ppMult (p1,p2,h), pk,h) ;
      showPp (cc ,p5, 0) ;
      ppExpand (cc, p5, h) ;
      showPp (cc ,p5, 0) ; fprintf (stderr, "\n") ;
    }
  ac_free (h) ;
  return 0 ;
} /* main */

/***********************************************************************/
/***********************************************************************/
