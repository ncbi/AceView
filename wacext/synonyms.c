/*
 * Take a dictionary of words
 * count the frequency of each letter
 * compute the distance as sum square (diff of number of a/b...z)
 * export all cases where the distance is < 3
 */

#include "ac.h"
#include "dict.h"
#include "aceio.h"
#include "utils.h"

#define LIMIT .80
#define DLIMIT 6
#define NALPHA 26

typedef struct sfStruct { int nam, cc[NALPHA*NALPHA*NALPHA] ;} SF ;
typedef struct spStruct { int a, b ; double angle ;} SP ;

/************************************************************/

static void usage (char *err)
{
  fprintf (stderr, "usage synonims -f file_of_words :: returns probable synonyms\n") ;
  exit (1) ;
} /* usage */

/************************************************************/
/* tokenise the input file into a dictionary of words */
static DICT *synParse (const char *fnam) 
{
  DICT *dict = dictCreate (10000) ;
  ACEIN fi = aceInCreateFromFile (fnam, "r", 0, 0) ;
  char *cp, *cq, cutter ;

  if (!fi)
    usage (messprintf ("cannot open fnam", fnam ? fnam : "NULL")) ;

  while (aceInCard (fi)) /* will close f */
    {
      if (1)
	{
	  while ((cp = aceInWordCut (fi, ".,;:?!#@$%*-+=\'\"/(){}[]| \n", &cutter)), cutter)
	    if (cp)
	      {
		cq = cp - 1 ;
		while (*++cq) *cq = ace_lower (*cq) ;
		dictAdd (dict, cp, 0) ;
	      }
	}
      else
	{
	  while ((cp = aceInWordCut (fi, "\n", 0)))
	    {
	      cq = cp - 1 ;
	      while (*++cq) *cq = ace_lower (*cq) ;
	      dictAdd (dict, cp, 0) ;
	    }
	}
    }
  fprintf (stderr, "// found %d words in file %s\n", dictMax (dict), fnam) ;
  return dict ;
} /* synParse */

/************************************************************/

static Array synLetterCount (DICT *dict)
{
  const char *ccp, *ccq ;
  int ii, i, j, k ;
  SF *sf ;

  Array aa = arrayCreate (10000, SF) ;

  for (ii = 1 ; ii <= dictMax (dict) ; ii++)
    {
      sf = arrayp (aa, ii, SF) ;
      sf->nam = ii ;
      ccp = dictName (dict, ii) - 1 ;
      while (*++ccp)
	{
	  ccq = ccp ;
	  i = *ccq++ - 'a' ;
	  j = *ccq++ - 'a' ;
	  k = *ccq++ - 'a' ;
	  if (i >= 0 && i < NALPHA &&
	      j >= 0 && j < NALPHA &&
	      k >= 0 && k < NALPHA)
	    sf->cc[i*NALPHA*NALPHA + j*NALPHA + k]++ ;
	}
    }
  return aa ;
} /* synLetterCount */

/************************************************************/

static Array synFindNeighbours  (DICT *dict, Array aa)
{
  int i, ii, jj, kk, n, na, nb, nd, nda, ndb, nn, nna, nnb ;
   double angle ;
   SP *sp ;
   SF *sfa, *sfb ;
   Array bb = arrayCreate (10000, SP) ;

  for (ii = kk = 0 ; ii < arrayMax (aa) ; ii++)
    {
      sfa = arrp (aa, ii, SF) ;
      for (jj = ii+1 ; jj < arrayMax (aa) ; jj++)
	{
	  sfb = arrp (aa, jj, SF) ;
	  nn = nd = nda = ndb = nna = nnb = 0 ;
	  for (i = 0 ; i < NALPHA*NALPHA*NALPHA ; i++)
	    {
	      na = sfa->cc[i] ;
	      nb = sfb->cc[i] ;
	      n = sfa->cc[i] - sfb->cc[i] ;
	      /* scalar product method */
	      {
		nn += na * nb ;
		nna += na * na ;
		nnb += nb * nb ;
	      }
	      
	      {
		nd += n > 0 ? n : -n ; 
		nda += na  ;
		ndb += nb ;
		if (nd > DLIMIT)
		  break ;
	      }
	    }
	  if (nna == 0 || nnb == 0)
	    continue ;
	  if (0) /* scalar product method */
	    {
	      angle = nn / sqrt (nna * nnb) ;
	      if (angle > LIMIT && angle < 1.0)
		{ /* register happy few */
		  sp = arrayp (bb, kk++, SP) ;
		  sp->a = sfa->nam ; sp->b = sfb->nam ; sp->angle = angle ;
		}
	    }
	  else
	    {
	      angle = nd ;
	      if (angle < DLIMIT && nda > 3 && ndb > 3 && angle > 0)
		{ /* register happy few */
		  sp = arrayp (bb, kk++, SP) ;
		  sp->a = sfa->nam ; sp->b = sfb->nam ; sp->angle = angle ;
		}
	    }
	}
    }
  fprintf (stderr, "// found %d synonym%s\n", kk, kk > 1 ? "s" : "") ;
  return bb ;
} /* synFindNeighbours */

/************************************************************/

static void synShow (DICT *dict, Array bb)
{
  int ii ;
  SP *sp ;
  
  for (ii = 0 ; ii < arrayMax (bb) ; ii++)
    {
      sp = arrp (bb, ii, SP) ;
      printf ("%s\t%s\t%.4f\n", dictName (dict, sp->a),  dictName (dict, sp->b), sp->angle) ;
    } 
  fprintf (stderr, "// found %d synonym%s\n", ii, ii > 1 ? "s" : "") ;
} /* synShow */

/************************************************************/

int main (int argc, const char **argv)
{
 const char *fnam = getArgV (&argc, argv, "-f") ;
 DICT *dict = 0 ;
 Array aa, bb ;

 if (!fnam) usage ("missing file name argument") ;

 dict = synParse (fnam) ;
 aa = synLetterCount (dict) ;
 bb = synFindNeighbours (dict, aa) ;
 synShow (dict, bb) ;

 return 0 ;
}

/************************************************************/
/************************************************************/
