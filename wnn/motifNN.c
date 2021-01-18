#include "ac.h"
#include "matrix.h"
#define NTYPES 12

typedef struct mfStruct {
  int NN ; /* number of motifs */
  MX motifMatrix ;
  Array found ;
  AC_HANDLE h ;
  BOOL gzi, gzo ;
  const char *motifFileName ;
  const char *readFileName ;
  const char *outFileName ;
  DICT *motifDict ;
  Array motifs ;
} MF ;


/****************************************************/

static void mfParseMotifs (MF *mf)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (motifFileName, mf->gzi, h) ;
  DICT *dict = dictHandleCreate (1024, h) ;
  DICT *mDict = mf->motifDict ;
  char *cp ;
  int n, state = 0 ;

  if (!ai)
    messcrash ("Cannot open file motifFile -m %s\n",motifFileName) ;
  
  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (!cp || !cp[1] ||  *cp == '#' || *cp == "/")
	continue ;
      if (state == 0 && cp[0] == '>')
	{
	  dictAdd (mDict, cp+1, &n) ;
	  state = 1 ;
	  continue ;
	}
      if (state == 1)
	{
	  array (lns, n, int) = strlen (cp) ;
	  
	  cp-- ;
	  while (*cp++)
	    switch (*cp)
	      {
	      }
	}	  
    }

  ac_free (h) ;
  return ;
} /* mfParseMotifs */

/****************************************************/
/* create words of length 2kk bits representing all k-mers */
static void mfInitKmers (MF *mf, int kk)
{
  int ii, iMax = 1 << 2 * kk ; /* number of motifs */
  Array aa ;

  mf->mask = iMax - 1 ;

  if (kk > 15)
    messcrash ("kk=%d > 15 in mfInitKmers", kk) ;

  aa = mf->motifs = arrayHandleCreate (iMax, int, h) ;
  for (ii = 0 ; ii < iMax ; ii++)
    {
      int x = 0, y, i = ii ;
      for (j = 0 ; j < kk ; j++)
	{
	  y = i & 0x3 ;
	  x  <<= 2 ; x |= y ; i >>= 2 ;
	}
      array (aa, ii, int) = x ;
    }
  return ;
} /* mfInitKmers */

/****************************************************/

static void mfInit (MF *mf)
{
  mf->motifDict = dictHandleCreate (1024, mf->h) ;
  return ;
} /* mfInit */

/****************************************************/
int main (int argc, const char **argv)
{
  AC_HANDLE h = ac_new_handle () ;
  MF mf ;
  AC_HANDLE h = ac_new_handle () ;
  memset (&mf, 0, sizeof (MF)) ;
  mf.h = h ;  
  mfInit (&mf) ;   /* init the 4x4 matrices */ 
  getCmdLineInt (&argc, argv, "-m", &(mf.motifFileName)) ;
  getCmdLineInt (&argc, argv, "-i", &(mf.readFileName)) ;
  getCmdLineInt (&argc, argv, "-o", &(mf.outFileName)) ;
  mf.gzi = getCmdLineOption (&argc, argv, "-gzi", 0);
  mf.gzo = getCmdLineOption (&argc, argv, "-gzo", 0);

  if (0)
    mfParseMotifs (&mf) ;
  if (1)
    mfInitKmers (mf, 7) ;

  ac_free (h) ;
  return 0 ;
} /* main */

/****************************************************/
/****************************************************/
